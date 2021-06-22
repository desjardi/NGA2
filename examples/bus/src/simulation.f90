!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use scalar_class,      only: scalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use datafile_class,    only: datafile
   implicit none
   private
   
   !> Incompressible flow and scalar solvers and corresponding time tracker and sgs model
   type(scalar), dimension(:), allocatable, public :: sc
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs
   
   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,intfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW,Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: resSC
   real(WP), dimension(:,:,:,:), allocatable :: SR
   
   !> Scalar source descriptor
   integer :: nsc
   real(WP), dimension(:,:), allocatable :: all_src_pos
   real(WP), dimension(3) :: src_pos
   real(WP) :: src_rad
   
   !> Fluid viscosity
   real(WP) :: visc
   
contains
   
   
   !> Localizes the vertical vent near the front of the bus (facing -z)
   function vertical_vent(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%x(i) .gt.bus%l_bus-bus%x_kneecurb.and.&
      &   pg%x(i+1).lt.bus%l_bus-bus%x_kneecurb+bus%w_vvent.and.&
      &   pg%y(j) .gt.bus%y_vvent.and.pg%y(j+1).lt.bus%y_vvent+bus%h_vvent.and.&
      &   pg%z(k) .gt.bus%w_bus-bus%w_seatcol-100.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k) .lt.bus%w_bus-bus%w_seatcol+100.0_WP*epsilon(1.0_WP)) isIn=.true.
   end function vertical_vent
   
   
   !> Localizes the 3 floor vents near the front of the bus (facing +y)
   function front_floor_vents(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      ! x and y specifications are the same, 3 regions for z
      if (pg%x(i)  .gt.bus%x_seatdriver+bus%l_seatrow*bus%nseat_driver-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%x(i+1).lt.bus%x_seatdriver+bus%l_seatrow*bus%nseat_driver+bus%l_fventfront+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%y(j)  .le.10.0_WP*epsilon(1.0_WP).and.pg%y(j+1).gt.10.0_WP*epsilon(1.0_WP).and.&
      & ((pg%z(k)  .gt.bus%w_ventwindow-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%w_ventwindow+bus%w_fventfront+10.0_WP*epsilon(1.0_WP)).or. &
      &  (pg%z(k)  .gt.bus%w_bus-bus%w_ventwindow-bus%w_fventfront-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%w_bus-bus%w_ventwindow+10.0_WP*epsilon(1.0_WP)).or. &
      &  (pg%z(k)  .gt.bus%w_bus-bus%w_ventwindow-bus%w_fventfront-bus%w_seatcol-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%w_bus-bus%w_ventwindow-bus%w_seatcol+10.0_WP*epsilon(1.0_WP)))) isIn=.true.
   end function front_floor_vents
   
   
   !> Localizes the back floor vent driver's side (facing +y)
   function back_floor_vent_driver(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%x(i)  .gt.bus%x_seatdriver.and.pg%x(i+1).lt.bus%x_seatdriver+bus%l_fventdriver.and.&
      &   pg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.pg%y(j+1).lt.bus%h_fventback+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).gt.bus%z_fventback+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k  ).lt.bus%z_fventback+10.0_WP*epsilon(1.0_WP)) isIn=.true.
   end function back_floor_vent_driver
   
   
   !> Localizes the back floor vent curb side (facing +y)
   function back_floor_vent_curb(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%x(i)  .gt.bus%x_seatcurb+2.0_WP*bus%l_seatrow.and.&
      &   pg%x(i+1).lt.bus%x_seatcurb+2.0_WP*bus%l_seatrow+bus%l_fventcurb.and.&
      &   pg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.pg%y(j+1).lt.bus%h_fventback+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k)  .gt.bus%w_bus-bus%z_fventback-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k)  .lt.bus%w_bus-bus%z_fventback+10.0_WP*epsilon(1.0_WP)) isIn=.true.
   end function back_floor_vent_curb
   
   
   !> Localizes the vent(s) at the lavatory (facing +x)
   function lavatory_vent(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%x(i)  .le.bus%x_lav+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%x(i+1).gt.bus%x_lav+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%y(j)  .gt.bus%h_seat+bus%h_bseat-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%y(j+1).lt.bus%h_seat+bus%h_bseat+bus%h_ventlav+10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k)  .gt.bus%w_bus-bus%w_seatcol+bus%w_seat-bus%w_ventlav.and.&
      &   pg%z(k+1).lt.bus%w_bus-bus%w_seatcol+bus%w_seat) isIn=.true.
   end function lavatory_vent
   
   
   !> Localizes the vents along the windows on both sides (facing +y)
   function window_vents(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      ! driver side
      if (pg%x(i)  .gt.bus%x_seatdriver.and.pg%x(i+1).lt.bus%l_bus-bus%x_closet.and.&
      &   pg%y(j)  .le.bus%y_ventwindow.and.pg%y(j+1).gt.bus%y_ventwindow.and.&
      &   pg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%w_ventwindow+10.0_WP*epsilon(1.0_WP)) isIn=.true.
      ! curb side
      if (pg%x(i)  .gt.bus%x_seatcurb.and.pg%x(i+1).lt.bus%l_bus-bus%x_kneecurb.and.&
      &   pg%y(j)  .le.bus%y_ventwindow.and.pg%y(j+1).gt.bus%y_ventwindow.and.&
      &   pg%z(k)  .gt.bus%w_bus-bus%w_ventwindow-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).le.bus%w_bus+10.0_WP*epsilon(1.0_WP)) isIn=.true.
   end function window_vents
   
   
   !> Localizes the intake vents of the parcel rack (facing -y)
   function rackintake_vents(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      ! same for x and y, two regions for z
      if (pg%x(i)  .gt.bus%x_vinrack-0.5_WP*bus%l_vinrack-100.0_WP*epsilon(1.0_WP).and.&
      &   pg%x(i+1).lt.bus%x_vinrack+0.5_WP*bus%l_vinrack+100.0_WP*epsilon(1.0_WP).and.&
      &   pg%y(j)  .gt.bus%h_rack-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%y(j)  .lt.bus%h_rack+10.0_WP*epsilon(1.0_WP).and.&
      & ((pg%z(k)  .gt.bus%min_length/real(bus%gres,WP)-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%min_length/real(bus%gres,WP)+10.0_WP*epsilon(1.0_WP)+bus%w_vinrack).or.&
      &  (pg%z(k)  .gt.bus%w_bus-bus%w_vinrack-bus%min_length/real(bus%gres,WP)-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%w_bus-bus%min_length/real(bus%gres,WP)+10.0_WP*epsilon(1.0_WP)))) isIn=.true.
   end function rackintake_vents
   
   
   !> Localizes the outlet vents of the parcel rack (facing -y)
   function rackoutlet_vents(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      use geometry, only: bus
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      integer :: n
      logical :: isIn
      real(WP):: x1face,adist
      isIn=.false.
      ! get factor to adjust spacing according to mesh
      adist = real(mod(nint(bus%l_seatrow/bus%i2m),nint(bus%min_length/bus%gres/bus%i2m)),WP)*bus%i2m
      ! same y
      if (pg%y(j)  .gt.bus%h_rack-10.0_WP*epsilon(1.0_WP).and.&
      &   pg%y(j)  .lt.bus%h_rack+10.0_WP*epsilon(1.0_WP)) isIn=.true.
      ! continue if proper y location
      if (isIn) then; isIn=.false.; else; return
      end if
      ! driver side, cycle through rows
      x1face = bus%x_seatdriver+real(bus%nseat_driver-1,WP)*bus%l_seatrow+bus%th_seat
      do n=1,bus%nseat_driver
         if (pg%z(k)  .gt.bus%z_ventface1-100.0_WP*epsilon(1.0_WP).and.&
         &   pg%z(k+1).lt.bus%z_ventface2+100.0_WP*epsilon(1.0_WP).and.&
         &   pg%x(i)  .gt.x1face+bus%x_v2facedriver-1000.0_WP*epsilon(1.0_WP) &
         &   -real(n-1,WP)*bus%l_ventrow+real(mod(n,2),WP)*adist .and.&
         &   pg%x(i+1).lt.x1face+bus%x_v2facedriver+1000.0_WP*epsilon(1.0_WP) &
         &   -real(n-1,WP)*bus%l_ventrow+real(mod(n,2),WP)*adist+bus%l_ventface) isIn=.true.
      end do
      ! curb side, front section
      x1face = bus%x_seatcurb2+real(bus%nseat_curb2-1,WP)*bus%l_seatrow+bus%th_seat
      do n=1,bus%nseat_curb2-1
         if (pg%z(k)  .gt.bus%w_bus-bus%z_ventface2-100.0_WP*epsilon(1.0_WP).and.&
         &   pg%z(k+1).lt.bus%w_bus-bus%z_ventface1+100.0_WP*epsilon(1.0_WP).and.&
         &   pg%x(i)  .gt.x1face+bus%x_v2facecurb-1000.0_WP*epsilon(1.0_WP) &
         &   -real(n-1,WP)*bus%l_ventrow+real(mod(n,2),WP)*adist.and.&
         &   pg%x(i+1).lt.x1face+bus%x_v2facecurb+1000.0_WP*epsilon(1.0_WP) &
         &   -real(n-1,WP)*bus%l_ventrow+real(mod(n,2),WP)*adist+bus%l_ventface) isIn=.true.
      end do
      ! curb side, backward seats
      if (pg%z(k)  .gt.bus%w_bus-bus%z_ventface2-100.0_WP*epsilon(1.0_WP).and.&
      &   pg%z(k+1).lt.bus%w_bus-bus%z_ventface1+100.0_WP*epsilon(1.0_WP).and.&
      &   pg%x(i)  .gt.bus%x_seatcurb2+bus%l_seat-bus%th_seat-bus%x_v2facecurb &
      &   -bus%l_ventface-1000.0_WP*epsilon(1.0_WP).and.  pg%x(i+1).lt.bus%x_seatcurb2 &
      &   +bus%l_seat-bus%th_seat-bus%x_v2facecurb+1000.0_WP*epsilon(1.0_WP)) isIn=.true.
      ! curb side, back section
      x1face = bus%x_seatcurb+real(bus%nseat_curb1-1,WP)*bus%l_seatrow+bus%th_seat
      do n=1,bus%nseat_curb1
         if (pg%z(k)  .gt.bus%w_bus-bus%z_ventface2-100.0_WP*epsilon(1.0_WP).and.&
         &   pg%z(k+1).lt.bus%w_bus-bus%z_ventface1+100.0_WP*epsilon(1.0_WP).and.&
         &   pg%x(i)  .gt.x1face+bus%x_v2facecurb-1000.0_WP*epsilon(1.0_WP) &
         &   -real(n-1,WP)*bus%l_ventrow+real(mod(n+1,2),WP)*adist.and.&
         &   pg%x(i+1).lt.x1face+bus%x_v2facecurb+1000.0_WP*epsilon(1.0_WP) &
         &   -real(n-1,WP)*bus%l_ventrow+real(mod(n+1,2),WP)*adist+bus%l_ventface) isIn=.true.
      end do
   end function rackoutlet_vents
   
   
   !> Functions that localizes the source of tracer scalar
   function scalar_src(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (sqrt((pg%xm(i)-src_pos(1))**2+(pg%ym(j)-src_pos(2))**2++(pg%zm(k)-src_pos(3))**2).le.src_rad) isIn=.true.
   end function scalar_src
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Start by reading in the source positions
      initialize_sources: block
         use param, only: param_getsize
         ! Figure out how many sources
         nsc=param_getsize('Source position')/3
         if (cfg%amRoot) print*,'Found',nsc,'sources!'
         ! Allocate storage for sources and read them
         if (nsc.gt.0) then
            allocate(all_src_pos(3,nsc))
            call param_read('Source position',all_src_pos)
            call param_read('Source radius',src_rad)
         end if
      end block initialize_sources
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='bustime')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max simulation time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Handle restart/saves here
      restart_and_save: block ! CAREFUL - WE NEED TO CREATE THE TIMETRACKER BEFORE THE EVENT
         character(len=str_medium) :: dir_restart
         ! Create event for saving restart files
         save_evt=event(time,'Restart output')
         call param_read('Restart output period',save_evt%tper)
         ! Check if we are restarting
         call param_read(tag='Restart from',val=dir_restart,short='r',default='')
         restarted=.false.; if (len_trim(dir_restart).gt.0) restarted=.true.
         if (restarted) then
            ! If we are, read the name of the directory
            call param_read('Restart from',dir_restart,'r')
            ! Read the datafile
            df=datafile(pg=cfg,fdata=trim(adjustl(dir_restart))//'/'//'data')
         else
            ! If we are not restarting, we will still need a datafile for saving restart files
            ! We're only saving the velocity/pressure here, not the tracers
            df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=6)
            df%valname(1)='t'
            df%valname(2)='dt'
            df%varname(1)='U'
            df%varname(2)='V'
            df%varname(3)='W'
            df%varname(4)='P'
            df%varname(5)='LM'
            df%varname(6)='MM'
         end if
      end block restart_and_save
      
      
      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
         if (restarted) then
            call df%pullval(name='t' ,val=time%t )
            call df%pullval(name='dt',val=time%dt)
            time%told=time%t-time%dt
         end if
      end block update_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_solver: block
         use ils_class, only: pcg_amg,pcg_pfmg,gmres
         use incomp_class, only: dirichlet,convective,neumann,clipped_neumann
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         !! -- Main HVAC system -- !!
         ! Intakes (outflows of fluid domain)
         call fs%add_bcond(name='vertical vent'    ,type=dirichlet,locator=vertical_vent         ,face='z',dir=+1,canCorrect=.false.)
         call fs%add_bcond(name='front floor vents',type=dirichlet,locator=front_floor_vents     ,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='back vent driver' ,type=dirichlet,locator=back_floor_vent_driver,face='z',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='back vent curb'   ,type=dirichlet,locator=back_floor_vent_curb  ,face='z',dir=+1,canCorrect=.false.)
         call fs%add_bcond(name='lavatory vent'    ,type=clipped_neumann,locator=lavatory_vent   ,face='x',dir=-1,canCorrect=.true. )
         ! Outputs (inflows of fluid domain)
         call fs%add_bcond(name='window vents'     ,type=dirichlet,locator=window_vents          ,face='y',dir=-1,canCorrect=.false.)
         !! -- Parcel rack system -- !!
         ! Intakes (outflows of fluid domain)
         call fs%add_bcond(name='rack intakes'     ,type=dirichlet,locator=rackintake_vents      ,face='y',dir=+1,canCorrect=.false.)
         ! Outputs (inflows of fluid domain)
         call fs%add_bcond(name='rack outlets'     ,type=dirichlet,locator=rackoutlet_vents      ,face='y',dir=+1,canCorrect=.false.)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         fs%psolv%maxlevel=16
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=gmres)
      end block create_solver
      
      
      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         ! Handle restart
         if (restarted) then
            call df%pullvar(name='LM',var=sgs%LM)
            call df%pullvar(name='MM',var=sgs%MM)
         end if
      end block create_sgs
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use geometry, only: bus
         use incomp_class, only: bcond
         type(bcond), pointer :: my_bc
         integer :: n,i,j,k
         ! Initial fields
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! Handle restart
         if (restarted) then
            call df%pullvar(name='U',var=fs%U)
            call df%pullvar(name='V',var=fs%V)
            call df%pullvar(name='W',var=fs%W)
            call df%pullvar(name='P',var=fs%P)
         end if
         ! Apply Dirichlet at the inflow vents and known outflow vents
         call fs%get_bcond('vertical vent',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%W(i,j,k)=bus%vel_vvent
         end do
         call fs%get_bcond('front floor vents',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%V(i,j,k)=bus%vel_fventfront !< Experimental data differs between the vents, may need to separate
         end do
         call fs%get_bcond('back vent driver',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%W(i,j,k)=bus%vel_fventdriver
         end do
         call fs%get_bcond('back vent curb',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%W(i,j,k)=bus%vel_fventcurb
         end do
         call fs%get_bcond('window vents',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%V(i,j,k)=bus%vel_ventwindow !< This is an average value, the measurements vary significantly
         end do
         call fs%get_bcond('rack intakes',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%V(i,j,k)=bus%vel_vinrack
         end do
         call fs%get_bcond('rack outlets',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            fs%V(i,j,k)=bus%vel_ventface !< This value has been balanced with the rack intake flowrate
         end do
         ! Apply all other boundary conditions (just the lavatory is neumann)
         call fs%apply_bcond(time%t,time%dt)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity
      
      
      ! Create nsc scalar solvers with bconds
      create_scalar: block
         use ils_class, only: gmres
         use scalar_class, only: bcond,dirichlet,neumann,quick
         type(bcond), pointer :: mybc
         integer :: i,j,k,n,ii
         character(len=2) :: id
         ! Allocate scalar solvers
         if (nsc.gt.0) allocate(sc(nsc))
         ! For each solver, prepare it completely
         do ii=1,nsc
            ! Prepare tracer name
            write(id,'(i2.2)') ii
            ! Create solver
            sc(ii)=scalar(cfg=cfg,scheme=quick,name='Tracker'//id)
            ! Assign density and diffusivity = viscosity
            sc(ii)%rho=fs%rho; sc(ii)%diff=visc
            ! Configure implicit scalar solver
            sc(ii)%implicit%maxit=fs%implicit%maxit; sc(ii)%implicit%rcvg=fs%implicit%rcvg
            ! Assign proper source location
            src_pos=all_src_pos(:,ii)
            ! Define boundary conditions
            call sc(ii)%add_bcond(name='source',type=dirichlet,locator=scalar_src,dir='+y')
            ! Setup the solver
            call sc(ii)%setup(implicit_ils=gmres)
            ! Initialize the scalar fields at zero except for source patch
            sc(ii)%SC=0.0_WP
            call sc(ii)%get_bcond('source',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               sc(ii)%SC(i,j,k)=1.0_WP
            end do
         end do
      end block create_scalar
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='bus')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         if (nsc.ge.1) call ens_out%add_scalar('tracer',sc(1)%SC)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: ii
         character(len=2) :: id
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create sc integral monitor
         intfile=monitor(fs%cfg%amRoot,'integral')
         call intfile%add_column(time%n,'Timestep number')
         call intfile%add_column(time%t,'Time')
         do ii=1,nsc
            call sc(ii)%get_max(); call sc(ii)%get_int()
            write(id,'(i2.2)') ii
            call intfile%add_column(sc(ii)%SCint,'SC'//id)
         end do
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      integer :: ii
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old scalar
         do ii=1,nsc
            sc(ii)%SCold=sc(ii)%SC
         end do
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Reset here fluid properties
         fs%visc=visc
         
         ! Turbulence modeling
         call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         resU=fs%rho
         call sgs%get_visc(dt=time%dtold,rho=resU,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         where (sgs%visc.lt.-fs%visc)
            sgs%visc=-fs%visc
         end where
         fs%visc=fs%visc+sgs%visc
         do ii=1,nsc
            sc(ii)%diff=fs%visc
         end do
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time scalar
            do ii=1,nsc
               sc(ii)%SC=0.5_WP*(sc(ii)%SC+sc(ii)%SCold)
            end do
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! ============= SCALAR SOLVER =======================
            resU=fs%rho*fs%U; resV=fs%rho*fs%V; resW=fs%rho*fs%W
            do ii=1,nsc
               ! Explicit calculation of drhoSC/dt from scalar equation
               call sc(ii)%get_drhoSCdt(drhoSCdt=resSC,rhoU=resU,rhoV=resV,rhoW=resW)
               ! Assemble explicit residual
               resSC=-2.0_WP*(sc(ii)%rho*sc(ii)%SC-sc(ii)%rho*sc(ii)%SCold)+time%dt*resSC
               ! Form implicit residual
               call sc(ii)%solve_implicit(dt=time%dt,resSC=resSC,rhoU=resU,rhoV=resV,rhoW=resW)
               ! Apply this residual
               sc(ii)%SC=2.0_WP*sc(ii)%SC-sc(ii)%SCold+resSC
               ! Apply other boundary conditions on the resulting field
               call sc(ii)%apply_bcond(time%t,time%dt)
            end do
            ! ===================================================
            
            
            ! ============ VELOCITY SOLVER ======================
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         do ii=1,nsc
            call sc(ii)%get_max(); call sc(ii)%get_int()
         end do
         call mfile%write()
         call cflfile%write()
         call intfile%write()
         
         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
               character(len=str_medium) :: dirname,timestamp
               ! Prefix for files
               dirname='restart_'; write(timestamp,'(es12.5)') time%t
               ! Prepare a new directory
               if (fs%cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
               ! Populate df and write it
               call df%pushval(name='t' ,val=time%t )
               call df%pushval(name='dt',val=time%dt)
               call df%pushvar(name='U' ,var=fs%U   )
               call df%pushvar(name='V' ,var=fs%V   )
               call df%pushvar(name='W' ,var=fs%W   )
               call df%pushvar(name='P' ,var=fs%P   )
               call df%pushvar(name='LM',var=sgs%LM )
               call df%pushvar(name='MM',var=sgs%MM )
               call df%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data')
            end block save_restart
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,resSC,Ui,Vi,Wi,SR)
      
   end subroutine simulation_final
   
   
end module simulation
