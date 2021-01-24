!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp,bcond
   use incomp_class,      only: dirichlet,convective,neumann,clipped_neumann
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single incompressible flow solver and corresponding time tracker
   type(incomp),      public :: fs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
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
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Create an incompressible flow solver with bconds
      create_solver: block
         use ils_class, only: pcg_amg,gmres,amg,pfmg
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',fs%visc)
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
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=gmres)
      end block create_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=fs%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max simulation time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use geometry, only: bus
         type(bcond), pointer :: my_bc
         integer :: n,i,j,k
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
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
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='bus')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('walls',fs%cfg%VF)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
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
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
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
         call mfile%write()
         call cflfile%write()
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
