!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg,Lv,IR,lpipe,ypipe,rpipe,dx,dy,dz
   use lowmach_class,     only: lowmach
   use sgsmodel_class,    only: sgsmodel
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single incompressible flow solver and corresponding time tracker
   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Provide datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   
   !> Equation of state and case conditions
   real(WP) :: pressure,pressure_old,Vtotal
   real(WP) :: MFR,Ain,rhoUin,fluid_mass,fluid_mass_old
   real(WP) :: Tinit,Tinlet,Twall,Twallold,Tavg
   logical  :: wall_losses
   real(WP), parameter :: Wmlr=44.01e-3_WP  ! kg/mol
   real(WP), parameter :: Rcst=8.314_WP     ! J/(mol.K)
   real(WP), parameter :: Cp=40.0_WP/Wmlr   ! ~40 J/(mol.K) from NIST, divided by Wmlr to get to kg
   real(WP), parameter :: Pr_turb=0.9_WP    ! For now, we're assuming a constant Prandtl number
   real(WP), parameter :: Cs=500.0_WP       ! Specific heat of steel in J/(kg.K)
   real(WP), parameter :: Msteel=268.5_WP   ! Mass of steel in the sled/basket assembly
   real(WP), parameter :: Lsteel=0.0762_WP  ! Thickness of the wall in m
   real(WP), parameter :: kwall=2.0_WP      ! W/(K.m) <- lower than steel because of Teflon
   real(WP), parameter :: tau_wall_in = 10.0_WP   ! Thermal timescale for wall heating by inside  (in s)
   real(WP), parameter :: tau_wall_out=100.0_WP   ! Thermal timescale for wall cooling by outside (in s)
   real(WP), parameter :: Tout=300.0_WP     ! Outside temp
   
   !> Info relevant to basket modeling
   real(WP), parameter :: bporo=0.4_WP      ! Assuming 40% porosity in the bag - depends on shape and arrangement of product
   real(WP), parameter :: bperm=1.0e-7_WP   ! Permeability - should be low, needs to be optimized
   real(WP), parameter :: Tperm=1.0e-6_WP   ! Thermal permeability - should be low, needs to be optimized
   real(WP), dimension(3), parameter :: bsize=[0.40_WP,0.05_WP,0.20_WP]  ! Dimensions of a product bag
   real(WP) :: be                                                        ! Smearing coeff for a product bag
   real(WP), dimension(:,:,:), allocatable :: epsp                       ! Volume fraction of product
   real(WP), dimension(:,:,:), allocatable :: epsf                       ! Volume fraction of fluid
   real(WP), dimension(:,:,:), allocatable :: Tprod,Tprodold             ! Temperature of product
   real(WP), parameter :: Cprod=3500.0_WP   ! In J/(kg.K) - this is for potatoes
   real(WP), parameter :: kprod=0.7_WP      ! In W/(m.K)  - this is a representative value for potatoes [0.545,0.957]
   real(WP), parameter :: rprod=1070.0_WP   ! Mass density of product
   real(WP), parameter :: Mbag=2.26796_WP   ! Weight of 1 bag of product
   integer , parameter :: nbag=55           ! Number of bags
   real(WP), parameter :: Mprod=nbag*Mbag   ! Weight of product
   
   ! Backup of viscosity and diffusivity
   real(WP), dimension(:,:,:), allocatable :: viscmol,diffmol
   
contains
   
   
   !> Function that provides volume fraction of bag at a point
   function bag_at_loc(bag_loc,mesh_loc) result(eps)
      implicit none
      real(WP), dimension(3), intent(in) :: bag_loc,mesh_loc
      real(WP) :: eps
      real(WP), dimension(3) :: dist
      dist=min(mesh_loc-(bag_loc-0.5_WP*bsize),(bag_loc+0.5_WP*bsize)-mesh_loc)
      eps=(0.5_WP+0.5_WP*tanh(dist(1)/be))*(0.5_WP+0.5_WP*tanh(dist(2)/be))*(0.5_WP+0.5_WP*tanh(dist(3)/be))
   end function bag_at_loc
   
   
   !> Function that localizes the top of the tube - y-face (v-vel)
   function vtube(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: r
      isIn=.false.
      if (abs(pg%xm(i)+dx-0.5_WP*lpipe).le.dx.and.pg%y(j)-0.5_WP*dy.lt.ypipe+rpipe.and.pg%y(j)+0.5_WP*dy.ge.ypipe+rpipe.and.abs(pg%zm(k)).lt.rpipe) isIn=.true.
      if (abs(pg%xm(i)                ).le.dx.and.pg%y(j)-0.5_WP*dy.lt.ypipe+rpipe.and.pg%y(j)+0.5_WP*dy.ge.ypipe+rpipe.and.abs(pg%zm(k)).lt.rpipe) isIn=.true.
      if (abs(pg%xm(i)-dx+0.5_WP*lpipe).le.dx.and.pg%y(j)-0.5_WP*dy.lt.ypipe+rpipe.and.pg%y(j)+0.5_WP*dy.ge.ypipe+rpipe.and.abs(pg%zm(k)).lt.rpipe) isIn=.true.
   end function vtube
   
   
   !> Function that localizes the top of the tube - scalar
   function sctube(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: r
      isIn=.false.
      if (abs(pg%xm(i)+dx-0.5_WP*lpipe).le.dx.and.pg%ym(j).lt.ypipe+rpipe.and.pg%ym(j)+dy.ge.ypipe+rpipe.and.abs(pg%zm(k)).lt.rpipe) isIn=.true.
      if (abs(pg%xm(i)                ).le.dx.and.pg%ym(j).lt.ypipe+rpipe.and.pg%ym(j)+dy.ge.ypipe+rpipe.and.abs(pg%zm(k)).lt.rpipe) isIn=.true.
      if (abs(pg%xm(i)-dx+0.5_WP*lpipe).le.dx.and.pg%ym(j).lt.ypipe+rpipe.and.pg%ym(j)+dy.ge.ypipe+rpipe.and.abs(pg%zm(k)).lt.rpipe) isIn=.true.
   end function sctube
   
   
   !> Function that localizes the left end of the tube
   !function left_of_tube(pg,i,j,k) result(isIn)
   !   use pgrid_class, only: pgrid
   !   implicit none
   !   class(pgrid), intent(in) :: pg
   !   integer, intent(in) :: i,j,k
   !   logical :: isIn
   !   real(WP) :: r
   !   isIn=.false.
   !   r=sqrt((pg%ym(j)-ypipe)**2+(pg%zm(k))**2)
   !   if (abs(pg%x(i)+0.5_WP*lpipe).lt.0.01_WP.and.r.lt.rpipe) isIn=.true.
   !end function left_of_tube
   
   
   !> Function that localizes the right end of the tube - x-face (u-vel)
   !function right_of_tube_vel(pg,i,j,k) result(isIn)
   !   use pgrid_class, only: pgrid
   !   implicit none
   !   class(pgrid), intent(in) :: pg
   !   integer, intent(in) :: i,j,k
   !   logical :: isIn
   !   real(WP) :: r
   !   isIn=.false.
   !   r=sqrt((pg%ym(j)-ypipe)**2+(pg%zm(k))**2)
   !   if (abs(pg%x(i)-0.5_WP*lpipe).lt.0.01_WP.and.r.lt.rpipe) isIn=.true.
   !end function right_of_tube_vel
   
   
   !> Function that localizes the right end of the tube - scalar
   !function right_of_tube_sc(pg,i,j,k) result(isIn)
   !   use pgrid_class, only: pgrid
   !   implicit none
   !   class(pgrid), intent(in) :: pg
   !   integer, intent(in) :: i,j,k
   !   logical :: isIn
   !   real(WP) :: r
   !   isIn=.false.
   !   r=sqrt((pg%ym(j)-ypipe)**2+(pg%zm(k))**2)
   !   if (abs(pg%x(i+1)-0.5_WP*lpipe).lt.0.01_WP.and.r.lt.rpipe) isIn=.true.
   !end function right_of_tube_sc
   
   
   !> Function that localizes the vessel walls
   function wall_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: r
      isIn=.false.
      r=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (pg%xm(i).le.-0.5_WP*Lv.or.pg%xm(i).ge.+0.5_WP*Lv.or.r.ge.IR) isIn=.true.
   end function wall_locator
   
   
   !> Define here our equation of state - rho(T,mass)
   subroutine get_rho(mass)
      use vdscalar_class, only: bcond
      implicit none
      real(WP), intent(in) :: mass
      type(bcond), pointer :: inflow
      integer :: i,j,k,n
      real(WP) :: one_over_T
      ! Integrate 1/T
      resSC=1.0_WP/sc%SC
      call sc%cfg%integrate(resSC,integral=one_over_T)
      ! Update pressure first
      pressure=mass*Rcst/(Wmlr*one_over_T)
      ! Update density in the domain
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               if (sc%cfg%VF(i,j,k).gt.0.0_WP) then
                  sc%rho(i,j,k)=pressure*Wmlr/(Rcst*sc%SC(i,j,k))
               else
                  sc%rho(i,j,k)=1.0_WP
               end if
            end do
         end do
      end do
      ! Also update the density in the bcond
      !call sc%get_bcond( 'left inflow',inflow)
      !do n=1,inflow%itr%no_
      !   i=inflow%itr%map(1,n); j=inflow%itr%map(2,n); k=inflow%itr%map(3,n)
      !   sc%rho(i,j,k)=pressure*Wmlr/(Rcst*sc%SC(i,j,k))
      !end do
      !call sc%get_bcond('right inflow',inflow)
      !do n=1,inflow%itr%no_
      !   i=inflow%itr%map(1,n); j=inflow%itr%map(2,n); k=inflow%itr%map(3,n)
      !   sc%rho(i,j,k)=pressure*Wmlr/(Rcst*sc%SC(i,j,k))
      !end do
      call sc%get_bcond('pipe inflow',inflow)
      do n=1,inflow%itr%no_
         i=inflow%itr%map(1,n); j=inflow%itr%map(2,n); k=inflow%itr%map(3,n)
         sc%rho(i,j,k)=pressure*Wmlr/(Rcst*sc%SC(i,j,k))
      end do
      ! Account for porosity here
      sc%rho=sc%rho*epsf
   end subroutine get_rho
   
   
   !> Calculate here our viscosity from local T and vessel pressure
   subroutine get_visc()
      implicit none
      real(WP), parameter :: A1=-1.146067e-01_WP
      real(WP), parameter :: A2=+6.978380e-07_WP
      real(WP), parameter :: A3=+3.976765e-10_WP
      real(WP), parameter :: A4=+6.336120e-02_WP
      real(WP), parameter :: A5=-1.166119e-02_WP
      real(WP), parameter :: A6=+7.142596e-04_WP
      real(WP), parameter :: A7=+6.519333e-06_WP
      real(WP), parameter :: A8=-3.567559e-01_WP
      real(WP), parameter :: A9=+3.180473e-02_WP
      integer :: i,j,k
      real(WP) :: Pbar,lnT,MUcp
      ! Pressure needs to be in bar and in the model range
      Pbar=max(min(1.0e-5_WP*pressure,1014.0_WP),75.0_WP)
      ! Loop over the entire domain
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               ! Log of temperature clipped to the model range
               lnT=log(min(max(sc%SC(i,j,k),305.0_WP),900.0_WP))
               ! Evaluate viscosity in cP
               MUcp=(A1+A2*Pbar+A3*Pbar**2+A4*lnT+A5*lnT**2+A6*lnT**3)/(1.0_WP+A7*Pbar+A8*lnT+A9*lnT**2)
               ! Set viscosity
               fs%visc(i,j,k)=MUcp/1000.0_WP
            end do
         end do
      end do
   end subroutine get_visc
   
   
   !> Calculate here our thermal conductivity from local T and rho
   subroutine get_cond()
      implicit none
      real(WP), parameter :: A1=-105.161_WP
      real(WP), parameter :: A2=+0.9007_WP
      real(WP), parameter :: A3=+0.0007_WP
      real(WP), parameter :: A4=+3.50e-15_WP
      real(WP), parameter :: A5=+3.76e-10_WP
      real(WP), parameter :: A6=+0.7500_WP
      real(WP), parameter :: A7=+0.0017_WP
      integer :: i,j,k
      real(WP) :: T,rho,cond
      ! Loop over the entire domain
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               ! Temperature clipped to the model range
               T=min(max(sc%SC(i,j,k),290.0_WP),800.0_WP)
               ! Density clipped to the model range
               rho=min(max(sc%rho(i,j,k)/epsf(i,j,k),1.0_WP),1200.0_WP)
               ! Evaluate conductivity in mW/(m.K)
               cond=(A1+A2*rho+A3*rho**2+A4*rho**3*T**3+A5*rho**4+A6*T+A7*T**2)/sqrt(T)
               ! Set heat diffusivity
               sc%diff(i,j,k)=0.001_WP*cond/Cp
            end do
         end do
      end do
   
end subroutine get_cond
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Handle restart/saves here
      restart_and_save: block
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
            ! Read the two datafiles and the name of the IRL file to read later
            df=datafile(pg=cfg,fdata=trim(adjustl(dir_restart))//'/'//'data')
         else
            ! If we are not restarting, we will still need datafiles for saving restart files
            df=datafile(pg=cfg,filename=trim(cfg%name),nval=4,nvar=9)
            df%valname(1)='t'
            df%valname(2)='dt'
            df%valname(3)='pressure'
            df%valname(4)='Twall'
            df%varname(1)='rhoU'
            df%varname(2)='rhoV'
            df%varname(3)='rhoW'
            df%varname(4)='P'
            df%varname(5)='rho_fs'
            df%varname(6)='rho'
            df%varname(7)='rhoold'
            df%varname(8)='SC'
            df%varname(9)='Tprod'
         end if
      end block restart_and_save
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         ! Flow solver
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Scalar solver
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
         ! Handle restart
         if (restarted) then
            call df%pullval(name='t' ,val=time%t )
            call df%pullval(name='dt',val=time%dt)
            time%told=time%t-time%dt
         end if
      end block initialize_timetracker
      
      
      ! Create an incompressible flow solver with bconds
      create_solver: block
         use ils_class,     only: pcg_pfmg,pcg_amg
         use lowmach_class, only: dirichlet
         real(WP) :: visc
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Define boundary conditions
         !call fs%add_bcond(name= 'left inflow',type=dirichlet,locator= left_of_tube    ,face='x',dir=+1,canCorrect=.false.)
         !call fs%add_bcond(name='right inflow',type=dirichlet,locator=right_of_tube_vel,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name= 'pipe inflow',type=dirichlet,locator=vtube,face='y',dir=-1,canCorrect=.false.)
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=pcg_pfmg)
      end block create_solver
      
      
      ! Allocate and initialize product volume fraction
      ! This is done ASAP but after fs%setup so we can directly modify the Laplacian
      prepare_eps: block
         integer :: i,j,k
         ! Set smearing thickness based on mesh size
         be=0.75_WP*cfg%min_meshsize
         ! Allocate volume fractions and permeability
         allocate(epsp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); epsp=0.0_WP
         allocate(epsf(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); epsf=0.0_WP
         ! Build volume fraction based on sum of bags
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  
                  ! First 7-bag test
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.0_WP, 0.0_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.5_WP, 0.0_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.5_WP, 0.0_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.0_WP,+0.2_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.0_WP,-0.2_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.0_WP, 0.0_WP,-0.3_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  !epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.0_WP, 0.0_WP,+0.3_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])
                  
                  ! 55-bag version from Amir
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP, 0.000_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! Center bag  1
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP, 0.000_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! X direction 2
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP, 0.000_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 3
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP, 0.000_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 4
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP, 0.000_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 5
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,+0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 6
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,-0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 7
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,+0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 8
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,+0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 9
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,-0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 10
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,-0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 11
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,+0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 12
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,+0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 13
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,-0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 14
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,-0.075_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 15
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,+0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 16
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,+0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 17
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,-0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 18
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,-0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 19
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,+0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 20
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,-0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 21
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,+0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 22
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,-0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 23
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,+0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 24
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,+0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 25
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,-0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 26
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,-0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 27
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,+0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 28
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,+0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 29
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,-0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 30
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,-0.150_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 31
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,+0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 32
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,+0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 33
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,-0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 34
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,-0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 35
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,+0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 36
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,-0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 37
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,+0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 38
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,+0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 39
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.45_WP,-0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 40
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.45_WP,-0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 41
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,+0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 42
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[ 0.00_WP,-0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 43
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,+0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 44
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,+0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 45
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,-0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 46
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,-0.225_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 47
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,+0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 48
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,+0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 49
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,-0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 50
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,-0.300_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 51
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,+0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 52
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,+0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 53
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[+0.90_WP,-0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 54
                  epsp(i,j,k)=epsp(i,j,k)+bag_at_loc(bag_loc=[-0.90_WP,-0.375_WP, 0.0_WP],mesh_loc=[cfg%xm(i),cfg%ym(j),cfg%zm(k)])  ! 55
                  
               end do
            end do
         end do
         ! Multiply by expected porosity in the bag
         epsp=epsp*bporo
         epsf=1.0_WP!-epsp     !< SET TO ONE TO REMOVE POROSITY EFFECT
         ! Setup the scaled Laplacian operator from incomp metrics: lap(*)=-vol*div(epsf*grad(*))
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  ! Set Laplacian
                  fs%psolv%opr(1,i,j,k)=fs%divp_x(1,i,j,k)*sum(fs%itpr_x(:,i+1,j,k)*epsf(i:i+1,j,k))*fs%divu_x(-1,i+1,j,k)+&
                  &                     fs%divp_x(0,i,j,k)*sum(fs%itpr_x(:,i  ,j,k)*epsf(i-1:i,j,k))*fs%divu_x( 0,i  ,j,k)+&
                  &                     fs%divp_y(1,i,j,k)*sum(fs%itpr_y(:,i,j+1,k)*epsf(i,j:j+1,k))*fs%divv_y(-1,i,j+1,k)+&
                  &                     fs%divp_y(0,i,j,k)*sum(fs%itpr_y(:,i,j  ,k)*epsf(i,j-1:j,k))*fs%divv_y( 0,i,j  ,k)+&
                  &                     fs%divp_z(1,i,j,k)*sum(fs%itpr_z(:,i,j,k+1)*epsf(i,j,k:k+1))*fs%divw_z(-1,i,j,k+1)+&
                  &                     fs%divp_z(0,i,j,k)*sum(fs%itpr_z(:,i,j,k  )*epsf(i,j,k-1:k))*fs%divw_z( 0,i,j,k  )
                  fs%psolv%opr(2,i,j,k)=fs%divp_x(1,i,j,k)*sum(fs%itpr_x(:,i+1,j,k)*epsf(i:i+1,j,k))*fs%divu_x( 0,i+1,j,k)
                  fs%psolv%opr(3,i,j,k)=fs%divp_x(0,i,j,k)*sum(fs%itpr_x(:,i  ,j,k)*epsf(i-1:i,j,k))*fs%divu_x(-1,i  ,j,k)
                  fs%psolv%opr(4,i,j,k)=fs%divp_y(1,i,j,k)*sum(fs%itpr_y(:,i,j+1,k)*epsf(i,j:j+1,k))*fs%divv_y( 0,i,j+1,k)
                  fs%psolv%opr(5,i,j,k)=fs%divp_y(0,i,j,k)*sum(fs%itpr_y(:,i,j  ,k)*epsf(i,j-1:j,k))*fs%divv_y(-1,i,j  ,k)
                  fs%psolv%opr(6,i,j,k)=fs%divp_z(1,i,j,k)*sum(fs%itpr_z(:,i,j,k+1)*epsf(i,j,k:k+1))*fs%divw_z( 0,i,j,k+1)
                  fs%psolv%opr(7,i,j,k)=fs%divp_z(0,i,j,k)*sum(fs%itpr_z(:,i,j,k  )*epsf(i,j,k-1:k))*fs%divw_z(-1,i,j,k  )
                  ! Scale it by the cell volume
                  fs%psolv%opr(:,i,j,k)=-fs%psolv%opr(:,i,j,k)*fs%cfg%vol(i,j,k)
               end do
            end do
         end do
         call fs%psolv%setup()
         ! Also allocate the temperature arrays
         allocate(Tprod(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Tprodold(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block prepare_eps
      
      
      ! Create a scalar solver
      create_scalar: block
         use ils_class,      only: gmres
         use vdscalar_class, only: dirichlet,quick
         real(WP) :: diffusivity
         ! Check if we want to model wall losses
         call param_read(tag='Wall temperature',val=Twall,default=-1.0_WP)
         wall_losses=.false.; if (Twall.gt.0.0_WP) wall_losses=.true.
         if (restarted) call df%pullval(name='Twall',val=Twall)
         ! Create scalar solver
         sc=vdscalar(cfg=cfg,scheme=quick,name='Temperature')
         ! Define boundary conditions
         !call sc%add_bcond(name= 'left inflow',type=dirichlet,locator= left_of_tube   )
         !call sc%add_bcond(name='right inflow',type=dirichlet,locator=right_of_tube_sc)
         call sc%add_bcond(name='pipe inflow',type=dirichlet,locator=sctube)
         if (wall_losses) call sc%add_bcond(name='wall',type=dirichlet,locator=wall_locator)
         ! Assign constant diffusivity
         call param_read('Dynamic diffusivity',diffusivity)
         sc%diff=diffusivity
         ! Configure implicit scalar solver
         sc%implicit%maxit=fs%implicit%maxit; sc%implicit%rcvg=fs%implicit%rcvg
         ! Setup the solver
         call sc%setup(implicit_ils=gmres)
      end block create_scalar
      
      
      ! Initialize our temperature field
      initialize_scalar: block
         use vdscalar_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Read in the temperature and pressure
         call param_read('Initial temperature',Tinit)
         call param_read('Inlet temperature',Tinlet)
         call param_read('Initial pressure',pressure)
         ! Handle restart
         if (restarted) then
            call df%pullval(name='pressure',val=pressure )
            call df%pullvar(name='SC'      ,var=sc%SC    )
            call df%pullvar(name='rho'     ,var=sc%rho   )
            call df%pullvar(name='rhoold'  ,var=sc%rhoold)
            call sc%cfg%integrate(sc%rho,integral=fluid_mass)
            call df%pullvar(name='Tprod'   ,var=Tprod    )
         else
            ! Uniform initial temperature and density
            sc%SC=Tinit
            Tprod=Tinit
            sc%rho=1.0_WP; where (sc%cfg%VF.gt.0.0_WP) sc%rho=pressure*Wmlr/(Rcst*Tinit)
            call sc%cfg%integrate(sc%rho,integral=fluid_mass)
            ! Account for porosity here
            sc%rho=sc%rho*epsf
            sc%rhoold=sc%rho
         end if
         ! Apply Dirichlet at the tube
         !call sc%get_bcond( 'left inflow',mybc)
         !do n=1,mybc%itr%no_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   sc%SC (i,j,k)=Tinlet
         !   sc%rho(i,j,k)=pressure*Wmlr/(Rcst*Tinlet)
         !end do
         !call sc%get_bcond('right inflow',mybc)
         !do n=1,mybc%itr%no_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   sc%SC(i,j,k)=Tinlet
         !   sc%rho(i,j,k)=pressure*Wmlr/(Rcst*Tinlet)
         !end do
         call sc%get_bcond('pipe inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc%SC(i,j,k)=Tinlet
            sc%rho(i,j,k)=pressure*Wmlr/(Rcst*Tinlet)
         end do
         if (wall_losses) then
            call sc%get_bcond('wall',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               sc%SC(i,j,k)=Twall
            end do
         end if
         ! Build rhoSC
         call sc%rho_multiply()
         ! Apply all other boundary conditions
         call sc%apply_bcond(time%t,time%dt)
         ! Compute fluid volume - this is the integral of epsf
         call sc%cfg%integrate(epsf,integral=Vtotal)
         ! Recompute pressure
         call sc%cfg%integrate(sc%rhoSC,integral=pressure); pressure=pressure*Rcst/(Wmlr*Vtotal)
      end block initialize_scalar
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools,     only: Pi
         use lowmach_class, only: bcond
         use parallel,      only: MPI_REAL_WP
         use mpi_f08,       only: MPI_ALLREDUCE,MPI_SUM
         type(bcond), pointer :: mybc
         integer :: n,i,j,k,ierr
         real(WP) :: myAin
         ! Handle restart
         if (restarted) then
            call df%pullvar(name='rhoU',var=fs%rhoU)
            call df%pullvar(name='rhoV',var=fs%rhoV)
            call df%pullvar(name='rhoW',var=fs%rhoW)
            call df%pullvar(name='P'   ,var=fs%P   )
         end if
         ! Read in MFR
         call param_read('Inlet MFR (kg/s)',MFR)
         ! Calculate inflow area
         myAin=0.0_WP
         !call fs%get_bcond('left inflow',mybc)
         !do n=1,mybc%itr%n_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   myAin=myAin+fs%cfg%dy(j)*fs%cfg%dz(k)
         !end do
         !call fs%get_bcond('right inflow',mybc)
         !do n=1,mybc%itr%n_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   myAin=myAin+fs%cfg%dy(j)*fs%cfg%dz(k)
         !end do
         call fs%get_bcond('pipe inflow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAin=myAin+fs%cfg%dx(i)*fs%cfg%dz(k)
         end do
         call MPI_ALLREDUCE(myAin,Ain,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
         ! Form inflow momentum
         rhoUin=MFR/Ain
         ! Apply Dirichlet at the tube
         !call fs%get_bcond('left inflow',mybc)
         !do n=1,mybc%itr%no_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   fs%rhoU(i,j,k)=-rhoUin
         !end do
         !call fs%get_bcond('right inflow',mybc)
         !do n=1,mybc%itr%no_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   fs%rhoU(i,j,k)=+rhoUin
         !end do
         call fs%get_bcond('pipe inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%rhoV(i,j,k)=+rhoUin
         end do
         ! Set density from scalar
         fs%rho=0.5_WP*(sc%rho+sc%rhoold)
         ! Handle restart
         if (restarted) call df%pullvar(name='rho_fs',var=fs%rho)
         ! Form momentum
         call fs%rho_divide
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resSC)
         call fs%get_div(drhodt=resSC)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         allocate(viscmol(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(diffmol(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pvessel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',sc%rho)
         call ens_out%add_scalar('temperature',sc%SC)
         call ens_out%add_scalar('visc',fs%visc)
         call ens_out%add_scalar('epsp',epsp)
         call ens_out%add_scalar('Tprod',Tprod)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         Tavg=sc%rhoSCint/sc%rhoint
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
         call mfile%add_column(sc%SCmax,'Tmax')
         call mfile%add_column(sc%SCmin,'Tmin')
         call mfile%add_column(sc%rhomax,'RHOmax')
         call mfile%add_column(sc%rhomin,'RHOmin')
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
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(sc%SCint,'SC integral')
         call consfile%add_column(Tavg,'Tavg')
         call consfile%add_column(sc%rhoint,'Mass')
         call consfile%add_column(pressure,'Pthermo')
         call consfile%add_column(Twall,'Twall')
         call consfile%write()
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
         
         ! Remember fluid mass, pressure, and wall temperature
         fluid_mass_old=fluid_mass
         pressure_old=pressure
         Twallold=Twall
         Tprodold=Tprod
         
         ! Remember old scalar
         sc%rhoold=sc%rho
         sc%SCold=sc%SC; sc%rhoSCold=sc%rhoSC
         
         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! ============ UPDATE PROPERTIES ====================
         call get_visc(); viscmol=fs%visc
         call get_cond(); diffmol=sc%diff
         
         ! Turbulence modeling
         call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         call sgs%get_visc(dt=time%dtold,rho=fs%rho,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
         where (sgs%visc.lt.-fs%visc)
            sgs%visc=-fs%visc
         end where
         fs%visc=fs%visc+sgs%visc
         sc%diff=sc%diff+sgs%visc/Pr_turb    !< This is an assumption of constant turbulent Pr number
         ! ===================================================
         
         
         ! ============ Estimate heat flux at wall ============
         heat_transfer: block
            use vdscalar_class, only: bcond
            type(bcond), pointer :: mybc
            integer :: i,j,k,n
            ! Update wall temperature based on empirical model
            Twall=Twallold+time%dt*(Tavg-Twallold)/tau_wall_in+time%dt*(Tout-Twallold)/tau_wall_out
            ! Update the boundary condition
            if (wall_losses) then
               call sc%get_bcond('wall',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k)=Twall
               end do
            end if
         end block heat_transfer
         ! ====================================================
         
         
         ! Account for heat transfer to product
         !product_heating: block
         !   integer :: i,j,k
         !   do k=fs%cfg%kmin_,fs%cfg%kmax_
         !      do j=fs%cfg%jmin_,fs%cfg%jmax_
         !         do i=fs%cfg%imin_,fs%cfg%imax_
         !            if (sc%cfg%VF(i,j,k).gt.0.0_WP) Tprod(i,j,k)=Tprod(i,j,k)+Cp*diffmol(i,j,k)*epsp(i,j,k)*(sc%SC(i,j,k)-Tprodold(i,j,k))/(rprod*bporo*bporo*Tperm*Cprod)
         !         end do
         !      end do
         !   end do
         !   call cfg%sync(Tprod)
         !end block product_heating
         
         
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            
            ! ============ VELOCITY SOLVER ======================
            ! Build n+1 density
            fs%rho=0.5_WP*(sc%rho+sc%rhoold)
            
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS - without grad(p) for now
            resSC=fs%P; fs%P=0.0_WP
            call fs%get_dmomdt(resU,resV,resW)
            fs%P=resSC
            
            ! Add epsf*grad(P)
            modified_pgrad_pred: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)-sum(fs%divu_x(:,i,j,k)*fs%P(i-1:i,j,k))*sum(fs%itpr_x(:,i,j,k)*epsf(i-1:i,j,k))
                        resV(i,j,k)=resV(i,j,k)-sum(fs%divv_y(:,i,j,k)*fs%P(i,j-1:j,k))*sum(fs%itpr_y(:,i,j,k)*epsf(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)-sum(fs%divw_z(:,i,j,k)*fs%P(i,j,k-1:k))*sum(fs%itpr_z(:,i,j,k)*epsf(i,j,k-1:k))
                     end do
                  end do
               end do
            end block modified_pgrad_pred
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Also add permeability
            permeability: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)-sum(fs%itpr_x(:,i,j,k)*epsf(i-1:i,j,k)*viscmol(i-1:i,j,k)*epsp(i-1:i,j,k))*fs%U(i,j,k)/(bporo*bperm)
                        resV(i,j,k)=resV(i,j,k)-sum(fs%itpr_y(:,i,j,k)*epsf(i,j-1:j,k)*viscmol(i,j-1:j,k)*epsp(i,j-1:j,k))*fs%V(i,j,k)/(bporo*bperm)
                        resW(i,j,k)=resW(i,j,k)-sum(fs%itpr_z(:,i,j,k)*epsf(i,j,k-1:k)*viscmol(i,j,k-1:k)*epsp(i,j,k-1:k))*fs%W(i,j,k)/(bporo*bperm)
                     end do
                  end do
               end do
            end block permeability
            
            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Get momentum
            call fs%rho_multiply()
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%tmid,time%dtmid)
            mom_bcond: block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               !call fs%get_bcond('left inflow',mybc)
               !do n=1,mybc%itr%no_
               !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               !   fs%rhoU(i,j,k)=-rhoUin
               !end do
               !call fs%get_bcond('right inflow',mybc)
               !do n=1,mybc%itr%no_
               !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               !   fs%rhoU(i,j,k)=+rhoUin
               !end do
               call fs%get_bcond('pipe inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%rhoV(i,j,k)=+rhoUin
               end do
            end block mom_bcond
            
            ! Solve Poisson equation
            call fs%correct_mfr()                          !< No outlet so this gets the MFR imbalance
            fluid_mass=fluid_mass_old-sum(fs%mfr)*time%dt  !< Update mass in system
            call get_rho(mass=fluid_mass)                  !< Adjust rho and pressure accordingly
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct momentum and rebuild velocity - pressure gradient is modified here due to porosity
            fs%P=fs%P+fs%psolv%sol
            modified_pgrad_corr: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%rhoU(i,j,k)=fs%rhoU(i,j,k)-time%dtmid*sum(fs%divu_x(:,i,j,k)*fs%psolv%sol(i-1:i,j,k))*sum(fs%itpr_x(:,i,j,k)*epsf(i-1:i,j,k))
                        fs%rhoV(i,j,k)=fs%rhoV(i,j,k)-time%dtmid*sum(fs%divv_y(:,i,j,k)*fs%psolv%sol(i,j-1:j,k))*sum(fs%itpr_y(:,i,j,k)*epsf(i,j-1:j,k))
                        fs%rhoW(i,j,k)=fs%rhoW(i,j,k)-time%dtmid*sum(fs%divw_z(:,i,j,k)*fs%psolv%sol(i,j,k-1:k))*sum(fs%itpr_z(:,i,j,k)*epsf(i,j,k-1:k))
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%rhoU)
               call fs%cfg%sync(fs%rhoV)
               call fs%cfg%sync(fs%rhoW)
            end block modified_pgrad_corr
            call fs%rho_divide
            
            ! Update div for monitoring
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            ! ===================================================
            
            
            ! ============= SCALAR SOLVER =======================
            ! Build mid-time scalar
            sc%SC=0.5_WP*(sc%SC+sc%SCold)
            
            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Add pressure term
            where (sc%cfg%VF.gt.0.0_WP) resSC=resSC+epsf/Cp*(pressure-pressure_old)/time%dt
            
            ! Account for heat transfer to product
            !heat_to_product: block
            !   integer :: i,j,k
            !   do k=fs%cfg%kmin_,fs%cfg%kmax_
            !      do j=fs%cfg%jmin_,fs%cfg%jmax_
            !         do i=fs%cfg%imin_,fs%cfg%imax_
            !            if (sc%cfg%VF(i,j,k).gt.0.0_WP) resSC(i,j,k)=resSC(i,j,k)-diffmol(i,j,k)*epsp(i,j,k)*(sc%SC(i,j,k)-Tprod(i,j,k))/(bporo*Tperm)
            !         end do
            !      end do
            !   end do
            !end block heat_to_product
            
            ! Temporarily redefine density to include sled/basket assembly
            sc%rho   =sc%rho   +Msteel*Cs/(Vtotal*Cp)
            sc%rhoold=sc%rhoold+Msteel*Cs/(Vtotal*Cp)
            
            ! Temporarily redefine density to include product as well
            sc%rho   =sc%rho   +Mprod*Cprod/(Vtotal*Cp)
            sc%rhoold=sc%rhoold+Mprod*Cprod/(Vtotal*Cp)
            
            ! Assemble explicit residual
            where (sc%cfg%VF.gt.0.0_WP) resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
            
            ! Form implicit residual - rho is modified to make volumetric heat losses implicit
            call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Set densities back to fluid only
            sc%rho   =sc%rho   -Msteel*Cs/(Vtotal*Cp)
            sc%rhoold=sc%rhoold-Msteel*Cs/(Vtotal*Cp)
            
            ! Set densities back to fluid only
            sc%rho   =sc%rho   -Mprod*Cprod/(Vtotal*Cp)
            sc%rhoold=sc%rhoold-Mprod*Cprod/(Vtotal*Cp)
            
            ! Apply this residual
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC
            
            ! Update density and pressure
            call get_rho(mass=fluid_mass)
            
            ! Multiply by density
            call sc%rho_multiply()
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         
         ! Recompute interpolated velocity
         call fs%interp_vel(Ui,Vi,Wi)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         Tavg=sc%rhoSCint/sc%rhoint
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
         
         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
               character(len=str_medium) :: dirname,timestamp
               ! Prefix for files
               dirname='restart_'; write(timestamp,'(es12.5)') time%t
               ! Prepare a new directory
               if (fs%cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
               ! Populate df and write it
               call df%pushval(name='t'       ,val=time%t   )
               call df%pushval(name='dt'      ,val=time%dt  )
               call df%pushval(name='pressure',val=pressure )
               call df%pushval(name='Twall'   ,val=Twall    )
               call df%pushvar(name='rhoU'    ,var=fs%rhoU  )
               call df%pushvar(name='rhoV'    ,var=fs%rhoV  )
               call df%pushvar(name='rhoW'    ,var=fs%rhoW  )
               call df%pushvar(name='P'       ,var=fs%P     )
               call df%pushvar(name='rho_fs'  ,var=fs%rho   )
               call df%pushvar(name='rho'     ,var=sc%rho   )
               call df%pushvar(name='rhoold'  ,var=sc%rhoold)
               call df%pushvar(name='SC'      ,var=sc%SC    )
               call df%pushvar(name='Tprod'   ,var=Tprod    )
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
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
