!> Definition for a channel class
module channel_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use ensight_class,     only: ensight
   use fft2d_class,       only: fft2d
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: channel
   
   !> Channel object
   type :: channel
      
      !> Input file
      type(inputfile) :: input

      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs    !< Incompressible flow solver
      type(fft2d)       :: ps    !< Pressure solver based on 2D FFT
      type(timetracker) :: time  !< Time info
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      
   contains
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type channel
   
   
contains
   
   
   !> Initialization of channel simulation
   subroutine init(this)
      implicit none
      class(channel), intent(inout) :: this
      
      ! Read in the input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='input.fluid')
      end block read_input
      

      ! Create the channel mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1))
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='channel')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create a multiphase flow solver with bconds
      create_flow_solver: block
         use incomp_class, only: dirichlet,clipped_neumann
         use incomp_class, only: bcond
         use random,       only: random_normal
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         real(WP) :: visc,Uin
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call this%input%read('Density',this%fs%rho)
         call this%input%read('Dynamic viscosity',visc); this%fs%visc=visc
         ! Define boundary conditions
         call this%fs%add_bcond(name='inflow', type=dirichlet      ,locator=xm_locator,face='x',dir=-1,canCorrect=.false.)
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=xp_locator,face='x',dir=+1,canCorrect=.true. )
         ! Configure pressure solver
         this%ps=fft2d(cfg=this%cfg,name='Pressure',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)
         ! Read inflow velocity
         call this%input%read('Inlet velocity',Uin)
         ! Random initial field
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  this%fs%U(i,j,k)=Uin+random_normal(m=0.0_WP,sd=0.01_WP*Uin)
                  this%fs%V(i,j,k)=    random_normal(m=0.0_WP,sd=0.01_WP*Uin)
                  this%fs%W(i,j,k)=    random_normal(m=0.0_WP,sd=0.01_WP*Uin)
               end do
            end do
         end do
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%V)
         call this%fs%cfg%sync(this%fs%W)
         ! Apply inflow condition
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=Uin
         end do
         ! Compute MFR through all boundary conditions
         call this%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block create_flow_solver
      

      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='fluid')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'fluid')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this,ls)
      use lss_class, only: lss
      implicit none
      class(channel), intent(inout) :: this
      class(lss), intent(inout) :: ls
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%rho*this%fs%U-this%fs%rho*this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%rho*this%fs%V-this%fs%rho*this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%rho*this%fs%W-this%fs%rho*this%fs%Wold)+this%time%dt*this%resW   
         
         ! Form implicit residuals
         !call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho
         
         ! Add momentum source term from direct forcing
         ibm_correction: block
            integer :: i,j,k
            this%resU=this%fs%rho
            call ls%get_source(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU)
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     this%fs%U(i,j,k)=this%fs%U(i,j,k)+ls%srcU(i,j,k)
                     this%fs%V(i,j,k)=this%fs%V(i,j,k)+ls%srcV(i,j,k)
                     this%fs%W(i,j,k)=this%fs%W(i,j,k)+ls%srcW(i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibm_correction
         
         ! Apply other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Output to ensight
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
         
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(channel), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   
end module channel_class