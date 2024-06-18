!> Definition for an hit class
module hit_class
   use precision,         only: WP
   use config_class,      only: config
   use fft3d_class,       only: fft3d
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: hit
   
   !> HIT object
   type :: hit
      !> Config
      type(config)      :: cfg   !< Mesh for solver
      !> Flow solver
      type(incomp)      :: fs    !< Incompressible flow solver
      type(fft3d)       :: ps    !< FFT-based linear solver
      type(timetracker) :: time  !< Time info
      !> Simulation monitor file
      type(monitor) :: mfile     !< General simulation monitoring
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:,:), allocatable :: SR                !< Strain rate tensor
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      !> Turbulence parameters
      real(WP) :: ti                            ! Turbulence intensity
      real(WP) :: visc,meanU,meanV,meanW
      real(WP) :: Urms_tgt,tke_tgt,eps_tgt      ! u',k, and dissipation rate
      real(WP) :: tko_tgt,eta_tgt               ! Kolmogorov time and length scales
      real(WP) :: Rel_tgt,Ret_tgt               ! Lambda and turbulent Reynolds numbers
      real(WP) :: tau_tgt                       ! Eddy turnover time
      real(WP) :: Urms,tke,eps,Ret,Rel,eta,ell  ! Current turbulence parameters (ell is large eddy size)
      !> Forcing constant
      real(WP) :: forcing
   contains
      procedure, private :: compute_stats          !< Turbulence information
      procedure :: init                            !< Initialize HIT simulation
      procedure :: step                            !< Advance HIT simulation by one time step
      procedure :: final                           !< Finalize HIT simulation
   end type hit
   
   
contains
   
   
   !> Compute turbulence stats (assumes rho=1)
   subroutine compute_stats(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      class(hit), intent(inout) :: this
      real(WP) :: myTKE,myEPS
      integer :: i,j,k,ierr
      ! Compute mean velocities
      call this%fs%cfg%integrate(A=this%fs%U,integral=this%meanU); this%meanU=this%meanU/this%fs%cfg%vol_total
      call this%fs%cfg%integrate(A=this%fs%V,integral=this%meanV); this%meanV=this%meanV/this%fs%cfg%vol_total
      call this%fs%cfg%integrate(A=this%fs%W,integral=this%meanW); this%meanW=this%meanW/this%fs%cfg%vol_total
      ! Compute strainrate and grad(U)
      call this%fs%get_strainrate(SR=this%SR)
      call this%fs%get_gradu(gradu=this%gradU)
      ! Compute current TKE and dissipation rate
      myTKE=0.0_WP
      myEPS=0.0_WP
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               myTKE=myTKE+0.5_WP*((this%fs%U(i,j,k)-this%meanU)**2+(this%fs%V(i,j,k)-this%meanV)**2+(this%fs%W(i,j,k)-this%meanW)**2)*this%fs%cfg%vol(i,j,k)
               myEPS=myEPS+2.0_WP*this%fs%cfg%vol(i,j,k)*(this%SR(1,i,j,k)**2+this%SR(2,i,j,k)**2+this%SR(3,i,j,k)**2+2.0_WP*(this%SR(4,i,j,k)**2+this%SR(5,i,j,k)**2+this%SR(6,i,j,k)**2))
            end do
         end do
      end do
      call MPI_ALLREDUCE(myTKE,this%tke,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%tke=this%tke/this%fs%cfg%vol_total
      call MPI_ALLREDUCE(myEPS,this%eps,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%eps=this%eps*this%visc/this%fs%cfg%vol_total
      ! Compute standard parameters for HIT
      this%Urms=sqrt(2.0_WP/3.0_WP*this%tke)
      this%Ret=this%tke**2.0_WP/(this%visc*this%eps)
      this%Rel=sqrt(20.0_WP*this%Ret/3.0_WP)
      this%eta=((this%visc)**3.0_WP/this%eps)**0.25_WP
      this%ell=(2.0_WP*this%tke/3.0_WP)**1.5_WP/this%eps
   end subroutine compute_stats
   
   
   !> Initialization of HIT simulation
   subroutine init(this,group,xend,visc)
      use mpi_f08, only: MPI_Group
      implicit none
      class(hit), intent(inout)   :: this
      type(MPI_Group), intent(in) :: group
      real(WP) :: xend,visc
      
      ! Create the HIT mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         real(WP), dimension(:), allocatable :: x,y
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: j,ny,nz
         real(WP) :: Ly,Lz
         ! Read in grid definition
         call param_read('Ly',Ly); call param_read('Lz',Lz); Ly=max(Ly,Lz)
         call param_read('ny',ny); call param_read('nz',nz); ny=max(ny,nz)
         allocate(y(ny+1)); allocate(x(ny+1))
         ! Create simple rectilinear grid in y and z
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         ! Same grid in x, but shifted so it ends at xend
         x=y-y(ny+1)+xend
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=y,xper=.true.,yper=.true.,zper=.true.,name='HIT')
         ! Read in partition
         call param_read('Partition',partition,short='p'); partition(1)=1
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config

      ! Initialize the work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SR  (1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Create a single-phase periodic flow solver
      create_flow_solver: block
         use mathtools, only: Pi
         use param,     only: param_read
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='NS solver')
         ! Set density to 1.0
         this%fs%rho=1.0_WP
         ! Set viscosity from Reynolds number
         this%visc=visc
         this%fs%visc=this%visc
         ! Prepare and configure pressure solver
         this%ps=fft3d(cfg=this%cfg,name='Pressure',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)
      end block create_flow_solver
      
      
      ! Prepare initial velocity field
      initialize_velocity: block
         use random,    only: random_normal
         use mathtools, only: Pi
         use param,     only: param_read,param_exists
         use messager,  only: log
         use string,    only: str_long
         character(str_long) :: message
         real(WP) :: max_forcing_estimate
         integer :: i,j,k
         ! Read in turbulence intensity for turbulence injection
         call param_read('Turbulence intensity',this%ti)
         ! Read in target Re_lambda and convert to target Urms
         call param_read('Target Re_lambda',this%Urms_tgt)
         this%Urms_tgt=this%visc/(3.0_WP*this%cfg%xL)*this%Urms_tgt**2
         ! Calculate other target quantities assuming l=0.2*xL
         this%tke_tgt=1.5_WP*this%Urms_tgt**2
         this%eps_tgt=5.0_WP*this%Urms_tgt**3/this%cfg%xL
         this%tko_tgt=sqrt(this%visc/this%eps_tgt)
         this%eta_tgt=(this%visc**3/this%eps_tgt)**(0.25_WP)
         this%Rel_tgt=sqrt(3.0_WP*this%Urms_tgt*this%cfg%xL/this%visc)
         this%Ret_tgt=this%tke_tgt**2/(this%eps_tgt*this%visc)
         this%tau_tgt=2.0_WP*this%tke_tgt/(3.0_WP*this%eps_tgt)
         ! Read in forcing parameter (we need dt<tau_tgt/forcing)
         max_forcing_estimate=3.0_WP*this%tau_tgt*this%Urms_tgt/this%cfg%min_meshsize
         call param_read('Forcing constant',this%forcing,default=max_forcing_estimate)
         ! Output nominal turbulence parameters
         if (this%cfg%amRoot) then
            write(message,'("Target turbulence parameters:")'); call log(message)
            write(message,'("[HIT setup] => Urms              =",es12.5)') this%Urms_tgt; call log(message)
            write(message,'("[HIT setup] => Re_lambda         =",es12.5)') this%Rel_tgt; call log(message)
            write(message,'("[HIT setup] => Re_turb           =",es12.5)') this%Ret_tgt; call log(message)
            write(message,'("[HIT setup] => Kolmogorov Lscale =",es12.5)') this%eta_tgt; call log(message)
            write(message,'("[HIT setup] => Kolmogorov Tscale =",es12.5)') this%tko_tgt; call log(message)
            write(message,'("[HIT setup] => Epsilon           =",es12.5)') this%eps_tgt; call log(message)
            write(message,'("[HIT setup] => Eddyturnover time =",es12.5)') this%tau_tgt; call log(message)
         end if
         ! Gaussian initial field
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  this%fs%U(i,j,k)=random_normal(m=0.0_WP,sd=this%Urms_tgt)
                  this%fs%V(i,j,k)=random_normal(m=0.0_WP,sd=this%Urms_tgt)
                  this%fs%W(i,j,k)=random_normal(m=0.0_WP,sd=this%Urms_tgt)
               end do
            end do
         end do
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%V)
         call this%fs%cfg%sync(this%fs%W)
         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call this%fs%cfg%integrate(A=this%fs%U,integral=this%meanU); this%meanU=this%meanU/this%fs%cfg%vol_total; this%fs%U=this%fs%U-this%meanU
         call this%fs%cfg%integrate(A=this%fs%V,integral=this%meanV); this%meanV=this%meanV/this%fs%cfg%vol_total; this%fs%V=this%fs%V-this%meanV
         call this%fs%cfg%integrate(A=this%fs%W,integral=this%meanW); this%meanW=this%meanW/this%fs%cfg%vol_total; this%fs%W=this%fs%W-this%meanW
         ! Project to ensure divergence-free
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%resU
         this%fs%V=this%fs%V-this%resV
         this%fs%W=this%fs%W-this%resW
         ! Calculate divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create monitoring file
      create_monitor: block
         ! Prepare some info about turbulence
         call this%fs%get_max()
         call this%compute_stats()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'hit')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%Ret,'Re_turb')
         call this%mfile%add_column(this%Rel,'Re_lambda')
         call this%mfile%add_column(this%Urms,'Urms')
         call this%mfile%add_column(this%TKE,'TKE')
         call this%mfile%add_column(this%EPS,'Epsilon')
         call this%mfile%add_column(this%ell,'Large eddy size')
         call this%mfile%add_column(this%eta,'Kolmogorov length')
         call this%mfile%write()
      end block create_monitor
      
      
   end subroutine init
   

   !> Take one time step with specified dt
   subroutine step(this,dt)
      implicit none
      class(hit), intent(inout) :: this
      real(WP), intent(in) :: dt
      
      ! Increment time based on provided dt
      this%time%dt=dt; call this%time%increment()
      
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
         this%resU=-2.0_WP*(this%fs%U-this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%V-this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%W-this%fs%Wold)+this%time%dt*this%resW
         
         ! Apply HIT forcing
         hit_forcing: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
            use parallel, only: MPI_REAL_WP
            real(WP) :: myTKE,A,myEPSp,EPSp
            integer :: i,j,k,ierr
            ! Calculate mean velocity
            call this%fs%cfg%integrate(A=this%fs%U,integral=this%meanU); this%meanU=this%meanU/this%fs%cfg%vol_total
            call this%fs%cfg%integrate(A=this%fs%V,integral=this%meanV); this%meanV=this%meanV/this%fs%cfg%vol_total
            call this%fs%cfg%integrate(A=this%fs%W,integral=this%meanW); this%meanW=this%meanW/this%fs%cfg%vol_total
            ! Calculate TKE and pseudo-EPS
            call this%fs%get_gradu(gradu=this%gradU)
            myTKE=0.0_WP; myEPSp=0.0_WP
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     myTKE =myTKE +0.5_WP*((this%fs%U(i,j,k)-this%meanU)**2+(this%fs%V(i,j,k)-this%meanV)**2+(this%fs%W(i,j,k)-this%meanW)**2)*this%fs%cfg%vol(i,j,k)
                     myEPSp=myEPSp+this%fs%cfg%vol(i,j,k)*(this%gradU(1,1,i,j,k)**2+this%gradU(1,2,i,j,k)**2+this%gradU(1,3,i,j,k)**2+&
                     &                                     this%gradU(2,1,i,j,k)**2+this%gradU(2,2,i,j,k)**2+this%gradU(2,3,i,j,k)**2+&
                     &                                     this%gradU(3,1,i,j,k)**2+this%gradU(3,2,i,j,k)**2+this%gradU(3,3,i,j,k)**2)
                  end do
               end do
            end do
            call MPI_ALLREDUCE(myTKE,this%tke,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%tke=this%tke/this%fs%cfg%vol_total
            call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); EPSp=EPSp*this%visc/this%fs%cfg%vol_total
            A=(EPSp-this%forcing*(this%tke-this%tke_tgt)/this%tau_tgt)/(2.0_WP*this%tke)
            this%resU=this%resU+A*this%time%dt*(this%fs%U-this%meanU)
            this%resV=this%resV+A*this%time%dt*(this%fs%V-this%meanV)
            this%resW=this%resW+A*this%time%dt*(this%fs%W-this%meanW)
         end block hit_forcing

         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Solve Poisson equation
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU
         this%fs%V=this%fs%V-this%time%dt*this%resV
         this%fs%W=this%fs%W-this%time%dt*this%resW
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute divergence
      call this%fs%get_div()
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%compute_stats()
      call this%mfile%write()
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(hit), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%gradU,this%SR)
      
   end subroutine final
   

end module hit_class