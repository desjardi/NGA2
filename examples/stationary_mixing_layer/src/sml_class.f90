!> Definition for an sml class
module sml_class
   use precision,         only: WP
   use config_class,      only: config
   use fft2d_class,       only: fft2d
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use monitor_class,     only: monitor
   use ensight_class,     only: ensight
   use event_class,       only: event
   implicit none
   private
   
   public :: sml
   
   !> SML object
   type :: sml
      !> Config
      type(config)      :: cfg   !< Mesh for solver
      !> Flow solver
      type(incomp)      :: fs    !< Incompressible flow solver
      type(fft2d)       :: ps    !< FFT-based linear solver
      type(timetracker) :: time  !< Time info
      !> Ensight postprocessing
      type(ensight)     :: ens_out
      type(event)       :: ens_evt
      !> Simulation monitor file
      type(monitor)     :: mfile,cflfile
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
      !> Fluid parameters
      real(WP) :: visc
      !> Mean velocity profile
      real(WP), dimension(:), allocatable :: meanU
   contains
      procedure :: init          !< Initialize SML simulation
      procedure :: step          !< Advance SML simulation by one time step
      procedure :: final         !< Finalize SML simulation
      procedure :: force         !< Force fluctuating flow using imposed mean flow
   end type sml
   
   
contains
   
   
   !> Initialization of SML simulation
   subroutine init(this)
      implicit none
      class(sml), intent(inout) :: this
      

      ! Create the SML mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1; x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx; end do
         do j=1,ny+1; y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly; end do
         do k=1,nz+1; z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz; end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='sml')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      

      ! Initialize the work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
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
      
      
      ! Create an incompressible flow solver with slip conditions top and bottom
      create_flow_solver: block
         use param, only: param_read
         use incomp_class, only: slip
         integer :: i,j,k
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='mixing_layer')
         ! Set density to 1.0
         this%fs%rho=1.0_WP
         ! Set viscosity from Reynolds number
         call param_read('Reynolds number',this%visc); this%visc=1.0_WP/this%visc
         this%fs%visc=this%visc
         ! Add slip conditions top and bottom
         call this%fs%add_bcond(name='ymslip',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         call this%fs%add_bcond(name='ypslip',type=slip,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         ! Create pressure solver
         this%ps=fft2d(cfg=this%cfg,name='Pressure',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)
      end block create_flow_solver
      
      
      ! Prepare random initial fluctuating velocity field
      initialize_velocity: block
         use random, only: random_normal
         integer :: i,j,k
         real(WP) :: initial_rms,Umean,Vmean,Wmean
         ! Gaussian initial field
         initial_rms=0.1_WP
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  this%fs%U(i,j,k)=random_normal(m=0.0_WP,sd=initial_rms)*exp(-this%cfg%ym(j)**2)
                  this%fs%V(i,j,k)=random_normal(m=0.0_WP,sd=initial_rms)*exp(-this%cfg%y (j)**2)
                  this%fs%W(i,j,k)=random_normal(m=0.0_WP,sd=initial_rms)*exp(-this%cfg%ym(j)**2)
               end do
            end do
         end do
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%V)
         call this%fs%cfg%sync(this%fs%W)
         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call this%fs%cfg%integrate(A=this%fs%U,integral=Umean); Umean=Umean/this%fs%cfg%vol_total; this%fs%U=this%fs%U-Umean
         call this%fs%cfg%integrate(A=this%fs%V,integral=Vmean); Vmean=Vmean/this%fs%cfg%vol_total; this%fs%V=this%fs%V-Vmean
         call this%fs%cfg%integrate(A=this%fs%W,integral=Wmean); Wmean=Wmean/this%fs%cfg%vol_total; this%fs%W=this%fs%W-Wmean
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
         ! Calculate cell-centered velocities and divergence
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         call this%fs%get_div()
      end block initialize_velocity
      
      ! Prepare mean velocity profile for forcing
      initialize_forcing: block
         integer :: j
         ! Allocate meanU
         allocate(this%meanU(this%cfg%jmino_:this%cfg%jmaxo_))
         ! Initialize to tanh profile
         do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
            this%meanU(j)=0.5_WP*tanh(0.5_WP*this%cfg%ym(j))
         end do
      end block initialize_forcing
      
      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='mixing_layer')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight

      
      ! Create monitoring file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation')
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
   

   !> Take one time step with specified dt
   subroutine step(this)
      implicit none
      class(sml), intent(inout) :: this
      
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
         this%resU=-2.0_WP*(this%fs%U-this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%V-this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%W-this%fs%Wold)+this%time%dt*this%resW
         
         ! Force using -div(u' meanU)-div(meanU u')
         call this%force()
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
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
      class(sml), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
   !> Function that localizes the top (y+) of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator
   
   
   !> Function that localizes the bottom (y-) of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator
   

   !> Subroutine that forces fluctuating velocity equation using imposed mean velocity
   subroutine force(this)
      implicit none
      class(sml), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Force u' by -div(meanU u')-div(u' meanU)
      do kk=this%cfg%kmin_,this%cfg%kmax_+1; do jj=this%cfg%jmin_,this%cfg%jmax_+1; do ii=this%cfg%imin_,this%cfg%imax_+1
         ! Fluxes on x-face
         i=ii-1; j=jj-1; k=kk-1; FX(i,j,k)=-this%fs%rho*this%meanU(j)*sum(this%fs%itpu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))&
         &                                 -this%fs%rho*sum(this%fs%itpu_x(:,i,j,k)*this%fs%U(i:i+1,j,k))*this%meanU(j)
         i=ii  ; j=jj  ; k=kk  ; FY(i,j,k)=-this%fs%rho*sum(this%fs%itpu_y(:,i,j,k)*this%meanU(j-1:j))*sum(this%fs%itpv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))
         i=ii  ; j=jj  ; k=kk  ; FZ(i,j,k)=-this%fs%rho*this%meanU(j)*sum(this%fs%itpw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))
      end do; end do; end do
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         this%resU(i,j,k)=this%resU(i,j,k)+this%time%dt*(sum(this%fs%divu_x(:,i,j,k)*FX(i-1:i,j,k))+sum(this%fs%divu_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%fs%divu_z(:,i,j,k)*FZ(i,j,k:k+1)))
      end do; end do; end do
      ! Force v' by -div(meanU u')-div(u' meanU)
      do kk=this%cfg%kmin_,this%cfg%kmax_+1; do jj=this%cfg%jmin_,this%cfg%jmax_+1; do ii=this%cfg%imin_,this%cfg%imax_+1
         i=ii  ; j=jj  ; k=kk  ; FX(i,j,k)=-this%fs%rho*sum(this%fs%itpv_x(:,i,j,k)*this%fs%V(i-1:i,j,k))*sum(this%fs%itpu_y(:,i,j,k)*this%meanU(j-1:j))
      end do; end do; end do
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         this%resV(i,j,k)=this%resV(i,j,k)+this%time%dt*sum(this%fs%divv_x(:,i,j,k)*FX(i:i+1,j,k))
      end do; end do; end do
      ! Force w' by -div(meanU u')-div(u' meanU)
      do kk=this%cfg%kmin_,this%cfg%kmax_+1; do jj=this%cfg%jmin_,this%cfg%jmax_+1; do ii=this%cfg%imin_,this%cfg%imax_+1
         i=ii  ; j=jj  ; k=kk  ; FX(i,j,k)=-this%fs%rho*sum(this%fs%itpw_x(:,i,j,k)*this%fs%W(i-1:i,j,k))*this%meanU(j)
      end do; end do; end do
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         this%resW(i,j,k)=this%resW(i,j,k)+this%time%dt*sum(this%fs%divw_x(:,i,j,k)*FX(i:i+1,j,k))
      end do; end do; end do
      
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      
   end subroutine force
   

end module sml_class