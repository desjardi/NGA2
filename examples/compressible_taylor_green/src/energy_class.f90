!> Internal energy solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes variable diffusivity and density.
module energy_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use linsol_class,   only: linsol
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: energy,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2                !< Dirichlet condition
   integer, parameter, public :: neumann=3                  !< Zero normal gradient
   
   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                          !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'     !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                       !< Bcond type
      integer :: dir                                        !< Bcond direction (1 to 6)
      type(iterator) :: itr                                 !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Internal energy solver object definition
   type :: energy
      
      ! This is our config
      class(config), pointer :: cfg                         !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_ENERGY'    !< Solver name (default=UNNAMED_ENERGY)
      
      ! Variable property fluid
      real(WP), dimension(:,:,:), allocatable :: diff       !< These is our dynamic diffusivity
      real(WP) :: Cv=1.0_WP
      
      ! Boundary condition list
      integer :: nbc                                        !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                      !< List of bcond for our solver
      
      ! Energy variables
      real(WP), dimension(:,:,:), allocatable :: E          !< Current E array
      real(WP), dimension(:,:,:), allocatable :: Eold       !< Old E array
      
      ! Implicit energy solver
      class(linsol), pointer :: implicit                    !< Iterative linear solver object for an implicit prediction of the energy residual
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itp_x,itp_y,itp_z        !< Interpolation of diffusivity to the face
      real(WP), dimension(:,:,:,:), allocatable :: itpe_x,itpe_y,itpe_z      !< 2nd order centered interpolation for energy
      real(WP), dimension(:,:,:,:), allocatable :: up2_xp,up2_yp,up2_zp      !< 2nd order upwind interpolation for energy (+)
      real(WP), dimension(:,:,:,:), allocatable :: up2_xm,up2_ym,up2_zm      !< 2nd order upwind interpolation for energy (+)
      real(WP), dimension(:,:,:,:), allocatable :: weno_xp,weno_yp,weno_zp   !< WENO interpolation for energy (+)
      real(WP), dimension(:,:,:,:), allocatable :: weno_xm,weno_ym,weno_zm   !< WENO interpolation for energy (-)
      real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z        !< Divergence of energy flux
      real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z        !< Gradient of temperature
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask        !< Integer array used for modifying E metrics
      
      ! Monitoring quantities
      real(WP) :: Emax,Emin,rhoEint                         !< Maximum and minimum, integral energy
      
   contains
      procedure :: print=>energy_print                      !< Output solver to the screen
      procedure :: initialize                               !< Initialize the energy solver
      procedure :: setup                                    !< Finish configuring the energy solver
      procedure :: add_bcond                                !< Add a boundary condition
      procedure :: get_bcond                                !< Get a boundary condition
      procedure :: apply_bcond                              !< Apply all boundary conditions
      procedure :: init_metrics                             !< Initialize metrics
      procedure :: adjust_metrics                           !< Adjust metrics
      procedure :: prepare_weno                             !< Prepare weno interpolation coefficients
      procedure :: get_drhoEdt                              !< Calculate drhoE/dt
      procedure :: get_max                                  !< Calculate stats of our energy
      procedure :: solve_implicit                           !< Solve for the energy residuals implicitly
   end type energy
   
   
contains
   
   
   !> Initialization for energy solver
   subroutine initialize(this,cfg,name)
      use messager, only: die
      implicit none
      class(energy), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to pgrid object
      this%cfg=>cfg
      
      ! Nullify bcond list
      this%nbc=0
      this%first_bc=>NULL()
      
      ! Allocate variables
      allocate(this%E   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%E   =0.0_WP
      allocate(this%Eold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Eold=0.0_WP
      allocate(this%diff(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%diff=0.0_WP
      
      ! Check current overlap
      if (this%cfg%no.lt.2) call die('[energy constructor] energy transport scheme requires larger overlap')
      
      ! Prepare default metrics
      call this%init_metrics()
      
      ! Prepare mask for E
      allocate(this%mask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%mask=0
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.           1) this%mask(:this%cfg%imin-1,:,:)=2
         if (this%cfg%iproc.eq.this%cfg%npx) this%mask(this%cfg%imax+1:,:,:)=2
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.           1) this%mask(:,:this%cfg%jmin-1,:)=2
         if (this%cfg%jproc.eq.this%cfg%npy) this%mask(:,this%cfg%jmax+1:,:)=2
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.           1) this%mask(:,:,:this%cfg%kmin-1)=2
         if (this%cfg%kproc.eq.this%cfg%npz) this%mask(:,:,this%cfg%kmax+1:)=2
      end if
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%cfg%VF(i,j,k).eq.0.0_WP) this%mask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%mask)
      
   end subroutine initialize
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      use mathtools, only: fv_itp_build
      implicit none
      class(energy), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference diffusivity interpolation coefficients to cell faces
      allocate(this%itp_x(-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< X-face-centered
      allocate(this%itp_y(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Y-face-centered
      allocate(this%itp_z(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Z-face-centered
      ! Create interpolation coefficients to cell face in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itp_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
            end do
         end do
      end do
      ! Create interpolation coefficients to cell face in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itp_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
            end do
         end do
      end do
      ! Create interpolation coefficients to cell face in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itp_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Allocate finite difference energy interpolation coefficients
      allocate(this%itpe_x(-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< X-face-centered
      allocate(this%itpe_y(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Y-face-centered
      allocate(this%itpe_z(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Z-face-centered
      ! Create energy interpolation coefficients to cell face in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpe_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
            end do
         end do
      end do
      ! Create energy interpolation coefficients to cell face in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpe_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
            end do
         end do
      end do
      ! Create energy interpolation coefficients to cell face in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpe_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Allocate 2nd order upwind energy interpolation coefficients
      allocate(this%up2_xp(-2:-1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%up2_xm( 0:+1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%up2_yp(-2:-1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%up2_ym( 0:+1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%up2_zp(-2:-1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      allocate(this%up2_zm( 0:+1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create 2nd order upwind energy interpolation coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Interpolation to x-face
               call fv_itp_build(n=2,x=this%cfg%x(i-2:i  ),xp=this%cfg%x(i),coeff=this%up2_xp(:,i,j,k))
               call fv_itp_build(n=2,x=this%cfg%x(i  :i+2),xp=this%cfg%x(i),coeff=this%up2_xm(:,i,j,k))
               ! Interpolation to y-face
               call fv_itp_build(n=2,x=this%cfg%y(j-2:j  ),xp=this%cfg%y(j),coeff=this%up2_yp(:,i,j,k))
               call fv_itp_build(n=2,x=this%cfg%y(j  :j+2),xp=this%cfg%y(j),coeff=this%up2_ym(:,i,j,k))
               ! Interpolation to z-face
               call fv_itp_build(n=2,x=this%cfg%z(k-2:k  ),xp=this%cfg%z(k),coeff=this%up2_zp(:,i,j,k))
               call fv_itp_build(n=2,x=this%cfg%z(k  :k+2),xp=this%cfg%z(k),coeff=this%up2_zm(:,i,j,k))
            end do
         end do
      end do
      
      ! Allocate WENO3 energy interpolation coefficients
      allocate(this%weno_xp(-2:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%weno_xm(-1:1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%weno_yp(-2:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%weno_ym(-1:1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%weno_zp(-2:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      allocate(this%weno_zm(-1:1,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      
      ! Allocate finite volume divergence operators
      allocate(this%div_x(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%div_y(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%div_z(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%div_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%div_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do
      
      ! Allocate finite difference energy gradient operators
      allocate(this%grd_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%grd_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%grd_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%grd_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of E in x from [xm,ym,zm] to [x,ym,zm]
               this%grd_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of E in y from [xm,ym,zm] to [xm,y,zm]
               this%grd_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of E in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(energy), intent(inout) :: this
      integer :: i,j,k
      
      ! Sync up masks
      call this%cfg%sync(this%mask)
      
      ! Adjust interpolation coefficients to cell faces - Neumann is used at walls only,
      ! leaving intact interpolation coefficients at Dirichlet boundaries (mask=2)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Linear interpolation in x
               if (this%mask(i,j,k).eq.0.and.this%mask(i-1,j,k).eq.1) this%itp_x(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).eq.1.and.this%mask(i-1,j,k).eq.0) this%itp_x(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in y
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j-1,k).eq.1) this%itp_y(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).eq.1.and.this%mask(i,j-1,k).eq.0) this%itp_y(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in z
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k-1).eq.1) this%itp_z(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).eq.1.and.this%mask(i,j,k-1).eq.0) this%itp_z(:,i,j,k)=[1.0_WP,0.0_WP]
            end do
         end do
      end do
      
      ! Adjust energy interpolation to reflect Dirichlet boundaries
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face
               if (this%mask(i-1,j,k).eq.2) then
                  this%up2_xm(:,i,j,k)=0.0_WP; this%up2_xm(-1,i,j,k)=1.0_WP
                  this%up2_xp(:,i,j,k)=0.0_WP; this%up2_xp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i  ,j,k).eq.2) then
                  this%up2_xm(:,i,j,k)=0.0_WP; this%up2_xm( 0,i,j,k)=1.0_WP
                  this%up2_xp(:,i,j,k)=0.0_WP; this%up2_xp( 0,i,j,k)=1.0_WP
               end if
               ! Y face
               if (this%mask(i,j-1,k).eq.2) then
                  this%up2_ym(:,i,j,k)=0.0_WP; this%up2_ym(-1,i,j,k)=1.0_WP
                  this%up2_yp(:,i,j,k)=0.0_WP; this%up2_yp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j  ,k).eq.2) then
                  this%up2_ym(:,i,j,k)=0.0_WP; this%up2_ym( 0,i,j,k)=1.0_WP
                  this%up2_yp(:,i,j,k)=0.0_WP; this%up2_yp( 0,i,j,k)=1.0_WP
               end if
               ! Z face
               if (this%mask(i,j,k-1).eq.2) then
                  this%up2_zm(:,i,j,k)=0.0_WP; this%up2_zm(-1,i,j,k)=1.0_WP
                  this%up2_zp(:,i,j,k)=0.0_WP; this%up2_zp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j,k  ).eq.2) then
                  this%up2_zm(:,i,j,k)=0.0_WP; this%up2_zm( 0,i,j,k)=1.0_WP
                  this%up2_zp(:,i,j,k)=0.0_WP; this%up2_zp( 0,i,j,k)=1.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to divergence of energy flux
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).gt.0) then
                  this%div_x(:,i,j,k)=0.0_WP
                  this%div_y(:,i,j,k)=0.0_WP
                  this%div_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell faces for walls (assume Neumann at wall) - NEEDS TO BE CHANGED AS WE MAY WANT ISOTHERMAL WALLS
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%mask(i,j,k).eq.1.or.this%mask(i-1,j,k).eq.1) this%grd_x(:,i,j,k)=0.0_WP     !< FD gradient in x of T
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j-1,k).eq.1) this%grd_y(:,i,j,k)=0.0_WP     !< FD gradient in y of T
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j,k-1).eq.1) this%grd_z(:,i,j,k)=0.0_WP     !< FD gradient in z of T
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%div_x=0.0_WP
         this%grd_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%div_y=0.0_WP
         this%grd_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%div_z=0.0_WP
         this%grd_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the energy solver now that bconds have been defined
   subroutine setup(this,implicit_solver)
      implicit none
      class(energy), intent(inout) :: this
      class(linsol), target, intent(in), optional :: implicit_solver
      
      ! Adjust metrics based on mask array
      call this%adjust_metrics()
      
      ! Prepare implicit solver if it had been provided
      if (present(implicit_solver)) then
         
         ! Point to implicit solver linsol object
         this%implicit=>implicit_solver

         ! Set 13-pt stencil map for the energy solver
         this%implicit%stc( 1,:)=[ 0, 0, 0]
         this%implicit%stc( 2,:)=[+1, 0, 0]
         this%implicit%stc( 3,:)=[-1, 0, 0]
         this%implicit%stc( 4,:)=[ 0,+1, 0]
         this%implicit%stc( 5,:)=[ 0,-1, 0]
         this%implicit%stc( 6,:)=[ 0, 0,+1]
         this%implicit%stc( 7,:)=[ 0, 0,-1]
         this%implicit%stc( 8,:)=[+2, 0, 0]
         this%implicit%stc( 9,:)=[-2, 0, 0]
         this%implicit%stc(10,:)=[ 0,+2, 0]
         this%implicit%stc(11,:)=[ 0,-2, 0]
         this%implicit%stc(12,:)=[ 0, 0,+2]
         this%implicit%stc(13,:)=[ 0, 0,-2]
         
         ! Set the diagonal to 1 to make sure all cells participate in solver
         this%implicit%opr(1,:,:,:)=1.0_WP
         
         ! Initialize the implicit velocity solver
         call this%implicit%init()
         
      else
         
         ! Point to implicit solver linsol object
         this%implicit=>NULL()
         
      end if
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,         only: lowercase
      use messager,       only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(energy), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: type
      procedure(locator_ftype) :: locator
      character(len=2), optional :: dir
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      if (present(dir)) then
         select case (lowercase(dir))
         case ('+x','x+','xp','px'); new_bc%dir=1
         case ('-x','x-','xm','mx'); new_bc%dir=2
         case ('+y','y+','yp','py'); new_bc%dir=3
         case ('-y','y-','ym','my'); new_bc%dir=4
         case ('+z','z+','zp','pz'); new_bc%dir=5
         case ('-z','z-','zm','mz'); new_bc%dir=6
         case default; call die('[energy add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[energy apply_bcond] Neumann requires a direction')
         new_bc%dir=0
      end if
      new_bc%itr=iterator(this%cfg,new_bc%name,locator,'c')
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         do n=1,new_bc%itr%n_
            i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
            this%mask(i,j,k)=2
         end do
      case (neumann)
         ! No modification - this assumes Neumann is only applied at walls or domain boundaries
      case default
         call die('[energy apply_bcond] Unknown bcond type')
      end select
      
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(energy), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[energy get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(energy), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
            case (dirichlet)           ! Apply Dirichlet conditions
               
               ! This is done by the user directly
               ! Unclear whether we want to do this within the solver...
               
            case (neumann)             ! Apply Neumann condition
               
               ! Implement based on bcond direction
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%E(i,j,k)=this%E(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[energy apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! Sync full field after all bcond
      call this%cfg%sync(this%E)
      
   end subroutine apply_bcond
   
   
   !> Prepare WENO stencil
   subroutine prepare_weno(this,e)
      use mathtools, only: fv_itp_build
      implicit none
      class(energy), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: e !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), parameter :: eps=1.0e-15_WP
      integer :: i,j,k
      real(WP) :: w
      ! Analyze passed E array and prepare interpolation coefficients
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X-face
               w=weno_weight((abs(e(i-1,j,k)-e(i-2,j,k))+eps)/(abs(e(i  ,j,k)-e(i-1,j,k))+eps))
               this%weno_xp(:,i,j,k)=(       w)*[this%up2_xp(-2,i,j,k),this%up2_xp(-1,i,j,k),0.0_WP               ]&
               &                    +(1.0_WP-w)*[0.0_WP               ,this%itpe_x(-1,i,j,k),this%itpe_x( 0,i,j,k)]
               w=weno_weight((abs(e(i+1,j,k)-e(i  ,j,k))+eps)/(abs(e(i  ,j,k)-e(i-1,j,k))+eps))
               this%weno_xm(:,i,j,k)=(       w)*[0.0_WP               ,this%up2_xm( 0,i,j,k),this%up2_xm(+1,i,j,k)]&
               &                    +(1.0_WP-w)*[this%itpe_x(-1,i,j,k),this%itpe_x( 0,i,j,k),0.0_WP               ]
               ! Y-face
               w=weno_weight((abs(e(i,j-1,k)-e(i,j-2,k))+eps)/(abs(e(i,j  ,k)-e(i,j-1,k))+eps))
               this%weno_yp(:,i,j,k)=(       w)*[this%up2_yp(-2,i,j,k),this%up2_yp(-1,i,j,k),0.0_WP               ]&
               &                    +(1.0_WP-w)*[0.0_WP               ,this%itpe_y(-1,i,j,k),this%itpe_y( 0,i,j,k)]
               w=weno_weight((abs(e(i,j+1,k)-e(i,j  ,k))+eps)/(abs(e(i,j  ,k)-e(i,j-1,k))+eps))
               this%weno_ym(:,i,j,k)=(       w)*[0.0_WP               ,this%up2_ym( 0,i,j,k),this%up2_ym(+1,i,j,k)]&
               &                    +(1.0_WP-w)*[this%itpe_y(-1,i,j,k),this%itpe_y( 0,i,j,k),0.0_WP               ]
               ! Z-face
               w=weno_weight((abs(e(i,j,k-1)-e(i,j,k-2))+eps)/(abs(e(i,j,k  )-e(i,j,k-1))+eps))
               this%weno_zp(:,i,j,k)=(       w)*[this%up2_zp(-2,i,j,k),this%up2_zp(-1,i,j,k),0.0_WP               ]&
               &                    +(1.0_WP-w)*[0.0_WP               ,this%itpe_z(-1,i,j,k),this%itpe_z( 0,i,j,k)]
               w=weno_weight((abs(e(i,j,k+1)-e(i,j,k  ))+eps)/(abs(e(i,j,k  )-e(i,j,k-1))+eps))
               this%weno_zm(:,i,j,k)=(       w)*[0.0_WP               ,this%up2_zm( 0,i,j,k),this%up2_zm(+1,i,j,k)]&
               &                    +(1.0_WP-w)*[this%itpe_z(-1,i,j,k),this%itpe_z( 0,i,j,k),0.0_WP               ]
            end do
         end do
      end do
   contains
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP ! Switching parameter
         real(WP), parameter :: delta=0.01_WP  ! Switching thickness
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight
   end subroutine prepare_weno
   
   
   !> Calculate the explicit rhoE time derivative given passed RHO/RHOold/rhoU/rhoV/rhoW
   subroutine get_drhoEdt(this,drhoEdt,rhoU,rhoV,rhoW)
      implicit none
      class(energy), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoEdt!< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoU   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoV   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoW   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ,E_
      ! Zero out drhoE/dt array
      drhoEdt=0.0_WP
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Allocate and compute mid-time energy
      allocate(E_(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); E_=0.5_WP*(this%E+this%Eold)
      ! Flux of rhoE
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%weno_xp(:,i,j,k)*E_(i-2:i  ,j,k)) &
               &         -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%weno_xm(:,i,j,k)*E_(i-1:i+1,j,k)) &
               &         +sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k))*sum(this%grd_x(:,i,j,k)*E_(i-1:i,j,k))/this%Cv
               FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%weno_yp(:,i,j,k)*E_(i,j-2:j  ,k)) &
               &         -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%weno_ym(:,i,j,k)*E_(i,j-1:j+1,k)) &
               &         +sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k))*sum(this%grd_y(:,i,j,k)*E_(i,j-1:j,k))/this%Cv
               FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%weno_zp(:,i,j,k)*E_(i,j,k-2:k  )) &
               &         -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%weno_zm(:,i,j,k)*E_(i,j,k-1:k+1)) &
               &         +sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k))*sum(this%grd_z(:,i,j,k)*E_(i,j,k-1:k))/this%Cv
            end do
         end do
      end do
      ! Time derivative of rhoE
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoEdt(i,j,k)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync residual
      call this%cfg%sync(drhoEdt)
      ! Deallocate work arrays
      deallocate(FX,FY,FZ,E_)
   end subroutine get_drhoEdt
   
   
   !> Calculate the int, min, and max of our E field
   subroutine get_max(this,rho)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(energy), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho
      integer :: ierr,i,j,k
      call this%cfg%integrate(rho*this%E,integral=this%rhoEint)
      this%Emax=-huge(1.0_WP)
      this%Emin=+huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip only walls
               if (this%mask(i,j,k).ne.1) then
                  this%Emax=max(this%E(i,j,k),this%Emax)
                  this%Emin=min(this%E(i,j,k),this%Emin)
               end if
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Emax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Emin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine get_max
   
   
   !> Solve for implicit energy residual
   subroutine solve_implicit(this,dt,resE,RHO,rhoU,rhoV,rhoW)
      implicit none
      class(energy), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resE   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: RHO    !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoU   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoV   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoW   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,std,sti
      
      ! If no implicit solver available, just divide by density and return
      if (.not.associated(this%implicit)) then
         resE=resE/RHO
         call this%cfg%sync(resE)
         return
      end if
      
      ! Prepare convective operator
      this%implicit%opr(1,:,:,:)=RHO; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=-2,0
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%div_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)+abs(rhoU(i+std,j,k)))*this%weno_xp(sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%div_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)+abs(rhoV(i,j+std,k)))*this%weno_yp(sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%div_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)+abs(rhoW(i,j,k+std)))*this%weno_zp(sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=-1,+1
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%div_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)-abs(rhoU(i+std,j,k)))*this%weno_xm(sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%div_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)-abs(rhoV(i,j+std,k)))*this%weno_ym(sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%div_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)-abs(rhoW(i,j,k+std)))*this%weno_zm(sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%div_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grd_x(-1,i+1,j,k)+&
               &                                                                this%div_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grd_x( 0,i  ,j,k)+&
               &                                                                this%div_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grd_y(-1,i,j+1,k)+&
               &                                                                this%div_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grd_y( 0,i,j  ,k)+&
               &                                                                this%div_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grd_z(-1,i,j,k+1)+&
               &                                                                this%div_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grd_z( 0,i,j,k  ))/this%Cv
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%div_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grd_x( 0,i+1,j,k))/this%Cv
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%div_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grd_x(-1,i  ,j,k))/this%Cv
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%div_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grd_y( 0,i,j+1,k))/this%Cv
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%div_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grd_y(-1,i,j  ,k))/this%Cv
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%div_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grd_z( 0,i,j,k+1))/this%Cv
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%div_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grd_z(-1,i,j,k  ))/this%Cv
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resE
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resE=this%implicit%sol
      
   end subroutine solve_implicit
   
   
   !> Print out info for energy solver
   subroutine energy_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(energy), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Energy solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine energy_print
   
   
end module energy_class
