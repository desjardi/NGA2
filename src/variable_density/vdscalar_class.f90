!> Variable density scalar solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes variable diffusivity and density.
module vdscalar_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ils_class,      only: ils
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: vdscalar,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2                !< Dirichlet condition
   integer, parameter, public :: neumann=3                  !< Zero normal gradient
   
   ! List of available advection schemes for scalar transport
   integer, parameter, public :: quick=1                    !< Quick scheme
   
   
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
   
   !> Variable density scalar solver object definition
   type :: vdscalar
      
      ! This is our config
      class(config), pointer :: cfg                         !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_VDSCALAR'  !< Solver name (default=UNNAMED_VDSCALAR)
      
      ! Variable property fluid
      real(WP), dimension(:,:,:), allocatable :: diff       !< These is our variable scalar dynamic diffusivity
      
      ! Boundary condition list
      integer :: nbc                                        !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                      !< List of bcond for our solver
      
      ! Scalar variable
      real(WP), dimension(:,:,:), allocatable :: rho        !< Density array
      real(WP), dimension(:,:,:), allocatable :: SC         !< SC array
      real(WP), dimension(:,:,:), allocatable :: rhoSC      !< rhoSC array
      
      ! Old scalar variable
      real(WP), dimension(:,:,:), allocatable :: rhoold     !< Old density array
      real(WP), dimension(:,:,:), allocatable :: SCold      !< SCold array
      real(WP), dimension(:,:,:), allocatable :: rhoSCold   !< rhoSCold array
      
      ! Implicit scalar solver
      type(ils) :: implicit                                 !< Iterative linear solver object for an implicit prediction of the scalar residual
      integer, dimension(:,:,:), allocatable :: stmap       !< Inverse map from stencil shift to index location
      
      ! Metrics
      integer :: scheme                                     !< Advection scheme for scalar
      integer :: nst                                        !< Scheme order (and elemental stencil size)
      integer :: stp1,stp2                                  !< Plus interpolation stencil extent for scalar advection
      integer :: stm1,stm2                                  !< Minus interpolation stencil extent for scalar advection
      real(WP), dimension(:,:,:,:), allocatable :: itpsc_xp,itpsc_yp,itpsc_zp        !< Plus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: itpsc_xm,itpsc_ym,itpsc_zm        !< Minus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: divsc_x ,divsc_y ,divsc_z         !< Divergence for SC
      real(WP), dimension(:,:,:,:), allocatable :: grdsc_x ,grdsc_y ,grdsc_z         !< Scalar gradient for SC
      real(WP), dimension(:,:,:,:), allocatable :: itp_x   ,itp_y   ,itp_z           !< Second order interpolation for SC diffusivity
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask        !< Integer array used for modifying SC metrics
      
      ! Monitoring quantities
      real(WP) :: SCmax,SCmin,SCint                         !< Maximum and minimum, integral scalar
      real(WP) :: rhomax,rhomin,rhoint                      !< Maximum and minimum, integral density
      real(WP) :: rhoSCmax,rhoSCmin,rhoSCint                !< Maximum and minimum, integral of rhoSC
      
   contains
      procedure :: print=>scalar_print                      !< Output solver to the screen
      procedure :: setup                                    !< Finish configuring the scalar solver
      procedure :: add_bcond                                !< Add a boundary condition
      procedure :: get_bcond                                !< Get a boundary condition
      procedure :: apply_bcond                              !< Apply all boundary conditions
      procedure :: init_metrics                             !< Initialize metrics
      procedure :: adjust_metrics                           !< Adjust metrics
      procedure :: get_drhoSCdt                             !< Calculate drhoSC/dt
      procedure :: get_max                                  !< Calculate maximum field values
      procedure :: get_int                                  !< Calculate integral field values
      procedure :: get_drhodt                               !< Calculate drhodt
      procedure :: solve_implicit                           !< Solve for the scalar residuals implicitly
      procedure :: rho_divide                               !< Divide rhoSC by rho to get SC
      procedure :: rho_multiply                             !< Multiply SC by rho to get rhoSC
   end type vdscalar
   
   
   !> Declare vdscalar solver constructor
   interface vdscalar
      procedure constructor
   end interface vdscalar
   
contains
   
   
   !> Default constructor for variable density scalar solver
   function constructor(cfg,scheme,name) result(self)
      use messager, only: die
      implicit none
      type(vdscalar) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: scheme
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate variables
      allocate(self%SC      (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SC      =0.0_WP
      allocate(self%SCold   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SCold   =0.0_WP
      allocate(self%rho     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho     =0.0_WP
      allocate(self%rhoold  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhoold  =0.0_WP
      allocate(self%rhoSC   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhoSC   =0.0_WP
      allocate(self%rhoSCold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhoSCold=0.0_WP
      allocate(self%diff    (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%diff    =0.0_WP
      
      ! Prepare advection scheme
      self%scheme=scheme
      select case (self%scheme)
      case (quick)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[scalar constructor] vdscalar transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=3
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case default
         call die('[scalar constructor] Unknown vdscalar transport scheme selected')
      end select
      
      ! Create implicit scalar solver object
      self%implicit=ils(cfg=self%cfg,name='Scalar',nst=1+6*abs(self%stp1))
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare mask for SC
      allocate(self%mask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%mask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%mask(:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(self%cfg%imax+1:,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%mask(:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(:,self%cfg%jmax+1:,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%mask(:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(:,:,self%cfg%kmax+1:)=2
      end if
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%mask)
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      use mathtools, only: fv_itp_build
      implicit none
      class(vdscalar), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference diffusivity interpolation coefficients
      allocate(this%itp_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create diffusivity interpolation coefficients to cell face
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%itp_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Allocate finite difference scalar interpolation coefficients
      allocate(this%itpsc_xp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itpsc_xm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itpsc_yp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itpsc_ym(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itpsc_zp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      allocate(this%itpsc_zm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create scalar interpolation coefficients to cell faces
      select case (this%scheme)
      case (quick)
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Interpolation to x-face
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itpsc_xp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itpsc_xm(:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itpsc_yp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itpsc_ym(:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itpsc_zp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itpsc_zm(:,i,j,k))
               end do
            end do
         end do
      end select
      
      ! Allocate finite volume divergence operators
      allocate(this%divsc_x(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%divsc_y(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%divsc_z(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%divsc_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divsc_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divsc_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(this%grdsc_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%grdsc_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%grdsc_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%grdsc_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in x from [xm,ym,zm] to [x,ym,zm]
               this%grdsc_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in y from [xm,ym,zm] to [xm,y,zm]
               this%grdsc_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(vdscalar), intent(inout) :: this
      integer :: i,j,k
      
      ! Sync up masks
      call this%cfg%sync(this%mask)
      
      ! Adjust interpolation coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Linear interpolation in x
               if (this%mask(i,j,k).eq.0.and.this%mask(i-1,j,k).gt.0) this%itp_x(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i-1,j,k).eq.0) this%itp_x(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in y
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j-1,k).gt.0) this%itp_y(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i,j-1,k).eq.0) this%itp_y(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in z
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k-1).gt.0) this%itp_z(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i,j,k-1).eq.0) this%itp_z(:,i,j,k)=[1.0_WP,0.0_WP]
            end do
         end do
      end do
      
      ! Adjust scalar interpolation to reflect Dirichlet boundaries
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face
               if (this%mask(i-1,j,k).eq.2) then
                  this%itpsc_xm(:,i,j,k)=0.0_WP; this%itpsc_xm(-1,i,j,k)=1.0_WP
                  this%itpsc_xp(:,i,j,k)=0.0_WP; this%itpsc_xp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i  ,j,k).eq.2) then
                  this%itpsc_xm(:,i,j,k)=0.0_WP; this%itpsc_xm( 0,i,j,k)=1.0_WP
                  this%itpsc_xp(:,i,j,k)=0.0_WP; this%itpsc_xp( 0,i,j,k)=1.0_WP
               end if
               ! Y face
               if (this%mask(i,j-1,k).eq.2) then
                  this%itpsc_ym(:,i,j,k)=0.0_WP; this%itpsc_ym(-1,i,j,k)=1.0_WP
                  this%itpsc_yp(:,i,j,k)=0.0_WP; this%itpsc_yp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j  ,k).eq.2) then
                  this%itpsc_ym(:,i,j,k)=0.0_WP; this%itpsc_ym( 0,i,j,k)=1.0_WP
                  this%itpsc_yp(:,i,j,k)=0.0_WP; this%itpsc_yp( 0,i,j,k)=1.0_WP
               end if
               ! Z face
               if (this%mask(i,j,k-1).eq.2) then
                  this%itpsc_zm(:,i,j,k)=0.0_WP; this%itpsc_zm(-1,i,j,k)=1.0_WP
                  this%itpsc_zp(:,i,j,k)=0.0_WP; this%itpsc_zp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j,k  ).eq.2) then
                  this%itpsc_zm(:,i,j,k)=0.0_WP; this%itpsc_zm( 0,i,j,k)=1.0_WP
                  this%itpsc_zp(:,i,j,k)=0.0_WP; this%itpsc_zp( 0,i,j,k)=1.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to SC divergence
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).gt.0) then
                  this%divsc_x(:,i,j,k)=0.0_WP
                  this%divsc_y(:,i,j,k)=0.0_WP
                  this%divsc_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell faces for walls (assume Neumann at wall)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%mask(i,j,k).eq.1.or.this%mask(i-1,j,k).eq.1) this%grdsc_x(:,i,j,k)=0.0_WP     !< FD gradient in x of SC
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j-1,k).eq.1) this%grdsc_y(:,i,j,k)=0.0_WP     !< FD gradient in y of SC
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j,k-1).eq.1) this%grdsc_z(:,i,j,k)=0.0_WP     !< FD gradient in z of SC
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%divsc_x=0.0_WP
         this%grdsc_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%divsc_y=0.0_WP
         this%grdsc_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%divsc_z=0.0_WP
         this%grdsc_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the variable density scalar solver now that bconds have been defined
   subroutine setup(this,implicit_ils)
      implicit none
      class(vdscalar), intent(inout) :: this
      integer, intent(in) :: implicit_ils
      integer :: count,st
      
      ! Adjust metrics based on mask array
      call this%adjust_metrics()
      
      ! Set dynamic stencil map for the scalar solver
      count=1; this%implicit%stc(count,:)=[0,0,0]
      do st=1,abs(this%stp1)
         count=count+1; this%implicit%stc(count,:)=[+st,0,0]
         count=count+1; this%implicit%stc(count,:)=[-st,0,0]
         count=count+1; this%implicit%stc(count,:)=[0,+st,0]
         count=count+1; this%implicit%stc(count,:)=[0,-st,0]
         count=count+1; this%implicit%stc(count,:)=[0,0,+st]
         count=count+1; this%implicit%stc(count,:)=[0,0,-st]
      end do
      
      ! Set the diagonal to 1 to make sure all cells participate in solver
      this%implicit%opr(1,:,:,:)=1.0_WP
      
      ! Initialize the implicit scalar solver
      call this%implicit%init(implicit_ils)
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(vdscalar), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: type
      external :: locator
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
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
         case default; call die('[vdscalar add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[vdscalar apply_bcond] Neumann requires a direction')
         new_bc%dir=0
      end if
      new_bc%itr=iterator(this%cfg,new_bc%name,locator)
      
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
         call die('[vdscalar apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(vdscalar), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[vdscalar get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vdscalar), intent(inout) :: this
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
                  this%SC(i,j,k)=this%SC(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  this%rhoSC(i,j,k)=this%rhoSC(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[vdscalar apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%SC)
         call this%cfg%sync(this%rhoSC)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit rhoSC time derivative based on rhoU/rhoV/rhoW
   subroutine get_drhoSCdt(this,drhoSCdt,rhoU,rhoV,rhoW)
      implicit none
      class(vdscalar), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoSCdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoU     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoV     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoW     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Flux of rhoSC
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%itpsc_xp(:,i,j,k)*this%SC(i+this%stp1:i+this%stp2,j,k)) &
               &         -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%itpsc_xm(:,i,j,k)*this%SC(i+this%stm1:i+this%stm2,j,k)) &
               &         +sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k))*sum(this%grdsc_x(:,i,j,k)*this%SC(i-1:i,j,k))
               ! Fluxes on y-face
               FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%itpsc_yp(:,i,j,k)*this%SC(i,j+this%stp1:j+this%stp2,k)) &
               &         -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%itpsc_ym(:,i,j,k)*this%SC(i,j+this%stm1:j+this%stm2,k)) &
               &         +sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k))*sum(this%grdsc_y(:,i,j,k)*this%SC(i,j-1:j,k))
               ! Fluxes on z-face
               FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%itpsc_zp(:,i,j,k)*this%SC(i,j,k+this%stp1:k+this%stp2)) &
               &         -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%itpsc_zm(:,i,j,k)*this%SC(i,j,k+this%stm1:k+this%stm2)) &
               &         +sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k))*sum(this%grdsc_z(:,i,j,k)*this%SC(i,j,k-1:k))
            end do
         end do
      end do
      ! Time derivative of rhoSC
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoSCdt(i,j,k)=sum(this%divsc_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &               sum(this%divsc_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &               sum(this%divsc_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      ! Sync residual
      call this%cfg%sync(drhoSCdt)
   end subroutine get_drhoSCdt
   
   
   !> Calculate the time derivative of rho
   subroutine get_drhodt(this,dt,drhodt)
      implicit none
      class(vdscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhodt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      drhodt=(this%rho-this%rhoold)/dt
   end subroutine get_drhodt
   
   
   !> Calculate the min and max of our SC field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vdscalar), intent(inout) :: this
      integer :: ierr,i,j,k
      real(WP) :: my_SCmax,my_SCmin,my_rhomax,my_rhomin,my_rhoSCmax,my_rhoSCmin
      my_SCmax   =-huge(1.0_WP)
      my_SCmin   =+huge(1.0_WP)
      my_rhomax  =-huge(1.0_WP)
      my_rhomin  =+huge(1.0_WP)
      my_rhoSCmax=-huge(1.0_WP)
      my_rhoSCmin=+huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip only walls
               if (this%mask(i,j,k).ne.1) then
                  my_SCmax   =max(   this%SC(i,j,k),my_SCmax   )
                  my_SCmin   =min(   this%SC(i,j,k),my_SCmin   )
                  my_rhomax  =max(  this%rho(i,j,k),my_rhomax  )
                  my_rhomin  =min(  this%rho(i,j,k),my_rhomin  )
                  my_rhoSCmax=max(this%rhoSC(i,j,k),my_rhoSCmax)
                  my_rhoSCmin=min(this%rhoSC(i,j,k),my_rhoSCmin)
               end if
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_SCmax   ,this%SCmax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_SCmin   ,this%SCmin   ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_rhomax  ,this%rhomax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_rhomin  ,this%rhomin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_rhoSCmax,this%rhoSCmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_rhoSCmin,this%rhoSCmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine get_max
   
   
   !> Calculate the integral of our SC field
   subroutine get_int(this)
      implicit none
      class(vdscalar), intent(inout) :: this
      ! Integrate SC
      call this%cfg%integrate(this%SC,integral=this%SCint)
      ! Integrate rho
      call this%cfg%integrate(this%rho,integral=this%rhoint)
      ! Integrate rhoSC
      call this%cfg%integrate(this%rhoSC,integral=this%rhoSCint)
   end subroutine get_int
   
   
   !> Divide rhoSC by rho to form SC
   subroutine rho_divide(this)
      implicit none
      class(vdscalar), intent(inout) :: this
      where (this%mask.eq.0) this%SC=this%rhoSC/this%rho
   end subroutine rho_divide
   
   
   !> Multiply SC by rho to form rhoSC
   subroutine rho_multiply(this)
      implicit none
      class(vdscalar), intent(inout) :: this
      integer :: i,j,k
      where (this%mask.eq.0) this%rhoSC=this%SC*this%rho
   end subroutine rho_multiply
   
   
   !> Solve for implicit vdscalar residual
   subroutine solve_implicit(this,dt,resSC,rhoU,rhoV,rhoW)
      implicit none
      class(vdscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resSC !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoU  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoV  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoW  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,sti,std
      
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1 ,i,j,k)=this%rho(i,j,k)
               this%implicit%opr(2:,i,j,k)=0.0_WP
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divsc_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)+abs(rhoU(i+std,j,k)))*this%itpsc_xp(sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divsc_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)+abs(rhoV(i,j+std,k)))*this%itpsc_yp(sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divsc_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)+abs(rhoW(i,j,k+std)))*this%itpsc_zp(sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divsc_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)-abs(rhoU(i+std,j,k)))*this%itpsc_xm(sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divsc_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)-abs(rhoV(i,j+std,k)))*this%itpsc_ym(sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divsc_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)-abs(rhoW(i,j,k+std)))*this%itpsc_zm(sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
      
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divsc_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grdsc_x(-1,i+1,j,k)+&
               &                                                                this%divsc_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grdsc_x( 0,i  ,j,k)+&
               &                                                                this%divsc_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grdsc_y(-1,i,j+1,k)+&
               &                                                                this%divsc_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grdsc_y( 0,i,j  ,k)+&
               &                                                                this%divsc_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grdsc_z(-1,i,j,k+1)+&
               &                                                                this%divsc_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grdsc_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divsc_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k))*this%grdsc_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divsc_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k))*this%grdsc_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divsc_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k))*this%grdsc_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divsc_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k))*this%grdsc_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divsc_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1))*this%grdsc_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divsc_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ))*this%grdsc_z(-1,i,j,k  ))
            end do
         end do
      end do
      
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resSC
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resSC=this%implicit%sol
      
      ! Sync up residual
      call this%cfg%sync(resSC)
      
   end subroutine solve_implicit
   
   
   !> Print out info for vdscalar solver
   subroutine scalar_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(vdscalar), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Variable density scalar solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine scalar_print
   
   
end module vdscalar_class
