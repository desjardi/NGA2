!> Incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density.
module incomp_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: incomp
   
   !> Incompressible solver object definition
   type :: incomp
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_INCOMP'  !< Solver name (default=UNNAMED_INCOMP)
      
      ! Constant property fluid
      real(WP) :: rho                                     !< This is our constant fluid density
      real(WP) :: visc                                    !< These is our constant fluid dynamic viscosity
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itpu_x ,itpu_y ,itpu_z   !< Interpolation for U
      real(WP), dimension(:,:,:,:), allocatable :: itpv_x ,itpv_y ,itpv_z   !< Interpolation for V
      real(WP), dimension(:,:,:,:), allocatable :: itpw_x ,itpw_y ,itpw_z   !< Interpolation for W
      real(WP), dimension(:,:,:,:), allocatable :: divp_x ,divp_y ,divp_z   !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: divu_x ,divu_y ,divu_z   !< Divergence for U-cell
      real(WP), dimension(:,:,:,:), allocatable :: divv_x ,divv_y ,divv_z   !< Divergence for V-cell
      real(WP), dimension(:,:,:,:), allocatable :: divw_x ,divw_y ,divw_z   !< Divergence for W-cell
      real(WP), dimension(:,:,:,:), allocatable :: gradu_x,gradu_y,gradu_z  !< Velocity gradient for U
      real(WP), dimension(:,:,:,:), allocatable :: gradv_x,gradv_y,gradv_z  !< Velocity gradient for V
      real(WP), dimension(:,:,:,:), allocatable :: gradw_x,gradw_y,gradw_z  !< Velocity gradient for W
      
      ! Boundary condition list
      !type(bcond), pointer :: bc=NULL()                   !<
      
      ! MAYBE SOME WORK ARRAYS TOO?
      
   contains
      procedure :: print=>incomp_print                    !< Output solver to the screen
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: mask_metrics                           !< Apply masks to metrics
   end type incomp
   
   
   !> Declare incompressible solver constructor
   interface incomp
      procedure constructor
   end interface incomp
   
contains
   
   
   !> Default constructor for incompressible flow solver
   function constructor(cfg,name) result(self)
      implicit none
      type(incomp) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Prepare metrics
      call self%init_metrics()
      call self%mask_metrics()
      
      ! Allocate flow variables
      allocate(self%U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%U=0.0_WP
      allocate(self%V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%V=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%W=0.0_WP
      allocate(self%P(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P=0.0_WP
      
   end function constructor
   
   
   !> Metric initialization (allocation and setup for 2nd order without BC/walls)
   subroutine init_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference velocity interpolation coefficients
      allocate(this%itpu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_-1,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%itpu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%itpu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  ))
      allocate(this%itpv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%itpv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_-1,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%itpv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  ))
      allocate(this%itpw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%itpw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%itpw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_-1))
      
      ! Create velocity interpolation coefficients in x
      do i=this%cfg%imino_  ,this%cfg%imaxo_-1
         ! Linear interpolation of U from [x,ym,zm] to [xm,ym,zm]
         this%itpu_x( 0,i,:,:)=this%cfg%dxi (i)*(this%cfg%x (i+1)-this%cfg%xm(i))
         this%itpu_x(+1,i,:,:)=this%cfg%dxi (i)*(this%cfg%xm(i+1)-this%cfg%x (i))
      end do
      do i=this%cfg%imino_+1,this%cfg%imaxo_
         ! Linear interpolation of V from [xm,y,zm] to [x,y,zm]
         this%itpv_x(-1,i,:,:)=this%cfg%dxmi(i)*(this%cfg%xm(i)-this%cfg%x (i  ))
         this%itpv_x( 0,i,:,:)=this%cfg%dxmi(i)*(this%cfg%x (i)-this%cfg%xm(i-1))
         ! Linear interpolation of W from [xm,ym,z] to [x,ym,z]
         this%itpw_x(-1,i,:,:)=this%cfg%dxmi(i)*(this%cfg%xm(i)-this%cfg%x (i  ))
         this%itpw_x( 0,i,:,:)=this%cfg%dxmi(i)*(this%cfg%x (i)-this%cfg%xm(i-1))
      end do
      
      ! Create velocity interpolation coefficients in y
      do j=this%cfg%jmino_  ,this%cfg%jmaxo_-1
         ! Linear interpolation of V from [xm,y,zm] to [xm,ym,zm]
         this%itpv_y( 0,:,j,:)=this%cfg%dyi (j)*(this%cfg%y (j+1)-this%cfg%ym(j))
         this%itpv_y(+1,:,j,:)=this%cfg%dyi (j)*(this%cfg%ym(j+1)-this%cfg%y (j))
      end do
      do j=this%cfg%jmino_+1,this%cfg%jmaxo_
         ! Linear interpolation of U from [x,ym,zm] to [x,y,zm]
         this%itpu_y(-1,:,j,:)=this%cfg%dymi(j)*(this%cfg%ym(j)-this%cfg%y (j  ))
         this%itpu_y( 0,:,j,:)=this%cfg%dymi(j)*(this%cfg%y (j)-this%cfg%ym(j-1))
         ! Linear interpolation of W from [xm,ym,z] to [xm,y,z]
         this%itpw_y(-1,:,j,:)=this%cfg%dymi(j)*(this%cfg%ym(j)-this%cfg%y (j  ))
         this%itpw_y( 0,:,j,:)=this%cfg%dymi(j)*(this%cfg%y (j)-this%cfg%ym(j-1))
      end do
      
      ! Create velocity interpolation coefficients in z
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_-1
         ! Linear interpolation of W from [xm,ym,z] to [xm,ym,zm]
         this%itpw_z( 0,:,:,k)=this%cfg%dzi (k)*(this%cfg%z (k+1)-this%cfg%zm(k))
         this%itpw_z(+1,:,:,k)=this%cfg%dzi (k)*(this%cfg%zm(k+1)-this%cfg%z (k))
      end do
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         ! Linear interpolation of U from [x,ym,zm] to [x,ym,z]
         this%itpu_z(-1,:,:,k)=this%cfg%dzmi(k)*(this%cfg%zm(k)-this%cfg%z (k  ))
         this%itpu_z( 0,:,:,k)=this%cfg%dzmi(k)*(this%cfg%z (k)-this%cfg%zm(k-1))
         ! Linear interpolation of V from [xm,y,zm] to [xm,y,z]
         this%itpv_z(-1,:,:,k)=this%cfg%dzmi(k)*(this%cfg%zm(k)-this%cfg%z (k  ))
         this%itpv_z( 0,:,:,k)=this%cfg%dzmi(k)*(this%cfg%z (k)-this%cfg%zm(k-1))
      end do
      
      ! Allocate finite volume divergence operators
      allocate(this%divp_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_-1,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divp_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_-1,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divp_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_-1))
      allocate(this%divu_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divu_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_-1,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divu_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_-1))
      allocate(this%divv_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_-1,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divv_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divv_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_-1))
      allocate(this%divw_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_-1,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divw_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_-1,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%divw_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  ))
      
      ! Create divergence coefficients in x
      do i=this%cfg%imino_  ,this%cfg%imaxo_-1
         ! FV divergence from [x,ym,zm] to [xm,ym,zm]
         this%divp_x( 0,i,:,:)=-this%cfg%dxi (i)
         this%divp_x(+1,i,:,:)=+this%cfg%dxi (i)
         ! FV divergence from [x,y,zm] to [xm,y,zm]
         this%divv_x( 0,i,:,:)=-this%cfg%dxi (i)
         this%divv_x(+1,i,:,:)=+this%cfg%dxi (i)
         ! FV divergence from [x,ym,z] to [xm,ym,z]
         this%divw_x( 0,i,:,:)=-this%cfg%dxi (i)
         this%divw_x(+1,i,:,:)=+this%cfg%dxi (i)
      end do
      do i=this%cfg%imino_+1,this%cfg%imaxo_
         ! FV divergence from [xm,ym,zm] to [x,ym,zm]
         this%divu_x(-1,i,:,:)=-this%cfg%dxmi(i)
         this%divu_x( 0,i,:,:)=+this%cfg%dxmi(i)
      end do
      
      ! Create divergence coefficients in y
      do j=this%cfg%jmino_  ,this%cfg%jmaxo_-1
         ! FV divergence from [xm,y,zm] to [xm,ym,zm]
         this%divp_y( 0,:,j,:)=-this%cfg%dyi (j)
         this%divp_y(+1,:,j,:)=+this%cfg%dyi (j)
         ! FV divergence from [x,y,zm] to [x,ym,zm]
         this%divu_y( 0,:,j,:)=-this%cfg%dyi (j)
         this%divu_y(+1,:,j,:)=+this%cfg%dyi (j)
         ! FV divergence from [xm,y,z] to [xm,ym,z]
         this%divw_y( 0,:,j,:)=-this%cfg%dyi (j)
         this%divw_y(+1,:,j,:)=+this%cfg%dyi (j)
      end do
      do j=this%cfg%jmino_+1,this%cfg%jmaxo_
         ! FV divergence from [xm,ym,zm] to [xm,y,zm]
         this%divv_y(-1,:,j,:)=-this%cfg%dymi(j)
         this%divv_y( 0,:,j,:)=+this%cfg%dymi(j)
      end do
      
      ! Create divergence coefficients in z
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_-1
         ! FV divergence from [xm,ym,z] to [xm,ym,zm]
         this%divp_z( 0,:,:,k)=-this%cfg%dzi (k)
         this%divp_z(+1,:,:,k)=+this%cfg%dzi (k)
         ! FV divergence from [x,ym,z] to [x,ym,zm]
         this%divu_z( 0,:,:,k)=-this%cfg%dzi (k)
         this%divu_z(+1,:,:,k)=+this%cfg%dzi (k)
         ! FV divergence from [xm,y,z] to [xm,y,zm]
         this%divv_z( 0,:,:,k)=-this%cfg%dzi (k)
         this%divv_z(+1,:,:,k)=+this%cfg%dzi (k)
      end do
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         ! FV divergence from [xm,ym,zm] to [xm,ym,z]
         this%divw_z(-1,:,:,k)=-this%cfg%dzmi(k)
         this%divw_z( 0,:,:,k)=+this%cfg%dzmi(k)
      end do
      
      ! Allocate finite different velocity gradient operators
      allocate(this%gradu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_-1,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%gradu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%gradu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  ))
      allocate(this%gradv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%gradv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_-1,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%gradv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  ))
      allocate(this%gradw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%gradw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  ))
      allocate(this%gradw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_-1))
      
      ! Create gradient coefficients in x
      do i=this%cfg%imino_+1,this%cfg%imaxo_
         ! FD gradient from [xm,y,zm] to [x,y,zm]
         this%gradv_x(-1,i,:,:)=-this%cfg%dxmi(i)
         this%gradv_x( 0,i,:,:)=+this%cfg%dxmi(i)
         ! FD gradient from [xm,ym,z] to [x,ym,z]
         this%gradw_x(-1,i,:,:)=-this%cfg%dxmi(i)
         this%gradw_x( 0,i,:,:)=+this%cfg%dxmi(i)
      end do
      do i=this%cfg%imino_  ,this%cfg%imaxo_-1
         ! FD gradient from [x,ym,zm] to [xm,ym,zm]
         this%gradu_x( 0,i,:,:)=-this%cfg%dxi (i)
         this%gradu_x(+1,i,:,:)=+this%cfg%dxi (i)
      end do
      
      ! Create gradient coefficients in y
      do j=this%cfg%jmino_+1,this%cfg%jmaxo_
         ! FD gradient from [x,ym,zm] to [x,y,zm]
         this%gradu_y(-1,:,j,:)=-this%cfg%dymi(j)
         this%gradu_y( 0,:,j,:)=+this%cfg%dymi(j)
         ! FD gradient from [xm,ym,z] to [xm,y,z]
         this%gradw_y(-1,:,j,:)=-this%cfg%dymi(j)
         this%gradw_y( 0,:,j,:)=+this%cfg%dymi(j)
      end do
      do j=this%cfg%jmino_  ,this%cfg%jmaxo_-1
         ! FD gradient from [xm,y,zm] to [xm,ym,zm]
         this%gradv_y( 0,:,j,:)=-this%cfg%dyi (j)
         this%gradv_y(+1,:,j,:)=+this%cfg%dyi (j)
      end do
      
      ! Create gradient coefficients in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         ! FD gradient from [x,ym,zm] to [x,ym,z]
         this%gradu_z(-1,:,:,k)=-this%cfg%dzmi(k)
         this%gradu_z( 0,:,:,k)=+this%cfg%dzmi(k)
         ! FD gradient from [xm,y,zm] to [xm,y,z]
         this%gradv_z(-1,:,:,k)=-this%cfg%dzmi(k)
         this%gradv_z( 0,:,:,k)=+this%cfg%dzmi(k)
      end do
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_-1
         ! FD gradient from [xm,ym,z] to [xm,ym,zm]
         this%gradw_z( 0,:,:,k)=-this%cfg%dzi (k)
         this%gradw_z(+1,:,:,k)=+this%cfg%dzi (k)
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment for masks
   subroutine mask_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      
      ! Adjust velocity interpolation coefficients for walls
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Check if [i,j,k] is a wall
               if (nint(this%cfg%mask(i,j,k)).eq.1) then
                  ! Zero out interpolators to [xm,ym,zm]
                  if (i.lt.this%cfg%imaxo_) this%itpu_x(:,i,j,k)=0.0_WP
                  if (j.lt.this%cfg%jmaxo_) this%itpv_y(:,i,j,k)=0.0_WP
                  if (k.lt.this%cfg%kmaxo_) this%itpw_z(:,i,j,k)=0.0_WP
                  ! Zero out interpolators to [xm,y,z]
                  if (j.gt.this%cfg%jmino_) this%itpw_y(:,i,j,k)=0.0_WP
                  if (k.gt.this%cfg%kmino_) this%itpv_z(:,i,j,k)=0.0_WP
                  ! Zero out interpolators to [x,ym,z]
                  if (k.gt.this%cfg%kmino_) this%itpu_z(:,i,j,k)=0.0_WP
                  if (i.gt.this%cfg%imino_) this%itpw_x(:,i,j,k)=0.0_WP
                  ! Zero out interpolators to [x,y,zm]
                  if (i.gt.this%cfg%imino_) this%itpv_x(:,i,j,k)=0.0_WP
                  if (j.gt.this%cfg%jmino_) this%itpu_y(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust divergence operators for walls
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Check if [i,j,k] is a wall
               if (nint(this%cfg%mask(i,j,k)).eq.1) then
                  ! Zero out divergence to [xm,ym,zm]
                  if (i.lt.this%cfg%imaxo_) this%divp_x(:,i,j,k)=0.0_WP
                  if (j.lt.this%cfg%jmaxo_) this%divp_y(:,i,j,k)=0.0_WP
                  if (k.lt.this%cfg%kmaxo_) this%divp_z(:,i,j,k)=0.0_WP
                  ! Zero out divergence to [x,ym,zm]
                  if (i.gt.this%cfg%imino_) this%divu_x(:,i,j,k)=0.0_WP
                  if (j.lt.this%cfg%jmaxo_) this%divu_y(:,i,j,k)=0.0_WP
                  if (k.lt.this%cfg%kmaxo_) this%divu_z(:,i,j,k)=0.0_WP
                  ! Zero out divergence to [xm,y,zm]
                  if (i.lt.this%cfg%imaxo_) this%divv_x(:,i,j,k)=0.0_WP
                  if (j.gt.this%cfg%jmino_) this%divv_y(:,i,j,k)=0.0_WP
                  if (k.lt.this%cfg%kmaxo_) this%divv_z(:,i,j,k)=0.0_WP
                  ! Zero out divergence to [xm,ym,z]
                  if (i.lt.this%cfg%imaxo_) this%divw_x(:,i,j,k)=0.0_WP
                  if (j.lt.this%cfg%jmaxo_) this%divw_y(:,i,j,k)=0.0_WP
                  if (k.gt.this%cfg%kmino_) this%divw_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust velocity gradient operators for walls
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Check if [i,j,k] is a wall
               if (nint(this%cfg%mask(i,j,k)).eq.1) then
                  ! Zero out gradients to [xm,ym,zm]
                  if (i.lt.this%cfg%imaxo_) this%gradu_x(:,i,j,k)=0.0_WP
                  if (j.lt.this%cfg%jmaxo_) this%gradv_y(:,i,j,k)=0.0_WP
                  if (k.lt.this%cfg%kmaxo_) this%gradw_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      
      
      ! Create gradient coefficients in x
      do i=this%cfg%imino_+1,this%cfg%imaxo_
         ! FD gradient from [xm,y,zm] to [x,y,zm]
         this%gradv_x(-1,i,:,:)=-this%cfg%dxmi(i)
         this%gradv_x( 0,i,:,:)=+this%cfg%dxmi(i)
         ! FD gradient from [xm,ym,z] to [x,ym,z]
         this%gradw_x(-1,i,:,:)=-this%cfg%dxmi(i)
         this%gradw_x( 0,i,:,:)=+this%cfg%dxmi(i)
      end do
      
      ! Create gradient coefficients in y
      do j=this%cfg%jmino_+1,this%cfg%jmaxo_
         ! FD gradient from [x,ym,zm] to [x,y,zm]
         this%gradu_y(-1,:,j,:)=-this%cfg%dymi(j)
         this%gradu_y( 0,:,j,:)=+this%cfg%dymi(j)
         ! FD gradient from [xm,ym,z] to [xm,y,z]
         this%gradw_y(-1,:,j,:)=-this%cfg%dymi(j)
         this%gradw_y( 0,:,j,:)=+this%cfg%dymi(j)
      end do
      
      ! Create gradient coefficients in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         ! FD gradient from [x,ym,zm] to [x,ym,z]
         this%gradu_z(-1,:,:,k)=-this%cfg%dzmi(k)
         this%gradu_z( 0,:,:,k)=+this%cfg%dzmi(k)
         ! FD gradient from [xm,y,zm] to [xm,y,z]
         this%gradv_z(-1,:,:,k)=-this%cfg%dzmi(k)
         this%gradv_z( 0,:,:,k)=+this%cfg%dzmi(k)
      end do
      
      
      
   end subroutine mask_metrics
   
   
   !> Print out info for incompressible flow solver
   subroutine incomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(incomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   density = ",es12.5)') this%rho
         write(output_unit,'(" > viscosity = ",es12.5)') this%visc
      end if
      
   end subroutine incomp_print
   
   
end module incomp_class
