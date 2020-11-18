!> Incompressible scalar solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes constant diffusivity and density.
module scalar_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ils_class,      only: ils
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: scalar,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   integer, parameter, public :: convective=4        !< Convective outflow condition
   integer, parameter, public :: clipped_neumann=5   !< Clipped Neumann condition (outflow only)
   
   ! List of available advection schemes for scalar transport
   integer, parameter, public :: quick=1             !< Quick scheme
   
   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      integer :: si,sj,sk                                 !< Index shift in the outward normal direction
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   
   !> Constant density scalar solver object definition
   type :: scalar
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_SCALAR'  !< Solver name (default=UNNAMED_SCALAR)
      
      ! Constant property fluid
      real(WP) :: rho                                     !< This is our constant fluid density
      real(WP) :: diff                                    !< These is our constant scalar dynamic diffusivity
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Scalar variable
      real(WP), dimension(:,:,:), allocatable :: SC       !< SC array
      
      ! Old scalar variable
      real(WP), dimension(:,:,:), allocatable :: SCold    !< SCold array
      
      ! Implicit scalar solver
      type(ils) :: implicit                               !< Iterative linear solver object for an implicit prediction of the scalar residual
      
      ! Metrics
      integer :: scheme                                   !< Advection scheme for scalar
      integer :: nst                                      !< Scheme order (and elemental stencil size)
      integer :: stp1,stp2                                !< Plus interpolation stencil extent for scalar advection
      integer :: stm1,stm2                                !< Minus interpolation stencil extent for scalar advection
      real(WP), dimension(:,:,:,:), allocatable :: itpsc_xp,itpsc_yp,itpsc_zp   !< Plus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: itpsc_xm,itpsc_ym,itpsc_zm   !< Minus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: divsc_x ,divsc_y ,divsc_z    !< Divergence for SC
      real(WP), dimension(:,:,:,:), allocatable :: grdsc_x ,grdsc_y ,grdsc_z    !< Scalar gradient for SC
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask                    !< Integer array used for modifying SC metrics
      
      ! Monitoring quantities
      real(WP) :: SCmax,SCmin                                             !< Maximum and minimum scalar
      
   contains
      procedure :: print=>scalar_print                    !< Output solver to the screen
      procedure :: setup                                  !< Finish configuring the scalar solver
      !procedure :: add_bcond                              !< Add a boundary condition
      !procedure :: get_bcond                              !< Get a boundary condition
      !procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: get_drhoSCdt                           !< Calculate drhoSC/dt
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: solve_implicit                         !< Solve for the scalar residuals implicitly
   end type scalar
   
   
   !> Declare scalar solver constructor
   interface scalar
      procedure constructor
   end interface scalar
   
contains
   
   
   !> Default constructor for scalar solver
   function constructor(cfg,scheme,name) result(self)
      implicit none
      type(scalar) :: self
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
      allocate(self%SC   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SC   =0.0_WP
      allocate(self%SCold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SCold=0.0_WP
      
      ! Prepare advection scheme
      self%scheme=scheme
      select case (self%scheme)
      case (quick)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[scalar constructor] Scalar transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=3
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case default
         call die('[scalar constructor] Unknown scalar transport scheme selected')
      end select
      
      ! Create implicit scalar solver object
      self%implicit=ils(cfg=self%cfg,name='Implicit scalar residual',nst=1+6*abs(self%stp1))
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare mask for SC
      allocate(self%mask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%mask=0
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%mask)
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%mask(self%cfg%imin  ,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(self%cfg%imax+1,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%mask(:,self%cfg%jmin  ,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(:,self%cfg%jmax+1,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%mask(:,:,self%cfg%kmin  )=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(:,:,self%cfg%kmax+1)=2
      end if
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      use mathtools, only: fv_itp_build
      implicit none
      class(scalar), intent(inout) :: this
      integer :: i,j,k
      
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
                  call fv_itp_build(n=3,x=this%cfg%x(i+stp1:i+stp2+1),xp=this%cfg%x(i),coeff=this%itpsc_xp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+stm1:i+stm2+1),xp=this%cfg%x(i),coeff=this%itpsc_xm(:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=3,x=this%cfg%y(j+stp1:j+stp2+1),xp=this%cfg%y(j),coeff=this%itpsc_yp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+stm1:j+stm2+1),xp=this%cfg%y(j),coeff=this%itpsc_ym(:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=3,x=this%cfg%z(k+stp1:k+stp2+1),xp=this%cfg%z(k),coeff=this%itpsc_zp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+stm1:k+stm2+1),xp=this%cfg%z(k),coeff=this%itpsc_zm(:,i,j,k))
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
   
   
   !> Metric adjustment accounting for bconds and walls
   subroutine adjust_metrics(this)
      implicit none
      class(scalar), intent(inout) :: this
      integer :: i,j,k
      real(WP) :: delta
      
      ! Sync up u/v/wmasks
      call this%cfg%sync(this%umask)
      call this%cfg%sync(this%vmask)
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
      ! Loop over the domain and adjust divergence for P cell
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%cfg%VF(i,j,k).eq.0.0_WP) then
                  this%divp_x(:,i,j,k)=0.0_WP
                  this%divp_y(:,i,j,k)=0.0_WP
                  this%divp_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to U metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (this%umask(i,j,k).gt.0) then
                  this%divu_x(:,i,j,k)=0.0_WP
                  this%divu_y(:,i,j,k)=0.0_WP
                  this%divu_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to V metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%vmask(i,j,k).gt.0) then
                  this%divv_x(:,i,j,k)=0.0_WP
                  this%divv_y(:,i,j,k)=0.0_WP
                  this%divv_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to W metrics
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%wmask(i,j,k).gt.0) then
                  this%divw_x(:,i,j,k)=0.0_WP
                  this%divw_y(:,i,j,k)=0.0_WP
                  this%divw_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! FD gradient in x of V from [xm,y ,zm]
               if (maxval(this%vmask(i-1:i,j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%vmask(i  ,j,k).eq.0) delta=delta+(this%cfg%xm(i)-this%cfg%x (i  ))
                  if (this%vmask(i-1,j,k).eq.0) delta=delta+(this%cfg%x (i)-this%cfg%xm(i-1))
                  if (delta.gt.0.0_WP) then
                     this%grdv_x(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdv_x(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in x of W from [xm,ym,z ]
               if (maxval(this%wmask(i-1:i,j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%wmask(i  ,j,k).eq.0) delta=delta+(this%cfg%xm(i)-this%cfg%x (i  ))
                  if (this%wmask(i-1,j,k).eq.0) delta=delta+(this%cfg%x (i)-this%cfg%xm(i-1))
                  if (delta.gt.0.0_WP) then
                     this%grdw_x(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdw_x(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! FD gradient in y of U from [x ,ym,zm]
               if (maxval(this%umask(i,j-1:j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%umask(i,j  ,k).eq.0) delta=delta+(this%cfg%ym(j)-this%cfg%y (j  ))
                  if (this%umask(i,j-1,k).eq.0) delta=delta+(this%cfg%y (j)-this%cfg%ym(j-1))
                  if (delta.gt.0.0_WP) then
                     this%grdu_y(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdu_y(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in y of W from [xm,ym,z ]
               if (maxval(this%wmask(i,j-1:j,k)).gt.0) then
                  delta=0.0_WP
                  if (this%wmask(i,j  ,k).eq.0) delta=delta+(this%cfg%ym(j)-this%cfg%y (j  ))
                  if (this%wmask(i,j-1,k).eq.0) delta=delta+(this%cfg%y (j)-this%cfg%ym(j-1))
                  if (delta.gt.0.0_WP) then
                     this%grdw_y(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdw_y(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! FD gradient in z of U from [x ,ym,zm]
               if (maxval(this%umask(i,j,k-1:k)).gt.0) then
                  delta=0.0_WP
                  if (this%umask(i,j,k  ).eq.0) delta=delta+(this%cfg%zm(k)-this%cfg%z (k  ))
                  if (this%umask(i,j,k-1).eq.0) delta=delta+(this%cfg%z (k)-this%cfg%zm(k-1))
                  if (delta.gt.0.0_WP) then
                     this%grdu_z(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdu_z(:,i,j,k)=0.0_WP
                  end if
               end if
               ! FD gradient in z of V from [xm,y ,zm]
               if (maxval(this%vmask(i,j,k-1:k)).gt.0) then
                  delta=0.0_WP
                  if (this%vmask(i,j,k  ).eq.0) delta=delta+(this%cfg%zm(k)-this%cfg%z (k  ))
                  if (this%vmask(i,j,k-1).eq.0) delta=delta+(this%cfg%z (k)-this%cfg%zm(k-1))
                  if (delta.gt.0.0_WP) then
                     this%grdv_z(:,i,j,k)=[-1.0_WP,+1.0_WP]/delta
                  else
                     this%grdv_z(:,i,j,k)=0.0_WP
                  end if
               end if
            end do
         end do
      end do
      
      ! Adjust interpolation coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! Linear interpolation in x of V from [xm,y ,zm]
               if (this%vmask(i,j,k).eq.0.and.this%vmask(i-1,j,k).gt.0) this%itpv_x(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%vmask(i,j,k).gt.0.and.this%vmask(i-1,j,k).eq.0) this%itpv_x(:,i,j,k)=[0.0_WP,1.0_WP]
               ! Linear interpolation in x of W from [xm,ym,z ]
               if (this%wmask(i,j,k).eq.0.and.this%wmask(i-1,j,k).gt.0) this%itpw_x(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%wmask(i,j,k).gt.0.and.this%wmask(i-1,j,k).eq.0) this%itpw_x(:,i,j,k)=[0.0_WP,1.0_WP]
            end do
         end do
      end do
      
      ! Adjust interpolation coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! Linear interpolation in y of U from [x ,ym,zm]
               if (this%umask(i,j,k).eq.0.and.this%umask(i,j-1,k).gt.0) this%itpu_y(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%umask(i,j,k).gt.0.and.this%umask(i,j-1,k).eq.0) this%itpu_y(:,i,j,k)=[0.0_WP,1.0_WP]
               ! Linear interpolation in y of W from [xm,ym,z ]
               if (this%wmask(i,j,k).eq.0.and.this%wmask(i,j-1,k).gt.0) this%itpw_y(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%wmask(i,j,k).gt.0.and.this%wmask(i,j-1,k).eq.0) this%itpw_y(:,i,j,k)=[0.0_WP,1.0_WP]
            end do
         end do
      end do
      
      ! Adjust interpolation coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               ! Linear interpolation in z of U from [x ,ym,zm]
               if (this%umask(i,j,k).eq.0.and.this%umask(i,j,k-1).gt.0) this%itpu_z(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%umask(i,j,k).gt.0.and.this%umask(i,j,k-1).eq.0) this%itpu_z(:,i,j,k)=[0.0_WP,1.0_WP]
               !  Linear interpolation in z of V from [xm,y ,zm]
               if (this%vmask(i,j,k).eq.0.and.this%vmask(i,j,k-1).gt.0) this%itpv_z(:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%vmask(i,j,k).gt.0.and.this%vmask(i,j,k-1).eq.0) this%itpv_z(:,i,j,k)=[0.0_WP,1.0_WP]
            end do
         end do
      end do
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the flow solver now that bconds have been defined
   subroutine setup(this,implicit_ils)
      implicit none
      class(scalar), intent(inout) :: this
      integer, intent(in) :: implicit_ils
      integer :: i,j,k,count,st
      
      ! Adjust metrics based on bcflag array
      call this%adjust_metrics()
      
      ! Set dynamic stencil map for the velocity solver
      count=1; this%implicit%stc(count,:)=[0,0,0]
      do st=1,+abs(this%stp1)
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
      call this%implicit%init_solver(implicit_ils)
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,dir,canCorrect,locator,inBalance)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(incomp), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      character(len=2), intent(in) :: dir
      logical,  intent(in) :: canCorrect
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
      logical, optional :: inBalance
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      new_bc%canCorrect=canCorrect
      select case (lowercase(dir))
      case ('+x','x+','xp','px')
         new_bc%dir=1
         new_bc%si=+1
         new_bc%sj= 0
         new_bc%sk= 0
      case ('-x','x-','xm','mx')
         new_bc%dir=2
         new_bc%si=-1
         new_bc%sj= 0
         new_bc%sk= 0
      case ('+y','y+','yp','py')
         new_bc%dir=3
         new_bc%si= 0
         new_bc%sj=+1
         new_bc%sk= 0
      case ('-y','y-','ym','my')
         new_bc%dir=4
         new_bc%si= 0
         new_bc%sj=-1
         new_bc%sk= 0
      case ('+z','z+','zp','pz')
         new_bc%dir=5
         new_bc%si= 0
         new_bc%sj= 0
         new_bc%sk=+1
      case ('-z','z-','zm','mz')
         new_bc%dir=6
         new_bc%si= 0
         new_bc%sj= 0
         new_bc%sk=-1
      case default; call die('[incomp add_bcond] Unknown bcond direction')
      end select
      new_bc%itr=iterator(this%cfg,new_bc%name,locator)
      if (present(inBalance)) then
         new_bc%inBalance=inBalance
      else
         new_bc%inBalance=.true.
      end if
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         ! Implement based on bcond direction
         select case (new_bc%dir)
         case (1) ! +x
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%umask(i+1,j,k)=2
            end do
         case (2) ! -x
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%umask(i,j,k)=2
            end do
         case (3) ! +y
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%vmask(i,j+1,k)=2
            end do
         case (4) ! -y
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%vmask(i,j,k)=2
            end do
         case (5) ! +z
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%wmask(i,j,k+1)=2
            end do
         case (6) ! -z
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%wmask(i,j,k)=2
            end do
         end select
         
      case (neumann)
         ! Not yet implemented
      case (clipped_neumann)
         ! Not yet implemented
      case (convective)
         ! Not yet implemented
      case default
         call die('[incomp apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      implicit none
      class(incomp), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(scalar), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n,ierr
      type(bcond), pointer :: my_bc
      real(WP) :: conv_vel,my_conv_vel
      
      ! ! First enfore zero velocity at walls
      ! do k=this%cfg%kmin_,this%cfg%kmax_
      !    do j=this%cfg%jmin_,this%cfg%jmax_
      !       do i=this%cfg%imin_,this%cfg%imax_
      !          if (minval(this%cfg%VF(i-1:i,j,k)).lt.10.0_WP*epsilon(1.0_WP)) this%U(i,j,k)=0.0_WP
      !          if (minval(this%cfg%VF(i,j-1:j,k)).lt.10.0_WP*epsilon(1.0_WP)) this%V(i,j,k)=0.0_WP
      !          if (minval(this%cfg%VF(i,j,k-1:k)).lt.10.0_WP*epsilon(1.0_WP)) this%W(i,j,k)=0.0_WP
      !       end do
      !    end do
      ! end do
      ! ! Sync fields
      ! call this%cfg%sync(this%U)
      ! call this%cfg%sync(this%V)
      ! call this%cfg%sync(this%W)
      
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
               select case (my_bc%dir)
               case (1) ! Neumann in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i+1,j,k)=this%U(i,j,k)
                     this%V(i+1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i+1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (2) ! Neumann in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i,j,k)=this%U(i+1,j,k)
                     this%V(i-1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i-1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (3) ! Neumann in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+1,k)=this%U(i:i+1,j,k)
                     this%V(i,j+1,k)=this%V(i,j,k)
                     this%W(i,j+1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (4) ! Neumann in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j-1,k)=this%U(i:i+1,j,k)
                     this%V(i,j,k)=this%V(i,j+1,k)
                     this%W(i,j-1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (5) ! Neumann in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k+1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k+1)=this%V(i,j:j+1,k)
                     this%W(i,j,k+1)=this%W(i,j,k)
                  end do
               case (6) ! Neumann in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k-1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k-1)=this%V(i,j:j+1,k)
                     this%W(i,j,k)=this%W(i,j,k+1)
                  end do
               end select
               
            case (clipped_neumann)     ! Apply clipped Neumann condition
               
               ! Implement based on bcond direction
               select case (my_bc%dir)
               case (1) ! Neumann in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i+1,j,k)=max(this%U(i,j,k),0.0_WP)
                     this%V(i+1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i+1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (2) ! Neumann in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i,j,k)=min(this%U(i+1,j,k),0.0_WP)
                     this%V(i-1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i-1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (3) ! Neumann in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+1,k)=this%U(i:i+1,j,k)
                     this%V(i,j+1,k)=max(this%V(i,j,k),0.0_WP)
                     this%W(i,j+1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (4) ! Neumann in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j-1,k)=this%U(i:i+1,j,k)
                     this%V(i,j,k)=min(this%V(i,j+1,k),0.0_WP)
                     this%W(i,j-1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (5) ! Neumann in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k+1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k+1)=this%V(i,j:j+1,k)
                     this%W(i,j,k+1)=max(this%W(i,j,k),0.0_WP)
                  end do
               case (6) ! Neumann in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k-1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k-1)=this%V(i,j:j+1,k)
                     this%W(i,j,k)=min(this%W(i,j,k+1),0.0_WP)
                  end do
               end select
            
            case (convective)   ! Apply convective condition
               
               ! Implement based on bcond direction
               select case (my_bc%dir)
               case (1) ! Neumann in +x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i+1,j,k)=this%U(i,j,k)
                     this%V(i+1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i+1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (2) ! Neumann in -x
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i,j,k)=this%U(i+1,j,k)
                     this%V(i-1,j:j+1,k)=this%V(i,j:j+1,k)
                     this%W(i-1,j,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (3) ! Neumann in +y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+1,k)=this%U(i:i+1,j,k)
                     this%V(i,j+1,k)=this%V(i,j,k)
                     this%W(i,j+1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (4) ! Neumann in -y
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j-1,k)=this%U(i:i+1,j,k)
                     this%V(i,j,k)=this%V(i,j+1,k)
                     this%W(i,j-1,k:k+1)=this%W(i,j,k:k+1)
                  end do
               case (5) ! Neumann in +z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k+1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k+1)=this%V(i,j:j+1,k)
                     this%W(i,j,k+1)=this%W(i,j,k)
                  end do
               case (6) ! Neumann in -z
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j,k-1)=this%U(i:i+1,j,k)
                     this%V(i,j:j+1,k-1)=this%V(i,j:j+1,k)
                     this%W(i,j,k)=this%W(i,j,k+1)
                  end do
               end select
               
               ! ! Get convective velocity
               ! my_conv_vel=0.0_WP
               ! do n=1,my_bc%itr%no_
               !    ! Get the indices
               !    i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
               !    ! Find maximum velocity
               !    my_conv_vel(1)=max(my_conv_vel(1),this%U(i+si,j+sj,k+sk))
               !    my_conv_vel(2)=max(my_conv_vel(2),this%V(i+si,j+sj,k+sk))
               !    my_conv_vel(3)=max(my_conv_vel(3),this%W(i+si,j+sj,k+sk))
               ! end do
               ! call MPI_ALLREDUCE(my_conv_vel,conv_vel,3,MPI_REAL_WP,MPI_MAX,my_bc%itr%comm,ierr)
               
               
            case default
               call die('[incomp apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%SC)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit rhoSC time derivative based on rhoU/rhoV/rhoW
   subroutine get_drhoSCdt(this,drhoSCdt,rhoU,rhoV,rhoW)
      implicit none
      class(scalar), intent(inout) :: this
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
               FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%itpsc_xp(:,i,j,k)*this%SC(i+stp1:i+stp2,j,k)) &
               &         -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%itpsc_xm(:,i,j,k)*this%SC(i+stm1:i+stm2,j,k)) &
               &         +this%diff*sum(this%grdsc_x(:,i,j,k)*this%SC(i-1:i,j,k))
               ! Fluxes on y-face
               FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%itpsc_yp(:,i,j,k)*this%SC(i,j+stp1:j+stp2,k)) &
               &         -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%itpsc_ym(:,i,j,k)*this%SC(i,j+stm1:j+stm2,k)) &
               &         +this%diff*sum(this%grdsc_y(:,i,j,k)*this%SC(i,j-1:j,k))
               ! Fluxes on z-face
               FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%itpsc_zp(:,i,j,k)*this%SC(i,j,k+stp1:k+stp2)) &
               &         -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%itpsc_zm(:,i,j,k)*this%SC(i,j,k+stm1:k+stm2)) &
               &         +this%diff*sum(this%grdsc_z(:,i,j,k)*this%SC(i,j,k-1:k))
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
   
   
   !> Calculate the min and max of our SC field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(scalar), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: my_SCmax,my_SCmin
      my_SCmax=maxval(this%SC); call MPI_ALLREDUCE(my_SCmax,this%SCmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_SCmin=minval(this%SC); call MPI_ALLREDUCE(my_SCmin,this%SCmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine get_max
   
   
   !> Solve for implicit scalar residual
   subroutine solve_implicit(this,dt,resSC,rho,rhoU,rhoV,rhoW)
      implicit none
      class(scalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resSC !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(in) :: rho
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rhoU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rhoV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rhoW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      
      ! Prepare operator
      this%implicit%opr(1,:,:,:)=rho; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%itpu_x( 0,i  ,j,k)*rhoUp+&
               &                                                                this%divu_x(-1,i,j,k)*this%itpu_x(+1,i-1,j,k)*rhoUm+&
               &                                                                this%divu_y(+1,i,j,k)*this%itpu_y(-1,i,j+1,k)*rhoVp+&
               &                                                                this%divu_y( 0,i,j,k)*this%itpu_y( 0,i,j  ,k)*rhoVm+&
               &                                                                this%divu_z(+1,i,j,k)*this%itpu_z(-1,i,j,k+1)*rhoWp+&
               &                                                                this%divu_z( 0,i,j,k)*this%itpu_z( 0,i,j,k  )*rhoWm)
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%itpu_x(+1,i  ,j,k)*rhoUp)
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+0.5_WP*dt*(this%divu_x(-1,i,j,k)*this%itpu_x( 0,i-1,j,k)*rhoUm)
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+0.5_WP*dt*(this%divu_y(+1,i,j,k)*this%itpu_y( 0,i,j+1,k)*rhoVp)
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+0.5_WP*dt*(this%divu_y( 0,i,j,k)*this%itpu_y(-1,i,j  ,k)*rhoVm)
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+0.5_WP*dt*(this%divu_z(+1,i,j,k)*this%itpu_z( 0,i,j,k+1)*rhoWp)
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+0.5_WP*dt*(this%divu_z( 0,i,j,k)*this%itpu_z(-1,i,j,k  )*rhoWm)
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%visc*this%grdu_x( 0,i  ,j,k)+&
               &                                                                this%divu_x(-1,i,j,k)*this%visc*this%grdu_x(+1,i-1,j,k)+&
               &                                                                this%divu_y(+1,i,j,k)*this%visc*this%grdu_y(-1,i,j+1,k)+&
               &                                                                this%divu_y( 0,i,j,k)*this%visc*this%grdu_y( 0,i,j  ,k)+&
               &                                                                this%divu_z(+1,i,j,k)*this%visc*this%grdu_z(-1,i,j,k+1)+&
               &                                                                this%divu_z( 0,i,j,k)*this%visc*this%grdu_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divu_x( 0,i,j,k)*this%visc*this%grdu_x(+1,i  ,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divu_x(-1,i,j,k)*this%visc*this%grdu_x( 0,i-1,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divu_y(+1,i,j,k)*this%visc*this%grdu_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divu_y( 0,i,j,k)*this%visc*this%grdu_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divu_z(+1,i,j,k)*this%visc*this%grdu_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divu_z( 0,i,j,k)*this%visc*this%grdu_z(-1,i,j,k  ))
            end do
         end do
      end do
      call this%implicit%update_solver()
      this%implicit%rhs=resU
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resU=this%implicit%sol
      
      ! Sync up residual
      call this%cfg%sync(resSC)
      
   end subroutine solve_implicit
   
   
   !> Print out info for scalar solver
   subroutine scalar_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(scalar), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Constant density scalar solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" > diffusivity = ",es12.5)') this%diff
      end if
      
   end subroutine scalar_print
   
   
end module scalar_class
