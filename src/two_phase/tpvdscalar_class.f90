!> Two-phase variable density scalar solver class:
!> Provides support for various BC, RHS calculation
!> Based on vfs geometric transport, hybridized with upwind for now.
module tpvdscalar_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   use linsol_class,   only: linsol
   implicit none
   private

   ! Expose type/constructor/methods
   public :: tpvdscalar,bcond

   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2                                    !< Dirichlet condition
   integer, parameter, public :: neumann=3                                      !< Zero normal gradient

   ! List of available advection schemes for scalar transport
   integer, parameter, public :: upwind=0                                       !< First order upwind scheme

   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                                              !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'                         !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                                           !< Bcond type
      integer :: dir                                                            !< Bcond direction (1 to 6)
      type(iterator) :: itr                                                     !< This is the iterator for the bcond
   end type bcond

   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))

   !> Two-phase scalar solver object definition
   type :: tpvdscalar

      ! This is our config
      class(config), pointer :: cfg                                             !< This is the config the solver is build for

      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_SCALAR'                        !< Solver name (default=UNNAMED_SCALAR)

      ! Variable property fluid
      integer, dimension(:), allocatable :: phase                               !< This is the phase for each scalar (0=liquid, 1=gas)
      real(WP), dimension(:,:,:,:), allocatable :: diff                         !< These is our constant+SGS dynamic diffusivity for the scalar

      ! Boundary condition list
      integer :: nbc                                                            !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                                          !< List of bcond for our solver

      ! Scalar variable
      integer :: nscalar                                                        !< Number of scalars
      character(len=str_medium), dimension(:), allocatable :: SCname            !< Names of scalars
      real(WP), dimension(:,:,:,:), allocatable :: rho                          !< Density array includes both phases: Liquid 0, Gas 1
      real(WP), dimension(:,:,:,:), allocatable :: SC                           !< SC array
      real(WP), dimension(:,:,:,:), allocatable :: rhoSC                        !< rhoSC array

      ! Old scalar variable
      real(WP), dimension(:,:,:,:), allocatable :: rhoold                       !< Old density array includes both phases: Liquid 0, Gas 1
      real(WP), dimension(:,:,:,:), allocatable :: SCold                        !< SCold array
      real(WP), dimension(:,:,:,:), allocatable :: rhoSCold                     !< rhoSCold array

      ! Implicit scalar solver
      class(linsol), pointer :: implicit                                        !< Iterative linear solver object for an implicit prediction of the scalar residual
      integer, dimension(:,:,:), allocatable :: stmap                           !< Inverse map from stencil shift to index location

      ! Metrics
      ! integer :: scheme                                                       !< Advection scheme for scalar
      ! integer :: nst                                                          !< Scheme order (and elemental stencil size)
      ! integer :: stp1,stp2                                                    !< Plus interpolation stencil extent for scalar advection
      ! integer :: stm1,stm2                                                    !< Minus interpolation stencil extent for scalar advection
      ! real(WP), dimension(:,:,:,:), allocatable :: itpsc_xp,itpsc_yp,itpsc_zp !< Plus interpolation for SC
      ! real(WP), dimension(:,:,:,:), allocatable :: itpsc_xm,itpsc_ym,itpsc_zm !< Minus interpolation for SC
      real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z            !< Divergence for SC
      real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z            !< Scalar gradient for SC
      real(WP), dimension(:,:,:,:), allocatable :: itp_x,itp_y,itp_z            !< Second order interpolation for SC diffusivity

      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask                            !< Integer array used for modifying SC metrics

      ! Monitoring quantities
      real(WP), dimension(:), allocatable :: SCmax,SCmin,SCint                  !< Maximum and minimum, integral scalar

   contains
      procedure :: initialize                                                   !< Initialization of the scalar solver
      procedure :: print=>tpvdscalar_print                                      !< Output solver to the screen
      procedure, private :: init_metrics                                        !< Initialize metrics
      procedure, private :: adjust_metrics                                      !< Adjust metrics
      procedure :: setup                                                        !< Finish configuring the scalar solver
      procedure :: add_bcond                                                    !< Add a boundary condition
      procedure :: get_bcond                                                    !< Get a boundary condition
      procedure :: apply_bcond                                                  !< Apply all boundary conditions
      procedure :: get_drhoSCdt                                                 !< Calculate drhoSC/dt
      procedure :: get_max                                                      !< Calculate maximum and integral field values
      procedure :: solve_implicit                                               !< Solve for the scalar residuals implicitly
   end type tpvdscalar


   !> Declare tpvdscalar solver constructor
   interface tpvdscalar
      procedure constructor
   end interface tpvdscalar

contains


   !> Default constructor for two-phase variable density scalar solver
   function constructor(self,cfg,nscalar,name) result(self)
      use messager, only: die
      implicit none
      type(tpvdscalar) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: nscalar
      character(len=*), optional :: name
      integer :: i,j,k

      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))

      ! Set the number of scalars
      self%nscalar=nscalar
      if (self%nscalar.le.0) call die('[tpvdscalar constructor] tpvdscalar object requires at least 1 scalar')

      ! Initialize scalar names
      allocate(self%SCname(1:self%nscalar))
      self%SCname='' ! User will set names

      ! Initialize scalar phase
      allocate(self%phase(1:self%nscalar))
      self%phase=0  ! User will set phase

      ! Point to pgrid object
      self%cfg=>cfg

      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()

      ! Allocate variables
      allocate(self%SC      (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,1:self%nscalar)); self%SC      =0.0_WP
      allocate(self%rho     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,0:1));            self%rho     =0.0_WP
      allocate(self%rhoold  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,0:1));            self%rhoold  =0.0_WP
      allocate(self%rhoSC   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,1:self%nscalar)); self%rhoSC   =0.0_WP
      allocate(self%rhoSCold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,1:self%nscalar)); self%rhoSCold=0.0_WP
      allocate(self%diff    (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,1:self%nscalar)); self%diff    =0.0_WP

      ! Check current overlap
      if (self%cfg%no.lt.1) call die('[tpvdscalar constructor] Scalar transport scheme requires larger overlap')

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

      ! Monitoring data needs to be allocated
      allocate(self%SCint(1:self%nscalar))
      allocate(self%SCmin(1:self%nscalar))
      allocate(self%SCmax(1:self%nscalar))

   end function constructor


   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      implicit none
      class(tpvdscalar), intent(inout) :: this
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

      ! Allocate finite difference scalar gradient operators
      allocate(this%grd_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%grd_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%grd_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%grd_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in x from [xm,ym,zm] to [x,ym,zm]
               this%grd_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in y from [xm,ym,zm] to [xm,y,zm]
               this%grd_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do

   end subroutine init_metrics


   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(tpvdscalar), intent(inout) :: this
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

      ! Loop over the domain and apply masked conditions to SC divergence
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

      ! Adjust gradient coefficients to cell faces for walls (assume Neumann at wall)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%mask(i,j,k).eq.1.or.this%mask(i-1,j,k).eq.1) this%grd_x(:,i,j,k)=0.0_WP     !< FD gradient in x of SC
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j-1,k).eq.1) this%grd_y(:,i,j,k)=0.0_WP     !< FD gradient in y of SC
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j,k-1).eq.1) this%grd_z(:,i,j,k)=0.0_WP     !< FD gradient in z of SC
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


   !> Finish setting up the scalar solver now that bconds have been defined
   subroutine setup(this,implicit_solver)
      use messager, only: die
      implicit none
      class(tpvdscalar), intent(inout) :: this
      class(linsol), target, intent(in), optional :: implicit_solver
      integer :: count

      ! Adjust metrics based on mask array
      call this%adjust_metrics()

      ! Prepare implicit solver if it had been provided
      if (present(implicit_solver)) then

         ! Point to implicit solver linsol object
         this%implicit=>implicit_solver

         ! Check implicit solver size
         if (this%implicit%nst.ne.7) call die('[tpvdscalar setup] Implicit solver needs nst=7')

         ! Set dynamic stencil map for the scalar solver - diffusion only
         count=      1; this%implicit%stc(count,:)=[ 0, 0, 0]
         count=count+1; this%implicit%stc(count,:)=[+1, 0, 0]
         count=count+1; this%implicit%stc(count,:)=[-1, 0, 0]
         count=count+1; this%implicit%stc(count,:)=[ 0,+1, 0]
         count=count+1; this%implicit%stc(count,:)=[ 0,-1, 0]
         count=count+1; this%implicit%stc(count,:)=[ 0, 0,+1]
         count=count+1; this%implicit%stc(count,:)=[ 0, 0,-1]

         ! Set the diagonal to 1 to make sure all cells participate in solver
         this%implicit%opr(1,:,:,:)=1.0_WP

         ! Initialize the implicit scalar solver
         call this%implicit%init()

      end if

   end subroutine setup


   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,         only: lowercase
      use messager,       only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(tpvdscalar), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
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
          case default; call die('[tpvdscalar add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[tpvdscalar apply_bcond] Neumann requires a direction')
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
         call die('[tpvdscalar apply_bcond] Unknown bcond type')
      end select

   end subroutine add_bcond


   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(tpvdscalar), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[tpvdscalar get_bcond] Boundary condition was not found')
   end subroutine get_bcond


   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpvdscalar), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n,nsc
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
               do nsc=1,this%nscalar
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%SC(i,j,k,nsc)=this%SC(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir),nsc)
                  end do
               end do

             case default
               call die('[tpvdscalar apply_bcond] Unknown bcond type')
            end select

         end if

         ! Move on to the next bcond
         my_bc=>my_bc%next

      end do

      ! Sync full fields after all bcond
      do nsc=1,this%nscalar
         call this%cfg%sync(this%SC(:,:,:,nsc))
      end do

   end subroutine apply_bcond


   !> Calculate the explicit rhoSC time derivative based on U/V/W
   subroutine get_drhoSCdt(this,drhoSCdt,U,V,W,VFold,VF,detailed_face_flux,dt)
      use irl_fortran_interface
      implicit none
      class(tpvdscalar), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: drhoSCdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)   , intent(in)  :: U        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)   , intent(in)  :: V        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)   , intent(in)  :: W        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)   , intent(in)  :: VFold    !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)   , intent(in)  :: VF       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      type(TagAccVM_SepVM_type), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:) :: detailed_face_flux !< Needs to be (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(in) :: dt  !< This is the time step size that was used to generate the detailed_face_flux geometric data
      type(SepVM_type) :: my_SepVM
      integer :: i,j,k,nsc,n
      real(WP), dimension(:,:,:),   allocatable :: FX,FY,FZ
      real(WP), dimension(:,:,:,:), allocatable :: grad
      real(WP) :: my_vol,SCm,SCp
      !real(WP), dimension(3) :: my_bar
      integer, dimension(3) :: ind
      ! Zero out dSC/dt array
      drhoSCdt=0.0_WP
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Allocate scalar gradient
      allocate(grad(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Work on each scalar
      do nsc=1,this%nscalar
         ! Reset fluxes and gradient to zero
         FX=0.0_WP; FY=0.0_WP; FZ=0.0_WP; grad=0.0_WP
         ! Calculate minmod-limited gradient of SC everywhere
         do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1
            do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1
               do i=this%cfg%imino_+1,this%cfg%imaxo_-1
                  ! No need to calculate gradient inside of wall cell
                  if (this%mask(i,j,k).eq.1) cycle
                  ! Get gradient
                  grad(1,i,j,k)=minmod((this%SC(i+1,j,k,nsc)-this%SC(i,j,k,nsc))*this%cfg%dxmi(i+1),(this%SC(i,j,k,nsc)-this%SC(i-1,j,k,nsc))*this%cfg%dxmi(i))
                  grad(2,i,j,k)=minmod((this%SC(i,j+1,k,nsc)-this%SC(i,j,k,nsc))*this%cfg%dymi(j+1),(this%SC(i,j,k,nsc)-this%SC(i,j-1,k,nsc))*this%cfg%dymi(j))
                  grad(3,i,j,k)=minmod((this%SC(i,j,k+1,nsc)-this%SC(i,j,k,nsc))*this%cfg%dzmi(k+1),(this%SC(i,j,k,nsc)-this%SC(i,j,k-1,nsc))*this%cfg%dzmi(k))
               end do
            end do
         end do
         call this%cfg%sync(grad)
         ! Convective flux of SC
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Flux on x-face
                  if (getSize(detailed_face_flux(1,i,j,k)).gt.0) then
                     ! Detailed geometric flux is available, use geometric fluxing
                     do n=0,getSize(detailed_face_flux(1,i,j,k))-1
                        ! Get cell index for nth object
                        ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux(1,i,j,k),n))
                        ! Get SepVM for nth object
                        call getSepVMAtIndex(detailed_face_flux(1,i,j,k),n,my_SepVM)
                        ! Extract volume for relevant phase
                        my_vol=getVolume(my_SepVM,this%phase(nsc))
                        ! Increment flux with first order estimate
                        FX(i,j,k)=FX(i,j,k)-my_vol*this%SCold(ind(1),ind(2),ind(3),nsc)
                        ! Second order correction
                        !my_bar=getCentroid(my_SepVM,this%phase(nsc))
                        !FX(i,j,k)=FX(i,j,k)-my_vol*(sum(grad(:,ii,jj,kk)*my_bar(:)-my_barold(:)))
                     end do
                     ! Scale by cell face area and time step size
                     FX(i,j,k)=FX(i,j,k)/(dt*this%cfg%dy(j)*this%cfg%dz(k))
                  else
                     ! No detailed geometric flux is available, use MUSCL flux
                     SCm=0.0_WP; if (VFold(i-1,j,k).ne.real(this%phase(nsc),WP)) SCm=this%SC(i-1,j,k,nsc)+0.5_WP*grad(1,i-1,j,k)*this%cfg%dx(i-1)
                     SCp=0.0_WP; if (VFold(i  ,j,k).ne.real(this%phase(nsc),WP)) SCp=this%SC(i  ,j,k,nsc)-0.5_WP*grad(1,i  ,j,k)*this%cfg%dx(i  )
                     FX(i,j,k)=-0.5_WP*(U(i,j,k)+abs(U(i,j,k)))*SCm-0.5_WP*(U(i,j,k)-abs(U(i,j,k)))*SCp
                  end if
                  ! Flux on y-face
                  if (getSize(detailed_face_flux(2,i,j,k)).gt.0) then
                     ! Detailed geometric flux is available, use geometric fluxing
                     do n=0,getSize(detailed_face_flux(2,i,j,k))-1
                        ! Get cell index for nth object
                        ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux(2,i,j,k),n))
                        ! Get SepVM for nth object
                        call getSepVMAtIndex(detailed_face_flux(2,i,j,k),n,my_SepVM)
                        ! Extract volume for relevant phase
                        my_vol=getVolume(my_SepVM,this%phase(nsc))
                        ! Increment flux with first order estimate
                        FY(i,j,k)=FY(i,j,k)-my_vol*this%SCold(ind(1),ind(2),ind(3),nsc)
                        ! Second order correction
                        !my_bar=getCentroid(my_SepVM,this%phase(nsc))
                        !FY(i,j,k)=FY(i,j,k)-my_vol*(sum(grad(:,ii,jj,kk)*my_bar(:)-my_barold(:)))
                     end do
                     ! Scale by cell face area and time step size
                     FY(i,j,k)=FY(i,j,k)/(dt*this%cfg%dx(i)*this%cfg%dz(k))
                  else
                     ! No detailed geometric flux is available, use MUSCL flux
                     SCm=0.0_WP; if (VFold(i,j-1,k).ne.real(this%phase(nsc),WP)) SCm=this%SC(i,j-1,k,nsc)+0.5_WP*grad(2,i,j-1,k)*this%cfg%dy(j-1)
                     SCp=0.0_WP; if (VFold(i,j  ,k).ne.real(this%phase(nsc),WP)) SCp=this%SC(i,j  ,k,nsc)-0.5_WP*grad(2,i,j  ,k)*this%cfg%dy(j  )
                     FY(i,j,k)=-0.5_WP*(V(i,j,k)+abs(V(i,j,k)))*SCm-0.5_WP*(V(i,j,k)-abs(V(i,j,k)))*SCp
                  end if
                  ! Flux on z-face
                  if (getSize(detailed_face_flux(3,i,j,k)).gt.0) then
                     ! Detailed geometric flux is available, use geometric fluxing
                     do n=0,getSize(detailed_face_flux(3,i,j,k))-1
                        ! Get cell index for nth object
                        ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux(3,i,j,k),n))
                        ! Get SepVM for nth object
                        call getSepVMAtIndex(detailed_face_flux(3,i,j,k),n,my_SepVM)
                        ! Extract volume for relevant phase
                        my_vol=getVolume(my_SepVM,this%phase(nsc))
                        ! Increment flux with first order estimate
                        FZ(i,j,k)=FZ(i,j,k)-my_vol*this%SCold(ind(1),ind(2),ind(3),nsc)
                        ! Second order correction
                        !my_bar=getCentroid(my_SepVM,this%phase(nsc))
                        !FZ(i,j,k)=FZ(i,j,k)-my_vol*(sum(grad(:,ii,jj,kk)*my_bar(:)-my_barold(:)))
                     end do
                     ! Scale by cell face area and time step size
                     FZ(i,j,k)=FZ(i,j,k)/(dt*this%cfg%dx(i)*this%cfg%dy(j))
                  else
                     ! No detailed geometric flux is available, use MUSCL flux
                     SCm=0.0_WP; if (VFold(i,j,k-1).ne.real(this%phase(nsc),WP)) SCm=this%SC(i,j,k-1,nsc)+0.5_WP*grad(3,i,j,k-1)*this%cfg%dz(k-1)
                     SCp=0.0_WP; if (VFold(i,j,k  ).ne.real(this%phase(nsc),WP)) SCp=this%SC(i,j,k  ,nsc)-0.5_WP*grad(3,i,j,k  )*this%cfg%dz(k  )
                     FZ(i,j,k)=-0.5_WP*(W(i,j,k)+abs(W(i,j,k)))*SCm-0.5_WP*(W(i,j,k)-abs(W(i,j,k)))*SCp
                  end if
               end do
            end do
         end do
         ! Diffusive flux of SC - needs to be made one-sided
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  FX(i,j,k)=FX(i,j,k)+sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k,nsc))*sum(this%grd_x(:,i,j,k)*this%SC(i-1:i,j,k,nsc))
                  FY(i,j,k)=FY(i,j,k)+sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k,nsc))*sum(this%grd_y(:,i,j,k)*this%SC(i,j-1:j,k,nsc))
                  FZ(i,j,k)=FZ(i,j,k)+sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k,nsc))*sum(this%grd_z(:,i,j,k)*this%SC(i,j,k-1:k,nsc))
               end do
            end do
         end do
         ! Time derivative of SC
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  drhoSCdt(i,j,k,nsc)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+&
                  &                sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+&
                  &                sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
               end do
            end do
         end do
         ! Sync residual
         call this%cfg%sync(drhoSCdt(:,:,:,nsc))
      end do
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ,grad)

   contains

      !> Minmod gradient
      function minmod(g1,g2) result(g)
         implicit none
         real(WP), intent(in) :: g1,g2
         real(WP) :: g
         if (g1*g2.le.0.0_WP) then
            g=0.0_WP
         else
            if (abs(g1).lt.abs(g2)) then
               g=g1
            else
               g=g2
            end if
         end if
      end function minmod

   end subroutine get_drhoSCdt


   !> Calculate the time derivative of rho
   subroutine get_drhodt(this,dt,drhodt)
      implicit none
      class(tpvdscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhodt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      drhodt=(this%rho-this%rhoold)/dt
   end subroutine get_drhodt


   !> Calculate the min, max, and int of our SC field
   subroutine get_max(this,VF)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpvdscalar), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: VF        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: ierr,nsc
      real(WP) :: my_SCmax,my_SCmin
      real(WP), dimension(:,:,:), allocatable :: tmp
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      do nsc=1,this%nscalar
         my_SCmax=maxval(this%SC(:,:,:,nsc)); call MPI_ALLREDUCE(my_SCmax,this%SCmax(nsc),1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         my_SCmin=minval(this%SC(:,:,:,nsc)); call MPI_ALLREDUCE(my_SCmin,this%SCmin(nsc),1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         if      (this%phase(nsc).eq.0) then ! Liquid scalar
            tmp=this%SC(:,:,:,nsc)*(       VF(:,:,:))
         else if (this%phase(nsc).eq.1) then ! Gas scalar
            tmp=this%SC(:,:,:,nsc)*(1.0_WP-VF(:,:,:))
         end if
         call this%cfg%integrate(A=tmp,integral=this%SCint(nsc))
      end do
      deallocate(tmp)
   end subroutine get_max


   !> Divide rhoSC by rho to form SC
   subroutine rho_divide(this)
      implicit none
      class(tpvdscalar), intent(inout) :: this
      where (this%mask.eq.0) this%SC=this%rhoSC/this%rho
   end subroutine rho_divide


   !> Multiply SC by rho to form rhoSC
   subroutine rho_multiply(this)
      implicit none
      class(tpvdscalar), intent(inout) :: this
      integer :: i,j,k
      where (this%mask.eq.0) this%rhoSC=this%SC*this%rho
   end subroutine rho_multiply


   !> Solve for implicit scalar residual
   subroutine solve_implicit(this,dt,resSC)
      implicit none
      class(tpvdscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar)
      integer :: i,j,k,nsc

      ! Apply implicit treatment for each scalar
      do nsc=1,this%nscalar

         ! Prepare diffusive operator
         this%implicit%opr(1,:,:,:)=1.0_WP; this%implicit%opr(2:,:,:,:)=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%div_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k,nsc))*this%grd_x(-1,i+1,j,k)+&
                  &                                                                this%div_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k,nsc))*this%grd_x( 0,i  ,j,k)+&
                  &                                                                this%div_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k,nsc))*this%grd_y(-1,i,j+1,k)+&
                  &                                                                this%div_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k,nsc))*this%grd_y( 0,i,j  ,k)+&
                  &                                                                this%div_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1,nsc))*this%grd_z(-1,i,j,k+1)+&
                  &                                                                this%div_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ,nsc))*this%grd_z( 0,i,j,k  ))
                  this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%div_x(+1,i,j,k)*sum(this%itp_x(:,i+1,j,k)*this%diff(i  :i+1,j,k,nsc))*this%grd_x( 0,i+1,j,k))
                  this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%div_x( 0,i,j,k)*sum(this%itp_x(:,i  ,j,k)*this%diff(i-1:i  ,j,k,nsc))*this%grd_x(-1,i  ,j,k))
                  this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%div_y(+1,i,j,k)*sum(this%itp_y(:,i,j+1,k)*this%diff(i,j  :j+1,k,nsc))*this%grd_y( 0,i,j+1,k))
                  this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%div_y( 0,i,j,k)*sum(this%itp_y(:,i,j  ,k)*this%diff(i,j-1:j  ,k,nsc))*this%grd_y(-1,i,j  ,k))
                  this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%div_z(+1,i,j,k)*sum(this%itp_z(:,i,j,k+1)*this%diff(i,j,k  :k+1,nsc))*this%grd_z( 0,i,j,k+1))
                  this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%div_z( 0,i,j,k)*sum(this%itp_z(:,i,j,k  )*this%diff(i,j,k-1:k  ,nsc))*this%grd_z(-1,i,j,k  ))
               end do
            end do
         end do

         ! Solve the linear system
         call this%implicit%setup()
         this%implicit%rhs=resSC(:,:,:,nsc)
         this%implicit%sol=0.0_WP
         call this%implicit%solve()
         resSC(:,:,:,nsc)=this%implicit%sol

         ! Sync up residual
         call this%cfg%sync(resSC(:,:,:,nsc))

      end do

   end subroutine solve_implicit


   !> Print out info for tpvdscalar solver
   subroutine tpvdscalar_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(tpvdscalar), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("Two-phase variable density scalar solver [",a,"] with [",i3,"] scalars for config [",a,"]")') trim(this%name),this%nscalar,trim(this%cfg%name)
      end if
   end subroutine tpvdscalar_print


end module tpvdscalar_class
