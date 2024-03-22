!> Two-phase incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density in each phase.
!> Interface is represented using VOF
module tpns_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use linsol_class,   only: linsol
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: tpns,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: wall=1              !< Dirichlet at zero condition
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   integer, parameter, public :: convective=4        !< Convective outflow condition
   integer, parameter, public :: clipped_neumann=5   !< Clipped Neumann condition (outflow only)
   integer, parameter, public :: slip=6              !< Free-slip condition

   ! List of available contact line models for this solver
   integer, parameter, public :: static_contact=1    !< Static contact line model
   
   ! List of available averaging strategies for viscosity
   integer, parameter, public :: harmonic_visc=1     !< Harmonically-averaged viscosity
   integer, parameter, public :: arithmetic_visc=2   !< Arithmetically-averaged viscosity
   
   ! Parameter for switching schemes around the interface
   real(WP), parameter :: rhoeps_coeff=1.0e-3_WP     !< Parameter for deciding when to switch to upwinded transport
   
   !> Boundary conditions for the two-phase solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      type(iterator) :: itr                               !< This is the iterator for the bcond - this identifies the (i,j,k)
      character(len=1) :: face                            !< Bcond face (x/y/z)
      integer :: dir                                      !< Bcond direction (+1,-1,0 for interior)
      real(WP) :: rdir                                    !< Bcond direction (real variable)
      logical :: canCorrect                               !< Can this bcond be corrected for global conservation?
   end type bcond
   
   !> Two-phase incompressible solver object definition
   type :: tpns
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_TPNS'    !< Solver name (default=UNNAMED_TPNS)
      
      ! Constant property fluids
      real(WP) :: contact_angle                           !< This is our static contact angle
      real(WP) :: sigma                                   !< This is our constant surface tension coefficient
      real(WP) :: rho_l,rho_g                             !< These are our constant densities in liquid and gas
      real(WP) :: visc_l,visc_g                           !< These is our constant dynamic viscosities in liquid and gas
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP            !< Acceleration of gravity
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      real(WP), dimension(:), allocatable :: mfr          !< MFR through each bcond
      real(WP), dimension(:), allocatable :: area         !< Area for each bcond
      real(WP) :: correctable_area                        !< Area of bcond that can be corrected
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Density and face density fields
      real(WP), dimension(:,:,:), allocatable :: rho      !< New density field
      real(WP), dimension(:,:,:), allocatable :: rho_U    !< New density field at U face
      real(WP), dimension(:,:,:), allocatable :: rho_V    !< New density field at V face
      real(WP), dimension(:,:,:), allocatable :: rho_W    !< New density field at W face
      
      ! Viscosity fields
      real(WP), dimension(:,:,:), allocatable :: visc     !< Viscosity field on P-cell
      real(WP), dimension(:,:,:), allocatable :: visc_xy  !< Viscosity field on U-cell
      real(WP), dimension(:,:,:), allocatable :: visc_yz  !< Viscosity field on V-cell
      real(WP), dimension(:,:,:), allocatable :: visc_zx  !< Viscosity field on W-cell
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: rhoU     !< U momentum array
      real(WP), dimension(:,:,:), allocatable :: rhoV     !< V momentum array
      real(WP), dimension(:,:,:), allocatable :: rhoW     !< W momentum array
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      real(WP), dimension(:,:,:), allocatable :: Pjx      !< Pressure jump to add to -dP/dx
      real(WP), dimension(:,:,:), allocatable :: Pjy      !< Pressure jump to add to -dP/dy
      real(WP), dimension(:,:,:), allocatable :: Pjz      !< Pressure jump to add to -dP/dz
      real(WP), dimension(:,:,:), allocatable :: dPjx     !< dPressure jump to add to -ddP/dx
      real(WP), dimension(:,:,:), allocatable :: dPjy     !< dPressure jump to add to -ddP/dy
      real(WP), dimension(:,:,:), allocatable :: dPjz     !< dPressure jump to add to -ddP/dz
      
      ! Old flow variables
      real(WP), dimension(:,:,:), allocatable :: rhoold   !< Old density field
      real(WP), dimension(:,:,:), allocatable :: rhoUold  !< rhoUold momentum array
      real(WP), dimension(:,:,:), allocatable :: rhoVold  !< rhoVold momentum array
      real(WP), dimension(:,:,:), allocatable :: rhoWold  !< rhoWold momentum array
      real(WP), dimension(:,:,:), allocatable :: Uold     !< Uold velocity array
      real(WP), dimension(:,:,:), allocatable :: Vold     !< Vold velocity array
      real(WP), dimension(:,:,:), allocatable :: Wold     !< Wold velocity array
      
      ! Flow divergence
      real(WP), dimension(:,:,:), allocatable :: div      !< Divergence array
      
      ! Pressure solver
      class(linsol), pointer :: psolv                     !< Iterative linear solver object for the pressure Poisson equation
      
      ! Implicit velocity solver
      class(linsol), pointer :: implicit                  !< Iterative linear solver object for an implicit prediction of the NS residual
      
      ! Metrics
      real(WP), dimension(:,:,:,:,:), allocatable :: itp_xy,itp_yz,itp_xz !< Interpolation for viscosity
      real(WP), dimension(:,:,:,:), allocatable :: itpr_x,itpr_y,itpr_z   !< Interpolation for density
      real(WP), dimension(:,:,:,:), allocatable :: itpu_x,itpu_y,itpu_z   !< Second order interpolation for U
      real(WP), dimension(:,:,:,:), allocatable :: itpv_x,itpv_y,itpv_z   !< Second order interpolation for V
      real(WP), dimension(:,:,:,:), allocatable :: itpw_x,itpw_y,itpw_z   !< Second order interpolation for W
      real(WP), dimension(:,:,:,:), allocatable :: hybu_x,hybu_y,hybu_z   !< Hybrid interpolation for U
      real(WP), dimension(:,:,:,:), allocatable :: hybv_x,hybv_y,hybv_z   !< Hybrid interpolation for V
      real(WP), dimension(:,:,:,:), allocatable :: hybw_x,hybw_y,hybw_z   !< Hybrid interpolation for W
      real(WP), dimension(:,:,:,:), allocatable :: divp_x,divp_y,divp_z   !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: divu_x,divu_y,divu_z   !< Divergence for U-cell
      real(WP), dimension(:,:,:,:), allocatable :: divv_x,divv_y,divv_z   !< Divergence for V-cell
      real(WP), dimension(:,:,:,:), allocatable :: divw_x,divw_y,divw_z   !< Divergence for W-cell
      real(WP), dimension(:,:,:,:), allocatable :: grdu_x,grdu_y,grdu_z   !< Velocity gradient for U
      real(WP), dimension(:,:,:,:), allocatable :: grdv_x,grdv_y,grdv_z   !< Velocity gradient for V
      real(WP), dimension(:,:,:,:), allocatable :: grdw_x,grdw_y,grdw_z   !< Velocity gradient for W
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable ::  mask                     !< Integer array used for modifying P metrics
      integer, dimension(:,:,:), allocatable :: umask                     !< Integer array used for modifying U metrics
      integer, dimension(:,:,:), allocatable :: vmask                     !< Integer array used for modifying V metrics
      integer, dimension(:,:,:), allocatable :: wmask                     !< Integer array used for modifying W metrics
      
      ! CFL numbers
      real(WP) :: CFLst                                                   !< Surface tension CFL
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                                    !< Convective CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                                    !< Viscous CFL numbers
      
      ! Monitoring quantities
      real(WP) :: Umax,Vmax,Wmax,Pmax,divmax                              !< Maximum velocity, pressure, divergence
      
   contains
      procedure :: print=>tpns_print                      !< Output solver to the screen
      procedure :: initialize                             !< Initialize the flow solver
      procedure :: setup                                  !< Finish configuring the flow solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: get_dmomdt                             !< Calculate dmom/dt
      procedure :: get_div                                !< Calculate velocity divergence
      procedure :: get_pgrad                              !< Calculate pressure gradient
      procedure :: update_laplacian                       !< Update the pressure Laplacian div(1/rho*grad(.))
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_strainrate                         !< Calculate deviatoric part of strain rate tensor
      procedure :: get_gradu                              !< Calculate velocity gradient tensor
      procedure :: get_vorticity                          !< Calculate vorticity tensor
      procedure :: get_mfr                                !< Calculate outgoing MFR through each bcond
      procedure :: correct_mfr                            !< Correct for mfr mismatch to ensure global conservation
      procedure :: shift_p                                !< Shift pressure to have zero average
      procedure :: update_density                         !< Calculate new density and momentum from vfs object
      procedure :: get_face_density                       !< Calculate face densities by interpolation
      procedure :: get_viscosity                          !< Calculate viscosity fields from subcell phasic volume data in a vfs object
      procedure :: solve_implicit                         !< Solve for the velocity residuals implicitly
      
      procedure :: addsrc_gravity                         !< Add gravitational body force
      procedure :: add_surface_tension_jump               !< Add surface tension jump
      procedure :: add_static_contact                     !< Add static contact line model to surface tension jump
      
   end type tpns
   

contains
   
   
   !> Initialization for two-phase incompressible flow solver
   subroutine initialize(this,cfg,name)
      implicit none
      class(tpns), intent(inout) :: this
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
      
      ! Allocate flow variables
      allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      allocate(this%P(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%P=0.0_WP
      allocate(this%Pjx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Pjx=0.0_WP
      allocate(this%Pjy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Pjy=0.0_WP
      allocate(this%Pjz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Pjz=0.0_WP
      allocate(this%dPjx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%dPjx=0.0_WP
      allocate(this%dPjy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%dPjy=0.0_WP
      allocate(this%dPjz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%dPjz=0.0_WP
      allocate(this%visc   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc   =0.0_WP
      allocate(this%visc_xy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_xy=0.0_WP
      allocate(this%visc_yz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_yz=0.0_WP
      allocate(this%visc_zx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_zx=0.0_WP
      
      ! Mass conservation data around which to build momentum conservation
      allocate(this%rho    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rho    =0.0_WP
      allocate(this%rho_U  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rho_U  =0.0_WP
      allocate(this%rho_V  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rho_V  =0.0_WP
      allocate(this%rho_W  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rho_W  =0.0_WP
      allocate(this%rhoold (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoold =0.0_WP
      allocate(this%rhoU   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoU   =0.0_WP
      allocate(this%rhoV   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoV   =0.0_WP
      allocate(this%rhoW   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoW   =0.0_WP
      allocate(this%rhoUold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoUold=0.0_WP
      allocate(this%rhoVold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoVold=0.0_WP
      allocate(this%rhoWold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoWold=0.0_WP
      
      ! Allocate flow divergence
      allocate(this%div(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%div=0.0_WP
      
      ! Allocate old flow variables
      allocate(this%Uold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Uold=0.0_WP
      allocate(this%Vold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Vold=0.0_WP
      allocate(this%Wold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Wold=0.0_WP
      
      ! Prepare default metrics
      call this%init_metrics()
      
      ! Prepare P-cell masks
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
      
      ! Prepare face mask for U
      allocate(this%umask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%umask=0
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.           1) this%umask(this%cfg%imin  ,:,:)=2
         if (this%cfg%iproc.eq.this%cfg%npx) this%umask(this%cfg%imax+1,:,:)=2
      end if
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (minval(this%cfg%VF(i-1:i,j,k)).eq.0.0_WP) this%umask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%umask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      
      ! Prepare face mask for V
      allocate(this%vmask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%vmask=0
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.           1) this%vmask(:,this%cfg%jmin  ,:)=2
         if (this%cfg%jproc.eq.this%cfg%npy) this%vmask(:,this%cfg%jmax+1,:)=2
      end if
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (minval(this%cfg%VF(i,j-1:j,k)).eq.0.0_WP) this%vmask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%vmask)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      
      ! Prepare face mask for W
      allocate(this%wmask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%wmask=0
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.           1) this%wmask(:,:,this%cfg%kmin  )=2
         if (this%cfg%kproc.eq.this%cfg%npz) this%wmask(:,:,this%cfg%kmax+1)=2
      end if
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (minval(this%cfg%VF(i,j,k-1:k)).eq.0.0_WP) this%wmask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
   end subroutine initialize
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP), dimension(-1:0) :: itpx,itpy,itpz
      
      ! Allocate hybrid interpolation coefficients - these are populated later
      allocate(this%hybu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybu_x=0.0_WP
      allocate(this%hybv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybv_y=0.0_WP
      allocate(this%hybw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybw_z=0.0_WP
      allocate(this%hybv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybv_x=0.0_WP
      allocate(this%hybw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybw_x=0.0_WP
      allocate(this%hybu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybu_y=0.0_WP
      allocate(this%hybw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)); this%hybw_y=0.0_WP
      allocate(this%hybu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)); this%hybu_z=0.0_WP
      allocate(this%hybv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)); this%hybv_z=0.0_WP
      
      ! Allocate finite difference density interpolation coefficients to cell faces
      allocate(this%itpr_x(-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< X-face-centered
      allocate(this%itpr_y(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Y-face-centered
      allocate(this%itpr_z(-1:0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Z-face-centered
      ! Create density interpolation coefficients to cell face in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpr_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
            end do
         end do
      end do
      ! Create density interpolation coefficients to cell face in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpr_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
            end do
         end do
      end do
      ! Create density interpolation coefficients to cell face in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpr_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Allocate finite difference viscosity interpolation coefficients
      allocate(this%itp_xy(-1:0,-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%itp_yz(-1:0,-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (yz)
      allocate(this%itp_xz(-1:0,-1:0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (zx)
      ! Create viscosity interpolation coefficients to cell edge
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! Prepare local 1D metrics
               itpx=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)]
               itpy=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)]
               itpz=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)]
               ! Combine for 2D interpolations
               do st1=-1,0
                  do st2=-1,0
                     this%itp_xy(st1,st2,i,j,k)=itpx(st1)*itpy(st2)
                     this%itp_yz(st1,st2,i,j,k)=itpy(st1)*itpz(st2)
                     this%itp_xz(st1,st2,i,j,k)=itpx(st1)*itpz(st2)
                  end do
               end do
            end do
         end do
      end do
      
      ! Allocate finite difference velocity interpolation coefficients
      allocate(this%itpu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%itpv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%itpw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%itpv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%itpw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%itpu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%itpw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (yz)
      allocate(this%itpu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%itpv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (yz)
      ! Create velocity interpolation coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%itpu_x(:,i,j,k)=[+0.5_WP,+0.5_WP] !< Linear interpolation in x of U from [x ,ym,zm]
               this%itpv_y(:,i,j,k)=[+0.5_WP,+0.5_WP] !< Linear interpolation in y of V from [xm,y ,zm]
               this%itpw_z(:,i,j,k)=[+0.5_WP,+0.5_WP] !< Linear interpolation in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpv_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x of V from [xm,y ,zm]
               this%itpw_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpu_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y of U from [x ,ym,zm]
               this%itpw_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpu_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z of U from [x ,ym,zm]
               this%itpv_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z of V from [xm,y ,zm]
            end do
         end do
      end do
      
      ! Allocate finite volume divergence operators
      allocate(this%divp_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divu_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divu_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divu_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divv_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divv_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divv_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divw_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (z)
      allocate(this%divw_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (z)
      allocate(this%divw_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Face-centered (z)
      ! Create divergence operator to cell center [xm,ym,zm] or tangent to cell face
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%divp_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divp_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divp_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               
               this%divu_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,y ,zm]
               this%divu_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,z ]
               
               this%divv_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,y ,zm]
               this%divv_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,z ]
               
               this%divw_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,z ]
               this%divw_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,z ]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [x ,ym,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%divu_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [xm,y ,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divv_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      ! Create divergence operator perpendicular to cell face [xm,ym,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divw_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(this%grdu_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdv_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdw_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%grdv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%grdw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%grdu_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (xy)
      allocate(this%grdw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Edge-centered (yz)
      allocate(this%grdu_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (zx)
      allocate(this%grdv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Edge-centered (yz)
      ! Create gradient coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%grdu_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of U from [x ,ym,zm]
               this%grdv_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of V from [xm,y ,zm]
               this%grdw_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in x
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%grdv_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of V from [xm,y ,zm]
               this%grdw_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in y
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdu_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of U from [x ,ym,zm]
               this%grdw_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge in z
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdu_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of U from [x ,ym,zm]
               this%grdv_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z of V from [xm,y ,zm]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls
   subroutine adjust_metrics(this)
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP) :: delta,mysum
      
      ! Sync up u/v/wmasks
      call this%cfg%sync(this%umask)
      call this%cfg%sync(this%vmask)
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
      ! I am assuming here that we do not really need to zero out wall cells
      ! as they could be used for Dirichlet (then the density needs to be available! could be problematic if we do not have an explicit BC for scalars, e.g. for a Couette flow)
      ! or outflow condition (then the density needs to be available but it should be directly calculated)
      ! or used for a real no-slip wall (then density is always multiplied by zero)
      ! Adjust density interpolation coefficients to cell faces in the presence of walls (only walls!)
      !do k=this%cfg%kmin_,this%cfg%kmax_+1
      !   do j=this%cfg%jmin_,this%cfg%jmax_+1
      !      do i=this%cfg%imin_,this%cfg%imax_+1
      !         ! Linear interpolation in x
      !         if (this%cfg%VF(i,j,k).eq.0.0_WP.and.this%cfg%VF(i-1,j,k).gt.0.0_WP) this%itpr_x(:,i,j,k)=[1.0_WP,0.0_WP]
      !         if (this%cfg%VF(i,j,k).gt.0.0_WP.and.this%cfg%VF(i-1,j,k).eq.0.0_WP) this%itpr_x(:,i,j,k)=[0.0_WP,1.0_WP]
      !         ! Linear interpolation in y
      !         if (this%cfg%VF(i,j,k).eq.0.0_WP.and.this%cfg%VF(i,j-1,k).gt.0.0_WP) this%itpr_y(:,i,j,k)=[1.0_WP,0.0_WP]
      !         if (this%cfg%VF(i,j,k).gt.0.0_WP.and.this%cfg%VF(i,j-1,k).eq.0.0_WP) this%itpr_y(:,i,j,k)=[0.0_WP,1.0_WP]
      !         ! Linear interpolation in z
      !         if (this%cfg%VF(i,j,k).eq.0.0_WP.and.this%cfg%VF(i,j,k-1).gt.0.0_WP) this%itpr_z(:,i,j,k)=[1.0_WP,0.0_WP]
      !         if (this%cfg%VF(i,j,k).gt.0.0_WP.and.this%cfg%VF(i,j,k-1).eq.0.0_WP) this%itpr_z(:,i,j,k)=[0.0_WP,1.0_WP]
      !      end do
      !   end do
      !end do
      
      ! Adjust interpolation coefficients to cell centers in the presence of walls (only walls!)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) this%itpu_x(:,i,j,k)=0.0_WP
               if (this%mask(i,j,k).eq.1) this%itpv_y(:,i,j,k)=0.0_WP
               if (this%mask(i,j,k).eq.1) this%itpw_z(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Adjust viscosity interpolation coefficients to cell edge in the presence of walls (only walls)
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               ! Zero out interpolation coefficients reaching in the walls
               do st1=-1,0
                  do st2=-1,0
                     if (this%mask(i+st1,j+st2,k).eq.1) this%itp_xy(st1,st2,i,j,k)=0.0_WP
                     if (this%mask(i,j+st1,k+st2).eq.1) this%itp_yz(st1,st2,i,j,k)=0.0_WP
                     if (this%mask(i+st1,j,k+st2).eq.1) this%itp_xz(st1,st2,i,j,k)=0.0_WP
                  end do
               end do
               ! Rescale to ensure sum(itp)=1
               mysum=sum(this%itp_xy(:,:,i,j,k)); if (mysum.gt.0.0_WP) this%itp_xy(:,:,i,j,k)=this%itp_xy(:,:,i,j,k)/mysum
               mysum=sum(this%itp_yz(:,:,i,j,k)); if (mysum.gt.0.0_WP) this%itp_yz(:,:,i,j,k)=this%itp_yz(:,:,i,j,k)/mysum
               mysum=sum(this%itp_xz(:,:,i,j,k)); if (mysum.gt.0.0_WP) this%itp_xz(:,:,i,j,k)=this%itp_xz(:,:,i,j,k)/mysum
            end do
         end do
      end do
      
      ! Loop over the domain and adjust divergence for P cell
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).gt.0) then
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
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%divp_x=0.0_WP
         this%divu_x=0.0_WP
         this%divv_x=0.0_WP
         this%divw_x=0.0_WP
         this%grdu_x=0.0_WP
         this%grdv_x=0.0_WP
         this%grdw_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%divp_y=0.0_WP
         this%divu_y=0.0_WP
         this%divv_y=0.0_WP
         this%divw_y=0.0_WP
         this%grdu_y=0.0_WP
         this%grdv_y=0.0_WP
         this%grdw_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%divp_z=0.0_WP
         this%divu_z=0.0_WP
         this%divv_z=0.0_WP
         this%divw_z=0.0_WP
         this%grdu_z=0.0_WP
         this%grdv_z=0.0_WP
         this%grdw_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the flow solver now that bconds have been defined
   subroutine setup(this,pressure_solver,implicit_solver)
      implicit none
      class(tpns), intent(inout) :: this
      class(linsol), target, intent(in) :: pressure_solver                      !< A pressure solver is required
      class(linsol), target, intent(in), optional :: implicit_solver            !< An implicit solver can be provided
      integer :: i,j,k

      ! Adjust metrics based on bcflag array
      call this%adjust_metrics()
      
      ! Point to pressure solver linsol object
      this%psolv=>pressure_solver
      
      ! Set 7-pt stencil map for the pressure solver
      this%psolv%stc(1,:)=[ 0, 0, 0]
      this%psolv%stc(2,:)=[+1, 0, 0]
      this%psolv%stc(3,:)=[-1, 0, 0]
      this%psolv%stc(4,:)=[ 0,+1, 0]
      this%psolv%stc(5,:)=[ 0,-1, 0]
      this%psolv%stc(6,:)=[ 0, 0,+1]
      this%psolv%stc(7,:)=[ 0, 0,-1]
      
      ! Setup the scaled Laplacian operator from incomp metrics: lap(*)=-vol*div(grad(*))
      ! Expectations is that this will be replaced later to lap(*)=-vol*div(grad(*)/rho)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Set Laplacian
               this%psolv%opr(1,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x(-1,i+1,j,k)+&
               &                       this%divp_x(0,i,j,k)*this%divu_x( 0,i  ,j,k)+&
               &                       this%divp_y(1,i,j,k)*this%divv_y(-1,i,j+1,k)+&
               &                       this%divp_y(0,i,j,k)*this%divv_y( 0,i,j  ,k)+&
               &                       this%divp_z(1,i,j,k)*this%divw_z(-1,i,j,k+1)+&
               &                       this%divp_z(0,i,j,k)*this%divw_z( 0,i,j,k  )
               this%psolv%opr(2,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x( 0,i+1,j,k)
               this%psolv%opr(3,i,j,k)=this%divp_x(0,i,j,k)*this%divu_x(-1,i  ,j,k)
               this%psolv%opr(4,i,j,k)=this%divp_y(1,i,j,k)*this%divv_y( 0,i,j+1,k)
               this%psolv%opr(5,i,j,k)=this%divp_y(0,i,j,k)*this%divv_y(-1,i,j  ,k)
               this%psolv%opr(6,i,j,k)=this%divp_z(1,i,j,k)*this%divw_z( 0,i,j,k+1)
               this%psolv%opr(7,i,j,k)=this%divp_z(0,i,j,k)*this%divw_z(-1,i,j,k  )
               ! Scale it by the cell volume
               this%psolv%opr(:,i,j,k)=-this%psolv%opr(:,i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do
      
      ! Initialize the pressure Poisson solver
      call this%psolv%init()
      call this%psolv%setup()
      
      ! Prepare implicit solver if it had been provided
      if (present(implicit_solver)) then
         
         ! Point to implicit solver linsol object
         this%implicit=>implicit_solver
         
         ! Set 7-pt stencil map for the velocity solver
         this%implicit%stc(1,:)=[ 0, 0, 0]
         this%implicit%stc(2,:)=[+1, 0, 0]
         this%implicit%stc(3,:)=[-1, 0, 0]
         this%implicit%stc(4,:)=[ 0,+1, 0]
         this%implicit%stc(5,:)=[ 0,-1, 0]
         this%implicit%stc(6,:)=[ 0, 0,+1]
         this%implicit%stc(7,:)=[ 0, 0,-1]
         
         ! Set the diagonal to 1 to make sure all cells participate in solver
         this%implicit%opr(1,:,:,:)=1.0_WP
         
         ! Initialize the implicit velocity solver
         call this%implicit%init()
         
      end if
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,face,dir,canCorrect)
      use string,   only: lowercase
      use messager, only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(tpns), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: type
      procedure(locator_ftype) :: locator
      character(len=1), intent(in) :: face
      integer, intent(in) :: dir
      logical, intent(in) :: canCorrect
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      select case (lowercase(face))
      case ('x'); new_bc%face='x'
      case ('y'); new_bc%face='y'
      case ('z'); new_bc%face='z'
      case default; call die('[tpns add_bcond] Unknown bcond face - expecting x, y, or z')
      end select
      new_bc%itr=iterator(pg=this%cfg,name=new_bc%name,locator=locator,face=new_bc%face)
      select case (dir) ! Outward-oriented
      case (+1); new_bc%dir=+1
      case (-1); new_bc%dir=-1
      case ( 0); new_bc%dir= 0
      case default; call die('[tpns add_bcond] Unknown bcond dir - expecting -1, +1, or 0')
      end select
      new_bc%rdir=real(new_bc%dir,WP)
      new_bc%canCorrect=canCorrect
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet) !< Dirichlet is set one face (i.e., velocit component) at the time
         select case (new_bc%face)
         case ('x')
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%umask(i,j,k)=2
            end do
         case ('y')
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%vmask(i,j,k)=2
            end do
         case ('z')
            do n=1,new_bc%itr%n_
               i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
               this%wmask(i,j,k)=2
            end do
         end select
         
      case (neumann) !< Neumann has to be at existing wall or at domain boundary!
      case (clipped_neumann)
      case (convective)
      case (slip)
      case default
         call die('[tpns apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(tpns), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[tpns get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n,stag
      type(bcond), pointer :: my_bc
      
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
               
            case (dirichlet)               !< Apply Dirichlet conditions
               
               ! This is done by the user directly
               ! Unclear whether we want to do this within the solver...
               
            case (neumann,clipped_neumann,slip) !< Apply Neumann condition to all 3 components
               ! Handle index shift due to staggering
               stag=min(my_bc%dir,0)
               ! Implement based on bcond direction
               select case (my_bc%face)
               case ('x')
                  stag=min(my_bc%dir,0)
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i     ,j    ,k    )=this%U(i-my_bc%dir     ,j    ,k    )
                     this%V(i+stag,j:j+1,k    )=this%V(i-my_bc%dir+stag,j:j+1,k    )
                     this%W(i+stag,j    ,k:k+1)=this%W(i-my_bc%dir+stag,j    ,k:k+1)
                     this%rhoU(i     ,j    ,k    )=this%rhoU(i-my_bc%dir     ,j    ,k    )
                     this%rhoV(i+stag,j:j+1,k    )=this%rhoV(i-my_bc%dir+stag,j:j+1,k    )
                     this%rhoW(i+stag,j    ,k:k+1)=this%rhoW(i-my_bc%dir+stag,j    ,k:k+1)
                  end do
               case ('y')
                  stag=min(my_bc%dir,0)
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+stag,k    )=this%U(i:i+1,j-my_bc%dir+stag,k    )
                     this%V(i    ,j     ,k    )=this%V(i    ,j-my_bc%dir     ,k    )
                     this%W(i    ,j+stag,k:k+1)=this%W(i    ,j-my_bc%dir+stag,k:k+1)
                     this%rhoU(i:i+1,j+stag,k    )=this%rhoU(i:i+1,j-my_bc%dir+stag,k    )
                     this%rhoV(i    ,j     ,k    )=this%rhoV(i    ,j-my_bc%dir     ,k    )
                     this%rhoW(i    ,j+stag,k:k+1)=this%rhoW(i    ,j-my_bc%dir+stag,k:k+1)
                  end do
               case ('z')
                  stag=min(my_bc%dir,0)
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j    ,k+stag)=this%U(i:i+1,j    ,k-my_bc%dir+stag)
                     this%V(i    ,j:j+1,k+stag)=this%V(i    ,j:j+1,k-my_bc%dir+stag)
                     this%W(i    ,j    ,k     )=this%W(i    ,j    ,k-my_bc%dir     )
                     this%rhoU(i:i+1,j    ,k+stag)=this%rhoU(i:i+1,j    ,k-my_bc%dir+stag)
                     this%rhoV(i    ,j:j+1,k+stag)=this%rhoV(i    ,j:j+1,k-my_bc%dir+stag)
                     this%rhoW(i    ,j    ,k     )=this%rhoW(i    ,j    ,k-my_bc%dir     )
                  end do
               end select
               ! If needed, clip
               if (my_bc%type.eq.clipped_neumann) then
                  select case (my_bc%face)
                  case ('x')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        if (this%U(i,j,k)*my_bc%rdir.lt.0.0_WP) then
                           this%U(i,j,k)=0.0_WP
                           this%rhoU(i,j,k)=0.0_WP
                        end if
                     end do
                  case ('y')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        if (this%V(i,j,k)*my_bc%rdir.lt.0.0_WP) then
                           this%V(i,j,k)=0.0_WP
                           this%rhoV(i,j,k)=0.0_WP
                        end if
                     end do
                  case ('z')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        if (this%W(i,j,k)*my_bc%rdir.lt.0.0_WP) then
                           this%W(i,j,k)=0.0_WP
                           this%rhoW(i,j,k)=0.0_WP
                        end if
                     end do
                  end select
               end if
               ! If needed, no penetration
               if (my_bc%type.eq.slip) then
                  select case (my_bc%face)
                  case ('x')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        this%U(i,j,k)=0.0_WP
                        this%rhoU(i,j,k)=0.0_WP
                     end do
                  case ('y')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        this%V(i,j,k)=0.0_WP
                        this%rhoV(i,j,k)=0.0_WP
                     end do
                  case ('z')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        this%W(i,j,k)=0.0_WP
                        this%rhoW(i,j,k)=0.0_WP
                     end do
                  end select
               end if
               
            case (convective)   ! Not implemented yet!
               
            case default
               call die('[tpns apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! Sync full fields after each bcond - this should be optimized
      call this%cfg%sync(this%U)
      call this%cfg%sync(this%V)
      call this%cfg%sync(this%W)
      call this%cfg%sync(this%rhoU)
      call this%cfg%sync(this%rhoV)
      call this%cfg%sync(this%rhoW)
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit momentum time derivative
   !> This assumes that rho, rhoold, rhoU/V/W have been updated already
   subroutine get_dmomdt(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      
      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
      ! Initialize the hybrid interpolator to centered scheme
      this%hybu_x=this%itpu_x; this%hybv_x=this%itpv_x; this%hybw_x=this%itpw_x
      this%hybu_y=this%itpu_y; this%hybv_y=this%itpv_y; this%hybw_y=this%itpw_y
      this%hybu_z=this%itpu_z; this%hybv_z=this%itpv_z; this%hybw_z=this%itpw_z
      
      ! Prepare hybrid interpolation scheme for mass and momentum
      ! do kk=this%cfg%kmin_,this%cfg%kmax_+1
      !    do jj=this%cfg%jmin_,this%cfg%jmax_+1
      !       do ii=this%cfg%imin_,this%cfg%imax_+1
      !          ! U-cell **********
      !          ! Fluxes on x-face
      !          i=ii-1; j=jj-1; k=kk-1
      !          if (abs(this%rho_Uold(i+1,j,k)-this%rho_Uold(i,j,k)).gt.rhoeps) then
      !             vel=sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybu_x(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybu_x(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! Fluxes on y-face
      !          i=ii; j=jj; k=kk
      !          if (abs(this%rho_Uold(i,j,k)-this%rho_Uold(i,j-1,k)).gt.rhoeps) then
      !             vel=sum(this%itpv_x(:,i,j,k)*this%V(i-1:i,j,k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybu_y(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybu_y(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! Fluxes on z-face
      !          i=ii; j=jj; k=kk
      !          if (abs(this%rho_Uold(i,j,k)-this%rho_Uold(i,j,k-1)).gt.rhoeps) then
      !             vel=sum(this%itpw_x(:,i,j,k)*this%W(i-1:i,j,k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybu_z(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybu_z(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! V-cell **********
      !          ! Fluxes on x-face
      !          i=ii; j=jj; k=kk
      !          if (abs(this%rho_Vold(i,j,k)-this%rho_Vold(i-1,j,k)).gt.rhoeps) then
      !             vel=sum(this%itpu_y(:,i,j,k)*this%U(i,j-1:j,k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybv_x(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybv_x(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! Fluxes on y-face
      !          i=ii-1; j=jj-1; k=kk-1
      !          if (abs(this%rho_Vold(i,j+1,k)-this%rho_Vold(i,j,k)).gt.rhoeps) then
      !             vel=sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybv_y(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybv_y(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! Fluxes on z-face
      !          i=ii; j=jj; k=kk
      !          if (abs(this%rho_Vold(i,j,k)-this%rho_Vold(i,j,k-1)).gt.rhoeps) then
      !             vel=sum(this%itpw_y(:,i,j,k)*this%W(i,j-1:j,k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybv_z(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybv_z(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! W-cell **********
      !          ! Fluxes on x-face
      !          i=ii; j=jj; k=kk
      !          if (abs(this%rho_Wold(i,j,k)-this%rho_Wold(i-1,j,k)).gt.rhoeps) then
      !             vel=sum(this%itpu_z(:,i,j,k)*this%U(i,j,k-1:k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybw_x(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybw_x(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! Fluxes on y-face
      !          i=ii; j=jj; k=kk
      !          if (abs(this%rho_Wold(i,j,k)-this%rho_Wold(i,j-1,k)).gt.rhoeps) then
      !             vel=sum(this%itpv_z(:,i,j,k)*this%V(i,j,k-1:k))
      !             if (vel.ge.0.0_WP) then
      !                this%hybw_y(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybw_y(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !          ! Fluxes on z-face
      !          i=ii-1; j=jj-1; k=kk-1
      !          if (abs(this%rho_Wold(i,j,k+1)-this%rho_Wold(i,j,k)).gt.rhoeps) then
      !             vel=sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1))
      !             if (vel.ge.0.0_WP) then
      !                this%hybw_z(:,i,j,k)=[1.0_WP,0.0_WP]
      !             else
      !                this%hybw_z(:,i,j,k)=[0.0_WP,1.0_WP]
      !             end if
      !          end if
      !       end do
      !    end do
      ! end do
      
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Flux of rhoU
      do kk=this%cfg%kmin_,this%cfg%kmax_+1
         do jj=this%cfg%jmin_,this%cfg%jmax_+1
            do ii=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               i=ii-1; j=jj-1; k=kk-1
               FX(i,j,k)=-sum(this%itpu_x(:,i,j,k)*this%rhoU(i:i+1,j,k))*sum(this%hybu_x(:,i,j,k)*this%U(i:i+1,j,k)) &
               &         +this%visc   (i,j,k)*(sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))+sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))) &
               &         -this%P(i,j,k)
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               FY(i,j,k)=-sum(this%itpv_x(:,i,j,k)*this%rhoV(i-1:i,j,k))*sum(this%hybu_y(:,i,j,k)*this%U(i,j-1:j,k)) &
               &         +this%visc_xy(i,j,k)*(sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))+sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               FZ(i,j,k)=-sum(this%itpw_x(:,i,j,k)*this%rhoW(i-1:i,j,k))*sum(this%hybu_z(:,i,j,k)*this%U(i,j,k-1:k)) &
               &         +this%visc_zx(i,j,k)*(sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))+sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoU
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoUdt(i,j,k)=sum(this%divu_x(:,i,j,k)*FX(i-1:i,j,k))+&
               &              sum(this%divu_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divu_z(:,i,j,k)*FZ(i,j,k:k+1))+this%Pjx(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoUdt)
      
      ! Flux of rhoV
      do kk=this%cfg%kmin_,this%cfg%kmax_+1
         do jj=this%cfg%jmin_,this%cfg%jmax_+1
            do ii=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               FX(i,j,k)=-sum(this%itpu_y(:,i,j,k)*this%rhoU(i,j-1:j,k))*sum(this%hybv_x(:,i,j,k)*this%V(i-1:i,j,k)) &
               &         +this%visc_xy(i,j,k)*(sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))+sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii-1; j=jj-1; k=kk-1
               FY(i,j,k)=-sum(this%itpv_y(:,i,j,k)*this%rhoV(i,j:j+1,k))*sum(this%hybv_y(:,i,j,k)*this%V(i,j:j+1,k)) &
               &         +this%visc   (i,j,k)*(sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))+sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))) &
               &         -this%P(i,j,k)
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               FZ(i,j,k)=-sum(this%itpw_y(:,i,j,k)*this%rhoW(i,j-1:j,k))*sum(this%hybv_z(:,i,j,k)*this%V(i,j,k-1:k)) &
               &         +this%visc_yz(i,j,k)*(sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))+sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k)))
            end do
         end do
      end do
      ! Time derivative of rhoV
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoVdt(i,j,k)=sum(this%divv_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divv_y(:,i,j,k)*FY(i,j-1:j,k))+&
               &              sum(this%divv_z(:,i,j,k)*FZ(i,j,k:k+1))+this%Pjy(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoVdt)
      
      ! Flux of rhoW
      do kk=this%cfg%kmin_,this%cfg%kmax_+1
         do jj=this%cfg%jmin_,this%cfg%jmax_+1
            do ii=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
               i=ii; j=jj; k=kk
               FX(i,j,k)=-sum(this%itpu_z(:,i,j,k)*this%rhoW(i,j,k-1:k))*sum(this%hybw_x(:,i,j,k)*this%W(i-1:i,j,k)) &
               &         +this%visc_zx(i,j,k)*(sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))+sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               FY(i,j,k)=-sum(this%itpv_z(:,i,j,k)*this%rhoW(i,j,k-1:k))*sum(this%hybw_y(:,i,j,k)*this%W(i,j-1:j,k)) &
               &         +this%visc_yz(i,j,k)*(sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))+sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii-1; j=jj-1; k=kk-1
               FZ(i,j,k)=-sum(this%itpw_z(:,i,j,k)*this%rhoW(i,j,k:k+1))*sum(this%hybw_z(:,i,j,k)*this%W(i,j,k:k+1)) &
               &         +this%visc   (i,j,k)*(sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))+sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))) &
               &         -this%P(i,j,k)
            end do
         end do
      end do
      ! Time derivative of rhoW
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoWdt(i,j,k)=sum(this%divw_x(:,i,j,k)*FX(i:i+1,j,k))+&
               &              sum(this%divw_y(:,i,j,k)*FY(i,j:j+1,k))+&
               &              sum(this%divw_z(:,i,j,k)*FZ(i,j,k-1:k))+this%Pjz(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(drhoWdt)
      
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      
   end subroutine get_dmomdt
   
   
   !> Update pressure Poisson operator
   subroutine update_laplacian(this)
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k,s1,s2
      ! Setup the scaled Laplacian operator from  metrics: lap(*)=-vol.div(grad(*)/rho)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Zero out Laplacian
               this%psolv%opr(:,i,j,k)=0.0_WP
               ! Tranverse the stencil and recompute Laplacian
               do s1=0,1
                  do s2=-1,0
                     this%psolv%opr(this%psolv%stmap(s1+s2,0,0),i,j,k)=this%psolv%opr(this%psolv%stmap(s1+s2,0,0),i,j,k)+this%divp_x(s1,i,j,k)*this%divu_x(s2,i+s1,j,k)/this%rho_U(i+s1,j,k)
                     this%psolv%opr(this%psolv%stmap(0,s1+s2,0),i,j,k)=this%psolv%opr(this%psolv%stmap(0,s1+s2,0),i,j,k)+this%divp_y(s1,i,j,k)*this%divv_y(s2,i,j+s1,k)/this%rho_V(i,j+s1,k)
                     this%psolv%opr(this%psolv%stmap(0,0,s1+s2),i,j,k)=this%psolv%opr(this%psolv%stmap(0,0,s1+s2),i,j,k)+this%divp_z(s1,i,j,k)*this%divw_z(s2,i,j,k+s1)/this%rho_W(i,j,k+s1)
                  end do
               end do
               ! Scale Laplacian by cell volume
               this%psolv%opr(:,i,j,k)=-this%psolv%opr(:,i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do
      ! Initialize the pressure Poisson solver
      call this%psolv%setup()
   end subroutine update_laplacian
   
   
   !> Calculate the velocity divergence based on U/V/W
   subroutine get_div(this)
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div(i,j,k)=sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k))+&
               &               sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k))+&
               &               sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(this%div)
   end subroutine get_div
   
   
   !> Add surface tension jump term using CSF
   subroutine add_surface_tension_jump(this,dt,div,vf,contact_model)
      use messager,  only: die
      use vfs_class, only: vfs
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), intent(inout) :: dt     !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      class(vfs), intent(inout) :: vf
      integer, intent(in), optional :: contact_model
      integer :: i,j,k,s1
      real(WP) :: mycurv,mysurf
      
      ! Store old jump
      this%DPjx=this%Pjx
      this%DPjy=this%Pjy
      this%DPjz=this%Pjz
      
      ! Calculate pressure jump
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face
               mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
               if (mysurf.gt.0.0_WP) then
                  mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
               else
                  mycurv=0.0_WP
               end if
               this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k))
               ! Y face
               mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
               if (mysurf.gt.0.0_WP) then
                  mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
               else
                  mycurv=0.0_WP
               end if
               this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*vf%VF(i,j-1:j,k))
               ! Z face
               mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
               if (mysurf.gt.0.0_WP) then
                  mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
               else
                  mycurv=0.0_WP
               end if
               this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k))
            end do
         end do
      end do
      
      ! Add wall contact force to pressure jump
      if (present(contact_model)) then
         select case (contact_model)
         case (static_contact)
            call this%add_static_contact(vf=vf)
         case default
            call die('[tpns: add_surface_tension_jump] Unknown contact model!')
         end select
      end if
      
      ! Compute jump of DP
      this%DPjx=this%Pjx-this%DPjx
      this%DPjy=this%Pjy-this%DPjy
      this%DPjz=this%Pjz-this%DPjz
      
      ! Add div(Pjump) to RP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do s1=0,1
                  div(i,j,k)=div(i,j,k)+dt*this%divp_x(s1,i,j,k)*this%DPjx(i+s1,j,k)/sum(this%itpr_x(:,i+s1,j,k)*this%rho(i+s1-1:i+s1,j,k))
                  div(i,j,k)=div(i,j,k)+dt*this%divp_y(s1,i,j,k)*this%DPjy(i,j+s1,k)/sum(this%itpr_y(:,i,j+s1,k)*this%rho(i,j+s1-1:j+s1,k))
                  div(i,j,k)=div(i,j,k)+dt*this%divp_z(s1,i,j,k)*this%DPjz(i,j,k+s1)/sum(this%itpr_z(:,i,j,k+s1)*this%rho(i,j,k+s1-1:k+s1))
               end do
            end do
         end do
      end do
      
   end subroutine add_surface_tension_jump
   
   
   !> Calculate the pressure gradient based on P
   subroutine get_pgrad(this,P,Pgradx,Pgrady,Pgradz)
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: P      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradx !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgrady !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradz !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               Pgradx(i,j,k)=sum(this%divu_x(:,i,j,k)*P(i-1:i,j,k))-this%dPjx(i,j,k)
               Pgrady(i,j,k)=sum(this%divv_y(:,i,j,k)*P(i,j-1:j,k))-this%dPjy(i,j,k)
               Pgradz(i,j,k)=sum(this%divw_z(:,i,j,k)*P(i,j,k-1:k))-this%dPjz(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(Pgradx)
      call this%cfg%sync(Pgrady)
      call this%cfg%sync(Pgradz)
   end subroutine get_pgrad
   
   
   !> Calculate the interpolated velocity, including overlap and ghosts
   subroutine interp_vel(this,Ui,Vi,Wi)
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate as far as possible each component
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               Ui(i,j,k)=sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k))
            end do
         end do
      end do
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_
               Vi(i,j,k)=sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k))
            end do
         end do
      end do
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               Wi(i,j,k)=sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1))
            end do
         end do
      end do
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) Ui(this%cfg%imaxo,:,:)=this%U(this%cfg%imaxo,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) Vi(:,this%cfg%jmaxo,:)=this%V(:,this%cfg%jmaxo,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) Wi(:,:,this%cfg%kmaxo)=this%W(:,:,this%cfg%kmaxo)
      ! Sync it
      call this%cfg%sync(Ui)
      call this%cfg%sync(Vi)
      call this%cfg%sync(Wi)
   end subroutine interp_vel
   
   
   !> Calculate the deviatoric part of the strain rate tensor from U/V/W
   !> 1: du/dx-div/3
   !> 2: dv/dy-div/3
   !> 3: dw/dz-div/3
   !> 4: (du/dy+dv/dx)/2
   !> 5: (dv/dz+dw/dy)/2
   !> 6: (dw/dx+du/dz)/2
   subroutine get_strainrate(this,SR)
      use messager, only: die
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: SR  !< Needs to be (1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      real(WP) :: div
      integer :: i,j,k
      
      ! Check SR's first dimension
	   if (size(SR,dim=1).ne.6) call die('[tpns get_strainrate] SR should be of size (1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
      ! Compute dudx, dvdy, and dwdz first
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               SR(1,i,j,k)=sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))
               SR(2,i,j,k)=sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))
               SR(3,i,j,k)=sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))
               div=sum(SR(1:3,i,j,k))/3.0_WP
               SR(1,i,j,k)=SR(1,i,j,k)-div
               SR(2,i,j,k)=SR(2,i,j,k)-div
               SR(3,i,j,k)=SR(3,i,j,k)-div
            end do
         end do
      end do
      
      ! Allocate velocity gradient components
	   allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
	   do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               dudy(i,j,k)=sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))
               dudz(i,j,k)=sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))
               dvdx(i,j,k)=sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))
               dvdz(i,j,k)=sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))
               dwdx(i,j,k)=sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))
               dwdy(i,j,k)=sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Interpolate off-diagonal components of the velocity gradient to the cell center and store strain rate
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               SR(4,i,j,k)=0.125_WP*(sum(dudy(i:i+1,j:j+1,k    ))+sum(dvdx(i:i+1,j:j+1,k    )))
               SR(5,i,j,k)=0.125_WP*(sum(dvdz(i    ,j:j+1,k:k+1))+sum(dwdy(i    ,j:j+1,k:k+1)))
               SR(6,i,j,k)=0.125_WP*(sum(dwdx(i:i+1,j    ,k:k+1))+sum(dudz(i:i+1,j    ,k:k+1)))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            SR(:,this%cfg%imin-1,:,:)=SR(:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) SR(:,this%cfg%imax+1,:,:)=SR(:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            SR(:,:,this%cfg%jmin-1,:)=SR(:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) SR(:,:,this%cfg%jmax+1,:)=SR(:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            SR(:,:,:,this%cfg%kmin-1)=SR(:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) SR(:,:,:,this%cfg%kmax+1)=SR(:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
	   do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) SR(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Sync it
	   call this%cfg%sync(SR)
      
      ! Deallocate velocity gradient storage
	   deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
      
   end subroutine get_strainrate

   
   !> Calculate the velocity gradient tensor from U/V/W
   !> Note that gradu(i,j)=duj/dxi
   subroutine get_gradu(this,gradu)
      use messager, only: die
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: gradu  !< Needs to be (1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      
      ! Check gradu's first two dimensions
	   if (size(gradu,dim=1).ne.3.or.size(gradu,dim=2).ne.3) call die('[tpns get_strainrate] gradu should be of size (1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
      ! Compute dudx, dvdy, and dwdz first
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               gradu(1,1,i,j,k)=sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))
               gradu(2,2,i,j,k)=sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))
               gradu(3,3,i,j,k)=sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))
            end do
         end do
      end do
      
      ! Allocate velocity gradient components
	   allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
	   do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               dudy(i,j,k)=sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))
               dudz(i,j,k)=sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))
               dvdx(i,j,k)=sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))
               dvdz(i,j,k)=sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))
               dwdx(i,j,k)=sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))
               dwdy(i,j,k)=sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Interpolate off-diagonal components of the velocity gradient to the cell center
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               gradu(2,1,i,j,k)=0.25_WP*sum(dudy(i:i+1,j:j+1,k))
               gradu(3,1,i,j,k)=0.25_WP*sum(dudz(i:i+1,j,k:k+1))
               gradu(1,2,i,j,k)=0.25_WP*sum(dvdx(i:i+1,j:j+1,k))
               gradu(3,2,i,j,k)=0.25_WP*sum(dvdz(i,j:j+1,k:k+1))
               gradu(1,3,i,j,k)=0.25_WP*sum(dwdx(i:i+1,j,k:k+1))
               gradu(2,3,i,j,k)=0.25_WP*sum(dwdy(i,j:j+1,k:k+1))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            gradu(:,:,this%cfg%imin-1,:,:)=gradu(:,:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) gradu(:,:,this%cfg%imax+1,:,:)=gradu(:,:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            gradu(:,:,:,this%cfg%jmin-1,:)=gradu(:,:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) gradu(:,:,:,this%cfg%jmax+1,:)=gradu(:,:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            gradu(:,:,:,:,this%cfg%kmin-1)=gradu(:,:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) gradu(:,:,:,:,this%cfg%kmax+1)=gradu(:,:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
	   do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) gradu(:,:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Sync it
	   call this%cfg%sync(gradu)
      
      ! Deallocate velocity gradient storage
	   deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
      
   end subroutine get_gradu
   
   
   !> Calculate vorticity vector
   subroutine get_vorticity(this,vort)
      use messager, only: die
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: vort  !< Needs to be (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      
      ! Check vort's first two dimensions
      if (size(vort,dim=1).ne.3) call die('[tpns get_vorticity] vort should be of size (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')

      ! Allocate velocity gradient components
      allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               dudy(i,j,k)=sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))
               dudz(i,j,k)=sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k))
               dvdx(i,j,k)=sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))
               dvdz(i,j,k)=sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k))
               dwdx(i,j,k)=sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))
               dwdy(i,j,k)=sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))
            end do
         end do
      end do

      ! Interpolate off-diagonal components of the velocity gradient to the cell center
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               vort(1,i,j,k)=0.25_WP*(sum(dwdy(i,j:j+1,k:k+1))-sum(dvdz(i,j:j+1,k:k+1)))
               vort(2,i,j,k)=0.25_WP*(sum(dudz(i:i+1,j,k:k+1))-sum(dwdx(i:i+1,j,k:k+1)))
               vort(3,i,j,k)=0.25_WP*(sum(dvdx(i:i+1,j:j+1,k))-sum(dudy(i:i+1,j:j+1,k)))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            vort(:,this%cfg%imin-1,:,:)=vort(:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) vort(:,this%cfg%imax+1,:,:)=vort(:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            vort(:,:,this%cfg%jmin-1,:)=vort(:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) vort(:,:,this%cfg%jmax+1,:)=vort(:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            vort(:,:,:,this%cfg%kmin-1)=vort(:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) vort(:,:,:,this%cfg%kmax+1)=vort(:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) vort(:,i,j,k)=0.0_WP
            end do
         end do
      end do

      ! Sync it
      call this%cfg%sync(vort)

      ! Deallocate velocity gradient storage
      deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)

   end subroutine get_vorticity

   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cflc,cfl)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFLc_x,my_CFLc_y,my_CFLc_z,my_CFLv_x,my_CFLv_y,my_CFLv_z,my_CFLst
      real(WP) :: max_nu
      
      ! Get surface tension CFL first
      my_CFLst=huge(1.0_WP)
      if (this%sigma.gt.0.0_WP) then
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  my_CFLst=min(my_CFLst,sqrt((this%rho_l+this%rho_g)*this%cfg%meshsize(i,j,k)**3.0_WP/(4.0_WP*Pi*this%sigma)))
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(my_CFLst,this%CFLst,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      this%CFLst=dt/this%CFLst
      
      ! Get largest kinematic viscosity
      max_nu=max(this%visc_l/this%rho_l,this%visc_g/this%rho_g)
      
      ! Set the CFLs to zero
      my_CFLc_x=0.0_WP; my_CFLc_y=0.0_WP; my_CFLc_z=0.0_WP
      my_CFLv_x=0.0_WP; my_CFLv_y=0.0_WP; my_CFLv_z=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFLc_x=max(my_CFLc_x,abs(this%U(i,j,k))*this%cfg%dxmi(i))
               my_CFLc_y=max(my_CFLc_y,abs(this%V(i,j,k))*this%cfg%dymi(j))
               my_CFLc_z=max(my_CFLc_z,abs(this%W(i,j,k))*this%cfg%dzmi(k))
               my_CFLv_x=max(my_CFLv_x,4.0_WP*max_nu*this%cfg%dxi(i)**2)
               my_CFLv_y=max(my_CFLv_y,4.0_WP*max_nu*this%cfg%dyi(j)**2)
               my_CFLv_z=max(my_CFLv_z,4.0_WP*max_nu*this%cfg%dzi(k)**2)
            end do
         end do
      end do
      my_CFLc_x=my_CFLc_x*dt; my_CFLc_y=my_CFLc_y*dt; my_CFLc_z=my_CFLc_z*dt
      my_CFLv_x=my_CFLv_x*dt; my_CFLv_y=my_CFLv_y*dt; my_CFLv_z=my_CFLv_z*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFLc_x,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_y,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_z,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_x,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_y,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_z,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Return the maximum convective + surface tension CFL
      cflc=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLst)
      
      ! If asked for, also return the maximum overall CFL
      if (present(CFL)) cfl =max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z,this%CFLst)
      
   end subroutine get_cfl
   
   
   !> Calculate the max of our fields
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: my_Umax,my_Vmax,my_Wmax,my_Pmax,my_divmax
      
      ! Set all to zero
      my_Umax=0.0_WP; my_Vmax=0.0_WP; my_Wmax=0.0_WP; my_Pmax=0.0_WP; my_divmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_Umax=max(my_Umax,abs(this%U(i,j,k)))
               my_Vmax=max(my_Vmax,abs(this%V(i,j,k)))
               my_Wmax=max(my_Wmax,abs(this%W(i,j,k)))
               if (this%cfg%VF(i,j,k).gt.0.0_WP) my_Pmax  =max(my_Pmax  ,abs(this%P(i,j,k)  ))
               if (this%cfg%VF(i,j,k).gt.0.0_WP) my_divmax=max(my_divmax,abs(this%div(i,j,k)))
            end do
         end do
      end do
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Umax  ,this%Umax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Vmax  ,this%Vmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Wmax  ,this%Wmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Pmax  ,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_divmax,this%divmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_max
   
   
   !> Compute MFR through all bcs
   subroutine get_mfr(this)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k,n,ibc,ierr
      type(bcond), pointer :: my_bc
      real(WP), dimension(:), allocatable :: my_mfr,my_area
      real(WP), dimension(:), allocatable :: canCorrect
      
      ! Ensure this%mfr is of proper size
      if (.not.allocated(this%mfr)) then
         allocate(this%mfr(this%nbc))
      else
         if (size(this%mfr).ne.this%nbc) then
            deallocate(this%mfr); allocate(this%mfr(this%nbc))
         end if
      end if
      
      ! Ensure this%area is of proper size
      if (.not.allocated(this%area)) then
         allocate(this%area(this%nbc))
      else
         if (size(this%area).ne.this%nbc) then
            deallocate(this%area); allocate(this%area(this%nbc))
         end if
      end if
      
      ! Allocate temp array for communication
      allocate(my_mfr(this%nbc))
      allocate(my_area(this%nbc))
      allocate(canCorrect(this%nbc))
      
      ! Traverse bcond list and integrate local outgoing MFR
      my_bc=>this%first_bc; ibc=1
      do while (associated(my_bc))
         
         ! Set zero local MFR and area
         my_mfr(ibc)=0.0_WP
         my_area(ibc)=0.0_WP
         if (my_bc%canCorrect) then
            canCorrect(ibc)=1.0_WP
         else
            canCorrect(ibc)=0.0_WP
         end if
         
         ! Only processes inside the bcond have a non-zero MFR
         if (my_bc%itr%amIn) then
            
            ! Implement based on bcond face and dir, loop over interior only
            select case (my_bc%face)
            case ('x')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  my_mfr(ibc)=my_mfr(ibc)+my_bc%rdir*this%U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                  my_area(ibc)=my_area(ibc)+this%cfg%dy(j)*this%cfg%dz(k)
               end do
            case ('y')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  my_mfr(ibc)=my_mfr(ibc)+my_bc%rdir*this%V(i,j,k)*this%cfg%dz(k)*this%cfg%dx(i)
                  my_area(ibc)=my_area(ibc)+this%cfg%dz(k)*this%cfg%dx(i)
               end do
            case ('z')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  my_mfr(ibc)=my_mfr(ibc)+my_bc%rdir*this%W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j)
                  my_area(ibc)=my_area(ibc)+this%cfg%dx(i)*this%cfg%dy(j)
               end do
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next; ibc=ibc+1
         
      end do
      
      ! Sum up all values
      call MPI_ALLREDUCE(my_mfr ,this%mfr ,this%nbc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_area,this%area,this%nbc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
      ! Compute the correctable area
      this%correctable_area=sum(this%area*canCorrect)
      
      ! Deallocate temp array
      deallocate(my_mfr,my_area,canCorrect)
      
   end subroutine get_mfr
   
   
   !> Correct MFR through correctable bconds
   subroutine correct_mfr(this)
      use mpi_f08, only: MPI_SUM
      implicit none
      class(tpns), intent(inout) :: this
      real(WP) :: mfr_error,vel_correction
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Evaluate MFR mismatch and velocity correction
      call this%get_mfr()
      mfr_error=sum(this%mfr)
      if (abs(mfr_error).lt.10.0_WP*epsilon(1.0_WP).or.abs(this%correctable_area).lt.10.0_WP*epsilon(1.0_WP)) return
      vel_correction=-mfr_error/(this%correctable_area)
      
      ! Traverse bcond list and correct bcond MFR
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside correctable bcond need to work
         if (my_bc%itr%amIn.and.my_bc%canCorrect) then
            
            ! Implement based on bcond direction, loop over all cell
            select case (my_bc%face)
            case ('x')
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%U(i,j,k)=this%U(i,j,k)+my_bc%rdir*vel_correction
               end do
            case ('y')
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%V(i,j,k)=this%V(i,j,k)+my_bc%rdir*vel_correction
               end do
            case ('z')
               do n=1,my_bc%itr%no_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%W(i,j,k)=this%W(i,j,k)+my_bc%rdir*vel_correction
               end do
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine correct_mfr
   
   
   !> Shift pressure to ensure zero average
   subroutine shift_p(this,pressure)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpns), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: pressure !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: vol_tot,pressure_tot,my_vol_tot,my_pressure_tot
      integer :: i,j,k,ierr
      
      ! Loop over domain and integrate volume and pressure
      my_vol_tot=0.0_WP
      my_pressure_tot=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_vol_tot     =my_vol_tot     +this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)
               my_pressure_tot=my_pressure_tot+this%cfg%vol(i,j,k)*this%cfg%VF(i,j,k)*pressure(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_vol_tot     ,vol_tot     ,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_pressure_tot,pressure_tot,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      pressure_tot=pressure_tot/vol_tot
      
      ! Shift the pressure
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%VF(i,j,k).gt.0.0_WP) pressure(i,j,k)=pressure(i,j,k)-pressure_tot
            end do
         end do
      end do
      call this%cfg%sync(pressure)
      
   end subroutine shift_p
   
   
   !> Solve for implicit velocity residual
   subroutine solve_implicit(this,dt,resU,resV,resW)
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP) :: rhoUp,rhoUm,rhoVp,rhoVm,rhoWp,rhoWm
      
      ! Solve implicit U problem
      this%implicit%opr(1,:,:,:)=sum(this%itpr_x(:,i,j,k)*this%rho(i-1:i,j,k)); this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=sum(this%itpu_x(:,i  ,j,k)*this%rhoU(i  :i+1,j,k))
               rhoUm=sum(this%itpu_x(:,i-1,j,k)*this%rhoU(i-1:i  ,j,k))
               rhoVp=sum(this%itpv_x(:,i,j+1,k)*this%rhoV(i-1:i,j+1,k))
               rhoVm=sum(this%itpv_x(:,i,j  ,k)*this%rhoV(i-1:i,j  ,k))
               rhoWp=sum(this%itpw_x(:,i,j,k+1)*this%rhoW(i-1:i,j,k+1))
               rhoWm=sum(this%itpw_x(:,i,j,k  )*this%rhoW(i-1:i,j,k  ))
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
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divu_x( 0,i,j,k)*2.0_WP*this%visc   (i  ,j,k)*this%grdu_x( 0,i  ,j,k)+&
               &                                                                this%divu_x(-1,i,j,k)*2.0_WP*this%visc   (i-1,j,k)*this%grdu_x(+1,i-1,j,k)+&
               &                                                                this%divu_y(+1,i,j,k)*       this%visc_xy(i,j+1,k)*this%grdu_y(-1,i,j+1,k)+&
               &                                                                this%divu_y( 0,i,j,k)*       this%visc_xy(i,j  ,k)*this%grdu_y( 0,i,j  ,k)+&
               &                                                                this%divu_z(+1,i,j,k)*       this%visc_zx(i,j,k+1)*this%grdu_z(-1,i,j,k+1)+&
               &                                                                this%divu_z( 0,i,j,k)*       this%visc_zx(i,j,k  )*this%grdu_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divu_x( 0,i,j,k)*2.0_WP*this%visc   (i  ,j,k)*this%grdu_x(+1,i  ,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divu_x(-1,i,j,k)*2.0_WP*this%visc   (i-1,j,k)*this%grdu_x( 0,i-1,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divu_y(+1,i,j,k)*       this%visc_xy(i,j+1,k)*this%grdu_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divu_y( 0,i,j,k)*       this%visc_xy(i,j  ,k)*this%grdu_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divu_z(+1,i,j,k)*       this%visc_zx(i,j,k+1)*this%grdu_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divu_z( 0,i,j,k)*       this%visc_zx(i,j,k  )*this%grdu_z(-1,i,j,k  ))
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resU
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resU=this%implicit%sol
      
      ! Solve implicit V problem
      this%implicit%opr(1,:,:,:)=sum(this%itpr_y(:,i,j,k)*this%rho(i,j-1:j,k)); this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=sum(this%itpu_y(:,i+1,j,k)*this%rhoU(i+1,j-1:j,k))
               rhoUm=sum(this%itpu_y(:,i  ,j,k)*this%rhoU(i  ,j-1:j,k))
               rhoVp=sum(this%itpv_y(:,i,j  ,k)*this%rhoV(i,j  :j+1,k))
               rhoVm=sum(this%itpv_y(:,i,j-1,k)*this%rhoV(i,j-1:j  ,k))
               rhoWp=sum(this%itpw_y(:,i,j,k+1)*this%rhoW(i,j-1:j,k+1))
               rhoWm=sum(this%itpw_y(:,i,j,k  )*this%rhoW(i,j-1:j,k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+0.5_WP*dt*(this%divv_x(+1,i,j,k)*this%itpv_x(-1,i+1,j,k)*rhoUp+&
               &                                                                this%divv_x( 0,i,j,k)*this%itpv_x( 0,i  ,j,k)*rhoUm+&
               &                                                                this%divv_y( 0,i,j,k)*this%itpv_y( 0,i,j  ,k)*rhoVp+&
               &                                                                this%divv_y(-1,i,j,k)*this%itpv_y(+1,i,j-1,k)*rhoVm+&
               &                                                                this%divv_z(+1,i,j,k)*this%itpv_z(-1,i,j,k+1)*rhoWp+&
               &                                                                this%divv_z( 0,i,j,k)*this%itpv_z( 0,i,j,k  )*rhoWm)
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+0.5_WP*dt*(this%divv_x(+1,i,j,k)*this%itpv_x( 0,i+1,j,k)*rhoUp)
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+0.5_WP*dt*(this%divv_x( 0,i,j,k)*this%itpv_x(-1,i  ,j,k)*rhoUm)
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+0.5_WP*dt*(this%divv_y( 0,i,j,k)*this%itpv_y(+1,i,j  ,k)*rhoVp)
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+0.5_WP*dt*(this%divv_y(-1,i,j,k)*this%itpv_y( 0,i,j-1,k)*rhoVm)
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+0.5_WP*dt*(this%divv_z(+1,i,j,k)*this%itpv_z( 0,i,j,k+1)*rhoWp)
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+0.5_WP*dt*(this%divv_z( 0,i,j,k)*this%itpv_z(-1,i,j,k  )*rhoWm)
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divv_x(+1,i,j,k)*       this%visc_xy(i+1,j,k)*this%grdv_x(-1,i+1,j,k)+&
               &                                                                this%divv_x( 0,i,j,k)*       this%visc_xy(i  ,j,k)*this%grdv_x( 0,i  ,j,k)+&
               &                                                                this%divv_y( 0,i,j,k)*2.0_WP*this%visc   (i,j  ,k)*this%grdv_y( 0,i,j  ,k)+&
               &                                                                this%divv_y(-1,i,j,k)*2.0_WP*this%visc   (i,j-1,k)*this%grdv_y(+1,i,j-1,k)+&
               &                                                                this%divv_z(+1,i,j,k)*       this%visc_yz(i,j,k+1)*this%grdv_z(-1,i,j,k+1)+&
               &                                                                this%divv_z( 0,i,j,k)*       this%visc_yz(i,j,k  )*this%grdv_z( 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divv_x(+1,i,j,k)*       this%visc_xy(i+1,j,k)*this%grdv_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divv_x( 0,i,j,k)*       this%visc_xy(i  ,j,k)*this%grdv_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divv_y( 0,i,j,k)*2.0_WP*this%visc   (i,j  ,k)*this%grdv_y(+1,i,j  ,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divv_y(-1,i,j,k)*2.0_WP*this%visc   (i,j-1,k)*this%grdv_y( 0,i,j-1,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divv_z(+1,i,j,k)*       this%visc_yz(i,j,k+1)*this%grdv_z( 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divv_z( 0,i,j,k)*       this%visc_yz(i,j,k  )*this%grdv_z(-1,i,j,k  ))
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resV
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resV=this%implicit%sol
      
      ! Solve implicit W problem
      this%implicit%opr(1,:,:,:)=sum(this%itpr_z(:,i,j,k)*this%rho(i,j,k-1:k)); this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=sum(this%itpu_z(:,i+1,j,k)*this%rhoU(i+1,j,k-1:k))
               rhoUm=sum(this%itpu_z(:,i  ,j,k)*this%rhoU(i  ,j,k-1:k))
               rhoVp=sum(this%itpv_z(:,i,j+1,k)*this%rhoV(i,j+1,k-1:k))
               rhoVm=sum(this%itpv_z(:,i,j  ,k)*this%rhoV(i,j  ,k-1:k))
               rhoWp=sum(this%itpw_z(:,i,j,k  )*this%rhoW(i,j,k  :k+1))
               rhoWm=sum(this%itpw_z(:,i,j,k-1)*this%rhoW(i,j,k-1:k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+0.5_WP*dt*(this%divw_x(+1,i,j,k)*this%itpw_x(-1,i+1,j,k)*rhoUp+&
               &                                                                this%divw_x( 0,i,j,k)*this%itpw_x( 0,i  ,j,k)*rhoUm+&
               &                                                                this%divw_y(+1,i,j,k)*this%itpw_y(-1,i,j+1,k)*rhoVp+&
               &                                                                this%divw_y( 0,i,j,k)*this%itpw_y( 0,i,j  ,k)*rhoVm+&
               &                                                                this%divw_z( 0,i,j,k)*this%itpw_z( 0,i,j,k  )*rhoWp+&
               &                                                                this%divw_z(-1,i,j,k)*this%itpw_z(+1,i,j,k-1)*rhoWm)
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+0.5_WP*dt*(this%divw_x(+1,i,j,k)*this%itpw_x( 0,i+1,j,k)*rhoUp)
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+0.5_WP*dt*(this%divw_x( 0,i,j,k)*this%itpw_x(-1,i  ,j,k)*rhoUm)
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+0.5_WP*dt*(this%divw_y(+1,i,j,k)*this%itpw_y( 0,i,j+1,k)*rhoVp)
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+0.5_WP*dt*(this%divw_y( 0,i,j,k)*this%itpw_y(-1,i,j  ,k)*rhoVm)
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+0.5_WP*dt*(this%divw_z( 0,i,j,k)*this%itpw_z(+1,i,j,k  )*rhoWp)
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+0.5_WP*dt*(this%divw_z(-1,i,j,k)*this%itpw_z( 0,i,j,k-1)*rhoWm)
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divw_x(+1,i,j,k)*       this%visc_zx(i+1,j,k)*this%grdw_x(-1,i+1,j,k)+&
               &                                                                this%divw_x( 0,i,j,k)*       this%visc_zx(i  ,j,k)*this%grdw_x( 0,i  ,j,k)+&
               &                                                                this%divw_y(+1,i,j,k)*       this%visc_yz(i,j+1,k)*this%grdw_y(-1,i,j+1,k)+&
               &                                                                this%divw_y( 0,i,j,k)*       this%visc_yz(i,j  ,k)*this%grdw_y( 0,i,j  ,k)+&
               &                                                                this%divw_z( 0,i,j,k)*2.0_WP*this%visc   (i,j,k  )*this%grdw_z( 0,i,j,k  )+&
               &                                                                this%divw_z(-1,i,j,k)*2.0_WP*this%visc   (i,j,k-1)*this%grdw_z(+1,i,j,k-1))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divw_x(+1,i,j,k)*       this%visc_zx(i+1,j,k)*this%grdw_x( 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divw_x( 0,i,j,k)*       this%visc_zx(i  ,j,k)*this%grdw_x(-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divw_y(+1,i,j,k)*       this%visc_yz(i,j+1,k)*this%grdw_y( 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divw_y( 0,i,j,k)*       this%visc_yz(i,j  ,k)*this%grdw_y(-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divw_z( 0,i,j,k)*2.0_WP*this%visc   (i,j,k  )*this%grdw_z(+1,i,j,k  ))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divw_z(-1,i,j,k)*2.0_WP*this%visc   (i,j,k-1)*this%grdw_z( 0,i,j,k-1))
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resW
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resW=this%implicit%sol
      
      ! Sync up all residuals
      call this%cfg%sync(resU)
      call this%cfg%sync(resV)
      call this%cfg%sync(resW)
      
   end subroutine solve_implicit
   

   !> Update density from VFS object
   subroutine update_density(this,dt,vf)
      use vfs_class, only: vfs
      use irl_fortran_interface, only: getVolumePtr
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), intent(in) :: dt
      class(vfs), intent(inout) :: vf
      integer :: i,j,k
      ! Update density and momentum from VFS object
      do k=vf%cfg%kmino_,vf%cfg%kmaxo_
         do j=vf%cfg%jmino_,vf%cfg%jmaxo_
            do i=vf%cfg%imino_,vf%cfg%imaxo_
               this%rho(i,j,k)=this%rho_l*vf%VF(i,j,k)+this%rho_g*(1.0_WP-vf%VF(i,j,k))
               this%rhoU(i,j,k)=this%rho_l*getVolumePtr(vf%face_flux(1,i,j,k),0)/(vf%cfg%dy(j)*vf%cfg%dz(k)*dt)+this%rho_g*getVolumePtr(vf%face_flux(1,i,j,k),1)/(vf%cfg%dy(j)*vf%cfg%dz(k)*dt)
               this%rhoV(i,j,k)=this%rho_l*getVolumePtr(vf%face_flux(2,i,j,k),0)/(vf%cfg%dz(k)*vf%cfg%dx(i)*dt)+this%rho_g*getVolumePtr(vf%face_flux(2,i,j,k),1)/(vf%cfg%dz(k)*vf%cfg%dx(i)*dt)
               this%rhoW(i,j,k)=this%rho_l*getVolumePtr(vf%face_flux(3,i,j,k),0)/(vf%cfg%dx(i)*vf%cfg%dy(j)*dt)+this%rho_g*getVolumePtr(vf%face_flux(3,i,j,k),1)/(vf%cfg%dx(i)*vf%cfg%dy(j)*dt)
            end do
         end do
      end do
      ! Synchronize boundaries
      call this%cfg%sync(this%rho)
      call this%cfg%sync(this%rhoU)
      call this%cfg%sync(this%rhoV)
      call this%cfg%sync(this%rhoW)
   end subroutine update_density

   
   !> Compute face densities by interpolation
   subroutine get_face_density(this)
      implicit none
      class(tpns), intent(inout) :: this
      integer :: i,j,k
      ! Calculate rho_U/V/W using interpolation
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%rho_U(i,j,k)=sum(this%itpr_x(:,i,j,k)*this%rho(i-1:i,j,k))
            end do
         end do
      end do
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%rho_V(i,j,k)=sum(this%itpr_y(:,i,j,k)*this%rho(i,j-1:j,k))
            end do
         end do
      end do
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%rho_W(i,j,k)=sum(this%itpr_z(:,i,j,k)*this%rho(i,j,k-1:k))
            end do
         end do
      end do
      ! Handle non-periodic borders
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%rho_U(this%cfg%imino,:,:)=this%rho(this%cfg%imino,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%rho_V(:,this%cfg%jmino,:)=this%rho(:,this%cfg%jmino,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%rho_W(:,:,this%cfg%kmino)=this%rho(:,:,this%cfg%kmino)
      ! Synchronize boundaries
      call this%cfg%sync(this%rho_U)
      call this%cfg%sync(this%rho_V)
      call this%cfg%sync(this%rho_W)
   end subroutine get_face_density
   
   
   !> Prepare viscosity arrays from vfs object
   subroutine get_viscosity(this,vf,strat)
      use vfs_class, only: vfs
      use messager,  only: die
      implicit none
      class(tpns), intent(inout) :: this
      class(vfs), intent(in) :: vf
      integer :: i,j,k,mystrat
      real(WP) :: liq_vol,gas_vol,tot_vol
      integer, optional :: strat
      ! Choose viscosity averaging strategy
      if (present(strat)) then
         mystrat=strat
      else
         mystrat=harmonic_visc
      end if
      ! Check what strategy should be used
      select case (mystrat)
      case (harmonic_visc)
         ! Compute harmonically-averaged staggered viscosities using subcell phasic volumes
         do k=this%cfg%kmino_+1,this%cfg%kmaxo_
            do j=this%cfg%jmino_+1,this%cfg%jmaxo_
               do i=this%cfg%imino_+1,this%cfg%imaxo_
                  ! VISC at [xm,ym,zm] - direct sum in x/y/z
                  liq_vol=sum(vf%Lvol(:,:,:,i,j,k))
                  gas_vol=sum(vf%Gvol(:,:,:,i,j,k))
                  tot_vol=gas_vol+liq_vol
                  this%visc(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
                  ! VISC_xy at [x,y,zm] - direct sum in z, staggered sum in x/y
                  liq_vol=sum(vf%Lvol(0,0,:,i,j,k))+sum(vf%Lvol(1,0,:,i-1,j,k))+sum(vf%Lvol(0,1,:,i,j-1,k))+sum(vf%Lvol(1,1,:,i-1,j-1,k))
                  gas_vol=sum(vf%Gvol(0,0,:,i,j,k))+sum(vf%Gvol(1,0,:,i-1,j,k))+sum(vf%Gvol(0,1,:,i,j-1,k))+sum(vf%Gvol(1,1,:,i-1,j-1,k))
                  tot_vol=gas_vol+liq_vol
                  this%visc_xy(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc_xy(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
                  ! VISC_yz at [xm,y,z] - direct sum in x, staggered sum in y/z
                  liq_vol=sum(vf%Lvol(:,0,0,i,j,k))+sum(vf%Lvol(:,1,0,i,j-1,k))+sum(vf%Lvol(:,0,1,i,j,k-1))+sum(vf%Lvol(:,1,1,i,j-1,k-1))
                  gas_vol=sum(vf%Gvol(:,0,0,i,j,k))+sum(vf%Gvol(:,1,0,i,j-1,k))+sum(vf%Gvol(:,0,1,i,j,k-1))+sum(vf%Gvol(:,1,1,i,j-1,k-1))
                  tot_vol=gas_vol+liq_vol
                  this%visc_yz(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc_yz(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
                  ! VISC_zx at [x,ym,z] - direct sum in y, staggered sum in z/x
                  liq_vol=sum(vf%Lvol(0,:,0,i,j,k))+sum(vf%Lvol(0,:,1,i,j,k-1))+sum(vf%Lvol(1,:,0,i-1,j,k))+sum(vf%Lvol(1,:,1,i-1,j,k-1))
                  gas_vol=sum(vf%Gvol(0,:,0,i,j,k))+sum(vf%Gvol(0,:,1,i,j,k-1))+sum(vf%Gvol(1,:,0,i-1,j,k))+sum(vf%Gvol(1,:,1,i-1,j,k-1))
                  tot_vol=gas_vol+liq_vol
                  this%visc_zx(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc_zx(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
               end do
            end do
         end do
      case (arithmetic_visc)
         ! Compute arithmetically-averaged staggered viscosities using subcell phasic volumes
         do k=this%cfg%kmino_+1,this%cfg%kmaxo_
            do j=this%cfg%jmino_+1,this%cfg%jmaxo_
               do i=this%cfg%imino_+1,this%cfg%imaxo_
                  ! VISC at [xm,ym,zm] - direct sum in x/y/z
                  liq_vol=sum(vf%Lvol(:,:,:,i,j,k))
                  gas_vol=sum(vf%Gvol(:,:,:,i,j,k))
                  tot_vol=gas_vol+liq_vol
                  this%visc(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc(i,j,k)=(this%visc_l*liq_vol+this%visc_g*gas_vol)/tot_vol
                  ! VISC_xy at [x,y,zm] - direct sum in z, staggered sum in x/y
                  liq_vol=sum(vf%Lvol(0,0,:,i,j,k))+sum(vf%Lvol(1,0,:,i-1,j,k))+sum(vf%Lvol(0,1,:,i,j-1,k))+sum(vf%Lvol(1,1,:,i-1,j-1,k))
                  gas_vol=sum(vf%Gvol(0,0,:,i,j,k))+sum(vf%Gvol(1,0,:,i-1,j,k))+sum(vf%Gvol(0,1,:,i,j-1,k))+sum(vf%Gvol(1,1,:,i-1,j-1,k))
                  tot_vol=gas_vol+liq_vol
                  this%visc_xy(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc_xy(i,j,k)=(this%visc_l*liq_vol+this%visc_g*gas_vol)/tot_vol
                  ! VISC_yz at [xm,y,z] - direct sum in x, staggered sum in y/z
                  liq_vol=sum(vf%Lvol(:,0,0,i,j,k))+sum(vf%Lvol(:,1,0,i,j-1,k))+sum(vf%Lvol(:,0,1,i,j,k-1))+sum(vf%Lvol(:,1,1,i,j-1,k-1))
                  gas_vol=sum(vf%Gvol(:,0,0,i,j,k))+sum(vf%Gvol(:,1,0,i,j-1,k))+sum(vf%Gvol(:,0,1,i,j,k-1))+sum(vf%Gvol(:,1,1,i,j-1,k-1))
                  tot_vol=gas_vol+liq_vol
                  this%visc_yz(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc_yz(i,j,k)=(this%visc_l*liq_vol+this%visc_g*gas_vol)/tot_vol
                  ! VISC_zx at [x,ym,z] - direct sum in y, staggered sum in z/x
                  liq_vol=sum(vf%Lvol(0,:,0,i,j,k))+sum(vf%Lvol(0,:,1,i,j,k-1))+sum(vf%Lvol(1,:,0,i-1,j,k))+sum(vf%Lvol(1,:,1,i-1,j,k-1))
                  gas_vol=sum(vf%Gvol(0,:,0,i,j,k))+sum(vf%Gvol(0,:,1,i,j,k-1))+sum(vf%Gvol(1,:,0,i-1,j,k))+sum(vf%Gvol(1,:,1,i-1,j,k-1))
                  tot_vol=gas_vol+liq_vol
                  this%visc_zx(i,j,k)=0.0_WP
                  if (tot_vol.gt.0.0_WP) this%visc_zx(i,j,k)=(this%visc_l*liq_vol+this%visc_g*gas_vol)/tot_vol
               end do
            end do
         end do
      case default
         call die('[tpns get_viscosity] Unknown viscosity averaging strategy')
      end select
      ! Synchronize boundaries - not really needed...
      call this%cfg%sync(this%visc)
      call this%cfg%sync(this%visc_xy)
      call this%cfg%sync(this%visc_yz)
      call this%cfg%sync(this%visc_zx)
   end subroutine get_viscosity
   
   
   !> Add gravity source term - assumes that rho has been updated before
   subroutine addsrc_gravity(this,resU,resV,resW)
      implicit none
      class(tpns), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+this%rho_U(i,j,k)*this%gravity(1)
               if (this%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+this%rho_V(i,j,k)*this%gravity(2)
               if (this%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+this%rho_W(i,j,k)*this%gravity(3)
            end do
         end do
      end do
   end subroutine addsrc_gravity
   
   
   !> Add a static contact line model
   subroutine add_static_contact(this,vf)
      use mathtools, only: normalize
      use vfs_class, only: vfs
      use irl_fortran_interface
      implicit none
      class(tpns), intent(inout) :: this
      class(vfs),  intent(in) :: vf
      integer :: i,j,k
      real(WP), dimension(3) :: nw
      real(WP), dimension(2) :: fvof
      real(WP), dimension(:,:,:), allocatable :: GFM
      real(WP) :: dd,mysurf,mycos
      real(WP) :: cos_contact_angle
      real(WP) :: sin_contact_angle
      real(WP) :: tan_contact_angle
      real(WP), parameter :: cfactor=1.0_WP
      
      ! Allocate and zero out binarized VF for GFM-style jump distribution
      allocate(GFM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); GFM=0.0_WP
      
      ! Prepare a GFM-based strategy
      GFM=real(nint(vf%VF),WP)
      
      ! Precalculate cos/sin/tan(contact angle)
      cos_contact_angle=cos(this%contact_angle)
      sin_contact_angle=sin(this%contact_angle)
      tan_contact_angle=tan(this%contact_angle)
      
      ! Loop over domain and identify cells that require contact angle model
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! Check if we have an interface in the vicinity of the x-face
               mysurf=abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
               if (mysurf.gt.0.0_WP) then
                  ! Compute the liquid area fractions from GFM
                  fvof=GFM(i-1:i,j,k)
                  ! Check for local wall configuration - wall in y-
                  if (this%umask(i,j,k).eq.0.and.this%mask(i,j-1,k).eq.1.and.this%mask(i-1,j-1,k).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,+1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     ! Compute the liquid area fractions from PLIC
                     !call getMoments(vf%polyface(2,i-1,j,k),vf%liquid_gas_interface(i-1,j,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(2,i-1,j,k)))
                     !call getMoments(vf%polyface(2,i  ,j,k),vf%liquid_gas_interface(i  ,j,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(2,i  ,j,k)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i-1,j,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i  ,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i  ,j,k)),nw))/mysurf
                     ! Add source term
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in y+
                  if (this%umask(i,j,k).eq.0.and.this%mask(i,j+1,k).eq.1.and.this%mask(i-1,j+1,k).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,-1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(2,i-1,j+1,k),vf%liquid_gas_interface(i-1,j,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(2,i-1,j+1,k)))
                     !call getMoments(vf%polyface(2,i  ,j+1,k),vf%liquid_gas_interface(i  ,j,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(2,i  ,j+1,k)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i-1,j,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i  ,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i  ,j,k)),nw))/mysurf
                     ! Add source term
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in z-
                  if (this%umask(i,j,k).eq.0.and.this%mask(i,j,k-1).eq.1.and.this%mask(i-1,j,k-1).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,0.0_WP,+1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(3,i-1,j,k),vf%liquid_gas_interface(i-1,j,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(3,i-1,j,k)))
                     !call getMoments(vf%polyface(3,i  ,j,k),vf%liquid_gas_interface(i  ,j,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(3,i  ,j,k)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i-1,j,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i  ,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i  ,j,k)),nw))/mysurf
                     ! Add source term
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in z+
                  if (this%umask(i,j,k).eq.0.and.this%mask(i,j,k+1).eq.1.and.this%mask(i-1,j,k+1).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,0.0_WP,-1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(3,i-1,j,k+1),vf%liquid_gas_interface(i-1,j,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(3,i-1,j,k+1)))
                     !call getMoments(vf%polyface(3,i  ,j,k+1),vf%liquid_gas_interface(i  ,j,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(3,i  ,j,k+1)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i-1,j,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i  ,j,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i  ,j,k)),nw))/mysurf
                     ! Add source term
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
               end if
               
               ! Check if we have an interface in the vicinity of the y-face
               mysurf=abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
               if (mysurf.gt.0.0_WP) then
                  ! Compute the liquid area fractions from GFM
                  fvof=GFM(i,j-1:j,k)
                  ! Check for local wall configuration - wall in x-
                  if (this%vmask(i,j,k).eq.0.and.this%mask(i-1,j,k).eq.1.and.this%mask(i-1,j-1,k).eq.1) then
                     ! Define wall
                     nw=[+1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(1,i,j-1,k),vf%liquid_gas_interface(i,j-1,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(1,i,j-1,k)))
                     !call getMoments(vf%polyface(1,i,j  ,k),vf%liquid_gas_interface(i,j  ,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(1,i,j  ,k)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j-1,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j  ,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j  ,k)),nw))/mysurf
                     ! Add source term
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in x+
                  if (this%vmask(i,j,k).eq.0.and.this%mask(i+1,j,k).eq.1.and.this%mask(i+1,j-1,k).eq.1) then
                     ! Define wall
                     nw=[-1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(1,i+1,j-1,k),vf%liquid_gas_interface(i,j-1,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(1,i+1,j-1,k)))
                     !call getMoments(vf%polyface(1,i+1,j  ,k),vf%liquid_gas_interface(i,j  ,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(1,i+1,j  ,k)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j-1,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j  ,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j  ,k)),nw))/mysurf
                     ! Add source term
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in z-
                  if (this%vmask(i,j,k).eq.0.and.this%mask(i,j,k-1).eq.1.and.this%mask(i,j-1,k-1).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,0.0_WP,+1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(3,i,j-1,k),vf%liquid_gas_interface(i,j-1,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(3,i,j-1,k)))
                     !call getMoments(vf%polyface(3,i,j  ,k),vf%liquid_gas_interface(i,j  ,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(3,i,j  ,k)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j-1,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j  ,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j  ,k)),nw))/mysurf
                     ! Add source term
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in z+
                  if (this%vmask(i,j,k).eq.0.and.this%mask(i,j,k+1).eq.1.and.this%mask(i,j-1,k+1).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,0.0_WP,-1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(3,i,j-1,k+1),vf%liquid_gas_interface(i,j-1,k),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(3,i,j-1,k+1)))
                     !call getMoments(vf%polyface(3,i,j  ,k+1),vf%liquid_gas_interface(i,j  ,k),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(3,i,j  ,k+1)))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j-1,k)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j  ,k)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j  ,k)),nw))/mysurf
                     ! Add source term
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
               end if
               
               ! Check if we have an interface in the vicinity of the z-face
               mysurf=abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
               if (mysurf.gt.0.0_WP) then
                  ! Compute the liquid area fractions from GFM
                  fvof=GFM(i,j,k-1:k)
                  ! Check for local wall configuration - wall in x-
                  if (this%wmask(i,j,k).eq.0.and.this%mask(i-1,j,k).eq.1.and.this%mask(i-1,j,k-1).eq.1) then
                     ! Define wall
                     nw=[+1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(1,i,j,k-1),vf%liquid_gas_interface(i,j,k-1),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(1,i,j,k-1)))
                     !call getMoments(vf%polyface(1,i,j,k  ),vf%liquid_gas_interface(i,j,k  ),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(1,i,j,k  )))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k-1)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j,k  )))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k  )),nw))/mysurf
                     ! Add source term
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in x+
                  if (this%wmask(i,j,k).eq.0.and.this%mask(i+1,j,k).eq.1.and.this%mask(i+1,j,k-1).eq.1) then
                     ! Define wall
                     nw=[-1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(1,i+1,j,k-1),vf%liquid_gas_interface(i,j,k-1),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(1,i+1,j,k-1)))
                     !call getMoments(vf%polyface(1,i+1,j,k  ),vf%liquid_gas_interface(i,j,k  ),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(1,i+1,j,k  )))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k-1)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j,k  )))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k  )),nw))/mysurf
                     ! Add source term
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in y-
                  if (this%wmask(i,j,k).eq.0.and.this%mask(i,j-1,k).eq.1.and.this%mask(i,j-1,k-1).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,+1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(2,i,j,k-1),vf%liquid_gas_interface(i,j,k-1),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(2,i,j,k-1)))
                     !call getMoments(vf%polyface(2,i,j,k  ),vf%liquid_gas_interface(i,j,k  ),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(2,i,j,k  )))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k-1)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j,k  )))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k  )),nw))/mysurf
                     ! Add source term
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
                  ! Check for local wall configuration - wall in y+
                  if (this%wmask(i,j,k).eq.0.and.this%mask(i,j+1,k).eq.1.and.this%mask(i,j+1,k-1).eq.1) then
                     ! Define wall
                     nw=[0.0_WP,-1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     ! Compute the liquid area fractions
                     !call getMoments(vf%polyface(2,i,j+1,k-1),vf%liquid_gas_interface(i,j,k-1),fvof(1)); fvof(1)=abs(fvof(1))/abs(calculateVolume(vf%polyface(2,i,j+1,k-1)))
                     !call getMoments(vf%polyface(2,i,j+1,k  ),vf%liquid_gas_interface(i,j,k  ),fvof(2)); fvof(2)=abs(fvof(2))/abs(calculateVolume(vf%polyface(2,i,j+1,k  )))
                     ! Surface-averaged local cos(CA)
                     mycos=(abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k-1)),nw)+&
                     &      abs(calculateVolume(vf%interface_polygon(1,i,j,k  )))*dot_product(calculateNormal(vf%interface_polygon(1,i,j,k  )),nw))/mysurf
                     ! Add source term
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*fvof(:))*(mycos-cos_contact_angle)/dd
                  end if
               end if
               
            end do
         end do
      end do
      
      ! Deallocate array
      deallocate(GFM)
      
   end subroutine add_static_contact
   
   
   !> Print out info for two-phase incompressible flow solver
   subroutine tpns_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(tpns), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Two-phase incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   liquid density = ",es12.5)') this%rho_l
         write(output_unit,'(" > liquid viscosity = ",es12.5)') this%visc_l
         write(output_unit,'(" >      gas density = ",es12.5)') this%rho_g
         write(output_unit,'(" >    gas viscosity = ",es12.5)') this%visc_g
      end if
      
   end subroutine tpns_print
   
   
end module tpns_class
