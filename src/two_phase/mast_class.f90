!> Two-phase compressible flow solver class:
!> Provides support for ...
!> Interface is represented using VOF
module mast_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ils_class,      only: ils
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: mast,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: wall=1              !< Dirichlet at zero condition
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   !integer, parameter, public :: convective=4        !< Convective outflow condition
   integer, parameter, public :: clipped_neumann=5   !< Clipped Neumann condition (outflow only)
   
   ! List of available contact line models for this solver
   integer, parameter, public :: static_contact=1    !< Static contact line model
   
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
   
   !> Two-phase compressible solver object definition ([M]ultiphase [A]ll-Mach [S]emi-Lagrangian [T]ransport)
   type :: mast
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_MAST'    !< Solver name (default=UNNAMED_MAST)
      
      ! Constant fluid properties
      real(WP) :: contact_angle                           !< This is our static contact angle
      real(WP) :: sigma                                   !< This is our constant surface tension coefficient

      ! Variable or constant fluid properties
      real(WP) :: visc_l0,visc_g0                         !< These are our constant/initial dynamic viscosities in liquid and gas
      real(WP) :: cv_l0,cv_g0                             !< These are our constant/initial specific heats (constant volume) in liquid and gas
      real(WP) :: kappa_l0,kappa_g0                       !< These are our constant/initial thermal conductivities in liquid and gas
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP            !< Acceleration of gravity
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      real(WP), dimension(:), allocatable :: mfr          !< MFR through each bcond
      real(WP), dimension(:), allocatable :: area         !< Area for each bcond
      real(WP) :: correctable_area                        !< Area of bcond that can be corrected
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver

      ! Note: naming convention of "U-cell", etc. is maintained, but this is a collocated solver.
      ! U-cells correspond to face velocities, which are interpolated and projected quantities.
      ! The cell-centered velocities, such as "Ui", correspond to the conserved mass and momentum and are not interpolated.
      
      ! Interpolated density fields
      real(WP), dimension(:,:,:), allocatable :: rho_U    !< Density field array on U-cell
      real(WP), dimension(:,:,:), allocatable :: rho_V    !< Density field array on V-cell
      real(WP), dimension(:,:,:), allocatable :: rho_W    !< Density field array on W-cell
      
      ! Viscosity fields
      real(WP), dimension(:,:,:), allocatable :: visc     !< Viscosity field on P-cell
      real(WP), dimension(:,:,:), allocatable :: visc_xy  !< Viscosity field on U-cell
      real(WP), dimension(:,:,:), allocatable :: visc_yz  !< Viscosity field on V-cell
      real(WP), dimension(:,:,:), allocatable :: visc_zx  !< Viscosity field on W-cell
      
      ! Flow variables - harmonized/mixture
      real(WP), dimension(:,:,:), allocatable :: RHO      !< density array
      real(WP), dimension(:,:,:), allocatable :: RHOSS2   !< bulk modulus array
      real(WP), dimension(:,:,:), allocatable :: rhoUi    !< U momentum array
      real(WP), dimension(:,:,:), allocatable :: rhoVi    !< V momentum array
      real(WP), dimension(:,:,:), allocatable :: rhoWi    !< W momentum array
      real(WP), dimension(:,:,:), allocatable :: U        !< U face-velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V face-velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W face-velocity array
      real(WP), dimension(:,:,:), allocatable :: Ui       !< U center-velocity array
      real(WP), dimension(:,:,:), allocatable :: Vi       !< V center-velocity array
      real(WP), dimension(:,:,:), allocatable :: Wi       !< W center-velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      real(WP), dimension(:,:,:), allocatable :: PA       !< Advected pressure array
      real(WP), dimension(:,:,:), allocatable :: Pjx      !< Pressure jump to add to -dP/dx
      real(WP), dimension(:,:,:), allocatable :: Pjy      !< Pressure jump to add to -dP/dy
      real(WP), dimension(:,:,:), allocatable :: Pjz      !< Pressure jump to add to -dP/dz
      real(WP), dimension(:,:,:), allocatable :: dPjx     !< dPressure jump to add to -ddP/dx
      real(WP), dimension(:,:,:), allocatable :: dPjy     !< dPressure jump to add to -ddP/dy
      real(WP), dimension(:,:,:), allocatable :: dPjz     !< dPressure jump to add to -ddP/dz
      ! Flow variables - individual phases
      real(WP), dimension(:,:,:), allocatable :: Grho,   Lrho    !< phase density arrays
      real(WP), dimension(:,:,:), allocatable :: GrhoU,  LrhoU   !< phase U momentum arrays
      real(WP), dimension(:,:,:), allocatable :: GrhoV,  LrhoV   !< phase V momentum arrays
      real(WP), dimension(:,:,:), allocatable :: GrhoW,  LrhoW   !< phase W momentum arrays
      real(WP), dimension(:,:,:), allocatable :: GrhoE,  LrhoE   !< phase energy arrays
      real(WP), dimension(:,:,:), allocatable :: GP,     LP      !< phase pressure arrays
      real(WP), dimension(:,:,:), allocatable :: GrhoSS2,LrhoSS2 !< phase bulk modulus arrays
      
      ! Old flow variables
      real(WP), dimension(:,:,:), allocatable :: RHOold
      real(WP), dimension(:,:,:), allocatable :: Uiold,   Viold,   Wiold
      real(WP), dimension(:,:,:), allocatable :: GrhoUold,GrhoVold,GrhoWold,LrhoUold,LrhoVold,LrhoWold
      real(WP), dimension(:,:,:), allocatable :: GrhoEold,LrhoEold
      real(WP), dimension(:,:,:), allocatable :: Grhoold, Lrhoold
      real(WP), dimension(:,:,:), allocatable :: GPold,   LPold
      
      ! Flux sum arrays
      real(WP), dimension(:,:,:), allocatable :: F_GrhoU,F_GrhoV,F_GrhoW !< Momentum fluxes (gas)
      real(WP), dimension(:,:,:), allocatable :: F_LrhoU,F_LrhoV,F_LrhoW !< Momentum fluxes (liquid)
      real(WP), dimension(:,:,:), allocatable :: F_Grho, F_Lrho          !< Density fluxes
      real(WP), dimension(:,:,:), allocatable :: F_GrhoE,F_LrhoE         !< Energy fluxes
      real(WP), dimension(:,:,:), allocatable :: F_GP,   F_LP            !< Pressure fluxes
      real(WP), dimension(:,:,:), allocatable :: F_VOL,  F_VF            !< Volume fluxes

      ! Hybrid advection
      integer,  dimension(:,:,:,:), allocatable :: sl_face ! < Flag for flux method switching
      ! Pressure relaxation
      real(WP), dimension(:,:,:), allocatable :: srcVF ! < Predicted volume exchange during advection

      ! Pressure solver
      type(ils) :: psolv                                  !< Iterative linear solver object for the pressure Helmholtz equation
      
      ! Implicit momentum solver
      type(ils) :: implicit                               !< Iterative linear solver object for an implicit prediction of the advection residual
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itpi_x,itpi_y,itpi_z   !< Interpolation for Ui/Vi/Wi
      real(WP), dimension(:,:,:,:), allocatable :: divp_x,divp_y,divp_z   !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: divu_x,divu_y,divu_z   !< Divergence for U-cell
      real(WP), dimension(:,:,:,:), allocatable :: divv_x,divv_y,divv_z   !< Divergence for V-cell
      real(WP), dimension(:,:,:,:), allocatable :: divw_x,divw_y,divw_z   !< Divergence for W-cell
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable ::  mask                     !< Integer array used for modifying P metrics
      integer, dimension(:,:,:), allocatable :: umask                     !< Integer array used for modifying U metrics
      integer, dimension(:,:,:), allocatable :: vmask                     !< Integer array used for modifying V metrics
      integer, dimension(:,:,:), allocatable :: wmask                     !< Integer array used for modifying W metrics
      
      ! CFL numbers
      real(WP) :: CFLst                                                   !< Surface tension CFL
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                                    !< Convective CFL numbers
      real(WP) :: CFLa_x,CFLa_y,CFLa_z                                    !< Acoustic CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                                    !< Viscous CFL numbers
      
      ! Monitoring quantities
      real(WP) :: Umax,Vmax,Wmax,Pmax,RHOmax                              !< Maximum velocity, pressure, density
      
   contains
      procedure :: print=>mast_print                      !< Output solver to the screen
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: add_static_contact                     !< Add static contact line model to surface tension jump

      ! For initialization of simulation
      procedure :: setup                                  !< Finish configuring the flow solver
      procedure :: add_bcond                              !< Add a boundary condition
      ! For beginning of timestep (before iterative loop)
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      ! For advection solve
      procedure :: advection_step
      ! For viscous/dissipative/body forces (will address later)
      ! For setting up pressure solve
      !procedure :: interp_vel                             !< Calculate interpolated velocity
      !procedure :: calcH_LHS                              !< Calculate the left-hand side of the pressure equation
      !procedure :: calcH_RHS                              !< Calculate the right-hand side of the pressure equation
      procedure :: add_surface_tension_jump               !< Add surface tension jump
      ! For pressure correction
      procedure :: get_pgrad                              !< Calculate pressure gradient
      ! For pressure relaxation

      ! Miscellaneous
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_max                                !< Calculate maximum field values
      
   end type mast
   
   
   !> Declare two-phase compressible solver constructor
   interface mast
      procedure constructor
   end interface mast
   
contains
   
   
   !> Default constructor for two-phase compressible flow solver
   function constructor(cfg,name) result(self)
      implicit none
      type(mast) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate flow variables
      allocate(self%RHO(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%RHO=0.0_WP
      allocate(self%RHOSS2(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%RHOSS2=0.0_WP
      allocate(self%rhoUi(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhoUi=0.0_WP
      allocate(self%rhoVi(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhoVi=0.0_WP
      allocate(self%rhoWi(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhoWi=0.0_WP
      allocate(self%U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%U=0.0_WP
      allocate(self%V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%V=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%W=0.0_WP
      allocate(self%Ui(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Ui=0.0_WP
      allocate(self%Vi(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Vi=0.0_WP
      allocate(self%Wi(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Wi=0.0_WP
      allocate(self%P(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P=0.0_WP
      allocate(self%PA(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%PA=0.0_WP
      allocate(self%Pjx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Pjx=0.0_WP
      allocate(self%Pjy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Pjy=0.0_WP
      allocate(self%Pjz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Pjz=0.0_WP
      allocate(self%dPjx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dPjx=0.0_WP
      allocate(self%dPjy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dPjy=0.0_WP
      allocate(self%dPjz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dPjz=0.0_WP
      allocate(self%rho_U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_U=0.0_WP
      allocate(self%rho_V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_V=0.0_WP
      allocate(self%rho_W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_W=0.0_WP
      allocate(self%visc   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc   =0.0_WP
      allocate(self%visc_xy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_xy=0.0_WP
      allocate(self%visc_yz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_yz=0.0_WP
      allocate(self%visc_zx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_zx=0.0_WP
      ! Two-phase
      allocate(self%Grho (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Grho =0.0_WP
      allocate(self%Lrho (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Lrho =0.0_WP
      allocate(self%GrhoU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoU=0.0_WP
      allocate(self%GrhoV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoV=0.0_WP
      allocate(self%GrhoW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoW=0.0_WP
      allocate(self%LrhoU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoU=0.0_WP
      allocate(self%LrhoV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoV=0.0_WP
      allocate(self%LrhoW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoW=0.0_WP
      allocate(self%GrhoE(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoE=0.0_WP
      allocate(self%LrhoE(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoE=0.0_WP
      allocate(self%GP   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GP   =0.0_WP
      allocate(self%LP   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LP   =0.0_WP

      
      ! Allocate old flow variables
      allocate(self%RHOold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%RHOold=0.0_WP
      allocate(self%Uiold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Uiold=0.0_WP
      allocate(self%Viold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Viold=0.0_WP
      allocate(self%Wiold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Wiold=0.0_WP
      allocate(self%GrhoUold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoUold=0.0_WP
      allocate(self%GrhoVold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoVold=0.0_WP
      allocate(self%GrhoWold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoWold=0.0_WP
      allocate(self%LrhoUold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoUold=0.0_WP
      allocate(self%LrhoVold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoVold=0.0_WP
      allocate(self%LrhoWold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoWold=0.0_WP
      allocate(self%GrhoEold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoEold=0.0_WP
      allocate(self%LrhoEold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoEold=0.0_WP
      allocate(self%GPold   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GPold   =0.0_WP
      allocate(self%LPold   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LPold   =0.0_WP
      
      ! Flux sum arrays need to be preallocated
      allocate(self%F_GrhoU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_GrhoU=0.0_WP
      allocate(self%F_GrhoV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_GrhoV=0.0_WP
      allocate(self%F_GrhoW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_GrhoW=0.0_WP
      allocate(self%F_LrhoU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_LrhoU=0.0_WP
      allocate(self%F_LrhoV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_LrhoV=0.0_WP
      allocate(self%F_LrhoW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_LrhoW=0.0_WP
      allocate(self%F_Grho (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_Grho =0.0_WP
      allocate(self%F_Lrho (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_Lrho =0.0_WP
      allocate(self%F_GrhoE(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_GrhoE=0.0_WP
      allocate(self%F_LrhoE(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_LrhoE=0.0_WP
      allocate(self%F_GP   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_GP   =0.0_WP
      allocate(self%F_LP   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_LP   =0.0_WP

      ! Hybrid advection
      allocate(self%sl_face(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3)); self%sl_face=0.0_WP
      ! Pressure relaxation
      allocate(self%srcVF(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcVF=0.0_WP
      
      ! Create pressure solver object
      self%psolv   =ils(cfg=self%cfg,name='Pressure')
      
      ! Create implicit velocity solver object
      self%implicit=ils(cfg=self%cfg,name='Momentum')
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare P-cell masks
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
      
      ! Prepare face mask for U
      allocate(self%umask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%umask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%umask(self%cfg%imin  ,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%umask(self%cfg%imax+1,:,:)=2
      end if
      do k=self%cfg%kmino_  ,self%cfg%kmaxo_
         do j=self%cfg%jmino_  ,self%cfg%jmaxo_
            do i=self%cfg%imino_+1,self%cfg%imaxo_
               if (minval(self%cfg%VF(i-1:i,j,k)).eq.0.0_WP) self%umask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%umask)
      if (.not.self%cfg%xper.and.self%cfg%iproc.eq.1) self%umask(self%cfg%imino,:,:)=self%umask(self%cfg%imino+1,:,:)
      
      ! Prepare face mask for V
      allocate(self%vmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%vmask=0
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%vmask(:,self%cfg%jmin  ,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%vmask(:,self%cfg%jmax+1,:)=2
      end if
      do k=self%cfg%kmino_  ,self%cfg%kmaxo_
         do j=self%cfg%jmino_+1,self%cfg%jmaxo_
            do i=self%cfg%imino_  ,self%cfg%imaxo_
               if (minval(self%cfg%VF(i,j-1:j,k)).eq.0.0_WP) self%vmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%vmask)
      if (.not.self%cfg%yper.and.self%cfg%jproc.eq.1) self%vmask(:,self%cfg%jmino,:)=self%vmask(:,self%cfg%jmino+1,:)
      
      ! Prepare face mask for W
      allocate(self%wmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%wmask=0
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%wmask(:,:,self%cfg%kmin  )=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%wmask(:,:,self%cfg%kmax+1)=2
      end if
      do k=self%cfg%kmino_+1,self%cfg%kmaxo_
         do j=self%cfg%jmino_  ,self%cfg%jmaxo_
            do i=self%cfg%imino_  ,self%cfg%imaxo_
               if (minval(self%cfg%VF(i,j,k-1:k)).eq.0.0_WP) self%wmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%wmask)
      if (.not.self%cfg%zper.and.self%cfg%kproc.eq.1) self%wmask(:,:,self%cfg%kmino)=self%wmask(:,:,self%cfg%kmino+1)
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      implicit none
      class(mast), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP), dimension(-1:0) :: itpx,itpy,itpz
      
      ! Allocate interpolation coefficients for center to cell face
      allocate(this%itpi_x(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itpi_y(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itpi_z(-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create interpolation coefficients to cell face
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%itpi_x(:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itpi_y(:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itpi_z(:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
            
      ! Allocate finite volume divergence operators
      allocate(this%divp_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      allocate(this%divp_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm] or tangent to cell face
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%divp_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divp_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divp_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do

      ! Divergence operators used in pressure gradients at faces
      allocate(this%divu_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (x)
      allocate(this%divv_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_+1:this%cfg%jmaxo_,this%cfg%kmino_  :this%cfg%kmaxo_)) !< Face-centered (y)
      allocate(this%divw_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_,this%cfg%jmino_  :this%cfg%jmaxo_,this%cfg%kmino_+1:this%cfg%kmaxo_)) !< Face-centered (z)
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
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls
   subroutine adjust_metrics(this)
      implicit none
      class(mast), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP) :: delta,mysum
      
      ! Sync up u/v/wmasks
      call this%cfg%sync(this%umask)
      call this%cfg%sync(this%vmask)
      call this%cfg%sync(this%wmask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%umask(this%cfg%imino,:,:)=this%umask(this%cfg%imino+1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%wmask(:,:,this%cfg%kmino)=this%wmask(:,:,this%cfg%kmino+1)
      
      ! Adjust Ui/Vi/Wi interpolation coefficients back to cell faces in the presence of walls (only walls!)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%umask(i,j,k).eq.1) this%itpi_x(:,i,j,k)=0.0_WP
               if (this%vmask(i,j,k).eq.1) this%itpi_y(:,i,j,k)=0.0_WP
               if (this%wmask(i,j,k).eq.1) this%itpi_z(:,i,j,k)=0.0_WP
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
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to V metrics
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%vmask(i,j,k).gt.0) then
                  this%divv_y(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to W metrics
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               if (this%wmask(i,j,k).gt.0) then
                  this%divw_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%divp_x=0.0_WP
         this%divu_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%divp_y=0.0_WP
         this%divv_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%divp_z=0.0_WP
         this%divw_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the flow solver now that bconds have been defined
   subroutine setup(this,pressure_ils,implicit_ils)
      implicit none
      class(mast), intent(inout) :: this
      integer, intent(in) :: pressure_ils
      integer, intent(in) :: implicit_ils
      
      ! Adjust metrics based on bcflag array
      call this%adjust_metrics()
      
      ! Set 7-pt stencil map for the pressure solver
      this%psolv%stc(1,:)=[ 0, 0, 0]
      this%psolv%stc(2,:)=[+1, 0, 0]
      this%psolv%stc(3,:)=[-1, 0, 0]
      this%psolv%stc(4,:)=[ 0,+1, 0]
      this%psolv%stc(5,:)=[ 0,-1, 0]
      this%psolv%stc(6,:)=[ 0, 0,+1]
      this%psolv%stc(7,:)=[ 0, 0,-1]
      
      ! Set the diagonal to VF to make sure all needed cells participate in solver
      this%psolv%opr(1,:,:,:)=this%cfg%VF
      
      ! Initialize the pressure Poisson solver
      call this%psolv%init(pressure_ils)
      
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
      call this%implicit%init(implicit_ils)
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,face,dir,canCorrect)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(mast), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: type
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
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
      case default; call die('[mast add_bcond] Unknown bcond face - expecting x, y, or z')
      end select
      new_bc%itr=iterator(pg=this%cfg,name=new_bc%name,locator=locator,face=new_bc%face)
      select case (dir) ! Outward-oriented
      case (+1); new_bc%dir=+1
      case (-1); new_bc%dir=-1
      case ( 0); new_bc%dir= 0
      case default; call die('[mast add_bcond] Unknown bcond dir - expecting -1, +1, or 0')
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
      !case (convective)
      case default
         call die('[mast apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(mast), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[mast get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      implicit none
      class(mast), intent(inout) :: this
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
               
            case (neumann,clipped_neumann) !< Apply Neumann condition to all 3 components
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
                  end do
               case ('y')
                  stag=min(my_bc%dir,0)
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j+stag,k    )=this%U(i:i+1,j-my_bc%dir+stag,k    )
                     this%V(i    ,j     ,k    )=this%V(i    ,j-my_bc%dir     ,k    )
                     this%W(i    ,j+stag,k:k+1)=this%W(i    ,j-my_bc%dir+stag,k:k+1)
                  end do
               case ('z')
                  stag=min(my_bc%dir,0)
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     this%U(i:i+1,j    ,k+stag)=this%U(i:i+1,j    ,k-my_bc%dir+stag)
                     this%V(i    ,j:j+1,k+stag)=this%V(i    ,j:j+1,k-my_bc%dir+stag)
                     this%W(i    ,j    ,k     )=this%W(i    ,j    ,k-my_bc%dir     )
                  end do
               end select
               ! If needed, clip
               if (my_bc%type.eq.clipped_neumann) then
                  select case (my_bc%face)
                  case ('x')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        if (this%U(i,j,k)*my_bc%rdir.lt.0.0_WP) this%U(i,j,k)=0.0_WP
                     end do
                  case ('y')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        if (this%V(i,j,k)*my_bc%rdir.lt.0.0_WP) this%V(i,j,k)=0.0_WP
                     end do
                  case ('z')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        if (this%W(i,j,k)*my_bc%rdir.lt.0.0_WP) this%W(i,j,k)=0.0_WP
                     end do
                  end select
               end if
               
            !case (convective)   ! Not implemented yet!
               
            case default
               call die('[mast apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%U)
         call this%cfg%sync(this%V)
         call this%cfg%sync(this%W)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond


   subroutine advection_step(this,dt,vf)
     use vfs_class, only: vfs, VFhi, VFlo
     implicit none
     class(mast), intent(inout) :: this   !< The two-phase flow solver
     class(vfs),  intent(inout) :: vf     !< The volume fraction solver
     real(WP),    intent(inout) :: dt     !< Timestep size over which to advance
     real(WP),   dimension(14)  :: flux   !< Passes flux to and from routines
     real(WP) :: Ga_i,Ga_nb,La_i,La_nb
     integer  :: i,j,k

     ! Initialize flux arrays
     this%F_VOL   = 0.0_WP; this%F_VF    = 0.0_WP
     this%F_Grho  = 0.0_WP; this%F_Lrho  = 0.0_WP
     this%F_GrhoU = 0.0_WP; this%F_LrhoU = 0.0_WP
     this%F_GrhoV = 0.0_WP; this%F_LrhoV = 0.0_WP
     this%F_GrhoW = 0.0_WP; this%F_LrhoW = 0.0_WP
     this%F_GrhoE = 0.0_WP; this%F_LrhoE = 0.0_WP
     this%F_GP    = 0.0_WP; this%F_LP    = 0.0_WP

     ! Initialize flux sum arrays
     this%F_VOL  =                                 this%cfg%vol
     this%F_VF   =                       vf%VFold *this%cfg%vol
     this%F_Grho =this%Grhoold *((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoE=this%GrhoEold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoU=this%GrhoUold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoV=this%GrhoVold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoW=this%GrhoWold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GP   =this%GPold   *((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_Lrho =this%Lrhoold *((       vf%VFold)*this%cfg%vol)
     this%F_LrhoE=this%LrhoEold*((       vf%VFold)*this%cfg%vol)
     this%F_LrhoU=this%LrhoUold*((       vf%VFold)*this%cfg%vol)
     this%F_LrhoV=this%LrhoVold*((       vf%VFold)*this%cfg%vol)
     this%F_LrhoW=this%LrhoWold*((       vf%VFold)*this%cfg%vol)
     this%F_LP   =this%LPold   *((       vf%VFold)*this%cfg%vol)

     !! ---------------------------------------!!
     !! 1. SL and TTSL flux calculations       !!
     !! ---------------------------------------!!
     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
              
              !! ---- LEFT X(I) FACE ---- !!
              select case(this%sl_face(i,j,k,1))
              case(1)!; call SL_advect  (flux,b_flux,fp,ffm              ,i,j,k,'x')
              case(0)!; call TTSL_advect(flux,b_flux,[x (i),ym(j),zm(k)],i,j,k,'x')
              end select
              call add_fluxes(flux,i,j,k,'x')

              !! ---- BOTTOM Y(J) FACE ---- !!
              select case(this%sl_face(i,j,k,2))
              case(1)!; call SL_advect  (flux,b_flux,fp,ffm              ,i,j,k,'y')
              case(0)!; call TTSL_advect(flux,b_flux,[xm(i),y (j),zm(k)],i,j,k,'y')
              end select
              call add_fluxes(flux,i,j,k,'y')

              !! ---- BACK Z(K) FACE ---- !!
              select case(this%sl_face(i,j,k,3))
              case(1)!; call SL_advect  (flux,b_flux,fp,ffm             ,i,j,k,'z')
              case(0)!; call TTSL_advect(flux,b_flux,[xm(i),ym(j),z (k)],i,j,k,'z')
              end select
              call add_fluxes(flux,i,j,k,'z')

           end do
        end do
     end do

     ! Zero the source term for volume fraction
     this%srcVF = 0.0_WP

     !! ---------------------------------------!!
     !! 2. VF, adv. pressure, and density      !!
     !! ---------------------------------------!!
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! Skip wall cells
              if (this%cfg%vol(i,j,k).eq.0.0_WP) cycle
              ! Update VF
              vf%VF(i,j,k) = this%F_VF(i,j,k)/this%F_VOL(i,j,k)
              ! Update phase density, pressure, and bulk moduli according to VF limits
              if (vf%VF(i,j,k).lt.VFlo) then
                 ! Gas only
                 vf%VF(i,j,k) = 0.0_WP
                 ! Calculate PA, bulkmod using primitively advected density
                 this%Grho (i,j,k)    = this%F_Grho (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%GP   (i,j,k)    = this%F_GP   (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 !this%GrhoSS2(i,j,k) = matm%EOS_gas(i,j,k,'K')
                 ! Zero quantities in opposite phase
                 this%Lrho (i,j,k)    = 0.0_WP
                 this%LP   (i,j,k)    = 0.0_WP
                 this%LrhoSS2(i,j,k) = 0.0_WP
                 ! Calculate density incorporating volume change
                 this%Grho (i,j,k)    = this%F_Grho (i,j,k)/this%cfg%vol(i,j,k)
              else if (vf%VF(i,j,k).gt.VFhi) then
                 ! Liquid only
                 vf%VF  (i,j,k) = 1.0_WP
                 ! Calculate PA, bulkmod using primitively advected density
                 this%Lrho (i,j,k)    = this%F_Lrho (i,j,k)/this%F_VF(i,j,k)
                 this%LP   (i,j,k)    = this%F_LP   (i,j,k)/this%F_VF(i,j,k)
                 !this%LrhoSS2(i,j,k) = matm%EOS_liq(i,j,k,'K')
                 ! Zero quantities in opposite phase
                 this%Grho (i,j,k)    = 0.0_WP
                 this%GP   (i,j,k)    = 0.0_WP
                 this%GrhoSS2(i,j,k) = 0.0_WP
                 ! Calculate density incorporating volume change
                 this%Lrho (i,j,k)    = this%F_Lrho (i,j,k)/this%cfg%vol(i,j,k)
              else
                 ! Primitive density, advected pressure, bulkmod for both phases
                 this%Grho (i,j,k)    = this%F_Grho (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%Lrho (i,j,k)    = this%F_Lrho (i,j,k)/(                  this%F_VF(i,j,k))
                 this%GP   (i,j,k)    = this%F_GP   (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%LP   (i,j,k)    = this%F_LP   (i,j,k)/(                  this%F_VF(i,j,k))
                 !this%GrhoSS2(i,j,k) = matm%EOS_gas(i,j,k,'K')
                 !this%LrhoSS2(i,j,k) = matm%EOS_liq(i,j,k,'K')
                 ! Store current VF
                 this%srcVF(i,j,k)    = vf%VF(i,j,k)
                 ! Get new VF according to quadratic source term that depends on compressibility
                 vf%VF(i,j,k) = VF_src_quad(this%GrhoSS2(i,j,k), this%LrhoSS2(i,j,k), &
                      this%F_VF(i,j,k), this%F_VOL(i,j,k), this%cfg%vol(i,j,k))
                 ! Update phase quantities according to VF limits (again)
                 if (vf%VF(i,j,k).lt.VFlo) then
                    vf%VF(i,j,k)        = 0.0_WP
                    this%srcVF(i,j,k)    = 0.0_WP
                    ! Zero density, PA, bulkmod, and dpde_rho of other phase
                    this%Lrho (i,j,k)    = 0.0_WP
                    this%LP   (i,j,k)    = 0.0_WP
                    this%LrhoSS2(i,j,k) = 0.0_WP
                    ! Calculate density, incorporating volume change
                    this%Grho (i,j,k) = this%F_Grho (i,j,k)/this%cfg%vol(i,j,k)
                 else if (vf%VF(i,j,k).gt.VFhi) then
                    vf%VF(i,j,k)        = 1.0_WP
                    this%srcVF(i,j,k)    = 0.0_WP
                    ! Zero density, PA, bulkmod, and dpde_rho of other phase
                    this%Grho (i,j,k)    = 0.0_WP
                    this%GP   (i,j,k)    = 0.0_WP
                    this%GrhoSS2(i,j,k) = 0.0_WP
                    ! Calculate density, incorporating volume change
                    this%Lrho (i,j,k) = this%F_Lrho (i,j,k)/this%cfg%vol(i,j,k)
                 else
                    ! Record srcVF
                    this%srcVF(i,j,k) = vf%VF(i,j,k)-this%srcVF(i,j,k)
                    ! PA, bulkmod stay unchanged
                    ! Calculate density, incorporating volume change
                    this%Grho (i,j,k) = this%F_Grho (i,j,k)/((1.0_WP-vf%VF(i,j,k))*this%cfg%vol(i,j,k))
                    this%Lrho (i,j,k) = this%F_Lrho (i,j,k)/(        vf%VF(i,j,k) *this%cfg%vol(i,j,k))
                 end if
              end if
           end do
        end do
     end do

     ! Boundaries for VF, phase density, phase pressure

     ! Mixture density
     this%RHO = vf%VF*this%Lrho + (1.0_WP-vf%VF)*this%Grho

     ! Initialize quantities
     this%implicit%opr(:,:,:,:) = 0.0_WP        ! zero operator
     this%implicit%opr(1,:,:,:) = this%cfg%vol  ! unity*volume diagonal
     
     !! ---------------------------------------!!
     !! 3. Implicit, centered momentum setup   !!
     !! ---------------------------------------!!
     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1

              !! ---- LEFT X(I) FACE ---- !!
              ! Advection term
              if (this%sl_face(i,j,k,1).eq.0) then
                 call SubCan_mom_RHS_flux(i,j,k,'x')

                 ! Calculate coefficients to go into alap array
                 call SubCan_LHS(i,j,k,'x',Ga_i,Ga_nb,La_i,La_nb)

                 if (.false.) then
                    ! Boundary conditions
                 else
                    ! Just like the fluxes, subtract from current cell (right from face), add to left cell (left from face)
                    ! Note that "neighbor" (nb) cell is left of current cell
                    ! Current cell (right from face)
                    this%implicit%opr(1,i  ,j,k) = this%implicit%opr(1,i  ,j,k) + (1.0_WP-vf%VF(i  ,j,k))*Ga_i + vf%VF(i  ,j,k)*La_i
                    this%implicit%opr(3,i  ,j,k) = this%implicit%opr(3,i  ,j,k) + (1.0_WP-vf%VF(i  ,j,k))*Ga_nb+ vf%VF(i  ,j,k)*La_nb
                    ! Left cell (left from face)
                    this%implicit%opr(1,i-1,j,k) = this%implicit%opr(1,i-1,j,k) + (1.0_WP-vf%VF(i-1,j,k))*Ga_i + vf%VF(i-1,j,k)*La_i
                    this%implicit%opr(2,i-1,j,k) = this%implicit%opr(2,i-1,j,k) + (1.0_WP-vf%VF(i-1,j,k))*Ga_nb+ vf%VF(i-1,j,k)*La_nb
                 end if
              end if

              !! ---- BOTTOM Y(J) FACE ---- !!
              ! Advection term
              if (this%sl_face(i,j,k,2).eq.0) then
                 call SubCan_mom_RHS_flux(i,j,k,'y')

                 ! Calculate coefficients to go into alap array
                 call SubCan_LHS(i,j,k,'y',Ga_i,Ga_nb,La_i,La_nb)

                 if (.false.) then
                    ! Boundary conditions
                 else
                    ! Current cell (above face)
                    this%implicit%opr(1,i,j  ,k) = this%implicit%opr(1,i,j  ,k) + (1.0_WP-vf%VF(i,j  ,k))*Ga_i + vf%VF(i,j  ,k)*La_i
                    this%implicit%opr(5,i,j  ,k) = this%implicit%opr(5,i,j  ,k) + (1.0_WP-vf%VF(i,j  ,k))*Ga_nb+ vf%VF(i,j  ,k)*La_nb
                    ! Left cell (below face)
                    this%implicit%opr(1,i,j-1,k) = this%implicit%opr(1,i,j-1,k) + (1.0_WP-vf%VF(i,j-1,k))*Ga_i + vf%VF(i,j-1,k)*La_i
                    this%implicit%opr(4,i,j-1,k) = this%implicit%opr(4,i,j-1,k) + (1.0_WP-vf%VF(i,j-1,k))*Ga_nb+ vf%VF(i,j-1,k)*La_nb
                 end if
              end if

              !! ---- BACK Z(K) FACE ---- !!
              ! Advection term
              if (this%sl_face(i,j,k,3).eq.0) then
                 call SubCan_mom_RHS_flux(i,j,k,'z')

                 ! Calculate coefficients to go into alap array
                 call SubCan_LHS(i,j,k,'z',Ga_i,Ga_nb,La_i,La_nb)

                 if (.false.) then
                    ! Boundary conditions
                 else
                    ! Current cell (in front of face)
                    this%implicit%opr(1,i,j,k  ) = this%implicit%opr(1,i,j,k  ) + (1.0_WP-vf%VF(i,j,k  ))*Ga_i + vf%VF(i,j,k  )*La_i
                    this%implicit%opr(7,i,j,k  ) = this%implicit%opr(5,i,j,k  ) + (1.0_WP-vf%VF(i,j,k  ))*Ga_nb+ vf%VF(i,j,k  )*La_nb
                    ! Left cell (behind face)
                    this%implicit%opr(1,i,j,k-1) = this%implicit%opr(1,i,j,k-1) + (1.0_WP-vf%VF(i,j,k-1))*Ga_i + vf%VF(i,j,k-1)*La_i
                    this%implicit%opr(6,i,j,k-1) = this%implicit%opr(6,i,j,k-1) + (1.0_WP-vf%VF(i,j,k-1))*Ga_nb+ vf%VF(i,j,k-1)*La_nb
                 end if
              end if

           end do
        end do
     end do

     ! Put momentum eq RHS into arrays
     this%rhoUi = (1.0_WP-vf%VF)*this%F_GrhoU + vf%VF*this%F_LrhoU
     this%rhoVi = (1.0_WP-vf%VF)*this%F_GrhoV + vf%VF*this%F_LrhoV
     this%rhoWi = (1.0_WP-vf%VF)*this%F_GrhoW + vf%VF*this%F_LrhoW

     ! Solve for x-momentum
     call this%implicit%setup()
     this%implicit%rhs=this%rhoUi
     this%implicit%sol=0.0_WP
     call this%implicit%solve()
     this%rhoUi=this%implicit%sol
     ! Solve for y-momentum
     call this%implicit%setup()
     this%implicit%rhs=this%rhoVi
     this%implicit%sol=0.0_WP
     call this%implicit%solve()
     this%rhoVi=this%implicit%sol
     ! Solve for z-momentum
     call this%implicit%setup()
     this%implicit%rhs=this%rhoWi
     this%implicit%sol=0.0_WP
     call this%implicit%solve()
     this%rhoWi=this%implicit%sol

     ! Boundary conditions for momentum

     ! Velocity
     this%Ui=this%rhoUi/this%RHO
     this%Vi=this%rhoVi/this%RHO
     this%Wi=this%rhoWi/this%RHO

     !! ---------------------------------------!!
     !! 4. Explicit KE fluxes and pres terms   !!
     !! ---------------------------------------!!
     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1

              !! ---- LEFT X(I) FACE ---- !!
              if (this%sl_face(i,j,k,1).eq.0) then
                 call SubCan_KE_RHS_flux(i,j,k,'x')
              end if

              !! ---- BOTTOM Y(J) FACE ---- !!
              if (this%sl_face(i,j,k,2).eq.0) then
                 call SubCan_KE_RHS_flux(i,j,k,'y')
              end if

              !! ---- BACK Z(K) FACE ---- !!
              if (this%sl_face(i,j,k,3).eq.0) then
                 call SubCan_KE_RHS_flux(i,j,k,'z')
              end if

           end do
        end do
     end do
     
     !! ---------------------------------------!!
     !! 5. Total energy calculation            !!
     !! ---------------------------------------!!
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! Cycle wall cells
              if (this%cfg%vol(i,j,k).eq.0.0_WP) cycle
              ! Update phase energy
              if      (vf%VF(i,j,k).lt.VFlo) then
                 this%GrhoE(i,j,k) = this%F_GrhoE(i,j,k)/this%cfg%vol(i,j,k)
                 this%LrhoE(i,j,k) = 0.0_WP
              else if (vf%VF(i,j,k).gt.VFhi) then
                 this%GrhoE(i,j,k) = 0.0_WP
                 this%LrhoE(i,j,k) = this%F_LrhoE(i,j,k)/this%cfg%vol(i,j,k)
              else
                 this%GrhoE(i,j,k)=this%F_GrhoE(i,j,k)/((1.0_WP-vf%VF(i,j,k))*this%cfg%vol(i,j,k))
                 this%LrhoE(i,j,k)=this%F_LrhoE(i,j,k)/(        vf%VF(i,j,k) *this%cfg%vol(i,j,k))
              end if
           end do
        end do
     end do

     ! Energy boundary conditions
  
   contains

     function VF_src_quad(gss2,lss2,fvf,fvl,vl) result(volfrac)
       real(WP), intent(in) :: gss2, lss2, fvf, fvl, vl
       real(WP) :: volfrac
       real(WP) :: quad_a, quad_b, quad_c, det
       ! Components of quadratic equation to solve for VF
       quad_a = (gss2-lss2)*vl
       quad_b = lss2*(vl+fvf)-gss2*(fvf+vl-fvl)
       quad_c = -lss2*fvf
       ! Solve for VF
       if (quad_a.eq.0.0_WP) then  ! Linear case
          volfrac = -quad_c/quad_b
       else ! Quadratic case, addition solution
          det = max(0.0_WP,quad_b**2-4.0_WP*quad_a*quad_c)
          volfrac = (-quad_b + sqrt(det))/(2.0_WP*quad_a)
       end if
     end function VF_src_quad

     subroutine add_fluxes(flux,i,j,k,dir)
       implicit none
       real(WP), dimension(14) :: flux
       character(len=1) :: dir
       integer :: i,j,k,st_i,st_j,st_k

       select case (trim(dir))
       case('x')
          st_i =-1; st_j = 0; st_k = 0
       case('y')
          st_i = 0; st_j =-1; st_k = 0
       case('z')
          st_i = 0; st_j = 0; st_k =-1
       end select

       ! Store update for right cell from face
       this%F_VOL  (i,j,k)=this%F_VOL  (i,j,k) +flux( 1)
       this%F_VF   (i,j,k)=this%F_VF   (i,j,k) +flux( 2)
       this%F_Grho (i,j,k)=this%F_Grho (i,j,k) +flux( 3)
       this%F_GrhoE(i,j,k)=this%F_GrhoE(i,j,k) +flux( 4)
       this%F_GrhoU(i,j,k)=this%F_GrhoU(i,j,k) +flux( 5)
       this%F_GrhoV(i,j,k)=this%F_GrhoV(i,j,k) +flux( 6)
       this%F_GrhoW(i,j,k)=this%F_GrhoW(i,j,k) +flux( 7)
       this%F_Lrho (i,j,k)=this%F_Lrho (i,j,k) +flux( 8)
       this%F_LrhoE(i,j,k)=this%F_LrhoE(i,j,k) +flux( 9)
       this%F_LrhoU(i,j,k)=this%F_LrhoU(i,j,k) +flux(10)
       this%F_LrhoV(i,j,k)=this%F_LrhoV(i,j,k) +flux(11)
       this%F_LrhoW(i,j,k)=this%F_LrhoW(i,j,k) +flux(12)
       this%F_GP   (i,j,k)=this%F_GP   (i,j,k) +flux(13)
       this%F_LP   (i,j,k)=this%F_LP   (i,j,k) +flux(14)
       ! Store update for left cell from face
       this%F_VOL  (i+st_i,j+st_j,k+st_k)=this%F_VOL  (i+st_i,j+st_j,k+st_k) -flux( 1)
       this%F_VF   (i+st_i,j+st_j,k+st_k)=this%F_VF   (i+st_i,j+st_j,k+st_k) -flux( 2)
       this%F_Grho (i+st_i,j+st_j,k+st_k)=this%F_Grho (i+st_i,j+st_j,k+st_k) -flux( 3)
       this%F_GrhoE(i+st_i,j+st_j,k+st_k)=this%F_GrhoE(i+st_i,j+st_j,k+st_k) -flux( 4)
       this%F_GrhoU(i+st_i,j+st_j,k+st_k)=this%F_GrhoU(i+st_i,j+st_j,k+st_k) -flux( 5)
       this%F_GrhoV(i+st_i,j+st_j,k+st_k)=this%F_GrhoV(i+st_i,j+st_j,k+st_k) -flux( 6)
       this%F_GrhoW(i+st_i,j+st_j,k+st_k)=this%F_GrhoW(i+st_i,j+st_j,k+st_k) -flux( 7)
       this%F_Lrho (i+st_i,j+st_j,k+st_k)=this%F_Lrho (i+st_i,j+st_j,k+st_k) -flux( 8)
       this%F_LrhoE(i+st_i,j+st_j,k+st_k)=this%F_LrhoE(i+st_i,j+st_j,k+st_k) -flux( 9)
       this%F_LrhoU(i+st_i,j+st_j,k+st_k)=this%F_LrhoU(i+st_i,j+st_j,k+st_k) -flux(10)
       this%F_LrhoV(i+st_i,j+st_j,k+st_k)=this%F_LrhoV(i+st_i,j+st_j,k+st_k) -flux(11)
       this%F_LrhoW(i+st_i,j+st_j,k+st_k)=this%F_LrhoW(i+st_i,j+st_j,k+st_k) -flux(12)
       this%F_GP   (i+st_i,j+st_j,k+st_k)=this%F_GP   (i+st_i,j+st_j,k+st_k) -flux(13)
       this%F_LP   (i+st_i,j+st_j,k+st_k)=this%F_LP   (i+st_i,j+st_j,k+st_k) -flux(14)
     end subroutine add_fluxes

     function SubCan_advect_mom_RHS(RHOi,RHOnb,oRHOi,oRHOnb,orhoUi,orhoUnb,rho_fterm) result(F_rhoU)
       implicit none
       real(WP)              :: F_rhoU
       real(WP), intent(in)  :: RHOi,RHOnb,oRHOi,oRHOnb,orhoUi,orhoUnb
       real(WP)              :: rhoU_Obar,rho_fterm
       real(WP)              :: beta_i, beta_nb

       ! variable starting with 'o' means old, or timestep n
       ! variable with no prefix means new, or timestep n+1
       ! variable ending with 'i' means current location/index
       ! variable ending with 'nb' means other, or neighbor cell
       ! in the direction nface

       ! rho_fterm includes velocity, area, timestep, sign, and flux velocity

       ! Time interpolation factors
       beta_i = 1.0_WP/(sqrt(RHOi )+sqrt(oRHOi ))/sqrt(oRHOi )
       beta_nb= 1.0_WP/(sqrt(RHOnb)+sqrt(oRHOnb))/sqrt(oRHOnb)
       ! Interpolate in space and time
       rhoU_Obar = 0.5_WP*(beta_i*orhoUi+beta_nb*orhoUnb)

       ! epsilon is to avoid dividing by zero. If the denominator is zero,
       ! the numerator should be too. Is this good coding practice?

       ! Calculate flux
       F_rhoU = rho_fterm*rhoU_Obar

       return
     end function SubCan_advect_mom_RHS

     subroutine SubCan_mom_RHS_flux(i,j,k,dir)
       character(len=1) :: dir
       integer :: i,j,k,st_i,st_j,st_k
       real(WP), dimension(14):: flux
       real(WP):: grho_fterm, lrho_fterm

       !print*,'begin mom RHS flux'
       select case (trim(dir))
       case('x')
          st_i =-1; st_j = 0; st_k = 0
          !grho_fterm = GrhoXf(i,j,k)
          !lrho_fterm = LrhoXf(i,j,k)
       case('y')
          st_i = 0; st_j =-1; st_k = 0
          !grho_fterm = GrhoYf(i,j,k)
          !lrho_fterm = LrhoYf(i,j,k)
       case('z')
          st_i = 0; st_j = 0; st_k =-1
          !grho_fterm = GrhoZf(i,j,k)
          !lrho_fterm = LrhoZf(i,j,k)
       end select

       !print*,'before flux calculations'
       ! Gas phase
       flux(5) = SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%GrhoUold(i,j,k),this%GrhoUold(i+st_i,j+st_j,k+st_k),&
            grho_fterm)
       flux(6) = SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%GrhoVold(i,j,k),this%GrhoVold(i+st_i,j+st_j,k+st_k),&
            grho_fterm)
       flux(7) = SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%GrhoWold(i,j,k),this%GrhoWold(i+st_i,j+st_j,k+st_k),&
            grho_fterm)
       ! Liquid phase
       flux(10)= SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%LrhoUold(i,j,k),this%LrhoUold(i+st_i,j+st_j,k+st_k),&
            lrho_fterm)
       flux(11)= SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%LrhoVold(i,j,k),this%LrhoVold(i+st_i,j+st_j,k+st_k),&
            lrho_fterm)
       flux(12)= SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%LrhoWold(i,j,k),this%LrhoWold(i+st_i,j+st_j,k+st_k),&
            lrho_fterm)

       ! Store update for right cell from face
       this%F_GrhoU(i,j,k)=this%F_GrhoU(i,j,k)  +flux(5)
       this%F_GrhoV(i,j,k)=this%F_GrhoV(i,j,k)  +flux(6)
       this%F_GrhoW(i,j,k)=this%F_GrhoW(i,j,k)  +flux(7)
       this%F_LrhoU(i,j,k)=this%F_LrhoU(i,j,k)  +flux(10)
       this%F_LrhoV(i,j,k)=this%F_LrhoV(i,j,k)  +flux(11)
       this%F_LrhoW(i,j,k)=this%F_LrhoW(i,j,k)  +flux(12)

       ! Store update for left cell from face
       this%F_GrhoU(i+st_i,j+st_j,k+st_k)=this%F_GrhoU(i+st_i,j+st_j,k+st_k) -flux(5)
       this%F_GrhoV(i+st_i,j+st_j,k+st_k)=this%F_GrhoV(i+st_i,j+st_j,k+st_k) -flux(6)
       this%F_GrhoW(i+st_i,j+st_j,k+st_k)=this%F_GrhoW(i+st_i,j+st_j,k+st_k) -flux(7)
       this%F_LrhoU(i+st_i,j+st_j,k+st_k)=this%F_LrhoU(i+st_i,j+st_j,k+st_k) -flux(10)
       this%F_LrhoV(i+st_i,j+st_j,k+st_k)=this%F_LrhoV(i+st_i,j+st_j,k+st_k) -flux(11)
       this%F_LrhoW(i+st_i,j+st_j,k+st_k)=this%F_LrhoW(i+st_i,j+st_j,k+st_k) -flux(12)

       return
     end subroutine SubCan_mom_RHS_flux

     subroutine SubCan_LHS(i,j,k,dir,Ga_i,Ga_nb,La_i,La_nb)
       character(len=1) :: dir
       real(WP)         :: Ga_i,Ga_nb,La_i,La_nb
       real(WP)         :: RHOi,RHOnb,oRHOi,oRHOnb
       integer          :: i,j,k,st_i,st_j,st_k
       real(WP)         :: zeta_i,zeta_nb
       real(WP)         :: grho_fterm, lrho_fterm

       ! variable starting with 'o' means old
       ! variable ending with 'i' means current location/index
       ! variable ending with 'nb' means other, or neighbor cell
       ! in the direction nface

       ! rho is the density

       ! this routine is only designed for single phase
       ! only one of grho_fterm or lrho_fterm should be nonzero
       ! this is the portion of the momentum flux that uses the new velocity

       select case (trim(dir))
       case('x')
          st_i =-1; st_j = 0; st_k = 0
          !grho_fterm = -GrhoXf(i,j,k)
          !lrho_fterm = -LrhoXf(i,j,k)
       case('y')
          st_i = 0; st_j =-1; st_k = 0
          !grho_fterm = -GrhoYf(i,j,k)
          !lrho_fterm = -LrhoYf(i,j,k)
       case('z')
          st_i = 0; st_j = 0; st_k =-1
          !grho_fterm = -GrhoZf(i,j,k)
          !lrho_fterm = -LrhoZf(i,j,k)
       end select
       ! Negative sign since this is on the left hand side

       ! Density assignment
       RHOi  = this%RHO   (i,j,k); RHOnb  = this%RHO   (i+st_i,j+st_j,k+st_k);
       oRHOi = this%RHOold(i,j,k); oRHOnb = this%RHOold(i+st_i,j+st_j,k+st_k);

       ! Flux calculation
       zeta_i = 1.0_WP/(sqrt(RHOi )*(sqrt(RHOi )+sqrt(oRHOi )))
       zeta_nb= 1.0_WP/(sqrt(RHOnb)*(sqrt(RHOnb)+sqrt(oRHOnb)))
       ! Interpolate in time
       Ga_i = zeta_i;      Ga_nb= zeta_nb
       ! Interpolate in space
       Ga_i = 0.5_WP*Ga_i; Ga_nb= 0.5_WP*Ga_nb
       ! Copy to liquid variables
       La_i = Ga_i;        La_nb= Ga_nb
       ! Multiply density flux term for each phase
       Ga_i = Ga_i *grho_fterm;  Ga_nb= Ga_nb*grho_fterm
       La_i = La_i *lrho_fterm;  La_nb= La_nb*lrho_fterm

       return
     end subroutine SubCan_LHS

     function SubCan_advect_KE_RHS(u_star_i,v_star_i,w_star_i,u_star_nb,v_star_nb,w_star_nb) result(F_rhok)
       implicit none
       real(WP)              :: F_rhok
       real(WP), intent(in)  :: u_star_i,u_star_nb,v_star_i,v_star_nb,w_star_i,w_star_nb
       real(WP)              :: k_star_i,k_star_nb
       real(WP)              :: k_starbar,u_starbar,v_starbar,w_starbar

       ! calculate time-interpolated kinetic energy
       k_star_i = 0.5_WP*(u_star_i **2 + v_star_i **2 + w_star_i **2)
       k_star_nb= 0.5_WP*(u_star_nb**2 + v_star_nb**2 + w_star_nb**2)

       ! spatially interpolate variables
       k_starbar = 0.5_WP*(k_star_i+k_star_nb)
       u_starbar = 0.5_WP*(u_star_i+u_star_nb)
       v_starbar = 0.5_WP*(v_star_i+v_star_nb)
       w_starbar = 0.5_WP*(w_star_i+w_star_nb)

       ! construct advected kinetic energy (k_startilde)
       F_rhok = -k_starbar + u_starbar**2 + v_starbar**2 + w_starbar**2
       ! This is multiplied by the density flux term later

       return
     end function SubCan_advect_KE_RHS

     subroutine SubCan_KE_RHS_flux(i,j,k,dir)
       character(len=1) :: dir
       integer :: i,j,k,st_i,st_j,st_k
       real(WP):: grho_fterm,lrho_fterm,f_gen
       real(WP), dimension(3) :: star_vel_i, star_vel_nb
       real(WP), dimension(14):: flux

       select case (trim(dir))
       case('x')
          st_i =-1; st_j = 0; st_k = 0
          !grho_fterm = GrhoXf(i,j,k)
          !lrho_fterm = LrhoXf(i,j,k)
       case('y')
          st_i = 0; st_j =-1; st_k = 0
          !grho_fterm = GrhoYf(i,j,k)
          !lrho_fterm = LrhoYf(i,j,k)
       case('z')
          st_i = 0; st_j = 0; st_k =-1
          !grho_fterm = GrhoZf(i,j,k)
          !lrho_fterm = LrhoZf(i,j,k)
       end select

       ! get starred velocity first
       call calculate_ustar(i     ,j     ,k     ,star_vel_i )
       call calculate_ustar(i+st_i,j+st_j,k+st_k,star_vel_nb)
       ! Kinetic energy
       f_gen = SubCan_advect_KE_RHS( &
            star_vel_i( 1),star_vel_i( 2),star_vel_i( 3), &
            star_vel_nb(1),star_vel_nb(2),star_vel_nb(3))
       ! Gas phase
       flux(4) = f_gen*grho_fterm
       ! Liquid phase
       flux(9) = f_gen*lrho_fterm

       ! Store update for right cell from face
       this%F_GrhoE(i,j,k)=this%F_GrhoE(i,j,k)  +flux(4)
       this%F_LrhoE(i,j,k)=this%F_LrhoE(i,j,k)  +flux(9)

       ! Store update for left cell from face
       this%F_GrhoE(i+st_i,j+st_j,k+st_k)=this%F_GrhoE(i+st_i,j+st_j,k+st_k) -flux(4)
       this%F_LrhoE(i+st_i,j+st_j,k+st_k)=this%F_LrhoE(i+st_i,j+st_j,k+st_k) -flux(9)

       return
     end subroutine SubCan_KE_RHS_flux

     subroutine calculate_ustar(i,j,k,star_vel)
       integer  :: i,j,k
       real(WP), dimension(3) :: star_vel

       ! Calculate starred velocity in each direction
       star_vel(1) = (sqrt(this%RHO(i,j,k))*this%Ui(i,j,k)+sqrt(this%RHOold(i,j,k))*this%Uiold(i,j,k)) / &
            (sqrt(this%RHO(i,j,k))+sqrt(this%RHOold(i,j,k)))
       star_vel(2) = (sqrt(this%RHO(i,j,k))*this%Vi(i,j,k)+sqrt(this%RHOold(i,j,k))*this%Viold(i,j,k)) / &
            (sqrt(this%RHO(i,j,k))+sqrt(this%RHOold(i,j,k)))
       star_vel(3) = (sqrt(this%RHO(i,j,k))*this%Wi(i,j,k)+sqrt(this%RHOold(i,j,k))*this%Wiold(i,j,k)) / &
            (sqrt(this%RHO(i,j,k))+sqrt(this%RHOold(i,j,k)))

     end subroutine calculate_ustar
       
   end subroutine advection_step
   
   !> Add surface tension jump term using CSF
   subroutine add_surface_tension_jump(this,dt,div,vf,contact_model)
      use messager,  only: die
      use vfs_class, only: vfs
      implicit none
      class(mast), intent(inout) :: this
      real(WP), intent(inout) :: dt     !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      class(vfs), intent(inout) :: vf
      integer, intent(in), optional :: contact_model
      integer :: i,j,k
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
            call die('[mast: add_surface_tension_jump] Unknown contact model!')
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
               div(i,j,k)=div(i,j,k)+dt*(sum(this%divp_x(:,i,j,k)*this%DPjx(i:i+1,j,k)/this%rho_U(i:i+1,j,k))&
               &                        +sum(this%divp_y(:,i,j,k)*this%DPjy(i,j:j+1,k)/this%rho_V(i,j:j+1,k))&
               &                        +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
            end do
         end do
      end do
      
   end subroutine add_surface_tension_jump
   
   
   !> Calculate the pressure gradient based on P
   subroutine get_pgrad(this,P,Pgradx,Pgrady,Pgradz)
      implicit none
      class(mast), intent(inout) :: this
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
      class(mast), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! ! Calculate as far as possible each component
      ! do k=this%cfg%kmino_,this%cfg%kmaxo_
      !    do j=this%cfg%jmino_,this%cfg%jmaxo_
      !       do i=this%cfg%imino_,this%cfg%imaxo_-1
      !          Ui(i,j,k)=sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k))
      !       end do
      !    end do
      ! end do
      ! do k=this%cfg%kmino_,this%cfg%kmaxo_
      !    do j=this%cfg%jmino_,this%cfg%jmaxo_-1
      !       do i=this%cfg%imino_,this%cfg%imaxo_
      !          Vi(i,j,k)=sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k))
      !       end do
      !    end do
      ! end do
      ! do k=this%cfg%kmino_,this%cfg%kmaxo_-1
      !    do j=this%cfg%jmino_,this%cfg%jmaxo_
      !       do i=this%cfg%imino_,this%cfg%imaxo_
      !          Wi(i,j,k)=sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1))
      !       end do
      !    end do
      ! end do
      ! ! Add last layer in each direction
      ! if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) Ui(this%cfg%imaxo,:,:)=this%U(this%cfg%imaxo,:,:)
      ! if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) Vi(:,this%cfg%jmaxo,:)=this%V(:,this%cfg%jmaxo,:)
      ! if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) Wi(:,:,this%cfg%kmaxo)=this%W(:,:,this%cfg%kmaxo)
      ! ! Sync it
      ! call this%cfg%sync(Ui)
      ! call this%cfg%sync(Vi)
      ! call this%cfg%sync(Wi)
   end subroutine interp_vel
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cflc,cfl)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(mast), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFLc_x,my_CFLc_y,my_CFLc_z,my_CFLv_x,my_CFLv_y,my_CFLv_z,my_CFLst
      real(WP) :: my_CFLa_x,my_CFLa_y,my_CFLa_z
      real(WP) :: max_nu
      
      ! Get surface tension CFL first
      my_CFLst=huge(1.0_WP)
      if (this%sigma.gt.0.0_WP) then
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  my_CFLst=min(my_CFLst,sqrt((this%Lrho(i,j,k)+this%Grho(i,j,k))*this%cfg%meshsize(i,j,k)**3.0_WP/(4.0_WP*Pi*this%sigma)))
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(my_CFLst,this%CFLst,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      this%CFLst=dt/this%CFLst
      
      ! Get largest kinematic viscosity
      !max_nu=max(this%visc_l/this%rho_l,this%visc_g/this%rho_g)
      
      ! Set the CFLs to zero
      my_CFLc_x=0.0_WP; my_CFLc_y=0.0_WP; my_CFLc_z=0.0_WP
      my_CFLv_x=0.0_WP; my_CFLv_y=0.0_WP; my_CFLv_z=0.0_WP
      my_CFLa_x=0.0_WP; my_CFLa_y=0.0_WP; my_CFLa_z=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFLc_x=max(my_CFLc_x,abs(this%U(i,j,k))*this%cfg%dxmi(i))
               my_CFLc_y=max(my_CFLc_y,abs(this%V(i,j,k))*this%cfg%dymi(j))
               my_CFLc_z=max(my_CFLc_z,abs(this%W(i,j,k))*this%cfg%dzmi(k))
               my_CFLa_x=max(my_CFLa_x,abs(this%U(i,j,k)+sqrt(this%RHOSS2(i,j,k)/this%RHO(i,j,k)))*this%cfg%dxmi(i))
               my_CFLa_y=max(my_CFLa_y,abs(this%V(i,j,k)+sqrt(this%RHOSS2(i,j,k)/this%RHO(i,j,k)))*this%cfg%dymi(j))
               my_CFLa_z=max(my_CFLa_z,abs(this%W(i,j,k)+sqrt(this%RHOSS2(i,j,k)/this%RHO(i,j,k)))*this%cfg%dzmi(k))
               my_CFLv_x=max(my_CFLv_x,4.0_WP*max_nu*this%cfg%dxi(i)**2)
               my_CFLv_y=max(my_CFLv_y,4.0_WP*max_nu*this%cfg%dyi(j)**2)
               my_CFLv_z=max(my_CFLv_z,4.0_WP*max_nu*this%cfg%dzi(k)**2)
            end do
         end do
      end do
      my_CFLc_x=my_CFLc_x*dt; my_CFLc_y=my_CFLc_y*dt; my_CFLc_z=my_CFLc_z*dt
      my_CFLa_x=my_CFLa_x*dt; my_CFLa_y=my_CFLa_y*dt; my_CFLa_z=my_CFLa_z*dt
      my_CFLv_x=my_CFLv_x*dt; my_CFLv_y=my_CFLv_y*dt; my_CFLv_z=my_CFLv_z*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFLc_x,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_y,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLc_z,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLa_x,this%CFLa_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLa_y,this%CFLa_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLa_z,this%CFLa_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_x,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_y,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_CFLv_z,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Return the maximum convective + surface tension CFL
      cflc=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLst)
      
      ! If asked for, also return the maximum overall CFL
      if (present(CFL)) cfl =max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z, &
           this%CFLst,this%CFLa_x,this%CFLa_y,this%CFLa_z)
      
   end subroutine get_cfl
   
   
   !> Calculate the max of our fields
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mast), intent(inout) :: this
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
               !if (this%cfg%VF(i,j,k).gt.0.0_WP) my_divmax=max(my_divmax,abs(this%div(i,j,k)))
            end do
         end do
      end do
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Umax  ,this%Umax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Vmax  ,this%Vmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Wmax  ,this%Wmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Pmax  ,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      !call MPI_ALLREDUCE(my_divmax,this%divmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_max
   
   
   !> Compute MFR through all bcs
   subroutine get_mfr(this)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mast), intent(inout) :: this
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
   
   !> Prepare viscosity arrays from vfs object
   subroutine get_viscosity(this,vf)
      use vfs_class, only: vfs
      implicit none
      class(mast), intent(inout) :: this
      class(vfs), intent(in) :: vf
      integer :: i,j,k
      real(WP) :: liq_vol,gas_vol,tot_vol
      ! ! Compute harmonically-averaged staggered viscosities using subcell phasic volumes
      ! do k=this%cfg%kmino_+1,this%cfg%kmaxo_
      !    do j=this%cfg%jmino_+1,this%cfg%jmaxo_
      !       do i=this%cfg%imino_+1,this%cfg%imaxo_
      !          ! VISC at [xm,ym,zm] - direct sum in x/y/z
      !          liq_vol=sum(vf%Lvol(:,:,:,i,j,k))
      !          gas_vol=sum(vf%Gvol(:,:,:,i,j,k))
      !          tot_vol=gas_vol+liq_vol
      !          this%visc(i,j,k)=0.0_WP
      !          if (tot_vol.gt.0.0_WP) this%visc(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
      !          ! VISC_xy at [x,y,zm] - direct sum in z, staggered sum in x/y
      !          liq_vol=sum(vf%Lvol(0,0,:,i,j,k))+sum(vf%Lvol(1,0,:,i-1,j,k))+sum(vf%Lvol(0,1,:,i,j-1,k))+sum(vf%Lvol(1,1,:,i-1,j-1,k))
      !          gas_vol=sum(vf%Gvol(0,0,:,i,j,k))+sum(vf%Gvol(1,0,:,i-1,j,k))+sum(vf%Gvol(0,1,:,i,j-1,k))+sum(vf%Gvol(1,1,:,i-1,j-1,k))
      !          tot_vol=gas_vol+liq_vol
      !          this%visc_xy(i,j,k)=0.0_WP
      !          if (tot_vol.gt.0.0_WP) this%visc_xy(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
      !          ! VISC_yz at [xm,y,z] - direct sum in x, staggered sum in y/z
      !          liq_vol=sum(vf%Lvol(:,0,0,i,j,k))+sum(vf%Lvol(:,1,0,i,j-1,k))+sum(vf%Lvol(:,0,1,i,j,k-1))+sum(vf%Lvol(:,1,1,i,j-1,k-1))
      !          gas_vol=sum(vf%Gvol(:,0,0,i,j,k))+sum(vf%Gvol(:,1,0,i,j-1,k))+sum(vf%Gvol(:,0,1,i,j,k-1))+sum(vf%Gvol(:,1,1,i,j-1,k-1))
      !          tot_vol=gas_vol+liq_vol
      !          this%visc_yz(i,j,k)=0.0_WP
      !          if (tot_vol.gt.0.0_WP) this%visc_yz(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
      !          ! VISC_zx at [x,ym,z] - direct sum in y, staggered sum in z/x
      !          liq_vol=sum(vf%Lvol(0,:,0,i,j,k))+sum(vf%Lvol(0,:,1,i,j,k-1))+sum(vf%Lvol(1,:,0,i-1,j,k))+sum(vf%Lvol(1,:,1,i-1,j,k-1))
      !          gas_vol=sum(vf%Gvol(0,:,0,i,j,k))+sum(vf%Gvol(0,:,1,i,j,k-1))+sum(vf%Gvol(1,:,0,i-1,j,k))+sum(vf%Gvol(1,:,1,i-1,j,k-1))
      !          tot_vol=gas_vol+liq_vol
      !          this%visc_zx(i,j,k)=0.0_WP
      !          if (tot_vol.gt.0.0_WP) this%visc_zx(i,j,k)=this%visc_g*this%visc_l/(this%visc_l*gas_vol/tot_vol+this%visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
      !       end do
      !    end do
      ! end do
      ! ! Synchronize boundaries - not really needed...
      ! call this%cfg%sync(this%visc)
      ! call this%cfg%sync(this%visc_xy)
      ! call this%cfg%sync(this%visc_yz)
      ! call this%cfg%sync(this%visc_zx)
   end subroutine get_viscosity
   
   
   ! !> Prepare old density arrays from vfs object
   ! subroutine get_olddensity(this,vf)
   !    use vfs_class, only: vfs
   !    implicit none
   !    class(mast), intent(inout) :: this
   !    class(vfs), intent(in) :: vf
   !    integer :: i,j,k
   !    real(WP) :: liq_vol,gas_vol,tot_vol
   !    ! Calculate rho_U/V/Wold using subcell phasic volumes
   !    do k=this%cfg%kmino_+1,this%cfg%kmaxo_
   !       do j=this%cfg%jmino_+1,this%cfg%jmaxo_
   !          do i=this%cfg%imino_+1,this%cfg%imaxo_
   !             ! U-cell
   !             liq_vol=sum(vf%Lvol(0,:,:,i,j,k))+sum(vf%Lvol(1,:,:,i-1,j,k))
   !             gas_vol=sum(vf%Gvol(0,:,:,i,j,k))+sum(vf%Gvol(1,:,:,i-1,j,k))
   !             tot_vol=gas_vol+liq_vol
   !             this%rho_Uold(i,j,k)=1.0_WP
   !             if (tot_vol.gt.0.0_WP) this%rho_Uold(i,j,k)=(liq_vol*this%rho_l+gas_vol*this%rho_g)/tot_vol
   !             ! V-cell
   !             liq_vol=sum(vf%Lvol(:,0,:,i,j,k))+sum(vf%Lvol(:,1,:,i,j-1,k))
   !             gas_vol=sum(vf%Gvol(:,0,:,i,j,k))+sum(vf%Gvol(:,1,:,i,j-1,k))
   !             tot_vol=gas_vol+liq_vol
   !             this%rho_Vold(i,j,k)=1.0_WP
   !             if (tot_vol.gt.0.0_WP) this%rho_Vold(i,j,k)=(liq_vol*this%rho_l+gas_vol*this%rho_g)/tot_vol
   !             ! W-cell
   !             liq_vol=sum(vf%Lvol(:,:,0,i,j,k))+sum(vf%Lvol(:,:,1,i,j,k-1))
   !             gas_vol=sum(vf%Gvol(:,:,0,i,j,k))+sum(vf%Gvol(:,:,1,i,j,k-1))
   !             tot_vol=gas_vol+liq_vol
   !             this%rho_Wold(i,j,k)=1.0_WP
   !             if (tot_vol.gt.0.0_WP) this%rho_Wold(i,j,k)=(liq_vol*this%rho_l+gas_vol*this%rho_g)/tot_vol
   !          end do
   !       end do
   !    end do
   !    ! Synchronize boundaries - not really needed...
   !    call this%cfg%sync(this%rho_Uold)
   !    call this%cfg%sync(this%rho_Vold)
   !    call this%cfg%sync(this%rho_Wold)
   ! end subroutine get_olddensity
   
   
   ! !> Add gravity source term - assumes that rho_[UVW] have been generated before
   ! subroutine addsrc_gravity(this,resU,resV,resW)
   !    implicit none
   !    class(mast), intent(inout) :: this
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
   !    integer :: i,j,k
   !    do k=this%cfg%kmin_,this%cfg%kmax_
   !       do j=this%cfg%jmin_,this%cfg%jmax_
   !          do i=this%cfg%imin_,this%cfg%imax_
   !             if (this%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+this%rho_U(i,j,k)*this%gravity(1)
   !             if (this%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+this%rho_V(i,j,k)*this%gravity(2)
   !             if (this%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+this%rho_W(i,j,k)*this%gravity(3)
   !          end do
   !       end do
   !    end do
   ! end subroutine addsrc_gravity
   
   
   !> Add a static contact line model
   subroutine add_static_contact(this,vf)
      use mathtools, only: normalize
      use vfs_class, only: vfs
      use irl_fortran_interface, only: calculateNormal,calculateVolume
      implicit none
      class(mast), intent(inout) :: this
      class(vfs),  intent(in) :: vf
      integer :: i,j,k
      real(WP), dimension(3) :: nw,mynorm
      real(WP) :: dd,mysurf
      real(WP) :: cos_contact_angle
      real(WP) :: sin_contact_angle
      real(WP) :: tan_contact_angle
      real(WP), dimension(:,:,:), allocatable :: GFM
      real(WP), parameter :: cfactor=1.0_WP
      
      ! Allocate and zero out binarized VF for GFM-style jump distribution
      allocate(GFM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); GFM=0.0_WP
      
      ! Prepare a GFM-based strategy
      GFM=real(nint(vf%VF),WP)
      
      ! Precalculate cos/sin/tan(contact angle)
      cos_contact_angle=cos(this%contact_angle)
      sin_contact_angle=sin(this%contact_angle)
      tan_contact_angle=tan(this%contact_angle)
      
      ! Loop over domain and identify cells that require contact angle model in GFM style
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! Check if we have an interface on the x-face then check walls
               mysurf=abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
               if (GFM(i,j,k).ne.GFM(i-1,j,k).and.mysurf.gt.0.0_WP) then
                  mynorm=normalize(abs(calculateVolume(vf%interface_polygon(1,i-1,j,k)))*calculateNormal(vf%interface_polygon(1,i-1,j,k))+&
                  &                abs(calculateVolume(vf%interface_polygon(1,i  ,j,k)))*calculateNormal(vf%interface_polygon(1,i  ,j,k)))
                  if (this%umask(i,j-1,k).eq.1) then
                     nw=[0.0_WP,+1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*GFM(i-1:i,j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%umask(i,j+1,k).eq.1) then
                     nw=[0.0_WP,-1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*GFM(i-1:i,j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%umask(i,j,k-1).eq.1) then
                     nw=[0.0_WP,0.0_WP,+1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*GFM(i-1:i,j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%umask(i,j,k+1).eq.1) then
                     nw=[0.0_WP,0.0_WP,-1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     this%Pjx(i,j,k)=this%Pjx(i,j,k)+this%sigma*sum(this%divu_x(:,i,j,k)*GFM(i-1:i,j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
               end if
               
               ! Check if we have an interface on the y-face then check walls
               mysurf=abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
               if (GFM(i,j,k).ne.GFM(i,j-1,k).and.mysurf.gt.0.0_WP) then
                  mynorm=normalize(abs(calculateVolume(vf%interface_polygon(1,i,j-1,k)))*calculateNormal(vf%interface_polygon(1,i,j-1,k))+&
                  &                abs(calculateVolume(vf%interface_polygon(1,i,j  ,k)))*calculateNormal(vf%interface_polygon(1,i,j  ,k)))
                  if (this%vmask(i-1,j,k).eq.1) then
                     nw=[+1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*GFM(i,j-1:j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%vmask(i+1,j,k).eq.1) then
                     nw=[-1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*GFM(i,j-1:j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%vmask(i,j,k-1).eq.1) then
                     nw=[0.0_WP,0.0_WP,+1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*GFM(i,j-1:j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%vmask(i,j,k+1).eq.1) then
                     nw=[0.0_WP,0.0_WP,-1.0_WP]; dd=cfactor*this%cfg%dz(k)
                     this%Pjy(i,j,k)=this%Pjy(i,j,k)+this%sigma*sum(this%divv_y(:,i,j,k)*GFM(i,j-1:j,k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
               end if
               
               ! Check if we have an interface on the z-face then check walls
               mysurf=abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))+abs(calculateVolume(vf%interface_polygon(1,i,j,k)))
               if (GFM(i,j,k).ne.GFM(i,j,k-1).and.mysurf.gt.0.0_WP) then
                  mynorm=normalize(abs(calculateVolume(vf%interface_polygon(1,i,j,k-1)))*calculateNormal(vf%interface_polygon(1,i,j,k-1))+&
                  &                abs(calculateVolume(vf%interface_polygon(1,i,j,k  )))*calculateNormal(vf%interface_polygon(1,i,j,k  )))
                  if (this%wmask(i-1,j,k).eq.1) then
                     nw=[+1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*GFM(i,j,k-1:k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%wmask(i+1,j,k).eq.1) then
                     nw=[-1.0_WP,0.0_WP,0.0_WP]; dd=cfactor*this%cfg%dx(i)
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*GFM(i,j,k-1:k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%wmask(i,j-1,k).eq.1) then
                     nw=[0.0_WP,+1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*GFM(i,j,k-1:k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
                  if (this%wmask(i,j+1,k).eq.1) then
                     nw=[0.0_WP,-1.0_WP,0.0_WP]; dd=cfactor*this%cfg%dy(j)
                     this%Pjz(i,j,k)=this%Pjz(i,j,k)+this%sigma*sum(this%divw_z(:,i,j,k)*GFM(i,j,k-1:k))*(dot_product(mynorm,nw)-cos_contact_angle)/dd
                  end if
               end if
               
            end do
         end do
      end do
      
      ! Deallocate array
      deallocate(GFM)
      
   end subroutine add_static_contact
   
   
   !> Print out info for two-phase incompressible flow solver
   subroutine mast_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(mast), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Two-phase compressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         !write(output_unit,'(" >   liquid density = ",es12.5)') this%rho_l
         !write(output_unit,'(" > liquid viscosity = ",es12.5)') this%visc_l
         !write(output_unit,'(" >      gas density = ",es12.5)') this%rho_g
         !write(output_unit,'(" >    gas viscosity = ",es12.5)') this%visc_g
      end if
      
   end subroutine mast_print
   
   
 end module mast_class
