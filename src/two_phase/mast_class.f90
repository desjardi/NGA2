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

   ! Buffer for labeling BCs
   character(len=str_medium), public :: bc_scope
   
   !> Boundary conditions for the two-phase solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      type(iterator) :: itr                               !< This is the iterator for the bcond - this identifies the (i,j,k)
      character(len=1) :: face                            !< Bcond face (x/y/z)
      integer :: dir                                      !< Bcond direction (+1,-1,0 for interior)
      real(WP) :: rdir                                    !< Bcond direction (real variable)
   end type bcond

   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Two-phase compressible solver object definition ([M]ultiphase [A]ll-Mach [S]emi-Lagrangian [T]ransport)
   type :: mast
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_MAST'    !< Solver name (default=UNNAMED_MAST)

      ! Solver parameters
      real(WP) :: shs_wt                                  !< Shock sensor weight (higher value --> more sensitive)
      
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
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver

      ! Ad-hoc fixing of variables
      logical :: energy_fix                               !< Turns on correction of bad energy values at the end of the timestep

      ! Note: naming convention of "U-cell", etc. is maintained, but this is a collocated solver.
      ! U-cells correspond to face velocities, which are interpolated and projected quantities.
      ! The cell-centered velocities, such as "Ui", correspond to the conserved mass and momentum and are not interpolated.
      
      ! Interpolated density fields
      real(WP), dimension(:,:,:), allocatable :: rho_U    !< Density field array on U-cell
      real(WP), dimension(:,:,:), allocatable :: rho_V    !< Density field array on V-cell
      real(WP), dimension(:,:,:), allocatable :: rho_W    !< Density field array on W-cell
      ! Interpolated pressure fields
      real(WP), dimension(:,:,:), allocatable :: P_U      !< Pressure field array on U-cell
      real(WP), dimension(:,:,:), allocatable :: P_V      !< Pressure field array on V-cell
      real(WP), dimension(:,:,:), allocatable :: P_W      !< Pressure field array on W-cell
      
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
      real(WP), dimension(:,:,:), allocatable :: GrhoE,  LrhoE   !< phase energy arrays
      real(WP), dimension(:,:,:), allocatable :: GP,     LP      !< phase pressure arrays
      real(WP), dimension(:,:,:), allocatable :: GrhoSS2,LrhoSS2 !< phase bulk modulus arrays
      
      ! Old flow variables
      real(WP), dimension(:,:,:), allocatable :: RHOold
      real(WP), dimension(:,:,:), allocatable :: Uiold,   Viold,   Wiold
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

      ! Density flux arrays
      real(WP), dimension(:,:,:,:), allocatable :: GrhoFf                !< Gas density flux (used for SC advection)
      real(WP), dimension(:,:,:,:), allocatable :: LrhoFf                !< Liquid density flux

      ! Scalar gradient arrays (for reconstruction)
      real(WP), dimension(:,:,:,:), allocatable :: gradGrho,gradGrhoE,gradGIE,gradGP !< Gradients: Gas scalars
      real(WP), dimension(:,:,:,:), allocatable :: gradGrhoU,gradGrhoV,gradGrhoW     !< Gradients: Gas momentum
      real(WP), dimension(:,:,:,:), allocatable :: gradLrho,gradLrhoE,gradLIE,gradLP !< Gradients: Liquid scalars
      real(WP), dimension(:,:,:,:), allocatable :: gradLrhoU,gradLrhoV,gradLrhoW     !< Gradients: Liquid momentum

      ! Hybrid advection
      integer,  dimension(:,:,:,:), allocatable :: sl_face ! < Flag for flux method switching
      ! Terms needed for pressure relaxation
      real(WP), dimension(:,:,:), allocatable :: srcVF   ! < Predicted volume exchange during advection
      real(WP), dimension(:,:,:), allocatable :: Hpjump  ! < Helmholtz pressure jump
      real(WP), dimension(:,:,:), allocatable :: dHpjump ! < Helmholtz pressure jump increment, also used for old Hpjump

      ! Temporary arrays for a few things
      real(WP), dimension(:,:,:), pointer :: tmp1,tmp2,tmp3
      
      ! Pressure solver
      type(ils) :: psolv                                  !< Iterative linear solver object for the pressure Helmholtz equation
      
      ! Implicit momentum solver
      type(ils) :: implicit                               !< Iterative linear solver object for an implicit prediction of the advection residual
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itpi_x,itpi_y,itpi_z   !< Interpolation fom cell center to face (for scalars, e.g. pressure)
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
      procedure :: apply_bcond                            !< Apply boundary conditions as specified
      procedure :: add_static_contact                     !< Add static contact line model to surface tension jump

      ! For initialization of simulation
      procedure :: setup                                  !< Finish configuring the flow solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: init_phase_bulkmod                     !< Calculate phase bulk moduli according to other variables
      ! For beginning of timestep (before iterative loop)
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: flag_sl                                !< Flag where SL scheme needs to be used
      procedure :: flow_reconstruct                       !< Calculate in-cell gradients for reconstructing flow variables
      procedure :: reinit_phase_pressure                  !< Calculates phase pressures according to other conserved variables
      
      ! For advection solve
      procedure :: advection_step                         !< Full, hybrid advection step
      ! For convenience in advection routines
      procedure :: GKEold, LKEold                         !< Calculate kinetic energy at timestep 'n'
      ! For viscous/dissipative/body forces (will address later)
      ! For setting up pressure solve
      procedure :: pressureproj_prepare                   !< Overall subroutine for pressure solve setup
      procedure :: interp_pressure_density                !< Calculate pressure and density interpolated to face
      procedure :: interp_vel_basic                       !< Calculate interpolated velocity
      procedure :: interp_vel_full                        !< Calculate interpolated velocity, with additional terms
      procedure :: terms_modified_face_velocity           !< Calculate additional terms for velocity interpolation
      procedure :: update_Helmholtz_LHS                   !< Calculate the left-hand side of the pressure equation
      procedure :: update_Helmholtz_RHS                   !< Calculate the right-hand side of the pressure equation
      procedure :: update_surface_tension_jump            !< Calculate surface tension jump
      procedure :: harmonize_advpressure_bulkmod          !< Calculate PA and RHOSS2 for Helmholtz equation
      ! For pressure correction
      procedure :: pressureproj_correct                   !< Overall subroutine for pressure solve correction
      ! For pressure relaxation
      procedure :: pressure_relax                         !< Loop to relax pressure, change VF and phase values, at end of timetep
      procedure :: pressure_relax_one                     !< Routine to perform mechanical relaxation for an individual cell

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
   function constructor(cfg,name,vf) result(self)
      use vfs_class, only: vfs
      implicit none
      type(mast) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      class(vfs),  intent(inout) :: vf     !< The volume fraction solver
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))

      ! Set fixes off initially
      self%energy_fix = .false.
      
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
      allocate(self%P_U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P_U=0.0_WP
      allocate(self%P_V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P_V=0.0_WP
      allocate(self%P_W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P_W=0.0_WP
      allocate(self%visc   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc   =0.0_WP
      allocate(self%visc_xy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_xy=0.0_WP
      allocate(self%visc_yz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_yz=0.0_WP
      allocate(self%visc_zx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_zx=0.0_WP
      ! Two-phase
      allocate(self%Grho (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Grho =0.0_WP
      allocate(self%Lrho (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Lrho =0.0_WP
      allocate(self%GrhoE(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoE=0.0_WP
      allocate(self%LrhoE(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoE=0.0_WP
      allocate(self%GP   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GP   =0.0_WP
      allocate(self%LP   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LP   =0.0_WP
      allocate(self%GrhoSS2(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoSS2=0.0_WP
      allocate(self%LrhoSS2(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoSS2=0.0_WP

      
      ! Allocate old flow variables
      allocate(self%RHOold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%RHOold=0.0_WP
      allocate(self%Uiold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Uiold=0.0_WP
      allocate(self%Viold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Viold=0.0_WP
      allocate(self%Wiold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Wiold=0.0_WP
      allocate(self%GrhoEold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GrhoEold=0.0_WP
      allocate(self%LrhoEold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LrhoEold=0.0_WP
      allocate(self%GPold   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%GPold   =0.0_WP
      allocate(self%LPold   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%LPold   =0.0_WP
      
      ! Flux sum arrays need to be preallocated
      allocate(self%F_VF (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_VF=0.0_WP
      allocate(self%F_VOL(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%F_VOL=0.0_WP
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

      ! Density flux arrays
      allocate(self%GrhoFf (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3)); self%GrhoFf =0.0_WP
      allocate(self%LrhoFf (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3)); self%LrhoFf =0.0_WP

      ! Gradients
      allocate(self%gradGrho (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGrho =0.0_WP
      allocate(self%gradGrhoE(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGrhoE=0.0_WP
      allocate(self%gradGIE  (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGIE  =0.0_WP
      allocate(self%gradGP   (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGP   =0.0_WP
      allocate(self%gradGrhoU(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGrhoU=0.0_WP
      allocate(self%gradGrhoV(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGrhoV=0.0_WP
      allocate(self%gradGrhoW(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradGrhoW=0.0_WP
      allocate(self%gradLrho (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLrho =0.0_WP
      allocate(self%gradLrhoE(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLrhoE=0.0_WP
      allocate(self%gradLIE  (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLIE  =0.0_WP
      allocate(self%gradLP   (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLP   =0.0_WP
      allocate(self%gradLrhoU(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLrhoU=0.0_WP
      allocate(self%gradLrhoV(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLrhoV=0.0_WP
      allocate(self%gradLrhoW(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%gradLrhoW=0.0_WP

      ! Hybrid advection
      allocate(self%sl_face(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3)); self%sl_face=0
      ! Pressure relaxation
      allocate(self%srcVF(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%srcVF=0.0_WP
      allocate(self%Hpjump(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Hpjump=0.0_WP
      allocate(self%dHpjump(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dHpjump=0.0_WP

      ! Arrays employed temporarily in a few places
      allocate(self%tmp1(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%tmp1   =0.0_WP
      allocate(self%tmp2(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%tmp2   =0.0_WP
      allocate(self%tmp3(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%tmp3   =0.0_WP
      
      ! Allocate vfs objects that are required for the MAST solver
      call vf%allocate_supplement()
      
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

      ! Adjust interpolation coefficients to cell faces - used for scalars
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Linear interpolation in x
               if (this%mask(i,j,k).eq.0.and.this%mask(i-1,j,k).gt.0) this%itpi_x(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i-1,j,k).eq.0) this%itpi_x(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in y               
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j-1,k).gt.0) this%itpi_y(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i,j-1,k).eq.0) this%itpi_y(:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in z
               if (this%mask(i,j,k).eq.0.and.this%mask(i,j,k-1).gt.0) this%itpi_z(:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(i,j,k).gt.0.and.this%mask(i,j,k-1).eq.0) this%itpi_z(:,i,j,k)=[1.0_WP,0.0_WP]
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
   subroutine add_bcond(this,name,type,locator,face,dir)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(mast), intent(inout) :: this
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
      character(len=1), intent(in) :: face
      integer, intent(in) :: dir
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
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1

      ! For this collocated solver, BCs specify cell-centered values and behavior,
      ! and the behavior of the face-centered values depends on the type of BC.
      ! (For now, dirichlet BCs always enforce the flowrate, i.e., mask face velocities,
      !  and clipped neumann BCs modify the sl_face flags to what works best at the outlet.)

      ! To work properly with n_ and locators, dirichlets work by specifying the face that is
      ! a part of the dirichlet, and the cell behind it (determined by the direction) is included.
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         do n=1,new_bc%itr%n_
            i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
            select case(face)
            case('x')
               ! Mask face
               this%umask(i,j,k)=2
               ! Mask cell
               this%mask(i+min(0,dir),j,k)=2
            case('y')
               ! Mask face
               this%vmask(i,j,k)=2
               ! Mask cell
               this%mask(i,j+min(0,dir),k)=2
            case('z')
               ! Mask face
               this%wmask(i,j,k)=2
               ! Mask cell
               this%mask(i,j,k+min(0,dir))=2
            end select
         end do
         
      case (neumann) !< Neumann has to be at existing wall or at domain boundary!
      case (clipped_neumann)
      !case (convective)
      case default
         call die('[mast add_bcond] Unknown bcond type')
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
   subroutine apply_bcond(this,dt,scope)
      use messager, only: die
      implicit none
      class(mast), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: i,j,k,n,stag,ii
      type(bcond), pointer :: my_bc
      character(len=str_medium) :: scope
      
      ! Scope
      ! Available options: density, phase_momentum, momentum, energy, velocity, flag
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
            case (dirichlet)               !< Apply Dirichlet conditions
               
               ! This is done by the user directly for conseerved quantities and face velocity

               ! Condition is necessary for sl_face
               if (trim(adjustl(scope)).eq.'flag') then
                  ! Usage of SL scheme is same for dirichlet and neumann
                  do n=1,my_bc%itr%n_
                     i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                     select case (my_bc%face)
                     case ('x')
                        do ii = 0,abs(2*my_bc%dir)
                           this%sl_face(i-ii*my_bc%dir,j,k,1) = 1
                        end do
                     case ('y')
                        do ii = 0,abs(2*my_bc%dir)
                           this%sl_face(i,j+ii*my_bc%dir,k,2) = 1
                        end do
                     case ('z')
                        do ii = 0,abs(2*my_bc%dir)
                           this%sl_face(i,j,k-ii*my_bc%dir,3) = 1
                        end do
                     end select
                  end do
               end if
               
            case (neumann,clipped_neumann) !< Apply Neumann condition to all 3 components
               ! Handle index shift due to staggering
               stag=min(my_bc%dir,0)
               ! Implement based on bcond direction
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  select case (trim(adjustl(scope)))
                  case('density')
                     this%Grho(i,j,k)=this%Grho(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                     this%Lrho(i,j,k)=this%Lrho(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  case('momentum')
                     this%rhoUi(i,j,k)=this%rhoUi(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                     this%rhoVi(i,j,k)=this%rhoVi(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                     this%rhoWi(i,j,k)=this%rhoWi(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  case('energy')
                     this%GrhoE(i,j,k)=this%GrhoE(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                     this%LrhoE(i,j,k)=this%LrhoE(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  case('velocity')
                     select case (my_bc%face)
                     case ('x')
                        this%U(i     ,j    ,k    )=this%U(i-my_bc%dir     ,j    ,k    )
                        this%V(i+stag,j:j+1,k    )=this%V(i-my_bc%dir+stag,j:j+1,k    )
                        this%W(i+stag,j    ,k:k+1)=this%W(i-my_bc%dir+stag,j    ,k:k+1)
                     case ('y')
                        this%U(i:i+1,j+stag,k    )=this%U(i:i+1,j-my_bc%dir+stag,k    )
                        this%V(i    ,j     ,k    )=this%V(i    ,j-my_bc%dir     ,k    )
                        this%W(i    ,j+stag,k:k+1)=this%W(i    ,j-my_bc%dir+stag,k:k+1)
                     case ('z')
                        this%U(i:i+1,j    ,k+stag)=this%U(i:i+1,j    ,k-my_bc%dir+stag)
                        this%V(i    ,j:j+1,k+stag)=this%V(i    ,j:j+1,k-my_bc%dir+stag)
                        this%W(i    ,j    ,k     )=this%W(i    ,j    ,k-my_bc%dir     )
                     end select
                  case('flag')
                     ! Extend 
                     select case (my_bc%face)
                     case ('x')
                        do ii = 0,abs(2*my_bc%dir)
                           this%sl_face(i-ii*my_bc%dir ,j     ,k     ,1     ) = 1
                        end do
                     case ('y')
                        do ii = 0,abs(2*my_bc%dir)
                           this%sl_face(i     ,j-ii*my_bc%dir ,k     ,2     ) = 1
                        end do
                     case ('z')
                        do ii = 0,abs(2*my_bc%dir)
                           this%sl_face(i     ,j     ,k-ii*my_bc%dir ,3     ) = 1
                        end do
                     end select
                  end select
               end do
               ! If needed, clip
               if (my_bc%type.eq.clipped_neumann.and.trim(adjustl(scope)).eq.'velocity') then
                  select case (my_bc%face)
                  case ('x')
                     if (this%U(i,j,k)*my_bc%rdir.lt.0.0_WP) this%U(i,j,k)=0.0_WP
                  case ('y')
                     if (this%V(i,j,k)*my_bc%rdir.lt.0.0_WP) this%V(i,j,k)=0.0_WP
                  case ('z')
                     if (this%W(i,j,k)*my_bc%rdir.lt.0.0_WP) this%W(i,j,k)=0.0_WP
                  end select
               end if
               
            !case (convective)   ! Not implemented yet!
               
            case default
               call die('[mast apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond (Needs to be done anyways)
         select case(trim(adjustl(scope)))
         case('density')
            call this%cfg%sync(this%Grho)
            call this%cfg%sync(this%Lrho)
         case('momentum')
            call this%cfg%sync(this%rhoUi)
            call this%cfg%sync(this%rhoVi)
            call this%cfg%sync(this%rhoWi)
         case('energy')
            call this%cfg%sync(this%GrhoE)
            call this%cfg%sync(this%LrhoE)
         case('velocity')
            call this%cfg%sync(this%U)
            call this%cfg%sync(this%V)
            call this%cfg%sync(this%W)
         case('flag')
            call this%cfg%sync(this%sl_face(:,:,:,1))
            call this%cfg%sync(this%sl_face(:,:,:,2))
            call this%cfg%sync(this%sl_face(:,:,:,3))
         end select
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond

   ! Function to more easily calculate gas kinetic energy
   function GKEold(this,i,j,k) result(val)
     implicit none
     class(mast), intent(inout) :: this
     real(WP) :: val
     integer, intent(in) :: i,j,k
     val = 0.5_WP*this%Grhoold(i,j,k)*(this%Uiold(i,j,k)**2+this%Viold(i,j,k)**2+this%Wiold(i,j,k)**2)
     return
   end function GKEold

   ! Function to more easily calculate liquid kinetic energy
   function LKEold(this,i,j,k) result(val)
     implicit none
     class(mast), intent(inout) :: this
     real(WP) :: val
     integer, intent(in) :: i,j,k
     val = 0.5_WP*this%Lrhoold(i,j,k)*(this%Uiold(i,j,k)**2+this%Viold(i,j,k)**2+this%Wiold(i,j,k)**2)
     return
   end function LKEold

   subroutine flag_sl(this,dt,vf)
     use vfs_class, only: vfs, VFhi, VFlo
     implicit none

     class(mast), intent(inout) :: this   !< The two-phase all-Mach flow solver
     class(vfs),  intent(inout) :: vf     !< The volume fraction solver
     real(WP), intent(in) :: dt           !< Passed to BC routine
     integer :: VF_check,shock_check,n_band,ni,nj,nk
     integer :: i,j,k

     n_band = 1 !ceiling(max_CFL)

     this%sl_face = 0
     ! 1 is for semi-Lagrangian, 0 is for centered scheme

     ! Loop over the domain and determine multiphase locations within band and shock locations
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              if (vf%VFold(i,j,k).lt.VFlo) then
                 ! In fully gas cells, check in band for change in phase
                 VF_check = 0
                 do ni=-n_band,n_band
                    do nj=-n_band,n_band
                       do nk=-n_band,n_band
                          if (vf%VFold(max(this%cfg%imino_,min(i+ni,this%cfg%imaxo_)), &
                               max(this%cfg%jmino_,min(j+nj,this%cfg%jmaxo_)), &
                               max(this%cfg%kmino_,min(k+nk,this%cfg%kmaxo_))).gt.VFlo .and. &
                               this%mask(max(this%cfg%imino_,min(i+ni,this%cfg%imaxo_)), &
                               max(this%cfg%jmino_,min(j+nj,this%cfg%jmaxo_)), &
                               max(this%cfg%kmino_,min(k+nk,this%cfg%kmaxo_))).ne.1) then
                             VF_check = 1
                          end if
                       end do
                    end do
                 end do
                 if (VF_check.eq.0) then
                    ! If change of phase is not found in band, check shock sensor
                    if (shock_sensor(i,j,k).gt.0) then
                       ! If shock is indicated, set values to negative, band will be extended around these cells
                       this%sl_face(i,j,k,1:3)=-1
                       this%sl_face(i+1,j,k,1)=-1; this%sl_face(i,j+1,k,2)=-1; this%sl_face(i,j,k+1,3)=-1;
                    end if
                 else
                    ! If change of phase is found in band, turn on SL fluxing
                    this%sl_face(i,j,k,1:3)=1
                    this%sl_face(i+1,j,k,1)=1; this%sl_face(i,j+1,k,2)=1; this%sl_face(i,j,k+1,3)=1;
                 end if
              elseif (vf%VFold(i,j,k).gt.VFhi) then
                 ! In fully liquid cells, check in band for change in phase
                 VF_check = 0
                 do ni=-n_band,n_band
                    do nj=-n_band,n_band
                       do nk=-n_band,n_band
                          if (vf%VFold(max(this%cfg%imino_,min(i+ni,this%cfg%imaxo_)), &
                               max(this%cfg%jmino_,min(j+nj,this%cfg%jmaxo_)), &
                               max(this%cfg%kmino_,min(k+nk,this%cfg%kmaxo_))).lt.VFhi .and. &
                               this%mask(max(this%cfg%imino_,min(i+ni,this%cfg%imaxo_)), &
                               max(this%cfg%jmino_,min(j+nj,this%cfg%jmaxo_)), &
                               max(this%cfg%kmino_,min(k+nk,this%cfg%kmaxo_))).ne.1) then
                             VF_check = 1
                          end if
                       end do
                    end do
                 end do
                 if (VF_check.eq.0) then
                    ! If change of phase is not found in band, check shock sensor
                    if (shock_sensor(i,j,k).gt.0.0_WP) then
                       ! If shock is indicated, set values to negative, band will be extended around these cells
                       this%sl_face(i,j,k,1:3)=-1
                       this%sl_face(i+1,j,k,1)=-1; this%sl_face(i,j+1,k,2)=-1; this%sl_face(i,j,k+1,3)=-1;
                    end if
                 else
                    ! If change of phase is found in band, turn on SL fluxing
                    this%sl_face(i,j,k,1:3)=1
                    this%sl_face(i+1,j,k,1)=1; this%sl_face(i,j+1,k,2)=1; this%sl_face(i,j,k+1,3)=1;
                 end if
              else
                 ! If cell is multiphase, turn on SL fluxing
                 this%sl_face(i,j,k,1:3)=1
                 this%sl_face(i+1,j,k,1)=1; this%sl_face(i,j+1,k,2)=1; this%sl_face(i,j,k+1,3)=1;
              end if
           end do
        end do
     end do

     ! BC for flag
     bc_scope = 'flag'
     call this%apply_bcond(dt,bc_scope)

     ! Loop over the domain and check within band containing shock locations
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              if (this%sl_face(i,j,k,1).eq.0.or.this%sl_face(i,j,k,2).eq.0.or.this%sl_face(i,j,k,3).eq.0) then
                 shock_check = 0
                 do ni=-n_band,n_band
                    do nj=-n_band,n_band
                       do nk=-n_band,n_band
                          if (minval(this%sl_face(max(this%cfg%imino_,min(i+ni,this%cfg%imaxo_)), &
                               max(this%cfg%jmino_,min(j+nj,this%cfg%jmaxo_)), &
                               max(this%cfg%kmino_,min(k+nk,this%cfg%kmaxo_)),:)).eq.-1) then
                             shock_check = 1
                          end if
                       end do
                    end do
                 end do
                 if (shock_check.eq.1) then
                    ! If shock is found in band, turn on SL fluxing, but only in faces where shock is not indicated
                    ! I want 0 to go to 1, 1 to stay at 1, and -1 to stay at -1
                    this%sl_face(i  ,j,k,1)=min(1,2*this%sl_face(i  ,j,k,1)+1)
                    this%sl_face(i+1,j,k,1)=min(1,2*this%sl_face(i+1,j,k,1)+1)
                    this%sl_face(i,j  ,k,2)=min(1,2*this%sl_face(i,j  ,k,2)+1)
                    this%sl_face(i,j+1,k,2)=min(1,2*this%sl_face(i,j+1,k,2)+1)
                    this%sl_face(i,j,k  ,3)=min(1,2*this%sl_face(i,j,k  ,3)+1)
                    this%sl_face(i,j,k+1,3)=min(1,2*this%sl_face(i,j,k+1,3)+1)
                 end if
              end if
           end do
        end do
     end do

     ! Make all the negatives positive
     this%sl_face = abs(this%sl_face)

     ! BCs again
     bc_scope = 'flag'
     call this%apply_bcond(dt,bc_scope)

   contains

     function shock_sensor(i,j,k) result(cb)
       real(WP) :: cb
       integer :: i,j,k

       ! --- Bhagatwala and Lele, modified --- !
       ! wt weights how important the acoustic scale is
       cb=abs(sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)) &
            + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)) &
            + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)) &
            )-min(this%cfg%dxi(i),this%cfg%dyi(j),this%cfg%dzi(k)) &
            / sqrt(this%RHOSS2(i,j,k)/this%RHO(i,j,k)) / this%shs_wt

     end function shock_sensor

   end subroutine flag_sl


   ! ==================================================== !
   ! Form conservative linear reconstruction in each cell !
   ! ==================================================== !
   subroutine flow_reconstruct(this,vf)
     use vfs_class, only: vfs
     implicit none
     class(mast), intent(inout) :: this   !< The two-phase all-Mach flow solver
     class(vfs),  intent(inout) :: vf     !< The volume fraction solver
     integer  :: i,j,k,ip,jp,kp,im,jm,km,idir
     real(WP) :: idp,idm !,idpg,idpl,idmg,idml

     ! Zero gradients
     this%gradGrho =0.0_WP
     this%gradGrhoE=0.0_WP
     this%gradGrhoU=0.0_WP
     this%gradGrhoV=0.0_WP
     this%gradGrhoW=0.0_WP
     this%gradLrho =0.0_WP
     this%gradLrhoE=0.0_WP
     this%gradLrhoU=0.0_WP
     this%gradLrhoV=0.0_WP
     this%gradLrhoW=0.0_WP

     this%gradGIE  =0.0_WP
     this%gradLIE  =0.0_WP

     this%gradGP   =0.0_WP
     this%gradLP   =0.0_WP

     ! Gradients for linear reconstruction
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! No need to calculate gradient inside of wall cell
              if (this%mask(i,j,k).eq.1) cycle
              ! Cycle through directions
              do idir=1,3
                 select case (idir)
                 case (1) ! X gradient
                    ip=i+1; jp=j; kp=k; idp=this%cfg%dxmi(i+1)
                    im=i-1; jm=j; km=k; idm=this%cfg%dxmi(i  )
                    ! idpg=1.0_WP/(vf%Gbary(1,i+1,j,k)-vf%Gbary(1,i,j,k))
                    ! idpl=1.0_WP/(vf%Lbary(1,i+1,j,k)-vf%Lbary(1,i,j,k))
                    ! idmg=1.0_WP/(vf%Gbary(1,i,j,k)-vf%Gbary(1,i-1,j,k))
                    ! idml=1.0_WP/(vf%Lbary(1,i,j,k)-vf%Lbary(1,i-1,j,k))
                 case (2) ! Y gradient
                    ip=i; jp=j+1; kp=k; idp=this%cfg%dymi(j+1)
                    im=i; jm=j-1; km=k; idm=this%cfg%dymi(j  )
                    ! idpg=1.0_WP/(vf%Gbary(2,i,j+1,k)-vf%Gbary(2,i,j,k))
                    ! idpl=1.0_WP/(vf%Lbary(2,i,j+1,k)-vf%Lbary(2,i,j,k))
                    ! idmg=1.0_WP/(vf%Gbary(2,i,j,k)-vf%Gbary(2,i,j-1,k))
                    ! idml=1.0_WP/(vf%Lbary(2,i,j,k)-vf%Lbary(2,i,j-1,k))
                 case (3) ! Z gradient
                    ip=i; jp=j; kp=k+1; idp=this%cfg%dzmi(k+1)
                    im=i; jm=j; km=k-1; idm=this%cfg%dzmi(k  )
                    ! idpg=1.0_WP/(vf%Gbary(3,i,j,k+1)-vf%Gbary(3,i,j,k))
                    ! idpl=1.0_WP/(vf%Lbary(3,i,j,k+1)-vf%Lbary(3,i,j,k))
                    ! idmg=1.0_WP/(vf%Gbary(3,i,j,k)-vf%Gbary(3,i,j,k+1))
                    ! idml=1.0_WP/(vf%Lbary(3,i,j,k)-vf%Lbary(3,i,j,k-1))
                 end select
                 this%gradGrho (idir,i,j,k)=mmgrad((this%Grho (ip,jp,kp)-this%Grho (i,j,k))*idp,&
                      (this%Grho (i,j,k)-this%Grho (im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradGrhoE(idir,i,j,k)=mmgrad((this%GrhoE(ip,jp,kp)-this%GrhoE(i,j,k))*idp,&
                      (this%GrhoE(i,j,k)-this%GrhoE(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradGrhoU(idir,i,j,k)=mmgrad((this%Grho(ip,jp,kp)*this%Ui(ip,jp,kp)-this%Grho(i,j,k)*this%Ui(i,j,k))*idp,&
                      (this%Grho(i,j,k)*this%Ui(i,j,k)-this%Grho(im,jm,km)*this%Ui(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradGrhoV(idir,i,j,k)=mmgrad((this%Grho(ip,jp,kp)*this%Vi(ip,jp,kp)-this%Grho(i,j,k)*this%Vi(i,j,k))*idp,&
                      (this%Grho(i,j,k)*this%Vi(i,j,k)-this%Grho(im,jm,km)*this%Vi(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradGrhoW(idir,i,j,k)=mmgrad((this%Grho(ip,jp,kp)*this%Wi(ip,jp,kp)-this%Grho(i,j,k)*this%Wi(i,j,k))*idp,&
                      (this%Grho(i,j,k)*this%Wi(i,j,k)-this%Grho(im,jm,km)*this%Wi(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)

                 this%gradLrho (idir,i,j,k)=mmgrad((this%Lrho (ip,jp,kp)-this%Lrho (i,j,k))*idp,&
                      (this%Lrho (i,j,k)-this%Lrho (im,jm,km))*idm)*real(floor(       vf%VF(i,j,k)),WP)
                 this%gradLrhoE(idir,i,j,k)=mmgrad((this%LrhoE(ip,jp,kp)-this%LrhoE(i,j,k))*idp,&
                      (this%LrhoE(i,j,k)-this%LrhoE(im,jm,km))*idm)*real(floor(       vf%VF(i,j,k)),WP)
                 this%gradLrhoU(idir,i,j,k)=mmgrad((this%Lrho(ip,jp,kp)*this%Ui(ip,jp,kp)-this%Lrho(i,j,k)*this%Ui(i,j,k))*idp,&
                      (this%Lrho(i,j,k)*this%Ui(i,j,k)-this%Lrho(im,jm,km)*this%Ui(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradLrhoV(idir,i,j,k)=mmgrad((this%Lrho(ip,jp,kp)*this%Vi(ip,jp,kp)-this%Lrho(i,j,k)*this%Vi(i,j,k))*idp,&
                      (this%Lrho(i,j,k)*this%Vi(i,j,k)-this%Lrho(im,jm,km)*this%Vi(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradLrhoW(idir,i,j,k)=mmgrad((this%Lrho(ip,jp,kp)*this%Wi(ip,jp,kp)-this%Lrho(i,j,k)*this%Wi(i,j,k))*idp,&
                      (this%Lrho(i,j,k)*this%Wi(i,j,k)-this%Lrho(im,jm,km)*this%Wi(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)                 

                 this%gradGIE  (idir,i,j,k)=mmgrad((this%GrhoE(ip,jp,kp)-this%GKEold(ip,jp,kp)-this%GrhoE(i,j,k) &
                      +this%GKEold(i,j,k))*idp,(this%GrhoE(i,j,k)-this%GKEold(i,j,k)-this%GrhoE(im,jm,km) &
                      +this%GKEold(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradLIE  (idir,i,j,k)=mmgrad((this%LrhoE(ip,jp,kp)-this%LKEold(ip,jp,kp)-this%LrhoE(i,j,k) &
                      +this%LKEold(i,j,k))*idp,(this%LrhoE(i,j,k)-this%LKEold(i,j,k)-this%LrhoE(im,jm,km) &
                      +this%LKEold(im,jm,km))*idm)*real(floor(       vf%VF(i,j,k)),WP)

                 this%gradGP   (idir,i,j,k)=mmgrad((this%GP(ip,jp,kp)-this%GP(i,j,k))*idp,&
                      (this%GP(i,j,k)-this%GP(im,jm,km))*idm)*real(floor(1.0_WP-vf%VF(i,j,k)),WP)
                 this%gradLP   (idir,i,j,k)=mmgrad((this%LP(ip,jp,kp)-this%LP(i,j,k))*idp,&
                      (this%LP(i,j,k)-this%LP(im,jm,km))*idm)*real(floor(       vf%VF(i,j,k)),WP)

                 ! this%gradGrho (idir,i,j,k)=mmgrad((Grho (ip,jp,kp)-Grho (i,j,k))*idpg,(Grho (i,j,k)-Grho (im,jm,km))*idmg)
                 ! this%gradGrhoE(idir,i,j,k)=mmgrad((GrhoE(ip,jp,kp)-GrhoE(i,j,k))*idpg,(GrhoE(i,j,k)-GrhoE(im,jm,km))*idmg)
                 ! this%gradGrhoU(idir,i,j,k)=mmgrad((GrhoU(ip,jp,kp)-GrhoU(i,j,k))*idpg,(GrhoU(i,j,k)-GrhoU(im,jm,km))*idmg)
                 ! this%gradGrhoV(idir,i,j,k)=mmgrad((GrhoV(ip,jp,kp)-GrhoV(i,j,k))*idpg,(GrhoV(i,j,k)-GrhoV(im,jm,km))*idmg)
                 ! this%gradGrhoW(idir,i,j,k)=mmgrad((GrhoW(ip,jp,kp)-GrhoW(i,j,k))*idpg,(GrhoW(i,j,k)-GrhoW(im,jm,km))*idmg)
                 ! this%gradLrho (idir,i,j,k)=mmgrad((Lrho (ip,jp,kp)-Lrho (i,j,k))*idpl,(Lrho (i,j,k)-Lrho (im,jm,km))*idml)
                 ! this%gradLrhoE(idir,i,j,k)=mmgrad((LrhoE(ip,jp,kp)-LrhoE(i,j,k))*idpl,(LrhoE(i,j,k)-LrhoE(im,jm,km))*idml)
                 ! this%gradLrhoU(idir,i,j,k)=mmgrad((LrhoU(ip,jp,kp)-LrhoU(i,j,k))*idpl,(LrhoU(i,j,k)-LrhoU(im,jm,km))*idml)
                 ! this%gradLrhoV(idir,i,j,k)=mmgrad((LrhoV(ip,jp,kp)-LrhoV(i,j,k))*idpl,(LrhoV(i,j,k)-LrhoV(im,jm,km))*idml)
                 ! this%gradLrhoW(idir,i,j,k)=mmgrad((LrhoW(ip,jp,kp)-LrhoW(i,j,k))*idpl,(LrhoW(i,j,k)-LrhoW(im,jm,km))*idml)

                 ! this%gradGIE(idir,i,j,k)=mmgrad((GrhoE(ip,jp,kp)-oldGKE(ip,jp,kp)-GrhoE(i,j,k)+oldGKE(i,j,k))*idpg,&
                 !                           (GrhoE(i,j,k)-oldGKE(i,j,k)-GrhoE(im,jm,km)+oldGKE(im,jm,km))*idmg)
                 ! this%gradLIE(idir,i,j,k)=mmgrad((LrhoE(ip,jp,kp)-oldLKE(ip,jp,kp)-LrhoE(i,j,k)+oldLKE(i,j,k))*idpl,&
                 !                           (LrhoE(i,j,k)-oldLKE(i,j,k)-LrhoE(im,jm,km)+oldLKE(im,jm,km))*idml)

                 ! this%gradGP(idir,i,j,k)=mmgrad((GP(ip,jp,kp)-GP(i,j,k))*idpg,(GP(i,j,k)-GP(im,jm,km))*idmg)
                 ! this%gradLP(idir,i,j,k)=mmgrad((LP(ip,jp,kp)-LP(i,j,k))*idpl,(LP(i,j,k)-LP(im,jm,km))*idml)
              end do
           end do
        end do
     end do

     ! Communication and BCs
     
   contains

     ! Minmod gradient
     function mmgrad(g1,g2) result(g)
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
     end function mmgrad

   end subroutine flow_reconstruct


   ! Full advection routine inside the inner loop
   subroutine advection_step(this,dt,vf,matmod)
     use vfs_class,  only: vfs, VFhi, VFlo
     use matm_class, only: matm
     use irl_fortran_interface, only : CapDod_type,TagAccVM_SepVM_type,new
     implicit none
     class(mast), intent(inout) :: this   !< The two-phase all-Mach flow solver
     class(vfs),  intent(inout) :: vf     !< The volume fraction solver
     class(matm), intent(inout) :: matmod !< The material models for this solver
     real(WP),    intent(inout) :: dt     !< Timestep size over which to advance
     real(WP),   dimension(14)  :: flux   !< Passes flux to and from routines
     real(WP), dimension(:,:,:), pointer :: PgradX,PgradY,PgradZ
     type(CapDod_type) :: fp              !< Object for flux polyhedron
     type(TagAccVM_SepVM_type) :: ffm     !< Object for flux moments
     real(WP), dimension(0:1) :: jump
     real(WP) :: Ga_i,Ga_nb,La_i,La_nb
     real(WP) :: rho_l,rho_r
     integer  :: i,j,k

     ! Initialize flux sum arrays
     this%F_VOL  =                                 this%cfg%vol
     this%F_VF   =                       vf%VFold *this%cfg%vol
     this%F_Grho =this%Grhoold *((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoE=this%GrhoEold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoU=this%Grhoold*this%Uiold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoV=this%Grhoold*this%Viold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GrhoW=this%Grhoold*this%Wiold*((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_GP   =this%GPold   *((1.0_WP-vf%VFold)*this%cfg%vol)
     this%F_Lrho =this%Lrhoold *((       vf%VFold)*this%cfg%vol)
     this%F_LrhoE=this%LrhoEold*((       vf%VFold)*this%cfg%vol)
     this%F_LrhoU=this%Lrhoold*this%Uiold*((       vf%VFold)*this%cfg%vol)
     this%F_LrhoV=this%Lrhoold*this%Viold*((       vf%VFold)*this%cfg%vol)
     this%F_LrhoW=this%Lrhoold*this%Wiold*((       vf%VFold)*this%cfg%vol)
     this%F_LP   =this%LPold   *((       vf%VFold)*this%cfg%vol)

     ! Allocate flux_polyhedron that will be used for fluxes
     call new(fp)
     ! Allocate face_flux_moment that will be used in fluxes
     call new(ffm)

     ! Designate the use of temporary arrays for pressure gradients
     PgradX => this%tmp1; PgradY => this%tmp2; PgradZ =>this%tmp3
     PgradX = 0.0_WP;     PgradY = 0.0_WP;     PgradZ = 0.0_WP


     !! ---------------------------------------!!
     !! 1. SL and TTSL flux calculations       !!
     !! ---------------------------------------!!
     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
              
              !! ---- LEFT X(I) FACE ---- !!
              select case(this%sl_face(i,j,k,1))
              case(1); call SL_advect  (flux,fp,ffm                                        ,i,j,k,'x')
              case(0); call TTSL_advect(flux,[this%cfg%x (i),this%cfg%ym(j),this%cfg%zm(k)],i,j,k,'x')
              end select
              call add_fluxes(flux,i,j,k,'x')

              !! ---- BOTTOM Y(J) FACE ---- !!
              select case(this%sl_face(i,j,k,2))
              case(1); call SL_advect  (flux,fp,ffm                                        ,i,j,k,'y')
              case(0); call TTSL_advect(flux,[this%cfg%xm(i),this%cfg%y (j),this%cfg%zm(k)],i,j,k,'y')
              end select
              call add_fluxes(flux,i,j,k,'y')

              !! ---- BACK Z(K) FACE ---- !!
              select case(this%sl_face(i,j,k,3))
              case(1); call SL_advect  (flux,fp,ffm                                        ,i,j,k,'z')
              case(0); call TTSL_advect(flux,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%z (k)],i,j,k,'z')
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
              ! Skip wall cells and masked BC cells
              if (this%mask(i,j,k).gt.0) cycle
              ! Update VF
              vf%VF(i,j,k) = this%F_VF(i,j,k)/this%F_VOL(i,j,k)
              ! Update phase density, pressure, and bulk moduli according to VF limits
              if (vf%VF(i,j,k).lt.VFlo) then
                 ! Gas only
                 vf%VF(i,j,k) = 0.0_WP
                 ! Calculate PA, bulkmod using primitively advected density
                 this%Grho (i,j,k)   = this%F_Grho (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%GP   (i,j,k)   = this%F_GP   (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%GrhoSS2(i,j,k) = matmod%EOS_gas   (i,j,k,'K')
                 ! Zero quantities in opposite phase
                 this%Lrho (i,j,k)   = 0.0_WP
                 this%LP   (i,j,k)   = 0.0_WP
                 this%LrhoSS2(i,j,k) = 0.0_WP
                 ! Calculate density incorporating volume change
                 this%Grho (i,j,k)   = this%F_Grho (i,j,k)/this%cfg%vol(i,j,k)
              else if (vf%VF(i,j,k).gt.VFhi) then
                 ! Liquid only
                 vf%VF  (i,j,k) = 1.0_WP
                 ! Calculate PA, bulkmod using primitively advected density
                 this%Lrho (i,j,k)   = this%F_Lrho (i,j,k)/this%F_VF(i,j,k)
                 this%LP   (i,j,k)   = this%F_LP   (i,j,k)/this%F_VF(i,j,k)
                 this%LrhoSS2(i,j,k) = matmod%EOS_liquid(i,j,k,'K')
                 ! Zero quantities in opposite phase
                 this%Grho (i,j,k)   = 0.0_WP
                 this%GP   (i,j,k)   = 0.0_WP
                 this%GrhoSS2(i,j,k) = 0.0_WP
                 ! Calculate density incorporating volume change
                 this%Lrho (i,j,k)   = this%F_Lrho (i,j,k)/this%cfg%vol(i,j,k)
              else
                 ! Primitive density, advected pressure, bulkmod for both phases
                 this%Grho (i,j,k)   = this%F_Grho (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%Lrho (i,j,k)   = this%F_Lrho (i,j,k)/(                  this%F_VF(i,j,k))
                 this%GP   (i,j,k)   = this%F_GP   (i,j,k)/(this%F_VOL(i,j,k)-this%F_VF(i,j,k))
                 this%LP   (i,j,k)   = this%F_LP   (i,j,k)/(                  this%F_VF(i,j,k))
                 this%GrhoSS2(i,j,k) = matmod%EOS_gas   (i,j,k,'K')
                 this%LrhoSS2(i,j,k) = matmod%EOS_liquid(i,j,k,'K')
                 ! Store current VF
                 this%srcVF(i,j,k)   = vf%VF(i,j,k)
                 ! Get new VF according to quadratic source term that depends on compressibility
                 vf%VF(i,j,k) = VF_src_quad(this%GrhoSS2(i,j,k), this%LrhoSS2(i,j,k), &
                      this%F_VF(i,j,k), this%F_VOL(i,j,k), this%cfg%vol(i,j,k))
                 ! Update phase quantities according to VF limits (again)
                 if (vf%VF(i,j,k).lt.VFlo) then
                    vf%VF(i,j,k)        = 0.0_WP
                    this%srcVF(i,j,k)   = 0.0_WP
                    ! Zero density, PA, bulkmod, and dpde_rho of other phase
                    this%Lrho (i,j,k)   = 0.0_WP
                    this%LP   (i,j,k)   = 0.0_WP
                    this%LrhoSS2(i,j,k) = 0.0_WP
                    ! Calculate density, incorporating volume change
                    this%Grho (i,j,k)   = this%F_Grho (i,j,k)/this%cfg%vol(i,j,k)
                 else if (vf%VF(i,j,k).gt.VFhi) then
                    vf%VF(i,j,k)        = 1.0_WP
                    this%srcVF(i,j,k)   = 0.0_WP
                    ! Zero density, PA, bulkmod, and dpde_rho of other phase
                    this%Grho (i,j,k)   = 0.0_WP
                    this%GP   (i,j,k)   = 0.0_WP
                    this%GrhoSS2(i,j,k) = 0.0_WP
                    ! Calculate density, incorporating volume change
                    this%Lrho (i,j,k)   = this%F_Lrho (i,j,k)/this%cfg%vol(i,j,k)
                 else
                    ! Record srcVF
                    this%srcVF(i,j,k)   = vf%VF(i,j,k)-this%srcVF(i,j,k)
                    ! PA, bulkmod stay unchanged
                    ! Calculate density, incorporating volume change
                    this%Grho (i,j,k)   = this%F_Grho (i,j,k)/((1.0_WP-vf%VF(i,j,k))*this%cfg%vol(i,j,k))
                    this%Lrho (i,j,k)   = this%F_Lrho (i,j,k)/(        vf%VF(i,j,k) *this%cfg%vol(i,j,k))
                 end if
              end if

           end do
        end do
     end do

     ! Boundaries for VF, phase density, phase pressure
     call vf%cfg%sync(vf%VF)
     
     ! Update phase interface to match VF field, etc.
     ! (Mimics the end of vf%advance())
     call vf%remove_flotsams()
     call vf%sync_and_clean_barycenters()
     call vf%advect_interface(dt,this%U,this%V,this%W)
     call vf%build_interface()
     !call vf%update_band()
     call vf%remove_sheets()
     call vf%polygonalize_interface()
     !call vf%distance_from_polygon()
     call vf%subcell_vol()
     call vf%get_curvature()
     ! Band and distance information should not be needed

     ! Boundary conditions density
     bc_scope = 'density'
     call this%apply_bcond(dt,bc_scope)
     ! Sync for phase pressure (boundaries are not used)
     call this%cfg%sync(this%GP)
     call this%cfg%sync(this%LP)

     ! Mixture density
     this%RHO = vf%VF*this%Lrho + (1.0_WP-vf%VF)*this%Grho

     ! Interpolate pressure (for predictor) and density (for predictor and Helmholtz equation operator)
     call this%interp_pressure_density(vf)
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
                    this%implicit%opr(1,i-1,j,k) = this%implicit%opr(1,i-1,j,k) - (1.0_WP-vf%VF(i-1,j,k))*Ga_nb- vf%VF(i-1,j,k)*La_nb
                    this%implicit%opr(2,i-1,j,k) = this%implicit%opr(2,i-1,j,k) - (1.0_WP-vf%VF(i-1,j,k))*Ga_i - vf%VF(i-1,j,k)*La_i
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
                    this%implicit%opr(1,i,j-1,k) = this%implicit%opr(1,i,j-1,k) - (1.0_WP-vf%VF(i,j-1,k))*Ga_nb- vf%VF(i,j-1,k)*La_nb
                    this%implicit%opr(4,i,j-1,k) = this%implicit%opr(4,i,j-1,k) - (1.0_WP-vf%VF(i,j-1,k))*Ga_i - vf%VF(i,j-1,k)*La_i
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
                    this%implicit%opr(1,i,j,k-1) = this%implicit%opr(1,i,j,k-1) - (1.0_WP-vf%VF(i,j,k-1))*Ga_nb- vf%VF(i,j,k-1)*La_nb
                    this%implicit%opr(6,i,j,k-1) = this%implicit%opr(6,i,j,k-1) - (1.0_WP-vf%VF(i,j,k-1))*Ga_i - vf%VF(i,j,k-1)*La_i
                 end if
              end if

           end do
        end do
     end do

     !! Pressure predictor - momentum
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! Pressure jump, gradient for x
              rho_l=sum(vf%Gvol(0,:,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(0,:,:,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(1,:,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(1,:,:,i,j,k))*this%Lrho(i,j,k)
              jump(0)=rho_l*this%Pjx(i  ,j,k)/(this%rho_U(i  ,j,k)*this%cfg%vol(i,j,k))
              jump(1)=rho_r*this%Pjx(i+1,j,k)/(this%rho_U(i+1,j,k)*this%cfg%vol(i,j,k))
              PgradX(i,j,k) = (this%P_U(i+1,j,k)-this%P_U(i,j,k))*this%cfg%dxi(i)-sum(jump)
              ! Pressure jump, gradient for y
              rho_l=sum(vf%Gvol(:,0,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,0,:,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(:,1,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,1,:,i,j,k))*this%Lrho(i,j,k)
              jump(0)=rho_l*this%Pjy(i,j  ,k)/(this%rho_V(i,j  ,k)*this%cfg%vol(i,j,k))
              jump(1)=rho_r*this%Pjy(i,j+1,k)/(this%rho_V(i,j+1,k)*this%cfg%vol(i,j,k))
              PgradY(i,j,k) = (this%P_V(i,j+1,k)-this%P_V(i,j,k))*this%cfg%dyi(j)-sum(jump)
              ! Pressure jump, gradient for z
              rho_l=sum(vf%Gvol(:,:,0,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,:,0,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(:,:,1,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,:,1,i,j,k))*this%Lrho(i,j,k)
              jump(0)=rho_l*this%Pjz(i,j,k  )/(this%rho_W(i,j,k  )*this%cfg%vol(i,j,k))
              jump(1)=rho_r*this%Pjz(i,j,k+1)/(this%rho_W(i,j,k+1)*this%cfg%vol(i,j,k))
              PgradZ(i,j,k) = (this%P_W(i,j,k+1)-this%P_W(i,j,k))*this%cfg%dzi(k)-sum(jump)
           end do
        end do
     end do


     ! Solve for x-momentum
     call this%implicit%setup()
     this%implicit%rhs=this%F_GrhoU+this%F_LrhoU-dt*this%cfg%vol*PgradX
     this%implicit%sol=this%rhoUi
     call this%implicit%solve()
     this%rhoUi=this%implicit%sol
     ! Solve for y-momentum
     call this%implicit%setup()
     this%implicit%rhs=this%F_GrhoV+this%F_LrhoV-dt*this%cfg%vol*PgradY
     this%implicit%sol=this%rhoVi
     call this%implicit%solve()
     this%rhoVi=this%implicit%sol
     ! Solve for z-momentum
     call this%implicit%setup()
     this%implicit%rhs=this%F_GrhoW+this%F_LrhoW-dt*this%cfg%vol*PgradZ
     this%implicit%sol=this%rhoWi
     call this%implicit%solve()
     this%rhoWi=this%implicit%sol

     ! Nullify pointers
     nullify(PgradX,PgradY,PgradZ)

     ! Boundary conditions for momentum
     bc_scope = 'momentum'
     call this%apply_bcond(dt,bc_scope)

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
              
              ! Cycle wall cells and masked BC cells
              if (this%mask(i,j,k).gt.0) cycle
              
              ! Incorporate pressure predictor terms
              this%F_GrhoE(i,j,k) = this%F_GrhoE(i,j,k) &
                   - dt*this%cfg%vol(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(this%Grho(i,j,k)/this%RHO(i,j,k)*&
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)*this%P_U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)*this%P_V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)*this%P_W(i,j,k:k+1)) &
                   - sum(this%Pjx (i:i+1,j,k)*this%U(i:i+1,j,k)) &
                   - sum(this%Pjy (i,j:j+1,k)*this%V(i,j:j+1,k)) &
                   - sum(this%Pjz (i,j,k:k+1)*this%W(i,j,k:k+1)) ) + &
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)) )*&
                   ( this%P(i,j,k)*(1.0_WP-this%Grho(i,j,k)/this%RHO(i,j,k)) &
                   - this%Hpjump(i,j,k)*(       vf%VF(i,j,k))) )

	      this%F_LrhoE(i,j,k) = this%F_LrhoE(i,j,k) &
                   - dt*this%cfg%vol(i,j,k)*(       vf%VF(i,j,k))*(this%Lrho(i,j,k)/this%RHO(i,j,k)*&
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)*this%P_U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)*this%P_V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)*this%P_W(i,j,k:k+1)) &
                   - sum(this%Pjx (i:i+1,j,k)*this%U(i:i+1,j,k)) &
                   - sum(this%Pjy (i,j:j+1,k)*this%V(i,j:j+1,k)) &
                   - sum(this%Pjz (i,j,k:k+1)*this%W(i,j,k:k+1)) ) + &
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)) )*&
                   ( this%P(i,j,k)*(1.0_WP-this%Lrho(i,j,k)/this%RHO(i,j,k)) &
                   + this%Hpjump(i,j,k)*(1.0_WP-vf%VF(i,j,k))) )              
              
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
     bc_scope = 'energy'
     call this%apply_bcond(dt,bc_scope)
  
   contains

     ! ================================================================= !
     ! Construct flux hexahedron, perform cutting, call flux calculation !
     ! ================================================================= !
     subroutine SL_advect(flux,a_flux_polyhedron,some_face_flux_moments,i,j,k,dir)
       character(len=1) :: dir
       integer :: i,j,k,idir
       real(WP), dimension(14) :: flux

       type(CapDod_type) :: a_flux_polyhedron
       type(TagAccVM_SepVM_type) :: some_face_flux_moments

       select case (trim(dir))
       case ('x')
          idir = 1
       case ('y')
          idir = 2
       case('z')
          idir = 3
       end select

       ! Prepare flux polyhedron
       call vf%fluxpoly_project_getmoments(i,j,k,dt,dir,this%U,this%V,this%W,&
            a_flux_polyhedron,vf%localized_separator_linkold(i,j,k),some_face_flux_moments)
       ! Calculate fluxes from volume moments
       call SL_getFaceFlux(some_face_flux_moments, flux)
       ! Store face density terms
       this%GrhoFf(i,j,k,idir) = flux(3)
       this%LrhoFf(i,j,k,idir) = flux(8)
       
       return
     end subroutine SL_advect

     ! ====================================================== !
     ! Given a flux hexahedron, calculate and return the flux !
     ! ====================================================== !
     subroutine SL_getFaceFlux(f_moments, flux)
       use irl_fortran_interface, only: getSize
       integer  :: ii,jj,kk,n
       integer  :: list_size,uniq_id
       integer,  dimension(3) :: ind
       real(WP) :: my_Gvol,my_Lvol
       real(WP), dimension(3) :: my_Gbary,my_Lbary
       real(WP), dimension(14)  :: flux
       type(TagAccVM_SepVM_type) :: f_moments
       logical  :: skip_flag

       !..... Using geometry, calculate fluxes .....!
       ! Initialize at zero
       flux = 0.0_WP
       ! Get number of tags (elements) in the f_moments
       list_size = getSize(f_moments)
       ! Loop through tags in the list
       do n = 0,list_size-1
          ! Get indices of current cell, volumes, and centroids
          call vf%fluxpoly_cell_getvolcentr(f_moments,n,ii,jj,kk, &
               my_Lbary,my_Gbary,my_Lvol,my_Gvol,skip_flag)

          ! Skip current cell if there is a reason
          if (skip_flag) cycle

          ! Store barycenter fluxes
          !b_flux(:,1) = b_flux(:,1) + my_Gvol*my_Gbary
          !b_flux(:,2) = b_flux(:,2) + my_Lvol*my_Lbary

          ! Bound barycenters by the cell dimensions
          my_Lbary(1) = max(this%cfg%x(ii),min(this%cfg%x(ii+1),my_Lbary(1)))
          my_Gbary(1) = max(this%cfg%x(ii),min(this%cfg%x(ii+1),my_Gbary(1)))
          my_Lbary(2) = max(this%cfg%y(jj),min(this%cfg%y(jj+1),my_Lbary(2)))
          my_Gbary(2) = max(this%cfg%y(jj),min(this%cfg%y(jj+1),my_Gbary(2)))
          my_Lbary(3) = max(this%cfg%z(kk),min(this%cfg%z(kk+1),my_Lbary(3)))
          my_Gbary(3) = max(this%cfg%z(kk),min(this%cfg%z(kk+1),my_Gbary(3)))

          ! Make barycenter relative to cell
          my_Lbary = my_Lbary - vf%Lbaryold(:,ii,jj,kk)
          my_Gbary = my_Gbary - vf%Gbaryold(:,ii,jj,kk)

          ! Add contribution to the flux in this tetrahedron
          flux( 1) = flux( 1) + my_Lvol + my_Gvol
          flux( 2) = flux( 2) + my_Lvol

          flux( 3) = flux( 3) + my_Gvol*(this%Grhoold (ii,jj,kk)+sum(this%gradGrho (:,ii,jj,kk)*my_Gbary(:)))
          flux( 4) = flux( 4) + my_Gvol*(this%GrhoEold(ii,jj,kk)+sum(this%gradGrhoE(:,ii,jj,kk)*my_Gbary(:)))
          flux( 5) = flux( 5) + my_Gvol*(this%Grhoold(ii,jj,kk)*this%Uiold(ii,jj,kk)+sum(this%gradGrhoU(:,ii,jj,kk)*my_Gbary(:)))
          flux( 6) = flux( 6) + my_Gvol*(this%Grhoold(ii,jj,kk)*this%Viold(ii,jj,kk)+sum(this%gradGrhoV(:,ii,jj,kk)*my_Gbary(:)))
          flux( 7) = flux( 7) + my_Gvol*(this%Grhoold(ii,jj,kk)*this%Wiold(ii,jj,kk)+sum(this%gradGrhoW(:,ii,jj,kk)*my_Gbary(:)))

          flux( 8) = flux( 8) + my_Lvol*(this%Lrhoold (ii,jj,kk)+sum(this%gradLrho (:,ii,jj,kk)*my_Lbary(:)))
          flux( 9) = flux( 9) + my_Lvol*(this%LrhoEold(ii,jj,kk)+sum(this%gradLrhoE(:,ii,jj,kk)*my_Lbary(:)))
          flux(10) = flux(10) + my_Lvol*(this%Lrhoold(ii,jj,kk)*this%Uiold(ii,jj,kk)+sum(this%gradLrhoU(:,ii,jj,kk)*my_Lbary(:)))
          flux(11) = flux(11) + my_Lvol*(this%Lrhoold(ii,jj,kk)*this%Viold(ii,jj,kk)+sum(this%gradLrhoV(:,ii,jj,kk)*my_Lbary(:)))
          flux(12) = flux(12) + my_Lvol*(this%Lrhoold(ii,jj,kk)*this%Wiold(ii,jj,kk)+sum(this%gradLrhoW(:,ii,jj,kk)*my_Lbary(:)))

          flux(13) = flux(13) + my_Gvol*(this%GPold   (ii,jj,kk)+sum(this%gradGP   (:,ii,jj,kk)*my_Gbary(:)))
          flux(14) = flux(14) + my_Lvol*(this%LPold   (ii,jj,kk)+sum(this%gradLP   (:,ii,jj,kk)*my_Lbary(:)))

       end do

       return
     end subroutine SL_getFaceFlux

     ! ===================================================== !
     ! Project the center of the face, call flux calculation !
     ! ===================================================== !
     subroutine TTSL_advect(flux,pt_f,i,j,k,dir)
       implicit none
       real(WP), dimension(3) :: pt_f,pt_p
       real(WP), dimension(14):: flux
       integer,  dimension(3) :: ind_p
       real(WP) :: vol_f,vol_check
       character(len=1) :: dir
       integer :: n,i,j,k,idir

       ! Trajectory-Tracing Semi-Lagrangian scheme

       ! Initialize flux at 0
       flux = 0.0_WP

       ! Get projected point
       pt_p = vf%project(pt_f,i,j,k,-dt,this%U,this%V,this%W)

       ! Get approximate barycenter of flux volume, assign to one phase
       ! b_flux(:,ceiling(oldVOF(i,j,k))+1) = 0.5_WP*(pt_p+pt_f)

       select case (trim(dir))
       case('x')
          idir = 1
          ! Flux volume
          vol_f = dt*this%U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
          ! Ratio of distance to check
          vol_check = abs(dt*this%U(i,j,k)*this%cfg%dxi(i))
       case('y')
          idir = 2
          ! Flux volume
          vol_f = dt*this%V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k)
          ! Ratio of distance to check
          vol_check = abs(dt*this%V(i,j,k)*this%cfg%dyi(j))
       case('z')
          idir = 3
          ! Flux volume
          vol_f = dt*this%W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j)
          ! Ratio of distance to check
          vol_check = abs(dt*this%W(i,j,k)*this%cfg%dzi(k))
       end select
       
       ! Calculate barycenter fluxes
       !b_flux = b_flux*vol_f
       
       ! Calculate fluxes (if there is any flux to calculate)
       if (vol_check.gt.epsilon(1.0_WP)) then
          ! Get initial indices
          ind_p = this%cfg%get_ijk_local(pt_p,[i,j,k])
          call TTSL_getFaceFlux(pt_p,pt_f,ind_p,vol_f,  flux)
       end if
       ! Store face density terms
       this%GrhoFf(i,j,k,idir) = flux(3)
       this%LrhoFf(i,j,k,idir) = flux(8)
       
       return
     end subroutine TTSL_advect

     ! ===================================================== !
     ! Calculate fluxes using amount of line segment in cell !
     ! ===================================================== !
     subroutine TTSL_getFaceFlux(pt_p,pt_f,ind_p,vol_f, flux)
       implicit none
       real(WP), dimension(3) :: pt_f,pt_p,pt_c,dx_f
       real(WP), dimension(14):: flux
       integer,  dimension(3) :: ind_p,ind_c
       real(WP) :: vol_f,l_f,l_c,l_d,l_s

       ! Total distance from projected face centroid to initial face centroid
       l_f = sqrt(sum((pt_p-pt_f)**2))
       ! Initial difference between integrated distance and total distance
       l_d = l_f
       ! Initial amount of length that has been summed over
       l_s = 0.0_WP
       ! Loop through cells while progressing from pt_p to pt_f
       do while (l_d.gt.epsilon(1.0_WP))
          ! Find intersection between flux line and mesh, get next point and indices
          call mesh_intersect(pt_p,pt_f,ind_p,   pt_c,ind_c)
          ! Find distance between intersection and origin face
          l_d = sqrt(sum((pt_c-pt_f)**2))
          ! Get distance between point p and intersection (current distance)
          l_c = l_f - l_d - l_s
          ! Calculate flux, multiply by distance ratio
          flux = flux + l_c/(l_f+tiny(1.0_WP))*vol_f*TTSL_getval(ind_p,0.5_WP*(pt_p+pt_c))
          ! Move point p to intersection
          pt_p = pt_c
          ! Move indices to next location
          ind_p = ind_c
          ! Add length that has been summed over
          l_s = l_s + l_c
       end do

     end subroutine TTSL_getFaceFlux
     
     ! Find value at point on line segment
     function TTSL_getval(ind,pt_i) result(f)

       real(WP), dimension(3),  intent(in)  :: pt_i
       integer,  dimension(3),  intent(in)  :: ind
       real(WP), dimension(14) :: f
       real(WP), dimension(3) :: dx_i
       integer :: ii,jj,kk

       ii = ind(1); jj = ind(2); kk = ind(3)
       f = 0.0_WP

       ! Avoid if in wall
       if (this%mask(ii,jj,kk).eq.1) return
       ! Avoid if beyond outflow boundary
       ! if (backflow_flux_flag(ii,jj,kk)) return

       ! normalized volume flux
       f(1) = 1.0_WP
       select case(ceiling(vf%VFold(ii,jj,kk)))
       case(0)
          ! displacement to barycenter
          dx_i = pt_i - vf%Gbaryold(:,ii,jj,kk)
          ! interpolated variables
          f(3) = this%Grhoold (ii,jj,kk)                       + sum(this%gradGrho(:,ii,jj,kk)*dx_i)
          f(4) = this%GrhoEold(ii,jj,kk)-this%GKEold(ii,jj,kk) + sum(this%gradGIE (:,ii,jj,kk)*dx_i)
          f(13)= this%GPold   (ii,jj,kk)                       + sum(this%gradGP  (:,ii,jj,kk)*dx_i)
       case(1)
          ! displacement to barycenter
          dx_i = pt_i - vf%Lbaryold(:,ii,jj,kk)
          ! interpolated variables
          f(2) = 1.0_WP
          f(8) = this%Lrhoold (ii,jj,kk)                       + sum(this%gradLrho(:,ii,jj,kk)*dx_i)
          f(9) = this%LrhoEold(ii,jj,kk)-this%LKEold(ii,jj,kk) + sum(this%gradLIE (:,ii,jj,kk)*dx_i)
          f(14)= this%LPold   (ii,jj,kk)                       + sum(this%gradLP  (:,ii,jj,kk)*dx_i)
       end select
     end function TTSL_getval

     ! Routine to find intesections with mesh and calculate distances
     subroutine mesh_intersect(pt_p,pt_f,ind_p, pt_c,ind_c)
       integer :: i,j,k,min_d
       real(WP) :: slope_factor
       real(WP), dimension(3) :: pt_p,pt_f,pt_c,vec
       integer,  dimension(3) :: ind_p,ind_c,ind_m
       real(WP), dimension(3) :: xint,yint,zint
       real(WP), dimension(4) :: dint
       logical :: p_flag

       ! Direction of intersection
       vec = pt_f - pt_p
       ! Identify indices of mesh planes that could be intersected
       ind_m = ind_p + nint(0.5_WP*(1.0_WP - sign(1.0_WP,-pt_f+pt_p)))
       ! Identify indices of potential next cells
       ind_c = ind_p - nint(sign(1.0_WP,-pt_f+pt_p))
       ! Signs are weird because I need sign(0) => -1
       ! Calculate intersection point for each mesh plane
       ! X-plane
       xint(1) = this%cfg%x(ind_m(1))                   ! x-coord at intersection
       slope_factor = (xint(1)-pt_f(1))/(vec(1)+tiny(1.0_WP))
       slope_factor = min(sqrt(huge(1.0_WP)),max(slope_factor,-sqrt(huge(1.0_WP))))
       xint(2) = pt_p(2) + vec(2)*(1.0_WP+slope_factor) ! y-coord at intersection
       xint(3) = pt_p(3) + vec(3)*(1.0_WP+slope_factor) ! z-coord at intersection
       ! Y-plane
       yint(2) = this%cfg%y(ind_m(2))                   ! y-coord at intersection
       slope_factor = (yint(2)-pt_f(2))/(vec(2)+tiny(1.0_WP))
       slope_factor = min(sqrt(huge(1.0_WP)),max(slope_factor,-sqrt(huge(1.0_WP))))
       yint(1) = pt_p(1) + vec(1)*(1.0_WP+slope_factor) ! x-coord at intersection
       yint(3) = pt_p(3) + vec(3)*(1.0_WP+slope_factor) ! z-coord at intersection
       ! Z-plane
       zint(3) = this%cfg%z(ind_m(3))                   ! z-coord at intersection
       slope_factor = (zint(3)-pt_f(3))/(vec(3)+tiny(1.0_WP))
       slope_factor = min(sqrt(huge(1.0_WP)),max(slope_factor,-sqrt(huge(1.0_WP))))
       zint(1) = pt_p(1) + vec(1)*(1.0_WP+slope_factor) ! x-coord at intersection
       zint(2) = pt_p(2) + vec(2)*(1.0_WP+slope_factor) ! y-coord at intersection
       ! If slope is zero, the equation should become a constant
       ! Calculate distance of line segment, find minimum
       dint(1) = sqrt(sum((pt_p-xint)**2)) ! distance to x intersection
       dint(2) = sqrt(sum((pt_p-yint)**2)) ! distance to y intersection
       dint(3) = sqrt(sum((pt_p-zint)**2)) ! distance to z intersection
       dint(4) = sqrt(sum((pt_p-pt_f)**2)) ! distance to origin face point
       min_d = minloc(dint,1)

       ! Update new point location, reset other indices to keep only correct new index
       select case (min_d)
       case(1) ! x intersection
          ind_c(2) = ind_p(2); ind_c(3) = ind_p(3)
          pt_c = xint
       case(2) ! y intersection
          ind_c(1) = ind_p(1); ind_c(3) = ind_p(3)
          pt_c = yint
       case(3) ! z intersection
          ind_c(1) = ind_p(1); ind_c(2) = ind_p(2)
          pt_c = zint
       case(4) ! origin point is closest
          ind_c = ind_p ! no need to change indices
          pt_c = pt_f
       end select

       return
     end subroutine mesh_intersect

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
          grho_fterm = this%GrhoFf(i,j,k,1)
          lrho_fterm = this%LrhoFf(i,j,k,1)
       case('y')
          st_i = 0; st_j =-1; st_k = 0
          grho_fterm = this%GrhoFf(i,j,k,2)
          lrho_fterm = this%LrhoFf(i,j,k,2)
       case('z')
          st_i = 0; st_j = 0; st_k =-1
          grho_fterm = this%GrhoFf(i,j,k,3)
          lrho_fterm = this%LrhoFf(i,j,k,3)
       end select

       !print*,'before flux calculations'
       ! Gas phase
       flux(5) = SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%Grhoold(i,j,k)*this%Uiold(i,j,k),this%Grhoold(i+st_i,j+st_j,k+st_k)*this%Uiold(i+st_i,j+st_j,k+st_k),&
            grho_fterm)
       flux(6) = SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%Grhoold(i,j,k)*this%Viold(i,j,k),this%Grhoold(i+st_i,j+st_j,k+st_k)*this%Viold(i+st_i,j+st_j,k+st_k),&
            grho_fterm)
       flux(7) = SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%Grhoold(i,j,k)*this%Wiold(i,j,k),this%Grhoold(i+st_i,j+st_j,k+st_k)*this%Wiold(i+st_i,j+st_j,k+st_k),&
            grho_fterm)
       ! Liquid phase
       flux(10)= SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%Lrhoold(i,j,k)*this%Uiold(i,j,k),this%Lrhoold(i+st_i,j+st_j,k+st_k)*this%Uiold(i+st_i,j+st_j,k+st_k),&
            lrho_fterm)
       flux(11)= SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%Lrhoold(i,j,k)*this%Viold(i,j,k),this%Lrhoold(i+st_i,j+st_j,k+st_k)*this%Viold(i+st_i,j+st_j,k+st_k),&
            lrho_fterm)
       flux(12)= SubCan_advect_mom_RHS(&
            this%RHO     (i,j,k),this%RHO     (i+st_i,j+st_j,k+st_k),&
            this%RHOold  (i,j,k),this%RHOold  (i+st_i,j+st_j,k+st_k),&
            this%Lrhoold(i,j,k)*this%Wiold(i,j,k),this%Lrhoold(i+st_i,j+st_j,k+st_k)*this%Wiold(i+st_i,j+st_j,k+st_k),&
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
          grho_fterm = -this%GrhoFf(i,j,k,1)
          lrho_fterm = -this%LrhoFf(i,j,k,1)
       case('y')
          st_i = 0; st_j =-1; st_k = 0
          grho_fterm = -this%GrhoFf(i,j,k,2)
          lrho_fterm = -this%LrhoFf(i,j,k,2)
       case('z')
          st_i = 0; st_j = 0; st_k =-1
          grho_fterm = -this%GrhoFf(i,j,k,3)
          lrho_fterm = -this%LrhoFf(i,j,k,3)
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
          grho_fterm = this%GrhoFf(i,j,k,1)
          lrho_fterm = this%LrhoFf(i,j,k,1)
       case('y')
          st_i = 0; st_j =-1; st_k = 0
          grho_fterm = this%GrhoFf(i,j,k,2)
          lrho_fterm = this%LrhoFf(i,j,k,2)
       case('z')
          st_i = 0; st_j = 0; st_k =-1
          grho_fterm = this%GrhoFf(i,j,k,3)
          lrho_fterm = this%LrhoFf(i,j,k,3)
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

   subroutine pressureproj_prepare(this,dt,vf,matmod,contact_model)
     use vfs_class, only : vfs
     use matm_class, only: matm
     implicit none
     class(mast), intent(inout) :: this
     class(vfs),  intent(in)    :: vf
     class(matm), intent(inout) :: matmod
     real(WP),    intent(in)    :: dt
     integer, intent(in), optional :: contact_model
     real(WP), dimension(:,:,:), pointer :: termU,termV,termW

     ! Designate the use of temporary arrays
     termU => this%tmp1; termV => this%tmp2; termW =>this%tmp3
     ! Calculate terms for modified velocity interpolation
     call this%terms_modified_face_velocity(vf,termU,termV,termW)
     ! Perform modified velocity interpolation
     call this%interp_vel_full(vf,dt,this%Ui,this%Vi,this%Wi,termU,termV,termW,this%U,this%V,this%W)
     ! Disassociate pointers at end of use
     nullify(termU,termV,termW)
     ! Address BCs of velocity field
     bc_scope = 'velocity'
     call this%apply_bcond(dt,bc_scope)
     ! Calculate mixed advected pressure and bulk modulus
     call this%harmonize_advpressure_bulkmod(vf,matmod)
     ! Calculate pressure jump
     call this%update_surface_tension_jump(vf,contact_model)
     ! Calculate LHS operator of Helmholtz equation
     call this%update_Helmholtz_LHS(dt)
     ! Calculate RHS of Helmholtz equation
     call this%update_Helmholtz_RHS(dt)
     ! Boundary conditions will be needed for pressure equation if dirichlet or fixed gradient is used.
     ! These boundary conditions could be applied by the user directly, though.
     ! (For Neumann pressure BCs, nothing needs to be done)

   end subroutine pressureproj_prepare

   subroutine pressureproj_correct(this,dt,vf,DP)
     use vfs_class, only : vfs
     implicit none
     class(mast), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     real(WP),    intent(in)    :: dt
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: DP
     real(WP), dimension(:,:,:), pointer :: DP_U,DP_V,DP_W
     real(WP) :: vol_l,vol_r,rho_l,rho_r
     real(WP), dimension(0:1) :: jumpx,jumpy,jumpz,rho_f
     integer :: i,j,k
     ! Set up temporary array
     DP_U=>this%tmp1; DP_V=>this%tmp2; DP_W=>this%tmp3

     ! Correct face velocity and calculate incremental face pressure
     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
               ! Update face pressure and density in X
               rho_f(1)=0.0_WP; vol_r=sum(vf%Gvol(0,:,:,i  ,j,k)+vf%Lvol(0,:,:,i  ,j,k))
               if (vol_r.gt.0.0_WP.and.this%mask(i  ,j,k).eq.0) &
                    rho_f(1)=(sum(vf%Gvol(0,:,:,i  ,j,k))*this%Grho(i  ,j,k)&
                             +sum(vf%Lvol(0,:,:,i  ,j,k))*this%Lrho(i  ,j,k))
               rho_f(0)=0.0_WP; vol_l=sum(vf%Gvol(1,:,:,i-1,j,k)+vf%Lvol(1,:,:,i-1,j,k))
               if (vol_l.gt.0.0_WP.and.this%mask(i-1,j,k).eq.0) &
                    rho_f(0)=(sum(vf%Gvol(1,:,:,i-1,j,k))*this%Grho(i-1,j,k)&
                             +sum(vf%Lvol(1,:,:,i-1,j,k))*this%Lrho(i-1,j,k))
               if (sum(rho_f).gt.0.0_WP) then
                  DP_U(i,j,k)=2.0_WP*sum(this%itpi_x(:,i,j,k)*DP(i-1:i,j,k))&
                       -sum(this%itpi_x(:,i,j,k)*rho_f*DP(i-1:i,j,k)) / &
                       (sum(this%itpi_x(:,i,j,k)*rho_f) + tiny(1.0_WP))
                  this%U(i,j,k)=this%U(i,j,k)-dt*(&
                       sum(this%divu_x(:,i,j,k)*DP(i-1:i,j,k))&
                       -this%DPjx(i,j,k))/this%rho_U(i,j,k)
               end if
               ! Update face pressure and density in Y
               rho_f(1)=0.0_WP; vol_r=sum(vf%Gvol(:,0,:,i,j  ,k)+vf%Lvol(:,0,:,i,j  ,k))
               if (vol_r.gt.0.0_WP.and.this%mask(i,j  ,k).eq.0) &
                    rho_f(1)=(sum(vf%Gvol(:,0,:,i,j  ,k))*this%Grho(i,j  ,k)&
                             +sum(vf%Lvol(:,0,:,i,j  ,k))*this%Lrho(i,j  ,k))
               rho_f(0)=0.0_WP; vol_l=sum(vf%Gvol(:,1,:,i,j-1,k)+vf%Lvol(:,1,:,i,j-1,k))
               if (vol_l.gt.0.0_WP.and.this%mask(i,j+1,k).eq.0) &
                    rho_f(0)=(sum(vf%Gvol(:,1,:,i,j-1,k))*this%Grho(i,j-1,k)&
                             +sum(vf%Lvol(:,1,:,i,j-1,k))*this%Lrho(i,j-1,k))
               if (sum(rho_f).gt.0.0_WP) then
                  DP_V(i,j,k)=2.0_WP*sum(this%itpi_y(:,i,j,k)*DP(i,j-1:j,k))&
                       -sum(this%itpi_y(:,i,j,k)*rho_f*DP(i,j-1:j,k)) / &
                       (sum(this%itpi_y(:,i,j,k)*rho_f) + tiny(1.0_WP))
                  this%V(i,j,k)=this%V(i,j,k)-dt*(&
                       sum(this%divv_y(:,i,j,k)*DP(i,j-1:j,k))&
                       -this%DPjy(i,j,k))/this%rho_V(i,j,k)
               end if
               ! Update face pressure and density in Z
               rho_f(1)=0.0_WP; vol_r=sum(vf%Gvol(:,:,0,i,j,k  )+vf%Lvol(:,:,0,i,j,k  ))
               if (vol_r.gt.0.0_WP.and.this%mask(i,j,k  ).eq.0) &
                    rho_f(1)=(sum(vf%Gvol(:,:,0,i,j,k  ))*this%Grho(i,j,k  )&
                             +sum(vf%Lvol(:,:,0,i,j,k  ))*this%Lrho(i,j,k  ))
               rho_f(0)=0.0_WP; vol_l=sum(vf%Gvol(:,:,1,i,j,k-1)+vf%Lvol(:,:,1,i,j,k-1))
               if (vol_l.gt.0.0_WP.and.this%mask(i,j,k+1).eq.0) &
                    rho_f(0)=(sum(vf%Gvol(:,:,1,i,j,k-1))*this%Grho(i,j,k-1)&
                             +sum(vf%Lvol(:,:,1,i,j,k-1))*this%Lrho(i,j,k-1))
               if (sum(rho_f).gt.0.0_WP) then
                  DP_W(i,j,k)=2.0_WP*sum(this%itpi_z(:,i,j,k)*DP(i,j,k-1:k))&
                       -sum(this%itpi_z(:,i,j,k)*rho_f*DP(i,j,k-1:k)) / &
                       (sum(this%itpi_z(:,i,j,k)*rho_f) + tiny(1.0_WP))
                  this%W(i,j,k)=this%W(i,j,k)-dt*(&
                       sum(this%divw_z(:,i,j,k)*DP(i,j,k-1:k))&
                       -this%DPjz(i,j,k))/this%rho_W(i,j,k)
               end if
            end do
         end do
      end do
      
      ! BCs
      bc_scope = 'velocity'
      call this%apply_bcond(dt,bc_scope)

     ! Correct cell-centered quantities
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! No need to calculate terms inside of wall cell or masked BC
              if (this%mask(i,j,k).gt.0) cycle

              ! r and l are with respect to the cell here, not a face
              rho_l=sum(vf%Gvol(0,:,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(0,:,:,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(1,:,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(1,:,:,i,j,k))*this%Lrho(i,j,k)
              jumpx(0)=rho_l*this%DPjx(i  ,j,k)/(this%rho_U(i  ,j,k)*this%cfg%vol(i,j,k))
              jumpx(1)=rho_r*this%DPjx(i+1,j,k)/(this%rho_U(i+1,j,k)*this%cfg%vol(i,j,k))
              rho_l=sum(vf%Gvol(:,0,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,0,:,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(:,1,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,1,:,i,j,k))*this%Lrho(i,j,k)
              jumpy(0)=rho_l*this%DPjy(i,j  ,k)/(this%rho_V(i,j  ,k)*this%cfg%vol(i,j,k))
              jumpy(1)=rho_r*this%DPjy(i,j+1,k)/(this%rho_V(i,j+1,k)*this%cfg%vol(i,j,k))
              rho_l=sum(vf%Gvol(:,:,0,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,:,0,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(:,:,1,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,:,1,i,j,k))*this%Lrho(i,j,k)
              jumpz(0)=rho_l*this%DPjz(i,j,k  )/(this%rho_W(i,j,k  )*this%cfg%vol(i,j,k))
              jumpz(1)=rho_r*this%DPjz(i,j,k+1)/(this%rho_W(i,j,k+1)*this%cfg%vol(i,j,k))

              ! Update momentum in X
              this%rhoUi(i,j,k)=this%rhoUi(i,j,k)-dt*((DP_U(i+1,j,k)-DP_U(i,j,k))*this%cfg%dxi(i)-sum(jumpx))
              ! Update momentum in Y
              this%rhoVi(i,j,k)=this%rhoVi(i,j,k)-dt*((DP_V(i,j+1,k)-DP_V(i,j,k))*this%cfg%dyi(j)-sum(jumpy))
              ! Update momentum in Z
              this%rhoWi(i,j,k)=this%rhoWi(i,j,k)-dt*((DP_W(i,j,k+1)-DP_W(i,j,k))*this%cfg%dzi(k)-sum(jumpz))

	      ! Update gas total energy
              this%GrhoE(i,j,k) = this%GrhoE(i,j,k) - dt*( this%Grho(i,j,k)/this%RHO(i,j,k)*&
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)*DP_U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)*DP_V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)*DP_W(i,j,k:k+1)) &
                   - sum(this%DPjx(i:i+1,j,k)*this%U(i:i+1,j,k)) &
                   - sum(this%DPjy(i,j:j+1,k)*this%V(i,j:j+1,k)) &
                   - sum(this%DPjz(i,j,k:k+1)*this%W(i,j,k:k+1)) ) + &
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)) )*(DP(i,j,k)* &
                   (real(ceiling(1.0_WP-vf%VF(i,j,k)),WP)-this%Grho(i,j,k)/this%RHO(i,j,k)) &
                   - (       vf%VF(i,j,k))*this%dHpjump(i,j,k) ) )

              ! Update liquid total energy
              this%LrhoE(i,j,k) = this%LrhoE(i,j,k) - dt*( this%Lrho(i,j,k)/this%RHO(i,j,k)*&
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)*DP_U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)*DP_V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)*DP_W(i,j,k:k+1)) &
                   - sum(this%DPjx(i:i+1,j,k)*this%U(i:i+1,j,k)) &
                   - sum(this%DPjy(i,j:j+1,k)*this%V(i,j:j+1,k)) &
                   - sum(this%DPjz(i,j,k:k+1)*this%W(i,j,k:k+1)) ) + &
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)) )*(DP(i,j,k)* &
                   (real(ceiling(       vf%VF(i,j,k)),WP)-this%Lrho(i,j,k)/this%RHO(i,j,k)) &
                   + (1.0_WP-vf%VF(i,j,k))*this%dHpjump(i,j,k) ) )

           end do
        end do
     end do

     ! Transfer current jump to old
     this%dHpjump = this%Hpjump

     ! BCs
     bc_scope = 'momentum'
     call this%apply_bcond(dt,bc_scope)
     bc_scope = 'energy'
     call this%apply_bcond(dt,bc_scope)

     ! Calculate other quantities
     this%Ui=this%rhoUi/this%RHO
     this%Vi=this%rhoVi/this%RHO
     this%Wi=this%rhoWi/this%RHO
     
     ! End use of temporary memory
     nullify(DP_U,DP_V,DP_W)

   end subroutine pressureproj_correct

   subroutine pressure_relax(this,vf,matmod)
     use vfs_class,  only: vfs, VFlo, VFhi
     use matm_class, only: matm
     implicit none
     class(mast), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     class(matm), intent(inout) :: matmod
     real(WP) :: KE
     logical  :: Gflag, Lflag
     integer  :: i,j,k

     ! Loop through domain and relax multiphase cells
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              ! Skip wall cells and masked BCs
              if (this%mask(i,j,k).gt.0) cycle
              ! Determine if multiphase cell
              if ((vf%VF(i,j,k).ge.VFlo).and.(vf%VF(i,j,k).le.VFhi)) then
                 if (.true.) then !(mult_iso) then
                    ! Mechanical relaxation
                    call this%pressure_relax_one(vf,matmod,i,j,k)
                 else
                    ! Thermal relaxation
                    ! call pressure_thermal_relax_one(i,j,k)
                 end if
                 ! Refresh temperature and temperature-dependent variables
                 ! call therm_refresh_one(i,j,k)
              end if
           end do
        end do
     end do

     ! Update phase interface to match VF field, etc.
     ! (Mimics the end of vf%advance())
     call vf%remove_flotsams()
     call vf%advect_interface(0.0_WP,this%U,this%V,this%W)
     call vf%build_interface()
     !call vf%update_band()
     call vf%remove_sheets()
     call vf%polygonalize_interface()
     !call vf%distance_from_polygon()
     call vf%subcell_vol()
     call vf%get_curvature()
     ! * Band and distance information should not be needed
     ! * Interface information won't be used before interface changes
     !   Just VOF and barycenter values
     ! * Relaxation is assumed to be instantaneous, so advect_interface happens over dt=0

     ! Correct bad energy values if permitted
     if (this%energy_fix) then
        Gflag = .false.; Lflag = .false.
        do k=this%cfg%kmino_,this%cfg%kmaxo_
           do j=this%cfg%jmino_,this%cfg%jmaxo_
              do i=this%cfg%imino_,this%cfg%imaxo_
                 ! Skip wall cells and masked BCs
                 if (this%mask(i,j,k).gt.0) cycle
                 ! Calculate kinematic KE
                 KE = 0.5_WP*(this%Ui(i,j,k)**2+this%Vi(i,j,k)**2+this%Wi(i,j,k)**2)
                 ! Alter energy if needed
                 call matmod%fix_energy(vf%VF(i,j,k),KE,this%P(i,j,k),this%sigma*vf%curv(i,j,k),i,j,k,Gflag,Lflag)
              end do
           end do
        end do

        if (Gflag.or.Lflag) then
           print*,'Energy fixed during timestep: Gas ',Gflag,'Liquid ',Lflag
        end if
     end if
                 
     return
   end subroutine pressure_relax

   subroutine pressure_relax_one(this,vf,matmod,i,j,k)
     use vfs_class, only: vfs, VFlo, VFhi
     use matm_class, only: matm
     implicit none
     class(mast), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     class(matm), intent(inout) :: matmod
     integer,     intent(in)    :: i,j,k
     real(WP) :: KE,Pint,buf_VF,detmnt
     real(WP) :: Peq,dVF,oVF,my_pjump
     real(WP) :: rGrhoE,rLrhoE,rGrho,rLrho,rGIE,rLIE
     real(WP), dimension(4) :: VF_terms
     real(WP), dimension(3) :: quad_terms
     real(WP), dimension(2) :: Pint_terms
     real(WP), parameter :: lelim = 0.1_WP
     real(WP), parameter :: gelim = 1.0_WP-lelim
     real(WP), parameter :: Plo = 1e-8_WP

     ! Calculate kinetic energy without the density
     KE  = 0.5_WP*(this%Ui(i,j,k)**2+this%Vi(i,j,k)**2+this%Wi(i,j,k)**2)

     ! Store old VF value
     oVF = vf%VF(i,j,k)
     ! Store pressure jump
     my_pjump = this%sigma*vf%curv(i,j,k)

     ! Get coefficients for quadratic equation
     call matmod%EOS_relax_quad(i,j,k,vf%VF(i,j,k),this%srcVF(i,j,k),my_pjump,&
          this%LrhoSS2(i,j,k),this%GrhoSS2(i,j,k),      VF_terms,quad_terms,Pint_terms)
     
     ! Check if solution can exist
     detmnt = quad_terms(2)**2-4.0_WP*quad_terms(1)*quad_terms(3)
     if (detmnt.lt.0.0_WP) then
        if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
             j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
             k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
           print*,'Ptherm relax could not be performed',i,j,k,'oVF',oVF, &
                'detmnt',detmnt,'Pliq',matmod%EOS_liquid(i,j,k,'p'),'Pgas',matmod%EOS_gas(i,j,k,'p')
        end if
        return
     end if
     
     ! Solve quadratic equation for equilibrium pressure
     Peq = (-quad_terms(2) + sqrt(detmnt))/(2.0_WP*quad_terms(1))
     
     ! Use equilibrium pressure to calculate the new VF,
     ! and check if it is not NaN. Let bounds be handled later.
     buf_VF = (VF_terms(1)*Peq+VF_terms(2))/(VF_terms(3)*Peq+VF_terms(4))
     if (buf_VF.ne.buf_VF.or.&
          (Peq+matmod%Pref_g-my_pjump.lt.Plo.and.oVF.le.gelim).or.&
          (Peq+matmod%Pref_l.lt.Plo.and.oVF.ge.lelim)) then
        if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
             j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
             k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
           print*,'P relax could not be performed',i,j,k,'oVF',oVF,'VF*',buf_VF, &
                'Peq',Peq,'Pliq',matmod%EOS_liquid(i,j,k,'p'),'Pgas',matmod%EOS_gas(i,j,k,'p')
        end if
        return
     end if

     ! Calculate how much VF changes from unrelaxed state
     dVF = buf_VF - (oVF-this%srcVF(i,j,k))
     ! Calculate interface pressure according to temporal interpolation and acoustic impedances
     Pint = Pint_terms(1)*Peq+Pint_terms(2)
     ! Use the VF increment and Pint to calculate the change of the internal energy
     ! Store "relaxed" values with prefix r
     if (buf_VF.le.VFhi.and.buf_VF.ge.VFlo) then
        ! Calculate future energy unless VF is out of bounds
        rLrhoE = (oVF*this%LrhoE(i,j,k) - Pint*dVF)/buf_VF
        rGrhoE = ((1.0_WP-oVF)*this%GrhoE(i,j,k) + Pint*dVF)/(1.0_WP-buf_VF)
        ! Adjust phase quantities now that the volume fraction has changed
        rLrho  = oVF*this%Lrho(i,j,k)/buf_VF
        rGrho  = (1.0_WP-oVF)*this%Grho(i,j,k)/(1.0_WP-buf_VF)
        ! Calculate relaxed internal energy of each phase
        rLIE   = rLrhoE-rLrho*KE
        rGIE   = rGrhoE-rGrho*KE
     end if

     ! Go through cases according to the values found
     ! Check for VF bounds, then for bad values
     if (buf_VF.gt.VFhi) then
        ! If equilibrium pressure or relaxed energy is negative, let liquid phase take gas energy
        vf%VF(i,j,k) = 1.0_WP
        this%LrhoE(i,j,k) = oVF*this%LrhoE(i,j,k) + (1.0_WP-oVF)*this%GrhoE(i,j,k)
        this%GrhoE(i,j,k) = 0.0_WP
        this%Lrho (i,j,k) = oVF*this%Lrho(i,j,k)
        this%Grho (i,j,k) = 0.0_WP
     else if (buf_VF.lt.VFlo) then
        vf%VF(i,j,k) = 0.0_WP
        this%GrhoE(i,j,k) = (1.0_WP-oVF)*this%GrhoE(i,j,k) + oVF*this%LrhoE(i,j,k)
        this%LrhoE(i,j,k) = 0.0_WP
        this%Grho (i,j,k) = (1.0_WP-oVF)*this%Grho(i,j,k)
        this%Lrho (i,j,k) = 0.0_WP
     else if (rGIE.lt.Plo.and.oVF.gt.gelim) then
        ! If equilibrium energy is negative, leave other values unchanged
        if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
             j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
             k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
           print*,'Gas phase eliminated in P relax',i,j,k,'oVF',oVF,'VF*',buf_VF, &
                'Peq-pjump',Peq-my_pjump,'Pliq',matmod%EOS_liquid(i,j,k,'p'), &
                'Grhoe0',this%GrhoE(i,j,k)-this%Grho(i,j,k)*KE,'Grho0',this%Grho(i,j,k)
        end if
        vf%VF(i,j,k) = 1.0_WP
        this%LrhoE(i,j,k) = this%LrhoE(i,j,k)
        this%GrhoE(i,j,k) = 0.0_WP
        this%Lrho (i,j,k) = this%Lrho(i,j,k)
        this%Grho (i,j,k) = 0.0_WP
     else if (rLIE.lt.Plo.and.oVF.lt.lelim) then
        ! If equilibrium pressure is negative, leave other values unchanged
        if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
             j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
             k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
           print*,'Liq phase eliminated in P relax',i,j,k,'oVF',oVF,'VF*',buf_VF, &
                'Peq',Peq,'Pgas',matmod%EOS_gas(i,j,k,'p'), &
                'Lrhoe0',this%LrhoE(i,j,k)-this%Lrho(i,j,k)*KE,'Lrho0',this%Lrho(i,j,k)
        end if
        vf%VF(i,j,k) = 0.0_WP
        this%GrhoE(i,j,k) = this%GrhoE(i,j,k)
        this%LrhoE(i,j,k) = 0.0_WP
        this%Grho (i,j,k) = this%Grho(i,j,k)
        this%Lrho (i,j,k) = 0.0_WP
     else if (Peq+matmod%Pref_g-my_pjump.lt.Plo.and.oVF.gt.gelim) then
        ! If equilibrium pressure is negative, let liquid phase take the energy
        if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
             j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
             k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
           print*,'Gas phase eliminated in P relax',i,j,k,'oVF',oVF,'VF*',buf_VF, &
                'Peq-pjump',Peq-my_pjump,'Pliq',matmod%EOS_liquid(i,j,k,'p'), &
                'Grhoe0',this%GrhoE(i,j,k)-this%Grho(i,j,k)*KE,'Grho0',this%Grho(i,j,k)
        end if
        vf%VF(i,j,k) = 1.0_WP
        this%LrhoE(i,j,k) = oVF*this%LrhoE(i,j,k) + max((1.0_WP-oVF)*this%GrhoE(i,j,k),0.0_WP)
        this%GrhoE(i,j,k) = 0.0_WP
        this%Lrho (i,j,k) = oVF*this%Lrho(i,j,k) ! + (1.0_WP-oVF)*Grho(i,j,k)
        this%Grho (i,j,k) = 0.0_WP
     else if (Peq+matmod%Pref_l.lt.Plo.and.oVF.lt.lelim) then
        ! If equilibrium pressure is negative, let gas phase take energy
        if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
             j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
             k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
           print*,'Liq phase eliminated in P relax',i,j,k,'oVF',oVF,'VF*',buf_VF, &
                'Peq',Peq,'Pgas',matmod%EOS_gas(i,j,k,'p'), &
                'Lrhoe0',this%LrhoE(i,j,k)-this%Lrho(i,j,k)*KE,'Lrho0',this%Lrho(i,j,k)
        end if
        vf%VF(i,j,k) = 0.0_WP
        this%GrhoE(i,j,k) = (1.0_WP-oVF)*this%GrhoE(i,j,k) + max(oVF*this%LrhoE(i,j,k),0.0_WP)
        this%LrhoE(i,j,k) = 0.0_WP
        this%Grho (i,j,k) = (1.0_WP-oVF)*this%Grho(i,j,k) ! + oVF*Lrho(i,j,k)
        this%Lrho (i,j,k) = 0.0_WP
     else
        ! Check for errant energy values
        if (min(rGIE,rLIE,Peq+matmod%Pref_g-my_pjump,Peq+matmod%Pref_l).lt.Plo) then
           if (i.ge.this%cfg%imin_.and.i.le.this%cfg%imax_.and. &
                j.ge.this%cfg%jmin_.and.j.le.this%cfg%jmax_.and. &
                k.ge.this%cfg%kmin_.and.k.le.this%cfg%kmax_) then
              print*,'P relax could not be performed',i,j,k,'oVF',oVF,'VF*',buf_VF, &
                   'Peq+Pref_l',Peq+matmod%Pref_l,'Peq+Pref_g-pjump',Peq+matmod%Pref_g-my_pjump, &
                   'rGrhoe',rGIE,'rLrhoe',rLIE
           end if
           return
        end if

        ! Set VF, rhoE, rho to relaxed values (made it past all the conditionals)
        vf%VF(i,j,k) = buf_VF
        this%LrhoE(i,j,k) = rLrhoE
        this%GrhoE(i,j,k) = rGrhoE
        this%Lrho (i,j,k) = rLrho
        this%Grho (i,j,k) = rGrho
     end if
     
     return
   end subroutine pressure_relax_one

   !> Combine phase pressures and bulk moduli according to mixture rules, add additional terms to PA
   subroutine harmonize_advpressure_bulkmod(this,vf,matmod)
     use vfs_class, only : vfs,VFlo,VFhi
     use matm_class, only: matm
     implicit none
     class(mast), intent(inout) :: this
     class(vfs), intent(in) :: vf
     class(matm), intent(in) :: matmod
     real(WP) :: rhoc2_lI,rhoc2_gI,denom
     integer :: i,j,k

     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
              ! Intermediate variables for mechanical equilibrium projection
              rhoc2_lI = this%LrhoSS2(i,j,k); rhoc2_gI = this%GrhoSS2(i,j,k)
              if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) call matmod%bulkmod_intf(i,j,k,rhoc2_lI,rhoc2_gI)
              denom = 1.0_WP/ &
                   ((       vf%VF(i,j,k))/(rhoc2_lI+tiny(1.0_WP)) &
                   +(1.0_WP-vf%VF(i,j,k))/(rhoc2_gI+tiny(1.0_WP)))
              ! Mixture bulk modulus
              this%RHOSS2(i,j,k) = (vf%VF(i,j,k) *(this%LrhoSS2(i,j,k)/(rhoc2_lI+tiny(1.0_WP))) &
                           +(1.0_WP-vf%VF(i,j,k))*(this%GrhoSS2(i,j,k)/(rhoc2_gI+tiny(1.0_WP))))*denom
              ! Mixture advected pressure with additional jump terms
              this%PA(i,j,k) = ((vf%VF(i,j,k) *(this%LP(i,j,k)/(rhoc2_lI+tiny(1.0_WP))) &
                   +(1.0_WP-vf%VF(i,j,k))*(this%GP(i,j,k)/(rhoc2_gI+tiny(1.0_WP)))) &
                   + vf%VF(i,j,k)*(1.0_WP-vf%VF(i,j,k))*(rhoc2_lI-rhoc2_gI)/&
                   (rhoc2_lI*rhoc2_gI+tiny(1.0_WP))*this%sigma*vf%curv(i,j,k) )*denom
           end do
        end do
     end do

   end subroutine harmonize_advpressure_bulkmod

   !> Add surface tension jump term using CSF
   subroutine update_surface_tension_jump(this,vf,contact_model)
      use messager,  only: die
      use vfs_class, only: vfs
      implicit none
      class(mast), intent(inout) :: this
      class(vfs), intent(in) :: vf
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

      ! Store pressure jump at cell centers
      this%Hpjump = this%sigma*vf%curv
      ! Get increment, while considering different phase values during iteration
      this%dHpjump = this%Hpjump - this%dHpjump*real(ceiling(1.0_WP-vf%VF)*ceiling(vf%VF),WP)
      
    end subroutine update_surface_tension_jump


    subroutine interp_pressure_density(this,vf)
      use vfs_class, only : vfs
      implicit none
      class(mast), intent(inout) :: this
      class(vfs),  intent(inout) :: vf
      integer :: i,j,k
      real(WP), dimension(0:1) :: rho_f
      real(WP) :: vol_l,vol_r
      ! Initialize
      this%P_U  =0.0_WP; this%P_V  =0.0_WP; this%P_W  =0.0_WP
      this%rho_U=1.0_WP; this%rho_V=1.0_WP; this%rho_W=1.0_WP

      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1

               ! Update face pressure and density in X
               rho_f(1)=0.0_WP; vol_r=sum(vf%Gvol(0,:,:,i  ,j,k)+vf%Lvol(0,:,:,i  ,j,k))
               if (vol_r.gt.0.0_WP) rho_f(1)=(sum(vf%Gvol(0,:,:,i  ,j,k))*this%Grho(i  ,j,k)&
                                             +sum(vf%Lvol(0,:,:,i  ,j,k))*this%Lrho(i  ,j,k))
               rho_f(0)=0.0_WP; vol_l=sum(vf%Gvol(1,:,:,i-1,j,k)+vf%Lvol(1,:,:,i-1,j,k))
               if (vol_l.gt.0.0_WP) rho_f(0)=(sum(vf%Gvol(1,:,:,i-1,j,k))*this%Grho(i-1,j,k)&
                                             +sum(vf%Lvol(1,:,:,i-1,j,k))*this%Lrho(i-1,j,k))
               if (sum(rho_f).gt.0.0_WP) then
                  this%P_U(i,j,k)=2.0_WP*sum(this%itpi_x(:,i,j,k)*this%P(i-1:i,j,k))&
                       -sum(this%itpi_x(:,i,j,k)*rho_f*this%P(i-1:i,j,k)) / &
                       (sum(this%itpi_x(:,i,j,k)*rho_f) + tiny(1.0_WP))
                  this%rho_U(i,j,k)=sum(rho_f)/(vol_l+vol_r)
               end if
               ! Update face pressure and density in Y
               rho_f(1)=0.0_WP; vol_r=sum(vf%Gvol(:,0,:,i,j  ,k)+vf%Lvol(:,0,:,i,j  ,k))
               if (vol_r.gt.0.0_WP) rho_f(1)=(sum(vf%Gvol(:,0,:,i,j  ,k))*this%Grho(i,j  ,k)&
                                             +sum(vf%Lvol(:,0,:,i,j  ,k))*this%Lrho(i,j  ,k))
               rho_f(0)=0.0_WP; vol_l=sum(vf%Gvol(:,1,:,i,j-1,k)+vf%Lvol(:,1,:,i,j-1,k))
               if (vol_l.gt.0.0_WP) rho_f(0)=(sum(vf%Gvol(:,1,:,i,j-1,k))*this%Grho(i,j-1,k)&
                                             +sum(vf%Lvol(:,1,:,i,j-1,k))*this%Lrho(i,j-1,k))
               if (sum(rho_f).gt.0.0_WP) then
                  this%P_V(i,j,k)=2.0_WP*sum(this%itpi_y(:,i,j,k)*this%P(i,j-1:j,k))&
                       -sum(this%itpi_y(:,i,j,k)*rho_f*this%P(i,j-1:j,k)) / &
                       (sum(this%itpi_y(:,i,j,k)*rho_f) + tiny(1.0_WP))
                  this%rho_V(i,j,k)=sum(rho_f)/(vol_l+vol_r)
               end if
               ! Update face pressure and density in Z
               rho_f(1)=0.0_WP; vol_r=sum(vf%Gvol(:,:,0,i,j,k  )+vf%Lvol(:,:,0,i,j,k  ))
               if (vol_r.gt.0.0_WP) rho_f(1)=(sum(vf%Gvol(:,:,0,i,j,k  ))*this%Grho(i,j,k  )&
                                             +sum(vf%Lvol(:,:,0,i,j,k  ))*this%Lrho(i,j,k  ))
               rho_f(0)=0.0_WP; vol_l=sum(vf%Gvol(:,:,1,i,j,k-1)+vf%Lvol(:,:,1,i,j,k-1))
               if (vol_l.gt.0.0_WP) rho_f(0)=(sum(vf%Gvol(:,:,1,i,j,k-1))*this%Grho(i,j,k-1)&
                                             +sum(vf%Lvol(:,:,1,i,j,k-1))*this%Lrho(i,j,k-1))
               if (sum(rho_f).gt.0.0_WP) then
                  this%P_W(i,j,k)=2.0_WP*sum(this%itpi_z(:,i,j,k)*this%P(i,j,k-1:k))&
                       -sum(this%itpi_z(:,i,j,k)*rho_f*this%P(i,j,k-1:k)) / &
                       (sum(this%itpi_z(:,i,j,k)*rho_f) + tiny(1.0_WP))
                  this%rho_W(i,j,k)=sum(rho_f)/(vol_l+vol_r)
               end if
               
            end do
         end do
      end do

    end subroutine interp_pressure_density

   !> Calculate the interpolated velocity, which is the velocity at the face
   subroutine interp_vel_basic(this,vf,Ui,Vi,Wi,U,V,W)
     use vfs_class, only : vfs
     implicit none
     class(mast), intent(in) :: this
     class(vfs),  intent(in) :: vf
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Ui
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Vi
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Wi
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: U
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: V
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: W
     integer  :: i,j,k
     real(WP) :: vol_r,vol_l,rho_r,rho_l


     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
              ! X face
              vol_r=sum(vf%Gvol(0,:,:,i  ,j,k)+vf%Lvol(0,:,:,i  ,j,k))
              vol_l=sum(vf%Gvol(1,:,:,i-1,j,k)+vf%Lvol(1,:,:,i-1,j,k))
              rho_r=sum(vf%Gvol(0,:,:,i  ,j,k))*this%Grho(i  ,j,k)+sum(vf%Lvol(0,:,:,i  ,j,k))*this%Lrho(i  ,j,k)
              rho_l=sum(vf%Gvol(1,:,:,i-1,j,k))*this%Grho(i-1,j,k)+sum(vf%Lvol(1,:,:,i-1,j,k))*this%Lrho(i-1,j,k)
              if (min(vol_l,vol_r).gt.0.0_WP) then
                 if (this%umask(i,j,k).ne.2) U(i,j,k)=(rho_l*Ui(i-1,j,k)+rho_r*Ui(i,j,k))/(rho_l+rho_r)
              else
                 U(i,j,k)=0.0_WP
              end if
              ! Y face
              vol_r=sum(vf%Gvol(:,0,:,i,j  ,k)+vf%Lvol(:,0,:,i,j  ,k))
              vol_l=sum(vf%Gvol(:,1,:,i,j-1,k)+vf%Lvol(:,1,:,i,j-1,k))
              rho_r=sum(vf%Gvol(:,0,:,i,j  ,k))*this%Grho(i,j  ,k)+sum(vf%Lvol(:,0,:,i,j  ,k))*this%Lrho(i,j  ,k)
              rho_l=sum(vf%Gvol(:,1,:,i,j-1,k))*this%Grho(i,j-1,k)+sum(vf%Lvol(:,1,:,i,j-1,k))*this%Lrho(i,j-1,k)
              if (min(vol_l,vol_r).gt.0.0_WP) then
                 if (this%vmask(i,j,k).ne.2) V(i,j,k)=(rho_l*Vi(i,j-1,k)+rho_r*Vi(i,j,k))/(rho_l+rho_r)
              else
                 V(i,j,k)=0.0_WP
              end if
              ! Z face
              vol_r=sum(vf%Gvol(:,:,0,i,j,k  )+vf%Lvol(:,:,0,i,j,k  ))
              vol_l=sum(vf%Gvol(:,:,1,i,j,k-1)+vf%Lvol(:,:,1,i,j,k-1))
              rho_r=sum(vf%Gvol(:,:,0,i,j,k  ))*this%Grho(i,j,k  )+sum(vf%Lvol(:,:,0,i,j,k  ))*this%Lrho(i,j,k  )
              rho_l=sum(vf%Gvol(:,:,1,i,j,k-1))*this%Grho(i,j,k-1)+sum(vf%Lvol(:,:,1,i,j,k-1))*this%Lrho(i,j,k-1)
              if (min(vol_l,vol_r).gt.0.0_WP) then
                 if (this%wmask(i,j,k).ne.2) W(i,j,k)=(rho_l*Wi(i,j,k-1)+rho_r*Wi(i,j,k))/(rho_l+rho_r)
              else
                 W(i,j,k)=0.0_WP
              end if
           end do
        end do
     end do

     ! Boundary conditions not included here (except through masks)

   end subroutine interp_vel_basic

   !> Calculate the interpolated velocity, intended for the pressure solve
   subroutine interp_vel_full(this,vf,dt,Ui,Vi,Wi,termU,termV,termW,U,V,W)
     use vfs_class, only : vfs
     implicit none
     class(mast), intent(in) :: this
     class(vfs),  intent(in) :: vf
     real(WP), intent(in)  :: dt
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Ui
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Vi
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Wi
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: termU
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: termV
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: termW
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: U
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: V
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: W
     integer  :: i,j,k
     real(WP) :: vol_r,vol_l,rho_r,rho_l
     real(WP) :: subtract_l,subtract_r,add
     real(WP) :: u_l,u_r


     do k=this%cfg%kmin_,this%cfg%kmax_+1
        do j=this%cfg%jmin_,this%cfg%jmax_+1
           do i=this%cfg%imin_,this%cfg%imax_+1
              ! X face
              vol_r=sum(vf%Gvol(0,:,:,i  ,j,k)+vf%Lvol(0,:,:,i  ,j,k))
              vol_l=sum(vf%Gvol(1,:,:,i-1,j,k)+vf%Lvol(1,:,:,i-1,j,k))
              rho_r=sum(vf%Gvol(0,:,:,i  ,j,k))*this%Grho(i  ,j,k)+sum(vf%Lvol(0,:,:,i  ,j,k))*this%Lrho(i  ,j,k)
              rho_l=sum(vf%Gvol(1,:,:,i-1,j,k))*this%Grho(i-1,j,k)+sum(vf%Lvol(1,:,:,i-1,j,k))*this%Lrho(i-1,j,k)
              if (min(vol_l,vol_r).gt.0.0_WP) then
                 if (this%sl_face(i,j,k,1).eq.1) then
                    subtract_l = dt*termU(i-1,j,k)
                    subtract_r = dt*termU(i  ,j,k)
                    add = -dt*sum(this%divu_x(:,i,j,k)*this%P(i-1:i,j,k))/this%rho_U(i,j,k)&
                          +dt*this%Pjx(i,j,k)/this%rho_U(i,j,k)
                 else
                    subtract_l = 0.0_WP; subtract_r = 0.0_WP; add = 0.0_WP
                 end if
                 u_l = Ui(i-1,j,k)+subtract_l; u_r = Ui(i,j,k)+subtract_r
                 if (this%umask(i,j,k).ne.2) U(i,j,k)=(rho_l*u_l+rho_r*u_r)/(rho_l+rho_r) + add
              else
                 U(i,j,k)=0.0_WP
              end if
              ! Y face
              vol_r=sum(vf%Gvol(:,0,:,i,j  ,k)+vf%Lvol(:,0,:,i,j  ,k))
              vol_l=sum(vf%Gvol(:,1,:,i,j-1,k)+vf%Lvol(:,1,:,i,j-1,k))
              rho_r=sum(vf%Gvol(:,0,:,i,j  ,k))*this%Grho(i,j  ,k)+sum(vf%Lvol(:,0,:,i,j  ,k))*this%Lrho(i,j  ,k)
              rho_l=sum(vf%Gvol(:,1,:,i,j-1,k))*this%Grho(i,j-1,k)+sum(vf%Lvol(:,1,:,i,j-1,k))*this%Lrho(i,j-1,k)
              if (min(vol_l,vol_r).gt.0.0_WP) then
                 if (this%sl_face(i,j,k,2).eq.1) then
                    subtract_l = dt*termV(i,j-1,k)
                    subtract_r = dt*termV(i,j  ,k)
                    add = -dt*sum(this%divv_y(:,i,j,k)*this%P(i,j-1:j,k))/this%rho_V(i,j,k)&
                          +dt*this%Pjy(i,j,k)/this%rho_V(i,j,k)
                 else
                    subtract_l = 0.0_WP; subtract_r = 0.0_WP; add = 0.0_WP
                 end if
                 u_l = Vi(i,j-1,k)+subtract_l; u_r = Vi(i,j,k)+subtract_r
                 if (this%vmask(i,j,k).ne.2) V(i,j,k)=(rho_l*u_l+rho_r*u_r)/(rho_l+rho_r) + add
              else
                 V(i,j,k)=0.0_WP
              end if
              ! Z face
              vol_r=sum(vf%Gvol(:,:,0,i,j,k  )+vf%Lvol(:,:,0,i,j,k  ))
              vol_l=sum(vf%Gvol(:,:,1,i,j,k-1)+vf%Lvol(:,:,1,i,j,k-1))
              rho_r=sum(vf%Gvol(:,:,0,i,j,k  ))*this%Grho(i,j,k  )+sum(vf%Lvol(:,:,0,i,j,k  ))*this%Lrho(i,j,k  )
              rho_l=sum(vf%Gvol(:,:,1,i,j,k-1))*this%Grho(i,j,k-1)+sum(vf%Lvol(:,:,1,i,j,k-1))*this%Lrho(i,j,k-1)
              if (min(vol_l,vol_r).gt.0.0_WP) then
                 if (this%sl_face(i,j,k,3).eq.1) then
                    subtract_l = dt*termW(i,j,k-1)
                    subtract_r = dt*termW(i,j,k  )
                    add = -dt*sum(this%divw_z(:,i,j,k)*this%P(i,j,k-1:k))/this%rho_W(i,j,k)&
                          +dt*this%Pjz(i,j,k)/this%rho_W(i,j,k)
                 else
                    subtract_l = 0.0_WP; subtract_r = 0.0_WP; add = 0.0_WP
                 end if
                 u_l = Wi(i,j,k-1)+subtract_l; u_r = Wi(i,j,k)+subtract_r
                 if (this%wmask(i,j,k).ne.2) W(i,j,k)=(rho_l*u_l+rho_r*u_r)/(rho_l+rho_r) + add
              else
                 W(i,j,k)=0.0_WP
              end if
           end do
        end do
     end do

     ! Boundary conditions not addressed here, except through masks

   end subroutine interp_vel_full
   
   ! Make terms to remove explicit pressure gradient from the cell-centered velocity before interpolation
   subroutine terms_modified_face_velocity(this,vf,termU,termV,termW)
     use vfs_class, only : vfs
     implicit none
     class(mast), intent(in) :: this
     class(vfs),  intent(in) :: vf
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: termU
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: termV
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: termW
     integer  :: i,j,k
     real(WP) :: vol_r,vol_l,rho_r,rho_l
     real(WP), dimension(0:1) :: jumpx,jumpy,jumpz

     ! Begin at zero every time
     termU = 0.0_WP
     termV = 0.0_WP
     termW = 0.0_WP

     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_

              ! No need to calculate terms inside of wall cell or masked BCs
              if (this%mask(i,j,k).gt.0) cycle

              ! Get contributions of pressure jump to faces of current cell
              ! r and l are with respect to the cell here, not a face
              rho_l=sum(vf%Gvol(0,:,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(0,:,:,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(1,:,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(1,:,:,i,j,k))*this%Lrho(i,j,k)
              jumpx(0)=rho_l*this%Pjx(i  ,j,k)/(this%rho_U(i  ,j,k)*this%cfg%vol(i,j,k))
              jumpx(1)=rho_r*this%Pjx(i+1,j,k)/(this%rho_U(i+1,j,k)*this%cfg%vol(i,j,k))
              rho_l=sum(vf%Gvol(:,0,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,0,:,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(:,1,:,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,1,:,i,j,k))*this%Lrho(i,j,k)
              jumpy(0)=rho_l*this%Pjy(i,j  ,k)/(this%rho_V(i,j  ,k)*this%cfg%vol(i,j,k))
              jumpy(1)=rho_r*this%Pjy(i,j+1,k)/(this%rho_V(i,j+1,k)*this%cfg%vol(i,j,k))
              rho_l=sum(vf%Gvol(:,:,0,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,:,0,i,j,k))*this%Lrho(i,j,k)
              rho_r=sum(vf%Gvol(:,:,1,i,j,k))*this%Grho(i,j,k)+sum(vf%Lvol(:,:,1,i,j,k))*this%Lrho(i,j,k)
              jumpz(0)=rho_l*this%Pjz(i,j,k  )/(this%rho_W(i,j,k  )*this%cfg%vol(i,j,k))
              jumpz(1)=rho_r*this%Pjz(i,j,k+1)/(this%rho_W(i,j,k+1)*this%cfg%vol(i,j,k))

              ! Combine contributions with cell-centered pressure gradient
              termU(i,j,k)=((this%P_U(i+1,j,k)-this%P_U(i,j,k))*this%cfg%dxi(i)-sum(jumpx))/this%RHO(i,j,k)
              termV(i,j,k)=((this%P_V(i,j+1,k)-this%P_V(i,j,k))*this%cfg%dyi(j)-sum(jumpy))/this%RHO(i,j,k)
              termW(i,j,k)=((this%P_W(i,j,k+1)-this%P_W(i,j,k))*this%cfg%dzi(k)-sum(jumpz))/this%RHO(i,j,k)
           end do
        end do
     end do

     ! Communicate
     call this%cfg%sync(termU)
     call this%cfg%sync(termV)
     call this%cfg%sync(termW)


   end subroutine terms_modified_face_velocity

   subroutine reinit_phase_pressure(this,vf,matmod)
     use vfs_class,  only: vfs
     use matm_class, only: matm
     implicit none
     class(mast), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     class(matm), intent(inout) :: matmod
     integer :: i,j,k
     
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              if (vf%VF(i,j,k).gt.0.0_WP) this%LP(i,j,k) = matmod%EOS_liquid(i,j,k,'p')
              if (vf%VF(i,j,k).lt.1.0_WP) this%GP(i,j,k) = matmod%EOS_gas(i,j,k,'p')
           end do
        end do
     end do

   end subroutine reinit_phase_pressure
   
   subroutine init_phase_bulkmod(this,vf,matmod)
     use vfs_class,  only: vfs
     use matm_class, only: matm
     implicit none
     class(mast), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     class(matm), intent(inout) :: matmod
     integer :: i,j,k
     
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              if (vf%VF(i,j,k).gt.0.0_WP) this%LrhoSS2(i,j,k) = matmod%EOS_liquid(i,j,k,'M')
              if (vf%VF(i,j,k).lt.1.0_WP) this%GrhoSS2(i,j,k) = matmod%EOS_gas(i,j,k,'M')
           end do
        end do
     end do

   end subroutine init_phase_bulkmod

   subroutine update_Helmholtz_LHS(this,dt)
     class(mast), intent(inout) :: this
     real(WP), intent(in)  :: dt
     integer :: i,j,k
     ! Work on this%psolv%opr

     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! Laplacian operator
              this%psolv%opr(1,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x(-1,i+1,j,k)/this%rho_U(i+1,j,k)+&
                   &                  this%divp_x(0,i,j,k)*this%divu_x( 0,i  ,j,k)/this%rho_U(i  ,j,k)+&
                   &                  this%divp_y(1,i,j,k)*this%divv_y(-1,i,j+1,k)/this%rho_V(i,j+1,k)+&
                   &                  this%divp_y(0,i,j,k)*this%divv_y( 0,i,j  ,k)/this%rho_V(i,j  ,k)+&
                   &                  this%divp_z(1,i,j,k)*this%divw_z(-1,i,j,k+1)/this%rho_W(i,j,k+1)+&
                   &                  this%divp_z(0,i,j,k)*this%divw_z( 0,i,j,k  )/this%rho_W(i,j,k  )
              this%psolv%opr(2,i,j,k)=this%divp_x(1,i,j,k)*this%divu_x( 0,i+1,j,k)/this%rho_U(i+1,j,k)
              this%psolv%opr(3,i,j,k)=this%divp_x(0,i,j,k)*this%divu_x(-1,i  ,j,k)/this%rho_U(i  ,j,k)
              this%psolv%opr(4,i,j,k)=this%divp_y(1,i,j,k)*this%divv_y( 0,i,j+1,k)/this%rho_V(i,j+1,k)
              this%psolv%opr(5,i,j,k)=this%divp_y(0,i,j,k)*this%divv_y(-1,i,j  ,k)/this%rho_V(i,j  ,k)
              this%psolv%opr(6,i,j,k)=this%divp_z(1,i,j,k)*this%divw_z( 0,i,j,k+1)/this%rho_W(i,j,k+1)
              this%psolv%opr(7,i,j,k)=this%divp_z(0,i,j,k)*this%divw_z(-1,i,j,k  )/this%rho_W(i,j,k  )
              ! Multiply by temporal term and sign
              this%psolv%opr(:,i,j,k) = -dt**2*this%psolv%opr(:,i,j,k)
              ! Diagonal term
              this%psolv%opr(1,i,j,k) = this%psolv%opr(1,i,j,k) + 1.0_WP/this%RHOSS2(i,j,k)

           end do
        end do
     end do
     
   end subroutine update_Helmholtz_LHS

   subroutine update_Helmholtz_RHS(this,dt)
     class(mast), intent(inout) :: this
     real(WP), intent(in)  :: dt
     integer :: i,j,k
     ! Work on this%psolv%rhs
     
     do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
           do i=this%cfg%imin_,this%cfg%imax_
              ! Advected pressure minus pressure at n
              this%psolv%rhs(i,j,k) = (this%PA(i,j,k) - this%P(i,j,k))/this%RHOSS2(i,j,k)

              ! Divergence
              this%psolv%rhs(i,j,k) = this%psolv%rhs(i,j,k) - dt*&
                   ( sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1)) )

              ! Pressure jump components
              this%psolv%rhs(i,j,k) = this%psolv%rhs(i,j,k) - dt**2*&
                   ( sum(this%divp_x(:,i,j,k)*this%DPjx(i:i+1,j,k)*this%rho_U(i:i+1,j,k)) &
                   + sum(this%divp_y(:,i,j,k)*this%DPjy(i,j:j+1,k)*this%rho_V(i,j:j+1,k)) &
                   + sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)*this%rho_W(i,j,k:k+1)) )

           end do
        end do
     end do

   end subroutine update_Helmholtz_RHS
   
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
