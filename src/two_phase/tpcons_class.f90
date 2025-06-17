!> Two-phase incompressible conservative flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density in each phase.
!> Interface is represented using VOF
module tpcons_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use linsol_class,   only: linsol
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: tpcons,bcond
   
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
   type :: tpcons
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_TPCONS'  !< Solver name (default=UNNAMED_TPCONS)
      
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
      
      ! Face densities
      real(WP), dimension(:,:,:), allocatable :: RHOX     !< Density at x-face
      real(WP), dimension(:,:,:), allocatable :: RHOY     !< Density at y-face
      real(WP), dimension(:,:,:), allocatable :: RHOZ     !< Density at z-face
      real(WP), dimension(:,:,:), allocatable :: RHOXold  !< Old density at x-face
      real(WP), dimension(:,:,:), allocatable :: RHOYold  !< Old density at y-face
      real(WP), dimension(:,:,:), allocatable :: RHOZold  !< Old density at z-face
      
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
      procedure :: print=>tpcons_print                    !< Output solver to the screen
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
      procedure :: get_gradU                              !< Calculate velocity gradient tensor
      procedure :: get_ugradu                             !< Calculate (u.grad)u vector
      procedure :: get_vorticity                          !< Calculate vorticity vector
      procedure :: get_mfr                                !< Calculate outgoing MFR through each bcond
      procedure :: correct_mfr                            !< Correct for mfr mismatch to ensure global conservation
      procedure :: shift_p                                !< Shift pressure to have zero average
      procedure :: update_faceRHO                         !< Calculate face densities from provided density
      procedure :: get_viscosity                          !< Calculate viscosity fields from subcell phasic volume data in a vfs object
      procedure :: solve_implicit                         !< Solve for the velocity residuals implicitly
      procedure :: addsrc_gravity                         !< Add gravitational body force
      procedure :: add_surface_tension_jump               !< Add surface tension jump
      procedure :: add_surface_tension_jump_thin          !< Add surface tension jump - thin region version
      procedure :: add_surface_tension_jump_twoVF         !< Add surface tension jump - two decomposed VF fields
      procedure :: add_static_contact                     !< Add static contact line model to surface tension jump
   end type tpcons
   
   
contains
   
   
   !> Initialization for two-phase incompressible flow solver
   subroutine initialize(this,cfg,name)
      implicit none
      class(tpcons), intent(inout) :: this
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
      allocate(this%U   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      allocate(this%P   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%P=0.0_WP
      allocate(this%Pjx (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Pjx=0.0_WP
      allocate(this%Pjy (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Pjy=0.0_WP
      allocate(this%Pjz (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Pjz=0.0_WP
      allocate(this%dPjx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%dPjx=0.0_WP
      allocate(this%dPjy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%dPjy=0.0_WP
      allocate(this%dPjz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%dPjz=0.0_WP
      allocate(this%visc   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc   =0.0_WP
      allocate(this%visc_xy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_xy=0.0_WP
      allocate(this%visc_yz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_yz=0.0_WP
      allocate(this%visc_zx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc_zx=0.0_WP
      
      ! Mass conservation data around which to build momentum/energy conservation
      allocate(this%RHOX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOX=0.0_WP
      allocate(this%RHOY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOY=0.0_WP
      allocate(this%RHOZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOZ=0.0_WP
      allocate(this%rhoU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoU=0.0_WP
      allocate(this%rhoV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoV=0.0_WP
      allocate(this%rhoW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%rhoW=0.0_WP
      
      ! Allocate flow divergence
      allocate(this%div(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%div=0.0_WP
      
      ! Allocate old flow variables
      allocate(this%Uold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Uold=0.0_WP
      allocate(this%Vold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Vold=0.0_WP
      allocate(this%Wold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Wold=0.0_WP
      allocate(this%RHOXold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOXold=0.0_WP
      allocate(this%RHOYold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOYold=0.0_WP
      allocate(this%RHOZold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOZold=0.0_WP
      
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
      class(tpcons), intent(inout) :: this
      integer :: i,j,k,st1,st2
      real(WP), dimension(-1:0) :: itpx,itpy,itpz
      
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
      class(tpcons), intent(inout) :: this
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
      class(tpcons), intent(inout) :: this
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
         
      else
         
         ! Point to implicit solver linsol object
         this%implicit=>NULL()
         
      end if
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,face,dir,canCorrect)
      use string,         only: lowercase
      use messager,       only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(tpcons), intent(inout) :: this
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
      case default; call die('[tpcons add_bcond] Unknown bcond face - expecting x, y, or z')
      end select
      new_bc%itr=iterator(pg=this%cfg,name=new_bc%name,locator=locator,face=new_bc%face)
      select case (dir) ! Outward-oriented
      case (+1); new_bc%dir=+1
      case (-1); new_bc%dir=-1
      case ( 0); new_bc%dir= 0
      case default; call die('[tpcons add_bcond] Unknown bcond dir - expecting -1, +1, or 0')
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
      case (dirichlet) !< Dirichlet is set one face (i.e., velocity component) at the time
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
         call die('[tpcons apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(tpcons), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[tpcons get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      implicit none
      class(tpcons), intent(inout) :: this
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
               ! If needed, no penetration
               if (my_bc%type.eq.slip) then
                  select case (my_bc%face)
                  case ('x')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        this%U(i,j,k)=0.0_WP
                     end do
                  case ('y')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        this%V(i,j,k)=0.0_WP
                     end do
                  case ('z')
                     do n=1,my_bc%itr%n_
                        i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                        this%W(i,j,k)=0.0_WP
                     end do
                  end select
               end if
               
            case (convective)   ! Not implemented yet!
               
            case default
               call die('[tpcons apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! Sync full fields after all bcond
      call this%cfg%sync(this%U)
      call this%cfg%sync(this%V)
      call this%cfg%sync(this%W)
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit time derivative of momentum given rhoU, U, and P
   subroutine get_dmomdt(this,drhoUdt,drhoVdt,drhoWdt)
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,jj,kk
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      
      ! Zero out drhoUVW/dt arrays
      drhoUdt=0.0_WP; drhoVdt=0.0_WP; drhoWdt=0.0_WP
      
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
               FX(i,j,k)=-sum(this%itpu_x(:,i,j,k)*this%rhoU(i:i+1,j,k))*sum(this%itpu_x(:,i,j,k)*this%U(i:i+1,j,k)) &
               &         +this%visc   (i,j,k)*(sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))+sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))) &
               &         -this%P(i,j,k)
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               FY(i,j,k)=-sum(this%itpv_x(:,i,j,k)*this%rhoV(i-1:i,j,k))*sum(this%itpu_y(:,i,j,k)*this%U(i,j-1:j,k)) &
               &         +this%visc_xy(i,j,k)*(sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k))+sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k)))
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               FZ(i,j,k)=-sum(this%itpw_x(:,i,j,k)*this%rhoW(i-1:i,j,k))*sum(this%itpu_z(:,i,j,k)*this%U(i,j,k-1:k)) &
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
               FX(i,j,k)=-sum(this%itpu_y(:,i,j,k)*this%rhoU(i,j-1:j,k))*sum(this%itpv_x(:,i,j,k)*this%V(i-1:i,j,k)) &
               &         +this%visc_xy(i,j,k)*(sum(this%grdv_x(:,i,j,k)*this%V(i-1:i,j,k))+sum(this%grdu_y(:,i,j,k)*this%U(i,j-1:j,k)))
               ! Fluxes on y-face
               i=ii-1; j=jj-1; k=kk-1
               FY(i,j,k)=-sum(this%itpv_y(:,i,j,k)*this%rhoV(i,j:j+1,k))*sum(this%itpv_y(:,i,j,k)*this%V(i,j:j+1,k)) &
               &         +this%visc   (i,j,k)*(sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))+sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))) &
               &         -this%P(i,j,k)
               ! Fluxes on z-face
               i=ii; j=jj; k=kk
               FZ(i,j,k)=-sum(this%itpw_y(:,i,j,k)*this%rhoW(i,j-1:j,k))*sum(this%itpv_z(:,i,j,k)*this%V(i,j,k-1:k)) &
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
               FX(i,j,k)=-sum(this%itpu_z(:,i,j,k)*this%rhoU(i,j,k-1:k))*sum(this%itpw_x(:,i,j,k)*this%W(i-1:i,j,k)) &
               &         +this%visc_zx(i,j,k)*(sum(this%grdw_x(:,i,j,k)*this%W(i-1:i,j,k))+sum(this%grdu_z(:,i,j,k)*this%U(i,j,k-1:k)))
               ! Fluxes on y-face
               i=ii; j=jj; k=kk
               FY(i,j,k)=-sum(this%itpv_z(:,i,j,k)*this%rhoV(i,j,k-1:k))*sum(this%itpw_y(:,i,j,k)*this%W(i,j-1:j,k)) &
               &         +this%visc_yz(i,j,k)*(sum(this%grdw_y(:,i,j,k)*this%W(i,j-1:j,k))+sum(this%grdv_z(:,i,j,k)*this%V(i,j,k-1:k)))
               ! Fluxes on z-face
               i=ii-1; j=jj-1; k=kk-1
               FZ(i,j,k)=-sum(this%itpw_z(:,i,j,k)*this%rhoW(i,j,k:k+1))*sum(this%itpw_z(:,i,j,k)*this%W(i,j,k:k+1)) &
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
      class(tpcons), intent(inout) :: this
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
                     this%psolv%opr(this%psolv%stmap(s1+s2,0,0),i,j,k)=this%psolv%opr(this%psolv%stmap(s1+s2,0,0),i,j,k)+this%divp_x(s1,i,j,k)*this%divu_x(s2,i+s1,j,k)/this%RHOX(i+s1,j,k)
                     this%psolv%opr(this%psolv%stmap(0,s1+s2,0),i,j,k)=this%psolv%opr(this%psolv%stmap(0,s1+s2,0),i,j,k)+this%divp_y(s1,i,j,k)*this%divv_y(s2,i,j+s1,k)/this%RHOY(i,j+s1,k)
                     this%psolv%opr(this%psolv%stmap(0,0,s1+s2),i,j,k)=this%psolv%opr(this%psolv%stmap(0,0,s1+s2),i,j,k)+this%divp_z(s1,i,j,k)*this%divw_z(s2,i,j,k+s1)/this%RHOZ(i,j,k+s1)
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
   subroutine get_div(this,src)
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), optional :: src !< Mass source term
      integer :: i,j,k
      ! Calculate divergence of velocity
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div(i,j,k)=sum(this%divp_x(:,i,j,k)*this%U(i:i+1,j,k))+&
               &               sum(this%divp_y(:,i,j,k)*this%V(i,j:j+1,k))+&
               &               sum(this%divp_z(:,i,j,k)*this%W(i,j,k:k+1))
            end do
         end do
      end do
      ! If present, account for mass source
      if (present(src)) then
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%div(i,j,k)=this%div(i,j,k)-src(i,j,k)
               end do
            end do
         end do
      end if
      ! Sync it
      call this%cfg%sync(this%div)
   end subroutine get_div
   
   
   !> Add surface tension jump term using CSF
   subroutine add_surface_tension_jump(this,dt,div,vf,contact_model)
      use messager,  only: die
      use vfs_class, only: vfs
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), intent(inout) :: dt     !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      class(vfs), intent(inout) :: vf
      integer, intent(in), optional :: contact_model
      integer :: i,j,k,s1
      real(WP) :: mycurv,mysurf
      
      ! Store old jump
      this%dPjx=this%Pjx
      this%dPjy=this%Pjy
      this%dPjz=this%Pjz
      
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
            call die('[tpcons: add_surface_tension_jump] Unknown contact model!')
         end select
      end if
      
      ! Compute jump of dP
      this%dPjx=this%Pjx-this%dPjx
      this%dPjy=this%Pjy-this%dPjy
      this%dPjz=this%Pjz-this%dPjz
      
      ! Add div(Pjump) to RP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do s1=0,1
                  div(i,j,k)=div(i,j,k)+dt*this%divp_x(s1,i,j,k)*this%dPjx(i+s1,j,k)/this%RHOX(i+s1,j,k)
                  div(i,j,k)=div(i,j,k)+dt*this%divp_y(s1,i,j,k)*this%dPjy(i,j+s1,k)/this%RHOY(i,j+s1,k)
                  div(i,j,k)=div(i,j,k)+dt*this%divp_z(s1,i,j,k)*this%dPjz(i,j,k+s1)/this%RHOZ(i,j,k+s1)
               end do
            end do
         end do
      end do
      
   end subroutine add_surface_tension_jump
   
   
   !> Add surface tension jump term using CSF
   !> Account for thin regions
   subroutine add_surface_tension_jump_thin(this,dt,div,vf)
      use messager,  only: die
      use vfs_class, only: vfs
      use irl_fortran_interface
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), intent(inout) :: dt     !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,n,np1,np2,ind,s1
      real(WP) :: mycurv,mysurf
      real(WP), dimension(3) :: n1,n2,nf
      real(WP), dimension(:,:,:), allocatable :: alpha  !< Modified volume fraction field
      
      ! Store old jump
      this%dPjx=this%Pjx
      this%dPjy=this%Pjy
      this%dPjz=this%Pjz
      
      ! Prepare modified volume fraction field that accounts for thin regions
      allocate(alpha(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (vf%thin_sensor(i,j,k).eq.1.0_WP) then
                  alpha(i,j,k)=1.0_WP !< Thin liquid region
               else if (vf%thin_sensor(i,j,k).eq.2.0_WP) then
                  alpha(i,j,k)=0.0_WP !< Thin gas region
               else
                  alpha(i,j,k)=vf%VF(i,j,k) !< Resolved region
               end if
            end do
         end do
      end do
      
      ! Calculate pressure jump
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Count my planes
               np1=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(vf%interface_polygon(n,i,j,k)).gt.0) np1=np1+1
               end do
               ! X face ===========================================
               mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
               if (mysurf.gt.0.0_WP) then
                  mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
               else
                  mycurv=0.0_WP
               end if
               ! Count neighbor's planes
               np2=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i-1,j,k))
                  if (getNumberOfVertices(vf%interface_polygon(n,i-1,j,k)).gt.0) np2=np2+1
               end do
               ! Assign correct curvature
               if (np1.eq.2.and.np2.eq.0) then
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                  nf=sign(1.0_WP,alpha(i-1,j,k)-0.5_WP)*[+1.0_WP,0.0_WP,0.0_WP]
                  ind=maxloc([dot_product(n1,nf),dot_product(n2,nf)],dim=1)
                  mycurv=vf%curv2p(ind,i,j,k)
               end if
               if (np1.eq.0.and.np2.eq.2) then
                  n1=calculateNormal(vf%interface_polygon(1,i-1,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i-1,j,k))
                  nf=sign(1.0_WP,alpha(i,j,k)-0.5_WP)*[-1.0_WP,0.0_WP,0.0_WP]
                  ind=maxloc([dot_product(n1,nf),dot_product(n2,nf)],dim=1)
                  mycurv=vf%curv2p(ind,i-1,j,k)
               end if
               ! Assemble pressure jump
               this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*alpha(i-1:i,j,k))
               ! Y face ===========================================
               mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
               if (mysurf.gt.0.0_WP) then
                  mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
               else
                  mycurv=0.0_WP
               end if
               ! Count neighbor's planes
               np2=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j-1,k))
                  if (getNumberOfVertices(vf%interface_polygon(n,i,j-1,k)).gt.0) np2=np2+1
               end do
               ! Assign correct curvature
               if (np1.eq.2.and.np2.eq.0) then
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                  nf=sign(1.0_WP,alpha(i,j-1,k)-0.5_WP)*[0.0_WP,+1.0_WP,0.0_WP]
                  ind=maxloc([dot_product(n1,nf),dot_product(n2,nf)],dim=1)
                  mycurv=vf%curv2p(ind,i,j,k)
               end if
               if (np1.eq.0.and.np2.eq.2) then
                  n1=calculateNormal(vf%interface_polygon(1,i,j-1,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j-1,k))
                  nf=sign(1.0_WP,alpha(i,j,k)-0.5_WP)*[0.0_WP,-1.0_WP,0.0_WP]
                  ind=maxloc([dot_product(n1,nf),dot_product(n2,nf)],dim=1)
                  mycurv=vf%curv2p(ind,i,j-1,k)
               end if
               ! Assemble pressure jump
               this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*alpha(i,j-1:j,k))
               ! Z face ===========================================
               mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
               if (mysurf.gt.0.0_WP) then
                  mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
               else
                  mycurv=0.0_WP
               end if
               ! Count neighbor's planes
               np2=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k-1))
                  if (getNumberOfVertices(vf%interface_polygon(n,i,j,k-1)).gt.0) np2=np2+1
               end do
               ! Assign correct curvature
               if (np1.eq.2.and.np2.eq.0) then
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                  nf=sign(1.0_WP,alpha(i,j,k-1)-0.5_WP)*[0.0_WP,0.0_WP,+1.0_WP]
                  ind=maxloc([dot_product(n1,nf),dot_product(n2,nf)],dim=1)
                  mycurv=vf%curv2p(ind,i,j,k)
               end if
               if (np1.eq.0.and.np2.eq.2) then
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k-1))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k-1))
                  nf=sign(1.0_WP,alpha(i,j,k)-0.5_WP)*[0.0_WP,0.0_WP,-1.0_WP]
                  ind=maxloc([dot_product(n1,nf),dot_product(n2,nf)],dim=1)
                  mycurv=vf%curv2p(ind,i,j,k-1)
               end if
               ! Assemble pressure jump
               this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*alpha(i,j,k-1:k))
            end do
         end do
      end do
      
      ! Deallocate alpha
      deallocate(alpha)
      
      ! Compute jump of DP
      this%dPjx=this%Pjx-this%dPjx
      this%dPjy=this%Pjy-this%dPjy
      this%dPjz=this%Pjz-this%dPjz
      
      ! Add div(Pjump) to RP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do s1=0,1
                  div(i,j,k)=div(i,j,k)+dt*this%divp_x(s1,i,j,k)*this%dPjx(i+s1,j,k)/this%RHOX(i+s1,j,k)
                  div(i,j,k)=div(i,j,k)+dt*this%divp_y(s1,i,j,k)*this%dPjy(i,j+s1,k)/this%RHOY(i,j+s1,k)
                  div(i,j,k)=div(i,j,k)+dt*this%divp_z(s1,i,j,k)*this%dPjz(i,j,k+s1)/this%RHOZ(i,j,k+s1)
               end do
            end do
         end do
      end do
      
   end subroutine add_surface_tension_jump_thin
   
   
   !> Add surface tension jump term using CSF
   !> Account for thin regions by using building two VF fields
   subroutine add_surface_tension_jump_twoVF(this,dt,div,vf)
      use messager,  only: die
      use vfs_class, only: vfs
      use irl_fortran_interface
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), intent(inout) :: dt     !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,ii,jj,kk,n,np1,np2,cn,s1
      integer, dimension(2,2) :: plane_ind
      real(WP), dimension(2,2) :: VF2p,curv2p,surf2p
      real(WP), dimension(3,2) :: pos
      real(WP), dimension(4,2,2) :: normals
      real(WP), dimension(4) :: plane
      type(RectCub_type) :: cell
      type(PlanarSep_type) :: lginterface
      
      ! Store old jump
      this%dPjx=this%Pjx
      this%dPjy=this%Pjy
      this%dPjz=this%Pjz
      
      ! Allocate IRL storage
      call new(cell)
      call new(lginterface)
      call setNumberOfPlanes(lginterface,1)
      
      ! Calculate pressure jump
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face ===========================================
               np1=0; np2=0; VF2p=0.0_WP; curv2p=0.0_WP; surf2p=0.0_WP; normals=0.0_WP; plane_ind=0
               ! Get info from cells on each side of the face
               ii=i-1; jj=j; kk=k; cn=1; np1=prepare_info()
               ii=i  ; jj=j; kk=k; cn=2; np2=prepare_info()
               ! Assemble pressure jump
               if ((np1.eq.0).and.(np2.eq.0)) then
                  this%Pjx(i,j,k)=0.0_WP
               else
                  this%Pjx(i,j,k)=this%sigma*sum(this%divu_x(:,i,j,k)*get_kalpha())
               end if
               ! Y face ===========================================
               np1=0; np2=0; VF2p=0.0_WP; curv2p=0.0_WP; surf2p=0.0_WP; normals=0.0_WP; plane_ind=0
               ! Get info from cells on each side of the face
               ii=i; jj=j-1; kk=k; cn=1; np1=prepare_info()
               ii=i; jj=j  ; kk=k; cn=2; np2=prepare_info()
               ! Assemble pressure jump
               if ((np1.eq.0).and.(np2.eq.0)) then
                  this%Pjy(i,j,k)=0.0_WP
               else
                  this%Pjy(i,j,k)=this%sigma*sum(this%divv_y(:,i,j,k)*get_kalpha())
               end if
               ! Z face ===========================================
               np1=0; np2=0; VF2p=0.0_WP; curv2p=0.0_WP; surf2p=0.0_WP; normals=0.0_WP; plane_ind=0
               ! Get info from cells on each side of the face
               ii=i; jj=j; kk=k-1; cn=1; np1=prepare_info()
               ii=i; jj=j; kk=k  ; cn=2; np2=prepare_info()
               ! Assemble pressure jump
               if ((np1.eq.0).and.(np2.eq.0)) then
                  this%Pjz(i,j,k)=0.0_WP
               else
                  this%Pjz(i,j,k)=this%sigma*sum(this%divw_z(:,i,j,k)*get_kalpha())
               end if
            end do
         end do
      end do
      
      ! Compute jump of DP
      this%dPjx=this%Pjx-this%dPjx
      this%dPjy=this%Pjy-this%dPjy
      this%dPjz=this%Pjz-this%dPjz
      
      ! Add div(Pjump) to RP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do s1=0,1
                  div(i,j,k)=div(i,j,k)+dt*this%divp_x(s1,i,j,k)*this%dPjx(i+s1,j,k)/this%RHOX(i+s1,j,k)
                  div(i,j,k)=div(i,j,k)+dt*this%divp_y(s1,i,j,k)*this%dPjy(i,j+s1,k)/this%RHOY(i,j+s1,k)
                  div(i,j,k)=div(i,j,k)+dt*this%divp_z(s1,i,j,k)*this%dPjz(i,j,k+s1)/this%RHOZ(i,j,k+s1)
               end do
            end do
         end do
      end do
      
   contains
      
      function prepare_info() result(np)
         implicit none
         integer :: np
         real(WP) :: myvol
         np=0; myvol=0.0_WP
         do n=1,getNumberOfPlanes(vf%liquid_gas_interface(ii,jj,kk))
            if (getNumberOfVertices(vf%interface_polygon(n,ii,jj,kk)).eq.0) cycle
            np=np+1; plane_ind(n,cn)=1
            call construct_2pt(cell,[vf%cfg%x(ii),vf%cfg%y(jj),vf%cfg%z(kk)],[vf%cfg%x(ii+1),vf%cfg%y(jj+1),vf%cfg%z(kk+1)])
            plane=getPlane(vf%liquid_gas_interface(ii,jj,kk),n-1)
            call setPlane(lginterface,0,plane(1:3),plane(4))
            call getNormMoments(cell,lginterface,myvol)
            VF2p(n,cn)=myvol/vf%cfg%vol(ii,jj,kk)
            normals(:,n,cn)=plane
            surf2p(n,cn)=abs(calculateVolume(vf%interface_polygon(n,ii,jj,kk)))
            curv2p(n,cn)=vf%curv2p(n,ii,jj,kk)
         end do
         pos(:,cn)=[vf%cfg%xm(ii),vf%cfg%ym(jj),vf%cfg%zm(kk)]
      end function  prepare_info
      
      function get_kalpha() result(kalpha)   
         use mathtools, only: normalize
         implicit none
         real(WP), dimension(3) :: pos_temp
         real(WP), dimension(2) :: kalpha
         integer, dimension(2,2) :: plane_check
         real(WP) :: surfsum,dist
         integer :: flag
         ! Merge left cell
         if (np1.eq.2) then
            if (dot_product(normals(1:3,1,1),normals(1:3,2,1)).gt.0.0_WP) then
               surfsum=sum(surf2p(:,1))
               normals(:,1,1)=(surf2p(1,1)*normals(:,1,1)+surf2p(2,1)*normals(:,2,1))/surfsum
               normals(:,1,1)=normals(:,1,1)/(norm2(normals(1:3,1,1))+tiny(1.0_WP))
               normals(:,2,1)=0.0_WP
               VF2p(:,1)=[sum(surf2p(:,1)*VF2p(:,1))/surfsum,0.0_WP]
               curv2p(:,1)=[sum(surf2p(:,1)*curv2p(:,1))/surfsum,0.0_WP]
               surf2p(:,1)=[surfsum,0.0_WP]
               plane_ind(:,1)=[1,0]
            end if
         end if
         ! Merge right cell
         if (np2.eq.2) then
            if (dot_product(normals(1:3,1,2),normals(1:3,2,2)).gt.0.0_WP) then
               surfsum = sum(surf2p(:,2))
               normals(:,1,2)=(surf2p(1,2)*normals(:,1,2)+surf2p(2,2)*normals(:,2,2))/surfsum
               normals(:,1,2)=normals(:,1,2)/(norm2(normals(1:3,1,2))+tiny(1.0_WP))
               normals(:,2,2)=0.0_WP
               VF2p(:,2)=[sum(surf2p(:,2)*VF2p(:,2))/surfsum,0.0_WP]
               curv2p(:,2)=[sum(surf2p(:,2)*curv2p(:,2))/surfsum,0.0_WP]
               surf2p(:,2)=[surfsum,0.0_WP]
               plane_ind(:,2)=[1,0]
            end if
         end if
         plane_check=0
         do ii=1,2
            if (plane_ind(ii,1).eq.0) cycle    ! If a plane doesn't exist, move on
            flag=0
            do jj=1,2
               if (plane_ind(jj,2).eq.0) cycle ! If a plane doesn't exist, move on
               if (dot_product(normals(1:3,ii,1),normals(1:3,jj,2)).gt.0.0_WP) then
                  ! IF two plane aligns, both planes get the surface area weighted curvature
                  flag=1
                  curv2p(ii,1)=(curv2p(ii,1)*surf2p(ii,1)+curv2p(jj,2)*surf2p(jj,2))/(surf2p(ii,1)+surf2p(jj,2))
                  curv2p(jj,2)=curv2p(ii,1)
                  plane_check(ii,1)=1; plane_check(jj,2)=1
               end if
            end do
            if (flag.eq.0) then ! NO PLANE with the same orientation is found
               dist=-dot_product(normals(1:3,ii,1),pos(:,2))+normals(4,ii,1)
               pos_temp=normalize(-dist*normals(1:3,ii,1))
               ! Set first volume of cell2 based on centroid projection
               if (plane_ind(1,2).ne.1.and.plane_ind(2,2).eq.1) then
                  if (dot_product(pos_temp,normals(1:3,ii,1)).le.0.0_WP) then
                     VF2p(1,2)=1.0_WP
                  else
                     VF2p(1,2)=0.0_WP
                  end if
                  curv2p(1,2)=curv2p(ii,1)
                  plane_check(ii,1)=1; plane_check(1,2)=1
               ! Set second volume of cell2 based on centroid projection
               else if (plane_ind(2,2).ne.1.and.plane_ind(1,2).eq.1) then
                  if (dot_product(pos_temp,normals(1:3,ii,1)).le.0.0_WP) then
                     VF2p(2,2)=1.0_WP
                  else
                     VF2p(2,2)=0.0_WP
                  end if
                  curv2p(2,2)=curv2p(ii,1)
                  plane_check(ii,1)=1; plane_check(2,2)=1
               ! If both cells are empty, pick the cell that hasn't been checked
               else if (plane_ind(2,2).ne.1.and.plane_ind(1,2).ne.1) then
                  if (plane_check(1,2).eq.0) then
                     if (dot_product(pos_temp,normals(1:3,ii,1)).le.0.0_WP) then
                        VF2p(1,2)=1.0_WP
                     else
                        VF2p(1,2)=0.0_WP
                     end if
                     curv2p(1,2)=curv2p(ii,1)
                     plane_check(ii,1)=1; plane_check(1,2)=1
                  else if (plane_check(2,2).eq.0) then
                     if (dot_product(pos_temp,normals(1:3,ii,1)).le.0.0_WP) then
                        VF2p(2,2)=1.0_WP
                     else
                        VF2p(2,2)=0.0_WP
                     end if
                     curv2p(2,2)=curv2p(ii,1)
                     plane_check(ii,1)=1; plane_check(2,2)=1
                  end if
               else ! don't know what to do, do nothing
               end if
            end if
         end do
         ! Next loop over all the planes from cell2 
         do ii=1,2
            ! If a plane doesn't exist or has been already accounted for, skip
            if ((plane_ind(ii,2).eq.0).or.(plane_check(ii,2).eq.1)) cycle
            dist=-dot_product(normals(1:3,ii,2),pos(:,1))+normals(4,ii,2)
            pos_temp=normalize(-dist*normals(1:3,ii,2))
            ! Set first volume of cell1 based on centroid projection
            if (plane_ind(1,1).ne.1.and.plane_ind(2,1).eq.1) then
               if (dot_product(pos_temp,normals(1:3,ii,2)).le.0.0_WP) then
                  VF2p(1,1)=1.0_WP
               else
                  VF2p(1,1)=0.0_WP
               end if
               curv2p(1,1)=curv2p(ii,2)
               plane_check(ii,2)=1; plane_check(1,1)=1
            ! Set second volume of cell1 based on centroid projection
            else if (plane_ind(2,1).ne.1.and.plane_ind(1,1).eq.1) then
               if (dot_product(pos_temp,normals(1:3,ii,2)).le.0.0_WP) then
                  VF2p(2,1)=1.0_WP
               else
                  VF2p(2,1)=0.0_WP
               end if
               curv2p(2,1)=curv2p(ii,2)
               plane_check(ii,2)=1; plane_check(2,1)=1
            ! If both cells are empty, pick the cell that hasn't been checked
            else if (plane_ind(1,1).ne.1.and.plane_ind(2,1).ne.1) then
               if (plane_check(1,1).eq.0) then
                  if (dot_product(pos_temp,normals(1:3,ii,2)).le.0.0_WP) then
                     VF2p(1,1)=1.0_WP
                  else
                     VF2p(1,1)=0.0_WP
                  end if
                  curv2p(1,1)=curv2p(ii,2)
                  plane_check(ii,2)=1; plane_check(1,1)=1
               else if (plane_check(2,1).eq.0) then
                  if (dot_product(pos_temp,normals(1:3,ii,2)).le.0.0_WP) then
                     VF2p(2,1)=1.0_WP
                  else
                     VF2p(2,1)=0.0_WP
                  end if
                  curv2p(2,1)=curv2p(ii,2)
                  plane_check(ii,2)=1; plane_check(2,1)=1
               end if
            else ! don't know what to do, do nothing
            end if
         end do
         kalpha=[sum(VF2p(:,1)*curv2p(:,1)),sum(VF2p(:,2)*curv2p(:,2))]
      end function get_kalpha
      
   end subroutine add_surface_tension_jump_twoVF
   
   
   !> Calculate the pressure gradient based on P
   subroutine get_pgrad(this,P,Pgradx,Pgrady,Pgradz)
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: P      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradx !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgrady !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pgradz !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      Pgradx=0.0_WP; Pgrady=0.0_WP; Pgradz=0.0_WP
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
      class(tpcons), intent(inout) :: this
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
      class(tpcons), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: SR  !< Needs to be (1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      real(WP) :: div
      integer :: i,j,k
      
      ! Check SR's first dimension
	   if (size(SR,dim=1).ne.6) call die('[tpcons get_strainrate] SR should be of size (1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
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

   
   !> Calculate the (u.grad)u vector from U/V/W
   subroutine get_ugradu(this,ugradu)
      use messager, only: die
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: ugradu  !< Needs to be (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      real(WP) :: Ui,gradUx,gradUy,gradUz
      real(WP) :: Vi,gradVx,gradVy,gradVz
      real(WP) :: Wi,gradWx,gradWy,gradWz
      
      ! Check ugradu's first dimension
	   if (size(ugradu,dim=1).ne.3) call die('[tpcons get_ugradu] gradu should be of size (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
      ! Allocate off-diagonal components of the velocity gradient
	   allocate(dudy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dudz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dvdz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(dwdy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate off-diagonal components of the velocity gradient at their natural locations with an extra cell for interpolation
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
      
      ! Assemble (u.grad)u at cell center
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Calculate velocity gradient directly or by interpolation
               gradUx=sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))
               gradUy=0.25_WP*sum(dudy(i:i+1,j:j+1,k))
               gradUz=0.25_WP*sum(dudz(i:i+1,j,k:k+1))
               gradVx=0.25_WP*sum(dvdx(i:i+1,j:j+1,k))
               gradVy=sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))
               gradVz=0.25_WP*sum(dvdz(i,j:j+1,k:k+1))
               gradWx=0.25_WP*sum(dwdx(i:i+1,j,k:k+1))
               gradWy=0.25_WP*sum(dwdy(i,j:j+1,k:k+1))
               gradWz=sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))
               ! Interpolate velocity to cell center
               Ui=0.5_WP*sum(this%U(i:i+1,j,k))
               Vi=0.5_WP*sum(this%V(i,j:j+1,k))
               Wi=0.5_WP*sum(this%W(i,j,k:k+1))
               ! Assemble (u.grad)u
               ugradu(1,i,j,k)=Ui*gradUx+Vi*gradUy+Wi*gradUz
               ugradu(2,i,j,k)=Ui*gradVx+Vi*gradVy+Wi*gradVz
               ugradu(3,i,j,k)=Ui*gradWx+Vi*gradWy+Wi*gradWz
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            ugradu(:,this%cfg%imin-1,:,:)=ugradu(:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) ugradu(:,this%cfg%imax+1,:,:)=ugradu(:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            ugradu(:,:,this%cfg%jmin-1,:)=ugradu(:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) ugradu(:,:,this%cfg%jmax+1,:)=ugradu(:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            ugradu(:,:,:,this%cfg%kmin-1)=ugradu(:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) ugradu(:,:,:,this%cfg%kmax+1)=ugradu(:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
	   do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) ugradu(:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Sync it
	   call this%cfg%sync(ugradu)
      
      ! Deallocate velocity gradient storage
	   deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
      
   end subroutine get_ugradu


   !> Calculate the velocity gradient tensor from U/V/W
   !> Note that gradU(i,j)=duj/dxi
   subroutine get_gradU(this,gradU)
      use messager, only: die
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: gradU  !< Needs to be (1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      
      ! Check gradU's first two dimensions
	   if (size(gradU,dim=1).ne.3.or.size(gradU,dim=2).ne.3) call die('[tpcons get_gradU] gradU should be of size (1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
      ! Compute dudx, dvdy, and dwdz first
	   do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               gradU(1,1,i,j,k)=sum(this%grdu_x(:,i,j,k)*this%U(i:i+1,j,k))
               gradU(2,2,i,j,k)=sum(this%grdv_y(:,i,j,k)*this%V(i,j:j+1,k))
               gradU(3,3,i,j,k)=sum(this%grdw_z(:,i,j,k)*this%W(i,j,k:k+1))
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
               gradU(2,1,i,j,k)=0.25_WP*sum(dudy(i:i+1,j:j+1,k))
               gradU(3,1,i,j,k)=0.25_WP*sum(dudz(i:i+1,j,k:k+1))
               gradU(1,2,i,j,k)=0.25_WP*sum(dvdx(i:i+1,j:j+1,k))
               gradU(3,2,i,j,k)=0.25_WP*sum(dvdz(i,j:j+1,k:k+1))
               gradU(1,3,i,j,k)=0.25_WP*sum(dwdx(i:i+1,j,k:k+1))
               gradU(2,3,i,j,k)=0.25_WP*sum(dwdy(i,j:j+1,k:k+1))
            end do
         end do
      end do
      
      ! Apply a Neumann condition in non-periodic directions
	   if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            gradU(:,:,this%cfg%imin-1,:,:)=gradU(:,:,this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) gradU(:,:,this%cfg%imax+1,:,:)=gradU(:,:,this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            gradU(:,:,:,this%cfg%jmin-1,:)=gradU(:,:,:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) gradU(:,:,:,this%cfg%jmax+1,:)=gradU(:,:,:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            gradU(:,:,:,:,this%cfg%kmin-1)=gradU(:,:,:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) gradU(:,:,:,:,this%cfg%kmax+1)=gradU(:,:,:,:,this%cfg%kmax)
      end if
      
      ! Ensure zero in walls
	   do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) gradU(:,:,i,j,k)=0.0_WP
            end do
         end do
      end do
      
      ! Sync it
	   call this%cfg%sync(gradU)
      
      ! Deallocate velocity gradient storage
	   deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
      
   end subroutine get_gradU
   
   
   !> Calculate vorticity vector
   subroutine get_vorticity(this,vort)
      use messager, only: die
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: vort  !< Needs to be (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      
      ! Check vort's first two dimensions
      if (size(vort,dim=1).ne.3) call die('[tpcons get_vorticity] vort should be of size (1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)')
      
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
      
      ! Interpolate components to the cell center and assemble vorticity
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
      class(tpcons), intent(inout) :: this
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
      if (present(CFL)) cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z,this%CFLst)
      
   end subroutine get_cfl
   
   
   !> Calculate the max of our fields
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpcons), intent(inout) :: this
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
      class(tpcons), intent(inout) :: this
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
   subroutine correct_mfr(this,src)
      use mpi_f08, only: MPI_SUM
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), optional :: src !< Mass source term
      real(WP) :: mfr_error,vel_correction,int
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Evaluate MFR mismatch and velocity correction
      call this%get_mfr()
      mfr_error=sum(this%mfr)
      if (present(src)) then
         ! Also account for provided source term
         call this%cfg%integrate_without_VF(src,int)
         mfr_error=mfr_error-int
      end if
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
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%U(i,j,k)=this%U(i,j,k)+my_bc%rdir*vel_correction
               end do
            case ('y')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%V(i,j,k)=this%V(i,j,k)+my_bc%rdir*vel_correction
               end do
            case ('z')
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%W(i,j,k)=this%W(i,j,k)+my_bc%rdir*vel_correction
               end do
            end select
            
         end if
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
      ! Sync full fields
      call this%cfg%sync(this%U)
      call this%cfg%sync(this%V)
      call this%cfg%sync(this%W)
      
   end subroutine correct_mfr
   
   
   !> Shift pressure to ensure zero average
   subroutine shift_p(this,pressure)
      implicit none
      class(tpcons), intent(in) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: pressure !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: pressure_tot
      integer :: i,j,k
      ! Compute volume-averaged pressure
      call this%cfg%integrate(A=pressure,integral=pressure_tot); pressure_tot=pressure_tot/this%cfg%fluid_vol
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
      class(tpcons), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP) :: rhoUp,rhoUm,rhoVp,rhoVm,rhoWp,rhoWm
      
      ! If no implicit solver available, just divide by density and return
      if (.not.associated(this%implicit)) then
         resU=resU/this%RHOX
         resV=resV/this%RHOY
         resW=resW/this%RHOZ
         call this%cfg%sync(resU)
         call this%cfg%sync(resV)
         call this%cfg%sync(resW)
         return
      end if
      
      ! Solve implicit U problem
      this%implicit%opr(1,:,:,:)=this%RHOX; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=sum(this%itpu_x(:,i  ,j,k)*this%rhoU(i  :i+1,j,k))
               rhoUm=sum(this%itpu_x(:,i-1,j,k)*this%rhoU(i-1:i  ,j,k))
               rhoVp=sum(this%itpv_x(:,i,j+1,k)*this%rhoV(i-1:i,j+1,k))
               rhoVm=sum(this%itpv_x(:,i,j  ,k)*this%rhoV(i-1:i,j  ,k))
               rhoWp=sum(this%itpw_x(:,i,j,k+1)*this%rhoW(i-1:i,j,k+1))
               rhoWm=sum(this%itpw_x(:,i,j,k  )*this%rhoW(i-1:i,j,k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*(this%divu_x( 0,i,j,k)*this%itpu_x( 0,i  ,j,k)*rhoUp+&
               &                                                         this%divu_x(-1,i,j,k)*this%itpu_x(+1,i-1,j,k)*rhoUm+&
               &                                                         this%divu_y(+1,i,j,k)*this%itpu_y(-1,i,j+1,k)*rhoVp+&
               &                                                         this%divu_y( 0,i,j,k)*this%itpu_y( 0,i,j  ,k)*rhoVm+&
               &                                                         this%divu_z(+1,i,j,k)*this%itpu_z(-1,i,j,k+1)*rhoWp+&
               &                                                         this%divu_z( 0,i,j,k)*this%itpu_z( 0,i,j,k  )*rhoWm)*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+dt*(this%divu_x( 0,i,j,k)*this%itpu_x(+1,i  ,j,k)*rhoUp)*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+dt*(this%divu_x(-1,i,j,k)*this%itpu_x( 0,i-1,j,k)*rhoUm)*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+dt*(this%divu_y(+1,i,j,k)*this%itpu_y( 0,i,j+1,k)*rhoVp)*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+dt*(this%divu_y( 0,i,j,k)*this%itpu_y(-1,i,j  ,k)*rhoVm)*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+dt*(this%divu_z(+1,i,j,k)*this%itpu_z( 0,i,j,k+1)*rhoWp)*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+dt*(this%divu_z( 0,i,j,k)*this%itpu_z(-1,i,j,k  )*rhoWm)*0.5_WP
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-dt*(this%divu_x( 0,i,j,k)*2.0_WP*this%visc   (i  ,j,k)*this%grdu_x( 0,i  ,j,k)+&
               &                                                         this%divu_x(-1,i,j,k)*2.0_WP*this%visc   (i-1,j,k)*this%grdu_x(+1,i-1,j,k)+&
               &                                                         this%divu_y(+1,i,j,k)*       this%visc_xy(i,j+1,k)*this%grdu_y(-1,i,j+1,k)+&
               &                                                         this%divu_y( 0,i,j,k)*       this%visc_xy(i,j  ,k)*this%grdu_y( 0,i,j  ,k)+&
               &                                                         this%divu_z(+1,i,j,k)*       this%visc_zx(i,j,k+1)*this%grdu_z(-1,i,j,k+1)+&
               &                                                         this%divu_z( 0,i,j,k)*       this%visc_zx(i,j,k  )*this%grdu_z( 0,i,j,k  ))*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-dt*(this%divu_x( 0,i,j,k)*2.0_WP*this%visc   (i  ,j,k)*this%grdu_x(+1,i  ,j,k))*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-dt*(this%divu_x(-1,i,j,k)*2.0_WP*this%visc   (i-1,j,k)*this%grdu_x( 0,i-1,j,k))*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-dt*(this%divu_y(+1,i,j,k)*       this%visc_xy(i,j+1,k)*this%grdu_y( 0,i,j+1,k))*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-dt*(this%divu_y( 0,i,j,k)*       this%visc_xy(i,j  ,k)*this%grdu_y(-1,i,j  ,k))*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-dt*(this%divu_z(+1,i,j,k)*       this%visc_zx(i,j,k+1)*this%grdu_z( 0,i,j,k+1))*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-dt*(this%divu_z( 0,i,j,k)*       this%visc_zx(i,j,k  )*this%grdu_z(-1,i,j,k  ))*0.5_WP
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resU
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resU=this%implicit%sol
      
      ! Solve implicit V problem
      this%implicit%opr(1,:,:,:)=this%RHOY; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=sum(this%itpu_y(:,i+1,j,k)*this%rhoU(i+1,j-1:j,k))
               rhoUm=sum(this%itpu_y(:,i  ,j,k)*this%rhoU(i  ,j-1:j,k))
               rhoVp=sum(this%itpv_y(:,i,j  ,k)*this%rhoV(i,j  :j+1,k))
               rhoVm=sum(this%itpv_y(:,i,j-1,k)*this%rhoV(i,j-1:j  ,k))
               rhoWp=sum(this%itpw_y(:,i,j,k+1)*this%rhoW(i,j-1:j,k+1))
               rhoWm=sum(this%itpw_y(:,i,j,k  )*this%rhoW(i,j-1:j,k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*(this%divv_x(+1,i,j,k)*this%itpv_x(-1,i+1,j,k)*rhoUp+&
               &                                                         this%divv_x( 0,i,j,k)*this%itpv_x( 0,i  ,j,k)*rhoUm+&
               &                                                         this%divv_y( 0,i,j,k)*this%itpv_y( 0,i,j  ,k)*rhoVp+&
               &                                                         this%divv_y(-1,i,j,k)*this%itpv_y(+1,i,j-1,k)*rhoVm+&
               &                                                         this%divv_z(+1,i,j,k)*this%itpv_z(-1,i,j,k+1)*rhoWp+&
               &                                                         this%divv_z( 0,i,j,k)*this%itpv_z( 0,i,j,k  )*rhoWm)*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+dt*(this%divv_x(+1,i,j,k)*this%itpv_x( 0,i+1,j,k)*rhoUp)*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+dt*(this%divv_x( 0,i,j,k)*this%itpv_x(-1,i  ,j,k)*rhoUm)*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+dt*(this%divv_y( 0,i,j,k)*this%itpv_y(+1,i,j  ,k)*rhoVp)*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+dt*(this%divv_y(-1,i,j,k)*this%itpv_y( 0,i,j-1,k)*rhoVm)*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+dt*(this%divv_z(+1,i,j,k)*this%itpv_z( 0,i,j,k+1)*rhoWp)*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+dt*(this%divv_z( 0,i,j,k)*this%itpv_z(-1,i,j,k  )*rhoWm)*0.5_WP
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-dt*(this%divv_x(+1,i,j,k)*       this%visc_xy(i+1,j,k)*this%grdv_x(-1,i+1,j,k)+&
               &                                                         this%divv_x( 0,i,j,k)*       this%visc_xy(i  ,j,k)*this%grdv_x( 0,i  ,j,k)+&
               &                                                         this%divv_y( 0,i,j,k)*2.0_WP*this%visc   (i,j  ,k)*this%grdv_y( 0,i,j  ,k)+&
               &                                                         this%divv_y(-1,i,j,k)*2.0_WP*this%visc   (i,j-1,k)*this%grdv_y(+1,i,j-1,k)+&
               &                                                         this%divv_z(+1,i,j,k)*       this%visc_yz(i,j,k+1)*this%grdv_z(-1,i,j,k+1)+&
               &                                                         this%divv_z( 0,i,j,k)*       this%visc_yz(i,j,k  )*this%grdv_z( 0,i,j,k  ))*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-dt*(this%divv_x(+1,i,j,k)*       this%visc_xy(i+1,j,k)*this%grdv_x( 0,i+1,j,k))*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-dt*(this%divv_x( 0,i,j,k)*       this%visc_xy(i  ,j,k)*this%grdv_x(-1,i  ,j,k))*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-dt*(this%divv_y( 0,i,j,k)*2.0_WP*this%visc   (i,j  ,k)*this%grdv_y(+1,i,j  ,k))*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-dt*(this%divv_y(-1,i,j,k)*2.0_WP*this%visc   (i,j-1,k)*this%grdv_y( 0,i,j-1,k))*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-dt*(this%divv_z(+1,i,j,k)*       this%visc_yz(i,j,k+1)*this%grdv_z( 0,i,j,k+1))*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-dt*(this%divv_z( 0,i,j,k)*       this%visc_yz(i,j,k  )*this%grdv_z(-1,i,j,k  ))*0.5_WP
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resV
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resV=this%implicit%sol
      
      ! Solve implicit W problem
      this%implicit%opr(1,:,:,:)=this%RHOZ; this%implicit%opr(2:,:,:,:)=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               rhoUp=sum(this%itpu_z(:,i+1,j,k)*this%rhoU(i+1,j,k-1:k))
               rhoUm=sum(this%itpu_z(:,i  ,j,k)*this%rhoU(i  ,j,k-1:k))
               rhoVp=sum(this%itpv_z(:,i,j+1,k)*this%rhoV(i,j+1,k-1:k))
               rhoVm=sum(this%itpv_z(:,i,j  ,k)*this%rhoV(i,j  ,k-1:k))
               rhoWp=sum(this%itpw_z(:,i,j,k  )*this%rhoW(i,j,k  :k+1))
               rhoWm=sum(this%itpw_z(:,i,j,k-1)*this%rhoW(i,j,k-1:k  ))
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)+dt*(this%divw_x(+1,i,j,k)*this%itpw_x(-1,i+1,j,k)*rhoUp+&
               &                                                         this%divw_x( 0,i,j,k)*this%itpw_x( 0,i  ,j,k)*rhoUm+&
               &                                                         this%divw_y(+1,i,j,k)*this%itpw_y(-1,i,j+1,k)*rhoVp+&
               &                                                         this%divw_y( 0,i,j,k)*this%itpw_y( 0,i,j  ,k)*rhoVm+&
               &                                                         this%divw_z( 0,i,j,k)*this%itpw_z( 0,i,j,k  )*rhoWp+&
               &                                                         this%divw_z(-1,i,j,k)*this%itpw_z(+1,i,j,k-1)*rhoWm)*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)+dt*(this%divw_x(+1,i,j,k)*this%itpw_x( 0,i+1,j,k)*rhoUp)*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)+dt*(this%divw_x( 0,i,j,k)*this%itpw_x(-1,i  ,j,k)*rhoUm)*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)+dt*(this%divw_y(+1,i,j,k)*this%itpw_y( 0,i,j+1,k)*rhoVp)*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)+dt*(this%divw_y( 0,i,j,k)*this%itpw_y(-1,i,j  ,k)*rhoVm)*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)+dt*(this%divw_z( 0,i,j,k)*this%itpw_z(+1,i,j,k  )*rhoWp)*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)+dt*(this%divw_z(-1,i,j,k)*this%itpw_z( 0,i,j,k-1)*rhoWm)*0.5_WP
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-dt*(this%divw_x(+1,i,j,k)*       this%visc_zx(i+1,j,k)*this%grdw_x(-1,i+1,j,k)+&
               &                                                         this%divw_x( 0,i,j,k)*       this%visc_zx(i  ,j,k)*this%grdw_x( 0,i  ,j,k)+&
               &                                                         this%divw_y(+1,i,j,k)*       this%visc_yz(i,j+1,k)*this%grdw_y(-1,i,j+1,k)+&
               &                                                         this%divw_y( 0,i,j,k)*       this%visc_yz(i,j  ,k)*this%grdw_y( 0,i,j  ,k)+&
               &                                                         this%divw_z( 0,i,j,k)*2.0_WP*this%visc   (i,j,k  )*this%grdw_z( 0,i,j,k  )+&
               &                                                         this%divw_z(-1,i,j,k)*2.0_WP*this%visc   (i,j,k-1)*this%grdw_z(+1,i,j,k-1))*0.5_WP
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-dt*(this%divw_x(+1,i,j,k)*       this%visc_zx(i+1,j,k)*this%grdw_x( 0,i+1,j,k))*0.5_WP
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-dt*(this%divw_x( 0,i,j,k)*       this%visc_zx(i  ,j,k)*this%grdw_x(-1,i  ,j,k))*0.5_WP
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-dt*(this%divw_y(+1,i,j,k)*       this%visc_yz(i,j+1,k)*this%grdw_y( 0,i,j+1,k))*0.5_WP
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-dt*(this%divw_y( 0,i,j,k)*       this%visc_yz(i,j  ,k)*this%grdw_y(-1,i,j  ,k))*0.5_WP
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-dt*(this%divw_z( 0,i,j,k)*2.0_WP*this%visc   (i,j,k  )*this%grdw_z(+1,i,j,k  ))*0.5_WP
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-dt*(this%divw_z(-1,i,j,k)*2.0_WP*this%visc   (i,j,k-1)*this%grdw_z( 0,i,j,k-1))*0.5_WP
            end do
         end do
      end do
      call this%implicit%setup()
      this%implicit%rhs=resW
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resW=this%implicit%sol
      
   end subroutine solve_implicit
   
   
   !> Compute RHOX/Y/Z from passed rho field
   subroutine update_faceRHO(this,rho)
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: rho !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate square root of face densities
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_; do j=this%cfg%jmino_  ,this%cfg%jmaxo_; do i=this%cfg%imino_+1,this%cfg%imaxo_
         this%RHOX(i,j,k)=sum(this%itpr_x(:,i,j,k)*rho(i-1:i,j,k))
      end do; end do; end do
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_; do j=this%cfg%jmino_+1,this%cfg%jmaxo_; do i=this%cfg%imino_  ,this%cfg%imaxo_
         this%RHOY(i,j,k)=sum(this%itpr_y(:,i,j,k)*rho(i,j-1:j,k))
      end do; end do; end do
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_; do j=this%cfg%jmino_  ,this%cfg%jmaxo_; do i=this%cfg%imino_  ,this%cfg%imaxo_
         this%RHOZ(i,j,k)=sum(this%itpr_z(:,i,j,k)*rho(i,j,k-1:k))
      end do; end do; end do
      ! Handle non-periodic borders
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%RHOX(this%cfg%imino,:,:)=rho(this%cfg%imino,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%RHOY(:,this%cfg%jmino,:)=rho(:,this%cfg%jmino,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%RHOZ(:,:,this%cfg%kmino)=rho(:,:,this%cfg%kmino)
      ! Synchronize boundaries
      call this%cfg%sync(this%RHOX)
      call this%cfg%sync(this%RHOY)
      call this%cfg%sync(this%RHOZ)
   end subroutine update_faceRHO
   
   
   !> Prepare viscosity arrays from vfs object
   subroutine get_viscosity(this,vf,strat)
      use vfs_class, only: vfs
      use messager,  only: die
      implicit none
      class(tpcons), intent(inout) :: this
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
         call die('[tpcons get_viscosity] Unknown viscosity averaging strategy')
      end select
      ! Synchronize boundaries - not really needed...
      call this%cfg%sync(this%visc)
      call this%cfg%sync(this%visc_xy)
      call this%cfg%sync(this%visc_yz)
      call this%cfg%sync(this%visc_zx)
   end subroutine get_viscosity
   
   
   !> Add gravity source term
   subroutine addsrc_gravity(this,resU,resV,resW)
      implicit none
      class(tpcons), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+this%gravity(1)*this%RHOX(i,j,k)
               if (this%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+this%gravity(2)*this%RHOY(i,j,k)
               if (this%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+this%gravity(3)*this%RHOZ(i,j,k)
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
      class(tpcons), intent(inout) :: this
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
   subroutine tpcons_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(tpcons), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Two-phase incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   liquid density = ",es12.5)') this%rho_l
         write(output_unit,'(" > liquid viscosity = ",es12.5)') this%visc_l
         write(output_unit,'(" >      gas density = ",es12.5)') this%rho_g
         write(output_unit,'(" >    gas viscosity = ",es12.5)') this%visc_g
      end if
      
   end subroutine tpcons_print
   
   
end module tpcons_class
