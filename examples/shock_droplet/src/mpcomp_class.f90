!> Multiphase compressible flow solver class:
!> Provides support for RHS calculation only
module mpcomp_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   use timer_class,  only: timer
   use irl_fortran_interface
   implicit none
   private
   
   ! Expose type
   public :: mpcomp,VFlo,VFhi
   
   ! Default parameters for volume fraction solver
   integer,  parameter :: half_edge=1                                !< Half-edge cutting (default)
   real(WP), parameter :: VFlo=1.0e-12_WP                            !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo                           !< Maximum VF value considered
   real(WP), parameter :: volume_epsilon_factor =1.0e-15_WP          !< Minimum volume  to consider for computational geometry (normalized by min_meshsize**3)
   real(WP), parameter :: surface_epsilon_factor=1.0e-15_WP          !< Minimum surface to consider for computational geometry (normalized by min_meshsize**2)
   real(WP), parameter :: iterative_distfind_tol=1.0e-12_WP          !< Tolerance for iterative plane distance finding
   
   !> Multiphase compressible solver object definition
   type :: mpcomp
      
      ! This is the config around which solver is built
      class(config), pointer :: cfg
      
      ! Solver name
      character(len=str_medium) :: name='UNNAMED_MPCOMP'
      
      ! Pointers to functions to evaluate P(RHO,E), T(RHO,P), and C(RHO,P)
      procedure(Pfunc_type), pointer, nopass :: getPL=>NULL()
      procedure(Tfunc_type), pointer, nopass :: getTL=>NULL()
      procedure(Cfunc_type), pointer, nopass :: getCL=>NULL()
      procedure(Sfunc_type), pointer, nopass :: getSL=>NULL()
      procedure(Pfunc_type), pointer, nopass :: getPG=>NULL()
      procedure(Tfunc_type), pointer, nopass :: getTG=>NULL()
      procedure(Cfunc_type), pointer, nopass :: getCG=>NULL()
      procedure(Sfunc_type), pointer, nopass :: getSG=>NULL()
      
      ! Pointer to subroutine for mixture cell relaxation
      procedure(relax_type), pointer, nopass :: relax=>NULL()
      
      ! Volume moments, interface, and semi-Lagrangian fluxes
      real(WP), dimension(:,:,:)  , allocatable :: VF,VFold
      real(WP), dimension(:,:,:,:), allocatable :: BL,BLold
      real(WP), dimension(:,:,:,:), allocatable :: BG,BGold
      type(PlanarSep_type), dimension(:,:,:), allocatable :: PLIC,PLICold
      
      ! Tag for semi-Lagrangian fluxing
      integer, dimension(:,:,:), allocatable :: iSL
      
      ! Conserved variables: 1=VF*RHOL, 2=(1-VF))*RHOG, 3=VF*RHOL*IL, 4=(1-VF)*RHOG*IL, 5=RHO*U, 6=RHO*V, 7=RHO*W
      integer :: nQ
      real(WP), dimension(:,:,:,:), allocatable :: Q,Qold
      
      ! Flow velocity
      real(WP), dimension(:,:,:), allocatable :: U,V,W
      
      ! Phasic densities
      real(WP), dimension(:,:,:), allocatable :: RHOL,RHOG,RHOLold,RHOGold
      
      ! Phasic internal energies
      real(WP), dimension(:,:,:), allocatable :: IL,IG,ILold,IGold
      
      ! Phasic pressures
      real(WP), dimension(:,:,:), allocatable :: PL,PG
      
      ! Phasic temperatures
      real(WP), dimension(:,:,:), allocatable :: TL,TG
      
      ! Mixture speed of sound
      real(WP), dimension(:,:,:), allocatable :: C
      
      ! Viscosities and heat diffusivity
      real(WP), dimension(:,:,:), allocatable :: VISCL,BETAL,DIFFL
      real(WP), dimension(:,:,:), allocatable :: VISCG,BETAG,DIFFG
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP
      
      ! Store mesh info
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      
      ! CFL numbers
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                    !< Convective CFL numbers
      real(WP) :: CFLa_x,CFLa_y,CFLa_z                    !< Acoustic CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                    !< Viscous CFL numbers
      
      ! Monitoring quantities for conserved variables
      real(WP) :: VFmin,VFmax,VFint
      real(WP), dimension(:), allocatable :: Qmin,Qmax,Qint
      real(WP) :: RHOKLint,RHOKGint
      real(WP) :: RHOSLint,RHOSGint
      
      ! Monitoring quantities for primitive variables
      real(WP) :: Umax,Vmax,Wmax                          !< Velocity stats
      real(WP) :: RHOLmin,RHOLmax                         !< Liquid density stats
      real(WP) :: RHOGmin,RHOGmax                         !< Gas    density stats
      real(WP) :: ILmin,ILmax                             !< Liquid internal energy stats
      real(WP) :: IGmin,IGmax                             !< Gas    internal energy stats
      real(WP) :: PLmin,PLmax                             !< Liquid pressure stats
      real(WP) :: PGmin,PGmax                             !< Gas    pressure stats
      real(WP) :: TLmin,TLmax                             !< Liquid temperature stats
      real(WP) :: TGmin,TGmax                             !< Gas    temperature stats
      
      ! Timers
      type(timer) :: trhs                                 !< Timer for RHS calculation (excluding tagged cells)
      type(timer) :: tsl                                  !< Timer for semi-Lagrangian transport
      type(timer) :: tplic                                !< Timer for PLIC reconstruction
      
      ! Volume fluxes, SL phasic fluxes, SL time step size
      type(SepVM_type),     dimension(:,:,:), allocatable :: FVx,FVy,FVz
      real(WP), dimension(:,:,:,:), allocatable :: SLFx,SLFy,SLFz
      real(WP) :: SLdt
      
      ! IRL-native data
      type(ByteBuffer_type) :: send_byte_buffer,recv_byte_buffer
      type(ObjServer_PlanarSep_type)  :: planar_separator_allocation,old_planar_separator_allocation
      type(ObjServer_PlanarLoc_type)  :: planar_localizer_allocation
      type(ObjServer_LocSepLink_type) :: localized_separator_link_allocation
      type(ObjServer_LocLink_type)    :: localizer_link_allocation
      type(PlanarLoc_type),  dimension(:,:,:), allocatable :: localizer
      type(LocSepLink_type), dimension(:,:,:), allocatable :: localized_separator_link
      type(LocLink_type),    dimension(:,:,:), allocatable :: localizer_link
      type(Poly_type),       dimension(:,:,:), allocatable :: interface_polygon
      
   contains
      procedure :: print=>mpcomp_print                    !< Output solver to the screen
      procedure :: initialize                             !< Initialize the flow solver
      procedure :: initialize_irl                         !< Initialize interface with interface reconstruction library
      procedure :: SLtag                                  !< Tag cell for semi-Lagrangian fluxing
      procedure :: SLstep                                 !< Perform an unsplit semi-Lagrangian transport step in tagged cells, store SL advection fluxes
      procedure :: SLincrement                            !< Increment Q by SLdt*rhs of all equations with rhs calculated using stored SL advection fluxes
      procedure :: rhs                                    !< Compute rhs of our equations using standard fluxes
      procedure :: build_interface                        !< Build interface from volume moments
      procedure :: get_primitive                          !< Calculate phasic and mixture primitive variables from conserved variables
      procedure :: apply_relax                            !< Apply user-provided relaxation model in interfacial cells
      procedure :: get_velocity                           !< Calculate velocity from momentum
      procedure :: get_ke                                 !< Calculate kinetic energy per unit mass from velocity
      procedure :: get_momentum                           !< Calculate momentum from velocity
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_info                               !< Calculate maximum field values
      procedure :: update_surfmesh                        !< Update a surfmesh object using current polygons
   end type mpcomp
   
   !> Interfaces for user-defined function
   abstract interface
      !> P=P(RHO,I)
      pure real(WP) function Pfunc_type(RHO,I)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: I
      end function Pfunc_type
      !> T=T(RHO,P)
      pure real(WP) function Tfunc_type(RHO,P)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: P
      end function Tfunc_type
      !> C=C(RHO,P)
      pure real(WP) function Cfunc_type(RHO,P)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: P
      end function Cfunc_type
      !> S=S(RHO,P)
      pure real(WP) function Sfunc_type(RHO,P)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: P
      end function Sfunc_type
      !> Pressure relaxation (acts only on the conserved quantities)
      subroutine relax_type(VF,Q)
         import :: WP
         implicit none
         real(WP),                intent(inout) :: VF
         real(WP), dimension(1:), intent(inout) :: Q
      end subroutine relax_type
   end interface
   
contains
   
   
   !> Initialization for compressible flow solver
   subroutine initialize(this,cfg,getPL,getCL,getPG,getCG,name)
      use messager, only: die
      implicit none
      class(mpcomp) :: this
      class(config), target, intent(in) :: cfg
      procedure(Pfunc_type)  :: getPL,getPG
      procedure(Cfunc_type)  :: getCL,getCG
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to config object
      this%cfg=>cfg
      
      ! Point to thermodynamic functions
      this%getPL=>getPL
      this%getCL=>getCL
      this%getPG=>getPG
      this%getCG=>getCG
      
      ! Check that config is uniform with at least 3 cells of overlap
      if (this%cfg%no.lt.3) call die('[mpcomp initialize] mpcomp solver requires at least 3 cells of overlap')
      if (.not.all([this%cfg%uniform_x,this%cfg%uniform_y,this%cfg%uniform_z])) call die('[mpcomp initialize] mpcomp solver requires a uniform mesh')
      
      ! Store constant cell size and its inverse, handle 2D conditions
      this%dx=this%cfg%dx(this%cfg%imin_); this%dxi=1.0_WP/this%dx; if (this%cfg%nx.eq.1) this%dxi=0.0_WP
      this%dy=this%cfg%dy(this%cfg%jmin_); this%dyi=1.0_WP/this%dy; if (this%cfg%ny.eq.1) this%dyi=0.0_WP
      this%dz=this%cfg%dz(this%cfg%kmin_); this%dzi=1.0_WP/this%dz; if (this%cfg%nz.eq.1) this%dzi=0.0_WP
      
      ! Allocate and zero out conserved variables
      this%nQ=7
      allocate(this%Q   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); this%Q   =0.0_WP
      allocate(this%Qold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); this%Qold=0.0_WP
      
      ! Allocate volume moments
      allocate(this%VF    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VF=0.0_WP
      allocate(this%BL(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BL=0.0_WP
      allocate(this%BG(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BG=0.0_WP
      allocate(this%VFold    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VFold=0.0_WP
      allocate(this%BLold(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BLold=0.0_WP
      allocate(this%BGold(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BGold=0.0_WP
      
      ! Allocate semi-Lagrangian fluxes
      allocate(this%SLFx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:4)); this%SLFx=0.0_WP
      allocate(this%SLFy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:4)); this%SLFy=0.0_WP
      allocate(this%SLFz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:4)); this%SLFz=0.0_WP
      
      ! Initialize Interface Reconstruction Library and its data
      call this%initialize_irl()
      
      ! Allocate semi-Lagrangian tag
      allocate(this%iSL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%iSL=0
      
      ! Conserved variables monitoring
      allocate(this%Qmin(1:this%nQ),this%Qmax(1:this%nQ),this%Qint(1:this%nQ))
      
      ! Flow variables
      allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      
      ! Phasic densities
      allocate(this%RHOL   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOL   =0.0_WP
      allocate(this%RHOG   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOG   =0.0_WP
      allocate(this%RHOLold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOLold=0.0_WP
      allocate(this%RHOGold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOGold=0.0_WP
      
      ! Phasic internal energies
      allocate(this%IL   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%IL   =0.0_WP
      allocate(this%IG   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%IG   =0.0_WP
      allocate(this%ILold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%ILold=0.0_WP
      allocate(this%IGold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%IGold=0.0_WP
      
      ! Phasic temperatures
      allocate(this%TL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%TL=0.0_WP
      allocate(this%TG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%TG=0.0_WP
      
      ! Phasic pressures
      allocate(this%PL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%PL=0.0_WP
      allocate(this%PG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%PG=0.0_WP
      
      ! Mixture speed of sound
      allocate(this%C(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%C=0.0_WP
      
      ! Fluid viscosities and heat diffusivity
      allocate(this%VISCL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VISCL=0.0_WP
      allocate(this%BETAL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BETAL=0.0_WP
      allocate(this%DIFFL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%DIFFL=0.0_WP
      allocate(this%VISCG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VISCG=0.0_WP
      allocate(this%BETAG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BETAG=0.0_WP
      allocate(this%DIFFG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%DIFFG=0.0_WP
      
   end subroutine initialize
   
   
   !> Initialize IRL interface
   subroutine initialize_irl(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k,n,tag
      real(WP), dimension(3,4) :: vert
      integer(IRL_LargeOffsetIndex_t) :: total_cells
      
      ! Transfer small constants to IRL
      call setVFBounds(VFlo)
      call setVFTolerance_IterativeDistanceFinding(iterative_distfind_tol)
      call setMinimumVolToTrack(volume_epsilon_factor*this%cfg%min_meshsize**3)
      call setMinimumSAToTrack(surface_epsilon_factor*this%cfg%min_meshsize**2)
      
      ! Set IRL's moment calculation method
      call getMoments_setMethod(half_edge)
      
      ! Allocate IRL storage for PLIC interface
      allocate(this%PLIC   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%PLICold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Allocate IRL storage for volume fluxes
      allocate(this%FVx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%FVy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%FVz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         call new(this%FVx(i,j,k))
         call new(this%FVy(i,j,k))
         call new(this%FVz(i,j,k))
      end do; end do; end do
      
      ! Allocate IRL arrays
      allocate(this%localizer               (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%localized_separator_link(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%localizer_link          (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%interface_polygon       (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Initialize size for IRL
      total_cells=int(this%cfg%nxo_,8)*int(this%cfg%nyo_,8)*int(this%cfg%nzo_,8)
      call new(this%planar_localizer_allocation,total_cells)
      call new(this%planar_separator_allocation,total_cells)
      call new(this%old_planar_separator_allocation,total_cells)
      call new(this%localized_separator_link_allocation,total_cells)
      call new(this%localizer_link_allocation,total_cells)
      
      ! Initialize arrays and setup linking
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         ! Transfer cell to IRL
         call new(this%localizer(i,j,k),this%planar_localizer_allocation)
         call setFromRectangularCuboid(this%localizer(i,j,k),[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
         ! PLIC interface
         call new(this%PLIC   (i,j,k),this%planar_separator_allocation   )
         call new(this%PLICold(i,j,k),this%old_planar_separator_allocation)
         ! PLIC+mesh with connectivity (i.e., link) - this uses PLICold
         call new(this%localized_separator_link(i,j,k),this%localized_separator_link_allocation,this%localizer(i,j,k),this%PLICold(i,j,k))
         ! Mesh with connectivity
         call new(this%localizer_link(i,j,k),this%localizer_link_allocation,this%localizer(i,j,k))
      end do; end do; end do
      
      ! Polygonal representation of the surface
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         call new(this%interface_polygon(i,j,k))
      end do; end do; end do
      
      ! Give each link a unique lexicographic tag (per processor)
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         tag=this%cfg%get_lexico_from_ijk([i,j,k])
         call setId(this%localized_separator_link(i,j,k),tag)
         call setId(this%localizer_link(i,j,k),tag)
      end do; end do; end do
      
      ! Set the connectivity
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         ! In the x- direction
         if (i.gt.this%cfg%imino_) then
            call setEdgeConnectivity(this%localized_separator_link(i,j,k),0,this%localized_separator_link(i-1,j,k))
            call setEdgeConnectivity(this%localizer_link(i,j,k),0,this%localizer_link(i-1,j,k))
         end if
         ! In the x+ direction
         if (i.lt.this%cfg%imaxo_) then
            call setEdgeConnectivity(this%localized_separator_link(i,j,k),1,this%localized_separator_link(i+1,j,k))
            call setEdgeConnectivity(this%localizer_link(i,j,k),1,this%localizer_link(i+1,j,k))
         end if
         ! In the y- direction
         if (j.gt.this%cfg%jmino_) then
            call setEdgeConnectivity(this%localized_separator_link(i,j,k),2,this%localized_separator_link(i,j-1,k))
            call setEdgeConnectivity(this%localizer_link(i,j,k),2,this%localizer_link(i,j-1,k))
         end if
         ! In the y+ direction
         if (j.lt.this%cfg%jmaxo_) then
            call setEdgeConnectivity(this%localized_separator_link(i,j,k),3,this%localized_separator_link(i,j+1,k))
            call setEdgeConnectivity(this%localizer_link(i,j,k),3,this%localizer_link(i,j+1,k))
         end if
         ! In the z- direction
         if (k.gt.this%cfg%kmino_) then
            call setEdgeConnectivity(this%localized_separator_link(i,j,k),4,this%localized_separator_link(i,j,k-1))
            call setEdgeConnectivity(this%localizer_link(i,j,k),4,this%localizer_link(i,j,k-1))
         end if
         ! In the z+ direction
         if (k.lt.this%cfg%kmaxo_) then
            call setEdgeConnectivity(this%localized_separator_link(i,j,k),5,this%localized_separator_link(i,j,k+1))
            call setEdgeConnectivity(this%localizer_link(i,j,k),5,this%localizer_link(i,j,k+1))
         end if
      end do; end do; end do
      
      ! Prepare byte storage for synchronization
      call new(this%send_byte_buffer)
      call new(this%recv_byte_buffer)
      
   end subroutine initialize_irl
   
   
   !> Tag cell for semi-Lagrangian fluxing based on proximity to 0<VOF<1
   !> Could become user-specified, use VOFold, or depend on current CFL...
   subroutine SLtag(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k,dir,n
      integer, dimension(3) :: ind
      ! Reset tag
      this%iSL=0
      ! First sweep to identify cells with interface
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1
            do i=this%cfg%imino_+1,this%cfg%imaxo_-1
               ! Flag all mixture cells
               if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) this%iSL(i,j,k)=1
               ! Check if cell-face is an interface
               do dir=1,3
                  do n=-1,+1,2
                     ind=[i,j,k]; ind(dir)=ind(dir)+n
                     if (this%VF(i,j,k).lt.VFlo.and.this%VF(ind(1),ind(2),ind(3)).gt.VFhi.or.&
                     &   this%VF(i,j,k).gt.VFhi.and.this%VF(ind(1),ind(2),ind(3)).lt.VFlo) this%iSL(i,j,k)=1
                  end do
               end do
            end do
         end do
      end do
      call this%cfg%sync(this%iSL)
      ! Second sweep to extend by one cell
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1; do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1; do i=this%cfg%imino_+1,this%cfg%imaxo_-1
         if (this%iSL(i,j,k).eq.0.and.any(this%iSL(i-1:i+1,j-1:j+1,k-1:k+1).eq.1)) this%iSL(i,j,k)=2
      end do; end do; end do
      call this%cfg%sync(this%iSL)
   end subroutine SLtag
   
   
   !> Perform an unsplit semi-Lagrangian transport step by dt in all cells tagged by this%SLtag>0
   !> Volume moments are updated and advection fluxes for phasic equations are computed
   !> Uses VFold, BLold, BGold, PLICold, RHOLold, RHOGold, ILold, IGold
   !> Vertex transport is done with RK2 using passed (U,V,W)
   subroutine SLstep(this,dt,U,V,W)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W
      real(IRL_double), dimension(3,9) :: face
      type(SepVM_type) :: my_SepVM
      type(CapDod_type) :: flux_polyhedron
      type(TagAccVM_SepVM_type) :: detailed_face_flux
      real(WP), dimension(3) :: lbar,gbar
      integer , dimension(3) :: ind
      real(WP) :: lvol,gvol,lrho,grho,lrhoi,grhoi
      real(WP) :: Lvolold,Gvolold,Lvolinc,Gvolinc,Lvolnew,Gvolnew
      integer :: i,j,k,n
      
      ! Start semi-Lagrangian timer
      call this%tsl%start()
      
      ! Remember dt used
      this%SLdt=dt
      
      ! Allocate flux polyhedron and detailed face flux
      call new(flux_polyhedron)
      call new(detailed_face_flux)
      
      ! Zero out semi-Lagrangian fluxes
      this%SLFx=0.0_WP; this%SLFy=0.0_WP; this%SLFz=0.0_WP
      
      ! Loop through all cell faces and get volume moment, mass, and internal energy fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! X flux
               if (maxval(this%iSL(i-1:i,j,k)).gt.0) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=project(face(:,1),-dt)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=project(face(:,2),-dt)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=project(face(:,3),-dt)
                  face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=project(face(:,4),-dt)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]; face(:,9)=project(face(:,9),-dt)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%dy*this%dz)
                  ! Build detailed geometric flux
                  call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),detailed_face_flux)
                  ! Build finite volume fluxes from detailed flux
                  lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP; lrho=0.0_WP; grho=0.0_WP; lrhoi=0.0_WP; grhoi=0.0_WP
                  do n=0,getSize(detailed_face_flux)-1
                     ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux,n))
                     call getSepVMAtIndex(detailed_face_flux,n,my_SepVM)
                     lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0); lrho=lrho+getVolume(my_SepVM,0)*this%RHOLold(ind(1),ind(2),ind(3)); lrhoi=lrhoi+getVolume(my_SepVM,0)*this%RHOLold(ind(1),ind(2),ind(3))*this%ILold(ind(1),ind(2),ind(3))
                     gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1); grho=grho+getVolume(my_SepVM,1)*this%RHOGold(ind(1),ind(2),ind(3)); grhoi=grhoi+getVolume(my_SepVM,1)*this%RHOGold(ind(1),ind(2),ind(3))*this%IGold(ind(1),ind(2),ind(3))
                  end do
                  call construct(this%FVx(i,j,k),[lvol,lbar,gvol,gbar])
                  this%SLFx(i,j,k,1:4)=-[lrho,grho,lrhoi,grhoi]/(dt*this%dy*this%dz)
                  ! Clear detailed flux
                  call clear(detailed_face_flux)
               end if
               
               ! Y flux
               if (maxval(this%iSL(i,j-1:j,k)).gt.0) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=project(face(:,1),-dt)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=project(face(:,2),-dt)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=project(face(:,3),-dt)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=project(face(:,4),-dt)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]; face(:,9)=project(face(:,9),-dt)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%dx*this%dz)
                  ! Build detailed geometric flux
                  call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),detailed_face_flux)
                  ! Build finite volume fluxes from detailed flux
                  lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP; lrho=0.0_WP; grho=0.0_WP; lrhoi=0.0_WP; grhoi=0.0_WP
                  do n=0,getSize(detailed_face_flux)-1
                     ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux,n))
                     call getSepVMAtIndex(detailed_face_flux,n,my_SepVM)
                     lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0); lrho=lrho+getVolume(my_SepVM,0)*this%RHOLold(ind(1),ind(2),ind(3)); lrhoi=lrhoi+getVolume(my_SepVM,0)*this%RHOLold(ind(1),ind(2),ind(3))*this%ILold(ind(1),ind(2),ind(3))
                     gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1); grho=grho+getVolume(my_SepVM,1)*this%RHOGold(ind(1),ind(2),ind(3)); grhoi=grhoi+getVolume(my_SepVM,1)*this%RHOGold(ind(1),ind(2),ind(3))*this%IGold(ind(1),ind(2),ind(3))
                  end do
                  call construct(this%FVy(i,j,k),[lvol,lbar,gvol,gbar])
                  this%SLFy(i,j,k,1:4)=-[lrho,grho,lrhoi,grhoi]/(dt*this%dx*this%dz)
                  ! Clear detailed flux
                  call clear(detailed_face_flux)
               end if
               
               ! Z flux
               if (maxval(this%iSL(i,j,k-1:k)).gt.0) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=project(face(:,1),-dt)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=project(face(:,2),-dt)
                  face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=project(face(:,3),-dt)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=project(face(:,4),-dt)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]; face(:,9)=project(face(:,9),-dt)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%dx*this%dy)
                  ! Build detailed geometric flux
                  call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),detailed_face_flux)
                  ! Build finite volume fluxes from detailed flux
                  lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP; lrho=0.0_WP; grho=0.0_WP; lrhoi=0.0_WP; grhoi=0.0_WP
                  do n=0,getSize(detailed_face_flux)-1
                     ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux,n))
                     call getSepVMAtIndex(detailed_face_flux,n,my_SepVM)
                     lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0); lrho=lrho+getVolume(my_SepVM,0)*this%RHOLold(ind(1),ind(2),ind(3)); lrhoi=lrhoi+getVolume(my_SepVM,0)*this%RHOLold(ind(1),ind(2),ind(3))*this%ILold(ind(1),ind(2),ind(3))
                     gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1); grho=grho+getVolume(my_SepVM,1)*this%RHOGold(ind(1),ind(2),ind(3)); grhoi=grhoi+getVolume(my_SepVM,1)*this%RHOGold(ind(1),ind(2),ind(3))*this%IGold(ind(1),ind(2),ind(3))
                  end do
                  call construct(this%FVz(i,j,k),[lvol,lbar,gvol,gbar])
                  this%SLFz(i,j,k,1:4)=-[lrho,grho,lrhoi,grhoi]/(dt*this%dx*this%dy)
                  ! Clear detailed flux
                  call clear(detailed_face_flux)
               end if
               
            end do
         end do
      end do
      
      ! Mass fluxes will be used to build momentum fluxes, they need to be extended by one cell on the left because of staggering
      call this%cfg%sync(this%SLFx(:,:,:,1)); call this%cfg%sync(this%SLFx(:,:,:,2)); if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%SLFx(this%cfg%imin-1,:,:,1:2)=this%SLFx(this%cfg%imin,:,:,1:2)
      call this%cfg%sync(this%SLFy(:,:,:,1)); call this%cfg%sync(this%SLFy(:,:,:,2)); if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%SLFy(:,this%cfg%jmin-1,:,1:2)=this%SLFy(:,this%cfg%jmin,:,1:2)
      call this%cfg%sync(this%SLFz(:,:,:,1)); call this%cfg%sync(this%SLFz(:,:,:,2)); if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%SLFz(:,:,this%cfg%kmin-1,1:2)=this%SLFz(:,:,this%cfg%kmin,1:2)
      
      ! Update volume moments in tagged cells
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Only work at interface
               if (this%iSL(i,j,k).eq.0) cycle
               ! Old liquid and gas volumes
               Lvolold=        this%VFold(i,j,k) *this%cfg%vol(i,j,k)
               Gvolold=(1.0_WP-this%VFold(i,j,k))*this%cfg%vol(i,j,k)
               ! Compute incoming liquid and gas volumes
               Lvolinc=-getVolumePtr(this%FVx(i+1,j,k),0)+getVolumePtr(this%FVx(i,j,k),0)&
               &       -getVolumePtr(this%FVy(i,j+1,k),0)+getVolumePtr(this%FVy(i,j,k),0)&
               &       -getVolumePtr(this%FVz(i,j,k+1),0)+getVolumePtr(this%FVz(i,j,k),0)
               Gvolinc=-getVolumePtr(this%FVx(i+1,j,k),1)+getVolumePtr(this%FVx(i,j,k),1)&
               &       -getVolumePtr(this%FVy(i,j+1,k),1)+getVolumePtr(this%FVy(i,j,k),1)&
               &       -getVolumePtr(this%FVz(i,j,k+1),1)+getVolumePtr(this%FVz(i,j,k),1)
               ! Compute new liquid and gas volumes
               Lvolnew=Lvolold+Lvolinc
               Gvolnew=Gvolold+Gvolinc
               ! Compute new liquid volume fraction
               this%VF(i,j,k)=Lvolnew/(Lvolnew+Gvolnew)
               ! Clip it within [VFlo,VFhi] and compute first order moments
               if (this%VF(i,j,k).le.VFlo) then
                  this%VF(i,j,k)=0.0_WP
                  this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               else if (this%VF(i,j,k).ge.VFhi) then
                  this%VF(i,j,k)=1.0_WP
                  this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               else
                  ! Compute old liquid barycenter
                  this%BL(:,i,j,k)=(this%BLold(:,i,j,k)*Lvolold-getCentroidPtr(this%FVx(i+1,j,k),0)+getCentroidPtr(this%FVx(i,j,k),0)&
                  &                                            -getCentroidPtr(this%FVy(i,j+1,k),0)+getCentroidPtr(this%FVy(i,j,k),0)&
                  &                                            -getCentroidPtr(this%FVz(i,j,k+1),0)+getCentroidPtr(this%FVz(i,j,k),0))/Lvolnew
                  ! Project it forward in time
                  this%BL(:,i,j,k)=project(this%BL(:,i,j,k),dt)
                  ! Compute old gas phase barycenter
                  this%BG(:,i,j,k)=(this%BGold(:,i,j,k)*Gvolold-getCentroidPtr(this%FVx(i+1,j,k),1)+getCentroidPtr(this%FVx(i,j,k),1)&
                  &                                            -getCentroidPtr(this%FVy(i,j+1,k),1)+getCentroidPtr(this%FVy(i,j,k),1)&
                  &                                            -getCentroidPtr(this%FVz(i,j,k+1),1)+getCentroidPtr(this%FVz(i,j,k),1))/Gvolnew
                  ! Project it forward in time
                  this%BG(:,i,j,k)=project(this%BG(:,i,j,k),dt)
               end if
            end do
         end do
      end do
      
      ! Synchronize VF and barycenters
      call this%cfg%sync(this%VF); call this%cfg%sync(this%BL); call this%cfg%sync(this%BG)
      
      ! Fix barycenter synchronization across periodic boundaries
      if (this%cfg%xper.and.this%cfg%iproc.eq.1           ) then; do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino,this%cfg%imin-1
         this%BL(1,i,j,k)=this%BL(1,i,j,k)-this%cfg%xL; this%BG(1,i,j,k)=this%BG(1,i,j,k)-this%cfg%xL
      end do; end do; end do; end if
      if (this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then; do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imax+1,this%cfg%imaxo
         this%BL(1,i,j,k)=this%BL(1,i,j,k)+this%cfg%xL; this%BG(1,i,j,k)=this%BG(1,i,j,k)+this%cfg%xL
      end do; end do; end do; end if
      if (this%cfg%yper.and.this%cfg%jproc.eq.1           ) then; do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino,this%cfg%jmin-1; do i=this%cfg%imino_,this%cfg%imaxo_
         this%BL(2,i,j,k)=this%BL(2,i,j,k)-this%cfg%yL; this%BG(2,i,j,k)=this%BG(2,i,j,k)-this%cfg%yL
      end do; end do; end do; end if
      if (this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then; do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmax+1,this%cfg%jmaxo; do i=this%cfg%imino_,this%cfg%imaxo_
         this%BL(2,i,j,k)=this%BL(2,i,j,k)+this%cfg%yL; this%BG(2,i,j,k)=this%BG(2,i,j,k)+this%cfg%yL
      end do; end do; end do; end if
      if (this%cfg%zper.and.this%cfg%kproc.eq.           1) then; do k=this%cfg%kmino,this%cfg%kmin-1; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         this%BL(3,i,j,k)=this%BL(3,i,j,k)-this%cfg%zL; this%BG(3,i,j,k)=this%BG(3,i,j,k)-this%cfg%zL
      end do; end do; end do; end if
      if (this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then; do k=this%cfg%kmax+1,this%cfg%kmaxo; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         this%BL(3,i,j,k)=this%BL(3,i,j,k)+this%cfg%zL; this%BG(3,i,j,k)=this%BG(3,i,j,k)+this%cfg%zL
      end do; end do; end do; end if
      
      ! Handle 2D barycenters
      if (this%cfg%nx.eq.1) then; do i=this%cfg%imino_,this%cfg%imaxo_; this%BL(1,i,:,:)=this%cfg%xm(i); this%BG(1,i,:,:)=this%cfg%xm(i); end do; end if
      if (this%cfg%ny.eq.1) then; do j=this%cfg%jmino_,this%cfg%jmaxo_; this%BL(2,:,j,:)=this%cfg%ym(j); this%BG(2,:,j,:)=this%cfg%ym(j); end do; end if
      if (this%cfg%nz.eq.1) then; do k=this%cfg%kmino_,this%cfg%kmaxo_; this%BL(3,:,:,k)=this%cfg%zm(k); this%BG(3,:,:,k)=this%cfg%zm(k); end do; end if
      
      ! Stop semi-Lagrangian timer
      call this%tsl%stop()
      
   contains
      !> Project function that moves a point p1 at (U,V,W) for mydt
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project
      !> Function that performs trilinear interpolation of staggered velocity to pos
      !> This version assumes a uniform mesh for maximum speed
      function interp_velocity(pos) result(vel)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer  :: ip,jp,kp
         real(WP) :: x0,y0,z0,wx1,wy1,wz1,wx2,wy2,wz2
         ! Perform tri-linear interpolation of U component
         ip=min(max(floor((pos(1)-this%cfg%x (this%cfg%imino))/this%dx)+this%cfg%imino,this%cfg%imino_),this%cfg%imaxo_-1); wx1=min(max((pos(1)-this%cfg%x (ip))*this%dxi,0.0_WP),1.0_WP); wx2=1.0_WP-wx1
         jp=min(max(floor((pos(2)-this%cfg%ym(this%cfg%jmino))/this%dy)+this%cfg%jmino,this%cfg%jmino_),this%cfg%jmaxo_-1); wy1=min(max((pos(2)-this%cfg%ym(jp))*this%dyi,0.0_WP),1.0_WP); wy2=1.0_WP-wy1
         kp=min(max(floor((pos(3)-this%cfg%zm(this%cfg%kmino))/this%dz)+this%cfg%kmino,this%cfg%kmino_),this%cfg%kmaxo_-1); wz1=min(max((pos(3)-this%cfg%zm(kp))*this%dzi,0.0_WP),1.0_WP); wz2=1.0_WP-wz1
         vel(1)=wz1*(wy1*(wx1*U(ip+1,jp+1,kp+1)+wx2*U(ip,jp+1,kp+1))+wy2*(wx1*U(ip+1,jp,kp+1)+wx2*U(ip,jp,kp+1)))+wz2*(wy1*(wx1*U(ip+1,jp+1,kp)+wx2*U(ip,jp+1,kp))+wy2*(wx1*U(ip+1,jp,kp)+wx2*U(ip,jp,kp)))
         ! Perform tri-linear interpolation of V component
         ip=min(max(floor((pos(1)-this%cfg%xm(this%cfg%imino))/this%dx)+this%cfg%imino,this%cfg%imino_),this%cfg%imaxo_-1); wx1=min(max((pos(1)-this%cfg%xm(ip))*this%dxi,0.0_WP),1.0_WP); wx2=1.0_WP-wx1
         jp=min(max(floor((pos(2)-this%cfg%y (this%cfg%jmino))/this%dy)+this%cfg%jmino,this%cfg%jmino_),this%cfg%jmaxo_-1); wy1=min(max((pos(2)-this%cfg%y (jp))*this%dyi,0.0_WP),1.0_WP); wy2=1.0_WP-wy1
         kp=min(max(floor((pos(3)-this%cfg%zm(this%cfg%kmino))/this%dz)+this%cfg%kmino,this%cfg%kmino_),this%cfg%kmaxo_-1); wz1=min(max((pos(3)-this%cfg%zm(kp))*this%dzi,0.0_WP),1.0_WP); wz2=1.0_WP-wz1
         vel(2)=wz1*(wy1*(wx1*V(ip+1,jp+1,kp+1)+wx2*V(ip,jp+1,kp+1))+wy2*(wx1*V(ip+1,jp,kp+1)+wx2*V(ip,jp,kp+1)))+wz2*(wy1*(wx1*V(ip+1,jp+1,kp)+wx2*V(ip,jp+1,kp))+wy2*(wx1*V(ip+1,jp,kp)+wx2*V(ip,jp,kp)))
         ! Perform tri-linear interpolation of W component
         ip=min(max(floor((pos(1)-this%cfg%xm(this%cfg%imino))/this%dx)+this%cfg%imino,this%cfg%imino_),this%cfg%imaxo_-1); wx1=min(max((pos(1)-this%cfg%xm(ip))*this%dxi,0.0_WP),1.0_WP); wx2=1.0_WP-wx1
         jp=min(max(floor((pos(2)-this%cfg%ym(this%cfg%jmino))/this%dy)+this%cfg%jmino,this%cfg%jmino_),this%cfg%jmaxo_-1); wy1=min(max((pos(2)-this%cfg%ym(jp))*this%dyi,0.0_WP),1.0_WP); wy2=1.0_WP-wy1
         kp=min(max(floor((pos(3)-this%cfg%z (this%cfg%kmino))/this%dz)+this%cfg%kmino,this%cfg%kmino_),this%cfg%kmaxo_-1); wz1=min(max((pos(3)-this%cfg%z (kp))*this%dzi,0.0_WP),1.0_WP); wz2=1.0_WP-wz1
         vel(3)=wz1*(wy1*(wx1*W(ip+1,jp+1,kp+1)+wx2*W(ip,jp+1,kp+1))+wy2*(wx1*W(ip+1,jp,kp+1)+wx2*W(ip,jp,kp+1)))+wz2*(wy1*(wx1*W(ip+1,jp+1,kp)+wx2*W(ip,jp+1,kp))+wy2*(wx1*W(ip+1,jp,kp)+wx2*W(ip,jp,kp)))
      end function interp_velocity
   end subroutine SLstep
   
   
   !> Add SL RHS for all equations
   subroutine SLincrement(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer  :: i,j,k,n
      real(WP), dimension(:,:,:), allocatable :: FUX,FUY,FUZ,FVX,FVY,FVZ,FWX,FWY,FWZ
      real(WP) :: div
      
      ! Start rhs timer
      call this%trhs%start()
      
      ! Allocate mixture momentum fluxes
      allocate(FUX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate cell-centered mixture momentum fluxes from SL mass fluxes with extra cell on the left due to staggering
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               FUX(i,j,k)=0.25_WP*sum(this%SLFx(i:i+1,j,k,1:2))*sum(this%U(i:i+1,j,k))
               FVY(i,j,k)=0.25_WP*sum(this%SLFy(i,j:j+1,k,1:2))*sum(this%V(i,j:j+1,k))
               FWZ(i,j,k)=0.25_WP*sum(this%SLFz(i,j,k:k+1,1:2))*sum(this%W(i,j,k:k+1))
            end do
         end do
      end do
      
      ! Calculate edge-centered mixture momentum fluxes from SL mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FUY(i,j,k)=0.25_WP*sum(this%SLFy(i-1:i,j,k,1:2))*sum(this%U(i,j-1:j,k))
               FUZ(i,j,k)=0.25_WP*sum(this%SLFz(i-1:i,j,k,1:2))*sum(this%U(i,j,k-1:k))
               FVX(i,j,k)=0.25_WP*sum(this%SLFx(i,j-1:j,k,1:2))*sum(this%V(i-1:i,j,k))
               FVZ(i,j,k)=0.25_WP*sum(this%SLFz(i,j-1:j,k,1:2))*sum(this%V(i,j,k-1:k))
               FWX(i,j,k)=0.25_WP*sum(this%SLFx(i,j,k-1:k,1:2))*sum(this%W(i-1:i,j,k))
               FWY(i,j,k)=0.25_WP*sum(this%SLFy(i,j,k-1:k,1:2))*sum(this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Increment all conserved variables using terms built from SL mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Phasic mass and internal energy advection
               this%Q(i,j,k,1:4)=this%Q(i,j,k,1:4)+this%SLdt*(this%dxi*(this%SLFx(i+1,j,k,1:4)-this%SLFx(i,j,k,1:4))+this%dyi*(this%SLFy(i,j+1,k,1:4)-this%SLFy(i,j,k,1:4))+this%dzi*(this%SLFz(i,j,k+1,1:4)-this%SLFz(i,j,k,1:4)))
               ! Pressure dilatation term
               if (this%iSL(i,j,k).gt.0) then
                  div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
                  this%Q(i,j,k,3)=this%Q(i,j,k,3)-this%SLdt*(       this%VF(i,j,k))*this%PL(i,j,k)*div
                  this%Q(i,j,k,4)=this%Q(i,j,k,4)-this%SLdt*(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)*div
               end if
               ! Mixture momentum advection
               this%Q(i,j,k,5)=this%Q(i,j,k,5)+this%SLdt*(this%dxi*(FUX(i  ,j,k)-FUX(i-1,j,k))+this%dyi*(FUY(i,j+1,k)-FUY(i,j  ,k))+this%dzi*(FUZ(i,j,k+1)-FUZ(i,j,k  )))
               this%Q(i,j,k,6)=this%Q(i,j,k,6)+this%SLdt*(this%dxi*(FVX(i+1,j,k)-FVX(i  ,j,k))+this%dyi*(FVY(i,j  ,k)-FVY(i,j-1,k))+this%dzi*(FVZ(i,j,k+1)-FVZ(i,j,k  )))
               this%Q(i,j,k,7)=this%Q(i,j,k,7)+this%SLdt*(this%dxi*(FWX(i+1,j,k)-FWX(i  ,j,k))+this%dyi*(FWY(i,j+1,k)-FWY(i,j  ,k))+this%dzi*(FWZ(i,j,k  )-FWZ(i,j,k-1)))
            end do
         end do
      end do
      
      ! Deallocate all flux arrays
      deallocate(FUX,FUY,FUZ,FVX,FVY,FVZ,FWX,FWY,FWZ)
      
      ! Synchronize all Q fields
      do n=1,this%nQ; call this%cfg%sync(this%Q(:,:,:,n)); end do
      
      ! Stop rhs timer
      call this%trhs%stop()
      
   end subroutine SLincrement
   
   
   !> Obtain non-SL RHS for all equations except volume
   subroutine rhs(this,dQdt)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: dQdt  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nVAR)
      real(WP), dimension(:,:,:,:), allocatable :: Fx,Fy,Fz
      real(WP), dimension(:,:,:)  , allocatable :: FUX,FUY,FUZ,FVX,FVY,FVZ,FWX,FWY,FWZ
      integer :: i,j,k,n
      real(WP) :: w,div
      real(WP), parameter :: eps=1.0e-15_WP
      real(WP), dimension(-2: 0) :: wenop
      real(WP), dimension(-1:+1) :: wenom
      
      ! Start rhs timer
      call this%trhs%start()
      
      ! Zero out RHS
      dQdt=0.0_WP
      
      ! ================================================================ !
      ! ======================== INVISID FLUXES ======================== !
      ! ================================================================ !
      
      ! Allocate phasic mass, phasic energy, and mixture momentum fluxes
      allocate(Fx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:4)); Fx=0.0_WP
      allocate(Fy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:4)); Fy=0.0_WP
      allocate(Fz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:4)); Fz=0.0_WP
      allocate(FUX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Calculate standard (i.e., non-SL) fluxes in the appropriate phase
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (maxval(this%iSL(i-1:i,j,k)).eq.0) then
                  if (this%VF(i,j,k).gt.0.5_WP) then
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOL(i-1,j,k)-this%RHOL(i-2,j,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOL(i+1,j,k)-this%RHOL(i  ,j,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fx(i,j,k,1)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wenop*this%RHOL(i-2:i  ,j,k))-0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wenom*this%RHOL(i-1:i+1,j,k))
                     ! Centered phasic mass flux
                     !Fx(i,j,k,1)=-this%U(i,j,k)*0.5_WP*sum(this%RHOL(i-1:i,j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IL(i-1,j,k)-this%IL(i-2,j,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IL(i+1,j,k)-this%IL(i  ,j,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fx(i,j,k,3)=0.5_WP*(Fx(i,j,k,1)-abs(-Fx(i,j,k,1)))*sum(wenop*this%IL(i-2:i,j,k))+0.5_WP*(Fx(i,j,k,1)+abs(-Fx(i,j,k,1)))*sum(wenom*this%IL(i-1:i+1,j,k))
                     ! Centered phasic internal energy flux
                     !Fx(i,j,k,3)=Fx(i,j,k,1)*0.5_WP*sum(this%IL(i-1:i,j,k))
                  else
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOG(i-1,j,k)-this%RHOG(i-2,j,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOG(i+1,j,k)-this%RHOG(i  ,j,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fx(i,j,k,2)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wenop*this%RHOG(i-2:i,j,k))-0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wenom*this%RHOG(i-1:i+1,j,k))
                     ! Centered phasic mass flux
                     !Fx(i,j,k,2)=-this%U(i,j,k)*0.5_WP*sum(this%RHOG(i-1:i,j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IG(i-1,j,k)-this%IG(i-2,j,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IG(i+1,j,k)-this%IG(i  ,j,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fx(i,j,k,4)=0.5_WP*(Fx(i,j,k,2)-abs(-Fx(i,j,k,2)))*sum(wenop*this%IG(i-2:i,j,k))+0.5_WP*(Fx(i,j,k,2)+abs(-Fx(i,j,k,2)))*sum(wenom*this%IG(i-1:i+1,j,k))
                     ! Centered phasic internal energy flux
                     !Fx(i,j,k,4)=Fx(i,j,k,2)*0.5_WP*sum(this%IG(i-1:i,j,k))
                  end if
               end if
               if (maxval(this%iSL(i,j-1:j,k)).eq.0) then
                  if (this%VF(i,j,k).gt.0.5_WP) then
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOL(i,j-1,k)-this%RHOL(i,j-2,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOL(i,j+1,k)-this%RHOL(i,j  ,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fy(i,j,k,1)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wenop*this%RHOL(i,j-2:j,k))-0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wenom*this%RHOL(i,j-1:j+1,k))
                     ! Centered phasic mass flux
                     !Fy(i,j,k,1)=-this%V(i,j,k)*0.5_WP*sum(this%RHOL(i,j-1:j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IL(i,j-1,k)-this%IL(i,j-2,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IL(i,j+1,k)-this%IL(i,j  ,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fy(i,j,k,3)=0.5_WP*(Fy(i,j,k,1)-abs(-Fy(i,j,k,1)))*sum(wenop*this%IL(i,j-2:j,k))+0.5_WP*(Fy(i,j,k,1)+abs(-Fy(i,j,k,1)))*sum(wenom*this%IL(i,j-1:j+1,k))
                     ! Centered phasic internal energy flux
                     !Fy(i,j,k,3)=Fy(i,j,k,1)*0.5_WP*sum(this%IL(i,j-1:j,k))
                  else
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOG(i,j-1,k)-this%RHOG(i,j-2,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOG(i,j+1,k)-this%RHOG(i,j  ,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fy(i,j,k,2)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wenop*this%RHOG(i,j-2:j,k))-0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wenom*this%RHOG(i,j-1:j+1,k))
                     ! Centered phasic mass flux
                     !Fy(i,j,k,2)=-this%V(i,j,k)*0.5_WP*sum(this%RHOG(i,j-1:j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IG(i,j-1,k)-this%IG(i,j-2,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IG(i,j+1,k)-this%IG(i,j  ,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fy(i,j,k,4)=0.5_WP*(Fy(i,j,k,2)-abs(-Fy(i,j,k,2)))*sum(wenop*this%IG(i,j-2:j,k))+0.5_WP*(Fy(i,j,k,2)+abs(-Fy(i,j,k,2)))*sum(wenom*this%IG(i,j-1:j+1,k))
                     ! Centered phasic internal energy flux
                     !Fy(i,j,k,4)=Fy(i,j,k,2)*0.5_WP*sum(this%IG(i,j-1:j,k))
                  end if
               end if
               if (maxval(this%iSL(i,j,k-1:k)).eq.0) then
                  if (this%VF(i,j,k).gt.0.5_WP) then
                     ! WENO mass flux
                     w=weno_weight((abs(this%RHOL(i,j,k-1)-this%RHOL(i,j,k-2))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOL(i,j,k+1)-this%RHOL(i,j,k  ))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fz(i,j,k,1)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wenop*this%RHOL(i,j,k-2:k))-0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wenom*this%RHOL(i,j,k-1:k+1))
                     ! Centered mass flux
                     !Fz(i,j,k,1)=-this%W(i,j,k)*0.5_WP*sum(this%RHOL(i,j,k-1:k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IL(i,j,k-1)-this%IL(i,j,k-2))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IL(i,j,k+1)-this%IL(i,j,k  ))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fz(i,j,k,3)=0.5_WP*(Fz(i,j,k,1)-abs(-Fz(i,j,k,1)))*sum(wenop*this%IL(i,j,k-2:k))+0.5_WP*(Fz(i,j,k,1)+abs(-Fz(i,j,k,1)))*sum(wenom*this%IL(i,j,k-1:k+1))
                     ! Centered internal energy flux
                     !Fz(i,j,k,3)=Fz(i,j,k,1)*0.5_WP*sum(this%IL(i,j,k-1:k))
                  else
                     ! WENO mass flux
                     w=weno_weight((abs(this%RHOG(i,j,k-1)-this%RHOG(i,j,k-2))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOG(i,j,k+1)-this%RHOG(i,j,k  ))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fz(i,j,k,2)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wenop*this%RHOG(i,j,k-2:k))-0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wenom*this%RHOG(i,j,k-1:k+1))
                     ! Centered mass flux
                     !Fz(i,j,k,2)=-this%W(i,j,k)*0.5_WP*sum(this%RHOG(i,j,k-1:k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IG(i,j,k-1)-this%IG(i,j,k-2))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IG(i,j,k+1)-this%IG(i,j,k  ))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     Fz(i,j,k,4)=0.5_WP*(Fz(i,j,k,2)-abs(-Fz(i,j,k,2)))*sum(wenop*this%IG(i,j,k-2:k))+0.5_WP*(Fz(i,j,k,2)+abs(-Fz(i,j,k,2)))*sum(wenom*this%IG(i,j,k-1:k+1))
                     ! Centered internal energy flux
                     !Fz(i,j,k,4)=Fz(i,j,k,2)*0.5_WP*sum(this%IG(i,j,k-1:k))
                  end if
               end if
            end do
         end do
      end do
      
      ! Mass fluxes will be used to build momentum fluxes, they need to be extended by one cell on the left because of staggering
      call this%cfg%sync(Fx(:,:,:,1)); call this%cfg%sync(Fx(:,:,:,2)); if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) Fx(this%cfg%imin-1,:,:,1:2)=Fx(this%cfg%imin,:,:,1:2)
      call this%cfg%sync(Fy(:,:,:,1)); call this%cfg%sync(Fy(:,:,:,2)); if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) Fy(:,this%cfg%jmin-1,:,1:2)=Fy(:,this%cfg%jmin,:,1:2)
      call this%cfg%sync(Fz(:,:,:,1)); call this%cfg%sync(Fz(:,:,:,2)); if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) Fz(:,:,this%cfg%kmin-1,1:2)=Fz(:,:,this%cfg%kmin,1:2)
      
      ! Calculate cell-centered mixture momentum fluxes from standard mass fluxes with extra cell on the left due to staggering
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               FUX(i,j,k)=0.25_WP*sum(Fx(i:i+1,j,k,1:2))*sum(this%U(i:i+1,j,k))-this%VF(i,j,k)*this%PL(i,j,k)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)
               FVY(i,j,k)=0.25_WP*sum(Fy(i,j:j+1,k,1:2))*sum(this%V(i,j:j+1,k))-this%VF(i,j,k)*this%PL(i,j,k)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)
               FWZ(i,j,k)=0.25_WP*sum(Fz(i,j,k:k+1,1:2))*sum(this%W(i,j,k:k+1))-this%VF(i,j,k)*this%PL(i,j,k)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)
            end do
         end do
      end do
      
      ! Calculate edge-centered mixture momentum fluxes from standard mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FUY(i,j,k)=0.25_WP*sum(Fy(i-1:i,j,k,1:2))*sum(this%U(i,j-1:j,k))
               FUZ(i,j,k)=0.25_WP*sum(Fz(i-1:i,j,k,1:2))*sum(this%U(i,j,k-1:k))
               FVX(i,j,k)=0.25_WP*sum(Fx(i,j-1:j,k,1:2))*sum(this%V(i-1:i,j,k))
               FVZ(i,j,k)=0.25_WP*sum(Fz(i,j-1:j,k,1:2))*sum(this%V(i,j,k-1:k))
               FWX(i,j,k)=0.25_WP*sum(Fx(i,j,k-1:k,1:2))*sum(this%W(i-1:i,j,k))
               FWY(i,j,k)=0.25_WP*sum(Fy(i,j,k-1:k,1:2))*sum(this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Assemble time derivative for conserved variables using terms built from standard mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Phasic mass and phasic internal energy advection
               dQdt(i,j,k,1:4)=this%dxi*(Fx(i+1,j,k,1:4)-Fx(i,j,k,1:4))+this%dyi*(Fy(i,j+1,k,1:4)-Fy(i,j,k,1:4))+this%dzi*(Fz(i,j,k+1,1:4)-Fz(i,j,k,1:4))
               ! Pressure dilatation term
               if (this%iSL(i,j,k).eq.0) then
                  div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
                  dQdt(i,j,k,3)=dQdt(i,j,k,3)-(       this%VF(i,j,k))*this%PL(i,j,k)*div
                  dQdt(i,j,k,4)=dQdt(i,j,k,4)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)*div
               end if
               ! Mixture momentum advection and pressure stress
               dQdt(i,j,k,5)=this%dxi*(FUX(i  ,j,k)-FUX(i-1,j,k))+this%dyi*(FUY(i,j+1,k)-FUY(i,j  ,k))+this%dzi*(FUZ(i,j,k+1)-FUZ(i,j,k  ))
               dQdt(i,j,k,6)=this%dxi*(FVX(i+1,j,k)-FVX(i  ,j,k))+this%dyi*(FVY(i,j  ,k)-FVY(i,j-1,k))+this%dzi*(FVZ(i,j,k+1)-FVZ(i,j,k  ))
               dQdt(i,j,k,7)=this%dxi*(FWX(i+1,j,k)-FWX(i  ,j,k))+this%dyi*(FWY(i,j+1,k)-FWY(i,j  ,k))+this%dzi*(FWZ(i,j,k  )-FWZ(i,j,k-1))
            end do
         end do
      end do
      
      ! ================================================================ !
      ! ======================== VISCOUS FLUXES ======================== !
      ! ================================================================ !
      
      
      
      ! Deallocate all flux arrays
      deallocate(Fx,Fy,Fz,FUX,FUY,FUZ,FVX,FVY,FVZ,FWX,FWY,FWZ)
      
      ! Synchronize all dQdt fields
      do n=1,this%nQ; call this%cfg%sync(dQdt(:,:,:,n)); end do
      
      ! Stop rhs timer
      call this%trhs%stop()
      
   contains
      !> WENO switch function
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP ! Switching parameter
         real(WP), parameter :: delta=0.01_WP  ! Switching thickness
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight
   end subroutine rhs
   
   
   !> Build the interface from volume moments using PLICnet
   subroutine build_interface(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      
      ! Start timer
      call this%tplic%start()
      
      ! Perfrom PLICnet reconstruction of the interface
      plicnet_reconstruct: block
         use mathtools, only: normalize
         use plicnet,   only: get_normal,reflect_moments
         integer(IRL_SignedIndex_t) :: i,j,k
         integer :: ind,ii,jj,kk
         real(IRL_double), dimension(0:2) :: normal
         real(IRL_double), dimension(0:188) :: moments
         integer :: direction
         logical :: flip
         real(IRL_double) :: m000,m100,m010,m001
         real(IRL_double), dimension(0:2) :: center
         real(IRL_double) :: initial_dist
         type(RectCub_type) :: cell
         ! Get a cell
         call new(cell)
         ! Traverse domain and reconstruct interface
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            ! Handle full cells differently
            if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
               call setNumberOfPlanes(this%PLIC(i,j,k),1)
               call setPlane(this%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               cycle
            end if
            ! Liquid-gas symmetry
            flip=.false.
            if (this%VF(i,j,k).ge.0.5_WP) flip=.true.
            m000=0; m100=0; m010=0; m001=0
            ! Construct neighborhood of volume moments
            if (flip) then
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=1.0_WP-this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%BG(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%BG(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%BG(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%BL(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%BL(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%BL(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  ! Calculate geometric moments of neighborhood
                  m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
               end do; end do; end do
            else
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%BL(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%BL(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%BL(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%BG(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%BG(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%BG(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  ! Calculate geometric moments of neighborhood
                  m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
               end do; end do; end do
            end if
            ! Calculate geometric center of neighborhood
            center=[m100,m010,m001]/m000
            ! Symmetry about Cartesian planes
            call reflect_moments(moments,center,direction)
            ! Get PLIC normal vector from neural network
            call get_normal(moments,normal)
            normal=normalize(normal)
            ! Rotate normal vector to original octant
            if (direction.eq.1) then
               normal(0)=-normal(0)
            else if (direction.eq.2) then
               normal(1)=-normal(1)
            else if (direction.eq.3) then
               normal(2)=-normal(2)
            else if (direction.eq.4) then
               normal(0)=-normal(0)
               normal(1)=-normal(1)
            else if (direction.eq.5) then
               normal(0)=-normal(0)
               normal(2)=-normal(2)
            else if (direction.eq.6) then
               normal(1)=-normal(1)
               normal(2)=-normal(2)
            else if (direction.eq.7) then
               normal(0)=-normal(0)
               normal(1)=-normal(1)
               normal(2)=-normal(2)
            end if
            if (.not.flip) then
               normal(0)=-normal(0)
               normal(1)=-normal(1)
               normal(2)=-normal(2)
            end if
            ! Locate PLIC plane in cell
            call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
            initial_dist=dot_product(normal,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
            call setNumberOfPlanes(this%PLIC(i,j,k),1)
            call setPlane(this%PLIC(i,j,k),0,normal,initial_dist)
            call matchVolumeFraction(cell,this%VF(i,j,k),this%PLIC(i,j,k))
         end do; end do; end do
      end block plicnet_reconstruct
      
      ! Synchronize PLIC interface across boundaries
      synchronize_interface: block
         integer :: i,j,k,ni
         real(WP), dimension(1:4) :: plane
         integer , dimension(2,3) :: send_range,recv_range
         ! Synchronize in x
         if (this%cfg%nx.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               call copy(this%PLIC(i,j,k),this%PLIC(this%cfg%imin,j,k))
            end do; end do; end do
         else
            ! Send minus
            send_range(1:2,1)=[this%cfg%imin_   ,this%cfg%imin_ +this%cfg%no-1]
            send_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            send_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            recv_range(1:2,1)=[this%cfg%imax_ +1,this%cfg%imaxo_              ]
            recv_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            recv_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            call sync_side(send_range,recv_range,0,-1)
            ! Send plus
            send_range(1:2,1)=[this%cfg%imax_ -this%cfg%no+1,this%cfg%imax_   ]
            send_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            send_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imin_ -1]
            recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            call sync_side(send_range,recv_range,0,+1)
         end if
         ! Synchronize in y
         if (this%cfg%ny.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               call copy(this%PLIC(i,j,k),this%PLIC(i,this%cfg%jmin,k))
            end do; end do; end do
         else
            ! Send minus side
            send_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            send_range(1:2,2)=[this%cfg%jmin_   ,this%cfg%jmin_ +this%cfg%no-1]
            send_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            recv_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            recv_range(1:2,2)=[this%cfg%jmax_ +1,this%cfg%jmaxo_              ]
            recv_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
            call sync_side(send_range,recv_range,1,-1)
            ! Send plus side
            send_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            send_range(1:2,2)=[this%cfg%jmax_ -this%cfg%no+1,this%cfg%jmax_   ]
            send_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmin_ -1]
            recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
            call sync_side(send_range,recv_range,1,+1)
         end if
         ! Synchronize in z
         if (this%cfg%nz.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               call copy(this%PLIC(i,j,k),this%PLIC(i,j,this%cfg%kmin))
            end do; end do; end do
         else
            ! Send minus side
            send_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            send_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            send_range(1:2,3)=[this%cfg%kmin_   ,this%cfg%kmin_ +this%cfg%no-1]
            recv_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
            recv_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
            recv_range(1:2,3)=[this%cfg%kmax_ +1,this%cfg%kmaxo_              ]
            call sync_side(send_range,recv_range,2,-1)
            ! Send plus side
            send_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            send_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            send_range(1:2,3)=[this%cfg%kmax_ -this%cfg%no+1,this%cfg%kmax_   ]
            recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
            recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
            recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmin_ -1]
            call sync_side(send_range,recv_range,2,+1)
         end if
         ! Fix plane position if we are periodic in x
         if (this%cfg%xper.and.this%cfg%iproc.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino,this%cfg%imin-1
               plane=getPlane(this%PLIC(i,j,k),0)
               plane(4)=plane(4)-plane(1)*this%cfg%xL
               call setPlane(this%PLIC(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         if (this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imax+1,this%cfg%imaxo
               plane=getPlane(this%PLIC(i,j,k),0)
               plane(4)=plane(4)+plane(1)*this%cfg%xL
               call setPlane(this%PLIC(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         ! Fix plane position if we are periodic in y
         if (this%cfg%yper.and.this%cfg%jproc.eq.1) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino,this%cfg%jmin-1; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%PLIC(i,j,k),0)
               plane(4)=plane(4)-plane(2)*this%cfg%yL
               call setPlane(this%PLIC(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         if (this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
            do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmax+1,this%cfg%jmaxo; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%PLIC(i,j,k),0)
               plane(4)=plane(4)+plane(2)*this%cfg%yL
               call setPlane(this%PLIC(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         ! Fix plane position if we are periodic in z
         if (this%cfg%zper.and.this%cfg%kproc.eq.1) then
            do k=this%cfg%kmino,this%cfg%kmin-1; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%PLIC(i,j,k),0)
               plane(4)=plane(4)-plane(3)*this%cfg%zL
               call setPlane(this%PLIC(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
         if (this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
            do k=this%cfg%kmax+1,this%cfg%kmaxo; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
               plane=getPlane(this%PLIC(i,j,k),0)
               plane(4)=plane(4)+plane(3)*this%cfg%zL
               call setPlane(this%PLIC(i,j,k),0,plane(1:3),plane(4))
            end do; end do; end do
         end if
      end block synchronize_interface
      
      ! Reset volume moments to match with reconstructed interface
      reset_moments: block
         integer :: i,j,k
         type(RectCub_type) :: cell
         type(SepVM_type) :: separated_volume_moments
         call new(cell)
         call new(separated_volume_moments)
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! Form the grid cell
            call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
            ! Cut it by the current interface
            call getNormMoments(cell,this%PLIC(i,j,k),separated_volume_moments)
            ! Recover relevant moments
            this%VF  (i,j,k)=getVolumePtr(separated_volume_moments,0)/this%cfg%vol(i,j,k)
            this%BL(:,i,j,k)= getCentroid(separated_volume_moments,0)
            this%BG(:,i,j,k)= getCentroid(separated_volume_moments,1)
            ! Clean up
            if (this%VF(i,j,k).lt.VFlo) then
               this%VF  (i,j,k)=0.0_WP
               this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end if
            if (this%VF(i,j,k).gt.VFhi) then
               this%VF  (i,j,k)=1.0_WP
               this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end if
         end do; end do; end do
      end block reset_moments
      
      ! Polygonize interface
      polygonize_interface: block
         integer :: i,j,k,n
         type(RectCub_type) :: cell
         real(WP), dimension(1:3,1:4) :: vert
         real(WP), dimension(1:3) :: norm
         ! Create a cell object
         call new(cell)
         ! Loop over full domain and form polygon
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! Zero out the polygon
            call zeroPolygon(this%interface_polygon(i,j,k))
            ! Create polygons for cells with interfaces, zero for those without
            if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) then
               call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               call getPoly(cell,this%PLIC(i,j,k),0,this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
         ! Find inferface between filled and empty cells on x-face
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_+1,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.and.this%VF(i-1,j,k).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i-1,j,k).lt.VFlo) then
               norm=[sign(1.0_WP,0.5_WP-this%VF(i,j,k)),0.0_WP,0.0_WP]
               call setNumberOfPlanes(this%PLIC(i,j,k),1)
               call setPlane(this%PLIC(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%x(i))
               vert(:,1)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k  )]
               vert(:,2)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k  )]
               vert(:,3)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k+1)]
               vert(:,4)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k+1)]
               call construct(this%interface_polygon(i,j,k),4,vert)
               call setPlaneOfExistence(this%interface_polygon(i,j,k),getPlane(this%PLIC(i,j,k),0))
               if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
         ! Find inferface between filled and empty cells on y-face
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_+1,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.and.this%VF(i,j-1,k).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i,j-1,k).lt.VFlo) then
               norm=[0.0_WP,sign(1.0_WP,0.5_WP-this%VF(i,j,k)),0.0_WP]
               call setNumberOfPlanes(this%PLIC(i,j,k),1)
               call setPlane(this%PLIC(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%y(j))
               vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k  )]
               vert(:,2)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k+1)]
               vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k+1)]
               vert(:,4)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k  )]
               call construct(this%interface_polygon(i,j,k),4,vert)
               call setPlaneOfExistence(this%interface_polygon(i,j,k),getPlane(this%PLIC(i,j,k),0))
               if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
         ! Find inferface between filled and empty cells on z-face
         do k=this%cfg%kmino_+1,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).lt.VFlo.and.this%VF(i,j,k-1).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i,j,k-1).lt.VFlo) then
               norm=[0.0_WP,0.0_WP,sign(1.0_WP,0.5_WP-this%VF(i,j,k))]
               call setNumberOfPlanes(this%PLIC(i,j,k),1)
               call setPlane(this%PLIC(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%z(k))
               vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k)]
               vert(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k)]
               vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k)]
               vert(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k)]
               call construct(this%interface_polygon(i,j,k),4,vert)
               call setPlaneOfExistence(this%interface_polygon(i,j,k),getPlane(this%PLIC(i,j,k),0))
               if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(i,j,k))
            end if
         end do; end do; end do
      end block polygonize_interface
      
      ! Stop timer
      call this%tplic%stop()
      
   contains
      
      !> Private procedure to perform communication across one boundary
      subroutine sync_side(a_send_range,a_recv_range,a_dimension,a_direction)
         implicit none
         integer, dimension(2,3), intent(in) :: a_send_range
         integer, dimension(2,3), intent(in) :: a_recv_range
         integer, intent(in) :: a_dimension
         integer, intent(in) :: a_direction
         integer :: i,j,k
         logical :: something_received
         ! Pack the buffer
         call resetBufferPointer(this%send_byte_buffer)
         call setSize(this%send_byte_buffer,int(0,8))
         do k=a_send_range(1,3),a_send_range(2,3); do j=a_send_range(1,2),a_send_range(2,2); do i=a_send_range(1,1),a_send_range(2,1)
            call serializeAndPack(this%PLIC(i,j,k),this%send_byte_buffer)
         end do; end do; end do
         ! Communicate
         call sync_ByteBuffer(this%send_byte_buffer,a_dimension,a_direction,this%recv_byte_buffer,something_received)
         ! If something was received, unpack it: traversal order is important and must be aligned with how the sent data was packed
         if (something_received) then
            call resetBufferPointer(this%recv_byte_buffer)
            do k=a_recv_range(1,3),a_recv_range(2,3); do j=a_recv_range(1,2),a_recv_range(2,2); do i=a_recv_range(1,1),a_recv_range(2,1)
               call unpackAndStore(this%PLIC(i,j,k),this%recv_byte_buffer)
            end do; end do; end do
         end if
      end subroutine sync_side
      
      !> Private procedure to communicate a package of bytes across one boundary
      subroutine sync_ByteBuffer(a_send_buffer,a_dimension,a_direction,a_receive_buffer,a_received_something)
         use mpi_f08
         implicit none
         type(ByteBuffer_type), intent(inout)  :: a_send_buffer !< Inout needed because it is preallocated
         integer, intent(in) :: a_dimension  !< Should be 0/1/2 for x/y/z
         integer, intent(in) :: a_direction  !< Should be -1 for left or +1 for right
         type(ByteBuffer_type), intent(inout) :: a_receive_buffer !< Inout needed because it is preallocated
         logical, intent(out) :: a_received_something
         type(MPI_Status) :: status
         integer :: isrc,idst,ierr
         integer(IRL_LargeOffsetIndex_t) :: my_size
         integer(IRL_LargeOffsetIndex_t) :: incoming_size
         integer :: my_size_small,incoming_size_small
         ! Figure out source and destination
         call MPI_CART_SHIFT(this%cfg%comm,a_dimension,a_direction,isrc,idst,ierr)
         ! Communicate sizes so that each processor knows what to expect in main communication
         my_size=getSize(a_send_buffer)
         call MPI_SENDRECV(my_size,1,MPI_INTEGER8,idst,0,incoming_size,1,MPI_INTEGER8,isrc,0,this%cfg%comm,status,ierr)
         ! Set size of recv buffer to appropriate size and perform send-receive
         if (isrc.ne.MPI_PROC_NULL) then
            a_received_something=.true.
            call setSize(a_receive_buffer,incoming_size)
         else
            a_received_something=.false.
            incoming_size=0
            call setSize(a_receive_buffer,int(1,8))
         end if
         ! Convert integers
         my_size_small=int(my_size,4)
         incoming_size_small=int(incoming_size,4)
         call MPI_SENDRECV(dataPtr(a_send_buffer),my_size_small,MPI_BYTE,idst,0,dataPtr(a_receive_buffer),incoming_size_small,MPI_BYTE,isrc,0,this%cfg%comm,status,ierr)
      end subroutine sync_ByteBuffer
      
   end subroutine build_interface
   
   
   !> Calculate all primitive and mixture variables from updated conserved variables
   subroutine get_primitive(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP) :: CL,CG
      integer :: i,j,k
      
      ! Get mixture velocity
      call this%get_velocity()
      
      ! Get other mixture and phasic primitive variables
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         
         ! Get liquid primitive variables
         if (this%VF(i,j,k).ge.VFlo) then
            this%RHOL(i,j,k)=this%Q (i,j,k,1)/this%VF(i,j,k)
            this%IL  (i,j,k)=this%Q (i,j,k,3)/this%Q (i,j,k,1)
            this%PL  (i,j,k)=this%getPL(this%RHOL(i,j,k),this%IL(i,j,k))
            CL              =this%getCL(this%RHOL(i,j,k),this%PL(i,j,k))
         else
            this%RHOL(i,j,k)=0.0_WP
            this%IL  (i,j,k)=0.0_WP
            this%PL  (i,j,k)=0.0_WP
            CL              =0.0_WP
         end if
         
         ! Get gas primitive variables
         if (this%VF(i,j,k).le.VFhi) then
            this%RHOG(i,j,k)=this%Q(i,j,k,2)/(1.0_WP-this%VF(i,j,k))
            this%IG  (i,j,k)=this%Q(i,j,k,4)/        this%Q (i,j,k,2)
            this%PG  (i,j,k)=this%getPG(this%RHOG(i,j,k),this%IG(i,j,k))
            CG              =this%getCG(this%RHOG(i,j,k),this%PG(i,j,k))
         else
            this%RHOG(i,j,k)=0.0_WP
            this%IG  (i,j,k)=0.0_WP
            this%PG  (i,j,k)=0.0_WP
            CG              =0.0_WP
         end if
         
         ! Rebuild conserved quantities to ensure consistency
         this%Q(i,j,k,1)=(       this%VF(i,j,k))*this%RHOL(i,j,k); this%Q(i,j,k,3)=this%Q(i,j,k,1)*this%IL(i,j,k)
         this%Q(i,j,k,2)=(1.0_WP-this%VF(i,j,k))*this%RHOG(i,j,k); this%Q(i,j,k,4)=this%Q(i,j,k,2)*this%IG(i,j,k)
         
         ! Get mixture speed of sound
         this%C(i,j,k)=sqrt((this%Q(i,j,k,1)*CL**2+this%Q(i,j,k,2)*CG**2)/sum(this%Q(i,j,k,1:2)))
         
      end do; end do; end do
      
      ! Finally rebuild momentum for consistency
      call this%get_momentum()
      
      ! Get phasic temperatures
      if (associated(this%getTL)) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).ge.VFlo) then
               this%TL(i,j,k)=this%getTL(this%RHOL(i,j,k),this%PL(i,j,k))
            else
               this%TL(i,j,k)=0.0_WP
            end if
         end do; end do; end do
      end if
      if (associated(this%getTG)) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            if (this%VF(i,j,k).le.VFhi) then
               this%TG(i,j,k)=this%getTG(this%RHOG(i,j,k),this%PG(i,j,k))   
            else
               this%TG(i,j,k)=0.0_WP
            end if
         end do; end do; end do
      end if
      
   end subroutine get_primitive
   
   
   !> Calculate velocity from momentum and mixture density, which is not recalculated
   subroutine get_velocity(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate mixture velocity as far as possible
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%U(i,j,k)=2.0_WP*this%Q(i,j,k,5)/sum(this%Q(i-1:i,j,k,1:2))
               this%V(i,j,k)=2.0_WP*this%Q(i,j,k,6)/sum(this%Q(i,j-1:j,k,1:2))
               this%W(i,j,k)=2.0_WP*this%Q(i,j,k,7)/sum(this%Q(i,j,k-1:k,1:2))
            end do
         end do
      end do
      ! Sync velocity
      call this%cfg%sync(this%U)
      call this%cfg%sync(this%V)
      call this%cfg%sync(this%W)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         this%U(this%cfg%imino,:,:)=this%Q(this%cfg%imino,:,:,5)/(this%Q(this%cfg%imino,:,:,1)+this%Q(this%cfg%imino,:,:,2))
         this%V(this%cfg%imino,:,:)=this%Q(this%cfg%imino,:,:,6)/(this%Q(this%cfg%imino,:,:,1)+this%Q(this%cfg%imino,:,:,2))
         this%W(this%cfg%imino,:,:)=this%Q(this%cfg%imino,:,:,7)/(this%Q(this%cfg%imino,:,:,1)+this%Q(this%cfg%imino,:,:,2))
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         this%U(:,this%cfg%jmino,:)=this%Q(:,this%cfg%jmino,:,5)/(this%Q(:,this%cfg%jmino,:,1)+this%Q(:,this%cfg%jmino,:,2))
         this%V(:,this%cfg%jmino,:)=this%Q(:,this%cfg%jmino,:,6)/(this%Q(:,this%cfg%jmino,:,1)+this%Q(:,this%cfg%jmino,:,2))
         this%W(:,this%cfg%jmino,:)=this%Q(:,this%cfg%jmino,:,7)/(this%Q(:,this%cfg%jmino,:,1)+this%Q(:,this%cfg%jmino,:,2))
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         this%U(:,:,this%cfg%kmino)=this%Q(:,:,this%cfg%kmino,5)/(this%Q(:,:,this%cfg%kmino,1)+this%Q(:,:,this%cfg%kmino,2))
         this%V(:,:,this%cfg%kmino)=this%Q(:,:,this%cfg%kmino,6)/(this%Q(:,:,this%cfg%kmino,1)+this%Q(:,:,this%cfg%kmino,2))
         this%W(:,:,this%cfg%kmino)=this%Q(:,:,this%cfg%kmino,7)/(this%Q(:,:,this%cfg%kmino,1)+this%Q(:,:,this%cfg%kmino,2))
      end if
   end subroutine get_velocity
   
   
   !> Calculate mixture kinetic energy per unit mass from pre-calculated mixture velocity
   !> Need to redo this better
   subroutine get_ke(this,KE)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: KE !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               KE(i,j,k)=0.5_WP*sum(this%U(i:i+1,j,k)**2+this%V(i,j:j+1,k)**2+this%W(i,j,k:k+1)**2)
            end do
         end do
      end do
      call this%cfg%sync(KE)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) KE(this%cfg%imaxo,:,:)=KE(this%cfg%imaxo-1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) KE(:,this%cfg%jmaxo,:)=KE(:,this%cfg%jmaxo-1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) KE(:,:,this%cfg%kmaxo)=KE(:,:,this%cfg%kmaxo-1)
   end subroutine get_ke
   
   
   !> Calculate momentum from velocity and mixture density
   subroutine get_momentum(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate momentum as far as possible
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%Q(i,j,k,5)=0.5_WP*sum(this%Q(i-1:i,j,k,1:2))*this%U(i,j,k)
               this%Q(i,j,k,6)=0.5_WP*sum(this%Q(i,j-1:j,k,1:2))*this%V(i,j,k)
               this%Q(i,j,k,7)=0.5_WP*sum(this%Q(i,j,k-1:k,1:2))*this%W(i,j,k)
            end do
         end do
      end do
      ! Sync momentum
      call this%cfg%sync(this%Q(:,:,:,5))
      call this%cfg%sync(this%Q(:,:,:,6))
      call this%cfg%sync(this%Q(:,:,:,7))
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         this%Q(this%cfg%imino,:,:,5)=(this%Q(this%cfg%imino,:,:,1)+this%Q(this%cfg%imino,:,:,2))*this%U(this%cfg%imino,:,:)
         this%Q(this%cfg%imino,:,:,6)=(this%Q(this%cfg%imino,:,:,1)+this%Q(this%cfg%imino,:,:,2))*this%V(this%cfg%imino,:,:)
         this%Q(this%cfg%imino,:,:,7)=(this%Q(this%cfg%imino,:,:,1)+this%Q(this%cfg%imino,:,:,2))*this%W(this%cfg%imino,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         this%Q(:,this%cfg%jmino,:,5)=(this%Q(:,this%cfg%jmino,:,1)+this%Q(:,this%cfg%jmino,:,2))*this%U(:,this%cfg%jmino,:)
         this%Q(:,this%cfg%jmino,:,6)=(this%Q(:,this%cfg%jmino,:,1)+this%Q(:,this%cfg%jmino,:,2))*this%V(:,this%cfg%jmino,:)
         this%Q(:,this%cfg%jmino,:,7)=(this%Q(:,this%cfg%jmino,:,1)+this%Q(:,this%cfg%jmino,:,2))*this%W(:,this%cfg%jmino,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         this%Q(:,:,this%cfg%kmino,5)=(this%Q(:,:,this%cfg%kmino,1)+this%Q(:,:,this%cfg%kmino,2))*this%U(:,:,this%cfg%kmino)
         this%Q(:,:,this%cfg%kmino,6)=(this%Q(:,:,this%cfg%kmino,1)+this%Q(:,:,this%cfg%kmino,2))*this%V(:,:,this%cfg%kmino)
         this%Q(:,:,this%cfg%kmino,7)=(this%Q(:,:,this%cfg%kmino,1)+this%Q(:,:,this%cfg%kmino,2))*this%W(:,:,this%cfg%kmino)
      end if
   end subroutine get_momentum
   
   
   !> Interpolate mixture velocity to cell-center, including overlap and ghosts
   subroutine interp_vel(this,Ui,Vi,Wi)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate interpolated velocity as far as possible
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               Ui(i,j,k)=0.5_WP*(this%U(i,j,k)+this%U(i+1,j,k))
               Vi(i,j,k)=0.5_WP*(this%V(i,j,k)+this%V(i,j+1,k))
               Wi(i,j,k)=0.5_WP*(this%W(i,j,k)+this%W(i,j,k+1))
            end do
         end do
      end do
      ! Sync interpolated velocity
      call this%cfg%sync(Ui)
      call this%cfg%sync(Vi)
      call this%cfg%sync(Wi)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         Ui(this%cfg%imaxo,:,:)=this%U(this%cfg%imaxo,:,:)
         Vi(this%cfg%imaxo,:,:)=this%V(this%cfg%imaxo,:,:)
         Wi(this%cfg%imaxo,:,:)=this%W(this%cfg%imaxo,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         Ui(:,this%cfg%jmaxo,:)=this%U(:,this%cfg%jmaxo,:)
         Vi(:,this%cfg%jmaxo,:)=this%V(:,this%cfg%jmaxo,:)
         Wi(:,this%cfg%jmaxo,:)=this%W(:,this%cfg%jmaxo,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         Ui(:,:,this%cfg%kmaxo)=this%U(:,:,this%cfg%kmaxo)
         Vi(:,:,this%cfg%kmaxo)=this%V(:,:,this%cfg%kmaxo)
         Wi(:,:,this%cfg%kmaxo)=this%W(:,:,this%cfg%kmaxo)
      end if
   end subroutine interp_vel
   
   
   !> Apply user-provided relaxation model to interfacial cells
   subroutine apply_relax(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k
      type(RectCub_type) :: cell
      type(SepVM_type)   :: separated_volume_moments
      ! If no relaxation model was provided, return
      if (.not.associated(this%relax)) return
      ! Allocate IRL objects
      call new(cell)
      call new(separated_volume_moments)
      ! Loop over the full domain
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         ! Ignore single-phase cells
         if (this%VF(i,j,k).le.VFlo.or.this%VF(i,j,k).ge.VFhi) cycle
         ! Apply user-provided relaxation model
         call this%relax(this%VF(i,j,k),this%Q(i,j,k,:))
         ! Adjust PLIC interface location
         call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
         call matchVolumeFraction(cell,this%VF(i,j,k),this%PLIC(i,j,k))
         ! Adjust volume moments to ensure consistency
         call getNormMoments(cell,this%PLIC(i,j,k),separated_volume_moments)
         this%VF  (i,j,k)=getVolumePtr(separated_volume_moments,0)/this%cfg%vol(i,j,k)
         this%BL(:,i,j,k)= getCentroid(separated_volume_moments,0)
         this%BG(:,i,j,k)= getCentroid(separated_volume_moments,1)
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF  (i,j,k)=0.0_WP
            this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
         end if
         if (this%VF(i,j,k).gt.VFhi) then
            this%VF  (i,j,k)=1.0_WP
            this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
         end if
      end do; end do; end do
   end subroutine apply_relax
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      integer :: ierr
      real(WP) :: maxvisc,maxC
      ! Compute convective+acoustic CFLs
      this%CFLc_x=maxval(abs(this%U)+abs(this%C))*dt*this%dxi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_y=maxval(abs(this%V)+abs(this%C))*dt*this%dyi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_z=maxval(abs(this%W)+abs(this%C))*dt*this%dzi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      ! Compute acoustic CFLs
      maxC=maxval(this%C); call MPI_ALLREDUCE(MPI_IN_PLACE,maxC,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLa_x=maxC*dt*this%dxi
      this%CFLa_y=maxC*dt*this%dyi
      this%CFLa_z=maxC*dt*this%dzi
      ! Compute viscous CFLs
      maxvisc=maxval((this%VF*this%VISCL+(1.0_WP-this%VF)*this%VISCG)/(this%Q(:,:,:,1)+this%Q(:,:,:,2))); call MPI_ALLREDUCE(MPI_IN_PLACE,maxvisc,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLv_x=4.0_WP*maxvisc*dt*this%dxi**2
      this%CFLv_y=4.0_WP*maxvisc*dt*this%dyi**2
      this%CFLv_z=4.0_WP*maxvisc*dt*this%dzi**2
      ! Return the maximum overall CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,&
      &       this%CFLa_x,this%CFLa_y,this%CFLa_z,&
      &       this%CFLv_x,this%CFLv_y,this%CFLv_z)
   end subroutine get_cfl
   
   
   !> Calculate info about our fields
   subroutine get_info(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: n,i,j,k,ierr
      real(WP), dimension(:,:,:), allocatable :: tmp
      
      ! Compute integrals and extrema of conserved variables
      call this%cfg%integrate(this%VF,integral=this%VFint)
      do n=1,this%nQ
         call this%cfg%integrate(this%Q(:,:,:,n),integral=this%Qint(n))
      end do
      this%VFmin=+huge(1.0_WP)
      this%VFmax=-huge(1.0_WP)
      this%Qmin=+huge(1.0_WP)
      this%Qmax=-huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         this%VFmin=min(this%VFmin,this%VF(i,j,k))
         this%VFmax=max(this%VFmax,this%VF(i,j,k))
         do n=1,this%nQ
            this%Qmin(n)=min(this%Qmin(n),this%Q(i,j,k,n))
            this%Qmax(n)=max(this%Qmax(n),this%Q(i,j,k,n))
         end do
      end do; end do; end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Qmin,this%nQ,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Qmax,this%nQ,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Also compute integral of phasic KE and phasic entropy
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      call this%get_ke(tmp); tmp=this%Q(:,:,:,1)*tmp; call this%cfg%integrate(tmp,integral=this%RHOKLint)
      call this%get_ke(tmp); tmp=this%Q(:,:,:,2)*tmp; call this%cfg%integrate(tmp,integral=this%RHOKGint)
      this%RHOSLint=0.0_WP
      if (associated(this%getSL)) then
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%VF(i,j,k).ge.VFlo) then
               tmp(i,j,k)=this%Q(i,j,k,1)*this%getSL(this%RHOL(i,j,k),this%PL(i,j,k))
            else
               tmp(i,j,k)=0.0_WP
            end if
         end do; end do; end do
         call this%cfg%integrate(tmp,integral=this%RHOSLint)
      end if
      this%RHOSGint=0.0_WP
      if (associated(this%getSG)) then
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%VF(i,j,k).le.VFhi) then
               tmp(i,j,k)=this%Q(i,j,k,2)*this%getSG(this%RHOG(i,j,k),this%PG(i,j,k))
            else
               tmp(i,j,k)=0.0_WP
            end if
         end do; end do; end do
         call this%cfg%integrate(tmp,integral=this%RHOSGint)
      end if
      deallocate(tmp)
      
      ! Calculate extrema of primitive fields
      this%RHOLmin=+huge(1.0_WP); this%RHOLmax=-huge(1.0_WP); this%RHOGmin=+huge(1.0_WP); this%RHOGmax=-huge(1.0_WP)
      this%ILmin  =+huge(1.0_WP); this%ILmax  =-huge(1.0_WP); this%IGmin  =+huge(1.0_WP); this%IGmax  =-huge(1.0_WP)
      this%PLmin  =+huge(1.0_WP); this%PLmax  =-huge(1.0_WP); this%PGmin  =+huge(1.0_WP); this%PGmax  =-huge(1.0_WP)
      this%TLmin  =+huge(1.0_WP); this%TLmax  =-huge(1.0_WP); this%TGmin  =+huge(1.0_WP); this%TGmax  =-huge(1.0_WP)
      this%Umax=0.0_WP; this%Vmax=0.0_WP; this%Wmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         if (this%VF(i,j,k).gt.VFlo) then
            this%RHOLmin=min(this%RHOLmin,this%RHOL(i,j,k)); this%RHOLmax=max(this%RHOLmax,this%RHOL(i,j,k))
            this%ILmin  =min(this%ILmin  ,this%IL  (i,j,k)); this%ILmax  =max(this%ILmax  ,this%IL  (i,j,k))
            this%PLmin  =min(this%PLmin  ,this%PL  (i,j,k)); this%PLmax  =max(this%PLmax  ,this%PL  (i,j,k))
            this%TLmin  =min(this%TLmin  ,this%TL  (i,j,k)); this%TLmax  =max(this%TLmax  ,this%TL  (i,j,k))
         end if
         if (this%VF(i,j,k).lt.VFhi) then
            this%RHOGmin=min(this%RHOGmin,this%RHOG(i,j,k)); this%RHOGmax=max(this%RHOGmax,this%RHOG(i,j,k))
            this%IGmin  =min(this%IGmin  ,this%IG  (i,j,k)); this%IGmax  =max(this%IGmax  ,this%IG  (i,j,k))
            this%PGmin  =min(this%PGmin  ,this%PG  (i,j,k)); this%PGmax  =max(this%PGmax  ,this%PG  (i,j,k))
            this%TGmin  =min(this%TGmin  ,this%TG  (i,j,k)); this%TGmax  =max(this%TGmax  ,this%TG  (i,j,k))
         end if
         this%Umax=max(this%Umax,abs(this%U(i,j,k)))
         this%Vmax=max(this%Vmax,abs(this%V(i,j,k)))
         this%Wmax=max(this%Wmax,abs(this%W(i,j,k)))
      end do; end do; end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOLmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOLmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOGmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOGmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%ILmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%ILmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%IGmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%IGmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%PLmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%PLmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%PGmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%PGmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%TLmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%TLmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%TGmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%TGmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax   ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_info
   
   
   !> Update a surfmesh object from our current polygons
   subroutine update_surfmesh(this,smesh)
      use surfmesh_class, only: surfmesh
      implicit none
      class(mpcomp),   intent(inout) :: this
      class(surfmesh), intent(inout) :: smesh
      integer :: i,j,k,n,shape,nv,np
      real(WP), dimension(3) :: tmp_vert
      ! Reset surface mesh storage
      call smesh%reset()
      ! First pass to count how many vertices and polygons are inside our processor
      nv=0; np=0
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         shape=getNumberOfVertices(this%interface_polygon(i,j,k))
         if (shape.gt.0) then
            nv=nv+shape
            np=np+1
         end if
      end do; end do; end do
      ! Reallocate storage and fill out arrays
      if (np.gt.0) then
         call smesh%set_size(nvert=nv,npoly=np)
         allocate(smesh%polyConn(smesh%nVert)) ! Also allocate naive connectivity
         nv=0; np=0
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            shape=getNumberOfVertices(this%interface_polygon(i,j,k))
            if (shape.gt.0) then
               ! Increment polygon counter
               np=np+1
               smesh%polySize(np)=shape
               ! Loop over its vertices and add them
               do n=1,shape
                  tmp_vert=getPt(this%interface_polygon(i,j,k),n-1)
                  ! Increment node counter
                  nv=nv+1
                  smesh%xVert(nv)=tmp_vert(1)
                  smesh%yVert(nv)=tmp_vert(2)
                  smesh%zVert(nv)=tmp_vert(3)
                  smesh%polyConn(nv)=nv
               end do
            end if
         end do; end do; end do
      else
         ! Add a zero-area triangle if this proc doesn't have one
         np=1; nv=3
         call smesh%set_size(nvert=nv,npoly=np)
         allocate(smesh%polyConn(smesh%nVert)) ! Also allocate naive connectivity
         smesh%xVert(1:3)=this%cfg%x(this%cfg%imin)
         smesh%yVert(1:3)=this%cfg%y(this%cfg%jmin)
         smesh%zVert(1:3)=this%cfg%z(this%cfg%kmin)
         smesh%polySize(1)=3
         smesh%polyConn(1:3)=[1,2,3]
      end if
   end subroutine update_surfmesh
   
   
   !> Print out info for mpcomp flow solver
   subroutine mpcomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(mpcomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("mpcomp solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine mpcomp_print
   
   
end module mpcomp_class
