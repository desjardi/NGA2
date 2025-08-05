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
      
      ! Conserved variables: 1=VF*RHOL, 2=(1-VF)*RHOG, 3=VF*RHOL*IL, 4=(1-VF)*RHOG*IL, 5=RHO*U, 6=RHO*V, 7=RHO*W
      integer :: nQ
      real(WP), dimension(:,:,:,:), allocatable :: Q,Qold
      
      ! Flow velocity
      real(WP), dimension(:,:,:), allocatable :: U,V,W
      
      ! Phasic densities
      real(WP), dimension(:,:,:), allocatable :: RHOL,RHOG,RHOLold,RHOGold
      
      ! Phasic internal energies
      real(WP), dimension(:,:,:), allocatable :: IL,IG,ILold,IGold
      
      ! Phasic pressures
      real(WP), dimension(:,:,:), allocatable :: PL,PG,PLold,PGold
      
      ! Phasic temperatures
      real(WP), dimension(:,:,:), allocatable :: TL,TG
      
      ! Mixture speed of sound
      real(WP), dimension(:,:,:), allocatable :: C
      
      ! Mixture viscosities and heat diffusivities
      real(WP), dimension(:,:,:), allocatable :: VISC,BETA,DIFF
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP
      
      ! Store mesh info
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,vol
      
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
      
      ! Semi-Lagrangian conservative variable increments
      real(WP), dimension(:,:,:,:), allocatable :: SLdQ
      
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
      procedure :: sync_volume_moments                    !< Helper subroutine that performs the sync for VF, BL, and BG
      procedure :: rhs                                    !< Compute rhs of our equations using standard fluxes
      procedure :: build_interface                        !< Build interface from volume moments
      procedure :: get_primitive                          !< Calculate phasic and mixture primitive variables from conserved variables
      procedure :: apply_relax                            !< Apply user-provided relaxation model in interfacial cells
      procedure :: get_viscartif                          !< Calculate artifical bulk kinematic viscosity
      procedure :: get_vreman                             !< Get kinematic eddy viscosity using Vreman's model
      procedure :: get_velocity                           !< Calculate velocity from momentum
      procedure :: get_ke                                 !< Calculate kinetic energy per unit mass from velocity
      procedure :: get_momentum                           !< Calculate momentum from velocity
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_info                               !< Calculate maximum field values
      procedure :: update_surfmesh                        !< Update a surfmesh object using current polygons
      procedure :: finalize                               !< Finalize mpcomp solver object
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
      !> Mixture relaxation (acts only on the conserved quantities)
      subroutine relax_type(VF,Q)
         import :: WP
         implicit none
         real(WP),                intent(inout) :: VF
         real(WP), dimension(1:), intent(inout) :: Q
      end subroutine relax_type
   end interface
   
contains
   
   
   !> Initialization for compressible flow solver
   subroutine initialize(this,cfg,name)
      use messager, only: die
      implicit none
      class(mpcomp) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to config object
      this%cfg=>cfg
      
      ! Check that config is uniform with at least 3 cells of overlap
      if (this%cfg%no.lt.3) call die('[mpcomp initialize] mpcomp solver requires at least 3 cells of overlap')
      if (.not.all([this%cfg%uniform_x,this%cfg%uniform_y,this%cfg%uniform_z])) call die('[mpcomp initialize] mpcomp solver requires a uniform mesh')
      
      ! Store constant cell size and its inverse, handle 2D conditions, store cell volume
      this%dx=this%cfg%dx(this%cfg%imin_); this%dxi=1.0_WP/this%dx; if (this%cfg%nx.eq.1) this%dxi=0.0_WP
      this%dy=this%cfg%dy(this%cfg%jmin_); this%dyi=1.0_WP/this%dy; if (this%cfg%ny.eq.1) this%dyi=0.0_WP
      this%dz=this%cfg%dz(this%cfg%kmin_); this%dzi=1.0_WP/this%dz; if (this%cfg%nz.eq.1) this%dzi=0.0_WP
      this%vol=this%dx*this%dy*this%dz
      
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
      
      ! Allocate semi-Lagrangian increments
      allocate(this%SLdQ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); this%SLdQ=0.0_WP
      
      ! Initialize Interface Reconstruction Library and its data
      call this%initialize_irl()
      
      ! Allocate semi-Lagrangian tag
      allocate(this%iSL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%iSL=0
      
      ! Conserved variables monitoring
      allocate(this%Qmin(1:this%nQ),this%Qmax(1:this%nQ),this%Qint(1:this%nQ))
      
      ! Flow velocity
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
      
      ! Phasic pressures
      allocate(this%PL   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%PL   =0.0_WP
      allocate(this%PG   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%PG   =0.0_WP
      allocate(this%PLold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%PLold=0.0_WP
      allocate(this%PGold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%PGold=0.0_WP
      
      ! Phasic temperatures
      allocate(this%TL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%TL=0.0_WP
      allocate(this%TG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%TG=0.0_WP
      
      ! Mixture speed of sound
      allocate(this%C(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%C=0.0_WP
      
      ! Mixture viscosities and heat diffusivity
      allocate(this%VISC(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VISC=0.0_WP
      allocate(this%BETA(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BETA=0.0_WP
      allocate(this%DIFF(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%DIFF=0.0_WP
      
      ! Create timers
      this%tplic=timer(comm=this%cfg%comm,name='PLIC')
      this%tsl  =timer(comm=this%cfg%comm,name='SemiLag')
      this%trhs =timer(comm=this%cfg%comm,name='RHS')
      
   end subroutine initialize
   
   
   !> Initialize IRL interface
   subroutine initialize_irl(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k,tag
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
            cell: do i=this%cfg%imino_+1,this%cfg%imaxo_-1
               ! Flag all obvious mixture cells
               if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) then; this%iSL(i,j,k)=1; cycle cell; end if
               ! We may have missed implicit interfaces, check those
               do dir=1,3; do n=-1,+1,2
                  ind=[i,j,k]; ind(dir)=ind(dir)+n; if (this%VF(i,j,k).lt.VFlo.and.this%VF(ind(1),ind(2),ind(3)).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(ind(1),ind(2),ind(3)).lt.VFlo) then; this%iSL(i,j,k)=1; cycle cell; end if
               end do; end do
            end do cell
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
   !> Uses VFold, BLold, BGold, PLICold, RHOLold, RHOGold, ILold, IGold, PLold, PGold
   !> Vertex transport is done with RK2 using passed (U,V,W)
   subroutine SLstep(this,dt,U,V,W)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W
      real(WP), dimension(3,9) :: face
      type(SepVM_type) :: my_SepVM
      type(CapDod_type) :: flux_polyhedron
      type(TagAccVM_SepVM_type) :: detailed_face_flux
      integer :: i,j,k,n
      integer , dimension(3) :: ind
      real(WP), dimension(3) :: Lbar,Gbar
      real(WP) :: Lvol,Gvol,Lmass,Gmass,VFold,div,flux,Lflux,Gflux
      real(WP), dimension(:,:,:,:), allocatable :: SLVx,SLVy,SLVz
      real(WP), dimension(:,:,:,:), allocatable :: SLQx,SLQy,SLQz
      real(WP), dimension(:,:,:,:), allocatable :: SLPx,SLPy,SLPz
      real(WP), dimension(:,:,:)  , allocatable :: dMX,dMY,dMZ
      real(WP), parameter :: Chybrid=-1.5_WP
      real(WP), parameter :: eps=1.0e-12_WP
      
      ! Start semi-Lagrangian timer
      call this%tsl%start()
      
      ! Reset SL increment to zero
      this%SLdQ=0.0_WP
      
      ! Allocate flux polyhedron and detailed face flux
      call new(flux_polyhedron)
      call new(detailed_face_flux)
      
      ! Allocate semi-Lagrangian volume fluxes
      allocate(SLVx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:8)); SLVx=0.0_WP
      allocate(SLVy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:8)); SLVy=0.0_WP
      allocate(SLVz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:8)); SLVz=0.0_WP
      
      ! Allocate semi-Lagrangian conserved variable fluxes
      allocate(SLQx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); SLQx=0.0_WP
      allocate(SLQy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); SLQy=0.0_WP
      allocate(SLQz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); SLQz=0.0_WP
      
      ! Allocate semi-Lagrangian pressure fluxes
      allocate(SLPx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:2)); SLPx=0.0_WP
      allocate(SLPy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:2)); SLPy=0.0_WP
      allocate(SLPz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:2)); SLPz=0.0_WP
      
      ! Loop through all cell faces and get volume moment, mass, and internal energy fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1; do j=this%cfg%jmin_,this%cfg%jmax_+1; do i=this%cfg%imin_,this%cfg%imax_+1
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
            call adjustCapToMatchVolume(flux_polyhedron,U(i,j,k)*dt*this%dy*this%dz)
            ! Build detailed geometric flux
            call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),detailed_face_flux)
            ! Traverse current detailed face flux and increment fluxes
            do n=0,getSize(detailed_face_flux)-1
               ! Get cell index
               ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux,n))
               ! Get separated volume moments
               call getSepVMAtIndex(detailed_face_flux,n,my_SepVM)
               ! Extract volume and mass data
               Lvol=getVolume(my_SepVM,0); Lbar=getCentroid(my_SepVM,0); Lmass=Lvol*this%RHOLold(ind(1),ind(2),ind(3))
               Gvol=getVolume(my_SepVM,1); Gbar=getCentroid(my_SepVM,1); Gmass=Gvol*this%RHOGold(ind(1),ind(2),ind(3))
               ! Increment volume flux
               SLVx(i,j,k,:)=SLVx(i,j,k,:)-[Lvol,Gvol,Lbar,Gbar]
               ! Increment conserved variable flux
               SLQx(i,j,k,:)=SLQx(i,j,k,:)-[Lmass,Gmass,Lmass*this%ILold(ind(1),ind(2),ind(3)),Gmass*this%IGold(ind(1),ind(2),ind(3)),0.0_WP,0.0_WP,0.0_WP]
               ! Increment pressure flux
               SLPx(i,j,k,:)=SLPx(i,j,k,:)-[Lvol*this%PLold(ind(1),ind(2),ind(3)),Gvol*this%PGold(ind(1),ind(2),ind(3))]
            end do
            SLQx(i,j,k,:)=SLQx(i,j,k,:)/(dt*this%dy*this%dz)
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
            call adjustCapToMatchVolume(flux_polyhedron,V(i,j,k)*dt*this%dz*this%dx)
            ! Build detailed geometric flux
            call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),detailed_face_flux)
            ! Traverse current detailed face flux and increment fluxes
            do n=0,getSize(detailed_face_flux)-1
               ! Get cell index
               ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux,n))
               ! Get separated volume moments
               call getSepVMAtIndex(detailed_face_flux,n,my_SepVM)
               ! Extract volume and mass data
               Lvol=getVolume(my_SepVM,0); Lbar=getCentroid(my_SepVM,0); Lmass=Lvol*this%RHOLold(ind(1),ind(2),ind(3))
               Gvol=getVolume(my_SepVM,1); Gbar=getCentroid(my_SepVM,1); Gmass=Gvol*this%RHOGold(ind(1),ind(2),ind(3))
               ! Increment volume flux
               SLVy(i,j,k,:)=SLVy(i,j,k,:)-[Lvol,Gvol,Lbar,Gbar]
               ! Increment conserved variable flux
               SLQy(i,j,k,:)=SLQy(i,j,k,:)-[Lmass,Gmass,Lmass*this%ILold(ind(1),ind(2),ind(3)),Gmass*this%IGold(ind(1),ind(2),ind(3)),0.0_WP,0.0_WP,0.0_WP]
               ! Increment pressure flux
               SLPy(i,j,k,:)=SLPy(i,j,k,:)-[Lvol*this%PLold(ind(1),ind(2),ind(3)),Gvol*this%PGold(ind(1),ind(2),ind(3))]
            end do
            SLQy(i,j,k,:)=SLQy(i,j,k,:)/(dt*this%dz*this%dx)
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
            call adjustCapToMatchVolume(flux_polyhedron,W(i,j,k)*dt*this%dx*this%dy)
            ! Build detailed geometric flux
            call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),detailed_face_flux)
            ! Traverse current detailed face flux and increment fluxes
            do n=0,getSize(detailed_face_flux)-1
               ! Get cell index
               ind=this%cfg%get_ijk_from_lexico(getTagForIndex(detailed_face_flux,n))
               ! Get separated volume moments
               call getSepVMAtIndex(detailed_face_flux,n,my_SepVM)
               ! Extract volume and mass data
               Lvol=getVolume(my_SepVM,0); Lbar=getCentroid(my_SepVM,0); Lmass=Lvol*this%RHOLold(ind(1),ind(2),ind(3))
               Gvol=getVolume(my_SepVM,1); Gbar=getCentroid(my_SepVM,1); Gmass=Gvol*this%RHOGold(ind(1),ind(2),ind(3))
               ! Increment volume flux
               SLVz(i,j,k,:)=SLVz(i,j,k,:)-[Lvol,Gvol,Lbar,Gbar]
               ! Increment conserved variable flux
               SLQz(i,j,k,:)=SLQz(i,j,k,:)-[Lmass,Gmass,Lmass*this%ILold(ind(1),ind(2),ind(3)),Gmass*this%IGold(ind(1),ind(2),ind(3)),0.0_WP,0.0_WP,0.0_WP]
               ! Increment pressure flux
               SLPz(i,j,k,:)=SLPz(i,j,k,:)-[Lvol*this%PLold(ind(1),ind(2),ind(3)),Gvol*this%PGold(ind(1),ind(2),ind(3))]
            end do
            SLQz(i,j,k,:)=SLQz(i,j,k,:)/(dt*this%dx*this%dy)
            ! Clear detailed flux
            call clear(detailed_face_flux)
         end if
      end do; end do; end do
      
      ! Mass fluxes will be used to build momentum fluxes, they need to be extended by one cell on the left because of staggering
      call this%cfg%sync(SLQx(:,:,:,1)); call this%cfg%sync(SLQx(:,:,:,2)); if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) SLQx(this%cfg%imin-1,:,:,1:2)=SLQx(this%cfg%imin,:,:,1:2)
      call this%cfg%sync(SLQy(:,:,:,1)); call this%cfg%sync(SLQy(:,:,:,2)); if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) SLQy(:,this%cfg%jmin-1,:,1:2)=SLQy(:,this%cfg%jmin,:,1:2)
      call this%cfg%sync(SLQz(:,:,:,1)); call this%cfg%sync(SLQz(:,:,:,2)); if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) SLQz(:,:,this%cfg%kmin-1,1:2)=SLQz(:,:,this%cfg%kmin,1:2)
      
      ! Update volume moments in tagged cells, also fix potential missing pressures
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Only work at interface
         if (this%iSL(i,j,k).eq.0) cycle
         ! Compute new liquid and gas volumes
         Lvol=(       this%VFold(i,j,k))*this%vol+SLVx(i+1,j,k,1)-SLVx(i,j,k,1)+SLVy(i,j+1,k,1)-SLVy(i,j,k,1)+SLVz(i,j,k+1,1)-SLVz(i,j,k,1)
         Gvol=(1.0_WP-this%VFold(i,j,k))*this%vol+SLVx(i+1,j,k,2)-SLVx(i,j,k,2)+SLVy(i,j+1,k,2)-SLVy(i,j,k,2)+SLVz(i,j,k+1,2)-SLVz(i,j,k,2)
         ! Compute new liquid volume fraction while remembering the current
         VFold=this%VF(i,j,k)
         this%VF(i,j,k)=Lvol/(Lvol+Gvol)
         ! Clip it within [VFlo,VFhi] and compute first order moments
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
            this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
            this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
         else
            ! Project old barycenters forward in time
            this%BL(:,i,j,k)=project((this%BLold(:,i,j,k)*(       this%VFold(i,j,k))*this%vol+SLVx(i+1,j,k,3:5)-SLVx(i,j,k,3:5)+SLVy(i,j+1,k,3:5)-SLVy(i,j,k,3:5)+SLVz(i,j,k+1,3:5)-SLVz(i,j,k,3:5))/Lvol,dt)
            this%BG(:,i,j,k)=project((this%BGold(:,i,j,k)*(1.0_WP-this%VFold(i,j,k))*this%vol+SLVx(i+1,j,k,6:8)-SLVx(i,j,k,6:8)+SLVy(i,j+1,k,6:8)-SLVy(i,j,k,6:8)+SLVz(i,j,k+1,6:8)-SLVz(i,j,k,6:8))/Gvol,dt)
         end if
         ! Compute new liquid and gas pressures in newly created cells
         if (VFold.lt.VFlo.and.this%VF(i,j,k).ge.VFlo) this%PL(i,j,k)=(SLPx(i+1,j,k,1)-SLPx(i,j,k,1)+SLPy(i,j+1,k,1)-SLPy(i,j,k,1)+SLPz(i,j,k+1,1)-SLPz(i,j,k,1))/Lvol
         if (VFold.gt.VFhi.and.this%VF(i,j,k).le.VFhi) this%PG(i,j,k)=(SLPx(i+1,j,k,2)-SLPx(i,j,k,2)+SLPy(i,j+1,k,2)-SLPy(i,j,k,2)+SLPz(i,j,k+1,2)-SLPz(i,j,k,2))/Gvol
      end do; end do; end do
      
      ! Deallocate volume, pressure flux arrays
      deallocate(SLVx,SLVy,SLVz,SLPx,SLPy,SLPz)
      
      ! Synchronize pressures
      call this%cfg%sync(this%PL); call this%cfg%sync(this%PG)
      
      ! Synchronize VF and barycenters, fix barycenter synchronization across periodic boundaries, and handle 2D barycenters
      call this%sync_volume_moments()
      
      ! ================================================================ !
      ! ======================== INVISID FLUXES ======================== !
      ! ================================================================ !
      
      ! Form the SL increment for phasic mass and energy
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Phasic mass and internal energy advection
               this%SLdQ(i,j,k,1:4)=dt*(this%dxi*(SLQx(i+1,j,k,1:4)-SLQx(i,j,k,1:4))+this%dyi*(SLQy(i,j+1,k,1:4)-SLQy(i,j,k,1:4))+this%dzi*(SLQz(i,j,k+1,1:4)-SLQz(i,j,k,1:4)))
               ! Pressure dilatation term
               if (this%iSL(i,j,k).gt.0) then
                  div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
                  this%SLdQ(i,j,k,3)=this%SLdQ(i,j,k,3)-dt*(       this%VF(i,j,k))*this%PL(i,j,k)*div
                  this%SLdQ(i,j,k,4)=this%SLdQ(i,j,k,4)-dt*(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)*div
               end if
            end do
         end do
      end do
      
      ! Synchronize phasic dQ fields
      do n=1,4; call this%cfg%sync(this%SLdQ(:,:,:,n)); end do
      
      ! Calculate normalized mass change on staggered cells
      allocate(dMX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); dMX=0.0_WP
      allocate(dMY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); dMY=0.0_WP
      allocate(dMZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); dMZ=0.0_WP
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_; do j=this%cfg%jmino_+1,this%cfg%jmaxo_; do i=this%cfg%imino_+1,this%cfg%imaxo_
         dMX(i,j,k)=sum(this%SLdQ(i-1:i,j,k,1:2))/sum(this%Qold(i-1:i,j,k,1:2)+this%SLdQ(i-1:i,j,k,1:2))
         dMY(i,j,k)=sum(this%SLdQ(i,j-1:j,k,1:2))/sum(this%Qold(i,j-1:j,k,1:2)+this%SLdQ(i,j-1:j,k,1:2))
         dMZ(i,j,k)=sum(this%SLdQ(i,j,k-1:k,1:2))/sum(this%Qold(i,j,k-1:k,1:2)+this%SLdQ(i,j,k-1:k,1:2))
      end do; end do; end do
      call this%cfg%sync(dMX); call this%cfg%sync(dMY); call this%cfg%sync(dMZ)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then; dMX(this%cfg%imino,:,:)=dMX(this%cfg%imino+1,:,:); dMY(this%cfg%imino,:,:)=dMY(this%cfg%imino+1,:,:); dMZ(this%cfg%imino,:,:)=dMZ(this%cfg%imino+1,:,:); end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then; dMX(:,this%cfg%jmino,:)=dMX(:,this%cfg%jmino+1,:); dMY(:,this%cfg%jmino,:)=dMY(:,this%cfg%jmino+1,:); dMZ(:,this%cfg%jmino,:)=dMZ(:,this%cfg%jmino+1,:); end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then; dMX(:,:,this%cfg%kmino)=dMX(:,:,this%cfg%kmino+1); dMY(:,:,this%cfg%kmino)=dMY(:,:,this%cfg%kmino+1); dMZ(:,:,this%cfg%kmino)=dMZ(:,:,this%cfg%kmino+1); end if
      
      ! Calculate mixture momentum fluxes from SL mass fluxes with extra cell on the left due to staggering
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               Lflux=0.5_WP*sum(SLQx(i:i+1,j,k,1)); Gflux=0.5_WP*sum(SLQx(i:i+1,j,k,2)); flux=Lflux+Gflux; SLQx(i,j,k,5)=flux*0.5_WP*sum(this%U(i:i+1,j,k)); if (minval(dMX(i:i+1,j,k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQx(i,j,k,5)=0.5_WP*(flux-abs(-flux))*this%U(i,j,k)+0.5_WP*(flux+abs(-flux))*this%U(i+1,j,k)
               Lflux=0.5_WP*sum(SLQy(i,j:j+1,k,1)); Gflux=0.5_WP*sum(SLQy(i,j:j+1,k,2)); flux=Lflux+Gflux; SLQy(i,j,k,6)=flux*0.5_WP*sum(this%V(i,j:j+1,k)); if (minval(dMY(i,j:j+1,k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQy(i,j,k,6)=0.5_WP*(flux-abs(-flux))*this%V(i,j,k)+0.5_WP*(flux+abs(-flux))*this%V(i,j+1,k)
               Lflux=0.5_WP*sum(SLQz(i,j,k:k+1,1)); Gflux=0.5_WP*sum(SLQz(i,j,k:k+1,2)); flux=Lflux+Gflux; SLQz(i,j,k,7)=flux*0.5_WP*sum(this%W(i,j,k:k+1)); if (minval(dMZ(i,j,k:k+1)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQz(i,j,k,7)=0.5_WP*(flux-abs(-flux))*this%W(i,j,k)+0.5_WP*(flux+abs(-flux))*this%W(i,j,k+1)
               ! Add pressure fluxes in SL cells
               if (this%iSL(i,j,k).gt.0) then
                  SLQx(i,j,k,5)=SLQx(i,j,k,5)-this%VF(i,j,k)*this%PL(i,j,k)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)
                  SLQy(i,j,k,6)=SLQy(i,j,k,6)-this%VF(i,j,k)*this%PL(i,j,k)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)
                  SLQz(i,j,k,7)=SLQz(i,j,k,7)-this%VF(i,j,k)*this%PL(i,j,k)-(1.0_WP-this%VF(i,j,k))*this%PG(i,j,k)
               end if
            end do
         end do
      end do
      
      ! Calculate edge-centered mixture momentum fluxes from SL mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               Lflux=0.5_WP*sum(SLQy(i-1:i,j,k,1)); Gflux=0.5_WP*sum(SLQy(i-1:i,j,k,2)); flux=Lflux+Gflux; SLQy(i,j,k,5)=flux*0.5_WP*sum(this%U(i,j-1:j,k)); if (minval(dMX(i,j-1:j,k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQy(i,j,k,5)=0.5_WP*(flux-abs(-flux))*this%U(i,j-1,k)+0.5_WP*(flux+abs(-flux))*this%U(i,j,k)
               Lflux=0.5_WP*sum(SLQz(i-1:i,j,k,1)); Gflux=0.5_WP*sum(SLQz(i-1:i,j,k,2)); flux=Lflux+Gflux; SLQz(i,j,k,5)=flux*0.5_WP*sum(this%U(i,j,k-1:k)); if (minval(dMX(i,j,k-1:k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQz(i,j,k,5)=0.5_WP*(flux-abs(-flux))*this%U(i,j,k-1)+0.5_WP*(flux+abs(-flux))*this%U(i,j,k)
               Lflux=0.5_WP*sum(SLQx(i,j-1:j,k,1)); Gflux=0.5_WP*sum(SLQx(i,j-1:j,k,2)); flux=Lflux+Gflux; SLQx(i,j,k,6)=flux*0.5_WP*sum(this%V(i-1:i,j,k)); if (minval(dMY(i-1:i,j,k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQx(i,j,k,6)=0.5_WP*(flux-abs(-flux))*this%V(i-1,j,k)+0.5_WP*(flux+abs(-flux))*this%V(i,j,k)
               Lflux=0.5_WP*sum(SLQz(i,j-1:j,k,1)); Gflux=0.5_WP*sum(SLQz(i,j-1:j,k,2)); flux=Lflux+Gflux; SLQz(i,j,k,6)=flux*0.5_WP*sum(this%V(i,j,k-1:k)); if (minval(dMY(i,j,k-1:k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQz(i,j,k,6)=0.5_WP*(flux-abs(-flux))*this%V(i,j,k-1)+0.5_WP*(flux+abs(-flux))*this%V(i,j,k)
               Lflux=0.5_WP*sum(SLQx(i,j,k-1:k,1)); Gflux=0.5_WP*sum(SLQx(i,j,k-1:k,2)); flux=Lflux+Gflux; SLQx(i,j,k,7)=flux*0.5_WP*sum(this%W(i-1:i,j,k)); if (minval(dMZ(i-1:i,j,k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQx(i,j,k,7)=0.5_WP*(flux-abs(-flux))*this%W(i-1,j,k)+0.5_WP*(flux+abs(-flux))*this%W(i,j,k)
               Lflux=0.5_WP*sum(SLQy(i,j,k-1:k,1)); Gflux=0.5_WP*sum(SLQy(i,j,k-1:k,2)); flux=Lflux+Gflux; SLQy(i,j,k,7)=flux*0.5_WP*sum(this%W(i,j-1:j,k)); if (minval(dMZ(i,j-1:j,k)).lt.Chybrid.and.min(abs(Lflux),abs(Gflux)).gt.eps*abs(flux)) SLQy(i,j,k,7)=0.5_WP*(flux-abs(-flux))*this%W(i,j-1,k)+0.5_WP*(flux+abs(-flux))*this%W(i,j,k)
            end do
         end do
      end do
      
      ! Form SL increment for mixture momentum
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%SLdQ(i,j,k,5)=dt*(this%dxi*(SLQx(i  ,j,k,5)-SLQx(i-1,j,k,5))+this%dyi*(SLQy(i,j+1,k,5)-SLQy(i,j  ,k,5))+this%dzi*(SLQz(i,j,k+1,5)-SLQz(i,j,k  ,5)))
               this%SLdQ(i,j,k,6)=dt*(this%dxi*(SLQx(i+1,j,k,6)-SLQx(i  ,j,k,6))+this%dyi*(SLQy(i,j  ,k,6)-SLQy(i,j-1,k,6))+this%dzi*(SLQz(i,j,k+1,6)-SLQz(i,j,k  ,6)))
               this%SLdQ(i,j,k,7)=dt*(this%dxi*(SLQx(i+1,j,k,7)-SLQx(i  ,j,k,7))+this%dyi*(SLQy(i,j+1,k,7)-SLQy(i,j  ,k,7))+this%dzi*(SLQz(i,j,k  ,7)-SLQz(i,j,k-1,7)))
            end do
         end do
      end do
      
      ! Deallocate staggered mass change
      deallocate(dMX,dMY,dMZ)
      
      ! ================================================================ !
      ! ======================== VISCOUS  FLUXES ======================= !
      ! ================================================================ !
      
      ! Zero out all fluxes
      SLQx=0.0_WP; SLQy=0.0_WP; SLQz=0.0_WP
      
      ! Compute cell-centered momentum viscous fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               if (this%iSL(i,j,k).gt.0) then
                  div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
                  SLQx(i,j,k,5)=(2.0_WP*this%VISC(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div)
                  SLQy(i,j,k,6)=(2.0_WP*this%VISC(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div)
                  SLQz(i,j,k,7)=(2.0_WP*this%VISC(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div)
               end if
            end do
         end do
      end do
      
      ! Compute edge-centered momentum viscous fluxes and corresponding viscous heating
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (maxval(this%iSL(i-1:i,j-1:j,k)).gt.0) then
                  SLQy(i,j,k,5)=0.25_WP*sum(this%VISC(i-1:i,j-1:j,k))*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k))); SLQx(i,j,k,6)=SLQy(i,j,k,5)
                  SLQz(i,j,k,3:4)=SLQy(i,j,k,5)*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
               end if
               if (maxval(this%iSL(i,j-1:j,k-1:k)).gt.0) then
                  SLQz(i,j,k,6)=0.25_WP*sum(this%VISC(i,j-1:j,k-1:k))*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k))); SLQy(i,j,k,7)=SLQz(i,j,k,6)
                  SLQx(i,j,k,3:4)=SLQz(i,j,k,6)*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
               end if
               if (maxval(this%iSL(i-1:i,j,k-1:k)).gt.0) then
                  SLQx(i,j,k,7)=0.25_WP*sum(this%VISC(i-1:i,j,k-1:k))*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1))); SLQz(i,j,k,5)=SLQx(i,j,k,7)
                  SLQy(i,j,k,3:4)=SLQx(i,j,k,7)*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
               end if
            end do
         end do
      end do
      
      ! Update the SL increment for all conserved variables
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Viscous momentum transport
               this%SLdQ(i,j,k,5)=this%SLdQ(i,j,k,5)+dt*(this%dxi*(SLQx(i  ,j,k,5)-SLQx(i-1,j,k,5))+this%dyi*(SLQy(i,j+1,k,5)-SLQy(i,j  ,k,5))+this%dzi*(SLQz(i,j,k+1,5)-SLQz(i,j,k  ,5)))
               this%SLdQ(i,j,k,6)=this%SLdQ(i,j,k,6)+dt*(this%dxi*(SLQx(i+1,j,k,6)-SLQx(i  ,j,k,6))+this%dyi*(SLQy(i,j  ,k,6)-SLQy(i,j-1,k,6))+this%dzi*(SLQz(i,j,k+1,6)-SLQz(i,j,k  ,6)))
               this%SLdQ(i,j,k,7)=this%SLdQ(i,j,k,7)+dt*(this%dxi*(SLQx(i+1,j,k,7)-SLQx(i  ,j,k,7))+this%dyi*(SLQy(i,j+1,k,7)-SLQy(i,j  ,k,7))+this%dzi*(SLQz(i,j,k  ,7)-SLQz(i,j,k-1,7)))
               ! Distribute viscous heating term based on VF
               this%SLdQ(i,j,k,3)=this%SLdQ(i,j,k,3)+dt*(       this%VF(i,j,k))*((SLQx(i,j,k,5)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+SLQy(i,j,k,6)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+SLQz(i,j,k,7)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+0.25_WP*sum(SLQz(i:i+1,j:j+1,k,3))+0.25_WP*sum(SLQx(i,j:j+1,k:k+1,3))+0.25_WP*sum(SLQy(i:i+1,j,k:k+1,3))))
               this%SLdQ(i,j,k,4)=this%SLdQ(i,j,k,4)+dt*(1.0_WP-this%VF(i,j,k))*((SLQx(i,j,k,5)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+SLQy(i,j,k,6)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+SLQz(i,j,k,7)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+0.25_WP*sum(SLQz(i:i+1,j:j+1,k,4))+0.25_WP*sum(SLQx(i,j:j+1,k:k+1,4))+0.25_WP*sum(SLQy(i:i+1,j,k:k+1,4))))
            end do
         end do
      end do
      
      ! Deallocate all flux arrays
      deallocate(SLQx,SLQy,SLQz)
      
      ! Synchronize all Q fields
      do n=1,this%nQ; call this%cfg%sync(this%SLdQ(:,:,:,n)); end do
      
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
         integer  :: ipc,jpc,kpc,ipu,jpv,kpw
         real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2
         real(WP) :: wxu1,wyv1,wzw1,wxu2,wyv2,wzw2
         ! Cell-centered indices and weights
         ipc=max(this%cfg%imino_,min(this%cfg%imaxo_-1,floor((pos(1)-this%cfg%xm(this%cfg%imino))*this%dxi)+this%cfg%imino)); wxc1=max(0.0_WP,min(1.0_WP,(pos(1)-this%cfg%xm(ipc))*this%dxi)); wxc2=1.0_WP-wxc1
         jpc=max(this%cfg%jmino_,min(this%cfg%jmaxo_-1,floor((pos(2)-this%cfg%ym(this%cfg%jmino))*this%dyi)+this%cfg%jmino)); wyc1=max(0.0_WP,min(1.0_WP,(pos(2)-this%cfg%ym(jpc))*this%dyi)); wyc2=1.0_WP-wyc1
         kpc=max(this%cfg%kmino_,min(this%cfg%kmaxo_-1,floor((pos(3)-this%cfg%zm(this%cfg%kmino))*this%dzi)+this%cfg%kmino)); wzc1=max(0.0_WP,min(1.0_WP,(pos(3)-this%cfg%zm(kpc))*this%dzi)); wzc2=1.0_WP-wzc1
         ! Face-centered indices and weights
         ipu=max(this%cfg%imino_,min(this%cfg%imaxo_-1,floor((pos(1)-this%cfg%x (this%cfg%imino))*this%dxi)+this%cfg%imino)); wxu1=max(0.0_WP,min(1.0_WP,(pos(1)-this%cfg%x (ipu))*this%dxi)); wxu2=1.0_WP-wxu1
         jpv=max(this%cfg%jmino_,min(this%cfg%jmaxo_-1,floor((pos(2)-this%cfg%y (this%cfg%jmino))*this%dyi)+this%cfg%jmino)); wyv1=max(0.0_WP,min(1.0_WP,(pos(2)-this%cfg%y (jpv))*this%dyi)); wyv2=1.0_WP-wyv1
         kpw=max(this%cfg%kmino_,min(this%cfg%kmaxo_-1,floor((pos(3)-this%cfg%z (this%cfg%kmino))*this%dzi)+this%cfg%kmino)); wzw1=max(0.0_WP,min(1.0_WP,(pos(3)-this%cfg%z (kpw))*this%dzi)); wzw2=1.0_WP-wzw1
         ! Tri-linear interpolation of staggered velocity
         vel(1)=wzc1*(wyc1*(wxu1*U(ipu+1,jpc+1,kpc+1)+wxu2*U(ipu,jpc+1,kpc+1))+wyc2*(wxu1*U(ipu+1,jpc,kpc+1)+wxu2*U(ipu,jpc,kpc+1)))+wzc2*(wyc1*(wxu1*U(ipu+1,jpc+1,kpc  )+wxu2*U(ipu,jpc+1,kpc  ))+wyc2*(wxu1*U(ipu+1,jpc,kpc  )+wxu2*U(ipu,jpc,kpc  )))
         vel(2)=wzc1*(wyv1*(wxc1*V(ipc+1,jpv+1,kpc+1)+wxc2*V(ipc,jpv+1,kpc+1))+wyv2*(wxc1*V(ipc+1,jpv,kpc+1)+wxc2*V(ipc,jpv,kpc+1)))+wzc2*(wyv1*(wxc1*V(ipc+1,jpv+1,kpc  )+wxc2*V(ipc,jpv+1,kpc  ))+wyv2*(wxc1*V(ipc+1,jpv,kpc  )+wxc2*V(ipc,jpv,kpc  )))
         vel(3)=wzw1*(wyc1*(wxc1*W(ipc+1,jpc+1,kpw+1)+wxc2*W(ipc,jpc+1,kpw+1))+wyc2*(wxc1*W(ipc+1,jpc,kpw+1)+wxc2*W(ipc,jpc,kpw+1)))+wzw2*(wyc1*(wxc1*W(ipc+1,jpc+1,kpw  )+wxc2*W(ipc,jpc+1,kpw  ))+wyc2*(wxc1*W(ipc+1,jpc,kpw  )+wxc2*W(ipc,jpc,kpw  )))
      end function interp_velocity
   end subroutine SLstep
   
   
   !> Sync volume moments (VF, BL, and BG)
   subroutine sync_volume_moments(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      integer :: i,j,k
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
   end subroutine sync_volume_moments
   
   
   !> Obtain non-SL RHS for all equations except volume
   subroutine rhs(this,dQdt)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: dQdt  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nVAR)
      real(WP), dimension(:,:,:,:), allocatable :: FQx,FQy,FQz
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
      
      ! Allocate fluxes of conserved variables
      allocate(FQx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); FQx=0.0_WP
      allocate(FQy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); FQy=0.0_WP
      allocate(FQz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); FQz=0.0_WP
      
      ! Calculate standard (i.e., non-SL) fluxes in the appropriate phase
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (maxval(this%iSL(i-1:i,j,k)).eq.0) then
                  if (this%VF(i,j,k).gt.0.5_WP) then
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOL(i-1,j,k)-this%RHOL(i-2,j,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOL(i+1,j,k)-this%RHOL(i  ,j,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQx(i,j,k,1)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wenop*this%RHOL(i-2:i  ,j,k))&
                     &            -0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wenom*this%RHOL(i-1:i+1,j,k))
                     ! Centered phasic mass flux
                     !FQx(i,j,k,1)=-this%U(i,j,k)*0.5_WP*sum(this%RHOL(i-1:i,j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IL(i-1,j,k)-this%IL(i-2,j,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IL(i+1,j,k)-this%IL(i  ,j,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQx(i,j,k,3)=0.5_WP*(FQx(i,j,k,1)-abs(-FQx(i,j,k,1)))*sum(wenop*this%IL(i-2:i  ,j,k))&
                     &           +0.5_WP*(FQx(i,j,k,1)+abs(-FQx(i,j,k,1)))*sum(wenom*this%IL(i-1:i+1,j,k))
                     ! Centered phasic internal energy flux
                     !FQx(i,j,k,3)=FQx(i,j,k,1)*0.5_WP*sum(this%IL(i-1:i,j,k))
                  else
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOG(i-1,j,k)-this%RHOG(i-2,j,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOG(i+1,j,k)-this%RHOG(i  ,j,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQx(i,j,k,2)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wenop*this%RHOG(i-2:i  ,j,k))&
                     &            -0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wenom*this%RHOG(i-1:i+1,j,k))
                     ! Centered phasic mass flux
                     !FQx(i,j,k,2)=-this%U(i,j,k)*0.5_WP*sum(this%RHOG(i-1:i,j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IG(i-1,j,k)-this%IG(i-2,j,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IG(i+1,j,k)-this%IG(i  ,j,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQx(i,j,k,4)=0.5_WP*(FQx(i,j,k,2)-abs(-FQx(i,j,k,2)))*sum(wenop*this%IG(i-2:i  ,j,k))&
                     &           +0.5_WP*(FQx(i,j,k,2)+abs(-FQx(i,j,k,2)))*sum(wenom*this%IG(i-1:i+1,j,k))
                     ! Centered phasic internal energy flux
                     !FQx(i,j,k,4)=FQx(i,j,k,2)*0.5_WP*sum(this%IG(i-1:i,j,k))
                  end if
               end if
               if (maxval(this%iSL(i,j-1:j,k)).eq.0) then
                  if (this%VF(i,j,k).gt.0.5_WP) then
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOL(i,j-1,k)-this%RHOL(i,j-2,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOL(i,j+1,k)-this%RHOL(i,j  ,k))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQy(i,j,k,1)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wenop*this%RHOL(i,j-2:j  ,k))&
                     &            -0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wenom*this%RHOL(i,j-1:j+1,k))
                     ! Centered phasic mass flux
                     !FQy(i,j,k,1)=-this%V(i,j,k)*0.5_WP*sum(this%RHOL(i,j-1:j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IL(i,j-1,k)-this%IL(i,j-2,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IL(i,j+1,k)-this%IL(i,j  ,k))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQy(i,j,k,3)=0.5_WP*(FQy(i,j,k,1)-abs(-FQy(i,j,k,1)))*sum(wenop*this%IL(i,j-2:j  ,k))&
                     &           +0.5_WP*(FQy(i,j,k,1)+abs(-FQy(i,j,k,1)))*sum(wenom*this%IL(i,j-1:j+1,k))
                     ! Centered phasic internal energy flux
                     !FQy(i,j,k,3)=FQy(i,j,k,1)*0.5_WP*sum(this%IL(i,j-1:j,k))
                  else
                     ! WENO phasic mass flux
                     w=weno_weight((abs(this%RHOG(i,j-1,k)-this%RHOG(i,j-2,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOG(i,j+1,k)-this%RHOG(i,j  ,k))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQy(i,j,k,2)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wenop*this%RHOG(i,j-2:j  ,k))&
                     &            -0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wenom*this%RHOG(i,j-1:j+1,k))
                     ! Centered phasic mass flux
                     !FQy(i,j,k,2)=-this%V(i,j,k)*0.5_WP*sum(this%RHOG(i,j-1:j,k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IG(i,j-1,k)-this%IG(i,j-2,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IG(i,j+1,k)-this%IG(i,j  ,k))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQy(i,j,k,4)=0.5_WP*(FQy(i,j,k,2)-abs(-FQy(i,j,k,2)))*sum(wenop*this%IG(i,j-2:j  ,k))&
                     &           +0.5_WP*(FQy(i,j,k,2)+abs(-FQy(i,j,k,2)))*sum(wenom*this%IG(i,j-1:j+1,k))
                     ! Centered phasic internal energy flux
                     !FQy(i,j,k,4)=FQy(i,j,k,2)*0.5_WP*sum(this%IG(i,j-1:j,k))
                  end if
               end if
               if (maxval(this%iSL(i,j,k-1:k)).eq.0) then
                  if (this%VF(i,j,k).gt.0.5_WP) then
                     ! WENO mass flux
                     w=weno_weight((abs(this%RHOL(i,j,k-1)-this%RHOL(i,j,k-2))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOL(i,j,k+1)-this%RHOL(i,j,k  ))+eps)/(abs(this%RHOL(i,j,k)-this%RHOL(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQz(i,j,k,1)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wenop*this%RHOL(i,j,k-2:k  ))&
                     &            -0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wenom*this%RHOL(i,j,k-1:k+1))
                     ! Centered mass flux
                     !FQz(i,j,k,1)=-this%W(i,j,k)*0.5_WP*sum(this%RHOL(i,j,k-1:k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IL(i,j,k-1)-this%IL(i,j,k-2))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IL(i,j,k+1)-this%IL(i,j,k  ))+eps)/(abs(this%IL(i,j,k)-this%IL(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQz(i,j,k,3)=0.5_WP*(FQz(i,j,k,1)-abs(-FQz(i,j,k,1)))*sum(wenop*this%IL(i,j,k-2:k  ))&
                     &           +0.5_WP*(FQz(i,j,k,1)+abs(-FQz(i,j,k,1)))*sum(wenom*this%IL(i,j,k-1:k+1))
                     ! Centered internal energy flux
                     !FQz(i,j,k,3)=FQz(i,j,k,1)*0.5_WP*sum(this%IL(i,j,k-1:k))
                  else
                     ! WENO mass flux
                     w=weno_weight((abs(this%RHOG(i,j,k-1)-this%RHOG(i,j,k-2))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%RHOG(i,j,k+1)-this%RHOG(i,j,k  ))+eps)/(abs(this%RHOG(i,j,k)-this%RHOG(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQz(i,j,k,2)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wenop*this%RHOG(i,j,k-2:k  ))&
                     &            -0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wenom*this%RHOG(i,j,k-1:k+1))
                     ! Centered mass flux
                     !FQz(i,j,k,2)=-this%W(i,j,k)*0.5_WP*sum(this%RHOG(i,j,k-1:k))
                     ! WENO phasic internal energy flux
                     w=weno_weight((abs(this%IG(i,j,k-1)-this%IG(i,j,k-2))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(this%IG(i,j,k+1)-this%IG(i,j,k  ))+eps)/(abs(this%IG(i,j,k)-this%IG(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     FQz(i,j,k,4)=0.5_WP*(FQz(i,j,k,2)-abs(-FQz(i,j,k,2)))*sum(wenop*this%IG(i,j,k-2:k  ))&
                     &           +0.5_WP*(FQz(i,j,k,2)+abs(-FQz(i,j,k,2)))*sum(wenom*this%IG(i,j,k-1:k+1))
                     ! Centered internal energy flux
                     !FQz(i,j,k,4)=FQz(i,j,k,2)*0.5_WP*sum(this%IG(i,j,k-1:k))
                  end if
               end if
            end do
         end do
      end do
      
      ! Mass fluxes will be used to build momentum fluxes, they need to be extended by one cell on the left because of staggering
      call this%cfg%sync(FQx(:,:,:,1)); call this%cfg%sync(FQx(:,:,:,2)); if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) FQx(this%cfg%imin-1,:,:,1:2)=FQx(this%cfg%imin,:,:,1:2)
      call this%cfg%sync(FQy(:,:,:,1)); call this%cfg%sync(FQy(:,:,:,2)); if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) FQy(:,this%cfg%jmin-1,:,1:2)=FQy(:,this%cfg%jmin,:,1:2)
      call this%cfg%sync(FQz(:,:,:,1)); call this%cfg%sync(FQz(:,:,:,2)); if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) FQz(:,:,this%cfg%kmin-1,1:2)=FQz(:,:,this%cfg%kmin,1:2)
      
      ! Calculate cell-centered mixture momentum fluxes from standard mass fluxes with extra cell on the left due to staggering
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               FQx(i,j,k,5)=0.25_WP*sum(FQx(i:i+1,j,k,1:2))*sum(this%U(i:i+1,j,k))
               FQy(i,j,k,6)=0.25_WP*sum(FQy(i,j:j+1,k,1:2))*sum(this%V(i,j:j+1,k))
               FQz(i,j,k,7)=0.25_WP*sum(FQz(i,j,k:k+1,1:2))*sum(this%W(i,j,k:k+1))
               ! Add pressure fluxes in non-SL cells
               if (this%iSL(i,j,k).eq.0.and.this%VF(i,j,k).ge.0.5_WP) then
                  FQx(i,j,k,5)=FQx(i,j,k,5)-this%PL(i,j,k)
                  FQy(i,j,k,6)=FQy(i,j,k,6)-this%PL(i,j,k)
                  FQz(i,j,k,7)=FQz(i,j,k,7)-this%PL(i,j,k)
               end if
               if (this%iSL(i,j,k).eq.0.and.this%VF(i,j,k).lt.0.5_WP) then
                  FQx(i,j,k,5)=FQx(i,j,k,5)-this%PG(i,j,k)
                  FQy(i,j,k,6)=FQy(i,j,k,6)-this%PG(i,j,k)
                  FQz(i,j,k,7)=FQz(i,j,k,7)-this%PG(i,j,k)
               end if
            end do
         end do
      end do
      
      ! Calculate edge-centered mixture momentum fluxes from standard mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FQy(i,j,k,5)=0.25_WP*sum(FQy(i-1:i,j,k,1:2))*sum(this%U(i,j-1:j,k))
               FQz(i,j,k,5)=0.25_WP*sum(FQz(i-1:i,j,k,1:2))*sum(this%U(i,j,k-1:k))
               FQx(i,j,k,6)=0.25_WP*sum(FQx(i,j-1:j,k,1:2))*sum(this%V(i-1:i,j,k))
               FQz(i,j,k,6)=0.25_WP*sum(FQz(i,j-1:j,k,1:2))*sum(this%V(i,j,k-1:k))
               FQx(i,j,k,7)=0.25_WP*sum(FQx(i,j,k-1:k,1:2))*sum(this%W(i-1:i,j,k))
               FQy(i,j,k,7)=0.25_WP*sum(FQy(i,j,k-1:k,1:2))*sum(this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Assemble time derivative for conserved variables using terms built from standard mass fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Phasic mass and phasic internal energy advection
               dQdt(i,j,k,1:4)=this%dxi*(FQx(i+1,j,k,1:4)-FQx(i,j,k,1:4))+this%dyi*(FQy(i,j+1,k,1:4)-FQy(i,j,k,1:4))+this%dzi*(FQz(i,j,k+1,1:4)-FQz(i,j,k,1:4))
               ! Mixture momentum advection and pressure stress
               dQdt(i,j,k,5)=this%dxi*(FQx(i  ,j,k,5)-FQx(i-1,j,k,5))+this%dyi*(FQy(i,j+1,k,5)-FQy(i,j  ,k,5))+this%dzi*(FQz(i,j,k+1,5)-FQz(i,j,k  ,5))
               dQdt(i,j,k,6)=this%dxi*(FQx(i+1,j,k,6)-FQx(i  ,j,k,6))+this%dyi*(FQy(i,j  ,k,6)-FQy(i,j-1,k,6))+this%dzi*(FQz(i,j,k+1,6)-FQz(i,j,k  ,6))
               dQdt(i,j,k,7)=this%dxi*(FQx(i+1,j,k,7)-FQx(i  ,j,k,7))+this%dyi*(FQy(i,j+1,k,7)-FQy(i,j  ,k,7))+this%dzi*(FQz(i,j,k  ,7)-FQz(i,j,k-1,7))
               ! Pressure dilatation term
               if (this%iSL(i,j,k).eq.0.and.this%VF(i,j,k).ge.0.5_WP) dQdt(i,j,k,3)=dQdt(i,j,k,3)-this%PL(i,j,k)*(this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k)))
               if (this%iSL(i,j,k).eq.0.and.this%VF(i,j,k).lt.0.5_WP) dQdt(i,j,k,4)=dQdt(i,j,k,4)-this%PG(i,j,k)*(this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k)))
            end do
         end do
      end do
      
      ! ================================================================ !
      ! ======================== VISCOUS  FLUXES ======================= !
      ! ================================================================ !
      
      ! Zero out fluxes
      FQx=0.0_WP; FQy=0.0_WP; FQz=0.0_WP
      
      ! Compute cell-centered momentum viscous fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               if (this%iSL(i,j,k).eq.0) then
                  div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
                  FQx(i,j,k,5)=2.0_WP*this%VISC(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div
                  FQy(i,j,k,6)=2.0_WP*this%VISC(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div
                  FQz(i,j,k,7)=2.0_WP*this%VISC(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div
               end if
            end do
         end do
      end do
      
      ! Compute edge-centered momentum viscous fluxes and corresponding viscous heating
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (maxval(this%iSL(i-1:i,j-1:j,k)).eq.0) then
                  FQy(i,j,k,5)=0.25_WP*sum(this%VISC(i-1:i,j-1:j,k))*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k))); FQx(i,j,k,6)=FQy(i,j,k,5)
                  FQz(i,j,k,3)=FQy(i,j,k,5)*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
               end if
               if (maxval(this%iSL(i,j-1:j,k-1:k)).eq.0) then
                  FQz(i,j,k,6)=0.25_WP*sum(this%VISC(i,j-1:j,k-1:k))*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k))); FQy(i,j,k,7)=FQz(i,j,k,6)
                  FQx(i,j,k,3)=FQz(i,j,k,6)*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
               end if
               if (maxval(this%iSL(i-1:i,j,k-1:k)).eq.0) then
                  FQx(i,j,k,7)=0.25_WP*sum(this%VISC(i-1:i,j,k-1:k))*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1))); FQz(i,j,k,5)=FQx(i,j,k,7)
                  FQy(i,j,k,3)=FQx(i,j,k,7)*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
               end if
            end do
         end do
      end do
      
      ! Assemble time derivative for conserved variables
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Viscous momentum transport
               dQdt(i,j,k,5)=dQdt(i,j,k,5)+this%dxi*(FQx(i  ,j,k,5)-FQx(i-1,j,k,5))+this%dyi*(FQy(i,j+1,k,5)-FQy(i,j  ,k,5))+this%dzi*(FQz(i,j,k+1,5)-FQz(i,j,k  ,5))
               dQdt(i,j,k,6)=dQdt(i,j,k,6)+this%dxi*(FQx(i+1,j,k,6)-FQx(i  ,j,k,6))+this%dyi*(FQy(i,j  ,k,6)-FQy(i,j-1,k,6))+this%dzi*(FQz(i,j,k+1,6)-FQz(i,j,k  ,6))
               dQdt(i,j,k,7)=dQdt(i,j,k,7)+this%dxi*(FQx(i+1,j,k,7)-FQx(i  ,j,k,7))+this%dyi*(FQy(i,j+1,k,7)-FQy(i,j  ,k,7))+this%dzi*(FQz(i,j,k  ,7)-FQz(i,j,k-1,7))
               ! Viscous heating term
               if (this%iSL(i,j,k).eq.0.and.this%VF(i,j,k).ge.0.5_WP) dQdt(i,j,k,3)=dQdt(i,j,k,3)+FQx(i,j,k,5)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+FQy(i,j,k,6)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+FQz(i,j,k,7)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+0.25_WP*sum(FQz(i:i+1,j:j+1,k,3))+0.25_WP*sum(FQx(i,j:j+1,k:k+1,3))+0.25_WP*sum(FQy(i:i+1,j,k:k+1,3))
               if (this%iSL(i,j,k).eq.0.and.this%VF(i,j,k).lt.0.5_WP) dQdt(i,j,k,4)=dQdt(i,j,k,4)+FQx(i,j,k,5)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+FQy(i,j,k,6)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+FQz(i,j,k,7)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+0.25_WP*sum(FQz(i:i+1,j:j+1,k,3))+0.25_WP*sum(FQx(i,j:j+1,k:k+1,3))+0.25_WP*sum(FQy(i:i+1,j,k:k+1,3))
            end do
         end do
      end do
      
      ! Deallocate flux arrays
      deallocate(FQx,FQy,FQz)
      
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
         integer :: ii,jj,kk
         real(IRL_double), dimension(0:2) :: normal
         real(IRL_double), dimension(0:188) :: moments
         integer :: direction,direction2
         logical :: flip
         real(IRL_double) :: m000,m100,m010,m001,temp
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
            flip=.false.; if (this%VF(i,j,k).ge.0.5_WP) flip=.true.
            m000=0.0_WP; m100=0.0_WP; m010=0.0_WP; m001=0.0_WP
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
            call reflect_moments(moments,center,direction,direction2)
            ! Get PLIC normal vector from neural network
            call get_normal(moments,normal)
            normal=normalize(normal)
            ! Rotate normal vector to original octant
            if (direction2.eq.1) then
               temp=normal(0)
               normal(0)=normal(1)
               normal(1)=temp
            else if (direction2.eq.2) then
               temp=normal(1)
               normal(1)=normal(2)
               normal(2)=temp
            else if (direction2.eq.3) then
               temp=normal(0)
               normal(0)=normal(2)
               normal(2)=temp
            else if (direction2.eq.4) then
               temp=normal(1)
               normal(1)=normal(2)
               normal(2)=temp
               temp=normal(0)
               normal(0)=normal(1)
               normal(1)=temp
            else if (direction2.eq.5) then
               temp=normal(0)
               normal(0)=normal(2)
               normal(2)=temp
               temp=normal(0)
               normal(0)=normal(1)
               normal(1)=temp
            end if
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
            ! Handle lower dimensions exactly
            if (this%cfg%nx.eq.1) normal(0)=0.0_WP
            if (this%cfg%ny.eq.1) normal(1)=0.0_WP
            if (this%cfg%nz.eq.1) normal(2)=0.0_WP
            normal=normalize(normal)
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
         integer :: i,j,k
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
         type(SepVM_type)   :: separated_volume_moments
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
         integer :: i,j,k
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
            this%RHOL(i,j,k)=this%Q(i,j,k,1)/this%VF(i,j,k)
            this%IL  (i,j,k)=this%Q(i,j,k,3)/this%Q (i,j,k,1)
            this%PL  (i,j,k)=this%getPL(this%RHOL(i,j,k),this%IL(i,j,k))
            CL              =this%getCL(this%RHOL(i,j,k),this%PL(i,j,k))
         else
            this%VF  (i,j,k)=0.0_WP
            this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
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
            this%VF  (i,j,k)=1.0_WP
            this%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            this%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
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
               Ui(i,j,k)=0.5_WP*sum(this%U(i:i+1,j,k))
               Vi(i,j,k)=0.5_WP*sum(this%V(i,j:j+1,k))
               Wi(i,j,k)=0.5_WP*sum(this%W(i,j,k:k+1))
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
         if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
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
         ! Update polygon for visualization
         call zeroPolygon(this%interface_polygon(i,j,k)); if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) call getPoly(cell,this%PLIC(i,j,k),0,this%interface_polygon(i,j,k))
      end do; end do; end do
   end subroutine apply_relax
   
   
   !> Get artifical bulk kinematic viscosity
   subroutine get_viscartif(this,dt,beta)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: beta
      integer :: i,j,k,si,sj,sk,n
      integer, parameter :: nfilter=1
      real(WP) :: max_beta,dudy,dudz,dvdx,dvdz,dwdx,dwdy,vort,grad_div
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP), parameter :: Cartif=2.0_WP
      real(WP), parameter :: Cartif_vort=100.0_WP
      real(WP), dimension(:,:,:), allocatable :: div
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]
      ! Calculate max beta permissible
      max_beta=max_cfl*min(this%dx**2,this%dy**2,this%dz**2)/(4.0_WP*dt)
      ! Zero out array
      beta=0.0_WP
      ! Compute velocity divergence
      allocate(div(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1; do j=this%cfg%jmino_,this%cfg%jmaxo_-1; do i=this%cfg%imino_,this%cfg%imaxo_-1
         div(i,j,k)=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
      end do; end do; end do
      call this%cfg%sync(div)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) div(this%cfg%imaxo,:,:)=div(this%cfg%imaxo-1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) div(:,this%cfg%jmaxo,:)=div(:,this%cfg%jmaxo-1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) div(:,:,this%cfg%kmaxo)=div(:,:,this%cfg%kmaxo-1)
      ! Compute artificial bulk viscosity based on gradU provided
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Only work in compression regions
         if (div(i,j,k).ge.0.0_WP) cycle
         ! Compute local vorticity
         dudy=0.25_WP*this%dyi*sum(this%U(i:i+1,j:j+1,k)-this%U(i:i+1,j-1:j,k))
         dudz=0.25_WP*this%dzi*sum(this%U(i:i+1,j,k:k+1)-this%U(i:i+1,j,k-1:k))
         dvdx=0.25_WP*this%dxi*sum(this%V(i:i+1,j:j+1,k)-this%V(i-1:i,j:j+1,k))
         dvdz=0.25_WP*this%dzi*sum(this%V(i,j:j+1,k:k+1)-this%V(i,j:j+1,k-1:k))
         dwdx=0.25_WP*this%dxi*sum(this%W(i:i+1,j,k:k+1)-this%W(i-1:i,j,k:k+1))
         dwdy=0.25_WP*this%dyi*sum(this%W(i,j:j+1,k:k+1)-this%W(i,j-1:j,k:k+1))
         vort=(dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2
         ! Compute |grad(div)|
         grad_div=max(abs(div(i+1,j,k)-div(i,j,k)),abs(div(i,j,k)-div(i-1,j,k)))*this%dx**2&
         &       +max(abs(div(i,j+1,k)-div(i,j,k)),abs(div(i,j,k)-div(i,j-1,k)))*this%dy**2&
         &       +max(abs(div(i,j,k+1)-div(i,j,k)),abs(div(i,j,k)-div(i,j,k-1)))*this%dz**2
         ! Estimate artificial kinematic viscosity using grad(div)
         beta(i,j,k)=Cartif*grad_div*div(i,j,k)**2/(div(i,j,k)**2+Cartif_vort*vort+1.0e-15_WP)
         ! Clip it so CFL<max_CFL
         beta(i,j,k)=min(beta(i,j,k),max_beta)
      end do; end do; end do
      call this%cfg%sync(beta)
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            beta(this%cfg%imin-1,:,:)=beta(this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) beta(this%cfg%imax+1,:,:)=beta(this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            beta(:,this%cfg%jmin-1,:)=beta(:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) beta(:,this%cfg%jmax+1,:)=beta(:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            beta(:,:,this%cfg%kmin-1)=beta(:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) beta(:,:,this%cfg%kmax+1)=beta(:,:,this%cfg%kmax)
      end if
      ! Filter beta
      do n=1,nfilter
         div=beta
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            beta(i,j,k)=0.0_WP
            do sk=-1,+1; do sj=-1,+1; do si=-1,+1
               beta(i,j,k)=beta(i,j,k)+filter(si)*filter(sj)*filter(sk)*div(i+si,j+sj,k+sk)
            end do; end do; end do
         end do; end do; end do
         call this%cfg%sync(beta)
         if (.not.this%cfg%xper) then
            if (this%cfg%iproc.eq.1)            beta(this%cfg%imin-1,:,:)=beta(this%cfg%imin,:,:)
            if (this%cfg%iproc.eq.this%cfg%npx) beta(this%cfg%imax+1,:,:)=beta(this%cfg%imax,:,:)
         end if
         if (.not.this%cfg%yper) then
            if (this%cfg%jproc.eq.1)            beta(:,this%cfg%jmin-1,:)=beta(:,this%cfg%jmin,:)
            if (this%cfg%jproc.eq.this%cfg%npy) beta(:,this%cfg%jmax+1,:)=beta(:,this%cfg%jmax,:)
         end if
         if (.not.this%cfg%zper) then
            if (this%cfg%kproc.eq.1)            beta(:,:,this%cfg%kmin-1)=beta(:,:,this%cfg%kmin)
            if (this%cfg%kproc.eq.this%cfg%npz) beta(:,:,this%cfg%kmax+1)=beta(:,:,this%cfg%kmax)
         end if
      end do
      ! Free up memory
      deallocate(div)
   end subroutine get_viscartif
   
   
   !> Get kinematic eddy viscosity using Vreman's model
   subroutine get_vreman(this,dt,visc)
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: visc
      real(WP), parameter :: Cs_ref=0.17_WP
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP) :: max_visc,A,B,C
      real(WP), dimension(1:3,1:3) :: beta,gradU
      real(WP), dimension(:,:,:), allocatable :: tmp
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]
      integer :: i,j,k,si,sj,sk,n
      integer, parameter :: nfilter=1
      ! Model constant is c=2.5*Cs_ref**2 - Vreman uses c=0.07 which corresponds to Cs_ref=0.17
      C=2.5_WP*Cs_ref**2
      ! Calculate max visc permissible
      max_visc=max_cfl*min(this%dx**2,this%dy**2,this%dz**2)/(4.0_WP*dt)
      ! Zero out array
      visc=0.0_WP
      ! Compute the eddy viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Compute velocity gradient tensor
         gradU(1,1)=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))
         gradU(2,1)=0.25_WP*this%dyi*sum(this%U(i:i+1,j:j+1,k)-this%U(i:i+1,j-1:j,k))
         gradU(3,1)=0.25_WP*this%dzi*sum(this%U(i:i+1,j,k:k+1)-this%U(i:i+1,j,k-1:k))
         gradU(1,2)=0.25_WP*this%dxi*sum(this%V(i:i+1,j:j+1,k)-this%V(i-1:i,j:j+1,k))
         gradU(2,2)=this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))
         gradU(3,2)=0.25_WP*this%dzi*sum(this%V(i,j:j+1,k:k+1)-this%V(i,j:j+1,k-1:k))
         gradU(1,3)=0.25_WP*this%dxi*sum(this%W(i:i+1,j,k:k+1)-this%W(i-1:i,j,k:k+1))
         gradU(2,3)=0.25_WP*this%dyi*sum(this%W(i,j:j+1,k:k+1)-this%W(i,j-1:j,k:k+1))
         gradU(3,3)=this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
         ! Compute A=gradu_ij*gradu_ij invariant
         A=sum(gradU**2)
         ! Compute beta_ij=dx_m*dx_m*gradu_mi*gradu_mj
         do sj=1,3; do si=1,3; beta(si,sj)=this%dx**2*gradU(1,si)*gradU(1,sj)+this%dy**2*gradU(2,si)*gradU(2,sj)+this%dz**2*gradU(3,si)*gradU(3,sj); end do; end do
         ! Compute B invariant
         B=beta(1,1)*beta(2,2)-beta(1,2)**2+beta(1,1)*beta(3,3)-beta(1,3)**2+beta(2,2)*beta(3,3)-beta(2,3)**2
         ! Assemble algebraic eddy viscosity model
         if (B.lt.1.0e-8_WP) then
            visc(i,j,k)=0.0_WP
         else
            visc(i,j,k)=C*sqrt(B/A)
         end if
         ! Clip it so CFL<max_CFL
         visc(i,j,k)=min(visc(i,j,k),max_visc)
      end do; end do; end do
      ! Synchronize visc
      call this%cfg%sync(visc)
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            visc(this%cfg%imin-1,:,:)=visc(this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) visc(this%cfg%imax+1,:,:)=visc(this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            visc(:,this%cfg%jmin-1,:)=visc(:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) visc(:,this%cfg%jmax+1,:)=visc(:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            visc(:,:,this%cfg%kmin-1)=visc(:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) visc(:,:,this%cfg%kmax+1)=visc(:,:,this%cfg%kmax)
      end if
      ! Filter visc
      do n=1,nfilter
         allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp=visc
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            visc(i,j,k)=0.0_WP
            do sk=-1,+1; do sj=-1,+1; do si=-1,+1
               visc(i,j,k)=visc(i,j,k)+filter(si)*filter(sj)*filter(sk)*tmp(i+si,j+sj,k+sk)
            end do; end do; end do
         end do; end do; end do
         call this%cfg%sync(visc)
         if (.not.this%cfg%xper) then
            if (this%cfg%iproc.eq.1)            visc(this%cfg%imin-1,:,:)=visc(this%cfg%imin,:,:)
            if (this%cfg%iproc.eq.this%cfg%npx) visc(this%cfg%imax+1,:,:)=visc(this%cfg%imax,:,:)
         end if
         if (.not.this%cfg%yper) then
            if (this%cfg%jproc.eq.1)            visc(:,this%cfg%jmin-1,:)=visc(:,this%cfg%jmin,:)
            if (this%cfg%jproc.eq.this%cfg%npy) visc(:,this%cfg%jmax+1,:)=visc(:,this%cfg%jmax,:)
         end if
         if (.not.this%cfg%zper) then
            if (this%cfg%kproc.eq.1)            visc(:,:,this%cfg%kmin-1)=visc(:,:,this%cfg%kmin)
            if (this%cfg%kproc.eq.this%cfg%npz) visc(:,:,this%cfg%kmax+1)=visc(:,:,this%cfg%kmax)
         end if
         deallocate(tmp)
      end do
   end subroutine get_vreman
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(mpcomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      integer  :: ierr
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
      maxvisc=maxval((this%VISC+this%BETA)/(this%Q(:,:,:,1)+this%Q(:,:,:,2))); call MPI_ALLREDUCE(MPI_IN_PLACE,maxvisc,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
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
         if (this%VF(i,j,k).ge.VFlo) then
            this%RHOLmin=min(this%RHOLmin,this%RHOL(i,j,k)); this%RHOLmax=max(this%RHOLmax,this%RHOL(i,j,k))
            this%ILmin  =min(this%ILmin  ,this%IL  (i,j,k)); this%ILmax  =max(this%ILmax  ,this%IL  (i,j,k))
            this%PLmin  =min(this%PLmin  ,this%PL  (i,j,k)); this%PLmax  =max(this%PLmax  ,this%PL  (i,j,k))
            this%TLmin  =min(this%TLmin  ,this%TL  (i,j,k)); this%TLmax  =max(this%TLmax  ,this%TL  (i,j,k))
         end if
         if (this%VF(i,j,k).le.VFhi) then
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
   
   
   !> Finalize mpcomp solver object
   subroutine finalize(this)
      implicit none
      class(mpcomp), intent(inout) :: this
      ! Nullify all function pointers
      nullify(this%getPL)
      nullify(this%getTL)
      nullify(this%getCL)
      nullify(this%getSL)
      nullify(this%getPG)
      nullify(this%getTG)
      nullify(this%getCG)
      nullify(this%getSG)
      nullify(this%relax)
      ! Deallocate allocatable arrays
      this%nQ=0
      if (allocated(this%VF))        deallocate(this%VF)
      if (allocated(this%VFold))     deallocate(this%VFold)
      if (allocated(this%BL))        deallocate(this%BL)
      if (allocated(this%BLold))     deallocate(this%BLold)
      if (allocated(this%BG))        deallocate(this%BG)
      if (allocated(this%BGold))     deallocate(this%BGold)
      if (allocated(this%Q))         deallocate(this%Q)
      if (allocated(this%Qold))      deallocate(this%Qold)
      if (allocated(this%U))         deallocate(this%U)
      if (allocated(this%V))         deallocate(this%V)
      if (allocated(this%W))         deallocate(this%W)
      if (allocated(this%RHOL))      deallocate(this%RHOL)
      if (allocated(this%RHOG))      deallocate(this%RHOG)
      if (allocated(this%RHOLold))   deallocate(this%RHOLold)
      if (allocated(this%RHOGold))   deallocate(this%RHOGold)
      if (allocated(this%IL))        deallocate(this%IL)
      if (allocated(this%IG))        deallocate(this%IG)
      if (allocated(this%ILold))     deallocate(this%ILold)
      if (allocated(this%IGold))     deallocate(this%IGold)
      if (allocated(this%PL))        deallocate(this%PL)
      if (allocated(this%PG))        deallocate(this%PG)
      if (allocated(this%PLold))     deallocate(this%PLold)
      if (allocated(this%PGold))     deallocate(this%PGold)
      if (allocated(this%TL))        deallocate(this%TL)
      if (allocated(this%TG))        deallocate(this%TG)
      if (allocated(this%C))         deallocate(this%C)
      if (allocated(this%VISC))      deallocate(this%VISC)
      if (allocated(this%BETA))      deallocate(this%BETA)
      if (allocated(this%DIFF))      deallocate(this%DIFF)
      if (allocated(this%Qmin))      deallocate(this%Qmin)
      if (allocated(this%Qmax))      deallocate(this%Qmax)
      if (allocated(this%Qint))      deallocate(this%Qint)
      if (allocated(this%SLdQ))      deallocate(this%SLdQ)
      if (allocated(this%iSL))       deallocate(this%iSL)
      ! Deallocate IRL data - should be automatically deleted by C++
      if (allocated(this%PLIC))      deallocate(this%PLIC)
      if (allocated(this%PLICold))   deallocate(this%PLICold)
      if (allocated(this%localizer)) deallocate(this%localizer)
      if (allocated(this%localized_separator_link)) deallocate(this%localized_separator_link)
      if (allocated(this%localizer_link)) deallocate(this%localizer_link)
      if (allocated(this%interface_polygon)) deallocate(this%interface_polygon)
      ! Nullify config pointer
      nullify(this%cfg)
      ! Finalize timers
      call this%trhs%finalize()
      call this%tsl%finalize()
      call this%tplic%finalize()
   end subroutine finalize
   
   
end module mpcomp_class
