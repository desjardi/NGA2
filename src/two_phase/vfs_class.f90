!> Volume fraction solver class:
!> Provides support for various BC, semi-Lagrangian geometric advancement,
!> curvature calculation, interface reconstruction.
module vfs_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   use irl_fortran_interface
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: vfs,bcond
   
   ! Also expose the min and max VF values
   public :: VFhi,VFlo
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available interface reconstructions schemes for VF
   integer, parameter, public :: lvira=1             !< LVIRA scheme
   integer, parameter, public :: elvira=2            !< ELVIRA scheme
   integer, parameter, public :: mof=3               !< MOF scheme
   integer, parameter, public :: wmof=4              !< Wide-MOF scheme
   integer, parameter, public :: r2p=5               !< R2P scheme
   integer, parameter, public :: swartz=6            !< Swartz scheme
   integer, parameter, public :: youngs=7            !< Youngs' scheme
   integer, parameter, public :: lvlset=8            !< Levelset-based scheme
   integer, parameter, public :: plicnet=9           !< PLICnet
   integer, parameter, public :: r2pnet=10           !< R2Pnet
   
   ! List of available interface transport schemes for VF
   integer, parameter, public :: flux=1             !< Flux-based geometric transport
   integer, parameter, public :: flux_storage=2     !< Flux-based geometric transport with storage of detailed face fluxes
   integer, parameter, public :: remap=3            !< Cell-based geometric transport (faster but fluxes are not available)
   integer, parameter, public :: remap_storage=4    !< Cell-based geometric transport with storage of detailed volume moments
   
   ! IRL cutting moment calculation method
   integer, parameter, public :: recursive_simplex=0 !< Recursive simplex cutting
   integer, parameter, public :: half_edge=1         !< Half-edge cutting (default)
   integer, parameter, public :: nonrecurs_simplex=2 !< Non-recursive simplex cutting
   
   ! Default parameters for volume fraction solver
   integer,  parameter :: nband=3                                 !< Number of cells around the interfacial cells on which localized work is performed
   integer,  parameter :: advect_band=1                           !< How far we do the transport
   integer,  parameter :: distance_band=2                         !< How far we build the distance
   integer,  parameter :: max_interface_planes=2                  !< Maximum number of interfaces allowed (2 for R2P)
   real(WP), parameter :: VFlo=1.0e-12_WP                         !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo                        !< Maximum VF value considered
   real(WP), parameter :: volume_epsilon_factor =1.0e-15_WP       !< Minimum volume  to consider for computational geometry (normalized by min_meshsize**3)
   real(WP), parameter :: surface_epsilon_factor=1.0e-15_WP       !< Minimum surface to consider for computational geometry (normalized by min_meshsize**2)
   real(WP), parameter :: iterative_distfind_tol=1.0e-12_WP       !< Tolerance for iterative plane distance finding
   
   !> Boundary conditions for the volume fraction solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Volume fraction solver object definition
   type :: vfs
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_VFS'     !< Solver name (default=UNNAMED_VFS)
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Volume fraction data
      real(WP), dimension(:,:,:), allocatable :: VF       !< VF array
      real(WP), dimension(:,:,:), allocatable :: VFold    !< VFold array
      
      ! Phase barycenter data
      real(WP), dimension(:,:,:,:), allocatable :: Lbary  !< Liquid barycenter
      real(WP), dimension(:,:,:,:), allocatable :: Gbary  !< Gas barycenter
      
      ! Subcell phasic volume fields
      real(WP), dimension(:,:,:,:,:,:), allocatable :: Lvol   !< Subcell liquid volume
      real(WP), dimension(:,:,:,:,:,:), allocatable :: Gvol   !< Subcell gas volume
      
      ! Surface density data
      real(WP), dimension(:,:,:), allocatable :: SD       !< Surface density array
      
      ! Distance level set
      real(WP) :: Gclip                                   !< Min/max distance
      real(WP), dimension(:,:,:), allocatable :: G        !< Distance level set array
      
      ! Curvature
      real(WP), dimension(:,:,:), allocatable :: curv     !< Interface mean curvature
      real(WP), dimension(:,:,:,:), allocatable :: curv2p !< Curvature for each interface
      
      ! Band strategy
      integer, dimension(:,:,:), allocatable :: band      !< Band to localize workload around the interface
      integer, dimension(:,:),   allocatable :: band_map  !< Unstructured band mapping
      integer, dimension(0:nband) :: band_count           !< Number of cells per band value
      
      ! Interface handling methods
      integer :: reconstruction_method                    !< Interface reconstruction method
      integer :: transport_method                         !< Interface transport method
      logical :: cons_correct=.true.                      !< Conservative correction (true by default)
      
      ! Flotsam removal parameter
      real(WP) :: flotsam_thld=0.0_WP                     !< Threshold VF parameter for flotsam removal (0.0=off by default)

      ! Parameters for SGS modeling of thin structures
      logical  :: two_planes                              !< Whether we're using a 2-plane reconstruction approach
      real(WP) :: twoplane_thld1=0.99_WP                  !< Average normal magnitude threshold for r2p to switch from one-plane to two-planes (purely local)
      real(WP) :: twoplane_thld2=0.5_WP                   !< Average normal magnitude threshold above which r2p switches to LVIRA (based on 3x3x3 stencil)
      real(WP) :: thin_thld_dotprod=-0.5_WP               !< Maximum dot product of two interface normals for their respective cells to be considered thin region cells
      real(WP) :: thin_thld_max=0.8_WP                    !< Maximum local thickness to be considered a thin region cell (as a factor of mesh size)
      real(WP) :: thin_thld_min=0.0_WP                    !< Minimum local thickness for thin structure removal (0.0=off, as a factor of mesh size)
      real(WP) :: edge_thld=0.80_WP                       !< Threshold for classifying edges
      
      ! Interface sensing variables
      real(WP), dimension(:,:,:), allocatable :: thin_sensor     !< Thin structure sensing (=1 is liquid, =2 is gas)
      real(WP), dimension(:,:,:), allocatable :: thickness       !< Local thickness of thin region
      real(WP), dimension(:,:,:), allocatable :: edge_sensor     !< Edge sensing (higher is edge)
      real(WP), dimension(:,:,:,:), allocatable :: edge_normal   !< Edge normal
      
      ! Curvature clipping parameter
      real(WP) :: maxcurv_times_mesh=1.0_WP               !< Clipping parameter for maximum curvature (classically set to 1, but is larger with r2p since we resolve more)
      
      ! IRL objects
      type(ByteBuffer_type) :: send_byte_buffer
      type(ByteBuffer_type) :: recv_byte_buffer
      type(ObjServer_PlanarSep_type)  :: planar_separator_allocation
      type(ObjServer_PlanarLoc_type)  :: planar_localizer_allocation
      type(ObjServer_LocSepLink_type) :: localized_separator_link_allocation
      type(ObjServer_LocLink_type)    :: localizer_link_allocation
      type(PlanarLoc_type),   dimension(:,:,:),   allocatable :: localizer
      type(PlanarSep_type),   dimension(:,:,:),   allocatable :: liquid_gas_interface
      type(LocSepLink_type),  dimension(:,:,:),   allocatable :: localized_separator_link
      type(ListVM_VMAN_type), dimension(:,:,:),   allocatable :: triangle_moments_storage
      type(LocLink_type),     dimension(:,:,:),   allocatable :: localizer_link
      type(Poly_type),        dimension(:,:,:,:), allocatable :: interface_polygon
      type(Poly_type),        dimension(:,:,:,:), allocatable :: polyface
      type(SepVM_type),       dimension(:,:,:,:), allocatable :: face_flux    !< Only stored if flux-based transport is used
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for enforcing bconds
      integer, dimension(:,:,:), allocatable :: vmask     !< Integer array used for enforcing bconds - for vertices
      
      ! Monitoring quantities
      real(WP) :: VFmax,VFmin,VFint,SDint                 !< Maximum, minimum, and integral volume fraction and surface density
      real(WP) :: flotsam_error                           !< Integral of flotsam removal error
      real(WP) :: thinstruct_error                        !< Integral of thin structure removal error
      
      ! Additional storage option for geometric transport of auxiliary quantities
      type(TagAccVM_SepVM_type), dimension(:,:,:,:), allocatable :: detailed_face_flux    !< Cell-decomposed face flux geometric data
      type(TagAccVM_SepVM_type), dimension(:,:,:),   allocatable :: detailed_remap        !< Cell-decomposed remapped cell geometric data
      
      ! Old arrays that are needed for the compressible MAST solver
      real(WP), dimension(:,:,:,:), allocatable :: Lbaryold  !< Liquid barycenter
      real(WP), dimension(:,:,:,:), allocatable :: Gbaryold  !< Gas barycenter
      type(PlanarSep_type),  dimension(:,:,:), allocatable :: liquid_gas_interfaceold
      type(LocSepLink_type), dimension(:,:,:), allocatable :: localized_separator_linkold
      type(ObjServer_PlanarSep_type)  :: planar_separatorold_allocation
      type(ObjServer_LocSepLink_type) :: localized_separator_linkold_allocation
      
   contains
      procedure :: initialize                             !< Initialize the vfs object
      procedure :: print=>vfs_print                       !< Output solver to the screen
      procedure :: initialize_irl                         !< Initialize the IRL objects
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: update_band                            !< Update the band info given the VF values
      procedure :: sync_interface                         !< Synchronize the IRL objects
      procedure :: remove_flotsams                        !< Remove flotsams manually
      procedure :: remove_thinstruct                      !< Remove thin structures manually
      procedure :: clean_irl_and_band                     !< After a manual change in VF (maybe due to transfer to drops), update IRL and band info
      procedure :: sync_and_clean_barycenters             !< Synchronize and clean up phasic barycenters
      procedure, private :: sync_side                     !< Synchronize the IRL objects across one side - another I/O helper
      procedure, private :: sync_ByteBuffer               !< Communicate byte packets across one side - another I/O helper
      procedure, private :: calculate_offset_to_planes    !< Helper routine for I/O
      procedure, private :: crude_phase_test              !< Helper function that rapidly assess if a mixed cell might be present
      procedure :: project                                !< Function that performs a Lagrangian projection of a vertex
      procedure :: read_interface                         !< Read an IRL interface from a file
      procedure :: write_interface                        !< Write an IRL interface to a file
      procedure :: advance                                !< Advance VF to next step
      procedure :: transport_flux                         !< Transport VF using geometric fluxing
      procedure :: transport_flux_storage                 !< Transport VF using geometric fluxing with storage
      procedure :: transport_remap                        !< Transport VF using geometric cell remap
      procedure :: transport_remap_storage                !< Transport VF using geometric cell remap with storage
      procedure :: advect_interface                       !< Advance IRL surface to next step
      procedure :: build_interface                        !< Reconstruct IRL interface from VF field
      procedure :: build_elvira                           !< ELVIRA reconstruction of the interface from VF field
      procedure :: build_lvira                            !< LVIRA reconstruction of the interface from VF field
      procedure :: build_mof                              !< MOF reconstruction of the interface from VF field
      procedure :: build_wmof                             !< Wide-MOF reconstruction of the interface from VF field
      procedure :: build_r2p                              !< R2P reconstruction of the interface from VF field
      procedure :: build_plicnet                          !< PLICnet reconstruction of the interface from VF and bary fields
      procedure :: build_r2pnet                           !< R2Pnet reconstruction of the interface
      procedure :: sense_interface                        !< Calculate various surface sensors
      procedure :: get_thickness                          !< Calculate multiphasic structure thickness
      procedure :: detect_thin_regions                    !< Detect thin regions
      procedure :: detect_edge_regions                    !< Detect edge regions
      procedure :: build_youngs                           !< Youngs' reconstruction of the interface from VF field
      !procedure :: build_lvlset                           !< LVLSET-based reconstruction of the interface from VF field
      procedure :: smooth_interface                       !< Interface smoothing based on Swartz idea
      procedure :: set_full_bcond                         !< Full liq/gas plane-setting for boundary cells - this is stair-stepped
      procedure :: polygonalize_interface                 !< Build a discontinuous polygonal representation of the IRL interface
      procedure :: distance_from_polygon                  !< Build a signed distance field from the polygonalized interface
      procedure :: subcell_vol                            !< Build subcell phasic volumes from reconstructed interface
      procedure :: reset_volume_moments                   !< Reconstruct volume moments from IRL interfaces
      procedure :: reset_moments                          !< Reconstruct first-order moments from IRL interfaces
      procedure :: update_surfmesh                        !< Update a surfmesh object using current polygons
      procedure :: update_surfmesh_nowall                 !< Update a surfmesh object using current polygons - do not show polygons in walls
      procedure :: get_curvature                          !< Compute curvature from IRL surface polygons
      procedure :: paraboloid_fit                         !< Perform local paraboloid fit of IRL surface using IRL barycenter data
      procedure :: paraboloid_integral_fit                !< Perform local paraboloid fit of IRL surface using surface-integrated IRL data
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: get_cfl                                !< Get CFL for the VF solver

      procedure :: allocate_supplement                    !< Allocate arrays and initialize objects needed for compressible MAST solver
      procedure :: copy_interface_to_old                  !< Copy interface variables at beginning of timestep
      procedure :: fluxpoly_project_getmoments            !< Project face to get flux volume and output its moments
      procedure :: fluxpoly_cell_getvolcentr              !< Get the volume and centroid during generalized SL advection

   end type vfs
   
   
contains
   
   
   !> Initialization for volume fraction solver
   subroutine initialize(this,cfg,reconstruction_method,transport_method,name)
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: reconstruction_method
      integer, intent(in), optional :: transport_method
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Set transport scheme
      if (present(transport_method)) then
         this%transport_method=transport_method
      else
         this%transport_method=flux ! Set flux-based geometric transport as default (could change at some point)
      end if
      
      ! Check that we have at least 3 overlap cells - we can push that to 2 with limited work!
      if (cfg%no.lt.3) call die('[vfs initialize] The config requires at least 3 overlap cells')
      
      ! Point to pgrid object
      this%cfg=>cfg
      
      ! Nullify bcond list
      this%nbc=0
      this%first_bc=>NULL()
      
      ! Allocate variables
      allocate(this%VF   (  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VF   =0.0_WP
      allocate(this%VFold(  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VFold=0.0_WP
      allocate(this%Lbary(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Lbary=0.0_WP
      allocate(this%Gbary(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Gbary=0.0_WP
      allocate(this%SD   (  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%SD   =0.0_WP
      allocate(this%G    (  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%G    =0.0_WP
      allocate(this%curv (  this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%curv =0.0_WP
      
      ! Set clipping distance
      this%Gclip=real(distance_band+1,WP)*this%cfg%min_meshsize
      
      ! Subcell phasic volumes
      allocate(this%Lvol(0:1,0:1,0:1,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Lvol=0.0_WP
      allocate(this%Gvol(0:1,0:1,0:1,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Gvol=0.0_WP
      
      ! Prepare the band arrays
      allocate(this%band(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%band=0
      if (allocated(this%band_map)) deallocate(this%band_map)
      
      ! Set reconstruction method
      select case (reconstruction_method)
      case (lvira,elvira,swartz,youngs,mof,wmof,plicnet)
         this%reconstruction_method=reconstruction_method
         this%two_planes=.false.
      case (r2p,r2pnet)
         this%reconstruction_method=reconstruction_method
         ! Allocate extra curvature storage
         this%two_planes=.true.
         allocate(this%curv2p(1:2,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%curv2p=0.0_WP
         ! Allocate extra sensors
         allocate(this%thin_sensor(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%thin_sensor=0.0_WP
         allocate(this%thickness  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%thickness  =0.0_WP
         allocate(this%edge_sensor(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%edge_sensor=0.0_WP
         allocate(this%edge_normal(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%edge_normal=0.0_WP
         ! By default, use thin structure removal
         this%thin_thld_min=1.0e-4_WP !< This removes any thin structure with thickness below dx/1000
         ! By default, use flotsam removal
         this%flotsam_thld=1.0e-3_WP  !< This considers any separated structure around dx/10 and below as bogus
         ! Also allow for larger curvatures to be calculated
         this%maxcurv_times_mesh=2.0_WP
      case default
         call die('[vfs initialize] Unknown interface reconstruction scheme.')
      end select
      
      ! Initialize IRL
      call this%initialize_irl()
      
      ! Prepare mask for VF
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
      
      ! Prepare mask for vertices
      allocate(this%vmask(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%vmask=0
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.           1) this%vmask(               :this%cfg%imin,:,:)=2
         if (this%cfg%iproc.eq.this%cfg%npx) this%vmask(this%cfg%imax+1:             ,:,:)=2
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.           1) this%vmask(:,               :this%cfg%jmin,:)=2
         if (this%cfg%jproc.eq.this%cfg%npy) this%vmask(:,this%cfg%jmax+1:             ,:)=2
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.           1) this%vmask(:,:,               :this%cfg%kmin)=2
         if (this%cfg%kproc.eq.this%cfg%npz) this%vmask(:,:,this%cfg%kmax+1:             )=2
      end if
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (minval(this%cfg%VF(i-1:i,j-1:j,k-1:k)).eq.0.0_WP) this%vmask(i,j,k)=1
            end do
         end do
      end do
      call this%cfg%sync(this%vmask)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) this%vmask(this%cfg%imino,:,:)=this%vmask(this%cfg%imino+1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) this%vmask(:,this%cfg%jmino,:)=this%vmask(:,this%cfg%jmino+1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) this%vmask(:,:,this%cfg%kmino)=this%vmask(:,:,this%cfg%kmino+1)
      
   end subroutine initialize
    

   !> Initialize the IRL representation of our interfaces
   subroutine initialize_irl(this)
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
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
      
      ! Allocate IRL arrays
      allocate(this%localizer               (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%liquid_gas_interface    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%localized_separator_link(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%triangle_moments_storage(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%localizer_link          (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%interface_polygon(1:max_interface_planes,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%polyface            (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Work arrays for flux-based transport
      select case (this%transport_method)
      case (flux)
         ! Arrays for storing face fluxes
         allocate(this%face_flux(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  do n=1,3
                     call new(this%face_flux(n,i,j,k))
                  end do
               end do
            end do
         end do
      case (flux_storage)
         ! Arrays for storing face fluxes and detailed flux geometry
         allocate(this%face_flux         (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%detailed_face_flux(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  do n=1,3
                     call new(this%face_flux(n,i,j,k))
                     call new(this%detailed_face_flux(n,i,j,k))
                  end do
               end do
            end do
         end do
      case (remap)
         ! No storage needed
      case (remap_storage)
         ! Array for storing detailed remapped cell geometry
         allocate(this%detailed_remap(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call new(this%detailed_remap(i,j,k))
               end do
            end do
         end do
      case default
         call die('[vfs initialize IRL] Unknown transport method')
      end select
      
      ! Precomputed face correction velocities
      !> allocate(face_correct_velocity(1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
      
      ! Initialize size for IRL
      total_cells=int(this%cfg%nxo_,8)*int(this%cfg%nyo_,8)*int(this%cfg%nzo_,8)
      call new(this%planar_localizer_allocation,total_cells)
      call new(this%planar_separator_allocation,total_cells)
      call new(this%localized_separator_link_allocation,total_cells)
      call new(this%localizer_link_allocation,total_cells)
      
      ! Initialize arrays and setup linking
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Transfer cell to IRL
               call new(this%localizer(i,j,k),this%planar_localizer_allocation)
               call setFromRectangularCuboid(this%localizer(i,j,k),[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               ! PLIC interface(s)
               call new(this%liquid_gas_interface(i,j,k),this%planar_separator_allocation)
               ! PLIC+mesh with connectivity (i.e., link)
               call new(this%localized_separator_link(i,j,k),this%localized_separator_link_allocation,this%localizer(i,j,k),this%liquid_gas_interface(i,j,k))
               ! For transport surface
               call new(this%triangle_moments_storage(i,j,k))
               ! Mesh with connectivity
               call new(this%localizer_link(i,j,k),this%localizer_link_allocation,this%localizer(i,j,k))
            end do
         end do
      end do
      
      ! Polygonal representation of the surface
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               do n=1,max_interface_planes
                  call new(this%interface_polygon(n,i,j,k))
               end do
            end do
         end do
      end do
      
      ! Polygonal representation of cell faces
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Polygonal representation of the x-face
               call new(this%polyface(1,i,j,k))
               vert(:,1)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k  )]
               vert(:,2)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k  )]
               vert(:,3)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k+1)]
               vert(:,4)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k+1)]
               call construct(this%polyface(1,i,j,k),4,vert)
               call setPlaneOfExistence(this%polyface(1,i,j,k),[1.0_WP,0.0_WP,0.0_WP,this%cfg%x(i)])
               ! Polygonal representation of the y-face
               call new(this%polyface(2,i,j,k))
               vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k  )]
               vert(:,2)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k+1)]
               vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k+1)]
               vert(:,4)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k  )]
               call construct(this%polyface(2,i,j,k),4,vert)
               call setPlaneOfExistence(this%polyface(2,i,j,k),[0.0_WP,1.0_WP,0.0_WP,this%cfg%y(j)])
               ! Polygonal representation of the z-face
               call new(this%polyface(3,i,j,k))
               vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k)]
               vert(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k)]
               vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k)]
               vert(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k)]
               call construct(this%polyface(3,i,j,k),4,vert)
               call setPlaneOfExistence(this%polyface(3,i,j,k),[0.0_WP,0.0_WP,1.0_WP,this%cfg%z(k)])
            end do
         end do
      end do
      
      ! Give each link a unique lexicographic tag (per processor)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               tag=this%cfg%get_lexico_from_ijk([i,j,k])
               call setId(this%localized_separator_link(i,j,k),tag)
               call setId(this%localizer_link(i,j,k),tag)
            end do
         end do
      end do
      
      ! Set the connectivity
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
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
            end do
         end do
      end do
      
      ! Prepare byte storage for synchronization
      call new(this%send_byte_buffer)
      call new(this%recv_byte_buffer)
      
   end subroutine initialize_irl
   

   ! Set up additional arrays needed for the compressible MAST solver
   subroutine allocate_supplement(this)
     implicit none
     class(vfs), intent(inout) :: this
     integer  :: i,j,k,tag
     integer(IRL_LargeOffsetIndex_t) :: total_cells
     
     ! Allocate arrays for storing old variables
     allocate(this%liquid_gas_interfaceold    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
     allocate(this%localized_separator_linkold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
     total_cells=int(this%cfg%nxo_,8)*int(this%cfg%nyo_,8)*int(this%cfg%nzo_,8)
     call new(this%planar_separatorold_allocation,total_cells)
     call new(this%localized_separator_linkold_allocation,total_cells)

     ! Initialize arrays and setup linking
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              ! PLIC interface(s)
              call new(this%liquid_gas_interfaceold(i,j,k),this%planar_separatorold_allocation)
              ! PLIC+mesh with connectivity (i.e., link)
              call new(this%localized_separator_linkold(i,j,k),this%localized_separator_linkold_allocation,this%localizer(i,j,k),this%liquid_gas_interfaceold(i,j,k))
           end do
        end do
     end do

     ! Give each link a unique lexicographic tag (per processor)
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              tag=this%cfg%get_lexico_from_ijk([i,j,k])
              call setId(this%localized_separator_linkold(i,j,k),tag)
           end do
        end do
     end do

     ! Set the connectivity
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              ! In the x- direction
              if (i.gt.this%cfg%imino_) then
                 call setEdgeConnectivity(this%localized_separator_linkold(i,j,k),0,this%localized_separator_linkold(i-1,j,k))
              end if
              ! In the x+ direction
              if (i.lt.this%cfg%imaxo_) then
                 call setEdgeConnectivity(this%localized_separator_linkold(i,j,k),1,this%localized_separator_linkold(i+1,j,k))
              end if
              ! In the y- direction
              if (j.gt.this%cfg%jmino_) then
                 call setEdgeConnectivity(this%localized_separator_linkold(i,j,k),2,this%localized_separator_linkold(i,j-1,k))
              end if
              ! In the y+ direction
              if (j.lt.this%cfg%jmaxo_) then
                 call setEdgeConnectivity(this%localized_separator_linkold(i,j,k),3,this%localized_separator_linkold(i,j+1,k))
              end if
              ! In the z- direction
              if (k.gt.this%cfg%kmino_) then
                 call setEdgeConnectivity(this%localized_separator_linkold(i,j,k),4,this%localized_separator_linkold(i,j,k-1))
              end if
              ! In the z+ direction
              if (k.lt.this%cfg%kmaxo_) then
                 call setEdgeConnectivity(this%localized_separator_linkold(i,j,k),5,this%localized_separator_linkold(i,j,k+1))
              end if
           end do
        end do
     end do


     ! Deallocate memory-intensive arrays that are not used
     deallocate(this%face_flux)

   end subroutine allocate_supplement
   
   
   ! Store old liquid_gas_interface at the beginning of every timestep
   subroutine copy_interface_to_old(this)
     implicit none
     class(vfs), intent(inout) :: this
     integer :: i,j,k

     ! Copy interface
     do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
           do i=this%cfg%imino_,this%cfg%imaxo_
              call copy(this%liquid_gas_interfaceold(i,j,k),this%liquid_gas_interface(i,j,k))
           end do
        end do
     end do

     ! Copy volume fraction
     this%VFold = this%VF

     ! Copy barycenters
     this%Gbaryold = this%Gbary
     this%Lbaryold = this%Lbary

     return
   end subroutine copy_interface_to_old

   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,         only: lowercase
      use messager,       only: die
      use iterator_class, only: locator_ftype
      implicit none
      class(vfs), intent(inout) :: this
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
         case default; call die('[vfs add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[vfs apply_bcond] Neumann requires a direction')
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
         call die('[vfs apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[vfs get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
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
                  this%VF(i,j,k)=this%VF(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[vfs apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%VF)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the new VF based on U/V/W and dt
   subroutine advance(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      
      ! First perform transport
      select case (this%transport_method)
      case (flux)
         call this%transport_flux(dt,U,V,W)
      case (flux_storage)
         call this%transport_flux_storage(dt,U,V,W)
      case (remap)
         call this%transport_remap(dt,U,V,W)
      case (remap_storage)
         call this%transport_remap_storage(dt,U,V,W)
      end select
      
      ! Advect interface polygons
      call this%advect_interface(dt,U,V,W)
      
      ! Remove flotsams and thin structures if needed
      call this%remove_flotsams()
      call this%remove_thinstruct()
      
      ! Synchronize and clean-up barycenter fields
      call this%sync_and_clean_barycenters()
      
      ! Update the band
      call this%update_band()
      
      ! Perform interface reconstruction from transported moments
      call this%build_interface()
      
      ! Create discontinuous polygon mesh from IRL interface
      call this%polygonalize_interface()
      
      ! Perform interface sensing
      if (this%two_planes) call this%sense_interface()
      
      ! Calculate distance from polygons
      call this%distance_from_polygon()
      
      ! Calculate subcell phasic volumes
      call this%subcell_vol()
      
      ! Calculate curvature
      call this%get_curvature()
      
      ! Reset moments to guarantee compatibility with interface reconstruction
      call this%reset_volume_moments()
      
   end subroutine advance
   

   !> Perform cell-based transport of VF based on U/V/W and dt
   subroutine transport_remap(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,index
      real(IRL_double), dimension(3,14) :: cell
      real(IRL_double), dimension(3, 9) :: face
      type(Poly24_type) :: remap_cell
      type(CapDod_type) :: remap_face
      real(WP) :: vol_now,crude_VF
      real(WP) :: lvol,gvol
      real(WP), dimension(3,2) :: bounding_pts
      integer, dimension(3,2) :: bb_indices
      type(SepVM_type) :: my_SepVM
      
      ! Allocate poly24 and capdod, as well as SepVM objects
      call new(remap_cell)
      call new(remap_face)
      call new(my_SepVM)
      
      ! Loop over the advection band and compute conservative cell-based remap using semi-Lagrangian algorithm
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Construct and project cell
         cell(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; if (this%vmask(i+1,j  ,k+1).ne.1) cell(:,1)=this%project(cell(:,1),i,j,k,-dt,U,V,W)
         cell(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; if (this%vmask(i+1,j  ,k  ).ne.1) cell(:,2)=this%project(cell(:,2),i,j,k,-dt,U,V,W)
         cell(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; if (this%vmask(i+1,j+1,k  ).ne.1) cell(:,3)=this%project(cell(:,3),i,j,k,-dt,U,V,W)
         cell(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; if (this%vmask(i+1,j+1,k+1).ne.1) cell(:,4)=this%project(cell(:,4),i,j,k,-dt,U,V,W)
         cell(:,5)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; if (this%vmask(i  ,j  ,k+1).ne.1) cell(:,5)=this%project(cell(:,5),i,j,k,-dt,U,V,W)
         cell(:,6)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; if (this%vmask(i  ,j  ,k  ).ne.1) cell(:,6)=this%project(cell(:,6),i,j,k,-dt,U,V,W)
         cell(:,7)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; if (this%vmask(i  ,j+1,k  ).ne.1) cell(:,7)=this%project(cell(:,7),i,j,k,-dt,U,V,W)
         cell(:,8)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; if (this%vmask(i  ,j+1,k+1).ne.1) cell(:,8)=this%project(cell(:,8),i,j,k,-dt,U,V,W)
         
         ! Correct volume of x- face
         face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=cell(:,5)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,6)
         face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=cell(:,7)
         face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,8)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
         cell(:,14)=getPt(remap_face,8)
         
         ! Correct volume of x+ face
         face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=cell(:,1)
         face(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,2)
         face(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=cell(:,3)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,4)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*U(i+1,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
         cell(:, 9)=getPt(remap_face,8)
         
         ! Correct volume of y- face
         face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=cell(:,2)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,6)
         face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=cell(:,5)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=cell(:,1)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
         cell(:,10)=getPt(remap_face,8)
         
         ! Correct volume of y+ face
         face(:,1)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=cell(:,3)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,6)=cell(:,7)
         face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,7)=cell(:,8)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,4)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*V(i,j+1,k)*this%cfg%dx(i)*this%cfg%dz(k))
         cell(:,12)=getPt(remap_face,8)
         
         ! Correct volume of z- face
         face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=cell(:,7)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,6)
         face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=cell(:,2)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=cell(:,3)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
         cell(:,11)=getPt(remap_face,8)
         
         ! Correct volume of z+ face
         face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,5)=cell(:,8)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,6)=cell(:,5)
         face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=cell(:,1)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,4)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*W(i,j,k+1)*this%cfg%dx(i)*this%cfg%dy(j))
         cell(:,13)=getPt(remap_face,8)
         
         ! Form remapped cell in IRL
         call construct(remap_cell,cell)
         
         ! Get bounding box for our remapped cell
         call getBoundingPts(remap_cell,bounding_pts(:,1),bounding_pts(:,2))
         bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
         bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
         
         ! Crudely check phase information for remapped cell and skip cells where nothing is changing
         crude_VF=this%crude_phase_test(bb_indices)
         if (crude_VF.ge.0.0_WP) cycle
         
         ! Need full geometric flux
         call getMoments(remap_cell,this%localized_separator_link(i,j,k),my_SepVM)

         ! Compute new liquid volume fraction
         lvol=getVolumePtr(my_SepVM,0)
         gvol=getVolumePtr(my_SepVM,1)
         this%VF(i,j,k)=lvol/(lvol+gvol)
         
         ! Only work on higher order moments if VF is in [VFlo,VFhi]
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
         else
            ! Get old phasic barycenters
            this%Lbary(:,i,j,k)=getCentroidPtr(my_SepVM,0)/lvol
            this%Gbary(:,i,j,k)=getCentroidPtr(my_SepVM,1)/gvol
            ! Project then forward in time
            this%Lbary(:,i,j,k)=this%project(this%Lbary(:,i,j,k),i,j,k,dt,U,V,W)
            this%Gbary(:,i,j,k)=this%project(this%Gbary(:,i,j,k),i,j,k,dt,U,V,W)
         end if
         
      end do

      ! Synchronize VF and barycenter fields
      call this%cfg%sync(this%VF)
      call this%sync_and_clean_barycenters()
      
   end subroutine transport_remap
   
   
   !> Perform cell-based transport of VF based on U/V/W and dt
   !> Include storage of detailed geometry of remapped cell
   subroutine transport_remap_storage(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,index,n
      real(IRL_double), dimension(3,14) :: cell
      real(IRL_double), dimension(3, 9) :: face
      type(Poly24_type) :: remap_cell
      type(CapDod_type) :: remap_face
      real(WP) :: vol_now,crude_VF
      real(WP) :: lvol,gvol
      real(WP), dimension(3,2) :: bounding_pts
      integer, dimension(3,2) :: bb_indices
      real(WP), dimension(3) :: lbar,gbar
      type(SepVM_type) :: my_SepVM
      
      ! Clean up detailed remap data
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               call clear(this%detailed_remap(i,j,k))
            end do
         end do
      end do

      ! Allocate poly24 and capdod, as well as SepVM objects
      call new(remap_cell)
      call new(remap_face)
      
      ! Loop over the advection band and compute conservative cell-based remap using semi-Lagrangian algorithm
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Construct and project cell
         cell(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; if (this%vmask(i+1,j  ,k+1).ne.1) cell(:,1)=this%project(cell(:,1),i,j,k,-dt,U,V,W)
         cell(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; if (this%vmask(i+1,j  ,k  ).ne.1) cell(:,2)=this%project(cell(:,2),i,j,k,-dt,U,V,W)
         cell(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; if (this%vmask(i+1,j+1,k  ).ne.1) cell(:,3)=this%project(cell(:,3),i,j,k,-dt,U,V,W)
         cell(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; if (this%vmask(i+1,j+1,k+1).ne.1) cell(:,4)=this%project(cell(:,4),i,j,k,-dt,U,V,W)
         cell(:,5)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; if (this%vmask(i  ,j  ,k+1).ne.1) cell(:,5)=this%project(cell(:,5),i,j,k,-dt,U,V,W)
         cell(:,6)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; if (this%vmask(i  ,j  ,k  ).ne.1) cell(:,6)=this%project(cell(:,6),i,j,k,-dt,U,V,W)
         cell(:,7)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; if (this%vmask(i  ,j+1,k  ).ne.1) cell(:,7)=this%project(cell(:,7),i,j,k,-dt,U,V,W)
         cell(:,8)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; if (this%vmask(i  ,j+1,k+1).ne.1) cell(:,8)=this%project(cell(:,8),i,j,k,-dt,U,V,W)
         
         ! Correct volume of x- face
         face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=cell(:,5)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,6)
         face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=cell(:,7)
         face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,8)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
         cell(:,14)=getPt(remap_face,8)
         
         ! Correct volume of x+ face
         face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=cell(:,1)
         face(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,2)
         face(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=cell(:,3)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,4)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*U(i+1,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
         cell(:, 9)=getPt(remap_face,8)
         
         ! Correct volume of y- face
         face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=cell(:,2)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,6)
         face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=cell(:,5)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=cell(:,1)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
         cell(:,10)=getPt(remap_face,8)
         
         ! Correct volume of y+ face
         face(:,1)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=cell(:,3)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,6)=cell(:,7)
         face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,7)=cell(:,8)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,4)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*V(i,j+1,k)*this%cfg%dx(i)*this%cfg%dz(k))
         cell(:,12)=getPt(remap_face,8)
         
         ! Correct volume of z- face
         face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=cell(:,7)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=cell(:,6)
         face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=cell(:,2)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=cell(:,3)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
         cell(:,11)=getPt(remap_face,8)
         
         ! Correct volume of z+ face
         face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,5)=cell(:,8)
         face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,6)=cell(:,5)
         face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=cell(:,1)
         face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=cell(:,4)
         face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
         call construct(remap_face,face)
         if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*W(i,j,k+1)*this%cfg%dx(i)*this%cfg%dy(j))
         cell(:,13)=getPt(remap_face,8)
         
         ! Form remapped cell in IRL
         call construct(remap_cell,cell)
         
         ! Get bounding box for our remapped cell
         call getBoundingPts(remap_cell,bounding_pts(:,1),bounding_pts(:,2))
         bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
         bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
         
         ! Crudely check phase information for remapped cell and skip cells where nothing is changing
         crude_VF=this%crude_phase_test(bb_indices)
         if (crude_VF.ge.0.0_WP) cycle
         
         ! Need full geometric flux
         call getMoments(remap_cell,this%localized_separator_link(i,j,k),this%detailed_remap(i,j,k))
         
         ! Rebuild face flux from detailed face flux
         lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
         do n=1,getSize(this%detailed_remap(i,j,k))
            call getSepVMAtIndex(this%detailed_remap(i,j,k),n-1,my_SepVM)
            lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
            gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
         end do
         
         ! Compute new liquid volume fraction
         this%VF(i,j,k)=lvol/(lvol+gvol)
         
         ! Only work on higher order moments if VF is in [VFlo,VFhi]
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
         else
            ! Get old phasic barycenters
            this%Lbary(:,i,j,k)=lbar/lvol
            this%Gbary(:,i,j,k)=gbar/gvol
            ! Project then forward in time
            this%Lbary(:,i,j,k)=this%project(this%Lbary(:,i,j,k),i,j,k,dt,U,V,W)
            this%Gbary(:,i,j,k)=this%project(this%Gbary(:,i,j,k),i,j,k,dt,U,V,W)
         end if
         
      end do

      ! Synchronize VF and barycenter fields
      call this%cfg%sync(this%VF)
      call this%sync_and_clean_barycenters()
      
   end subroutine transport_remap_storage
   
   
   !> Perform flux-based transport of VF based on U/V/W and dt
   subroutine transport_flux(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,index
      real(IRL_double), dimension(3,9) :: face
      type(CapDod_type) :: flux_polyhedron
      real(WP) :: Lvolold,Gvolold
      real(WP) :: Lvolinc,Gvolinc
      real(WP) :: Lvolnew,Gvolnew
      real(WP) :: vol_now,crude_VF
      real(WP), dimension(3) :: ctr_now
      real(WP), dimension(3,2) :: bounding_pts
      integer, dimension(3,2) :: bb_indices
      
      ! Allocate
      call new(flux_polyhedron)
      
      ! Reset face fluxes to crude estimate (just needs to be valid for volume away from interface)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%band(i,j,k).lt.0) then
                  call construct(this%face_flux(1,i,j,k),[dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k),[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(2,i,j,k),[dt*V(i,j,k)*this%cfg%dz(k)*this%cfg%dx(i),[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(3,i,j,k),[dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j),[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               else
                  call construct(this%face_flux(1,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k),0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(2,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*V(i,j,k)*this%cfg%dz(k)*this%cfg%dx(i),0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(3,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j),0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! X flux
               if (minval(abs(this%band(i-1:i,j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  if (this%cons_correct) call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(1,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(1,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Y flux
               if (minval(abs(this%band(i,j-1:j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  if (this%cons_correct) call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(2,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(2,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Z flux
               if (minval(abs(this%band(i,j,k-1:k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j+1,k  ).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  if (this%cons_correct) call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(3,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(3,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
            end do
         end do
      end do
      
      ! Compute transported moments
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle
         
         ! Old liquid and gas volumes
         Lvolold=        this%VFold(i,j,k) *this%cfg%vol(i,j,k)
         Gvolold=(1.0_WP-this%VFold(i,j,k))*this%cfg%vol(i,j,k)
         
         ! Compute incoming liquid and gas volumes
         Lvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),0)+getVolumePtr(this%face_flux(1,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),0)+getVolumePtr(this%face_flux(2,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),0)+getVolumePtr(this%face_flux(3,i,j,k),0)
         Gvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),1)+getVolumePtr(this%face_flux(1,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),1)+getVolumePtr(this%face_flux(2,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),1)+getVolumePtr(this%face_flux(3,i,j,k),1)
         
         ! Compute new liquid and gas volumes
         Lvolnew=Lvolold+Lvolinc
         Gvolnew=Gvolold+Gvolinc
         
         ! Compute new liquid volume fraction
         this%VF(i,j,k)=Lvolnew/(Lvolnew+Gvolnew)
         
         ! Only work on higher order moments if VF is in [VFlo,VFhi]
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
         else
            ! Compute old phase barycenters
            this%Lbary(:,i,j,k)=(this%Lbary(:,i,j,k)*Lvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),0)+getCentroidPtr(this%face_flux(1,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),0)+getCentroidPtr(this%face_flux(2,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),0)+getCentroidPtr(this%face_flux(3,i,j,k),0))/Lvolnew
            this%Gbary(:,i,j,k)=(this%Gbary(:,i,j,k)*Gvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),1)+getCentroidPtr(this%face_flux(1,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),1)+getCentroidPtr(this%face_flux(2,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),1)+getCentroidPtr(this%face_flux(3,i,j,k),1))/Gvolnew
            ! Project forward in time
            this%Lbary(:,i,j,k)=this%project(this%Lbary(:,i,j,k),i,j,k,dt,U,V,W)
            this%Gbary(:,i,j,k)=this%project(this%Gbary(:,i,j,k),i,j,k,dt,U,V,W)
         end if
      end do
      
      ! Synchronize VF and barycenter fields
      call this%cfg%sync(this%VF)
      call this%sync_and_clean_barycenters()
      
   end subroutine transport_flux


   !> Perform flux-based transport of VF based on U/V/W and dt
   !> Include storage of detailed fluxes
   subroutine transport_flux_storage(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,index,n
      real(IRL_double), dimension(3,9) :: face
      type(CapDod_type) :: flux_polyhedron
      real(WP) :: Lvolold,Gvolold
      real(WP) :: Lvolinc,Gvolinc
      real(WP) :: Lvolnew,Gvolnew
      real(WP) :: vol_now,crude_VF
      real(WP) :: lvol,gvol
      real(WP), dimension(3) :: ctr_now,lbar,gbar
      real(WP), dimension(3,2) :: bounding_pts
      integer, dimension(3,2) :: bb_indices
      type(SepVM_type) :: my_SepVM
      
      ! Allocate
      call new(flux_polyhedron)
      
      ! Reset face fluxes to crude estimate (just needs to be valid for volume away from interface)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%band(i,j,k).lt.0) then
                  call construct(this%face_flux(1,i,j,k),[dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k),[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(2,i,j,k),[dt*V(i,j,k)*this%cfg%dz(k)*this%cfg%dx(i),[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(3,i,j,k),[dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j),[0.0_WP,0.0_WP,0.0_WP],0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               else
                  call construct(this%face_flux(1,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k),0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(2,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*V(i,j,k)*this%cfg%dz(k)*this%cfg%dx(i),0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
                  call construct(this%face_flux(3,i,j,k),[0.0_WP,[0.0_WP,0.0_WP,0.0_WP],dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j),0.0_WP,[0.0_WP,0.0_WP,0.0_WP]])
               end if
               ! Also empty out detailed fluxes
               call clear(this%detailed_face_flux(1,i,j,k))
               call clear(this%detailed_face_flux(2,i,j,k))
               call clear(this%detailed_face_flux(3,i,j,k))
            end do
         end do
      end do
      
      ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! X flux
               if (minval(abs(this%band(i-1:i,j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  if (this%cons_correct) call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%detailed_face_flux(1,i,j,k))
                     ! Rebuild face flux from detailed face flux
                     lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
                     do n=0,getSize(this%detailed_face_flux(1,i,j,k))-1
                        call getSepVMAtIndex(this%detailed_face_flux(1,i,j,k),n,my_SepVM)
                        lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
                        gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
                     end do
                     call construct(this%face_flux(1,i,j,k),[lvol,lbar,gvol,gbar])
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(1,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Y flux
               if (minval(abs(this%band(i,j-1:j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  if (this%cons_correct) call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%detailed_face_flux(2,i,j,k))
                     ! Rebuild face flux from detailed face flux
                     lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
                     do n=0,getSize(this%detailed_face_flux(2,i,j,k))-1
                        call getSepVMAtIndex(this%detailed_face_flux(2,i,j,k),n,my_SepVM)
                        lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
                        gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
                     end do
                     call construct(this%face_flux(2,i,j,k),[lvol,lbar,gvol,gbar])
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(2,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Z flux
               if (minval(abs(this%band(i,j,k-1:k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j+1,k  ).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  if (this%cons_correct) call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%detailed_face_flux(3,i,j,k))
                     ! Rebuild face flux from detailed face flux
                     lvol=0.0_WP; gvol=0.0_WP; lbar=0.0_WP; gbar=0.0_WP
                     do n=0,getSize(this%detailed_face_flux(3,i,j,k))-1
                        call getSepVMAtIndex(this%detailed_face_flux(3,i,j,k),n,my_SepVM)
                        lvol=lvol+getVolume(my_SepVM,0); lbar=lbar+getCentroid(my_SepVM,0)
                        gvol=gvol+getVolume(my_SepVM,1); gbar=gbar+getCentroid(my_SepVM,1)
                     end do
                     call construct(this%face_flux(3,i,j,k),[lvol,lbar,gvol,gbar])
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(3,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
            end do
         end do
      end do
      
      ! Compute transported moments
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle
         
         ! Old liquid and gas volumes
         Lvolold=        this%VFold(i,j,k) *this%cfg%vol(i,j,k)
         Gvolold=(1.0_WP-this%VFold(i,j,k))*this%cfg%vol(i,j,k)
         
         ! Compute incoming liquid and gas volumes
         Lvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),0)+getVolumePtr(this%face_flux(1,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),0)+getVolumePtr(this%face_flux(2,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),0)+getVolumePtr(this%face_flux(3,i,j,k),0)
         Gvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),1)+getVolumePtr(this%face_flux(1,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),1)+getVolumePtr(this%face_flux(2,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),1)+getVolumePtr(this%face_flux(3,i,j,k),1)
         
         ! Compute new liquid and gas volumes
         Lvolnew=Lvolold+Lvolinc
         Gvolnew=Gvolold+Gvolinc
         
         ! Compute new liquid volume fraction
         this%VF(i,j,k)=Lvolnew/(Lvolnew+Gvolnew)
         
         ! Only work on higher order moments if VF is in [VFlo,VFhi]
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
         else
            ! Compute old phase barycenters
            this%Lbary(:,i,j,k)=(this%Lbary(:,i,j,k)*Lvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),0)+getCentroidPtr(this%face_flux(1,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),0)+getCentroidPtr(this%face_flux(2,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),0)+getCentroidPtr(this%face_flux(3,i,j,k),0))/Lvolnew
            this%Gbary(:,i,j,k)=(this%Gbary(:,i,j,k)*Gvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),1)+getCentroidPtr(this%face_flux(1,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),1)+getCentroidPtr(this%face_flux(2,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),1)+getCentroidPtr(this%face_flux(3,i,j,k),1))/Gvolnew
            ! Project forward in time
            this%Lbary(:,i,j,k)=this%project(this%Lbary(:,i,j,k),i,j,k,dt,U,V,W)
            this%Gbary(:,i,j,k)=this%project(this%Gbary(:,i,j,k),i,j,k,dt,U,V,W)
         end if
      end do
      
      ! Synchronize VF and barycenter fields
      call this%cfg%sync(this%VF)
      call this%sync_and_clean_barycenters()
      
   end subroutine transport_flux_storage
   
   
   !> Project a single face to get flux polyhedron and get its moments
   subroutine fluxpoly_project_getmoments(this,i,j,k,dt,dir,U,V,W,a_flux_polyhedron,a_locseplink,some_face_flux_moments)
     implicit none
     class(vfs), intent(inout)    :: this
     integer,          intent(in) :: i,j,k    !< Index 
     real(WP),         intent(in) :: dt       !< Timestep size over which to advance
     character(len=1), intent(in) :: dir      !< Orientation of current face
     type(CapDod_type)         :: a_flux_polyhedron
     type(LocSepLink_type)     :: a_locseplink
     type(TagAccVM_SepVM_type) :: some_face_flux_moments
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
     real(IRL_double), dimension(3,9) :: face !< Points on face being projected to form flux volume
     real(WP) :: vol_f                        !< Consistent volume according to face velocity and mesh

     ! Project points, form flux volume, set directional parameters
     select case(trim(dir))
     case('x')
        ! Construct and project left x(i) face
        face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]
        face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]
        face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]
        face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]
        face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W)
        face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W)
        face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W)
        face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W)
        face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
        face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
        ! Limit projection in the case of walls
        if (this%vmask(i  ,j  ,k+1).eq.1) face(:,5)=face(:,1)
        if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
        if (this%vmask(i  ,j+1,k  ).eq.1) face(:,7)=face(:,3)
        if (this%vmask(i  ,j+1,k+1).eq.1) face(:,8)=face(:,4)
        ! Calculate consistent volume
        vol_f = dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
     case ('y')
        ! Construct and project bottom y(j) face
        face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]
        face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]
        face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]
        face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]
        face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W)
        face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W)
        face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W)
        face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W)
        face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
        face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
        ! Limit projection in the case of walls
        if (this%vmask(i+1,j  ,k  ).eq.1) face(:,5)=face(:,1)
        if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
        if (this%vmask(i  ,j  ,k+1).eq.1) face(:,7)=face(:,3)
        if (this%vmask(i+1,j  ,k+1).eq.1) face(:,8)=face(:,4)
        ! Calculate consistent volume
        vol_f = dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k)
     case('z')
        ! Construct and project bottom z(k) face
        face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]
        face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]
        face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]
        face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]
        face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W)
        face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W)
        face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W)
        face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W)
        face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
        face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
        ! Limit projection in the case of walls
        if (this%vmask(i  ,j+1,k  ).eq.1) face(:,5)=face(:,1)
        if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
        if (this%vmask(i+1,j  ,k  ).eq.1) face(:,7)=face(:,3)
        if (this%vmask(i+1,j+1,k  ).eq.1) face(:,8)=face(:,4)
        ! Calculate consistent volume
        vol_f = dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j)
     end select

     ! Form flux polyhedron
     call construct(a_flux_polyhedron,face)
     ! Add volume-consistent correction
     call adjustCapToMatchVolume(a_flux_polyhedron,vol_f)
     ! Get all volumetric fluxes
     call getNormMoments(a_flux_polyhedron,a_locseplink,some_face_flux_moments)

   end subroutine fluxpoly_project_getmoments
   
   
   !> From the moments object in a single cell, get the volume and centroid out of IRL
   subroutine fluxpoly_cell_getvolcentr(this,f_moments,n,ii,jj,kk,my_Lbary,my_Gbary,my_Lvol,my_Gvol,skip_flag)
     implicit none
     class(vfs), intent(inout) :: this
     type(TagAccVM_SepVM_type) :: f_moments
     integer, intent(in) :: n
     integer  :: ii,jj,kk,localizer_id
     integer, dimension(3) :: ind
     type(SepVM_type) :: my_SepVM
     real(WP) :: my_Gvol,my_Lvol
     real(WP), dimension(3) :: my_Gbary,my_Lbary
     logical  :: skip_flag
     
     skip_flag = .false.
     ! Get unique id of current cell
     localizer_id = getTagForIndex(f_moments,n)
     ! Convert unique id to indices ii,jj,kk
     ind=this%cfg%get_ijk_from_lexico(localizer_id)
     ii = ind(1); jj = ind(2); kk = ind(3)
     
     ! If inside wall, nothing should be added to flux
     if (this%mask(ii,jj,kk).eq.1) then
        skip_flag = .true.
        return
     end if
     
     ! If bringing material back from beyond outflow boundary, skip
     !if (backflow_flux_flag(ii,jj,kk)) cycle

     ! Get barycenter and volume of tets in current cell
     call getSepVMAtIndex(f_moments,n,my_SepVM)
     my_Lbary = getCentroid(my_SepVM, 0)
     my_Gbary = getCentroid(my_SepVM, 1)
     my_Lvol  = getVolume(my_SepVM, 0)
     my_Gvol  = getVolume(my_SepVM, 1)
     
   end subroutine fluxpoly_cell_getvolcentr
   
   
   !> Remove likely flotsams
   subroutine remove_flotsams(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,ierr
      real(WP) :: FSlo,FShi,myerror
      ! Do not do anything if VFflot<=0.0_WP
      if (this%flotsam_thld.le.0.0_WP) return
      ! Build lo and hi values
      FSlo=this%flotsam_thld
      FShi=1.0_WP-this%flotsam_thld
      ! Reset error monitoring
      this%flotsam_error=0.0_WP
      ! Loop inside and remove flotsams
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            oloop: do i=this%cfg%imin_,this%cfg%imax_
               ! Handle liquid flotsams
               if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).lt.FSlo) then
                  ! Check for fuller neighbors
                  do kk=k-1,k+1
                     do jj=j-1,j+1
                        do ii=i-1,i+1
                           if (i.eq.ii.and.j.eq.jj.and.k.eq.kk) cycle
                           if (this%VF(ii,jj,kk).ge.VFlo) cycle oloop
                        end do
                     end do
                  end do
                  ! None was found, we have an isolated flotsam
                  this%flotsam_error=this%flotsam_error+(this%VF(i,j,k)-0.0_WP)*this%cfg%vol(i,j,k)
                  this%VF(i,j,k)=0.0_WP
               end if
               ! Handle gas flotsams
               if (this%VF(i,j,k).le.VFhi.and.this%VF(i,j,k).gt.FShi) then
                  ! Check for fuller neighbors
                  do kk=k-1,k+1
                     do jj=j-1,j+1
                        do ii=i-1,i+1
                           if (i.eq.ii.and.j.eq.jj.and.k.eq.kk) cycle
                           if (this%VF(ii,jj,kk).le.VFhi) cycle oloop
                        end do
                     end do
                  end do
                  ! None was found, we have an isolated flotsam
                  this%flotsam_error=this%flotsam_error+(this%VF(i,j,k)-1.0_WP)*this%cfg%vol(i,j,k)
                  this%VF(i,j,k)=1.0_WP
               end if
            end do oloop
         end do
      end do
      ! Synchronize VF field
      call this%cfg%sync(this%VF)
      ! Gather error
      call MPI_ALLREDUCE(this%flotsam_error,myerror,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      this%flotsam_error=myerror
   end subroutine remove_flotsams
   
   
   !> Remove thin structures below a specified thickness
   subroutine remove_thinstruct(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: newVF,myerror
      ! Do not do anything if VFsheet<=0.0_WP
      if (this%thin_thld_min.le.0.0_WP) return
      ! Reset error monitoring
      this%thinstruct_error=0.0_WP
      ! First compute thickness based on advected surface and volume moments
      call this%get_thickness()
      ! Remove thin structures below cut-off
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%thickness(i,j,k).gt.0.0_WP.and.this%thickness(i,j,k).lt.this%thin_thld_min*this%cfg%meshsize(i,j,k)) then
                  newVF=real(nint(this%VF(i,j,k)),WP)
                  this%thinstruct_error=this%thinstruct_error+(this%VF(i,j,k)-newVF)*this%cfg%vol(i,j,k)
                  this%VF(i,j,k)=newVF
               end if
            end do
         end do
      end do
      ! Synchronize VF field
      call this%cfg%sync(this%VF)
      ! Gather error
      call MPI_ALLREDUCE(this%thinstruct_error,myerror,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      this%thinstruct_error=myerror
   end subroutine remove_thinstruct
   
   
   !> Measure local thickness of multiphasic structure
   subroutine get_thickness(this)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: tmp
      integer :: i,j,k,ii,jj,kk
      real(WP) :: lvol,gvol,area
      ! Reset thickness
      this%thickness=0.0_WP
      ! First compute thickness based on current surface and volume moments (SD and VF)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%mask(i,j,k).ne.0) cycle
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Extract thickness estimate from local phasic volumes and surface area
               lvol=0.0_WP; gvol=0.0_WP; area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  lvol=lvol+(       this%VF(ii,jj,kk))*this%cfg%vol(ii,jj,kk)
                  gvol=gvol+(1.0_WP-this%VF(ii,jj,kk))*this%cfg%vol(ii,jj,kk)
                  area=area+        this%SD(ii,jj,kk) *this%cfg%vol(ii,jj,kk)
               end do; end do; end do
               if (area.gt.0.0_WP) this%thickness(i,j,k)=2.0_WP*min(lvol,gvol)/area
            end do
         end do
      end do
      call this%cfg%sync(this%thickness)
      ! Filter thickness
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%mask(i,j,k).ne.0) cycle
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Surface-average thickness
               area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  area      =area      +this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
                  tmp(i,j,k)=tmp(i,j,k)+this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%thickness(ii,jj,kk)
               end do; end do; end do
               if (area.gt.0.0_WP) tmp(i,j,k)=tmp(i,j,k)/area
            end do
         end do
      end do
      call this%cfg%sync(tmp)
      this%thickness=tmp
      deallocate(tmp)
   end subroutine get_thickness
   
   
   !> Clean up after VF change removal
   subroutine clean_irl_and_band(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,n
      ! Loop everywhere and remove leftover IRL objects
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%VF(i,j,k).lt.VFlo) then
                  ! Pure gas moments
                  this%VF(i,j,k)=0.0_WP
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  ! Provide default interface
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  ! Zero out the polygons
                  do n=1,max_interface_planes
                     call zeroPolygon(this%interface_polygon(n,i,j,k))
                  end do
               else if (this%VF(i,j,k).gt.VFhi) then
                  ! Pure liquid moments
                  this%VF(i,j,k)=1.0_WP
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  ! Provide default interface
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  ! Zero out the polygons
                  do n=1,max_interface_planes
                     call zeroPolygon(this%interface_polygon(n,i,j,k))
                  end do
               end if
            end do
         end do
      end do
      ! Update the band
      call this%update_band()
   end subroutine clean_irl_and_band
   
   
   !> Synchronize and clean up barycenter fields
   subroutine sync_and_clean_barycenters(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      ! Clean up barycenters everywhere - SD too...
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%SD(i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      ! Synchronize barycenters
      call this%cfg%sync(this%Lbary)
      call this%cfg%sync(this%Gbary)
      ! Fix barycenter synchronization across periodic boundaries
      if (this%cfg%xper.and.this%cfg%iproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino,this%cfg%imin-1
                  this%Lbary(1,i,j,k)=this%Lbary(1,i,j,k)-this%cfg%xL
                  this%Gbary(1,i,j,k)=this%Gbary(1,i,j,k)-this%cfg%xL
               end do
            end do
         end do
      end if
      if (this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imax+1,this%cfg%imaxo
                  this%Lbary(1,i,j,k)=this%Lbary(1,i,j,k)+this%cfg%xL
                  this%Gbary(1,i,j,k)=this%Gbary(1,i,j,k)+this%cfg%xL
               end do
            end do
         end do
      end if
      if (this%cfg%yper.and.this%cfg%jproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino,this%cfg%jmin-1
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%Lbary(2,i,j,k)=this%Lbary(2,i,j,k)-this%cfg%yL
                  this%Gbary(2,i,j,k)=this%Gbary(2,i,j,k)-this%cfg%yL
               end do
            end do
         end do
      end if
      if (this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmax+1,this%cfg%jmaxo
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%Lbary(2,i,j,k)=this%Lbary(2,i,j,k)+this%cfg%yL
                  this%Gbary(2,i,j,k)=this%Gbary(2,i,j,k)+this%cfg%yL
               end do
            end do
         end do
      end if
      if (this%cfg%zper.and.this%cfg%kproc.eq.1) then
         do k=this%cfg%kmino,this%cfg%kmin-1
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%Lbary(3,i,j,k)=this%Lbary(3,i,j,k)-this%cfg%zL
                  this%Gbary(3,i,j,k)=this%Gbary(3,i,j,k)-this%cfg%zL
               end do
            end do
         end do
      end if
      if (this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         do k=this%cfg%kmax+1,this%cfg%kmaxo
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%Lbary(3,i,j,k)=this%Lbary(3,i,j,k)+this%cfg%zL
                  this%Gbary(3,i,j,k)=this%Gbary(3,i,j,k)+this%cfg%zL
               end do
            end do
         end do
      end if
      ! Handle 2D barycenters
      if (this%cfg%nx.eq.1) then
         do i=this%cfg%imino_,this%cfg%imaxo_
            this%Lbary(1,i,:,:)=this%cfg%xm(i)
            this%Gbary(1,i,:,:)=this%cfg%xm(i)
         end do
      end if
      if (this%cfg%ny.eq.1) then
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            this%Lbary(2,:,j,:)=this%cfg%ym(j)
            this%Gbary(2,:,j,:)=this%cfg%ym(j)
         end do
      end if
      if (this%cfg%nz.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            this%Lbary(3,:,:,k)=this%cfg%zm(k)
            this%Gbary(3,:,:,k)=this%cfg%zm(k)
         end do
      end if
   end subroutine sync_and_clean_barycenters
   
   
   !> Lagrangian advection of the IRL surface using U,V,W and dt
   subroutine advect_interface(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(in) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ii,n,t,vt
      integer :: list_size,localizer_id
      type(TagAccListVM_VMAN_type) :: accumulated_moments_from_tri
      type(ListVM_VMAN_type) :: moments_list_from_tri
      type(DivPoly_type) :: divided_polygon
      type(Tri_type) :: triangle
      real(IRL_double), dimension(1:4) :: plane_data
      integer, dimension(3) :: ind
      real(IRL_double), dimension(1:3,1:3) :: tri_vert
      type(VMAN_type) :: volume_moments_and_normal
      
      ! Clear moments from before
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               call clear(this%triangle_moments_storage(i,j,k))
            end do
         end do
      end do
      
      ! Allocate IRL data
      call new(accumulated_moments_from_tri)
      call new(moments_list_from_tri)
      call new(divided_polygon)
      call new(triangle)
      call new(volume_moments_and_normal)
      
      ! Loop over domain to forward transport interface
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               
               ! Skip if no interface
               if (this%VFold(i,j,k).lt.VFlo.or.this%VFold(i,j,k).gt.VFhi) cycle
               
               ! Construct triangulation of each interface plane
               do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  
                  ! Skip planes outside of the cell
                  if (getNumberOfVertices(this%interface_polygon(n,i,j,k)).eq.0) cycle
                  
                  ! Get DividedPolygon from the plane
                  call constructFromPolygon(divided_polygon,this%interface_polygon(n,i,j,k))
                  
                  ! Check if point ordering correct, flip if not
                  plane_data=getPlane(this%liquid_gas_interface(i,j,k),n-1)
                  if (abs(1.0_WP-dot_product(calculateNormal(divided_polygon),plane_data(1:3))).gt.1.0_WP) call reversePtOrdering(divided_polygon)
                  
                  ! Loop over triangles from DividedPolygon
                  do t=1,getNumberOfSimplicesInDecomposition(divided_polygon)
                     ! Get the triangle
                     call getSimplexFromDecomposition(divided_polygon,t-1,triangle)
                     ! Forward project triangle vertices
                     tri_vert=getVertices(triangle)
                     do vt=1,3
                        tri_vert(:,vt)=this%project(tri_vert(:,vt),i,j,k,dt,U,V,W)
                     end do
                     call construct(triangle,tri_vert)
                     call calculateAndSetPlaneOfExistence(triangle)
                     ! Cut it by the mesh
                     call getMoments(triangle,this%localizer_link(i,j,k),accumulated_moments_from_tri)
                     ! Loop through each cell and append to triangle_moments_storage in each cell
                     list_size=getSize(accumulated_moments_from_tri)
                     do ii=1,list_size
                        localizer_id=getTagForIndex(accumulated_moments_from_tri,ii-1)
                        ind=this%cfg%get_ijk_from_lexico(localizer_id)
                        call getListAtIndex(accumulated_moments_from_tri,ii-1,moments_list_from_tri)
                        call append(this%triangle_moments_storage(ind(1),ind(2),ind(3)),moments_list_from_tri)
                     end do
                  end do
                  
               end do
               
            end do
         end do
      end do

      ! Recompute surface density from advected interface
      this%SD=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do ii=0,getSize(this%triangle_moments_storage(i,j,k))-1
                  call getMoments(this%triangle_moments_storage(i,j,k),ii,volume_moments_and_normal)
                  this%SD(i,j,k)=this%SD(i,j,k)+getVolume(volume_moments_and_normal)
               end do
               this%SD(i,j,k)=this%SD(i,j,k)/this%cfg%vol(i,j,k)
            end do
         end do
      end do
      call this%cfg%sync(this%SD)
      
   end subroutine advect_interface
   
   
   !> Band update from VF dataset
   subroutine update_band(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,dir,n,index
      integer, dimension(3) :: ind
      integer :: ibmin_,ibmax_,jbmin_,jbmax_,kbmin_,kbmax_
      integer, dimension(0:nband) :: band_size
      
      ! Loop extents
      ibmin_=this%cfg%imin_; if (this%cfg%iproc.eq.1           .and..not.this%cfg%xper) ibmin_=this%cfg%imin-1
      ibmax_=this%cfg%imax_; if (this%cfg%iproc.eq.this%cfg%npx.and..not.this%cfg%xper) ibmax_=this%cfg%imax+1
      jbmin_=this%cfg%jmin_; if (this%cfg%jproc.eq.1           .and..not.this%cfg%yper) jbmin_=this%cfg%jmin-1
      jbmax_=this%cfg%jmax_; if (this%cfg%jproc.eq.this%cfg%npy.and..not.this%cfg%yper) jbmax_=this%cfg%jmax+1
      kbmin_=this%cfg%kmin_; if (this%cfg%kproc.eq.1           .and..not.this%cfg%zper) kbmin_=this%cfg%kmin-1
      kbmax_=this%cfg%kmax_; if (this%cfg%kproc.eq.this%cfg%npz.and..not.this%cfg%zper) kbmax_=this%cfg%kmax+1
      
      ! Reset band
      this%band=(nband+1)*int(sign(1.0_WP,this%VF-0.5_WP))
      
      ! First sweep to identify cells with interface
      do k=kbmin_,kbmax_
         do j=jbmin_,jbmax_
            do i=ibmin_,ibmax_
               ! Skip *real* wall cells
               if (this%mask(i,j,k).eq.1) cycle
               ! Identify cells with interface
               if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) this%band(i,j,k)=0
               ! Check if cell-face is an interface
               do dir=1,3
                  do n=-1,+1,2
                     ind=[i,j,k]; ind(dir)=ind(dir)+n
                     if (this%mask(ind(1),ind(2),ind(3)).ne.1) then
                        if (this%VF(i,j,k).lt.VFlo.and.this%VF(ind(1),ind(2),ind(3)).gt.VFhi.or.&
                        &   this%VF(i,j,k).gt.VFhi.and.this%VF(ind(1),ind(2),ind(3)).lt.VFlo) this%band(i,j,k)=0
                     end if
                  end do
               end do
            end do
         end do
      end do
      call this%cfg%sync(this%band)
      
      ! Sweep to identify the bands up to nband
      do n=1,nband
         ! For each band
         do k=kbmin_,kbmax_
            do j=jbmin_,jbmax_
               do i=ibmin_,ibmax_
                  ! Skip wall cells
                  if (this%mask(i,j,k).eq.1) cycle
                  ! Work on one band at a time
                  if (abs(this%band(i,j,k)).gt.n) then
                     ! Loop over 26 neighbors
                     do kk=k-1,k+1
                        do jj=j-1,j+1
                           do ii=i-1,i+1
                              ! Skip wall cells
                              if (this%mask(ii,jj,kk).eq.1) cycle
                              ! Extend the band
                              if (abs(this%band(ii,jj,kk)).eq.n-1) this%band(i,j,k)=int(sign(real(n,WP),this%VF(i,j,k)-0.5_WP))
                           end do
                        end do
                     end do
                  end if
               end do
            end do
         end do
         call this%cfg%sync(this%band)
      end do
      
      ! Count the number of cells in each band value
      band_size=0
      do k=kbmin_,kbmax_
         do j=jbmin_,jbmax_
            do i=ibmin_,ibmax_
               if (abs(this%band(i,j,k)).le.nband) band_size(abs(this%band(i,j,k)))=band_size(abs(this%band(i,j,k)))+1
            end do
         end do
      end do
      
      ! Rebuild the unstructured mapping
      if (allocated(this%band_map)) deallocate(this%band_map); allocate(this%band_map(3,sum(band_size)))
      this%band_count=0
      do k=kbmin_,kbmax_
         do j=jbmin_,jbmax_
            do i=ibmin_,ibmax_
               if (abs(this%band(i,j,k)).le.nband) then
                  this%band_count(abs(this%band(i,j,k)))=this%band_count(abs(this%band(i,j,k)))+1
                  index=sum(band_size(0:abs(this%band(i,j,k))-1))+this%band_count(abs(this%band(i,j,k)))
                  this%band_map(:,index)=[i,j,k]
               end if
            end do
         end do
      end do
      
   end subroutine update_band
   
   
   !> Reconstruct an IRL interface from the VF field distribution
   subroutine build_interface(this)
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
      ! Reconstruct interface - will need to support various methods
      select case (this%reconstruction_method)
      case (elvira); call this%build_elvira()
      case (lvira) ; call this%build_lvira()
      case (mof)   ; call this%build_mof()
      case (wmof)  ; call this%build_wmof()
      case (r2p)   ; call this%build_r2p()
      case (swartz)
         call this%build_lvira()
         call this%smooth_interface()
      case (youngs); call this%build_youngs()
      !case (lvlset); call this%build_lvlset()
      case (plicnet); call this%build_plicnet()
      case (r2pnet); call this%build_r2pnet()
      case default; call die('[vfs build interface] Unknown interface reconstruction scheme')
      end select
   end subroutine build_interface
   

   !> Compute interface sensors
   subroutine sense_interface(this)
      implicit none
      class(vfs), intent(inout) :: this
      ! Update local thickness
      call this%get_thickness()
      ! Identify thin regions
      call this%detect_thin_regions()
      ! Identify edge regions
      call this%detect_edge_regions()
   end subroutine sense_interface
   
   
   !> Detect thin regions of the interface
   subroutine detect_thin_regions(this)
      use mathtools, only: normalize
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,dim,dir,ni
      integer , dimension(3) :: pos
      real(WP), dimension(3) :: n1,n2,c1,c2
      real(WP), dimension(:,:,:), allocatable :: mysensor
      real(WP) :: a1,a2
      ! Default value is 0
      this%thin_sensor=0.0_WP
      ! First pass to handle 2-plane cells
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip wall/bcond cells
               if (this%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Detect thin regions from local polygon data
               n1=calculateNormal(this%interface_polygon(1,i,j,k))
               if (getNumberOfVertices(this%interface_polygon(2,i,j,k)).gt.0) then
                  n2=calculateNormal(this%interface_polygon(2,i,j,k))
                  ! Check normal orientation to identify thin regions
                  if (dot_product(n1,n2).lt.this%thin_thld_dotprod) then
                     ! Check if liquid or gas
                     c1=calculateCentroid(this%interface_polygon(1,i,j,k))
                     c2=calculateCentroid(this%interface_polygon(2,i,j,k))
                     if (dot_product(c2-c1,n2).gt.0.0_WP) then
                        ! Thin liquid region
                        this%thin_sensor(i,j,k)=1.0_WP
                     else
                        ! Thin gas region
                        this%thin_sensor(i,j,k)=2.0_WP
                     end if
                  end if
               end if
            end do
         end do
      end do
      ! Second pass to extend sensor
      allocate(mysensor(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); mysensor=this%thin_sensor
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Sensor is still zero, check direct neighbors
               if (this%thin_sensor(i,j,k).eq.0.0_WP) then
                  do dim=1,3
                     do dir=-1,1,2
                        pos=0; pos(dim)=dir; ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                        if (this%thin_sensor(ii,jj,kk).gt.0.0_WP) mysensor(i,j,k)=this%thin_sensor(ii,jj,kk)
                     end do
                  end do
               end if
            end do
         end do
      end do
      this%thin_sensor=mysensor
      call this%cfg%sync(this%thin_sensor)
      deallocate(mysensor)
      ! Final pass to check single-plane cells
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Sensor is still zero and 1-plane cell, check direct neighbors
               if (this%thin_sensor(i,j,k).eq.0.0_WP.and.getNumberOfPlanes(this%liquid_gas_interface(i,j,k)).eq.1) then
                  n1=calculateNormal(this%interface_polygon(1,i,j,k))
                  do dim=1,3
                     do dir=-1,1,2
                        pos=0; pos(dim)=dir; ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                        ! Skip full cells
                        if (this%VF(ii,jj,kk).lt.VFlo.or.this%VF(ii,jj,kk).gt.VFhi) cycle
                        ! Check normal orientation to identify thin regions
                        ! If neighbor has two planes, then surface-average its normals and centroids
                        if (getNumberOfVertices(this%interface_polygon(2,ii,jj,kk)).gt.0) then
                           a1=calculateVolume(this%interface_polygon(1,ii,jj,kk))/this%cfg%meshsize(ii,jj,kk)
                           a2=calculateVolume(this%interface_polygon(2,ii,jj,kk))/this%cfg%meshsize(ii,jj,kk)
                           n2=normalize(a1*calculateNormal(this%interface_polygon(1,ii,jj,kk))&
                           &           +a2*calculateNormal(this%interface_polygon(2,ii,jj,kk)))
                           if (dot_product(n1,n2).lt.this%thin_thld_dotprod) then
                              c1=calculateCentroid(this%interface_polygon(1,i ,j ,k ))
                              c2=(a1*calculateCentroid(this%interface_polygon(1,ii,jj,kk))&
                              &  +a2*calculateCentroid(this%interface_polygon(2,ii,jj,kk)))/(a1+a2)
                              ! Check if liquid or gas
                              if (dot_product(c2-c1,n2).gt.0.0_WP.and.dot_product(c1-c2,n1).gt.0.0_WP) then
                                 this%thin_sensor(i,j,k)=1.0_WP
                              else if (dot_product(c2-c1,n2).lt.0.0_WP.and.dot_product(c1-c2,n1).lt.0.0_WP) then
                                 this%thin_sensor(i,j,k)=2.0_WP
                              else
                                 this%thin_sensor(i,j,k)=this%thin_sensor(ii,jj,kk) ! what if it hasn't been assigned a thin sensor value yet?                          
                              end if
                           end if
                        else
                           n2=calculateNormal(this%interface_polygon(1,ii,jj,kk))
                           if (dot_product(n1,n2).lt.this%thin_thld_dotprod) then
                              ! Check if liquid or gas
                              c1=calculateCentroid(this%interface_polygon(1,i ,j ,k ))
                              c2=calculateCentroid(this%interface_polygon(1,ii,jj,kk))
                              ! Check if liquid or gas
                              if (dot_product(c2-c1,n2).gt.0.0_WP.and.dot_product(c1-c2,n1).gt.0.0_WP) then
                                 this%thin_sensor(i,j,k)=1.0_WP
                              else if (dot_product(c2-c1,n2).lt.0.0_WP.and.dot_product(c1-c2,n1).lt.0.0_WP) then
                                 this%thin_sensor(i,j,k)=2.0_WP
                              else
                                 this%thin_sensor(i,j,k)=this%thin_sensor(ii,jj,kk) ! what if it hasn't been assigned a thin sensor value yet?                          
                              end if
                           end if
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end do
      call this%cfg%sync(this%thin_sensor)
      ! Finally, ensure thin regions have small enough thickness
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%thickness(i,j,k).gt.this%thin_thld_max*this%cfg%meshsize(i,j,k)) this%thin_sensor(i,j,k)=0.0_WP
            end do
         end do
      end do
   end subroutine detect_thin_regions
   
   
   !> Detect edge-like regions of the interface
   subroutine detect_edge_regions(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,ni
      real(WP) :: fvol,myvol,volume_sensor,surface_sensor
      real(WP), dimension(3) :: fbary,mybary
      real(WP), dimension(:,:,:)  , allocatable :: s_tmp
      real(WP), dimension(:,:,:,:), allocatable :: v_tmp
      real(WP) :: surface_area
      ! Default value is 0
      this%edge_sensor=0.0_WP
      this%edge_normal=0.0_WP
      ! Traverse domain and compute sensors
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells
               if (this%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Compute filtered volume barycenter sensor
               fvol=0.0_WP; fbary=0.0_WP
               do kk=k-2,k+2
                  do jj=j-2,j+2
                     do ii=i-2,i+2
                        myvol=this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        fvol =fvol +myvol
                        fbary=fbary+myvol*this%Lbary(:,ii,jj,kk)
                     end do
                  end do
               end do
               fbary=fbary/fvol
               volume_sensor=norm2(fbary-this%Lbary(:,i,j,k))/this%cfg%meshsize(i,j,k)
               ! Compute filtered surface barycenter sensor
               fvol=0.0_WP; mybary=0.0_WP
               do ni=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%interface_polygon(ni,i,j,k)).ne.0) then
                     myvol =abs(calculateVolume(this%interface_polygon(ni,i,j,k)))
                     fvol  =fvol  +myvol
                     mybary=mybary+myvol*calculateCentroid(this%interface_polygon(ni,i,j,k))
                  end if
               end do
               mybary=mybary/fvol
               fvol=0.0_WP; fbary=0.0_WP
               do kk=k-2,k+2
                  do jj=j-2,j+2
                     do ii=i-2,i+2
                        do ni=1,getNumberOfPlanes(this%liquid_gas_interface(ii,jj,kk))
                           if (getNumberOfVertices(this%interface_polygon(ni,ii,jj,kk)).ne.0) then
                              myvol=abs(calculateVolume(this%interface_polygon(ni,ii,jj,kk)))
                              fvol =fvol +myvol
                              fbary=fbary+myvol*calculateCentroid(this%interface_polygon(ni,ii,jj,kk))
                           end if
                        end do
                     end do
                  end do
               end do
               fbary=fbary/fvol
               surface_sensor=norm2(fbary-mybary)/this%cfg%meshsize(i,j,k)
               ! Aggregate into an edge sensor
               this%edge_sensor(i,j,k)=volume_sensor*surface_sensor
               ! Clip based on thickness
               if (this%thickness(i,j,k).gt.this%thin_thld_max*this%cfg%meshsize(i,j,k)) this%edge_sensor(i,j,k)=0.0_WP
               ! Finally, store edge orientation
               if (this%edge_sensor(i,j,k).gt.0.0_WP) then
                  this%edge_normal(:,i,j,k)=(fbary-mybary)/(norm2(fbary-mybary)+epsilon(1.0_WP))
               end if
            end do
         end do
      end do
      ! Communicate
      call this%cfg%sync(this%edge_sensor)
      call this%cfg%sync(this%edge_normal)
      ! Apply an extra step of surface smoothing to our edge info
      allocate(s_tmp(    this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); s_tmp=0.0_WP
      allocate(v_tmp(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); v_tmp=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%mask(i,j,k).ne.0) cycle
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Surface-averaged normal magnitude
               surface_area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  surface_area  =  surface_area+this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
                  s_tmp  (i,j,k)=s_tmp  (i,j,k)+this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%edge_sensor  (ii,jj,kk)
                  v_tmp(:,i,j,k)=v_tmp(:,i,j,k)+this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%edge_normal(:,ii,jj,kk)
               end do; end do; end do
               if (surface_area.gt.0.0_WP) then
                  s_tmp  (i,j,k)=s_tmp  (i,j,k)/surface_area
                  v_tmp(:,i,j,k)=v_tmp(:,i,j,k)/surface_area
                  v_tmp(:,i,j,k)=v_tmp(:,i,j,k)/(norm2(v_tmp(:,i,j,k))+epsilon(1.0_WP))
               end if
            end do
         end do
      end do
      call this%cfg%sync(s_tmp); this%edge_sensor=s_tmp; deallocate(s_tmp)
      call this%cfg%sync(v_tmp); this%edge_normal=v_tmp; deallocate(v_tmp)
   end subroutine detect_edge_regions
   
   
   !> ELVIRA reconstruction of a planar interface in mixed cells
   subroutine build_elvira(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer(IRL_SignedIndex_t) :: i,j,k
      integer :: ind,ii,jj,kk
      type(ELVIRANeigh_type) :: neighborhood
      type(RectCub_type), dimension(0:26) :: neighborhood_cells
      real(IRL_double), dimension(0:26) :: liquid_volume_fraction
      
      ! Give ourselves a ELVIRA neighborhood of 27 cells
      call new(neighborhood)
      do i=0,26
         call new(neighborhood_cells(i))
      end do
      call setSize(neighborhood,27)
      ind=0
      do k=-1,+1
         do j=-1,+1
            do i=-1,+1
               call setMember(neighborhood,neighborhood_cells(ind),liquid_volume_fraction(ind),i,j,k)
               ind=ind+1
            end do
         end do
      end do
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if

               ! Set neighborhood_cells and liquid_volume_fraction to current correct values
               ind=0
               do kk=k-1,k+1
                  do jj=j-1,j+1
                     do ii=i-1,i+1
                        ! Build the cell
                        call construct_2pt(neighborhood_cells(ind),[this%cfg%x(ii),this%cfg%y(jj),this%cfg%z(kk)],[this%cfg%x(ii+1),this%cfg%y(jj+1),this%cfg%z(kk+1)])
                        ! Assign volume fraction
                        liquid_volume_fraction(ind)=this%VF(ii,jj,kk)
                        ! Increment counter
                        ind=ind+1
                     end do
                  end do
               end do
               ! Perform the reconstruction
               call reconstructELVIRA3D(neighborhood,this%liquid_gas_interface(i,j,k))
            end do
         end do
      end do
      
      ! Synchronize across boundaries
      call this%sync_interface()
      
   end subroutine build_elvira

   
   !> Smoothing of an IRL interface based on Swartz-like algorithm
   subroutine smooth_interface(this)
      use mathtools, only: cross_product,normalize,Pi,qrotate
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,ierr,ite,count
      real(WP) :: myres,res,mag
      real(WP), dimension(3) :: mynorm,mybary,bary,norm,newnorm,r
      real(WP), dimension(4) :: plane,q
      type(RectCub_type) :: cell
      !real(WP), parameter :: norm_threshold=0.85_WP ! About 30 degrees
      real(WP), parameter :: norm_threshold=0.0_WP ! 0 degrees
      real(WP), parameter :: maxres=1.0e-6_WP
      integer , parameter :: maxite=5
      
      ! Allocate cell
      call new(cell)
      
      ! Iterate until convergence criterion is met
      res=huge(1.0_WP); ite=0
      do while (res.ge.maxres.and.ite.lt.maxite)
         
         ! Create discontinuous polygon mesh from IRL interface
         call this%polygonalize_interface()
         
         ! Traverse domain and form new normal
         myres=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  ! Skip wall/bcond cells
                  if (this%mask(i,j,k).ne.0) cycle
                  ! Skip cells without interface
                  if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
                  ! Compute polygon barycenter and normal
                  mybary=calculateCentroid(this%interface_polygon(1,i,j,k))
                  mynorm=calculateNormal  (this%interface_polygon(1,i,j,k))
                  ! Loop over our neighbors and form new normal
                  newnorm=0.0_WP; count=0
                  do kk=k-1,k+1
                     do jj=j-1,j+1
                        do ii=i-1,i+1
                           ! Skip stencil center
                           if (ii.eq.i.and.jj.eq.j.and.kk.eq.k) cycle
                           ! Skip cells without polygons
                           if (getNumberOfVertices(this%interface_polygon(1,ii,jj,kk)).eq.0) cycle
                           ! Increment neighbor counter
                           count=count+1
                           ! Compute polygon barycenter and normal
                           bary=calculateCentroid(this%interface_polygon(1,ii,jj,kk))-mybary
                           norm=calculateNormal  (this%interface_polygon(1,ii,jj,kk))
                           ! Skip polygons with normal too different from ours
                           if (dot_product(mynorm,norm).lt.norm_threshold) cycle
                           ! Build a quaternion to rotate Pi/2 around r axis
                           r=cross_product(bary,mynorm)
                           q(1)=cos(0.25_WP*Pi); q(2:4)=sin(0.25_WP*Pi)*normalize(r)
                           ! Increment our normal estimate using a barycenter-based normal
                           newnorm=newnorm+qrotate(v=bary,q=q)
                        end do
                     end do
                  end do
                  ! Set minimum number of neighbors to 1
                  if (count.eq.0) cycle
                  ! Normalize new normal vector
                  newnorm=normalize(newnorm)
                  ! Adjust plane orientation (not position yet)
                  plane=getPlane(this%liquid_gas_interface(i,j,k),0)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,newnorm,plane(4))
                  call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
                  call matchVolumeFraction(cell,this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
                  ! Monitor convergence
                  myres=max(myres,1.0_WP-dot_product(mynorm,newnorm))
               end do
            end do
         end do
         
         ! Collect maximum residual and increment iteration counter
         call MPI_ALLREDUCE(myres,res,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); ite=ite+1
         if (this%cfg%amRoot) print*,'ite=',ite,'residual=',res
         
         ! Synchronize across boundaries
         call this%sync_interface()
         
      end do
      
   end subroutine smooth_interface
   
   
   !> Youngs' algorithm for reconstructing a planar interface
   subroutine build_youngs(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if

               ! Apply Youngs' method to get normal
               

            end do
         end do
      end do

      ! Synchronize across boundaries
      call this%sync_interface()

   end subroutine build_youngs
   
   
   !> LVIRA reconstruction of a planar interface in mixed cells
   subroutine build_lvira(this)
      use mathtools, only: normalize
      implicit none
      class(vfs), intent(inout) :: this
      integer(IRL_SignedIndex_t) :: i,j,k
      integer :: ind,ii,jj,kk,icenter
      type(LVIRANeigh_RectCub_type) :: neighborhood
      type(RectCub_type), dimension(0:26) :: neighborhood_cells
      real(IRL_double)  , dimension(0:26) :: liquid_volume_fraction
      real(IRL_double), dimension(3) :: initial_norm
      real(IRL_double) :: initial_dist
      
      ! Give ourselves an LVIRA neighborhood of 27 cells
      call new(neighborhood)
      do i=0,26
         call new(neighborhood_cells(i))
      end do
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if
               
               ! Set neighborhood_cells and liquid_volume_fraction to current correct values
               ind=0
               do kk=k-1,k+1
                  do jj=j-1,j+1
                     do ii=i-1,i+1
                        ! Skip true wall cells - bconds can be used here
                        if (this%mask(ii,jj,kk).eq.1) cycle
                        ! Add cell to neighborhood
                        call addMember(neighborhood,neighborhood_cells(ind),liquid_volume_fraction(ind))
                        ! Build the cell
                        call construct_2pt(neighborhood_cells(ind),[this%cfg%x(ii),this%cfg%y(jj),this%cfg%z(kk)],[this%cfg%x(ii+1),this%cfg%y(jj+1),this%cfg%z(kk+1)])
                        ! Assign volume fraction
                        liquid_volume_fraction(ind)=this%VF(ii,jj,kk)
                        ! Trap and set stencil center
                        if (ii.eq.i.and.jj.eq.j.and.kk.eq.k) then
                           icenter=ind
                           call setCenterOfStencil(neighborhood,icenter)
                        end if
                        ! Increment counter
                        ind=ind+1
                     end do
                  end do
               end do
               
               ! Formulate initial guess
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               initial_norm=normalize(this%Gbary(:,i,j,k)-this%Lbary(:,i,j,k))
               initial_dist=dot_product(initial_norm,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
               call setPlane(this%liquid_gas_interface(i,j,k),0,initial_norm,initial_dist)
               call matchVolumeFraction(neighborhood_cells(icenter),this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
               
               ! Perform the reconstruction
               call reconstructLVIRA3D(neighborhood,this%liquid_gas_interface(i,j,k))
               
               ! Clean up neighborhood
               call emptyNeighborhood(neighborhood)
               
            end do
         end do
      end do
      
      ! Synchronize across boundaries
      call this%sync_interface()
      
   end subroutine build_lvira
   
   
   !> MOF reconstruction of a planar interface in mixed cells
   subroutine build_mof(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer(IRL_SignedIndex_t) :: i,j,k
      type(RectCub_type) :: my_cell
      type(SepVM_type)   :: separated_volume_moments
      
      ! Storage for a cell and corresponding separated_volume_moments
      call new(my_cell)
      call new(separated_volume_moments)
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if
               
               ! Perform MoF reconstruction
               call construct_2pt(my_cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               call construct(separated_volume_moments,[this%VF(i,j,k)*this%cfg%vol(i,j,k),this%Lbary(:,i,j,k),(1.0_WP-this%VF(i,j,k))*this%cfg%vol(i,j,k),this%Gbary(:,i,j,k)])
               call reconstructMOF3D(my_cell,separated_volume_moments,this%liquid_gas_interface(i,j,k))
               
            end do
         end do
      end do
      
      ! Synchronize across boundaries
      call this%sync_interface()
      
   end subroutine build_mof


   !> Wide-MOF reconstruction of a planar interface in mixed cells
   subroutine build_wmof(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: ii,jj,kk
      real(WP) :: lvol,gvol
      real(WP), dimension(3) :: lbary,gbary
      integer(IRL_SignedIndex_t) :: i,j,k
      type(RectCub_type) :: my_cell
      type(SepVM_type)   :: separated_volume_moments
      
      ! Storage for a cell and corresponding separated_volume_moments
      call new(my_cell)
      call new(separated_volume_moments)
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if
               
               ! Compute moments on a 3x3x3 stencil
               lvol=0.0_WP; gvol=0.0_WP; lbary=0.0_WP; gbary=0.0_WP
               do kk=k-1,k+1
                  do jj=j-1,j+1
                     do ii=i-1,i+1
                        lvol =lvol +this%cfg%vol(ii,jj,kk)*        this%VF(ii,jj,kk)
                        gvol =gvol +this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        lbary=lbary+this%cfg%vol(ii,jj,kk)*        this%VF(ii,jj,kk) *this%Lbary(:,ii,jj,kk)
                        gbary=gbary+this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))*this%Gbary(:,ii,jj,kk)
                     end do
                  end do
               end do
               lbary=lbary/lvol; gbary=gbary/gvol
               
               ! Perform MoF reconstruction with these wider moments
               call construct_2pt(my_cell,[this%cfg%x(i-1),this%cfg%y(j-1),this%cfg%z(k-1)],[this%cfg%x(i+2),this%cfg%y(j+2),this%cfg%z(k+2)])
               call construct(separated_volume_moments,[lvol,lbary,gvol,gbary])
               call reconstructMOF3D(my_cell,separated_volume_moments,this%liquid_gas_interface(i,j,k))
               
               ! Reset distance to match VOF in central cell
               call construct_2pt(my_cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               call matchVolumeFraction(my_cell,this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
               
            end do
         end do
      end do
      
      ! Synchronize across boundaries
      call this%sync_interface()
      
   end subroutine build_wmof

   
   !> R2P reconstruction of a planar interface in mixed cells
   subroutine build_r2p(this)
      use mathtools, only: normalize
      implicit none
      class(vfs), intent(inout) :: this
      integer(IRL_SignedIndex_t) :: i,j,k
      integer :: ind,ii,jj,kk,icenter
      type(LVIRANeigh_RectCub_type) :: nh_lvr
      type(R2PNeigh_RectCub_type)   :: nh_r2p
      type(RectCub_type), dimension(0:26) :: neighborhood_cells
      real(IRL_double)  , dimension(0:26) :: liquid_volume_fraction
      type(SepVM_type)  , dimension(0:26) :: separated_volume_moments
      type(VMAN_type) :: volume_moments_and_normal
      
      !type(R2PWeighting_type) :: r2p_weight
      real(WP) :: surface_area,area!,l2g_weight
      real(WP), dimension(3) :: surface_norm
      real(WP), dimension(:,:,:), allocatable :: surf_norm_mag,tmp
      
      real(IRL_double), dimension(3) :: initial_norm
      real(IRL_double) :: initial_dist
      logical :: is_wall
      
      ! Get storage for volume moments and normal
      call new(volume_moments_and_normal)

      ! Get r2p object for optimization weights
      !call new(r2p_weight)
      
      ! Give ourselves an R2P and an LVIRA neighborhood of 27 cells along with separated volume moments
      call new(nh_r2p)
      call new(nh_lvr)
      do i=0,26
         call new(neighborhood_cells(i))
         call new(separated_volume_moments(i))
      end do
      
      ! Compute magnitude of the surface-averaged normal vector
      allocate(surf_norm_mag(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); surf_norm_mag=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%mask(i,j,k).ne.0) cycle
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Extract average normal magnitude from neighborhood surface moments
               surface_area=0.0_WP; surface_norm=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  do ind=0,getSize(this%triangle_moments_storage(ii,jj,kk))-1
                     call getMoments(this%triangle_moments_storage(ii,jj,kk),ind,volume_moments_and_normal)
                     surface_area=surface_area+getVolume(volume_moments_and_normal)
                     surface_norm=surface_norm+getNormal(volume_moments_and_normal)
                  end do
               end do; end do; end do
               if (surface_area.gt.0.0_WP) surf_norm_mag(i,j,k)=norm2(surface_norm/surface_area)
            end do
         end do
      end do
      call this%cfg%sync(surf_norm_mag)
      
      ! Apply an extra step of surface smoothing to our normal magnitude
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond/full cells
               if (this%mask(i,j,k).ne.0) cycle
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
               ! Surface-averaged normal magnitude
               surface_area=0.0_WP
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  surface_area=surface_area+this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
                  tmp(i,j,k)  =tmp(i,j,k)  +this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*surf_norm_mag(ii,jj,kk)
               end do; end do; end do
               if (surface_area.gt.0.0_WP) tmp(i,j,k)=tmp(i,j,k)/surface_area
            end do
         end do
      end do
      call this%cfg%sync(tmp); surf_norm_mag=tmp; deallocate(tmp)
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if
               
               ! If a wall is in our neighborhood, apply LVIRA
               is_wall=.false.
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  if (this%mask(ii,jj,kk).eq.1) is_wall=.true.
               end do; end do; end do
               if (is_wall) then
                  ! Set neighborhood_cells and liquid_volume_fraction to current correct values
                  ind=0; call emptyNeighborhood(nh_lvr)
                  do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                     ! Skip true wall cells - bconds can be used here
                     if (this%mask(ii,jj,kk).eq.1) cycle
                     ! Add cell to neighborhood
                     call addMember(nh_lvr,neighborhood_cells(ind),liquid_volume_fraction(ind))
                     ! Build the cell
                     call construct_2pt(neighborhood_cells(ind),[this%cfg%x(ii),this%cfg%y(jj),this%cfg%z(kk)],[this%cfg%x(ii+1),this%cfg%y(jj+1),this%cfg%z(kk+1)])
                     ! Assign volume fraction
                     liquid_volume_fraction(ind)=this%VF(ii,jj,kk)
                     ! Trap and set stencil center
                     if (ii.eq.i.and.jj.eq.j.and.kk.eq.k) then
                        icenter=ind
                        call setCenterOfStencil(nh_lvr,icenter)
                     end if
                     ! Increment counter
                     ind=ind+1
                  end do; end do; end do
                  ! Formulate initial guess
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  initial_norm=normalize(this%Gbary(:,i,j,k)-this%Lbary(:,i,j,k))
                  initial_dist=dot_product(initial_norm,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
                  call setPlane(this%liquid_gas_interface(i,j,k),0,initial_norm,initial_dist)
                  call matchVolumeFraction(neighborhood_cells(icenter),this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
                  ! Perform the reconstruction
                  call reconstructLVIRA3D(nh_lvr,this%liquid_gas_interface(i,j,k))
                  ! Done with that cell
                  cycle
               end if
               
               ! If the neighborhood normals are sufficiently consistent, just use LVIRA
               if (surf_norm_mag(i,j,k).gt.this%twoplane_thld2) then
                  ! Build LVIRA neighborhood
                  ind=0; call emptyNeighborhood(nh_lvr)
                  do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                     call addMember(nh_lvr,neighborhood_cells(ind),liquid_volume_fraction(ind))
                     call construct_2pt(neighborhood_cells(ind),[this%cfg%x(ii),this%cfg%y(jj),this%cfg%z(kk)],[this%cfg%x(ii+1),this%cfg%y(jj+1),this%cfg%z(kk+1)])
                     liquid_volume_fraction(ind)=this%VF(ii,jj,kk)
                     if (ii.eq.i.and.jj.eq.j.and.kk.eq.k) then
                        icenter=ind
                        call setCenterOfStencil(nh_lvr,icenter)
                     end if
                     ind=ind+1
                  end do; end do; end do
                  ! Formulate initial guess
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  initial_norm=normalize(this%Gbary(:,i,j,k)-this%Lbary(:,i,j,k))
                  initial_dist=dot_product(initial_norm,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
                  call setPlane(this%liquid_gas_interface(i,j,k),0,initial_norm,initial_dist)
                  call matchVolumeFraction(neighborhood_cells(icenter),this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
                  ! Perform the reconstruction
                  call reconstructLVIRA3D(nh_lvr,this%liquid_gas_interface(i,j,k))
                  ! Done with that cell
                  cycle
               end if
               
               ! Prepare R2P data
               ind=0; call emptyNeighborhood(nh_r2p)
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  call addMember(nh_r2p,neighborhood_cells(ind),separated_volume_moments(ind))
                  call construct_2pt(neighborhood_cells(ind),[this%cfg%x(ii),this%cfg%y(jj),this%cfg%z(kk)],[this%cfg%x(ii+1),this%cfg%y(jj+1),this%cfg%z(kk+1)])
                  call construct(separated_volume_moments(ind),[this%VF(ii,jj,kk)*this%cfg%vol(ii,jj,kk),this%Lbary(:,ii,jj,kk),(1.0_WP-this%VF(ii,jj,kk))*this%cfg%vol(ii,jj,kk),this%Gbary(:,ii,jj,kk)])
                  if (ii.eq.i.and.jj.eq.j.and.kk.eq.k) then
                     icenter=ind
                     call setCenterOfStencil(nh_r2p,icenter)
                  end if
                  ind=ind+1
               end do; end do; end do
               
               ! Generate initial guess for R2P based on availability of in-cell surface data
               surface_area=0.0_WP
               do ind=0,getSize(this%triangle_moments_storage(i,j,k))-1
                  call getMoments(this%triangle_moments_storage(i,j,k),ind,volume_moments_and_normal)
                  surface_area=surface_area+getVolume(volume_moments_and_normal)
               end do
               if (surface_area.gt.surface_epsilon_factor*this%cfg%meshsize(i,j,k)**2) then
                  ! Local normals are available, reconstruction from surface data
                  call reconstructAdvectedNormals(this%triangle_moments_storage(i,j,k),nh_r2p,this%twoplane_thld1,this%liquid_gas_interface(i,j,k))
                  if (getNumberOfPlanes(this%liquid_gas_interface(i,j,k)).eq.1) then
                     call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                     initial_norm=normalize(this%Gbary(:,i,j,k)-this%Lbary(:,i,j,k))
                     initial_dist=dot_product(initial_norm,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
                     call setPlane(this%liquid_gas_interface(i,j,k),0,initial_norm,initial_dist)
                     call matchVolumeFraction(neighborhood_cells(icenter),this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
                  end if
                  call setSurfaceArea(nh_r2p,surface_area)
               else
                  ! No interface was advected in our cell, use MoF
                  call reconstructMOF3D(neighborhood_cells(icenter),separated_volume_moments(icenter),this%liquid_gas_interface(i,j,k))
                  call setSurfaceArea(nh_r2p,getSA(neighborhood_cells(icenter),this%liquid_gas_interface(i,j,k)))
               end if
               
               ! Perform R2P reconstruction
               !l2g_weight=0.5_WP
               !if (this%VF(i,j,k).lt.0.1_WP) l2g_weight=1.0_WP
               !if (this%VF(i,j,k).gt.0.9_WP) l2g_weight=0.0_WP
               !l2g_weight=min(max(0.5_WP+1.25_WP*(0.5_WP-vf_nbr),0.0_WP),1.0_WP)
               !call setImportances(r2p_weight,[0.0_WP,l2g_weight,1.0_WP,-1.0_WP])
               !call setImportances(r2p_weight,[0.0_WP,l2g_weight,1.0_WP,0.0_WP])
               call reconstructR2P3D(nh_r2p,this%liquid_gas_interface(i,j,k))!,r2p_weight)
               
            end do
         end do
      end do
      
      ! Synchronize across boundaries
      call this%sync_interface()

      ! Deallocate
      deallocate(surf_norm_mag)
      
   end subroutine build_r2p


   !> Machine learning reconstruction of a planar interface in mixed cells
   subroutine build_plicnet(this)
      use mathtools, only: normalize
      use plicnet,   only: get_normal,reflect_moments
      implicit none
      class(vfs), intent(inout) :: this
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
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
               if (this%mask(i,j,k).ne.0) cycle
               ! Handle full cells differently
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                  cycle
               end if
               ! Liquid-gas symmetry
               flip=.false.
               if (this%VF(i,j,k).ge.0.5_WP) flip=.true.
               m000=0; m100=0; m010=0; m001=0
               ! Construct neighborhood of volume moments
               if (flip) then
                  do kk=k-1,k+1
                     do jj=j-1,j+1
                        do ii=i-1,i+1
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=1.0_WP-this%VF(ii,jj,kk)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                           ! Calculate geometric moments of neighborhood
                           m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                        end do
                     end do
                  end do
               else
                  do kk=k-1,k+1
                     do jj=j-1,j+1
                        do ii=i-1,i+1
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=this%VF(ii,jj,kk)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                           moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                           ! Calculate geometric moments of neighborhood
                           m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                        end do
                     end do
                  end do
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
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               call setPlane(this%liquid_gas_interface(i,j,k),0,normal,initial_dist)
               call matchVolumeFraction(cell,this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
            end do
         end do
      end do
      ! Synchronize across boundaries
      call this%sync_interface()
   end subroutine build_plicnet
   

   !> Hybrid PLICnet-R2P reconstruction of a planar interface in mixed cells
   subroutine build_r2pnet(this)
      use mathtools, only: normalize
      use plicnet,   only: get_normal,reflect_moments
      implicit none
      class(vfs), intent(inout) :: this
      integer(IRL_SignedIndex_t) :: i,j,k
      integer :: ind,ii,jj,kk,icenter
      type(R2PNeigh_RectCub_type)   :: nh_r2p
      type(RectCub_type), dimension(0:26) :: neighborhood_cells
      real(IRL_double)  , dimension(0:26) :: liquid_volume_fraction
      type(SepVM_type)  , dimension(0:26) :: separated_volume_moments
      type(VMAN_type) :: volume_moments_and_normal
      
      real(WP) :: surface_area,area
      real(WP), dimension(3) :: surface_norm
      real(WP), dimension(:,:,:), allocatable :: surf_norm_mag,tmp
      
      real(IRL_double), dimension(3) :: initial_norm
      real(IRL_double) :: initial_dist
      logical :: is_wall

      real(IRL_double), dimension(0:2) :: normal
      real(IRL_double), dimension(0:188) :: moments
      integer :: direction
      logical :: flip
      real(IRL_double) :: m000,m100,m010,m001
      real(IRL_double), dimension(0:2) :: center
      type(RectCub_type) :: cell
      
      ! Get storage for volume moments and normal
      call new(volume_moments_and_normal)
      call new(cell)
      
      ! Give ourselves an R2P neighborhood of 27 cells along with separated volume moments
      call new(nh_r2p)
      do i=0,26
         call new(neighborhood_cells(i))
         call new(separated_volume_moments(i))
      end do
      
      ! Compute magnitude of the surface-averaged normal vector
      allocate(surf_norm_mag(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); surf_norm_mag=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Skip wall/bcond/full cells
         if (this%mask(i,j,k).ne.0) cycle
         if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
         ! Extract average normal magnitude from neighborhood surface moments
         surface_area=0.0_WP; surface_norm=0.0_WP
         do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
            do ind=0,getSize(this%triangle_moments_storage(ii,jj,kk))-1
               call getMoments(this%triangle_moments_storage(ii,jj,kk),ind,volume_moments_and_normal)
               surface_area=surface_area+getVolume(volume_moments_and_normal)
               surface_norm=surface_norm+getNormal(volume_moments_and_normal)
            end do
         end do; end do; end do
         if (surface_area.gt.0.0_WP) surf_norm_mag(i,j,k)=norm2(surface_norm/surface_area)
      end do; end do; end do
      call this%cfg%sync(surf_norm_mag)
      
      ! Apply an extra step of surface smoothing to our normal magnitude
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Skip wall/bcond/full cells
         if (this%mask(i,j,k).ne.0) cycle
         if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) cycle
         ! Surface-averaged normal magnitude
         surface_area=0.0_WP
         do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
            surface_area=surface_area+this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)
            tmp(i,j,k)  =tmp(i,j,k)  +this%SD(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*surf_norm_mag(ii,jj,kk)
         end do; end do; end do
         if (surface_area.gt.0.0_WP) tmp(i,j,k)=tmp(i,j,k)/surface_area
      end do; end do; end do
      call this%cfg%sync(tmp); surf_norm_mag=tmp; deallocate(tmp)
      
      ! Traverse domain and reconstruct interface
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle
         
         ! Handle full cells differently
         if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
            call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
            call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
            cycle
         end if
         
         ! If a wall is in our neighborhood, apply PLICNET
         !!! ALSO FORCING PLICNET IN MY NOZZLE HERE
         is_wall=.false.
         do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
            if (this%mask(ii,jj,kk).eq.1) is_wall=.true.
         end do; end do; end do
         if (is_wall.or.this%cfg%xm(i).lt.0.0001_WP) then
            ! PLICNET
            ! Liquid-gas symmetry
            flip=.false.; if (this%VF(i,j,k).ge.0.5_WP) flip=.true.
            m000=0; m100=0; m010=0; m001=0
            ! Construct neighborhood of volume moments
            if (flip) then
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=1.0_WP-this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  ! Calculate geometric moments of neighborhood
                  m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
               end do; end do; end do
            else
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
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
            call get_normal(moments,normal); normal=normalize(normal)
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
            call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
            call setPlane(this%liquid_gas_interface(i,j,k),0,normal,initial_dist)
            call matchVolumeFraction(cell,this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
            ! Done with that cell
            cycle
         end if
         
         ! If the neighborhood normals are sufficiently consistent, just use PLICNET
         if (surf_norm_mag(i,j,k).gt.this%twoplane_thld2) then
            ! PLICNET
            ! Liquid-gas symmetry
            flip=.false.; if (this%VF(i,j,k).ge.0.5_WP) flip=.true.
            m000=0; m100=0; m010=0; m001=0
            ! Construct neighborhood of volume moments
            if (flip) then
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=1.0_WP-this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  ! Calculate geometric moments of neighborhood
                  m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                  m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
               end do; end do; end do
            else
               do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k)))=this%VF(ii,jj,kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(this%Lbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(this%Lbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(this%Lbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(this%Gbary(1,ii,jj,kk)-this%cfg%xm(ii))/this%cfg%dx(ii)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(this%Gbary(2,ii,jj,kk)-this%cfg%ym(jj))/this%cfg%dy(jj)
                  moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(this%Gbary(3,ii,jj,kk)-this%cfg%zm(kk))/this%cfg%dz(kk)
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
            call get_normal(moments,normal); normal=normalize(normal)
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
            call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
            call setPlane(this%liquid_gas_interface(i,j,k),0,normal,initial_dist)
            call matchVolumeFraction(cell,this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
            ! Done with that cell
            cycle
         end if
         
         ! Prepare R2P data
         ind=0; call emptyNeighborhood(nh_r2p)
         do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
            call addMember(nh_r2p,neighborhood_cells(ind),separated_volume_moments(ind))
            call construct_2pt(neighborhood_cells(ind),[this%cfg%x(ii),this%cfg%y(jj),this%cfg%z(kk)],[this%cfg%x(ii+1),this%cfg%y(jj+1),this%cfg%z(kk+1)])
            call construct(separated_volume_moments(ind),[this%VF(ii,jj,kk)*this%cfg%vol(ii,jj,kk),this%Lbary(:,ii,jj,kk),(1.0_WP-this%VF(ii,jj,kk))*this%cfg%vol(ii,jj,kk),this%Gbary(:,ii,jj,kk)])
            if (ii.eq.i.and.jj.eq.j.and.kk.eq.k) then
               icenter=ind
               call setCenterOfStencil(nh_r2p,icenter)
            end if
            ind=ind+1
         end do; end do; end do
         
         ! Generate initial guess for R2P based on availability of in-cell surface data
         surface_area=0.0_WP
         do ind=0,getSize(this%triangle_moments_storage(i,j,k))-1
            call getMoments(this%triangle_moments_storage(i,j,k),ind,volume_moments_and_normal)
            surface_area=surface_area+getVolume(volume_moments_and_normal)
         end do
         if (surface_area.gt.surface_epsilon_factor*this%cfg%meshsize(i,j,k)**2) then
            ! Local normals are available, reconstruction from surface data
            call reconstructAdvectedNormals(this%triangle_moments_storage(i,j,k),nh_r2p,this%twoplane_thld1,this%liquid_gas_interface(i,j,k))
            if (getNumberOfPlanes(this%liquid_gas_interface(i,j,k)).eq.1) then
               call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
               initial_norm=normalize(this%Gbary(:,i,j,k)-this%Lbary(:,i,j,k))
               initial_dist=dot_product(initial_norm,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
               call setPlane(this%liquid_gas_interface(i,j,k),0,initial_norm,initial_dist)
               call matchVolumeFraction(neighborhood_cells(icenter),this%VF(i,j,k),this%liquid_gas_interface(i,j,k))
            end if
            call setSurfaceArea(nh_r2p,surface_area)
         else
            ! No interface was advected in our cell, use MoF
            call reconstructMOF3D(neighborhood_cells(icenter),separated_volume_moments(icenter),this%liquid_gas_interface(i,j,k))
            call setSurfaceArea(nh_r2p,getSA(neighborhood_cells(icenter),this%liquid_gas_interface(i,j,k)))
         end if
         
         ! Perform R2P reconstruction
         call reconstructR2P3D(nh_r2p,this%liquid_gas_interface(i,j,k))
         
      end do; end do; end do
      
      ! Synchronize across boundaries
      call this%sync_interface()
      
      ! Deallocate
      deallocate(surf_norm_mag)
      
   end subroutine build_r2pnet

   
   !> Set all domain boundaries to full liquid/gas based on VOF value
   subroutine set_full_bcond(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      ! In X-
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino,this%cfg%imin-1
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               end do
            end do
         end do
      end if
      ! In X+
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imax+1,this%cfg%imaxo
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               end do
            end do
         end do
      end if
      ! In Y-
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino,this%cfg%jmin-1
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               end do
            end do
         end do
      end if
      ! In Y+
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmax+1,this%cfg%jmaxo
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               end do
            end do
         end do
      end if
      ! In Z-
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         do k=this%cfg%kmino,this%cfg%kmin-1
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               end do
            end do
         end do
      end if
      ! In Z+
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         do k=this%cfg%kmax+1,this%cfg%kmaxo
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               end do
            end do
         end do
      end if
   end subroutine set_full_bcond
   
   
   !> Polygonalization of the IRL interface (calculates SD at the same time)
   !> Here, only mask=1 is skipped (i.e., real walls), so bconds should be handled
   subroutine polygonalize_interface(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,n
      real(WP) :: tsd
      type(RectCub_type) :: cell
      real(WP), dimension(1:3,1:4) :: vert
      real(WP), dimension(1:3) :: norm
      
      ! Create a cell object
      call new(cell)
      
      ! Loop over full domain and form polygon
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Zero out the polygons
               do n=1,max_interface_planes
                  call zeroPolygon(this%interface_polygon(n,i,j,k))
               end do
               ! Skip wall cells only here
               if (this%mask(i,j,k).eq.1) cycle
               ! Create polygons for cells with interfaces, zero for those without
               if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) then
                  call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
                  do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                     call getPoly(cell,this%liquid_gas_interface(i,j,k),n-1,this%interface_polygon(n,i,j,k))
                  end do
               end if
            end do
         end do
      end do
      
      ! Find inferface between filled and empty cells on x-face
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               if (this%VF(i,j,k).lt.VFlo.and.this%VF(i-1,j,k).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i-1,j,k).lt.VFlo) then
                  if (this%mask(i,j,k).eq.1.or.this%mask(i-1,j,k).eq.1) cycle
                  norm=[sign(1.0_WP,0.5_WP-this%VF(i,j,k)),0.0_WP,0.0_WP]
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%x(i))
                  vert(:,1)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k  )]
                  vert(:,2)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k  )]
                  vert(:,3)=[this%cfg%x(i),this%cfg%y(j+1),this%cfg%z(k+1)]
                  vert(:,4)=[this%cfg%x(i),this%cfg%y(j  ),this%cfg%z(k+1)]
                  call construct(this%interface_polygon(1,i,j,k),4,vert)
                  call setPlaneOfExistence(this%interface_polygon(1,i,j,k),getPlane(this%liquid_gas_interface(i,j,k),0))
                  if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(1,i,j,k))
               end if
            end do
         end do
      end do
      
      ! Find inferface between filled and empty cells on y-face
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%VF(i,j,k).lt.VFlo.and.this%VF(i,j-1,k).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i,j-1,k).lt.VFlo) then
                  if (this%mask(i,j,k).eq.1.or.this%mask(i,j-1,k).eq.1) cycle
                  norm=[0.0_WP,sign(1.0_WP,0.5_WP-this%VF(i,j,k)),0.0_WP]
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%y(j))
                  vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k  )]
                  vert(:,2)=[this%cfg%x(i  ),this%cfg%y(j),this%cfg%z(k+1)]
                  vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k+1)]
                  vert(:,4)=[this%cfg%x(i+1),this%cfg%y(j),this%cfg%z(k  )]
                  call construct(this%interface_polygon(1,i,j,k),4,vert)
                  call setPlaneOfExistence(this%interface_polygon(1,i,j,k),getPlane(this%liquid_gas_interface(i,j,k),0))
                  if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(1,i,j,k))
               end if
            end do
         end do
      end do
      
      ! Find inferface between filled and empty cells on z-face
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%VF(i,j,k).lt.VFlo.and.this%VF(i,j,k-1).gt.VFhi.or.this%VF(i,j,k).gt.VFhi.and.this%VF(i,j,k-1).lt.VFlo) then
                  if (this%mask(i,j,k).eq.1.or.this%mask(i,j,k-1).eq.1) cycle
                  norm=[0.0_WP,0.0_WP,sign(1.0_WP,0.5_WP-this%VF(i,j,k))]
                  call setNumberOfPlanes(this%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%liquid_gas_interface(i,j,k),0,norm,sign(1.0_WP,0.5_WP-this%VF(i,j,k))*this%cfg%z(k))
                  vert(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k)]
                  vert(:,2)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k)]
                  vert(:,3)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k)]
                  vert(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k)]
                  call construct(this%interface_polygon(1,i,j,k),4,vert)
                  call setPlaneOfExistence(this%interface_polygon(1,i,j,k),getPlane(this%liquid_gas_interface(i,j,k),0))
                  if (this%VF(i,j,k).gt.VFhi) call reversePtOrdering(this%interface_polygon(1,i,j,k))
               end if
            end do
         end do
      end do
      
      ! Now compute surface area divided by cell volume
      this%SD=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) cycle
               tsd=0.0_WP
               do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%interface_polygon(n,i,j,k)).gt.0) then
                     tsd=tsd+abs(calculateVolume(this%interface_polygon(n,i,j,k)))
                  end if
               end do
               this%SD(i,j,k)=tsd/this%cfg%vol(i,j,k)
            end do
         end do
      end do
      
   end subroutine polygonalize_interface
   
   
   !> Update a surfmesh object from our current polygons
   subroutine update_surfmesh(this,smesh)
      use surfmesh_class, only: surfmesh
      implicit none
      class(vfs), intent(inout) :: this
      class(surfmesh), intent(inout) :: smesh
      integer :: i,j,k,n,shape,nv,np,nplane
      real(WP), dimension(3) :: tmp_vert
      
      ! Reset surface mesh storage
      call smesh%reset()
      
      ! First pass to count how many vertices and polygons are inside our processor
      nv=0; np=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do nplane=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  shape=getNumberOfVertices(this%interface_polygon(nplane,i,j,k))
                  if (shape.gt.0) then
                     nv=nv+shape
                     np=np+1
                  end if
               end do
            end do
         end do
      end do
      
      ! Reallocate storage and fill out arrays
      if (np.gt.0) then
         call smesh%set_size(nvert=nv,npoly=np)
         allocate(smesh%polyConn(smesh%nVert)) ! Also allocate naive connectivity
         nv=0; np=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(this%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        smesh%polySize(np)=shape
                        ! Loop over its vertices and add them
                        do n=1,shape
                           tmp_vert=getPt(this%interface_polygon(nplane,i,j,k),n-1)
                           ! Increment node counter
                           nv=nv+1
                           smesh%xVert(nv)=tmp_vert(1)
                           smesh%yVert(nv)=tmp_vert(2)
                           smesh%zVert(nv)=tmp_vert(3)
                           smesh%polyConn(nv)=nv
                        end do
                     end if
                  end do
               end do
            end do
         end do
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


   !> Update a surfmesh object from our current polygons - near-empty cells are not shown
   subroutine update_surfmesh_nowall(this,smesh,threshold)
      use surfmesh_class, only: surfmesh
      implicit none
      class(vfs), intent(inout) :: this
      class(surfmesh), intent(inout) :: smesh
      real(WP), optional :: threshold
      real(WP) :: VFclip
      integer :: i,j,k,n,shape,nv,np,nplane
      real(WP), dimension(3) :: tmp_vert
      
      ! Handle VF threshold
      if (present(threshold)) then
         VFclip=threshold
      else
         VFclip=2.0_WP*epsilon(1.0_WP)
      end if

      ! Reset surface mesh storage
      call smesh%reset()
      
      ! First pass to count how many vertices and polygons are inside our processor
      nv=0; np=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%VF(i,j,k).lt.VFclip) cycle ! Skip cells below VF threshold
               do nplane=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  shape=getNumberOfVertices(this%interface_polygon(nplane,i,j,k))
                  if (shape.gt.0) then
                     nv=nv+shape
                     np=np+1
                  end if
               end do
            end do
         end do
      end do
      
      ! Reallocate storage and fill out arrays
      if (np.gt.0) then
         call smesh%set_size(nvert=nv,npoly=np)
         allocate(smesh%polyConn(smesh%nVert)) ! Also allocate naive connectivity
         nv=0; np=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%cfg%VF(i,j,k).lt.VFclip) cycle ! Skip cells below VF threshold
                  do nplane=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(this%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        smesh%polySize(np)=shape
                        ! Loop over its vertices and add them
                        do n=1,shape
                           tmp_vert=getPt(this%interface_polygon(nplane,i,j,k),n-1)
                           ! Increment node counter
                           nv=nv+1
                           smesh%xVert(nv)=tmp_vert(1)
                           smesh%yVert(nv)=tmp_vert(2)
                           smesh%zVert(nv)=tmp_vert(3)
                           smesh%polyConn(nv)=nv
                        end do
                     end if
                  end do
               end do
            end do
         end do
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
      
   end subroutine update_surfmesh_nowall
   
   
   !> Calculate distance from polygonalized interface inside the band
   !> Domain edges are not done here
   subroutine distance_from_polygon(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: ni,i,j,k,ii,jj,kk,index
      real(IRL_double), dimension(3) :: pos,nearest_pt
      
      ! First reset distance
      this%G=huge(1.0_WP)
      
      ! Loop over 1/2-band
      do index=1,sum(this%band_count(0:distance_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         ! Get cell centroid location
         pos=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
         ! Loop over neighboring polygons and compute distance
         do kk=k-2,k+2
            do jj=j-2,j+2
               do ii=i-2,i+2
                  do ni=1,getNumberOfPlanes(this%liquid_gas_interface(ii,jj,kk))
                     if (getNumberOfVertices(this%interface_polygon(ni,ii,jj,kk)).ne.0) then
                        nearest_pt=calculateNearestPtOnSurface(this%interface_polygon(ni,ii,jj,kk),pos)
                        nearest_pt=pos-nearest_pt
                        this%G(i,j,k)=min(this%G(i,j,k),dot_product(nearest_pt,nearest_pt))
                     end if
                  end do
               end do
            end do
         end do
         this%G(i,j,k)=sqrt(this%G(i,j,k))
         ! Only need to consult planes in own cell to know sign
         ! Even "empty" cells have one plane, which is really far away from it..
         if (.not.isPtInt(pos,this%liquid_gas_interface(i,j,k))) this%G(i,j,k)=-this%G(i,j,k)
      end do
      
      ! Clip distance field and sign it properly
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (abs(this%G(i,j,k)).gt.this%Gclip) this%G(i,j,k)=sign(this%Gclip,this%VF(i,j,k)-0.5_WP)
            end do
         end do
      end do
      
      ! Sync boundaries
      call this%cfg%sync(this%G)
      
   end subroutine distance_from_polygon
   
   
   !> Calculate subcell phasic volumes from reconstructed interface
   !> Here, only mask=1 is skipped (i.e., real walls), so bconds are handled
   subroutine subcell_vol(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk
      real(WP), dimension(0:2) :: subx,suby,subz
      type(RectCub_type) :: cell
      type(SepVM_type) :: separated_volume_moments
      
      ! Allocate IRL objects for moment calculation
      call new(cell)
      call new(separated_volume_moments)
      
      ! Compute subcell liquid and gas information
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Deal with walls only - we do compute inside bconds here
               if (this%mask(i,j,k).eq.1) then
                  this%Lvol(:,:,:,i,j,k)=0.0_WP
                  this%Gvol(:,:,:,i,j,k)=0.0_WP
                  cycle
               end if
               ! Deal with other cells
               if (this%VF(i,j,k).gt.VFhi) then
                  this%Lvol(:,:,:,i,j,k)=0.125_WP*this%cfg%vol(i,j,k)
                  this%Gvol(:,:,:,i,j,k)=0.0_WP
               else if (this%VF(i,j,k).lt.VFlo) then
                  this%Lvol(:,:,:,i,j,k)=0.0_WP
                  this%Gvol(:,:,:,i,j,k)=0.125_WP*this%cfg%vol(i,j,k)
               else
                  ! Prepare subcell extent
                  subx=[this%cfg%x(i),this%cfg%xm(i),this%cfg%x(i+1)]
                  suby=[this%cfg%y(j),this%cfg%ym(j),this%cfg%y(j+1)]
                  subz=[this%cfg%z(k),this%cfg%zm(k),this%cfg%z(k+1)]
                  ! Loop over sub-cells
                  do kk=0,1
                     do jj=0,1
                        do ii=0,1
                           call construct_2pt(cell,[subx(ii),suby(jj),subz(kk)],[subx(ii+1),suby(jj+1),subz(kk+1)])
                           call getNormMoments(cell,this%liquid_gas_interface(i,j,k),separated_volume_moments)
                           this%Lvol(ii,jj,kk,i,j,k)=getVolume(separated_volume_moments,0)
                           this%Gvol(ii,jj,kk,i,j,k)=getVolume(separated_volume_moments,1)
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      
   end subroutine subcell_vol
   
   
   !> Reset volumetric moments based on reconstructed interface
   !> NGA finishes this with a comm - I removed it as it does not seem useful
   subroutine reset_volume_moments(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      type(RectCub_type) :: cell
      type(SepVM_type) :: separated_volume_moments
      
      ! Calculate volume moments and store
      call new(cell)
      call new(separated_volume_moments)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Handle pure wall cells
               if (this%mask(i,j,k).eq.1) then
                  this%VF(i,j,k)     =0.0_WP
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  cycle
               end if
               ! Form the grid cell
               call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
               ! Cut it by the current interface(s)
               call getNormMoments(cell,this%liquid_gas_interface(i,j,k),separated_volume_moments)
               ! Recover relevant moments
               this%VF(i,j,k)     =getVolumePtr(separated_volume_moments,0)/this%cfg%vol(i,j,k)
               this%Lbary(:,i,j,k)=getCentroid(separated_volume_moments,0)
               this%Gbary(:,i,j,k)=getCentroid(separated_volume_moments,1)
               ! Clean up
               if (this%VF(i,j,k).lt.VFlo) then
                  this%VF(i,j,k)     =0.0_WP
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               end if
               if (this%VF(i,j,k).gt.VFhi) then
                  this%VF(i,j,k)     =1.0_WP
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               end if
            end do
         end do
      end do
      
      ! NGA had comms here: unclear to me why it would be necessary...
      ! Synchronize VF field
      !call this%cfg%sync(this%VF)
      ! Synchronize and clean-up barycenter fields
      !call this%sync_and_clean_barycenters()
      
   end subroutine reset_volume_moments
   

   ! Reset only moments, leave VF unchanged
   subroutine reset_moments(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      type(RectCub_type) :: cell
      type(SepVM_type) :: separated_volume_moments
      
      ! Calculate volume moments and store
      call new(cell)
      call new(separated_volume_moments)
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Handle pure wall cells
               if (this%mask(i,j,k).eq.1) then
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  cycle
               end if
               ! Reset first order moments based on VOF value
               if (this%VF(i,j,k).lt.VFlo) then
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               else if (this%VF(i,j,k).gt.VFhi) then
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               else
                  ! Form the grid cell
                  call construct_2pt(cell,[this%cfg%x(i),this%cfg%y(j),this%cfg%z(k)],[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
                  ! Cut it by the current interface(s)
                  call getNormMoments(cell,this%liquid_gas_interface(i,j,k),separated_volume_moments)
                  ! Recover relevant moments
                  this%Lbary(:,i,j,k)=getCentroid(separated_volume_moments,0)
                  this%Gbary(:,i,j,k)=getCentroid(separated_volume_moments,1)
               end if
            end do
         end do
      end do
      
   end subroutine reset_moments
   
   
   !> Compute curvature from a least squares fit of the IRL surface
   subroutine get_curvature(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,n
      real(WP), dimension(max_interface_planes) :: mycurv,mysurf
      real(WP), dimension(max_interface_planes,3) :: mynorm
      real(WP), dimension(3) :: csn,sn
      ! Reset curvature
      this%curv=0.0_WP
      if (this%two_planes) this%curv2p=0.0_WP
      ! Traverse interior domain and compute curvature in cells with polygons
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Zero out curvature and surface storage
               mycurv=0.0_WP; mysurf=0.0_WP; mynorm=0.0_WP
               ! Get a curvature for each plane
               do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  ! Skip empty polygon
                  if (getNumberOfVertices(this%interface_polygon(n,i,j,k)).eq.0) cycle
                  ! Perform LSQ PLIC barycenter fitting to get curvature
                  !call this%paraboloid_fit(i,j,k,n,mycurv(n))
                  ! Perform PLIC surface fitting to get curvature
                  call this%paraboloid_integral_fit(i,j,k,n,mycurv(n))
                  ! Also store surface and normal
                  mysurf(n)  =abs(calculateVolume(this%interface_polygon(n,i,j,k)))
                  mynorm(n,:)=    calculateNormal(this%interface_polygon(n,i,j,k))
               end do
               ! Oriented-surface-average curvature
               !csn=0.0_WP; sn=0.0_WP
               !do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
               !   csn=csn+mysurf(n)*mynorm(n,:)*mycurv(n)
               !   sn = sn+mysurf(n)*mynorm(n,:)
               !end do
               !if (dot_product(sn,sn).gt.10.0_WP*tiny(1.0_WP)) this%curv(i,j,k)=dot_product(csn,sn)/dot_product(sn,sn)
               ! Surface-averaged curvature
               if (sum(mysurf).gt.0.0_WP) this%curv(i,j,k)=sum(mysurf*mycurv)/sum(mysurf)
               ! Curvature of largest surface
               !if (mysurf(maxloc(mysurf,1)).gt.0.0_WP) this%curv(i,j,k)=mycurv(maxloc(mysurf,1))
               ! Largest curvature
               !this%curv(i,j,k)=mycurv(maxloc(abs(mycurv),1))
               ! Smallest curvature
               !if (getNumberOfPlanes(this%liquid_gas_interface(i,j,k)).eq.2) then
               !   this%curv(i,j,k)=mycurv(minloc(abs(mycurv),1))
               !else
               !   this%curv(i,j,k)=mycurv(1)
               !end if
               ! Clip curvature - may not be needed if we select polygons carefully
               this%curv(i,j,k)=max(min(this%curv(i,j,k),this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k)),-this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k))
               ! Also store 2-plane curvature if needed
               if (this%two_planes) this%curv2p(:,i,j,k)=max(min(mycurv,this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k)),-this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k))
               ! Model edge curvature at 1/thickness
               !if (this%two_planes.and.this%edge_sensor(i,j,k).gt.this%edge_thld) this%curv2p(:,i,j,k)=sign(1.0_WP/this%thickness(i,j,k),0.5_WP-this%VF(i,j,k))
            end do
         end do
      end do
      ! Synchronize boundaries
      call this%cfg%sync(this%curv)
      if (this%two_planes) call this%cfg%sync(this%curv2p)
   end subroutine get_curvature
   
   
   !> Perform local paraboloid fit of IRL surface in pointwise sense
   subroutine paraboloid_fit(this,i,j,k,iplane,mycurv)
      use mathtools, only: normalize,cross_product
      implicit none
      ! In/out variables
      class(vfs), intent(inout) :: this
      integer,  intent(in)  :: i,j,k,iplane
      real(WP), intent(out) :: mycurv
      ! Variables used to process the polygonal surface
      real(WP), dimension(3) :: pref,nref,tref,sref
      real(WP), dimension(3) :: ploc,nloc
      real(WP), dimension(3) :: buf
      real(WP) :: surf,ww
      integer :: n,ii,jj,kk,ndata,info
      ! Storage for least squares problem
      real(WP), dimension(125,6) :: A=0.0_WP
      real(WP), dimension(125)   :: b=0.0_WP
      real(WP), dimension(6)     :: sol
      real(WP), dimension(200)   :: work
      ! Curvature evaluation
      real(WP) :: dF_dt,dF_ds,ddF_dtdt,ddF_dsds,ddF_dtds
      
      ! Store polygon centroid - this is our reference point
      pref=calculateCentroid(this%interface_polygon(iplane,i,j,k))
      
      ! Create local basis from polygon normal
      nref=calculateNormal(this%interface_polygon(iplane,i,j,k))
      select case (maxloc(abs(nref),1))
      case (1); tref=normalize([+nref(2),-nref(1),0.0_WP])
      case (2); tref=normalize([0.0_WP,+nref(3),-nref(2)])
      case (3); tref=normalize([-nref(3),0.0_WP,+nref(1)])
      end select; sref=cross_product(nref,tref)
      
      ! Collect all data
      ndata=0
      do kk=k-2,k+2
         do jj=j-2,j+2
            do ii=i-2,i+2
               
               ! Skip the cell if it's a true wall
               if (this%mask(ii,jj,kk).eq.1) cycle
               
               ! Check all planes
               do n=1,getNumberOfPlanes(this%liquid_gas_interface(ii,jj,kk))
                  
                  ! Skip empty polygon
                  if (getNumberOfVertices(this%interface_polygon(n,ii,jj,kk)).eq.0) cycle
                  
                  ! Get local polygon normal
                  nloc=calculateNormal(this%interface_polygon(n,ii,jj,kk))
                  
                  ! Store triangle centroid, and surface
                  ploc=    calculateCentroid(this%interface_polygon(n,ii,jj,kk))
                  surf=abs(calculateVolume  (this%interface_polygon(n,ii,jj,kk)))/this%cfg%meshsize(i,j,k)**2
                  
                  ! Transform polygon barycenter to a local coordinate system
                  buf=(ploc-pref)/this%cfg%meshsize(i,j,k); ploc=[dot_product(buf,nref),dot_product(buf,tref),dot_product(buf,sref)]
                  
                  ! Distance from ref point AND projected surface weighting (clipped to ensure positivity)
                  ww=surf*max(dot_product(nloc,nref),0.0_WP)*wgauss(sqrt(dot_product(ploc,ploc)),2.5_WP)
                  
                  ! If we have data, add it to the LS problem
                  if (ww.gt.0.0_WP) then
                     ! Increment counter
                     ndata=ndata+1
                     ! Store least squares matrix and RHS
                     A(ndata,1)=sqrt(ww)*1.0_WP
                     A(ndata,2)=sqrt(ww)*ploc(2)
                     A(ndata,3)=sqrt(ww)*ploc(3)
                     A(ndata,4)=sqrt(ww)*0.5_WP*ploc(2)*ploc(2)
                     A(ndata,5)=sqrt(ww)*0.5_WP*ploc(3)*ploc(3)
                     A(ndata,6)=sqrt(ww)*1.0_WP*ploc(2)*ploc(3)
                     b(ndata  )=sqrt(ww)*ploc(1)
                  end if
                  
               end do
               
            end do
         end do
      end do
      
      ! Solve for surface as n=F(t,s)=b1+b2*t+b3*s+b4*t^2+b5*s^2+b6*t*s using Lapack
      call dgels('N',ndata,6,1,A,125,b,125,work,200,info); sol=b(1:6)
      
      ! Get the curvature at (t,s)=(0,0)
      dF_dt=sol(2)+sol(4)*0.0_WP+sol(6)*0.0_WP; ddF_dtdt=sol(4)
      dF_ds=sol(3)+sol(5)*0.0_WP+sol(6)*0.0_WP; ddF_dsds=sol(5)
      ddF_dtds=sol(6)
      mycurv=-((1.0_WP+dF_dt**2)*ddF_dsds-2.0_WP*dF_dt*dF_ds*ddF_dtds+(1.0_WP+dF_ds**2)*ddF_dtdt)/(1.0_WP+dF_dt**2+dF_ds**2)**(1.5_WP)
      mycurv=mycurv/this%cfg%meshsize(i,j,k)
      
   contains
      
      ! Some weighting function - h=0.75 looks okay
      real(WP) function wkernel(d,h)
         implicit none
         real(WP), intent(in) :: d,h
         wkernel=(1.0_WP+(d/h)**2)**(-1.4_WP)
      end function wkernel
      
      ! Tri-cubic Weighting function - h=2 looks okay
      real(WP) function tricubic(d,h)
         implicit none
         real(WP), intent(in) :: d,h
         if (d.ge.h) then
            tricubic=0.0_WP
         else
            tricubic=(1.0_WP-(d/h)**3)**3
         end if
      end function tricubic
      
      ! Quasi-Gaussian weighting function - h=2.5 looks okay
      real(WP) function wgauss(d,h)
         implicit none
         real(WP), intent(in) :: d,h
         if (d.ge.h) then
            wgauss=0.0_WP
         else
            wgauss=(1.0_WP+4.0_WP*d/h)*(1.0_WP-d/h)**4
         end if
      end function wgauss
      
   end subroutine paraboloid_fit


   !> Perform local paraboloid fit of IRL surface in integral sense
   subroutine paraboloid_integral_fit(this,i,j,k,iplane,mycurv)
      use mathtools, only: normalize,cross_product
      implicit none
      ! In/out variables
      class(vfs), intent(inout) :: this
      integer,  intent(in)  :: i,j,k,iplane
      real(WP), intent(out) :: mycurv
      ! Variables used to process the polygons
      real(WP), dimension(3) :: pref,nref,tref,sref
      real(WP), dimension(3) :: vert1,vert2,ploc,nloc
      real(WP), dimension(3) :: buf,reconst_plane_coeffs
      integer :: nplane,shape,n,ii,jj,kk,ai,aj
      real(WP), dimension(6) :: integrals
      real(WP) :: xv,xvn,yv,yvn,ww,b_dot_sum
      ! Storage for symmetric problem
      real(WP), dimension(6,6) :: A
      integer , dimension(6)   :: ipiv
      real(WP), dimension(6)   :: b
      real(WP), dimension(6)   :: sol
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      integer  :: lwork,info
      ! Curvature evaluation
      real(WP) :: dF_dt,dF_ds,ddF_dtdt,ddF_dsds,ddF_dtds
      
      ! Store polygon centroid - this is our reference point
      pref=calculateCentroid(this%interface_polygon(iplane,i,j,k))
      
      ! Create local basis from polygon normal
      nref=calculateNormal(this%interface_polygon(iplane,i,j,k))
      select case (maxloc(abs(nref),1))
      case (1); tref=normalize([+nref(2),-nref(1),0.0_WP])
      case (2); tref=normalize([0.0_WP,+nref(3),-nref(2)])
      case (3); tref=normalize([-nref(3),0.0_WP,+nref(1)])
      end select; sref=cross_product(nref,tref)
      
      ! Collect all data
      A=0.0_WP
      b=0.0_WP
      do kk=k-2,k+2
         do jj=j-2,j+2
            do ii=i-2,i+2
               
               ! Skip the cell if it's a true wall
               if (this%mask(ii,jj,kk).eq.1) cycle
               
               ! Check all planes
               do nplane=1,getNumberOfPlanes(this%liquid_gas_interface(ii,jj,kk))
                  
                  ! Skip empty polygon
                  shape=getNumberOfVertices(this%interface_polygon(nplane,ii,jj,kk))
                  if (shape.eq.0) cycle
                  
                  ! Get local polygon normal and skip if normal is not aligned with center polygon normal
                  nloc=calculateNormal(this%interface_polygon(nplane,ii,jj,kk))
                  if (dot_product(nloc,nref).le.0.0_WP) cycle
                  
                  ! Get local polygon centroid
                  ploc=calculateCentroid(this%interface_polygon(nplane,ii,jj,kk))
                  
                  ! Transform normal and centroid to a local coordinate system
                  buf=(ploc-pref)/this%cfg%meshsize(i,j,k); ploc=[dot_product(buf,tref),dot_product(buf,sref),dot_product(buf,nref)]
                  buf=nloc; nloc=[dot_product(buf,tref),dot_product(buf,sref),dot_product(buf,nref)]
                  
                  ! Get plane coefficients
                  reconst_plane_coeffs(1)=-dot_product(nloc,ploc)
                  reconst_plane_coeffs(2)=nloc(1)
                  reconst_plane_coeffs(3)=nloc(2)
                  reconst_plane_coeffs=reconst_plane_coeffs/(-nloc(3))
                  
                  ! Get integrals
                  integrals=0.0_WP
                  b_dot_sum=0.0_WP
                  do n=1,shape
                     vert1=getPt(this%interface_polygon(nplane,ii,jj,kk),n-1)
                     vert2=getPt(this%interface_polygon(nplane,ii,jj,kk),modulo(n,shape))
                     ! Transform vertices to a local coordinate system
                     buf=(vert1-pref)/this%cfg%meshsize(i,j,k); vert1=[dot_product(buf,tref),dot_product(buf,sref),dot_product(buf,nref)]
                     buf=(vert2-pref)/this%cfg%meshsize(i,j,k); vert2=[dot_product(buf,tref),dot_product(buf,sref),dot_product(buf,nref)]
                     ! Add to area integral
                     xv=vert1(1); xvn=vert2(1); yv=vert1(2); yvn=vert2(2)
                     integrals = integrals + [&
                     (xv*yvn - xvn*yv) / 2.0_WP, &
                     (xv + xvn)*(xv*yvn - xvn*yv) / 6.0_WP, &
                     (yv + yvn)*(xv*yvn - xvn*yv) / 6.0_WP, &
                     (xv + xvn)*(xv**2 + xvn**2)*(yvn - yv) / 12.0_WP, &
                     (yvn - yv)*(3.0_WP*xv**2*yv + xv**2*yvn + 2.0_WP*xv*xvn*yv + 2.0_WP*xv*xvn*yvn + xvn**2*yv + 3.0_WP*xvn**2*yvn)/24.0_WP, &
                     (xv - xvn)*(yv + yvn)*(yv**2 + yvn**2) / 12.0_WP]
                  end do
                  b_dot_sum=b_dot_sum+dot_product(reconst_plane_coeffs,integrals(1:3))
                  
                  ! Get weighting
                  ww=wgauss(sqrt(dot_product(ploc,ploc)),2.5_WP)
                  
                  ! Add to symmetric matrix and RHS
                  do aj=1,6
                     do ai=1,aj
                        A(ai,aj)=A(ai,aj)+ww*integrals(ai)*integrals(aj)
                     end do
                  end do
                  b=b+ww*integrals*b_dot_sum
                  
               end do
            end do
         end do
      end do
      
      ! Query optimal work array size then solve for paraboloid as n=F(t,s)=b1+b2*t+b3*s+b4*t^2+b5*t*s+b6*s^2
      call dsysv('U',6,1,A,6,ipiv,b,6,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      call dsysv('U',6,1,A,6,ipiv,b,6,work,lwork,info); sol=b(1:6); deallocate(work)
      
      ! Get the curvature at (t,s)=(0,0)
      dF_dt=sol(2)+2.0_WP*sol(4)*0.0_WP+sol(5)*0.0_WP; ddF_dtdt=2.0_WP*sol(4)
      dF_ds=sol(3)+2.0_WP*sol(6)*0.0_WP+sol(5)*0.0_WP; ddF_dsds=2.0_WP*sol(6)
      ddF_dtds=sol(5)
      mycurv=-((1.0_WP+dF_dt**2)*ddF_dsds-2.0_WP*dF_dt*dF_ds*ddF_dtds+(1.0_WP+dF_ds**2)*ddF_dtdt)/(1.0_WP+dF_dt**2+dF_ds**2)**(1.5_WP)
      mycurv=mycurv/this%cfg%meshsize(i,j,k)
      
   contains
      
      ! Quasi-Gaussian weighting function - h=2.5 looks okay
      real(WP) function wgauss(d,h)
         implicit none
         real(WP), intent(in) :: d,h
         if (d.ge.h) then
            wgauss=0.0_WP
         else
            wgauss=(1.0_WP+4.0_WP*d/h)*(1.0_WP-d/h)**4
         end if
      end function wgauss
      
   end subroutine paraboloid_integral_fit
   
   
   !> Private function to rapidly assess if a mixed cell is possible
   pure function crude_phase_test(this,b_ind) result(crude_phase)
      implicit none
      class(vfs), intent(in) :: this
      integer, dimension(3,2), intent(in) :: b_ind
      real(WP) :: crude_phase
      integer :: i,j,k
      ! Check if originating cell is mixed, continue if not
      crude_phase=this%VF(b_ind(1,2),b_ind(2,2),b_ind(3,2))
      if (crude_phase.ge.VFlo.and.crude_phase.le.VFhi) then
         ! Already have a mixed cell, we need the full geometry
         crude_phase=-1.0_WP; return
      end if
      ! Check cells in bounding box
      do k=b_ind(3,1),b_ind(3,2)
         do j=b_ind(2,1),b_ind(2,2)
            do i=b_ind(1,1),b_ind(1,2)
               if (this%VF(i,j,k).ne.crude_phase) then
                  ! We could have changed phase, we need the full geometry
                  crude_phase=-1.0_WP; return
               end if
            end do
         end do
      end do
      ! Ensure proper values
      if (crude_phase.gt.VFhi) then
         crude_phase=1.0_WP
      else if (crude_phase.lt.VFlo) then
         crude_phase=0.0_WP
      end if
   end function crude_phase_test
   
   
   !> Private function that performs a Lagrangian projection of a vertex p1 to position p2
   !> using the provided velocity U/V/W, time step dt, and a guess of the i/j/k
   function project(this,p1,i,j,k,dt,U,V,W) result(p2)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), dimension(3), intent(in) :: p1
      integer,                intent(in) :: i,j,k
      real(WP), intent(in) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(3) :: p2
      ! Explicit RK4
      !real(WP), dimension(3) :: v1,v2,v3,v4
      !v1=this%cfg%get_velocity(p1             ,i,j,k,U,V,W)
      !v2=this%cfg%get_velocity(p1+0.5_WP*dt*v1,i,j,k,U,V,W)
      !v3=this%cfg%get_velocity(p1+0.5_WP*dt*v2,i,j,k,U,V,W)
      !v4=this%cfg%get_velocity(p1+       dt*v3,i,j,k,U,V,W)
      !p2=p1+dt/6.0_WP*(v1+2.0_WP*v2+2.0_WP*v3+v4)
      ! For implicit RK2
      real(WP), dimension(3) :: p2old,v1
      real(WP) :: tolerance
      integer :: iter
      p2=p1
      tolerance=(1.0e-3_WP*this%cfg%min_meshsize)*(1.0e-3_WP*this%cfg%min_meshsize)
      do iter=1,10
         v1=this%cfg%get_velocity(0.5_WP*(p1+p2),i,j,k,U,V,W)
         p2old=p2
         p2=p1+dt*v1
         if (dot_product(p2-p2old,p2-p2old).lt.tolerance) exit
      end do
   end function project
   
   
   !> Calculate the min/max/int of our VF field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      integer :: ierr
      real(WP) :: my_VFmax,my_VFmin
      my_VFmax=maxval(this%VF); call MPI_ALLREDUCE(my_VFmax,this%VFmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_VFmin=minval(this%VF); call MPI_ALLREDUCE(my_VFmin,this%VFmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call this%cfg%integrate(this%VF,integral=this%VFint)
      call this%cfg%integrate(this%SD,integral=this%SDint)
   end subroutine get_max
   
   
   !> Write an IRL interface to a file
   subroutine write_interface(this,filename)
      use mpi_f08
      use messager, only: die
      use parallel, only: info_mpiio
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: filename
      logical :: file_is_there
      integer :: i,j,k,ind,ierr
      type(MPI_File) :: ifile
      type(MPI_Status) :: status
      integer(kind=MPI_OFFSET_KIND) :: disp
      integer :: size_to_write
      integer, dimension(3) :: dims
      integer,                        dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: number_of_planes
      integer(kind=MPI_OFFSET_KIND),  dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: offset_to_planes
      integer,                        dimension(this%cfg%ny_*this%cfg%nz ) :: array_of_block_lengths
      integer(kind=MPI_ADDRESS_KIND), dimension(this%cfg%ny_*this%cfg%nz_) :: array_of_displacements
      type(MPI_Datatype) :: MPI_OFFSET_ARRAY_TYPE
      type(ByteBuffer_type) :: byte_buffer
      
      ! Open the file
      inquire(file=trim(filename),exist=file_is_there)
      if (file_is_there.and.this%cfg%amRoot) call MPI_FILE_DELETE(trim(filename),info_mpiio,ierr)
      call MPI_FILE_OPEN(this%cfg%comm,trim(adjustl(filename)),IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[vfs interface write] Problem encountered while opening IRL data file: '//trim(filename))
      
      ! Write dimensions in header
      if (this%cfg%amRoot) then
         dims=[this%cfg%nx,this%cfg%ny,this%cfg%nz]
         call MPI_FILE_WRITE(ifile,dims,3,MPI_INTEGER,status,ierr)
      end if
      
      ! Calculate and store number of planes in each cell
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               number_of_planes(i,j,k)=getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
            end do
         end do
      end do
      
      ! Write out number of planes in each cell
      disp=int(4,8)*int(3,8) !< Only 3 int(4) - would need two more r(8) if we add time and dt
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_INTEGER,this%cfg%Iview,'native',info_mpiio,ierr)
      call MPI_FILE_WRITE_ALL(ifile,number_of_planes(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_INTEGER,status,ierr)
      
      ! Calculate the offset to each plane, needed for reading
      call this%calculate_offset_to_planes(number_of_planes,offset_to_planes)
      
      ! Make custom offset vector type for offsets
      ind=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            ind=ind+1
            array_of_block_lengths(ind)=int((offset_to_planes(this%cfg%imax_,j,k)-offset_to_planes(this%cfg%imin_,j,k)),4)+4+number_of_planes(this%cfg%imax_,j,k)*4*8+8
            array_of_displacements(ind)=int(offset_to_planes(this%cfg%imin_,j,k),MPI_ADDRESS_KIND)
         end do
      end do
      call MPI_TYPE_CREATE_HINDEXED(this%cfg%ny_*this%cfg%nz_,array_of_block_lengths,array_of_displacements,MPI_BYTE,MPI_OFFSET_ARRAY_TYPE,ierr)
      call MPI_TYPE_COMMIT(MPI_OFFSET_ARRAY_TYPE,ierr)
      
      ! Write out the actual PlanarSeps as packed bytes
      call new(byte_buffer)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call serializeAndPack(this%liquid_gas_interface(i,j,k),byte_buffer)
            end do
         end do
      end do
      disp=disp+int(4*this%cfg%nx*this%cfg%ny*this%cfg%nz,MPI_OFFSET_KIND)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_BYTE,MPI_OFFSET_ARRAY_TYPE,'native',info_mpiio,ierr)
      size_to_write=int(getSize(byte_buffer),4)
      call MPI_FILE_WRITE_ALL(ifile,dataPtr(byte_buffer),size_to_write,MPI_BYTE,status,ierr)
      
      ! Close file
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Free the type
      call MPI_TYPE_FREE(MPI_OFFSET_ARRAY_TYPE,ierr)
      
   end subroutine write_interface
   
   
   !> Read an IRL interface from a file
   subroutine read_interface(this,filename)
      use mpi_f08
      use messager, only: die
      use parallel, only: info_mpiio
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: i,j,k,ind,ierr
      type(MPI_File) :: ifile
      type(MPI_Status) :: status
      integer(kind=MPI_OFFSET_KIND) :: disp,size_to_read
      integer, dimension(3) :: dims
      integer,                        dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: number_of_planes
      integer(kind=MPI_OFFSET_KIND),  dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: offset_to_planes
      integer,                        dimension(this%cfg%ny_*this%cfg%nz ) :: array_of_block_lengths
      integer(kind=MPI_ADDRESS_KIND), dimension(this%cfg%ny_*this%cfg%nz_) :: array_of_displacements
      type(MPI_Datatype) :: MPI_OFFSET_ARRAY_TYPE
      type(ByteBuffer_type) :: byte_buffer
      
      ! Open the file
      call MPI_FILE_OPEN(this%cfg%comm,trim(adjustl(filename)),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[vfs interface read] Problem encountered while reading IRL data file: '//trim(filename))
      
      ! Read dimensions from header
      call MPI_FILE_READ_ALL(ifile,dims,4,MPI_INTEGER,status,ierr)
      
      ! Throw error if size mismatch
      if ((dims(1).ne.this%cfg%nx).or.(dims(2).ne.this%cfg%ny).or.(dims(3).ne.this%cfg%nz)) then
         if (this%cfg%amRoot) then
            print*, '    grid size = ',this%cfg%nx,this%cfg%ny,this%cfg%nz
            print*, 'IRL file size = ',dims(1),dims(2),dims(3)
         end if
         call die('[vfs interface read] The size of the interface file does not correspond to the grid')
      end if
      
      ! Read in number of planes
      call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_INTEGER,this%cfg%Iview,'native',info_mpiio,ierr)
      call MPI_FILE_READ_ALL(ifile,number_of_planes(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_INTEGER,status,ierr)
      
      ! Fill in ghost cells for number of planes
      call this%cfg%sync(number_of_planes)
      
      ! Calculate the offset to each plane, needed for reading
      call this%calculate_offset_to_planes(number_of_planes,offset_to_planes)
      
      ! Make custom offset vector type for offsets
      ind=0
      size_to_read=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            ind=ind+1
            array_of_block_lengths(ind)=int((offset_to_planes(this%cfg%imax_,j,k)-offset_to_planes(this%cfg%imin_,j,k)),4)+4+number_of_planes(this%cfg%imax_,j,k)*4*8+8
            size_to_read=size_to_read+int(array_of_block_lengths(ind),MPI_OFFSET_KIND)
            array_of_displacements(ind)=int(offset_to_planes(this%cfg%imin_,j,k),MPI_ADDRESS_KIND)
         end do
      end do
      call MPI_TYPE_CREATE_HINDEXED(this%cfg%ny_*this%cfg%nz_,array_of_block_lengths,array_of_displacements,MPI_BYTE,MPI_OFFSET_ARRAY_TYPE,ierr)
      call MPI_TYPE_COMMIT(MPI_OFFSET_ARRAY_TYPE,ierr)
      
      ! Read in the bytes and pack in to buffer, then loop through and unpack to PlanarSep
      call new(byte_buffer); call setSize(byte_buffer,size_to_read)
      if (size_to_read.ne.int(size_to_read,4)) call die('[vfs read interface] Cannot read that much data using the current I/O strategy') !< I/O WILL CRASH FOR IRL DATA >2Go/PROCESS
      disp=disp+int(4*this%cfg%nx*this%cfg%ny*this%cfg%nz,MPI_OFFSET_KIND)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_BYTE,MPI_OFFSET_ARRAY_TYPE,'native',info_mpiio,ierr)
      call MPI_FILE_READ_ALL(ifile,dataPtr(byte_buffer),int(size_to_read,4),MPI_BYTE,status,ierr)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call unpackAndStore(this%liquid_gas_interface(i,j,k),byte_buffer)
            end do
         end do
      end do
      
      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Free the type
      call MPI_TYPE_FREE(MPI_OFFSET_ARRAY_TYPE,ierr)
      
      ! Communicate interfaces
      call this%sync_interface()
      
   end subroutine read_interface
   
   
   !> Find byte offset for I/O of interface
   subroutine calculate_offset_to_planes(this,number_of_planes,offset_to_planes)
      use mpi_f08
      implicit none
      class(vfs), intent(in) :: this
      integer,                       dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: number_of_planes !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer(kind=MPI_OFFSET_KIND), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: offset_to_planes !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      integer(IRL_LargeOffsetIndex_t) :: nbytes
      integer :: isrc,idst,ierr
      type(MPI_Status) :: status
      
      ! Zero the offset to planes
      offset_to_planes=int(0,MPI_OFFSET_KIND)
      
      ! Calculate offsets in x direction
      isrc=this%cfg%xrank-1
      idst=this%cfg%xrank+1
      ! If leftmost processor, calculate offsets
      if (this%cfg%iproc.eq.1) then
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               nbytes=0
               do i=this%cfg%imin_,this%cfg%imax_+1
                  offset_to_planes(i,j,k)=nbytes
                  nbytes=nbytes+int(4,MPI_OFFSET_KIND)+int(number_of_planes(i,j,k),MPI_OFFSET_KIND)*int(4,MPI_OFFSET_KIND)*int(8,MPI_OFFSET_KIND)+int(8,MPI_OFFSET_KIND)
               end do
            end do
         end do
         if (this%cfg%iproc.ne.this%cfg%npx) call MPI_SEND(offset_to_planes(this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%ny_*this%cfg%nz_,MPI_INTEGER8,idst,0,this%cfg%xcomm,ierr)
      else
         ! Receive from the left processor
         call MPI_RECV(offset_to_planes(this%cfg%imin_-1,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%ny_*this%cfg%nz_,MPI_INTEGER8,isrc,0,this%cfg%xcomm,status,ierr)
         ! Calculate my offsets now that I know where my neighbor ended
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               nbytes=offset_to_planes(this%cfg%imin_-1,j,k)
               do i=this%cfg%imin_,this%cfg%imax_+1
                  offset_to_planes(i,j,k)=nbytes
                  nbytes=nbytes+int(4,MPI_OFFSET_KIND)+int(number_of_planes(i,j,k),MPI_OFFSET_KIND)*int(4,MPI_OFFSET_KIND)*int(8,MPI_OFFSET_KIND)+int(8,MPI_OFFSET_KIND)
               end do
            end do
         end do
         ! If not rightmost processor, send to next processor on the right
         if (this%cfg%iproc.ne.this%cfg%npx) call MPI_SEND(offset_to_planes(this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%ny_*this%cfg%nz_,MPI_INTEGER8,idst,0,this%cfg%xcomm,ierr)
      end if
      
      ! Calculate offsets in y direction
      if (this%cfg%iproc.eq.this%cfg%npx) then
         isrc=this%cfg%yrank-1
         idst=this%cfg%yrank+1
         ! If bottom processor, calculate offsets
         if (this%cfg%jproc.eq.1) then
            do k=this%cfg%kmin_,this%cfg%kmax_
               nbytes=0
               do j=this%cfg%jmin_,this%cfg%jmax_
                  offset_to_planes(this%cfg%imax_+1,j,k)=offset_to_planes(this%cfg%imax_+1,j,k)+nbytes
                  nbytes=offset_to_planes(this%cfg%imax_+1,j,k)
               end do
            end do
            if (this%cfg%jproc.ne.this%cfg%npy) call MPI_SEND(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%nz_,MPI_INTEGER8,idst,0,this%cfg%ycomm,ierr)
         else
            ! Receive from the bottom processor
            call MPI_RECV(offset_to_planes(this%cfg%imax_+1,this%cfg%jmin_-1,this%cfg%kmin_:this%cfg%kmax_),this%cfg%nz_,MPI_INTEGER8,isrc,0,this%cfg%ycomm,status,ierr)
            ! Calculate my offsets now that I know where my neighbor ended
            do k=this%cfg%kmin_,this%cfg%kmax_
               nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmin_-1,k)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  offset_to_planes(this%cfg%imax_+1,j,k)=offset_to_planes(this%cfg%imax_+1,j,k)+nbytes
                  nbytes=offset_to_planes(this%cfg%imax_+1,j,k)
               end do
            end do
            ! Send to the top processor, sent to next processor above
            if (this%cfg%jproc.ne.this%cfg%npy) call MPI_SEND(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%nz_,MPI_INTEGER8,idst,0,this%cfg%ycomm,ierr)
         end if
      end if
      
      ! Calculate offsets in z direction
      if (this%cfg%iproc.eq.this%cfg%npx.and.this%cfg%jproc.eq.this%cfg%npy) then
         isrc=this%cfg%zrank-1
         idst=this%cfg%zrank+1
         ! If first processor in z, calculate offsets
         if (this%cfg%kproc.eq.1) then
            nbytes=0
            do k=this%cfg%kmin_,this%cfg%kmax_
               offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)+nbytes
               nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)
            end do
            if (this%cfg%kproc.ne.this%cfg%npz) call MPI_SEND(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmax_),1,MPI_INTEGER8,idst,0,this%cfg%zcomm,ierr)
         else
            ! Receive from the previous processor
            call MPI_RECV(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_-1),1,MPI_INTEGER8,isrc,0,this%cfg%zcomm,status,ierr)
            ! Calculate my offsets now that I know where my neighbor ended
            nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_-1)
            do k=this%cfg%kmin_,this%cfg%kmax_
               offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)+nbytes
               nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)
            end do
            ! If not last processor in z, send to the next proc in z
            if (this%cfg%kproc.ne.this%cfg%npz) call MPI_SEND(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmax_),1,MPI_INTEGER8,idst,0,this%cfg%zcomm,ierr)
         end if
      end if
      
      ! Now need to unravel all this to update all of the offsets
      ! Be informed and add the j-offsets
      call MPI_BCAST(offset_to_planes(this%cfg%imax_+1,this%cfg%jmin_-1:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),(this%cfg%ny_+1)*this%cfg%nz_,MPI_INTEGER8,this%cfg%npx-1,this%cfg%xcomm,ierr)
      do k=this%cfg%kmin_,this%cfg%kmax_
         nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmin_-1,k)
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               offset_to_planes(i,j,k)=offset_to_planes(i,j,k)+nbytes
            end do
            nbytes=offset_to_planes(this%cfg%imax_+1,j,k)
         end do
      end do
      ! Be informed and add the z-offsets
      if (this%cfg%jproc.eq.this%cfg%npy) call MPI_BCAST(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_-1:this%cfg%kmax_),this%cfg%nz_+1,MPI_INTEGER8,this%cfg%npx-1,this%cfg%xcomm,ierr)
      call MPI_BCAST(offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_-1:this%cfg%kmax_),this%cfg%nz_+1,MPI_INTEGER8,this%cfg%npy-1,this%cfg%ycomm,ierr)
      nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,this%cfg%kmin_-1)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               offset_to_planes(i,j,k)=offset_to_planes(i,j,k)+nbytes
            end do
         end do
         nbytes=offset_to_planes(this%cfg%imax_+1,this%cfg%jmax_,k)
      end do
      
   end subroutine calculate_offset_to_planes
   
   
   !> Synchronize IRL objects across processors
   subroutine sync_interface(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ni
      real(WP), dimension(1:4) :: plane
      integer , dimension(2,3) :: send_range,recv_range
      ! Synchronize in x
      if (this%cfg%nx.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call copy(this%liquid_gas_interface(i,j,k),this%liquid_gas_interface(this%cfg%imin,j,k))
               end do
            end do
         end do
      else
         ! Send minus
         send_range(1:2,1)=[this%cfg%imin_   ,this%cfg%imin_ +this%cfg%no-1]
         send_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
         send_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
         recv_range(1:2,1)=[this%cfg%imax_ +1,this%cfg%imaxo_              ]
         recv_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
         recv_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
         call this%sync_side(send_range,recv_range,0,-1)
         ! Send plus
         send_range(1:2,1)=[this%cfg%imax_ -this%cfg%no+1,this%cfg%imax_   ]
         send_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
         send_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
         recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imin_ -1]
         recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
         recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
         call this%sync_side(send_range,recv_range,0,+1)
      end if
      ! Synchronize in y
      if (this%cfg%ny.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call copy(this%liquid_gas_interface(i,j,k),this%liquid_gas_interface(i,this%cfg%jmin,k))
               end do
            end do
         end do
      else
         ! Send minus side
         send_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
         send_range(1:2,2)=[this%cfg%jmin_   ,this%cfg%jmin_ +this%cfg%no-1]
         send_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
         recv_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
         recv_range(1:2,2)=[this%cfg%jmax_ +1,this%cfg%jmaxo_              ]
         recv_range(1:2,3)=[this%cfg%kmino_  ,this%cfg%kmaxo_              ]
         call this%sync_side(send_range,recv_range,1,-1)
         ! Send plus side
         send_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
         send_range(1:2,2)=[this%cfg%jmax_ -this%cfg%no+1,this%cfg%jmax_   ]
         send_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
         recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
         recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmin_ -1]
         recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmaxo_  ]
         call this%sync_side(send_range,recv_range,1,+1)
      end if
      ! Synchronize in z
      if (this%cfg%nz.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  call copy(this%liquid_gas_interface(i,j,k),this%liquid_gas_interface(i,j,this%cfg%kmin))
               end do
            end do
         end do
      else
         ! Send minus side
         send_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
         send_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
         send_range(1:2,3)=[this%cfg%kmin_   ,this%cfg%kmin_ +this%cfg%no-1]
         recv_range(1:2,1)=[this%cfg%imino_  ,this%cfg%imaxo_              ]
         recv_range(1:2,2)=[this%cfg%jmino_  ,this%cfg%jmaxo_              ]
         recv_range(1:2,3)=[this%cfg%kmax_ +1,this%cfg%kmaxo_              ]
         call this%sync_side(send_range,recv_range,2,-1)
         ! Send plus side
         send_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
         send_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
         send_range(1:2,3)=[this%cfg%kmax_ -this%cfg%no+1,this%cfg%kmax_   ]
         recv_range(1:2,1)=[this%cfg%imino_              ,this%cfg%imaxo_  ]
         recv_range(1:2,2)=[this%cfg%jmino_              ,this%cfg%jmaxo_  ]
         recv_range(1:2,3)=[this%cfg%kmino_              ,this%cfg%kmin_ -1]
         call this%sync_side(send_range,recv_range,2,+1)
      end if
      ! Fix plane position if we are periodic in x
      if (this%cfg%xper.and.this%cfg%iproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino,this%cfg%imin-1
                  do ni=0,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))-1
                     plane=getPlane(this%liquid_gas_interface(i,j,k),ni)
                     plane(4)=plane(4)-plane(1)*this%cfg%xL
                     call setPlane(this%liquid_gas_interface(i,j,k),ni,plane(1:3),plane(4))
                  end do
               end do
            end do
         end do
      end if
      if (this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imax+1,this%cfg%imaxo
                  do ni=0,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))-1
                     plane=getPlane(this%liquid_gas_interface(i,j,k),ni)
                     plane(4)=plane(4)+plane(1)*this%cfg%xL
                     call setPlane(this%liquid_gas_interface(i,j,k),ni,plane(1:3),plane(4))
                  end do
               end do
            end do
         end do
      end if
      ! Fix plane position if we are periodic in y
      if (this%cfg%yper.and.this%cfg%jproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino,this%cfg%jmin-1
               do i=this%cfg%imino_,this%cfg%imaxo_
                  do ni=0,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))-1
                     plane=getPlane(this%liquid_gas_interface(i,j,k),ni)
                     plane(4)=plane(4)-plane(2)*this%cfg%yL
                     call setPlane(this%liquid_gas_interface(i,j,k),ni,plane(1:3),plane(4))
                  end do
               end do
            end do
         end do
      end if
      if (this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmax+1,this%cfg%jmaxo
               do i=this%cfg%imino_,this%cfg%imaxo_
                  do ni=0,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))-1
                     plane=getPlane(this%liquid_gas_interface(i,j,k),ni)
                     plane(4)=plane(4)+plane(2)*this%cfg%yL
                     call setPlane(this%liquid_gas_interface(i,j,k),ni,plane(1:3),plane(4))
                  end do
               end do
            end do
         end do
      end if
      ! Fix plane position if we are periodic in z
      if (this%cfg%zper.and.this%cfg%kproc.eq.1) then
         do k=this%cfg%kmino,this%cfg%kmin-1
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  do ni=0,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))-1
                     plane=getPlane(this%liquid_gas_interface(i,j,k),ni)
                     plane(4)=plane(4)-plane(3)*this%cfg%zL
                     call setPlane(this%liquid_gas_interface(i,j,k),ni,plane(1:3),plane(4))
                  end do
               end do
            end do
         end do
      end if
      if (this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         do k=this%cfg%kmax+1,this%cfg%kmaxo
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  do ni=0,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))-1
                     plane=getPlane(this%liquid_gas_interface(i,j,k),ni)
                     plane(4)=plane(4)+plane(3)*this%cfg%zL
                     call setPlane(this%liquid_gas_interface(i,j,k),ni,plane(1:3),plane(4))
                  end do
               end do
            end do
         end do
      end if
   end subroutine sync_interface
   
   
   !> Private procedure to perform communication across one boundary
   subroutine sync_side(this,a_send_range,a_recv_range,a_dimension,a_direction)
      implicit none
      class(vfs), intent(inout) :: this
      integer, dimension(2,3), intent(in) :: a_send_range
      integer, dimension(2,3), intent(in) :: a_recv_range
      integer, intent(in) :: a_dimension
      integer, intent(in) :: a_direction
      integer :: i,j,k
      logical :: something_received
      ! Pack the buffer
      call resetBufferPointer(this%send_byte_buffer)
      call setSize(this%send_byte_buffer,int(0,8))
      do k=a_send_range(1,3),a_send_range(2,3)
         do j=a_send_range(1,2),a_send_range(2,2)
            do i=a_send_range(1,1),a_send_range(2,1)
               call serializeAndPack(this%liquid_gas_interface(i,j,k),this%send_byte_buffer)
            end do
         end do
      end do
      ! Communicate
      call this%sync_ByteBuffer(this%send_byte_buffer,a_dimension,a_direction,this%recv_byte_buffer,something_received)
      ! If something was received, unpack it: traversal order is important and must be aligned with how the sent data was packed
      if (something_received) then
         call resetBufferPointer(this%recv_byte_buffer)
         do k=a_recv_range(1,3),a_recv_range(2,3)
            do j=a_recv_range(1,2),a_recv_range(2,2)
               do i=a_recv_range(1,1),a_recv_range(2,1)
                  call unpackAndStore(this%liquid_gas_interface(i,j,k),this%recv_byte_buffer)
               end do
            end do
         end do
      end if
   end subroutine sync_side
   
   
   !> Private procedure to communicate a package of bytes across one boundary
   subroutine sync_ByteBuffer(this,a_send_buffer,a_dimension,a_direction,a_receive_buffer,a_received_something)
      use mpi_f08
      implicit none
      class(vfs), intent(inout) :: this
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
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,U,V,W,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFL
      
      ! Set the CFL to zero
      my_CFL=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFL=max(my_CFL,abs(U(i,j,k))*this%cfg%dxmi(i))
               my_CFL=max(my_CFL,abs(V(i,j,k))*this%cfg%dymi(j))
               my_CFL=max(my_CFL,abs(W(i,j,k))*this%cfg%dzmi(k))
            end do
         end do
      end do
      my_CFL=my_CFL*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFL,cfl,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_cfl
   
   
   !> Print out info for vf solver
   subroutine vfs_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(vfs), intent(in) :: this
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Volume fraction solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
   end subroutine vfs_print
   

end module vfs_class
