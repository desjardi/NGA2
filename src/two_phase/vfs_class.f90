!> Volume fraction solver class:
!> Provides support for various BC, semi-Lagrangian geometric advancement,
!> curvature calculation, interface reconstruction.
module vfs_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   use surfmesh_class, only: surfmesh
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
   !integer, parameter, public :: elvira=2            !< ELVIRA scheme
   !integer, parameter, public :: mof=3               !< MOF scheme
   !integer, parameter, public :: r2p=4               !< R2P scheme
   !integer, parameter, public :: swartz=5            !< Swartz scheme
   
   ! IRL cutting moment calculation method
   integer, parameter, public :: recursive_simplex=0 !< Recursive simplex cutting
   integer, parameter, public :: half_edge=1         !< Half-edge cutting
   integer, parameter, public :: nonrecurs_simplex=2 !< Non-recursive simplex cutting
   
   ! Default parameters for volume fraction solver
   integer,  parameter :: nband=3                                 !< Number of cells around the interfacial cells on which localized work is performed
   integer,  parameter :: advect_band=2                           !< How far we do the transport
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
      
      ! Subcell VF field
      real(WP), dimension(:,:,:,:,:,:), allocatable :: subvf   !< Subcell liquid volume fraction
      
      ! Surface density data
      real(WP), dimension(:,:,:), allocatable :: SD       !< Surface density array
      
      ! Distance level set
      real(WP) :: Gclip                                   !< Min/max distance
      real(WP), dimension(:,:,:), allocatable :: G        !< Distance level set array
      
      ! Curvature
      real(WP), dimension(:,:,:), allocatable :: curv     !< Interface mean curvature
      
      ! Band strategy
      integer, dimension(:,:,:), allocatable :: band      !< Band to localize workload around the interface
      integer, dimension(:,:),   allocatable :: band_map  !< Unstructured band mapping
      integer, dimension(nband) :: band_count             !< Number of cells per band value
      
      ! Interface reconstruction method
      integer :: reconstruction_method                    !< Interface reconstruction method
      
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
      type(SepVM_type),       dimension(:,:,:,:), allocatable :: face_flux
      
      ! More basic representation of the surface grid
      type(surfmesh) :: surfgrid
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for enforcing bconds
      
      ! Monitoring quantities
      real(WP) :: VFmax,VFmin,VFint                       !< Maximum, minimum, and integral volume fraction
      
   contains
      procedure :: print=>vfs_print                       !< Output solver to the screen
      procedure :: initialize_irl                         !< Initialize the IRL objects
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: update_band                            !< Update the band info given the VF values
      procedure :: sync_interface                         !< Synchronize the IRL objects
      procedure :: sync_and_clean_barycenters             !< Synchronize and clean up phasic barycenters
      procedure, private :: sync_side                     !< Synchronize the IRL objects across one side - another I/O helper
      procedure, private :: sync_ByteBuffer               !< Communicate byte packets across one side - another I/O helper
      procedure, private :: calculate_offset_to_planes    !< Helper routine for I/O
      procedure, private :: crude_phase_test              !< Helper function that rapidly assess if a mixed cell might be present
      procedure, private :: project                       !< Helper function that performs a Lagrangian projection of a vertex
      procedure, private :: get_velocity                  !< Helper function that interpolates a velocity field to a point
      procedure :: read_interface                         !< Read an IRL interface from a file
      procedure :: write_interface                        !< Write an IRL interface to a file
      procedure :: advance                                !< Advance VF to next step
      procedure :: advect_interface                       !< Advance IRL surface to next step
      procedure :: build_interface                        !< Reconstruct IRL interface from VF field
      procedure :: build_lvira                            !< LVIRA reconstruction of the interface from VF field
      procedure :: polygonalize_interface                 !< Build a discontinuous polygonal representation of the IRL interface
      procedure :: distance_from_polygon                  !< Build a signed distance field from the polygonalized interface
      procedure :: subcell_vf                             !< Build subcell volume fraction from reconstructed interface
      procedure :: reset_volume_moments                   !< Reconstruct volume moments from IRL interfaces
      procedure :: update_surfgrid                        !< Create a simple surface mesh from the IRL polygons
      procedure :: get_curvature                          !< Compute curvature from IRL surface polygons
      procedure :: paraboloid_fit                         !< Perform local paraboloid fit of IRL surface
      procedure :: get_max                                !< Calculate maximum field values
      procedure :: get_cfl                                !< Get CFL for the VF solver
   end type vfs
   
   
   !> Declare volume fraction solver constructor
   interface vfs
      procedure constructor
   end interface vfs
   
contains
   
   
   !> Default constructor for volume fraction solver
   function constructor(cfg,reconstruction_method,name) result(self)
      use messager, only: die
      implicit none
      type(vfs) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: reconstruction_method
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Check that we have at least 2 overlap cells
      if (cfg%no.lt.2) call die('[vfs constructor] The config requires at least 2 overlap cells')
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate variables
      allocate(self%VF   (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF   =0.0_WP
      allocate(self%VFold(  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFold=0.0_WP
      allocate(self%Lbary(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Lbary=0.0_WP
      allocate(self%Gbary(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Gbary=0.0_WP
      allocate(self%SD   (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SD   =0.0_WP
      allocate(self%G    (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%G    =0.0_WP
      allocate(self%curv (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%curv =0.0_WP
      
      ! Set clipping distance
      self%Gclip=real(distance_band+1,WP)*self%cfg%min_meshsize
      
      ! Subcell volume fraction
      allocate(self%subvf(0:1,0:1,0:1,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%subvf=0.0_WP
      
      ! Prepare the band arrays
      allocate(self%band(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%band=0
      if (allocated(self%band_map)) deallocate(self%band_map)
      
      ! Set reconstruction method
      self%reconstruction_method=reconstruction_method
      
      ! Initialize IRL
      call self%initialize_irl()
      
      ! Also initialize a surface grid
      self%surfgrid=surfmesh(name='plic')
      
      ! Prepare mask for VF
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
      
   end function constructor
   
   
   !> Initialize the IRL representation of our interfaces
   subroutine initialize_irl(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,n,tag
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
      
      ! Work arrays for face fluxes
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
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,dir,locator)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      character(len=2), intent(in) :: dir
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      select case (lowercase(dir))
      case ('+x','x+','xp','px'); new_bc%dir=1
      case ('-x','x-','xm','mx'); new_bc%dir=2
      case ('+y','y+','yp','py'); new_bc%dir=3
      case ('-y','y-','ym','my'); new_bc%dir=4
      case ('+z','z+','zp','pz'); new_bc%dir=5
      case ('-z','z-','zm','mz'); new_bc%dir=6
      case default; call die('[vfs add_bcond] Unknown bcond direction')
      end select
      new_bc%itr=iterator(this%cfg,new_bc%name,locator)
      
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
         ! Not yet implemented
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
      
      ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! X flux
               if (minval(abs(this%band(i-1:i,j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W)
                  face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_indices(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_indices(bounding_pts(:,2),[i,j,k])
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
                  face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_indices(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_indices(bounding_pts(:,2),[i,j,k])
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
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W)
                  face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_indices(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_indices(bounding_pts(:,2),[i,j,k])
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
      do index=1,sum(this%band_count(1:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
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
      
      ! Synchronize VF field
      call this%cfg%sync(this%VF)
      
      ! Synchronize and clean-up barycenter fields
      call this%sync_and_clean_barycenters()
      
      ! Advect interface polygons
      call this%advect_interface(dt,U,V,W)
      
      ! Update the band
      call this%update_band()
      
      ! Perform interface reconstruction from transported moments
      call this%build_interface()
      
      ! Create discontinuous polygon mesh from IRL interface
      call this%polygonalize_interface()
      
      ! Calculate distance from polygons
      call this%distance_from_polygon()
      
      ! Calculate subcell volume fraction
      call this%subcell_vf()
      
      ! Calculate curvature
      call this%get_curvature()
      
      ! Reset moments to guarantee compatibility with interface reconstruction
      call this%reset_volume_moments()
      
   end subroutine advance
   
   
   !> Synchronize and clean up barycenter fields
   subroutine sync_and_clean_barycenters(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      ! Clean up barycenters everywhere
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%VF(i,j,k).lt.VFlo.or.this%VF(i,j,k).gt.VFhi) then
                  this%Lbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  this%Gbary(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
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
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
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
      
   end subroutine advect_interface
   
   
   !> Band update from VF dataset
   subroutine update_band(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,dir,n,index
      integer, dimension(3) :: ind
      integer :: ibmin_,ibmax_,jbmin_,jbmax_,kbmin_,kbmax_
      integer, dimension(nband) :: band_size
      
      ! Loop extents
      ibmin_=this%cfg%imin_; if (this%cfg%iproc.eq.1           .and.this%cfg%xper) ibmin_=this%cfg%imin-1
      ibmax_=this%cfg%imax_; if (this%cfg%iproc.eq.this%cfg%npx.and.this%cfg%xper) ibmax_=this%cfg%imax+1
      jbmin_=this%cfg%jmin_; if (this%cfg%jproc.eq.1           .and.this%cfg%yper) jbmin_=this%cfg%jmin-1
      jbmax_=this%cfg%jmax_; if (this%cfg%jproc.eq.this%cfg%npy.and.this%cfg%yper) jbmax_=this%cfg%jmax+1
      kbmin_=this%cfg%kmin_; if (this%cfg%kproc.eq.1           .and.this%cfg%zper) kbmin_=this%cfg%kmin-1
      kbmax_=this%cfg%kmax_; if (this%cfg%kproc.eq.this%cfg%npz.and.this%cfg%zper) kbmax_=this%cfg%kmax+1
      
      ! Reset band
      this%band=(nband+1)*int(sign(1.0_WP,this%VF-0.5_WP))
      
      ! First sweep to identify cells with interface
      do k=kbmin_,kbmax_
         do j=jbmin_,jbmax_
            do i=ibmin_,ibmax_
               ! Skip *real* wall cells
               if (this%mask(i,j,k).eq.1) cycle
               ! Identify cells with interface
               if (this%VF(i,j,k).ge.VFlo.and.this%VF(i,j,k).le.VFhi) this%band(i,j,k)=int(sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
               ! Check if cell-face is an interface
               do dir=1,3
                  do n=-1,+1,2
                     ind=[i,j,k]; ind(dir)=ind(dir)+n
                     if (this%mask(ind(1),ind(2),ind(3)).ne.1) then
                        if (this%VF(i,j,k).lt.VFlo.and.this%VF(ind(1),ind(2),ind(3)).gt.VFhi.or.&
                        &   this%VF(i,j,k).gt.VFhi.and.this%VF(ind(1),ind(2),ind(3)).lt.VFlo) this%band(i,j,k)=int(sign(1.0_WP,this%VF(i,j,k)-0.5_WP))
                     end if
                  end do
               end do
            end do
         end do
      end do
      call this%cfg%sync(this%band)
      
      ! Sweep to identify the bands up to nband
      do n=2,nband
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
               if (this%band(i,j,k).le.nband.and.this%band(i,j,k).ne.-nband-1) band_size(abs(this%band(i,j,k)))=band_size(abs(this%band(i,j,k)))+1
            end do
         end do
      end do
      
      ! Rebuild the unstructured mapping
      if (allocated(this%band_map)) deallocate(this%band_map); allocate(this%band_map(3,sum(band_size)))
      this%band_count=0
      do k=kbmin_,kbmax_
         do j=jbmin_,jbmax_
            do i=ibmin_,ibmax_
               if (this%band(i,j,k).le.nband.and.this%band(i,j,k).ne.-nband-1) then
                  this%band_count(abs(this%band(i,j,k)))=this%band_count(abs(this%band(i,j,k)))+1
                  index=sum(band_size(1:abs(this%band(i,j,k))-1))+this%band_count(abs(this%band(i,j,k)))
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
      case (lvira); call this%build_lvira()
      case default; call die('[vfs build interface] Unknown interface reconstruction scheme')
      end select
   end subroutine build_interface
   
   
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
      
      ! Finally, update the basic unstructured surface mesh representation of our polygons
      call this%update_surfgrid()
      
   end subroutine polygonalize_interface
   
   
   !> Create a basic surface mesh from our IRL polygons
   subroutine update_surfgrid(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k,n,shape,nv,np,nplane
      real(WP), dimension(3) :: tmp_vert
      
      ! Reset surface mesh storage
      call this%surfgrid%reset()
      
      ! First pass to count how many vertices and polygones are inside our processor
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
         call this%surfgrid%set_size(nvert=nv,npoly=np)
         allocate(this%surfgrid%polyConn(this%surfgrid%nVert)) ! Also allocate naive connectivity
         nv=0; np=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(this%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        this%surfgrid%polySize(np)=shape
                        ! Loop over its vertices and add them
                        do n=1,shape
                           tmp_vert=getPt(this%interface_polygon(nplane,i,j,k),n-1)
                           ! Increment node counter
                           nv=nv+1
                           this%surfgrid%xVert(nv)=tmp_vert(1)
                           this%surfgrid%yVert(nv)=tmp_vert(2)
                           this%surfgrid%zVert(nv)=tmp_vert(3)
                           this%surfgrid%polyConn(nv)=nv
                        end do
                     end if
                  end do
               end do
            end do
         end do
      !else
      !   ! Add a zero-area triangle if this proc doesn't have one
      !   np=1; nv=3
      !   call this%surfgrid%set_size(nvert=nv,npoly=np)
      !   allocate(this%surfgrid%polyConn(this%surfgrid%nVert)) ! Also allocate naive connectivity
      !   this%surfgrid%xVert(1:3)=this%cfg%x(this%cfg%imin)
      !   this%surfgrid%yVert(1:3)=this%cfg%y(this%cfg%jmin)
      !   this%surfgrid%zVert(1:3)=this%cfg%z(this%cfg%kmin)
      !   this%surfgrid%polySize(1)=3
      !   this%surfgrid%polyConn(1:3)=[1,2,3]
      end if
      
   end subroutine update_surfgrid
   
   
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
      do index=1,sum(this%band_count(1:distance_band))
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
   
   
   !> Calculate subcell volume fraction from reconstructed interface
   !> Here, only mask=1 is skipped (i.e., real walls), so bconds are handled
   subroutine subcell_vf(this)
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
                  this%subvf(:,:,:,i,j,k)=0.0_WP
                  cycle
               end if
               ! Deal with other cells
               if (this%VF(i,j,k).gt.VFhi) then
                  this%subvf(:,:,:,i,j,k)=1.0_WP
               else if (this%VF(i,j,k).lt.VFlo) then
                  this%subvf(:,:,:,i,j,k)=0.0_WP
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
                           this%subvf(ii,jj,kk,i,j,k)=8.0_WP*getVolume(separated_volume_moments,0)/this%cfg%vol(i,j,k)
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      
   end subroutine subcell_vf
   
   
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
   
   
   !> Compute curvature from a least squares fit of the IRL surface
   subroutine get_curvature(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer :: i,j,k
      real(WP) :: mycurv
      ! Reset curvature
      this%curv=0.0_WP
      ! Traverse interior domain and compute curvature in cells with polygons
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Skip cells without polygon
               if (getNumberOfVertices(this%interface_polygon(1,i,j,k)).eq.0) cycle
               ! Perform LSQ PLIC barycenter fitting
               call this%paraboloid_fit(i,j,k,mycurv)
               ! Store clipped curvature
               this%curv(i,j,k)=max(min(mycurv,1.0_WP/this%cfg%meshsize(i,j,k)),-1.0_WP/this%cfg%meshsize(i,j,k))
            end do
         end do
      end do
      ! Synchronize boundaries
      call this%cfg%sync(this%curv)
   end subroutine get_curvature
   
   
   !> Perform local paraboloid fit of IRL surface
   subroutine paraboloid_fit(this,i,j,k,mycurv)
      use mathtools, only: normalize,cross_product
      implicit none
      ! In/out variables
      class(vfs), intent(inout) :: this
      integer,  intent(in)  :: i,j,k
      real(WP), intent(out) :: mycurv
      ! Variables used to process the polygonal surface
      real(WP), dimension(3) :: pref,nref,tref,sref
      real(WP), dimension(3) :: ploc,nloc
      real(WP), dimension(3) :: buf
      real(WP) :: surf,ww
      integer :: ii,jj,kk,ndata,info
      ! Storage for least squares problem
      real(WP), dimension(125,6) :: A=0.0_WP
      real(WP), dimension(125)   :: b=0.0_WP
      real(WP), dimension(6)     :: sol
      real(WP), dimension(200)   :: work
      ! Curvature evaluation
      real(WP) :: dF_dt,dF_ds,ddF_dtdt,ddF_dsds,ddF_dtds
      
      ! Store polygon centroid - this is our reference point
      pref=calculateCentroid(this%interface_polygon(1,i,j,k))
      
      ! Create local basis from polygon normal
      nref=calculateNormal(this%interface_polygon(1,i,j,k))
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
               
               ! Skip the cell if there is no polygon in it
               if (getNumberOfVertices(this%interface_polygon(1,ii,jj,kk)).eq.0) cycle
               
               ! Skip the cell if it's a true wall
               if (this%mask(ii,jj,kk).eq.1) cycle
               
               ! Get local polygon normal
               nloc=calculateNormal(this%interface_polygon(1,ii,jj,kk))
               
               ! Store triangle centroid, and surface
               ploc=    calculateCentroid(this%interface_polygon(1,ii,jj,kk))
               surf=abs(calculateVolume  (this%interface_polygon(1,ii,jj,kk)))/this%cfg%meshsize(i,j,k)**2
               
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
      
      ! Solve for surface as n=F(t,s)=b1+b2*t+b3*s+b4*t^2+b5*s^2+b6*t*s using Lapack
      call dgels('N',ndata,6,1,A,125,b,125,work,200,info); sol=b(1:6)
      
      ! Get the curvature at (t,s)=(0,0)
      dF_dt=sol(2)+sol(4)*0.0_WP+sol(6)*0.0_WP; ddF_dtdt=sol(4)
      dF_ds=sol(3)+sol(5)*0.0_WP+sol(6)*0.0_WP; ddF_dsds=sol(5)
      ddF_dtds=sol(6)
      mycurv=-((1+dF_dt**2)*ddF_dsds-2.0_WP*dF_dt*dF_ds*ddF_dtds+(1.0_WP+dF_ds**2)*ddF_dtdt)/(1.0_WP+dF_dt**2+dF_ds**2)**(1.5_WP)
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
      !v1=this%get_velocity(p1             ,i,j,k,U,V,W)
      !v2=this%get_velocity(p1+0.5_WP*dt*v1,i,j,k,U,V,W)
      !v3=this%get_velocity(p1+0.5_WP*dt*v2,i,j,k,U,V,W)
      !v4=this%get_velocity(p1+       dt*v3,i,j,k,U,V,W)
      !p2=p1+dt/6.0_WP*(v1+2.0_WP*v2+2.0_WP*v3+v4)
      ! For implicit RK2
      real(WP), dimension(3) :: p2old,v1
      real(WP) :: tolerance
      integer :: iter
      p2=p1
      tolerance=(1.0e-3_WP*this%cfg%min_meshsize)*(1.0e-3_WP*this%cfg%min_meshsize)
      do iter=1,10
         v1=this%get_velocity(0.5_WP*(p1+p2),i,j,k,U,V,W)
         p2old=p2
         p2=p1+dt*v1
         if (dot_product(p2-p2old,p2-p2old).lt.tolerance) exit
      end do
   end function project
   
   
   !> Private function that performs an trilinear interpolation of the provided velocity U,V,W
   !> to the provided position pos in the vicinity of cell i0,j0,k0
   function get_velocity(this,pos,i0,j0,k0,U,V,W) result(vel)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), dimension(3) :: vel
      real(WP), dimension(3), intent(in) :: pos
      integer, intent(in) :: i0,j0,k0
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP) :: wx1,wy1,wz1
      real(WP) :: wx2,wy2,wz2
      ! Interpolate U velocity ------------------------------
      ! Find right i index
      i=max(min(this%cfg%imaxo_-1,i0),this%cfg%imino_)
      do while (pos(1)-this%cfg%x (i  ).lt.0.0_WP.and.i  .gt.this%cfg%imino_); i=i-1; end do
      do while (pos(1)-this%cfg%x (i+1).ge.0.0_WP.and.i+1.lt.this%cfg%imaxo_); i=i+1; end do
      ! Find right j index
      j=max(min(this%cfg%jmaxo_-1,j0),this%cfg%jmino_)
      do while (pos(2)-this%cfg%ym(j  ).lt.0.0_WP.and.j  .gt.this%cfg%jmino_); j=j-1; end do
      do while (pos(2)-this%cfg%ym(j+1).ge.0.0_WP.and.j+1.lt.this%cfg%jmaxo_); j=j+1; end do
      ! Find right k index
      k=max(min(this%cfg%kmaxo_-1,k0),this%cfg%kmino_)
      do while (pos(3)-this%cfg%zm(k  ).lt.0.0_WP.and.k  .gt.this%cfg%kmino_); k=k-1; end do
      do while (pos(3)-this%cfg%zm(k+1).ge.0.0_WP.and.k+1.lt.this%cfg%kmaxo_); k=k+1; end do
      ! Prepare tri-linear interpolation coefficients
      wx1=(pos(1)-this%cfg%x (i))/(this%cfg%x (i+1)-this%cfg%x (i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%cfg%ym(j))/(this%cfg%ym(j+1)-this%cfg%ym(j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%cfg%zm(k))/(this%cfg%zm(k+1)-this%cfg%zm(k)); wz2=1.0_WP-wz1
      ! Tri-linear interpolation of U
      vel(1)=wz1*(wy1*(wx1*U(i+1,j+1,k+1)  + &
      &                wx2*U(i  ,j+1,k+1)) + &
      &           wy2*(wx1*U(i+1,j  ,k+1)  + &
      &                wx2*U(i  ,j  ,k+1)))+ &
      &      wz2*(wy1*(wx1*U(i+1,j+1,k  )  + &
      &                wx2*U(i  ,j+1,k  )) + &
      &           wy2*(wx1*U(i+1,j  ,k  )  + &
      &                wx2*U(i  ,j  ,k  )))
      ! Interpolate V velocity ------------------------------
      ! Find right i index
      i=max(min(this%cfg%imaxo_-1,i0),this%cfg%imino_)
      do while (pos(1)-this%cfg%xm(i  ).lt.0.0_WP.and.i  .gt.this%cfg%imino_); i=i-1; end do
      do while (pos(1)-this%cfg%xm(i+1).ge.0.0_WP.and.i+1.lt.this%cfg%imaxo_); i=i+1; end do
      ! Find right j index
      j=max(min(this%cfg%jmaxo_-1,j0),this%cfg%jmino_)
      do while (pos(2)-this%cfg%y (j  ).lt.0.0_WP.and.j  .gt.this%cfg%jmino_); j=j-1; end do
      do while (pos(2)-this%cfg%y (j+1).ge.0.0_WP.and.j+1.lt.this%cfg%jmaxo_); j=j+1; end do
      ! Find right k index
      k=max(min(this%cfg%kmaxo_-1,k0),this%cfg%kmino_)
      do while (pos(3)-this%cfg%zm(k  ).lt.0.0_WP.and.k  .gt.this%cfg%kmino_); k=k-1; end do
      do while (pos(3)-this%cfg%zm(k+1).ge.0.0_WP.and.k+1.lt.this%cfg%kmaxo_); k=k+1; end do
      ! Prepare tri-linear interpolation coefficients
      wx1=(pos(1)-this%cfg%xm(i))/(this%cfg%xm(i+1)-this%cfg%xm(i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%cfg%y (j))/(this%cfg%y (j+1)-this%cfg%y (j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%cfg%zm(k))/(this%cfg%zm(k+1)-this%cfg%zm(k)); wz2=1.0_WP-wz1
      ! Tri-linear interpolation of V
      vel(2)=wz1*(wy1*(wx1*V(i+1,j+1,k+1)  + &
      &                wx2*V(i  ,j+1,k+1)) + &
      &           wy2*(wx1*V(i+1,j  ,k+1)  + &
      &                wx2*V(i  ,j  ,k+1)))+ &
      &      wz2*(wy1*(wx1*V(i+1,j+1,k  )  + &
      &                wx2*V(i  ,j+1,k  )) + &
      &           wy2*(wx1*V(i+1,j  ,k  )  + &
      &                wx2*V(i  ,j  ,k  )))
      ! Interpolate W velocity ------------------------------
      ! Find right i index
      i=max(min(this%cfg%imaxo_-1,i0),this%cfg%imino_)
      do while (pos(1)-this%cfg%xm(i  ).lt.0.0_WP.and.i  .gt.this%cfg%imino_); i=i-1; end do
      do while (pos(1)-this%cfg%xm(i+1).ge.0.0_WP.and.i+1.lt.this%cfg%imaxo_); i=i+1; end do
      ! Find right j index
      j=max(min(this%cfg%jmaxo_-1,j0),this%cfg%jmino_)
      do while (pos(2)-this%cfg%ym(j  ).lt.0.0_WP.and.j  .gt.this%cfg%jmino_); j=j-1; end do
      do while (pos(2)-this%cfg%ym(j+1).ge.0.0_WP.and.j+1.lt.this%cfg%jmaxo_); j=j+1; end do
      ! Find right k index
      k=max(min(this%cfg%kmaxo_-1,k0),this%cfg%kmino_)
      do while (pos(3)-this%cfg%z (k  ).lt.0.0_WP.and.k  .gt.this%cfg%kmino_); k=k-1; end do
      do while (pos(3)-this%cfg%z (k+1).ge.0.0_WP.and.k+1.lt.this%cfg%kmaxo_); k=k+1; end do
      ! Prepare tri-linear interpolation coefficients
      wx1=(pos(1)-this%cfg%xm(i))/(this%cfg%xm(i+1)-this%cfg%xm(i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%cfg%ym(j))/(this%cfg%ym(j+1)-this%cfg%ym(j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%cfg%z (k))/(this%cfg%z (k+1)-this%cfg%z (k)); wz2=1.0_WP-wz1
      ! Tri-linear interpolation of W
      vel(3)=wz1*(wy1*(wx1*W(i+1,j+1,k+1)  + &
      &                wx2*W(i  ,j+1,k+1)) + &
      &           wy2*(wx1*W(i+1,j  ,k+1)  + &
      &                wx2*W(i  ,j  ,k+1)))+ &
      &      wz2*(wy1*(wx1*W(i+1,j+1,k  )  + &
      &                wx2*W(i  ,j+1,k  )) + &
      &           wy2*(wx1*W(i+1,j  ,k  )  + &
      &                wx2*W(i  ,j  ,k  )))
      return
   end function get_velocity
   
   
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
      integer                       , dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: number_of_planes
      integer(kind=MPI_OFFSET_KIND) , dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: offset_to_planes
      integer                       , dimension(this%cfg%ny_*this%cfg%nz) :: array_of_block_lengths
      integer(kind=MPI_ADDRESS_KIND), dimension(this%cfg%ny_*this%cfg%nz_):: array_of_displacements
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
      disp=int(4,8)*int(3,8) !< Only 3 int(4) - would need to more r(8) if we add time and dt
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_INTEGER,this%cfg%Iview,'native',info_mpiio,ierr)
      call MPI_FILE_WRITE_ALL(ifile,number_of_planes(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_INTEGER,status,ierr)
      
      ! Calculate the offset to each plane, needed for reading
      call this%calculate_offset_to_planes(number_of_planes,offset_to_planes)
      
      ! Make custom offset vector type for offsets
      ind=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            ind=ind+1
            array_of_block_lengths(ind)=int((offset_to_planes(this%cfg%imax_,j,k)-offset_to_planes(this%cfg%imin_,j,k)),4)+int(4,4)+int(number_of_planes(this%cfg%imax_,j,k),4)*int(4,4)*int(8,4)+int(8,4)
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
      disp=disp+int(4,MPI_OFFSET_KIND)*int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)
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
      integer(kind=MPI_OFFSET_KIND) :: disp,size_to_read_big
      integer :: size_to_read
      integer, dimension(3) :: dims
      integer                       , dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: number_of_planes
      integer(kind=MPI_OFFSET_KIND) , dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: offset_to_planes
      integer                       , dimension(this%cfg%ny_*this%cfg%nz) :: array_of_block_lengths
      integer(kind=MPI_ADDRESS_KIND), dimension(this%cfg%ny_*this%cfg%nz_):: array_of_displacements
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
      size_to_read_big=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            ind=ind+1
            array_of_block_lengths(ind)=int((offset_to_planes(this%cfg%imax_,j,k)-offset_to_planes(this%cfg%imin_,j,k)),4)+int(4,4)+int(number_of_planes(this%cfg%imax_,j,k),4)*int(4,4)*int(8,4)+int(8,4)
            size_to_read_big=size_to_read_big+int(array_of_block_lengths(ind),MPI_OFFSET_KIND)
            array_of_displacements(ind)=int(offset_to_planes(this%cfg%imin_,j,k),MPI_ADDRESS_KIND)
         end do
      end do
      call MPI_TYPE_CREATE_HINDEXED(this%cfg%ny_*this%cfg%nz_,array_of_block_lengths,array_of_displacements,MPI_BYTE,MPI_OFFSET_ARRAY_TYPE,ierr)
      call MPI_TYPE_COMMIT(MPI_OFFSET_ARRAY_TYPE,ierr)
      
      ! Read in the bytes and pack in to buffer, then loop through and unpack to PlanarSep
      disp=disp+int(4,MPI_OFFSET_KIND)*int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)
      call new(byte_buffer)
      call setSize(byte_buffer,size_to_read_big)
      if (size_to_read_big.gt.huge(size_to_read)) call die('[vfs read interface] Cannot read that much data using the current I/O strategy') !< I/O WILL CRASH FOR IRL DATA >2Go/PROCESS
      size_to_read=int(size_to_read_big)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_BYTE,MPI_OFFSET_ARRAY_TYPE,'native',info_mpiio,ierr)
      call MPI_FILE_READ_ALL(ifile,dataPtr(byte_buffer),size_to_read,MPI_BYTE,status,ierr)
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
      offset_to_planes = int(0,MPI_OFFSET_KIND)
      
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
      ! Fix plane posistion if we are periodic in x
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
      ! Fix plane posistion if we are periodic in z
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