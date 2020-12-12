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
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available interface reconstructions schemes for VF
   integer, parameter, public :: lvira=1             !< LVIRA scheme
   !integer, parameter, public :: elvira=2            !< ELVIRA scheme
   !integer, parameter, public :: mof=3               !< MOF scheme
   !integer, parameter, public :: r2p=4               !< R2P scheme
   !integer, parameter, public :: swartz=5            !< Swartz scheme
   
   
   ! Default parameters for volume fraction solver
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
      procedure :: sync_interface                         !< Synchronize the IRL objects
      procedure, private :: sync_side                     !< Synchronize the IRL objects across one side - another I/O helper
      procedure, private :: sync_ByteBuffer               !< Communicate byte packets across one side - another I/O helper
      procedure, private :: calculate_offset_to_planes    !< Helper routine for I/O
      procedure :: read_interface                         !< Read an IRL interface from a file
      procedure :: write_interface                        !< Write an IRL interface to a file
      procedure :: advance                                !< Advance VF to next step
      procedure :: build_interface                        !< Reconstruct IRL interface from VF field
      procedure :: build_lvira                            !< LVIRA reconstruction of the interface from VF field
      !procedure :: reset_moments                          !< Reconstruct volume moments from IRL interfaces
      procedure :: get_max                                !< Calculate maximum field values
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
      
      ! Set reconstruction method
      self%reconstruction_method=reconstruction_method
      
      ! Initialize IRL
      call self%initialize_irl()
      
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
               tag=(i-this%cfg%imino_)+(j-this%cfg%jmino_)*this%cfg%nxo_+(k-this%cfg%kmino_)*this%cfg%nxo_*this%cfg%nyo_
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
      
      ! Create polygon mesh based on IRL interface
      !call multiphase_plic_poly
      
      ! Calculate distance from polygons field
      !call multiphase_plic_distance
      
      ! Calculate subcell volume fractions
      !call multiphase_plic_subcell
      
      ! Calculate curvature
      !call multiphase_curvature_calc
      
   end subroutine build_interface
   
   
   !> LVIRA reconstruction of a planar interface in mixed cells
   subroutine build_lvira(this)
      implicit none
      class(vfs), intent(inout) :: this
      integer(IRL_SignedIndex_t) :: i,j,k
      integer :: ind,ii,jj,kk,icenter
      type(LVIRANeigh_RectCub_type) :: neighborhood
      type(RectCub_type), dimension(0:26) :: neighborhood_cells
      real(IRL_double)  , dimension(0:26) :: liquid_volume_fraction
      real(IRL_double), dimension(3) :: initial_normal
      real(IRL_double) :: initial_distance
      
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
                        liquid_volume_fraction(ind)=this%cfg%VF(ii,jj,kk)
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
               initial_normal   = this%Gbary(:,i,j,k) - this%Lbary(:,i,j,k)
               initial_distance = sqrt(sum(initial_normal**2.0_WP))
               initial_normal   = initial_normal/(initial_distance+tiny(1.0_WP))
               initial_distance = dot_product(initial_normal,[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)])
               call setPlane(this%liquid_gas_interface(i,j,k),0,initial_normal,initial_distance)
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
   
   
   !> Write an IRL interface to a file
   subroutine write_interface(this,filename)
      use mpi_f08
      use messager, only: die
      use parallel, only: info_mpiio,MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: filename
      logical :: file_is_there
      integer :: i,j,k,ind,ierr,var
      integer, dimension(3) :: gsizes,lsizes,start
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
      use parallel, only: info_mpiio,MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: i,j,k,ind,ierr,var
      integer, dimension(3) :: gsizes,lsizes,start
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
      ! Fix plane posistion if we are periodic in y
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
      type(ByteBuffer_type), intent(in)  :: a_send_buffer
      integer, intent(in) :: a_dimension  !< Should be 0/1/2 for x/y/z
      integer, intent(in) :: a_direction  !< Should be -1 for left or +1 for right
      type(ByteBuffer_type), intent(out) :: a_receive_buffer
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
   
   
   !> Calculate the new VF based on U/V/W
   subroutine advance(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(in) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Flux of VF
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
      !         FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%itpsc_xp(:,i,j,k)*this%SC(i+this%stp1:i+this%stp2,j,k)) &
      !         &         -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%itpsc_xm(:,i,j,k)*this%SC(i+this%stm1:i+this%stm2,j,k)) &
      !         &         +this%diff*sum(this%grdsc_x(:,i,j,k)*this%SC(i-1:i,j,k))
               ! Fluxes on y-face
      !         FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%itpsc_yp(:,i,j,k)*this%SC(i,j+this%stp1:j+this%stp2,k)) &
      !         &         -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%itpsc_ym(:,i,j,k)*this%SC(i,j+this%stm1:j+this%stm2,k)) &
      !         &         +this%diff*sum(this%grdsc_y(:,i,j,k)*this%SC(i,j-1:j,k))
               ! Fluxes on z-face
      !         FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%itpsc_zp(:,i,j,k)*this%SC(i,j,k+this%stp1:k+this%stp2)) &
      !         &         -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%itpsc_zm(:,i,j,k)*this%SC(i,j,k+this%stm1:k+this%stm2)) &
      !         &         +this%diff*sum(this%grdsc_z(:,i,j,k)*this%SC(i,j,k-1:k))
            end do
         end do
      end do
      ! Update VF
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
      !         drhoSCdt(i,j,k)=sum(this%divsc_x(:,i,j,k)*FX(i:i+1,j,k))+&
      !         &               sum(this%divsc_y(:,i,j,k)*FY(i,j:j+1,k))+&
      !         &               sum(this%divsc_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      ! Sync result
      call this%cfg%sync(this%VF)
   end subroutine advance
   
   
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
