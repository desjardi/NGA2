!> Connected component labeling class:
!> Provides support for identifying and manipulating Lagrangian objects from a vfs solution
module ccl_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   use irl_fortran_interface
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ccl
   
   !> Object purely for building meta structures for now - holdover from NGA
   type :: struct_type
      ! Location WRT mesh
      integer, dimension(:,:), pointer :: node => null()
      integer :: nnode
      integer :: id
      integer, dimension(3) :: per
      type(struct_type), pointer :: next => null()
   end type struct_type
   
   !> Meta-structure object
   type :: meta_struct_type
      integer :: id
      real(WP) :: vol
      real(WP) :: x
      real(WP) :: y
      real(WP) :: z
      real(WP) :: u
      real(WP) :: v
      real(WP) :: w
      real(WP), dimension(3) :: lengths
      real(WP), dimension(3,3) :: axes
   end type meta_struct_type
   
   !> Structure object
   type :: struct
      integer :: id                                       !< ID of struct
      integer :: parent                                   !< ID of parent struct
      integer :: rank = 0                                 !< Upper bound on height of tree; only used for roots
      integer :: counter = 0                              !< Counter for filling node array
      integer :: nnode                                    !< Number of cells contained in struct
      integer, dimension(:,:), allocatable :: node        !< List of cells contained in struct, dimension(1:nnode,3)
      integer, dimension(3) :: per = 0                    !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct
   
   !> Film object
   type, extends(struct) :: film
      integer :: phase                                    !< Film phase; 1: liquid, 2: gas
      integer :: type                                     !< Film type; 0: droplet, 1: ligament, 2: sheet
      integer, dimension(2) :: adjacent_structs           !< IDs of structures adjacent to gas film
      real(WP) :: min_thickness                           !< Global minimum film thickness
   end type film
   
   !> CCL object definition
   type :: ccl
      
      ! This is our config
      class(config), pointer :: cfg                          !< This is the config object for which the CCL is built
      
      ! This is the name of the CCL
      character(len=str_medium) :: name='UNNAMED_CCL'        !< Solver name (default=UNNAMED_CCL)
      
      ! Volume fraction information
      real(WP), dimension(:,:,:), pointer :: VF              !< Volume fraction array
      
      ! Interface polygon information
      type(Poly_type), dimension(:,:,:,:), pointer :: poly   !< Array of IRL interface polygons (n,i,j,k)
      
      ! Maximum number of PLIC interfaces per cell
      integer :: max_interface_planes                        !< Number of planar interfaces per cell (0=VF-only, 1=PLIC, 2=R2P, etc.)
      
      ! Structure and films
      ! type(struct), dimension(:), allocatable :: struct_list
      ! type(film)  , dimension(:), allocatable :: film_list
      type(struct), dimension(:), pointer :: struct_list => null()
      type(film)  , dimension(:), pointer :: film_list => null()
      
      ! Linked-list struct and meta_struct
      type(struct_type), pointer :: first_struct => null()   !< Pointer to first struct in linked list
      type(struct_type), pointer :: my_struct => null()      !< Pointer to subsequent struct in linked list
      ! list of meta-structures' ID tags
      integer, dimension(:), pointer :: meta_structures => null()
      ! List of meta-structures
      type(meta_struct_type), dimension(:), pointer :: meta_structures_list => null()
      ! List of output variables names
      integer :: meta_structures_nname
      character(len=str_medium), dimension(:), pointer :: meta_structures_name => null()
      ! Time info
      integer :: nout_time
      real(WP), dimension(:), pointer :: out_time => null()
      
      ! CCL selection parameters
      real(WP) :: VFlo=1.0e-10_WP                            !< Minimum VF value considered for a structure to exist
      real(WP) :: dot_threshold=-0.5_WP                      !< Maximum dot product of two interface normals for their respective cells to be considered film cells
      real(WP) :: thickness_cutoff=0.5_WP                    !< Maximum film thickness as fraction of meshsize
      
      ! Feature counts
      integer :: n_struct, n_struct_max !, n_border_struct, n_border_struct_max
      integer :: n_film, n_film_max !, n_border_film, n_border_film_max
      integer :: n_meta_struct
      
      ! Global tag offset
      integer :: id_offset
      integer :: sync_offset, film_sync_offset
      
      ! Local equivalence array (array version of film/struct(:)%parent)
      integer, dimension(:), allocatable :: struct_map_
      integer, dimension(:), allocatable :: film_map_
      
      ! Global equivalence array
      ! parent: unique to each proc; only used as input in allreduce and allgather
      ! parent_all: maximum of all procs (allreduce output)
      ! parent_own: each proc contributes parents for its own structures only (allgather output)
      integer, dimension(:), allocatable :: parent,parent_all,parent_own
      integer, dimension(:), allocatable :: film_parent,film_parent_all,film_parent_own
      
      ! Periodicity array
      integer, dimension(:,:), allocatable :: per_,per
      
      ! Output of the CCL
      integer, dimension(:,:,:), allocatable :: id             !< ID of the structure that contains the cell
      integer, dimension(:,:,:), allocatable :: film_id        !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: film_phase     !< Phase of the film cell - 0/1/2 for none/liquid/gas
      real(WP),dimension(:,:,:), allocatable :: film_thickness !< Local thickness of the film cell
      integer, dimension(:,:,:), allocatable :: film_type      !< Local film type - 0: droplet, 1: ligament, 2: sheet
      real(WP),dimension(:,:,:), allocatable :: film_edge      !< Some measure of whether the film cell is an edge cell 
      real(WP),dimension(:,:,:,:), allocatable :: film_edge_normal !< Outward facing normal of film edge 

      ! Work arrays
      integer, dimension(:,:,:,:), allocatable :: idp          !< ID of the structure that contains the cell
      integer, dimension(:,:,:,:), allocatable :: film_idp     !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: film_pair      !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: border_id      !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: film_border_id !< ID of the film that contains the cell
      real(WP),dimension(:,:,:), allocatable :: SD             !< Surface density array
      
   contains
      procedure :: build_lists
      procedure :: deallocate_lists
      procedure, private :: label
      procedure, private :: calculate_film_thickness
      procedure, private :: is_film_cell_upper
      procedure, private :: struct_sync
      procedure, private :: is_connected
      procedure, private :: film_sync
      procedure, private :: film_is_connected
      procedure, private :: struct_label_update
      procedure, private :: struct_final
      procedure, private :: meta_structures_sort
      procedure, private :: meta_structures_stats
      procedure :: film_classify
      procedure :: get_min_thickness
      procedure :: sort_by_thickness
      procedure, private :: kill_struct
   end type ccl
   
   
   !> Declare CCL algorithm constructor
   interface ccl
      procedure constructor
   end interface ccl
   
contains
   
   
   !> Default constructor for CCL algorithm
   function constructor(cfg,name) result(self)
      use messager, only: die
      implicit none
      type(ccl) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      ! logical :: file_is_there,i
      integer :: ierr
      
      ! Set the name for the object
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to cfg object
      self%cfg=>cfg
      
      ! Allocate id arrays
      allocate(self%id     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%id=0
      allocate(self%film_id(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_id=0
      
      ! Allocate film phase array
      allocate(self%film_phase(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_phase=0
      
      ! Allocate film thickness array
      allocate(self%film_thickness(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_thickness=0.0_WP
      
      ! Allocate film type array
      allocate(self%film_type(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_type=0
      
      ! Allocate film edge array
      allocate(self%film_edge(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_edge=0.0_WP
      allocate(self%film_edge_normal(1:3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_edge_normal=0.0_WP

      ! Allocate periodicity work arrays
      allocate(self%idp     (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%idp=0
      allocate(self%film_idp(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_idp=0
     
      ! Allocate film work arrays
      allocate(self%film_pair (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_pair=0
      
      ! Allocate border id arrays
      allocate(self%border_id     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%border_id=0
      allocate(self%film_border_id(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_border_id=0
      
      ! Allocate surface density array
      allocate(self%SD(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SD=0.0_WP
      
      ! Variable names
      self%meta_structures_nname = 20
      allocate(self%meta_structures_name(self%meta_structures_nname))
      self%meta_structures_name(1) = 'Tag'
      self%meta_structures_name(2) = 'Volume'
      self%meta_structures_name(3) = 'X'
      self%meta_structures_name(4) = 'Y'
      self%meta_structures_name(5) = 'Z'
      self%meta_structures_name(6) = 'U'
      self%meta_structures_name(7) = 'V'
      self%meta_structures_name(8) = 'W'
      self%meta_structures_name(9) = 'L1'
      self%meta_structures_name(10) = 'L2'
      self%meta_structures_name(11) = 'L3'
      self%meta_structures_name(12) = 'V11'
      self%meta_structures_name(13) = 'V12'
      self%meta_structures_name(14) = 'V13'
      self%meta_structures_name(15) = 'V21'
      self%meta_structures_name(16) = 'V22'
      self%meta_structures_name(17) = 'V23'
      self%meta_structures_name(18) = 'V31'
      self%meta_structures_name(19) = 'V32'
      self%meta_structures_name(20) = 'V33'
      ! Open the file
      ! inquire(file='structs/times',exist=file_is_there)
      call MPI_BARRIER(self%cfg%comm,ierr)
      
      ! if (file_is_there) then
      !    ! Read the file
      !    call parser_parsefile('structs/times')
      !    ! Get the time
      !    call parser_getsize('time values',self%nout_time)
      !    allocate(self%out_time(self%nout_time))
      !    call parser_read('time values',self%out_time)
      !    ! Remove future time
      !    future: do i=1,size(self%out_time)
      !       if (self%out_time(i).GE.time*0.99999_WP) then
      !          self%nout_time = i-1
      !          exit future
      !       end if
      !    end do future
      ! else
      ! ! Create directory
      ! if (self%cfg%amRoot) call execute_command_line('mkdir -p structs')
      ! Set the time
      self%nout_time = 0
      allocate(self%out_time(1))
      ! end if
   end function constructor
   
   !> Build lists of structures and films
   subroutine build_lists(this,VF,poly,U,V,W)
      implicit none
      class(ccl), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: VF      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      type(Poly_type), dimension(:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in), optional:: poly !< Needs to be (1:2,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in), optional :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,n
      real(WP) :: tsd
      
      ! Point to volume fraction field
      this%VF(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)=>VF
      
      ! Point to polygon object
      if (present(poly)) then
         this%poly(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)=>poly
         ! Now compute surface area divided by cell volume
         this%SD=0.0_WP
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (VF(i,j,k).eq.0.0_WP) cycle
                  tsd=0.0_WP
                  do n=1,this%max_interface_planes
                     if (getNumberOfVertices(this%poly(n,i,j,k)).gt.0) then
                        tsd=tsd+abs(calculateVolume(this%poly(n,i,j,k)))
                     end if
                  end do
                  this%SD(i,j,k)=tsd/this%cfg%vol(i,j,k)
               end do
            end do
         end do
      else
         this%max_interface_planes = 0
         this%n_film = 0
         this%n_film_max = 0
      end if
      
      ! Label features locally
      call this%label()
      
      ! Synchronize labels across procs
      call this%struct_sync()
      
      if (this%max_interface_planes.eq.2) call this%film_sync()
      
      call this%struct_label_update()
      if (present(U).and.present(V).and.present(W)) call this%struct_final(U,V,W)
      
   end subroutine build_lists
   
   !> Deallocate lists of structures and films
   subroutine deallocate_lists(this)
      implicit none
      class(ccl), intent(inout) :: this
      
      ! Deallocate arrays
      call this%kill_struct()
      
      ! Deallocate list of meta_structures
      if (associated(this%meta_structures).and.associated(this%meta_structures_list)) then
         deallocate(this%meta_structures)
         deallocate(this%meta_structures_list)
         nullify(this%meta_structures)
         nullify(this%meta_structures_list)
      end if
      
   end subroutine deallocate_lists
   
   !> Build local lists of structures and films
   subroutine label(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
      implicit none
      class(ccl), intent(inout) :: this
      integer :: npp,idd,i,j,k,ierr,dim,ii,jj,kk
      integer :: max_structs, max_films
      ! Using PLIC normal information
      logical :: is_contiguous = .true., has_normal, use_normal
      integer, dimension(3) :: pos
      real(WP), dimension(3) :: n1, n2, c1, c2
      logical :: is_film = .false., is_two_plane_film = .false.
      ! Only if two-plane cells are used
      real(WP), dimension(3) :: c22, n22
      
      ! Initialize id (for tagging)
      this%id = 0
      this%film_id = 0
      this%film_phase = 0
      
      ! Reset film thickness
      this%film_thickness=0.0_WP
      
      ! Initialize work arrays
      this%idp = 0
      this%film_idp = 0
      this%film_pair = 0
      this%border_id = 0
      this%film_border_id = 0
      !   n_border_struct = 0 ! counter for local number of structs on processor boundaries
      !   n_border_film   = 0
      
      ! Compute some useful values
      npp = int((this%cfg%nx*this%cfg%ny*this%cfg%nz)/(this%cfg%nproc)) ! number points per processor
      max_structs = ceiling(npp/2.0_WP)
      this%id_offset = npp*(this%cfg%rank) ! global tag offset; unique for each proc
      
      ! Allocate union-find data structure for structs
      allocate(this%struct_list(this%id_offset+1:this%id_offset+max_structs))
      
      ! initialize this%struct_list
      this%struct_list%parent = 0
      this%struct_list%rank = 0
      this%struct_list%counter = 0
      this%struct_list%nnode = 0
      
      ! initialize the global tag to this%id_offset - our "0"
      idd = this%id_offset
      if (this%max_interface_planes.gt.0) then
         ! Allocate union-find data structure for films
         allocate(this%film_list(this%id_offset+1:this%id_offset+max_structs))
         ! Initialize this%film_list
         this%film_list%parent = 0
         this%film_list%rank = 0
         this%film_list%counter = 0
         this%film_list%nnode = 0
         ! outer loop over domain to find "corner" of structure
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  
                  ! Find untagged point in liquid phase
                  if (this%VF(i,j,k).ge.this%VFlo) then
                     
                     is_two_plane_film = .false.
                     has_normal = getNumberOfVertices(this%poly(1,i,j,k)).gt.0
                     if (has_normal) then
                        n1 = calculateNormal(this%poly(1,i,j,k))
                        c1 = calculateCentroid(this%poly(1,i,j,k))
                        ! Rudimentary two-plane cell treatment
                        if (getNumberOfVertices(this%poly(2,i,j,k)).ne.0) then
                           n2 = calculateNormal(this%poly(2,i,j,k))
                           c2 = calculateCentroid(this%poly(2,i,j,k))
                           is_two_plane_film = (dot_product(n1,n2).lt.this%dot_threshold)
                           
                           if (is_two_plane_film) then
                              if (dot_product(c2-c1,n2).gt.0.0_WP) then
                                 this%film_phase(i,j,k) = 1 ! Liquid film
                              else
                                 this%film_phase(i,j,k) = 2 ! Gas film
                              end if
                              if (this%film_phase(i,j,k).eq.1) then ! if liquid film
                                 do dim = 1,3 ! need parallel treatment
                                    pos = 0
                                    pos(dim) = -1
                                    ii = i + pos(1)
                                    jj = j + pos(2)
                                    kk = k + pos(3)
                                    if (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0) then
                                       this%film_phase(ii,jj,kk) = 1
                                    end if
                                    if (this%id(ii,jj,kk).gt.0) then ! Liquid neighbors should already be labeled
                                       if (this%id(i,j,k).ne.0) then
                                          this%id(i,j,k) = union(this%struct_list,this%id(i,j,k),this%id(ii,jj,kk))
                                       else
                                          ! Propagate neighbor label
                                          this%id(i,j,k) = this%id(ii,jj,kk)
                                       end if
                                    end if
                                 end do ! dim = 1,3
                                 if (this%id(i,j,k).eq.0) then
                                    ! Create new label
                                    idd = idd + 1
                                    this%id(i,j,k) = idd
                                    this%struct_list(idd)%parent = idd
                                 end if
                                 ! Periodicity
                                 if (this%cfg%xper .and. i.eq.this%cfg%imax) this%struct_list(this%id(i,j,k))%per(1) = 1
                                 if (this%cfg%yper .and. j.eq.this%cfg%jmax) this%struct_list(this%id(i,j,k))%per(2) = 1
                                 if (this%cfg%zper .and. k.eq.this%cfg%kmax) this%struct_list(this%id(i,j,k))%per(3) = 1
                                 this%idp(:,i,j,k) = this%struct_list(this%id(i,j,k))%per
                              else ! if gas film
                                 do dim = 1,3 ! need parallel treatment
                                    pos = 0
                                    pos(dim) = -1
                                    ii = i + pos(1)
                                    jj = j + pos(2)
                                    kk = k + pos(3)
                                    if (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0) then
                                       this%film_phase(ii,jj,kk) = 2
                                    end if
                                 end do
                              end if
                              do dim = 1,3
                                 pos = 0
                                 pos(dim) = +1
                                 ii = i + pos(1)
                                 jj = j + pos(2)
                                 kk = k + pos(3)
                                 if (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0) then
                                    this%film_phase(ii,jj,kk) = this%film_phase(i,j,k)
                                 end if
                              end do
                           end if
                        end if ! if (getNumberOfVertices(this%poly(2,i,j,k)).ne.0)
                     end if ! if (has_normal)
                     if (.not.is_two_plane_film) then
                        do dim = 1,3
                           pos = 0
                           pos(dim) = -1
                           ii = i + pos(1)
                           jj = j + pos(2)
                           kk = k + pos(3)
                           
                           is_contiguous = .true.
                           is_film = .false.
                           use_normal = has_normal .and. (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0)
                           if (use_normal) then
                              n2 = calculateNormal(this%poly(1,ii,jj,kk))
                              c2 = calculateCentroid(this%poly(1,ii,jj,kk))
                              ! If neighbor is two-plane cell
                              if (getNumberOfVertices(this%poly(2,ii,jj,kk)).ne.0) then
                                 c22 = calculateCentroid(this%poly(2,ii,jj,kk))
                                 n22 = calculateNormal(this%poly(2,ii,jj,kk))
                                 is_contiguous = dot_product((c22-c2),n22).gt.0.0_WP ! .true. if liquid film
                                 is_film = .true.
                              else
                                 ! If neighbor is one-plane cell
                                 is_contiguous = (dot_product(c2-c1,n2).ge.0.0_WP).or.(dot_product(c1-c2,n1).ge.0.0_WP)
                                 is_film = (dot_product(n1,n2).lt.this%dot_threshold)
                              end if
                              
                              if (is_film) then
                                 if (is_contiguous) then
                                    ! Liquid film
                                    this%film_phase(i,j,k) = 1
                                    this%film_phase(ii,jj,kk) = 1
                                 else
                                    ! Gas film
                                    this%film_phase(i,j,k) = 2
                                    this%film_phase(ii,jj,kk) = 2
                                    this%film_pair(i,j,k) = this%cfg%get_lexico_from_ijk([ii,jj,kk])
                                    this%film_pair(ii,jj,kk) = this%cfg%get_lexico_from_ijk([i,j,k])
                                 end if ! is_contiguous
                              else
                                 ! is_contiguous = is_contiguous.or.this%VF(i,j,k).gt.0.1_WP.or.this%VF(ii,jj,kk).gt.0.1_WP
                                 ! is_contiguous = is_contiguous.or.this%VF(ii,jj,kk).gt.0.1_WP
                                 ! is_contiguous = is_contiguous.or.norm2(c2-c1).le.0.5_WP*(this%cfg%meshsize(i,j,k)+this%cfg%meshsize(ii,jj,kk))
                                 is_contiguous = is_contiguous.or.norm2(c2-c1).le.1.0_WP*norm2([this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]-[this%cfg%xm(ii),this%cfg%ym(jj),this%cfg%zm(kk)])
                              end if ! is_film
                              
                           end if ! use_normal
                           if (this%id(ii,jj,kk).gt.0 .and. is_contiguous) then ! Liquid neighbors should already be labeled
                              if (this%id(i,j,k).ne.0) then
                                 this%id(i,j,k) = union(this%struct_list,this%id(i,j,k),this%id(ii,jj,kk))
                              else
                                 ! Propagate neighbor label
                                 this%id(i,j,k) = this%id(ii,jj,kk)
                              end if
                           end if
                        end do ! dim = 1,3
                        if (this%id(i,j,k).eq.0) then
                           ! Create new label
                           idd = idd + 1
                           this%id(i,j,k) = idd
                           this%struct_list(idd)%parent = idd
                        end if
                        ! Periodicity
                        if (this%cfg%xper .and. i.eq.this%cfg%imax) this%struct_list(this%id(i,j,k))%per(1) = 1
                        if (this%cfg%yper .and. j.eq.this%cfg%jmax) this%struct_list(this%id(i,j,k))%per(2) = 1
                        if (this%cfg%zper .and. k.eq.this%cfg%kmax) this%struct_list(this%id(i,j,k))%per(3) = 1
                        this%idp(:,i,j,k) = this%struct_list(this%id(i,j,k))%per
                     end if ! .not.is_two_plane_film
                  end if ! this%VF >= this%VFlo
               end do ! k
            end do ! j
         end do ! i
         
         ! Check upper boundaries for film cells
         ! this%cfg%imax_
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (this%VF(this%cfg%imax_,j,k).ge.this%VFlo .and. this%film_phase(this%cfg%imax_,j,k).eq.0) then
                  this%film_phase(this%cfg%imax_,j,k) = this%is_film_cell_upper(this%cfg%imax_,j,k,1)
               end if
            end do
         end do
         ! this%cfg%jmax_
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (this%VF(i,this%cfg%jmax_,k).ge.this%VFlo .and. this%film_phase(i,this%cfg%jmax_,k).eq.0) then
                  this%film_phase(i,this%cfg%jmax_,k) = this%is_film_cell_upper(i,this%cfg%jmax_,k,2)
               end if
            end do
         end do
         ! this%cfg%kmax_
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               if (this%VF(i,j,this%cfg%kmax_).ge.this%VFlo .and. this%film_phase(i,j,this%cfg%kmax_).eq.0) then
                  this%film_phase(i,j,this%cfg%kmax_) = this%is_film_cell_upper(i,j,this%cfg%kmax_,3)
               end if
            end do
         end do
         
         ! Number of this%struct_list graph nodes
         max_structs = idd-this%id_offset
         
         ! Reset idd
         idd = this%id_offset
         ! Tag films
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  ! Find untagged point in film
                  if (this%film_phase(i,j,k).gt.0) then
                     if (this%calculate_film_thickness(i,j,k).lt.this%thickness_cutoff*this%cfg%meshsize(i,j,k)) then
                        do dim = 1,3
                           pos = 0
                           pos(dim) = -1
                           ii = i + pos(1)
                           jj = j + pos(2)
                           kk = k + pos(3)
                           if (this%film_id(ii,jj,kk).gt.0) then ! Liquid neighbors should already be labeled
                              if (this%film_id(i,j,k).ne.0) then
                                 this%film_id(i,j,k) = union(this%film_list,this%film_id(i,j,k),this%film_id(ii,jj,kk))
                              else
                                 ! Propagate neighbor label
                                 this%film_id(i,j,k) = this%film_id(ii,jj,kk)
                              end if
                           end if
                           if (this%film_id(i,j,k).eq.0) then
                              ! Create new label
                              idd = idd + 1
                              this%film_id(i,j,k) = idd
                              this%film_list(idd)%parent = idd
                           end if
                        end do
                        ! Periodicity
                        if (this%cfg%xper .and. i.eq.this%cfg%imax) this%film_list(this%film_id(i,j,k))%per(1) = 1
                        if (this%cfg%yper .and. j.eq.this%cfg%jmax) this%film_list(this%film_id(i,j,k))%per(2) = 1
                        if (this%cfg%zper .and. k.eq.this%cfg%kmax) this%film_list(this%film_id(i,j,k))%per(3) = 1
                        this%film_idp(:,i,j,k) = this%film_list(this%film_id(i,j,k))%per
                     else
                        this%film_phase(i,j,k) = 0
                     end if
                  end if
               end do ! k
            end do ! j
         end do ! i
         ! Update ghost cells
         call this%cfg%sync(this%film_phase)
         call this%cfg%sync(this%film_thickness)
         
         ! Number of this%film_list graph nodes
         max_films = idd-this%id_offset
         
         ! Collapse tree
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  
                  if (this%film_id(i,j,k).gt.0) then
                     this%film_id(i,j,k) = find(this%film_list,this%film_id(i,j,k))
                     this%film_list(this%film_id(i,j,k))%nnode = this%film_list(this%film_id(i,j,k))%nnode + 1
                     this%film_idp(1,i,j,k) = max(this%film_idp(1,i,j,k),this%film_list(this%film_id(i,j,k))%per(1))
                     this%film_idp(2,i,j,k) = max(this%film_idp(2,i,j,k),this%film_list(this%film_id(i,j,k))%per(2))
                     this%film_idp(3,i,j,k) = max(this%film_idp(3,i,j,k),this%film_list(this%film_id(i,j,k))%per(3))
                     this%film_list(this%film_id(i,j,k))%per = this%film_idp(:,i,j,k)
                  end if
                  
               end do ! k
            end do ! j
         end do ! i

         ! Allocate this%film_list%node array, calculate this%n_film
         this%n_film = 0
         do i=this%id_offset+1,this%id_offset+max_films
            if (this%film_list(i)%nnode.gt.0) then
               allocate(this%film_list(i)%node(3,this%film_list(i)%nnode))
               this%n_film = this%n_film + 1
               ! if (sum(this%film_list(i)%per).gt.0) n_border_struct = n_border_struct + 1
            end if
         end do

         ! total films or parts of films  across all procs
         call MPI_ALLREDUCE(this%n_film,this%n_film_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         ! call MPI_ALLREDUCE(this%n_film_border_struct,this%n_film_border_struct_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)

         this%film_sync_offset = this%cfg%rank*this%n_film_max

         allocate(this%film_map_(this%film_sync_offset+1:this%film_sync_offset+this%n_film_max)); this%film_map_ = 0
         idd = this%film_sync_offset
         do i=this%id_offset+1,this%id_offset+max_films
            ! Only if this%film_list has nodes
            if (this%film_list(i)%nnode.eq.0) cycle
            idd = idd + 1
            this%film_map_(idd) = i
            this%film_list(i)%id = idd
         end do
         call fill_border_id(this%film_list,this%film_id,this%film_border_id)

         ! Fill this%film_list%node array
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_

                  if (this%film_phase(i,j,k).gt.0) then
                     this%film_list(this%film_id(i,j,k))%counter = this%film_list(this%film_id(i,j,k))%counter + 1
                     this%film_list(this%film_id(i,j,k))%node(1,this%film_list(this%film_id(i,j,k))%counter) = i
                     this%film_list(this%film_id(i,j,k))%node(2,this%film_list(this%film_id(i,j,k))%counter) = j
                     this%film_list(this%film_id(i,j,k))%node(3,this%film_list(this%film_id(i,j,k))%counter) = k
                  end if
                  
               end do ! i
            end do ! j
         end do ! k
      
         ! Film phase
         film_phase: block
            integer :: m
            do m=this%film_sync_offset+1,this%film_sync_offset+this%n_film
               i = this%film_list(this%film_map_(m))%node(1,1)
               j = this%film_list(this%film_map_(m))%node(2,1)
               k = this%film_list(this%film_map_(m))%node(3,1)
               this%film_list(this%film_map_(m))%phase = this%film_phase(i,j,k)
            end do
         end block film_phase

      else ! if no interface polygons given
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  
                  ! Find untagged point in liquid phase
                  if (this%VF(i,j,k).ge.this%VFlo) then
                     
                     do dim = 1,3
                        pos = 0
                        pos(dim) = -1
                        ii = i + pos(1)
                        jj = j + pos(2)
                        kk = k + pos(3)
                        
                        if (this%id(ii,jj,kk).gt.0) then ! Liquid neighbors should already be labeled
                           if (this%id(i,j,k).ne.0) then
                              this%id(i,j,k) = union(this%struct_list,this%id(i,j,k),this%id(ii,jj,kk))
                           else
                              ! Propagate neighbor label
                              this%id(i,j,k) = this%id(ii,jj,kk)
                           end if
                        end if
                     end do ! dim = 1,3
                     if (this%id(i,j,k).eq.0) then
                        ! Create new label
                        idd = idd + 1
                        this%id(i,j,k) = idd
                        this%struct_list(idd)%parent = idd
                     end if
                     ! Periodicity
                     if (this%cfg%xper .and. i.eq.this%cfg%imax) this%struct_list(this%id(i,j,k))%per(1) = 1
                     if (this%cfg%yper .and. j.eq.this%cfg%jmax) this%struct_list(this%id(i,j,k))%per(2) = 1
                     if (this%cfg%zper .and. k.eq.this%cfg%kmax) this%struct_list(this%id(i,j,k))%per(3) = 1
                     this%idp(:,i,j,k) = this%struct_list(this%id(i,j,k))%per
                  end if ! this%VF >= this%VFlo
               end do ! k
            end do ! j
         end do ! i
         
         ! Number of this%struct_list graph nodes
         max_structs = idd-this%id_offset
         
      end if
      
      ! Collapse tree
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               
               if (this%id(i,j,k).gt.0) then
                  this%id(i,j,k) = find(this%struct_list,this%id(i,j,k))
                  this%struct_list(this%id(i,j,k))%nnode = this%struct_list(this%id(i,j,k))%nnode + 1
                  this%idp(1,i,j,k) = max(this%idp(1,i,j,k),this%struct_list(this%id(i,j,k))%per(1))
                  this%idp(2,i,j,k) = max(this%idp(2,i,j,k),this%struct_list(this%id(i,j,k))%per(2))
                  this%idp(3,i,j,k) = max(this%idp(3,i,j,k),this%struct_list(this%id(i,j,k))%per(3))
                  this%struct_list(this%id(i,j,k))%per = this%idp(:,i,j,k)
               end if
               
            end do ! k
         end do ! j
      end do ! i
      ! At this point, this%struct_list(this%id(i,j,k))%nnode > 0 only if this%id(i,j,k) is the root of a structure tree
      ! Allocate this%struct_list%node array, calculate this%n_struct
      this%n_struct = 0 ! counter for local number of structs
      do i=this%id_offset+1,this%id_offset+max_structs
         if (this%struct_list(i)%nnode.gt.0) then
            allocate(this%struct_list(i)%node(3,this%struct_list(i)%nnode))
            this%n_struct = this%n_struct + 1
            ! if (sum(this%struct_list(i)%per).gt.0) n_border_struct = n_border_struct + 1
         end if
      end do

      ! total structures or parts of structures across all procs
      call MPI_ALLREDUCE(this%n_struct,this%n_struct_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      ! call MPI_ALLREDUCE(n_border_struct,n_border_struct_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      
      !! Whenever we switch to variable struct number per proc
      ! call MPI_ALLGATHER(this%n_struct, 1, MPI_INTEGER, this%n_struct_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      ! call MPI_ALLGATHER(n_border_struct, 1, MPI_INTEGER, n_border_struct_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      ! call MPI_ALLGATHER(this%n_film, 1, MPI_INTEGER, this%n_film_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      ! call MPI_ALLGATHER(this%n_film_border_struct, 1, MPI_INTEGER, this%n_film_border_struct_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      !   this%sync_offset = sum(this%n_struct_per_processor(0:this%cfg%rank)) ???? not sum(this%n_struct_per_processor(0:this%cfg%rank)) - this%n_struct_per_processor(this%cfg%rank)??
      !   this%sync_offset = sum(n_border_struct_per_processor(0:this%cfg%rank))
      !   this%film_sync_offset = sum(this%n_film_per_processor(0:this%cfg%rank))
      !   this%film_sync_offset = sum(this%n_film_border_struct_per_processor(0:this%cfg%rank))
      !   allocate(this%struct_map_(this%sync_offset+1:this%sync_offset+this%n_struct)); this%struct_map_ = 0
      !   allocate(this%film_map_(this%film_sync_offset+1:this%film_sync_offset+this%n_film)); this%film_map_ = 0
      
      this%sync_offset = this%cfg%rank*this%n_struct_max
      ! Could allocate this%struct_map_(max_structs) and integrate this loop with above
      ! allocate(parent_    (this%sync_offset+1:this%sync_offset+this%n_struct_max)); parent_ = 0
      allocate(this%struct_map_(this%sync_offset+1:this%sync_offset+this%n_struct_max)); this%struct_map_ = 0
      idd = this%sync_offset
      do i=this%id_offset+1,this%id_offset+max_structs
         ! Only if this%struct_list has nodes
         if (this%struct_list(i)%nnode.eq.0) cycle
         idd = idd + 1
         this%struct_map_(idd) = i
         this%struct_list(i)%id = idd
      end do

      
      call fill_border_id(this%struct_list,this%id,this%border_id)
      
      ! Fill this%struct_list%node array
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%id(i,j,k).gt.0) then
                  this%struct_list(this%id(i,j,k))%counter = this%struct_list(this%id(i,j,k))%counter + 1
                  this%struct_list(this%id(i,j,k))%node(1,this%struct_list(this%id(i,j,k))%counter) = i
                  this%struct_list(this%id(i,j,k))%node(2,this%struct_list(this%id(i,j,k))%counter) = j
                  this%struct_list(this%id(i,j,k))%node(3,this%struct_list(this%id(i,j,k))%counter) = k
               end if
            end do ! i
         end do ! j
      end do ! k

      ! Copy this%struct_list to my_struct
      call copy_to_my_struct
      
      return
      
   contains
      
      ! function find(a_struct,a_x) result(a_y)
      !    implicit none
      !    type(struct), dimension(this%id_offset+1:), intent(inout) :: a_struct
      !    integer :: a_x
      !    integer :: a_y
      !    ! integer :: a_z
      
      !    a_y = a_x
      !    ! if(irank.eq.2) print *,"irank",irank,"a_y",a_y,"parent",a_struct(a_y)%parent
      !    ! find root of x with path halving
      !    do while ( a_y.ne.a_struct(a_y)%parent )
      !       a_struct(a_y)%parent = a_struct(a_struct(a_y)%parent)%parent
      !       a_y = a_struct(a_y)%parent
      !    end do
      
      !    ! ! alternate way of collapsing tree
      !    ! do while (a_struct(a_x)%parent.ne.a_x)
      !    !    a_z = a_struct(a_x)%parent
      !    !    a_struct(a_x)%parent = a_y
      !    !    a_x = a_z
      !    ! end do
      
      !    return
      ! end function find
      
      ! This function points the parent to root and returns that root
      recursive function find(a_struct,a_x) result(a_y)
         implicit none
         class(struct), dimension(this%id_offset+1:), intent(inout) :: a_struct
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.a_struct(a_y)%parent ) then
            a_struct(a_y)%parent = find(a_struct,a_struct(a_y)%parent)
            a_y = a_struct(a_y)%parent
         end if
         return
      end function find
      
      ! This function joins two branches at their roots
      ! At the moment, it joints them at the root of a_y and returns that root
      function union(a_struct, a_x, a_y) result(a_z)
         implicit none
         class(struct), dimension(this%id_offset+1:), intent(inout) :: a_struct
         integer, intent(in) :: a_x
         integer, intent(in) :: a_y
         integer :: a_z
         
         ! if a_x tree smaller/lower rank than a_y tree
         a_z = find(a_struct,a_y)
         a_struct(find(a_struct,a_x))%parent = a_z
         
         return
      end function union
      
      ! ! Union-by-rank
      ! function union(a_struct, a_x, a_y) result(a_z)
      !    implicit none
      !    class(struct), dimension(this%id_offset+1:), intent(inout) :: a_struct
      !    ! integer, intent(inout) :: a_x
      !    ! integer, intent(inout) :: a_y
      !    integer :: a_z
      
      !    ! if a_x tree smaller/lower rank than a_y tree
      !    a_x = find(a_struct,a_x)
      !    a_y = find(a_struct,a_y)
      !    if (a_x.eq.a_y) a_z = a_x; return
      !    ! if (a_struct(a_x)%rank.ge.a_struct(a_y)%rank) then
      !    !    a_z = a_x
      !    !    a_x = a_y
      !    !    a_y = a_z
      !    ! end if
      !    ! a_struct(a_x)%parent = a_y
      !    ! if (a_struct(a_x)%rank.eq.a_struct(a_y)%rank) then
      !    !    a_struct(a_y)%rank = a_struct(a_y)%rank + 1
      !    ! end if
      !    if (a_struct(a_x)%rank.gt.a_struct(a_y)%rank) then
      !       a_struct(a_y)%parent = a_x
      !    else
      !       a_struct(a_x)%parent = a_y
      !       if (a_struct(a_x)%rank.eq.a_struct(a_y)%rank) then
      !          a_struct(a_y)%rank = a_struct(a_y)%rank + 1
      !       end if
      !    end if
      !    return
      ! end function union
      
      ! Fill border cell id array with compact IDs
      subroutine fill_border_id(a_list_type,a_id_array,a_border_array)
         implicit none
         class(struct), dimension(this%id_offset+1:), intent(in) :: a_list_type
         integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(in)   :: a_id_array
         integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(out)  :: a_border_array
         ! this%cfg%imin_
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (a_id_array(this%cfg%imin_,j,k).gt.0) then
                  a_border_array(this%cfg%imin_,j,k) = a_list_type(a_id_array(this%cfg%imin_,j,k))%id
               end if
            end do
         end do
         
         ! this%cfg%imax_
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (a_id_array(this%cfg%imax_,j,k).gt.0) then
                  a_border_array(this%cfg%imax_,j,k) = a_list_type(a_id_array(this%cfg%imax_,j,k))%id
               end if
            end do
         end do
         
         ! this%cfg%jmin_
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (a_id_array(i,this%cfg%jmin_,k).gt.0) then
                  a_border_array(i,this%cfg%jmin_,k) = a_list_type(a_id_array(i,this%cfg%jmin_,k))%id
               end if
            end do
         end do
         
         ! this%cfg%jmax_
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (a_id_array(i,this%cfg%jmax_,k).gt.0) then
                  a_border_array(i,this%cfg%jmax_,k) = a_list_type(a_id_array(i,this%cfg%jmax_,k))%id
               end if
            end do
         end do
         
         ! this%cfg%kmin_
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               if (a_id_array(i,j,this%cfg%kmin_).gt.0) then
                  a_border_array(i,j,this%cfg%kmin_) = a_list_type(a_id_array(i,j,this%cfg%kmin_))%id
               end if
            end do
         end do
         
         ! this%cfg%kmax_
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               if (a_id_array(i,j,this%cfg%kmax_).gt.0) then
                  a_border_array(i,j,this%cfg%kmax_) = a_list_type(a_id_array(i,j,this%cfg%kmax_))%id
               end if
            end do
         end do
      end subroutine fill_border_id
      
      subroutine copy_to_my_struct
         implicit none
         
         do idd = this%id_offset+1,this%id_offset+max_structs
            ! Only allocate if this%struct_list has nodes
            if (this%struct_list(idd)%nnode.eq.0) cycle
            ! Allocate next structure
            if(.not.associated(this%first_struct)) then
               allocate(this%first_struct)
               this%my_struct => this%first_struct
            else
               allocate(this%my_struct%next)
               this%my_struct => this%my_struct%next
            end if
            this%my_struct%nnode = this%struct_list(idd)%nnode
            this%my_struct%id = idd
            this%my_struct%per = this%struct_list(idd)%per
            this%my_struct%node => this%struct_list(idd)%node
            nullify(this%my_struct%next)
         end do
      end subroutine copy_to_my_struct
      
   end subroutine label
   
   ! Calculate film thickness
   function calculate_film_thickness(this,i,j,k) result(local_thickness)
      implicit none
      class(ccl), intent(inout) :: this
      integer, intent(in) :: i,j,k
      real(WP) :: local_thickness
      real(WP) :: SD_local_sum, VOF_local_sum
      integer :: ii,jj,kk
      
      SD_local_sum = 0.0_WP
      VOF_local_sum = 0.0_WP
      if (this%film_phase(i,j,k).eq.2) then
         do kk = k-1,k+1
            do jj = j-1,j+1
               do ii = i-1,i+1
                  VOF_local_sum = VOF_local_sum + (1.0_WP-this%VF(ii,jj,kk))
               end do
            end do
         end do
      else
         do kk = k-1,k+1
            do jj = j-1,j+1
               do ii = i-1,i+1
                  VOF_local_sum = VOF_local_sum + this%VF(ii,jj,kk)
               end do
            end do
         end do
      end if
      do kk = k-1,k+1
         do jj = j-1,j+1
            do ii = i-1,i+1
               SD_local_sum = SD_local_sum + this%SD(ii,jj,kk)
            end do
         end do
      end do
      local_thickness = 2.0_WP*VOF_local_sum/(SD_local_sum+epsilon(1.0_WP))
      this%film_thickness(i,j,k) = local_thickness
      return
   end function calculate_film_thickness

   ! Is the cell above (i,j,k) in the dim direction a film cell?
   function is_film_cell_upper(this,i,j,k,dim) result(film_phase)
      implicit none
      class(ccl), intent(inout) :: this
      integer, intent(in) :: i,j,k,dim
      integer :: ii,jj,kk
      logical :: is_film, use_normal, is_contiguous
      integer, dimension(3) :: pos
      real(WP), dimension(3) :: nref, nloc, pref, ploc
      integer :: film_phase
      
      pos = 0
      pos(dim) = +1
      ii = i + pos(1)
      jj = j + pos(2)
      kk = k + pos(3)
      film_phase = 0
      use_normal = getNumberOfVertices(this%poly(1,i,j,k)).gt.0 .and. (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0)
      if (use_normal) then
         ! If neighbor is a two-plane cell
         if (getNumberOfVertices(this%poly(2,ii,jj,kk)).ne.0) then
            nref = calculateNormal(this%poly(1,ii,jj,kk))
            nloc = calculateNormal(this%poly(2,ii,jj,kk))
            pref = calculateCentroid(this%poly(1,ii,jj,kk))
            ploc = calculateCentroid(this%poly(2,ii,jj,kk))
            is_contiguous = dot_product(ploc-pref,nloc).gt.0.0_WP ! .true. if liquid film
            is_film = .true.
         ! If self is a two-plane cell
         elseif (getNumberOfVertices(this%poly(2,i,j,k)).ne.0) then
            nref = calculateNormal(this%poly(1,i,j,k))
            nloc = calculateNormal(this%poly(2,i,j,k))
            pref = calculateCentroid(this%poly(1,i,j,k))
            ploc = calculateCentroid(this%poly(2,i,j,k))
            is_contiguous = dot_product(ploc-pref,nloc).gt.0.0_WP ! .true. if liquid film
            is_film = .true.
         ! If both are one-plane cells
         else
            nref = calculateNormal(this%poly(1,i,j,k))
            nloc = calculateNormal(this%poly(1,ii,jj,kk))
            pref = calculateCentroid(this%poly(1,i,j,k))
            ploc = calculateCentroid(this%poly(1,ii,jj,kk))
            is_contiguous = (dot_product(ploc-pref,nloc).ge.0.0_WP).or.(dot_product(pref-ploc,nref).ge.0.0_WP)
            is_film = (dot_product(nref,nloc).lt.this%dot_threshold)
         end if
         if (is_film) then
            if (is_contiguous) then
               ! Liquid film
               film_phase = 1
            else
               ! Gas film
               film_phase = 2
               this%film_pair(i,j,k) = this%cfg%get_lexico_from_ijk([ii,jj,kk])
            end if ! is_contiguous
         end if ! is_film
      end if
      return
   end function is_film_cell_upper
   
   !> Synchronize struct labels across all procs
   subroutine struct_sync(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
      implicit none
      class(ccl), intent(inout) :: this
      integer :: i,j,k,n,ii,jj,kk,stop_,stop_global,counter,ierr
      integer :: find_parent,find_parent_own
      
      ! Eventually replace with sum(n_struct_per_processor(1:irank)) or border_struct equivalent
      ! Also all "fill parent" loops and allreduces
      allocate(this%parent(this%cfg%nproc*this%n_struct_max)); this%parent = 0
      allocate(this%parent_all(this%cfg%nproc*this%n_struct_max)); this%parent_all = 0
      allocate(this%parent_own(this%cfg%nproc*this%n_struct_max)); this%parent_own = 0
      
      ! Fill parent with selves
      do i= 1,this%cfg%nproc*this%n_struct_max
         this%parent(i) = i
      end do
      
      ! update ghost cells
      call this%cfg%sync(this%border_id)
      
      ! imin_
      if (this%cfg%imin_.ne.this%cfg%imin) then
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%border_id(this%cfg%imin_,j,k).gt.0.and.this%border_id(this%cfg%imin_-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%border_id(this%cfg%imin_,j,k)).eq.this%border_id(this%cfg%imin_,j,k)) then
                     ! call union(this%border_id(this%cfg%imin_,j,k),this%border_id(this%cfg%imin_-1,j,k))
                     ! is_contiguous = is_connected(this%cfg%imin_,j,k,1)
                     if (this%is_connected(this%cfg%imin_,j,k,1)) this%parent(this%border_id(this%cfg%imin_,j,k)) = this%border_id(this%cfg%imin_-1,j,k)
                  elseif (find(this%border_id(this%cfg%imin_,j,k)).ne.find(this%border_id(this%cfg%imin_-1,j,k))) then
                     ! call union(this%border_id(this%cfg%imin_-1,j,k),this%parent(this%border_id(this%cfg%imin_,j,k)))
                     if (this%is_connected(this%cfg%imin_,j,k,1)) then
                        if (this%border_id(this%cfg%imin_-1,j,k).gt.this%parent(this%border_id(this%cfg%imin_,j,k))) then
                           call union(this%border_id(this%cfg%imin_-1,j,k),this%parent(this%border_id(this%cfg%imin_,j,k)))
                        else
                           call union(this%parent(this%border_id(this%cfg%imin_,j,k)),this%border_id(this%cfg%imin_-1,j,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%jmin_
      if (this%cfg%jmin_.ne.this%cfg%jmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%border_id(i,this%cfg%jmin_,k).gt.0.and.this%border_id(i,this%cfg%jmin_-1,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%border_id(i,this%cfg%jmin_,k)).eq.this%border_id(i,this%cfg%jmin_,k)) then
                     ! call union(this%border_id(i,this%cfg%jmin_,k),this%border_id(i,this%cfg%jmin_-1,k))
                     if (this%is_connected(i,this%cfg%jmin_,k,2)) this%parent(this%border_id(i,this%cfg%jmin_,k)) = this%border_id(i,this%cfg%jmin_-1,k)
                  elseif (find(this%border_id(i,this%cfg%jmin_,k)).ne.find(this%border_id(i,this%cfg%jmin_-1,k))) then
                     ! call union(this%border_id(i,this%cfg%jmin_-1,k),this%parent(this%border_id(i,this%cfg%jmin_,k)))
                     if (this%is_connected(i,this%cfg%jmin_,k,2)) then
                        if (this%border_id(i,this%cfg%jmin_-1,k).gt.this%parent(this%border_id(i,this%cfg%jmin_,k))) then
                           call union(this%border_id(i,this%cfg%jmin_-1,k),this%parent(this%border_id(i,this%cfg%jmin_,k)))
                        else
                           call union(this%parent(this%border_id(i,this%cfg%jmin_,k)),this%border_id(i,this%cfg%jmin_-1,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%kmin_
      if (this%cfg%kmin_.ne.this%cfg%kmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               ! If the neighboring ghost cell is filled
               if (this%border_id(i,j,this%cfg%kmin_).gt.0.and.this%border_id(i,j,this%cfg%kmin_-1).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%border_id(i,j,this%cfg%kmin_)).eq.this%border_id(i,j,this%cfg%kmin_)) then
                     ! call union(this%border_id(i,j,this%cfg%kmin_),this%border_id(i,j,this%cfg%kmin_-1))
                     if (this%is_connected(i,j,this%cfg%kmin_,3)) this%parent(this%border_id(i,j,this%cfg%kmin_)) = this%border_id(i,j,this%cfg%kmin_-1)
                  elseif (find(this%border_id(i,j,this%cfg%kmin_)).ne.find(this%border_id(i,j,this%cfg%kmin_-1))) then
                     ! call union(this%border_id(i,j,this%cfg%kmin_-1),this%parent(this%border_id(i,j,this%cfg%kmin_)))
                     if (this%is_connected(i,j,this%cfg%kmin_,3)) then
                        if (this%border_id(i,j,this%cfg%kmin_-1).gt.this%parent(this%border_id(i,j,this%cfg%kmin_))) then
                           call union(this%border_id(i,j,this%cfg%kmin_-1),this%parent(this%border_id(i,j,this%cfg%kmin_)))
                        else
                           call union(this%parent(this%border_id(i,j,this%cfg%kmin_)),this%border_id(i,j,this%cfg%kmin_-1))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! initialize global stop criterion
      stop_global = 1
      
      ! initialize a counter
      counter = 0
      
      do while (stop_global.ne.0)
         
         ! Initialize local stop flag
         stop_ = 0
         
         this%parent_own = this%parent
         ! Set self-parents to huge(1)
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent(i).eq.i) this%parent(i) = huge(1)
         end do
         
         call MPI_ALLREDUCE(this%parent,this%parent_all,this%cfg%nproc*this%n_struct_max,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
         
         ! Set self-parents back to selves
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent_all(i).eq.huge(1)) this%parent_all(i) = i
            ! if (this%parent_own(i).eq.huge(1)) this%parent_own(i) = i
         end do
         
         ! Flatten trees - is this necessary?
         do i= 1,this%cfg%nproc*this%n_struct_max
            this%parent_all(i) = find_all(i)
            this%parent_own(i) = find_own(i)
         end do
         
         ! Start with final this%parent array being equal to this%parent_all
         this%parent = this%parent_all
         
         ! Count how many global updates we do
         counter = counter + 1
         
         ! Reconcile conflicts between this%parent_all and this%parent_own
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent_own(i).ne.i) then
               find_parent_own = find(this%parent_own(i))
               find_parent = find(this%parent(i))
               if (find_parent_own.gt.find_parent) then
                  call union(this%parent_own(i),this%parent(i))
                  stop_ = 1
               else if (find_parent.gt.find_parent_own) then
                  call union(this%parent(i),this%parent_own(i))
                  stop_ = 1
               end if
            end if
         end do
         ! Check if we did some changes
         call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         
      end do ! do while (stop_global.ne.0) excluding domain boundaries
      
      ! Update this%struct_list%parent and point all parents to root
      do i=this%sync_offset+1,this%sync_offset+this%n_struct
         this%struct_list(this%struct_map_(i))%parent = find(this%parent(i))
      end do
      
      ! Update id array with compacted and syncrhonized ids
      do i=this%sync_offset+1,this%sync_offset+this%n_struct
         do n=1,this%struct_list(this%struct_map_(i))%nnode
            ii = this%struct_list(this%struct_map_(i))%node(1,n)
            jj = this%struct_list(this%struct_map_(i))%node(2,n)
            kk = this%struct_list(this%struct_map_(i))%node(3,n)
            this%id(ii,jj,kk) = this%struct_list(this%struct_map_(i))%parent
         end do
      end do
      
      ! Update periodicity array across processors
      call struct_sync_per
      
      !! Update domain boundaries
      ! this%cfg%imin
      if (this%cfg%imin_.eq.this%cfg%imin) then
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%id(this%cfg%imin,j,k).gt.0.and.this%border_id(this%cfg%imin-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%id(this%cfg%imin,j,k)).eq.this%id(this%cfg%imin,j,k)) then
                     if (this%is_connected(this%cfg%imin,j,k,1)) call union(this%id(this%cfg%imin,j,k),this%border_id(this%cfg%imin-1,j,k))
                     ! this%parent(this%id(this%cfg%imin,j,k)) = this%border_id(this%cfg%imin-1,j,k)
                  elseif (find(this%id(this%cfg%imin,j,k)).ne.find(this%border_id(this%cfg%imin-1,j,k))) then
                     if (this%is_connected(this%cfg%imin,j,k,1)) then
                        if (this%border_id(this%cfg%imin-1,j,k).gt.this%parent(this%id(this%cfg%imin,j,k))) then
                           call union(this%border_id(this%cfg%imin-1,j,k),this%parent(this%id(this%cfg%imin,j,k)))
                        else
                           call union(this%parent(this%id(this%cfg%imin,j,k)),this%border_id(this%cfg%imin-1,j,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%jmin
      if (this%cfg%jmin_.eq.this%cfg%jmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the neighboring ghost cell is filled
               if (this%id(i,this%cfg%jmin,k).gt.0.and.this%border_id(i,this%cfg%jmin-1,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%id(i,this%cfg%jmin,k)).eq.this%id(i,this%cfg%jmin,k)) then
                     if (this%is_connected(i,this%cfg%jmin,k,2)) call union(this%id(i,this%cfg%jmin,k),this%border_id(i,this%cfg%jmin-1,k))
                     ! this%parent(this%id(i,this%cfg%jmin,k)) = this%border_id(i,this%cfg%jmin-1,k)
                  elseif (find(this%id(i,this%cfg%jmin,k)).ne.find(this%border_id(i,this%cfg%jmin-1,k))) then
                     ! call union(this%border_id(i,this%cfg%jmin-1,k),this%parent(this%id(i,this%cfg%jmin,k)))
                     if (this%is_connected(i,this%cfg%jmin,k,2)) then
                        if (this%border_id(i,this%cfg%jmin-1,k).gt.this%parent(this%id(i,this%cfg%jmin,k))) then
                           call union(this%border_id(i,this%cfg%jmin-1,k),this%parent(this%id(i,this%cfg%jmin,k)))
                        else
                           call union(this%parent(this%id(i,this%cfg%jmin,k)),this%border_id(i,this%cfg%jmin-1,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%kmin
      if (this%cfg%kmin_.eq.this%cfg%kmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               ! If the neighboring ghost cell is filled
               if (this%id(i,j,this%cfg%kmin).gt.0.and.this%border_id(i,j,this%cfg%kmin-1).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%id(i,j,this%cfg%kmin)).eq.this%id(i,j,this%cfg%kmin)) then
                     if (this%is_connected(i,j,this%cfg%kmin,3)) call union(this%id(i,j,this%cfg%kmin),this%border_id(i,j,this%cfg%kmin-1))
                     ! this%parent(this%id(i,j,this%cfg%kmin)) = this%border_id(i,j,this%cfg%kmin-1)
                  elseif (find(this%id(i,j,this%cfg%kmin)).ne.find(this%border_id(i,j,this%cfg%kmin-1))) then
                     ! call union(this%border_id(i,j,this%cfg%kmin-1),this%parent(this%id(i,j,this%cfg%kmin)))
                     if (this%is_connected(i,j,this%cfg%kmin,3)) then
                        if (this%border_id(i,j,this%cfg%kmin-1).gt.this%parent(this%id(i,j,this%cfg%kmin))) then
                           call union(this%border_id(i,j,this%cfg%kmin-1),this%parent(this%id(i,j,this%cfg%kmin)))
                        else
                           call union(this%parent(this%id(i,j,this%cfg%kmin)),this%border_id(i,j,this%cfg%kmin-1))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! initialize global stop criterion
      stop_global = 1
      
      ! initialize a counter
      counter = 0
      
      do while (stop_global.ne.0)
         
         ! Initialize local stop flag
         stop_ = 0
         
         this%parent_own = this%parent
         ! Set self-parents to huge(1)
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent(i).eq.i) this%parent(i) = huge(1)
         end do
         
         call MPI_ALLREDUCE(this%parent,this%parent_all,this%cfg%nproc*this%n_struct_max,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
         
         ! Set self-parents back to selves
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent_all(i).eq.huge(1)) this%parent_all(i) = i
            ! if (this%parent_own(i).eq.huge(1)) this%parent_own(i) = i
         end do
         
         ! Flatten trees - is this necessary?
         do i= 1,this%cfg%nproc*this%n_struct_max
            this%parent_all(i) = find_all_2(i,i)
            this%parent_own(i) = find_own(i)
         end do
         
         ! Start with final this%parent array being equal to this%parent_all
         this%parent = this%parent_all
         
         ! Count how many global updates we do
         counter = counter + 1
         
         ! Reconcile conflicts between this%parent_all and this%parent_own
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent_own(i).ne.i) then
               find_parent_own = find(this%parent_own(i))
               find_parent = find(this%parent(i))
               if (find_parent_own.gt.find_parent) then
                  call union(this%parent_own(i),this%parent(i))
                  stop_ = 1
               else if (find_parent.gt.find_parent_own) then
                  call union(this%parent(i),this%parent_own(i))
                  stop_ = 1
               end if
            end if
         end do
         ! Check if we did some changes
         call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         
      end do ! do while (stop_global.ne.0) including domain boundaries
      
      ! Update this%struct_list%parent and point all parents to root
      do i=this%sync_offset+1,this%sync_offset+this%n_struct
         this%struct_list(this%struct_map_(i))%parent = find(this%parent(i))
      end do
      
      ! Update id array with compacted and syncrhonized ids
      do i=this%sync_offset+1,this%sync_offset+this%n_struct
         do n=1,this%struct_list(this%struct_map_(i))%nnode
            ii = this%struct_list(this%struct_map_(i))%node(1,n)
            jj = this%struct_list(this%struct_map_(i))%node(2,n)
            kk = this%struct_list(this%struct_map_(i))%node(3,n)
            this%id(ii,jj,kk) = this%struct_list(this%struct_map_(i))%parent
         end do
      end do
      ! Update ghost cells
      call this%cfg%sync(this%id)
      
      ! Update this%my_struct%id
      ! start marching thru list, starting at first_struct
      this%my_struct => this%first_struct
      
      do while(associated(this%my_struct))
         ii = this%my_struct%node(1,1)
         jj = this%my_struct%node(2,1)
         kk = this%my_struct%node(3,1)
         ! Update structure id with id cell value
         this%my_struct%id  = this%id(ii,jj,kk)
         
         ! Go to next structure
         this%my_struct => this%my_struct%next
         
      end do
      
      return
      
   contains
      
      ! ! For debugging only
      ! function find_only(a_x) result(a_y)
      !    implicit none
      !    integer :: a_x
      !    integer :: a_y
      
      !    a_y = a_x
      !    ! find root of x with path compression
      !    do while ( a_y.ne.this%parent(a_y) )
      !       a_y = this%parent(a_y)
      !    end do
      !    return
      ! end function find_only
      
      ! This function points the this%parent to root and returns that root
      recursive function find(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%parent(a_y) ) then
            this%parent(a_y) = find(this%parent(a_y))
            a_y = this%parent(a_y)
         end if
         return
      end function find
      
      ! Not the same as union in struct_build
      ! This subroutine joins two branches at their roots
      ! At the moment, it joints them at a_y and returns it
      ! but does not return the root
      subroutine union(a_x, a_y)
         implicit none
         integer,intent(in) :: a_x
         integer,intent(in) :: a_y
         
         ! tbd: if a_x tree smaller/lower rank than a_y tree
         this%parent(find(a_x)) = find(a_y)
         
         return
      end subroutine union
      
      ! For this%parent_all array
      ! This function points the this%parent to root and returns that root
      
      recursive function find_all(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%parent_all(a_y) ) then
            this%parent_all(a_y) = find_all(this%parent_all(a_y))
            a_y = this%parent_all(a_y)
         end if
         return
      end function find_all
      
      ! Version that stops at the completion of a cycle
      recursive function find_all_2(a_x,a_starting_node) result(a_y)
         implicit none
         integer :: a_x,a_starting_node
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%parent_all(a_y) ) then
            if ( this%parent_all(a_y).eq.a_starting_node ) then
               a_y = this%parent_all(a_y)
               return
            else
               this%parent_all(a_y) = find_all_2(this%parent_all(a_y),a_starting_node)
               a_y = this%parent_all(a_y)
            end if
         end if
         return
      end function find_all_2
      
      ! For this%parent_own array
      ! This function points the this%parent to root and returns that root
      recursive function find_own(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%parent_own(a_y) ) then
            this%parent_own(a_y) = find_own(this%parent_own(a_y))
            a_y = this%parent_own(a_y)
         end if
         return
      end function find_own
      
      ! This function points the parent to root and returns that root
      recursive function film_find(a_parent,a_x) result(a_y)
         implicit none
         integer, dimension(:) :: a_parent
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.a_parent(a_y) ) then
            a_parent(a_y) = film_find(a_parent,a_parent(a_y))
            a_y = a_parent(a_y)
         end if
         return
      end function film_find
      
      !> Synchronize struct periodicity across all procs
      subroutine struct_sync_per
         implicit none
         
         ! Allocate local and global perodicity arrays
         allocate(this%per_(1:3,this%sync_offset+1:this%sync_offset+this%n_struct_max)); this%per_ = 0
         allocate(this%per (1:3,this%cfg%nproc*this%n_struct_max)); this%per = 0
         
         ! Fill per_ array
         do i=this%sync_offset+1,this%sync_offset+this%n_struct
            this%per_(:,i) = this%struct_list(this%struct_map_(i))%per
         end do
         
         ! Communitcate per
         call MPI_ALLGATHER(this%per_(1,:),this%n_struct_max,MPI_INTEGER,this%per(1,:),this%n_struct_max,MPI_INTEGER,this%cfg%comm,ierr)
         call MPI_ALLGATHER(this%per_(2,:),this%n_struct_max,MPI_INTEGER,this%per(2,:),this%n_struct_max,MPI_INTEGER,this%cfg%comm,ierr)
         call MPI_ALLGATHER(this%per_(3,:),this%n_struct_max,MPI_INTEGER,this%per(3,:),this%n_struct_max,MPI_INTEGER,this%cfg%comm,ierr)
        
         ! Update parent per
         do i=1,this%cfg%nproc*this%n_struct_max
            this%per(:,this%parent(i)) = max(this%per(:,this%parent(i)),this%per(:,i))
         end do
         
         do i=this%sync_offset+1,this%sync_offset+this%n_struct
            do n=1,this%struct_list(this%struct_map_(i))%nnode
               ii = this%struct_list(this%struct_map_(i))%node(1,n)
               jj = this%struct_list(this%struct_map_(i))%node(2,n)
               kk = this%struct_list(this%struct_map_(i))%node(3,n)
               this%idp(:,ii,jj,kk) = this%per(:,this%id(ii,jj,kk))
            end do
         end do
         
         return
         
      end subroutine struct_sync_per
      
   end subroutine struct_sync
   
   
   ! Use PLIC normals to determine if two cells are in the same struct
   function is_connected(this,i,j,k,dim)
      implicit none
      class(ccl), intent(inout) :: this
      integer, intent(in) :: i,j,k,dim
      integer :: ii,jj,kk
      logical :: is_connected, use_normal
      integer, dimension(3) :: pos
      real(WP), dimension(3) :: nref, nloc, pref, ploc
      
      pos = 0
      pos(dim) = -1
      ii = i + pos(1)
      jj = j + pos(2)
      kk = k + pos(3)
      is_connected = .true.
      if (this%max_interface_planes.eq.0) return
      use_normal = getNumberOfVertices(this%poly(1,i,j,k)).gt.0 .and. (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0)
      if (use_normal) then
         ! If neighbor is a two-plane cell
         if (getNumberOfVertices(this%poly(2,ii,jj,kk)).ne.0) then
            nref = calculateNormal(this%poly(1,ii,jj,kk))
            nloc = calculateNormal(this%poly(2,ii,jj,kk))
            pref = calculateCentroid(this%poly(1,ii,jj,kk))
            ploc = calculateCentroid(this%poly(2,ii,jj,kk))
            is_connected = dot_product(ploc-pref,nloc).gt.0.0_WP ! .true. if liquid film
            ! If self is a two-plane cell
         elseif (getNumberOfVertices(this%poly(2,i,j,k)).ne.0) then
            nref = calculateNormal(this%poly(1,i,j,k))
            nloc = calculateNormal(this%poly(2,i,j,k))
            pref = calculateCentroid(this%poly(1,i,j,k))
            ploc = calculateCentroid(this%poly(2,i,j,k))
            is_connected = dot_product(ploc-pref,nloc).gt.0.0_WP ! .true. if liquid film
            ! If both are one-plane cells
         else
            nref = calculateNormal(this%poly(1,i,j,k))
            nloc = calculateNormal(this%poly(1,ii,jj,kk))
            pref = calculateCentroid(this%poly(1,i,j,k))
            ploc = calculateCentroid(this%poly(1,ii,jj,kk))
            is_connected = (dot_product(ploc-pref,nloc).ge.0.0_WP).or.(dot_product(pref-ploc,nref).ge.0.0_WP)
         end if
      end if
      return
   end function is_connected
      
   
   !> Synchronize film labels across all procs
   subroutine film_sync(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
      implicit none
      class(ccl), intent(inout) :: this
      integer :: i,j,k,n,ii,jj,kk,stop_,stop_global,counter,ierr
      integer :: find_parent,find_parent_own
      
      ! Eventually replace with sum(n_struct_per_processor(1:irank)) or border_struct equivalent
      ! Also all "fill parent" loops and allreduces
      allocate(this%film_parent(this%cfg%nproc*this%n_film_max)); this%film_parent = 0
      allocate(this%film_parent_all(this%cfg%nproc*this%n_film_max)); this%film_parent_all = 0
      allocate(this%film_parent_own(this%cfg%nproc*this%n_film_max)); this%film_parent_own = 0
      
      ! Fill parent with selves
      do i= 1,this%cfg%nproc*this%n_film_max
         this%film_parent(i) = i
      end do
      
      ! update ghost cells
      call this%cfg%sync(this%film_border_id)
      
      ! imin_
      if (this%cfg%imin_.ne.this%cfg%imin) then
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%film_border_id(this%cfg%imin_,j,k).gt.0.and.this%film_border_id(this%cfg%imin_-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%film_parent(this%film_border_id(this%cfg%imin_,j,k)).eq.this%film_border_id(this%cfg%imin_,j,k)) then
                     ! call union(this%film_border_id(this%cfg%imin_,j,k),this%film_border_id(this%cfg%imin_-1,j,k))
                     ! is_contiguous = this%film_is_connected(this%cfg%imin_,j,k,1)
                     ! if (is_contiguous) parent(this%film_border_id(this%cfg%imin_,j,k)) = this%film_border_id(this%cfg%imin_-1,j,k)
                     if (this%film_is_connected(this%cfg%imin_,j,k,1)) this%film_parent(this%film_border_id(this%cfg%imin_,j,k)) = this%film_border_id(this%cfg%imin_-1,j,k)
                  elseif (find(this%film_border_id(this%cfg%imin_,j,k)).ne.find(this%film_border_id(this%cfg%imin_-1,j,k))) then
                     ! call union(this%film_border_id(this%cfg%imin_-1,j,k),this%film_parent(this%film_border_id(this%cfg%imin_,j,k)))
                     if (this%film_is_connected(this%cfg%imin_,j,k,1)) then
                        if (this%film_border_id(this%cfg%imin_-1,j,k).gt.this%film_parent(this%film_border_id(this%cfg%imin_,j,k))) then
                           call union(this%film_border_id(this%cfg%imin_-1,j,k),this%film_parent(this%film_border_id(this%cfg%imin_,j,k)))
                        else
                           call union(this%film_parent(this%film_border_id(this%cfg%imin_,j,k)),this%film_border_id(this%cfg%imin_-1,j,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%jmin_
      if (this%cfg%jmin_.ne.this%cfg%jmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%film_border_id(i,this%cfg%jmin_,k).gt.0.and.this%film_border_id(i,this%cfg%jmin_-1,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%film_parent(this%film_border_id(i,this%cfg%jmin_,k)).eq.this%film_border_id(i,this%cfg%jmin_,k)) then
                     ! call union(this%film_border_id(i,this%cfg%jmin_,k),this%film_border_id(i,this%cfg%jmin_-1,k))
                     if (this%film_is_connected(i,this%cfg%jmin_,k,2)) this%film_parent(this%film_border_id(i,this%cfg%jmin_,k)) = this%film_border_id(i,this%cfg%jmin_-1,k)
                  elseif (find(this%film_border_id(i,this%cfg%jmin_,k)).ne.find(this%film_border_id(i,this%cfg%jmin_-1,k))) then
                     ! call union(this%film_border_id(i,this%cfg%jmin_-1,k),this%film_parent(this%film_border_id(i,this%cfg%jmin_,k)))
                     if (this%film_is_connected(i,this%cfg%jmin_,k,2)) then
                        if (this%film_border_id(i,this%cfg%jmin_-1,k).gt.this%film_parent(this%film_border_id(i,this%cfg%jmin_,k))) then
                           call union(this%film_border_id(i,this%cfg%jmin_-1,k),this%film_parent(this%film_border_id(i,this%cfg%jmin_,k)))
                        else
                           call union(this%film_parent(this%film_border_id(i,this%cfg%jmin_,k)),this%film_border_id(i,this%cfg%jmin_-1,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%kmin_
      if (this%cfg%kmin_.ne.this%cfg%kmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               ! If the neighboring ghost cell is filled
               if (this%film_border_id(i,j,this%cfg%kmin_).gt.0.and.this%film_border_id(i,j,this%cfg%kmin_-1).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%film_parent(this%film_border_id(i,j,this%cfg%kmin_)).eq.this%film_border_id(i,j,this%cfg%kmin_)) then
                     ! call union(this%film_border_id(i,j,this%cfg%kmin_),this%film_border_id(i,j,this%cfg%kmin_-1))
                     if (this%film_is_connected(i,j,this%cfg%kmin_,3)) this%film_parent(this%film_border_id(i,j,this%cfg%kmin_)) = this%film_border_id(i,j,this%cfg%kmin_-1)
                  elseif (find(this%film_border_id(i,j,this%cfg%kmin_)).ne.find(this%film_border_id(i,j,this%cfg%kmin_-1))) then
                     ! call union(this%film_border_id(i,j,this%cfg%kmin_-1),this%film_parent(this%film_border_id(i,j,this%cfg%kmin_)))
                     if (this%film_is_connected(i,j,this%cfg%kmin_,3)) then
                        if (this%film_border_id(i,j,this%cfg%kmin_-1).gt.this%film_parent(this%film_border_id(i,j,this%cfg%kmin_))) then
                           call union(this%film_border_id(i,j,this%cfg%kmin_-1),this%film_parent(this%film_border_id(i,j,this%cfg%kmin_)))
                        else
                           call union(this%film_parent(this%film_border_id(i,j,this%cfg%kmin_)),this%film_border_id(i,j,this%cfg%kmin_-1))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! initialize global stop criterion
      stop_global = 1
      
      ! initialize a counter
      counter = 0
      
      do while (stop_global.ne.0)
         
         ! Initialize local stop flag
         stop_ = 0
         
         this%film_parent_own = this%film_parent
         ! Set self-parents to huge(1)
         do i= 1,this%cfg%nproc*this%n_film_max
            if (this%film_parent(i).eq.i) this%film_parent(i) = huge(1)
         end do
         
         call MPI_ALLREDUCE(this%film_parent,this%film_parent_all,this%cfg%nproc*this%n_film_max,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
         
         ! Set self-parents back to selves
         do i= 1,this%cfg%nproc*this%n_film_max
            if (this%film_parent_all(i).eq.huge(1)) this%film_parent_all(i) = i
            ! if (this%film_parent_own(i).eq.huge(1)) this%film_parent_own(i) = i
         end do
         
         ! Flatten trees - is this necessary?
         do i= 1,this%cfg%nproc*this%n_film_max
            this%film_parent_all(i) = find_all(i)
            this%film_parent_own(i) = find_own(i)
         end do
         
         ! Start with final this%film_parent array being equal to this%film_parent_all
         this%film_parent = this%film_parent_all
         
         ! Count how many global updates we do
         counter = counter + 1
         
         ! Reconcile conflicts between this%film_parent_all and this%film_parent_own
         do i= 1,this%cfg%nproc*this%n_film_max
            if (this%film_parent_own(i).ne.i) then
               find_parent_own = find(this%film_parent_own(i))
               find_parent = find(this%film_parent(i))
               if (find_parent_own.gt.find_parent) then
                  call union(this%film_parent_own(i),this%film_parent(i))
                  stop_ = 1
               else if (find_parent.gt.find_parent_own) then
                  call union(this%film_parent(i),this%film_parent_own(i))
                  stop_ = 1
               end if
            end if
         end do
         ! Check if we did some changes
         call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         
      end do ! do while (stop_global.ne.0) excluding domain boundaries
      
      ! Update this%film_list%film_parent and point all parents to root
      do i=this%film_sync_offset+1,this%film_sync_offset+this%n_film
         this%film_list(this%film_map_(i))%parent = find(this%film_parent(i))
      end do
      
      ! Update id array with compacted and syncrhonized ids
      do i=this%film_sync_offset+1,this%film_sync_offset+this%n_film
         do n=1,this%film_list(this%film_map_(i))%nnode
            ii = this%film_list(this%film_map_(i))%node(1,n)
            jj = this%film_list(this%film_map_(i))%node(2,n)
            kk = this%film_list(this%film_map_(i))%node(3,n)
            this%film_id(ii,jj,kk) = this%film_list(this%film_map_(i))%parent
            ! this%film_id(ii,jj,kk) = this%film_list(this%film_id(ii,jj,kk))%parent
         end do
      end do
      
      !   call film_sync_per ! only if we need to keep track of film periodicity
      
      !! Update domain boundaries
      ! this%cfg%imin
      if (this%cfg%imin_.eq.this%cfg%imin) then
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%film_id(this%cfg%imin,j,k).gt.0.and.this%film_border_id(this%cfg%imin-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%film_parent(this%film_id(this%cfg%imin,j,k)).eq.this%film_id(this%cfg%imin,j,k)) then
                     if (this%film_is_connected(this%cfg%imin,j,k,1)) call union(this%film_id(this%cfg%imin,j,k),this%film_border_id(this%cfg%imin-1,j,k))
                     ! this%film_parent(this%film_id(this%cfg%imin,j,k)) = this%film_border_id(this%cfg%imin-1,j,k)
                  elseif (find(this%film_id(this%cfg%imin,j,k)).ne.find(this%film_border_id(this%cfg%imin-1,j,k))) then
                     if (this%film_is_connected(this%cfg%imin,j,k,1)) then
                        if (this%film_border_id(this%cfg%imin-1,j,k).gt.this%film_parent(this%film_id(this%cfg%imin,j,k))) then
                           call union(this%film_border_id(this%cfg%imin-1,j,k),this%film_parent(this%film_id(this%cfg%imin,j,k)))
                        else
                           call union(this%film_parent(this%film_id(this%cfg%imin,j,k)),this%film_border_id(this%cfg%imin-1,j,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%jmin
      if (this%cfg%jmin_.eq.this%cfg%jmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the neighboring ghost cell is filled
               if (this%film_id(i,this%cfg%jmin,k).gt.0.and.this%film_border_id(i,this%cfg%jmin-1,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%film_parent(this%film_id(i,this%cfg%jmin,k)).eq.this%film_id(i,this%cfg%jmin,k)) then
                     if (this%film_is_connected(i,this%cfg%jmin,k,2)) call union(this%film_id(i,this%cfg%jmin,k),this%film_border_id(i,this%cfg%jmin-1,k))
                     ! this%film_parent(this%film_id(i,this%cfg%jmin,k)) = this%film_border_id(i,this%cfg%jmin-1,k)
                  elseif (find(this%film_id(i,this%cfg%jmin,k)).ne.find(this%film_border_id(i,this%cfg%jmin-1,k))) then
                     ! call union(this%film_border_id(i,this%cfg%jmin-1,k),this%film_parent(this%film_id(i,this%cfg%jmin,k)))
                     if (this%film_is_connected(i,this%cfg%jmin,k,2)) then
                        if (this%film_border_id(i,this%cfg%jmin-1,k).gt.this%film_parent(this%film_id(i,this%cfg%jmin,k))) then
                           call union(this%film_border_id(i,this%cfg%jmin-1,k),this%film_parent(this%film_id(i,this%cfg%jmin,k)))
                        else
                           call union(this%film_parent(this%film_id(i,this%cfg%jmin,k)),this%film_border_id(i,this%cfg%jmin-1,k))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! this%cfg%kmin
      if (this%cfg%kmin_.eq.this%cfg%kmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               ! If the neighboring ghost cell is filled
               if (this%film_id(i,j,this%cfg%kmin).gt.0.and.this%film_border_id(i,j,this%cfg%kmin-1).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%film_parent(this%film_id(i,j,this%cfg%kmin)).eq.this%film_id(i,j,this%cfg%kmin)) then
                     if (this%film_is_connected(i,j,this%cfg%kmin,3)) call union(this%film_id(i,j,this%cfg%kmin),this%film_border_id(i,j,this%cfg%kmin-1))
                     ! this%film_parent(this%film_id(i,j,this%cfg%kmin)) = this%film_border_id(i,j,this%cfg%kmin-1)
                  elseif (find(this%film_id(i,j,this%cfg%kmin)).ne.find(this%film_border_id(i,j,this%cfg%kmin-1))) then
                     ! call union(this%film_border_id(i,j,this%cfg%kmin-1),this%film_parent(this%film_id(i,j,this%cfg%kmin)))
                     if (this%film_is_connected(i,j,this%cfg%kmin,3)) then
                        if (this%film_border_id(i,j,this%cfg%kmin-1).gt.this%film_parent(this%film_id(i,j,this%cfg%kmin))) then
                           call union(this%film_border_id(i,j,this%cfg%kmin-1),this%film_parent(this%film_id(i,j,this%cfg%kmin)))
                        else
                           call union(this%film_parent(this%film_id(i,j,this%cfg%kmin)),this%film_border_id(i,j,this%cfg%kmin-1))
                        end if
                     end if
                  end if
               end if
            end do
         end do
      end if
      
      ! initialize global stop criterion
      stop_global = 1
      
      ! initialize a counter
      counter = 0
      
      do while (stop_global.ne.0)
         
         ! Initialize local stop flag
         stop_ = 0
         
         this%film_parent_own = this%film_parent
         ! Set self-parents to huge(1)
         do i= 1,this%cfg%nproc*this%n_film_max
            if (this%film_parent(i).eq.i) this%film_parent(i) = huge(1)
         end do
         
         call MPI_ALLREDUCE(this%film_parent,this%film_parent_all,this%cfg%nproc*this%n_film_max,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
         
         ! Set self-parents back to selves
         do i= 1,this%cfg%nproc*this%n_film_max
            if (this%film_parent_all(i).eq.huge(1)) this%film_parent_all(i) = i
         end do
         
         ! Flatten trees - is this necessary?
         do i= 1,this%cfg%nproc*this%n_film_max
            this%film_parent_all(i) = find_all_2(i,i)
            this%film_parent_own(i) = find_own(i)
         end do
         
         ! Start with final this%film_parent array being equal to this%film_parent_all
         this%film_parent = this%film_parent_all
         
         ! Count how many global updates we do
         counter = counter + 1
         
         ! Reconcile conflicts between this%film_parent_all and this%film_parent_own
         do i= 1,this%cfg%nproc*this%n_film_max
            if (this%film_parent_own(i).ne.i) then
               find_parent_own = find(this%film_parent_own(i))
               find_parent = find(this%film_parent(i))
               if (find_parent_own.gt.find_parent) then
                  call union(this%film_parent_own(i),this%film_parent(i))
                  stop_ = 1
               else if (find_parent.gt.find_parent_own) then
                  call union(this%film_parent(i),this%film_parent_own(i))
                  stop_ = 1
               end if
            end if
         end do
         ! Check if we did some changes
         call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         
      end do ! do while (stop_global.ne.0) including domain boundaries
      
      ! Update this%film_list%film_parent and point all parents to root
      do i=this%film_sync_offset+1,this%film_sync_offset+this%n_film
         this%film_list(this%film_map_(i))%parent = find(this%film_parent(i))
      end do
      
      ! Update id array with compacted and syncrhonized ids
      do i=this%film_sync_offset+1,this%film_sync_offset+this%n_film
         do n=1,this%film_list(this%film_map_(i))%nnode
            ii = this%film_list(this%film_map_(i))%node(1,n)
            jj = this%film_list(this%film_map_(i))%node(2,n)
            kk = this%film_list(this%film_map_(i))%node(3,n)
            this%film_id(ii,jj,kk) = this%film_list(this%film_map_(i))%parent
         end do
      end do
      ! Update ghost cells
      call this%cfg%sync(this%film_id)
      
      return
      
   contains
      
      ! ! For debugging only
      ! function find_only(a_x) result(a_y)
      !    implicit none
      !    integer :: a_x
      !    integer :: a_y
      
      !    a_y = a_x
      !    ! find root of x with path compression
      !    do while ( a_y.ne.this%film_parent(a_y) )
      !       a_y = this%film_parent(a_y)
      !    end do
      !    return
      ! end function find_only
      
      ! This function points the this%film_parent to root and returns that root
      recursive function find(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%film_parent(a_y) ) then
            this%film_parent(a_y) = find(this%film_parent(a_y))
            a_y = this%film_parent(a_y)
         end if
         return
      end function find
      
      ! Not the same as union in struct_build
      ! This subroutine joins two branches at their roots
      ! At the moment, it joints them at a_y and returns it
      ! but does not return the root
      subroutine union(a_x, a_y)
         implicit none
         integer,intent(in) :: a_x
         integer,intent(in) :: a_y
         
         ! tbd: if a_x tree smaller/lower rank than a_y tree
         this%film_parent(find(a_x)) = find(a_y)
         
         return
      end subroutine union
      
      ! For this%film_parent_all array
      ! This function points the this%film_parent to root and returns that root
      recursive function find_all(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%film_parent_all(a_y) ) then
            this%film_parent_all(a_y) = find_all(this%film_parent_all(a_y))
            a_y = this%film_parent_all(a_y)
         end if
         return
      end function find_all
      
      ! Version that stops at the completion of a cycle
      recursive function find_all_2(a_x,a_starting_node) result(a_y)
         implicit none
         integer :: a_x,a_starting_node
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%film_parent_all(a_y) ) then
            if ( this%film_parent_all(a_y).eq.a_starting_node ) then
               a_y = this%film_parent_all(a_y)
               return
            else
               this%film_parent_all(a_y) = find_all_2(this%film_parent_all(a_y),a_starting_node)
               a_y = this%film_parent_all(a_y)
            end if
         end if
         return
      end function find_all_2
      
      ! For this%film_parent_own array
      ! This function points the this%film_parent to root and returns that root
      recursive function find_own(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.this%film_parent_own(a_y) ) then
            this%film_parent_own(a_y) = find_own(this%film_parent_own(a_y))
            a_y = this%film_parent_own(a_y)
         end if
         return
      end function find_own
      
      ! This function points the parent to root and returns that root
      recursive function film_find(a_parent,a_x) result(a_y)
         implicit none
         integer, dimension(:) :: a_parent
         integer :: a_x
         integer :: a_y
         
         a_y = a_x
         ! find root of x with path compression
         if ( a_y.ne.a_parent(a_y) ) then
            a_parent(a_y) = film_find(a_parent,a_parent(a_y))
            a_y = a_parent(a_y)
         end if
         return
      end function film_find
      
      ! ! Synchronize struct periodicity across all procs
      ! subroutine film_sync_per
      !    implicit none
      
      !    ! Allocate local and global perodicity arrays
      !    allocate(this%per_(1:3,this%film_sync_offset+1:this%film_sync_offset+this%n_film_max)); this%per_ = 0
      !    allocate(this%per (1:3,this%cfg%nproc*this%n_film_max)); this%per = 0
      
      !    ! Fill per_ array
      !    do i=this%film_sync_offset+1,this%film_sync_offset+this%n_film
      !       this%per_(:,i) = this%film_list(this%film_map_(i))%per
      !    end do
      
      !    ! Communitcate per
      !    call MPI_ALLGATHER(this%per_(1,:),this%n_film_max,MPI_INTEGER,this%per(1,:),this%n_film_max,MPI_INTEGER,this%cfg%comm,ierr)
      !    call MPI_ALLGATHER(this%per_(2,:),this%n_film_max,MPI_INTEGER,this%per(2,:),this%n_film_max,MPI_INTEGER,this%cfg%comm,ierr)
      !    call MPI_ALLGATHER(this%per_(3,:),this%n_film_max,MPI_INTEGER,this%per(3,:),this%n_film_max,MPI_INTEGER,this%cfg%comm,ierr)
      !    ! Update parent per
      !    do i=1,this%cfg%nproc*this%n_film_max
      !       this%per(:,this%film_parent(i)) = max(this%per(:,this%film_parent(i)),this%per(:,i))
      !    end do
      
      !    do i=this%film_sync_offset+1,this%film_sync_offset+this%n_film
      !       do n=1,this%film_list(this%film_map_(i))%nnode
      !          ii = this%film_list(this%film_map_(i))%node(1,n)
      !          jj = this%film_list(this%film_map_(i))%node(2,n)
      !          kk = this%film_list(this%film_map_(i))%node(3,n)
      !          this%film_idp(:,ii,jj,kk) = this%per(:,this%film_id(ii,jj,kk))
      !       end do
      !    end do
      
      !    return
      
      ! end subroutine film_sync_per
      
   end subroutine film_sync
   
   
   ! Use PLIC normals to determine if two cells are in the same film
   ! Fow now, always connected
   function film_is_connected(this,i,j,k,dim)
      implicit none
      class(ccl), intent(inout) :: this
      integer, intent(in) :: i,j,k,dim
      integer :: ii,jj,kk
      logical :: film_is_connected!, use_normal
      integer, dimension(3) :: pos
      ! real(WP), dimension(3) :: nref, nloc, pref, ploc
      
      ! pos = 0
      ! pos(dim) = -1
      ! ii = i + pos(1)
      ! jj = j + pos(2)
      ! kk = k + pos(3)
      film_is_connected = .true.
      ! use_normal = getNumberOfVertices(this%poly(1,i,j,k)).gt.0 .and. (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0)
      ! if (use_normal) then
      !    ! If neighbor is a two-plane cell
      !    if (getNumberOfVertices(this%poly(2,ii,jj,kk)).ne.0) then
      !       nref = calculateNormal(this%poly(1,ii,jj,kk))
      !       nloc = calculateNormal(this%poly(2,ii,jj,kk))
      !       pref = calculateCentroid(this%poly(1,ii,jj,kk))
      !       ploc = calculateCentroid(this%poly(2,ii,jj,kk))
   !       film_is_connected = dot_product(ploc-pref,nloc).gt.0.0_WP ! .true. if liquid film
      !    ! If self is a two-plane cell
      !    elseif (getNumberOfVertices(this%poly(2,i,j,k)).ne.0) then
      !       nref = calculateNormal(this%poly(1,i,j,k))
      !       nloc = calculateNormal(this%poly(2,i,j,k))
      !       pref = calculateCentroid(this%poly(1,i,j,k))
      !       ploc = calculateCentroid(this%poly(2,i,j,k))
   !       film_is_connected = dot_product(ploc-pref,nloc).gt.0.0_WP ! .true. if liquid film
      !    ! If both are one-plane cells
      !    else
      !       nref = calculateNormal(this%poly(1,i,j,k))
      !       nloc = calculateNormal(this%poly(1,ii,jj,kk))
      !       pref = calculateCentroid(this%poly(1,i,j,k))
      !       ploc = calculateCentroid(this%poly(1,ii,jj,kk))
   !       film_is_connected = (dot_product(ploc-pref,nloc).ge.0.0_WP).or.(dot_product(pref-ploc,nref).ge.0.0_WP)
      !    end if
      ! end if
      return
   end function film_is_connected
   
   
   !> Final update of struct labels in global array
   subroutine struct_label_update(this)
      implicit none
      class(ccl), intent(inout) :: this
      integer :: i,ii,jj,kk,n_node,tag_final
      
      ! Start at first_struct
      this%my_struct => this%first_struct
      
      ! init local vars
      n_node = 0
      tag_final = 0
      
      ! ! re-init global tag array
      ! this%id = 0
      
      ! loop through list
      do while(associated(this%my_struct))
         n_node = this%my_struct%nnode
         ! tag_final = this%my_struct%id
         ! loop over nodes
         do i=1,n_node
            ! Update tag
            ii= this%my_struct%node(1,i)
            jj= this%my_struct%node(2,i)
            kk= this%my_struct%node(3,i)
            ! this%id(ii,jj,kk) = tag_final
            ! Update periodicity
            this%my_struct%per(1) = max(this%my_struct%per(1),this%idp(1,ii,jj,kk))
            this%my_struct%per(2) = max(this%my_struct%per(2),this%idp(2,ii,jj,kk))
            this%my_struct%per(3) = max(this%my_struct%per(3),this%idp(3,ii,jj,kk))
         end do
         
         this%my_struct => this%my_struct%next
         
      end do
      
      return
   end subroutine struct_label_update
   
   
   !> Build array of meta_structures
   subroutine struct_final(this,U,V,W)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      implicit none
      class(ccl), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,ierr
      integer, dimension(:), pointer :: buf => null()
      
      ! Allocate global meta_structures array
      this%n_meta_struct = this%n_struct_max * this%cfg%nproc
      allocate(buf(1:this%n_meta_struct),this%meta_structures(1:this%n_meta_struct))
      
      ! Point to head element
      this%my_struct => this%first_struct
      
      ! Offset to part of array belonging to local proc
      i = this%cfg%rank*this%n_struct_max + 1
      
      ! Gather id of local structs
      buf = 0
      do while(associated(this%my_struct))
         buf(i) = this%my_struct%id
         this%my_struct => this%my_struct%next
         i = i+1
      end do
      
      ! Collect data
      call MPI_ALLREDUCE(buf,this%meta_structures,this%n_meta_struct,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
      
      ! Sort, purge to single list of unique ID tags
      call meta_structures_sort(this)
      
      ! Compute stats
      call meta_structures_stats(this,U,V,W)
      
      ! Clean up
      deallocate(buf)
      nullify(buf)
      
      return
   end subroutine struct_final
   
   
   !> Sort list of ID tags, purge extra elements, eliminate duplicates
   subroutine meta_structures_sort(this)
      implicit none
      class(ccl), intent(inout) :: this
      integer, dimension(:), pointer :: meta_structures_tmp => null()
      integer :: i,n_meta_struct_tmp
      
      ! Sort first
      call quick_sort(this%meta_structures)
      
      ! Compact list
      n_meta_struct_tmp = 1
      i = 1
      do while (i.le.this%n_meta_struct)
         if (this%meta_structures(i).eq.0) then
            i = i+1
            this%meta_structures(n_meta_struct_tmp) = this%meta_structures(i)
         elseif (this%meta_structures(i).gt.this%meta_structures(n_meta_struct_tmp)) then
            n_meta_struct_tmp = n_meta_struct_tmp + 1
            this%meta_structures(n_meta_struct_tmp) = this%meta_structures(i)
            i = i+1
         else
            i = i+1
         end if
      end do
      
      ! Resize this%meta_structures array
      this%n_meta_struct = n_meta_struct_tmp
      allocate(meta_structures_tmp(this%n_meta_struct))
      meta_structures_tmp = this%meta_structures(1:this%n_meta_struct)
      deallocate(this%meta_structures)
      nullify(this%meta_structures)
      this%meta_structures => meta_structures_tmp
      return
      
   contains
      
      recursive subroutine quick_sort(A)
         implicit none
         integer, dimension(:)  :: A
         integer :: imark
         
         if (size(A).gt.1) then
            call qs_partition(A,imark)
            call quick_sort(A(:imark-1))
            call quick_sort(A(imark:))
         end if
         
         return
      end subroutine quick_sort
      
      subroutine qs_partition(A,marker)
         implicit none
         integer,  dimension(:) :: A
         integer, intent(out) :: marker
         integer :: ii,jj
         integer :: itmp
         integer :: x
         
         x = A(1)
         ii = 0
         jj = size(A) + 1
         
         do
            jj = jj-1
            do
               if (A(jj).le.x) exit
               jj = jj-1
            end do
            ii = ii+1
            do
               if (A(ii).ge.x) exit
               ii = ii+1
            end do
            if (ii.lt.jj) then
               ! Exchange A(ii) and A(jj)
               itmp = A(ii)
               A(ii) = A(jj)
               A(jj) = itmp
            else if (ii.eq.jj) then
               marker = ii+1
               return
            else
               marker = ii
               return
            endif
         end do
         
         return
      end subroutine qs_partition
   end subroutine meta_structures_sort
   
   
   !> Compute stats for full meta-structures
   subroutine meta_structures_stats(this,U,V,W)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(ccl), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(:), pointer :: x_cg,y_cg,z_cg
      real(WP), dimension(:), pointer :: vol_struct_,vol_struct
      real(WP), dimension(:), pointer :: x_vol_,x_vol,y_vol_,y_vol,z_vol_,z_vol
      real(WP), dimension(:), pointer :: u_vol_,u_vol,v_vol,v_vol_,w_vol,w_vol_
      real(WP), dimension(:,:,:), pointer :: Imom_,Imom
      integer  :: i,j,ii,jj,kk,ierr!,iunit,var
      integer  :: per_x,per_y,per_z
      real(WP) :: xtmp,ytmp,ztmp
      ! character(len=str_medium) :: filename
      ! Eigenvalues/eigenvectors
      real(WP), dimension(3,3) :: A
      real(WP), dimension(3) :: d
      integer , parameter :: order = 3
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      integer  :: lwork,info,n
      
      ! Query optimal work array size
      call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))

      ! allocate / initialize temps arrays for computation
      allocate(vol_struct(1:this%n_meta_struct),vol_struct_(1:this%n_meta_struct))
      allocate(x_vol(1:this%n_meta_struct),x_vol_(1:this%n_meta_struct))
      allocate(y_vol(1:this%n_meta_struct),y_vol_(1:this%n_meta_struct))
      allocate(z_vol(1:this%n_meta_struct),z_vol_(1:this%n_meta_struct))
      allocate(u_vol(1:this%n_meta_struct),u_vol_(1:this%n_meta_struct))
      allocate(v_vol(1:this%n_meta_struct),v_vol_(1:this%n_meta_struct))
      allocate(w_vol(1:this%n_meta_struct),w_vol_(1:this%n_meta_struct))
      allocate(Imom(1:this%n_meta_struct,3,3),Imom_(1:this%n_meta_struct,3,3))
      allocate(x_cg(1:this%n_meta_struct),y_cg(1:this%n_meta_struct),z_cg(1:this%n_meta_struct))
      
      ! initialize temp arrays for computation
      vol_struct(:) = 0.0_WP;vol_struct_(:)= 0.0_WP
      x_vol(:) = 0.0_WP;x_vol_(:) = 0.0_WP
      y_vol(:) = 0.0_WP;y_vol_(:) = 0.0_WP
      z_vol(:) = 0.0_WP;z_vol_(:) = 0.0_WP
      u_vol(:) = 0.0_WP;u_vol_(:) = 0.0_WP
      v_vol(:) = 0.0_WP;v_vol_(:) = 0.0_WP
      w_vol(:) = 0.0_WP;w_vol_(:) = 0.0_WP
      Imom(:,:,:) = 0.0_WP;Imom_(:,:,:) = 0.0_WP;
      
      ! fill in data (except moments of inertia)
      do i = 1,this%n_meta_struct
         this%my_struct => this%first_struct
         do while(associated(this%my_struct))
            if (this%my_struct%id.eq.this%meta_structures(i)) then
               
               ! Periodicity
               per_x = this%my_struct%per(1)
               per_y = this%my_struct%per(2)
               per_z = this%my_struct%per(3)
               
               do j = 1,this%my_struct%nnode
                  
                  ! Indices of struct node
                  ii = this%my_struct%node(1,j)
                  jj = this%my_struct%node(2,j)
                  kk = this%my_struct%node(3,j)
                  
                  ! Location of struct node
                  xtmp = this%cfg%xm(ii)-per_x*this%cfg%xL
                  ytmp = this%cfg%ym(jj)-per_y*this%cfg%yL
                  ztmp = this%cfg%zm(kk)-per_z*this%cfg%zL
                  
                  ! Volume
                  vol_struct_(i) = vol_struct_(i) + this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  
                  ! Center of gravity
                  x_vol_(i) = x_vol_(i) + xtmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  y_vol_(i) = y_vol_(i) + ytmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  z_vol_(i) = z_vol_(i) + ztmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  
                  ! Average gas velocity inside struct
                  u_vol_(i) = u_vol_(i) + U(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  v_vol_(i) = v_vol_(i) + V(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  w_vol_(i) = w_vol_(i) + W(ii,jj,kk)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  
               end do
            end if
            this%my_struct => this%my_struct%next
         end do
      end do
      
      ! Sum parallel stats
      call MPI_ALLREDUCE(vol_struct_,vol_struct,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(x_vol_,x_vol,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_vol_,y_vol,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_vol_,z_vol,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(u_vol_,u_vol,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(v_vol_,v_vol,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(w_vol_,w_vol,this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      
      ! fill in moments of inertia
      do i = 1,this%n_meta_struct
         this%my_struct => this%first_struct
         do while(associated(this%my_struct))
            if (this%my_struct%id.eq.this%meta_structures(i)) then
               
               ! Periodicity
               per_x = this%my_struct%per(1)
               per_y = this%my_struct%per(2)
               per_z = this%my_struct%per(3)
               
               do j = 1,this%my_struct%nnode
                  
                  ! Indices of struct node
                  ii = this%my_struct%node(1,j)
                  jj = this%my_struct%node(2,j)
                  kk = this%my_struct%node(3,j)
                  
                  ! Location of struct node
                  xtmp = this%cfg%xm(ii)-per_x*this%cfg%xL-x_vol(i)/vol_struct(i)
                  ytmp = this%cfg%ym(jj)-per_y*this%cfg%yL-y_vol(i)/vol_struct(i)
                  ztmp = this%cfg%zm(kk)-per_z*this%cfg%zL-z_vol(i)/vol_struct(i)
                  
                  ! Moment of Inertia
                  Imom_(i,1,1) = Imom_(i,1,1) + (ytmp**2 + ztmp**2)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  Imom_(i,2,2) = Imom_(i,2,2) + (xtmp**2 + ztmp**2)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  Imom_(i,3,3) = Imom_(i,3,3) + (xtmp**2 + ytmp**2)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  
                  Imom_(i,1,2) = Imom_(i,1,2) - xtmp*ytmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  Imom_(i,1,3) = Imom_(i,1,3) - xtmp*ztmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  Imom_(i,2,3) = Imom_(i,2,3) - ytmp*ztmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                  
               end do
            end if
            this%my_struct => this%my_struct%next
         end do
      end do
      
      ! Sum parallel stat on Imom
      do i=1,3
         do j=1,3
            call MPI_ALLREDUCE(Imom_(:,i,j),Imom(:,i,j),this%n_meta_struct,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         end do
      end do
      
      ! Allocate and initalize final list of meta_structures
      allocate(this%meta_structures_list(1:this%n_meta_struct))
      this%meta_structures_list(:)%id = 0
      this%meta_structures_list(:)%vol= 0.0_WP
      this%meta_structures_list(:)%x  = 0.0_WP
      this%meta_structures_list(:)%y  = 0.0_WP
      this%meta_structures_list(:)%z  = 0.0_WP
      this%meta_structures_list(:)%u  = 0.0_WP
      this%meta_structures_list(:)%v  = 0.0_WP
      this%meta_structures_list(:)%w  = 0.0_WP
      
      ! Store data
      do i=1,this%n_meta_struct
         
         ! Center of gravity
         x_cg(i) = x_vol(i)/vol_struct(i)
         y_cg(i) = y_vol(i)/vol_struct(i)
         z_cg(i) = z_vol(i)/vol_struct(i)
         
         ! Periodicity: transport back inside domain if needed
         if (x_cg(i).lt.this%cfg%x(this%cfg%imin)) x_cg(i) = x_cg(i)+this%cfg%xL
         if (y_cg(i).lt.this%cfg%y(this%cfg%jmin)) y_cg(i) = y_cg(i)+this%cfg%yL
         if (z_cg(i).lt.this%cfg%z(this%cfg%kmin)) z_cg(i) = z_cg(i)+this%cfg%zL
         
         ! tag, volume, location and velocity
         this%meta_structures_list(i)%id = this%meta_structures(i)
         this%meta_structures_list(i)%vol= vol_struct(i)
         this%meta_structures_list(i)%x  = x_cg(i)
         this%meta_structures_list(i)%y  = y_cg(i)
         this%meta_structures_list(i)%z  = z_cg(i)
         this%meta_structures_list(i)%u  = u_vol(i)/vol_struct(i)
         this%meta_structures_list(i)%v  = v_vol(i)/vol_struct(i)
         this%meta_structures_list(i)%w  = w_vol(i)/vol_struct(i)
         
         ! Eigenvalues/eigenvectors of moments of inertia tensor
         A = Imom(i,:,:)
         n = 3
         ! On exit, A contains eigenvectors, and d contains eigenvalues in ascending order
         call dsyev('V','U',n,A,n,d,work,lwork,info)
         ! Get rid of very small negative values (due to machine accuracy)
         d = max(0.0_WP,d)
         ! Store characteristic lengths
         this%meta_structures_list(i)%lengths(1) = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol_struct(i))
         this%meta_structures_list(i)%lengths(2) = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol_struct(i))
         this%meta_structures_list(i)%lengths(3) = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol_struct(i))
         ! Zero out length in 3rd dimension if 2D
         if (this%cfg%nx.eq.1.or.this%cfg%ny.eq.1.or.this%cfg%nz.eq.1) this%meta_structures_list(i)%lengths(3)=0.0_WP
         ! Store principal axes
         this%meta_structures_list(i)%axes(:,:) = A
         
      end do
      
      ! Write output
      !if (this%cfg%amRoot) then
      !
      !   ! Open file with correct index
      !   filename = "struct_stat."
      !   write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%nout_time
      !   open(iunit,file=filename,form="formatted",iostat=ierr,status="REPLACE")
      !
      !   ! Header
      !   write(iunit,'(10000a20)') (trim(adjustl(this%meta_structures_name(var))),var=1,this%meta_structures_nname)
      !
      !   ! Data
      !   do i=1,this%n_meta_struct
      !      write(iunit,'(I20,10000ES20.12)')       this%meta_structures_list(i)%id,        this%meta_structures_list(i)%vol,       &
      !      this%meta_structures_list(i)%x,         this%meta_structures_list(i)%y,         this%meta_structures_list(i)%z,         &
      !      this%meta_structures_list(i)%u,         this%meta_structures_list(i)%v,         this%meta_structures_list(i)%w,         &
      !      this%meta_structures_list(i)%lengths(1),this%meta_structures_list(i)%lengths(2),this%meta_structures_list(i)%lengths(3),&
      !      this%meta_structures_list(i)%axes(1,1), this%meta_structures_list(i)%axes(1,2), this%meta_structures_list(i)%axes(1,3), &
      !      this%meta_structures_list(i)%axes(2,1), this%meta_structures_list(i)%axes(2,2), this%meta_structures_list(i)%axes(2,3), &
      !      this%meta_structures_list(i)%axes(3,1), this%meta_structures_list(i)%axes(3,2), this%meta_structures_list(i)%axes(3,3)
      !   end do
      !   close(iunit)
      !end if
      
      ! Deallocate arrays
      deallocate(vol_struct_,x_vol_,y_vol_,z_vol_,u_vol_,v_vol_,w_vol_,Imom_)
      deallocate(vol_struct,x_vol,y_vol,z_vol,u_vol,v_vol,w_vol,Imom)
      deallocate(x_cg,y_cg,z_cg)
      nullify(vol_struct_,x_vol_,y_vol_,z_vol_,u_vol_,v_vol_,w_vol_,Imom_)
      nullify(vol_struct,x_vol,y_vol,z_vol,u_vol,v_vol,w_vol,Imom)
      nullify(x_cg,y_cg,z_cg)
      
      return
   end subroutine meta_structures_stats
   

   !> Classify film by shape
   subroutine film_classify(this,Lbary,Gbary)
      implicit none
      class(ccl), intent(inout) :: this
      real(WP), dimension(:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Lbary  !< Liquid barycenter
      real(WP), dimension(:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Gbary  !< Gas barycenter
      integer  :: m,n,i,j,k,ii,jj,kk
      real(WP) :: xtmp, ytmp, ztmp
      real(WP) :: x_vol,y_vol,z_vol,vol_total
      ! Local moment of inertia tensor
      real(WP), dimension(3,3) :: Imom
      ! Eigenvalues/eigenvectors
      real(WP), dimension(3)     :: d
      integer , parameter        :: order = 3
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      integer  :: lwork,info
      ! ! Characteristic lengths
      ! real(WP) :: l1,l2,l3
      real(WP), parameter :: ratio = 1.5_WP
      ! Edge detection
      real(WP), dimension(3) :: c_local, c_filter

      ! Query optimal work array size
      call dsyev('V','U',order,Imom,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))

      this%film_type = 0
      this%film_edge = 0.0_WP
      this%film_edge_normal = 0.0_WP
      do m=this%film_sync_offset+1,this%film_sync_offset+this%n_film ! Loops over film segments contained locally
         if (this%film_list(this%film_map_(m))%phase.eq.1) then ! Liquid film
            do n=1,this%film_list(this%film_map_(m))%nnode ! Loops over cells within local film segment
               i = this%film_list(this%film_map_(m))%node(1,n)
               j = this%film_list(this%film_map_(m))%node(2,n)
               k = this%film_list(this%film_map_(m))%node(3,n)
               vol_total = 0.0_WP
               x_vol = 0.0_WP; y_vol = 0.0_WP; z_vol = 0.0_WP
               Imom(:,:) = 0.0_WP
               do kk = k-2,k+2
                  do jj = j-2,j+2
                     do ii = i-2,i+2
                        
                        ! Location of film node
                        xtmp = Lbary(1,ii,jj,kk)
                        ytmp = Lbary(2,ii,jj,kk)
                        ztmp = Lbary(3,ii,jj,kk)
                        
                        ! Volume
                        vol_total = vol_total + this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        
                        ! Center of gravity
                        x_vol = x_vol + xtmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        y_vol = y_vol + ytmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        z_vol = z_vol + ztmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        
                     end do
                  end do
               end do
               do kk = k-2,k+2
                  do jj = j-2,j+2
                     do ii = i-2,i+2
                        
                        ! Location of film node
                        xtmp = Lbary(1,ii,jj,kk)-x_vol/vol_total
                        ytmp = Lbary(2,ii,jj,kk)-y_vol/vol_total
                        ztmp = Lbary(3,ii,jj,kk)-z_vol/vol_total
                        
                        Imom(1,1) = Imom(1,1) + (ytmp**2 + ztmp**2)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        Imom(2,2) = Imom(2,2) + (xtmp**2 + ztmp**2)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        Imom(3,3) = Imom(3,3) + (xtmp**2 + ytmp**2)*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        
                        Imom(1,2) = Imom(1,2) - xtmp*ytmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        Imom(1,3) = Imom(1,3) - xtmp*ztmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        Imom(2,3) = Imom(2,3) - ytmp*ztmp*this%cfg%vol(ii,jj,kk)*this%VF(ii,jj,kk)
                        
                     end do
                  end do
               end do
               ! On exit, Imom contains eigenvectors, and d contains eigenvalues in ascending order
               call dsyev('V','U',order,Imom,order,d,work,lwork,info)
               ! Get rid of very small negative values (due to machine accuracy)
               d = max(0.0_WP,d)
               ! ! Calculate characteristic lengths assuming ellipsoid
               ! l1 = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol_total)
               ! l2 = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol_total)
               ! l3 = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol_total)
               this%film_type(i,j,k) = 0
               if (d(3).gt.(ratio*d(1))) this%film_type(i,j,k) = this%film_type(i,j,k) + 1
               if (d(3).gt.(ratio*d(2))) this%film_type(i,j,k) = this%film_type(i,j,k) + 1
               ! if (l1.gt.(ratio*l2)) film_type(i,j,k) = film_type(i,j,k) + 1
               ! if (l1.gt.(ratio*l3)) film_type(i,j,k) = film_type(i,j,k) + 1

               ! We get edge detection for free - need to distinguish between sheet and ligament edge
               ! this%film_edge(i,j,k) = l1/l2
               this%film_edge(i,j,k) = d(2)/d(1)

               ! Calculate edge normal
               c_local = Lbary(:,i,j,k)
               c_filter = [x_vol,y_vol,z_vol]/vol_total
               this%film_edge_normal(:,i,j,k) = (c_local-c_filter)/norm2(c_local-c_filter)
            end do
         else ! Gas film
            do n=1,this%film_list(this%film_map_(m))%nnode
               i = this%film_list(this%film_map_(m))%node(1,n)
               j = this%film_list(this%film_map_(m))%node(2,n)
               k = this%film_list(this%film_map_(m))%node(3,n)
               vol_total = 0.0_WP
               x_vol = 0.0_WP; y_vol = 0.0_WP; z_vol = 0.0_WP
               Imom(:,:) = 0.0_WP
               do kk = k-2,k+2
                  do jj = j-2,j+2
                     do ii = i-2,i+2
                        
                        ! Location of film node
                        xtmp = Gbary(1,ii,jj,kk)
                        ytmp = Gbary(2,ii,jj,kk)
                        ztmp = Gbary(3,ii,jj,kk)

                        ! Volume
                        vol_total = vol_total + this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))

                        ! Center of gravity
                        x_vol = x_vol + xtmp*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        y_vol = y_vol + ytmp*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        z_vol = z_vol + ztmp*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
            
                     end do
                  end do
               end do
               do kk = k-2,k+2
                  do jj = j-2,j+2
                     do ii = i-2,i+2
                        
                        ! Location of film node
                        xtmp = Gbary(1,ii,jj,kk)-x_vol/vol_total
                        ytmp = Gbary(2,ii,jj,kk)-y_vol/vol_total
                        ztmp = Gbary(3,ii,jj,kk)-z_vol/vol_total

                        Imom(1,1) = Imom(1,1) + (ytmp**2 + ztmp**2)*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        Imom(2,2) = Imom(2,2) + (xtmp**2 + ztmp**2)*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        Imom(3,3) = Imom(3,3) + (xtmp**2 + ytmp**2)*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
      
                        Imom(1,2) = Imom(1,2) - xtmp*ytmp*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        Imom(1,3) = Imom(1,3) - xtmp*ztmp*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        Imom(2,3) = Imom(2,3) - ytmp*ztmp*this%cfg%vol(ii,jj,kk)*(1.0_WP-this%VF(ii,jj,kk))
                        
                     end do
                  end do
               end do
               ! On exit, Imom contains eigenvectors, and d contains eigenvalues in ascending order
               call dsyev('V','U',order,Imom,order,d,work,lwork,info)
               ! Get rid of very small negative values (due to machine accuracy)
               d = max(0.0_WP,d)
               ! ! Calculate characteristic lengths assuming ellipsoid
               ! l1 = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol_total)
               ! l2 = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol_total)
               ! l3 = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol_total)
               this%film_type(i,j,k) = 0
               if (d(3).gt.(ratio*d(1))) this%film_type(i,j,k) = this%film_type(i,j,k) + 1
               if (d(3).gt.(ratio*d(2))) this%film_type(i,j,k) = this%film_type(i,j,k) + 1
               ! if (l1.gt.(ratio*l2)) film_type(i,j,k) = film_type(i,j,k) + 1
               ! if (l1.gt.(ratio*l3)) film_type(i,j,k) = film_type(i,j,k) + 1
            end do
         end if ! Film phase
      end do
      ! Update ghost cells
      call this%cfg%sync(this%film_type)
      call this%cfg%sync(this%film_edge)
      call this%cfg%sync(this%film_edge_normal)

   end subroutine film_classify
   
   
   !> Find the minimum thickness
   subroutine get_min_thickness(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(ccl), intent(inout) :: this
      integer  :: id,m,n,i,j,k,ierr
      real(WP), dimension(1:this%cfg%nproc*this%n_film_max) :: min_thickness_,min_thickness

      if (this%n_film_max.eq.0) return ! If there are no films globally
      min_thickness_=huge(1.0_WP)
      do m=this%film_sync_offset+1,this%film_sync_offset+this%n_film ! Loops over film segments contained locally
         id=this%film_list(this%film_map_(m))%parent
         do n=1,this%film_list(this%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%film_list(this%film_map_(m))%node(1,n)
            j=this%film_list(this%film_map_(m))%node(2,n)
            k=this%film_list(this%film_map_(m))%node(3,n)
            min_thickness_(id)=min(min_thickness_(id),this%film_thickness(i,j,k))
         end do
      end do
      call MPI_ALLREDUCE(min_thickness_,min_thickness,this%cfg%nproc*this%n_film_max,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      do m=this%film_sync_offset+1,this%film_sync_offset+this%n_film ! Loops over film segments contained locally
         id=this%film_list(this%film_map_(m))%parent
         this%film_list(this%film_map_(m))%min_thickness = min_thickness(id)
      end do
   end subroutine get_min_thickness
   
   
   !> Sort film indices by increasing thickness
   subroutine sort_by_thickness(this)
      implicit none
      class(ccl), intent(inout) :: this
      integer  :: m,n
      real(WP), dimension(  :), allocatable :: mythick
      ! Nothing to do if there is no local film
      if (this%n_film.eq.0) return
      ! Build an auxiliary array with thickness for sorting for each film
      do m=this%film_sync_offset+1,this%film_sync_offset+this%n_film
         ! Allocate and copy 1D thickness and node map
         allocate(mythick(    this%film_list(this%film_map_(m))%nnode))
         do n=1,this%film_list(this%film_map_(m))%nnode
            mythick(n)=this%film_thickness(this%film_list(this%film_map_(m))%node(1,n),&
            &                              this%film_list(this%film_map_(m))%node(2,n),&
            &                              this%film_list(this%film_map_(m))%node(3,n))
         end do
         ! Sort it
         call quick_sort_by_thickness(mythick,this%film_list(this%film_map_(m))%node)
         ! Deallocate thickness and map
         deallocate(mythick)
      end do
   contains
      ! Thickness sorting
      recursive subroutine quick_sort_by_thickness(thick,nodes)
         implicit none
         real(WP), dimension(:)   :: thick
         integer , dimension(:,:) :: nodes
         integer :: imark
         if (size(thick).gt.1) then
            call quick_sort_by_thickness_partition(thick,nodes,imark)
            call quick_sort_by_thickness(thick(     :imark-1),nodes(:,     :imark-1))
            call quick_sort_by_thickness(thick(imark:       ),nodes(:,imark:       ))
         end if
      end subroutine quick_sort_by_thickness
      subroutine quick_sort_by_thickness_partition(thick,nodes,marker)
         implicit none
         real(WP), dimension(  :) :: thick
         integer , dimension(:,:) :: nodes
         integer , intent(out)    :: marker
         integer :: i,j
         integer, dimension(3) :: i3tmp
         real(WP) :: dtmp,x
         x=thick(1)
         i=0; j=size(thick)+1
         do
            j=j-1
            do
               if (thick(j).le.x) exit
               j=j-1
            end do
            i=i+1
            do
               if (thick(i).ge.x) exit
               i=i+1
            end do
            if (i.lt.j) then
               dtmp =thick(  i); thick(  i)=thick(  j); thick(  j)=dtmp
               i3tmp=nodes(:,i); nodes(:,i)=nodes(:,j); nodes(:,j)=i3tmp
            else if (i.eq.j) then
               marker=i+1
               return
            else
               marker=i
               return
            endif
         end do
      end subroutine quick_sort_by_thickness_partition
   end subroutine sort_by_thickness


   !> Deallocate local structures
   subroutine kill_struct(this)
      implicit none
      class(ccl), intent(inout) :: this
      type(struct_type), pointer :: current => null()
      type(struct_type), pointer :: next => null()
      
      current => this%first_struct
      do while (associated(current))
         next => current%next
         if (associated(current%node)) then
            deallocate(current%node)
            nullify(current%node)
         end if
         deallocate(current)
         nullify(current)
         current => next
      end do
      
      ! Finish clean-up
      nullify(this%first_struct)
      nullify(this%my_struct)
      
      deallocate(this%struct_list)
      deallocate(this%struct_map_)
      if (this%max_interface_planes.gt.0) then
         deallocate(this%film_list)
         deallocate(this%film_map_)
      end if
      deallocate(this%parent)
      deallocate(this%parent_all)
      deallocate(this%parent_own)
      deallocate(this%per_)
      deallocate(this%per)
      if (this%max_interface_planes.eq.2) then
         deallocate(this%film_parent)
         deallocate(this%film_parent_all)
         deallocate(this%film_parent_own)
      end if
      
      return
   end subroutine kill_struct
   
   
end module ccl_class
