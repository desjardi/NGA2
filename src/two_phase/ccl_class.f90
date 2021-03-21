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
   public :: ccl ! make private/public more granular
   
   ! Default parameters for film CCL
   real(WP), parameter :: dot_threshold=-0.5_WP                      !< Maximum dot product of two interface normals for their respective cells to be considered film cells

   !> Structure object - could move to new module
   type :: struct
      ! 
      integer :: id                                       !< ID of struct
      integer :: parent                                   !< ID of parent struct
      integer :: rank = 0                                 !< Upper bound on height of tree; only used for roots
      integer :: counter = 0                              !< Counter for filling node array
      integer :: nnode != 1                               !< Number of cells contained in struct
      integer, dimension(:,:), allocatable :: node        !< List of cells contained in struct, dimension(1:nnode,3)
      integer, dimension(3) :: per = 0                    !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct

   !> Film object - could move to new module
   type, extends(struct) :: film
      integer :: phase                                    !< Film phase; 1 if liquid, 2 if gas
      integer, dimension(2) :: adjacent_structs           !< IDs of structures adjacent to gas film
   end type film

   !> CCL object definition
   type :: ccl
      
      ! This is our config
      class(config), pointer :: cfg                          !< This is the config object for which the CCL is built
      
      ! This is the name of the CCL
      character(len=str_medium) :: name='UNNAMED_CCL'        !< Solver name (default=UNNAMED_CCL)

      ! Volume fraction information
      real(WP), dimension(:,:,:), pointer :: VF

      ! Interface polygon information
      type(Poly_type), dimension(:,:,:,:), pointer :: poly   !< Array of IRL interface polygons (n,i,j,k)

      ! Maximum number of PLIC interfaces per cell
      integer :: max_interface_planes                        !< Number of planar interfaces per cell (0=VF-only, 1=PLIC, 2=R2P, etc)
      
      ! Structure and films
      ! type(struct), dimension(:), allocatable :: struct_list
      ! type(film)  , dimension(:), allocatable :: film_list
      type(struct), dimension(:), pointer :: struct_list => null()
      type(film)  , dimension(:), pointer :: film_list => null()

      ! CCL selection parameters - *read from input by user*
      real(WP) :: VFlo=1.0e-10_WP                            !< Minimum VF value considered for a structure to exist
      
      ! Feature counts
      integer :: n_struct, n_struct_max !, n_border_struct, n_border_struct_max
      integer :: n_film, n_film_max !, n_border_film, n_border_film_max

      ! Global tag offset
      integer :: id_offset
      integer :: synch_offset, film_synch_offset

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
      integer, dimension(:,:,:), allocatable :: id           !< ID of the structure that contains the cell
      integer, dimension(:,:,:), allocatable :: film_id      !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: film_phase   !< Phase of the film cell - 0/1/2/3 for none/liquid/gas/both

      ! Work arrays
      integer, dimension(:,:,:,:), allocatable :: idp          !< ID of the structure that contains the cell
      integer, dimension(:,:,:,:), allocatable :: film_idp     !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: film_pair      !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: border_id      !< ID of the film that contains the cell
      integer, dimension(:,:,:), allocatable :: film_border_id !< ID of the film that contains the cell

   contains
      procedure :: build_lists
      ! procedure, private :: is_connected
      ! procedure, private :: is_film_cell_upper
      ! procedure, private :: find
      ! procedure, private :: safe_find
      ! procedure, private :: union
      procedure, private :: label
      procedure, private :: struct_synch
      ! procedure, private :: struct_synch_per
      ! procedure, private :: film_synch
      ! procedure, private :: film_synch_per
      procedure, private :: kill_struct
   end type ccl
   
   
   !> Declare CCL algorithm constructor
   interface ccl
      procedure constructor
      ! procedure construct_struct_only
   end interface ccl
   
contains
   
   
   !> Default constructor for CCL algorithm
   function constructor(cfg,name) result(self)
      use messager, only: die
      implicit none
      type(ccl) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name

      ! Set the name for the object
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to cfg object
      self%cfg=>cfg
      
      ! Allocate id arrays
      allocate(self%id     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%id=0
      allocate(self%film_id(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_id=0

      ! Allocate film phase array
      allocate(self%film_phase(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_phase=0

      ! Allocate border id arrays
      allocate(self%border_id     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%border_id=0
      allocate(self%film_border_id(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_border_id=0

      ! Allocate film work arrays
      allocate(self%film_pair (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_pair=0

      ! Allocate periodicity work arrays
      ! allocate(self%idp     (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%idp=0
      ! allocate(self%film_idp(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%film_idp=0
      allocate(self%idp     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3)); self%idp=0
      allocate(self%film_idp(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3)); self%film_idp=0

   end function constructor
   
   !> Build lists of structures and films
   subroutine build_lists(this,VF,poly)
      implicit none
      class(ccl), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: VF
      type(Poly_type), dimension(:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in), optional:: poly

      ! Point to volume fraction field
      this%VF(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)=>VF

      ! Point to polygon object
      if (present(poly)) then
         this%poly(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:)=>poly
         print *,"max planes", this%max_interface_planes
         ! max_interface_planes = 2 - should user set this value?
      else
         this%max_interface_planes = 0
      end if

      ! Label features locally
      call this%label()

      ! Synchronize labels across procs
      call this%struct_synch()

      ! if (this%max_interface_planes.eq.2) call this%film_synch()

      ! Deallocate arrays
      call this%kill_struct()

   end subroutine build_lists

   !> Build local lists of structures and films
   subroutine label(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
      implicit none
      class(ccl), intent(inout) :: this
      integer :: npp,idd,i,j,k,ierr,dim,ii,jj,kk,cutoff_ind
      integer :: max_structs, max_films
      ! Using PLIC normal information
      logical :: is_contiguous = .true., has_normal, use_normal
      integer, dimension(3) :: pos
      real(WP), dimension(3) :: n1, n2, c1, c2
      logical :: is_film = .false., is_two_plane_film = .false.
      ! Only if two-plane cells are used
      real(WP), dimension(3) :: c22, n22

      ! Initialize id (for tagging) (is this redundant???)
      this%id = 0
      this%film_id = 0
      this%film_phase = 0

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
      this%id_offset = npp*(this%cfg%rank)         ! global tag offset; unique for each proc
      ! Allocate union-find data structure
      allocate(this%struct_list(this%id_offset+1:this%id_offset+max_structs))
      allocate(this%film_list(this%id_offset+1:this%id_offset+max_structs))
      ! initialize this%struct_list...
      !   this%struct_list(:)%volume = 0.0_WP
      !   this%struct_list(:)%per = 0
      this%struct_list%parent = 0
      this%struct_list%rank = 0
      this%struct_list%counter = 0
      this%struct_list%nnode = 0
      this%film_list%parent = 0
      this%film_list%rank = 0
      this%film_list%counter = 0
      this%film_list%nnode = 0
      ! print *,"0th check random id",this%id(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_),this%id(this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_)
      ! print *,"VFlo",this%VFlo,"mid VF",this%VF(this%cfg%imax_-20,this%cfg%jmax_-20,this%cfg%kmax_-20),"max_ VF",this%VF(this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_)
      ! initialize the global tag to this%id_offset - our "0"
      idd = this%id_offset
      if (this%max_interface_planes.gt.0) then
         ! outer loop over domain to find "corner" of structure
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
         
                  ! Find untagged point in liquid phase
                  if (this%VF(i,j,k).ge.this%VFlo) then

                     ! is_contiguous = .true.
                     has_normal = getNumberOfVertices(this%poly(1,i,j,k)).gt.0
                     if (has_normal) n1 = calculateNormal(this%poly(1,i,j,k))
                     is_two_plane_film = .false.
                     ! Rudimentary two-plane cell treatment
                     if (this%max_interface_planes.eq.2) then
                        if (getNumberOfVertices(this%poly(1,i,j,k)).ne.0.and.getNumberOfVertices(this%poly(2,i,j,k)).ne.0) then
                           n2 = calculateNormal(this%poly(2,i,j,k))
                           c1 = calculateCentroid(this%poly(1,i,j,k))
                           c2 = calculateCentroid(this%poly(2,i,j,k))
                           is_two_plane_film = (dot_product(n1,n2).lt.dot_threshold) ! .true. if normals point in opposing directions
                           ! if (.not.is_two_plane_film) print *,"ijk",i,j,k," not two plane dot product", dot_product(n1,n2)
                           ! if (is_two_plane_film) print *,"ijk",i,j,k," is two plane dot product", dot_product(n1,n2)
                           if (is_two_plane_film) then
                              if (dot_product(c2-c1,n2).gt.0.0_WP) then
                                 this%film_phase(i,j,k) = 1 ! Liquid film
                              else
                                 this%film_phase(i,j,k) = 2 ! Gas film
                                 ! this%id(i,j,k) = 0
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
                                    end if ! getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0
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
                                 this%idp(i,j,k,:) = this%struct_list(this%id(i,j,k))%per
                              else ! if gas film 
                                 do dim = 1,3 ! need parallel treatment
                                    pos = 0
                                    pos(dim) = -1
                                    ii = i + pos(1)
                                    jj = j + pos(2)
                                    kk = k + pos(3)
                                    ! if (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0) then
                                    if (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0.or.getNumberOfVertices(this%poly(2,ii,jj,kk)).gt.0) then
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
                                 if (getNumberOfVertices(this%poly(1,ii,jj,kk)).gt.0.or.getNumberOfVertices(this%poly(2,ii,jj,kk)).gt.0) then
                                    this%film_phase(ii,jj,kk) = this%film_phase(i,j,k)
                                 end if
                              end do
                           end if
                        end if
                     end if ! (this%max_interface_planes.eq.2)
                     if (.not.is_two_plane_film) then
                        if (has_normal) c1 = calculateCentroid(this%poly(1,i,j,k))
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
                              ! if (this%max_interface_planes.eq.2) then
                              if (getNumberOfVertices(this%poly(1,ii,jj,kk)).ne.0.and.getNumberOfVertices(this%poly(2,ii,jj,kk)).ne.0) then
                                 c22 = calculateCentroid(this%poly(2,ii,jj,kk))
                                 n22 = calculateNormal(this%poly(2,ii,jj,kk))
                                 is_contiguous = dot_product((c22-c2),n22).gt.0.0_WP ! .true. if liquid film
                                 is_film = .true.
                                 ! end if 
                              else
                                 ! If neighbor is one-plane cell
                                 is_contiguous = (dot_product(c2-c1,n2).ge.0.0_WP).or.(dot_product(c1-c2,n1).ge.0.0_WP)
                                 ! if (.not.is_contiguous) print *,"ijkVF",i,j,k,this%VF(i,j,k),"iijjkkVF",ii,jj,kk,this%VF(ii,jj,kk),"not contig: dot 1",dot_product(c2-c1,n2),"dot 2",dot_product(c1-c2,n1)
                                 is_film = (dot_product(n1,n2).lt.dot_threshold) ! .true. if normals point in opposing directions
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
                              end if ! is_film

                           end if ! use_normal
                           if (this%id(ii,jj,kk).gt.0 .and. is_contiguous) then ! Liquid neighbors should already be labeled
                              if (this%id(i,j,k).ne.0) then
                                 this%id(i,j,k) = union(this%struct_list,this%id(i,j,k),this%id(ii,jj,kk))
                                 ! if (is_film) print *,"irank",irank,"ijk",i,j,k,"iijjkk",ii,jj,kk,"film union","is_contiguous",is_contiguous,"dim dot",n1(dim)*n2(dim),"dot",dot_product(n1,n2)
                              else
                                 ! Propagate neighbor label
                                 this%id(i,j,k) = this%id(ii,jj,kk)
                                 ! if (is_film) print *,"irank",irank,"ijk",i,j,k,"iijjkk",ii,jj,kk,"film propagate"
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
                        this%idp(i,j,k,:) = this%struct_list(this%id(i,j,k))%per
                     end if ! .not.is_two_plane_film
                  end if ! this%VF >= this%VFlo
               end do ! k
            end do ! j
         end do ! i
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
                              ! if (is_film) print *,"irank",irank,"ijk",i,j,k,"iijjkk",ii,jj,kk,"film union","is_contiguous",is_contiguous,"dim dot",n1(dim)*n2(dim),"dot",dot_product(n1,n2)
                           else
                              ! Propagate neighbor label
                              this%id(i,j,k) = this%id(ii,jj,kk)
                              ! if (is_film) print *,"irank",irank,"ijk",i,j,k,"iijjkk",ii,jj,kk,"film propagate"
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
                     this%idp(i,j,k,:) = this%struct_list(this%id(i,j,k))%per
                  end if ! this%VF >= this%VFlo
               end do ! k
            end do ! j
         end do ! i
      end if
      ! print *,"first check random id",this%id(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_),this%id(this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_)

      ! Check upper boundaries for film cells
      ! this%cfg%imax_
      do j= this%cfg%jmin_,this%cfg%jmax_
         do k= this%cfg%kmin_,this%cfg%kmax_
            if (this%VF(this%cfg%imax_,j,k).ge.this%VFlo .and. this%film_phase(this%cfg%imax_,j,k).eq.0) then
               this%film_phase(this%cfg%imax_,j,k) = is_film_cell_upper(this%cfg%imax_,j,k,1)
            end if
         end do
      end do
      ! this%cfg%jmax_
      do i= this%cfg%imin_,this%cfg%imax_
         do k= this%cfg%kmin_,this%cfg%kmax_
            if (this%VF(i,this%cfg%jmax_,k).ge.this%VFlo .and. this%film_phase(i,this%cfg%jmax_,k).eq.0) then
               this%film_phase(i,this%cfg%jmax_,k) = is_film_cell_upper(i,this%cfg%jmax_,k,2)
            end if
         end do
      end do
      ! this%cfg%kmax_
      do i= this%cfg%imin_,this%cfg%imax_
         do j= this%cfg%jmin_,this%cfg%jmax_
            if (this%VF(i,j,this%cfg%kmax_).ge.this%VFlo .and. this%film_phase(i,j,this%cfg%kmax_).eq.0) then
               this%film_phase(i,j,this%cfg%kmax_) = is_film_cell_upper(i,j,this%cfg%kmax_,3)
            end if
         end do
      end do

      ! Number of this%struct_list graph nodes
      max_structs = idd-this%id_offset

      ! Reset idd
      idd = this%id_offset
      ! Tag film structs
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Find untagged point in film
               if (this%film_phase(i,j,k).gt.0) then
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
                  this%film_idp(i,j,k,:) = this%film_list(this%film_id(i,j,k))%per
               end if
            end do ! k
         end do ! j
      end do ! i

      ! Number of this%film_list graph nodes
      max_films = idd-this%id_offset

      ! Collapse tree
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_

               if (this%id(i,j,k).gt.0) then
                  this%id(i,j,k) = find(this%struct_list,this%id(i,j,k))
                  this%struct_list(this%id(i,j,k))%nnode = this%struct_list(this%id(i,j,k))%nnode + 1
                  this%idp(i,j,k,1) = max(this%idp(i,j,k,1),this%struct_list(this%id(i,j,k))%per(1)) 
                  this%idp(i,j,k,2) = max(this%idp(i,j,k,2),this%struct_list(this%id(i,j,k))%per(2)) 
                  this%idp(i,j,k,3) = max(this%idp(i,j,k,3),this%struct_list(this%id(i,j,k))%per(3)) 
                  this%struct_list(this%id(i,j,k))%per = this%idp(i,j,k,:)
               end if
               if (this%film_id(i,j,k).gt.0) then
                  this%film_id(i,j,k) = find(this%film_list,this%film_id(i,j,k))
                  this%film_list(this%film_id(i,j,k))%nnode = this%film_list(this%film_id(i,j,k))%nnode + 1
                  this%film_idp(i,j,k,1) = max(this%film_idp(i,j,k,1),this%film_list(this%film_id(i,j,k))%per(1)) 
                  this%film_idp(i,j,k,2) = max(this%film_idp(i,j,k,2),this%film_list(this%film_id(i,j,k))%per(2)) 
                  this%film_idp(i,j,k,3) = max(this%film_idp(i,j,k,3),this%film_list(this%film_id(i,j,k))%per(3)) 
                  this%film_list(this%film_id(i,j,k))%per = this%film_idp(i,j,k,:)
               end if
               
            end do ! k
         end do ! j
      end do ! i
      ! At this point, this%struct_list(this%id(i,j,k))%nnode > 0 only if this%id(i,j,k) is the root of a structure tree
      ! Allocate this%struct_list%node array, calculate this%n_struct
      this%n_struct = 0 ! counter for local number of structs
      do i=this%id_offset+1,this%id_offset+max_structs
         if (this%struct_list(i)%nnode.gt.0) then
            allocate(this%struct_list(i)%node(this%struct_list(i)%nnode,3))
            this%n_struct = this%n_struct + 1
            ! if (sum(this%struct_list(i)%per).gt.0) n_border_struct = n_border_struct + 1
         end if
      end do
      ! Allocate this%film_list%node array, calculate this%n_film
      this%n_film = 0
      do i=this%id_offset+1,this%id_offset+max_films
         if (this%film_list(i)%nnode.gt.0) then
            allocate(this%film_list(i)%node(this%film_list(i)%nnode,3))
            this%n_film = this%n_film + 1
            ! if (sum(this%film_list(i)%per).gt.0) n_border_struct = n_border_struct + 1
         end if
      end do
      ! total structures or parts of structures across all procs
      call MPI_ALLREDUCE(this%n_struct,this%n_struct_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      ! call MPI_ALLREDUCE(n_border_struct,n_border_struct_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(this%n_film,this%n_film_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      ! call MPI_ALLREDUCE(this%n_film_border_struct,this%n_film_border_struct_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)

      !! Whenever we switch to variable struct number per proc
      ! call MPI_ALLGATHER(this%n_struct, 1, MPI_INTEGER, this%n_struct_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      ! call MPI_ALLGATHER(n_border_struct, 1, MPI_INTEGER, n_border_struct_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      ! call MPI_ALLGATHER(this%n_film, 1, MPI_INTEGER, this%n_film_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      ! call MPI_ALLGATHER(this%n_film_border_struct, 1, MPI_INTEGER, this%n_film_border_struct_per_processor, 1, MPI_INTEGER, this%cfg%comm, ierr)
      !   this%synch_offset = sum(this%n_struct_per_processor(0:this%cfg%rank)) ???? not sum(this%n_struct_per_processor(0:this%cfg%rank)) - this%n_struct_per_processor(this%cfg%rank)??
      !   this%synch_offset = sum(n_border_struct_per_processor(0:this%cfg%rank))
      !   this%film_synch_offset = sum(this%n_film_per_processor(0:this%cfg%rank))
      !   this%film_synch_offset = sum(this%n_film_border_struct_per_processor(0:this%cfg%rank))
      !   allocate(this%struct_map_(this%synch_offset+1:this%synch_offset+this%n_struct)); this%struct_map_ = 0
      !   allocate(this%film_map_(this%film_synch_offset+1:this%film_synch_offset+this%n_film)); this%film_map_ = 0

      this%synch_offset = this%cfg%rank*this%n_struct_max
      this%film_synch_offset = this%cfg%rank*this%n_film_max
      ! Could allocate this%struct_map_(max_structs) and integrate this loop with above
      ! allocate(parent_    (this%synch_offset+1:this%synch_offset+this%n_struct_max)); parent_ = 0
      allocate(this%struct_map_(this%synch_offset+1:this%synch_offset+this%n_struct_max)); this%struct_map_ = 0
      allocate(this%film_map_(this%film_synch_offset+1:this%film_synch_offset+this%n_film_max)); this%film_map_ = 0
      idd = this%synch_offset
      do i=this%id_offset+1,this%id_offset+max_structs
      ! Only if this%struct_list has nodes
         if (this%struct_list(i)%nnode.eq.0) cycle
         idd = idd + 1
         this%struct_map_(idd) = i
         this%struct_list(i)%id = idd
      end do
      idd = this%film_synch_offset
      do i=this%id_offset+1,this%id_offset+max_films
      ! Only if this%film_list has nodes
         if (this%film_list(i)%nnode.eq.0) cycle
         idd = idd + 1
         this%film_map_(idd) = i
         this%film_list(i)%id = idd
      end do

      ! print *,"irank",irank,"start this%struct_map_ allgather"
      call fill_border_id(this%struct_list,this%id,this%border_id)
      call fill_border_id(this%film_list,this%film_id,this%film_border_id)

      ! Fill this%struct_list%node array
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
            if (this%id(i,j,k).gt.0) then
               this%struct_list(this%id(i,j,k))%counter = this%struct_list(this%id(i,j,k))%counter + 1
               this%struct_list(this%id(i,j,k))%node(this%struct_list(this%id(i,j,k))%counter,1) = i
               this%struct_list(this%id(i,j,k))%node(this%struct_list(this%id(i,j,k))%counter,2) = j
               this%struct_list(this%id(i,j,k))%node(this%struct_list(this%id(i,j,k))%counter,3) = k
            end if
            if (this%film_phase(i,j,k).gt.0) then
               this%film_list(this%film_id(i,j,k))%counter = this%film_list(this%film_id(i,j,k))%counter + 1
               this%film_list(this%film_id(i,j,k))%node(this%film_list(this%film_id(i,j,k))%counter,1) = i
               this%film_list(this%film_id(i,j,k))%node(this%film_list(this%film_id(i,j,k))%counter,2) = j
               this%film_list(this%film_id(i,j,k))%node(this%film_list(this%film_id(i,j,k))%counter,3) = k
            end if

            end do ! i
         end do ! j
      end do ! k
      print *,"this%n_struct",this%n_struct,"size"," max_structs",max_structs,"random id",this%id(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_),this%id(this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_)
      ! print *,"irank",irank,"start copy this%struct_list"
      ! ! Copy this%struct_list to my_struct
      ! call copy_to_my_struct

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

      ! ----------------------------------------------------------- !
      ! Is the cell above (i,j,k) in the dim direction a film cell? !
      ! ----------------------------------------------------------- !
      function is_film_cell_upper(i,j,k,dim) result(film_phase)
         implicit none
         integer, intent(in) :: i,j,k,dim
         integer :: ii,jj,kk
         logical :: is_film, use_normal, is_contiguous
         integer, dimension(3) :: pos
         real, dimension(3) :: nref, nloc, pref, ploc
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
               is_film = (dot_product(nref,nloc).lt.dot_threshold)
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

      ! Fill border cell id array with compact IDs
      subroutine fill_border_id(a_list_type,a_id_array,a_border_array)
         implicit none
         class(struct), dimension(this%id_offset+1:), intent(in) :: a_list_type
         integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(in)   :: a_id_array
         integer, dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(out)  :: a_border_array
         ! this%cfg%imin_
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               if (a_id_array(this%cfg%imin_,j,k).gt.0) then ! why is this necessary
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

   end subroutine label

   ! ----------------------------------------- !
   ! Synchronize struct tags across all procs  !
   ! ----------------------------------------- !
   subroutine struct_synch(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
      implicit none
      class(ccl), intent(inout) :: this
      integer :: i,j,k,n,ii,jj,kk,stop_,stop_global,counter,ierr,print_rank=2
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
      ! if (irank.eq.print_rank) print *,"Starting imin_"
      ! imin_
      if (this%cfg%imin_.ne.this%cfg%imin) then
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%border_id(this%cfg%imin_,j,k).gt.0.and.this%border_id(this%cfg%imin_-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%border_id(this%cfg%imin_,j,k)).eq.this%border_id(this%cfg%imin_,j,k)) then
                     ! if (irank.eq.print_rank) print *,"create parentage"
                     ! call union(this%border_id(this%cfg%imin_,j,k),this%border_id(this%cfg%imin_-1,j,k))
                     ! is_contiguous = is_connected(this%cfg%imin_,j,k,1)
                     ! if (is_contiguous) parent(this%border_id(this%cfg%imin_,j,k)) = this%border_id(this%cfg%imin_-1,j,k)
                     if (is_connected(this%cfg%imin_,j,k,1)) this%parent(this%border_id(this%cfg%imin_,j,k)) = this%border_id(this%cfg%imin_-1,j,k)
                  elseif (find(this%border_id(this%cfg%imin_,j,k)).ne.find(this%border_id(this%cfg%imin_-1,j,k))) then
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()"
                     ! if (irank.eq.print_rank) print *,"parent of ghost/border",this%parent(this%border_id(this%cfg%imin_-1,j,k)),this%parent(this%border_id(this%cfg%imin_,j,k)) 
                     ! if (irank.eq.print_rank) print *,"root   of ghost/border",find_only(this%border_id(this%cfg%imin_-1,j,k)),find_only(this%border_id(this%cfg%imin_,j,k)) 
                     ! call union(this%border_id(this%cfg%imin_-1,j,k),this%parent(this%border_id(this%cfg%imin_,j,k)))
                     if (is_connected(this%cfg%imin_,j,k,1)) then
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
      ! if (irank.eq.2) print *,"this%synch_offset",this%synch_offset
      ! if (irank.eq.2) print *,"irank",irank,"border id ",this%border_id(this%cfg%imin_,:,:)
      ! if (irank.eq.2) print *,"irank",irank,"ghost id  ",this%border_id(this%cfg%imin_-1,:,:)
      ! if (irank.eq.2) print *,"irank",irank,"this%parent", this%parent
   
      ! if (irank.eq.print_rank) print *,"Starting this%cfg%jmin_"
      ! this%cfg%jmin_
      if (this%cfg%jmin_.ne.this%cfg%jmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%border_id(i,this%cfg%jmin_,k).gt.0.and.this%border_id(i,this%cfg%jmin_-1,k).gt.0) then
                  ! if (irank.eq.print_rank) print *,"ghost border",this%border_id(i,this%cfg%jmin_-1,k),this%border_id(i,this%cfg%jmin_,k)
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%border_id(i,this%cfg%jmin_,k)).eq.this%border_id(i,this%cfg%jmin_,k)) then 
                     ! if (irank.eq.print_rank) print *,"create parentage"
                     ! call union(this%border_id(i,this%cfg%jmin_,k),this%border_id(i,this%cfg%jmin_-1,k))
                     if (is_connected(i,this%cfg%jmin_,k,2)) this%parent(this%border_id(i,this%cfg%jmin_,k)) = this%border_id(i,this%cfg%jmin_-1,k)
                  elseif (find(this%border_id(i,this%cfg%jmin_,k)).ne.find(this%border_id(i,this%cfg%jmin_-1,k))) then
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()" 
                  ! if (irank.eq.print_rank) print *,"structure borders two structures - union()" 
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()" 
                     ! call union(this%border_id(i,this%cfg%jmin_-1,k),this%parent(this%border_id(i,this%cfg%jmin_,k)))
                     if (is_connected(i,this%cfg%jmin_,k,2)) then
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
      ! if (irank.eq.2) print *,"irank",irank,"this%parent", this%parent
   
      ! if (irank.eq.print_rank) print *,"Starting this%cfg%kmin_"
      ! this%cfg%kmin_
      if (this%cfg%kmin_.ne.this%cfg%kmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               ! If the neighboring ghost cell is filled
               if (this%border_id(i,j,this%cfg%kmin_).gt.0.and.this%border_id(i,j,this%cfg%kmin_-1).gt.0) then
                  ! if (irank.eq.print_rank) print *,"ghost border",this%border_id(i,j,this%cfg%kmin_-1),this%border_id(i,j,this%cfg%kmin_)
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%border_id(i,j,this%cfg%kmin_)).eq.this%border_id(i,j,this%cfg%kmin_)) then 
                     ! if (irank.eq.print_rank) print *,"create parentage"
                     ! call union(this%border_id(i,j,this%cfg%kmin_),this%border_id(i,j,this%cfg%kmin_-1))
                     if (is_connected(i,j,this%cfg%kmin_,3)) this%parent(this%border_id(i,j,this%cfg%kmin_)) = this%border_id(i,j,this%cfg%kmin_-1)
                  elseif (find(this%border_id(i,j,this%cfg%kmin_)).ne.find(this%border_id(i,j,this%cfg%kmin_-1))) then
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()" 
                     ! if (irank.eq.print_rank) print *,"parent of ghost/border",this%parent(this%border_id(i,j,this%cfg%kmin_-1)),this%parent(this%border_id(i,j,this%cfg%kmin_)) 
                     ! if (irank.eq.print_rank) print *,"root   of ghost/border",find_only(this%border_id(i,j,this%cfg%kmin_-1)),find_only(this%border_id(i,j,this%cfg%kmin_)) 
                     ! call union(this%border_id(i,j,this%cfg%kmin_-1),this%parent(this%border_id(i,j,this%cfg%kmin_)))
                     if (is_connected(i,j,this%cfg%kmin_,3)) then 
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
      ! if (irank.eq.2) print *,"irank",irank,"this%parent", this%parent
   
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
   
         ! print *,"counter",counter,"irank",irank,"after mpi"
         ! Set self-parents back to selves
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent_all(i).eq.huge(1)) this%parent_all(i) = i
            ! if (this%parent_own(i).eq.huge(1)) this%parent_own(i) = i
           end do
         ! if (irank.eq.print_rank) print *,"counter",counter,"irank",irank,"this%parent_all before flatten",this%parent_all
         ! print *,"counter",counter,"irank",irank,"this%parent_all before flatten",this%parent_all
         ! if (irank.eq.print_rank) print *,"counter",counter,"irank",irank,"this%parent_own before flatten",this%parent_own
         ! print *,"counter",counter,"irank",irank,"this%parent_own before flatten",this%parent_own
         ! Flatten trees - is this necessary?
         do i= 1,this%cfg%nproc*this%n_struct_max
            ! print *,"counter",counter,"irank",irank,"flatten all tree i",i
            this%parent_all(i) = find_all(i)
            ! print *,"counter",counter,"irank",irank,"flatten own tree i",i
            this%parent_own(i) = find_own(i)
                  ! if (irank.eq.print_rank) print *,"counter",counter,"irank",irank,"i",i,"this%parent_all i",this%parent_all
         end do
         ! print *,"counter",counter,"irank",irank,"after flatten"
         ! ! Fill this%parent with selves again for printing
         ! do i= 1,this%cfg%nproc*this%n_struct_max
         !    this%parent(i) = i
         ! end do
         ! if (irank.eq.print_rank) print *, "irank",irank,"    struct_id         ",this%parent
         ! print *,"counter",counter,"irank",irank,"this%parent_all before find",this%parent_all
         ! print *,"counter",counter,"irank",irank,"before reconciliation"
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
                  ! print *,"counter",counter,"irank",irank,"reconciliation union between",this%parent_own(i),this%parent(i),"find_own(this%parent_own) find(this%parent) find(this%parent_own)",find_own(this%parent_own(i)),find(this%parent(i)),find(this%parent_own(i))
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
   
      ! if (irank.eq.1) print *,'counter',counter
      ! Update this%struct_list%parent and point all parents to root
      do i=this%synch_offset+1,this%synch_offset+this%n_struct
         this%struct_list(this%struct_map_(i))%parent = find(this%parent(i))
      end do
         ! if (irank.eq.print_rank) print *,"irank",irank,"this%parent     after find ",this%parent
         ! Update id array with compacted and syncrhonized ids
      do i=this%synch_offset+1,this%synch_offset+this%n_struct
         do n=1,this%struct_list(this%struct_map_(i))%nnode
            ii = this%struct_list(this%struct_map_(i))%node(n,1)
            jj = this%struct_list(this%struct_map_(i))%node(n,2)
            kk = this%struct_list(this%struct_map_(i))%node(n,3)
            this%id(ii,jj,kk) = this%struct_list(this%struct_map_(i))%parent
            ! this%id(ii,jj,kk) = this%struct_list(this%id(ii,jj,kk))%parent
         end do
      end do
      ! print *,'rank',irank,'Done with tree collapse'
   
     call struct_synch_per
   !   print *,"irank",irank,'id before boundary update',this%id(11,11,:)
   !   print *,"irank",irank,'this%parent before boundary update',this%parent
   !   print *,"irank",irank,'this%border_id before boundary update [border,ghost]',this%border_id(4,11,11),this%border_id(3,11,11)
   !   call MPI_BARRIER(this%cfg%comm,ierr)
      !! Update domain boundaries
      ! this%cfg%imin
      if (this%cfg%imin_.eq.this%cfg%imin) then
         do j= this%cfg%jmin_,this%cfg%jmax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the border cell and the neighboring ghost cell are filled
               if (this%id(this%cfg%imin,j,k).gt.0.and.this%border_id(this%cfg%imin-1,j,k).gt.0) then
                  ! if (irank.eq.print_rank) print *,"ijk",this%cfg%imin,j,k
                  ! if (irank.eq.print_rank) print *,"i ghost border",this%border_id(this%cfg%imin-1,j,k),this%id(this%cfg%imin,j,k)
                  ! if (irank.eq.print_rank) print *,"i border parent",this%parent(this%id(this%cfg%imin,j,k))
                  ! if (irank.eq.print_rank) print *,"this%parent",this%parent
                  ! if (irank.eq.print_rank) print *,"find border root",find(this%id(this%cfg%imin,j,k))
                  ! if (irank.eq.print_rank) print *,"find ghost root",find(this%border_id(this%cfg%imin-1,j,k))
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%id(this%cfg%imin,j,k)).eq.this%id(this%cfg%imin,j,k)) then
                     ! if (irank.eq.print_rank) print *,"create parentage"
                     if (is_connected(this%cfg%imin,j,k,1)) call union(this%id(this%cfg%imin,j,k),this%border_id(this%cfg%imin-1,j,k))
                     ! this%parent(this%id(this%cfg%imin,j,k)) = this%border_id(this%cfg%imin-1,j,k)
                  elseif (find(this%id(this%cfg%imin,j,k)).ne.find(this%border_id(this%cfg%imin-1,j,k))) then
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()"
                     ! if (irank.eq.print_rank) print *,"parent of ghost/border",this%parent(this%border_id(this%cfg%imin-1,j,k)),this%parent(this%border_id(this%cfg%imin,j,k)) 
                     ! if (irank.eq.print_rank) print *,"root   of ghost/border",find_only(this%border_id(this%cfg%imin-1,j,k)),find_only(this%border_id(this%cfg%imin,j,k)) 
                     if (is_connected(this%cfg%imin,j,k,1)) then
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
      ! if (irank.eq.2) print *,"this%synch_offset",this%synch_offset
      ! if (irank.eq.2) print *,"irank",irank,"border id ",this%border_id(this%cfg%imin,:,:)
      ! if (irank.eq.2) print *,"irank",irank,"ghost id  ",this%border_id(this%cfg%imin-1,:,:)
      ! if (irank.eq.2) print *,"irank",irank,"this%parent", this%parent
   
      ! if (irank.eq.print_rank) print *,"Starting this%cfg%jmin"
      ! print *,"irank",irank,'starting this%cfg%jmin'
      ! this%cfg%jmin
      if (this%cfg%jmin_.eq.this%cfg%jmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do k= this%cfg%kmin_,this%cfg%kmax_
               ! If the neighboring ghost cell is filled
               if (this%id(i,this%cfg%jmin,k).gt.0.and.this%border_id(i,this%cfg%jmin-1,k).gt.0) then
                  ! if (irank.eq.print_rank) print *,"j ghost border",this%border_id(i,this%cfg%jmin-1,k),this%id(i,this%cfg%jmin,k)
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%id(i,this%cfg%jmin,k)).eq.this%id(i,this%cfg%jmin,k)) then 
                     ! if (irank.eq.print_rank) print *,"create parentage"
                     if (is_connected(i,this%cfg%jmin,k,2)) call union(this%id(i,this%cfg%jmin,k),this%border_id(i,this%cfg%jmin-1,k))
                     ! this%parent(this%id(i,this%cfg%jmin,k)) = this%border_id(i,this%cfg%jmin-1,k)
                  elseif (find(this%id(i,this%cfg%jmin,k)).ne.find(this%border_id(i,this%cfg%jmin-1,k))) then
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()" 
                     ! call union(this%border_id(i,this%cfg%jmin-1,k),this%parent(this%id(i,this%cfg%jmin,k)))
                     if (is_connected(i,this%cfg%jmin,k,2)) then
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
      ! if (irank.eq.2) print *,"irank",irank,"this%parent", this%parent
   
      ! if (irank.eq.print_rank) print *,"Starting this%cfg%kmin"
      ! print *,"irank",irank,'starting this%cfg%kmin'
      ! this%cfg%kmin
      if (this%cfg%kmin_.eq.this%cfg%kmin) then
         do i= this%cfg%imin_,this%cfg%imax_
            do j= this%cfg%jmin_,this%cfg%jmax_
               ! If the neighboring ghost cell is filled
               if (this%id(i,j,this%cfg%kmin).gt.0.and.this%border_id(i,j,this%cfg%kmin-1).gt.0) then
                  ! if (irank.eq.print_rank) print *,"k ghost border",this%border_id(i,j,this%cfg%kmin-1),this%id(i,j,this%cfg%kmin)
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (this%parent(this%id(i,j,this%cfg%kmin)).eq.this%id(i,j,this%cfg%kmin)) then 
                     ! if (irank.eq.print_rank) print *,"create parentage"
                     if (is_connected(i,j,this%cfg%kmin,3)) call union(this%id(i,j,this%cfg%kmin),this%border_id(i,j,this%cfg%kmin-1))
                     ! this%parent(this%id(i,j,this%cfg%kmin)) = this%border_id(i,j,this%cfg%kmin-1)
                  elseif (find(this%id(i,j,this%cfg%kmin)).ne.find(this%border_id(i,j,this%cfg%kmin-1))) then
                     ! if (irank.eq.print_rank) print *,"structure borders two structures - union()" 
                     ! if (irank.eq.print_rank) print *,"parent of ghost/border",this%parent(this%border_id(i,j,this%cfg%kmin-1)),this%parent(this%id(i,j,this%cfg%kmin)) 
                     ! if (irank.eq.print_rank) print *,"root   of ghost/border",find_only(this%border_id(i,j,this%cfg%kmin-1)),find_only(this%id(i,j,this%cfg%kmin)) 
                     ! call union(this%border_id(i,j,this%cfg%kmin-1),this%parent(this%id(i,j,this%cfg%kmin)))
                     if (is_connected(i,j,this%cfg%kmin,3)) then
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
      ! if (irank.eq.print_rank) print *,"irank",irank,"this%parent after boundary synch", this%parent
   
      ! initialize global stop criterion
      stop_global = 1
   
      ! initialize a counter
      counter = 0
      ! print *,"irank",irank,'starting do while 2'
   
      do while (stop_global.ne.0)
   
         ! Initialize local stop flag
         stop_ = 0
   
         this%parent_own = this%parent
         ! Set self-parents to huge(1)
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent(i).eq.i) this%parent(i) = huge(1)
         end do
   
         call MPI_ALLREDUCE(this%parent,this%parent_all,this%cfg%nproc*this%n_struct_max,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
   
         ! print *,"counter",counter,"irank",irank,"after mpi"
         ! Set self-parents back to selves
         do i= 1,this%cfg%nproc*this%n_struct_max
            if (this%parent_all(i).eq.huge(1)) this%parent_all(i) = i
            ! if (this%parent_own(i).eq.huge(1)) this%parent_own(i) = i
         end do
         ! if (irank.eq.print_rank) print *,"counter",counter,"irank",irank,"this%parent_all before flatten",this%parent_all
         ! print *,"counter",counter,"irank",irank,"this%parent_all before flatten",this%parent_all
         ! if (irank.eq.print_rank) print *,"counter",counter,"irank",irank,"this%parent_own before flatten",this%parent_own
         ! print *,"counter",counter,"irank",irank,"this%parent_own before flatten",this%parent_own
         ! Flatten trees - is this necessary?
         do i= 1,this%cfg%nproc*this%n_struct_max
            ! print *,"counter",counter,"irank",irank,"flatten all tree i",i
            this%parent_all(i) = find_all_2(i,i)
            ! print *,"counter",counter,"irank",irank,"flatten own tree i",i
            this%parent_own(i) = find_own(i)
            ! if (irank.eq.print_rank) print *,"counter",counter,"irank",irank,"i",i,"this%parent_all i",this%parent_all
         end do
         ! print *,"counter",counter,"irank",irank,"after flatten"
         ! ! Fill this%parent with selves again for printing
         ! do i= 1,this%cfg%nproc*this%n_struct_max
         !    this%parent(i) = i
         ! end do
         ! if (irank.eq.print_rank) print *, "irank",irank,"    struct_id         ",this%parent
         ! print *,"counter",counter,"irank",irank,"this%parent_all before find",this%parent_all
         ! print *,"counter",counter,"irank",irank,"before reconciliation"
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
   
      ! if (irank.eq.1) print *,'counter',counter
      ! Update this%struct_list%parent and point all parents to root
      do i=this%synch_offset+1,this%synch_offset+this%n_struct
            this%struct_list(this%struct_map_(i))%parent = find(this%parent(i))
            ! print *,'this%struct_list%parent',this%struct_list(this%struct_map_(i))%parent
      end do
      ! if (irank.eq.print_rank) print *,"irank",irank,"this%parent     after find ",this%parent
      ! print *,'id before update',this%id(11,11,:)
      ! Update id array with compacted and syncrhonized ids
      do i=this%synch_offset+1,this%synch_offset+this%n_struct
         do n=1,this%struct_list(this%struct_map_(i))%nnode
            ii = this%struct_list(this%struct_map_(i))%node(n,1)
            jj = this%struct_list(this%struct_map_(i))%node(n,2)
            kk = this%struct_list(this%struct_map_(i))%node(n,3)
            this%id(ii,jj,kk) = this%struct_list(this%struct_map_(i))%parent
         end do
      end do
      ! print *,'rank',irank,'Done with tree collapse'
   
      ! ! Update my_struct%id
      ! ! start marching thru list, starting at first_struct
      ! my_struct => first_struct
      
      ! do while(associated(my_struct))
      !    ii = my_struct%node(1,1)
      !    jj = my_struct%node(1,2)
      !    kk = my_struct%node(1,3)
      !    ! Update structure id with id cell value
      !    my_struct%id  = this%id(ii,jj,kk)
         
      !    ! Go to next structure
      !    my_struct => my_struct%next
      
      ! end do
   
      return
   
   contains
   
      ! For debugging only
      function find_only(a_x) result(a_y)
         implicit none
         integer :: a_x
         integer :: a_y
   
         a_y = a_x
         ! find root of x with path compression
         do while ( a_y.ne.this%parent(a_y) )
            a_y = this%parent(a_y)
         end do
         return
      end function find_only
   
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

      ! --------------------------------- !
      ! Use PLIC normals to determine if  !
      ! two cells are in the same struct  !
      ! --------------------------------- !
      function is_connected(i,j,k,dim)
         implicit none
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

      ! ------------------------------------------------ !
      ! Synchronize struct periodicity across all procs  !
      ! ------------------------------------------------ !
      subroutine struct_synch_per
         implicit none
         ! integer :: i,j,k,n,ii,jj,kk,stop_,stop_global,counter,ierr
         integer :: flag

         ! Allocate local and global perodicity arrays
         allocate(this%per_(this%synch_offset+1:this%synch_offset+this%n_struct_max,1:3)); this%per_ = 0
         allocate(this%per (this%cfg%nproc*this%n_struct_max,1:3)); this%per = 0 

         ! Fill per_ array
         do i=this%synch_offset+1,this%synch_offset+this%n_struct
            this%per_(i,:) = this%struct_list(this%struct_map_(i))%per
         end do

         ! Communitcate per
         call MPI_ALLGATHER(this%per_(:,1),this%n_struct_max,MPI_INTEGER,this%per(:,1),this%n_struct_max,MPI_INTEGER,this%cfg%comm,ierr)
         call MPI_ALLGATHER(this%per_(:,2),this%n_struct_max,MPI_INTEGER,this%per(:,2),this%n_struct_max,MPI_INTEGER,this%cfg%comm,ierr)
         call MPI_ALLGATHER(this%per_(:,3),this%n_struct_max,MPI_INTEGER,this%per(:,3),this%n_struct_max,MPI_INTEGER,this%cfg%comm,ierr)
         ! Update parent per
         do i=1,this%cfg%nproc*this%n_struct_max
            this%per(this%parent(i),:) = max(this%per(this%parent(i),:),this%per(i,:))
         end do

         do i=this%synch_offset+1,this%synch_offset+this%n_struct
            do n=1,this%struct_list(this%struct_map_(i))%nnode
               ii = this%struct_list(this%struct_map_(i))%node(n,1)
               jj = this%struct_list(this%struct_map_(i))%node(n,2)
               kk = this%struct_list(this%struct_map_(i))%node(n,3)
               this%idp(ii,jj,kk,:) = this%per(this%id(ii,jj,kk),:)
            end do
         end do

         return

      end subroutine struct_synch_per
   
   end subroutine struct_synch

   ! --------------------------- !
   ! Deallocate local structures !
   ! --------------------------- !
   subroutine kill_struct(this)
      implicit none
      class(ccl), intent(inout) :: this

      deallocate(this%struct_list)
      deallocate(this%film_list)
      deallocate(this%struct_map_)
      deallocate(this%film_map_)
      deallocate(this%parent)
      deallocate(this%parent_all)
      deallocate(this%parent_own)
      deallocate(this%per_)
      deallocate(this%per)
      ! deallocate(this%film_parent)
      ! deallocate(this%film_parent_all)
      ! deallocate(this%film_parent_own)


      return
   end subroutine kill_struct

end module ccl_class
