!> Connected component labeling class: identifies Lagrangian objects from a Eulerian logical field
!> and provides unstructured mapping to traverse these objects
module amrcclabel_class
   use precision,       only: WP
   use string,          only: str_medium
   use amrdata_class,   only: amrdata
   use amrsolver_class, only: amrsolver
   use amrgrid_class,   only: amrgrid
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: amrcclabel,make_label_ftype,same_label_ftype,stats_type
   
   
   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%
   
   !> Structure object
   type :: struct_type
      integer :: parent                                   !< ID of parent struct
      integer :: n_                                       !< Number of local cells contained in struct
      integer, dimension(3) :: per                        !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct_type
   
   !> Statistics object
   type :: stats_type
      integer :: id                    !< ID of structure
      real(WP) :: vol                  !< Volme of structure
      real(WP), dimension(3) :: com    !< Center of mass of structure
   end type stats_type
   
   !> amrcclabel object definition
   type :: amrcclabel
      character(len=str_medium) :: name = 'UNNAMED_CCLABEL'
      ! ID of the structure that contains each cell
      type(amrdata) :: id
      ! Array of structures
      integer :: nstruct
      type(struct_type), dimension(:), allocatable :: struct
      ! Ghost cells
      integer :: nover=1
      ! Associated amr grid 
      class(amrgrid), pointer, private :: amr => null()
      ! Temporary arrays for interlevel sync
      type(amrdata) :: tmp_id,tmp_conflict
   contains
      procedure :: initialize
      procedure :: build
      procedure :: empty
      procedure :: compute_stats
      procedure :: finalize
   end type amrcclabel
   
   !> Type of the make_label function used to generate a structure
   interface
      logical function make_label_ftype(pdata,lo,i,j,k)
         use precision,    only: WP
         real(WP), dimension(:,:,:,:), intent(in) :: pdata
         integer, dimension(3), intent(in) :: lo
         integer, intent(in) :: i,j,k
      end function make_label_ftype
   end interface
   
   !> Type of the same_label function used to connect two structures
   interface
      logical function same_label_ftype(pdata,lo,i,j,k,ii,jj,kk)
         use precision,    only: WP
         real(WP), dimension(:,:,:,:), intent(in) :: pdata
         integer, dimension(3), intent(in) :: lo
         integer, intent(in) :: i,j,k,ii,jj,kk
      end function same_label_ftype
   end interface
   
contains
   
   
   !> Initialization for amrcclabel class
   subroutine initialize(this,amr,name)
      use amrdata_class, only: interp_none
      implicit none
      class(amrcclabel) :: this
      class(amrgrid), target, intent(in) :: amr 
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Point cclabel to amr grid
      this%amr => amr
      ! Allocate and initialize ID array
      call this%id%initialize(amr,name='id',ncomp=1,ng=this%nover,interp=interp_none);! this%id%parent=>this
      call this%id%register() ! Update with regriding
      call this%id%setval(val=0.0_WP)
      ! Allocate temporary arrays for interlevel sync
      call this%tmp_id%initialize(amr,name='tmp_id',ncomp=1,ng=this%nover)
      call this%tmp_conflict%initialize(amr,name='tmp_conflict',ncomp=1,ng=this%nover)
      call this%tmp_id%register() ! Update with regriding
      call this%tmp_conflict%register() ! Update with regriding
      ! Zero structures
      this%nstruct=0
   end subroutine initialize
   
   
   !> Build structure using the user-set test functions
   subroutine build(this,make_label,same_label,coarse_make_label,coarse_same_label,data)
      use amrdata_class,    only: amrdata
      use amrdata_class, only: interp_none
      implicit none
      class(amrcclabel), intent(inout) :: this
      procedure(make_label_ftype) :: make_label,coarse_make_label
      procedure(same_label_ftype) :: same_label,coarse_same_label
      type(amrdata), intent(in) :: data
      type(amrdata) :: idp
      integer :: nstruct_,stmin,stmax
      integer, dimension(:), allocatable :: parent             !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_all         !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_own         !< Resolving structure id across procs

      ! Initialized id to zero on all levels
      call this%id%setval(0.0_WP)

      ! Build CCL on finest level
      call build_lvl(data%amr%maxlvl,make_label,same_label)


      ! Create unique IDs for each structure on coarser levels
      build_coarser: block
         integer :: lvl
         integer, dimension(3) :: ref_ratio
         do lvl = data%amr%maxlvl-1, 0, -1   ! finest-1 → coarsest

            ref_ratio(1)=data%amr%rrefx(lvl)
            ref_ratio(2)=data%amr%rrefy(lvl)
            ref_ratio(3)=data%amr%rrefz(lvl)

            ! Implemented in C to get access to additional functions
            call restrict_unique_id( &
               this%id%mf(lvl),     & ! coarse
               this%id%mf(lvl+1),   & ! fine
               ref_ratio,           &
               this%amr%geom(lvl+1) )

            ! Build CCL on coarse level
            call build_lvl(lvl,coarse_make_label,coarse_same_label)

         end do
      end block build_coarser

      ! testing_end_build: block 
      !    integer :: lvl
      !    do lvl = 0,data%amr%maxlvl
      !       call print_ids(lvl,"after build")
      !    end do
      ! end block testing_end_build

   contains

      !> Build structure on a level using user-set test functions
      subroutine build_lvl(lvl,make_label,same_label)
         integer, intent(in) :: lvl
         procedure(make_label_ftype) :: make_label
         procedure(same_label_ftype) :: same_label
         logical :: finest
         integer :: nstruct_work

         ! Set finest logical
         finest=.false.
         if (lvl.eq.this%amr%maxlvl) finest=.true.

         ! Start by cleaning up
         call this%empty()
         
         ! Then allocate struct to a default size
         nstruct_=0
         allocate(this%struct(min_struct_size))
         this%struct(:)%parent=0
         this%struct(:)%per(1)=0
         this%struct(:)%per(2)=0
         this%struct(:)%per(3)=0
         this%struct(:)%n_=0

         ! Add any ids from finer levels to struct array
         previous_ids: block 
            use mpi_f08, only: MPI_ALLREDUCE,MPI_INTEGER,MPI_MAX,MPI_IN_PLACE
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k,ierr
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
            ! Only do if on coarser level
            if (finest) exit previous_ids
            ! Set structure counter to not overwrite any existing structures
            nstruct_=this%nstruct
            
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data arrays
               pid=>this%id%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pid(i,j,k,1).gt.0.5_WP) then
                     call add_existing(nint(pid(i,j,k,1)))
                  end if
               end do; end do; end do
            end do
         end block previous_ids

         ! Perform a first pass to build proc-local structures and corresponding tree
         first_pass: block
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k
            integer :: ii,jj,kk,dim
            integer, dimension(3) :: pos
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid,pdata

            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data arrays
               pid=>this%id%mf(lvl)%dataptr(mfi)
               ! pidp=>idp%mf(lvl)%dataptr(mfi)
               pdata=>data%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Find next cell in a structure
                  if (make_label(pdata,lbound(pdata),i,j,k)) then
                     ! Loop through one-sided neighbors
                     do dim=1,3
                        pos=0; pos(dim)=-1
                        ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                        ! Check if neighbor is labeled
                        if (pid(ii,jj,kk,1).gt.0.5_WP) then
                           ! Neighbor is labeled, but are we?
                           if (pid(i,j,k,1).gt.0.5_WP) then
                              ! We already have a label, perform a union of both labels
                              if (same_label(pdata,lbound(pdata),i,j,k,ii,jj,kk)) then
                                 pid(i,j,k,1)=union_struct(nint(pid(i,j,k,1)),nint(pid(ii,jj,kk,1)))
                              end if
                           else
                              ! We don't have a label, check if we take the neighbor's label
                              if (same_label(pdata,lbound(pdata),i,j,k,ii,jj,kk)) then
                                 pid(i,j,k,1)=pid(ii,jj,kk,1)
                              else
                                 pid(i,j,k,1)=add()
                              end if
                           end if
                        end if
                     end do
                     ! If no neighbor was labeled, we need a new structure
                     if (pid(i,j,k,1).eq.0) then 
                        pid(i,j,k,1)=add()
                     end if
                     ! ! Identify periodicity cases
                     ! if (this%amr%xper.and.i.eq.this%pg%imax) this%struct(pid(i,j,k,1))%per(1)=1
                     ! if (this%amr%yper.and.j.eq.this%pg%jmax) this%struct(pid(i,j,k,1))%per(2)=1
                     ! if (this%amr%zper.and.k.eq.this%pg%kmax) this%struct(pid(i,j,k,1))%per(3)=1
                     ! pidp(i,j,k,:)=this%struct(pid(i,j,k,1))%per
                  end if
               end do; end do; end do
            end do
         end block first_pass

         ! Now collapse the tree, count the cells and resolve periodicity in each structure
         collapse_tree: block
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid,pidp
            ! Loop over tiles
            call data%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pid=>this%id%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pid(i,j,k,1).gt.0.5_WP) then
                     pid(i,j,k,1)=rootify_struct(nint(pid(i,j,k,1)))
                     this%struct(nint(pid(i,j,k,1)))%n_=this%struct(nint(pid(i,j,k,1)))%n_+1
                     ! pidp(i,j,k,1)=max(nint(pidp(i,j,k,1)),this%struct(nint(pid(i,j,k,1)))%per(1))
                     ! pidp(i,j,k,2)=max(nint(pidp(i,j,k,2)),this%struct(nint(pid(i,j,k,1)))%per(2))
                     ! pidp(i,j,k,3)=max(nint(pidp(i,j,k,3)),this%struct(nint(pid(i,j,k,1)))%per(3))
                     ! this%struct(nint(pid(i,j,k,1)))%per=nint(pidp(:,i,j,k))
                  end if
               end do; end do; end do
            end do
         end block collapse_tree
         
         ! Compact structure array
         compact_tree: block
            use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER,MPI_MAX
            integer :: i,j,k,n,ierr
            integer, dimension(:), allocatable :: my_nstruct,all_nstruct,idmap
            type(struct_type), dimension(:), allocatable :: tmp
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
            
            ! If not finest just compute the number of structures
            if (.not.finest) then
               nstruct_=0
               do n=1,size(this%struct,dim=1)
                  if (this%struct(n)%n_.gt.0) nstruct_=n
               end do
               call MPI_ALLREDUCE(nstruct_,nstruct_work,1,MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)
            else
               ! Count exact number of local structures
               nstruct_=0
               do n=1,size(this%struct,dim=1)
                  if (this%struct(n)%n_.gt.0) nstruct_=nstruct_+1
               end do
               ! Gather this info to ensure unique index
               allocate( my_nstruct(0:this%amr%nproc-1)); my_nstruct=0; my_nstruct(this%amr%rank)=nstruct_
               allocate(all_nstruct(0:this%amr%nproc-1)); call MPI_ALLREDUCE(my_nstruct,all_nstruct,this%amr%nproc,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
               stmin=1
               if (this%amr%rank.gt.0) stmin=stmin+sum(all_nstruct(0:this%amr%rank-1))
               nstruct_work=sum(all_nstruct)
               deallocate(my_nstruct,all_nstruct)
               stmax=stmin+nstruct_-1
               ! Generate an index map
               allocate(idmap(1:size(this%struct,dim=1))); idmap=0
               nstruct_=0
               do n=1,size(this%struct,dim=1)
                  if (this%struct(n)%n_.gt.0) then
                     nstruct_=nstruct_+1
                     idmap(n)=stmin+nstruct_-1
                  end if
               end do
               ! Update id array to new index
               update_id: block
                  use amrex_amr_module, only: amrex_mfiter,amrex_box
                  type(amrex_mfiter) :: mfi
                  type(amrex_box) :: bx
                  ! Loop over tiles
                  call this%amr%mfiter_build(lvl,mfi)
                  do while (mfi%next())
                     ! Get pointers to data
                     pid=>this%id%mf(lvl)%dataptr(mfi)
                     ! Perform local loop
                     bx=mfi%tilebox()
                     do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        if (pid(i,j,k,1).gt.0.5_WP) then
                           pid(i,j,k,1)=idmap(nint(pid(i,j,k,1)))
                        end if
                     end do; end do; end do  
                  end do
               end block update_id
               deallocate(idmap)
               ! Finish compacting and renumbering
               allocate(tmp(stmin:stmax))
               nstruct_=0
               do n=1,size(this%struct,dim=1)
                  if (this%struct(n)%n_.gt.0) then
                     nstruct_=nstruct_+1
                     tmp(stmin+nstruct_-1)=this%struct(n)
                  end if
               end do
               call move_alloc(tmp,this%struct)
            end if
         end block compact_tree

         ! Interprocessor treatment of our structures
         interproc_handling: block
            use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k
            integer :: ii,jj,kk,dim
            integer, dimension(3) :: pos
            integer ::stop_global,stop_,counter,n,m,ierr,find_parent,find_parent_own
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid,pidp,pdata
            ! Allocate to total number of structures
            allocate(parent    (nstruct_work)); parent    =0
            allocate(parent_all(nstruct_work)); parent_all=0
            allocate(parent_own(nstruct_work)); parent_own=0
            ! Fill global lineage with selves
            do n=1,nstruct_work
               parent(n)=n
            end do
            ! Synchronize id array
            call this%id%sync()
            ! Loop over cells and check for connections across periodic boundaries, storing parent connections in parent array
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pid=>this%id%mf(lvl)%dataptr(mfi)
               pdata=>data%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Only work with labeled cells 
                  if (pid(i,j,k,1).lt.0.5_WP) cycle
                  ! Loop through one-sided neighbors 
                  do dim=1,3
                     pos=0; pos(dim)=-1
                     ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                     if (pid(ii,jj,kk,1).lt.0.5_WP) cycle
                     ! Check if we should connect these two cells
                     if (same_label(pdata,lbound(pdata),i,j,k,ii,jj,kk)) then
                        ! Update parent array to reflect connection
                        call union_parent(nint(pid(i,j,k,1)),nint(pid(ii,jj,kk,1)))
                     end if
                  end do
               end do; end do; end do
            end do
         
            ! Initialize global stop criterion and counter
            stop_global=1
            counter=0
            ! Resolve lineage
            do while (stop_global.ne.0)
               ! Initialize local stop flag
               stop_=0
               ! Remember own parents
               parent_own=parent
               ! Set self-parents to huge(1)
               do n=1,nstruct_work
                  if (parent(n).eq.n) parent(n)=huge(1)
               end do
               ! Take global min
               call MPI_ALLREDUCE(parent,parent_all,nstruct_work,MPI_INTEGER,MPI_MIN,this%amr%comm,ierr)
               ! Set self-parents back to selves
               do n=1,nstruct_work
                  if (parent_all(n).eq.huge(1)) parent_all(n)=n
               end do
               ! Flatten trees
               do n=1,nstruct_work
                  parent_all(n)=find_all(n)
                  parent_own(n)=find_own(n)
               end do
               ! Start with final parent array being equal to parent_all
               parent=parent_all
               ! Increment counter
               counter=counter+1
               ! Reconcile conflicts between parent_all and parent_own
               do n=1,nstruct_work
                  if (parent_own(n).ne.n) then
                     find_parent_own=rootify_parent(parent_own(n))
                     find_parent    =rootify_parent(parent(n))
                     if (find_parent_own.ne.find_parent) then
                        call union_parent(find_parent,find_parent_own)
                        stop_=1
                     end if
                  end if
               end do
               ! Check if we did some changes
               call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)
            end do
            ! Update this%struct%parent by pointing all parents to root and update id
            ! do n=stmin,stmax
            !    this%struct(n)%parent=rootify_parent(parent(n))
            !    do m=1,this%struct(n)%n_
            !       this%id(this%struct(n)%map(m)%i,this%struct(n)%map(m)%j,this%struct(n)%map(m)%k)=this%struct(n)%parent
            !    end do
            ! end do
            ! Loop over tiles
            call data%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data arrays
               pid=>this%id%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pid(i,j,k,1).gt.0.5_WP) then
                     pid(i,j,k,1)=rootify_parent(parent(nint(pid(i,j,k,1))))
                  end if
               end do; end do; end do
            end do
            ! Update ghost cells with new ids
            call this%id%sync()
         end block interproc_handling

         ! Now we need to compact the data based on id only if on finest level
         renumber_ids: block
            use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER,MPI_IN_PLACE
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k,n,nn,ierr,count
            integer, dimension(:), allocatable :: idmap,counter
            type(struct_type), dimension(:), allocatable :: tmp
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
            ! Only renumber of finest level
            if (.not.finest) exit renumber_ids
            ! Prepare global id map
            allocate(   idmap(1:nstruct_work));    idmap=0
            ! Traverse id array and tag used id values
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pid=>this%id%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pid(i,j,k,1).gt.0) idmap(pid(i,j,k,1))=1
               end do; end do; end do
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,idmap,nstruct_work,MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)
            ! Count number of used structures, set nstruct, and create map
            this%nstruct=sum(idmap)
            count=0
            do n=1,size(idmap,dim=1)
               if (idmap(n).gt.0) then
                  count=count+1
                  idmap(n)=count
               end if
            end do
            ! Rename all structures
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pid=>this%id%mf(lvl)%dataptr(mfi)
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pid(i,j,k,1).gt.0) pid(i,j,k,1)=idmap(pid(i,j,k,1))
               end do; end do; end do
            end do
         end block renumber_ids

         ! Sync final ids
         call this%id%sync()

         ! Release scratch
         call idp%finalize()

         ! Deallocate arrays
         deallocate(parent,parent_all,parent_own)

      end subroutine build_lvl

      !> Debug function to print id's that exist on a level
      subroutine print_ids(lvl,msg)
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         implicit none
         integer, intent(in) :: lvl
         character(len=*), intent(in) :: msg
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer, parameter :: max_id = 100000   ! adjust as needed
         logical :: seen(0:max_id) 
         integer :: count(0:max_id) 
         integer :: id,i,j,k,root

         seen = .false.
         count = 0
         ! Loop over tiles
         call data%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pid => this%id%mf(lvl)%dataptr(mfi)
            bx = mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               id = nint(pid(i,j,k,1))
               if (id <= max_id) then
                  seen(id) = .true.
                  count(id) = count(id) + 1
               else
                  print *, "Warning: ID ", id, " exceeds max_id ", max_id
               end if
            end do; end do; end do
         end do
         ! Collect and print unique IDs
         communicate: block
                  use mpi_f08,   only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_Logical,MPI_LOR, MPI_INTEGER, MPI_SUM
                  integer :: ierr
            call MPI_AllREDUCE(MPI_IN_PLACE,  seen, max_id+1, MPI_LOGICAL, MPI_LOR, this%amr%comm, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, count, max_id+1, MPI_INTEGER, MPI_SUM, this%amr%comm, ierr)
         end block communicate
         if (this%amr%amRoot) then
            print *, "Unique IDs on level ", lvl,' ',msg
            do id=0,max_id
               ! if (seen(id).and.id.gt.0) then
               !    print *,'rootifying on ',id
               !    root = rootify_struct(id)
               ! else
               !    root = 0
               ! end if
               if (seen(id)) print *, 'id = ',id,' count = ',count(id)!, ' root =',root
            end do
         end if
      end subroutine print_ids
            
      
      !> This recursive function that points the lineage of a structure to its root and returns that root
      recursive function rootify_struct(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.this%struct(y)%parent) then
            this%struct(y)%parent=rootify_struct(this%struct(y)%parent)
            y=this%struct(y)%parent
         end if
      end function rootify_struct
      
      !> This function joins two structures at their roots (the smallest root is chosen and returned)
      function union_struct(x,y) result(rmin)
         implicit none
         integer, intent(in) :: x,y
         integer :: rx,ry,rmin,rmax
         rx=rootify_struct(x); ry=rootify_struct(y)
         rmin=min(rx,ry); rmax=max(rx,ry)
         this%struct(rmax)%parent=rmin
      end function union_struct
      
      !> This function adds one new root while dynamically handling storage space
      function add() result(x)
         implicit none
         integer :: x
         integer :: size_now,size_new
         type(struct_type), dimension(:), allocatable :: tmp
         ! Check if there is enough room for storing a new structure
         size_now=size(this%struct,dim=1)
         if (nstruct_.eq.size_now) then
            size_new=nint(real(size_now,WP)*coeff_up)
            allocate(tmp(size_new))
            tmp(1:nstruct_)=this%struct
            tmp(nstruct_+1:)%parent=0
            tmp(nstruct_+1:)%per(1)=0
            tmp(nstruct_+1:)%per(2)=0
            tmp(nstruct_+1:)%per(3)=0
            tmp(nstruct_+1:)%n_=0
            call move_alloc(tmp,this%struct)
         end if
         ! Add new root
         nstruct_=nstruct_+1
         this%struct(nstruct_)%parent=nstruct_
         this%struct(nstruct_)%per=0
         this%struct(nstruct_)%n_=0
         x=nstruct_
      end function add

      !> This subroutine adds an existing root while dynamically handling storage space
      subroutine add_existing(id)
         implicit none
         integer, intent(in) :: id
         integer :: x
         integer :: size_now,size_new
         type(struct_type), dimension(:), allocatable :: tmp
         ! Check if there is enough room for storing a new structure
         size_now=size(this%struct,dim=1)
         if (id.gt.size_now) then
            size_new=id
            allocate(tmp(size_new))
            tmp(1:nstruct_)=this%struct
            tmp(nstruct_+1:)%parent=0
            tmp(nstruct_+1:)%per(1)=0
            tmp(nstruct_+1:)%per(2)=0
            tmp(nstruct_+1:)%per(3)=0
            tmp(nstruct_+1:)%n_=0
            call move_alloc(tmp,this%struct)
         end if
         ! Add new root if doesn't already exist
         if (this%struct(id)%parent.ne.id) then
            nstruct_=nstruct_+1
            this%struct(id)%parent=id
            this%struct(id)%per=0
            this%struct(id)%n_=0
         end if
      end subroutine add_existing
      
      !> This recursive function points global parent to root and returns that root
      recursive function rootify_parent(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.parent(y)) then
            parent(y)=rootify_parent(parent(y))
            y=parent(y)
         end if
      end function rootify_parent
      
      !> This function joins two branches at their roots (the smallest root is chosen)
      subroutine union_parent(x,y)
         implicit none
         integer, intent(in) :: x,y
         integer :: rx,ry,rmin,rmax
         rx=rootify_parent(x); ry=rootify_parent(y); rmin=min(rx,ry); rmax=max(rx,ry)
         parent(rmax)=rmin
      end subroutine union_parent
      
      !> For parent_all array: this function points the parent to root and returns that root
      recursive function find_all(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.parent_all(y)) then
            parent_all(y)=find_all(parent_all(y))
            y=parent_all(y)
         end if
      end function find_all
      
      !> Version of previous function that stops at the completion of a cycle
      recursive function find_all_2(x,x0) result(y)
         implicit none
         integer, intent(in) :: x,x0
         integer :: y
         y=x
         if (y.ne.parent_all(y)) then
            if (parent_all(y).eq.x0) then
               y=parent_all(y)
               return
            else
               parent_all(y)=find_all_2(parent_all(y),x0)
               y=parent_all(y)
            end if
         end if
      end function find_all_2
      
      !> For parent_own array: this function points the parent to root and returns that root
      recursive function find_own(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.parent_own(y)) then
            parent_own(y)=find_own(parent_own(y))
            y=parent_own(y)
         end if
      end function find_own

      subroutine restrict_unique_id(cmf, fmf, ratio, geom)
         use amrex_multifab_module, only : amrex_multifab,amrex_multifab_build,amrex_multifab_destroy
         use amrex_amr_module,   only : amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy
         use amrex_amr_module,   only : amrex_box, amrex_long, amrex_geometry
         use amrex_boxarray_module,   only : amrex_boxarray, amrex_boxarray_build, amrex_boxarray_destroy
         implicit none
         type(amrex_multifab), intent(inout) :: cmf
         type(amrex_multifab), intent(in)    :: fmf
         integer,              intent(in)    :: ratio(3)
         type(amrex_geometry), intent(in)    :: geom

         type(amrex_multifab) :: fine_tmp
         type(amrex_boxarray) :: fba
         type(amrex_mfiter)   :: mfi
         type(amrex_box)      :: bx

         real(WP), contiguous, pointer :: cp(:,:,:,:) => null()
         real(WP), contiguous, pointer :: fp(:,:,:,:) => null()

         integer(amrex_long) :: nb, n
         integer, allocatable :: bxs(:,:,:)   ! (2, 3, nboxes) — lo/hi, dim, box index

         integer, dimension(3) :: clo,chi,flo,fhi

         ! Build refined boxarray by scaling each coarse box's lo/hi
         nb = cmf%ba%nboxes()
         allocate(bxs(2, 3, nb))
         do n = 1, nb
            bx = cmf%ba%get_box(int(n-1))   ! get_box is 0-indexed on the C side
            bxs(1,:,n) = bx%lo * ratio
            bxs(2,:,n) = (bx%hi + 1) * ratio - 1
         end do
         call amrex_boxarray_build(fba, bxs)
         deallocate(bxs)

         call amrex_multifab_build(fine_tmp, fba, cmf%dm, 1, 0)
         call amrex_boxarray_destroy(fba)

         call fine_tmp%setval(0.0_WP)
         call fine_tmp%parallel_copy(fmf, geom)   

         call amrex_mfiter_build(mfi, cmf)
         do while (mfi%next())
            bx =  mfi%validbox()
            cp => cmf%dataptr(mfi)
            fp => fine_tmp%dataptr(mfi)

            clo = [lbound(cp,1), lbound(cp,2), lbound(cp,3)]
            chi = [ubound(cp,1), ubound(cp,2), ubound(cp,3)]
            flo = [lbound(fp,1), lbound(fp,2), lbound(fp,3)]
            fhi = [ubound(fp,1), ubound(fp,2), ubound(fp,3)]


            call restrict_kernel(cp(:,:,:,1), clo, chi, &
                     fp(:,:,:,1), flo, fhi, &
                     bx%lo, bx%hi, ratio)

            nullify(cp, fp)
         end do
         call amrex_mfiter_destroy(mfi)
         call amrex_multifab_destroy(fine_tmp)

      end subroutine restrict_unique_id

      !---------------------------------------------------------------------------
      ! Private kernel — operates on a single patch
      !---------------------------------------------------------------------------
      subroutine restrict_kernel(crse, clo, chi, fine, flo, fhi, lo, hi, ratio)
         implicit none
         integer,  intent(in)    :: clo(3), chi(3)
         integer,  intent(in)    :: flo(3), fhi(3)
         integer,  intent(in)    :: lo(3), hi(3), ratio(3)
         real(WP), intent(inout) :: crse(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
         real(WP), intent(in)    :: fine(flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3))

         integer :: i,  j,  k
         integer :: ii, jj, kk
         integer :: id_val, id_store
         logical :: found, conflict

         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)

               found    = .false.
               conflict = .false.
               id_store = 0

               do kk = k*ratio(3), k*ratio(3) + ratio(3) - 1
               do jj = j*ratio(2), j*ratio(2) + ratio(2) - 1
               do ii = i*ratio(1), i*ratio(1) + ratio(1) - 1

                  if (abs(fine(ii,jj,kk)) > 0.5_WP) then
                     id_val = nint(fine(ii,jj,kk))
                     if (.not. found) then
                           id_store = id_val
                           found    = .true.
                     else if (id_val /= id_store) then
                           conflict = .true.
                     end if
                  end if

               end do
               end do
               end do

               if (found .and. .not. conflict) then
                  crse(i,j,k) = real(id_store, WP)
               else
                  crse(i,j,k) = 0.0_WP
               end if

         end do
         end do
         end do

      end subroutine restrict_kernel
      
   end subroutine build
   
   
   !> Empty structure info
   subroutine empty(this)
      implicit none
      class(amrcclabel), intent(inout) :: this
      integer :: n
      ! Deallocate structure array
      if (allocated(this%struct)) deallocate(this%struct)
   end subroutine empty


   !> Compute common statistics for structures 
   !> identified by id in this%id and weighted by VF array 
   subroutine compute_stats(this,VF,stats)
      use amrex_fort_module,     only : amrex_spacedim
      use amrex_multifab_module, only : amrex_multifab, amrex_mfiter, &
                                        amrex_mfiter_build, amrex_mfiter_destroy
      use amrex_box_module,      only : amrex_box
      use amrex_boxarray_module, only : amrex_boxarray, amrex_boxarray_build, amrex_boxarray_destroy
      use amrex_geometry_module, only : amrex_geometry
      use amrex_box_module,      only : amrex_box
      use amrex_amr_module,      only : amrex_long
      use amrex_parallel_module, only : amrex_parallel_reduce_sum
      implicit none
      class(amrcclabel) :: this
      type(amrdata), intent(in) :: VF
      type(stats_type), allocatable, dimension(:), intent(out) :: stats
      real(WP), allocatable, dimension(:)   :: vol_map
      real(WP), allocatable, dimension(:,:) :: com_map
      logical :: id_seen(this%nstruct)
      type(amrex_mfiter)   :: mfi
      type(amrex_box)      :: bx
      type(amrex_boxarray) :: fine_ba_crse
      type(amrex_box) :: pt_box
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pid,pVF
      real(WP) :: dx(3), cell_vol, prob_lo(3)
      real(WP) :: xc, yc, zc
      integer  :: lo(3), hi(3), ilo(3), ihi(3)
      integer  :: i, j, k, lvl, id_val
      real(WP) :: VF_val
      integer(amrex_long) :: nb, n
      integer, allocatable :: bxs(:,:,:)
      integer :: ratio(3)

      allocate(vol_map(1:this%nstruct))
      allocate(com_map(1:this%nstruct,3))

      vol_map = 0.0_WP
      com_map = 0.0_WP
      prob_lo = [this%amr%xlo, this%amr%ylo, this%amr%zlo]

      do lvl = 0, this%amr%maxlvl

         dx(1)    = this%amr%dx(lvl)
         dx(2)    = this%amr%dy(lvl)
         dx(3)    = this%amr%dz(lvl)
         cell_vol = this%amr%cell_vol(lvl)

         if (lvl < this%amr%maxlvl) then
               ratio = [this%amr%rrefx(lvl), this%amr%rrefy(lvl), this%amr%rrefz(lvl)]
               nb = this%id%mf(lvl+1)%ba%nboxes()
               allocate(bxs(2, 3, nb))
               do n = 1, nb
                  bx = this%id%mf(lvl+1)%ba%get_box(int(n-1))
                  bxs(1,:,n) = bx%lo / ratio
                  bxs(2,:,n) = bx%hi / ratio
               end do
               call amrex_boxarray_build(fine_ba_crse, bxs)
               deallocate(bxs)
         end if

         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data arrays
            pid => this%id%mf(lvl)%dataptr(mfi)   
            pVF =>      VF%mf(lvl)%dataptr(mfi)
            bx  = mfi%validbox()
            lo  = bx%lo
            hi  = bx%hi
            
            ilo = [lbound(pid,1), lbound(pid,2), lbound(pid,3)]
            ihi = [ubound(pid,1), ubound(pid,2), ubound(pid,3)]

            accumulate_stats: block 
               real(WP) :: id_arr(ilo(1):ihi(1), ilo(2):ihi(2), ilo(3):ihi(3))

               do k = lo(3), hi(3)
               do j = lo(2), hi(2)
               do i = lo(1), hi(1)

                  id_val = nint(pid(i,j,k,1))
                  VF_val =      pVF(i,j,k,1)
                  if (id_val <= 0) cycle

                  if (lvl < this%amr%maxlvl) then
                     pt_box%lo = [i, j, k]
                     pt_box%hi = [i, j, k]
                     if (fine_ba_crse%intersects(pt_box)) cycle
                  end if

                  xc = prob_lo(1) + (real(i, WP) + 0.5_WP) * dx(1)
                  yc = prob_lo(2) + (real(j, WP) + 0.5_WP) * dx(2)
                  zc = prob_lo(3) + (real(k, WP) + 0.5_WP) * dx(3)

                  vol_map(id_val  ) = vol_map(id_val  ) + cell_vol * VF_val 
                  com_map(id_val,1) = com_map(id_val,1) + cell_vol * VF_val * xc
                  com_map(id_val,2) = com_map(id_val,2) + cell_vol * VF_val * yc
                  com_map(id_val,3) = com_map(id_val,3) + cell_vol * VF_val * zc

               end do
               end do
               end do
            end block accumulate_stats

            nullify(pid)
            nullify(pVF)
         end do
         call amrex_mfiter_destroy(mfi)

         if (lvl < this%amr%maxlvl) call amrex_boxarray_destroy(fine_ba_crse)

      end do

      call amrex_parallel_reduce_sum(vol_map, this%nstruct)
      call amrex_parallel_reduce_sum(com_map(:,1), this%nstruct)
      call amrex_parallel_reduce_sum(com_map(:,2), this%nstruct)
      call amrex_parallel_reduce_sum(com_map(:,3), this%nstruct)

      allocate(stats(this%nstruct))
      do i = 1, this%nstruct
         stats(i)%id  = i
         stats(i)%vol = vol_map(i)
         stats(i)%com = com_map(i,:)/vol_map(i)
      end do

   end subroutine compute_stats
   
   
   !> Finalize CCL object
   subroutine finalize(this)
      implicit none
      class(amrcclabel), intent(inout) :: this
      call this%empty()
      call this%id%finalize()
      ! nullify(this%pg)
      this%name='UNNAMED_CCL'
   end subroutine finalize
   
   
end module amrcclabel_class
