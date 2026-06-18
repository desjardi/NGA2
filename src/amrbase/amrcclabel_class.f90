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
    integer  :: id           ! Structure ID
    real(WP) :: vol          ! Liquid volume
    real(WP) :: com(3)       ! Center of mass
    real(WP) :: vel(3)       ! Volume-weighted liquid velocity
    real(WP) :: moi(3,3)     ! Moment of inertia tensor
    real(WP) :: Deq          ! Equivalent sphere diameter
    real(WP) :: gvel(3)      ! Volume-weighted surrounding gas velocity
    real(WP) :: weber        ! Weber number
    logical  :: remove       ! .true. if structure touches domain boundary
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
      call this%id%reset() ! Update with current grids
      call this%id%setval(0.0_WP)
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
            use mpi_f08, only: MPI_ALLREDUCE,MPI_INTEGER,MPI_MAX
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k
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
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
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
            integer ::stop_global,stop_,counter,n,ierr,find_parent,find_parent_own
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid,pdata
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
            integer :: i,j,k,n,ierr,count
            integer, dimension(:), allocatable :: idmap
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
                  if (pid(i,j,k,1).gt.0) idmap(nint(pid(i,j,k,1)))=1
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
                  if (pid(i,j,k,1).gt.0) pid(i,j,k,1)=idmap(nint(pid(i,j,k,1)))
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
      ! Deallocate structure array
      if (allocated(this%struct)) deallocate(this%struct)
   end subroutine empty


   !> Compute common statistics for structures 
   !> identified by id in this%id
   subroutine compute_stats(this, VF, Q, rhoG, sigma, stats)
      use amrex_fort_module,     only : amrex_spacedim
      use amrex_multifab_module, only : amrex_multifab, amrex_mfiter, &
                                       amrex_mfiter_build, amrex_mfiter_destroy
      use amrex_box_module,      only : amrex_box
      use amrex_boxarray_module, only : amrex_boxarray, amrex_boxarray_build, amrex_boxarray_destroy
      use amrex_parallel_module, only : amrex_parallel_reduce_sum
      use amrex_amr_module, only : amrex_long
      use mathtools,             only : pi
      implicit none
      class(amrcclabel) :: this
      type(amrdata), intent(in) :: VF     ! Volume fraction
      type(amrdata), intent(in) :: Q      ! Cell-centred velocity (3 components)
      real(WP),      intent(in) :: rhoG   ! Gas density
      real(WP),      intent(in) :: sigma  ! Surface tension coefficient
      type(stats_type), allocatable, dimension(:), intent(out) :: stats

      ! Accumulator arrays
      real(WP), allocatable :: vol_map(:)      ! (nstruct)
      real(WP), allocatable :: com_map(:,:)    ! (nstruct, 3)
      real(WP), allocatable :: vel_map(:,:)    ! (nstruct, 3) liquid velocity
      real(WP), allocatable :: moi_map(:,:,:)  ! (nstruct, 3, 3)
      real(WP), allocatable :: gvel_map(:,:)   ! (nstruct, 3) gas velocity
      real(WP), allocatable :: gwt_map(:)      ! (nstruct) gas velocity weights
      real(WP), allocatable :: rem_map(:)      ! (nstruct) boundary flag (real for MPI reduce)

      type(amrex_mfiter)   :: mfi
      type(amrex_box)      :: bx, pt_box
      type(amrex_boxarray) :: fine_ba_crse

      real(WP), pointer, contiguous :: pid(:,:,:,:), pVF(:,:,:,:), pQ(:,:,:,:)
      real(WP) :: dx(3), cell_vol, prob_lo(3)
      real(WP) :: xc, yc, zc, VF_val
      real(WP) :: xr, yr, zr, x0, y0, z0
      real(WP) :: slip_vel, Deq
      integer  :: lo(3), hi(3), ilo(3), ihi(3)
      integer  :: i, j, k, lvl, id_val, ns
      integer  :: ratio(3)
      integer, parameter :: nlayer = 2         ! cells near domain face to flag
      integer(amrex_long) :: nb, nn
      integer, allocatable :: bxs(:,:,:)

      allocate(vol_map (this%nstruct));         vol_map  = 0.0_WP
      allocate(com_map (this%nstruct, 3));      com_map  = 0.0_WP
      allocate(vel_map (this%nstruct, 3));      vel_map  = 0.0_WP
      allocate(moi_map (this%nstruct, 3, 3));   moi_map  = 0.0_WP
      allocate(gvel_map(this%nstruct, 3));      gvel_map = 0.0_WP
      allocate(gwt_map (this%nstruct));         gwt_map  = 0.0_WP
      allocate(rem_map (this%nstruct));         rem_map  = 0.0_WP

      prob_lo = [this%amr%xlo, this%amr%ylo, this%amr%zlo]

      ! =========================================================================
      ! PASS 1: volume, CoM, liquid velocity, boundary removal flag
      ! =========================================================================
      do lvl = 0, this%amr%clvl()

         dx(1)    = this%amr%dx(lvl)
         dx(2)    = this%amr%dy(lvl)
         dx(3)    = this%amr%dz(lvl)
         cell_vol = this%amr%cell_vol(lvl)

         if (lvl < this%amr%clvl()) then
            ratio = [this%amr%rrefx(lvl), this%amr%rrefy(lvl), this%amr%rrefz(lvl)]
            nb = this%id%mf(lvl+1)%ba%nboxes()
            allocate(bxs(2, 3, nb))
            do nn = 1, nb
               bx = this%id%mf(lvl+1)%ba%get_box(int(nn-1))
               bxs(1,:,nn) = bx%lo / ratio
               bxs(2,:,nn) = bx%hi / ratio
            end do
            call amrex_boxarray_build(fine_ba_crse, bxs)
            deallocate(bxs)
         end if

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            pid => this%id%mf(lvl)%dataptr(mfi)
            pVF =>      VF%mf(lvl)%dataptr(mfi)
            pQ  =>       Q%mf(lvl)%dataptr(mfi)
            bx  = mfi%validbox()
            lo  = bx%lo;  hi = bx%hi

            pass1: block
               do k = lo(3), hi(3)
               do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  id_val = nint(pid(i,j,k,1))
                  VF_val =      pVF(i,j,k,1)
                  if (id_val <= 0 .or. id_val > this%nstruct) cycle

                  if (lvl < this%amr%clvl()) then
                        pt_box%lo = [i,j,k];  pt_box%hi = [i,j,k]
                        if (fine_ba_crse%intersects(pt_box)) cycle
                  end if

                  xc = prob_lo(1) + (real(i, WP) + 0.5_WP) * dx(1)
                  yc = prob_lo(2) + (real(j, WP) + 0.5_WP) * dx(2)
                  zc = prob_lo(3) + (real(k, WP) + 0.5_WP) * dx(3)

                  vol_map(id_val)   = vol_map(id_val)   + cell_vol * VF_val
                  com_map(id_val,1) = com_map(id_val,1) + cell_vol * VF_val * xc
                  com_map(id_val,2) = com_map(id_val,2) + cell_vol * VF_val * yc
                  com_map(id_val,3) = com_map(id_val,3) + cell_vol * VF_val * zc
                  vel_map(id_val,1) = vel_map(id_val,1) + cell_vol * VF_val * pQ(i,j,k,1)
                  vel_map(id_val,2) = vel_map(id_val,2) + cell_vol * VF_val * pQ(i,j,k,2)
                  vel_map(id_val,3) = vel_map(id_val,3) + cell_vol * VF_val * pQ(i,j,k,3)

                  ! Flag if within nlayer cells of any domain face
                  if (xc < this%amr%xlo + nlayer*dx(1) .or. &
                        xc > this%amr%xhi - nlayer*dx(1) .or. &
                        yc < this%amr%ylo + nlayer*dx(2) .or. &
                        yc > this%amr%yhi - nlayer*dx(2) .or. &
                        zc < this%amr%zlo + nlayer*dx(3) .or. &
                        zc > this%amr%zhi - nlayer*dx(3)) then
                        rem_map(id_val) = 1.0_WP
                  end if
               end do
               end do
               end do
            end block pass1

            nullify(pid, pVF, pQ)
         end do
         call amrex_mfiter_destroy(mfi)

         if (lvl < this%amr%clvl()) call amrex_boxarray_destroy(fine_ba_crse)
      end do

      ! Reduce pass 1
      call amrex_parallel_reduce_sum(vol_map,      this%nstruct)
      call amrex_parallel_reduce_sum(com_map(:,1), this%nstruct)
      call amrex_parallel_reduce_sum(com_map(:,2), this%nstruct)
      call amrex_parallel_reduce_sum(com_map(:,3), this%nstruct)
      call amrex_parallel_reduce_sum(vel_map(:,1), this%nstruct)
      call amrex_parallel_reduce_sum(vel_map(:,2), this%nstruct)
      call amrex_parallel_reduce_sum(vel_map(:,3), this%nstruct)
      call amrex_parallel_reduce_sum(rem_map,      this%nstruct)

      ! Normalize CoM and velocity
      do ns = 1, this%nstruct
         if (vol_map(ns) > 0.0_WP) then
               com_map(ns,:) = com_map(ns,:) / vol_map(ns)
               vel_map(ns,:) = vel_map(ns,:) / vol_map(ns)
         end if
      end do

      ! =========================================================================
      ! PASS 2: moment of inertia + surrounding gas velocity
      ! Ghost cells on id%mf are filled by sync() at end of build(),
      ! so neighbor IDs are visible across FAB boundaries.
      ! =========================================================================
      do lvl = 0, this%amr%clvl()

         dx(1)    = this%amr%dx(lvl)
         dx(2)    = this%amr%dy(lvl)
         dx(3)    = this%amr%dz(lvl)
         cell_vol = this%amr%cell_vol(lvl)

         if (lvl < this%amr%clvl()) then
               ratio = [this%amr%rrefx(lvl), this%amr%rrefy(lvl), this%amr%rrefz(lvl)]
               nb = this%id%mf(lvl+1)%ba%nboxes()
               allocate(bxs(2, 3, nb))
               do nn = 1, nb
                  bx = this%id%mf(lvl+1)%ba%get_box(int(nn-1))
                  bxs(1,:,nn) = bx%lo / ratio
                  bxs(2,:,nn) = bx%hi / ratio
               end do
               call amrex_boxarray_build(fine_ba_crse, bxs)
               deallocate(bxs)
         end if

         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
               pid => this%id%mf(lvl)%dataptr(mfi)
               pVF =>      VF%mf(lvl)%dataptr(mfi)
               pQ  =>       Q%mf(lvl)%dataptr(mfi)
               bx  = mfi%validbox()
               lo  = bx%lo;  hi = bx%hi
               ilo = [lbound(pid,1), lbound(pid,2), lbound(pid,3)]
               ihi = [ubound(pid,1), ubound(pid,2), ubound(pid,3)]

               pass2: block
                  real(WP) :: id_arr(ilo(1):ihi(1), ilo(2):ihi(2), ilo(3):ihi(3))
                  integer  :: unique_ids(6), n_unique, d, nbid, ii, jj, kk
                  integer, dimension(3,6) :: off

                  ! Copy with correct AMReX bounds so neighbor lookup works
                  id_arr = pid(:,:,:,1)

                  off(:,1)=[1,0,0]; off(:,2)=[-1,0,0]
                  off(:,3)=[0,1,0]; off(:,4)=[0,-1,0]
                  off(:,5)=[0,0,1]; off(:,6)=[0,0,-1]

                  do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     id_val = nint(id_arr(i,j,k))
                     VF_val = pVF(i,j,k,1)

                     if (lvl < this%amr%clvl()) then
                           pt_box%lo = [i,j,k];  pt_box%hi = [i,j,k]
                           if (fine_ba_crse%intersects(pt_box)) cycle
                     end if

                     xc = prob_lo(1) + (real(i, WP) + 0.5_WP) * dx(1)
                     yc = prob_lo(2) + (real(j, WP) + 0.5_WP) * dx(2)
                     zc = prob_lo(3) + (real(k, WP) + 0.5_WP) * dx(3)

                     ! --- Moment of inertia for liquid cells ---
                     if (id_val > 0 .and. id_val <= this%nstruct) then
                           x0 = com_map(id_val,1);  xr = xc - x0
                           y0 = com_map(id_val,2);  yr = yc - y0
                           z0 = com_map(id_val,3);  zr = zc - z0
                           moi_map(id_val,1,1) = moi_map(id_val,1,1) + cell_vol*VF_val*(yr**2+zr**2)
                           moi_map(id_val,2,2) = moi_map(id_val,2,2) + cell_vol*VF_val*(zr**2+xr**2)
                           moi_map(id_val,3,3) = moi_map(id_val,3,3) + cell_vol*VF_val*(xr**2+yr**2)
                           moi_map(id_val,1,2) = moi_map(id_val,1,2) - cell_vol*VF_val*(xr*yr)
                           moi_map(id_val,1,3) = moi_map(id_val,1,3) - cell_vol*VF_val*(xr*zr)
                           moi_map(id_val,2,3) = moi_map(id_val,2,3) - cell_vol*VF_val*(yr*zr)
                     end if

                     ! --- Gas velocity: accumulate gas cell adjacent to structures ---
                     ! Ghost-cell-filled id_arr lets us see neighbor IDs across FABs.
                     ! Collect unique structure IDs from 6-connected neighbors to avoid
                     ! double-counting a gas cell that touches multiple cells of the
                     ! same structure.
                     if (VF_val < 0.5_WP) then
                           unique_ids = 0;  n_unique = 0
                           do d = 1, 6
                              ii = i + off(1,d)
                              jj = j + off(2,d)
                              kk = k + off(3,d)
                              if (ii < ilo(1) .or. ii > ihi(1) .or. &
                                 jj < ilo(2) .or. jj > ihi(2) .or. &
                                 kk < ilo(3) .or. kk > ihi(3)) cycle
                              nbid = nint(id_arr(ii,jj,kk))
                              if (nbid <= 0 .or. nbid > this%nstruct) cycle
                              if (any(unique_ids(1:n_unique) == nbid)) cycle
                              n_unique = n_unique + 1
                              unique_ids(n_unique) = nbid
                           end do
                           do d = 1, n_unique
                              nbid = unique_ids(d)
                              gvel_map(nbid,1) = gvel_map(nbid,1) + cell_vol*(1.0_WP-VF_val)*pQ(i,j,k,1)
                              gvel_map(nbid,2) = gvel_map(nbid,2) + cell_vol*(1.0_WP-VF_val)*pQ(i,j,k,2)
                              gvel_map(nbid,3) = gvel_map(nbid,3) + cell_vol*(1.0_WP-VF_val)*pQ(i,j,k,3)
                              gwt_map(nbid)    = gwt_map(nbid)    + cell_vol*(1.0_WP-VF_val)
                           end do
                     end if

                  end do
                  end do
                  end do
               end block pass2

               nullify(pid, pVF, pQ)
         end do
         call amrex_mfiter_destroy(mfi)

         if (lvl < this%amr%clvl()) call amrex_boxarray_destroy(fine_ba_crse)
      end do

      ! Reduce pass 2
      call amrex_parallel_reduce_sum(moi_map(:,1,1), this%nstruct)
      call amrex_parallel_reduce_sum(moi_map(:,2,2), this%nstruct)
      call amrex_parallel_reduce_sum(moi_map(:,3,3), this%nstruct)
      call amrex_parallel_reduce_sum(moi_map(:,1,2), this%nstruct)
      call amrex_parallel_reduce_sum(moi_map(:,1,3), this%nstruct)
      call amrex_parallel_reduce_sum(moi_map(:,2,3), this%nstruct)
      call amrex_parallel_reduce_sum(gvel_map(:,1),  this%nstruct)
      call amrex_parallel_reduce_sum(gvel_map(:,2),  this%nstruct)
      call amrex_parallel_reduce_sum(gvel_map(:,3),  this%nstruct)
      call amrex_parallel_reduce_sum(gwt_map,        this%nstruct)

      ! =========================================================================
      ! Pack results
      ! =========================================================================
      allocate(stats(this%nstruct))
      do ns = 1, this%nstruct
         stats(ns)%id  = ns
         stats(ns)%vol = vol_map(ns)
         stats(ns)%com = com_map(ns,:)
         stats(ns)%vel = vel_map(ns,:)
         stats(ns)%remove = (rem_map(ns) > 0.0_WP)

         ! Fill symmetric off-diagonal components
         stats(ns)%moi      = moi_map(ns,:,:)
         stats(ns)%moi(2,1) = moi_map(ns,1,2)
         stats(ns)%moi(3,1) = moi_map(ns,1,3)
         stats(ns)%moi(3,2) = moi_map(ns,2,3)

         ! Equivalent diameter from liquid volume
         Deq = ((vol_map(ns) * 6.0_WP) / pi)**(1.0_WP/3.0_WP)
         stats(ns)%Deq = Deq

         ! Surrounding gas velocity
         if (gwt_map(ns) > 0.0_WP) then
               stats(ns)%gvel = gvel_map(ns,:) / gwt_map(ns)
         else
               stats(ns)%gvel = 0.0_WP
         end if

         ! Weber number: rhoG * |slip|^2 * Deq / sigma
         slip_vel = sqrt(sum((stats(ns)%gvel - stats(ns)%vel)**2))
         if (sigma > 0.0_WP) then
               stats(ns)%weber = rhoG * slip_vel**2 * Deq / sigma
         else
               stats(ns)%weber = 0.0_WP
         end if
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
