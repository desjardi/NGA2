!> TODO
! - restict seems to be working
! - Now need to update cells on a coarse level that are completely liquid - i think.  Test with level = 3 or 4?


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
   public :: amrcclabel,make_label_ftype,same_label_ftype
   
   
   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%
   
   !> Structure object
   type :: struct_type
      integer :: parent                                   !< ID of parent struct
      integer :: n_                                       !< Number of local cells contained in struct
      integer, dimension(3) :: per                        !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct_type
   
   
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

      testing_finest: block 
         integer :: lvl
         do lvl = 0,data%amr%maxlvl
            call print_ids(lvl,"after build_lvl(finest)")
         end do
      end block testing_finest

      ! Create unique IDs for each structure on coarser levels
      restrict_unique_id: block
         use amrex_interface, only: amrmfab_restrict_unique_id         
         integer :: lvl
         integer, dimension(3) :: ref_ratio
         do lvl = data%amr%maxlvl-1, 0, -1   ! finest-1 → coarsest
            if (this%amr%amRoot) print *,'Restricting to level ',lvl

            ref_ratio(1)=data%amr%rrefx(lvl)
            ref_ratio(2)=data%amr%rrefy(lvl)
            ref_ratio(3)=data%amr%rrefz(lvl)

            ! Implemented in C to get access to additional functions
            call amrmfab_restrict_unique_id( &
               this%id%mf(lvl),     & ! coarse
               this%id%mf(lvl+1),   & ! fine
               ref_ratio )

            testing_after_restrict: block 
               integer :: lvl
               do lvl = 0,data%amr%maxlvl
                  call print_ids(lvl,"after restrict")
               end do
            end block testing_after_restrict

            ! Build CCL on coarse level
            call build_lvl(lvl,coarse_make_label,coarse_same_label)

         end do
      end block restrict_unique_id

      testing_end_build: block 
         integer :: lvl
         do lvl = 0,data%amr%maxlvl
            call print_ids(lvl,"after build")
         end do
      end block testing_end_build

   contains

      !> Build structure on a level using user-set test functions
      subroutine build_lvl(lvl,make_label,same_label)
         integer, intent(in) :: lvl
         procedure(make_label_ftype) :: make_label
         procedure(same_label_ftype) :: same_label

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

         ! Add any ids from coarser levels to struct array
         previous_ids: block 
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            integer :: i,j,k
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
            ! Loop over tiles
            call data%amr%mfiter_build(lvl,mfi)
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
            call data%amr%mfiter_build(lvl,mfi)
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
                                 ! print *,'Using neighbor''s label'
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
                     ! if (nint(pid(i,j,k,1)).ne.1) then
                     !    print *,' collapsing ',pid(i,j,k,1),' into ',rootify_struct(nint(pid(i,j,k,1)))
                     ! end if
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
            use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
            integer :: i,j,k,n,ierr
            integer, dimension(:), allocatable :: my_nstruct,all_nstruct,idmap
            type(struct_type), dimension(:), allocatable :: tmp
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pid
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
            this%nstruct=sum(all_nstruct)
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
            allocate(parent    (this%nstruct)); parent    =0
            allocate(parent_all(this%nstruct)); parent_all=0
            allocate(parent_own(this%nstruct)); parent_own=0
            ! Fill global lineage with selves
            do n=1,this%nstruct
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
               do n=1,this%nstruct
                  if (parent(n).eq.n) parent(n)=huge(1)
               end do
               ! Take global min
               call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%amr%comm,ierr)
               ! Set self-parents back to selves
               do n=1,this%nstruct
                  if (parent_all(n).eq.huge(1)) parent_all(n)=n
               end do
               ! Flatten trees
               do n=1,this%nstruct
                  parent_all(n)=find_all(n)
                  parent_own(n)=find_own(n)
               end do
               ! Start with final parent array being equal to parent_all
               parent=parent_all
               ! Increment counter
               counter=counter+1
               ! Reconcile conflicts between parent_all and parent_own
               do n=1,this%nstruct
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
         integer :: id,i,j,k

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
               if (seen(id)) print *, 'id = ',id,' count = ',count(id)
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
      
   end subroutine build
   
   
   !> Empty structure info
   subroutine empty(this)
      implicit none
      class(amrcclabel), intent(inout) :: this
      integer :: n
      ! Deallocate structure array
      if (allocated(this%struct)) deallocate(this%struct)
      ! Zero structures
      this%nstruct=0
   end subroutine empty
   
   
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
