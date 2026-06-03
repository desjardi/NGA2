!> TODO
! - How to deal with what used to be domain boundaries 
!   if (this%pg%imin_.ne.this%pg%imin) then ! ??????????????
! - How to deal with VF? How to access values?
! - How should map be represented? If it should.
!    old - i,j,k
!    new - level, tile(?), i,j,k

!> Connected component labeling class: identifies Lagrangian objects from a Eulerian logical field
!> and provides unstructured mapping to traverse these objects
module amrcclabel_class
   use precision,       only: WP
   use string,          only: str_medium
   use amrdata_class,   only: amrdata
   use amrsolver_class, only: amrsolver
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: amrcclabel,make_label_ftype,same_label_ftype
   
   
   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%

   !> Map object
   type :: map_type
      integer :: lvl !< AMR level
      integer :: fab !< mfi%index()
      integer, dimension(:), allocatable :: i,j,k !< Cell index
   end type map_type
   
   
   !> Structure object
   type :: struct_type
      integer :: parent                                   !< ID of parent struct
      integer :: n_                                       !< Number of local cells contained in struct
      type(map_type), dimension(:), allocatable :: map         !< List of cells contained in struct
      integer, dimension(3) :: per                        !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct_type
   
   
   !> amrcclabel object definition
   type, extends(amrsolver) :: amrcclabel
      ! ID of the structure that contains each cell
      type(amrdata) :: id
      ! Periodicity treatement
      type(amrdata) :: idp
      ! Array of structures
      integer :: nstruct
      type(struct_type), dimension(:), allocatable :: struct
      ! Ghost cells
      integer :: nover=1
   contains
      procedure :: initialize
      procedure :: build
      procedure :: empty
      procedure :: finalize
   end type amrcclabel
   
   !> Type of the make_label function used to generate a structure
   interface
      logical function make_label_ftype(pVF,i,j,k)
         use precision,    only: WP
         real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pVF
         integer, intent(in) :: i,j,k
      end function make_label_ftype
   end interface
   
   !> Type of the same_label function used to connect two structures
   interface
      logical function same_label_ftype(pVF,i,j,k,ii,jj,kk)
         use precision,    only: WP
         real(WP), dimension(:,:,:,:), contiguous, intent(in) :: pVF
         integer, intent(in) :: i,j,k,ii,jj,kk
      end function same_label_ftype
   end interface
   
   
contains
   
   
   !> Initialization for amrcclabel class
   subroutine initialize(this,name)
      implicit none
      class(amrcclabel) :: this
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Allocate and initialize ID array
      call this%id%initialize(amr,name='id',ncomp=1,ng=this%nover); this%id%parent=>this
      call this%id%setval(val=0)
      ! Allocate and initialize periodicity array
      call this%idp%initialize(amr,name='idp',ncomp=3,ng=this%nover); this%idp%parent=>this
      ! Zero structures
      this%nstruct=0
   end subroutine initialize
   
   
   !> Build structure using the user-set test functions
   subroutine build(this,make_label,same_label)
      use amrdata_class,    only: amrdata
      implicit none
      class(amrcclabel), intent(inout) :: this
      procedure(make_label_ftype) :: make_label
      procedure(same_label_ftype) :: same_label
      integer :: nstruct_,stmin,stmax
      integer, dimension(:), allocatable :: parent             !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_all         !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_own         !< Resolving structure id across procs
      
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
      
      ! Allocate periodicity work array
      call this%idp%setval(val=0)
      
      ! Perform a first pass to build proc-local structures and corresponding tree
      first_pass: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         integer :: lvl,i,j,k
         integer :: ii,jj,kk,dim
         integer :: fab 
         integer, dimension(3) :: pos
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF

         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pid=>this%id%mf(lvl)%dataptr(mfi)
               pidp=>this%idp%mf(lvl)%dataptr(mfi)
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               ! Only work on finest level for now
               if (lvl.ne.this%amr%finest_level()) cycle
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Find next cell in a structure
                  if (make_label(pVF,i,j,k)) then
                     ! Loop through one-sided neighbors
                     do dim=1,3
                        pos=0; pos(dim)=-1
                        ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                        ! Check if neighbor is labeled
                        if (pid(ii,jj,kk,1).gt.0) then
                           ! Neighbor is labeled, but are we?
                           if (pid(i,j,k,1).ne.0) then
                              ! We already have a label, perform a union of both labels
                              if (same_label(pVF,i,j,k,ii,jj,kk)) then
                                 pid(i,j,k,1)=union_struct(pid(i,j,k,1),pid(ii,jj,kk,1))
                              end if
                           else
                              ! We don't have a label, check if we take the neighbor's label
                              if (same_label(pVF,i,j,k,ii,jj,kk)) then
                                 pid(i,j,k,1)=pid(ii,jj,kk,1)
                              else
                                 pid(i,j,k,1)=add()
                              end if
                           end if
                        end if
                     end do
                     ! If no neighbor was labeled, we need a new structure
                     if (pid(i,j,k,1).eq.0) pid(i,j,k,1)=add()
                     ! Identify periodicity cases
                     if (this%pg%xper.and.i.eq.this%pg%imax) this%struct(pid(i,j,k,1))%per(1)=1
                     if (this%pg%yper.and.j.eq.this%pg%jmax) this%struct(pid(i,j,k,1))%per(2)=1
                     if (this%pg%zper.and.k.eq.this%pg%kmax) this%struct(pid(i,j,k,1))%per(3)=1
                     pidp(i,j,k,:)=this%struct(pid(i,j,k,1))%per
                  end if
               end do; end do; end do
            end do
         end do
      end block first_pass
      
      ! Now collapse the tree, count the cells and resolve periodicity in each structure
      collapse_tree: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         integer :: i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         do lvl=0,this%amr%clvl()
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pid=>this%id%mf(lvl)%dataptr(mfi)
               pidp=>this%idp%mf(lvl)%dataptr(mfi)
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               ! Only work on finest level for now
               if (lvl.ne.this%amr%finest_level()) cycle
               ! Perform local loop
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pid(i,j,k,1).gt.0) then
                     pid(i,j,k,1)=rootify_struct(pid(i,j,k,1))
                     this%struct(pid(i,j,k,1))%n_=this%struct(pid(i,j,k,1))%n_+1
                     pidp(i,j,k,1)=max(pidp(1,i,j,k),this%struct(pid(i,j,k,1))%per(1))
                     pidp(i,j,k,2)=max(pidp(2,i,j,k),this%struct(pid(i,j,k,1))%per(2))
                     pidp(i,j,k,3)=max(pidp(3,i,j,k),this%struct(pid(i,j,k,1))%per(3))
                     this%struct(pid(i,j,k,1))%per=pidp(:,i,j,k)
                  end if
               end do; end do; end do
            end do
         end do
      end block collapse_tree
      
      ! Compact structure array
      compact_tree: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
         integer :: i,j,k,n,ierr
         integer, dimension(:), allocatable :: my_nstruct,all_nstruct,idmap
         type(struct_type), dimension(:), allocatable :: tmp
         ! Count exact number of local structures
         nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%n_.gt.0) nstruct_=nstruct_+1
         end do
         ! Gather this info to ensure unique index
         allocate( my_nstruct(0:this%pg%nproc-1)); my_nstruct=0; my_nstruct(this%pg%rank)=nstruct_
         allocate(all_nstruct(0:this%pg%nproc-1)); call MPI_ALLREDUCE(my_nstruct,all_nstruct,this%pg%nproc,MPI_INTEGER,MPI_SUM,this%pg%comm,ierr)
         stmin=1
         if (this%pg%rank.gt.0) stmin=stmin+sum(all_nstruct(0:this%pg%rank-1))
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
            integer :: lvl
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            ! Traverse levels ! Only work on finest level for now
            do lvl=this%amr%finest_level()
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
               allocate(tmp(stmin+nstruct_-1)%map(3,tmp(stmin+nstruct_-1)%n_))
            end if
         end do
         call move_alloc(tmp,this%struct)
      end block compact_tree
      
      ! ! Fill out the node map
      ! node_map: block
      !    use amrex_amr_module, only: amrex_mfiter,amrex_box
      !    integer :: i,j,k
      !    integer, dimension(:), allocatable :: counter
      !    integer :: lvl
      !    type(amrex_mfiter) :: mfi
      !    type(amrex_box) :: bx
      !    allocate(counter(stmin:stmax)); counter=0
      !    ! Traverse levels ! Only work on finest level for now
      !    do lvl=this%amr%finest_level()
      !       ! Loop over tiles
      !       call this%amr%mfiter_build(lvl,mfi)
      !       do while (mfi%next())
      !          ! Get pointers to data
      !          pid=>this%id%mf(lvl)%dataptr(mfi)
      !          ! Perform local loop
      !          bx=mfi%tilebox()
      !          do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
      !             if (pid(i,j,k,1).gt.0) then
      !                counter(pid(i,j,k,1))=counter(pid(i,j,k,1))+1
      !                this%struct(pid(i,j,k,1))%map(:,counter(pid(i,j,k,1))))=[i,j,k]
      !             end if
      !          end do; end do; end do
      !    deallocate(counter)
      ! end block node_map
      
      ! ! Interprocessor treatment of our structures
      ! interproc_handling: block
      !    use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
      !    integer :: i,j,k,stop_global,stop_,counter,n,m,ierr,find_parent,find_parent_own
      !    ! Allocate to total number of structures
      !    allocate(parent    (this%nstruct)); parent    =0
      !    allocate(parent_all(this%nstruct)); parent_all=0
      !    allocate(parent_own(this%nstruct)); parent_own=0
      !    ! Fill global lineage with selves
      !    do n=1,this%nstruct
      !       parent(n)=n
      !    end do
      !    ! Synchronize id array
      !    call sync_lvl(this%id,this%amr%finest_level())
      !    ! Handle imin_ border
      !    if (this%pg%imin_.ne.this%pg%imin) then ! ??????????????
      !       ! Traverse levels ! Only work on finest level for now
      !       do lvl=this%amr%finest_level()
      !          ! Loop over tiles
      !          call this%amr%mfiter_build(lvl,mfi)
      !          do while (mfi%next())
      !             ! Get pointers to data
      !             pid=>this%id%mf(lvl)%dataptr(mfi)
      !             pVF=>this%VF%mf(lvl)%dataptr(mfi)
      !             ! Perform local loop
      !             bx=mfi%tilebox()
      !             do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
      !                if (pid(bx%lo(1),j,k,1).gt.0.and(pid(bx%lo(1)-1,j,k,1).gt.0)) then
      !                   if (same_label(pVF(bx%lo(1),j,k),pVF(bx%lo(1)-1,j,k))) call union_parent(pid(bx%lo(1),j,k),pid(bx%lo(1)-1,j,k))
      !                end if
      !             end do; end do
      !          end if
      !       end do; end do
      !    end if
      !    ! Handle jmin_ border
      !    if (this%pg%jmin_.ne.this%pg%jmin) then ! ?????????????            
      !       ! Traverse levels ! Only work on finest level for now
      !       do lvl=this%amr%finest_level()
      !          ! Loop over tiles
      !          call this%amr%mfiter_build(lvl,mfi)
      !          do while (mfi%next())
      !             ! Get pointers to data
      !             pid=>this%id%mf(lvl)%dataptr(mfi)
      !             pVF=>this%VF%mf(lvl)%dataptr(mfi)
      !             ! Perform local loop
      !             bx=mfi%tilebox()
      !             do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
      !                if (pid(i,bx%lo(2),k,1).gt.0.and(pid(i,bx%lo(2)-1,k,1).gt.0)) then
      !                   if (same_label(pVF(i,bx%lo(2),k),pVF(i,bx%lo(2)-1,k))) call union_parent(pid(i,bx%lo(2),k),pid(i,bx%lo(2)-1,k))
      !                end if
      !             end do; end do
      !          end if
      !       end do; end do
      !    end if
      !    ! Handle kmin_ border
      !    if (this%pg%kmin_.ne.this%pg%kmin) then ! ?????????????
      !       ! Traverse levels ! Only work on finest level for now
      !       do lvl=this%amr%finest_level()
      !          ! Loop over tiles
      !          call this%amr%mfiter_build(lvl,mfi)
      !          do while (mfi%next())
      !             ! Get pointers to data
      !             pid=>this%id%mf(lvl)%dataptr(mfi)
      !             pVF=>this%VF%mf(lvl)%dataptr(mfi)
      !             ! Perform local loop
      !             bx=mfi%tilebox()
      !             do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
      !                if (pid(i,j,bx%lo(3),1).gt.0.and(pid(i,j,bx%lo(3)-1,1).gt.0)) then
      !                   if (same_label(pVF(i,j,bx%lo(3)),pVF(i,j,bx%lo(3)-1))) call union_parent(pid(i,j,bx%lo(3)),pid(i,j,bx%lo(3)-1))
      !                end if
      !             end do; end do
      !          end if
      !       end do; end do
      !    end if
      !    ! Initialize global stop criterion and counter
      !    stop_global=1
      !    counter=0
      !    ! Resolve lineage
      !    do while (stop_global.ne.0)
      !       ! Initialize local stop flag
      !       stop_=0
      !       ! Remember own parents
      !       parent_own=parent
      !       ! Set self-parents to huge(1)
      !       do n=1,this%nstruct
      !          if (parent(n).eq.n) parent(n)=huge(1)
      !       end do
      !       ! Take global min
      !       call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%pg%comm,ierr)
      !       ! Set self-parents back to selves
      !       do n=1,this%nstruct
      !          if (parent_all(n).eq.huge(1)) parent_all(n)=n
      !       end do
      !       ! Flatten trees
      !       do n=1,this%nstruct
      !          parent_all(n)=find_all(n)
      !          parent_own(n)=find_own(n)
      !       end do
      !       ! Start with final parent array being equal to parent_all
      !       parent=parent_all
      !       ! Increment counter
      !       counter=counter+1
      !       ! Reconcile conflicts between parent_all and parent_own
      !       do n=1,this%nstruct
      !          if (parent_own(n).ne.n) then
      !             find_parent_own=rootify_parent(parent_own(n))
      !             find_parent    =rootify_parent(parent(n))
      !             if (find_parent_own.ne.find_parent) then
      !                call union_parent(find_parent,find_parent_own)
      !                stop_=1
      !             end if
      !          end if
      !       end do
      !       ! Check if we did some changes
      !       call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%pg%comm,ierr)
      !    end do
      !    ! Update this%struct%parent by pointing all parents to root and update id
      !    do n=stmin,stmax
      !       this%struct(n)%parent=rootify_parent(parent(n))
      !       do m=1,this%struct(n)%n_
      !          this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=this%struct(n)%parent
      !       end do
      !    end do
      ! end block interproc_handling
      
      ! ! Update periodicity array across processors
      ! periodicity_update: block
      !    use mpi_f08, only: MPI_ALLGATHER,MPI_MAX,MPI_INTEGER
      !    integer, dimension(:,:), allocatable :: ownper,allper
      !    integer :: n,m,ierr
      !    ! Allocate local and global perodicity arrays
      !    allocate(ownper(1:3,this%nstruct)); ownper=0
      !    allocate(allper(1:3,this%nstruct)); allper=0
      !    ! Fill ownper array
      !    do n=stmin,stmax
      !       ownper(:,n)=this%struct(n)%per
      !    end do
      !    ! Communicate per
      !    call MPI_ALLREDUCE(ownper,allper,3*this%nstruct,MPI_INTEGER,MPI_MAX,this%pg%comm,ierr)
      !    ! Update parent per
      !    do n=1,this%nstruct
      !       allper(:,parent(n))=max(allper(:,parent(n)),allper(:,n))
      !    end do
      !    ! Update idp array
      !    do n=stmin,stmax
      !       do m=1,this%struct(n)%n_
      !          idp(:,this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=allper(:,this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m)))
      !       end do
      !    end do
      !    ! Clean up
      !    deallocate(ownper,allper)
      ! end block periodicity_update
      
      ! ! One more pass for domain boundaries
      ! boundary_handling: block
      !    use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
      !    integer :: i,j,k,stop_global,stop_,counter,n,m,ierr,find_parent,find_parent_own
      !    ! Handle imin border
      !    if (this%pg%imin_.eq.this%pg%imin) then
      !       do k=this%pg%kmin_,this%pg%kmax_; do j=this%pg%jmin_,this%pg%jmax_
      !          if (this%id(this%pg%imin_,j,k).gt.0.and.this%id(this%pg%imin_-1,j,k).gt.0) then
      !             if (same_label(this%pg%imin_,j,k,this%pg%imin_-1,j,k)) call union_parent(this%id(this%pg%imin_,j,k),this%id(this%pg%imin_-1,j,k))
      !          end if
      !       end do; end do
      !    end if
      !    ! Handle jmin border
      !    if (this%pg%jmin_.eq.this%pg%jmin) then
      !       do k=this%pg%kmin_,this%pg%kmax_; do i=this%pg%imin_,this%pg%imax_
      !          if (this%id(i,this%pg%jmin_,k).gt.0.and.this%id(i,this%pg%jmin_-1,k).gt.0) then
      !             if (same_label(i,this%pg%jmin_,k,i,this%pg%jmin_-1,k)) call union_parent(this%id(i,this%pg%jmin_,k),this%id(i,this%pg%jmin_-1,k))
      !          end if
      !       end do; end do
      !    end if
      !    ! Handle kmin border
      !    if (this%pg%kmin_.eq.this%pg%kmin) then
      !       do j=this%pg%jmin_,this%pg%jmax_; do i=this%pg%imin_,this%pg%imax_
      !          if (this%id(i,j,this%pg%kmin_).gt.0.and.this%id(i,j,this%pg%kmin_-1).gt.0) then
      !             if (same_label(i,j,this%pg%kmin_,i,j,this%pg%kmin_-1)) call union_parent(this%id(i,j,this%pg%kmin_),this%id(i,j,this%pg%kmin_-1))
      !          end if
      !       end do; end do
      !    end if
      !    ! Initialize global stop criterion and counter
      !    stop_global=1
      !    counter=0
      !    ! Resolve lineage
      !    do while (stop_global.ne.0)
      !       ! Initialize local stop flag
      !       stop_=0
      !       ! Remember own parents
      !       parent_own=parent
      !       ! Set self-parents to huge(1)
      !       do n=1,this%nstruct
      !          if (parent(n).eq.n) parent(n)=huge(1)
      !       end do
      !       ! Take global min
      !       call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%pg%comm,ierr)
      !       ! Set self-parents back to selves
      !       do n=1,this%nstruct
      !          if (parent_all(n).eq.huge(1)) parent_all(n)=n
      !       end do
      !       ! Flatten trees
      !       do n=1,this%nstruct
      !          parent_all(n)=find_all_2(n,n)
      !          parent_own(n)=find_own(n)
      !       end do
      !       ! Start with final parent array being equal to parent_all
      !       parent=parent_all
      !       ! Increment counter
      !       counter=counter+1
      !       ! Reconcile conflicts between parent_all and parent_own
      !       do n=1,this%nstruct
      !          if (parent_own(n).ne.n) then
      !             find_parent_own=rootify_parent(parent_own(n))
      !             find_parent    =rootify_parent(parent(n))
      !             if (find_parent_own.ne.find_parent) then
      !                call union_parent(find_parent,find_parent_own)
      !                stop_=1
      !             end if
      !          end if
      !       end do
      !       ! Check if we did some changes
      !       call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%pg%comm,ierr)
      !    end do
      !    ! Update this%struct%parent and point all parents to root and update id
      !    do n=stmin,stmax
      !       this%struct(n)%parent=rootify_parent(parent(n))
      !       do m=1,this%struct(n)%n_
      !          this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=this%struct(n)%parent
      !       end do
      !    end do
      !    ! Update ghost cells
      !    call this%pg%sync(this%id)
      !    ! Clean up parent info
      !    deallocate(parent,parent_all,parent_own)
      ! end block boundary_handling
      
      ! ! Now we need to compact the data based on id only
      ! compact_struct: block
      !    use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
      !    integer :: i,j,k,n,nn,ierr,count
      !    integer, dimension(:), allocatable :: my_idmap,idmap,counter
      !    type(struct_type), dimension(:), allocatable :: tmp
      !    ! Prepare global id map
      !    allocate(my_idmap(1:this%nstruct)); my_idmap=0
      !    allocate(   idmap(1:this%nstruct));    idmap=0
      !    ! Traverse id array and tag used id values
      !    do k=this%pg%kmin_,this%pg%kmax_; do j=this%pg%jmin_,this%pg%jmax_; do i=this%pg%imin_,this%pg%imax_
      !       if (this%id(i,j,k).gt.0) my_idmap(this%id(i,j,k))=1
      !    end do; end do; end do
      !    call MPI_ALLREDUCE(my_idmap,idmap,this%nstruct,MPI_INTEGER,MPI_MAX,this%pg%comm,ierr)
      !    deallocate(my_idmap)
      !    ! Count number of used structures and create the map
      !    this%nstruct=sum(idmap)
      !    count=0
      !    do n=1,size(idmap,dim=1)
      !       if (idmap(n).gt.0) then
      !          count=count+1
      !          idmap(n)=count
      !       end if
      !    end do
      !    ! Rename all structures
      !    do k=this%pg%kmin_,this%pg%kmax_; do j=this%pg%jmin_,this%pg%jmax_; do i=this%pg%imin_,this%pg%imax_
      !       if (this%id(i,j,k).gt.0) this%id(i,j,k)=idmap(this%id(i,j,k))
      !    end do; end do; end do
      !    call this%pg%sync(this%id)
      !    ! Allocate temporary storage for structure
      !    allocate(tmp(this%nstruct))
      !    allocate(counter(this%nstruct)); counter=0
      !    do k=this%pg%kmin_,this%pg%kmax_; do j=this%pg%jmin_,this%pg%jmax_; do i=this%pg%imin_,this%pg%imax_
      !       if (this%id(i,j,k).gt.0) counter(this%id(i,j,k))=counter(this%id(i,j,k))+1
      !    end do; end do; end do
      !    do n=1,this%nstruct
      !       tmp(n)%parent=n
      !       tmp(n)%per=0
      !       tmp(n)%n_=counter(n)
      !       allocate(tmp(n)%map(1:3,1:tmp(n)%n_))
      !    end do
      !    ! Transfer periodicity info
      !    do n=stmin,stmax
      !       if (idmap(n).gt.0) then
      !          tmp(idmap(n))%per=this%struct(n)%per
      !       end if
      !    end do
      !    deallocate(idmap)
      !    ! Store the map
      !    counter=0
      !    do k=this%pg%kmin_,this%pg%kmax_; do j=this%pg%jmin_,this%pg%jmax_; do i=this%pg%imin_,this%pg%imax_
      !       if (this%id(i,j,k).gt.0) then
      !          counter(this%id(i,j,k))=counter(this%id(i,j,k))+1
      !          tmp(this%id(i,j,k))%map(:,counter(this%id(i,j,k)))=[i,j,k]
      !       end if
      !    end do; end do; end do
      !    deallocate(counter)
      !    ! Transfer allocation
      !    call move_alloc(tmp,this%struct)
      !    ! Final pass to fix periodicity info
      !    do n=1,this%nstruct
      !       do nn=1,this%struct(n)%n_
      !          i=this%struct(n)%map(1,nn)
      !          j=this%struct(n)%map(2,nn)
      !          k=this%struct(n)%map(3,nn)
      !          this%struct(n)%per(1)=max(this%struct(n)%per(1),idp(1,i,j,k))
      !          this%struct(n)%per(2)=max(this%struct(n)%per(2),idp(2,i,j,k))
      !          this%struct(n)%per(3)=max(this%struct(n)%per(3),idp(3,i,j,k))
      !       end do
      !    end do
      !    deallocate(idp)
      ! end block compact_struct
      
      ! Extra QOL step to ensure that id=1 is always the largest structure in terms of number of cells
      rename_largest_structure: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER,MPI_IN_PLACE
         integer :: ierr,bigid,i,j,k
         integer, dimension(:), allocatable :: ncells
         type(struct_type) :: tmp
         ! Skip if no structure was found
         if (this%nstruct.eq.0) exit rename_largest_structure
         ! Loop over all structures and count total number of cells to find ID of largest structure
         allocate(ncells(1:this%nstruct)); ncells=this%struct(:)%n_
         call MPI_ALLREDUCE(MPI_IN_PLACE,ncells,this%nstruct,MPI_INTEGER,MPI_SUM,this%pg%comm,ierr)
         bigid=maxloc(ncells,1)
         deallocate(ncells)
         ! Swap structures
         tmp=this%struct(1); this%struct(1)=this%struct(bigid); this%struct(bigid)=tmp
         do k=this%pg%kmino_,this%pg%kmaxo_; do j=this%pg%jmino_,this%pg%jmaxo_; do i=this%pg%imino_,this%pg%imaxo_
            if (this%id(i,j,k).eq.1) then; this%id(i,j,k)=bigid; else if (this%id(i,j,k).eq.bigid) then; this%id(i,j,k)=1; end if
         end do; end do; end do
      end block rename_largest_structure
      
      
   contains
      
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
            size_new=int(real(size_now,WP)*coeff_up)
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
      class(cclabel), intent(inout) :: this
      integer :: n
      ! Loop over all structures and deallocate maps
      if (allocated(this%struct)) then
         do n=1,size(this%struct,dim=1)
            if (allocated(this%struct(n)%map)) deallocate(this%struct(n)%map)
         end do
         ! Deallocate structure array
         deallocate(this%struct)
      end if
      ! Zero structures
      this%nstruct=0
      ! Reset id to zero
      this%id=0
   end subroutine empty
   
   
   !> Finalize CCL object
   subroutine finalize(this)
      implicit none
      class(cclabel), intent(inout) :: this
      call this%empty()
      if (allocated(this%id)) deallocate(this%id)
      nullify(this%pg)
      this%name='UNNAMED_CCL'
   end subroutine finalize
   
   
end module cclabel_class
