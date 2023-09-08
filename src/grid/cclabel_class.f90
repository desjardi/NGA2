!> Connected component labeling class: identifies Lagrangian objects from a Eulerian logical field
!> and provides unstructured mapping to traverse these objects
module cclabel_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   implicit none
   private
   

   ! Expose type/constructor/methods
   public :: cclabel
   

   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%
   

   !> Structure object
   type :: struct_type
      integer :: parent                                   !< ID of parent struct
      integer :: n_                                       !< Number of local cells contained in struct
      integer, dimension(:,:), allocatable :: map         !< List of cells contained in struct, dimension(1:3,1:n_)
      integer, dimension(3) :: per                        !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct_type
   

   !> cclabel object definition
   type :: cclabel
      ! This is our config
      class(config), pointer :: cfg
      ! This is the name of the CCL
      character(len=str_medium) :: name='UNNAMED_CCL'
      ! Eulerian field of tagged cells - input set by user
      logical, dimension(:,:,:), allocatable :: tagged
      ! ID of the structure that contains each cell
      integer, dimension(:,:,:), allocatable :: id
      ! Array of structures
      integer :: nstruct_,nstruct
      type(struct_type), dimension(:), allocatable :: struct
      integer :: stmin,stmax
   contains
      procedure :: initialize
      procedure :: build
      procedure :: empty
   end type cclabel
   

contains
   
   
   !> Initialization for cclabel class
   subroutine initialize(this,cfg,name)
      implicit none
      class(cclabel) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Point to cfg object
      this%cfg=>cfg
      ! Allocate and initialize logical array
      allocate(this%tagged(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%tagged=.false.
      ! Allocate and initialize ID array
      allocate(this%id(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%id=0
      ! Zero structures
      this%nstruct_=0; this%nstruct=0
   end subroutine initialize
   

   !> Build structure using the user-set tagged field
   subroutine build(this)
      implicit none
      class(cclabel), intent(inout) :: this
      integer, dimension(:,:,:,:), allocatable :: idp          !< Periodicity treatment
      integer, dimension(:), allocatable :: parent             !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_all         !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_own         !< Resolving structure id across procs
      
      ! Allocate periodicity work array
      allocate(idp(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); idp=0
      
      ! Perform a first pass to build proc-local structures and corresponding tree
      first_pass: block
         integer :: i,j,k,ii,jj,kk,dim,dir
         integer, dimension(3) :: pos
         ! Perform local loop
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            ! Find next tagged cell
            if (this%tagged(i,j,k)) then
               ! Loop through one-sided neighbors
               do dim=1,3
                  pos=0; pos(dim)=-1
                  ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                  ! Check if neighbor is labeled
                  if (this%id(ii,jj,kk).gt.0) then
                     ! Neighbor is labeled, but are we?
                     if (this%id(i,j,k).ne.0) then
                        ! We already have a label, perform a union of both labels
                        this%id(i,j,k)=union(this%id(i,j,k),this%id(ii,jj,kk))
                     else
                        ! We don't have a label, take the neighbor's label
                        this%id(i,j,k)=this%id(ii,jj,kk)
                     end if
                  end if
               end do
               ! If no neighbor was labeled, we need a new structure
               if (this%id(i,j,k).eq.0) this%id(i,j,k)=add()
               ! Identify periodicity cases
               if (this%cfg%xper.and.i.eq.this%cfg%imax) this%struct(this%id(i,j,k))%per(1)=1
               if (this%cfg%yper.and.j.eq.this%cfg%jmax) this%struct(this%id(i,j,k))%per(2)=1
               if (this%cfg%zper.and.k.eq.this%cfg%kmax) this%struct(this%id(i,j,k))%per(3)=1
               idp(:,i,j,k)=this%struct(this%id(i,j,k))%per
            end if
         end do; end do; end do
      end block first_pass

      ! Now collapse the tree, count the cells and resolve periodicity in each structure
      collapse_tree: block
         integer :: i,j,k
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%id(i,j,k).gt.0) then
               this%id(i,j,k)=rootify(this%id(i,j,k))
               this%struct(this%id(i,j,k))%n_=this%struct(this%id(i,j,k))%n_+1
               idp(1,i,j,k)=max(idp(1,i,j,k),this%struct(this%id(i,j,k))%per(1))
               idp(2,i,j,k)=max(idp(2,i,j,k),this%struct(this%id(i,j,k))%per(2))
               idp(3,i,j,k)=max(idp(3,i,j,k),this%struct(this%id(i,j,k))%per(3))
               this%struct(this%id(i,j,k))%per=idp(:,i,j,k)
            end if
         end do; end do; end do
      end block collapse_tree
      
      ! Compact structure array
      compact_tree: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
         integer :: i,j,k,n,ierr,max_nstruct
         integer, dimension(:), allocatable :: my_nstruct,all_nstruct,idmap
         type(struct_type), dimension(:), allocatable :: tmp
         ! Count exact number of local structures
         this%nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%n_.gt.0) this%nstruct_=this%nstruct_+1
         end do
         ! Gather this info to ensure unique index
         allocate( my_nstruct(0:this%cfg%nproc-1)); my_nstruct=0; my_nstruct(this%cfg%rank)=this%nstruct_
         allocate(all_nstruct(0:this%cfg%nproc-1)); call MPI_ALLREDUCE(my_nstruct,all_nstruct,this%cfg%nproc,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
         this%stmin=1
         if (this%cfg%rank.gt.0) this%stmin=this%stmin+sum(all_nstruct(0:this%cfg%rank-1))
         this%nstruct=sum(all_nstruct)
         deallocate(my_nstruct,all_nstruct)
         this%stmax=this%stmin+this%nstruct_-1
         ! Generate an index map
         allocate(idmap(1:size(this%struct,dim=1))); idmap=0
         this%nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%n_.gt.0) then
               this%nstruct_=this%nstruct_+1
               idmap(n)=this%stmin+this%nstruct_-1
            end if
         end do
         ! Update id array to new index
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%id(i,j,k).gt.0) this%id(i,j,k)=idmap(this%id(i,j,k))
         end do; end do; end do
         deallocate(idmap)
         ! Finish compacting and renumbering
         allocate(tmp(this%stmin:this%stmax))
         this%nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%n_.gt.0) then
               this%nstruct_=this%nstruct_+1
               tmp(this%stmin+this%nstruct_-1)=this%struct(n)
               allocate(tmp(this%stmin+this%nstruct_-1)%map(3,tmp(this%stmin+this%nstruct_-1)%n_))
            end if
         end do
         call move_alloc(tmp,this%struct)
      end block compact_tree
      
      ! Fill out the node map
      node_map: block
         integer :: i,j,k
         integer, dimension(:), allocatable :: counter
         allocate(counter(this%stmin:this%stmax)); counter=0
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%id(i,j,k).gt.0) then
               counter(this%id(i,j,k))=counter(this%id(i,j,k))+1
               this%struct(this%id(i,j,k))%map(:,counter(this%id(i,j,k)))=[i,j,k]
            end if
         end do; end do; end do
         deallocate(counter)
      end block node_map
      
      ! Interprocessor treatment of our structures
      interproc_handling: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
         integer :: i,j,k,stop_global,stop_,counter,n,m,ierr,find_parent,find_parent_own
         ! Allocate to total number of structures
         allocate(parent    (this%nstruct)); parent    =0
         allocate(parent_all(this%nstruct)); parent_all=0
         allocate(parent_own(this%nstruct)); parent_own=0
         ! Fill parent with selves
         do i=1,this%nstruct
            parent(i)=i
         end do
         ! Synchronize id array
         call this%cfg%sync(this%id)
         ! Handle imin_ border
         if (this%cfg%imin_.ne.this%cfg%imin) then
            do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_
               ! If both border cell and neighboring ghost cell are filled
               if (this%id(this%cfg%imin_,j,k).gt.0.and.this%id(this%cfg%imin_-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (parent(this%id(this%cfg%imin_,j,k)).eq.this%id(this%cfg%imin_,j,k)) then
                     ! Directly set the lineage
                     parent(this%id(this%cfg%imin_,j,k))=this%id(this%cfg%imin_-1,j,k)
                  else if (find(this%id(this%cfg%imin_,j,k)).ne.find(this%id(this%cfg%imin_-1,j,k))) then
                     if (this%id(this%cfg%imin_-1,j,k).gt.parent(this%id(this%cfg%imin_,j,k))) then
                        parent(find(this%id(this%cfg%imin_-1,j,k)))=find(parent(this%id(this%cfg%imin_,j,k)))
                     else
                        parent(find(parent(this%id(this%cfg%imin_,j,k))))=find(this%id(this%cfg%imin_-1,j,k))
                     end if
                  end if
               end if
            end do; end do
         end if
         ! Handle jmin_ border
         if (this%cfg%jmin_.ne.this%cfg%jmin) then
            do k=this%cfg%kmin_,this%cfg%kmax_; do i=this%cfg%imin_,this%cfg%imax_
               ! If both border cell and neighboring ghost cell are filled
               if (this%id(i,this%cfg%jmin_,k).gt.0.and.this%id(i,this%cfg%jmin_-1,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (parent(this%id(i,this%cfg%jmin_,k)).eq.this%id(i,this%cfg%jmin_,k)) then
                     ! Directly set the lineage
                     parent(this%id(i,this%cfg%jmin_,k))=this%id(i,this%cfg%jmin_-1,k)
                  elseif (find(this%id(i,this%cfg%jmin_,k)).ne.find(this%id(i,this%cfg%jmin_-1,k))) then
                     if (this%id(i,this%cfg%jmin_-1,k).gt.parent(this%id(i,this%cfg%jmin_,k))) then
                        parent(find(this%id(i,this%cfg%jmin_-1,k)))=find(parent(this%id(i,this%cfg%jmin_,k)))
                     else
                        parent(find(parent(this%id(i,this%cfg%jmin_,k))))=find(this%id(i,this%cfg%jmin_-1,k))
                     end if
                  end if
               end if
            end do; end do
         end if
         ! Handle jmin_ border
         if (this%cfg%kmin_.ne.this%cfg%kmin) then
            do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
               ! If both border cell and neighboring ghost cell are filled
               if (this%id(i,j,this%cfg%kmin_).gt.0.and.this%id(i,j,this%cfg%kmin_-1).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (parent(this%id(i,j,this%cfg%kmin_)).eq.this%id(i,j,this%cfg%kmin_)) then
                     ! Directly set the lineage
                     parent(this%id(i,j,this%cfg%kmin_))=this%id(i,j,this%cfg%kmin_-1)
                  elseif (find(this%id(i,j,this%cfg%kmin_)).ne.find(this%id(i,j,this%cfg%kmin_-1))) then
                     if (this%id(i,j,this%cfg%kmin_-1).gt.parent(this%id(i,j,this%cfg%kmin_))) then
                        parent(find(this%id(i,j,this%cfg%kmin_-1)))=find(parent(this%id(i,j,this%cfg%kmin_)))
                     else
                        parent(find(parent(this%id(i,j,this%cfg%kmin_))))=find(this%id(i,j,this%cfg%kmin_-1))
                     end if
                  end if
               end if
            end do; end do
         end if
         ! Initialize global stop criterion and counter
         stop_global=1
         counter=0
         ! Resolve lineage
         do while (stop_global.ne.0)
            ! Initialize local stop flag
            stop_=0
            ! Remember own parents
            parent_own=parent
            ! Set self-parents to huge(1) and take global min
            do n=1,this%nstruct
               if (parent(n).eq.n) parent(n)=huge(1)
            end do
            call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
            ! Set self-parents back to selves
            do n=1,this%nstruct
               if (parent_all(n).eq.huge(1)) parent_all(n)=n
            end do
            ! Flatten trees - is this necessary?
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
                  find_parent_own=find(parent_own(n))
                  find_parent    =find(parent(n))
                  if (find_parent_own.gt.find_parent) then
                     parent(find(parent_own(n)))=find(parent(n))
                     stop_=1
                  else if (find_parent.gt.find_parent_own) then
                     parent(find(parent(n)))=find(parent_own(n))
                     stop_=1
                  end if
               end if
            end do
            ! Check if we did some changes
            call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         end do
         ! Update this%struct%parent and point all parents to root and update id
         do n=this%stmin,this%stmax
            this%struct(n)%parent=find(parent(n))
            do m=1,this%struct(n)%n_
               this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=this%struct(n)%parent
            end do
         end do
      end block interproc_handling
      
      ! Update periodicity array across processors
      periodicity_update: block
         use mpi_f08, only: MPI_ALLGATHER,MPI_MAX,MPI_INTEGER
         integer, dimension(:,:), allocatable :: ownper,allper
         integer :: n,m,ierr
         ! Allocate local and global perodicity arrays
         allocate(ownper(1:3,this%nstruct)); ownper=0
         allocate(allper(1:3,this%nstruct)); allper=0
         ! Fill ownper array
         do n=this%stmin,this%stmax
            ownper(:,n)=this%struct(n)%per
         end do
         ! Communicate per
         call MPI_ALLREDUCE(ownper,allper,3*this%nstruct,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         ! Update parent per
         do n=1,this%nstruct
            allper(:,parent(n))=max(allper(:,parent(n)),allper(:,n))
         end do
         ! Update idp array
         do n=this%stmin,this%stmax
            do m=1,this%struct(n)%n_
               idp(:,this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=allper(:,this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m)))
            end do
         end do
      end block periodicity_update

      ! One more pass for domain boundaries
      !boundary_handling: block
      !   use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
      !   integer :: i,j,k,stop_global,stop_,counter,n,m,ierr,find_parent,find_parent_own
      !   ! Handle imin border
      !   if (this%cfg%imin_.eq.this%cfg%imin) then
      !      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_
      !         ! If both border cell and neighboring ghost cell are filled
      !         if (this%id(this%cfg%imin,j,k).gt.0.and.this%id(this%cfg%imin-1,j,k).gt.0) then
      !            ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
      !            if (parent(this%id(this%cfg%imin,j,k)).eq.this%id(this%cfg%imin,j,k)) then
      !               parent(find(this%id(this%cfg%imin,j,k)))=find(this%id(this%cfg%imin-1,j,k))
      !            else if (find(this%id(this%cfg%imin,j,k)).ne.find(this%id(this%cfg%imin-1,j,k))) then
      !               if (this%id(this%cfg%imin-1,j,k).gt.parent(this%id(this%cfg%imin,j,k))) then
      !                  parent(find(this%id(this%cfg%imin-1,j,k)))=find(parent(this%id(this%cfg%imin,j,k)))
      !               else
      !                  parent(find(parent(this%id(this%cfg%imin,j,k))))=find(this%id(this%cfg%imin-1,j,k))
      !               end if
      !            end if
      !         end if
      !      end do; end do
      !   end if
      !   ! Handle jmin border
      !   if (this%cfg%jmin_.eq.this%cfg%jmin) then
      !      do k=this%cfg%kmin_,this%cfg%kmax_; do i=this%cfg%imin_,this%cfg%imax_
      !         ! If both border cell and neighboring ghost cell are filled
      !         if (this%id(i,this%cfg%jmin,k).gt.0.and.this%id(i,this%cfg%jmin-1,k).gt.0) then
      !            ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
      !            if (parent(this%id(i,this%cfg%jmin,k)).eq.this%id(i,this%cfg%jmin,k)) then
      !               parent(find(this%id(i,this%cfg%jmin,k)))=find(this%id(i,this%cfg%jmin-1,k))
      !            elseif (find(this%id(i,this%cfg%jmin,k)).ne.find(this%id(i,this%cfg%jmin-1,k))) then
      !               if (this%id(i,this%cfg%jmin-1,k).gt.parent(this%id(i,this%cfg%jmin,k))) then
      !                  parent(find(this%id(i,this%cfg%jmin-1,k)))=find(parent(this%id(i,this%cfg%jmin,k)))
      !               else
      !                  parent(find(parent(this%id(i,this%cfg%jmin,k))))=find(this%id(i,this%cfg%jmin-1,k))
      !               end if
      !            end if
      !         end if
      !      end do; end do
      !   end if
      !   ! Handle kmin border
      !   if (this%cfg%kmin_.eq.this%cfg%kmin) then
      !      do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
      !         ! If both border cell and neighboring ghost cell are filled
      !         if (this%id(i,j,this%cfg%kmin).gt.0.and.this%id(i,j,this%cfg%kmin-1).gt.0) then
      !            ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
      !            if (parent(this%id(i,j,this%cfg%kmin)).eq.this%id(i,j,this%cfg%kmin)) then
      !               parent(find(this%id(i,j,this%cfg%kmin)))=find(this%id(i,j,this%cfg%kmin-1))
      !            elseif (find(this%id(i,j,this%cfg%kmin)).ne.find(this%id(i,j,this%cfg%kmin-1))) then
      !               if (this%id(i,j,this%cfg%kmin-1).gt.parent(this%id(i,j,this%cfg%kmin))) then
      !                  parent(find(this%id(i,j,this%cfg%kmin-1)))=find(parent(this%id(i,j,this%cfg%kmin)))
      !               else
      !                  parent(find(parent(this%id(i,j,this%cfg%kmin))))=find(this%id(i,j,this%cfg%kmin-1))
      !               end if
      !            end if
      !         end if
      !      end do; end do
      !   end if
      !   ! Initialize global stop criterion and counter
      !   stop_global=1
      !   counter=0
      !   ! Resolve lineage
      !   do while (stop_global.ne.0)
      !      ! Initialize local stop flag
      !      stop_=0
      !      ! Remember own parents
      !      parent_own=parent
      !      ! Set self-parents to huge(1) and take global min
      !      do n=1,this%nstruct
      !         if (parent(n).eq.n) parent(n)=huge(1)
      !      end do
      !      call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%cfg%comm,ierr)
      !      ! Set self-parents back to selves
      !      do n=1,this%nstruct
      !         if (parent_all(n).eq.huge(1)) parent_all(n)=n
      !      end do
      !      ! Flatten trees - is this necessary?
      !      do n=1,this%nstruct
      !         parent_all(n)=find_all_2(n,n)
      !         parent_own(n)=find_own(n)
      !      end do
      !      ! Start with final parent array being equal to parent_all
      !      parent=parent_all
      !      ! Increment counter
      !      counter=counter+1
      !      ! Reconcile conflicts between parent_all and parent_own
      !      do n=1,this%nstruct
      !         if (parent_own(n).ne.n) then
      !            find_parent_own=find(parent_own(n))
      !            find_parent    =find(parent(n))
      !            if (find_parent_own.gt.find_parent) then
      !               parent(find(parent_own(n)))=find(parent(n))
      !               stop_=1
      !            else if (find_parent.gt.find_parent_own) then
      !               parent(find(parent(n)))=find(parent_own(n))
      !               stop_=1
      !            end if
      !         end if
      !      end do
      !      ! Check if we did some changes
      !      call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      !   end do
      !   ! Update this%struct%parent and point all parents to root and update id
      !   do n=this%stmin,this%stmax
      !      this%struct(n)%parent=find(parent(n))
      !      do m=1,this%struct(n)%n_
      !         this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=this%struct(n)%parent
      !      end do
      !   end do
      !   ! Update ghost cells
      !   call this%cfg%sync(this%id)
      !end block boundary_handling
      
      ! Now we need to compact the data based on id only
      !compact_struct: block
      !   use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
      !   integer :: i,j,k,n,ierr,count
      !   integer, dimension(:), allocatable :: my_idmap,idmap
      !   type(struct_type), dimension(:), allocatable :: tmp
      !   ! Prepare global id map
      !   allocate(my_idmap(1:this%nstruct)); my_idmap=0
      !   allocate(   idmap(1:this%nstruct));    idmap=0
      !   ! Traverse id array and tag used id values
      !   do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
      !      if (this%id(i,j,k).gt.0) my_idmap(this%id(i,j,k))=1
      !   end do; end do; end do
      !   call MPI_ALLREDUCE(my_idmap,idmap,this%nstruct,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      !   deallocate(my_idmap)
      !   ! Count number of used structures and create the map
      !   this%nstruct=sum(idmap)
      !   count=0
      !   do n=1,size(idmap,dim=1)
      !      if (idmap(n).gt.0) then
      !         count=count+1
      !         idmap(n)=count
      !      end if
      !   end do
      !   ! Rename all structures
      !   do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
      !      if (this%id(i,j,k).gt.0) this%id(i,j,k)=idmap(this%id(i,j,k))
      !   end do; end do; end do
      !   call this%cfg%sync(this%id)
      !   deallocate(idmap)
      !   ! Allocate temporary storage for structure
      !   allocate(tmp(this%nstruct))
      !   allocate(idmap(this%nstruct)); idmap=0
      !   do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
      !      if (this%id(i,j,k).gt.0) idmap(this%id(i,j,k))=idmap(this%id(i,j,k))+1
      !   end do; end do; end do
      !   print*,idmap
      !   stop
      !   do n=1,this%nstruct
      !      tmp(n)%parent=n
      !      tmp(n)%n_=idmap(n)
      !      allocate(tmp(n)%map(1:3,1:tmp(n)%n_))
      !      tmp(n)%per=0
      !   end do
      !   ! Store the map
      !   idmap=0
      !   do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
      !      if (this%id(i,j,k).gt.0) idmap(this%id(i,j,k))=idmap(this%id(i,j,k))+1
      !      tmp(this%id(i,j,k))%map(:,idmap(this%id(i,j,k)))=[i,j,k]
      !   end do; end do; end do
      !   deallocate(idmap)
      !   ! Transfer allocation
      !   call move_alloc(tmp,this%struct)
      !end block compact_struct
      
   contains
      
      !> This recursive function points the parent to root and returns that root
      recursive function rootify(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.this%struct(y)%parent) then
            this%struct(y)%parent=rootify(this%struct(y)%parent)
            y=this%struct(y)%parent
         end if
      end function rootify
      
      !> This function joins two branches at their roots (the root of y is chosen and returned)
      function union(x,y) result(z)
         implicit none
         integer, intent(in) :: x,y
         integer :: z
         z=rootify(y)
         this%struct(rootify(x))%parent=z
      end function union
      
      !> This function adds one new root while dynamically handling storage space
      function add() result(x)
         implicit none
         integer :: x
         integer :: size_now,size_new
         type(struct_type), dimension(:), allocatable :: tmp
         ! Check size of struct
         if (.not.allocated(this%struct)) then
            ! Allocate minimum storage
            allocate(this%struct(min_struct_size))
            ! Add first root
            this%nstruct_=1
            this%struct(this%nstruct_)%parent=1
            this%struct(this%nstruct_)%per=0
            this%struct(this%nstruct_)%n_=0
         else
            ! Check if there is enough room for storing a new structure
            size_now=size(this%struct,dim=1)
            if (this%nstruct_.eq.size_now) then
               size_new=int(real(size_now,WP)*coeff_up)
               allocate(tmp(size_new))
               tmp(1:this%nstruct_)=this%struct
               tmp(this%nstruct_+1:)%parent=0
               tmp(this%nstruct_+1:)%n_=0
               call move_alloc(tmp,this%struct)
            end if
            ! Add new root
            this%nstruct_=this%nstruct_+1
            this%struct(this%nstruct_)%parent=this%nstruct_
            this%struct(this%nstruct_)%per=0
            this%struct(this%nstruct_)%n_=0
         end if
         x=this%nstruct_
      end function add
      
      !> This recursive function points parent to root and returns that root
      recursive function find(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.parent(y)) then
            parent(y)=find(parent(y))
            y=parent(y)
         end if
      end function find

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
      recursive function find_all_2(x,starting_node) result(y)
         implicit none
         integer, intent(in) :: x,starting_node
         integer :: y
         y=x
         if (y.ne.parent_all(y)) then
            if (parent_all(y).eq.starting_node) then
               y=parent_all(y)
               return
            else
               parent_all(y)=find_all_2(parent_all(y),starting_node)
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
      do n=this%stmin,this%stmax
         if (allocated(this%struct(n)%map)) deallocate(this%struct(n)%map)
      end do
      ! Deallocate structure array
      deallocate(this%struct)
      ! Zero structures
      this%nstruct_=0; this%nstruct=0
   end subroutine empty
   
   
end module cclabel_class
