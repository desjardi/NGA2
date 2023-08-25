!> Connected component labeling class: identifies Lagrangian objects from a Eulerian logical field
!> and provides unstructured mapping to traverse these objects
module ccl_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ccl
   
   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%

   !> Structure object
   type :: struct_type
      integer :: parent                                   !< ID of parent struct
      integer :: rank    ! NEEDED?                        !< Upper bound on height of tree; only used for roots
      integer :: counter ! NEEDED?                        !< Counter for filling node array
      integer :: nnode                                    !< Number of cells contained in struct
      integer, dimension(:,:), allocatable :: node        !< List of cells contained in struct, dimension(1:nnode,3)
      integer, dimension(3) :: per                        !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct_type
   
   !> CCL object definition
   type :: ccl
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
   end type ccl
   
contains
   
   !> Initialization for CCL class
   subroutine initialize(this,cfg,name)
      implicit none
      class(ccl) :: this
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
   end subroutine

   !> Build structure using the user-set
   subroutine build(this)
      implicit none
      class(ccl), intent(inout) :: this
      integer, dimension(:,:,:,:), allocatable :: idp          !< Periodicity treatment
      integer, dimension(:), allocatable :: parent             !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_all         !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_own         !< Resolving structure id across procs

      ! Allocate periodicity work array
      allocate(idp(3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); idp=0
      
      ! Perform a first pass to build proc-local structures and corresponding tree
      first_pass: block
         integer :: i,j,k,dim,dir
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
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%id(i,j,k).gt.0) then
               this%id(i,j,k)=rootify(this%struct,this%id(i,j,k))
               this%struct(this%id(i,j,k))%nnode=this%struct(this%id(i,j,k))%nnode+1
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
         ! Count exact number of local structures
         this%nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%nnode.gt.0) this%nstruct_=this%nstruct_+1
         end do
         ! Gather this info to ensure unique index
         allocate(my_nstruct(0:this%cfg%nproc-1)); my_nstruct=0; my_nstruct(this%cfg%rank)=this%nstruct_
         allocate(all_nstruct(0:this%cfg%nproc-1)); call MPI_ALLREDUCE(my_nstruct,all_nstruct,this%cfg%nproc,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
         this%stmin=1
         if (this%cfg%rank.gt.0) this%stmin=this%stmin+sum(all_nstruct(0:this%cfg%rank-1))
         this%nstruct=sum(all_nstruct)
         deallocate(my_nstruct,all_nstruct)
         this%stmax=this%stmin+this%nstruct_-1
         ! Generate an index map
         allocate(idmap(1,size(this%struct,dim=1))); idmap=0
         this%nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%nnode.gt.0) then
               this%nstruct_=this%nstruct_+1
               idmap(n)=this%stmin+this%nstruct_-1
            end if
         end do
         ! Update id array to new index
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            if (this%id(i,j,k).gt.0) this%id(i,j,k)=idmap(this%id(i,j,k))
         end do; end do; end do
         deallocate(idmap)
         ! Finish compact and renumbering
         allocate(tmp(this%stmin:this%stmax))
         this%nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%nnode.gt.0) then
               this%nstruct_=this%nstruct_+1
               tmp(this%stmin+this%nstruct_-1)=this%struct(n)
               allocate(tmp(this%stmin+this%nstruct_-1)%node(3,tmp(this%stmin+this%nstruct_-1)%nnode))
            end if
         end do
         call move_alloc(tmp,this%struct)
      end block compact_tree
      
      ! Fill out the node map
      node_map: block
         integer :: i,j,k
         integer, dimension(:), allocatable :: counter
         allocate(counter(this%stmin:this%stmax)); counter=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%id(i,j,k).gt.0) then
                     counter(this%id(i,j,k))=counter(this%id(i,j,k))+1
                     this%struct(this%id(i,j,k))%node(:,counter(this%id(i,j,k)))=[i,j,k]
                  end if
               end do
            end do
         end do
         deallocate(counter)
      end block node_map
      
      ! Interprocessor treatment of our structures
      interproc_handling: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
         integer :: i,j,k,stop_global,stop_,counter,n,m,ierr
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
                        parent(rootify2(this%id(this%cfg%imin_-1,j,k)))=rootify2(parent(this%id(this%cfg%imin_,j,k)))
                     else
                        parent(rootify2(parent(this%id(this%cfg%imin_,j,k))))=rootify2(this%id(this%cfg%imin_-1,j,k))
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
                        parent(rootify2(this%id(i,this%cfg%jmin_-1,k)))=rootify2(parent(this%id(i,this%cfg%jmin_,k)))
                     else
                        parent(rootify2(parent(this%id(i,this%cfg%jmin_,k))))=rootify2(this%id(i,this%cfg%jmin_-1,k))
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
                        parent(rootify2(this%id(i,j,this%cfg%kmin_-1)))=rootify2(parent(this%id(i,j,this%cfg%kmin_)))
                     else
                        parent(rootify2(parent(this%id(i,j,this%cfg%kmin_))))=rootify2(this%id(i,j,this%cfg%kmin_-1))
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
                     parent(rootify2(parent_own(n)))=rootify2(parent(n))
                     stop_=1
                  else if (find_parent.gt.find_parent_own) then
                     parent(rootify2(parent(n)))=rootify2(parent_own(n))
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
            do m=1,this%struct(n)%nnode
               this%id(this%struct(n)%node(1,m),this%struct(n)%node(2,m),this%struct(n)%node(3,m))=this%struct(n)%parent
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
         call MPI_ALLREDUCE(ownper(1,:),allper(1,:),this%nstruct,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(ownper(2,:),allper(2,:),this%nstruct,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(ownper(3,:),allper(3,:),this%nstruct,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         ! Update parent per
         do n=1,this%nstruct
            allper(:,parent(n))=max(allper(:,parent(n)),allper(:,n))
         end do
         ! Update idp array
         do n=this%stmin,this%stmax
            do m=1,this%struct(n)%nnode
               this%idp(:,this%struct(n)%node(1,m),this%struct(n)%node(2,m),this%struct(n)%node(3,m))=allper(:,this%id(this%struct(n)%node(1,m),this%struct(n)%node(2,m),this%struct(n)%node(3,m)))
            end do
         end do
      end block periodicity_update

      ! Final pass for domain boundaries
      boundary_handling: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
         integer :: i,j,k,stop_global,stop_,counter,n,m,ierr
         ! Handle imin border
         if (this%cfg%imin_.eq.this%cfg%imin) then
            do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_
               ! If both border cell and neighboring ghost cell are filled
               if (this%id(this%cfg%imin,j,k).gt.0.and.this%id(this%cfg%imin-1,j,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (parent(this%id(this%cfg%imin,j,k)).eq.this%id(this%cfg%imin,j,k)) then
                     parent(rootify2(this%id(this%cfg%imin,j,k)))=rootify2(this%id(this%cfg%imin-1,j,k))
                  else if (find(this%id(this%cfg%imin,j,k)).ne.find(this%id(this%cfg%imin-1,j,k))) then
                     if (this%id(this%cfg%imin-1,j,k).gt.parent(this%id(this%cfg%imin,j,k))) then
                        parent(rootify2(this%id(this%cfg%imin-1,j,k)))=rootify2(parent(this%id(this%cfg%imin,j,k)))
                     else
                        parent(rootify2(parent(this%id(this%cfg%imin,j,k))))=rootify2(this%id(this%cfg%imin-1,j,k))
                     end if
                  end if
               end if
            end do; end do
         end if
         ! Handle jmin border
         if (this%cfg%jmin_.eq.this%cfg%jmin) then
            do k=this%cfg%kmin_,this%cfg%kmax_; do i=this%cfg%imin_,this%cfg%imax_
               ! If both border cell and neighboring ghost cell are filled
               if (this%id(i,this%cfg%jmin,k).gt.0.and.this%id(i,this%cfg%jmin-1,k).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (parent(this%id(i,this%cfg%jmin,k)).eq.this%id(i,this%cfg%jmin,k)) then
                     parent(rootify2(this%id(i,this%cfg%jmin,k)))=rootify2(this%id(i,this%cfg%jmin-1,k))
                  elseif (find(this%id(i,this%cfg%jmin,k)).ne.find(this%id(i,this%cfg%jmin-1,k))) then
                     if (this%id(i,this%cfg%jmin-1,k).gt.parent(this%id(i,this%cfg%jmin,k))) then
                        parent(rootify2(this%id(i,this%cfg%jmin-1,k)))=rootify2(parent(this%id(i,this%cfg%jmin,k)))
                     else
                        parent(rootify2(parent(this%id(i,this%cfg%jmin,k))))=rootify2(this%id(i,this%cfg%jmin-1,k))
                     end if
                  end if
               end if
            end do; end do
         end if
         ! Handle kmin border
         if (this%cfg%kmin_.eq.this%cfg%kmin) then
            do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
               ! If both border cell and neighboring ghost cell are filled
               if (this%id(i,j,this%cfg%kmin).gt.0.and.this%id(i,j,this%cfg%kmin-1).gt.0) then
                  ! If the border cell id already has a parent id, then union() the ghost cell id with the border cell id parent
                  if (parent(this%id(i,j,this%cfg%kmin)).eq.this%id(i,j,this%cfg%kmin)) then
                     parent(rootify2(this%id(i,j,this%cfg%kmin)))=rootify2(this%id(i,j,this%cfg%kmin-1))
                  elseif (find(this%id(i,j,this%cfg%kmin)).ne.find(this%id(i,j,this%cfg%kmin-1))) then
                     if (this%id(i,j,this%cfg%kmin-1).gt.parent(this%id(i,j,this%cfg%kmin))) then
                        parent(rootify2(this%id(i,j,this%cfg%kmin-1)))=rootify2(parent(this%id(i,j,this%cfg%kmin)))
                     else
                        parent(rootify2(parent(this%id(i,j,this%cfg%kmin))))=rootify2(this%id(i,j,this%cfg%kmin-1))
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
               parent_all(n)=find_all_2(n,n)
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
                     parent(rootify2(parent_own(n)))=rootify2(parent(n))
                     stop_=1
                  else if (find_parent.gt.find_parent_own) then
                     parent(rootify2(parent(n)))=rootify2(parent_own(n))
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
            do m=1,this%struct(n)%nnode
               this%id(this%struct(n)%node(1,m),this%struct(n)%node(2,m),this%struct(n)%node(3,m))=this%struct(n)%parent
            end do
         end do
         ! Update ghost cells
         call this%cfg%sync(this%id)
      end block boundary_handling
      
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
         integer :: size_now
         type(struct_type), dimension(:), allocatable :: tmp
         ! Check size of struct
         if (.not.allocated(this%struct)) then
            ! Allocate minimum storage
            allocate(this%struct(min_struct_size))
            ! Add first root
            this%nstruct_=1
            this%struct(this%nstruct_)%parent=1
            this%struct(this%nstruct_)%per=0
            this%struct(this%nstruct_)%nnode=0
         else
            ! Check if there is enough room for storing a new structure
            size_now=size(this%struct,dim=1)
            if (this%nstruct_.eq.size_now) then
               size_new=int(real(size_now,WP)*coeff_up)
               allocate(tmp(size_new))
               tmp(1:this%nstruct_)=this%struct
               tmp(this%nstruct_+1:)%parent=0
               tmp(this%nstruct_+1:)%nnode=0
               call move_alloc(tmp,this%struct)
            end if
            ! Add new root
            this%nstruct_=this%nstruct_+1
            this%struct(this%nstruct_)%parent=this%nstruct_
            this%struct(this%nstruct_)%per=0
            this%struct(this%nstruct_)%nnode=0
         end if
         x=this%nstruct_
      end function add
      
      !> This recursive function points parent to root and returns that root
      recursive function rootify2(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.parent(y)) then
            parent(y)=rootify2(parent(y))
            y=parent(y)
         end if
      end function rootify2

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
   
   
   !> Synchronize struct labels across all procs
   subroutine struct_sync(this)
      
   contains
      
      ! This function points the this%parent to root and returns that root
      recursive function find(a_x) result(a_y)
         implicit none
         integer :: a_x,a_y
         a_y=a_x
         if (a_y.ne.this%parent(a_y)) then
            this%parent(a_y)=find(this%parent(a_y))
            a_y=this%parent(a_y)
         end if
      end function find
      
      ! This subroutine joins two branches at their roots *it's not the same as union in struct_build*
      ! At the moment, it joints them at a_y and returns it but does not return the root
      subroutine union(a_x, a_y)
         implicit none
         integer, intent(in) :: a_x,a_y
         ! tbd: if a_x tree smaller/lower rank than a_y tree
         this%parent(find(a_x))=find(a_y)
      end subroutine union
      
      ! For this%parent_all array: this function points the this%parent to root and returns that root
      recursive function find_all(a_x) result(a_y)
         implicit none
         integer :: a_x,a_y
         a_y=a_x
         if (a_y.ne.this%parent_all(a_y)) then
            this%parent_all(a_y)=find_all(this%parent_all(a_y))
            a_y=this%parent_all(a_y)
         end if
      end function find_all
      
      ! Version that stops at the completion of a cycle
      recursive function find_all_2(a_x,a_starting_node) result(a_y)
         implicit none
         integer :: a_x,a_starting_node,a_y
         a_y=a_x
         if (a_y.ne.this%parent_all(a_y)) then
            if (this%parent_all(a_y).eq.a_starting_node) then
               a_y=this%parent_all(a_y)
               return
            else
               this%parent_all(a_y)=find_all_2(this%parent_all(a_y),a_starting_node)
               a_y=this%parent_all(a_y)
            end if
         end if
      end function find_all_2
      
      ! For this%parent_own array: this function points the this%parent to root and returns that root
      recursive function find_own(a_x) result(a_y)
         implicit none
         integer :: a_x,a_y
         a_y=a_x
         if (a_y.ne.this%parent_own(a_y)) then
            this%parent_own(a_y)=find_own(this%parent_own(a_y))
            a_y=this%parent_own(a_y)
         end if
      end function find_own
      
      
      
   end subroutine struct_sync
   
   
   !> Final update of struct labels in global array
   subroutine struct_label_update(this)
      implicit none
      class(ccl), intent(inout) :: this
      integer :: i,ii,jj,kk,n_node,tag_final
      
      ! Start at first_struct
      this%my_struct=>this%first_struct
      
      ! Initialize local vars
      n_node=0
      tag_final=0
      
      ! ! re-init global tag array
      ! this%id = 0
      
      ! Traverse structure list
      do while (associated(this%my_struct))
         n_node=this%my_struct%nnode
         ! Loop over nodes
         do i=1,n_node
            ! Update tag
            ii=this%my_struct%node(1,i)
            jj=this%my_struct%node(2,i)
            kk=this%my_struct%node(3,i)
            ! this%id(ii,jj,kk) = tag_final
            ! Update periodicity
            this%my_struct%per(1)=max(this%my_struct%per(1),this%idp(1,ii,jj,kk))
            this%my_struct%per(2)=max(this%my_struct%per(2),this%idp(2,ii,jj,kk))
            this%my_struct%per(3)=max(this%my_struct%per(3),this%idp(3,ii,jj,kk))
         end do
         this%my_struct=>this%my_struct%next
      end do
      
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
               this%film_edge_normal(:,i,j,k) = (c_local-c_filter)/(norm2(c_local-c_filter)+tiny(1.0_WP))
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
