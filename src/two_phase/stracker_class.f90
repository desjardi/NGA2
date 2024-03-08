!> Two-phase structure tracking class
!> Based on cclabel object with greometric transport of index for persistent tracking of structures
module stracker_class
   use precision, only: WP
   use vfs_class, only: vfs
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: stracker,make_label_ftype
   
   
   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%
   
   
   !> Structure object
   type :: struct_type
      integer :: id                                       !< Persistent ID of struct
      integer :: n_                                       !< Number of local cells contained in struct
      integer, dimension(:,:), allocatable :: map         !< List of cells contained in struct, dimension(1:3,1:n_)
      integer, dimension(3) :: per                        !< Periodicity array - per(dim)=1 if structure is periodic in dim direction
   end type struct_type
   
   
   !> Structure tracker object definition
   type :: stracker
      ! Pointer to our VF solver
      class(vfs), pointer :: vf
      ! This is the name of the stracker
      character(len=str_medium) :: name='UNNAMED_STRACKER'
      ! Phase containing the tracked structures
      integer :: phase                !< 0 is liquid, 1 is gas
      ! ID of the structure that contains each cell
      integer, dimension(:,:,:), allocatable :: id
      ! Array of structures
      integer :: nstruct
      type(struct_type), dimension(:), allocatable :: struct
      ! Old ID of the structure that contains each cell
      integer, dimension(:,:,:), allocatable :: id_old
      ! Array of structures
      integer :: nstruct_old
      type(struct_type), dimension(:), allocatable :: struct_old
   contains
      procedure :: initialize         !< Initialization of stracker based on a provided VFS object
      procedure, private :: build     !< Private cclabel step without persistent id
      procedure :: advance            !< Perform cclabel step with persistent id
   end type stracker
   
   
   !> Type of the make_label function used to generate a structure - user-provided
   abstract interface
      logical function make_label_ftype(ind1,ind2,ind3)
         integer, intent(in) :: ind1,ind2,ind3
      end function make_label_ftype
   end interface
   
   
contains
   
   
   !> Structure tracker initialization
   subroutine initialize(this,vf,phase,make_label,name)
      use messager,  only: die
      use vfs_class, only: remap_storage
      implicit none
      class(stracker), intent(inout) :: this
      class(vfs), target, intent(in) :: vf
      integer, intent(in) :: phase
      procedure(make_label_ftype) :: make_label
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Set the phase
      this%phase=phase
      ! Point to our vfs object
      this%vf=>vf
      ! Check vfs object uses cell remap with storage
      if (this%vf%transport_method.ne.remap_storage) call die('[stracker initialize] Our vfs needs to use full cell remap with storage')
      ! Allocate and initialize ID array
      allocate(this%id    (this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%id    =0
      allocate(this%id_old(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%id_old=0
      ! Zero structures
      this%nstruct    =0
      this%nstruct_old=0
      ! Perform a first CCL
      call this%build(make_label)
   end subroutine initialize
   

   !> Structure tracking step
   subroutine advance(this,make_label)
      implicit none
      class(stracker), intent(inout) :: this
      procedure(make_label_ftype) :: make_label
      integer, dimension(:,:,:), allocatable :: id_remap
      
      ! Copy old structure data - may not be needed
      copy_to_old: block
         integer :: n
         this%id_old=this%id
         this%nstruct_old=this%nstruct
         if (allocated(this%struct_old)) then
            do n=1,size(this%struct_old,dim=1)
               if (allocated(this%struct_old(n)%map)) deallocate(this%struct_old(n)%map)
            end do
            deallocate(this%struct_old)
         end if
         allocate(this%struct_old(this%nstruct_old))
         do n=1,size(this%struct,dim=1)
            this%struct_old(n)%id =this%struct(n)%id
            this%struct_old(n)%per=this%struct(n)%per
            this%struct_old(n)%n_ =this%struct(n)%n_
            allocate(this%struct_old(n)%map(this%struct_old(n)%n_))
            this%struct_old(n)%map=this%struct(n)%map
         end do
      end block copy_to_old
      
      ! Remap id using VF geometry data to identify merge events
      remap_id: block
         use irl_fortran_interface
         integer :: i,j,k,n
         integer, dimension(3) :: ind
         type(SepVM_type) :: my_SepVM
         integer :: my_id
         ! Allocate remapped id array
         allocate(id_remap(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); id_remap=0
         ! Perform remapping step
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            ! Check if detailed remap geometry is available
            if (getSize(this%vf%detailed_remap(i,j,k)).gt.0) then
               ! It is, so analyze it for merge events
               do n=0,getSize(this%vf%detailed_remap(i,j,k))-1
                  ! Get SepVM for nth object
                  call getSepVMAtIndex(this%vf%detailed_remap(i,j,k),n,my_SepVM)
                  ! Verify volume of tracked cell is non-zero
                  if (getVolume(my_SepVM,this%phase).eq.0.0_WP) cycle
                  ! Get cell index for nth object
                  ind=this%vf%cfg%get_ijk_from_lexico(getTagForIndex(this%vf%detailed_remap(i,j,k),n))
                  ! Get the id
                  my_id=this%id(ind(1),ind(2),ind(3))
               end do
               ! Get remapped id
               id_remap(i,j,k)=my_id
            end if
         end do; end do; end do
      end block remap_id

      ! Perform a CCL build
      call this%build(make_label)

      ! Resolve persistent id
      new_id: block
      
      end block new_id
      
   end subroutine advance

   
   !> Build structure using the user-set test function
   subroutine build(this,make_label)
      implicit none
      class(stracker), intent(inout) :: this
      procedure(make_label_ftype) :: make_label
      integer :: nstruct_,stmin,stmax
      integer, dimension(:,:,:,:), allocatable :: idp          !< Periodicity treatment
      integer, dimension(:), allocatable :: parent             !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_all         !< Resolving structure id across procs
      integer, dimension(:), allocatable :: parent_own         !< Resolving structure id across procs
      
      ! Cleanup of storage before doing CCL
      cleanup: block
         integer :: n
         ! First loop over all structures and deallocate maps
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
         ! Then allocate struct to a default size
         nstruct_=0
         allocate(this%struct(min_struct_size))
         this%struct(:)%id    =0
         this%struct(:)%per(1)=0
         this%struct(:)%per(2)=0
         this%struct(:)%per(3)=0
         this%struct(:)%n_=0
      end block cleanup
      
      ! Allocate periodicity work array
      allocate(idp(3,this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); idp=0
      
      ! Perform a first pass to build proc-local structures and corresponding tree
      first_pass: block
         integer :: i,j,k,ii,jj,kk,dim
         integer, dimension(3) :: pos
         ! Perform local loop
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            ! Find next cell in a structure
            if (make_label(i,j,k)) then
               ! Loop through one-sided neighbors
               do dim=1,3
                  pos=0; pos(dim)=-1
                  ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                  ! Check if neighbor is labeled
                  if (this%id(ii,jj,kk).gt.0) then
                     ! Neighbor is labeled, but are we?
                     if (this%id(i,j,k).ne.0) then
                        ! We already have a label, perform a union of both labels
                        if (same_label(i,j,k,ii,jj,kk)) this%id(i,j,k)=union_struct(this%id(i,j,k),this%id(ii,jj,kk))
                     else
                        ! We don't have a label, check if we take the neighbor's label
                        if (same_label(i,j,k,ii,jj,kk)) then
                           this%id(i,j,k)=this%id(ii,jj,kk)
                        else
                           this%id(i,j,k)=add()
                        end if
                     end if
                  end if
               end do
               ! If no neighbor was labeled, we need a new structure
               if (this%id(i,j,k).eq.0) this%id(i,j,k)=add()
               ! Identify periodicity cases
               if (this%vf%cfg%xper.and.i.eq.this%vf%cfg%imax) this%struct(this%id(i,j,k))%per(1)=1
               if (this%vf%cfg%yper.and.j.eq.this%vf%cfg%jmax) this%struct(this%id(i,j,k))%per(2)=1
               if (this%vf%cfg%zper.and.k.eq.this%vf%cfg%kmax) this%struct(this%id(i,j,k))%per(3)=1
               idp(:,i,j,k)=this%struct(this%id(i,j,k))%per
            end if
         end do; end do; end do
      end block first_pass
      
      ! Now collapse the tree, count the cells and resolve periodicity in each structure
      collapse_tree: block
         integer :: i,j,k
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            if (this%id(i,j,k).gt.0) then
               this%id(i,j,k)=rootify_struct(this%id(i,j,k))
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
         integer :: i,j,k,n,ierr
         integer, dimension(:), allocatable :: my_nstruct,all_nstruct,idmap
         type(struct_type), dimension(:), allocatable :: tmp
         ! Count exact number of local structures
         nstruct_=0
         do n=1,size(this%struct,dim=1)
            if (this%struct(n)%n_.gt.0) nstruct_=nstruct_+1
         end do
         ! Gather this info to ensure unique index
         allocate( my_nstruct(0:this%vf%cfg%nproc-1)); my_nstruct=0; my_nstruct(this%vf%cfg%rank)=nstruct_
         allocate(all_nstruct(0:this%vf%cfg%nproc-1)); call MPI_ALLREDUCE(my_nstruct,all_nstruct,this%vf%cfg%nproc,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
         stmin=1
         if (this%vf%cfg%rank.gt.0) stmin=stmin+sum(all_nstruct(0:this%vf%cfg%rank-1))
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
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            if (this%id(i,j,k).gt.0) this%id(i,j,k)=idmap(this%id(i,j,k))
         end do; end do; end do
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
      
      ! Fill out the node map
      node_map: block
         integer :: i,j,k
         integer, dimension(:), allocatable :: counter
         allocate(counter(stmin:stmax)); counter=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
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
         ! Fill global lineage with selves
         do n=1,this%nstruct
            parent(n)=n
         end do
         ! Synchronize id array
         call this%vf%cfg%sync(this%id)
         ! Handle imin_ border
         if (this%vf%cfg%imin_.ne.this%vf%cfg%imin) then
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               if (this%id(this%vf%cfg%imin_,j,k).gt.0.and.this%id(this%vf%cfg%imin_-1,j,k).gt.0) then
                  if (same_label(this%vf%cfg%imin_,j,k,this%vf%cfg%imin_-1,j,k)) call union_parent(this%id(this%vf%cfg%imin_,j,k),this%id(this%vf%cfg%imin_-1,j,k))
               end if
            end do; end do
         end if
         ! Handle jmin_ border
         if (this%vf%cfg%jmin_.ne.this%vf%cfg%jmin) then
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               if (this%id(i,this%vf%cfg%jmin_,k).gt.0.and.this%id(i,this%vf%cfg%jmin_-1,k).gt.0) then
                  if (same_label(i,this%vf%cfg%jmin_,k,i,this%vf%cfg%jmin_-1,k)) call union_parent(this%id(i,this%vf%cfg%jmin_,k),this%id(i,this%vf%cfg%jmin_-1,k))
               end if
            end do; end do
         end if
         ! Handle kmin_ border
         if (this%vf%cfg%kmin_.ne.this%vf%cfg%kmin) then
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               if (this%id(i,j,this%vf%cfg%kmin_).gt.0.and.this%id(i,j,this%vf%cfg%kmin_-1).gt.0) then
                  if (same_label(i,j,this%vf%cfg%kmin_,i,j,this%vf%cfg%kmin_-1)) call union_parent(this%id(i,j,this%vf%cfg%kmin_),this%id(i,j,this%vf%cfg%kmin_-1))
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
            ! Set self-parents to huge(1)
            do n=1,this%nstruct
               if (parent(n).eq.n) parent(n)=huge(1)
            end do
            ! Take global min
            call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%vf%cfg%comm,ierr)
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
            call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%vf%cfg%comm,ierr)
         end do
         ! Update this%struct%id by pointing all parents to root and update id
         do n=stmin,stmax
            this%struct(n)%id=rootify_parent(parent(n))
            do m=1,this%struct(n)%n_
               this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=this%struct(n)%id
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
         do n=stmin,stmax
            ownper(:,n)=this%struct(n)%per
         end do
         ! Communicate per
         call MPI_ALLREDUCE(ownper,allper,3*this%nstruct,MPI_INTEGER,MPI_MAX,this%vf%cfg%comm,ierr)
         ! Update parent per
         do n=1,this%nstruct
            allper(:,parent(n))=max(allper(:,parent(n)),allper(:,n))
         end do
         ! Update idp array
         do n=stmin,stmax
            do m=1,this%struct(n)%n_
               idp(:,this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=allper(:,this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m)))
            end do
         end do
         ! Clean up
         deallocate(ownper,allper)
      end block periodicity_update
      
      ! One more pass for domain boundaries
      boundary_handling: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_INTEGER
         integer :: i,j,k,stop_global,stop_,counter,n,m,ierr,find_parent,find_parent_own
         ! Handle imin border
         if (this%vf%cfg%imin_.eq.this%vf%cfg%imin) then
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               if (this%id(this%vf%cfg%imin_,j,k).gt.0.and.this%id(this%vf%cfg%imin_-1,j,k).gt.0) then
                  if (same_label(this%vf%cfg%imin_,j,k,this%vf%cfg%imin_-1,j,k)) call union_parent(this%id(this%vf%cfg%imin_,j,k),this%id(this%vf%cfg%imin_-1,j,k))
               end if
            end do; end do
         end if
         ! Handle jmin border
         if (this%vf%cfg%jmin_.eq.this%vf%cfg%jmin) then
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               if (this%id(i,this%vf%cfg%jmin_,k).gt.0.and.this%id(i,this%vf%cfg%jmin_-1,k).gt.0) then
                  if (same_label(i,this%vf%cfg%jmin_,k,i,this%vf%cfg%jmin_-1,k)) call union_parent(this%id(i,this%vf%cfg%jmin_,k),this%id(i,this%vf%cfg%jmin_-1,k))
               end if
            end do; end do
         end if
         ! Handle kmin border
         if (this%vf%cfg%kmin_.eq.this%vf%cfg%kmin) then
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               if (this%id(i,j,this%vf%cfg%kmin_).gt.0.and.this%id(i,j,this%vf%cfg%kmin_-1).gt.0) then
                  if (same_label(i,j,this%vf%cfg%kmin_,i,j,this%vf%cfg%kmin_-1)) call union_parent(this%id(i,j,this%vf%cfg%kmin_),this%id(i,j,this%vf%cfg%kmin_-1))
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
            ! Set self-parents to huge(1)
            do n=1,this%nstruct
               if (parent(n).eq.n) parent(n)=huge(1)
            end do
            ! Take global min
            call MPI_ALLREDUCE(parent,parent_all,this%nstruct,MPI_INTEGER,MPI_MIN,this%vf%cfg%comm,ierr)
            ! Set self-parents back to selves
            do n=1,this%nstruct
               if (parent_all(n).eq.huge(1)) parent_all(n)=n
            end do
            ! Flatten trees
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
                  find_parent_own=rootify_parent(parent_own(n))
                  find_parent    =rootify_parent(parent(n))
                  if (find_parent_own.ne.find_parent) then
                     call union_parent(find_parent,find_parent_own)
                     stop_=1
                  end if
               end if
            end do
            ! Check if we did some changes
            call MPI_ALLREDUCE(stop_,stop_global,1,MPI_INTEGER,MPI_MAX,this%vf%cfg%comm,ierr)
         end do
         ! Update this%struct%id and point all parents to root and update id
         do n=stmin,stmax
            this%struct(n)%id=rootify_parent(parent(n))
            do m=1,this%struct(n)%n_
               this%id(this%struct(n)%map(1,m),this%struct(n)%map(2,m),this%struct(n)%map(3,m))=this%struct(n)%id
            end do
         end do
         ! Update ghost cells
         call this%vf%cfg%sync(this%id)
         ! Clean up parent info
         deallocate(parent,parent_all,parent_own)
      end block boundary_handling
      
      ! Now we need to compact the data based on id only
      compact_struct: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
         integer :: i,j,k,n,nn,ierr,count
         integer, dimension(:), allocatable :: my_idmap,idmap,counter
         type(struct_type), dimension(:), allocatable :: tmp
         ! Prepare global id map
         allocate(my_idmap(1:this%nstruct)); my_idmap=0
         allocate(   idmap(1:this%nstruct));    idmap=0
         ! Traverse id array and tag used id values
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            if (this%id(i,j,k).gt.0) my_idmap(this%id(i,j,k))=1
         end do; end do; end do
         call MPI_ALLREDUCE(my_idmap,idmap,this%nstruct,MPI_INTEGER,MPI_MAX,this%vf%cfg%comm,ierr)
         deallocate(my_idmap)
         ! Count number of used structures and create the map
         this%nstruct=sum(idmap)
         count=0
         do n=1,size(idmap,dim=1)
            if (idmap(n).gt.0) then
               count=count+1
               idmap(n)=count
            end if
         end do
         ! Rename all structures
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            if (this%id(i,j,k).gt.0) this%id(i,j,k)=idmap(this%id(i,j,k))
         end do; end do; end do
         call this%vf%cfg%sync(this%id)
         ! Allocate temporary storage for structure
         allocate(tmp(this%nstruct))
         allocate(counter(this%nstruct)); counter=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            if (this%id(i,j,k).gt.0) counter(this%id(i,j,k))=counter(this%id(i,j,k))+1
         end do; end do; end do
         do n=1,this%nstruct
            tmp(n)%id=n
            tmp(n)%n_=counter(n)
            allocate(tmp(n)%map(1:3,1:tmp(n)%n_))
         end do
         ! Transfer periodicity info
         do n=stmin,stmax
            if (idmap(n).gt.0) then
               tmp(idmap(n))%per=this%struct(n)%per
            end if
         end do
         deallocate(idmap)
         ! Store the map
         counter=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            if (this%id(i,j,k).gt.0) then
               counter(this%id(i,j,k))=counter(this%id(i,j,k))+1
               tmp(this%id(i,j,k))%map(:,counter(this%id(i,j,k)))=[i,j,k]
            end if
         end do; end do; end do
         deallocate(counter)
         ! Transfer allocation
         call move_alloc(tmp,this%struct)
         ! Final pass to fix periodicity info
         do n=1,this%nstruct
            do nn=1,this%struct(n)%n_
               i=this%struct(n)%map(1,nn)
               j=this%struct(n)%map(2,nn)
               k=this%struct(n)%map(3,nn)
               this%struct(n)%per(1)=max(this%struct(n)%per(1),idp(1,i,j,k))
               this%struct(n)%per(2)=max(this%struct(n)%per(2),idp(2,i,j,k))
               this%struct(n)%per(3)=max(this%struct(n)%per(3),idp(3,i,j,k))
            end do
         end do
         deallocate(idp)
      end block compact_struct
      
   contains
      
      !> This recursive function that points the lineage of a structure to its root and returns that root
      recursive function rootify_struct(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y
         y=x
         if (y.ne.this%struct(y)%id) then
            this%struct(y)%id=rootify_struct(this%struct(y)%id)
            y=this%struct(y)%id
         end if
      end function rootify_struct
      
      !> This function joins two structures at their roots (the smallest root is chosen and returned)
      function union_struct(x,y) result(rmin)
         implicit none
         integer, intent(in) :: x,y
         integer :: rx,ry,rmin,rmax
         rx=rootify_struct(x); ry=rootify_struct(y)
         rmin=min(rx,ry); rmax=max(rx,ry)
         this%struct(rmax)%id=rmin
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
            tmp(nstruct_+1:)%id    =0
            tmp(nstruct_+1:)%per(1)=0
            tmp(nstruct_+1:)%per(2)=0
            tmp(nstruct_+1:)%per(3)=0
            tmp(nstruct_+1:)%n_=0
            call move_alloc(tmp,this%struct)
         end if
         ! Add new root
         nstruct_=nstruct_+1
         this%struct(nstruct_)%id =nstruct_
         this%struct(nstruct_)%per=0
         this%struct(nstruct_)%n_ =0
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
   
   
end module stracker_class