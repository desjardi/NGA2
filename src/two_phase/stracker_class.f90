!> Two-phase structure tracking class
!> Based on cclabel object with greometric transport of index for persistent tracking of structures
module stracker_class
   use precision, only: WP
   use string,    only: str_medium
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
      ! Merging events
      integer :: nmerge
      integer, dimension(:,:), allocatable :: merge
      ! Merging2 events
      integer :: nmerge2
      integer, dimension(:,:), allocatable :: merge2
      ! ID of the structure that contains each cell
      integer, dimension(:,:,:), allocatable :: id
      ! Array of structures
      integer :: nstruct
      type(struct_type), dimension(:), allocatable :: struct
      ! Old ID of the structure that contains each cell
      integer, dimension(:,:,:), allocatable :: id_old
      ! Remapped old ID of the structure that contains each cell
      integer, dimension(:,:,:), allocatable :: id_rmp
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
      allocate(this%id_rmp(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%id_rmp=0
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
      
      ! Reset merge event storage
      reset_merge: block
         this%nmerge=0
         if (allocated(this%merge)) deallocate(this%merge)
         allocate(this%merge(1:2,1:min_struct_size)); this%merge=0
         this%nmerge2=0
         if (allocated(this%merge2)) deallocate(this%merge2)
         allocate(this%merge2(1:2,1:min_struct_size)); this%merge2=0
      end block reset_merge
      
      ! Remap id using VF geometry data to identify merge events
      remap_id: block
         use irl_fortran_interface
         integer, dimension(3) :: ind
         type(SepVM_type) :: my_SepVM
         real(WP), dimension(0:1) :: vols
         integer, dimension(:), allocatable :: ids
         integer :: i,j,k,n,nn
         integer :: nid,nobj,my_id
         ! Allocate ids storage
         nobj=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            nobj=max(count,getSize(this%vf%detailed_remap(i,j,k)))
         end do; end do; end do
         allocate(ids(nobj))
         ! Perform remapping step
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_; do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_; do i=this%vf%cfg%imin_,this%vf%cfg%imax_
            ! Check detailed remap geometry to get nobj
            nobj=getSize(this%vf%detailed_remap(i,j,k))
            ! Skip cell if no detailed remap geometry is available
            if (nobj.eq.0) cycle
            ! Zero out ids and reset counter
            nid=0; ids=0
            ! Initialize id_rmp
            this%id_rmp(i,j,k)=0
            ! Traverse every object and check for phase presence
            obj_loop: do n=1,nobj
               ! Get SepVM for nth object
               call getSepVMAtIndex(this%vf%detailed_remap(i,j,k),n-1,my_SepVM)
               ! Verify volume fraction for our phase of interest is non-zero
               vols(0)=getVolume(my_SepVM,0); vols(1)=getVolume(my_SepVM,1)
               if (vols(this%phase)/sum(vols).lt.VFlo) cycle
               ! Get cell index for nth object
               ind=this%vf%cfg%get_ijk_from_lexico(getTagForIndex(this%vf%detailed_remap(i,j,k),n-1))
               ! Get the local id
               my_id=this%id(ind(1),ind(2),ind(3))
               ! If my_id is zero, cycle
               if (my_id.eq.0) cycle obj_loop
               ! If my_id has already been encountered, cycle
               do nn=1,nid
                  if (ids(nn).eq.my_id) cycle obj_loop
               end do
               ! Increment the ids array
               nid=nid+1; ids(nid)=my_id
            end do obj_loop
            ! If no ids were encountered, done
            if (nid.eq.0) cycle
            ! Set id_rmp to smallest non-zero id encountered
            this%id_rmp(i,j,k)=minval(ids(1:nid))
            ! Create max(nid-1,0) merging events - this should be based on a binary search tree instead
            do n=1,nid
               call add_merge(this%id_rmp(i,j,k),ids(n))
            end do
         end do; end do; end do
         ! Deallocate ids
         deallocate(ids)
      end block remap_id

      ! Gather all merge events
      gather_merge: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_INTEGER
         integer, dimension(:), allocatable :: nmerge_proc
         integer, dimension(:,:), allocatable :: merge_all
         integer :: ierr,nmerge,n
         ! Sum all merge events
         allocate(nmerge_proc(1:this%vf%cfg%nproc)); nmerge_proc=0; nmerge_proc(this%vf%cfg%rank+1)=this%nmerge
         call MPI_ALLREDUCE(MPI_IN_PLACE,nmerge_proc,this%vf%cfg%nproc,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
         nmerge=sum(nmerge_proc)
         ! Gather all merge events
         allocate(merge_all(1:2,1:nmerge)); merge_all=0; merge_all(:,sum(nmerge_proc(1:this%vf%cfg%rank))+1:sum(nmerge_proc(1:this%vf%cfg%rank+1)))=this%merge(:,1:this%nmerge)
         call MPI_ALLREDUCE(MPI_IN_PLACE,merge_all,2*nmerge,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
         ! Reset merge data
         this%nmerge=0; deallocate(this%merge)
         allocate(this%merge(1:2,1:min_struct_size)); this%merge=0
         ! Compress merge_all array
         do n=1,nmerge
            call add_merge(merge_all(1,n),merge_all(2,n))
         end do
         ! Deallocate
         deallocate(nmerge_proc,merge_all)
      end block gather_merge

      ! Execute all merges
      execute_merge: block
         integer :: i,j,k,n,m
         ! Traverse merge list
         do n=1,this%nmerge
            ! Loop over full domain and update id based on that merge event
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_; do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_; do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
               if (this%id_rmp(i,j,k).eq.this%merge(2,n)) this%id_rmp(i,j,k)=this%merge(1,n)
            end do; end do; end do
            ! Loop over remaining merge events and update id based on that merge event
            do m=1,this%merge
               if (this%merge(1,m).eq.this%merge(2,n)) this%merge(1,m)=this%merge(1,n)
               if (this%merge(2,m).eq.this%merge(2,n)) this%merge(2,m)=this%merge(1,n)
               i=minval(this%merge(:,m)); j=maxval(this%merge(:,m)); this%merge(:,m)=[i,j]
            end do
         end do
      end block execute_merge
      
      ! Perform a CCL build - same_label is false if id_rmp are different and non-zero
      call this%build(make_label)
      
      ! Resolve persistent id:
      ! For each new structure, build a list of all distinct id_rmp that are contained
      ! Zeros are ignored here - if more than one non-zero id_rmp value, then another
      ! merging step is needed
      new_id: block
         integer :: n,nn,nnn,nid,my_id,nobj
         integer, dimension(:), allocatable :: ids
         ! Allocate maximum storage for id_rmp values
         nobj=0
         do n=1,this%nstruct
            nobj=max(nobj,this%struct(n)%n_)
         end do
         allocate(ids(nobj))
         ! Traverse each new structure and gather list of potential id
         do n=1,this%nstruct
            ! Zero out ids and reset counter
            nid=0; ids=0
            ! Loop over local cells in the structure
            str_loop: do nn=1,this%struct(n)%n_
               ! Local cell id_rmp value
               my_id=this%id_rmp(this%struct(n)%map(1,nn),this%struct(n)%map(2,nn),this%struct(n)%map(3,nn))
               ! Ignore zeros
               if (my_id.eq.0) cycle
               ! If my_id has already been encountered, cycle
               do nnn=1,nid
                  if (ids(nnn).eq.my_id) cycle str_loop
               end do
               ! Increment the ids array
               nid=nid+1; ids(nid)=my_id
            end do str_loop
            ! If no ids were encountered, done - we may need a brand new id
            if (nid.eq.0) cycle
            ! Set id_rmp to smallest non-zero id encountered
            this%struct(n)%id=minval(ids(1:nid))
            ! Create max(nid-1,0) merging events - this should be based on a binary search tree instead
            do nnn=1,nid
               call add_merge2(this%struct(n)%id,ids(nnn))
            end do
         end do
      end block new_id

      ! Gather all merge2 events
      gather_merge2: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_INTEGER
         integer, dimension(:), allocatable :: nmerge2_proc
         integer, dimension(:,:), allocatable :: merge2_all
         integer :: ierr,nmerge2,n
         ! Sum all merge2 events
         allocate(nmerge2_proc(1:this%vf%cfg%nproc)); nmerge2_proc=0; nmerge2_proc(this%vf%cfg%rank+1)=this%nmerge2
         call MPI_ALLREDUCE(MPI_IN_PLACE,nmerge2_proc,this%vf%cfg%nproc,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
         nmerge2=sum(nmerge2_proc)
         ! Gather all merge2 events
         allocate(merge2_all(1:2,1:nmerge2)); merge2_all=0; merge2_all(:,sum(nmerge2_proc(1:this%vf%cfg%rank))+1:sum(nmerge2_proc(1:this%vf%cfg%rank+1)))=this%merge2(:,1:this%nmerge2)
         call MPI_ALLREDUCE(MPI_IN_PLACE,merge2_all,2*nmerge2,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)
         ! Reset merge2 data
         this%nmerge2=0; deallocate(this%merge2)
         allocate(this%merge2(1:2,1:min_struct_size)); this%merge2=0
         ! Compress merge2_all array
         do n=1,nmerge2
            call add_merge2(merge2_all(1,n),merge2_all(2,n))
         end do
         ! Deallocate
         deallocate(nmerge2_proc,merge2_all)
      end block gather_merge2

      ! Execute all merge2
      execute_merge2: block
         integer :: i,j,n,m
         ! Traverse merge2 list
         do n=1,this%nmerge2
            ! Loop over structures and update id based on merge2 events
            do nn=1,this%nstruct
               if (this%struct(nn)%id.eq.this%merge2(2,n)) this%struct(nn)%id=this%merge2(1,n)
            end do
            ! Loop over remaining merge events and update id based on that merge2 event
            do m=1,this%merge2
               if (this%merge2(1,m).eq.this%merge2(2,n)) this%merge2(1,m)=this%merge2(1,n)
               if (this%merge2(2,m).eq.this%merge2(2,n)) this%merge2(2,m)=this%merge2(1,n)
               i=minval(this%merge2(:,m)); j=maxval(this%merge2(:,m)); this%merge2(:,m)=[i,j]
            end do
         end do
      end block execute_merge2
      
      ! Now deal with splits - case where two structs share the same id
      
      
      ! Finally, one last pass to give unique id to structs with id=0
      
      
   contains
      
      !> Subroutine that adds a new merge event (id1<id2 required)
      subroutine add_merge(id1,id2)
         implicit none
         integer, intent(in) :: id1,id2
         integer :: n,size_now,size_new
         integer, dimension(:,:), allocatable :: tmp
         ! Skip self-merge
         if (id1.eq.id2) return
         ! Skip redundant merge
         do n=1,this%nmerge
            if (all(this%merge(:,n).eq.[id1,id2])) return
         end do
         ! Create new merge event
         size_now=size(this%merge,dim=2)
         if (this%nmerge.eq.size_now) then
            size_new=int(real(size_now,WP)*coeff_up)
            allocate(tmp(1:2,1:size_new))
            tmp(:,1:this%nmerge)=this%merge
            tmp(:,this%nmerge+1:)=0
            call move_alloc(tmp,this%merge)
         end if
         this%nmerge=this%nmerge+1
         this%merge(:,this%nmerge)=[id1,id2]
      end subroutine add_merge

      !> Subroutine that adds a new merge2 event (id1<id2 required)
      subroutine add_merge2(id1,id2)
         implicit none
         integer, intent(in) :: id1,id2
         integer :: n,size_now,size_new
         integer, dimension(:,:), allocatable :: tmp
         ! Skip self-merge2
         if (id1.eq.id2) return
         ! Skip redundant merge2
         do n=1,this%nmerge2
            if (all(this%merge2(:,n).eq.[id1,id2])) return
         end do
         ! Create new merge2 event
         size_now=size(this%merge2,dim=2)
         if (this%nmerge2.eq.size_now) then
            size_new=int(real(size_now,WP)*coeff_up)
            allocate(tmp(1:2,1:size_new))
            tmp(:,1:this%nmerge2)=this%merge2
            tmp(:,this%nmerge2+1:)=0
            call move_alloc(tmp,this%merge2)
         end if
         this%nmerge2=this%nmerge2+1
         this%merge2(:,this%nmerge2)=[id1,id2]
      end subroutine add_merge2
      
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
            tmp(n)%id=0 !< This is different from CCLABEL, our id is initialized to zero here
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
      
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         ! True by default
         same_label=.true.
         ! Exception is if pair carries distinct non-zero remapped ids
         if (this%id_rmp(i1,j1,k1)*this%id_rmp(i2,j2,k2).gt.0.and.this%id_rmp(i1,j1,k1).ne.this%id_rmp(i2,j2,k2)) same_label=.false.
      end function same_label

   end subroutine build
   
   
end module stracker_class