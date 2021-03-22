!> Coupler concept is defined here: it takes in two pgrid
!> objects and builds the communication and interpolation
!> layer to exchange data between them.
module coupler_class
   use precision,      only: WP
   use string,         only: str_medium
   use pgrid_class,    only: pgrid
   use mpi_f08
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: coupler
   
   !> Coupler object definition
   type :: coupler
      ! This is the name of the coupler
      character(len=str_medium) :: name='UNNAMED_CPL'     !< Coupler name (default=UNNAMED_CPL)
      ! This is our communication information
      type(MPI_Comm) :: comm                              !< Intracommunicator over the union of both groups
      type(MPI_Group) :: sgrp,dgrp,grp                    !< Source and destination groups and their union
      integer :: nproc                                    !< Number of processors
      integer :: rank                                     !< Processor grid rank
      logical :: amRoot                                   !< Am I root for the coupler?
      integer :: sroot                                    !< Rank of src grid root on union group
      integer :: droot                                    !< Rank of dst grid root on union group
      ! These are our two pgrids
      type(pgrid), pointer :: src=>NULL()                 !< Source grid
      type(pgrid), pointer :: dst=>NULL()                 !< Destination grid
      ! Logicals to help us know if we have received a src or dst grid
      logical :: got_src=.false.                          !< Were we given a src grid
      logical :: got_dst=.false.                          !< Were we given a dst grid
      ! Rank map for dst grid
      integer :: dnproc,dnpx,dnpy,dnpz                    !< Destination grid partitioning
      integer, dimension(:,:,:), allocatable :: rankmap   !< Processor coordinate to union group rank map
   contains
      procedure :: initialize                             !< Routine that prepares all interpolation metrics from src to dst
      procedure :: set_src                                !< Routine that sets the source grid
      procedure :: set_dst                                !< Routine that sets the destination grid
   end type coupler
   
   
   !> Declare coupler constructor
   interface coupler
      procedure construct_from_two_groups
   end interface coupler
   
contains
   
   
   !> Coupler constructor from two groups
   function construct_from_two_groups(src_grp,dst_grp,name) result(self)
      use messager, only: die
      use parallel, only: comm
      implicit none
      type(coupler) :: self
      type(MPI_Group), intent(in) :: src_grp,dst_grp
      character(len=*), intent(in) :: name
      integer, dimension(1) :: rankin,rankout
      integer :: ierr
      
      ! Set name for the coupler
      self%name=trim(adjustl(name))
      
      ! Build group union
      self%sgrp=src_grp
      self%dgrp=dst_grp
      call MPI_GROUP_UNION(self%sgrp,self%dgrp,self%grp,ierr)
      
      ! Gather some info for communication
      call MPI_GROUP_SIZE(self%grp,self%nproc,ierr)
      if (self%nproc.eq.0) call die('[coupler constructor] Somehow the union of both groups is of size zero')
      call MPI_GROUP_RANK(self%grp,self%rank ,ierr)
      if (self%rank.eq.MPI_UNDEFINED) call die('[coupler constructor] All processors that call the constructor need to be in one of the two groups')
      
      ! Create intracommunicator for the new group
      call MPI_COMM_CREATE_GROUP(comm,self%grp,0,self%comm,ierr)
      
      ! Find roots for both grids on the shared communicator
      rankin=0; call MPI_GROUP_TRANSLATE_RANKS(self%sgrp,1,rankin,self%grp,rankout,ierr); self%sroot=rankout(1)
      rankin=0; call MPI_GROUP_TRANSLATE_RANKS(self%dgrp,1,rankin,self%grp,rankout,ierr); self%droot=rankout(1)
      
      ! Set coupler root to src root
      self%amRoot=(self%rank.eq.self%sroot)
      
   end function construct_from_two_groups
   
   
   !> Set the source grid - to be called by processors in src_group
   subroutine set_src(this,pg)
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      ! Point to the grid
      this%src=>pg
      ! Set a flag
      this%got_src=.true.
   end subroutine set_src
   
   
   !> Set the destination grid - to be called by processors in dst_group
   subroutine set_dst(this,pg)
      implicit none
      class(coupler), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      ! Point to the grid
      this%dst=>pg
      ! Set a flag
      this%got_dst=.true.
   end subroutine set_dst
   
   
   !> Prepare interpolation metrics from src to dst
   subroutine initialize(this)
      implicit none
      class(coupler), intent(inout) :: this
      
      
      ! First step is to make destination grid available to all
      share_grid: block
         use sgrid_class, only: sgrid
         use parallel,    only: MPI_REAL_WP
         character(len=str_medium) :: simu_name
         real(WP), dimension(:), allocatable :: x
         real(WP), dimension(:), allocatable :: y
         real(WP), dimension(:), allocatable :: z
         logical :: xper,yper,zper
         integer :: no,nx,ny,nz,coord,ierr
         
         ! Destination root process extracts its own sgrid
         if (this%rank.eq.this%droot) then
            simu_name=this%dst%name
            coord=this%dst%coordsys
            xper=this%dst%xper
            yper=this%dst%yper
            zper=this%dst%zper
            nx=this%dst%nx
            ny=this%dst%ny
            nz=this%dst%nz
            no=this%dst%no
            this%dnproc=this%dst%nproc
            this%dnpx=this%dst%npx
            this%dnpy=this%dst%npy
            this%dnpz=this%dst%npz
         end if
         
         ! Then it broadcasts it to our group
         call MPI_BCAST(simu_name,len(simu_name),MPI_CHARACTER,this%droot,this%comm,ierr)
         call MPI_BCAST(coord    ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(xper     ,1             ,MPI_LOGICAL  ,this%droot,this%comm,ierr)
         call MPI_BCAST(yper     ,1             ,MPI_LOGICAL  ,this%droot,this%comm,ierr)
         call MPI_BCAST(zper     ,1             ,MPI_LOGICAL  ,this%droot,this%comm,ierr)
         call MPI_BCAST(nx       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(ny       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(nz       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(no       ,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpx,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpy,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpz,1             ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnproc,1           ,MPI_INTEGER  ,this%droot,this%comm,ierr)
         
         ! Allocate x/y/z, fill it, and bcast
         allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1))
         if (this%rank.eq.this%droot) then
            x(1:nx+1)=this%dst%x(this%dst%imin:this%dst%imax+1)
            y(1:ny+1)=this%dst%y(this%dst%jmin:this%dst%jmax+1)
            z(1:nz+1)=this%dst%z(this%dst%kmin:this%dst%kmax+1)
         end if
         call MPI_BCAST(x,nx+1,MPI_REAL_WP,this%droot,this%comm,ierr)
         call MPI_BCAST(y,ny+1,MPI_REAL_WP,this%droot,this%comm,ierr)
         call MPI_BCAST(z,nz+1,MPI_REAL_WP,this%droot,this%comm,ierr)
         
         ! Finish creating the sgrid
         if (.not.this%got_dst) then
            allocate(this%dst)
            this%dst%sgrid=sgrid(coord,no,x,y,z,xper,yper,zper,trim(adjustl(simu_name)))
         end if
         
         ! Deallocate
         deallocate(x,y,z)
         
      end block share_grid
      
      
      ! Second step is to make destination partition map available to all
      share_partition: block
         integer :: ierr
         integer, dimension(:), allocatable :: diproc,djproc,dkproc
         
         ! Destination root process extracts partition
         if (this%rank.eq.this%droot) then
            this%dnproc=this%dst%nproc
            this%dnpx=this%dst%npx
            this%dnpy=this%dst%npy
            this%dnpz=this%dst%npz
         end if
         
         ! Broadcast it to our group
         call MPI_BCAST(this%dnpx,  1,MPI_INTEGER,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpy,  1,MPI_INTEGER,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnpz,  1,MPI_INTEGER,this%droot,this%comm,ierr)
         call MPI_BCAST(this%dnproc,1,MPI_INTEGER,this%droot,this%comm,ierr)
         
         ! Allocate the destination rankmap
         allocate(this%rankmap(this%dnpx,this%dnpy,this%dnpz))
         
         ! Prepare communication arrays
         allocate(diproc(this%nproc),djproc(this%nproc),dkproc(this%nproc))
         
         ! Provide a default i/j/kproc to processors without dst grid
         if (.not.this%got_dst) then
            this%dst%iproc=-1
            this%dst%jproc=-1
            this%dst%kproc=-1
         end if
         
         ! Allgather the rank->(iproc,jproc,kproc) info
         call MPI_ALLGATHER(this%dst%iproc,1,MPI_INTEGER,diproc,1,MPI_INTEGER,this%dst%comm,ierr)
         call MPI_ALLGATHER(this%dst%jproc,1,MPI_INTEGER,djproc,1,MPI_INTEGER,this%dst%comm,ierr)
         call MPI_ALLGATHER(this%dst%kproc,1,MPI_INTEGER,dkproc,1,MPI_INTEGER,this%dst%comm,ierr)
         
         ! Finally, flip the rankmap data
         
         
         ! Deallocate communication arrays
         deallocate(diproc,djproc,dkproc)
         
      end block share_partition
      
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (this%amRoot) then
            write(message,'("Coupler [",a,"] from pgrid [",a,"] to pgrid [",a,"]")') trim(this%name),trim(this%src%name),trim(this%dst%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
      
   end subroutine initialize
   
   
end module coupler_class
