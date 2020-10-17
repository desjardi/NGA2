!> Definition of a partitioned grid class in NGA.
!> @todo Add other parallelization strategies
module pgrid_class
   use precision,   only: WP
   use sgrid_class, only: sgrid
   use string,      only: str_medium
   use mpi_f08,     only: MPI_COMM,MPI_GROUP,MPI_Datatype
   implicit none
   private
   
   !> Default parallelization strategy
   character(len=str_medium), parameter :: defstrat='fewest_dir'
   integer, parameter :: defmincell=4
   
   ! Expose type/constructor/methods
   public :: pgrid
   
   !> Partitioned grid type
   type, extends(sgrid) :: pgrid
      ! Parallelization information
      type(MPI_Group) :: group              !< Grid group
      type(MPI_Comm) :: comm                !< Grid communicator
      type(MPI_Datatype) :: view            !< Local to global array mapping info - real(WP)
      integer :: nproc                      !< Number of processors
      integer :: rank                       !< Processor grid rank
      logical :: amRoot                     !< Am I grid root?
      logical :: amIn                       !< Am I working in this grid?
      integer :: npx,npy,npz                !< Number of processors per direction
      integer :: iproc,jproc,kproc          !< Coordinates location of processor
      ! Local grid size
      integer :: nx_,ny_,nz_                !< Local grid size in x/y/z
      ! Local grid size with overlap
      integer :: nxo_,nyo_,nzo_             !< Local grid size in x/y/z with overlap
      ! Local index bounds
      integer :: imin_,imax_                !< Domain-decomposed index bounds in x
      integer :: jmin_,jmax_                !< Domain-decomposed index bounds in y
      integer :: kmin_,kmax_                !< Domain-decomposed index bounds in z
      ! Local index bounds with overlap
      integer :: imino_,imaxo_              !< Domain-decomposed index bounds in x with overlap
      integer :: jmino_,jmaxo_              !< Domain-decomposed index bounds in y with overlap
      integer :: kmino_,kmaxo_              !< Domain-decomposed index bounds in z with overlap
      
      !> Communication buffers
      real(WP), dimension(:,:,:), allocatable, private :: syncbuf_x1,syncbuf_x2
      real(WP), dimension(:,:,:), allocatable, private :: syncbuf_y1,syncbuf_y2
      real(WP), dimension(:,:,:), allocatable, private :: syncbuf_z1,syncbuf_z2
      
   contains
      procedure, private :: init_mpi=>pgrid_init_mpi           !< Prepare the MPI environment for the pgrid
      procedure, private :: domain_decomp=>pgrid_domain_decomp !< Perform the domain decomposition
      
      procedure :: allprint=>pgrid_allprint                    !< Output grid to screen - blocking and requires all procs...
      procedure :: print   =>pgrid_print                       !< Output grid to screen
      procedure :: log     =>pgrid_log                         !< Output grid info to log
      procedure :: sync    =>pgrid_sync,pgrid_sync_no          !< Commmunicate inner and periodic boundaries
   end type pgrid
   
   
   !> Declare partitioned grid constructor
   interface pgrid
      procedure construct_pgrid_from_sgrid
      procedure construct_pgrid_from_file
   end interface pgrid
   
contains
   
   
   !> Partitioned grid constructor from file
   function construct_pgrid_from_file(no,file,grp,decomp) result(self)
      use string,   only: lowercase
      use monitor,  only: die
      use param,    only: verbose
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      implicit none
      
      type(pgrid) :: self                               !< Parallel grid
      integer, intent(in) :: no                         !< Overlap size
      character(len=*), intent(in) :: file              !< Grid file
      type(MPI_Group), intent(in) :: grp                !< Processor group for parallelization
      integer, dimension(3), intent(in) :: decomp       !< Desired decomposition
      integer :: ierr
      character(len=str_medium) :: simu_name
      integer :: nx,ny,nz,coord
      real(WP), dimension(:), allocatable :: x
      real(WP), dimension(:), allocatable :: y
      real(WP), dimension(:), allocatable :: z
      logical :: xper,yper,zper
      
      ! Initialize MPI environment
      self%group=grp; call self%init_mpi
      
      ! Nothing more to do if the processor is not inside
      if (.not.self%amIn) return
      
      ! Root process can now read in the grid
      if (self%amRoot) then
         self%sgrid=sgrid(no,file)
         simu_name=self%name
         coord=self%coordsys
         xper=self%xper
         yper=self%yper
         zper=self%zper
         nx=self%nx
         ny=self%ny
         nz=self%nz
      end if
      
      ! Broadcast it to our group
      call MPI_BCAST(simu_name,len(simu_name),MPI_CHARACTER,0,self%comm,ierr)
      call MPI_BCAST(coord    ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      call MPI_BCAST(xper     ,1             ,MPI_LOGICAL  ,0,self%comm,ierr)
      call MPI_BCAST(yper     ,1             ,MPI_LOGICAL  ,0,self%comm,ierr)
      call MPI_BCAST(zper     ,1             ,MPI_LOGICAL  ,0,self%comm,ierr)
      call MPI_BCAST(nx       ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      call MPI_BCAST(ny       ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      call MPI_BCAST(nz       ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      
      ! Allocate x/y/z, fill it, and bc
      allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1))
      if (self%amRoot) then
         x(1:nx+1)=self%x(self%imin:self%imax+1)
         y(1:ny+1)=self%y(self%jmin:self%jmax+1)
         z(1:nz+1)=self%z(self%kmin:self%kmax+1)
      end if
      call MPI_BCAST(x,nx+1,MPI_REAL_WP,0,self%comm,ierr)
      call MPI_BCAST(y,ny+1,MPI_REAL_WP,0,self%comm,ierr)
      call MPI_BCAST(z,nz+1,MPI_REAL_WP,0,self%comm,ierr)
      
      ! Finish creating the sgrid
      if (.not.self%amRoot) self%sgrid=sgrid(coord,no,x,y,z,xper,yper,zper,trim(adjustl(simu_name)))
      
      ! Deallocate
      deallocate(x,y,z)
      
      ! Perform actual domain decomposition of grid
      call self%domain_decomp(decomp)
      
      ! If verbose run, log and or print grid
      if (verbose.gt.0) call self%log
      if (verbose.gt.1) call self%print
      if (verbose.gt.2) call self%allprint
      
   end function construct_pgrid_from_file
   
   
   !> Partitioned grid constructor from sgrid
   function construct_pgrid_from_sgrid(grid,grp,decomp) result(self)
      use string,  only: lowercase
      use monitor, only: die
      use param,   only: verbose
      implicit none
      include 'mpif.h'
      
      type(pgrid) :: self                               !< Parallel grid
      
      type(sgrid), intent(in) :: grid                   !< Base grid
      type(MPI_Group), intent(in) :: grp                !< Processor group for parallelization
      integer, dimension(3) :: decomp                   !< Requested domain decomposition
      
      ! Initialize MPI environment
      self%group=grp; call self%init_mpi
      
      ! Nothing more to do if the processor is not inside
      if (.not.self%amIn) return
      
      ! Assign base grid data
      self%sgrid=grid
      
      ! Perform actual domain decomposition of grid
      call self%domain_decomp(decomp)
      
      ! If verbose run, log and or print grid
      if (verbose.gt.0) call self%log
      if (verbose.gt.1) call self%print
      if (verbose.gt.2) call self%allprint
      
   end function construct_pgrid_from_sgrid
   
   
   !> Prepares the MPI environment for the pgrid
   subroutine pgrid_init_mpi(self)
      use parallel, only: comm
      use monitor , only: die
      use mpi_f08
      implicit none
      class(pgrid) :: self
      integer :: ierr
      ! Get group size, and intracommunicator
      call MPI_GROUP_SIZE(self%group,self%nproc,ierr)
      if (self%nproc.eq.0) call die('[pgrid constructor] Non-empty group is required!')
      call MPI_COMM_CREATE(comm,self%group,self%comm,ierr)
      ! Get rank and amIn info
      call MPI_GROUP_RANK(self%group,self%rank,ierr)
      self%amIn=(self%rank.ne.MPI_UNDEFINED)
      ! Handle processors that are not part of the group
      if (.not.self%amIn) then
         self%amRoot=.false.
         self%iproc=0; self%nx_=0; self%imin_=0; self%imax_=0; self%nxo_=0; self%imino_=0; self%imaxo_=0
         self%jproc=0; self%ny_=0; self%jmin_=0; self%jmax_=0; self%nyo_=0; self%jmino_=0; self%jmaxo_=0
         self%kproc=0; self%nz_=0; self%kmin_=0; self%kmax_=0; self%nzo_=0; self%kmino_=0; self%kmaxo_=0
      end if
   end subroutine pgrid_init_mpi
   
   
   !> Prepares the domain decomposition of the pgrid
   subroutine pgrid_domain_decomp(self,decomp)
      use monitor , only: die
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      implicit none
      class(pgrid) :: self
      integer, dimension(3), intent(in) :: decomp
      integer :: ierr,q,r
      type(MPI_Comm) :: tmp_comm
      integer, parameter :: ndims=3
      logical, parameter :: reorder=.true.
      integer, dimension(3) :: coords
      integer, dimension(3) :: gsizes,lsizes,lstart
      
      ! Store decomposition on the grid
      self%npx=decomp(1); self%npy=decomp(2); self%npz=decomp(3)
      
      ! Perform sanity check of the decomposition
      if (self%npx*self%npy*self%npz.ne.self%nproc) call die('[pgrid constructor] Parallel decomposition is improper')
      
      ! Give cartesian layout to intracommunicator
      call MPI_CART_CREATE(self%comm,ndims,[self%npx,self%npy,self%npz],[self%xper,self%yper,self%zper],reorder,tmp_comm,ierr); self%comm=tmp_comm
      call MPI_COMM_RANK  (self%comm,self%rank,ierr)
      call MPI_CART_COORDS(self%comm,self%rank,ndims,coords,ierr)
      self%iproc=coords(1)+1; self%jproc=coords(2)+1; self%kproc=coords(3)+1
      self%amRoot=(self%rank.eq.0)
      
      ! Perform decomposition in x
      q=self%nx/self%npx; r=mod(self%nx,self%npx)
      if (self%iproc.le.r) then
         self%nx_  =q+1
         self%imin_=self%imin+(self%iproc-1)*self%nx_
      else
         self%nx_  =q
         self%imin_=self%imin+(self%iproc-1)*self%nx_+r
      end if
      self%imax_ =self%imin_+self%nx_-1
      self%nxo_  =self%nx_+2*self%no
      self%imino_=self%imin_-self%no
      self%imaxo_=self%imax_+self%no
      
      ! Perform decomposition in y
      q=self%ny/self%npy; r=mod(self%ny,self%npy)
      if (self%jproc.le.r) then
         self%ny_  =q+1
         self%jmin_=self%jmin+(self%jproc-1)*self%ny_
      else
         self%ny_  =q
         self%jmin_=self%jmin+(self%jproc-1)*self%ny_+r
      end if
      self%jmax_ =self%jmin_+self%ny_-1
      self%nyo_  =self%ny_+2*self%no
      self%jmino_=self%jmin_-self%no
      self%jmaxo_=self%jmax_+self%no
      
      ! Perform decomposition in z
      q=self%nz/self%npz; r=mod(self%nz,self%npz)
      if (self%kproc.le.r) then
         self%nz_  =q+1
         self%kmin_=self%kmin+(self%kproc-1)*self%nz_
      else
         self%nz_  =q
         self%kmin_=self%kmin+(self%kproc-1)*self%nz_+r
      end if
      self%kmax_ =self%kmin_+self%nz_-1
      self%nzo_  =self%nz_+2*self%no
      self%kmino_=self%kmin_-self%no
      self%kmaxo_=self%kmax_+self%no
      
      ! We also need to prepare communication buffers
      allocate(self%syncbuf_x1(self%no,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_))
      allocate(self%syncbuf_x2(self%no,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_))
      allocate(self%syncbuf_y1(self%imino_:self%imaxo_,self%no,self%kmino_:self%kmaxo_))
      allocate(self%syncbuf_y2(self%imino_:self%imaxo_,self%no,self%kmino_:self%kmaxo_))
      allocate(self%syncbuf_z1(self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%no))
      allocate(self%syncbuf_z2(self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%no))
      
      ! Finally, we need to define a proper MPI-I/O view
      gsizes=[self%nx ,self%ny ,self%nz ]
      lsizes=[self%nx_,self%ny_,self%nz_]
      lstart=[self%imin_-self%imin,self%jmin_-self%jmin,self%kmin_-self%kmin]
      call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_REAL_WP,self%view,ierr)
      call MPI_TYPE_COMMIT(self%view,ierr)
      
   end subroutine pgrid_domain_decomp
   
   
   !> Print out partitioned grid info to screen
   !> This is a slow and blocking routine for debugging only!
   subroutine pgrid_allprint(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use parallel, only: rank
      implicit none
      class(pgrid), intent(in) :: this
      integer :: i,ierr
      ! Return all non-involved procs
      if (.not.this%amIn) return
      ! Root writes header
      call pgrid_print(this)
      ! Each proc provides info on his involvement in the grid
      do i=0,this%nproc-1
         ! Block for clean output
         call MPI_BARRIER(this%comm,ierr)
         ! Output info
         if (this%rank.eq.i) then
            write(output_unit,'(" --> Rank ",i0,"(",i0,") -> [",i0,",",i0,",",i0,"] owns [",i0,",",i0,"]x[",i0,",",i0,"]x[",i0,",",i0,"]")') &
            this%rank,rank,this%iproc,this%jproc,this%kproc,this%imin_,this%imax_,this%jmin_,this%jmax_,this%kmin_,this%kmax_
         end if
      end do
   end subroutine pgrid_allprint
   
   
   !> Cheap print of partitioned grid info to screen
   subroutine pgrid_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use sgrid_class, only: cartesian,cylindrical,spherical
      implicit none
      class(pgrid), intent(in) :: this
      if (this%amRoot) then
         select case (this%coordsys)
         case (cartesian);   write(output_unit,'("Partitioned cartesian grid [",a,"] on ",i0," processes")') trim(this%name),this%nproc
         case (cylindrical); write(output_unit,'("Partitioned cylindrical grid [",a,"]")') trim(this%name)
         case (spherical);   write(output_unit,'("Partitioned spherical grid [",a,"]")') trim(this%name)
         end select
         write(output_unit,'(" >   overlap = ",i0)') this%no
         write(output_unit,'(" >      size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
         write(output_unit,'(" >    decomp = ",i0,"x",i0,"x",i0)') this%npx,this%npy,this%npz
         write(output_unit,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
         write(output_unit,'(" >   uniform = ",l1,"x",l1,"x",l1)') this%uniform_x,this%uniform_y,this%uniform_z
         write(output_unit,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper
         write(output_unit,'(" >    volume = ",es12.5)') this%vol_total
      end if
   end subroutine pgrid_print
   
   
   !> Cheap print of partitioned grid info to log
   subroutine pgrid_log(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use monitor,     only: log
      use string,      only: str_long
      use sgrid_class, only: cartesian,cylindrical,spherical
      implicit none
      class(pgrid), intent(in) :: this
      character(len=str_long) :: message
      if (this%amRoot) then
         select case (this%coordsys)
         case (cartesian);   write(message,'("Partitioned cartesian grid [",a,"] on ",i0," processes")') trim(this%name),this%nproc; call log(message)
         case (cylindrical); write(message,'("Partitioned cylindrical grid [",a,"]")') trim(this%name); call log(message)
         case (spherical);   write(message,'("Partitioned spherical grid [",a,"]")') trim(this%name); call log(message)
         end select
         write(message,'(" >   overlap = ",i0)') this%no; call log(message)
         write(message,'(" >      size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz; call log(message)
         write(message,'(" >    decomp = ",i0,"x",i0,"x",i0)') this%npx,this%npy,this%npz; call log(message)
         write(message,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1); call log(message)
         write(message,'(" >   uniform = ",l1,"x",l1,"x",l1)') this%uniform_x,this%uniform_y,this%uniform_z; call log(message)
         write(message,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper; call log(message)
         write(message,'(" >    volume = ",es12.5)') this%vol_total; call log(message)
      end if
   end subroutine pgrid_log
   
   
   !> Synchronization of overlap cells - uses full no
   !> This routine assumes that the default overlap size is used
   !> It allows the use of pre-allocated buffers for speed
   subroutine pgrid_sync(this,A)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pgrid) :: this
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: A !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr,isize,i,j,k
      
      ! Work in x - is it 2D or 3D?
      if (this%nx.eq.1) then
         ! Direct copy if 2D
         do i=this%imax_+1,this%imaxo_
            A(i,:,:)=A(this%imin_,:,:)
         end do
         do i=this%imino_,this%imin_-1
            A(i,:,:)=A(this%imin_,:,:)
         end do
      else
         isize=(this%no)*(this%nyo_)*(this%nzo_)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,0,-1,isrc,idst,ierr)
         this%syncbuf_x1=A(this%imin_:this%imin_+this%no-1,:,:)
         call MPI_SENDRECV(this%syncbuf_x1,isize,MPI_REAL_WP,idst,0,this%syncbuf_x2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imax_+1:this%imaxo_,:,:)=this%syncbuf_x2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
         this%syncbuf_x1=A(this%imax_-this%no+1:this%imax_,:,:)
         call MPI_SENDRECV(this%syncbuf_x1,isize,MPI_REAL_WP,idst,0,this%syncbuf_x2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imino_:this%imin_-1,:,:)=this%syncbuf_x2
      end if
      
      ! Work in y - is it 2D or 3D?
      if (this%ny.eq.1) then
         ! Direct copy if 2D
         do j=this%jmax_+1,this%jmaxo_
            A(:,j,:)=A(:,this%jmin_,:)
         end do
         do j=this%jmino_,this%jmin_-1
            A(:,j,:)=A(:,this%jmin_,:)
         end do
      else
         isize=(this%nxo_)*(this%no)*(this%nzo_)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,1,-1,isrc,idst,ierr)
         this%syncbuf_y1=A(:,this%jmin_:this%jmin_+this%no-1,:)
         call MPI_SENDRECV(this%syncbuf_y1,isize,MPI_REAL_WP,idst,0,this%syncbuf_y2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmax_+1:this%jmaxo_,:)=this%syncbuf_y2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
         this%syncbuf_y1=A(:,this%jmax_-this%no+1:this%jmax_,:)
         call MPI_SENDRECV(this%syncbuf_y1,isize,MPI_REAL_WP,idst,0,this%syncbuf_y2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmino_:this%jmin_-1,:)=this%syncbuf_y2
      end if
      
      ! Work in z - is it 2D or 3D?
      if (this%nz.eq.1) then
         ! Direct copy if 2D
         do k=this%kmax_+1,this%kmaxo_
            A(:,:,k)=A(:,:,this%kmin_)
         end do
         do j=this%kmino_,this%kmin_-1
            A(:,:,k)=A(:,:,this%kmin_)
         end do
      else
         isize=(this%nxo_)*(this%nyo_)*(this%no)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
         this%syncbuf_z1=A(:,:,this%kmin_:this%kmin_+this%no-1)
         call MPI_SENDRECV(this%syncbuf_z1,isize,MPI_REAL_WP,idst,0,this%syncbuf_z2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmax_+1:this%kmaxo_)=this%syncbuf_z2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
         this%syncbuf_z1=A(:,:,this%kmax_-this%no+1:this%kmax_)
         call MPI_SENDRECV(this%syncbuf_z1,isize,MPI_REAL_WP,idst,0,this%syncbuf_z2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmino_:this%kmin_-1)=this%syncbuf_z2
      end if
      
   end subroutine pgrid_sync
   
   
   !> Synchronization of overlap cells
   !> This version is capable of handling any overlap size
   subroutine pgrid_sync_no(this,A,no)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pgrid) :: this
      integer, intent(in) :: no
      real(WP), dimension(this%imin_-no:,this%jmin_-no:,this%kmin_-no:), intent(inout) :: A !< Needs to be (imin_-no:imax_+no,jmin_-no:jmax_+no,kmin_-no:kmax_+no)
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr,isize,i,j,k
      real(WP), dimension(:,:,:), allocatable :: buf1,buf2
      
      ! Work in x - is it 2D or 3D?
      if (this%nx.eq.1) then
         ! Direct copy if 2D
         do i=this%imax_+1,this%imax_+no
            A(i,:,:)=A(this%imin_,:,:)
         end do
         do i=this%imin_-no,this%imin_-1
            A(i,:,:)=A(this%imin_,:,:)
         end do
      else
         isize=(no)*(this%ny_+2*no)*(this%nz_+2*no)
         allocate(buf1(no,this%ny_+2*no,this%nz_+2*no))
         allocate(buf2(no,this%ny_+2*no,this%nz_+2*no))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,0,-1,isrc,idst,ierr)
         buf1=A(this%imin_:this%imin_+no-1,:,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imax_+1:this%imax_+no,:,:)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
         buf1=A(this%imax_-no+1:this%imax_,:,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imin_-no:this%imin_-1,:,:)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
      ! Work in y - is it 2D or 3D?
      if (this%ny.eq.1) then
         ! Direct copy if 2D
         do j=this%jmax_+1,this%jmax_+no
            A(:,j,:)=A(:,this%jmin_,:)
         end do
         do j=this%jmin_-no,this%jmin_-1
            A(:,j,:)=A(:,this%jmin_,:)
         end do
      else
         isize=(this%nx_+2*no)*(no)*(this%nz_+2*no)
         allocate(buf1(this%nx_+2*no,no,this%nz_+2*no))
         allocate(buf2(this%nx_+2*no,no,this%nz_+2*no))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,1,-1,isrc,idst,ierr)
         buf1=A(:,this%jmin_:this%jmin_+no-1,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmax_+1:this%jmax_+no,:)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
         buf1=A(:,this%jmax_-no+1:this%jmax_,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmin_-no:this%jmin_-1,:)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
      ! Work in z - is it 2D or 3D?
      if (this%nz.eq.1) then
         ! Direct copy if 2D
         do k=this%kmax_+1,this%kmax_+no
            A(:,:,k)=A(:,:,this%kmin_)
         end do
         do j=this%kmin_-no,this%kmin_-1
            A(:,:,k)=A(:,:,this%kmin_)
         end do
      else
         isize=(this%nx_+2*no)*(this%ny_+2*no)*(no)
         allocate(buf1(this%nx_+2*no,this%ny_+2*no,no))
         allocate(buf2(this%nx_+2*no,this%ny_+2*no,no))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
         buf1=A(:,:,this%kmin_:this%kmin_+no-1)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmax_+1:this%kmax_+no)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
         buf1=A(:,:,this%kmax_-no+1:this%kmax_)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmin_-no:this%kmin_-1)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
   end subroutine pgrid_sync_no
   
   
end module pgrid_class
