!> Definition of a partitioned grid class in NGA.
!> All processors that call the constructor need to be in the group passed by the user.
module pgrid_class
   use precision,   only: WP
   use sgrid_class, only: sgrid
   use string,      only: str_medium
   use mpi_f08
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
      type(MPI_Group)    :: group                             !< Group of processors working on the pgrid
      type(MPI_Comm)     :: comm  =MPI_COMM_NULL              !< Communicator for our group
      type(MPI_Datatype) :: view  =MPI_DATATYPE_NULL          !< Local to global array mapping info - real(WP)
      type(MPI_Datatype) :: Iview =MPI_DATATYPE_NULL          !< Local to global array mapping info - integer
      type(MPI_Datatype) :: SPview=MPI_DATATYPE_NULL          !< Local to global array mapping info - real(SP)
      integer :: nproc                      !< Number of processors
      integer :: rank                       !< Processor grid rank
      logical :: amRoot                     !< Am I grid root?
      integer :: npx                        !< Number of processors in x
      integer :: npy                        !< Number of processors in y
      integer :: npz                        !< Number of processors in z
      integer :: iproc=0                    !< Coordinate location of processor in x
      integer :: jproc=0                    !< Coordinate location of processor in y
      integer :: kproc=0                    !< Coordinate location of processor in z
      type(MPI_Comm) :: xcomm=MPI_COMM_NULL !< 1D x-communicator
      type(MPI_Comm) :: ycomm=MPI_COMM_NULL !< 1D y-communicator
      type(MPI_Comm) :: zcomm=MPI_COMM_NULL !< 1D z-communicator
      integer :: xrank=MPI_UNDEFINED        !< 1D rank for xcomm
      integer :: yrank=MPI_UNDEFINED        !< 1D rank for ycomm
      integer :: zrank=MPI_UNDEFINED        !< 1D rank for zcomm
      type(MPI_Comm) :: xycomm=MPI_COMM_NULL!< 2D xy-communicator
      type(MPI_Comm) :: yzcomm=MPI_COMM_NULL!< 2D yz-communicator
      type(MPI_Comm) :: zxcomm=MPI_COMM_NULL!< 2D zx-communicator
      integer :: xyrank=MPI_UNDEFINED       !< 2D rank for xycomm
      integer :: yzrank=MPI_UNDEFINED       !< 2D rank for yzcomm
      integer :: zxrank=MPI_UNDEFINED       !< 2D rank for zxcomm
      
      ! Local grid size
      integer :: nx_=0                      !< Local grid size in x
      integer :: ny_=0                      !< Local grid size in y
      integer :: nz_=0                      !< Local grid size in z
      ! Local grid size with overlap
      integer :: nxo_=0                     !< Local grid size in x with overlap
      integer :: nyo_=0                     !< Local grid size in y with overlap
      integer :: nzo_=0                     !< Local grid size in z with overlap
      ! Local index bounds
      integer :: imin_=0                    !< Domain-decomposed index lower bound in x
      integer :: imax_=0                    !< Domain-decomposed index upper bound in x
      integer :: jmin_=0                    !< Domain-decomposed index lower bound in y
      integer :: jmax_=0                    !< Domain-decomposed index upper bound in y
      integer :: kmin_=0                    !< Domain-decomposed index lower bound in z
      integer :: kmax_=0                    !< Domain-decomposed index upper bound in z
      ! Local index bounds with overlap
      integer :: imino_=0                   !< Domain-decomposed index lower bound in x with overlap
      integer :: imaxo_=0                   !< Domain-decomposed index upper bound in x with overlap
      integer :: jmino_=0                   !< Domain-decomposed index lower bound in y with overlap
      integer :: jmaxo_=0                   !< Domain-decomposed index upper bound in y with overlap
      integer :: kmino_=0                   !< Domain-decomposed index lower bound in z with overlap
      integer :: kmaxo_=0                   !< Domain-decomposed index upper bound in z with overlap
      
      ! Index to processor coordinate map
      integer, dimension(:), allocatable :: xcoord                              !< Conversion from grid index to processor coordinate in x
      integer, dimension(:), allocatable :: ycoord                              !< Conversion from grid index to processor coordinate in y
      integer, dimension(:), allocatable :: zcoord                              !< Conversion from grid index to processor coordinate in z
      
      !> Communication buffers
      real(WP), dimension(:,:,:), allocatable, private ::  syncbuf_x1, syncbuf_x2
      real(WP), dimension(:,:,:), allocatable, private ::  syncbuf_y1, syncbuf_y2
      real(WP), dimension(:,:,:), allocatable, private ::  syncbuf_z1, syncbuf_z2
      integer , dimension(:,:,:), allocatable, private :: isyncbuf_x1,isyncbuf_x2
      integer , dimension(:,:,:), allocatable, private :: isyncbuf_y1,isyncbuf_y2
      integer , dimension(:,:,:), allocatable, private :: isyncbuf_z1,isyncbuf_z2
      
   contains
      procedure, private :: init_mpi=>pgrid_init_mpi                            !< Prepare the MPI environment for the pgrid
      procedure, private :: domain_decomp=>pgrid_domain_decomp                  !< Perform the domain decomposition
      procedure :: allprint=>pgrid_allprint                                     !< Output grid to screen - blocking and requires all procs...
      procedure :: print   =>pgrid_print                                        !< Output grid to screen
      procedure :: log     =>pgrid_log                                          !< Output grid info to log
      generic :: sync=>pgrid_rsync,pgrid_rsync_array,pgrid_rsync_tensor,pgrid_rsync_no,pgrid_isync,pgrid_isync_no    !< Commmunicate inner and periodic boundaries - generic
      procedure, private :: pgrid_isync,pgrid_isync_no                          !< Commmunicate inner and periodic boundaries for integer
      procedure, private :: pgrid_rsync,pgrid_rsync_no                          !< Commmunicate inner and periodic boundaries for real(WP)
      procedure, private :: pgrid_rsync_array                                   !< Commmunicate inner and periodic boundaries for arrays of real(WP) of the form (:,i,j,k)
	  procedure, private :: pgrid_rsync_tensor                                  !< Commmunicate inner and periodic boundaries for tensors of real(WP) of the form (:,:,i,j,k)
      generic :: syncsum=>pgrid_rsyncsum                                        !< Summation across inner and periodic boundaries - generic
      procedure, private :: pgrid_rsyncsum                                      !< Summation inner and periodic boundaries for real(WP)
      procedure :: get_rank                                                     !< Function that returns rank of processor that contains provided indices
      procedure :: get_ijk_local                                                !< Function that returns closest mesh indices to a provided position - local to processor subdomain
      procedure :: get_ijk_global                                               !< Function that returns closest mesh indices to a provided position - global over full pgrid
      procedure :: get_ijk_from_lexico,get_lexico_from_ijk                      !< Functions that convert a lexicographic index to (i,j,k) and vice-versa
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
      use messager, only: die
      use param,    only: verbose
      use parallel, only: MPI_REAL_WP
      implicit none
      type(pgrid) :: self                               !< Parallel grid
      integer, intent(in) :: no                         !< Overlap size
      character(len=*), intent(in) :: file              !< Grid file
      type(MPI_Group), intent(in) :: grp                !< MPI group
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
      
      ! Check and store decomposition
      if (product(decomp).ne.self%nproc) call die('[pgrid constructor] Parallel decomposition is improper')
      self%npx=decomp(1); self%npy=decomp(2); self%npz=decomp(3)
      
      ! Perform actual domain decomposition of grid
      call self%domain_decomp()
      
      ! If verbose run, log and or print grid
      if (verbose.gt.0) call self%log
      if (verbose.gt.1) call self%print
      if (verbose.gt.2) call self%allprint
      
   end function construct_pgrid_from_file
   
   
   !> Partitioned grid constructor from sgrid
   function construct_pgrid_from_sgrid(grid,grp,decomp) result(self)
      use string,   only: lowercase
      use messager, only: die
      use param,    only: verbose
      implicit none
      type(pgrid) :: self                               !< Parallel grid
      type(sgrid), intent(in) :: grid                   !< Base grid
      type(MPI_Group), intent(in) :: grp                !< MPI group
      integer, dimension(3), intent(in) :: decomp       !< Requested domain decomposition
      
      ! Initialize MPI environment
      self%group=grp; call self%init_mpi
      
      ! Copy base grid data
      self%sgrid=grid
      
      ! Check and store decomposition
      if (product(decomp).ne.self%nproc) call die('[pgrid constructor] Parallel decomposition is improper')
      self%npx=decomp(1); self%npy=decomp(2); self%npz=decomp(3)
      
      ! Perform actual domain decomposition of grid
      call self%domain_decomp()
      
      ! If verbose run, log and or print grid
      if (verbose.gt.0) call self%log
      if (verbose.gt.1) call self%print
      if (verbose.gt.2) call self%allprint
      
   end function construct_pgrid_from_sgrid
   
   
   !> Prepares the MPI environment for the pgrid
   subroutine pgrid_init_mpi(self)
      use parallel, only: comm
      use messager, only: die
      implicit none
      class(pgrid), intent(inout) :: self
      integer :: ierr
      ! Get group size
      call MPI_GROUP_SIZE(self%group,self%nproc,ierr)
      if (self%nproc.eq.0) call die('[pgrid constructor] A non-empty group is required')
      ! Get rank
      call MPI_GROUP_RANK(self%group,self%rank,ierr)
      ! Get a root process
      self%amRoot=(self%rank.eq.0)
      ! Test if a processors was not in the group
      if (self%rank.eq.MPI_UNDEFINED) call die('[pgrid constructor] All processors that call the constructor need to be in the group')
      ! Create intracommunicator for the group
      call MPI_COMM_CREATE_GROUP(comm,self%group,0,self%comm,ierr)
   end subroutine pgrid_init_mpi
   
   
   !> Prepares the domain decomposition of the pgrid
   subroutine pgrid_domain_decomp(self)
      use messager, only: die
      use parallel, only: MPI_REAL_WP,MPI_REAL_SP
      implicit none
      class(pgrid), intent(inout) :: self
      integer :: ierr,q,r
      type(MPI_Comm) :: tmp_comm
      integer, parameter :: ndims=3
      logical, parameter :: reorder=.true.
      logical, dimension(3) :: dir
      integer, dimension(3) :: coords
      integer, dimension(3) :: gsizes,lsizes,lstart
      
      ! Give cartesian layout to intracommunicator
      call MPI_CART_CREATE(self%comm,ndims,[self%npx,self%npy,self%npz],[self%xper,self%yper,self%zper],reorder,tmp_comm,ierr); self%comm=tmp_comm
      call MPI_COMM_RANK  (self%comm,self%rank,ierr)
      call MPI_CART_COORDS(self%comm,self%rank,ndims,coords,ierr)
      self%iproc=coords(1)+1; self%jproc=coords(2)+1; self%kproc=coords(3)+1
      self%amRoot=(self%rank.eq.0)
      
      ! Create 1D communicators
      dir=[.true.,.false.,.false.]
      call MPI_CART_SUB(self%comm,dir,self%xcomm,ierr)
      call MPI_COMM_RANK(self%xcomm,self%xrank,ierr)
      dir=[.false.,.true.,.false.]
      call MPI_CART_SUB(self%comm,dir,self%ycomm,ierr)
      call MPI_COMM_RANK(self%ycomm,self%yrank,ierr)
      dir=[.false.,.false.,.true.]
      call MPI_CART_SUB(self%comm,dir,self%zcomm,ierr)
      call MPI_COMM_RANK(self%zcomm,self%zrank,ierr)
      
      ! Create 2D communicators
      dir=[.true.,.true.,.false.]
      call MPI_CART_SUB(self%comm,dir,self%xycomm,ierr)
      call MPI_COMM_RANK(self%xycomm,self%xyrank,ierr)
      dir=[.false.,.true.,.true.]
      call MPI_CART_SUB(self%comm,dir,self%yzcomm,ierr)
      call MPI_COMM_RANK(self%yzcomm,self%yzrank,ierr)
      dir=[.true.,.false.,.true.]
      call MPI_CART_SUB(self%comm,dir,self%zxcomm,ierr)
      call MPI_COMM_RANK(self%zxcomm,self%zxrank,ierr)
      
      ! Perform decomposition in x
      q=self%nx/self%npx; r=mod(self%nx,self%npx)
      self%imin_ =self%imin+ coords(1)   *q+min(coords(1)  ,r)
      self%imax_ =self%imin+(coords(1)+1)*q+min(coords(1)+1,r)-1
      self%nx_   =self%imax_-self%imin_+1
      self%nxo_  =self%nx_+2*self%no
      self%imino_=self%imin_-self%no
      self%imaxo_=self%imax_+self%no
      
      ! Perform decomposition in y
      q=self%ny/self%npy; r=mod(self%ny,self%npy)
      self%jmin_ =self%jmin+ coords(2)   *q+min(coords(2)  ,r)
      self%jmax_ =self%jmin+(coords(2)+1)*q+min(coords(2)+1,r)-1
      self%ny_   =self%jmax_-self%jmin_+1
      self%nyo_  =self%ny_+2*self%no
      self%jmino_=self%jmin_-self%no
      self%jmaxo_=self%jmax_+self%no
      
      ! Perform decomposition in z
      q=self%nz/self%npz; r=mod(self%nz,self%npz)
      self%kmin_ =self%kmin+ coords(3)   *q+min(coords(3)  ,r)
      self%kmax_ =self%kmin+(coords(3)+1)*q+min(coords(3)+1,r)-1
      self%nz_   =self%kmax_-self%kmin_+1
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
      allocate(self%isyncbuf_x1(self%no,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_))
      allocate(self%isyncbuf_x2(self%no,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_))
      allocate(self%isyncbuf_y1(self%imino_:self%imaxo_,self%no,self%kmino_:self%kmaxo_))
      allocate(self%isyncbuf_y2(self%imino_:self%imaxo_,self%no,self%kmino_:self%kmaxo_))
      allocate(self%isyncbuf_z1(self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%no))
      allocate(self%isyncbuf_z2(self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%no))
      
      ! We need to define a proper MPI-I/O view
      gsizes=[self%nx ,self%ny ,self%nz ]
      lsizes=[self%nx_,self%ny_,self%nz_]
      lstart=[self%imin_-self%imin,self%jmin_-self%jmin,self%kmin_-self%kmin]
      call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_REAL_WP,self%view,ierr)
      call MPI_TYPE_COMMIT(self%view,ierr)
      call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_REAL_SP,self%SPview,ierr)
      call MPI_TYPE_COMMIT(self%SPview,ierr)
      call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_INTEGER,self%Iview,ierr)
      call MPI_TYPE_COMMIT(self%Iview,ierr)
      
      ! Finally, create x/y/zcoord array for rapid finding of processor cartesian coordinates
      allocate(self%xcoord(self%imino:self%imaxo)); self%xcoord=0
      allocate(self%ycoord(self%jmino:self%jmaxo)); self%ycoord=0
      allocate(self%zcoord(self%kmino:self%kmaxo)); self%zcoord=0
      find_coords: block
         integer :: i,j,k
         q=self%nx/self%npx; r=mod(self%nx,self%npx)
         do i=self%imino,self%imaxo
            do while (i.ge.self%imin+(self%xcoord(i)+1)*q+min(self%xcoord(i)+1,r).and.self%xcoord(i).lt.self%npx-1)
               self%xcoord(i)=self%xcoord(i)+1
            end do
         end do
         q=self%ny/self%npy; r=mod(self%ny,self%npy)
         do j=self%jmino,self%jmaxo
            do while (j.ge.self%jmin+(self%ycoord(j)+1)*q+min(self%ycoord(j)+1,r).and.self%ycoord(j).lt.self%npy-1)
               self%ycoord(j)=self%ycoord(j)+1
            end do
         end do
         q=self%nz/self%npz; r=mod(self%nz,self%npz)
         do k=self%kmino,self%kmaxo
            do while (k.ge.self%kmin+(self%zcoord(k)+1)*q+min(self%zcoord(k)+1,r).and.self%zcoord(k).lt.self%npz-1)
               self%zcoord(k)=self%zcoord(k)+1
            end do
         end do
      end block find_coords
      
   end subroutine pgrid_domain_decomp
   
   
   !> Print out partitioned grid info to screen
   !> This is a slow and blocking routine for debugging only!
   subroutine pgrid_allprint(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use parallel, only: rank
      implicit none
      class(pgrid), intent(in) :: this
      integer :: i,ierr
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
      use messager,    only: log
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
   subroutine pgrid_rsync(this,A)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pgrid), intent(inout) :: this
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
         do k=this%kmino_,this%kmin_-1
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
      
   end subroutine pgrid_rsync
   
   
   !> Synchronization of overlap cells - uses full no
   !> This routine assumes that the default overlap size is used
   !> It allows the use of pre-allocated buffers for speed
   subroutine pgrid_isync(this,A)
      implicit none
      class(pgrid), intent(inout) :: this
      integer, dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: A !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
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
         this%isyncbuf_x1=A(this%imin_:this%imin_+this%no-1,:,:)
         call MPI_SENDRECV(this%isyncbuf_x1,isize,MPI_INTEGER,idst,0,this%isyncbuf_x2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imax_+1:this%imaxo_,:,:)=this%isyncbuf_x2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
         this%isyncbuf_x1=A(this%imax_-this%no+1:this%imax_,:,:)
         call MPI_SENDRECV(this%isyncbuf_x1,isize,MPI_INTEGER,idst,0,this%isyncbuf_x2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imino_:this%imin_-1,:,:)=this%isyncbuf_x2
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
         this%isyncbuf_y1=A(:,this%jmin_:this%jmin_+this%no-1,:)
         call MPI_SENDRECV(this%isyncbuf_y1,isize,MPI_INTEGER,idst,0,this%isyncbuf_y2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmax_+1:this%jmaxo_,:)=this%isyncbuf_y2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
         this%isyncbuf_y1=A(:,this%jmax_-this%no+1:this%jmax_,:)
         call MPI_SENDRECV(this%isyncbuf_y1,isize,MPI_INTEGER,idst,0,this%isyncbuf_y2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmino_:this%jmin_-1,:)=this%isyncbuf_y2
      end if
      
      ! Work in z - is it 2D or 3D?
      if (this%nz.eq.1) then
         ! Direct copy if 2D
         do k=this%kmax_+1,this%kmaxo_
            A(:,:,k)=A(:,:,this%kmin_)
         end do
         do k=this%kmino_,this%kmin_-1
            A(:,:,k)=A(:,:,this%kmin_)
         end do
      else
         isize=(this%nxo_)*(this%nyo_)*(this%no)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
         this%isyncbuf_z1=A(:,:,this%kmin_:this%kmin_+this%no-1)
         call MPI_SENDRECV(this%isyncbuf_z1,isize,MPI_INTEGER,idst,0,this%isyncbuf_z2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmax_+1:this%kmaxo_)=this%isyncbuf_z2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
         this%isyncbuf_z1=A(:,:,this%kmax_-this%no+1:this%kmax_)
         call MPI_SENDRECV(this%isyncbuf_z1,isize,MPI_INTEGER,idst,0,this%isyncbuf_z2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmino_:this%kmin_-1)=this%isyncbuf_z2
      end if
      
   end subroutine pgrid_isync
   
   
   !> Synchronization of overlap cells
   !> This version is capable of handling any overlap size
   subroutine pgrid_rsync_no(this,A,no)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pgrid), intent(in) :: this
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
         do k=this%kmin_-no,this%kmin_-1
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
      
   end subroutine pgrid_rsync_no
   
   
   !> Synchronization of overlap cells
   !> This version is capable of handling an array of the shape (:,i,j,k)
   subroutine pgrid_rsync_array(this,A)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pgrid), intent(in) :: this
      real(WP), dimension(1:,this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: A !< Needs to be (:,imin_-no:imax_+no,jmin_-no:jmax_+no,kmin_-no:kmax_+no)
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr,isize,i,j,k,dim
      real(WP), dimension(:,:,:,:), allocatable :: buf1,buf2
      
      ! Get first dimension
      dim=size(A,DIM=1)
      
      ! Work in x - is it 2D or 3D?
      if (this%nx.eq.1) then
         ! Direct copy if 2D
         do i=this%imax_+1,this%imaxo_
            A(:,i,:,:)=A(:,this%imin_,:,:)
         end do
         do i=this%imino_,this%imin_-1
            A(:,i,:,:)=A(:,this%imin_,:,:)
         end do
      else
         isize=dim*(this%no)*(this%nyo_)*(this%nzo_)
         allocate(buf1(dim,this%no,this%nyo_,this%nzo_))
         allocate(buf2(dim,this%no,this%nyo_,this%nzo_))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,0,-1,isrc,idst,ierr)
         buf1=A(:,this%imin_:this%imin_+this%no-1,:,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%imax_+1:this%imaxo_,:,:)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
         buf1=A(:,this%imax_-this%no+1:this%imax_,:,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%imino_:this%imin_-1,:,:)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
      ! Work in y - is it 2D or 3D?
      if (this%ny.eq.1) then
         ! Direct copy if 2D
         do j=this%jmax_+1,this%jmaxo_
            A(:,:,j,:)=A(:,:,this%jmin_,:)
         end do
         do j=this%jmino_,this%jmin_-1
            A(:,:,j,:)=A(:,:,this%jmin_,:)
         end do
      else
         isize=dim*(this%nxo_)*(this%no)*(this%nzo_)
         allocate(buf1(dim,this%nxo_,this%no,this%nzo_))
         allocate(buf2(dim,this%nxo_,this%no,this%nzo_))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,1,-1,isrc,idst,ierr)
         buf1=A(:,:,this%jmin_:this%jmin_+this%no-1,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%jmax_+1:this%jmaxo_,:)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
         buf1=A(:,:,this%jmax_-this%no+1:this%jmax_,:)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%jmino_:this%jmin_-1,:)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
      ! Work in z - is it 2D or 3D?
      if (this%nz.eq.1) then
         ! Direct copy if 2D
         do k=this%kmax_+1,this%kmaxo_
            A(:,:,:,k)=A(:,:,:,this%kmin_)
         end do
         do k=this%kmino_,this%kmin_-1
            A(:,:,:,k)=A(:,:,:,this%kmin_)
         end do
      else
         isize=dim*(this%nxo_)*(this%nyo_)*(this%no)
         allocate(buf1(dim,this%nxo_,this%nyo_,this%no))
         allocate(buf2(dim,this%nxo_,this%nyo_,this%no))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
         buf1=A(:,:,:,this%kmin_:this%kmin_+this%no-1)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,:,this%kmax_+1:this%kmaxo_)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
         buf1=A(:,:,:,this%kmax_-this%no+1:this%kmax_)
         call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,:,this%kmino_:this%kmin_-1)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
   end subroutine pgrid_rsync_array


   !> Synchronization of overlap cells
   !> This version is capable of handling a tensor of the shape (:,:,i,j,k)
   subroutine pgrid_rsync_tensor(this,A)
	  use parallel, only: MPI_REAL_WP
	  implicit none
	  class(pgrid), intent(in) :: this
	  real(WP), dimension(1:,1:,this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: A !< Needs to be (:,:,imin_-no:imax_+no,jmin_-no:jmax_+no,kmin_-no:kmax_+no)
	  type(MPI_Status) :: status
	  integer :: isrc,idst,ierr,isize,i,j,k,dim1,dim2
	  real(WP), dimension(:,:,:,:,:), allocatable :: buf1,buf2
	  
	  ! Get first two dimensions
	  dim1=size(A,DIM=1)
	  dim2=size(A,DIM=2)
	  
	  ! Work in x - is it 2D or 3D?
	  if (this%nx.eq.1) then
	     ! Direct copy if 2D
	     do i=this%imax_+1,this%imaxo_
		    A(:,:,i,:,:)=A(:,:,this%imin_,:,:)
	     end do
	     do i=this%imino_,this%imin_-1
		    A(:,:,i,:,:)=A(:,:,this%imin_,:,:)
	     end do
	  else
	     isize=dim1*dim2*(this%no)*(this%nyo_)*(this%nzo_)
	     allocate(buf1(dim1,dim2,this%no,this%nyo_,this%nzo_))
	     allocate(buf2(dim1,dim2,this%no,this%nyo_,this%nzo_))
	     ! Send left buffer to left neighbour
	     call MPI_CART_SHIFT(this%comm,0,-1,isrc,idst,ierr)
	     buf1=A(:,:,this%imin_:this%imin_+this%no-1,:,:)
	     call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
	     if (isrc.ne.MPI_PROC_NULL) A(:,:,this%imax_+1:this%imaxo_,:,:)=buf2
	     ! Send right buffer to right neighbour
	     call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
	     buf1=A(:,:,this%imax_-this%no+1:this%imax_,:,:)
	     call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
	     if (isrc.ne.MPI_PROC_NULL) A(:,:,this%imino_:this%imin_-1,:,:)=buf2
	     ! Deallocate
	     deallocate(buf1,buf2)
	  end if
	  
	  ! Work in y - is it 2D or 3D?
	  if (this%ny.eq.1) then
	     ! Direct copy if 2D
	     do j=this%jmax_+1,this%jmaxo_
		    A(:,:,:,j,:)=A(:,:,:,this%jmin_,:)
	     end do
	     do j=this%jmino_,this%jmin_-1
		    A(:,:,:,j,:)=A(:,:,:,this%jmin_,:)
	     end do
	  else
	     isize=dim1*dim2*(this%nxo_)*(this%no)*(this%nzo_)
	     allocate(buf1(dim1,dim2,this%nxo_,this%no,this%nzo_))
	     allocate(buf2(dim1,dim2,this%nxo_,this%no,this%nzo_))
	     ! Send left buffer to left neighbour
	     call MPI_CART_SHIFT(this%comm,1,-1,isrc,idst,ierr)
	     buf1=A(:,:,:,this%jmin_:this%jmin_+this%no-1,:)
	     call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
	     if (isrc.ne.MPI_PROC_NULL) A(:,:,:,this%jmax_+1:this%jmaxo_,:)=buf2
	     ! Send right buffer to right neighbour
	     call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
	     buf1=A(:,:,:,this%jmax_-this%no+1:this%jmax_,:)
	     call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
	     if (isrc.ne.MPI_PROC_NULL) A(:,:,:,this%jmino_:this%jmin_-1,:)=buf2
	     ! Deallocate
	     deallocate(buf1,buf2)
	  end if
	  
	  ! Work in z - is it 2D or 3D?
	  if (this%nz.eq.1) then
	     ! Direct copy if 2D
	     do k=this%kmax_+1,this%kmaxo_
		    A(:,:,:,:,k)=A(:,:,:,:,this%kmin_)
	     end do
	     do k=this%kmino_,this%kmin_-1
		    A(:,:,:,:,k)=A(:,:,:,:,this%kmin_)
	     end do
	  else
	     isize=dim1*dim2*(this%nxo_)*(this%nyo_)*(this%no)
	     allocate(buf1(dim1,dim2,this%nxo_,this%nyo_,this%no))
	     allocate(buf2(dim1,dim2,this%nxo_,this%nyo_,this%no))
	     ! Send left buffer to left neighbour
	     call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
	     buf1=A(:,:,:,:,this%kmin_:this%kmin_+this%no-1)
	     call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
	     if (isrc.ne.MPI_PROC_NULL) A(:,:,:,:,this%kmax_+1:this%kmaxo_)=buf2
	     ! Send right buffer to right neighbour
	     call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
	     buf1=A(:,:,:,:,this%kmax_-this%no+1:this%kmax_)
	     call MPI_SENDRECV(buf1,isize,MPI_REAL_WP,idst,0,buf2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
	     if (isrc.ne.MPI_PROC_NULL) A(:,:,:,:,this%kmino_:this%kmin_-1)=buf2
	     ! Deallocate
	     deallocate(buf1,buf2)
	  end if
	  
   end subroutine pgrid_rsync_tensor
   
   
   !> Synchronization of overlap cells for integer
   !> This version is capable of handling any overlap size
   subroutine pgrid_isync_no(this,A,no)
      implicit none
      class(pgrid), intent(in) :: this
      integer, intent(in) :: no
      integer, dimension(this%imin_-no:,this%jmin_-no:,this%kmin_-no:), intent(inout) :: A !< Needs to be (imin_-no:imax_+no,jmin_-no:jmax_+no,kmin_-no:kmax_+no)
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr,isize,i,j,k
      integer, dimension(:,:,:), allocatable :: buf1,buf2
      
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
         call MPI_SENDRECV(buf1,isize,MPI_INTEGER,idst,0,buf2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imax_+1:this%imax_+no,:,:)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
         buf1=A(this%imax_-no+1:this%imax_,:,:)
         call MPI_SENDRECV(buf1,isize,MPI_INTEGER,idst,0,buf2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
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
         call MPI_SENDRECV(buf1,isize,MPI_INTEGER,idst,0,buf2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmax_+1:this%jmax_+no,:)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
         buf1=A(:,this%jmax_-no+1:this%jmax_,:)
         call MPI_SENDRECV(buf1,isize,MPI_INTEGER,idst,0,buf2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
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
         do k=this%kmin_-no,this%kmin_-1
            A(:,:,k)=A(:,:,this%kmin_)
         end do
      else
         isize=(this%nx_+2*no)*(this%ny_+2*no)*(no)
         allocate(buf1(this%nx_+2*no,this%ny_+2*no,no))
         allocate(buf2(this%nx_+2*no,this%ny_+2*no,no))
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
         buf1=A(:,:,this%kmin_:this%kmin_+no-1)
         call MPI_SENDRECV(buf1,isize,MPI_INTEGER,idst,0,buf2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmax_+1:this%kmax_+no)=buf2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
         buf1=A(:,:,this%kmax_-no+1:this%kmax_)
         call MPI_SENDRECV(buf1,isize,MPI_INTEGER,idst,0,buf2,isize,MPI_INTEGER,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmin_-no:this%kmin_-1)=buf2
         ! Deallocate
         deallocate(buf1,buf2)
      end if
      
   end subroutine pgrid_isync_no
   
   
   !> Synchronization by summation of overlap cells - uses full no
   !> This routine assumes that the default overlap size is used
   !> It allows the use of pre-allocated buffers for speed
   subroutine pgrid_rsyncsum(this,A)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pgrid), intent(inout) :: this
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: A !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr,isize,i,j,k
      
      ! Work in x - is it 2D or 3D?
      if (this%nx.eq.1) then
         ! Sum along x
         do i=this%imax_+1,this%imaxo_
            A(this%imin,:,:)=A(this%imin,:,:)+A(i,:,:)
         end do
         do i=this%imino_,this%imin_-1
            A(this%imin,:,:)=A(this%imin,:,:)+A(i,:,:)
         end do
      else
         isize=(this%no)*(this%nyo_)*(this%nzo_)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,0,-1,isrc,idst,ierr)
         this%syncbuf_x1=A(this%imino_:this%imin_-1,:,:)
         call MPI_SENDRECV(this%syncbuf_x1,isize,MPI_REAL_WP,idst,0,this%syncbuf_x2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imax_-this%no+1:this%imax_,:,:)=A(this%imax_-this%no+1:this%imax_,:,:)+this%syncbuf_x2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,0,+1,isrc,idst,ierr)
         this%syncbuf_x1=A(this%imax_+1:this%imaxo_,:,:)
         call MPI_SENDRECV(this%syncbuf_x1,isize,MPI_REAL_WP,idst,0,this%syncbuf_x2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(this%imin_:this%imin_+this%no-1,:,:)=A(this%imin_:this%imin_+this%no-1,:,:)+this%syncbuf_x2
      end if
      
      ! Work in y - is it 2D or 3D?
      if (this%ny.eq.1) then
         ! Sum along x
         do j=this%jmax_+1,this%jmaxo_
            A(:,this%jmin,:)=A(:,this%jmin,:)+A(:,j,:)
         end do
         do j=this%jmino_,this%jmin_-1
            A(:,this%jmin,:)=A(:,this%jmin,:)+A(:,j,:)
         end do
      else
         isize=(this%nxo_)*(this%no)*(this%nzo_)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,1,-1,isrc,idst,ierr)
         this%syncbuf_y1=A(:,this%jmino_:this%jmin_-1,:)
         call MPI_SENDRECV(this%syncbuf_y1,isize,MPI_REAL_WP,idst,0,this%syncbuf_y2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmax_-this%no+1:this%jmax_,:)=A(:,this%jmax_-this%no+1:this%jmax_,:)+this%syncbuf_y2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,1,+1,isrc,idst,ierr)
         this%syncbuf_y1=A(:,this%jmax_+1:this%jmaxo_,:)
         call MPI_SENDRECV(this%syncbuf_y1,isize,MPI_REAL_WP,idst,0,this%syncbuf_y2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,this%jmin_:this%jmin_+this%no-1,:)=A(:,this%jmin_:this%jmin_+this%no-1,:)+this%syncbuf_y2
      end if
      
      ! Work in z - is it 2D or 3D?
      if (this%nz.eq.1) then
         ! Sum along z
         do k=this%kmax_+1,this%kmaxo_
            A(:,:,this%kmin)=A(:,:,this%kmin)+A(:,:,k)
         end do
         do k=this%kmino_,this%kmin_-1
            A(:,:,this%kmin)=A(:,:,this%kmin)+A(:,:,k)
         end do
      else
         isize=(this%nxo_)*(this%nyo_)*(this%no)
         ! Send left buffer to left neighbour
         call MPI_CART_SHIFT(this%comm,2,-1,isrc,idst,ierr)
         this%syncbuf_z1=A(:,:,this%kmino_:this%kmin_-1)
         call MPI_SENDRECV(this%syncbuf_z1,isize,MPI_REAL_WP,idst,0,this%syncbuf_z2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmax_-this%no+1:this%kmax_)=A(:,:,this%kmax_-this%no+1:this%kmax_)+this%syncbuf_z2
         ! Send right buffer to right neighbour
         call MPI_CART_SHIFT(this%comm,2,+1,isrc,idst,ierr)
         this%syncbuf_z1=A(:,:,this%kmax_+1:this%kmaxo_)
         call MPI_SENDRECV(this%syncbuf_z1,isize,MPI_REAL_WP,idst,0,this%syncbuf_z2,isize,MPI_REAL_WP,isrc,0,this%comm,status,ierr)
         if (isrc.ne.MPI_PROC_NULL) A(:,:,this%kmin_:this%kmin_+this%no-1)=A(:,:,this%kmin_:this%kmin_+this%no-1)+this%syncbuf_z2
      end if
      
      ! Follow by a sync step
      call this%pgrid_rsync(A)

   end subroutine pgrid_rsyncsum

   
   !> Returns the closest local indices "ind" to the provided position "pos" with initial guess "ind_guess"
   function get_ijk_local(this,pos,ind_guess) result(ind)
      implicit none
      class(pgrid), intent(in) :: this
      real(WP), dimension(3), intent(in) :: pos
      integer,  dimension(3), intent(in) :: ind_guess
      integer,  dimension(3) :: ind
      ! X direction
      ind(1)=ind_guess(1)
      do while (pos(1).gt.this%x(ind(1)+1).and.ind(1).lt.this%imaxo_); ind(1)=ind(1)+1; end do
      do while (pos(1).lt.this%x(ind(1)  ).and.ind(1).gt.this%imino_); ind(1)=ind(1)-1; end do
      ! Y direction
      ind(2)=ind_guess(2)
      do while (pos(2).gt.this%y(ind(2)+1).and.ind(2).lt.this%jmaxo_); ind(2)=ind(2)+1; end do
      do while (pos(2).lt.this%y(ind(2)  ).and.ind(2).gt.this%jmino_); ind(2)=ind(2)-1; end do
      ! Z direction
      ind(3)=ind_guess(3)
      do while (pos(3).gt.this%z(ind(3)+1).and.ind(3).lt.this%kmaxo_); ind(3)=ind(3)+1; end do
      do while (pos(3).lt.this%z(ind(3)  ).and.ind(3).gt.this%kmino_); ind(3)=ind(3)-1; end do
   end function get_ijk_local
   
   
   !> Returns the closest global indices "ind" to the provided position "pos" with initial guess "ind_guess"
   function get_ijk_global(this,pos,ind_guess) result(ind)
      implicit none
      class(pgrid), intent(in) :: this
      real(WP), dimension(3), intent(in) :: pos
      integer,  dimension(3), intent(in) :: ind_guess
      integer,  dimension(3) :: ind
      ! X direction
      ind(1)=ind_guess(1)
      do while (pos(1).gt.this%x(ind(1)+1).and.ind(1).lt.this%imaxo); ind(1)=ind(1)+1; end do
      do while (pos(1).lt.this%x(ind(1)  ).and.ind(1).gt.this%imino); ind(1)=ind(1)-1; end do
      ! Y direction
      ind(2)=ind_guess(2)
      do while (pos(2).gt.this%y(ind(2)+1).and.ind(2).lt.this%jmaxo); ind(2)=ind(2)+1; end do
      do while (pos(2).lt.this%y(ind(2)  ).and.ind(2).gt.this%jmino); ind(2)=ind(2)-1; end do
      ! Z direction
      ind(3)=ind_guess(3)
      do while (pos(3).gt.this%z(ind(3)+1).and.ind(3).lt.this%kmaxo); ind(3)=ind(3)+1; end do
      do while (pos(3).lt.this%z(ind(3)  ).and.ind(3).gt.this%kmino); ind(3)=ind(3)-1; end do
   end function get_ijk_global
   
   
   !> Returns the rank that contains provided indices
   function get_rank(this,ind) result(rank)
      implicit none
      class(pgrid), intent(in) :: this
      integer, dimension(3), intent(in) :: ind
      integer :: rank,ierr
      call MPI_CART_RANK(this%comm,[this%xcoord(ind(1)),this%ycoord(ind(2)),this%zcoord(ind(3))],rank,ierr)
   end function get_rank
   
   
   !> Function that returns an (i,j,k) index from a lexicographic index
   pure function get_ijk_from_lexico(this,lexico) result(ijk)
      implicit none
      class(pgrid), intent(in) :: this
      integer, intent(in) :: lexico
      integer, dimension(3) :: ijk
      ijk(3)=lexico/(this%nxo_*this%nyo_)
      ijk(2)=(lexico-this%nxo_*this%nyo_*ijk(3))/this%nxo_
      ijk(1)= lexico-this%nxo_*this%nyo_*ijk(3)-this%nxo_*ijk(2)
      ijk=ijk+[this%imino_,this%jmino_,this%kmino_]
   end function get_ijk_from_lexico
   
   
   !> Function that returns a lexicographic index from an (i,j,k) index
   pure function get_lexico_from_ijk(this,ijk) result(lexico)
      implicit none
      class(pgrid), intent(in) :: this
      integer, dimension(3), intent(in) :: ijk
      integer :: lexico
      lexico=(ijk(1)-this%imino_)+(ijk(2)-this%jmino_)*this%nxo_+(ijk(3)-this%kmino_)*this%nxo_*this%nyo_
   end function get_lexico_from_ijk
   
   
end module pgrid_class
