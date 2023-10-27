!> Parallel data file concept is defined here: given a partitioned grid and
!> an I/O partition, it provides efficient parallel I/O access to a data file
module pardata_class
   use precision,     only: WP
   use string,        only: str_short,str_medium
   use pgrid_class,   only: pgrid
   use coupler_class, only: coupler
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: pardata   
   

   !> pardata object definition as an array of values and variables
   !> The choice is made here to *copy* the data to simplify the process
   type :: pardata
      ! pardata has two partitioned grid
      type(pgrid), pointer     :: pg                                  !< Original partitioned grid
      type(pgrid), allocatable :: pg_io                               !< Partitioned grid for I/O
      ! I/O partition info
      logical :: ingrp_io
      integer, dimension(3) :: partition_io
      ! Couplers for transfering data for reading and writing
      type(coupler), allocatable :: wcpl,rcpl
      ! Filename for read/write operations
      character(len=str_medium) :: filename
      ! A pardata stores scalar values
      integer :: nval                                                 !< Number of scalar values
      character(len=str_short), dimension(:), allocatable :: valname  !< Name of scalar values
      real(WP), dimension(:), allocatable :: val                      !< Values
      ! A pardata stores real(WP) variables of pgrid-compatible size
      integer :: nvar                                                 !< Number of field variables
      character(len=str_short), dimension(:), allocatable :: varname  !< Name of field variables
      real(WP), dimension(:,:,:,:), allocatable :: var                !< Variables (without overlap!)
      real(WP), dimension(:,:,:,:), allocatable :: var_io             !< Variables for I/O (only allocated when reading/writing)
   contains
      procedure :: initialize                     !< Initialization for pardata
      procedure, private :: prep_iomap            !< IO group/pgrid/map
      procedure, private :: findval               !< Function that returns val index if name is found, zero otherwise
      procedure, private :: findvar               !< Function that returns var index if name is found, zero otherwise
      procedure :: pushval=>pardata_pushval       !< Push data to pardata
      procedure :: pullval=>pardata_pullval       !< Pull data from pardata
      procedure :: pushvar=>pardata_pushvar       !< Push data to pardata
      procedure :: pullvar=>pardata_pullvar       !< Pull data from pardata
      procedure :: write=>pardata_write           !< Parallel write a pardata object to disk
      procedure :: print=>pardata_print           !< Print out debugging info to screen
      procedure :: log  =>pardata_log             !< Print out debugging info to log
      final     :: destructor                     !< Destructor for pardata
   end type pardata
   
   
contains
   
   
   !> Destructor for pardata object
   subroutine destructor(this)
      implicit none
      type(pardata) :: this
      if (allocated(this%valname)) deallocate(this%valname)
      if (allocated(this%val    )) deallocate(this%val    )
      if (allocated(this%varname)) deallocate(this%varname)
      if (allocated(this%var    )) deallocate(this%var    )
   end subroutine destructor
   
   
   !> Prepare the io group/pgrid/map
   subroutine prep_iomap(this,iopartition)
      use messager, only: die
      use mpi_f08,  only: MPI_Group_range_incl
      implicit none
      class(pardata) :: this
      integer, dimension(3), intent(in) :: iopartition
      type(MPI_Group) :: iogrp
      integer, dimension(3,1) :: iorange
      integer :: ierr
      
      ! Store IO partition
      this%partition_io=iopartition
      
      ! Basic checks on the partition
      if (product(this%partition_io).le.0) call die('[pardata constructor] At least one process must be assigned to IO')
      if (product(this%partition_io).gt.this%pg%nproc) call die('[pardata constructor] Number of IO processes must be less or equal to pgrid processes')
      
      ! Create IO MPI group using the first ranks (could explore other strategies here)
      iorange(:,1)=[0,product(this%partition_io)-1,1]
      call MPI_Group_range_incl(this%pg%group,1,iorange,iogrp,ierr)
      
      ! Create IO logical
      this%ingrp_io=.false.
      if (this%pg%rank.le.product(this%partition_io)-1) this%ingrp_io=.true.
      
      ! Create IO pgrid
      if (this%ingrp_io) this%pg_io=pgrid(this%pg%sgrid,iogrp,this%partition_io)
      
      ! Prepare coupler for writing
      this%wcpl=coupler(src_grp=this%pg%group,dst_grp=this%pg_io%group,name='write')
      call this%wcpl%set_src(this%pg); if (this%ingrp_io) call this%wcpl%set_dst(this%pg_io)
      call this%wcpl%initialize()
      
      ! Prepare coupler for reading
      this%rcpl=coupler(src_grp=this%pg_io%group,dst_grp=this%pg%group,name='read')
      if (this%ingrp_io) call this%rcpl%set_src(this%pg_io); call this%rcpl%set_dst(this%pg)
      call this%rcpl%initialize()
      
   end subroutine prep_iomap
   
   
   !> Initialization for a pardata object
   subroutine initialize(this,filename,nval,nvar,iopartition)
      implicit none
      class(pardata) :: this
      class(pgrid), target, intent(in) :: pg
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nval,nvar
      integer, dimension(3), intent(in) :: iopartition
      
      ! Link the partitioned grid and store the filename
      this%pg=>pg
      this%filename=trim(adjustl(filename))
      
      ! Prepare the IO-specialized partition and map
      call this%prepare_iomap(iopartition)
      
      ! Allocate storage and store names for vals
      this%nval=nval
      allocate(this%valname(this%nval))
      this%valname=''
      allocate(this%val(this%nval))
      this%val=0.0_WP
      
      ! Allocate storage and store names for vars
      this%nvar=nvar
      allocate(this%varname(this%nvar))
      this%varname=''
      allocate(this%var(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,this%pg%kmin_:this%pg%kmax_,this%nvar))
      this%var=0.0_WP
      
   end subroutine initialize
   
   
   !> Constructor for a pardata object based on an existing file
   function pardata_from_file(pg,fdata) result(self)
      use param,    only: verbose
      use messager, only: die
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      type(pardata) :: self
      character(len=*), intent(in) :: fdata
      class(pgrid), target, intent(in) :: pg
      integer :: ierr,n
      type(MPI_File) :: ifile
      type(MPI_Status) :: status
      integer, dimension(5) :: dims
      integer(kind=MPI_OFFSET_KIND) :: disp
      
      ! Link the partitioned grid and store the filename
      self%pg=>pg
      self%filename=trim(adjustl(fdata))
      
      ! Prepare the IO-specialized partition and map
      call self%prepare_iomap(iopartition)

      ! First open the file in parallel
      call MPI_FILE_OPEN(self%pg%comm,trim(self%filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[datafile constructor] Problem encountered while reading data file: '//trim(self%filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,dims,5,MPI_INTEGER,status,ierr)
      
      ! Throw error if size mismatch
      if ((dims(1).ne.self%pg%nx).or.(dims(2).ne.self%pg%ny).or.(dims(3).ne.self%pg%nz)) then
         if (self%pg%amRoot) then
            print*,'grid size = ',self%pg%nx,self%pg%ny,self%pg%nz
            print*,'data size = ',dims(1),dims(2),dims(3)
         end if
         call die('[datafile constructor] Size of datafile ['//trim(self%filename)//'] does not match pgrid ['//self%pg%name//']')
      end if
      if (dims(4).lt.0) call die('[datafile constructor] Number of values cannot be negative')
      if (dims(5).lt.0) call die('[datafile constructor] Number of variables cannot be negative')
      
      ! Allocate necessary storage
      self%nval=dims(4); allocate(self%valname(self%nval)); allocate(self%val(self%nval))
      self%nvar=dims(5); allocate(self%varname(self%nvar)); allocate(self%var(self%pg%imin_:self%pg%imax_,self%pg%jmin_:self%pg%jmax_,self%pg%kmin_:self%pg%kmax_,self%nvar))
      
      ! Read the names for all the values
      do n=1,self%nval
         call MPI_FILE_READ_ALL(ifile,self%valname(n),str_short,MPI_CHARACTER,status,ierr)
      end do
      
      ! Read all the values
      do n=1,self%nval
         call MPI_FILE_READ_ALL(ifile,self%val(n),1,MPI_REAL_WP,status,ierr)
      end do
      
      ! Read the names for all the variables
      do n=1,self%nvar
         call MPI_FILE_READ_ALL(ifile,self%varname(n),str_short,MPI_CHARACTER,status,ierr)
      end do
      
      ! Read the data for all the variables
      call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      do n=1,self%nvar
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,self%pg%view,'native',info_mpiio,ierr)
         call MPI_FILE_READ_ALL(ifile,self%var(:,:,:,n),self%pg%nx_*self%pg%ny_*self%pg%nz_,MPI_REAL_WP,status,ierr)
         disp=disp+int(self%pg%nx,MPI_OFFSET_KIND)*int(self%pg%ny,MPI_OFFSET_KIND)*int(self%pg%nz,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
      end do
      
      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! If verbose run, log and or print grid
      if (verbose.gt.1) call self%log('Read')
      if (verbose.gt.2) call self%print('Read')
      
   end function pardata_from_file
   
   
   !> Write a pardata object to a file
   subroutine pardata_write(this,fdata)
      use param,    only: verbose
      use messager, only: die
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      class(pardata), intent(in) :: this
      character(len=*), optional :: fdata
      integer :: ierr,n,iunit
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: disp
      character(len=str_medium) :: filename
      
      ! Choose the filename
      if (present(fdata)) then
         filename=trim(adjustl(fdata))
      else
         filename=trim(adjustl(this%filename))
      end if
      
      ! Root serial-writes the file header
      if (this%pg%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[datafile write] Problem encountered while serial-opening data file: '//trim(filename))
         ! Dimensions
         write(iunit) this%pg%nx,this%pg%ny,this%pg%nz,this%nval,this%nvar
         ! Name of all values
         do n=1,this%nval
            write(iunit) this%valname(n)
         end do
         ! Write all values
         do n=1,this%nval
            write(iunit) this%val(n)
         end do
         ! Name of all variables
         do n=1,this%nvar
            write(iunit) this%varname(n)
         end do
         ! Close the file
         close(iunit)
      end if
      
      ! The rest is done in parallel
      call MPI_FILE_OPEN(this%pg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[datafile write] Problem encountered while parallel-opening data file: '//trim(filename))
      
      ! Get current position
      call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      
      ! Write all variables
      do n=1,this%nvar
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,this%pg%view,'native',info_mpiio,ierr)
         call MPI_FILE_WRITE_ALL(ifile,this%var(:,:,:,n),this%pg%nx_*this%pg%ny_*this%pg%nz_,MPI_REAL_WP,status,ierr)
         disp=disp+int(this%pg%nx,MPI_OFFSET_KIND)*int(this%pg%ny,MPI_OFFSET_KIND)*int(this%pg%nz,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
      end do
      
      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! If verbose run, log and or print grid
      if (verbose.gt.1) call this%log('Wrote')
      if (verbose.gt.2) call this%print('Wrote')
      
   end subroutine pardata_write
   
   
   !> Print pardata content to the screen
   subroutine pardata_print(this,text)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(pardata), intent(in) :: this
      character(*), intent(in) :: text
      if (this%pg%amRoot) then
         write(output_unit,'(a," pardata [",a,"] for partitioned grid [",a,"]")') trim(text),trim(this%filename),trim(this%pg%name)
         write(output_unit,'(" >      nval = ",i0)') this%nval
         write(output_unit,'(" > val names = ",1000(a,1x))') this%valname
         write(output_unit,'(" >      nvar = ",i0)') this%nvar
         write(output_unit,'(" > var names = ",1000(a,1x))') this%varname
         write(output_unit,'(" > IO partition is ",i0,"x",i0,"x",i0)') this%pg_io%npx,this%pg_io%npy,this%pg_io%npz
      end if
   end subroutine pardata_print
   
   
   !> Print pardata content to the log
   subroutine pardata_log(this,text)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(pardata), intent(in) :: this
      character(*), intent(in) :: text
      character(len=str_long) :: message
      if (this%pg%amRoot) then
         write(message,'(a," pardata [",a,"] for partitioned grid [",a,"]")') trim(text),trim(this%filename),trim(this%pg%name); call log(message)
         write(message,'(" >      nval = ",i0)') this%nval; call log(message)
         write(message,'(" > val names = ",1000(a,1x))') this%valname; call log(message)
         write(message,'(" >      nvar = ",i0)') this%nvar; call log(message)
         write(message,'(" > var names = ",1000(a,1x))') this%varname; call log(message)
         write(message,'(" > IO partition is ",i0,"x",i0,"x",i0)') this%pg_io%npx,this%pg_io%npy,this%pg_io%npz; call log(message)
      end if
   end subroutine pardata_log
   
   
   !> Index finding for val
   function findval(this,name) result(ind)
      implicit none
      class(pardata), intent(in) :: this
      integer :: ind
      character(len=*), intent(in) :: name
      integer :: n
      ind=0
      do n=1,this%nval
         if (this%valname(n).eq.name) ind=n
      end do
   end function findval
   
   
   !> Push data to a val
   subroutine pardata_pushval(this,name,val)
      use messager, only: die
      implicit none
      class(pardata), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(in) :: val
      integer :: n
      n=this%findval(name)
      if (n.gt.0) then
         this%val(n)=val
      else
         call die('[pardata pushval] Val does not exist in the data file: '//name)
      end if
   end subroutine pardata_pushval
   
   
   !> Pull data from a val
   subroutine pardata_pullval(this,name,val)
      use messager, only: die
      implicit none
      class(pardata), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(out) :: val
      integer :: n
      n=this%findval(name)
      if (n.gt.0) then
         val=this%val(n)
      else
         call die('[pardata pullval] Val does not exist in the data file: '//name)
      end if
   end subroutine pardata_pullval
   
   
   !> Index finding for var
   function findvar(this,name) result(ind)
      implicit none
      class(pardata), intent(in) :: this
      integer :: ind
      character(len=*), intent(in) :: name
      integer :: n
      ind=0
      do n=1,this%nvar
         if (this%varname(n).eq.name) ind=n
      end do
   end function findvar
   
   
   !> Push data to a var
   subroutine pardata_pushvar(this,name,var)
      use messager, only: die
      implicit none
      class(pardata), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), intent(in) :: var !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n
      n=this%findvar(name)
      if (n.gt.0) then
         this%var(:,:,:,n)=var(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,this%pg%kmin_:this%pg%kmax_)
      else
         call die('[pardata pushvar] Var does not exist in the data file: '//name)
      end if
   end subroutine pardata_pushvar
   
   
   !> Pull data from a var and synchronize it
   subroutine pardata_pullvar(this,name,var)
      use messager, only: die
      implicit none
      class(pardata), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), intent(out) :: var !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n
      n=this%findvar(name)
      if (n.gt.0) then
         var(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,this%pg%kmin_:this%pg%kmax_)=this%var(:,:,:,n)
         call this%pg%sync(var)
      else
         call die('[pardata pullvar] Var does not exist in the data file: '//name)
      end if
   end subroutine pardata_pullvar
   
   
end module pardata_class
