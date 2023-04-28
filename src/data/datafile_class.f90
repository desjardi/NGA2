!> Datafile concept is defined here: given a partitioned grid,
!> it provides parallel I/O access to a data file
module datafile_class
   use precision,   only: WP
   use string,      only: str_short,str_medium
   use pgrid_class, only: pgrid
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: datafile
   
   
   !> Datafile object definition as an array of values and variables
   !> The choice is made here to *copy* the data to simplify reading.
   type :: datafile
      ! A datafile works in a partitioned grid
      type(pgrid), pointer :: pg                                      !< Partitioned grid for the data I/O
      ! A datafile has a filename
      character(len=str_medium) :: filename                           !< Name of datafile to read/write
      ! A datafile stores scalar values
      integer :: nval                                                 !< Number of scalar values
      character(len=str_short), dimension(:), allocatable :: valname  !< Name of scalar values
      real(WP), dimension(:), allocatable :: val                      !< Values
      ! A datafile stores real(WP) variables of pgrid-compatible size
      integer :: nvar                                                 !< Number of field variables
      character(len=str_short), dimension(:), allocatable :: varname  !< Name of field variables
      real(WP), dimension(:,:,:,:), allocatable :: var                !< Variables (without overlap!)
   contains
      procedure, private :: findval              !< Function that returns val index if name is found, zero otherwise
      procedure, private :: findvar              !< Function that returns var index if name is found, zero otherwise
      procedure :: pushval=>datafile_pushval     !< Push data to datafile
      procedure :: pullval=>datafile_pullval     !< Pull data from datafile
      procedure :: pushvar=>datafile_pushvar     !< Push data to datafile
      procedure :: pullvar=>datafile_pullvar     !< Pull data from datafile
      procedure :: write  =>datafile_write       !< Parallel write a datafile object to disk
      procedure :: print  =>datafile_print       !< Print out debugging info to screen
      procedure :: log    =>datafile_log         !< Print out debugging info to log
      final     :: destructor                    !< Destructor for datafile
   end type datafile
   
   
   !> Declare datafile constructor
   interface datafile
      procedure datafile_from_args
      procedure datafile_from_file
   end interface datafile
   
   
contains
   
   
   !> Destructor for datafile object
   subroutine destructor(this)
      implicit none
      type(datafile) :: this
      if (allocated(this%valname)) deallocate(this%valname)
      if (allocated(this%val    )) deallocate(this%val    )
      if (allocated(this%varname)) deallocate(this%varname)
      if (allocated(this%var    )) deallocate(this%var    )
   end subroutine destructor
   
   
   !> Constructor for an empty datafile object
   function datafile_from_args(pg,filename,nval,nvar) result(self)
      use messager, only: die
      implicit none
      type(datafile) :: self
      class(pgrid), target, intent(in) :: pg
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nval,nvar
      
      ! Link the partitioned grid and store the filename
      self%pg=>pg
      self%filename=trim(adjustl(filename))
      
      ! Allocate storage and store names for vals
      self%nval=nval
      allocate(self%valname(self%nval))
      self%valname=''
      allocate(self%val(self%nval))
      self%val=0.0_WP
      
      ! Allocate storage and store names for vars
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname=''
      allocate(self%var(self%pg%imin_:self%pg%imax_,self%pg%jmin_:self%pg%jmax_,self%pg%kmin_:self%pg%kmax_,self%nvar))
      self%var=0.0_WP
      
   end function datafile_from_args
   
   
   !> Constructor for a datafile object based on an existing file
   function datafile_from_file(pg,fdata) result(self)
      use param,    only: verbose
      use messager, only: die
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      type(datafile) :: self
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
      
   end function datafile_from_file
   
   
   !> Write a datafile object to a file
   subroutine datafile_write(this,fdata)
      use param,    only: verbose
      use messager, only: die
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      class(datafile), intent(in) :: this
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
      
   end subroutine datafile_write
   
   
   !> Print datafile content to the screen
   subroutine datafile_print(this,text)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(datafile), intent(in) :: this
      character(*), intent(in) :: text
      if (this%pg%amRoot) then
         write(output_unit,'(a," datafile [",a,"] for partitioned grid [",a,"]")') trim(text),trim(this%filename),trim(this%pg%name)
         write(output_unit,'(" >      nval = ",i0)') this%nval
         write(output_unit,'(" > val names = ",1000(a,1x))') this%valname
         write(output_unit,'(" >      nvar = ",i0)') this%nvar
         write(output_unit,'(" > var names = ",1000(a,1x))') this%varname
      end if
   end subroutine datafile_print
   
   
   !> Print datafile content to the log
   subroutine datafile_log(this,text)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(datafile), intent(in) :: this
      character(*), intent(in) :: text
      character(len=str_long) :: message
      if (this%pg%amRoot) then
         write(message,'(a," datafile [",a,"] for partitioned grid [",a,"]")') trim(text),trim(this%filename),trim(this%pg%name); call log(message)
         write(message,'(" >      nval = ",i0)') this%nval; call log(message)
         write(message,'(" > val names = ",1000(a,1x))') this%valname; call log(message)
         write(message,'(" >      nvar = ",i0)') this%nvar; call log(message)
         write(message,'(" > var names = ",1000(a,1x))') this%varname; call log(message)
      end if
   end subroutine datafile_log
   
   
   !> Index finding for val
   function findval(this,name) result(ind)
      implicit none
      class(datafile), intent(in) :: this
      integer :: ind
      character(len=*), intent(in) :: name
      integer :: n
      ind=0
      do n=1,this%nval
         if (this%valname(n).eq.name) ind=n
      end do
   end function findval
   
   
   !> Push data to a val
   subroutine datafile_pushval(this,name,val)
      use messager, only: die
      implicit none
      class(datafile), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(in) :: val
      integer :: n
      n=this%findval(name)
      if (n.gt.0) then
         this%val(n)=val
      else
         call die('[datafile pushval] Val does not exist in the datafile: '//name)
      end if
   end subroutine datafile_pushval
   
   
   !> Pull data from a val
   subroutine datafile_pullval(this,name,val)
      use messager, only: die
      implicit none
      class(datafile), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(out) :: val
      integer :: n
      n=this%findval(name)
      if (n.gt.0) then
         val=this%val(n)
      else
         call die('[datafile pullval] Val does not exist in the datafile: '//name)
      end if
   end subroutine datafile_pullval
   
   
   !> Index finding for var
   function findvar(this,name) result(ind)
      implicit none
      class(datafile), intent(in) :: this
      integer :: ind
      character(len=*), intent(in) :: name
      integer :: n
      ind=0
      do n=1,this%nvar
         if (this%varname(n).eq.name) ind=n
      end do
   end function findvar
   
   
   !> Push data to a var
   subroutine datafile_pushvar(this,name,var)
      use messager, only: die
      implicit none
      class(datafile), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), intent(in) :: var !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n
      n=this%findvar(name)
      if (n.gt.0) then
         this%var(:,:,:,n)=var(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,this%pg%kmin_:this%pg%kmax_)
      else
         call die('[datafile pushvar] Var does not exist in the datafile: '//name)
      end if
   end subroutine datafile_pushvar
   
   
   !> Pull data from a var and synchronize it
   subroutine datafile_pullvar(this,name,var)
      use messager, only: die
      implicit none
      class(datafile), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), intent(out) :: var !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n
      n=this%findvar(name)
      if (n.gt.0) then
         var(this%pg%imin_:this%pg%imax_,this%pg%jmin_:this%pg%jmax_,this%pg%kmin_:this%pg%kmax_)=this%var(:,:,:,n)
         call this%pg%sync(var)
      else
         call die('[datafile pullvar] Var does not exist in the datafile: '//name)
      end if
   end subroutine datafile_pullvar
   
   
end module datafile_class
