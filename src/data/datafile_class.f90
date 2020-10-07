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
   
   !> Datafile object definition as an extension of pgrid
   type :: datafile
      ! A datafile works in a partitioned grid
      type(pgrid), pointer :: pg                                      !< Partitioned grid for the data I/O
      ! A datafile has a filename
      character(len=str_medium) :: filename                           !< Name of datafile to read/write
      ! A datafile stores scalar values
      integer :: nval                                                 !< Number of scalar values
      character(len=str_short), dimension(:), allocatable :: valname  !< Name of scalar values
      real(WP), dimension(:), allocatable :: val                      !< Values
      ! A datafile stores real(WP) variables pgrid-compatible size
      integer :: nvar                                                 !< Number of field variables
      character(len=str_short), dimension(:), allocatable :: varname  !< Name of field variables
      real(WP), dimension(:,:,:,:), allocatable :: var                !< Variables
   contains
      procedure :: push=>datafile_push       !< Push data to datafile
      procedure :: pull=>datafile_pull       !< Pull data from datafile
      !procedure :: read =>datafile_read      !< Parallel read a datafile object to disk
      !procedure :: write=>datafile_write     !< Parallel write a datafile object to disk
      procedure :: print=>datafile_print     !< Print out debugging info to screen
      procedure :: log=>datafile_log         !< Print out debugging info to log
   end type datafile
   
   
   !> Declare datafile constructor
   interface datafile
      procedure datafile_from_file
   end interface datafile
   
   
contains
   
   
   !> Constructor for a datafile object
   function datafile_from_file(pg,fdata) result(self)
      use monitor,  only: die
      use parallel, only: mpi_info,MPI_REAL_WP
      use mpi
      implicit none
      type(datafile) :: self
      character(len=*), intent(in) :: fdata
      class(pgrid), target, intent(in) :: pg
      integer :: ierr,ifile,n,data_size,view
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(5) :: dims
      integer(kind=MPI_OFFSET_KIND) :: disp,full_size
      
      ! Link the partitioned grid and store the filename
      self%pg=>pg
      self%filename=trim(adjustl(fdata))
      
      ! First open the file in parallel
      call MPI_FILE_OPEN(self%pg%comm,trim(self%filename),MPI_MODE_RDONLY,mpi_info,ifile,ierr)
      if (ierr.ne.0) call die('[datafile constructor] Problem encountered while reading data file: '//trim(self%filename))
      
      ! Read file header first
      call MPI_FILE_READ_ALL(ifile,dims,5,MPI_INTEGER,status,ierr)
      
      ! Throw error if size mismatch
      if ((dims(1).ne.self%pg%nx).or.(dims(2).ne.self%pg%ny).or.(dims(3).ne.self%pg%nz)) then
         if (self%pg%amRoot) then
            print*, 'grid size = ',self%pg%nx,self%pg%ny,self%pg%nz
            print*, 'data size = ',dims(1),dims(2),dims(3)
         end if
         call die('[datafile constructor] Size of datafile ['//trim(self%filename)//'] does not match pgrid ['//self%pg%name//']')
      end if
      if (dims(4).lt.0) call die('[datafile constructor] Number of values cannot be negative')
      if (dims(5).lt.0) call die('[datafile constructor] Number of variables cannot be negative')
      
      ! Allocate necessary storage
      self%nval=dims(4); allocate(self%valname(self%nval)); allocate(self%val(self%nval))
      self%nvar=dims(5); allocate(self%varname(self%nvar)); allocate(self%var(self%nvar,self%pg%imin_:self%pg%imax_,self%pg%jmin_:self%pg%jmax_,self%pg%kmin_:self%pg%kmax_))
      
      ! Read the name and data for all the values
      do n=1,self%nval
         call MPI_FILE_READ_ALL(ifile,self%valname(n),str_short,MPI_CHARACTER,status,ierr)
         call MPI_FILE_READ_ALL(ifile,self%val(n)    ,1        ,MPI_REAL_WP  ,status,ierr)
      end do
      
      ! Size of local and full arrays
      data_size=self%pg%nx_*self%pg%ny_*self%pg%nz_
      full_size=int(self%pg%nx,MPI_OFFSET_KIND)*int(self%pg%ny,MPI_OFFSET_KIND)*int(self%pg%nz,MPI_OFFSET_KIND)
      
      ! Read the name and data for all the variables
      do n=1,self%nvar
         call MPI_FILE_READ_ALL(ifile,self%varname(n),str_short,MPI_CHARACTER,status,ierr)
         disp=int(5*4+self%nval*(str_short+WP),MPI_OFFSET_KIND)+int(n-1,MPI_OFFSET_KIND)*(full_size+WP)
         call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,view,'native',mpi_info,ierr)
         call MPI_FILE_READ_ALL(ifile,self%var(n,:,:,:),data_size,MPI_REAL_WP,status,ierr)
      end do
      
      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
      
   end function datafile_from_file
   
   
   !> Print datafile content to the screen
   subroutine datafile_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(datafile), intent(in) :: this
      if (this%pg%amRoot) then
         write(output_unit,'("Datafile [",a,"] for partitioned grid [",a,"]")') trim(this%filename),trim(this%pg%name)
         write(output_unit,'(" >      nval = ",i0)') this%nval
         write(output_unit,'(" > val names = ",1000(a,x))') this%valname
         write(output_unit,'(" >      nvar = ",i0)') this%nvar
         write(output_unit,'(" > var names = ",1000(a,x))') this%varname
      end if
   end subroutine datafile_print
   
   
   !> Print datafile content to the log
   subroutine datafile_log(this)
      use monitor,  only: log
      use string,   only: str_long
      implicit none
      class(datafile), intent(in) :: this
      character(len=str_long) :: message
      if (this%pg%amRoot) then
         write(message,'("Datafile [",a,"] for partitioned grid [",a,"]")') trim(this%filename),trim(this%pg%name); call log(message)
         write(message,'(" >      nval = ",i0)') this%nval; call log(message)
         write(message,'(" > val names = ",1000(a,x))') this%valname; call log(message)
         write(message,'(" >      nvar = ",i0)') this%nvar; call log(message)
         write(message,'(" > var names = ",1000(a,x))') this%varname; call log(message)
      end if
   end subroutine datafile_log
   
   
end module datafile_class
