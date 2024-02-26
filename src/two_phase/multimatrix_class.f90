!> Multi matrix class: loads weights and biases of neural network into vectors and matrices
module multimatrix_class
    use precision,    only: WP
    use string,       only: str_medium
    use config_class, only: config
    implicit none

    private
 
    public :: multimatrix
 
    ! Vectors
    type :: marr1D
       real(WP), dimension(:), allocatable :: vector
       character(:), allocatable :: name
    end type marr1D
 
    ! Matrices
    type :: marr2D
       real(WP), dimension(:,:), allocatable :: matrix
       character(:), allocatable :: name
    end type marr2D
 
    type :: multimatrix
       class(config), pointer :: cfg
       character(len=str_medium) :: name='UNNAMED_MULTIMATRIX'
       character(len=str_medium) :: filename
       integer :: nvector
       type(marr1D), dimension(:), allocatable :: vectors
       integer :: nmatrix
       type(marr2D), dimension(:), allocatable :: matrices

    contains
       procedure :: print=>multimatrix_print
 
    end type multimatrix
 
    interface multimatrix
       procedure constructor
    end interface multimatrix
 
 contains
    function constructor(cfg,fdata,name) result(self)
       use messager, only: die
       use parallel, only: info_mpiio,MPI_REAL_WP
       use mpi_f08
       implicit none
       type(multimatrix)                        :: self
       class(config), target, intent(in)        :: cfg
       character(len=*)  , intent(in)           :: fdata
       character(len=*)  , intent(in), optional :: name
       integer :: ierr
       type(MPI_File) :: ifile
       type(MPI_Status) :: status
       integer :: imatrix,ivector,n_row,n_col,n_char,i

       self%cfg=>cfg
 
       if (present(name)) self%name=trim(adjustl(name))
 
       self%filename=trim(adjustl(fdata))
 
       call MPI_FILE_OPEN(self%cfg%comm,trim(self%filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
       if (ierr.ne.0) call die('[multimatrix constructor] Problem encountered while reading file: '//trim(self%filename))
 
       ! Read number of bias vectors
       call MPI_FILE_READ_ALL(ifile,self%nvector,1,MPI_INTEGER,status,ierr)
       if (self%nvector.le.0) call die('[multimatrix constructor] multimatrix object requires at least 1 vector')
 
       allocate(self%vectors(1:self%nvector))
 
       ! Read in the bias vectors
       do ivector=1,self%nvector
          call MPI_FILE_READ_ALL(ifile,n_char,1,MPI_INTEGER,status,ierr)
          allocate(character(n_char) :: self%vectors(ivector)%name)
          call MPI_FILE_READ_ALL(ifile,self%vectors(ivector)%name,n_char,MPI_CHARACTER,status,ierr)
          call MPI_FILE_READ_ALL(ifile,n_row,1,MPI_INTEGER,status,ierr)
          allocate(self%vectors(ivector)%vector(n_row))
          call MPI_FILE_READ_ALL(ifile,self%vectors(ivector)%vector,n_row,MPI_REAL_WP,status,ierr)
       end do
 
       ! Read in number of weight matrices
       call MPI_FILE_READ_ALL(ifile,self%nmatrix,1,MPI_INTEGER,status,ierr)
       if (self%nmatrix.le.0) call die('[multimatrix constructor] multimatrix object requires at least 1 matrix')
 
       allocate(self%matrices(1:self%nmatrix))
 
       ! Read in the weight matrices
       do imatrix=1,self%nmatrix
          call MPI_FILE_READ_ALL(ifile,n_char,1,MPI_INTEGER,status,ierr)
          allocate(character(n_char) :: self%matrices(imatrix)%name)
          call MPI_FILE_READ_ALL(ifile,self%matrices(imatrix)%name,n_char,MPI_CHARACTER,status,ierr)
          call MPI_FILE_READ_ALL(ifile,n_row,1,MPI_INTEGER,status,ierr)
          call MPI_FILE_READ_ALL(ifile,n_col,1,MPI_INTEGER,status,ierr)
          allocate(self%matrices(imatrix)%matrix(n_row,n_col))
          call MPI_FILE_READ_ALL(ifile,self%matrices(imatrix)%matrix,n_row*n_col,MPI_REAL_WP,status,ierr)
       end do
 
       call MPI_FILE_CLOSE(ifile,ierr)
    end function constructor
 
    !> Print Neural Network information
    subroutine multimatrix_print(this)
       use, intrinsic :: iso_fortran_env, only: output_unit
       implicit none
       class(multimatrix), intent(in) :: this
       
       if (this%cfg%amRoot) then
          write(output_unit,'("Multi matrix object [",a,"] containing [",i5,"] vectors and [",i5,"] matrices was read from file [",a,"].")') trim(this%name),this%nvector,this%nmatrix,trim(this%filename)
       end if
    end subroutine multimatrix_print
 
 end module multimatrix_class