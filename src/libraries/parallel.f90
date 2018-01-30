!> Parallel operations - wrappers for a number of recurring MPI operations
!> Fortran notation is used (rank goes from 1 to nproc)
module parallel
  use precision, only: WP
  use string
  implicit none
  private
  !> Number of processors
  integer :: nproc
  !> Rank of this processor
  integer :: rank
  !> Rank of root process
  integer :: root
  !> This is our communicator
  integer :: comm
  !> These are MPI type for our precision
  integer :: MPI_REAL_WP,MPI_REAL_SP,MPI_COMPLEX_WP
  include 'mpif.h'

  public :: parallel_max,parallel_min,parallel_sum,parallel_bc
  public :: parallel_init,parallel_final,parallel_time,parallel_kill
  public :: nproc,rank,root,comm
  
  !> Compute the global maximum of anything
  interface parallel_max
     module procedure parallel_max_real_3d
     module procedure parallel_max_real_0d
     module procedure parallel_max_int_3d
     module procedure parallel_max_int_0d
  end interface parallel_max
  !> Compute the global minimum of anything
  interface parallel_min
     module procedure parallel_min_real_3d
     module procedure parallel_min_real_0d
     module procedure parallel_min_int_3d
     module procedure parallel_min_int_0d
  end interface parallel_min
  !> Compute the global sum of anything
  interface parallel_sum
     module procedure parallel_sum_int_0d
     module procedure parallel_sum_int_1d
     module procedure parallel_sum_real_0d
     module procedure parallel_sum_real_1d
     module procedure parallel_sum_real_2d
     module procedure parallel_sum_real_3d
  end interface parallel_sum
  !> Perform broadcast
  interface parallel_bc
     module procedure parallel_bc_char
     module procedure parallel_bc_int_0d
     module procedure parallel_bc_int_1d
     module procedure parallel_bc_int_2d
     module procedure parallel_bc_int_3d
     module procedure parallel_bc_real_0d
     module procedure parallel_bc_real_1d
     module procedure parallel_bc_real_2d
     module procedure parallel_bc_real_3d
  end interface parallel_bc
  
contains
  
  
  !> MPI MAX
  subroutine parallel_max_real_3d(A,B)
    implicit none
    real(WP), dimension(:,:,:), intent(in) :: A
    real(WP), intent(out) :: B
    real(WP) :: C
    integer :: ierr
    C = maxval(A)
    call MPI_ALLREDUCE(C,B,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
  end subroutine parallel_max_real_3d
  subroutine parallel_max_int_3d(A,B)
    implicit none
    integer, dimension(:,:,:), intent(in)  :: A
    integer, intent(out) :: B
    integer :: C
    integer :: ierr
    C = maxval(A)
    call MPI_ALLREDUCE(C,B,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine parallel_max_int_3d
  subroutine parallel_max_real_0d(A,B)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
  end subroutine parallel_max_real_0d
  subroutine parallel_max_int_0d(A,B)
    implicit none
    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine parallel_max_int_0d
  
  
  !> MPI MIN
  subroutine parallel_min_real_3d(A,B)
    implicit none
    real(WP), dimension(:,:,:), intent(in) :: A
    real(WP), intent(out) :: B
    real(WP) :: C
    integer :: ierr
    C = minval(A)
    call MPI_ALLREDUCE(C,B,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
  end subroutine parallel_min_real_3d
  subroutine parallel_min_int_3d(A,B)
    implicit none
    integer, dimension(:,:,:), intent(in) :: A
    integer, intent(out) :: B
    integer :: C
    integer :: ierr
    C = minval(A)
    call MPI_ALLREDUCE(C,B,1,MPI_INTEGER,MPI_MIN,comm,ierr)
  end subroutine parallel_min_int_3d
  subroutine parallel_min_real_0d(A,B)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
  end subroutine parallel_min_real_0d
  subroutine parallel_min_int_0d(A,B)
    implicit none
    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_INTEGER,MPI_MIN,comm,ierr)
  end subroutine parallel_min_int_0d
  
  
  !> MPI SUM
  subroutine parallel_sum_real_0d(A,B)
    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_real_0d
  subroutine parallel_sum_int_0d(A,B)
    implicit none
    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_int_0d
  subroutine parallel_sum_real_1d(A,B)
    implicit none
    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_REAL_WP,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_real_1d
  subroutine parallel_sum_int_1d(A,B)
    implicit none
    integer, dimension(:), intent(in)  :: A
    integer, dimension(:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_INTEGER,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_int_1d
  subroutine parallel_sum_real_2d(A,B)
    implicit none
    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_REAL_WP,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_real_2d
  subroutine parallel_sum_int_2d(A,B)
    implicit none
    integer, dimension(:,:), intent(in)  :: A
    integer, dimension(:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_INTEGER,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_int_2d
  subroutine parallel_sum_real_3d(A,B)
    implicit none
    real(WP), dimension(:,:,:), intent(in)  :: A
    real(WP), dimension(:,:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_REAL_WP,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_real_3d
  subroutine parallel_sum_int_3d(A,B)
    implicit none
    integer, dimension(:,:,:), intent(in)  :: A
    integer, dimension(:,:,:), intent(out) :: B
    integer :: ierr
    call MPI_ALLREDUCE(A,B,size(A),MPI_INTEGER,MPI_SUM,comm,ierr)
  end subroutine parallel_sum_int_3d
  
  
  !> MPI_BCAST
  subroutine parallel_bc_char(A)
    implicit none
    integer :: ierr
    character(len=*) :: A
    call MPI_BCAST(A,len(A),MPI_CHARACTER,root-1,comm,ierr)
  end subroutine parallel_bc_char
  subroutine parallel_bc_int_0d(A)
    implicit none
    integer :: ierr
    integer :: A
    call MPI_BCAST(A,1,MPI_INTEGER,root-1,comm,ierr)
  end subroutine parallel_bc_int_0d
  subroutine parallel_bc_int_1d(A)
    implicit none
    integer :: ierr
    integer, dimension(:) :: A
    call MPI_BCAST(A,size(A),MPI_INTEGER,root-1,comm,ierr)
  end subroutine parallel_bc_int_1d
  subroutine parallel_bc_int_2d(A)
    implicit none
    integer :: ierr
    integer, dimension(:,:) :: A
    call MPI_BCAST(A,size(A),MPI_INTEGER,root-1,comm,ierr)
  end subroutine parallel_bc_int_2d
  subroutine parallel_bc_int_3d(A)
    implicit none
    integer :: ierr
    integer, dimension(:,:,:) :: A
    call MPI_BCAST(A,size(A),MPI_INTEGER,root-1,comm,ierr)
  end subroutine parallel_bc_int_3d
  subroutine parallel_bc_real_0d(A)
    implicit none
    integer :: ierr
    real(WP) :: A
    call MPI_BCAST(A,1,MPI_REAL_WP,root-1,comm,ierr)
  end subroutine parallel_bc_real_0d
  subroutine parallel_bc_real_1d(A)
    implicit none
    integer :: ierr
    real(WP), dimension(:) :: A
    call MPI_BCAST(A,size(A),MPI_REAL_WP,root-1,comm,ierr)
  end subroutine parallel_bc_real_1d
  subroutine parallel_bc_real_2d(A)
    implicit none
    integer :: ierr
    real(WP), dimension(:,:) :: A
    call MPI_BCAST(A,size(A),MPI_REAL_WP,root-1,comm,ierr)
  end subroutine parallel_bc_real_2d
  subroutine parallel_bc_real_3d(A)
    implicit none
    integer :: ierr
    real(WP), dimension(:,:,:) :: A
    call MPI_BCAST(A,size(A),MPI_REAL_WP,root-1,comm,ierr)
  end subroutine parallel_bc_real_3d
  
  
  !> Initialize the parallel environment
  subroutine parallel_init
    use precision, only: SP
    implicit none
    integer :: ierr
    integer :: size_real,size_dp,size_complex,size_complex_dp
    
    ! Initialize a first basic MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr) 
    rank=rank+1
    root=1
    
    ! Set MPI working precision - WP
    call MPI_TYPE_SIZE(MPI_REAL,size_real,ierr)
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,size_dp,ierr)
    if (WP.eq.size_real) then
       MPI_REAL_WP=MPI_REAL
    else if (WP.eq.size_dp) then
       MPI_REAL_WP=MPI_DOUBLE_PRECISION
    else
       call parallel_kill('Error in parallel_init: no WP equivalent in MPI')
    end if
    
    ! Set MPI single precision
    call MPI_TYPE_SIZE(MPI_REAL,size_real,ierr)
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,size_dp,ierr)
    if (SP.eq.size_real) then
       MPI_REAL_SP=MPI_REAL
    else if (SP.eq.size_dp) then
       MPI_REAL_SP=MPI_DOUBLE_PRECISION
    else
       call parallel_kill('Error in parallel_init: no SP equivalent in MPI')
    end if
    
    ! Set MPI working precision (for complex types)
    call MPI_TYPE_SIZE(MPI_COMPLEX,size_complex,ierr)
    call MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,size_complex_dp,ierr)
    if (2*WP.eq.size_complex) then
       MPI_COMPLEX_WP=MPI_COMPLEX
    else if (2*WP.eq.size_complex_dp) then
       MPI_COMPLEX_WP=MPI_DOUBLE_COMPLEX
    else
       call parallel_kill('Error in parallel_init: no complex WP equivalent in MPI')
    end if
    
    ! Create shorthand for MPI_COMM_WORLD
    comm=MPI_COMM_WORLD
    
  end subroutine parallel_init
  
  
  !> Finalize MPI environment
  subroutine parallel_final
    implicit none
    integer :: ierr
    call MPI_FINALIZE(ierr)
  end subroutine parallel_final
  
  
  !> Wrapper for time inquiry
  function parallel_time() result(wtime)
    implicit none
    real(WP) :: wtime
    wtime=MPI_WTIME()
  end function parallel_time
  
  
  !> MPI KILL
  subroutine parallel_kill(error_text)
    use iso_fortran_env, only: error_unit
    implicit none
    integer :: ierr
    character(len=*), intent(in), optional :: error_text
    ! Specify who sends the abort signal and what it means
    if (present(error_text)) then
       write(error_unit,'(a,i3,a,/,a)') '[',rank,'] initiated general abort signal due to the following error:',trim(error_text)
    else
       write(error_unit,'(a,i3,a)') '[',rank,'] initiated general abort signal due to an unspecified error'
    end if
    ! Call general abort
    call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
  end subroutine parallel_kill
  
  
end module parallel
