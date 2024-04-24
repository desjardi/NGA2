!> Basic parallel environment - global communicator, rank, root...
!> Standard notation is used (rank goes from 0 to nproc-1)
module parallel
   use precision, only: WP
   use mpi_f08
   implicit none
   private
   
   !> Start-up wallclock time
   real(WP), public, protected :: wtinit
   !> Number of processors
   integer, public, protected :: nproc
   !> Rank of this processor
   integer, public, protected :: rank
   !> Am I the global root process?
   logical, public, protected :: amRoot
   !> This is our global communicator
   type(MPI_Comm), public, protected :: comm
   !> This is our global group
   type(MPI_Group), public, protected :: group
   
   !> These are MPI type for our precision
   type(MPI_Datatype), public, protected :: MPI_REAL_WP,MPI_REAL_SP,MPI_COMPLEX_WP
   type(MPI_Datatype), public, protected :: MPI_2REAL_WP

   !> I/O info
   type(MPI_Info), public, protected :: info_mpiio
   
   public :: parallel_init,parallel_final,parallel_time,parallel_kill
   
contains
   
   
   !> Initialize the parallel environment
   subroutine parallel_init
      use precision, only: SP,WP
      implicit none
      integer :: ierr
      integer :: size_real,size_dp,size_complex,size_complex_dp,size_2real,size_2dp
      
      ! Initialize a first basic MPI environment
      call MPI_INIT(ierr)
      wtinit=MPI_WTIME()
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_GROUP(MPI_COMM_WORLD,group,ierr)
      
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
      
      ! Set MPI working precision for pairs - WP
      call MPI_TYPE_SIZE(MPI_2REAL,size_2real,ierr)
      call MPI_TYPE_SIZE(MPI_2DOUBLE_PRECISION,size_2dp,ierr)
      if (2*WP.eq.size_2real) then
         MPI_2REAL_WP=MPI_2REAL
      else if (2*WP.eq.size_2dp) then
         MPI_2REAL_WP=MPI_2DOUBLE_PRECISION
      else
         call parallel_kill('Error in parallel_init: no WP equivalent in MPI')
      end if
      
      ! Create shorthand for MPI_COMM_WORLD
      comm=MPI_COMM_WORLD
      
      ! Am I the global root?
      amRoot=(rank.eq.0)
      
      ! Also initialize I/O info
      call MPI_INFO_CREATE(info_mpiio,ierr)
      info_mpiio=MPI_INFO_NULL
      !call MPI_INFO_SET(info_mpiio,'romio_ds_read','automatic',ierr)
      !call MPI_INFO_SET(info_mpiio,'romio_cb_read','enable',ierr)
      
   end subroutine parallel_init
   
   
   !> Finalize MPI environment
   subroutine parallel_final
      implicit none
      integer :: ierr
      call MPI_FINALIZE(ierr)
   end subroutine parallel_final
   
   
   !> Wrapper for time inquiry
   function parallel_time() result(wtime)
      use precision, only: WP
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
         write(error_unit,'("[",i0,"] requested a general abort with message: ",a)') rank,trim(error_text)
      else
         write(error_unit,'("[",i0,"] requested a general abort with no message")') rank
      end if
      ! Call general abort
      call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
   end subroutine parallel_kill
   
   
end module parallel
