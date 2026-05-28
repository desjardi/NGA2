!> Module providing support for resource usage tracking
module resource_tracker
   use iso_c_binding
   use precision, only: WP
   implicit none
   private
   public :: getRSS
   public :: maxRSS,minRSS,avgRSS
   
   !> Max, min, and avg RSS
   real(WP) :: maxRSS
   real(WP) :: minRSS
   real(WP) :: avgRSS
   
   !> Structure to hold resource usage information
   type, bind(C) :: rusage_struct
      integer(c_long) :: ru_utime_tv_sec,ru_utime_tv_usec
      integer(c_long) :: ru_stime_tv_sec,ru_stime_tv_usec
      integer(c_long) :: ru_maxrss,ru_ixrss,ru_idrss,ru_isrss
      integer(c_long) :: ru_minflt,ru_majflt
      integer(c_long) :: ru_nswap
      integer(c_long) :: ru_inblock,ru_oublock
      integer(c_long) :: ru_msgsnd,ru_msgrcv
      integer(c_long) :: ru_nsignals
      integer(c_long) :: ru_nvcsw,ru_nivcsw
   end type rusage_struct
   
   !> Interface to getrusage system call
   interface
      function getrusage(usage,rusage) bind(C,name='getrusage')
         import :: c_int,c_ptr
         integer(c_int), value :: usage
         type(c_ptr), value :: rusage
         integer(c_int) :: getrusage
      end function getrusage
   end interface
   
contains
   
   !> Get current RSS
   function get_process_rss() result(RSS)
      real(WP) :: RSS
      type(rusage_struct), target :: ru
      integer :: ierr
      ierr=getrusage(0_c_int,c_loc(ru))
      RSS=real(ru%ru_maxrss,WP)/1024.0_WP/1024.0_WP !< Convert RSS to MB on MacOS and GB on Linux
   end function get_process_rss
   
   
   !> Get max, min, and avg RSS
   subroutine getRSS()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
      use parallel, only: MPI_REAL_WP,comm,nproc
      real(WP) :: RSS
      integer :: ierr
      ! Each process gets its own memory usage
      RSS=get_process_rss()
      ! Get global max, min, and avg
      call MPI_ALLREDUCE(RSS,maxRSS,1,MPI_REAL_WP,MPI_MAX,comm,ierr)
      call MPI_ALLREDUCE(RSS,minRSS,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
      call MPI_ALLREDUCE(RSS,avgRSS,1,MPI_REAL_WP,MPI_SUM,comm,ierr); avgRSS=avgRSS/real(nproc,WP)
   end subroutine getRSS
   
   
end module resource_tracker