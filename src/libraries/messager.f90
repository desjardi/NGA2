!> Module handling standard messaging to the screen
module messager
   implicit none
   private
   
   ! Expose messaging subroutines
   public :: die,warn,log
   
   !> Handle to log file
   integer, private :: logfile
   
   ! Expose initialization and finalization for messaging
   public :: messager_init,messager_final
   
contains
   
   
   !> Initialize the messager module
   subroutine messager_init
      use parallel, only: amRoot
      implicit none
      integer :: ierr
      ! Only the root process does the following
      if (amRoot) then
         ! Create the monitor directory
         call execute_command_line('mkdir -p monitor')
         ! Open log file
         open(newunit=logfile,file='monitor/log',form='formatted',iostat=ierr,status='replace')
         call log('******************** NGA APPLICATION STARTED ********************')
      end if
   end subroutine messager_init
   
   
   !> Finalize messaging by closing log file
   subroutine messager_final
      use parallel, only: amRoot
      implicit none
      call log('******************** NGA APPLICATION COMPLETE *******************')
      if (amRoot) close(logfile)
   end subroutine messager_final
   
   
   !> This routine outputs the current YMDHS date in a readeable fashion.
   subroutine timestamp(tstamp)
      implicit none
      character(len=*), intent(out) :: tstamp !< This is the timestamp returned by the routine - should be 22 long
      integer :: d,h,m,n,s,y
      integer, dimension(8) :: values
      character(len=3), parameter, dimension(12) :: month=(/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
      ! Get date and time and make output more readable
      call date_and_time(values=values)
      y=values(1); m=values(2); d=values(3)
      h=values(5); n=values(6); s=values(7)
      ! Generate output
      write(tstamp,'(i2,1x,a3,1x,i4,1x,a1,1x,i2,a1,i2.2,a1,i2.2)') d,month(m),y,'@',h,':',n,':',s
   end subroutine timestamp
   
   
   !> Subroutine that violently interrupts code execution
   subroutine die(error_text)
      use parallel, only: parallel_kill
      implicit none
      character(len=*), intent(in) :: error_text
      call log("<KILL> "//trim(error_text))
      call parallel_kill(trim(error_text))
   end subroutine die
   
   
   !> Subroutine that warns of issue
   subroutine warn(warn_text)
      implicit none
      character(len=*), intent(in) :: warn_text
      call log('<WARN> '//trim(warn_text))
   end subroutine warn
   
   
   !> Monitor action in a log file
   subroutine log(text)
      use parallel, only: amRoot
      implicit none
      character(len=*), intent(in) :: text !< Text to be printed to the log file
      character(len=22) :: tstamp
      if (amRoot) then
         call timestamp(tstamp)
         write(logfile,'(a22,1x,a2,1x,a)') tstamp,'=>',trim(text)
         call flush(logfile)
      end if
   end subroutine log
   
   
end module messager
