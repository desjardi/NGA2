!> Module handling standard output to the screen
!> or to text files for monitoring purposes.
!> @todo Handle column type automatically via interfacing
module monitor
   use precision, only: WP
   use string, only: str_medium,str_long
   implicit none
   private
   
   ! General error handling
   public :: die,warn,log
   
   ! Monitoring to files/screen
   public :: monitor_dump_timestep,monitor_dump_iteration,monitor_dump_screen
   public :: monitor_init,monitor_final,monitor_create_file
   public :: monitor_set_header,monitor_set_value,monitor_set_array
   
   ! Timing tools
   public :: timer_add,timer_start,timer_stop
   
   !> Log file
   integer :: logfile
   
   !> Preset some length and formats for the columns
   integer,          parameter :: col_len=14
   character(len=*), parameter :: f1='(a12)'
   character(len=*), parameter :: f2='(i12)'
   character(len=*), parameter :: f3='(ES12.5)'
   
   ! File numbers
   integer :: nfiles  !< Total number of files
   integer :: iscreen !< Latest screen datafile
   
   !> Definition of type 'mfile'
   type :: mfile
      ! Data to dump
      integer :: ncol   !< Number of columns to output
      character(len=str_medium), dimension(:), allocatable :: col_head
      character(len=str_medium), dimension(:), allocatable :: col_type
      real(WP),                  dimension(:), allocatable :: val
      ! File information
      integer :: iunit
      character(len=str_medium) :: filename
      logical :: isnew  !< True: header needs to be written, False: no need to write header
      ! Type of monitor file
      character(len=str_medium) :: type   !< type='screen':    monitoring to the screen (only latest definition is used)
                                          !! type='timestep':  file dumped once per timestep
                                          !! type='iteration': file dumped once per subiteration
   end type mfile
   
   !> Array of mfiles
   type(mfile), dimension(:), allocatable :: mfiles
   
   !> Definition of a timer
   type :: timer_type
      character(len=str_medium) :: name
      real(WP) :: time
      real(WP) :: time_in
      logical  :: started
   end type timer_type
   
   !> Array of timers and their number
   type(timer_type), dimension(:), allocatable :: timers
   integer :: ntimers
   
   !> Timing for full timestep
   real(WP) :: full_time_in
   
   !> Array of values to transfer to monitor
   real(WP), dimension(:), allocatable :: mval
   
contains
   
   
   !> Write the header of the file
   subroutine monitor_first_dump(myfile)
      implicit none
      integer, intent(in) :: myfile
      integer :: icol,offset,index1
      character(len=col_len)    :: col
      character(len=str_long)   :: header
      character(len=str_medium) :: buffer
      logical :: twolines
      ! Create the header
      select case(trim(mfiles(myfile)%type))
      case ('screen')
         write(header(1+0*col_len:),f1) 'Step'
         write(header(1+1*col_len:),f1) 'Time'
         offset=2*col_len
      case ('timestep')
         write(header(1+0*col_len:),f1) 'Step'
         write(header(1+1*col_len:),f1) 'Time'
         offset=2*col_len
      case ('iteration')
         write(header(1+0*col_len:),f1) 'Step'
         write(header(1+1*col_len:),f1) 'Time'
         write(header(1+2*col_len:),f1) 'Niter'
         offset=3*col_len
      end select
      ! Extract the first line, detect if a second line is needed
      twolines=.false.
      do icol=1,mfiles(myfile)%ncol
         read(mfiles(myfile)%col_head(icol),'(a)') buffer
         index1=index(trim(buffer),' ')
         if (index1.ne.0.and.index1.lt.col_len-1) then
            twolines=.true.
            read(buffer(1:index1),f1) col
         else
            read(buffer,f1) col
         end if
         write(header(offset+1+(icol-1)*col_len:),f1) trim(col)
      end do
      ! Dump the first line of header
      write(mfiles(myfile)%iunit,'(a)') trim(header)
      ! Dump second line if needed
      if (twolines) then
         header=''
         do icol=1,mfiles(myfile)%ncol
            read (mfiles(myfile)%coL_head(icol),'(a)') buffer
            index1=index(trim(buffer),' ')
            if (index1.ne.0.and.index1.lt.col_len-1) then
               read(buffer(index1:),f1) col
            else
               col=''
            end if
            write(header(offset+1+(icol-1)*col_len:),f1) trim(col)
         end do
         write(mfiles(myfile)%iunit,'(a)') trim(header)
      end if
      ! File is not new anymore
      mfiles(myfile)%isnew=.false.
   end subroutine monitor_first_dump
   
   
   !> Dump the values to the files at each timestep
   subroutine monitor_dump_timestep(ntime,time)
      use parallel, only: amroot
      implicit none
      integer,  intent(in) :: ntime
      real(WP), intent(in) :: time
      integer :: icol,offset,i
      character(len=str_long) :: line
      logical, save :: first=.true.
      ! Handle timer monitoring here
      if (first) then
         call timer_monitor_init
         first=.false.
      end if
      call timer_monitor
      ! Only the root process does something
      if (.not.amroot) return
      ! Loop over all files
      do i=1,nfiles
         ! Test if we need to dump the values
         if (trim(mfiles(i)%type).eq.'timestep') then
            ! Check if we need to open the file and write the header
            if (mfiles(i)%isnew) call monitor_first_dump(i)
            ! Start the line to dump with time step info
            write(line(1+0*col_len:),f2) ntime
            write(line(1+1*col_len:),f3)  time
            offset=2*col_len
            ! Add all variables passed to monitor
            do icol=1,mfiles(i)%ncol
               select case(mfiles(i)%col_type(icol))
               case('i')
                  write(line(offset+1+(icol-1)*col_len:),f2) int(mfiles(i)%val(icol))
               case('r')
                  write(line(offset+1+(icol-1)*col_len:),f3) real(mfiles(i)%val(icol),WP)
               end select
            end do
            ! Dump the line
            write(mfiles(i)%iunit,'(a)') trim(line)
            ! Flush file
            call flush(mfiles(i)%iunit)
         end if
      end do
   end subroutine monitor_dump_timestep
   
   
   !> Dump the values to the files at each subiteration
   subroutine monitor_dump_iteration(ntime,time,niter)
      use parallel, only: amroot
      implicit none
      integer,  intent(in) :: ntime,niter
      real(WP), intent(in) :: time
      integer :: icol,offset,i
      character(len=str_long) :: line
      ! Only the root process does something
      if (.not.amroot) return
      ! Loop over all files
      do i=1,nfiles
         ! Test if we need to dump the values
         if (trim(mfiles(i)%type).eq.'iteration') then
            ! Check if we need to open the file and write the header
            if (mfiles(i)%isnew) call monitor_first_dump(i)
            ! Start the line to dump with time step and iteration info
            write(line(1+0*col_len:),f2) ntime
            write(line(1+1*col_len:),f3) time
            write(line(1+2*col_len:),f2) niter
            offset=3*col_len
            ! Add all variables passed to monitor
            do icol=1,mfiles(i)%ncol
               select case(mfiles(i)%col_type(icol))
               case('i')
                  write(line(offset+1+(icol-1)*col_len:),f2) int(mfiles(i)%val(icol))
               case('r')
                  write(line(offset+1+(icol-1)*col_len:),f3) real(mfiles(i)%val(icol),WP)
               end select
            end do
            ! Dump the line
            write(mfiles(i)%iunit,'(a)') trim(line)
            ! Flush file
            call flush(mfiles(i)%iunit)
         end if
      end do
   end subroutine monitor_dump_iteration
   
   
   !> Print some values to the screen
   subroutine monitor_dump_screen(ntime,time)
      use parallel, only: amroot
      implicit none
      integer,  intent(in) :: ntime
      real(WP), intent(in) :: time
      integer :: icol,offset
      character(len=str_long) :: line
      ! Only the root process does something
      if (.not.amroot) return
      ! Work only on last screen output defined
      if (mfiles(iscreen)%isnew) call monitor_first_dump(iscreen)
      ! Start the line to dump with time step and iteration info
      write(line(1+0*col_len:),f2) ntime
      write(line(1+1*col_len:),f3) time
      offset=2*col_len
      ! Add all variables passed to monitor
      do icol=1,mfiles(iscreen)%ncol
         select case(mfiles(iscreen)%col_type(icol))
         case('i')
            write(line(offset+1+(icol-1)*col_len:),f2) int(mfiles(iscreen)%val(icol))
         case('r')
            write(line(offset+1+(icol-1)*col_len:),f3) real(mfiles(iscreen)%val(icol),WP)
         end select
      end do
      ! Dump the line
      write(mfiles(iscreen)%iunit,'(a)') trim(line)
      ! Flush file
      call flush(mfiles(iscreen)%iunit)
   end subroutine monitor_dump_screen
   
   
   !> Initialize the monitor module
   subroutine monitor_init
      use parallel, only: amroot,parallel_time
      implicit none
      integer :: ierr
      ! Ensure everything is closed and empty
      nfiles=0; iscreen=0; ntimers=0
      ! Only the root process does the following
      if (amroot) then
         ! Create the monitor directory
         call execute_command_line('mkdir -p monitor')
         ! Open log file
         open(newunit=logfile,file='monitor/log',form='formatted',iostat=ierr,status='REPLACE')
         call log('******************** NGA APPLICATION STARTED ********************')
      end if
      ! Initialize timing info as well
      full_time_in=parallel_time()
   end subroutine monitor_init
   
   
   !> Finalize monitoring by closing all files
   subroutine monitor_final
      use parallel, only: amroot
      implicit none
      integer :: i
      do i=1,nfiles
         ! Deallocate all arrays
         if (allocated(mfiles(i)%col_head)) deallocate(mfiles(i)%col_head)
         if (allocated(mfiles(i)%col_type)) deallocate(mfiles(i)%col_type)
         if (allocated(mfiles(i)%val))      deallocate(mfiles(i)%val)
         ! Close files
         if (amroot.and.i.ne.iscreen) close(mfiles(i)%iunit)
      end do
      ! Close logfiles
      call log('******************** NGA APPLICATION COMPLETE *******************')
      if (amroot) close(logfile)
      ! Deallocate mfiles
      if (allocated(mfiles)) deallocate(mfiles)
      ! Zero files
      nfiles=0; iscreen=0
      ! Clean up timers
      ntimers=0
      if (allocated(mval))   deallocate(mval)
      if (allocated(timers)) deallocate(timers)
   end subroutine monitor_final
   
   
   !> Create a new file to monitor
   subroutine monitor_create_file(filename,ncol,type)
      use parallel, only: amroot
      use string, only: lowercase
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: ncol
      character(len=*), intent(in) :: type
      type(mfile), dimension(:), allocatable :: temp
      integer :: ierr,i
      ! Check that the file does not already exit
      do i=1,nfiles
         if (trim(mfiles(i)%filename).eq.trim(filename)) call die('[monitor_create_file] File has already been defined: '//trim(filename))
      end do
      ! Make room for new file
      allocate(temp(nfiles+1))
      if (nfiles.gt.0) temp(1:nfiles)=mfiles
      call move_alloc(temp,mfiles)
      ! Add the file
      nfiles=nfiles+1
      mfiles(nfiles)%filename=trim(filename)
      mfiles(nfiles)%type    =lowercase(trim(type))
      mfiles(nfiles)%ncol    =ncol
      ! Deal with screen output
      if (trim(type).eq.'screen') then
         if (iscreen.ne.0) call warn('[monitor_create_file] Screen output is being redefined!')
         iscreen=nfiles
      end if
      ! Allocate the arrays
      allocate(mfiles(nfiles)%col_head(ncol))
      allocate(mfiles(nfiles)%col_type(ncol))
      allocate(mfiles(nfiles)%val(ncol))
      ! Set file as new
      mfiles(nfiles)%isnew=.true.
      ! Open the file
      select case (trim(mfiles(nfiles)%type))
      case ('screen')
         ! Set the unit to std out
         if (amroot) mfiles(nfiles)%iunit=output_unit
      case ('timestep','iteration')
         ! Open the file
         if (amroot) open(newunit=mfiles(nfiles)%iunit,file='monitor/'//trim(mfiles(nfiles)%filename),form='formatted',iostat=ierr,status='REPLACE')
      case default
         call die('[monitor_create_file] Unknown type for file: '//trim(filename))
      end select
   end subroutine monitor_create_file
   
   
   ! Search for the file and set the allocatable
   function monitor_get_file(filename) result(ind)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: ind
      do ind=1,nfiles
         if (trim(mfiles(ind)%filename).eq.trim(filename)) return
      end do
      call die('[monitor_get_file] Could not find following file: '//trim(filename))
   end function monitor_get_file
   
   
   !> Set the header of the file
   subroutine monitor_set_header(filename,icol,header,col_type)
      use string, only: lowercase
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: icol
      character(len=*), intent(in) :: header
      character(len=1), intent(in) :: col_type
      integer :: ind
      ! Get the right file
      ind=monitor_get_file(filename)
      ! Test if enough columns
      if (icol.gt.mfiles(ind)%ncol) call die('[monitor_set_header] Column index too large for file '//trim(mfiles(ind)%filename))
      mfiles(ind)%col_head(icol)=trim(header)
      ! Test if column type is acceptable
      select case (lowercase(col_type))
      case ('i')
         mfiles(ind)%col_type(icol)='i'
      case ('r')
         mfiles(ind)%col_type(icol)='r'
      case default
         call die('[monitor_set_header] Column type unknown for file '//trim(mfiles(ind)%filename))
      end select
   end subroutine monitor_set_header
   
   
   !> Set the values to dump in the file
   subroutine monitor_set_value(filename,icol,val)
      implicit none
      character(len=*), intent(in) :: filename
      integer,  intent(in) :: icol !< Column number for the new value
      real(WP), intent(in) :: val  !< Value to dump in file
      integer :: ind
      ! Get the right file
      ind=monitor_get_file(filename)
      ! Put the data at the right place
      if (icol.gt.mfiles(ind)%ncol) then
         call die('[monitor_set_value] Column index too large for file '//trim(mfiles(ind)%filename))
      else
         mfiles(ind)%val(icol)=val
      end if
   end subroutine monitor_set_value
   
   
   !> Set all the values to dump in the file simultaneously
   subroutine monitor_set_array(filename,val)
      implicit none
      character(len=*), intent(in) :: filename
      real(WP), dimension(:), intent(in) :: val
      integer :: ind
      ! Get the right file
      ind=monitor_get_file(filename)
      ! Put the data at the right place
      if (size(val).gt.mfiles(ind)%ncol) then
         call die('[monitor_set_array] Array size too large for file '//trim(mfiles(ind)%filename))
      else
         mfiles(ind)%val=val
      end if
   end subroutine monitor_set_array
   
   
   !> Find a timer of a given name
   function timer_find(name) result(ind)
      implicit none
      character(len=*), intent(in) :: name
      integer :: ind
      ! Loop over timers and find right one
      do ind=1,ntimers
         if (trim(timers(ind)%name).eq.trim(name)) return
      end do
      ! Failed to find the timer
      call die('[timer_find] Unknown timer:'//trim(name))
   end function timer_find
   
   
   !> Create a new timer object
   subroutine timer_add(name)
      implicit none
      character(len=*), intent(in) :: name
      type(timer_type), dimension(:), allocatable :: temp
      integer :: i
      ! Check that the name does not already exist
      do i=1,ntimers
         if (trim(timers(i)%name).eq.trim(name)) call die('[timer_add] Attempted to add an already-existing timer:'//trim(name))
      end do
      ! Expend the storage
      allocate(temp(ntimers+1))
      if (ntimers.gt.0) temp(1:size(timers))=timers
      call move_alloc(temp,timers)
      ! Add the new timer
      ntimers=ntimers+1
      timers(ntimers)%name=trim(name)
      timers(ntimers)%time=0.0_WP
      timers(ntimers)%time_in=0.0_WP
      timers(ntimers)%started=.false.
   end subroutine timer_add
   
   
   !> Start the timer
   subroutine timer_start(name)
      use parallel, only: parallel_time
      implicit none
      character(len=*), intent(in) :: name
      integer :: i
      ! Find the timer
      i=timer_find(name)
      ! Check if started already
      if (timers(i)%started) call warn('[timer_start] Timer has already been started: '//trim(name))
      ! Get the time and mark timer as started
      timers(i)%time_in=parallel_time()
      timers(i)%started=.true.
   end subroutine timer_start
   
   
   !> Stop the timer
   subroutine timer_stop(name)
      use parallel, only: parallel_time
      implicit none
      character(len=*), intent(in) :: name
      integer :: i
      ! Find the timer
      i=timer_find(name)
      ! Check timer was previously started
      if (.not.timers(i)%started) call warn('[timer_stop] Timer was not started: '//trim(name))
      ! Get time and stop timer
      timers(i)%time=timers(i)%time+parallel_time()-timers(i)%time_in
      timers(i)%time_in=0.0_WP
      timers(i)%started=.false.
   end subroutine timer_stop
   
   
   !> Get all the timers and prepare for monitoring
   subroutine timer_monitor_init
      implicit none
      integer :: i
      ! Allocate the array of data to transfer to monitor
      allocate(mval(2*ntimers+3))
      ! Create the file to monitor
      call monitor_create_file('timer',2*ntimers+3,'timestep')
      call monitor_set_header('timer',1,'Total [s]','r')
      do i=1,ntimers
         call monitor_set_header('timer',i*2+0,trim(timers(i)%name)//' [s]','r')
         call monitor_set_header('timer',i*2+1,trim(timers(i)%name)//' [%]','r')
      end do
      call monitor_set_header('timer',2*ntimers+2,'Rest [s]','r')
      call monitor_set_header('timer',2*ntimers+3,'Rest [%]','r')
   end subroutine timer_monitor_init
   
   
   !> Monitor the timers
   subroutine timer_monitor
      use parallel, only: parallel_time
      implicit none
      integer :: i
      real(WP) :: full_time_out,rest_time
      ! Get time for full time step
      full_time_out=parallel_time()
      ! Create the array of values
      rest_time=0.0_WP
      mval(1)=full_time_out-full_time_in
      do i=1,ntimers
         rest_time=rest_time+timers(i)%time
         mval(2*i+0)=timers(i)%time
         mval(2*i+1)=100.0_WP*timers(i)%time/mval(1)
      end do
      mval(2*ntimers+2)=mval(1)-rest_time
      mval(2*ntimers+3)=100.0_WP*(mval(1)-rest_time)/mval(1)
      ! Transfer data to monitor
      call monitor_set_array('timer',mval)
      ! Reset all timers
      full_time_in=full_time_out
      do i=1,ntimers
         timers(i)%time=0.0_WP
         timers(i)%time_in=0.0_WP
         timers(i)%started=.false.
      end do
   end subroutine timer_monitor
   
   
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
      use parallel, only: amroot
      implicit none
      character(len=*), intent(in) :: text !< Text to be printed to the log file
      character(len=22) :: tstamp
      if (amroot) then
         call timestamp(tstamp)
         write(logfile,'(a22,1x,a2,1x,a)') tstamp,'=>',trim(text)
         call flush(logfile)
      end if
   end subroutine log
   
   
end module monitor
