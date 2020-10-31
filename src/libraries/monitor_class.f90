!> Module handling standard output to the screen
!> or to text files for monitoring purposes.
module monitor_class
   use precision, only: WP
   use string,    only: str_medium,str_long
   implicit none
   private
   
   ! Expose general error handling subroutines
   public :: die,warn,log
   
   !> Handle to log file
   integer, private :: logfile
   
   ! Expose constructor, initialization, and finalization for monitor
   public :: monitor,monitor_init,monitor_final
   
   !> Preset some length and formats for the columns
   integer,          parameter :: col_len=14
   character(len=*), parameter :: f1='( a12  )'
   character(len=*), parameter :: f2='( i12  )'
   character(len=*), parameter :: f3='(es12.5)'
   
   !> Type for data list
   type :: mdata
      type(mdata), pointer :: next
      character(len=str_medium) :: name
      integer , pointer :: iptr
      real(WP), pointer :: rptr
   end type mdata
   
   !> Definition of monitor object
   type :: mfile
      ! Data to dump
      integer :: ncol                                                    !< Number of columns to output
      type(mdata), pointer :: data_first                                 !< Data to dump
      ! File information
      integer :: iunit                                                   !< File handle
      character(len=str_medium) :: name                                  !< File name
      ! First dump or not
      logical :: isfirst                                                 !< Is it our first time dumping this file?
   contains
      procedure :: add_column=>add_column_real,add_column_integer        !< Add a column to the monitor file
      procedure :: write_data                                            !< Writes the content of the monitor object to a file
      procedure :: write_head                                            !< Writes the header of the monitor object to a file
   end type mfile
   
   
   ! Constructor for monitor object
   interface monitor
      procedure constructor
   end interface monitor
   
   
contains
   
   
   !> Default constructor for monitor object
   function constructor(name) result(self)
      use parallel, only: amRoot
      implicit none
      type(monitor) :: self
      character(len=*), intent(in) :: name
      ! Set the name of the monitor file and root opens it
      self%name=trim(adjustl(name))
      if (amRoot) open(newunit=self%iunit,file='monitor/'//trim(self%name),form='formatted',iostat=ierr,status='replace')
      ! Set number of columns to zero for now
      self%ncol=0
      self%data_first=>NULL()
      ! We haven't yet dumped the file
      self%isfirst=.true.
   end function constructor
   
   
   !> Add a column to the monitor file - real version
   subroutine add_column_real(this,value,name)
      implicit none
      class(monitor), intent(inout) :: this
      real(WP), target, intent(in) :: value
      character(len=*), intent(in) :: name
      type(mdata), pointer :: new_data,last_data
      ! Create a new real data entry
      allocate(new_data)
      new_data%next=>NULL()
      new_data%name=trim(adjustl(name))
      new_data%iptr=>NULL()
      new_data%rptr=>value
      ! Add it to the end of the list
      if (.not.associated(this%first_data)) then
         this%first_data=>new_data
      else
         ! Move to the end of the list
         last_data=>this%first_data
         do while (associated(last_data%next))
            last_data=>last_data%next
         end do
         last_data%next=>new_data
      end if
      ! Increment list size
      this%size=this%size+1
   end subroutine add_column_real
   
   
   !> Add a column to the monitor file - integer version
   subroutine add_column_integer(this,value,name)
      implicit none
      class(monitor), intent(inout) :: this
      integer, target, intent(in) :: value
      character(len=*), intent(in) :: name
      type(mdata), pointer :: new_data,last_data
      ! Create a new real data entry
      allocate(new_data)
      new_data%next=>NULL()
      new_data%name=trim(adjustl(name))
      new_data%iptr=>value
      new_data%rptr=>NULL()
      ! Add it to the end of the list
      if (.not.associated(this%first_data)) then
         this%first_data=>new_data
      else
         ! Move to the end of the list
         last_data=>this%first_data
         do while (associated(last_data%next))
            last_data=>last_data%next
         end do
         last_data%next=>new_data
      end if
      ! Increment list size
      this%size=this%size+1
   end subroutine add_column_integer
   
   
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
   
   
   !> Initialize the monitor module
   subroutine monitor_init
      use parallel, only: amRoot
      implicit none
      ! Only the root process does the following
      if (amRoot) then
         ! Create the monitor directory
         call execute_command_line('mkdir -p monitor')
         ! Open log file
         open(newunit=logfile,file='monitor/log',form='formatted',iostat=ierr,status='replace')
         call log('******************** NGA APPLICATION STARTED ********************')
      end if
   end subroutine monitor_init
   
   
   !> Finalize monitoring by closing all files
   subroutine monitor_final
      use parallel, only: amRoot
      implicit none
      call log('******************** NGA APPLICATION COMPLETE *******************')
      if (amRoot) close(logfile)
   end subroutine monitor_final
   
   
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
   
   
end module monitor
