!> Module handling standard output to the screen
!> or to text files for monitoring purposes.
module monitor_class
   use precision, only: WP
   use string,    only: str_medium,str_long
   implicit none
   private
   
   
   ! Expose constructor for monitor object
   public :: monitor
   
   
   !> Preset some length and formats for the columns
   integer,          parameter :: col_len=14
   character(len=*), parameter :: aformat='(a12)'
   character(len=*), parameter :: iformat='(i12)'
   character(len=*), parameter :: rformat='(es12.5)'
   
   
   !> Type for column list
   type :: column
      type(column), pointer :: next
      character(len=str_medium) :: name
      integer , pointer :: iptr
      real(WP), pointer :: rptr
   end type column
   
   
   !> Definition of monitor object
   type :: monitor
      ! Assign a process responsible for writing
      logical :: amRoot                                                  !< Root is in charge of writing
      ! Data to dump
      integer :: ncol                                                    !< Number of columns to output
      type(column), pointer :: first_col                                 !< Data to dump
      ! File information
      integer :: iunit                                                   !< File handle
      character(len=str_medium) :: name                                  !< File name
      ! First dump or not
      logical :: isfirst                                                 !< Is it our first time dumping this file?
   contains
      generic :: add_column=>add_column_real,add_column_integer          !< Add a column to the monitor file
      procedure, private :: add_column_real,add_column_integer
      procedure :: write                                                 !< Writes the content of the monitor object to a file
      !procedure :: write_header                                          !< Writes the header of the monitor object to a file
   end type monitor
   
   
   ! Constructor for monitor object
   interface monitor
      procedure constructor
   end interface monitor
   
   
contains
   
   
   !> Default constructor for monitor object
   function constructor(amRoot,name) result(self)
      implicit none
      type(monitor) :: self
      logical, intent(in) :: amRoot
      character(len=*), intent(in) :: name
      integer :: ierr
      ! Set root process
      self%amRoot=amRoot
      ! Set the name of the monitor file and open it
      self%name=trim(adjustl(name))
      if (self%amRoot) open(newunit=self%iunit,file='monitor/'//trim(self%name),form='formatted',iostat=ierr,status='replace')
      ! Set number of columns to zero for now
      self%ncol=0
      self%first_col=>NULL()
      ! We haven't yet dumped the file
      self%isfirst=.true.
   end function constructor
   
   
   !> Write out monitor file
   subroutine write(this)
      implicit none
      class(monitor), intent(inout) :: this
      type(column), pointer :: my_col
      character(len=str_long)   :: line
      character(len=str_medium) :: buffer
      character(len=col_len)    :: col
      integer :: icol,index1
      logical :: twolines
      
      ! Only root works here
      if (.not.this%amRoot) return
      
      ! Check if we need to write out the header
      if (this%isfirst) then
         ! Extract the first line, detect if a second line is needed
         twolines=.false.
         ! Loop over columns
         my_col=>this%first_col
         icol=0
         do while (associated(my_col))
            ! Read that column header - first line
            read(my_col%name,'(a)') buffer
            index1=index(trim(buffer),' ')
            if (index1.ne.0.and.index1.lt.col_len-1) then
               twolines=.true.
               read(buffer(1:index1),aformat) col
            else
               read(buffer,aformat) col
            end if
            ! Write to buffer line
            write(line(1+icol*col_len:),aformat) trim(col)
            ! Increment column counter
            icol=icol+1
            ! Move to the next column
            my_col=>my_col%next
         end do
         ! Dump the first line of header
         write(this%iunit,'(a)') trim(line)
         ! Dump second line if needed
         if (twolines) then
            line=''
            ! Loop over columns
            my_col=>this%first_col
            icol=0
            do while (associated(my_col))
               ! Read that column header - 2nd line
               read(my_col%name,'(a)') buffer
               index1=index(trim(buffer),' ')
               if (index1.ne.0.and.index1.lt.col_len-1) then
                  read(buffer(index1+1:),aformat) col
               else
                  col=''
               end if
               ! Write to buffer line
               write(line(1+icol*col_len:),aformat) trim(col)
               ! Increment column counter
               icol=icol+1
               ! Move to the next column
               my_col=>my_col%next
            end do
            ! Dump the first line of header
            write(this%iunit,'(a)') trim(line)
         end if
         ! No need to dump header again in the future
         this%isfirst=.false.
      end if
      
      ! Write out all columns to an temporary line
      my_col=>this%first_col
      icol=0
      do while (associated(my_col))
         ! Write that column - integer or real
         if (associated(my_col%iptr)) write(line(1+icol*col_len:),iformat) my_col%iptr
         if (associated(my_col%rptr)) write(line(1+icol*col_len:),rformat) my_col%rptr
         ! Increment column counter
         icol=icol+1
         ! Move to the next column
         my_col=>my_col%next
      end do
      
      ! Dump the line tp the file
      write(this%iunit,'(a)') trim(line)
      
      ! Flush file (the flush statement is part of the 2003 standard)
      flush(this%iunit)
      
   end subroutine write
   
   
   !> Add a column to the monitor file - real version
   subroutine add_column_real(this,value,name)
      implicit none
      class(monitor), intent(inout) :: this
      real(WP), target, intent(in) :: value
      character(len=*), intent(in) :: name
      type(column), pointer :: new_col,last_col
      ! Create a new real column entry
      allocate(new_col)
      new_col%next=>NULL()
      new_col%name=trim(adjustl(name))
      new_col%iptr=>NULL()
      new_col%rptr=>value
      ! Add it to the end of the list
      if (.not.associated(this%first_col)) then
         this%first_col=>new_col
      else
         ! Move to the end of the list
         last_col=>this%first_col
         do while (associated(last_col%next))
            last_col=>last_col%next
         end do
         last_col%next=>new_col
      end if
      ! Increment list size
      this%ncol=this%ncol+1
   end subroutine add_column_real
   
   
   !> Add a column to the monitor file - integer version
   subroutine add_column_integer(this,value,name)
      implicit none
      class(monitor), intent(inout) :: this
      integer, target, intent(in) :: value
      character(len=*), intent(in) :: name
      type(column), pointer :: new_col,last_col
      ! Create a new integer column entry
      allocate(new_col)
      new_col%next=>NULL()
      new_col%name=trim(adjustl(name))
      new_col%iptr=>value
      new_col%rptr=>NULL()
      ! Add it to the end of the list
      if (.not.associated(this%first_col)) then
         this%first_col=>new_col
      else
         ! Move to the end of the list
         last_col=>this%first_col
         do while (associated(last_col%next))
            last_col=>last_col%next
         end do
         last_col%next=>new_col
      end if
      ! Increment list size
      this%ncol=this%ncol+1
   end subroutine add_column_integer
   
   
end module monitor_class
