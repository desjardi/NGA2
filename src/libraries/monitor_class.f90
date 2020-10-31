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
   character(len=*), parameter :: f1='( a12  )'
   character(len=*), parameter :: f2='( i12  )'
   character(len=*), parameter :: f3='(es12.5)'
   
   
   !> Type for column list
   type :: column
      type(column), pointer :: next
      character(len=str_medium) :: name
      integer , pointer :: iptr
      real(WP), pointer :: rptr
   end type column
   
   
   !> Definition of monitor object
   type :: monitor
      ! Data to dump
      integer :: ncol                                                    !< Number of columns to output
      type(column), pointer :: first_col                                 !< Data to dump
      ! File information
      integer :: iunit                                                   !< File handle
      character(len=str_medium) :: name                                  !< File name
      ! First dump or not
      logical :: isfirst                                                 !< Is it our first time dumping this file?
   contains
      procedure :: add_column=>add_column_real,add_column_integer        !< Add a column to the monitor file
      procedure :: write                                                 !< Writes the content of the monitor object to a file
      !procedure :: write_header                                          !< Writes the header of the monitor object to a file
   end type monitor
   
   
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
      integer :: ierr
      ! Set the name of the monitor file and root opens it
      self%name=trim(adjustl(name))
      if (amRoot) open(newunit=self%iunit,file='monitor/'//trim(self%name),form='formatted',iostat=ierr,status='replace')
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
      
      ! Check if we need to write out the header
      if (this%isfirst) then
         !call this%write_header()
         this%isfirst=.false.
      end if
      
      ! Write out our data
      
      
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
