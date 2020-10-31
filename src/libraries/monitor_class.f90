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
   
   
   !> Type for data list
   type :: mdata
      type(mdata), pointer :: next
      character(len=str_medium) :: name
      integer , pointer :: iptr
      real(WP), pointer :: rptr
   end type mdata
   
   
   !> Definition of monitor object
   type :: monitor
      ! Data to dump
      integer :: ncol                                                    !< Number of columns to output
      type(mdata), pointer :: first_data                                 !< Data to dump
      ! File information
      integer :: iunit                                                   !< File handle
      character(len=str_medium) :: name                                  !< File name
      ! First dump or not
      logical :: isfirst                                                 !< Is it our first time dumping this file?
   contains
      procedure :: add_column=>add_column_real,add_column_integer        !< Add a column to the monitor file
      !procedure :: write_data                                            !< Writes the content of the monitor object to a file
      !procedure :: write_head                                            !< Writes the header of the monitor object to a file
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
      self%first_data=>NULL()
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
      this%ncol=this%ncol+1
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
      this%ncol=this%ncol+1
   end subroutine add_column_integer
   
   
end module monitor_class
