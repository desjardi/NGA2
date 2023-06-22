module event_class
   use precision,         only: WP
   use string,            only: str_medium
   use timetracker_class, only: timetracker
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: event
   
   !> Safety parameter to avoid barely missing an occurence due to round off
   real(WP), parameter :: dtsafe=1.0e-6_WP
   
   !> Event object
   type :: event
      character(len=str_medium)   :: name='UNNAMED_EVENT'   !< Name for event
      class(timetracker), pointer :: time                   !< Timetracker for event
      integer  :: nper                                      !< Period in elapsed number of time steps
      real(WP) :: tper                                      !< Period in elapsed time
   contains
      procedure :: occurs                                   !< Check if event is occuring
   end type event
   
   
   !> Declare event constructor
   interface event
      procedure constructor
   end interface event
   
   
contains
   
   
   !> Constructor for event object
   function constructor(time,name) result(self)
      implicit none
      type(event) :: self
      class(timetracker), target, intent(in) :: time
      character(len=*), optional :: name
      ! Point to timetracker object
      self%time=>time
      ! Set the event name
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 periods for event
      self%nper=0
      self%tper=0.0_WP
   end function constructor
   
   
   !> Occurence check
   logical function occurs(this)
      implicit none
      class(event), intent(in) :: this
      ! Assume not occuring
      occurs=.false.
      ! Go through standard occurence tests
      if (this%nper.gt.0) then
         if (mod(this%time%n,this%nper).eq.0) occurs=.true.
      end if
      if (this%tper.gt.0.0_WP) then
         if (mod(this%time%t+dtsafe*this%time%dt,this%tper).lt.this%time%dt*(1.0_WP-dtsafe)) occurs=.true.
      end if
   end function occurs
   
   
end module event_class
