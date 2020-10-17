!> Boundary condition class
!> Combines an iterator and an action
module bcond_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: bcond
   
   !> Boundary condition object definition
   type :: bcond
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      class(config), pointer :: cfg                       !< This is the config the bcond is defined for
      type(iterator) :: itr                               !< This is the iterator
   contains
      procedure :: print=>bcond_print                     !< Output bcond to the screen
   end type bcond
   
   
   !> Declare bcond constructor
   interface bcond
      procedure constructor
   end interface bcond
   
contains
   
   
   !> Default constructor for bcond
   function constructor(cfg,name) result(self)
      implicit none
      type(bcond) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      
      
   end function constructor
   
   
   
   !> Print out info for boundary condition
   subroutine bcond_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(bcond), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Boundary condition [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine bcond_print
   
   
end module bcond_class
