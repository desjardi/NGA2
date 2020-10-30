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
      type(iterator) :: itr                               !< This is the iterator for the bcond
      
      class(config), pointer :: cfg                       !< This is the config the bcond is defined for
      
      integer :: type                                     !< Boundary condition type
      
   contains
      procedure :: print=>bcond_print                     !< Output bcond to the screen
   end type bcond
   
   
   !> Declare bcond constructor
   interface bcond
      procedure constructor
   end interface bcond
   
contains
   
   
   !> Default constructor for bcond
   function constructor(cfg,name,type,bc_locator) result(self)
      implicit none
      type(bcond) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      interface
         logical function bc_locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function bc_locator
      end interface
      integer, intent(in) :: type
      
      ! Link to the config
      self%cfg=>cfg
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Create iterator
      self%itr=iterator(cfg,name,bc_locator)
      
      ! Store type
      self%type
      
   end function constructor
   
   
   
   !> Print out info for boundary condition
   subroutine bcond_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(bcond), intent(in) :: this
      
      ! Output
      if (this%itr%amRoot) then
         write(output_unit,'("Boundary condition [",a,"] with iterator [",a,"]")') trim(this%name),trim(this%itr%name)
      end if
      
   end subroutine bcond_print
   
   
end module bcond_class
