!> Single-grid config concept is defined here.
!> The intention is to have a concept that mimics the original NGA,
!> which will be evolved into mgconfig (multi-grid), and agconfig (adaptive-grid)
module config_class
   use bgrid_class, only: bgrid
   use pgrid_class, only: pgrid
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: config
   
   !> Config object definition
   !> A "config" is essentially NGA's old config concept, but enhanced a bit.
   !> Currently, it contains a pgrid, periodicity info, 
   type, extends(pgrid) :: config
      
      ! Metrics
      
      ! Wall geometry
      ! to do...
      
      ! Boundary conditions
      ! to do....
      
   contains
      procedure :: print=>config_print     !< Output configuration information to the screen
   end type config
   
   
   !> Declare single-grid config constructor
   interface config
      procedure construct_from_bgrid
   end interface config
   
   
contains
   
   
   !> Single-grid config constructor
   function construct_from_bgrid(grid) result(self)
      use parallel, only: group
      implicit none
      type(config) :: self
      type(bgrid), intent(in) :: grid
      ! Create a parallel grid with our entire group
      self%pgrid=pgrid(grid,group)
   end function construct_from_bgrid
   
   
   !> Cheap print of sgconfig info to screen
   subroutine config_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(config), intent(in) :: this
      call this%print
   end subroutine config_print
   
end module config_class
