!> Single-grid config concept is defined here.
!> The intention is to have a concept that mimics the original NGA
module config_class
   use precision,   only: WP
   use bgrid_class, only: bgrid
   use pgrid_class, only: pgrid
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: config
   
   !> Config object definition
   !> A "config" is essentially NGA's old config concept, but enhanced a bit.
   !> Currently, it contains a pgrid, periodicity info, more metrics, case name,
   !> geometry info, and boundary conditions
   type, extends(pgrid) :: config
      ! Some more metrics
      real(WP), dimension(:,:,:), allocatable :: vol           !< Local cell volume
      real(WP), dimension(:,:,:), allocatable :: meshsize      !< Local effective cell size
      real(WP) :: min_meshsize                                 !< Global minimum mesh size
      ! Wall geometry - mask=0 is fluid, mask=1 is wall
      integer,  dimension(:,:,:), allocatable :: mask          !< Masking info
      ! Boundary conditions
      integer :: xper,yper,zper                                !< Periodicity in x/y/z
      integer :: nbound                                        !< Number of boundary conditions
      type(bcond), dimension(:), pointer :: bc                 !< Storage array for boundary conditions
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
