!> Single-grid config concept is defined here.
!> The intention is to have a concept that mimics the original NGA,
!> which will be evolved into mgconfig (multi-grid), and agconfig (adaptive-grid)
module sgconfig_class
   use pgrid_class, only: pgrid
   use bgrid_class, only: bgrid
   implicit none
   private
   
   !> Single-grid configuration
   type :: sgconfig
      
      ! Periodicity info
      logical :: xper,yper,zper              !< Periodicity in x/y/z
      
      ! Grid info
      type(pgrid) :: grid                    !< Config's parallelized grid
      
      ! Wall geometry
      ! to do...
      
      ! Boundary conditions
      ! to do....
      
   contains
      procedure :: print=>sgconfig_print     !< Output configuration information to the screen
   end type sgconfig
   
   
   !> Declare single-grid config constructor
   interface sgconfig
      procedure constructor
   end interface sgconfig
   
   
contains
   
   
   !> Single-grid config constructor
   function constructor(grid,xper,yper,zper) result(self)
      use parallel, only: group
      implicit none
      type(sgconfig) :: self
      type(bgrid), intent(in) :: grid
      logical, intent(in) :: xper,yper,zper
      ! Create a parallel grid with our entire group
      self%grid=pgrid(grid,group)
      ! Set periodicity
      self%xper=xper; self%yper=yper; self%zper=zper
   end function constructor
   
   
   !> Cheap print of sgconfig info to screen
   subroutine sgconfig_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(sgconfig), intent(in) :: this
      call this%grid%print
   end subroutine sgconfig_print
   
end module sgconfig_class
