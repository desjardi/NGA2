!> Various definitions and tools for handling NGA grids
!> @todo Provide a flexible multi-grid environment
!> @todo Provide a flexible parallelization strategy
module grid_class
  use precision, only: WP
  implicit none
  private

  public :: grid_type,constructor
  
  !> Definition of a generic grid type
  type :: grid_type
     !> Communicator assigned to grid
     !integer :: comm
     !> Number of cells per direction
     integer :: nx
     integer :: ny
     integer :: nz
     !> Physical bounds
     !real(WP) :: xlo,xhi
     !real(WP) :: ylo,yhi
     !real(WP) :: zlo,zhi
   contains
     private
     !> Simple output to screen for debugging
     procedure, public :: print=>grid_print
  end type grid_type
  
contains
  
  
  !> Constructor for grid object
  function constructor(nx,ny,nz)
    implicit none
    type(grid_type) :: constructor
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    !call constructor%init_grid(nx,ny,nz)
    constructor%nx=nx
    constructor%ny=ny
    constructor%nz=nz
  end function constructor
  
  
  !> Print out a grid to screen
  subroutine grid_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    use parallel, only: amroot
    implicit none
    class(grid_type), intent(in) :: this
    ! This assumes all processors know of all grids!
    if (amroot) then
       write(output_unit,'("Mesh size is ",i4,"x",i4,"x",i4)') this%nx,this%ny,this%nz
    end if
  end subroutine grid_print
  
  
end module grid_class
