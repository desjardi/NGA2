module partmesh_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: partmesh
   
   !> Particle mesh object
   type :: partmesh
      character(len=str_medium) :: name='PARTMESH'         !< Name for the particle mesh
      integer :: n                                         !< Number of particles
      real(WP), dimension(:,:), allocatable :: pos         !< Position of the particles
      real(WP), dimension(:),   allocatable :: d           !< Diameter of the particles
   contains
      procedure :: reset                                   !< Reset particle mesh to zero size
      procedure :: set_size                                !< Set particle mesh to provided size
   end type partmesh
   
   
   !> Declare particle mesh constructor
   interface partmesh
      procedure constructor
   end interface partmesh
   
   
contains
   
   
   !> Constructor for particle mesh object
   function constructor(name) result(self)
      implicit none
      type(partmesh) :: self
      character(len=*), optional :: name
      ! Set the name of the surface mesh
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 size
      self%n=0
   end function constructor
   
   
   !> Reset mesh storage
   subroutine reset(this)
      implicit none
      class(partmesh), intent(inout) :: this
      this%n=0
      if (allocated(this%pos)) deallocate(this%pos)
      if (allocated(this%d))   deallocate(this%d)
   end subroutine reset
   
   
   ! Set mesh storage size
   subroutine set_size(this,size)
      implicit none
      class(partmesh), intent(inout) :: this
      integer, intent(in) :: size
      this%n=size
      allocate(this%pos(3,this%n))
      allocate(this%d  (  this%n))
   end subroutine set_size
   
   
end module partmesh_class
