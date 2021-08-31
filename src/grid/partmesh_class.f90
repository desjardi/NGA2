module partmesh_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: partmesh
   
   !> Particle mesh object
   type :: partmesh
      character(len=str_medium) :: name='PARTMESH'                      !< Name for the particle mesh
      integer :: n                                                      !< Number of particles
      real(WP), dimension(:,:), allocatable :: pos                      !< Position of the particles
      integer :: nvar                                                   !< Number of particle variables stored
      real(WP), dimension(:,:), allocatable :: var                      !< Particle variable storage
      character(len=str_medium), dimension(:), allocatable :: varname   !< Name of particle variable fields
   contains
      procedure :: reset                                   !< Reset particle mesh to zero size
      procedure :: set_size                                !< Set particle mesh to provided size
      procedure :: add_vars                                !< Add particle variables
   end type partmesh
   
   
   !> Declare particle mesh constructor
   interface partmesh
      procedure constructor
   end interface partmesh
   
   
contains
   
   
   !> Constructor for particle mesh object
   function constructor(nvar,name) result(self)
      implicit none
      type(partmesh) :: self
      integer, intent(in) :: nvar
      character(len=*), optional :: name
      ! Set the name of the particle mesh
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 particles
      self%n=0
      ! Initialize additional variables
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname='' !< Users will set the name themselves
   end function constructor
   
   
   !> Reset mesh storage
   subroutine reset(this)
      implicit none
      class(partmesh), intent(inout) :: this
      this%n=0
      if (allocated(this%pos)) deallocate(this%pos)
      if (allocated(this%var)) deallocate(this%var)
   end subroutine reset
   
   
   ! Set mesh storage size
   subroutine set_size(this,size)
      implicit none
      class(partmesh), intent(inout) :: this
      integer, intent(in) :: size
      this%n=size
      allocate(this%pos(        3,this%n))
      allocate(this%var(this%nvar,this%n))
   end subroutine set_size
   
   
   ! Add particle variables
   subroutine add_vars(this,nvars_to_add)
      implicit none
      class(partmesh), intent(inout) :: this
      integer, intent(in) :: nvars_to_add
      this%nvar=this%nvar+nvars_to_add
      if (allocated(this%var)) then
         deallocate(this%var)
         allocate(this%var(this%nvar,this%n))
      end if
      if (allocated(this%varname)) then
         deallocate(this%varname)
         allocate(this%varname(this%nvar))
      end if
      this%varname='' !< Users will set the names themselves
   end subroutine add_vars

end module partmesh_class
