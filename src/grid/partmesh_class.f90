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
      integer :: nvar                                                   !< Number of particle scalars stored
      integer :: nvec                                                   !< Number of particle vectors stored
      real(WP), dimension(:,:), allocatable :: var                      !< Particle scalar storage
      real(WP), dimension(:,:,:), allocatable :: vec                    !< Particle vector storage
      character(len=str_medium), dimension(:), allocatable :: varname   !< Name of particle scalar fields
      character(len=str_medium), dimension(:), allocatable :: vecname   !< Name of particle vector fields
   contains
      procedure :: reset                                   !< Reset particle mesh to zero size
      procedure :: set_size                                !< Set particle mesh to provided size
      procedure :: finalize                                !< Finalize particle mesh
   end type partmesh
   
   
   !> Declare particle mesh constructor
   interface partmesh
      procedure constructor
   end interface partmesh
   
   
contains
   
   
   !> Constructor for particle mesh object
   function constructor(nvar,nvec,name) result(self)
      implicit none
      type(partmesh) :: self
      integer, intent(in) :: nvar,nvec
      character(len=*), optional :: name
      ! Set the name of the particle mesh
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 particles
      self%n=0
      ! Initialize additional variables
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname='' !< Users will set the name themselves
      self%nvec=nvec
      allocate(self%vecname(self%nvec))
      self%vecname='' !< Users will set the name themselves
   end function constructor
   
   
   !> Reset mesh storage
   subroutine reset(this)
      implicit none
      class(partmesh), intent(inout) :: this
      this%n=0
      if (allocated(this%pos)) deallocate(this%pos)
      if (allocated(this%var)) deallocate(this%var)
      if (allocated(this%vec)) deallocate(this%vec)
   end subroutine reset
   
   
   ! Set mesh storage size
   subroutine set_size(this,size)
      implicit none
      class(partmesh), intent(inout) :: this
      integer, intent(in) :: size
      this%n=size
      allocate(this%pos(        3          ,this%n)); this%pos=0.0_WP
      allocate(this%vec(        3,this%nvec,this%n)); this%vec=0.0_WP
      allocate(this%var(          this%nvar,this%n)); this%var=0.0_WP
   end subroutine set_size
   
   
   !> Finalize partmesh
   subroutine finalize(this)
      implicit none
      class(partmesh), intent(inout) :: this
      call this%reset()
      if (allocated(this%varname)) deallocate(this%varname); this%nvar=0
      if (allocated(this%vecname)) deallocate(this%vecname); this%nvec=0
      this%name='PARTMESH'
   end subroutine finalize
   
   
end module partmesh_class
