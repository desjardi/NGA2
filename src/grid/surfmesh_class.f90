module surfmesh_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: surfmesh
   
   !> Surface mesh object
   type :: surfmesh
      character(len=str_medium) :: name='UNNAMED_SURFMESH' !< Name for the surface mesh
      integer :: nVert                                     !< Number of vertices
      real(WP), dimension(:), allocatable :: xVert         !< X position of the vertices - size=nVert
      real(WP), dimension(:), allocatable :: yVert         !< Y position of the vertices - size=nVert
      real(WP), dimension(:), allocatable :: zVert         !< Z position of the vertices - size=nVert
      integer :: nPoly                                     !< Number of polygons
      integer,  dimension(:), allocatable :: polySize      !< Size of polygons - size=nPoly
      integer,  dimension(:), allocatable :: polyConn      !< Connectivity - size=sum(polySize)
      integer :: nvar                                                   !< Number of surface variables stored
      real(WP), dimension(:,:), allocatable :: var                      !< Surface variable storage
      character(len=str_medium), dimension(:), allocatable :: varname   !< Name of surface variable fields
   contains
      procedure :: reset                                   !< Reset surface mesh to zero size
      procedure :: set_size                                !< Set surface mesh to provided size
   end type surfmesh
   
   
   !> Declare surface mesh constructor
   interface surfmesh
      procedure constructor
   end interface surfmesh
   
   
contains
   
   
   !> Constructor for surface mesh object
   function constructor(nvar,name) result(self)
      implicit none
      type(surfmesh) :: self
      integer, intent(in) :: nvar
      character(len=*), optional :: name
      ! Set the name of the surface mesh
      if (present(name)) self%name=trim(adjustl(name))
      ! Default to 0 size
      self%nVert=0
      self%nPoly=0
      ! Initialize additional variables
      self%nvar=nvar
      allocate(self%varname(self%nvar))
      self%varname='' !< Users will set the name themselves
   end function constructor
   
   
   !> Reset mesh storage
   subroutine reset(this)
      implicit none
      class(surfmesh), intent(inout) :: this
      this%nPoly=0; this%nVert=0
      if (allocated(this%xVert))    deallocate(this%xvert)
      if (allocated(this%yVert))    deallocate(this%yvert)
      if (allocated(this%zVert))    deallocate(this%zvert)
      if (allocated(this%polySize)) deallocate(this%polySize)
      if (allocated(this%polyConn)) deallocate(this%polyConn)
      if (allocated(this%var))      deallocate(this%var)
   end subroutine reset
   
   
   ! Set mesh storage size - leave connectivity alone
   subroutine set_size(this,nvert,npoly)
      implicit none
      class(surfmesh), intent(inout) :: this
      integer, intent(in) :: nvert,npoly
      this%nPoly=npoly; this%nVert=nvert
      allocate(this%xVert   (this%nVert))
      allocate(this%yVert   (this%nVert))
      allocate(this%zVert   (this%nVert))
      allocate(this%polySize(this%nPoly))
      allocate(this%var     (this%nvar,this%nPoly))
   end subroutine set_size
   
   
end module surfmesh_class
