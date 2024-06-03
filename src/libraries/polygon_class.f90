!> 2D polygon storage and signed distance measurement
module polygon_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private
   
   public :: polygon
   
   !> Polygon object
   type :: polygon
      character(len=str_medium) :: name='UNNAMED_POLYGON'   !< Name for polygon
      integer :: nvert                                      !< Number of vertices
      real(WP), dimension(:,:), allocatable :: vert         !< Ordered list of vertices
   contains
      procedure :: initialize                               !< Initialize a polygon object
      procedure :: get_distance                             !< Function that returns the signed distance to the polygon
   end type polygon
   
   
contains
   
   
   ! Initialization of polygon object
   subroutine initialize(this,nvert,name)
      use messager, only: die
      implicit none
      class(polygon), intent(inout) :: this
      integer, intent(in) :: nvert
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Set the number of vertices
      if (nvert.gt.1) then
         this%nvert=nvert
      else
         call die('[polygon initialize] At least 2 vertices are needed')
      end if
      ! Allocate storage for vertices
      allocate(this%vert(1:2,1:this%nvert)); this%vert=0.0_WP
   end subroutine initialize


   ! Function that returns the signed distance to a polygon object
   function get_distance(this,x) result(d)
      implicit none
      class(polygon), intent(in) :: this
      real(WP), dimension(2), intent(in) :: x
      real(WP) :: d,s
      integer :: i,i2,wn
      logical, dimension(3) :: cond
      real(WP), dimension(2) :: p,e,v
      ! Initialize distance and sign
      d=dot_product(x-this%vert(:,1),x-this%vert(:,1))
      s=1.0_WP
      ! Loop over each segment of polygon to calculate distance
      do i=1,this%nvert
         ! Handle cyclicity
         i2=i+1; if (i2.gt.this%nvert) i2=1
         ! Get distance
         e=this%vert(:,i2)-this%vert(:,i)
         v=x-this%vert(:,i)
         p=v-e*min(max(dot_product(v,e)/dot_product(e,e),0.0_WP),1.0_WP)
         d=min(d,dot_product(p,p))
         ! Get winding number
         cond=[(x(2).ge.this%vert(2,i)),(x(2).lt.this%vert(2,i2)),(e(1)*v(2).gt.e(2)*v(1))]
         if (all(cond).or.all(.not.(cond))) s=-s
      end do
      ! Return signed distance
      d=s*sqrt(d)
   end function get_distance
   
   
end module polygon_class
