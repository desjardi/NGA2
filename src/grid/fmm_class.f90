!> Fast marching method class:
!> Provides support for creating a signed distance field from a vfs solution
module fmm_class
   use precision,   only: WP
   use string,      only: str_medium
   use pgrid_class, only: pgrid
   implicit none
   private

   ! Expose type/constructor/methods
   public :: fmm, distance_init_ftype

   !> fmm object definition 
   type :: fmm 
      ! This is our pgrid
      class(pgrid), pointer :: pg
      ! This is the name of the CCL
      character(len=str_medium) :: name='UNNAMED_FFM'
      ! Distance to the interface 
      real(WP), dimension(:,:,:), allocatable :: dist
      ! Tag that is true for cells with defined distance 
      logical, dimension(:,:,:), allocatable :: tag
   contains
      procedure :: initialize
      procedure :: build
   end type fmm

   !> Type of the interface function used to set initial distance
   interface
      subroutine distance_init_ftype(ind1,ind2,ind3,dist,tag)
         use precision,   only: WP
         integer, intent(in) :: ind1,ind2,ind3
         real(WP), intent(out) :: dist
         logical, intent(out) :: tag
      end subroutine distance_init_ftype
   end interface
   
contains

   !> Initialize the fmm class 
   subroutine initialize(this,pg,name)
      implicit none
      class(fmm), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Point to pgrid object
      this%pg=>pg
      ! Allocate and initialize distance array
      allocate(this%dist(this%pg%imino_:this%pg%imaxo_,this%pg%jmino_:this%pg%jmaxo_,this%pg%kmino_:this%pg%kmaxo_)); this%dist=0
      allocate(this%tag(this%pg%imino_:this%pg%imaxo_,this%pg%jmino_:this%pg%jmaxo_,this%pg%kmino_:this%pg%kmaxo_)); this%tag=.false.
   end subroutine initialize

   !> Build the distance field using the user-set distance_init function
   subroutine build(this,distance_init)
      implicit none
      class(fmm), intent(inout) :: this
      procedure(distance_init_ftype) :: distance_init
      integer :: i,j,k
      ! Loop over the grid and set the distance
      do k=this%pg%kmino_,this%pg%kmaxo_
         do j=this%pg%jmino_,this%pg%jmaxo_
            do i=this%pg%imino_,this%pg%imaxo_
               call distance_init(i,j,k,this%dist(i,j,k),this%tag(i,j,k))
            end do
         end do
      end do
   end subroutine build
   
end module fmm_class