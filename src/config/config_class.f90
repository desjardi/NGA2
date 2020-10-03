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
      logical :: xper,yper,zper                                !< Periodicity in x/y/z
      integer :: nbound                                        !< Number of boundary conditions
      !type(bcond), dimension(:), pointer :: bc                 !< Storage array for boundary conditions
   contains
      procedure :: print=>config_print                         !< Output configuration information to the screen
   end type config
   
   
   !> Declare single-grid config constructor
   interface config
      procedure construct_from_bgrid
   end interface config
   
   
contains
   
   
   !> Single-grid config constructor
   function construct_from_bgrid(grid,grp,per) result(self)
      use parallel, only: parallel_min
      implicit none
      type(config) :: self
      type(bgrid), intent(in) :: grid
      integer, intent(in) :: grp
      logical, dimension(3), intent(in) :: per
      integer :: i,j,k
      integer :: powx,powy,powz
      
      ! Create a parallel grid with the provided group
      self%pgrid=pgrid(grid,grp,per)
      
      ! Store periodicity
      self%xper=per(1); self%yper=per(2); self%zper=per(3)
      
      ! Prepare cell volume
      allocate(self%vol(self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_))
      do k=self%kmino_,self%kmaxo_
         do j=self%jmino_,self%jmaxo_
            do i=self%imino_,self%imaxo_
               self%vol(i,j,k)=self%dx(i)*self%dy(j)*self%dz(k)
            end do
         end do
      end do
      
      ! Prepare cell size and small cell size
      powx=1; if (self%nx.eq.1) powx=0
      powy=1; if (self%ny.eq.1) powy=0
      powz=1; if (self%nz.eq.1) powz=0
      allocate(self%meshsize(self%imino_:self%imaxo_,self%jmino_:self%jmaxo_,self%kmino_:self%kmaxo_))
      do k=self%kmino_,self%kmaxo_
         do j=self%jmino_,self%jmaxo_
            do i=self%imino_,self%imaxo_
               self%meshsize(i,j,k)=(self%dx(i)**powx*self%dy(j)**powy*self%dz(k)**powz)**(1.0_WP/real(powx+powy+powz,WP))
            end do
         end do
      end do
      self%min_meshsize=parallel_min(self%meshsize,self%comm)
      
   end function construct_from_bgrid
   
   
   !> Cheap print of sgconfig info to screen
   subroutine config_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(config), intent(in) :: this
      call this%print
   end subroutine config_print
   
end module config_class
