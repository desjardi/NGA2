!> Definition of a basic grid in NGA
module bgrid_class
   use precision, only: WP
   use string, only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: bgrid
   
   ! Reference index
   integer, parameter :: refindex=1
   
   !> Basic grid type: nx is the number of cells, x represents the nx+1 cell faces
   type :: bgrid
      ! General grid info
      character(len=str_medium) :: name='UNDEF'     !< Grid name (default=UNDEF)
      integer :: nx,ny,nz                           !< Grid size in x/y/z
      logical :: xper,yper,zper                     !< Grid periodicity in x/y/z
      real(WP), dimension(:), allocatable :: x,y,z  !< Mesh in x/y/z
      ! Global index bounds
      integer :: imin,imax,jmin,jmax,kmin,kmax      !< Full domain index bounds in x/y/z
   contains
      procedure :: print=>bgrid_print               !< Output grid to screen
      !> Will need to learn to read and write bgrid!
   end type bgrid
   
   !> Declare basic grid constructor
   interface bgrid
      module procedure constructor
   end interface bgrid
   
contains
   
   !> Constructor for a basic grid object
   function constructor(n,per,name) result(self)
      use monitor, only: die
      implicit none
      type(bgrid) :: self
      integer, dimension(3), intent(in) :: n
      logical, dimension(3), intent(in) :: per
      character(len=*), optional :: name
      ! Check sizes
      if (minval(n).le.0) call die('[bgrid constructor] Grid size has to be larger than zero')
      ! Set sizes, periodicity, and bounds
      self%nx=n(1); self%xper=per(1); self%imin=refindex; self%imax=self%imin+self%nx
      self%ny=n(2); self%yper=per(2); self%jmin=refindex; self%jmax=self%jmin+self%ny
      self%nz=n(3); self%zper=per(3); self%kmin=refindex; self%kmax=self%kmin+self%nz
      ! Allocate x/y/z arrays
      allocate(self%x(self%imin:self%imax+1),self%y(self%jmin:self%jmax+1),self%z(self%kmin:self%kmax+1))
      ! Give it a name if one was provided
      if (present(name)) then
         self%name=trim(adjustl(name))
      end if
   end function constructor
   
   !> Print out a basic grid to screen
   subroutine bgrid_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use parallel, only: amRoot
      implicit none
      class(bgrid), intent(in) :: this
      if (amRoot) then
         write(output_unit,'("Basic grid ",a)') trim(this%name)
         write(output_unit,'(" >   size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
         write(output_unit,'(" > extent = [",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]")') &
         this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
      end if
   end subroutine bgrid_print
   
end module bgrid_class
