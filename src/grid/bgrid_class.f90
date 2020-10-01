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
      real(WP), dimension(:), allocatable :: x,y,z  !< Mesh in x/y/z
      ! Global index bounds
      integer :: imin,imax                          !< Full domain index bounds in x
      integer :: jmin,jmax                          !< Full domain index bounds in y
      integer :: kmin,kmax                          !< Full domain index bounds in z
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
      function constructor(nx,ny,nz,name) result(self)
         use monitor, only: die
         implicit none
         type(bgrid) :: self
         integer, intent(in) :: nx,ny,nz
         character(len=*), optional :: name
         ! Check sizes
         if (min(nx,ny,nz).le.0) call die('[bgrid constructor] Grid size has to be larger than zero')
         ! Set sizes
         self%nx=nx; self%ny=ny; self%nz=nz
         ! Set bounds
         self%imin=refindex; self%imax=self%imin+self%nx
         self%jmin=refindex; self%jmax=self%jmin+self%ny
         self%kmin=refindex; self%kmax=self%kmin+self%nz
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
         use parallel, only: amroot
         implicit none
         class(bgrid), intent(in) :: this
         if (amroot) then
            write(output_unit,'("Basic grid ",a)') trim(this%name)
            write(output_unit,'(" >   size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
            write(output_unit,'(" > extent = [",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]")') &
            this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
         end if
      end subroutine bgrid_print
      
   end module bgrid_class
