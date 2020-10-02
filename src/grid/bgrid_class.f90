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
   
   !> Basic grid type
   !> Contains grid size, index extent, and actual mesh - all global 1D info
   type :: bgrid
      ! Grid name
      character(len=str_medium) :: name='UNDEF'        !< Grid name (default=UNDEF)
      ! Mesh dimensions
      integer :: nx,ny,nz                              !< Grid size in x/y/z
      integer :: imin,imax                             !< Domain index bounds in x
      integer :: jmin,jmax                             !< Domain index bounds in y
      integer :: kmin,kmax                             !< Domain index bounds in z
      ! Mesh dimensions with overlap
      integer :: no                                    !< Grid overlap - this overlap can be large
      integer :: nxo,nyo,nzo                           !< Extended grid size in x/y/z
      integer :: imino,imaxo                           !< Extended domain index bounds in x
      integer :: jmino,jmaxo                           !< Extended domain index bounds in y
      integer :: kmino,kmaxo                           !< Extended domain index bounds in z
      ! Mesh data
      real(WP), dimension(:), allocatable :: x,y,z     !< Vertex mesh in x/y/z
      real(WP), dimension(:), allocatable :: xm,ym,zm  !< Centroid mesh in x/y/z
      ! Cell sizes
      real(WP), dimension(:), allocatable :: dx,dxm    !< Cell sizes in x
      real(WP), dimension(:), allocatable :: dy,dym    !< Cell sizes in y
      real(WP), dimension(:), allocatable :: dz,dzm    !< Cell sizes in z
      ! Inverse cell sizes
      real(WP), dimension(:), allocatable :: dxi,dxmi  !< Inverse of cell sizes in x
      real(WP), dimension(:), allocatable :: dyi,dymi  !< Inverse of cell sizes in y
      real(WP), dimension(:), allocatable :: dzi,dzmi  !< Inverse of cell sizes in z
      ! Total length
      real(WP) :: xL,yL,zL                             !< Total domain size
      real(WP) :: vol_total                            !< Total domain volume
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
   function constructor(no,nx,x,ny,y,nz,z,name) result(self)
      use monitor, only: die
      implicit none
      
      type(bgrid) :: self
      
      integer,  intent(in) :: no,nx,ny,nz
      real(WP), intent(in), dimension(nx+1) :: x
      real(WP), intent(in), dimension(ny+1) :: y
      real(WP), intent(in), dimension(nz+1) :: z
      character(len=*), optional :: name
      
      integer :: i,j,k
      
      ! Check sizes
      if (min(nx,ny,nz).le.0) call die('[bgrid constructor] Grid size has to be larger than zero')
      
      ! Set overlap size if it is appropriate
      if (no.lt.0) call die('[bgrid constructor] Overlap size cannot be negative')
      self%no=no
      
      ! Set index sizes and bounds
      self%nx=nx; self%imin=refindex; self%imax=self%imin+self%nx
      self%ny=ny; self%jmin=refindex; self%jmax=self%jmin+self%ny
      self%nz=nz; self%kmin=refindex; self%kmax=self%kmin+self%nz
      
      ! Set extended index sizes and bounds
      self%nxo=self%nx+2*self%no; self%imino=self%imin-self%no; self%imaxo=self%imax+self%no
      self%nyo=self%ny+2*self%no; self%jmino=self%jmin-self%no; self%jmaxo=self%jmax+self%no
      self%nzo=self%nz+2*self%no; self%kmino=self%kmin-self%no; self%kmaxo=self%kmax+self%no
      
      ! Allocate mesh positions based on extended sizes
      allocate(self%x(self%imino:self%imaxo+1)); allocate(self%xm(self%imino:self%imaxo))
      allocate(self%y(self%jmino:self%jmaxo+1)); allocate(self%ym(self%jmino:self%jmaxo))
      allocate(self%z(self%kmino:self%kmaxo+1)); allocate(self%zm(self%kmino:self%kmaxo))
      
      ! Store mesh position
      self%x(imin:imax+1)=x
      self%y(jmin:jmax+1)=y
      self%z(kmin:kmax+1)=z
      
      ! Extend mesh over the overlapping region - this is naive but sufficient
      do i=self%imin-1,self%imino,-1
         self%x(i)=self%x(i+1)-(self%x(self%imin+1)-self%x(self%imin))
      end do
      do i=self%imax+2,self%imaxo+1
         self%x(i)=self%x(i-1)+(self%x(self%imax+1)-self%x(self%imax))
      end do
      do j=self%jmin-1,self%jmino,-1
         self%y(j)=self%y(j+1)-(self%y(self%jmin+1)-self%y(self%jmin))
      end do
      do j=self%jmax+2,self%jmaxo+1
         self%y(j)=self%y(j-1)+(self%y(self%jmax+1)-self%y(self%jmax))
      end do
      do k=self%kmin-1,self%kmino,-1
         self%z(k)=self%z(k+1)-(self%z(self%kmin+1)-self%z(self%kmin))
      end do
      do k=self%kmax+2,self%kmaxo+1
         self%z(k)=self%z(k-1)+(self%z(self%kmax+1)-self%z(self%kmax))
      end do
      
      ! Create centroid mesh
      do i=self%imino,self%imaxo
         self%xm(i)=0.5_WP*(self%x(i)+self%x(i+1))
      end do
      do j=self%jmino,self%jmaxo
         self%ym(j)=0.5_WP*(self%y(j)+self%y(j+1))
      end do
      do k=self%kmino,self%kmaxo
         self%zm(k)=0.5_WP*(self%z(k)+self%z(k+1))
      end do
      
      ! Allocate cell sizes and inverse cell sizes
      allocate(self%dx (self%imino:self%imaxo)); allocate(self%dxm (self%imino+1:self%imaxo))
      allocate(self%dy (self%jmino:self%jmaxo)); allocate(self%dym (self%jmino+1:self%jmaxo))
      allocate(self%dz (self%kmino:self%kmaxo)); allocate(self%dzm (self%kmino+1:self%kmaxo))
      allocate(self%dxi(self%imino:self%imaxo)); allocate(self%dxmi(self%imino+1:self%imaxo))
      allocate(self%dyi(self%jmino:self%jmaxo)); allocate(self%dymi(self%jmino+1:self%jmaxo))
      allocate(self%dzi(self%kmino:self%kmaxo)); allocate(self%dzmi(self%kmino+1:self%kmaxo))
      
      ! Prepare the cell size
      do i=self%imino,self%imaxo
         self%dx(i)=self%x(i+1)-self%x(i)
      end do
      do j=self%jmino,self%jmaxo
         self%dy(j)=self%y(j+1)-self%y(j)
      end do
      do k=self%kmino,self%kmaxo
         self%dz(k)=self%z(k+1)-self%z(k)
      end do
      
      ! Check that the mesh is proper
      if (minval(self%dx).le.0.0_WP) call die('[bgrid constructor] Mesh in x is not monotonically increasing')
      if (minval(self%dy).le.0.0_WP) call die('[bgrid constructor] Mesh in y is not monotonically increasing')
      if (minval(self%dz).le.0.0_WP) call die('[bgrid constructor] Mesh in z is not monotonically increasing')
      
      ! Prepare the inverse of cell size
      do i=self%imino,self%imaxo
         self%dxi(i)=1.0_WP/self%dx(i)
      end do
      do j=self%jmino,self%jmaxo
         self%dyi(j)=1.0_WP/self%dy(j)
      end do
      do k=self%kmino,self%kmaxo
         self%dzi(k)=1.0_WP/self%dz(k)
      end do
      
      ! Prepare the centroid distance
      do i=self%imino+1,self%imaxo
         self%dxm(i)=self%xm(i)-self%xm(i-1)
      end do
      do j=self%jmino+1,self%jmaxo
         self%dym(j)=self%ym(j)-self%ym(j-1)
      end do
      do k=self%kmino+1,self%kmaxo
         self%dzm(k)=self%zm(k)-self%zm(k-1)
      end do
      
      ! Prepare the inverse of centroid distance
      do i=self%imino+1,self%imaxo
         self%dxmi(i)=1.0_WP/self%dxm(i)
      end do
      do j=self%jmino+1,self%jmaxo
         self%dymi(j)=1.0_WP/self%dym(j)
      end do
      do k=self%kmino+1,self%kmaxo
         self%dzmi(k)=1.0_WP/self%dzm(k)
      end do
      
      ! Compute total length and total volume
      self%xL=self%x(imax+1)-self%x(imin)
      self%yL=self%y(jmax+1)-self%y(jmin)
      self%zL=self%z(kmax+1)-self%z(kmin)
      self%vol_total=self%xL*self%yL*self%zL
      
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
         write(output_unit,'(" > overlap = ",i0)') this%no
         write(output_unit,'(" >    size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
         write(output_unit,'(" >  extent = [",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]")') &
         this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
         write(output_unit,'(" >  volume = ",es12.6)') this%vol_total
      end if
   end subroutine bgrid_print
   
end module bgrid_class
