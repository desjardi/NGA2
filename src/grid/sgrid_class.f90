!> Definition of a serial grid in NGA - this is essentially a OOP remaster of NGA's original config concept
module sgrid_class
   use precision, only: WP
   use string, only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: sgrid
   
   ! Reference index
   integer, parameter :: refindex=1
   
   ! Coordinate systems
   integer, parameter, public :: cartesian  =0
   integer, parameter, public :: cylindrical=1
   integer, parameter, public :: spherical  =2
   
   !> Serial grid type:
   !> Contains grid size, index extent, overlap and periodicity,
   !> as well as the actual 1D mesh along with 1D metrics
   type :: sgrid
      ! Grid name
      character(len=str_medium) :: name='UNNAMED_GRID' !< Grid name (default=UNNAMED_GRID)
      ! Coordinate system
      integer :: coordsys                              !< Coordinate system
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
      ! Mesh uniformity
      logical :: uniform_x,uniform_y,uniform_z         !< Mesh uniformity in x/y/z
      ! Grid periodicity
      logical :: xper,yper,zper                        !< Periodicity in x/y/z
   contains
      procedure :: print=>sgrid_print                  !< Output grid to screen
      procedure :: log  =>sgrid_log                    !< Output grid to log file
      procedure :: write=>sgrid_write                  !< Output grid to grid file
   end type sgrid
   
   !> Declare basic grid constructor
   interface sgrid
      module procedure construct_from_file
      module procedure construct_from_args
   end interface sgrid
   
contains
   
   
   !> Constructor for a serial grid object based on NGA2 new grid file
   function construct_from_file(no,file,name) result(self)
      use messager, only: die
      implicit none
      type(sgrid) :: self
      integer, intent(in) :: no
      character(len=*), intent(in) :: file
      character(len=*), optional :: name
      integer :: iunit,ierr
      character(len=str_medium) :: simu_name
      integer :: nx,ny,nz,coord
      real(WP), dimension(:), allocatable :: x
      real(WP), dimension(:), allocatable :: y
      real(WP), dimension(:), allocatable :: z
      logical :: xper,yper,zper
      
      ! Open the file, read grid info, close the file
      open(newunit=iunit,file=trim(adjustl(file)),form='unformatted',status='old',access='stream',iostat=ierr)
      if (ierr.ne.0) call die('[sgrid constructor] Could not open file: '//trim(file))
      read(iunit) simu_name,coord,xper,yper,zper,nx,ny,nz
      allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1)); read(iunit) x,y,z
      close(iunit)
      
      ! Update name if another one was provided
      if (present(name)) simu_name=trim(adjustl(name))
      
      ! Use arg-based constructor now that all info are known
      self=construct_from_args(coord,no,x,y,z,xper,yper,zper,trim(simu_name))
      
      ! Deallocate 1D arrays
      deallocate(x,y,z)
      
   end function construct_from_file
   
   
   !> Constructor for a serial grid object
   function construct_from_args(coord,no,x,y,z,xper,yper,zper,name) result(self)
      use messager, only: die
      use param,    only: verbose
      implicit none
      type(sgrid) :: self
      integer,  intent(in) :: coord
      integer,  intent(in) :: no
      real(WP), dimension(:), intent(in) :: x
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(:), intent(in) :: z
      logical, intent(in) :: xper,yper,zper
      character(len=*), optional :: name
      
      integer :: nx,ny,nz
      integer :: i,j,k
      
      ! Obtain sizes from arrays and check them
      nx=size(x)-1; ny=size(y)-1; nz=size(z)-1
      if (min(nx,ny,nz).le.0) call die('[sgrid constructor] Grid size has to be larger than zero')
      
      ! Check coordinate system
      select case (coord)
      case (cartesian)
         self%coordsys=cartesian
      case (cylindrical)
         self%coordsys=cylindrical; call die('[sgrid constructor] Cylindrical coordinates are not implemented yet')
      case (spherical)
         self%coordsys=spherical; call die('[sgrid constructor] Spherical coordinates are not implemented yet')
      case default
         call die('[sgrid constructor] The coordinate system is unknown')
      end select
      
      ! Set overlap size if it is appropriate
      if (no.lt.0) call die('[sgrid constructor] Overlap size cannot be negative')
      self%no=no
      
      ! Store periodicity
      self%xper=xper; self%yper=yper; self%zper=zper
      
      ! Give it a name if one was provided
      if (present(name)) then
         self%name=trim(adjustl(name))
      end if
      
      ! Set index sizes and bounds
      self%nx=nx; self%imin=refindex; self%imax=self%imin+self%nx-1
      self%ny=ny; self%jmin=refindex; self%jmax=self%jmin+self%ny-1
      self%nz=nz; self%kmin=refindex; self%kmax=self%kmin+self%nz-1
      
      ! Set extended index sizes and bounds
      self%nxo=self%nx+2*self%no; self%imino=self%imin-self%no; self%imaxo=self%imax+self%no
      self%nyo=self%ny+2*self%no; self%jmino=self%jmin-self%no; self%jmaxo=self%jmax+self%no
      self%nzo=self%nz+2*self%no; self%kmino=self%kmin-self%no; self%kmaxo=self%kmax+self%no
      
      ! Allocate mesh positions based on extended sizes
      allocate(self%x(self%imino:self%imaxo+1)); allocate(self%xm(self%imino:self%imaxo))
      allocate(self%y(self%jmino:self%jmaxo+1)); allocate(self%ym(self%jmino:self%jmaxo))
      allocate(self%z(self%kmino:self%kmaxo+1)); allocate(self%zm(self%kmino:self%kmaxo))
      
      ! Store mesh position
      self%x(self%imin:self%imax+1)=x
      self%y(self%jmin:self%jmax+1)=y
      self%z(self%kmin:self%kmax+1)=z
      
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
      if (minval(self%dx).le.0.0_WP) call die('[sgrid constructor] Mesh in x is not monotonically increasing')
      if (minval(self%dy).le.0.0_WP) call die('[sgrid constructor] Mesh in y is not monotonically increasing')
      if (minval(self%dz).le.0.0_WP) call die('[sgrid constructor] Mesh in z is not monotonically increasing')
      
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
      self%xL=self%x(self%imax+1)-self%x(self%imin)
      self%yL=self%y(self%jmax+1)-self%y(self%jmin)
      self%zL=self%z(self%kmax+1)-self%z(self%kmin)
      self%vol_total=self%xL*self%yL*self%zL
      
      ! Check mesh uniformity - using xL*epsilon to test equality here
      self%uniform_x=.false.; if (abs(maxval(self%dx)-minval(self%dx)).lt.self%xL*10.0_WP*epsilon(1.0_WP)) self%uniform_x=.true.
      self%uniform_y=.false.; if (abs(maxval(self%dy)-minval(self%dy)).lt.self%yL*10.0_WP*epsilon(1.0_WP)) self%uniform_y=.true.
      self%uniform_z=.false.; if (abs(maxval(self%dz)-minval(self%dz)).lt.self%zL*10.0_WP*epsilon(1.0_WP)) self%uniform_z=.true.
      
      ! If verbose run, log and or print grid
      if (verbose.gt.2) call self%log
      if (verbose.gt.3) call self%print
      
   end function construct_from_args
   
   
   !> Print out a serial grid to screen
   subroutine sgrid_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use parallel, only: amRoot
      implicit none
      class(sgrid), intent(in) :: this
      if (amRoot) then
         select case (this%coordsys)
         case (cartesian);   write(output_unit,'("Serial cartesian grid [",a,"]")') trim(this%name)
         case (cylindrical); write(output_unit,'("Serial cylindrical grid [",a,"]")') trim(this%name)
         case (spherical);   write(output_unit,'("Serial spherical grid [",a,"]")') trim(this%name)
         end select
         write(output_unit,'(" >   overlap = ",i0)') this%no
         write(output_unit,'(" >      size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
         write(output_unit,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
         write(output_unit,'(" >   uniform = ",l1,"x",l1,"x",l1)') this%uniform_x,this%uniform_y,this%uniform_z
         write(output_unit,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper
         write(output_unit,'(" >    volume = ",es12.5)') this%vol_total
      end if
   end subroutine sgrid_print
   
   
   !> Print out a serial grid to log
   subroutine sgrid_log(this)
      use parallel, only: amRoot
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(sgrid), intent(in) :: this
      character(len=str_long) :: message
      if (amRoot) then
         select case (this%coordsys)
         case (cartesian);   write(message,'("Serial cartesian grid [",a,"]")') trim(this%name); call log(message)
         case (cylindrical); write(message,'("Serial cylindrical grid [",a,"]")') trim(this%name); call log(message)
         case (spherical);   write(message,'("Serial spherical grid [",a,"]")') trim(this%name); call log(message)
         end select
         write(message,'(" >   overlap = ",i0)') this%no; call log(message)
         write(message,'(" >      size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz; call log(message)
         write(message,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1); call log(message)
         write(message,'(" >   uniform = ",l1,"x",l1,"x",l1)') this%uniform_x,this%uniform_y,this%uniform_z; call log(message)
         write(message,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper; call log(message)
         write(message,'(" >    volume = ",es12.5)') this%vol_total; call log(message)
      end if
   end subroutine sgrid_log
   
   
   !> Output a serial grid object to a file
   subroutine sgrid_write(this,file)
      use messager, only: die,log
      use param,    only: verbose
      implicit none
      class(sgrid), intent(in) :: this
      character(len=*), intent(in) :: file
      integer :: iunit,ierr
      
      ! Open the file, write grid info, close the file
      open(newunit=iunit,file=trim(adjustl(file)),form='unformatted',status='replace',access='stream',iostat=ierr)
      if (ierr.ne.0) call die('[sgrid write] Could not open file: '//trim(adjustl(file)))
      write(iunit) this%name,this%coordsys,this%xper,this%yper,this%zper,this%nx,this%ny,this%nz,this%x(this%imin:this%imax+1),this%y(this%jmin:this%jmax+1),this%z(this%kmin:this%kmax+1)
      close(iunit)
      
      ! Log the action
      if (verbose.gt.0) call log('Grid file written: '//trim(adjustl(file)))
      
   end subroutine sgrid_write
   
end module sgrid_class
