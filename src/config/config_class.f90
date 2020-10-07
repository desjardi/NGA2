!> Single-grid config concept is defined here: this is a partitioned grid
!> as well as geometry (i.e., masks, Gib, or similar)
module config_class
   use precision,      only: WP
   use pgrid_class,    only: pgrid
   use datafile_class, only: datafile,datafile_from_file
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: config
   
   !> Config object definition as an extension of pgrid
   type, extends(pgrid) :: config
      ! Some more metrics
      real(WP), dimension(:,:,:), allocatable :: vol           !< Local cell volume
      real(WP), dimension(:,:,:), allocatable :: meshsize      !< Local effective cell size
      real(WP) :: min_meshsize                                 !< Global minimum mesh size
      ! Geometry
      real(WP), dimension(:,:,:), allocatable :: mask          !< Masking info (mask=0 is fluid, mask=1 is wall)
   contains
      procedure :: print=>config_print                         !< Output configuration information to the screen
      procedure, private :: prep=>config_prep                  !< Finish preparing config after the partitioned grid is loaded
      procedure :: write=>config_write                         !< Write out config files: grid and geometry
   end type config
   
   
   !> Declare single-grid config constructor
   interface config
      procedure construct_from_sgrid
      procedure construct_from_file
   end interface config
   
   
contains
   
   
   !> Single-grid config constructor from a serial grid
   function construct_from_sgrid(grp,decomp,grid) result(self)
      use sgrid_class, only: sgrid
      use string,      only: str_medium
      implicit none
      type(config) :: self
      type(sgrid), intent(in) :: grid
      integer, intent(in) :: grp
      integer, dimension(3), intent(in) :: decomp
      ! Create a partitioned grid with the provided group and decomposition
      self%pgrid=pgrid(grid,grp,decomp)
      ! Finish preparing the config
      call self%prep
   end function construct_from_sgrid
   
   
   !> Single-grid config constructor from NGA grid file
   function construct_from_file(grp,decomp,no,fgrid) result(self)
      implicit none
      type(config) :: self
      character(*), intent(in) :: fgrid
      integer, intent(in) :: grp,no
      integer, dimension(3), intent(in) :: decomp
      ! Create a partitioned grid with the provided group and decomposition
      self%pgrid=pgrid(no,fgrid,grp,decomp)
      ! Finish preparing the config
      call self%prep
   end function construct_from_file
   
   
   !> Prepare a config once the pgrid is set
   subroutine config_prep(this)
      use mpi
      use parallel, only: MPI_REAL_WP
      implicit none
      class(config) :: this
      integer :: i,j,k,ierr
      integer :: powx,powy,powz
      
      ! Prepare cell volume
      allocate(this%vol(this%imino_:this%imaxo_,this%jmino_:this%jmaxo_,this%kmino_:this%kmaxo_))
      do k=this%kmino_,this%kmaxo_
         do j=this%jmino_,this%jmaxo_
            do i=this%imino_,this%imaxo_
               this%vol(i,j,k)=this%dx(i)*this%dy(j)*this%dz(k)
            end do
         end do
      end do
      
      ! Prepare cell size and small cell size
      powx=1; if (this%nx.eq.1) powx=0
      powy=1; if (this%ny.eq.1) powy=0
      powz=1; if (this%nz.eq.1) powz=0
      allocate(this%meshsize(this%imino_:this%imaxo_,this%jmino_:this%jmaxo_,this%kmino_:this%kmaxo_))
      do k=this%kmino_,this%kmaxo_
         do j=this%jmino_,this%jmaxo_
            do i=this%imino_,this%imaxo_
               this%meshsize(i,j,k)=(this%dx(i)**powx*this%dy(j)**powy*this%dz(k)**powz)**(1.0_WP/real(powx+powy+powz,WP))
            end do
         end do
      end do
      call MPI_ALLREDUCE(minval(this%meshsize(this%imin_:this%imax_,this%jmin_:this%jmax_,this%kmin_:this%kmax_)),this%min_meshsize,1,MPI_REAL_WP,MPI_MIN,this%comm,ierr)
      
      ! Allocate wall geometry - assume all fluid until told otherwise
      allocate(this%mask(this%imino_:this%imaxo_,this%jmino_:this%jmaxo_,this%kmino_:this%kmaxo_))
      this%mask=0.0_WP
      
   end subroutine config_prep
   
   
   !> Cheap print of config info to screen
   subroutine config_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(config), intent(in) :: this
      call this%print
   end subroutine config_print
   
   
   !> Output a config object to a file
   subroutine config_write(this,file)
      implicit none
      class(config) :: this
      character(len=*), intent(in) :: file
      ! Write out the grid
      call this%pgrid%write(trim(adjustl(file))//".grid")
   end subroutine config_write
   
   
end module config_class
