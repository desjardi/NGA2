!> Single-grid config concept is defined here: this is a partitioned grid
!> as well as geometry (i.e., masks, Gib, or similar)
module config_class
   use precision,   only: WP
   use pgrid_class, only: pgrid
   use mpi_f08,     only: MPI_Group
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
      procedure :: write=>config_write                         !< Write out config files: grid and geometry
      procedure :: maskupdate                                  !< Takes in a mask array and updates the config consequently
      procedure, private :: prep=>config_prep                  !< Finish preparing config after the partitioned grid is loaded
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
      type(MPI_Group), intent(in) :: grp
      integer, dimension(3), intent(in) :: decomp
      ! Create a partitioned grid with the provided group and decomposition
      self%pgrid=pgrid(grid,grp,decomp)
      ! Finish preparing the config
      call self%prep
   end function construct_from_sgrid
   
   
   !> Single-grid config constructor from NGA grid file
   function construct_from_file(grp,decomp,no,fgrid,fgeom) result(self)
      implicit none
      type(config) :: self
      type(MPI_Group), intent(in) :: grp
      integer, dimension(3), intent(in) :: decomp
      integer, intent(in) :: no
      character(*), intent(in) :: fgrid
      character(*), intent(in) :: fgeom
      ! Create a partitioned grid with the provided group and decomposition
      self%pgrid=pgrid(no,fgrid,grp,decomp)
      ! Finish preparing the config
      call self%prep
      ! Read in the mask info and update the masks
      read_mask: block
         use datafile_class, only: datafile
         type(datafile) :: geomfile
         geomfile=datafile(self,fgeom)
         call geomfile%pullvar('mask',self%mask)
         call self%maskupdate()
      end block read_mask
   end function construct_from_file
   
   
   !> Prepare a config once the pgrid is set
   subroutine config_prep(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN
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
   
   
   !> Updates mask info for the config
   subroutine maskupdate(this)
      implicit none
      class(config) :: this
      ! Communicate mask data
      call this%sync(this%mask)
      ! Apply Neumann in all non-periodic directions
      
   end subroutine maskupdate
   
   
   !> Cheap print of config info to screen
   subroutine config_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(config), intent(in) :: this
      call this%print
   end subroutine config_print
   
   
   !> Output a config object to a file: both grid and geom
   !> Expect a parallel call here
   subroutine config_write(this,file)
      implicit none
      class(config) :: this
      character(len=*), intent(in) :: file
      ! Root process writes out the grid
      if (this%amRoot) call this%sgrid%write(trim(adjustl(file))//'.grid')
      ! All processes write out the geometry
      write_mask: block
         use datafile_class, only: datafile
         type(datafile) :: geomfile
         geomfile=datafile(this,trim(this%name),0,2)
         geomfile%varname(1)='mask2'
         call geomfile%pushvar('mask2',this%mask)
         geomfile%varname(2)='mask'
         call geomfile%pushvar('mask',this%mask)
         call geomfile%write(trim(adjustl(file))//'.geom')
      end block write_mask
   end subroutine config_write
   
   
end module config_class
