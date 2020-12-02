!> Single-grid config concept is defined here: this is a partitioned grid
!> as well as geometry (i.e., masks, Gib, or similar)
module config_class
   use precision,   only: WP
   use pgrid_class, only: pgrid
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
      real(WP), dimension(:,:,:), allocatable :: VF            !< Volume fraction info (VF=1 is fluid, VF=0 is wall)
      
   contains
      procedure :: print=>config_print                         !< Output configuration information to the screen
      procedure :: write=>config_write                         !< Write out config files: grid and geometry
      procedure, private :: prep=>config_prep                  !< Finish preparing config after the partitioned grid is loaded
      procedure :: VF_extend                                   !< Extend VF array into the non-periodic domain overlaps
      procedure :: integrate                                   !< Integrate a variable on config while accounting for VF
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
      use mpi_f08,     only: MPI_Group
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
      use mpi_f08, only: MPI_Group
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
      ! Read in the VF info
      read_VF: block
         use datafile_class, only: datafile
         type(datafile) :: geomfile
         ! Access the file
         geomfile=datafile(self,fgeom)
         ! Get the VF array
         call geomfile%pullvar('VF',self%VF)
         ! Sync up the VF array
         call self%sync(self%VF)
         ! Perform an extension in the overlap
         call self%VF_extend()
      end block read_VF
   end function construct_from_file
   
   
   !> Prepare a config once the pgrid is set
   subroutine config_prep(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(config), intent(inout) :: this
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
      allocate(this%VF(this%imino_:this%imaxo_,this%jmino_:this%jmaxo_,this%kmino_:this%kmaxo_)); this%VF=1.0_WP
      
   end subroutine config_prep
   
   
   !> Extend VF array into the non-periodic domain overlaps
   subroutine VF_extend(this)
      implicit none
      class(config), intent(inout) :: this
      integer :: i,j,k
      if (.not.this%xper) then
         if (this%iproc.eq.1) then
            do i=this%imino,this%imin-1
               this%VF(i,:,:)=this%VF(this%imin,:,:)
            end do
         else if (this%iproc.eq.this%npx) then
            do i=this%imax+1,this%imaxo
               this%VF(i,:,:)=this%VF(this%imax,:,:)
            end do
         end if
      end if
      if (.not.this%yper) then
         if (this%jproc.eq.1) then
            do j=this%jmino,this%jmin-1
               this%VF(:,j,:)=this%VF(:,this%jmin,:)
            end do
         else if (this%jproc.eq.this%npy) then
            do j=this%jmax+1,this%jmaxo
               this%VF(:,j,:)=this%VF(:,this%jmax,:)
            end do
         end if
      end if
      if (.not.this%zper) then
         if (this%kproc.eq.1) then
            do k=this%kmino,this%kmin-1
               this%VF(:,:,k)=this%VF(:,:,this%kmin)
            end do
         else if (this%kproc.eq.this%npz) then
            do j=this%kmax+1,this%kmaxo
               this%VF(:,:,k)=this%VF(:,:,this%kmax)
            end do
         end if
      end if
   end subroutine VF_extend
   
   
   !> Calculate the integral of a field on the config while accounting for VF
   subroutine integrate(this,A,integral)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(config), intent(inout) :: this
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(in) :: A      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: integral
      integer :: i,j,k,ierr
      real(WP) :: my_int,my_vol,volume
      my_int=0.0_WP
      my_vol=0.0_WP
      do k=this%kmin_,this%kmax_
         do j=this%jmin_,this%jmax_
            do i=this%imin_,this%imax_
               my_int=my_int+this%vol(i,j,k)*this%VF(i,j,k)*A(i,j,k)
               my_vol=my_vol+this%vol(i,j,k)*this%VF(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_int,integral,1,MPI_REAL_WP,MPI_SUM,this%comm,ierr)
      call MPI_ALLREDUCE(my_vol,volume  ,1,MPI_REAL_WP,MPI_SUM,this%comm,ierr)
      integral=integral/volume
   end subroutine integrate
   
   
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
      class(config), intent(in) :: this
      character(len=*), intent(in) :: file
      ! Root process writes out the grid
      if (this%amRoot) call this%pgrid%write(trim(adjustl(file))//'.grid')
      ! All processes write out the geometry
      write_VF: block
         use datafile_class, only: datafile
         type(datafile) :: geomfile
         geomfile=datafile(this,trim(this%name),0,1)
         geomfile%varname(1)='VF'
         call geomfile%pushvar('VF',this%VF)
         call geomfile%write(trim(adjustl(file))//'.geom')
      end block write_VF
   end subroutine config_write
   
   
end module config_class
