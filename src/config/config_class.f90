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
      procedure :: set_scalar                                  !< Subroutine that extrapolates a provided scalar value at a point to a pgrid field
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
      real(WP) :: my_int
      my_int=0.0_WP
      do k=this%kmin_,this%kmax_
         do j=this%jmin_,this%jmax_
            do i=this%imin_,this%imax_
               my_int=my_int+this%vol(i,j,k)*this%VF(i,j,k)*A(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_int,integral,1,MPI_REAL_WP,MPI_SUM,this%comm,ierr)
   end subroutine integrate


   !> Subroutine that performs trilinear extrapolation of provided value Sp
   !> at provided position pos near cell i0,j0,k0 to provided scalar field S
   subroutine set_scalar(this,Sp,pos,i0,j0,k0,S)
      implicit none
      class(config), intent(in) :: this
      real(WP), intent(in) :: Sp
      real(WP), dimension(3), intent(in) :: pos
      integer, intent(in) :: i0,j0,k0
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: S     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ni,nj,nk
      real(WP) :: wx1,wy1,wz1
      real(WP) :: wx2,wy2,wz2
      real(WP), dimension(2,2,2) :: ww
      ! Find right i index
      i=max(min(this%imaxo_-1,i0),this%imino_)
      do while (pos(1)-this%xm(i  ).lt.0.0_WP.and.i  .gt.this%imino_); i=i-1; end do
      do while (pos(1)-this%xm(i+1).ge.0.0_WP.and.i+1.lt.this%imaxo_); i=i+1; end do
      ! Find right j index
      j=max(min(this%jmaxo_-1,j0),this%jmino_)
      do while (pos(2)-this%ym(j  ).lt.0.0_WP.and.j  .gt.this%jmino_); j=j-1; end do
      do while (pos(2)-this%ym(j+1).ge.0.0_WP.and.j+1.lt.this%jmaxo_); j=j+1; end do
      ! Find right k index
      k=max(min(this%kmaxo_-1,k0),this%kmino_)
      do while (pos(3)-this%zm(k  ).lt.0.0_WP.and.k  .gt.this%kmino_); k=k-1; end do
      do while (pos(3)-this%zm(k+1).ge.0.0_WP.and.k+1.lt.this%kmaxo_); k=k+1; end do
      ! For exact conservation, all information needs to be placed *inside* the domain
      if (.not.this%xper) i=max(this%imin,min(this%imax-1,i))
      if (.not.this%yper) j=max(this%jmin,min(this%jmax-1,j))
      if (.not.this%zper) k=max(this%kmin,min(this%kmax-1,k))
      ! Prepare tri-linear extrapolation coefficients
      wx1=(pos(1)-this%xm(i))/(this%xm(i+1)-this%xm(i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%ym(j))/(this%ym(j+1)-this%ym(j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%zm(k))/(this%zm(k+1)-this%zm(k)); wz2=1.0_WP-wz1
      ! Combine into array for easy rescaling
      ww(1,1,1)=wx2*wy2*wz2
      ww(2,1,1)=wx1*wy2*wz2
      ww(1,2,1)=wx2*wy1*wz2
      ww(2,2,1)=wx1*wy1*wz2
      ww(1,1,2)=wx2*wy2*wz1
      ww(2,1,2)=wx1*wy2*wz1
      ww(1,2,2)=wx2*wy1*wz1
      ww(2,2,2)=wx1*wy1*wz1
      ! Apply Neumann condition at walls
      do nk=0,1
         do nj=0,1
            do ni=0,1
               if (this%VF(i+ni,j+nj,k+nk).eq.0.0_WP) ww(1+ni,1+nj,1+nk)=0.0_WP
            end do
         end do
      end do
      ww=ww/sum(ww)
      ! Tri-linear extrapolation of Sp onto S
      S(i:i+1,j:j+1,k:k+1)=S(i:i+1,j:j+1,k:k+1)+ww*Sp
   end subroutine set_scalar
   
   
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
