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
      real(WP) :: fluid_vol                                    !< Total fluid volume in the config
      
   contains
      procedure :: print=>config_print                         !< Output configuration information to the screen
      procedure :: write=>config_write                         !< Write out config files: grid and geometry
      procedure, private :: prep=>config_prep                  !< Finish preparing config after the partitioned grid is loaded
      procedure :: calc_fluid_vol                              !< Compute fluid_vol after VF has been set
      procedure :: VF_extend                                   !< Extend VF array into the non-periodic domain overlaps
      procedure :: integrate                                   !< Integrate a variable on config while accounting for VF
      procedure :: integrate_without_VF                        !< Integrate a variable on config while ignoring VF
      procedure :: set_scalar                                  !< Subroutine that extrapolates a provided scalar value at a point to a pgrid field
      procedure :: get_velocity                                !< Function that interpolates a provided velocity field staggered on the pgrid to a point
      procedure :: get_scalar                                  !< Function that interpolates a provided scalar field centered on the pgrid to a point
      procedure :: maximum                                     !< Find global max of variable on config for VF>0
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
         ! Update fluid_vol
         call self%calc_fluid_vol()
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
      this%fluid_vol=this%vol_total
      
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
   

   !> Calculate the total fluid volume by integrating VF
   subroutine calc_fluid_vol(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(config), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: buf
      buf=0.0_WP
      do k=this%kmin_,this%kmax_
         do j=this%jmin_,this%jmax_
            do i=this%imin_,this%imax_
               buf=buf+this%vol(i,j,k)*this%VF(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(buf,this%fluid_vol,1,MPI_REAL_WP,MPI_SUM,this%comm,ierr)
   end subroutine calc_fluid_vol
   
   
   !> Find global max of variable on config for VF>0
   subroutine maximum(this,A,globalmax)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(config), intent(inout) :: this
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(in) :: A      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: globalmax
      integer :: i,j,k,ierr
      real(WP) :: my_max
      my_max=0.0_WP
      do k=this%kmin_,this%kmax_
         do j=this%jmin_,this%jmax_
            do i=this%imin_,this%imax_
               if (this%VF(i,j,k).gt.0.0_WP) my_max=max(my_max,A(i,j,k))
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_max,globalmax,1,MPI_REAL_WP,MPI_MAX,this%comm,ierr)
   end subroutine maximum

   
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


   !> Calculate the integral of a field on the config while ignoring VF
   subroutine integrate_without_VF(this,A,integral)
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
               my_int=my_int+this%vol(i,j,k)*A(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_int,integral,1,MPI_REAL_WP,MPI_SUM,this%comm,ierr)
   end subroutine integrate_without_VF
   
   
   !> Subroutine that performs trilinear interpolation of provided scalar field S
   !> to provided position pos near cell i0,j0,k0, stores it in Sp 
   function get_scalar(this,pos,i0,j0,k0,S,bc) result(Sp)
      use messager, only: die
      implicit none
      class(config), intent(in) :: this
      real(WP) :: Sp
      real(WP), dimension(3), intent(in) :: pos
      integer, intent(in) :: i0,j0,k0
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: S     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      character(len=1), intent(in) :: bc    !< Supports n for Neumann, d for Dirichlet, 0 for zero Dirichlet
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
      ! Apply bc
      select case (bc)
      case ('n','N')
         ! Apply Neumann condition at walls
         do nk=0,1
            do nj=0,1
               do ni=0,1
                  if (this%VF(i+ni,j+nj,k+nk).eq.0.0_WP) ww(1+ni,1+nj,1+nk)=0.0_WP
               end do
            end do
         end do
         ww=ww/(sum(ww)+epsilon(1.0_WP))
      case ('0')
         ! Apply zero Dirichlet at walls
         do nk=0,1
            do nj=0,1
               do ni=0,1
                  if (this%VF(i+ni,j+nj,k+nk).eq.0.0_WP) ww(1+ni,1+nj,1+nk)=0.0_WP
               end do
            end do
         end do
      case ('d','D')
         ! Apply Dirichlet at walls
      case default
         ! Not recognized
         call die('[cfg get_scalar] Boundary condition not recognized')
      end select
      ! Tri-linear interpolation of S
      Sp=sum(ww*S(i:i+1,j:j+1,k:k+1))
   end function get_scalar

   
   !> Subroutine that performs trilinear extrapolation of provided value Sp
   !> at provided position pos near cell i0,j0,k0 to provided scalar field S
   subroutine set_scalar(this,Sp,pos,i0,j0,k0,S,bc)
      use messager, only: die
      implicit none
      class(config), intent(in) :: this
      real(WP), intent(in) :: Sp
      real(WP), dimension(3), intent(in) :: pos
      integer, intent(in) :: i0,j0,k0
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(inout) :: S     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      character(len=1), intent(in) :: bc    !< Supports n for Neumann, d for Dirichlet, 0 for zero Dirichlet
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
      !if (.not.this%xper) i=max(this%imin,min(this%imax-1,i))
      !if (.not.this%yper) j=max(this%jmin,min(this%jmax-1,j))
      !if (.not.this%zper) k=max(this%kmin,min(this%kmax-1,k))
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
      ! Apply bc
      select case (bc)
      case ('n','N')
         ! Apply Neumann condition at walls
         do nk=0,1
            do nj=0,1
               do ni=0,1
                  if (this%VF(i+ni,j+nj,k+nk).eq.0.0_WP) ww(1+ni,1+nj,1+nk)=0.0_WP
               end do
            end do
         end do
         ww=ww/(sum(ww)+epsilon(1.0_WP))
      case ('0')
         ! Apply zero Dirichlet at walls
         do nk=0,1
            do nj=0,1
               do ni=0,1
                  if (this%VF(i+ni,j+nj,k+nk).eq.0.0_WP) ww(1+ni,1+nj,1+nk)=0.0_WP
               end do
            end do
         end do
      case ('d','D')
         ! Apply Dirichlet at walls
      case default
         ! Not recognized
         call die('[cfg set_scalar] Boundary condition not recognized')
      end select
      ! Tri-linear extrapolation of Sp onto S
      S(i:i+1,j:j+1,k:k+1)=S(i:i+1,j:j+1,k:k+1)+ww*Sp
   end subroutine set_scalar
   
   
   !> Function that performs trilinear interpolation of provided
   !> velocity U,V,W to provided position pos near cell i0,j0,k0
   function get_velocity(this,pos,i0,j0,k0,U,V,W) result(vel)
      implicit none
      class(config), intent(in) :: this
      real(WP), dimension(3) :: vel
      real(WP), dimension(3), intent(in) :: pos
      integer, intent(in) :: i0,j0,k0
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%imino_:,this%jmino_:,this%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP) :: wx1,wy1,wz1
      real(WP) :: wx2,wy2,wz2
      ! Interpolate U velocity ------------------------------
      ! Find right i index
      i=max(min(this%imaxo_-1,i0),this%imino_)
      do while (pos(1)-this%x (i  ).lt.0.0_WP.and.i  .gt.this%imino_); i=i-1; end do
      do while (pos(1)-this%x (i+1).ge.0.0_WP.and.i+1.lt.this%imaxo_); i=i+1; end do
      ! Find right j index
      j=max(min(this%jmaxo_-1,j0),this%jmino_)
      do while (pos(2)-this%ym(j  ).lt.0.0_WP.and.j  .gt.this%jmino_); j=j-1; end do
      do while (pos(2)-this%ym(j+1).ge.0.0_WP.and.j+1.lt.this%jmaxo_); j=j+1; end do
      ! Find right k index
      k=max(min(this%kmaxo_-1,k0),this%kmino_)
      do while (pos(3)-this%zm(k  ).lt.0.0_WP.and.k  .gt.this%kmino_); k=k-1; end do
      do while (pos(3)-this%zm(k+1).ge.0.0_WP.and.k+1.lt.this%kmaxo_); k=k+1; end do
      ! Prepare tri-linear interpolation coefficients
      wx1=(pos(1)-this%x (i))/(this%x (i+1)-this%x (i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%ym(j))/(this%ym(j+1)-this%ym(j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%zm(k))/(this%zm(k+1)-this%zm(k)); wz2=1.0_WP-wz1
      ! Tri-linear interpolation of U
      vel(1)=wz1*(wy1*(wx1*U(i+1,j+1,k+1)  + &
      &                wx2*U(i  ,j+1,k+1)) + &
      &           wy2*(wx1*U(i+1,j  ,k+1)  + &
      &                wx2*U(i  ,j  ,k+1)))+ &
      &      wz2*(wy1*(wx1*U(i+1,j+1,k  )  + &
      &                wx2*U(i  ,j+1,k  )) + &
      &           wy2*(wx1*U(i+1,j  ,k  )  + &
      &                wx2*U(i  ,j  ,k  )))
      ! Interpolate V velocity ------------------------------
      ! Find right i index
      i=max(min(this%imaxo_-1,i0),this%imino_)
      do while (pos(1)-this%xm(i  ).lt.0.0_WP.and.i  .gt.this%imino_); i=i-1; end do
      do while (pos(1)-this%xm(i+1).ge.0.0_WP.and.i+1.lt.this%imaxo_); i=i+1; end do
      ! Find right j index
      j=max(min(this%jmaxo_-1,j0),this%jmino_)
      do while (pos(2)-this%y (j  ).lt.0.0_WP.and.j  .gt.this%jmino_); j=j-1; end do
      do while (pos(2)-this%y (j+1).ge.0.0_WP.and.j+1.lt.this%jmaxo_); j=j+1; end do
      ! Find right k index
      k=max(min(this%kmaxo_-1,k0),this%kmino_)
      do while (pos(3)-this%zm(k  ).lt.0.0_WP.and.k  .gt.this%kmino_); k=k-1; end do
      do while (pos(3)-this%zm(k+1).ge.0.0_WP.and.k+1.lt.this%kmaxo_); k=k+1; end do
      ! Prepare tri-linear interpolation coefficients
      wx1=(pos(1)-this%xm(i))/(this%xm(i+1)-this%xm(i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%y (j))/(this%y (j+1)-this%y (j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%zm(k))/(this%zm(k+1)-this%zm(k)); wz2=1.0_WP-wz1
      ! Tri-linear interpolation of V
      vel(2)=wz1*(wy1*(wx1*V(i+1,j+1,k+1)  + &
      &                wx2*V(i  ,j+1,k+1)) + &
      &           wy2*(wx1*V(i+1,j  ,k+1)  + &
      &                wx2*V(i  ,j  ,k+1)))+ &
      &      wz2*(wy1*(wx1*V(i+1,j+1,k  )  + &
      &                wx2*V(i  ,j+1,k  )) + &
      &           wy2*(wx1*V(i+1,j  ,k  )  + &
      &                wx2*V(i  ,j  ,k  )))
      ! Interpolate W velocity ------------------------------
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
      do while (pos(3)-this%z (k  ).lt.0.0_WP.and.k  .gt.this%kmino_); k=k-1; end do
      do while (pos(3)-this%z (k+1).ge.0.0_WP.and.k+1.lt.this%kmaxo_); k=k+1; end do
      ! Prepare tri-linear interpolation coefficients
      wx1=(pos(1)-this%xm(i))/(this%xm(i+1)-this%xm(i)); wx2=1.0_WP-wx1
      wy1=(pos(2)-this%ym(j))/(this%ym(j+1)-this%ym(j)); wy2=1.0_WP-wy1
      wz1=(pos(3)-this%z (k))/(this%z (k+1)-this%z (k)); wz2=1.0_WP-wz1
      ! Tri-linear interpolation of W
      vel(3)=wz1*(wy1*(wx1*W(i+1,j+1,k+1)  + &
      &                wx2*W(i  ,j+1,k+1)) + &
      &           wy2*(wx1*W(i+1,j  ,k+1)  + &
      &                wx2*W(i  ,j  ,k+1)))+ &
      &      wz2*(wy1*(wx1*W(i+1,j+1,k  )  + &
      &                wx2*W(i  ,j+1,k  )) + &
      &           wy2*(wx1*W(i+1,j  ,k  )  + &
      &                wx2*W(i  ,j  ,k  )))
      return
   end function get_velocity

   
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
