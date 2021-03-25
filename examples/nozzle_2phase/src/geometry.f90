!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   public :: geometry_init
   
   !> Single config files
   type(config), public :: cfg !< This is the atomizing flow config
   
   ! Some nozzle geometry data
   real(WP), parameter :: nozzle_end   =-0.09130000_WP     !< back end of the nozzle plenum
   real(WP), parameter :: xinj_dist    =-0.06623329454_WP  !< x location of the gas port axis
   real(WP), parameter :: inj_norm_diam=+0.016002_WP       !< diameter of the gas injection ports
   real(WP) :: rli0  !< Liquid pipe wall inner diameter at xmin, if xmin>nozzle_end
   real(WP) :: rlo0  !< Liquid pipe wall outer diameter at xmin, if xmin>nozzle_end
   real(WP) :: rgi0  !< Gas annular wall pipe inner diameter at xmin, if xmin>nozzle_end
   real(WP) :: rgo0  !< Gas annular wall pipe outer diameter at xmin, if xmin>nozzle_end
   public :: xinj_dist,inj_norm_diam,rli0,rlo0,rgi0,rgo0 ! These will be needed by the bconds later
   
contains
   
   ! -------------------------------------------------- !
   ! Bisection routine                                  !
   ! gets an array and its size as well as a position x !
   ! Returns the index between 1 and nx                 !
   ! x(iloc)<=xloc<x(iloc+1)                            !
   ! Assuming xarray is monotonically increasing        !
   ! -------------------------------------------------- !
   subroutine bisection(xloc,iloc,x,nx)
      implicit none
      real(WP), intent(in ) :: xloc
      integer,  intent(out) :: iloc
      integer,  intent(in ) :: nx
      real(WP), dimension(1:nx), intent(in) :: x
      integer :: il,im,iu
      ! Take care of outside points
      if (xloc.lt.x(1)) then
         iloc=1
      else if (xloc.ge.x(nx)) then
         iloc=nx-1
      else
         ! Initialize lower and upper limits
         il=1
         iu=nx
         ! While not done
         do while (iu-il.gt.1)
            ! Compute a mid-point
            im=(iu+il)/2
            ! Replace lower of upper limit as appropriate
            if (xloc.ge.x(im)) then
               il=im
            else
               iu=im
            end if
         end do
         ! Return
         iloc = il
      end if
   end subroutine bisection
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      real(WP), dimension(   2) :: x_swh1,swh1    ! dli
      real(WP), dimension( 889) :: x_swh2,swh2    ! dlo
      real(WP), dimension(3302) :: x_swh3,swh3    ! dgi
      real(WP), dimension(5053) :: x_swh4,swh4    ! dgo
      
      
      ! Read in the nozzle CAD data first
      read_cad: block
         use mpi_f08,  only: MPI_BCAST
         use parallel, only: comm,amRoot,MPI_REAL_WP
         integer :: ierr,iunit,i
         ! Only the global root process reads the CAD files
         if (amRoot) then
            open(newunit=iunit,file='cad/wallrli.txt',form='formatted',status='old',access='stream',iostat=ierr); do i=1,   2; read(iunit,'(2F12.8)') x_swh1(i),swh1(i); end do; close(iunit)
            open(newunit=iunit,file='cad/wallrlo.txt',form='formatted',status='old',access='stream',iostat=ierr); do i=1, 889; read(iunit,'(2F12.8)') x_swh2(i),swh2(i); end do; close(iunit)
            open(newunit=iunit,file='cad/wallrgi.txt',form='formatted',status='old',access='stream',iostat=ierr); do i=1,3302; read(iunit,'(2F12.8)') x_swh3(i),swh3(i); end do; close(iunit)
            open(newunit=iunit,file='cad/wallrgo.txt',form='formatted',status='old',access='stream',iostat=ierr); do i=1,5053; read(iunit,'(2F12.8)') x_swh4(i),swh4(i); end do; close(iunit)
         end if
         ! Then the root broadcasts
         call MPI_BCAST(x_swh1,   2,MPI_REAL_WP,0,comm,ierr); call MPI_BCAST(swh1,   2,MPI_REAL_WP,0,comm,ierr)
         call MPI_BCAST(x_swh2, 889,MPI_REAL_WP,0,comm,ierr); call MPI_BCAST(swh2, 889,MPI_REAL_WP,0,comm,ierr)
         call MPI_BCAST(x_swh3,3302,MPI_REAL_WP,0,comm,ierr); call MPI_BCAST(swh3,3302,MPI_REAL_WP,0,comm,ierr)
         call MPI_BCAST(x_swh4,5053,MPI_REAL_WP,0,comm,ierr); call MPI_BCAST(swh4,5053,MPI_REAL_WP,0,comm,ierr)
      end block read_cad
      
      
      ! Create second grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xmax
         real(WP), dimension(:), allocatable :: x,y,z
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         call param_read('xmax',xmax)
         x=x-x(nx+1)+xmax
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object with overlap=3 for tpns
         if (nz.eq.1) then
            grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='nozzle_exterior')
         else
            grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='nozzle_exterior')
         end if
      end block create_grid
      
      
      ! Create config2 from grid2 on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create masks for config2
      create_nozzle_exit: block
         use messager, only: warn
         integer :: i,j,k
         real(WP) :: xloc,rad,run,rise,hyp
         real(WP), dimension(:,:), allocatable :: swh
         
         ! Basic check
         if (cfg%xm(cfg%imin  ).lt.xinj_dist+0.5_WP*inj_norm_diam) call warn('[geometry init] cfg extends beyond the nozzle ports!')
         if (cfg%x (cfg%imax+1).lt.0.0_WP) call warn('[geometry init] cfg does not extend beyond the nozzle tip!')
         
         ! Interpolate CAD on our mesh
         allocate(swh(4,cfg%imino_:cfg%imaxo_)); swh=0.0_WP
         do i=cfg%imino_,cfg%imaxo_
            xloc=cfg%xm(i)
            ! Skip if past nozzle exit on the right
            if (xloc.gt.0.0_WP) cycle
            ! Skip if past nozzle backend on the left
            if (cfg%xm(i).lt.nozzle_end) swh(:,i)=0.0_WP
            ! Interpolate switch height
            call bisection(xloc,k,x_swh1,   2)
            swh(1,i)=(xloc-x_swh1(k))/(x_swh1(k+1)-x_swh1(k))*swh1(k+1)+(x_swh1(k+1)-xloc)/(x_swh1(k+1)-x_swh1(k))*swh1(k)
            call bisection(xloc,k,x_swh2, 890)
            swh(2,i)=(xloc-x_swh2(k))/(x_swh2(k+1)-x_swh2(k))*swh2(k+1)+(x_swh2(k+1)-xloc)/(x_swh2(k+1)-x_swh2(k))*swh2(k)
            call bisection(xloc,k,x_swh3,3302)
            swh(3,i)=(xloc-x_swh3(k))/(x_swh3(k+1)-x_swh3(k))*swh3(k+1)+(x_swh3(k+1)-xloc)/(x_swh3(k+1)-x_swh3(k))*swh3(k)
            call bisection(xloc,k,x_swh4,5053)
            swh(4,i)=(xloc-x_swh4(k))/(x_swh4(k+1)-x_swh4(k))*swh4(k+1)+(x_swh4(k+1)-xloc)/(x_swh4(k+1)-x_swh4(k))*swh4(k)
         end do
         
         ! If the domain does not extend past the back on the left, store the diameters for use by the boundary conditions
         if (cfg%xm(cfg%imin).gt.nozzle_end) then
            xloc=cfg%xm(cfg%imin)
            call bisection(xloc,k,x_swh1,   2)
            rli0=(xloc-x_swh1(k))/(x_swh1(k+1)-x_swh1(k))*swh1(k+1)+(x_swh1(k+1)-xloc)/(x_swh1(k+1)-x_swh1(k))*swh1(k)
            call bisection(xloc,k,x_swh2, 890)
            rlo0=(xloc-x_swh2(k))/(x_swh2(k+1)-x_swh2(k))*swh2(k+1)+(x_swh2(k+1)-xloc)/(x_swh2(k+1)-x_swh2(k))*swh2(k)
            call bisection(xloc,k,x_swh3,3302)
            rgi0=(xloc-x_swh3(k))/(x_swh3(k+1)-x_swh3(k))*swh3(k+1)+(x_swh3(k+1)-xloc)/(x_swh3(k+1)-x_swh3(k))*swh3(k)
            call bisection(xloc,k,x_swh4,5053)
            rgo0=(xloc-x_swh4(k))/(x_swh4(k+1)-x_swh4(k))*swh4(k+1)+(x_swh4(k+1)-xloc)/(x_swh4(k+1)-x_swh4(k))*swh4(k)
         else
            rli0=0.0_WP
            rlo0=0.0_WP
            rgi0=0.0_WP
            rgo0=0.0_WP
         end if
         
         ! Write VF array - includes the overlap here
         cfg%VF=1.0_WP
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  rad=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
                  ! Skip if past nozzle exit
                  if (cfg%xm(i).gt.0.0_WP) cycle
                  ! Skip if past nozzle end
                  if (cfg%xm(i).lt.nozzle_end) then
                     cfg%VF(i,j,k)=0.0_WP
                     cycle
                  end if
                  ! first nozzle wall
                  if (swh(1,i).le.rad.and.rad.lt.swh(2,i)) cfg%VF(i,j,k)=0.0_WP
                  ! second nozzle wall
                  if (swh(3,i).le.rad.and.rad.lt.swh(4,i)) cfg%VF(i,j,k)=0.0_WP
                  ! Carve out injector ports
                  if (swh(3,i).le.rad.and.rad.lt.swh(4,i)) then
                     ! Check injector x-y plane
                     hyp =norm2([cfg%xm(i),cfg%ym(j),cfg%zm(k)]-[xinj_dist,0.0_WP,0.0_WP])
                     rise=abs(cfg%zm(k))
                     run =sqrt(hyp**2-rise**2)
                     if (run.le.inj_norm_diam/2.0_WP) cfg%VF(i,j,k)=1.0_WP
                     ! Check injector x-z plane
                     hyp =norm2([cfg%xm(i),cfg%ym(j),cfg%zm(k)]-[xinj_dist,0.0_WP,0.0_WP])
                     rise=abs(cfg%ym(j))
                     run =sqrt(hyp**2-rise**2)
                     if (run.le.inj_norm_diam/2.0_WP) cfg%VF(i,j,k)=1.0_WP
                  end if
               end do
            end do
         end do
         
         ! Deallocate CAD on our mesh
         deallocate(swh)
         
      end block create_nozzle_exit
      
      
   end subroutine geometry_init
   
   
end module geometry
