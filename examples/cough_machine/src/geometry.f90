!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   public :: geometry_init,t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   
   !> Three configs:
   !> Block 1 (b1) -> turbulent duct, generate inflow velocity into b2
   type(config), target, public :: cfg1
   !> Block 2 (b2) -> zoomed in on mouth, two-phase with VOF
   type(config), target, public :: cfg2
   !> Block 3 (b3) -> zoomed out, Euler-Lagrange
   type(config), target, public :: cfg3
   
   ! Virtual mouth dimensions
   real(WP), parameter :: t_wall =5.0e-3_WP
   ! real(WP), parameter :: L_mouth=5.0e-2_WP
   real(WP), parameter :: H_mouth=1.0e-2_WP
   real(WP), parameter :: W_mouth=1.0e-2_WP
   ! real(WP), parameter :: L_film =4.0e-2_WP
   ! real(WP), parameter :: H_film =1.0e-3_WP
   ! real(WP), parameter :: W_film =1.0e-2_WP
   ! real(WP), parameter :: L_lip  =5.0e-3_WP 
   
   ! Dimensions for extended mouth
   real(WP), parameter :: H_film =2.8e-3_WP
   real(WP), parameter :: L_lip  =1.0e-1_WP
   real(WP), parameter :: L_film =2.8_WP
   real(WP), parameter :: t_base =1.0e-1_WP
   real(WP), parameter :: W_film =3.0e-2_WP
   real(WP), parameter :: L_mouth=1.0e-1_WP
   
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use param, only: param_read
      implicit none
      
      
      ! Create config for block 1
      create_block1: block
         use sgrid_class, only: cartesian,sgrid
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,stretch
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('1 Lx',Lx); call param_read('1 nx',nx); allocate(x(nx+1))
         call param_read('1 Ly',Ly); call param_read('1 ny',ny); allocate(y(ny+1)); call param_read('1 Stretching',stretch)
         call param_read('1 Lz',Lz); call param_read('1 nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid in x and z, tanh-stretched grid in y
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            ! y(j)=0.5_WP*Ly*tanh(stretch*(2.0_WP*real(j-1,WP)/real(ny,WP)-1.0_WP))/tanh(stretch)
            y(j)=0.5_WP*Ly*tanh(stretch*(2.0_WP*real(j-1,WP)/real(ny,WP)-1.0_WP))/tanh(stretch)+(Ly/2.0_WP)
         end do
         do k=1,nz+1
            z(k)=0.5_WP*Lz*tanh(stretch*(2.0_WP*real(k-1,WP)/real(nz,WP)-1.0_WP))/tanh(stretch)
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.false.,name='duct')
         ! Read in partition
         call param_read('1 Partition',partition,short='p')
         ! Create partitioned grid
         cfg1=config(grp=group,decomp=partition,grid=grid)
         ! Create masks for this config
         cfg1%VF=0.0_WP
         cfg1%VF(cfg1%imin_:cfg1%imax_,cfg1%jmin_:cfg1%jmax_,cfg1%kmin_:cfg1%kmax_)=1.0_WP
         call cfg1%sync(cfg1%VF)
      end block create_block1
      
      ! Create config for block 2
      create_block2: block
         use sgrid_class, only: cartesian,sgrid
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz,unx,uny,unz
         real(WP) :: Lx,Ly,Lz,Sx,Sy,Sz,stretch
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('2 Uniform Lx',Lx); call param_read('2 Uniform nx',unx); call param_read('2 Total nx',nx); allocate(x(nx+1))
         call param_read('2 Uniform Ly',Ly); call param_read('2 Uniform ny',uny); call param_read('2 Total ny',ny); allocate(y(ny+1))
         call param_read('2 Uniform Lz',Lz); call param_read('2 Uniform nz',unz); call param_read('2 Total nz',nz); allocate(z(nz+1))
         call param_read('2 Stretching x',Sx)
         call param_read('2 Stretching y',Sy)
         call param_read('2 Stretching z',Sz)
         ! Create simple rectilinear grid
         do i=1,unx+1
            x(i)=real(i-1,WP)/real(unx,WP)*Lx
         end do
         do j=1+(ny-uny)/2,1+(ny+uny)/2
            y(j)=real(j-(1+(ny-uny)/2),WP)/real(uny,WP)*(Ly+t_base)-t_base
         end do
         do k=1+(nz-unz)/2,1+(nz+unz)/2
            z(k)=real(k-(1+(nz-unz)/2),WP)/real(unz,WP)*Lz-0.5_WP*Lz
         end do
         ! Grid for square duct domain
         ! do i=1,unx+1
         !    x(i)=real(i-1,WP)/real(unx,WP)*Lx-L_mouth
         ! end do
         ! do j=1+(ny-uny)/2,1+(ny+uny)/2
         !    y(j)=real(j-(1+(ny-uny)/2),WP)/real(uny,WP)*Ly-0.5_WP*Ly+0.5_WP*H_mouth
         ! end do
         ! do k=1+(nz-unz)/2,1+(nz+unz)/2
         !    z(k)=real(k-(1+(nz-unz)/2),WP)/real(unz,WP)*Lz-0.5_WP*Lz
         ! end do
         ! ! Add stretching
         ! do i=unx+2,nx+1
         !    x(i)=(1.0_WP+Sx)*x(i-1)-Sx*x(i-2)
         ! end do
         ! do j=(ny+uny)/2+2,ny+1
         !    y(j)=(1.0_WP+Sy)*y(j-1)-Sy*y(j-2)
         ! end do
         ! do j=(ny-uny)/2,1,-1
         !    y(j)=(1.0_WP+Sy)*y(j+1)-Sy*y(j+2)
         ! end do
         ! do k=(nz+unz)/2+2,nz+1
         !    z(k)=(1.0_WP+Sz)*z(k-1)-Sz*z(k-2)
         ! end do
         ! do k=(nz-unz)/2,1,-1
         !    z(k)=(1.0_WP+Sz)*z(k+1)-Sz*z(k+2)
         ! end do
         ! General serial grid object with overlap=3 for tpns
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cough_machine_in')
         ! Read in partition
         call param_read('2 Partition',partition,short='p')
         ! Create partitioned grid
         cfg2=config(grp=group,decomp=partition,grid=grid)
         ! Create walls - extended mouth
         cfg2%VF=1.0_WP
         do k=cfg2%kmino_,cfg2%kmaxo_
            do j=cfg2%jmino_,cfg2%jmaxo_
               do i=cfg2%imino_,cfg2%imaxo_
                  ! Zero out base cells
                  if (cfg2%ym(j).lt.0.0_WP) cfg2%VF(i,j,k)=0.0_WP
                  ! Carve out tray for liquid film
                  if (cfg2%xm(i).gt.L_lip.and.cfg2%xm(i).lt.L_lip+L_film.and.cfg2%ym(j).lt.0.0_WP.and.cfg2%ym(j).gt.-H_film) cfg2%VF(i,j,k)=1.0_WP
               end do
            end do
         end do
         ! ! Create walls - duct geometry
         ! cfg2%VF=1.0_WP
         ! do k=cfg2%kmino_,cfg2%kmaxo_
         !    do j=cfg2%jmino_,cfg2%jmaxo_
         !       do i=cfg2%imino_,cfg2%imaxo_
         !          ! Zero out wall cells including inside mouth
         !          if (cfg2%xm(i).lt.0.0_WP.and.abs(cfg2%zm(k)).lt.0.5_WP*W_mouth+t_wall.and.cfg2%ym(j).gt.-t_wall.and.cfg2%ym(j).lt.H_mouth+t_wall) cfg2%VF(i,j,k)=0.0_WP
         !          ! Carve out inside of mouth
         !          if (cfg2%xm(i).lt.0.0_WP.and.abs(cfg2%zm(k)).lt.0.5_WP*W_mouth.and.cfg2%ym(j).gt.0.0_WP.and.cfg2%ym(j).lt.H_mouth) cfg2%VF(i,j,k)=1.0_WP
         !          ! Carve out tray for liquid film
         !          if (cfg2%xm(i).lt.-L_lip.and.cfg2%xm(i).gt.-L_lip-L_film.and.abs(cfg2%zm(k)).lt.0.5_WP*W_film.and.cfg2%ym(j).lt.0.0_WP.and.cfg2%ym(j).gt.-H_film) cfg2%VF(i,j,k)=1.0_WP
         !       end do
         !    end do
         ! end do
      end block create_block2
      
      
      ! Create config for block 3
      create_block3: block
         use sgrid_class, only: cartesian,sgrid
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz,unx,uny,unz
         real(WP) :: Lx,Ly,Lz,Sx,Sy,Sz
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('3 Uniform Lx',Lx); call param_read('3 Uniform nx',unx); call param_read('3 Total nx',nx); allocate(x(nx+1))
         call param_read('3 Uniform Ly',Ly); call param_read('3 Uniform ny',uny); call param_read('3 Total ny',ny); allocate(y(ny+1))
         call param_read('3 Uniform Lz',Lz); call param_read('3 Uniform nz',unz); call param_read('3 Total nz',nz); allocate(z(nz+1))
         call param_read('3 Stretching x',Sx)
         call param_read('3 Stretching y',Sy)
         call param_read('3 Stretching z',Sz)
         ! Create simple rectilinear grid in the uniform region
         do i=1,unx+1
            x(i)=real(i-1,WP)/real(unx,WP)*Lx-L_mouth
         end do
         do j=1+(ny-uny)/2,1+(ny+uny)/2
            y(j)=real(j-(1+(ny-uny)/2),WP)/real(uny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1+(nz-unz)/2,1+(nz+unz)/2
            z(k)=real(k-(1+(nz-unz)/2),WP)/real(unz,WP)*Lz-0.5_WP*Lz
         end do
         ! Add stretching
         do i=unx+2,nx+1
            x(i)=(1.0_WP+Sx)*x(i-1)-Sx*x(i-2)
         end do
         do j=(ny+uny)/2+2,ny+1
            y(j)=(1.0_WP+Sy)*y(j-1)-Sy*y(j-2)
         end do
         do j=(ny-uny)/2,1,-1
            y(j)=(1.0_WP+Sy)*y(j+1)-Sy*y(j+2)
         end do
         do k=(nz+unz)/2+2,nz+1
            z(k)=(1.0_WP+Sz)*z(k-1)-Sz*z(k-2)
         end do
         do k=(nz-unz)/2,1,-1
            z(k)=(1.0_WP+Sz)*z(k+1)-Sz*z(k+2)
         end do
         ! General serial grid object with overlap=2 for Euler-Lagrange solver
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cough_machine_out')
         ! Read in partition
         call param_read('2 Partition',partition,short='p')
         ! Create partitioned grid
         cfg3=config(grp=group,decomp=partition,grid=grid)
         ! Create walls
         cfg3%VF=1.0_WP
      end block create_block3
      
      
   end subroutine geometry_init
   
   
end module geometry
