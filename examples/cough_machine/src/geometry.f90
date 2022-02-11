!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   public :: geometry_init,t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   
   !> Two configs:
   !> Block 1 (b1) -> zoomed in on mouth, two-phase with VOF
   type(config), target, public :: cfg1
   !> Block 2 (b2) -> zoomed out, Euler-Lagrange
   type(config), target, public :: cfg2
   
   ! Virtual mouth dimensions
   real(WP), parameter :: t_wall =5.0e-3_WP
   real(WP), parameter :: L_mouth=5.0e-2_WP
   real(WP), parameter :: H_mouth=1.0e-2_WP
   real(WP), parameter :: W_mouth=1.0e-2_WP
   real(WP), parameter :: L_film =4.0e-2_WP
   real(WP), parameter :: H_film =1.0e-3_WP
   real(WP), parameter :: W_film =1.0e-2_WP
   real(WP), parameter :: L_lip  =5.0e-3_WP
   
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use param, only: param_read
      implicit none
      
      
      ! First create config for block 1
      create_block1: block
         use sgrid_class, only: cartesian,sgrid
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('1 Lx',Lx); call param_read('1 nx',nx); allocate(x(nx+1))
         call param_read('1 Ly',Ly); call param_read('1 ny',ny); allocate(y(ny+1))
         call param_read('1 Lz',Lz); call param_read('1 nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-L_mouth
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly+0.5_WP*H_mouth
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object with overlap=3 for tpns
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cough_machine_in')
         ! Read in partition
         call param_read('1 Partition',partition,short='p')
         ! Create partitioned grid
         cfg1=config(grp=group,decomp=partition,grid=grid)
         ! Create walls
         cfg1%VF=1.0_WP
         do k=cfg1%kmino_,cfg1%kmaxo_
            do j=cfg1%jmino_,cfg1%jmaxo_
               do i=cfg1%imino_,cfg1%imaxo_
                  ! Zero out wall cells including inside mouth
                  if (cfg1%xm(i).lt.0.0_WP.and.abs(cfg1%zm(k)).lt.0.5_WP*W_mouth+t_wall.and.cfg1%ym(j).gt.-t_wall.and.cfg1%ym(j).lt.H_mouth+t_wall) cfg1%VF(i,j,k)=0.0_WP
                  ! Carve out inside of mouth
                  if (cfg1%xm(i).lt.0.0_WP.and.abs(cfg1%zm(k)).lt.0.5_WP*W_mouth.and.cfg1%ym(j).gt.0.0_WP.and.cfg1%ym(j).lt.H_mouth) cfg1%VF(i,j,k)=1.0_WP
                  ! Carve out tray for liquid film
                  if (cfg1%xm(i).lt.-L_lip.and.cfg1%xm(i).gt.-L_lip-L_film.and.abs(cfg1%zm(k)).lt.0.5_WP*W_film.and.cfg1%ym(j).lt.0.0_WP.and.cfg1%ym(j).gt.-H_film) cfg1%VF(i,j,k)=1.0_WP
               end do
            end do
         end do
      end block create_block1
      
      
      ! Second create config for block 2
      create_block2: block
         use sgrid_class, only: cartesian,sgrid
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('2 Lx',Lx); call param_read('2 nx',nx); allocate(x(nx+1))
         call param_read('2 Ly',Ly); call param_read('2 ny',ny); allocate(y(ny+1))
         call param_read('2 Lz',Lz); call param_read('2 nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-L_mouth
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object with overlap=2 for Euler-Lagrange solver
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cough_machine_out')
         ! Read in partition
         call param_read('2 Partition',partition,short='p')
         ! Create partitioned grid
         cfg2=config(grp=group,decomp=partition,grid=grid)
         ! Create walls
         cfg2%VF=1.0_WP
      end block create_block2
      
      
   end subroutine geometry_init
   
   
end module geometry
