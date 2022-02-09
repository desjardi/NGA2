!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   public :: geometry_init,t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   
   !> Two config
   type(config), public :: cfg_in,cfg_out
   
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
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid_in,grid_out
      
      
      ! First create inner computational domain
      ! ======================================-
      ! Create inner grid from input params
      create_grid_inner: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         ! Read in grid definition
         call param_read('Inner Lx',Lx); call param_read('Inner nx',nx); allocate(x(nx+1))
         call param_read('Inner Ly',Ly); call param_read('Inner ny',ny); allocate(y(ny+1))
         call param_read('Inner Lz',Lz); call param_read('Inner nz',nz); allocate(z(nz+1))
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
         grid_in=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cough_machine_in')
      end block create_grid_inner
      
      
      ! Create config from grid on our entire group
      create_cfg_inner: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Inner partition',partition,short='p')
         ! Create partitioned grid
         cfg_in=config(grp=group,decomp=partition,grid=grid_in)
      end block create_cfg_inner
      
      
      ! Create masks for config
      create_walls_inner: block
         integer :: i,j,k
         cfg_in%VF=1.0_WP
         do k=cfg_in%kmino_,cfg_in%kmaxo_
            do j=cfg_in%jmino_,cfg_in%jmaxo_
               do i=cfg_in%imino_,cfg_in%imaxo_
                  ! Zero out wall cells including inside mouth
                  if (cfg_in%xm(i).lt.0.0_WP.and.abs(cfg_in%zm(k)).lt.0.5_WP*W_mouth+t_wall.and.cfg_in%ym(j).gt.-t_wall.and.cfg_in%ym(j).lt.H_mouth+t_wall) cfg_in%VF(i,j,k)=0.0_WP
                  ! Carve out inside of mouth
                  if (cfg_in%xm(i).lt.0.0_WP.and.abs(cfg_in%zm(k)).lt.0.5_WP*W_mouth.and.cfg_in%ym(j).gt.0.0_WP.and.cfg_in%ym(j).lt.H_mouth) cfg_in%VF(i,j,k)=1.0_WP
                  ! Carve out tray for liquid film
                  if (cfg_in%xm(i).lt.-L_lip.and.cfg_in%xm(i).gt.-L_lip-L_film.and.abs(cfg_in%zm(k)).lt.0.5_WP*W_film.and.cfg_in%ym(j).lt.0.0_WP.and.cfg_in%ym(j).gt.-H_film) cfg_in%VF(i,j,k)=1.0_WP
               end do
            end do
         end do
      end block create_walls_inner
      
      
      ! Second create outer computational domain
      ! ======================================-
      ! Create outer grid from input params
      create_grid_outer: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         ! Read in grid definition
         call param_read('Outer Lx',Lx); call param_read('Outer nx',nx); allocate(x(nx+1))
         call param_read('Outer Ly',Ly); call param_read('Outer ny',ny); allocate(y(ny+1))
         call param_read('Outer Lz',Lz); call param_read('Outer nz',nz); allocate(z(nz+1))
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
         grid_out=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cough_machine_out')
      end block create_grid_outer
      
      
      ! Create config from grid on our entire group
      create_cfg_outer: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Outer partition',partition,short='p')
         ! Create partitioned grid
         cfg_out=config(grp=group,decomp=partition,grid=grid_out)
      end block create_cfg_outer
      
      
      ! Create masks for config
      create_walls_outer: block
         cfg_out%VF=1.0_WP
      end block create_walls_outer
      
      
   end subroutine geometry_init
   
   
end module geometry
