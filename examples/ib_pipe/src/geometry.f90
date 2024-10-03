!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg
   
   !> Pipe diameter
   real(WP), public :: D
   
   public :: geometry_init
   
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Pipe length',Lx)
         call param_read('Pipe diameter',D)
         call param_read('ny',ny); allocate(y(ny+1))
         call param_read('nx',nx); allocate(x(nx+1))
         call param_read('nz',nz); allocate(z(nz+1))
         
         dx=Lx/real(nx,WP)
         no=4
         if (ny.gt.1) then
            Ly=D+real(2*no,WP)*D/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=D+real(2*no,WP)*D/real(ny-2*no,WP)
         else
            Lz=dx
         end if
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='pipe')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
         use ibconfig_class, only: bigot,sharp
         integer :: i,j,k
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  cfg%Gib(i,j,k)=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)-0.5_WP*D
               end do
            end do
         end do
          ! Get normal vector
         call cfg%calculate_normal()
         ! Get VF field
         call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry