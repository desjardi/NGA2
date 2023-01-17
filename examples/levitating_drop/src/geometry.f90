!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read,param_exists
      use parallel,    only: amRoot
      use mathtools,   only: pi
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,ny_uni
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         real(WP) :: ry, y_uni, dy, dy_tmp, y_tmp
         
         ! Read in grid definition: can be 2D or 3D
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny)
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid in x and z
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         ! Create grid in y, with stretching if specified
         dy = Ly/ny
         if (param_exists('Stretching ratio in y')) then
            call param_read('Stretching ratio in y',ry)
            call param_read('Stretching begins at y',y_uni)

            ! Determine actual number of points to be used
            ny = 0
            y_tmp = 0.0_WP
            dy_tmp = dy
            do while (y_tmp.lt.0.5_WP*Ly)
               ny = ny+2 ! Symmetric
               if (y_tmp.gt.y_uni) dy_tmp = dy_tmp*ry
               y_tmp = y_tmp + dy_tmp
            end do
            ny_uni = ceiling(y_uni/dy)
            if (amRoot) print*,'ny with stretched mesh',ny
            allocate(y(ny+1))
            y(ny/2+1) = 0.0_WP
            do j=ny/2+2,ny/2+ny_uni+1
               y(j) = y(j-1)+dy
            end do
            do j=ny/2+ny_uni+2,ny+1
               dy = dy*ry
               y(j) = y(j-1)+dy
            end do
            do j=1,ny/2
               y(j) = -y(ny+2-j)
            end do
         else
            ! Uniform mesh
            allocate(y(ny+1))
            do j=1,ny+1
               y(j)=real(j-1,WP)*dy-0.5_WP*Ly
            end do
         end if
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='LevDrop')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
   end subroutine geometry_init
   
   
end module geometry
