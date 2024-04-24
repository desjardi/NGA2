!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   !> Mesh info
   integer :: nx
   real(WP) :: Lx

   public :: geometry_init,nx,Lx
   
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
         integer :: i
         real(WP), dimension(:), allocatable :: x
         ! Read in grid definition
         call param_read('Lx',Lx)
         call param_read('nx',nx)
         allocate(x(nx+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=x,z=x,xper=.true.,yper=.true.,zper=.true.,name='HIT')
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
      
      ! Create masks for this config
      create_walls: block
         cfg%VF=1.0_WP
      end block create_walls
      
   end subroutine geometry_init
   
   
end module geometry
