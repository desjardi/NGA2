!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private

   public :: geometry_init

   !> Config:
   type(config), target, public :: cfg1


contains

   !> Initialization of problem geometry
   subroutine geometry_init
      use param, only: param_read
      implicit none

      ! First create config for block 1
      create_block1: block
         use sgrid_class, only: cartesian,sgrid
         use parallel,    only: group
         use random,      only: random_uniform,random_initialize
         use mathtools,   only: twoPi
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,rand,sx,sy
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('1 Lx',Lx); call param_read('1 nx',nx); allocate(x(nx+1))
         call param_read('1 Ly',Ly); call param_read('1 ny',ny); allocate(y(ny+1))
         call param_read('1 Lz',Lz); call param_read('1 nz',nz); allocate(z(nz+1))
         ! Read in streching coefficient
         call param_read('Stretching x',sx)
         call param_read('Stretching y',sy)
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
         ! General serial grid object with overlap=2 for Euler-Lagrange solver
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='block1')
         ! Read in partition
         call param_read('1 Partition',partition,short='p')
         ! Create partitioned grid
         cfg1=config(grp=group,decomp=partition,grid=grid)
         ! Apply stretching
         !if (sx.gt.0.0_WP) then
         !   x = x + sx*sin(1.0_WP*twoPi*x/Lx)
         !else if (sx.lt.0.0_WP) then
         !   call random_initialize
         !   do i=2,nx+1
         !      call random_number(rand)
         !      rand=1.0_WP+sx-2.0_WP*sx*rand
         !      x(i) = x(i-1)+rand*(x(i)-x(i-1))
         !   end do
         !   x=x+0.5_WP*Lx
         !   rand=x(nx+1)
         !   x=x*Lx/rand
         !   x=x-0.5_WP*Lx
         !end if
         !if (sy.gt.0.0_WP) then
         !   y = y + sy*sin(1.0_WP*twoPi*y/Ly)
         !else if (sy.lt.0.0_WP) then
         !   call random_initialize
         !   do j=2,ny+1
         !      call random_number(rand)
         !      rand=1.0_WP+sy-2.0_WP*sy*rand
         !      y(j) = y(j-1)+rand*(y(j)-y(j-1))
         !   end do
         !   y=y+0.5_WP*Ly
         !   rand=y(ny+1)
         !   y=y*Ly/rand
         !   y=y-0.5_WP*Ly
         !end if
      end block create_block1

   end subroutine geometry_init

end module geometry
