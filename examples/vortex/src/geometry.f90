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
         integer  :: i,j,k,nx,ny,nz
         real(WP)                            :: Lx,Ly,Lz,rand,sx,sy
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('1 Lx',Lx); call param_read('1 nx',nx); allocate(x(nx+1))
         call param_read('1 Ly',Ly); call param_read('1 ny',ny); allocate(y(ny+1))
         call param_read('1 Lz',Lz); call param_read('1 nz',nz); allocate(z(nz+1))
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
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.false.,name='block1')
         ! Read in partition
         call param_read('1 Partition',partition,short='p')
         ! Create partitioned grid
         cfg1=config(grp=group,decomp=partition,grid=grid)
         ! Create masks for this config
         cfg1%VF=1.0_WP
      end block create_block1

   end subroutine geometry_init

end module geometry
