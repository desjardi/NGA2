!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg

   !> Tank radius
   real(WP), public :: r_tank
   real(WP), public :: height_tank

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
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)*Lx/real(nx,WP)-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)*Ly/real(ny,WP)
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)*Lz/real(nz,WP)-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='Stirred Tank')
         
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
         use param,          only: param_read
         integer :: i,j,k
         real(WP) :: xm, ym, zm, tmp1, tmp2
         real(WP) :: Ly

         ! Read in tank radius
         call param_read('tank radius',r_tank)
         call param_read('tank height',height_tank)
         call param_read('Ly',Ly)

         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  xm=cfg%xm(i); ym=cfg%ym(j); zm=cfg%zm(k)
                  tmp1 = sqrt(xm**2 + zm**2) - r_tank ! distance to wall
                  if (ym .gt. 0.5_WP*(Ly+height_tank)) then
                     tmp2 = ym - 0.5_WP*(Ly+height_tank) ! distance to top
                  else if (ym .lt. 0.5_WP*(Ly-height_tank)) then
                     tmp2 = 0.5_WP*(Ly-height_tank) - ym
                  else
                     tmp2 = min(ym-0.5_WP*(Ly+height_tank),0.5_WP*(Ly-height_tank)-ym)
                  end if
                  cfg%Gib(i,j,k)=max(tmp1,tmp2)
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