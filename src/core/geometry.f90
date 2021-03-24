!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg,cfg2
   
   public :: geometry_init
   
contains
   
   
   !> Function that assists the creation of a sphere iterator
   function sphere_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (sqrt(pg%xm(i)**2+pg%ym(j)**2+pg%zm(k)**2).lt.0.005_WP) then
         if (pg%ym(j).gt.0.0_WP) isIn=.true.
      end if
   end function sphere_locator
   
   
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
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='DropletSpreading')
         
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
         real(WP) :: hole_size,hole_dist,hole_depth
         integer :: i,j,k
         
         ! Read in wall definitions
         call param_read('Hole size',hole_size)
         call param_read('Hole dist',hole_dist)
         call param_read('Hole depth',hole_depth)
         
         ! Write VF array - includes the overlap here
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  if (cfg%ym(j).gt.0.0_WP) then
                     ! Above the plate
                     cfg%VF(i,j,k)=1.0_WP
                  else if (cfg%ym(j).lt.-hole_depth) then
                     ! Below the plate
                     cfg%VF(i,j,k)=1.0_WP
                  else
                     ! This is the plate
                     cfg%VF(i,j,k)=0.0_WP
                     ! Now perforate it
                     if ((0.5_WP*hole_dist-abs(modulo(cfg%xm(i),hole_dist)-0.5_WP*hole_dist)).lt.0.5_WP*hole_size.and.&
                     &   (0.5_WP*hole_dist-abs(modulo(cfg%zm(k),hole_dist)-0.5_WP*hole_dist)).lt.0.5_WP*hole_size) cfg%VF(i,j,k)=1.0_WP
                  end if
               end do
            end do
         end do
         
      end block create_walls
      
      
      ! Test output of a config to the disk
      call cfg%write('test')
      
      
      ! Create a config from grid and config
      create_cfg2: block
         use parallel, only: group
         use string,   only: str_medium
         character(len=str_medium) :: fgrid
         character(len=str_medium) :: fgeom
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Read in the grid/geom file names
         call param_read('Grid file',fgrid)
         call param_read('Geom file',fgeom)
         
         ! Create the cfg by reading in the file
         cfg2=config(grp=group,decomp=partition,no=2,fgrid=fgrid,fgeom=fgeom)
         
      end block create_cfg2
      
      
      ! Attempt to create an iterator of a sphere in the center
      test_itr: block
         use iterator_class, only: iterator
         type(iterator) :: sphere_itr
         
         ! Create an unstructured iterator and print it out
         sphere_itr=iterator(cfg2,'sphere',sphere_locator)
         call sphere_itr%print()
         
      end block test_itr
      
      
   end subroutine geometry_init
   
   
end module geometry
