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
            y(j)=real(j-1,WP)/real(ny,WP)*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='FallingDrop')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         use mpi_f08,  only: MPI_Group,MPI_Group_range_incl
         integer, dimension(3) :: partition
         type(MPI_Group) :: mygrp
         integer :: ierr
         integer, dimension(3,1) :: myrange
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create a group that macthes the partition size
         mygrp=group
         !myrange(:,1)=[0,partition(1)*partition(2)*partition(3)-1,1]
         !call MPI_Group_range_incl(group,1,myrange,mygrp,ierr)
         
         ! Create partitioned grid
         cfg=config(grp=mygrp,decomp=partition,grid=grid)
         call cfg%allprint()
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
         ! Put walls all around
         cfg%VF=0.0_WP
         cfg%VF(cfg%imin_:cfg%imax_,cfg%jmin_:cfg%jmax_,cfg%kmin_:cfg%kmax_)=1.0_WP
         call cfg%sync(cfg%VF)
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
