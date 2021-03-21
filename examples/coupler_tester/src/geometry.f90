!> Various definitions and tools for initializing NGA2 config
module geometry
   use mpi_f08,      only: MPI_Group
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Public declarations
   public :: geometry_init
   public :: cfg1,grp1,isInGrp1
   public :: cfg2,grp2,isInGrp2
   
   !> Two groups and partitions, along with logicals
   integer, dimension(3) :: partition1,partition2
   logical :: isInGrp1,isInGrp2
   type(MPI_Group) :: grp1,grp2
   
   !> These are the two configs
   type(config) :: cfg1,cfg2
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use param,    only: param_read
      use parallel, only: comm,group,nproc,rank
      use mpi_f08,  only: MPI_Group,MPI_Group_range_incl
      implicit none
      integer, dimension(3,1) :: grange
      integer :: ierr
      
      ! We start by reading in the two partitions
      call param_read('1 Partition',partition1)
      call param_read('2 Partition',partition2)
      
      ! Create an MPI group along with logical for the first grid on the lowest ranks
      grange(:,1)=[0,product(partition1)-1,1]
      call MPI_Group_range_incl(group,1,grange,grp1,ierr)
      isInGrp1=.false.; if (rank.le.product(partition1)-1) isInGrp1=.true.
      
      ! Create an MPI group along with logical for the second grid on the highest ranks
      grange(:,1)=[nproc-product(partition2),nproc-1,1]
      call MPI_Group_range_incl(group,1,grange,grp2,ierr)
      isInGrp2=.false.; if (rank.ge.nproc-product(partition2)) isInGrp2=.true.
      
      ! Initialize grid 1
      if (isInGrp1) then
         
         create_cfg1: block
            use param,       only: param_read
            use sgrid_class, only: cartesian,sgrid
            type(sgrid) :: grid
            integer :: i,j,k
            integer,  dimension(3) :: ncell
            real(WP), dimension(3) :: dsize
            real(WP), dimension(3) :: dcent
            real(WP), dimension(:), allocatable :: x,y,z
            
            ! Read in grid definition
            call param_read('1 Number of cells',ncell)
            call param_read('1 Domain size'    ,dsize)
            call param_read('1 Domain center'  ,dcent)
            allocate(x(ncell(1)+1))
            allocate(y(ncell(2)+1))
            allocate(z(ncell(3)+1))
            
            ! Create simple rectilinear grid
            do i=1,ncell(1)+1
               x(i)=real(i-1,WP)/real(ncell(1),WP)*dsize(1)-0.5_WP*dsize(1)+dcent(1)
            end do
            do j=1,ncell(2)+1
               y(j)=real(j-1,WP)/real(ncell(2),WP)*dsize(2)-0.5_WP*dsize(2)+dcent(2)
            end do
            do k=1,ncell(3)+1
               z(k)=real(k-1,WP)/real(ncell(3),WP)*dsize(3)-0.5_WP*dsize(3)+dcent(3)
            end do
            
            ! General serial grid object
            grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='grid1')
            
            ! Use it to create a config
            cfg1=config(grp=grp1,decomp=partition1,grid=grid)
            
            ! Do not place walls
            cfg1%VF=1.0_WP
            
         end block create_cfg1
         
      end if
      
      ! Initialize grid 2
      if (isInGrp2) then
         
         create_cfg2: block
            use param,       only: param_read
            use sgrid_class, only: cartesian,sgrid
            type(sgrid) :: grid
            integer :: i,j,k
            integer,  dimension(3) :: ncell
            real(WP), dimension(3) :: dsize
            real(WP), dimension(3) :: dcent
            real(WP), dimension(:), allocatable :: x,y,z
            
            ! Read in grid definition
            call param_read('2 Number of cells',ncell)
            call param_read('2 Domain size'    ,dsize)
            call param_read('2 Domain center'  ,dcent)
            allocate(x(ncell(1)+1))
            allocate(y(ncell(2)+1))
            allocate(z(ncell(3)+1))
            
            ! Create simple rectilinear grid
            do i=1,ncell(1)+1
               x(i)=real(i-1,WP)/real(ncell(1),WP)*dsize(1)-0.5_WP*dsize(1)+dcent(1)
            end do
            do j=1,ncell(2)+1
               y(j)=real(j-1,WP)/real(ncell(2),WP)*dsize(2)-0.5_WP*dsize(2)+dcent(2)
            end do
            do k=1,ncell(3)+1
               z(k)=real(k-1,WP)/real(ncell(3),WP)*dsize(3)-0.5_WP*dsize(3)+dcent(3)
            end do
            
            ! General serial grid object
            grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='grid2')
            
            ! Use it to create a config
            cfg2=config(grp=grp2,decomp=partition2,grid=grid)
            
            ! Do not place walls
            cfg2%VF=1.0_WP
            
         end block create_cfg2
         
      end if
      
   end subroutine geometry_init
   
   
end module geometry
