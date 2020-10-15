!> Various definitions and tools for handling NGA grids
!> @todo Provide a flexible multi-grid environment
!> @todo Provide a flexible parallelization strategy
module geometry
   use config_class,  only: config
   use ensight_class, only: ensight
   implicit none
   private
   
   !> Single config
   type(config) :: cfg,cfg2
   
   !> Ensight output
   type(ensight) :: ens_out
   
   !> Main config
   !type(config) :: cfg
   
   public :: geometry_init
   
contains
   
   function sphere_locator(pg,i1,i2,i3) result(isIn)
      use pgrid_class, only: pgrid
      use precision,   only: WP
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i1,i2,i3
      logical :: isIn
      isIn=.false.
      if (sqrt(pg%xm(i1)**2+pg%ym(i2)**2+pg%zm(i3)**2).lt.0.005_WP) then
         if (pg%ym(i2).gt.0.0_WP) isIn=.true.
      end if
   end function sphere_locator
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use precision,      only: WP
      use string,         only: str_medium
      use param,          only: param_read
      use parallel,       only: group
      use sgrid_class,    only: sgrid,cartesian
      use datafile_class, only: datafile
      implicit none
      ! Test case initialization
      integer, dimension(3) :: partition
      type(sgrid) :: grid
      integer :: i,j,k
      integer :: nx,ny,nz
      real(WP), dimension(:), allocatable :: x,y,z
      real(WP) :: Lx,Ly,Lz,hole_size,hole_dist
      ! We also test the creation of a grid/geom file set
      character(len=str_medium) :: fgrid
      character(len=str_medium) :: fgeom
      
      ! Create a grid from input params
      call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
      call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
      call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
      do i=1,nx+1
         x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
      end do
      do j=1,ny+1
         y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
      end do
      do k=1,nz+1
         z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
      end do
      grid=sgrid(cartesian,1,x,y,z,.true.,.false.,.true.,name='DropletSpreading')
      
      ! Create a config from that grid on our entire group
      call param_read('Partition',partition,short='p')
      cfg=config(group,partition,grid)
      
      ! Create masks for this config
      call param_read('Hole size',hole_size)
      call param_read('Hole dist',hole_dist)
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               if (cfg%ym(j).gt.0.0_WP) then
                  ! Above the plate
                  cfg%mask(i,j,k)=0.0_WP
               else
                  ! This is the plate
                  cfg%mask(i,j,k)=1.0_WP
                  ! Now perforate it
                  if (modulo(cfg%xm(i)-0.5_WP*hole_size,hole_dist).lt.hole_size.and.modulo(cfg%zm(k)-0.5_WP*hole_size,hole_dist).lt.hole_size) cfg%mask(i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      call cfg%maskupdate()
      
      ! Print it back out
      call cfg%write('test')
      
      ! Create a config from grid and config
      call param_read('Grid file',fgrid)
      call param_read('Geom file',fgeom)
      cfg2=config(group,partition,1,fgrid,fgeom)
      
      ! Ensight output
      ens_out=ensight(cfg2,'test')
      
      ! Attempt to create an iterator of a sphere in the center
      test_itr: block
         use itr_class, only: itr
         type(itr) :: sphere_itr
         sphere_itr=itr(cfg,'sphere',sphere_locator)
         call sphere_itr%print()
      end block test_itr
      
      
   end subroutine geometry_init
   
   
end module geometry




   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   ! !> Read in geometry file
   ! subroutine geometry_read_from_file(fgeom)
   !   use hdf5
   !   implicit none
   !   character(len=*) :: fgeom
   !   ! Check that the file exists
   
   ! end subroutine geometry_read_from_file
   
   ! !> Write a geometry file
   ! subroutine geometry_write_to_file(fgeom)
   !   use hdf5
   !   implicit none
   !   character(len=*) :: fgeom                      !< File name
   !   integer(HID_T) :: file_id                      !< File identifier
   !   integer(HID_T) :: group_id                     !< Group identifier
   !   integer(HID_T) :: space_id                     !< Dataspace identifier
   !   integer(HID_T) :: dset_id                      !< Dataset identifier
   
   !   integer(HSIZE_T), dimension(2) :: dims=(/4,6/) !< Dataset dimensions
   !   integer :: rank=2                              !< Dataset rank
   
   !   integer, dimension(4,6) :: data
   
   !   integer :: error                               !< Error flag
   
   !   !character(len=4), parameter :: dsetname="dset" !< Dataset name
   
   !   !integer, dimension(4,6) :: dset_data,data_out  !< Data buffers
   !   !integer(HSIZE_T), dimension(2) :: data_dims
   !   !integer :: i,j
   
   !   ! Initialize FORTRAN interface
   !   call h5open_f(error)
   
   !   ! Create a new file using default properties
   !   call h5fcreate_f(fgeom,H5F_ACC_TRUNC_F,file_id,error)
   
   !   ! Create a grid directory
   !   call h5gcreate_f(file_id,'grid',group_id,error)
   !   ! Open the directory
   !   call h5gopen_f(file_id,'grid',group_id,error)
   
   !   ! Create the data space for the dataset
   !   call h5screate_simple_f(rank,dims,space_id,error)
   !   ! Create the dataset in group "grid" with default properties
   !   call h5dcreate_f(group_id,'mytest',H5T_NATIVE_INTEGER,space_id,dset_id,error)
   !   ! Write the dataset
   !   call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,data,dims,error)
   
   !   ! Close dataspace
   !   call h5sclose_f(space_id,error)
   !   ! Close the dataset
   !   call h5dclose_f(dset_id,error)
   
   !   ! Close the directory
   !   call h5gclose_f(group_id,error)
   
   !   ! Close the file
   !   call h5fclose_f(file_id,error)
   
   !   ! Close FORTRAN interface
   !   call h5close_f(error)
   
   
   
   
   
   !   ! Create the dataspace
   !   !call h5screate_simple_f(rank,dims,dspace_id,error)
   
   !   ! Create the dataset with default properties
   !   !call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
   
   !   ! End access to the dataset and release resources used by it
   !   !call h5dclose_f(dset_id,error)
   
   !   ! Terminate access to the data space
   !   !call h5sclose_f(dspace_id,error)
   
   !   ! Close the file
   !   !call h5fclose_f(file_id,error)
   
   !   ! Close FORTRAN interface
   !   !call h5close_f(error)
   
   
   
   
   
   !   ! Re-initialize FORTRAN interface
   !   !call h5open_f(error)
   
   !   ! Initialize the dset_data array
   !   !do i=1,4
   !   !   do j=1,6
   !   !      dset_data(i,j)=(i-1)*6+j
   !   !   end do
   !   !end do
   
   !   ! Open an existing file
   !   !call h5fopen_f(fgeom,H5F_ACC_RDWR_F,file_id,error)
   
   !   ! Open an existing dataset
   !   !call h5dopen_f(file_id,dsetname,dset_id,error)
   
   !   ! Write the dataset
   !   !data_dims(1)=4
   !   !data_dims(2)=6
   !   !call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,dset_data,data_dims,error)
   
   !   ! Read the dataset
   !   !call h5dread_f(dset_id,H5T_NATIVE_INTEGER,data_out,data_dims,error)
   
   !   ! Close the dataset
   !   !call h5dclose_f(dset_id,error)
   
   !   ! Close the file
   !   !call h5fclose_f(file_id,error)
   
   !   ! Close FORTRAN interface
   !   !call h5close_f(error)
   
   ! end subroutine geometry_write_to_file
   
   
   
   
   
