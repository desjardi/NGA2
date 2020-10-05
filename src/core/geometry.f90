!> Various definitions and tools for handling NGA grids
!> @todo Provide a flexible multi-grid environment
!> @todo Provide a flexible parallelization strategy
module geometry
   use precision,    only: WP
   use sgrid_class,  only: sgrid,cartesian
   use pgrid_class,  only: pgrid
   !use config_class, only: config
   implicit none
   private
   
   !> Array of grids
   integer :: ngrid
   type(sgrid) :: grid
   type(pgrid), dimension(:), allocatable :: pg
   
   !> Main config
   !type(config) :: cfg
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use string,   only: str_medium
      use param,    only: param_read
      use parallel, only: group,nproc
      implicit none
      integer :: i,ierr,n,grid_group
      integer, dimension(3) :: range
      character(len=str_medium) :: fgrid
      
      ! Create a config from a grid file
      call param_read('Grid file to read',fgrid,short='g')
      grid=sgrid(3,fgrid)
      
      ! Allocate two partitioned grids for testing
      allocate(pg(2))
      
      ! First one based on sgrid
      pg(1)=pgrid(grid,group,[1,1,4])
      
      ! Second one based on grid file
      pg(2)=pgrid(2,fgrid,group,[2,1,2])
      
      ! Try to use HDF5 to create a file
      !call param_read('Grid file to read',fgeom,short='g'); call geometry_write_to_file(fgeom)
      
      ! The logic here would be to either call some mesh creation, or read in a geometry file
      ! Here, we'll start with reading a geometry file
      !call param_read('Grid file to read',fgeom,short='g'); call geometry_read_from_file(fgeom)
      !print*,'fgeom file =',fgeom
      
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
   
   
   
   
   
