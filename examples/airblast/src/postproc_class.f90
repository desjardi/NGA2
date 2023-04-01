!> Definition for a postproc class
module postproc_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use timetracker_class, only: timetracker
   implicit none
   private
   
   public :: postproc
   
   !> postproc object
   type :: postproc
      !> Input file for the simulation
      type(inputfile) :: input
      !> Config
      type(config) :: cfg
      !> Time info
      type(timetracker) :: time  
      !> Ensight postprocessing
      type(partmesh) :: pmesh    !< Particle mesh for core output
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      !> Data arrays
      real(WP), dimension(:,:,:), allocatable :: VOF
   contains
      procedure :: analyze
   end type postproc
      
contains
   
   
   !> Analysis of atom simulation
   subroutine analyze(this,flag)
      use parallel, only: amRoot
      implicit none
      class(postproc), intent(inout) :: this
      logical :: flag
      
      ! Switch flag to true
      flag=.true.
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_postproc')
      
      ! Create the config
      create_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='postproc')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      ! Allocate work arrays
      allocate_data: block
         allocate(this%VOF(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_data
      
      
      
      
      
      
      
      ! Create partmesh object for core output
      create_pmesh: block
         this%pmesh=partmesh(nvar=0,nvec=0,name='core')
         ! Transfer stuff here
      end block create_pmesh
      
      ! Perform Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='postproc')
         ! Add variables to output
         call this%ens_out%add_particle('core',this%pmesh)
         ! Output to ensight
         call this%ens_out%write_data(time=0.0_WP)
      end block create_ensight
      
   end subroutine analyze
   

end module postproc_class