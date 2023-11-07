!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use config_class,  only: config
   use pardata_class, only: pardata
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Give ourselves a few test arrays
   real(WP), dimension(:,:,:), allocatable :: Var1,Var2,Var3
   real(WP), dimension(:,:,:), allocatable :: test
   
   !> Config object
   type(config) :: cfg

   !> IO objects
   type(pardata) :: rf,wf
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      
      ! Create config
      create_cfg: block
         use parallel,    only: group
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         real(WP), dimension(:), allocatable :: x
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,nx
         ! Read in grid definition
         call param_read('nx',nx,default=100)
         allocate(x(nx+1))
         ! Create simple rectilinear grid in y and z
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=x,z=x,xper=.true.,yper=.true.,zper=.true.,name='io_tester')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create some test fields
      create_fields: block
         use mathtools, only: twoPi
         integer :: i,j,k
         ! Allocate arrays
         allocate(Var1(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Var2(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Var3(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Initialize to solid body rotation
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  Var1(i,j,k)=1.0_WP
                  Var2(i,j,k)=2.0_WP
                  Var3(i,j,k)=cos(twoPi*cfg%xm(i))*cos(twoPi*cfg%ym(j))*cos(twoPi*cfg%zm(k))
               end do
            end do
         end do
      end block create_fields
      
      
      ! Create pardata object
      create_pardata_write: block
         use param, only: param_read
         integer, dimension(3) :: iopartition
         ! Read in partition
         call param_read('Write partition',iopartition,short='w')
         ! Create pardata object with that IO partition
         call wf%initialize(pg=cfg,iopartition=iopartition,filename='test.file',nval=2,nvar=3)
         ! Assign names to values/variables
         wf%valname=['val1','val2']
         wf%varname=['var1','var2','var3']
         ! Push data to pardata
         call wf%push(name='val1',val=1.0_WP)
         call wf%push(name='val2',val=2.0_WP)
         call wf%push(name='var1',var=Var1)
         call wf%push(name='var2',var=Var2)
         call wf%push(name='var3',var=Var3)
         ! Write out the file
         call wf%write()
      end block create_pardata_write
      
      
      ! Create another pardata object for reading
      create_pardata_read: block
         use param, only: param_read
         integer, dimension(3) :: iopartition
         real(WP) :: val
         ! Read in partition
         call param_read('Read partition',iopartition,short='r')
         ! Create pardata by reading the file
         call rf%initialize(pg=cfg,iopartition=iopartition,fdata='test.file')
         ! Compare our data
         allocate(test(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         call rf%pull(name='val1',val=val ); if (cfg%amRoot) print*,'val1 error',abs(val-1.0_WP)
         call rf%pull(name='val2',val=val ); if (cfg%amRoot) print*,'val2 error',abs(val-2.0_WP)
         call rf%pull(name='var1',var=test); call cfg%integrate(A=abs(test-Var1),integral=val); if (cfg%amRoot) print*,'var1 error',val
         call rf%pull(name='var2',var=test); call cfg%integrate(A=abs(test-Var2),integral=val); if (cfg%amRoot) print*,'var2 error',val
         call rf%pull(name='var3',var=test); call cfg%integrate(A=abs(test-Var3),integral=val); if (cfg%amRoot) print*,'var3 error',val
         ! Try writing again to another name
         call rf%write(fdata='another.file')
      end block create_pardata_read
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
   end subroutine simulation_run
   !> Finalize the NGA2 simulation
   subroutine simulation_final
   end subroutine simulation_final
end module simulation
