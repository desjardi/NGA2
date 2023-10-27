!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use config_class,  only: config
   use cclabel_class, only: cclabel
   use ensight_class, only: ensight
   implicit none
   private
   
   !> Config and ccl
   type(config)  :: cfg
   type(cclabel) :: ccl
   
   !> A level set array
   real(WP), dimension(:,:,:), allocatable :: G
   
   !> Global array of particles
   integer :: np
   real(WP), dimension(:,:), allocatable :: xp
   real(WP), dimension(:),   allocatable :: dp

   !> Ensight output to visualize results
   type(ensight) :: ens
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      

      ! Create the mesh
      create_mesh: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,nx
         ! Read in grid definition
         call param_read('nx',nx); allocate(x(nx+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=x,z=x,xper=.true.,yper=.true.,zper=.true.,name='ccl_test')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_mesh
      

      ! Create the particles
      create_particles: block
         use param,    only: param_read
         use random,   only: random_uniform
         use mpi_f08,  only: MPI_BCAST
         use parallel, only: MPI_REAL_WP
         integer :: n,ierr
         real(WP) :: dmin,dmax
         ! Read in number of particles
         call param_read('np',np)
         ! Read in diameter range
         call param_read('dmin',dmin)
         call param_read('dmax',dmax)
         ! Allocate position and diameter arrays
         allocate(xp(3,np),dp(np))
         ! Root process randomly initializes values
         if (cfg%amRoot) then
            do n=1,np
               xp(:,n)=[random_uniform(cfg%x(cfg%imin),cfg%x(cfg%imax+1)),&
               &        random_uniform(cfg%y(cfg%jmin),cfg%y(cfg%jmax+1)),&
               &        random_uniform(cfg%z(cfg%kmin),cfg%z(cfg%kmax+1))]
               dp(n)=random_uniform(dmin,dmax)
            end do
         end if
         ! Broadcast from root
         call MPI_BCAST(xp,3*np,MPI_REAL_WP,0,cfg%comm,ierr)
         call MPI_BCAST(dp,1*np,MPI_REAL_WP,0,cfg%comm,ierr)
      end block create_particles
      

      ! Create distance field
      create_distance: block
         integer :: n,i,j,k
         real(WP), dimension(3) :: pos
         ! Allocate array
         allocate(G(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); G=huge(1.0_WP)
         ! Measure distance from each particle
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Loop over each particle
                  do n=1,np
                     pos(1)=cfg%xm(i)-xp(1,n); if (pos(1).gt.0.5_WP*cfg%xL) pos(1)=pos(1)-cfg%xL; if (pos(1).lt.-0.5_WP*cfg%xL) pos(1)=pos(1)+cfg%xL
                     pos(2)=cfg%ym(j)-xp(2,n); if (pos(2).gt.0.5_WP*cfg%yL) pos(2)=pos(2)-cfg%yL; if (pos(2).lt.-0.5_WP*cfg%yL) pos(2)=pos(2)+cfg%yL
                     pos(3)=cfg%zm(k)-xp(3,n); if (pos(3).gt.0.5_WP*cfg%zL) pos(3)=pos(3)-cfg%zL; if (pos(3).lt.-0.5_WP*cfg%zL) pos(3)=pos(3)+cfg%zL
                     G(i,j,k)=min(G(i,j,k),norm2(pos)-0.5_WP*dp(n))
                  end do
               end do
            end do
         end do
      end block create_distance
      

      ! Create CCL
      create_ccl: block
         integer :: i,j,k
         ! Initialize CCL
         call ccl%initialize(cfg=cfg,name='ccl_test')
         ! Tag the cells
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  if (G(i,j,k).lt.0.0_WP) ccl%tagged(i,j,k)=.true.
               end do
            end do
         end do
         ! Perform CCL
         call ccl%build()
      end block create_ccl
      
      ! Ensight output
      ensight_output: block
         ens=ensight(cfg=cfg,name='ccl_test')
         call ens%add_scalar('G',G)
         call ens%add_scalar('id',ccl%id)
         call ens%write_data(0.0_WP)
      end block ensight_output
      
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
   end subroutine simulation_run
   !> Finalize the NGA2 simulation
   subroutine simulation_final
   end subroutine simulation_final
end module simulation