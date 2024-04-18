!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,      only: WP
   use film_class,     only: film
   use config_class,   only: config
   use vfs_class,      only: vfs
   use surfmesh_class, only: surfmesh
   use ensight_class,  only: ensight
   use event_class,    only: event
   
   implicit none
   private
   
   !> Film simulation
   type(film) :: dns
   
   !> R2P post-processing
   type(config)   :: cfg
   type(vfs)      :: vf
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: pproc_evt
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   !> Initialization of our simulation
   subroutine simulation_init
      
      ! Initialize film simulation
      call dns%init()
      
      ! Create a coarse grid config
      coarse_grid: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz,coarse_factor
         real(WP) :: Lx,Ly,Lz
         ! Coarsening factor
         call param_read('Coarsening factor',coarse_factor)
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); nx=nx/coarse_factor; allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); ny=ny/coarse_factor; allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); nz=nz/coarse_factor; allocate(z(nz+1))
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
         ! General coarse serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='coarseFilm')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block coarse_grid
      
      ! Create a coarse VF solver
      coarse_vof: block
         use vfs_class, only: r2p,remap
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=r2p,transport_method=remap,name='coarseVOF')
      end block coarse_vof

      ! Create a coarse surfmesh object for interface polygon output
      coarse_smesh: block
         smesh=surfmesh(nvar=5,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         smesh%varname(3)='edge_sensor'
         smesh%varname(4)='thin_sensor'
         smesh%varname(5)='thickness'
      end block coarse_smesh
      
      ! Also create ensight output
      coarse_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='coarse')
         ! Add variables to output
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_surface('r2p',smesh)
      end block coarse_ensight
      
      ! Finally, create event for postprocessing
      pproc_event: block
         use param, only: param_read
         pproc_evt=event(time=dns%time,name='Postproc')
         call param_read('Postproc period',pproc_evt%tper)
      end block pproc_event

      ! Perform postprocessing
      if (pproc_evt%occurs()) call pproc()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Film drives overall time integration
      do while (.not.dns%time%done())

         ! Advance film simulation
         call dns%step()
         
         ! Perform postprocessing
         if (pproc_evt%occurs()) call pproc()
         
      end do
      
   end subroutine simulation_run


   !> Postprocess phasic moments from DNS
   !> Feed to coarse mesh for R2P reconstruction
   subroutine pproc()
      implicit none
      
      ! Integrate volume moments from fine to coarse mesh
      
      
      ! Integrate surface moments from fine to coarse mesh
      
      
      ! Perform R2P reconstruction on coarse mesh
      update_vf: block
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block update_vf

      ! Update smesh data
      update_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         smesh%var(1,:)=1.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                        smesh%var(2,np)=vf%curv2p(nplane,i,j,k)
                        smesh%var(3,np)=vf%edge_sensor(i,j,k)
                        smesh%var(4,np)=vf%thin_sensor(i,j,k)
                        smesh%var(5,np)=vf%thickness  (i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block update_smesh
      
      ! Perform ensight output
      call ens_out%write_data(dns%time%t)
      
   end subroutine pproc
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize film simulation
      call dns%final()
      
   end subroutine simulation_final
   
   
end module simulation
