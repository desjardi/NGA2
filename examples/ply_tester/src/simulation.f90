!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   implicit none
   private
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   
   public :: simulation_init,simulation_run,simulation_final
   
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use string, only: str_medium
      use param,  only: param_read
      implicit none
      character(len=str_medium) :: plyfile
      
      ! Read in ply file name
      call param_read('ply filename',plyfile)
      
      ! Create surface mesh
      smesh=surfmesh(plyfile=plyfile,nvar=0,name='ply')
      
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='plytester')
      
      ! Add surfmesh to output
      call ens_out%add_surface('ply',smesh)
      
      ! Output to ensight
      call ens_out%write_data(time=0.0_WP)
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
   end subroutine simulation_final
   
end module simulation
