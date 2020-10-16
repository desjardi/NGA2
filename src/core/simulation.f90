!> Various definitions and tools for running an NGA2 simulation
module simulation
   use incomp_class, only: incomp
   use geometry,     only: cfg
   implicit none
   private
   
   !> Single flow solver
   type(incomp), public :: fs
   
   public :: simulation_init,simulation_run
   
contains
   
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Create an incompressible flow solver
      create_solver: block
         ! Create solver
         fs=incomp(cfg,'Bob')
         ! Assign constant fluid properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',fs%visc)
         ! Prepare metrics
         
         ! Check solver object
         call fs%print()
      end block create_solver
      
      
      
   end subroutine simulation_init
   
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      
      
      
   end subroutine simulation_run
   
   
end module simulation
