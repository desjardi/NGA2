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
      use param,     only: param_read
      use precision, only: WP
      use ils_class, only: rbgs
      implicit none
      
      ! Create an incompressible flow solver
      create_solver: block
         ! Create solver
         fs=incomp(cfg,'Bob')
         ! Assign constant fluid properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',fs%visc)
         ! Configure pressure solver
         fs%psolv%maxit=100
         fs%psolv%acvg=1.0e-4_WP
         fs%psolv%rcvg=1.0e-4_WP
         fs%psolv%method=rbgs
         
         ! Check solver objects
         call fs%print()
         call fs%psolv%print()
         
      end block create_solver
      
      
      
   end subroutine simulation_init
   
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      
      
      
   end subroutine simulation_run
   
   
end module simulation
