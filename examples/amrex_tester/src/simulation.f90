!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use amrex_amr_module, only: amrex_init,amrex_parallel_ioprocessor,amrex_finalize
      implicit none
      
      ! Play with some AMReX functionalities here
      call amrex_init()
      if (amrex_parallel_ioprocessor()) print*,'Hello from AMReX!'
      call amrex_finalize()
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
   end subroutine simulation_final
   
   
end module simulation
