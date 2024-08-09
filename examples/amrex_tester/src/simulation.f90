!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use amrcore_class, only: amrcore
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Give ourselves an amrcore object
   type(amrcore) :: amr

contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      ! Create an amrcore object
      call amr%initialize()
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      use amrcore_class, only: finalize_amrex
      implicit none
      call amr%finalize()
      call finalize_amrex()
   end subroutine simulation_final
   
   
end module simulation
