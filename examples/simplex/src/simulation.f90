!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use simplex_class, only: simplex
   implicit none
   private
   
   !> Simplex simulation
   type(simplex) :: spx
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize simplex simulation
      call spx%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Simplex drives overall time integration
      do while (.not.spx%time%done())
         ! Advance simplex simulation
         call spx%step()
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize simplex simulation
      call spx%final()
      
   end subroutine simulation_final
   
   
end module simulation