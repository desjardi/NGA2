!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use sml_class, only: sml
   implicit none
   private
   
   !> SML simulation
   type(sml) :: mixing
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize mixing layer
      call mixing%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.mixing%time%done())

         ! Advance mixing layer simulation
         call mixing%step()

      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize mixing layer simulation
      call mixing%final()
      
   end subroutine simulation_final
   
   
end module simulation