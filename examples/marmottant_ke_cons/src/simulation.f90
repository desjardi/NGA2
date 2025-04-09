!> Various definitions and tools for running an NGA2 simulation
module simulation
   use coaxialjet_class, only: coaxialjet
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final
   
   !> Jet simulation
   type(coaxialjet) :: jet
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize jet simulation
      call jet%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Jet drives overall time integration
      do while (.not.jet%time%done())
         
         ! Advance jet simulation
         call jet%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize jet simulation
      call jet%final()
      
   end subroutine simulation_final
   
   
end module simulation
