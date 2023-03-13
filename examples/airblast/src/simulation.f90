!> Various definitions and tools for running an NGA2 simulation
module simulation
   !use nozzle_class, only: nozzle
   use atom_class,   only: atom
   implicit none
   private
   
   !> Injector simulation
   !type(nozzle) :: injector
   
   !> Atomization simulation
   type(atom) :: atomization
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize injector simulation
      !call injector%init()

      ! Initialize atomization simulation
      call atomization%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Perform global time integration
      !do while (.not.injector%time%done())
      do while (.not.atomization%time%done())
         
         ! Advance injector simulation
         !call injector%step()
         
         ! Advance atomization simulation
         call atomization%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize injector simulation
      !call injector%final()
      
      ! Finalize atomization simulation
      call atomization%final()
      
   end subroutine simulation_final
   

end module simulation
