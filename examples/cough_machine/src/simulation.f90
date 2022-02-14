!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,    only: WP
   use block1_class, only: block1
   use block2_class, only: block2
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Block 1 and 2 objects
   type(block1) :: b1
   type(block2) :: b2
   
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use geometry, only: cfg1,cfg2
      implicit none
      
      ! Initialize both blocks
      b1%cfg=>cfg1; call b1%init()
      b2%cfg=>cfg2; call b2%init()
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration - block 1 is the main driver here
      do while (.not.b1%time%done())
         
         ! Advance block 1
         call b1%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays for both blocks
      call b1%final()
      call b2%final()
      
   end subroutine simulation_final
   
   
end module simulation
