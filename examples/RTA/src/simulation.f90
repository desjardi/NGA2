!> Various definitions and tools for running an NGA2 simulation
module simulation
   use hit_class, only: hit
   !use rta_class, only: rta
   implicit none
   private
   
   !> HIT simulation
   type(hit) :: turb
   
   !> Rayleigh-Taylor atomization simulation
   !type(rta) :: atom
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use parallel, only: group
      implicit none
      
      ! Initialize atomization simulation
      !call atom%init()

      ! Initialize HIT simulation
      call turb%init(group=group)
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      !do while (.not.atom%time%done())
      do while (.not.turb%time%done())
         
         ! Advance atomization simulation
         !call atom%step()
         
         ! Advance HIT simulation
         call turb%step(turb%time%dt)
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize atomization simulation
      !call atom%final()
      
      ! Finalize HIT simulation
      call turb%final()
      
   end subroutine simulation_final
   

end module simulation