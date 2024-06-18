!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,             only: WP
   use film_retraction_class, only: film_retraction
   implicit none
   private
   
   !> Film retraction simulation
   type(film_retraction) :: film
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize film simulation
      call film%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Film drives overall time integration
      do while (.not.film%time%done())
         
         ! Advance film simulation
         call film%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize film simulation
      call film%final()
      
   end subroutine simulation_final
   
   
end module simulation