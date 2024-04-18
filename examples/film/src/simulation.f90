!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,  only: WP
   use film_class, only: film
   implicit none
   private
   
   !> Film simulation
   type(film) :: dns
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: hit_group
      
      ! Initialize film simulation
      call dns%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Film drives overall time integration
      do while (.not.dns%time%done())

         ! Advance film simulation
         call dns%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize film simulation
      call dns%final()
      
   end subroutine simulation_final
   
   
end module simulation
