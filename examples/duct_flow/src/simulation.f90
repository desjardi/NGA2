!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,             only: str_medium
   use precision,          only: WP
   use block1_class,       only: block1
   use event_class,        only: event
   use timetracker_class,  only: timetracker
   use datafile_class,     only: datafile
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> Block 1 object1
   type(block1) :: b1

contains


   !> Initialization of problem solver
   subroutine simulation_init
      use geometry, only: cfg1
      use param, only: param_read
      implicit none

      ! Blocks
      b1%cfg=>cfg1; call b1%init()
      
   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      integer :: i,j,k
      character(len=str_medium) :: dirname,timestamp

      ! Perform time integration 
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

   end subroutine simulation_final


end module simulation
