!> Various definitions and tools for running an NGA2 simulation
module simulation
   !use beamimpact_class,   only: beamimpact
   !use fracturetest_class, only: fracturetest
   use object_class,  only: object
   use channel_class, only: channel
   implicit none
   private
   
   !> Solid simulation
   !type(beamimpact) :: solid
   !type(fracturetest) :: solid
   type(object) :: solid
   
   !> Fluid simulation
   type(channel) :: fluid
   
   public :: simulation_init,simulation_run,simulation_final
   
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize solid simulation
      call solid%init()

      ! Initialize fluid simulation
      call fluid%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Solid drives overall time integration
      do while (.not.solid%time%done())
         
         ! Advance solid simulation
         call solid%step()
         
         ! Advance fluid simulation
         call fluid%step(ls=solid%ls)
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize solid simulation
      call solid%final()

      ! Finalize fluid simulation
      call fluid%final()
      
   end subroutine simulation_final
   
   
end module simulation
