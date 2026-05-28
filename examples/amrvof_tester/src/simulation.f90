!> Simulation hooks for AMRVOF Tester
module simulation
   use param, only: param_read
   use messager, only: log, die
   implicit none

contains

   !> Initialization hook
   subroutine simulation_init()
      call log("Simulation Init: AMRVOF Tester")
   end subroutine simulation_init

   !> Run selected test(s)
   subroutine simulation_run()
      use mod_test_geometry, only: test_geometry
      use mod_test_amrvof,   only: test_amrvof
      character(len=32) :: test_name

      call log("Simulation Run: AMRVOF Tester")

      ! Read test selection from input file
      call param_read('Test', test_name, short='t', default='geometry')

      select case (trim(test_name))
       case ('geometry')
         call test_geometry()
       case ('advect')
         call test_amrvof()
       case ('all')
         call test_geometry()
         call test_amrvof()
       case default
         call die('Unknown test: '//trim(test_name)//'. Valid: geometry, advect, all')
      end select
   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      call log("Simulation Final: Done")
   end subroutine simulation_final

end module simulation
