!> Simulation hooks for AMR Tester
!> Uses param module for test selection via input file or command line
module simulation
   use param, only: param_init, param_final, param_read, param_exists
   use messager, only: log, die
   implicit none

contains

   !> Initialization hook
   subroutine simulation_init()
      call log("Simulation Init: AMR Tester")
      call param_init()
   end subroutine simulation_init

   !> Run selected test(s)
   subroutine simulation_run()
      use mod_test_amrdata,   only: test_amrdata
      use mod_test_amrscalar, only: test_amrscalar
      use mod_test_amrmg,     only: test_amrmg
      use mod_test_project,   only: test_project
      use mod_test_amrincomp, only: test_amrincomp
      character(len=32) :: test_name

      call log("Simulation Run: AMR Tester")

      ! Read test selection from input file or command line
      ! Usage: ./amr_tester.*.exe --Test=amrmg
      !    or: ./amr_tester.*.exe -i input  (with "Test: amrmg" in input file)
      call param_read('Test', test_name, short='t', default='all')

      select case (trim(test_name))
       case ('amrdata')
         call test_amrdata()
       case ('amrscalar')
         call test_amrscalar()
       case ('amrmg')
         call test_amrmg()
       case ('project')
         call test_project()
       case ('amrincomp')
         call test_amrincomp()
       case ('all')
         call test_amrdata()
         call test_amrscalar()
         call test_amrmg()
         call test_project()
         call test_amrincomp()
       case default
         call die('Unknown test: '//trim(test_name)//'. Valid: amrdata, amrscalar, amrmg, project, amrincomp, all')
      end select
   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      call log("Simulation Final: Done")
      call param_final()
   end subroutine simulation_final

end module simulation
