!> This is a typical NGA driver.
!> @version 1.0
!> @author O. Desjardins
program nga
  use random, only: random_init,random_final
  use option, only: option_init,option_final
  use parallel, only: parallel_init,parallel_final
  implicit none
  
  ! Code initialization ===================
  
  ! Initialize basic parallel environment
  call parallel_init
  
  ! Initialize command line options
  call option_init
  
  ! Initialize RNG
  call random_init
  
  ! Initialization
  !call monitor_init
  !call timing_init
  !call parser_init
  !call io_init
  
  ! Initialize geometry, data, and solvers
  !call geometry_init
  !call data_init
  !call simulation_init


  
  ! =======================================
  
  
  
  ! Run flow solver
  !call simulation_run
  

  
  
  ! Code termination ======================
  ! Terminate solvers, data, and geometry
  !call simulation_final
  !call data_final
  !call geometry_final
  
  ! Terminate I/O, parser, timing, and monitor
  !call io_final
  !call parser_final
  !call timing_final
  !call monitor_final
  
  ! Clean up RNG
  call random_final
  
  ! Clean up options
  call option_final
  
  ! Terminate basic parallel environment
  call parallel_final
  
  ! =======================================
  
end program nga
