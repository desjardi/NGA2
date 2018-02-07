!> This is a typical NGA driver.
!> @version 1.0
!> @author O. Desjardins
program nga
  use random  , only: random_init,random_final
  use option  , only: option_init,option_final,option_isinvoked,option_read,option_print
  use parser  , only: parser_init,parser_final
  use parallel, only: parallel_init,parallel_final,amroot
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  logical :: verbose !< Verbosity of this run
  logical :: is_there
  
  ! First initialize a basic parallel environment
  call parallel_init
  
  ! Framework initialization ==============
  ! Initialize command line parser
  call option_init
  ! Check if help option is invoked
  call option_isinvoked(is_there,shn='h',lgn='help')
  if (is_there) then
     if (amroot) then
        write(output_unit,'(a)') ' ____________________________________________________________ '
        write(output_unit,'(a)') '|                    _   _  ____    _                        |'
        write(output_unit,'(a)') '|                   | \ | |/ ___|  / \                       |'
        write(output_unit,'(a)') '|                   |  \| | |  _  / _ \                      |'
        write(output_unit,'(a)') '|                   | |\  | |_| |/ ___ \                     |'
        write(output_unit,'(a)') '|                   |_| \_|\____/_/   \_\                    |'
        write(output_unit,'(a)') '|                                                            |'
        write(output_unit,'(a)') '| Multiphase Reactive Turbulent Flow Solver by O. Desjardins |'
        write(output_unit,'(a)') '|                   ---------------------                    |'
        write(output_unit,'(a)') '|     The following format is expected when calling NGA:     |'
        write(output_unit,'(a)') '|                                                            |'
        write(output_unit,'(a)') '|     > nga.().exe  [-X value] [--XXX=value]  filename       |'
        write(output_unit,'(a)') '|       executable | short and long options | input file     |'
        write(output_unit,'(a)') '|                  |-- value is optional! --|                |'
        write(output_unit,'(a)') '|____________________________________________________________|'
     end if
     ! Clean exit
     call parallel_final; stop
  end if
  ! Check if we are verbose
  call option_isinvoked(verbose,shn='v',lgn='verbose')
  if (verbose.and.amroot) call option_print
  
  ! Initialize file parser
  call parser_init
  
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
  !call timing_final
  !call monitor_final
  
  ! Clean up RNG
  call random_final

  ! Clean up file parser
  call parser_final
  
  ! Clean up options
  call option_final
  
  ! =======================================
  
  ! Terminate the basic parallel environment
  call parallel_final
  
end program nga
