!> This is a typical NGA driver.
!> @version 1.0
!> @author O. Desjardins
program nga
   use random,     only: random_initialize
   use param,      only: param_init,param_final
   use parallel,   only: parallel_init,parallel_final
   use messager,   only: messager_init,messager_final
   use geometry,   only: geometry_init
   use simulation, only: simulation_init,simulation_run,simulation_final
   implicit none
   
   
   ! =======================================
   ! Initialization sequence ===============
   ! =======================================
   ! Initialize parallel environment =======
   call parallel_init
   ! Initialize messaging capabilities ====
   call messager_init
   ! Initialize user interaction ===========
   call param_init
   ! Initialize random number generator ====
   call random_initialize
   ! =======================================
   
   
   
   
   ! =======================================
   ! Setup the problem grid ================
   ! =======================================
   call geometry_init
   ! =======================================
   
   
   
   
   
   
   ! =======================================
   ! Setup the solver ======================
   ! =======================================
   call simulation_init
   ! =======================================
   
   
   
   
   ! =======================================
   ! Run the solver ========================
   ! =======================================
   call simulation_run
   ! =======================================
   
   
   
   
   ! =======================================
   ! Finalize the solver ===================
   ! =======================================
   call simulation_final
   ! =======================================
   
   
   
   
   ! Code termination ======================
   ! Terminate solvers, data, and geometry
   !call simulation_final
   !call data_final
   !call geometry_final
   
   ! Terminate I/O, parser, timing, and monitor
   !call io_final
   
   
   
   
   
   ! =======================================
   ! Termination sequence ==================
   ! =======================================
   ! Clean up user parameters ==============
   call param_final
   ! Clean up monitoring tools =============
   call messager_final
   ! Clean parallel exit ===================
   call parallel_final
   ! =======================================
   
   
end program nga
