!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp,dirichlet
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single incompressible flow solver and corresponding time tracker
   type(incomp),      public :: fs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: dudt,dvdt,dwdt,div
   
   
contains
   
   !> Function that localizes the top of the domain
   function yplus_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yplus_locator
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use precision, only: WP
      use param,     only: param_read
      use ils_class, only: rbgs,amg
      implicit none
      
      
      ! Allocate work arrays
      allocate(dudt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(dvdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(dwdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(div (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      
      
      ! Create an incompressible flow solver
      create_solver: block
         ! Create solver
         fs=incomp(cfg,'Bob')
         ! Assign constant fluid properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',fs%visc)
         ! Configure pressure solver
         fs%psolv%maxit=100
         fs%psolv%acvg=1.0e-4_WP
         fs%psolv%rcvg=1.0e-4_WP
         ! Initialize solver
         call fs%psolv%init_solver(amg)
         ! Check solver objects
         call fs%print()
      end block create_solver
      
      
      ! Initialize boundary conditions
      initialize_bc: block
         call fs%add_bcond('inflow',dirichlet,yplus_locator)
      end block initialize_bc
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(fs%cfg%amRoot)
         call param_read('Time step size',time%dt)
      end block initialize_timetracker
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         fs%U=0.0_WP
         fs%V=1.0_WP
         fs%W=0.0_WP
      end block initialize_velocity
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg,'test')
         ! Create event for Ensight output
         ens_evt=event(time,'Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
         call ens_out%add_vector('dveldt',dudt,dvdt,dwdt)
         call ens_out%add_scalar('div',div)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Create monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         ! Add time info first
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         ! Write it out
         call mfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform explicit Euler time integration
      do while (.not.time%done())
         
         ! Increment time
         call time%adjust_dt()
         call time%increment()
         
         ! Explicit Euler advancement of NS
         call fs%get_dmomdt(dudt,dvdt,dwdt)
         fs%U=fs%U+time%dt*dudt/fs%rho
         fs%V=fs%V+time%dt*dvdt/fs%rho
         fs%W=fs%W+time%dt*dwdt/fs%rho
         
         ! Solve Poisson equation
         call fs%get_div(div)
         fs%psolv%rhs=-fs%cfg%vol*div*fs%rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         
         ! Correct velocity
         call fs%get_pgrad(fs%psolv%sol,dudt,dvdt,dwdt)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*dudt/fs%rho
         fs%V=fs%V-time%dt*dvdt/fs%rho
         fs%W=fs%W-time%dt*dwdt/fs%rho
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Write out monitor file
         call mfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(dudt,dvdt,dwdt,div)
      
   end subroutine simulation_final
   
   
end module simulation
