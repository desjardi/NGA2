!> Various definitions and tools for running an NGA2 simulation
module simulation
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
   
   public :: simulation_init,simulation_run
   
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
         fs%V=0.0_WP
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
         call ens_out%add_scalar('Pressure',fs%P)
         call ens_out%add_vector('Velocity',fs%U,fs%V,fs%W)
      end block create_ensight
      
      
      ! Try to use the pressure solver
      test_pressure_solver: block
         ! Create a scaled RHS and output it
         fs%psolv%rhs=0.0_WP
         if (fs%cfg%jproc.eq.         1) fs%psolv%rhs(:,fs%cfg%jmin_,:)=+1.0_WP
         if (fs%cfg%jproc.eq.fs%cfg%npy) fs%psolv%rhs(:,fs%cfg%jmax_,:)=-1.0_WP
         fs%psolv%rhs=-fs%cfg%vol*fs%psolv%rhs
         call ens_out%add_scalar('RHS',fs%psolv%rhs)
         ! Set initial guess to zero
         fs%psolv%sol=0.0_WP
         ! Call the solver
         call fs%psolv%solve()
         ! Copy back to pressure
         fs%P=fs%psolv%sol
         ! Output to ensight
         !if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block test_pressure_solver
      
      
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
      use precision, only: WP
      implicit none
      real(WP), dimension(:,:,:), allocatable :: dudt,dvdt,dwdt
      
      ! Allocate work arrays
      allocate(dudt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(dvdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(dwdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      
      call ens_out%add_vector('dveldt',dudt,dvdt,dwdt)
      
      ! Perform explicit Euler time integration
      do while (.not.time%done())
         ! Increment time
         call time%adjust_dt()
         call time%increment()
         ! Evaluate velocity rate of change
         call fs%get_dmomdt(dudt,dvdt,dwdt)
         ! Explicit Euler advancement
         fs%U=fs%U+time%dt*dudt/fs%rho
         fs%V=fs%V+time%dt*dvdt/fs%rho
         fs%W=fs%W+time%dt*dwdt/fs%rho
         ! Output to ensight
         !if (ens_evt%occurs()) call ens_out%write_data(time%t)
         call ens_out%write_data(time%t)
         ! Write out monitor file
         call mfile%write()
         
         stop
         
      end do
      
      ! Deallocate work arrays
      deallocate(dudt,dvdt,dwdt)
      
   end subroutine simulation_run
   
   
end module simulation
