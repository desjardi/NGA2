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
   real(WP), dimension(:,:,:), allocatable :: dudt,dvdt,dwdt
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
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
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%acvg); fs%psolv%rcvg=fs%psolv%acvg
         ! Initialize solver
         call fs%psolv%init_solver(amg)
         ! Check solver objects
         call fs%print()
      end block create_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dudt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(dvdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(dwdt(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize boundary conditions
      initialize_bc: block
         call fs%add_bcond('inflow',dirichlet,yplus_locator)
      end block initialize_bc
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(fs%cfg%amRoot)
         call param_read('Time step size',time%dtmax)
         time%dt=time%dtmax
      end block initialize_timetracker
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         fs%U=0.0_WP
         fs%V=1.0_WP
         fs%W=0.0_WP
         call fs%apply_bcond(time%t)
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
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
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('dveldt',dudt,dvdt,dwdt)
         call ens_out%add_scalar('div',fs%div)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Create monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         ! Add time info first
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         ! Write it out
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call mfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform explicit Euler time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply BCs
         call fs%apply_bcond(time%t)
         
         ! Explicit Euler advancement of NS
         call fs%get_dmomdt(dudt,dvdt,dwdt)
         fs%U=fs%U+time%dt*dudt/fs%rho
         fs%V=fs%V+time%dt*dvdt/fs%rho
         fs%W=fs%W+time%dt*dwdt/fs%rho
         
         ! Solve Poisson equation
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         
         ! Correct velocity
         call fs%get_pgrad(fs%psolv%sol,dudt,dvdt,dwdt)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*dudt/fs%rho
         fs%V=fs%V-time%dt*dvdt/fs%rho
         fs%W=fs%W-time%dt*dwdt/fs%rho
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
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
      deallocate(dudt,dvdt,dwdt,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
