!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single low Mach flow solver and scalar solver and corresponding time tracker
   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Create an incompressible flow solver with bconds
      create_solver: block
         use ils_class,     only: rbgs,amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,smg,pfmg
         use lowmach_class, only: dirichlet,convective,neumann,clipped_neumann
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=pfmg)
      end block create_solver
      
      
      ! Create a scalar solver
      create_scalar: block
         use ils_class,      only: gmres
         use vdscalar_class, only: dirichlet,neumann,quick
         ! Create scalar solver
         sc=vdscalar(cfg=cfg,scheme=quick,name='Temperature')
         ! Configure implicit scalar solver
         sc%implicit%maxit=fs%implicit%maxit; sc%implicit%rcvg=fs%implicit%rcvg
         ! Setup the solver
         call sc%setup(implicit_ils=gmres)
      end block create_scalar
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         ! Flow solver
         allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         ! Scalar solver
         allocate(resSC(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=fs%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         type(bcond), pointer :: inflow
         integer :: n,i,j,k
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%rho=1.0_WP
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity
      
      
      ! Initialize our temperature field
      initialize_scalar: block
         use vdscalar_class, only: bcond
         type(bcond), pointer :: inflow
         integer :: n,i,j,k
         ! Zero initial field
         sc%SC=0.0_WP
         ! Apply all other boundary conditions
         call sc%apply_bcond(time%t,time%dt)
      end block initialize_scalar
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='test')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('temperature',sc%SC)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(sc%SCmax,'Tmax')
         call mfile%add_column(sc%SCmin,'Tmin')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
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
         
         ! Remember old scalar
         sc%rhoold=sc%rho
         sc%SCold =sc%SC
         
         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ============= SCALAR SOLVER =======================
            ! Build mid-time scalar
            sc%SC=0.5_WP*(sc%SC+sc%SCold)
            
            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Assemble explicit residual
            resSC=-2.0_WP*(sc%rho*sc%SC-sc%rhoold*sc%SCold)+time%dt*resSC
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Apply this residual
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! SOMEWHERE WE NEED TO UPDATE THE DENSITY AND THE VISCOSITY AND THE DIFFUSIVITY
            
            ! ============ VELOCITY SOLVER ======================
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rhoU-fs%rhoUold)+time%dt*resU
            resV=-2.0_WP*(fs%rhoV-fs%rhoVold)+time%dt*resV
            resW=-2.0_WP*(fs%rhoW-fs%rhoWold)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dt*resU
            fs%rhoV=fs%rhoV-time%dt*resV
            fs%rhoW=fs%rhoW-time%dt*resW
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call mfile%write()
         call cflfile%write()
         
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
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
