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
   type(monitor) :: mfile,cflfile,consfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Equation of state
   real(WP) :: Tmin,Tmax,fluid_mass
   
contains
   
   
   !> Define here our equation of state - rho(T,mass)
   subroutine get_rho(mass)
      implicit none
      real(WP), intent(in) :: mass
      integer :: i,j,k
      real(WP) :: one_over_T
      ! Integrate 1/T
      resSC=1.0_WP/sc%SC
      call sc%cfg%integrate(resSC,integral=one_over_T)
      ! Calculate density
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               sc%rho(i,j,k)=mass/(sc%SC(i,j,k)*one_over_T)
            end do
         end do
      end do
   end subroutine get_rho
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read in the EOS info
      !call param_read('Min density',rho_min)
      !call param_read('Max density',rho_max)
      
      ! Create an incompressible flow solver with bconds
      create_solver: block
         use ils_class, only: pcg_amg,pfmg
         real(WP) :: visc
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
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
         real(WP) :: diffusivity
         ! Create scalar solver
         sc=vdscalar(cfg=cfg,scheme=quick,name='Temperature')
         ! Assign constant diffusivity
         call param_read('Dynamic diffusivity',diffusivity)
         sc%diff=diffusivity
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
      
      
      ! Initialize our temperature field
      initialize_scalar: block
         integer :: j
         ! Read in the temperature and mass
         call param_read('Min temperature',Tmin)
         call param_read('Max temperature',Tmax)
         call param_read('Total mass',fluid_mass)
         ! Stratified initial field
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            sc%SC(:,j,:)=Tmin+0.5_WP*(tanh(-fs%cfg%ym(j)/0.05_WP)+1.0_WP)*(Tmax-Tmin)
         end do
         ! Build corresponding density
         call get_rho(mass=fluid_mass)
      end block initialize_scalar
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Set density from scalar
         fs%rho=sc%rho
         ! Form momentum
         call fs%rho_multiply
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         call fs%interp_vel(Ui,Vi,Wi)
         resSC=0.0_WP
         call fs%get_div(drhodt=resSC)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity
      
      
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
         call ens_out%add_scalar('density',sc%rho)
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
         call sc%get_int()
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
         call mfile%add_column(sc%rhomax,'RHOmax')
         call mfile%add_column(sc%rhomin,'RHOmin')
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
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(sc%SCint,'SC integral')
         call consfile%add_column(sc%rhoint,'RHO integral')
         call consfile%add_column(sc%rhoSCint,'rhoSC integral')
         call consfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
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
            resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Apply this residual
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ UPDATE PROPERTIES ====================
            ! Backup rhoSC
            resSC=sc%rho*sc%SC
            ! Update density
            call get_rho(mass=fluid_mass)
            ! Rescale scalar for conservation
            !sc%SC=resSC/sc%rho
            ! UPDATE THE VISCOSITY
            ! UPDATE THE DIFFUSIVITY
            ! ===================================================
            
            ! ============ VELOCITY SOLVER ======================
            
            ! Build n+1 density
            fs%rho=0.5_WP*(sc%rho+sc%rhoold)
            
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions and update momentum
            call fs%apply_bcond(time%tmid,time%dtmid)
            call fs%rho_multiply()
            call fs%apply_bcond(time%tmid,time%dtmid)
            
            ! Solve Poisson equation
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%correct_mfr(drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resSC)
         call fs%get_div(drhodt=resSC)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
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
