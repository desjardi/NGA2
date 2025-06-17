!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use hypre_str_class,   only: hypre_str
   use compress_class,    only: compress
   use energy_class,      only: energy
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Single compressible flow solver and scalar solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(ddadi),       public :: ss
   type(compress),    public :: fs
   type(energy),      public :: sc
   type(timetracker), public :: time
   
   !> SGS modeling
   logical        :: use_sgsb
   type(sgsmodel) :: sgsb
   logical        :: use_sgss
   type(sgsmodel) :: sgss
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   real(WP) :: RHOcvg=0.0_WP,RHOtol=1.0e-4_WP
   
   !> Monitoring of conservation
   type(monitor) :: consfile
   real(WP) :: rhoKEint,dkdt,dedt,sdiss,ddiss,sdiss_sgss,ddiss_sgss,ddiss_sgsb,pdil
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resE
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,C2,Ma
   
   !> Equation of state and flow conditions
   real(WP) :: Gamma,Mach,Reynolds,Prandtl,bulk2shear
   
contains
   
   
   !> Obtain density from equation of state
   subroutine get_rho()
      implicit none
      integer :: i,j,k
      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               fs%RHO(i,j,k)=fs%P(i,j,k)/(sc%E(i,j,k)*(Gamma-1.0_WP))
            end do
         end do
      end do
   end subroutine get_rho
   
   
   !> Obtain speed of sound squared from equation of state
   subroutine get_c2()
      implicit none
      integer :: i,j,k
      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               C2(i,j,k)=Gamma*fs%P(i,j,k)/fs%RHO(i,j,k)
            end do
         end do
      end do
   end subroutine get_c2
   
   
   !> Sutherland's law for viscosity as a function of temperature
   subroutine get_visc()
      implicit none
      integer :: i,j,k
      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               fs%viscs(i,j,k)=Reynolds**(-1.0_WP)*(1.4042_WP*(sc%E(i,j,k)/sc%Cv)**1.5_WP)/(sc%E(i,j,k)/sc%Cv+0.4042_WP)
            end do
         end do
      end do
   end subroutine get_visc
   
   
   !> Post-process conservation properties
   subroutine analyze_conservation()
      implicit none
      ! Get integral of KE and dk/dt
      calculate_KE: block
         integer :: i,j,k
         do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
            resU(i,j,k)=0.5_WP*sum(fs%itpu_x(:,i,j,k)*(fs%sRHOX(i:i+1,j,k)*fs%U(i:i+1,j,k))**2)&
            &          +0.5_WP*sum(fs%itpv_y(:,i,j,k)*(fs%sRHOY(i,j:j+1,k)*fs%V(i,j:j+1,k))**2)&
            &          +0.5_WP*sum(fs%itpw_z(:,i,j,k)*(fs%sRHOZ(i,j,k:k+1)*fs%W(i,j,k:k+1))**2)
            resV(i,j,k)=0.5_WP*sum(fs%itpu_x(:,i,j,k)*(fs%sRHOXold(i:i+1,j,k)*fs%Uold(i:i+1,j,k))**2)&
            &          +0.5_WP*sum(fs%itpv_y(:,i,j,k)*(fs%sRHOYold(i,j:j+1,k)*fs%Vold(i,j:j+1,k))**2)&
            &          +0.5_WP*sum(fs%itpw_z(:,i,j,k)*(fs%sRHOZold(i,j,k:k+1)*fs%Wold(i,j,k:k+1))**2)
         end do; end do; end do
         call cfg%integrate(resU,integral=rhoKEint)
         call cfg%integrate(resV,integral=dkdt); dkdt=(rhoKEint-dkdt)/time%dt
      end block calculate_KE
      ! Get de/dt
      calculate_e: block
         real(WP) :: rhoe,rhoeold
         resU=fs%RHO*sc%E; call cfg%integrate(resU,integral=rhoe)
         resU=fs%RHOold*sc%Eold; call cfg%integrate(resU,integral=rhoeold)
         dedt=(rhoe-rhoeold)/time%dt
      end block calculate_e
      ! Calculate dissipations
      calculate_dissipations: block
         ! Get resolved+sgs dilatational and solenoidal dissipations
         call fs%get_dilatational_dissipation(diss=resW); call cfg%integrate(resW,integral=ddiss)
         call fs%get_solenoidal_dissipation  (diss=resW); call cfg%integrate(resW,integral=sdiss)
         ! Remember viscosities
         resU=fs%viscs
         resV=fs%viscb
         ! Get sgs dilatational dissipation from sgsb
         if (use_sgsb) then
            fs%viscs=0.0_WP; fs%viscb=sgsb%visc
            call fs%get_dilatational_dissipation(diss=resW); call cfg%integrate(resW,integral=ddiss_sgsb)
         else
            ddiss_sgsb=0.0_WP
         end if
         ! Get sgs dilatational and solenoidal dissipations from sgss
         if (use_sgss) then
            fs%viscs=sgss%visc; fs%viscb=0.0_WP
            call fs%get_dilatational_dissipation(diss=resW); call cfg%integrate(resW,integral=ddiss_sgss)
            call fs%get_solenoidal_dissipation  (diss=resW); call cfg%integrate(resW,integral=sdiss_sgss)
         else
            ddiss_sgss=0.0_WP
            sdiss_sgss=0.0_WP
         end if
         ! Compute resolved dissipations
         ddiss=ddiss-ddiss_sgsb-ddiss_sgss
         sdiss=sdiss-sdiss_sgss
         ! Set viscosities back
         fs%viscs=resU
         fs%viscb=resV
      end block calculate_dissipations
      ! Calculate pressure dilatation
      calculate_pdil: block
         resU=fs%P; call fs%get_pdil(P=resU,Pdil=resV); call cfg%integrate(resV,integral=pdil)
      end block calculate_pdil
   end subroutine analyze_conservation
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Prepare EoS and flow conditions
      initialize_eos: block
         call param_read('Gamma'   ,gamma   )
         call param_read('Mach'    ,Mach    )
         call param_read('Reynolds',Reynolds)
         call param_read('Prandtl' ,Prandtl )
         call param_read('Bulk to shear ratio',bulk2shear)
      end block initialize_eos
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resE(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(C2  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ma  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax)
         time%itmin=2
         call param_read('Relaxation coeff',time%relax,default=1.0_WP)
      end block initialize_timetracker
      
      ! Create a compressible flow solver
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Compressible NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=18
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_velocity_solver
      
      ! Create an energy solver
      create_energy: block
         real(WP) :: Pr
         ! Create energy solver
         call sc%initialize(cfg=cfg,name='Energy')
         ! Set Cv such that T=1
         sc%Cv=1.0_WP/(Gamma*(Gamma-1.0_WP)*Mach**2)
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Energy',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_energy
      
      ! Initialize our initial conditions
      initial_conditions: block
         use mathtools, only: Pi,twoPi
         integer :: i,j,k
         ! Initialize density field using tanh and corresponding energy from EOS
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  fs%P(i,j,k)=1.0_WP/(Gamma*Mach**2)+(cos(2.0_WP*fs%cfg%xm(i))+cos(2.0_WP*fs%cfg%ym(j)))*(cos(2.0_WP*fs%cfg%zm(k))+2.0_WP)/16.0_WP
                  sc%E(i,j,k)=1.0_WP/((Gamma-1.0_WP)*Gamma*Mach**2)
                  fs%RHO(i,j,k)=fs%P(i,j,k)/(sc%E(i,j,k)*(Gamma-1.0_WP))
                  fs%U(i,j,k)=+sin(fs%cfg%x (i))*cos(fs%cfg%ym(j))*cos(fs%cfg%zm(k))
                  fs%V(i,j,k)=-cos(fs%cfg%xm(i))*sin(fs%cfg%y (j))*cos(fs%cfg%zm(k))
                  fs%W(i,j,k)=0.0_WP
               end do
            end do
         end do
         ! Get face RHO
         call fs%update_faceRHO()
         sc%Eold=sc%E
         fs%RHOold=fs%RHO
         fs%sRHOXold=fs%sRHOX
         fs%sRHOYold=fs%sRHOY
         fs%sRHOZold=fs%sRHOZ
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         ! Calculate cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
         ! Compute local Mach number
         call get_c2(); Ma=sqrt((Ui**2+Vi**2+Wi**2)/C2)
      end block initial_conditions
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use bSGS model',use_sgsb); if (use_sgsb) sgsb=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         call param_read('Use sSGS model',use_sgss); if (use_sgss) sgss=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         if (use_sgsb.or.use_sgss) allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block create_sgs
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='TaylorGreen')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('energy',sc%E)
         if (use_sgsb) call ens_out%add_scalar('sgsb',sgsb%visc)
         if (use_sgss) call ens_out%add_scalar('sgss',sgss%visc)
         call ens_out%add_scalar('Mach',Ma)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max(fs%RHO)
         call analyze_conservation()
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
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(sc%Emax,'Emax')
         call mfile%add_column(sc%Emin,'Emin')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(   RHOcvg,'RHOcvg')
         call mfile%add_column(  time%it,'Subiter')
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
         call consfile%add_column(cfg%vol_total,'Volume')
         call consfile%add_column(fs%RHOint,'Mass')
         call consfile%add_column(sc%rhoEint,'Internal energy')
         call consfile%add_column(rhoKEint,'Kinetic energy')
         call consfile%add_column(dkdt,'dkdt')
         call consfile%add_column(dedt,'dedt')
         call consfile%add_column(sdiss,'sdiss')
         call consfile%add_column(ddiss,'ddiss')
         call consfile%add_column(ddiss_sgsb,'ddiss_sgsb')
         call consfile%add_column(sdiss_sgss,'sdiss_sgss')
         call consfile%add_column(ddiss_sgss,'ddiss_sgss')
         call consfile%add_column(pdil,'pdil')
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
         
         ! Remember old energy, velocity, and density
         sc%Eold=sc%E
         fs%RHOold=fs%RHO
         fs%sRHOXold=fs%sRHOX
         fs%sRHOYold=fs%sRHOY
         fs%sRHOZold=fs%sRHOZ
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         
         ! ============= RECOMPUTE VISCOSITY =================
         call get_visc(); fs%viscb=bulk2shear*fs%viscs
         ! ===================================================
         
         ! ============= SUBGRID SCALE MODELING ==============
         if (use_sgsb.or.use_sgss) then
            sgs_modeling: block
               use sgsmodel_class, only: vreman,localartif
               call fs%get_gradU(gradU)
               if (use_sgsb) then
                  call sgsb%get_visc(type=localartif,dt=time%dt,rho=fs%RHO,gradu=gradU); fs%viscb=fs%viscb+sgsb%visc
               end if
               if (use_sgss) then
                  call sgss%get_visc(type=vreman    ,dt=time%dt,rho=fs%RHO,gradu=gradU); fs%viscs=fs%viscs+sgss%visc
               end if
            end block sgs_modeling
         end if
         ! ===================================================
         
         ! ============= RECOMPUTE DIFFUSIVITY ===============
         sc%diff=sc%Cv*Gamma*fs%viscs/Prandtl
         ! ===================================================
         
         ! ============= PREPARE WENO SCHEMES ================
         call fs%prepare_weno(rho=fs%RHO)
         call sc%prepare_weno(e=sc%E)
         ! ===================================================

         ! ============= INITIAL GUESS FOR RHO ===============
         call fs%predict_rho(dt=time%dt)
         ! ===================================================
         
         ! Perform sub-iterations until RHO is sufficiently converged
         RHOcvg=huge(1.0_WP); time%it=0
         do while (RHOcvg.gt.RHOtol.and.time%it.lt.time%itmax.or.time%it.lt.time%itmin)
            
            ! ============= ENERGY SOLVER =======================
            ! Explicit calculation of drhoE/dt from energy equation
            call sc%get_drhoEdt(resE,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Add pressure dilatation term
            call fs%get_pdil(P=fs%P,Pdil=resV); resE=resE+resV
            
            ! Add viscous heating term
            call fs%get_visc_heating(visc_heating=resU); resE=resE+resU
            
            ! Assemble explicit residual
            resE=time%dt*resE-(fs%RHO*sc%E-fs%RHOold*sc%Eold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resE,fs%RHO,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Apply this residual
            sc%E=sc%E+resE*time%relax
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ UPDATE DENSITY AND C2 ================
            ! Remember RHO
            resE=fs%RHO
            
            ! Calculate RHO predictor
            call get_rho(); call fs%update_faceRHO()
            
            ! Compute speed of sound squared
            call get_c2()
            ! ===================================================
            
            ! ============ VELOCITY SOLVER ======================
            ! Compute Umid
            call fs%get_Umid()
            
            ! Get rhoU/rhoV/rhoW
            call fs%rho_multiply()
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=time%dt*resU-(fs%sRHOX**2*fs%U-fs%sRHOXold**2*fs%Uold)
            resV=time%dt*resV-(fs%sRHOY**2*fs%V-fs%sRHOYold**2*fs%Vold)
            resW=time%dt*resW-(fs%sRHOZ**2*fs%W-fs%sRHOZold**2*fs%Wold)
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals with under-relaxation
            fs%U=fs%U+resU*time%relax
            fs%V=fs%V+resV*time%relax
            fs%W=fs%W+resW*time%relax
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ PRESSURE SOLVER ======================
            ! Compute Umid
            call fs%get_Umid()
            ! Get rhoU/rhoV/rhoW
            call fs%rho_multiply()
            ! Solve Poisson equation
            call fs%update_laplacian(dt=time%dt,c2=C2)
            call fs%get_div(dt=time%dt)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            ! Correct pressure
            fs%P=fs%P+fs%psolv%sol
            ! Correct density
            fs%RHO=fs%RHO+fs%psolv%sol/C2; call fs%update_faceRHO()
            ! Correct mass flux
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%rhoU=fs%rhoU-time%dt*resU*((1.0_WP-fs%theta)*fs%sRHOXold**2+fs%theta*fs%sRHOX**2)/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%rhoV=fs%rhoV-time%dt*resV*((1.0_WP-fs%theta)*fs%sRHOYold**2+fs%theta*fs%sRHOY**2)/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%rhoW=fs%rhoW-time%dt*resW*((1.0_WP-fs%theta)*fs%sRHOZold**2+fs%theta*fs%sRHOZ**2)/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            ! Recover Umid
            call fs%rho_divide()
            ! Recover U
            call fs%get_U()
            ! Also update internal energy
            call fs%get_pdil(P=fs%psolv%sol,Pdil=resU); sc%E=sc%E+time%dt*resU/fs%RHO
            ! ===================================================
            
            ! ============= CVG CHECKING ========================
            residual_check: block
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
               use parallel, only: MPI_REAL_WP
               integer :: ierr
               RHOcvg=maxval(abs(resE-fs%RHO)/fs%RHO); call MPI_ALLREDUCE(MPI_IN_PLACE,RHOcvg,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
            end block residual_check
            ! ===================================================
            
            ! Increment sub-iteration
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
         
         ! Compute Mach number
         call get_c2(); Ma=sqrt((Ui**2+Vi**2+Wi**2)/C2)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max(fs%RHO)
         call analyze_conservation()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(resE,resU,resV,resW,Ui,Vi,Wi,Ma)
      if (use_sgsb.or.use_sgss) deallocate(gradU)
   end subroutine simulation_final
   
   
end module simulation
