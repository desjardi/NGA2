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
   logical        :: use_sgs
   type(sgsmodel) :: sgs
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   real(WP) :: RHOcvg=0.0_WP,RHOtol=1.0e-6_WP
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resE
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,C2,Ma
   
   !> Equation of state and shock properties
   real(WP) :: Gamma
   real(WP) :: M2,p2,rho2,a2,u2
   real(WP) :: M1,p1,rho1,a1,u1

   !> Fluid properties
   real(WP) :: Reynolds=1.0_WP,Prandtl,bulk2shear
   
   !> Postprocessing of dissipation
   type(monitor) :: dissfile
   real(WP), dimension(:,:,:), allocatable :: sdiss,ddiss
   real(WP) :: sdiss_int,ddiss_int
   
contains
   
   
   !> Function that returns a smooth Heaviside representation of exact shock
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock
   
   
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
      ! Constant kinematic viscosity
      fs%viscs=fs%RHO
      ! Constant dynamic viscosity
      !fs%visc=1.0_WP
   end subroutine get_visc
   
   
   !> Calculate dissipations
   subroutine get_dissipation()
      implicit none
      ! Get full dilatational and solenoidal dissipations
      call fs%get_dilatational_dissipation(diss=ddiss); call cfg%integrate(ddiss,integral=ddiss_int); ddiss_int=ddiss_int/(cfg%yL*cfg%zL)
      call fs%get_solenoidal_dissipation  (diss=sdiss); call cfg%integrate(sdiss,integral=sdiss_int); sdiss_int=sdiss_int/(cfg%yL*cfg%zL)
   end subroutine get_dissipation
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Prepare shock data
      initialize_shock: block
         use mathtools, only: Pi
         use string,    only: str_long
         use messager,  only: log
         character(str_long) :: message
         ! Read in gamma
         call param_read('Gamma',gamma)
         ! Read in minimal shock information
         call param_read('Mach',M1)
         ! Generate preshock conditions
         rho1=1.0_WP
         p1=0.25_WP*rho1/gamma*((gamma+1.0_WP)*M1/(M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
         a1=sqrt(gamma*p1/rho1)
         ! Generate pre-shock velocity
         u1=M1*a1
         ! Also generate post-shock conditions
         p2=p1*(2.0_WP*gamma/(gamma+1.0_WP)*(M1**2-1.0_WP)+1.0_WP)
         rho2=rho1*(gamma+1.0_WP)*M1**2/((gamma-1.0_WP)*M1**2+2.0_WP)
         u2=u1*rho1/rho2
         a2=sqrt(gamma*p2/rho2)
         M2=u2/a2
         ! Output shock info
         if (cfg%amRoot) then
            write(message,'("[Pre -shock conditions] =>  rho1=",es12.5)')  rho1; call log(message)
            write(message,'("[Pre -shock conditions] =>    p1=",es12.5)')    p1; call log(message)
            write(message,'("[Pre -shock conditions] =>    u1=",es12.5)')    u1; call log(message)
            write(message,'("[Pre -shock conditions] =>    a1=",es12.5)')    a1; call log(message)
            write(message,'("[Pre -shock conditions] =>    M1=",es12.5)')    M1; call log(message)
            write(message,'("[Post-shock conditions] =>  rho2=",es12.5)')  rho2; call log(message)
            write(message,'("[Post-shock conditions] =>    p2=",es12.5)')    p2; call log(message)
            write(message,'("[Post-shock conditions] =>    u2=",es12.5)')    u2; call log(message)
            write(message,'("[Post-shock conditions] =>    a2=",es12.5)')    a2; call log(message)
            write(message,'("[Post-shock conditions] =>    M2=",es12.5)')    M2; call log(message)
         end if
      end block initialize_shock
      
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
         allocate(sdiss(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(ddiss(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
      
      ! Create a compressible flow solver with bconds
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use compress_class,  only: dirichlet,clipped_neumann
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Compressible NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         ! Assign viscosity for Re_1=1
         call param_read('Bulk to shear ratio',bulk2shear)
         ! Define boundary conditions
         !call fs%add_bcond(name='xm',type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         !call fs%add_bcond(name='xp',type=dirichlet      ,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         !call fs%add_bcond(name='xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
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
      
      ! Create a energy solver
      create_energy: block
         use energy_class, only: neumann,dirichlet
         real(WP) :: Pr
         ! Create energy solver
         call sc%initialize(cfg=cfg,name='Energy')
         ! Set Cv such that T1=1
         sc%Cv=p1/(rho1*(Gamma-1.0_WP))
         call param_read('Prandtl',Prandtl)
         ! Define boundary conditions
         !call sc%add_bcond(name='xm',type=dirichlet,locator=xm_locator_sc,dir='-x')
         !call sc%add_bcond(name='xp',type=dirichlet,locator=xp_locator   ,dir='+x')
         !call sc%add_bcond(name='xp',type=neumann  ,locator=xp_locator   ,dir='+x')
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Energy',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_energy
      
      ! Initialize our initial conditions
      initial_conditions: block
         integer :: i,j,k
         ! Initialize all fields
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Setup normal shock at t=0
                  fs%RHO(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i)+0.40_WP,delta=1.13_WP)  ! Shift rho a bit
                  fs%U  (i,j,k)=u1  +(u2  -u1  )*Hshock(cfg%x (i)+2.35_WP,delta=1.10_WP)  ! Shift u a bit
                  fs%P  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i)+0.02_WP,delta=1.00_WP)
                  ! Corresponding internal energy
                  sc%E(i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
               end do
            end do
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Get face RHO
         fs%RHOold=fs%RHO
         call fs%prepare_weno(fs%RHO)
         call sc%prepare_weno(sc%E)
         call fs%update_faceRHO()
         fs%sRHOXold=fs%sRHOX
         fs%sRHOYold=fs%sRHOY
         fs%sRHOZold=fs%sRHOZ
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         fs%Pold=fs%P
         call fs%update_faceRHO()
         ! Get Umid and mass flux
         call fs%get_Umid()
         call fs%rho_multiply()
         ! Calculate cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
         ! Compute local Mach number
         call get_c2(); Ma=sqrt((Ui**2+Vi**2+Wi**2)/C2)
      end block initial_conditions
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',use_sgs)
         if (use_sgs) then
            allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         end if
      end block create_sgs
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='staticshock')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('energy',sc%E)
         call ens_out%add_scalar('sdiss',sdiss)
         call ens_out%add_scalar('ddiss',ddiss)
         if (use_sgs) call ens_out%add_scalar('viscb',sgs%visc)
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
         call get_dissipation()
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
         call mfile%add_column(sc%rhoEint,'rhoEint')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(fs%RHOint,'RHOint')
         call mfile%add_column(   RHOcvg,'RHOcvg')
         call mfile%add_column(  time%it,'Subiter')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create dissipation monitor
         dissfile=monitor(fs%cfg%amRoot,'dissipation')
         call dissfile%add_column(time%n,'Timestep number')
         call dissfile%add_column(time%t,'Time')
         call dissfile%add_column(time%dt,'Timestep size')
         call dissfile%add_column(sdiss_int,'Solenoidal dissipation')
         call dissfile%add_column(ddiss_int,'Dilatational dissipation')
         call dissfile%write()
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
   
   contains
      
      !> Function that localizes x- boundary
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator
      
      !> Function that localizes x- boundary for a scalar
      function xm_locator_sc(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin-1) isIn=.true.
      end function xm_locator_sc
      
      !> Function that localizes x+ boundary
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
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
         fs%Pold=fs%P
         
         ! Apply time-varying Dirichlet conditions
         
         ! ============= RECOMPUTE VISCOSITY =================
         call get_visc(); fs%viscb=bulk2shear*fs%viscs
         sc%diff=sc%Cv*Gamma*fs%viscs/Prandtl
         ! ===================================================
         
         ! ============= SUBGRID SCALE MODELING ==============
         if (use_sgs) then
            sgs_modeling: block
               use sgsmodel_class, only: vreman,localartif
               call fs%get_gradU(gradU)
               call sgs%get_visc(type=localartif,dt=time%dt,rho=fs%RHO,gradu=gradU)
               fs%viscb=sgs%visc
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
         !call fs%predict_rho(dt=time%dt)
         ! ===================================================
         
         ! Perform sub-iterations until RHO is sufficiently converged
         RHOcvg=huge(1.0_WP); time%it=0
         do while (RHOcvg.gt.RHOtol.and.time%it.lt.time%itmax.or.time%it.lt.time%itmin)
            
            ! ============= ENERGY SOLVER =======================
            ! Explicit calculation of drhoE/dt from energy equation
            call sc%get_drhoEdt(resE,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Add pressure dilatation term
            call fs%get_pdil(P=fs%P,Pdil=resU); resE=resE+resU
            
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
            
            ! Compute speed of sound squared at mid time
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
         call get_dissipation()
         call mfile%write()
         call dissfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(resE,resU,resV,resW,Ui,Vi,Wi,Ma,ddiss,sdiss)
      if (use_sgs) deallocate(gradU)
   end subroutine simulation_final
   
   
end module simulation
