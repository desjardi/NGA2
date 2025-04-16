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
   use iterator_class,    only: iterator
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
   real(WP) :: RHOcvg=0.0_WP,RHOtol=1.0e-4_WP
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resE
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,C2,Ma
   
   !> Equation of state and shock properties
   real(WP) :: Gamma
   real(WP) :: Ms,Xs0,theta,delta,sintheta,costheta
   real(WP) :: M2,p2,rho2,a2,u2,u2x,u2y,p1,rho1,a1,us
   
   !> Sponge zone
   real(WP) :: yspg
   type(iterator) :: spg
   
contains
   
   
   !> Function that returns a smooth Heaviside representation of exact shock
   function Hshock(x,y,t) result(H)
      real(WP), intent(in)  :: x,y,t
      real(WP) :: H
      H=0.5_WP+0.5_WP*tanh((sintheta*y-costheta*(x-Xs0-u2x*t))/delta)
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
         call param_read('Shock Mach',Ms)
         call param_read('Shock angle',theta); theta=theta*Pi/180.0_WP; costheta=cos(theta); sintheta=sin(theta)
         call param_read('Shock position',Xs0)
         call param_read('Shock thickness',delta)
         ! Generate preshock conditions
         rho1=1.0_WP                                                   ! This is our reference density
         p1=0.25_WP*rho1/gamma*((gamma+1.0_WP)*Ms/(Ms**2-1.0_WP))**2   ! We choose p1 to ensure that u2=1
         a1=sqrt(gamma*p1/rho1)
         ! Generate shock velocity
         us=Ms*a1
         ! Also generate post-shock conditions
         p2=p1*(2.0_WP*gamma*Ms**2-(gamma-1.0_WP))/(gamma+1.0_WP)
         rho2=rho1*(gamma+1.0_WP)*Ms**2/((gamma-1.0_WP)*Ms**2+2.0_WP)
         u2=us*(1.0_WP-rho1/rho2)   ! Should come out to 1
         a2=sqrt(gamma*p2/rho2)
         M2=u2/a2
         u2x=+u2*costheta
         u2y=-u2*sintheta
         ! Output shock info
         if (cfg%amRoot) then
            write(message,'("[     Shock conditions] =>    Ms=",es12.5)')    Ms; call log(message)
            write(message,'("[     Shock conditions] =>    Vs=",es12.5)')    us; call log(message)
            write(message,'("[     Shock conditions] => theta=",es12.5)') theta; call log(message)
            write(message,'("[Pre -shock conditions] =>  rho1=",es12.5)')  rho1; call log(message)
            write(message,'("[Pre -shock conditions] =>    p1=",es12.5)')    p1; call log(message)
            write(message,'("[Pre -shock conditions] =>    a1=",es12.5)')    a1; call log(message)
            write(message,'("[Post-shock conditions] =>  rho2=",es12.5)')  rho2; call log(message)
            write(message,'("[Post-shock conditions] =>    p2=",es12.5)')    p2; call log(message)
            write(message,'("[Post-shock conditions] =>    a2=",es12.5)')    a2; call log(message)
            write(message,'("[Post-shock conditions] =>    u2=",es12.5)')    u2; call log(message)
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
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax)
         time%itmin=2
      end block initialize_timetracker
      
      ! Create a compressible flow solver with bconds
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Compressible NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         fs%viscb=0.01_WP
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
         use energy_class, only: neumann
         real(WP) :: Pr
         ! Create energy solver
         call sc%initialize(cfg=cfg,name='Energy')
         sc%diff=0.01_WP
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Energy',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_energy
      
      ! Initialize our initial conditions
      initial_conditions: block
         integer  :: i,j,k
         ! Initialize all fields
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Setup normal shock at t=0
                  fs%RHO(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i),cfg%ym(j),0.0_WP)
                  fs%P  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i),cfg%ym(j),0.0_WP)
                  fs%U  (i,j,k)=        u2x     *Hshock(cfg%x (i),cfg%ym(j),0.0_WP)
                  fs%V  (i,j,k)=        u2y     *Hshock(cfg%xm(i),cfg%y (j),0.0_WP)
                  ! Corresponding internal energy
                  sc%E(i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
               end do
            end do
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Get Umid/Vmid/Wmid
         call fs%update_sRHO(); fs%RHOold=fs%RHO; fs%sRHOXold=fs%sRHOX; fs%sRHOYold=fs%sRHOY; fs%sRHOZold=fs%sRHOZ
         call fs%get_Umid()
         ! Calculate cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
         ! Compute local Mach number
         call get_c2(); Ma=sqrt((Ui**2+Vi**2+Wi**2)/C2)
      end block initial_conditions
      
      ! Create sponge zone
      create_sponge: block
         call param_read('Y sponge',yspg)
         spg=iterator(cfg,'Sponge zone',sponge_locator)
      end block create_sponge

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
         ens_out=ensight(cfg=cfg,name='shockreflection')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('energy',sc%E)
         call ens_out%add_scalar('viscb',fs%viscb)
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
      
      !> Function that localizes the y+ sponge layer
      function sponge_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (pg%ym(j).ge.yspg) isIn=.true.
      end function sponge_locator
      
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
         !top_update: block
         !   integer :: i,j,k
         !   if (cfg%jproc.eq.cfg%npy) then
         !      do k=cfg%kmino_,cfg%kmaxo_
         !         do j=cfg%jmax+1,cfg%jmaxo
         !            do i=cfg%imino_,cfg%imaxo_
         !               ! Setup normal shock at current time
         !               fs%RHO(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i),cfg%ym(j),time%t)
         !               fs%P  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i),cfg%ym(j),time%t)
         !               fs%U  (i,j,k)=        u2x     *Hshock(cfg%x (i),cfg%ym(j),time%t)
         !               fs%V  (i,j,k)=        u2y     *Hshock(cfg%xm(i),cfg%y (j),time%t)
         !               ! Corresponding internal energy
         !               sc%E(i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
         !            end do
         !         end do
         !      end do
         !   end if
         !   if (cfg%jproc.eq.1) then
         !      do k=cfg%kmino_,cfg%kmaxo_
         !         do j=cfg%jmino,cfg%jmin-1
         !            do i=cfg%imino_,cfg%imaxo_
         !               ! Setup normal shock at current time
         !               fs%RHO(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i),cfg%ym(j),time%t)
         !               fs%P  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i),cfg%ym(j),time%t)
         !               fs%U  (i,j,k)=        u2x     *Hshock(cfg%x (i),cfg%ym(j),time%t)
         !               fs%V  (i,j,k)=        u2y     *Hshock(cfg%xm(i),cfg%y (j),time%t)
         !               ! Corresponding internal energy
         !               sc%E(i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
         !            end do
         !         end do
         !      end do
         !   end if
         !end block top_update
         
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
            sc%E=sc%E+resE
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ UPDATE DENSITY AND C2 ================
            ! Remember RHO
            resE=fs%RHO
            
            ! Calculate RHO predictor
            call get_rho(); call fs%update_sRHO()
            
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
            
            ! Apply these residuals
            fs%U=fs%U+resU
            fs%V=fs%V+resV
            fs%W=fs%W+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! Apply time-varying sponge conditions
            !sponge_update: block
            !   integer :: i,j,k,n
            !   do n=1,spg%no_
            !      i=spg%map(1,n)
            !      j=spg%map(2,n)
            !      k=spg%map(3,n)
            !      ! Apply exact shock solution
            !      fs%RHO(i,j,k)=rho1+(rho2-rho1)*Hshock(cfg%xm(i),cfg%ym(j),time%t)
            !      fs%P  (i,j,k)=p1  +(p2  -p1  )*Hshock(cfg%xm(i),cfg%ym(j),time%t)
            !      fs%U  (i,j,k)=        u2x     *Hshock(cfg%x (i),cfg%ym(j),time%t)
            !      fs%V  (i,j,k)=        u2y     *Hshock(cfg%xm(i),cfg%y (j),time%t)
            !      ! Corresponding internal energy
            !      sc%E(i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
            !   end do
            !   call fs%update_sRHO()
            !   call get_c2()
            !end block sponge_update
            
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
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%RHO=fs%RHO+fs%psolv%sol/C2; call fs%update_sRHO()
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
         call mfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(resE,resU,resV,resW,Ui,Vi,Wi,Ma,gradU)
   end subroutine simulation_final
   
   
end module simulation
