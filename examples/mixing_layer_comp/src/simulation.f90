!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use hypre_str_class,   only: hypre_str
   use compress_class,    only: compress
   use energy_class,      only: energy
   use accelerator_class, only: accelerator
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
   type(accelerator), public :: accel
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   real(WP) :: RHOcvg=0.0_WP,RHOtol=1.0e-3_WP
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resE
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,C2
   
   !> Equation of state
   real(WP) :: Gamma,Mach
   
   !> Initial perturbations
   integer :: nwaveX,nwaveZ
   real(WP) :: wamp
   real(WP), dimension(:), allocatable :: wnumbX,wshiftX,wnumbZ,wshiftZ
   
contains
   
   
   !> Obtain density from equation of state
   subroutine get_rho()
      implicit none
      integer :: i,j,k
      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               !fs%RHO(i,j,k)=fs%P(i,j,k)/(sc%E(i,j,k)*(Gamma-1.0_WP))
               fs%RHO(i,j,k)=1.0_WP/(sc%E(i,j,k)*Gamma*Mach**2*(Gamma-1.0_WP))
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
               !C2(i,j,k)=Gamma*fs%P(i,j,k)/fs%RHO(i,j,k)
               C2(i,j,k)=huge(1.0_WP)
            end do
         end do
      end do
   end subroutine get_c2

   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read in EOS parameters
      initialize_eos: block
         ! EOS parameters
         call param_read('Gamma',Gamma)
         call param_read('Mach number',Mach)
      end block initialize_eos
      
      ! Initialize disturbances
      initialize_disturbances: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         use parallel,  only: MPI_REAL_WP
         use mpi_f08,   only: MPI_BCAST
         use string,    only: str_long
         use messager,  only: log
         character(str_long) :: message
         integer :: ierr,n
         ! Prepare interface disturbance in X
         call param_read('Wave amplitude',wamp)
         call param_read('NwaveX',NwaveX,default=-1)
         if (NwaveX.gt.0) then
            allocate(wnumbX(nwaveX),wshiftX(nwaveX))
            call param_read('wnumbX',wnumbX)
            call param_read('wshiftX',wshiftX)
         else
            nwaveX=6
            allocate(wnumbX(nwaveX),wshiftX(nwaveX))
            wnumbX=[3.0_WP,4.0_WP,5.0_WP,6.0_WP,7.0_WP,8.0_WP]*twoPi/cfg%xL
            if (cfg%amRoot) then
               do n=1,nwaveX
                  wshiftX(n)=random_uniform(lo=-0.5_WP*cfg%xL,hi=+0.5_WP*cfg%xL)
               end do
            end if
            call MPI_BCAST(wshiftX,nwaveX,MPI_REAL_WP,0,cfg%comm,ierr)
         end if
         ! Prepare interface disturbance in Z
         call param_read('NwaveZ',NwaveZ,default=-1)
         if (NwaveZ.gt.0) then
            allocate(wnumbZ(nwaveZ),wshiftZ(nwaveZ))
            call param_read('wnumbZ',wnumbZ)
            call param_read('wshiftZ',wshiftZ)
         else
            nwaveZ=6
            allocate(wnumbZ(nwaveZ),wshiftZ(nwaveZ))
            wnumbZ=[3.0_WP,4.0_WP,5.0_WP,6.0_WP,7.0_WP,8.0_WP]*twoPi/cfg%zL
            if (cfg%amRoot) then
               do n=1,nwaveZ
                  wshiftZ(n)=random_uniform(lo=-0.5_WP*cfg%zL,hi=+0.5_WP*cfg%zL)
               end do
            end if
            call MPI_BCAST(wshiftZ,nwaveZ,MPI_REAL_WP,0,cfg%comm,ierr)
         end if
         ! Handle 2D case
         if (cfg%nz.eq.1) wnumbZ=0.0_WP
         ! Print out initial disturbance
         if (cfg%amRoot) then
            write(message,'("[Initial conditions] =>  NwaveX =",i6)')              nwaveX; call log(message)
            write(message,'("[Initial conditions] =>  wnumbX =",1000(es12.5,x))')  wnumbX; call log(message)
            write(message,'("[Initial conditions] => wshiftX =",1000(es12.5,x))') wshiftX; call log(message)
            write(message,'("[Initial conditions] =>  NwaveZ =",i6)')              nwaveZ; call log(message)
            write(message,'("[Initial conditions] =>  wnumbZ =",1000(es12.5,x))')  wnumbZ; call log(message)
            write(message,'("[Initial conditions] => wshiftZ =",1000(es12.5,x))') wshiftZ; call log(message)
         end if
      end block initialize_disturbances
      
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
      end block allocate_work_arrays
      
      ! Initialize convergence accelerator
      initialize_accelerator: block
         call accel%initialize(cfg=cfg,storage_size=6,name='Density solver')
         accel%beta=0.75_WP ! Under-relaxation improves convergence
      end block initialize_accelerator

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
         use compress_class,  only: slip
         real(WP) :: Re
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Compressible NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         ! Assign constant viscosity based on Reynolds number
         call param_read('Reynolds number',Re); fs%visc=1.0_WP/Re
         ! Define boundary conditions
         call fs%add_bcond(name='yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call fs%add_bcond(name='ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
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
         use energy_class, only: dirichlet,neumann
         real(WP) :: Pr
         ! Create energy solver
         call sc%initialize(cfg=cfg,name='Energy')
         ! Define boundary conditions
         call sc%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call sc%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         ! Assign constant diffusivity
         sc%Cv=1.0_WP/(Gamma*(Gamma-1.0_WP)*Mach**2)
         ! Assign constant diffusivity
         call param_read('Prandtl number',Pr); sc%diff=sc%Cv*Gamma*fs%visc/Pr
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Energy',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_energy
      
      ! Initialize our initial conditions
      initial_conditions: block
         integer :: i,j,k,nX,nZ
         real(WP) :: y,delta,rho1,rho2
         ! Set pressure from Mach number
         fs%P=1.0_WP/(Gamma*Mach**2)
         ! Read in density profile
         call param_read('Density thickness',delta,default=1.0_WP)
         call param_read('Density top',rho1,default=1.0_WP)
         call param_read('Density bottom',rho2,default=1.0_WP)
         ! Initialize density field using tanh and corresponding energy from EOS
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  fs%RHO(i,j,k)=0.5_WP*(rho1+rho2)+0.5_WP*(rho2-rho1)*tanh(y/delta)
                  sc%E(i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
               end do
            end do
         end do
         ! Initialize velocity field
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Distance to the nominal interface position
                  y=cfg%ym(j)
                  ! Perturb interface position
                  do nX=1,nwaveX
                     do nZ=1,nwaveZ
                        y=y+wamp*cos(wnumbX(nX)*(cfg%x(i)-wshiftX(nX)))*cos(wnumbZ(nZ)*(cfg%zm(k)-wshiftZ(nZ)))
                     end do
                  end do
                  ! Hyperbolic tangent
                  fs%U(i,j,k)=0.5_WP*tanh(y)
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
      end block initial_conditions
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='mixing')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('energy',sc%E)
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
      
      !> Function that localizes y- boundary
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      !> Function that localizes y- boundary
      function ym_locator_sc(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin-1) isIn=.true.
      end function ym_locator_sc
      
      !> Function that localizes y+ boundary
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
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
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Perform sub-iterations until RHO is sufficiently converged
         RHOcvg=huge(1.0_WP); time%it=0; call accel%restart(ig=fs%RHO)
         do while (RHOcvg.gt.RHOtol.and.time%it.lt.time%itmax.or.time%it.lt.time%itmin)
            
            ! ============= ENERGY SOLVER =======================
            ! Explicit calculation of drhoE/dt from energy equation
            call sc%get_drhoEdt(resE,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Add pressure dilatation term, storing it
            call fs%get_pdil(resU); resE=resE+resU
            
            ! Assemble explicit residual
            resE=time%dt*resE-(fs%RHO*sc%E-fs%RHOold*sc%Eold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resE,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Apply this residual
            sc%E=sc%E+resE
            
            ! Remember pressure dilatation in resE
            resE=resU
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            ! ============ UPDATE DENSITY AND C2 ================
            update_rho: block
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
               use parallel, only: MPI_REAL_WP
               integer :: ierr
               ! Remember RHO
               resU=fs%RHO
               ! Calculate RHO using our convergence accelerator
               call get_rho(); call accel%increment(fs%RHO); call fs%update_sRHO()
               ! Check RHO convergence
               RHOcvg=maxval(abs(resU-fs%RHO)/fs%RHO); call MPI_ALLREDUCE(MPI_IN_PLACE,RHOcvg,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
               ! Compute speed of sound squared
               call get_c2()
            end block update_rho
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
            
            ! ============ PRESSURE SOLVER ======================
            ! Compute Umid
            call fs%get_Umid()
            ! Get rhoU/rhoV/rhoW
            call fs%rho_multiply()
            ! Correct MFR
            call fs%correct_mfr(dt=time%dt)
            ! Solve Poisson equation
            call fs%update_laplacian(dt=time%dt,c2=C2)
            call fs%get_div(dt=time%dt)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dt*resU*((1.0_WP-fs%theta)*fs%sRHOXold**2+fs%theta*fs%sRHOX**2)/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%rhoV=fs%rhoV-time%dt*resV*((1.0_WP-fs%theta)*fs%sRHOYold**2+fs%theta*fs%sRHOY**2)/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%rhoW=fs%rhoW-time%dt*resW*((1.0_WP-fs%theta)*fs%sRHOZold**2+fs%theta*fs%sRHOZ**2)/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            ! Recover Umid
            call fs%rho_divide()
            ! Recover U
            call fs%get_U()
            ! Also update internal energy (resE contains the predicted pressure dilatation term)
            sc%E=sc%E-time%dt*resE/fs%RHO
            call fs%get_pdil(resE)
            sc%E=sc%E+time%dt*resE/fs%RHO
            ! ===================================================
            
            ! Increment sub-iteration
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
         
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
      deallocate(resE,resU,resV,resW,Ui,Vi,Wi)
   end subroutine simulation_final
   
   
end module simulation
