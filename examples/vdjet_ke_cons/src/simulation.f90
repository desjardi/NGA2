!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use hypre_str_class,   only: hypre_str
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single low Mach flow solver and scalar solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs,ss
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
   
   !> Equation of state
   real(WP) :: Zst,rho0,rho1,rhost
   
   !> Inlet parameters
   real(WP) :: Zjet,Zcof
   real(WP) :: Djet,Dcof
   real(WP) :: Ujet,Ucof
   real(WP) :: thick
   
contains
   
   
   !> Obtain density from equation of state based on Burke-Schumann
   subroutine get_rho()
      implicit none
      integer :: i,j,k
      real(WP) :: Zclip
      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               ! Clip mixture fraction between - and 1
               Zclip=min(max(sc%SC(i,j,k),0.0_WP),1.0_WP)
               ! Pure mixing flamelet
               fs%RHO(i,j,k)=rho0*rho1/(rho1*(1.0_WP-Zclip)+rho0*Zclip)
               ! Burke-Schumann flamelet
               !if (Zclip.le.Zst) then
               !   fs%RHO(i,j,k)=Zst*rho0*rhost/(rhost*(Zst-Zclip)+rho0*Zclip)
               !else
               !   fs%RHO(i,j,k)=(1.0_WP-Zst)*rhost*rho1/(rho1*(1.0_WP-Zclip)+rhost*(Zclip-Zst))
               !end if
            end do
         end do
      end do
   end subroutine get_rho
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read in EOS parameters
      initialize_eos: block
         call param_read('rho0',rho0)
         call param_read('rho1',rho1)
         call param_read('rhost',rhost)
         call param_read('Zst',Zst)
         if (Zst.le.0.0_WP) Zst=0.0_WP+epsilon(Zst)
         if (Zst.ge.1.0_WP) Zst=1.0_WP-epsilon(Zst)
      end block initialize_eos
      
      ! Read in inlet parameters
      initialize_inlet: block
         call param_read('Z jet',Zjet)
         call param_read('D jet',Djet)
         call param_read('U jet',Ujet)
         call param_read('Z coflow',Zcof)
         call param_read('D coflow',Dcof)
         call param_read('U coflow',Ucof)
         call param_read('Thickness',thick,default=cfg%xL/real(cfg%nx,WP))
      end block initialize_inlet
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker
      
      ! Create a scalar solver
      create_scalar: block
         use vdscalar_class, only: dirichlet,neumann
         real(WP) :: diffusivity
         ! Create scalar solver
         call sc%initialize(cfg=cfg,name='MixFrac')
         ! Add slight backward bias to CN scheme
         sc%theta=sc%theta+1.0e-2_WP
         ! Define boundary conditions
         call sc%add_bcond(name='inflow' ,type=dirichlet,locator=xm_locator_sc)
         call sc%add_bcond(name='outflow',type=neumann  ,locator=xp_locator,dir='+x')
         ! Assign constant diffusivity
         call param_read('Dynamic diffusivity',diffusivity)
         sc%diff=diffusivity
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=7)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_scalar
      
      ! Create a low-Mach flow solver with bconds
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use lowmach_class,   only: dirichlet,clipped_neumann,slip
         real(WP) :: visc
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Variable density low Mach NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         !call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         !call fs%add_bcond(name='yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         !call fs%add_bcond(name='ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         !call fs%add_bcond(name='zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         !call fs%add_bcond(name='zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
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
      
      ! Initialize our mixture fraction field
      initialize_scalar: block
         use vdscalar_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc
         real(WP) :: radius
         ! Zero initial field
         sc%SC=0.0_WP
         ! Apply BCs
         call sc%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            radius=norm2([sc%cfg%ym(j),sc%cfg%zm(k)]-[0.0_WP,0.0_WP])
            sc%SC(i,j,k)=(Zcof+(Zjet-Zcof)*0.5_WP*(1.0_WP+tanh((0.5_WP*Djet-radius)/thick)))*0.5_WP*(1.0_WP+tanh((0.5_WP*Dcof-radius)/thick))
         end do
      end block initialize_scalar
      
      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc
         real(WP) :: radius
         ! Initialize density from scalar
         call get_rho(); call fs%update_sRHO()
         fs%RHOold=fs%RHO; fs%sRHOXold=fs%sRHOX; fs%sRHOYold=fs%sRHOY; fs%sRHOZold=fs%sRHOZ
         ! Initialize velocity field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply inlet condition
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            radius=norm2([fs%cfg%ym(j),fs%cfg%zm(k)]-[0.0_WP,0.0_WP])
            fs%U(i,j,k)=(Ucof+(Ujet-Ucof)*0.5_WP*(1.0_WP+tanh((0.5_WP*Djet-radius)/thick)))*0.5_WP*(1.0_WP+tanh((0.5_WP*Dcof-radius)/thick))
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Get Umid/Vmid/Wmid
         call fs%get_Umid()
         ! Calculate cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(dt=time%dt)
      end block initialize_velocity
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='vdjet')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('mixfrac',sc%SC)
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
         call mfile%add_column(sc%SCmax,'SCmax')
         call mfile%add_column(sc%SCmin,'SCmin')
         call mfile%add_column(sc%rhoSCint,'rhoSCint')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(fs%RHOint,'RHOint')
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
      
      !> Function that localizes y+ boundary
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      !> Function that localizes z- boundary
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      !> Function that localizes z+ boundary
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      !> Function that localizes the x+ boundary
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
      !> Function that localizes the x- boundary
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator
      
      !> Function that localizes the x- boundary for SC
      function xm_locator_sc(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin-1) isIn=.true.
      end function xm_locator_sc
      
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
         
         ! Remember old scalar, velocity, and density
         sc%SCold=sc%SC
         fs%RHOold=fs%RHO
         fs%sRHOXold=fs%sRHOX
         fs%sRHOYold=fs%sRHOY
         fs%sRHOZold=fs%sRHOZ
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ============= SCALAR SOLVER =======================
            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Assemble explicit residual
            resSC=time%dt*resSC-(fs%RHO*sc%SC-fs%RHOold*sc%SCold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Apply this residual
            sc%SC=sc%SC+resSC
            
            ! Apply other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            dirichlet_scalar: block
               use vdscalar_class, only: bcond
               integer :: n,i,j,k
               type(bcond), pointer :: mybc
               real(WP) :: radius
               call sc%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  radius=norm2([sc%cfg%ym(j),sc%cfg%zm(k)]-[0.0_WP,0.0_WP])
                  sc%SC(i,j,k)=(Zcof+(Zjet-Zcof)*0.5_WP*(1.0_WP+tanh((0.5_WP*Djet-radius)/thick)))*0.5_WP*(1.0_WP+tanh((0.5_WP*Dcof-radius)/thick))
               end do
            end block dirichlet_scalar
            ! ===================================================
            
            ! ============ UPDATE DENSITY========================
            call get_rho(); call fs%update_sRHO()
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
            dirichlet_velocity: block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               real(WP) :: radius
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  radius=norm2([fs%cfg%ym(j),fs%cfg%zm(k)]-[0.0_WP,0.0_WP])
                  fs%U(i,j,k)=(Ucof+(Ujet-Ucof)*0.5_WP*(1.0_WP+tanh((0.5_WP*Djet-radius)/thick)))*0.5_WP*(1.0_WP+tanh((0.5_WP*Dcof-radius)/thick))
               end do
            end block dirichlet_velocity
            ! ===================================================
            
            ! ============ PRESSURE SOLVER ======================
            ! Compute Umid
            call fs%get_Umid()
            ! Get rhoU/rhoV/rhoW
            call fs%rho_multiply()
            ! Correct MFR
            call fs%correct_mfr(dt=time%dt)
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%get_div(dt=time%dt)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dt*resU*((1.0_WP-fs%theta)*fs%sRHOXold**2+fs%theta*fs%sRHOX**2)/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%rhoV=fs%rhoV-time%dt*resV*((1.0_WP-fs%theta)*fs%sRHOYold**2+fs%theta*fs%sRHOY**2)/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%rhoW=fs%rhoW-time%dt*resW*((1.0_WP-fs%theta)*fs%sRHOZold**2+fs%theta*fs%sRHOZ**2)/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            ! Recover Umid
            call fs%rho_divide()
            ! Recover U
            call fs%get_U()
            ! ===================================================
            
            ! Increment sub-iteration
            time%it=time%it+1
            
         end do
         
         ! Could consider doing final scalar step here
         !call sc%get_drhoSCdt(resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
         !resSC=time%dt*resSC-(fs%RHO*sc%SC-fs%RHOold*sc%SCold)
         !call sc%solve_implicit(time%dt,resSC,fs%RHO,fs%RHOold,fs%rhoU,fs%rhoV,fs%rhoW)
         !sc%SC=sc%SC+resSC
         !call sc%apply_bcond(time%t,time%dt)
         
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
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi)
   end subroutine simulation_final
   
   
end module simulation
