!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use matm_class,        only: matm
   use lpt_class,         only: lpt
   use partmesh_class,    only: partmesh
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use string,            only: str_medium
   implicit none
   private

   !> Single two-phase flow solver, volume fraction solver, material model, and LPT solver set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(lpt),         public :: lp
   type(timetracker), public :: time

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,cvgfile
   type(monitor) :: lptfile

   public :: simulation_init,simulation_run,simulation_final

   !> Problem definition
   character(len=str_medium) :: phase
   integer :: relax_model

   !> Work arrays and fluid properties
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resE
   real(WP), dimension(:,:,:), allocatable :: dRHOdt

 contains

   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
        time=timetracker(amRoot=cfg%amRoot)
        call param_read('Max timestep size',time%dtmax)
        call param_read('Max cfl number',time%cflmax)
        call param_read('Max time',time%tmax)
        call param_read('Max steps',time%nmax)
        time%dt=time%dtmax
        time%itmax=2
      end block initialize_timetracker


      ! Allocate work arrays
      allocate_work_arrays: block
        allocate(dRHOdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(resE  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize our LPT
      initialize_lpt: block
        use random, only: random_uniform
        use mathtools,  only: Pi
        use parallel,   only: amRoot
        use lpt_class,  only: drag_model
        real(WP) :: dp,VFp,L,x0,vol,volp
         integer :: i,np
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',drag_model,default='Tenneti')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Get particle volume fraction from the input
         call param_read('Particle volume fraction',VFp)
         ! Get curtain width from the input
         call param_read('Particle curtain width',L)
         ! Get starting position
         x0=0.5_WP*(lp%cfg%x(lp%cfg%imax+1)-lp%cfg%x(lp%cfg%imin)-L)
         ! Get number of particles
         vol=L*lp%cfg%yL*lp%cfg%zL
         volp=Pi*dp**3/6.0_WP
         np=int(VFp * vol / volp)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Root process initializes np particles randomly
         if (lp%cfg%amRoot) then
            call lp%resize(np)
            do i=1,np
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! Assign random position in the domain
               lp%p(i)%pos=[random_uniform(x0,x0+L),&
                    &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
                    &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
               if (lp%cfg%nz.eq.1) lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin)
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               ! Give zero collision force
               lp%p(i)%col=0.0_WP
               ! Give it temperature
               lp%p(i)%T=293.15_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         lp%Tcol=5.0_WP*time%dt
         if (amRoot) then
            print*,"===== Particle Setup Description ====="
            print*,'Number of particles', np
            print*,'Mean volume fraction',vfp
         end if
       end block initialize_lpt


       ! Create partmesh object for Lagrangian particle output
       create_pmesh: block
         pmesh=partmesh(nvar=0,name='lpt')
         call lp%update_partmesh(pmesh)
       end block create_pmesh


       ! Initialize our VOF solver and field
       create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,elvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
         call param_read('Phase',phase)
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  select case(trim(phase))
                  case('gas')
                     vf%VF(i,j,k)=0.0_WP
                  case('liquid')
                     vf%VF(i,j,k)=1.0_WP
                  case default
                  end select
                  vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Boundary conditions on VF are built into the mast solver
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set initial interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
       end block create_and_initialize_vof


       ! Create a compressible two-phase flow solver
       create_and_initialize_flow_solver: block
         use mast_class, only: clipped_neumann,dirichlet,bc_scope,bcond,mech_egy_mech_hhz
         use matm_class, only : air
         use ils_class,  only: pcg_bbox,gmres_smg
         use mathtools,  only: Pi
         use parallel,   only: amRoot
         integer :: i,j,k,n
         real(WP), dimension(3) :: xyz
         real(WP) :: gamm,Pref,visc,hdff
         real(WP) :: xshock,vshock,relshockvel
         real(WP) :: rho0, P0, rho1, P1, Ma1, Ma, rho, P, U
         type(bcond), pointer :: mybc
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get EOS parameters from input
         call param_read('Gamma',gamm)
         ! Register equations of state
         select case(trim(phase))
         case('gas')
            call matmod%register_idealgas('gas',gamm)
         case('liquid')
            call param_read('Pref', Pref)
            call matmod%register_stiffenedgas('liquid',gamm,Pref)
         end select
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         ! Assign constant viscosity, also heat diffusion
         call param_read('Dynamic viscosity',visc)
         call param_read('Thermal conductivity',hdff)
         select case(trim(phase))
         case('gas')
            call matmod%register_diffusion_thermo_models(viscmodel_gas=air, viscconst_gas=visc, &
                 viscconst_liquid=0.0_WP, hdffconst_gas=hdff, hdffconst_liquid=0.0_WP)
            fs%Lrho = 0.0_WP
            fs%LrhoE= 0.0_WP
         case('liquid')
            call matmod%register_diffusion_thermo_models(viscconst_gas=0.0_WP, viscconst_liquid=visc, &
                 hdffconst_gas=0.0_WP, hdffconst_liquid=hdff)
            fs%Grho = 0.0_WP
         end select
         ! Read in surface tension coefficient
         fs%sigma=0.0_WP
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_bbox,implicit_ils=gmres_smg)
         ! Initially 0 velocity in y and z
         fs%Vi = 0.0_WP; fs%Wi = 0.0_WP
         ! Zero face velocities as well for the sake of dirichlet boundaries
         fs%V = 0.0_WP; fs%W = 0.0_WP
         ! Set up initial thermo properties (0)
         ! Default is from Meng & Colonius (2015),
         ! which is basically the same as Igra & Takayama and Terashima & Tryggvason
         call param_read('Pre-shock density',rho0,default=1.204_WP)
         call param_read('Pre-shock pressure',P0,default=1.01325e5_WP)
         call param_read('Mach number of shock',Ma,default=1.47_WP)
         call param_read('Shock location',xshock)
         ! Use shock relations to get post-shock numbers (1)
         P1 = P0 * (2.0_WP*gamm*Ma**2 - (gamm-1.0_WP)) / (gamm+1.0_WP)
         rho1 = rho0 * (Ma**2 * (gamm+1.0_WP) / ((gamm-1.0_WP)*Ma**2 + 2.0_WP))
         ! Calculate post-shock Mach number
         Ma1 = sqrt(((gamm-1.0_WP)*(Ma**2)+2.0_WP)/(2.0_WP*gamm*(Ma**2)-(gamm-1.0_WP)))
         ! Calculate post-shock velocity
         vshock = -Ma1 * sqrt(gamm*P1/rho1) + Ma*sqrt(gamm*P0/rho0)
         ! Velocity at which shock moves
         relshockvel = -rho1*vshock/(rho0-rho1)

         if (amRoot) then
            print*,"===== Fluid Setup Description ====="
            print*,'Mach number', Ma
            print*,'Pre-shock:  Density',rho0,'Pressure',P0
            print*,'Post-shock: Density',rho1,'Pressure',P1,'Velocity',vshock
            print*,'Shock velocity', relshockvel
         end if

         ! Initialize gas phase quantities
         do i=fs%cfg%imino_,fs%cfg%imaxo_
            ! pressure, velocity, use matmod for energy
            if (fs%cfg%x(i).lt.xshock) then
               rho = rho1
               U = vshock
               P = P1
            else
               rho = rho0
               U = 0.0_WP
               P = P0
            end if
            select case(trim(phase))
            case('gas')
               fs%Grho(i,:,:)  = rho
               fs%Ui(i,:,:)    = U
               fs%GrhoE(i,:,:) = matmod%EOS_energy(P,rho,U,0.0_WP,0.0_WP,'gas')
            case('liquid')
               fs%Lrho(i,:,:)  = rho
               fs%Ui(i,:,:)    = U
               fs%LrhoE(i,:,:) = matmod%EOS_energy(P,rho,U,0.0_WP,0.0_WP,'liquid')
            end select
         end do

         ! Define boundary conditions - initialized values are intended dirichlet values too, for the cell centers
         call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1)

         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         ! Apply face BC - inflow
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=vshock
         end do
         ! Apply face BC - outflow
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt,bc_scope)

         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         ! Perform initial pressure relax
         relax_model = mech_egy_mech_hhz
         call fs%pressure_relax(vf,matmod,relax_model)
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA

      end block create_and_initialize_flow_solver


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='ShockTube')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('T',fs%Tmptr)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('PA',fs%PA)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_scalar('Bulkmod',fs%RHOSS2)
         call ens_out%add_scalar('divU',fs%divU)
         call ens_out%add_scalar('Viscosity',fs%visc)
         call ens_out%add_scalar('Mach',fs%Mach)
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_scalar('epsp',lp%VF)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call lp%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Tmax,'Tmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'cvg')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(fs%impl_it_x,'Impl_x iteration')
         call cvgfile%add_column(fs%impl_rerr_x,'Impl_x error')
         call cvgfile%add_column(fs%impl_it_y,'Impl_y iteration')
         call cvgfile%add_column(fs%impl_rerr_y,'Impl_y error')
         call cvgfile%add_column(fs%implicit%it,'Impl_z iteration')
         call cvgfile%add_column(fs%implicit%rerr,'Impl_z error')
         call cvgfile%add_column(fs%psolv%it,'Pressure iteration')
         call cvgfile%add_column(fs%psolv%rerr,'Pressure error')
         ! Create LPT monitor
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%VFmean,'VFp_mean')
         call lptfile%add_column(lp%VFmax,'VFp_max')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
       end block create_monitor


   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()

         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Collide and advance particles
         call lp%collide(dt=time%dtmid)
         resU=fs%rho
         call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=fs%visc)

         ! Update density based on particle volume fraction
         fs%rho=resU*(1.0_WP-lp%VF)
         dRHOdt=(fs%RHO-fs%RHOold)/time%dtmid

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Add momentum source term from lpt (cell-centered)
            add_lpt_src: block
              integer :: i,j,k
              resU=0.0_WP; resV=0.0_WP; resW=0.0_WP; resE=0.0_WP
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       resU(i,j,k)=resU(i,j,k)+lp%srcU(i,j,k)
                       resV(i,j,k)=resV(i,j,k)+lp%srcV(i,j,k)
                       resW(i,j,k)=resW(i,j,k)+lp%srcW(i,j,k)
                       resE(i,j,k)=resE(i,j,k)+lp%srcE(i,j,k)
                    end do
                 end do
              end do
            end block add_lpt_src

            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Viscous step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod,srcU=resU,srcV=resV,srcW=resW,srcE=resE)

            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            fs%psolv%sol=fs%PA-fs%P
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)

            ! Record convergence monitor
            call cvgfile%write()
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do

         ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              integer :: i
              call lp%update_partmesh(pmesh)
              do i=1,lp%np_
                 pmesh%var(1,i)=0.5_WP*lp%p(i)%d
              end do
            end block update_pmesh
            call fs%get_viz()
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call lptfile%write()

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
      deallocate(resU,resV,resW,resE,dRHOdt)

   end subroutine simulation_final

end module simulation
