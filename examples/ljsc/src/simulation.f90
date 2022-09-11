!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use matm_class,        only: matm,air,water
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single two-phase flow solver, volume fraction solver, and material model set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(timetracker), public :: time

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Problem definition
   real(WP) :: djet
   real(WP), dimension(:), allocatable :: xjet
   integer :: relax_model, nwall

contains

   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imin-1) isIn=.true.
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
   
   !> Function that localizes the top (y+) of the domain
   function top_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain
   
   !> Function that localizes the back (z-) of the domain
   function back_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (k.eq.pg%kmin-1) isIn=.true.
   end function back_of_domain
   
   !> Function that localizes the front (z+) of the domain
   function front_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (k.eq.pg%kmax+1) isIn=.true.
   end function front_of_domain
   
   !> Function that localizes the jet(s) initial location
   function jet(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     real(WP), dimension(3) :: xyz
     logical :: isIn
     isIn=.false.
     xyz(1)=pg%xm(i); xyz(2)=0.0_WP; xyz(3)=pg%zm(k)
     if (j.le.pg%jmin-1+nwall.and.levelset_halfdrop(xyz,0.0_WP).gt.0.0_WP) isIn=.true.
   end function jet
   
   !> Function that localizes the walls surrounding the jets
   function wall(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.le.pg%jmin-1+nwall.and.(.not.jet(pg,i,j,k))) isIn=.true.
   end function wall
   
   !> Function that localizes the jet(s) BCs at edge of domain
   function jet_bdy(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     real(WP), dimension(3) :: xyz
     logical :: isIn
     isIn=.false.
     xyz(1)=pg%xm(i); xyz(2)=0.0_WP; xyz(3)=pg%zm(k)
     if (j.eq.pg%jmin-1.and.(jet(pg,i,j,k))) isIn=.true.
   end function jet_bdy

   !> Function that defines a level set function for the start of a jet at wall
   function levelset_halfdrop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      integer  :: n
      ! Loop jets
      G = -huge(1.0_WP)
      do n = 1,size(xjet)
        G=max(G,0.5_WP*djet-sqrt((xyz(1)-xjet(n))**2+max(0.0_WP,xyz(2))**2+xyz(3)**2))
      end do
   end function levelset_halfdrop

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_getsize
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
      
      ! Set up walls before solvers are initialized
      create_walls: block
        integer :: i,j,k,njet
        ! Initialize liquid jet(s)
        call param_read('Jet diameter',djet)
        njet = param_getsize('Jet location')
        allocate(xjet(njet))
        call param_read('Jet location',xjet)
        ! Number of wall cells
        call param_read('Wall cells in domain', nwall, default=1)
        do k=cfg%kmino_,cfg%kmaxo_
          do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
              if (wall(cfg%pgrid,i,j,k)) then
                cfg%VF(i,j,k)=0.0_WP
              end if
            end do
          end do
        end do
        
      end block create_walls


      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: elvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with lvira reconstruction
         vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
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
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_halfdrop,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  ! Round up to fully liquid within wall (accuracy of mdot improves with resolution)
                  if (vf%VF(i,j,k).gt.0.0_WP.and.vf%cfg%ym(j).lt.0.0_WP) vf%VF(i,j,k)=1.0_WP
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
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
         use mast_class, only: clipped_neumann,dirichlet,bc_scope,bcond,thermmech_egy_mech_hhz
         use ils_class,  only: pcg_bbox,gmres_smg
         use mathtools,  only: Pi
         use parallel,   only: amRoot
         integer :: i,j,k,n
         real(WP), dimension(3) :: xyz
         real(WP) :: gamm_l,Pref_l,q_l,b_l,gamm_g
         real(WP) :: ST, Ma, Ptot, Ttot, mdot, BL, vcross
         real(WP) :: T_cf, P_cf, U_cf, RHO_cf
         real(WP) :: T_j , P_j , V_j , RHO_j
         type(bcond), pointer :: mybc
         real(WP), parameter :: R = 287.0_WP ! Gas constant for air
         
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get EOS parameters from input
         ! Paramters from Kuhn and Desjardins (2021): gamm_l = 1.19,
         ! Pref_l = 7.028e8, q_l = -1.178e6, b_l = 6.61e-4, gamm_g = 1.4
         call param_read('Liquid gamma',gamm_l)
         call param_read('Liquid Pref', Pref_l)
         call param_read('Liquid q', q_l)
         call param_read('Liquid b', b_l)
         call param_read('Gas gamma',gamm_g)
         ! Register equations of state
         call matmod%register_NobleAbelstiffenedgas('liquid',gamm_l,Pref_l,q_l,b_l)
         call matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         ! Use built-in temperature-dependent models for diffusion and specific heat
         call matmod%register_diffusion_thermo_models(viscmodel_gas=air, viscmodel_liquid=water, &
                                                      hdffmodel_gas=air, hdffmodel_liquid=water, &
                                                      sphtmodel_gas=air, sphtmodel_liquid=water)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit momentum solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_bbox,implicit_ils=gmres_smg)

         
         ! Read crossflow parameters
         call param_read('Mach number',Ma)
         call param_read('Stagnation pressure',Ptot)
         call param_read('Stagnation temperature',Ttot)
         call param_read('Crossflow BL thickness',BL)
         ! Get crossflow quantities
         T_cf = Ttot / (1+(gamm_g-1.0_WP)/2.0_WP*Ma**2)
         P_cf = Ptot * (1.0_WP+(gamm_g-1.0_WP)/2.0_WP*Ma**2)**(gamm_g/(1.0_WP-gamm_g))
         U_cf = Ma * sqrt(gamm_g*R*T_cf)
         RHO_cf = P_cf/R/Ttot * (1+(gamm_g-1.0_WP)/2.0_WP*Ma**2)
         
         ! Read jet parameters
         call param_read('Mass flowrate',mdot)
         call param_read('Liquid temperature',T_j)
         call param_read('Liquid density',RHO_j)
         ! Get jet quantities
         V_j = mdot / RHO_j / (pi * 0.25_WP * djet**2)
         P_j = matmod%EOS_pressure(T_j,RHO_j,matmod%spec_heat_liquid(T_j),'liquid')
         
         ! Print useful information
         if (amRoot) then
           print*,"===== Problem Setup Description ====="
           print*,'Mach number', Ma, 'Number of jets', size(xjet), 'Mass flowrate', mdot
           print*,'Stagnation: Pressure', Ptot, 'Temperature', Ttot
           print*,'Crossflow:  Pressure', P_cf,'Density',RHO_cf,'Velocity',U_cf
           print*,'Jet:        Pressure',  P_j,'Density', RHO_j,'Velocity', V_j
           print*,'Re_BL',RHO_cf*U_cf*BL/matmod%viscosity_air(T_cf), &
                  'Re_D',RHO_cf*U_cf*djet/matmod%viscosity_air(T_cf), &
                  'Re_jet',RHO_j*V_j*djet/matmod%viscosity_water(T_j)
           if (fs%sigma.gt.0.0_WP) then
             print*,'We',RHO_cf*U_cf**2*djet/fs%sigma, &
                    'We_eff', (2.0_WP + (gamm_g-1.0_WP)*Ma**2) / & 
                    ((gamm_g-1.0_WP)*Ma**2) * (RHO_cf*U_cf**2*djet/fs%sigma)
           else
             print*,'We      Infinity'
           end if
         end if
         
         !! -- Populate crossflow data -- !!
         ! Initially 0 velocity in y and z
         fs%Vi = 0.0_WP; fs%Wi = 0.0_WP
         ! Zero face velocities as well for the sake of dirichlet boundaries
         fs%V = 0.0_WP; fs%W = 0.0_WP
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           if (j.le.nwall) then
             vcross = 0.0_WP
             RHO_cf = P_cf/R/Ttot
           else
             vcross = U_cf * min(1.0_WP, (fs%cfg%ym(j)/BL)**(1.0_WP/7.0_WP))
             RHO_cf = P_cf/R/Ttot/(1+(gamm_g-1.0_WP)/2.0_WP*Ma**2*(1.0_WP-(vcross/U_cf)**2))*(1+(gamm_g-1.0_WP)/2.0_WP*Ma**2)
           end if
           ! pressure, velocity, use matmod for energy
           fs%Grho(:,j,:) = RHO_cf
           fs%Ui(:,j,:) = vcross
           fs%GrhoE(:,j,:) = matmod%EOS_energy(P_cf,RHO_cf,vcross,0.0_WP,0.0_WP,'gas')
         end do
         
         !! -- Populate jet data -- !!
         fs%Lrho = RHO_j
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
           do j=fs%cfg%jmino_,fs%cfg%jmaxo_
             do i=fs%cfg%imino_,fs%cfg%imaxo_
               ! Populate liquid velocity (use uniform for now)
               if (vf%VF(i,j,k).gt.0.0_WP) then
                 fs%Ui(i,j,k) = 0.0_WP
                 fs%Vi(i,j,k) = V_j
                 fs%LrhoE(i,j,k) = matmod%EOS_energy(P_j,RHO_j,0.0_WP,V_j,0.0_WP,'liquid')
               end if
             end do
           end do
         end do
         
         ! Define boundary conditions - initialized values are intended dirichlet values too, for the cell centers
         call fs%add_bcond(name='crossflow',type=dirichlet      ,locator=left_of_domain ,celldir='xm')
         call fs%add_bcond(name='jet'      ,type=dirichlet      ,locator=jet_bdy        ,celldir='ym')
         call fs%add_bcond(name='outflow'  ,type=clipped_neumann,locator=right_of_domain,celldir='xp')
         call fs%add_bcond(name='out_yp'   ,type=clipped_neumann,locator=top_of_domain  ,celldir='yp')
         if (fs%cfg%nz.gt.1) then
           call fs%add_bcond(name='out_zm'   ,type=clipped_neumann,locator=back_of_domain ,celldir='zm')
           call fs%add_bcond(name='out_zp'   ,type=clipped_neumann,locator=front_of_domain,celldir='zp')
         end if

         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         ! Apply face BC - air inflow
         call fs%get_bcond('crossflow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i:i+1,j,k)=fs%Ui(i,j,k)
         end do
         ! Apply face BC - water inflow
         call fs%get_bcond('jet',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j:j+1,k)=fs%Vi(i,j,k)
         end do
         ! Apply face BC - outflows
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt,bc_scope)

         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         ! Perform initial pressure relax
         relax_model = thermmech_egy_mech_hhz
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
         ens_out=ensight(cfg=cfg,name='LiqJetinSSCrossflow')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('PA',fs%PA)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_scalar('Bulkmod',fs%RHOSS2)
         call ens_out%add_scalar('Temperature',fs%Tmptr)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
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
         call cflfile%add_column(fs%CFLst,'STension CFL')
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

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)
            
            ! Viscous step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

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
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
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

      ! Deallocate work arrays - none

   end subroutine simulation_final

end module simulation
