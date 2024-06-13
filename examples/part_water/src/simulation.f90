!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lpt_class,         only: lpt
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Two-phase flow solver, volume fraction solver, LPT solver, time tracker
   type(lpt),         public :: lp
   type(hypre_str),   public :: ps
   !type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: sfile,pfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: srcU,srcV,srcW
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,rho
   
   !> Problem definition
   real(WP) :: depth
   
   !> Max timestep size for LPT
   real(WP) :: lp_dt,lp_dt_max
   
contains
   
   
   !> Function that defines a level set function for a pool of water
   function levelset_pool(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Add the pool
      G=xyz(1)-depth
   end function levelset_pool
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(srcU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rho (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our LPT
      initialize_lpt: block
         use random, only: random_uniform
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Initialize with zero particles
         call lp%resize(0)
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Maximum timestep size used for particles
         call param_read('Particle timestep size',lp_dt_max,default=huge(1.0_WP))
         lp_dt=lp_dt_max
         ! Set collision timescale
         lp%tau_col=5.0_WP*time%dt
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n)
         call param_read('Wall restitution',lp%e_w)
         call param_read('Friction coefficient',lp%mu_f)
         ! Injection parameters
         call param_read('Particle mass flow rate',lp%mfr)
         call param_read('Particle velocity',lp%inj_vel)
         call param_read('Particle mean diameter',lp%inj_dmean)
         call param_read('Particle standard deviation',lp%inj_dsd,default=0.0_WP)
         call param_read('Particle min diameter',lp%inj_dmin,default=tiny(1.0_WP))
         call param_read('Particle max diameter',lp%inj_dmax,default=huge(1.0_WP))
         call param_read('Particle diameter shift',lp%inj_dshift,default=0.0_WP)
         if (lp%inj_dsd.le.epsilon(1.0_WP)) then
            lp%inj_dmin=lp%inj_dmean
            lp%inj_dmax=lp%inj_dmean
         end if
         call param_read('Particle inject diameter',lp%inj_d)
         lp%inj_pos(1)=lp%cfg%x(lp%cfg%imin)+lp%inj_dmax
         lp%inj_pos(2:3)=0.0_WP
         lp%inj_T=300.0_WP
      end block initialize_lpt
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: plicnet,VFhi,VFlo,remap
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=remap,name='VOF')
         ! Initialize to a pool
         call param_read('Pool depth',depth)
         depth=cfg%xL-depth
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_pool,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
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
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
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
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         !vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)!,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=1,name='lpt')
         pmesh%varname(1)='radius'
         pmesh%vecname(1)='velocity'
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(1,i)=0.5_WP*lp%p(i)%d
            pmesh%vec(:,1,i)=lp%p(i)%vel
         end do
      end block create_pmesh
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh

      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=fs%cfg,name='part_water')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_surface('plic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call lp%get_max()
         ! Create simulation monitor
         sfile=monitor(fs%cfg%amRoot,'simulation')
         call sfile%add_column(time%n,'Timestep number')
         call sfile%add_column(time%t,'Time')
         call sfile%add_column(time%dt,'Timestep size')
         call sfile%add_column(time%cfl,'Maximum CFL')
         call sfile%add_column(fs%Umax,'Umax')
         call sfile%add_column(fs%Vmax,'Vmax')
         call sfile%add_column(fs%Wmax,'Wmax')
         call sfile%add_column(fs%Pmax,'Pmax')
         call sfile%add_column(vf%VFmax,'VOF maximum')
         call sfile%add_column(vf%VFmin,'VOF minimum')
         call sfile%add_column(vf%VFint,'VOF integral')
         call sfile%add_column(fs%divmax,'Maximum divergence')
         call sfile%add_column(fs%psolv%it,'Pressure iteration')
         call sfile%add_column(fs%psolv%rerr,'Pressure error')
         call sfile%write()
         ! Create particle monitor
         pfile=monitor(amroot=lp%cfg%amRoot,name='particles')
         call pfile%add_column(time%n,'Timestep number')
         call pfile%add_column(time%t,'Time')
         call pfile%add_column(lp_dt,'Particle dt')
         call pfile%add_column(lp%np,'Particle number')
         call pfile%add_column(lp%np_new,'Npart new')
         call pfile%add_column(lp%np_out,'Npart removed')
         call pfile%add_column(lp%ncol,'Particle collisions')
         call pfile%add_column(lp%VFmax,'Max VF')
         call pfile%add_column(lp%Umin,'Particle Umin')
         call pfile%add_column(lp%Umax,'Particle Umax')
         call pfile%add_column(lp%Vmin,'Particle Vmin')
         call pfile%add_column(lp%Vmax,'Particle Vmax')
         call pfile%add_column(lp%Wmin,'Particle Wmin')
         call pfile%add_column(lp%Wmax,'Particle Wmax')
         call pfile%add_column(lp%dmin,'Particle dmin')
         call pfile%add_column(lp%dmax,'Particle dmax')
         call pfile%write()
         ! Create CFL monitor
         cflfile=monitor(lp%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(lp%CFLp_x,'Particle xCFL')
         call cflfile%add_column(lp%CFLp_y,'Particle yCFL')
         call cflfile%add_column(lp%CFLp_z,'Particle zCFL')
         call cflfile%add_column(lp%CFL_col,'Collision CFL')
         call cflfile%write()
      end block create_monitor
      

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
      implicit none
      real(WP) :: cfl
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Perform particle injection, collisions, and advancement
         
         call lp%collide(dt=time%dt)
         resU=vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g
         call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=fs%visc,srcU=srcU,srcV=srcV,srcW=srcW)
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Particle update
         lpt_step: block
            real(WP) :: dt_done,mydt
            real(WP), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3
            ! Allocate src storage
            allocate(tmp1(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp1=0.0_WP
            allocate(tmp2(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp2=0.0_WP
            allocate(tmp3(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp3=0.0_WP
            ! Get fluid stress
            !call fs%get_div_stress(resU,resV,resW)
            resU=0.0_WP; resV=0.0_WP; resW=0.0_WP
            ! Zero-out LPT source terms
            srcU=0.0_WP; srcV=0.0_WP; srcW=0.0_WP
            ! Sub-iterate
            call lp%get_cfl(lp_dt,cflc=cfl,cfl=cfl)
            if (cfl.gt.0.0_WP) lp_dt=min(lp_dt*time%cflmax/cfl,lp_dt_max)
            dt_done=0.0_WP
            do while (dt_done.lt.time%dtmid)
               ! Decide the timestep size
               mydt=min(lp_dt,time%dtmid-dt_done)
               ! Inject, collide and advance particles
               call lp%inject (dt=mydt,avoid_overlap=.true.)
               call lp%collide(dt=mydt)
               rho=vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g
               call lp%advance(dt=mydt,U=fs%U,V=fs%V,W=fs%W,rho=rho,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,srcU=tmp1,srcV=tmp2,srcW=tmp3)
               srcU=srcU+tmp1
               srcV=srcV+tmp2
               srcW=srcW+tmp3
               ! Increment
               dt_done=dt_done+mydt
            end do
            ! Update density based on particle volume fraction
            !fs%rho=rho*(1.0_WP-lp%VF)
            !dRHOdt=(fs%RHO-fs%RHOold)/time%dtmid
            ! Deallocate
            deallocate(tmp1,tmp2,tmp3)
         end block lpt_step
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Add momentum source term from lpt
            add_lpt_src: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
                  if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcU(i-1:i,j,k))
                  if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcV(i,j-1:j,k))
                  if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcW(i,j,k-1:k))
               end do; end do; end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
            end block add_lpt_src
            
            ! Form implicit residuals
            !call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho_U
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho_V
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho_W
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(1,i)=0.5_WP*lp%p(i)%d
                  pmesh%vec(:,1,i)=lp%p(i)%vel
               end do
            end block update_pmesh
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call lp%get_max()
         call sfile%write()
         call pfile%write()
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
      deallocate(srcU,srcV,srcW,resU,resV,resW,Ui,Vi,Wi,rho)
      
   end subroutine simulation_final
   
end module simulation
