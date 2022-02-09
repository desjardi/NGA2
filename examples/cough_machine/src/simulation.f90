!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg_in,t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Two-phase incompressible flow solver, VF solver with CCL, and corresponding time tracker and sgs model
   type(tpns),        public :: fs_in
   type(vfs),         public :: vf_in
   type(timetracker), public :: time_in
   type(sgsmodel),    public :: sgs_in
   
   !> Ensight postprocessing
   type(ensight)  :: ens_out_in
   type(event)    :: ens_evt_in
   
   !> Simulation monitor file
   type(monitor) :: mfile_in,cflfile_in
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU_in,resV_in,resW_in
   real(WP), dimension(:,:,:), allocatable :: Ui_in,Vi_in,Wi_in
   real(WP), dimension(:,:,:,:), allocatable :: SR_in
   
   !> Inflow parameters
   real(WP) :: Uin,delta,Urand
   
contains
   
   
   !> Function that localizes the left domain boundary, inside the mouth
   function left_boundary_mouth(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j).gt.0.0_WP.and.pg%ym(j).lt.H_mouth.and.abs(pg%zm(k)).lt.0.5_WP*W_mouth) isIn=.true.
   end function left_boundary_mouth
   
   
   !> Function that localizes the rightmost domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Zoomed in simulation setup
      ! ==========================
      
      ! Allocate work arrays for cfg
      allocate_work_arrays_in: block
         allocate(resU_in(cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
         allocate(resV_in(cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
         allocate(resW_in(cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
         allocate(Ui_in  (cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
         allocate(Vi_in  (cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
         allocate(Wi_in  (cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
         allocate(SR_in(6,cfg_in%imino_:cfg_in%imaxo_,cfg_in%jmino_:cfg_in%jmaxo_,cfg_in%kmino_:cfg_in%kmaxo_))
      end block allocate_work_arrays_in
      
      
      ! Initialize time tracker
      initialize_timetracker_in: block
         time_in=timetracker(cfg_in%amRoot,name='cough_machine_in')
         call param_read('Max timestep size',time_in%dtmax)
         call param_read('Max cfl number',time_in%cflmax)
         call param_read('Max time',time_in%tmax)
         time_in%dt=time_in%dtmax
         time_in%itmax=2
      end block initialize_timetracker_in
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: lvira,r2p
         integer :: i,j,k
         ! Create a VOF solver with LVIRA
         vf_in=vfs(cfg=cfg_in,reconstruction_method=lvira,name='VOF')
         ! Create a VOF solver with R2P
         !vf=vfs(cfg=cfg,reconstruction_method=r2p,name='VOF')
         ! Initialize to flat interface in liquid tray
         do k=vf_in%cfg%kmino_,vf_in%cfg%kmaxo_
            do j=vf_in%cfg%jmino_,vf_in%cfg%jmaxo_
               do i=vf_in%cfg%imino_,vf_in%cfg%imaxo_
                  if (vf_in%cfg%xm(i).lt.-L_lip.and.vf_in%cfg%xm(i).gt.-L_lip-L_film.and.abs(vf_in%cfg%zm(k)).lt.0.5_WP*W_film.and.vf_in%cfg%ym(j).lt.0.0_WP.and.vf_in%cfg%ym(j).gt.-H_film) then
                     vf_in%VF(i,j,k)=1.0_WP
                  else
                     vf_in%VF(i,j,k)=0.0_WP
                  end if
                  vf_in%Lbary(:,i,j,k)=[vf_in%cfg%xm(i),vf_in%cfg%ym(j),vf_in%cfg%zm(k)]
                  vf_in%Gbary(:,i,j,k)=[vf_in%cfg%xm(i),vf_in%cfg%ym(j),vf_in%cfg%zm(k)]
               end do
            end do
         end do
         ! Update the band
         call vf_in%update_band()
         ! Perform interface reconstruction from VOF field
         call vf_in%build_interface()
         ! Set interface planes at the boundaries
         call vf_in%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf_in%polygonalize_interface()
         ! Calculate distance from polygons
         call vf_in%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf_in%subcell_vol()
         ! Calculate curvature
         call vf_in%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf_in%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver with bconds
      create_solver_in: block
         use tpns_class, only: dirichlet,clipped_neumann,neumann
         use ils_class,  only: pcg_pfmg,gmres_amg
         use mathtools,  only: Pi
         ! Create a two-phase flow solver
         fs_in=tpns(cfg=cfg_in,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs_in%visc_l)
         call param_read('Gas dynamic viscosity'   ,fs_in%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs_in%rho_l)
         call param_read('Gas density'   ,fs_in%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs_in%sigma)
         call param_read('Static contact angle',fs_in%contact_angle)
         fs_in%contact_angle=fs_in%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs_in%gravity)
         ! Inflow on the left
         call fs_in%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_mouth)
         ! Outflow on the right
         call fs_in%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs_in%psolv%maxit)
         call param_read('Pressure tolerance',fs_in%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs_in%implicit%maxit)
         call param_read('Implicit tolerance',fs_in%implicit%rcvg)
         ! Setup the solver
         !fs%psolv%maxlevel=10
         call fs_in%setup(pressure_ils=gmres_amg,implicit_ils=gmres_amg)
      end block create_solver_in
      
      
      ! Initialize our velocity field
      initialize_velocity_in: block
         use tpns_class, only: bcond
         use random,     only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! Zero initial field
         fs_in%U=0.0_WP; fs_in%V=0.0_WP; fs_in%W=0.0_WP
         ! Apply Dirichlet at inlet
         call param_read('Gas velocity',Uin)
         call param_read('Gas thickness',delta)
         call param_read('Gas perturbation',Urand)
         call fs_in%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs_in%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(fs_in%cfg%zm(k)))/delta)*tanh(2.0_WP*fs_in%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-fs_in%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
         end do
         ! Apply all other boundary conditions
         call fs_in%apply_bcond(time_in%t,time_in%dt)
         ! Compute MFR through all boundary conditions
         call fs_in%get_mfr()
         ! Adjust MFR for global mass balance
         call fs_in%correct_mfr()
         ! Compute cell-centered velocity
         call fs_in%interp_vel(Ui_in,Vi_in,Wi_in)
         ! Compute divergence
         call fs_in%get_div()
      end block initialize_velocity_in
      
      
      ! Create an LES model
      create_sgs_in: block
         sgs_in=sgsmodel(cfg=fs_in%cfg,umask=fs_in%umask,vmask=fs_in%vmask,wmask=fs_in%wmask)
      end block create_sgs_in
      
      
      ! Add Ensight output
      create_ensight_in: block
         ! Create Ensight output from cfg
         ens_out_in=ensight(cfg_in,'cough_in')
         ! Create event for Ensight output
         ens_evt_in=event(time_in,'Ensight output')
         call param_read('Ensight output period',ens_evt_in%tper)
         ! Add variables to output
         call ens_out_in%add_vector('velocity',Ui_in,Vi_in,Wi_in)
         call ens_out_in%add_scalar('VOF',vf_in%VF)
         call ens_out_in%add_scalar('curvature',vf_in%curv)
         !call ens_out_in%add_scalar('visc_t',sgs_in%visc)
         ! Output to ensight
         if (ens_evt_in%occurs()) call ens_out_in%write_data(time_in%t)
      end block create_ensight_in
      
      
      ! Create a monitor file
      create_monitor_in: block
         ! Prepare some info about fields
         call fs_in%get_cfl(time_in%dt,time_in%cfl)
         call fs_in%get_max()
         call vf_in%get_max()
         ! Create simulation monitor
         mfile_in=monitor(fs_in%cfg%amRoot,'simulation')
         call mfile_in%add_column(time_in%n,'Timestep number')
         call mfile_in%add_column(time_in%t,'Time')
         call mfile_in%add_column(time_in%dt,'Timestep size')
         call mfile_in%add_column(time_in%cfl,'Maximum CFL')
         call mfile_in%add_column(fs_in%Umax,'Umax')
         call mfile_in%add_column(fs_in%Vmax,'Vmax')
         call mfile_in%add_column(fs_in%Wmax,'Wmax')
         call mfile_in%add_column(fs_in%Pmax,'Pmax')
         call mfile_in%add_column(vf_in%VFmax,'VOF maximum')
         call mfile_in%add_column(vf_in%VFmin,'VOF minimum')
         call mfile_in%add_column(vf_in%VFint,'VOF integral')
         call mfile_in%add_column(fs_in%divmax,'Maximum divergence')
         call mfile_in%add_column(fs_in%psolv%it,'Pressure iteration')
         call mfile_in%add_column(fs_in%psolv%rerr,'Pressure error')
         call mfile_in%write()
         ! Create CFL monitor
         cflfile_in=monitor(fs_in%cfg%amRoot,'cfl')
         call cflfile_in%add_column(time_in%n,'Timestep number')
         call cflfile_in%add_column(time_in%t,'Time')
         call cflfile_in%add_column(fs_in%CFLc_x,'Convective xCFL')
         call cflfile_in%add_column(fs_in%CFLc_y,'Convective yCFL')
         call cflfile_in%add_column(fs_in%CFLc_z,'Convective zCFL')
         call cflfile_in%add_column(fs_in%CFLv_x,'Viscous xCFL')
         call cflfile_in%add_column(fs_in%CFLv_y,'Viscous yCFL')
         call cflfile_in%add_column(fs_in%CFLv_z,'Viscous zCFL')
         call cflfile_in%write()
      end block create_monitor_in
      
      
      ! Zoomed out simulation setup
      ! ===========================
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: static_contact
      implicit none
      
      ! Perform time integration - the zoomed in solver is the main driver here
      do while (.not.time_in%done())
         
         ! Increment time
         call fs_in%get_cfl(time_in%dt,time_in%cfl)
         call time_in%adjust_dt()
         call time_in%increment()
         
         ! Remember old VOF
         vf_in%VFold=vf_in%VF
         
         ! Remember old velocity
         fs_in%Uold=fs_in%U
         fs_in%Vold=fs_in%V
         fs_in%Wold=fs_in%W
         
         ! Reapply Dirichlet at inlet
         reapply_dirichlet_in: block
            use tpns_class, only: bcond
            use random,     only: random_uniform
            type(bcond), pointer :: mybc
            integer  :: n,i,j,k
            call fs_in%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs_in%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(fs_in%cfg%zm(k)))/delta)*tanh(2.0_WP*fs_in%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-fs_in%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
            end do
         end block reapply_dirichlet_in
         
         ! Prepare old staggered density (at n)
         call fs_in%get_olddensity(vf=vf_in)
         
         ! VOF solver step
         call vf_in%advance(dt=time_in%dt,U=fs_in%U,V=fs_in%V,W=fs_in%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs_in%get_viscosity(vf=vf_in)
         
         ! Turbulence modeling - only work with gas properties here
         sgsmodel_in: block
            integer :: i,j,k
            call fs_in%get_strainrate(Ui=Ui_in,Vi=Vi_in,Wi=Wi_in,SR=SR_in)
            resU_in=fs_in%rho_g
            call sgs_in%get_visc(dt=time_in%dtold,rho=resU_in,Ui=Ui_in,Vi=Vi_in,Wi=Wi_in,SR=SR_in)
            where (sgs_in%visc.lt.-fs_in%visc_g)
               sgs_in%visc=-fs_in%visc_g
            end where
            do k=fs_in%cfg%kmino_+1,fs_in%cfg%kmaxo_
               do j=fs_in%cfg%jmino_+1,fs_in%cfg%jmaxo_
                  do i=fs_in%cfg%imino_+1,fs_in%cfg%imaxo_
                     fs_in%visc(i,j,k)   =fs_in%visc(i,j,k)   +sgs_in%visc(i,j,k)
                     fs_in%visc_xy(i,j,k)=fs_in%visc_xy(i,j,k)+sum(fs_in%itp_xy(:,:,i,j,k)*sgs_in%visc(i-1:i,j-1:j,k))
                     fs_in%visc_yz(i,j,k)=fs_in%visc_yz(i,j,k)+sum(fs_in%itp_yz(:,:,i,j,k)*sgs_in%visc(i,j-1:j,k-1:k))
                     fs_in%visc_zx(i,j,k)=fs_in%visc_zx(i,j,k)+sum(fs_in%itp_xz(:,:,i,j,k)*sgs_in%visc(i-1:i,j,k-1:k))
                  end do
               end do
            end do
         end block sgsmodel_in
         
         ! Perform sub-iterations
         do while (time_in%it.le.time_in%itmax)
            
            ! Build mid-time velocity
            fs_in%U=0.5_WP*(fs_in%U+fs_in%Uold)
            fs_in%V=0.5_WP*(fs_in%V+fs_in%Vold)
            fs_in%W=0.5_WP*(fs_in%W+fs_in%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs_in%prepare_advection_upwind(dt=time_in%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs_in%get_dmomdt(resU_in,resV_in,resW_in)
            
            ! Add momentum source terms
            call fs_in%addsrc_gravity(resU_in,resV_in,resW_in)
            
            ! Assemble explicit residual
            resU_in=-2.0_WP*fs_in%rho_U*fs_in%U+(fs_in%rho_Uold+fs_in%rho_U)*fs_in%Uold+time_in%dt*resU_in
            resV_in=-2.0_WP*fs_in%rho_V*fs_in%V+(fs_in%rho_Vold+fs_in%rho_V)*fs_in%Vold+time_in%dt*resV_in
            resW_in=-2.0_WP*fs_in%rho_W*fs_in%W+(fs_in%rho_Wold+fs_in%rho_W)*fs_in%Wold+time_in%dt*resW_in
            
            ! Form implicit residuals
            call fs_in%solve_implicit(time_in%dt,resU_in,resV_in,resW_in)
            
            ! Apply these residuals
            fs_in%U=2.0_WP*fs_in%U-fs_in%Uold+resU_in
            fs_in%V=2.0_WP*fs_in%V-fs_in%Vold+resV_in
            fs_in%W=2.0_WP*fs_in%W-fs_in%Wold+resW_in
            
            ! Apply other boundary conditions
            call fs_in%apply_bcond(time_in%t,time_in%dt)
            
            ! Solve Poisson equation
            call fs_in%update_laplacian()
            call fs_in%correct_mfr()
            call fs_in%get_div()
            call fs_in%add_surface_tension_jump(dt=time_in%dt,div=fs_in%div,vf=vf_in,contact_model=static_contact)
            fs_in%psolv%rhs=-fs_in%cfg%vol*fs_in%div/time_in%dt
            fs_in%psolv%sol=0.0_WP
            call fs_in%psolv%solve()
            call fs_in%shift_p(fs_in%psolv%sol)
            
            ! Correct velocity
            call fs_in%get_pgrad(fs_in%psolv%sol,resU_in,resV_in,resW_in)
            fs_in%P=fs_in%P+fs_in%psolv%sol
            fs_in%U=fs_in%U-time_in%dt*resU_in/fs_in%rho_U
            fs_in%V=fs_in%V-time_in%dt*resV_in/fs_in%rho_V
            fs_in%W=fs_in%W-time_in%dt*resW_in/fs_in%rho_W
            
            ! Increment sub-iteration counter
            time_in%it=time_in%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs_in%interp_vel(Ui_in,Vi_in,Wi_in)
         call fs_in%get_div()
         
         ! Output to ensight
         if (ens_evt_in%occurs()) call ens_out_in%write_data(time_in%t)
         
         ! Perform and output monitoring
         call fs_in%get_max()
         call vf_in%get_max()
         call mfile_in%write()
         call cflfile_in%write()
         
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
      deallocate(resU_in,resV_in,resW_in,Ui_in,Vi_in,Wi_in,SR_in)
      
   end subroutine simulation_final
   
   
end module simulation
