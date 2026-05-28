!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use iterator_class,    only: iterator
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use resource_tracker,  only: getRSS
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Implicit solver
   logical     :: use_implicit !< Is an implicit solver used?
   type(ddadi) :: vs           !< DDADI solver for velocity   
   
   !> SGS modeling
   logical        :: use_sgs   !< Is an LES model used?
   type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
   
   !> Iterator for VOF removal
   type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
   real(WP)       :: vof_removed        !< Integral of VOF removed
   integer        :: nlayer=4           !< Size of buffer layer for VOF removal
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,memfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use mms_geom,  only: initialize_volume_moments
         use vfs_class, only: VFlo,plicnet,remap
         integer :: i,j,k
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=remap,name='VOF')
         ! Create the PLIC interface
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_; do j=vf%cfg%jmino_,vf%cfg%jmaxo_; do i=vf%cfg%imino_,vf%cfg%imaxo_
            ! Initialize liquid volume to zero and set corresponding PLIC
            vf%VF(i,j,k)=0.0_WP; vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]; vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
            call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1); call setPlane(vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,vf%VF(i,j,k)-0.5_WP))
            ! Not set volume moments for a liquid jet
            call initialize_volume_moments(lo=[vf%cfg%x(i),vf%cfg%y(j),vf%cfg%z(k)],hi=[vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)],&
            &                              levelset=levelset_jet,time=0.0_WP,level=5,VFlo=VFlo,VF=vf%VF(i,j,k),BL=vf%Lbary(:,i,j,k),BG=vf%Gbary(:,i,j,k))
         end do; end do; end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
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
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         vof_removal_layer=iterator(cfg,'VOF removal',vof_removal_layer_locator)
         vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create a two-phase flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: slip
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Read in flow conditions
         fs%rho_l=1.0_WP
         call param_read('Density ratio'  ,fs%rho_g);  fs%rho_g =1.0_WP/fs%rho_g
         call param_read('Reynolds number',fs%visc_l); fs%visc_l=1.0_WP/fs%visc_l
         call param_read('Viscosity ratio',fs%visc_g); fs%visc_g=fs%visc_l/fs%visc_g
         call param_read('Weber number'   ,fs%sigma);  fs%sigma =1.0_WP/fs%sigma
         ! Add slip conditions on the sides
         call fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         call fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         if (cfg%nz.gt.1) then
            call fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.false.,locator=zp_locator)
            call fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.false.,locator=zm_locator)
         end if
         ! Configure pressure solver
			ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Check if we want to use an implicit solver
         call param_read('Use implicit solver',use_implicit)
         if (use_implicit) then
            ! Configure implicit velocity solver
            vs=ddadi(cfg=cfg,name='Velocity',nst=7)
            ! Setup the solver
            call fs%setup(pressure_solver=ps,implicit_solver=vs)
         else
            ! Setup the solver
            call fs%setup(pressure_solver=ps)
         end if
      end block create_flow_solver
      
      
      ! Initialize velocity field
      initialize_velocity: block
         use random,    only: random_uniform
         use mathtools, only: Pi
         integer :: i,j,k
         real(WP) :: amp
         ! Initialize with powerlaw profile in the liquid jet normalized to Ubulk=1.0
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
            fs%U(i,j,k)=max(0.0_WP,(1.0_WP-sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)/0.5_WP)**(1.0_WP/7.0_WP))
         end do; end do; end do
         resU=vf%VF*fs%U; call cfg%integrate(A=resU,integral=amp); fs%U=0.25_WP*fs%U*cfg%xL*Pi/amp
         ! Add random disturbances
         call param_read('Perturbation amplitude',amp,default=0.05_WP)
         do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
            if (sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2).lt.0.5_WP) then
               fs%U(i,j,k)=fs%U(i,j,k)+random_uniform(lo=-amp,hi=+amp)
               fs%V(i,j,k)=            random_uniform(lo=-amp,hi=+amp)
               fs%W(i,j,k)=            random_uniform(lo=-amp,hi=+amp)
            end if
         end do; end do; end do
         call cfg%sync(fs%U); call cfg%sync(fs%V); call cfg%sync(fs%W)
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',use_sgs)
         if (use_sgs) sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='temp_jet')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_surface('plic',smesh)
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
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(vof_removed,'VOF removed')
         call mfile%add_column(vf%SDint,'SD integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
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
      end block create_monitor
      
      
      ! Track memory usage
      create_memory_monitor: block
         use resource_tracker, only: maxRSS,minRSS,avgRSS
         call getRSS()
         memfile=monitor(cfg%amRoot,'memory')
         call memfile%add_column(time%n,'Timestep number')
         call memfile%add_column(time%t,'Time')
         call memfile%add_column(maxRSS,'Maximum RSS')
         call memfile%add_column(minRSS,'Minimum RSS')
         call memfile%add_column(avgRSS,'Average RSS')
         call memfile%write()
      end block create_memory_monitor
      
      
   contains
      
      
      !> Function that defines a level set function for a cylindrical jet interface
      function levelset_jet(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-sqrt(xyz(2)**2+xyz(3)**2)
      end function levelset_jet
      
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      
      !> Function that localizes the back (z+) of the domain
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      
      !> Function that localizes the front (z-) of the domain
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.le.pg%jmin+nlayer.or.j.ge.pg%jmax-nlayer.or.k.le.pg%kmin+nlayer.or.k.ge.pg%kmax-nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())!.and.amp.lt.0.1_WP*vf%cfg%yL.and.time%t.lt.20.0_WP*tau)
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)
         
         ! Turbulence modeling
         if (use_sgs) then
            sgs_modeling: block
               use sgsmodel_class, only: vreman
               integer :: i,j,k
               resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF)
               call fs%get_gradu(gradU)
               call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
               do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_; do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_; do i=fs%cfg%imino_+1,fs%cfg%imaxo_
                  fs%visc   (i,j,k)=fs%visc   (i,j,k)+sgs%visc(i,j,k)
                  fs%visc_xy(i,j,k)=fs%visc_xy(i,j,k)+sum(fs%itp_xy(:,:,i,j,k)*sgs%visc(i-1:i,j-1:j,k))
                  fs%visc_yz(i,j,k)=fs%visc_yz(i,j,k)+sum(fs%itp_yz(:,:,i,j,k)*sgs%visc(i,j-1:j,k-1:k))
                  fs%visc_zx(i,j,k)=fs%visc_zx(i,j,k)+sum(fs%itp_xz(:,:,i,j,k)*sgs%visc(i-1:i,j,k-1:k))
               end do; end do; end do
            end block sgs_modeling
         end if
         
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
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
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
         
         ! Remove VOF at edge of domain
         remove_vof: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            use parallel, only: MPI_REAL_WP
            integer :: n,i,j,k,ierr
            vof_removed=0.0_WP
            do n=1,vof_removal_layer%no_
               i=vof_removal_layer%map(1,n)
               j=vof_removal_layer%map(2,n)
               k=vof_removal_layer%map(3,n)
               if (n.le.vof_removal_layer%n_) vof_removed=vof_removed+cfg%vol(i,j,k)*vf%VF(i,j,k)
               vf%VF(i,j,k)=0.0_WP
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,vof_removed,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
            call vf%clean_irl_and_band()
         end block remove_vof
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         
         ! Monitor memory usage
         call getRSS(); call memfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays
      deallocate(gradU,resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
