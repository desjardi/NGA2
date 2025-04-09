!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use iterator_class,    only: iterator
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final
   
   !> Flow solver objects
   type(hypre_str),   public :: ps     !< Structured Hypre linear solver for pressure
   type(ddadi),       public :: vs     !< DDADI solver for velocity
   type(tpns),        public :: fs     !< Two-phase flow solver
   type(vfs),         public :: vf     !< Volume fraction solver
   type(timetracker), public :: time   !< Time info
   
   !> SGS modeling
   logical        :: use_sgs   !< Is an LES model used?
   type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh     !< Surface mesh for interface
   type(ensight)  :: ens_out   !< Ensight output for flow variables
   type(event)    :: ens_evt   !< Event trigger for Ensight output
   
   !> Monitoring files
   type(monitor) :: mfile      !< General simulation monitoring
   type(monitor) :: cflfile    !< CFL monitoring
   
   !> Iterator for VOF removal
   type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
   integer        :: nlayer=4           !< Size of buffer layer for VOF removal
   real(WP)       :: vof_removed        !< Integral of VOF removed
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities (based on Umid)
   real(WP), dimension(:,:,:,:), allocatable :: vel               !< Other cell-centered velocity (based on U)
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      real(WP), dimension(3) :: center
      real(WP) :: radius,depth
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(vel (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: plicnet,VFhi,VFlo,flux
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')
         ! Initialize to a droplet and a pool
         radius=0.5_WP
         call param_read('Droplet center',center)
         call param_read('Pool depth',depth)
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_falling_drop,0.0_WP,amr_ref_lvl)
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
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         vof_removal_layer=iterator(cfg,'VOF removal',vof_removal_layer_locator)
         vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create a two-phase flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         real(WP) :: Re,Fr,We,r,m
         ! Create flow solver
         call fs%initialize(cfg=cfg,name='Two-phase NS')
         ! Add slight backward bias to CN scheme
         fs%theta=fs%theta+1.0e-2_WP
         ! Read in adimensional parameters
         call param_read('Froude number',Fr)     ! Fr=U^2/(g*D)=1/g
         call param_read('Weber number',We)      ! We=rho_l*U^2*D/sigma=1/sigma
         call param_read('Reynolds number',Re)   ! Re=rho_l*U*D/mu_l=1/mu_l
         call param_read('Density ratio',r)      ! r=rho_l/rho_g=1/rho_g
         call param_read('Viscosity ratio',m)    ! m=mu_l/mu_g
         ! Assign fluid properties to each phase
         fs%gravity=[0.0_WP,-Fr**(-1.0_WP),0.0_WP]
         fs%sigma=We**(-1.0_WP)
         fs%rho_l=1.0_WP
         fs%rho_g=fs%rho_l/r
         fs%visc_l=Re**(-1.0_WP)
         fs%visc_g=fs%visc_l/m
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=24
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Generate initial conditions for velocity
      initialize_velocity: block
         use vfs_class, only: VFlo
         integer :: i,j,k
         ! Initialize density
         resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_sRHO(rho=resU)
         ! Initialize velocity
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  if (cfg%ym(j).gt.depth.and.maxval(vf%VF(i,j-1:j,k)).gt.VFlo) fs%V(i,j,k)=-1.0_WP
               end do
            end do
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Copy to Umid and make it solenoidal
         call fs%get_Umid()
         call fs%correct_mfr()
         call fs%update_laplacian()
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%Umid=fs%Umid-resU/(fs%sRHOX**2); fs%U=fs%Umid
         fs%Vmid=fs%Vmid-resV/(fs%sRHOY**2); fs%V=fs%Vmid
         fs%Wmid=fs%Wmid-resW/(fs%sRHOZ**2); fs%W=fs%Wmid
         call fs%get_U()
         ! Calculate cell-centered velocities and divergence
         call fs%interp_velmid(Ui,Vi,Wi)
         call fs%interp_vel(vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',use_sgs)
         if (use_sgs) then
            allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         end if
      end block create_sgs
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='FallingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('othervel',vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
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
      
   contains
      
      !> Function that defines a level set function for a falling drop problem
      function levelset_falling_drop(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         ! Create the droplet
         G=radius-sqrt(sum((xyz-center)**2))
         ! Add the pool
         G=max(G,depth-xyz(2))
      end function levelset_falling_drop
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.ge.pg%jmax-nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocities and sRHOs
         fs%Uold=fs%U; fs%sRHOxold=fs%sRHOx
         fs%Vold=fs%V; fs%sRHOyold=fs%sRHOy
         fs%Wold=fs%W; fs%sRHOzold=fs%sRHOz
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! VOF equation ====================================================
            ! Advance VOF equation
            vf%VF=vf%VFold
            if (time%it.eq.time%itmax) then   
               call vf%advance(dt=time%dt,U=fs%Umid,V=fs%Vmid,W=fs%Wmid)
            else
               call vf%advance_tmp(dt=time%dt,U=fs%Umid,V=fs%Vmid,W=fs%Wmid)
            end if
            
            ! Update sqrt(face density) and momentum vector
            resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_sRHO(rho=resU)
            fs%rhoU=fs%rho_l*vf%UFl(1,:,:,:)+fs%rho_g*vf%UFg(1,:,:,:)
            fs%rhoV=fs%rho_l*vf%UFl(2,:,:,:)+fs%rho_g*vf%UFg(2,:,:,:)
            fs%rhoW=fs%rho_l*vf%UFl(3,:,:,:)+fs%rho_g*vf%UFg(3,:,:,:)
            
            ! Prepare new staggered viscosity (at n+1)
            call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
            
            ! Turbulence modeling
            if (use_sgs) then
               sgs_modeling: block
                  use sgsmodel_class, only: vreman
                  integer :: i,j,k
                  resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF)
                  call fs%get_gradUmid(gradU)
                  call sgs%get_visc(type=vreman,dt=time%dt,rho=resU,gradu=gradU)
                  do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_
                     do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_
                        do i=fs%cfg%imino_+1,fs%cfg%imaxo_
                           fs%visc(i,j,k)   =fs%visc(i,j,k)   +sgs%visc(i,j,k)
                           fs%visc_xy(i,j,k)=fs%visc_xy(i,j,k)+sum(fs%itp_xy(:,:,i,j,k)*sgs%visc(i-1:i,j-1:j,k))
                           fs%visc_yz(i,j,k)=fs%visc_yz(i,j,k)+sum(fs%itp_yz(:,:,i,j,k)*sgs%visc(i,j-1:j,k-1:k))
                           fs%visc_zx(i,j,k)=fs%visc_zx(i,j,k)+sum(fs%itp_xz(:,:,i,j,k)*sgs%visc(i-1:i,j,k-1:k))
                        end do
                     end do
                  end do
               end block sgs_modeling
            end if
            
            ! Momentum equation ===============================================
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-(fs%U*fs%sRHOX**2-fs%Uold*fs%sRHOXold**2)+time%dt*resU
            resV=-(fs%V*fs%sRHOY**2-fs%Vold*fs%sRHOYold**2)+time%dt*resV
            resW=-(fs%W*fs%sRHOZ**2-fs%Wold*fs%sRHOZold**2)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Compute predictor U
            fs%U=fs%U+resU
            fs%V=fs%V+resV
            fs%W=fs%W+resW
            
            ! Sync and apply boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Poisson equation ================================================
            ! Compute Umid from U and Uold
            call fs%get_Umid()
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct pressure and Umid
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%Umid=fs%Umid-time%dt*resU/((fs%sRHOX+fs%sRHOXold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOX)
            fs%Vmid=fs%Vmid-time%dt*resV/((fs%sRHOY+fs%sRHOYold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOY)
            fs%Wmid=fs%Wmid-time%dt*resW/((fs%sRHOZ+fs%sRHOZold*(1.0_WP-fs%theta)/fs%theta)*fs%sRHOZ)
            
            ! Regenerate U from Umid and Uold
            call fs%get_U()
            
            ! Increment sub-iteration counter =================================
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_velmid(Ui,Vi,Wi)
         call fs%interp_vel(vel(1,:,:,:),vel(2,:,:,:),vel(3,:,:,:))
         call fs%get_div()
         
         ! Remove VOF at edge of domain
         remove_vof: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            use parallel, only: MPI_REAL_WP
            integer :: n,i,j,k,ierr
            ! Adjust VOF field
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
            ! Also adjust density
            resU=fs%rho_l*vf%VF+fs%rho_g*(1.0_WP-vf%VF); call fs%update_sRHO(rho=resU)
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
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      if (use_sgs) deallocate(gradU)
      
   end subroutine simulation_final
   
   
end module simulation