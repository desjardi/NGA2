!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,wheight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use fene_class,        only: fene
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs,ss
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(fene),        public :: nn
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: resSC,SCtmp
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU

   !> Problem definition
   real(WP) :: radius,Ujet
   
contains


   !> Function that defines a level set function for a falling drop problem
   function levelset_jet(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=radius-sqrt(xyz(1)**2+xyz(3)**2)
   end function levelset_jet
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(SCtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize the jet coming into the domain
         call param_read('Jet radius',radius)
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Initialize to circle at the top
                  if (j.gt.vf%cfg%jmax) then
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
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_jet,0.0_WP,amr_ref_lvl)
                     vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                     if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                        vf%Lbary(:,i,j,k)=v_cent
                        vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                     else
                        vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                        vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     end if
                  else
                     ! Inside, just set to zero
                     vf%VF(i,j,k)=0.0_WP
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
         use hypre_str_class, only: pcg_pfmg
         use tpns_class,      only: bcond,dirichlet,clipped_neumann
         type(bcond), pointer :: mybc
         real(WP) :: myr
         integer :: n,i,j,k
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
         ! Dirichlet inflow at the top
         call fs%add_bcond(name='bc_yp',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
         ! Outflow on the sides
         call fs%add_bcond(name='bc_xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fs%add_bcond(name='bc_xm',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='bc_zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call fs%add_bcond(name='bc_zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Setup jet inflow here
         call param_read('Jet velocity',Ujet)
         call fs%get_bcond('bc_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myr=sqrt(fs%cfg%xm(i)**2+fs%cfg%zm(k)**2)/radius
            fs%V(i,j,k)=min(0.0_WP,-2.0_WP*Ujet*(1.0_WP-myr**2))
         end do
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      

      ! Create a FENE model 
      create_fene: block 
         use multiscalar_class, only: bquick
         use fene_class,        only: fenecr,oldroydb,clipped_fenecr
         integer :: i,j,k
         ! Create FENE model solver
         nn=fene(cfg=cfg,model=fenecr,scheme=bquick,name='FENE')
         ! Assign unity density for simplicity
         nn%rho=1.0_WP
         ! Maximum extensibility of polymer chain
         call param_read('Maximum polymer extensibility',nn%Lmax)
         ! Relaxation time for polymer
         call param_read('Polymer relaxation time',nn%trelax)
         ! Polymer viscosity at zero strain rate
         call param_read('Polymer viscosity',nn%visc)
         ! Powerlaw coefficient in Carreau model
         call param_read('Carreau powerlaw',nn%ncoeff)
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='scalar',nst=13)
         ! Setup the solver
         call nn%setup(implicit_solver=ss)
         ! Initialize conformation tensor to identity
         nn%SC(:,:,:,1)=1.0_WP !< Cxx
         nn%SC(:,:,:,4)=1.0_WP !< Cyy
         nn%SC(:,:,:,6)=1.0_WP !< Czz
      end block create_fene
      
      
      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='jetBox')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('viscosity',fs%visc)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         do nsc=1,nn%nscalar
            call ens_out%add_scalar(trim(nn%SCname(nsc)),nn%SC(:,:,:,nsc))
         end do
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: nsc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call nn%get_max()
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
         ! Create scalar monitor
         scfile=monitor(nn%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         call scfile%add_column(nn%visc_pmax,'Maximum visc_p')
         call scfile%add_column(nn%visc_pmin,'Minimum visc_p')
         do nsc=1,nn%nscalar
            call scfile%add_column(nn%SCmin(nsc),trim(nn%SCname(nsc))//'_min')
            call scfile%add_column(nn%SCmax(nsc),trim(nn%SCname(nsc))//'_max')
         end do
         call scfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
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
         
         ! Remember old scalars
         nn%SCold=nn%SC
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Calculate grad(U)
         call fs%get_gradU(gradU)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! ============= SCALAR SOLVER =======================
            
            ! Reset interpolation metrics to QUICK scheme
            call nn%metric_reset()
            
            ! Build mid-time scalar
            nn%SC=0.5_WP*(nn%SC+nn%SCold)
            
            ! Explicit calculation of drhoSC/dt from scalar equation
            call nn%get_drhoSCdt(resSC,fs%Uold,fs%Vold,fs%Wold)
            
            ! Perform bquick procedure
            bquick: block
               integer :: i,j,k
               logical, dimension(:,:,:), allocatable :: flag
               ! Allocate work array
               allocate(flag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               ! Assemble explicit residual
               resSC=-2.0_WP*(nn%SC-nn%SCold)+time%dt*resSC
               ! Apply it to get explicit scalar prediction
               SCtmp=2.0_WP*nn%SC-nn%SCold+resSC
               ! Check cells that require bquick
               do k=nn%cfg%kmino_,nn%cfg%kmaxo_
                  do j=nn%cfg%jmino_,nn%cfg%jmaxo_
                     do i=nn%cfg%imino_,nn%cfg%imaxo_
                        if (SCtmp(i,j,k,1).le.0.0_WP.or.SCtmp(i,j,k,4).le.0.0_WP.or.SCtmp(i,j,k,6).le.0.0_WP.or.&
                        &   SCtmp(i,j,k,1)+SCtmp(i,j,k,4)+SCtmp(i,j,k,6).ge.nn%Lmax**2) then
                           flag(i,j,k)=.true.
                        else
                           flag(i,j,k)=.false.
                        end if
                     end do
                  end do
               end do
               ! Adjust metrics
               call nn%metric_adjust(SCtmp,flag)
               ! Clean up
               deallocate(flag)
               ! Recompute drhoSC/dt
               call nn%get_drhoSCdt(resSC,fs%Uold,fs%Vold,fs%Wold)
            end block bquick
            
            ! Add fene sources
            call nn%addsrc_CgradU(gradU,resSC)
            call nn%addsrc_relax(resSC,time%dt)
            
            ! Assemble explicit residual
            resSC=-2.0_WP*(nn%SC-nn%SCold)+time%dt*resSC
            
            ! Form implicit residual
            call nn%solve_implicit(time%dt,resSC,fs%Uold,fs%Vold,fs%Wold)
            
            ! Update scalars
            nn%SC=2.0_WP*nn%SC-nn%SCold+resSC
            
            ! Force the gas scalar to identity
            gas_scalar_forcing: block
               integer :: i,j,k
               do k=nn%cfg%kmino_,nn%cfg%kmaxo_
                  do j=nn%cfg%jmino_,nn%cfg%jmaxo_
                     do i=nn%cfg%imino_,nn%cfg%imaxo_
                        if (nn%mask(i,j,k).eq.0) then
                           nn%SC(i,j,k,1)=vf%VF(i,j,k)*nn%SC(i,j,k,1)+(1.0_WP-vf%VF(i,j,k))*1.0_WP
                           nn%SC(i,j,k,2)=vf%VF(i,j,k)*nn%SC(i,j,k,2)
                           nn%SC(i,j,k,3)=vf%VF(i,j,k)*nn%SC(i,j,k,3)
                           nn%SC(i,j,k,4)=vf%VF(i,j,k)*nn%SC(i,j,k,4)+(1.0_WP-vf%VF(i,j,k))*1.0_WP
                           nn%SC(i,j,k,5)=vf%VF(i,j,k)*nn%SC(i,j,k,5)
                           nn%SC(i,j,k,6)=vf%VF(i,j,k)*nn%SC(i,j,k,6)+(1.0_WP-vf%VF(i,j,k))*1.0_WP
                        end if
                     end do
                  end do
               end do
            end block gas_scalar_forcing
            
            ! Apply all other boundary conditions on the resulting field
            call nn%apply_bcond(time%t,time%dt)
            ! ===================================================
            
            
            ! ============ VELOCITY SOLVER ======================
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Include shear-thinning effect here by adjusting viscosity based on mid-time strain-rate
            ! fs%visc_l is the solvent viscosity, nn%visc is the zero strainrate polymer viscosity
            shear_thinning: block
               integer :: i,j,k
               real(WP) :: liq_vol,gas_vol,tot_vol
               real(WP) :: visc_l
               real(WP), dimension(:,:,:,:), allocatable :: SR
               ! Allocate SR array
               allocate(SR(1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               ! Calculate strain rate
               call fs%get_strainrate(SR)
               ! Update polymer viscosity using Carreau model
               call nn%update_visc_p(SR)
               ! Handle mixture viscosity
               do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_
                  do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_
                     do i=fs%cfg%imino_+1,fs%cfg%imaxo_
                        ! VISC at [xm,ym,zm] - direct sum in x/y/z
                        liq_vol=sum(vf%Lvol(:,:,:,i,j,k))
                        gas_vol=sum(vf%Gvol(:,:,:,i,j,k))
                        tot_vol=gas_vol+liq_vol
                        visc_l=fs%visc_l+nn%visc_p(i,j,k)
                        fs%visc(i,j,k)=0.0_WP; if (tot_vol.gt.0.0_WP) fs%visc(i,j,k)=(visc_l*liq_vol+fs%visc_g*gas_vol)/tot_vol
                        ! VISC_xy at [x,y,zm] - direct sum in z, staggered sum in x/y
                        liq_vol=sum(vf%Lvol(0,0,:,i,j,k))+sum(vf%Lvol(1,0,:,i-1,j,k))+sum(vf%Lvol(0,1,:,i,j-1,k))+sum(vf%Lvol(1,1,:,i-1,j-1,k))
                        gas_vol=sum(vf%Gvol(0,0,:,i,j,k))+sum(vf%Gvol(1,0,:,i-1,j,k))+sum(vf%Gvol(0,1,:,i,j-1,k))+sum(vf%Gvol(1,1,:,i-1,j-1,k))
                        tot_vol=gas_vol+liq_vol
                        visc_l=fs%visc_l+sum(fs%itp_xy(:,:,i,j,k)*nn%visc_p(i-1:i,j-1:j,k))
                        fs%visc_xy(i,j,k)=0.0_WP; if (tot_vol.gt.0.0_WP) fs%visc_xy(i,j,k)=(visc_l*liq_vol+fs%visc_g*gas_vol)/tot_vol
                        ! VISC_yz at [xm,y,z] - direct sum in x, staggered sum in y/z
                        liq_vol=sum(vf%Lvol(:,0,0,i,j,k))+sum(vf%Lvol(:,1,0,i,j-1,k))+sum(vf%Lvol(:,0,1,i,j,k-1))+sum(vf%Lvol(:,1,1,i,j-1,k-1))
                        gas_vol=sum(vf%Gvol(:,0,0,i,j,k))+sum(vf%Gvol(:,1,0,i,j-1,k))+sum(vf%Gvol(:,0,1,i,j,k-1))+sum(vf%Gvol(:,1,1,i,j-1,k-1))
                        tot_vol=gas_vol+liq_vol
                        visc_l=fs%visc_l+sum(fs%itp_yz(:,:,i,j,k)*nn%visc_p(i,j-1:j,k-1:k))
                        fs%visc_yz(i,j,k)=0.0_WP; if (tot_vol.gt.0.0_WP) fs%visc_yz(i,j,k)=(visc_l*liq_vol+fs%visc_g*gas_vol)/tot_vol
                        ! VISC_zx at [x,ym,z] - direct sum in y, staggered sum in z/x
                        liq_vol=sum(vf%Lvol(0,:,0,i,j,k))+sum(vf%Lvol(0,:,1,i,j,k-1))+sum(vf%Lvol(1,:,0,i-1,j,k))+sum(vf%Lvol(1,:,1,i-1,j,k-1))
                        gas_vol=sum(vf%Gvol(0,:,0,i,j,k))+sum(vf%Gvol(0,:,1,i,j,k-1))+sum(vf%Gvol(1,:,0,i-1,j,k))+sum(vf%Gvol(1,:,1,i-1,j,k-1))
                        tot_vol=gas_vol+liq_vol
                        visc_l=fs%visc_l+sum(fs%itp_xz(:,:,i,j,k)*nn%visc_p(i-1:i,j,k-1:k))
                        fs%visc_zx(i,j,k)=0.0_WP; if (tot_vol.gt.0.0_WP) fs%visc_zx(i,j,k)=(visc_l*liq_vol+fs%visc_g*gas_vol)/tot_vol
                     end do
                  end do
               end do
               ! Deallocate SR array
               deallocate(SR)
            end block shear_thinning
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Add polymer stress term
            polymer_stress: block
               integer :: i,j,k,n
               real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
               real(WP), dimension(:,:,:,:), allocatable :: stress
               ! Allocate work arrays
               allocate(stress(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
               allocate(Txy   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(Tyz   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(Tzx   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               ! Calculate the polymer relaxation
               stress=0.0_WP; call nn%addsrc_relax(stress,time%dt)
               ! Build liquid stress tensor
               do n=1,6
                  stress(:,:,:,n)=-nn%visc_p(:,:,:)*vf%VF*stress(:,:,:,n)
               end do
               ! Interpolate tensor components to cell edges
               do k=cfg%kmin_,cfg%kmax_+1
                  do j=cfg%jmin_,cfg%jmax_+1
                     do i=cfg%imin_,cfg%imax_+1
                        Txy(i,j,k)=sum(fs%itp_xy(:,:,i,j,k)*stress(i-1:i,j-1:j,k,2))
                        Tyz(i,j,k)=sum(fs%itp_yz(:,:,i,j,k)*stress(i,j-1:j,k-1:k,5))
                        Tzx(i,j,k)=sum(fs%itp_xz(:,:,i,j,k)*stress(i-1:i,j,k-1:k,3))
                     end do
                  end do
               end do
               ! Add divergence of stress to residual
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%divu_x(:,i,j,k)*stress(i-1:i,j,k,1))&
                        &                                                +sum(fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     &
                        &                                                +sum(fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
                        if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
                        &                                                +sum(fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))&
                        &                                                +sum(fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
                        if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
                        &                                                +sum(fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     &                  
                        &                                                +sum(fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6))        
                     end do
                  end do
               end do
               ! Clean up
               deallocate(stress,Txy,Tyz,Tzx)
            end block polymer_stress

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
            
            ! Solve Poisson equation - pinned version
            call fs%update_laplacian(pinpoint=[fs%cfg%imin,fs%cfg%jmin,fs%cfg%kmin])
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            if (cfg%amRoot) fs%psolv%rhs(cfg%imin,cfg%jmin,cfg%kmin)=0.0_WP
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
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call nn%get_max()
         call mfile%write()
         call cflfile%write()
         call scfile%write()

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
      deallocate(resSC,SCtmp,gradU)
      
   end subroutine simulation_final
   
   
   !> Function that localizes the y+ side of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes the x+ side of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1.and.pg%ym(j).gt.wheight) isIn=.true.
   end function xp_locator
   
   
   !> Function that localizes the x- side of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j).gt.wheight) isIn=.true.
   end function xm_locator
   
   
   !> Function that localizes the z+ side of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1.and.pg%ym(j).gt.wheight) isIn=.true.
   end function zp_locator
   
   
   !> Function that localizes the z- side of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin.and.pg%ym(j).gt.wheight) isIn=.true.
   end function zm_locator
   
   
end module simulation