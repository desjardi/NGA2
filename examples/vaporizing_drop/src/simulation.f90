
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Lz
   use hypre_str_class,   only: hypre_str
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use evap_class,        only: evap
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps,psL
   type(tpns),        public :: fs,fsL
   type(vfs),         public :: vf
   type(evap),        public :: evp
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,mfileL,cflfile,evpfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Ui_L,Vi_L,Wi_L
   
   !> Problem definition
   real(WP), dimension(3) :: center
   real(WP) :: radius
   real(WP) :: mdotdp
   real(WP) :: rad_drop
   
contains


   !> Function that defines a level set function for a drop problem
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=radius-sqrt(sum((xyz-center)**2))
   end function levelset_drop


   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x- boundary for scalar fields
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y- boundary for scalar fields
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


   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator


   !> Function that localizes z- boundary for scalar fields
   function zm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin-1) isIn=.true.
   end function zm_locator_sc
   
   
   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui_L(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi_L(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi_L(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub-iterations',time%itmax)
         time%dt=time%dtmax
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,flux,neumann
         use mathtools, only: Pi
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux,name='VOF')
         vf%cons_correct=.false.
         ! Boundary conditinos
         call vf%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call vf%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call vf%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call vf%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         ! call vf%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         ! call vf%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Initialize to a droplet
         call param_read('Droplet center',center)
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_drop,0.0_WP,amr_ref_lvl)
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
         ! Droplet size
         call vf%get_max()
         rad_drop=sqrt(vf%VFint/(Lz*Pi))
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         use tpns_class,      only: clipped_neumann,bcond
         type(bcond), pointer :: mybc
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
         call param_read('Static contact angle',fs%contact_angle)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Boundary conditions
         call fs%add_bcond(name='xm',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fs%add_bcond(name='ym',type=clipped_neumann,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         ! call fs%add_bcond(name='zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! call fs%add_bcond(name='zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=16
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver


      ! Create a two-phase flow solver for the divergence-free liquid velocity
      create_div_free_flow_solver: block
         use tpns_class, only: clipped_neumann,bcond
         use hypre_str_class, only: pcg_pfmg2
         type(bcond), pointer :: mybc
         ! Create flow solver
         fsL=tpns(cfg=cfg,name='Liquid NS')
         ! Assign constant viscosity to each phase
         fsL%visc_l=fs%visc_l
         fsL%visc_g=fs%visc_g
         ! Assign constant density to each phase
         fsL%rho_l=fs%rho_l
         fsL%rho_g=fs%rho_g
         ! Assign surface tension coefficient
         fsL%sigma=fs%sigma
         fsL%contact_angle=fs%contact_angle
         ! Assign acceleration of gravity
         fsL%gravity=fs%gravity
         ! Boundary conditions
         call fsL%add_bcond(name='xm',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fsL%add_bcond(name='xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fsL%add_bcond(name='ym',type=clipped_neumann,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call fsL%add_bcond(name='yp',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         ! call fsL%add_bcond(name='zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! call fsL%add_bcond(name='zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         ! Configure pressure solver
         psL=hypre_str(cfg=cfg,name='Liquid P',method=pcg_pfmg2,nst=7)
         psL%maxlevel=ps%maxlevel
         psL%maxit=ps%maxit
         psL%rcvg=ps%rcvg
         ! Setup the solver
         call fsL%setup(pressure_solver=psL)
         ! Zero initial field
         fsL%U=0.0_WP; fsL%V=0.0_WP; fsL%W=0.0_WP
         ! Apply boundary conditions
         call fsL%apply_bcond(time%t,time%dt)
         ! Calculate cell-centered velocities and divergence
         call fsL%interp_vel(Ui_L,Vi_L,Wi_L)
         call fsL%get_div()
      end block create_div_free_flow_solver
      

      ! Create and initialize an evp object
      create_evp: block
         ! Create the object
         call evp%initialize(cfg=cfg,vf=vf,itp_x=fs%itpr_x,itp_y=fs%itpr_y,itp_z=fs%itpr_z,div_x=fs%divp_x,div_y=fs%divp_y,div_z=fs%divp_z,name='liquid gas pc')
         call param_read('Mass flux tolerence',     evp%mflux_tol)
         call param_read('Evaporation mass flux',   mdotdp)
         call param_read('Max pseudo timestep size',evp%pseudo_time%dtmax)
         call param_read('Max pseudo cfl number',   evp%pseudo_time%cflmax)
         call param_read('Max pseudo time steps',   evp%pseudo_time%nmax)
         evp%pseudo_time%dt=evp%pseudo_time%dtmax
         ! Get densities from the flow solver
         evp%rho_l=fs%rho_l
         evp%rho_g=fs%rho_g
         ! Interface jump condition
         where ((vf%VF.gt.0.0_WP).and.(vf%VF.lt.1.0_WP))
            evp%mdotdp=mdotdp
         else where
            evp%mdotdp=0.0_WP
         end where
         ! Get the volumetric evaporation mass flux
         call evp%get_mflux()
         ! Initialize the liquid and gas side mass fluxes
         call evp%init_mflux()
      end block create_evp
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='VaporizingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         call ens_out%add_scalar('mflux',evp%mflux)
         call ens_out%add_scalar('evp_div',evp%div_src)
         call ens_out%add_scalar('mdotdp',evp%mdotdp)
         call ens_out%add_scalar('mfluxL',evp%mfluxL)
         call ens_out%add_scalar('mfluxG',evp%mfluxG)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_vector('vel_itf',evp%U_itf,evp%V_itf,evp%W_itf)
         call ens_out%add_vector('vel_L',Ui_L,Vi_L,Wi_L)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call fsL%get_max()
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
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%add_column(rad_drop,'Droplet raduis')
         call mfile%write()
         ! Create simulation monitor for liquid
         mfileL=monitor(fsL%cfg%amRoot,'simulation_liquid')
         call mfileL%add_column(time%n,'Timestep number')
         call mfileL%add_column(time%t,'Time')
         call mfileL%add_column(time%dt,'Timestep size')
         call mfileL%add_column(time%cfl,'Maximum CFL')
         call mfileL%add_column(fsL%Umax,'Umax')
         call mfileL%add_column(fsL%Vmax,'Vmax')
         call mfileL%add_column(fsL%Wmax,'Wmax')
         call mfileL%add_column(fsL%Pmax,'Pmax')
         call mfileL%add_column(fsL%divmax,'Maximum divergence')
         call mfileL%add_column(fsL%psolv%it,'Pressure iteration')
         call mfileL%add_column(fsL%psolv%rerr,'Pressure error')
         call mfileL%write()
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
         ! Create evaporation monitor
         evpfile=monitor(evp%cfg%amRoot,'evaporation')
         call evpfile%add_column(time%n,'Timestep number')
         call evpfile%add_column(time%t,'Time')
         call evpfile%add_column(evp%pseudo_time%dt,'Pseudo time step')
         call evpfile%add_column(evp%pseudo_time%cfl,'Maximum pseudo CFL')
         call evpfile%add_column(evp%pseudo_time%n,'No. pseudo steps')
         call evpfile%add_column(evp%mflux_int,'mflux int')
         call evpfile%add_column(evp%mfluxL_int,'shifted mfluxL int')
         call evpfile%add_column(evp%mfluxG_int,'shifted mfluxG int')
         call evpfile%add_column(evp%mfluxL_int_err,'mfluxL int err')
         call evpfile%add_column(evp%mfluxG_int_err,'mfluxG int err')
         call evpfile%add_column(evp%mfluxL_err,'max mfluxL err')
         call evpfile%add_column(evp%mfluxG_err,'max mfluxG err')
         call evpfile%add_column(evp%Upcmax,'Upc_max')
         call evpfile%add_column(evp%Vpcmax,'Vpc_max')
         call evpfile%add_column(evp%Wpcmax,'Wpc_max')
         call evpfile%write()
      end block create_monitor
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact,harmonic_visc
      use mathtools, only: Pi
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
         fs%Uold=fs%U; fsL%Uold=fsL%U
         fs%Vold=fs%V; fsL%Vold=fsL%V
         fs%Wold=fs%W; fsL%Wold=fsL%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf); call fsL%get_olddensity(vf=vf)
         
         ! Interface jump conditions
         where ((vf%VF.gt.0.0_WP).and.(vf%VF.lt.1.0_WP))
            evp%mdotdp=mdotdp
         else where
            evp%mdotdp=0.0_WP
         end where

         ! Get the volumetric evaporation mass flux
         call evp%get_mflux()

         ! Get the interface velocity
         call evp%get_vel_pc()
         evp%U_itf=fsL%U-evp%vel_pc(:,:,:,1)
         evp%V_itf=fsL%V-evp%vel_pc(:,:,:,2)
         evp%W_itf=fsL%W-evp%vel_pc(:,:,:,3)

         ! VOF solver step
         call vf%advance(dt=time%dt,U=evp%U_itf,V=evp%V_itf,W=evp%W_itf)
         call vf%apply_bcond(time%t,time%dt)

         ! Shift the evaporation mass flux
         call evp%shift_mflux()

         ! Get the phase-change induced divergence
         call evp%get_div()

         ! Advance flow
         advance_flow: block
            ! Prepare new staggered viscosity (at n+1)
            call fs%get_viscosity(vf=vf,strat=harmonic_visc)
            
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
               
               ! Add momentum mass fluxes
               call fs%addsrc_gravity(resU,resV,resW)
               
               ! Assemble explicit residual
               resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
               resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
               resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
               
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho_U
               fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho_V
               fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho_W
               
               ! Apply other boundary conditions
               call fs%apply_bcond(time%t,time%dt)
               
               ! Solve Poisson equation
               call fs%update_laplacian()
               call fs%correct_mfr(src=evp%div_src)
               call fs%get_div(src=evp%div_src)
               ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
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
            call fs%get_div(src=evp%div_src)

         end block advance_flow

         ! Advance flow for the divergence-free velocity field
         advance_div_free_liquid: block
            integer  :: it_L,itmax_L

            ! Prepare new staggered viscosity (at n+1)
            call fsL%get_viscosity(vf=vf,strat=harmonic_visc)
            
            ! Perform sub-iterations
            itmax_L=2
            it_L=1            
            do while (it_L.le.itmax_L)
               
               ! Build mid-time velocity
               fsL%U=0.5_WP*(fsL%U+fsL%Uold)
               fsL%V=0.5_WP*(fsL%V+fsL%Vold)
               fsL%W=0.5_WP*(fsL%W+fsL%Wold)
               
               ! Preliminary mass and momentum transport step at the interface
               call fsL%prepare_advection_upwind(dt=time%dt)
               
               ! Explicit calculation of drho*u/dt from NS
               call fsL%get_dmomdt(resU,resV,resW)
               
               ! Add momentum mass fluxes
               call fsL%addsrc_gravity(resU,resV,resW)
               
               ! Assemble explicit residual
               resU=-2.0_WP*fsL%rho_U*fsL%U+(fsL%rho_Uold+fsL%rho_U)*fsL%Uold+time%dt*resU
               resV=-2.0_WP*fsL%rho_V*fsL%V+(fsL%rho_Vold+fsL%rho_V)*fsL%Vold+time%dt*resV
               resW=-2.0_WP*fsL%rho_W*fsL%W+(fsL%rho_Wold+fsL%rho_W)*fsL%Wold+time%dt*resW
               
               ! Apply these residuals
               fsL%U=2.0_WP*fsL%U-fsL%Uold+resU/fsL%rho_U
               fsL%V=2.0_WP*fsL%V-fsL%Vold+resV/fsL%rho_V
               fsL%W=2.0_WP*fsL%W-fsL%Wold+resW/fsL%rho_W
               
               ! Apply other boundary conditions
               call fsL%apply_bcond(time%t,time%dt)
               
               ! Solve Poisson equation
               call fsL%update_laplacian()
               call fsL%correct_mfr()
               call fsL%get_div()
               ! call fsL%add_surface_tension_jump(dt=time%dt,div=fsL%div,vf=vf,contact_model=static_contact)
               call fsL%add_surface_tension_jump(dt=time%dt,div=fsL%div,vf=vf)
               fsL%psolv%rhs=-fsL%cfg%vol*fsL%div/time%dt
               fsL%psolv%sol=0.0_WP
               call fsL%psolv%solve()
               call fsL%shift_p(fsL%psolv%sol)
               
               ! Correct velocity
               call fsL%get_pgrad(fsL%psolv%sol,resU,resV,resW)
               fsL%P=fsL%P+fsL%psolv%sol
               fsL%U=fsL%U-time%dt*resU/fsL%rho_U
               fsL%V=fsL%V-time%dt*resV/fsL%rho_V
               fsL%W=fsL%W-time%dt*resW/fsL%rho_W
               
               ! Increment sub-iteration counter
               it_L=it_L+1
               
            end do
            
            ! Recompute interpolated velocity and divergence
            call fsL%interp_vel(Ui_L,Vi_L,Wi_L)
            call fsL%get_div()

         end block advance_div_free_liquid
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max(); call fsL%get_max()
         call vf%get_max()
         call evp%get_max_vel_pc()
         rad_drop=sqrt(vf%VFint/(Lz*Pi))
         call mfile%write(); call mfileL%write()
         call cflfile%write()
         call evpfile%write()
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,Ui_L,Vi_L,Wi_L)

   end subroutine simulation_final
   
   
end module simulation