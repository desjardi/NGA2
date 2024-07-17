!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use tpscalar_class,    only: tpscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   !type(ddadi),       public :: vs
   type(ddadi),       public :: ss
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(tpscalar),    public :: sc
   type(timetracker), public :: time,pseudo_time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out,ens_mdot_out
   type(event)    :: ens_evt,ens_mdot_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,mdotfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: VFgradX,VFgradY,VFgradZ
   real(WP), dimension(:,:,:),   allocatable :: mdot
   real(WP), dimension(:,:,:),   allocatable :: mdotL,mdotL_old,resmdotL,mdotL_err_field
   real(WP), dimension(:,:,:),   allocatable :: mdotG,mdotG_old,resmdotG
   
   !> Problem definition
   real(WP), dimension(3) :: center
   real(WP) :: radius,depth
   integer  :: iTl,iTg
   real(WP) :: mdot_int,mdot_err,mdot_tol
   real(WP) :: mdotL_int,mdotL_err,mdotL_int_err
   real(WP) :: mdotG_int,mdotG_err,mdotG_int_err
   real(WP) :: mdot_ens_time
   real(WP) :: evp_src
   
contains


   !> Function that defines a level set function for a drop problem
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=radius-sqrt(sum((xyz-center)**2))
      ! Add the pool
      ! G=max(G,depth-xyz(2))
   end function levelset_drop
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resSC    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:2))
         allocate(VFgradX  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFgradX=0.0_WP
         allocate(VFgradY  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFgradY=0.0_WP
         allocate(VFgradZ  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFgradZ=0.0_WP
         allocate(mdot     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mdot =0.0_WP
         allocate(mdotL    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mdotL=0.0_WP
         allocate(mdotG    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mdotG=0.0_WP
         allocate(mdotL_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mdotL_old=0.0_WP
         allocate(mdotG_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mdotG_old=0.0_WP
         allocate(resmdotL (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resmdotL=0.0_WP
         allocate(resmdotG (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resmdotG=0.0_WP
         allocate(mdotL_err_field(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mdotL_err_field=0.0_WP
         allocate(resU     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui       (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi       (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi       (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,remap,flux_storage,neumann
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux_storage,name='VOF')
         ! call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=remap,name='VOF')
         vf%cons_correct=.false.
         ! Initialize to a droplet and a pool
         !center=[0.0_WP,0.05_WP,0.0_WP]
         call param_read('Droplet center',center)
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
         depth =0.02_WP
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
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
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
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=10
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
      
      
      ! Create a one-sided scalar solver
      create_scalar: block
         use param, only: param_read
         use tpscalar_class, only: Lphase,Gphase
         use hypre_str_class, only: pcg_pfmg2,pfmg,gmres_pfmg
         integer :: i,j,k
         real(WP) :: Ldiff,Gdiff
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=2,name='tpscalar')
         ! Initialize the phase specific VOF
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Initialize the phase specific density
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         ! Temperature on the liquid and gas sides
         sc%SCname=[  'Tl',  'Tg']; iTl=1; iTg=2
         sc%phase =[Lphase,Gphase]
         ! Read diffusivity
         call param_read('Liquid thermal diffusivity',Ldiff)
         sc%diff(:,:,:,iTl)=Ldiff
         call param_read('Gas thermal diffusivity',Gdiff)
         sc%diff(:,:,:,iTg)=Gdiff
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=7)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize scalar fields
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Liquid scalar
                  if (vf%VF(i,j,k).gt.0.0_WP) then
                     ! We are in the liquid
                     if (cfg%ym(j).gt.depth+cfg%dy(j)) then
                        ! We are above the pool
                        sc%SC(i,j,k,iTl)=1.0_WP
                     else
                        ! We are in the pool
                        sc%SC(i,j,k,iTl)=2.0_WP
                     end if
                  end if
                  ! Gas scalar
                  if (vf%VF(i,j,k).lt.1.0_WP) then
                     ! We are in the gas
                     sc%SC(i,j,k,iTg)=(cfg%ym(j)-depth)/(cfg%yL-depth)
                  end if
               end do
            end do
         end do
      end block create_scalar
      
      ! Create a framework for shifting the mass source term
      initialize_mdot: block
         ! Initialize a pseudo time tracker
         pseudo_time=timetracker(amRoot=cfg%amRoot,name='Pseudo',print_info=.false.)
         call param_read('Max pseudo timestep size',pseudo_time%dtmax)
         call param_read('Max pseudo cfl number',pseudo_time%cflmax)
         call param_read('Max pseudo time steps',pseudo_time%nmax)
         call param_read('Tolerence',mdot_tol)
         call param_read('Evaporation surface mass rate',evp_src)
         pseudo_time%dt=pseudo_time%dtmax
         ! Initialize the evaporation mass source term and errors
         where ((vf%VF.gt.0.0_WP).and.(vf%VF.lt.1.0_WP)); mdot=evp_src; else where; mdot=0.0_WP; end where
         mdotL=mdot; mdotG=mdot
         mdot_err=0.0_WP; mdotL_err=0.0_WP; mdotG_err=0.0_WP
         mdotL_int_err=0.0_WP; mdotG_int_err=0.0_WP
         ! Integral
         call cfg%integrate(mdot,mdot_int)
         call cfg%integrate(mdotL,mdotL_int)
         call cfg%integrate(mdotG,mdotG_int)
      end block initialize_mdot
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         logical  :: ens_for_mdot

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
         do nsc=1,sc%nscalar
           call ens_out%add_scalar(trim(sc%SCname(nsc)),sc%SC(:,:,:,nsc))
         end do
         call ens_out%add_scalar('mdot' ,mdot)
         call ens_out%add_scalar('mdotL',mdotL)
         call ens_out%add_scalar('mdotG',mdotG)
         call ens_out%add_scalar('mdotL_err_field',mdotL_err_field)
         call ens_out%add_vector('VFgrad',VFgradx,VFgrady,VFgradz)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Creat ensight for mdot evolution
         call param_read('Ensight for mdot evolution',ens_for_mdot)
         if (ens_for_mdot) then
            call param_read('When to write mdot ensight',mdot_ens_time)
            ens_mdot_out=ensight(cfg=cfg,name='mdotEvolution')
            ens_mdot_evt=event(time=pseudo_time,name='mdot ensight')
            call param_read('Ensight output period for mdot',ens_mdot_evt%tper)
            call ens_mdot_out%add_scalar('VOF'  ,vf%VF)
            call ens_mdot_out%add_scalar('mdot' ,mdot)
            call ens_mdot_out%add_scalar('mdotL',mdotL)
            call ens_mdot_out%add_scalar('mdotG',mdotG)
            call ens_mdot_out%add_scalar('mdotL_err_field',mdotL_err_field)
            if (ens_mdot_evt%occurs()) call ens_mdot_out%write_data(pseudo_time%t)
         else
            mdot_ens_time=-1.0_WP
         end if
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: nsc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call sc%get_max(VF=vf%VF)
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
         scfile=monitor(sc%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         do nsc=1,sc%nscalar
           call scfile%add_column(sc%SCmin(nsc),trim(sc%SCname(nsc))//'_min')
           call scfile%add_column(sc%SCmax(nsc),trim(sc%SCname(nsc))//'_max')
           call scfile%add_column(sc%SCint(nsc),trim(sc%SCname(nsc))//'_int')
         end do
         call scfile%write()
         ! Create mdot monitor
         mdotfile=monitor(sc%cfg%amRoot,'mdot')
         call mdotfile%add_column(time%n,'Timestep number')
         call mdotfile%add_column(time%t,'Time')
         call mdotfile%add_column(pseudo_time%t,'Pseudo time')
         call mdotfile%add_column(pseudo_time%dt,'Pseudo time step')
         call mdotfile%add_column(pseudo_time%cfl,'Maximum pseudo CFL')
         call mdotfile%add_column(pseudo_time%n,'No. pseudo steps')
         call mdotfile%add_column(mdot_int,'mdot int')
         call mdotfile%add_column(mdotL_int,'shifted mdotL int')
         call mdotfile%add_column(mdotG_int,'shifted mdotG int')
         call mdotfile%add_column(mdotL_int_err,'mdotL int err')
         call mdotfile%add_column(mdotG_int_err,'mdotG int err')
         call mdotfile%add_column(mdotL_err,'max mdotL err')
         call mdotfile%add_column(mdotG_err,'max mdotG err')
         call mdotfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact,harmonic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old SC
         sc%SCold =sc%SC
         sc%PVFold=sc%PVF
         
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

         ! Update the evaporation source term
         where ((vf%VF.gt.0.0_WP).and.(vf%VF.lt.1.0_WP)); mdot=evp_src; else where; mdot=0.0_WP; end where
         
         ! ! Now transport our phase-specific scalars
         ! advance_scalar: block
         !    use tpscalar_class, only: Lphase,Gphase
         !    integer :: nsc
            
         !    ! Get the phas-specific VOF
         !    sc%PVF(:,:,:,Lphase)=vf%VF
         !    sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
            
         !    ! Explicit calculation of dSC/dt from advective term
         !    call sc%get_dSCdt(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,VFold=vf%VFold,VF=vf%VF,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)
            
         !    ! Advance advection
         !    do nsc=1,sc%nscalar
         !       where (sc%mask.eq.0.and.sc%PVF(:,:,:,sc%phase(nsc)).gt.0.0_WP) sc%SC(:,:,:,nsc)=((sc%PVFold(:,:,:,sc%phase(nsc)))*sc%SCold(:,:,:,nsc)+time%dt*resSC(:,:,:,nsc))/(sc%PVF(:,:,:,sc%phase(nsc)))
         !       where (sc%PVF(:,:,:,sc%phase(nsc)).eq.0.0_WP) sc%SC(:,:,:,nsc)=0.0_WP
         !    end do

         !    ! Apply the mass/energy transfer source term for the species/temperature
         !    do nsc=1,sc%nscalar
         !       where (sc%PVF(:,:,:,sc%phase(nsc)).gt.0.0_WP.and.sc%PVF(:,:,:,sc%phase(nsc)).lt.1.0_WP) sc%SC(:,:,:,nsc)=sc%SC(:,:,:,nsc)-mdot/sc%Prho(sc%phase(nsc))*vf%SD*time%dt
         !    end do
               
         !    ! Advance diffusion
         !    call sc%solve_implicit(time%dt,sc%SC)
            
         !    ! Apply boundary conditions
         !    call sc%apply_bcond(time%t,time%dt)
               
         ! end block advance_scalar         
               
         ! Shift the mass source term away from the interface
         shift_mdot: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
            use parallel, only: MPI_REAL_WP
            use irl_fortran_interface, only: calculateNormal,getNumberOfVertices
            real(WP), dimension(:,:,:), allocatable :: ccVFgradX,ccVFgradY,ccVFgradZ
            integer  :: ierr,i,j,k
            real(WP) :: my_mdot_err
            real(WP), dimension(3) :: n1,n2
            
            ! Allocate memory for the cell-centered VOF gradient
            allocate(ccVFgradX(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            allocate(ccVFgradY(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            allocate(ccVFgradZ(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            
            ! Get the scaled gradient of VOF
            do k=vf%cfg%kmino_,vf%cfg%kmaxo_
               do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                  do i=vf%cfg%imino_,vf%cfg%imaxo_
                     n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                        n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                        n1=0.5_WP*(n1+n2)
                     end if
                     ccVFgradX(i,j,k)=-n1(1)
                     ccVFgradY(i,j,k)=-n1(2)
                     ccVFgradZ(i,j,k)=-n1(3)
                  end do
               end do
            end do
            call sc%cellVec_to_face(ccf_x=ccVFgradX,ccf_y=ccVFgradY,ccf_z=ccVFgradZ,fcf_x=VFgradX,fcf_y=VFgradY,fcf_z=VFgradZ)
            
            ! Deallocate the unused
            deallocate(ccVFgradX,ccVFgradY,ccVFgradZ)
            
            ! Get the CFL based on the gradient of the VOF
            call sc%get_cfl(VFgradX,VFgradY,VFgradZ,pseudo_time%dt,pseudo_time%cfl)
            
            ! Reset the pseudo time
            call pseudo_time%reset()
            
            ! Adjust the pseudo time step
            call pseudo_time%adjust_dt()
            
            ! Initialize the source terms on the liquid and gas sides
            mdotL=mdot; mdotG=mdot

            ! Move the evaporation source term away from the interface
            do while (.not.pseudo_time%done())
               
               ! Remember old mdot
               mdotL_old=mdotL
               mdotG_old=mdotG
               
               ! Increment pseudo time
               call pseudo_time%increment()
               
               ! Assemble explicit residual
               call sc%get_dmdotdtau( VFgradX, VFgradY, VFgradZ,mdotL_old,resmdotL)
               call sc%get_dmdotdtau(-VFgradX,-VFgradY,-VFgradZ,mdotG_old,resmdotG)
               
               ! Apply these residuals
               mdotL=mdotL_old+pseudo_time%dt*resmdotL
               mdotG=mdotG_old+pseudo_time%dt*resmdotG

               ! Output to ensight
               mdotL_err_field=mdotL-mdotL_old
               if (mdot_ens_time.eq.time%t.and.ens_mdot_evt%occurs()) then
                  call ens_mdot_out%write_data(pseudo_time%t)
               end if

               ! Calculate the error on the liquid side
               my_mdot_err=maxval(abs(mdotL_err_field))
               call MPI_ALLREDUCE(my_mdot_err,mdotL_err,1,MPI_REAL_WP,MPI_Max,sc%cfg%comm,ierr)
               ! Calculate the error on the gas side
               my_mdot_err=maxval(abs(mdotG-mdotG_old))
               call MPI_ALLREDUCE(my_mdot_err,mdotG_err,1,MPI_REAL_WP,MPI_Max,sc%cfg%comm,ierr)
               ! Check convergence
               mdot_err=max(mdotL_err,mdotG_err)
               if (mdot_err.lt.mdot_tol) exit

            end do

            ! Integral of mdot
            call cfg%integrate(mdot,mdot_int)
            call cfg%integrate(mdotL,mdotL_int)
            call cfg%integrate(mdotG,mdotG_int)
            mdotL_int_err=abs(mdotL_int-mdot_int)
            mdotG_int_err=abs(mdotG_int-mdot_int)

         end block shift_mdot

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
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
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
            ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*(fs%div+(mdotG/fs%rho_g-mdotL/fs%rho_l)*vf%SD)/time%dt ! Evaporation mass source term is taken into account here
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
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call sc%get_max(VF=vf%VF)
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call mdotfile%write()
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,resSC,mdot,mdotL,mdotL_old,mdotG,mdotG_old,resmdotL,resmdotG)
   end subroutine simulation_final
   
   
end module simulation