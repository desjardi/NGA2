!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps,es
   !type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Problem definition
   real(WP) :: radius,q0
   real(WP) :: eps_g,eps_l
   
   !> EHD variables
   real(WP), dimension(:,:,:), allocatable :: Ex,Ey,Ez
   real(WP), dimension(:,:,:), allocatable :: Exi,Eyi,Ezi
   real(WP), dimension(:,:,:), allocatable :: phi
   real(WP), dimension(:,:,:), allocatable :: q
   
contains
   
   
   !> Function that defines a level set function for a drop
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,theta,amp,mode
      ! Create a spherical droplet
      !G=radius-sqrt(sum(xyz**2))
      ! Create a deformed droplet
      theta=atan2(xyz(2),xyz(1))
      amp=0.3_WP
      mode=2.0_WP
      G=radius*(1.0_WP+amp*0.5_WP*cos(mode*theta))-sqrt(sum(xyz**2))
   end function levelset_drop
   
   
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
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,remap
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=remap,name='VOF')
         !vf%cons_correct=.false.
         ! Initialize to a droplet
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
      end block create_and_initialize_vof
      
      
      ! Prepare EHD model
      create_ehd_solver: block
         use hypre_str_class, only: pcg_pfmg2
         integer :: i,j,k
         ! Read in charge density and permittivities
         call param_read('Charge density',q0)
         call param_read('Liquid permittivity',eps_l)
         call param_read('Gas permittivity',eps_g)
         ! Allocate EHD variables
         allocate(Ex (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Ex =0.0_WP
         allocate(Ey (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Ey =0.0_WP
         allocate(Ez (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Ez =0.0_WP
         allocate(Exi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Exi=0.0_WP
         allocate(Eyi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Eyi=0.0_WP
         allocate(Ezi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Ezi=0.0_WP
         allocate(q  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); q  =0.0_WP
         allocate(phi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); phi=0.0_WP
         ! Initialize charge density
         q=q0*vf%VF
         ! Prepare Poisson solver
         es=hypre_str(cfg=cfg,name='Electric',method=pcg_pfmg2,nst=7)
         es%maxlevel=10
         call param_read('Electric iteration',es%maxit)
         call param_read('Electric tolerance',es%rcvg)
         ! Set 7-pt stencil map for the electric Poisson solver
         es%stc(1,:)=[ 0, 0, 0]
         es%stc(2,:)=[+1, 0, 0]
         es%stc(3,:)=[-1, 0, 0]
         es%stc(4,:)=[ 0,+1, 0]
         es%stc(5,:)=[ 0,-1, 0]
         es%stc(6,:)=[ 0, 0,+1]
         es%stc(7,:)=[ 0, 0,-1]
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  ! Set operator to div(grad(.))
                  es%opr(1,i,j,k)=-cfg%dxi(i)*cfg%dxmi(i+1)&
                  &               -cfg%dxi(i)*cfg%dxmi(i  )&
                  &               -cfg%dyi(j)*cfg%dymi(j+1)&
                  &               -cfg%dyi(j)*cfg%dymi(j  )&
                  &               -cfg%dzi(k)*cfg%dzmi(k+1)&
                  &               -cfg%dzi(k)*cfg%dzmi(k  )
                  es%opr(2,i,j,k)=+cfg%dxi(i)*cfg%dxmi(i+1)
                  es%opr(3,i,j,k)=+cfg%dxi(i)*cfg%dxmi(i  )
                  es%opr(4,i,j,k)=+cfg%dyi(j)*cfg%dymi(j+1)
                  es%opr(5,i,j,k)=+cfg%dyi(j)*cfg%dymi(j  )
                  es%opr(6,i,j,k)=+cfg%dzi(k)*cfg%dzmi(k+1)
                  es%opr(7,i,j,k)=+cfg%dzi(k)*cfg%dzmi(k  )
                  ! Scale it by the cell volume
                  es%opr(:,i,j,k)=-es%opr(:,i,j,k)*cfg%vol(i,j,k)
               end do
            end do
         end do
         call es%init()
         call es%setup()
      end block create_ehd_solver
      
      
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
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='EHDdrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('q',q)
         call ens_out%add_scalar('phi',phi)
         call ens_out%add_vector('E',Exi,Eyi,Ezi)
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
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
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
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
         ! Solve EHD equations
         solve_ehd: block
            integer :: i,j,k
            real(WP) :: liq_vol,gas_vol,tot_vol
            real(WP) :: eps_xm,eps_xp,eps_ym,eps_yp,eps_zm,eps_zp
            ! Update charge density
            q=q0*vf%VF
            ! Recompute electric Poisson operator
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     ! Get permittivity at x- face
                     liq_vol=sum(vf%Lvol(0,:,:,i,j,k))+sum(vf%Lvol(1,:,:,i-1,j,k))
                     gas_vol=sum(vf%Gvol(0,:,:,i,j,k))+sum(vf%Gvol(1,:,:,i-1,j,k))
                     tot_vol=gas_vol+liq_vol
                     eps_xm=0.0_WP; if (tot_vol.gt.0.0_WP) eps_xm=eps_g*eps_l/(eps_l*gas_vol/tot_vol+eps_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     ! Get permittivity at x+ face
                     liq_vol=sum(vf%Lvol(0,:,:,i+1,j,k))+sum(vf%Lvol(1,:,:,i,j,k))
                     gas_vol=sum(vf%Gvol(0,:,:,i+1,j,k))+sum(vf%Gvol(1,:,:,i,j,k))
                     tot_vol=gas_vol+liq_vol
                     eps_xp=0.0_WP; if (tot_vol.gt.0.0_WP) eps_xp=eps_g*eps_l/(eps_l*gas_vol/tot_vol+eps_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     ! Get permittivity at y- face
                     liq_vol=sum(vf%Lvol(:,0,:,i,j,k))+sum(vf%Lvol(:,1,:,i,j-1,k))
                     gas_vol=sum(vf%Gvol(:,0,:,i,j,k))+sum(vf%Gvol(:,1,:,i,j-1,k))
                     tot_vol=gas_vol+liq_vol
                     eps_ym=0.0_WP; if (tot_vol.gt.0.0_WP) eps_ym=eps_g*eps_l/(eps_l*gas_vol/tot_vol+eps_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     ! Get permittivity at y+ face
                     liq_vol=sum(vf%Lvol(:,0,:,i,j+1,k))+sum(vf%Lvol(:,1,:,i,j,k))
                     gas_vol=sum(vf%Gvol(:,0,:,i,j+1,k))+sum(vf%Gvol(:,1,:,i,j,k))
                     tot_vol=gas_vol+liq_vol
                     eps_yp=0.0_WP; if (tot_vol.gt.0.0_WP) eps_yp=eps_g*eps_l/(eps_l*gas_vol/tot_vol+eps_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     ! Get permittivity at z- face
                     liq_vol=sum(vf%Lvol(:,:,0,i,j,k))+sum(vf%Lvol(:,:,1,i,j,k-1))
                     gas_vol=sum(vf%Gvol(:,:,0,i,j,k))+sum(vf%Gvol(:,:,1,i,j,k-1))
                     tot_vol=gas_vol+liq_vol
                     eps_zm=0.0_WP; if (tot_vol.gt.0.0_WP) eps_zm=eps_g*eps_l/(eps_l*gas_vol/tot_vol+eps_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     ! Get permittivity at z+ face
                     liq_vol=sum(vf%Lvol(:,:,0,i,j,k+1))+sum(vf%Lvol(:,:,1,i,j,k))
                     gas_vol=sum(vf%Gvol(:,:,0,i,j,k+1))+sum(vf%Gvol(:,:,1,i,j,k))
                     tot_vol=gas_vol+liq_vol
                     eps_zp=0.0_WP; if (tot_vol.gt.0.0_WP) eps_zp=eps_g*eps_l/(eps_l*gas_vol/tot_vol+eps_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     ! Set operator to div(eps*grad(.))
                     es%opr(1,i,j,k)=-cfg%dxi(i)*cfg%dxmi(i+1)*eps_xp&
                     &               -cfg%dxi(i)*cfg%dxmi(i  )*eps_xm&
                     &               -cfg%dyi(j)*cfg%dymi(j+1)*eps_yp&
                     &               -cfg%dyi(j)*cfg%dymi(j  )*eps_ym&
                     &               -cfg%dzi(k)*cfg%dzmi(k+1)*eps_zp&
                     &               -cfg%dzi(k)*cfg%dzmi(k  )*eps_zm
                     es%opr(2,i,j,k)=+cfg%dxi(i)*cfg%dxmi(i+1)*eps_xp
                     es%opr(3,i,j,k)=+cfg%dxi(i)*cfg%dxmi(i  )*eps_xm
                     es%opr(4,i,j,k)=+cfg%dyi(j)*cfg%dymi(j+1)*eps_yp
                     es%opr(5,i,j,k)=+cfg%dyi(j)*cfg%dymi(j  )*eps_ym
                     es%opr(6,i,j,k)=+cfg%dzi(k)*cfg%dzmi(k+1)*eps_zp
                     es%opr(7,i,j,k)=+cfg%dzi(k)*cfg%dzmi(k  )*eps_zm
                     ! Scale it by the cell volume
                     es%opr(:,i,j,k)=-es%opr(:,i,j,k)*cfg%vol(i,j,k)
                     ! Apply Dirichlet boundary conditions here - grounded box
                     if (cfg%VF(i+1,j,k).eq.0.0_WP) es%opr(2,i,j,k)=0.0_WP
                     if (cfg%VF(i-1,j,k).eq.0.0_WP) es%opr(3,i,j,k)=0.0_WP
                     if (cfg%VF(i,j+1,k).eq.0.0_WP) es%opr(4,i,j,k)=0.0_WP
                     if (cfg%VF(i,j-1,k).eq.0.0_WP) es%opr(5,i,j,k)=0.0_WP
                     if (cfg%VF(i,j,k+1).eq.0.0_WP) es%opr(6,i,j,k)=0.0_WP
                     if (cfg%VF(i,j,k-1).eq.0.0_WP) es%opr(7,i,j,k)=0.0_WP
                  end do
               end do
            end do
            ! Update electric potential
            es%rhs=fs%cfg%vol*q
            es%sol=0.0_WP
            call es%solve()
            phi=es%sol
            call cfg%sync(phi)
            ! Update electric field
            do k=cfg%kmino_+1,cfg%kmaxo_
               do j=cfg%jmino_+1,cfg%jmaxo_
                  do i=cfg%imino_+1,cfg%imaxo_
                     Ex(i,j,k)=-(phi(i,j,k)-phi(i-1,j,k))*cfg%dxmi(i)
                     Ey(i,j,k)=-(phi(i,j,k)-phi(i,j-1,k))*cfg%dymi(j)
                     Ez(i,j,k)=-(phi(i,j,k)-phi(i,j,k-1))*cfg%dzmi(k)
                  end do
               end do
            end do
            call cfg%sync(Ex)
            call cfg%sync(Ey)
            call cfg%sync(Ez)
            ! Update cell-centered electric field
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_-1
                     Exi(i,j,k)=0.5_WP*sum(Ex(i:i+1,j,k))
                  end do
               end do
            end do
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_-1
                  do i=cfg%imino_,cfg%imaxo_
                     Eyi(i,j,k)=0.5_WP*sum(Ey(i,j:j+1,k))
                  end do
               end do
            end do
            do k=cfg%kmino_,cfg%kmaxo_-1
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     Ezi(i,j,k)=0.5_WP*sum(Ez(i,j,k:k+1))
                  end do
               end do
            end do
            if (.not.cfg%xper.and.cfg%iproc.eq.cfg%npx) Exi(cfg%imaxo,:,:)=Ex(cfg%imaxo,:,:)
            if (.not.cfg%yper.and.cfg%jproc.eq.cfg%npy) Eyi(:,cfg%jmaxo,:)=Ey(:,cfg%jmaxo,:)
            if (.not.cfg%zper.and.cfg%kproc.eq.cfg%npz) Ezi(:,:,cfg%kmaxo)=Ez(:,:,cfg%kmaxo)
            call cfg%sync(Exi)
            call cfg%sync(Eyi)
            call cfg%sync(Ezi)
         end block solve_ehd

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
            
            ! Add Coulomb force
            add_coulomb: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+(fs%rho_U(i,j,k)-fs%rho_g)/(fs%rho_l-fs%rho_g)*q0*Ex(i,j,k)
                        if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+(fs%rho_V(i,j,k)-fs%rho_g)/(fs%rho_l-fs%rho_g)*q0*Ey(i,j,k)
                        if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+(fs%rho_W(i,j,k)-fs%rho_g)/(fs%rho_l-fs%rho_g)*q0*Ez(i,j,k)
                     end do
                  end do
               end do
            end block add_coulomb
            
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
      
   end subroutine simulation_final
   
   
end module simulation