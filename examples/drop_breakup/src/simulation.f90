!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use ccl_class,         only: ccl
   use lpt_class,         only: lpt
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Two-phase incompressible flow solver, VF solver with CCL, and corresponding time tracker and sgs model
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs
   type(ccl),         public :: cc
   type(lpt),         public :: lp
   
   !> Provide two datafiles and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   character(len=str_medium) :: irl_file
   character(len=str_medium) :: lpt_file
   logical :: restarted
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,sprayfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   
   !> Problem definition
   real(WP), dimension(3) :: center,radii
   
   !> Transfer model parameters
   real(WP) :: filmthickness_over_dx  =5.0e-1_WP
   real(WP) :: min_filmthickness      =1.0e-3_WP
   real(WP) :: diam_over_filmthickness=1.0e+1_WP
   real(WP) :: max_eccentricity       =5.0e-1_WP
   real(WP) :: d_threshold            =1.0e-1_WP
   
   !> SGS surface tension model
   real(WP), dimension(:,:,:), allocatable :: sgsSTx,sgsSTy,sgsSTz
   
contains
   
   
   !> Function that localizes the left (x-) of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator
   
   
   !> Function that localizes the right (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   
   !> Function that defines a level set function for a spherical droplet
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP-sqrt(((xyz(1)-center(1))/radii(1))**2+((xyz(2)-center(2))/radii(2))**2+((xyz(3)-center(3))/radii(3))**2)
   end function levelset_drop
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Handle restart/saves here
      restart_and_save: block
         character(len=str_medium) :: dir_restart
         ! CAREFUL - WE NEED TO CREATE THE TIMETRACKER BEFORE THE EVENT !
         time=timetracker(cfg%amRoot,name='drop break-up')
         ! Create event for saving restart files
         save_evt=event(time,'Restart output')
         call param_read('Restart output period',save_evt%tper)
         ! Check if we are restarting
         call param_read(tag='Restart from',val=dir_restart,short='r',default='')
         restarted=.false.; if (len_trim(dir_restart).gt.0) restarted=.true.
         if (restarted) then
            ! If we are, read the name of the directory
            call param_read('Restart from',dir_restart,'r')
            ! Read the datafile and the name of the IRL file to read later
            df=datafile(pg=cfg,fdata=trim(adjustl(dir_restart))//'/'//'data')
            irl_file=trim(adjustl(dir_restart))//'/'//'data.irl'
            lpt_file=trim(adjustl(dir_restart))//'/'//'data.lpt'
         else
            ! If we are not restarting, we will still need datafiles for saving restart files
            df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=10)
            df%valname(1)='t'
            df%valname(2)='dt'
            df%varname(1)='U'
            df%varname(2)='V'
            df%varname(3)='W'
            df%varname(4)='P'
            df%varname(5)='Pjx'
            df%varname(6)='Pjy'
            df%varname(7)='Pjz'
            df%varname(8)='LM'
            df%varname(9)='MM'
            df%varname(10)='VF'
         end if
      end block restart_and_save
      
      
      ! Allocate work arrays for cfg
      allocate_work_arrays: block
         allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
      end block allocate_work_arrays
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         !time=timetracker(cfg%amRoot,name='drop break-up')   !< This is moved up for restarts!
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
         ! Handle restart
         if (restarted) then
            call df%pullval(name='t' ,val=time%t )
            call df%pullval(name='dt',val=time%dt)
            time%told=time%t-time%dt
         end if
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: VFlo,VFhi,lvira,r2p,art,swartz
         integer :: i,j,k,si,sj,sk,n
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with LVIRA
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Create a VOF solver with R2P
         !vf=vfs(cfg=cfg,reconstruction_method=r2p,name='VOF')
         !vf%VFflot =1.0e-4_WP !< Enables flotsam removal
         !vf%VFsheet=1.0e-2_WP !< Enables sheet removal
         ! Create a VOF solver with ART
         !vf=vfs(cfg=cfg,reconstruction_method=art,name='VOF')
         ! Initialize to droplet
         call param_read('Droplet center',center)
         call param_read('Droplet radii',radii)
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
         ! Handle restart - using IRL data - tested but not working
         !if (restarted) then
         !   ! Get the IRL interface
         !   call vf%read_interface(filename=trim(irl_file))
         !   ! Reset moments to guarantee compatibility with interface reconstruction
         !   call vf%reset_volume_moments()
         !end if
         ! Handle restart - using VF data
         if (restarted) call df%pullvar(name='VF',var=vf%VF)
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
      
      
      ! Create a two-phase flow solver with bconds
      create_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,neumann
         use hypre_str_class,  only: pcg_pfmg2
         ! Create a two-phase flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity'   ,fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density'   ,fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Define boundary conditions: gas inflow on the left and outflow on the right
         call fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools,  only: pi
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         real(WP) :: Uin
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Handle restart
         if (restarted) then
            call df%pullvar(name='U'  ,var=fs%U  )
            call df%pullvar(name='V'  ,var=fs%V  )
            call df%pullvar(name='W'  ,var=fs%W  )
            call df%pullvar(name='P'  ,var=fs%P  )
            call df%pullvar(name='Pjx',var=fs%Pjx)
            call df%pullvar(name='Pjy',var=fs%Pjy)
            call df%pullvar(name='Pjz',var=fs%Pjz)
         end if
         ! Apply Dirichlet at inflow
         call param_read('Gas velocity',Uin)
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=Uin
         end do
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
      
      
      ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         cc=ccl(cfg=cfg,name='CCL')
         cc%max_interface_planes=2
         cc%VFlo=VFlo
         cc%dot_threshold=-0.5_WP
         cc%thickness_cutoff=filmthickness_over_dx
         ! Perform CCL step
         call cc%build_lists(VF=vf%VF,poly=vf%interface_polygon,U=fs%U,V=fs%V,W=fs%W)
         call cc%film_classify(Lbary=vf%Lbary,Gbary=vf%Gbary)
         call cc%deallocate_lists()
      end block create_and_initialize_ccl
      
      
      ! Create a Lagrangian spray tracker
      create_lpt: block
         ! Create the solver
         lp=lpt(cfg=cfg,name='spray')
         ! Get particle density from the flow solver
         lp%rho=fs%rho_l
         ! Handle restarts
         if (restarted) call lp%read(filename=trim(lpt_file))
      end block create_lpt
      
      
      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         ! Handle restart
         if (restarted) then
            call df%pullvar(name='LM',var=sgs%LM)
            call df%pullvar(name='MM',var=sgs%MM)
         end if
      end block create_sgs
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='nplane'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         smesh%var(1,:)=1.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         ! Include an extra variable for droplet diameter
         pmesh=partmesh(nvar=1,nvec=0,name='lpt')
         pmesh%varname(1)='diameter'
         ! Transfer particles to pmesh
         call lp%update_partmesh(pmesh)
         ! Also populate diameter variable
         do i=1,lp%np_
            pmesh%var(1,i)=lp%p(i)%d
         end do
      end block create_pmesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg,'drop_breakup')
         ! Create event for Ensight output
         ens_evt=event(time,'Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curv',vf%curv)
         call ens_out%add_scalar('filmThickness',cc%film_thickness)
         call ens_out%add_surface('vofplic',smesh)
         call ens_out%add_particle('spray',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call lp%get_max()
         ! Create simulation monitor
         mfile=monitor(amroot=fs%cfg%amRoot,name='simulation')
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
         cflfile=monitor(amroot=fs%cfg%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(fs%CFLst,'ST CFL')
         call cflfile%write()
         ! Create a spray monitor
         sprayfile=monitor(amroot=lp%cfg%amRoot,name='spray')
         call sprayfile%add_column(time%n,'Timestep number')
         call sprayfile%add_column(time%t,'Time')
         call sprayfile%add_column(time%dt,'Timestep size')
         call sprayfile%add_column(lp%np,'Droplet number')
         call sprayfile%add_column(lp%Umin, 'Umin')
         call sprayfile%add_column(lp%Umax, 'Umax')
         call sprayfile%add_column(lp%Umean,'Umean')
         call sprayfile%add_column(lp%Vmin, 'Vmin')
         call sprayfile%add_column(lp%Vmax, 'Vmax')
         call sprayfile%add_column(lp%Vmean,'Vmean')
         call sprayfile%add_column(lp%Wmin, 'Wmin')
         call sprayfile%add_column(lp%Wmax, 'Wmax')
         call sprayfile%add_column(lp%Wmean,'Wmean')
         call sprayfile%add_column(lp%dmin, 'dmin')
         call sprayfile%add_column(lp%dmax, 'dmax')
         call sprayfile%add_column(lp%dmean,'dmean')
         call sprayfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration - the second solver is the main driver here
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Advance our spray
         resU=fs%rho_g; resV=fs%visc_g
         call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV)
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)
         
         ! Turbulence modeling - only work with gas properties here
         sgsmodel: block
            use sgsmodel_class, only: vreman
            integer :: i,j,k
            call fs%get_gradu(gradU)
            call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,Ui=Ui,Vi=Vi,Wi=Wi,gradU=gradU)
            where (sgs%visc.lt.-fs%visc_g)
               sgs%visc=-fs%visc_g
            end where
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
         end block sgsmodel
         
         ! ! SGS modeling of surface tension force
         ! STmodel_prep: block
         !    use mathtools, only: normalize
         !    integer :: i,j,k,ii,jj,kk
         !    real(WP) :: thickness,STmag,fvol,sgscurv,fSD,xtmp,ytmp,ztmp
         !    real(WP), dimension(3)   :: fbary
         !    real(WP), dimension(3,3) :: Imom
         !    real(WP), parameter :: ligam_threshold=1.3_WP
         !    real(WP), parameter :: thick_threshold=1.0_WP
         !    ! LAPACK Eigenvalues/eigenvectors
         !    integer :: info
         !    real(WP), dimension(3)     :: d
         !    integer , parameter        :: order=3
         !    integer , parameter        :: lwork=102 ! dsyev optimal length (nb+2)*order, where block size nb=32
         !    real(WP), dimension(lwork) :: work
         !    ! Zero out sgs term
         !    sgsSTx=0.0_WP; sgsSTy=0.0_WP; sgsSTz=0.0_WP
         !    ! Loop over the domain and find ligament/edge cells
         !    do k=vf%cfg%kmin_,vf%cfg%kmax_
         !       do j=vf%cfg%jmin_,vf%cfg%jmax_
         !          do i=vf%cfg%imin_,vf%cfg%imax_
         !             ! Work only if interface is present
         !             if (vf%SD(i,j,k).eq.0.0_WP) cycle
         !             ! Compute thickness
         !             fvol=0.0_WP; fSD=0.0_WP
         !             do kk=k-1,k+1
         !                do jj=j-1,j+1
         !                   do ii=i-1,i+1
         !                      fvol=fvol+vf%VF(ii,jj,kk)
         !                      fSD =fSD +vf%SD(ii,jj,kk)
         !                   end do
         !                end do
         !             end do
         !             thickness=2.0_WP*fvol/fSD
         !             ! Calculate barycenter
         !             fvol=0.0_WP; fbary=0.0_WP; Imom=0.0_WP
         !             do kk=k-1,k+1
         !                do jj=j-1,j+1
         !                   do ii=i-1,i+1
         !                      fvol =fvol +vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                      fbary=fbary+vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)*vf%Lbary(:,ii,jj,kk)
         !                   end do
         !                end do
         !             end do
         !             fbary=fbary/fvol
         !             ! Second pass to calculate moment of inertia
         !             do kk=k-1,k+1
         !                do jj=j-1,j+1
         !                   do ii=i-1,i+1
         !                      ! Shifted local barycenter
         !                      xtmp=vf%Lbary(1,ii,jj,kk)-fbary(1)
         !                      ytmp=vf%Lbary(2,ii,jj,kk)-fbary(2)
         !                      ztmp=vf%Lbary(3,ii,jj,kk)-fbary(3)
         !                      ! Moment of inertia
         !                      Imom(1,1)=Imom(1,1)+(ytmp**2+ztmp**2)*vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                      Imom(2,2)=Imom(2,2)+(xtmp**2+ztmp**2)*vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                      Imom(3,3)=Imom(3,3)+(xtmp**2+ytmp**2)*vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                      Imom(1,2)=Imom(1,2)-xtmp*ytmp*vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                      Imom(1,3)=Imom(1,3)-xtmp*ztmp*vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                      Imom(2,3)=Imom(2,3)-ytmp*ztmp*vf%cfg%vol(ii,jj,kk)*vf%VF(ii,jj,kk)
         !                   end do
         !                end do
         !             end do
         !             ! On exit, Imom contains eigenvectors, and d contains eigenvalues in ascending order
         !             call dsyev('V','U',order,Imom,order,d,work,lwork,info); d=max(0.0_WP,d)
         !             ! Only apply model for ligaments/edges that are thin enough
         !             if (d(2).gt.ligam_threshold*sqrt(d(1)*d(3)).and.thickness.lt.thick_threshold*vf%cfg%meshsize(i,j,k)) then
         !                ! First compute ST force direction
         !                fbary=normalize(fbary-vf%Lbary(:,i,j,k))
         !                ! Estimate SGS curvature
         !                sgscurv=min(max(1.0_WP/thickness-abs(vf%curv(i,j,k)),0.0_WP),1.0_WP/(fs%CFLst*vf%cfg%meshsize(i,j,k)))
         !                ! Compute force magnitude
         !                STmag=fs%sigma*vf%SD(i,j,k)*sgscurv
         !                ! Apply it
         !                sgsSTx(i,j,k)=STmag*fbary(1)
         !                sgsSTy(i,j,k)=STmag*fbary(2)
         !                sgsSTz(i,j,k)=STmag*fbary(3)
         !             end if
         !          end do
         !       end do
         !    end do
         !    ! Sync it
         !    call vf%cfg%sync(sgsSTx)
         !    call vf%cfg%sync(sgsSTy)
         !    call vf%cfg%sync(sgsSTz)
         ! end block STmodel_prep
         
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
            
            ! ! Add sgs ST model
            ! STmodel_add: block
            !    integer :: i,j,k
            !    do k=vf%cfg%kmin_,vf%cfg%kmax_
            !       do j=vf%cfg%jmin_,vf%cfg%jmax_
            !          do i=vf%cfg%imin_,vf%cfg%imax_
            !             if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%itpi_x(:,i,j,k)*sgsSTx(i-1:i,j,k))
            !             if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%itpi_y(:,i,j,k)*sgsSTy(i,j-1:j,k))
            !             if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%itpi_z(:,i,j,k)*sgsSTz(i,j,k-1:k))
            !          end do
            !       end do
            !    end do
            ! end block STmodel_add
            
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
         
         ! Perform volume-fraction-to-droplet transfer
         call transfer_vf_to_drops()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Also populate nplane variable
               smesh%var(1,:)=1.0_WP
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            ! Update partmesh object
            update_pmesh: block
               integer :: i
               ! Transfer particles to pmesh
               call lp%update_partmesh(pmesh)
               ! Also populate diameter variable
               do i=1,lp%np_
                  pmesh%var(1,i)=lp%p(i)%d
               end do
            end block update_pmesh
            ! Perform ensight output
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call lp%get_max()
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         call sprayfile%write()
         
         ! After we're done clip all VOF at the exit area
         vf_side_clipping: block
            integer :: i,j,k
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     if (i.ge.vf%cfg%imax-5) vf%VF(i,j,k)=0.0_WP
                  end do
               end do
            end do
         end block vf_side_clipping
         
         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
               character(len=str_medium) :: dirname,timestamp
               ! Prefix for files
               dirname='restart_'; write(timestamp,'(es12.5)') time%t
               ! Prepare a new directory
               if (fs%cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
               ! Populate df and write it
               call df%pushval(name=  't',val=time%t )
               call df%pushval(name= 'dt',val=time%dt)
               call df%pushvar(name=  'U',var=fs%U   )
               call df%pushvar(name=  'V',var=fs%V   )
               call df%pushvar(name=  'W',var=fs%W   )
               call df%pushvar(name=  'P',var=fs%P   )
               call df%pushvar(name='Pjx',var=fs%Pjx )
               call df%pushvar(name='Pjy',var=fs%Pjy )
               call df%pushvar(name='Pjz',var=fs%Pjz )
               call df%pushvar(name= 'LM',var=sgs%LM )
               call df%pushvar(name= 'MM',var=sgs%MM )
               call df%pushvar(name= 'VF',var=vf%VF  )
               call df%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data')
               ! Also output IRL interface
               call vf%write_interface(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.irl')
               ! Also output particles
               call lp%write(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.lpt')
            end block save_restart
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Transfer vf to drops
   subroutine transfer_vf_to_drops()
      implicit none
      
      ! Perform a first pass with simplest CCL
      call cc%build_lists(VF=vf%VF,U=fs%U,V=fs%V,W=fs%W)
      
      ! Loop through identified detached structs and remove those that are spherical enough
      remove_struct: block
         use mathtools, only: pi
         integer :: m,n,l,i,j,k,np
         real(WP) :: lmin,lmax,eccentricity,diam
         
         ! Loops over film segments contained locally
         do m=1,cc%n_meta_struct
            
            ! Test if sphericity is compatible with transfer
            lmin=cc%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=cc%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=cc%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
            if (eccentricity.gt.max_eccentricity) cycle
            
            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if (diam.eq.0.0_WP.or.diam.gt.d_threshold) cycle
            
            ! Create drop from available liquid volume - only one root does that
            if (cc%cfg%amRoot) then
               ! Make room for new drop
               np=lp%np_+1; call lp%resize(np)
               ! Add the drop
               lp%p(np)%id  =int(0,8)                                                                                 !< Give id (maybe based on break-up model?)
               lp%p(np)%dt  =0.0_WP                                                                                   !< Let the drop find it own integration time
               lp%p(np)%Acol=0.0_WP                                                                                   !< Give zero collision force (axial)
               lp%p(np)%Tcol=0.0_WP                                                                                   !< Give zero collision force (tangential)
               lp%p(np)%d   =diam                                                                                     !< Assign diameter to account for full volume
               lp%p(np)%pos =[cc%meta_structures_list(m)%x,cc%meta_structures_list(m)%y,cc%meta_structures_list(m)%z] !< Place the drop at the liquid barycenter
               lp%p(np)%vel =[cc%meta_structures_list(m)%u,cc%meta_structures_list(m)%v,cc%meta_structures_list(m)%w] !< Assign mean structure velocity as drop velocity
               lp%p(np)%ind =lp%cfg%get_ijk_global(lp%p(np)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])                !< Place the drop in the proper cell for the lp%cfg
               lp%p(np)%flag=0                                                                                        !< Activate it
               ! Increment particle counter
               lp%np_=np
            end if
            
            ! Find local structs with matching id
            do n=cc%sync_offset+1,cc%sync_offset+cc%n_struct
               if (cc%struct_list(cc%struct_map_(n))%parent.ne.cc%meta_structures_list(m)%id) cycle
               ! Remove liquid in meta-structure cells
               do l=1,cc%struct_list(cc%struct_map_(n))%nnode ! Loops over cells within local
                  i=cc%struct_list(cc%struct_map_(n))%node(1,l)
                  j=cc%struct_list(cc%struct_map_(n))%node(2,l)
                  k=cc%struct_list(cc%struct_map_(n))%node(3,l)
                  ! Remove liquid in that cell
                  vf%VF(i,j,k)=0.0_WP
               end do
            end do
            
         end do
         
      end block remove_struct
      
      ! Sync VF and clean up IRL and band
      call vf%cfg%sync(vf%VF)
      call vf%clean_irl_and_band()
      
      ! Clean up CCL
      call cc%deallocate_lists()
      
      ! Perform more detailed CCL in a second pass
      cc%max_interface_planes=2
      call cc%build_lists(VF=vf%VF,poly=vf%interface_polygon,U=fs%U,V=fs%V,W=fs%W)
      call cc%get_min_thickness()
      call cc%sort_by_thickness()
      
      ! Loop through identified films and remove those that are thin enough
      remove_film: block
         use mathtools, only: pi
         integer :: m,n,i,j,k,np,ip,np_old
         real(WP) :: Vt,Vl,Hl,Vd
         
         ! Loops over film segments contained locally
         do m=cc%film_sync_offset+1,cc%film_sync_offset+cc%n_film
            
            ! Skip non-liquid films
            if (cc%film_list(cc%film_map_(m))%phase.ne.1) cycle
            
            ! Skip films that are still thick enough
            if (cc%film_list(cc%film_map_(m))%min_thickness.gt.min_filmthickness) cycle
            
            ! We are still here: transfer the film to drops
            Vt=0.0_WP      ! Transferred volume
            Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
            np_old=lp%np_  ! Remember old number of particles
            do n=1,cc%film_list(cc%film_map_(m))%nnode ! Loops over cells within local film segment
               i=cc%film_list(cc%film_map_(m))%node(1,n)
               j=cc%film_list(cc%film_map_(m))%node(2,n)
               k=cc%film_list(cc%film_map_(m))%node(3,n)
               ! Increment liquid volume to remove
               Vl=Vl+vf%VF(i,j,k)*vf%cfg%vol(i,j,k)
               ! Estimate drop size based on local film thickness in current cell
               Hl=max(cc%film_thickness(i,j,k),min_filmthickness)
               Vd=pi/6.0_WP*(diam_over_filmthickness*Hl)**3
               ! Create drops from available liquid volume
               do while (Vl-Vd.gt.0.0_WP)
                  ! Make room for new drop
                  np=lp%np_+1; call lp%resize(np)
                  ! Add the drop
                  lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                  lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  lp%p(np)%Acol=0.0_WP                                     !< Collision acceleration                                              !< Give zero collision force (axial)
                  lp%p(np)%Tcol=0.0_WP                                     !< Collision torque                                              !< Give zero collision force (tangential)
                  lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                  lp%p(np)%pos =vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
                  lp%p(np)%vel =fs%cfg%get_velocity(pos=lp%p(np)%pos,i0=i,j0=j,k0=k,U=fs%U,V=fs%V,W=fs%W) !< Interpolate local cell velocity as drop velocity
                  lp%p(np)%ind =lp%cfg%get_ijk_global(lp%p(np)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                  lp%p(np)%flag=0                                          !< Activate it
                  ! Increment particle counter
                  lp%np_=np
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
               end do
               ! Remove liquid in that cell
               vf%VF(i,j,k)=0.0_WP
            end do
            
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
               ! Add one last drop for remaining liquid volume
               np=lp%np_+1; call lp%resize(np)
               ! Add the drop
               lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
               lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               lp%p(np)%Acol=0.0_WP                                     !< Collision acceleration
               lp%p(np)%Tcol=0.0_WP                                     !< Collision torque
               lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
               lp%p(np)%pos =vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
               lp%p(np)%vel =fs%cfg%get_velocity(pos=lp%p(np)%pos,i0=i,j0=j,k0=k,U=fs%U,V=fs%V,W=fs%W) !< Interpolate local cell velocity as drop velocity
               lp%p(np)%ind =lp%cfg%get_ijk_global(lp%p(np)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
               lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               lp%np_=np
            else ! Some particles were created, make them all larger
               do ip=np_old+1,lp%np_
                  lp%p(ip)%d=lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
               end do
            end if
         end do
         
      end block remove_film
      
      ! Sync VF and clean up IRL and band
      call vf%cfg%sync(vf%VF)
      call vf%clean_irl_and_band()
      
      ! Clean up CCL
      call cc%deallocate_lists()
      
      ! Resync the spray
      call lp%sync()
      
   end subroutine transfer_vf_to_drops
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)
      
   end subroutine simulation_final
   
   
end module simulation
