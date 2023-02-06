!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg1,cfg2,cfg3
   use geometry,          only: xinj_dist,inj_norm_diam,rli0,rlo0,rgi0,rlo0_spray,rgi0_spray
   use incomp_class,      only: incomp
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use ccl_class,         only: ccl
   use lpt_class,         only: lpt
   use coupler_class,     only: coupler
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
   
   !> Single-phase incompressible flow solver, and corresponding time tracker and sgs model
   type(incomp),      public :: fs1
   type(timetracker), public :: time1
   type(sgsmodel),    public :: sgs1
   
   ! Couplers between 1 and 2
   type(coupler),     public :: cpl12x,cpl12y,cpl12z
   
   !> Two-phase incompressible flow solver, VF solver with CCL, and corresponding time tracker and sgs model
   type(tpns),        public :: fs2
   type(vfs),         public :: vf2
   type(timetracker), public :: time2
   type(sgsmodel),    public :: sgs2
   type(ccl),         public :: cc2
   
   ! Couplers between 2 and 3
   type(coupler),     public :: cpl23x,cpl23y,cpl23z
   
   !> Single-phase incompressible flow solver with lpt, and corresponding time tracker and sgs model
   type(incomp),      public :: fs3
   type(lpt),         public :: lp3
   type(timetracker), public :: time3
   type(sgsmodel),    public :: sgs3
   
   !> Provide three datafiles and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df1,df2,df3
   character(len=str_medium) :: irl_file
   character(len=str_medium) :: lpt_file
   logical :: restarted
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh2
   type(partmesh) :: pmesh3
   type(ensight)  :: ens_out1
   type(ensight)  :: ens_out2
   type(ensight)  :: ens_out3
   type(event)    :: ens_evt1
   type(event)    :: ens_evt2
   type(event)    :: ens_evt3
   
   !> Simulation monitor file
   type(monitor) :: mfile1,cflfile1
   type(monitor) :: mfile2,cflfile2
   type(monitor) :: mfile3,cflfile3,sprayfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU1,resV1,resW1
   real(WP), dimension(:,:,:), allocatable :: Ui1,Vi1,Wi1
   real(WP), dimension(:,:,:,:), allocatable :: SR1
   real(WP), dimension(:,:,:), allocatable :: resU2,resV2,resW2
   real(WP), dimension(:,:,:), allocatable :: Ui2,Vi2,Wi2
   real(WP), dimension(:,:,:,:), allocatable :: SR2
   real(WP), dimension(:,:,:), allocatable :: resU3,resV3,resW3
   real(WP), dimension(:,:,:), allocatable :: Ui3,Vi3,Wi3
   real(WP), dimension(:,:,:,:), allocatable :: SR3
   real(WP), dimension(:,:,:), allocatable :: U2on3,V2on3,W2on3
   
   !> Transfer model parameters
   real(WP) :: filmthickness_over_dx  =5.0e-1_WP
   real(WP) :: min_filmthickness      =1.0e-6_WP
   real(WP) :: diam_over_filmthickness=1.0e+1_WP
   real(WP) :: max_eccentricity       =5.0e-1_WP
   real(WP) :: d_threshold            =1.0e-3_WP
   
contains
   
   
   !> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary
   
   !> Function that localizes injector at -y
   function inj1(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) ::hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-z plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[xinj_dist,0.0_WP,0.0_WP])
      rise=abs(pg%ym(j))
      run =sqrt(hyp**2-rise**2)
      if (run.le.inj_norm_diam/2.0_WP.and.j.eq.pg%jmin) isIn=.true.
   end function inj1

   !> Function that localizes injector at +y
   function inj2(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) ::hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-z plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[xinj_dist,0.0_WP,0.0_WP])
      rise=abs(pg%ym(j))
      run =sqrt(hyp**2-rise**2)
      if (run.le.inj_norm_diam/2.0_WP.and.j.eq.pg%jmax+1) isIn=.true.
   end function inj2

   !> Function that localizes injector at -z
   function inj3(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) ::hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-y plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[xinj_dist,0.0_WP,0.0_WP])
      rise=abs(pg%zm(k))
      run =sqrt(hyp**2-rise**2)
      if (run.le.inj_norm_diam/2.0_WP.and.k.eq.pg%kmin) isIn=.true.
   end function inj3

   !> Function that localizes injector at +z
   function inj4(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) ::hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-y plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[xinj_dist,0.0_WP,0.0_WP] )
      rise=abs(pg%zm(k))
      run =sqrt(hyp**2-rise**2)
      if (run.le.inj_norm_diam/2.0_WP.and.k.eq.pg%kmax+1) isIn=.true.
   end function inj4
   
   
   !> Function that localizes liquid stream at -x
   function liq_inj(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: rad
      isIn=.false.
      rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (rad.lt.rli0.and.i.eq.pg%imin) isIn=.true.
   end function liq_inj
   
   
   !> Function that localizes gas stream at -x
   function gas_inj(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: rad
      isIn=.false.
      rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (rad.ge.rlo0.and.rad.lt.rgi0.and.i.eq.pg%imin) isIn=.true.
   end function gas_inj
   
   !> Function that localizes gas stream at -x
   function gas_inj_spray(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: rad
      isIn=.false.
      rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (rad.ge.rlo0_spray.and.rad.lt.rgi0_spray.and.i.eq.pg%imin) isIn=.true.
   end function gas_inj_spray
   
   
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
   
   
   !> Function that localizes the top (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator
   
   
   !> Function that localizes the bottom (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Handle restart/saves here
      restart_and_save: block
         character(len=str_medium) :: dir_restart
         ! CAREFUL - WE NEED TO CREATE THE TIMETRACKER BEFORE THE EVENT !
         time2=timetracker(cfg2%amRoot,name='nozzle_exterior')
         ! Create event for saving restart files
         save_evt=event(time2,'Restart output')
         call param_read('Restart output period',save_evt%tper)
         ! Check if we are restarting
         call param_read(tag='Restart from',val=dir_restart,short='r',default='')
         restarted=.false.; if (len_trim(dir_restart).gt.0) restarted=.true.
         if (restarted) then
            ! If we are, read the name of the directory
            call param_read('Restart from',dir_restart,'r')
            ! Read the two datafiles and the name of the IRL file to read later
            df1=datafile(pg=cfg1,fdata=trim(adjustl(dir_restart))//'/'//'data.1')
            df2=datafile(pg=cfg2,fdata=trim(adjustl(dir_restart))//'/'//'data.2')
            df3=datafile(pg=cfg3,fdata=trim(adjustl(dir_restart))//'/'//'data.3')
            irl_file=trim(adjustl(dir_restart))//'/'//'data.irl'
            lpt_file=trim(adjustl(dir_restart))//'/'//'data.lpt'
         else
            ! If we are not restarting, we will still need datafiles for saving restart files
            df1=datafile(pg=cfg1,filename=trim(cfg1%name),nval=2,nvar=9)
            df1%valname(1)='t'
            df1%valname(2)='dt'
            df1%varname(1)='U'
            df1%varname(2)='V'
            df1%varname(3)='W'
            df1%varname(4)='P'
            df1%varname(5)='Uold'
            df1%varname(6)='Vold'
            df1%varname(7)='Wold'
            df1%varname(8)='LM'
            df1%varname(9)='MM'
            df2=datafile(pg=cfg2,filename=trim(cfg2%name),nval=2,nvar=10)
            df2%valname(1)='t'
            df2%valname(2)='dt'
            df2%varname(1)='U'
            df2%varname(2)='V'
            df2%varname(3)='W'
            df2%varname(4)='P'
            df2%varname(5)='Pjx'
            df2%varname(6)='Pjy'
            df2%varname(7)='Pjz'
            df2%varname(8)='LM'
            df2%varname(9)='MM'
            df2%varname(10)='VF'
            df3=datafile(pg=cfg3,filename=trim(cfg3%name),nval=2,nvar=6)
            df3%valname(1)='t'
            df3%valname(2)='dt'
            df3%varname(1)='U'
            df3%varname(2)='V'
            df3%varname(3)='W'
            df3%varname(4)='P'
            df3%varname(5)='LM'
            df3%varname(6)='MM'
         end if
      end block restart_and_save
      
      
      ! ###############################################
      ! ######## PREPARE SOLVER COUPLING HERE #########
      ! ###############################################
      coupler_prep: block
         use parallel, only: group
         ! Create nozzle-to-atomization couplers
         cpl12x=coupler(src_grp=group,dst_grp=group,name='nozzle_to_atom_x'); call cpl12x%set_src(cfg1,'x'); call cpl12x%set_dst(cfg2,'x'); call cpl12x%initialize()
         cpl12y=coupler(src_grp=group,dst_grp=group,name='nozzle_to_atom_y'); call cpl12y%set_src(cfg1,'y'); call cpl12y%set_dst(cfg2,'y'); call cpl12y%initialize()
         cpl12z=coupler(src_grp=group,dst_grp=group,name='nozzle_to_atom_z'); call cpl12z%set_src(cfg1,'z'); call cpl12z%set_dst(cfg2,'z'); call cpl12z%initialize()
         ! Create atomization-to-spray couplers
         cpl23x=coupler(src_grp=group,dst_grp=group,name='atom_to_spray__x'); call cpl23x%set_src(cfg2,'x'); call cpl23x%set_dst(cfg3,'x'); call cpl23x%initialize()
         cpl23y=coupler(src_grp=group,dst_grp=group,name='atom_to_spray__y'); call cpl23y%set_src(cfg2,'y'); call cpl23y%set_dst(cfg3,'y'); call cpl23y%initialize()
         cpl23z=coupler(src_grp=group,dst_grp=group,name='atom_to_spray__z'); call cpl23z%set_src(cfg2,'z'); call cpl23z%set_dst(cfg3,'z'); call cpl23z%initialize()
      end block coupler_prep
      
      
      ! *****************************************************************************
      ! **************************** INITIALIZE SOLVER 2 ****************************
      ! *****************************************************************************
      
      ! Allocate work arrays for cfg2
      allocate_work_arrays2: block
         allocate(resU2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         allocate(resV2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         allocate(resW2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         allocate(Ui2  (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         allocate(Vi2  (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         allocate(Wi2  (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
         allocate(SR2(6,cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_))
      end block allocate_work_arrays2
      
      
      ! Initialize time tracker
      initialize_timetracker2: block
         !time2=timetracker(cfg2%amRoot,name='nozzle_exterior')   !< This is moved up for restarts!
         call param_read('2 Max timestep size',time2%dtmax)
         call param_read('2 Max cfl number',time2%cflmax)
         time2%dt=time2%dtmax
         time2%itmax=2
         ! Handle restart
         if (restarted) then
            call df2%pullval(name='t' ,val=time2%t )
            call df2%pullval(name='dt',val=time2%dt)
            time2%told=time2%t-time2%dt
         end if
      end block initialize_timetracker2
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof2: block
         use vfs_class, only: lvira,r2p
         integer :: i,j,k
         real(WP) :: xloc,rad
         ! Create a VOF solver with LVIRA
         !vf2=vfs(cfg=cfg2,reconstruction_method=lvira,name='VOF')
         ! Create a VOF solver with R2P
         vf2=vfs(cfg=cfg2,reconstruction_method=r2p,name='VOF')
         !vf2%VFflot =1.0e-4_WP !< Enables flotsam removal
         !vf2%VFsheet=1.0e-2_WP !< Enables sheet removal
         ! Initialize to flat interface in liquid needle
         xloc=-0.001_WP !< Interface initially just inside the nozzle
         do k=vf2%cfg%kmino_,vf2%cfg%kmaxo_
            do j=vf2%cfg%jmino_,vf2%cfg%jmaxo_
               do i=vf2%cfg%imino_,vf2%cfg%imaxo_
                  rad=sqrt(vf2%cfg%ym(j)**2+vf2%cfg%zm(k)**2)
                  if (vf2%cfg%xm(i).lt.xloc.and.rad.le.rli0) then
                     vf2%VF(i,j,k)=1.0_WP
                  else
                     vf2%VF(i,j,k)=0.0_WP
                  end if
                  vf2%Lbary(:,i,j,k)=[vf2%cfg%xm(i),vf2%cfg%ym(j),vf2%cfg%zm(k)]
                  vf2%Gbary(:,i,j,k)=[vf2%cfg%xm(i),vf2%cfg%ym(j),vf2%cfg%zm(k)]
               end do
            end do
         end do
         ! Handle restart - using IRL data - tested but not working
         !if (restarted) then
         !   ! Get the IRL interface
         !   call vf2%read_interface(filename=trim(irl_file))
         !   ! Reset moments to guarantee compatibility with interface reconstruction
         !   call vf2%reset_volume_moments()
         !end if
         ! Handle restart - using VF data
         if (restarted) call df2%pullvar(name='VF',var=vf2%VF)
         ! Update the band
         call vf2%update_band()
         ! Perform interface reconstruction from VOF field
         call vf2%build_interface()
         ! Set interface planes at the boundaries
         call vf2%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf2%polygonalize_interface()
         ! Calculate distance from polygons
         call vf2%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf2%subcell_vol()
         ! Calculate curvature
         call vf2%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf2%reset_volume_moments()
      end block create_and_initialize_vof2
      
      
      ! Create a two-phase flow solver with bconds
      create_solver2: block
         use tpns_class, only: dirichlet,clipped_neumann,neumann
         use ils_class,  only: pcg_pfmg
         use mathtools,  only: Pi
         ! Instanciate the solver
         fs2=tpns(cfg=cfg2,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs2%visc_l)
         call param_read('Gas dynamic viscosity'   ,fs2%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs2%rho_l)
         call param_read('Gas density'   ,fs2%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs2%sigma)
         call param_read('Static contact angle',fs2%contact_angle)
         fs2%contact_angle=fs2%contact_angle*Pi/180.0_WP
         ! Define direction gas/liquid stream boundary conditions
         call fs2%add_bcond(name='gas_inj',type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=gas_inj)
         call fs2%add_bcond(name='liq_inj',type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=liq_inj)
         ! Neumann on the side
         call fs2%add_bcond(name='bc_yp'  ,type=clipped_neumann,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call fs2%add_bcond(name='bc_ym'  ,type=clipped_neumann,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         call fs2%add_bcond(name='bc_zp'  ,type=clipped_neumann,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call fs2%add_bcond(name='bc_zm'  ,type=clipped_neumann,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         ! Outflow on the right
         call fs2%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs2%psolv%maxit)
         call param_read('Pressure tolerance',fs2%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs2%implicit%maxit)
         call param_read('Implicit tolerance',fs2%implicit%rcvg)
         ! Setup the solver
         fs2%psolv%maxlevel=16
         call fs2%setup(pressure_ils=pcg_pfmg,implicit_ils=pcg_pfmg)
      end block create_solver2
      
      
      ! Initialize our velocity field
      initialize_velocity2: block
         use mathtools,  only: pi
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         real(WP) :: rad,Uports
         real(WP) :: Q_SLPM,Q_SI,Aport
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         fs2%U=0.0_WP; fs2%V=0.0_WP; fs2%W=0.0_WP
         ! Handle restart
         if (restarted) then
            call df2%pullvar(name='U'  ,var=fs2%U  )
            call df2%pullvar(name='V'  ,var=fs2%V  )
            call df2%pullvar(name='W'  ,var=fs2%W  )
            call df2%pullvar(name='P'  ,var=fs2%P  )
            call df2%pullvar(name='Pjx',var=fs2%Pjx)
            call df2%pullvar(name='Pjy',var=fs2%Pjy)
            call df2%pullvar(name='Pjz',var=fs2%Pjz)
         end if
         ! Apply Dirichlet at liquid injector port
         call param_read('Liquid flow rate (SLPM)',Q_SLPM)
         Q_SI=Q_SLPM*SLPM2SI
         Aport=pi*rli0**2
         Uports=Q_SI/Aport
         call fs2%get_bcond('liq_inj',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs2%U(i,j,k)=Uports
         end do
         ! Apply all other boundary conditions
         call fs2%apply_bcond(time2%t,time2%dt)
         ! Compute MFR through all boundary conditions
         call fs2%get_mfr()
         ! Adjust MFR for global mass balance
         call fs2%correct_mfr()
         ! Compute cell-centered velocity
         call fs2%interp_vel(Ui2,Vi2,Wi2)
         ! Compute divergence
         call fs2%get_div()
      end block initialize_velocity2
      
      
      ! Create a connected-component labeling object
      create_and_initialize_ccl2: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         cc2=ccl(cfg=cfg2,name='CCL')
         cc2%max_interface_planes=2
         cc2%VFlo=VFlo
         cc2%dot_threshold=-0.5_WP
         cc2%thickness_cutoff=filmthickness_over_dx
         ! Perform CCL step
         call cc2%build_lists(VF=vf2%VF,poly=vf2%interface_polygon,U=fs2%U,V=fs2%V,W=fs2%W)
         call cc2%film_classify(Lbary=vf2%Lbary,Gbary=vf2%Gbary)
         call cc2%deallocate_lists()
      end block create_and_initialize_ccl2
      
      
      ! Create an LES model
      create_sgs2: block
         sgs2=sgsmodel(cfg=fs2%cfg,umask=fs2%umask,vmask=fs2%vmask,wmask=fs2%wmask)
         ! Handle restart
         if (restarted) then
            call df2%pullvar(name='LM',var=sgs2%LM)
            call df2%pullvar(name='MM',var=sgs2%MM)
         end if
      end block create_sgs2
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh2: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh2=surfmesh(nvar=1,name='plic')
         smesh2%varname(1)='nplane'
         ! Transfer polygons to smesh
         call vf2%update_surfmesh(smesh2)
         ! Also populate nplane variable
         smesh2%var(1,:)=1.0_WP
         np=0
         do k=vf2%cfg%kmin_,vf2%cfg%kmax_
            do j=vf2%cfg%jmin_,vf2%cfg%jmax_
               do i=vf2%cfg%imin_,vf2%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf2%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf2%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh2%var(1,np)=real(getNumberOfPlanes(vf2%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh2
      
      
      ! Add Ensight output
      create_ensight2: block
         ! Create Ensight output from cfg2
         ens_out2=ensight(cfg2,'atom')
         ! Create event for Ensight output
         ens_evt2=event(time2,'Ensight output')
         call param_read('2 Ensight output period',ens_evt2%tper)
         ! Add variables to output
         call ens_out2%add_vector('velocity',Ui2,Vi2,Wi2)
         call ens_out2%add_scalar('VOF',vf2%VF)
         call ens_out2%add_scalar('curvature',vf2%curv)
         call ens_out2%add_scalar('visc_t',sgs2%visc)
         call ens_out2%add_scalar('filmID',cc2%film_id)
         call ens_out2%add_scalar('filmThickness',cc2%film_thickness)
         call ens_out2%add_surface('vofplic',smesh2)
         ! Output to ensight
         if (ens_evt2%occurs()) call ens_out2%write_data(time2%t)
      end block create_ensight2
      
      
      ! Create a monitor file
      create_monitor2: block
         ! Prepare some info about fields
         call fs2%get_cfl(time2%dt,time2%cfl)
         call fs2%get_max()
         call vf2%get_max()
         ! Create simulation monitor
         mfile2=monitor(fs2%cfg%amRoot,'simulation2')
         call mfile2%add_column(time2%n,'Timestep number')
         call mfile2%add_column(time2%t,'Time')
         call mfile2%add_column(time2%dt,'Timestep size')
         call mfile2%add_column(time2%cfl,'Maximum CFL')
         call mfile2%add_column(fs2%Umax,'Umax')
         call mfile2%add_column(fs2%Vmax,'Vmax')
         call mfile2%add_column(fs2%Wmax,'Wmax')
         call mfile2%add_column(fs2%Pmax,'Pmax')
         call mfile2%add_column(vf2%VFmax,'VOF maximum')
         call mfile2%add_column(vf2%VFmin,'VOF minimum')
         call mfile2%add_column(vf2%VFint,'VOF integral')
         call mfile2%add_column(fs2%divmax,'Maximum divergence')
         call mfile2%add_column(fs2%psolv%it,'Pressure iteration')
         call mfile2%add_column(fs2%psolv%rerr,'Pressure error')
         call mfile2%write()
         ! Create CFL monitor
         cflfile2=monitor(fs2%cfg%amRoot,'cfl2')
         call cflfile2%add_column(time2%n,'Timestep number')
         call cflfile2%add_column(time2%t,'Time')
         call cflfile2%add_column(fs2%CFLc_x,'Convective xCFL')
         call cflfile2%add_column(fs2%CFLc_y,'Convective yCFL')
         call cflfile2%add_column(fs2%CFLc_z,'Convective zCFL')
         call cflfile2%add_column(fs2%CFLv_x,'Viscous xCFL')
         call cflfile2%add_column(fs2%CFLv_y,'Viscous yCFL')
         call cflfile2%add_column(fs2%CFLv_z,'Viscous zCFL')
         call cflfile2%write()
      end block create_monitor2
      
      
      ! *****************************************************************************
      ! **************************** INITIALIZE SOLVER 1 ****************************
      ! *****************************************************************************
      
      ! Allocate work arrays for cfg1 solver
      allocate_work_arrays1: block
         allocate(resU1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         allocate(resV1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         allocate(resW1(cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         allocate(Ui1  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         allocate(Vi1  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         allocate(Wi1  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
         allocate(SR1(6,cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_))
      end block allocate_work_arrays1
      
      
      ! Initialize time tracker
      initialize_timetracker1: block
         time1=timetracker(cfg1%amRoot,name='nozzle_interior')
         call param_read('1 Max timestep size',time1%dtmax)
         call param_read('1 Max cfl number',time1%cflmax)
         time1%dt=time1%dtmax
         time1%itmax=2
         ! Handle restart
         if (restarted) then
            call df1%pullval(name='t' ,val=time1%t )
            call df1%pullval(name='dt',val=time1%dt)
            time1%told=time1%t-time1%dt
         end if
      end block initialize_timetracker1
      
      
      ! Create an incompressible flow solver with bconds
      create_solver1: block
         use incomp_class, only: dirichlet,clipped_neumann
         use ils_class,    only: pcg_pfmg,pcg_amg
         ! Create flow solver
         fs1=incomp(cfg=cfg1,name='Incompressible NS')
         ! Set the flow properties
         fs1%rho=fs2%rho_g
         fs1%visc=fs2%visc_g
         ! Define gas port boundary conditions
         call fs1%add_bcond(name='inj1',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=inj1)
         call fs1%add_bcond(name='inj2',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=inj2)
         call fs1%add_bcond(name='inj3',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=inj3)
         call fs1%add_bcond(name='inj4',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=inj4)
         ! Outflow on the right
         call fs1%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs1%psolv%maxit)
         call param_read('Pressure tolerance',fs1%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs1%implicit%maxit)
         call param_read('Implicit tolerance',fs1%implicit%rcvg)
         ! Setup the solver
         fs1%psolv%maxlevel=16
         call fs1%setup(pressure_ils=pcg_amg,implicit_ils=pcg_pfmg)
      end block create_solver1
      
      
      ! Initialize our velocity field
      initialize_velocity1: block
         use mathtools, only: pi
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         real(WP) :: rad,Uports
         real(WP) :: Q_SLPM,Q_SI,Aport
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         fs1%U=0.0_WP; fs1%V=0.0_WP; fs1%W=0.0_WP
         ! Handle restart
         if (restarted) then
            call df1%pullvar(name='U'   ,var=fs1%U   )
            call df1%pullvar(name='V'   ,var=fs1%V   )
            call df1%pullvar(name='W'   ,var=fs1%W   )
            call df1%pullvar(name='P'   ,var=fs1%P   )
            call df1%pullvar(name='Uold',var=fs1%Uold)
            call df1%pullvar(name='Vold',var=fs1%Vold)
            call df1%pullvar(name='Wold',var=fs1%Wold)
         end if
         ! Read in mean gas
         call param_read('Gas flow rate (SLPM)',Q_SLPM)
         Q_SI=Q_SLPM*SLPM2SI
         Aport=pi/4.0_WP*inj_norm_diam**2
         Uports=Q_SI/(4.0_WP*Aport) ! fix here conversion
         ! Apply Dirichlet at 4 injector ports
         call fs1%get_bcond('inj1',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs1%V(i,j,k)=+Uports
         end do
         call fs1%get_bcond('inj2',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs1%V(i,j,k)=-Uports
         end do
         call fs1%get_bcond('inj3',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs1%W(i,j,k)=+Uports
         end do
         call fs1%get_bcond('inj4',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs1%W(i,j,k)=-Uports
         end do
         ! Apply all other boundary conditions
         call fs1%apply_bcond(time1%t,time1%dt)
         ! Compute MFR through all boundary conditions
         call fs1%get_mfr()
         ! Adjust MFR for global mass balance
         call fs1%correct_mfr()
         ! Compute cell-centered velocity
         call fs1%interp_vel(Ui1,Vi1,Wi1)
         ! Compute divergence
         call fs1%get_div()
      end block initialize_velocity1
      
      
      ! Create an LES model
      create_sgs1: block
         sgs1=sgsmodel(cfg=fs1%cfg,umask=fs1%umask,vmask=fs1%vmask,wmask=fs1%wmask)
         ! Handle restart
         if (restarted) then
            call df1%pullvar(name='LM',var=sgs1%LM)
            call df1%pullvar(name='MM',var=sgs1%MM)
         end if
      end block create_sgs1
      
      
      ! Add Ensight output
      create_ensight1: block
         ! Create Ensight output from cfg
         ens_out1=ensight(cfg1,'nozzle')
         ! Create event for Ensight output
         ens_evt1=event(time1,'Ensight output')
         call param_read('1 Ensight output period',ens_evt1%tper)
         ! Add variables to output
         call ens_out1%add_vector('velocity',Ui1,Vi1,Wi1)
         call ens_out1%add_scalar('visc_t',sgs1%visc)
         ! Output to ensight
         if (ens_evt1%occurs()) call ens_out1%write_data(time1%t)
      end block create_ensight1
      
      
      ! Create a monitor file
      create_monitor1: block
         ! Prepare some info about fields
         call fs1%get_cfl(time1%dt,time1%cfl)
         call fs1%get_max()
         ! Create simulation monitor
         mfile1=monitor(fs1%cfg%amRoot,'simulation1')
         call mfile1%add_column(time1%n,'Timestep number')
         call mfile1%add_column(time1%t,'Time')
         call mfile1%add_column(time1%dt,'Timestep size')
         call mfile1%add_column(time1%cfl,'Maximum CFL')
         call mfile1%add_column(fs1%Umax,'Umax')
         call mfile1%add_column(fs1%Vmax,'Vmax')
         call mfile1%add_column(fs1%Wmax,'Wmax')
         call mfile1%add_column(fs1%Pmax,'Pmax')
         call mfile1%add_column(fs1%divmax,'Maximum divergence')
         call mfile1%add_column(fs1%psolv%it,'Pressure iteration')
         call mfile1%add_column(fs1%psolv%rerr,'Pressure error')
         call mfile1%write()
         ! Create CFL monitor
         cflfile1=monitor(fs1%cfg%amRoot,'cfl1')
         call cflfile1%add_column(time1%n,'Timestep number')
         call cflfile1%add_column(time1%t,'Time')
         call cflfile1%add_column(fs1%CFLc_x,'Convective xCFL')
         call cflfile1%add_column(fs1%CFLc_y,'Convective yCFL')
         call cflfile1%add_column(fs1%CFLc_z,'Convective zCFL')
         call cflfile1%add_column(fs1%CFLv_x,'Viscous xCFL')
         call cflfile1%add_column(fs1%CFLv_y,'Viscous yCFL')
         call cflfile1%add_column(fs1%CFLv_z,'Viscous zCFL')
         call cflfile1%write()
      end block create_monitor1
      
      
      ! *****************************************************************************
      ! **************************** INITIALIZE SOLVER 3 ****************************
      ! *****************************************************************************
      
      ! Allocate work arrays for cfg3 solver
      allocate_work_arrays3: block
         allocate(resU3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(resV3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(resW3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(Ui3  (cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(Vi3  (cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(Wi3  (cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(SR3(6,cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(U2on3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(V2on3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
         allocate(W2on3(cfg3%imino_:cfg3%imaxo_,cfg3%jmino_:cfg3%jmaxo_,cfg3%kmino_:cfg3%kmaxo_))
      end block allocate_work_arrays3
      
      
      ! Initialize time tracker
      initialize_timetracker3: block
         time3=timetracker(cfg3%amRoot,name='spray_region')
         call param_read('3 Max timestep size',time3%dtmax)
         call param_read('3 Max cfl number',time3%cflmax)
         time3%dt=time3%dtmax
         time3%itmax=2
         ! Handle restart
         if (restarted) then
            call df3%pullval(name='t' ,val=time3%t )
            call df3%pullval(name='dt',val=time3%dt)
            time3%told=time3%t-time3%dt
         end if
      end block initialize_timetracker3
      
      
      ! Create an incompressible flow solver with bconds
      create_solver3: block
         use incomp_class, only: dirichlet,clipped_neumann
         use ils_class,    only: pcg_pfmg,pcg_amg
         ! Create flow solver
         fs3=incomp(cfg=cfg3,name='Incompressible NS')
         ! Set the flow properties
         fs3%rho=fs2%rho_g
         fs3%visc=fs2%visc_g
         ! Define direction gas/liquid stream boundary conditions
         call fs3%add_bcond(name='gas_inj',type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=gas_inj_spray)
         ! Neumann on the side
         call fs3%add_bcond(name='bc_yp'  ,type=clipped_neumann,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call fs3%add_bcond(name='bc_ym'  ,type=clipped_neumann,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         call fs3%add_bcond(name='bc_zp'  ,type=clipped_neumann,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call fs3%add_bcond(name='bc_zm'  ,type=clipped_neumann,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         ! Outflow on the right
         call fs3%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs3%psolv%maxit)
         call param_read('Pressure tolerance',fs3%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs3%implicit%maxit)
         call param_read('Implicit tolerance',fs3%implicit%rcvg)
         ! Setup the solver
         fs3%psolv%maxlevel=18
         call fs3%setup(pressure_ils=pcg_amg,implicit_ils=pcg_pfmg)
      end block create_solver3
      
      
      ! Initialize our velocity field
      initialize_velocity3: block
         use mathtools, only: pi
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         real(WP) :: rad,Uports
         real(WP) :: Q_SLPM,Q_SI,Aport
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         fs3%U=0.0_WP; fs3%V=0.0_WP; fs3%W=0.0_WP
         ! Handle restart
         if (restarted) then
            call df3%pullvar(name='U'   ,var=fs3%U   )
            call df3%pullvar(name='V'   ,var=fs3%V   )
            call df3%pullvar(name='W'   ,var=fs3%W   )
            call df3%pullvar(name='P'   ,var=fs3%P   )
         end if
         ! Read in mean gas flow rate
         call param_read('Gas flow rate (SLPM)',Q_SLPM)
         Q_SI=Q_SLPM*SLPM2SI
         Aport=pi*rgi0_spray**2-pi*rlo0_spray**2
         Uports=Q_SI/Aport
         call fs3%get_bcond('gas_inj',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs3%U(i,j,k)=Uports
         end do
         ! Apply all other boundary conditions
         call fs3%apply_bcond(time3%t,time3%dt)
         ! Compute MFR through all boundary conditions
         call fs3%get_mfr()
         ! Adjust MFR for global mass balance
         call fs3%correct_mfr()
         ! Compute cell-centered velocity
         call fs3%interp_vel(Ui3,Vi3,Wi3)
         ! Compute divergence
         call fs3%get_div()
      end block initialize_velocity3
      
      
      ! Create an LES model
      create_sgs3: block
         sgs3=sgsmodel(cfg=fs3%cfg,umask=fs3%umask,vmask=fs3%vmask,wmask=fs3%wmask)
         ! Handle restart
         if (restarted) then
            call df3%pullvar(name='LM',var=sgs3%LM)
            call df3%pullvar(name='MM',var=sgs3%MM)
         end if
      end block create_sgs3
      
      
      ! Initialize our Lagrangian spray solver
      initialize_lpt3: block
         ! Create solver
         lp3=lpt(cfg=cfg3,name='spray')
         ! Get droplet density from the input
         call param_read('Liquid density',lp3%rho)
         ! Handle restarts
         if (restarted) call lp3%read(filename=trim(lpt_file))
      end block initialize_lpt3
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh3: block
         integer :: i
         ! Include an extra variable for droplet diameter
         pmesh3=partmesh(nvar=1,nvec=0,name='lpt')
         pmesh3%varname(1)='diameter'
         ! Transfer particles to pmesh
         call lp3%update_partmesh(pmesh3)
         ! Also populate diameter variable
         do i=1,lp3%np_
            pmesh3%var(1,i)=lp3%p(i)%d
         end do
      end block create_pmesh3
      
      
      ! Add Ensight output
      create_ensight3: block
         ! Create Ensight output from cfg
         ens_out3=ensight(cfg3,'spray')
         ! Create event for Ensight output
         ens_evt3=event(time3,'Ensight output')
         call param_read('3 Ensight output period',ens_evt3%tper)
         ! Add variables to output
         call ens_out3%add_vector('velocity',Ui3,Vi3,Wi3)
         call ens_out3%add_scalar('visc_t',sgs3%visc)
         call ens_out3%add_particle('spray',pmesh3)
         ! Output to ensight
         if (ens_evt3%occurs()) call ens_out3%write_data(time3%t)
      end block create_ensight3
      
      
      ! Create a monitor file
      create_monitor3: block
         ! Prepare some info about fields
         call fs3%get_cfl(time3%dt,time3%cfl)
         call fs3%get_max()
         call lp3%get_max()
         ! Create simulation monitor
         mfile3=monitor(fs3%cfg%amRoot,'simulation3')
         call mfile3%add_column(time3%n,'Timestep number')
         call mfile3%add_column(time3%t,'Time')
         call mfile3%add_column(time3%dt,'Timestep size')
         call mfile3%add_column(time3%cfl,'Maximum CFL')
         call mfile3%add_column(fs3%Umax,'Umax')
         call mfile3%add_column(fs3%Vmax,'Vmax')
         call mfile3%add_column(fs3%Wmax,'Wmax')
         call mfile3%add_column(fs3%Pmax,'Pmax')
         call mfile3%add_column(fs3%divmax,'Maximum divergence')
         call mfile3%add_column(fs3%psolv%it,'Pressure iteration')
         call mfile3%add_column(fs3%psolv%rerr,'Pressure error')
         call mfile3%write()
         ! Create CFL monitor
         cflfile3=monitor(fs3%cfg%amRoot,'cfl3')
         call cflfile3%add_column(time3%n,'Timestep number')
         call cflfile3%add_column(time3%t,'Time')
         call cflfile3%add_column(fs3%CFLc_x,'Convective xCFL')
         call cflfile3%add_column(fs3%CFLc_y,'Convective yCFL')
         call cflfile3%add_column(fs3%CFLc_z,'Convective zCFL')
         call cflfile3%add_column(fs3%CFLv_x,'Viscous xCFL')
         call cflfile3%add_column(fs3%CFLv_y,'Viscous yCFL')
         call cflfile3%add_column(fs3%CFLv_z,'Viscous zCFL')
         call cflfile3%write()
         ! Create a spray monitor
         sprayfile=monitor(amroot=lp3%cfg%amRoot,name='spray')
         call sprayfile%add_column(time3%n,'Timestep number')
         call sprayfile%add_column(time3%t,'Time')
         call sprayfile%add_column(time3%dt,'Timestep size')
         call sprayfile%add_column(lp3%np,'Droplet number')
         call sprayfile%add_column(lp3%Umin, 'Umin')
         call sprayfile%add_column(lp3%Umax, 'Umax')
         call sprayfile%add_column(lp3%Umean,'Umean')
         call sprayfile%add_column(lp3%Vmin, 'Vmin')
         call sprayfile%add_column(lp3%Vmax, 'Vmax')
         call sprayfile%add_column(lp3%Vmean,'Vmean')
         call sprayfile%add_column(lp3%Wmin, 'Wmin')
         call sprayfile%add_column(lp3%Wmax, 'Wmax')
         call sprayfile%add_column(lp3%Wmean,'Wmean')
         call sprayfile%add_column(lp3%dmin, 'dmin')
         call sprayfile%add_column(lp3%dmax, 'dmax')
         call sprayfile%add_column(lp3%dmean,'dmean')
         call sprayfile%write()
      end block create_monitor3
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use parallel,   only: MPI_REAL_WP
      use mpi_f08,    only: MPI_ALLREDUCE,MPI_SUM
      use tpns_class, only: bcond,static_contact
      implicit none
      integer :: n,i,j,k,ierr
      type(bcond), pointer :: mybc
      character(len=str_medium) :: dirname,timestamp
      
      ! Perform time integration - the second solver is the main driver here
      do while (.not.time2%done())
         
         ! Increment time
         call fs2%get_cfl(time2%dt,time2%cfl)
         call time2%adjust_dt()
         call time2%increment()
         
         ! ###############################################
         ! ############ ADVANCE SOLVER 1 HERE ############
         ! ###############################################
         
         ! Advance until we've caught up
         do while (time1%t.lt.time2%t)
            
            ! Increment time
            call fs1%get_cfl(time1%dt,time1%cfl)
            call time1%adjust_dt()
            call time1%increment()
            
            ! Remember old velocity
            fs1%Uold=fs1%U
            fs1%Vold=fs1%V
            fs1%Wold=fs1%W
            
            ! Apply time-varying Dirichlet conditions
            ! This is where time-dpt Dirichlet would be enforced
            
            ! Reset here fluid properties
            fs1%visc=fs2%visc_g
            
            ! Turbulence modeling
            call fs1%get_strainrate(Ui=Ui1,Vi=Vi1,Wi=Wi1,SR=SR1)
            resU1=fs1%rho
            call sgs1%get_visc(dt=time1%dtold,rho=resU1,Ui=Ui1,Vi=Vi1,Wi=Wi1,SR=SR1)
            where (sgs1%visc.lt.-fs1%visc)
               sgs1%visc=-fs1%visc
            end where
            fs1%visc=fs1%visc+sgs1%visc
            
            ! Perform sub-iterations
            do while (time1%it.le.time1%itmax)
               
               ! Build mid-time velocity
               fs1%U=0.5_WP*(fs1%U+fs1%Uold)
               fs1%V=0.5_WP*(fs1%V+fs1%Vold)
               fs1%W=0.5_WP*(fs1%W+fs1%Wold)
               
               ! Explicit calculation of drho*u/dt from NS
               call fs1%get_dmomdt(resU1,resV1,resW1)
               
               ! Assemble explicit residual
               resU1=-2.0_WP*(fs1%rho*fs1%U-fs1%rho*fs1%Uold)+time1%dt*resU1
               resV1=-2.0_WP*(fs1%rho*fs1%V-fs1%rho*fs1%Vold)+time1%dt*resV1
               resW1=-2.0_WP*(fs1%rho*fs1%W-fs1%rho*fs1%Wold)+time1%dt*resW1
               
               ! Form implicit residuals
               call fs1%solve_implicit(time1%dt,resU1,resV1,resW1)
               
               ! Apply these residuals
               fs1%U=2.0_WP*fs1%U-fs1%Uold+resU1
               fs1%V=2.0_WP*fs1%V-fs1%Vold+resV1
               fs1%W=2.0_WP*fs1%W-fs1%Wold+resW1
               
               ! Apply other boundary conditions on the resulting fields
               call fs1%apply_bcond(time1%t,time1%dt)
               
               ! Solve Poisson equation
               call fs1%correct_mfr()
               call fs1%get_div()
               fs1%psolv%rhs=-fs1%cfg%vol*fs1%div*fs1%rho/time1%dt
               fs1%psolv%sol=0.0_WP
               call fs1%psolv%solve()
               
               ! Correct velocity
               call fs1%get_pgrad(fs1%psolv%sol,resU1,resV1,resW1)
               fs1%P=fs1%P+fs1%psolv%sol
               fs1%U=fs1%U-time1%dt*resU1/fs1%rho
               fs1%V=fs1%V-time1%dt*resV1/fs1%rho
               fs1%W=fs1%W-time1%dt*resW1/fs1%rho
               
               ! Increment sub-iteration counter
               time1%it=time1%it+1
               
            end do
            
            ! Recompute interpolated velocity and divergence
            call fs1%interp_vel(Ui1,Vi1,Wi1)
            call fs1%get_div()
            
            ! Output to ensight
            if (ens_evt1%occurs()) call ens_out1%write_data(time1%t)
            
            ! Perform and output monitoring
            call fs1%get_max()
            call mfile1%write()
            call cflfile1%write()
            
         end do
         
         
         ! ###############################################
         ! ######## PERFORM ONE-WAY COUPLING HERE ########
         ! ###############################################
         first_coupling: block
            real(WP) :: wt
            ! At this point, we have ensured that t1>t2, therefore t2 is in [t1old,t1], so interpolate solver 1's velocity in time to t2
            wt=(time2%t-time1%told)/time1%dt
            resU1=wt*fs1%U+(1.0_WP-wt)*fs1%Uold
            resV1=wt*fs1%V+(1.0_WP-wt)*fs1%Vold
            resW1=wt*fs1%W+(1.0_WP-wt)*fs1%Wold
            ! Exchange data using cpl12x/y/z couplers
            call cpl12x%push(resU1); call cpl12x%transfer(); call cpl12x%pull(resU2)
            call cpl12y%push(resV1); call cpl12y%transfer(); call cpl12y%pull(resV2)
            call cpl12z%push(resW1); call cpl12z%transfer(); call cpl12z%pull(resW2)
            ! Apply time-varying Dirichlet conditions
            call fs2%get_bcond('gas_inj',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs2%U(i  ,j,k)=resU2(i,j,k)
               fs2%V(i-1,j,k)=resV2(i,j,k)
               fs2%W(i-1,j,k)=resW2(i,j,k)
            end do
         end block first_coupling
         
         
         ! ###############################################
         ! ############ ADVANCE SOLVER 2 HERE ############
         ! ###############################################
         
         ! Remember old VOF
         vf2%VFold=vf2%VF
         
         ! Remember old velocity
         fs2%Uold=fs2%U
         fs2%Vold=fs2%V
         fs2%Wold=fs2%W
         
         ! Prepare old staggered density (at n)
         call fs2%get_olddensity(vf=vf2)
         
         ! VOF solver step
         call vf2%advance(dt=time2%dt,U=fs2%U,V=fs2%V,W=fs2%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs2%get_viscosity(vf=vf2)
         
         ! Turbulence modeling - only work with gas properties here
         call fs2%get_strainrate(Ui=Ui2,Vi=Vi2,Wi=Wi2,SR=SR2)
         resU2=fs2%rho_g
         call sgs2%get_visc(dt=time2%dtold,rho=resU2,Ui=Ui2,Vi=Vi2,Wi=Wi2,SR=SR2)
         where (sgs2%visc.lt.-fs2%visc_g)
            sgs2%visc=-fs2%visc_g
         end where
         do k=fs2%cfg%kmino_+1,fs2%cfg%kmaxo_
            do j=fs2%cfg%jmino_+1,fs2%cfg%jmaxo_
               do i=fs2%cfg%imino_+1,fs2%cfg%imaxo_
                  fs2%visc(i,j,k)   =fs2%visc(i,j,k)   +sgs2%visc(i,j,k)
                  fs2%visc_xy(i,j,k)=fs2%visc_xy(i,j,k)+sum(fs2%itp_xy(:,:,i,j,k)*sgs2%visc(i-1:i,j-1:j,k))
                  fs2%visc_yz(i,j,k)=fs2%visc_yz(i,j,k)+sum(fs2%itp_yz(:,:,i,j,k)*sgs2%visc(i,j-1:j,k-1:k))
                  fs2%visc_zx(i,j,k)=fs2%visc_zx(i,j,k)+sum(fs2%itp_xz(:,:,i,j,k)*sgs2%visc(i-1:i,j,k-1:k))
               end do
            end do
         end do
         
         ! Perform sub-iterations
         do while (time2%it.le.time2%itmax)
            
            ! Build mid-time velocity
            fs2%U=0.5_WP*(fs2%U+fs2%Uold)
            fs2%V=0.5_WP*(fs2%V+fs2%Vold)
            fs2%W=0.5_WP*(fs2%W+fs2%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs2%prepare_advection_upwind(dt=time2%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs2%get_dmomdt(resU2,resV2,resW2)
            
            ! Assemble explicit residual
            resU2=-2.0_WP*fs2%rho_U*fs2%U+(fs2%rho_Uold+fs2%rho_U)*fs2%Uold+time2%dt*resU2
            resV2=-2.0_WP*fs2%rho_V*fs2%V+(fs2%rho_Vold+fs2%rho_V)*fs2%Vold+time2%dt*resV2
            resW2=-2.0_WP*fs2%rho_W*fs2%W+(fs2%rho_Wold+fs2%rho_W)*fs2%Wold+time2%dt*resW2
            
            ! Form implicit residuals
            call fs2%solve_implicit(time2%dt,resU2,resV2,resW2)
            
            ! Apply these residuals
            fs2%U=2.0_WP*fs2%U-fs2%Uold+resU2
            fs2%V=2.0_WP*fs2%V-fs2%Vold+resV2
            fs2%W=2.0_WP*fs2%W-fs2%Wold+resW2
            
            ! Apply other boundary conditions
            call fs2%apply_bcond(time2%t,time2%dt)
            
            ! Solve Poisson equation
            call fs2%update_laplacian()
            call fs2%correct_mfr()
            call fs2%get_div()
            call fs2%add_surface_tension_jump(dt=time2%dt,div=fs2%div,vf=vf2,contact_model=static_contact)
            fs2%psolv%rhs=-fs2%cfg%vol*fs2%div/time2%dt
            fs2%psolv%sol=0.0_WP
            call fs2%psolv%solve()
            call fs2%shift_p(fs2%psolv%sol)
            
            ! Correct velocity
            call fs2%get_pgrad(fs2%psolv%sol,resU2,resV2,resW2)
            fs2%P=fs2%P+fs2%psolv%sol
            fs2%U=fs2%U-time2%dt*resU2/fs2%rho_U
            fs2%V=fs2%V-time2%dt*resV2/fs2%rho_V
            fs2%W=fs2%W-time2%dt*resW2/fs2%rho_W
            
            ! Increment sub-iteration counter
            time2%it=time2%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs2%interp_vel(Ui2,Vi2,Wi2)
         call fs2%get_div()
         
         ! Output to ensight
         if (ens_evt2%occurs()) then
            ! Update surfmesh object
            update_smesh2: block
               use irl_fortran_interface
               integer :: nplane,np
               ! Transfer polygons to smesh
               call vf2%update_surfmesh(smesh2)
               ! Also populate nplane variable
               smesh2%var(1,:)=1.0_WP
               np=0
               do k=vf2%cfg%kmin_,vf2%cfg%kmax_
                  do j=vf2%cfg%jmin_,vf2%cfg%jmax_
                     do i=vf2%cfg%imin_,vf2%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf2%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf2%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh2%var(1,np)=real(getNumberOfPlanes(vf2%liquid_gas_interface(i,j,k)),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh2
            ! Perform ensight output
            call ens_out2%write_data(time2%t)
         end if
         
         ! Perform and output monitoring
         call fs2%get_max()
         call vf2%get_max()
         call mfile2%write()
         call cflfile2%write()
         
         
         ! ###############################################
         ! ####### TRANSFER DROPLETS FROM 2->3 HERE ######
         ! ###############################################
         ! Perform CCL and transfer
         call transfer_vf_to_drops()
         ! After we're done clip all VOF at the exit area and along the sides - hopefully nothing's left
         do k=fs2%cfg%kmino_,fs2%cfg%kmaxo_
            do j=fs2%cfg%jmino_,fs2%cfg%jmaxo_
               do i=fs2%cfg%imino_,fs2%cfg%imaxo_
                  if (i.ge.vf2%cfg%imax-5) vf2%VF(i,j,k)=0.0_WP
                  if (j.ge.vf2%cfg%jmax-5) vf2%VF(i,j,k)=0.0_WP
                  if (j.le.vf2%cfg%jmin+5) vf2%VF(i,j,k)=0.0_WP
                  if (k.ge.vf2%cfg%kmax-5) vf2%VF(i,j,k)=0.0_WP
                  if (k.le.vf2%cfg%kmin+5) vf2%VF(i,j,k)=0.0_WP
               end do
            end do
         end do
         
         
         ! ###############################################
         ! ############ ADVANCE SOLVER 3 HERE ############
         ! ###############################################
         
         ! Advance until we've caught up
         do while (time3%t.lt.time2%t)
            
            
            ! ######################################
            ! ######## PREPARE NUDGING HERE ########
            ! ######################################
            ! Exchange data using cpl23x/y/z couplers and the most recent velocity
            U2on3=0.0_WP; call cpl23x%push(fs2%U); call cpl23x%transfer(); call cpl23x%pull(U2on3)
            V2on3=0.0_WP; call cpl23y%push(fs2%V); call cpl23y%transfer(); call cpl23y%pull(V2on3)
            W2on3=0.0_WP; call cpl23z%push(fs2%W); call cpl23z%transfer(); call cpl23z%pull(W2on3)
            
            ! Increment time
            call fs3%get_cfl(time3%dt,time3%cfl)
            call time3%adjust_dt()
            call time3%increment()
            
            ! Advance particles by full dt
            resU3=fs3%rho; resV3=fs2%visc_g
            call lp3%advance(dt=time3%dt,U=fs3%U,V=fs3%V,W=fs3%W,rho=resU3,visc=resV3)
            
            ! Remember old velocity
            fs3%Uold=fs3%U
            fs3%Vold=fs3%V
            fs3%Wold=fs3%W
            
            ! Apply time-varying Dirichlet conditions
            ! This is where time-dpt Dirichlet would be enforced
            
            ! Reset here fluid properties
            fs3%visc=fs2%visc_g
            
            ! Turbulence modeling
            call fs3%get_strainrate(Ui=Ui3,Vi=Vi3,Wi=Wi3,SR=SR3)
            resU3=fs3%rho
            call sgs3%get_visc(dt=time3%dtold,rho=resU3,Ui=Ui3,Vi=Vi3,Wi=Wi3,SR=SR3)
            where (sgs3%visc.lt.-fs3%visc)
               sgs3%visc=-fs3%visc
            end where
            fs3%visc=fs3%visc+sgs3%visc
            
            ! Perform sub-iterations
            do while (time3%it.le.time3%itmax)
               
               ! Build mid-time velocity
               fs3%U=0.5_WP*(fs3%U+fs3%Uold)
               fs3%V=0.5_WP*(fs3%V+fs3%Vold)
               fs3%W=0.5_WP*(fs3%W+fs3%Wold)
               
               ! Explicit calculation of drho*u/dt from NS
               call fs3%get_dmomdt(resU3,resV3,resW3)
               
               ! Assemble explicit residual
               resU3=-2.0_WP*(fs3%rho*fs3%U-fs3%rho*fs3%Uold)+time3%dt*resU3
               resV3=-2.0_WP*(fs3%rho*fs3%V-fs3%rho*fs3%Vold)+time3%dt*resV3
               resW3=-2.0_WP*(fs3%rho*fs3%W-fs3%rho*fs3%Wold)+time3%dt*resW3
               
               ! Add nudging term here
               nudge: block
                  real(WP) :: xcoord,ycoord,zcoord
                  do k=fs3%cfg%kmin_,fs3%cfg%kmax_
                     do j=fs3%cfg%jmin_,fs3%cfg%jmax_
                        do i=fs3%cfg%imin_,fs3%cfg%imax_
                           if (fs3%umask(i,j,k).eq.0) then
                              xcoord=max((fs2%cfg%x(fs2%cfg%imax+1)-    fs3%cfg%x (i) )/(fs2%cfg%x(fs2%cfg%imax+1)),0.0_WP)
                              ycoord=max((        0.5_WP*fs2%cfg%yL-abs(fs3%cfg%ym(j)))/(0.5_WP*fs2%cfg%yL        ),0.0_WP)
                              zcoord=max((        0.5_WP*fs2%cfg%zL-abs(fs3%cfg%zm(k)))/(0.5_WP*fs2%cfg%zL        ),0.0_WP)
                              resU3(i,j,k)=resU3(i,j,k)+(U2on3(i,j,k)-fs3%U(i,j,k))*(xcoord*ycoord*zcoord)**2
                           end if
                           if (fs3%vmask(i,j,k).eq.0) then
                              xcoord=max((fs2%cfg%x(fs2%cfg%imax+1)-    fs3%cfg%xm(i) )/(fs2%cfg%x(fs2%cfg%imax+1)),0.0_WP)
                              ycoord=max((        0.5_WP*fs2%cfg%yL-abs(fs3%cfg%y (j)))/(0.5_WP*fs2%cfg%yL        ),0.0_WP)
                              zcoord=max((        0.5_WP*fs2%cfg%zL-abs(fs3%cfg%zm(k)))/(0.5_WP*fs2%cfg%zL        ),0.0_WP)
                              resV3(i,j,k)=resV3(i,j,k)+(V2on3(i,j,k)-fs3%V(i,j,k))*(xcoord*ycoord*zcoord)**2
                           end if
                           if (fs3%wmask(i,j,k).eq.0) then
                              xcoord=max((fs2%cfg%x(fs2%cfg%imax+1)-    fs3%cfg%xm(i) )/(fs2%cfg%x(fs2%cfg%imax+1)),0.0_WP)
                              ycoord=max((        0.5_WP*fs2%cfg%yL-abs(fs3%cfg%ym(j)))/(0.5_WP*fs2%cfg%yL        ),0.0_WP)
                              zcoord=max((        0.5_WP*fs2%cfg%zL-abs(fs3%cfg%z (k)))/(0.5_WP*fs2%cfg%zL        ),0.0_WP)
                              resW3(i,j,k)=resW3(i,j,k)+(W2on3(i,j,k)-fs3%W(i,j,k))*(xcoord*ycoord*zcoord)**2
                           end if
                        end do
                     end do
                  end do
               end block nudge
               
               ! Form implicit residuals
               call fs3%solve_implicit(time3%dt,resU3,resV3,resW3)
               
               ! Apply these residuals
               fs3%U=2.0_WP*fs3%U-fs3%Uold+resU3
               fs3%V=2.0_WP*fs3%V-fs3%Vold+resV3
               fs3%W=2.0_WP*fs3%W-fs3%Wold+resW3
               
               ! Apply other boundary conditions on the resulting fields
               call fs3%apply_bcond(time3%t,time3%dt)
               
               ! Solve Poisson equation
               call fs3%correct_mfr()
               call fs3%get_div()
               fs3%psolv%rhs=-fs3%cfg%vol*fs3%div*fs3%rho/time3%dt
               fs3%psolv%sol=0.0_WP
               call fs3%psolv%solve()
               
               ! Correct velocity
               call fs3%get_pgrad(fs3%psolv%sol,resU3,resV3,resW3)
               fs3%P=fs3%P+fs3%psolv%sol
               fs3%U=fs3%U-time3%dt*resU3/fs3%rho
               fs3%V=fs3%V-time3%dt*resV3/fs3%rho
               fs3%W=fs3%W-time3%dt*resW3/fs3%rho
               
               ! Increment sub-iteration counter
               time3%it=time3%it+1
               
            end do
            
            ! Recompute interpolated velocity and divergence
            call fs3%interp_vel(Ui3,Vi3,Wi3)
            call fs3%get_div()
            
            ! Output to ensight
            if (ens_evt3%occurs()) then
               ! Update partmesh object
               update_pmesh3: block
                  ! Transfer particles to pmesh
                  call lp3%update_partmesh(pmesh3)
                  ! Also populate diameter variable
                  do i=1,lp3%np_
                     pmesh3%var(1,i)=lp3%p(i)%d
                  end do
               end block update_pmesh3
               ! Perform ensight output
               call ens_out3%write_data(time3%t)
            end if
            
            ! Perform and output monitoring
            call fs3%get_max()
            call lp3%get_max()
            call mfile3%write()
            call cflfile3%write()
            call sprayfile%write()
            
         end do
         
         
         ! ###############################################
         ! ############## PERFORM I/O HERE ###############
         ! ###############################################
         
         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            ! Prefix for files
            dirname='restart_'; write(timestamp,'(es12.5)') time2%t
            ! Prepare a new directory
            if (fs1%cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
            ! Populate df1 and write it
            call df1%pushval(name=   't',val=time1%t )
            call df1%pushval(name=  'dt',val=time1%dt)
            call df1%pushvar(name=   'U',var=fs1%U   )
            call df1%pushvar(name=   'V',var=fs1%V   )
            call df1%pushvar(name=   'W',var=fs1%W   )
            call df1%pushvar(name=   'P',var=fs1%P   )
            call df1%pushvar(name='Uold',var=fs1%Uold)
            call df1%pushvar(name='Vold',var=fs1%Vold)
            call df1%pushvar(name='Wold',var=fs1%Wold)
            call df1%pushvar(name=  'LM',var=sgs1%LM )
            call df1%pushvar(name=  'MM',var=sgs1%MM )
            call df1%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.1')
            ! Populate df2 and write it
            call df2%pushval(name=  't',val=time2%t )
            call df2%pushval(name= 'dt',val=time2%dt)
            call df2%pushvar(name=  'U',var=fs2%U   )
            call df2%pushvar(name=  'V',var=fs2%V   )
            call df2%pushvar(name=  'W',var=fs2%W   )
            call df2%pushvar(name=  'P',var=fs2%P   )
            call df2%pushvar(name='Pjx',var=fs2%Pjx )
            call df2%pushvar(name='Pjy',var=fs2%Pjy )
            call df2%pushvar(name='Pjz',var=fs2%Pjz )
            call df2%pushvar(name= 'LM',var=sgs2%LM )
            call df2%pushvar(name= 'MM',var=sgs2%MM )
            call df2%pushvar(name= 'VF',var=vf2%VF  )
            call df2%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.2')
            ! Also output IRL interface
            call vf2%write_interface(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.irl')
            ! Populate df3 and write it
            call df3%pushval(name=   't',val=time3%t )
            call df3%pushval(name=  'dt',val=time3%dt)
            call df3%pushvar(name=   'U',var=fs3%U   )
            call df3%pushvar(name=   'V',var=fs3%V   )
            call df3%pushvar(name=   'W',var=fs3%W   )
            call df3%pushvar(name=   'P',var=fs3%P   )
            call df3%pushvar(name=  'LM',var=sgs3%LM )
            call df3%pushvar(name=  'MM',var=sgs3%MM )
            call df3%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.3')
            ! Also output particles
            call lp3%write(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.lpt')
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Transfer vf to drops
   subroutine transfer_vf_to_drops()
      implicit none
      
      ! Perform a first pass with simplest CCL
      call cc2%build_lists(VF=vf2%VF,U=fs2%U,V=fs2%V,W=fs2%W)
      
      ! Loop through identified detached structs and remove those that are spherical enough
      remove_struct: block
         use mathtools, only: pi
         integer :: m,n,l,i,j,k,np
         real(WP) :: lmin,lmax,eccentricity,diam
         
         ! Loops over film segments contained locally
         do m=1,cc2%n_meta_struct
            
            ! Test if sphericity is compatible with transfer
            lmin=cc2%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=cc2%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=cc2%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
            if (eccentricity.gt.max_eccentricity) cycle
            
            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*cc2%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if (diam.eq.0.0_WP.or.diam.gt.d_threshold) cycle
            
            ! Create drop from available liquid volume - only one root does that
            if (cc2%cfg%amRoot) then
               ! Make room for new drop
               np=lp3%np_+1; call lp3%resize(np)
               ! Add the drop
               lp3%p(np)%id  =int(0,8)                                                                                      !< Give id (maybe based on break-up model?)
               lp3%p(np)%dt  =0.0_WP                                                                                        !< Let the drop find it own integration time
               lp3%p(np)%col =0.0_WP                                                                                        !< Give zero collision force               
               lp3%p(np)%d   =diam                                                                                          !< Assign diameter to account for full volume
               lp3%p(np)%pos =[cc2%meta_structures_list(m)%x,cc2%meta_structures_list(m)%y,cc2%meta_structures_list(m)%z]   !< Place the drop at the liquid barycenter
               lp3%p(np)%vel =[cc2%meta_structures_list(m)%u,cc2%meta_structures_list(m)%v,cc2%meta_structures_list(m)%w]   !< Assign mean structure velocity as drop velocity
               lp3%p(np)%ind =lp3%cfg%get_ijk_global(lp3%p(np)%pos,[lp3%cfg%imin,lp3%cfg%jmin,lp3%cfg%kmin])                !< Place the drop in the proper cell for the lp%cfg
               lp3%p(np)%flag=0                                                                                             !< Activate it
               ! Increment particle counter
               lp3%np_=np
            end if
            
            ! Find local structs with matching id
            do n=cc2%sync_offset+1,cc2%sync_offset+cc2%n_struct
               if (cc2%struct_list(cc2%struct_map_(n))%parent.ne.cc2%meta_structures_list(m)%id) cycle
               ! Remove liquid in meta-structure cells
               do l=1,cc2%struct_list(cc2%struct_map_(n))%nnode ! Loops over cells within local
                  i=cc2%struct_list(cc2%struct_map_(n))%node(1,l)
                  j=cc2%struct_list(cc2%struct_map_(n))%node(2,l)
                  k=cc2%struct_list(cc2%struct_map_(n))%node(3,l)
                  ! Remove liquid in that cell
                  vf2%VF(i,j,k)=0.0_WP
               end do
            end do
            
         end do
         
      end block remove_struct
      
      ! Sync VF and clean up IRL and band
      call vf2%cfg%sync(vf2%VF)
      call vf2%clean_irl_and_band()
      
      ! Clean up CCL
      call cc2%deallocate_lists()
      
      ! Perform more detailed CCL in a second pass
      cc2%max_interface_planes=2
      call cc2%build_lists(VF=vf2%VF,poly=vf2%interface_polygon,U=fs2%U,V=fs2%V,W=fs2%W)
      call cc2%get_min_thickness()
      call cc2%sort_by_thickness()
      
      ! Loop through identified films and remove those that are thin enough
      remove_film: block
         use mathtools, only: pi
         integer :: m,n,i,j,k,np,ip,np_old
         real(WP) :: Vt,Vl,Hl,Vd
         
         ! Loops over film segments contained locally
         do m=cc2%film_sync_offset+1,cc2%film_sync_offset+cc2%n_film
            
            ! Skip non-liquid films
            if (cc2%film_list(cc2%film_map_(m))%phase.ne.1) cycle
            
            ! Skip films that are still thick enough
            if (cc2%film_list(cc2%film_map_(m))%min_thickness.gt.min_filmthickness) cycle
            
            ! We are still here: transfer the film to drops
            Vt=0.0_WP       ! Transferred volume
            Vl=0.0_WP       ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
            np_old=lp3%np_  ! Remember old number of particles
            do n=1,cc2%film_list(cc2%film_map_(m))%nnode ! Loops over cells within local film segment
               i=cc2%film_list(cc2%film_map_(m))%node(1,n)
               j=cc2%film_list(cc2%film_map_(m))%node(2,n)
               k=cc2%film_list(cc2%film_map_(m))%node(3,n)
               ! Increment liquid volume to remove
               Vl=Vl+vf2%VF(i,j,k)*vf2%cfg%vol(i,j,k)
               ! Estimate drop size based on local film thickness in current cell
               Hl=max(cc2%film_thickness(i,j,k),min_filmthickness)
               Vd=pi/6.0_WP*(diam_over_filmthickness*Hl)**3
               ! Create drops from available liquid volume
               do while (Vl-Vd.gt.0.0_WP)
                  ! Make room for new drop
                  np=lp3%np_+1; call lp3%resize(np)
                  ! Add the drop
                  lp3%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                  lp3%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  lp3%p(np)%col =0.0_WP                                     !< Give zero collision force
                  lp3%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                  lp3%p(np)%pos =vf2%Lbary(:,i,j,k)                         !< Place the drop at the liquid barycenter
                  lp3%p(np)%vel =fs2%cfg%get_velocity(pos=lp3%p(np)%pos,i0=i,j0=j,k0=k,U=fs2%U,V=fs2%V,W=fs2%W) !< Interpolate local cell velocity as drop velocity
                  lp3%p(np)%ind =lp3%cfg%get_ijk_global(lp3%p(np)%pos,[lp3%cfg%imin,lp3%cfg%jmin,lp3%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                  lp3%p(np)%flag=0                                          !< Activate it
                  ! Increment particle counter
                  lp3%np_=np
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
               end do
               ! Remove liquid in that cell
               vf2%VF(i,j,k)=0.0_WP
            end do
            
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
               ! Add one last drop for remaining liquid volume
               np=lp3%np_+1; call lp3%resize(np)
               ! Add the drop
               lp3%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
               lp3%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               lp3%p(np)%Acol=0.0_WP                                     !< Give zero collision force (axial)
               lp3%p(np)%Tcol=0.0_WP                                     !< Give zero collision force (tangential)
               lp3%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
               lp3%p(np)%pos =vf2%Lbary(:,i,j,k)                         !< Place the drop at the liquid barycenter
               lp3%p(np)%vel =fs2%cfg%get_velocity(pos=lp3%p(np)%pos,i0=i,j0=j,k0=k,U=fs2%U,V=fs2%V,W=fs2%W) !< Interpolate local cell velocity as drop velocity
               lp3%p(np)%ind =lp3%cfg%get_ijk_global(lp3%p(np)%pos,[lp3%cfg%imin,lp3%cfg%jmin,lp3%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
               lp3%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               lp3%np_=np
            else ! Some particles were created, make them all larger
               do ip=np_old+1,lp3%np_
                  lp3%p(ip)%d=lp3%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
               end do
            end if
         end do
         
      end block remove_film
      
      ! Sync VF and clean up IRL and band
      call vf2%cfg%sync(vf2%VF)
      call vf2%clean_irl_and_band()
      
      ! Clean up CCL
      call cc2%deallocate_lists()
      
      ! Resync the spray
      call lp3%sync()
      
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
      deallocate(resU1,resV1,resW1,Ui1,Vi1,Wi1,SR1)
      deallocate(resU2,resV2,resW2,Ui2,Vi2,Wi2,SR2)
      deallocate(resU3,resV3,resW3,Ui3,Vi3,Wi3,SR3,U2on3,V2on3,W2on3)
      
   end subroutine simulation_final
   
   
end module simulation
