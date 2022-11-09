!> Definition for a block 2 simulation: zoomed in cough machine
module block2_class
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   use config_class,      only: config
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use ccl_class,         only: ccl
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   use object_timer,      only: objtimer
   implicit none
   private

   public :: block2,filmthickness_over_dx,min_filmthickness,diam_over_filmthickness,max_eccentricity,d_threshold

   !> Block 2 object
   type :: block2
      class(config), pointer :: cfg                                        !< Pointer to config
      type(tpns) :: fs                                                     !< Two-phase incompressible flow solver
      type(vfs) :: vf                                                      !< VF solver
      type(ccl) :: cc2                                                     !< Connected component labeling
      type(timetracker) :: time                                            !< Time tracker
      type(objtimer) :: timer                                              !< Method timer
      type(sgsmodel) ::  sgs                                               !< SGS model
      type(surfmesh) :: smesh                                              !< Surfmesh 
      type(ensight) :: ens_out                                             !< Ensight output
      type(event) :: ens_evt                                               !< Ensight output event
      type(monitor) :: mfile,cflfile,volfile,timerfile,timersummaryfile    !< Monitor files
      type(datafile) :: df                                                 !< Datafile for restart
      !> Private work arrays
      real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:), allocatable :: SR
      !> Nudging region
      real(WP) :: nudge_trans
      real(WP) :: nudge_xmin,nudge_xmax
      real(WP) :: nudge_ymin,nudge_ymax
      real(WP) :: nudge_zmin,nudge_zmax
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block2

   !> Transfer model parameters
   real(WP) :: filmthickness_over_dx  =5.0e-1_WP   !Model parameter

   real(WP) :: min_filmthickness      =1.0e-7_WP   !Model parameter
   real(WP) :: diam_over_filmthickness=1.0e+1_WP   !Model parameter
   
   real(WP) :: max_eccentricity       =2.0e-1_WP   !Saliva specific?
   real(WP) :: d_threshold            =1.0e-3_WP   ! Uperbound of respiratory droplets generated in oral cavity 

   !> Inflow parameters
   real(WP) :: Uin,delta,Urand,Uco,CPFR

   !> Parameter ratios
   real(WP) :: rhog_rhol,mug_mul
   
   !> Liquid viscosity
   real(WP) :: visc_l

   !> Logical constant for evaluating restart
   logical :: restart_test
   
   
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
   
   
   !> Function that localizes the left domain boundary, outside the mouth
   function left_boundary_coflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.(pg%ym(j).le.-t_wall.or.pg%ym(j).ge.H_mouth+t_wall.or.abs(pg%zm(k)).ge.0.5_WP*W_mouth+t_wall)) isIn=.true.
   end function left_boundary_coflow
   
   
   !> Function that localizes the rightmost domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary

   !> Function that localizes the left domain boundary at imin
   function left_boundary_u_inflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j).gt.0.0_WP) isIn=.true.
   end function left_boundary_u_inflow

   !> Function that localizes the left domain boundary at imin-1
   function left_boundary_vw_inflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1.and.pg%ym(j).gt.0.0_WP) isIn=.true.
   end function left_boundary_vw_inflow

   !> Function that localizes the rightmost domain boundary
   function right_boundary_outflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1.and.pg%ym(j).gt.0.0_WP) isIn=.true.
   end function right_boundary_outflow
   
   !> Function that calcuates velocity at current time 
   function inflowVelocity(time,CPFR,H,W) result(UCPFR)
      real(WP), intent(in)   :: time,CPFR,H,W
      real(WP)               :: UCPFR,CEV,PVT,a1,b1,c1,a2,b2,c2,tau,M
      ! Model parameters
      CEV=0.20_WP*CPFR-4e-5_WP    ! cough expiratory volume
      PVT=2.85_WP*CPFR+0.07_WP    ! Peak velocity time
      a1=1.68_WP
      b1=3.34_WP
      c1=0.43_WP
      a2=(CEV/(PVT*CPFR))-a1
      b2=((-2.16_WP*CEV)/(PVT*CPFR))+10.46_WP
      c2=(1.8_WP/(b2-1.0_WP))
      ! Dimensionless time
      tau=time/PVT
      ! Dimensionless flow rate
      if (tau.eq.0.0_WP) then
         M=0.0_WP
      else if (tau.lt.1.2_WP.and.tau.gt.0.0_WP) then 
         M=(a1*tau**(b1-1.0_WP)*exp(-tau/c1))/(gamma(b1)*c1**b1)    
      else if (tau.ge.1.2_WP) then 
         M=(a1*tau**(b1-1.0_WP)*exp(-tau/c1))/(gamma(b1)*c1**b1)+(a2*(tau-1.2_WP)**(b2-1.0_WP)*exp(-(tau-1.2_WP)/c2))/(gamma(b2)*c2**b2)   
      end if
      ! Inflow velocity 
      UCPFR=(M*CPFR)/(H*W)
   end function inflowVelocity
   
   !> Function that localizes the top domain boundary
   ! function top_boundaryV(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (j.eq.pg%jmax+1) isIn=.true.
   ! end function top_boundaryV
   
   
   !> Function that localizes the bottom domain boundary
   ! function bottom_boundaryV(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (j.eq.pg%jmin) isIn=.true.
   ! end function bottom_boundaryV
   
   
   !> Function that localizes the top domain boundary
   ! function top_boundaryU(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (j.eq.pg%jmax+1) isIn=.true.
   ! end function top_boundaryU
   
   
   !> Function that localizes the bottom domain boundary
   ! function bottom_boundaryU(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (j.eq.pg%jmin-1) isIn=.true.
   ! end function bottom_boundaryU
   
   
   !> Initialization of block 2
   subroutine init(b,restart_test)
      use param, only: param_read
      implicit none
      class(block2), intent(inout) :: b
      logical,       intent(in) :: restart_test
   
      ! Allocate work arrays for cfg
      allocate_work_arrays: block
         allocate(b%resU(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resV(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resW(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Ui  (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Vi  (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Wi  (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%SR(6,b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
      end block allocate_work_arrays

      ! Initialize time tracker
      initialize_timetracker: block
         !b%time=timetracker(b%cfg%amRoot,name='cough_machine_in')
         call param_read('2 Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         call param_read('Max time',b%time%tmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
         ! Handle restart
         if (restart_test) then
            call b%df%pullval(name='t' ,val=b%time%t )
            call b%df%pullval(name='dt',val=b%time%dt)
            b%time%told=b%time%t-b%time%dt
         end if
      end block initialize_timetracker

      ! Initalize object time tracker
      initialize_objtimer: block
         b%timer=objtimer(b%cfg%amRoot,name='cough_in_timer')
      end block initialize_objtimer

      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: lvira,r2p
         integer :: i,j,k
         ! Create a VOF solver with LVIRA
         !b%vf=vfs(cfg=b%cfg,reconstruction_method=lvira,name='VOF')
         ! Create a VOF solver with R2P
         b%vf=vfs(cfg=b%cfg,reconstruction_method=r2p,name='VOF')
         ! Initialize to flat interface in liquid tray
         do k=b%vf%cfg%kmino_,b%vf%cfg%kmaxo_
            do j=b%vf%cfg%jmino_,b%vf%cfg%jmaxo_
               do i=b%vf%cfg%imino_,b%vf%cfg%imaxo_
                  ! if (b%vf%cfg%xm(i).lt.-L_lip.and.b%vf%cfg%xm(i).gt.-L_lip-L_film.and.abs(b%vf%cfg%zm(k)).lt.0.5_WP*W_film.and.b%vf%cfg%ym(j).lt.0.0_WP.and.b%vf%cfg%ym(j).gt.-H_film) then !original 
                  if (b%vf%cfg%xm(i).gt.L_lip.and.b%vf%cfg%xm(i).lt.L_lip+L_film.and.b%vf%cfg%ym(j).lt.0.0_WP.and.b%vf%cfg%ym(j).gt.-H_film) then
                     b%vf%VF(i,j,k)=1.0_WP
                  else
                     b%vf%VF(i,j,k)=0.0_WP
                  end if
                  b%vf%Lbary(:,i,j,k)=[b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]
                  b%vf%Gbary(:,i,j,k)=[b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Handle restart - using VF data
         if (restart_test) call b%df%pullvar(name='VF',var=b%vf%VF)
         ! Update the band
         call b%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call b%vf%build_interface()
         ! Set interface planes at the boundaries
         call b%vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call b%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call b%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call b%vf%subcell_vol()
         ! Calculate curvature
         call b%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call b%vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a two-phase flow solver with bconds
      create_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,neumann
         use ils_class,  only: pcg_pfmg,gmres_amg,gmres_pilut !Didn't work for pressure_ils: gmres_pilut, pcg_pfmg
         use mathtools,  only: Pi
         integer  :: i,j,k
         ! Create a two-phase flow solver
         b%fs=tpns(cfg=b%cfg,name='Two-phase NS')
         ! Assign gas viscosity and density
         call param_read('Gas dynamic viscosity',b%fs%visc_g)
         call param_read('Gas density'          ,b%fs%rho_g)
         ! Assign liquid viscosity and density based on ratios
         call param_read('Viscosity ratio',mug_mul)
         b%fs%visc_l=b%fs%visc_g/mug_mul
         call param_read('Density ratio',  rhog_rhol)
         b%fs%rho_l=b%fs%rho_g/rhog_rhol
         ! Assign viscosity to each phase
         ! call param_read('Liquid dynamic viscosity',visc_l)
         ! b%fs%visc_l=visc_l
         ! call param_read('Gas dynamic viscosity'   ,b%fs%visc_g)
         ! Assign constant density to each phase
         ! call param_read('Liquid density',b%fs%rho_l)
         ! call param_read('Gas density'   ,b%fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',b%fs%sigma)
         call param_read('Static contact angle',b%fs%contact_angle)
         b%fs%contact_angle=b%fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',b%fs%gravity)
         ! Inflow on the left
         call b%fs%add_bcond(name='u_inflow' ,type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_u_inflow )
         call b%fs%add_bcond(name='vw_inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_vw_inflow)
         ! call b%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_mouth)
         ! call b%fs%add_bcond(name='coflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_coflow)
         ! Outflow on the right
         call b%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary_outflow)
         ! Configure pressure solver
         call param_read('Pressure iteration',b%fs%psolv%maxit)
         call param_read('Pressure tolerance',b%fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',b%fs%implicit%maxit)
         call param_read('Implicit tolerance',b%fs%implicit%rcvg)
         ! Setup the solver
         b%fs%psolv%maxlevel=24
         call b%fs%setup(pressure_ils=gmres_amg,implicit_ils=gmres_amg)
      end block create_solver


      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         use random,     only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! Zero initial field
         b%fs%U=0.0_WP; b%fs%V=0.0_WP; b%fs%W=0.0_WP
         ! Handle restart
         if (restart_test) then
            call b%df%pullvar(name='U'  ,var=b%fs%U  )
            call b%df%pullvar(name='V'  ,var=b%fs%V  )
            call b%df%pullvar(name='W'  ,var=b%fs%W  )
            call b%df%pullvar(name='P'  ,var=b%fs%P  )
         end if
         ! Gas velocity parameters
         ! call param_read('Gas velocity',Uin)
         ! call param_read('Peak flow rate',CPFR)
         ! Uin=inflowVelocity(b%time%t,CPFR,H_mouth,W_mouth)
         ! call param_read('Gas thickness',delta)
         ! call param_read('Gas perturbation',Urand)
         ! call b%fs%get_bcond('inflow',mybc)
         ! ! Apply Dirichlet at inlet
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    !b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
         !    b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)
         ! end do
         ! Apply coflow around inlet geometry
         !call param_read('Gas coflow',Uco)
         ! Uco=0.10_WP*Uin
         ! call b%fs%get_bcond('coflow',mybc)
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    b%fs%U(i,j,k)=Uco
         ! end do
         ! Apply all other boundary conditions
         call b%fs%apply_bcond(b%time%t,b%time%dt)
         ! Compute MFR through all boundary conditions
         call b%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call b%fs%correct_mfr()
         ! Compute cell-centered velocity
         call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
         ! Compute divergence
         call b%fs%get_div()
      end block initialize_velocity

      ! Create a connected-component labeling object
      create_and_initialize_ccl2: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         b%cc2=ccl(cfg=b%cfg,name='CCL')
         b%cc2%max_interface_planes=2
         b%cc2%VFlo=VFlo
         b%cc2%dot_threshold=-0.5_WP
         b%cc2%thickness_cutoff=filmthickness_over_dx
         ! Perform CCL step
         call b%cc2%build_lists(VF=b%vf%VF,poly=b%vf%interface_polygon,U=b%fs%U,V=b%fs%V,W=b%fs%W)
         call b%cc2%film_classify(Lbary=b%vf%Lbary,Gbary=b%vf%Gbary)
         call b%cc2%deallocate_lists()
      end block create_and_initialize_ccl2

      ! Create an LES model
      create_sgs: block
         b%sgs=sgsmodel(cfg=b%fs%cfg,umask=b%fs%umask,vmask=b%fs%vmask,wmask=b%fs%wmask)
         ! Handle restart
         if (restart_test) then
            call b%df%pullvar(name='LM',var=b%sgs%LM)
            call b%df%pullvar(name='MM',var=b%sgs%MM)
         end if
      end block create_sgs

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         b%smesh=surfmesh(nvar=1,name='plic')
         b%smesh%varname(1)='nplane'
         ! Transfer polygons to smesh
         call b%vf%update_surfmesh(b%smesh)
         ! Also populate nplane variable
         b%smesh%var(1,:)=1.0_WP
         np=0
         do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
            do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
               do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(b%cfg,'cough_in')
         ! Create event for Ensight output
         b%ens_evt=event(b%time,'Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call b%ens_out%add_scalar('VOF',b%vf%VF)
         call b%ens_out%add_scalar('curvature',b%vf%curv)
         call b%ens_out%add_scalar('visc_t',b%sgs%visc)
         call b%ens_out%add_surface('vofplic',b%smesh)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         call b%vf%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation2')
         call b%mfile%add_column(b%time%n,'Timestep number')
         call b%mfile%add_column(b%time%t,'Time')
         call b%mfile%add_column(b%time%dt,'Timestep size')
         call b%mfile%add_column(b%time%cfl,'Maximum CFL')
         call b%mfile%add_column(b%fs%Umax,'Umax')
         call b%mfile%add_column(b%fs%Vmax,'Vmax')
         call b%mfile%add_column(b%fs%Wmax,'Wmax')
         call b%mfile%add_column(b%fs%Pmax,'Pmax')
         call b%mfile%add_column(b%vf%VFmax,'VOF maximum')
         call b%mfile%add_column(b%vf%VFmin,'VOF minimum')
         call b%mfile%add_column(b%vf%VFint,'VOF integral')
         call b%mfile%add_column(b%fs%divmax,'Maximum divergence')
         call b%mfile%add_column(b%fs%psolv%it,'Pressure iteration')
         call b%mfile%add_column(b%fs%psolv%rerr,'Pressure error')
         call b%mfile%add_column(Uin,'Inflow Velocity')
         call b%mfile%write()
         ! Create CFL monitor
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl2')
         call b%cflfile%add_column(b%time%n,'Timestep number')
         call b%cflfile%add_column(b%time%t,'Time')
         call b%cflfile%add_column(b%fs%CFLc_x,'Convective xCFL')
         call b%cflfile%add_column(b%fs%CFLc_y,'Convective yCFL')
         call b%cflfile%add_column(b%fs%CFLc_z,'Convective zCFL')
         call b%cflfile%add_column(b%fs%CFLv_x,'Viscous xCFL')
         call b%cflfile%add_column(b%fs%CFLv_y,'Viscous yCFL')
         call b%cflfile%add_column(b%fs%CFLv_z,'Viscous zCFL')
         call b%cflfile%write()
         ! Create 0 volume struct monitor
         b%volfile=monitor(b%fs%cfg%amRoot,'Zero_Vol_Struct')
         call b%volfile%add_column(b%time%n,'Timestep number')
         call b%volfile%add_column(b%time%t,'Time')
         call b%volfile%add_column(b%cc2%zero_struct_id,'Structure ID')
         call b%volfile%add_column(b%cc2%zero_struct_vol,'Structure volume')
         call b%volfile%write()
         ! Create object time tracker monitor
         b%timerfile=monitor(b%fs%cfg%amRoot,'cough_in_timers')
         call b%timerfile%add_column(b%time%n,'Timestep number')
         call b%timerfile%add_column(b%time%t,'Simulation Time')
         call b%timerfile%add_column(b%timer%vf_wt,'VF_advance Wall Time')
         call b%timerfile%add_column(b%timer%sgs_wt,'sgs_visc Wall Time')
         call b%timerfile%add_column(b%timer%implicit_wt,'imp_solv Wall Time')
         call b%timerfile%add_column(b%timer%pressure_wt,'pres_solv Wall Time')
         call b%timerfile%add_column(b%timer%step_wt,'time_step Wall Time')
         call b%timerfile%write()
         ! Create object time and cost summary monitor
         b%timersummaryfile=monitor(b%fs%cfg%amRoot,'cough_in_timer_summary')
         call b%timersummaryfile%add_column(b%time%n,'Timestep number')
         call b%timersummaryfile%add_column(b%time%t,'Simulation Time')
         call b%timersummaryfile%add_column(b%timer%vf_wt_total,'VF_advance Total Hours')
         call b%timersummaryfile%add_column(b%timer%vf_core_hours,'VF_advance Core Hours')
         call b%timersummaryfile%add_column(b%timer%sgs_wt_total,'sgs_visc Total Hours')
         call b%timersummaryfile%add_column(b%timer%sgs_core_hours,'sgs_visc Core Hours')
         call b%timersummaryfile%add_column(b%timer%implicit_wt_total,'imp_solv WT Hours')
         call b%timersummaryfile%add_column(b%timer%implicit_core_hours,'imp_solv Core Hours')
         call b%timersummaryfile%add_column(b%timer%pressure_wt_total,'pres_solv WT Hours')
         call b%timersummaryfile%add_column(b%timer%pressure_core_hours,'pres_solv Core Hours')
         call b%timersummaryfile%add_column(b%timer%step_wt_total,'time_step WT Hours')
         call b%timersummaryfile%add_column(b%timer%step_core_hours,'time_step Core Hours')
         call b%timersummaryfile%write()
      end block create_monitor


   end subroutine init


   !> Take a time step with block 2
   subroutine step(b,Unudge,Vnudge,Wnudge)
      use tpns_class, only: static_contact
      use mpi,        only: mpi_wtime
      implicit none
      class(block2), intent(inout) :: b
      real(WP), dimension(b%cfg%imino_:,b%cfg%jmino_:,b%cfg%kmino_:), intent(inout) :: Unudge     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(b%cfg%imino_:,b%cfg%jmino_:,b%cfg%kmino_:), intent(inout) :: Vnudge     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(b%cfg%imino_:,b%cfg%jmino_:,b%cfg%kmino_:), intent(inout) :: Wnudge     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: starttime,endtime

      ! Start time step timer
      starttime=mpi_wtime()

      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! Apply time-varying Dirichlet conditions
      reapply_dirichlet: block
         use tpns_class, only: bcond
         use random,     only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! U velocity inflow from duct
         call b%fs%get_bcond('u_inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            b%fs%U(i,j,k)=Unudge(i,j,k)
         end do
         ! V and W velocity inflow from duct
         call b%fs%get_bcond('vw_inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            b%fs%V(i,j,k)=Vnudge(i,j,k)
            b%fs%W(i,j,k)=Wnudge(i,j,k)
         end do
         ! Reapply Dirichlet at inlet
         ! Uin=inflowVelocity(b%time%t,CPFR,H_mouth,W_mouth)
         ! U velocity inflow from duct
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    b%fs%U(i,j,k)=Uin!*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
         !    b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)
         ! end do
         ! Reapply coflow around inlet geometry
         ! Uco=0.10_WP*Uin
         ! call b%fs%get_bcond('coflow',mybc)
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    b%fs%U(i,j,k)=Uco
         ! end do
      end block reapply_dirichlet

      ! Calculate SR
      call b%fs%get_strainrate(Ui=b%Ui,Vi=b%Vi,Wi=b%Wi,SR=b%SR)

      ! ! Model non-Newtonian fluid
      ! nonewt: block
      !    integer :: i,j,k
      !    real(WP) :: SRmag
      !    real(WP), parameter :: C=1.137e-3_WP
      !    real(WP), parameter :: n=0.3_WP 
      !    ! Update viscosity
      !    do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
      !       do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
      !          do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
      !             SRmag=sqrt(b%SR(1,i,j,k)**2+b%SR(2,i,j,k)**2+b%SR(3,i,j,k)**2+2.0_WP*(b%SR(4,i,j,k)**2+b%SR(5,i,j,k)**2+b%SR(6,i,j,k)**2))
      !             SRmag=max(SRmag,1000.0_WP**(1.0_WP/(n-1.0_WP)))
      !             b%fs%visc_l(i,j,k)=C*SRmag**(n-1.0_WP)
      !          end do
      !       end do
      !    end do
      !    call b%fs%cfg%sync(b%fs%visc_l)
      ! end block nonewt

      ! Remember old VOF
      b%vf%VFold=b%vf%VF

      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W
      
      ! Prepare old staggered density (at n)
      call b%fs%get_olddensity(vf=b%vf)
      
      ! VOF solver step
      call b%vf%advance(dt=b%time%dt,U=b%fs%U,V=b%fs%V,W=b%fs%W)

      ! Prepare new staggered constant viscosity (at n+1)
      call b%fs%get_viscosity(vf=b%vf)

      ! Turbulence modeling - only work with gas properties here
      sgs_model: block
         integer :: i,j,k
         b%resU=b%fs%rho_g
         call b%sgs%get_visc(dt=b%time%dtold,rho=b%resU,Ui=b%Ui,Vi=b%Vi,Wi=b%Wi,SR=b%SR)
         where (b%sgs%visc.lt.-b%fs%visc_g)
            b%sgs%visc=-b%fs%visc_g
         end where
         do k=b%fs%cfg%kmino_+1,b%fs%cfg%kmaxo_
            do j=b%fs%cfg%jmino_+1,b%fs%cfg%jmaxo_
               do i=b%fs%cfg%imino_+1,b%fs%cfg%imaxo_
                  b%fs%visc(i,j,k)   =b%fs%visc(i,j,k)   +b%sgs%visc(i,j,k)
                  b%fs%visc_xy(i,j,k)=b%fs%visc_xy(i,j,k)+sum(b%fs%itp_xy(:,:,i,j,k)*b%sgs%visc(i-1:i,j-1:j,k))
                  b%fs%visc_yz(i,j,k)=b%fs%visc_yz(i,j,k)+sum(b%fs%itp_yz(:,:,i,j,k)*b%sgs%visc(i,j-1:j,k-1:k))
                  b%fs%visc_zx(i,j,k)=b%fs%visc_zx(i,j,k)+sum(b%fs%itp_xz(:,:,i,j,k)*b%sgs%visc(i-1:i,j,k-1:k))
               end do
            end do
         end do
      end block sgs_model

      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)

         ! Build mid-time velocity
         b%fs%U=0.5_WP*(b%fs%U+b%fs%Uold)
         b%fs%V=0.5_WP*(b%fs%V+b%fs%Vold)
         b%fs%W=0.5_WP*(b%fs%W+b%fs%Wold)

         ! Preliminary mass and momentum transport step at the interface
         call b%fs%prepare_advection_upwind(dt=b%time%dt)

         ! Explicit calculation of drho*u/dt from NS
         call b%fs%get_dmomdt(b%resU,b%resV,b%resW)

         ! Add momentum source terms
         call b%fs%addsrc_gravity(b%resU,b%resV,b%resW)

         ! Assemble explicit residual
         b%resU=-2.0_WP*b%fs%rho_U*b%fs%U+(b%fs%rho_Uold+b%fs%rho_U)*b%fs%Uold+b%time%dt*b%resU
         b%resV=-2.0_WP*b%fs%rho_V*b%fs%V+(b%fs%rho_Vold+b%fs%rho_V)*b%fs%Vold+b%time%dt*b%resV
         b%resW=-2.0_WP*b%fs%rho_W*b%fs%W+(b%fs%rho_Wold+b%fs%rho_W)*b%fs%Wold+b%time%dt*b%resW

         ! Form implicit residuals
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)
         
         ! Apply these residuals
         b%fs%U=2.0_WP*b%fs%U-b%fs%Uold+b%resU
         b%fs%V=2.0_WP*b%fs%V-b%fs%Vold+b%resV
         b%fs%W=2.0_WP*b%fs%W-b%fs%Wold+b%resW

         ! Apply other boundary conditions
         call b%fs%apply_bcond(b%time%t,b%time%dt)

         ! Solve Poisson equation
         call b%fs%update_laplacian()
         call b%fs%correct_mfr()
         call b%fs%get_div()
         call b%fs%add_surface_tension_jump(dt=b%time%dt,div=b%fs%div,vf=b%vf,contact_model=static_contact)
         b%fs%psolv%rhs=-b%fs%cfg%vol*b%fs%div/b%time%dt
         b%fs%psolv%sol=0.0_WP
         call b%fs%psolv%solve()
         call b%fs%shift_p(b%fs%psolv%sol)

         ! Correct velocity
         call b%fs%get_pgrad(b%fs%psolv%sol,b%resU,b%resV,b%resW)
         b%fs%P=b%fs%P+b%fs%psolv%sol
         b%fs%U=b%fs%U-b%time%dt*b%resU/b%fs%rho_U
         b%fs%V=b%fs%V-b%time%dt*b%resV/b%fs%rho_V
         b%fs%W=b%fs%W-b%time%dt*b%resW/b%fs%rho_W

         ! Increment sub-iteration counter
         b%time%it=b%time%it+1

      end do

      ! Recompute interpolated velocity and divergence
      call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
      call b%fs%get_div()

      ! End time steo timer
      endtime=mpi_wtime()

      ! Wall time spent in current time step
      b%cfg%step_wt=endtime-starttime

      ! Output to ensight
      if (b%ens_evt%occurs()) then 
         ! Update surfmesh object
         update_smesh: block
         use irl_fortran_interface
         integer :: nplane,np,i,j,k
         ! Transfer polygons to smesh
         call b%vf%update_surfmesh(b%smesh)
         ! Also populate nplane variable
         b%smesh%var(1,:)=1.0_WP
         np=0
         do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
            do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
               do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block update_smesh
      ! Perform ensight output 
         call b%ens_out%write_data(b%time%t)
      end if

      ! Update object time trackers
      call b%timer%vf_advance_timer(b%cfg,b%vf)
      call b%timer%sgs_visc_timer(b%cfg,b%sgs)
      call b%timer%implicit_timer(b%cfg,b%fs)
      call b%timer%pressure_timer(b%cfg,b%fs)
      call b%timer%step_timer(b%cfg)

      ! Perform and output monitoring
      call b%fs%get_max()
      call b%vf%get_max()
      call b%mfile%write()
      call b%cflfile%write()
      call b%timerfile%write()
      call b%timersummaryfile%write()

   end subroutine step


   !> Finalize b2 simulation
   subroutine final(b)
      implicit none
      class(block2), intent(inout) :: b

      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR)

   end subroutine final


end module block2_class
