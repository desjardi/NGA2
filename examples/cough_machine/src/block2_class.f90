!> Definition for a block 2 simulation: zoomed out cough machine
module block2_class
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   use config_class,      only: config
   use incomp_class,      only: incomp
   use hypre_str_class,   only: hypre_str
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

   public :: block2


   !> Block 2 object
   type :: block2
      class(config), pointer :: cfg                 !< Pointer to config
      type(incomp) :: fs                            !< Single-phase incompressible flow solver
      type(hypre_str) :: ps                         !< Unstructured HYPRE pressure solver
      type(hypre_str) :: is                         !< Unstructured HYPRE implicit solver
      type(lpt)       :: lp                         !< Lagrangian particle tracking
      type(partmesh)  :: pmesh                      !< Partmesh for Lagrangian particle output
      type(timetracker) :: time                     !< Time tracker
      type(sgsmodel) ::  sgs                        !< SGS model
      type(ensight) :: ens_out                      !< Ensight output
      type(event) :: ens_evt                        !< Ensight output event
      type(monitor) :: mfile,cflfile,sprayfile      !< Monitor files
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


   !> Gas viscosity
   real(WP) :: visc_g


contains

   !> Function that localizes the left domain boundary
   function left_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_boundary

   !> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary

   !> Function that localizes cough stream at L_mouth
   function gas_inj_spray(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.L_mouth.and.pg%ym(j).gt.0.0_WP.and.pg%ym(j).lt.H_mouth.and.abs(pg%zm(k)).lt.0.5_WP*W_mouth) isIn=.true.
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


   !> Initialization of block 2
   subroutine init(b)
      use param, only: param_read
      implicit none
      class(block2), intent(inout) :: b


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
         b%time=timetracker(b%cfg%amRoot,name='cough_machine_out')
         call param_read('2 Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         call param_read('Max time',b%time%tmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker


      ! Create a single-phase flow solver with bconds
      create_solver: block
         use incomp_class,    only: dirichlet,clipped_neumann,neumann
         use hypre_str_class, only: pcg_pfmg
         ! Create a single-phase flow solver
         b%fs=incomp(cfg=b%cfg,name='Single-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Gas dynamic viscosity',visc_g)
         b%fs%visc=visc_g
         ! Assign constant density to each phase
         call param_read('Gas density',b%fs%rho)
         ! Define direction gas/liquid stream boundary conditions (inflow from block1 to block2)
         !call b%fs%add_bcond(name='gas_inj',type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=gas_inj_spray)
         ! Neumann on block 1 interface
         call b%fs%add_bcond(name='b1_interface'  ,type=neumann,face='x',dir=-1,canCorrect=.true. ,locator=left_boundary)
         ! Neumann on the sides
         call b%fs%add_bcond(name='bc_yp'  ,type=neumann,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call b%fs%add_bcond(name='bc_ym'  ,type=neumann,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         call b%fs%add_bcond(name='bc_zp'  ,type=neumann,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call b%fs%add_bcond(name='bc_zm'  ,type=neumann,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         ! Outflow on the right
         call b%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary)
         ! Prepare and configure pressure solver
         b%ps=hypre_str(cfg=b%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         b%ps%maxlevel=20
         call param_read('Pressure iteration',b%ps%maxit)
         call param_read('Pressure tolerance',b%ps%rcvg)
         ! Prepare and configure implicit solver
         b%is=hypre_str(cfg=b%cfg,name='Implicit',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',b%is%maxit)
         call param_read('Implicit tolerance',b%is%rcvg)
         ! Setup the solver
         call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%is)
      end block create_solver


      ! Initialize our velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         use random,     only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! Zero initial field
         b%fs%U=0.0_WP; b%fs%V=0.0_WP; b%fs%W=0.0_WP
         ! Apply Dirichlet at inlet
         !call b%fs%get_bcond('gas_inj',mybc)
         !do n=1,mybc%itr%no_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
         !   ! U velocity in x/i
         !   ! V velocity in y/j
         !   ! W velocity in z/k
         !end do
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


      ! Create an LES model
      create_sgs: block
         b%sgs=sgsmodel(cfg=b%fs%cfg,umask=b%fs%umask,vmask=b%fs%vmask,wmask=b%fs%wmask)
      end block create_sgs

      ! Initialize our Lagrangian spray solver
      initialize_lpt: block
         ! Create solver
         b%lp=lpt(cfg=b%cfg,name='spray')
         ! Get droplet density from the input
         call param_read('Liquid density',b%lp%rho)
      end block initialize_lpt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         ! Include an extra variable for droplet diameter
         b%pmesh=partmesh(nvar=1,name='lpt')
         b%pmesh%varname(1)='diameter'
         ! Transfer particles to pmesh
         call b%lp%update_partmesh(b%pmesh)
         ! Also populate diameter variable
         do i=1,b%lp%np_
            b%pmesh%var(1,i)=b%lp%p(i)%d
         end do
      end block create_pmesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(b%cfg,'cough_out')
         ! Create event for Ensight output
         b%ens_evt=event(b%time,'Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call b%ens_out%add_scalar('visc_t',b%sgs%visc)
         call b%ens_out%add_particle('spray',b%pmesh)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
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
         call b%mfile%add_column(b%fs%divmax,'Maximum divergence')
         call b%mfile%add_column(b%fs%psolv%it,'Pressure iteration')
         call b%mfile%add_column(b%fs%psolv%rerr,'Pressure error')
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
         ! Create a spray monitor
         b%sprayfile=monitor(amroot=b%lp%cfg%amRoot,name='spray')
         call b%sprayfile%add_column(b%time%n,'Timestep number')
         call b%sprayfile%add_column(b%time%t,'Time')
         call b%sprayfile%add_column(b%time%dt,'Timestep size')
         call b%sprayfile%add_column(b%lp%np,'Droplet number')
         call b%sprayfile%add_column(b%lp%Umin, 'Umin')
         call b%sprayfile%add_column(b%lp%Umax, 'Umax')
         call b%sprayfile%add_column(b%lp%Umean,'Umean')
         call b%sprayfile%add_column(b%lp%Vmin, 'Vmin')
         call b%sprayfile%add_column(b%lp%Vmax, 'Vmax')
         call b%sprayfile%add_column(b%lp%Vmean,'Vmean')
         call b%sprayfile%add_column(b%lp%Wmin, 'Wmin')
         call b%sprayfile%add_column(b%lp%Wmax, 'Wmax')
         call b%sprayfile%add_column(b%lp%Wmean,'Wmean')
         call b%sprayfile%add_column(b%lp%dmin, 'dmin')
         call b%sprayfile%add_column(b%lp%dmax, 'dmax')
         call b%sprayfile%add_column(b%lp%dmean,'dmean')
         call b%sprayfile%write()
      end block create_monitor


   end subroutine init


   !> Take a time step with block 2
   subroutine step(b,Unudge,Vnudge,Wnudge)
      implicit none
      class(block2), intent(inout) :: b
      real(WP), dimension(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%jmaxo_), intent(inout) :: Unudge     !< Updated to (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%jmaxo_), intent(inout) :: Vnudge     !< Updated to (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%jmaxo_), intent(inout) :: Wnudge     !< Updated to (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)

      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! Advance particles by full dt
      b%resU=b%fs%rho; b%resV=b%fs%visc
      call b%lp%advance(dt=b%time%dt,U=b%fs%U,V=b%fs%V,W=b%fs%W,rho=b%resU,visc=b%resV)

      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W

      ! Apply time-varying Dirichlet conditions
      ! This is where time-dpt Dirichlet would be enforced

      ! Reset here gas viscosity
      b%fs%visc=visc_g

      ! Turbulence modeling
      call b%fs%get_strainrate(Ui=b%Ui,Vi=b%Vi,Wi=b%Wi,SR=b%SR)
      b%resU=b%fs%rho
      call b%sgs%get_visc(dt=b%time%dtold,rho=b%resU,Ui=b%Ui,Vi=b%Vi,Wi=b%Wi,SR=b%SR)
      where (b%sgs%visc.lt.-b%fs%visc)
         b%sgs%visc=-b%fs%visc
      end where
      b%fs%visc=b%fs%visc+b%sgs%visc

      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)

         ! Build mid-time velocity
         b%fs%U=0.5_WP*(b%fs%U+b%fs%Uold)
         b%fs%V=0.5_WP*(b%fs%V+b%fs%Vold)
         b%fs%W=0.5_WP*(b%fs%W+b%fs%Wold)

         ! Explicit calculation of drho*u/dt from NS
         call b%fs%get_dmomdt(b%resU,b%resV,b%resW)

         ! Assemble explicit residual
         b%resU=-2.0_WP*(b%fs%rho*b%fs%U-b%fs%rho*b%fs%Uold)+b%time%dt*b%resU
         b%resV=-2.0_WP*(b%fs%rho*b%fs%V-b%fs%rho*b%fs%Vold)+b%time%dt*b%resV
         b%resW=-2.0_WP*(b%fs%rho*b%fs%W-b%fs%rho*b%fs%Wold)+b%time%dt*b%resW

         ! Add nudging term here
         nudge: block
            real(WP), dimension(3) :: pos
            real(WP) :: dist
            integer :: i,j,k
            do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
               do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
                  do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                     if (b%fs%umask(i,j,k).eq.0) then
                        pos=[b%fs%cfg%x(i),b%fs%cfg%ym(j),b%fs%cfg%zm(k)]
                        dist=min(b%nudge_xmax-pos(1),pos(1)-b%nudge_xmin,&
                        &        b%nudge_ymax-pos(2),pos(2)-b%nudge_ymin,&
                        &        b%nudge_zmax-pos(3),pos(3)-b%nudge_zmin)/b%nudge_trans
                        dist=min(max(dist,0.0_WP),1.0_WP)
                        b%resU(i,j,k)=b%resU(i,j,k)+(Unudge(i,j,k)-b%fs%U(i,j,k))*dist**2
                     end if
                     if (b%fs%vmask(i,j,k).eq.0) then
                        pos=[b%fs%cfg%xm(i),b%fs%cfg%y(j),b%fs%cfg%zm(k)]
                        dist=min(b%nudge_xmax-pos(1),pos(1)-b%nudge_xmin,&
                        &        b%nudge_ymax-pos(2),pos(2)-b%nudge_ymin,&
                        &        b%nudge_zmax-pos(3),pos(3)-b%nudge_zmin)/b%nudge_trans
                        dist=min(max(dist,0.0_WP),1.0_WP)
                        b%resV(i,j,k)=b%resV(i,j,k)+(Vnudge(i,j,k)-b%fs%V(i,j,k))*dist**2
                     end if
                     if (b%fs%wmask(i,j,k).eq.0) then
                        pos=[b%fs%cfg%xm(i),b%fs%cfg%ym(j),b%fs%cfg%z(k)]
                        dist=min(b%nudge_xmax-pos(1),pos(1)-b%nudge_xmin,&
                        &        b%nudge_ymax-pos(2),pos(2)-b%nudge_ymin,&
                        &        b%nudge_zmax-pos(3),pos(3)-b%nudge_zmin)/b%nudge_trans
                        dist=min(max(dist,0.0_WP),1.0_WP)
                        b%resW(i,j,k)=b%resW(i,j,k)+(Wnudge(i,j,k)-b%fs%W(i,j,k))*dist**2
                     end if
                  end do
               end do
            end do
         end block nudge

         ! Form implicit residuals
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)

         ! Apply these residuals
         b%fs%U=2.0_WP*b%fs%U-b%fs%Uold+b%resU
         b%fs%V=2.0_WP*b%fs%V-b%fs%Vold+b%resV
         b%fs%W=2.0_WP*b%fs%W-b%fs%Wold+b%resW

         ! Apply other boundary conditions on the resulting fields
         call b%fs%apply_bcond(b%time%t,b%time%dt)

         ! Solve Poisson equation
         call b%fs%correct_mfr()
         call b%fs%get_div()
         b%ps%rhs=-b%fs%cfg%vol*b%fs%div*b%fs%rho/b%time%dt
         b%ps%sol=0.0_WP
         call b%ps%solve()

         ! Correct velocity
         call b%fs%get_pgrad(b%ps%sol,b%resU,b%resV,b%resW)
         b%fs%P=b%fs%P+b%ps%sol
         b%fs%U=b%fs%U-b%time%dt*b%resU/b%fs%rho
         b%fs%V=b%fs%V-b%time%dt*b%resV/b%fs%rho
         b%fs%W=b%fs%W-b%time%dt*b%resW/b%fs%rho

         ! Increment sub-iteration counter
         b%time%it=b%time%it+1

      end do

      ! Recompute interpolated velocity and divergence
      call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
      call b%fs%get_div()

      ! Output to ensight
      if (b%ens_evt%occurs()) then
         ! Update partmesh object
         update_pmesh: block
            integer :: i
            ! Transfer particles to pmesh
            call b%lp%update_partmesh(b%pmesh)
            ! Also populate diameter variable
            do i=1,b%lp%np_
               b%pmesh%var(1,i)=b%lp%p(i)%d
            end do
         end block update_pmesh
         ! Perform ensight output
         call b%ens_out%write_data(b%time%t)
      end if

      ! Perform and output monitoring
      call b%fs%get_max()
      call b%lp%get_max()
      call b%mfile%write()
      call b%cflfile%write()
      call b%sprayfile%write()

   end subroutine step


   !> Finalize b1 simulation
   subroutine final(b)
      implicit none
      class(block2), intent(inout) :: b

      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR)

   end subroutine final


end module block2_class
