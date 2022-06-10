!> Definition for a block 1 simulation
module block1_class
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: t_wall,L_mouth,H_mouth,W_mouth
   use config_class,      only: config
   use incomp_class,      only: incomp
   use hypre_uns_class,   only: hypre_uns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private

   public :: block1


   !> Block 1 object
   type :: block1
      class(config), pointer :: cfg                               !< Pointer to config
      type(incomp) :: fs                                          !< Single-phase incompressible flow solver
      type(hypre_uns) :: ps                                       !< Unstructured HYPRE pressure solver
      type(hypre_uns) :: is                                       !< Unstructured HYPRE implicit solver
      type(timetracker) :: time                                   !< Time tracker
      type(sgsmodel) ::  sgs                                      !< SGS model
      type(ensight) :: ens_out                                    !< Ensight output
      type(event) :: ens_evt                                      !< Ensight output event
      type(monitor) :: mfile,cflfile                              !< Monitor files

      !> Private work arrays
      real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:), allocatable :: SR
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block1

   !> Gas viscosity
   real(WP) :: visc_g
   
   !> Gas Density
   real(WP) :: rho_g
   
   !> Inflow parameters
   real(WP) :: Uin,delta,Urand,Uco,CPFR
   
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


   !> Function that calcuates velocity at current time 
   function inflowVelocity(time,CPFR,H,W) result(UCPFR)
      real(WP), intent(in)   :: time,CPFR,H,W
      real(WP)               :: UCPFR
      real(WP)               :: CEV,PVT,a1,b1,c1,a2,b2,c2,tau,M
      class(config), pointer :: cfg 
      !Model parameters
      CEV=0.20_WP*CPFR-4e-5_WP    !cough expiratory volume
      PVT=2.85_WP*CPFR+0.07_WP    !Peak velocity time
      a1=1.68_WP
      b1=3.34_WP
      c1=0.43_WP
      a2=(CEV/(PVT*CPFR))-a1
      b2=((-2.16_WP*CEV)/(PVT*CPFR))+10.46_WP
      c2=(1.8_WP/(b2-1.0_WP))
      !Dimensionless time
      tau=time/PVT
      !Dimensionless flow rate
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
   
   
   !> Initialization of block 1
   subroutine init(b)
      use param, only: param_read
      implicit none
      class(block1), intent(inout) :: b

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
         b%time=timetracker(b%cfg%amRoot,name='duct_flow')
         call param_read('Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         call param_read('Max time',b%time%tmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker

      ! Create a single-phase flow solver with bconds
      create_solver: block
         use incomp_class,    only: dirichlet,clipped_neumann,neumann
         use hypre_uns_class, only: pcg_amg,gmres_amg
         ! Create a single-phase flow solver
         b%fs=incomp(cfg=b%cfg,name='Single-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Gas dynamic viscosity',visc_g)
         b%fs%visc=visc_g
         ! Assign constant density to each phase
         call param_read('Gas density',b%fs%rho)
         ! Inflow on the left
         call b%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_mouth)
         call b%fs%add_bcond(name='coflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_coflow)
         ! Apply clipped Neumann on the right
         call b%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary)
         ! Prepare and configure pressure solver
         b%ps=hypre_uns(cfg=b%cfg,name='Pressure',method=pcg_amg,nst=7)
         b%ps%maxlevel=22
         call param_read('Pressure iteration',b%ps%maxit)
         call param_read('Pressure tolerance',b%ps%rcvg)
         ! Prepare and configure implicit solver
         b%is=hypre_uns(cfg=b%cfg,name='Implicit',method=gmres_amg,nst=7)
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
         !call param_read('Gas velocity',Uin)
         call param_read('Peak flow rate',CPFR)
         Uin=inflowVelocity(b%time%t,CPFR,H_mouth,W_mouth)
         call param_read('Gas thickness',delta)
         call param_read('Gas perturbation',Urand)
         call b%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            !b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)
            b%fs%U(i,j,k)=Uin
         end do
         !Apply coflow around inlet geometry
         !call param_read('Gas coflow',Uco)
         Uco=0.10_WP*Uin
         call b%fs%get_bcond('coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            b%fs%U(i,j,k)=Uco
         end do
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

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(b%cfg,'duct_flow')
         ! Create event for Ensight output
         b%ens_evt=event(b%time,'Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call b%ens_out%add_scalar('visc_t',b%sgs%visc)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation1')
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
         call b%mfile%add_column(Uin,'Inflow Velocity')
         call b%mfile%write()
         ! Create CFL monitor
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl1')
         call b%cflfile%add_column(b%time%n,'Timestep number')
         call b%cflfile%add_column(b%time%t,'Time')
         call b%cflfile%add_column(b%fs%CFLc_x,'Convective xCFL')
         call b%cflfile%add_column(b%fs%CFLc_y,'Convective yCFL')
         call b%cflfile%add_column(b%fs%CFLc_z,'Convective zCFL')
         call b%cflfile%add_column(b%fs%CFLv_x,'Viscous xCFL')
         call b%cflfile%add_column(b%fs%CFLv_y,'Viscous yCFL')
         call b%cflfile%add_column(b%fs%CFLv_z,'Viscous zCFL')
         call b%cflfile%write()
      end block create_monitor


   end subroutine init


   !> Take a time step with block 1
   subroutine step(b)
      implicit none
      class(block1), intent(inout) :: b

      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! Apply time-varying Dirichlet conditions
      reapply_dirichlet: block
         use incomp_class, only: bcond
         use random,     only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! Reapply Dirichlet at inlet
         Uin=inflowVelocity(b%time%t,CPFR,H_mouth,W_mouth)
         call b%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            !b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
            b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)
            b%fs%U(i,j,k)=Uin
         end do
         ! Reapply coflow around inlet geometry
         Uco=0.10_WP*Uin
         call b%fs%get_bcond('coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            b%fs%U(i,j,k)=Uco
         end do
      end block reapply_dirichlet

      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W

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

         ! Form implicit residuals
         !if (b%cfg%rank.eq.1) print *, "B2 - pre solve_implicit, time step: ", b%time%n," iterations: ",b%time%it
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)
         !if (b%cfg%rank.eq.1) print *, "B2 - post solve_implicit, time step: ", b%time%n," iterations: ",b%time%it

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
         !if (b%cfg%rank.eq.1) print *, "B2 - pre ps%sol, time step: ", b%time%n," iterations: ",b%time%it
         call b%ps%solve()
         !if (b%cfg%rank.eq.1) print *, "B2 - post ps%sol, time step: ", b%time%n," iterations: ",b%time%it

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
         ! Perform ensight output
         call b%ens_out%write_data(b%time%t)
      end if

      ! Perform and output monitoring
      call b%fs%get_max()
      call b%mfile%write()
      call b%cflfile%write()

   end subroutine step


   !> Finalize b1 simulation
   subroutine final(b)
      implicit none
      class(block1), intent(inout) :: b

      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR)

   end subroutine final


end module block1_class
