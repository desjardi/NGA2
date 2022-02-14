!> Definition for a block 2 simulation: zoomed out cough machine
module block2_class
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   use config_class,      only: config
   use incomp_class,      only: incomp
   use hypre_uns_class,   only: hypre_uns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private
   
   
   public :: block2
   
   
   !> Block 2 object
   type :: block2
      class(config), pointer :: cfg       !< Pointer to config
      type(incomp) :: fs                  !< Single-phase incompressible flow solver
      type(hypre_uns) :: ps               !< Unstructured HYPRE pressure solver
      type(hypre_uns) :: is               !< Unstructured HYPRE implicit solver
      type(timetracker) :: time           !< Time tracker
      type(sgsmodel) ::  sgs              !< SGS model
      type(ensight) :: ens_out            !< Ensight output
      type(event) :: ens_evt              !< Ensight output event
      type(monitor) :: mfile,cflfile      !< Monitor files
      !> Private work arrays
      real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:), allocatable :: SR
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block2
   
   
   !> Gas viscosity
   real(WP) :: visc_g
   
   
contains
   
   
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
         call param_read('Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         call param_read('Max time',b%time%tmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker
      
      
      ! Create a single-phase flow solver with bconds
      create_solver: block
         use incomp_class,    only: dirichlet,clipped_neumann,neumann
         use hypre_uns_class, only: gmres_amg
         ! Create a single-phase flow solver
         b%fs=incomp(cfg=b%cfg,name='Single-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Gas dynamic viscosity',visc_g)
         b%fs%visc=visc_g
         ! Assign constant density to each phase
         call param_read('Gas density',b%fs%rho)
         ! Inflow on the left
         !call b%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_mouth)
         ! Outflow on the right
         !call b%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=right_boundary)
         ! Prepare and configure pressure solver
         b%ps=hypre_uns(cfg=b%cfg,name='Pressure',method=gmres_amg,nst=7)
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
         use tpns_class, only: bcond
         use random,     only: random_uniform
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! Zero initial field
         b%fs%U=0.0_WP; b%fs%V=0.0_WP; b%fs%W=0.0_WP
         ! Apply Dirichlet at inlet
         !call param_read('Gas velocity',Uin)
         !call param_read('Gas thickness',delta)
         !call param_read('Gas perturbation',Urand)
         !call b%fs%get_bcond('inflow',mybc)
         !do n=1,mybc%itr%no_
         !   i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !   b%fs%U(i,j,k)=Uin*tanh(2.0_WP*(0.5_WP*W_mouth-abs(b%fs%cfg%zm(k)))/delta)*tanh(2.0_WP*b%fs%cfg%ym(j)/delta)*tanh(2.0_WP*(H_mouth-b%fs%cfg%ym(j))/delta)+random_uniform(-Urand,Urand)
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
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(b%cfg,'cough_in')
         ! Create event for Ensight output
         b%ens_evt=event(b%time,'Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation')
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
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl')
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
   
   
   !> Take a time step with block 2
   subroutine step(b)
      use tpns_class, only: static_contact
      implicit none
      class(block2), intent(inout) :: b
      
      
      
   end subroutine step
   
   
   !> Finalize b1 simulation
   subroutine final(b)
      implicit none
      class(block2), intent(inout) :: b
      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR)
   end subroutine final
   
   
end module block2_class
