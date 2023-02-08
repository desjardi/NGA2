!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use lpt_class,         only: lpt
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Only get a LPT solver and corresponding time tracker
  type(lpt),         public :: lp
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile

  public :: simulation_init,simulation_run,simulation_final

  !> Fluid phase arrays
  real(WP), dimension(:,:,:), allocatable :: U,V,W
  real(WP), dimension(:,:,:), allocatable :: rho,visc

contains


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none


    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max cfl number',time%cflmax)
      time%dt=time%dtmax
      time%itmax=1
    end block initialize_timetracker


    ! Initialize our LPT
    initialize_lpt: block
      use random, only: random_uniform
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Set gravity
      call param_read('Gravity',lp%gravity)
      ! Set filter scale to 3.5*dx
      lp%filter_width=3.5_WP*cfg%min_meshsize
      ! Turn off drag
      lp%drag_model='none'
      ! Initialize with zero particles
      call lp%resize(0)
      ! Get initial particle volume fraction
      call lp%update_VF()
      ! Set collision timescale
      lp%tau_col=5.0_WP*time%dt
      ! Set coefficient of restitution
      call param_read('Coefficient of restitution',lp%e_n)
      call param_read('Wall restitution',lp%e_w)
      call param_read('Friction coefficient',lp%mu_f)
      ! Injection parameters
      call param_read('Particle mass flow rate',lp%mfr)
      call param_read('Particle velocity',lp%inj_vel)
      call param_read('Particle mean diameter',lp%inj_dmean)
      call param_read('Particle standard deviation',lp%inj_dsd,default=0.0_WP)
      call param_read('Particle min diameter',lp%inj_dmin,default=tiny(1.0_WP))
      call param_read('Particle max diameter',lp%inj_dmax,default=huge(1.0_WP))
      call param_read('Particle diameter shift',lp%inj_dshift,default=0.0_WP)
      if (lp%inj_dsd.le.epsilon(1.0_WP)) then
         lp%inj_dmin=lp%inj_dmean
         lp%inj_dmax=lp%inj_dmean
      end if
      call param_read('Particle inject diameter',lp%inj_d)
      lp%inj_pos(1)=lp%cfg%x(lp%cfg%imin)+lp%inj_dmax
      lp%inj_pos(2:3)=0.0_WP
      lp%inj_T=300.0_WP
    end block initialize_lpt


    ! Initialize quiescent fluid
    initialize_fluid: block
      use mathtools, only: twoPi
      integer :: i,j,k
      real(WP) :: rhof,viscf
      ! Allocate arrays
      allocate(rho (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
      allocate(visc(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
      allocate(U(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
      allocate(V(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
      allocate(W(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
      U=0.0_WP; V=0.0_WP; W=0.0_WP
      ! Set constant density and viscosity
      call param_read('Density',rhof); rho=rhof
      call param_read('Viscosity',viscf); visc=viscf
    end block initialize_fluid


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=2,nvec=2,name='lpt')
      pmesh%varname(1)='id'
      pmesh%varname(2)='diameter'
      pmesh%vecname(1)='velocity'
      pmesh%vecname(2)='ang_vel'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=real(lp%p(i)%id,WP)
         pmesh%var(2,i)=lp%p(i)%d
         pmesh%vec(:,1,i)=lp%p(i)%vel
         pmesh%vec(:,2,i)=lp%p(i)%angVel
      end do
    end block create_pmesh


    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=lp%cfg,name='particles')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',U,V,W)
      call ens_out%add_scalar('epsp',lp%VF)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
      call lp%get_max()
      ! Create simulation monitor
      mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(lp%np,'Particle number')
      call mfile%add_column(lp%np_new,'Npart new')
      call mfile%add_column(lp%np_out,'Npart removed')
      call mfile%add_column(lp%ncol,'Particle collisions')
      call mfile%add_column(lp%VFmax,'Max VF')
      call mfile%add_column(lp%Umin,'Particle Umin')
      call mfile%add_column(lp%Umax,'Particle Umax')
      call mfile%add_column(lp%Vmin,'Particle Vmin')
      call mfile%add_column(lp%Vmax,'Particle Vmax')
      call mfile%add_column(lp%Wmin,'Particle Wmin')
      call mfile%add_column(lp%Wmax,'Particle Wmax')
      call mfile%add_column(lp%dmin,'Particle dmin')
      call mfile%add_column(lp%dmax,'Particle dmax')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(lp%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(lp%CFLp_x,'Particle xCFL')
      call cflfile%add_column(lp%CFLp_y,'Particle yCFL')
      call cflfile%add_column(lp%CFLp_z,'Particle zCFL')
      call cflfile%add_column(lp%CFL_col,'Collision CFL')
      call cflfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Inject particles
       call lp%inject(dt=time%dt,avoid_overlap=.true.)

       ! Collide particles
       call lp%collide(dt=time%dt)

       ! Advance particles by dt
       U=0.0_WP; V=0.0_WP; W=0.0_WP
       call lp%advance(dt=time%dt,U=U,V=V,W=W,rho=rho,visc=visc,stress_x=U,stress_y=V,stress_z=W)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            do i=1,lp%np_
               pmesh%var(1,i)=real(lp%p(i)%id,WP)
               pmesh%var(2,i)=lp%p(i)%d
               pmesh%vec(:,1,i)=lp%p(i)%vel
               pmesh%vec(:,2,i)=lp%p(i)%angVel
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call lp%get_max()
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
    deallocate(rho,visc,U,V,W)

  end subroutine simulation_final

end module simulation
