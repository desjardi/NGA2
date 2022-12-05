!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use lpt_class,         only: lpt
  use lowmach_class,     only: lowmach
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get an LPT solver, a lowmach solver, and corresponding time tracker
  type(lowmach),     public :: fs
  type(lpt),         public :: lp
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile
  type(monitor) :: lptfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays and fluid properties
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,dRHOdt
  real(WP) :: visc,rho,inlet_velocity

contains


   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none


    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max time',time%tmax)
      time%dt=time%dtmax
      time%itmax=2
    end block initialize_timetracker


    ! Create a low Mach flow solver with bconds
    create_flow_solver: block
      use ils_class,     only: pcg_pfmg,pcg_amg
      use lowmach_class, only: dirichlet,clipped_neumann
      ! Create flow solver
      fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
      ! Define boundary conditions
      call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
      call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )

      ! Assign constant density
      call param_read('Density',rho); fs%rho=rho
      ! Assign constant viscosity
      call param_read('Dynamic viscosity',visc); fs%visc=visc
      ! Assign acceleration of gravity
      call param_read('Gravity',fs%gravity)
      ! Configure pressure solver
      call param_read('Pressure iteration',fs%psolv%maxit)
      call param_read('Pressure tolerance',fs%psolv%rcvg)
      ! Configure implicit velocity solver
      call param_read('Implicit iteration',fs%implicit%maxit)
      call param_read('Implicit tolerance',fs%implicit%rcvg)
      ! Setup the solver
      call fs%setup(pressure_ils=pcg_amg,implicit_ils=pcg_pfmg)
    end block create_flow_solver


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(dRHOdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays


    ! Initialize our LPT solver
    initialize_lpt: block
      use random, only: random_uniform
      use mathtools, only: Pi
      real(WP) :: dp,Hbed,VFavg,Tp,Lpart,Lparty,Lpartz,Volp
      integer :: i,ix,iy,iz,np,npx,npy,npz
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Get particle diameter from the input
      call param_read('Particle diameter',dp)
      ! Set filter scale to 3.5*dx
      lp%filter_width=3.5_WP*cfg%min_meshsize

      ! Root process initializes particles uniformly
      call param_read('Bed height',Hbed)
      call param_read('Avg. volume fraction',VFavg)
      call param_read('Particle temperature',Tp,default=298.15_WP)
      if (lp%cfg%amRoot) then
         ! Particle volume
         Volp = Pi/6.0_WP*dp**3
         ! Get number of particles
         Lpart = (Volp/VFavg)**(1.0_WP/3.0_WP)
         npx=int(Hbed/Lpart)
         npy = int(cfg%yL/Lpart)
         Lparty = cfg%yL/real(npy,WP)
         npz = int(cfg%zL/Lpart)
         Lpartz = cfg%zL/real(npz,WP)
         np = npx*npy*npz
         call lp%resize(np)
         ! Distribute particles
         do i=1,np
            ! Give position
            ix = (i-1)/(npy*npz)
            iy = (i-1-npy*npz*ix)/npz
            iz = i-1-npy*npz*ix-npz*iy
            lp%p(i)%pos(1) = lp%cfg%x(lp%cfg%imin)+(real(ix,WP)+0.5_WP)*Lpart
            lp%p(i)%pos(2) = lp%cfg%y(lp%cfg%jmin)+(real(iy,WP)+0.5_WP)*Lparty
            lp%p(i)%pos(3) = lp%cfg%z(lp%cfg%kmin)+(real(iz,WP)+0.5_WP)*Lpartz

            ! Give id
            lp%p(i)%id=int(i,8)
            ! Set the diameter
            lp%p(i)%d=dp
            ! Set the temperature
            lp%p(i)%T=Tp
            ! Give zero velocity
            lp%p(i)%vel=0.0_WP
            ! Give zero collision force
            lp%p(i)%col=0.0_WP
            ! Give zero dt
            lp%p(i)%dt=0.0_WP
            ! Locate the particle on the mesh
            lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
            ! Activate the particle
            lp%p(i)%flag=0
         end do
      end if
      call lp%sync()

      ! Get initial particle volume fraction
      call lp%update_VF()
      ! Set collision timescale
      call param_read('Collision timescale',lp%Tcol,default=5.0_WP*time%dt)
      lp%Tcol=5.0_WP*time%dt
      ! Set coefficient of restitution
      call param_read('Coefficient of restitution',lp%e_n)
      ! Set gravity
      call param_read('Gravity',lp%gravity)
    end block initialize_lpt


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=1,name='lpt')
      pmesh%varname(1)='radius'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=0.5_WP*lp%p(i)%d
      end do
    end block create_pmesh


    ! Initialize our velocity field
    initialize_velocity: block
      use lowmach_class, only: bcond
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      ! Zero initial field
      fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
      ! Set inflow velocity/momentum
      call param_read('Inlet velocity',inlet_velocity)
      call fs%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         fs%rhoU(i,j,k)=inlet_velocity
         fs%U(i,j,k)   =inlet_velocity
      end do
      ! Set density from particle volume fraction
      fs%rho=rho*(1.0_WP-lp%VF)
      ! Form momentum
      call fs%rho_multiply
      ! Apply all other boundary conditions
      call fs%apply_bcond(time%t,time%dt)
      call fs%interp_vel(Ui,Vi,Wi)
      call fs%get_div(drhodt=dRHOdt)
      ! Compute MFR through all boundary conditions
      call fs%get_mfr()
    end block initialize_velocity


    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='fluidized_bed')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('epsp',lp%VF)
      call ens_out%add_scalar('pressure',fs%P)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight


    ! Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_max()
      call lp%get_max()
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
      call mfile%add_column(fs%divmax,'Maximum divergence')
      call mfile%add_column(fs%psolv%it,'Pressure iteration')
      call mfile%add_column(fs%psolv%rerr,'Pressure error')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(fs%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
      call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
      call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
      call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
      call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
      call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
      call cflfile%write()
      ! Create LPT monitor
      lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
      call lptfile%add_column(time%n,'Timestep number')
      call lptfile%add_column(time%t,'Time')
      call lptfile%add_column(lp%VFmean,'VFp mean')
      call lptfile%add_column(lp%VFmax,'VFp max')
      call lptfile%add_column(lp%np,'Particle number')
      call lptfile%add_column(lp%Umin,'Particle Umin')
      call lptfile%add_column(lp%Umax,'Particle Umax')
      call lptfile%add_column(lp%Vmin,'Particle Vmin')
      call lptfile%add_column(lp%Vmax,'Particle Vmax')
      call lptfile%add_column(lp%Wmin,'Particle Wmin')
      call lptfile%add_column(lp%Wmax,'Particle Wmax')
      call lptfile%add_column(lp%dmin,'Particle dmin')
      call lptfile%add_column(lp%dmax,'Particle dmax')
      call lptfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    use mathtools, only: twoPi
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call fs%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Remember old density, velocity, and momentum
       fs%rhoold=fs%rho
       fs%Uold=fs%U; fs%rhoUold=fs%rhoU
       fs%Vold=fs%V; fs%rhoVold=fs%rhoV
       fs%Wold=fs%W; fs%rhoWold=fs%rhoW

       ! Collide and advance particles
       call lp%collide(dt=time%dtmid)
       resU=rho
       call lp%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=fs%visc)

       ! Update density based on particle volume fraction
       fs%rho=rho*(1.0_WP-lp%VF)
       dRHOdt=(fs%RHO-fs%RHOold)/time%dtmid

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)

          ! Build mid-time velocity and momentum
          fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
          fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
          fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

          ! Explicit calculation of drho*u/dt from NS
          call fs%get_dmomdt(resU,resV,resW)

          ! Add momentum source terms
          call fs%addsrc_gravity(resU,resV,resW)

          ! Assemble explicit residual
          resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
          resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
          resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

          ! Add momentum source term from lpt
          add_lpt_src: block
            integer :: i,j,k
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*lp%srcU(i-1:i,j,k))
                     resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*lp%srcV(i,j-1:j,k))
                     resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*lp%srcW(i,j,k-1:k))
                  end do
               end do
            end do
          end block add_lpt_src

          ! Form implicit residuals
          call fs%solve_implicit(time%dtmid,resU,resV,resW)

          ! Apply these residuals
          fs%U=2.0_WP*fs%U-fs%Uold+resU
          fs%V=2.0_WP*fs%V-fs%Vold+resV
          fs%W=2.0_WP*fs%W-fs%Wold+resW

          ! Apply other boundary conditions and update momentum
          call fs%apply_bcond(time%tmid,time%dtmid)
          call fs%rho_multiply()
          call fs%apply_bcond(time%tmid,time%dtmid)

          ! Reset Dirichlet BCs
          dirichlet_velocity: block
            use lowmach_class, only: bcond
            type(bcond), pointer :: mybc
            integer :: n,i,j,k
            call fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs%rhoU(i,j,k)=inlet_velocity
               fs%U(i,j,k)   =inlet_velocity
            end do
          end block dirichlet_velocity

          ! Solve Poisson equation
          call fs%correct_mfr(drhodt=dRHOdt)
          call fs%get_div(drhodt=dRHOdt)
          fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
          fs%psolv%sol=0.0_WP
          call fs%psolv%solve()
          call fs%shift_p(fs%psolv%sol)

          ! Correct momentum and rebuild velocity
          call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
          fs%P=fs%P+fs%psolv%sol
          fs%rhoU=fs%rhoU-time%dtmid*resU
          fs%rhoV=fs%rhoV-time%dtmid*resV
          fs%rhoW=fs%rhoW-time%dtmid*resW
          call fs%rho_divide

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute interpolated velocity and divergence
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div(drhodt=dRHOdt)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            do i=1,lp%np_
               pmesh%var(1,i)=0.5_WP*lp%p(i)%d
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call fs%get_max()
       call lp%get_max()
       call mfile%write()
       call cflfile%write()
       call lptfile%write()

    end do

    ! Ouput particle bed
    !call lp%write(filename='part.file')

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
    deallocate(resU,resV,resW,Ui,Vi,Wi,dRHOdt)

  end subroutine simulation_final

end module simulation
