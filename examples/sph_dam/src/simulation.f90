!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use sph_class,         only: sph
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private
  
  !> Get an SPH solver and corresponding time tracker
  type(sph),         public :: ss
  type(timetracker), public :: time
  
  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt
  
  !> Simulation monitor file
  type(monitor) :: mfile,cflfile
  
  public :: simulation_init,simulation_run,simulation_final

contains


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none


    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max time',time%tmax)
      call param_read('Max cfl number',time%cflmax)
      time%dt=time%dtmax
    end block initialize_timetracker


    ! Initialize SPH solver
    initialize_sph: block
      use random, only: random_uniform
      use mathtools, only: twoPi
      real(WP) :: Hbed,Lpx,Lpy,dx,dy
      integer :: i,j,n,np,npx,npy
      ! Create solver
      ss=sph(cfg=cfg,name='SPH')
      ! Set gravity
      call param_read('Gravity',ss%gravity)
      ! Set fluid properties
      call param_read('Viscosity',ss%visc)
      call param_read('Density',ss%rho)
      call param_read('Gamma',ss%gamm)
      call param_read('Min sound speed',ss%cmin)
      ! Set kernel properties
      call param_read('Smoothing length',ss%h)
      call param_read('Cutoff radius',ss%rmax)
      ! Particle configuration
      call param_read('Lpx',Lpx)
      call param_read('Lpy',Lpy)
      call param_read('Npx',npx); dx=Lpx/real(npx,WP)
      call param_read('Npy',npy); dy=Lpy/real(npy,WP)
      ! Root process distributes particles
      if (ss%cfg%amRoot) then
         np = npx*npy
         call ss%resize(np)
         ! Distribute particles
         n=0
         do j=1,npy
            do i=1,npx
               n=n+1
               ! Give position
               ss%p(n)%pos(1)=real(i,WP)*dx-0.5_WP*Lpx
               ss%p(n)%pos(2)=ss%cfg%y(ss%cfg%jmin)+real(j,WP)*dy
               ss%p(n)%pos(3)=0.0_WP
               ! Give id
               ss%p(n)%id=n
               ! Set the volume
               ss%p(n)%vol=dx*dy
               ! Give zero velocity
               ss%p(n)%vel=0.0_WP
               ! Give zero RHS
               ss%p(n)%dxdt=0.0_WP
               ss%p(n)%dudt=0.0_WP
               ss%p(n)%dRhoDt=0.0_WP
               ! Set density
               ss%p(n)%rho=ss%rho
               ! Locate the particle on the mesh
               ss%p(n)%ind=ss%cfg%get_ijk_global(ss%p(n)%pos,[ss%cfg%imin,ss%cfg%jmin,ss%cfg%kmin])
               ! Activate the particle
               ss%p(n)%flag=0
            end do
         end do
      end if
      call ss%sync()

      ! Initialize solver now that parameters are set
      call ss%sph_init()

    end block initialize_sph


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=3,nvec=1,name='lpt')
      pmesh%varname(1)='diameter'
      pmesh%varname(2)='density'
      pmesh%varname(3)='pressure'
      pmesh%vecname(1)='velocity'
      call ss%update_partmesh(pmesh)
      do i=1,ss%np_
         pmesh%var(1,i)=ss%p(i)%dp
         pmesh%var(2,i)=ss%p(i)%rho
         pmesh%var(3,i)=ss%p(i)%P
         pmesh%vec(:,1,i)=ss%p(i)%vel
      end do
    end block create_pmesh
    
    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='sph')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_particle('particles',pmesh)
      call ens_out%add_scalar('Wdist',ss%wdist)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight

    
    ! Create monitor filea
    create_monitor: block
      ! Prepare some info about fields
      call ss%get_cfl(time%dt,time%cfl)
      call ss%get_max()
      ! Create simulation monitor
      mfile=monitor(ss%cfg%amRoot,'simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(time%cfl,'Maximum CFL')
      call mfile%add_column(ss%np,'Particle number')
      call mfile%add_column(ss%Umax,'Umax')
      call mfile%add_column(ss%Vmax,'Vmax')
      call mfile%add_column(ss%Wmax,'Wmax')
      call mfile%add_column(ss%RHOmax,'RHOmax')
      call mfile%add_column(ss%dRHO,'dRHO (%)')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(ss%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(ss%CFLp_x,'xCFL')
      call cflfile%add_column(ss%CFLp_y,'yCFL')
      call cflfile%add_column(ss%CFLp_z,'zCFL')
      call cflfile%add_column(ss%CFLp_c,'cCFL')
      call cflfile%add_column(ss%CFLp_f,'fCFL')
      call cflfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call ss%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Advance SPH solver
       call ss%advance(dt=time%dtmid)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call ss%update_partmesh(pmesh)
            do i=1,ss%np_
               pmesh%var(1,i)=ss%p(i)%dp
               pmesh%var(2,i)=ss%p(i)%rho
               pmesh%var(3,i)=ss%p(i)%P
               pmesh%vec(:,1,i)=ss%p(i)%vel
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call ss%get_max()
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
    ! timetracker

    ! Deallocate work arrays

  end subroutine simulation_final

end module simulation
