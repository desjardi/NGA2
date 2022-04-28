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
   type(monitor) :: mfile
   
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
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker
      

      ! Initialize our LPT
      initialize_lpt: block
         use random, only: random_uniform
         real(WP) :: dp
         integer :: i
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Root process initializes 1000 particles randomly
         if (lp%cfg%amRoot) then
            call lp%resize(1000)
            do i=1,1000
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! Assign random position in the domain
               lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
               &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
               &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
               lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin)
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
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         lp%Tcol=5.0_WP*time%dt
      end block initialize_lpt
      
      
      ! Test particle I/O
      !test_lpt_io: block
      !   ! Write it out
      !   call lp%write(filename='part.file')
      !   ! Read it back up
      !   call lp%read(filename='part.file')
      !end block test_lpt_io
      
      
      ! Prepare our fluid phase info based on a Taylor vortex
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
         ! Initialize to solid body rotation
         do k=lp%cfg%kmino_,lp%cfg%kmaxo_
            do j=lp%cfg%jmino_,lp%cfg%jmaxo_
               do i=lp%cfg%imino_,lp%cfg%imaxo_
                  ! Solid body rotation
                  !U(i,j,k)=-twoPi*lp%cfg%ym(j)
                  !V(i,j,k)=+twoPi*lp%cfg%xm(i)
                  !W(i,j,k)=0.0_WP
                  ! Taylor-Green vortex
                  U(i,j,k)=+cos(twoPi*lp%cfg%x (i))*sin(twoPi*lp%cfg%ym(j))
                  V(i,j,k)=-sin(twoPi*lp%cfg%xm(i))*cos(twoPi*lp%cfg%y (j))
                  W(i,j,k)=0.0_WP
               end do
            end do
         end do
         ! Set constant density and viscosity
         call param_read('Density',rhof); rho=rhof
         call param_read('Viscosity',viscf); visc=viscf
      end block initialize_fluid
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         pmesh=partmesh(nvar=0,name='lpt')
         call lp%update_partmesh(pmesh)
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
         call ens_out%add_vector('source',lp%srcU,lp%srcV,lp%srcW)
         call ens_out%add_scalar('epsp',lp%VF)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call lp%get_max()
         ! Create simulation monitor
         mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(lp%np,'Particle number')
         call mfile%add_column(lp%VFmean,'Mean VF')
         call mfile%add_column(lp%Umin,'Particle Umin')
         call mfile%add_column(lp%Umax,'Particle Umax')
         call mfile%add_column(lp%Vmin,'Particle Vmin')
         call mfile%add_column(lp%Vmax,'Particle Vmax')
         call mfile%add_column(lp%Wmin,'Particle Wmin')
         call mfile%add_column(lp%Wmax,'Particle Wmax')
         call mfile%add_column(lp%dmin,'Particle dmin')
         call mfile%add_column(lp%dmax,'Particle dmax')
         call mfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call time%increment()
         
         ! Collide particles
         call lp%collide(dt=time%dt)
         
         ! Advance particles by dt
         call lp%advance(dt=time%dt,U=U,V=V,W=W,rho=rho,visc=visc)
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call lp%get_max()
         call mfile%write()
         
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
