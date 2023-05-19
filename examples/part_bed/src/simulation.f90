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
   
   !> Fluid info
   real(WP), dimension(:,:,:), allocatable :: U,V,W
   real(WP), dimension(:,:,:), allocatable :: rho,visc
   
   !> Particle barycenter and reference diameter
   real(WP) :: mean_y
   real(WP) :: dp
   
contains
   
   
   !> Calculation of average particle position
   subroutine calc_barycenter
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: n,ierr
      integer :: my_number_particles,number_particles
      real(WP) :: my_mean_y
      ! Loop over local particles
      my_number_particles=0
      my_mean_y=0.0_WP
      do n=1,lp%np_
         my_number_particles=my_number_particles+1
         my_mean_y=my_mean_y+lp%p(n)%pos(2)
      end do
      call MPI_ALLREDUCE(my_number_particles,number_particles,1,MPI_INTEGER,MPI_SUM,lp%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_mean_y,mean_y,1,MPI_REAL_WP,MPI_SUM,lp%cfg%comm,ierr)
      mean_y=mean_y/real(max(number_particles,1),WP)
   end subroutine calc_barycenter
   
   
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
         time%itmax=1
      end block initialize_timetracker
      

      ! Initialize our LPT solver
      initialize_lpt: block
         use random, only: random_uniform
         integer :: i,np
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         
         ! Read in file
         !call lp%read(filename='bed.file')

         ! Root process initializes np particles randomly
         call param_read('Number of particles',np)
         if (lp%cfg%amRoot) then
            call lp%resize(np)
            do i=1,np
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! Assign random position in the domain
               lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
               &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
               &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
               if (cfg%nz.eq.1) lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin_)
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               ! Give zero collision force
               lp%p(i)%Acol=0.0_WP
               lp%p(i)%Tcol=0.0_WP
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
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n,default=0.7_WP)
         call param_read('Wall restitution',lp%e_w,default=lp%e_n)
         call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         
      end block initialize_lpt
      

      ! Prepare our fluid phase info based on a Taylor vortex
      initialize_fluid: block
         real(WP) :: rhof,viscf
         ! Allocate arrays
         allocate(rho (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); rho =1.0_WP
         allocate(visc(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); visc=0.0_WP
         allocate(U   (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); U   =0.0_WP
         allocate(V   (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); V   =0.0_WP
         allocate(W   (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); W   =0.0_WP
         ! Set constant density and viscosity
         call param_read('Density'  ,rhof ); rho =rhof
         call param_read('Viscosity',viscf); visc=viscf
      end block initialize_fluid
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=2,name='lpt')
         pmesh%varname(1)='radius'
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='Fcol'
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(  1,i)=0.5_WP*lp%p(i)%d
            pmesh%vec(:,1,i)=lp%p(i)%vel
            pmesh%vec(:,2,i)=lp%p(i)%Acol
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
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call lp%get_max()
         call calc_barycenter()
         ! Create simulation monitor
         mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(lp%np,'Particle number')
         call mfile%add_column(mean_y,'Mean y')
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
      use mathtools, only: twoPi
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
            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(  1,i)=0.5_WP*lp%p(i)%d
                  pmesh%vec(:,1,i)=lp%p(i)%vel
                  pmesh%vec(:,2,i)=lp%p(i)%Acol
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call lp%get_max()
         call calc_barycenter()
         call mfile%write()
         
      end do
      
      ! Output particle bed
      call lp%write(filename='part.file')
      
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
