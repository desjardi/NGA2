!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Only get a LPT solver and corresponding time tracker
   type(lpt),         public :: lp
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Velocity arrays
   real(WP), dimension(:,:,:), allocatable :: U,V,W
   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Initialize our LPT
      initialize_lpt: block
         lp=lpt(cfg=cfg,name='LPT')
         call param_read('Particle density',lp%rho)
      end block initialize_lpt
      
      
      ! Prepare our velocity field based on a Taylor vortex
      initialize_velocity: block
         use mathtools, only: twoPi
         integer :: i,j,k
         ! Allocate arrays
         allocate(U(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         allocate(V(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         allocate(W(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         ! Initialize to solid body rotation
         do k=lp%cfg%kmino_,lp%cfg%kmaxo_
            do j=lp%cfg%jmino_,lp%cfg%jmaxo_
               do i=lp%cfg%imino_,lp%cfg%imaxo_
                  U(i,j,k)=-twoPi*lp%cfg%ym(j)
                  V(i,j,k)=+twoPi*lp%cfg%xm(i)
                  W(i,j,k)=0.0_WP
               end do
            end do
         end do
      end block initialize_velocity
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=lp%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='particles')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',lp%pmesh)
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
         
         ! Advance particles by dt
         call lp%advance(dt=time%dt,U=U,V=V,W=W)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
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
      deallocate(U,V,W)
      
   end subroutine simulation_final
   
end module simulation
