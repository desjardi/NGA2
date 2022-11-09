!> Definition for a block 1 simulation: zoomed in cough machine
module block1_class
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: t_wall,L_mouth,H_mouth,W_mouth,L_film,H_film,W_film,L_lip
   use config_class,      only: config
   use incomp_class,      only: incomp
   use hypre_uns_class,   only: hypre_uns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   use object_timer,      only: objtimer
   implicit none
   private

   public :: block1

   !> block 1 object
   type :: block1
      class(config), pointer :: cfg                                        !< Pointer to config
      type(incomp) :: fs                                                   !< Single-phase incompressible flow solver
      type(hypre_uns) :: ps                                                !< Unstructured HYPRE pressure solver
      type(hypre_uns) :: is                                                !< Unstructured HYPRE implicit solver
      type(timetracker) :: time                                            !< Time tracker
      type(objtimer) :: timer                                              !< Method timer
      type(sgsmodel) ::  sgs                                               !< SGS model
      type(ensight) :: ens_out                                             !< Ensight output
      type(event) :: ens_evt                                               !< Ensight output event
      type(monitor) :: mfile,cflfile,timerfile,timersummaryfile,forcefile  !< Monitor files
      type(datafile) :: df                                                 !< Datafile for restart
      !> Private work arrays
      real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:), allocatable :: SR
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block1

   !> Fluid viscosity
   real(WP) :: visc

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW
   
contains

   !> Initialization of block 1
   subroutine init(b,restart_test)
      use param, only: param_read
      implicit none
      class(block1), intent(inout) :: b
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
         b%time=timetracker(b%cfg%amRoot,name='duct')
         call param_read('1 Max timestep size',b%time%dtmax)
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
         b%timer=objtimer(b%cfg%amRoot,name='duct_timer')
      end block initialize_objtimer

      ! Create a single-phase flow solver with bconds
      create_and_initialize_flow_solver: block
         use hypre_uns_class, only: pcg_amg,gmres_amg
         use mathtools, only: twoPi
         integer :: i,j,k
         real(WP) :: amp,vel
         ! Create a single-phase flow solver
         b%fs=incomp(cfg=b%cfg,name='Single-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Gas dynamic viscosity',visc); b%fs%visc=visc
         ! Assign constant density to each phase
         call param_read('Gas density',b%fs%rho)
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
         ! Handle restart
         if (restart_test) then
            call b%df%pullvar(name='U'  ,var=b%fs%U  )
            call b%df%pullvar(name='V'  ,var=b%fs%V  )
            call b%df%pullvar(name='W'  ,var=b%fs%W  )
            call b%df%pullvar(name='P'  ,var=b%fs%P  )
         end if
         ! Initialize velocity based on specified bulk
         call param_read('Ubulk',Ubulk)
         call param_read('Wbulk',Wbulk)
         where (b%fs%umask.eq.0) b%fs%U=Ubulk
         where (b%fs%wmask.eq.0) b%fs%W=Wbulk
         meanU=Ubulk
         meanW=Wbulk
         ! To facilitate transition
         call param_read('Perturbation',amp)
         vel=sqrt(Ubulk**2+Wbulk**2)
         do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
            do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
               do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
                  if (b%fs%umask(i,j,k).eq.0) b%fs%U(i,j,k)=b%fs%U(i,j,k)+amp*vel*cos(8.0_WP*twoPi*b%fs%cfg%zm(k)/b%fs%cfg%zL)
                  if (b%fs%wmask(i,j,k).eq.0) b%fs%W(i,j,k)=b%fs%W(i,j,k)+amp*vel*cos(8.0_WP*twoPi*b%fs%cfg%xm(i)/b%fs%cfg%xL)
               end do
            end do
         end do
         ! Calculate cell-centered velocities and divergence
         call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
         call b%fs%get_div()
      end block create_and_initialize_flow_solver

      ! Create an LES model
      create_sgs: block
         b%sgs=sgsmodel(cfg=b%fs%cfg,umask=b%fs%umask,vmask=b%fs%vmask,wmask=b%fs%wmask)
         ! Handle restart
         if (restart_test) then
            call b%df%pullvar(name='LM',var=b%sgs%LM)
            call b%df%pullvar(name='MM',var=b%sgs%MM)
         end if
      end block create_sgs

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(b%cfg,'duct')
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
         ! Create object time tracker monitor
         b%timerfile=monitor(b%fs%cfg%amRoot,'duct_timers')
         call b%timerfile%add_column(b%time%n,'Timestep number')
         call b%timerfile%add_column(b%time%t,'Simulation Time')
         call b%timerfile%add_column(b%timer%sgs_wt,'sgs_visc Wall Time')
         call b%timerfile%add_column(b%timer%implicit_wt,'imp_solv Wall Time')
         call b%timerfile%add_column(b%timer%pressure_wt,'pres_solv Wall Time')
         call b%timerfile%add_column(b%timer%step_wt,'time_step Wall Time')
         call b%timerfile%write()
         ! Create object time and cost summary monitor
         b%timersummaryfile=monitor(b%fs%cfg%amRoot,'duct_timer_summary')
         call b%timersummaryfile%add_column(b%time%n,'Timestep number')
         call b%timersummaryfile%add_column(b%time%t,'Simulation Time')
         call b%timersummaryfile%add_column(b%timer%sgs_wt_total,'sgs_visc Total Hours')
         call b%timersummaryfile%add_column(b%timer%sgs_core_hours,'sgs_visc Core Hours')
         call b%timersummaryfile%add_column(b%timer%implicit_wt_total,'imp_solv WT Hours')
         call b%timersummaryfile%add_column(b%timer%implicit_core_hours,'imp_solv Core Hours')
         call b%timersummaryfile%add_column(b%timer%pressure_wt_total,'pres_solv WT Hours')
         call b%timersummaryfile%add_column(b%timer%pressure_core_hours,'pres_solv Core Hours')
         call b%timersummaryfile%add_column(b%timer%step_wt_total,'time_step WT Hours')
         call b%timersummaryfile%add_column(b%timer%step_core_hours,'time_step Core Hours')
         call b%timersummaryfile%write()
         ! Create forcing monitor
         b%forcefile=monitor(b%fs%cfg%amRoot,'forcing')
         call b%forcefile%add_column(b%time%n,'Timestep number')
         call b%forcefile%add_column(b%time%t,'Time')
         call b%forcefile%add_column(meanU,'Bulk U')
         call b%forcefile%add_column(meanW,'Bulk W')
         call b%forcefile%write()
      end block create_monitor


   end subroutine init


   !> Take a time step with block 1
   subroutine step(b)
      use mpi,        only: mpi_wtime
      implicit none
      class(block1), intent(inout) :: b
      real(WP) :: starttime,endtime

      ! Start time step timer
      starttime=mpi_wtime()

      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W

      ! Reset here gas viscosity
      b%fs%visc=visc

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

         ! Add body forcing
         forcing: block
            use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
            use parallel, only: MPI_REAL_WP
            integer :: i,j,k,ierr
            real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
            myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
            do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
               do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
                  do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                     if (b%fs%umask(i,j,k).eq.0) then
                        myU   =myU   +b%fs%cfg%dxm(i)*b%fs%cfg%dy(j)*b%fs%cfg%dz(k)*(2.0_WP*b%fs%U(i,j,k)-b%fs%Uold(i,j,k))
                        myUvol=myUvol+b%fs%cfg%dxm(i)*b%fs%cfg%dy(j)*b%fs%cfg%dz(k)
                     end if
                     if (b%fs%wmask(i,j,k).eq.0) then
                        myW   =myW   +b%fs%cfg%dx(i)*b%fs%cfg%dy(j)*b%fs%cfg%dzm(k)*(2.0_WP*b%fs%W(i,j,k)-b%fs%Wold(i,j,k))
                        myWvol=myWvol+b%fs%cfg%dx(i)*b%fs%cfg%dy(j)*b%fs%cfg%dzm(k)
                     end if
                  end do
               end do
            end do
            call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,b%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,b%fs%cfg%comm,ierr); meanU=meanU/Uvol
            where (b%fs%umask.eq.0) b%resU=b%resU+Ubulk-meanU
            call MPI_ALLREDUCE(myWvol,Wvol ,1,MPI_REAL_WP,MPI_SUM,b%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(myW   ,meanW,1,MPI_REAL_WP,MPI_SUM,b%fs%cfg%comm,ierr); meanW=meanW/Wvol
            where (b%fs%wmask.eq.0) b%resW=b%resW+Wbulk-meanW
         end block forcing

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

      ! End time steo timer
      endtime=mpi_wtime()

      ! Wall time spent in current time step
      b%cfg%step_wt=endtime-starttime

      ! Output to ensight
      if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)

      ! Update object time trackers
      call b%timer%sgs_visc_timer(b%cfg,b%sgs)
      call b%timer%implicit_timer(b%cfg,b%fs)
      call b%timer%pressure_timer(b%cfg,b%fs)
      call b%timer%step_timer(b%cfg)

      ! Perform and output monitoring
      call b%fs%get_max()
      call b%mfile%write()
      call b%cflfile%write()
      call b%timerfile%write()
      call b%timersummaryfile%write()

   end subroutine step


   !> Finalize b1 simulation
   subroutine final(b)
      implicit none
      class(block1), intent(inout) :: b

      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR)

   end subroutine final


end module block1_class
