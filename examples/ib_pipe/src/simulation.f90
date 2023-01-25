!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use hypre_str_class,   only: hypre_str
  use incomp_class,      only: incomp
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private

  !> Get an an incompressible solver, pressure solver, and corresponding time tracker
  type(incomp),      public :: fs
  type(hypre_str),   public :: ps
  type(hypre_str),   public :: vs
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile

  public :: simulation_init,simulation_run,simulation_final

  !> Work arrays
  real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
  real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
  real(WP), dimension(:,:,:), allocatable :: G,VF
  real(WP) :: Ubulk,Umean,Dpipe

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
      time%itmax=2
    end block initialize_timetracker


    ! Create an incompressible flow solver without bconds
    create_flow_solver: block
      use hypre_str_class, only: pcg_pfmg,gmres_pfmg
      real(WP) :: visc
      ! Create flow solver
      fs=incomp(cfg=cfg,name='Incompressible NS')
      ! Set the flow properties
      call param_read('Density',fs%rho)
      call param_read('Dynamic viscosity',visc); fs%visc=visc
      ! Configure pressure solver
      ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
      ps%maxlevel=10
      call param_read('Pressure iteration',ps%maxit)
      call param_read('Pressure tolerance',ps%rcvg)
      ! Configure implicit velocity solver
      vs=hypre_str(cfg=cfg,name='Velocity',method=gmres_pfmg,nst=7)
      call param_read('Implicit iteration',vs%maxit)
      call param_read('Implicit tolerance',vs%rcvg)
      ! Setup the solver
      call fs%setup(pressure_solver=ps,implicit_solver=vs)
    end block create_flow_solver


    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(G       (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(VF      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
    end block allocate_work_arrays


    ! Initialize our velocity field
    initialize_velocity: block
      use mathtools, only: twoPi
      use random, only: random_uniform
      integer :: i,j,k
      real(WP) :: amp
      ! Initial fields
      call param_read('Pipe diameter',Dpipe)
      call param_read('Bulk velocity',Ubulk); Umean=Ubulk
      call param_read('Fluctuation amp',amp,default=0.0_WP)
      fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
      ! For faster transition
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               fs%U(i,j,k)=Ubulk*(1.0_WP+random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp))
               fs%V(i,j,k)=Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)
               fs%U(i,j,k)=fs%U(i,j,k)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)*cos(8.0_WP*twoPi*fs%cfg%ym(j)/fs%cfg%yL)
               fs%V(i,j,k)=fs%V(i,j,k)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
            end do
         end do
      end do
      call fs%cfg%sync(fs%U)
      call fs%cfg%sync(fs%V)
      call fs%cfg%sync(fs%W)
      ! Compute cell-centered velocity
      call fs%interp_vel(Ui,Vi,Wi)
      ! Compute divergence
      call fs%get_div()
    end block initialize_velocity


    ! Initialize IBM fields
    initialize_ibm: block
      integer :: i,j,k
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               VF(i,j,k)=get_VF(i,j,k,'SC')
               G(i,j,k)=0.5_WP*Dpipe-sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
            end do
         end do
      end do
      call fs%cfg%sync(VF)
      call fs%cfg%sync(G)
    end block initialize_ibm

    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg,name='pipe')
      ! Create event for Ensight output
      ens_evt=event(time=time,name='Ensight output')
      call param_read('Ensight output period',ens_evt%tper)
      ! Add variables to output
      call ens_out%add_vector('velocity',Ui,Vi,Wi)
      call ens_out%add_scalar('levelset',G)
      call ens_out%add_scalar('ibm_vf',VF)
      call ens_out%add_scalar('pressure',fs%P)
      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)
    end block create_ensight

    ! Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_max()
      ! Create simulation monitor
      mfile=monitor(fs%cfg%amRoot,'simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(time%cfl,'Maximum CFL')
      call mfile%add_column(Umean,'Umean')
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
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call fs%get_cfl(time%dt,time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Remember old velocity
       fs%Uold=fs%U
       fs%Vold=fs%V
       fs%Wold=fs%W

       ! Perform sub-iterations
       do while (time%it.le.time%itmax)

          ! Build mid-time velocity
          fs%U=0.5_WP*(fs%U+fs%Uold)
          fs%V=0.5_WP*(fs%V+fs%Vold)
          fs%W=0.5_WP*(fs%W+fs%Wold)

          ! Explicit calculation of drho*u/dt from NS
          call fs%get_dmomdt(resU,resV,resW)

          ! Assemble explicit residual
          resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
          resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
          resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW

          ! Add body forcing
          forcing: block
            use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
            use parallel, only: MPI_REAL_WP
            integer :: i,j,k,ierr
            real(WP) :: myU,myUvol,Uvol,VFx
            myU=0.0_WP; myUvol=0.0_WP
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     VFx=get_VF(i,j,k,'U')
                     myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*VFx*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
                     myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*VFx
                  end do
               end do
            end do
            call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(myU   ,Umean,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Umean=Umean/Uvol
            resU=resU+Ubulk-Umean
          end block forcing

          ! Form implicit residuals
          call fs%solve_implicit(time%dt,resU,resV,resW)

          ! Apply these residuals
          fs%U=2.0_WP*fs%U-fs%Uold+resU
          fs%V=2.0_WP*fs%V-fs%Vold+resV
          fs%W=2.0_WP*fs%W-fs%Wold+resW

          ! Apply direct forcing to enforce BC at the pipe walls
          ibm_correction: block
            integer :: i,j,k
            real(WP) :: VFx,VFy,VFz
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     VF(i,j,k)=get_VF(i,j,k,'SC')
                     VFx      =get_VF(i,j,k,'U')
                     VFy      =get_VF(i,j,k,'V')
                     VFz      =get_VF(i,j,k,'W')
                     fs%U(i,j,k)=VFx*fs%U(i,j,k)
                     fs%V(i,j,k)=VFy*fs%V(i,j,k)
                     fs%W(i,j,k)=VFz*fs%W(i,j,k)
                  end do
               end do
            end do
            call fs%cfg%sync(fs%U)
            call fs%cfg%sync(fs%V)
            call fs%cfg%sync(fs%W)
            call fs%cfg%sync(VF)
          end block ibm_correction

          ! Apply other boundary conditions on the resulting fields
          call fs%apply_bcond(time%t,time%dt)

          ! Solve Poisson equation
          call fs%correct_mfr()
          call fs%get_div()
          fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
          fs%psolv%sol=0.0_WP
          call fs%psolv%solve()
          call fs%shift_p(fs%psolv%sol)

          ! Correct velocity
          call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
          fs%P=fs%P+fs%psolv%sol
          fs%U=fs%U-time%dt*resU/fs%rho
          fs%V=fs%V-time%dt*resV/fs%rho
          fs%W=fs%W-time%dt*resW/fs%rho
          
          ! Increment sub-iteration counter
          time%it=time%it+1

          ! Increment sub-iteration counter
          time%it=time%it+1

       end do

       ! Recompute interpolated velocity and divergence
       call fs%interp_vel(Ui,Vi,Wi)
       call fs%get_div()

       ! Output to ensight
       if (ens_evt%occurs()) call ens_out%write_data(time%t)

       ! Perform and output monitoring
       call fs%get_max()
       call mfile%write()
       call cflfile%write()

    end do

  end subroutine simulation_run


  !> Get volume fraction for direct forcing
  function get_VF(i,j,k,dir) result(VF)
    implicit none
    integer, intent(in)    :: i,j,k
    character(len=*)       :: dir
    real(WP)               :: VF
    real(WP)               :: r,eta,lam,delta,VFx,VFy,VFz
    real(WP), dimension(3) :: norm
    select case(trim(dir))
    case('U')
       delta=(fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k))**(1.0_WP/3.0_WP)
       r=sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)+epsilon(1.0_WP)
       norm(1)=0.0_WP; norm(2)=fs%cfg%ym(j)/r; norm(3)=fs%cfg%zm(k)/r
    case('V')
       delta=(fs%cfg%dx(i)*fs%cfg%dym(j)*fs%cfg%dz(k))**(1.0_WP/3.0_WP)
       r=sqrt(fs%cfg%y(j)**2+fs%cfg%zm(k)**2)+epsilon(1.0_WP)
       norm(1)=0.0_WP; norm(2)=fs%cfg%y(j)/r; norm(3)=fs%cfg%zm(k)/r
    case('W')
       delta=(fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k))**(1.0_WP/3.0_WP)
       r=sqrt(fs%cfg%ym(j)**2+fs%cfg%z(k)**2)+epsilon(1.0_WP)
       norm(1)=0.0_WP; norm(2)=fs%cfg%ym(j)/r; norm(3)=fs%cfg%z(k)/r
    case default
       delta=(fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k))**(1.0_WP/3.0_WP)
       r=sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
       norm(1)=0.0_WP; norm(2)=fs%cfg%ym(j)/r; norm(3)=fs%cfg%zm(k)/r
    end select
    lam=sum(abs(norm)); eta=0.065_WP*(1.0_WP-lam**2)+0.39_WP
    VF=0.5_WP*(1.0_WP-tanh((r-0.5_WP*Dpipe)/(sqrt(2.0_WP)*lam*eta*delta+epsilon(1.0_WP))))
  end function get_VF


  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Get rid of all objects - need destructors
    ! monitor
    ! ensight
    ! timetracker

    ! Deallocate work arrays
    deallocate(resU,resV,resW,Ui,Vi,Wi)

  end subroutine simulation_final

end module simulation
