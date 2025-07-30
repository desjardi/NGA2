!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single-phase incompressible flow solver and corresponding time tracker
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(sgsmodel),    public :: sgs

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,forcefile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU

   !> Fluid viscosity
   real(WP) :: visc

   !> Channel forcing
   real(WP) :: Ubulk
   real(WP) :: meanU,meanW

 contains


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
      end block allocate_work_arrays



      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker


      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         use mpi_f08,  only: MPI_SUM,MPI_MAX,MPI_ALLREDUCE
         use parallel, only: MPI_REAL_WP
         integer :: i,j,k,ierr
         real(WP) :: amp,vel,geo_fac,my_Ly,Ly
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
        ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initialize velocity based on specified bulk
         call param_read('Uh',Ubulk)
         ! Rescale bulk velocity to enforce correct mean velocity at x=0
         my_Ly=0.0_WP
         i=fs%cfg%imin; k=fs%cfg%kmin
         if (fs%cfg%imin_.le.i.and.fs%cfg%imax_.ge.i.and.&
              fs%cfg%kmin_.le.k.and.fs%cfg%kmax_.ge.k) then
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               my_Ly=my_Ly+fs%cfg%dy(j)*fs%cfg%VF(i,j,k)
            end do
         end if
         call MPI_ALLREDUCE(my_Ly,Ly,1,MPI_REAL_WP,MPI_SUM,fs%cfg%ycomm,ierr)
         call MPI_ALLREDUCE(Ly,my_Ly,1,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr); Ly=my_Ly
         geo_fac=Ly*fs%cfg%xL*fs%cfg%zL/fs%cfg%fluid_vol
         Ubulk=Ubulk*geo_fac
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         where (fs%umask.eq.0) fs%U=Ubulk
         meanU=Ubulk
         meanW=0.0_WP
         ! To facilitate transition
         call param_read('Perturbation',amp)
         vel=sqrt(Ubulk**2)
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=fs%U(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)
                  if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
               end do
            end do
         end do
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=fs%U(i,j,k)+random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)*fs%U(i,j,k)
                  if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)+random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)*fs%W(i,j,k)
               end do
            end do
         end do
         if (fs%cfg%nx.eq.1) fs%U=0.0_WP
         if (fs%cfg%ny.eq.1) fs%V=0.0_WP
         if (fs%cfg%nz.eq.1) fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver


      ! Create an LES model
      create_sgs: block
        sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='channel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('levelset',cfg%Gib)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('viscosity',fs%visc)
         call ens_out%add_scalar('visc_sgs',sgs%visc)
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
         ! Create forcing monitor
         forcefile=monitor(fs%cfg%amRoot,'forcing')
         call forcefile%add_column(time%n,'Timestep number')
         call forcefile%add_column(time%t,'Time')
         call forcefile%add_column(meanU,'Bulk U')
         call forcefile%add_column(meanW,'Bulk W')
         call forcefile%write()
      end block create_monitor


   end subroutine simulation_init


   !> Time integrate our problem
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

         ! Turbulence modeling
         sgs_modeling: block
           use sgsmodel_class, only: dynamic_smag,vreman
           resU=fs%rho
           call fs%get_gradu(gradU)
           call fs%interp_vel(Ui,Vi,Wi)
           call fs%get_strainrate(SR)
           call sgs%get_visc(type=dynamic_smag,dt=time%dtold,rho=resU,gradu=gradU,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
           !sgs%visc=sgs%visc*fs%cfg%vf
           where (cfg%Gib.lt.0.0_WP) sgs%visc=0.0_WP
           fs%visc=visc+sgs%visc
         end block sgs_modeling

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
               real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol,VFx,VFz
               myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        VFx=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                        if (VFx.le.0.5_WP) cycle
                        myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*VFx*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
                        myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*VFx
                        VFz=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                        myW   =myW   +fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)*VFz*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))
                        myWvol=myWvol+fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)*VFz
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               where (fs%umask.eq.0) resU=resU+fs%rho*(Ubulk-meanU)
               call MPI_ALLREDUCE(myWvol,Wvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myW   ,meanW,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanW=meanW/Wvol
               where (fs%wmask.eq.0) resW=resW-fs%rho*meanW
            end block forcing

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
              integer :: i,j,k
              real(WP) :: VFx,VFy,VFz
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       VFx=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                       VFy=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                       VFz=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                       fs%U(i,j,k)=VFx*fs%U(i,j,k)
                       fs%V(i,j,k)=VFy*fs%V(i,j,k)
                       fs%W(i,j,k)=VFz*fs%W(i,j,k)
                    end do
                 end do
              end do
              call fs%cfg%sync(fs%U)
              call fs%cfg%sync(fs%V)
              call fs%cfg%sync(fs%W)
            end block ibforcing

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
         call forcefile%write()

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
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradU)

   end subroutine simulation_final





end module simulation
