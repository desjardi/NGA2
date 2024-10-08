!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,D
   use fft3d_class,       only: fft3d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp) :: fs
   type(fft3d)  :: ps
   type(timetracker) :: time
   
   !> Implicit solver
   logical :: use_implicit
   type(ddadi) :: vs
   
   !> SGS model
   logical :: use_sgs
   type(sgsmodel) :: sgs
   
   !> Ensight postprocessing
   type(event) :: ens_evt
   type(ensight) :: ens_out
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP) :: visc,Ubulk,meanU,bforce
   
   !> Provide a pardata objects for restarts
   logical :: restarted
   type(event) :: save_evt
   type(pardata) :: df

   
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
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Configure pressure solver
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Check if implicit velocity solver is used
         call param_read('Use implicit solver',use_implicit)
         if (use_implicit) then
            ! Configure implicit solver
            vs=ddadi(cfg=cfg,name='Velocity',nst=7)
            ! Finish flow solver setup
            call fs%setup(pressure_solver=ps,implicit_solver=vs)
         else
            ! Finish flow solver setup
            call fs%setup(pressure_solver=ps)
         end if
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: amp,VFx,VFy,VFz
         ! Initial fields
         call param_read('Bulk velocity',Ubulk)
         fs%U=Ubulk; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         meanU=Ubulk
         bforce=0.0_WP
         ! For faster transition
         call param_read('Fluctuation amp',amp,default=0.0_WP)
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  ! Add fluctuations for faster transition
                  fs%U(i,j,k)=fs%U(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)*cos(8.0_WP*twoPi*fs%cfg%ym(j)/fs%cfg%yL)
                  fs%V(i,j,k)=fs%V(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  fs%W(i,j,k)=fs%W(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  ! Remove values in the wall
                  VFx=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                  VFy=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                  VFz=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                  fs%U(i,j,k)=fs%U(i,j,k)*VFx
                  fs%V(i,j,k)=fs%V(i,j,k)*VFy
                  fs%W(i,j,k)=fs%W(i,j,k)*VFz
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
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',use_sgs)
         if (use_sgs) sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs


      ! Handle restart here
      perform_restart: block
         use string,  only: str_medium
         use filesys, only: makedir,isdir
         character(len=str_medium) :: filename
         ! Create event for saving restart files
         save_evt=event(time,'Restart output')
         call param_read('Restart output period',save_evt%tper)
         ! Check if a restart file was provided
         call param_read('Restart from',filename,default='')
         restarted=.false.; if (len_trim(filename).gt.0) restarted=.true.
         ! Restart the simulation
         if (restarted) then
            ! Read in the file
            call df%initialize(pg=cfg,iopartition=[cfg%npx,cfg%npy,cfg%npz],fdata=trim(filename))
            ! Put the data at the right place
            call df%pull(name='U',var=fs%U)
            call df%pull(name='V',var=fs%V)
            call df%pull(name='W',var=fs%W)
            call df%pull(name='P',var=fs%P)
            ! Update cell-centered velocity
            call fs%interp_vel(Ui,Vi,Wi)
            ! Update divergence
            call fs%get_div()
            ! Also update time
            call df%pull(name='t' ,val=time%t )
            call df%pull(name='dt',val=time%dt)
            time%told=time%t-time%dt
         else
            ! Prepare a new directory for storing files for restart
            if (cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! If we are not restarting, we will still need a datafile for saving restart files
            call df%initialize(pg=cfg,iopartition=[cfg%npx,cfg%npy,cfg%npz],filename=trim(cfg%name),nval=2,nvar=4)
            df%valname=['dt','t ']; df%varname=['U','V','W','P']
         end if
      end block perform_restart
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pipe')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('Gib',cfg%Gib)
         call ens_out%add_scalar('pressure',fs%P)
         if (use_sgs) call ens_out%add_scalar('visc_sgs',sgs%visc)
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
         call mfile%add_column(meanU,'Bulk U')
         call mfile%add_column(bforce,'Body force')
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
         
         ! Turbulence modeling
         if (use_sgs) then
            sgs_modeling: block
               use sgsmodel_class, only: vreman
               resU=fs%rho
               call fs%get_gradu(gradU)
               call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
               where (cfg%Gib.gt.0.0_WP) sgs%visc=0.0_WP
               fs%visc=visc+sgs%visc
            end block sgs_modeling
         end if
         
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
               use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
               use parallel, only: MPI_REAL_WP
               integer :: i,j,k,ierr
               real(WP) :: Uvol,VFx
               Uvol=0.0_WP; meanU=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        VFx=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                        if (VFx.le.0.5_WP) cycle
                        meanU=meanU+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*VFx*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
                        Uvol =Uvol +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*VFx
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               resU=resU+fs%rho*(Ubulk-meanU)
               bforce=fs%rho*(Ubulk-meanU)/time%dt
            end block forcing
            
            ! Finish update
            if (use_implicit) then
               ! Form implicit residuals
               call fs%solve_implicit(time%dt,resU,resV,resW)
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU
               fs%V=2.0_WP*fs%V-fs%Vold+resV
               fs%W=2.0_WP*fs%W-fs%Wold+resW
            else
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
               fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
               fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            end if
            
            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%U(i,j,k)=fs%U(i,j,k)*sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                        fs%V(i,j,k)=fs%V(i,j,k)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                        fs%W(i,j,k)=fs%W(i,j,k)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
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

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
               use string, only: str_medium
               character(len=str_medium) :: timestamp
               ! Prefix for files
               write(timestamp,'(es12.5)') time%t
               ! Populate df and write it
               call df%push(name='t' ,val=time%t )
               call df%push(name='dt',val=time%dt)
               call df%push(name='U' ,var=fs%U   )
               call df%push(name='V' ,var=fs%V   )
               call df%push(name='W' ,var=fs%W   )
               call df%push(name='P' ,var=fs%P   )
               call df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            end block save_restart
         end if
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)
      
   end subroutine simulation_final
   
end module simulation
