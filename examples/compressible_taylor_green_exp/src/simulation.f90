!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fastcomp_class,    only: fastcomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Fast compressible flow solver and corresponding time tracker
   type(fastcomp),    public :: fs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   !> Monitoring of conservation
   type(monitor) :: consfile
   real(WP) :: RHOKint,RHOSint,DilDiss,SolDiss
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: dRHOdt,dRHOUdt,dRHOVdt,dRHOWdt,dRHOEdt
   real(WP), dimension(:,:,:)  , allocatable :: Ui,Vi,Wi,C,Ma
   
   !> Equation of state and flow conditions
   real(WP) :: Gamma,Mach,Reynolds,Prandtl
   
contains
   
   
   !> Sutherland's law for viscosity as a function of temperature
   subroutine get_visc()
      implicit none
      integer :: i,j,k
      do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
         fs%visc(i,j,k)=Reynolds**(-1.0_WP)*(1.4042_WP*(fs%E(i,j,k)/fs%Cv)**1.5_WP)/(fs%E(i,j,k)/fs%Cv+0.4042_WP)
      end do; end do; end do
   end subroutine get_visc
   
   
   !> Post-process conservation properties
   subroutine analyze_conservation()
      implicit none
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: tmp,XY,YZ,ZX
      ! Allocate temporary storage
      allocate(tmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(XY (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(YZ (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(ZX (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      ! Compute kinetic energy
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               tmp(i,j,k)=0.25_WP*(sum(fs%RHOU(i:i+1,j,k)*fs%U(i:i+1,j,k))+sum(fs%RHOV(i,j:j+1,k)*fs%V(i,j:j+1,k))+sum(fs%RHOW(i,j,k:k+1)*fs%W(i,j,k:k+1)))
            end do
         end do
      end do
      call cfg%integrate(tmp,integral=RHOKint)
      ! Compute entropy
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               tmp(i,j,k)=fs%Cv*log(fs%P(i,j,k)/fs%RHO(i,j,k)**Gamma)
            end do
         end do
      end do
      call cfg%integrate(tmp,integral=RHOSint)
      ! Calculate dilatational dissipation
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               tmp(i,j,k)=(fs%beta(i,j,k)+4.0_WP/3.0_WP*fs%visc(i,j,k))*(fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k)))**2
            end do
         end do
      end do
      call cfg%integrate(tmp,integral=DilDiss)
      ! Calculate solenoidal dissipation
      do k=cfg%kmin_,cfg%kmax_+1
         do j=cfg%jmin_,cfg%jmax_+1
            do i=cfg%imin_,cfg%imax_+1
               ! (dvdx-dudy)^2 at xy edge
               XY(i,j,k)=0.25_WP*sum(fs%visc(i-1:i,j-1:j,k))*(fs%dxi*(fs%V(i,j,k)-fs%V(i-1,j,k))-fs%dyi*(fs%U(i,j,k)-fs%U(i,j-1,k)))**2
               ! (dwdy-dvdz)^2 at yz edge
               YZ(i,j,k)=0.25_WP*sum(fs%visc(i,j-1:j,k-1:k))*(fs%dyi*(fs%W(i,j,k)-fs%W(i,j-1,k))-fs%dzi*(fs%V(i,j,k)-fs%V(i,j,k-1)))**2
               ! (dudz-dwdx)^2 at zx edge
               ZX(i,j,k)=0.25_WP*sum(fs%visc(i-1:i,j,k-1:k))*(fs%dzi*(fs%U(i,j,k)-fs%U(i,j,k-1))-fs%dxi*(fs%W(i,j,k)-fs%W(i-1,j,k)))**2
            end do
         end do
      end do
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               tmp(i,j,k)=0.25_WP*sum(XY(i:i+1,j:j+1,k))+0.25_WP*sum(YZ(i,j:j+1,k:k+1))+0.25_WP*sum(ZX(i:i+1,j,k:k+1))
            end do
         end do
      end do
      call cfg%integrate(tmp,integral=SolDiss)
      ! Deallocate storage
      deallocate(tmp,XY,YZ,ZX)
   end subroutine analyze_conservation


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Prepare EoS and flow conditions
      initialize_eos: block
         call param_read('Gamma'   ,gamma   )
         call param_read('Mach'    ,Mach    )
         call param_read('Reynolds',Reynolds)
         call param_read('Prandtl' ,Prandtl )
      end block initialize_eos
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dRHOdt (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:4))
         allocate(dRHOUdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:4))
         allocate(dRHOVdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:4))
         allocate(dRHOWdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:4))
         allocate(dRHOEdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:4))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(C (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ma(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker
      
      ! Create a fast compressible flow solver
      create_velocity_solver: block
         call fs%initialize(cfg=cfg,name='Compressible NS')
         fs%Cv=1.0_WP/(Gamma*(Gamma-1.0_WP)*Mach**2) ! Specific heat corresponding to normalized temperature
      end block create_velocity_solver
      
      ! Prepare initial conditions
      initial_conditions: block
         use mathtools, only: Pi,twoPi
         !use random,    only: random_uniform
         integer :: i,j,k
         ! Initialize primary variables
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Standard TGV
                  fs%P  (i,j,k)=1.0_WP/(Gamma*Mach**2)+(cos(2.0_WP*fs%cfg%xm(i))+cos(2.0_WP*fs%cfg%ym(j)))*(cos(2.0_WP*fs%cfg%zm(k))+2.0_WP)/16.0_WP
                  fs%E  (i,j,k)=1.0_WP/((Gamma-1.0_WP)*Gamma*Mach**2)
                  fs%RHO(i,j,k)=fs%P(i,j,k)/(fs%E(i,j,k)*(Gamma-1.0_WP))
                  fs%U  (i,j,k)=+sin(fs%cfg%x (i))*cos(fs%cfg%ym(j))*cos(fs%cfg%zm(k))
                  fs%V  (i,j,k)=-cos(fs%cfg%xm(i))*sin(fs%cfg%y (j))*cos(fs%cfg%zm(k))
                  fs%W  (i,j,k)=0.0_WP
                  ! Uniform flow
                  !fs%U(i,j,k)=1.0_WP+random_uniform(lo=-1.0e-14_WP,hi=+1.0e-14_WP)
                  !fs%V(i,j,k)=0.0_WP+random_uniform(lo=-1.0e-14_WP,hi=+1.0e-14_WP)
                  !fs%W(i,j,k)=0.0_WP+random_uniform(lo=-1.0e-14_WP,hi=+1.0e-14_WP)
                  !fs%P  (i,j,k)=1.0_WP/(Gamma*Mach**2)
                  !fs%RHO(i,j,k)=1000.0_WP
                  !fs%E  (i,j,k)=fs%P(i,j,k)/(fs%RHO(i,j,k)*(Gamma-1.0_WP))
               end do
            end do
         end do
         ! Initialize conserved variables
         call fs%get_momentum(); fs%RHOE=fs%RHO*fs%E; fs%T=fs%E/fs%Cv
         ! Interpolate velcoity, compute speed of sound and local Mach number
         call fs%interp_vel(Ui,Vi,Wi)
         C=sqrt(Gamma*(Gamma-1.0_WP)*fs%E)
         Ma=sqrt((Ui**2+Vi**2+Wi**2))/C
      end block initial_conditions
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='TaylorGreen')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('density' ,fs%RHO)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('energy'  ,fs%E)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('Mach',Ma)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(dt=time%dt,C=C,cfl=time%cfl)
         call fs%get_info()
         call analyze_conservation()
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
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(fs%Emax,'Emax')
         call mfile%add_column(fs%Emin,'Emin')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(cfg%vol_total,'Volume')
         call consfile%add_column(fs%RHOint ,'Mass')
         call consfile%add_column(fs%RHOUint,'U Momentum')
         call consfile%add_column(fs%RHOVint,'V Momentum')
         call consfile%add_column(fs%RHOWint,'W Momentum')
         call consfile%add_column(fs%RHOEint,'Internal energy')
         call consfile%add_column(   RHOKint,'Kinetic energy')
         call consfile%add_column(   RHOSint,'Entropy')
         call consfile%add_column(DilDiss,'Dilatation dissipation')
         call consfile%add_column(SolDiss,'Solenoidal dissipation')
         call consfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,C=C,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember conserved variables
         fs%RHOold =fs%RHO
         fs%RHOUold=fs%RHOU
         fs%RHOVold=fs%RHOV
         fs%RHOWold=fs%RHOW
         fs%RHOEold=fs%RHOE

         ! Recalculate viscosity and heat diffusivity
         call get_visc()
         fs%diff=fs%Cv*Gamma*fs%visc/Prandtl
         
         ! First RK step ====================================================================================
         ! Get RHS (primitive variables have been computed before)
         call fs%get_RHS(dRHOdt(:,:,:,1),dRHOUdt(:,:,:,1),dRHOVdt(:,:,:,1),dRHOWdt(:,:,:,1),dRHOEdt(:,:,:,1))
         ! Increment conserved variables
         fs%RHO =fs%RHOold +0.5_WP*time%dt*dRHOdt (:,:,:,1)
         fs%RHOU=fs%RHOUold+0.5_WP*time%dt*dRHOUdt(:,:,:,1)
         fs%RHOV=fs%RHOVold+0.5_WP*time%dt*dRHOVdt(:,:,:,1)
         fs%RHOW=fs%RHOWold+0.5_WP*time%dt*dRHOWdt(:,:,:,1)
         fs%RHOE=fs%RHOEold+0.5_WP*time%dt*dRHOEdt(:,:,:,1)
         ! Recompute primitive variables
         call fs%get_velocity(); fs%E=fs%RHOE/fs%RHO; fs%P=fs%RHOE*(Gamma-1.0_WP); fs%T=fs%E/fs%Cv
         
         ! Second RK step ===================================================================================
         ! Get RHS
         call fs%get_RHS(dRHOdt(:,:,:,2),dRHOUdt(:,:,:,2),dRHOVdt(:,:,:,2),dRHOWdt(:,:,:,2),dRHOEdt(:,:,:,2))
         ! Increment conserved variables
         fs%RHO =fs%RHOold +0.5_WP*time%dt*dRHOdt (:,:,:,2)
         fs%RHOU=fs%RHOUold+0.5_WP*time%dt*dRHOUdt(:,:,:,2)
         fs%RHOV=fs%RHOVold+0.5_WP*time%dt*dRHOVdt(:,:,:,2)
         fs%RHOW=fs%RHOWold+0.5_WP*time%dt*dRHOWdt(:,:,:,2)
         fs%RHOE=fs%RHOEold+0.5_WP*time%dt*dRHOEdt(:,:,:,2)
         ! Recompute primitive variables
         call fs%get_velocity(); fs%E=fs%RHOE/fs%RHO; fs%P=fs%RHOE*(Gamma-1.0_WP); fs%T=fs%E/fs%Cv
         
         ! Third RK step ====================================================================================
         ! Get RHS
         call fs%get_RHS(dRHOdt(:,:,:,3),dRHOUdt(:,:,:,3),dRHOVdt(:,:,:,3),dRHOWdt(:,:,:,3),dRHOEdt(:,:,:,3))
         ! Increment conserved variables
         fs%RHO =fs%RHOold +1.0_WP*time%dt*dRHOdt (:,:,:,3)
         fs%RHOU=fs%RHOUold+1.0_WP*time%dt*dRHOUdt(:,:,:,3)
         fs%RHOV=fs%RHOVold+1.0_WP*time%dt*dRHOVdt(:,:,:,3)
         fs%RHOW=fs%RHOWold+1.0_WP*time%dt*dRHOWdt(:,:,:,3)
         fs%RHOE=fs%RHOEold+1.0_WP*time%dt*dRHOEdt(:,:,:,3)
         ! Recompute primitive variables
         call fs%get_velocity(); fs%E=fs%RHOE/fs%RHO; fs%P=fs%RHOE*(Gamma-1.0_WP); fs%T=fs%E/fs%Cv
         
         ! Fourth RK step ===================================================================================
         ! Get RHS
         call fs%get_RHS(dRHOdt(:,:,:,4),dRHOUdt(:,:,:,4),dRHOVdt(:,:,:,4),dRHOWdt(:,:,:,4),dRHOEdt(:,:,:,4))
         ! Increment conserved variables
         fs%RHO =fs%RHOold +time%dt/6.0_WP*(dRHOdt (:,:,:,1)+2.0_WP*dRHOdt (:,:,:,2)+2.0_WP*dRHOdt (:,:,:,3)+dRHOdt (:,:,:,4))
         fs%RHOU=fs%RHOUold+time%dt/6.0_WP*(dRHOUdt(:,:,:,1)+2.0_WP*dRHOUdt(:,:,:,2)+2.0_WP*dRHOUdt(:,:,:,3)+dRHOUdt(:,:,:,4))
         fs%RHOV=fs%RHOVold+time%dt/6.0_WP*(dRHOVdt(:,:,:,1)+2.0_WP*dRHOVdt(:,:,:,2)+2.0_WP*dRHOVdt(:,:,:,3)+dRHOVdt(:,:,:,4))
         fs%RHOW=fs%RHOWold+time%dt/6.0_WP*(dRHOWdt(:,:,:,1)+2.0_WP*dRHOWdt(:,:,:,2)+2.0_WP*dRHOWdt(:,:,:,3)+dRHOWdt(:,:,:,4))
         fs%RHOE=fs%RHOEold+time%dt/6.0_WP*(dRHOEdt(:,:,:,1)+2.0_WP*dRHOEdt(:,:,:,2)+2.0_WP*dRHOEdt(:,:,:,3)+dRHOEdt(:,:,:,4))
         ! Recompute primitive variables
         call fs%get_velocity(); fs%E=fs%RHOE/fs%RHO; fs%P=fs%RHOE*(Gamma-1.0_WP); fs%T=fs%E/fs%Cv
         
         ! Interpolate velcoity, compute speed of sound and local Mach number
         call fs%interp_vel(Ui,Vi,Wi)
         C=sqrt(Gamma*(Gamma-1.0_WP)*fs%E)
         Ma=sqrt((Ui**2+Vi**2+Wi**2))/C
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_info()
         call analyze_conservation()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(dRHOdt,dRHOUdt,dRHOVdt,dRHOWdt,dRHOEdt,Ui,Vi,Wi,C,Ma)
   end subroutine simulation_final
   
   
end module simulation