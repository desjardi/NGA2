!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mpcomp_class,      only: mpcomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Multiphase compressible flow solver and corresponding time tracker
   type(mpcomp),      public :: fs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   !> Monitoring of conservation
   type(monitor) :: consfile
   real(WP) :: RHOKint,RHOSLint,RHOSGint,DilDiss,SolDiss
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma
   
   !> Equation of state and flow conditions
   real(WP) :: Mach
   real(WP) :: GammaL,PinfL,CvL,RHOL,GammaG,PinfG,CvG,RHOG
   
   !> Drop radius and center
   real(WP) :: radius=1.0_WP
   real(WP), dimension(3) :: center=[0.0_WP,0.0_WP,0.0_WP]
   
contains
   
   
   !> Sutherland's law for viscosity as a function of temperature
   !subroutine get_visc()
   !   implicit none
   !   integer :: i,j,k
   !   do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
   !      fs%visc(i,j,k)=Reynolds**(-1.0_WP)*(1.4042_WP*fs%T(i,j,k)**1.5_WP)/(fs%T(i,j,k)+0.4042_WP)
   !   end do; end do; end do
   !end subroutine get_visc
   
   
   !> P=EOS(RHO,E) for liquid
   pure real(WP) function get_PL(RHO,E)
      implicit none
      real(WP), intent(in) :: RHO,E
      get_PL=RHO*E*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> T=f(RHO,P) for liquid
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> C=f(RHO,P) for liquid
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(GammaL*(P+PinfL)/RHO)
   end function get_CL
   
   
   !> P=EOS(RHO,E) for gas
   pure real(WP) function get_PG(RHO,E)
      implicit none
      real(WP), intent(in) :: RHO,E
      get_PG=RHO*E*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> T=f(RHO,P) for gas
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> C=f(RHO,P) for gas
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(GammaG*(P+PinfG)/RHO)
   end function get_CG
   
   
   !> Post-process conservation properties
   subroutine analyze_conservation()
      implicit none
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: tmp,XY,YZ,ZX
      ! Allocate temporary storage
      allocate(tmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      !allocate(XY (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      !allocate(YZ (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      !allocate(ZX (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      ! Compute kinetic energy
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         tmp(i,j,k)=0.25_WP*(sum(fs%Q(i:i+1,j,k,5)*fs%U(i:i+1,j,k))+sum(fs%Q(i,j:j+1,k,6)*fs%V(i,j:j+1,k))+sum(fs%Q(i,j,k:k+1,7)*fs%W(i,j,k:k+1)))
      end do; end do; end do
      call cfg%integrate(tmp,integral=RHOKint)
      ! Compute liquid entropy
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         if (fs%VF(i,j,k).gt.0.0_WP) then
            tmp(i,j,k)=fs%VF(i,j,k)*fs%RHOL(i,j,k)*CvL*log((fs%PL(i,j,k)+PinfL)/fs%RHOL(i,j,k)**GammaL)
         else
            tmp(i,j,k)=0.0_WP
         end if
      end do; end do; end do
      call cfg%integrate(tmp,integral=RHOSLint)
      ! Compute gas entropy
      do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
         if (fs%VF(i,j,k).lt.1.0_WP) then
            tmp(i,j,k)=(1.0_WP-fs%VF(i,j,k))*fs%RHOG(i,j,k)*CvG*log((fs%PG(i,j,k)+PinfG)/fs%RHOG(i,j,k)**GammaG)
         else
            tmp(i,j,k)=0.0_WP
         end if
      end do; end do; end do
      call cfg%integrate(tmp,integral=RHOSGint)
      ! Calculate dilatational dissipation
      !do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
      !   tmp(i,j,k)=(fs%beta(i,j,k)+4.0_WP/3.0_WP*fs%visc(i,j,k))*(fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k)))**2
      !end do; end do; end do
      !call cfg%integrate(tmp,integral=DilDiss)
      ! Calculate solenoidal dissipation
      !do k=cfg%kmin_,cfg%kmax_+1; do j=cfg%jmin_,cfg%jmax_+1; do i=cfg%imin_,cfg%imax_+1
         ! (dvdx-dudy)^2 at xy edge
      !   XY(i,j,k)=0.25_WP*sum(fs%visc(i-1:i,j-1:j,k))*(fs%dxi*(fs%V(i,j,k)-fs%V(i-1,j,k))-fs%dyi*(fs%U(i,j,k)-fs%U(i,j-1,k)))**2
         ! (dwdy-dvdz)^2 at yz edge
      !   YZ(i,j,k)=0.25_WP*sum(fs%visc(i,j-1:j,k-1:k))*(fs%dyi*(fs%W(i,j,k)-fs%W(i,j-1,k))-fs%dzi*(fs%V(i,j,k)-fs%V(i,j,k-1)))**2
      !   ! (dudz-dwdx)^2 at zx edge
      !   ZX(i,j,k)=0.25_WP*sum(fs%visc(i-1:i,j,k-1:k))*(fs%dzi*(fs%U(i,j,k)-fs%U(i,j,k-1))-fs%dxi*(fs%W(i,j,k)-fs%W(i-1,j,k)))**2
      !end do; end do; end do
      !do k=cfg%kmin_,cfg%kmax_; do j=cfg%jmin_,cfg%jmax_; do i=cfg%imin_,cfg%imax_
      !   tmp(i,j,k)=0.25_WP*sum(XY(i:i+1,j:j+1,k))+0.25_WP*sum(YZ(i,j:j+1,k:k+1))+0.25_WP*sum(ZX(i:i+1,j,k:k+1))
      !end do; end do; end do
      !call cfg%integrate(tmp,integral=SolDiss)
      ! Deallocate storage
      deallocate(tmp)!,XY,YZ,ZX)
   end subroutine analyze_conservation
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Prepare EoS and flow conditions
      initialize_eos: block
         call param_read('GammaL'  ,GammaL  )
         call param_read('PinfL'   ,PinfL   )
         call param_read('RHOL'    ,RHOL    )
         call param_read('GammaG'  ,GammaG  )
         call param_read('PinfG'   ,PinfG   )
         call param_read('RHOG'    ,RHOG    )
         call param_read('Mach'    ,Mach    )
         CvL=1.0_WP/(GammaL*(GammaL-1.0_WP)*Mach**2)
         CvG=1.0_WP/(GammaG*(GammaG-1.0_WP)*Mach**2)
      end block initialize_eos
      
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
         call fs%initialize(cfg=cfg,getPL=get_PL,getTL=get_TL,getCL=get_CL,&
         &                          getPG=get_PG,getTG=get_TG,getCG=get_CG,name='Compressible NS')
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dQdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ma(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Prepare initial conditions
      initial_conditions: block
         use mms_geom,     only: initialize_volume_moments
         use mpcomp_class, only: VFlo
         use random,       only: random_uniform
         integer :: i,j,k
         ! Initialize primary variables
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Volume moments for a droplet
                  call initialize_volume_moments(lo=[cfg%x(i  ),cfg%y(j  ),cfg%z(k  )],hi=[cfg%x(i+1),cfg%y(j+1),cfg%z(k+1)],&
                  levelset=levelset_slab,time=0.0_WP,level=4,VFlo=VFlo,VF=fs%VF(i,j,k),BL=fs%BL(:,i,j,k),BG=fs%BG(:,i,j,k))
                  ! Mixture velocity
                  fs%U(i,j,k)=+Mach*sin(cfg%x (i))*cos(cfg%ym(j))*cos(cfg%zm(k))
                  fs%V(i,j,k)=-Mach*cos(cfg%xm(i))*sin(cfg%y (j))*cos(cfg%zm(k))
                  fs%W(i,j,k)=0.0_WP
                  ! Liquid variables
                  if (fs%VF(i,j,k).gt.0.0_WP) then
                     fs%RHOL(i,j,k)=RHOL
                     fs%PL  (i,j,k)=1.0_WP/GammaL+fs%RHOL(i,j,k)*Mach**2/16.0_WP*(cos(2.0_WP*cfg%xm(i))+cos(2.0_WP*cfg%ym(j)))*(cos(2.0_WP*cfg%zm(k))+2.0_WP)
                     fs%EL  (i,j,k)=(fs%PL(i,j,k)+GammaL*PinfL)/(fs%RHOL(i,j,k)*(GammaL-1.0_WP))
                  end if
                  ! Gas variables
                  if (fs%VF(i,j,k).lt.1.0_WP) then
                     fs%RHOG(i,j,k)=RHOG
                     fs%PG  (i,j,k)=1.0_WP/GammaG+fs%RHOG(i,j,k)*Mach**2/16.0_WP*(cos(2.0_WP*cfg%xm(i))+cos(2.0_WP*cfg%ym(j)))*(cos(2.0_WP*cfg%zm(k))+2.0_WP)
                     fs%EG  (i,j,k)=(fs%PG(i,j,k)+GammaG*PinfG)/(fs%RHOG(i,j,k)*(GammaG-1.0_WP))
                  end if
               end do
            end do
         end do
         ! Build PLIC interface
         call fs%build_interface()
         ! Initialize conserved variables
         fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
         fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
         fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%EL
         fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%EG
         call fs%get_momentum()
         ! Rebuild primitive variables
         call fs%get_primitive()
         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
      end block initial_conditions
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='TaylorGreen')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',fs%VF)
         call ens_out%add_scalar('RHO',fs%RHO)
         call ens_out%add_scalar('E',fs%E)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('Mach',Ma)
         ! Create surface mesh for PLIC
         smesh=surfmesh(nvar=0,name='plic')
         call fs%update_surfmesh(smesh)
         call ens_out%add_surface('plic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
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
         call mfile%add_column(fs%RHOLmax,'max(RHOL)')
         call mfile%add_column(fs%RHOLmin,'min(RHOL)')
         call mfile%add_column(fs%ELmax  ,'max(EL)'  )
         call mfile%add_column(fs%ELmin  ,'min(EL)'  )
         call mfile%add_column(fs%PLmax  ,'max(PL)'  )
         call mfile%add_column(fs%PLmin  ,'min(PL)'  )
         call mfile%add_column(fs%TLmax  ,'max(TL)'  )
         call mfile%add_column(fs%TLmin  ,'min(TL)'  )
         call mfile%add_column(fs%RHOGmax,'max(RHOG)')
         call mfile%add_column(fs%RHOGmin,'min(RHOG)')
         call mfile%add_column(fs%EGmax  ,'max(EG)'  )
         call mfile%add_column(fs%EGmin  ,'min(EG)'  )
         call mfile%add_column(fs%PGmax  ,'max(PG)'  )
         call mfile%add_column(fs%PGmin  ,'min(PG)'  )
         call mfile%add_column(fs%TGmax  ,'max(TG)'  )
         call mfile%add_column(fs%TGmin  ,'min(TG)'  )
         call mfile%add_column(fs%VFmax  ,'VFmax'    )
         call mfile%add_column(fs%VFmin  ,'VFmin'    )
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
         call consfile%add_column(fs%VFint  ,'Volume')
         call consfile%add_column(fs%Qint(1),'Liquid mass')
         call consfile%add_column(fs%Qint(2),'Gas mass')
         call consfile%add_column(fs%Qint(3),'Liquid energy')
         call consfile%add_column(fs%Qint(4),'Gas energy')
         call consfile%add_column(RHOSLint,'Liquid Entropy')
         call consfile%add_column(RHOSGint,'Gas Entropy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(RHOKint,'Kinetic energy')
         !call consfile%add_column(DilDiss,'Dilatation dissipation')
         !call consfile%add_column(SolDiss,'Solenoidal dissipation')
         call consfile%write()
      end block create_monitor
      
   contains
      !> Function that defines a level set function for a sphere
      function levelset_drop(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=radius-sqrt(sum((xyz-center)**2))
      end function levelset_drop
      !> Function that defines a level set function for a slab
      function levelset_slab(xyz,t) result(G)
         use mathtools, only: Pi
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=1.0_WP-abs(xyz(1)-Pi)
      end function levelset_slab
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember conserved variables
         fs%Qold=fs%Q
         
         ! Remember phasic quantities
         fs%RHOLold=fs%RHOL; fs%ELold=fs%EL
         fs%RHOGold=fs%RHOG; fs%EGold=fs%EG
         
         ! Remember volume moments and interface
         fs%VFold=fs%VF
         fs%BLold=fs%BL
         fs%BGold=fs%BG
         copy_plic_to_old: block
            use irl_fortran_interface, only: copy
            integer :: i,j,k
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
               call copy(fs%PLICold(i,j,k),fs%PLIC(i,j,k))
            end do; end do; end do
         end block copy_plic_to_old
         
         ! Tag cells for semi-Lagrangian transport
         call fs%SLtag()
         
         ! Perform first semi-Lagrangian transport step =====================================================
         call fs%SLstep(dt=0.5_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         !call fs%build_interface()
         
         ! First RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,1))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Second RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,2))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Perform second semi-Lagrangian transport step ====================================================
         call fs%SLstep(dt=1.0_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         call fs%build_interface()
         
         ! Third RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,3))
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Fourth RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,4))
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Increment Q with SL terms
         call fs%SLincrement()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)
         
         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call fs%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
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
      deallocate(dQdt,Ui,Vi,Wi,Ma)
   end subroutine simulation_final
   
   
end module simulation
