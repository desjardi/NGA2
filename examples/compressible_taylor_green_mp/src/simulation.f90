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
   type(monitor) :: mfile,cflfile,consfile
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ma
   
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
   
   
   !> P=EOS(RHO,I) for liquid
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
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
   !> S=f(RHO,P) for liquid
   pure real(WP) function get_SL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_SL=CvL*log((P+PinfL)/RHO**GammaL)
   end function get_SL
   
   
   !> P=EOS(RHO,I) for gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
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
   !> S=f(RHO,P) for gas
   pure real(WP) function get_SG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_SG=CvG*log((P+PinfG)/RHO**GammaG)
   end function get_SG
   
   
   !> Mechanical relaxation model
   subroutine P_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1)-0.5_WP*((Q(5)/(Q(1)+Q(2)))**2+(Q(6)/(Q(1)+Q(2)))**2+(Q(7)/(Q(1)+Q(2)))**2))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2)-0.5_WP*((Q(5)/(Q(1)+Q(2)))**2+(Q(6)/(Q(1)+Q(2)))**2+(Q(7)/(Q(1)+Q(2)))**2))
      ! Handle limit cases
      if (PL.le.-PinfL) then; VF=0.0_WP; Q(4)=Q(3)+Q(4); Q(3)=0.0_WP; return; end if
      if (PG.le.-PinfG) then; VF=1.0_WP; Q(3)=Q(3)+Q(4); Q(4)=0.0_WP; return; end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)**2
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG)**2
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      coeffL=(GammaL-1.0_WP)*Pint+2.0_WP*GammaL*PinfL
      coeffG=(GammaG-1.0_WP)*Pint+2.0_WP*GammaG*PinfG
      a=1.0_WP+GammaG*VF+GammaL*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+GammaG)*VF*PL-(1.0_WP+GammaL)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Get equilibrium volume fraction
      VFeq=VF*((gammaL-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+gammaL)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax


   !> Pseudo single-phase mechanical relaxation
   !subroutine P_relax(VF,Q)
   !   implicit none
   !   real(WP),                intent(inout) :: VF
   !   real(WP), dimension(1:), intent(inout) :: Q
   !   real(WP) :: E
   !   VF=Q(1)/(Q(1)+Q(2))
   !   E=Q(3)+Q(4)
   !   Q(3)=E*VF
   !   Q(4)=E*(1.0_WP-VF)
   !end subroutine P_relax
   
   
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
      
      ! Create multipgase compressible flow solver
      create_velocity_solver: block
         ! Initialize solver with required thermodynamic functions
         call fs%initialize(cfg=cfg,getPL=get_PL,getCL=get_CL,getPG=get_PG,getCG=get_CG,name='Compressible NS')
         ! Provide pressure relaxation model
         fs%Prelax=>P_relax
         ! Provide entropy calculation functions
         fs%getSL=>get_SL; fs%getSG=>get_SG
         ! Provide temperature calculation functions
         !fs%getTL=>get_TL; fs%getTG=>get_TG
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dQdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
         allocate(Ma  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Prepare initial conditions
      initial_conditions: block
         use mms_geom,     only: initialize_volume_moments
         use mpcomp_class, only: VFlo
         use random,       only: random_uniform
         integer :: i,j,k
         ! Calculate domain center
         center=[cfg%x(cfg%imin),cfg%y(cfg%jmin),cfg%z(cfg%kmin)]+0.5_WP*[cfg%x(cfg%imax+1)-cfg%x(cfg%imin),cfg%y(cfg%jmax+1)-cfg%y(cfg%jmin),cfg%z(cfg%kmax+1)-cfg%z(cfg%kmin)]
         ! Initialize primary variables
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Volume moments for a droplet
                  call initialize_volume_moments(lo=[cfg%x(i),cfg%y(j),cfg%z(k)],hi=[cfg%x(i+1),cfg%y(j+1),cfg%z(k+1)],&
                  levelset=levelset_slab,time=0.0_WP,level=4,VFlo=VFlo,VF=fs%VF(i,j,k),BL=fs%BL(:,i,j,k),BG=fs%BG(:,i,j,k))
                  ! Mixture velocity
                  fs%U(i,j,k)=+Mach*sin(cfg%xm(i))*cos(cfg%ym(j))*cos(cfg%zm(k))
                  fs%V(i,j,k)=-Mach*cos(cfg%xm(i))*sin(cfg%ym(j))*cos(cfg%zm(k))
                  fs%W(i,j,k)=0.0_WP
                  ! Liquid variables
                  if (fs%VF(i,j,k).gt.0.0_WP) then
                     fs%RHOL(i,j,k)=RHOL
                     fs%PL  (i,j,k)=1.0_WP/GammaL+fs%RHOL(i,j,k)*Mach**2/16.0_WP*(cos(2.0_WP*cfg%xm(i))+cos(2.0_WP*cfg%ym(j)))*(cos(2.0_WP*cfg%zm(k))+2.0_WP)
                     fs%IL  (i,j,k)=(fs%PL(i,j,k)+GammaL*PinfL)/(fs%RHOL(i,j,k)*(GammaL-1.0_WP))
                  end if
                  ! Gas variables
                  if (fs%VF(i,j,k).lt.1.0_WP) then
                     fs%RHOG(i,j,k)=RHOG
                     fs%PG  (i,j,k)=1.0_WP/GammaG+fs%RHOG(i,j,k)*Mach**2/16.0_WP*(cos(2.0_WP*cfg%xm(i))+cos(2.0_WP*cfg%ym(j)))*(cos(2.0_WP*cfg%zm(k))+2.0_WP)
                     fs%IG  (i,j,k)=(fs%PG(i,j,k)+GammaG*PinfG)/(fs%RHOG(i,j,k)*(GammaG-1.0_WP))
                  end if
               end do
            end do
         end do
         ! Build PLIC interface
         call fs%build_interface()
         ! Build phasic total energies
         where (fs%VF.gt.0.0_WP) fs%EL=fs%IL+0.5_WP*(fs%U**2+fs%V**2+fs%W**2)
         where (fs%VF.lt.1.0_WP) fs%EG=fs%IG+0.5_WP*(fs%U**2+fs%V**2+fs%W**2)
         ! Initialize conserved variables
         fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
         fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
         fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%EL
         fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%EG
         fs%RHO=fs%Q(:,:,:,1)+fs%Q(:,:,:,2)
         fs%Q(:,:,:,5)=fs%RHO*fs%U
         fs%Q(:,:,:,6)=fs%RHO*fs%V
         fs%Q(:,:,:,7)=fs%RHO*fs%W
         ! Communicate conserved variables (not needed in general, but allows 2D runs without changing loop above...)
         do i=1,fs%nQ; call fs%cfg%sync(fs%Q(:,:,:,i)); end do
         ! Rebuild primitive variables
         call fs%get_primitive()
         ! Compute local Mach number
         Ma=sqrt(fs%U**2+fs%V**2+fs%W**2)/fs%C
      end block initial_conditions
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='TaylorGreen')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
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
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%RHOKLint,'Liquid KE')
         call consfile%add_column(fs%RHOKGint,'Gas KE')
         call consfile%add_column(fs%RHOSLint,'Liquid entropy')
         call consfile%add_column(fs%RHOSGint,'Gas entropy')
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
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=1.0_WP-abs(xyz(1)-center(1))
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
         fs%RHOLold=fs%RHOL; fs%ILold=fs%IL
         fs%RHOGold=fs%RHOG; fs%IGold=fs%IG
         
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
         ! Enforce mechanical equilibrium
         call fs%relax_pressure()
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Compute local Mach number
         Ma=sqrt(fs%U**2+fs%V**2+fs%W**2)/fs%C
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call fs%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(dQdt,Ma)
   end subroutine simulation_final
   
   
end module simulation
