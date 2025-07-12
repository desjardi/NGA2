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
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc
   
   !> Equations of state
   real(WP) :: PinfL,GammaL,CvL
   real(WP) :: PinfG,GammaG,CvG
   
   !> Flow parameters
   real(WP) :: Ms,Xs
   real(WP) :: rho1,p1,u1,M1
   real(WP) :: rho2,p2,u2,M2
   real(WP) :: rho_ratio,c_ratio
   real(WP) :: rhoL,ML
   real(WP) :: ReG,viscG,viscL,visc_ratio
   
contains
   
   
   !> Sutherland's law for viscosity as a function of temperature
   !subroutine get_visc()
   !   implicit none
   !   integer :: i,j,k
   !   do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
   !      fs%visc(i,j,k)=Reynolds**(-1.0_WP)*(1.4042_WP*fs%T(i,j,k)**1.5_WP)/(fs%T(i,j,k)+0.4042_WP)
   !   end do; end do; end do
   !end subroutine get_visc
   
   
   !> Function that returns a smooth Heaviside of thickness delta
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      ! Goes from 0 to 1 as x goes from begative to positive
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock
   
   
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
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2))
      ! Handle limit cases - should mass/energy be tranasfered or lost? - this should probably never happen...
      if (PL.le.-PinfL) then
         print*,"****************** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP; return
      end if
      if (PG.le.-PinfG) then
         print*,"****************** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; return
      end if
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
      ! Last debugging check... Probably should never happen...
      if (Peq.lt.-PinfG) print*,"****************** NEGATIVE PRESSURE! - time",time%t,"VFeq",VFeq,"Peq",Peq
   end subroutine P_relax
   
   
   !> Thermo-mechanical relaxation model
   subroutine PT_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2))
      ! Handle limit cases - should mass/energy be tranasfered or lost? - this should probably never happen...
      if (PL.le.-PinfL) then
         print*,"****************** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP; return
      end if
      if (PG.le.-PinfG) then
         print*,"****************** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; return
      end if
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
      ! ================= Second step for thermal relaxation =================
      ! Setup quadratic problem
      a=Q(1)*CvL+Q(2)*CvG
      b=Q(1)*CvL*(GammaL*PinfL+PinfG)+Q(2)*CvG*(GammaG*PinfG+PinfL)-sum(Q(3:4))*(Q(1)*CvL*(GammaL-1.0_WP)+Q(2)*CvG*(GammaG-1.0_WP))
      d=(Q(1)*CvL*GammaL+Q(2)*CvG*GammaG)*PinfL*PinfG-sum(Q(3:4))*(Q(1)*CvL*(GammaL-1.0_WP)*PinfG+Q(2)*CvG*(GammaG-1.0_WP)*PinfL)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Get equilibrium volume fraction
      VFeq=Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)/(Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)+Q(2)*CvG*(GammaG-1.0_WP)*(Peq+PinfL))
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-PinfL); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,-PinfG); end if
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*(Peq+GammaL*PinfL)/(GammaL-1.0_WP)
      Q(4)=(1.0_WP-VFeq)*(Peq+GammaG*PinfG)/(GammaG-1.0_WP)
      VF=VFeq
      ! Last debugging check... Probably should never happen...
      if (Peq.lt.-PinfG) print*,"****************** NEGATIVE PRESSURE! - time",time%t,"VFeq",VFeq,"Peq",Peq
   end subroutine PT_relax
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Initialize eos and flow parameters
      initialize_parameters: block
         use string,   only: str_long
         use messager, only: log
         character(str_long) :: message
         ! Set PinfG to zero
         PinfG=0.0_WP
         ! Read in Gammas
         call param_read('Liquid gamma',GammaL)
         call param_read('Gas gamma'   ,GammaG)
         ! Read in shock Mach number and location
         call param_read('Shock Mach number',Ms)
         call param_read('Shock location',Xs)
         ! First generate static shock with normalized pre-shock conditions
         M1=Ms
         rho1=1.0_WP
         rho2=rho1*(GammaG+1.0_WP)*M1**2/((GammaG-1.0_WP)*M1**2+2.0_WP)
         p1=0.25_WP*rho1/GammaG*((GammaG+1.0_WP)*M1/(M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
         p2=p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(M1**2-1.0_WP)+1.0_WP)
         u1=M1*sqrt(GammaG*p1/rho1)
         u2=u1*rho1/rho2
         ! Now shift frame of reference to obtain moving shock
         u2=abs(u2-u1); M2=u2/sqrt(GammaG*p2/rho2); u1=0.0_WP; M1=u1/sqrt(GammaG*p1/rho1)
         ! Read in density ratio and use it to set liquid density
         call param_read('Density ratio',rho_ratio); rhoL=rho_ratio*rho1
         ! Read in liquid Mach number and use it to set PinfL
         !call param_read('Liquid Mach number',ML)
         !PinfL=(u2/ML)**2*rhoL/GammaL-p1
         !c_ratio=sqrt(GammaL*(p1+PinfL)/rhoL)/sqrt(GammaG*p1/rho1)
         ! Read in sound speed ratio and use it to set PinfL
         call param_read('Sound speed ratio',c_ratio)
         PinfL=p1*(rho_ratio*c_ratio**2*GammaG/GammaL-1.0_WP)
         ML=u2/sqrt(GammaL*(p1+PinfL)/rhoL)
         ! Set heat capacities corresponding to a normalized pre-shock and liquid temperature
         CvL=(p1+PinfL)/(rhoL*(GammaL-1.0_WP))
         CvG=(p1+PinfG)/(rho1*(GammaG-1.0_WP))
         ! Viscous parameters
         call param_read('Gas Reynolds number',ReG); viscG=rho1*1.0_WP*u2/ReG 
         call param_read('Viscosity ratio',visc_ratio); viscL=visc_ratio*viscG
         ! Output case info
         if (cfg%amRoot) then
            write(message,'("[Liquid EOS] => Gamma=",es12.5)') GammaL; call log(message)
            write(message,'("[Liquid EOS] =>  Pinf=",es12.5)')  PinfL; call log(message)
            write(message,'("[Liquid EOS] =>    Cv=",es12.5)')    CvL; call log(message)
            write(message,'("[Gas EOS]    => Gamma=",es12.5)') GammaG; call log(message)
            write(message,'("[Gas EOS]    =>    Cv=",es12.5)')    CvG; call log(message)
            write(message,'("[Shock Mach number]     =>     Ms=",es12.5)')     Ms; call log(message)
            write(message,'("[Pre -shock conditions] =>   rho1=",es12.5)')   rho1; call log(message)
            write(message,'("[Pre -shock conditions] =>     p1=",es12.5)')     p1; call log(message)
            write(message,'("[Pre -shock conditions] =>     u1=",es12.5)')     u1; call log(message)
            write(message,'("[Pre -shock conditions] =>     M1=",es12.5)')     M1; call log(message)
            write(message,'("[Post-shock conditions] =>   rho2=",es12.5)')   rho2; call log(message)
            write(message,'("[Post-shock conditions] =>     p2=",es12.5)')     p2; call log(message)
            write(message,'("[Post-shock conditions] =>     u2=",es12.5)')     u2; call log(message)
            write(message,'("[Post-shock conditions] =>     M2=",es12.5)')     M2; call log(message)
            write(message,'("[Liquid Mach number] =>        ML=",es12.5)')     Ml; call log(message)
            write(message,'("[Density ratio]      => rhoL/rho1=",es12.5)') rho_ratio; call log(message)
            write(message,'("[Sound speed ratio]  =>     cl/c1=",es12.5)')   c_ratio; call log(message)
            write(message,'("[Gas Reynolds]     =>     ReG=",es12.5)')        ReG; call log(message)
            write(message,'("[Viscosity ratio]  => muL/muG=",es12.5)') visc_ratio; call log(message)
            write(message,'("[Gas    viscosity] =>     muG=",es12.5)')      viscG; call log(message)
            write(message,'("[Liquid viscosity] =>     muL=",es12.5)')      viscL; call log(message)
         end if
      end block initialize_parameters
      
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
         ! Provide relaxation model
         fs%relax=>P_relax
         ! Provide entropy calculation functions
         fs%getSL=>get_SL; fs%getSG=>get_SG
         ! Provide temperature calculation functions
         fs%getTL=>get_TL; fs%getTG=>get_TG
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dQdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
         allocate(beta(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(visc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ma(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Prepare initial conditions
      initial_conditions: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use mms_geom,              only: initialize_volume_moments
         use mpcomp_class,          only: VFlo
         integer :: i,j,k
         ! Initialize primary variables
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Initialize liquid volume to zero and set corresponding PLIC
                  fs%VF(i,j,k)=0.0_WP; fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]; fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  call setNumberOfPlanes(fs%PLIC(i,j,k),1); call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,fs%VF(i,j,k)-0.5_WP))
                  ! Not set volume moments for a droplet or a slab
                  call initialize_volume_moments(lo=[cfg%x(i),cfg%y(j),cfg%z(k)],hi=[cfg%x(i+1),cfg%y(j+1),cfg%z(k+1)],&
                  levelset=levelset_drop,time=0.0_WP,level=4,VFlo=VFlo,VF=fs%VF(i,j,k),BL=fs%BL(:,i,j,k),BG=fs%BG(:,i,j,k))
                  ! Initialize mixture velocity to normal shock
                  fs%U(i,j,k)=u2*Hshock(Xs-cfg%x(i),delta=0.5_WP*fs%dx)
                  fs%V(i,j,k)=0.0_WP
                  fs%W(i,j,k)=0.0_WP
                  ! Gas variables
                  if (fs%VF(i,j,k).lt.1.0_WP) then
                     fs%RHOG(i,j,k)=rho1+(rho2-rho1)*Hshock(Xs-cfg%xm(i),delta=0.5_WP*fs%dx)
                     fs%PG  (i,j,k)=p1  +(p2  -p1  )*Hshock(Xs-cfg%xm(i),delta=0.5_WP*fs%dx)
                     fs%IG  (i,j,k)=(fs%PG(i,j,k)+GammaG*PinfG)/(fs%RHOG(i,j,k)*(GammaG-1.0_WP))
                  end if
                  ! Liquid variables
                  if (fs%VF(i,j,k).gt.0.0_WP) then
                     fs%RHOL(i,j,k)=rhoL
                     fs%PL  (i,j,k)=p1
                     fs%IL  (i,j,k)=(fs%PL(i,j,k)+GammaL*PinfL)/(fs%RHOL(i,j,k)*(GammaL-1.0_WP))
                  end if
               end do
            end do
         end do
         ! Build PLIC interface
         call fs%build_interface()
         ! Initialize conserved variables
         fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
         fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
         fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
         fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG
         call fs%get_momentum()
         ! Communicate conserved variables (not needed in general, but allows 2D runs without changing loop above...)
         do i=1,fs%nQ; call fs%cfg%sync(fs%Q(:,:,:,i)); end do
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
         ens_out=ensight(cfg=cfg,name='ShockDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',fs%VF)
         call ens_out%add_scalar('RHOL',fs%RHOL)
         call ens_out%add_scalar('RHOG',fs%RHOG)
         call ens_out%add_scalar('IL',fs%IL)
         call ens_out%add_scalar('IG',fs%IG)
         call ens_out%add_scalar('PL',fs%PL)
         call ens_out%add_scalar('PG',fs%PG)
         call ens_out%add_scalar('Mach',Ma)
         call ens_out%add_scalar('beta',beta)
         call ens_out%add_scalar('visc',visc)
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
         call mfile%add_column(fs%ILmax  ,'max(IL)'  )
         call mfile%add_column(fs%ILmin  ,'min(IL)'  )
         call mfile%add_column(fs%PLmax  ,'max(PL)'  )
         call mfile%add_column(fs%PLmin  ,'min(PL)'  )
         call mfile%add_column(fs%TLmax  ,'max(TL)'  )
         call mfile%add_column(fs%TLmin  ,'min(TL)'  )
         call mfile%add_column(fs%RHOGmax,'max(RHOG)')
         call mfile%add_column(fs%RHOGmin,'min(RHOG)')
         call mfile%add_column(fs%IGmax  ,'max(IG)'  )
         call mfile%add_column(fs%IGmin  ,'min(IG)'  )
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
      !> Level set function for a sphere of unity diameter centered at (0,0,0)
      function levelset_drop(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-sqrt(sum(xyz**2))
      end function levelset_drop
      !> Level set function for a slab of unity width centered at x=0
      function levelset_slab(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-abs(xyz(1))
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
         fs%RHOLold=fs%RHOL; fs%ILold=fs%IL; fs%PLold=fs%PL
         fs%RHOGold=fs%RHOG; fs%IGold=fs%IG; fs%PGold=fs%PG
         
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
         
         ! Prepare SGS viscosity models
         call fs%get_viscartif(dt=time%dt,beta=beta)
         visc=0.0_WP !call fs%get_vreman(dt=time%dt,visc=visc)
         mixture_viscosity: block
            integer  :: i,j,k
            real(WP) :: Lvof,Lrho,Gvof,Grho
            real(WP) :: Lvisc,Gvisc,Lbeta,Gbeta
            real(WP), parameter :: eps=1.0e-15_WP
            do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_-1; do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_-1; do i=fs%cfg%imino_+1,fs%cfg%imaxo_-1
               ! Create smooth mass info distribution
               Lvof=sum(       fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
               Gvof=sum(1.0_WP-fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
               Lrho=sum(       fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,1))/(Lvof+eps)
               Grho=sum(       fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,2))/(Gvof+eps)
               ! Harmonic average of VISC
               Lvisc=Lrho*viscL; Gvisc=Grho*viscG; fs%VISC(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lvisc,eps)+Gvof/max(Gvisc,eps))
               ! Harmonic average of BETA
               Lbeta=Lrho*beta(i,j,k); Gbeta=Grho*beta(i,j,k); fs%BETA(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lbeta,eps)+Gvof/max(Gbeta,eps))
            end do; end do; end do
         end block mixture_viscosity
         
         ! Perform first semi-Lagrangian transport step =====================================================
         call fs%SLstep(dt=0.5_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         !call fs%build_interface()
         
         ! First RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,1))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Second RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,2))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
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
         fs%Q=fs%Q+fs%SLdQ
         ! Recompute primitive variables
         call fs%get_primitive()
         
         ! Fourth RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,4))
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
         ! Apply user-provided relaxation model
         call fs%apply_relax()
         ! Recompute primitive variables
         call fs%get_primitive()
         ! Apply Neumann condition at the outflow
         neumann_outflow: block
            use irl_fortran_interface, only: setPlane
            integer :: i,j,k
            ! Apply clipped Neumann on primitive variables
            if (fs%cfg%iproc.eq.fs%cfg%npx) then
               do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imax+1,fs%cfg%imaxo
                  ! Copy primitive variables
                  fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imax,j,k)
                  fs%PL  (i,j,k)=fs%PL  (fs%cfg%imax,j,k)
                  fs%IL  (i,j,k)=fs%IL  (fs%cfg%imax,j,k)
                  fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imax,j,k)
                  fs%PG  (i,j,k)=fs%PG  (fs%cfg%imax,j,k)
                  fs%IG  (i,j,k)=fs%IG  (fs%cfg%imax,j,k)
                  fs%U  (i,j,k)=max(fs%U(fs%cfg%imax,j,k),0.0_WP)
                  fs%V   (i,j,k)=fs%V   (fs%cfg%imax,j,k)
                  fs%W   (i,j,k)=fs%W   (fs%cfg%imax,j,k)
                  fs%VF  (i,j,k)=fs%VF  (fs%cfg%imax,j,k)
                  ! Also adjust interface data
                  call setPlane(fs%PLIC(i,j,k),0,[1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%cfg%dx(i)*fs%VF(i,j,k))
                  fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
               end do; end do; end do
            end if
            ! Rebuild conserved quantities
            fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
            fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
            fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
            fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG
            call fs%get_momentum()
         end block neumann_outflow
         
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
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(dQdt,Ui,Vi,Wi,Ma,beta,visc)
   end subroutine simulation_final
   
   
end module simulation
