!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,       only: WP
   use shockdrop_class, only: shockdrop
   use ffshock_class,   only: ffshock
   use coupler_class,   only: coupler
   use event_class,     only: event
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Shock-drop simulation
   type(shockdrop) :: sd
   
   !> Far-field shock simulation
   type(ffshock) :: ff
   
   !> Couplers between domains
   type(coupler) :: sd2ff
   type(coupler) :: ff2sd
   
   !> Ensight output event
   type(event) :: ens_evt
   
   !> Remeshing event
   type(event) :: remesh_evt
   
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
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! ================ Handle gas flotsams ================
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
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
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=VF*((gammaL-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+gammaL)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
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
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)/(Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)+Q(2)*CvG*(GammaG-1.0_WP)*(Peq+PinfL))
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-PinfL); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,-PinfG); end if
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*(Peq+GammaL*PinfL)/(GammaL-1.0_WP)
      Q(4)=(1.0_WP-VFeq)*(Peq+GammaG*PinfG)/(GammaG-1.0_WP)
      VF=VFeq
   end subroutine PT_relax
   
   
   !> Solver initialization
   subroutine simulation_init
      implicit none
      
      ! Initialize eos and flow parameters - all cores
      initialize_parameters: block
         use string,   only: str_long
         use messager, only: log
         use parallel, only: amRoot
         use param,    only: param_read
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
         if (amRoot) then
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
      
      ! Setup shock-drop simulation - all cores
      setup_sd: block
         use param,    only: param_read
         use parallel, only: group
         integer , dimension(3) :: meshsize,partition
         real(WP), dimension(3) :: X0
         real(WP) :: dx
         ! Read in mesh size and desired partition
         call param_read('Shock-drop dx',dx)
         call param_read('Shock-drop nx',meshsize)
         call param_read('Shock-drop partition',partition)
         X0=-0.5_WP*real(meshsize,WP)*dx !< This assumes that the domain is centered on (0,0,0)
         ! Initialize the shock-drop solver
         call sd%init(dx=dx,meshsize=meshsize,startloc=X0,group=group,partition=partition)
         ! Provide relaxation and thermodynamic models
         sd%fs%relax=>P_relax
         sd%fs%getPL=>get_PL; sd%fs%getCL=>get_CL; sd%fs%getSL=>get_SL; sd%fs%getTL=>get_TL
         sd%fs%getPG=>get_PG; sd%fs%getCG=>get_CG; sd%fs%getSG=>get_SG; sd%fs%getTG=>get_TG
         ! Set time integration parameters
         call param_read('Shock-drop dt',sd%time%dtmax); sd%time%dt=sd%time%dtmax
         call param_read('Shock-drop CFL',sd%time%cflmax)
         call param_read('Max time',sd%time%tmax)
         ! We finally need to transfer our viscosities explicitly...
         sd%cst_viscL=viscL; sd%cst_viscG=viscG
      end block setup_sd
      
      ! Generate initial conditions for shock-drop problem
      initialize_sd: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use mms_geom,              only: initialize_volume_moments
         use mpcomp_class,          only: VFlo
         integer :: i,j,k
         ! Initialize primary variables
         do k=sd%cfg%kmino_,sd%cfg%kmaxo_; do j=sd%cfg%jmino_,sd%cfg%jmaxo_; do i=sd%cfg%imino_,sd%cfg%imaxo_
            ! Initialize liquid volume to zero and set corresponding PIC
            sd%fs%VF(i,j,k)=0.0_WP; sd%fs%BL(:,i,j,k)=[sd%fs%cfg%xm(i),sd%fs%cfg%ym(j),sd%fs%cfg%zm(k)]; sd%fs%BG(:,i,j,k)=[sd%fs%cfg%xm(i),sd%fs%cfg%ym(j),sd%fs%cfg%zm(k)]
            call setNumberOfPlanes(sd%fs%PLIC(i,j,k),1); call setPlane(sd%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,sd%fs%VF(i,j,k)-0.5_WP))
            ! Not set volume moments for a droplet or a slab
            call initialize_volume_moments(lo=[sd%fs%cfg%x(i),sd%fs%cfg%y(j),sd%fs%cfg%z(k)],hi=[sd%fs%cfg%x(i+1),sd%fs%cfg%y(j+1),sd%fs%cfg%z(k+1)],&
            levelset=levelset_drop,time=0.0_WP,level=5,VFlo=VFlo,VF=sd%fs%VF(i,j,k),BL=sd%fs%BL(:,i,j,k),BG=sd%fs%BG(:,i,j,k))
            ! Initialize mixture velocity to normal shock
            sd%fs%U(i,j,k)=u2*Hshock(Xs-sd%fs%cfg%x(i),delta=0.5_WP*sd%fs%dx)
            sd%fs%V(i,j,k)=0.0_WP
            sd%fs%W(i,j,k)=0.0_WP
            ! Gas variables
            if (sd%fs%VF(i,j,k).lt.1.0_WP) then
               sd%fs%RHOG(i,j,k)=rho1+(rho2-rho1)*Hshock(Xs-sd%fs%cfg%xm(i),delta=0.5_WP*sd%fs%dx)
               sd%fs%PG  (i,j,k)=p1  +(p2  -p1  )*Hshock(Xs-sd%fs%cfg%xm(i),delta=0.5_WP*sd%fs%dx)
               sd%fs%IG  (i,j,k)=(sd%fs%PG(i,j,k)+GammaG*PinfG)/(sd%fs%RHOG(i,j,k)*(GammaG-1.0_WP))
            end if
            ! Liquid variables
            if (sd%fs%VF(i,j,k).gt.0.0_WP) then
               sd%fs%RHOL(i,j,k)=rhoL
               sd%fs%PL  (i,j,k)=p1
               sd%fs%IL  (i,j,k)=(sd%fs%PL(i,j,k)+GammaL*PinfL)/(sd%fs%RHOL(i,j,k)*(GammaL-1.0_WP))
            end if
         end do; end do; end do
         ! Build PLIC interface
         call sd%fs%build_interface()
         ! Initialize conserved variables
         sd%fs%Q(:,:,:,1)=        sd%fs%VF *sd%fs%RHOL
         sd%fs%Q(:,:,:,2)=(1.0_WP-sd%fs%VF)*sd%fs%RHOG
         sd%fs%Q(:,:,:,3)= sd%fs%Q(:,:,:,1)*sd%fs%IL
         sd%fs%Q(:,:,:,4)= sd%fs%Q(:,:,:,2)*sd%fs%IG
         call sd%fs%get_momentum()
         ! Communicate conserved variables (not needed in general, but allows 2D runs without changing loop above...)
         do i=1,sd%fs%nQ; call sd%fs%cfg%sync(sd%fs%Q(:,:,:,i)); end do
         ! Rebuild primitive variables
         call sd%fs%get_primitive()
         ! Interpolate velocity
         call sd%fs%interp_vel(sd%Ui,sd%Vi,sd%Wi)
         ! Compute local Mach number
         sd%Ma=sqrt(sd%Ui**2+sd%Vi**2+sd%Wi**2)/sd%fs%C
         ! Perform monitoring
         call sd%output_monitor()
      end block initialize_sd
      
      ! Setup far-field shock simulation - all cores
      setup_ff: block
         use param,    only: param_read
         use parallel, only: group
         integer , dimension(3) :: meshsize,partition
         real(WP), dimension(3) :: X0
         real(WP) :: dx
         ! Read in mesh size and desired partition
         call param_read('Farfield dx',dx)
         call param_read('Farfield nx',meshsize)
         call param_read('Farfield partition',partition)
         X0=-0.5_WP*real(meshsize,WP)*dx      !< This assumes that the domain is centered on (0,0,0)
         call param_read('Farfield X0',X0(1)) !< This shifts the domain in x based on user input
         ! Initialize the farfield solver
         call ff%init(dx=dx,meshsize=meshsize,startloc=X0,group=group,partition=partition)
         ! Provide thermodynamic model
         ff%fs%getP=>get_PG; ff%fs%getC=>get_CG; ff%fs%getS=>get_SG; ff%fs%getT=>get_TG
         ! We finally need to transfer our viscosity explicitly...
         ff%cst_visc=viscG
      end block setup_ff
      
      ! Generate initial conditions for far-field shock problem
      initialize_ff: block
         integer :: i,j,k
         ! Initialize primary variables to normal shock
         do k=ff%cfg%kmino_,ff%cfg%kmaxo_; do j=ff%cfg%jmino_,ff%cfg%jmaxo_; do i=ff%cfg%imino_,ff%cfg%imaxo_
            ff%fs%U(i,j,k)  =u2*Hshock(Xs-ff%fs%cfg%x(i),delta=0.5_WP*ff%fs%dx)
            ff%fs%V(i,j,k)  =0.0_WP
            ff%fs%W(i,j,k)  =0.0_WP
            ff%fs%Q(i,j,k,1)=rho1+(rho2-rho1)*Hshock(Xs-ff%fs%cfg%xm(i),delta=0.5_WP*ff%fs%dx)
            ff%fs%P(i,j,k)  =p1  +(p2  -p1  )*Hshock(Xs-ff%fs%cfg%xm(i),delta=0.5_WP*ff%fs%dx)
            ff%fs%I(i,j,k)  =(ff%fs%P(i,j,k)+GammaG*PinfG)/(ff%fs%Q(i,j,k,1)*(GammaG-1.0_WP))
         end do; end do; end do
         ! Initialize conserved variables
         ff%fs%Q(:,:,:,2)=ff%fs%Q(:,:,:,1)*ff%fs%I
         call ff%fs%get_momentum()
         ! Rebuild primitive variables
         call ff%fs%get_primitive()
         ! Interpolate velocity
         call ff%fs%interp_vel(ff%Ui,ff%Vi,ff%Wi)
         ! Compute local Mach number
         ff%Ma=sqrt(ff%Ui**2+ff%Vi**2+ff%Wi**2)/ff%fs%C
         ! Perform monitoring
         call ff%output_monitor()
      end block initialize_ff
      
      ! Create couplers
      create_couplers: block
         use parallel, only: group
         sd2ff=coupler(src_grp=group,dst_grp=group,name='sd2ff'); call sd2ff%set_src(sd%cfg); call sd2ff%set_dst(ff%cfg); call sd2ff%initialize()
         ff2sd=coupler(src_grp=group,dst_grp=group,name='ff2sd'); call ff2sd%set_src(ff%cfg); call ff2sd%set_dst(sd%cfg); call ff2sd%initialize()
      end block create_couplers
      
      ! Initialize Ensight output event and perform initial Ensight output
      initialize_ensight: block
         use param, only: param_read
         ens_evt=event(time=sd%time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         if (ens_evt%occurs()) then
            call sd%output_ensight(t=sd%time%t)
            call ff%output_ensight(t=sd%time%t)
         end if
      end block initialize_ensight
      
      ! Initialize remeshing event
      initialize_remeshing: block
         use param, only: param_read
         remesh_evt=event(time=sd%time,name='Remeshing')
         call param_read('Remeshing period',remesh_evt%tper)
      end block initialize_remeshing
      
   contains
      !> Level set function for a sphere of unity diameter centered at (0,0,0)
      function levelset_drop(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-sqrt(sum(xyz**2))
      end function levelset_drop
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Shock-drop drives overall time integration
      do while (.not.sd%time%done())
         
         ! Handle coupling
         call couple_sd2ff()
         call couple_ff2sd()
         
         ! Advance shock-drop simulation
         call sd%step()
         
         ! Advance farfield simulation
         call ff%step(dt=sd%time%dt)
         
         ! Perform monitoring
         call sd%output_monitor()
         call ff%output_monitor()

         ! Remesh sd
         if (remesh_evt%occurs()) call remesh()
         
         ! Perform Ensight output
         if (ens_evt%occurs()) then
            call sd%output_ensight(t=sd%time%t)
            call ff%output_ensight(t=sd%time%t)
         end if
         
      end do
      
   end subroutine simulation_run
   

   !> Remesh sd to follow the drop
   subroutine remesh()
      implicit none
      type(shockdrop) :: sdnew
      type(coupler) :: sdnew2ff
      type(coupler) :: ff2sdnew
      
      ! Setup new shock-drop simulation - all cores
      setup_sdnew: block
         use param,    only: param_read
         use parallel, only: group
         real(WP), dimension(3) :: X0
         ! Shift domain based on core barycenter
         X0=[sd%cfg%x(sd%cfg%imin),sd%cfg%y(sd%cfg%jmin),sd%cfg%z(sd%cfg%kmin)]+sd%fs%dx*real([int(sd%Xcore/sd%fs%dx),int(sd%Ycore/sd%fs%dy),int(sd%Zcore/sd%fs%dz)],WP)
         ! Initialize the shock-drop solver
         call sdnew%init(dx=sd%fs%dx,meshsize=[sd%cfg%nx,sd%cfg%ny,sd%cfg%nz],startloc=X0,group=group,partition=[sd%cfg%npx,sd%cfg%npy,sd%cfg%npz])
         ! Provide relaxation and thermodynamic models
         sdnew%fs%relax=>sd%fs%relax
         sdnew%fs%getPL=>sd%fs%getPL; sdnew%fs%getCL=>sd%fs%getCL; sdnew%fs%getSL=>sd%fs%getSL; sdnew%fs%getTL=>sd%fs%getTL
         sdnew%fs%getPG=>sd%fs%getPG; sdnew%fs%getCG=>sd%fs%getCG; sdnew%fs%getSG=>sd%fs%getSG; sdnew%fs%getTG=>sd%fs%getTG
         ! Set time integration parameters
         sdnew%time%dtmax =sd%time%dtmax
         sdnew%time%dt    =sd%time%dt
         sdnew%time%cflmax=sd%time%cflmax
         sdnew%time%tmax  =sd%time%tmax
         ! We finally need to transfer our viscosities explicitly...
         sdnew%cst_viscL=sd%cst_viscL; sdnew%cst_viscG=sd%cst_viscG
      end block setup_sdnew
      
      ! Create new couplers
      setup_new_couplers: block
         use parallel, only: group
         sdnew2ff=coupler(src_grp=group,dst_grp=group,name='sd2ff'); call sdnew2ff%set_src(sdnew%cfg); call sdnew2ff%set_dst(ff%cfg); call sdnew2ff%initialize()
         ff2sdnew=coupler(src_grp=group,dst_grp=group,name='ff2sd'); call ff2sdnew%set_src(ff%cfg); call ff2sdnew%set_dst(sdnew%cfg); call ff2sdnew%initialize()
      end block setup_new_couplers
      
      ! Initialize sdnew using ff
      initialize_sdnew_from_ff: block
         integer :: i,j,k,n
         call ff2sdnew%push(ff%fs%Q(:,:,:,1)); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,2))
         call ff2sdnew%push(ff%fs%Q(:,:,:,2)); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,4))
         call ff2sdnew%push(ff%Ui);            call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%Ui)
         call ff2sdnew%push(ff%Vi);            call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%Vi)
         call ff2sdnew%push(ff%Wi);            call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%Wi)
         do k=sdnew%cfg%kmino_+1,sdnew%cfg%kmaxo_; do j=sdnew%cfg%jmino_+1,sdnew%cfg%jmaxo_; do i=sdnew%cfg%imino_+1,sdnew%cfg%imaxo_
            sdnew%fs%Q(i,j,k,5)=0.5_WP*sum(sdnew%fs%Q(i-1:i,j,k,2)*sdnew%Ui(i-1:i,j,k))
            sdnew%fs%Q(i,j,k,6)=0.5_WP*sum(sdnew%fs%Q(i,j-1:j,k,2)*sdnew%Vi(i,j-1:j,k))
            sdnew%fs%Q(i,j,k,7)=0.5_WP*sum(sdnew%fs%Q(i,j,k-1:k,2)*sdnew%Wi(i,j,k-1:k))
         end do; end do; end do
         do n=1,sdnew%fs%nQ; call sdnew%cfg%sync(sdnew%fs%Q(:,:,:,n)); end do
      end block initialize_sdnew_from_ff
      
      ! Initialize sdnew using sd
      initialize_sdnew_from_sd: block
         use parallel, only: group
         integer :: i,j,k,n
         type(coupler) :: sd2sdnew
         real(WP), dimension(:,:,:), allocatable :: tmp
         ! Create new coupler
         sd2sdnew=coupler(src_grp=group,dst_grp=group,name='sd2sd'); call sd2sdnew%set_src(sd%cfg); call sd2sdnew%set_dst(sdnew%cfg); call sd2sdnew%initialize()
         ! Allocate tmp array for transfer
         allocate(tmp(sdnew%fs%cfg%imino_:sdnew%fs%cfg%imaxo_,sdnew%fs%cfg%jmino_:sdnew%fs%cfg%jmaxo_,sdnew%fs%cfg%kmino_:sdnew%fs%cfg%kmaxo_))
         ! Transfer Q(1-7) - since the mesh is the same, we can safely ignore staggering here
         do n=1,7
            tmp=0.0_WP; call sd2sdnew%push(sd%fs%Q(:,:,:,n)); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp)
            do k=sdnew%cfg%kmino_,sdnew%cfg%kmaxo_; do j=sdnew%cfg%jmino_,sdnew%cfg%jmaxo_; do i=sdnew%cfg%imino_,sdnew%cfg%imaxo_
               if (sdnew%fs%cfg%xm(i).gt.sd%fs%cfg%x(sd%fs%cfg%imin).and.sdnew%fs%cfg%xm(i).lt.sd%fs%cfg%x(sd%fs%cfg%imax+1).and.&
               &   sdnew%fs%cfg%ym(j).gt.sd%fs%cfg%y(sd%fs%cfg%jmin).and.sdnew%fs%cfg%ym(j).lt.sd%fs%cfg%y(sd%fs%cfg%jmax+1).and.&
               &   sdnew%fs%cfg%zm(k).gt.sd%fs%cfg%z(sd%fs%cfg%kmin).and.sdnew%fs%cfg%zm(k).lt.sd%fs%cfg%z(sd%fs%cfg%kmax+1)) sdnew%fs%Q(i,j,k,n)=tmp(i,j,k)
            end do; end do; end do
         end do
         ! Transfer VOF
         tmp=0.0_WP; call sd2sdnew%push(sd%fs%VF); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp)
         do k=sdnew%cfg%kmino_,sdnew%cfg%kmaxo_; do j=sdnew%cfg%jmino_,sdnew%cfg%jmaxo_; do i=sdnew%cfg%imino_,sdnew%cfg%imaxo_
            if (sdnew%fs%cfg%xm(i).gt.sd%fs%cfg%x(sd%fs%cfg%imin).and.sdnew%fs%cfg%xm(i).lt.sd%fs%cfg%x(sd%fs%cfg%imax+1).and.&
            &   sdnew%fs%cfg%ym(j).gt.sd%fs%cfg%y(sd%fs%cfg%jmin).and.sdnew%fs%cfg%ym(j).lt.sd%fs%cfg%y(sd%fs%cfg%jmax+1).and.&
            &   sdnew%fs%cfg%zm(k).gt.sd%fs%cfg%z(sd%fs%cfg%kmin).and.sdnew%fs%cfg%zm(k).lt.sd%fs%cfg%z(sd%fs%cfg%kmax+1)) sdnew%fs%VF(i,j,k)=tmp(i,j,k)
         end do; end do; end do
         ! Transfer barycenters
         do n=1,3
            tmp=0.0_WP; call sd2sdnew%push(sd%fs%BG(n,:,:,:)); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp)
            do k=sdnew%cfg%kmino_,sdnew%cfg%kmaxo_; do j=sdnew%cfg%jmino_,sdnew%cfg%jmaxo_; do i=sdnew%cfg%imino_,sdnew%cfg%imaxo_
               if (sdnew%fs%cfg%xm(i).gt.sd%fs%cfg%x(sd%fs%cfg%imin).and.sdnew%fs%cfg%xm(i).lt.sd%fs%cfg%x(sd%fs%cfg%imax+1).and.&
               &   sdnew%fs%cfg%ym(j).gt.sd%fs%cfg%y(sd%fs%cfg%jmin).and.sdnew%fs%cfg%ym(j).lt.sd%fs%cfg%y(sd%fs%cfg%jmax+1).and.&
               &   sdnew%fs%cfg%zm(k).gt.sd%fs%cfg%z(sd%fs%cfg%kmin).and.sdnew%fs%cfg%zm(k).lt.sd%fs%cfg%z(sd%fs%cfg%kmax+1)) sdnew%fs%BG(n,i,j,k)=tmp(i,j,k)
            end do; end do; end do
         end do
         do n=1,3
            tmp=0.0_WP; call sd2sdnew%push(sd%fs%BL(n,:,:,:)); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp)
            do k=sdnew%cfg%kmino_,sdnew%cfg%kmaxo_; do j=sdnew%cfg%jmino_,sdnew%cfg%jmaxo_; do i=sdnew%cfg%imino_,sdnew%cfg%imaxo_
               if (sdnew%fs%cfg%xm(i).gt.sd%fs%cfg%x(sd%fs%cfg%imin).and.sdnew%fs%cfg%xm(i).lt.sd%fs%cfg%x(sd%fs%cfg%imax+1).and.&
               &   sdnew%fs%cfg%ym(j).gt.sd%fs%cfg%y(sd%fs%cfg%jmin).and.sdnew%fs%cfg%ym(j).lt.sd%fs%cfg%y(sd%fs%cfg%jmax+1).and.&
               &   sdnew%fs%cfg%zm(k).gt.sd%fs%cfg%z(sd%fs%cfg%kmin).and.sdnew%fs%cfg%zm(k).lt.sd%fs%cfg%z(sd%fs%cfg%kmax+1)) sdnew%fs%BL(n,i,j,k)=tmp(i,j,k)
            end do; end do; end do
         end do
         ! Build PLIC interface
         call sdnew%fs%build_interface()
         ! Communicate conserved variables just to be safe
         do n=1,sdnew%fs%nQ; call sdnew%fs%cfg%sync(sdnew%fs%Q(:,:,:,i)); end do
         ! Rebuild primitive variables
         call sdnew%fs%get_primitive()
         ! Interpolate velocity
         call sdnew%fs%interp_vel(sdnew%Ui,sdnew%Vi,sdnew%Wi)
         ! Compute local Mach number
         sdnew%Ma=sqrt(sdnew%Ui**2+sdnew%Vi**2+sdnew%Wi**2)/sdnew%fs%C
         ! Monitor
         call sdnew%output_monitor()
         ! Free memory
         deallocate(tmp)
      end block initialize_sdnew_from_sd
      
      ! Finally, transfer allocations
      transfer_allocation: block
         ! Copy couplers
         sdnew2ff=sd2ff
         ff2sdnew=ff2sd
         ! Copy sd
         sd=sdnew
      end block transfer_allocation
      
   end subroutine remesh
   
   !> Coupling from sd to ff
   subroutine couple_sd2ff()
      implicit none
      integer  :: i,j,k,n
      real(WP) :: coeff,lambda,strength
      real(WP), dimension(:,:,:,:), allocatable :: Q
      
      ! Sponge parameters
      lambda=2.5_WP*ff%fs%dx
      strength=0.25_WP
      
      ! Allocate storage for transfered variables
      allocate(Q(ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_,0:sd%fs%nQ)); Q=0.0_WP
      
      ! Exchange data using coupler
      call sd2ff%push(sd%fs%VF); call sd2ff%transfer(); call sd2ff%pull(Q(:,:,:,0))
      do n=1,4
         call sd2ff%push(sd%fs%Q(:,:,:,n)); call sd2ff%transfer(); call sd2ff%pull(Q(:,:,:,n))
      end do
      call sd2ff%push(sd%Ui); call sd2ff%transfer(); call sd2ff%pull(Q(:,:,:,5))
      call sd2ff%push(sd%Vi); call sd2ff%transfer(); call sd2ff%pull(Q(:,:,:,6))
      call sd2ff%push(sd%Wi); call sd2ff%transfer(); call sd2ff%pull(Q(:,:,:,7))
      
      ! Compute nudging increment
      do k=ff%cfg%kmino_,ff%cfg%kmaxo_; do j=ff%cfg%jmino_,ff%cfg%jmaxo_; do i=ff%cfg%imino_,ff%cfg%imaxo_
         ! Cell-centered nudging coefficient
         coeff=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[ff%cfg%xm(i),ff%cfg%ym(j),ff%cfg%zm(k)],delta=lambda)
         if (Q(i,j,k,0).gt.0.0_WP) coeff=0.0_WP
         ! Compute cell-centered momentum increment
         Q(i,j,k,5)=coeff*(Q(i,j,k,2)*Q(i,j,k,5)-ff%fs%Q(i,j,k,1)*ff%Ui(i,j,k))
         Q(i,j,k,6)=coeff*(Q(i,j,k,2)*Q(i,j,k,6)-ff%fs%Q(i,j,k,1)*ff%Vi(i,j,k))
         Q(i,j,k,7)=coeff*(Q(i,j,k,2)*Q(i,j,k,7)-ff%fs%Q(i,j,k,1)*ff%Wi(i,j,k))
         ! Compute RHO increment
         Q(i,j,k,1)=0.0_WP
         Q(i,j,k,2)=coeff*(Q(i,j,k,2)-ff%fs%Q(i,j,k,1))
         ! Compute RHO*I increment
         Q(i,j,k,3)=0.0_WP
         Q(i,j,k,4)=coeff*(Q(i,j,k,4)-ff%fs%Q(i,j,k,2))
      end do; end do; end do
      
      ! Second pass to apply forcing
      do k=ff%cfg%kmino_+1,ff%cfg%kmaxo_; do j=ff%cfg%jmino_+1,ff%cfg%jmaxo_; do i=ff%cfg%imino_+1,ff%cfg%imaxo_
         ff%fs%Q(i,j,k,1)=ff%fs%Q(i,j,k,1)+Q(i,j,k,2)
         ff%fs%Q(i,j,k,2)=ff%fs%Q(i,j,k,2)+Q(i,j,k,4)
         ff%fs%Q(i,j,k,3)=ff%fs%Q(i,j,k,3)+0.5_WP*sum(Q(i-1:i,j,k,5))
         ff%fs%Q(i,j,k,4)=ff%fs%Q(i,j,k,4)+0.5_WP*sum(Q(i,j-1:j,k,6))
         ff%fs%Q(i,j,k,5)=ff%fs%Q(i,j,k,5)+0.5_WP*sum(Q(i,j,k-1:k,7))
      end do; end do; end do
      
      ! Communicate conserved variables
      do n=1,ff%fs%nQ; call ff%cfg%sync(ff%fs%Q(:,:,:,n)); end do
      
      ! Recompute primitive variables
      call ff%fs%get_primitive()
      
      ! Free memory
      deallocate(Q)
      
   contains
      !> Function that calculates the signed distance between to domains
      real(WP) function sponge_forcing(inner,pos,delta)
         use pgrid_class, only: pgrid
         use mathtools,   only: Pi
         implicit none
         type(pgrid), intent(in) :: inner
         real(WP), dimension(3), intent(in) :: pos
         real(WP), intent(in) :: delta
         real(WP) :: dx,dy,dz,dx_in,dy_in,dz_in
         logical :: is_out_x,is_out_y,is_out_z
         ! X direction
         if (inner%nx.gt.1) then
            dx=max(inner%x(inner%imin)-pos(1),0.0_WP,pos(1)-inner%x(inner%imax+1))
            dx_in=min(inner%x(inner%imax+1)-pos(1),pos(1)-inner%x(inner%imin))
            is_out_x=(pos(1).lt.inner%x(inner%imin).or.pos(1).gt.inner%x(inner%imax+1))
         else
            dx=0.0_WP
            dx_in=huge(1.0_WP)
            is_out_x=.false.
         end if
         ! Y direction
         if (inner%ny.gt.1) then
            dy=max(inner%y(inner%jmin)-pos(2),0.0_WP,pos(2)-inner%y(inner%jmax+1))
            dy_in=min(inner%y(inner%jmax+1)-pos(2),pos(2)-inner%y(inner%jmin))
            is_out_y=(pos(2).lt.inner%y(inner%jmin).or.pos(2).gt.inner%y(inner%jmax+1))
         else
            dy=0.0_WP
            dy_in=huge(1.0_WP)
            is_out_y=.false.
         end if
         ! Z direction
         if (inner%nz.gt.1) then
            dz=max(inner%z(inner%kmin)-pos(3),0.0_WP,pos(3)-inner%z(inner%kmax+1))
            dz_in=min(inner%z(inner%kmax+1)-pos(3),pos(3)-inner%z(inner%kmin))
            is_out_z=(pos(3).lt.inner%z(inner%kmin).or.pos(3).gt.inner%z(inner%kmax+1))
         else
            dz=0.0_WP
            dz_in=huge(1.0_WP)
            is_out_z=.false.
         end if
         ! Signed distance
         if (is_out_x.or.is_out_y.or.is_out_z) then
            sponge_forcing=-sqrt(dx**2+dy**2+dz**2)
         else
            sponge_forcing=min(dx_in,dy_in,dz_in)
         end if
         ! Return a smooth forcing coefficient in [0,1]
         sponge_forcing=sponge_forcing/delta-1.0_WP
         if (sponge_forcing.le.0.0_WP) then
            sponge_forcing=0.0_WP
         else if (sponge_forcing.ge.1.0_WP) then
            sponge_forcing=1.0_WP
         else
            sponge_forcing=1.0_WP-cos(0.5_WP*Pi*sponge_forcing)**4
         end if
      end function sponge_forcing
   end subroutine couple_sd2ff
   
   
   !> Coupling from ff to sd
   subroutine couple_ff2sd()
      implicit none
      integer  :: i,j,k,n
      real(WP) :: coeff,lambda,strength
      real(WP), dimension(:,:,:,:), allocatable :: Q
      
      ! Sponge parameters
      lambda=5.0_WP*sd%fs%dx
      strength=0.1_WP
      
      ! Allocate storage for transfered variables
      allocate(Q(sd%fs%cfg%imino_:sd%fs%cfg%imaxo_,sd%fs%cfg%jmino_:sd%fs%cfg%jmaxo_,sd%fs%cfg%kmino_:sd%fs%cfg%kmaxo_,1:ff%fs%nQ)); Q=0.0_WP
      
      ! Exchange data using coupler
      do n=1,2
         call ff2sd%push(ff%fs%Q(:,:,:,n)); call ff2sd%transfer(); call ff2sd%pull(Q(:,:,:,n))
      end do
      call ff2sd%push(ff%Ui); call ff2sd%transfer(); call ff2sd%pull(Q(:,:,:,3))
      call ff2sd%push(ff%Vi); call ff2sd%transfer(); call ff2sd%pull(Q(:,:,:,4))
      call ff2sd%push(ff%Wi); call ff2sd%transfer(); call ff2sd%pull(Q(:,:,:,5))
      
      ! Compute nudging increment
      do k=sd%cfg%kmino_,sd%cfg%kmaxo_; do j=sd%cfg%jmino_,sd%cfg%jmaxo_; do i=sd%cfg%imino_,sd%cfg%imaxo_
         ! Cell-centered nudging coefficient
         coeff=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[sd%cfg%xm(i),sd%cfg%ym(j),sd%cfg%zm(k)],delta=lambda)
         ! Compute cell-centered momentum increment
         Q(i,j,k,3)=coeff*((Q(i,j,k,1)+sd%fs%Q(i,j,k,1))*Q(i,j,k,3)-sum(sd%fs%Q(i,j,k,1:2))*sd%Ui(i,j,k))
         Q(i,j,k,4)=coeff*((Q(i,j,k,1)+sd%fs%Q(i,j,k,1))*Q(i,j,k,4)-sum(sd%fs%Q(i,j,k,1:2))*sd%Vi(i,j,k))
         Q(i,j,k,5)=coeff*((Q(i,j,k,1)+sd%fs%Q(i,j,k,1))*Q(i,j,k,5)-sum(sd%fs%Q(i,j,k,1:2))*sd%Wi(i,j,k))
         ! Compute RHO increment
         Q(i,j,k,1)=coeff*(Q(i,j,k,1)-sd%fs%Q(i,j,k,2))
         ! Compute RHO*I increment
         Q(i,j,k,2)=coeff*(Q(i,j,k,2)-sd%fs%Q(i,j,k,4))
      end do; end do; end do
      
      ! Second pass to apply forcing
      do k=sd%cfg%kmino_+1,sd%cfg%kmaxo_; do j=sd%cfg%jmino_+1,sd%cfg%jmaxo_; do i=sd%cfg%imino_+1,sd%cfg%imaxo_
         sd%fs%Q(i,j,k,1)=sd%fs%Q(i,j,k,1)
         sd%fs%Q(i,j,k,2)=sd%fs%Q(i,j,k,2)+Q(i,j,k,1)
         sd%fs%Q(i,j,k,3)=sd%fs%Q(i,j,k,3)
         sd%fs%Q(i,j,k,4)=sd%fs%Q(i,j,k,4)+Q(i,j,k,2)
         sd%fs%Q(i,j,k,5)=sd%fs%Q(i,j,k,5)+0.5_WP*sum(Q(i-1:i,j,k,3))
         sd%fs%Q(i,j,k,6)=sd%fs%Q(i,j,k,6)+0.5_WP*sum(Q(i,j-1:j,k,4))
         sd%fs%Q(i,j,k,7)=sd%fs%Q(i,j,k,7)+0.5_WP*sum(Q(i,j,k-1:k,5))
      end do; end do; end do
      
      ! Communicate conserved variables
      do n=1,sd%fs%nQ; call sd%cfg%sync(sd%fs%Q(:,:,:,n)); end do
      
      ! Recompute primitive variables
      call sd%fs%get_primitive()
      
      ! Free memory
      deallocate(Q)
      
   contains
      !> Function that calculates the signed distance between to domains
      real(WP) function sponge_forcing(inner,pos,delta)
         use pgrid_class, only: pgrid
         use mathtools,   only: Pi
         implicit none
         type(pgrid), intent(in) :: inner
         real(WP), dimension(3), intent(in) :: pos
         real(WP), intent(in) :: delta
         real(WP) :: dx,dy,dz,dx_in,dy_in,dz_in
         logical :: is_out_x,is_out_y,is_out_z
         ! X direction
         if (inner%nx.gt.1) then ! Only force in -x
            dx=max(inner%x(inner%imin)-pos(1),0.0_WP)!,pos(1)-inner%x(inner%imax+1))
            dx_in=pos(1)-inner%x(inner%imin)!min(pos(1)-inner%x(inner%imin),inner%x(inner%imax+1)-pos(1))
            is_out_x=(pos(1).lt.inner%x(inner%imin))!.or.pos(1).gt.inner%x(inner%imax+1))
         else
            dx=0.0_WP
            dx_in=huge(1.0_WP)
            is_out_x=.false.
         end if
         ! Y direction
         if (inner%ny.gt.1) then
            dy=max(inner%y(inner%jmin)-pos(2),0.0_WP,pos(2)-inner%y(inner%jmax+1))
            dy_in=min(inner%y(inner%jmax+1)-pos(2),pos(2)-inner%y(inner%jmin))
            is_out_y=(pos(2).lt.inner%y(inner%jmin).or.pos(2).gt.inner%y(inner%jmax+1))
         else
            dy=0.0_WP
            dy_in=huge(1.0_WP)
            is_out_y=.false.
         end if
         ! Z direction
         if (inner%nz.gt.1) then
            dz=max(inner%z(inner%kmin)-pos(3),0.0_WP,pos(3)-inner%z(inner%kmax+1))
            dz_in=min(inner%z(inner%kmax+1)-pos(3),pos(3)-inner%z(inner%kmin))
            is_out_z=(pos(3).lt.inner%z(inner%kmin).or.pos(3).gt.inner%z(inner%kmax+1))
         else
            dz=0.0_WP
            dz_in=huge(1.0_WP)
            is_out_z=.false.
         end if
         ! Signed distance
         if (is_out_x.or.is_out_y.or.is_out_z) then
            sponge_forcing=-sqrt(dx**2+dy**2+dz**2)
         else
            sponge_forcing=min(dx_in,dy_in,dz_in)
         end if
         ! Return a smooth forcing coefficient in [0,1]
         sponge_forcing=sponge_forcing/delta
         if (sponge_forcing.ge.1.0_WP) then
            sponge_forcing=0.0_WP
         else if (sponge_forcing.ge.0.0_WP) then
            sponge_forcing=(0.5_WP*(1.0_WP+cos(Pi*sponge_forcing)))**4
         else
            sponge_forcing=1.0_WP
         end if
      end function sponge_forcing
   end subroutine couple_ff2sd
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      call sd%final()
      call ff%final()
   end subroutine simulation_final
   
   
end module simulation
