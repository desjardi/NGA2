!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,       only: WP
   use shockdrop_class, only: shockdrop
   use event_class,     only: event
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Shock-drop simulation
   type(shockdrop) :: sd
   
   !> Ensight output event
   type(event) :: ens_evt
   
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
         sd%viscL=viscL; sd%viscG=viscG
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
      
      ! Initialize Ensight output event and perform initial Ensight output
      initialize_ensight: block
         use param, only: param_read
         ens_evt=event(time=sd%time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         if (ens_evt%occurs()) then
            call sd%output_ensight(t=sd%time%t)
         end if
      end block initialize_ensight
      
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
         ! Advance shock-drop simulation
         call sd%step()
         ! Perform monitoring
         call sd%output_monitor()
         ! Perform Ensight output
         if (ens_evt%occurs()) then
            call sd%output_ensight(t=sd%time%t)
         end if
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize shock-drop simulation
      call sd%final()
      
   end subroutine simulation_final
   
   
end module simulation
