!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use timetracker_class, only: timetracker
   use shockdrop_class,   only: shockdrop
   use ffshock_class,     only: ffshock
   use coupler_class,     only: coupler
   use event_class,       only: event
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final
   
   !> Track time from here
   type(timetracker) :: time
   
   !> Shock-drop simulation - pointer since we will dynamically remesh
   type(shockdrop), pointer :: sd=>null()
   
   !> Far-field shock simulation - this is static
   type(ffshock) :: ff
   
   !> Couplers between domains - pointers since we will dynamically remesh
   type(coupler), pointer :: sd2ff=>null()
   type(coupler), pointer :: ff2sd=>null()
   
   !> Ensight output event
   type(event) :: ens_evt
   
   !> Remeshing event
   type(event) :: remesh_evt
   real(WP), dimension(3) :: Lmin,Lmax
   real(WP) :: Lmargin
   
   !> Maximum mesh size
   integer :: max_nx,max_ny,max_nz
   
   !> Equations of state
   real(WP) :: PinfL,GammaL,CvL
   real(WP) :: PinfG,GammaG,CvG
   
   !> Spherical harmonics perturbation parameters
   integer :: nsh_modes=0
   integer , dimension(:), allocatable :: l_modes,m_modes
   real(WP), dimension(:), allocatable :: amp_modes,phase_modes
   
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
   
   
   !> Mechanical relaxation model (implicit)
   subroutine P_relax_implicit(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: a,b,d,d1,d0,Peq,VFeq,invG1G,invG1L,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! Handle gas flotsams
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Setup quadratic problem
      invG1G=1.0_WP/(GammaG-1.0_WP); invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=PinfL*GammaL*invG1L; d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G; facL=invG1G+VF
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      d=d0*(facG-Q(4))-Q(3)*facG
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Adjust conserved quantities
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_implicit
   
   
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
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
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
      
      ! Initialize spherical harmonics modes - all cores
      initialize_sph_modes: block
         use param,     only: param_read,param_exists
         use parallel,  only: amRoot,MPI_REAL_WP,comm
         use string,    only: str_long
         use messager,  only: log
         use mpi_f08,   only: MPI_BCAST,MPI_INTEGER
         use random,    only: random_initialize,random_uniform
         use mathtools, only: twoPi
         character(str_long) :: message
         integer  :: i,ierr
         real(WP) :: r
         ! Read number of modes (or set default)
         call param_read('SphHarm Nmode',nsh_modes,default=0)
         if (nsh_modes.gt.0) then
            allocate(l_modes(nsh_modes),m_modes(nsh_modes),amp_modes(nsh_modes),phase_modes(nsh_modes))
            if (param_exists('SphHarm l'))     call param_read('SphHarm l',    l_modes)
            if (param_exists('SphHarm m'))     call param_read('SphHarm m',    m_modes)
            if (param_exists('SphHarm amp'))   call param_read('SphHarm amp',  amp_modes)
            if (param_exists('SphHarm phase')) call param_read('SphHarm phase',phase_modes)
         else
            ! Generate random modes if not provided
            nsh_modes=8
            allocate(l_modes(nsh_modes),m_modes(nsh_modes),amp_modes(nsh_modes),phase_modes(nsh_modes))
            if (amRoot) then
               call random_initialize()
               do i=1,nsh_modes
                  l_modes    (i)=1+i
                  m_modes    (i)=i-nsh_modes/2
                  amp_modes  (i)=1.0e-3_WP
                  phase_modes(i)=random_uniform(lo=0.0_WP,hi=twoPi)
               end do
            end if
            call MPI_BCAST(l_modes    ,nsh_modes,MPI_INTEGER,0,comm,ierr)
            call MPI_BCAST(m_modes    ,nsh_modes,MPI_INTEGER,0,comm,ierr)
            call MPI_BCAST(amp_modes  ,nsh_modes,MPI_REAL_WP,0,comm,ierr)
            call MPI_BCAST(phase_modes,nsh_modes,MPI_REAL_WP,0,comm,ierr)
         end if
         ! Log the selected modes
         if (amRoot) then
            write(message,'("[SphHarm] Nmode  =",i6)')         nsh_modes; call log(message)
            write(message,'("[SphHarm] l      =",1000(i4,x))') l_modes  ; call log(message)
            write(message,'("[SphHarm] m      =",1000(i4,x))') m_modes  ; call log(message)
            write(message,'("[SphHarm] amp    =",1000(es12.5,x))') amp_modes  ; call log(message)
            write(message,'("[SphHarm] phase  =",1000(es12.5,x))') phase_modes; call log(message)
         end if
      end block initialize_sph_modes
      
      ! Initialize time tracker
      initialize_timetracker: block
         use parallel, only: amRoot
         use param,    only: param_read
         ! Create time tracker object
         time=timetracker(amRoot=amRoot,name='Global')
         ! Set time integration parameters
         call param_read('Shock-drop dt',time%dtmax); time%dt=time%dtmax
         call param_read('Shock-drop CFL',time%cflmax)
         call param_read('Max time',time%tmax)
      end block initialize_timetracker
      
      ! Setup shock-drop simulation - all cores
      setup_sd: block
         use param,    only: param_read
         use parallel, only: group
         real(WP), dimension(3) :: X0
         integer , dimension(3) :: meshsize,partition
         real(WP) :: dx
         ! Read in mesh size and desired partition
         call param_read('Shock-drop dx',dx)
         call param_read('Shock-drop margin',Lmargin)
         call param_read('Shock-drop max nx',max_nx,default=0)
         call param_read('Shock-drop max ny',max_ny,default=0)
         call param_read('Shock-drop max nz',max_nz,default=0)
         call param_read('Shock-drop partition',partition)
         ! Set initial domain of size (2D)^3 centered on (0,0,0)
         meshsize=min(nint([(1.0_WP+2.0_WP*Lmargin)/dx,(1.0_WP+2.0_WP*Lmargin)/dx,(1.0_WP+2.0_WP*Lmargin)/dx]),merge([max_nx,max_ny,max_nz],huge(1),[max_nx,max_ny,max_nz].gt.0))
         X0=-0.5_WP*real(meshsize,WP)*dx
         ! Allocate and initialize the shock-drop solver
         allocate(sd); call sd%initialize(dx=dx,meshsize=meshsize,startloc=X0,group=group,partition=partition,continue_monitor=.false.)
         ! Provide relaxation and thermodynamic models
         sd%fs%relax=>P_relax_implicit
         sd%fs%getPL=>get_PL; sd%fs%getCL=>get_CL; sd%fs%getSL=>get_SL; sd%fs%getTL=>get_TL
         sd%fs%getPG=>get_PG; sd%fs%getCG=>get_CG; sd%fs%getSG=>get_SG; sd%fs%getTG=>get_TG
         ! We need to transfer our viscosities explicitly...
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
         call ff%initialize(dx=dx,meshsize=meshsize,startloc=X0,group=group,partition=partition)
         ! Provide thermodynamic model
         ff%fs%getP=>get_PG; ff%fs%getC=>get_CG; ff%fs%getS=>get_SG; ff%fs%getT=>get_TG
         ! We need to transfer our viscosity explicitly...
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
         ! Initialize lpt thermodynamic variables
         ff%lp%Cp=GammaL*CvL !< Incorrect for stiffened gas...
         ff%lp%rho=rhoL      !< Should be variable...
         ! Perform monitoring
         call ff%output_monitor()
      end block initialize_ff
      
      ! Create couplers
      create_couplers: block
         use parallel, only: group
         allocate(sd2ff); sd2ff=coupler(src_grp=group,dst_grp=group,name='sd2ff'); call sd2ff%set_src(sd%cfg); call sd2ff%set_dst(ff%cfg); call sd2ff%initialize()
         allocate(ff2sd); ff2sd=coupler(src_grp=group,dst_grp=group,name='ff2sd'); call ff2sd%set_src(ff%cfg); call ff2sd%set_dst(sd%cfg); call ff2sd%initialize()
      end block create_couplers
      
      ! Initialize Ensight output event and perform initial Ensight output
      initialize_ensight: block
         use param, only: param_read
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         if (ens_evt%occurs()) then
            call sd%output_ensight(t=time%t)
            call ff%output_ensight(t=time%t)
         end if
      end block initialize_ensight
      
      ! Initialize remeshing event
      initialize_remeshing: block
         use param, only: param_read
         remesh_evt=event(time=time,name='Remeshing')
         call param_read('Remeshing period',remesh_evt%tper)
         call sd%meshfile%write()
      end block initialize_remeshing
      
   contains
      !> Level set function for a perturbed sphere of unity diameter centered at (0,0,0)
      function levelset_drop(xyz,t) result(G)
         use mathtools, only: spherical_harmonic
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G,r,theta,phi,perturb
         integer :: i
         ! Spherical coordinates
         r=sqrt(sum(xyz**2))
         if (r.gt.1.0e-12_WP) then
            theta= acos(xyz(3)/r)
            phi  =atan2(xyz(2),xyz(1))
         else
            theta=0.0_WP
            phi  =0.0_WP
         end if
         ! Compute perturbation
         perturb=0.0_WP
         do i=1,nsh_modes
            perturb=perturb+amp_modes(i)*spherical_harmonic(l_modes(i),m_modes(i),theta,phi+phase_modes(i))
         end do
         ! Level set function for a sphere with radius 0.5 and perturbation
         G=0.5_WP+perturb-r
      end function levelset_drop
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Overall time integration
      do while (.not.time%done())
         
         ! Adjust time step size using sd CFL info
         call sd%fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Handle coupling
         call couple_sd2ff()
         call couple_ff2sd()
         
         ! Advance shock-drop simulation
         call sd%step(dt=time%dt)
         
         ! Advance farfield simulation
         call ff%step(dt=time%dt)
         
         ! Transfer drops
         call transfer()
         
         ! Perform monitoring
         call sd%output_monitor()
         call ff%output_monitor()
         
         ! Perform Ensight output
         if (ens_evt%occurs()) then
            call sd%output_ensight(t=time%t)
            call ff%output_ensight(t=time%t)
         end if
         
         ! Remesh sd
         if (remesh_evt%occurs()) call remesh()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Remesh sd to follow the drop
   subroutine remesh()
      use, intrinsic :: iso_fortran_env, only: output_unit
      use messager, only: log
      use parallel, only: amRoot
      use string,   only: str_long
      implicit none
      character(len=str_long) :: message
      type(shockdrop), pointer :: sdnew
      type(coupler)  , pointer :: sdnew2ff
      type(coupler)  , pointer :: ff2sdnew
      
      ! Setup new shock-drop simulation - all cores
      setup_sdnew: block
         use parallel, only: group,MPI_REAL_WP
         use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
         real(WP), dimension(3) :: X0,X1
         integer , dimension(3) :: meshsize,partition
         real(WP) :: dx
         integer  :: ierr
         ! Use same dx and partition as sd
         dx=sd%fs%dx
         partition=[sd%cfg%npx,sd%cfg%npy,sd%cfg%npz]
         ! Get meshsize that encompasses current liquid extent
         X0(1)=dx*floor  ((Lmin(1)-Lmargin)/dx); if (sd%cfg%nx.eq.1) X0(1)=Lmin(1)
         X0(2)=dx*floor  ((Lmin(2)-Lmargin)/dx); if (sd%cfg%ny.eq.1) X0(2)=Lmin(2)
         X0(3)=dx*floor  ((Lmin(3)-Lmargin)/dx); if (sd%cfg%nz.eq.1) X0(3)=Lmin(3)
         X1(1)=dx*ceiling((Lmax(1)+Lmargin)/dx); if (sd%cfg%nx.eq.1) X1(1)=Lmax(1)
         X1(2)=dx*ceiling((Lmax(2)+Lmargin)/dx); if (sd%cfg%ny.eq.1) X1(2)=Lmax(2)
         X1(3)=dx*ceiling((Lmax(3)+Lmargin)/dx); if (sd%cfg%nz.eq.1) X1(3)=Lmax(3)
         meshsize=nint((X1-X0)/dx)
         ! Ensure consistency across ranks
         call MPI_BCAST(meshsize,3,MPI_INTEGER,0,sd%cfg%comm,ierr)
         call MPI_BCAST(X0      ,3,MPI_REAL_WP,0,sd%cfg%comm,ierr)
         ! Initialize the shock-drop solver
         allocate(sdnew); call sdnew%initialize(dx=dx,meshsize=meshsize,startloc=X0,group=group,partition=partition,continue_monitor=.true.)
         ! Provide relaxation and thermodynamic models
         sdnew%fs%relax=>sd%fs%relax
         sdnew%fs%getPL=>sd%fs%getPL; sdnew%fs%getCL=>sd%fs%getCL; sdnew%fs%getSL=>sd%fs%getSL; sdnew%fs%getTL=>sd%fs%getTL
         sdnew%fs%getPG=>sd%fs%getPG; sdnew%fs%getCG=>sd%fs%getCG; sdnew%fs%getSG=>sd%fs%getSG; sdnew%fs%getTG=>sd%fs%getTG
         ! We need to transfer our viscosities explicitly...
         sdnew%cst_viscL=sd%cst_viscL; sdnew%cst_viscG=sd%cst_viscG
         ! Inform sdnew's timetracker of our current time, but leave n unchanged to make remeshing obvious
         sdnew%time%t=time%t
      end block setup_sdnew
      
      ! Create new couplers
      setup_new_couplers: block
         use parallel, only: group
         allocate(sdnew2ff); sdnew2ff=coupler(src_grp=group,dst_grp=group,name='sd2ff'); call sdnew2ff%set_src(sdnew%cfg); call sdnew2ff%set_dst(ff%cfg); call sdnew2ff%initialize()
         allocate(ff2sdnew); ff2sdnew=coupler(src_grp=group,dst_grp=group,name='ff2sd'); call ff2sdnew%set_src(ff%cfg); call ff2sdnew%set_dst(sdnew%cfg); call ff2sdnew%initialize()
      end block setup_new_couplers
      
      ! Initialize all sdnew to gas including in ghost cells
      initialize_to_gas: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         integer :: i,j,k
         ! Let us set PLIC explicitly to gas everywhere and provide corresponding volume moments
         do k=sdnew%cfg%kmino_,sdnew%cfg%kmaxo_; do j=sdnew%cfg%jmino_,sdnew%cfg%jmaxo_; do i=sdnew%cfg%imino_,sdnew%cfg%imaxo_
            sdnew%fs%VF(i,j,k)=0.0_WP; sdnew%fs%BL(:,i,j,k)=[sdnew%fs%cfg%xm(i),sdnew%fs%cfg%ym(j),sdnew%fs%cfg%zm(k)]; sdnew%fs%BG(:,i,j,k)=[sdnew%fs%cfg%xm(i),sdnew%fs%cfg%ym(j),sdnew%fs%cfg%zm(k)]
            call setNumberOfPlanes(sdnew%fs%PLIC(i,j,k),1); call setPlane(sdnew%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,sdnew%fs%VF(i,j,k)-0.5_WP))
         end do; end do; end do
      end block initialize_to_gas
      
      ! Initialize sdnew using ff
      initialize_sdnew_from_ff: block
         integer :: n
         ! Transfer data
         call ff2sdnew%push(ff%fs%Q(:,:,:,1),loc='c'); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,2),loc='c')
         call ff2sdnew%push(ff%fs%Q(:,:,:,2),loc='c'); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,4),loc='c')
         call ff2sdnew%push(ff%fs%Q(:,:,:,3),loc='x'); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,5),loc='x')
         call ff2sdnew%push(ff%fs%Q(:,:,:,4),loc='y'); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,6),loc='y')
         call ff2sdnew%push(ff%fs%Q(:,:,:,5),loc='z'); call ff2sdnew%transfer(); call ff2sdnew%pull(sdnew%fs%Q(:,:,:,7),loc='z')
         ! Communicate conserved variables
         do n=1,sdnew%fs%nQ; call sdnew%cfg%sync(sdnew%fs%Q(:,:,:,n)); end do
         ! Rebuild primitive variables
         call sdnew%fs%get_primitive()
      end block initialize_sdnew_from_ff
      
      ! Initialize sdnew using sd
      initialize_sdnew_from_sd: block
         use parallel, only: group
         integer :: n
         type(coupler) :: sd2sdnew
         real(WP), dimension(:,:,:), allocatable :: tmp1,tmp2
         ! Create new coupler
         sd2sdnew=coupler(src_grp=group,dst_grp=group,name='sd2sd'); call sd2sdnew%set_src(sd%cfg); call sd2sdnew%set_dst(sdnew%cfg); call sd2sdnew%initialize()
         ! Transfer Q(1-7)
         call sd2sdnew%push(sd%fs%Q(:,:,:,1),loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,1),loc='c')
         call sd2sdnew%push(sd%fs%Q(:,:,:,2),loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,2),loc='c')
         call sd2sdnew%push(sd%fs%Q(:,:,:,3),loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,3),loc='c')
         call sd2sdnew%push(sd%fs%Q(:,:,:,4),loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,4),loc='c')
         call sd2sdnew%push(sd%fs%Q(:,:,:,5),loc='x'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,5),loc='x')
         call sd2sdnew%push(sd%fs%Q(:,:,:,6),loc='y'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,6),loc='y')
         call sd2sdnew%push(sd%fs%Q(:,:,:,7),loc='z'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%Q(:,:,:,7),loc='z')
         ! Transfer VF
         call sd2sdnew%push(sd%fs%VF,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(sdnew%fs%VF,loc='c')
         ! Transfer barycenters
         allocate(tmp1(   sd%fs%cfg%imino_:   sd%fs%cfg%imaxo_,   sd%fs%cfg%jmino_:   sd%fs%cfg%jmaxo_,   sd%fs%cfg%kmino_:   sd%fs%cfg%kmaxo_))
         allocate(tmp2(sdnew%fs%cfg%imino_:sdnew%fs%cfg%imaxo_,sdnew%fs%cfg%jmino_:sdnew%fs%cfg%jmaxo_,sdnew%fs%cfg%kmino_:sdnew%fs%cfg%kmaxo_))
         tmp1=sd%fs%BG(1,:,:,:); call sd2sdnew%push(tmp1,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp2,loc='c'); sdnew%fs%BG(1,:,:,:)=tmp2
         tmp1=sd%fs%BG(2,:,:,:); call sd2sdnew%push(tmp1,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp2,loc='c'); sdnew%fs%BG(2,:,:,:)=tmp2
         tmp1=sd%fs%BG(3,:,:,:); call sd2sdnew%push(tmp1,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp2,loc='c'); sdnew%fs%BG(3,:,:,:)=tmp2
         tmp1=sd%fs%BL(1,:,:,:); call sd2sdnew%push(tmp1,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp2,loc='c'); sdnew%fs%BL(1,:,:,:)=tmp2
         tmp1=sd%fs%BL(2,:,:,:); call sd2sdnew%push(tmp1,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp2,loc='c'); sdnew%fs%BL(2,:,:,:)=tmp2
         tmp1=sd%fs%BL(3,:,:,:); call sd2sdnew%push(tmp1,loc='c'); call sd2sdnew%transfer(); call sd2sdnew%pull(tmp2,loc='c'); sdnew%fs%BL(3,:,:,:)=tmp2
         deallocate(tmp1,tmp2)
         ! Destroy coupler
         call sd2sdnew%finalize()
         ! Communicate conserved variables
         do n=1,sdnew%fs%nQ; call sdnew%fs%cfg%sync(sdnew%fs%Q(:,:,:,n)); end do
         ! Also sync volume moments
         call sdnew%fs%sync_volume_moments()
         ! Build PLIC interface
         call sdnew%fs%build_interface()
         ! Rebuild primitive variables
         call sdnew%fs%get_primitive()
         ! Re-apply boundary conditions
         call sdnew%apply_bconds()
         ! Interpolate velocity
         call sdnew%fs%interp_vel(sdnew%Ui,sdnew%Vi,sdnew%Wi)
         ! Compute local Mach number
         sdnew%Ma=sqrt(sdnew%Ui**2+sdnew%Vi**2+sdnew%Wi**2)/sdnew%fs%C
      end block initialize_sdnew_from_sd
      
      ! Finally, transfer allocation
      transfer_allocation: block
         ! Finalize and free up couplers, point to new ones
         call sd2ff%finalize(); deallocate(sd2ff); sd2ff=>sdnew2ff
         call ff2sd%finalize(); deallocate(ff2sd); ff2sd=>ff2sdnew
         ! Finalize and free up sd, point to new one
         call sd%finalize(); deallocate(sd); sd=>sdnew
      end block transfer_allocation
      
      ! Monitor mesh size
      call sd%meshfile%write()
      
      ! Some messaging
      if (amRoot) then
         write(message,'(">> Remeshed shockdrop domain: ",i0,"x",i0,"x",i0," from X0=[",es12.5,"x",es12.5,"x",es12.5,"]")') sd%cfg%nx,sd%cfg%ny,sd%cfg%nz,sd%cfg%x(sd%cfg%imin),sd%cfg%y(sd%cfg%jmin),sd%cfg%z(sd%cfg%kmin)
         call log(trim(message))
         write(output_unit,*) trim(message)
      end if
      
   end subroutine remesh
   
   
   !> Coupling from sd to ff
   subroutine couple_sd2ff()
      implicit none
      integer  :: i,j,k,n
      real(WP) :: coeffc,coeffx,coeffy,coeffz,lambda,strength
      real(WP), dimension(:,:,:), allocatable :: VOFtmp,RHOtmp,Itmp,Utmp,Vtmp,Wtmp
      
      ! Sponge parameters
      lambda=2.5_WP*ff%fs%dx
      strength=0.25_WP
      
      ! Allocate storage for transfered variables
      allocate(VOFtmp(ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_)); VOFtmp=0.0_WP
      allocate(RHOtmp(ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_)); RHOtmp=0.0_WP
      allocate(Itmp  (ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_)); Itmp  =0.0_WP
      allocate(Utmp  (ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_)); Utmp  =0.0_WP
      allocate(Vtmp  (ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_)); Vtmp  =0.0_WP
      allocate(Wtmp  (ff%fs%cfg%imino_:ff%fs%cfg%imaxo_,ff%fs%cfg%jmino_:ff%fs%cfg%jmaxo_,ff%fs%cfg%kmino_:ff%fs%cfg%kmaxo_)); Wtmp  =0.0_WP
      
      ! Nudge ff using sd data
      call sd2ff%push(sd%fs%VF  ,loc='c'); call sd2ff%transfer(); call sd2ff%pull(VOFtmp,loc='c')
      call sd2ff%push(sd%fs%RHOG,loc='c'); call sd2ff%transfer(); call sd2ff%pull(RHOtmp,loc='c')
      call sd2ff%push(sd%fs%IG  ,loc='c'); call sd2ff%transfer(); call sd2ff%pull(Itmp  ,loc='c')
      call sd2ff%push(sd%fs%U   ,loc='x'); call sd2ff%transfer(); call sd2ff%pull(Utmp  ,loc='x')
      call sd2ff%push(sd%fs%V   ,loc='y'); call sd2ff%transfer(); call sd2ff%pull(Vtmp  ,loc='y')
      call sd2ff%push(sd%fs%W   ,loc='z'); call sd2ff%transfer(); call sd2ff%pull(Wtmp  ,loc='z')
      do k=ff%cfg%kmino_,ff%cfg%kmaxo_; do j=ff%cfg%jmino_,ff%cfg%jmaxo_; do i=ff%cfg%imino_,ff%cfg%imaxo_
         ! Compute cell-centered nudging coefficient
         coeffc=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[ff%cfg%xm(i),ff%cfg%ym(j),ff%cfg%zm(k)],delta=lambda)
         ! Only nudge in pure gas regions
         if (VOFtmp(i,j,k).gt.0.0_WP) coeffc=0.0_WP
         ! Apply nudging on density and internal energy
         ff%fs%Q(i,j,k,1)=ff%fs%Q(i,j,k,1)+coeffc*(RHOtmp(i,j,k)            -ff%fs%Q(i,j,k,1))
         ff%fs%Q(i,j,k,2)=ff%fs%Q(i,j,k,2)+coeffc*(RHOtmp(i,j,k)*Itmp(i,j,k)-ff%fs%Q(i,j,k,2))
         ! Compute face-centered nudging coefficients
         coeffx=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[ff%cfg%x (i),ff%cfg%ym(j),ff%cfg%zm(k)],delta=lambda)
         coeffy=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[ff%cfg%xm(i),ff%cfg%y (j),ff%cfg%zm(k)],delta=lambda)
         coeffz=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[ff%cfg%xm(i),ff%cfg%ym(j),ff%cfg%z (k)],delta=lambda)
         ! Apply nudging on velocity
         ff%fs%U(i,j,k)=ff%fs%U(i,j,k)+coeffx*(Utmp(i,j,k)-ff%fs%U(i,j,k))
         ff%fs%V(i,j,k)=ff%fs%V(i,j,k)+coeffy*(Vtmp(i,j,k)-ff%fs%V(i,j,k))
         ff%fs%W(i,j,k)=ff%fs%W(i,j,k)+coeffz*(Wtmp(i,j,k)-ff%fs%W(i,j,k))
      end do; end do; end do
      
      ! Recompute momentum
      call ff%fs%get_momentum()
      
      ! Recompute primitive variables
      call ff%fs%get_primitive()
      
      ! Free up memory
      deallocate(VOFtmp,RHOtmp,Itmp,Utmp,Vtmp,Wtmp)
      
   contains
      !> Function that calculates the signed distance between two domains
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
      real(WP) :: coeffc,coeffx,coeffy,coeffz,lambda,strength
      real(WP), dimension(:,:,:), allocatable :: RHOtmp,RHOItmp,Utmp,Vtmp,Wtmp
      
      ! Sponge parameters
      lambda=5.0_WP*sd%fs%dx
      strength=0.125_WP
      
      ! Allocate storage for transfered variables
      allocate(RHOtmp (sd%fs%cfg%imino_:sd%fs%cfg%imaxo_,sd%fs%cfg%jmino_:sd%fs%cfg%jmaxo_,sd%fs%cfg%kmino_:sd%fs%cfg%kmaxo_)); RHOtmp =0.0_WP
      allocate(RHOItmp(sd%fs%cfg%imino_:sd%fs%cfg%imaxo_,sd%fs%cfg%jmino_:sd%fs%cfg%jmaxo_,sd%fs%cfg%kmino_:sd%fs%cfg%kmaxo_)); RHOItmp=0.0_WP
      allocate(Utmp   (sd%fs%cfg%imino_:sd%fs%cfg%imaxo_,sd%fs%cfg%jmino_:sd%fs%cfg%jmaxo_,sd%fs%cfg%kmino_:sd%fs%cfg%kmaxo_)); Utmp   =0.0_WP
      allocate(Vtmp   (sd%fs%cfg%imino_:sd%fs%cfg%imaxo_,sd%fs%cfg%jmino_:sd%fs%cfg%jmaxo_,sd%fs%cfg%kmino_:sd%fs%cfg%kmaxo_)); Vtmp   =0.0_WP
      allocate(Wtmp   (sd%fs%cfg%imino_:sd%fs%cfg%imaxo_,sd%fs%cfg%jmino_:sd%fs%cfg%jmaxo_,sd%fs%cfg%kmino_:sd%fs%cfg%kmaxo_)); Wtmp   =0.0_WP
      
      ! Nudge sd using ff data
      call ff2sd%push(ff%fs%Q(:,:,:,1),loc='c'); call ff2sd%transfer(); call ff2sd%pull(RHOtmp ,loc='c')
      call ff2sd%push(ff%fs%Q(:,:,:,2),loc='c'); call ff2sd%transfer(); call ff2sd%pull(RHOItmp,loc='c')
      call ff2sd%push(ff%fs%U         ,loc='x'); call ff2sd%transfer(); call ff2sd%pull(Utmp   ,loc='x')
      call ff2sd%push(ff%fs%V         ,loc='y'); call ff2sd%transfer(); call ff2sd%pull(Vtmp   ,loc='y')
      call ff2sd%push(ff%fs%W         ,loc='z'); call ff2sd%transfer(); call ff2sd%pull(Wtmp   ,loc='z')
      do k=sd%cfg%kmino_,sd%cfg%kmaxo_; do j=sd%cfg%jmino_,sd%cfg%jmaxo_; do i=sd%cfg%imino_,sd%cfg%imaxo_
         ! Compute cell-centered nudging coefficient
         coeffc=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[sd%cfg%xm(i),sd%cfg%ym(j),sd%cfg%zm(k)],delta=lambda)
         ! Apply nudging on gas density and gas internal energy
         sd%fs%Q(i,j,k,2)=sd%fs%Q(i,j,k,2)+coeffc*((1.0_WP-sd%fs%VF(i,j,k))*RHOtmp (i,j,k)-sd%fs%Q(i,j,k,2))
         sd%fs%Q(i,j,k,4)=sd%fs%Q(i,j,k,4)+coeffc*((1.0_WP-sd%fs%VF(i,j,k))*RHOItmp(i,j,k)-sd%fs%Q(i,j,k,4))
         ! Compute face-centered nudging coefficients
         coeffx=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[sd%cfg%x (i),sd%cfg%ym(j),sd%cfg%zm(k)],delta=lambda)
         coeffy=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[sd%cfg%xm(i),sd%cfg%y (j),sd%cfg%zm(k)],delta=lambda)
         coeffz=strength*sponge_forcing(inner=sd%cfg%pgrid,pos=[sd%cfg%xm(i),sd%cfg%ym(j),sd%cfg%z (k)],delta=lambda)
         ! Apply nudging on velocity
         sd%fs%U(i,j,k)=sd%fs%U(i,j,k)+coeffx*(Utmp(i,j,k)-sd%fs%U(i,j,k))
         sd%fs%V(i,j,k)=sd%fs%V(i,j,k)+coeffy*(Vtmp(i,j,k)-sd%fs%V(i,j,k))
         sd%fs%W(i,j,k)=sd%fs%W(i,j,k)+coeffz*(Wtmp(i,j,k)-sd%fs%W(i,j,k))
      end do; end do; end do
      
      ! Recompute momentum
      call sd%fs%get_momentum()
      
      ! Recompute primitive variables
      call sd%fs%get_primitive()
      
      ! Free up memory
      deallocate(RHOtmp,RHOItmp,Utmp,Vtmp,Wtmp)
      
   contains
      !> Function that calculates the signed distance between two domains
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
            dx_in=min(pos(1)-inner%x(inner%imin),inner%x(inner%imax+1)-pos(1))
            is_out_x=(pos(1).lt.inner%x(inner%imin).or.pos(1).gt.inner%x(inner%imax+1))
            ! Only force in -x
            !dx=max(inner%x(inner%imin)-pos(1),0.0_WP)!,pos(1)-inner%x(inner%imax+1))
            !dx_in=pos(1)-inner%x(inner%imin)!min(pos(1)-inner%x(inner%imin),inner%x(inner%imax+1)-pos(1))
            !is_out_x=(pos(1).lt.inner%x(inner%imin))!.or.pos(1).gt.inner%x(inner%imax+1))
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
   
   
   !> Transfer drops
   subroutine transfer()
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MIN,MPI_MAX
      use parallel,  only: MPI_REAL_WP,amRoot
      use mathtools, only: Pi
      use messager,  only: log
      use string,    only: str_long
      use irl_fortran_interface, only: setPlane,zeroPolygon
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      integer :: i,j,k,n,m,ierr,nremoved,ncreated,np
      real(WP), dimension(:)  , allocatable :: Vd,Md,Pd
      real(WP), dimension(:,:), allocatable :: Bd,Ud
      real(WP), parameter :: vol_coeff=10.0_WP
      real(WP), parameter :: diameter_threshold=1.0e-2_WP
      real(WP), dimension(3) :: edgelo,edgehi
      character(len=str_long) :: message
      
      ! Update sd's CCL
      call sd%ccl%build(make_label,same_label)
      
      ! Allocate volume, mass, pressure, barycenter, and velocity arrays
      allocate(Vd(    1:sd%ccl%nstruct)); Vd=0.0_WP
      allocate(Md(    1:sd%ccl%nstruct)); Md=0.0_WP
      allocate(Pd(    1:sd%ccl%nstruct)); Pd=0.0_WP
      allocate(Bd(1:3,1:sd%ccl%nstruct)); Bd=0.0_WP
      allocate(Ud(1:3,1:sd%ccl%nstruct)); Ud=0.0_WP
      
      ! Loop over individual structures and compute structure properties
      do n=1,sd%ccl%nstruct
         ! Loop over cells in structure and accumulate data
         do m=1,sd%ccl%struct(n)%n_
            ! Get cell index
            i=sd%ccl%struct(n)%map(1,m); j=sd%ccl%struct(n)%map(2,m); k=sd%ccl%struct(n)%map(3,m)
            ! Increment volume, mass, pressure, barycenter, and momentum
            Vd  (n)=Vd  (n)+sd%fs%cfg%vol(i,j,k)*sd%fs%VF(i,j,k)
            Md  (n)=Md  (n)+sd%fs%cfg%vol(i,j,k)*sd%fs%Q (i,j,k,1)
            Pd  (n)=Pd  (n)+sd%fs%cfg%vol(i,j,k)*sd%fs%VF(i,j,k)*sd%fs%PL(i,j,k)
            Bd(:,n)=Bd(:,n)+sd%fs%cfg%vol(i,j,k)*sd%fs%Q (i,j,k,1)*sd%fs%BL(:,i,j,k)
            Ud(1,n)=Ud(1,n)+sd%fs%cfg%vol(i,j,k)*sd%fs%VF(i,j,k)*sd%Ui(i,j,k)
            Ud(2,n)=Ud(2,n)+sd%fs%cfg%vol(i,j,k)*sd%fs%VF(i,j,k)*sd%Vi(i,j,k)
            Ud(3,n)=Ud(3,n)+sd%fs%cfg%vol(i,j,k)*sd%fs%VF(i,j,k)*sd%Wi(i,j,k)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,Vd,  sd%ccl%nstruct,MPI_REAL_WP,MPI_SUM,sd%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Md,  sd%ccl%nstruct,MPI_REAL_WP,MPI_SUM,sd%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Pd,  sd%ccl%nstruct,MPI_REAL_WP,MPI_SUM,sd%fs%cfg%comm,ierr); Pd=Pd/Vd
      call MPI_ALLREDUCE(MPI_IN_PLACE,Bd,3*sd%ccl%nstruct,MPI_REAL_WP,MPI_SUM,sd%fs%cfg%comm,ierr); do n=1,sd%ccl%nstruct; Bd(:,n)=Bd(:,n)/Md(n); end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,Ud,3*sd%ccl%nstruct,MPI_REAL_WP,MPI_SUM,sd%fs%cfg%comm,ierr); do n=1,sd%ccl%nstruct; Ud(:,n)=Ud(:,n)/Vd(n); end do
      
      ! Loop again to calculate extent of large enough drops - will be used for remeshing
      Lmin=[+huge(1.0_WP),+huge(1.0_WP),+huge(1.0_WP)]; Lmax=[-huge(1.0_WP),-huge(1.0_WP),-huge(1.0_WP)]
      do n=1,sd%ccl%nstruct
         ! Ensure structure is large enough
         if (Vd(n).le.vol_coeff*sd%fs%vol) cycle
         ! Loop over cells in structure and increment extent
         do m=1,sd%ccl%struct(n)%n_
            ! Get cell index
            i=sd%ccl%struct(n)%map(1,m); j=sd%ccl%struct(n)%map(2,m); k=sd%ccl%struct(n)%map(3,m)
            ! Update liquid extent
            Lmin=min(Lmin,[sd%cfg%x(i  ),sd%cfg%y(j  ),sd%cfg%z(k  )])
            Lmax=max(Lmax,[sd%cfg%x(i+1),sd%cfg%y(j+1),sd%cfg%z(k+1)])
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,Lmin,3,MPI_REAL_WP,MPI_MIN,sd%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Lmax,3,MPI_REAL_WP,MPI_MAX,sd%fs%cfg%comm,ierr)
      
      ! Transfer is based on closeness to edge of sd domain
      edgelo=[sd%fs%cfg%x(sd%fs%cfg%imin  ),sd%fs%cfg%y(sd%fs%cfg%jmin  ),sd%fs%cfg%z(sd%fs%cfg%kmin  )]
      edgehi=[sd%fs%cfg%x(sd%fs%cfg%imax+1),sd%fs%cfg%y(sd%fs%cfg%jmax+1),sd%fs%cfg%z(sd%fs%cfg%kmax+1)]
      ! If we're about to remesh sd, use extent of new sd domain
      if (remesh_evt%occurs()) then
         edgelo(1)=sd%fs%dx*floor  ((Lmin(1)-Lmargin)/sd%fs%dx); if (sd%cfg%nx.eq.1) edgelo(1)=Lmin(1)
         edgelo(2)=sd%fs%dx*floor  ((Lmin(2)-Lmargin)/sd%fs%dx); if (sd%cfg%ny.eq.1) edgelo(2)=Lmin(2)
         edgelo(3)=sd%fs%dx*floor  ((Lmin(3)-Lmargin)/sd%fs%dx); if (sd%cfg%nz.eq.1) edgelo(3)=Lmin(3)
         edgehi(1)=sd%fs%dx*ceiling((Lmax(1)+Lmargin)/sd%fs%dx); if (sd%cfg%nx.eq.1) edgehi(1)=Lmax(1)
         edgehi(2)=sd%fs%dx*ceiling((Lmax(2)+Lmargin)/sd%fs%dx); if (sd%cfg%ny.eq.1) edgehi(2)=Lmax(2)
         edgehi(3)=sd%fs%dx*ceiling((Lmax(3)+Lmargin)/sd%fs%dx); if (sd%cfg%nz.eq.1) edgehi(3)=Lmax(3)
      end if
      ! Add a buffer of 10 cells
      edgelo=edgelo+10.0_WP*sd%fs%dx
      edgehi=edgehi-10.0_WP*sd%fs%dx
      
      ! Now traverse all structures again to perform transfer
      nremoved=0; ncreated=0
      do n=1,sd%ccl%nstruct
         ! Check if volume is small enough for transfer
         if (Vd(n).gt.vol_coeff*sd%fs%vol) cycle
         ! Check if drop is close to the edge of the sd domain
         if (close_to_edge(Bd(:,n))) then
            ! Increment counter
            nremoved=nremoved+1
            ! Perform transfer
            if (ff%cfg%amRoot) then
               ! Make room for new drop
               np=ff%lp%np_+1; call ff%lp%resize(np)
               ! Give the drop an id
               ff%lp%p(np)%id  =int(1,8)
               ! Assign equivalent diameter from the volume
               ff%lp%p(np)%d   =(6.0_WP*Vd(n)/Pi)**(1.0_WP/3.0_WP)
               ! Assign droplet temperature from mass and EOS
               ff%lp%p(np)%T   =get_TL(RHO=Md(n)/Vd(n),P=Pd(n))
               ! Place the drop at the liquid barycenter
               ff%lp%p(np)%pos =Bd(:,n)
               ! Assign mean structure velocity as drop velocity
               ff%lp%p(np)%vel =Ud(:,n)
               ! Give zero angular velocity
               ff%lp%p(np)%angVel=0.0_WP
               ! Give zero collision force
               ff%lp%p(np)%Acol=0.0_WP
               ff%lp%p(np)%Tcol=0.0_WP
               ! Let the drop find it own integration time
               ff%lp%p(np)%dt  =0.0_WP
               ! Localize the drop on the ff mesh
               ff%lp%p(np)%ind =ff%lp%cfg%get_ijk_global(ff%lp%p(np)%pos,[ff%lp%cfg%imin,ff%lp%cfg%jmin,ff%lp%cfg%kmin])
               ! Activate it if it is large enough
               ff%lp%p(np)%flag=1
               if (ff%lp%p(np)%d.gt.diameter_threshold*ff%fs%dx) then
                  ff%lp%p(np)%flag=0
                  ncreated=ncreated+1
               end if
               ! Increment particle counter
               ff%lp%np_=np
            end if
            ! Remove drop from VOF representation
            do m=1,sd%ccl%struct(n)%n_
               ! Get cell index
               i=sd%ccl%struct(n)%map(1,m); j=sd%ccl%struct(n)%map(2,m); k=sd%ccl%struct(n)%map(3,m)
               ! Zero out Q(1) and Q(3)
               sd%fs%Q(i,j,k,1)=0.0_WP
               sd%fs%Q(i,j,k,3)=0.0_WP
               ! Rescale Q(2) and Q(4)
               if (sd%fs%VF(i,j,k).lt.0.99_WP) then
                  sd%fs%Q(i,j,k,2)=sd%fs%Q(i,j,k,2)/(1.0_WP-sd%fs%VF(i,j,k))
                  sd%fs%Q(i,j,k,4)=sd%fs%Q(i,j,k,4)/(1.0_WP-sd%fs%VF(i,j,k))
               else
                  ! Too little gas to trust local properties, rebuild from neighbors instead
                  sd%fs%Q(i,j,k,2)=sum((1.0_WP-sd%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1))*sd%fs%RHOG(i-1:i+1,j-1:j+1,k-1:k+1))/sum((1.0_WP-sd%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)))
                  sd%fs%Q(i,j,k,4)=sum((1.0_WP-sd%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1))*sd%fs%RHOG(i-1:i+1,j-1:j+1,k-1:k+1)*sd%fs%IG(i-1:i+1,j-1:j+1,k-1:k+1))/sum((1.0_WP-sd%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)))
               end if
               ! Remove liquid from volume moments
               sd%fs%VF  (i,j,k)=0.0_WP
               sd%fs%BL(:,i,j,k)=[sd%fs%cfg%xm(i),sd%fs%cfg%ym(j),sd%fs%cfg%zm(k)]
               sd%fs%BG(:,i,j,k)=[sd%fs%cfg%xm(i),sd%fs%cfg%ym(j),sd%fs%cfg%zm(k)]
               ! Adjust PLIC interface
               call setPlane(sd%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,sd%fs%VF(i,j,k)-0.5_WP))
               ! Remove polygon
               call zeroPolygon(sd%fs%interface_polygon(i,j,k))
            end do
         end if
      end do
      ! Resync the spray
      call ff%lp%sync()
      ! Communicate conserved variables and volume moments
      call sd%fs%cfg%sync(sd%fs%Q(:,:,:,1))
      call sd%fs%cfg%sync(sd%fs%Q(:,:,:,2))
      call sd%fs%cfg%sync(sd%fs%Q(:,:,:,3))
      call sd%fs%cfg%sync(sd%fs%Q(:,:,:,4))
      call sd%fs%sync_volume_moments()
      ! Synchronize interface
      call sd%fs%sync_interface()
      ! Rebuild momentum
      call sd%fs%get_momentum()
      ! Rebuild primitive variables
      call sd%fs%get_primitive()
      
      ! Deallocate memory
      deallocate(Vd,Md,Pd,Bd,Ud)
      
      ! Some messaging
      if (amRoot.and.nremoved.gt.0) then
         write(message,'(">> Transfered ",i0," structures from shockdrop, created ",i0," particles in farfield")') nremoved,ncreated 
         call log(trim(message))
         write(output_unit,*) trim(message)
      end if
      
   contains
      !> Function that identifies cells that need a label
      logical function make_label(i1,j1,k1)
         implicit none
         integer, intent(in) :: i1,j1,k1
         if (sd%fs%VF(i1,j1,k1).gt.0.0_WP) then; make_label=.true.; else; make_label=.false.; end if
      end function make_label
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
      !> Function that test closeness of a point X0 to edge of sd domain
      logical function close_to_edge(X0)
         implicit none
         real(WP), dimension(3), intent(in) :: X0
         close_to_edge=.false.
         ! X direction
         if (sd%fs%cfg%nx.gt.1) then
            if (X0(1).lt.edgelo(1)) close_to_edge=.true.
            if (X0(1).gt.edgehi(1)) close_to_edge=.true.
         end if
         ! Y direction
         if (sd%fs%cfg%ny.gt.1) then
            if (X0(2).lt.edgelo(2)) close_to_edge=.true.
            if (X0(2).gt.edgehi(2)) close_to_edge=.true.
         end if
         ! Z direction
         if (sd%fs%cfg%nz.gt.1) then
            if (X0(3).lt.edgelo(3)) close_to_edge=.true.
            if (X0(3).gt.edgehi(3)) close_to_edge=.true.
         end if
      end function close_to_edge
   end subroutine transfer
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      if (associated(sd)) then
         call sd%finalize()
         deallocate(sd)
      end if
      if (associated(sd2ff)) then
         call sd2ff%finalize()
         deallocate(sd2ff)
      end if
      if (associated(ff2sd)) then
         call ff2sd%finalize()
         deallocate(ff2sd)
      end if
      call ff%finalize()
      call ens_evt%finalize()
      call remesh_evt%finalize()
   end subroutine simulation_final
   
   
end module simulation
