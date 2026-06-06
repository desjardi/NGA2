!> SG-liquid + ideal-gas-mixture relaxation: mechanical (p), thermal (pT),
!> chemical/phase-change (pTg) with non-condensable gas
!> Liquid: class(stiffened_gas) (accepts SG or NASG via inheritance)
!> Gas: class(igmix) (multi-species mixture; ns=1 reduces to pure vapor)
!> Assumes vapor is the transported gas species at Q(7+liq%ns+indV-1); air (or other carrier) is the implicit ns-th species
module relax_igmix_sg_class
   use precision,             only: WP
   use messager,              only: die
   use thermorelax_class,     only: thermorelax
   use stiffened_gas_class,   only: stiffened_gas
   use igmix_class,           only: igmix
   implicit none
   private

   public :: relax_igmix_sg
   public :: Prelax,PTrelax,PTgrelax
   public :: Mv,Ma

   !> Molar masses of vapor and air [kg/mol]
   real(WP), parameter :: Mv=0.0180153_WP
   real(WP), parameter :: Ma=0.02897_WP

   !> Model enum
   integer, parameter :: Prelax  =1   !< Mechanical only
   integer, parameter :: PTrelax =2   !< Mechanical + thermal
   integer, parameter :: PTgrelax=3   !< Mechanical + thermal + chemical/phase change

   type, extends(thermorelax) :: relax_igmix_sg
      class(stiffened_gas), pointer :: liq => null()
      class(igmix),         pointer :: gas => null()
      !> Species indices for vapor and air in the gas mixture
      integer  :: indV =0
      integer  :: indA =0
      !> Saturation curve coefficients (computed in initialize from material properties)
      real(WP) :: AS=0.0_WP,BS=0.0_WP,CS=0.0_WP,DS=0.0_WP,ES=0.0_WP
      !> Convergence tolerances
      real(WP) :: p_tol     =1.0e-5_WP
      real(WP) :: Yv_tol    =1.0e-5_WP
      real(WP) :: Yv_tol_abs=1.0e-8_WP
      real(WP) :: Tsat_tol  =1.0e-5_WP
      real(WP) :: rho_tol   =1.0e-5_WP
      real(WP) :: rhoe_tol  =1.0e-5_WP
      real(WP) :: F1_tol    =1.0e-5_WP
      real(WP) :: F2_tol    =1.0e-5_WP
      !> Iteration limits
      integer  :: Tsat_itmax=40
      integer  :: NR_itmax  =40
      !> Dispatch
      integer  :: model=Prelax
   contains
      procedure :: initialize
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: get_T_lvg
      procedure :: get_coeffs_lv
      procedure :: get_p_eq
      procedure :: pTsat
      procedure :: dpTsatdT
      procedure :: dpTsatdp_lv
      procedure :: dpTsatdlnp
      procedure :: get_Tsat
      procedure :: get_pvsat
      procedure :: get_xv
   end type relax_igmix_sg

contains

   !> Store EOS pointers, vapor/air species indices, and compute saturation-curve coefficients
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_sg),         intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      real(WP) :: cpV,cvV,RV
      this%liq=>liq
      this%gas=>gas
      this%indV=indV
      this%indA=indA
      ! Vapor properties from the gas mixture (flat-array access)
      cvV=gas%cv(indV)
      cpV=gas%cp(indV)
      RV =cpV-cvV
      ! Saturation curve coefficients (integrated Clapeyron form)
      this%AS=(liq%cp-cpV+gas%qp(indV)-liq%qp)/RV
      this%BS=(liq%q -gas%q (indV))            /RV
      this%CS=(cpV-liq%cp)                     /RV
      this%DS=(liq%cp-liq%cv)                  /RV
      ! ES stays 0 for SG; NASG override sets ES=liq%b/RV after calling parent initialize
   end subroutine initialize

   !> Dispatch via model. Encapsulates mixture-cell gate.
   subroutine apply(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      ! Mixture cells only
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP) return
      ! Dispatch on model
      select case (this%model)
      case (Prelax);   call this%p_relax  (dt,VF,Q,Pjump)
      case (PTrelax);  call this%pT_relax (dt,VF,Q,Pjump)
      case (PTgrelax); call this%pTg_relax(dt,VF,Q,Pjump)
      case default;    call die('[relax_igmix_sg apply] unknown model')
      end select
   end subroutine apply

   !> Mechanical relaxation (Pelanti quadratic). Has clipping for unphysical phasic pressures.
   subroutine p_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP), dimension(:), allocatable   :: y
      real(WP) :: PL,PG,ZL,ZG,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! Vapor mass fraction (from Q layout)
      if (Q(2).gt.0.0_WP) then
         Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
      else
         Yv=0.0_WP
      end if
      ! Gas composition vector (vapor + carrier-by-closure)
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      ! Inline mixture coefficients
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=Q(1)/(       VF),e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
      ! Hard clipping for unphysical phasic pressure (cavitation, collapse)
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(7+this%liq%ns+this%indV-1)=Yv*Q(2)
         deallocate(y); return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(7+this%liq%ns+this%indV-1)=0.0_WP
         deallocate(y); return
      end if
      ! Phasic acoustic impedances (rho*c)
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(p=PL,rho=Q(1)/(       VF),y=[1.0_WP])
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(p=PG,rho=Q(2)/(1.0_WP-VF),y=y)
      ! Interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Quadratic for Peq
      coeffL=(this%liq%gamma-1.0_WP)*Pint+2.0_WP*this%liq%gamma*this%liq%pinf
      coeffG=(gammaG       -1.0_WP)*Pint
      a=1.0_WP+gammaG*VF+this%liq%gamma*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+gammaG)*VF*PL-(1.0_WP+this%liq%gamma)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      Peq =(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      VFeq=VF*((this%liq%gamma-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+this%liq%gamma)*Peq+coeffL)
      ! Update conservatives
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
      deallocate(y)
   end subroutine p_relax

   !> Mechanical + thermal relaxation. Calls p_relax first, then enforces TL=TG.
   subroutine pT_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: y
      real(WP) :: a,b,d,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! Step 1: mechanical
      call this%p_relax(dt,VF,Q,Pjump)
      ! Step 2: thermal
      if (Q(2).gt.0.0_WP) then
         Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Quadratic for the equilibrium pressure under PL=PG, TL=TG (SG liquid + IG mixture gas)
      a=Q(1)*this%liq%cv+Q(2)*cvG
      b=this%liq%q*this%liq%cv*(this%liq%gamma-1.0_WP)*Q(1)**2+qG*cvG*(gammaG-1.0_WP)*Q(2)**2+&
      &  Q(1)*this%liq%cv*this%liq%gamma*this%liq%pinf+Q(2)*cvG*this%liq%pinf                +&
      &  Q(1)*Q(2)*(this%liq%q*cvG*(gammaG-1.0_WP)+qG*this%liq%cv*(this%liq%gamma-1.0_WP))   -&
      &  sum(Q(3:4))*(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)+Q(2)*cvG*(gammaG-1.0_WP))
      d=cvG*(gammaG-1.0_WP)*this%liq%pinf*(qG*Q(2)**2+this%liq%q*Q(1)*Q(2)-sum(Q(3:4))*Q(2))
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y); return
      end if
      VFeq=Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq &
      &   /(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq+Q(2)*cvG*(gammaG-1.0_WP)*(Peq+this%liq%pinf))
      ! Clamp
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-this%liq%pinf); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,0.0_WP);          end if
      ! Update conservatives
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=Q(1)/max(VFeq,tiny(1.0_WP)),y=[1.0_WP])
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=Q(2)/max(1.0_WP-VFeq,tiny(1.0_WP)),y=y)
      VF=VFeq
      deallocate(y)
   end subroutine pT_relax

   !> Mechanical + thermal + chemical (phase change) relaxation.
   !> Nucleates a tiny opposite phase in metastable pure cells; calls pT_relax;
   !> runs pure-phase admissibility tests; falls back to LV (pure water) or LVG (with non-condensable) Newton solves.
   subroutine pTg_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: Q0,Qin,y
      real(WP) :: VF0,VFin,p,T,Yv
      real(WP) :: rho0,rhoe0,rhoA0
      real(WP) :: RHOL,RHOG
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP), parameter :: p_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,pv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: fd_eps=1.0e-7_WP,F_line_search_tol=0.3_WP
      logical :: chem_relax,nucleated
      allocate(Qin(size(Q))); Qin=Q
      VFin=VF
      nucleated=.false.
      ! Nucleation: seed a tiny opposite phase in metastable pure-ish cells so pTg starts well-conditioned
      nucleation: block
         real(WP), parameter :: VF_nuc=1.0e-7_WP
         real(WP) :: rhoL_nuc,pL_nuc,TL_nuc,pv_sat,rhoV_nuc,eV_nuc
         real(WP) :: rhoG_nuc,pG_nuc,TG_nuc,Yv_nuc,xv_nuc,pv_nuc,Tsat_nuc
         real(WP) :: rhoL_new,eL_new,drho,de
         real(WP) :: y_nuc(this%gas%ns)
         logical  :: conv_nuc
         integer  :: Tsat_it_nuc
         if (VF.ge.1.0_WP-VF_nuc) then
            ! Near-pure liquid: check cavitation
            rhoL_nuc=Q(1)/VF
            pL_nuc=this%liq%get_p_from_rho_e(rho=rhoL_nuc,e=Q(3)/Q(1),y=[1.0_WP])
            TL_nuc=this%liq%get_T_from_p_rho(p=pL_nuc,rho=rhoL_nuc,y=[1.0_WP])
            if (pL_nuc.le.-this%liq%pinf.or.TL_nuc.le.0.0_WP) return
            pv_sat=this%get_pvsat(pL_nuc,TL_nuc)
            if (pv_sat.le.pL_nuc) return  ! stable pure liquid
            ! Superheated liquid: nucleate tiny vapor
            y_nuc=0.0_WP; y_nuc(this%indV)=1.0_WP
            rhoV_nuc=this%gas%get_rho_from_p_T(p=pL_nuc,T=TL_nuc,y=y_nuc)
            eV_nuc  =this%gas%get_e_from_p_T  (p=pL_nuc,T=TL_nuc,y=y_nuc)
            drho=VF_nuc*rhoV_nuc; de=drho*eV_nuc
            Q(1)=Q(1)-drho; Q(2)=Q(2)+drho
            Q(3)=Q(3)-de;   Q(4)=Q(4)+de
            Q(7+this%liq%ns+this%indV-1)=Q(7+this%liq%ns+this%indV-1)+drho
            VF=1.0_WP-VF_nuc
            nucleated=.true.
         else if (VF.le.VF_nuc) then
            ! Near-pure gas: check condensation
            if (Q(2).le.0.0_WP) return
            Yv_nuc=Q(7+this%liq%ns+this%indV-1)/Q(2); Yv_nuc=max(Yvmin,min(Yvmax,Yv_nuc))
            if (Yv_nuc.le.Yv_dry) return
            y_nuc(this%indV)=Yv_nuc; y_nuc(this%indA)=1.0_WP-Yv_nuc
            rhoG_nuc=Q(2)/max(1.0_WP-VF,tiny(1.0_WP))
            pG_nuc=this%gas%get_p_from_rho_e(rho=rhoG_nuc,e=Q(4)/Q(2),y=y_nuc)
            TG_nuc=this%gas%get_T_from_p_rho(p=pG_nuc,rho=rhoG_nuc,y=y_nuc)
            if (pG_nuc.le.0.0_WP.or.TG_nuc.le.0.0_WP) return
            xv_nuc=this%get_xv(Yv_nuc); pv_nuc=xv_nuc*pG_nuc
            if (pv_nuc.le.p_eps) return
            call this%get_Tsat(pG_nuc,pv_nuc,TG_nuc,Tsat_nuc,conv_nuc,Tsat_it_nuc)
            if (.not.conv_nuc) return
            if (TG_nuc.ge.Tsat_nuc) return  ! stable pure vapor/gas
            ! Supersaturated: nucleate tiny liquid
            rhoL_new=this%liq%get_rho_from_p_T(p=pG_nuc,T=TG_nuc,y=[1.0_WP])
            eL_new  =this%liq%get_e_from_p_T  (p=pG_nuc,T=TG_nuc,y=[1.0_WP])
            drho=VF_nuc*rhoL_new
            drho=min(drho,0.5_WP*Q(7+this%liq%ns+this%indV-1),0.5_WP*Q(2))
            if (drho.le.0.0_WP) return
            de=drho*eL_new
            Q(1)=Q(1)+drho; Q(2)=Q(2)-drho
            Q(3)=Q(3)+de;   Q(4)=Q(4)-de
            Q(7+this%liq%ns+this%indV-1)=Q(7+this%liq%ns+this%indV-1)-drho
            VF=drho/rhoL_new
            nucleated=.true.
         end if
      end block nucleation
      ! Steps 1+2: mechanical + thermal
      call this%pT_relax(dt,VF,Q,Pjump)
      ! Step 3: chemical (phase change)
      if (Q(2).gt.0.0_WP) then
         Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Recover p, T from dominant phase
      if (VF.gt.0.5_WP) then
         p=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
         T=this%liq%get_T_from_p_rho(p=p,rho=Q(1)/VF,y=[1.0_WP])
      else
         p=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
         T=this%gas%get_T_from_p_rho(p=p,rho=Q(2)/(1.0_WP-VF),y=y)
      end if
      ! Conserve totals
      VF0=VF
      allocate(Q0(size(Q))); Q0=Q
      rho0 =sum(Q0(1:2))
      rhoe0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      ! Pure-phase admissibility tests (Caze et al. trigger)
      pure_phase_bounds: block
         real(WP), parameter :: rhoA_pure=1.0e-12_WP
         real(WP) :: rhoL_pure,eL_pure,pL_pure,TL_pure,Tsat_pure
         real(WP) :: rhoV_pure,eV_pure,pV_pure,TV_pure
         real(WP) :: rhoV0,Yv_gas,pG_pure,TG_pure,xv_gas,pv_gas
         integer  :: Tsat_it_pure
         logical  :: conv_pure
         if (rho0.le.0.0_WP.or.rhoe0.le.0.0_WP) exit pure_phase_bounds
         if (rhoA0/rho0.le.rhoA_pure) then
            ! No air content: literal pure-liquid / pure-vapor tests
            rhoL_pure=rho0; eL_pure=rhoe0/rho0
            pL_pure=this%liq%get_p_from_rho_e(rho=rhoL_pure,e=eL_pure,y=[1.0_WP])
            TL_pure=this%liq%get_T_from_p_rho(p=pL_pure,rho=rhoL_pure,y=[1.0_WP])
            if (pL_pure.gt.p_eps.and.TL_pure.gt.0.0_WP) then
               call this%get_Tsat(pL_pure,pL_pure,TL_pure,Tsat_pure,conv_pure,Tsat_it_pure)
               if (conv_pure.and.TL_pure.le.Tsat_pure*(1.0_WP+this%Tsat_tol)) then
                  VF=1.0_WP
                  Q(1)=rho0;     Q(2)=0.0_WP
                  Q(3)=rhoe0;    Q(4)=0.0_WP
                  Q(7+this%liq%ns+this%indV-1)=0.0_WP
                  call dealloc(); return
               end if
            end if
            rhoV_pure=rho0; eV_pure=rhoe0/rho0
            y=0.0_WP; y(this%indV)=1.0_WP
            pV_pure=this%gas%get_p_from_rho_e(rho=rhoV_pure,e=eV_pure,y=y)
            TV_pure=this%gas%get_T_from_p_rho(p=pV_pure,rho=rhoV_pure,y=y)
            if (pV_pure.gt.p_eps.and.TV_pure.gt.0.0_WP) then
               call this%get_Tsat(pV_pure,pV_pure,TV_pure,Tsat_pure,conv_pure,Tsat_it_pure)
               if (conv_pure.and.TV_pure.ge.Tsat_pure*(1.0_WP-this%Tsat_tol)) then
                  VF=0.0_WP
                  Q(1)=0.0_WP;   Q(2)=rho0
                  Q(3)=0.0_WP;   Q(4)=rhoe0
                  Q(7+this%liq%ns+this%indV-1)=rho0
                  call dealloc(); return
               end if
            end if
         else
            ! With non-condensable air: only check all-water-as-vapor gas state
            rhoV0=rho0-rhoA0
            if (rhoV0.le.0.0_WP) exit pure_phase_bounds
            Yv_gas=rhoV0/rho0; Yv_gas=max(Yvmin,min(Yvmax,Yv_gas))
            y=0.0_WP; y(this%indV)=Yv_gas; y(this%indA)=1.0_WP-Yv_gas
            pG_pure=this%gas%get_p_from_rho_e(rho=rho0,e=rhoe0/rho0,y=y)
            TG_pure=this%gas%get_T_from_p_rho(p=pG_pure,rho=rho0,y=y)
            if (pG_pure.gt.p_eps.and.TG_pure.gt.0.0_WP) then
               xv_gas=this%get_xv(Yv_gas); pv_gas=xv_gas*pG_pure
               if (pv_gas.gt.p_eps) then
                  call this%get_Tsat(pG_pure,pv_gas,TG_pure,Tsat_pure,conv_pure,Tsat_it_pure)
                  if (conv_pure.and.TG_pure.ge.Tsat_pure*(1.0_WP-this%Tsat_tol)) then
                     VF=0.0_WP
                     Q(1)=0.0_WP;   Q(2)=rho0
                     Q(3)=0.0_WP;   Q(4)=rhoe0
                     Q(7+this%liq%ns+this%indV-1)=rhoV0
                     call dealloc(); return
                  end if
               end if
            end if
         end if
      end block pure_phase_bounds
      ! Activation: decide if chemical relaxation should fire
      chem_relax=activate_chem(p,T,Yv)
      if (.not.chem_relax) then
         call restore(); call dealloc(); return
      end if
      ! Solve chemical equilibrium for p, T, Yv (without modifying Q yet)
      if (Yv.gt.Yv_pure) then
         Yv=Yvmax
         y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
         cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
         qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
         gammaG=cpG/cvG
         call solve_lv(p,T,chem_relax)
      else
         call solve_lvg(p,T,Yv,chem_relax)
      end if
      if (.not.chem_relax) then
         if (nucleated) then
            VF=VFin; Q=Qin
         else
            call restore()
         end if
         call dealloc(); return
      end if
      ! Update Q with the converged equilibrium state
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      RHOL=this%liq%get_rho_from_p_T(p=p,T=T,y=[1.0_WP])
      RHOG=this%gas%get_rho_from_p_T(p=p,T=T,y=y)
      VF=(rho0-RHOG)/(RHOL-RHOG)
      if (VF.lt.0.0_WP) then; VF=0.0_WP; p=max(p,-this%liq%pinf); end if
      if (VF.gt.1.0_WP) then; VF=1.0_WP; p=max(p,0.0_WP);          end if
      Q(1)=(       VF)*RHOL
      Q(2)=(1.0_WP-VF)*RHOG
      Q(3)=Q(1)*this%liq%get_e_from_p_T(p=p,T=T,y=[1.0_WP])
      Q(4)=Q(2)*this%gas%get_e_from_p_T(p=p,T=T,y=y)
      Q(7+this%liq%ns+this%indV-1)=Q(2)*Yv
      if (.not.check_cons()) then
         call restore(); call dealloc(); return
      end if
      call dealloc()
   contains
      subroutine restore()
         VF=VF0; Q=Q0
      end subroutine restore
      subroutine dealloc()
         if (allocated(Q0)) deallocate(Q0)
         if (allocated(y))  deallocate(y)
         if (allocated(Qin))deallocate(Qin)
      end subroutine dealloc
      logical function check_pv(pv_)
         real(WP), intent(in) :: pv_
         check_pv=(pv_.gt.p_eps)
      end function check_pv
      logical function check_cons()
         real(WP) :: re,ee
         re=(sum(Q(1:2))-rho0)/rho0
         ee=(sum(Q(3:4))-rhoe0)/rhoe0
         check_cons=(abs(re).le.this%rho_tol).and.(abs(ee).le.this%rhoe_tol)
      end function check_cons
      logical function activate_chem(p_,T_,Yv_)
         real(WP), intent(in)    :: p_,T_
         real(WP), intent(inout) :: Yv_
         real(WP) :: xv,pv_,Fsat
         activate_chem=.false.
         xv=this%get_xv(Yv_); pv_=xv*p_
         if ((Yv_.le.Yv_dry).or.(pv_.le.pv_dry).or.(.not.check_pv(pv_))) then
            ! Dry / ill-conditioned edge: seed Yv from saturation at current state
            pv_=exp(this%AS+(this%BS+this%ES*p_)/T_)*T_**this%CS*(p_+this%liq%pinf)**this%DS
            if (.not.check_pv(pv_)) return
            if (pv_.ge.p_) then
               Yv_=sqrt(0.0001_WP)
            else
               xv=pv_/p_; Yv_=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
            end if
            Yv_=max(Yvmin,min(Yvmax,Yv_))
         else
            ! Direct saturation residual
            Fsat=this%pTsat(p_,pv_,T_)
            if (abs(Fsat).lt.this%F1_tol) return
         end if
         activate_chem=.true.
      end function activate_chem
      real(WP) function get_T_lv(ap,bp,dp)
         real(WP), intent(in) :: ap,bp,dp
         get_T_lv=(-bp+sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
      end function get_T_lv
      real(WP) function get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
         real(WP), intent(in) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         get_dTdp_lv=(ap*(-dbpdp+(bp*dbpdp-2.0_WP*(dapdp*dp+ap*ddpdp))/sqrt(bp**2-4.0_WP*ap*dp))-dapdp*(-bp+sqrt(bp**2-4.0_WP*ap*dp)))/(2.0_WP*ap**2)
      end function get_dTdp_lv
      real(WP) function rhoe_res_lvg(p_,T_,Yv_)
         real(WP), intent(in) :: p_,T_,Yv_
         rhoe_res_lvg=(rho0*(1.0_WP-Yv_)-rhoA0)*this%liq%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP])+&
         &             rhoA0*this%gas%get_e_from_p_T(p=p_,T=T_,y=[Yv_,1.0_WP-Yv_])           -&
         &             rhoe0*(1.0_WP-Yv_)
      end function rhoe_res_lvg
      real(WP) function dlnxvdYv(Yv_)
         real(WP), intent(in) :: Yv_
         real(WP) :: Ys,den
         Ys=max(Yv_,Yvmin+fd_eps)
         den=Ys*Ma+(1.0_WP-Ys)*Mv
         dlnxvdYv=1.0_WP/Ys-(Ma-Mv)/den
      end function dlnxvdYv
      !> LV solve: damped Newton in ln(p)
      subroutine solve_lv(p_eq,T_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq
         logical,  intent(out)   :: conv
         real(WP) :: pOld,lnpOld,p_try,T_try
         real(WP) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         real(WP) :: dTdp,dTdlnp,dF1dlnp
         real(WP) :: F1,F1_try,dlnp_nr,p_err,alpha
         integer  :: it
         logical  :: accepted
         conv=.false.; p_err=10.0_WP*this%p_tol
         do it=1,this%NR_itmax
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp); if (T_eq.le.0.0_WP) return
            dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp); dTdlnp=p_eq*dTdp
            F1=this%pTsat(p_eq,p_eq,T_eq); dF1dlnp=this%dpTsatdlnp(p_eq,T_eq,dTdlnp)
            if (abs(dF1dlnp).lt.1.0e-30_WP) exit
            dlnp_nr=-F1/dF1dlnp
            dlnp_nr=max(log(0.5_WP),min(log(1.5_WP),dlnp_nr))
            pOld=p_eq; lnpOld=log(pOld); alpha=1.0_WP; accepted=.false.
            do while (alpha.gt.1.0e-8_WP)
               p_try=exp(lnpOld+alpha*dlnp_nr)
               if (p_try.le.p_eps) then; alpha=0.5_WP*alpha; cycle; end if
               call this%get_coeffs_lv(p_try,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
               T_try=get_T_lv(ap,bp,dp)
               if (T_try.le.0.0_WP) then; alpha=0.5_WP*alpha; cycle; end if
               F1_try=this%pTsat(p_try,p_try,T_try)
               if (abs(F1_try).lt.abs(F1)) then
                  p_eq=p_try; T_eq=T_try; accepted=.true.; exit
               end if
               alpha=0.5_WP*alpha
            end do
            if (.not.accepted) exit
            p_err=abs(log(p_eq/pOld))
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp); if (T_eq.le.0.0_WP) return
            F1=this%pTsat(p_eq,p_eq,T_eq)
            if ((p_err.lt.this%p_tol).and.(abs(F1).lt.this%F1_tol)) then
               conv=.true.; exit
            end if
         end do
         if (.not.conv) return
         call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
         T_eq=get_T_lv(ap,bp,dp)
      end subroutine solve_lv
      !> LVG solve: 2x2 damped Newton in (ln(p), Yv) with step limiter
      subroutine solve_lvg(p_eq,T_eq,Yv_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         logical,  intent(out)   :: conv
         real(WP) :: xv,pv
         real(WP) :: F1,F2,F2p,F2Y,dF1dlnp,dF1dYv,dF2dlnp,dF2dYv,detJ
         real(WP) :: lnp_pert,p_pert,Yv_pert,T_pert,xv_pert,pv_pert,dTdlnp,dTdYv,dlnp_nr,dYv_nr
         real(WP) :: pOld,YvOld,lnpOld,p_err,Yv_err
         real(WP) :: alpha,res0,res_try
         real(WP) :: p_try,Yv_try,T_try,xv_try,pv_try,F1_try,F2_try
         real(WP) :: Yv_max_phys,Yv_hi
         integer  :: it,lsit
         logical  :: accepted
         Yv_max_phys=1.0_WP-(rhoA0/rho0)
         Yv_hi=min(Yvmax,Yv_max_phys)
         conv=.false.
         p_err=10.0_WP*this%p_tol; Yv_err=10.0_WP*this%Yv_tol
         do it=1,this%NR_itmax
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0); if (T_eq.le.0.0_WP) return
            xv=this%get_xv(Yv_eq); pv=xv*p_eq
            if (.not.check_pv(pv)) return
            F1=this%pTsat(p_eq,pv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            res0=sqrt(F1**2+F2**2)
            ! d/dlnp via forward difference
            lnp_pert=log(p_eq)+fd_eps; p_pert=exp(lnp_pert); pv_pert=xv*p_pert
            if (.not.check_pv(pv_pert)) return
            T_pert=this%get_T_lvg(p_pert,Yv_eq,rho0,rhoA0); if (T_pert.le.0.0_WP) return
            F2p=rhoe_res_lvg(p_pert,T_pert,Yv_eq)/rhoe0
            dTdlnp=(T_pert-T_eq)/fd_eps
            dF1dlnp=this%dpTsatdlnp(p_eq,T_eq,dTdlnp)
            dF2dlnp=(F2p-F2)/fd_eps
            ! d/dYv via forward (or backward) difference, staying inside domain
            Yv_pert=Yv_eq+fd_eps
            if (Yv_pert.gt.Yv_hi-fd_eps) Yv_pert=Yv_eq-fd_eps
            if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
            if ((Yv_pert.le.Yvmin+fd_eps).or.(Yv_pert.ge.Yv_hi-fd_eps)) return
            xv_pert=this%get_xv(Yv_pert); pv_pert=xv_pert*p_eq
            if (.not.check_pv(pv_pert)) return
            T_pert=this%get_T_lvg(p_eq,Yv_pert,rho0,rhoA0); if (T_pert.le.0.0_WP) return
            F2Y=rhoe_res_lvg(p_eq,T_pert,Yv_pert)/rhoe0
            dTdYv=(T_pert-T_eq)/(Yv_pert-Yv_eq)
            dF1dYv=this%dpTsatdT(p_eq,T_eq)*dTdYv-dlnxvdYv(Yv_eq)
            dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
            ! 2x2 Newton
            detJ=dF1dlnp*dF2dYv-dF1dYv*dF2dlnp
            if (abs(detJ).lt.1.0e-30_WP) exit
            dlnp_nr=-( dF2dYv *F1-dF1dYv  *F2)/detJ
            dYv_nr =-(-dF2dlnp*F1+dF1dlnp*F2)/detJ
            ! Direction-preserving step limiter
            step_limit: block
               real(WP) :: ms,lnp_up,lnp_dn
               ms=1.0_WP
               lnp_up=log(1.5_WP); lnp_dn=log(0.5_WP)
               if (dlnp_nr.gt.lnp_up) ms=min(ms,lnp_up/dlnp_nr)
               if (dlnp_nr.lt.lnp_dn) ms=min(ms,lnp_dn/dlnp_nr)
               if (dYv_nr.gt.0.0_WP) then
                  if (Yv_eq+dYv_nr.ge.Yv_hi-fd_eps) ms=min(ms,0.9_WP*(Yv_hi-fd_eps-Yv_eq)/dYv_nr)
               else if (dYv_nr.lt.0.0_WP) then
                  if (Yv_eq+dYv_nr.le.Yvmin+fd_eps) ms=min(ms,0.9_WP*(Yv_eq-Yvmin-fd_eps)/abs(dYv_nr))
               end if
               ms=max(0.0_WP,min(1.0_WP,ms))
               dlnp_nr=dlnp_nr*ms; dYv_nr=dYv_nr*ms
            end block step_limit
            ! Damped update with line search
            pOld=p_eq; YvOld=Yv_eq; lnpOld=log(pOld); alpha=1.0_WP; lsit=0
            if ((abs(F1).lt.F_line_search_tol).and.(abs(F2).lt.F_line_search_tol)) then
               p_eq=exp(lnpOld+dlnp_nr); Yv_eq=YvOld+dYv_nr
            else
               accepted=.false.
               do while (alpha.gt.1.0e-8_WP)
                  lsit=lsit+1
                  p_try=exp(lnpOld+alpha*dlnp_nr); Yv_try=YvOld+alpha*dYv_nr
                  if (p_try.le.p_eps) then; alpha=0.5_WP*alpha; cycle; end if
                  if ((Yv_try.le.Yvmin+fd_eps).or.(Yv_try.ge.Yv_hi-fd_eps)) then
                     alpha=0.5_WP*alpha; cycle
                  end if
                  if ((rho0*(1.0_WP-Yv_try)-rhoA0).le.0.0_WP) then
                     alpha=0.5_WP*alpha; cycle
                  end if
                  T_try=this%get_T_lvg(p_try,Yv_try,rho0,rhoA0)
                  if (T_try.le.0.0_WP) then; alpha=0.5_WP*alpha; cycle; end if
                  xv_try=this%get_xv(Yv_try); pv_try=xv_try*p_try
                  if (.not.check_pv(pv_try)) then; alpha=0.5_WP*alpha; cycle; end if
                  F1_try=this%pTsat(p_try,pv_try,T_try)
                  F2_try=rhoe_res_lvg(p_try,T_try,Yv_try)/rhoe0
                  res_try=sqrt(F1_try**2+F2_try**2)
                  if (res_try.lt.res0) then
                     p_eq=p_try; Yv_eq=Yv_try; T_eq=T_try; accepted=.true.; exit
                  end if
                  alpha=0.5_WP*alpha
               end do
               if (.not.accepted) exit
            end if
            p_err=abs(log(p_eq/pOld))
            Yv_err=abs(Yv_eq-YvOld)
            xv=this%get_xv(Yv_eq); pv=xv*p_eq
            if (.not.check_pv(pv)) return
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0); if (T_eq.le.0.0_WP) return
            F1=this%pTsat(p_eq,pv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            if ((p_err.lt.this%p_tol).and.Yv_err.lt.this%Yv_tol_abs+this%Yv_tol*max(abs(YvOld),abs(Yv_eq)).and.(abs(F1).lt.this%F1_tol).and.(abs(F2).lt.this%F2_tol)) then
               conv=.true.; exit
            end if
         end do
      end subroutine solve_lvg
   end subroutine pTg_relax

   !> p-T saturation residual (general form: ES=0 for SG, ES=b/RV for NASG)
   real(WP) function pTsat(this,pl_,pv_,T_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: pl_,pv_,T_
      pTsat=this%AS+(this%BS+this%ES*pl_)/T_+this%CS*log(T_)+this%DS*log(pl_+this%liq%pinf)-log(pv_)
   end function pTsat

   real(WP) function dpTsatdT(this,pl_,T_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: pl_,T_
      dpTsatdT=-(this%BS+this%ES*pl_)/T_**2+this%CS/T_
   end function dpTsatdT

   real(WP) function dpTsatdp_lv(this,p_,T_,dTdp_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdp_
      dpTsatdp_lv=this%dpTsatdT(p_,T_)*dTdp_+this%ES/T_+this%DS/(p_+this%liq%pinf)-1.0_WP/p_
   end function dpTsatdp_lv

   real(WP) function dpTsatdlnp(this,p_,T_,dTdlnp_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdlnp_
      dpTsatdlnp=this%dpTsatdT(p_,T_)*dTdlnp_+this%ES*p_/T_+this%DS*p_/(p_+this%liq%pinf)-1.0_WP
   end function dpTsatdlnp

   !> Safeguarded Newton on the saturation curve for Tsat(pl, pv)
   subroutine get_Tsat(this,pl_,pv_,Tguess,Tsat,conv,Tsat_it)
      implicit none
      class(relax_igmix_sg), intent(inout) :: this
      real(WP), intent(in)  :: pl_,pv_,Tguess
      real(WP), intent(out) :: Tsat
      logical,  intent(out) :: conv
      integer,  intent(out) :: Tsat_it
      real(WP) :: Tlo,Thi,Told,Tnew,Flo,Fhi,Fold,Fnew,dFold
      integer  :: it,expand_it
      conv=.false.; Tsat_it=0
      Tlo=250.0_WP; Thi=900.0_WP
      Flo=this%pTsat(pl_,pv_,Tlo)
      Fhi=this%pTsat(pl_,pv_,Thi)
      expand_it=0
      do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
         if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
            Tlo=max(1.0_WP,0.8_WP*Tlo); Flo=this%pTsat(pl_,pv_,Tlo)
         else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
            Thi=1.2_WP*Thi;             Fhi=this%pTsat(pl_,pv_,Thi)
         else
            exit
         end if
         expand_it=expand_it+1
      end do
      if (Flo*Fhi.gt.0.0_WP) return
      Tsat=max(Tlo,min(Thi,Tguess))
      do it=1,this%Tsat_itmax
         Told=Tsat
         Fold=this%pTsat(pl_,pv_,Told); dFold=this%dpTsatdT(pl_,Told)
         if (abs(dFold).gt.tiny(1.0_WP)) then
            Tnew=Told-Fold/dFold
         else
            Tnew=0.5_WP*(Tlo+Thi)
         end if
         if ((Tnew.ne.Tnew).or.(Tnew.le.Tlo).or.(Tnew.ge.Thi)) Tnew=0.5_WP*(Tlo+Thi)
         Fnew=this%pTsat(pl_,pv_,Tnew)
         if (Fnew.ne.Fnew) then
            Tnew=0.5_WP*(Tlo+Thi); Fnew=this%pTsat(pl_,pv_,Tnew)
         end if
         if (Flo*Fnew.le.0.0_WP) then
            Thi=Tnew; Fhi=Fnew
         else
            Tlo=Tnew; Flo=Fnew
         end if
         Tsat_it=it; Tsat=Tnew
         if ((abs((Tnew-Told)/max(abs(Told),tiny(1.0_WP))).lt.this%Tsat_tol).or.(abs(Fnew).lt.this%F1_tol)) then
            conv=.true.; return
         end if
      end do
   end subroutine get_Tsat

   real(WP) function get_pvsat(this,pl_,T_)
      implicit none
      class(relax_igmix_sg), intent(inout) :: this
      real(WP), intent(in)  :: pl_,T_
      get_pvsat=exp(this%AS+(this%BS+this%ES*pl_)/T_)*T_**this%CS*(pl_+this%liq%pinf)**this%DS
   end function get_pvsat

   real(WP) function get_xv(this,Yv_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: Yv_
      get_xv=Yv_*Ma/(Yv_*Ma+(1.0_WP-Yv_)*Mv)
   end function get_xv

   !> Equilibrium T from energy conservation in liquid-vapor-gas mixture (SG form: no co-volume)
   real(WP) function get_T_lvg(this,p_,Yv_,rho0,rhoA0)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: p_,Yv_,rho0,rhoA0
      get_T_lvg=(1.0_WP-Yv_)/((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)+&
      &           rhoA0*((this%gas%gamma(this%indV)-1.0_WP)*this%gas%cv(this%indV)*Yv_+ &
      &                  (this%gas%gamma(this%indA)-1.0_WP)*this%gas%cv(this%indA)*(1.0_WP-Yv_))/p_)
   end function get_T_lvg

   !> Quadratic coefficients for the equilibrium-T equation (SG form: PinfG=0)
   subroutine get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
      implicit none
      class(relax_igmix_sg), intent(in)  :: this
      real(WP), intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      ap=rho0*this%liq%cv*cvG*((GammaG-1.0_WP)*(p_eq+this%liq%gamma*this%liq%pinf)-(this%liq%gamma-1.0_WP)*p_eq)
      bp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv*p_eq-(GammaG-1.0_WP)*cvG*(p_eq+this%liq%pinf))              +&
      &  rho0*((GammaG-1.0_WP)*cvG*this%liq%q*(p_eq+this%liq%pinf)-(this%liq%gamma-1.0_WP)*this%liq%cv*qG*p_eq) +&
      &  cvG*p_eq*(p_eq+this%liq%pinf)-this%liq%cv*p_eq*(p_eq+this%liq%gamma*this%liq%pinf)
      dp=(qG-this%liq%q)*(p_eq+this%liq%pinf)*p_eq
      dapdp=rho0*this%liq%cv*cvG*(GammaG-this%liq%gamma)
      dbpdp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(GammaG-1.0_WP)*cvG)                                     +&
      &     rho0*((GammaG-1.0_WP)*cvG*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qG)                        +&
      &     cvG*(2.0_WP*p_eq+this%liq%pinf)-this%liq%cv*(2.0_WP*p_eq+this%liq%gamma*this%liq%pinf)
      ddpdp=(qG-this%liq%q)*(2.0_WP*p_eq+this%liq%pinf)
   end subroutine get_coeffs_lv

   !> SG-form energy-conserving equilibrium pressure at given VF (used by NASG override)
   real(WP) function get_p_eq(this,VF_,Q0_,qG_,gammaG_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP),               intent(in) :: VF_
      real(WP), dimension(:), intent(in) :: Q0_
      real(WP),               intent(in) :: qG_,gammaG_
      get_p_eq=(sum(Q0_(3:4))-Q0_(1)*this%liq%q-VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0_(2)*qG_)/&
      &        (VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG_-1.0_WP))
   end function get_p_eq

end module relax_igmix_sg_class
