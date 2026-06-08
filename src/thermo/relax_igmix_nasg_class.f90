!> NASG-liquid + ideal-gas-mixture relaxation (extends SG-IG mixture with co-volume b corrections)
!> Inherits pTg_relax unchanged from parent; virtual dispatch into the overridden helpers
!> (relax_p, pT_relax, get_T_lvg, get_coeffs_lv, get_p_eq) does the right thing.
!> References Pelanti 2022 https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
module relax_igmix_nasg_class
   use precision,               only: WP
   use relax_igmix_sg_class,    only: relax_igmix_sg
   use stiffened_gas_class,     only: stiffened_gas
   use nasg_class,              only: nasg
   use igmix_class,             only: igmix
   implicit none
   private

   public :: relax_igmix_nasg

   type, extends(relax_igmix_sg) :: relax_igmix_nasg
      !> Typed pointer for direct b access (avoids select type in hot path)
      type(nasg), pointer :: liq_nasg=>null()
      !> matm-style p/pT relaxation parameters (mirror relax_ig_nasg)
      real(WP) :: RHOGmin=1.0e-2_WP   !< Skip mechanical relax when gas density falls below this
      real(WP) :: phist  =1.0_WP      !< Temporal weighting on equilibrium pressure
      real(WP) :: phi0   =0.0_WP      !< Temporal weighting on interface pressure
   contains
      procedure :: initialize
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: get_T_lvg
      procedure :: get_coeffs_lv
      procedure :: get_p_eq
   end type relax_igmix_nasg

contains

   !> Initialize: call parent (sets liq, gas, indV/A, AS-DS, ES=0), then set NASG-specific fields
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_nasg),       intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      real(WP) :: cpV,cvV,RV
      ! Parent: sets this%liq=>liq, this%gas=>gas, indV/A, saturation coeffs AS-DS, ES=0
      call this%relax_igmix_sg%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      ! Set typed pointer and override ES (select type to extract type(nasg) from class(stiffened_gas))
      select type (liq)
      type is (nasg)
         this%liq_nasg=>liq
         cpV=gas%cp(indV); cvV=gas%cv(indV); RV=cpV-cvV
         this%ES=liq%b/RV
      end select
   end subroutine initialize

   !> Mechanical relaxation: matm-style quadratic (PL=Peq, PG=Peq-Pjump), generalized to an
   !> ideal-gas MIXTURE via mass-fraction-weighted (gammaG,cvG,qG) from the frozen composition,
   !> + NASG liquid (co-volume b). Reduces exactly to relax_ig_nasg%p_relax at single-component
   !> gas. Composition (Q(iVQ)) is untouched (no mass transfer in p-relax).
   subroutine p_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_nasg), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ,bL,rhoL,rhoG
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG,gL,pinfL,qL,Yv
      real(WP) :: y(this%gas%ns)
      integer  :: iVQ
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.this%RHOGmin) return
      ! Frozen gas composition -> mixture-effective ideal-gas parameters
      iVQ=7+this%liq%ns+this%indV-1
      Yv=Q(iVQ)/Q(2)
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Liquid parameters (b via the typed nasg pointer; EOS calls dispatch to nasg)
      bL=this%liq_nasg%b; gL=this%liq%gamma; pinfL=this%liq%pinf; qL=this%liq%q
      rhoL=Q(1)/(VF); rhoG=Q(2)/(1.0_WP-VF)
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=rhoG,e=Q(4)/Q(2),y=y)
      ! No cavitation model: leave sub-vacuum / over-packed cells untouched (avoid NaN impedances)
      if (PL.le.-pinfL.or.PG.le.0.0_WP.or.1.0_WP-bL*rhoL.le.0.0_WP) return
      ! Phasic acoustic impedances (rho*c)
      ZL=rhoL*this%liq%get_c_from_p_rho(p=PL,rho=rhoL,y=[1.0_WP])
      ZG=rhoG*this%gas%get_c_from_p_rho(p=PG,rho=rhoG,y=y)
      cJ=ZL/(ZG+ZL)
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Quadratic (pinf_g=0; liquid co-volume b in n1,n0; mixture gas gamma/q)
      n1=VF*this%phist+Q(1)*bL/(gL-1.0_WP)
      n0=VF*(this%phi0*Pint-this%phist*cJ*Pjump)+Q(3)-Q(1)*qL+Q(1)*bL*gL/(gL-1.0_WP)*pinfL
      d1=this%phist+1.0_WP/(gL-1.0_WP)
      d0=this%phi0*Pint-this%phist*cJ*Pjump+gL/(gL-1.0_WP)*pinfL
      a=d1*(1.0_WP/(gammaG-1.0_WP)+this%phist*VF)+n1*(-1.0_WP/(gammaG-1.0_WP)-this%phist)
      b=d1*(-Pjump/(gammaG-1.0_WP)+Q(2)*qG-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n1*(Pjump/(gammaG-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)+d0*(1.0_WP/(gammaG-1.0_WP)+this%phist*VF)+n0*(-1.0_WP/(gammaG-1.0_WP)-this%phist)
      d=d0*(-Pjump/(gammaG-1.0_WP)+Q(2)*qG-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n0*(Pjump/(gammaG-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g=0)
      if (Peq.le.-pinfL.or.Peq-Pjump.le.0.0_WP) return
      ! Equilibrium VF (strict bounds)
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) return
      ! Update Q with p*dV work at the relaxed interface pressure
      Q(3)=Q(3)-(this%phist*Peq+this%phi0*Pint-this%phist*cJ*Pjump)*(VFeq-VF)
      Q(4)=Q(4)+(this%phist*Peq+this%phi0*Pint-this%phist*cJ*Pjump)*(VFeq-VF)
      VF=VFeq
   end subroutine p_relax

   !> Mechanical + thermal relaxation: matm-style exact quadratic (PL=Peq, PG=Peq-Pjump, TL=TG),
   !> ideal-gas MIXTURE (mass-fraction-weighted gammaG,cvG,qG) + NASG liquid (co-volume via the
   !> (1-b*Q(1)) factor). Reduces exactly to relax_ig_nasg%pT_relax at single-component gas.
   subroutine pT_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_nasg), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP) :: a,b,d,Peq,VFeq,Eth
      real(WP) :: cv1,cv2,g1,g2,pinf,R1,R2,bL,ombm,cpG,qG,Yv
      real(WP) :: y(this%gas%ns)
      integer  :: iVQ
      ! Mechanical relaxation first
      call this%p_relax(dt,VF,Q,Pjump)
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) return
      ! Frozen gas composition -> mixture-effective ideal-gas parameters
      iVQ=7+this%liq%ns+this%indV-1
      Yv=Q(iVQ)/Q(2)
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cv2=sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG=sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      g2 =cpG/cv2
      ! Liquid shorthands (b via typed nasg pointer)
      cv1=this%liq%cv; g1=this%liq%gamma; pinf=this%liq%pinf
      R1=cv1*(g1-1.0_WP); R2=cv2*(g2-1.0_WP)
      bL=this%liq_nasg%b
      ombm=1.0_WP-bL*Q(1)             ! (1 - m1*b) co-volume factor
      if (ombm.le.0.0_WP) return      ! liquid past co-volume packing limit -> do nothing
      ! Thermal internal energy (formation energies removed); invariant under thermal relax
      Eth=Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*qG
      ! Quadratic for liquid equilibrium pressure Peq (gas pressure Peq-Pjump), TL=TG, pinf_g=0
      a=ombm*(Q(1)*cv1+Q(2)*cv2)
      b=ombm*(Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf-Pjump*(Q(1)*cv1+Q(2)*cv2))-Eth*(Q(1)*R1+Q(2)*R2)
      d=-ombm*Pjump*(Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf)+Pjump*Eth*Q(1)*R1-Eth*Q(2)*R2*pinf
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g=0): liquid p=Peq>-pinf, gas p=Peq-Pjump>0
      if (Peq.le.-this%liq%pinf.or.Peq-Pjump.le.0.0_WP) return
      ! Equilibrium VF (strict bounds); co-volume floor b*Q(1) appears naturally
      VFeq=bL*Q(1)+ombm*Q(1)*R1*(Peq-Pjump)/(Q(1)*R1*(Peq-Pjump)+Q(2)*R2*(Peq+pinf))
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) return
      ! Update Q with the new equilibrium state (co-volume in liquid rhoe, formation energies re-added)
      Q(3)=(VFeq-bL*Q(1))*(Peq+g1*pinf)/(g1-1.0_WP)+Q(1)*this%liq%q
      Q(4)=(1.0_WP-VFeq)*(Peq-Pjump  )/(g2-1.0_WP)+Q(2)*qG
      VF=VFeq
   end subroutine pT_relax

   !> NASG-form energy-conserving equilibrium pressure
   real(WP) function get_p_eq(this,VF_,Q0_,qG_,gammaG_) result(p_eq)
      implicit none
      class(relax_igmix_nasg), intent(in) :: this
      real(WP),                intent(in) :: VF_
      real(WP), dimension(:),  intent(in) :: Q0_
      real(WP),                intent(in) :: qG_,gammaG_
      real(WP) :: one_brho
      one_brho=1.0_WP-this%liq_nasg%b*Q0_(1)/VF_
      p_eq=(sum(Q0_(3:4))-Q0_(1)*this%liq%q-one_brho*VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0_(2)*qG_)/&
      &    (one_brho*VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG_-1.0_WP))
   end function get_p_eq

   !> Equilibrium T from energy conservation (NASG form: includes co-volume b correction)
   real(WP) function get_T_lvg(this,p_,Yv_,rho0,rhoA0) result(T)
      implicit none
      class(relax_igmix_nasg), intent(in) :: this
      real(WP),                intent(in) :: p_,Yv_,rho0,rhoA0
      T=(1.0_WP-Yv_-this%liq_nasg%b*(rho0*(1.0_WP-Yv_)-rhoA0))/                                                       &
      & ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)                          +  &
      &   rhoA0*((this%gas%gamma(this%indV)-1.0_WP)*this%gas%cv(this%indV)*Yv_                                     +  &
      &          (this%gas%gamma(this%indA)-1.0_WP)*this%gas%cv(this%indA)*(1.0_WP-Yv_))/p_)
   end function get_T_lvg

   !> Quadratic coefficients for the equilibrium-T equation (NASG form: includes co-volume terms)
   subroutine get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
      implicit none
      class(relax_igmix_nasg), intent(in)  :: this
      real(WP),                intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP),                intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      real(WP) :: cvV_,gammaV_,qV_
      cvV_   =this%gas%cv   (this%indV)
      gammaV_=this%gas%gamma(this%indV)
      qV_    =this%gas%q    (this%indV)
      ap=rho0*this%liq%cv*cvV_*((gammaV_-this%liq%gamma)*p_eq+this%liq%gamma*(gammaV_-1.0_WP)*this%liq%pinf)
      bp=(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq**2                                                    +&
      &  (this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)                                                          -&
      &   this%liq%gamma*this%liq%cv)+rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+&
      &   rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_))*p_eq                                    +&
      &   (gammaV_-1.0_WP)*cvV_*this%liq%pinf*(rho0*this%liq%q-rhoe0)
      dp=p_eq*(p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
      dapdp=rho0*this%liq%cv*cvV_*(gammaV_-this%liq%gamma)
      dbpdp=2.0_WP*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq                                             +&
      &     this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%gamma*this%liq%cv)                            +&
      &     rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)                          +&
      &     rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_)
      ddpdp=(2.0_WP*p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
   end subroutine get_coeffs_lv

end module relax_igmix_nasg_class
