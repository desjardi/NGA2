!> Thermodynamic relaxation for ideal-gas / stiffened-gas pair
!> Provides p_relax (mechanical) and pT_relax (mechanical+thermal) steps
module relax_ig_sg_class
   use precision,           only: WP
   use thermorelax_class,   only: thermorelax
   use ideal_gas_class,     only: ideal_gas
   use stiffened_gas_class, only: stiffened_gas
   implicit none
   private

   public :: relax_ig_sg
   public :: Prelax,PTrelax

   integer, parameter :: Prelax =1    !< Mechanical relaxation only
   integer, parameter :: PTrelax=2    !< Mechanical + thermal relaxation

   type, extends(thermorelax) :: relax_ig_sg
      class(ideal_gas),     pointer :: gas=>null()
      class(stiffened_gas), pointer :: liq=>null()
      integer  :: model  =Prelax
      real(WP) :: RHOGmin=1.0e-2_WP   !< Skip mechanical relax when gas density falls below this
      real(WP) :: phist  =1.0_WP      !< Temporal weighting on equilibrium pressure
      real(WP) :: phi0   =0.0_WP      !< Temporal weighting on interface pressure
   contains
      procedure :: initialize
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
   end type relax_ig_sg

contains

   !> Initialize: store EOS pointers
   subroutine initialize(this,gas,liq)
      implicit none
      class(relax_ig_sg),           intent(inout) :: this
      class(ideal_gas),     target, intent(in)    :: gas
      class(stiffened_gas), target, intent(in)    :: liq
      this%gas=>gas; this%liq=>liq
   end subroutine initialize

   !> Apply: gate on mixture cells, dispatch via model
   subroutine apply(this,dt,VF,Q,Pjump)
      use messager, only: die
      implicit none
      class(relax_ig_sg),     intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      ! Mixture cells only
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP) return
      ! Dispatch
      select case (this%model)
      case (Prelax);  call this%p_relax (dt,VF,Q,Pjump)
      case (PTrelax); call this%pT_relax(dt,VF,Q,Pjump)
      case default; call die('[relax_ig_sg apply] unknown model')
      end select
   end subroutine apply

   !> Mechanical relaxation: solve quadratic for Peq where PL+Pjump=PG=Peq,
   !> then adjust VF and phasic internal energies via p*dV work exchange.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7).
   subroutine p_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_ig_sg),     intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.this%RHOGmin) return
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=Q(1)/(       VF),e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])
      ! Phasic acoustic impedances
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(p=PL,rho=Q(1)/(       VF),y=[1.0_WP])**2
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(p=PG,rho=Q(2)/(1.0_WP-VF),y=[1.0_WP])**2
      cJ=ZL/(ZG+ZL)
      ! Interface pressure model
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Quadratic setup (pinf_g = 0)
      n1=VF*this%phist
      n0=VF*(this%phi0*Pint-this%phist*cJ*Pjump)+Q(3)
      d1=this%phist+1.0_WP/(this%liq%gamma-1.0_WP)
      d0=this%phi0*Pint-this%phist*cJ*Pjump+this%liq%gamma/(this%liq%gamma-1.0_WP)*this%liq%pinf
      a=d1*(1.0_WP/(this%gas%gamma-1.0_WP)+this%phist*VF)+n1*(-1.0_WP/(this%gas%gamma-1.0_WP)-this%phist)
      b=d1*(-Pjump/(this%gas%gamma-1.0_WP)-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n1*(Pjump/(this%gas%gamma-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)+d0*(1.0_WP/(this%gas%gamma-1.0_WP)+this%phist*VF)+n0*(-1.0_WP/(this%gas%gamma-1.0_WP)-this%phist)
      d=d0*(-Pjump/(this%gas%gamma-1.0_WP)-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n0*(Pjump/(this%gas%gamma-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g = 0)
      if (Peq.le.-this%liq%pinf.or.Peq-Pjump.le.0.0_WP) return
      ! Equilibrium VF (strict bounds)
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) return
      ! Update Q with p*dV work
      Q(3)=Q(3)-(this%phi0*Pint+this%phist*Peq)*(VFeq-VF)
      Q(4)=Q(4)+(this%phi0*Pint+this%phist*Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine p_relax

   !> Thermal relaxation: cascades into p_relax first, then enforces TL=TG.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7).
   !> PJUMP IS NOT PROPERLY ACCOUNTED FOR HERE!
   subroutine pT_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_ig_sg),     intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP) :: a,b,d,Peq,VFeq,Etot
      ! Mechanical relaxation first
      call this%p_relax(dt,VF,Q,Pjump)
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) return
      ! Total internal energy is invariant under thermal relax
      Etot=Q(3)+Q(4)
      ! Quadratic setup (pinf_g = 0)
      a=Q(1)*this%liq%cv+Q(2)*this%gas%cv
      b=Q(1)*this%liq%cv*this%liq%gamma*this%liq%pinf+Q(2)*this%gas%cv*this%liq%pinf &
      & -Etot*(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)+Q(2)*this%gas%cv*(this%gas%gamma-1.0_WP))
      d=-Etot*Q(2)*this%gas%cv*(this%gas%gamma-1.0_WP)*this%liq%pinf
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g = 0)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) return
      ! Equilibrium VF (strict bounds)
      VFeq=Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq &
      &   /(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq+Q(2)*this%gas%cv*(this%gas%gamma-1.0_WP)*(Peq+this%liq%pinf))
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) return
      ! Update Q with the new equilibrium state
      Q(3)=(       VFeq)*(Peq+this%liq%gamma*this%liq%pinf)/(this%liq%gamma-1.0_WP)
      Q(4)=(1.0_WP-VFeq)*Peq/(this%gas%gamma-1.0_WP)
      VF=VFeq
   end subroutine pT_relax

end module relax_ig_sg_class
