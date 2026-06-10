!> Thermodynamic relaxation for ideal-gas / stiffened-gas pair
!> Provides p_relax (mechanical) and pT_relax (mechanical+thermal) steps
module relax_ig_sg_class
   use precision,           only: WP
   use thermorelax_class,   only: thermorelax,RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_VACUUM_GAS,RELAX_DEGENERATE
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
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use messager, only: die
      implicit none
      class(relax_ig_sg),     intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      ! Mixture cells only
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Dispatch
      select case (this%model)
      case (Prelax);  call this%p_relax (dt,VF,Q,Pjump,ierr)
      case (PTrelax); call this%pT_relax(dt,VF,Q,Pjump,ierr)
      case default; call die('[relax_ig_sg apply] unknown model')
      end select
   end subroutine apply

   !> Mechanical relaxation: solve quadratic for the liquid equilibrium pressure Peq
   !> (gas pressure Peq-Pjump, i.e. PL-PG=Pjump), then adjust VF and phasic internal
   !> energies via p*dV work exchange at the relaxed interface pressure.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7).
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_ig_sg),     intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.this%RHOGmin) then; if (present(ierr)) ierr=RELAX_VACUUM_GAS; return; end if
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=Q(1)/(       VF),e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])
      ! Return if phasic pressures are out of bound
      if (PL.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (PG.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Phasic acoustic impedances (rho*c)
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(p=PL,rho=Q(1)/(       VF),y=[1.0_WP])
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(p=PG,rho=Q(2)/(1.0_WP-VF),y=[1.0_WP])
      cJ=ZL/(ZG+ZL)
      ! Interface pressure model
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Quadratic setup (pinf_g = 0)
      n1=VF*this%phist
      n0=VF*(this%phi0*Pint-this%phist*cJ*Pjump)+Q(3)-Q(1)*this%liq%q
      d1=this%phist+1.0_WP/(this%liq%gamma-1.0_WP)
      d0=this%phi0*Pint-this%phist*cJ*Pjump+this%liq%gamma/(this%liq%gamma-1.0_WP)*this%liq%pinf
      a=d1*(1.0_WP/(this%gas%gamma-1.0_WP)+this%phist*VF)+n1*(-1.0_WP/(this%gas%gamma-1.0_WP)-this%phist)
      b=d1*(-Pjump/(this%gas%gamma-1.0_WP)+Q(2)*this%gas%q-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n1*(Pjump/(this%gas%gamma-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)+d0*(1.0_WP/(this%gas%gamma-1.0_WP)+this%phist*VF)+n0*(-1.0_WP/(this%gas%gamma-1.0_WP)-this%phist)
      d=d0*(-Pjump/(this%gas%gamma-1.0_WP)+Q(2)*this%gas%q-Q(4)+VF*(this%phi0*Pint-this%phist*cJ*Pjump))+n0*(Pjump/(this%gas%gamma-1.0_WP)-this%phi0*Pint+this%phist*cJ*Pjump)
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g = 0)
      if (Peq.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Equilibrium VF (strict bounds)
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Update Q with p*dV work at the relaxed interface pressure
      Q(3)=Q(3)-(this%phist*Peq+this%phi0*Pint-this%phist*cJ*Pjump)*(VFeq-VF)
      Q(4)=Q(4)+(this%phist*Peq+this%phi0*Pint-this%phist*cJ*Pjump)*(VFeq-VF)
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   end subroutine p_relax

   !> Thermal relaxation: cascades into p_relax first, then enforces TL=TG with
   !> mechanical equilibrium PL=Peq, PG=Peq-Pjump (PL-PG=Pjump). Closed-form quadratic.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7).
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_ig_sg),     intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: a,b,d,Peq,VFeq,Eth
      real(WP) :: cv1,cv2,g1,g2,pinf,R1,R2
      ! Mechanical relaxation first (its own ierr is not propagated; pT owns the final verdict)
      call this%p_relax(dt,VF,Q,Pjump)
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Shorthands
      cv1=this%liq%cv; cv2=this%gas%cv
      g1=this%liq%gamma; g2=this%gas%gamma; pinf=this%liq%pinf
      R1=cv1*(g1-1.0_WP); R2=cv2*(g2-1.0_WP)
      ! Thermal internal energy (formation energies removed); invariant under thermal relax
      Eth=Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*this%gas%q
      ! Quadratic for liquid equilibrium pressure Peq (gas pressure Peq-Pjump), TL=TG, pinf_g=0
      a=Q(1)*cv1+Q(2)*cv2
      b=Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf-Pjump*(Q(1)*cv1+Q(2)*cv2)-Eth*(Q(1)*R1+Q(2)*R2)
      d=-Pjump*(Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf)+Pjump*Eth*Q(1)*R1-Eth*Q(2)*R2*pinf
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g=0): liquid p=Peq>-pinf, gas p=Peq-Pjump>0
      if (Peq.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Equilibrium VF (strict bounds)
      VFeq=Q(1)*R1*(Peq-Pjump)/(Q(1)*R1*(Peq-Pjump)+Q(2)*R2*(Peq+pinf))
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Update Q with the new equilibrium state (formation energies re-added)
      Q(3)=(       VFeq)*(Peq+g1*pinf)/(g1-1.0_WP)+Q(1)*this%liq%q
      Q(4)=(1.0_WP-VFeq)*(Peq-Pjump  )/(g2-1.0_WP)+Q(2)*this%gas%q
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pT_relax

end module relax_ig_sg_class
