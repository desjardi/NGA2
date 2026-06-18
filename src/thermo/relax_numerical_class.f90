!> EOS-agnostic mechanical (pressure) relaxation.
!> Drives the two phases to a common pressure (PL=Peq, PG=Peq-Pjump) by a small 2x2
!> Newton in (VFeq, Peq), exchanging energy through p*dV work at the matm interface
!> pressure. Works for ANY class(material) pair and reduces exactly to the NASG
!> analytical p_relax quadratic. Conserves phasic masses Q(1:2), total internal
!> energy Q(3)+Q(4), momentum Q(5:).
!>
!> The Jacobian is ANALYTIC, built from the EOS sound speed c and Gruneisen Gamma:
!>   (dp/drho)_e = c^2 - p*Gamma/rho ,  (dp/de)_rho = rho*Gamma
!> so per Newton iteration each phase is evaluated at ONE (rho,e) state (get_p +
!> get_c + get_gruneisen) -- no finite-difference perturbations. With an EOS whose
!> evaluation is expensive AND caches by (rho,e) (e.g. a CoolProp-backed material),
!> the three queries hit one flash, so it is ~1 flash/iteration instead of FD's ~3
!> (FD perturbs to 3 distinct states). Iteration count and accuracy are unchanged vs
!> FD; the gain is purely flashes-per-iteration. The achievable rtol floor is set by
!> the EOS, not the Jacobian (e.g. stiffened-gas pinf cancellation limits NASG to
!> ~1e-11 relative; an EOS without pinf does better).
module relax_numerical_class
   use precision,         only: WP
   use material_class,    only: material
   use thermorelax_class, only: thermorelax,RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_VACUUM_GAS,RELAX_DEGENERATE,RELAX_SINGULAR
   implicit none
   private

   public :: prelax_numerical
   public :: relax_numerical

   !> Thermorelax wrapper around prelax_numerical: an EOS-agnostic MECHANICAL (p) relaxation
   !> dispatchable from the amrmpcomp apply path. Holds the two phase materials + the p-relax
   !> options; apply() forwards to the prelax_numerical kernel. Pure phases only (single-component
   !> liquid & gas); pT/pTg relaxation and multicomponent composition extraction are not yet done.
   type, extends(thermorelax) :: relax_numerical
      class(material), pointer :: liq=>null(), gas=>null()
      real(WP) :: RHOGmin=1.0e-2_WP   !< skip near-pure-liquid cells (gas density floor)
      real(WP) :: phist  =1.0_WP      !< interface-pressure weight on Peq
      real(WP) :: phi0   =0.0_WP      !< interface-pressure weight on Pint
      real(WP) :: rtol   =1.0e-6_WP   !< Newton relative tolerance
      integer  :: itmax  =50          !< Newton iteration cap
   contains
      procedure :: initialize
      procedure :: apply
   end type relax_numerical

contains

   !> Store the two phase materials (any class(material))
   subroutine initialize(this,liq,gas)
      implicit none
      class(relax_numerical),  intent(inout) :: this
      class(material), target, intent(in)    :: liq,gas
      this%liq=>liq; this%gas=>gas
   end subroutine initialize

   !> Apply: EOS-agnostic mechanical relaxation of one cell (pure phases -> y=[1]).
   !> The kernel does the mixture-cell / soundness gating and sets ierr.
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_numerical), intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: yL(this%liq%ns),yG(this%gas%ns)
      integer  :: k,Yl_lo,Yg_lo
      ! Phasic compositions from the amrmpcomp Q-layout: liquid species at slots 8.., gas species
      ! at 7+liq%ns.. (Q(slot)=phase_mass*Y_k, last species by closure). Pure phase -> [1].
      yL=0.0_WP; yL(1)=1.0_WP
      yG=0.0_WP; yG(1)=1.0_WP
      if (Q(1).gt.0.0_WP.and.Q(2).gt.0.0_WP) then
         if (this%liq%ns.gt.1) then
            Yl_lo=8
            do k=1,this%liq%ns-1; yL(k)=Q(Yl_lo+k-1)/Q(1); end do
            yL(this%liq%ns)=1.0_WP-sum(yL(1:this%liq%ns-1))
         end if
         if (this%gas%ns.gt.1) then
            Yg_lo=7+this%liq%ns
            do k=1,this%gas%ns-1; yG(k)=Q(Yg_lo+k-1)/Q(2); end do
            yG(this%gas%ns)=1.0_WP-sum(yG(1:this%gas%ns-1))
         end if
      end if
      call prelax_numerical(this%liq,this%gas,yL,yG,VF,Q,Pjump, &
      &    phist=this%phist,phi0=this%phi0,itmax=this%itmax,rtol=this%rtol,RHOGmin=this%RHOGmin,ierr=ierr)
   end subroutine apply

   !> Numerical pressure relaxation of a single mixture cell.
   !>   liq/gas : the two phase materials (any class(material))
   !>   yL/yG   : phasic compositions passed to the EOS calls (pure phase -> [1.0])
   !>   VF      : liquid volume fraction (inout)
   !>   Q       : conserved vector; uses Q(1:4) = [mL, mG, alphaL*rhoL*eL, alphaG*rhoG*eG]
   !>   Pjump   : interface pressure jump PL-PG (surface tension); 0 if none
   !>   phist/phi0 : interface-pressure weights (defaults 1, 0 -> matm/Pelanti)
   !>   RHOGmin : skip near-pure-liquid cells with gas density below this (default 1e-2)
   !>   niter   : (out) Newton iterations taken
   !>   ierr    : (out) relaxation status (see thermorelax_class RELAX_* codes; 0=ok)
   !> Soundness is tested via get_c<=0: the material EOS clamp returns a non-positive sound
   !> speed for an unphysical (rho,e), so cL<=0 flags a bad liquid (co-volume / sub-vacuum)
   !> and cG<=0 a bad gas. The cell is left untouched (not committed) on any non-ok status.
   subroutine prelax_numerical(liq,gas,yL,yG,VF,Q,Pjump,phist,phi0,itmax,rtol,RHOGmin,niter,ierr)
      implicit none
      class(material),        intent(in)    :: liq,gas
      real(WP), dimension(:), intent(in)    :: yL,yG
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP), optional,     intent(in)    :: phist,phi0,rtol,RHOGmin
      integer,  optional,     intent(in)    :: itmax
      integer,  optional,     intent(out)   :: niter,ierr
      real(WP) :: ph,p0,tol,rgmin,VF0
      real(WP) :: rhoL,rhoG,PL,PG,cL0,cG0,ZL,ZG,Pint,cJ,pI
      real(WP) :: x(2),R(2),Jm(2,2),detJ,VFeq,Peq
      real(WP) :: rhoLf,rhoGf,eLf,eGf,cLf,cGf
      integer  :: nit,it,ni
      ! Options (defaults reproduce the deployed NASG p_relax; rtol 1e-6 -> ~2 iters)
      ph =1.0_WP;     if (present(phist)) ph =phist
      p0 =0.0_WP;     if (present(phi0 )) p0 =phi0
      tol=1.0e-6_WP;  if (present(rtol )) tol=rtol
      nit=50;         if (present(itmax)) nit=itmax
      rgmin=1.0e-2_WP;if (present(RHOGmin)) rgmin=RHOGmin
      if (present(niter)) niter=0
      ! Skip degenerate cells
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      VF0=VF
      ! Initial phasic state; acoustic impedances and interface pressure are frozen during the relax
      rhoL=Q(1)/VF; rhoG=Q(2)/(1.0_WP-VF)
      ! Skip near-pure-liquid cells (gas density too low)
      if (rhoG.lt.rgmin) then; if (present(ierr)) ierr=RELAX_VACUUM_GAS; return; end if
      PL=liq%get_p_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=yL)
      PG=gas%get_p_from_rho_e(rho=rhoG,e=Q(4)/Q(2),y=yG)
      cL0=liq%get_c_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=yL)
      cG0=gas%get_c_from_rho_e(rho=rhoG,e=Q(4)/Q(2),y=yG)
      ! Initial-state soundness (clamped get_c returns 0 for an unphysical state)
      if (cL0.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (cG0.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS;    return; end if
      ZL=rhoL*cL0; ZG=rhoG*cG0
      cJ  =ZL/(ZG+ZL)
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! 2x2 Newton in x = (VFeq, Peq), analytic Jacobian; seed at current VF and liquid pressure
      x=[VF,PL]; ni=0
      do it=1,nit
         call resjac(x,R,Jm)
         if (maxval(abs(R)).lt.tol*max(abs(x(2)),1.0_WP)) exit
         ni=ni+1
         detJ=Jm(1,1)*Jm(2,2)-Jm(1,2)*Jm(2,1)
         if (abs(detJ).lt.tiny(1.0_WP)) then; if (present(ierr)) ierr=RELAX_SINGULAR; return; end if
         x=x-[ Jm(2,2)*R(1)-Jm(1,2)*R(2), -Jm(2,1)*R(1)+Jm(1,1)*R(2)]/detJ
         if (x(1).le.0.0_WP.or.x(1).ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if   ! VFeq left the physical interval
      end do
      if (present(niter)) niter=ni
      if (maxval(abs(R)).ge.tol*max(abs(x(2)),1.0_WP)) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if   ! no convergence
      VFeq=x(1); Peq=x(2)
      if (VFeq.le.0.0_WP.or.VFeq.ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Converged-state soundness: re-evaluate both phases at the committed (rho,e); leave cell if unphysical
      pI=ph*Peq+p0*Pint-ph*cJ*Pjump
      rhoLf=Q(1)/VFeq; rhoGf=Q(2)/(1.0_WP-VFeq)
      eLf=(Q(3)-pI*(VFeq-VF0))/Q(1); eGf=(Q(4)+pI*(VFeq-VF0))/Q(2)
      cLf=liq%get_c_from_rho_e(rho=rhoLf,e=eLf,y=yL)
      cGf=gas%get_c_from_rho_e(rho=rhoGf,e=eGf,y=yG)
      if (cLf.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (cGf.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS;    return; end if
      ! Commit: p*dV work at the relaxed interface pressure
      Q(3)=Q(3)-pI*(VFeq-VF0)
      Q(4)=Q(4)+pI*(VFeq-VF0)
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   contains
      !> Residuals (liquid at Peq, gas at Peq-Pjump) and the analytic 2x2 Jacobian.
      !> Each phase is evaluated at ONE (rho,e): get_p, get_c, get_gruneisen.
      subroutine resjac(xx,RR,JJ)
         real(WP), intent(in)  :: xx(2)
         real(WP), intent(out) :: RR(2),JJ(2,2)
         real(WP) :: s,pe,pIl,dWp,rL,rG,eL,eG,PLv,PGv,cL,cG,GL,GG,pLr,pLe,pGr,pGe
         s=xx(1); pe=xx(2)
         pIl=ph*pe+p0*Pint-ph*cJ*Pjump            ! interface pressure (linear in Peq)
         dWp=ph*(s-VF0)                            ! d(work)/d(Peq)   ; d(work)/ds = pIl
         rL =Q(1)/s;                rG =Q(2)/(1.0_WP-s)
         eL =(Q(3)-pIl*(s-VF0))/Q(1);  eG =(Q(4)+pIl*(s-VF0))/Q(2)
         PLv=liq%get_p_from_rho_e(rho=rL,e=eL,y=yL)
         PGv=gas%get_p_from_rho_e(rho=rG,e=eG,y=yG)
         cL =liq%get_c_from_rho_e(rho=rL,e=eL,y=yL)
         cG =gas%get_c_from_rho_e(rho=rG,e=eG,y=yG)
         GL =liq%get_gruneisen_from_rho_e(rho=rL,e=eL,y=yL)
         GG =gas%get_gruneisen_from_rho_e(rho=rG,e=eG,y=yG)
         pLe=rL*GL; pLr=cL*cL-PLv*GL/rL           ! (dp/de)_rho , (dp/drho)_e
         pGe=rG*GG; pGr=cG*cG-PGv*GG/rG
         RR(1)=PLv-pe
         RR(2)=PGv-(pe-Pjump)
         JJ(1,1)=pLr*(-rL/s)            + pLe*(-pIl/Q(1))      ! dR1/ds
         JJ(1,2)=pLe*(-dWp/Q(1)) - 1.0_WP                     ! dR1/dPeq
         JJ(2,1)=pGr*( rG/(1.0_WP-s))  + pGe*( pIl/Q(2))      ! dR2/ds
         JJ(2,2)=pGe*( dWp/Q(2)) - 1.0_WP                     ! dR2/dPeq
      end subroutine resjac
   end subroutine prelax_numerical

end module relax_numerical_class
