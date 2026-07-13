!> Case-specific "safe" relaxation model for the amrcomp_impact problem.
!> Extends the stock relax_ig_nasg quadratic relaxation with the case's single-cell
!> equilibration policy. Per mixture cell, apply() runs, in order:
!>  1) absorb (seed culling): sub-resolution gas at the high-pressure extreme is absorbed
!>     into the liquid rather than sustained — a gas packet above diss_P is supercritical
!>     and mixes into the liquid (conserves cell totals exactly; the cell becomes pure
!>     liquid and the solver's pure-cell snap completes the PLIC reset). The low-pressure
!>     counterpart (vacuum kill) was removed 2026-07-13: measured not load-bearing
!>     (swap+floor+pool cover its population), and the near-vacuum seeds it culled are
!>     physical-topology candidates (pocket growth) that phase change should adjudicate.
!>  2) propose VF: the stock quadratic relaxation (Prelax/PTrelax/PThybrid dispatch,
!>     VFratmax intact); its return code only feeds the ledger
!>  3) set energies: the fixed-VF energy swap runs UNCONDITIONALLY — whatever VF step 2
!>     produced (converged, VFratmax-clamped, or unchanged on failure), the unique energy
!>     split with PL-PG=Pjump is imposed at frozen VF and masses (conserves phasic masses,
!>     total energy, momentum; a no-op to roundoff where the quadratic fully converged;
!>     completes the equilibration where it was clamped or failed; also the vacuum-runaway
!>     cutoff). Only a non-positive phase mass exits untouched (counted as stuck).
!>  4) floor: if the shared pressure sits below any user limit (PL,TL,PG,TG), it is raised
!>     to the binding limit before the energies are written (all four limits rise
!>     monotonically with the shared pressure, so the target is a direct computation)
!> The per-phase (Pmin,Tmin) corner rescue lives solver-side in clean_Q (EOS-generic,
!> all levels); the same limits are duplicated here for the floor, and the case sets
!> both pairs together. Every action is ledgered in acc (rank-local cumulative; dm/dE
!> scaled by the vol member, valid because apply() only runs on the finest level);
!> the case reduces acc across ranks for its own monitoring.
module safe_relax_class
   use precision,           only: WP
   use thermorelax_class,   only: RELAX_OK,RELAX_DEGENERATE
   use relax_ig_nasg_class, only: relax_ig_nasg,Prelax,PTrelax,PThybrid
   use amrvof_class,        only: VFlo,VFhi
   implicit none
   private

   public :: safe_relax

   type, extends(relax_ig_nasg) :: safe_relax
      ! Floor limits: minimal (P,T) per phase at the post-relaxation equilibrium
      ! (same values as the solver's clean_Q rescue limits; defaults never fire)
      real(WP) :: Pmin_liq=-1.0e30_WP   !< Liquid pressure limit (e.g. max sustainable tension)
      real(WP) :: Tmin_liq=-1.0_WP      !< Liquid temperature limit
      real(WP) :: Pmin_gas=-1.0e30_WP   !< Gas pressure limit (e.g. ~saturation pressure)
      real(WP) :: Tmin_gas=-1.0_WP      !< Gas temperature limit
      ! Absorb threshold (default off)
      real(WP) :: diss_P=1.0e30_WP      !< Dissolution: absorb gas where phasic gas pressure exceeds this (1e30=off)
      ! Cell volume for ledger units (set by the case; apply() runs on the finest level only)
      real(WP) :: vol=1.0_WP
      ! Ledger: rank-local cumulative accumulators (1-2 dissolution n/dm; 3 quadratic proposal
      ! succeeded; 4 proposal failed, completed by the swap alone; 5-6 floor n/dE; 7 untouched
      ! cells, non-positive phase mass). The model only counts; reduction across ranks is the
      ! monitoring code's business (the case reduces acc itself)
      real(WP), dimension(7) :: acc=0.0_WP
   contains
      procedure :: apply
   end type safe_relax

contains

   !> Apply the safe relaxation policy: absorb -> equilibrate -> floor
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(safe_relax),      intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: TL,TG,ombm,Asum,Eth,PL,Ptar
      integer :: ier
      ! Only mixture cells; the solver's pure-cell snap owns the rest
      if (VF.lt.VFlo.or.VF.gt.VFhi) then
         if (present(ierr)) ierr=RELAX_DEGENERATE
         return
      end if
      ! Absorb (seed culling, high-pressure extreme): a gas packet above diss_P is
      ! supercritical and mixes into the liquid; conserves cell totals exactly, the cell
      ! becomes pure liquid, and the caller's pure-cell snap completes the PLIC reset
      if (this%diss_P.lt.1.0e30_WP.and.Q(2).gt.0.0_WP.and.Q(4).gt.0.0_WP) then
         if (this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP]).gt.this%diss_P) then
            this%acc(1)=this%acc(1)+1.0_WP; this%acc(2)=this%acc(2)+Q(2)*this%vol
            Q(1)=Q(1)+Q(2); Q(3)=Q(3)+Q(4)
            Q(2)=0.0_WP; Q(4)=0.0_WP
            VF=1.0_WP
            if (present(ierr)) ierr=RELAX_OK
            return
         end if
      end if
      ! Stage 1 — propose VF: stock quadratic relaxation (Prelax/PTrelax/PThybrid dispatch).
      ! Its return code only feeds the ledger; the energy split is set unconditionally below.
      ier=RELAX_OK
      select case (this%model)
      case (Prelax);  call this%p_relax (dt,VF,Q,Pjump,ier)
      case (PTrelax); call this%pT_relax(dt,VF,Q,Pjump,ier)
      case (PThybrid)
         if (any(Q(1:4).le.0.0_WP)) then
            call this%p_relax(dt,VF,Q,Pjump,ier)
         else
            TL=this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
            TG=this%gas%get_T_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])
            if (TL.gt.0.0_WP.and.TG.gt.0.0_WP.and.max(TG/TL,TL/TG).gt.this%Tratmax) then
               call this%pT_relax(dt,VF,Q,Pjump,ier)
            else
               call this%p_relax(dt,VF,Q,Pjump,ier)
            end if
         end if
      end select
      ! Ledger the proposal outcome; a non-positive phase mass is the only untouched exit
      if (ier.eq.RELAX_OK) then
         this%acc(3)=this%acc(3)+1.0_WP
      else if (Q(1).le.0.0_WP.or.Q(2).le.0.0_WP) then
         this%acc(7)=this%acc(7)+1.0_WP
         if (present(ierr)) ierr=ier
         return
      else
         this%acc(4)=this%acc(4)+1.0_WP
      end if
      ! The proposal may exit at the VF bounds; the solver's pure-cell snap owns those
      if (VF.lt.VFlo.or.VF.gt.VFhi) then
         if (present(ierr)) ierr=RELAX_OK
         return
      end if
      ! Stage 2 — set energies (unconditional fixed-VF swap): whatever VF stage 1 produced
      ! (converged, VFratmax-clamped, or unchanged on failure), impose the unique energy
      ! split with PL-PG=Pjump at frozen VF and masses. Conserves phasic masses, total
      ! energy, momentum; a no-op to roundoff where the quadratic fully converged; completes
      ! the equilibration where it was clamped or failed (also the vacuum-runaway cutoff).
      ! Co-volume factor VF*(1-b*rhoL) clamped consistently with the nasg accessors.
      ombm=VF*max(1.0_WP-this%liq%b*Q(1)/VF,1.0_WP-this%liq%brhomax)
      Asum=ombm/(this%liq%gamma-1.0_WP)+(1.0_WP-VF)/(this%gas%gamma-1.0_WP)
      Eth=Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*this%gas%q
      PL=(Eth-ombm*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)+(1.0_WP-VF)*Pjump/(this%gas%gamma-1.0_WP))/Asum
      ! Stage 3 — floor: if the shared pressure sits below any user limit (PL,TL,PG,TG),
      ! raise it to the binding limit (all four rise monotonically with the shared pressure);
      ! the energy added is Asum*(Ptar-PL). Ledgered.
      Ptar=-1.0e30_WP
      if (this%Pmin_liq.gt.-1.0e29_WP) Ptar=max(Ptar,this%Pmin_liq)
      if (this%Pmin_gas.gt.-1.0e29_WP) Ptar=max(Ptar,this%Pmin_gas+Pjump)
      if (this%Tmin_liq.gt.0.0_WP) Ptar=max(Ptar,this%liq%get_p_from_rho_T(rho=Q(1)/VF,T=this%Tmin_liq,y=[1.0_WP]))
      if (this%Tmin_gas.gt.0.0_WP) Ptar=max(Ptar,this%gas%get_p_from_rho_T(rho=Q(2)/(1.0_WP-VF),T=this%Tmin_gas,y=[1.0_WP])+Pjump)
      if (PL.lt.Ptar) then
         this%acc(5)=this%acc(5)+1.0_WP
         this%acc(6)=this%acc(6)+Asum*(Ptar-PL)*this%vol
         PL=Ptar
      end if
      ! Write the split
      Q(3)=ombm*(PL+this%liq%gamma*this%liq%pinf)/(this%liq%gamma-1.0_WP)+Q(1)*this%liq%q
      Q(4)=(1.0_WP-VF)*(PL-Pjump)/(this%gas%gamma-1.0_WP)+Q(2)*this%gas%q
      if (present(ierr)) ierr=RELAX_OK
   end subroutine apply

end module safe_relax_class
