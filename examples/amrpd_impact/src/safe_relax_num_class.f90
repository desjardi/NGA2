!> Case-specific "safe" relaxation model for the amrcomp_impact problem, EOS-AGNOSTIC:
!> the numerical mechanical relaxation (relax_numerical) plus the case's single-cell
!> equilibration policy, working for ANY class(material) pair (pure phases). Replaces
!> the NASG-analytic safe_relax. Thermal relaxation (pT/PThybrid) is deliberately
!> dropped: measured near-inert on this case (relax_therm census 2026-06-24, Peq_therm
!> ~= Peq_mech with ~no gas-ward energy motion), and interphase heat exchange is
!> carried physically by the solver's phasic conduction.
!> Per mixture cell, apply() runs, in order:
!>  1) absorb (seed culling): sub-resolution gas at the high-pressure extreme is
!>     absorbed into the liquid -- a gas packet above diss_P is supercritical and
!>     mixes into the liquid (conserves cell totals exactly; the cell becomes pure
!>     liquid and the solver's pure-cell snap completes the PLIC reset)
!>  2) propose VF: FULL pT relaxation (pTrelax_numerical: PL-PG=Pjump and TL=TG) where
!>     the phasic temperature contrast exceeds Tratmax -- the conservative in-cell
!>     quench of superheated sub-resolution wisps -- and the parent EOS-agnostic
!>     Newton p-relax (relax_numerical%apply) otherwise or as fallback; return codes
!>     only feed the ledger. The per-call phase-volume change is then clamped to a
!>     factor VFratmax -- post-hoc clamping is consistent here because stage 3
!>     re-splits the conserved energy at the final VF regardless.
!>  3) set energies: the fixed-VF energy split runs UNCONDITIONALLY -- whatever VF
!>     stage 2 produced (converged, clamped, or unchanged on failure), the unique
!>     energy split with PL-PG=Pjump is imposed at frozen VF and masses by solving
!>       VF*rhoe_L(PL,rhoL) + (1-VF)*rhoe_G(PL-Pjump,rhoG) = Etot
!>     for PL by Newton with the analytic slope Asum = VF/GammaL + (1-VF)/GammaG
!>     from get_gruneisen. For any Gruneisen-form EOS (ideal gas, SG, NASG, MG --
!>     p affine in e at fixed rho) one step is EXACT; otherwise it is a short
!>     fixed-point. Conserves phasic masses, total energy, momentum; a no-op to
!>     roundoff where the proposal fully converged; completes the equilibration
!>     where it was clamped or failed (also the vacuum-runaway cutoff). Only a
!>     non-positive phase mass exits untouched (counted as stuck).
!>  4) floor: if the shared pressure sits below any user limit (PL,TL,PG,TG), it is
!>     raised to the binding limit before the energies are written (all four limits
!>     rise monotonically with the shared pressure); the energy added is
!>     Asum*(Ptar-PL), ledgered.
!> The per-phase (Pmin,Tmin) corner rescue lives solver-side in clean_Q (EOS-generic,
!> all levels); the same limits are duplicated here for the floor, and the case sets
!> both pairs together. Ledger acc(1:7) keeps the safe_relax slot meanings (1-2
!> dissolution n/dm; 3 proposal succeeded; 4 proposal failed, completed by the split;
!> 5-6 floor n/dE; 7 untouched, non-positive phase mass), so the case's rescue
!> monitor works unchanged.
module safe_relax_num_class
   use precision,             only: WP
   use material_class,        only: material
   use thermorelax_class,     only: RELAX_OK,RELAX_DEGENERATE
   use relax_numerical_class, only: relax_numerical,pTrelax_numerical
   use amrvof_class,          only: VFlo,VFhi
   implicit none
   private

   public :: safe_relax_num

   type, extends(relax_numerical) :: safe_relax_num
      ! Floor limits: minimal (P,T) per phase at the post-relaxation equilibrium
      ! (same values as the solver's clean_Q rescue limits; defaults never fire)
      real(WP) :: Pmin_liq=-1.0e30_WP   !< Liquid pressure limit (e.g. max sustainable tension)
      real(WP) :: Tmin_liq=-1.0_WP      !< Liquid temperature limit
      real(WP) :: Pmin_gas=-1.0e30_WP   !< Gas pressure limit (e.g. ~saturation pressure)
      real(WP) :: Tmin_gas=-1.0_WP      !< Gas temperature limit
      ! Absorb threshold (default off)
      real(WP) :: diss_P=1.0e30_WP      !< Dissolution: absorb gas where phasic gas pressure exceeds this (1e30=off)
      real(WP) :: diss_VFmin=0.9_WP     !< Absorb only LIQUID-DOMINANT cells (VF>=this). Absorbing a gas-dominant cell
                                        !< dilutes its liquid below the corner density and the clean_Q guard re-masses it:
                                        !< a mass-fabrication feedback (measured 2026-07-16: LiqMass x2.6, dt 5e-6, blowup)
      ! pT-hybrid dispatch
      real(WP) :: Tratmax=10.0_WP       !< Phasic temperature contrast max(TG/TL,TL/TG) above which the proposal is the
                                        !< FULL pT relaxation (TL=TG quench of sub-resolution wisps in intimate thermal
                                        !< contact -- conservative, in-cell; the PThybrid role the relax_therm census
                                        !< never measured: pT leaves Peq unchanged but caps the temperature split);
                                        !< 1e30 = mechanical-only
      real(WP) :: Tliftmax=1.5_WP       !< Direction gate: run pT only if the capacity-weighted equilibrium estimate
                                        !< T* = (mL*cvL*TL+mG*cvG*TG)/(mL*cvL+mG*cvG) stays below Tliftmax*TL -- pT must
                                        !< QUENCH the hot minority toward the liquid, never drag minority liquid to the
                                        !< gas temperature (an MG liquid at gas T equilibrates by flash expansion to
                                        !< rho~14-26: unphysical hot "liquid" that seeds dP hotspots; that regime is
                                        !< really evaporation and belongs to pTg phase change)
      ! Proposal clamp
      real(WP) :: VFratmax=10.0_WP      !< Max per-call phase-volume change factor (post-hoc clamp on stage 2)
      ! Cell volume for ledger units (set by the case; apply() runs on the finest level only)
      real(WP) :: vol=1.0_WP
      ! Ledger: rank-local cumulative accumulators (see header; the case reduces acc itself)
      real(WP), dimension(7) :: acc=0.0_WP
   contains
      procedure :: initialize
      procedure :: apply
   end type safe_relax_num

contains

   !> Initialize: parent stores the material pointers; the ladder is pure-phase only
   subroutine initialize(this,liq,gas)
      use messager, only: die
      implicit none
      class(safe_relax_num),   intent(inout) :: this
      class(material), target, intent(in)    :: liq,gas
      call this%relax_numerical%initialize(liq=liq,gas=gas)
      if (liq%ns.ne.1.or.gas%ns.ne.1) call die('[safe_relax_num initialize] pure phases only (ns=1): the ladder passes y=[1]')
   end subroutine initialize

   !> Apply the safe relaxation policy: absorb -> propose -> split -> floor
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(safe_relax_num),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP) :: VF0,rhoL,rhoG,Etot,PL,Ptar,QL,QG,GL,GG,Asum,dPL,TL,TG,cvL,cvG,Tstar
      logical :: didpT
      integer :: ier,k
      ! Only mixture cells; the solver's pure-cell snap owns the rest
      if (VF.lt.VFlo.or.VF.gt.VFhi) then
         if (present(ierr)) ierr=RELAX_DEGENERATE
         return
      end if
      ! Absorb (seed culling, high-pressure extreme): a gas packet above diss_P is
      ! supercritical and mixes into the liquid; conserves cell totals exactly, the cell
      ! becomes pure liquid, and the caller's pure-cell snap completes the PLIC reset.
      ! LIQUID-DOMINANT cells only (VF>=diss_VFmin): absorbing a gas-dominant cell dilutes
      ! its liquid to rhoL*VF and the corner guard re-masses it -- a fabrication feedback.
      if (this%diss_P.lt.1.0e30_WP.and.VF.ge.this%diss_VFmin &
      &   .and.Q(2).gt.0.0_WP.and.Q(4).gt.0.0_WP) then
         if (this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP]).gt.this%diss_P) then
            this%acc(1)=this%acc(1)+1.0_WP; this%acc(2)=this%acc(2)+Q(2)*this%vol
            Q(1)=Q(1)+Q(2); Q(3)=Q(3)+Q(4)
            Q(2)=0.0_WP; Q(4)=0.0_WP
            VF=1.0_WP
            if (present(ierr)) ierr=RELAX_OK
            return
         end if
      end if
      ! Stage 1 -- propose VF: FULL pT relaxation where the phasic temperature contrast
      ! exceeds Tratmax (instantaneous TL=TG quench of sub-resolution wisps in intimate
      ! thermal contact -- conservative and in-cell, structurally incapable of the
      ! absorb/guard fabrication feedback), mechanical-only p-relax otherwise, with
      ! p-relax as the fallback when pT fails. Return codes only feed the ledger; the
      ! per-call phase-volume change is then clamped to a factor VFratmax (bounds
      ! consistent for VFratmax>=1; the split below re-derives the energies)
      VF0=VF
      ier=RELAX_OK
      didpT=.false.
      if (this%Tratmax.lt.1.0e30_WP.and.all(Q(1:4).gt.0.0_WP)) then
         TL=this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
         TG=this%gas%get_T_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])
         if (TL.gt.0.0_WP.and.TG.gt.0.0_WP) then
            if (max(TG/TL,TL/TG).gt.this%Tratmax) then
               ! Direction gate: capacity-weighted equilibrium estimate must not lift the
               ! liquid temperature by more than Tliftmax (quench semantics; see field doc)
               cvL=this%liq%get_cv_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
               cvG=this%gas%get_cv_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])
               Tstar=(Q(1)*cvL*TL+Q(2)*cvG*TG)/(Q(1)*cvL+Q(2)*cvG)
               if (Tstar.lt.this%Tliftmax*TL) then
                  call pTrelax_numerical(this%liq,this%gas,[1.0_WP],[1.0_WP],VF,Q,Pjump, &
                  &    itmax=this%itmax,rtol=this%rtol,RHOGmin=this%RHOGmin,ierr=ier)
                  didpT=(ier.eq.RELAX_OK)
               end if
            end if
         end if
      end if
      if (.not.didpT) call this%relax_numerical%apply(dt,VF,Q,Pjump,ier)
      VF=max(1.0_WP-this%VFratmax*(1.0_WP-VF0),VF0/this%VFratmax, &
      &      min(VF,1.0_WP-(1.0_WP-VF0)/this%VFratmax,this%VFratmax*VF0))
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
      ! Stage 2 -- set energies (unconditional fixed-VF split): impose the unique energy
      ! split with PL-PG=Pjump at frozen VF and masses. Newton in PL on
      ! E(PL) = VF*rhoe_L(PL) + (1-VF)*rhoe_G(PL-Pjump) = Etot, analytic slope
      ! Asum = dE/dPL = VF/GammaL + (1-VF)/GammaG; exact in one step for Gruneisen-form
      ! EOS (rhoe affine in p at fixed rho), short fixed-point otherwise.
      rhoL=Q(1)/VF; rhoG=Q(2)/(1.0_WP-VF); Etot=Q(3)+Q(4)
      PL=this%liq%get_p_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=[1.0_WP])
      Asum=0.0_WP
      do k=1,5
         QL=VF*this%liq%get_rhoe_from_p_rho(p=PL,rho=rhoL,y=[1.0_WP])
         QG=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=PL-Pjump,rho=rhoG,y=[1.0_WP])
         GL=this%liq%get_gruneisen_from_rho_e(rho=rhoL,e=QL/Q(1),y=[1.0_WP])
         GG=this%gas%get_gruneisen_from_rho_e(rho=rhoG,e=QG/Q(2),y=[1.0_WP])
         if (GL.le.0.0_WP.or.GG.le.0.0_WP) then   ! non-Gruneisen state: keep the stage-1 result
            if (present(ierr)) ierr=ier
            return
         end if
         Asum=VF/GL+(1.0_WP-VF)/GG
         dPL=(Etot-(QL+QG))/Asum
         PL=PL+dPL
         if (abs(dPL).le.this%rtol*max(abs(PL),1.0_WP)) exit
      end do
      ! Stage 3 -- floor: if the shared pressure sits below any user limit (PL,TL,PG,TG),
      ! raise it to the binding limit (all four rise monotonically with the shared
      ! pressure); the energy added is Asum*(Ptar-PL). Ledgered.
      Ptar=-1.0e30_WP
      if (this%Pmin_liq.gt.-1.0e29_WP) Ptar=max(Ptar,this%Pmin_liq)
      if (this%Pmin_gas.gt.-1.0e29_WP) Ptar=max(Ptar,this%Pmin_gas+Pjump)
      if (this%Tmin_liq.gt.0.0_WP) Ptar=max(Ptar,this%liq%get_p_from_rho_T(rho=rhoL,T=this%Tmin_liq,y=[1.0_WP]))
      if (this%Tmin_gas.gt.0.0_WP) Ptar=max(Ptar,this%gas%get_p_from_rho_T(rho=rhoG,T=this%Tmin_gas,y=[1.0_WP])+Pjump)
      if (PL.lt.Ptar) then
         this%acc(5)=this%acc(5)+1.0_WP
         this%acc(6)=this%acc(6)+Asum*(Ptar-PL)*this%vol
         PL=Ptar
      end if
      ! Write the split
      Q(3)=VF*this%liq%get_rhoe_from_p_rho(p=PL,rho=rhoL,y=[1.0_WP])
      Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=PL-Pjump,rho=rhoG,y=[1.0_WP])
      if (present(ierr)) ierr=RELAX_OK
   end subroutine apply

end module safe_relax_num_class
