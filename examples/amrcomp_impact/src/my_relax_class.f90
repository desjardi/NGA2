!> Impact-specific extension of relax_ig_sg
!> Adds:
!>   - Cavitation check for pure-liquid cells (energy injection when PL < -0.9*pinf)
!>   - Naive air-dissolution clip after mechanical relax (when Peq > Peq_diss)
module my_relax_class
   use precision,         only: WP
   use relax_ig_sg_class, only: relax_ig_sg,Prelax,PTrelax
   implicit none
   private

   public :: my_relax

   type, extends(relax_ig_sg) :: my_relax
      real(WP) :: PL_cav  =-0.9_WP             !< Cavitation target pressure as a factor of pinf
      real(WP) :: Peq_diss=200.0_WP            !< Dissolution clip threshold on the equilibrium pressure
   contains
      procedure :: apply
      procedure :: p_relax
   end type my_relax

contains

   !> Apply: cavitation branch for pure-liquid cells, then dispatch to overridden p_relax/pT_relax
   subroutine apply(this,dt,VF,Q,Pjump)
      use messager, only: die
      implicit none
      class(my_relax),        intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP) :: PL
      ! Cavitation territory: pure-liquid cells
      if (VF.ge.1.0_WP) then
         if (Q(1).le.0.0_WP.or.Q(3).le.0.0_WP) return
         PL=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
         if (PL.lt.this%PL_cav*this%liq%pinf) then
            print*,'Warning: cavitation detected, injecting energy'
            Q(3)=VF*(this%PL_cav*this%liq%pinf+this%liq%gamma*this%liq%pinf)/(this%liq%gamma-1.0_WP)
         end if
         return
      end if
      ! Mixture cells only
      if (VF.le.0.0_WP) return
      ! Dispatch (this%p_relax hits the override polymorphically since this is class(my_relax) here)
      select case (this%model)
      case (Prelax);  call this%p_relax (dt,VF,Q,Pjump)
      case (PTrelax); call this%pT_relax(dt,VF,Q,Pjump)
      case default; call die('[my_relax apply] unknown model')
      end select
   end subroutine apply

   !> p_relax: run parent's mechanical relax, then apply dissolution clip if Peq exceeds threshold
   subroutine p_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(my_relax),        intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      real(WP) :: Peq
      ! Run parent's mechanical relax
      call this%relax_ig_sg%p_relax(dt,VF,Q,Pjump)
      ! Post-relax dissolution check: compute equilibrium pressure from updated state
      Peq=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
      if (Peq.gt.this%Peq_diss) then
         print*,'Warning: Pressure relaxation hit upper limit, applying naive model for air dissolution into water'
         Q(1)=Q(1)/VF
         Q(2)=0.0_WP
         Q(3)=Q(3)/VF
         Q(4)=0.0_WP
         VF  =1.0_WP
      end if
   end subroutine p_relax

end module my_relax_class
