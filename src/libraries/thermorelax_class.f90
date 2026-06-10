!> Abstract thermodynamic relaxation model
module thermorelax_class
   use precision, only: WP
   implicit none
   private
   public :: thermorelax
   !> Relaxation status codes, returned via the optional ierr argument of apply
   !> (and the p_relax/pT_relax/pTg_relax workers). RELAX_OK=0 means relaxed cleanly.
   public :: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_VACUUM_GAS,RELAX_DEGENERATE,RELAX_SINGULAR
   integer, parameter :: RELAX_OK        =0   !< Relaxed successfully; both phases sound
   integer, parameter :: RELAX_FAILED    =1   !< No convergence / no real equilibrium found
   integer, parameter :: RELAX_BAD_LIQUID=2   !< Unphysical liquid (cL<=0: past co-volume or p<=-pinf)
   integer, parameter :: RELAX_BAD_GAS   =3   !< Unphysical gas (cG<=0)
   integer, parameter :: RELAX_VACUUM_GAS=4   !< Near-vacuum gas (rhoG<RHOGmin)
   integer, parameter :: RELAX_DEGENERATE=5   !< Degenerate cell (Q(1:4)<=0 or VF outside (0,1))
   integer, parameter :: RELAX_SINGULAR  =6   !< Singular Jacobian (numerical relax only)
   type, abstract :: thermorelax
   contains
      procedure(relax_iface), deferred :: apply
   end type thermorelax
   abstract interface
      subroutine relax_iface(this,dt,VF,Q,Pjump,ierr)
         import :: thermorelax,WP
         class(thermorelax),     intent(inout) :: this
         real(WP),               intent(in)    :: dt
         real(WP),               intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
         real(WP),               intent(in)    :: Pjump
         integer,  optional,     intent(out)   :: ierr
      end subroutine relax_iface
   end interface
end module thermorelax_class
