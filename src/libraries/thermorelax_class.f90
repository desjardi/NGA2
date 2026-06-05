!> Abstract thermodynamic relaxation model
module thermorelax_class
   use precision, only: WP
   implicit none
   private
   public :: thermorelax
   type, abstract :: thermorelax
   contains
      procedure(relax_iface), deferred :: apply
   end type thermorelax
   abstract interface
      subroutine relax_iface(this,dt,VF,Q,Pjump)
         import :: thermorelax,WP
         class(thermorelax),     intent(inout) :: this
         real(WP),               intent(in)    :: dt
         real(WP),               intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
         real(WP),               intent(in)    :: Pjump
      end subroutine relax_iface
   end interface
end module thermorelax_class
