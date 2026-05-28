!> Abstract base class for AMR solvers
!> Provides common data (amr pointer, name) and enforces deferred interfaces for checkpoint I/O and solver info
module amrsolver_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrio_class,      only: amrio
   implicit none
   private

   public :: amrsolver

   !> Abstract solver base - extended by concrete solvers (amrscalar, etc.)
   type, abstract :: amrsolver
      ! Pointer to amrgrid
      class(amrgrid), pointer :: amr => null()
      ! Solver name for logging
      character(len=str_medium) :: name = 'UNNAMED_SOLVER'

   contains
      ! Deferred - concrete solvers must implement
      procedure(info_iface), deferred :: get_info
      procedure(chkpt_iface), deferred :: register_checkpoint
      procedure(restore_iface), deferred :: restore_checkpoint
   end type amrsolver

   !> Interface for get_info
   abstract interface
      subroutine info_iface(this)
         import :: amrsolver
         class(amrsolver), intent(inout) :: this
      end subroutine info_iface
   end interface

   !> Interface for register_checkpoint
   abstract interface
      subroutine chkpt_iface(this, io)
         import :: amrsolver, amrio
         class(amrsolver), intent(inout) :: this
         class(amrio), intent(inout) :: io
      end subroutine chkpt_iface
   end interface

   !> Interface for restore_checkpoint
   abstract interface
      subroutine restore_iface(this, io, dirname, time)
         import :: amrsolver, amrio, WP
         class(amrsolver), intent(inout) :: this
         class(amrio), intent(inout) :: io
         character(len=*), intent(in) :: dirname
         real(WP), intent(in) :: time
      end subroutine restore_iface
   end interface

end module amrsolver_class

