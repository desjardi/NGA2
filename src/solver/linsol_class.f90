!> Abstract linear solver concept is defined here:
!> given a config, it provides parallel solvers
module linsol_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: linsol
   
   
   !> linsol object definition
   type, abstract :: linsol
      
      ! A linear solver works for a config
      type(config), pointer :: cfg                                    !< Config for linsol
      
      ! An linsol has a name
      character(len=str_medium) :: name                               !< Name of solver
      
      ! Specific method to use
      integer  :: method                                              !< Solution method (we assume all solvers will have some sort of options)
      
      ! An linsol has a stencil size
      integer  :: nst                                                 !< Stencil size in 3D
      integer, dimension(:,:),   allocatable :: stc                   !< Stencil map in 3D (from stencil entry to i/j/k shift)
      integer, dimension(:,:,:), allocatable :: stmap                 !< Inverse stencil map in 3D (from i/j/k shift to stencil entry)
      
      ! An linsol stores the linear operator and rhs
      real(WP), dimension(:,:,:,:), allocatable :: opr                !< Linear operator
      real(WP), dimension(:,:,:),   allocatable :: rhs                !< RHS
      real(WP), dimension(:,:,:),   allocatable :: sol                !< Solution
      
      ! Is the solver setup?
      logical :: setup_done                                           !< Check whether the solver has been setup
      
      ! Convergence criteria
      integer  :: maxit                                               !< Maximum number of iterations allowed
      real(WP) :: rcvg                                                !< Desired relative convergence criterion
      real(WP) :: acvg                                                !< Desired absolue convergence criterion
      
      ! Current convergence info
      integer  :: it                                                  !< Current number of iterations
      real(WP) :: rerr                                                !< Current relative error
      real(WP) :: aerr                                                !< Current absolute error
      
   contains
      procedure(in_noarg_interface   ), deferred :: print_short             !< One-line printing of solver status
      procedure(in_noarg_interface   ), deferred :: print                   !< Long-form printing of solver status
      procedure(in_noarg_interface   ), deferred :: log                     !< Long-form logging of solver status
      procedure(inout_noarg_interface), deferred :: init                    !< Grid and stencil initialization - done once for the grid and stencil
      procedure(inout_noarg_interface), deferred :: setup                   !< Solver setup (every time the operator changes)
      procedure(inout_noarg_interface), deferred :: solve                   !< Execute solver (assumes new RHS and initial guess at every call)
      procedure(inout_noarg_interface), deferred :: destroy                 !< Solver destruction (every time the operator changes)
   end type linsol
   
   
   !> Interface
   abstract interface
      subroutine in_noarg_interface(this)
         import linsol
         class(linsol), intent(in) :: this
      end subroutine in_noarg_interface
      subroutine inout_noarg_interface(this)
         import linsol
         class(linsol), intent(inout) :: this
      end subroutine inout_noarg_interface
   end interface
   
   
end module linsol_class
