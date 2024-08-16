!> Incompressible scalar solver class built onto an amrcore
module amrscalar_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrconfig_class,  only: amrconfig
   use amrex_amr_module, only: amrex_multifab
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrscalar
   
   !> Constant density scalar solver object definition
   type :: amrscalar
      
      ! This is our amrconfig
      class(amrconfig), pointer :: amr                          !< This is the amrconfig the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_AMRSCALAR'     !< Solver name (default=UNNAMED_AMRSCALAR)
      
      ! Scalar and old scalar variables
      type(amrex_multifab), dimension(:), allocatable ::  SC    !< SC multifab array
      type(amrex_multifab), dimension(:), allocatable ::  SCold !< Old SC multifab array
      
      ! Number of components and ghost cells
      integer :: ncomp=1
      integer :: nover=0
      
      ! Monitoring quantities
      real(WP) :: SCmax,SCmin,SCint                             !< Maximum and minimum, integral scalar
      
   contains
      procedure :: initialize                                   !< Initialize scalar solver
      procedure :: finalize                                     !< Finalize scalar solver
      procedure :: clr_lvl
   end type amrscalar
   
   
contains
   
   
   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,name)
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrconfig), target, intent(in) :: amr
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to amrconfig object
      this%amr=>amr
      
      ! Allocate variables
      allocate(this%SC   (0:this%amr%nlvl-1))
      allocate(this%SCold(0:this%amr%nlvl-1))
      
   end subroutine initialize


   !> Clear solver data at level lev - should this be bind(c)?
   subroutine clr_lvl(this,lev)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in), value :: lev
      call amrex_multifab_destroy(this%SC   (lev))
      call amrex_multifab_destroy(this%SCold(lev))
   end subroutine clr_lvl
   
   
   !> Finalization for amrscalar solver
   subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: n
      do n=0,this%amr%nlvl-1
         call amrex_multifab_destroy(this%SC   (n))
         call amrex_multifab_destroy(this%SCold(n))
      end do
      deallocate(this%SC,this%SCold)
   end subroutine finalize
   
   
end module amrscalar_class
