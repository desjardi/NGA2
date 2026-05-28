!> Flux register class wrapping AMReX FluxRegister
!> Designed to be managed by amrgrid Registry
module amrflux_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_associated
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_fluxregister,amrex_multifab,&
   &                           amrex_boxarray,amrex_distromap,amrex_spacedim,&
   &                           amrex_real,amrex_geometry
   use amrex_interface,  only: amrfluxreg_build,amrfluxreg_destroy
   implicit none
   private

   public :: amrflux

   !> Flux register object wrapping array of FluxRegisters (one per fine level)
   type :: amrflux
      ! Underlying AMReX objects (array 1:maxlvl, no level 0)
      type(amrex_fluxregister), dimension(:), allocatable :: fr
      ! Metadata
      character(len=str_medium) :: name='UNNAMED_AMRFLUX'
      integer :: ncomp=1   !< Number of components
   contains
      procedure :: initialize    !< Initialize amrflux with parameters
      procedure :: finalize      !< Finalize amrflux
      procedure :: reset_level   !< Reset flux register at a fine level
      procedure :: clear_level   !< Clear flux register at a fine level
      ! Wrapper methods that forward to underlying fluxregister
      procedure :: crseinit
      procedure :: fineadd
      procedure :: reflux
      procedure :: setval
   end type amrflux

contains

   !> Initialize amrflux with parameters
   !> Allocates the fr array
   subroutine initialize(this, maxlvl, name, ncomp)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: maxlvl              !< Max level (from amr%maxlvl)
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      ! Set metadata
      this%name = name
      this%ncomp = ncomp
      ! Allocate fr array (flux registers only exist for levels >= 1)
      if (.not.allocated(this%fr)) allocate(this%fr(1:maxlvl))
   end subroutine initialize

   !> Finalize flux register object
   subroutine finalize(this)
      class(amrflux), intent(inout) :: this
      integer :: i
      if (allocated(this%fr)) then
         do i=1,ubound(this%fr,1)
            if (c_associated(this%fr(i)%p)) call amrfluxreg_destroy(this%fr(i)%p)
         end do
         deallocate(this%fr)
      end if
   end subroutine finalize

   !> Reset flux register at a fine level
   subroutine reset_level(this,lvl,ba,dm,ref_ratio)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      integer, intent(in) :: ref_ratio(3)
      ! Destroy if already built
      if (c_associated(this%fr(lvl)%p)) call amrfluxreg_destroy(this%fr(lvl)%p)
      ! Build new register with per-direction ref_ratio
      this%fr(lvl)%owner = .false.  ! We manage lifecycle via our C wrapper
      this%fr(lvl)%flev  = lvl
      call amrfluxreg_build(this%fr(lvl)%p, ba%p, dm%p, ref_ratio, lvl, this%ncomp)
   end subroutine reset_level

   !> Clear flux register at a specific level
   subroutine clear_level(this,lvl)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      if (c_associated(this%fr(lvl)%p)) call amrfluxreg_destroy(this%fr(lvl)%p)
      this%fr(lvl)%p = c_null_ptr
   end subroutine clear_level

   !> Initialize coarse-side fluxes (resets the register)
   subroutine crseinit(this,lvl,fluxes,scale)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
      real(amrex_real), intent(in) :: scale
      call this%fr(lvl)%crseinit(fluxes,scale)
   end subroutine crseinit

   !> Add fine-side fluxes
   subroutine fineadd(this,lvl,fluxes,scale)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
      real(amrex_real), intent(in) :: scale
      call this%fr(lvl)%fineadd(fluxes,scale)
   end subroutine fineadd

   !> Apply correction to coarse-level data
   subroutine reflux(this,lvl,mf,scale,geom)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_multifab), intent(inout) :: mf
      real(amrex_real), intent(in) :: scale
      type(amrex_geometry), intent(in) :: geom
      call this%fr(lvl)%reflux(mf,scale)
   end subroutine reflux

   !> Set all values in the register
   subroutine setval(this,lvl,val)
      class(amrflux), intent(inout) :: this
      integer, intent(in) :: lvl
      real(amrex_real), intent(in) :: val
      call this%fr(lvl)%setval(val)
   end subroutine setval

end module amrflux_class
