!> Data class wrapping AMReX MultiFab
!> Designed to be managed by amrgrid Registry
module amrdata_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrex_amr_module, only: amrex_multifab,amrex_boxarray,amrex_distromap,&
   &                           amrex_geometry,amrex_interp_pc,amrex_bc_foextrap,&
   &                           amrex_interp_cell_cons,amrex_interp_face_linear,amrex_interp_node_bilinear
   implicit none
   private

   public :: amrdata
   public :: amrdata_on_init,amrdata_on_coarse,amrdata_on_remake,amrdata_on_clear,amrdata_fillbc
   public :: amrex_interp_pc,amrex_bc_foextrap,amrex_interp_cell_cons,amrex_interp_face_linear,amrex_interp_node_bilinear
   public :: default_fillbc

   ! Special interpolation modes for amrdata
   integer, parameter, public :: amrex_interp_none   = -1  !< Workspace: allocate but don't fill
   integer, parameter, public :: amrex_interp_reinit = -2  !< Reinit: always call user_init on regrid

   ! Forward declaration for interfaces
   type :: amrdata
      ! The underlying AMReX objects (array of MultiFabs, one per level)
      type(amrex_multifab), dimension(:), allocatable :: mf
      ! Pointer to parent amrgrid (set by initialize)
      class(amrgrid), pointer :: amr => null()
      ! Pointer to parent solver (set by solver, for BC context access)
      class(*), pointer :: parent => null()
      ! Metadata
      character(len=str_medium) :: name='UNNAMED_AMRDATA'
      integer :: ncomp=1                                   !< Number of components (default=1)
      integer :: ng=0                                      !< Number of ghost cells (default=0)
      logical :: nodal(3) = [.false., .false., .false.]    !< false=cell, true=vertex in that direction
      integer :: interp=amrex_interp_cell_cons             !< Interpolation method
      integer, dimension(:,:), allocatable :: lo_bc,hi_bc  !< Boundary conditions: lo_bc(3,ncomp), hi_bc(3,ncomp)
      ! Cache level index for fillbc callback
      integer :: fill_lvl_cache=-1
      ! Callback pointers (set to defaults in initialize)
      procedure(on_init_iface),   pointer, pass :: on_init   => null()
      procedure(on_coarse_iface), pointer, pass :: on_coarse => null()
      procedure(on_remake_iface), pointer, pass :: on_remake => null()
      procedure(on_clear_iface),  pointer, pass :: on_clear  => null()
      procedure(fillbc_iface),    pointer, pass :: fillbc    => null()
      ! User-provided initialization callback
      procedure(on_init_iface),   pointer, pass :: user_init => null()
   contains
      ! Lifecycle methods
      procedure :: initialize       !< Initialize amrdata with amrgrid and parameters
      procedure :: finalize         !< Finalize the amrdata object
      procedure :: register         !< Register regrid callbacks with amrgrid
      ! Level management
      procedure :: reset_level      !< Regenerates mfab at level
      procedure :: clear_level      !< Clears mfab at level
      ! Fill operations
      procedure :: fill_from_coarse !< Interpolate from coarse level only (single level)
      procedure :: fill_lvl         !< Fill ghost cells at single level
      procedure :: fill             !< Fill ghost cells on all levels
      procedure :: fill_mfab        !< Fill into a target MultiFab (single level)
      procedure :: sync_lvl         !< Lightweight same-level ghost sync (single level)
      procedure :: sync             !< Lightweight ghost sync on all levels
      procedure :: syncsum_lvl      !< Same-level ghost-to-valid accumulation (single level)
      procedure :: syncsum          !< Ghost-to-valid accumulation on all levels
      procedure :: average_down     !< Average from finest to coarsest level
      procedure :: average_downto   !< Average level lvl+1 down to level lvl
      procedure :: sum_down         !< Restrict-SUM from finest to coarsest level
      procedure :: sum_downto       !< Restrict-SUM level lvl+1 into level lvl
      ! Scalar operations (Y = op(Y, scalar))
      procedure :: setval           !< Y = val
      procedure :: plus             !< Y = Y + val
      procedure :: mult             !< Y = Y * val
      procedure :: clip             !< Y = clip(Y, minval, maxval)
      ! Binary operations (Y = op(Y, X))
      procedure :: add              !< Y = Y + X
      procedure :: subtract         !< Y = Y - X
      procedure :: multiply         !< Y = Y * X (element-wise)
      procedure :: divide           !< Y = Y / X (element-wise)
      procedure :: copy             !< Y = X
      ! BLAS-like operations
      procedure :: saxpy            !< Y = aX + Y
      procedure :: lincomb          !< Y = aX1 + bX2
      ! Reduction operations (require explicit lvl)
      procedure :: get_min          !< min value at level
      procedure :: get_max          !< max value at level
      procedure :: get_sum          !< sum at level
      procedure :: norm0            !< L-infinity norm at level
      procedure :: norm1            !< L1 norm at level
      procedure :: norm2            !< L2 norm at level
      ! Vector operations
      procedure :: get_magnitude    !< M = sqrt(srcX(compX)²+srcY(compY)²+srcZ(compZ)²), comps default to 1
      ! Iteration helper
      procedure :: mfiter_build     !< Build MFIter from this data's MultiFab
      ! Cloning
      procedure :: clone            !< Clone amrdata with optional alternate DM per level
   end type amrdata

   !> Abstract interface for on_init callback
   abstract interface
      subroutine on_init_iface(this, lvl, time, ba, dm)
         import :: amrdata, amrex_boxarray, amrex_distromap, WP
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine on_init_iface
   end interface

   !> Abstract interface for on_coarse callback
   abstract interface
      subroutine on_coarse_iface(this, lvl, time, ba, dm)
         import :: amrdata, amrex_boxarray, amrex_distromap, WP
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine on_coarse_iface
   end interface

   !> Abstract interface for on_remake callback
   abstract interface
      subroutine on_remake_iface(this, lvl, time, ba, dm)
         import :: amrdata, amrex_boxarray, amrex_distromap, WP
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine on_remake_iface
   end interface

   !> Abstract interface for on_clear callback
   abstract interface
      subroutine on_clear_iface(this, lvl)
         import :: amrdata
         class(amrdata), intent(inout) :: this
         integer, intent(in) :: lvl
      end subroutine on_clear_iface
   end interface

   !> Abstract interface for fillbc callback (matches AMReX order)
   abstract interface
      subroutine fillbc_iface(this, mf, scomp, ncomp, time, geom)
         import :: amrdata, amrex_multifab, amrex_geometry, WP
         class(amrdata), intent(inout) :: this
         type(amrex_multifab), intent(inout) :: mf
         integer, intent(in) :: scomp, ncomp
         real(WP), intent(in) :: time
         type(amrex_geometry), intent(in) :: geom
      end subroutine fillbc_iface
   end interface

contains

   !> Initialize amrdata with amrgrid reference and parameters
   !> This sets the amr pointer, allocates mf array, and configures BCs
   subroutine initialize(this, amr, name, ncomp, ng, nodal, interp)
      use amrex_amr_module, only: amrex_bc_int_dir
      use amrex_amr_module, only: amrex_interp_face_linear,amrex_interp_node_bilinear
      implicit none
      class(amrdata), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      integer, intent(in) :: ng
      logical, intent(in), optional :: nodal(3)  !< Data location: false=cell, true=vertex
      integer, intent(in), optional :: interp    !< Interpolation method override
      ! Store pointer to amrgrid
      this%amr => amr
      ! Set metadata
      this%name = name
      this%ncomp = ncomp
      this%ng = ng
      ! Set nodal if provided
      if (present(nodal)) this%nodal = nodal
      ! Set interpolation: user override or auto-select based on nodal
      if (present(interp)) then
         this%interp = interp
      else if (any(this%nodal)) then
         if (all(this%nodal)) then
            this%interp = amrex_interp_node_bilinear  ! Node-centered
         else
            this%interp = amrex_interp_face_linear    ! Face-centered
         end if
      end if
      ! Allocate mf array
      if (.not.allocated(this%mf)) allocate(this%mf(0:amr%maxlvl))
      ! Allocate and set default BCs to handle periodicity only
      if (.not.allocated(this%lo_bc)) allocate(this%lo_bc(3,ncomp))
      if (.not.allocated(this%hi_bc)) allocate(this%hi_bc(3,ncomp))
      this%lo_bc = amrex_bc_int_dir
      this%hi_bc = amrex_bc_int_dir
      ! Set default callbacks
      this%on_init   => default_on_init
      this%on_coarse => default_on_coarse
      this%on_remake => default_on_remake
      this%on_clear  => default_on_clear
      this%fillbc    => default_fillbc
   end subroutine initialize

   !> Finalize the amrdata object
   subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy
      class(amrdata), intent(inout) :: this
      integer :: i
      ! Destroy all MultiFabs
      if (allocated(this%mf)) then
         do i=0,ubound(this%mf,1)
            call amrex_multifab_destroy(this%mf(i))
         end do
         deallocate(this%mf)
      end if
      ! Deallocate BC arrays
      if (allocated(this%lo_bc)) deallocate(this%lo_bc)
      if (allocated(this%hi_bc)) deallocate(this%hi_bc)
      ! Reset pointers
      nullify(this%amr)
      nullify(this%parent)
      nullify(this%on_init)
      nullify(this%on_coarse)
      nullify(this%on_remake)
      nullify(this%on_clear)
      nullify(this%fillbc)
      nullify(this%user_init)
      ! Reset scalars to defaults
      this%name = 'UNNAMED_AMRDATA'
      this%ncomp = 1
      this%ng = 0
      this%nodal = [.false., .false., .false.]
      this%interp = amrex_interp_cell_cons
   end subroutine finalize

   !> Register regrid callbacks with amrgrid
   subroutine register(this)
      use iso_c_binding, only: c_loc
      use messager, only: die
      class(amrdata), target, intent(inout) :: this
      if (.not.associated(this%amr)) call die('[amrdata%register] armdata%initialize must be called first')
      select type (this)
       type is (amrdata)
         call this%amr%add_on_init  (amrdata_on_init,   c_loc(this))
         call this%amr%add_on_coarse(amrdata_on_coarse, c_loc(this))
         call this%amr%add_on_remake(amrdata_on_remake, c_loc(this))
         call this%amr%add_on_clear (amrdata_on_clear,  c_loc(this))
      end select
   end subroutine register

   !> Reset mfab on a level given new BoxArray and DistroMap
   subroutine reset_level(this,lvl,ba,dm)
      use amrex_amr_module, only: amrex_multifab_build,amrex_multifab_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      call amrex_multifab_destroy(this%mf(lvl))
      call amrex_multifab_build(this%mf(lvl),ba,dm,this%ncomp,this%ng,this%nodal)
   end subroutine reset_level

   !> Destroy mfab on a level
   subroutine clear_level(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      call amrex_multifab_destroy(this%mf(lvl))
   end subroutine clear_level

   !> Build MFIter from this data's MultiFab (ensures correct layout for staggered data)
   subroutine mfiter_build(this, lvl, mfi, tiling)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(inout) :: mfi
      logical, intent(in), optional :: tiling
      logical :: use_tiling
      use_tiling = this%amr%default_tiling
      if (present(tiling)) use_tiling = tiling
      call amrex_mfiter_build(mfi, this%mf(lvl), tiling=use_tiling)
   end subroutine mfiter_build

   !> Default fillbc: applies amrex_filcc using lo_bc/hi_bc
   !> NOTE: amrex_filcc is designed for cell-centered data only.
   subroutine default_fillbc(this, mf, scomp, ncomp, time, geom)
      use amrex_amr_module, only: amrex_filcc, amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer, dimension(4) :: plo, phi
      ! Skip if all directions periodic
      if (this%amr%xper .and. this%amr%yper .and. this%amr%zper) return
      ! Loop over boxes and apply filcc
      call amrex_mfiter_build(mfi, mf, tiling=.false.)
      do while (mfi%next())
         p => mf%dataptr(mfi)
         if (.not.geom%domain%contains(p)) then
            plo = lbound(p); phi = ubound(p)
            call amrex_filcc(p, plo, phi, geom%domain%lo, geom%domain%hi, geom%dx, &
            &   geom%get_physical_location(plo), this%lo_bc, this%hi_bc)
         end if
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine default_fillbc

   !> Default on_init: reset level and initialize based on interp mode
   subroutine default_on_init(this, lvl, time, ba, dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level
      call this%reset_level(lvl, ba, dm)
      ! Handle different modes
      select case (this%interp)
       case (amrex_interp_none)
         ! Workspace mode: just allocate, don't fill
       case default
         ! Standard or reinit mode: zero out, then call user_init
         call this%mf(lvl)%setval(0.0_WP)
      end select
      ! User-provided initialization (called for all modes except none without user_init)
      if (associated(this%user_init)) call this%user_init(lvl, time, ba, dm)
   end subroutine default_on_init

   !> Default on_coarse: reset level and fill based on interp mode
   subroutine default_on_coarse(this, lvl, time, ba, dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level
      call this%reset_level(lvl, ba, dm)
      ! Handle different modes
      select case (this%interp)
       case (amrex_interp_none)
         ! Workspace mode: just allocate, don't fill
       case (amrex_interp_reinit)
         ! Reinit mode: call user_init instead of interpolating
         if (associated(this%user_init)) call this%user_init(lvl, time, ba, dm)
       case default
         ! Standard interpolation: fill from coarse
         call this%fill_from_coarse(lvl, time)
      end select
   end subroutine default_on_coarse

   !> Default on_remake: handle level relayout based on interp mode
   subroutine default_on_remake(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_multifab_build, amrex_multifab_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab) :: mf_tmp
      ! Handle different modes
      select case (this%interp)
       case (amrex_interp_none)
         ! Workspace mode: just reallocate
         call this%reset_level(lvl, ba, dm)
       case (amrex_interp_reinit)
         ! Reinit mode: reallocate and call user_init
         call this%reset_level(lvl, ba, dm)
         if (associated(this%user_init)) call this%user_init(lvl, time, ba, dm)
       case default
         ! Standard interpolation: FillPatch old data into new layout
         ! Build temp MultiFab with new layout (0 ghost cells for FillPatch)
         call amrex_multifab_build(mf_tmp, ba, dm, this%ncomp, 0, this%nodal)
         ! Fill temp from old data via FillPatch
         call this%fill_mfab(mf_tmp, lvl, time)
         ! Reset level and copy from temp
         call this%reset_level(lvl, ba, dm)
         call this%mf(lvl)%copy(mf_tmp, 1, 1, this%ncomp, 0)
         ! Destroy temp
         call amrex_multifab_destroy(mf_tmp)
      end select
   end subroutine default_on_remake

   !> Default on_clear: just clear the level
   subroutine default_on_clear(this, lvl)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%clear_level(lvl)
   end subroutine default_on_clear

   !> Dispatch callback for fillbc
   subroutine amrdata_fillbc(ctx, mf_ptr, scomp, ncomp, time, geom_ptr) bind(c)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_double, c_int
      use amrex_amr_module, only: amrex_geometry, amrex_multifab
      type(c_ptr), value, intent(in) :: ctx
      type(c_ptr), value, intent(in) :: mf_ptr
      integer(c_int), value, intent(in) :: scomp
      integer(c_int), value, intent(in) :: ncomp
      real(c_double), value, intent(in) :: time
      type(c_ptr), value, intent(in) :: geom_ptr
      type(amrdata), pointer :: this
      type(amrex_multifab) :: mf
      type(amrex_geometry) :: geom
      call c_f_pointer(ctx, this)
      mf = mf_ptr
      geom = geom_ptr
      ! Call the fillbc callback
      call this%fillbc(mf, int(scomp), int(ncomp), real(time, WP), geom)
   end subroutine amrdata_fillbc

   !> Dispatch callback for on_init
   subroutine amrdata_on_init(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
   end subroutine amrdata_on_init

   !> Dispatch callback for on_coarse
   subroutine amrdata_on_coarse(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrdata_on_coarse

   !> Dispatch callback for on_remake
   subroutine amrdata_on_remake(ctx, lvl, time, ba, dm)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrdata_on_remake

   !> Dispatch callback for on_clear
   subroutine amrdata_on_clear(ctx, lvl)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrdata), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrdata_on_clear

   !> Fill fine level from coarse only (for creating new fine levels)
   subroutine fill_from_coarse(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillcoarsepatch
      implicit none
      class(amrdata), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      ! Only valid for fine levels
      if (lvl .lt. 1) return
      ! Get context pointer
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(amrdata_fillbc)
      ! Call C++ wrapper
      this%fill_lvl_cache=lvl ! Cache current level
      call amrmfab_fillcoarsepatch(this%mf(lvl), time, this%mf(lvl-1), &
      &   this%amr%geom(lvl-1), this%amr%geom(lvl), data_ctx, bc_dispatch_ptr, &
      &   1, 1, this%ncomp, [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)], this%interp, this%lo_bc, this%hi_bc, this%ncomp)
   end subroutine fill_from_coarse

   !> Fill ghost cells and coarse-fine boundary data at a single level
   !> For lvl=0: fills ghost cells using boundary conditions
   !> For lvl>0: interpolates from coarser level and fills ghosts
   subroutine fill_lvl(this, lvl, time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrdata), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get context pointer (this amrdata)
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(amrdata_fillbc)
      ! Call appropriate FillPatch (scomp/dcomp use 1-indexed Fortran convention)
      this%fill_lvl_cache=lvl ! Cache current level
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(this%mf(0), t_old, this%mf(0), &
         &   t_new, this%mf(0), this%amr%geom(0), data_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, this%ncomp)
      else
         call amrmfab_fillpatch_two(this%mf(lvl), t_old, this%mf(lvl-1), &
         &   t_new, this%mf(lvl-1), this%amr%geom(lvl-1), &
         &   t_old, this%mf(lvl), t_new, this%mf(lvl), this%amr%geom(lvl), &
         &   data_ctx, bc_dispatch_ptr, time, 1, 1, this%ncomp, &
         &   [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)], this%interp, this%lo_bc, this%hi_bc, this%ncomp)
      end if
      ! For nodal/face data: reconcile shared valid faces
      if (any(this%nodal)) call this%mf(lvl)%override_sync(this%amr%geom(lvl))
   end subroutine fill_lvl

   !> Fill ghost cells and coarse-fine boundary data on all levels
   subroutine fill(this,time,lbase)
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in) :: time
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=lb,this%amr%clvl()
         call this%fill_lvl(lvl,time)
      end do
   end subroutine fill

   !> Fill into a target MultiFab from this amrdata (for regrid callbacks)
   subroutine fill_mfab(this,dest,lvl,time)
      use iso_c_binding, only: c_loc, c_funloc, c_funptr, c_ptr
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface, only: amrmfab_fillpatch_single, amrmfab_fillpatch_two
      implicit none
      class(amrdata), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: dest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP) :: t_old, t_new
      type(c_ptr) :: data_ctx
      type(c_funptr) :: bc_dispatch_ptr
      t_old = time - 1.0e200_WP
      t_new = time
      ! Get context pointer
      select type (this)
       type is (amrdata)
         data_ctx = c_loc(this)
      end select
      bc_dispatch_ptr = c_funloc(amrdata_fillbc)
      ! Call appropriate FillPatch
      this%fill_lvl_cache=lvl ! Cache current level
      if (lvl .eq. 0) then
         call amrmfab_fillpatch_single(dest, t_old, this%mf(0), &
         &   t_new, this%mf(0), this%amr%geom(0), data_ctx, bc_dispatch_ptr, &
         &   time, 1, 1, this%ncomp)
      else
         call amrmfab_fillpatch_two(dest, t_old, this%mf(lvl-1), &
         &   t_new, this%mf(lvl-1), this%amr%geom(lvl-1), &
         &   t_old, this%mf(lvl), t_new, this%mf(lvl), this%amr%geom(lvl), &
         &   data_ctx, bc_dispatch_ptr, time, 1, 1, this%ncomp, &
         &   [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)], this%interp, this%lo_bc, this%hi_bc, this%ncomp)
      end if
      ! For nodal/face data: reconcile shared valid faces
      if (any(this%nodal)) call dest%override_sync(this%amr%geom(lvl))
   end subroutine fill_mfab

   !> Lightweight same-level ghost sync at a single level (no C/F, no BCs except periodic)
   !> For nodal/face data, also reconciles shared valid faces
   subroutine sync_lvl(this, lvl)
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      ! For nodal/face data: reconcile shared valid faces first
      if (any(this%nodal)) call this%mf(lvl)%override_sync(this%amr%geom(lvl))
      ! Then fill ghosts from valid cells (same level only)
      call this%mf(lvl)%fill_boundary(this%amr%geom(lvl))
   end subroutine sync_lvl

   !> Lightweight same-level ghost sync on all levels
   subroutine sync(this,lbase)
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=lb,this%amr%clvl()
         call this%sync_lvl(lvl)
      end do
   end subroutine sync

   !> Same-level ghost-to-valid accumulation at a single level (SumBoundary, no C/F)
   subroutine syncsum_lvl(this,lvl)
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%mf(lvl)%sum_boundary(this%amr%geom(lvl))
   end subroutine syncsum_lvl

   !> Same-level ghost-to-valid accumulation on all levels
   subroutine syncsum(this,lbase)
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=lb,this%amr%clvl()
         call this%syncsum_lvl(lvl)
      end do
   end subroutine syncsum

   !> Average down from finest level to lbase (ensures level consistency)
   !> Simply calls average_downto in a loop from finest to coarsest
   subroutine average_down(this,lbase)
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      ! Loop from finest to coarsest
      do lvl=this%amr%clvl()-1,lb,-1
         call this%average_downto(lvl)
      end do
   end subroutine average_down

   !> Average level lvl+1 down to level lvl
   !> - nodal_count=0 (cell): amrmfab_average_down_cell
   !> - nodal_count=1 (face): amrmfab_average_down_face
   !> - nodal_count=2 (edge): amrmfab_average_down_edge
   !> - nodal_count=3 (node): amrmfab_average_down_node
   subroutine average_downto(this,lvl)
      use messager, only: die
      use amrex_interface, only: amrmfab_average_down_cell,amrmfab_average_down_face, &
      &                          amrmfab_average_down_edge,amrmfab_average_down_node
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Check that level is valid
      if (lvl.lt.0.or.lvl.ge.this%amr%clvl()) call die('[amrdata average_downto] invalid level provided')
      ! Pass geometry for periodic fix-up
      select case (count(this%nodal))
       case (0) ! Cell-centered
         call amrmfab_average_down_cell(fmf=this%mf(lvl+1),cmf=this%mf(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
       case (1) ! Face-centered
         call amrmfab_average_down_face(fmf=this%mf(lvl+1),cmf=this%mf(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
       case (2) ! Edge-centered
         call amrmfab_average_down_edge(fmf=this%mf(lvl+1),cmf=this%mf(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
       case (3) ! Node-centered
         call amrmfab_average_down_node(fmf=this%mf(lvl+1),cmf=this%mf(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
      end select
   end subroutine average_downto

   !> Restrict-SUM all fine deposits (valid+ghost) into the coarse level.
   !> Mirrors AMReX's sumFineToCrseNodal. lvl is the coarse destination; lvl+1 is the source.
   subroutine sum_downto(this,lvl)
      use amrex_interface, only: amrmfab_sum_downto
      use messager, only: die
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      if (lvl.lt.0.or.lvl.ge.this%amr%clvl()) call die('[amrdata sum_downto] invalid level')
      call amrmfab_sum_downto(this%mf(lvl+1),this%mf(lvl),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl),fgeom=this%amr%geom(lvl+1))
   end subroutine sum_downto

   !> Restrict-SUM from finest to coarsest, looping lvl=clvl()-1 down to 0
   subroutine sum_down(this,lbase)
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=this%amr%clvl()-1,lb,-1
         call this%sum_downto(lvl)
      end do
   end subroutine sum_down

   ! ============================================================================
   ! SCALAR OPERATIONS
   ! ============================================================================

   !> Set all values to val
   !> @param val Value to set
   !> @param lvl Optional: single level to operate on
   !> @param lbase Optional: lowest level (operates on lbase to clvl)
   subroutine setval(this, val, lvl, lbase)
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in) :: val
      integer, intent(in), optional :: lvl, lbase
      integer :: l, l0, l1
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      do l = l0, l1
         call this%mf(l)%setval(val)
      end do
   end subroutine setval

   !> Add scalar: Y = Y + val
   subroutine plus(this, val, lvl, lbase, comp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in) :: val
      integer, intent(in), optional :: lvl, lbase, comp, ncomp, nghost
      integer :: l, l0, l1, ic, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      ic = 1; if (present(comp)) ic = comp
      nc = this%ncomp; if (present(ncomp)) nc = ncomp
      ng = this%ng; if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%plus(val, ic, nc, ng)
      end do
   end subroutine plus

   !> Multiply by scalar: Y = Y * val
   subroutine mult(this, val, lvl, lbase, comp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in) :: val
      integer, intent(in), optional :: lvl, lbase, comp, ncomp, nghost
      integer :: l, l0, l1, ic, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      ic = 1; if (present(comp)) ic = comp
      nc = this%ncomp; if (present(ncomp)) nc = ncomp
      ng = this%ng; if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%mult(val, ic, nc, ng)
      end do
   end subroutine mult

   ! ============================================================================
   ! BINARY OPERATIONS
   ! ============================================================================

   !> Add: Y = Y + X
   subroutine add(this, src, lvl, lbase, srccomp, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrdata), intent(in) :: src
      integer, intent(in), optional :: lvl, lbase, srccomp, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc = 1; if (present(srccomp)) sc = srccomp
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, src%ncomp); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, src%ng); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%add(src%mf(l), sc, dc, nc, ng)
      end do
   end subroutine add

   !> Subtract: Y = Y - X
   subroutine subtract(this, src, lvl, lbase, srccomp, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrdata), intent(in) :: src
      integer, intent(in), optional :: lvl, lbase, srccomp, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc = 1; if (present(srccomp)) sc = srccomp
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, src%ncomp); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, src%ng); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%subtract(src%mf(l), sc, dc, nc, ng)
      end do
   end subroutine subtract

   !> Multiply element-wise: Y = Y * X
   subroutine multiply(this, src, lvl, lbase, srccomp, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrdata), intent(in) :: src
      integer, intent(in), optional :: lvl, lbase, srccomp, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc = 1; if (present(srccomp)) sc = srccomp
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, src%ncomp); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, src%ng); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%multiply(src%mf(l), sc, dc, nc, ng)
      end do
   end subroutine multiply

   !> Divide element-wise: Y = Y / X
   subroutine divide(this, src, lvl, lbase, srccomp, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrdata), intent(in) :: src
      integer, intent(in), optional :: lvl, lbase, srccomp, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc = 1; if (present(srccomp)) sc = srccomp
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, src%ncomp); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, src%ng); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%divide(src%mf(l), sc, dc, nc, ng)
      end do
   end subroutine divide

   !> Copy: Y = X
   subroutine copy(this, src, lvl, lbase, srccomp, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrdata), intent(in) :: src
      integer, intent(in), optional :: lvl, lbase, srccomp, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc = 1; if (present(srccomp)) sc = srccomp
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, src%ncomp); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, src%ng); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%copy(src%mf(l), sc, dc, nc, ng)
      end do
   end subroutine copy

   !> Clip values: Y = max(cliplo, min(cliphi, Y))
   subroutine clip(this, cliplo, cliphi, lvl, lbase, comp, ncomp, nghost)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in), optional :: cliplo, cliphi
      integer, intent(in), optional :: lvl, lbase, comp, ncomp, nghost
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: lo, hi
      integer :: i, j, k, n, l, l0, l1, ic, nc, ng
      if (.not.associated(this%amr)) return
      lo = -huge(1.0_WP); if (present(cliplo)) lo = cliplo
      hi = +huge(1.0_WP); if (present(cliphi)) hi = cliphi
      call get_level_range(this, lvl, lbase, l0, l1)
      ic = 1; if (present(comp)) ic = comp
      nc = this%ncomp; if (present(ncomp)) nc = ncomp
      ng = this%ng; if (present(nghost)) ng = nghost
      do l = l0, l1
         call amrex_mfiter_build(mfi, this%mf(l), tiling=.true.)
         do while (mfi%next())
            p => this%mf(l)%dataptr(mfi)
            bx = mfi%growntilebox(ng)
            do n = ic, ic+nc-1
               do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  p(i,j,k,n) = max(lo, min(hi, p(i,j,k,n)))
               end do; end do; end do
            end do
         end do
         call amrex_mfiter_destroy(mfi)
      end do
   end subroutine clip

   ! ============================================================================
   ! BLAS-LIKE OPERATIONS
   ! ============================================================================

   !> SAXPY: Y = a*X + Y
   subroutine saxpy(this, a, src, lvl, lbase, srccomp, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in) :: a
      class(amrdata), intent(in) :: src
      integer, intent(in), optional :: lvl, lbase, srccomp, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc = 1; if (present(srccomp)) sc = srccomp
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, src%ncomp); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, src%ng); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%saxpy(a, src%mf(l), sc, dc, nc, ng)
      end do
   end subroutine saxpy

   !> Linear combination: Y = a*X1 + b*X2
   subroutine lincomb(this, a, src1, b, src2, lvl, lbase, srccomp1, srccomp2, dstcomp, ncomp, nghost)
      implicit none
      class(amrdata), intent(inout) :: this
      real(WP), intent(in) :: a, b
      class(amrdata), intent(in) :: src1, src2
      integer, intent(in), optional :: lvl, lbase, srccomp1, srccomp2, dstcomp, ncomp, nghost
      integer :: l, l0, l1, sc1, sc2, dc, nc, ng
      if (.not.associated(this%amr)) return
      call get_level_range(this, lvl, lbase, l0, l1)
      sc1 = 1; if (present(srccomp1)) sc1 = srccomp1
      sc2 = 1; if (present(srccomp2)) sc2 = srccomp2
      dc = 1; if (present(dstcomp)) dc = dstcomp
      nc = min(this%ncomp, min(src1%ncomp, src2%ncomp)); if (present(ncomp)) nc = ncomp
      ng = min(this%ng, min(src1%ng, src2%ng)); if (present(nghost)) ng = nghost
      do l = l0, l1
         call this%mf(l)%lincomb(a, src1%mf(l), sc1, b, src2%mf(l), sc2, dc, nc, ng)
      end do
   end subroutine lincomb

   ! ============================================================================
   ! REDUCTION OPERATIONS (require explicit lvl)
   ! ============================================================================

   !> Get minimum value at level
   function get_min(this, lvl, comp) result(val)
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in) :: lvl
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      val = this%mf(lvl)%min(ic)
   end function get_min

   !> Get maximum value at level
   function get_max(this, lvl, comp) result(val)
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in) :: lvl
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      val = this%mf(lvl)%max(ic)
   end function get_max

   !> Get sum at level (uses sum_unique for nodal data to avoid double-counting)
   function get_sum(this, lvl, comp) result(val)
      use amrex_interface, only: amrmfab_sum_unique
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in) :: lvl
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      ! Use sum_unique for nodal data to avoid double-counting at periodic boundaries
      if (any(this%nodal)) then
         val = amrmfab_sum_unique(this%mf(lvl), this%amr%geom(lvl), ic)
      else
         val = this%mf(lvl)%sum(ic)
      end if
   end function get_sum

   !> Get L-infinity norm (max abs) at level
   function norm0(this, lvl, comp) result(val)
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in) :: lvl
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      val = this%mf(lvl)%norm0(ic)
   end function norm0

   !> Get L1 norm at level
   function norm1(this, lvl, comp) result(val)
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in) :: lvl
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      val = this%mf(lvl)%norm1(ic)
   end function norm1

   !> Get L2 norm at level
   function norm2(this, lvl, comp) result(val)
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in) :: lvl
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      val = this%mf(lvl)%norm2(ic)
   end function norm2

   !> Compute magnitude: this = sqrt(srcX² + srcY² + srcZ²)
   subroutine get_magnitude(this,srcX,srcY,srcZ,compX,compY,compZ,nghost)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
      class(amrdata), intent(inout) :: this
      class(amrdata), intent(in) :: srcX,srcY,srcZ
      integer, intent(in), optional :: compX,compY,compZ
      integer, intent(in), optional :: nghost
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pM,pX,pY,pZ
      integer :: i,j,k,lvl,ng,cx,cy,cz
      cx=1; if (present(compX)) cx=compX
      cy=1; if (present(compY)) cy=compY
      cz=1; if (present(compZ)) cz=compZ
      ng=min(this%ng,srcX%ng,srcY%ng,srcZ%ng); if (present(nghost)) ng=nghost
      do lvl=0,this%amr%clvl()
         call amrex_mfiter_build(mfi,this%mf(lvl),tiling=.true.)
         do while (mfi%next())
            ! Get pointers to the data
            pM=>this%mf(lvl)%dataptr(mfi)
            pX=>srcX%mf(lvl)%dataptr(mfi)
            pY=>srcY%mf(lvl)%dataptr(mfi)
            pZ=>srcZ%mf(lvl)%dataptr(mfi)
            ! Loop over grown tile
            bx=mfi%growntilebox(ng)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pM(i,j,k,1)=sqrt(pX(i,j,k,cx)**2+pY(i,j,k,cy)**2+pZ(i,j,k,cz)**2)
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end do
   end subroutine get_magnitude

   ! ============================================================================
   ! CLONING
   ! ============================================================================

   !> Clone this amrdata into dest, optionally using alternate DMs per level
   !> All class members are fully populated. MultiFabs are built and parallel_copied.
   subroutine clone(this,dest,dm)
      use amrex_amr_module, only: amrex_multifab_build
      use messager, only: die
      implicit none
      class(amrdata), intent(in) :: this
      type(amrdata), intent(inout) :: dest
      type(amrex_distromap), dimension(0:), intent(in), optional :: dm
      integer :: lvl
      ! Copy all metadata
      dest%amr   =>this%amr
      dest%parent=>this%parent
      dest%name  =trim(this%name)//'_clone'
      dest%ncomp =this%ncomp
      dest%ng    =this%ng
      dest%nodal =this%nodal
      dest%interp=this%interp
      dest%fill_lvl_cache=this%fill_lvl_cache
      ! Copy BCs
      if (allocated(this%lo_bc)) allocate(dest%lo_bc,source=this%lo_bc)
      if (allocated(this%hi_bc)) allocate(dest%hi_bc,source=this%hi_bc)
      ! Copy callbacks
      dest%on_init  =>this%on_init
      dest%on_coarse=>this%on_coarse
      dest%on_remake=>this%on_remake
      dest%on_clear =>this%on_clear
      dest%fillbc   =>this%fillbc
      dest%user_init=>this%user_init
      ! Build MultiFabs on alternate or same DM, parallel_copy data
      if (present(dm)) then
         if (size(dm).lt.size(this%mf)) call die('[amrdata clone] dm array too small for number of levels')
      end if
      allocate(dest%mf(lbound(this%mf,1):ubound(this%mf,1)))
      do lvl=lbound(this%mf,1),ubound(this%mf,1)
         if (present(dm)) then; call amrex_multifab_build(dest%mf(lvl),this%amr%ba(lvl),         dm(lvl),this%ncomp,this%ng,this%nodal)
         else;                  call amrex_multifab_build(dest%mf(lvl),this%amr%ba(lvl),this%amr%dm(lvl),this%ncomp,this%ng,this%nodal)
         end if
         call dest%mf(lvl)%parallel_copy(this%mf(lvl),1,1,this%ncomp,this%ng,this%ng,this%amr%geom(lvl))
      end do
   end subroutine clone

   ! ============================================================================
   ! HELPER ROUTINES
   ! ============================================================================

   !> Compute level range from optional lvl and lbase arguments
   subroutine get_level_range(this, lvl, lbase, l0, l1)
      implicit none
      class(amrdata), intent(in) :: this
      integer, intent(in), optional :: lvl, lbase
      integer, intent(out) :: l0, l1
      if (present(lvl)) then
         l0 = lvl; l1 = lvl
      else if (present(lbase)) then
         l0 = lbase; l1 = this%amr%clvl()
      else
         l0 = 0; l1 = this%amr%clvl()
      end if
   end subroutine get_level_range

end module amrdata_class
