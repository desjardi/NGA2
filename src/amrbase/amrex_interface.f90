!> Centralized C interfaces for all AMReX wrappers in NGA2
!> These wrap the C++ functions in amrex_interface.cpp
!> Provides: AmrCore lifecycle, dispatchers, queries, and FillPatch operations
module amrex_interface
   use iso_c_binding, only: c_ptr,c_funptr,c_double,c_int
   use precision !< Force dep.py to track this module
   implicit none
   private

   !=====================================================================
   ! AmrCore Lifecycle
   !=====================================================================
   public :: amrcore_create
   public :: amrcore_destroy
   public :: amrcore_set_owner

   !=====================================================================
   ! AmrCore Dispatcher Setters
   !=====================================================================
   public :: amrcore_set_on_init_dispatch
   public :: amrcore_set_on_coarse_dispatch
   public :: amrcore_set_on_remake_dispatch
   public :: amrcore_set_on_clear_dispatch
   public :: amrcore_set_on_tag_dispatch
   public :: amrcore_set_on_postregrid_dispatch
   public :: amrcore_set_on_cost_dispatch
   public :: amrcore_set_cost_strategy
   public :: amrdm_make_knapsack,amrdm_make_sfc

   !=====================================================================
   ! AmrCore Grid Operations
   !=====================================================================
   public :: amrcore_init_from_scratch
   public :: amrcore_regrid

   !=====================================================================
   ! AmrCore Queries
   !=====================================================================
   public :: amrcore_finest_level
   public :: amrcore_get_geometry
   public :: amrcore_get_ref_ratio
   public :: amrcore_get_ref_ratio_xyz
   public :: amrcore_get_boxarray
   public :: amrcore_get_distromap

   !=====================================================================
   ! FillPatch Operations
   !=====================================================================
   public :: amrmfab_fillpatch_single
   public :: amrmfab_fillpatch_two
   public :: amrmfab_fillcoarsepatch
   public :: amrmfab_fillcoarsepatch_faces  ! 3-component
   public :: amrmfab_fillpatch_two_faces    ! 3-component

   !=====================================================================
   ! Plotfile Output
   !=====================================================================
   public :: amrplotfile_write_native
   public :: amrplotfile_write_hdf5
   public :: amrplotfile_read_time

   !=====================================================================
   ! Checkpoint I/O (VisMF)
   !=====================================================================
   public :: amrmfab_vismf_write
   public :: amrmfab_vismf_read
   public :: amrcheckpoint_prebuild_dirs
   public :: amrcheckpoint_mfab_prefix
   public :: amrvismf_set_noutfiles
   public :: amrvismf_get_noutfiles
   public :: amrcore_build_level
   public :: amrcore_set_finest_level

   !=====================================================================
   ! MLMG Utilities
   !=====================================================================
   public :: amrmlmg_get_niters
   public :: amrmlmg_get_fluxes
   public :: amrmlmg_dot_composite
   public :: amrpoisson_build
   public :: amrabeclap_build

   !=====================================================================
   ! MultiFab Averaging (unified API - not in AMReX Fortran interface)
   !=====================================================================
   public :: amrmfab_average_down_cell  ! Cell-centered
   public :: amrmfab_average_down_face  ! Face-centered (nodal in 1 dir)
   public :: amrmfab_average_down_edge  ! Edge-centered (nodal in 2 dirs)
   public :: amrmfab_average_down_node  ! Node-centered (nodal in 3 dirs)
   public :: amrmfab_sum_downto         ! Restrict-SUM fine deposits into coarse level (lvl+1 → lvl)
   public :: amrmfab_interp_from_coarse ! PCInterp coarse→fine for deposit processing
   public :: amrmfab_compute_divergence ! Compute div(u) from face velocities
   public :: amrmfab_sum_unique         ! Sum for face/nodal data (no double-counting)
   public :: amrmask_make_fine          ! Create mask for cells covered by finer level

   !=====================================================================
   ! Per-direction wrappers (bypass AMReX scalar-only Fortran interfaces)
   !=====================================================================
   public :: amrfluxreg_build           ! FluxRegister build with IntVect ratio
   public :: amrfluxreg_destroy         ! FluxRegister destroy
   public :: amrmfab_average_down_faces ! average_down_faces with IntVect ratio
   public :: amrlinop_set_coarse_fine_bc ! linop set_coarse_fine_bc with IntVect ratio

   interface

      !------------------------------------------------------------------
      ! AmrCore Lifecycle
      !------------------------------------------------------------------

      !> Create an NGA2AmrCore object
      type(c_ptr) function amrcore_create() bind(c)
         import :: c_ptr
      end function amrcore_create

      !> Destroy an NGA2AmrCore object
      subroutine amrcore_destroy(core) bind(c)
         import :: c_ptr
         type(c_ptr), value :: core
      end subroutine amrcore_destroy

      !> Set owner pointer for callbacks
      subroutine amrcore_set_owner(core,owner) bind(c)
         import :: c_ptr
         type(c_ptr), value :: core,owner
      end subroutine amrcore_set_owner

      !------------------------------------------------------------------
      ! AmrCore Dispatcher Setters
      !------------------------------------------------------------------

      !> Set MakeNewLevelFromScratch dispatcher
      subroutine amrcore_set_on_init_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_init_dispatch

      !> Set MakeNewLevelFromCoarse dispatcher
      subroutine amrcore_set_on_coarse_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_coarse_dispatch

      !> Set RemakeLevel dispatcher
      subroutine amrcore_set_on_remake_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_remake_dispatch

      !> Set ClearLevel dispatcher
      subroutine amrcore_set_on_clear_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_clear_dispatch

      !> Set ErrorEst (tagging) dispatcher
      subroutine amrcore_set_on_tag_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_tag_dispatch

      !> Set post-regrid dispatcher
      subroutine amrcore_set_on_postregrid_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_postregrid_dispatch

      !> Set cost callback dispatcher for load balancing
      subroutine amrcore_set_on_cost_dispatch(core,f) bind(c)
         import :: c_ptr,c_funptr
         type(c_ptr), value :: core
         type(c_funptr), value :: f
      end subroutine amrcore_set_on_cost_dispatch

      !> Set load balancing strategy (0=SFC, 1=KnapSack)
      subroutine amrcore_set_cost_strategy(core,strat) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: core
         integer(c_int), value :: strat
      end subroutine amrcore_set_cost_strategy

      !> Build a KnapSack-optimized DistributionMapping from per-box costs
      subroutine amrdm_make_knapsack_c(dm,costs,nboxes) bind(c,name='amrdm_make_knapsack')
         import :: c_ptr,c_int,c_double
         type(c_ptr) :: dm
         real(c_double) :: costs(*)
         integer(c_int), value :: nboxes
      end subroutine amrdm_make_knapsack_c

      !> Build an SFC-optimized DistributionMapping from per-box costs + BoxArray
      subroutine amrdm_make_sfc_c(dm,costs,nboxes,ba) bind(c,name='amrdm_make_sfc')
         import :: c_ptr,c_int,c_double
         type(c_ptr) :: dm
         real(c_double) :: costs(*)
         integer(c_int), value :: nboxes
         type(c_ptr), value :: ba
      end subroutine amrdm_make_sfc_c

      !> Destroy a heap-allocated DistributionMapping
      subroutine amrdm_destroy_c(dm) bind(c,name='amrdm_destroy')
         import :: c_ptr
         type(c_ptr), value :: dm
      end subroutine amrdm_destroy_c

      !------------------------------------------------------------------
      ! AmrCore Grid Operations
      !------------------------------------------------------------------

      !> Initialize grid from scratch
      subroutine amrcore_init_from_scratch(core,t) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: core
         real(c_double), value :: t
      end subroutine amrcore_init_from_scratch

      !> Perform regrid operation
      subroutine amrcore_regrid(core,blvl,t) bind(c)
         import :: c_ptr,c_int,c_double
         type(c_ptr), value :: core
         integer(c_int), value :: blvl
         real(c_double), value :: t
      end subroutine amrcore_regrid

      !------------------------------------------------------------------
      ! AmrCore Queries
      !------------------------------------------------------------------

      !> Get current finest level
      integer(c_int) function amrcore_finest_level(core) bind(c)
         import :: c_int,c_ptr
         type(c_ptr), value :: core
      end function amrcore_finest_level

      !> Get geometry at a level
      subroutine amrcore_get_geometry(geom,lvl,core) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), intent(out) :: geom
         integer(c_int), value :: lvl
         type(c_ptr), value :: core
      end subroutine amrcore_get_geometry

      !> Get refinement ratios
      subroutine amrcore_get_ref_ratio(ref_ratio,core) bind(c)
         import :: c_ptr
         integer, dimension(*), intent(inout) :: ref_ratio
         type(c_ptr), value :: core
      end subroutine amrcore_get_ref_ratio

      !> Get per-direction refinement ratios
      subroutine amrcore_get_ref_ratio_xyz(rrefx,rrefy,rrefz,core) bind(c)
         import :: c_ptr
         integer, dimension(*), intent(inout) :: rrefx,rrefy,rrefz
         type(c_ptr), value :: core
      end subroutine amrcore_get_ref_ratio_xyz

      !> Get BoxArray at a level
      subroutine amrcore_get_boxarray(barray,lev,core) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), intent(out) :: barray
         integer(c_int), value :: lev
         type(c_ptr), value :: core
      end subroutine amrcore_get_boxarray

      !> Get DistributionMapping at a level
      subroutine amrcore_get_distromap(dmap,lev,core) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), intent(out) :: dmap
         integer(c_int), value :: lev
         type(c_ptr), value :: core
      end subroutine amrcore_get_distromap

      !------------------------------------------------------------------
      ! FillPatch Operations
      !------------------------------------------------------------------

      !> FillPatch for level 0 (single level, physical BCs only)
      subroutine amrmfab_fillpatch_single_c(mf,t_old,mf_old,t_new,mf_new, &
      &   geom,solver_ctx,bc_dispatch,time,scomp,dcomp,ncomp) &
      &   bind(c, name='amrmfab_fillpatch_single')
         import :: c_ptr,c_funptr,c_double,c_int
         type(c_ptr), value :: mf,mf_old,mf_new,geom,solver_ctx
         type(c_funptr), value :: bc_dispatch
         real(c_double), value :: t_old,t_new,time
         integer(c_int), value :: scomp,dcomp,ncomp
      end subroutine amrmfab_fillpatch_single_c

      !> FillPatch for fine levels (two-level interpolation + BCs)
      subroutine amrmfab_fillpatch_two_c(mf,t_old_c,mf_old_c,t_new_c,mf_new_c,geom_c, &
      &   t_old_f,mf_old_f,t_new_f,mf_new_f,geom_f,solver_ctx,bc_dispatch, &
      &   time,scomp,dcomp,ncomp,ref_ratio,interp_type,lo_bc,hi_bc,nbc) &
      &   bind(c, name='amrmfab_fillpatch_two')
         import :: c_ptr,c_funptr,c_double,c_int
         type(c_ptr), value :: mf,mf_old_c,mf_new_c,geom_c
         type(c_ptr), value :: mf_old_f,mf_new_f,geom_f,solver_ctx
         type(c_funptr), value :: bc_dispatch
         real(c_double), value :: t_old_c,t_new_c,t_old_f,t_new_f,time
         integer(c_int), value :: scomp,dcomp,ncomp,interp_type,nbc
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), intent(in) :: lo_bc(*),hi_bc(*)
      end subroutine amrmfab_fillpatch_two_c

      !> FillCoarsePatch - fill fine level from coarse only (for new levels)
      subroutine amrmfab_fillcoarsepatch_c(mf_f,time,mf_c,geom_c,geom_f,solver_ctx, &
      &   bc_dispatch,scomp,dcomp,ncomp,ref_ratio,interp_type,lo_bc,hi_bc,nbc) &
      &   bind(c, name='amrmfab_fillcoarsepatch')
         import :: c_ptr,c_funptr,c_double,c_int
         type(c_ptr), value :: mf_f,mf_c,geom_c,geom_f,solver_ctx
         type(c_funptr), value :: bc_dispatch
         real(c_double), value :: time
         integer(c_int), value :: scomp,dcomp,ncomp,interp_type,nbc
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), intent(in) :: lo_bc(*),hi_bc(*)
      end subroutine amrmfab_fillcoarsepatch_c


      !> FillCoarsePatch for 3-component face-centered velocity
      subroutine amrmfab_fillcoarsepatch_faces_c(mf_u, mf_v, mf_w, time, &
      &   cmf_u, cmf_v, cmf_w, geom_c, geom_f, &
      &   ctx_u, ctx_v, ctx_w, bc_u, bc_v, bc_w, &
      &   scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc) &
      &   bind(c, name='amrmfab_fillcoarsepatch_faces')
         import :: c_ptr, c_funptr, c_double, c_int
         type(c_ptr), value :: mf_u, mf_v, mf_w
         type(c_ptr), value :: cmf_u, cmf_v, cmf_w
         type(c_ptr), value :: geom_c, geom_f
         type(c_ptr), value :: ctx_u, ctx_v, ctx_w
         type(c_funptr), value :: bc_u, bc_v, bc_w
         real(c_double), value :: time
         integer(c_int), value :: scomp, dcomp, ncomp, interp_type
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), intent(in) :: lo_bc(*), hi_bc(*)
      end subroutine amrmfab_fillcoarsepatch_faces_c

      !> FillPatch for 3-component face-centered velocity from two levels
      subroutine amrmfab_fillpatch_two_faces_c(mf_u, mf_v, mf_w, time, &
      &   t_old_c, mf_old_c_u, mf_old_c_v, mf_old_c_w, &
      &   t_new_c, mf_new_c_u, mf_new_c_v, mf_new_c_w, geom_c, &
      &   t_old_f, mf_old_f_u, mf_old_f_v, mf_old_f_w, &
      &   t_new_f, mf_new_f_u, mf_new_f_v, mf_new_f_w, geom_f, &
      &   ctx_u, ctx_v, ctx_w, bc_u, bc_v, bc_w, &
      &   scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc) &
      &   bind(c, name='amrmfab_fillpatch_two_faces')
         import :: c_ptr, c_funptr, c_double, c_int
         type(c_ptr), value :: mf_u, mf_v, mf_w
         type(c_ptr), value :: mf_old_c_u, mf_old_c_v, mf_old_c_w
         type(c_ptr), value :: mf_new_c_u, mf_new_c_v, mf_new_c_w
         type(c_ptr), value :: mf_old_f_u, mf_old_f_v, mf_old_f_w
         type(c_ptr), value :: mf_new_f_u, mf_new_f_v, mf_new_f_w
         type(c_ptr), value :: geom_c, geom_f
         type(c_ptr), value :: ctx_u, ctx_v, ctx_w
         type(c_funptr), value :: bc_u, bc_v, bc_w
         real(c_double), value :: time, t_old_c, t_new_c, t_old_f, t_new_f
         integer(c_int), value :: scomp, dcomp, ncomp, interp_type
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), intent(in) :: lo_bc(*), hi_bc(*)
      end subroutine amrmfab_fillpatch_two_faces_c

      !====================================================================
      ! Plotfile Output
      !====================================================================

      !> Write multi-level plotfile in native AMReX format
      subroutine amrplotfile_write_native(name, nlevels, mf_ptrs, varnames, ncomp, &
      &   geom_ptrs, time, level_steps, ref_ratios) bind(c)
         import :: c_ptr,c_double,c_int
         character(kind=1), intent(in) :: name(*)
         integer(c_int), value :: nlevels, ncomp
         type(c_ptr), intent(in) :: mf_ptrs(*), geom_ptrs(*)
         type(c_ptr), intent(in) :: varnames(*)
         real(c_double), value :: time
         integer(c_int), intent(in) :: level_steps(*), ref_ratios(*)
      end subroutine amrplotfile_write_native

      !> Write multi-level plotfile in HDF5 format (if available)
      subroutine amrplotfile_write_hdf5(name, nlevels, mf_ptrs, varnames, ncomp, &
      &   geom_ptrs, time, level_steps, ref_ratios, compression) bind(c)
         import :: c_ptr,c_double,c_int
         character(kind=1), intent(in) :: name(*)
         integer(c_int), value :: nlevels, ncomp
         type(c_ptr), intent(in) :: mf_ptrs(*), geom_ptrs(*)
         type(c_ptr), intent(in) :: varnames(*)
         real(c_double), value :: time
         integer(c_int), intent(in) :: level_steps(*), ref_ratios(*)
         character(kind=1), intent(in) :: compression(*)
      end subroutine amrplotfile_write_hdf5

      !> Read time attribute from HDF5 plotfile (returns -1 if file not found)
      function amrplotfile_read_time(filename) result(time) bind(c)
         import :: c_double
         character(kind=1), intent(in) :: filename(*)
         real(c_double) :: time
      end function amrplotfile_read_time

      !====================================================================
      ! Checkpoint I/O (VisMF)
      !====================================================================

      !> Write MultiFab to disk using VisMF
      subroutine amrmfab_vismf_write(mf, path) bind(c)
         import :: c_ptr
         type(c_ptr), value :: mf
         character(kind=1), intent(in) :: path(*)
      end subroutine amrmfab_vismf_write

      !> Read MultiFab from disk using VisMF
      subroutine amrmfab_vismf_read(mf, path) bind(c)
         import :: c_ptr
         type(c_ptr), value :: mf
         character(kind=1), intent(in) :: path(*)
      end subroutine amrmfab_vismf_read

      !> Create directory hierarchy for checkpoints
      subroutine amrcheckpoint_prebuild_dirs(dirname, subdir_prefix, nlevels) bind(c)
         import :: c_int
         character(kind=1), intent(in) :: dirname(*), subdir_prefix(*)
         integer(c_int), value :: nlevels
      end subroutine amrcheckpoint_prebuild_dirs

      !> Get MultiFab file prefix for a level
      subroutine amrcheckpoint_mfab_prefix(result, result_len, lev, dirname, &
      &   level_prefix, mfab_name) bind(c)
         import :: c_int
         character(kind=1), intent(out) :: result(*)
         integer(c_int), value :: result_len, lev
         character(kind=1), intent(in) :: dirname(*), level_prefix(*), mfab_name(*)
      end subroutine amrcheckpoint_mfab_prefix

      !> Set number of output files for VisMF (controls I/O aggregation)
      subroutine amrvismf_set_noutfiles(nfiles) bind(c)
         import :: c_int
         integer(c_int), value :: nfiles
      end subroutine amrvismf_set_noutfiles

      !> Get current number of output files setting
      integer(c_int) function amrvismf_get_noutfiles() bind(c)
         import :: c_int
      end function amrvismf_get_noutfiles

      !> Build a single level from pre-built BoxArray and DistributionMapping
      subroutine amrcore_build_level(core, lev, time, ba_ptr, dm_ptr) bind(c)
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: core
         integer(c_int), value :: lev
         real(c_double), value :: time
         type(c_ptr), value :: ba_ptr
         type(c_ptr), value :: dm_ptr
      end subroutine amrcore_build_level

      !> Set finest_level on AmrCore
      subroutine amrcore_set_finest_level(core, finest_level) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: core
         integer(c_int), value :: finest_level
      end subroutine amrcore_set_finest_level

      !====================================================================
      ! MLMG Utilities (not available in AMReX Fortran interface)
      !====================================================================

      !> Get number of iterations from last MLMG solve
      integer(c_int) function amrmlmg_get_niters(mlmg) bind(c)
         import :: c_int,c_ptr
         type(c_ptr), value :: mlmg
      end function amrmlmg_get_niters

      !> Get face-centered fluxes from MLMG using solver's C/F stencils
      !> For (alpha*A - beta*div(B*grad))phi = rhs, flux = -B*grad(phi)
      subroutine amrmlmg_get_fluxes(mlmg, sol_mfs, flux_x, flux_y, flux_z, nlevs) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: mlmg
         type(c_ptr), intent(in) :: sol_mfs(*)
         type(c_ptr), intent(in) :: flux_x(*), flux_y(*), flux_z(*)
         integer(c_int), value :: nlevs
      end subroutine amrmlmg_get_fluxes

      !> Composite masked dot product over all AMR levels
      !> Excludes coarse cells covered by finer levels (no double counting)
      !>   mf1_ptrs, mf2_ptrs - arrays of MultiFab pointers (one per level)
      !>   ba_fine_ptrs       - BoxArray pointers (ba_fine_ptrs(lev+1) used for mask at lev)
      !>   ref_ratios         - flattened (rx,ry,rz)*nlevs array
      !>   nlevs              - total number of AMR levels
      real(c_double) function amrmlmg_dot_composite(mf1_ptrs, mf2_ptrs, ba_fine_ptrs, ref_ratios, nlevs) bind(c)
         import :: c_ptr, c_int, c_double
         type(c_ptr), intent(in) :: mf1_ptrs(*), mf2_ptrs(*), ba_fine_ptrs(*)
         integer(c_int), intent(in) :: ref_ratios(*)
         integer(c_int), value :: nlevs
      end function amrmlmg_dot_composite

      !====================================================================
      ! Per-direction wrappers (bypass AMReX scalar-only Fortran interfaces)
      !====================================================================

      subroutine amrfluxreg_build(fr_ptr, ba_ptr, dm_ptr, ref_ratio, fine_lev, ncomp) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), intent(out) :: fr_ptr
         type(c_ptr), value :: ba_ptr, dm_ptr
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), value :: fine_lev, ncomp
      end subroutine amrfluxreg_build

      subroutine amrfluxreg_destroy(fr_ptr) bind(c)
         import :: c_ptr
         type(c_ptr), value :: fr_ptr
      end subroutine amrfluxreg_destroy

      subroutine amrmfab_average_down_faces(fine_x, fine_y, fine_z, &
         crse_x, crse_y, crse_z, geom, scomp, ncomp, ref_ratio) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_x, fine_y, fine_z
         type(c_ptr), value :: crse_x, crse_y, crse_z
         type(c_ptr), value :: geom
         integer(c_int), value :: scomp, ncomp
         integer(c_int), intent(in) :: ref_ratio(3)
      end subroutine amrmfab_average_down_faces

      subroutine amrlinop_set_coarse_fine_bc(linop, crse_mf, ref_ratio) bind(c)
         import :: c_ptr, c_int
         type(c_ptr), value :: linop, crse_mf
         integer(c_int), intent(in) :: ref_ratio(3)
      end subroutine amrlinop_set_coarse_fine_bc

      !> Build MLPoisson
      subroutine amrpoisson_build_c(linop, nlevels, geom, ba, dm, semicoarsen, hidden_direction) bind(c, name='amrpoisson_build_c')
         import :: c_ptr, c_int
         type(c_ptr), intent(out) :: linop
         integer(c_int), value :: nlevels, semicoarsen, hidden_direction
         type(c_ptr), intent(in) :: geom(*), ba(*), dm(*)
      end subroutine amrpoisson_build_c

      !> Build MLABecLaplacian
      subroutine amrabeclap_build_c(linop, nlevels, geom, ba, dm, semicoarsen, hidden_direction) bind(c, name='amrabeclap_build_c')
         import :: c_ptr, c_int
         type(c_ptr), intent(out) :: linop
         integer(c_int), value :: nlevels, semicoarsen, hidden_direction
         type(c_ptr), intent(in) :: geom(*), ba(*), dm(*)
      end subroutine amrabeclap_build_c

      !====================================================================
      ! MultiFab Averaging (C bindings - unified API for all 4 types)
      ! Signature: (fine, crse, geom, ratio, ngcrse)
      ! geom can be null (c_null_ptr) to skip periodic fix-up
      !====================================================================

      subroutine amrmfab_average_down_cell_c(fine_mf, crse_mf, crse_geom, ref_ratio, ngcrse) &
         bind(c, name='amrmfab_average_down_cell')
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_mf, crse_mf, crse_geom
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), value :: ngcrse
      end subroutine amrmfab_average_down_cell_c

      subroutine amrmfab_average_down_face_c(fine_mf, crse_mf, crse_geom, ref_ratio, ngcrse) &
         bind(c, name='amrmfab_average_down_face')
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_mf, crse_mf, crse_geom
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), value :: ngcrse
      end subroutine amrmfab_average_down_face_c

      subroutine amrmfab_average_down_edge_c(fine_mf, crse_mf, crse_geom, ref_ratio, ngcrse) &
         bind(c, name='amrmfab_average_down_edge')
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_mf, crse_mf, crse_geom
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), value :: ngcrse
      end subroutine amrmfab_average_down_edge_c

      subroutine amrmfab_average_down_node_c(fine_mf, crse_mf, crse_geom, ref_ratio, ngcrse) &
         bind(c, name='amrmfab_average_down_node')
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_mf, crse_mf, crse_geom
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), value :: ngcrse
      end subroutine amrmfab_average_down_node_c

      !> Restrict-SUM all fine deposits into the coarse level (lvl+1 → lvl)
      subroutine amrmfab_sum_downto_c(fine_mf, crse_mf, crse_geom, fine_geom, ref_ratio) &
         bind(c, name='amrmfab_sum_downto')
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_mf, crse_mf, crse_geom, fine_geom
         integer(c_int), intent(in) :: ref_ratio(3)
      end subroutine amrmfab_sum_downto_c

      !> PCInterp coarse→fine for deposit processing (no-op BCs, 0-indexed scomp)
      subroutine amrmfab_interp_from_coarse_c(fine_mf, crse_mf, crse_geom, fine_geom, scomp, ncomp, ref_ratio) &
         bind(c, name='amrmfab_interp_from_coarse')
         import :: c_ptr, c_int
         type(c_ptr), value :: fine_mf, crse_mf, crse_geom, fine_geom
         integer(c_int), value :: scomp, ncomp
         integer(c_int), intent(in) :: ref_ratio(3)
      end subroutine amrmfab_interp_from_coarse_c

      !> Compute divergence of face-centered velocity into cell-centered MultiFab
      subroutine amrmfab_compute_divergence_c(divu, umac_x, umac_y, umac_z, geom) &
         bind(c, name='amrmfab_compute_divergence')
         import :: c_ptr
         type(c_ptr), value :: divu, umac_x, umac_y, umac_z, geom
      end subroutine amrmfab_compute_divergence_c

      !> Sum MultiFab avoiding double-counting at shared nodes/faces (for periodic)
      real(c_double) function amrmfab_sum_unique_c(mf, geom, comp) &
         bind(c, name='amrmfab_sum_unique')
         import :: c_ptr, c_double, c_int
         type(c_ptr), value :: mf, geom
         integer(c_int), value :: comp
      end function amrmfab_sum_unique_c

      !> Create mask identifying cells covered by finer level
      subroutine amrmask_make_fine_c(mask, ba_fine, ref_ratio, covered_val, notcovered_val) &
         bind(c, name='amrmask_make_fine')
         import :: c_ptr, c_int
         type(c_ptr), value :: mask, ba_fine
         integer(c_int), intent(in) :: ref_ratio(3)
         integer(c_int), value :: covered_val, notcovered_val
      end subroutine amrmask_make_fine_c

   end interface

contains

   !====================================================================
   ! MultiFab Averaging Wrappers (Unified API)
   ! All 4 types have identical signature:
   !   (fmf, cmf, rr, cgeom, ngcrse)
   ! cgeom and ngcrse are optional. If cgeom provided, FillBoundary is
   ! called in C++ for periodic ghost fix-up.
   !====================================================================

   !> Average down cell-centered MultiFab
   subroutine amrmfab_average_down_cell(fmf, cmf, rr, cgeom, ngcrse)
      use iso_c_binding, only: c_null_ptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(in) :: fmf
      type(amrex_multifab), intent(inout) :: cmf
      integer, intent(in) :: rr(3)
      type(amrex_geometry), intent(in), optional :: cgeom
      integer, intent(in), optional :: ngcrse
      integer :: ng
      ng = 0; if (present(ngcrse)) ng = ngcrse
      if (present(cgeom)) then
         call amrmfab_average_down_cell_c(fmf%p, cmf%p, cgeom%p, rr, ng)
      else
         call amrmfab_average_down_cell_c(fmf%p, cmf%p, c_null_ptr, rr, ng)
      end if
   end subroutine amrmfab_average_down_cell

   !> Average down face-centered MultiFab (nodal in 1 dir, cell in 2)
   subroutine amrmfab_average_down_face(fmf, cmf, rr, cgeom, ngcrse)
      use iso_c_binding, only: c_null_ptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(in) :: fmf
      type(amrex_multifab), intent(inout) :: cmf
      integer, intent(in) :: rr(3)
      type(amrex_geometry), intent(in), optional :: cgeom
      integer, intent(in), optional :: ngcrse
      integer :: ng
      ng = 0; if (present(ngcrse)) ng = ngcrse
      if (present(cgeom)) then
         call amrmfab_average_down_face_c(fmf%p, cmf%p, cgeom%p, rr, ng)
      else
         call amrmfab_average_down_face_c(fmf%p, cmf%p, c_null_ptr, rr, ng)
      end if
   end subroutine amrmfab_average_down_face

   !> Average down edge-centered MultiFab (nodal in 2 dirs, cell in 1)
   subroutine amrmfab_average_down_edge(fmf, cmf, rr, cgeom, ngcrse)
      use iso_c_binding, only: c_null_ptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(in) :: fmf
      type(amrex_multifab), intent(inout) :: cmf
      integer, intent(in) :: rr(3)
      type(amrex_geometry), intent(in), optional :: cgeom
      integer, intent(in), optional :: ngcrse
      integer :: ng
      ng = 0; if (present(ngcrse)) ng = ngcrse
      if (present(cgeom)) then
         call amrmfab_average_down_edge_c(fmf%p, cmf%p, cgeom%p, rr, ng)
      else
         call amrmfab_average_down_edge_c(fmf%p, cmf%p, c_null_ptr, rr, ng)
      end if
   end subroutine amrmfab_average_down_edge

   !> Average down node-centered MultiFab (nodal in all dirs)
   subroutine amrmfab_average_down_node(fmf, cmf, rr, cgeom, ngcrse)
      use iso_c_binding, only: c_null_ptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(in) :: fmf
      type(amrex_multifab), intent(inout) :: cmf
      integer, intent(in) :: rr(3)
      type(amrex_geometry), intent(in), optional :: cgeom
      integer, intent(in), optional :: ngcrse
      integer :: ng
      ng = 0; if (present(ngcrse)) ng = ngcrse
      if (present(cgeom)) then
         call amrmfab_average_down_node_c(fmf%p, cmf%p, cgeom%p, rr, ng)
      else
         call amrmfab_average_down_node_c(fmf%p, cmf%p, c_null_ptr, rr, ng)
      end if
   end subroutine amrmfab_average_down_node

   !> Restrict-SUM fine deposits (valid+ghost) into the coarse level.
   !> Delegates to AMReX's sum_fine_to_coarse; requires nGrow % ratio == 0.
   !> Call after SumBoundary; average_down afterward fixes double-counted cells.
   subroutine amrmfab_sum_downto(fmf,cmf,rr,cgeom,fgeom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(in)    :: fmf
      type(amrex_multifab), intent(inout) :: cmf
      integer,              intent(in)    :: rr(3)
      type(amrex_geometry), intent(in)    :: cgeom
      type(amrex_geometry), intent(in)    :: fgeom
      call amrmfab_sum_downto_c(fmf%p, cmf%p, cgeom%p, fgeom%p, rr)
   end subroutine amrmfab_sum_downto

   !> PCInterp coarse→fine interpolation for deposit processing.
   !> Propagates coarse-level deposits into fine-covered cells (no-op BCs).
   subroutine amrmfab_interp_from_coarse(fmf,cmf,rr,cgeom,fgeom,scomp,ncomp)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: fmf
      type(amrex_multifab), intent(in)    :: cmf
      integer,              intent(in)    :: rr(3)
      type(amrex_geometry), intent(in)    :: cgeom
      type(amrex_geometry), intent(in)    :: fgeom
      integer,              intent(in)    :: scomp,ncomp
      call amrmfab_interp_from_coarse_c(fmf%p, cmf%p, cgeom%p, fgeom%p, scomp-1, ncomp, rr)
   end subroutine amrmfab_interp_from_coarse

   !> Compute divergence of face-centered velocity into cell-centered MultiFab
   subroutine amrmfab_compute_divergence(divu, umac_x, umac_y, umac_z, geom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: divu
      type(amrex_multifab), intent(in) :: umac_x, umac_y, umac_z
      type(amrex_geometry), intent(in) :: geom
      call amrmfab_compute_divergence_c(divu%p, umac_x%p, umac_y%p, umac_z%p, geom%p)
   end subroutine amrmfab_compute_divergence

   !> Sum MultiFab avoiding double-counting at shared nodes/faces (for periodic)
   function amrmfab_sum_unique(mf, geom, comp) result(val)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      use precision, only: WP
      type(amrex_multifab), intent(in) :: mf
      type(amrex_geometry), intent(in) :: geom
      integer, intent(in), optional :: comp
      real(WP) :: val
      integer :: ic
      ic = 1; if (present(comp)) ic = comp
      val = amrmfab_sum_unique_c(mf%p, geom%p, ic)
   end function amrmfab_sum_unique

   !> Create integer mask identifying cells covered by finer level
   !> mask must already be built with same BA/DM as coarse level
   !> After call: mask(i,j,k) = notcovered_val where valid, covered_val where fine exists
   subroutine amrmask_make_fine(mask, ba_fine, ref_ratio, covered_val, notcovered_val)
      use amrex_amr_module, only: amrex_imultifab, amrex_boxarray
      type(amrex_imultifab), intent(inout) :: mask
      type(amrex_boxarray), intent(in) :: ba_fine
      integer, intent(in) :: ref_ratio(3)
      integer, intent(in) :: covered_val, notcovered_val
      call amrmask_make_fine_c(mask%p, ba_fine%p, ref_ratio, covered_val, notcovered_val)
   end subroutine amrmask_make_fine

   !> Fill coarse patch for 3-component face-centered velocity
   subroutine amrmfab_fillcoarsepatch_faces(mf_u, mf_v, mf_w, time, &
   &   cmf_u, cmf_v, cmf_w, geom_c, geom_f, &
   &   ctx_u, ctx_v, ctx_w, bc_u, bc_v, bc_w, &
   &   scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc)
      use iso_c_binding, only: c_ptr, c_funptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: mf_u, mf_v, mf_w
      type(amrex_multifab), intent(in) :: cmf_u, cmf_v, cmf_w
      type(amrex_geometry), intent(in) :: geom_c, geom_f
      type(c_ptr), intent(in) :: ctx_u, ctx_v, ctx_w
      type(c_funptr), intent(in) :: bc_u, bc_v, bc_w
      real(8), intent(in) :: time
      integer, intent(in) :: scomp, dcomp, ncomp, interp_type
      integer, intent(in) :: ref_ratio(3)
      integer, intent(in) :: lo_bc(*), hi_bc(*)
      call amrmfab_fillcoarsepatch_faces_c(mf_u%p, mf_v%p, mf_w%p, time, &
      &   cmf_u%p, cmf_v%p, cmf_w%p, geom_c%p, geom_f%p, &
      &   ctx_u, ctx_v, ctx_w, bc_u, bc_v, bc_w, &
      &   scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc)
   end subroutine amrmfab_fillcoarsepatch_faces

   !> FillPatch for 3-component face-centered velocity from two levels
   subroutine amrmfab_fillpatch_two_faces(mf_u, mf_v, mf_w, time, &
   &   t_old_c, mf_old_c_u, mf_old_c_v, mf_old_c_w, &
   &   t_new_c, mf_new_c_u, mf_new_c_v, mf_new_c_w, geom_c, &
   &   t_old_f, mf_old_f_u, mf_old_f_v, mf_old_f_w, &
   &   t_new_f, mf_new_f_u, mf_new_f_v, mf_new_f_w, geom_f, &
   &   ctx_u, ctx_v, ctx_w, bc_u, bc_v, bc_w, &
   &   scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc)
      use iso_c_binding, only: c_ptr, c_funptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: mf_u, mf_v, mf_w
      type(amrex_multifab), intent(in) :: mf_old_c_u, mf_old_c_v, mf_old_c_w
      type(amrex_multifab), intent(in) :: mf_new_c_u, mf_new_c_v, mf_new_c_w
      type(amrex_multifab), intent(in) :: mf_old_f_u, mf_old_f_v, mf_old_f_w
      type(amrex_multifab), intent(in) :: mf_new_f_u, mf_new_f_v, mf_new_f_w
      type(amrex_geometry), intent(in) :: geom_c, geom_f
      type(c_ptr), intent(in) :: ctx_u, ctx_v, ctx_w
      type(c_funptr), intent(in) :: bc_u, bc_v, bc_w
      real(8), intent(in) :: time, t_old_c, t_new_c, t_old_f, t_new_f
      integer, intent(in) :: scomp, dcomp, ncomp, interp_type
      integer, intent(in) :: ref_ratio(3)
      integer, intent(in) :: lo_bc(*), hi_bc(*)
      call amrmfab_fillpatch_two_faces_c(mf_u%p, mf_v%p, mf_w%p, time, &
      &   t_old_c, mf_old_c_u%p, mf_old_c_v%p, mf_old_c_w%p, &
      &   t_new_c, mf_new_c_u%p, mf_new_c_v%p, mf_new_c_w%p, geom_c%p, &
      &   t_old_f, mf_old_f_u%p, mf_old_f_v%p, mf_old_f_w%p, &
      &   t_new_f, mf_new_f_u%p, mf_new_f_v%p, mf_new_f_w%p, geom_f%p, &
      &   ctx_u, ctx_v, ctx_w, bc_u, bc_v, bc_w, &
      &   scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc)
   end subroutine amrmfab_fillpatch_two_faces

   !> FillPatch for level 0 (single level, physical BCs only)
   subroutine amrmfab_fillpatch_single(mf, t_old, mf_old, t_new, mf_new, &
   &   geom, solver_ctx, bc_dispatch, time, scomp, dcomp, ncomp)
      use iso_c_binding, only: c_ptr, c_funptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: mf
      type(amrex_multifab), intent(in) :: mf_old, mf_new
      type(amrex_geometry), intent(in) :: geom
      type(c_ptr), intent(in) :: solver_ctx
      type(c_funptr), intent(in) :: bc_dispatch
      real(8), intent(in) :: t_old, t_new, time
      integer, intent(in) :: scomp, dcomp, ncomp
      call amrmfab_fillpatch_single_c(mf%p, t_old, mf_old%p, t_new, mf_new%p, &
      &   geom%p, solver_ctx, bc_dispatch, time, scomp, dcomp, ncomp)
   end subroutine amrmfab_fillpatch_single

   !> FillPatch for fine levels (two-level interpolation + BCs)
   subroutine amrmfab_fillpatch_two(mf, t_old_c, mf_old_c, t_new_c, mf_new_c, geom_c, &
   &   t_old_f, mf_old_f, t_new_f, mf_new_f, geom_f, solver_ctx, bc_dispatch, &
   &   time, scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc, nbc)
      use iso_c_binding, only: c_ptr, c_funptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: mf
      type(amrex_multifab), intent(in) :: mf_old_c, mf_new_c, mf_old_f, mf_new_f
      type(amrex_geometry), intent(in) :: geom_c, geom_f
      type(c_ptr), intent(in) :: solver_ctx
      type(c_funptr), intent(in) :: bc_dispatch
      real(8), intent(in) :: t_old_c, t_new_c, t_old_f, t_new_f, time
      integer, intent(in) :: scomp, dcomp, ncomp, interp_type, nbc
      integer, intent(in) :: ref_ratio(3)
      integer, intent(in) :: lo_bc(*), hi_bc(*)
      call amrmfab_fillpatch_two_c(mf%p, t_old_c, mf_old_c%p, t_new_c, mf_new_c%p, &
      &   geom_c%p, t_old_f, mf_old_f%p, t_new_f, mf_new_f%p, geom_f%p, &
      &   solver_ctx, bc_dispatch, time, scomp, dcomp, ncomp, ref_ratio, &
      &   interp_type, lo_bc, hi_bc, nbc)
   end subroutine amrmfab_fillpatch_two

   !> FillCoarsePatch - fill fine level from coarse only (for new levels)
   subroutine amrmfab_fillcoarsepatch(mf_f, time, mf_c, geom_c, geom_f, solver_ctx, &
   &   bc_dispatch, scomp, dcomp, ncomp, ref_ratio, interp_type, lo_bc, hi_bc, nbc)
      use iso_c_binding, only: c_ptr, c_funptr
      use amrex_amr_module, only: amrex_multifab, amrex_geometry
      type(amrex_multifab), intent(inout) :: mf_f
      type(amrex_multifab), intent(in) :: mf_c
      type(amrex_geometry), intent(in) :: geom_c, geom_f
      type(c_ptr), intent(in) :: solver_ctx
      type(c_funptr), intent(in) :: bc_dispatch
      real(8), intent(in) :: time
      integer, intent(in) :: scomp, dcomp, ncomp, interp_type, nbc
      integer, intent(in) :: ref_ratio(3)
      integer, intent(in) :: lo_bc(*), hi_bc(*)
      call amrmfab_fillcoarsepatch_c(mf_f%p, time, mf_c%p, geom_c%p, geom_f%p, &
      &   solver_ctx, bc_dispatch, scomp, dcomp, ncomp, ref_ratio, interp_type, &
      &   lo_bc, hi_bc, nbc)
   end subroutine amrmfab_fillcoarsepatch

   !> Build MLPoisson operator
   subroutine amrpoisson_build(linop, nlevels, geom, ba, dm, semicoarsen, hidden_direction)
      use iso_c_binding, only: c_ptr
      use amrex_amr_module, only: amrex_geometry, amrex_boxarray, amrex_distromap
      use amrex_linear_solver_module, only: amrex_poisson
      type(amrex_poisson), intent(inout) :: linop
      integer, intent(in) :: nlevels
      type(amrex_geometry), intent(in) :: geom(0:)
      type(amrex_boxarray), intent(in) :: ba(0:)
      type(amrex_distromap), intent(in) :: dm(0:)
      logical, intent(in), optional :: semicoarsen
      integer, intent(in), optional :: hidden_direction
      type(c_ptr) :: gp(0:nlevels-1), bp(0:nlevels-1), dp(0:nlevels-1)
      integer :: lev, sc, hd
      sc = 0; if (present(semicoarsen)) sc = merge(1,0,semicoarsen)
      hd = -1; if (present(hidden_direction)) hd = hidden_direction
      do lev = 0, nlevels-1
         gp(lev) = geom(lev)%p
         bp(lev) = ba(lev)%p
         dp(lev) = dm(lev)%p
      end do
      linop%owner = .true.
      call amrpoisson_build_c(linop%p, nlevels, gp, bp, dp, sc, hd)
   end subroutine amrpoisson_build

   !> Build MLABecLaplacian operator
   subroutine amrabeclap_build(linop, nlevels, geom, ba, dm, semicoarsen, hidden_direction)
      use iso_c_binding, only: c_ptr
      use amrex_amr_module, only: amrex_geometry, amrex_boxarray, amrex_distromap
      use amrex_linear_solver_module, only: amrex_abeclaplacian
      type(amrex_abeclaplacian), intent(inout) :: linop
      integer, intent(in) :: nlevels
      type(amrex_geometry), intent(in) :: geom(0:)
      type(amrex_boxarray), intent(in) :: ba(0:)
      type(amrex_distromap), intent(in) :: dm(0:)
      logical, intent(in), optional :: semicoarsen
      integer, intent(in), optional :: hidden_direction
      type(c_ptr) :: gp(0:nlevels-1), bp(0:nlevels-1), dp(0:nlevels-1)
      integer :: lev, sc, hd
      sc = 0; if (present(semicoarsen)) sc = merge(1,0,semicoarsen)
      hd = -1; if (present(hidden_direction)) hd = hidden_direction
      do lev = 0, nlevels-1
         gp(lev) = geom(lev)%p
         bp(lev) = ba(lev)%p
         dp(lev) = dm(lev)%p
      end do
      linop%owner = .true.
      call amrabeclap_build_c(linop%p, nlevels, gp, bp, dp, sc, hd)
   end subroutine amrabeclap_build

   !====================================================================
   ! Cost-weighted DistributionMapping Factories
   !====================================================================

   !> Build a KnapSack-optimized DM from per-box costs
   subroutine amrdm_make_knapsack(dm,costs,nboxes)
      use amrex_amr_module, only: amrex_distromap
      type(amrex_distromap), intent(out) :: dm
      integer, intent(in) :: nboxes
      real(WP), intent(in) :: costs(nboxes)
      call amrdm_make_knapsack_c(dm%p,costs,nboxes)
      dm%owner=.true.
   end subroutine amrdm_make_knapsack

   !> Build an SFC-optimized DM from per-box costs + BoxArray
   subroutine amrdm_make_sfc(dm,costs,nboxes,ba)
      use amrex_amr_module, only: amrex_distromap,amrex_boxarray
      type(amrex_distromap), intent(out) :: dm
      integer, intent(in) :: nboxes
      real(WP), intent(in) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      call amrdm_make_sfc_c(dm%p,costs,nboxes,ba%p)
      dm%owner=.true.
   end subroutine amrdm_make_sfc

end module amrex_interface

