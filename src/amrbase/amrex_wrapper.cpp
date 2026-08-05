// AMReX Interface for NGA2
// C++ bridge between NGA2 Fortran and AMReX C++
// Provides amrcore_* and amrmfab_* functions callable from Fortran

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

namespace nga2 {

//=============================================================================
// Callback function types - these are bind(c) Fortran procedures
// They receive the owner pointer and handle dispatching in Fortran
//
// Notes on AMReX callback parameters:
// - ngrow in ErrorEst: AMReX always passes 0, not needed
// - initial in regrid: Unused in AmrCore implementation, dropped
//=============================================================================
extern "C" {
// Level callbacks: owner, level, time, boxarray ptr, distromap ptr
using LevelDispatcher = void (*)(void *owner, int lev, double time,
                                 const amrex::BoxArray *,
                                 const amrex::DistributionMapping *);
// Clear callback: owner, level
using ClearDispatcher = void (*)(void *owner, int lev);
// Tag callback: owner, level, tagboxarray ptr, time
// Note: SET/CLEAR values are constants (char(2)/char(0)) - defined in Fortran
using TagDispatcher = void (*)(void *owner, int lev, amrex::TagBoxArray *,
                               double time);
// Post-regrid callback: owner, lbase, time
// Called after regrid completes; lbase is the base level for regridding
using PostRegridDispatcher = void (*)(void *owner, int lbase, double time);
// Cost callback: owner fills per-box cost array for load balancing
// costs array is pre-filled with 1.0 (uniform); callback overwrites with actual costs
// ba_ptr points to the new BoxArray being distributed
using CostDispatcher = void (*)(void *owner, int lev, int nboxes, double *costs,
                                const void *ba_ptr);
}

//=============================================================================
// NGA2AmrCore - extends AmrCore directly (bypasses FAmrCore which
// enforces isotropic refinement ratios)
//=============================================================================
class NGA2AmrCore : public amrex::AmrCore {
public:
  void *owner = nullptr;

  LevelDispatcher on_init_dispatch = nullptr;
  LevelDispatcher on_coarse_dispatch = nullptr;
  LevelDispatcher on_remake_dispatch = nullptr;
  ClearDispatcher on_clear_dispatch = nullptr;
  TagDispatcher on_tag_dispatch = nullptr;
  PostRegridDispatcher on_postregrid_dispatch = nullptr;
  CostDispatcher on_cost_dispatch = nullptr;
  int cost_strategy = 0;  // 0=SFC (AMReX default), 1=KnapSack

  // Public override of regrid to call post-regrid dispatch after base class
  // Note: 'initial' parameter from AMReX is unused, we drop it
  void regrid(int lbase, amrex::Real time, bool initial = false) override {
    amrex::AmrCore::regrid(lbase, time, initial);
    if (on_postregrid_dispatch && owner) {
      on_postregrid_dispatch(owner, lbase, time);
    }
  }

  // Override MakeDistributionMap for cost-aware load balancing
  // If a cost callback is registered, calls it to get per-box costs,
  // then distributes using the selected strategy.
  // If no callback, falls back to default (SFC by cell count).
  //   cost_strategy: 0 = SFC (default), 1 = KnapSack
  amrex::DistributionMapping
  MakeDistributionMap(int lev, amrex::BoxArray const &ba) override {
    if (on_cost_dispatch && owner) {
      int nboxes = static_cast<int>(ba.size());
      amrex::Vector<amrex::Real> costs(nboxes, 1.0);
      on_cost_dispatch(owner, lev, nboxes, costs.data(), &ba);
      switch (cost_strategy) {
      case 0:
        return amrex::DistributionMapping::makeSFC(costs, ba);
      case 1:
        return amrex::DistributionMapping::makeKnapSack(costs);
      default:
        amrex::Abort("NGA2AmrCore::MakeDistributionMap: unknown cost_strategy="
                     + std::to_string(cost_strategy)
                     + " (valid: 0=SFC, 1=KnapSack)");
      }
    }
    return amrex::AmrCore::MakeDistributionMap(lev, ba);
  }

  // Public wrapper for MakeNewLevelFromScratch (protected in AmrCore)
  // Needed by amrcore_read_grids to build levels from checkpoint
  void BuildLevelFromScratch(int lev, amrex::Real time,
                             const amrex::BoxArray &ba,
                             const amrex::DistributionMapping &dm) {
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    MakeNewLevelFromScratch(lev, time, ba, dm);
  }
protected:
  void MakeNewLevelFromScratch(int lev, amrex::Real time,
                               const amrex::BoxArray &ba,
                               const amrex::DistributionMapping &dm) override {
    if (on_init_dispatch && owner) {
      on_init_dispatch(owner, lev, time, &ba, &dm);
    }
  }

  void MakeNewLevelFromCoarse(int lev, amrex::Real time,
                              const amrex::BoxArray &ba,
                              const amrex::DistributionMapping &dm) override {
    if (on_coarse_dispatch && owner) {
      on_coarse_dispatch(owner, lev, time, &ba, &dm);
    }
  }

  void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray &ba,
                   const amrex::DistributionMapping &dm) override {
    if (on_remake_dispatch && owner) {
      on_remake_dispatch(owner, lev, time, &ba, &dm);
    }
  }

  void ClearLevel(int lev) override {
    if (on_clear_dispatch && owner) {
      on_clear_dispatch(owner, lev);
    }
  }

  // Note: ngrow parameter from AMReX is always 0, we drop it
  // SET/CLEAR are TagBox::SET=char(2) and TagBox::CLEAR=char(0)
  void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time,
                int /*ngrow*/) override {
    if (on_tag_dispatch && owner) {
      on_tag_dispatch(owner, lev, &tags, time);
    }
  }
};

//=============================================================================
// FillPatch BC Callback - receives solver context, can access all fields
//=============================================================================

// Fortran BC callback: ctx first, then AMReX order (mf, scomp, ncomp, time,
// geom)
using FillPatchBCDispatcher = void (*)(void *solver_ctx, amrex::MultiFab *mf,
                                       int scomp, int ncomp, amrex::Real time,
                                       const amrex::Geometry *geom);

// Functor that AMReX FillPatch calls; forwards to Fortran with context
class NGA2BCFunctor {
public:
  void *solver_ctx = nullptr;
  FillPatchBCDispatcher bc_dispatch = nullptr;
  const amrex::Geometry *geom = nullptr;

  NGA2BCFunctor() = default;
  NGA2BCFunctor(void *ctx, FillPatchBCDispatcher f, const amrex::Geometry *g)
      : solver_ctx(ctx), bc_dispatch(f), geom(g) {}

  // Called by AMReX FillPatch when physical BCs needed
  void operator()(amrex::MultiFab &mf, int icomp, int ncomp,
                  amrex::IntVect const & /*nghost*/, amrex::Real time,
                  int /*bccomp*/) {
    if (bc_dispatch && solver_ctx) {
      bc_dispatch(solver_ctx, &mf, icomp, ncomp, time, geom);
    }
  }
};

} // namespace nga2

//=============================================================================
// C Interface for Fortran - all functions use amrcore_* prefix
//=============================================================================
extern "C" {

//-----------------------------------------------------------------------------
// Lifecycle
//-----------------------------------------------------------------------------
void *amrcore_create() { return new nga2::NGA2AmrCore(); }

void amrcore_destroy(void *core) {
  delete static_cast<nga2::NGA2AmrCore *>(core);
}

//-----------------------------------------------------------------------------
// Set Owner
//-----------------------------------------------------------------------------
void amrcore_set_owner(void *core, void *owner) {
  static_cast<nga2::NGA2AmrCore *>(core)->owner = owner;
}

//-----------------------------------------------------------------------------
// Set Dispatchers
//-----------------------------------------------------------------------------
void amrcore_set_on_init_dispatch(void *core, nga2::LevelDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_init_dispatch = f;
}

void amrcore_set_on_coarse_dispatch(void *core, nga2::LevelDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_coarse_dispatch = f;
}

void amrcore_set_on_remake_dispatch(void *core, nga2::LevelDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_remake_dispatch = f;
}

void amrcore_set_on_clear_dispatch(void *core, nga2::ClearDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_clear_dispatch = f;
}

void amrcore_set_on_tag_dispatch(void *core, nga2::TagDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_tag_dispatch = f;
}

void amrcore_set_on_postregrid_dispatch(void *core,
                                        nga2::PostRegridDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_postregrid_dispatch = f;
}

void amrcore_set_on_cost_dispatch(void *core, nga2::CostDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_cost_dispatch = f;
}

void amrcore_set_cost_strategy(void *core, int strategy) {
  static_cast<nga2::NGA2AmrCore *>(core)->cost_strategy = strategy;
}

//-----------------------------------------------------------------------------
// Grid Operations
//-----------------------------------------------------------------------------
void amrcore_init_from_scratch(void *core, double time) {
  static_cast<nga2::NGA2AmrCore *>(core)->InitFromScratch(time);
}

void amrcore_regrid(void *core, int base_level, double time) {
  static_cast<nga2::NGA2AmrCore *>(core)->regrid(base_level, time);
}

//-----------------------------------------------------------------------------
// Level Queries
//-----------------------------------------------------------------------------
int amrcore_finest_level(void *core) {
  return static_cast<nga2::NGA2AmrCore *>(core)->finestLevel();
}

int amrcore_max_level(void *core) {
  return static_cast<nga2::NGA2AmrCore *>(core)->maxLevel();
}

void amrcore_set_finest_level(void *core, int lev) {
  static_cast<nga2::NGA2AmrCore *>(core)->SetFinestLevel(lev);
}

//-----------------------------------------------------------------------------
// Geometry Access
//-----------------------------------------------------------------------------
void amrcore_get_geometry(void **geom_ptr, int lev, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  *geom_ptr = const_cast<amrex::Geometry *>(&(amr->Geom(lev)));
}

void amrcore_get_ref_ratio(int *ref_ratio, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  for (int lev = 0; lev < amr->maxLevel(); ++lev) {
    ref_ratio[lev] = amr->refRatio(lev)[0];
  }
}

void amrcore_get_ref_ratio_xyz(int *rrefx, int *rrefy, int *rrefz, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  for (int lev = 0; lev < amr->maxLevel(); ++lev) {
    rrefx[lev] = amr->refRatio(lev)[0];
    rrefy[lev] = amr->refRatio(lev)[1];
    rrefz[lev] = amr->refRatio(lev)[2];
  }
}

//-----------------------------------------------------------------------------
// BoxArray and DistributionMapping Access
//-----------------------------------------------------------------------------
void amrcore_get_boxarray(void **ba_ptr, int lev, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  *ba_ptr = const_cast<amrex::BoxArray *>(&(amr->boxArray(lev)));
}

void amrcore_get_distromap(void **dm_ptr, int lev, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  *dm_ptr =
      const_cast<amrex::DistributionMapping *>(&(amr->DistributionMap(lev)));
}

//-----------------------------------------------------------------------------
// Standalone cost-weighted DistributionMapping factories
//   - amrdm_make_knapsack: KnapSack algorithm (cost vector only)
//   - amrdm_make_sfc:      Space-filling curve (cost vector + BoxArray)
//   - amrdm_destroy:       Free a heap-allocated DM
//-----------------------------------------------------------------------------
void amrdm_make_knapsack(void **dm_out, double *costs, int nboxes) {
  amrex::Vector<amrex::Real> rcost(costs, costs + nboxes);
  auto *dm = new amrex::DistributionMapping(
      amrex::DistributionMapping::makeKnapSack(rcost));
  *dm_out = dm;
}

void amrdm_make_sfc(void **dm_out, double *costs, int nboxes, void *ba_ptr) {
  amrex::Vector<amrex::Real> rcost(costs, costs + nboxes);
  auto *ba = static_cast<amrex::BoxArray *>(ba_ptr);
  auto *dm = new amrex::DistributionMapping(
      amrex::DistributionMapping::makeSFC(rcost, *ba));
  *dm_out = dm;
}

void amrdm_destroy(void *dm) {
  delete static_cast<amrex::DistributionMapping *>(dm);
}

//=============================================================================
// MultiFab Operations - amrmfab_* prefix
//=============================================================================

//-----------------------------------------------------------------------------
// FillPatch for Level 0 (single level, just physical BCs)
// Note: scomp/dcomp are 1-indexed (Fortran convention), converted to 0-indexed
// here
//-----------------------------------------------------------------------------
void amrmfab_fillpatch_single(void *mf_ptr, double time_old, void *mf_old_ptr,
                              double time_new, void *mf_new_ptr, void *geom_ptr,
                              void *solver_ctx,
                              nga2::FillPatchBCDispatcher bc_dispatch,
                              double time, int scomp, int dcomp, int ncomp,
                              int nghost) {
  auto *mf = static_cast<amrex::MultiFab *>(mf_ptr);
  auto *mf_old = static_cast<amrex::MultiFab *>(mf_old_ptr);
  auto *mf_new = static_cast<amrex::MultiFab *>(mf_new_ptr);
  auto *geom = static_cast<amrex::Geometry *>(geom_ptr);

  amrex::Vector<amrex::MultiFab *> smf = {mf_old, mf_new};
  amrex::Vector<amrex::Real> stime = {time_old, time_new};

  nga2::NGA2BCFunctor bc_functor(solver_ctx, bc_dispatch, geom);

  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  // Use explicit nghost if non-negative, otherwise use mf.nGrowVect()
  amrex::IntVect ng = (nghost >= 0) ? amrex::IntVect(nghost) : mf->nGrowVect();
  amrex::FillPatchSingleLevel(*mf, ng, time, smf, stime, scomp - 1, dcomp - 1,
                              ncomp, *geom, bc_functor, 0);
}
//-----------------------------------------------------------------------------
// FillPatch for Fine Levels (two levels: coarse + fine)
// Note: scomp/dcomp are 1-indexed (Fortran convention), converted to 0-indexed
// here
//-----------------------------------------------------------------------------
void amrmfab_fillpatch_two(void *mf_ptr, double time_old_c, void *mf_old_c_ptr,
                           double time_new_c, void *mf_new_c_ptr,
                           void *geom_c_ptr, double time_old_f,
                           void *mf_old_f_ptr, double time_new_f,
                           void *mf_new_f_ptr, void *geom_f_ptr,
                           void *solver_ctx,
                           nga2::FillPatchBCDispatcher bc_dispatch, double time,
                           int scomp, int dcomp, int ncomp, int *ref_ratio,
                           int interp_type, int *lo_bc, int *hi_bc, int nbc,
                           int nghost) {
  auto *mf = static_cast<amrex::MultiFab *>(mf_ptr);
  auto *mf_old_c = static_cast<amrex::MultiFab *>(mf_old_c_ptr);
  auto *mf_new_c = static_cast<amrex::MultiFab *>(mf_new_c_ptr);
  auto *mf_old_f = static_cast<amrex::MultiFab *>(mf_old_f_ptr);
  auto *mf_new_f = static_cast<amrex::MultiFab *>(mf_new_f_ptr);
  auto *geom_c = static_cast<amrex::Geometry *>(geom_c_ptr);
  auto *geom_f = static_cast<amrex::Geometry *>(geom_f_ptr);

  // Look up interpolator from type constant (matches AMReX Fortran constants)
  // See amrex_interpolater_module in AMReX_interpolater_mod.F90
  amrex::Interpolater *interp = nullptr;
  switch (interp_type) {
  case 0: // amrex_interp_pc
    interp = &amrex::pc_interp;
    break;
  case 1: // amrex_interp_node_bilinear
    interp = &amrex::node_bilinear_interp;
    break;
  case 2: // amrex_interp_cell_bilinear
    interp = &amrex::cell_bilinear_interp;
    break;
  case 3: // amrex_interp_quadratic
    interp = &amrex::quadratic_interp;
    break;
  case 4: // amrex_interp_lincc
    interp = &amrex::lincc_interp;
    break;
  case 6: // amrex_interp_protected
    interp = &amrex::protected_interp;
    break;
  case 7: // amrex_interp_quartic
    interp = &amrex::quartic_interp;
    break;
  case 8: // amrex_interp_face_divfree
    interp = &amrex::face_divfree_interp;
    break;
  case 9: // amrex_interp_face_linear
    interp = &amrex::face_linear_interp;
    break;
  case 5: // amrex_interp_cell_cons
  default:
    interp = &amrex::cell_cons_interp;
    break;
  }

  amrex::Vector<amrex::MultiFab *> cmf = {mf_old_c, mf_new_c};
  amrex::Vector<amrex::Real> ctime = {time_old_c, time_new_c};
  amrex::Vector<amrex::MultiFab *> fmf = {mf_old_f, mf_new_f};
  amrex::Vector<amrex::Real> ftime = {time_old_f, time_new_f};

  // Build BCRec from lo_bc/hi_bc arrays
  amrex::Vector<amrex::BCRec> bcs(ncomp);
  for (int i = 0; i < ncomp; ++i) {
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::low),
               lo_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::high),
               hi_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::low),
               lo_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::high),
               hi_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::low),
               lo_bc[i * 3 + 2]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::high),
               hi_bc[i * 3 + 2]);
  }

  nga2::NGA2BCFunctor bc_functor_c(solver_ctx, bc_dispatch, geom_c);
  nga2::NGA2BCFunctor bc_functor_f(solver_ctx, bc_dispatch, geom_f);

  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));

  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  // Use explicit nghost if non-negative, otherwise use mf.nGrowVect()
  amrex::IntVect ng = (nghost >= 0) ? amrex::IntVect(nghost) : mf->nGrowVect();
  amrex::FillPatchTwoLevels(*mf, ng, time, cmf, ctime, fmf, ftime, scomp - 1,
                            dcomp - 1, ncomp, *geom_c, *geom_f, bc_functor_c, 0,
                            bc_functor_f, 0, ratio, interp, bcs, 0);
}

//-----------------------------------------------------------------------------
// FillCoarsePatch - fill fine level from coarse only (no fine data needed)
// Note: scomp/dcomp are 1-indexed (Fortran convention), converted to 0-indexed
// here
//-----------------------------------------------------------------------------
void amrmfab_fillcoarsepatch(void *mf_f_ptr, double time, void *mf_c_ptr,
                             void *geom_c_ptr, void *geom_f_ptr,
                             void *solver_ctx,
                             nga2::FillPatchBCDispatcher bc_dispatch, int scomp,
                             int dcomp, int ncomp, int *ref_ratio,
                             int interp_type, int *lo_bc, int *hi_bc, int nbc) {
  auto *mf_f = static_cast<amrex::MultiFab *>(mf_f_ptr);
  auto *mf_c = static_cast<amrex::MultiFab *>(mf_c_ptr);
  auto *geom_c = static_cast<amrex::Geometry *>(geom_c_ptr);
  auto *geom_f = static_cast<amrex::Geometry *>(geom_f_ptr);

  // Look up interpolator (see amrex_interpolater_module)
  amrex::Interpolater *interp = nullptr;
  switch (interp_type) {
  case 0: // amrex_interp_pc
    interp = &amrex::pc_interp;
    break;
  case 1: // amrex_interp_node_bilinear
    interp = &amrex::node_bilinear_interp;
    break;
  case 2: // amrex_interp_cell_bilinear
    interp = &amrex::cell_bilinear_interp;
    break;
  case 3: // amrex_interp_quadratic
    interp = &amrex::quadratic_interp;
    break;
  case 4: // amrex_interp_lincc
    interp = &amrex::lincc_interp;
    break;
  case 6: // amrex_interp_protected
    interp = &amrex::protected_interp;
    break;
  case 7: // amrex_interp_quartic
    interp = &amrex::quartic_interp;
    break;
  case 8: // amrex_interp_face_divfree
    interp = &amrex::face_divfree_interp;
    break;
  case 9: // amrex_interp_face_linear
    interp = &amrex::face_linear_interp;
    break;
  case 5: // amrex_interp_cell_cons
  default:
    interp = &amrex::cell_cons_interp;
    break;
  }

  // Build BCRec
  amrex::Vector<amrex::BCRec> bcs(ncomp);
  for (int i = 0; i < ncomp; ++i) {
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::low),
               lo_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::high),
               hi_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::low),
               lo_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::high),
               hi_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::low),
               lo_bc[i * 3 + 2]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::high),
               hi_bc[i * 3 + 2]);
  }

  nga2::NGA2BCFunctor bc_functor_c(solver_ctx, bc_dispatch, geom_c);
  nga2::NGA2BCFunctor bc_functor_f(solver_ctx, bc_dispatch, geom_f);

  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));

  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  // Use InterpFromCoarseLevel which doesn't require fine-level source data
  amrex::InterpFromCoarseLevel(*mf_f, time, *mf_c, scomp - 1, dcomp - 1, ncomp,
                               *geom_c, *geom_f, bc_functor_c, 0, bc_functor_f,
                               0, ratio, interp, bcs, 0);
}

//-----------------------------------------------------------------------------
// MultiFab sum_unique - avoids double counting at shared nodes/faces
//-----------------------------------------------------------------------------

double amrmfab_sum_unique(void *mf_ptr, void *geom_ptr, int comp) {
  auto *mf = static_cast<amrex::MultiFab *>(mf_ptr);
  auto *geom = static_cast<amrex::Geometry *>(geom_ptr);
  // Build periodicity from geometry
  amrex::Periodicity period = geom->periodicity();
  // Use sum_unique with 0-indexed component
  return mf->sum_unique(comp - 1, false, period);
}

//-----------------------------------------------------------------------------
// Make fine mask - identifies cells covered by finer level
//-----------------------------------------------------------------------------

void amrmask_make_fine(void *mask_ptr, void *ba_fine_ptr, int *ref_ratio,
                       int covered_val, int notcovered_val) {
  auto *mask = static_cast<amrex::iMultiFab *>(mask_ptr);
  auto *ba_fine = static_cast<amrex::BoxArray *>(ba_fine_ptr);
  amrex::IntVect rr(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  // Use makeFineMask to fill the mask
  // covered_val = value where fine grids exist
  // notcovered_val = value where no fine grids
  auto fine_mask = amrex::makeFineMask(*mask, *ba_fine, rr, notcovered_val, covered_val);
  // Copy result into provided mask
  amrex::iMultiFab::Copy(*mask, fine_mask, 0, 0, 1, 0);
}

//-----------------------------------------------------------------------------
// Plotfile Output (Native and HDF5)
//-----------------------------------------------------------------------------


void amrplotfile_write_native(const char *name, int nlevels, void **mf_ptrs,
                              const char **varnames, int ncomp,
                              void **geom_ptrs, double time, int *level_steps,
                              int *ref_ratios) {
  // Convert void* arrays to AMReX types
  amrex::Vector<const amrex::MultiFab *> mf_vec(nlevels);
  amrex::Vector<amrex::Geometry> geom_vec(nlevels);
  amrex::Vector<int> steps_vec(nlevels);
  amrex::Vector<amrex::IntVect> rr_vec(nlevels - 1);
  amrex::Vector<std::string> varname_vec(ncomp);

  for (int lev = 0; lev < nlevels; ++lev) {
    mf_vec[lev] = static_cast<amrex::MultiFab *>(mf_ptrs[lev]);
    geom_vec[lev] = *static_cast<amrex::Geometry *>(geom_ptrs[lev]);
    steps_vec[lev] = level_steps[lev];
  }
  for (int lev = 0; lev < nlevels - 1; ++lev) {
    rr_vec[lev] = amrex::IntVect(
        AMREX_D_DECL(ref_ratios[3*lev+0], ref_ratios[3*lev+1], ref_ratios[3*lev+2]));
  }
  for (int i = 0; i < ncomp; ++i) {
    varname_vec[i] = std::string(varnames[i]);
  }

  // Pre-delete existing directory to avoid AMReX ".old" rename
  std::string dirstr(name);
  if (amrex::ParallelDescriptor::IOProcessor()) {
    if (amrex::FileExists(dirstr)) {
      amrex::FileSystem::RemoveAll(dirstr);
    }
  }
  amrex::ParallelDescriptor::Barrier();

  amrex::WriteMultiLevelPlotfile(dirstr, nlevels, mf_vec,
                                 varname_vec, geom_vec, time, steps_vec,
                                 rr_vec);
}

#ifdef AMREX_USE_HDF5
void amrplotfile_write_hdf5(const char *name, int nlevels, void **mf_ptrs,
                            const char **varnames, int ncomp, void **geom_ptrs,
                            double time, int *level_steps, int *ref_ratios,
                            const char *compression) {
  // Convert void* arrays to AMReX types
  amrex::Vector<const amrex::MultiFab *> mf_vec(nlevels);
  amrex::Vector<amrex::Geometry> geom_vec(nlevels);
  amrex::Vector<int> steps_vec(nlevels);
  amrex::Vector<amrex::IntVect> rr_vec(nlevels - 1);
  amrex::Vector<std::string> varname_vec(ncomp);

  for (int lev = 0; lev < nlevels; ++lev) {
    mf_vec[lev] = static_cast<amrex::MultiFab *>(mf_ptrs[lev]);
    geom_vec[lev] = *static_cast<amrex::Geometry *>(geom_ptrs[lev]);
    steps_vec[lev] = level_steps[lev];
  }
  for (int lev = 0; lev < nlevels - 1; ++lev) {
    rr_vec[lev] = amrex::IntVect(
        AMREX_D_DECL(ref_ratios[3*lev+0], ref_ratios[3*lev+1], ref_ratios[3*lev+2]));
  }
  for (int i = 0; i < ncomp; ++i) {
    varname_vec[i] = std::string(varnames[i]);
  }

  std::string comp_str = compression ? std::string(compression) : "";
  amrex::WriteMultiLevelPlotfileHDF5(std::string(name), nlevels, mf_vec,
                                     varname_vec, geom_vec, time, steps_vec,
                                     rr_vec, comp_str);
}
#endif

//-----------------------------------------------------------------------------
// Checkpoint I/O (VisMF)
//-----------------------------------------------------------------------------
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

// Write a MultiFab to disk using VisMF
void amrmfab_vismf_write(void *mf, const char *path) {
  amrex::VisMF::Write(*static_cast<amrex::MultiFab *>(mf), std::string(path));
}

// Read a MultiFab from disk using VisMF
// NOTE: The MultiFab must already be defined with correct BoxArray/DistroMap
void amrmfab_vismf_read(void *mf, const char *path) {
  amrex::VisMF::Read(*static_cast<amrex::MultiFab *>(mf), std::string(path));
}

// Create directory hierarchy for checkpoints
void amrcheckpoint_prebuild_dirs(const char *dirname, const char *subdir_prefix,
                                 int nlevels) {
  // Only I/O processor removes existing directory to avoid MPI race condition
  std::string dirstr(dirname);
  if (amrex::ParallelDescriptor::IOProcessor()) {
    if (amrex::FileExists(dirstr)) {
      amrex::FileSystem::RemoveAll(dirstr);
    }
  }
  amrex::ParallelDescriptor::Barrier(); // Wait for removal to complete
  amrex::PreBuildDirectorHierarchy(dirstr, std::string(subdir_prefix), nlevels,
                                   true);
}

// Get MultiFab file prefix for a level (e.g., "chk00100/Level_0/phi")
// Returns the full path that can be used with VisMF::Write/Read
void amrcheckpoint_mfab_prefix(char *result, int result_len, int lev,
                               const char *dirname, const char *level_prefix,
                               const char *mfab_name) {
  std::string prefix = amrex::MultiFabFileFullPrefix(lev, std::string(dirname),
                                                     std::string(level_prefix),
                                                     std::string(mfab_name));
  // Copy to result buffer
  size_t len = std::min(prefix.size(), static_cast<size_t>(result_len - 1));
  std::copy(prefix.begin(), prefix.begin() + len, result);
  result[len] = '\0';
}

// Set number of output files for VisMF (controls I/O aggregation)
// nfiles=1 gives single file output; nfiles=nranks gives file-per-rank
void amrvismf_set_noutfiles(int nfiles) { amrex::VisMF::SetNOutFiles(nfiles); }

// Get current number of output files setting
int amrvismf_get_noutfiles() { return amrex::VisMF::GetNOutFiles(); }

// Build a single level from pre-built BoxArray and DistributionMapping
// Called from Fortran after reading box data and constructing BA/DM
// Fires on_init callbacks (MakeNewLevelFromScratch)
void amrcore_build_level(void *core, int lev, double time,
                         void *ba_ptr, void *dm_ptr) {
  auto *amrcore = static_cast<nga2::NGA2AmrCore *>(core);
  auto &ba = *static_cast<amrex::BoxArray *>(ba_ptr);
  auto &dm = *static_cast<amrex::DistributionMapping *>(dm_ptr);
  amrcore->BuildLevelFromScratch(lev, time, ba, dm);
}
//=============================================================================
// HDF5 Plotfile Utilities
//=============================================================================
// Read time from a native AMReX plotfile directory (reads Header text file)
// Header format: version, ncomp, var names (ncomp lines), spacedim, time, ...
// Returns -1.0 if not found or unreadable
static double read_time_from_native(const char *dirname) {
  std::string header_path = std::string(dirname) + "/Header";
  std::ifstream ifs(header_path);
  if (!ifs.good()) return -1.0;

  // Line 1: version string
  std::string line;
  std::getline(ifs, line);

  // Line 2: number of components
  int ncomp = 0;
  ifs >> ncomp;
  std::getline(ifs, line); // consume newline

  // Skip ncomp variable name lines
  for (int i = 0; i < ncomp; ++i) {
    std::getline(ifs, line);
  }

  // Next line: spacedim
  int sdim = 0;
  ifs >> sdim;

  // Next line: time
  double time = -1.0;
  ifs >> time;

  return ifs.good() ? time : -1.0;
}

// Read the time from an AMReX plotfile (native directory or HDF5 file)
// Tries native Header first, then HDF5 if available
// Returns -1.0 if unreadable
double amrplotfile_read_time(const char *filename) {
  // Try native format first (plotfile is a directory with Header)
  double time = read_time_from_native(filename);
  if (time >= 0.0) return time;

#ifdef AMREX_USE_HDF5
#include <hdf5.h>
  // Try HDF5 format
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) return -1.0;

  if (H5Aexists(file_id, "time") > 0) {
    hid_t attr_id = H5Aopen(file_id, "time", H5P_DEFAULT);
    if (attr_id >= 0) {
      double time_arr[1];
      H5Aread(attr_id, H5T_NATIVE_DOUBLE, time_arr);
      time = time_arr[0];
      H5Aclose(attr_id);
    }
  }
  H5Fclose(file_id);
  return time;
#else
  return -1.0;
#endif
}

//=============================================================================
// MLMG Utilities (not available in AMReX Fortran interface)
//=============================================================================

// Get number of iterations from last solve
int amrmlmg_get_niters(void *mlmg) {
  return static_cast<amrex::MLMG *>(mlmg)->getNumIters();
}

// Get face-centered fluxes from MLMG solution using solver's C/F stencils
// For (alpha*A - beta*div(B*grad))phi = rhs, flux = -B*grad(phi)
//   mlmg - pointer to MLMG object (must have valid solution from solve)
//   sol_mfs - array of solution MultiFab pointers (one per level)
//   flux_x, flux_y, flux_z - arrays of flux MultiFab pointers (one per level)
//   nlevs - number of AMR levels
void amrmlmg_get_fluxes(void *mlmg, void **sol_mfs, void **flux_x,
                        void **flux_y, void **flux_z, int nlevs) {
  auto *mg = static_cast<amrex::MLMG *>(mlmg);

  // Build solution vector
  amrex::Vector<amrex::MultiFab *> sol(nlevs);
  for (int lev = 0; lev < nlevs; ++lev) {
    sol[lev] = static_cast<amrex::MultiFab *>(sol_mfs[lev]);
  }

  // Build flux vector of arrays (one Array<MultiFab*,3> per level)
  amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fluxes(nlevs);
  for (int lev = 0; lev < nlevs; ++lev) {
    fluxes[lev][0] = static_cast<amrex::MultiFab *>(flux_x[lev]);
    fluxes[lev][1] = static_cast<amrex::MultiFab *>(flux_y[lev]);
    fluxes[lev][2] = static_cast<amrex::MultiFab *>(flux_z[lev]);
  }

  // Call MLMG getFluxes with solution - uses solver's C/F stencils
  mg->getFluxes(fluxes, sol, amrex::MLMG::Location::FaceCenter);
}

// Composite dot product for PCG with AMR: sums all valid cells at all levels.
// Covered coarse cells are included because comp_residual's average_down
// overwrites them with averaged fine values — they participate in the actual
// residual and must be accounted for in convergence checking.
//   mf1_ptrs, mf2_ptrs  - arrays of MultiFab pointers (one per level, 0-indexed)
//   ba_fine_ptrs, ref_ratios - reserved (unused, kept for ABI compatibility)
//   nlevs                - number of AMR levels
double amrmlmg_dot_composite(void **mf1_ptrs, void **mf2_ptrs,
                              void **ba_fine_ptrs, int *ref_ratios, int nlevs) {
  amrex::Real result = 0.0;
  for (int lev = 0; lev < nlevs; ++lev) {
    auto *mf1 = static_cast<amrex::MultiFab *>(mf1_ptrs[lev]);
    auto *mf2 = static_cast<amrex::MultiFab *>(mf2_ptrs[lev]);
    result += amrex::MultiFab::Dot(*mf1, 0, *mf2, 0, 1, 0);
  }
  return static_cast<double>(result);
}

//=============================================================================
// MultiFab Averaging Utilities (Unified API)
// All 4 types (cell, face, edge, node) have the same signature pattern:
//   (fine_mf, crse_mf, crse_geom, ref_ratio, ngcrse)
// When crse_geom is provided (non-null), FillBoundary is called for periodic
// ghost fix-up. ngcrse controls how many ghost cells to average into.
//=============================================================================

// Average down cell-centered MultiFab
// If crse_geom is provided, uses ParallelCopy with periodicity to propagate
// averaged valid cells to periodic partners, then FillBoundary for ghost cells.
void amrmfab_average_down_cell(void *fine_mf, void *crse_mf, void *crse_geom,
                               int *ref_ratio, int ngcrse) {
  auto *fmf = static_cast<amrex::MultiFab *>(fine_mf);
  auto *cmf = static_cast<amrex::MultiFab *>(crse_mf);
  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  int ncomp = cmf->nComp();

  // Always average into temp first, then copy to crse with periodicity
  amrex::BoxArray crse_fine_BA = fmf->boxArray();
  crse_fine_BA.coarsen(ratio);
  amrex::MultiFab ctmp(crse_fine_BA, fmf->DistributionMap(), ncomp, ngcrse);

  for (amrex::MFIter mfi(ctmp, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    amrex::Box const &bx = mfi.growntilebox(ngcrse);
    auto const &crsearr = ctmp.array(mfi);
    auto const &finearr = fmf->const_array(mfi);
    amrex::ParallelFor(
        bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          amrex::amrex_avgdown(i, j, k, n, crsearr, finearr, 0, 0, ratio);
        });
  }

  // Copy averaged data to crse - with periodicity if geometry provided
  if (crse_geom) {
    auto *cg = static_cast<amrex::Geometry *>(crse_geom);
    // ParallelCopy with periodicity propagates valid cells to periodic partners
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse, cg->periodicity());
    // FillBoundary syncs ghost cells from valid cells
    cmf->FillBoundary(cg->periodicity());
  } else {
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse);
  }
}

// Restrict-SUM fine deposits (valid+ghost) into the coarse level with ADD semantics.
// Delegates to AMReX's sum_fine_to_coarse (AMReX_MultiFabUtil.H), which:
//   - requires nGrow % ratio == 0 (i.e., nover must be a multiple of refinement ratio)
//   - iterates over growntilebox(nGrow/ratio) to capture ghost deposits at C/F boundaries
//   - uses ParallelCopy(..., nGrow, IntVect(0), ..., ADD) to merge into coarse valid cells
// Call after SumBoundary at the fine level; average_down afterward fixes double-counted cells.
extern "C" void amrmfab_sum_downto(void *fine_mf_ptr, void *crse_mf_ptr,
                                   void *crse_geom_ptr, void *fine_geom_ptr,
                                   const int *ref_ratio) {
    auto *fmf   = static_cast<amrex::MultiFab *>(fine_mf_ptr);
    auto *cmf   = static_cast<amrex::MultiFab *>(crse_mf_ptr);
    auto *cgeom = static_cast<amrex::Geometry *>(crse_geom_ptr);
    auto *fgeom = static_cast<amrex::Geometry *>(fine_geom_ptr);
    amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
    amrex::sum_fine_to_coarse(*fmf, *cmf, 0, fmf->nComp(), ratio, *cgeom, *fgeom);
}

// Interpolate coarse-level cell-centered data onto a fine-level MultiFab using
// piecewise-constant (PCInterp) interpolation with no-op physical BCs.
// Used to propagate coarse particle deposits into fine-covered cells so that
// average_down preserves them.  Mirrors AMReX AssignDensity's InterpFromCoarseLevel call.
// scomp is 0-indexed (C convention).
extern "C" void amrmfab_interp_from_coarse(void *fine_mf_ptr, void *crse_mf_ptr,
                                           void *crse_geom_ptr, void *fine_geom_ptr,
                                           int scomp, int ncomp, const int *ref_ratio) {
    auto *fmf   = static_cast<amrex::MultiFab *>(fine_mf_ptr);
    auto *cmf   = static_cast<amrex::MultiFab *>(crse_mf_ptr);
    auto *cgeom = static_cast<amrex::Geometry *>(crse_geom_ptr);
    auto *fgeom = static_cast<amrex::Geometry *>(fine_geom_ptr);
    amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));

    int lo_bc[] = {amrex::BCType::int_dir, amrex::BCType::int_dir, amrex::BCType::int_dir};
    int hi_bc[] = {amrex::BCType::int_dir, amrex::BCType::int_dir, amrex::BCType::int_dir};
    amrex::Vector<amrex::BCRec> bcs(ncomp, amrex::BCRec(lo_bc, hi_bc));
    amrex::PCInterp mapper;
    amrex::PhysBCFunctNoOp cbc, fbc;

    amrex::InterpFromCoarseLevel(*fmf, 0.0, *cmf, scomp, scomp, ncomp,
                                 *cgeom, *fgeom, cbc, 0, fbc, 0,
                                 ratio, &mapper, bcs, 0);
}

// Average down face-centered MultiFab (nodal in 1 dir, cell in 2)
// If crse_geom is provided, uses ParallelCopy with periodicity to propagate
// averaged valid cells to periodic partners, then FillBoundary for ghost cells.
void amrmfab_average_down_face(void *fine_mf, void *crse_mf, void *crse_geom,
                               int *ref_ratio, int ngcrse) {
  auto *fmf = static_cast<amrex::MultiFab *>(fine_mf);
  auto *cmf = static_cast<amrex::MultiFab *>(crse_mf);
  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  int ncomp = cmf->nComp();

  // Average into temp, then copy to crse with periodicity
  amrex::BoxArray crse_fine_BA = amrex::coarsen(fmf->boxArray(), ratio);
  amrex::MultiFab ctmp(crse_fine_BA, fmf->DistributionMap(), ncomp, ngcrse);
  amrex::average_down_faces(*fmf, ctmp, ratio, ngcrse);

  // Copy averaged data to crse - with periodicity if geometry provided
  if (crse_geom) {
    auto *cg = static_cast<amrex::Geometry *>(crse_geom);
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse, cg->periodicity());
    cmf->FillBoundary(cg->periodicity());
    cmf->OverrideSync(cg->periodicity());
  } else {
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse);
    cmf->OverrideSync();
  }
}

// Average down edge-centered MultiFab (nodal in 2 dirs, cell in 1)
// If crse_geom is provided, uses ParallelCopy with periodicity to propagate
// averaged valid cells to periodic partners, then FillBoundary for ghost cells.
void amrmfab_average_down_edge(void *fine_mf, void *crse_mf, void *crse_geom,
                               int *ref_ratio, int ngcrse) {
  auto *fmf = static_cast<amrex::MultiFab *>(fine_mf);
  auto *cmf = static_cast<amrex::MultiFab *>(crse_mf);
  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  int ncomp = cmf->nComp();

  // Average into temp, then copy to crse with periodicity
  amrex::BoxArray crse_fine_BA = amrex::coarsen(fmf->boxArray(), ratio);
  amrex::MultiFab ctmp(crse_fine_BA, fmf->DistributionMap(), ncomp, ngcrse);
  amrex::average_down_edges(*fmf, ctmp, ratio, ngcrse);

  // Copy averaged data to crse - with periodicity if geometry provided
  if (crse_geom) {
    auto *cg = static_cast<amrex::Geometry *>(crse_geom);
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse, cg->periodicity());
    cmf->FillBoundary(cg->periodicity());
    cmf->OverrideSync(cg->periodicity());
  } else {
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse);
    cmf->OverrideSync();
  }
}

// Average down node-centered MultiFab (nodal in all dirs)
// If crse_geom is provided, uses ParallelCopy with periodicity to propagate
// averaged valid cells to periodic partners, then FillBoundary for ghost cells.
void amrmfab_average_down_node(void *fine_mf, void *crse_mf, void *crse_geom,
                               int *ref_ratio, int ngcrse) {
  auto *fmf = static_cast<amrex::MultiFab *>(fine_mf);
  auto *cmf = static_cast<amrex::MultiFab *>(crse_mf);
  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  int ncomp = cmf->nComp();

  // Average into temp, then copy to crse with periodicity
  amrex::BoxArray crse_fine_BA = amrex::coarsen(fmf->boxArray(), ratio);
  amrex::MultiFab ctmp(crse_fine_BA, fmf->DistributionMap(), ncomp, ngcrse);
  amrex::average_down_nodal(*fmf, ctmp, ratio, ngcrse);

  // Copy averaged data to crse - with periodicity if geometry provided
  if (crse_geom) {
    auto *cg = static_cast<amrex::Geometry *>(crse_geom);
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse, cg->periodicity());
    cmf->FillBoundary(cg->periodicity());
    cmf->OverrideSync(cg->periodicity());
  } else {
    cmf->ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse);
    cmf->OverrideSync();
  }
}

// Compute divergence of face-centered velocity into cell-centered MultiFab
// Uses AMReX's computeDivergence from MultiFabUtil
void amrmfab_compute_divergence(void *divu, void *umac_x, void *umac_y,
                                void *umac_z, void *geom) {
  auto *div = static_cast<amrex::MultiFab *>(divu);
  auto const *ux = static_cast<amrex::MultiFab const *>(umac_x);
  auto const *uy = static_cast<amrex::MultiFab const *>(umac_y);
  auto const *uz = static_cast<amrex::MultiFab const *>(umac_z);
  auto const *g = static_cast<amrex::Geometry const *>(geom);

  amrex::Array<amrex::MultiFab const *, AMREX_SPACEDIM> umac = {ux, uy, uz};
  amrex::computeDivergence(*div, umac, *g);
}

// =============================================================================
// 3-Component Face FillPatch Wrappers
// These wrappers call AMReX's coupled FillPatchTwoLevels/InterpFromCoarseLevel
// for face-centered data using the specified face interpolator.
// Note: Only face_divfree_interp (type=8) is currently supported for the
// coupled 3-component array versions - other face interps don't have interp_arr
// =============================================================================

// FillPatch from coarse level only (used during MakeNewLevelFromCoarse)
// Takes 3 MultiFabs (U,V,W) and fills them with interpolated coarse data
// interp_type: 8 = face_divfree (required for coupled divfree interpolation)

void amrmfab_fillcoarsepatch_faces(
    void *mf_u, void *mf_v, void *mf_w, double time, void *cmf_u, void *cmf_v,
    void *cmf_w, void *geom_c_ptr, void *geom_f_ptr, void *ctx_u, void *ctx_v,
    void *ctx_w, nga2::FillPatchBCDispatcher bc_u,
    nga2::FillPatchBCDispatcher bc_v, nga2::FillPatchBCDispatcher bc_w,
    int scomp, int dcomp, int ncomp, int *ref_ratio, int interp_type, int *lo_bc,
    int *hi_bc) {
  auto *u = static_cast<amrex::MultiFab *>(mf_u);
  auto *v = static_cast<amrex::MultiFab *>(mf_v);
  auto *w = static_cast<amrex::MultiFab *>(mf_w);
  auto *cu = static_cast<amrex::MultiFab *>(cmf_u);
  auto *cv = static_cast<amrex::MultiFab *>(cmf_v);
  auto *cw = static_cast<amrex::MultiFab *>(cmf_w);
  auto *geom_c = static_cast<amrex::Geometry *>(geom_c_ptr);
  auto *geom_f = static_cast<amrex::Geometry *>(geom_f_ptr);

  amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> mf_arr = {u, v, w};
  amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> cmf_arr = {cu, cv, cw};

  // Build BCRecs per dimension (each dimension has ncomp components)
  amrex::Array<amrex::Vector<amrex::BCRec>, AMREX_SPACEDIM> bcs;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = 0; i < ncomp; ++i) {
      int base = d * 3;
      bcs[d].emplace_back(lo_bc[base], lo_bc[base + 1], lo_bc[base + 2],
                          hi_bc[base], hi_bc[base + 1], hi_bc[base + 2]);
    }
  }

  // BC functors per dimension
  nga2::NGA2BCFunctor bc_functor_c_u(ctx_u, bc_u, geom_c);
  nga2::NGA2BCFunctor bc_functor_c_v(ctx_v, bc_v, geom_c);
  nga2::NGA2BCFunctor bc_functor_c_w(ctx_w, bc_w, geom_c);
  nga2::NGA2BCFunctor bc_functor_f_u(ctx_u, bc_u, geom_f);
  nga2::NGA2BCFunctor bc_functor_f_v(ctx_v, bc_v, geom_f);
  nga2::NGA2BCFunctor bc_functor_f_w(ctx_w, bc_w, geom_f);

  amrex::Array<nga2::NGA2BCFunctor, AMREX_SPACEDIM> cbc = {
      bc_functor_c_u, bc_functor_c_v, bc_functor_c_w};
  amrex::Array<nga2::NGA2BCFunctor, AMREX_SPACEDIM> fbc = {
      bc_functor_f_u, bc_functor_f_v, bc_functor_f_w};

  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));

  // Call with appropriate interpolator based on type
  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  if (interp_type == 8) {
    amrex::InterpFromCoarseLevel(mf_arr, time, cmf_arr, scomp - 1, dcomp - 1,
                                 ncomp, *geom_c, *geom_f, cbc, 0, fbc, 0, ratio,
                                 &amrex::face_divfree_interp, bcs, 0);
  } else if (interp_type == 9) {
    amrex::InterpFromCoarseLevel(mf_arr, time, cmf_arr, scomp - 1, dcomp - 1,
                                 ncomp, *geom_c, *geom_f, cbc, 0, fbc, 0, ratio,
                                 &amrex::face_linear_interp, bcs, 0);
  } else if (interp_type == 10) {
    amrex::InterpFromCoarseLevel(mf_arr, time, cmf_arr, scomp - 1, dcomp - 1,
                                 ncomp, *geom_c, *geom_f, cbc, 0, fbc, 0, ratio,
                                 &amrex::face_cons_linear_interp, bcs, 0);
  } else {
    amrex::Abort("amrmfab_fillcoarsepatch_faces: unsupported interp_type");
  }
}

// FillPatch from two levels (used during regular FillPatch operations)
// Takes 3 destination MultiFabs and source data from both coarse and fine
// levels
void amrmfab_fillpatch_two_faces(
    void *mf_u, void *mf_v, void *mf_w, double time,
    // Coarse level: old and new states
    double time_old_c, void *mf_old_c_u, void *mf_old_c_v, void *mf_old_c_w,
    double time_new_c, void *mf_new_c_u, void *mf_new_c_v, void *mf_new_c_w,
    void *geom_c_ptr,
    // Fine level: old and new states
    double time_old_f, void *mf_old_f_u, void *mf_old_f_v, void *mf_old_f_w,
    double time_new_f, void *mf_new_f_u, void *mf_new_f_v, void *mf_new_f_w,
    void *geom_f_ptr,
    // BC callbacks and contexts
    void *ctx_u, void *ctx_v, void *ctx_w, nga2::FillPatchBCDispatcher bc_u,
    nga2::FillPatchBCDispatcher bc_v, nga2::FillPatchBCDispatcher bc_w,
    int scomp, int dcomp, int ncomp, int *ref_ratio, int interp_type, int *lo_bc,
    int *hi_bc, int nghost) {
  auto *u = static_cast<amrex::MultiFab *>(mf_u);
  auto *v = static_cast<amrex::MultiFab *>(mf_v);
  auto *w = static_cast<amrex::MultiFab *>(mf_w);
  auto *geom_c = static_cast<amrex::Geometry *>(geom_c_ptr);
  auto *geom_f = static_cast<amrex::Geometry *>(geom_f_ptr);

  // Destination array
  amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> mf_arr = {u, v, w};

  // Build coarse source arrays
  amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> cmf(2);
  cmf[0] = {static_cast<amrex::MultiFab *>(mf_old_c_u),
            static_cast<amrex::MultiFab *>(mf_old_c_v),
            static_cast<amrex::MultiFab *>(mf_old_c_w)};
  cmf[1] = {static_cast<amrex::MultiFab *>(mf_new_c_u),
            static_cast<amrex::MultiFab *>(mf_new_c_v),
            static_cast<amrex::MultiFab *>(mf_new_c_w)};
  amrex::Vector<amrex::Real> ctime = {time_old_c, time_new_c};

  // Build fine source arrays
  amrex::Vector<amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM>> fmf(2);
  fmf[0] = {static_cast<amrex::MultiFab *>(mf_old_f_u),
            static_cast<amrex::MultiFab *>(mf_old_f_v),
            static_cast<amrex::MultiFab *>(mf_old_f_w)};
  fmf[1] = {static_cast<amrex::MultiFab *>(mf_new_f_u),
            static_cast<amrex::MultiFab *>(mf_new_f_v),
            static_cast<amrex::MultiFab *>(mf_new_f_w)};
  amrex::Vector<amrex::Real> ftime = {time_old_f, time_new_f};

  // Build BCRecs per dimension
  amrex::Array<amrex::Vector<amrex::BCRec>, AMREX_SPACEDIM> bcs;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    for (int i = 0; i < ncomp; ++i) {
      int base = d * 3;
      bcs[d].emplace_back(lo_bc[base], lo_bc[base + 1], lo_bc[base + 2],
                          hi_bc[base], hi_bc[base + 1], hi_bc[base + 2]);
    }
  }

  // BC functors per dimension
  nga2::NGA2BCFunctor bc_functor_c_u(ctx_u, bc_u, geom_c);
  nga2::NGA2BCFunctor bc_functor_c_v(ctx_v, bc_v, geom_c);
  nga2::NGA2BCFunctor bc_functor_c_w(ctx_w, bc_w, geom_c);
  nga2::NGA2BCFunctor bc_functor_f_u(ctx_u, bc_u, geom_f);
  nga2::NGA2BCFunctor bc_functor_f_v(ctx_v, bc_v, geom_f);
  nga2::NGA2BCFunctor bc_functor_f_w(ctx_w, bc_w, geom_f);

  amrex::Array<nga2::NGA2BCFunctor, AMREX_SPACEDIM> cbc = {
      bc_functor_c_u, bc_functor_c_v, bc_functor_c_w};
  amrex::Array<nga2::NGA2BCFunctor, AMREX_SPACEDIM> fbc = {
      bc_functor_f_u, bc_functor_f_v, bc_functor_f_w};

  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));

  // Use explicit nghost if non-negative, otherwise use mf[0]->nGrowVect()
  amrex::IntVect ng = (nghost >= 0) ? amrex::IntVect(nghost) : u->nGrowVect();

  // Call with appropriate interpolator based on type
  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  if (interp_type == 8) {
    amrex::FillPatchTwoLevels(mf_arr, ng, time, cmf, ctime, fmf, ftime,
                              scomp - 1, dcomp - 1, ncomp, *geom_c, *geom_f,
                              cbc, 0, fbc, 0, ratio,
                              &amrex::face_divfree_interp, bcs, 0);
  } else if (interp_type == 9) {
    amrex::FillPatchTwoLevels(mf_arr, ng, time, cmf, ctime, fmf, ftime,
                              scomp - 1, dcomp - 1, ncomp, *geom_c, *geom_f,
                              cbc, 0, fbc, 0, ratio,
                              &amrex::face_linear_interp, bcs, 0);
  } else if (interp_type == 10) {
    amrex::FillPatchTwoLevels(mf_arr, ng, time, cmf, ctime, fmf, ftime,
                              scomp - 1, dcomp - 1, ncomp, *geom_c, *geom_f,
                              cbc, 0, fbc, 0, ratio,
                              &amrex::face_cons_linear_interp, bcs, 0);
  } else {
    amrex::Abort("amrmfab_fillpatch_two_faces: unsupported interp_type");
  }
}

// =====================================================================
// Per-direction wrappers for AMReX APIs whose Fortran interfaces
// only accept scalar ref_ratio. These call the C++ APIs directly
// with IntVect for full per-direction support.
// =====================================================================

// FluxRegister: build with per-direction ref_ratio
void amrfluxreg_build(void **fr_ptr, void *ba_ptr, void *dm_ptr,
                      int *ref_ratio, int fine_lev, int ncomp) {
  auto *ba = static_cast<amrex::BoxArray *>(ba_ptr);
  auto *dm = static_cast<amrex::DistributionMapping *>(dm_ptr);
  amrex::IntVect rr(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  auto *fr = new amrex::FluxRegister(*ba, *dm, rr, fine_lev, ncomp);
  *fr_ptr = static_cast<void *>(fr);
}

// FluxRegister: destroy
void amrfluxreg_destroy(void *fr_ptr) {
  auto *fr = static_cast<amrex::FluxRegister *>(fr_ptr);
  delete fr;
}

// average_down_faces with per-direction ref_ratio
// Takes 3 fine + 3 coarse MultiFab pointers (x,y,z) + geometry + ratio
void amrmfab_average_down_faces(void *fine_x, void *fine_y, void *fine_z,
                                void *crse_x, void *crse_y, void *crse_z,
                                void *geom_ptr, int scomp, int ncomp,
                                int *ref_ratio) {
  auto *fx = static_cast<amrex::MultiFab *>(fine_x);
  auto *fy = static_cast<amrex::MultiFab *>(fine_y);
  auto *fz = static_cast<amrex::MultiFab *>(fine_z);
  auto *cx = static_cast<amrex::MultiFab *>(crse_x);
  auto *cy = static_cast<amrex::MultiFab *>(crse_y);
  auto *cz = static_cast<amrex::MultiFab *>(crse_z);
  auto *geom = static_cast<amrex::Geometry *>(geom_ptr);
  amrex::IntVect rr(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  amrex::Array<amrex::MultiFab const *, AMREX_SPACEDIM> fine{
      AMREX_D_DECL(fx, fy, fz)};
  amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> crse{
      AMREX_D_DECL(cx, cy, cz)};
  amrex::average_down_faces(fine, crse, rr, *geom);
}

// linop set_coarse_fine_bc with per-direction ref_ratio
void amrlinop_set_coarse_fine_bc(void *linop_ptr, void *crse_mf_ptr,
                                 int *ref_ratio) {
  auto *linop = static_cast<amrex::MLLinOp *>(linop_ptr);
  auto *crse = static_cast<amrex::MultiFab const *>(crse_mf_ptr);
  amrex::IntVect rr(AMREX_D_DECL(ref_ratio[0], ref_ratio[1], ref_ratio[2]));
  linop->setCoarseFineBC(crse, rr);
}

// =============================================================================
// Linear Operator Builders
// Mirror AMReX's amrex_fi_new_poisson/abeclaplacian with additional options:
//   semicoarsen:      0=off, nonzero=on (anisotropic MG coarsening)
//   hidden_direction: -1=none, 0/1/2=ignore that direction (quasi-2D)
// =============================================================================

// Build MLPoisson
void amrpoisson_build_c(amrex::MLLinOp *&linop, int nlevels,
                              const amrex::Geometry *geom[],
                              const amrex::BoxArray *ba[],
                              const amrex::DistributionMapping *dm[],
                              int semicoarsen, int hidden_direction) {
  amrex::LPInfo info;
  info.setAgglomeration(true);
  info.setConsolidation(true);
  info.setMaxCoarseningLevel(30);
  if (semicoarsen != 0) {
    info.setSemicoarsening(true);
    info.setMaxSemicoarseningLevel(100);
  }
  if (hidden_direction >= 0) {
    info.setHiddenDirection(hidden_direction);
  }
  amrex::Vector<amrex::Geometry> g;
  amrex::Vector<amrex::BoxArray> b;
  amrex::Vector<amrex::DistributionMapping> d;
  for (int i = 0; i < nlevels; ++i) {
    g.push_back(*geom[i]);
    b.push_back(*ba[i]);
    d.push_back(*dm[i]);
  }
  auto *poisson = new amrex::MLPoisson(g, b, d, info);
  linop = static_cast<amrex::MLLinOp *>(poisson);
}

// Build MLABecLaplacian
void amrabeclap_build_c(amrex::MLLinOp *&linop, int nlevels,
                              const amrex::Geometry *geom[],
                              const amrex::BoxArray *ba[],
                              const amrex::DistributionMapping *dm[],
                              int semicoarsen, int hidden_direction) {
  amrex::LPInfo info;
  info.setAgglomeration(true);
  info.setConsolidation(true);
  info.setMaxCoarseningLevel(30);
  if (semicoarsen != 0) {
    info.setSemicoarsening(true);
    info.setMaxSemicoarseningLevel(100);
  }
  if (hidden_direction >= 0) {
    info.setHiddenDirection(hidden_direction);
  }
  amrex::Vector<amrex::Geometry> g;
  amrex::Vector<amrex::BoxArray> b;
  amrex::Vector<amrex::DistributionMapping> d;
  for (int i = 0; i < nlevels; ++i) {
    g.push_back(*geom[i]);
    b.push_back(*ba[i]);
    d.push_back(*dm[i]);
  }
  auto *abeclap = new amrex::MLABecLaplacian(g, b, d, info);
  linop = static_cast<amrex::MLLinOp *>(abeclap);
}

//-----------------------------------------------------------------------------
// ParallelCopy with ADD semantics (for deposit accumulation across DMs)
//-----------------------------------------------------------------------------

void amrmfab_parallel_add(void *dst_ptr, void *src_ptr,
                          int srccomp, int dstcomp, int ncomp,
                          int srcng, int dstng, void *geom_ptr) {
  auto *dst  = static_cast<amrex::MultiFab *>(dst_ptr);
  auto *src  = static_cast<amrex::MultiFab *>(src_ptr);
  auto *geom = static_cast<amrex::Geometry *>(geom_ptr);
  dst->ParallelCopy(*src, srccomp, dstcomp, ncomp, srcng, dstng,
                    geom->periodicity(), amrex::FabArrayBase::ADD);
}

} // extern "C"