// amrlpt_wrapper.cpp
// C++ bridge between Fortran amrlpt_class and AMReX NeighborParticleContainer
// Particle struct: 14 extra reals (d, vel[3], angVel[3], Acol[3], Tcol[3], dt)
//                  1 extra int   (flag)
// All functions callable from Fortran via bind(C)
//
// Uses NeighborParticleContainer (backward-compatible superset of AmrParticleContainer):
//   - fillNeighbors / clearNeighbors: ghost particles for collision detection
//   - buildNeighborList: explicit pair lists for DEM / peridynamics

#include <AMReX_NeighborParticleContainer.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

// -----------------------------------------------------------------------
// Particle container type: 14 extra reals, 1 extra int
// Full in-memory layout per particle (double precision):
//   pos[3]     (24 B)  -- AMReX base
//   rdata[14]  (112 B) -- d, vel[3], angVel[3], Acol[3], Tcol[3], dt
//   idcpu      (8 B)   -- AMReX packed id+cpu
//   idata[1]   (4 B)   -- flag
// -----------------------------------------------------------------------
#define AMRLPT_NREAL 14
#define AMRLPT_NINT   1

namespace {
    // NeighborParticleContainer is a backward-compatible superset of AmrParticleContainer.
    // Neighbor particles have the same layout as primary particles (NNeighborReal=NREAL, NNeighborInt=NINT).
    using FPC = NeighborParticleContainer<AMRLPT_NREAL, AMRLPT_NINT>;
    using FPT = FPC::ParticleType;
}

extern "C" {

// -----------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------

void amrlpt_new_pc(FPC*& pc, void* amrcore_raw)
{
    // amrcore_raw is a NGA2AmrCore* (is-a AmrCore*) stored as void* via Fortran c_ptr.
    // Single public inheritance => same address; safe to cast directly.
    pc = new FPC(static_cast<AmrCore*>(amrcore_raw));
}

void amrlpt_delete_pc(FPC* pc)
{
    delete pc;
}

// -----------------------------------------------------------------------
// Redistribution (AMR-aware particle sorting to finest covering level)
// -----------------------------------------------------------------------

void amrlpt_redistribute(FPC* pc, int lev_min, int lev_max, int ng)
{
    pc->Redistribute(lev_min, lev_max, ng);
}

// -----------------------------------------------------------------------
// Neighbor/ghost particles for collision detection and short-range interactions
// fillNeighbors: communicate particles within ngrow cells into neighbor buffer
// clearNeighbors: release neighbor buffer
// -----------------------------------------------------------------------

void amrlpt_fill_neighbors(FPC* pc, int ngrow)
{
    pc->fillNeighbors(ngrow);
}

void amrlpt_clear_neighbors(FPC* pc)
{
    pc->clearNeighbors();
}

// -----------------------------------------------------------------------
// Neighbor particle access per tile (read-only ghost copies).
// Returns pointer to neighbor particle array + count.
// Call after amrlpt_fill_neighbors; neighbor particles may overlap
// valid particles from adjacent tiles.
// -----------------------------------------------------------------------

void amrlpt_get_neighbor_particles_mfi(FPC* pc, int lev, MFIter* mfi,
                                        FPT*& dp, long long& np)
{
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& neighbors = pc->GetNeighbors(lev, grid, tile);
    np = static_cast<long long>(neighbors.numParticles());
    dp = (np > 0) ? neighbors.GetArrayOfStructs().data() : nullptr;
}

// -----------------------------------------------------------------------
// Neighbor list (explicit pair list) for DEM / peridynamics.
// buildNeighborList: build pairs with |r_i - r_j| < rcrit using
//   cell-linked-list search over neighbor particles.
// Pairs are stored as flat int arrays (2*npairs): [i0,j0, i1,j1, ...]
// where i is index into valid particles and j into neighbor particles.
// -----------------------------------------------------------------------

void amrlpt_build_neighbor_list(FPC* pc, double rcrit)
{
    const double rcrit2 = rcrit * rcrit;
    auto check_pair = [rcrit2](const FPT& p1, const FPT& p2) -> bool {
        double d2 = 0.0;
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
            d2 += (p1.pos(dim) - p2.pos(dim)) * (p1.pos(dim) - p2.pos(dim));
        return d2 < rcrit2;
    };
    // BuildNeighborList populates GetNeighborList per tile
    bool sort = false;
    pc->buildNeighborList(check_pair, sort);
}

void amrlpt_get_neighbor_list_mfi(FPC* pc, int lev, MFIter* mfi,
                                   const int*& pairs, long long& npairs)
{
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& nl = pc->getNeighborList(lev, grid, tile);
    npairs = static_cast<long long>(nl.size() / 2);  // each pair is (i,j)
    pairs  = (npairs > 0) ? nl.dataPtr() : nullptr;
}

// -----------------------------------------------------------------------
// Particle access via MFIter (grid+tile indices come from the MFIter)
// Returns pointer to contiguous particle array + count for this tile
// -----------------------------------------------------------------------

void amrlpt_get_particles_mfi(FPC* pc, int lev, MFIter* mfi,
                               FPT*& dp, long long& np)
{
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& plev = pc->GetParticles(lev);
    auto it = plev.find(std::make_pair(grid, tile));
    if (it != plev.end()) {
        auto& ptile = it->second;
        np = ptile.numParticles();
        dp = (np > 0) ? ptile.GetArrayOfStructs().data() : nullptr;
    } else {
        np = 0;
        dp = nullptr;
    }
}

void amrlpt_num_particles_mfi(FPC* pc, int lev, MFIter* mfi, long long& np)
{
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& plev = pc->GetParticles(lev);
    auto it = plev.find(std::make_pair(grid, tile));
    np = (it != plev.end()) ? it->second.numParticles() : 0;
}

// -----------------------------------------------------------------------
// Add a single particle to a tile (by explicit grid+tile index)
// -----------------------------------------------------------------------

void amrlpt_add_particle_i(FPC* pc, int lev, int grid, int tile, FPT* p)
{
    auto& plev = pc->GetParticles(lev);
    plev[std::make_pair(grid, tile)].push_back(*p);
}

// -----------------------------------------------------------------------
// Global unique particle ID counter and current CPU rank
// -----------------------------------------------------------------------

void amrlpt_get_next_id(long long& id)
{
    id = FPT::NextID();
}

void amrlpt_set_next_id(long long id)
{
    FPT::NextID(id);
}

void amrlpt_get_cpu(int& cpu)
{
    cpu = ParallelDescriptor::MyProc();
}

// -----------------------------------------------------------------------
// ID/CPU accessors for the packed idcpu field (private in Fortran struct)
// -----------------------------------------------------------------------

void amrlpt_get_particle_id(long long& id, const FPT* p)
{
    id = p->id();
}

void amrlpt_set_particle_id(long long id, FPT* p)
{
    p->id() = id;
}

void amrlpt_get_particle_cpu(int& cpu, const FPT* p)
{
    cpu = p->cpu();
}

void amrlpt_set_particle_cpu(int cpu, FPT* p)
{
    p->cpu() = cpu;
}

void amrlpt_particle_is_valid(int& valid, const FPT* p)
{
    valid = p->id().is_valid() ? 1 : 0;
}

// -----------------------------------------------------------------------
// Total particle count (global, across all ranks and levels)
// -----------------------------------------------------------------------

void amrlpt_total_np(FPC* pc, long long& np)
{
    np = pc->TotalNumberOfParticles();
}

// -----------------------------------------------------------------------
// Checkpoint: write particles to fullpath (e.g. "restart/part_1.00E+00").
// Splits at last '/' to satisfy AMReX Checkpoint(dir, name) API.
// Creates parent directory hierarchy via PreBuildDirectorHierarchy.
// -----------------------------------------------------------------------

void amrlpt_write(FPC* pc, const char* fullpath, int is_chk)
{
    std::string path(fullpath);
    while (!path.empty() && path.back() == '/') path.pop_back();
    auto pos = path.rfind('/');
    std::string parent = (pos != std::string::npos) ? path.substr(0, pos) : std::string(".");
    std::string leaf   = (pos != std::string::npos) ? path.substr(pos+1) : path;
    amrex::PreBuildDirectorHierarchy(parent, "", 1, true);
    pc->Checkpoint(parent, leaf, is_chk != 0);
}

// -----------------------------------------------------------------------
// Visualization plotfile with selective component output.
// write_real[AMRLPT_NREAL]: bitmask (1=write, 0=skip) for each extra real.
// write_int[AMRLPT_NINT]:   bitmask for each extra int.
// Component names are hardcoded here to match the Fortran part struct layout:
//   rdata[0..13]: d, vx, vy, vz, wx, wy, wz, ax, ay, az, tx, ty, tz, dt
//   idata[0]:     flag
// Position (pos[3]) is always written by AMReX (baked into particle format).
// -----------------------------------------------------------------------

void amrlpt_write_plotfile(FPC* pc, const char* basedir, const char* pname,
                            const int* write_real, const int* write_int)
{
    static const std::vector<std::string> rnames = {
        "d", "vx", "vy", "vz", "wx", "wy", "wz",
        "ax", "ay", "az", "tx", "ty", "tz", "dt"
    };
    static const std::vector<std::string> inames = { "flag" };

    Vector<int> wr(write_real, write_real + AMRLPT_NREAL);
    Vector<int> wi(write_int,  write_int  + AMRLPT_NINT);

    pc->WritePlotFile(std::string(basedir), std::string(pname),
                      wr, wi, rnames, inames);
}

// -----------------------------------------------------------------------
// Restart: read back a particle checkpoint written by amrlpt_write.
// -----------------------------------------------------------------------

void amrlpt_read(FPC* pc, const char* dir, const char* name)
{
    pc->Restart(std::string(dir), std::string(name));
}

} // extern "C"
