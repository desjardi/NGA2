// amrpd_wrapper.cpp
// C++ bridge between Fortran amrpd_class and AMReX particle containers.
//
// Two containers:
//   - AMRPDPC: NeighborParticleContainer<13,1> for solid particles
//               * 13 extra reals: vel[3], F_bond[3], F_fluid[3], mw, dil, damage, nb0
//               * 1  extra int  : flag
//               * setEnableInverse(true) so fillNeighbors records inverse_tags,
//                 enabling sumNeighbors for ghost-to-owner force/state reductions
//   - AMRPDBC: ParticleContainer<4,5> for bonds (no neighbor machinery)
//               * 4 extra reals: d0, w, damage, hist1 (reserved)
//               * 5 extra ints : id_lo_lo, id_lo_hi, id_hi_lo, id_hi_hi, alive
//
// Bonds are one-sided: each physical bond is stored once, with bond.pos anchored
// at the position of its lower-GID endpoint. AMReX's position-based Redistribute
// then keeps bond and lower-GID particle co-located without an explicit GID->rank
// map. The higher-GID endpoint is accessed via the particle ghost layer.

#include <AMReX_NeighborParticles.H>
#include <AMReX_NeighborList.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <fstream>
#include <iomanip>

using namespace amrex;

// Layout constants -- MUST match Fortran AMRPD_NREAL_PART/NINT_PART/NREAL_BOND/NINT_BOND
#define AMRPD_NREAL_PART 13
#define AMRPD_NINT_PART   1
#define AMRPD_NREAL_BOND  4
#define AMRPD_NINT_BOND   5

namespace {

// -----------------------------------------------------------------------
// Particle container (NeighborParticleContainer with sumNeighbors enabled)
// -----------------------------------------------------------------------
class AMRPDPC : public NeighborParticleContainer<AMRPD_NREAL_PART, AMRPD_NINT_PART>
{
public:
    using Base = NeighborParticleContainer<AMRPD_NREAL_PART, AMRPD_NINT_PART>;
    using ParticleType = Base::ParticleType;

    explicit AMRPDPC(AmrCore* amrcore)
        : Base(amrcore->GetParGDB(), 1)
    {
        // Record (orig_grid, orig_tile, orig_index, orig_level) for each ghost
        // during fillNeighbors so sumNeighbors can scatter-add back to owners
        setEnableInverse(true);
    }

    // Allow changing ghost cell width before fillNeighbors() (uniform per call)
    void setNeighborCells(int n) { m_num_neighbor_cells = n; }

    // CSR access to the neighbor list for a tile (built by buildNeighborList).
    // m_neighbor_list is protected on the base class; expose it here.
    NeighborList<ParticleType>& getNeighborList(int lev, int grid, int tile)
    {
        return m_neighbor_list[lev][std::make_pair(grid, tile)];
    }
};

// -----------------------------------------------------------------------
// Bond container (plain ParticleContainer -- bonds don't have neighbors)
// -----------------------------------------------------------------------
class AMRPDBC : public ParticleContainer<AMRPD_NREAL_BOND, AMRPD_NINT_BOND>
{
public:
    using Base = ParticleContainer<AMRPD_NREAL_BOND, AMRPD_NINT_BOND>;
    using ParticleType = Base::ParticleType;

    explicit AMRPDBC(AmrCore* amrcore)
        : Base(amrcore->GetParGDB())
    {}
};

using PCP = AMRPDPC;
using PTP = PCP::ParticleType;
using PCB = AMRPDBC;
using PTB = PCB::ParticleType;

} // namespace

extern "C" {

// -----------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------

void amrpd_new_pcp(PCP*& pc, void* amrcore_raw)
{
    pc = new PCP(static_cast<AmrCore*>(amrcore_raw));
}

void amrpd_delete_pcp(PCP* pc)
{
    delete pc;
}

void amrpd_new_pcb(PCB*& pc, void* amrcore_raw)
{
    pc = new PCB(static_cast<AmrCore*>(amrcore_raw));
}

void amrpd_delete_pcb(PCB* pc)
{
    delete pc;
}

// -----------------------------------------------------------------------
// Redistribute (AMR-aware position-based migration)
// -----------------------------------------------------------------------

void amrpd_redistribute_p(PCP* pc, int lev_min, int lev_max, int ng)
{
    pc->Redistribute(lev_min, lev_max, ng);
    // Invalidate cached neighbor/inverse-tag state. After Redistribute
    // (especially across a regrid where the BA/DM changes), inverse_tags
    // still reference grid indices in the OLD DistributionMap. Any
    // subsequent sumNeighbors call would read ParticleDistributionMap(lev)
    // out-of-bounds. clearNeighbors() drops the stale state so the next
    // fillNeighbors() rebuilds fresh against the new BA/DM. amrlpt doesn't
    // need this because it never enables inverse_tags; we do.
    pc->clearNeighbors();
}

void amrpd_redistribute_b(PCB* pc, int lev_min, int lev_max, int ng)
{
    pc->Redistribute(lev_min, lev_max, ng);
}

// -----------------------------------------------------------------------
// Ghost (neighbor) particle exchange on the PARTICLE container
//
// fillNeighbors_radius ghosts particles within physical distance search_radius
// of tile boundaries (much less data than cell-based when r << dx).
// sumNeighbors performs the inverse: adds ghost-particle rdata/idata back into
// the corresponding real particles on their owner ranks (intra- and inter-rank).
// Requires setEnableInverse(true) at construction (already set in AMRPDPC ctor).
// -----------------------------------------------------------------------

void amrpd_fill_neighbors_radius(PCP* pc, double r)
{
    pc->fillNeighbors(amrex::Real(r));
}

void amrpd_clear_neighbors(PCP* pc)
{
    pc->clearNeighbors();
}

// Refresh existing ghost data with current owner values (topology unchanged).
// Cheaper than fillNeighbors -- skips mask rebuild and tag re-cache. Use to
// propagate updated owner state (e.g., theta after compute_dilatation) into
// the ghost slots that the next kernel will read.
void amrpd_update_neighbors(PCP* pc)
{
    pc->updateNeighbors();
}

// -----------------------------------------------------------------------
// Neighbor list (explicit pair list) for short-range contact.
// buildNeighborList: must be called AFTER fillNeighbors. Cell-binned, gives a
// per-tile CSR (offsets, list) of (i, j) pairs with |r_i - r_j| < rcrit.
//   check_pair criterion: |r_i - r_j|^2 < rcrit^2
//   offsets[0..np]: prefix sums (unsigned int)
//   list[0..ntot-1]: neighbor indices into combined real+ghost particle array
//   np: number of valid particles (driver loop runs i in [0, np-1])
//   ntot: total neighbor entries across all i
// -----------------------------------------------------------------------

void amrpd_build_neighbor_list(PCP* pc, double rcrit)
{
    const double rcrit2 = rcrit * rcrit;
    auto check_pair = [rcrit2](const PTP& p1, const PTP& p2) -> bool {
        double d2 = 0.0;
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
            d2 += (p1.pos(dim) - p2.pos(dim)) * (p1.pos(dim) - p2.pos(dim));
        return d2 < rcrit2;
    };
    pc->buildNeighborList(check_pair, amrex::Real(rcrit), false);
}

void amrpd_get_neighbor_list_mfi(PCP* pc, int lev, MFIter* mfi,
                                  const unsigned int*& offsets,
                                  const unsigned int*& list,
                                  long long& np, long long& ntot)
{
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& nl = pc->getNeighborList(lev, grid, tile);
    np   = static_cast<long long>(nl.numParticles());
    ntot = static_cast<long long>(nl.GetList().size());
    offsets = (np   > 0) ? nl.GetOffsets().dataPtr() : nullptr;
    list    = (ntot > 0) ? nl.GetList().dataPtr()    : nullptr;
}

// Reduce ghost rdata[rs..rs+rn-1] and idata[is..is+in-1] back to owners
void amrpd_sum_neighbors(PCP* pc, int rs, int rn, int is, int in)
{
    pc->sumNeighbors(rs, rn, is, in);
}

// -----------------------------------------------------------------------
// MFIter accessors -- particles (real, ghost, combined)
// -----------------------------------------------------------------------

void amrpd_get_particles_mfi(PCP* pc, int lev, MFIter* mfi,
                              PTP*& dp, long long& np)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) {
        np = 0; dp = nullptr; return;
    }
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& plev = pc->GetParticles(lev);
    auto it = plev.find(std::make_pair(grid, tile));
    if (it != plev.end()) {
        auto& ptile = it->second;
        np = static_cast<long long>(ptile.numRealParticles());
        dp = (np > 0) ? ptile.GetArrayOfStructs().data() : nullptr;
    } else {
        np = 0; dp = nullptr;
    }
}

void amrpd_get_all_particles_mfi(PCP* pc, int lev, MFIter* mfi,
                                  PTP*& dp, long long& np_total, long long& np_valid)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) {
        np_total = 0; np_valid = 0; dp = nullptr; return;
    }
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& plev = pc->GetParticles(lev);
    auto it = plev.find(std::make_pair(grid, tile));
    if (it != plev.end()) {
        auto& ptile = it->second;
        np_valid = static_cast<long long>(ptile.numRealParticles());
        np_total = static_cast<long long>(ptile.numTotalParticles());
        dp = (np_total > 0) ? ptile.GetArrayOfStructs().data() : nullptr;
    } else {
        np_total = 0; np_valid = 0; dp = nullptr;
    }
}

// Returns a writable pointer to the ghost-particle slice of this tile. This
// targets AMReX's internal neighbor buffer (pc->GetNeighbors), NOT the ghost
// copies appended to the ptile by fillNeighbors. The two hold identical data
// right after fillNeighbors (parallel arrays, identical ordering), but only
// the internal buffer is what sumNeighbors reads from when reducing back to
// owners — writes to the ptile's appended ghost slots would be discarded.
void amrpd_get_ghosts_mfi(PCP* pc, int lev, MFIter* mfi,
                           PTP*& dp, long long& ng)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) {
        ng = 0; dp = nullptr; return;
    }
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& nbuf = pc->GetNeighbors(lev, grid, tile);
    ng = static_cast<long long>(nbuf.numParticles());
    dp = (ng > 0) ? nbuf.GetArrayOfStructs().data() : nullptr;
}

// -----------------------------------------------------------------------
// MFIter accessor -- bonds
// -----------------------------------------------------------------------

void amrpd_get_bonds_mfi(PCB* pc, int lev, MFIter* mfi,
                          PTB*& dp, long long& nb)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) {
        nb = 0; dp = nullptr; return;
    }
    const int grid = mfi->index();
    const int tile = mfi->LocalTileIndex();
    auto& plev = pc->GetParticles(lev);
    auto it = plev.find(std::make_pair(grid, tile));
    if (it != plev.end()) {
        auto& ptile = it->second;
        nb = static_cast<long long>(ptile.numRealParticles());
        dp = (nb > 0) ? ptile.GetArrayOfStructs().data() : nullptr;
    } else {
        nb = 0; dp = nullptr;
    }
}

// -----------------------------------------------------------------------
// Add a single particle or bond to a specific (level, grid, tile)
// -----------------------------------------------------------------------

void amrpd_add_particle_i(PCP* pc, int lev, int grid, int tile, PTP* p)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) return;
    auto& plev = pc->GetParticles(lev);
    plev[std::make_pair(grid, tile)].push_back(*p);
}

void amrpd_add_bond_i(PCB* pc, int lev, int grid, int tile, PTB* b)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) return;
    auto& plev = pc->GetParticles(lev);
    plev[std::make_pair(grid, tile)].push_back(*b);
}

// -----------------------------------------------------------------------
// Bulk append an array of particles to the container at level 0.
// Collective: all ranks must call. Ranks with nothing to add pass n=0
// and raw=nullptr. Each appended particle gets a fresh AMReX (id,cpu)
// and AddParticlesAtLevel internally redistributes by spatial position.
// -----------------------------------------------------------------------

void amrpd_append_particles(PCP* pc, const void* raw, long long n)
{
    PCP::ParticleTileType ptile;

    if (n > 0 && raw != nullptr) {
        const PTP* src = reinterpret_cast<const PTP*>(raw);
        ptile.resize(static_cast<int>(n));
        auto& aos = ptile.GetArrayOfStructs();
        for (long long i = 0; i < n; ++i) {
            PTP p = src[i];                         // copy pos + rdata + idata
            p.id()  = PTP::NextID();                // unique AMReX identity
            p.cpu() = ParallelDescriptor::MyProc();
            aos[static_cast<int>(i)] = p;
        }
    }

    pc->AddParticlesAtLevel(ptile, 0);              // collective: level 0 redistribute inside
}

// Bulk append of bonds (mirror of amrpd_append_particles for the bond container).
// Collective: all ranks must call. Each appended bond gets a fresh AMReX (id,cpu).
// AddParticlesAtLevel only redistributes within the single level passed; we
// follow with a full-hierarchy Redistribute() so bonds migrate to the finest
// AMR level containing their pos (= position of the lower-GID endpoint), keeping
// the bond co-located with its lower-GID owner across all refinement levels.
void amrpd_append_bonds(PCB* pc, const void* raw, long long n)
{
    PCB::ParticleTileType ptile;

    if (n > 0 && raw != nullptr) {
        const PTB* src = reinterpret_cast<const PTB*>(raw);
        ptile.resize(static_cast<int>(n));
        auto& aos = ptile.GetArrayOfStructs();
        for (long long i = 0; i < n; ++i) {
            PTB b = src[i];
            b.id()  = PTB::NextID();
            b.cpu() = ParallelDescriptor::MyProc();
            aos[static_cast<int>(i)] = b;
        }
    }

    pc->AddParticlesAtLevel(ptile, 0);
    pc->Redistribute();
}

// -----------------------------------------------------------------------
// BoxArray / DistributionMap accessors (for post-regrid sync)
// -----------------------------------------------------------------------

void amrpd_get_particle_boxarray_p(PCP* pc, int lev, void** ba_ptr)
{
    *ba_ptr = const_cast<amrex::BoxArray*>(&(pc->ParticleBoxArray(lev)));
}

void amrpd_get_particle_distromap_p(PCP* pc, int lev, void** dm_ptr)
{
    *dm_ptr = const_cast<amrex::DistributionMapping*>(&(pc->ParticleDistributionMap(lev)));
}

void amrpd_set_particle_boxarray_p(PCP* pc, int lev, void* ba_ptr)
{
    pc->SetParticleBoxArray(lev, *static_cast<amrex::BoxArray*>(ba_ptr));
}

void amrpd_set_particle_distromap_p(PCP* pc, int lev, void* dm_ptr)
{
    pc->SetParticleDistributionMap(lev, *static_cast<amrex::DistributionMapping*>(dm_ptr));
}

void amrpd_get_particle_boxarray_b(PCB* pc, int lev, void** ba_ptr)
{
    *ba_ptr = const_cast<amrex::BoxArray*>(&(pc->ParticleBoxArray(lev)));
}

void amrpd_get_particle_distromap_b(PCB* pc, int lev, void** dm_ptr)
{
    *dm_ptr = const_cast<amrex::DistributionMapping*>(&(pc->ParticleDistributionMap(lev)));
}

void amrpd_set_particle_boxarray_b(PCB* pc, int lev, void* ba_ptr)
{
    pc->SetParticleBoxArray(lev, *static_cast<amrex::BoxArray*>(ba_ptr));
}

void amrpd_set_particle_distromap_b(PCB* pc, int lev, void* dm_ptr)
{
    pc->SetParticleDistributionMap(lev, *static_cast<amrex::DistributionMapping*>(dm_ptr));
}

// -----------------------------------------------------------------------
// ID and CPU counters / accessors
// -----------------------------------------------------------------------

void amrpd_get_next_id_p(long long& id) { id = PTP::NextID(); }
void amrpd_set_next_id_p(long long id)  { PTP::NextID(id); }
void amrpd_get_next_id_b(long long& id) { id = PTB::NextID(); }
void amrpd_set_next_id_b(long long id)  { PTB::NextID(id); }

void amrpd_get_cpu(int& cpu) { cpu = ParallelDescriptor::MyProc(); }

void amrpd_get_particle_id(long long& id, const PTP* p) { id = p->id(); }
void amrpd_set_particle_id(long long id, PTP* p)        { p->id() = id; }
void amrpd_get_particle_cpu(int& cpu, const PTP* p)     { cpu = p->cpu(); }
void amrpd_set_particle_cpu(int cpu, PTP* p)            { p->cpu() = cpu; }

// Unique 64-bit key composed from (id, cpu) for the Fortran-side GID hash.
// Encoding is internal — only needs stability across one run. 40 bits for
// id (up to 1T particles per rank) and 24 bits for cpu (16M ranks); both
// well above any realistic simulation scale.
void amrpd_get_particle_idcpu(long long& key, const PTP* p)
{
    key = (static_cast<long long>(p->cpu()) << 40) |
          (static_cast<long long>(p->id()) & 0xFFFFFFFFFFLL);
}

void amrpd_get_bond_id(long long& id, const PTB* b)     { id = b->id(); }
void amrpd_set_bond_id(long long id, PTB* b)            { b->id() = id; }
void amrpd_get_bond_cpu(int& cpu, const PTB* b)         { cpu = b->cpu(); }
void amrpd_set_bond_cpu(int cpu, PTB* b)                { b->cpu() = cpu; }

// -----------------------------------------------------------------------
// Total counts (global, across all ranks and levels)
// -----------------------------------------------------------------------

void amrpd_total_np(PCP* pc, long long& np) { np = pc->TotalNumberOfParticles(); }
void amrpd_total_nb(PCB* pc, long long& nb) { nb = pc->TotalNumberOfParticles(); }

// -----------------------------------------------------------------------
// Checkpoint I/O (one subdirectory per container under a shared checkpoint
// directory: <dir>/particles/ and <dir>/bonds/). Caller is responsible for
// creating <dir>; AMReX's Checkpoint() builds the subdirectory.
// -----------------------------------------------------------------------

// Split a single fullpath into (parent, leaf) the same way amrlpt does, so the
// Fortran caller passes one composed path argument per call.
static inline void split_path(const char* fullpath,
                               std::string& parent, std::string& leaf)
{
    std::string path(fullpath);
    while (!path.empty() && path.back() == '/') path.pop_back();
    auto pos = path.rfind('/');
    parent = (pos != std::string::npos) ? path.substr(0, pos) : std::string(".");
    leaf   = (pos != std::string::npos) ? path.substr(pos+1) : path;
}

void amrpd_checkpoint_p(PCP* pc, const char* fullpath)
{
    std::string parent, leaf;
    split_path(fullpath, parent, leaf);
    pc->Checkpoint(parent, leaf, true);
}

void amrpd_restart_p(PCP* pc, const char* fullpath)
{
    std::string parent, leaf;
    split_path(fullpath, parent, leaf);
    pc->Restart(parent, leaf);
}

void amrpd_checkpoint_b(PCB* pc, const char* fullpath)
{
    std::string parent, leaf;
    split_path(fullpath, parent, leaf);
    pc->Checkpoint(parent, leaf, true);
}

void amrpd_restart_b(PCB* pc, const char* fullpath)
{
    std::string parent, leaf;
    split_path(fullpath, parent, leaf);
    pc->Restart(parent, leaf);
}

// -----------------------------------------------------------------------
// Visualization plotfile for the PARTICLE container (bonds not written; we
// can add a separate bond writer if and when we need bond visualization).
//
// write_real[AMRPD_NREAL_PART]: bitmask (1=write, 0=skip) per extra real.
// write_int[AMRPD_NINT_PART]:   bitmask per extra int.
// Component names are hardcoded here to match the Fortran part struct layout:
//   rdata[0..2]   vel      -> vx,  vy,  vz
//   rdata[3..5]   F_bond   -> fbx, fby, fbz
//   rdata[6..8]   F_fluid  -> ffx, ffy, ffz
//   rdata[9]      mw
//   rdata[10]     dil
//   rdata[11]     damage
//   rdata[12]     nb0
//   idata[0]      flag
// Position (pos[3]) is always written by AMReX (baked into the particle format).
// -----------------------------------------------------------------------

void amrpd_write_plotfile(PCP* pc, const char* basedir, const char* pname,
                           const int* write_real, const int* write_int, double time)
{
    static const Vector<std::string> rnames = {
        "vx", "vy", "vz",
        "fbx", "fby", "fbz",
        "ffx", "ffy", "ffz",
        "mw", "dil", "damage", "nb0"
    };
    static const Vector<std::string> inames = { "flag" };

    Vector<int> wr(write_real, write_real + AMRPD_NREAL_PART);
    Vector<int> wi(write_int,  write_int  + AMRPD_NINT_PART);

    pc->WritePlotFile(std::string(basedir), std::string(pname),
                      wr, wi, rnames, inames);

    // Store simulation time inside the plotfile directory (IOProcessor only).
    // Used by amrpd_read_plotfile_time on restart to recover the time series.
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream tf(std::string(basedir) + "/time");
        tf << std::setprecision(17) << time << '\n';
    }
}

// -----------------------------------------------------------------------
// Read back the simulation time stored by amrpd_write_plotfile.
// Returns -1.0 if the file does not exist or cannot be read.
// -----------------------------------------------------------------------

double amrpd_read_plotfile_time(const char* basedir)
{
    std::ifstream tf(std::string(basedir) + "/time");
    double t = -1.0;
    if (tf.good()) tf >> t;
    return t;
}

} // extern "C"
