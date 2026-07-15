// amrpd_wrapper.cpp
// C++ bridge between Fortran amrpd_class and AMReX particle containers.
//
// Two containers:
//   - AMRPDPC: plain ParticleContainer<15,1> holding the grid-side copy of
//     a pdsolver solid (positions/velocities/damage) for deposits, interp,
//     tagging, viz, and particle checkpoints. No physics, no ghosts.
//
// Bonds are one-sided: each physical bond is stored once, with bond.pos anchored
// at the position of its lower-GID endpoint. AMReX's position-based Redistribute
// then keeps bond and lower-GID particle co-located without an explicit GID->rank
// map. The higher-GID endpoint is accessed via the particle ghost layer.

#include <AMReX_Particles.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <fstream>
#include <iomanip>

using namespace amrex;

// Layout constants -- MUST match Fortran AMRPD_NREAL_PART/NINT_PART/NREAL_BOND/NINT_BOND
#define AMRPD_NREAL_PART 15
#define AMRPD_NINT_PART   1

namespace {

// -----------------------------------------------------------------------
// Particle container: plain AMReX ParticleContainer -- the grid-side face of
// a pdsolver solid. No neighbor/ghost machinery: physics lives in pdsolver.
// -----------------------------------------------------------------------
class AMRPDPC : public ParticleContainer<AMRPD_NREAL_PART, AMRPD_NINT_PART>
{
public:
    using Base = ParticleContainer<AMRPD_NREAL_PART, AMRPD_NINT_PART>;
    using ParticleType = Base::ParticleType;
    explicit AMRPDPC(AmrCore* amrcore) : Base(amrcore->GetParGDB()) {}
};

using PCP = AMRPDPC;
using PTP = PCP::ParticleType;

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



// -----------------------------------------------------------------------
// Redistribute (AMR-aware position-based migration)
// -----------------------------------------------------------------------

void amrpd_redistribute_p(PCP* pc, int lev_min, int lev_max, int ng)
{
    pc->Redistribute(lev_min, lev_max, ng);
}


// -----------------------------------------------------------------------




// -----------------------------------------------------------------------




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



// -----------------------------------------------------------------------
// MFIter accessor -- bonds
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// Add a single particle or bond to a specific (level, grid, tile)
// -----------------------------------------------------------------------

void amrpd_add_particle_i(PCP* pc, int lev, int grid, int tile, PTP* p)
{
    if (lev >= static_cast<int>(pc->GetParticles().size())) return;
    auto& plev = pc->GetParticles(lev);
    plev[std::make_pair(grid, tile)].push_back(*p);
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

// Append particles PRESERVING caller-provided identities: gids[i] is the RAW
// packed AMReX idcpu (valid bit | id | cpu) exactly as part_gid exposes it --
// restored verbatim, no decode. Used to rebuild the grid-side face from a
// pdsolver checkpoint, where identities must match the solver's node gids for
// exchange routing. NextID is bumped past the local max id so any later
// append cannot collide.
void amrpd_append_particles_gid(PCP* pc, const void* raw, long long n,
                                 const long long* gids)
{
    PCP::ParticleTileType ptile;
    long long maxid = 0;
    if (n > 0 && raw != nullptr) {
        const PTP* src = reinterpret_cast<const PTP*>(raw);
        ptile.resize(static_cast<int>(n));
        auto& aos = ptile.GetArrayOfStructs();
        for (long long i = 0; i < n; ++i) {
            PTP p = src[i];
            p.m_idcpu = static_cast<uint64_t>(gids[i]);
            const long long idv = p.id();
            if (idv > maxid) maxid = idv;
            aos[static_cast<int>(i)] = p;
        }
    }
    if (maxid >= PTP::NextID()) PTP::NextID(maxid + 1);
    pc->AddParticlesAtLevel(ptile, 0);
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





// -----------------------------------------------------------------------
// ID and CPU counters / accessors
// -----------------------------------------------------------------------

void amrpd_get_next_id_p(long long& id) { id = PTP::NextID(); }
void amrpd_set_next_id_p(long long id)  { PTP::NextID(id); }

void amrpd_get_cpu(int& cpu) { cpu = ParallelDescriptor::MyProc(); }

void amrpd_get_particle_id(long long& id, const PTP* p) { id = p->id(); }
void amrpd_set_particle_id(long long id, PTP* p)        { p->id() = id; }
void amrpd_get_particle_cpu(int& cpu, const PTP* p)     { cpu = p->cpu(); }
void amrpd_set_particle_cpu(int cpu, PTP* p)            { p->cpu() = cpu; }



// -----------------------------------------------------------------------
// Total counts (global, across all ranks and levels)
// -----------------------------------------------------------------------

void amrpd_total_np(PCP* pc, long long& np) { np = pc->TotalNumberOfParticles(); }

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
//   rdata[13]     td2
//   rdata[14]     td2a
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
        "mw", "dil", "damage", "nb0", "td2", "td2a"
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
