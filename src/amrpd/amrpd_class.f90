!> AMR-aware peridynamics solver
!>
!> Particles in NeighborParticleContainer<11,1>; bonds in ParticleContainer<4,5>.
!> Bonds are one-sided: each physical bond stored once, anchored at its lower-GID
!> endpoint's position so AMReX position-based Redistribute co-locates bond with
!> owner. Higher-GID endpoint accessed via particle ghost layer.
!> Force and dilatation contributions to the higher-GID (ghost) endpoint are
!> reduced back to owners via AMReX sumNeighbors (inverse fillNeighbors).
!>
!> Skeleton: compilable shell with stubbed physics kernels.
module amrpd_class
   use precision,     only: WP,I8
   use string,        only: str_medium
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use iso_c_binding
   implicit none
   private

   ! Public exports
   public :: amrpd,part,bond

   ! Particle motion-control bit flags (composed by bit-OR into idata[0])
   ! Standard cases:
   !   Free particle       = PART_MOVES + PART_INTEGRATES + PART_BONDS  (= 7)
   !   Clamped fixed       = PART_BONDS                                  (= 4)
   !   Velocity-prescribed = PART_MOVES + PART_BONDS                     (= 5)
   !   Witness/probe       = PART_MOVES                                  (= 1)
   !   Inactive (recycled) = PART_IS_DEAD                                (= 0)
   integer(c_int), parameter, public :: PART_IS_DEAD    = 0    !< Inactive; recycled out by AMReX
   integer(c_int), parameter, public :: PART_MOVES      = 1    !< pos += dt*vel during advance
   integer(c_int), parameter, public :: PART_INTEGRATES = 2    !< vel += (dt/2)*acc during Verlet half-kick
   integer(c_int), parameter, public :: PART_BONDS      = 4    !< Eligible for bond-network participation

   ! Domain boundary-condition flags (per face, set on lo_bc(d)/hi_bc(d))
   integer, parameter, public :: AMRPD_OPEN = 0    !< Particle leaves the domain (dropped by Redistribute)
   integer, parameter, public :: AMRPD_WALL = 1    !< Particle reflects off the domain face (position+velocity)

   ! Struct layout constants -- MUST match #defines in amrpd_wrapper.cpp
   integer, parameter, public :: AMRPD_NREAL_PART = 13
   integer, parameter, public :: AMRPD_NINT_PART  = 1
   integer, parameter, public :: AMRPD_NREAL_BOND = 4
   integer, parameter, public :: AMRPD_NINT_BOND  = 5

   ! Component indices into particle rdata, for use with amrpd_sum_neighbors
   ! (0-indexed to match AMReX's C++ convention)
   integer(c_int), parameter, public :: AMRPD_RC_VEL    = 0   !< vel[0..2]
   integer(c_int), parameter, public :: AMRPD_RC_FBOND  = 3   !< F_bond[0..2]  -- reduced
   integer(c_int), parameter, public :: AMRPD_RC_FFLUID = 6   !< F_fluid[0..2]
   integer(c_int), parameter, public :: AMRPD_RC_MW     = 9   !< mw            -- reduced
   integer(c_int), parameter, public :: AMRPD_RC_DIL    = 10  !< dil           -- reduced
   integer(c_int), parameter, public :: AMRPD_RC_DAMAGE = 11  !< damage        -- reduced
   integer(c_int), parameter, public :: AMRPD_RC_NB0    = 12  !< nb0           -- reduced (reference bond count)

   !> Solid particle struct -- must match C++ Particle<13,1> memory layout:
   !> pos[3], rdata[13], idcpu, idata[1]
   type, bind(C), public :: part
      real(c_double) :: pos(3)                  !< AMReX-managed position
      real(c_double) :: vel(3)                  !< rdata[0..2]
      real(c_double) :: F_bond(3)               !< rdata[3..5]  -- reduced via sumNeighbors
      real(c_double) :: F_fluid(3)              !< rdata[6..8]
      real(c_double) :: mw                      !< rdata[9]     -- reduced via sumNeighbors
      real(c_double) :: dil                     !< rdata[10]    -- reduced via sumNeighbors
      real(c_double) :: damage                  !< rdata[11]    -- reduced via sumNeighbors (broken-bond fraction in [0,1])
      real(c_double) :: nb0                     !< rdata[12]    -- reduced via sumNeighbors (reference bond count, stamped once at bond_init)
      integer(c_int64_t), private :: idcpu      !< AMReX packed id+cpu
      integer(c_int) :: flag                    !< idata[0]: PART_ALIVE or PART_IS_DEAD
   end type part

   !> Bond struct -- must match C++ Particle<4,5> memory layout:
   !> pos[3], rdata[5], idcpu, idata[5]
   !> pos is the position of the LOWER-GID endpoint (ownership invariant).
   type, bind(C), public :: bond
      real(c_double) :: pos(3)                  !< AMReX-managed; = pos(lower-GID endpoint)
      real(c_double) :: d0                      !< rdata[0]: reference distance
      real(c_double) :: w                       !< rdata[1]: cached influence weight
      real(c_double) :: damage                  !< rdata[2]: scalar damage (0=intact, 1=broken)
      real(c_double) :: hist1                   !< rdata[3]: packed periodic image offset (n_x,n_y,n_z) of the higher endpoint
      real(c_double) :: e_v                     !< rdata[4]: inelastic (Maxwell) deviatoric bond stretch
      integer(c_int64_t), private :: idcpu      !< AMReX packed id+cpu of this bond
      integer(c_int) :: id_lo_lo                !< idata[0]: low  32 bits of lower-GID endpoint idcpu
      integer(c_int) :: id_lo_hi                !< idata[1]: high 32 bits
      integer(c_int) :: id_hi_lo                !< idata[2]: low  32 bits of higher-GID endpoint idcpu
      integer(c_int) :: id_hi_hi                !< idata[3]: high 32 bits
      integer(c_int) :: alive                   !< idata[4]: 0=broken, 1=intact
   end type bond


   !> C interface bindings to amrpd_wrapper.cpp
   interface

      ! Lifecycle
      subroutine amrpd_new_pcp(pc,amrcore) bind(c)
         import :: c_ptr
         type(c_ptr) :: pc
         type(c_ptr), value :: amrcore
      end subroutine
      subroutine amrpd_delete_pcp(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine
      subroutine amrpd_new_pcb(pc,amrcore) bind(c)
         import :: c_ptr
         type(c_ptr) :: pc
         type(c_ptr), value :: amrcore
      end subroutine
      subroutine amrpd_delete_pcb(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      ! Redistribute
      subroutine amrpd_redistribute_p(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine
      subroutine amrpd_redistribute_b(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine

      ! Particle ghost/neighbor exchange and reduction
      subroutine amrpd_fill_neighbors_radius(pc,r) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: r
      end subroutine
      subroutine amrpd_clear_neighbors(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine
      subroutine amrpd_update_neighbors(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine
      subroutine amrpd_build_neighbor_list(pc,rcrit) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: rcrit
      end subroutine
      subroutine amrpd_get_neighbor_list_mfi(pc,lev,mfi,offsets,list,np,ntot) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: offsets,list
         integer(c_int64_t) :: np,ntot
      end subroutine
      subroutine amrpd_sum_neighbors(pc,rs,rn,is,in) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: rs,rn,is,in
      end subroutine

      ! MFIter accessors -- particles
      subroutine amrpd_get_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine
      subroutine amrpd_get_all_particles_mfi(pc,lev,mfi,dp,np_total,np_valid) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np_total,np_valid
      end subroutine
      subroutine amrpd_get_ghosts_mfi(pc,lev,mfi,dp,ng) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: ng
      end subroutine

      ! MFIter accessor -- bonds
      subroutine amrpd_get_bonds_mfi(pc,lev,mfi,dp,nb) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: nb
      end subroutine

      ! Single-element insertion (initialization)
      subroutine amrpd_add_particle_i(pc,lev,grid,tile,p) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,p
         integer(c_int), value :: lev,grid,tile
      end subroutine
      subroutine amrpd_add_bond_i(pc,lev,grid,tile,b) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,b
         integer(c_int), value :: lev,grid,tile
      end subroutine

      ! Bulk append at level 0 (collective; ranks with n=0 pass raw=NULL)
      subroutine amrpd_append_particles(pc,raw,n) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         type(c_ptr), value :: raw
         integer(c_int64_t), value :: n
      end subroutine
      subroutine amrpd_append_bonds(pc,raw,n) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         type(c_ptr), value :: raw
         integer(c_int64_t), value :: n
      end subroutine

      ! BoxArray / DistributionMap accessors
      subroutine amrpd_get_particle_boxarray_p(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: ba
      end subroutine
      subroutine amrpd_get_particle_distromap_p(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: dm
      end subroutine
      subroutine amrpd_set_particle_boxarray_p(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: ba
      end subroutine
      subroutine amrpd_set_particle_distromap_p(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: dm
      end subroutine
      subroutine amrpd_get_particle_boxarray_b(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: ba
      end subroutine
      subroutine amrpd_get_particle_distromap_b(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: dm
      end subroutine
      subroutine amrpd_set_particle_boxarray_b(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: ba
      end subroutine
      subroutine amrpd_set_particle_distromap_b(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: dm
      end subroutine

      ! ID/CPU counters and accessors
      subroutine amrpd_get_next_id_p(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine
      subroutine amrpd_set_next_id_p(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t), value :: id
      end subroutine
      subroutine amrpd_get_next_id_b(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine
      subroutine amrpd_set_next_id_b(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t), value :: id
      end subroutine
      subroutine amrpd_get_cpu(cpu) bind(c)
         import :: c_int
         integer(c_int) :: cpu
      end subroutine
      subroutine amrpd_get_particle_id(id,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t) :: id
         type(c_ptr), value :: p
      end subroutine
      subroutine amrpd_set_particle_id(id,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t), value :: id
         type(c_ptr), value :: p
      end subroutine
      subroutine amrpd_get_particle_cpu(cpu,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int) :: cpu
         type(c_ptr), value :: p
      end subroutine
      subroutine amrpd_set_particle_cpu(cpu,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int), value :: cpu
         type(c_ptr), value :: p
      end subroutine
      subroutine amrpd_get_bond_id(id,b) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t) :: id
         type(c_ptr), value :: b
      end subroutine
      subroutine amrpd_set_bond_id(id,b) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t), value :: id
         type(c_ptr), value :: b
      end subroutine
      subroutine amrpd_get_bond_cpu(cpu,b) bind(c)
         import :: c_int,c_ptr
         integer(c_int) :: cpu
         type(c_ptr), value :: b
      end subroutine
      subroutine amrpd_set_bond_cpu(cpu,b) bind(c)
         import :: c_int,c_ptr
         integer(c_int), value :: cpu
         type(c_ptr), value :: b
      end subroutine

      ! Unique 64-bit GID key composed from (id, cpu) of a particle, used to
      ! build the gid_hash for fast bond-endpoint resolution.
      subroutine amrpd_get_particle_idcpu(key,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t) :: key
         type(c_ptr), value :: p
      end subroutine

      ! Global counts
      subroutine amrpd_total_np(pc,np) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: np
      end subroutine
      subroutine amrpd_total_nb(pc,nb) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: nb
      end subroutine

      ! Checkpoint I/O: caller composes fullpath (e.g. <dirname>/particles)
      subroutine amrpd_checkpoint_p(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char) :: path(*)
      end subroutine
      subroutine amrpd_restart_p(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char) :: path(*)
      end subroutine
      subroutine amrpd_checkpoint_b(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char) :: path(*)
      end subroutine
      subroutine amrpd_restart_b(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char) :: path(*)
      end subroutine

   end interface


   !> Peridynamics solver type
   type :: amrpd

      !> Associated AMR grid
      class(amrgrid), pointer :: amr => null()

      !> Opaque AMReX container handles
      type(c_ptr) :: pcp = c_null_ptr           !< Particle container (NeighborParticleContainer<11,1>)
      type(c_ptr) :: pcb = c_null_ptr           !< Bond     container (ParticleContainer<4,5>)

      !> Solver name
      character(len=str_medium) :: name = 'UNNAMED_AMRPD'

      !> Global counts
      integer(I8) :: np = 0                     !< Global particle count
      integer(I8) :: nb = 0                     !< Global bond count (total, alive + broken)
      integer(I8) :: nb_broken = 0              !< Global count of broken (alive=0) bonds

      !> Local count + load-balance metrics across ranks
      integer(I8) :: np_loc = 0                 !< This rank's particle count
      integer(I8) :: np_min = 0                 !< Min particle count across ranks
      integer(I8) :: np_max = 0                 !< Max particle count across ranks
      real(WP)    :: np_eff = 0.0_WP            !< Load efficiency = mean/max across ranks

      integer(I8) :: nb_loc = 0                 !< This rank's bond count
      integer(I8) :: nb_min = 0                 !< Min bond count across ranks
      integer(I8) :: nb_max = 0                 !< Max bond count across ranks
      real(WP)    :: nb_eff = 0.0_WP            !< Load efficiency for bonds

      !> Material/physical parameters
      real(WP) :: elastic_modulus = 0.0_WP      !< Young's modulus
      real(WP) :: poisson_ratio   = 0.0_WP      !< Poisson's ratio
      real(WP) :: rho             = 0.0_WP      !< Material density
      real(WP) :: crit_energy     = 0.0_WP      !< Critical energy release rate G_c
      real(WP) :: s0              = huge(1.0_WP)!< Critical bond stretch (set by bond_init from G_c if >0; huge() = no damage)
      real(WP) :: tau             = huge(1.0_WP)!< Maxwell deviatoric relaxation time (huge = purely elastic, no viscoplastic flow)
      real(WP) :: visc_lambda     = 1.0_WP      !< SLS relaxing fraction [0,1] (1 = pure Maxwell/full flow; <1 keeps long-term elastic stiffness)
      real(WP) :: fail_stretch    = huge(1.0_WP)!< Direct failure-stretch override (huge = use G_c-derived s0; finite = ductile, decoupled from G_c)
      real(WP) :: dV              = 0.0_WP      !< Element (representative) volume

      !> Short-range contact (soft-sphere model ported from amrlpt%collide).
      !> Contact duration tau_col defaults to 5*dt (as stiff as integrable) but
      !> can be user-overridden by setting tau_col > 0. CFLc = dt/tau_col is
      !> reported as a diagnostic but does NOT constrain dt.
      real(WP) :: contact_dist    = 0.0_WP      !< d_c: contact threshold (default 0.9 * dV^(1/3))
      real(WP) :: tau_col         = 0.0_WP      !< Collision duration (<=0 -> auto = 5*dt each step)
      real(WP) :: e_n             = 0.7_WP      !< Normal restitution coefficient (particle-particle)
      real(WP) :: e_w             = 0.7_WP      !< Normal restitution coefficient (wall / IB)
      real(WP) :: clip_col        = 0.2_WP      !< Overlap clip fraction of d_eff

      !> Bonding parameters
      real(WP) :: delta           = 0.0_WP      !< Reference horizon (bonding distance)
      !> Ghost-layer search radius used by every fill_ghosts call (bond_init,
      !> compute_dilatation, compute_force, ...). Must be >= delta. Should also
      !> be >= the largest bond length any bond can reach before it breaks, so
      !> stretched bonds always have both endpoints accessible via the ghost
      !> layer. Default placeholder: 1.5 * delta. Will eventually be
      !> (1 + max_stretch) * delta + safety once damage is wired (M5).
      real(WP) :: search_radius = 0.0_WP

      !> Gravitational acceleration
      real(WP), dimension(3) :: gravity = 0.0_WP

      !> Maximum AMR level particles are allowed on (cap passed to AMReX
      !> Redistribute as lev_max). Particles span levels [0, maxlvl] and
      !> AMReX places each at the finest level covering its position.
      !> Defaults to amr%maxlvl in initialize.
      integer :: maxlvl = 0

      !> Overlap (ghost cell) width. Must be a multiple of the AMR refinement
      !> ratio (2 in standard AMReX setups) because amrdata's process_deposit
      !> uses sum_fine_to_coarse, which asserts nGrow % ratio == 0. Matches
      !> amrlpt's default.
      integer :: nover = 2

      !> Per-face domain BCs (default: open on all faces).
      !> Override with AMRPD_WALL for hard reflection. Periodic faces (set on
      !> amrgrid via xper/yper/zper) skip this logic — AMReX wraps automatically
      !> during Redistribute.
      integer, dimension(3) :: lo_bc = AMRPD_OPEN
      integer, dimension(3) :: hi_bc = AMRPD_OPEN

      !> Knapsack rebalancing of the particle (and bond) containers.
      !> When .true., post_regrid builds a knapsack DM weighted by per-box
      !> particle count and applies it to BOTH containers, decoupling solid-
      !> particle load balance from the fluid grid's DM. Solid bodies typically
      !> concentrate in a small region, so this matters whenever np_max/np_mean
      !> across ranks gets large under the fluid-driven DM.
      !> Bonds inherit the same DM so the lower-GID ownership invariant survives
      !> redistribute (bonds stay co-located with their owner particles).
      logical :: rebalance = .false.

      !> Particle volume fraction on the Eulerian AMR mesh. Cell-centered scalar
      !> (one component), one ghost layer. Updated each advance step by
      !> update_VF: trilinear deposition of each particle's volume dV onto the
      !> 8 surrounding cell centers, then average-down + optional smoothing.
      !> Drives the AMR tagging callback when VF_tag > 0.
      type(amrdata) :: VF
      real(WP)     :: VF_tag       = -1.0_WP    !< Refinement threshold (<=0 disables VF-driven tagging)
      real(WP)     :: filter_width =  0.0_WP    !< Gaussian-equivalent filter width for VF; 0 disables
      real(WP)     :: VFmin=0.0_WP,VFmax=0.0_WP,VFmean=0.0_WP   !< VF statistics

      !> Optional user-supplied tagging callback. Called AFTER the built-in VF
      !> tagging. Use to add custom refinement criteria (e.g., damage > 0.3).
      procedure(pd_tagging_iface), pointer, pass :: user_pd_tagging => null()

      !> Monitoring info
      real(WP) :: Umin=0.0_WP,Umax=0.0_WP,Umean=0.0_WP
      real(WP) :: Vmin=0.0_WP,Vmax=0.0_WP,Vmean=0.0_WP
      real(WP) :: Wmin=0.0_WP,Wmax=0.0_WP,Wmean=0.0_WP
      real(WP) :: CFLp=0.0_WP                                  !< convective: max(|v_d|) * dt / dp -- binds, limit 0.1 (scaled by 5)
      real(WP) :: CFLe=0.0_WP                                  !< elastic-wave: c_p * dt / dp -- binds, limit 0.5
      real(WP) :: CFLc=0.0_WP                                  !< contact: dt / tau_col -- diagnostic only, does not bind dt
      real(WP) :: CFLv=0.0_WP                                  !< viscous: dt / tau -- diagnostic only (exponential relaxation is unconditionally stable)

   contains
      ! Lifecycle
      procedure :: initialize
      procedure :: finalize
      ! Container utilities
      procedure :: redistribute
      procedure :: fill_ghosts
      procedure :: clear_ghosts
      procedure :: update_ghosts
      procedure :: build_neighbor_list
      procedure, private :: get_neighbor_list
      procedure :: sum_ghosts_force
      procedure :: sum_ghosts_mw
      procedure :: sum_ghosts_dil
      procedure :: sum_ghosts_damage
      procedure :: sum_ghosts_nb0
      procedure :: get_info             !< Global counts + min/max/mean velocities + load-balance metrics
      procedure :: log_box_loads        !< Per-tile (rank, lev, grid, tile, np, nb) dump to stdout
      procedure :: set_particle_ba_p
      procedure :: set_particle_dm_p
      procedure :: set_particle_ba_b
      procedure :: set_particle_dm_b
      ! Particle population
      procedure :: append
      procedure :: append_bonds
      ! MFIter helpers (particle container's BA/DM)
      procedure :: mfiter_build
      procedure :: mfiter_destroy
      procedure :: get_particles
      procedure :: get_all_particles
      procedure :: get_ghosts
      procedure :: get_bonds
      ! AMR callbacks
      procedure :: post_regrid
      procedure :: tagging
      ! Particle volume fraction + AMR tagging
      procedure :: update_VF                         !< Compute VF from particle positions (trilinear deposit)
      procedure :: process_deposit                   !< Post-process a deposited field (extensive -> intensive + C/F transfers; public: also used on driver-deposited fields)
      procedure :: filter                            !< Explicit-diffusion smoothing of a cell-centered amrdata (public: also used on driver-deposited fields)
      ! Physics -- STUBBED in skeleton
      procedure :: bond_init
      procedure :: compute_dilatation
      procedure :: compute_force
      procedure :: compute_contact
      procedure :: interp                !< Trilinear cell-centered interpolation (used by compute_contact for IB)
      procedure :: advance
      procedure :: get_cfl
      ! Checkpoint I/O
      procedure :: write
      procedure :: read
      ! Diagnostics
      procedure :: print
   end type amrpd


   !> Abstract interface for user-overridable tagging callback. Invoked AFTER
   !> the built-in VF-based tagging by the registered AMReX tagging dispatch.
   abstract interface
      subroutine pd_tagging_iface(solver,lvl,time,tags)
         import :: amrpd,c_ptr,WP
         class(amrpd), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine pd_tagging_iface
   end interface


contains


   ! ============================================================================
   ! DISPATCHERS (module-level) -- recover concrete amrpd type from c_ptr ctx
   ! ============================================================================

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrpd_postregrid_dispatch(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrpd), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrpd_postregrid_dispatch

   !> Dispatch tagging: calls type-bound method, then user override (if any)
   subroutine amrpd_tagging_dispatch(ctx,lvl,time,tags)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrpd), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_pd_tagging)) call this%user_pd_tagging(lvl,time,tags)
   end subroutine amrpd_tagging_dispatch


   ! ============================================================================
   ! LIFECYCLE
   ! ============================================================================

   !> Initialize amrpd solver: create particle and bond containers, register AMR callbacks
   subroutine initialize(this,amr,name)
      use amrex_amr_module, only: amrex_bc_foextrap
      implicit none
      class(amrpd), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), optional :: name
      ! Set solver name
      if (present(name)) this%name = trim(adjustl(name))
      ! Point to associated AMR grid
      this%amr => amr
      ! Default level cap: allow particles up to the AMR grid's max refinement
      this%maxlvl = amr%maxlvl
      ! Create AMReX particle and bond containers
      call amrpd_new_pcp(this%pcp,this%amr%amrcore)
      call amrpd_new_pcb(this%pcb,this%amr%amrcore)
      ! Particle volume fraction field (cell-centered, 1 ghost layer; foextrap
      ! on non-periodic faces matches amrlpt's convention)
      call this%VF%initialize(amr=amr,name='VF',ncomp=1,ng=this%nover); call this%VF%register()
      if (.not.this%amr%xper) then; this%VF%lo_bc(1,1)=amrex_bc_foextrap; this%VF%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.this%amr%yper) then; this%VF%lo_bc(2,1)=amrex_bc_foextrap; this%VF%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.this%amr%zper) then; this%VF%lo_bc(3,1)=amrex_bc_foextrap; this%VF%hi_bc(3,1)=amrex_bc_foextrap; end if
      ! Register AMR callbacks (post_regrid + tagging) so containers stay in sync
      select type (this)
       type is (amrpd)
         call this%amr%add_postregrid(amrpd_postregrid_dispatch,c_loc(this))
         call this%amr%add_tagging   (amrpd_tagging_dispatch,   c_loc(this))
      end select
      ! Print solver info
      call this%print()
   end subroutine initialize

   !> Finalize: destroy containers and release amr pointer
   subroutine finalize(this)
      implicit none
      class(amrpd), intent(inout) :: this
      ! Drop user tagging hook
      nullify(this%user_pd_tagging)
      ! Tear down the volume-fraction field
      call this%VF%finalize()
      if (c_associated(this%pcb)) then
         call amrpd_delete_pcb(this%pcb); this%pcb = c_null_ptr
      end if
      if (c_associated(this%pcp)) then
         call amrpd_delete_pcp(this%pcp); this%pcp = c_null_ptr
      end if
      nullify(this%amr)
   end subroutine finalize


   ! ============================================================================
   ! CONTAINER UTILITIES
   ! ============================================================================

   !> Redistribute both particle and bond containers across ranks. Particles
   !> can land on any currently-defined level (AMReX picks the finest covering
   !> level for each particle's position). lev_max=-1 tells AMReX to use the
   !> AmrCore's current finestLevel(), which avoids the
   !>   `Assertion lev_max <= finestLevel()' failed
   !> crash when redistribute is called before all maxlvl levels exist (e.g.,
   !> right after init_from_scratch). Matches amrlpt's default.
   subroutine redistribute(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_redistribute_p(this%pcp,0,-1,0)
      call amrpd_redistribute_b(this%pcb,0,-1,0)
   end subroutine redistribute

   !> Fill particle ghost layer to a given physical radius (in lieu of cell count)
   subroutine fill_ghosts(this,radius)
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP), intent(in) :: radius
      call amrpd_fill_neighbors_radius(this%pcp,real(radius,c_double))
   end subroutine fill_ghosts

   !> Clear particle ghost layer
   subroutine clear_ghosts(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_clear_neighbors(this%pcp)
   end subroutine clear_ghosts

   !> Refresh existing ghost data with current owner values (topology unchanged).
   !> Use after a kernel updates owner state (e.g., compute_dilatation -> theta)
   !> so subsequent kernels reading ghost rdata see the new values.
   subroutine update_ghosts(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_update_neighbors(this%pcp)
   end subroutine update_ghosts

   !> Build a per-tile CSR neighbor list of particle pairs with |r_i - r_j| < rcrit.
   !> Must be called AFTER fill_ghosts. Cell-binned by AMReX using the existing
   !> ghost mask; rcrit is the predicate radius. Consumed by get_neighbor_list.
   subroutine build_neighbor_list(this,rcrit)
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP), intent(in) :: rcrit
      call amrpd_build_neighbor_list(this%pcp,real(rcrit,c_double))
   end subroutine build_neighbor_list

   !> Get the CSR neighbor-list pointers (offsets, list) for the current tile.
   !> Iteration:  for i=1..np, do k=off(i)+1,off(i+1); j=lst(k)+1; ...
   !> Indices in lst are 0-based and index the COMBINED real+ghost particle
   !> array (same range as get_all_particles).
   subroutine get_neighbor_list(this,lvl,mfi,off,lst)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      integer(c_int32_t), dimension(:), pointer, intent(out) :: off,lst
      type(c_ptr) :: dp_off,dp_lst
      integer(c_int64_t) :: np_c,ntot_c
      call amrpd_get_neighbor_list_mfi(this%pcp,lvl,mfi%p,dp_off,dp_lst,np_c,ntot_c)
      if (np_c.gt.0_c_int64_t) then
         call c_f_pointer(dp_off,off,[np_c+1_c_int64_t])
      else
         nullify(off)
      end if
      if (ntot_c.gt.0_c_int64_t) then
         call c_f_pointer(dp_lst,lst,[ntot_c])
      else
         nullify(lst)
      end if
   end subroutine get_neighbor_list

   !> Reduce F_bond contributions from ghost slots back to owners (inverse fillNeighbors).
   !> F_bond occupies rdata[3..5]; no idata reduction.
   subroutine sum_ghosts_force(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_sum_neighbors(this%pcp,AMRPD_RC_FBOND,3,0,0)
   end subroutine sum_ghosts_force

   !> Reduce mw contributions from ghost slots back to owners. mw is reference-
   !> only (stamped once in bond_init); keep it on its own reducer so subsequent
   !> ghost-buffer activity for other fields doesn't re-touch it.
   subroutine sum_ghosts_mw(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_sum_neighbors(this%pcp,AMRPD_RC_MW,1,0,0)
   end subroutine sum_ghosts_mw

   !> Reduce dil contributions from ghost slots back to owners. dil is
   !> recomputed every step in compute_dilatation.
   subroutine sum_ghosts_dil(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_sum_neighbors(this%pcp,AMRPD_RC_DIL,1,0,0)
   end subroutine sum_ghosts_dil

   !> Reduce damage-fraction contributions from ghost slots back to owners.
   subroutine sum_ghosts_damage(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_sum_neighbors(this%pcp,AMRPD_RC_DAMAGE,1,0,0)
   end subroutine sum_ghosts_damage

   !> Reduce nb0 (reference bond count) contributions from ghost slots back to
   !> owners. Reference-only -- called once at bond_init alongside sum_ghosts_mw.
   subroutine sum_ghosts_nb0(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_sum_neighbors(this%pcp,AMRPD_RC_NB0,1,0,0)
   end subroutine sum_ghosts_nb0

   !> Compute global counts, per-rank load metrics, and min/max/mean velocity
   !> statistics in one collective pass. Mirrors amrlpt's get_info pattern: walks
   !> owned particles once via MFIter accumulating local stats, then does the
   !> global reductions (counts via SUM, min/max via MIN/MAX, etc).
   !>
   !> Populates on this:
   !>   np / np_loc / np_min / np_max / np_eff    (and nb analogues for bonds)
   !>   Umin/Umax/Umean, Vmin/Vmax/Vmean, Wmin/Wmax/Wmean
   subroutine get_info(this)
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP) :: my_Usum,my_Vsum,my_Wsum,safe_np
      integer(c_int64_t) :: nbtot

      ! Init per-rank accumulators
      this%np_loc=0_I8
      this%Umin= huge(1.0_WP); this%Umax=-huge(1.0_WP); my_Usum=0.0_WP
      this%Vmin= huge(1.0_WP); this%Vmax=-huge(1.0_WP); my_Vsum=0.0_WP
      this%Wmin= huge(1.0_WP); this%Wmax=-huge(1.0_WP); my_Wsum=0.0_WP

      ! Per-rank loop: count and accumulate velocity stats over all levels
      local_pass: block
         use amrex_amr_module, only: amrex_mfiter
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         integer(I8) :: np_,n
         integer :: lvl
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  this%np_loc=this%np_loc+1_I8
                  this%Umin=min(this%Umin,p(n)%vel(1)); this%Umax=max(this%Umax,p(n)%vel(1)); my_Usum=my_Usum+p(n)%vel(1)
                  this%Vmin=min(this%Vmin,p(n)%vel(2)); this%Vmax=max(this%Vmax,p(n)%vel(2)); my_Vsum=my_Vsum+p(n)%vel(2)
                  this%Wmin=min(this%Wmin,p(n)%vel(3)); this%Wmax=max(this%Wmax,p(n)%vel(3)); my_Wsum=my_Wsum+p(n)%vel(3)
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do
      end block local_pass

      ! Global reductions
      global_reduce: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8
         use parallel, only: MPI_REAL_WP
         integer :: ierr
         ! Particle counts and load balance: seed min/max with this rank's np_loc,
         ! then reduce; the global np is the SUM of np_loc across ranks.
         this%np_min=this%np_loc; this%np_max=this%np_loc; this%np=this%np_loc; this%np_eff=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_min,1,MPI_INTEGER8,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_max,1,MPI_INTEGER8,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np,    1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
         if (this%np_max.gt.0_I8) this%np_eff=real(this%np,WP)/real(this%np_max,WP)/real(this%amr%nproc,WP)
         ! Velocity min/max/mean
         safe_np=real(max(this%np,1_I8),WP)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,my_Usum  ,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr); this%Umean=my_Usum/safe_np
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,my_Vsum  ,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr); this%Vmean=my_Vsum/safe_np
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,my_Wsum  ,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr); this%Wmean=my_Wsum/safe_np
      end block global_reduce

      ! Bond counts: global total via AMReX TotalNumberOfParticles. Per-rank
      ! pass counts local bonds and broken bonds in one sweep; global min/max/eff
      ! and broken total via reduce.
      call amrpd_total_nb(this%pcb,nbtot); this%nb=int(nbtot,I8)
      count_bonds: block
         use amrex_amr_module, only: amrex_mfiter
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8
         type(amrex_mfiter) :: mfi
         type(bond), dimension(:), pointer :: b
         integer(I8) :: nb_tile,ib
         integer :: lvl,ierr
         this%nb_loc=0_I8; this%nb_broken=0_I8
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_bonds(lvl,mfi,b,nb_tile)
               this%nb_loc=this%nb_loc+nb_tile
               do ib=1_I8,nb_tile
                  if (b(ib)%alive.eq.0) this%nb_broken=this%nb_broken+1_I8
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do
         this%nb_min=this%nb_loc; this%nb_max=this%nb_loc; this%nb_eff=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nb_min,1,MPI_INTEGER8,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nb_max,1,MPI_INTEGER8,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nb_broken,1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
         if (this%nb_max.gt.0_I8) this%nb_eff=real(this%nb,WP)/real(this%nb_max,WP)/real(this%amr%nproc,WP)
      end block count_bonds
   end subroutine get_info

   !> Per-tile (rank, lev, grid_index, tile_index, np, nb) dump to stdout.
   !> Output interleaves across ranks; grep "[BOXLOAD]" and sort to view.
   !> Use as a one-shot diagnostic after a regrid to see balance directly.
   subroutine log_box_loads(this)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      type(bond), dimension(:), pointer :: b
      integer(I8) :: np_,nb_
      integer :: lvl,gi,ti
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            gi=mfi%grid_index(); ti=mfi%local_tile_index()
            call this%get_particles(lvl,mfi,p,np_)
            call this%get_bonds(lvl,mfi,b,nb_)
            write(*,'("[BOXLOAD] rank=",i0," lev=",i0," grid=",i0," tile=",i0," np=",i0," nb=",i0)') &
            &  this%amr%rank,lvl,gi,ti,np_,nb_
         end do
         call this%mfiter_destroy(mfi)
      end do
   end subroutine log_box_loads

   !> Set particle-container BoxArray for a given level
   subroutine set_particle_ba_p(this,lvl,ba)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      call amrpd_set_particle_boxarray_p(this%pcp,lvl,ba%p)
   end subroutine set_particle_ba_p

   !> Set particle-container DistributionMapping for a given level
   subroutine set_particle_dm_p(this,lvl,dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_distromap), intent(in) :: dm
      call amrpd_set_particle_distromap_p(this%pcp,lvl,dm%p)
   end subroutine set_particle_dm_p

   !> Set bond-container BoxArray for a given level
   subroutine set_particle_ba_b(this,lvl,ba)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      call amrpd_set_particle_boxarray_b(this%pcb,lvl,ba%p)
   end subroutine set_particle_ba_b

   !> Set bond-container DistributionMapping for a given level
   subroutine set_particle_dm_b(this,lvl,dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_distromap), intent(in) :: dm
      call amrpd_set_particle_distromap_b(this%pcb,lvl,dm%p)
   end subroutine set_particle_dm_b

   !> Build an MFIter over the particle container's BA/DM at level lvl.
   subroutine mfiter_build(this,lvl,mfi,tiling)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_mfiter_build
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(out) :: mfi
      logical, intent(in), optional :: tiling
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      logical :: use_tiling
      use_tiling=.false.; if (present(tiling)) use_tiling=tiling
      call amrpd_get_particle_boxarray_p (this%pcp,lvl,ba%p)
      call amrpd_get_particle_distromap_p(this%pcp,lvl,dm%p)
      call amrex_mfiter_build(mfi,ba,dm,tiling=use_tiling)
   end subroutine mfiter_build

   !> Destroy an MFIter built via mfiter_build.
   subroutine mfiter_destroy(this,mfi)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_destroy
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrex_mfiter), intent(inout) :: mfi
      call amrex_mfiter_destroy(mfi)
   end subroutine mfiter_destroy

   !> Return a Fortran pointer to the valid particle array on the current tile
   !> (ghost particles excluded). np is the number of valid particles.
   subroutine get_particles(this,lvl,mfi,p,np)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_c
      call amrpd_get_particles_mfi(this%pcp,lvl,mfi%p,dp,np_c)
      np=int(np_c,I8)
      if (np.gt.0_I8) then
         call c_f_pointer(dp,p,[np])
      else
         nullify(p)
      end if
   end subroutine get_particles

   !> Return a Fortran pointer to the COMBINED (valid + ghost) particle array on
   !> the current tile. np_total = valid + ghost; np_valid = valid count only.
   !> Ghosts are at indices [np_valid+1 .. np_total]. Use after fill_ghosts.
   subroutine get_all_particles(this,lvl,mfi,p,np_total,np_valid)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np_total,np_valid
      type(c_ptr) :: dp
      integer(c_int64_t) :: nt_c,nv_c
      call amrpd_get_all_particles_mfi(this%pcp,lvl,mfi%p,dp,nt_c,nv_c)
      np_total=int(nt_c,I8); np_valid=int(nv_c,I8)
      if (np_total.gt.0_I8) then
         call c_f_pointer(dp,p,[np_total])
      else
         nullify(p)
      end if
   end subroutine get_all_particles

   !> Return a Fortran pointer to the writable ghost-particle slice for the
   !> current tile. The pointer targets AMReX's internal neighbor buffer (the
   !> source that sum_ghosts_mw / sum_ghosts_dil / sum_ghosts_force read from when reducing
   !> back to owners) — NOT the ghost copies appended to the ptile by
   !> fill_ghosts. The two hold identical data right after fill_ghosts, but
   !> only writes to this buffer are picked up by the reduction.
   !> ng matches (np_total - np_valid) from get_all_particles.
   subroutine get_ghosts(this,lvl,mfi,pg,ng)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: pg
      integer(I8), intent(out) :: ng
      type(c_ptr) :: dp
      integer(c_int64_t) :: ng_c
      call amrpd_get_ghosts_mfi(this%pcp,lvl,mfi%p,dp,ng_c)
      ng=int(ng_c,I8)
      if (ng.gt.0_I8) then
         call c_f_pointer(dp,pg,[ng])
      else
         nullify(pg)
      end if
   end subroutine get_ghosts

   !> Return a Fortran pointer to the bond array on the current tile.
   subroutine get_bonds(this,lvl,mfi,b,nb)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(bond), dimension(:), pointer, intent(out) :: b
      integer(I8), intent(out) :: nb
      type(c_ptr) :: dp
      integer(c_int64_t) :: nb_c
      call amrpd_get_bonds_mfi(this%pcb,lvl,mfi%p,dp,nb_c)
      nb=int(nb_c,I8)
      if (nb.gt.0_I8) then
         call c_f_pointer(dp,b,[nb])
      else
         nullify(b)
      end if
   end subroutine get_bonds

   !> Bulk-append Fortran particle array into the container at level 0.
   !> Collective: every rank must call; ranks with nothing to add pass n=0.
   !> AMReX assigns unique (id,cpu) to each appended particle and
   !> AddParticlesAtLevel internally redistributes by position.
   subroutine append(this,plist,n)
      use messager, only: die
      implicit none
      class(amrpd), intent(inout) :: this
      type(part), dimension(:), allocatable, target, intent(in) :: plist
      integer(I8), intent(in) :: n
      type(c_ptr) :: raw
      raw=c_null_ptr
      if (n.gt.0_I8.and.allocated(plist)) then
         if (int(size(plist),I8).lt.n) call die('[amrpd append] plist array smaller than n')
         raw=c_loc(plist(1))
      end if
      call amrpd_append_particles(this%pcp,raw,int(n,c_int64_t))
   end subroutine append

   !> Bulk-append Fortran bond array into the bond container at level 0.
   !> Collective: every rank must call; ranks with nothing to add pass n=0.
   !> AMReX assigns each bond a unique (id,cpu) and AddParticlesAtLevel
   !> internally redistributes by spatial position (= position of the lower-GID
   !> endpoint, per the bond ownership invariant).
   subroutine append_bonds(this,blist,n)
      use messager, only: die
      implicit none
      class(amrpd), intent(inout) :: this
      type(bond), dimension(:), allocatable, target, intent(in) :: blist
      integer(I8), intent(in) :: n
      type(c_ptr) :: raw
      raw=c_null_ptr
      if (n.gt.0_I8.and.allocated(blist)) then
         if (int(size(blist),I8).lt.n) call die('[amrpd append_bonds] blist array smaller than n')
         raw=c_loc(blist(1))
      end if
      call amrpd_append_bonds(this%pcb,raw,int(n,c_int64_t))
   end subroutine append_bonds


   ! ============================================================================
   ! AMR CALLBACKS
   ! ============================================================================

   !> Post-regrid: re-sync particle and bond containers' BA/DM, redistribute.
   !> When rebalance is enabled, additionally build a knapsack DM weighted by
   !> per-box BOND counts (bonds dominate compute_dilatation/compute_force cost,
   !> so balancing bond load is what actually matters). Applied to BOTH containers
   !> so bonds stay co-located with their lower-GID owners. Particles span levels
   !> [0, maxlvl] and the BA/DM sync walks all current levels so the containers
   !> track AMR.
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      integer :: lvl

      ! Re-sync both containers to the fluid grid's BA/DM. This sets the
      ! AmrParGDB's per-level particle BA/DM caches. Keeping these in sync
      ! enables the dual-DM load-balancing path (knapsack rebalance below)
      ! and matches the pattern used in amrlpt.
      sync_to_amrgrid: block
         do lvl=0,this%amr%clvl()
            call this%set_particle_ba_p(lvl,this%amr%get_boxarray(lvl))
            call this%set_particle_dm_p(lvl,this%amr%get_distromap(lvl))
            call this%set_particle_ba_b(lvl,this%amr%get_boxarray(lvl))
            call this%set_particle_dm_b(lvl,this%amr%get_distromap(lvl))
         end do
      end block sync_to_amrgrid

      ! First settle particles onto the freshly-synced BA/DM before we can
      ! count them per box for the knapsack weighting below
      call this%redistribute()

      ! Optionally retarget BOTH containers' DM via knapsack on particle count.
      ! Same DM applied to particles and bonds → lower-GID ownership preserved.
      if (this%rebalance) then
         knapsack_rebalance: block
            use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_distromap_destroy,amrex_mfiter
            use amrex_interface,  only: amrdm_make_knapsack
            use parallel,         only: MPI_REAL_WP
            use mpi_f08,          only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            integer :: lvl_,nboxes,ierr
            real(WP), dimension(:), allocatable :: costs
            type(amrex_boxarray)  :: ba
            type(amrex_distromap) :: new_dm
            type(amrex_mfiter)    :: mfi
            type(bond), dimension(:), pointer :: b
            integer(I8) :: nb_
            do lvl_=0,this%amr%clvl()
               ! Per-box bond count on this rank (cost metric for knapsack)
               ba=this%amr%get_boxarray(lvl_); nboxes=int(ba%nboxes())
               allocate(costs(nboxes)); costs=0.0_WP
               call this%mfiter_build(lvl_,mfi)
               do while (mfi%next())
                  call this%get_bonds(lvl_,mfi,b,nb_)
                  costs(mfi%grid_index()+1)=real(nb_,WP)
               end do
               call this%mfiter_destroy(mfi)
               ! Global sum so every rank sees the full cost vector
               call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
               ! Knapsack DM weighted by bond count, applied to both containers
               call amrdm_make_knapsack(new_dm,costs,nboxes)
               call this%set_particle_dm_p(lvl_,new_dm)
               call this%set_particle_dm_b(lvl_,new_dm)
               call amrex_distromap_destroy(new_dm)
               deallocate(costs)
            end do
            ! Redistribute onto the new (rebalanced) DM
            call this%redistribute()
         end block knapsack_rebalance
      end if

      ! Recompute particle VF on the (new) mesh. AMReX fills new fine cells via
      ! amrdata's coarse-to-fine interpolation during regrid, but that's a
      ! guess; a fresh deposit from the redistributed particles is the truth.
      ! Matches amrlpt%post_regrid.
      call this%update_VF()
   end subroutine post_regrid

   !> Tag cells for refinement where particle VF exceeds VF_tag. Disabled if
   !> VF_tag <= 0. Mirrors amrlpt%tagging.
   subroutine tagging(this,lvl,time,tags)
      use amrex_amr_module, only: amrex_tagboxarray,amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use amrgrid_class,    only: SETtag
      implicit none
      class(amrpd), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrex_tagboxarray) :: tba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer :: i,j,k
      ! Skip if VF tagging is disabled
      if (this%VF_tag.le.0.0_WP) return
      ! Resolve tagboxarray pointer
      tba=tags
      ! Loop over tiles and tag
      call amrex_mfiter_build(mfi,this%VF%mf(lvl))
      do while (mfi%next())
         tagarr=>tba%dataPtr(mfi)
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            if (pVF(i,j,k,1).gt.this%VF_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine tagging

   !> Compute particle volume fraction on the Eulerian AMR mesh: trilinear
   !> deposit of each particle's volume dV onto the 8 surrounding cell centers,
   !> then convert extensive->intensive (divide by cell_vol), propagate across
   !> C/F boundaries, average down, fill ghosts, and optionally smooth.
   subroutine update_VF(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_multifab,amrex_multifab_build,amrex_multifab_destroy
      use amrex_distromap_module, only: operator(.eq.)
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer(I8) :: np_,i
      integer :: lvl,ii,jj,kk
      real(WP) :: dxi,dyi,dzi,wx,wy,wz,Vp
      type(amrex_multifab) :: tmpVF
      logical :: dual_dm

      ! Zero VF on all levels
      call this%VF%setval(0.0_WP)
      Vp=this%dV

      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)

         ! Particles may live on a different DM than VF -- deposit on the
         ! particle DM into a scratch mfab, then parallel_copy to VF
         dual_dm=(.not.(this%VF%mf(lvl)%dm.eq.get_pdm()))
         if (dual_dm) then
            call amrex_multifab_build(mf=tmpVF,ba=this%amr%ba(lvl),dm=get_pdm(),nc=this%VF%mf(lvl)%ncomp(),ng=this%VF%mf(lvl)%nghost(),nodal=this%VF%nodal)
            call tmpVF%setval(0.0_WP)
         end if

         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            if (dual_dm) then
               pVF=>tmpVF%dataptr(mfi)
            else
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
            end if
            call this%get_particles(lvl,mfi,p,np_)
            do i=1_I8,np_
               if (p(i)%flag.eq.PART_IS_DEAD) cycle
               ii=floor((p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ii,WP)
               jj=floor((p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jj,WP)
               kk=floor((p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kk,WP)
               ! Clamp 8-cell stencil at WALL faces (mirrors amrlpt)
               if (this%lo_bc(1).eq.AMRPD_WALL.and.ii  .lt.this%amr%geom(lvl)%domain%lo(1)) then; ii=this%amr%geom(lvl)%domain%lo(1)  ; wx=0.0_WP; end if
               if (this%hi_bc(1).eq.AMRPD_WALL.and.ii+1.gt.this%amr%geom(lvl)%domain%hi(1)) then; ii=this%amr%geom(lvl)%domain%hi(1)-1; wx=1.0_WP; end if
               if (this%lo_bc(2).eq.AMRPD_WALL.and.jj  .lt.this%amr%geom(lvl)%domain%lo(2)) then; jj=this%amr%geom(lvl)%domain%lo(2)  ; wy=0.0_WP; end if
               if (this%hi_bc(2).eq.AMRPD_WALL.and.jj+1.gt.this%amr%geom(lvl)%domain%hi(2)) then; jj=this%amr%geom(lvl)%domain%hi(2)-1; wy=1.0_WP; end if
               if (this%lo_bc(3).eq.AMRPD_WALL.and.kk  .lt.this%amr%geom(lvl)%domain%lo(3)) then; kk=this%amr%geom(lvl)%domain%lo(3)  ; wz=0.0_WP; end if
               if (this%hi_bc(3).eq.AMRPD_WALL.and.kk+1.gt.this%amr%geom(lvl)%domain%hi(3)) then; kk=this%amr%geom(lvl)%domain%hi(3)-1; wz=1.0_WP; end if
               pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)=pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)+Vp*reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])
            end do
         end do
         call this%mfiter_destroy(mfi)

         if (dual_dm) then
            call this%VF%mf(lvl)%parallel_copy(tmpVF,1,1,this%VF%mf(lvl)%ncomp(),this%VF%mf(lvl)%nghost(),this%VF%mf(lvl)%nghost(),this%amr%geom(lvl))
            call amrex_multifab_destroy(tmpVF)
         end if
      end do

      ! Convert extensive -> intensive (VF) and reconcile across C/F
      call this%process_deposit(this%VF)
      ! Fill ghost cells via amrdata's standard machinery
      call this%VF%fill(time=0.0_WP)
      ! Optional smoothing (zero filter_width disables)
      call this%filter(this%VF)

   contains

      !> Helper: get the particle distribution map for this level. Returns the
      !> DM that the particle container is currently using (which may differ
      !> from the Eulerian DM when rebalance reweights it).
      function get_pdm() result(dm)
         use amrex_distromap_module, only: amrex_distromap
         type(amrex_distromap) :: dm
         call amrpd_get_particle_distromap_p(this%pcp,lvl,dm%p)
         dm%owner=.false.
      end function get_pdm
   end subroutine update_VF

   !> Post-process an extensive deposit (sum of particle volumes per cell) into
   !> an intensive field (VF = sum/cell_vol), with cross-level transfers to
   !> avoid double-counting on covered cells. Verbatim port of amrlpt's
   !> process_deposit.
   subroutine process_deposit(this,A)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_sum_downto,amrmfab_interp_from_coarse
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      type(amrex_multifab), dimension(:), allocatable :: tmp
      integer :: lvl
      ! Convert extensive deposits to intensive
      do lvl=0,this%amr%clvl()
         call A%mf(lvl)%mult(1.0_WP/this%amr%cell_vol(lvl),1,A%ncomp,A%ng)
      end do
      ! Scratch mfabs for coarse->fine interpolation
      allocate(tmp(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,tmp(lvl),ncomp=A%ncomp,nover=0); call tmp(lvl)%setval(0.0_WP)
      end do
      ! Forward pass (coarse to fine)
      do lvl=0,this%amr%clvl()
         call A%syncsum_lvl(lvl)
         if (lvl.lt.this%amr%clvl()) then
            call amrmfab_interp_from_coarse(tmp(lvl+1),A%mf(lvl),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl),fgeom=this%amr%geom(lvl+1),scomp=1,ncomp=A%ncomp)
         end if
         if (lvl.gt.0) then
            call amrmfab_sum_downto(A%mf(lvl),A%mf(lvl-1),[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),fgeom=this%amr%geom(lvl))
         end if
         call A%mf(lvl)%add(tmp(lvl),1,1,A%ncomp,0)
      end do
      ! Backward pass: average down to fix double-counted covered cells
      do lvl=this%amr%clvl()-1,0,-1
         call A%average_downto(lvl)
      end do
      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(tmp(lvl))
      end do
      deallocate(tmp)
   end subroutine process_deposit

   !> Explicit-diffusion (Gaussian-equivalent) smoothing of a cell-centered
   !> amrdata field. Skips if filter_width<=mesh-size. Verbatim port of amrlpt.
   subroutine filter(this,A)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      real(WP) :: alpha,alpha_step,dxi,dyi,dzi
      integer  :: nstep,n,nc,lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pA,pFx,pFy,pFz

      alpha=max(this%filter_width**2-this%amr%min_meshsize(this%amr%clvl())**2,0.0_WP)/(16.0_WP*log(2.0_WP))
      if (alpha.le.0.0_WP) return

      nstep=ceiling(6.0_WP*alpha/this%amr%min_meshsize(this%amr%clvl())**2)
      alpha_step=alpha/real(nstep,WP)

      allocate(Fx(0:this%amr%maxlvl),Fy(0:this%amr%maxlvl),Fz(0:this%amr%maxlvl))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=A%ncomp,nover=0,atface=[.true., .false.,.false.]); call Fx(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.true., .false.]); call Fy(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
      end do

      do n=1,nstep
         do lvl=0,this%amr%clvl()
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call amrex_mfiter_build(mfi,A%mf(lvl),tiling=.false.)
            do while (mfi%next())
               pA =>A%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(1)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFx(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i-1,j,k,nc))*dxi
               end do; end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFy(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i,j-1,k,nc))*dyi
               end do; end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFz(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i,j,k-1,nc))*dzi
               end do; end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
         end do
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         end do
         do lvl=0,this%amr%clvl()
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call amrex_mfiter_build(mfi,A%mf(lvl),tiling=.false.)
            do while (mfi%next())
               pA =>A%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%tilebox()
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pA(i,j,k,nc)=pA(i,j,k,nc)+dxi*(pFx(i+1,j,k,nc)-pFx(i,j,k,nc))+dyi*(pFy(i,j+1,k,nc)-pFy(i,j,k,nc))+dzi*(pFz(i,j,k+1,nc)-pFz(i,j,k,nc))
               end do; end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
         end do
         call A%average_down()
         call A%fill(time=0.0_WP)
      end do

      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(Fx(lvl))
         call amrex_multifab_destroy(Fy(lvl))
         call amrex_multifab_destroy(Fz(lvl))
      end do
      deallocate(Fx,Fy,Fz)
   end subroutine filter


   ! ============================================================================
   ! PHYSICS -- STUBBED IN SKELETON
   ! ============================================================================

   !> Build initial bonds and stamp reference state.
   !>
   !> Phase 1 (detect): fill_ghosts(delta), then for each owned particle i pair
   !>   with every (owned + ghost) particle j within delta. Keep the pair only if
   !>   i has the lower-GID (the lower-GID endpoint owns the bond). Append the
   !>   bond at pos(i) with d0 = |x_j - x_i|, w = omega(d0, delta), damage=0,
   !>   alive=1. The two endpoints' AMReX idcpu values are split into 2 x int32
   !>   and stored in the bond's idata.
   !>
   !> Phase 2 (m_w stamp): walk bonds, scatter w*d0^2*dV to both endpoints' mw
   !>   slot. Owned-endpoint writes go directly to the ptile real slots; ghost-
   !>   endpoint writes go to the AMReX-internal neighbor (aux) buffer because
   !>   sumNeighbors reads its reduction inputs from there, not from the ptile's
   !>   appended ghost slots. sum_ghosts_mw then reduces ghost contributions
   !>   back to owners. mw is reference-only — set once here, never updated again.
   !>
   !> Influence function is a host-associated pure function (contained below) so
   !> the compiler can inline it. To switch to the gaussian (lss) form, comment
   !> the constant line and uncomment the gaussian one.
   subroutine bond_init(this)
      use amrex_amr_module, only: amrex_mfiter
      use amrpd_hash_class, only: gid_hash
      use string,           only: str_long
      use messager,         only: log
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP) :: r2_cut,K_bulk

      r2_cut=this%delta**2

      ! Critical bond stretch for brittle damage, derived from G_c
      ! (Silling-Askari 3D): s0 = sqrt(5 * G_c / (9 * K * delta)).
      ! If user did not provide G_c, leave s0 at its huge() default (no damage).
      if (this%crit_energy.gt.0.0_WP) then
         K_bulk=this%elastic_modulus/(3.0_WP*(1.0_WP-2.0_WP*this%poisson_ratio))
         this%s0=sqrt(5.0_WP*this%crit_energy/(9.0_WP*K_bulk*this%delta))
      end if
      ! Direct failure-stretch override: decouples rupture from the brittle G_c
      ! value (large -> ductile). Wins over the G_c-derived s0 when set.
      if (this%fail_stretch.lt.huge(1.0_WP)) this%s0=this%fail_stretch

      ! Search radius MUST cover the largest a LIVE bond can stretch (= failure
      ! stretch s0): a bond stretched beyond that has already ruptured, so its
      ! partner is never needed. This keeps every live bond's partner inside the
      ! ghost layer under large (ductile) deformation. Enlarge only; never shrink
      ! (and the radius is fixed for the run -- the AMReX neighbor mask is sized
      ! on the first fill_ghosts below). 1.2 = one-step drift margin.
      if (this%s0.lt.huge(1.0_WP)) then
         this%search_radius=max(this%search_radius,(1.0_WP+this%s0)*this%delta*1.2_WP)
      else if (this%tau.lt.huge(1.0_WP)) then
         call log('[amrpd] WARNING: viscoelastic flow (finite tau) with no failure stretch (s0=huge) &
         &-- bonds may stretch beyond search_radius and be silently dropped; set Critical energy or Failure stretch.')
      end if

      ! Short-range contact: only contact_dist gets defaulted here. The
      ! collision duration tau_col is set fresh each step inside compute_contact
      ! as 5*dt -- "as stiff as the current timestep allows" -- so it doesn't
      ! need a one-shot default. Other fields (e_n, e_w, clip_col) take their
      ! declared defaults.
      if (this%contact_dist.le.0.0_WP) this%contact_dist=0.9_WP*this%dV**(1.0_WP/3.0_WP)

      ! Fill ghost particles within search_radius of each tile boundary. We
      ! always use the same search_radius for every fill_ghosts call — AMReX's
      ! neighbor mask is sized on the FIRST fill_neighbors call and not rebuilt
      ! on subsequent calls with larger radius, so a single fixed radius across
      ! the run avoids out-of-bounds mask accesses.
      call this%fill_ghosts(radius=this%search_radius)

      ! Phase 1: detect candidate bonds, accumulate them in a local scratch
      ! buffer, then bulk-append. Each bond is created exactly once, by the rank
      ! that owns the bond's lower-GID endpoint.
      detect_bonds: block
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         type(bond), dimension(:), allocatable, target :: blist
         integer(I8) :: np_total,np_valid,i,j
         integer(I8) :: nb_local,ncap
         integer :: lvl,nx,ny,nz
         integer(c_int64_t) :: key_i,key_j
         integer(c_int) :: parts_lo(2),parts_hi(2)
         real(WP) :: r2,dist,Lx,Ly,Lz

         ncap=1024_I8
         allocate(blist(ncap))
         nb_local=0_I8
         Lx=this%amr%xhi-this%amr%xlo; Ly=this%amr%yhi-this%amr%ylo; Lz=this%amr%zhi-this%amr%zlo

         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_all_particles(lvl,mfi,p,np_total,np_valid)
               ! Outer loop: only OWNED particles initiate. Inner: all (owned+ghost).
               do i=1_I8,np_valid
                  key_i=p(i)%idcpu
                  do j=1_I8,np_total
                     if (i.eq.j) cycle
                     r2=sum((p(j)%pos-p(i)%pos)**2)
                     if (r2.gt.r2_cut) cycle
                     key_j=p(j)%idcpu
                     ! Lower-GID owns the bond. Skip when i is strictly higher;
                     ! equal ids are kept so a particle bonds to its own periodic
                     ! images (self-image bonds when period < horizon).
                     if (key_i.gt.key_j) cycle
                     ! Grow buffer if needed
                     nb_local=nb_local+1_I8
                     if (nb_local.gt.ncap) then
                        grow: block
                           type(bond), dimension(:), allocatable, target :: tmp
                           ncap=2_I8*ncap
                           allocate(tmp(ncap))
                           tmp(1:size(blist))=blist
                           call move_alloc(tmp,blist)
                        end block grow
                     end if
                     ! Stamp the bond. hist1 packs the higher endpoint's periodic
                     ! image offset (n_x,n_y,n_z) so the exact bonded image is
                     ! reconstructed at force time instead of guessed by proximity.
                     nx=0; ny=0; nz=0
                     if (this%amr%xper) nx=floor((p(j)%pos(1)-this%amr%xlo)/Lx)
                     if (this%amr%yper) ny=floor((p(j)%pos(2)-this%amr%ylo)/Ly)
                     if (this%amr%zper) nz=floor((p(j)%pos(3)-this%amr%zlo)/Lz)
                     dist=sqrt(r2)
                     blist(nb_local)%pos    =p(i)%pos
                     blist(nb_local)%d0     =dist
                     blist(nb_local)%w      =w(dist,this%delta)
                     blist(nb_local)%damage =0.0_WP
                     blist(nb_local)%hist1  =real((nx+128)+(ny+128)*256+(nz+128)*65536,WP)
                     blist(nb_local)%e_v    =0.0_WP
                     parts_lo=transfer(key_i,parts_lo)
                     parts_hi=transfer(key_j,parts_hi)
                     blist(nb_local)%id_lo_lo=parts_lo(1)
                     blist(nb_local)%id_lo_hi=parts_lo(2)
                     blist(nb_local)%id_hi_lo=parts_hi(1)
                     blist(nb_local)%id_hi_hi=parts_hi(2)
                     blist(nb_local)%alive  =1
                  end do
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do

         ! Bulk-append all collected bonds (collective; AddParticlesAtLevel
         ! redistributes by position so they land at their lower-GID's tile)
         call this%append_bonds(blist,nb_local)
         deallocate(blist)
      end block detect_bonds
      call this%clear_ghosts()

      ! Refresh global bond count for logging
      call this%get_info()

      ! Phase 2: stamp reference-state quantities m_w (weighted volume) AND nb0
      ! (raw bond count) on every particle. Both are reference-only: set once
      ! here, never updated by damage. Single bond sweep scatters w*d0^2*dV to
      ! mw and +1 to nb0 at each endpoint; sum_ghosts reduces ghost contributions.
      compute_mw: block
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p,pg
         type(bond), dimension(:), pointer :: b
         type(gid_hash) :: hash
         integer(I8) :: np_total,np_valid,ng,nb_tile,n,ib
         integer :: lvl,np_int,lid_lo,lid_hi
         integer(c_int64_t), allocatable :: keys(:)
         integer(c_int64_t) :: key_lo,key_hi
         integer(c_int) :: parts(2)
         real(WP) :: contrib

         ! Re-fill ghosts; the bond container's redistribute may have invalidated
         ! the particle ghost layer (bonds and particles share fillNeighbors state
         ! only loosely — safest to refresh)
         call this%fill_ghosts(radius=this%search_radius)

         ! Zero mw and nb0 on every owned particle (ptile real slots) AND every
         ! ghost particle in the aux buffer (what sum_ghosts_* reads).
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_valid)
               do n=1_I8,np_valid
                  p(n)%mw =0.0_WP
                  p(n)%nb0=0.0_WP
               end do
               call this%get_ghosts(lvl,mfi,pg,ng)
               do n=1_I8,ng
                  pg(n)%mw =0.0_WP
                  pg(n)%nb0=0.0_WP
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do

         ! Walk bonds; scatter mw to both endpoints. Build hash over the
         ! COMBINED array (real + ghost) for GID->LID lookup. Writes to ghost
         ! endpoints (LID > np_valid) target the AUX buffer so sumNeighbors picks
         ! them up; writes to owned endpoints (LID <= np_valid) go directly to p.
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_all_particles(lvl,mfi,p,np_total,np_valid)
               call this%get_ghosts(lvl,mfi,pg,ng)
               ! Build per-tile GID->LID hash from the combined array
               np_int=int(np_total)
               allocate(keys(np_int))
               do n=1_I8,np_total
                  keys(n)=p(n)%idcpu
               end do
               call hash%build(np_int,keys)
               deallocate(keys)
               ! Iterate this tile's bonds
               call this%get_bonds(lvl,mfi,b,nb_tile)
               do ib=1_I8,nb_tile
                  ! Reconstruct endpoint keys from packed idata
                  parts(1)=b(ib)%id_lo_lo; parts(2)=b(ib)%id_lo_hi
                  key_lo=transfer(parts,0_c_int64_t)
                  parts(1)=b(ib)%id_hi_lo; parts(2)=b(ib)%id_hi_hi
                  key_hi=transfer(parts,0_c_int64_t)
                  ! Resolve LIDs via the hash (1-based in p; > np_valid = ghost).
                  ! Periodic-image disambiguation NOT needed here: mw/nb0
                  ! contributions are direction-invariant scalars, so scattering
                  ! to any periodic image of an endpoint reduces back to the
                  ! same owner via sum_ghosts_mw / sum_ghosts_nb0.
                  lid_lo=hash%lookup(key_lo)
                  lid_hi=hash%lookup(key_hi)
                  if (lid_lo.lt.1.or.lid_hi.lt.1) cycle   ! shouldn't happen
                  ! Contributions: mw += w*|xi|^2*V (weighted volume), nb0 += 1
                  ! (raw bond count, used to normalize the damage fraction).
                  contrib=b(ib)%w*b(ib)%d0**2*this%dV
                  if (lid_lo.le.int(np_valid)) then
                     p(lid_lo)%mw =p(lid_lo)%mw +contrib
                     p(lid_lo)%nb0=p(lid_lo)%nb0+1.0_WP
                  else
                     pg(lid_lo-int(np_valid))%mw =pg(lid_lo-int(np_valid))%mw +contrib
                     pg(lid_lo-int(np_valid))%nb0=pg(lid_lo-int(np_valid))%nb0+1.0_WP
                  end if
                  ! Self-image bonds (id_lo==id_hi) scatter to lo only; the
                  ! opposite-side self-image bond supplies the hi contribution.
                  if (key_lo.ne.key_hi) then
                     if (lid_hi.le.int(np_valid)) then
                        p(lid_hi)%mw =p(lid_hi)%mw +contrib
                        p(lid_hi)%nb0=p(lid_hi)%nb0+1.0_WP
                     else
                        pg(lid_hi-int(np_valid))%mw =pg(lid_hi-int(np_valid))%mw +contrib
                        pg(lid_hi-int(np_valid))%nb0=pg(lid_hi-int(np_valid))%nb0+1.0_WP
                     end if
                  end if
               end do
               call hash%finalize()
            end do
            call this%mfiter_destroy(mfi)
         end do

         ! Reduce ghost-slot mw and nb0 contributions (from aux buffer) back to owners
         call this%sum_ghosts_mw()
         call this%sum_ghosts_nb0()
      end block compute_mw
      call this%clear_ghosts()


      ! Log bond count
      log_bonds: block
         character(len=str_long) :: message
         if (this%amr%amRoot) then
            write(message,'("[",a,"] bond_init: ",i0," bonds created (horizon=",es12.5,")")') &
            &     trim(this%name),this%nb,this%delta
            call log(message)
         end if
      end block log_bonds

   contains

      !> Influence function w(d, h). Host-associated pure function — compiler
      !> inlines it where called above. Pick one of the two bodies:
      !>   constant (Peridigm default)  : omega = 1
      !>   gaussian/quartic (lss form) : (1+4d/h)(1-d/h)^4, zero for d>=h
      pure function w(d,h) result(omega)
         real(WP), intent(in) :: d,h
         real(WP) :: omega
         omega=1.0_WP                                                                            ! constant (Peridigm)
         ! omega=merge(0.0_WP,(1.0_WP+4.0_WP*d/h)*(1.0_WP-d/h)**4,d.ge.h)                        ! gaussian (lss)
      end function w

   end subroutine bond_init

   !> Compute LPS dilatation theta on every owned particle:
   !>   theta_i = (3 / m_w_i) * sum_{j in H_i} w_ij * |xi_ij| * e_ij * V_j
   !> with e_ij = |x_j - x_i| - |X_j - X_i| (current minus reference bond length).
   !>
   !> Sweep bonds; for each alive bond, compute the current length using the
   !> particle pointer (combined real+ghost view) and scatter w*d0*e*dV to
   !> BOTH endpoints' dil slot. Owned-endpoint writes go to p; ghost-endpoint
   !> writes go to the aux buffer pg (so sum_ghosts_dil can reduce them).
   !> Finally divide owned dil by m_w and multiply by 3.
   !>
   !> Caller must have called fill_ghosts(search_radius) with the CURRENT
   !> particle positions before invoking this routine.
   subroutine compute_dilatation(this)
      use amrex_amr_module, only: amrex_mfiter
      use amrpd_hash_class, only: gid_hash
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p,pg
      type(bond), dimension(:), pointer :: b
      type(gid_hash) :: hash
      integer(I8) :: np_total,np_valid,ng,nb_tile,n,ib
      integer :: lvl,np_int,lid_lo,lid_hi,nx,ny,nz,ip
      integer(c_int64_t), allocatable :: keys(:)
      integer(c_int64_t) :: key_lo,key_hi
      integer(c_int) :: parts(2)
      real(WP) :: dx,dy,dz,curr_len,e_bond,contrib,Lx,Ly,Lz
      real(WP), dimension(3) :: xlo,xhi

      Lx=this%amr%xhi-this%amr%xlo; Ly=this%amr%yhi-this%amr%ylo; Lz=this%amr%zhi-this%amr%zlo

      ! Zero dil on every owned particle AND every ghost in the aux buffer.
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_valid)
            do n=1_I8,np_valid
               p(n)%dil=0.0_WP
            end do
            call this%get_ghosts(lvl,mfi,pg,ng)
            do n=1_I8,ng
               pg(n)%dil=0.0_WP
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do

      ! Walk bonds; scatter w*d0*e*dV to both endpoints
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_all_particles(lvl,mfi,p,np_total,np_valid)
            call this%get_ghosts(lvl,mfi,pg,ng)
            ! Build per-tile GID->LID hash from the combined array
            np_int=int(np_total)
            allocate(keys(np_int))
            do n=1_I8,np_total
               keys(n)=p(n)%idcpu
            end do
            call hash%build(np_int,keys)
            deallocate(keys)
            ! Iterate this tile's bonds
            call this%get_bonds(lvl,mfi,b,nb_tile)
            do ib=1_I8,nb_tile
               if (b(ib)%alive.eq.0) cycle
               ! Reconstruct endpoint keys from packed idata
               parts(1)=b(ib)%id_lo_lo; parts(2)=b(ib)%id_lo_hi
               key_lo=transfer(parts,0_c_int64_t)
               parts(1)=b(ib)%id_hi_lo; parts(2)=b(ib)%id_hi_hi
               key_hi=transfer(parts,0_c_int64_t)
               ! Resolve endpoint LIDs (any copy; the scatter is folded to owners
               ! by sum_ghosts). The exact bonded image of the higher endpoint is
               ! reconstructed from the stored periodic offset, not guessed.
               lid_lo=hash%lookup(key_lo)
               lid_hi=hash%lookup(key_hi)
               if (lid_lo.lt.1.or.lid_hi.lt.1) cycle
               ip=nint(b(ib)%hist1)
               nx=mod(ip,256)-128; ny=mod(ip/256,256)-128; nz=ip/65536-128
               xlo=canon(p(lid_lo)%pos)
               xhi=canon(p(lid_hi)%pos)
               xhi(1)=xhi(1)+real(nx,WP)*Lx; xhi(2)=xhi(2)+real(ny,WP)*Ly; xhi(3)=xhi(3)+real(nz,WP)*Lz
               ! Current deformed bond length and extension
               dx=xhi(1)-xlo(1); dy=xhi(2)-xlo(2); dz=xhi(3)-xlo(3)
               curr_len=sqrt(dx*dx+dy*dy+dz*dz)
               e_bond=curr_len-b(ib)%d0
               contrib=b(ib)%w*b(ib)%d0*e_bond*this%dV
               if (lid_lo.le.int(np_valid)) then
                  p(lid_lo)%dil=p(lid_lo)%dil+contrib
               else
                  pg(lid_lo-int(np_valid))%dil=pg(lid_lo-int(np_valid))%dil+contrib
               end if
               ! Self-image bonds scatter to lo only (see compute_mw rationale)
               if (key_lo.ne.key_hi) then
                  if (lid_hi.le.int(np_valid)) then
                     p(lid_hi)%dil=p(lid_hi)%dil+contrib
                  else
                     pg(lid_hi-int(np_valid))%dil=pg(lid_hi-int(np_valid))%dil+contrib
                  end if
               end if
            end do
            call hash%finalize()
         end do
         call this%mfiter_destroy(mfi)
      end do

      ! Reduce ghost-slot dil contributions back to owners
      call this%sum_ghosts_dil()

      ! Finalize: dil_i = 3 * (accumulated sum) / m_w_i  (owned only)
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_valid)
            do n=1_I8,np_valid
               if (p(n)%mw.gt.0.0_WP) then
                  p(n)%dil=3.0_WP*p(n)%dil/p(n)%mw
               else
                  p(n)%dil=0.0_WP
               end if
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do

   contains

      !> Wrap a position into the base domain [lo,hi) along periodic directions.
      !> Used with the per-bond stored offset to reconstruct the exact bonded
      !> periodic image (host-associated Lx/Ly/Lz).
      function canon(pos) result(c)
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: c
         c=pos
         if (this%amr%xper) c(1)=pos(1)-Lx*floor((pos(1)-this%amr%xlo)/Lx)
         if (this%amr%yper) c(2)=pos(2)-Ly*floor((pos(2)-this%amr%ylo)/Ly)
         if (this%amr%zper) c(3)=pos(3)-Lz*floor((pos(3)-this%amr%zlo)/Lz)
      end function canon

   end subroutine compute_dilatation

   !> Compute LPS bond forces and accumulate F_bond (force/volume) on every
   !> owned particle. LPS pairwise force density at endpoint i toward j
   !> (Silling 2007, matches Peridigm/elastic.cxx):
   !>   t_i = w/m_w_i * [ (3K - 5mu) * theta_i * |xi|  +  15 mu * e ]
   !> with e = |x_j - x_i| - |X_j - X_i| (full bond extension, NOT decomposed).
   !> Newton's-3rd-law pair force on lo from bond:
   !>   f_lo = (t_lo + t_hi) * dV * M_hat,  M_hat = (x_hi - x_lo)/|x_hi - x_lo|
   !>   f_hi = -f_lo
   !>
   !> Bonds with stretch e/d0 > s0 are irreversibly broken (alive=0, damage=1)
   !> and contribute zero force from that step onward. m_w stays at its initial
   !> value (reference quantity), so damage softens the material naturally.
   !>
   !> Caller must have (in order):
   !>   1) called fill_ghosts(search_radius) with current positions
   !>   2) called compute_dilatation
   !>   3) called update_ghosts so aux-buffer dil/mw match current owner values
   subroutine compute_force(this,dt)
      use amrex_amr_module, only: amrex_mfiter
      use amrpd_hash_class, only: gid_hash
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p,pg
      type(bond), dimension(:), pointer :: b
      type(gid_hash) :: hash
      integer(I8) :: np_total,np_valid,ng,nb_tile,n,ib
      integer :: lvl,np_int,lid_lo,lid_hi,nx,ny,nz,ip
      integer(c_int64_t), allocatable :: keys(:)
      integer(c_int64_t) :: key_lo,key_hi
      integer(c_int) :: parts(2)
      real(WP) :: K_bulk,mu_shear,coef_vol,coef_dev,e_d_avg,decay
      real(WP) :: dx,dy,dz,curr_len,e_bond,Lx,Ly,Lz
      real(WP), dimension(3) :: xlo,xhi
      real(WP) :: t_lo,t_hi,pair_mag
      real(WP) :: fx,fy,fz,mhat_x,mhat_y,mhat_z

      Lx=this%amr%xhi-this%amr%xlo; Ly=this%amr%yhi-this%amr%ylo; Lz=this%amr%zhi-this%amr%zlo

      ! Elastic moduli from (E, nu)
      K_bulk  =this%elastic_modulus/(3.0_WP*(1.0_WP-2.0_WP*this%poisson_ratio))
      mu_shear=this%elastic_modulus/(2.0_WP*(1.0_WP+this%poisson_ratio))
      coef_vol=3.0_WP*K_bulk                       ! volumetric (elastic)
      coef_dev=15.0_WP*mu_shear                    ! deviatoric (carries Maxwell relaxation; e_v=0 recovers LPS)

      ! Zero F_bond on every owned particle AND every ghost in the aux buffer.
      ! Also zero damage on GHOSTS only (owned damage is accumulated state and
      ! must persist; ghost slots carry only this step's break increments,
      ! which will be reduced to owners by sum_ghosts_damage at the end).
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_valid)
            do n=1_I8,np_valid
               p(n)%F_bond=0.0_WP
            end do
            call this%get_ghosts(lvl,mfi,pg,ng)
            do n=1_I8,ng
               pg(n)%F_bond=0.0_WP
               pg(n)%damage=0.0_WP
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do

      ! Walk bonds; compute pair force, scatter to both endpoints
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_all_particles(lvl,mfi,p,np_total,np_valid)
            call this%get_ghosts(lvl,mfi,pg,ng)
            np_int=int(np_total)
            allocate(keys(np_int))
            do n=1_I8,np_total
               keys(n)=p(n)%idcpu
            end do
            call hash%build(np_int,keys)
            deallocate(keys)
            call this%get_bonds(lvl,mfi,b,nb_tile)
            do ib=1_I8,nb_tile
               if (b(ib)%alive.eq.0) cycle
               parts(1)=b(ib)%id_lo_lo; parts(2)=b(ib)%id_lo_hi
               key_lo=transfer(parts,0_c_int64_t)
               parts(1)=b(ib)%id_hi_lo; parts(2)=b(ib)%id_hi_hi
               key_hi=transfer(parts,0_c_int64_t)
               ! Resolve LIDs and reconstruct the exact bonded image from the
               ! stored periodic offset (see compute_dilatation)
               lid_lo=hash%lookup(key_lo)
               lid_hi=hash%lookup(key_hi)
               if (lid_lo.lt.1.or.lid_hi.lt.1) cycle
               ip=nint(b(ib)%hist1)
               nx=mod(ip,256)-128; ny=mod(ip/256,256)-128; nz=ip/65536-128
               xlo=canon(p(lid_lo)%pos)
               xhi=canon(p(lid_hi)%pos)
               xhi(1)=xhi(1)+real(nx,WP)*Lx; xhi(2)=xhi(2)+real(ny,WP)*Ly; xhi(3)=xhi(3)+real(nz,WP)*Lz
               ! Deformed bond vector, length, extension
               dx=xhi(1)-xlo(1); dy=xhi(2)-xlo(2); dz=xhi(3)-xlo(3)
               curr_len=sqrt(dx*dx+dy*dy+dz*dz)
               if (curr_len.le.0.0_WP) cycle
               e_bond=curr_len-b(ib)%d0
               ! Brittle damage: irreversibly break the bond if its stretch
               ! exceeds the critical stretch s0 (derived from G_c in bond_init,
               ! or huge() = disabled if G_c was not provided). Broken bonds
               ! contribute zero force from this step onward. Per-particle
               ! damage is INCREMENTED here (at break time) by 1/nb0 on BOTH
               ! endpoints -- accumulated on the particle's own rdata, so it
               ! survives any subsequent redistribute and stays monotonic.
               ! Ghost-side increments are reduced to their owners via
               ! sum_ghosts_damage at the end of this routine.
               if (e_bond.gt.this%s0*b(ib)%d0) then
                  b(ib)%alive =0
                  b(ib)%damage=1.0_WP
                  ! Lower-GID is owned on this rank (bond is co-located with it)
                  if (p(lid_lo)%nb0.gt.0.0_WP) p(lid_lo)%damage=p(lid_lo)%damage+1.0_WP/p(lid_lo)%nb0
                  ! Upper-GID (skip for self-image bonds; lo already counted it)
                  if (key_lo.ne.key_hi) then
                     if (lid_hi.le.int(np_valid)) then
                        if (p(lid_hi)%nb0.gt.0.0_WP) p(lid_hi)%damage=p(lid_hi)%damage+1.0_WP/p(lid_hi)%nb0
                     else
                        if (pg(lid_hi-int(np_valid))%nb0.gt.0.0_WP) &
                        &  pg(lid_hi-int(np_valid))%damage=pg(lid_hi-int(np_valid))%damage+1.0_WP/pg(lid_hi-int(np_valid))%nb0
                     end if
                  end if
                  cycle
               end if
               ! Force-density magnitudes along M_hat at each endpoint
               ! (uses each owner's own theta and m_w; e is symmetric)
               ! Volumetric part elastic; deviatoric extension e_d=e-theta*d0/3 carries
               ! the Maxwell inelastic stretch e_v. e_v=0 -> identical to the LPS form.
               if (p(lid_lo)%mw.gt.0.0_WP) then
                  t_lo=b(ib)%w/p(lid_lo)%mw*(coef_vol*p(lid_lo)%dil*b(ib)%d0 &
                  &    +coef_dev*(e_bond-p(lid_lo)%dil*b(ib)%d0/3.0_WP-this%visc_lambda*b(ib)%e_v))
               else
                  t_lo=0.0_WP
               end if
               if (p(lid_hi)%mw.gt.0.0_WP) then
                  t_hi=b(ib)%w/p(lid_hi)%mw*(coef_vol*p(lid_hi)%dil*b(ib)%d0 &
                  &    +coef_dev*(e_bond-p(lid_hi)%dil*b(ib)%d0/3.0_WP-this%visc_lambda*b(ib)%e_v))
               else
                  t_hi=0.0_WP
               end if
               ! Unit bond vector and pair-force-per-volume (Newton's 3rd)
               mhat_x=dx/curr_len; mhat_y=dy/curr_len; mhat_z=dz/curr_len
               pair_mag=(t_lo+t_hi)*this%dV
               fx=pair_mag*mhat_x; fy=pair_mag*mhat_y; fz=pair_mag*mhat_z
               ! Scatter +f to lo, -f to hi (route ghost writes to aux buffer)
               if (lid_lo.le.int(np_valid)) then
                  p(lid_lo)%F_bond(1)=p(lid_lo)%F_bond(1)+fx
                  p(lid_lo)%F_bond(2)=p(lid_lo)%F_bond(2)+fy
                  p(lid_lo)%F_bond(3)=p(lid_lo)%F_bond(3)+fz
               else
                  pg(lid_lo-int(np_valid))%F_bond(1)=pg(lid_lo-int(np_valid))%F_bond(1)+fx
                  pg(lid_lo-int(np_valid))%F_bond(2)=pg(lid_lo-int(np_valid))%F_bond(2)+fy
                  pg(lid_lo-int(np_valid))%F_bond(3)=pg(lid_lo-int(np_valid))%F_bond(3)+fz
               end if
               ! Self-image bonds apply +f to lo only; the opposite-side self-image
               ! bond supplies the reaction (no -f onto the same owner DOF).
               if (key_lo.ne.key_hi) then
                  if (lid_hi.le.int(np_valid)) then
                     p(lid_hi)%F_bond(1)=p(lid_hi)%F_bond(1)-fx
                     p(lid_hi)%F_bond(2)=p(lid_hi)%F_bond(2)-fy
                     p(lid_hi)%F_bond(3)=p(lid_hi)%F_bond(3)-fz
                  else
                     pg(lid_hi-int(np_valid))%F_bond(1)=pg(lid_hi-int(np_valid))%F_bond(1)-fx
                     pg(lid_hi-int(np_valid))%F_bond(2)=pg(lid_hi-int(np_valid))%F_bond(2)-fy
                     pg(lid_hi-int(np_valid))%F_bond(3)=pg(lid_hi-int(np_valid))%F_bond(3)-fz
                  end if
               end if
               ! Maxwell relaxation of the bond's inelastic deviatoric stretch
               ! (owner-local; tau=huge -> e_v frozen -> purely elastic LPS).
               ! Exact exponential integration (Peridigm-style) -> unconditionally
               ! stable for any dt, so no viscous CFL constraint.
               if (this%tau.gt.0.0_WP.and.this%tau.lt.huge(1.0_WP)) then
                  e_d_avg=e_bond-0.5_WP*(p(lid_lo)%dil+p(lid_hi)%dil)*b(ib)%d0/3.0_WP
                  decay=exp(-dt/this%tau)
                  b(ib)%e_v=e_d_avg*(1.0_WP-decay)+b(ib)%e_v*decay
               end if
            end do
            call hash%finalize()
         end do
         call this%mfiter_destroy(mfi)
      end do

      ! Reduce ghost-slot F_bond AND damage contributions back to owners.
      ! (Damage was accumulated on ghost slots at bond-break time above; this
      ! routes those increments to the upper-GID's owner.)
      call this%sum_ghosts_force()
      call this%sum_ghosts_damage()

   contains

      !> Wrap a position into the base domain [lo,hi) along periodic directions
      !> (see compute_dilatation). Host-associated Lx/Ly/Lz.
      function canon(pos) result(c)
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: c
         c=pos
         if (this%amr%xper) c(1)=pos(1)-Lx*floor((pos(1)-this%amr%xlo)/Lx)
         if (this%amr%yper) c(2)=pos(2)-Ly*floor((pos(2)-this%amr%ylo)/Ly)
         if (this%amr%zper) c(3)=pos(3)-Lz*floor((pos(3)-this%amr%zlo)/Lz)
      end function canon

   end subroutine compute_force

   !> Short-range contact: soft-sphere penalty + damping, ported from
   !> amrlpt%collide. Handles three kinds of interactions:
   !>   - axis-aligned WALL faces of the domain (lo_bc/hi_bc == AMRPD_WALL)
   !>   - optional immersed boundary via Gib (level-set, +ve inside fluid)
   !>   - particle-particle pairs from a cell-binned neighbor list
   !>
   !> Contact duration is fixed at 5*dt every step -- "as stiff as the current
   !> timestep allows". Not user-tunable, not CFL-constrained: dt is set by the
   !> bond (wave) and particle (convective) CFLs alone; the contact spring then
   !> rides on whatever dt comes out.
   !>
   !> Smart features carried over from amrlpt%col_force:
   !>   - k/eta from (tau_col, e_n/e_w): k = (pi^2+ln(e)^2)/tau^2, eta = -2 ln(e)/tau
   !>   - early-engagement inflation r_influ = min(|rnv|*dt, 0.2*d_eff)
   !>   - overlap clipping delta_n <= clip_col*d_eff
   !>   - normal force f_n = (-m_eff*k*delta_n - m_eff*eta*rnv) * n12
   !>
   !> Walls/IB use e_w (m2=inf => m_eff=m1, d_eff=0.5*contact_dist).
   !> Particle-particle uses e_n (m_eff=m1/2, d_eff=contact_dist).
   !>
   !> Particle-particle uses one-sided pattern (each owned i accumulates from
   !> its neighbors; j picks up its share when iterated as i). No scatter, no
   !> sum_ghosts_force. Adds to F_bond as force/volume. Caller must have called
   !> fill_ghosts upstream.
   subroutine compute_contact(this,dt,Gib,Gibcomp)
      use amrex_amr_module, only: amrex_mfiter
      use amrdata_class,    only: amrdata
      use mathtools,        only: Pi
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in), optional :: Gib
      integer, intent(in), optional :: Gibcomp
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pG
      integer(c_int32_t), dimension(:), pointer :: nbr_off,nbr_lst
      integer(I8) :: np_total,np_valid
      integer :: lvl,gc
      integer(c_int32_t) :: i1,j_c,k
      real(WP) :: tau_col,k_n,eta_n,k_w,eta_w,m1
      real(WP) :: dxm,dym,dzm,d_eff_w
      real(WP), dimension(3) :: r1,v1,r2,v2,f_local

      ! Bail out if contact disabled or restitution out of range
      if (this%contact_dist.le.0.0_WP.or.dt.le.0.0_WP) return
      if (this%e_n.le.0.0_WP.or.this%e_w.le.0.0_WP) return

      ! Contact duration: user override (this%tau_col > 0) or auto = 5*dt.
      ! Auto means contact is as stiff as the current dt allows -- contact
      ! rides on whatever dt the wave/convective CFLs set.
      if (this%tau_col.gt.0.0_WP) then
         tau_col=this%tau_col
      else
         tau_col=5.0_WP*dt
      end if

      ! Spring + damping coefficients (particle-particle and wall/IB use
      ! separate restitution coefficients; same collision duration tau_col)
      k_n  =(Pi**2+log(this%e_n)**2)/tau_col**2
      eta_n=-2.0_WP*log(this%e_n)/tau_col
      k_w  =(Pi**2+log(this%e_w)**2)/tau_col**2
      eta_w=-2.0_WP*log(this%e_w)/tau_col

      ! IB component selector and per-particle (uniform) mass
      gc=1; if (present(Gibcomp)) gc=Gibcomp
      m1=this%rho*this%dV

      ! Wall/IB use d_eff = 0.5*contact_dist (single particle radius)
      d_eff_w=0.5_WP*this%contact_dist

      ! Build a separate (contact) neighbor list with rcrit slightly inflated
      ! to accommodate the r_influ early-engagement bump (capped at 0.2*d_eff,
      ! so 1.2*contact_dist is safe).
      call this%build_neighbor_list(rcrit=1.2_WP*this%contact_dist)

      do lvl=0,this%amr%clvl()
         dxm=this%amr%dx(lvl); dym=this%amr%dy(lvl); dzm=this%amr%dz(lvl)
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_all_particles(lvl,mfi,p,np_total,np_valid)
            call this%get_neighbor_list(lvl,mfi,nbr_off,nbr_lst)
            if (present(Gib)) pG=>Gib%mf(lvl)%dataptr(mfi)

            do i1=1,int(np_valid,c_int32_t)
               r1=p(i1)%pos; v1=p(i1)%vel
               f_local=0.0_WP
               v2=0.0_WP

               ! Wall collisions on axis-aligned faces flagged AMRPD_WALL.
               ! Virtual partner sits on the wall directly normal to particle.
               if (this%lo_bc(1).eq.AMRPD_WALL) then; r2=[this%amr%xlo,r1(2),r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2); end if
               if (this%hi_bc(1).eq.AMRPD_WALL) then; r2=[this%amr%xhi,r1(2),r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2); end if
               if (this%lo_bc(2).eq.AMRPD_WALL) then; r2=[r1(1),this%amr%ylo,r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2); end if
               if (this%hi_bc(2).eq.AMRPD_WALL) then; r2=[r1(1),this%amr%yhi,r1(3)]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2); end if
               if (this%lo_bc(3).eq.AMRPD_WALL) then; r2=[r1(1),r1(2),this%amr%zlo]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2); end if
               if (this%hi_bc(3).eq.AMRPD_WALL) then; r2=[r1(1),r1(2),this%amr%zhi]; call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2); end if

               ! Optional IB collision (virtual partner at IB surface)
               if (present(Gib)) then
                  ib_block: block
                     real(WP) :: d_ib,nrm
                     real(WP), dimension(3) :: pos_p,pos_m,n_out
                     d_ib=this%interp(lvl,r1,pG,gc)
                     pos_p=[r1(1)+0.5_WP*dxm,r1(2),r1(3)]; pos_m=[r1(1)-0.5_WP*dxm,r1(2),r1(3)]
                     n_out(1)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dxm
                     pos_p=[r1(1),r1(2)+0.5_WP*dym,r1(3)]; pos_m=[r1(1),r1(2)-0.5_WP*dym,r1(3)]
                     n_out(2)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dym
                     pos_p=[r1(1),r1(2),r1(3)+0.5_WP*dzm]; pos_m=[r1(1),r1(2),r1(3)-0.5_WP*dzm]
                     n_out(3)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dzm
                     nrm=norm2(n_out)+epsilon(1.0_WP); n_out=n_out/nrm
                     r2=r1-d_ib*n_out
                     call apply_col(k_w,eta_w,d_eff_w,m1,r2,v2)
                  end block ib_block
               end if

               ! Particle-particle collisions via cell-binned neighbor list
               if (associated(nbr_off)) then
                  do k=nbr_off(i1)+1,nbr_off(i1+1)
                     j_c=nbr_lst(k)+1
                     if (j_c.eq.i1) cycle
                     call apply_col(k_n,eta_n,this%contact_dist,0.5_WP*m1,p(j_c)%pos,p(j_c)%vel)
                  end do
               end if

               ! Accumulate as force/volume into F_bond (matches bond force units)
               p(i1)%F_bond=p(i1)%F_bond+f_local/this%dV
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do

   contains

      !> Soft-sphere normal force from virtual partner (r2_in, v2_in) onto i1.
      !> Mirrors amrlpt%col_force minus the rotational/friction terms (amrpd
      !> particles have no angular state). r1, v1, dt, this%clip_col are host-
      !> associated. f_local accumulates the force in [N].
      subroutine apply_col(kk,ee,d_eff,m_eff,r2_in,v2_in)
         real(WP), intent(in) :: kk,ee,d_eff,m_eff
         real(WP), dimension(3), intent(in) :: r2_in,v2_in
         real(WP) :: d12,rnv,r_influ,delta_n
         real(WP), dimension(3) :: n12,v12,f_n
         d12=norm2(r2_in-r1)
         if (d12.lt.10.0_WP*epsilon(d12)) return   ! self-overlap guard
         n12=(r2_in-r1)/d12
         v12=v1-v2_in
         rnv=dot_product(v12,n12)
         r_influ=min(abs(rnv)*dt,0.2_WP*d_eff)
         delta_n=min(d_eff+r_influ-d12,this%clip_col*d_eff)
         if (delta_n.le.0.0_WP) return
         f_n=(-m_eff*kk*delta_n-m_eff*ee*rnv)*n12
         f_local=f_local+f_n
      end subroutine apply_col

   end subroutine compute_contact

   !> Trilinear cell-centered interpolation of a multifab data array at a 3D
   !> position. Copied from amrlpt%interp -- self-contained, used by
   !> compute_contact for IB normal-distance evaluation.
   function interp(this,lvl,pos,arr,comp) result(val)
      implicit none
      class(amrpd), intent(in) :: this
      integer, intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: arr
      integer, intent(in) :: comp
      real(WP) :: val,wx,wy,wz
      integer  :: ii,jj,kk
      ii=floor((pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP); wx=(pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP-real(ii,WP)
      jj=floor((pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP); wy=(pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP-real(jj,WP)
      kk=floor((pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP); wz=(pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP-real(kk,WP)
      val=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*arr(ii  ,jj  ,kk  ,comp) &
      &  +        wx *(1.0_WP-wy)*(1.0_WP-wz)*arr(ii+1,jj  ,kk  ,comp) &
      &  +(1.0_WP-wx)*        wy *(1.0_WP-wz)*arr(ii  ,jj+1,kk  ,comp) &
      &  +        wx *        wy *(1.0_WP-wz)*arr(ii+1,jj+1,kk  ,comp) &
      &  +(1.0_WP-wx)*(1.0_WP-wy)*        wz *arr(ii  ,jj  ,kk+1,comp) &
      &  +        wx *(1.0_WP-wy)*        wz *arr(ii+1,jj  ,kk+1,comp) &
      &  +(1.0_WP-wx)*        wy *        wz *arr(ii  ,jj+1,kk+1,comp) &
      &  +        wx *        wy *        wz *arr(ii+1,jj+1,kk+1,comp)
   end function interp

   !> Velocity-Verlet step with full LPS bond pipeline + soft-contact walls.
   !>   1) Half-kick: vel += (dt/2) * (g + F_bond/rho + F_fluid/rho)   [PART_INTEGRATES]
   !>   2) Drift:     pos += dt * vel                                  [PART_MOVES]
   !>      WALL faces are NOT hard-reflected here -- soft contact in
   !>      compute_contact pushes particles back. OPEN faces let particles
   !>      drift out (Redistribute drops them).
   !>   3) Re-anchor bond.pos to lower-GID particle's current position so bonds
   !>      travel with their owning particle under the position-based redistribute.
   !>      Requires a rank-wide GID->position lookup since drift may move the
   !>      lower-GID across tiles on this rank.
   !>   4) Redistribute particles + bonds
   !>   5) fill_ghosts -> compute_dilatation -> update_ghosts -> compute_force
   !>      (damage accumulated at break time inside compute_force) ->
   !>      compute_contact (walls + IB + p-p)
   !>   6) Half-kick: vel += (dt/2) * acc_new                          [PART_INTEGRATES]
   subroutine advance(this,dt,Gib,Gibcomp)
      use amrdata_class, only: amrdata
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in), optional :: Gib       !< Optional IB level-set for contact (+ve inside fluid)
      integer, intent(in), optional :: Gibcomp         !< Component selector into Gib (defaults to 1)
      real(WP) :: rho_inv

      ! Cache invariants used by both half-kick blocks below
      rho_inv=1.0_WP/this%rho

      ! First half-kick (vel += dt/2 * acc) and drift (pos += dt * vel).
      ! WALL/OPEN BCs are NOT enforced here -- compute_contact handles WALL via
      ! soft-sphere penalty, and OPEN faces let particles drift out for
      ! Redistribute to drop. Loops over all AMR levels.
      first_halfkick_and_drift: block
         use amrex_amr_module, only: amrex_mfiter
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         integer(I8) :: np_,n
         integer :: lvl
         real(WP), dimension(3) :: acc
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  acc=this%gravity+(p(n)%F_bond+p(n)%F_fluid)*rho_inv
                  if (iand(p(n)%flag,PART_INTEGRATES).ne.0) then
                     p(n)%vel=p(n)%vel+0.5_WP*dt*acc
                  end if
                  if (iand(p(n)%flag,PART_MOVES).ne.0) then
                     p(n)%pos=p(n)%pos+dt*p(n)%vel
                  end if
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do
      end block first_halfkick_and_drift

      ! Re-anchor each bond's pos to the CURRENT position of its lower-GID
      ! endpoint, so that AMReX's position-based Redistribute keeps every bond
      ! co-located with its owning particle. The lower-GID particle is still on
      ! THIS rank (pre-redistribute invariant: bond and lower-GID share a rank),
      ! but drift may have moved it to a different tile -- so we build a rank-
      ! wide (GID -> position) lookup across all owned tiles before walking bonds.
      reanchor_bonds: block
         use amrex_amr_module, only: amrex_mfiter
         use amrpd_hash_class, only: gid_hash
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         type(bond), dimension(:), pointer :: b
         type(gid_hash) :: rhash
         integer(c_int64_t), allocatable :: gids(:)
         real(WP), allocatable :: poses(:,:)
         integer(I8) :: np_,nb_tile,n,ib,idx,np_owned
         integer :: lvl,lid
         integer(c_int64_t) :: key
         integer(c_int) :: parts(2)
         ! Count total owned particles on this rank
         np_owned=0_I8
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_)
               np_owned=np_owned+np_
            end do
            call this%mfiter_destroy(mfi)
         end do
         if (np_owned.gt.0_I8) then
            allocate(gids(np_owned),poses(3,np_owned))
            idx=0_I8
            do lvl=0,this%amr%clvl()
               call this%mfiter_build(lvl,mfi)
               do while (mfi%next())
                  call this%get_particles(lvl,mfi,p,np_)
                  do n=1_I8,np_
                     idx=idx+1_I8
                     gids(idx)=p(n)%idcpu
                     poses(:,idx)=p(n)%pos
                  end do
               end do
               call this%mfiter_destroy(mfi)
            end do
            call rhash%build(int(np_owned),gids)
            ! Walk bonds across all tiles, update pos to current lower-GID pos
            do lvl=0,this%amr%clvl()
               call this%mfiter_build(lvl,mfi)
               do while (mfi%next())
                  call this%get_bonds(lvl,mfi,b,nb_tile)
                  do ib=1_I8,nb_tile
                     if (b(ib)%alive.eq.0) cycle
                     parts(1)=b(ib)%id_lo_lo; parts(2)=b(ib)%id_lo_hi
                     key=transfer(parts,0_c_int64_t)
                     lid=rhash%lookup(key)
                     if (lid.lt.1) cycle   ! lower-GID not on this rank -- shouldn't happen
                     b(ib)%pos=poses(:,lid)
                  end do
               end do
               call this%mfiter_destroy(mfi)
            end do
            call rhash%finalize()
            deallocate(gids,poses)
         end if
      end block reanchor_bonds

      ! Settle particles AND bonds onto their new (rank,tile) -- handles periodic
      ! wrap and drops particles that left the domain through an OPEN face.
      ! Bonds follow their lower-GID particle because bond.pos was just re-anchored.
      call this%redistribute()

      ! Rebuild ghost layer (topology changed by redistribute), then recompute
      ! the bond state: dilatation first (each owner needs its own theta), then
      ! refresh ghost-slot values so compute_force can read neighboring owners'
      ! theta and m_w, then bond forces.
      call this%fill_ghosts(radius=this%search_radius)
      call this%compute_dilatation()
      call this%update_ghosts()
      call this%compute_force(dt)
      call this%compute_contact(dt=dt,Gib=Gib,Gibcomp=Gibcomp)
      call this%clear_ghosts()

      ! Refresh particle volume fraction on the Eulerian mesh so AMReX sees a
      ! fresh field if it considers a regrid before the next advance call.
      call this%update_VF()

      ! Second half-kick using freshly recomputed acceleration
      second_halfkick: block
         use amrex_amr_module, only: amrex_mfiter
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         integer(I8) :: np_,n
         integer :: lvl
         real(WP), dimension(3) :: acc
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  if (iand(p(n)%flag,PART_INTEGRATES).ne.0) then
                     acc=this%gravity+(p(n)%F_bond+p(n)%F_fluid)*rho_inv
                     p(n)%vel=p(n)%vel+0.5_WP*dt*acc
                  end if
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do
      end block second_halfkick
   end subroutine advance

   !> Compute the per-step CFL number from two BINDING sources, both using
   !> particle spacing dp = dV^(1/3) (the bond-network length scale, NOT the
   !> AMReX mesh dx -- waves propagate through bonds, not through cells):
   !>   - convective:   CFLp = max(|v_d|) * dt / dp    -- limit 0.1
   !>   - elastic:      CFLe = c_p        * dt / dp    -- limit 0.5
   !>
   !> Also reports a diagnostic contact CFL:
   !>   - contact:      CFLc = dt / tau_col            -- NOT binding
   !> tau_col is either user-set (this%tau_col > 0) or auto-derived as 5*dt;
   !> in the latter case CFLc is identically 0.2 (= 1/5).
   !>
   !> Convective limit is tighter than wave. To bind against the same 'Max CFL'
   !> threshold (typical 0.5), CFLp is scaled by 5x in the returned cfl. Raw
   !> values are stored on the type for monitor output.
   subroutine get_cfl(this,dt,cfl)
      implicit none
      class(amrpd), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      real(WP) :: K_bulk,mu_shear,c_p,dp_inv,tau_eff
      real(WP), parameter :: CFL_scale_conv=5.0_WP   !< convective scale (raw limit 0.1 -> 0.5)

      this%CFLp=0.0_WP; this%CFLe=0.0_WP; this%CFLc=0.0_WP

      ! P-wave speed and 1/dp (uniform particle spacing through bond network)
      K_bulk  =this%elastic_modulus/(3.0_WP*(1.0_WP-2.0_WP*this%poisson_ratio))
      mu_shear=this%elastic_modulus/(2.0_WP*(1.0_WP+this%poisson_ratio))
      c_p     =sqrt((K_bulk+4.0_WP*mu_shear/3.0_WP)/this%rho)
      dp_inv  =1.0_WP/this%dV**(1.0_WP/3.0_WP)

      ! Wave CFL is uniform (no per-particle dependence)
      this%CFLe=c_p*dp_inv

      ! Per-rank max convective CFL across all levels' owned particles
      local_max: block
         use amrex_amr_module, only: amrex_mfiter
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         integer(I8) :: np_,n
         integer :: lvl
         real(WP) :: vmax
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  vmax=max(abs(p(n)%vel(1)),abs(p(n)%vel(2)),abs(p(n)%vel(3)))
                  this%CFLp=max(this%CFLp,vmax*dp_inv)
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do
      end block local_max

      ! Global max across ranks, then scale by dt to get the actual CFL
      global_max: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: ierr
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp=this%CFLp*dt
         this%CFLe=this%CFLe*dt   ! uniform; no reduction needed
      end block global_max

      ! Contact CFL: diagnostic only, NOT included in returned cfl. tau_col is
      ! the user override if > 0, else the auto default = 5*dt (matches
      ! compute_contact). With the auto default, CFLc is identically 0.2.
      if (this%tau_col.gt.0.0_WP) then
         tau_eff=this%tau_col
      else
         tau_eff=5.0_WP*dt
      end if
      this%CFLc=dt/tau_eff

      ! Viscous (Maxwell) relaxation: dt/tau -- DIAGNOSTIC ONLY. The relaxation
      ! uses exact exponential integration (unconditionally stable), so it does
      ! NOT bind dt; reported for monitoring how resolved tau is.
      this%CFLv=0.0_WP
      if (this%tau.gt.0.0_WP.and.this%tau.lt.huge(1.0_WP)) this%CFLv=dt/this%tau
      ! Combine BINDING constraints only. Convective has tighter raw limit
      ! (0.1) than wave (0.5); scale by 5 to bind against the same Max CFL.
      cfl=max(CFL_scale_conv*this%CFLp,this%CFLe)
   end subroutine get_cfl



   ! ============================================================================
   ! CHECKPOINT I/O
   ! ============================================================================

   !> Write checkpoint for both containers under a shared checkpoint directory:
   !>   <dirname>/particles/   (AMReX particle Checkpoint)
   !>   <dirname>/bonds/       (AMReX bond Checkpoint)
   !> Caller is responsible for creating <dirname> (typically via io%write under
   !> the same directory).
   subroutine write(this,dirname)
      implicit none
      class(amrpd), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrpd_checkpoint_p(this%pcp,trim(dirname)//'/particles'//c_null_char)
      call amrpd_checkpoint_b(this%pcb,trim(dirname)//'/bonds'//c_null_char)
   end subroutine write

   !> Restore both containers from a checkpoint directory written by write.
   !> The amrgrid must already have been rebuilt via amr%init_from_checkpoint
   !> before calling this. Syncs both containers' BA/DM to the restored grid
   !> before AMReX Restart so particles land on the right ranks, then
   !> redistributes and refreshes counters.
   subroutine read(this,dirname)
      implicit none
      class(amrpd), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      integer :: lvl
      ! Sync both containers to the restored grid's BA/DM
      do lvl=0,this%amr%clvl()
         call this%set_particle_ba_p(lvl,this%amr%get_boxarray(lvl))
         call this%set_particle_dm_p(lvl,this%amr%get_distromap(lvl))
         call this%set_particle_ba_b(lvl,this%amr%get_boxarray(lvl))
         call this%set_particle_dm_b(lvl,this%amr%get_distromap(lvl))
      end do
      ! AMReX Restart on each container
      call amrpd_restart_p(this%pcp,trim(dirname)//'/particles'//c_null_char)
      call amrpd_restart_b(this%pcb,trim(dirname)//'/bonds'//c_null_char)
      ! Settle and refresh
      call this%redistribute()
      call this%get_info()
   end subroutine read


   ! ============================================================================
   ! DIAGNOSTICS
   ! ============================================================================

   !> Log solver info on root
   subroutine print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(amrpd), intent(inout) :: this
      character(len=str_long) :: message
      if (this%amr%amRoot) then
         write(message,'("AMRPD solver [",a,"] on AMR grid [",a,"]")') trim(this%name),trim(this%amr%name)
         if (verbose .gt. 1) write(output_unit,'(a)') trim(message)
         if (verbose .gt. 0) call log(message)
      end if
   end subroutine print


end module amrpd_class
