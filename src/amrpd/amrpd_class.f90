!> Grid-side extension of the PD solver: an amrpd IS a pdsolver (extends it)
!> plus an AMReX particle mirror of its nodes (the "face"), distributed by
!> position over the fluid grid's boxes so every grid-facing operation runs
!> where the cells live: volume-fraction and velocity deposits, field
!> interpolation at particle positions (F_fluid), VF-driven AMR tagging,
!> and plotfile visualization.
!>
!> Usage tiers (amrpd EXTENDS pdsolver: one object is the solid solver AND
!> its grid face -- assign material/damage/contact fields, then call handoff):
!>   1. pdsolver alone      -- grid-free solid dynamics (no visualization)
!>   2. amrpd               -- adds viz, solid VF on the mesh, AMR refinement,
!>                             seeding, and the face<->solver exchange
!>   3. ... + a flow solver -- two-way FSI (driver deposits IB forcing; the
!>                             fluid load returns via exchange_solid)
module amrpd_class
   use precision,      only: WP,I8
   use string,         only: str_medium
   use amrgrid_class,  only: amrgrid
   use amrdata_class,  only: amrdata
   use pdsolver_class, only: pdsolver,pd_partition,PD_WALL,PDC_IS_DEAD
   use iso_c_binding
   implicit none
   private

   ! Public exports
   public :: amrpd,part,part_gid

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

   ! Struct layout constants -- MUST match #defines in amrpd_wrapper.cpp
   integer, parameter, public :: AMRPD_NREAL_PART = 15
   integer, parameter, public :: AMRPD_NINT_PART  = 1
   integer, parameter, public :: AMRPD_NREAL_BOND = 4
   integer, parameter, public :: AMRPD_NINT_BOND  = 5

   !> Solid particle struct -- must match C++ Particle<15,1> memory layout:
   !> pos[3], rdata[15], idcpu, idata[1]. Physics lives in pdsolver; only
   !> pos/vel/F_fluid/damage are meaningful here (rest is legacy layout).
   type, bind(C), public :: part
      real(c_double) :: pos(3)                  !< AMReX-managed position
      real(c_double) :: vel(3)                  !< rdata[0..2]
      real(c_double) :: F_bond(3)               !< rdata[3..5]  (unused)
      real(c_double) :: F_fluid(3)              !< rdata[6..8]  fluid load interpolated by the driver
      real(c_double) :: mw                      !< rdata[9]     (unused)
      real(c_double) :: dil                     !< rdata[10]    (unused)
      real(c_double) :: damage                  !< rdata[11]    broken-bond fraction in [0,1]
      real(c_double) :: nb0                     !< rdata[12]    (unused)
      real(c_double) :: td2                     !< rdata[13]    (unused)
      real(c_double) :: td2a                    !< rdata[14]    (unused)
      integer(c_int64_t), private :: idcpu      !< AMReX packed id+cpu
      integer(c_int) :: flag                    !< idata[0]: PART_* flags at seeding; owner routing tag (7+8*owner) after handoff
   end type part



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

      ! Redistribute
      subroutine amrpd_redistribute_p(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine


      ! MFIter accessors -- particles
      subroutine amrpd_get_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine


      ! Single-element insertion (initialization)
      subroutine amrpd_add_particle_i(pc,lev,grid,tile,p) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,p
         integer(c_int), value :: lev,grid,tile
      end subroutine

      ! Bulk append at level 0 (collective; ranks with n=0 pass raw=NULL)
      subroutine amrpd_append_particles(pc,raw,n) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         type(c_ptr), value :: raw
         integer(c_int64_t), value :: n
      end subroutine
      subroutine amrpd_append_particles_gid(pc,raw,n,gids) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         type(c_ptr), value :: raw
         integer(c_int64_t), value :: n
         type(c_ptr), value :: gids
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

      ! ID/CPU counters and accessors
      subroutine amrpd_get_next_id_p(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine
      subroutine amrpd_set_next_id_p(id) bind(c)
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


      ! Global counts
      subroutine amrpd_total_np(pc,np) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: np
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

   end interface


   !> Grid-side extension of the PD solver: an amrpd IS a pdsolver, plus an
   !> AMReX particle mirror of its nodes (the "face") distributed by position
   !> over the fluid grid's boxes for deposits, interpolation, tagging, and viz
   type, extends(pdsolver) :: amrpd

      !> Associated AMR grid
      class(amrgrid), pointer :: amr => null()

      !> Opaque AMReX container handle
      type(c_ptr) :: pcp = c_null_ptr           !< Particle container

      !> Face load-balance metrics across ranks (grid decomposition; the
      !> solver's own Morton partition is balanced by construction)
      integer(I8) :: np_loc = 0                 !< This rank's face particle count
      integer(I8) :: np_min = 0                 !< Min across ranks
      integer(I8) :: np_max = 0                 !< Max across ranks
      real(WP)    :: np_eff = 0.0_WP            !< Load efficiency = mean/max

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


      !> Particle volume fraction on the Eulerian AMR mesh. Cell-centered scalar
      !> (one component), one ghost layer. Updated each advance step by
      !> update_VF: trilinear deposition of each particle's volume dV onto the
      !> 8 surrounding cell centers, then average-down + optional smoothing.
      !> Drives the AMR tagging callback when VF_tag > 0.
      type(amrdata) :: VF
      real(WP)     :: VF_tag       = -1.0_WP    !< Refinement threshold (<=0 disables VF-driven tagging)
      real(WP)     :: filter_width =  0.0_WP    !< Gaussian-equivalent filter width for VF; 0 disables
      real(WP)     :: VF_snap      =  0.1_WP    !< Deposit-moire amplitude: VF is rescaled by 1/(1-VF_snap) and clipped,
                                                !< so a fully packed interior reads exactly 1. 0 disables.
      real(WP)     :: VFmin=0.0_WP,VFmax=0.0_WP,VFmean=0.0_WP   !< VF statistics

      !> Optional user-supplied tagging callback. Called AFTER the built-in VF
      !> tagging. Use to add custom refinement criteria (e.g., damage > 0.3).
      procedure(pd_tagging_iface), pointer, pass :: user_pd_tagging => null()

   contains
      ! Lifecycle
      procedure :: initialize
      procedure :: finalize
      ! Container utilities
      procedure :: redistribute
      procedure :: get_info             !< Global counts + min/max/mean velocities + load-balance metrics
      procedure :: set_particle_ba_p
      procedure :: set_particle_dm_p
      ! Particle population
      procedure :: append
      procedure :: append_with_gids
      ! Solver coupling (grid face <-> contained pdsolver)
      procedure :: handoff
      procedure :: exchange_solid
      procedure :: rebuild_face
      ! MFIter helpers (particle container's BA/DM)
      procedure :: mfiter_build
      procedure :: mfiter_destroy
      procedure :: get_particles
      ! AMR callbacks
      procedure :: post_regrid
      procedure :: tagging
      ! Particle volume fraction + AMR tagging
      procedure :: update_VF                         !< Compute VF from particle positions (trilinear deposit)
      procedure :: process_deposit                   !< Post-process a deposited field (extensive -> intensive + C/F transfers; public: also used on driver-deposited fields)
      procedure :: filter                            !< Explicit-diffusion smoothing of a cell-centered amrdata (public: also used on driver-deposited fields)
      ! Physics -- STUBBED in skeleton
      procedure :: interp                !< Trilinear cell-centered interpolation (used by compute_contact for IB)
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


   !> Public accessor for a particle's unique 64-bit GID key (its AMReX idcpu).
   !> The idcpu component is private to protect the AMReX layout; the graph-core
   !> handoff and parity tooling need the key for gid-matched state exchange.
   function part_gid(p) result(gid)
      implicit none
      type(part), intent(in) :: p
      integer(I8) :: gid
      gid=p%idcpu
   end function part_gid


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

   !> Initialize amrpd solver: create particle container, register AMR callbacks
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
      ! Default deposit-smoothing width
      this%filter_width = 2.0_WP*this%amr%min_meshsize(this%amr%maxlvl)
      ! Create AMReX particle container
      call amrpd_new_pcp(this%pcp,this%amr%amrcore)
      ! Stamp the solver's domain geometry from the grid (material, damage,
      ! and contact fields are driver-assigned)
      this%Ldom=[amr%xhi-amr%xlo,amr%yhi-amr%ylo,amr%zhi-amr%zlo]
      this%per=[amr%xper,amr%yper,amr%zper]
      this%collapsed=[amr%nx.eq.1,amr%ny.eq.1,amr%nz.eq.1]
      this%dom_lo=[amr%xlo,amr%ylo,amr%zlo]
      this%dom_hi=[amr%xhi,amr%yhi,amr%zhi]
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

   !> Finalize: destroy face and container, then the parent solver
   subroutine finalize(this)
      implicit none
      class(amrpd), intent(inout) :: this
      ! Drop user tagging hook
      nullify(this%user_pd_tagging)
      ! Tear down the volume-fraction field
      call this%VF%finalize()
      if (c_associated(this%pcp)) then
         call amrpd_delete_pcp(this%pcp); this%pcp = c_null_ptr
      end if
      nullify(this%amr)
      ! Tear down the solver state
      call this%pdsolver%finalize()
   end subroutine finalize


   ! ============================================================================
   ! CONTAINER UTILITIES
   ! ============================================================================

   !> Redistribute the particle container across ranks (finest covering level
   !> per position; lev_max=-1 avoids the lev_max>finestLevel() assert before
   !> all levels exist). Matches amrlpt's default.
   subroutine redistribute(this)
      implicit none
      class(amrpd), intent(inout) :: this
      call amrpd_redistribute_p(this%pcp,0,-1,0)
   end subroutine redistribute












   !> Collective info: runs the solver's get_info (np, velocity extrema, bond
   !> censuses, timers), then adds the FACE load-balance metrics -- particles
   !> per rank across the grid decomposition (the solver's own partition is
   !> balanced by construction; this measures the position-based mirror).
   subroutine get_info(this)
      implicit none
      class(amrpd), intent(inout) :: this
      integer(I8) :: np_sum

      ! Solver-side info
      call this%pdsolver%get_info()

      ! Face census: live particles owned by this rank across all levels
      local_pass: block
         use amrex_amr_module, only: amrex_mfiter
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         integer(I8) :: np_,n
         integer :: lvl
         this%np_loc=0_I8
         do lvl=0,this%amr%clvl()
            call this%mfiter_build(lvl,mfi)
            do while (mfi%next())
               call this%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  this%np_loc=this%np_loc+1_I8
               end do
            end do
            call this%mfiter_destroy(mfi)
         end do
      end block local_pass

      ! Load-balance reductions
      global_reduce: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8
         integer :: ierr
         this%np_min=this%np_loc; this%np_max=this%np_loc; np_sum=this%np_loc; this%np_eff=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_min,1,MPI_INTEGER8,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_max,1,MPI_INTEGER8,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,np_sum,     1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
         if (this%np_max.gt.0_I8) this%np_eff=real(np_sum,WP)/real(this%np_max,WP)/real(this%amr%nproc,WP)
      end block global_reduce

   end subroutine get_info


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

   !> Bulk-append with PRESERVED identities: each particle takes the (id,cpu)
   !> packed in gids (the part_gid key). Used to rebuild the grid-side face
   !> from a pdsolver checkpoint so identities match the solver's node gids.
   !> Collective; AddParticlesAtLevel redistributes by position.
   subroutine append_with_gids(this,plist,n,gids)
      use messager, only: die
      implicit none
      class(amrpd), intent(inout) :: this
      type(part), dimension(:), allocatable, target, intent(in) :: plist
      integer(I8), intent(in) :: n
      integer(I8), dimension(:), target, intent(in) :: gids
      type(c_ptr) :: raw,graw
      raw=c_null_ptr; graw=c_null_ptr
      if (n.gt.0_I8.and.allocated(plist)) then
         if (int(size(plist),I8).lt.n) call die('[amrpd append_with_gids] plist smaller than n')
         raw=c_loc(plist(1)); graw=c_loc(gids(1))
      end if
      call amrpd_append_particles_gid(this%pcp,raw,int(n,c_int64_t),graw)
   end subroutine append_with_gids


   ! ============================================================================
   ! SOLVER COUPLING (grid face <-> contained pdsolver)
   ! ============================================================================

   !> Hand the seeded face population to the solver: extract nodes, Morton-
   !> partition them (balanced, motion-invariant), detect families, and stamp
   !> each face particle's flag with its solver owner rank (flag = 7+8*owner).
   !> The driver must have assigned the material/damage/contact fields first.
   subroutine handoff(this)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8), allocatable :: gids(:),rgid(:)
      real(WP), allocatable :: pos(:,:),vel(:,:),voll(:),rpos(:,:),rvel(:,:),rvol(:)
      integer, allocatable :: flags(:),owner(:),rflag(:)
      integer(I8) :: np_,n
      integer :: lvl,nn,i,nr
      ! Extract this rank's owned particles
      nn=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            nn=nn+int(np_)
         end do
         call this%mfiter_destroy(mfi)
      end do
      allocate(gids(max(nn,1)),pos(3,max(nn,1)),vel(3,max(nn,1)),flags(max(nn,1)),voll(max(nn,1)),owner(max(nn,1)))
      i=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               i=i+1
               gids(i) =part_gid(p(n))
               pos(:,i)=p(n)%pos
               vel(:,i)=p(n)%vel
               flags(i)=p(n)%flag
               voll(i) =this%dV
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Balanced static partition of the reference configuration
      call pd_partition(nn,gids,pos,vel,flags,voll,owner,nr,rgid,rpos,rvel,rflag,rvol)
      call this%set_nodes(nr,rgid,rpos,rvel,rflag,rvol)
      call this%detect_families()
      ! Stamp owner routing tags on the grid face (same walk order as the
      ! extraction above, so owner(i) lines up)
      i=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               i=i+1
               p(n)%flag=7+8*owner(i)
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      deallocate(gids,pos,vel,flags,voll,owner,rgid,rpos,rvel,rflag,rvol)
   end subroutine handoff

   !> Face <-> solver exchange (collective; call once per coupling step): push
   !> each face particle's interpolated F_fluid to its solver owner; pull back
   !> the owner's current (pos, vel, damage, alive). Dead solver nodes (exited
   !> an open face) get an outside-domain position written back, so the next
   !> redistribute drops the tombstone.
   subroutine exchange_solid(this)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrpd), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8), allocatable :: mgid(:)
      integer, allocatable :: mown(:)
      real(WP), allocatable :: mff(:,:),mpos(:,:),mvel(:,:),mdmg(:),malive(:)
      integer(I8) :: np_,n
      integer :: lvl,nm,i
      ! Count live face particles
      nm=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               nm=nm+1
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      allocate(mgid(max(nm,1)),mown(max(nm,1)),mff(3,max(nm,1)))
      allocate(mpos(3,max(nm,1)),mvel(3,max(nm,1)),mdmg(max(nm,1)),malive(max(nm,1)))
      ! Pack (gid, owner tag, F_fluid)
      i=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               i=i+1
               mgid(i)=part_gid(p(n))
               mown(i)=p(n)%flag/8
               mff(:,i)=p(n)%F_fluid
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Collective round-trip with the solver
      call this%exchange(nm,mgid,mown,mff,mpos,mvel,mdmg,malive)
      ! Write the solver state back onto the grid face (same walk order)
      i=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               i=i+1
               p(n)%pos=mpos(:,i)
               p(n)%vel=mvel(:,i)
               p(n)%damage=mdmg(i)
               if (malive(i).lt.0.5_WP) p(n)%flag=PART_IS_DEAD
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      deallocate(mgid,mown,mff,mpos,mvel,mdmg,malive)
   end subroutine exchange_solid

   !> Rebuild the grid face from the (restored) solver state: one particle per
   !> owned live node with PRESERVED identity (gid -> id,cpu), current
   !> position/velocity/damage, and the owner routing tag = this rank.
   subroutine rebuild_face(this)
      use parallel, only: rank
      implicit none
      class(amrpd), intent(inout) :: this
      type(part), dimension(:), allocatable, target :: plist
      integer(I8), dimension(:), allocatable, target :: gl
      integer(I8) :: n
      integer :: i,m
      allocate(plist(max(this%nown,1)),gl(max(this%nown,1)))
      m=0
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle   ! exited nodes stay dead
         m=m+1
         plist(m)%pos    =this%y(:,i)
         plist(m)%vel    =this%v(:,i)
         plist(m)%F_bond =0.0_WP
         plist(m)%F_fluid=0.0_WP
         plist(m)%mw     =0.0_WP
         plist(m)%dil    =0.0_WP
         plist(m)%damage =this%damage(i)
         plist(m)%nb0    =0.0_WP
         plist(m)%td2    =0.0_WP
         plist(m)%td2a   =0.0_WP
         plist(m)%flag   =7+8*rank
         gl(m)=this%gid(i)
      end do
      n=int(m,I8)
      call this%append_with_gids(plist,n,gl)
      call this%redistribute()
      deallocate(plist,gl)
   end subroutine rebuild_face


   ! ============================================================================
   ! AMR CALLBACKS
   ! ============================================================================

   !> Post-regrid: re-sync particle and bond containers' BA/DM, redistribute.
   !> Particles span levels
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
      ! and matches the pattern used in amrlpt.
      sync_to_amrgrid: block
         do lvl=0,this%amr%clvl()
            call this%set_particle_ba_p(lvl,this%amr%get_boxarray(lvl))
            call this%set_particle_dm_p(lvl,this%amr%get_distromap(lvl))
         end do
      end block sync_to_amrgrid

      ! First settle particles onto the freshly-synced BA/DM before we can
      ! count them per box for the knapsack weighting below
      call this%redistribute()

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
               if (this%lo_bc(1).eq.PD_WALL.and.ii  .lt.this%amr%geom(lvl)%domain%lo(1)) then; ii=this%amr%geom(lvl)%domain%lo(1)  ; wx=0.0_WP; end if
               if (this%hi_bc(1).eq.PD_WALL.and.ii+1.gt.this%amr%geom(lvl)%domain%hi(1)) then; ii=this%amr%geom(lvl)%domain%hi(1)-1; wx=1.0_WP; end if
               if (this%lo_bc(2).eq.PD_WALL.and.jj  .lt.this%amr%geom(lvl)%domain%lo(2)) then; jj=this%amr%geom(lvl)%domain%lo(2)  ; wy=0.0_WP; end if
               if (this%hi_bc(2).eq.PD_WALL.and.jj+1.gt.this%amr%geom(lvl)%domain%hi(2)) then; jj=this%amr%geom(lvl)%domain%hi(2)-1; wy=1.0_WP; end if
               if (this%lo_bc(3).eq.PD_WALL.and.kk  .lt.this%amr%geom(lvl)%domain%lo(3)) then; kk=this%amr%geom(lvl)%domain%lo(3)  ; wz=0.0_WP; end if
               if (this%hi_bc(3).eq.PD_WALL.and.kk+1.gt.this%amr%geom(lvl)%domain%hi(3)) then; kk=this%amr%geom(lvl)%domain%hi(3)-1; wz=1.0_WP; end if
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
      ! Snap out the deposit moire: the particle lattice is incommensurate with the
      ! grid (and moves), so a fully packed interior deposits VF slightly below 1.
      ! Rescale+clip so it reads exactly 1; continuous, so no jump is introduced.
      if (this%VF_snap.gt.0.0_WP) then
         do lvl=0,this%amr%clvl()
            call this%VF%mf(lvl)%mult(1.0_WP/(1.0_WP-this%VF_snap),1,1,this%VF%ng)
         end do
         call this%VF%clip(0.0_WP,1.0_WP)
         call this%VF%fill(time=0.0_WP)
      end if

   contains

      !> Helper: get the particle distribution map for this level. Returns the
      !> DM that the particle container is currently using (which may differ
      !> from the Eulerian DM).
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
   ! FIELD INTERPOLATION
   ! ============================================================================



   !> Trilinear cell-centered interpolation of a multifab data array at a 3D
   !> position (drivers use it to sample fields at particle positions).
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





   ! ============================================================================
   ! CHECKPOINT I/O
   ! ============================================================================

   !> Write the particle-container checkpoint under <dirname>/particles.
   !> Caller is responsible for creating <dirname>.
   subroutine write(this,dirname)
      implicit none
      class(amrpd), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrpd_checkpoint_p(this%pcp,trim(dirname)//'/particles'//c_null_char)
   end subroutine write

   !> Restore the particle container from a checkpoint written by write. The
   !> amrgrid must already be rebuilt; syncs BA/DM, restarts, redistributes.
   subroutine read(this,dirname)
      implicit none
      class(amrpd), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      integer :: lvl
      ! Sync the container to the restored grid's BA/DM
      do lvl=0,this%amr%clvl()
         call this%set_particle_ba_p(lvl,this%amr%get_boxarray(lvl))
         call this%set_particle_dm_p(lvl,this%amr%get_distromap(lvl))
      end do
      ! AMReX Restart
      call amrpd_restart_p(this%pcp,trim(dirname)//'/particles'//c_null_char)
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
