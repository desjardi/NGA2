!> AMR-aware Lagrangian particle tracking solver
!> Mirrors lpt_class capabilities on an AMReX AMR hierarchy.
!> Particle communication and sorting handled by NeighborParticleContainer<14,1>.
module amrlpt_class
   use precision,     only: WP,I8
   use string,        only: str_medium
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use timer_class,   only: timer
   use iso_c_binding
   implicit none
   private

   ! Public exports
   public :: amrlpt,part

   ! Boundary condition flags for collision detection
   integer, parameter, public :: AMRLPT_OPEN=0  !< No wall collision (open or periodic boundary)
   integer, parameter, public :: AMRLPT_WALL=1  !< Hard wall collision at domain boundary

   ! Particle struct layout constants (must match #define in amrlpt_wrapper.cpp)
   integer, parameter, public :: AMRLPT_NREAL=14  !< extra reals per particle
   integer, parameter, public :: AMRLPT_NINT =1   !< extra ints  per particle

   ! Particle flags
   integer(c_int), parameter, public :: PART_MOVES    =1
   integer(c_int), parameter, public :: PART_COLLIDES =2
   integer(c_int), parameter, public :: PART_EXCHANGES=4
   integer(c_int), parameter, public :: PART_IS_DEAD  =0

   ! Particle struct -- must match C++ Particle<14,1> memory layout exactly:
   ! - pos[3]    (pos, managed by AMReX)
   ! - rdata[14] (d, vel[3], angVel[3], Acol[3], Tcol[3], dt)
   ! - idcpu     (packed id+cpu, private)
   ! - idata[1]  (flag)
   type, bind(C), public :: part
      !> AMReX position (physical coordinates)
      real(c_double) :: pos(3)
      !> Extra reals (NStructReal=14)
      real(c_double) :: d               !< Particle diameter
      real(c_double) :: vel(3)          !< Particle velocity
      real(c_double) :: angVel(3)       !< Angular velocity
      real(c_double) :: Acol(3)         !< Collision acceleration
      real(c_double) :: Tcol(3)         !< Collision torque
      real(c_double) :: dt              !< Particle sub-timestep
      !> Packed id (39 bits, sign=valid) + cpu (24 bits)  [AMReX internal]
      integer(c_int64_t), private :: idcpu
      !> Extra int (NStructInt=1)
      integer(c_int) :: flag            !< 0=active, 1=mark for removal
   end type part

   !> C interface: bind(C) declarations to amrlpt_wrapper.cpp
   interface

      subroutine amrlpt_new_pc(pc, amrcore) bind(c)
         import :: c_ptr
         type(c_ptr) :: pc
         type(c_ptr), value :: amrcore
      end subroutine

      subroutine amrlpt_delete_pc(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_redistribute(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine

      subroutine amrlpt_fill_neighbors(pc,ngrow) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: ngrow
      end subroutine

      subroutine amrlpt_fill_neighbors_radius(pc,search_radius) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: search_radius
      end subroutine

      subroutine amrlpt_clear_neighbors(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_get_neighbor_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_build_neighbor_list(pc,rcrit) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: rcrit
      end subroutine

      subroutine amrlpt_get_neighbor_list_mfi(pc,lev,mfi,offsets,list,np,ntot) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: offsets
         type(c_ptr) :: list
         integer(c_int64_t) :: np
         integer(c_int64_t) :: ntot
      end subroutine

      subroutine amrlpt_get_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_get_all_particles_mfi(pc,lev,mfi,dp,np_total,np_valid) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np_total
         integer(c_int64_t) :: np_valid
      end subroutine

      subroutine amrlpt_add_particle_i(pc,lev,grid,tile,p) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,p
         integer(c_int), value :: lev,grid,tile
      end subroutine

      subroutine amrlpt_get_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine

      subroutine amrlpt_set_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t), value :: id
      end subroutine

      subroutine amrlpt_get_cpu(cpu) bind(c)
         import :: c_int
         integer(c_int) :: cpu
      end subroutine

      subroutine amrlpt_set_particle_id(id,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t), value :: id
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_set_particle_cpu(cpu,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int), value :: cpu
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_particle_is_valid(valid,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int) :: valid
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_total_np(pc,np) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_write(pc,path,is_chk) bind(c)
         import :: c_ptr,c_char,c_int
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
         integer(c_int), value :: is_chk
      end subroutine

      subroutine amrlpt_read(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
      end subroutine

      subroutine amrlpt_append_particles(pc,raw,n) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         type(c_ptr), value :: raw
         integer(c_int64_t), value :: n
      end subroutine

      subroutine amrlpt_get_particle_boxarray(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: ba
      end subroutine

      subroutine amrlpt_get_particle_distromap(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: dm
      end subroutine

      subroutine amrlpt_set_particle_boxarray(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: ba
      end subroutine

      subroutine amrlpt_set_particle_distromap(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: dm
      end subroutine

   end interface

   !> AMR LPT solver type
   type :: amrlpt

      !> Associated AMR grid
      class(amrgrid), pointer :: amr=>null()

      !> AMReX NeighborParticleContainer<14,1> opaque handle
      type(c_ptr) :: pc=c_null_ptr

      !> Solver name
      character(len=str_medium) :: name='UNNAMED_AMRLPT'

      !> User-provided injection callback
      procedure(lpt_inject_iface ), pointer, pass :: inject=>null()
      procedure(lpt_tagging_iface), pointer, pass :: user_lpt_tagging=>null()

      !> Global particle count
      integer(I8) :: np=0

      !> Physics parameters
      real(WP) :: rho                                   !< Particle material density
      real(WP), dimension(3) :: gravity=0.0_WP          !< Acceleration of gravity
      integer :: nstep=1                                !< Substeps per timestep
      character(len=str_medium) :: drag_model='Tenneti' !< Drag model

      !> Collision parameters
      real(WP) :: tau_col=1.0e-3_WP             !< Collision time scale
      real(WP) :: e_n=1.0_WP                    !< Normal restitution
      real(WP) :: e_w=1.0_WP                    !< Wall restitution
      real(WP) :: mu_f=0.0_WP                   !< Friction coefficient
      real(WP) :: clip_col=0.2_WP               !< Max overlap fraction

      !> CFL numbers
      real(WP) :: CFLp_x=0.0_WP,CFLp_y=0.0_WP,CFLp_z=0.0_WP
      real(WP) :: CFL_col=0.0_WP

      !> Monitoring data
      real(WP) :: dmin,dmax,dmean,dvar
      real(WP) :: Umin,Umax,Umean,Uvar
      real(WP) :: Vmin,Vmax,Vmean,Vvar
      real(WP) :: Wmin,Wmax,Wmean,Wvar
      real(WP) :: VFmin,VFmax,VFmean,VFvar
      integer  :: np_new=0,np_out=0
      real(WP) :: Vp_new=0.0_WP,Vp_out=0.0_WP,Vp_tot=0.0_WP
      integer  :: np_new_loc=0,np_out_loc=0
      real(WP) :: Vp_new_loc=0.0_WP,Vp_out_loc=0.0_WP
      integer  :: ncol=0

      !> Particle distribution across ranks
      integer(I8) :: np_loc=0           !< This rank's particle count
      integer(I8) :: np_min=0           !< Min particle count across ranks
      integer(I8) :: np_max=0           !< Max particle count across ranks
      real(WP)    :: np_eff=0.0_WP      !< Load efficiency (mean/max)

      !> Overlap size
      integer :: nover=2

      !> Filter width
      real(WP) :: filter_width=0.0_WP

      !> Particle volume fraction
      type(amrdata) :: VF

      !> VF threshold for tagging; disabled if <=0
      real(WP) :: VF_tag=-1.0_WP

      !> Two-way coupling source terms
      type(amrdata) :: src

      !> Boundary condition flags for collision detection
      integer, dimension(3) :: lo_bc=AMRLPT_OPEN
      integer, dimension(3) :: hi_bc=AMRLPT_OPEN

      !> Particle sub-stepping state
      real(WP) :: t=0.0_WP            !< Current particle time
      real(WP) :: dt=huge(1.0_WP)     !< Persistent CFL-limited sub-step size
      real(WP) :: dtmax=huge(1.0_WP)  !< Maximum allowed particle sub-step
      real(WP) :: cflmax=huge(1.0_WP) !< CFL limit for particle sub-stepping

      !> Knapsack load balancing
      logical :: rebalance=.false.    !< Enable knapsack DM for particles

      !> Per-phase wall-clock timers
      type(timer) :: tmr_coll         !< Collision phase total
      type(timer) :: tmr_coll_        !< Collision phase w/o MPI
      type(timer) :: tmr_fill         !< Ghost fill time
      type(timer) :: tmr_nbl          !< Neighbor list build time
      type(timer) :: tmr_step         !< Particle stepping phase
      type(timer) :: tmr_step_        !< Particle stepping phase w/o MPI
      type(timer) :: tmr_vf           !< Volume fraction computation
      type(timer) :: tmr_src          !< Source post-processing

   contains
      ! Lifecycle callbacks
      procedure :: post_regrid            !< Post-regrid callback
      procedure :: tagging                !< Tag cells for refinement based on VF
      procedure :: get_cost               !< Compute per-box particle cost for load balancing
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Physics procedures
      procedure :: collide                !< Soft-sphere collision model
      procedure :: advance                !< Advance particle ODEs one timestep
      procedure :: advance_to             !< CFL-substepped advance to a target time
      procedure, private :: step          !< Single sub-step: MFIter loop + redistribute
      procedure :: update_VF              !< Compute particle volume fraction field
      procedure :: get_cfl                !< Compute particle CFL numbers
      ! Utilities
      procedure :: redistribute           !< Call AMReX redistribute
      procedure :: fill_ghosts            !< Fill ghost particle buffer
      procedure :: fill_ghosts_radius     !< Fill ghosts with position-based filtering
      procedure :: clear_ghosts           !< Release ghost particle buffer
      procedure :: build_neighbor_list    !< Build explicit pair list within rcrit
      procedure :: get_np                 !< Update global particle count
      procedure, private :: get_particles     !< Get valid particle array for MFIter tile
      procedure, private :: get_ghosts        !< Get ghost particle array for MFIter tile
      procedure, private :: get_all_particles !< Get combined valid+ghost array for MFIter tile
      procedure, private :: get_neighbor_list !< Get CSR neighbor list for MFIter tile
      procedure :: interp                 !< Trilinear cell-centered interpolation
      procedure :: interp_face_velocities !< Trilinear face-centered interpolation
      procedure, private :: process_deposit   !< Post-process deposit: intensive conversion + F↔C transfers
      procedure, private :: filter        !< Explicit diffusion filter
      procedure :: gather_region          !< Allgather particles within a bounding box
      procedure :: append                 !< Append a Fortran part array
      procedure :: mfiter_build           !< Build MFIter for particle grid
      procedure :: mfiter_destroy         !< Destroy MFIter
      procedure :: set_particle_ba        !< Set particle BoxArray for a level
      procedure :: set_particle_dm        !< Set particle DistributionMapping for a level
      procedure :: get_particle_dm        !< Get particle DistributionMapping for a level
      ! Print solver info
      procedure :: get_info
      procedure :: print
      ! Checkpoint I/O
      procedure :: read
      procedure :: write
   end type amrlpt

   !> Abstract interface for user-provided particle injection callback
   abstract interface
      subroutine lpt_inject_iface(this,dt)
         import :: amrlpt,WP
         class(amrlpt), intent(inout) :: this
         real(WP), intent(in) :: dt
      end subroutine lpt_inject_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine lpt_tagging_iface(solver,lvl,time,tags)
         import :: amrlpt,c_ptr,WP
         class(amrlpt), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine lpt_tagging_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrlpt type
   ! ============================================================================

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrlpt_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrlpt), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrlpt_postregrid

   !> Dispatch tagging: calls type-bound method then user callback
   subroutine amrlpt_tagging(ctx,lvl,time,tags)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrlpt), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_lpt_tagging)) call this%user_lpt_tagging(lvl,time,tags)
   end subroutine amrlpt_tagging

   !> Dispatch get_cost: calls type-bound method
   subroutine amrlpt_get_cost(ctx,lvl,nboxes,costs,ba)
      use amrex_amr_module, only: amrex_boxarray
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrlpt), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrlpt_get_cost

   ! ============================================================================
   ! LIFECYCLE CALLBACKS
   ! ============================================================================

   !> Post-regrid callback: redistribute particles
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time

      ! Force resync particle container with new ba/dm
      resync: block
         integer :: lvl
         do lvl=0,this%amr%clvl()
            call this%set_particle_ba(lvl,this%amr%get_boxarray(lvl))
            call this%set_particle_dm(lvl,this%amr%get_distromap(lvl))
         end do
      end block resync

      ! Redistribute particles onto new grid
      call this%redistribute()

      ! Optionally rebalance particles with knapsack
      if (this%rebalance) then
         rebalance_particles: block
            use amrex_amr_module, only: amrex_distromap,amrex_distromap_destroy,amrex_mfiter,amrex_boxarray
            use amrex_interface,  only: amrdm_make_knapsack
            use parallel, only: MPI_REAL_WP
            use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            integer :: lvl,nboxes,ierr
            real(WP), dimension(:), allocatable :: costs
            type(amrex_boxarray) :: ba
            type(amrex_distromap) :: new_dm
            type(amrex_mfiter) :: mfi
            type(part), dimension(:), pointer :: p
            integer(I8) :: np_
            ! Loop over levels
            do lvl=0,this%amr%clvl()
               ! Count number of boxes at that level and prepare costs array
               ba=this%amr%get_boxarray(lvl); nboxes=int(ba%nboxes())
               allocate(costs(nboxes));costs=0.0_WP
               ! Count particles per box via MFIter
               call this%mfiter_build(lvl,mfi)
               do while (mfi%next())
                  call this%get_particles(lvl,mfi,p,np_)
                  costs(mfi%grid_index()+1)=real(np_,WP)
               end do
               call this%mfiter_destroy(mfi)
               ! Global sum
               call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
               ! Build knapsack DM from costs
               call amrdm_make_knapsack(new_dm,costs,nboxes)
               ! Apply to particle container
               call this%set_particle_dm(lvl,new_dm)
               ! Clean up
               call amrex_distromap_destroy(new_dm)
               deallocate(costs)
            end do
            ! Redistribute particles onto new DM
            call this%redistribute()
         end block rebalance_particles
      end if

      ! Recompute VF from particle positions
      call this%update_VF()

   end subroutine post_regrid

   !> Tag cells for refinement where particle VF exceeds VF_tag
   subroutine tagging(this,lvl,time,tags)
      use amrex_amr_module, only: amrex_tagboxarray,amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use amrgrid_class, only: SETtag
      implicit none
      class(amrlpt), intent(inout) :: this
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

   !> Estimate per-box cost for load balancing based on particle count
   !> This assumes that particles don't change level much...
   subroutine get_cost(this,lvl,nboxes,costs,ba)
      use amrex_amr_module, only: amrex_boxarray,amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: old_bx,new_bx
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: n,m,ierr
      integer, dimension(3) :: ijk
      real(WP) :: dxi,dyi,dzi
      ! Inverse cell size at this level
      dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
      ! Zero costs, then accumulate particle counts
      costs=0.0_WP
      call this%mfiter_build(lvl,mfi)
      do while (mfi%next())
         old_bx=mfi%tilebox()
         call this%get_particles(lvl,mfi,p,np_)
         do m=1,np_
            if (p(m)%flag.eq.PART_IS_DEAD) cycle
            ! Compute cell index
            ijk(1)=floor((p(m)%pos(1)-this%amr%xlo)*dxi)
            ijk(2)=floor((p(m)%pos(2)-this%amr%ylo)*dyi)
            ijk(3)=floor((p(m)%pos(3)-this%amr%zlo)*dzi)
            ! Find which proposed new box owns this particle
            do n=1,nboxes
               new_bx=ba%get_box(n-1)
               if (.not.old_bx%intersects(new_bx)) cycle
               if (new_bx%contains(ijk)) then
                  costs(n)=costs(n)+1.0_WP; exit
               end if
            end do
         end do
      end do
      call this%mfiter_destroy(mfi)
      ! Global sum
      call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine get_cost

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize amrlpt solver
   subroutine initialize(this,amr,name)
      use amrdata_class, only: amrex_bc_foextrap
      implicit none
      class(amrlpt), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), optional :: name
      ! Assign solver name
      if (present(name)) this%name=trim(adjustl(name))
      ! Point to amr
      this%amr=>amr
      ! Create particle container
      call amrlpt_new_pc(this%pc,this%amr%amrcore)
      ! Initialize VF field
      call this%VF%initialize(amr=amr,name='VF',ncomp=1,ng=this%nover); call this%VF%register()
      if (.not.this%amr%xper) then; this%VF%lo_bc(1,1)=amrex_bc_foextrap; this%VF%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.this%amr%yper) then; this%VF%lo_bc(2,1)=amrex_bc_foextrap; this%VF%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.this%amr%zper) then; this%VF%lo_bc(3,1)=amrex_bc_foextrap; this%VF%hi_bc(3,1)=amrex_bc_foextrap; end if
      ! Initialize source terms
      call this%src%initialize(amr=amr,name='src',ncomp=3,ng=this%nover); call this%src%register()
      if (.not.this%amr%xper) then; this%src%lo_bc(1,:)=amrex_bc_foextrap; this%src%hi_bc(1,:)=amrex_bc_foextrap; end if
      if (.not.this%amr%yper) then; this%src%lo_bc(2,:)=amrex_bc_foextrap; this%src%hi_bc(2,:)=amrex_bc_foextrap; end if
      if (.not.this%amr%zper) then; this%src%lo_bc(3,:)=amrex_bc_foextrap; this%src%hi_bc(3,:)=amrex_bc_foextrap; end if
      ! Register post-regrid callback
      select type (this)
       type is (amrlpt)
         call this%amr%add_postregrid(amrlpt_postregrid,c_loc(this))
         call this%amr%add_tagging   (amrlpt_tagging,   c_loc(this))
         call this%amr%set_get_cost  (amrlpt_get_cost,  c_loc(this))
      end select
      ! Initialize per-phase timers (unsynchronized for imbalance detection)
      this%tmr_coll =timer(comm=this%amr%comm,name='coll', sync=.false.)
      this%tmr_coll_=timer(comm=this%amr%comm,name='coll_',sync=.false.)
      this%tmr_fill =timer(comm=this%amr%comm,name='fill', sync=.false.)
      this%tmr_nbl  =timer(comm=this%amr%comm,name='nbl',  sync=.false.)
      this%tmr_step =timer(comm=this%amr%comm,name='step', sync=.false.)
      this%tmr_step_=timer(comm=this%amr%comm,name='step_',sync=.false.)
      this%tmr_vf   =timer(comm=this%amr%comm,name='VF',   sync=.false.)
      this%tmr_src  =timer(comm=this%amr%comm,name='src',  sync=.false.)
      ! Print out info
      call this%print()
   end subroutine initialize

   !> Finalize: destroy particle container and release grid pointer
   subroutine finalize(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call this%VF%finalize()
      call this%src%finalize()
      nullify(this%inject)
      nullify(this%user_lpt_tagging)
      call amrlpt_delete_pc(this%pc)
      this%pc=c_null_ptr
      nullify(this%amr)
      call this%tmr_coll%finalize()
      call this%tmr_coll_%finalize()
      call this%tmr_fill%finalize()
      call this%tmr_nbl%finalize()
      call this%tmr_step%finalize()
      call this%tmr_step_%finalize()
      call this%tmr_vf%finalize()
      call this%tmr_src%finalize()
   end subroutine finalize

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Soft-sphere collision model: computes Acol, Tcol on every valid particle with flag PART_COLLIDES set
   subroutine collide(this,dt,Gib,Gibcomp)
      use amrex_amr_module, only: amrex_mfiter,amrex_multifab,amrex_multifab_build,amrex_multifab_destroy,amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use amrdata_class, only: amrdata
      use mathtools, only: Pi,cross_product
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in), optional :: Gib
      integer, intent(in), optional :: Gibcomp
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pG
      integer(c_int32_t), dimension(:), pointer :: nbr_off,nbr_lst
      integer(I8) :: np_v,np_all
      integer :: lvl,no,i1,k
      integer(c_int32_t) :: j_c
      real(WP) :: d1,m1,d2,m2
      real(WP), dimension(3) :: r1,v1,w1,r2,v2,w2
      real(WP) :: k_coeff,eta_coeff,k_coeff_w,eta_coeff_w
      real(WP) :: dx,dy,dz
      logical :: hit,dual_dm
      integer :: gc
      type(amrex_multifab) :: tmpGib

      ! Start collision timer
      call this%tmr_coll%start()

      ! Precompute spring/damping coefficients
      k_coeff=(Pi**2+log(this%e_n)**2)/this%tau_col**2
      eta_coeff=-2.0_WP*log(this%e_n)/this%tau_col
      k_coeff_w=(Pi**2+log(this%e_w)**2)/this%tau_col**2
      eta_coeff_w=-2.0_WP*log(this%e_w)/this%tau_col

      ! Resolve Gib component
      gc=1; if (present(Gibcomp)) gc=Gibcomp

      ! Reset collision counter
      this%ncol=0

      ! Ghost fill and neighbor list radius: dmax+r_influ_max=1.2*dmax
      call this%tmr_fill%start()
      call this%fill_ghosts_radius(1.2_WP*this%dmax)
      call this%tmr_fill%stop()
      call this%tmr_nbl%start()
      call this%build_neighbor_list(rcrit=1.2_WP*this%dmax)
      call this%tmr_nbl%stop()

      ! Loop over levels
      call this%tmr_coll_%start()
      do lvl=0,this%amr%clvl()

         ! Get mesh size
         dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)

         ! Check if we need dual-grid bridge for Gib
         if (present(Gib)) then
            dual_dm=(.not.(Gib%mf(lvl)%dm.eq.this%get_particle_dm(lvl)))
            if (dual_dm) then
               call amrex_multifab_build(mf=tmpGib,ba=this%amr%ba(lvl),dm=this%get_particle_dm(lvl),nc=Gib%mf(lvl)%ncomp(),ng=Gib%mf(lvl)%nghost(),nodal=Gib%nodal)
               call tmpGib%parallel_copy(Gib%mf(lvl),1,1,Gib%mf(lvl)%ncomp(),Gib%mf(lvl)%nghost(),Gib%mf(lvl)%nghost(),this%amr%geom(lvl))
            end if
         end if

         ! Loop over tiles
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())

            ! Get combined valid+ghost particle array and neighbor list
            call this%get_all_particles(lvl,mfi,p,np_v,np_all)
            call this%get_neighbor_list(lvl,mfi,nbr_off,nbr_lst)

            ! Gib array pointer for IB collision
            if (present(Gib)) then
               if (dual_dm) then
                  pG=>tmpGib%dataptr(mfi)
               else
                  pG=>Gib%mf(lvl)%dataptr(mfi)
               end if
            end if

            ! Loop over valid particles
            do i1=1,np_v

               ! Only particles with PART_COLLIDES
               if (IAND(p(i1)%flag,PART_COLLIDES).eq.0) cycle

               ! Zero collision acceleration and torque
               p(i1)%Acol=0.0_WP; p(i1)%Tcol=0.0_WP

               ! Set first particle properties
               r1=p(i1)%pos; v1=p(i1)%vel; w1=p(i1)%angVel
               d1=p(i1)%d;   m1=this%rho*Pi/6.0_WP*d1**3

               ! Wall and IB collisions: r2=contact point, d2=0, v2=w2=0, and m2=huge so that m_eff=m1
               v2=0.0_WP; w2=0.0_WP; d2=0.0_WP; m2=huge(m1)
               if (this%lo_bc(1).eq.AMRLPT_WALL) then; r2=[this%amr%xlo,r1(2),r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%hi_bc(1).eq.AMRLPT_WALL) then; r2=[this%amr%xhi,r1(2),r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%lo_bc(2).eq.AMRLPT_WALL) then; r2=[r1(1),this%amr%ylo,r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%hi_bc(2).eq.AMRLPT_WALL) then; r2=[r1(1),this%amr%yhi,r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%lo_bc(3).eq.AMRLPT_WALL) then; r2=[r1(1),r1(2),this%amr%zlo]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%hi_bc(3).eq.AMRLPT_WALL) then; r2=[r1(1),r1(2),this%amr%zhi]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (present(Gib)) then
                  ib_col: block
                     real(WP) :: d_ib,buf
                     real(WP), dimension(3) :: pos_p,pos_m,n12_out
                     ! Signed distance and outward IB normal nabla G/|nabla G| (from IB into fluid)
                     d_ib=this%interp(lvl,r1,pG,gc)
                     pos_p=[r1(1)+0.5_WP*dx,r1(2),r1(3)]; pos_m=[r1(1)-0.5_WP*dx,r1(2),r1(3)]; n12_out(1)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dx
                     pos_p=[r1(1),r1(2)+0.5_WP*dy,r1(3)]; pos_m=[r1(1),r1(2)-0.5_WP*dy,r1(3)]; n12_out(2)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dy
                     pos_p=[r1(1),r1(2),r1(3)+0.5_WP*dz]; pos_m=[r1(1),r1(2),r1(3)-0.5_WP*dz]; n12_out(3)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dz
                     buf=norm2(n12_out)+epsilon(1.0_WP); n12_out=n12_out/buf; r2=r1-d_ib*n12_out
                     call col_force(k_coeff_w,eta_coeff_w)
                  end block ib_col
               end if

               ! Particle-particle collisions
               if (associated(nbr_off)) then
                  do k=nbr_off(i1)+1,nbr_off(i1+1)
                     ! Get neighbor index (note the 1-shift)
                     j_c=nbr_lst(k)+1
                     ! Skip self-collision
                     if (j_c.eq.i1) cycle
                     ! Skip if neighbor doesn't collide
                     if (IAND(p(j_c)%flag,PART_COLLIDES).eq.0) cycle
                     ! Set neighbor properties
                     r2=p(j_c)%pos; v2=p(j_c)%vel; w2=p(j_c)%angVel
                     d2=p(j_c)%d;   m2=this%rho*Pi/6.0_WP*d2**3
                     ! Compute collision force
                     call col_force(k_coeff,eta_coeff,hit)
                     if (hit) this%ncol=this%ncol+1
                  end do
               end if

               ! Suppress forces in collapsed directions
               if (this%amr%nx.eq.1) then; p(i1)%Acol(1)=0.0_WP; p(i1)%Tcol(2)=0.0_WP; p(i1)%Tcol(3)=0.0_WP; end if
               if (this%amr%ny.eq.1) then; p(i1)%Acol(2)=0.0_WP; p(i1)%Tcol(1)=0.0_WP; p(i1)%Tcol(3)=0.0_WP; end if
               if (this%amr%nz.eq.1) then; p(i1)%Acol(3)=0.0_WP; p(i1)%Tcol(1)=0.0_WP; p(i1)%Tcol(2)=0.0_WP; end if

            end do

         end do
         call this%mfiter_destroy(mfi)

         ! Clean up temporary Gib multifab if used
         if (dual_dm.and.present(Gib)) call amrex_multifab_destroy(tmpGib)

      end do
      call this%tmr_coll_%stop()

      ! Clear ghosts
      call this%clear_ghosts()

      ! Sum collision count (each pair counted once per partner => divide by 2)
      reduce_collision_count: block
         use mpi_f08, only: MPI_IN_PLACE,MPI_INTEGER,MPI_SUM
         integer :: ierr
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%ncol,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         this%ncol=this%ncol/2
      end block reduce_collision_count

      ! Stop collision timer
      call this%tmr_coll%stop()

   contains

      !> Soft-sphere collision force between (r1,v1,w1,d1,m1) and virtual partner (r2,v2,w2,d2,m2).
      !> n12=(r2-r1)/|r2-r1| computed internally => points from particle toward partner.
      !> Wall/IB: set r2=contact point, d2=0, v2=w2=0, m2=huge => m_eff->m1, d_eff=0.5*d1.
      subroutine col_force(kk,ee,collided)
         real(WP), intent(in)  :: kk,ee
         logical,  intent(out), optional :: collided
         real(WP) :: d12,d_eff,rnv,r_influ,delta_n,rtv,m_eff
         real(WP), dimension(3) :: n12,v12,t12,f_n,f_t
         real(WP), parameter :: aclipnorm=1.0e-6_WP,acliptan=1.0e-9_WP,rcliptan=0.05_WP
         ! No collision yet
         if (present(collided)) collided=.false.
         ! Get distance
         d12=norm2(r2-r1)
         ! Skip if particles are too close - likely self-collision
         if (d12.lt.10.0_WP*epsilon(d12)) return
         ! Get normal
         n12=(r2-r1)/d12
         ! Get effective diameter, relative velocity, and relative normal velocity
         d_eff=0.5_WP*(d1+d2); v12=v1-v2; rnv=dot_product(v12,n12)
         ! Get influence radius
         r_influ=min(abs(rnv)*dt,0.2_WP*d_eff)
         ! Get overlap
         delta_n=min(d_eff+r_influ-d12,this%clip_col*d_eff)
         ! Done if no overlap
         if (delta_n.le.0.0_WP) return
         ! Collision detected
         if (present(collided)) collided=.true.
         ! Get effective mass
         m_eff=m1*m2/(m1+m2)
         ! Get tangential velocity
         t12=v12-rnv*n12+cross_product(0.5_WP*(d1*w1+d2*w2),n12)
         ! Get normal force
         f_n=(-m_eff*kk*delta_n-m_eff*ee*rnv)*n12
         ! Get tangential force
         f_t=0.0_WP
         if (this%mu_f.gt.0.0_WP) then
            rtv=sqrt(sum(t12*t12))
            if (rnv*dt/d_eff.gt.aclipnorm) then; if (   rtv/rnv  .lt.rcliptan) rtv=0.0_WP
            else;                                if (dt*rtv/d_eff.lt.acliptan) rtv=0.0_WP; end if
            if (rtv.gt.0.0_WP) f_t=-this%mu_f*norm2(f_n)*t12/rtv
         end if
         ! Increment accelerations on p(i1)
         p(i1)%Acol=p(i1)%Acol+(f_n+f_t)/m1
         p(i1)%Tcol=p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t/m1)
      end subroutine col_force

   end subroutine collide

   !> Advance all particles by dt using Euler-midpoint substepping
   !> U,V,W: velocity amrdata (staggered MAC or collocated, with ghosts filled)
   !> rho,visc: cell-centered amrdata (with ghosts filled)
   !> Ucomp/Vcomp/Wcomp: component to use from each velocity field (optional, default 1)
   !> rhocomp/visccomp: component to use from rho/visc fields (optional, default 1)
   !> cst_rho/cst_visc: constant scalar alternatives to rho/visc amrdata
   subroutine advance(this,dt,U,Ucomp,V,Vcomp,W,Wcomp,rho,rhocomp,visc,visccomp,cst_rho,cst_visc)
      use amrex_amr_module, only: amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use amrex_interface, only: amrmfab_parallel_add
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in) :: U,V,W
      type(amrdata), intent(in), optional :: rho,visc
      real(WP),      intent(in), optional :: cst_rho,cst_visc
      integer,       intent(in), optional :: Ucomp,Vcomp,Wcomp,rhocomp,visccomp
      type(amrex_distromap), dimension(:), allocatable :: pdms
      logical :: any_dual_dm
      type(amrdata) :: U_pdm,V_pdm,W_pdm,rho_pdm,visc_pdm,VF_pdm,src_pdm
      integer :: lvl

      ! Build per-level particle DM array and check for dual-grid
      allocate(pdms(0:this%amr%clvl()))
      any_dual_dm=.false.
      do lvl=0,this%amr%clvl()
         pdms(lvl)=this%get_particle_dm(lvl)
         if (.not.(this%amr%dm(lvl).eq.pdms(lvl))) any_dual_dm=.true.
      end do

      ! Clone fields onto particle DMs if needed
      if (any_dual_dm) then
         call U%clone(U_pdm,dm=pdms)
         call V%clone(V_pdm,dm=pdms)
         call W%clone(W_pdm,dm=pdms)
         if (present(rho))  call rho%clone(rho_pdm,dm=pdms)
         if (present(visc)) call visc%clone(visc_pdm,dm=pdms)
         call this%VF%clone(VF_pdm,dm=pdms)
         call this%src%clone(src_pdm,dm=pdms)
      end if

      ! Reset sources
      call this%src%setval(0.0_WP)
      if (any_dual_dm) call src_pdm%setval(0.0_WP)

      ! Perform one step
      if (any_dual_dm) then
         if (present(rho).and.present(visc)) then
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,visc=visc_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp)
         else if (present(rho)) then
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp)
         else if (present(visc)) then
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,visc=visc_pdm,cst_rho=cst_rho,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,visccomp=visccomp)
         else
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp)
         end if
      else
         call this%step(dt=dt,U=U,V=V,W=W,VF=this%VF,src=this%src,rho=rho,visc=visc,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp)
      end if

      ! Writeback src_pdm to src and clean up clones
      if (any_dual_dm) then
         do lvl=0,this%amr%clvl()
            call amrmfab_parallel_add(this%src%mf(lvl),src_pdm%mf(lvl),1,1,this%src%mf(lvl)%ncomp(),0,0,this%amr%geom(lvl))
         end do
         call U_pdm%finalize()
         call V_pdm%finalize()
         call W_pdm%finalize()
         if (present(rho))  call rho_pdm%finalize()
         if (present(visc)) call visc_pdm%finalize()
         call VF_pdm%finalize()
         call src_pdm%finalize()
      end if
      deallocate(pdms)

      ! Recompute particle volume fraction
      call this%update_VF()

      ! Process accumulated sources
      call this%tmr_src%start()
      process_sources: block
         ! Process deposit
         call this%process_deposit(this%src)
         ! Fill ghost cells
         call this%src%fill(time=0.0_WP)
         ! Filter
         call this%filter(this%src)
      end block process_sources
      call this%tmr_src%stop()

      ! Log particle advance
      log_particle_advance: block
         use iso_fortran_env, only: output_unit
         use string,          only: str_long
         use param,           only: verbose
         use messager,        only: log
         character(len=str_long) :: message
         if (verbose.gt.0) then
            write(message,'(" [",a,"] LPT advance | Np=",I0," | dt=",ES9.3)') trim(this%name),this%np,dt
            if (verbose.gt.1.and.this%amr%amRoot) write(output_unit,'(a)') trim(message)
            call log(message)
         end if
      end block log_particle_advance
      
   end subroutine advance

   !> Advance particles to a target time using CFL-limited sub-steps:
   !> Sub-step size is adapted from this%dt (persistent) using this%cflmax
   !> On exit, this%t=t_target and this%dt holds the last adapted sub-step
   subroutine advance_to(this,time,do_collide,U,V,W,rho,visc,cst_rho,cst_visc,Ucomp,Vcomp,Wcomp,rhocomp,visccomp,Gib,Gibcomp)
      use amrex_amr_module, only: amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use amrex_interface, only: amrmfab_parallel_add
      use messager, only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: time
      logical,  intent(in) :: do_collide
      type(amrdata), intent(in) :: U,V,W
      type(amrdata), intent(in), optional :: rho,visc
      real(WP),      intent(in), optional :: cst_rho,cst_visc
      integer,       intent(in), optional :: Ucomp,Vcomp,Wcomp,rhocomp,visccomp
      type(amrdata), intent(in), optional :: Gib
      integer,       intent(in), optional :: Gibcomp
      real(WP) :: dt_done,mydt,cfl
      integer  :: n_sub,lvl
      type(amrex_distromap), dimension(:), allocatable :: pdms
      logical :: any_dual_dm
      type(amrdata) :: U_pdm,V_pdm,W_pdm,rho_pdm,visc_pdm,Gib_pdm,VF_pdm,src_pdm

      ! Validate time target
      if (time.le.this%t) return

      ! CFL-adapt the persistent sub-step size
      call this%get_cfl(this%dt,cflc=cfl,cfl=cfl)
      if (cfl.gt.0.0_WP) this%dt=min(this%dt*this%cflmax/cfl,this%dtmax)

      ! Build per-level particle DM array and check for dual-grid
      allocate(pdms(0:this%amr%clvl()))
      any_dual_dm=.false.
      do lvl=0,this%amr%clvl()
         pdms(lvl)=this%get_particle_dm(lvl)
         if (.not.(this%amr%dm(lvl).eq.pdms(lvl))) any_dual_dm=.true.
      end do

      ! Clone fields onto particle DMs if needed (ONE TIME before sub-step loop)
      if (any_dual_dm) then
         call U%clone(U_pdm,dm=pdms)
         call V%clone(V_pdm,dm=pdms)
         call W%clone(W_pdm,dm=pdms)
         if (present(rho))  call rho%clone(rho_pdm,dm=pdms)
         if (present(visc)) call visc%clone(visc_pdm,dm=pdms)
         if (present(Gib))  call Gib%clone(Gib_pdm,dm=pdms)
         call this%VF%clone(VF_pdm,dm=pdms)
         call this%src%clone(src_pdm,dm=pdms)
      end if

      ! Reset sources
      call this%src%setval(0.0_WP)
      if (any_dual_dm) call src_pdm%setval(0.0_WP)

      ! Sub-step loop: inject, collide, step, accumulate src
      dt_done=0.0_WP; n_sub=0
      do while (dt_done.lt.(time-this%t)-epsilon(1.0_WP))
         ! Select mydt
         mydt=min(this%dt,time-this%t-dt_done)
         ! Inject if needed
         if (associated(this%inject)) call this%inject(dt=mydt)
         ! Collide if needed
         if (do_collide) then
            if (any_dual_dm.and.present(Gib)) then
               call this%collide(dt=mydt,Gib=Gib_pdm,Gibcomp=Gibcomp)
            else
               call this%collide(dt=mydt,Gib=Gib,Gibcomp=Gibcomp)
            end if
         end if
         ! Step particles
         if (any_dual_dm) then
            if (present(rho).and.present(visc)) then
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,visc=visc_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp)
            else if (present(rho)) then
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp)
            else if (present(visc)) then
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,visc=visc_pdm,cst_rho=cst_rho,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,visccomp=visccomp)
            else
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp)
            end if
         else
            call this%step(dt=mydt,U=U,V=V,W=W,VF=this%VF,src=this%src,rho=rho,visc=visc,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp)
         end if
         ! Increment time and sub-step counter
         dt_done=dt_done+mydt; n_sub=n_sub+1
      end do

      ! Writeback src_pdm to src and clean up clones
      if (any_dual_dm) then
         do lvl=0,this%amr%clvl()
            call amrmfab_parallel_add(this%src%mf(lvl),src_pdm%mf(lvl),1,1,this%src%mf(lvl)%ncomp(),0,0,this%amr%geom(lvl))
         end do
         call U_pdm%finalize()
         call V_pdm%finalize()
         call W_pdm%finalize()
         if (present(rho))  call rho_pdm%finalize()
         if (present(visc)) call visc_pdm%finalize()
         if (present(Gib))  call Gib_pdm%finalize()
         call VF_pdm%finalize()
         call src_pdm%finalize()
      end if
      deallocate(pdms)

      ! Recompute particle volume fraction
      call this%update_VF()

      ! Process accumulated sources
      call this%tmr_src%start()
      process_sources: block
         ! Process deposit
         call this%process_deposit(this%src)
         ! Fill ghost cells
         call this%src%fill(time=0.0_WP)
         ! Filter
         call this%filter(this%src)
      end block process_sources
      call this%tmr_src%stop()

      ! Advance particle time
      this%t=time

      ! Log particle advance
      log_particle_advance: block
         use iso_fortran_env, only: output_unit
         use string,          only: str_long
         use param,           only: verbose
         use messager,        only: log
         character(len=str_long) :: message
         if (verbose.gt.0) then
            write(message,'(" [",a,"] ",i0," LPT sub-steps | Np=",I0," | dt=",ES9.3," | CFL=",F6.4)') trim(this%name),n_sub,this%np,this%dt,cfl
            if (verbose.gt.1.and.this%amr%amRoot) write(output_unit,'(a)') trim(message)
            call log(message)
         end if
      end block log_particle_advance

   end subroutine advance_to

   !> Step all particles on all AMR levels by dt, deposit momentum source
   !> All amrdata arguments must be on the particle DM (caller ensures via clone)
   subroutine step(this,dt,U,V,W,VF,src,rho,visc,cst_rho,cst_visc,Ucomp,Vcomp,Wcomp,rhocomp,visccomp)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mathtools, only: Pi
      use messager,  only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in) :: U,V,W,VF
      type(amrdata), intent(inout) :: src
      type(amrdata), intent(in), optional :: rho,visc
      real(WP),      intent(in), optional :: cst_rho,cst_visc
      integer,       intent(in), optional :: Ucomp,Vcomp,Wcomp,rhocomp,visccomp

      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRho,pVisc,pVolFrac,pSrc
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: lvl,i,uc,vc,wc,rhoc,viscc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,mydt,dt_done,Ip
      real(WP), dimension(3) :: acc,dmom
      type(part) :: myp,pold
      logical :: is_stag

      ! Start step timer
      call this%tmr_step%start()

      ! Resolve optional component indices
      uc=1; if (present(Ucomp)) uc=Ucomp
      vc=1; if (present(Vcomp)) vc=Vcomp
      wc=1; if (present(Wcomp)) wc=Wcomp
      rhoc =1; if (present(rhocomp))  rhoc =rhocomp
      viscc=1; if (present(visccomp)) viscc=visccomp

      ! Validate: each of rho and visc must be provided in exactly one form
      if (.not.present(rho) .and..not.present(cst_rho))  call die('[amrlpt step] rho or cst_rho required')
      if (.not.present(visc).and..not.present(cst_visc)) call die('[amrlpt step] visc or cst_visc required')

      ! Check velocity nodal locations
      check_velocity: block
         logical, dimension(3) :: nU,nV,nW
         nU=U%nodal; nV=V%nodal; nW=W%nodal
         if (all(nU.eqv.[.true.,.false.,.false.]).and. &
         &   all(nV.eqv.[.false.,.true.,.false.]).and. &
         &   all(nW.eqv.[.false.,.false.,.true.])) then
            is_stag=.true.
         else if (.not.any(nU).and..not.any(nV).and..not.any(nW)) then
            is_stag=.false.
         else
            call die('[amrlpt step] U/V/W must be staggered (face-centered) or collocated (cell-centered)')
         end if
      end block check_velocity

      ! Loop over all AMR levels
      call this%tmr_step_%start()
      do lvl=0,this%amr%clvl()

         ! Get mesh size
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz

         ! MFIter over the level
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())

            ! Get pointers to data
            pU      =>U%mf(lvl)%dataptr(mfi)
            pV      =>V%mf(lvl)%dataptr(mfi)
            pW      =>W%mf(lvl)%dataptr(mfi)
            if (present(rho))  pRho =>rho%mf(lvl)%dataptr(mfi)
            if (present(visc)) pVisc=>visc%mf(lvl)%dataptr(mfi)
            pVolFrac=>VF%mf(lvl)%dataptr(mfi)
            pSrc    =>src%mf(lvl)%dataptr(mfi)

            ! Get particles on this tile
            call this%get_particles(lvl=lvl,mfi=mfi,p=p,np=np_)

            ! Loop over local particles
            do i=1,np_
               ! Skip particles that are not moving or exchanging
               if (IAND(p(i)%flag,PART_MOVES+PART_EXCHANGES).eq.0) cycle
               ! Create copy of particle
               myp=p(i)
               ! Time-integrate until dt_done=dt
               dt_done=0.0_WP
               do while (dt_done.lt.dt)
                  mydt=min(myp%dt,dt-dt_done)
                  if (mydt.le.0.0_WP) mydt=dt-dt_done
                  ! Remember the particle
                  pold=myp
                  ! Precompute moment of inertia for a sphere
                  Ip=0.1_WP*myp%d**2
                  ! Advance with Euler prediction
                  acc=get_rhs()
                  if (IAND(myp%flag,PART_MOVES).ne.0) then
                     myp%pos=pold%pos+0.5_WP*mydt*myp%vel
                     myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity+myp%Acol)
                     myp%angVel=pold%angVel+0.5_WP*mydt*myp%Tcol/Ip
                  end if
                  acc=get_rhs()
                  if (IAND(myp%flag,PART_MOVES).ne.0) then
                     myp%pos=pold%pos+mydt*myp%vel
                     myp%vel=pold%vel+mydt*(acc+this%gravity+myp%Acol)
                     myp%angVel=pold%angVel+mydt*myp%Tcol/Ip
                  end if

                  ! DEBUG: catch bad deposit
                  if (floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP).lt.lbound(pSrc,1) .or. &
                  &   floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP)+1.gt.ubound(pSrc,1) .or. &
                  &   myp%pos(1).ne.myp%pos(1)) then
                     write(*,'(A,I0,A,I0)') 'lvl=',lvl,' i=',i
                     write(*,'(A,3ES18.10)') 'pold%pos=',pold%pos
                     write(*,'(A,3ES18.10)') ' myp%pos=',myp%pos
                     write(*,'(A,3ES18.10)') ' myp%vel=',myp%vel
                     write(*,'(A,3ES18.10)') '     acc=',acc
                     write(*,'(A,3ES18.10)') '    Acol=',myp%Acol
                     write(*,'(A,ES18.10)')  '      dt=',mydt
                     write(*,'(A,I0,A,I0)')  'pSrc x bounds: ',lbound(pSrc,1),' to ',ubound(pSrc,1)
                     error stop 'BAD DEPOSIT'
                  end if

                  ! Transfer back to the mesh
                  if (IAND(myp%flag,PART_EXCHANGES).ne.0) then
                     dmom=mydt*acc*this%rho*Pi/6.0_WP*myp%d**3
                     call deposit(val=-dmom)
                  end if
                  ! Increment
                  dt_done=dt_done+mydt
               end do
               ! Track escape from non-periodic boundaries
               if ((.not.this%amr%xper.and.(myp%pos(1).lt.this%amr%xlo.or.myp%pos(1).gt.this%amr%xhi)).or. &
               &   (.not.this%amr%yper.and.(myp%pos(2).lt.this%amr%ylo.or.myp%pos(2).gt.this%amr%yhi)).or. &
               &   (.not.this%amr%zper.and.(myp%pos(3).lt.this%amr%zlo.or.myp%pos(3).gt.this%amr%zhi))) then
                  this%np_out_loc=this%np_out_loc+1
                  this%Vp_out_loc=this%Vp_out_loc+Pi/6.0_WP*myp%d**3
               end if
               ! Write back
               p(i)=myp
            end do
         end do
         call this%mfiter_destroy(mfi)

      end do
      call this%tmr_step_%stop()

      ! Redistribute particles
      call this%redistribute()

      ! Stop step timer
      call this%tmr_step%stop()

   contains

      !> Calculate rhs of particle equations of motion
      function get_rhs() result(acc)
         real(WP), dimension(3) :: acc,fvel
         real(WP) :: pVF,fVF,frho,fvisc,Re,tau,corr,b1,b2
         real(WP) :: wx,wy,wz,fx,fy,fz
         integer  :: ic,jc,kc,ix,iy,iz
         ! Cell-centered stencil indices and weights
         ic=floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(myp%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
         jc=floor((myp%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(myp%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
         kc=floor((myp%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(myp%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
         ! Velocity interpolation
         if (is_stag) then
            ! Face-centered indices and weights
            ix=floor((myp%pos(1)-this%amr%xlo)*dxi); fx=(myp%pos(1)-this%amr%xlo)*dxi-real(ix,WP)
            iy=floor((myp%pos(2)-this%amr%ylo)*dyi); fy=(myp%pos(2)-this%amr%ylo)*dyi-real(iy,WP)
            iz=floor((myp%pos(3)-this%amr%zlo)*dzi); fz=(myp%pos(3)-this%amr%zlo)*dzi-real(iz,WP)
            fvel(1)=(1.0_WP-fx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix,jc,kc,uc)+fx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix+1,jc,kc,uc)+(1.0_WP-fx)*wy*(1.0_WP-wz)*pU(ix,jc+1,kc,uc)+fx*wy*(1.0_WP-wz)*pU(ix+1,jc+1,kc,uc)+(1.0_WP-fx)*(1.0_WP-wy)*wz*pU(ix,jc,kc+1,uc)+fx*(1.0_WP-wy)*wz*pU(ix+1,jc,kc+1,uc)+(1.0_WP-fx)*wy*wz*pU(ix,jc+1,kc+1,uc)+fx*wy*wz*pU(ix+1,jc+1,kc+1,uc)
            fvel(2)=(1.0_WP-wx)*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic,iy,kc,vc)+wx*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic+1,iy,kc,vc)+(1.0_WP-wx)*fy*(1.0_WP-wz)*pV(ic,iy+1,kc,vc)+wx*fy*(1.0_WP-wz)*pV(ic+1,iy+1,kc,vc)+(1.0_WP-wx)*(1.0_WP-fy)*wz*pV(ic,iy,kc+1,vc)+wx*(1.0_WP-fy)*wz*pV(ic+1,iy,kc+1,vc)+(1.0_WP-wx)*fy*wz*pV(ic,iy+1,kc+1,vc)+wx*fy*wz*pV(ic+1,iy+1,kc+1,vc)
            fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic,jc,iz,wc)+wx*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic+1,jc,iz,wc)+(1.0_WP-wx)*wy*(1.0_WP-fz)*pW(ic,jc+1,iz,wc)+wx*wy*(1.0_WP-fz)*pW(ic+1,jc+1,iz,wc)+(1.0_WP-wx)*(1.0_WP-wy)*fz*pW(ic,jc,iz+1,wc)+wx*(1.0_WP-wy)*fz*pW(ic+1,jc,iz+1,wc)+(1.0_WP-wx)*wy*fz*pW(ic,jc+1,iz+1,wc)+wx*wy*fz*pW(ic+1,jc+1,iz+1,wc)
         else
            fvel(1)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ic,jc,kc,uc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ic+1,jc,kc,uc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pU(ic,jc+1,kc,uc)+wx*wy*(1.0_WP-wz)*pU(ic+1,jc+1,kc,uc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pU(ic,jc,kc+1,uc)+wx*(1.0_WP-wy)*wz*pU(ic+1,jc,kc+1,uc)+(1.0_WP-wx)*wy*wz*pU(ic,jc+1,kc+1,uc)+wx*wy*wz*pU(ic+1,jc+1,kc+1,uc)
            fvel(2)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pV(ic,jc,kc,vc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pV(ic+1,jc,kc,vc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pV(ic,jc+1,kc,vc)+wx*wy*(1.0_WP-wz)*pV(ic+1,jc+1,kc,vc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pV(ic,jc,kc+1,vc)+wx*(1.0_WP-wy)*wz*pV(ic+1,jc,kc+1,vc)+(1.0_WP-wx)*wy*wz*pV(ic,jc+1,kc+1,vc)+wx*wy*wz*pV(ic+1,jc+1,kc+1,vc)
            fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pW(ic,jc,kc,wc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pW(ic+1,jc,kc,wc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pW(ic,jc+1,kc,wc)+wx*wy*(1.0_WP-wz)*pW(ic+1,jc+1,kc,wc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pW(ic,jc,kc+1,wc)+wx*(1.0_WP-wy)*wz*pW(ic+1,jc,kc+1,wc)+(1.0_WP-wx)*wy*wz*pW(ic,jc+1,kc+1,wc)+wx*wy*wz*pW(ic+1,jc+1,kc+1,wc)
         end if
         ! Density and viscosity
         if (present(cst_rho)) then
            frho=cst_rho
         else
            frho=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pRho(ic,jc,kc,rhoc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pRho(ic+1,jc,kc,rhoc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pRho(ic,jc+1,kc,rhoc)+wx*wy*(1.0_WP-wz)*pRho(ic+1,jc+1,kc,rhoc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pRho(ic,jc,kc+1,rhoc)+wx*(1.0_WP-wy)*wz*pRho(ic+1,jc,kc+1,rhoc)+(1.0_WP-wx)*wy*wz*pRho(ic,jc+1,kc+1,rhoc)+wx*wy*wz*pRho(ic+1,jc+1,kc+1,rhoc)
         end if
         if (present(cst_visc)) then
            fvisc=cst_visc
         else
            fvisc=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pVisc(ic,jc,kc,viscc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pVisc(ic+1,jc,kc,viscc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pVisc(ic,jc+1,kc,viscc)+wx*wy*(1.0_WP-wz)*pVisc(ic+1,jc+1,kc,viscc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pVisc(ic,jc,kc+1,viscc)+wx*(1.0_WP-wy)*wz*pVisc(ic+1,jc,kc+1,viscc)+(1.0_WP-wx)*wy*wz*pVisc(ic,jc+1,kc+1,viscc)+wx*wy*wz*pVisc(ic+1,jc+1,kc+1,viscc)
         end if
         fvisc=fvisc+epsilon(1.0_WP)
         ! Volume fraction
         pVF=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pVolFrac(ic,jc,kc,1)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pVolFrac(ic+1,jc,kc,1)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pVolFrac(ic,jc+1,kc,1)+wx*wy*(1.0_WP-wz)*pVolFrac(ic+1,jc+1,kc,1)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pVolFrac(ic,jc,kc+1,1)+wx*(1.0_WP-wy)*wz*pVolFrac(ic+1,jc,kc+1,1)+(1.0_WP-wx)*wy*wz*pVolFrac(ic,jc+1,kc+1,1)+wx*wy*wz*pVolFrac(ic+1,jc+1,kc+1,1)
         fVF=max(1.0_WP-pVF,epsilon(1.0_WP))
         ! Drag correction factor
         select case(trim(this%drag_model))
         case('None','none')
            corr=epsilon(1.0_WP)
         case('Stokes')
            corr=1.0_WP
         case('Schiller-Naumann','Schiller Naumann','SN')
            Re=frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            corr=1.0_WP+0.15_WP*Re**(0.687_WP)
         case('Tenneti') ! Tenneti and Subramaniam (2011)
            Re=fVF*frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            b1=5.81_WP*pVF/fVF**3+0.48_WP*pVF**(1.0_WP/3.0_WP)/fVF**4
            b2=pVF**3*Re*(0.95_WP+0.61_WP*pVF**3/fVF**2)
            corr=fVF*((1.0_WP+0.15_WP*Re**(0.687_WP))/fVF**3+b1+b2)
         case('Beetstra') ! Beetstra et al. (2007)
            Re=fVF*frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            b1=10.0_WP*pVF/fVF**2+fVF**2*(1.0_WP+1.5_WP*sqrt(pVF))
            b2=0.413_WP/24.0_WP*Re/fVF**2*(1.0_WP/fVF+3.0_WP*fVF*pVF+8.4_WP*Re**(-0.343_WP))/(1.0_WP+10.0_WP**(3.0_WP*pVF)*Re**(2.0_WP*fVF-2.5_WP))
            corr=b1+b2
         case default
            corr=1.0_WP
         end select
         ! Particle response time
         tau=this%rho*myp%d**2/(18.0_WP*fvisc*corr)
         ! Return acceleration and update particle timestep size
         acc=(fvel-myp%vel)/tau
         myp%dt=tau/real(this%nstep,WP)
      end function get_rhs

      !> Deposit a particle source term back to the mesh
      subroutine deposit(val)
         real(WP), dimension(3), intent(in) :: val
         real(WP) :: wx,wy,wz
         integer  :: ic,jc,kc
         ! Cell-centered stencil indices and weights
         ic=floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(myp%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
         jc=floor((myp%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(myp%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
         kc=floor((myp%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(myp%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
         ! Clamp stencil at wall boundaries
         if (this%lo_bc(1).eq.AMRLPT_WALL.and.ic  .lt.this%amr%geom(lvl)%domain%lo(1)) then; ic=this%amr%geom(lvl)%domain%lo(1)  ; wx=0.0_WP; end if
         if (this%hi_bc(1).eq.AMRLPT_WALL.and.ic+1.gt.this%amr%geom(lvl)%domain%hi(1)) then; ic=this%amr%geom(lvl)%domain%hi(1)-1; wx=1.0_WP; end if
         if (this%lo_bc(2).eq.AMRLPT_WALL.and.jc  .lt.this%amr%geom(lvl)%domain%lo(2)) then; jc=this%amr%geom(lvl)%domain%lo(2)  ; wy=0.0_WP; end if
         if (this%hi_bc(2).eq.AMRLPT_WALL.and.jc+1.gt.this%amr%geom(lvl)%domain%hi(2)) then; jc=this%amr%geom(lvl)%domain%hi(2)-1; wy=1.0_WP; end if
         if (this%lo_bc(3).eq.AMRLPT_WALL.and.kc  .lt.this%amr%geom(lvl)%domain%lo(3)) then; kc=this%amr%geom(lvl)%domain%lo(3)  ; wz=0.0_WP; end if
         if (this%hi_bc(3).eq.AMRLPT_WALL.and.kc+1.gt.this%amr%geom(lvl)%domain%hi(3)) then; kc=this%amr%geom(lvl)%domain%hi(3)-1; wz=1.0_WP; end if
         ! Add to source
         pSrc(ic:ic+1,jc:jc+1,kc:kc+1,1)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,1)+reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])*val(1)
         pSrc(ic:ic+1,jc:jc+1,kc:kc+1,2)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,2)+reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])*val(2)
         pSrc(ic:ic+1,jc:jc+1,kc:kc+1,3)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,3)+reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])*val(3)
      end subroutine deposit

   end subroutine step

   !> Update particle volume fraction field based on our current particles
   subroutine update_VF(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_multifab_build,amrex_multifab_destroy,amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer(I8) :: np_
      integer :: lvl,i,ii,jj,kk
      real(WP) :: dxi,dyi,dzi,Vp,wx,wy,wz
      type(amrex_multifab) :: tmpVF
      logical :: dual_dm

      ! Start VF timer
      call this%tmr_vf%start()

      ! Zero VF on all levels
      call this%VF%setval(0.0_WP)
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         ! Get inverse of mesh size
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Check if we need dual-grid parallel_copy bridge
         dual_dm=(.not.(this%VF%mf(lvl)%dm.eq.this%get_particle_dm(lvl)))
         if (dual_dm) then
            call amrex_multifab_build(mf=tmpVF,ba=this%amr%ba(lvl),dm=this%get_particle_dm(lvl),nc=this%VF%mf(lvl)%ncomp(),ng=this%VF%mf(lvl)%nghost(),nodal=this%VF%nodal)
            call tmpVF%setval(0.0_WP)
         end if
         ! Loop over tiles
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointer to data
            if (dual_dm) then
               pVF=>tmpVF%dataptr(mfi)
            else
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
            end if
            ! Loop over particles
            call this%get_particles(lvl=lvl,mfi=mfi,p=p,np=np_)
            do i=1,np_
               ! Skip particles that don't exchange with the fluid
               if (IAND(p(i)%flag,PART_EXCHANGES).eq.0) cycle
               ! Get particle volume
               Vp=Pi/6.0_WP*p(i)%d**3
               ! Get indices and weights
               ii=floor((p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ii,WP)
               jj=floor((p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jj,WP)
               kk=floor((p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kk,WP)
               ! Clamp stencil at wall boundaries
               if (this%lo_bc(1).eq.AMRLPT_WALL.and.ii  .lt.this%amr%geom(lvl)%domain%lo(1)) then; ii=this%amr%geom(lvl)%domain%lo(1)  ; wx=0.0_WP; end if
               if (this%hi_bc(1).eq.AMRLPT_WALL.and.ii+1.gt.this%amr%geom(lvl)%domain%hi(1)) then; ii=this%amr%geom(lvl)%domain%hi(1)-1; wx=1.0_WP; end if
               if (this%lo_bc(2).eq.AMRLPT_WALL.and.jj  .lt.this%amr%geom(lvl)%domain%lo(2)) then; jj=this%amr%geom(lvl)%domain%lo(2)  ; wy=0.0_WP; end if
               if (this%hi_bc(2).eq.AMRLPT_WALL.and.jj+1.gt.this%amr%geom(lvl)%domain%hi(2)) then; jj=this%amr%geom(lvl)%domain%hi(2)-1; wy=1.0_WP; end if
               if (this%lo_bc(3).eq.AMRLPT_WALL.and.kk  .lt.this%amr%geom(lvl)%domain%lo(3)) then; kk=this%amr%geom(lvl)%domain%lo(3)  ; wz=0.0_WP; end if
               if (this%hi_bc(3).eq.AMRLPT_WALL.and.kk+1.gt.this%amr%geom(lvl)%domain%hi(3)) then; kk=this%amr%geom(lvl)%domain%hi(3)-1; wz=1.0_WP; end if
               ! Add to VF
               pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)=pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)+Vp*reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])
            end do
         end do
         call amrex_mfiter_destroy(mfi)
         ! If dual DM, parallel_copy VF back to Euler DM
         if (dual_dm) then
            call this%VF%mf(lvl)%parallel_copy(tmpVF,1,1,this%VF%mf(lvl)%ncomp(),this%VF%mf(lvl)%nghost(),this%VF%mf(lvl)%nghost(),this%amr%geom(lvl))
            call amrex_multifab_destroy(tmpVF)
         end if
      end do
      ! Post-process deposit
      call this%process_deposit(this%VF)
      ! Fill ghost cells
      call this%VF%fill(time=0.0_WP)
      ! Filter
      call this%filter(this%VF)

      ! Stop VF timer
      call this%tmr_vf%stop()

   end subroutine update_VF

   !> CFL based on particle velocities relative to their local cell size
   subroutine get_cfl(this,dt,cflc,cfl)
      use amrex_amr_module, only: amrex_mfiter
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional    :: cfl
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: lvl,n,ierr
      ! Initialize
      this%CFLp_x=0.0_WP; this%CFLp_y=0.0_WP; this%CFLp_z=0.0_WP; this%CFL_col=0.0_WP
      ! Loop over levels and tiles
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Loop over particles
            call this%get_particles(lvl,mfi,p,np_)
            do n=1,np_
               ! Skip dead particles
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               ! Get CFL numbers
               this%CFLp_x =max(this%CFLp_x ,abs(p(n)%vel(1))/this%amr%dx(lvl))
               this%CFLp_y =max(this%CFLp_y ,abs(p(n)%vel(2))/this%amr%dy(lvl))
               this%CFLp_z =max(this%CFLp_z ,abs(p(n)%vel(3))/this%amr%dz(lvl))
               this%CFL_col=max(this%CFL_col,norm2(p(n)%vel)/p(n)%d)
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Global MPI max, scale by dt (10 dt for collision CFL to ensure max CFL of 0.1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp_x, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp_x=this%CFLp_x*dt
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp_y, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp_y=this%CFLp_y*dt
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp_z, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp_z=this%CFLp_z*dt
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFL_col,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFL_col=10.0_WP*this%CFL_col*dt
      ! Return convective CFL and optionally overall CFL
      cflc=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)
      if (present(cfl)) cfl=max(cflc,this%CFL_col)
   end subroutine get_cfl

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Redistribute particles
   subroutine redistribute(this,minlvl,maxlvl,nover)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in), optional :: minlvl,maxlvl,nover
      integer :: lmin,lmax,no
      lmin= 0; if (present(minlvl)) lmin=minlvl
      lmax=-1; if (present(maxlvl)) lmax=maxlvl
      no  = 0; if (present(nover))  no  =nover
      call amrlpt_redistribute(this%pc,lmin,lmax,no)
      call amrlpt_total_np(this%pc,this%np)
   end subroutine redistribute

   !> Fill ghost particle buffer within no cells
   subroutine fill_ghosts(this,no)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: no
      call amrlpt_fill_neighbors(this%pc,no)
   end subroutine fill_ghosts

   !> Fill ghost particles, but only communicate particles within radius of tile boundaries
   subroutine fill_ghosts_radius(this,radius)
      use iso_c_binding, only: c_double
      implicit none
      class(amrlpt), intent(inout) :: this
      real(c_double), intent(in) :: radius
      call amrlpt_fill_neighbors_radius(this%pc,radius)
   end subroutine fill_ghosts_radius

   !> Release ghost particle buffer
   subroutine clear_ghosts(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call amrlpt_clear_neighbors(this%pc)
   end subroutine clear_ghosts

   !> Build explicit pair list within interaction radius rcrit
   subroutine build_neighbor_list(this,rcrit)
      use iso_c_binding, only: c_double
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: rcrit
      call amrlpt_build_neighbor_list(this%pc,real(rcrit,c_double))
   end subroutine build_neighbor_list

   !> Update global particle count
   subroutine get_np(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer(c_int64_t) :: np_c
      call amrlpt_total_np(this%pc,np_c)
      this%np=np_c
   end subroutine get_np

   !> Get pointer to particle array for a given tile
   subroutine get_particles(this,lvl,mfi,p,np)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_c
      call amrlpt_get_particles_mfi(this%pc,lvl,mfi%p,dp,np_c)
      np=np_c
      if (np.gt.0) then
         call c_f_pointer(dp,p,[np])
      else
         nullify(p)
      end if
   end subroutine get_particles

   !> Get pointer to ghost particle array for a given tile (call after fill_ghosts)
   subroutine get_ghosts(this,lvl,mfi,g,ng)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: g
      integer(I8), intent(out) :: ng
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_c
      call amrlpt_get_neighbor_particles_mfi(this%pc,lvl,mfi%p,dp,np_c)
      ng=np_c
      if (ng.gt.0) then
         call c_f_pointer(dp,g,[ng])
      else
         nullify(g)
      end if
   end subroutine get_ghosts

   !> Get combined valid+ghost particle array for a tile (call after fill_ghosts)
   subroutine get_all_particles(this,lvl,mfi,p,np_v,np_all)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np_v,np_all
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_total_c,np_valid_c
      call amrlpt_get_all_particles_mfi(this%pc,lvl,mfi%p,dp,np_total_c,np_valid_c)
      np_v=np_valid_c; np_all=np_total_c
      if (np_total_c.gt.0) then
         call c_f_pointer(dp,p,[np_total_c])
      else
         nullify(p)
      end if
   end subroutine get_all_particles

   !> Get neighbor list pointers for a given tile (call after build_neighbor_list)
   subroutine get_neighbor_list(this,lvl,mfi,off,lst)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      integer(c_int32_t), dimension(:), pointer, intent(out) :: off,lst
      type(c_ptr) :: dp_off,dp_lst
      integer(c_int64_t) :: np_c,ntot_c
      call amrlpt_get_neighbor_list_mfi(this%pc,lvl,mfi%p,dp_off,dp_lst,np_c,ntot_c)
      if (ntot_c.gt.0) then
         call c_f_pointer(dp_off,off,[np_c+1_8])
         call c_f_pointer(dp_lst,lst,[ntot_c])
      else
         nullify(off); nullify(lst)
      end if
   end subroutine get_neighbor_list

   !> Interpolate from cell-centered data to a point
   function interp(this,lvl,pos,arr,comp) result(val)
      implicit none
      class(amrlpt), intent(in) :: this
      integer,  intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: arr
      integer,  intent(in) :: comp
      real(WP) :: val,wx,wy,wz
      integer  :: ii,jj,kk
      ! Get indices and weights
      ii=floor((pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP); wx=(pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP-real(ii,WP)
      jj=floor((pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP); wy=(pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP-real(jj,WP)
      kk=floor((pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP); wz=(pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP-real(kk,WP)
      ! Trilinear interpolation
      val=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*arr(ii  ,jj  ,kk  ,comp) &
      &  +        wx *(1.0_WP-wy)*(1.0_WP-wz)*arr(ii+1,jj  ,kk  ,comp) &
      &  +(1.0_WP-wx)*        wy *(1.0_WP-wz)*arr(ii  ,jj+1,kk  ,comp) &
      &  +        wx *        wy *(1.0_WP-wz)*arr(ii+1,jj+1,kk  ,comp) &
      &  +(1.0_WP-wx)*(1.0_WP-wy)*        wz *arr(ii  ,jj  ,kk+1,comp) &
      &  +        wx *(1.0_WP-wy)*        wz *arr(ii+1,jj  ,kk+1,comp) &
      &  +(1.0_WP-wx)*        wy *        wz *arr(ii  ,jj+1,kk+1,comp) &
      &  +        wx *        wy *        wz *arr(ii+1,jj+1,kk+1,comp)
   end function interp

   !> Interpolate face velocities to a point
   function interp_face_velocities(this,lvl,pos,pU,uc,pV,vc,pW,wc) result(fvel)
      implicit none
      class(amrlpt), intent(in) :: this
      integer,  intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: pU,pV,pW
      integer,  intent(in) :: uc,vc,wc
      real(WP), dimension(3) :: fvel
      real(WP) :: dxi,dyi,dzi,wx,wy,wz,fx,fy,fz
      integer  :: ic,jc,kc,ix,iy,iz
      ! Get cell-centered indices and weights
      dxi=1.0_WP/this%amr%dx(lvl); ic=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
      dyi=1.0_WP/this%amr%dy(lvl); jc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
      dzi=1.0_WP/this%amr%dz(lvl); kc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
      ! Get face-centered indices and weights
      ix=floor((pos(1)-this%amr%xlo)*dxi); fx=(pos(1)-this%amr%xlo)*dxi-real(ix,WP)
      iy=floor((pos(2)-this%amr%ylo)*dyi); fy=(pos(2)-this%amr%ylo)*dyi-real(iy,WP)
      iz=floor((pos(3)-this%amr%zlo)*dzi); fz=(pos(3)-this%amr%zlo)*dzi-real(iz,WP)
      ! U: staggered in x — cell indices jc/wy, kc/wz; face index ix/fx
      fvel(1)=(1.0_WP-fx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix  ,jc  ,kc  ,uc)+fx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix+1,jc  ,kc  ,uc) &
      &      +(1.0_WP-fx)*        wy *(1.0_WP-wz)*pU(ix  ,jc+1,kc  ,uc)+fx*        wy *(1.0_WP-wz)*pU(ix+1,jc+1,kc  ,uc) &
      &      +(1.0_WP-fx)*(1.0_WP-wy)*        wz *pU(ix  ,jc  ,kc+1,uc)+fx*(1.0_WP-wy)*        wz *pU(ix+1,jc  ,kc+1,uc) &
      &      +(1.0_WP-fx)*        wy *        wz *pU(ix  ,jc+1,kc+1,uc)+fx*        wy *        wz *pU(ix+1,jc+1,kc+1,uc)
      ! V: staggered in y — cell indices ic/wx, kc/wz; face index iy/fy
      fvel(2)=(1.0_WP-wx)*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic  ,iy  ,kc  ,vc)+wx*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic+1,iy  ,kc  ,vc) &
      &      +(1.0_WP-wx)*        fy *(1.0_WP-wz)*pV(ic  ,iy+1,kc  ,vc)+wx*        fy *(1.0_WP-wz)*pV(ic+1,iy+1,kc  ,vc) &
      &      +(1.0_WP-wx)*(1.0_WP-fy)*        wz *pV(ic  ,iy  ,kc+1,vc)+wx*(1.0_WP-fy)*        wz *pV(ic+1,iy  ,kc+1,vc) &
      &      +(1.0_WP-wx)*        fy *        wz *pV(ic  ,iy+1,kc+1,vc)+wx*        fy *        wz *pV(ic+1,iy+1,kc+1,vc)
      ! W: staggered in z — cell indices ic/wx, jc/wy; face index iz/fz
      fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic  ,jc  ,iz  ,wc)+wx*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic+1,jc  ,iz  ,wc) &
      &      +(1.0_WP-wx)*        wy *(1.0_WP-fz)*pW(ic  ,jc+1,iz  ,wc)+wx*        wy *(1.0_WP-fz)*pW(ic+1,jc+1,iz  ,wc) &
      &      +(1.0_WP-wx)*(1.0_WP-wy)*        fz *pW(ic  ,jc  ,iz+1,wc)+wx*(1.0_WP-wy)*        fz *pW(ic+1,jc  ,iz+1,wc) &
      &      +(1.0_WP-wx)*        wy *        fz *pW(ic  ,jc+1,iz+1,wc)+wx*        wy *        fz *pW(ic+1,jc+1,iz+1,wc)
   end function interp_face_velocities

   !> Post-process a deposited extensive field into an intensive composite field
   subroutine process_deposit(this,A)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_build,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_sum_downto,amrmfab_interp_from_coarse
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      type(amrex_multifab), dimension(:), allocatable :: tmp
      integer :: lvl
      ! Convert extensive deposits to intensive
      do lvl=0,this%amr%clvl()
         call A%mf(lvl)%mult(1.0_WP/this%amr%cell_vol(lvl),1,A%ncomp,A%ng)
      end do
      ! Allocate scratch MFs for coarse→fine interpolation
      allocate(tmp(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,tmp(lvl),ncomp=A%ncomp,nover=0); call tmp(lvl)%setval(0.0_WP)
      end do
      ! Forward pass (coarse to fine)
      do lvl=0,this%amr%clvl()
         ! SumBoundary: ghost→valid at this level
         call A%syncsum_lvl(lvl)
         ! InterpFromCoarseLevel: propagate coarse deposit to fine tmp
         if (lvl.lt.this%amr%clvl()) then
            call amrmfab_interp_from_coarse(tmp(lvl+1),A%mf(lvl),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl),fgeom=this%amr%geom(lvl+1),scomp=1,ncomp=A%ncomp)
         end if
         ! Sum_fine_to_coarse: fine→coarse ADD
         if (lvl.gt.0) then
            call amrmfab_sum_downto(A%mf(lvl),A%mf(lvl-1),[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),fgeom=this%amr%geom(lvl))
         end if
         ! Add interpolated coarse contribution to this level
         call A%mf(lvl)%add(tmp(lvl),1,1,A%ncomp,0)
      end do
      ! Backward pass: average_down fixes double-counting in covered cells
      do lvl=this%amr%clvl()-1,0,-1
         call A%average_downto(lvl)
      end do
      ! Clean up scratch MFs
      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(tmp(lvl))
      end do
      deallocate(tmp)
   end subroutine process_deposit

   !> Explicit diffusion filter for a cell-centered amrdata field
   subroutine filter(this,A)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      real(WP) :: alpha,alpha_step
      integer  :: nstep,n,nc,lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pA,pFx,pFy,pFz
      real(WP) :: dxi,dyi,dzi

      ! Determine filter coefficient from finest level mesh size
      alpha=max(this%filter_width**2-this%amr%min_meshsize(this%amr%clvl())**2,0.0_WP)/(16.0_WP*log(2.0_WP))
      if (alpha.le.0.0_WP) return

      ! Number of explicit sub-steps for stability (dt < dx^2/6 in 3D)
      nstep=ceiling(6.0_WP*alpha/this%amr%min_meshsize(this%amr%clvl())**2)
      alpha_step=alpha/real(nstep,WP)

      ! Allocate face flux mfabs (no ghost needed for fluxes)
      allocate(Fx(0:this%amr%maxlvl),Fy(0:this%amr%maxlvl),Fz(0:this%amr%maxlvl))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=A%ncomp,nover=0,atface=[.true., .false.,.false.]); call Fx(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.true., .false.]); call Fy(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
      end do

      ! Explicit sub-steps
      do n=1,nstep
         ! Compute diffusive fluxes at each level
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
         ! Average down face fluxes for C/F conservation (finest→coarsest)
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         end do
         ! Apply divergence to update A
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
         ! Restore coarse valid cells covered by fine
         call A%average_down()
         ! Fill ghosts at all levels
         call A%fill(time=0.0_WP)
      end do

      ! Destroy flux mfabs
      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(Fx(lvl))
         call amrex_multifab_destroy(Fy(lvl))
         call amrex_multifab_destroy(Fz(lvl))
      end do
      deallocate(Fx,Fy,Fz)

   end subroutine filter

   !> Allgather all particles whose position lies inside [lo,hi] across all processors
   !> Returns a newly allocated array pgather(1:ngather) containing copies from every rank
   subroutine gather_region(this,lo,hi,pgather,ngather)
      use mpi_f08, only: MPI_ALLGATHER,MPI_ALLGATHERV,MPI_BYTE,MPI_INTEGER
      use amrex_amr_module, only: amrex_mfiter
      use iso_c_binding, only: c_sizeof
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), dimension(3), intent(in) :: lo,hi
      type(part), dimension(:), allocatable, intent(out) :: pgather
      integer(I8), intent(out) :: ngather
      ! Local
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      type(part), dimension(:), allocatable :: pbuf
      type(part) :: dummy
      integer(I8) :: np_,m
      integer :: lvl,nlocal,ntotal,part_bytes,ierr,i
      integer, dimension(:), allocatable :: rcounts,rdisps
      ! Size of one particle in bytes
      part_bytes=int(c_sizeof(dummy))
      ! Pass 1: count local particles in box
      nlocal=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while(mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do m=1,np_
               if (p(m)%pos(1).ge.lo(1).and.p(m)%pos(1).le.hi(1).and. &
               &   p(m)%pos(2).ge.lo(2).and.p(m)%pos(2).le.hi(2).and. &
               &   p(m)%pos(3).ge.lo(3).and.p(m)%pos(3).le.hi(3)) nlocal=nlocal+1
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Pack local matches
      allocate(pbuf(max(nlocal,1)))
      nlocal=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while(mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do m=1,np_
               if (p(m)%pos(1).ge.lo(1).and.p(m)%pos(1).le.hi(1).and. &
               &   p(m)%pos(2).ge.lo(2).and.p(m)%pos(2).le.hi(2).and. &
               &   p(m)%pos(3).ge.lo(3).and.p(m)%pos(3).le.hi(3)) then
                  nlocal=nlocal+1
                  pbuf(nlocal)=p(m)
               end if
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Allgatherv: exchange counts
      allocate(rcounts(this%amr%nproc),rdisps(this%amr%nproc))
      call MPI_ALLGATHER(nlocal,1,MPI_INTEGER,rcounts,1,MPI_INTEGER,this%amr%comm,ierr)
      ntotal=sum(rcounts)
      ! Convert counts and displacements to bytes
      rdisps(1)=0
      do i=2,this%amr%nproc
         rdisps(i)=rdisps(i-1)+rcounts(i-1)*part_bytes
      end do
      rcounts=rcounts*part_bytes
      ! Allgatherv particles
      allocate(pgather(max(ntotal,1)))
      call MPI_ALLGATHERV(pbuf,nlocal*part_bytes,MPI_BYTE,pgather,rcounts,rdisps,MPI_BYTE,this%amr%comm,ierr)
      ngather=int(ntotal,I8)
      deallocate(pbuf,rcounts,rdisps)
   end subroutine gather_region

   !> Append a Fortran part array to the particle container (collective)
   subroutine append(this,pnew,n)
      use messager, only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      type(part), dimension(:), allocatable, target, intent(in) :: pnew
      integer(I8), intent(in) :: n
      type(c_ptr) :: raw
      raw=c_null_ptr
      if (n.gt.0_I8.and.allocated(pnew)) then
         if (int(size(pnew),I8).lt.n) call die('[amrlpt append] pnew array smaller than n')
         raw=c_loc(pnew(1))
      end if
      call amrlpt_append_particles(this%pc,raw,int(n,c_int64_t))
   end subroutine append

   !> Build MFIter from the particle grid's BA/DM
   subroutine mfiter_build(this,lvl,mfi,tiling)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_mfiter_build
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(out) :: mfi
      logical, intent(in), optional :: tiling
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      logical :: use_tiling
      use_tiling=.false.; if (present(tiling)) use_tiling=tiling
      call amrlpt_get_particle_boxarray (this%pc,lvl,ba%p)
      call amrlpt_get_particle_distromap(this%pc,lvl,dm%p)
      call amrex_mfiter_build(mfi,ba,dm,tiling=use_tiling)
   end subroutine mfiter_build

   !> Destroy MFIter
   subroutine mfiter_destroy(this,mfi)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_destroy
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter), intent(inout) :: mfi
      call amrex_mfiter_destroy(mfi)
   end subroutine mfiter_destroy

   !> Set particle-specific BoxArray for a given level
   subroutine set_particle_ba(this,lvl,ba)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      call amrlpt_set_particle_boxarray(this%pc,lvl,ba%p)
   end subroutine set_particle_ba

   !> Set particle-specific DistributionMapping for a given level
   subroutine set_particle_dm(this,lvl,dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_distromap), intent(in) :: dm
      call amrlpt_set_particle_distromap(this%pc,lvl,dm%p)
   end subroutine set_particle_dm

   !> Get particle-specific DistributionMapping for a given level
   function get_particle_dm(this,lvl) result(dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrlpt), intent(in) :: this
      integer, intent(in) :: lvl
      type(amrex_distromap) :: dm
      call amrlpt_get_particle_distromap(this%pc,lvl,dm%p)
      dm%owner=.false.
   end function get_particle_dm

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Compute particle statistics: np, d/vel min/max/mean/var, Vp_tot, and VF field stats
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(amrex_box)    :: bx
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), pointer :: pVF
      integer(I8) :: np_
      integer :: lvl,n,i,j,k,ierr
      real(WP) :: d_sum,d_sq,vx_sum,vx_sq,vy_sum,vy_sq,vz_sum,vz_sq,inv_np
      ! Init per-rank accumulators
      this%np=0; this%Vp_tot=0.0_WP
      this%dmin=huge(1.0_WP); this%dmax=-huge(1.0_WP);  d_sum=0.0_WP;  d_sq=0.0_WP
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); vx_sum=0.0_WP; vx_sq=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); vy_sum=0.0_WP; vy_sq=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); vz_sum=0.0_WP; vz_sq=0.0_WP
      ! Loop over all AMR levels and tiles
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               this%np=this%np+1
               this%Vp_tot=this%Vp_tot+Pi/6.0_WP*p(n)%d**3
               this%dmin=min(this%dmin,p(n)%d);      this%dmax=max(this%dmax,p(n)%d);       d_sum= d_sum+p(n)%d;       d_sq= d_sq+p(n)%d**2
               this%Umin=min(this%Umin,p(n)%vel(1)); this%Umax=max(this%Umax,p(n)%vel(1)); vx_sum=vx_sum+p(n)%vel(1); vx_sq=vx_sq+p(n)%vel(1)**2
               this%Vmin=min(this%Vmin,p(n)%vel(2)); this%Vmax=max(this%Vmax,p(n)%vel(2)); vy_sum=vy_sum+p(n)%vel(2); vy_sq=vy_sq+p(n)%vel(2)**2
               this%Wmin=min(this%Wmin,p(n)%vel(3)); this%Wmax=max(this%Wmax,p(n)%vel(3)); vz_sum=vz_sum+p(n)%vel(3); vz_sq=vz_sq+p(n)%vel(3)**2
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Compute per-rank particle distribution stats
      this%np_loc=this%np; this%np_min=this%np; this%np_max=this%np; this%np_eff=0.0_WP
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_min,1,MPI_INTEGER8,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_max,1,MPI_INTEGER8,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np,    1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
      if (this%np_max.gt.0) this%np_eff=real(this%np,WP)/real(this%np_max,WP)/real(this%amr%nproc,WP)
      ! Reduce other statistics
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_tot,1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sum,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sq,       1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sum,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sq,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sum,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sq,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sum,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sq,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      ! Derive mean and variance from single-pass accumulation
      if (this%np.gt.0) then
         inv_np=1.0_WP/real(this%np,WP)
         this%dmean= d_sum*inv_np; this%dvar=max(0.0_WP, d_sq*inv_np-this%dmean**2)
         this%Umean=vx_sum*inv_np; this%Uvar=max(0.0_WP,vx_sq*inv_np-this%Umean**2)
         this%Vmean=vy_sum*inv_np; this%Vvar=max(0.0_WP,vy_sq*inv_np-this%Vmean**2)
         this%Wmean=vz_sum*inv_np; this%Wvar=max(0.0_WP,vz_sq*inv_np-this%Wmean**2)
      end if
      ! Volume fraction statistics
      vf_stats_block: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface,  only: amrmask_make_fine
         type(amrex_imultifab) :: mask
         integer,  dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP) :: var_sum
         ! Initialize stats
         this%VFmean=this%VF%get_sum(lvl=0)/real(this%amr%nx*this%amr%ny*this%amr%nz,WP)
         var_sum=0.0_WP; this%VFmin=huge(1.0_WP); this%VFmax=-huge(1.0_WP)
         ! Loop over levels
         do lvl=0,this%amr%clvl()
            ! Build fine mask
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            ! Loop over tiles
            call amrex_mfiter_build(mfi,this%VF%mf(lvl))
            do while (mfi%next())
               ! Get pointer to data
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over tile
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then
                     if (pMask(i,j,k,1).eq.0) cycle
                  end if
                  ! Accumulate statistics
                  this%VFmin=min(this%VFmin,pVF(i,j,k,1)); this%VFmax=max(this%VFmax,pVF(i,j,k,1))
                  var_sum=var_sum+(pVF(i,j,k,1)-this%VFmean)**2*this%amr%cell_vol(lvl)
               end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,var_sum,   1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%VFvar=max(0.0_WP,var_sum/((this%amr%xhi-this%amr%xlo)*(this%amr%yhi-this%amr%ylo)*(this%amr%zhi-this%amr%zlo)))
      end block vf_stats_block

      ! Reduce local accumulators, publish to monitoring fields, then zero for next interval
      reduce_local_stats: block
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_out_loc,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_out_loc,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_new_loc,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_new_loc,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%np_out=this%np_out_loc; this%Vp_out=this%Vp_out_loc
         this%np_new=this%np_new_loc; this%Vp_new=this%Vp_new_loc
         this%np_out_loc=0; this%Vp_out_loc=0.0_WP
         this%np_new_loc=0; this%Vp_new_loc=0.0_WP
      end block reduce_local_stats

      ! Compute timing statistics across ranks
      call this%tmr_coll%get_stats()
      call this%tmr_coll_%get_stats()
      call this%tmr_fill%get_stats()
      call this%tmr_nbl%get_stats()
      call this%tmr_step%get_stats()
      call this%tmr_step_%get_stats()
      call this%tmr_vf%get_stats()
      call this%tmr_src%get_stats()

      ! Reset timers
      call this%tmr_coll%reset()
      call this%tmr_coll_%reset()
      call this%tmr_fill%reset()
      call this%tmr_nbl%reset()
      call this%tmr_step%reset()
      call this%tmr_step_%reset()
      call this%tmr_vf%reset()
      call this%tmr_src%reset()

   end subroutine get_info

   !> Print solver info
   subroutine print(this)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(amrlpt), intent(in) :: this
      character(len=str_long) :: message
      call log('AMR Lagrangian particle solver: '//trim(this%name))
      write(message,'("  Particle density  : ",ES12.5)') this%rho
      call log(trim(message))
      write(message,'("  Drag model        : ",a)') trim(this%drag_model)
      call log(trim(message))
      write(message,'("  ODE substeps      : ",i0)') this%nstep
      call log(trim(message))
      call log('  Grid: '//trim(this%amr%name))
   end subroutine print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Write AMReX checkpoint for particles to dirname
   subroutine write(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_write(this%pc,trim(dirname)//c_null_char,1_c_int)
   end subroutine write

   !> Read AMReX checkpoint for particles from dirname
   subroutine read(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_read(this%pc,trim(dirname)//c_null_char)
      call amrlpt_redistribute(this%pc,0,-1,0)
      call amrlpt_total_np(this%pc,this%np)
      call this%update_VF()
   end subroutine read

end module amrlpt_class
