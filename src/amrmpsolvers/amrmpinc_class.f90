!> AMR collocated incompressible multiphase solver class
!> Inherits from amrvof_class
module amrmpinc_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrmg_class,      only: amrmg
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_interp_face_divfree
   use amrvof_class,     only: amrvof,VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
   implicit none
   private

   ! Expose type
   public :: amrmpinc,VFlo,VFhi,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER

   !> AMR collocated incompressible multiphase solver type
   type, extends(amrvof) :: amrmpinc

      ! User-configurable callbacks
      procedure(mpinc_init_iface),    pointer, pass :: user_mpinc_init   =>null()
      procedure(mpinc_tagging_iface), pointer, pass :: user_mpinc_tagging=>null()
      procedure(mpinc_bc_iface),      pointer, pass :: user_mpinc_bc     =>null()

      ! Flow data
      type(amrdata) :: U,V,W            !< Face velocities
      type(amrdata) :: Uold,Vold,Wold   !< Old face velocities
      type(amrdata) :: UVW              !< Cell-centered velocity
      type(amrdata) :: UVWold           !< Old cell-centered velocity
      type(amrdata) :: P                !< Pressure (cell-centered)
      type(amrdata) :: div              !< Divergence (cell-centered)

      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rhoL=1.0_WP           !< Constant liquid density
      real(WP) :: rhoG=1.0_WP           !< Constant gas density
      type(amrdata) :: visc             !< Variable dynamic viscosity
      real(WP) :: sigma=0.0_WP          !< Surface tension coefficient

      ! Monitoring quantities
      real(WP) :: Umax=0.0_WP           !< Max U velocity
      real(WP) :: Vmax=0.0_WP           !< Max V velocity
      real(WP) :: Wmax=0.0_WP           !< Max W velocity
      real(WP) :: Pmax=0.0_WP           !< Max pressure
      real(WP) :: divmax=0.0_WP         !< Maximum divergence
      real(WP) :: rhoUint=0.0_WP        !< Integral of rho*U
      real(WP) :: rhoVint=0.0_WP        !< Integral of rho*V
      real(WP) :: rhoWint=0.0_WP        !< Integral of rho*W
      real(WP) :: rhoKint=0.0_WP        !< Integral of rho*K
      real(WP) :: CFLc_x=0.0_WP, CFLc_y=0.0_WP, CFLc_z=0.0_WP  !< Convective CFL
      real(WP) :: CFLv_x=0.0_WP, CFLv_y=0.0_WP, CFLv_z=0.0_WP  !< Viscous CFL
      real(WP) :: CFLst=0.0_WP          !< Surface tension CFL
      real(WP) :: CFL=0.0_WP            !< Maximum CFL

      ! Velocity interpolation method (for C/F interface fills)
      integer :: interp_vel=amrex_interp_face_divfree

      ! Cost for mixed cells in load balancing (1.0 = same as pure, higher = more expensive)
      real(WP) :: SLcost=100.0_WP

   contains
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Lifecycle callbacks
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      procedure :: tagging
      procedure :: get_cost
      ! Staggered velocity fills
      procedure :: fill_velocity_lvl          !< Fill face velocity ghosts at single level
      procedure :: fill_velocity              !< Fill face velocity ghosts on all levels
      procedure :: fill_velocity_from_coarse  !< Fill face velocity from coarse
      procedure :: sync_velocity_lvl          !< Sync face velocity ghosts at single level
      procedure :: sync_velocity              !< Sync face velocity ghosts on all levels
      procedure :: fill_velocity_mfab         !< Fill dest MultiFabs for regridding
      procedure :: average_down_velocity      !< Average down face velocity for C/F consistency
      procedure :: average_down_velocity_to   !< Average down face velocity for single level
      ! Utilities
      procedure :: store_old                  !< Store current state to old state
      procedure :: get_div                    !< Compute divergence (assumes velocity ghosts filled)
      procedure :: interp_vel_to_face         !< Interpolate cell-centered velocity to face
      procedure :: prepare_psolver            !< Prepare pressure solver with new densities
      procedure :: correct_both_velocities    !< Correct both face and cell-centered velocities
      procedure :: add_surface_tension        !< Add surface tension increment to both velocities
      procedure, private :: apply_face_fluxes !< Apply pre-built face fluxes to both face and cell-centered velocities
      ! Physics procedures
      procedure :: get_dUVWdt                 !< Compute d(UVW)/dt
      procedure :: get_dmomdt                 !< Compute d(mom)/dt
      procedure :: add_vreman                 !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                    !< Compute CFL numbers
      procedure :: correct_outflow            !< Correct outflow for global mass conservation
      ! Print solver info
      procedure :: get_info
      procedure :: print => amrmpinc_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrmpinc

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine mpinc_init_iface(solver,lvl,time,ba,dm)
         import :: amrmpinc,WP,amrex_boxarray,amrex_distromap
         class(amrmpinc), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine mpinc_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine mpinc_tagging_iface(solver,lvl,time,tags)
         import :: amrmpinc,c_ptr,WP
         class(amrmpinc), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine mpinc_tagging_iface
   end interface

   !> Abstract interface for user-provided velocity BC callback
   !> Called for ext_dir faces; user fills the boundary box with their own values
   abstract interface
      subroutine mpinc_bc_iface(solver,lvl,time,face,bx,comp,p)
         import :: amrmpinc,amrex_box,WP
         class(amrmpinc), intent(in) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< 'U', 'V', or 'W'
         real(WP), dimension(:,:,:,:), pointer, intent(inout) :: p
      end subroutine mpinc_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrcinc type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrmpinc_on_init(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_mpinc_init)) call this%user_mpinc_init(lvl,time,ba,dm)
   end subroutine amrmpinc_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrmpinc_on_coarse(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrmpinc_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrmpinc_on_remake(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrmpinc_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrmpinc_on_clear(ctx,lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrmpinc_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrmpinc_tagging(ctx,lvl,time,tags)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_mpinc_tagging)) call this%user_mpinc_tagging(lvl,time,tags)
   end subroutine amrmpinc_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrmpinc_postregrid(ctx,lbase,time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrmpinc_postregrid

   !> Dispatch cost: calls type-bound method
   subroutine amrmpinc_get_cost(ctx,lvl,nboxes,costs,ba)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrmpinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrmpinc_get_cost

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the incompressible solver
   subroutine initialize(this,amr,name)
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrmg_class,      only: amrmg_varcoef
      implicit none
      class(amrmpinc), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Initialize amrvof parent without callback registration
      this%amrvof%skip_registration=.true.
      call this%amrvof%initialize(amr,name)

      ! Initialize staggered velocity
      call this%U%initialize(amr,name='U',ncomp=1,ng=this%nover,nodal=[.true. ,.false.,.false.]); this%U%parent=>this
      call this%V%initialize(amr,name='V',ncomp=1,ng=this%nover,nodal=[.false.,.true. ,.false.]); this%V%parent=>this
      call this%W%initialize(amr,name='W',ncomp=1,ng=this%nover,nodal=[.false.,.false.,.true. ]); this%W%parent=>this
      call this%Uold%initialize(amr,name='Uold',ncomp=1,ng=this%nover,nodal=[.true. ,.false.,.false.]); this%Uold%parent=>this
      call this%Vold%initialize(amr,name='Vold',ncomp=1,ng=this%nover,nodal=[.false.,.true. ,.false.]); this%Vold%parent=>this
      call this%Wold%initialize(amr,name='Wold',ncomp=1,ng=this%nover,nodal=[.false.,.false.,.true. ]); this%Wold%parent=>this

      ! Set velocity fillbc callbacks to shared internal handler
      this%U%fillbc=>velocity_fillbc
      this%V%fillbc=>velocity_fillbc
      this%W%fillbc=>velocity_fillbc

      ! Initialize collocated velocities
      call this%UVW%initialize   (amr,name='UVW'   ,ncomp=3,ng=this%nover); this%UVW%parent   =>this
      call this%UVWold%initialize(amr,name='UVWold',ncomp=3,ng=this%nover); this%UVWold%parent=>this

      ! Set UVW fillbc callback to internal handler
      this%UVW%fillbc=>UVW_fillbc

      ! Initialize pressure with Neumann BCs
      call this%P%initialize(amr,name='P',ncomp=1,ng=this%nover); this%P%parent=>this
      if (.not.amr%xper) then; this%P%lo_bc(1,1)=amrex_bc_foextrap; this%P%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%P%lo_bc(2,1)=amrex_bc_foextrap; this%P%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%P%lo_bc(3,1)=amrex_bc_foextrap; this%P%hi_bc(3,1)=amrex_bc_foextrap; end if

      ! Initialize divergence (no ghosts needed)
      call this%div%initialize(amr,name='div',ncomp=1,ng=0); this%div%parent=>this

      ! Initialize viscosity with Neumann BCs
      call this%visc%initialize(amr,name='visc',ncomp=1,ng=this%nover); this%visc%parent=>this
      if (.not.amr%xper) then; this%visc%lo_bc(1,1)=amrex_bc_foextrap; this%visc%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%visc%lo_bc(2,1)=amrex_bc_foextrap; this%visc%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%visc%lo_bc(3,1)=amrex_bc_foextrap; this%visc%hi_bc(3,1)=amrex_bc_foextrap; end if

      ! Initialize pressure solver
      call this%psolver%initialize(amr,type=amrmg_varcoef); this%psolver%beta=-1.0_WP

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      select type (this)
       type is (amrmpinc)
         call this%amr%add_on_init   (amrmpinc_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrmpinc_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrmpinc_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrmpinc_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrmpinc_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrmpinc_postregrid,c_loc(this))
         call this%amr%set_get_cost  (amrmpinc_get_cost,  c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the incompressible solver
   subroutine finalize(this)
      implicit none
      class(amrmpinc), intent(inout) :: this
      call this%UVW%finalize()
      call this%UVWold%finalize()
      call this%U%finalize()
      call this%V%finalize()
      call this%W%finalize()
      call this%Uold%finalize()
      call this%Vold%finalize()
      call this%Wold%finalize()
      call this%P%finalize()
      call this%div%finalize()
      call this%visc%finalize()
      call this%psolver%finalize()
      nullify(this%user_mpinc_init)
      nullify(this%user_mpinc_tagging)
      nullify(this%user_mpinc_bc)
      call this%amrvof%finalize()
   end subroutine finalize

   ! ============================================================================
   ! LIFECYCLE CALLBACKS
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles VF/VFold + CL/CG/PLIC multifabs
      call this%amrvof%on_init(lvl,time,ba,dm)
      ! Reset level layouts
      call this%UVW%reset_level(lvl,ba,dm)
      call this%UVWold%reset_level(lvl,ba,dm)
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%div%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      ! Set to zero
      call this%U%setval(val=0.0_WP,lvl=lvl)
      call this%V%setval(val=0.0_WP,lvl=lvl)
      call this%W%setval(val=0.0_WP,lvl=lvl)
      call this%Uold%setval(val=0.0_WP,lvl=lvl)
      call this%Vold%setval(val=0.0_WP,lvl=lvl)
      call this%Wold%setval(val=0.0_WP,lvl=lvl)
      call this%UVW%setval(val=0.0_WP,lvl=lvl)
      call this%UVWold%setval(val=0.0_WP,lvl=lvl)
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%div%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using divergence-free interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles VF interpolation and CL/CG/PLIC rebuild
      call this%amrvof%on_coarse(lvl,time,ba,dm)
      ! Collocated velocity
      call this%UVW%on_coarse(lvl,time,ba,dm)
      call this%UVWold%reset_level(lvl,ba,dm)
      ! Face velocity: allocate then fill with divergence-free interpolation
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%fill_velocity_from_coarse(lvl,time)
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      ! Pressure uses default on_coarse
      call this%P%on_coarse(lvl,time,ba,dm)
      ! Divergence and viscosity just need to be reset
      call this%div%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles VF remake and CL/CG/PLIC parallel_copy
      call this%amrvof%on_remake(lvl,time,ba,dm)
      ! Collocated velocity
      call this%UVW%on_remake(lvl,time,ba,dm)
      call this%UVWold%reset_level(lvl,ba,dm)
      ! Face velocity: fill with div-free interpolation
      face_vel_remake: block
         use amrex_amr_module, only: amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
         type(amrex_multifab) :: Utmp,Vtmp,Wtmp
         ! Build temp MultiFabs with new layout (0 ghost cells for FillPatch)
         call amrex_multifab_build(Utmp,ba,dm,1,0,this%U%nodal)
         call amrex_multifab_build(Vtmp,ba,dm,1,0,this%V%nodal)
         call amrex_multifab_build(Wtmp,ba,dm,1,0,this%W%nodal)
         ! Fill temps from old data via coupled FillPatch
         call this%fill_velocity_mfab(Utmp,Vtmp,Wtmp,lvl,time)
         ! Reset levels and copy from temps
         call this%U%reset_level(lvl,ba,dm)
         call this%V%reset_level(lvl,ba,dm)
         call this%W%reset_level(lvl,ba,dm)
         call this%U%mf(lvl)%copy(Utmp,1,1,1,0)
         call this%V%mf(lvl)%copy(Vtmp,1,1,1,0)
         call this%W%mf(lvl)%copy(Wtmp,1,1,1,0)
         ! Destroy temps
         call amrex_multifab_destroy(Utmp)
         call amrex_multifab_destroy(Vtmp)
         call amrex_multifab_destroy(Wtmp)
      end block face_vel_remake
      ! Reset old velocities
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      ! Pressure remake
      call this%P%on_remake(lvl,time,ba,dm)
      ! Divergence and viscosity just need geometry
      call this%div%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles VF/VFold + CL/CG/PLIC
      call this%amrvof%on_clear(lvl)
      ! Flow solver data
      call this%UVW%clear_level(lvl)
      call this%UVWold%clear_level(lvl)
      call this%U%clear_level(lvl)
      call this%V%clear_level(lvl)
      call this%W%clear_level(lvl)
      call this%Uold%clear_level(lvl)
      call this%Vold%clear_level(lvl)
      call this%Wold%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%div%clear_level(lvl)
      call this%visc%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles VF average-down + fill
      call this%amrvof%post_regrid(lbase,time)
      ! Average down UVW, face velocities, and pressure
      call this%UVW%average_down(lbase)
      call this%average_down_velocity(lbase)
      call this%P%average_down(lbase)
      ! Fill ghosts
      call this%UVW%fill(time,lbase)
      call this%fill_velocity(time,lbase)
      call this%P%fill(time,lbase)
   end subroutine post_regrid

   !> Tag cells near interface with regrid_buffer layer growth
   subroutine tagging(this,lvl,time,tags)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      ! Parent handles VF tagging
      call this%amrvof%tagging(lvl,time,tags)
   end subroutine tagging

   !> Estimate per-box costs for load balancing
   !> Cost based on number of mixed cells vs pure cells
   subroutine get_cost(this,lvl,nboxes,costs,ba)
      use iso_c_binding, only: c_associated
      use amrex_amr_module, only: amrex_boxarray,amrex_box,amrex_mfiter,&
      &                           amrex_mfiter_build,amrex_mfiter_destroy,amrex_intersection
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: old_bx,new_bx,isect
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer :: n,i,j,k,ierr
      ! Guard: if VF data doesn't exist yet, return uniform costs
      if (.not.allocated(this%VF%mf)) then; costs=1.0_WP; return; end if
      if (.not.c_associated(this%VF%mf(lvl)%p)) then; costs=1.0_WP; return; end if
      ! Coarser levels are never mixed
      if (lvl.lt.this%amr%clvl()) then
         do n=1,nboxes
            new_bx=ba%get_box(n-1)
            costs(n)=real(new_bx%numpts(),WP)
         end do
         return
      end if
      ! At finest level, count mixed cells per new box from local old data
      costs=0.0_WP
      call amrex_mfiter_build(mfi,this%VF%mf(lvl),tiling=.false.)
      do while (mfi%next())
         old_bx=mfi%tilebox()
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         do n=1,nboxes
            new_bx=ba%get_box(n-1)  ! 0-indexed
            if (.not.old_bx%intersects(new_bx)) cycle
            isect=amrex_intersection(old_bx,new_bx)
            do k=isect%lo(3),isect%hi(3); do j=isect%lo(2),isect%hi(2); do i=isect%lo(1),isect%hi(1)
               if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
                  costs(n)=costs(n)+this%SLcost
               else
                  costs(n)=costs(n)+1.0_WP
               end if
            end do; end do; end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! Allreduce: sum partial mixed-cell counts across ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine get_cost

   ! ============================================================================
   ! Staggered velocity fills
   ! ============================================================================

   !> Average down MAC velocity for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles face-centered averaging correctly
   subroutine average_down_velocity_to(this,lvl)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%average_downto(lvl)
      call this%V%average_downto(lvl)
      call this%W%average_downto(lvl)
   end subroutine average_down_velocity_to

   !> Average down MAC velocity from finest to lbase
   !> Simply calls average_down_velocity_to in a loop
   !> @param lbase Optional: lowest level to average down to (default 0)
   subroutine average_down_velocity(this,lbase)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=this%amr%clvl()-1,lb,-1
         call this%average_down_velocity_to(lvl)
      end do
   end subroutine average_down_velocity

   !> Fill velocity ghost cells at a single level using divergence-free interpolation
   subroutine fill_velocity_lvl(this,lvl,time)
      use iso_c_binding, only: c_loc,c_funloc,c_funptr,c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single,amrmfab_fillpatch_two_faces
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrmpinc), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u,ctx_v,ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr(3),lo_bc(9),hi_bc(9)
      real(WP) :: t_old,t_new

      ! Get contexts for each velocity component
      ctx_u=c_loc(this%U);ctx_v=c_loc(this%V);ctx_w=c_loc(this%W)
      bc_dispatch=c_funloc(amrdata_fillbc)
      t_old=time-1.0e200_WP
      t_new=time

      ! Store current work level
      this%U%fill_lvl_cache=lvl
      this%V%fill_lvl_cache=lvl
      this%W%fill_lvl_cache=lvl

      if (lvl.eq.0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(this%U%mf(0),time,this%U%mf(0),time,this%U%mf(0),this%amr%geom(0),ctx_u,bc_dispatch,time,1,1,1)
         call amrmfab_fillpatch_single(this%V%mf(0),time,this%V%mf(0),time,this%V%mf(0),this%amr%geom(0),ctx_v,bc_dispatch,time,1,1,1)
         call amrmfab_fillpatch_single(this%W%mf(0),time,this%W%mf(0),time,this%W%mf(0),this%amr%geom(0),ctx_w,bc_dispatch,time,1,1,1)
      else
         ! Build combined BC array: [U_x,U_y,U_z, V_x,V_y,V_z, W_x,W_y,W_z]
         lo_bc(1:3)=this%U%lo_bc(:,1)
         lo_bc(4:6)=this%V%lo_bc(:,1)
         lo_bc(7:9)=this%W%lo_bc(:,1)
         hi_bc(1:3)=this%U%hi_bc(:,1)
         hi_bc(4:6)=this%V%hi_bc(:,1)
         hi_bc(7:9)=this%W%hi_bc(:,1)
         rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)]
         ! Call 3-component divfree FillPatch from two levels
         call amrmfab_fillpatch_two_faces( &
         &   this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl),time, &
         &   t_old,this%U%mf(lvl-1),this%V%mf(lvl-1),this%W%mf(lvl-1), &
         &   t_new,this%U%mf(lvl-1),this%V%mf(lvl-1),this%W%mf(lvl-1), &
         &   this%amr%geom(lvl-1), &
         &   t_old,this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl), &
         &   t_new,this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl), &
         &   this%amr%geom(lvl), &
         &   ctx_u,ctx_v,ctx_w,bc_dispatch,bc_dispatch,bc_dispatch, &
         &   1,1,1,rr,this%interp_vel,lo_bc,hi_bc)
      end if
      ! Reconcile shared face values at box boundaries
      call this%U%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%V%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%W%mf(lvl)%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_lvl

   !> Fill velocity ghost cells on all levels
   subroutine fill_velocity(this,time,lbase)
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: time
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=lb,this%amr%clvl()
         call this%fill_velocity_lvl(lvl,time)
      end do
   end subroutine fill_velocity

   !> Fill new fine level velocity from coarse using divergence-free interpolation
   !> Used during MakeNewLevelFromCoarse (creation of new fine levels)
   subroutine fill_velocity_from_coarse(this,lvl,time)
      use iso_c_binding, only: c_loc,c_funloc,c_funptr,c_ptr
      use amrex_interface, only: amrmfab_fillcoarsepatch_faces
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrmpinc), target, intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u,ctx_v,ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr(3),lo_bc(9),hi_bc(9)
      ! Get contexts for each velocity component
      ctx_u=c_loc(this%U);ctx_v=c_loc(this%V);ctx_w=c_loc(this%W)
      bc_dispatch=c_funloc(amrdata_fillbc)
      ! Build combined BC array: [U_x,U_y,U_z, V_x,V_y,V_z, W_x,W_y,W_z]
      lo_bc(1:3)=this%U%lo_bc(:,1)
      lo_bc(4:6)=this%V%lo_bc(:,1)
      lo_bc(7:9)=this%W%lo_bc(:,1)
      hi_bc(1:3)=this%U%hi_bc(:,1)
      hi_bc(4:6)=this%V%hi_bc(:,1)
      hi_bc(7:9)=this%W%hi_bc(:,1)
      rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)]
      ! Call 3-component divfree FillCoarsePatch
      call amrmfab_fillcoarsepatch_faces( &
      &   this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl),time, &
      &   this%U%mf(lvl-1),this%V%mf(lvl-1),this%W%mf(lvl-1), &
      &   this%amr%geom(lvl-1),this%amr%geom(lvl), &
      &   ctx_u,ctx_v,ctx_w,bc_dispatch,bc_dispatch,bc_dispatch, &
      &   1,1,1,rr,this%interp_vel,lo_bc,hi_bc)
      ! Reconcile shared face values at box boundaries
      call this%U%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%V%mf(lvl)%override_sync(this%amr%geom(lvl))
      call this%W%mf(lvl)%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_from_coarse

   !> Sync velocity ghost cells at a single level (lightweight, no C/F interpolation)
   subroutine sync_velocity_lvl(this,lvl)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%sync_lvl(lvl)
      call this%V%sync_lvl(lvl)
      call this%W%sync_lvl(lvl)
   end subroutine sync_velocity_lvl

   !> Sync velocity ghost cells on all levels
   subroutine sync_velocity(this,lbase)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=lb,this%amr%clvl()
         call this%sync_velocity_lvl(lvl)
      end do
   end subroutine sync_velocity

   !> Fill destination MultiFabs with velocity using divergence-free interpolation
   !> Used during regridding (on_remake) to fill new layout MultiFabs
   subroutine fill_velocity_mfab(this,Udest,Vdest,Wdest,lvl,time)
      use iso_c_binding, only: c_loc,c_funloc,c_funptr,c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single,amrmfab_fillpatch_two_faces
      use amrex_amr_module, only: amrex_multifab
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrmpinc), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: Udest,Vdest,Wdest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr) :: ctx_u,ctx_v,ctx_w
      type(c_funptr) :: bc_dispatch
      integer :: rr(3),lo_bc(9),hi_bc(9)
      real(WP) :: t_old,t_new

      ! Get contexts for each velocity component
      ctx_u=c_loc(this%U); ctx_v=c_loc(this%V); ctx_w=c_loc(this%W)
      bc_dispatch=c_funloc(amrdata_fillbc)
      t_old=time-1.0e200_WP
      t_new=time

      ! Store current work level
      this%U%fill_lvl_cache=lvl
      this%V%fill_lvl_cache=lvl
      this%W%fill_lvl_cache=lvl

      if (lvl .eq. 0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(Udest,t_old,this%U%mf(0),t_new,this%U%mf(0),this%amr%geom(0),ctx_u,bc_dispatch,time,1,1,1)
         call amrmfab_fillpatch_single(Vdest,t_old,this%V%mf(0),t_new,this%V%mf(0),this%amr%geom(0),ctx_v,bc_dispatch,time,1,1,1)
         call amrmfab_fillpatch_single(Wdest,t_old,this%W%mf(0),t_new,this%W%mf(0),this%amr%geom(0),ctx_w,bc_dispatch,time,1,1,1)
      else
         ! Build combined BC array
         lo_bc(1:3)=this%U%lo_bc(:,1)
         lo_bc(4:6)=this%V%lo_bc(:,1)
         lo_bc(7:9)=this%W%lo_bc(:,1)
         hi_bc(1:3)=this%U%hi_bc(:,1)
         hi_bc(4:6)=this%V%hi_bc(:,1)
         hi_bc(7:9)=this%W%hi_bc(:,1)
         rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)]
         ! Call 3-component divfree FillPatch
         call amrmfab_fillpatch_two_faces( &
         &   Udest,Vdest,Wdest,time, &
         &   t_old,this%U%mf(lvl-1),this%V%mf(lvl-1),this%W%mf(lvl-1), &
         &   t_new,this%U%mf(lvl-1),this%V%mf(lvl-1),this%W%mf(lvl-1), &
         &   this%amr%geom(lvl-1), &
         &   t_old,this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl), &
         &   t_new,this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl), &
         &   this%amr%geom(lvl), &
         &   ctx_u,ctx_v,ctx_w,bc_dispatch,bc_dispatch,bc_dispatch, &
         &   1,1,1,rr,this%interp_vel,lo_bc,hi_bc)
      end if
      ! Reconcile shared face values at box boundaries
      call Udest%override_sync(this%amr%geom(lvl))
      call Vdest%override_sync(this%amr%geom(lvl))
      call Wdest%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_mfab

   !> Velocity boundary condition callback (shared by U, V, W)
   !> Handles staggering-aware BC fills for face-centered velocity data
   !> - ext_dir: calls user_bc callback for user-controlled values
   !> - foextrap: copies from interior (Neumann, zero gradient)
   !> - reflect_even/odd: symmetry/anti-symmetry
   subroutine velocity_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_geometry,&
      &                           amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_even,amrex_bc_reflect_odd
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp,ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      class(amrmpinc), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: dlo(3),dhi(3),flo(3),fhi(3)
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: lvl
      character(len=1) :: comp

      ! Get point to solver
      select type (s=>this%parent)
       class is (amrmpinc)
         solver=>s
      end select

      ! Get domain bounds, component ID, and level
      dlo=geom%domain%lo
      dhi=geom%domain%hi
      comp=this%name(1:1)  ! 'U', 'V', or 'W'
      lvl=this%fill_lvl_cache

      ! Compute staggered domain bounds for this component
      ! Cell domain: [dlo, dhi]. Face domains extend by 1 in staggered direction.
      flo=dlo; fhi=dhi
      select case (comp)
       case ('U'); fhi(1)=dhi(1)+1
       case ('V'); fhi(2)=dhi(2)+1
       case ('W'); fhi(3)=dhi(3)+1
      end select

      ! Loop over FABs
      call amrex_mfiter_build(mfi,mf,tiling=.false.)
      do while (mfi%next())
         p=>mf%dataptr(mfi)
         ilo=lbound(p,1); ihi=ubound(p,1)
         jlo=lbound(p,2); jhi=ubound(p,2)
         klo=lbound(p,3); khi=ubound(p,3)
         ! Skip if FAB entirely within staggered domain
         if (ilo.ge.flo(1).and.ihi.le.fhi(1).and.jlo.ge.flo(2).and.jhi.le.fhi(2).and.klo.ge.flo(3).and.khi.le.fhi(3)) cycle
         ! X-LOW BOUNDARY
         if (.not.solver%amr%xper.and.ilo.lt.flo(1)) call apply_vel_bc(bnd=flo(1),bctype=this%lo_bc(1,1),face=1)
         ! X-HIGH BOUNDARY
         if (.not.solver%amr%xper.and.ihi.gt.fhi(1)) call apply_vel_bc(bnd=fhi(1),bctype=this%hi_bc(1,1),face=2)
         ! Y-LOW BOUNDARY
         if (.not.solver%amr%yper.and.jlo.lt.flo(2)) call apply_vel_bc(bnd=flo(2),bctype=this%lo_bc(2,1),face=3)
         ! Y-HIGH BOUNDARY
         if (.not.solver%amr%yper.and.jhi.gt.fhi(2)) call apply_vel_bc(bnd=fhi(2),bctype=this%hi_bc(2,1),face=4)
         ! Z-LOW BOUNDARY
         if (.not.solver%amr%zper.and.klo.lt.flo(3)) call apply_vel_bc(bnd=flo(3),bctype=this%lo_bc(3,1),face=5)
         ! Z-HIGH BOUNDARY
         if (.not.solver%amr%zper.and.khi.gt.fhi(3)) call apply_vel_bc(bnd=fhi(3),bctype=this%hi_bc(3,1),face=6)
      end do
      call amrex_mfiter_destroy(mfi)

   contains

      !> Apply BC at boundary face
      !> face: 1=xlo, 2=xhi, 3=ylo, 4=yhi, 5=zlo, 6=zhi
      !> Derives direction and side from face
      !> For NORMAL component (e.g., U in x): fills boundary face + ghosts
      !> For TANGENT component (e.g., V in x): fills ghosts only
      subroutine apply_vel_bc(bnd,bctype,face)
         implicit none
         integer, intent(in) :: bnd,bctype,face
         type(amrex_box) :: bc_bx
         integer :: i,j,k,dir,fill_edge,src_from,toff
         integer :: slo(3),shi(3)
         logical :: is_normal,is_lo
         ! Derive direction and side from face
         dir=(face+1)/2
         is_lo=mod(face,2).eq.1
         ! Check if this is a normal component (U in x, V in y, W in z)
         is_normal=(comp.eq.'U'.and.dir.eq.1).or.(comp.eq.'V'.and.dir.eq.2).or.(comp.eq.'W'.and.dir.eq.3)
         ! Compute fill edge and source index
         ! Normal: fill includes boundary face, source from first interior
         ! Tangent: fill ghosts only, source from boundary
         if (is_normal) then
            fill_edge=bnd; src_from=bnd+merge(1,-1,is_lo)
         else
            fill_edge=bnd+merge(-1,1,is_lo); src_from=bnd
         end if
         ! Compute slab bounds for fills (ext_dir, foextrap)
         slo=[ilo,jlo,klo]; shi=[ihi,jhi,khi]
         if (is_lo) then; shi(dir)=fill_edge; else; slo(dir)=fill_edge; end if
         ! Tangent mirror offset: -1 for lo, +1 for hi
         toff=merge(-1,1,is_lo)
         select case (bctype)
          case (amrex_bc_ext_dir)
            ! User-controlled: fill via user callback
            if (associated(solver%user_mpinc_bc)) then
               bc_bx=amrex_box(slo,shi)
               call solver%user_mpinc_bc(lvl=lvl,time=time,face=face,bx=bc_bx,comp=comp,p=p)
            else
               ! Default to zero if no callback provided
               do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1)
                  p(i,j,k,1)=0.0_WP
               end do; end do; end do
            end if
          case (amrex_bc_foextrap)
            ! Neumann: zero-gradient extrapolation
            if (is_normal) then
               ! Normal component: boundary face is a solver DOF (set by correct_outflow/projection)
               ! Only fill ghosts beyond it, copying from the boundary face value
               if (is_lo) then; shi(dir)=bnd-1; else; slo(dir)=bnd+1; end if
               select case (dir)
                case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(bnd,j,k,1); end do; end do; end do
                case (2); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,bnd,k,1); end do; end do; end do
                case (3); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,j,bnd,1); end do; end do; end do
               end select
            else
                ! Tangent component: fill ghosts from boundary
                select case (dir)
                 case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(src_from,j,k,1); end do; end do; end do
                 case (2); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,src_from,k,1); end do; end do; end do
                 case (3); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,j,src_from,1); end do; end do; end do
                end select
            end if
          case (amrex_bc_reflect_even)
            ! Symmetry: F(-n) = F(n) - ghosts only
            ! Adjust slab to ghosts only
            if (is_lo) then; shi(dir)=bnd-1; else; slo(dir)=bnd+1; end if
            if (is_normal) then
               ! Normal: mirror across wall face (source=2*bnd-ii)
               select case (dir)
                case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(2*bnd-i,j,k,1); end do; end do; end do
                case (2); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,2*bnd-j,k,1); end do; end do; end do
                case (3); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,j,2*bnd-k,1); end do; end do; end do
               end select
            else
               ! Tangent: mirror around half-cell (source=2*bnd-ii+toff)
               select case (dir)
                case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(2*bnd-i+toff,j,k,1); end do; end do; end do
                case (2); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,2*bnd-j+toff,k,1); end do; end do; end do
                case (3); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=p(i,j,2*bnd-k+toff,1); end do; end do; end do
               end select
            end if
          case (amrex_bc_reflect_odd)
            ! Anti-symmetry: F(-n) = -F(n)
            ! Adjust slab to ghosts only
            if (is_lo) then; shi(dir)=bnd-1; else; slo(dir)=bnd+1; end if
            if (is_normal) then
               ! Zero the wall face
               select case (dir)
                case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); p(bnd,j,k,1)=0.0_WP; end do; end do
                case (2); do k=slo(3),shi(3); do i=slo(1),shi(1); p(i,bnd,k,1)=0.0_WP; end do; end do
                case (3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,bnd,1)=0.0_WP; end do; end do
               end select
               ! Mirror across wall face (source=2*bnd-ii)
               select case (dir)
                case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=-p(2*bnd-i,j,k,1); end do; end do; end do
                case (2); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=-p(i,2*bnd-j,k,1); end do; end do; end do
                case (3); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=-p(i,j,2*bnd-k,1); end do; end do; end do
               end select
            else
               ! Tangent: mirror around half-cell (source=2*bnd-ii+toff)
               select case (dir)
                case (1); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=-p(2*bnd-i+toff,j,k,1); end do; end do; end do
                case (2); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=-p(i,2*bnd-j+toff,k,1); end do; end do; end do
                case (3); do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1); p(i,j,k,1)=-p(i,j,2*bnd-k+toff,1); end do; end do; end do
               end select
            end if
         end select
      end subroutine apply_vel_bc

   end subroutine velocity_fillbc

   ! ============================================================================
   ! Collocated velocity boundary conditions
   ! ============================================================================

   !> Internal fillbc for UVW - calls default_fillbc first, then user_bc for ext_dir faces
   subroutine UVW_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,&
      &                           amrex_geometry,amrex_multifab,amrex_bc_ext_dir
      use amrdata_class, only: default_fillbc
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp,ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bc_bx
      class(amrmpinc), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: ilo,ihi,jlo,jhi,klo,khi,n
      integer, dimension(3) :: dlo,dhi
      integer :: lvl
      character(len=1), parameter :: cname(3)=['U','V','W']

      ! First apply default BC handling (foextrap,hoextrap,reflect,etc.)
      call default_fillbc(this,mf,scomp,ncomp,time,geom)

      ! Access parent solver
      select type (s=>this%parent)
       class is (amrmpinc)
         solver=>s
      end select

      ! Check if user callback exists
      if (.not.associated(solver%user_mpinc_bc)) return

      ! Get domain bounds
      dlo=geom%domain%lo
      dhi=geom%domain%hi

      ! Get current level
      lvl=this%fill_lvl_cache

      ! Loop over FABs and apply user_bc for ext_dir faces
      call amrex_mfiter_build(mfi,mf,tiling=.false.)
      do while (mfi%next())
         p=>mf%dataptr(mfi)
         ilo=lbound(p,1); ihi=ubound(p,1)
         jlo=lbound(p,2); jhi=ubound(p,2)
         klo=lbound(p,3); khi=ubound(p,3)
         ! X-LOW (face=1)
         if (ilo.lt.dlo(1)) then
            bc_bx=amrex_box([ilo,jlo,klo],[dlo(1)-1,jhi,khi])
            do n=1,3; if (this%lo_bc(1,n).eq.amrex_bc_ext_dir) call solver%user_mpinc_bc(lvl=lvl,time=time,face=1,bx=bc_bx,comp=cname(n),p=p); end do
         end if
         ! X-HIGH (face=2)
         if (ihi.gt.dhi(1)) then
            bc_bx=amrex_box([dhi(1)+1,jlo,klo],[ihi,jhi,khi])
            do n=1,3; if (this%hi_bc(1,n).eq.amrex_bc_ext_dir) call solver%user_mpinc_bc(lvl=lvl,time=time,face=2,bx=bc_bx,comp=cname(n),p=p); end do
         end if
         ! Y-LOW (face=3)
         if (jlo.lt.dlo(2)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,dlo(2)-1,khi])
            do n=1,3; if (this%lo_bc(2,n).eq.amrex_bc_ext_dir) call solver%user_mpinc_bc(lvl=lvl,time=time,face=3,bx=bc_bx,comp=cname(n),p=p); end do
         end if
         ! Y-HIGH (face=4)
         if (jhi.gt.dhi(2)) then
            bc_bx=amrex_box([ilo,dhi(2)+1,klo],[ihi,jhi,khi])
            do n=1,3; if (this%hi_bc(2,n).eq.amrex_bc_ext_dir) call solver%user_mpinc_bc(lvl=lvl,time=time,face=4,bx=bc_bx,comp=cname(n),p=p); end do
         end if
         ! Z-LOW (face=5)
         if (klo.lt.dlo(3)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,jhi,dlo(3)-1])
            do n=1,3; if (this%lo_bc(3,n).eq.amrex_bc_ext_dir) call solver%user_mpinc_bc(lvl=lvl,time=time,face=5,bx=bc_bx,comp=cname(n),p=p); end do
         end if
         ! Z-HIGH (face=6)
         if (khi.gt.dhi(3)) then
            bc_bx=amrex_box([ilo,jlo,dhi(3)+1],[ihi,jhi,khi])
            do n=1,3; if (this%hi_bc(3,n).eq.amrex_bc_ext_dir) call solver%user_mpinc_bc(lvl=lvl,time=time,face=6,bx=bc_bx,comp=cname(n),p=p); end do
         end if
      end do
      call amrex_mfiter_destroy(mfi)

   end subroutine UVW_fillbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Copy current state to old state
   subroutine store_old(this)
      implicit none
      class(amrmpinc), intent(inout) :: this
      ! Store amrvof state
      call this%amrvof%store_old()
      ! Store velocity state
      call this%UVWold%copy(src=this%UVW)
      call this%Uold%copy(src=this%U)
      call this%Vold%copy(src=this%V)
      call this%Wold%copy(src=this%W)
   end subroutine store_old

   !> Compute divergence of velocity into internal div field, update divmax
   !> Uses composite fine masking so covered coarse cells don't pollute divmax
   subroutine get_div(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmfab_compute_divergence,amrmask_make_fine
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pDiv
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      ! Use our wrapper to amrex's 2nd order staggered divergence
      do lvl=0,this%amr%clvl()
         call amrmfab_compute_divergence(this%div%mf(lvl),this%U%mf(lvl),this%V%mf(lvl),this%W%mf(lvl),this%amr%geom(lvl))
      end do
      ! Update divmax using composite fine masking
      this%divmax=0.0_WP
      do lvl=0,this%amr%clvl()
         ! Build fine mask for this level (if not finest)
         if (lvl.lt.this%amr%clvl()) then
            call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
            call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
         end if
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            pDiv=>this%div%mf(lvl)%dataptr(mfi)
            if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
               this%divmax=max(this%divmax,abs(pDiv(i,j,k,1)))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%divmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
   end subroutine get_div

   !> Interpolate cell-centered UVW to face U,V,W using density-weighting
   subroutine interp_vel_to_face(this)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW,pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pSubVF
      real(WP) :: rho_lo,rho_hi
      ! Traverse levels
      do lvl=0,this%amr%clvl()
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            ! Get X-face velocity
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               rho_lo=this%rhoL*pVF(i-1,j,k,1)+this%rhoG*(1.0_WP-pVF(i-1,j,k,1))
               rho_hi=this%rhoL*pVF(i  ,j,k,1)+this%rhoG*(1.0_WP-pVF(i  ,j,k,1))
               if (lvl.eq.this%amr%maxlvl) then
                  rho_lo=this%rhoL*pSubVF(i-1,j,k,2)+this%rhoG*(1.0_WP-pSubVF(i-1,j,k,2))
                  rho_hi=this%rhoL*pSubVF(i  ,j,k,1)+this%rhoG*(1.0_WP-pSubVF(i  ,j,k,1))
               end if
               pU(i,j,k,1)=(rho_lo*pUVW(i-1,j,k,1)+rho_hi*pUVW(i,j,k,1))/(rho_lo+rho_hi)
            end do; end do; end do
            ! Get Y-face velocity
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               rho_lo=this%rhoL*pVF(i,j-1,k,1)+this%rhoG*(1.0_WP-pVF(i,j-1,k,1))
               rho_hi=this%rhoL*pVF(i,j  ,k,1)+this%rhoG*(1.0_WP-pVF(i,j  ,k,1))
               if (lvl.eq.this%amr%maxlvl) then
                  rho_lo=this%rhoL*pSubVF(i,j-1,k,4)+this%rhoG*(1.0_WP-pSubVF(i,j-1,k,4))
                  rho_hi=this%rhoL*pSubVF(i,j  ,k,3)+this%rhoG*(1.0_WP-pSubVF(i,j  ,k,3))
               end if
               pV(i,j,k,1)=(rho_lo*pUVW(i,j-1,k,2)+rho_hi*pUVW(i,j,k,2))/(rho_lo+rho_hi)
            end do; end do; end do
            ! Get Z-face velocity
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               rho_lo=this%rhoL*pVF(i,j,k-1,1)+this%rhoG*(1.0_WP-pVF(i,j,k-1,1))
               rho_hi=this%rhoL*pVF(i,j,k  ,1)+this%rhoG*(1.0_WP-pVF(i,j,k  ,1))
               if (lvl.eq.this%amr%maxlvl) then
                  rho_lo=this%rhoL*pSubVF(i,j,k-1,6)+this%rhoG*(1.0_WP-pSubVF(i,j,k-1,6))
                  rho_hi=this%rhoL*pSubVF(i,j,k  ,5)+this%rhoG*(1.0_WP-pSubVF(i,j,k  ,5))
               end if
               pW(i,j,k,1)=(rho_lo*pUVW(i,j,k-1,3)+rho_hi*pUVW(i,j,k,3))/(rho_lo+rho_hi)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine interp_vel_to_face

   !> Prepare variable-coefficient pressure solver using face densities
   subroutine prepare_psolver(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_multifab
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP) :: VF_f
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      type(amrex_multifab), dimension(:), allocatable :: Bx,By,Bz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pSubVF,pBx,pBy,pBz
      ! Allocate temporary face coefficient mfabs
      allocate(Bx(0:this%amr%clvl()),By(0:this%amr%clvl()),Bz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Bx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,By(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Bz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Fill 1/rho_face at all levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pBx=>Bx(lvl)%dataptr(mfi)
            pBy=>By(lvl)%dataptr(mfi)
            pBz=>Bz(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            ! X-faces
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               VF_f=0.5_WP*(pVF(i-1,j,k,1)+pVF(i,j,k,1)); if (lvl.eq.this%amr%maxlvl) VF_f=0.5_WP*(pSubVF(i-1,j,k,2)+pSubVF(i,j,k,1))
               pBx(i,j,k,1)=1.0_WP/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
            end do; end do; end do
            ! Y-faces
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               VF_f=0.5_WP*(pVF(i,j-1,k,1)+pVF(i,j,k,1)); if (lvl.eq.this%amr%maxlvl) VF_f=0.5_WP*(pSubVF(i,j-1,k,4)+pSubVF(i,j,k,3))
               pBy(i,j,k,1)=1.0_WP/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
            end do; end do; end do
            ! Z-faces
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               VF_f=0.5_WP*(pVF(i,j,k-1,1)+pVF(i,j,k,1)); if (lvl.eq.this%amr%maxlvl) VF_f=0.5_WP*(pSubVF(i,j,k-1,6)+pSubVF(i,j,k,5))
               pBz(i,j,k,1)=1.0_WP/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Rebuild operator
      call this%psolver%setup(bcoef_x=Bx,bcoef_y=By,bcoef_z=Bz)
      ! Destroy temporary mfabs
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Bx(lvl))
         call this%amr%mfab_destroy(By(lvl))
         call this%amr%mfab_destroy(Bz(lvl))
      end do
   end subroutine prepare_psolver

   !> Add (-scale*pressure gradient/rho) to both U/V/W and UVW velocities. Two flavors:
   !>   phi present -> direct path: use explicit stencil that reads phi ghost cells directly (for predictor with fs%P)
   !>   phi absent  -> MLMG path:   use psolver internal fluxes (for projection with dP)
   !> Cell-center correction averages the face gradients back to cell center
   subroutine correct_both_velocities(this,scale,phi)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrdata), intent(in), optional :: phi
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP,pFx,pFy,pFz,pVF,pSubVF
      real(WP) :: dxi,dyi,dzi,VF_f
      integer :: lvl,i,j,k
      ! Build temp face mfabs to store pressure fluxes
      allocate(Fx(0:this%amr%clvl()),Fy(0:this%amr%clvl()),Fz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Compute -pressure gradient/rho at faces
      if (present(phi)) then
         ! Use provided phi and its ghosts cells
         do lvl=0,this%amr%clvl()
            dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pP =>phi%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  VF_f=0.5_WP*(pVF(i-1,j,k,1)+pVF(i,j,k,1)); if (lvl.eq.this%amr%maxlvl) VF_f=0.5_WP*(pSubVF(i-1,j,k,2)+pSubVF(i,j,k,1))
                  pFx(i,j,k,1)=-(pP(i,j,k,1)-pP(i-1,j,k,1))*dxi/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
               end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  VF_f=0.5_WP*(pVF(i,j-1,k,1)+pVF(i,j,k,1)); if (lvl.eq.this%amr%maxlvl) VF_f=0.5_WP*(pSubVF(i,j-1,k,4)+pSubVF(i,j,k,3))
                  pFy(i,j,k,1)=-(pP(i,j,k,1)-pP(i,j-1,k,1))*dyi/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
               end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  VF_f=0.5_WP*(pVF(i,j,k-1,1)+pVF(i,j,k,1)); if (lvl.eq.this%amr%maxlvl) VF_f=0.5_WP*(pSubVF(i,j,k-1,6)+pSubVF(i,j,k,5))
                  pFz(i,j,k,1)=-(pP(i,j,k,1)-pP(i,j,k-1,1))*dzi/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      else
         ! Use psolver's solution and its internal ghosts
         call this%psolver%get_fluxes(Fx,Fy,Fz)
      end if
      ! Delegate face+cell application to apply_face_fluxes method
      call this%apply_face_fluxes(scale,Fx,Fy,Fz)
      ! Destroy temps
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Fx(lvl))
         call this%amr%mfab_destroy(Fy(lvl))
         call this%amr%mfab_destroy(Fz(lvl))
      end do
   end subroutine correct_both_velocities

   !> Add CSF surface-tension velocity increment
   !> Builds face fluxes 1/rho*sigma*kappa*grad(VF) at maxlvl, avg_down, then calls apply_face_fluxes
   subroutine add_surface_tension(this,scale)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      use amrex_interface, only: amrmfab_average_down_face
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrex_multifab), dimension(:), allocatable :: STFx,STFy,STFz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSTFx,pSTFy,pSTFz,pVF,pSubVF,pCurv,pSD
      real(WP) :: dxi,dyi,dzi,VF_f,mysurf,mycurv
      integer :: lvl,i,j,k
      ! Guard: no surface tension or clvl<maxlvl
      if (this%sigma.eq.0.0_WP.or.this%amr%clvl().lt.this%amr%maxlvl) return
      ! Build temp face flux mfabs
      allocate(STFx(0:this%amr%clvl()),STFy(0:this%amr%clvl()),STFz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,STFx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.]); call STFx(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,STFy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.]); call STFy(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,STFz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ]); call STFz(lvl)%setval(0.0_WP)
      end do
      ! Compute ST face fluxes at finest level
      lvl=this%amr%maxlvl
      dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF   =>this%VF%mf(lvl)%dataptr(mfi)
         pSubVF=>this%subVF%dataptr(mfi)
         pCurv =>this%curv%dataptr(mfi)
         pSD   =>this%SD%dataptr(mfi)
         pSTFx =>STFx(lvl)%dataptr(mfi)
         pSTFy =>STFy(lvl)%dataptr(mfi)
         pSTFz =>STFz(lvl)%dataptr(mfi)
         ! X-faces
         bx=mfi%nodaltilebox(1)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            mycurv=0.0_WP
            mysurf=sum(pSD(i-1:i,j,k,1)); if (mysurf.gt.0.0_WP) mycurv=sum(pSD(i-1:i,j,k,1)*pCurv(i-1:i,j,k,1))/mysurf
            VF_f=0.5_WP*(pSubVF(i-1,j,k,2)+pSubVF(i,j,k,1))
            pSTFx(i,j,k,1)=this%sigma*mycurv*(pVF(i,j,k,1)-pVF(i-1,j,k,1))*dxi/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
         end do; end do; end do
         ! Y-faces
         bx=mfi%nodaltilebox(2)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            mycurv=0.0_WP
            mysurf=sum(pSD(i,j-1:j,k,1)); if (mysurf.gt.0.0_WP) mycurv=sum(pSD(i,j-1:j,k,1)*pCurv(i,j-1:j,k,1))/mysurf
            VF_f=0.5_WP*(pSubVF(i,j-1,k,4)+pSubVF(i,j,k,3))
            pSTFy(i,j,k,1)=this%sigma*mycurv*(pVF(i,j,k,1)-pVF(i,j-1,k,1))*dyi/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
         end do; end do; end do
         ! Z-faces
         bx=mfi%nodaltilebox(3)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            mycurv=0.0_WP
            mysurf=sum(pSD(i,j,k-1:k,1)); if (mysurf.gt.0.0_WP) mycurv=sum(pSD(i,j,k-1:k,1)*pCurv(i,j,k-1:k,1))/mysurf
            VF_f=0.5_WP*(pSubVF(i,j,k-1,6)+pSubVF(i,j,k,5))
            pSTFz(i,j,k,1)=this%sigma*mycurv*(pVF(i,j,k,1)-pVF(i,j,k-1,1))*dzi/(this%rhoL*VF_f+this%rhoG*(1.0_WP-VF_f))
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Average down face fluxes from finest to coarser levels
      do lvl=this%amr%clvl()-1,0,-1
         call amrmfab_average_down_face(fmf=STFx(lvl+1),cmf=STFx(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_face(fmf=STFy(lvl+1),cmf=STFy(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
         call amrmfab_average_down_face(fmf=STFz(lvl+1),cmf=STFz(lvl),rr=[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl))
      end do
      ! Apply fluxes
      call this%apply_face_fluxes(scale,STFx,STFy,STFz)
      ! Destroy temps
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(STFx(lvl))
         call this%amr%mfab_destroy(STFy(lvl))
         call this%amr%mfab_destroy(STFz(lvl))
      end do
      deallocate(STFx,STFy,STFz)
   end subroutine add_surface_tension

   !> Apply pre-built face fluxes to both face and cell-centered velocities
   !> scale * Fx/Fy/Fz is saxpy'd onto U/V/W, then density-weighted to UVW
   subroutine apply_face_fluxes(this,scale,Fx,Fy,Fz)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrex_multifab), intent(in) :: Fx(0:),Fy(0:),Fz(0:)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pUVW,pVF,pSubVF
      real(WP) :: rho,rhoLo,rhoHi
      integer :: lvl,i,j,k
      do lvl=0,this%amr%clvl()
         ! Face velocities
         call this%U%mf(lvl)%saxpy(scale,Fx(lvl),1,1,1,0)
         call this%V%mf(lvl)%saxpy(scale,Fy(lvl),1,1,1,0)
         call this%W%mf(lvl)%saxpy(scale,Fz(lvl),1,1,1,0)
         ! Cell-centred velocities: density-weighted average of face fluxes
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
               ! X
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,1)+this%rhoG*(1.0_WP-pSubVF(i,j,k,1))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,2)+this%rhoG*(1.0_WP-pSubVF(i,j,k,2))
               pUVW(i,j,k,1)=pUVW(i,j,k,1)+scale*0.5_WP*(rhoLo*pFx(i,j,k,1)+rhoHi*pFx(i+1,j,k,1))/rho
               ! Y
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,3)+this%rhoG*(1.0_WP-pSubVF(i,j,k,3))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,4)+this%rhoG*(1.0_WP-pSubVF(i,j,k,4))
               pUVW(i,j,k,2)=pUVW(i,j,k,2)+scale*0.5_WP*(rhoLo*pFy(i,j,k,1)+rhoHi*pFy(i,j+1,k,1))/rho
               ! Z
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,5)+this%rhoG*(1.0_WP-pSubVF(i,j,k,5))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,6)+this%rhoG*(1.0_WP-pSubVF(i,j,k,6))
               pUVW(i,j,k,3)=pUVW(i,j,k,3)+scale*0.5_WP*(rhoLo*pFz(i,j,k,1)+rhoHi*pFz(i,j,k+1,1))/rho
            end do; end do; end do
            ! Non-periodic boundary cells (only one face flux available)
            if (.not.this%amr%xper.and.bx%lo(1).eq.this%amr%geom(lvl)%domain%lo(1)) then
               i=this%amr%geom(lvl)%domain%lo(1); do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,2)+this%rhoG*(1.0_WP-pSubVF(i,j,k,2))
                  pUVW(i,j,k,1)=pUVW(i,j,k,1)+scale*0.5_WP*rhoHi/rho*pFx(i+1,j,k,1)
               end do; end do
            end if
            if (.not.this%amr%xper.and.bx%hi(1).eq.this%amr%geom(lvl)%domain%hi(1)) then
               i=this%amr%geom(lvl)%domain%hi(1); do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,1)+this%rhoG*(1.0_WP-pSubVF(i,j,k,1))
                  pUVW(i,j,k,1)=pUVW(i,j,k,1)+scale*0.5_WP*rhoLo/rho*pFx(i,  j,k,1)
               end do; end do
            end if
            if (.not.this%amr%yper.and.bx%lo(2).eq.this%amr%geom(lvl)%domain%lo(2)) then
               j=this%amr%geom(lvl)%domain%lo(2); do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,4)+this%rhoG*(1.0_WP-pSubVF(i,j,k,4))
                  pUVW(i,j,k,2)=pUVW(i,j,k,2)+scale*0.5_WP*rhoHi/rho*pFy(i,j+1,k,1)
               end do; end do
            end if
            if (.not.this%amr%yper.and.bx%hi(2).eq.this%amr%geom(lvl)%domain%hi(2)) then
               j=this%amr%geom(lvl)%domain%hi(2); do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,3)+this%rhoG*(1.0_WP-pSubVF(i,j,k,3))
                  pUVW(i,j,k,2)=pUVW(i,j,k,2)+scale*0.5_WP*rhoLo/rho*pFy(i,j,  k,1)
               end do; end do
            end if
            if (.not.this%amr%zper.and.bx%lo(3).eq.this%amr%geom(lvl)%domain%lo(3)) then
               k=this%amr%geom(lvl)%domain%lo(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,6)+this%rhoG*(1.0_WP-pSubVF(i,j,k,6))
                  pUVW(i,j,k,3)=pUVW(i,j,k,3)+scale*0.5_WP*rhoHi/rho*pFz(i,j,k+1,1)
               end do; end do
            end if
            if (.not.this%amr%zper.and.bx%hi(3).eq.this%amr%geom(lvl)%domain%hi(3)) then
               k=this%amr%geom(lvl)%domain%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,5)+this%rhoG*(1.0_WP-pSubVF(i,j,k,5))
                  pUVW(i,j,k,3)=pUVW(i,j,k,3)+scale*0.5_WP*rhoLo/rho*pFz(i,j,k,  1)
               end do; end do
            end if
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine apply_face_fluxes

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Compute momentum advection and viscous terms for all levels, returning d(UVW)/dt
   !> No pressure gradient, user can add it in the main loop
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dUVWdt(this,dUVWdt,dt,time)
      implicit none
      class(amrmpinc), intent(inout) :: this
      class(amrdata), intent(inout) :: dUVWdt                        ! Output: velocity RHS (cell-centered)
      real(WP), intent(in) :: dt,time

      ! Get dmomdt
      call this%get_dmomdt(dmomdt=dUVWdt,dt=dt,time=time)

      ! Transform dmomdt to dUVWdt
      transform_to_velocity: block
         use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: lvl,i,j,k
         real(WP) :: rho,rho_old
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pRhs,pVF,pUVWold,pVFold
         do lvl=0,this%amr%clvl()
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               pRhs   =>dUVWdt%mf(lvl)%dataptr(mfi)
               pUVWold=>this%UVWold%mf(lvl)%dataptr(mfi)
               pVFold =>this%VFold%mf(lvl)%dataptr(mfi)
               pVF    =>this%VF%mf(lvl)%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho    =this%rhoL*pVF   (i,j,k,1)+this%rhoG*(1.0_WP-pVF   (i,j,k,1))
                  rho_old=this%rhoL*pVFold(i,j,k,1)+this%rhoG*(1.0_WP-pVFold(i,j,k,1))
                  pRhs(i,j,k,:)=pRhs(i,j,k,:)/rho+(rho_old-rho)/(rho*dt)*pUVWold(i,j,k,:)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block transform_to_velocity

   end subroutine get_dUVWdt

   !> Compute momentum advection and viscous terms for all levels, returning d(mom)/dt
   !> No pressure gradient, user can add it in the main loop
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dmomdt(this,dmomdt,dt,time)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      implicit none
      class(amrmpinc), intent(inout) :: this
      class(amrdata), intent(inout) :: dmomdt                        ! Output: momentum RHS (cell-centered)
      real(WP), intent(in) :: dt,time
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz ! Flux mfabs
      type(amrex_multifab) :: Vx,Vy,Vz
      type(amrex_multifab) :: band
      ! Shared variables for internal functions
      real(WP) :: dx,dy,dz,dxi,dyi,dzi                               ! Needed for SL transport
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW  ! Velocity used for project
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLICold  ! PLICold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVWold   ! UVWold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold    ! VFold used in tet2flux_plic
      logical :: crossed_plic ! Used in tet2flux/tet2flux_plic

      ! Build transport band at finest level to localize SL computation
      call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=band,ncomp=1,nover=1)
      call this%build_band(lvl=this%amr%clvl(),VF=this%VFold%mf(this%amr%clvl()),band=band,nband=2)

      ! Allocate all fluxes
      define_fluxes: block
         integer :: lvl
         ! Face-centered conserved variable fluxes (3 components: rhoU, rhoV, rhoW)
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_build(lvl=lvl,mfab=Fx(lvl),ncomp=3,nover=0,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl=lvl,mfab=Fy(lvl),ncomp=3,nover=0,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl=lvl,mfab=Fz(lvl),ncomp=3,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
         end do
         ! Volume moment fluxes at finest level (8 components: Lvol,Gvol,Lbar,Gbar)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vx,ncomp=8,nover=0,atface=[.true. ,.false.,.false.]); call Vx%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vy,ncomp=8,nover=0,atface=[.false.,.true. ,.false.]); call Vy%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vz,ncomp=8,nover=0,atface=[.false.,.false.,.true. ]); call Vz%setval(0.0_WP)
      end block define_fluxes

      ! Phase 1a: Semi-Lagrangian fluxes at finest level
      semilagrangian_fluxes: block
         use amrvof_geometry, only: tet_sign,tet_map,correct_flux_poly
         integer :: lvl,i,j,k,n,nn
         real(WP), dimension(3,9) :: face
         real(WP), dimension(3,4) :: tet
         integer , dimension(3,4) :: ijk
         integer , dimension(3,9) :: fijk
         real(WP), dimension(:,:,:,:), allocatable :: proj
         real(WP), dimension(8) :: Vflux
         real(WP), dimension(3) :: Qflux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pFx,pFy,pFz,pUVW
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx,nbx
         ! Skip if clvl < maxlvl
         if (this%amr%clvl().lt.this%amr%maxlvl) exit semilagrangian_fluxes
         ! Get finest level info
         lvl=this%amr%maxlvl
         dx=this%amr%dx(lvl); dxi=1.0_WP/this%amr%dx(lvl)
         dy=this%amr%dy(lvl); dyi=1.0_WP/this%amr%dy(lvl)
         dz=this%amr%dz(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         ! Loop over finest level tiles
         call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
         do while (mfi%next())
            ! Get data pointers: PLICold, Qold, VFold, band, velocity, fluxes
            pPLICold=>this%PLICold%dataptr(mfi)
            pUVWold =>this%UVWold%mf(lvl)%dataptr(mfi)
            pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
            pBand   =>band%dataptr(mfi)
            pUVW    =>this%UVW%mf(lvl)%dataptr(mfi)
            pU      =>this%U%mf(lvl)%dataptr(mfi)
            pV      =>this%V%mf(lvl)%dataptr(mfi)
            pW      =>this%W%mf(lvl)%dataptr(mfi)
            pVx     =>Vx%dataptr(mfi)
            pVy     =>Vy%dataptr(mfi)
            pVz     =>Vz%dataptr(mfi)
            pFx     =>Fx(lvl)%dataptr(mfi)
            pFy     =>Fy(lvl)%dataptr(mfi)
            pFz     =>Fz(lvl)%dataptr(mfi)
            ! Remap vertices in the band via RK2
            nbx=mfi%nodaltilebox()
            allocate(proj(3,nbx%lo(1):nbx%hi(1),nbx%lo(2):nbx%hi(2),nbx%lo(3):nbx%hi(3)))
            do k=nbx%lo(3),nbx%hi(3); do j=nbx%lo(2),nbx%hi(2); do i=nbx%lo(1),nbx%hi(1)
               if (maxval(pBand(i-1:i,j-1:j,k-1:k,1)).gt.0.0_WP) proj(:,i,j,k)=project([this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k,WP)*dz],-dt)
            end do; end do; end do
            ! X-fluxes
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i-1:i,j,k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,5)=proj(:,i,j  ,k  )
               face(:,2)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=proj(:,i,j  ,k+1)
               face(:,3)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,7)=proj(:,i,j+1,k+1)
               face(:,4)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=proj(:,i,j+1,k  )
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dy*dz*pU(i,j,k,1))
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(1,nn)=merge(i-1,i,pU(i,j,k,1).gt.0.0_WP); end do
               ! Are we crossing plic?
               crossed_plic=.false.
               ! Decompose into tets, cut, and accumulate
               pVx(i,j,k,:)=0.0_WP
               pFx(i,j,k,:)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVx(i,j,k,:)=pVx(i,j,k,:)+tet_sign(tet)*Vflux
                  pFx(i,j,k,:)=pFx(i,j,k,:)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFx(i,j,k,:)=-pFx(i,j,k,:)/(dt*dy*dz)
               ! Switch to dissipation-free momentum flux for BB-pure regions
               if (.not.crossed_plic) pFx(i,j,k,:)=-(this%rhoL*pVx(i,j,k,1)+this%rhoG*pVx(i,j,k,2))/(dt*dy*dz)*0.5_WP*(pUVW(i-1,j,k,:)+pUVW(i,j,k,:))
            end do; end do; end do
            ! Y-fluxes
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,5)=proj(:,i+1,j,k+1)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=proj(:,i  ,j,k+1)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,7)=proj(:,i  ,j,k  )
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=proj(:,i+1,j,k  )
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dz*dx*pV(i,j,k,1))
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(2,nn)=merge(j-1,j,pV(i,j,k,1).gt.0.0_WP); end do
               ! Are we crossing plic?
               crossed_plic=.false.
               ! Decompose into tets, cut, and accumulate
               pVy(i,j,k,:)=0.0_WP
               pFy(i,j,k,:)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVy(i,j,k,:)=pVy(i,j,k,:)+tet_sign(tet)*Vflux
                  pFy(i,j,k,:)=pFy(i,j,k,:)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFy(i,j,k,:)=-pFy(i,j,k,:)/(dt*dz*dx)
               ! Switch to dissipation-free momentum flux for BB-pure regions
               if (.not.crossed_plic) pFy(i,j,k,:)=-(this%rhoL*pVy(i,j,k,1)+this%rhoG*pVy(i,j,k,2))/(dt*dz*dx)*0.5_WP*(pUVW(i,j-1,k,:)+pUVW(i,j,k,:))
            end do; end do; end do
            ! Z-fluxes
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,5)=proj(:,i+1,j  ,k)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,6)=proj(:,i  ,j  ,k)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,7)=proj(:,i  ,j+1,k)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,8)=proj(:,i+1,j+1,k)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dx*dy*pW(i,j,k,1))
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(3,nn)=merge(k-1,k,pW(i,j,k,1).gt.0.0_WP); end do
               ! Are we crossing plic?
               crossed_plic=.false.
               ! Decompose into tets, cut, and accumulate
               pVz(i,j,k,:)=0.0_WP
               pFz(i,j,k,:)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVz(i,j,k,:)=pVz(i,j,k,:)+tet_sign(tet)*Vflux
                  pFz(i,j,k,:)=pFz(i,j,k,:)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFz(i,j,k,:)=-pFz(i,j,k,:)/(dt*dx*dy)
               ! Switch to dissipation-free momentum flux for BB-pure regions
               if (.not.crossed_plic) pFz(i,j,k,:)=-(this%rhoL*pVz(i,j,k,1)+this%rhoG*pVz(i,j,k,2))/(dt*dx*dy)*0.5_WP*(pUVW(i,j,k-1,:)+pUVW(i,j,k,:))
            end do; end do; end do
            ! Deallocate proj for this tile
            deallocate(proj)
         end do
         call this%amr%mfiter_destroy(mfi)
      end block semilagrangian_fluxes

      ! Phase 1b: Finite volume fluxes for all levels (Euler fluxes skip band cells at finest level)
      finitevolume_fluxes: block
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW,pFx,pFy,pFz,pBand,pVF,pVisc
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Intentional masking
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: div,visc_f,mass_flux
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         integer :: lvl,i,j,k
         logical :: in_band
         do lvl=0,this%amr%clvl()
            ! Grid spacings for this level
            dx=this%amr%dx(lvl); dxi=1.0_WP/this%amr%dx(lvl)
            dy=this%amr%dy(lvl); dyi=1.0_WP/this%amr%dy(lvl)
            dz=this%amr%dz(lvl); dzi=1.0_WP/this%amr%dz(lvl)
            ! Loop over tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               pUVW =>this%UVW%mf(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pFx  =>Fx(lvl)%dataptr(mfi)
               pFy  =>Fy(lvl)%dataptr(mfi)
               pFz  =>Fz(lvl)%dataptr(mfi)
               if (lvl.eq.this%amr%clvl()) pBand=>band%dataptr(mfi)
               ! X-fluxes
               fbx=mfi%nodaltilebox(1)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Check if in band
                  if (lvl.eq.this%amr%clvl()) then; in_band=maxval(pBand(i-1:i,j,k,1)).gt.0.0_WP; else; in_band=.false.; end if
                  ! Outside band, compute finite volume Euler fluxes
                  if (.not.in_band) then
                     mass_flux=-merge(this%rhoL,this%rhoG,pVF(i,j,k,1).gt.0.5_WP)*pU(i,j,k,1)
                     pFx(i,j,k,:)=mass_flux*0.5_WP*(pUVW(i-1,j,k,:)+pUVW(i,j,k,:))
                  end if
                  ! Velocity gradients at x-face
                  gradU(1,1)=dxi*(pUVW(i,j,k,1)-pUVW(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*sum(pUVW(i-1:i,j+1,k,1)-pUVW(i-1:i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*sum(pUVW(i-1:i,j,k+1,1)-pUVW(i-1:i,j,k-1,1))
                  gradU(1,2)=dxi*(pUVW(i,j,k,2)-pUVW(i-1,j,k,2))
                  gradU(2,2)=0.25_WP*dyi*sum(pUVW(i-1:i,j+1,k,2)-pUVW(i-1:i,j-1,k,2))
                  gradU(3,2)=0.25_WP*dzi*sum(pUVW(i-1:i,j,k+1,2)-pUVW(i-1:i,j,k-1,2))
                  gradU(1,3)=dxi*(pUVW(i,j,k,3)-pUVW(i-1,j,k,3))
                  gradU(2,3)=0.25_WP*dyi*sum(pUVW(i-1:i,j+1,k,3)-pUVW(i-1:i,j-1,k,3))
                  gradU(3,3)=0.25_WP*dzi*sum(pUVW(i-1:i,j,k+1,3)-pUVW(i-1:i,j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscosity at x-face
                  visc_f=2.0_WP*product(pVisc(i-1:i,j,k,1))/(sum(pVisc(i-1:i,j,k,1))+tiny(1.0_WP))
                  ! Viscous stress at x-face
                  pFx(i,j,k,1)=pFx(i,j,k,1)+visc_f*(gradU(1,1)+gradU(1,1))-2.0_WP/3.0_WP*visc_f*div
                  pFx(i,j,k,2)=pFx(i,j,k,2)+visc_f*(gradU(2,1)+gradU(1,2))
                  pFx(i,j,k,3)=pFx(i,j,k,3)+visc_f*(gradU(3,1)+gradU(1,3))
               end do; end do; end do
               ! Y-fluxes
               fbx=mfi%nodaltilebox(2)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Check if in band
                  if (lvl.eq.this%amr%clvl()) then; in_band=maxval(pBand(i,j-1:j,k,1)).gt.0.0_WP; else; in_band=.false.; end if
                  ! Outside band, compute finite volume Euler fluxes
                  if (.not.in_band) then
                     mass_flux=-merge(this%rhoL,this%rhoG,pVF(i,j,k,1).gt.0.5_WP)*pV(i,j,k,1)
                     pFy(i,j,k,:)=mass_flux*0.5_WP*(pUVW(i,j-1,k,:)+pUVW(i,j,k,:))
                  end if
                  ! Velocity gradients at y-face
                  gradU(1,1)=0.25_WP*dxi*sum(pUVW(i+1,j-1:j,k,1)-pUVW(i-1,j-1:j,k,1))
                  gradU(2,1)=dyi*(pUVW(i,j,k,1)-pUVW(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*sum(pUVW(i,j-1:j,k+1,1)-pUVW(i,j-1:j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*sum(pUVW(i+1,j-1:j,k,2)-pUVW(i-1,j-1:j,k,2))
                  gradU(2,2)=dyi*(pUVW(i,j,k,2)-pUVW(i,j-1,k,2))
                  gradU(3,2)=0.25_WP*dzi*sum(pUVW(i,j-1:j,k+1,2)-pUVW(i,j-1:j,k-1,2))
                  gradU(1,3)=0.25_WP*dxi*sum(pUVW(i+1,j-1:j,k,3)-pUVW(i-1,j-1:j,k,3))
                  gradU(2,3)=dyi*(pUVW(i,j,k,3)-pUVW(i,j-1,k,3))
                  gradU(3,3)=0.25_WP*dzi*sum(pUVW(i,j-1:j,k+1,3)-pUVW(i,j-1:j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscosity at y-face
                  visc_f=2.0_WP*product(pVisc(i,j-1:j,k,1))/(sum(pVisc(i,j-1:j,k,1))+tiny(1.0_WP))
                  ! Viscous stress at y-face
                  pFy(i,j,k,1)=pFy(i,j,k,1)+visc_f*(gradU(1,2)+gradU(2,1))
                  pFy(i,j,k,2)=pFy(i,j,k,2)+visc_f*(gradU(2,2)+gradU(2,2))-2.0_WP/3.0_WP*visc_f*div
                  pFy(i,j,k,3)=pFy(i,j,k,3)+visc_f*(gradU(3,2)+gradU(2,3))
               end do; end do; end do
               ! Z-fluxes
               fbx=mfi%nodaltilebox(3)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Check if in band
                  if (lvl.eq.this%amr%clvl()) then; in_band=maxval(pBand(i,j,k-1:k,1)).gt.0.0_WP; else; in_band=.false.; end if
                  ! Outside band, compute finite volume Euler fluxes
                  if (.not.in_band) then
                     mass_flux=-merge(this%rhoL,this%rhoG,pVF(i,j,k,1).gt.0.5_WP)*pW(i,j,k,1)
                     pFz(i,j,k,:)=mass_flux*0.5_WP*(pUVW(i,j,k-1,:)+pUVW(i,j,k,:))
                  end if
                  ! Velocity gradients at z-face
                  gradU(1,1)=0.25_WP*dxi*sum(pUVW(i+1,j,k-1:k,1)-pUVW(i-1,j,k-1:k,1))
                  gradU(2,1)=0.25_WP*dyi*sum(pUVW(i,j+1,k-1:k,1)-pUVW(i,j-1,k-1:k,1))
                  gradU(3,1)=dzi*(pUVW(i,j,k,1)-pUVW(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*sum(pUVW(i+1,j,k-1:k,2)-pUVW(i-1,j,k-1:k,2))
                  gradU(2,2)=0.25_WP*dyi*sum(pUVW(i,j+1,k-1:k,2)-pUVW(i,j-1,k-1:k,2))
                  gradU(3,2)=dzi*(pUVW(i,j,k,2)-pUVW(i,j,k-1,2))
                  gradU(1,3)=0.25_WP*dxi*sum(pUVW(i+1,j,k-1:k,3)-pUVW(i-1,j,k-1:k,3))
                  gradU(2,3)=0.25_WP*dyi*sum(pUVW(i,j+1,k-1:k,3)-pUVW(i,j-1,k-1:k,3))
                  gradU(3,3)=dzi*(pUVW(i,j,k,3)-pUVW(i,j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscosity at z-face
                  visc_f=2.0_WP*product(pVisc(i,j,k-1:k,1))/(sum(pVisc(i,j,k-1:k,1))+tiny(1.0_WP))
                  ! Viscous stress at z-face
                  pFz(i,j,k,1)=pFz(i,j,k,1)+visc_f*(gradU(1,3)+gradU(3,1))
                  pFz(i,j,k,2)=pFz(i,j,k,2)+visc_f*(gradU(2,3)+gradU(3,2))
                  pFz(i,j,k,3)=pFz(i,j,k,3)+visc_f*(gradU(3,3)+gradU(3,3))-2.0_WP/3.0_WP*visc_f*div
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block finitevolume_fluxes

      ! Phase 2: Average down all fluxes for C/F conservation
      c_f_consistency: block
         use amrex_interface, only: amrmfab_average_down_face
         integer :: lvl
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         end do
      end block c_f_consistency

      ! Phase 3: Compute divergence and source terms for all levels, update VF/bary at band
      divergence_and_sources: block
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: lvl,i,j,k
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pCLold,pCGold
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pRhs,pFx,pFy,pFz,pVx,pVy,pVz,pBand
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold  ! Intentional masking
         real(WP) :: vol
         real(WP) :: Lvol_old,Lvol_new,Lvol_flux
         real(WP) :: Gvol_old,Gvol_new,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Lbar_new,Lbar_flux
         real(WP), dimension(3) :: Gbar_old,Gbar_new,Gbar_flux
         do lvl=0,this%amr%clvl()
            ! Grid spacings for this level
            dx=this%amr%dx(lvl); dxi=1.0_WP/dx
            dy=this%amr%dy(lvl); dyi=1.0_WP/dy
            dz=this%amr%dz(lvl); dzi=1.0_WP/dz
            vol=this%amr%cell_vol(lvl)
            ! Loop over tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               pRhs =>dmomdt%mf(lvl)%dataptr(mfi)
               pFx  =>Fx(lvl)%dataptr(mfi)
               pFy  =>Fy(lvl)%dataptr(mfi)
               pFz  =>Fz(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               ! Extra pointers at finest level
               if (lvl.eq.this%amr%maxlvl) then
                  pBand =>band%dataptr(mfi)
                  pVx   =>Vx%dataptr(mfi)
                  pVy   =>Vy%dataptr(mfi)
                  pVz   =>Vz%dataptr(mfi)
                  pVF   =>this%VF%mf(lvl)%dataptr(mfi)
                  pVFold=>this%VFold%mf(lvl)%dataptr(mfi)
                  pCL   =>this%CL%dataptr(mfi)
                  pCG   =>this%CG%dataptr(mfi)
                  pCLold=>this%CLold%dataptr(mfi)
                  pCGold=>this%CGold%dataptr(mfi)
               end if
               ! Loop over interior
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! VF/barycenter update at band cells (finest level only)
                  if (lvl.eq.this%amr%maxlvl) then
                     ! Work on band cells only
                     if (pBand(i,j,k,1).gt.0.0_WP) then
                        ! Old phasic moments
                        Lvol_old=(       pVFold(i,j,k,1))*vol
                        Gvol_old=(1.0_WP-pVFold(i,j,k,1))*vol
                        Lbar_old=pCLold(i,j,k,1:3)
                        Gbar_old=pCGold(i,j,k,1:3)
                        ! Net volume flux (outflow positive) from SL volume moments
                        Lvol_flux=pVx(i+1,j,k, 1 )-pVx(i,j,k, 1 )+pVy(i,j+1,k, 1 )-pVy(i,j,k, 1 )+pVz(i,j,k+1, 1 )-pVz(i,j,k, 1 )
                        Gvol_flux=pVx(i+1,j,k, 2 )-pVx(i,j,k, 2 )+pVy(i,j+1,k, 2 )-pVy(i,j,k, 2 )+pVz(i,j,k+1, 2 )-pVz(i,j,k, 2 )
                        Lbar_flux=pVx(i+1,j,k,3:5)-pVx(i,j,k,3:5)+pVy(i,j+1,k,3:5)-pVy(i,j,k,3:5)+pVz(i,j,k+1,3:5)-pVz(i,j,k,3:5)
                        Gbar_flux=pVx(i+1,j,k,6:8)-pVx(i,j,k,6:8)+pVy(i,j+1,k,6:8)-pVy(i,j,k,6:8)+pVz(i,j,k+1,6:8)-pVz(i,j,k,6:8)
                        ! New phasic volumes
                        Lvol_new=Lvol_old-Lvol_flux
                        Gvol_new=Gvol_old-Gvol_flux
                        ! New VF and default barycenters
                        pVF(i,j,k,1)=Lvol_new/max(Lvol_new+Gvol_new,tiny(1.0_WP))
                        pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                        pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                        ! Clip and update barycenters
                        if (pVF(i,j,k,1).lt.VFlo) then
                           pVF(i,j,k,1)=0.0_WP
                        else if (pVF(i,j,k,1).gt.VFhi) then
                           pVF(i,j,k,1)=1.0_WP
                        else
                           ! Update barycenters from moment conservation and project forward
                           if (Lvol_new/max(Lvol_new+Gvol_new,tiny(1.0_WP)).gt.vol_eps) then; Lbar_new=(Lbar_old*Lvol_old-Lbar_flux)/Lvol_new; pCL(i,j,k,1:3)=project(Lbar_new,dt); end if
                           if (Gvol_new/max(Lvol_new+Gvol_new,tiny(1.0_WP)).gt.vol_eps) then; Gbar_new=(Gbar_old*Gvol_old-Gbar_flux)/Gvol_new; pCG(i,j,k,1:3)=project(Gbar_new,dt); end if
                        end if
                     end if
                  end if
                  ! Divergence of momentum flux
                  pRhs(i,j,k,:)=dxi*(pFx(i+1,j,k,:)-pFx(i,j,k,:))+dyi*(pFy(i,j+1,k,:)-pFy(i,j,k,:))+dzi*(pFz(i,j,k+1,:)-pFz(i,j,k,:))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block divergence_and_sources

      ! Sync and apply BC
      call this%fill(lvl=this%amr%clvl(),time=time)

      ! Cleanup temporary mfabs
      cleanup: block
         integer :: lvl
         call this%amr%mfab_destroy(band)
         call this%amr%mfab_destroy(Vx)
         call this%amr%mfab_destroy(Vy)
         call this%amr%mfab_destroy(Vz)
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_destroy(Fx(lvl))
            call this%amr%mfab_destroy(Fy(lvl))
            call this%amr%mfab_destroy(Fz(lvl))
         end do
      end block cleanup

   contains

      !> Recursive subroutine that cuts a tet by grid planes to compute volume and Q fluxes
      recursive subroutine tet2flux(mytet,myind,myVflux,myQflux)
         use amrvof_geometry, only: cut_side,cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert
         implicit none
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  dimension(3,4), intent(in) :: myind
         real(WP), dimension(8),  intent(out) :: myVflux
         real(WP), dimension(3),  intent(out) :: myQflux
         integer :: dir,cut_ind,icase,n1,n2,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         integer,  dimension(3,8,2) :: vert_ind
         real(WP) :: mu,my_vol
         real(WP), dimension(3,4) :: newtet
         integer,  dimension(3,4) :: newind
         real(WP), dimension(3) :: a,b,c
         real(WP), dimension(8) :: subVflux
         real(WP), dimension(3) :: subQflux
         real(WP) :: xcut,ycut,zcut
         
         myVflux=0.0_WP
         myQflux=0.0_WP
         
         ! Determine if tet spans multiple cells and needs cutting
         if (maxval(myind(1,:))-minval(myind(1,:)).gt.0) then
            dir=1; cut_ind=maxval(myind(1,:))
            xcut=this%amr%xlo+real(cut_ind,WP)*dx
            dd(:)=mytet(1,:)-xcut
         else if (maxval(myind(2,:))-minval(myind(2,:)).gt.0) then
            dir=2; cut_ind=maxval(myind(2,:))
            ycut=this%amr%ylo+real(cut_ind,WP)*dy
            dd(:)=mytet(2,:)-ycut
         else if (maxval(myind(3,:))-minval(myind(3,:)).gt.0) then
            dir=3; cut_ind=maxval(myind(3,:))
            zcut=this%amr%zlo+real(cut_ind,WP)*dz
            dd(:)=mytet(3,:)-zcut
         else
            ! All vertices in same cell - cut by PLIC and return
            call tet2flux_plic(mytet,myind(1,1),myind(2,1),myind(3,1),myVflux,myQflux)
            return
         end if
         
         ! Find cut case (1-indexed: 1-16)
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
         &    +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
         &    +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
         &    +8*int(0.5_WP+sign(0.5_WP,dd(4)))
         
         ! Copy vertices and indices
         do n1=1,4
            vert(:,n1)=mytet(:,n1)
            vert_ind(:,n1,1)=myind(:,n1)
            vert_ind(:,n1,2)=myind(:,n1)
            vert_ind(dir,n1,1)=min(vert_ind(dir,n1,1),cut_ind-1)
            vert_ind(dir,n1,2)=max(vert_ind(dir,n1,1),cut_ind)
         end do
         
         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
            vert_ind(1,4+n1,1)=floor((vert(1,4+n1)-this%amr%xlo)*dxi)
            vert_ind(2,4+n1,1)=floor((vert(2,4+n1)-this%amr%ylo)*dyi)
            vert_ind(3,4+n1,1)=floor((vert(3,4+n1)-this%amr%zlo)*dzi)
            vert_ind(:,4+n1,1)=max(vert_ind(:,4+n1,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,1)=min(vert_ind(:,4+n1,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,2)=vert_ind(:,4+n1,1)
            vert_ind(dir,4+n1,1)=cut_ind-1
            vert_ind(dir,4+n1,2)=cut_ind
         end do
         
         ! Create and process sub-tets
         do n1=1,cut_ntets(icase)
            do n2=1,4
               newtet(:,n2)=vert(:,cut_vtet(n2,n1,icase))
               newind(:,n2)=vert_ind(:,cut_vtet(n2,n1,icase),cut_side(n1,icase))
            end do
            a=newtet(:,1)-newtet(:,4)
            b=newtet(:,2)-newtet(:,4)
            c=newtet(:,3)-newtet(:,4)
            my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            if (my_vol.lt.VFlo*dx*dy*dz) cycle
            call tet2flux(newtet,newind,subVflux,subQflux)
            myVflux=myVflux+subVflux
            myQflux=myQflux+subQflux
         end do
         
      end subroutine tet2flux

      !> Cut tet by PLIC and compute volume + conserved variable fluxes
      subroutine tet2flux_plic(mytet,i0,j0,k0,myVflux,myQflux)
         use amrvof_geometry, only: cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert,cut_nntet,tet_vol
         !use messager, only: die
         implicit none
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  intent(in) :: i0,j0,k0
         real(WP), dimension(8),  intent(out) :: myVflux
         real(WP), dimension(3),  intent(out) :: myQflux
         integer :: icase,n1,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         real(WP), dimension(3) :: a,b,c,bary,normal,bary_tot
         real(WP) :: mu,my_vol,dist,VF0,vol_tot

         ! Zero out flux arrays
         myVflux=0.0_WP
         myQflux=0.0_WP

         ! Check indices are within PLICold bounds
         !if (i0.lt.lbound(pPLICold,1).or.i0.gt.ubound(pPLICold,1).or. &
         !    j0.lt.lbound(pPLICold,2).or.j0.gt.ubound(pPLICold,2).or. &
         !    k0.lt.lbound(pPLICold,3).or.k0.gt.ubound(pPLICold,3)) then
         !   call die('[tet2flux_plic] Index out of bounds - check CFL or ghost cells')
         !end if
         
         ! Get old VF for this cell
         VF0=pVFold(i0,j0,k0,1)
         
         ! Tet volume and barycenter
         vol_tot=abs(tet_vol(mytet))
         bary_tot=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
         
         ! Pure cell shortcut
         if (pPLICold(i0,j0,k0,4).gt.+1.0e9_WP) then
            ! Pure liquid
            myVflux( 1 )=vol_tot
            myVflux(3:5)=vol_tot*bary_tot
            ! Q flux: pure liquid momentum
            myQflux=vol_tot*this%rhoL*pUVWold(i0,j0,k0,1:3)
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas
            myVflux( 2 )=vol_tot
            myVflux(6:8)=vol_tot*bary_tot
            ! Q flux: pure gas momentum
            myQflux=vol_tot*this%rhoG*pUVWold(i0,j0,k0,1:3)
            return
         end if

         ! If we get here, we ARE cutting by a PLIC plane
         crossed_plic=.true.
         
         ! Get PLIC from this cell
         normal=pPLICold(i0,j0,k0,1:3)
         dist=pPLICold(i0,j0,k0,4)
         
         ! Compute signed distance to plane for each vertex
         dd(1)=normal(1)*mytet(1,1)+normal(2)*mytet(2,1)+normal(3)*mytet(3,1)-dist
         dd(2)=normal(1)*mytet(1,2)+normal(2)*mytet(2,2)+normal(3)*mytet(3,2)-dist
         dd(3)=normal(1)*mytet(1,3)+normal(2)*mytet(2,3)+normal(3)*mytet(3,3)-dist
         dd(4)=normal(1)*mytet(1,4)+normal(2)*mytet(2,4)+normal(3)*mytet(3,4)-dist
         
         ! Find cut case
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
         &    +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
         &    +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
         &    +8*int(0.5_WP+sign(0.5_WP,dd(4)))
         
         ! Copy vertices
         vert(:,1:4)=mytet(:,1:4)
         
         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
         end do

         ! Cut the minority phase (safer as we subtract small from large)
         if (VF0.gt.0.5_WP) then
            ! Liquid is dominant → compute gas directly
            do n1=1,cut_nntet(icase)-1
               a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
               bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
               &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
               myVflux( 2 )=myVflux( 2 )+my_vol
               myVflux(6:8)=myVflux(6:8)+my_vol*bary
            end do
            ! Liquid = total - gas
            myVflux( 1 )=vol_tot-myVflux( 2 )
            myVflux(3:5)=vol_tot*bary_tot-myVflux(6:8)
         else
            ! Gas is dominant → compute liquid directly
            do n1=cut_ntets(icase),cut_nntet(icase),-1
               a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
               bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
               &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
               myVflux( 1 )=myVflux( 1 )+my_vol
               myVflux(3:5)=myVflux(3:5)+my_vol*bary
            end do
            ! Gas = total - liquid
            myVflux( 2 )=vol_tot-myVflux( 1 )
            myVflux(6:8)=vol_tot*bary_tot-myVflux(3:5)
         end if

         ! Compute Q flux from UVWold
         myQflux=(this%rhoL*myVflux(1)+this%rhoG*myVflux(2))*pUVWold(i0,j0,k0,1:3)

      end subroutine tet2flux_plic

      !> RK2 vertex projection back in time
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project

      !> Trilinear interpolation of staggered velocity - uses pU,pV,pW
      function interp_velocity(pos) result(vel)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer  :: ipc,jpc,kpc   ! Cell-centered indices
         integer  :: ipu,jpv,kpw   ! Face-centered indices
         real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2  ! Cell-centered weights
         real(WP) :: wxu1,wyv1,wzw1,wxu2,wyv2,wzw2  ! Face-centered weights
         ! Compute raw indices
         ipc=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP)
         jpc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP)
         kpc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP)
         ipu=floor((pos(1)-this%amr%xlo)*dxi)
         jpv=floor((pos(2)-this%amr%ylo)*dyi)
         kpw=floor((pos(3)-this%amr%zlo)*dzi)
         ! Clamp to array bounds
         !ipu=max(lbound(pU,1),min(ubound(pU,1)-1,ipu))
         !jpc=max(lbound(pU,2),min(ubound(pU,2)-1,jpc))
         !kpc=max(lbound(pU,3),min(ubound(pU,3)-1,kpc))
         !ipc=max(lbound(pV,1),min(ubound(pV,1)-1,ipc))
         !jpv=max(lbound(pV,2),min(ubound(pV,2)-1,jpv))
         !kpw=max(lbound(pW,3),min(ubound(pW,3)-1,kpw))
         ! Cell-centered weights
         wxc1=(pos(1)-(this%amr%xlo+(real(ipc,WP)+0.5_WP)*dx))*dxi
         wyc1=(pos(2)-(this%amr%ylo+(real(jpc,WP)+0.5_WP)*dy))*dyi
         wzc1=(pos(3)-(this%amr%zlo+(real(kpc,WP)+0.5_WP)*dz))*dzi
         wxc1=max(0.0_WP,min(1.0_WP,wxc1)); wxc2=1.0_WP-wxc1
         wyc1=max(0.0_WP,min(1.0_WP,wyc1)); wyc2=1.0_WP-wyc1
         wzc1=max(0.0_WP,min(1.0_WP,wzc1)); wzc2=1.0_WP-wzc1
         ! Face-centered weights
         wxu1=(pos(1)-(this%amr%xlo+real(ipu,WP)*dx))*dxi
         wyv1=(pos(2)-(this%amr%ylo+real(jpv,WP)*dy))*dyi
         wzw1=(pos(3)-(this%amr%zlo+real(kpw,WP)*dz))*dzi
         wxu1=max(0.0_WP,min(1.0_WP,wxu1)); wxu2=1.0_WP-wxu1
         wyv1=max(0.0_WP,min(1.0_WP,wyv1)); wyv2=1.0_WP-wyv1
         wzw1=max(0.0_WP,min(1.0_WP,wzw1)); wzw2=1.0_WP-wzw1
         ! U at x-faces: face-centered in x, cell-centered in y,z
         vel(1)=wzc1*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc+1,1)+wxu2*pU(ipu,jpc+1,kpc+1,1)) +&
         &            wyc2*(wxu1*pU(ipu+1,jpc  ,kpc+1,1)+wxu2*pU(ipu,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc  ,1)+wxu2*pU(ipu,jpc+1,kpc  ,1)) +&
         &            wyc2*(wxu1*pU(ipu+1,jpc  ,kpc  ,1)+wxu2*pU(ipu,jpc  ,kpc  ,1)))
         ! V at y-faces: cell-centered in x, face-centered in y, cell-centered in z
         vel(2)=wzc1*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc+1,1)+wxc2*pV(ipc,jpv+1,kpc+1,1)) +&
         &            wyv2*(wxc1*pV(ipc+1,jpv  ,kpc+1,1)+wxc2*pV(ipc,jpv  ,kpc+1,1)))+&
         &      wzc2*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc  ,1)+wxc2*pV(ipc,jpv+1,kpc  ,1)) +&
         &            wyv2*(wxc1*pV(ipc+1,jpv  ,kpc  ,1)+wxc2*pV(ipc,jpv  ,kpc  ,1)))
         ! W at z-faces: cell-centered in x,y, face-centered in z
         vel(3)=wzw1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw+1,1)+wxc2*pW(ipc,jpc+1,kpw+1,1)) +&
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpw+1,1)+wxc2*pW(ipc,jpc  ,kpw+1,1)))+&
         &      wzw2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw  ,1)+wxc2*pW(ipc,jpc+1,kpw  ,1)) +&
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpw  ,1)+wxc2*pW(ipc,jpc  ,kpw  ,1)))
      end function interp_velocity

   end subroutine get_dmomdt

   !> Add Vreman SGS eddy viscosity to this%visc
   !> Assumes velocity ghosts are filled. User must reset visc
   !> to molecular value before calling this routine.
   subroutine add_vreman(this,dt,Cs)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: visc_t
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pUVW,pVisc,pVF
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,Cmodel,Aij,Bij
      real(WP), dimension(1:3,1:3) :: gradU,betaij
      integer :: lvl,i,j,k,si,sj
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2

      ! Model constant: c=2.5*Cs**2 (Vreman uses c=0.07 ~ Cs=0.17)
      if (present(Cs)) then; Cmodel=2.5_WP*Cs**2; else; Cmodel=2.5_WP*0.17_WP**2; end if

      ! Loop over levels
      do lvl=0,this%amr%clvl()
         ! Grid spacings
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz
         ! Max visc from CFL
         max_visc=max_cfl*this%amr%min_meshsize(lvl)**2/(4.0_WP*dt)
         ! Build temp multifab for eddy viscosity (1 ghost for filtering)
         call this%amr%mfab_build(lvl=lvl,mfab=visc_t,ncomp=1,nover=1); call visc_t%setval(0.0_WP)
         ! Phase 1: Compute kinematic eddy viscosity
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            pVisc_t=>visc_t%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute cell-centered velocity gradient tensor
               gradU(1,1)=0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
               gradU(2,1)=0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
               gradU(3,1)=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
               gradU(1,2)=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
               gradU(2,2)=0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
               gradU(3,2)=0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
               gradU(1,3)=0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
               gradU(2,3)=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
               gradU(3,3)=0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
               ! A=gradU_ij*gradU_ij invariant
               Aij=sum(gradU**2)
               ! beta_ij=dx_m^2*gradU_mi*gradU_mj
               do sj=1,3; do si=1,3
                  betaij(si,sj)=dx**2*gradU(1,si)*gradU(1,sj)+dy**2*gradU(2,si)*gradU(2,sj)+dz**2*gradU(3,si)*gradU(3,sj)
               end do; end do
               ! B invariant
               Bij=betaij(1,1)*betaij(2,2)-betaij(1,2)**2+betaij(1,1)*betaij(3,3)-betaij(1,3)**2+betaij(2,2)*betaij(3,3)-betaij(2,3)**2
               ! Assemble eddy viscosity
               if (Bij.gt.0.0_WP) pVisc_t(i,j,k,1)=Cmodel*sqrt(Bij/Aij)
               ! Clip to CFL limit
               pVisc_t(i,j,k,1)=min(pVisc_t(i,j,k,1),max_visc)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Phase 2: Filter visc_t
         call this%amr%mfab_filter(lvl=lvl,mfab=visc_t,npass=nfilter)
         ! Phase 3: Convert to dynamic viscosity and add to this%visc
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)/(pVF(i,j,k,1)/this%rhoL+(1.0_WP-pVF(i,j,k,1))/this%rhoG)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Destroy temp multifab
         call amrex_multifab_destroy(visc_t)
      end do
   end subroutine add_vreman

   !> Compute CFL numbers (convective and viscous)
   subroutine get_cfl(this,dt,cfl,cflc)
      use mathtools, only: Pi
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP), intent(out), optional :: cflc
      integer :: lvl
      real(WP) :: Umax_lvl,Vmax_lvl,Wmax_lvl,viscmax
      ! Reset CFLs
      this%CFLc_x=0.0_WP; this%CFLc_y=0.0_WP; this%CFLc_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      this%CFLst=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         Umax_lvl=this%U%norm0(lvl=lvl)
         Vmax_lvl=this%V%norm0(lvl=lvl)
         Wmax_lvl=this%W%norm0(lvl=lvl)
         ! Max viscosity
         get_viscmax: block
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            use parallel, only: MPI_REAL_WP
            use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc,pVF
            integer :: i,j,k,ierr
            viscmax=0.0_WP
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               ! Get data pointers
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               ! Loop over interior tiles
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  viscmax=max(viscmax,pVisc(i,j,k,1)/(pVF(i,j,k,1)*this%rhoL+(1.0_WP-pVF(i,j,k,1))*this%rhoG))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            call MPI_ALLREDUCE(MPI_IN_PLACE,viscmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         end block get_viscmax
         ! Convective CFL
         if (this%amr%nx.gt.1) this%CFLc_x=max(this%CFLc_x,dt*Umax_lvl/this%amr%dx(lvl))
         if (this%amr%ny.gt.1) this%CFLc_y=max(this%CFLc_y,dt*Vmax_lvl/this%amr%dy(lvl))
         if (this%amr%nz.gt.1) this%CFLc_z=max(this%CFLc_z,dt*Wmax_lvl/this%amr%dz(lvl))
         ! Viscous CFL (explicit stability: dt < dx^2 / (4*nu))
         if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt/this%amr%dx(lvl)**2)
         if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt/this%amr%dy(lvl)**2)
         if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt/this%amr%dz(lvl)**2)
      end do
      ! Surface-tension CFL (capillary wave stability criterion)
      if (this%sigma.gt.0.0_WP.and.this%amr%clvl().eq.this%amr%maxlvl) this%CFLst=dt/sqrt((this%rhoL+this%rhoG)*this%amr%min_meshsize(this%amr%maxlvl)**3/(4.0_WP*Pi*this%sigma))
      ! Compute max overall CFL
      this%CFL=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z,this%CFLst)
      ! Return max overall CFL
      cfl=this%CFL
      ! Optionally return max convective CFL
      if (present(cflc)) cflc=max(this%CFLc_x,this%CFLc_y,this%CFLc_z)
   end subroutine get_cfl

   !> Correct outflow velocity to ensure global mass conservation
   !> Scans all 6 domain faces: ext_dir faces contribute fixed flux,
   !> foextrap faces are correctable. Correction is distributed uniformly
   !> over all foextrap faces, optionally weighted by VF (fluid volume fraction).
   !> Uses composite integration with fine_mask to avoid double-counting across AMR levels.
   subroutine correct_outflow(this,VF)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      use messager, only: die
      use amrex_amr_module, only: amrex_mfiter,amrex_bc_foextrap,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmask_make_fine
      implicit none
      class(amrmpinc), intent(inout) :: this
      class(amrdata), intent(in), optional :: VF
      ! Face classification: ftype(dir,side) where dir=1,2,3 and side=1(lo),2(hi)
      integer, parameter :: SKIP=0,FIXED=1,CORR=2,INTEGRATE=1,CORRECT=2
      integer, dimension(3,2) :: ftype
      logical, dimension(3) :: per
      logical :: has_corr
      ! Composite integration
      real(WP) :: Qflux,Aout,Ucorr,dA,VFf
      integer :: lvl,i,j,k,dir,side,ierr,dlo(3),dhi(3),bc
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pVF
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      ! Sanity check
      if (present(VF)) then; if (VF%ng.lt.1) call die('[amrmpinc correct_outflow] VF must have at least 1 ghost cell'); end if
      ! Classify faces by dir and side
      per=[this%amr%xper,this%amr%yper,this%amr%zper]
      ftype=SKIP; has_corr=.false.
      do dir=1,3; do side=1,2
         if (per(dir)) cycle
         ! Get BC type on normal component for this dir/side
         select case (dir)
          case (1); bc=merge(this%U%lo_bc(1,1),this%U%hi_bc(1,1),side.eq.1)
          case (2); bc=merge(this%V%lo_bc(2,1),this%V%hi_bc(2,1),side.eq.1)
          case (3); bc=merge(this%W%lo_bc(3,1),this%W%hi_bc(3,1),side.eq.1)
         end select
         ftype(dir,side)=merge(CORR,FIXED,bc.eq.amrex_bc_foextrap)
         if (ftype(dir,side).eq.CORR) has_corr=.true.
      end do; end do
      if (.not.has_corr) return
      ! Pass 1: Integrate all face fluxes
      Qflux=0.0_WP; Aout=0.0_WP
      call composite_loop(INTEGRATE)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Qflux,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Aout, 1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      ! Compute correction velocity
      if (Aout.le.0.0_WP) return
      Ucorr=-Qflux/Aout
      ! Pass 2: Apply correction to foextrap faces
      call composite_loop(CORRECT)

   contains

      !> Composite loop over all levels with fine masking
      !> mode=INTEGRATE: accumulate Qflux/Aout; mode=CORRECT: apply Ucorr
      subroutine composite_loop(mode)
         implicit none
         integer, intent(in) :: mode
         do lvl=0,this%amr%clvl()
            dlo=this%amr%geom(lvl)%domain%lo
            dhi=this%amr%geom(lvl)%domain%hi
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
            do while (mfi%next())
               bx=mfi%tilebox()
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               if (present(VF)) pVF=>VF%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               do dir=1,3; do side=1,2
                  if (mode.eq.INTEGRATE.and.ftype(dir,side).eq.SKIP) cycle
                  if (mode.eq.CORRECT.and.ftype(dir,side).ne.CORR) cycle
                  call process_face(dir,side,mode)
               end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
      end subroutine composite_loop

      !> Process a single domain boundary face
      !> dir=1,2,3 (x,y,z); side=1(lo),2(hi)
      !> mode=INTEGRATE: accumulate flux into Qflux/Aout
      !> mode=CORRECT: apply Ucorr correction to velocity
      subroutine process_face(dir,side,mode)
         implicit none
         integer, intent(in) :: dir,side,mode
         integer :: bnd,ci
         real(WP) :: sgn
         ! Outward normal sign: lo face -> -1, hi face -> +1
         sgn=merge(-1.0_WP,1.0_WP,side.eq.1)
         ! Staggered boundary index and adjacent interior cell
         if (side.eq.1) then; bnd=dlo(dir); ci=dlo(dir); else; bnd=dhi(dir)+1; ci=dhi(dir); end if
         ! Check tile ownership (tile must contain the interior cell)
         if (bx%lo(dir).gt.ci.or.bx%hi(dir).lt.ci) return
         ! Face area (product of transverse mesh spacings)
         select case (dir)
          case (1); dA=this%amr%dy(lvl)*this%amr%dz(lvl)
          case (2); dA=this%amr%dx(lvl)*this%amr%dz(lvl)
          case (3); dA=this%amr%dx(lvl)*this%amr%dy(lvl)
         end select
         ! Loop over transverse indices
         select case (dir)
          case (1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
               if (lvl.lt.this%amr%clvl()) then; if (pMask(ci,j,k,1).eq.0) cycle; end if
               VFf=1.0_WP; if (present(VF)) VFf=0.5_WP*sum(pVF(bnd-1:bnd,j,k,1))
               if (mode.eq.INTEGRATE) then
                  Qflux=Qflux+sgn*pU(bnd,j,k,1)*VFf*dA
                  if (ftype(dir,side).eq.CORR) Aout=Aout+VFf*dA
               else
                  pU(bnd,j,k,1)=pU(bnd,j,k,1)+sgn*Ucorr*VFf
               end if
            end do; end do
          case (2)
            do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
               if (lvl.lt.this%amr%clvl()) then; if (pMask(i,ci,k,1).eq.0) cycle; end if
               VFf=1.0_WP; if (present(VF)) VFf=0.5_WP*sum(pVF(i,bnd-1:bnd,k,1))
               if (mode.eq.INTEGRATE) then
                  Qflux=Qflux+sgn*pV(i,bnd,k,1)*VFf*dA
                  if (ftype(dir,side).eq.CORR) Aout=Aout+VFf*dA
               else
                  pV(i,bnd,k,1)=pV(i,bnd,k,1)+sgn*Ucorr*VFf
               end if
            end do; end do
          case (3)
            do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,ci,1).eq.0) cycle; end if
               VFf=1.0_WP; if (present(VF)) VFf=0.5_WP*sum(pVF(i,j,bnd-1:bnd,1))
               if (mode.eq.INTEGRATE) then
                  Qflux=Qflux+sgn*pW(i,j,bnd,1)*VFf*dA
                  if (ftype(dir,side).eq.CORR) Aout=Aout+VFf*dA
               else
                  pW(i,j,bnd,1)=pW(i,j,bnd,1)+sgn*Ucorr*VFf
               end if
            end do; end do
         end select
      end subroutine process_face

   end subroutine correct_outflow

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Get solver information: min/max velocity, min/max pressure, divergence, momentum, TKE
   subroutine get_info(this)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer :: lvl

      ! Get VF info from parent
      call this%amrvof%get_info()

      ! First compute divergence (this updates divmax)
      call this%get_div()

      ! Initialize min/max values
      this%Umax=-huge(1.0_WP)
      this%Vmax=-huge(1.0_WP)
      this%Wmax=-huge(1.0_WP)
      this%Pmax=-huge(1.0_WP)

      ! Loop over all levels for min/max
      do lvl=0,this%amr%clvl()
         this%Umax=max(this%Umax,this%U%norm0(lvl=lvl),this%UVW%norm0(lvl=lvl,comp=1))
         this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl),this%UVW%norm0(lvl=lvl,comp=2))
         this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl),this%UVW%norm0(lvl=lvl,comp=3))
         this%Pmax=max(this%Pmax,this%P%norm0(lvl=lvl))
      end do

      ! Integrate momentum and kinetic energy
      get_momentum_and_kinetic_energy: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         integer :: i,j,k,ierr
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW,pVF
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP) :: rho
         ! Uses composite integration with fine masking to avoid double-counting
         this%rhoUint=0.0_WP
         this%rhoVint=0.0_WP
         this%rhoWint=0.0_WP
         this%rhoKint=0.0_WP
         do lvl=0,this%amr%clvl()
            ! Build fine mask for this level (if not finest)
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            ! Loop over all cells
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over tile
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then
                     if (pMask(i,j,k,1).eq.0) cycle
                  end if
                  ! Get local density
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  ! Accumulate momentum
                  this%rhoUint=this%rhoUint+rho*pUVW(i,j,k,1)*this%amr%cell_vol(lvl)
                  this%rhoVint=this%rhoVint+rho*pUVW(i,j,k,2)*this%amr%cell_vol(lvl)
                  this%rhoWint=this%rhoWint+rho*pUVW(i,j,k,3)*this%amr%cell_vol(lvl)
                  ! Accumulate kinetic energy
                  this%rhoKint=this%rhoKint+0.5_WP*rho*(pUVW(i,j,k,1)**2+pUVW(i,j,k,2)**2+pUVW(i,j,k,3)**2)*this%amr%cell_vol(lvl)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoUint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoVint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoWint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_momentum_and_kinetic_energy

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrmpinc_print(this)
      use messager, only: log
      use string, only: str_long
      implicit none
      class(amrmpinc), intent(in) :: this
      character(len=str_long) :: message
      call log("Incompressible multiphase collocated solver: "//trim(this%name))
      write(message,'("  rhoL = ",ES12.5)') this%rhoL; call log(trim(message))
      write(message,'("  rhoG = ",ES12.5)') this%rhoG; call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrmpinc_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrmpinc), intent(inout) :: this
      class(amrio), intent(inout) :: io
      ! VOF data is registered via parent
      call this%amrvof%register_checkpoint(io)
      ! Register flow data
      call io%add_data(this%UVW,'UVW')
      call io%add_data(this%U,'U')
      call io%add_data(this%V,'V')
      call io%add_data(this%W,'W')
      call io%add_data(this%P,'P')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrmpinc), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      ! VOF data is restored via parent
      call this%amrvof%restore_checkpoint(io,dirname,time)
      ! Restore flow data
      call io%read_data(dirname,this%UVW,'UVW')
      call io%read_data(dirname,this%U,'U')
      call io%read_data(dirname,this%V,'V')
      call io%read_data(dirname,this%W,'W')
      call io%read_data(dirname,this%P,'P')
      ! Fill ghost cells (VisMF reads valid data only)
      call this%UVW%fill(time=time)
      call this%fill_velocity(time=time)
      call this%P%fill(time=time)
   end subroutine restore_checkpoint

end module amrmpinc_class
