!> AMR Incompressible solver class
!> Provides projection and basic operations for constant-density flow
module amrincomp_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrmg_class,      only: amrmg
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_interp_face_divfree
   implicit none
   private

   ! Expose type
   public :: amrincomp

   !> AMR Incompressible solver type
   type, extends(amrsolver) :: amrincomp
      ! User-configurable callbacks
      procedure(incomp_init_iface),    pointer, pass :: user_init   =>null()
      procedure(incomp_tagging_iface), pointer, pass :: user_tagging=>null()
      procedure(incomp_bc_iface),      pointer, pass :: user_bc     =>null()

      ! Flow data
      type(amrdata) :: U,V,W            !< Current face velocities
      type(amrdata) :: Uold,Vold,Wold   !< Old face velocities
      type(amrdata) :: P                !< Pressure (cell-centered)
      type(amrdata) :: div              !< Divergence (cell-centered)

      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rho=1.0_WP            !< Constant density
      type(amrdata) :: visc             !< variable dynamic viscosity

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
      real(WP) :: CFL=0.0_WP            !< Maximum CFL

      ! Overlap size
      integer :: nover=1

      ! Velocity interpolation method (for C/F interface fills)
      integer :: interp_vel=amrex_interp_face_divfree

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
      ! Staggered velocity fills
      procedure :: fill_velocity_lvl         !< Fill velocity ghosts at single level
      procedure :: fill_velocity             !< Fill velocity ghosts on all levels
      procedure :: fill_velocity_from_coarse !< Fill velocity from coarse
      procedure :: sync_velocity_lvl         !< Sync velocity ghosts at single level
      procedure :: sync_velocity             !< Sync velocity ghosts on all levels
      procedure :: fill_velocity_mfab        !< Fill dest MultiFabs for regridding
      procedure :: average_down_velocity     !< Average down MAC velocity for C/F consistency
      procedure :: average_down_velocity_to  !< Average down MAC velocity for single level
      ! Utilities
      procedure :: get_div                   !< Compute divergence (assumes velocity ghosts filled)
      procedure :: correct_velocity          !< Correct face velocity with pressure gradient
      ! Physics procedures
      procedure :: get_dmomdt                !< Compute momentum advection RHS
      procedure :: add_vreman                !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                   !< Compute CFL numbers
      procedure :: correct_outflow           !< Correct outflow for global mass conservation
      ! Print solver info
      procedure :: get_info
      procedure :: print => amrincomp_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrincomp

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine incomp_init_iface(solver,lvl,time,ba,dm)
         import :: amrincomp,WP,amrex_boxarray,amrex_distromap
         class(amrincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine incomp_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine incomp_tagging_iface(solver,lvl,time,tags)
         import :: amrincomp,c_ptr,WP
         class(amrincomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine incomp_tagging_iface
   end interface

   !> Abstract interface for user-provided velocity BC callback
   !> Called for ext_dir faces; user fills the boundary box with their own values
   abstract interface
      subroutine incomp_bc_iface(solver,lvl,time,face,bx,comp,p)
         import :: amrincomp,amrex_box,WP
         class(amrincomp), intent(in) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< 'U', 'V', or 'W'
         real(WP), dimension(:,:,:,:), pointer, intent(inout) :: p
      end subroutine incomp_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrincomp type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrincomp_on_init(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(lvl,time,ba,dm)
   end subroutine amrincomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrincomp_on_coarse(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrincomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrincomp_on_remake(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrincomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrincomp_on_clear(ctx,lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrincomp_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrincomp_tagging(ctx,lvl,time,tags)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx,this)
      if (associated(this%user_tagging)) call this%user_tagging(lvl,time,tags)
   end subroutine amrincomp_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrincomp_postregrid(ctx,lbase,time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrincomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrincomp_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the incompressible solver
   subroutine initialize(this,amr,name)
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrmg_class,      only: amrmg_cstcoef
      implicit none
      class(amrincomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name=trim(name)
      else
         this%name='UNNAMED_INCOMP'
      end if

      ! Store amrgrid pointer
      this%amr=>amr

      ! Initialize staggered velocity
      call this%U%initialize(amr,name='U',ncomp=1,ng=this%nover,nodal=[.true. ,.false.,.false.]); this%U%parent=>this
      call this%V%initialize(amr,name='V',ncomp=1,ng=this%nover,nodal=[.false.,.true. ,.false.]); this%V%parent=>this
      call this%W%initialize(amr,name='W',ncomp=1,ng=this%nover,nodal=[.false.,.false.,.true. ]); this%W%parent=>this
      call this%Uold%initialize(amr,name='Uold',ncomp=1,ng=this%nover,nodal=[.true. ,.false.,.false.]); this%Uold%parent=>this
      call this%Vold%initialize(amr,name='Vold',ncomp=1,ng=this%nover,nodal=[.false.,.true. ,.false.]); this%Vold%parent=>this
      call this%Wold%initialize(amr,name='Wold',ncomp=1,ng=this%nover,nodal=[.false.,.false.,.true. ]); this%Wold%parent=>this

      ! Set velocity fillbc callbacks to shared handler
      this%U%fillbc=>velocity_fillbc
      this%V%fillbc=>velocity_fillbc
      this%W%fillbc=>velocity_fillbc

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
      call this%psolver%initialize(amr,type=amrmg_cstcoef)

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      select type (this)
       type is (amrincomp)
         call this%amr%add_on_init   (amrincomp_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrincomp_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrincomp_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrincomp_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrincomp_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrincomp_postregrid,c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the incompressible solver
   subroutine finalize(this)
      implicit none
      class(amrincomp), intent(inout) :: this
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
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
      nullify(this%user_bc)
   end subroutine finalize

   ! ============================================================================
   ! LIFECYCLE CALLBACKS
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
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
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%div%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using divergence-free interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Velocity: allocate then fill with divergence-free interpolation
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%fill_velocity_from_coarse(lvl,time)
      ! Old velocity just needs geometry
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      ! Pressure on coarse
      call this%P%on_coarse(lvl,time,ba,dm)
      ! Divergence and viscosity just need geometry
      call this%div%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
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
      ! Old velocity just needs geometry
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
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
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
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Average down velocity and pressure
      call this%average_down_velocity(lbase)
      call this%P%average_down(lbase)
      ! Fill ghosts
      call this%fill_velocity(time,lbase)
      call this%P%fill(time,lbase)
      ! Rebuild pressure solver operators for new grid
      call this%psolver%setup()
   end subroutine post_regrid

   ! ============================================================================
   ! Staggered velocity fills
   ! ============================================================================

   !> Average down MAC velocity for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles face-centered averaging correctly
   subroutine average_down_velocity_to(this,lvl)
      implicit none
      class(amrincomp), intent(inout) :: this
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
      class(amrincomp), intent(inout) :: this
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
      class(amrincomp), target, intent(inout) :: this
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
      class(amrincomp), intent(inout) :: this
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
      class(amrincomp), target, intent(inout) :: this
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
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%sync_lvl(lvl)
      call this%V%sync_lvl(lvl)
      call this%W%sync_lvl(lvl)
   end subroutine sync_velocity_lvl

   !> Sync velocity ghost cells on all levels
   subroutine sync_velocity(this,lbase)
      implicit none
      class(amrincomp), intent(inout) :: this
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
      class(amrincomp), target, intent(inout) :: this
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
      class(amrincomp), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: dlo(3),dhi(3),flo(3),fhi(3)
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: lvl
      character(len=1) :: comp

      ! Get point to solver
      select type (s=>this%parent)
       class is (amrincomp)
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
            if (associated(solver%user_bc)) then
               bc_bx=amrex_box(slo,shi)
               call solver%user_bc(lvl=lvl,time=time,face=face,bx=bc_bx,comp=comp,p=p)
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
   ! UTILITIES
   ! ============================================================================

   !> Compute divergence of velocity into internal div field, update divmax
   !> Uses composite fine masking so covered coarse cells don't pollute divmax
   subroutine get_div(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmfab_compute_divergence,amrmask_make_fine
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrincomp), intent(inout) :: this
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

   !> Add (-scale*pressure gradient) to U/V/W face velocities. Two paths:
   !>   phi present -> direct path: use explicit stencil that reads phi ghost cells directly (for predictor with fs%P)
   !>   phi absent  -> MLMG path:   use psolver internal fluxes (use for projection with dP)
   subroutine correct_velocity(this,scale,phi)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrdata), intent(in), optional :: phi
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP,pFx,pFy,pFz
      real(WP) :: dxi,dyi,dzi
      integer :: lvl,i,j,k
      ! Build temp face mfabs to store pressure gradient
      allocate(Fx(0:this%amr%clvl()),Fy(0:this%amr%clvl()),Fz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Compute -pressure gradient at faces
      if (present(phi)) then
         ! Direct path: differentiate provided phi, using its ghosts cells
         do lvl=0,this%amr%clvl()
            dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pP =>phi%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFx(i,j,k,1)=-(pP(i,j,k,1)-pP(i-1,j,k,1))*dxi
               end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFy(i,j,k,1)=-(pP(i,j,k,1)-pP(i,j-1,k,1))*dyi
               end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFz(i,j,k,1)=-(pP(i,j,k,1)-pP(i,j,k-1,1))*dzi
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      else
         ! MLMG path: use psolver's C/F-consistent internal fluxes
         call this%psolver%get_fluxes(Fx,Fy,Fz)
      end if
      ! Apply to face velocities
      do lvl=0,this%amr%clvl()
         call this%U%mf(lvl)%saxpy(scale,Fx(lvl),1,1,1,0)
         call this%V%mf(lvl)%saxpy(scale,Fy(lvl),1,1,1,0)
         call this%W%mf(lvl)%saxpy(scale,Fz(lvl),1,1,1,0)
      end do
      ! Destroy temps
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Fx(lvl))
         call this%amr%mfab_destroy(Fy(lvl))
         call this%amr%mfab_destroy(Fz(lvl))
      end do
   end subroutine correct_velocity

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Compute momentum advection and viscous terms for all levels
   !> No pressure gradient, user can add it in the main loop
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dmomdt(this,U,V,W,drhoUdt,drhoVdt,drhoWdt)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_destroy,amrex_mfiter,amrex_box
      use amrex_interface,  only: amrmfab_average_down_cell,amrmfab_average_down_edge
      implicit none
      class(amrincomp), intent(inout) :: this
      class(amrdata), intent(inout) :: U,V,W                          !< Velocity state (face-centered)
      class(amrdata), intent(inout) :: drhoUdt,drhoVdt,drhoWdt        !< Output: momentum RHS (face-centered)
      ! Flux MultiFabs (9 total: 3 CC, 2 xy-edge, 2 xz-edge, 2 yz-edge)
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FUx,FUy,FUz
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FVx,FVy,FVz
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: FWx,FWy,FWz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: lvl,i,j,k
      real(WP) :: dxi,dyi,dzi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFUx,pFUy,pFUz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFVx,pFVy,pFVz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFWx,pFWy,pFWz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pdUdt,pdVdt,pdWdt
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc

      ! Compute fluxes on all levels
      do lvl=0,this%amr%clvl()
         ! Get mesh size
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Build momentum flux MultiFabs
         ! FUx, FVy, FWz: cell-centered (diagonal fluxes)
         call this%amr%mfab_build(lvl,FUx(lvl),ncomp=1,nover=1,atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl,FVy(lvl),ncomp=1,nover=1,atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl,FWz(lvl),ncomp=1,nover=1,atface=[.false.,.false.,.false.])
         ! FUy, FVx: xy-edge (cross-fluxes)
         call this%amr%mfab_build(lvl,FUy(lvl),ncomp=1,nover=0,atface=[.true.,.true.,.false.])
         call this%amr%mfab_build(lvl,FVx(lvl),ncomp=1,nover=0,atface=[.true.,.true.,.false.])
         ! FUz, FWx: xz-edge (cross-fluxes)
         call this%amr%mfab_build(lvl,FUz(lvl),ncomp=1,nover=0,atface=[.true.,.false.,.true.])
         call this%amr%mfab_build(lvl,FWx(lvl),ncomp=1,nover=0,atface=[.true.,.false.,.true.])
         ! FVz, FWy: yz-edge (cross-fluxes)
         call this%amr%mfab_build(lvl,FVz(lvl),ncomp=1,nover=0,atface=[.false.,.true.,.true.])
         call this%amr%mfab_build(lvl,FWy(lvl),ncomp=1,nover=0,atface=[.false.,.true.,.true.])
         ! MFIter loop: compute all 9 fluxes
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            ! Cell-centered tile
            bx=mfi%tilebox()
            ! Get pointers to data
            pU=>U%mf(lvl)%dataptr(mfi)
            pV=>V%mf(lvl)%dataptr(mfi)
            pW=>W%mf(lvl)%dataptr(mfi)
            pFUx=>FUx(lvl)%dataptr(mfi)
            pFUy=>FUy(lvl)%dataptr(mfi)
            pFUz=>FUz(lvl)%dataptr(mfi)
            pFVx=>FVx(lvl)%dataptr(mfi)
            pFVy=>FVy(lvl)%dataptr(mfi)
            pFVz=>FVz(lvl)%dataptr(mfi)
            pFWx=>FWx(lvl)%dataptr(mfi)
            pFWy=>FWy(lvl)%dataptr(mfi)
            pFWz=>FWz(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            ! Diagonal fluxes
            do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+1
               pFUx(i,j,k,1)=-0.25_WP*this%rho*sum(pU(i:i+1,j,k,1))**2+2.0_WP*pVisc(i,j,k,1)*(pU(i+1,j,k,1)-pU(i,j,k,1))*dxi
               pFVy(i,j,k,1)=-0.25_WP*this%rho*sum(pV(i,j:j+1,k,1))**2+2.0_WP*pVisc(i,j,k,1)*(pV(i,j+1,k,1)-pV(i,j,k,1))*dyi
               pFWz(i,j,k,1)=-0.25_WP*this%rho*sum(pW(i,j,k:k+1,1))**2+2.0_WP*pVisc(i,j,k,1)*(pW(i,j,k+1,1)-pW(i,j,k,1))*dzi
            end do; end do; end do
            ! xy-edge (FUy, FVx): nodal in x,y; cell in z -> [lo,hi] in z; [lo,hi+1] in x,y
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)+1; do i=bx%lo(1),bx%hi(1)+1
               pFUy(i,j,k,1)=-0.25_WP*this%rho*sum(pV(i-1:i,j,k,1))*sum(pU(i,j-1:j,k,1))+0.25_WP*sum(pVisc(i-1:i,j-1:j,k,1))*((pU(i,j,k,1)-pU(i,j-1,k,1))*dyi+(pV(i,j,k,1)-pV(i-1,j,k,1))*dxi)
               pFVx(i,j,k,1)=pFUy(i,j,k,1)
            end do; end do; end do
            ! yz-edge (FVz, FWy): nodal in y,z; cell in x -> [lo,hi] in x; [lo,hi+1] in y,z
            do k=bx%lo(3),bx%hi(3)+1; do j=bx%lo(2),bx%hi(2)+1; do i=bx%lo(1),bx%hi(1)
               pFVz(i,j,k,1)=-0.25_WP*this%rho*sum(pW(i,j-1:j,k,1))*sum(pV(i,j,k-1:k,1))+0.25_WP*sum(pVisc(i,j-1:j,k-1:k,1))*((pV(i,j,k,1)-pV(i,j,k-1,1))*dzi+(pW(i,j,k,1)-pW(i,j-1,k,1))*dyi)
               pFWy(i,j,k,1)=pFVz(i,j,k,1)
            end do; end do; end do
            ! zx-edge (FWx, FUz): nodal in z,x; cell in y -> [lo,hi] in y; [lo,hi+1] in z,x
            do k=bx%lo(3),bx%hi(3)+1; do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)+1
               pFWx(i,j,k,1)=-0.25_WP*this%rho*sum(pU(i,j,k-1:k,1))*sum(pW(i-1:i,j,k,1))+0.25_WP*sum(pVisc(i-1:i,j,k-1:k,1))*((pW(i,j,k,1)-pW(i-1,j,k,1))*dxi+(pU(i,j,k,1)-pU(i,j,k-1,1))*dzi)
               pFUz(i,j,k,1)=pFWx(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do

      ! Average down fluxes (fine -> coarse) for conservation
      do lvl=this%amr%clvl(),1,-1
         ! Cell-centered fluxes
         call amrmfab_average_down_cell(fmf=FUx(lvl),cmf=FUx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_cell(fmf=FVy(lvl),cmf=FVy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_cell(fmf=FWz(lvl),cmf=FWz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         ! Edge-centered fluxes
         call amrmfab_average_down_edge(fmf=FUy(lvl),cmf=FUy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_edge(fmf=FVx(lvl),cmf=FVx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_edge(fmf=FUz(lvl),cmf=FUz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_edge(fmf=FWx(lvl),cmf=FWx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_edge(fmf=FVz(lvl),cmf=FVz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         call amrmfab_average_down_edge(fmf=FWy(lvl),cmf=FWy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
      end do

      ! Compute divergence to get momentum RHS
      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            ! Get pointers to data
            pFUx=>FUx(lvl)%dataptr(mfi)
            pFUy=>FUy(lvl)%dataptr(mfi)
            pFUz=>FUz(lvl)%dataptr(mfi)
            pFVx=>FVx(lvl)%dataptr(mfi)
            pFVy=>FVy(lvl)%dataptr(mfi)
            pFVz=>FVz(lvl)%dataptr(mfi)
            pFWx=>FWx(lvl)%dataptr(mfi)
            pFWy=>FWy(lvl)%dataptr(mfi)
            pFWz=>FWz(lvl)%dataptr(mfi)
            pdUdt=>drhoUdt%mf(lvl)%dataptr(mfi)
            pdVdt=>drhoVdt%mf(lvl)%dataptr(mfi)
            pdWdt=>drhoWdt%mf(lvl)%dataptr(mfi)
            ! U-momentum RHS at x-faces: -d(FUx)/dx - d(FUy)/dy - d(FUz)/dz
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pdUdt(i,j,k,1)=dxi*(pFUx(i,j,k,1)-pFUx(i-1,j,k,1))+dyi*(pFUy(i,j+1,k,1)-pFUy(i,j,k,1))+dzi*(pFUz(i,j,k+1,1)-pFUz(i,j,k,1))
            end do; end do; end do
            ! V-momentum RHS at y-faces: -d(FVx)/dx - d(FVy)/dy - d(FVz)/dz
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pdVdt(i,j,k,1)=dxi*(pFVx(i+1,j,k,1)-pFVx(i,j,k,1))+dyi*(pFVy(i,j,k,1)-pFVy(i,j-1,k,1))+dzi*(pFVz(i,j,k+1,1)-pFVz(i,j,k,1))
            end do; end do; end do
            ! W-momentum RHS at z-faces: -d(FWx)/dx - d(FWy)/dy - d(FWz)/dz
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pdWdt(i,j,k,1)=dxi*(pFWx(i+1,j,k,1)-pFWx(i,j,k,1))+dyi*(pFWy(i,j+1,k,1)-pFWy(i,j,k,1))+dzi*(pFWz(i,j,k,1)-pFWz(i,j,k-1,1))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do

      ! Cleanup flux MultiFabs
      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(FUx(lvl))
         call amrex_multifab_destroy(FUy(lvl))
         call amrex_multifab_destroy(FUz(lvl))
         call amrex_multifab_destroy(FVx(lvl))
         call amrex_multifab_destroy(FVy(lvl))
         call amrex_multifab_destroy(FVz(lvl))
         call amrex_multifab_destroy(FWx(lvl))
         call amrex_multifab_destroy(FWy(lvl))
         call amrex_multifab_destroy(FWz(lvl))
      end do

   end subroutine get_dmomdt

   !> Add Vreman SGS eddy viscosity to this%visc
   !> Assumes velocity ghosts are filled. User must reset visc
   !> to molecular value before calling this routine.
   subroutine add_vreman(this,dt,Cs)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: visc_t,scratch
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pScratch,pU,pV,pW,pVisc
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,Cmodel,Aij,Bij
      real(WP), dimension(1:3,1:3) :: gradU,betaij
      integer :: lvl,i,j,k,si,sj,sk,n
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]

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
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pVisc_t=>visc_t%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Cell-centered velocity gradient tensor
               gradU(1,1)=        dxi*   (pU(i+1,j,k,1)      -pU(i,j,k,1)        )
               gradU(2,1)=0.25_WP*dyi*sum(pU(i:i+1,j:j+1,k,1)-pU(i:i+1,j-1:j,k,1))
               gradU(3,1)=0.25_WP*dzi*sum(pU(i:i+1,j,k:k+1,1)-pU(i:i+1,j,k-1:k,1))
               gradU(1,2)=0.25_WP*dxi*sum(pV(i:i+1,j:j+1,k,1)-pV(i-1:i,j:j+1,k,1))
               gradU(2,2)=        dyi*   (pV(i,j+1,k,1)      -pV(i,j,k,1)        )
               gradU(3,2)=0.25_WP*dzi*sum(pV(i,j:j+1,k:k+1,1)-pV(i,j:j+1,k-1:k,1))
               gradU(1,3)=0.25_WP*dxi*sum(pW(i:i+1,j,k:k+1,1)-pW(i-1:i,j,k:k+1,1))
               gradU(2,3)=0.25_WP*dyi*sum(pW(i,j:j+1,k:k+1,1)-pW(i,j-1:j,k:k+1,1))
               gradU(3,3)=        dzi*   (pW(i,j,k+1,1)      -pW(i,j,k,1)        )
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
         call this%amr%mfab_build(lvl=lvl,mfab=scratch,ncomp=1,nover=1); call scratch%setval(0.0_WP)
         do n=1,nfilter
            call scratch%setval(0.0_WP)
            call scratch%copy(srcmf=visc_t,srccomp=1,dstcomp=1,nc=1,ng=0)
            call this%amr%mfab_validextrap(lvl=lvl,mfab=scratch)
            call scratch%fill_boundary(this%amr%geom(lvl))
            call this%amr%mfab_foextrap(lvl=lvl,mfab=scratch)
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               pScratch=>scratch%dataptr(mfi)
               pVisc_t=>visc_t%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pVisc_t(i,j,k,1)=0.0_WP
                  do sk=-1,+1; do sj=-1,+1; do si=-1,+1
                     pVisc_t(i,j,k,1)=pVisc_t(i,j,k,1)+filter(si)*filter(sj)*filter(sk)*pScratch(i+si,j+sj,k+sk,1)
                  end do; end do; end do
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
         call amrex_multifab_destroy(scratch)
         ! Fill visc_t ghosts after filtering
         call this%amr%mfab_validextrap(lvl=lvl,mfab=visc_t)
         call visc_t%fill_boundary(this%amr%geom(lvl))
         call this%amr%mfab_foextrap(lvl=lvl,mfab=visc_t)
         ! Phase 3: Convert to dynamic viscosity and add to this%visc
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)*this%rho
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Destroy temp multifab
         call amrex_multifab_destroy(visc_t)
      end do
      ! Average down and fill
      call this%visc%average_down(); call this%visc%fill(time=0.0_WP)
   end subroutine add_vreman

   !> Compute CFL numbers (convective and viscous)
   subroutine get_cfl(this,dt,cfl,cflc)
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP), intent(out), optional :: cflc
      integer :: lvl
      real(WP) :: Umax_lvl,Vmax_lvl,Wmax_lvl
      ! Reset CFLs
      this%CFLc_x=0.0_WP; this%CFLc_y=0.0_WP; this%CFLc_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         Umax_lvl=this%U%norm0(lvl=lvl)
         Vmax_lvl=this%V%norm0(lvl=lvl)
         Wmax_lvl=this%W%norm0(lvl=lvl)
         ! Convective CFL
         if (this%amr%nx.gt.1) this%CFLc_x=max(this%CFLc_x,dt*Umax_lvl/this%amr%dx(lvl))
         if (this%amr%ny.gt.1) this%CFLc_y=max(this%CFLc_y,dt*Vmax_lvl/this%amr%dy(lvl))
         if (this%amr%nz.gt.1) this%CFLc_z=max(this%CFLc_z,dt*Wmax_lvl/this%amr%dz(lvl))
         ! Viscous CFL (explicit stability: dt < dx^2 / (4*nu))
         if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*this%visc%norm0(lvl=lvl)*dt/(this%rho*this%amr%dx(lvl)**2))
         if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*this%visc%norm0(lvl=lvl)*dt/(this%rho*this%amr%dy(lvl)**2))
         if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*this%visc%norm0(lvl=lvl)*dt/(this%rho*this%amr%dz(lvl)**2))
      end do
      ! Compute max overall CFL
      this%CFL=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
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
      class(amrincomp), intent(inout) :: this
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
      if (present(VF)) then; if (VF%ng.lt.1) call die('[amrincomp correct_outflow] VF must have at least 1 ghost cell'); end if
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
      class(amrincomp), intent(inout) :: this
      integer :: lvl

      ! First compute divergence (this updates divmax)
      call this%get_div()

      ! Initialize min/max values
      this%Umax=-huge(1.0_WP)
      this%Vmax=-huge(1.0_WP)
      this%Wmax=-huge(1.0_WP)
      this%Pmax=-huge(1.0_WP)

      ! Loop over all levels for min/max
      do lvl=0,this%amr%clvl()
         this%Umax=max(this%Umax,this%U%norm0(lvl=lvl))
         this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl))
         this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl))
         this%Pmax=max(this%Pmax,this%P%norm0(lvl=lvl))
      end do

      ! Momentum integrals (rho * U * dV, summed over cells at level 0)
      this%rhoUint=this%rho*this%U%get_sum(lvl=0)*this%amr%cell_vol(0)
      this%rhoVint=this%rho*this%V%get_sum(lvl=0)*this%amr%cell_vol(0)
      this%rhoWint=this%rho*this%W%get_sum(lvl=0)*this%amr%cell_vol(0)

      ! Kinetic energy integral: 0.5 * rho * (Uc^2 + Vc^2 + Wc^2) * dV
      get_kinetic_energy: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         integer :: i,j,k,ierr
         real(WP) :: Uc,Vc,Wc
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         ! Uses composite integration with fine masking to avoid double-counting
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
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over tile
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then
                     if (pMask(i,j,k,1).eq.0) cycle
                  end if
                  ! Interpolate face velocities to cell center
                  Uc=0.5_WP*(pU(i,j,k,1)+pU(i+1,j,k,1))
                  Vc=0.5_WP*(pV(i,j,k,1)+pV(i,j+1,k,1))
                  Wc=0.5_WP*(pW(i,j,k,1)+pW(i,j,k+1,1))
                  ! Accumulate kinetic energy
                  this%rhoKint=this%rhoKint+0.5_WP*this%rho*(Uc**2+Vc**2+Wc**2)*this%amr%cell_vol(lvl)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_kinetic_energy

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrincomp_print(this)
      use messager, only: log
      use string, only: str_long
      implicit none
      class(amrincomp), intent(in) :: this
      character(len=str_long) :: message
      call log("Incompressible solver: "//trim(this%name))
      write(message,'("  rho = ",ES12.5)') this%rho
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrincomp_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%U,'U')
      call io%add_data(this%V,'V')
      call io%add_data(this%W,'W')
      call io%add_data(this%P,'P')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrincomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      call io%read_data(dirname,this%U,'U')
      call io%read_data(dirname,this%V,'V')
      call io%read_data(dirname,this%W,'W')
      call io%read_data(dirname,this%P,'P')
      ! Fill ghost cells (VisMF reads valid data only)
      call this%fill_velocity(time=time)
      call this%P%fill(time=time)
   end subroutine restore_checkpoint

end module amrincomp_class
