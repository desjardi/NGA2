!> AMR flow solver class: base class for all single-phase flow solvers
!> Provides face-centered velocity storage, divergence calculation,
!> convective CFL calculation, and outflow correction
!> Provides storage for conserved quantity Q(nQ)
module amrflow_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrdata_class,    only: amrdata,interp_face_div,interp_ublin
   use amrsolver_class,  only: amrsolver
   implicit none
   private

   ! Expose type
   public :: amrflow

   !> AMR flow solver type
   type, extends(amrsolver) :: amrflow

      ! Overlap size
      integer :: nover=1

      ! Face velocities (current and old)
      type(amrdata) :: U,Uold
      type(amrdata) :: V,Vold
      type(amrdata) :: W,Wold

      ! Velocity interpolation method
      integer :: interp_vel=interp_face_div

      ! Divergence
      type(amrdata) :: div

      ! Conserved quantity Q with nQ components (current and old)
      type(amrdata) :: Q,Qold
      integer :: nQ=0

      ! Q interpolation method
      integer :: interp_Q=interp_ublin

      ! Monitoring quantities
      real(WP) :: Umax=0.0_WP,Vmax=0.0_WP,Wmax=0.0_WP,divmax=0.0_WP
      real(WP), dimension(:), allocatable :: Qmax,Qmin,Qint

      ! Convective CFLs
      real(WP) :: CFLc_x=0.0_WP,CFLc_y=0.0_WP,CFLc_z=0.0_WP,CFLc=0.0_WP
      
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
      ! Face velocity fills
      procedure :: fill_velocity_lvl         !< Fill face velocity ghosts at single level
      procedure :: fill_velocity             !< Fill face velocity ghosts on all levels
      procedure :: fill_velocity_from_coarse !< Fill face velocity from coarse
      procedure :: sync_velocity_lvl         !< Sync face velocity ghosts at single level
      procedure :: sync_velocity             !< Sync face velocity ghosts on all levels
      procedure :: fill_velocity_mfab        !< Fill dest MultiFabs for regridding
      procedure :: average_down_velocity     !< Average down face velocity for C/F consistency
      procedure :: average_down_velocity_to  !< Average down face velocity for single level
      ! Boundary condition hooks (children should override)
      procedure :: apply_velbc               !< Velocity BC hook
      procedure :: apply_Qbc                 !< Q BC hook
      ! Utilities
      procedure :: get_div                   !< Compute divergence (assumes velocity ghosts filled)
      procedure :: get_cflc                  !< Compute convective CFL
      procedure :: correct_outflow           !< Correct outflow for global mass conservation
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrflow_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrflow

contains

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the flow solver
   subroutine initialize(this,amr,name)
      use amrgrid_class, only: amrgrid
      implicit none
      class(amrflow), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name
      ! Set name
      if (present(name)) then
         this%name=trim(name)
      else
         this%name='UNNAMED_AMRFLOW'
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
      ! Set velocity fillbc callbacks to shared internal handler
      this%U%fillbc=>velocity_fillbc
      this%V%fillbc=>velocity_fillbc
      this%W%fillbc=>velocity_fillbc
      ! Initialize divergence
      call this%div%initialize(amr,name='div',ncomp=1,ng=0); this%div%parent=>this
      ! Initialize Q
      if (this%nQ.gt.0) then
         ! Initialize conserved variables
         call this%Q%initialize   (amr,name='Q'   ,ncomp=this%nQ,ng=this%nover,interp=this%interp_Q); this%Q%parent   =>this
         call this%Qold%initialize(amr,name='Qold',ncomp=this%nQ,ng=this%nover,interp=this%interp_Q); this%Qold%parent=>this
         ! Set Q fillbc callback to internal handler
         this%Q%fillbc=>Q_fillbc
         ! Allocate min/max/integral arrays
         allocate(this%Qmax(this%nQ)); this%Qmax=0.0_WP
         allocate(this%Qmin(this%nQ)); this%Qmin=0.0_WP
         allocate(this%Qint(this%nQ)); this%Qint=0.0_WP
      end if
      ! Print solver info
      call this%print()
   end subroutine initialize

   !> Finalize the flow solver
   subroutine finalize(this)
      implicit none
      class(amrflow), intent(inout) :: this
      call this%U%finalize()
      call this%V%finalize()
      call this%W%finalize()
      call this%Uold%finalize()
      call this%Vold%finalize()
      call this%Wold%finalize()
      call this%div%finalize()
      if (this%nQ.gt.0) then
         call this%Q%finalize()
         call this%Qold%finalize()
         deallocate(this%Qmax,this%Qmin,this%Qint)
      end if
      nullify(this%amr)
   end subroutine finalize

   ! ============================================================================
   ! LIFECYCLE CALLBACKS
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrflow), intent(inout) :: this
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
      call this%div%reset_level(lvl,ba,dm)
      ! Set to zero
      call this%U%setval(val=0.0_WP,lvl=lvl)
      call this%V%setval(val=0.0_WP,lvl=lvl)
      call this%W%setval(val=0.0_WP,lvl=lvl)
      call this%Uold%setval(val=0.0_WP,lvl=lvl)
      call this%Vold%setval(val=0.0_WP,lvl=lvl)
      call this%Wold%setval(val=0.0_WP,lvl=lvl)
      call this%div%setval(val=0.0_WP,lvl=lvl)
      ! Conserved quantity Q
      if (this%nQ.gt.0) then
         call this%Q%reset_level(lvl,ba,dm)
         call this%Qold%reset_level(lvl,ba,dm)
         call this%Q%setval(val=0.0_WP,lvl=lvl)
         call this%Qold%setval(val=0.0_WP,lvl=lvl)
      end if
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using divergence-free interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Face velocity: allocate then fill with divergence-free interpolation
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%fill_velocity_from_coarse(lvl,time)
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      ! Divergence just needs to be reset
      call this%div%reset_level(lvl,ba,dm)
      ! Conserved quantity Q
      if (this%nQ.gt.0) then
         call this%Q%on_coarse(lvl,time,ba,dm)
         call this%Qold%reset_level(lvl,ba,dm)
      end if
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Face velocity: fill with div-free interpolation
      face_vel_remake: block
         use amrex_amr_module, only: amrex_multifab_build,amrex_multifab
         type(amrex_multifab) :: Utmp,Vtmp,Wtmp
         ! Build temp MultiFabs with new layout
         call amrex_multifab_build(Utmp,ba,dm,1,this%U%ng,this%U%nodal)
         call amrex_multifab_build(Vtmp,ba,dm,1,this%V%ng,this%V%nodal)
         call amrex_multifab_build(Wtmp,ba,dm,1,this%W%ng,this%W%nodal)
         ! Fill temps from old data via coupled FillPatch
         call this%fill_velocity_mfab(Utmp,Vtmp,Wtmp,lvl,time,ng=0)
         ! Transfer ownership via pointer swap
         call this%U%mf(lvl)%move(Utmp)
         call this%V%mf(lvl)%move(Vtmp)
         call this%W%mf(lvl)%move(Wtmp)
      end block face_vel_remake
      ! Reset old velocities
      call this%Uold%reset_level(lvl,ba,dm)
      call this%Vold%reset_level(lvl,ba,dm)
      call this%Wold%reset_level(lvl,ba,dm)
      ! Divergence just needs to be reset
      call this%div%reset_level(lvl,ba,dm)
      ! Conserved quantity Q
      if (this%nQ.gt.0) then
         call this%Q%on_remake(lvl,time,ba,dm)
         call this%Qold%reset_level(lvl,ba,dm)
      end if
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Clear velocity
      call this%U%clear_level(lvl)
      call this%V%clear_level(lvl)
      call this%W%clear_level(lvl)
      call this%Uold%clear_level(lvl)
      call this%Vold%clear_level(lvl)
      call this%Wold%clear_level(lvl)
      ! Clear divergence
      call this%div%clear_level(lvl)
      ! Conserved quantity Q
      if (this%nQ.gt.0) then
         call this%Q%clear_level(lvl)
         call this%Qold%clear_level(lvl)
      end if
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Average down face velocities and fill ghosts
      call this%average_down_velocity(lbase)
      call this%fill_velocity(time,lbase)
      ! Average down Q and fill ghosts
      if (this%nQ.gt.0) then
         call this%Q%average_down(lbase)
         call this%Q%fill(time,lbase)
      end if
   end subroutine post_regrid

   ! ============================================================================
   ! Staggered velocity fills
   ! ============================================================================

   !> Average down MAC velocity for a single level (lvl+1 -> lvl)
   !> Uses amrdata infrastructure which handles face-centered averaging correctly
   subroutine average_down_velocity_to(this,lvl)
      implicit none
      class(amrflow), intent(inout) :: this
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
      class(amrflow), intent(inout) :: this
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
      class(amrflow), target, intent(inout) :: this
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
      ! Fill face velocities
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
      class(amrflow), intent(inout) :: this
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
      class(amrflow), target, intent(inout) :: this
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

   !> Sync velocity ghost cells at a single level (no C/F interpolation)
   subroutine sync_velocity_lvl(this,lvl)
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%U%sync_lvl(lvl)
      call this%V%sync_lvl(lvl)
      call this%W%sync_lvl(lvl)
   end subroutine sync_velocity_lvl

   !> Sync velocity ghost cells on all levels
   subroutine sync_velocity(this,lbase)
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in), optional :: lbase
      integer :: lvl,lb
      lb=0; if (present(lbase)) lb=lbase
      do lvl=lb,this%amr%clvl()
         call this%sync_velocity_lvl(lvl)
      end do
   end subroutine sync_velocity

   !> Fill destination MultiFabs with velocity using divergence-free interpolation
   !> Used during regridding (on_remake) to fill new layout MultiFabs
   subroutine fill_velocity_mfab(this,Udest,Vdest,Wdest,lvl,time,ng)
      use iso_c_binding, only: c_loc,c_funloc,c_funptr,c_ptr
      use amrex_interface, only: amrmfab_fillpatch_single,amrmfab_fillpatch_two_faces
      use amrex_amr_module, only: amrex_multifab
      use amrdata_class, only: amrdata_fillbc
      implicit none
      class(amrflow), target, intent(inout) :: this
      type(amrex_multifab), intent(inout) :: Udest,Vdest,Wdest
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in), optional :: ng
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
      ! Fill face velocities
      if (lvl .eq. 0) then
         ! Level 0: single-level fill (just physical BCs)
         call amrmfab_fillpatch_single(Udest,t_old,this%U%mf(0),t_new,this%U%mf(0),this%amr%geom(0),ctx_u,bc_dispatch,time,1,1,1,nghost=ng)
         call amrmfab_fillpatch_single(Vdest,t_old,this%V%mf(0),t_new,this%V%mf(0),this%amr%geom(0),ctx_v,bc_dispatch,time,1,1,1,nghost=ng)
         call amrmfab_fillpatch_single(Wdest,t_old,this%W%mf(0),t_new,this%W%mf(0),this%amr%geom(0),ctx_w,bc_dispatch,time,1,1,1,nghost=ng)
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
         &   1,1,1,rr,this%interp_vel,lo_bc,hi_bc,nghost=ng)
      end if
      ! Reconcile shared face values at box boundaries
      call Udest%override_sync(this%amr%geom(lvl))
      call Vdest%override_sync(this%amr%geom(lvl))
      call Wdest%override_sync(this%amr%geom(lvl))
   end subroutine fill_velocity_mfab

   ! ============================================================================
   ! Boundary conditions
   ! ============================================================================

   !> Internal velocity boundary condition callback (shared by U,V,W)
   !> Handles staggering-aware BC fills for face-centered velocity data
   !> - ext_dir: calls user_bc callback for user-controlled values
   !> - foextrap: copies from interior (Neumann, zero gradient)
   !> - reflect_even/odd: symmetry/anti-symmetry
   subroutine velocity_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_geometry,&
      &                           amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_even,amrex_bc_reflect_odd
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp,ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      class(amrflow), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: dlo(3),dhi(3),flo(3),fhi(3)
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: lvl
      character(len=1) :: comp

      ! Get point to solver
      select type (s=>this%parent)
       class is (amrflow)
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
         integer :: i,j,k,dir,fill_edge,src_from,toff
         integer, dimension(3) :: slo,shi
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
            ! User-controlled: fill via overridable hook
            call solver%apply_velbc(lvl=lvl,time=time,face=face,bx=amrex_box(slo,shi),comp=comp,p=p)
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

   !> Internal conserved quantity boundary condition callback
   !> Calls default_fillbc first, then user_bc for ext_dir faces
   subroutine Q_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,&
      &                           amrex_geometry,amrex_multifab,amrex_bc_ext_dir
      use amrdata_class, only: default_fillbc
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp,ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      class(amrflow), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer, dimension(3) :: dlo,dhi
      integer :: lvl

      ! First apply default BC handling (foextrap,hoextrap,reflect,etc.)
      call default_fillbc(this,mf,scomp,ncomp,time,geom)

      ! Access parent solver
      select type (s=>this%parent)
       class is (amrflow)
         solver=>s
      end select



      ! Get domain bounds and level
      dlo=geom%domain%lo
      dhi=geom%domain%hi
      lvl=this%fill_lvl_cache

      ! Loop over FABs and apply user_bc for ext_dir faces
      call amrex_mfiter_build(mfi,mf,tiling=.false.)
      do while (mfi%next())
         p=>mf%dataptr(mfi)
         ilo=lbound(p,1); ihi=ubound(p,1)
         jlo=lbound(p,2); jhi=ubound(p,2)
         klo=lbound(p,3); khi=ubound(p,3)
         ! X-LOW (face=1)
         if (any(this%lo_bc(1,:).eq.amrex_bc_ext_dir).and.ilo.lt.dlo(1)) call solver%apply_Qbc(lvl=lvl,time=time,face=1,bx=amrex_box([ilo,jlo,klo],[dlo(1)-1,jhi,khi]),p=p)
         ! X-HIGH (face=2)
         if (any(this%hi_bc(1,:).eq.amrex_bc_ext_dir).and.ihi.gt.dhi(1)) call solver%apply_Qbc(lvl=lvl,time=time,face=2,bx=amrex_box([dhi(1)+1,jlo,klo],[ihi,jhi,khi]),p=p)
         ! Y-LOW (face=3)
         if (any(this%lo_bc(2,:).eq.amrex_bc_ext_dir).and.jlo.lt.dlo(2)) call solver%apply_Qbc(lvl=lvl,time=time,face=3,bx=amrex_box([ilo,jlo,klo],[ihi,dlo(2)-1,khi]),p=p)
         ! Y-HIGH (face=4)
         if (any(this%hi_bc(2,:).eq.amrex_bc_ext_dir).and.jhi.gt.dhi(2)) call solver%apply_Qbc(lvl=lvl,time=time,face=4,bx=amrex_box([ilo,dhi(2)+1,klo],[ihi,jhi,khi]),p=p)
         ! Z-LOW (face=5)
         if (any(this%lo_bc(3,:).eq.amrex_bc_ext_dir).and.klo.lt.dlo(3)) call solver%apply_Qbc(lvl=lvl,time=time,face=5,bx=amrex_box([ilo,jlo,klo],[ihi,jhi,dlo(3)-1]),p=p)
         ! Z-HIGH (face=6)
         if (any(this%hi_bc(3,:).eq.amrex_bc_ext_dir).and.khi.gt.dhi(3)) call solver%apply_Qbc(lvl=lvl,time=time,face=6,bx=amrex_box([ilo,jlo,dhi(3)+1],[ihi,jhi,khi]),p=p)
      end do
      call amrex_mfiter_destroy(mfi)

   end subroutine Q_fillbc

   !> Default velocity BC hook: no-op, children should override to fill ext_dir faces
   subroutine apply_velbc(this,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
   end subroutine apply_velbc

   !> Default Q BC hook: no-op, children should override to fill ext_dir faces
   subroutine apply_Qbc(this,lvl,time,face,bx,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrflow), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
   end subroutine apply_Qbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Compute divergence of velocity into internal div field, update divmax
   !> Uses composite fine masking so covered coarse cells don't pollute divmax
   subroutine get_div(this)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmfab_compute_divergence,amrmask_make_fine
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrflow), intent(inout) :: this
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

   !> Compute convective CFL numbers
   subroutine get_cflc(this,dt)
      implicit none
      class(amrflow), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl
      ! Reset CFLs
      this%CFLc_x=0.0_WP; this%CFLc_y=0.0_WP; this%CFLc_z=0.0_WP
      ! Compute convective CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         if (this%amr%nx.gt.1) this%CFLc_x=max(this%CFLc_x,dt*this%U%norm0(lvl=lvl)/this%amr%dx(lvl))
         if (this%amr%ny.gt.1) this%CFLc_y=max(this%CFLc_y,dt*this%V%norm0(lvl=lvl)/this%amr%dy(lvl))
         if (this%amr%nz.gt.1) this%CFLc_z=max(this%CFLc_z,dt*this%W%norm0(lvl=lvl)/this%amr%dz(lvl))
      end do
      ! Compute max overall CFL
      this%CFLc=max(this%CFLc_x,this%CFLc_y,this%CFLc_z)
   end subroutine get_cflc

   !> Correct outflow velocity to ensure global mass conservation
   !> Scans all 6 domain faces: ext_dir faces contribute fixed flux,
   !> foextrap faces are correctable. Correction is distributed uniformly
   !> over all foextrap faces, optionally weighted by VF (fluid volume fraction).
   !> Uses composite integration with fine_mask to avoid double-counting across AMR levels.
   subroutine correct_outflow(this,VF)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      use messager, only: die
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_bc_foextrap,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmask_make_fine
      implicit none
      class(amrflow), intent(inout) :: this
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
      if (present(VF)) then; if (VF%ng.lt.1) call die('[amrflow correct_outflow] VF must have at least 1 ghost cell'); end if
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

   !> Get solver information: min/max velocity and divergence, Q stats
   subroutine get_info(this)
      implicit none
      class(amrflow), intent(inout) :: this
      integer :: lvl,n
      ! First compute divergence (this updates divmax)
      call this%get_div()
      ! Initialize min/max values
      this%Umax=-huge(1.0_WP)
      this%Vmax=-huge(1.0_WP)
      this%Wmax=-huge(1.0_WP)
      if (this%nQ.gt.0) then
         this%Qmax=-huge(1.0_WP)
         this%Qmin=+huge(1.0_WP)
      end if
      ! Loop over all levels for min/max
      do lvl=0,this%amr%clvl()
         this%Umax=max(this%Umax,this%U%norm0(lvl=lvl))
         this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl))
         this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl))
         ! Extrema of conserved variables
         do n=1,this%nQ
            this%Qmin(n)=min(this%Qmin(n),this%Q%get_min(lvl=lvl,comp=n))
            this%Qmax(n)=max(this%Qmax(n),this%Q%get_max(lvl=lvl,comp=n))
         end do
      end do
      ! Integrate Q at base level
      do n=1,this%nQ
         this%Qint(n)=this%Q%get_sum(lvl=0,comp=n)*this%amr%cell_vol(0)
      end do
   end subroutine get_info

   !> Print solver info to screen
   subroutine amrflow_print(this)
      use messager, only: log
      use string, only: str_long,itoa
      implicit none
      class(amrflow), intent(in) :: this
      character(len=str_long) :: message
      call log("Flow solver: "//trim(this%name))
      call log("  Grid: "//trim(this%amr%name))
      call log("    nQ: "//itoa(this%nQ))
   end subroutine amrflow_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrflow), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%U,'U')
      call io%add_data(this%V,'V')
      call io%add_data(this%W,'W')
      if (this%nQ.gt.0) call io%add_data(this%Q,'Q')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrflow), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      call io%read_data(dirname,this%U,'U')
      call io%read_data(dirname,this%V,'V')
      call io%read_data(dirname,this%W,'W')
      if (this%nQ.gt.0) call io%read_data(dirname,this%Q,'Q')
      ! Fill ghost cells (VisMF reads valid data only)
      call this%fill_velocity(time=time)
      if (this%nQ.gt.0) call this%Q%fill(time=time)
   end subroutine restore_checkpoint

end module amrflow_class
