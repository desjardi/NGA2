!> AMR compressible multiphase solver class
!> Inherits from amrmpflow_class
module amrmpcomp_class
   use iso_c_binding,    only: c_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use amrdata_class,    only: amrdata
   use amrmpflow_class,  only: amrmpflow
   use amrmg_class,      only: amrmg
   use amrvof_class,     only: VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter
   implicit none
   private

   ! Expose type and constants
   public :: amrmpcomp,VFlo,VFhi,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER

   !> AMR compressible multiphase solver type
   type, extends(amrmpflow) :: amrmpcomp

      ! User-configurable callbacks
      procedure(mpcomp_init_iface),    pointer, pass :: user_init   =>null()
      procedure(mpcomp_tagging_iface), pointer, pass :: user_tagging=>null()
      procedure(mpcomp_bc_iface),      pointer, pass :: user_bc     =>null()
      procedure(mpcomp_vofbc_iface),   pointer, pass :: user_vofbc  =>null()

      ! Liquid equation of state function pointers: PL=PL(rho,I), CL=CL(rho,P), TL=TL(rho,P)
      procedure(eos_P_iface), pointer, nopass :: getPL=>null()
      procedure(eos_C_iface), pointer, nopass :: getCL=>null()
      procedure(eos_T_iface), pointer, nopass :: getTL=>null()

      ! Gas equation of state function pointers: PG=PG(rho,I), CG=CG(rho,P), TG=TG(rho,P)
      procedure(eos_P_iface), pointer, nopass :: getPG=>null()
      procedure(eos_C_iface), pointer, nopass :: getCG=>null()
      procedure(eos_T_iface), pointer, nopass :: getTG=>null()

      ! Pointer to subroutine for mixture cell relaxation
      procedure(relax_iface), pointer, nopass :: relax=>null()

      ! Pressure solver for pressure projection
      logical :: use_projection=.false.
      type(amrmg) :: psolver

      ! Phasic primitive variables
      type(amrdata) :: RHOL,RHOG         !< Phasic densities
      type(amrdata) :: IL,IG             !< Phasic internal energies
      type(amrdata) :: PL,PG             !< Phasic pressures
      type(amrdata) :: TL,TG             !< Phasic temperatures

      ! Mixture variables
      type(amrdata) :: UVW               !< Cell-centered mixture velocity
      type(amrdata) :: C                 !< Speed of sound
      type(amrdata) :: visc              !< Dynamic viscosity
      type(amrdata) :: beta              !< Bulk viscosity
      type(amrdata) :: diff              !< Heat diffusivity
      real(WP) :: sigma                  !< Surface tension coefficient

      ! CFL numbers
      real(WP) :: CFLst=0.0_WP                               !< Surface tension
      real(WP) :: CFLp=0.0_WP                                !< Pressure+convection
      real(WP) :: CFLa_x=0.0_WP,CFLa_y=0.0_WP,CFLa_z=0.0_WP  !< Acoustic
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous

      ! Monitoring quantities
      real(WP) :: RHOLmin=0.0_WP,RHOLmax=0.0_WP,RHOGmin=0.0_WP,RHOGmax=0.0_WP
      real(WP) :: ILmin=0.0_WP,ILmax=0.0_WP,IGmin=0.0_WP,IGmax=0.0_WP
      real(WP) :: PLmin=0.0_WP,PLmax=0.0_WP,PGmin=0.0_WP,PGmax=0.0_WP
      real(WP) :: TLmin=0.0_WP,TLmax=0.0_WP,TGmin=0.0_WP,TGmax=0.0_WP
      real(WP) :: Cmin=0.0_WP,Cmax=0.0_WP
      real(WP) :: dPmax=0.0_WP
      real(WP) :: rhoKint=0.0_WP

      ! Minimum density for stability
      real(WP) :: rho_floor=1.0e-10_WP

      ! Load distribution diagnostics
      ! Per-rank timing
      real(WP) :: wt_prim =0.0_WP       !< Get_primitive loops
      real(WP) :: wt_dQdt =0.0_WP       !< Full get_dQdt
      real(WP) :: wt_fv   =0.0_WP       !< FV flux loops
      real(WP) :: wt_div  =0.0_WP       !< Divergence + source loops
      real(WP) :: wt_relax=0.0_WP       !< Apply_relax
      real(WP) :: wt_visc =0.0_WP       !< Viscosity models
      ! Reduced timing
      real(WP) :: wtmax_prim =0.0_WP, wtmin_prim =0.0_WP
      real(WP) :: wtmax_dQdt =0.0_WP, wtmin_dQdt =0.0_WP
      real(WP) :: wtmax_fv   =0.0_WP, wtmin_fv   =0.0_WP
      real(WP) :: wtmax_div  =0.0_WP, wtmin_div  =0.0_WP
      real(WP) :: wtmax_relax=0.0_WP, wtmin_relax=0.0_WP
      real(WP) :: wtmax_visc =0.0_WP, wtmin_visc =0.0_WP

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
      ! BC overrides
      procedure :: apply_velbc=>mpcomp_apply_velbc
      procedure :: apply_Qbc=>mpcomp_apply_Qbc
      procedure :: apply_vofbc=>mpcomp_apply_vofbc
      ! Utilities
      procedure :: get_face_velocity          !< Update face velocities from cell-centered data
      procedure :: add_phasic_pressure        !< Add phasic pressure term to face velocities, cell-centered momentum, and phasic internal energies
      procedure :: prepare_psolver            !< Prepare Helmholtz pressure solver
      procedure :: add_pressure_correction    !< Add mixture pressure correction to face velocities, cell-centered momentum, and phasic internal energies
      procedure :: add_surface_tension        !< Add surface tension increment consistently to face velocities and cell-centered momentum
      procedure, private :: apply_face_fluxes !< Apply pre-built face fluxes to face velocities and cell-centered momentum
      ! Physics methods
      procedure :: get_primitive
      procedure :: get_conserved
      procedure :: get_dQdt
      procedure :: build_plic
      procedure :: clean_Q
      procedure :: apply_relax
      procedure :: add_viscartif
      procedure :: add_vreman
      procedure :: get_cfl
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrmpcomp_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrmpcomp

   !> Abstract interface for user-provided on_init callback
   abstract interface
      subroutine mpcomp_init_iface(solver,lvl,time,ba,dm)
         import :: amrmpcomp,WP,amrex_boxarray,amrex_distromap
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine mpcomp_init_iface
   end interface

   !> Abstract interface for user-provided tagging callback
   abstract interface
      subroutine mpcomp_tagging_iface(solver,lvl,time,tags)
         import :: amrmpcomp,c_ptr,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine mpcomp_tagging_iface
   end interface

   !> Abstract interface for user-provided BC callback
   abstract interface
      subroutine mpcomp_bc_iface(solver,lvl,time,face,bx,comp,p)
         use amrex_amr_module, only: amrex_box
         import :: amrmpcomp,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< Can be 'U','V','W','Q'
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      end subroutine mpcomp_bc_iface
   end interface

   !> Abstract interface for user-provided VOF BC callback
   abstract interface
      subroutine mpcomp_vofbc_iface(solver,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
         import :: amrmpcomp,amrex_box,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl,face
         real(WP), intent(in) :: time
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      end subroutine mpcomp_vofbc_iface
   end interface

   !> Abstract interface for EoS: P=P(rho,I)
   abstract interface
      pure real(WP) function eos_P_iface(rho,I)
         import :: WP
         real(WP), intent(in) :: rho
         real(WP), intent(in) :: I
      end function eos_P_iface
   end interface

   !> Abstract interface for EoS: C=C(rho,P)
   abstract interface
      pure real(WP) function eos_C_iface(rho,P)
         import :: WP
         real(WP), intent(in) :: rho
         real(WP), intent(in) :: P
      end function eos_C_iface
   end interface

   !> Abstract interface for EoS: T=T(rho,P)
   abstract interface
      pure real(WP) function eos_T_iface(rho,P)
         import :: WP
         real(WP), intent(in) :: rho
         real(WP), intent(in) :: P
      end function eos_T_iface
   end interface

   !> Abstract interface for pressure relaxation callback
   abstract interface
      subroutine relax_iface(VF,Q,Pjump)
         import :: WP
         real(WP), intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
         real(WP), intent(in) :: Pjump
      end subroutine relax_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrmpcomp_on_init(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrmpcomp_on_coarse(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrmpcomp_on_remake(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrmpcomp_on_clear(ctx,lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrmpcomp_on_clear

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrmpcomp_postregrid(ctx,lbase,time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrmpcomp_postregrid

   !> Dispatch tagging: calls type-bound method then user callback
   subroutine amrmpcomp_tagging(ctx,lvl,time,tags)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_tagging)) call this%user_tagging(lvl,time,tags)
   end subroutine amrmpcomp_tagging

   !> Dispatch cost: calls type-bound method
   subroutine amrmpcomp_get_cost(ctx,lvl,nboxes,costs,ba)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrmpcomp_get_cost

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the compressible solver
   subroutine initialize(this,amr,name)
      use amrgrid_class,    only: amrgrid
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrmg_class,      only: amrmg_varcoef
      implicit none
      class(amrmpcomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Initialize amrmpflow parent with 7 conserved components and at least 2 ghost cells
      this%nQ=7; this%nover=max(this%nover,2)
      call this%amrmpflow%initialize(amr,name); call this%set_parent()

      ! Initialize mixture velocity
      call this%UVW%initialize(amr,name='UVW',ncomp=3,ng=this%nover); this%UVW%parent=>this

      ! Initialize phasic primitive variables
      call this%RHOL%initialize(amr,name='RHOL',ncomp=1,ng=this%nover); this%RHOL%parent=>this
      call this%RHOG%initialize(amr,name='RHOG',ncomp=1,ng=this%nover); this%RHOG%parent=>this
      call this%IL%initialize  (amr,name='IL'  ,ncomp=1,ng=this%nover); this%IL%parent  =>this
      call this%IG%initialize  (amr,name='IG'  ,ncomp=1,ng=this%nover); this%IG%parent  =>this
      call this%PL%initialize  (amr,name='PL'  ,ncomp=1,ng=this%nover); this%PL%parent  =>this
      call this%PG%initialize  (amr,name='PG'  ,ncomp=1,ng=this%nover); this%PG%parent  =>this
      call this%TL%initialize  (amr,name='TL'  ,ncomp=1,ng=this%nover); this%TL%parent  =>this
      call this%TG%initialize  (amr,name='TG'  ,ncomp=1,ng=this%nover); this%TG%parent  =>this

      ! Initialize mixture properties
      call this%C%initialize(amr,name='C',ncomp=1,ng=this%nover); this%C%parent=>this

      ! Initialize physical properties (Neumann BCs on those)
      call this%visc%initialize(amr,name='visc',ncomp=1,ng=this%nover); this%visc%parent=>this
      call this%beta%initialize(amr,name='beta',ncomp=1,ng=this%nover); this%beta%parent=>this
      call this%diff%initialize(amr,name='diff',ncomp=1,ng=this%nover); this%diff%parent=>this
      if (.not.amr%xper) then
         this%visc%lo_bc(1,1)=amrex_bc_foextrap; this%visc%hi_bc(1,1)=amrex_bc_foextrap
         this%beta%lo_bc(1,1)=amrex_bc_foextrap; this%beta%hi_bc(1,1)=amrex_bc_foextrap
         this%diff%lo_bc(1,1)=amrex_bc_foextrap; this%diff%hi_bc(1,1)=amrex_bc_foextrap
      end if
      if (.not.amr%yper) then
         this%visc%lo_bc(2,1)=amrex_bc_foextrap; this%visc%hi_bc(2,1)=amrex_bc_foextrap
         this%beta%lo_bc(2,1)=amrex_bc_foextrap; this%beta%hi_bc(2,1)=amrex_bc_foextrap
         this%diff%lo_bc(2,1)=amrex_bc_foextrap; this%diff%hi_bc(2,1)=amrex_bc_foextrap
      end if
      if (.not.amr%zper) then
         this%visc%lo_bc(3,1)=amrex_bc_foextrap; this%visc%hi_bc(3,1)=amrex_bc_foextrap
         this%beta%lo_bc(3,1)=amrex_bc_foextrap; this%beta%hi_bc(3,1)=amrex_bc_foextrap
         this%diff%lo_bc(3,1)=amrex_bc_foextrap; this%diff%hi_bc(3,1)=amrex_bc_foextrap
      end if

      ! Initialize pressure solver if requested
      if (this%use_projection) call this%psolver%initialize(amr=amr,type=amrmg_varcoef)

      ! Register callbacks with amrgrid
      select type (this)
       type is (amrmpcomp)
         call this%amr%add_on_init   (amrmpcomp_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrmpcomp_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrmpcomp_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrmpcomp_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrmpcomp_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrmpcomp_postregrid,c_loc(this))
         call this%amr%set_get_cost  (amrmpcomp_get_cost,  c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the compressible multiphase solver
   subroutine finalize(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      ! Pressure solver
      call this%psolver%finalize()
      ! Velocity
      call this%UVW%finalize()
      ! Phasic primitives
      call this%RHOL%finalize(); call this%RHOG%finalize()
      call this%IL%finalize(); call this%IG%finalize()
      call this%PL%finalize(); call this%PG%finalize()
      call this%TL%finalize(); call this%TG%finalize()
      ! Mixture properties
      call this%C%finalize()
      ! Physical properties
      call this%visc%finalize(); call this%beta%finalize(); call this%diff%finalize()
      ! Nullify pointers
      nullify(this%user_init); nullify(this%user_tagging); nullify(this%user_bc); nullify(this%user_vofbc)
      nullify(this%getPL); nullify(this%getCL); nullify(this%getTL)
      nullify(this%getPG); nullify(this%getCG); nullify(this%getTG)
      nullify(this%relax)
      ! Finalize parent
      call this%amrmpflow%finalize()
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_init(lvl,time,ba,dm)
      ! Reset mixture velocity
      call this%UVW%reset_level(lvl,ba,dm)
      ! Reset phasic primitives
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm); call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm); call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm); call this%TG%reset_level(lvl,ba,dm)
      ! Reset mixture properties
      call this%C%reset_level(lvl,ba,dm)
      ! Reset physical properties
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      ! Zero out everything
      call this%UVW%setval(val=0.0_WP,lvl=lvl)
      call this%RHOL%setval(val=0.0_WP,lvl=lvl); call this%RHOG%setval(val=0.0_WP,lvl=lvl)
      call this%IL%setval(val=0.0_WP,lvl=lvl); call this%IG%setval(val=0.0_WP,lvl=lvl)
      call this%PL%setval(val=0.0_WP,lvl=lvl); call this%PG%setval(val=0.0_WP,lvl=lvl)
      call this%TL%setval(val=0.0_WP,lvl=lvl); call this%TG%setval(val=0.0_WP,lvl=lvl)
      ! Reset mixture properties
      call this%C%setval(val=0.0_WP,lvl=lvl)
      ! Reset physical properties
      call this%visc%setval(val=0.0_WP,lvl=lvl)
      call this%beta%setval(val=0.0_WP,lvl=lvl)
      call this%diff%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using conservative interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_coarse(lvl,time,ba,dm)
      ! Auxiliary / derived quantities (just reset, will be recomputed)
      call this%UVW%reset_level(lvl,ba,dm)
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm); call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm); call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm); call this%TG%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using conservative interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_remake(lvl,time,ba,dm)
      ! Auxiliary / derived quantities (just reset, will be recomputed)
      call this%UVW%reset_level(lvl,ba,dm)
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm); call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm); call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm); call this%TG%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_clear(lvl)
      ! Auxiliary / derived quantities
      call this%UVW%clear_level(lvl)
      call this%RHOL%clear_level(lvl); call this%RHOG%clear_level(lvl)
      call this%IL%clear_level(lvl); call this%IG%clear_level(lvl)
      call this%PL%clear_level(lvl); call this%PG%clear_level(lvl)
      call this%TL%clear_level(lvl); call this%TG%clear_level(lvl)
      call this%C%clear_level(lvl)
      call this%visc%clear_level(lvl)
      call this%beta%clear_level(lvl)
      call this%diff%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles face velocities, conserved quantities, and VOF-related quantities
      call this%amrmpflow%post_regrid(lbase,time)
      ! Rebuild primitives
      call this%get_primitive(this%Q)
   end subroutine post_regrid

   ! ============================================================================
   ! BOUNDARY CONDITIONS
   ! ============================================================================

   !> Velocity BC override: forward to user_bc with U/V/W component name
   subroutine mpcomp_apply_velbc(this,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp=comp,p=p)
   end subroutine mpcomp_apply_velbc

   !> Q BC override: forward to user_bc with comp='Q'
   subroutine mpcomp_apply_Qbc(this,lvl,time,face,bx,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp='Q',p=p)
   end subroutine mpcomp_apply_Qbc

   !> VOF BC override: forward to user_vofbc
   subroutine mpcomp_apply_vofbc(this,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      if (associated(this%user_vofbc)) call this%user_vofbc(lvl=lvl,time=time,face=face,bx=bx,pVF=pVF,pCL=pCL,pCG=pCG,pPLIC=pPLIC)
   end subroutine mpcomp_apply_vofbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Update face velocity from VF and Q using sub-cell density-weighting (subVF and primitives must be up-to-date)
   subroutine get_face_velocity(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pUVW,pRHOL,pRHOG,pSubVF
      real(WP) :: rhoLo,rhoHi
      ! Traverse levels
      do lvl=0,this%amr%clvl()
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pU   =>this%U%mf(lvl)%dataptr(mfi)
            pV   =>this%V%mf(lvl)%dataptr(mfi)
            pW   =>this%W%mf(lvl)%dataptr(mfi)
            pUVW =>this%UVW%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            ! Get X-face velocity
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               rhoLo=sum(pQ(i-1,j,k,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i-1,j,k,1)*pSubVF(i-1,j,k,2)+pRHOG(i-1,j,k,1)*(1.0_WP-pSubVF(i-1,j,k,2))
                  rhoHi=pRHOL(i  ,j,k,1)*pSubVF(i  ,j,k,1)+pRHOG(i  ,j,k,1)*(1.0_WP-pSubVF(i  ,j,k,1))
               end if
               pU(i,j,k,1)=(rhoLo*pUVW(i-1,j,k,1)+rhoHi*pUVW(i,j,k,1))/(rhoLo+rhoHi)
            end do; end do; end do
            ! Get Y-face velocity
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               rhoLo=sum(pQ(i,j-1,k,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i,j-1,k,1)*pSubVF(i,j-1,k,4)+pRHOG(i,j-1,k,1)*(1.0_WP-pSubVF(i,j-1,k,4))
                  rhoHi=pRHOL(i,j  ,k,1)*pSubVF(i,j  ,k,3)+pRHOG(i,j  ,k,1)*(1.0_WP-pSubVF(i,j  ,k,3))
               end if
               pV(i,j,k,1)=(rhoLo*pUVW(i,j-1,k,2)+rhoHi*pUVW(i,j,k,2))/(rhoLo+rhoHi)
            end do; end do; end do
            ! Get Z-face velocity
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               rhoLo=sum(pQ(i,j,k-1,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i,j,k-1,1)*pSubVF(i,j,k-1,6)+pRHOG(i,j,k-1,1)*(1.0_WP-pSubVF(i,j,k-1,6))
                  rhoHi=pRHOL(i,j,k  ,1)*pSubVF(i,j,k  ,5)+pRHOG(i,j,k  ,1)*(1.0_WP-pSubVF(i,j,k  ,5))
               end if
               pW(i,j,k,1)=(rhoLo*pUVW(i,j,k-1,3)+rhoHi*pUVW(i,j,k,3))/(rhoLo+rhoHi)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_face_velocity

   !> Add explicit phasic pressure gradient to face velocities, cell-centered momentum, and phasic internal energies
   !> Uses this%PL and this%PG directly
   !> Optional gravity(1:3) adds a constant face acceleration alongside -grad(Pmix)/rho
   subroutine add_phasic_pressure(this,scale,gravity,mask)
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: scale
      real(WP), dimension(3), intent(in), optional :: gravity
      type(amrdata), intent(in), optional :: mask
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pQ,pVF,pSubVF,pUVW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pPL,pPG,pRHOL,pRHOG,pMask
      real(WP) :: dxi,dyi,dzi,rhoLo,rhoHi,div,coeff
      integer :: lvl,i,j,k
      ! Build temp face mfabs to store pressure fluxes
      allocate(Fx(0:this%amr%clvl()),Fy(0:this%amr%clvl()),Fz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Compute -grad(Pmix)/rho at faces using phasic pressures
      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pPL  =>this%PL%mf(lvl)%dataptr(mfi)
            pPG  =>this%PG%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pFx  =>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            ! X-faces
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rhoLo=sum(pQ(i-1,j,k,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i-1,j,k,1)*pSubVF(i-1,j,k,2)+pRHOG(i-1,j,k,1)*(1.0_WP-pSubVF(i-1,j,k,2))
                  rhoHi=pRHOL(i  ,j,k,1)*pSubVF(i  ,j,k,1)+pRHOG(i  ,j,k,1)*(1.0_WP-pSubVF(i  ,j,k,1))
               end if
               pFx(i,j,k,1)=-2.0_WP*((pVF(i  ,j,k,1)*pPL(i  ,j,k,1)+(1.0_WP-pVF(i  ,j,k,1))*pPG(i  ,j,k,1))&
               &                    -(pVF(i-1,j,k,1)*pPL(i-1,j,k,1)+(1.0_WP-pVF(i-1,j,k,1))*pPG(i-1,j,k,1)))*dxi/max(rhoLo+rhoHi,this%rho_floor)
               if (present(gravity)) pFx(i,j,k,1)=pFx(i,j,k,1)+gravity(1)
            end do; end do; end do
            ! Y-faces
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rhoLo=sum(pQ(i,j-1,k,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i,j-1,k,1)*pSubVF(i,j-1,k,4)+pRHOG(i,j-1,k,1)*(1.0_WP-pSubVF(i,j-1,k,4))
                  rhoHi=pRHOL(i,j  ,k,1)*pSubVF(i,j  ,k,3)+pRHOG(i,j  ,k,1)*(1.0_WP-pSubVF(i,j  ,k,3))
               end if
               pFy(i,j,k,1)=-2.0_WP*((pVF(i,j  ,k,1)*pPL(i,j  ,k,1)+(1.0_WP-pVF(i,j  ,k,1))*pPG(i,j  ,k,1))&
               &                    -(pVF(i,j-1,k,1)*pPL(i,j-1,k,1)+(1.0_WP-pVF(i,j-1,k,1))*pPG(i,j-1,k,1)))*dyi/max(rhoLo+rhoHi,this%rho_floor)
               if (present(gravity)) pFy(i,j,k,1)=pFy(i,j,k,1)+gravity(2)
            end do; end do; end do
            ! Z-faces
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rhoLo=sum(pQ(i,j,k-1,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i,j,k-1,1)*pSubVF(i,j,k-1,6)+pRHOG(i,j,k-1,1)*(1.0_WP-pSubVF(i,j,k-1,6))
                  rhoHi=pRHOL(i,j,k  ,1)*pSubVF(i,j,k  ,5)+pRHOG(i,j,k  ,1)*(1.0_WP-pSubVF(i,j,k  ,5))
               end if
               pFz(i,j,k,1)=-2.0_WP*((pVF(i,j,k  ,1)*pPL(i,j,k  ,1)+(1.0_WP-pVF(i,j,k  ,1))*pPG(i,j,k  ,1))&
               &                    -(pVF(i,j,k-1,1)*pPL(i,j,k-1,1)+(1.0_WP-pVF(i,j,k-1,1))*pPG(i,j,k-1,1)))*dzi/max(rhoLo+rhoHi,this%rho_floor)
               if (present(gravity)) pFz(i,j,k,1)=pFz(i,j,k,1)+gravity(3)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Enforce flux consistency between levels
      do lvl=this%amr%clvl(),1,-1
         call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
      end do
      ! Apply face acceleration and momentum update
      call this%apply_face_fluxes(scale,Fx,Fy,Fz)
      ! Add phasic pressure-dilatation: -P*div(U)
      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ =>this%Q%mf(lvl)%dataptr(mfi)
            pU =>this%U%mf(lvl)%dataptr(mfi)
            pV =>this%V%mf(lvl)%dataptr(mfi)
            pW =>this%W%mf(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pPL=>this%PL%mf(lvl)%dataptr(mfi)
            pPG=>this%PG%mf(lvl)%dataptr(mfi)
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            if (present(mask)) pMask=>mask%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               div=dxi*(pU(i+1,j,k,1)-pU(i,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k,1))
               coeff=1.0_WP; if (present(mask)) coeff=pMask(i,j,k,1)
               pQ(i,j,k,3)=pQ(i,j,k,3)-scale*coeff*(       pVF(i,j,k,1))*pPL(i,j,k,1)*div
               pQ(i,j,k,4)=pQ(i,j,k,4)-scale*coeff*(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)*div
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temps
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Fx(lvl))
         call this%amr%mfab_destroy(Fy(lvl))
         call this%amr%mfab_destroy(Fz(lvl))
      end do
      deallocate(Fx,Fy,Fz)
   end subroutine add_phasic_pressure

   !> Prepare variable-coefficient pressure solver using face densities
   subroutine prepare_psolver(this,dt)
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(:), allocatable :: AA,BBx,BBy,BBz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pAA,pBBx,pBBy,pBBz,pQ,pC,pSubVF,pRHOL,pRHOG
      real(WP) :: rhoLo,rhoHi
      ! Allocate temporary face coefficient mfabs
      allocate(AA(0:this%amr%clvl()),BBx(0:this%amr%clvl()),BBy(0:this%amr%clvl()),BBz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,AA (lvl),ncomp=1,nover=0,atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl,BBx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,BBy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,BBz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Fill 1/rho_face at all levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pC   =>this%C%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pAA  =>AA(lvl)%dataptr(mfi)
            pBBx =>BBx(lvl)%dataptr(mfi)
            pBBy =>BBy(lvl)%dataptr(mfi)
            pBBz =>BBz(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            ! Cell-centered
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pAA(i,j,k,1)=1.0_WP/(max(sum(pQ(i,j,k,1:2)),this%rho_floor)*pC(i,j,k,1)**2)
            end do; end do; end do
            ! X-faces
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rhoLo=sum(pQ(i-1,j,k,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i-1,j,k,1)*pSubVF(i-1,j,k,2)+pRHOG(i-1,j,k,1)*(1.0_WP-pSubVF(i-1,j,k,2))
                  rhoHi=pRHOL(i  ,j,k,1)*pSubVF(i  ,j,k,1)+pRHOG(i  ,j,k,1)*(1.0_WP-pSubVF(i  ,j,k,1))
               end if
               pBBx(i,j,k,1)=2.0_WP/max(rhoLo+rhoHi,this%rho_floor)
            end do; end do; end do
            ! Y-faces
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rhoLo=sum(pQ(i,j-1,k,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i,j-1,k,1)*pSubVF(i,j-1,k,4)+pRHOG(i,j-1,k,1)*(1.0_WP-pSubVF(i,j-1,k,4))
                  rhoHi=pRHOL(i,j  ,k,1)*pSubVF(i,j  ,k,3)+pRHOG(i,j  ,k,1)*(1.0_WP-pSubVF(i,j  ,k,3))
               end if
               pBBy(i,j,k,1)=2.0_WP/max(rhoLo+rhoHi,this%rho_floor)
            end do; end do; end do
            ! Z-faces
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rhoLo=sum(pQ(i,j,k-1,1:2)); rhoHi=sum(pQ(i,j,k,1:2))
               if (lvl.eq.this%amr%maxlvl) then
                  rhoLo=pRHOL(i,j,k-1,1)*pSubVF(i,j,k-1,6)+pRHOG(i,j,k-1,1)*(1.0_WP-pSubVF(i,j,k-1,6))
                  rhoHi=pRHOL(i,j,k  ,1)*pSubVF(i,j,k  ,5)+pRHOG(i,j,k  ,1)*(1.0_WP-pSubVF(i,j,k  ,5))
               end if
               pBBz(i,j,k,1)=2.0_WP/max(rhoLo+rhoHi,this%rho_floor)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Enforce coefficient consistency between levels
      do lvl=this%amr%clvl(),1,-1
         call amrmfab_average_down_face(fmf=BBx(lvl),cmf=BBx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=BBy(lvl),cmf=BBy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=BBz(lvl),cmf=BBz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
      end do
      ! Rebuild operator
      this%psolver%alpha=-1.0_WP/dt**2
      this%psolver%beta =-1.0_WP
      call this%psolver%setup(acoef=AA,bcoef_x=BBx,bcoef_y=BBy,bcoef_z=BBz)
      ! Destroy temporary mfabs
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(AA (lvl))
         call this%amr%mfab_destroy(BBx(lvl))
         call this%amr%mfab_destroy(BBy(lvl))
         call this%amr%mfab_destroy(BBz(lvl))
      end do
      deallocate(AA,BBx,BBy,BBz)
   end subroutine prepare_psolver

   !> Add pressure correction term to face velocities and cell-centered momentum and phasic internal energies
   !> Assumes mixture pressure only, uses psolver's internal fluxes
   !> Optional mask argument is used for IB masking
   subroutine add_pressure_correction(this,scale,mask)
      use amrex_amr_module, only: amrex_multifab
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrdata), intent(in), optional :: mask
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pQ,pP,pU,pV,pW,pVF,pPL,pPG,pMask
      real(WP) :: dxi,dyi,dzi,div,coeff,crossterm,Pmix
      integer :: lvl,i,j,k
      ! Build temp face mfabs to store pressure fluxes
      allocate(Fx(0:this%amr%clvl()),Fy(0:this%amr%clvl()),Fz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Get face fluxes from MLMG solver
      call this%psolver%get_fluxes(Fx,Fy,Fz)
      ! Apply face acceleration and momentum update
      call this%apply_face_fluxes(scale,Fx,Fy,Fz)
      ! Add phasic pressure-dilatation from correction
      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ =>this%Q%mf(lvl)%dataptr(mfi)
            pU =>this%U%mf(lvl)%dataptr(mfi)
            pV =>this%V%mf(lvl)%dataptr(mfi)
            pW =>this%W%mf(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pPL=>this%PL%mf(lvl)%dataptr(mfi)
            pPG=>this%PG%mf(lvl)%dataptr(mfi)
            pP =>this%psolver%sol%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            if (present(mask)) pMask=>mask%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               coeff=1.0_WP; if (present(mask)) coeff=pMask(i,j,k,1)
               div=dxi*(pU(i+1,j,k,1)-pU(i,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k,1))
               Pmix=pVF(i,j,k,1)*pPL(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)
               crossterm=-coeff*scale**2*Pmix*(dxi*(pFx(i+1,j,k,1)-pFx(i,j,k,1))+dyi*(pFy(i,j+1,k,1)-pFy(i,j,k,1))+dzi*(pFz(i,j,k+1,1)-pFz(i,j,k,1)))
               pQ(i,j,k,3)=pQ(i,j,k,3)-(       pVF(i,j,k,1))*(coeff*scale*pP(i,j,k,1)*div-crossterm)
               pQ(i,j,k,4)=pQ(i,j,k,4)-(1.0_WP-pVF(i,j,k,1))*(coeff*scale*pP(i,j,k,1)*div-crossterm)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temps
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Fx(lvl))
         call this%amr%mfab_destroy(Fy(lvl))
         call this%amr%mfab_destroy(Fz(lvl))
      end do
      deallocate(Fx,Fy,Fz)
   end subroutine add_pressure_correction

   !> Add CSF surface-tension term to face velocities and cell-centered momentum
   !> Builds face fluxes 1/rho*sigma*kappa*grad(VF) at maxlvl, avg_down, then calls apply_face_fluxes
   subroutine add_surface_tension(this,scale)
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface, only: amrmfab_average_down_face
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrex_multifab), dimension(:), allocatable :: STFx,STFy,STFz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSTFx,pSTFy,pSTFz,pVF,pQ,pSubVF,pCurv,pSD,pRHOL,pRHOG
      real(WP) :: dxi,dyi,dzi,mysurf,mycurv,rhoLo,rhoHi
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
         pQ    =>this%Q%mf(lvl)%dataptr(mfi)
         pSubVF=>this%subVF%dataptr(mfi)
         pCurv =>this%curv%dataptr(mfi)
         pSD   =>this%SD%dataptr(mfi)
         pRHOL =>this%RHOL%mf(lvl)%dataptr(mfi)
         pRHOG =>this%RHOG%mf(lvl)%dataptr(mfi)
         pSTFx =>STFx(lvl)%dataptr(mfi)
         pSTFy =>STFy(lvl)%dataptr(mfi)
         pSTFz =>STFz(lvl)%dataptr(mfi)
         ! X-faces
         bx=mfi%nodaltilebox(1)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            mycurv=0.0_WP
            mysurf=sum(pSD(i-1:i,j,k,1)); if (mysurf.gt.0.0_WP) mycurv=sum(pSD(i-1:i,j,k,1)*pCurv(i-1:i,j,k,1))/mysurf
            rhoLo=pRHOL(i-1,j,k,1)*pSubVF(i-1,j,k,2)+pRHOG(i-1,j,k,1)*(1.0_WP-pSubVF(i-1,j,k,2))
            rhoHi=pRHOL(i  ,j,k,1)*pSubVF(i  ,j,k,1)+pRHOG(i  ,j,k,1)*(1.0_WP-pSubVF(i  ,j,k,1))
            pSTFx(i,j,k,1)=2.0_WP*this%sigma*mycurv*(pVF(i,j,k,1)-pVF(i-1,j,k,1))*dxi/max(rhoLo+rhoHi,this%rho_floor)
         end do; end do; end do
         ! Y-faces
         bx=mfi%nodaltilebox(2)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            mycurv=0.0_WP
            mysurf=sum(pSD(i,j-1:j,k,1)); if (mysurf.gt.0.0_WP) mycurv=sum(pSD(i,j-1:j,k,1)*pCurv(i,j-1:j,k,1))/mysurf
            rhoLo=pRHOL(i,j-1,k,1)*pSubVF(i,j-1,k,4)+pRHOG(i,j-1,k,1)*(1.0_WP-pSubVF(i,j-1,k,4))
            rhoHi=pRHOL(i,j  ,k,1)*pSubVF(i,j  ,k,3)+pRHOG(i,j  ,k,1)*(1.0_WP-pSubVF(i,j  ,k,3))
            pSTFy(i,j,k,1)=2.0_WP*this%sigma*mycurv*(pVF(i,j,k,1)-pVF(i,j-1,k,1))*dyi/max(rhoLo+rhoHi,this%rho_floor)
         end do; end do; end do
         ! Z-faces
         bx=mfi%nodaltilebox(3)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            mycurv=0.0_WP
            mysurf=sum(pSD(i,j,k-1:k,1)); if (mysurf.gt.0.0_WP) mycurv=sum(pSD(i,j,k-1:k,1)*pCurv(i,j,k-1:k,1))/mysurf
            rhoLo=pRHOL(i,j,k-1,1)*pSubVF(i,j,k-1,6)+pRHOG(i,j,k-1,1)*(1.0_WP-pSubVF(i,j,k-1,6))
            rhoHi=pRHOL(i,j,k  ,1)*pSubVF(i,j,k  ,5)+pRHOG(i,j,k  ,1)*(1.0_WP-pSubVF(i,j,k  ,5))
            pSTFz(i,j,k,1)=2.0_WP*this%sigma*mycurv*(pVF(i,j,k,1)-pVF(i,j,k-1,1))*dzi/max(rhoLo+rhoHi,this%rho_floor)
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

   !> Apply face-centered acceleration to face velocities and cell-centered momentum
   subroutine apply_face_fluxes(this,scale,Fx,Fy,Fz)
      use amrex_amr_module, only: amrex_multifab,amrex_bc_reflect_odd
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrex_multifab), intent(in) :: Fx(0:),Fy(0:),Fz(0:)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pQ,pSubVF,pRHOL,pRHOG
      real(WP) :: rho,rhoLo,rhoHi
      integer :: lvl,i,j,k
      do lvl=0,this%amr%clvl()
         ! Face velocities: direct saxpy
         call this%U%mf(lvl)%saxpy(scale,Fx(lvl),1,1,1,0)
         call this%V%mf(lvl)%saxpy(scale,Fy(lvl),1,1,1,0)
         call this%W%mf(lvl)%saxpy(scale,Fz(lvl),1,1,1,0)
         ! Cell-centered momentum: density-weighted average of face accelerations
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rho=sum(pQ(i,j,k,1:2))
               ! X-momentum
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=pRHOL(i,j,k,1)*pSubVF(i,j,k,1)+pRHOG(i,j,k,1)*(1.0_WP-pSubVF(i,j,k,1))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=pRHOL(i,j,k,1)*pSubVF(i,j,k,2)+pRHOG(i,j,k,1)*(1.0_WP-pSubVF(i,j,k,2))
               pQ(i,j,k,5)=pQ(i,j,k,5)+scale*0.5_WP*(rhoLo*pFx(i,j,k,1)+rhoHi*pFx(i+1,j,k,1))
               if (.not.this%amr%xper) then
                  if (i.eq.this%amr%geom(lvl)%domain%lo(1).and.this%U%lo_bc(1,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,5)=pQ(i,j,k,5)+scale*0.5_WP*rhoLo*pFx(i+1,j,k,1)
                  if (i.eq.this%amr%geom(lvl)%domain%hi(1).and.this%U%hi_bc(1,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,5)=pQ(i,j,k,5)+scale*0.5_WP*rhoHi*pFx(i  ,j,k,1)
               end if
               ! Y-momentum
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=pRHOL(i,j,k,1)*pSubVF(i,j,k,3)+pRHOG(i,j,k,1)*(1.0_WP-pSubVF(i,j,k,3))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=pRHOL(i,j,k,1)*pSubVF(i,j,k,4)+pRHOG(i,j,k,1)*(1.0_WP-pSubVF(i,j,k,4))
               pQ(i,j,k,6)=pQ(i,j,k,6)+scale*0.5_WP*(rhoLo*pFy(i,j,k,1)+rhoHi*pFy(i,j+1,k,1))
               if (.not.this%amr%yper) then
                  if (j.eq.this%amr%geom(lvl)%domain%lo(2).and.this%V%lo_bc(2,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,6)=pQ(i,j,k,6)+scale*0.5_WP*rhoLo*pFy(i,j+1,k,1)
                  if (j.eq.this%amr%geom(lvl)%domain%hi(2).and.this%V%hi_bc(2,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,6)=pQ(i,j,k,6)+scale*0.5_WP*rhoHi*pFy(i,j  ,k,1)
               end if
               ! Z-momentum
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=pRHOL(i,j,k,1)*pSubVF(i,j,k,5)+pRHOG(i,j,k,1)*(1.0_WP-pSubVF(i,j,k,5))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=pRHOL(i,j,k,1)*pSubVF(i,j,k,6)+pRHOG(i,j,k,1)*(1.0_WP-pSubVF(i,j,k,6))
               pQ(i,j,k,7)=pQ(i,j,k,7)+scale*0.5_WP*(rhoLo*pFz(i,j,k,1)+rhoHi*pFz(i,j,k+1,1))
               if (.not.this%amr%zper) then
                  if (k.eq.this%amr%geom(lvl)%domain%lo(3).and.this%W%lo_bc(3,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,7)=pQ(i,j,k,7)+scale*0.5_WP*rhoLo*pFz(i,j,k+1,1)
                  if (k.eq.this%amr%geom(lvl)%domain%hi(3).and.this%W%hi_bc(3,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,7)=pQ(i,j,k,7)+scale*0.5_WP*rhoHi*pFz(i,j,k  ,1)
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine apply_face_fluxes

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Calculate primitive variables from conserved variables
   !> Q layout: (1) VF*rhoL, (2) (1-VF)*rhoG, (3) VF*rhoL*IL, (4) (1-VF)*rhoG*IG, (5) rhoU, (6) rhoV, (7) rhoW
   subroutine get_primitive(this,Q)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08, only: MPI_Wtime
      use messager, only: die
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(in) :: Q
      integer :: lvl,i,j,k
      real(WP) :: t0
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pUVW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG,pC
      real(WP) :: irho,CL,CG
      ! Start timer
      t0=MPI_Wtime()
      ! Check passed Q is as expected
      if (Q%ncomp.ne.7) call die('[amrmpcomp get_primitive] Q must have 7 components')
      if (Q%ng.lt.this%nover) call die('[amrmpcomp get_primitive] Q must have at least nover ghost cells')
      ! Check EoS functions are set
      if (.not.associated(this%getPL)) call die('[amrmpcomp get_primitive] getPL not set')
      if (.not.associated(this%getCL)) call die('[amrmpcomp get_primitive] getCL not set')
      if (.not.associated(this%getTL)) call die('[amrmpcomp get_primitive] getTL not set')
      if (.not.associated(this%getPG)) call die('[amrmpcomp get_primitive] getPG not set')
      if (.not.associated(this%getCG)) call die('[amrmpcomp get_primitive] getCG not set')
      if (.not.associated(this%getTG)) call die('[amrmpcomp get_primitive] getTG not set')
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>Q%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pUVW =>this%UVW%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL  =>this%IL%mf(lvl)%dataptr(mfi)
            pIG  =>this%IG%mf(lvl)%dataptr(mfi)
            pPL  =>this%PL%mf(lvl)%dataptr(mfi)
            pPG  =>this%PG%mf(lvl)%dataptr(mfi)
            pTL  =>this%TL%mf(lvl)%dataptr(mfi)
            pTG  =>this%TG%mf(lvl)%dataptr(mfi)
            pC   =>this%C%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute mixture velocity from momentum
               irho=1.0_WP/max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
               pUVW(i,j,k,1)=pQ(i,j,k,5)*irho
               pUVW(i,j,k,2)=pQ(i,j,k,6)*irho
               pUVW(i,j,k,3)=pQ(i,j,k,7)*irho
               ! Get liquid primitive variables
               if (pVF(i,j,k,1).ge.VFlo.and.pQ(i,j,k,1).gt.0.0_WP.and.pQ(i,j,k,3).gt.0.0_WP) then
                  pRHOL(i,j,k,1)=pQ(i,j,k,1)/pVF(i,j,k,1)
                  pIL  (i,j,k,1)=pQ(i,j,k,3)/pQ(i,j,k,1)
                  pPL  (i,j,k,1)=this%getPL(pRHOL(i,j,k,1),pIL(i,j,k,1))
                  pTL  (i,j,k,1)=this%getTL(pRHOL(i,j,k,1),pPL(i,j,k,1))
                  CL            =this%getCL(pRHOL(i,j,k,1),pPL(i,j,k,1))
               else
                  pRHOL(i,j,k,1)=0.0_WP
                  pIL  (i,j,k,1)=0.0_WP
                  pPL  (i,j,k,1)=0.0_WP
                  pTL  (i,j,k,1)=0.0_WP
                  CL            =0.0_WP
               end if
               ! Get gas primitive variables
               if (pVF(i,j,k,1).le.VFhi.and.pQ(i,j,k,2).gt.0.0_WP.and.pQ(i,j,k,4).gt.0.0_WP) then
                  pRHOG(i,j,k,1)=pQ(i,j,k,2)/(1.0_WP-pVF(i,j,k,1))
                  pIG  (i,j,k,1)=pQ(i,j,k,4)/pQ(i,j,k,2)
                  pPG  (i,j,k,1)=this%getPG(pRHOG(i,j,k,1),pIG(i,j,k,1))
                  pTG  (i,j,k,1)=this%getTG(pRHOG(i,j,k,1),pPG(i,j,k,1))
                  CG            =this%getCG(pRHOG(i,j,k,1),pPG(i,j,k,1))
               else
                  pRHOG(i,j,k,1)=0.0_WP
                  pIG  (i,j,k,1)=0.0_WP
                  pPG  (i,j,k,1)=0.0_WP
                  pTG  (i,j,k,1)=0.0_WP
                  CG            =0.0_WP
               end if
               ! Get mixture speed of sound
               pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*CL**2+pQ(i,j,k,2)*CG**2)*irho)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! End timer
      this%wt_prim=this%wt_prim+(MPI_Wtime()-t0)
   end subroutine get_primitive

   !> Calculate conserved variables from primitive variables without pressure
   !> Rebuilds Q from VF, phasic densities, phasic energies, and mixture velocity
   subroutine get_conserved(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pUVW,pRHOL,pRHOG,pIL,pIG
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pUVW =>this%UVW%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL  =>this%IL%mf(lvl)%dataptr(mfi)
            pIG  =>this%IG%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pQ(i,j,k,1)=(       pVF(i,j,k,1))*pRHOL(i,j,k,1)
               pQ(i,j,k,2)=(1.0_WP-pVF(i,j,k,1))*pRHOG(i,j,k,1)
               pQ(i,j,k,3)=pQ(i,j,k,1)*pIL(i,j,k,1)
               pQ(i,j,k,4)=pQ(i,j,k,2)*pIG(i,j,k,1)
               pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pUVW(i,j,k,1)
               pQ(i,j,k,6)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pUVW(i,j,k,2)
               pQ(i,j,k,7)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pUVW(i,j,k,3)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_conserved

   !> Calculate dQdt using current this%Q and primitives
   subroutine get_dQdt(this,dQdt,dt,time)
      use amrex_amr_module, only: amrex_multifab
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: dQdt
      real(WP), intent(in) :: dt,time
      real(WP) :: t0,t1
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz
      type(amrex_multifab) :: Vx,Vy,Vz
      type(amrex_multifab) :: band
      ! Shared variables for internal functions
      real(WP) :: dx,dy,dz,dxi,dyi,dzi                              ! Needed for SL transport
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Velocity used for project
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLICold ! PLICold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQold    ! Qold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold   ! VFold used in tet2flux_plic
      logical :: crossed_plic ! Used in tet2flux/tet2flux_plic
      ! Start full routine timer
      t0=MPI_Wtime()

      ! Build transport band at finest level to localize SL computation
      call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=band,ncomp=1,nover=1)
      call this%build_band(lvl=this%amr%clvl(),VF=this%VFold%mf(this%amr%clvl()),band=band,nband=2)

      ! Allocate all fluxes
      define_fluxes: block
         integer :: lvl
         ! Face-centered conserved variable fluxes (7 components)
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_build(lvl=lvl,mfab=Fx(lvl),ncomp=7,nover=0,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl=lvl,mfab=Fy(lvl),ncomp=7,nover=0,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl=lvl,mfab=Fz(lvl),ncomp=7,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
         end do
         ! Volume moment fluxes at finest level (8 components: Lvol,Gvol,Lbar,Gbar)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vx,ncomp=8,nover=0,atface=[.true. ,.false.,.false.]); call Vx%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vy,ncomp=8,nover=0,atface=[.false.,.true. ,.false.]); call Vy%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%clvl(),mfab=Vz,ncomp=8,nover=0,atface=[.false.,.false.,.true. ]); call Vz%setval(0.0_WP)
      end block define_fluxes

      ! Phase 1a: Semi-Lagrangian fluxes at finest level
      t1=MPI_Wtime()
      semilagrangian_fluxes: block
         use amrvof_geometry, only: tet_sign,tet_map,correct_flux_poly
         integer :: lvl,i,j,k,n,nn
         real(WP), dimension(3,9) :: face
         real(WP), dimension(3,4) :: tet
         integer , dimension(3,4) :: ijk
         integer , dimension(3,9) :: fijk
         real(WP), dimension(:,:,:,:), allocatable :: proj
         real(WP), dimension(8) :: Vflux
         real(WP), dimension(7) :: Qflux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pFx,pFy,pFz
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx,nbx
         real(WP) :: rhoLo,rhoHi
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
            pQold   =>this%Qold%mf(lvl)%dataptr(mfi)
            pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
            pBand   =>band%dataptr(mfi)
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
               pVx(i,j,k,1:8)=0.0_WP
               pFx(i,j,k,1:7)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVx(i,j,k,1:8)=pVx(i,j,k,1:8)+tet_sign(tet)*Vflux
                  pFx(i,j,k,1:7)=pFx(i,j,k,1:7)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFx(i,j,k,1:7)=-pFx(i,j,k,1:7)/(dt*dy*dz)
               ! Switch to dissipation-free momentum flux for BB-pure regions
               if (.not.crossed_plic) then
                  rhoLo=max(pQold(i-1,j,k,1)+pQold(i-1,j,k,2),this%rho_floor)
                  rhoHi=max(pQold(i  ,j,k,1)+pQold(i  ,j,k,2),this%rho_floor)
                  pFx(i,j,k,5)=sum(pFx(i,j,k,1:2))*0.5_WP*(pQold(i-1,j,k,5)/rhoLo+pQold(i,j,k,5)/rhoHi)
                  pFx(i,j,k,6)=sum(pFx(i,j,k,1:2))*0.5_WP*(pQold(i-1,j,k,6)/rhoLo+pQold(i,j,k,6)/rhoHi)
                  pFx(i,j,k,7)=sum(pFx(i,j,k,1:2))*0.5_WP*(pQold(i-1,j,k,7)/rhoLo+pQold(i,j,k,7)/rhoHi)
               end if
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
               pVy(i,j,k,1:8)=0.0_WP
               pFy(i,j,k,1:7)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVy(i,j,k,1:8)=pVy(i,j,k,1:8)+tet_sign(tet)*Vflux
                  pFy(i,j,k,1:7)=pFy(i,j,k,1:7)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFy(i,j,k,1:7)=-pFy(i,j,k,1:7)/(dt*dz*dx)
               ! Switch to dissipation-free momentum flux for BB-pure regions
               if (.not.crossed_plic) then
                  rhoLo=max(pQold(i,j-1,k,1)+pQold(i,j-1,k,2),this%rho_floor)
                  rhoHi=max(pQold(i,j  ,k,1)+pQold(i,j  ,k,2),this%rho_floor)
                  pFy(i,j,k,5)=sum(pFy(i,j,k,1:2))*0.5_WP*(pQold(i,j-1,k,5)/rhoLo+pQold(i,j,k,5)/rhoHi)
                  pFy(i,j,k,6)=sum(pFy(i,j,k,1:2))*0.5_WP*(pQold(i,j-1,k,6)/rhoLo+pQold(i,j,k,6)/rhoHi)
                  pFy(i,j,k,7)=sum(pFy(i,j,k,1:2))*0.5_WP*(pQold(i,j-1,k,7)/rhoLo+pQold(i,j,k,7)/rhoHi)
               end if
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
               pVz(i,j,k,1:8)=0.0_WP
               pFz(i,j,k,1:7)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  pVz(i,j,k,1:8)=pVz(i,j,k,1:8)+tet_sign(tet)*Vflux
                  pFz(i,j,k,1:7)=pFz(i,j,k,1:7)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFz(i,j,k,1:7)=-pFz(i,j,k,1:7)/(dt*dx*dy)
               ! Switch to dissipation-free momentum flux for BB-pure regions
               if (.not.crossed_plic) then
                  rhoLo=max(pQold(i,j,k-1,1)+pQold(i,j,k-1,2),this%rho_floor)
                  rhoHi=max(pQold(i,j,k  ,1)+pQold(i,j,k  ,2),this%rho_floor)
                  pFz(i,j,k,5)=sum(pFz(i,j,k,1:2))*0.5_WP*(pQold(i,j,k-1,5)/rhoLo+pQold(i,j,k,5)/rhoHi)
                  pFz(i,j,k,6)=sum(pFz(i,j,k,1:2))*0.5_WP*(pQold(i,j,k-1,6)/rhoLo+pQold(i,j,k,6)/rhoHi)
                  pFz(i,j,k,7)=sum(pFz(i,j,k,1:2))*0.5_WP*(pQold(i,j,k-1,7)/rhoLo+pQold(i,j,k,7)/rhoHi)
               end if
            end do; end do; end do
            ! Deallocate proj for this tile
            deallocate(proj)
         end do
         call this%amr%mfiter_destroy(mfi)
      end block semilagrangian_fluxes

      ! Lower dissipation by blending SL with centered momentum fluxes
      reduce_dissipation: block
         integer :: lvl,i,j,k
         type(amrex_multifab) :: blend
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBlend,pBand,pFx,pFy,pFz
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx,bx
         real(WP) :: rho_old,drho,rhoLo,rhoHi,coeff
         real(WP), dimension(3) :: FC
         ! Skip if clvl < maxlvl
         if (this%amr%clvl().lt.this%amr%maxlvl) exit reduce_dissipation
         ! Get finest level info
         lvl=this%amr%maxlvl
         dx=this%amr%dx(lvl); dxi=1.0_WP/this%amr%dx(lvl)
         dy=this%amr%dy(lvl); dyi=1.0_WP/this%amr%dy(lvl)
         dz=this%amr%dz(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         ! Pass 1: build antidiffusive blend weight (1=full centered, 0=full SL)
         call this%amr%mfab_build(lvl=lvl,mfab=blend,ncomp=1,nover=1); call blend%setval(1.0_WP)
         call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
         do while (mfi%next())
            pQold =>this%Qold%mf(lvl)%dataptr(mfi)
            pBand =>band%dataptr(mfi)
            pFx   =>Fx(lvl)%dataptr(mfi)
            pFy   =>Fy(lvl)%dataptr(mfi)
            pFz   =>Fz(lvl)%dataptr(mfi)
            pBlend=>blend%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBlend(i,j,k,1)=1.0_WP
               if (pBand(i,j,k,1).le.0.0_WP) cycle
               rho_old=max(pQold(i,j,k,1)+pQold(i,j,k,2),this%rho_floor)
               drho=dt*sum(dxi*(pFx(i+1,j,k,1:2)-pFx(i,j,k,1:2))+dyi*(pFy(i,j+1,k,1:2)-pFy(i,j,k,1:2))+dzi*(pFz(i,j,k+1,1:2)-pFz(i,j,k,1:2)))
               pBlend(i,j,k,1)=min(1.0_WP,1.0_WP+drho/rho_old)**2
               !pBlend(i,j,k,1)=min(1.0_WP,max(rho_old+drho,0.0_WP)/max(-drho,this%rho_floor))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call blend%fill_boundary(this%amr%geom(lvl))
         ! Pass 2: F(5:7)=F_SL+coeff*(F_centered-F_SL)
         call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
         do while (mfi%next())
            pQold =>this%Qold%mf(lvl)%dataptr(mfi)
            pBand =>band%dataptr(mfi)
            pFx   =>Fx(lvl)%dataptr(mfi)
            pFy   =>Fy(lvl)%dataptr(mfi)
            pFz   =>Fz(lvl)%dataptr(mfi)
            pBlend=>blend%dataptr(mfi)
            ! X-fluxes
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               if (maxval(pBand(i-1:i,j,k,1)).eq.0.0_WP) cycle
               rhoLo=max(pQold(i-1,j,k,1)+pQold(i-1,j,k,2),this%rho_floor)
               rhoHi=max(pQold(i  ,j,k,1)+pQold(i  ,j,k,2),this%rho_floor)
               FC=sum(pFx(i,j,k,1:2))*0.5_WP*(pQold(i-1,j,k,5:7)/rhoLo+pQold(i,j,k,5:7)/rhoHi)
               coeff=min(pBlend(i-1,j,k,1),pBlend(i,j,k,1))
               pFx(i,j,k,5:7)=pFx(i,j,k,5:7)+coeff*(FC-pFx(i,j,k,5:7))
            end do; end do; end do
            ! Y-fluxes
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               rhoLo=max(pQold(i,j-1,k,1)+pQold(i,j-1,k,2),this%rho_floor)
               rhoHi=max(pQold(i,j  ,k,1)+pQold(i,j  ,k,2),this%rho_floor)
               FC=sum(pFy(i,j,k,1:2))*0.5_WP*(pQold(i,j-1,k,5:7)/rhoLo+pQold(i,j,k,5:7)/rhoHi)
               coeff=min(pBlend(i,j-1,k,1),pBlend(i,j,k,1))
               pFy(i,j,k,5:7)=pFy(i,j,k,5:7)+coeff*(FC-pFy(i,j,k,5:7))
            end do; end do; end do
            ! Z-fluxes
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               rhoLo=max(pQold(i,j,k-1,1)+pQold(i,j,k-1,2),this%rho_floor)
               rhoHi=max(pQold(i,j,k  ,1)+pQold(i,j,k  ,2),this%rho_floor)
               FC=sum(pFz(i,j,k,1:2))*0.5_WP*(pQold(i,j,k-1,5:7)/rhoLo+pQold(i,j,k,5:7)/rhoHi)
               coeff=min(pBlend(i,j,k-1,1),pBlend(i,j,k,1))
               pFz(i,j,k,5:7)=pFz(i,j,k,5:7)+coeff*(FC-pFz(i,j,k,5:7))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Cleanup
         call this%amr%mfab_destroy(blend)
      end block reduce_dissipation
      this%wt_sl=this%wt_sl+(MPI_Wtime()-t1)
      
      ! Phase 1b: Finite volume fluxes for all levels (Euler fluxes skip band cells at finest level)
      t1=MPI_Wtime()
      finitevolume_fluxes: block
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pFx,pFy,pFz,pBand
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pTL,pTG,pIL,pIG,pVF,pUVW!,pPL,pPG
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc,pBeta,pDiff
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Intentional masking
         real(WP), dimension(-2: 0) :: wenop
         real(WP), dimension(-1:+1) :: wenom
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: w,div,visc_f,beta_f
         real(WP), parameter :: eps=1.0e-15_WP
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
               pQ   =>this%Q%mf(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               !pPL  =>this%PL%mf(lvl)%dataptr(mfi)
               !pPG  =>this%PG%mf(lvl)%dataptr(mfi)
               pTL  =>this%TL%mf(lvl)%dataptr(mfi)
               pTG  =>this%TG%mf(lvl)%dataptr(mfi)
               pIL  =>this%IL%mf(lvl)%dataptr(mfi)
               pIG  =>this%IG%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pDiff=>this%diff%mf(lvl)%dataptr(mfi)
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
                     ! WENO liquid mass and energy fluxes
                     if (any(pVF(i-1:i,j,k,1).ge.VFlo)) then
                        w=weno_weight((abs(pQ(i-1,j,k,1)-pQ(i-2,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i+1,j,k,1)-pQ(i  ,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,1)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wenop*pQ(i-2:i  ,j,k,1)) &
                        &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wenom*pQ(i-1:i+1,j,k,1))
                        w=weno_weight((abs(pIL(i-1,j,k,1)-pIL(i-2,j,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIL(i+1,j,k,1)-pIL(i  ,j,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,3)=0.5_WP*(pFx(i,j,k,1)-abs(pFx(i,j,k,1)))*sum(wenop*pIL(i-2:i  ,j,k,1)) &
                        &           +0.5_WP*(pFx(i,j,k,1)+abs(pFx(i,j,k,1)))*sum(wenom*pIL(i-1:i+1,j,k,1))
                     end if
                     ! WENO gas mass and energy fluxes
                     if (any(pVF(i-1:i,j,k,1).le.VFhi)) then
                        w=weno_weight((abs(pQ(i-1,j,k,2)-pQ(i-2,j,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i-1,j,k,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i+1,j,k,2)-pQ(i  ,j,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i-1,j,k,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,2)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wenop*pQ(i-2:i  ,j,k,2)) &
                        &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wenom*pQ(i-1:i+1,j,k,2))
                        w=weno_weight((abs(pIG(i-1,j,k,1)-pIG(i-2,j,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIG(i+1,j,k,1)-pIG(i  ,j,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,4)=0.5_WP*(pFx(i,j,k,2)-abs(pFx(i,j,k,2)))*sum(wenop*pIG(i-2:i  ,j,k,1)) &
                        &           +0.5_WP*(pFx(i,j,k,2)+abs(pFx(i,j,k,2)))*sum(wenom*pIG(i-1:i+1,j,k,1))
                     end if
                     ! Momentum fluxes
                     pFx(i,j,k,5)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pUVW(i-1:i,j,k,1))
                     pFx(i,j,k,6)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pUVW(i-1:i,j,k,2))
                     pFx(i,j,k,7)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pUVW(i-1:i,j,k,3))
                  end if
                  ! Add pressure stress
                  !pFx(i,j,k,5)=pFx(i,j,k,5)-0.5_WP*sum(pVF(i-1:i,j,k,1)*pPL(i-1:i,j,k,1)+(1.0_WP-pVF(i-1:i,j,k,1))*pPG(i-1:i,j,k,1))
                  ! Velocity gradients at x-face
                  gradU(1,1)=dxi*(pUVW(i,j,k,1)-pUVW(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*(pUVW(i-1,j+1,k,1)-pUVW(i-1,j-1,k,1)+pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*(pUVW(i-1,j,k+1,1)-pUVW(i-1,j,k-1,1)+pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
                  gradU(1,2)=dxi*(pUVW(i,j,k,2)-pUVW(i-1,j,k,2))
                  gradU(2,2)=0.25_WP*dyi*(pUVW(i-1,j+1,k,2)-pUVW(i-1,j-1,k,2)+pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
                  gradU(3,2)=0.25_WP*dzi*(pUVW(i-1,j,k+1,2)-pUVW(i-1,j,k-1,2)+pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
                  gradU(1,3)=dxi*(pUVW(i,j,k,3)-pUVW(i-1,j,k,3))
                  gradU(2,3)=0.25_WP*dyi*(pUVW(i-1,j+1,k,3)-pUVW(i-1,j-1,k,3)+pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
                  gradU(3,3)=0.25_WP*dzi*(pUVW(i-1,j,k+1,3)-pUVW(i-1,j,k-1,3)+pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscosities at x-face
                  visc_f=2.0_WP*product(pVisc(i-1:i,j,k,1))/(sum(pVisc(i-1:i,j,k,1))+tiny(1.0_WP))
                  beta_f=2.0_WP*product(pBeta(i-1:i,j,k,1))/(sum(pBeta(i-1:i,j,k,1))+tiny(1.0_WP))
                  ! Viscous stress at x-face
                  pFx(i,j,k,5)=pFx(i,j,k,5)+visc_f*(gradU(1,1)+gradU(1,1))+(beta_f-2.0_WP/3.0_WP*visc_f)*div
                  pFx(i,j,k,6)=pFx(i,j,k,6)+visc_f*(gradU(2,1)+gradU(1,2))
                  pFx(i,j,k,7)=pFx(i,j,k,7)+visc_f*(gradU(3,1)+gradU(1,3))
                  ! Phasic heat diffusion flux (pure cells only)
                  if (all(pVF(i-1:i,j,k,1).gt.VFhi)) pFx(i,j,k,3)=pFx(i,j,k,3)+0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pTL(i,j,k,1)-pTL(i-1,j,k,1))
                  if (all(pVF(i-1:i,j,k,1).lt.VFlo)) pFx(i,j,k,4)=pFx(i,j,k,4)+0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pTG(i,j,k,1)-pTG(i-1,j,k,1))
               end do; end do; end do
               ! Y-fluxes
               fbx=mfi%nodaltilebox(2)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Check if in band
                  if (lvl.eq.this%amr%clvl()) then; in_band=maxval(pBand(i,j-1:j,k,1)).gt.0.0_WP; else; in_band=.false.; end if
                  ! Outside band, compute finite volume Euler fluxes
                  if (.not.in_band) then
                     ! WENO liquid mass and energy fluxes
                     if (any(pVF(i,j-1:j,k,1).ge.VFlo)) then
                        w=weno_weight((abs(pQ(i,j-1,k,1)-pQ(i,j-2,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j+1,k,1)-pQ(i,j  ,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,1)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wenop*pQ(i,j-2:j  ,k,1)) &
                        &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wenom*pQ(i,j-1:j+1,k,1))
                        w=weno_weight((abs(pIL(i,j-1,k,1)-pIL(i,j-2,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIL(i,j+1,k,1)-pIL(i,j  ,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,3)=0.5_WP*(pFy(i,j,k,1)-abs(pFy(i,j,k,1)))*sum(wenop*pIL(i,j-2:j  ,k,1)) &
                        &           +0.5_WP*(pFy(i,j,k,1)+abs(pFy(i,j,k,1)))*sum(wenom*pIL(i,j-1:j+1,k,1))
                     end if
                     ! WENO gas mass and energy fluxes
                     if (any(pVF(i,j-1:j,k,1).le.VFhi)) then
                        w=weno_weight((abs(pQ(i,j-1,k,2)-pQ(i,j-2,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j-1,k,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j+1,k,2)-pQ(i,j  ,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j-1,k,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,2)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wenop*pQ(i,j-2:j  ,k,2)) &
                        &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wenom*pQ(i,j-1:j+1,k,2))
                        w=weno_weight((abs(pIG(i,j-1,k,1)-pIG(i,j-2,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIG(i,j+1,k,1)-pIG(i,j  ,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,4)=0.5_WP*(pFy(i,j,k,2)-abs(pFy(i,j,k,2)))*sum(wenop*pIG(i,j-2:j  ,k,1)) &
                        &           +0.5_WP*(pFy(i,j,k,2)+abs(pFy(i,j,k,2)))*sum(wenom*pIG(i,j-1:j+1,k,1))
                     end if
                     ! Momentum fluxes
                     pFy(i,j,k,5)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pUVW(i,j-1:j,k,1))
                     pFy(i,j,k,6)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pUVW(i,j-1:j,k,2))
                     pFy(i,j,k,7)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pUVW(i,j-1:j,k,3))
                  end if
                  ! Add pressure stress
                  !pFy(i,j,k,6)=pFy(i,j,k,6)-0.5_WP*sum(pVF(i,j-1:j,k,1)*pPL(i,j-1:j,k,1)+(1.0_WP-pVF(i,j-1:j,k,1))*pPG(i,j-1:j,k,1))
                  ! Velocity gradients at y-face
                  gradU(1,1)=0.25_WP*dxi*(pUVW(i+1,j-1,k,1)-pUVW(i-1,j-1,k,1)+pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
                  gradU(2,1)=dyi*(pUVW(i,j,k,1)-pUVW(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*(pUVW(i,j-1,k+1,1)-pUVW(i,j-1,k-1,1)+pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*(pUVW(i+1,j-1,k,2)-pUVW(i-1,j-1,k,2)+pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
                  gradU(2,2)=dyi*(pUVW(i,j,k,2)-pUVW(i,j-1,k,2))
                  gradU(3,2)=0.25_WP*dzi*(pUVW(i,j-1,k+1,2)-pUVW(i,j-1,k-1,2)+pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
                  gradU(1,3)=0.25_WP*dxi*(pUVW(i+1,j-1,k,3)-pUVW(i-1,j-1,k,3)+pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
                  gradU(2,3)=dyi*(pUVW(i,j,k,3)-pUVW(i,j-1,k,3))
                  gradU(3,3)=0.25_WP*dzi*(pUVW(i,j-1,k+1,3)-pUVW(i,j-1,k-1,3)+pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscosities at y-face
                  visc_f=2.0_WP*product(pVisc(i,j-1:j,k,1))/(sum(pVisc(i,j-1:j,k,1))+tiny(1.0_WP))
                  beta_f=2.0_WP*product(pBeta(i,j-1:j,k,1))/(sum(pBeta(i,j-1:j,k,1))+tiny(1.0_WP))
                  ! Viscous stress at y-face
                  pFy(i,j,k,5)=pFy(i,j,k,5)+visc_f*(gradU(1,2)+gradU(2,1))
                  pFy(i,j,k,6)=pFy(i,j,k,6)+visc_f*(gradU(2,2)+gradU(2,2))+(beta_f-2.0_WP/3.0_WP*visc_f)*div
                  pFy(i,j,k,7)=pFy(i,j,k,7)+visc_f*(gradU(3,2)+gradU(2,3))
                  ! Phasic heat diffusion flux (pure cells only)
                  if (all(pVF(i,j-1:j,k,1).gt.VFhi)) pFy(i,j,k,3)=pFy(i,j,k,3)+0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pTL(i,j,k,1)-pTL(i,j-1,k,1))
                  if (all(pVF(i,j-1:j,k,1).lt.VFlo)) pFy(i,j,k,4)=pFy(i,j,k,4)+0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pTG(i,j,k,1)-pTG(i,j-1,k,1))
               end do; end do; end do
               ! Z-fluxes
               fbx=mfi%nodaltilebox(3)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Check if in band
                  if (lvl.eq.this%amr%clvl()) then; in_band=maxval(pBand(i,j,k-1:k,1)).gt.0.0_WP; else; in_band=.false.; end if
                  ! Outside band, compute finite volume Euler fluxes
                  if (.not.in_band) then
                     ! WENO liquid mass and energy fluxes
                     if (any(pVF(i,j,k-1:k,1).ge.VFlo)) then
                        w=weno_weight((abs(pQ(i,j,k-1,1)-pQ(i,j,k-2,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j,k+1,1)-pQ(i,j,k  ,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,1)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wenop*pQ(i,j,k-2:k  ,1)) &
                        &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wenom*pQ(i,j,k-1:k+1,1))
                        w=weno_weight((abs(pIL(i,j,k-1,1)-pIL(i,j,k-2,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIL(i,j,k+1,1)-pIL(i,j,k  ,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,3)=0.5_WP*(pFz(i,j,k,1)-abs(pFz(i,j,k,1)))*sum(wenop*pIL(i,j,k-2:k  ,1)) &
                        &           +0.5_WP*(pFz(i,j,k,1)+abs(pFz(i,j,k,1)))*sum(wenom*pIL(i,j,k-1:k+1,1))
                     end if
                     ! WENO gas mass and energy fluxes
                     if (any(pVF(i,j,k-1:k,1).le.VFhi)) then
                        w=weno_weight((abs(pQ(i,j,k-1,2)-pQ(i,j,k-2,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j,k-1,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j,k+1,2)-pQ(i,j,k  ,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j,k-1,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,2)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wenop*pQ(i,j,k-2:k  ,2)) &
                        &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wenom*pQ(i,j,k-1:k+1,2))
                        w=weno_weight((abs(pIG(i,j,k-1,1)-pIG(i,j,k-2,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIG(i,j,k+1,1)-pIG(i,j,k  ,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,4)=0.5_WP*(pFz(i,j,k,2)-abs(pFz(i,j,k,2)))*sum(wenop*pIG(i,j,k-2:k  ,1)) &
                        &           +0.5_WP*(pFz(i,j,k,2)+abs(pFz(i,j,k,2)))*sum(wenom*pIG(i,j,k-1:k+1,1))
                     end if
                     ! Momentum fluxes
                     pFz(i,j,k,5)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pUVW(i,j,k-1:k,1))
                     pFz(i,j,k,6)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pUVW(i,j,k-1:k,2))
                     pFz(i,j,k,7)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pUVW(i,j,k-1:k,3))
                  end if
                  ! Add pressure stress
                  !pFz(i,j,k,7)=pFz(i,j,k,7)-0.5_WP*sum(pVF(i,j,k-1:k,1)*pPL(i,j,k-1:k,1)+(1.0_WP-pVF(i,j,k-1:k,1))*pPG(i,j,k-1:k,1))
                  ! Velocity gradients at z-face
                  gradU(1,1)=0.25_WP*dxi*(pUVW(i+1,j,k-1,1)-pUVW(i-1,j,k-1,1)+pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*(pUVW(i,j+1,k-1,1)-pUVW(i,j-1,k-1,1)+pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
                  gradU(3,1)=dzi*(pUVW(i,j,k,1)-pUVW(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*(pUVW(i+1,j,k-1,2)-pUVW(i-1,j,k-1,2)+pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
                  gradU(2,2)=0.25_WP*dyi*(pUVW(i,j+1,k-1,2)-pUVW(i,j-1,k-1,2)+pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
                  gradU(3,2)=dzi*(pUVW(i,j,k,2)-pUVW(i,j,k-1,2))
                  gradU(1,3)=0.25_WP*dxi*(pUVW(i+1,j,k-1,3)-pUVW(i-1,j,k-1,3)+pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
                  gradU(2,3)=0.25_WP*dyi*(pUVW(i,j+1,k-1,3)-pUVW(i,j-1,k-1,3)+pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
                  gradU(3,3)=dzi*(pUVW(i,j,k,3)-pUVW(i,j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Viscosities at z-face
                  visc_f=2.0_WP*product(pVisc(i,j,k-1:k,1))/(sum(pVisc(i,j,k-1:k,1))+tiny(1.0_WP))
                  beta_f=2.0_WP*product(pBeta(i,j,k-1:k,1))/(sum(pBeta(i,j,k-1:k,1))+tiny(1.0_WP))
                  ! Viscous stress at z-face
                  pFz(i,j,k,5)=pFz(i,j,k,5)+visc_f*(gradU(1,3)+gradU(3,1))
                  pFz(i,j,k,6)=pFz(i,j,k,6)+visc_f*(gradU(2,3)+gradU(3,2))
                  pFz(i,j,k,7)=pFz(i,j,k,7)+visc_f*(gradU(3,3)+gradU(3,3))+(beta_f-2.0_WP/3.0_WP*visc_f)*div
                  ! Phasic heat diffusion flux (pure cells only)
                  if (all(pVF(i,j,k-1:k,1).gt.VFhi)) pFz(i,j,k,3)=pFz(i,j,k,3)+0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pTL(i,j,k,1)-pTL(i,j,k-1,1))
                  if (all(pVF(i,j,k-1:k,1).lt.VFlo)) pFz(i,j,k,4)=pFz(i,j,k,4)+0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pTG(i,j,k,1)-pTG(i,j,k-1,1))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block finitevolume_fluxes
      this%wt_fv=this%wt_fv+(MPI_Wtime()-t1)

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
      t1=MPI_Wtime()
      divergence_and_sources: block
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: lvl,i,j,k
         real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs,pFx,pFy,pFz,pBand
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFold  ! Intentional masking
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pVisc,pBeta,pUVW!,pPL,pPG
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVx,pVy,pVz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pCL,pCG,pCLold,pCGold
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: div,vol
         real(WP) :: Lvol_old,Lvol_new,Lvol_flux
         real(WP) :: Gvol_old,Gvol_new,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Lbar_new,Lbar_flux
         real(WP), dimension(3) :: Gbar_old,Gbar_new,Gbar_flux
         do lvl=0,this%amr%clvl()
            ! Grid spacings for this level
            dx=this%amr%dx(lvl); dxi=1.0_WP/dx
            dy=this%amr%dy(lvl); dyi=1.0_WP/dy
            dz=this%amr%dz(lvl); dzi=1.0_WP/dz
            vol=dx*dy*dz
            ! Loop over tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               rhs  =>dQdt%mf(lvl)%dataptr(mfi)
               pFx  =>Fx(lvl)%dataptr(mfi)
               pFy  =>Fy(lvl)%dataptr(mfi)
               pFz  =>Fz(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pQ   =>this%Q%mf(lvl)%dataptr(mfi)
               !pPL  =>this%PL%mf(lvl)%dataptr(mfi)
               !pPG  =>this%PG%mf(lvl)%dataptr(mfi)
               pUVW =>this%UVW%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               ! Extra pointers at finest level
               if (lvl.eq.this%amr%maxlvl) then
                  pBand =>band%dataptr(mfi)
                  pVx   =>Vx%dataptr(mfi)
                  pVy   =>Vy%dataptr(mfi)
                  pVz   =>Vz%dataptr(mfi)
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
                        pVF(i,j,k,1)=Lvol_new/(Lvol_new+Gvol_new)
                        pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                        pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                        ! Clip and update barycenters
                        if (pVF(i,j,k,1).lt.VFlo) then
                           pVF(i,j,k,1)=0.0_WP
                        else if (pVF(i,j,k,1).gt.VFhi) then
                           pVF(i,j,k,1)=1.0_WP
                        else
                           ! Update barycenters from moment conservation and project forward
                           if (Lvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then; Lbar_new=(Lbar_old*Lvol_old-Lbar_flux)/Lvol_new; pCL(i,j,k,1:3)=project(Lbar_new,dt); end if
                           if (Gvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then; Gbar_new=(Gbar_old*Gvol_old-Gbar_flux)/Gvol_new; pCG(i,j,k,1:3)=project(Gbar_new,dt); end if
                        end if
                     end if
                  end if
                  ! Divergence of conserved variable fluxes (7 components)
                  rhs(i,j,k,:)=dxi*(pFx(i+1,j,k,:)-pFx(i,j,k,:))+dyi*(pFy(i,j+1,k,:)-pFy(i,j,k,:))+dzi*(pFz(i,j,k+1,:)-pFz(i,j,k,:))
                  ! Velocity gradients at cell center
                  gradU(1,1)=0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
                  gradU(2,1)=0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
                  gradU(3,1)=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
                  gradU(1,2)=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
                  gradU(2,2)=0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
                  gradU(3,2)=0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
                  gradU(1,3)=0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
                  gradU(2,3)=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
                  gradU(3,3)=0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Pressure dilatation: split by VF between phasic energies - discontinuous
                  !rhs(i,j,k,3)=rhs(i,j,k,3)-(       pVF(i,j,k,1))*pPL(i,j,k,1)*div
                  !rhs(i,j,k,4)=rhs(i,j,k,4)-(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)*div
                  ! Viscous heating: τ:∇U, split by VF between phasic energies
                  rhs(i,j,k,3)=rhs(i,j,k,3)+(       pVF(i,j,k,1))*( &
                  & (2.0_WP*pVisc(i,j,k,1)*gradU(1,1)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(1,1) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(2,2)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(2,2) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(3,3)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(3,3) &
                  &+pVisc(i,j,k,1)*(gradU(2,1)+gradU(1,2))*(gradU(2,1)+gradU(1,2)) &
                  &+pVisc(i,j,k,1)*(gradU(3,1)+gradU(1,3))*(gradU(3,1)+gradU(1,3)) &
                  &+pVisc(i,j,k,1)*(gradU(3,2)+gradU(2,3))*(gradU(3,2)+gradU(2,3)))
                  rhs(i,j,k,4)=rhs(i,j,k,4)+(1.0_WP-pVF(i,j,k,1))*( &
                  & (2.0_WP*pVisc(i,j,k,1)*gradU(1,1)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(1,1) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(2,2)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(2,2) &
                  &+(2.0_WP*pVisc(i,j,k,1)*gradU(3,3)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(3,3) &
                  &+pVisc(i,j,k,1)*(gradU(2,1)+gradU(1,2))*(gradU(2,1)+gradU(1,2)) &
                  &+pVisc(i,j,k,1)*(gradU(3,1)+gradU(1,3))*(gradU(3,1)+gradU(1,3)) &
                  &+pVisc(i,j,k,1)*(gradU(3,2)+gradU(2,3))*(gradU(3,2)+gradU(2,3)))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block divergence_and_sources
      this%wt_div=this%wt_div+(MPI_Wtime()-t1)

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

      ! Sync and apply BC
      call this%fill(lvl=this%amr%clvl(),time=time)

      ! Stop full routine timer
      this%wt_dQdt=this%wt_dQdt+(MPI_Wtime()-t0)
   contains

      !> WENO switch function
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP
         real(WP), parameter :: delta=0.01_WP
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight

      !> Recursive subroutine that cuts a tet by grid planes to compute volume and Q fluxes
      recursive subroutine tet2flux(mytet,myind,myVflux,myQflux)
         use amrvof_geometry, only: cut_side,cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert
         implicit none
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  dimension(3,4), intent(in) :: myind
         real(WP), dimension(8),  intent(out) :: myVflux
         real(WP), dimension(7),  intent(out) :: myQflux
         integer :: dir,cut_ind,icase,n1,n2,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         integer,  dimension(3,8,2) :: vert_ind
         real(WP) :: mu,my_vol
         real(WP), dimension(3,4) :: newtet
         integer,  dimension(3,4) :: newind
         real(WP), dimension(3) :: a,b,c
         real(WP), dimension(8) :: subVflux
         real(WP), dimension(7) :: subQflux
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
         use messager, only: die
         implicit none
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  intent(in) :: i0,j0,k0
         real(WP), dimension(8),  intent(out) :: myVflux
         real(WP), dimension(7),  intent(out) :: myQflux
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
            ! Q flux: all mass is liquid
            myQflux=vol_tot*pQold(i0,j0,k0,:)
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas
            myVflux( 2 )=vol_tot
            myVflux(6:8)=vol_tot*bary_tot
            ! Q flux: all mass is gas
            myQflux=vol_tot*pQold(i0,j0,k0,:)
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

         ! Compute Q flux from Qold (guard may be needed at C/F boundaries)
         if (VF0.ge.VFlo) then
            myQflux(1)=myVflux(1)*pQold(i0,j0,k0,1)/VF0
            myQflux(3)=myVflux(1)*pQold(i0,j0,k0,3)/VF0
         end if
         if (VF0.le.VFhi) then
            myQflux(2)=myVflux(2)*pQold(i0,j0,k0,2)/(1.0_WP-VF0)
            myQflux(4)=myVflux(2)*pQold(i0,j0,k0,4)/(1.0_WP-VF0)
         end if
         myQflux(5:7)=sum(myQflux(1:2))*pQold(i0,j0,k0,5:7)/max(sum(pQold(i0,j0,k0,1:2)),this%rho_floor)
         
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

   end subroutine get_dQdt

   !> Override build_plic to include Q clean-up
   subroutine build_plic(this,time)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: time
      ! Call parent build_plic
      call this%amrmpflow%build_plic(time)
      ! Clean up Q
      call this%clean_Q()
   end subroutine build_plic

   !> Clean up Q in pure cells
   subroutine clean_Q(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      ! Traverse all levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pQ =>this%Q%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if      (pVF(i,j,k,1).lt.VFlo) then; pQ(i,j,k,1)=0.0_WP; pQ(i,j,k,3)=0.0_WP
               else if (pVF(i,j,k,1).gt.VFhi) then; pQ(i,j,k,2)=0.0_WP; pQ(i,j,k,4)=0.0_WP
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine clean_Q

   !> Apply relaxation to mixture cells
   subroutine apply_relax(this,time)
      use mpi_f08, only: MPI_Wtime
      use amrvof_geometry, only: get_plane_dist,cut_hex_vol
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl,i,j,k
      real(WP) :: t0
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pCL,pCG,pCurv,pPLIC
      logical :: oldmix,newmix
      real(WP) :: dx,dy,dz,cell_vol,vol_liq,vol_gas
      real(WP), dimension(3) :: lo,hi,bary_liq,bary_gas
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      ! If no relaxation model was provided, return
      if (.not.associated(this%relax)) return
      ! Return if clvl<maxlvl
      if (this%amr%clvl().lt.this%amr%maxlvl) return
      ! Start timer
      t0=MPI_Wtime()
      ! Apply relaxation on finest level only (mixture cells are always at finest)
      lvl=this%amr%maxlvl
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      cell_vol=this%amr%cell_vol(lvl)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pQ   =>this%Q%mf(lvl)%dataptr(mfi)
         pCL  =>this%CL%dataptr(mfi)
         pCG  =>this%CG%dataptr(mfi)
         pCurv=>this%curv%dataptr(mfi)
         pPLIC=>this%plic%dataptr(mfi)
         ! Loop over valid cells
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)

            ! Check if mixture cell prior to relaxation
            oldmix=(pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi)
            ! Apply user-provided relaxation model (modifies VF and Q)
            call this%relax(VF=pVF(i,j,k,1),Q=pQ(i,j,k,:),Pjump=this%sigma*pCurv(i,j,k,1))
            ! Check if mixture cell after relaxation
            newmix=(pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi)

            ! If not mixture cells, clean up and cycle
            if (.not.newmix) then
               if (pVF(i,j,k,1).lt.VFlo) then
                  ! Pure gas
                  pVF(i,j,k,1)=0.0_WP
                  pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pQ(i,j,k,1)=0.0_WP
                  pQ(i,j,k,3)=0.0_WP
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,-1.0e10_WP]
               else if (pVF(i,j,k,1).gt.VFhi) then
                  ! Pure liquid
                  pVF(i,j,k,1)=1.0_WP
                  pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pQ(i,j,k,2)=0.0_WP
                  pQ(i,j,k,4)=0.0_WP
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,+1.0e10_WP]
               end if
               cycle
            end if

            ! If mixture cell, post-process PLIC and barycenters
            if (.not.oldmix) pPLIC(i,j,k,1:3)=[1.0_WP,0.0_WP,0.0_WP]
            ! Adjust PLIC plane to match new VF
            lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            ! Reposition plane: keep normal, adjust distance for new VF
            pPLIC(i,j,k,4)=get_plane_dist(pPLIC(i,j,k,1:3),lo,hi,pVF(i,j,k,1))
            ! Recompute barycenters from adjusted PLIC
            hex(:,1)=[lo(1),lo(2),lo(3)]
            hex(:,2)=[hi(1),lo(2),lo(3)]
            hex(:,3)=[hi(1),hi(2),lo(3)]
            hex(:,4)=[lo(1),hi(2),lo(3)]
            hex(:,5)=[lo(1),lo(2),hi(3)]
            hex(:,6)=[hi(1),lo(2),hi(3)]
            hex(:,7)=[hi(1),hi(2),hi(3)]
            hex(:,8)=[lo(1),hi(2),hi(3)]
            plane=pPLIC(i,j,k,:)
            call cut_hex_vol(hex,plane,vol_liq,vol_gas,bary_liq,bary_gas)
            pVF(i,j,k,1)=vol_liq/cell_vol
            pCL(i,j,k,1:3)=bary_liq
            pCG(i,j,k,1:3)=bary_gas

         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Sync and apply BC
      call this%fill(lvl=this%amr%maxlvl,time=time)
      call this%Q%fill(time=time)
      ! End timer
      this%wt_relax=this%wt_relax+(MPI_Wtime()-t0)
   end subroutine apply_relax

   !> Add artificial bulk viscosity to this%beta and this%visc
   subroutine add_viscartif(this,dt,Cartif,Cvisc)
      use amrsgs, only: get_viscartif
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cartif
      real(WP), intent(in), optional :: Cvisc
      real(WP) :: myCvisc,t0,my_beta
      type(amrdata) :: beta_t
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBeta_t,pBeta,pVisc,pVF,pRHOL,pRHOG
      ! Start timer
      t0=MPI_Wtime()
      ! Set shear viscosity constant
      if (present(Cvisc)) then; myCvisc=Cvisc; else; myCvisc=0.0_WP; end if
      ! Create temp amrdata
      call beta_t%initialize(amr=this%amr,name='beta_t',ncomp=1,ng=this%nover); call beta_t%reset()
      ! Compute kinematic artificial bulk viscosity into temp
      call get_viscartif(dt=dt,visc=beta_t,U=this%UVW,V=this%UVW,W=this%UVW,Ucomp=1,Vcomp=2,Wcomp=3,C=this%C,Cartif=Cartif)
      ! Convert kinematic to dynamic via harmonic-averaged density, add to beta and visc
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pBeta_t=>beta_t%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               my_beta=pBeta_t(i,j,k,1)/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),this%rho_floor)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),this%rho_floor))
               pBeta(i,j,k,1)=pBeta(i,j,k,1)+my_beta
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+my_beta*myCvisc
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temp amrdata
      call beta_t%finalize()
      ! End timer
      this%wt_visc=this%wt_visc+(MPI_Wtime()-t0)
   end subroutine add_viscartif

   !> Add Vreman SGS eddy viscosity to this%visc
   subroutine add_vreman(this,dt,Cs)
      use amrsgs, only: get_vreman
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrdata) :: visc_t
      real(WP) :: t0
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pVisc,pVF,pRHOL,pRHOG
      ! Start timer
      t0=MPI_Wtime()
      ! Create temp amrdata
      call visc_t%initialize(amr=this%amr,name='visc_t',ncomp=1,ng=this%nover); call visc_t%reset()
      ! Compute kinematic eddy viscosity into temp
      call get_vreman(dt=dt,visc=visc_t,U=this%UVW,V=this%UVW,W=this%UVW,Ucomp=1,Vcomp=2,Wcomp=3,Cs=Cs)
      ! Convert kinematic to dynamic via harmonic-averaged density, add to visc
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),this%rho_floor)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),this%rho_floor))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temp amrdata
      call visc_t%finalize()
      ! End timer
      this%wt_visc=this%wt_visc+(MPI_Wtime()-t0)
   end subroutine add_vreman

   !> Calculate CFL numbers
   subroutine get_cfl(this,dt,cfl)
      use mathtools, only: Pi
      use parallel,  only: MPI_REAL_WP
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: lvl,i,j,k,ierr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pUVW,pVisc,pBeta,pDiff,pVF,pPL,pPG,pC,pTL,pTG
      real(WP) :: dxi,dyi,dzi,rho,viscmax,conv,pgrad
      real(WP) :: Pmix_ip,Pmix_im,Pmix_jp,Pmix_jm,Pmix_kp,Pmix_km
      ! Get convective CFL from parent
      call this%amrmpflow%get_cflc(dt=dt)
      ! Reset child CFLs
      this%CFLp=0.0_WP; this%CFLst=0.0_WP
      this%CFLa_x=0.0_WP; this%CFLa_y=0.0_WP; this%CFLa_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute CFL at each level
      do lvl=0,this%amr%clvl()
         ! Get mesh spacing
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get pointers to data
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pUVW =>this%UVW%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pDiff=>this%diff%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pPL  =>this%PL%mf(lvl)%dataptr(mfi)
            pPG  =>this%PG%mf(lvl)%dataptr(mfi)
            pC   =>this%C%mf(lvl)%dataptr(mfi)
            pTL  =>this%TL%mf(lvl)%dataptr(mfi)
            pTG  =>this%TG%mf(lvl)%dataptr(mfi)
            ! Loop over cells
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Get density
               rho=max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
               ! Viscous CFL
               viscmax=max(pVisc(i,j,k,1)/rho, &
               &           pBeta(i,j,k,1)/rho, &
               &           pDiff(i,j,k,1)*(pVF(i,j,k,1)*pTL(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*pTG(i,j,k,1))/max(pQ(i,j,k,3)+pQ(i,j,k,4),this%rho_floor))
               if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt*dxi**2)
               if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt*dyi**2)
               if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt*dzi**2)
               ! Convective+pressure CFL
               conv=abs(pUVW(i,j,k,1))*dxi+abs(pUVW(i,j,k,2))*dyi+abs(pUVW(i,j,k,3))*dzi
               pgrad=0.0_WP
               if (this%amr%nx.gt.1) then
                  Pmix_ip=pVF(i+1,j,k,1)*pPL(i+1,j,k,1)+(1.0_WP-pVF(i+1,j,k,1))*pPG(i+1,j,k,1)
                  Pmix_im=pVF(i-1,j,k,1)*pPL(i-1,j,k,1)+(1.0_WP-pVF(i-1,j,k,1))*pPG(i-1,j,k,1)
                  pgrad=pgrad+abs(Pmix_ip-Pmix_im)*0.5_WP*dxi**2
               end if
               if (this%amr%ny.gt.1) then
                  Pmix_jp=pVF(i,j+1,k,1)*pPL(i,j+1,k,1)+(1.0_WP-pVF(i,j+1,k,1))*pPG(i,j+1,k,1)
                  Pmix_jm=pVF(i,j-1,k,1)*pPL(i,j-1,k,1)+(1.0_WP-pVF(i,j-1,k,1))*pPG(i,j-1,k,1)
                  pgrad=pgrad+abs(Pmix_jp-Pmix_jm)*0.5_WP*dyi**2
               end if
               if (this%amr%nz.gt.1) then
                  Pmix_kp=pVF(i,j,k+1,1)*pPL(i,j,k+1,1)+(1.0_WP-pVF(i,j,k+1,1))*pPG(i,j,k+1,1)
                  Pmix_km=pVF(i,j,k-1,1)*pPL(i,j,k-1,1)+(1.0_WP-pVF(i,j,k-1,1))*pPG(i,j,k-1,1)
                  pgrad=pgrad+abs(Pmix_kp-Pmix_km)*0.5_WP*dzi**2
               end if
               pgrad=pgrad/rho
               this%CFLp=max(this%CFLp,0.5_WP*dt*(conv+sqrt(conv**2+4.0_WP*pgrad)))
               ! Acoustic CFL
               if (this%amr%nx.gt.1) this%CFLa_x=max(this%CFLa_x,(abs(pUVW(i,j,k,1))+pC(i,j,k,1))*dt*dxi)
               if (this%amr%ny.gt.1) this%CFLa_y=max(this%CFLa_y,(abs(pUVW(i,j,k,2))+pC(i,j,k,1))*dt*dyi)
               if (this%amr%nz.gt.1) this%CFLa_z=max(this%CFLa_z,(abs(pUVW(i,j,k,3))+pC(i,j,k,1))*dt*dzi)
               ! Surface tension CFL (only at finest level, near interface)
               if (this%sigma.gt.0.0_WP.and.lvl.eq.this%amr%maxlvl.and.pVF(i,j,k,1).gt.VFlo.and.pVF(i,j,k,1).lt.VFhi) then
                  this%CFLst=max(this%CFLst,dt/sqrt(rho*this%amr%min_meshsize(lvl)**3/(4.0_WP*Pi*this%sigma)))
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Reduce across ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLst ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_x,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_y,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_z,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z,this%CFLst)
      if (.not.this%use_projection) cfl=max(cfl,this%CFLa_x,this%CFLa_y,this%CFLa_z)
   end subroutine get_cfl

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Calculate monitoring info
   subroutine get_info(this)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MAX,MPI_MIN
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,n,ierr
      real(WP) :: dV

      ! Use parent's method first
      call this%amrmpflow%get_info()

      ! VF-conditional phasic extrema
      phasic_extrema: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG
         integer :: i,j,k
         ! Initialize extrema
         this%RHOLmin=huge(1.0_WP); this%RHOLmax=-huge(1.0_WP); this%RHOGmin=huge(1.0_WP); this%RHOGmax=-huge(1.0_WP)
         this%ILmin=huge(1.0_WP); this%ILmax=-huge(1.0_WP); this%IGmin=huge(1.0_WP); this%IGmax=-huge(1.0_WP)
         this%PLmin=huge(1.0_WP); this%PLmax=-huge(1.0_WP); this%PGmin=huge(1.0_WP); this%PGmax=-huge(1.0_WP)
         this%TLmin=huge(1.0_WP); this%TLmax=-huge(1.0_WP); this%TGmin=huge(1.0_WP); this%TGmax=-huge(1.0_WP)
         this%Cmin=huge(1.0_WP); this%Cmax=-huge(1.0_WP)
         this%dPmax=0.0_WP
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Velocity norm 0
            this%Umax=max(this%Umax,this%UVW%norm0(lvl=lvl,comp=1))
            this%Vmax=max(this%Vmax,this%UVW%norm0(lvl=lvl,comp=2))
            this%Wmax=max(this%Wmax,this%UVW%norm0(lvl=lvl,comp=3))
            ! Extrema of mixture speed of sound
            this%Cmin=min(this%Cmin,this%C%get_min(lvl=lvl)); this%Cmax=max(this%Cmax,this%C%get_max(lvl=lvl))
            ! Build fine mask for this level (if not finest)
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            ! Manual loops for discontinuous variables
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get data pointers
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
               pIL=>this%IL%mf(lvl)%dataptr(mfi); pIG=>this%IG%mf(lvl)%dataptr(mfi)
               pPL=>this%PL%mf(lvl)%dataptr(mfi); pPG=>this%PG%mf(lvl)%dataptr(mfi)
               pTL=>this%TL%mf(lvl)%dataptr(mfi); pTG=>this%TG%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over interior tiles
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
                  ! Liquid variables
                  if (pVF(i,j,k,1).ge.VFlo) then
                     this%RHOLmin=min(this%RHOLmin,pRHOL(i,j,k,1)); this%RHOLmax=max(this%RHOLmax,pRHOL(i,j,k,1))
                     this%ILmin  =min(this%ILmin  ,pIL  (i,j,k,1)); this%ILmax  =max(this%ILmax  ,pIL  (i,j,k,1))
                     this%PLmin  =min(this%PLmin  ,pPL  (i,j,k,1)); this%PLmax  =max(this%PLmax  ,pPL  (i,j,k,1))
                     this%TLmin  =min(this%TLmin  ,pTL  (i,j,k,1)); this%TLmax  =max(this%TLmax  ,pTL  (i,j,k,1))
                  end if
                  ! Gas variables
                  if (pVF(i,j,k,1).le.VFhi) then
                     this%RHOGmin=min(this%RHOGmin,pRHOG(i,j,k,1)); this%RHOGmax=max(this%RHOGmax,pRHOG(i,j,k,1))
                     this%IGmin  =min(this%IGmin  ,pIG  (i,j,k,1)); this%IGmax  =max(this%IGmax  ,pIG  (i,j,k,1))
                     this%PGmin  =min(this%PGmin  ,pPG  (i,j,k,1)); this%PGmax  =max(this%PGmax  ,pPG  (i,j,k,1))
                     this%TGmin  =min(this%TGmin  ,pTG  (i,j,k,1)); this%TGmax  =max(this%TGmax  ,pTG  (i,j,k,1))
                  end if
                  ! Pressure gap in mixed cells
                  if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
                     this%dPmax=max(this%dPmax,abs(pPL(i,j,k,1)-pPG(i,j,k,1)))
                  end if
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         ! Reduce across MPI ranks
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOLmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOLmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%ILmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%ILmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%PLmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%PLmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%TLmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%TLmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOGmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOGmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%IGmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%IGmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%PGmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%PGmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%TGmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%TGmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%dPmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      end block phasic_extrema

      ! Kinetic energy integral: 0.5 * rho * (U^2 + V^2 + W^2) * dV
      ! Uses composite integration with fine masking to avoid double-counting
      get_rhoKint: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pUVW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         integer :: i,j,k
         this%rhoKint=0.0_WP
         do lvl=0,this%amr%clvl()
            ! Get cell volume
            dV=this%amr%cell_vol(lvl)
            ! Build fine mask for this level (if not finest)
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               bx=mfi%tilebox()
               ! Get pointers to data
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
               ! Get pointer to fine mask
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over cells
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
                  ! Accumulate kinetic energy (rho = Q(1)+Q(2) = VF*rhoL + (1-VF)*rhoG)
                  this%rhoKint=this%rhoKint+0.5_WP*(pQ(i,j,k,1)+pQ(i,j,k,2))*(pUVW(i,j,k,1)**2+pUVW(i,j,k,2)**2+pUVW(i,j,k,3)**2)*dV
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         ! Reduce across MPI ranks
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_rhoKint

      ! Reduce per-rank timing to min/max across ranks
      call MPI_ALLREDUCE(this%wt_prim, this%wtmax_prim, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_prim, this%wtmin_prim, 1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_dQdt, this%wtmax_dQdt, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_dQdt, this%wtmin_dQdt, 1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_fv,   this%wtmax_fv,   1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_fv,   this%wtmin_fv,   1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_div,  this%wtmax_div,  1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_div,  this%wtmin_div,  1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_relax,this%wtmax_relax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_relax,this%wtmin_relax,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc, this%wtmax_visc, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc, this%wtmin_visc, 1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      ! Reset per-rank timing accumulators for next interval
      this%wt_prim=0.0_WP; this%wt_dQdt=0.0_WP; this%wt_fv=0.0_WP; this%wt_div=0.0_WP; this%wt_relax=0.0_WP; this%wt_visc=0.0_WP

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrmpcomp_print(this)
      use messager, only: log
      implicit none
      class(amrmpcomp), intent(in) :: this
      call log("Compressible Multiphase solver: "//trim(this%name))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrmpcomp_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      ! Register face velocities, conserved variables, and VOF data via parent
      call this%amrmpflow%register_checkpoint(io)
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      ! Restore face velocities and conserved variables via parent
      call this%amrmpflow%restore_checkpoint(io,dirname,time)
      ! Rebuild primitive variables
      call this%get_primitive(this%Q)
   end subroutine restore_checkpoint

end module amrmpcomp_class
