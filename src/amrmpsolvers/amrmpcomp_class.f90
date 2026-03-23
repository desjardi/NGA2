!> AMR compressible multiphase solver class
!> Inherits from amrvof_class
module amrmpcomp_class
   use iso_c_binding,    only: c_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap
   use amrvof_class,     only: amrvof,VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
   implicit none
   private

   ! Expose type and constants
   public :: amrmpcomp,VFlo,VFhi,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER

   !> AMR compressible multiphase solver type
   type, extends(amrvof) :: amrmpcomp

      ! User-configurable callbacks
      procedure(mpcomp_init_iface   ), pointer, pass :: user_mpcomp_init   =>null()
      procedure(mpcomp_tagging_iface), pointer, pass :: user_mpcomp_tagging=>null()
      procedure(mpcomp_bc_iface     ), pointer, pass :: user_mpcomp_bc     =>null()

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

      ! Conserved variables (1: VF*rhoL, 2: (1-VF)*rhoG, 3: VF*rhoL*IL, 4: (1-VF)*rhoG*IG, 5: rhoU, 6: rhoV, 7: rhoW)
      type(amrdata) :: Q,Qold

      ! Primitive variables: mixture velocity
      type(amrdata) :: U,V,W

      ! Phasic primitive variables
      type(amrdata) :: RHOL,RHOG           !< Phasic densities
      type(amrdata) :: IL,IG               !< Phasic internal energies
      type(amrdata) :: PL,PG               !< Phasic pressures
      type(amrdata) :: TL,TG               !< Phasic temperatures

      ! Mixture properties
      type(amrdata) :: C                   !< Speed of sound

      ! Physical properties
      type(amrdata) :: visc              !< Dynamic viscosity
      type(amrdata) :: beta              !< Bulk viscosity
      type(amrdata) :: diff              !< Heat diffusivity

      ! CFL numbers
      real(WP) :: CFLc_x=0.0_WP,CFLc_y=0.0_WP,CFLc_z=0.0_WP  !< Convective
      real(WP) :: CFLa_x=0.0_WP,CFLa_y=0.0_WP,CFLa_z=0.0_WP  !< Acoustic
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous

      ! Monitoring quantities
      real(WP) :: Umax=0.0_WP,Vmax=0.0_WP,Wmax=0.0_WP
      real(WP) :: RHOLmin=0.0_WP,RHOLmax=0.0_WP,RHOGmin=0.0_WP,RHOGmax=0.0_WP
      real(WP) :: ILmin=0.0_WP,ILmax=0.0_WP,IGmin=0.0_WP,IGmax=0.0_WP
      real(WP) :: PLmin=0.0_WP,PLmax=0.0_WP,PGmin=0.0_WP,PGmax=0.0_WP
      real(WP) :: TLmin=0.0_WP,TLmax=0.0_WP,TGmin=0.0_WP,TGmax=0.0_WP
      real(WP) :: Cmin=0.0_WP,Cmax=0.0_WP
      real(WP) :: dPmax=0.0_WP
      real(WP), dimension(7) :: Qint=0.0_WP,Qmin=0.0_WP,Qmax=0.0_WP
      real(WP) :: rhoKint=0.0_WP

      ! Minimum density for stability
      real(WP) :: rho_floor=1.0e-10_WP

      ! Cost for mixed cells in load balancing (1.0 = same as pure, higher = more expensive)
      real(WP) :: SLcost=100.0_WP

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
      procedure :: tagging
      procedure :: get_cost
      ! Physics methods
      procedure :: store_old
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

   !> Abstract interface for user init callback
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

   !> Abstract interface for tagging callback
   abstract interface
      subroutine mpcomp_tagging_iface(solver,lvl,time,tags)
         import :: amrmpcomp,c_ptr,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine mpcomp_tagging_iface
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
   ! abstract interface
   !    subroutine relax_iface(VF,Q,Peq,dVF,dQ)
   !       import :: WP
   !       real(WP), intent(in) :: VF
   !       real(WP), dimension(:), intent(in) :: Q
   !       real(WP), intent(inout) :: Peq
   !       real(WP), intent(out), optional :: dVF
   !       real(WP), dimension(:), intent(out), optional :: dQ
   !    end subroutine relax_iface
   ! end interface

   !> Abstract interface for pressure relaxation callback
   abstract interface
      subroutine relax_iface(VF,Q)
         import :: WP
         real(WP), intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
      end subroutine relax_iface
   end interface

   !> Abstract interface for user BC callback
   abstract interface
      subroutine mpcomp_bc_iface(solver,lvl,time,face,bx,pQ)
         use amrex_amr_module, only: amrex_box
         import :: amrmpcomp,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      end subroutine mpcomp_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrmpcomp type
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
      if (associated(this%user_mpcomp_init)) call this%user_mpcomp_init(lvl,time,ba,dm)
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
      if (associated(this%user_mpcomp_tagging)) call this%user_mpcomp_tagging(lvl,time,tags)
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
      use amrex_amr_module, only: amrex_bc_foextrap
      implicit none
      class(amrmpcomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name
      
      ! Initialize amrvof parent without callback registration
      this%amrvof%skip_registration=.true.
      call this%amrvof%initialize(amr,name)
      
      ! Initialize conserved variables Q (7 components for multiphase)
      call this%Q%initialize   (amr,name='Q'   ,ncomp=7,ng=this%nover); this%Q%parent   =>this
      call this%Qold%initialize(amr,name='Qold',ncomp=7,ng=this%nover); this%Qold%parent=>this
      
      ! Initialize mixture velocity
      call this%U%initialize(amr,name='U',ncomp=1,ng=this%nover); this%U%parent=>this
      call this%V%initialize(amr,name='V',ncomp=1,ng=this%nover); this%V%parent=>this
      call this%W%initialize(amr,name='W',ncomp=1,ng=this%nover); this%W%parent=>this
      
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

      ! Set Q fillbc callback to internal handler
      this%Q%fillbc=>Q_fillbc

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
      ! Conserved variables
      call this%Q%finalize(); call this%Qold%finalize()
      ! Velocity
      call this%U%finalize(); call this%V%finalize(); call this%W%finalize()
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
      nullify(this%user_mpcomp_init); nullify(this%user_mpcomp_tagging); nullify(this%user_mpcomp_bc)
      nullify(this%getPL); nullify(this%getCL); nullify(this%getTL)
      nullify(this%getPG); nullify(this%getCG); nullify(this%getTG)
      nullify(this%relax)
      ! Finalize parent
      call this%amrvof%finalize()
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles VF/VFold + CL/CG/PLIC multifabs
      call this%amrvof%on_init(lvl,time,ba,dm)
      ! Reset level layouts for conserved variables
      call this%Q%reset_level(lvl,ba,dm); call this%Qold%reset_level(lvl,ba,dm)
      ! Reset mixture velocity
      call this%U%reset_level(lvl,ba,dm); call this%V%reset_level(lvl,ba,dm); call this%W%reset_level(lvl,ba,dm)
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
      call this%Q%setval(val=0.0_WP,lvl=lvl); call this%Qold%setval(val=0.0_WP,lvl=lvl)
      call this%U%setval(val=0.0_WP,lvl=lvl); call this%V%setval(val=0.0_WP,lvl=lvl); call this%W%setval(val=0.0_WP,lvl=lvl)
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
      ! Parent handles VF interpolation and CL/CG/PLIC rebuild
      call this%amrvof%on_coarse(lvl,time,ba,dm)
      ! Conserved variables
      call this%Q%on_coarse(lvl,time,ba,dm)
      call this%Qold%reset_level(lvl,ba,dm)
      ! Auxiliary / derived quantities (just reset, will be recomputed)
      call this%U%reset_level(lvl,ba,dm); call this%V%reset_level(lvl,ba,dm); call this%W%reset_level(lvl,ba,dm)
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
      ! Parent handles VF remake and CL/CG/PLIC parallel_copy
      call this%amrvof%on_remake(lvl,time,ba,dm)
      ! Conserved variables
      call this%Q%on_remake(lvl,time,ba,dm)
      call this%Qold%reset_level(lvl,ba,dm)
      ! Auxiliary / derived quantities (just reset, will be recomputed)
      call this%U%reset_level(lvl,ba,dm); call this%V%reset_level(lvl,ba,dm); call this%W%reset_level(lvl,ba,dm)
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
      ! Parent handles VF/VFold + CL/CG/PLIC
      call this%amrvof%on_clear(lvl)
      ! Conserved variables
      call this%Q%clear_level(lvl); call this%Qold%clear_level(lvl)
      ! Auxiliary / derived quantities
      call this%U%clear_level(lvl); call this%V%clear_level(lvl); call this%W%clear_level(lvl)
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
      ! Parent handles VF average-down + fill
      call this%amrvof%post_regrid(lbase,time)
      ! Average down conserved variables Q for C/F consistency
      call this%Q%average_down(lbase)
      ! Fill Q ghosts and rebuild primitives
      call this%Q%fill(time)
      call this%get_primitive(this%Q)
   end subroutine post_regrid

   !> Tag cells near interface with regrid_buffer layer growth
   subroutine tagging(this,lvl,time,tags)
      implicit none
      class(amrmpcomp), intent(inout) :: this
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
      class(amrmpcomp), intent(inout) :: this
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

   !> Internal fillbc for Q - calls default_fillbc first, then user_mpcomp_bc for ext_dir faces
   subroutine Q_fillbc(this,mf,scomp,ncomp,time,geom)
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
      class(amrmpcomp), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer, dimension(3) :: dlo,dhi
      integer :: lvl
      
      ! First apply default BC handling (foextrap,hoextrap,reflect,etc.)
      call default_fillbc(this,mf,scomp,ncomp,time,geom)
      
      ! Access parent solver
      select type (s=>this%parent)
       class is (amrmpcomp)
         solver=>s
      end select
      
      ! Check if user callback exists
      if (.not.associated(solver%user_mpcomp_bc)) return

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
         if (this%lo_bc(1,1).eq.amrex_bc_ext_dir .and. ilo.lt.dlo(1)) then
            bc_bx=amrex_box([ilo,jlo,klo],[dlo(1)-1,jhi,khi])
            call solver%user_mpcomp_bc(lvl=lvl,time=time,face=1,bx=bc_bx,pQ=p)
         end if
         ! X-HIGH (face=2)
         if (this%hi_bc(1,1).eq.amrex_bc_ext_dir .and. ihi.gt.dhi(1)) then
            bc_bx=amrex_box([dhi(1)+1,jlo,klo],[ihi,jhi,khi])
            call solver%user_mpcomp_bc(lvl=lvl,time=time,face=2,bx=bc_bx,pQ=p)
         end if
         ! Y-LOW (face=3)
         if (this%lo_bc(2,1).eq.amrex_bc_ext_dir .and. jlo.lt.dlo(2)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,dlo(2)-1,khi])
            call solver%user_mpcomp_bc(lvl=lvl,time=time,face=3,bx=bc_bx,pQ=p)
         end if
         ! Y-HIGH (face=4)
         if (this%hi_bc(2,1).eq.amrex_bc_ext_dir .and. jhi.gt.dhi(2)) then
            bc_bx=amrex_box([ilo,dhi(2)+1,klo],[ihi,jhi,khi])
            call solver%user_mpcomp_bc(lvl=lvl,time=time,face=4,bx=bc_bx,pQ=p)
         end if
         ! Z-LOW (face=5)
         if (this%lo_bc(3,1).eq.amrex_bc_ext_dir .and. klo.lt.dlo(3)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,jhi,dlo(3)-1])
            call solver%user_mpcomp_bc(lvl=lvl,time=time,face=5,bx=bc_bx,pQ=p)
         end if
         ! Z-HIGH (face=6)
         if (this%hi_bc(3,1).eq.amrex_bc_ext_dir .and. khi.gt.dhi(3)) then
            bc_bx=amrex_box([ilo,jlo,dhi(3)+1],[ihi,jhi,khi])
            call solver%user_mpcomp_bc(lvl=lvl,time=time,face=6,bx=bc_bx,pQ=p)
         end if
      end do
      call amrex_mfiter_destroy(mfi)

   end subroutine Q_fillbc

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Copy current state to old state
   subroutine store_old(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      ! Store interface state
      call this%amrvof%store_old()
      ! Store conserved variables
      call this%Qold%copy(src=this%Q)
   end subroutine store_old

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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG,pC
      real(WP) :: rho_inv,CL,CG
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
            pU   =>this%U%mf(lvl)%dataptr(mfi)
            pV   =>this%V%mf(lvl)%dataptr(mfi)
            pW   =>this%W%mf(lvl)%dataptr(mfi)
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
               rho_inv=1.0_WP/max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
               pU(i,j,k,1)=pQ(i,j,k,5)*rho_inv
               pV(i,j,k,1)=pQ(i,j,k,6)*rho_inv
               pW(i,j,k,1)=pQ(i,j,k,7)*rho_inv
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
               pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*CL**2+pQ(i,j,k,2)*CG**2)*rho_inv)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! End timer
      this%wt_prim=this%wt_prim+(MPI_Wtime()-t0)
   end subroutine get_primitive

   !> Calculate conserved variables from primitive variables
   !> Rebuilds Q from VF, phasic densities, phasic energies, and mixture velocity
   subroutine get_conserved(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pU,pV,pW,pRHOL,pRHOG,pIL,pIG
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ   =>this%Q%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            pU   =>this%U%mf(lvl)%dataptr(mfi)
            pV   =>this%V%mf(lvl)%dataptr(mfi)
            pW   =>this%W%mf(lvl)%dataptr(mfi)
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
               pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pU(i,j,k,1)
               pQ(i,j,k,6)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pV(i,j,k,1)
               pQ(i,j,k,7)=(pQ(i,j,k,1)+pQ(i,j,k,2))*pW(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_conserved

   !> Calculate dQdt from passed Q
   subroutine get_dQdt(this,Q,dQdt,dt,time)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: Q
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

      ! First build primitive variables from Q
      call this%get_primitive(Q)

      ! Build transport band at finest level to localize SL computation
      call this%amr%mfab_build(lvl=this%amr%maxlvl,mfab=band,ncomp=1,nover=1)
      call this%build_band(lvl=this%amr%maxlvl,VF=this%VFold%mf(this%amr%maxlvl),band=band,nband=2)

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
         call this%amr%mfab_build(lvl=this%amr%maxlvl,mfab=Vx,ncomp=8,nover=0,atface=[.true. ,.false.,.false.]); call Vx%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%maxlvl,mfab=Vy,ncomp=8,nover=0,atface=[.false.,.true. ,.false.]); call Vy%setval(0.0_WP)
         call this%amr%mfab_build(lvl=this%amr%maxlvl,mfab=Vz,ncomp=8,nover=0,atface=[.false.,.false.,.true. ]); call Vz%setval(0.0_WP)
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
         real(WP) :: vel
         real(WP), dimension(8) :: Vflux
         real(WP), dimension(7) :: Qflux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pFx,pFy,pFz
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx,nbx
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
               ! Build face velocity
               vel=0.5_WP*sum(pU(i-1:i,j,k,1))
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,5)=proj(:,i,j  ,k  )
               face(:,2)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=proj(:,i,j  ,k+1)
               face(:,3)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,7)=proj(:,i,j+1,k+1)
               face(:,4)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=proj(:,i,j+1,k  )
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dy*dz*vel)
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(1,nn)=merge(i-1,i,vel.gt.0.0_WP); end do
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
                  pFx(i,j,k,5)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pU(i-1:i,j,k,1))
                  pFx(i,j,k,6)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pV(i-1:i,j,k,1))
                  pFx(i,j,k,7)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pW(i-1:i,j,k,1))
               end if
            end do; end do; end do
            ! Y-fluxes
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               ! Build face velocity
               vel=0.5_WP*sum(pV(i,j-1:j,k,1))
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,5)=proj(:,i+1,j,k+1)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=proj(:,i  ,j,k+1)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,7)=proj(:,i  ,j,k  )
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=proj(:,i+1,j,k  )
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dz*dx*vel)
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(2,nn)=merge(j-1,j,vel.gt.0.0_WP); end do
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
                  pFy(i,j,k,5)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pU(i,j-1:j,k,1))
                  pFy(i,j,k,6)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pV(i,j-1:j,k,1))
                  pFy(i,j,k,7)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pW(i,j-1:j,k,1))
               end if
            end do; end do; end do
            ! Z-fluxes
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               ! Build face velocity
               vel=0.5_WP*sum(pW(i,j,k-1:k,1))
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,5)=proj(:,i+1,j  ,k)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,6)=proj(:,i  ,j  ,k)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,7)=proj(:,i  ,j+1,k)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,8)=proj(:,i+1,j+1,k)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dx*dy*vel)
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(3,nn)=merge(k-1,k,vel.gt.0.0_WP); end do
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
                  pFz(i,j,k,5)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pU(i,j,k-1:k,1))
                  pFz(i,j,k,6)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pV(i,j,k-1:k,1))
                  pFz(i,j,k,7)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pW(i,j,k-1:k,1))
               end if
            end do; end do; end do
            ! Deallocate proj for this tile
            deallocate(proj)
         end do
         call this%amr%mfiter_destroy(mfi)
      end block semilagrangian_fluxes
      this%wt_sl=this%wt_sl+(MPI_Wtime()-t1)
      
      ! Phase 1b: Finite volume fluxes for all levels (Euler fluxes skip band cells at finest level)
      t1=MPI_Wtime()
      finitevolume_fluxes: block
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pFx,pFy,pFz,pBand
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pPL,pPG,pTL,pTG,pIL,pIG,pVF
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc,pBeta,pDiff
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW ! Intentional masking
         real(WP), dimension(-2: 0) :: wenop
         real(WP), dimension(-1:+1) :: wenom
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: w,div,vel,visc_f,beta_f,irho_f
         real(WP), dimension(-2:+1) :: Pmix
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
               pQ   =>Q%mf(lvl)%dataptr(mfi)
               pU   =>this%U%mf(lvl)%dataptr(mfi)
               pV   =>this%V%mf(lvl)%dataptr(mfi)
               pW   =>this%W%mf(lvl)%dataptr(mfi)
               pVF  =>this%VF%mf(lvl)%dataptr(mfi)
               pPL  =>this%PL%mf(lvl)%dataptr(mfi)
               pPG  =>this%PG%mf(lvl)%dataptr(mfi)
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
                     ! Face velocity
                     vel=0.5_WP*sum(pU(i-1:i,j,k,1))
                     ! Rhie-Chow correction
                     irho_f=2.0_WP/max(sum(pQ(i-1:i,j,k,1:2)),this%rho_floor)
                     Pmix=pVF(i-2:i+1,j,k,1)*pPL(i-2:i+1,j,k,1)+(1.0_WP-pVF(i-2:i+1,j,k,1))*pPG(i-2:i+1,j,k,1)
                     vel=vel+0.25_WP*dt*dxi*irho_f*(Pmix(1)-3.0_WP*Pmix(0)+3.0_WP*Pmix(-1)-Pmix(-2))
                     ! WENO liquid mass and energy fluxes
                     if (any(pVF(i-1:i,j,k,1).ge.VFlo)) then
                        w=weno_weight((abs(pQ(i-1,j,k,1)-pQ(i-2,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i+1,j,k,1)-pQ(i  ,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i-2:i  ,j,k,1)) &
                        &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i-1:i+1,j,k,1))
                        w=weno_weight((abs(pIL(i-1,j,k,1)-pIL(i-2,j,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIL(i+1,j,k,1)-pIL(i  ,j,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,3)=0.5_WP*(pFx(i,j,k,1)-abs(pFx(i,j,k,1)))*sum(wenop*pIL(i-2:i  ,j,k,1)) &
                        &           +0.5_WP*(pFx(i,j,k,1)+abs(pFx(i,j,k,1)))*sum(wenom*pIL(i-1:i+1,j,k,1))
                     end if
                     ! WENO gas mass and energy fluxes
                     if (any(pVF(i-1:i,j,k,1).le.VFhi)) then
                        w=weno_weight((abs(pQ(i-1,j,k,2)-pQ(i-2,j,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i-1,j,k,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i+1,j,k,2)-pQ(i  ,j,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i-1,j,k,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,2)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i-2:i  ,j,k,2)) &
                        &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i-1:i+1,j,k,2))
                        w=weno_weight((abs(pIG(i-1,j,k,1)-pIG(i-2,j,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIG(i+1,j,k,1)-pIG(i  ,j,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFx(i,j,k,4)=0.5_WP*(pFx(i,j,k,2)-abs(pFx(i,j,k,2)))*sum(wenop*pIG(i-2:i  ,j,k,1)) &
                        &           +0.5_WP*(pFx(i,j,k,2)+abs(pFx(i,j,k,2)))*sum(wenom*pIG(i-1:i+1,j,k,1))
                     end if
                     ! Momentum fluxes
                     pFx(i,j,k,5)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pU(i-1:i,j,k,1))
                     pFx(i,j,k,6)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pV(i-1:i,j,k,1))
                     pFx(i,j,k,7)=sum(pFx(i,j,k,1:2))*0.5_WP*sum(pW(i-1:i,j,k,1))
                  end if
                  ! Add pressure stress: P_face = VF*PL + (1-VF)*PG
                  pFx(i,j,k,5)=pFx(i,j,k,5)-0.5_WP*sum(pVF(i-1:i,j,k,1)*pPL(i-1:i,j,k,1)+(1.0_WP-pVF(i-1:i,j,k,1))*pPG(i-1:i,j,k,1))
                  ! Velocity gradients at x-face
                  gradU(1,1)=dxi*(pU(i,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*(pU(i-1,j+1,k,1)-pU(i-1,j-1,k,1)+pU(i,j+1,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*(pU(i-1,j,k+1,1)-pU(i-1,j,k-1,1)+pU(i,j,k+1,1)-pU(i,j,k-1,1))
                  gradU(1,2)=dxi*(pV(i,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=0.25_WP*dyi*(pV(i-1,j+1,k,1)-pV(i-1,j-1,k,1)+pV(i,j+1,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=0.25_WP*dzi*(pV(i-1,j,k+1,1)-pV(i-1,j,k-1,1)+pV(i,j,k+1,1)-pV(i,j,k-1,1))
                  gradU(1,3)=dxi*(pW(i,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=0.25_WP*dyi*(pW(i-1,j+1,k,1)-pW(i-1,j-1,k,1)+pW(i,j+1,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=0.25_WP*dzi*(pW(i-1,j,k+1,1)-pW(i-1,j,k-1,1)+pW(i,j,k+1,1)-pW(i,j,k-1,1))
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
                     ! Face velocity
                     vel=0.5_WP*sum(pV(i,j-1:j,k,1))
                     ! Rhie-Chow correction
                     irho_f=2.0_WP/max(sum(pQ(i,j-1:j,k,1:2)),this%rho_floor)
                     Pmix=pVF(i,j-2:j+1,k,1)*pPL(i,j-2:j+1,k,1)+(1.0_WP-pVF(i,j-2:j+1,k,1))*pPG(i,j-2:j+1,k,1)
                     vel=vel+0.25_WP*dt*dyi*irho_f*(Pmix(1)-3.0_WP*Pmix(0)+3.0_WP*Pmix(-1)-Pmix(-2))
                     ! WENO liquid mass and energy fluxes
                     if (any(pVF(i,j-1:j,k,1).ge.VFlo)) then
                        w=weno_weight((abs(pQ(i,j-1,k,1)-pQ(i,j-2,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j+1,k,1)-pQ(i,j  ,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j-2:j  ,k,1)) &
                        &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j-1:j+1,k,1))
                        w=weno_weight((abs(pIL(i,j-1,k,1)-pIL(i,j-2,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIL(i,j+1,k,1)-pIL(i,j  ,k,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,3)=0.5_WP*(pFy(i,j,k,1)-abs(pFy(i,j,k,1)))*sum(wenop*pIL(i,j-2:j  ,k,1)) &
                        &           +0.5_WP*(pFy(i,j,k,1)+abs(pFy(i,j,k,1)))*sum(wenom*pIL(i,j-1:j+1,k,1))
                     end if
                     ! WENO gas mass and energy fluxes
                     if (any(pVF(i,j-1:j,k,1).le.VFhi)) then
                        w=weno_weight((abs(pQ(i,j-1,k,2)-pQ(i,j-2,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j-1,k,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j+1,k,2)-pQ(i,j  ,k,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j-1,k,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,2)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j-2:j  ,k,2)) &
                        &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j-1:j+1,k,2))
                        w=weno_weight((abs(pIG(i,j-1,k,1)-pIG(i,j-2,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIG(i,j+1,k,1)-pIG(i,j  ,k,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFy(i,j,k,4)=0.5_WP*(pFy(i,j,k,2)-abs(pFy(i,j,k,2)))*sum(wenop*pIG(i,j-2:j  ,k,1)) &
                        &           +0.5_WP*(pFy(i,j,k,2)+abs(pFy(i,j,k,2)))*sum(wenom*pIG(i,j-1:j+1,k,1))
                     end if
                     ! Momentum fluxes
                     pFy(i,j,k,5)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pU(i,j-1:j,k,1))
                     pFy(i,j,k,6)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pV(i,j-1:j,k,1))
                     pFy(i,j,k,7)=sum(pFy(i,j,k,1:2))*0.5_WP*sum(pW(i,j-1:j,k,1))
                  end if
                  ! Add pressure stress: P_face = VF*PL + (1-VF)*PG
                  pFy(i,j,k,6)=pFy(i,j,k,6)-0.5_WP*sum(pVF(i,j-1:j,k,1)*pPL(i,j-1:j,k,1)+(1.0_WP-pVF(i,j-1:j,k,1))*pPG(i,j-1:j,k,1))
                  ! Velocity gradients at y-face
                  gradU(1,1)=0.25_WP*dxi*(pU(i+1,j-1,k,1)-pU(i-1,j-1,k,1)+pU(i+1,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=dyi*(pU(i,j,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*(pU(i,j-1,k+1,1)-pU(i,j-1,k-1,1)+pU(i,j,k+1,1)-pU(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*(pV(i+1,j-1,k,1)-pV(i-1,j-1,k,1)+pV(i+1,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=dyi*(pV(i,j,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=0.25_WP*dzi*(pV(i,j-1,k+1,1)-pV(i,j-1,k-1,1)+pV(i,j,k+1,1)-pV(i,j,k-1,1))
                  gradU(1,3)=0.25_WP*dxi*(pW(i+1,j-1,k,1)-pW(i-1,j-1,k,1)+pW(i+1,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=dyi*(pW(i,j,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=0.25_WP*dzi*(pW(i,j-1,k+1,1)-pW(i,j-1,k-1,1)+pW(i,j,k+1,1)-pW(i,j,k-1,1))
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
                     ! Face velocity
                     vel=0.5_WP*sum(pW(i,j,k-1:k,1))
                     ! Rhie-Chow correction
                     irho_f=2.0_WP/max(sum(pQ(i,j,k-1:k,1:2)),this%rho_floor)
                     Pmix=pVF(i,j,k-2:k+1,1)*pPL(i,j,k-2:k+1,1)+(1.0_WP-pVF(i,j,k-2:k+1,1))*pPG(i,j,k-2:k+1,1)
                     vel=vel+0.25_WP*dt*dzi*irho_f*(Pmix(1)-3.0_WP*Pmix(0)+3.0_WP*Pmix(-1)-Pmix(-2))
                     ! WENO liquid mass and energy fluxes
                     if (any(pVF(i,j,k-1:k,1).ge.VFlo)) then
                        w=weno_weight((abs(pQ(i,j,k-1,1)-pQ(i,j,k-2,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j,k+1,1)-pQ(i,j,k  ,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j,k-2:k  ,1)) &
                        &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j,k-1:k+1,1))
                        w=weno_weight((abs(pIL(i,j,k-1,1)-pIL(i,j,k-2,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIL(i,j,k+1,1)-pIL(i,j,k  ,1))+eps)/(abs(pIL(i,j,k,1)-pIL(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,3)=0.5_WP*(pFz(i,j,k,1)-abs(pFz(i,j,k,1)))*sum(wenop*pIL(i,j,k-2:k  ,1)) &
                        &           +0.5_WP*(pFz(i,j,k,1)+abs(pFz(i,j,k,1)))*sum(wenom*pIL(i,j,k-1:k+1,1))
                     end if
                     ! WENO gas mass and energy fluxes
                     if (any(pVF(i,j,k-1:k,1).le.VFhi)) then
                        w=weno_weight((abs(pQ(i,j,k-1,2)-pQ(i,j,k-2,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j,k-1,2))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pQ(i,j,k+1,2)-pQ(i,j,k  ,2))+eps)/(abs(pQ(i,j,k,2)-pQ(i,j,k-1,2))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,2)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j,k-2:k  ,2)) &
                        &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j,k-1:k+1,2))
                        w=weno_weight((abs(pIG(i,j,k-1,1)-pIG(i,j,k-2,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pIG(i,j,k+1,1)-pIG(i,j,k  ,1))+eps)/(abs(pIG(i,j,k,1)-pIG(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                        pFz(i,j,k,4)=0.5_WP*(pFz(i,j,k,2)-abs(pFz(i,j,k,2)))*sum(wenop*pIG(i,j,k-2:k  ,1)) &
                        &           +0.5_WP*(pFz(i,j,k,2)+abs(pFz(i,j,k,2)))*sum(wenom*pIG(i,j,k-1:k+1,1))
                     end if
                     ! Momentum fluxes
                     pFz(i,j,k,5)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pU(i,j,k-1:k,1))
                     pFz(i,j,k,6)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pV(i,j,k-1:k,1))
                     pFz(i,j,k,7)=sum(pFz(i,j,k,1:2))*0.5_WP*sum(pW(i,j,k-1:k,1))
                  end if
                  ! Add pressure stress: P_face = VF*PL + (1-VF)*PG
                  pFz(i,j,k,7)=pFz(i,j,k,7)-0.5_WP*sum(pVF(i,j,k-1:k,1)*pPL(i,j,k-1:k,1)+(1.0_WP-pVF(i,j,k-1:k,1))*pPG(i,j,k-1:k,1))
                  ! Velocity gradients at z-face
                  gradU(1,1)=0.25_WP*dxi*(pU(i+1,j,k-1,1)-pU(i-1,j,k-1,1)+pU(i+1,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*(pU(i,j+1,k-1,1)-pU(i,j-1,k-1,1)+pU(i,j+1,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=dzi*(pU(i,j,k,1)-pU(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*(pV(i+1,j,k-1,1)-pV(i-1,j,k-1,1)+pV(i+1,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=0.25_WP*dyi*(pV(i,j+1,k-1,1)-pV(i,j-1,k-1,1)+pV(i,j+1,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=dzi*(pV(i,j,k,1)-pV(i,j,k-1,1))
                  gradU(1,3)=0.25_WP*dxi*(pW(i+1,j,k-1,1)-pW(i-1,j,k-1,1)+pW(i+1,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=0.25_WP*dyi*(pW(i,j+1,k-1,1)-pW(i,j-1,k-1,1)+pW(i,j+1,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=dzi*(pW(i,j,k,1)-pW(i,j,k-1,1))
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
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPL,pPG,pQ,pVisc,pBeta
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
               pPL  =>this%PL%mf(lvl)%dataptr(mfi)
               pPG  =>this%PG%mf(lvl)%dataptr(mfi)
               pQ   =>this%Q%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               ! Extra pointers at finest level
               if (lvl.eq.this%amr%maxlvl) then
                  pBand   =>band%dataptr(mfi)
                  pVx     =>Vx%dataptr(mfi)
                  pVy     =>Vy%dataptr(mfi)
                  pVz     =>Vz%dataptr(mfi)
                  pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
                  pCL     =>this%CL%dataptr(mfi)
                  pCG     =>this%CG%dataptr(mfi)
                  pCLold  =>this%CLold%dataptr(mfi)
                  pCGold  =>this%CGold%dataptr(mfi)
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
                  gradU(1,1)=0.5_WP*dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))
                  gradU(2,1)=0.5_WP*dyi*(pU(i,j+1,k,1)-pU(i,j-1,k,1))
                  gradU(3,1)=0.5_WP*dzi*(pU(i,j,k+1,1)-pU(i,j,k-1,1))
                  gradU(1,2)=0.5_WP*dxi*(pV(i+1,j,k,1)-pV(i-1,j,k,1))
                  gradU(2,2)=0.5_WP*dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))
                  gradU(3,2)=0.5_WP*dzi*(pV(i,j,k+1,1)-pV(i,j,k-1,1))
                  gradU(1,3)=0.5_WP*dxi*(pW(i+1,j,k,1)-pW(i-1,j,k,1))
                  gradU(2,3)=0.5_WP*dyi*(pW(i,j+1,k,1)-pW(i,j-1,k,1))
                  gradU(3,3)=0.5_WP*dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1))
                  div=gradU(1,1)+gradU(2,2)+gradU(3,3)
                  ! Pressure dilatation: split by VF between phasic energies - discontinuous
                  rhs(i,j,k,3)=rhs(i,j,k,3)-(       pVF(i,j,k,1))*pPL(i,j,k,1)*div
                  rhs(i,j,k,4)=rhs(i,j,k,4)-(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)*div
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
      call this%fill(lvl=this%amr%maxlvl,time=time)

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

      !> Trilinear interpolation of collocated velocity - uses pU,pV,pW
      function interp_velocity(pos) result(vel)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer  :: ipc,jpc,kpc
         real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2
         ! All cell-centered
         ipc=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP)
         jpc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP)
         kpc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP)
         ! Clamp to array bounds
         !ipc=max(lbound(pU,1),min(ubound(pU,1)-1,ipc))
         !jpc=max(lbound(pU,2),min(ubound(pU,2)-1,jpc))
         !kpc=max(lbound(pU,3),min(ubound(pU,3)-1,kpc))
         ! Cell-centered weights
         wxc1=(pos(1)-(this%amr%xlo+(real(ipc,WP)+0.5_WP)*dx))*dxi
         wyc1=(pos(2)-(this%amr%ylo+(real(jpc,WP)+0.5_WP)*dy))*dyi
         wzc1=(pos(3)-(this%amr%zlo+(real(kpc,WP)+0.5_WP)*dz))*dzi
         wxc1=max(0.0_WP,min(1.0_WP,wxc1)); wxc2=1.0_WP-wxc1
         wyc1=max(0.0_WP,min(1.0_WP,wyc1)); wyc2=1.0_WP-wyc1
         wzc1=max(0.0_WP,min(1.0_WP,wzc1)); wzc2=1.0_WP-wzc1
         vel(1)=wzc1*(wyc1*(wxc1*pU(ipc+1,jpc+1,kpc+1,1)+wxc2*pU(ipc,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxc1*pU(ipc+1,jpc  ,kpc+1,1)+wxc2*pU(ipc,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxc1*pU(ipc+1,jpc+1,kpc  ,1)+wxc2*pU(ipc,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxc1*pU(ipc+1,jpc  ,kpc  ,1)+wxc2*pU(ipc,jpc  ,kpc  ,1)))
         vel(2)=wzc1*(wyc1*(wxc1*pV(ipc+1,jpc+1,kpc+1,1)+wxc2*pV(ipc,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxc1*pV(ipc+1,jpc  ,kpc+1,1)+wxc2*pV(ipc,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxc1*pV(ipc+1,jpc+1,kpc  ,1)+wxc2*pV(ipc,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxc1*pV(ipc+1,jpc  ,kpc  ,1)+wxc2*pV(ipc,jpc  ,kpc  ,1)))
         vel(3)=wzc1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpc+1,1)+wxc2*pW(ipc,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpc+1,1)+wxc2*pW(ipc,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpc  ,1)+wxc2*pW(ipc,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpc  ,1)+wxc2*pW(ipc,jpc  ,kpc  ,1)))
      end function interp_velocity

   end subroutine get_dQdt

   !> Override build_plic to include Q clean-up
   subroutine build_plic(this,time)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: time
      ! Call parent build_plic
      call this%amrvof%build_plic(time)
      ! Clean up Q
      call this%clean_Q()
   end subroutine build_plic

   !> Clean up Q in pure cells
   subroutine clean_Q(this)
      use amrex_amr_module, only: amrex_mfiter
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
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl,i,j,k
      real(WP) :: t0
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pCL,pCG
      real(WP) :: dx,dy,dz
      ! If no relaxation model was provided, return
      if (.not.associated(this%relax)) return
      ! Start timer
      t0=MPI_Wtime()
      ! Apply relaxation on finest level only (mixture cells are always at finest)
      lvl=this%amr%maxlvl
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pQ   =>this%Q%mf(lvl)%dataptr(mfi)
         pCL  =>this%CL%dataptr(mfi)
         pCG  =>this%CG%dataptr(mfi)
         ! Loop over all cells
         bx=mfi%growntilebox(this%nover)
         ! Loop over valid cells
         !bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Only relax mixture cells
            if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) cycle
            ! Apply user-provided relaxation model (modifies VF and Q)
            call this%relax(VF=pVF(i,j,k,1),Q=pQ(i,j,k,:))
            ! Ensure consistency with modified VF
            if (pVF(i,j,k,1).lt.VFlo) then
               ! Pure liquid
               pVF(i,j,k,1)=0.0_WP
               pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pQ(i,j,k,1)=0.0_WP
               pQ(i,j,k,3)=0.0_WP
            else if (pVF(i,j,k,1).gt.VFhi) then
               ! Pure gas
               pVF(i,j,k,1)=1.0_WP
               pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pQ(i,j,k,2)=0.0_WP
               pQ(i,j,k,4)=0.0_WP
            end if
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Sync and apply BC
      !call this%fill(lvl=this%amr%maxlvl,time=time)
      !call this%Q%fill(time=time)
      ! End timer
      this%wt_relax=this%wt_relax+(MPI_Wtime()-t0)
   end subroutine apply_relax

   !> Apply relaxation to mixture cells
   ! subroutine apply_relax(this,time)
   !    use amrex_amr_module, only: amrex_mfiter,amrex_box
   !    use mpi_f08, only: MPI_Wtime
   !    implicit none
   !    class(amrmpcomp), intent(inout) :: this
   !    real(WP), intent(in) :: time
   !    integer :: lvl,i,j,k,n
   !    real(WP) :: t0,Peq,dVF,Peq_avg,VF_avg,dVF_avg
   !    real(WP), dimension(7) :: dQ,Q_avg,dQ_avg
   !    type(amrex_mfiter) :: mfi
   !    type(amrex_box) :: bx
   !    real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pCL,pCG
   !    real(WP) :: dx,dy,dz
   !    ! If no relaxation model was provided, return
   !    if (.not.associated(this%relax)) return
   !    ! Start timer
   !    t0=MPI_Wtime()
   !    ! Apply relaxation on finest level only (mixture cells are always at finest)
   !    lvl=this%amr%maxlvl
   !    dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
   !    call this%amr%mfiter_build(lvl,mfi)
   !    do while (mfi%next())
   !       ! Get pointers to data
   !       pVF  =>this%VF%mf(lvl)%dataptr(mfi)
   !       pQ   =>this%Q%mf(lvl)%dataptr(mfi)
   !       pCL  =>this%CL%dataptr(mfi)
   !       pCG  =>this%CG%dataptr(mfi)
   !       ! Loop over all cells
   !       bx=mfi%growntilebox(this%nover)
   !       ! Loop over valid cells
   !       !bx=mfi%tilebox()
   !       do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
   !          ! Only relax mixture cells
   !          if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) cycle
   !          ! Apply relaxation

   !          ! Get equilibrium pressure
   !          Peq=-huge(1.0_WP); call this%relax(VF=pVF(i,j,k,1),Q=pQ(i,j,k,:),Peq=Peq)
   !          if (Peq.eq.-huge(1.0_WP)) cycle
   !          ! Get increments
   !          call this%relax(VF=pVF(i,j,k,1),Q=pQ(i,j,k,:),Peq=Peq,dVF=dVF,dQ=dQ)
   !          ! Diagnostic: check supercell feasibility
   !          VF_avg=sum(pVF(i-1:i+1,j-1:j+1,k-1:k+1,1))/27.0_WP
   !          do n=1,7; Q_avg(n)=sum(pQ(i-1:i+1,j-1:j+1,k-1:k+1,n))/27.0_WP; end do
   !          Peq_avg=-huge(1.0_WP); call this%relax(VF=VF_avg,Q=Q_avg,Peq=Peq_avg)
   !          if (Peq_avg.ne.-huge(1.0_WP)) then
   !             call this%relax(VF=pVF(i,j,k,1),Q=pQ(i,j,k,:),Peq=Peq_avg,dVF=dVF_avg,dQ=dQ_avg)
   !             if (pQ(i,j,k,3)+dQ_avg(3).le.0.0_WP.or.pQ(i,j,k,4)+dQ_avg(4).le.0.0_WP.or.pVF(i,j,k,1)+dVF_avg.le.0.0_WP.or.pVF(i,j,k,1)+dVF_avg.ge.1.0_WP) then
   !                print '(A,3I5,6(A,ES10.3),4(A,L1))', ' SC ',i,j,k, &
   !                & ' VF=',pVF(i,j,k,1),' Peq=',Peq,' Pavg=',Peq_avg, &
   !                & ' dVF=',dVF_avg,' Q3n=',pQ(i,j,k,3)+dQ_avg(3),' Q4n=',pQ(i,j,k,4)+dQ_avg(4), &
   !                & ' liqE=',pQ(i,j,k,3)+dQ_avg(3).gt.0.0_WP, &
   !                & ' gasE=',pQ(i,j,k,4)+dQ_avg(4).gt.0.0_WP, &
   !                & ' liqV=',pVF(i,j,k,1)+dVF_avg.gt.0.0_WP, &
   !                & ' gasV=',pVF(i,j,k,1)+dVF_avg.lt.1.0_WP
   !             end if
   !          end if
   !          ! Apply increments
   !          pVF(i,j,k,1)=pVF(i,j,k,1)+dVF
   !          pQ(i,j,k,3)=pQ(i,j,k,3)+dQ(3)
   !          pQ(i,j,k,4)=pQ(i,j,k,4)+dQ(4)
   !          ! Ensure consistency with modified VF
   !          if (pVF(i,j,k,1).lt.VFlo) then
   !             ! Pure liquid
   !             pVF(i,j,k,1)=0.0_WP
   !             pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
   !             pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
   !             pQ(i,j,k,1)=0.0_WP
   !             pQ(i,j,k,3)=0.0_WP
   !          else if (pVF(i,j,k,1).gt.VFhi) then
   !             ! Pure gas
   !             pVF(i,j,k,1)=1.0_WP
   !             pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
   !             pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
   !             pQ(i,j,k,2)=0.0_WP
   !             pQ(i,j,k,4)=0.0_WP
   !          end if
   !       end do; end do; end do
   !    end do
   !    call this%amr%mfiter_destroy(mfi)
   !    ! Sync and apply BC
   !    call this%fill(lvl=this%amr%maxlvl,time=time)
   !    call this%Q%fill(time=time)
   !    ! End timer
   !    this%wt_relax=this%wt_relax+(MPI_Wtime()-t0)
   ! end subroutine apply_relax
   
   !> Add artificial bulk viscosity to this%beta
   subroutine add_viscartif(this,dt,Cartif)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cartif
      ! Local variables
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: beta_t,scratch
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBeta_t,pScratch,pU,pV,pW,pC,pBeta,pVF,pRHOL,pRHOG
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_beta,myCartif,t0
      real(WP) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy,vort,grad_div
      integer :: lvl,i,j,k
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2
      
      ! Start timer
      t0=MPI_Wtime()

      ! Set model constant
      if (present(Cartif)) then; myCartif=Cartif; else; myCartif=5.0_WP; end if
      
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         
         ! Grid spacings
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz
         
         ! Max beta from CFL
         max_beta=max_cfl*this%amr%min_meshsize(lvl)**2/(4.0_WP*dt)
         
         ! Build temp multifabs
         call this%amr%mfab_build(lvl=lvl,mfab=scratch,ncomp=1,nover=1); call scratch%setval(0.0_WP)
         call this%amr%mfab_build(lvl=lvl,mfab=beta_t,ncomp=1,nover=this%nover); call beta_t%setval(0.0_WP)
         
         ! Phase 1: Compute divergence into scratch
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pScratch=>scratch%dataptr(mfi)
            bx=mfi%growntilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pScratch(i,j,k,1)=0.5_WP*(dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1)))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         
         ! Phase 2: Compute kinematic beta
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pC=>this%C%mf(lvl)%dataptr(mfi)
            pScratch=>scratch%dataptr(mfi)
            pBeta_t=>beta_t%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Only work in compression regions
               if (pScratch(i,j,k,1).ge.0.0_WP) cycle
               ! Compute local vorticity
               dudy=0.5_WP*dyi*(pU(i,j+1,k,1)-pU(i,j-1,k,1))
               dudz=0.5_WP*dzi*(pU(i,j,k+1,1)-pU(i,j,k-1,1))
               dvdx=0.5_WP*dxi*(pV(i+1,j,k,1)-pV(i-1,j,k,1))
               dvdz=0.5_WP*dzi*(pV(i,j,k+1,1)-pV(i,j,k-1,1))
               dwdx=0.5_WP*dxi*(pW(i+1,j,k,1)-pW(i-1,j,k,1))
               dwdy=0.5_WP*dyi*(pW(i,j+1,k,1)-pW(i,j-1,k,1))
               vort=(dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2
               ! Compute |grad(div)|
               grad_div=max(abs(pScratch(i+1,j,k,1)-pScratch(i,j,k,1)),abs(pScratch(i,j,k,1)-pScratch(i-1,j,k,1)))*dx**2 &
               &       +max(abs(pScratch(i,j+1,k,1)-pScratch(i,j,k,1)),abs(pScratch(i,j,k,1)-pScratch(i,j-1,k,1)))*dy**2 &
               &       +max(abs(pScratch(i,j,k+1,1)-pScratch(i,j,k,1)),abs(pScratch(i,j,k,1)-pScratch(i,j,k-1,1)))*dz**2
               ! Floor vorticity with sound speed
               vort=max(vort,(0.05_WP*pC(i,j,k,1)/this%amr%min_meshsize(lvl))**2)
               ! Compute beta
               pBeta_t(i,j,k,1)=myCartif*grad_div*min(4.0_WP/3.0_WP*pScratch(i,j,k,1)**2/(pScratch(i,j,k,1)**2+vort+1.0e-15_WP),1.0_WP)
               ! Clip to max
               pBeta_t(i,j,k,1)=min(pBeta_t(i,j,k,1),max_beta)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy scratch
         call amrex_multifab_destroy(scratch)

         ! Phase 3: Filter beta_t
         call this%amr%mfab_filter(lvl=lvl,mfab=beta_t,npass=nfilter)

         ! Phase 4: Convert to dynamic viscosity via harmonic averaging and add to this%beta
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pBeta_t=>beta_t%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBeta(i,j,k,1)=pBeta(i,j,k,1)+pBeta_t(i,j,k,1)/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),this%rho_floor)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),this%rho_floor))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy temp multifab
         call amrex_multifab_destroy(beta_t)
         
      end do

      ! End timer
      this%wt_visc=this%wt_visc+(MPI_Wtime()-t0)

   end subroutine add_viscartif

   !> Add Vreman SGS eddy viscosity to this%visc
   subroutine add_vreman(this,dt,Cs)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      ! Local variables
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: visc_t
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pU,pV,pW,pVisc,pVF,pRHOL,pRHOG
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,Cmodel,Aij,Bij,t0
      real(WP), dimension(1:3,1:3) :: gradU,betaij
      integer :: lvl,i,j,k,si,sj
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2
      
      ! Start timer
      t0=MPI_Wtime()

      ! Model constant: c=2.5*Cs**2 (Vreman uses c=0.07 which corresponds to Cs=0.17)
      if (present(Cs)) then; Cmodel=2.5_WP*Cs**2; else; Cmodel=2.5_WP*0.17_WP**2; end if
      
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         
         ! Grid spacings
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz
         
         ! Max visc from CFL
         max_visc=max_cfl*this%amr%min_meshsize(lvl)**2/(4.0_WP*dt)
         
         ! Build temp multifab for eddy viscosity at this level (nover ghost cells)
         call this%amr%mfab_build(lvl=lvl,mfab=visc_t,ncomp=1,nover=this%nover); call visc_t%setval(0.0_WP)
         
         ! Phase 1: Compute kinematic eddy viscosity
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get data pointers
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pVisc_t=>visc_t%dataptr(mfi)
            ! Loop over interior tiles
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute cell-centered velocity gradient tensor
               gradU(1,1)=0.5_WP*dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))
               gradU(2,1)=0.5_WP*dyi*(pU(i,j+1,k,1)-pU(i,j-1,k,1))
               gradU(3,1)=0.5_WP*dzi*(pU(i,j,k+1,1)-pU(i,j,k-1,1))
               gradU(1,2)=0.5_WP*dxi*(pV(i+1,j,k,1)-pV(i-1,j,k,1))
               gradU(2,2)=0.5_WP*dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))
               gradU(3,2)=0.5_WP*dzi*(pV(i,j,k+1,1)-pV(i,j,k-1,1))
               gradU(1,3)=0.5_WP*dxi*(pW(i+1,j,k,1)-pW(i-1,j,k,1))
               gradU(2,3)=0.5_WP*dyi*(pW(i,j+1,k,1)-pW(i,j-1,k,1))
               gradU(3,3)=0.5_WP*dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1))
               ! Compute A=gradU_ij*gradU_ij invariant
               Aij=sum(gradU**2)
               ! Compute beta_ij=dx_m^2*gradU_mi*gradU_mj
               do sj=1,3; do si=1,3
                  betaij(si,sj)=dx**2*gradU(1,si)*gradU(1,sj)+dy**2*gradU(2,si)*gradU(2,sj)+dz**2*gradU(3,si)*gradU(3,sj)
               end do; end do
               ! Compute B invariant
               Bij=betaij(1,1)*betaij(2,2)-betaij(1,2)**2+betaij(1,1)*betaij(3,3)-betaij(1,3)**2+betaij(2,2)*betaij(3,3)-betaij(2,3)**2
               ! Assemble eddy viscosity
               if (Bij.gt.0.0_WP) then
                  pVisc_t(i,j,k,1)=Cmodel*sqrt(Bij/Aij)
               end if
               ! Clip to CFL limit
               pVisc_t(i,j,k,1)=min(pVisc_t(i,j,k,1),max_visc)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Phase 2: Filter visc_t
         call this%amr%mfab_filter(lvl=lvl,mfab=visc_t,npass=nfilter)

         ! Phase 3: Convert to dynamic viscosity via harmonic averaging and add to this%visc
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%dataptr(mfi)
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

         ! Destroy temp multifab
         call amrex_multifab_destroy(visc_t)
         
      end do

      ! End timer
      this%wt_visc=this%wt_visc+(MPI_Wtime()-t0)

   end subroutine add_vreman

   !> Calculate CFL numbers
   subroutine get_cfl(this,dt,cfl)
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      integer :: lvl
      real(WP) :: Umax,Vmax,Wmax,Cmax,viscmax
      ! Reset CFLs
      this%CFLc_x=0.0_WP; this%CFLc_y=0.0_WP; this%CFLc_z=0.0_WP
      this%CFLa_x=0.0_WP; this%CFLa_y=0.0_WP; this%CFLa_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         ! Max velocity
         Umax=this%U%norm0(lvl=lvl)
         Vmax=this%V%norm0(lvl=lvl)
         Wmax=this%W%norm0(lvl=lvl)
         ! Max speed of sound
         Cmax=this%C%norm0(lvl=lvl)
         ! Max viscosities (TODO: use phasic viscosities viscL/RHOL, viscG/RHOG instead of mixture)
         get_viscmax: block
            use amrex_amr_module, only: amrex_mfiter,amrex_box
            use parallel, only: MPI_REAL_WP
            use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
            type(amrex_mfiter) :: mfi
            type(amrex_box) :: bx
            real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVisc,pBeta,pDiff
            integer :: i,j,k,ierr
            real(WP) :: rho
            viscmax=0.0_WP
            call this%amr%mfiter_build(lvl,mfi)
            do while(mfi%next())
               ! Get data pointers
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pDiff=>this%diff%mf(lvl)%dataptr(mfi)
               ! Loop over interior tiles
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
                  viscmax=max(viscmax,pVisc(i,j,k,1)/rho,pBeta(i,j,k,1)/rho,pDiff(i,j,k,1)/rho)   ! This is incorrect for heat diffusion!
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            call MPI_ALLREDUCE(MPI_IN_PLACE,viscmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         end block get_viscmax
         ! Convective+acoustic
         if (this%amr%nx.gt.1) this%CFLc_x=max(this%CFLc_x,(Umax+Cmax)*dt/this%amr%dx(lvl))
         if (this%amr%ny.gt.1) this%CFLc_y=max(this%CFLc_y,(Vmax+Cmax)*dt/this%amr%dy(lvl))
         if (this%amr%nz.gt.1) this%CFLc_z=max(this%CFLc_z,(Wmax+Cmax)*dt/this%amr%dz(lvl))
         ! Acoustic
         if (this%amr%nx.gt.1) this%CFLa_x=max(this%CFLa_x,Cmax*dt/this%amr%dx(lvl))
         if (this%amr%ny.gt.1) this%CFLa_y=max(this%CFLa_y,Cmax*dt/this%amr%dy(lvl))
         if (this%amr%nz.gt.1) this%CFLa_z=max(this%CFLa_z,Cmax*dt/this%amr%dz(lvl))
         ! Viscous
         if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt/this%amr%dx(lvl)**2)
         if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt/this%amr%dy(lvl)**2)
         if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt/this%amr%dz(lvl)**2)
      end do
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLa_x,this%CFLa_y,this%CFLa_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
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

      ! Get VF info from parent
      call this%amrvof%get_info()

      ! VF-conditional phasic extrema
      phasic_extrema: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG
         integer :: i,j,k
         ! Initialize extrema
         this%Umax=0.0_WP; this%Vmax=0.0_WP; this%Wmax=0.0_WP
         this%RHOLmin=huge(1.0_WP); this%RHOLmax=-huge(1.0_WP); this%RHOGmin=huge(1.0_WP); this%RHOGmax=-huge(1.0_WP)
         this%ILmin=huge(1.0_WP); this%ILmax=-huge(1.0_WP); this%IGmin=huge(1.0_WP); this%IGmax=-huge(1.0_WP)
         this%PLmin=huge(1.0_WP); this%PLmax=-huge(1.0_WP); this%PGmin=huge(1.0_WP); this%PGmax=-huge(1.0_WP)
         this%TLmin=huge(1.0_WP); this%TLmax=-huge(1.0_WP); this%TGmin=huge(1.0_WP); this%TGmax=-huge(1.0_WP)
         this%Cmin=huge(1.0_WP); this%Cmax=-huge(1.0_WP)
         this%dPmax=0.0_WP
         this%Qmin=huge(1.0_WP); this%Qmax=-huge(1.0_WP)
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Velocity norm 0
            this%Umax=max(this%Umax,this%U%norm0(lvl=lvl))
            this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl))
            this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl))
            ! Extrema of mixture speed of sound
            this%Cmin=min(this%Cmin,this%C%get_min(lvl=lvl)); this%Cmax=max(this%Cmax,this%C%get_max(lvl=lvl))
            ! Extrema of conserved variables
            do n=1,this%Q%ncomp
               this%Qmin(n)=min(this%Qmin(n),this%Q%get_min(lvl=lvl,comp=n))
               this%Qmax(n)=max(this%Qmax(n),this%Q%get_max(lvl=lvl,comp=n))
            end do
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

      ! Conserved integrals at base level
      dV=this%amr%cell_vol(0)
      do n=1,this%Q%ncomp
         this%Qint(n)=this%Q%get_sum(lvl=0,comp=n)*dV
      end do

      ! Kinetic energy integral: 0.5 * rho * (U^2 + V^2 + W^2) * dV
      ! Uses composite integration with fine masking to avoid double-counting
      get_rhoKint: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW
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
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               ! Get pointer to fine mask
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over cells
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
                  ! Accumulate kinetic energy (rho = Q(1)+Q(2) = VF*rhoL + (1-VF)*rhoG)
                  this%rhoKint=this%rhoKint+0.5_WP*(pQ(i,j,k,1)+pQ(i,j,k,2))*(pU(i,j,k,1)**2+pV(i,j,k,1)**2+pW(i,j,k,1)**2)*dV
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         ! Reduce across MPI ranks
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_rhoKint

      ! Reduce per-rank timing to min/max across ranks
      call MPI_ALLREDUCE(this%wt_prim,      this%wtmax_prim,      1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_prim,      this%wtmin_prim,      1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_dQdt,      this%wtmax_dQdt,      1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_dQdt,      this%wtmin_dQdt,      1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_fv,        this%wtmax_fv,        1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_fv,        this%wtmin_fv,        1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_div,       this%wtmax_div,       1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_div,       this%wtmin_div,       1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_relax,     this%wtmax_relax,     1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_relax,     this%wtmin_relax,     1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc,      this%wtmax_visc,      1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc,      this%wtmin_visc,      1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
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
      ! VOF data is registered via parent
      call this%amrvof%register_checkpoint(io)
      ! Register flow data
      call io%add_data(this%Q,'Q')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      ! VOF data is restored via parent
      call this%amrvof%restore_checkpoint(io,dirname,time)
      ! Restore flow data
      call io%read_data(dirname,this%Q,'Q')
      ! Fill ghost cells as io reads valid data only
      call this%Q%fill(time=time)
      ! Rebuild primitive variables
      call this%get_primitive(this%Q)
   end subroutine restore_checkpoint

end module amrmpcomp_class
