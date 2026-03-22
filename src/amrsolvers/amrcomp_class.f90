!> AMR Compressible solver class
!> Port of fastcomp_class to AMReX infrastructure
module amrcomp_class
   use iso_c_binding,    only: c_ptr
   use precision,        only: WP
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrsolver_class,  only: amrsolver
   use amrex_amr_module, only: amrex_boxarray,amrex_distromap
   implicit none
   private

   ! Expose type
   public :: amrcomp

   !> AMR Compressible solver type
   type, extends(amrsolver) :: amrcomp
      ! User-configurable callbacks
      procedure(comp_init_iface), pointer, pass :: user_init=>null()
      procedure(comp_tagging_iface), pointer, pass :: user_tagging=>null()
      procedure(comp_bc_iface), pointer, pass :: user_bc=>null()

      ! Equation of state function pointers: P=P(rho,I), C=C(rho,P), T=T(rho,P)
      procedure(eos_P_iface), pointer, nopass :: getP=>null()
      procedure(eos_C_iface), pointer, nopass :: getC=>null()
      procedure(eos_T_iface), pointer, nopass :: getT=>null()

      ! Conserved variables (1: rho, 2: rhoU, 3: rhoV, 4: rhoW, 5: rhoI)
      type(amrdata) :: Q,Qold

      ! Primitive variables (velocities, internal energy, and pressure)
      type(amrdata) :: U,V,W,I,P

      ! Temperature
      type(amrdata) :: T

      ! Speed of sound
      type(amrdata) :: C

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
      real(WP) :: Imin=0.0_WP,Imax=0.0_WP
      real(WP) :: Pmin=0.0_WP,Pmax=0.0_WP
      real(WP) :: Tmin=0.0_WP,Tmax=0.0_WP
      real(WP) :: Cmin=0.0_WP,Cmax=0.0_WP
      real(WP), dimension(5) :: Qint=0.0_WP,Qmin=0.0_WP,Qmax=0.0_WP
      real(WP) :: rhoKint=0.0_WP

      ! Minimum density for stability
      real(WP) :: rho_floor=1.0e-10_WP

      ! Number of overlap cells (2 for WENO3 stencil)
      integer :: nover=2

   contains
      procedure :: initialize
      procedure :: finalize
      ! Lifecycle callbacks
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      ! Physics
      procedure :: get_primitive
      procedure :: get_conserved
      procedure :: get_dQdt
      procedure :: add_viscartif
      procedure :: add_vreman
      procedure :: get_cfl
      ! Print
      procedure :: get_info
      procedure :: print=>amrcomp_print
      ! Checkpoint I/O (deferred from amrsolver)
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrcomp

   !> Abstract interface for user init callback
   abstract interface
      subroutine comp_init_iface(solver,lvl,time,ba,dm)
         import :: amrcomp,WP,amrex_boxarray,amrex_distromap
         class(amrcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine comp_init_iface
   end interface

   !> Abstract interface for tagging callback
   abstract interface
      subroutine comp_tagging_iface(solver,lvl,time,tags)
         import :: amrcomp,c_ptr,WP
         class(amrcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine comp_tagging_iface
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

   !> Abstract interface for user BC callback
   abstract interface
      subroutine comp_bc_iface(solver,lvl,time,face,bx,pQ)
         use amrex_amr_module, only: amrex_box
         import :: amrcomp,WP
         class(amrcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      end subroutine comp_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrcomp_on_init(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(lvl,time,ba,dm)
   end subroutine amrcomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrcomp_on_coarse(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrcomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrcomp_on_remake(ctx,lvl,time,ba,dm)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrcomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrcomp_on_clear(ctx,lvl)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrcomp_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrcomp_tagging(ctx,lvl,time,tags)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      if (associated(this%user_tagging)) call this%user_tagging(lvl,time,tags)
   end subroutine amrcomp_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrcomp_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrcomp_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the compressible solver
   subroutine initialize(this,amr,name)
      use iso_c_binding, only: c_loc
      use amrex_amr_module, only: amrex_bc_foextrap
      implicit none
      class(amrcomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Set name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'UNNAMED_COMP'
      end if
      this%amr => amr

      ! Initialize conserved variables Q
      call this%Q%initialize   (amr,name='Q'   ,ncomp=5,ng=this%nover); this%Q%parent   =>this
      call this%Qold%initialize(amr,name='Qold',ncomp=5,ng=this%nover); this%Qold%parent=>this

      ! Initialize primitive/derived variables
      call this%U%initialize(amr,name='U',ncomp=1,ng=this%nover); this%U%parent=>this
      call this%V%initialize(amr,name='V',ncomp=1,ng=this%nover); this%V%parent=>this
      call this%W%initialize(amr,name='W',ncomp=1,ng=this%nover); this%W%parent=>this
      call this%I%initialize(amr,name='I',ncomp=1,ng=this%nover); this%I%parent=>this
      call this%P%initialize(amr,name='P',ncomp=1,ng=this%nover); this%P%parent=>this
      call this%T%initialize(amr,name='T',ncomp=1,ng=this%nover); this%T%parent=>this
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
       type is (amrcomp)
         call this%amr%add_on_init   (amrcomp_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrcomp_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrcomp_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrcomp_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrcomp_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrcomp_postregrid,c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the compressible solver
   subroutine finalize(this)
      implicit none
      class(amrcomp), intent(inout) :: this
      call this%Q%finalize(); call this%Qold%finalize()
      call this%U%finalize(); call this%V%finalize(); call this%W%finalize()
      call this%I%finalize(); call this%P%finalize(); call this%T%finalize(); call this%C%finalize()
      call this%visc%finalize(); call this%beta%finalize(); call this%diff%finalize()
      nullify(this%amr); nullify(this%user_init); nullify(this%user_tagging)
      nullify(this%getP); nullify(this%getC); nullify(this%getT)
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%Q%reset_level(lvl,ba,dm); call this%Qold%reset_level(lvl,ba,dm)
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%I%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%T%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      ! Zero out
      call this%Q%setval(val=0.0_WP,lvl=lvl); call this%Qold%setval(val=0.0_WP,lvl=lvl)
      call this%U%setval(val=0.0_WP,lvl=lvl)
      call this%V%setval(val=0.0_WP,lvl=lvl)
      call this%W%setval(val=0.0_WP,lvl=lvl)
      call this%I%setval(val=0.0_WP,lvl=lvl)
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%T%setval(val=0.0_WP,lvl=lvl)
      call this%C%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl)
      call this%beta%setval(val=0.0_WP,lvl=lvl)
      call this%diff%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using conservative interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Use conservative interpolation from coarse for Q
      call this%Q%on_coarse(lvl,time,ba,dm)
      ! Qold just needs to be reset
      call this%Qold%reset_level(lvl,ba,dm)
      ! Derived variables are just reset
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%I%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%T%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using conservative interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Use on_remake from amrdata for Q
      call this%Q%on_remake(lvl,time,ba,dm)
      ! Qold just needs to be reset
      call this%Qold%reset_level(lvl,ba,dm)
      ! Derived variables are just reset
      call this%U%reset_level(lvl,ba,dm)
      call this%V%reset_level(lvl,ba,dm)
      call this%W%reset_level(lvl,ba,dm)
      call this%I%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%T%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%Q%clear_level(lvl)
      call this%Qold%clear_level(lvl)
      call this%U%clear_level(lvl)
      call this%V%clear_level(lvl)
      call this%W%clear_level(lvl)
      call this%I%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%T%clear_level(lvl)
      call this%C%clear_level(lvl)
      call this%visc%clear_level(lvl)
      call this%beta%clear_level(lvl)
      call this%diff%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      integer :: lvl
      ! Average down conserved variables for C/F consistency
      do lvl=this%amr%clvl()-1,lbase,-1
         call this%Q%average_downto(lvl)
      end do
      ! Fill ghosts
      call this%Q%fill(time)
      ! Rebuild primitive variables
      call this%get_primitive(this%Q)
   end subroutine post_regrid

   !> Internal fillbc for Q - calls default_fillbc first, then user_bc for ext_dir faces
   subroutine Q_fillbc(this,mf,scomp,ncomp,time,geom)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,&
      &                           amrex_box,amrex_geometry,amrex_multifab,amrex_bc_ext_dir
      use amrdata_class, only: default_fillbc
      implicit none
      class(amrdata), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: mf
      integer, intent(in) :: scomp, ncomp
      real(WP), intent(in) :: time
      type(amrex_geometry), intent(in) :: geom
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bc_bx
      class(amrcomp), pointer :: solver
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: dlo(3),dhi(3)
      integer :: lvl
      
      ! First apply default BC handling (foextrap,hoextrap,reflect,etc.)
      call default_fillbc(this,mf,scomp,ncomp,time,geom)
      
      ! Access parent solver
      select type (s=>this%parent)
       class is (amrcomp)
         solver=>s
      end select
      
      ! Check if user callback exists
      if (.not.associated(solver%user_bc)) return

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
         if (this%lo_bc(1,1).eq.amrex_bc_ext_dir .and. ilo.lt.dlo(1)) then
            bc_bx=amrex_box([ilo,jlo,klo],[dlo(1)-1,jhi,khi])
            call solver%user_bc(lvl=lvl,time=time,face=1,bx=bc_bx,pQ=p)
         end if

         ! X-HIGH (face=2)
         if (this%hi_bc(1,1).eq.amrex_bc_ext_dir .and. ihi.gt.dhi(1)) then
            bc_bx=amrex_box([dhi(1)+1,jlo,klo],[ihi,jhi,khi])
            call solver%user_bc(lvl=lvl,time=time,face=2,bx=bc_bx,pQ=p)
         end if

         ! Y-LOW (face=3)
         if (this%lo_bc(2,1).eq.amrex_bc_ext_dir .and. jlo.lt.dlo(2)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,dlo(2)-1,khi])
            call solver%user_bc(lvl=lvl,time=time,face=3,bx=bc_bx,pQ=p)
         end if

         ! Y-HIGH (face=4)
         if (this%hi_bc(2,1).eq.amrex_bc_ext_dir .and. jhi.gt.dhi(2)) then
            bc_bx=amrex_box([ilo,dhi(2)+1,klo],[ihi,jhi,khi])
            call solver%user_bc(lvl=lvl,time=time,face=4,bx=bc_bx,pQ=p)
         end if

         ! Z-LOW (face=5)
         if (this%lo_bc(3,1).eq.amrex_bc_ext_dir .and. klo.lt.dlo(3)) then
            bc_bx=amrex_box([ilo,jlo,klo],[ihi,jhi,dlo(3)-1])
            call solver%user_bc(lvl=lvl,time=time,face=5,bx=bc_bx,pQ=p)
         end if

         ! Z-HIGH (face=6)
         if (this%hi_bc(3,1).eq.amrex_bc_ext_dir .and. khi.gt.dhi(3)) then
            bc_bx=amrex_box([ilo,jlo,dhi(3)+1],[ihi,jhi,khi])
            call solver%user_bc(lvl=lvl,time=time,face=6,bx=bc_bx,pQ=p)
         end if

      end do
      call amrex_mfiter_destroy(mfi)

   end subroutine Q_fillbc

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Calculate primitive variables from conserved variables
   subroutine get_primitive(this,Q)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use messager, only: die
      implicit none
      class(amrcomp), intent(inout) :: this
      type(amrdata), intent(in) :: Q
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pI,pP,pT,pC
      real(WP) :: rho_inv
      ! Check passed Q is as expected
      if (Q%ncomp.ne.5) call die('[amrcomp get_primitive] Q must have 5 components')
      if (Q%ng.lt.this%nover) call die('[amrcomp get_primitive] Q must have at least nover ghost cells')
      ! Check EoS functions are set
      if (.not.associated(this%getP)) call die('[amrcomp get_primitive] getP not set')
      if (.not.associated(this%getC)) call die('[amrcomp get_primitive] getC not set')
      if (.not.associated(this%getT)) call die('[amrcomp get_primitive] getT not set')
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>Q%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pI=>this%I%mf(lvl)%dataptr(mfi)
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pT=>this%T%mf(lvl)%dataptr(mfi)
            pC=>this%C%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute velocity from momentum
               rho_inv=1.0_WP/max(pQ(i,j,k,1),this%rho_floor)
               pU(i,j,k,1)=pQ(i,j,k,2)*rho_inv
               pV(i,j,k,1)=pQ(i,j,k,3)*rho_inv
               pW(i,j,k,1)=pQ(i,j,k,4)*rho_inv
               ! Compute internal energy per unit mass
               pI(i,j,k,1)=pQ(i,j,k,5)*rho_inv
               ! Compute pressure via EoS: P = P(rho, I)
               pP(i,j,k,1)=this%getP(rho=pQ(i,j,k,1),I=pI(i,j,k,1))
               ! Compute speed of sound via EoS: C = C(rho, P)
               pC(i,j,k,1)=this%getC(rho=pQ(i,j,k,1),P=pP(i,j,k,1))
               ! Compute temperature via EoS: T = T(rho, P)
               pT(i,j,k,1)=this%getT(rho=pQ(i,j,k,1),P=pP(i,j,k,1))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_primitive

   !> Calculate conserved variables from primitive variables
   subroutine get_conserved(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      class(amrcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pI
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%growntilebox(this%nover)
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pI=>this%I%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pQ(i,j,k,2)=pQ(i,j,k,1)*pU(i,j,k,1)
               pQ(i,j,k,3)=pQ(i,j,k,1)*pV(i,j,k,1)
               pQ(i,j,k,4)=pQ(i,j,k,1)*pW(i,j,k,1)
               pQ(i,j,k,5)=pQ(i,j,k,1)*pI(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_conserved

   !> Calculate dQdt from passed Q
   subroutine get_dQdt(this,Q,dQdt,time)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_box
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: Q
      type(amrdata), intent(inout) :: dQdt
      real(WP), intent(in) :: time
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pP,pI,rhs,pFx,pFy,pFz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc,pBeta,pDiff,pT
      real(WP), dimension(-2: 0) :: wenop
      real(WP), dimension(-1:+1) :: wenom
      real(WP), dimension(1:3,1:3) :: gradU
      real(WP) :: w,dxi,dyi,dzi,div,vel
      real(WP), parameter :: eps=1.0e-15_WP
      integer :: lvl,i,j,k

      ! First build primitive variables from Q
      call this%get_primitive(Q)
      
      ! Phase 1: Compute fluxes for all levels
      do lvl=0,this%amr%clvl()
         
         ! Grid spacings for this level
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         
         ! Build face-centered flux MultiFabs for this level
         call this%amr%mfab_build(lvl=lvl,mfab=Fx(lvl),ncomp=5,nover=0,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl=lvl,mfab=Fy(lvl),ncomp=5,nover=0,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl=lvl,mfab=Fz(lvl),ncomp=5,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
         
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())

            ! Get data pointers
            pQ=>Q%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pI=>this%I%mf(lvl)%dataptr(mfi)
            pT=>this%T%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pDiff=>this%diff%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi)
            pFy=>Fy(lvl)%dataptr(mfi)
            pFz=>Fz(lvl)%dataptr(mfi)
            
            ! X-fluxes
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Face velocity
               vel=0.5_WP*sum(pU(i-1:i,j,k,1))
               ! WENO mass flux
               w=weno_weight((abs(pQ(i-1,j,k,1)-pQ(i-2,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(pQ(i+1,j,k,1)-pQ(i  ,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
               pFx(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i-2:i  ,j,k,1)) &
               &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i-1:i+1,j,k,1))
               ! Momentum fluxes with pressure stress
               pFx(i,j,k,2)=pFx(i,j,k,1)*0.5_WP*sum(pU(i-1:i,j,k,1))-0.5_WP*sum(pP(i-1:i,j,k,1))
               pFx(i,j,k,3)=pFx(i,j,k,1)*0.5_WP*sum(pV(i-1:i,j,k,1))
               pFx(i,j,k,4)=pFx(i,j,k,1)*0.5_WP*sum(pW(i-1:i,j,k,1))
               ! WENO internal energy flux
               w=weno_weight((abs(pI(i-1,j,k,1)-pI(i-2,j,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(pI(i+1,j,k,1)-pI(i  ,j,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
               pFx(i,j,k,5)=0.5_WP*(pFx(i,j,k,1)-abs(pFx(i,j,k,1)))*sum(wenop*pI(i-2:i  ,j,k,1)) &
               &           +0.5_WP*(pFx(i,j,k,1)+abs(pFx(i,j,k,1)))*sum(wenom*pI(i-1:i+1,j,k,1))
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
               ! Viscous stress at x-face (added to momentum fluxes)
               pFx(i,j,k,2)=pFx(i,j,k,2)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(1,1)+gradU(1,1))+0.5_WP*(sum(pBeta(i-1:i,j,k,1))-2.0_WP/3.0_WP*sum(pVisc(i-1:i,j,k,1)))*div
               pFx(i,j,k,3)=pFx(i,j,k,3)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(2,1)+gradU(1,2))
               pFx(i,j,k,4)=pFx(i,j,k,4)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(3,1)+gradU(1,3))
               ! Heat diffusion flux
               pFx(i,j,k,5)=pFx(i,j,k,5)+0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pT(i,j,k,1)-pT(i-1,j,k,1))
            end do; end do; end do
            
            ! Y-fluxes
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Face velocity
               vel=0.5_WP*sum(pV(i,j-1:j,k,1))
               ! WENO mass flux
               w=weno_weight((abs(pQ(i,j-1,k,1)-pQ(i,j-2,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(pQ(i,j+1,k,1)-pQ(i,j  ,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
               pFy(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j-2:j  ,k,1)) &
               &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j-1:j+1,k,1))
               ! Momentum fluxes with pressure stress
               pFy(i,j,k,2)=pFy(i,j,k,1)*0.5_WP*sum(pU(i,j-1:j,k,1))
               pFy(i,j,k,3)=pFy(i,j,k,1)*0.5_WP*sum(pV(i,j-1:j,k,1))-0.5_WP*sum(pP(i,j-1:j,k,1))
               pFy(i,j,k,4)=pFy(i,j,k,1)*0.5_WP*sum(pW(i,j-1:j,k,1))
               ! WENO internal energy flux
               w=weno_weight((abs(pI(i,j-1,k,1)-pI(i,j-2,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(pI(i,j+1,k,1)-pI(i,j  ,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
               pFy(i,j,k,5)=0.5_WP*(pFy(i,j,k,1)-abs(pFy(i,j,k,1)))*sum(wenop*pI(i,j-2:j  ,k,1)) &
               &           +0.5_WP*(pFy(i,j,k,1)+abs(pFy(i,j,k,1)))*sum(wenom*pI(i,j-1:j+1,k,1))
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
               ! Viscous stress at y-face (added to momentum fluxes)
               pFy(i,j,k,2)=pFy(i,j,k,2)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(1,2)+gradU(2,1))
               pFy(i,j,k,3)=pFy(i,j,k,3)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(2,2)+gradU(2,2))+0.5_WP*(sum(pBeta(i,j-1:j,k,1))-2.0_WP/3.0_WP*sum(pVisc(i,j-1:j,k,1)))*div
               pFy(i,j,k,4)=pFy(i,j,k,4)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(3,2)+gradU(2,3))
               ! Heat diffusion flux
               pFy(i,j,k,5)=pFy(i,j,k,5)+0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pT(i,j,k,1)-pT(i,j-1,k,1))
            end do; end do; end do
            
            ! Z-fluxes
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Face velocity
               vel=0.5_WP*sum(pW(i,j,k-1:k,1))
               ! WENO mass flux
               w=weno_weight((abs(pQ(i,j,k-1,1)-pQ(i,j,k-2,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(pQ(i,j,k+1,1)-pQ(i,j,k  ,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
               pFz(i,j,k,1)=-0.5_WP*(vel+abs(vel))*sum(wenop*pQ(i,j,k-2:k,1)) &
               &            -0.5_WP*(vel-abs(vel))*sum(wenom*pQ(i,j,k-1:k+1,1))
               ! Momentum fluxes with pressure stress
               pFz(i,j,k,2)=pFz(i,j,k,1)*0.5_WP*sum(pU(i,j,k-1:k,1))
               pFz(i,j,k,3)=pFz(i,j,k,1)*0.5_WP*sum(pV(i,j,k-1:k,1))
               pFz(i,j,k,4)=pFz(i,j,k,1)*0.5_WP*sum(pW(i,j,k-1:k,1))-0.5_WP*sum(pP(i,j,k-1:k,1))
               ! WENO internal energy flux
               w=weno_weight((abs(pI(i,j,k-1,1)-pI(i,j,k-2,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(pI(i,j,k+1,1)-pI(i,j,k  ,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
               pFz(i,j,k,5)=0.5_WP*(pFz(i,j,k,1)-abs(pFz(i,j,k,1)))*sum(wenop*pI(i,j,k-2:k  ,1)) &
               &           +0.5_WP*(pFz(i,j,k,1)+abs(pFz(i,j,k,1)))*sum(wenom*pI(i,j,k-1:k+1,1))
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
               ! Viscous stress at z-face (added to momentum fluxes)
               pFz(i,j,k,2)=pFz(i,j,k,2)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(1,3)+gradU(3,1))
               pFz(i,j,k,3)=pFz(i,j,k,3)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(2,3)+gradU(3,2))
               pFz(i,j,k,4)=pFz(i,j,k,4)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(3,3)+gradU(3,3))+0.5_WP*(sum(pBeta(i,j,k-1:k,1))-2.0_WP/3.0_WP*sum(pVisc(i,j,k-1:k,1)))*div
               ! Heat diffusion flux
               pFz(i,j,k,5)=pFz(i,j,k,5)+0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pT(i,j,k,1)-pT(i,j,k-1,1))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

      end do

      ! Phase 2: Average down all fluxes for C/F conservation
      do lvl=this%amr%clvl(),1,-1
         call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
      end do
         
      ! Phase 3: Compute divergence for all levels
      do lvl=0,this%amr%clvl()

         ! Grid spacings
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())

            ! Get pointers to data
            rhs=>dQdt%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi)
            pFy=>Fy(lvl)%dataptr(mfi)
            pFz=>Fz(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)

            ! Loop over interior
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Advection
               rhs(i,j,k,:)=dxi*(pFx(i+1,j,k,:)-pFx(i,j,k,:))+dyi*(pFy(i,j+1,k,:)-pFy(i,j,k,:))+dzi*(pFz(i,j,k+1,:)-pFz(i,j,k,:))
               ! Pressure dilatation
               rhs(i,j,k,5)=rhs(i,j,k,5)-pP(i,j,k,1)*(0.5_WP*dxi*(pU(i+1,j,k,1)-pU(i-1,j,k,1))+0.5_WP*dyi*(pV(i,j+1,k,1)-pV(i,j-1,k,1))+0.5_WP*dzi*(pW(i,j,k+1,1)-pW(i,j,k-1,1)))
               ! Viscous heating: compute cell-centered gradU and stress tensor
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
               ! τ:∇U = τ_ij * gradU(i,j)
               rhs(i,j,k,5)=rhs(i,j,k,5) &
               & +(2.0_WP*pVisc(i,j,k,1)*gradU(1,1)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(1,1) &
               & +(2.0_WP*pVisc(i,j,k,1)*gradU(2,2)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(2,2) &
               & +(2.0_WP*pVisc(i,j,k,1)*gradU(3,3)+(pBeta(i,j,k,1)-2.0_WP/3.0_WP*pVisc(i,j,k,1))*div)*gradU(3,3) &
               & +pVisc(i,j,k,1)*(gradU(2,1)+gradU(1,2))*(gradU(2,1)+gradU(1,2)) &
               & +pVisc(i,j,k,1)*(gradU(3,1)+gradU(1,3))*(gradU(3,1)+gradU(1,3)) &
               & +pVisc(i,j,k,1)*(gradU(3,2)+gradU(2,3))*(gradU(3,2)+gradU(2,3))
            end do; end do; end do

         end do
         call this%amr%mfiter_destroy(mfi)
      end do

      ! Cleanup flux mfabs
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Fx(lvl))
         call this%amr%mfab_destroy(Fy(lvl))
         call this%amr%mfab_destroy(Fz(lvl))
      end do
      
   contains
      !> WENO switch function
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP
         real(WP), parameter :: delta=0.01_WP
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight
   end subroutine get_dQdt

   !> Add artificial bulk viscosity to this%beta
   subroutine add_viscartif(this,dt,Cartif)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cartif
      ! Local variables
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: beta_t,scratch
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBeta_t,pScratch,pU,pV,pW,pC,pBeta,pQ
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_beta,myCartif
      real(WP) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy,vort,grad_div
      integer :: lvl,i,j,k
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2

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

         ! Phase 4: Convert to dynamic viscosity and add to this%beta
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pBeta_t=>beta_t%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBeta(i,j,k,1)=pBeta(i,j,k,1)+pBeta_t(i,j,k,1)*pQ(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy temp multifab
         call amrex_multifab_destroy(beta_t)
         
      end do

      ! Average down and fill
      call this%beta%average_down(); call this%beta%fill(time=0.0_WP)

   end subroutine add_viscartif

   !> Add Vreman SGS eddy viscosity to this%visc
   subroutine add_vreman(this,dt,Cs)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      ! Local variables
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: visc_t
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pU,pV,pW,pVisc,pQ
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,Cmodel,Aij,Bij
      real(WP), dimension(1:3,1:3) :: gradU,betaij
      integer :: lvl,i,j,k,si,sj
      ! Parameters
      real(WP), parameter :: max_cfl=0.5_WP
      integer, parameter :: nfilter=2

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

         ! Phase 3: Convert to dynamic viscosity and add to this%visc
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)*pQ(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)

         ! Destroy temp multifab
         call amrex_multifab_destroy(visc_t)
         
      end do

      ! Average down and fill
      call this%visc%average_down(); call this%visc%fill(time=0.0_WP)

   end subroutine add_vreman

   !> Calculate CFL numbers
   subroutine get_cfl(this,dt,cfl)
      class(amrcomp), intent(inout) :: this
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
         ! Max viscosities
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
                  rho=max(pQ(i,j,k,1),this%rho_floor)
                  viscmax=max(viscmax,pVisc(i,j,k,1)/rho,pBeta(i,j,k,1)/rho,pDiff(i,j,k,1)/rho) ! Heat diffusion cfl is incorrect
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

   !> Calculate monitoring info
   subroutine get_info(this)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer :: lvl,n
      real(WP) :: dV

      ! Compute extrema across all levels
      this%Umax=0.0_WP; this%Vmax=0.0_WP; this%Wmax=0.0_WP
      this%Imin=huge(1.0_WP); this%Imax=-huge(1.0_WP)
      this%Pmin=huge(1.0_WP); this%Pmax=-huge(1.0_WP)
      this%Tmin=huge(1.0_WP); this%Tmax=-huge(1.0_WP)
      this%Cmin=huge(1.0_WP); this%Cmax=-huge(1.0_WP)
      this%Qmin=huge(1.0_WP); this%Qmax=-huge(1.0_WP)
      do lvl=0,this%amr%clvl()
         ! Velocity norm 0
         this%Umax=max(this%Umax,this%U%norm0(lvl=lvl))
         this%Vmax=max(this%Vmax,this%V%norm0(lvl=lvl))
         this%Wmax=max(this%Wmax,this%W%norm0(lvl=lvl))
         ! Extrema of internal energy, pressure, and temperature
         this%Imin=min(this%Imin,this%I%get_min(lvl=lvl)); this%Imax=max(this%Imax,this%I%get_max(lvl=lvl))
         this%Pmin=min(this%Pmin,this%P%get_min(lvl=lvl)); this%Pmax=max(this%Pmax,this%P%get_max(lvl=lvl))
         this%Tmin=min(this%Tmin,this%T%get_min(lvl=lvl)); this%Tmax=max(this%Tmax,this%T%get_max(lvl=lvl))
         ! Extrema of speed of sound
         this%Cmin=min(this%Cmin,this%C%get_min(lvl=lvl)); this%Cmax=max(this%Cmax,this%C%get_max(lvl=lvl))
         ! Extrema of conserved variables
         do n=1,this%Q%ncomp
            this%Qmin(n)=min(this%Qmin(n),this%Q%get_min(lvl=lvl,comp=n))
            this%Qmax(n)=max(this%Qmax(n),this%Q%get_max(lvl=lvl,comp=n))
         end do
      end do

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
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         integer :: i,j,k,ierr
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
                  ! Accumulate kinetic energy
                  this%rhoKint=this%rhoKint+0.5_WP*pQ(i,j,k,1)*(pU(i,j,k,1)**2+pV(i,j,k,1)**2+pW(i,j,k,1)**2)*dV
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         ! Reduce across MPI ranks
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_rhoKint
      
   end subroutine get_info

   !> Print solver info to screen
   subroutine amrcomp_print(this)
      use messager, only: log
      implicit none
      class(amrcomp), intent(in) :: this
      call log("Compressible solver: "//trim(this%name))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrcomp_print

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%Q,'Q')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      call io%read_data(dirname,this%Q,'Q')
      ! Fill ghost cells (IO reads valid data only)
      call this%Q%fill(time=time)
   end subroutine restore_checkpoint
end module amrcomp_class
