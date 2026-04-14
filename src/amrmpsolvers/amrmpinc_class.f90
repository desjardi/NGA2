!> AMR collocated incompressible multiphase solver class
!> Inherits from amrvof_class
module amrmpinc_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrdata_class,    only: amrdata
   use amrmpflow_class,  only: amrmpflow
   use amrmg_class,      only: amrmg
   use amrvof_class,     only: VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter
   implicit none
   private

   ! Expose type
   public :: amrmpinc,VFlo,VFhi,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER

   !> AMR collocated incompressible multiphase solver type
   type, extends(amrmpflow) :: amrmpinc

      ! User-configurable callbacks
      procedure(mpinc_init_iface),    pointer, pass :: user_init   =>null()
      procedure(mpinc_tagging_iface), pointer, pass :: user_tagging=>null()
      procedure(mpinc_bc_iface),      pointer, pass :: user_bc     =>null()
      procedure(mpinc_vofbc_iface),   pointer, pass :: user_vofbc  =>null()

      ! Pressure
      type(amrdata) :: P

      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rhoL=1.0_WP           !< Constant liquid density
      real(WP) :: rhoG=1.0_WP           !< Constant gas density
      type(amrdata) :: visc             !< Variable dynamic viscosity
      real(WP) :: sigma=0.0_WP          !< Surface tension coefficient

      ! Monitoring quantities
      real(WP) :: Pmax=0.0_WP           !< Max pressure
      real(WP) :: rhoUint=0.0_WP        !< Integral of rho*U
      real(WP) :: rhoVint=0.0_WP        !< Integral of rho*V
      real(WP) :: rhoWint=0.0_WP        !< Integral of rho*W
      real(WP) :: rhoKint=0.0_WP        !< Integral of rho*K
      real(WP) :: CFLv_x=0.0_WP         !< Viscous CFL in x
      real(WP) :: CFLv_y=0.0_WP         !< Viscous CFL in y
      real(WP) :: CFLv_z=0.0_WP         !< Viscous CFL in z
      real(WP) :: CFLst=0.0_WP          !< Surface tension CFL
      real(WP) :: CFL=0.0_WP            !< Maximum CFL

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
      procedure :: apply_velbc=>mpinc_apply_velbc
      procedure :: apply_Qbc=>mpinc_apply_Qbc
      procedure :: apply_vofbc=>mpinc_apply_vofbc
      ! Utilities
      procedure :: get_face_velocity          !< Interpolate cell-centered velocity to face
      procedure :: prepare_psolver            !< Prepare pressure solver with new densities
      procedure :: add_pressure               !< Add pressure term consistently to face and cell-centered velocities
      procedure :: add_surface_tension        !< Add surface tension increment consistently to face and cell-centered velocities
      procedure, private :: apply_face_fluxes !< Apply pre-built face fluxes to both face and cell-centered velocities
      ! Physics procedures
      procedure :: get_dQdt                   !< Compute rate of change of Q=UVW
      procedure :: add_vreman                 !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                    !< Compute CFL numbers
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrmpinc_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrmpinc

   !> Abstract interface for user-provided on_init callback
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

   !> Abstract interface for user-provided tagging callback
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
   abstract interface
      subroutine mpinc_bc_iface(solver,lvl,time,face,bx,comp,p)
         import :: amrmpinc,amrex_box,WP
         class(amrmpinc), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< Can be 'U','V','W','Q'
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      end subroutine mpinc_bc_iface
   end interface

   !> Abstract interface for user-provided VOF BC callback
   abstract interface
      subroutine mpinc_vofbc_iface(solver,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
         import :: amrmpinc,amrex_box,WP
         class(amrmpinc), intent(inout) :: solver
         integer, intent(in) :: lvl,face
         real(WP), intent(in) :: time
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      end subroutine
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
      if (associated(this%user_init)) call this%user_init(lvl,time,ba,dm)
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
      if (associated(this%user_tagging)) call this%user_tagging(lvl,time,tags)
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
      use amrgrid_class,    only: amrgrid
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrmg_class,      only: amrmg_varcoef
      implicit none
      class(amrmpinc), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Initialize amrmpflow parent with 3 conserved components and at least 2 ghost cells
      this%nQ=3; this%nover=max(this%nover,2)
      call this%amrmpflow%initialize(amr,name); call this%set_parent()

      ! Initialize pressure with Neumann BCs
      call this%P%initialize(amr,name='P',ncomp=1,ng=this%nover); this%P%parent=>this
      if (.not.amr%xper) then; this%P%lo_bc(1,1)=amrex_bc_foextrap; this%P%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%P%lo_bc(2,1)=amrex_bc_foextrap; this%P%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%P%lo_bc(3,1)=amrex_bc_foextrap; this%P%hi_bc(3,1)=amrex_bc_foextrap; end if

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
      call this%P%finalize()
      call this%visc%finalize()
      call this%psolver%finalize()
      nullify(this%user_init)
      nullify(this%user_tagging)
      nullify(this%user_bc)
      nullify(this%user_vofbc)
      call this%amrmpflow%finalize()
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
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_init(lvl,time,ba,dm)
      ! Reset level layouts
      call this%P%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      ! Set to zero
      call this%P%setval(val=0.0_WP,lvl=lvl)
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
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_coarse(lvl,time,ba,dm)
      ! Pressure uses default on_coarse
      call this%P%on_coarse(lvl,time,ba,dm)
      ! Viscosity just needs to be reset
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
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_remake(lvl,time,ba,dm)
      ! Pressure remake
      call this%P%on_remake(lvl,time,ba,dm)
      ! Viscosity just needs to be reset
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles face velocities, divergence, conserved quantities, and VOF-related quantities
      call this%amrmpflow%on_clear(lvl)
      ! Pressure
      call this%P%clear_level(lvl)
      ! Viscosity
      call this%visc%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles face velocities, conserved quantities, and VOF-related quantities
      call this%amrmpflow%post_regrid(lbase,time)
      ! Average down pressure and fill ghosts
      call this%P%average_down(lbase)
      call this%P%fill(time,lbase)
   end subroutine post_regrid

   ! ============================================================================
   ! BOUNDARY CONDITIONS
   ! ============================================================================

   !> Velocity BC override: forward to user_bc with U/V/W component name
   subroutine mpinc_apply_velbc(this,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp=comp,p=p)
   end subroutine mpinc_apply_velbc

   !> Q BC override: forward to user_bc with comp='Q'
   subroutine mpinc_apply_Qbc(this,lvl,time,face,bx,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp='Q',p=p)
   end subroutine mpinc_apply_Qbc

   !> VOF BC override: forward to user_vofbc
   subroutine mpinc_apply_vofbc(this,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      if (associated(this%user_vofbc)) call this%user_vofbc(lvl=lvl,time=time,face=face,bx=bx,pVF=pVF,pCL=pCL,pCG=pCG,pPLIC=pPLIC)
   end subroutine mpinc_apply_vofbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Update face velocity from Q using density-weighting
   subroutine get_face_velocity(this)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pSubVF
      real(WP) :: rho_lo,rho_hi
      ! Traverse levels
      do lvl=0,this%amr%clvl()
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
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
               pU(i,j,k,1)=(rho_lo*pQ(i-1,j,k,1)+rho_hi*pQ(i,j,k,1))/(rho_lo+rho_hi)
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
               pV(i,j,k,1)=(rho_lo*pQ(i,j-1,k,2)+rho_hi*pQ(i,j,k,2))/(rho_lo+rho_hi)
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
               pW(i,j,k,1)=(rho_lo*pQ(i,j,k-1,3)+rho_hi*pQ(i,j,k,3))/(rho_lo+rho_hi)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_face_velocity

   !> Prepare variable-coefficient pressure solver using face densities
   subroutine prepare_psolver(this)
      use amrex_amr_module, only: amrex_multifab
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

   !> Add (-scale*pressure gradient/rho) to both U/V/W and Q=UVW velocities. Two flavors:
   !>   phi present -> direct path: use explicit stencil that reads phi ghost cells directly (for predictor with fs%P)
   !>   phi absent  -> MLMG path:   use psolver internal fluxes (for projection with dP)
   !> Cell-center correction averages the face gradients back to cell center
   subroutine add_pressure(this,scale,phi)
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface,  only: amrmfab_average_down_face
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
         ! Enforce flux consistency between levels
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
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
      deallocate(Fx,Fy,Fz)
   end subroutine add_pressure

   !> Add CSF surface-tension velocity increment to both U/V/W and Q=UVW velocities
   !> Builds face fluxes 1/rho*sigma*kappa*grad(VF) at maxlvl, avg_down, then calls apply_face_fluxes
   subroutine add_surface_tension(this,scale)
      use amrex_amr_module, only: amrex_multifab
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
   !> scale * Fx/Fy/Fz is saxpy'd onto U/V/W, then density-weighted to Q=UVW
   subroutine apply_face_fluxes(this,scale,Fx,Fy,Fz)
      use amrex_amr_module, only: amrex_multifab
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrex_multifab), intent(in) :: Fx(0:),Fy(0:),Fz(0:)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pQ,pVF,pSubVF
      real(WP) :: rho,rhoLo,rhoHi
      integer :: lvl,i,j,k
      do lvl=0,this%amr%clvl()
         ! Face velocities
         call this%U%mf(lvl)%saxpy(scale,Fx(lvl),1,1,1,0)
         call this%V%mf(lvl)%saxpy(scale,Fy(lvl),1,1,1,0)
         call this%W%mf(lvl)%saxpy(scale,Fz(lvl),1,1,1,0)
         ! Cell-centered velocities: density-weighted average of face fluxes
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pSubVF=>this%subVF%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
               ! X
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,1)+this%rhoG*(1.0_WP-pSubVF(i,j,k,1))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,2)+this%rhoG*(1.0_WP-pSubVF(i,j,k,2))
               pQ(i,j,k,1)=pQ(i,j,k,1)+scale*0.5_WP*(rhoLo*pFx(i,j,k,1)+rhoHi*pFx(i+1,j,k,1))/rho
               ! Y
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,3)+this%rhoG*(1.0_WP-pSubVF(i,j,k,3))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,4)+this%rhoG*(1.0_WP-pSubVF(i,j,k,4))
               pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.5_WP*(rhoLo*pFy(i,j,k,1)+rhoHi*pFy(i,j+1,k,1))/rho
               ! Z
               rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,5)+this%rhoG*(1.0_WP-pSubVF(i,j,k,5))
               rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,6)+this%rhoG*(1.0_WP-pSubVF(i,j,k,6))
               pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.5_WP*(rhoLo*pFz(i,j,k,1)+rhoHi*pFz(i,j,k+1,1))/rho
            end do; end do; end do
            ! Non-periodic boundary cells (only one face flux available)
            if (.not.this%amr%xper.and.bx%lo(1).eq.this%amr%geom(lvl)%domain%lo(1)) then
               i=this%amr%geom(lvl)%domain%lo(1); do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,2)+this%rhoG*(1.0_WP-pSubVF(i,j,k,2))
                  pQ(i,j,k,1)=pQ(i,j,k,1)+scale*0.5_WP*rhoHi/rho*pFx(i+1,j,k,1)
               end do; end do
            end if
            if (.not.this%amr%xper.and.bx%hi(1).eq.this%amr%geom(lvl)%domain%hi(1)) then
               i=this%amr%geom(lvl)%domain%hi(1); do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,1)+this%rhoG*(1.0_WP-pSubVF(i,j,k,1))
                  pQ(i,j,k,1)=pQ(i,j,k,1)+scale*0.5_WP*rhoLo/rho*pFx(i,  j,k,1)
               end do; end do
            end if
            if (.not.this%amr%yper.and.bx%lo(2).eq.this%amr%geom(lvl)%domain%lo(2)) then
               j=this%amr%geom(lvl)%domain%lo(2); do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,4)+this%rhoG*(1.0_WP-pSubVF(i,j,k,4))
                  pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.5_WP*rhoHi/rho*pFy(i,j+1,k,1)
               end do; end do
            end if
            if (.not.this%amr%yper.and.bx%hi(2).eq.this%amr%geom(lvl)%domain%hi(2)) then
               j=this%amr%geom(lvl)%domain%hi(2); do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,3)+this%rhoG*(1.0_WP-pSubVF(i,j,k,3))
                  pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.5_WP*rhoLo/rho*pFy(i,j,  k,1)
               end do; end do
            end if
            if (.not.this%amr%zper.and.bx%lo(3).eq.this%amr%geom(lvl)%domain%lo(3)) then
               k=this%amr%geom(lvl)%domain%lo(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoHi=rho; if (lvl.eq.this%amr%maxlvl) rhoHi=this%rhoL*pSubVF(i,j,k,6)+this%rhoG*(1.0_WP-pSubVF(i,j,k,6))
                  pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.5_WP*rhoHi/rho*pFz(i,j,k+1,1)
               end do; end do
            end if
            if (.not.this%amr%zper.and.bx%hi(3).eq.this%amr%geom(lvl)%domain%hi(3)) then
               k=this%amr%geom(lvl)%domain%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho=this%rhoL*pVF(i,j,k,1)+this%rhoG*(1.0_WP-pVF(i,j,k,1))
                  rhoLo=rho; if (lvl.eq.this%amr%maxlvl) rhoLo=this%rhoL*pSubVF(i,j,k,5)+this%rhoG*(1.0_WP-pSubVF(i,j,k,5))
                  pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.5_WP*rhoLo/rho*pFz(i,j,k,  1)
               end do; end do
            end if
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine apply_face_fluxes

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Compute dQ/dt for all levels without pressure term (user can add it via add_pressure)
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dQdt(this,dQdt,dt,time)
      use amrex_amr_module, only: amrex_multifab
      implicit none
      class(amrmpinc), intent(inout) :: this
      class(amrdata), intent(inout) :: dQdt                        ! Output: momentum RHS (cell-centered)
      real(WP), intent(in) :: dt,time
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz ! Flux mfabs
      type(amrex_multifab) :: Vx,Vy,Vz
      type(amrex_multifab) :: band
      ! Shared variables for internal functions
      real(WP) :: dx,dy,dz,dxi,dyi,dzi                               ! Needed for SL transport
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW  ! Velocity used for project
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPLICold  ! PLICold used in tet2flux_plic
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQold     ! Qold used in tet2flux_plic
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
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pBand,pVx,pVy,pVz,pFx,pFy,pFz,pQ
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
            pQold   =>this%Qold%mf(lvl)%dataptr(mfi)
            pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
            pBand   =>band%dataptr(mfi)
            pQ      =>this%Q%mf(lvl)%dataptr(mfi)
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
               if (.not.crossed_plic) pFx(i,j,k,:)=-(this%rhoL*pVx(i,j,k,1)+this%rhoG*pVx(i,j,k,2))/(dt*dy*dz)*0.5_WP*(pQ(i-1,j,k,:)+pQ(i,j,k,:))
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
               if (.not.crossed_plic) pFy(i,j,k,:)=-(this%rhoL*pVy(i,j,k,1)+this%rhoG*pVy(i,j,k,2))/(dt*dz*dx)*0.5_WP*(pQ(i,j-1,k,:)+pQ(i,j,k,:))
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
               if (.not.crossed_plic) pFz(i,j,k,:)=-(this%rhoL*pVz(i,j,k,1)+this%rhoG*pVz(i,j,k,2))/(dt*dx*dy)*0.5_WP*(pQ(i,j,k-1,:)+pQ(i,j,k,:))
            end do; end do; end do
            ! Deallocate proj for this tile
            deallocate(proj)
         end do
         call this%amr%mfiter_destroy(mfi)
      end block semilagrangian_fluxes

      ! Phase 1b: Finite volume fluxes for all levels (Euler fluxes skip band cells at finest level)
      finitevolume_fluxes: block
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pFx,pFy,pFz,pBand,pVF,pVisc
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
               pQ   =>this%Q%mf(lvl)%dataptr(mfi)
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
                     pFx(i,j,k,:)=mass_flux*0.5_WP*(pQ(i-1,j,k,:)+pQ(i,j,k,:))
                  end if
                  ! Velocity gradients at x-face
                  gradU(1,1)=dxi*(pQ(i,j,k,1)-pQ(i-1,j,k,1))
                  gradU(2,1)=0.25_WP*dyi*sum(pQ(i-1:i,j+1,k,1)-pQ(i-1:i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*sum(pQ(i-1:i,j,k+1,1)-pQ(i-1:i,j,k-1,1))
                  gradU(1,2)=dxi*(pQ(i,j,k,2)-pQ(i-1,j,k,2))
                  gradU(2,2)=0.25_WP*dyi*sum(pQ(i-1:i,j+1,k,2)-pQ(i-1:i,j-1,k,2))
                  gradU(3,2)=0.25_WP*dzi*sum(pQ(i-1:i,j,k+1,2)-pQ(i-1:i,j,k-1,2))
                  gradU(1,3)=dxi*(pQ(i,j,k,3)-pQ(i-1,j,k,3))
                  gradU(2,3)=0.25_WP*dyi*sum(pQ(i-1:i,j+1,k,3)-pQ(i-1:i,j-1,k,3))
                  gradU(3,3)=0.25_WP*dzi*sum(pQ(i-1:i,j,k+1,3)-pQ(i-1:i,j,k-1,3))
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
                     pFy(i,j,k,:)=mass_flux*0.5_WP*(pQ(i,j-1,k,:)+pQ(i,j,k,:))
                  end if
                  ! Velocity gradients at y-face
                  gradU(1,1)=0.25_WP*dxi*sum(pQ(i+1,j-1:j,k,1)-pQ(i-1,j-1:j,k,1))
                  gradU(2,1)=dyi*(pQ(i,j,k,1)-pQ(i,j-1,k,1))
                  gradU(3,1)=0.25_WP*dzi*sum(pQ(i,j-1:j,k+1,1)-pQ(i,j-1:j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*sum(pQ(i+1,j-1:j,k,2)-pQ(i-1,j-1:j,k,2))
                  gradU(2,2)=dyi*(pQ(i,j,k,2)-pQ(i,j-1,k,2))
                  gradU(3,2)=0.25_WP*dzi*sum(pQ(i,j-1:j,k+1,2)-pQ(i,j-1:j,k-1,2))
                  gradU(1,3)=0.25_WP*dxi*sum(pQ(i+1,j-1:j,k,3)-pQ(i-1,j-1:j,k,3))
                  gradU(2,3)=dyi*(pQ(i,j,k,3)-pQ(i,j-1,k,3))
                  gradU(3,3)=0.25_WP*dzi*sum(pQ(i,j-1:j,k+1,3)-pQ(i,j-1:j,k-1,3))
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
                     pFz(i,j,k,:)=mass_flux*0.5_WP*(pQ(i,j,k-1,:)+pQ(i,j,k,:))
                  end if
                  ! Velocity gradients at z-face
                  gradU(1,1)=0.25_WP*dxi*sum(pQ(i+1,j,k-1:k,1)-pQ(i-1,j,k-1:k,1))
                  gradU(2,1)=0.25_WP*dyi*sum(pQ(i,j+1,k-1:k,1)-pQ(i,j-1,k-1:k,1))
                  gradU(3,1)=dzi*(pQ(i,j,k,1)-pQ(i,j,k-1,1))
                  gradU(1,2)=0.25_WP*dxi*sum(pQ(i+1,j,k-1:k,2)-pQ(i-1,j,k-1:k,2))
                  gradU(2,2)=0.25_WP*dyi*sum(pQ(i,j+1,k-1:k,2)-pQ(i,j-1,k-1:k,2))
                  gradU(3,2)=dzi*(pQ(i,j,k,2)-pQ(i,j,k-1,2))
                  gradU(1,3)=0.25_WP*dxi*sum(pQ(i+1,j,k-1:k,3)-pQ(i-1,j,k-1:k,3))
                  gradU(2,3)=0.25_WP*dyi*sum(pQ(i,j+1,k-1:k,3)-pQ(i,j-1,k-1:k,3))
                  gradU(3,3)=dzi*(pQ(i,j,k,3)-pQ(i,j,k-1,3))
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
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pdQdt,pFx,pFy,pFz,pVx,pVy,pVz,pBand
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
               pdQdt=>dQdt%mf(lvl)%dataptr(mfi)
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
                  pdQdt(i,j,k,:)=dxi*(pFx(i+1,j,k,:)-pFx(i,j,k,:))+dyi*(pFy(i,j+1,k,:)-pFy(i,j,k,:))+dzi*(pFz(i,j,k+1,:)-pFz(i,j,k,:))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block divergence_and_sources

      ! Transform momentum rhs to dQdt
      transform_to_velocity: block
         use amrex_amr_module, only: amrex_multifab
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         integer :: lvl,i,j,k
         real(WP) :: rho,rho_old
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pdQdt,pVF,pQold,pVFold
         do lvl=0,this%amr%clvl()
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               pdQdt =>dQdt%mf(lvl)%dataptr(mfi)
               pQold =>this%Qold%mf(lvl)%dataptr(mfi)
               pVFold=>this%VFold%mf(lvl)%dataptr(mfi)
               pVF   =>this%VF%mf(lvl)%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  rho    =this%rhoL*pVF   (i,j,k,1)+this%rhoG*(1.0_WP-pVF   (i,j,k,1))
                  rho_old=this%rhoL*pVFold(i,j,k,1)+this%rhoG*(1.0_WP-pVFold(i,j,k,1))
                  pdQdt(i,j,k,:)=pdQdt(i,j,k,:)/rho+(rho_old-rho)/(rho*dt)*pQold(i,j,k,:)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block transform_to_velocity

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
            myQflux=vol_tot*this%rhoL*pQold(i0,j0,k0,1:3)
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas
            myVflux( 2 )=vol_tot
            myVflux(6:8)=vol_tot*bary_tot
            ! Q flux: pure gas momentum
            myQflux=vol_tot*this%rhoG*pQold(i0,j0,k0,1:3)
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

         ! Compute Q flux from Qold
         myQflux=(this%rhoL*myVflux(1)+this%rhoG*myVflux(2))*pQold(i0,j0,k0,1:3)

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

   !> Add Vreman SGS eddy viscosity to this%visc: assumes velocity ghosts are filled
   !> User must reset visc to molecular value before calling this routine
   subroutine add_vreman(this,dt,Cs)
      use amrsgs, only: get_vreman
      implicit none
      class(amrmpinc), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrdata) :: visc_t
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pVisc,pVF
      ! Create temp amrdata
      call visc_t%initialize(amr=this%amr,name='visc_t',ncomp=1,ng=this%nover); call visc_t%reset()
      ! Compute kinematic eddy viscosity into temp
      call get_vreman(dt=dt,visc=visc_t,U=this%Q,V=this%Q,W=this%Q,Ucomp=1,Vcomp=2,Wcomp=3,Cs=Cs)
      ! Add rho*visc_t to dynamic viscosity
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)/(pVF(i,j,k,1)/this%rhoL+(1.0_WP-pVF(i,j,k,1))/this%rhoG)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temp amrdata
      call visc_t%finalize()
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
      real(WP) :: viscmax
      ! Get convective CFL from parent
      call this%amrmpflow%get_cflc(dt=dt)
      ! Reset CFLs
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP; this%CFLst=0.0_WP
      ! Compute viscous CFL at each level (explicit stability: dt < dx^2 / (4*nu))
      do lvl=0,this%amr%clvl()
         ! Max viscosity
         get_viscmax: block
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

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Get solver information: min/max velocity, min/max pressure, divergence, momentum, TKE
   subroutine get_info(this)
      implicit none
      class(amrmpinc), intent(inout) :: this
      integer :: lvl

      ! Use parent's method first
      call this%amrmpflow%get_info()

      ! Initialize min/max values
      this%Pmax=-huge(1.0_WP)

      ! Loop over all levels for min/max
      do lvl=0,this%amr%clvl()
         this%Umax=max(this%Umax,this%Q%norm0(lvl=lvl,comp=1))
         this%Vmax=max(this%Vmax,this%Q%norm0(lvl=lvl,comp=2))
         this%Wmax=max(this%Wmax,this%Q%norm0(lvl=lvl,comp=3))
         this%Pmax=max(this%Pmax,this%P%norm0(lvl=lvl))
      end do

      ! Integrate momentum and kinetic energy
      get_momentum_and_kinetic_energy: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         integer :: i,j,k,ierr
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF
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
               pQ =>this%Q%mf(lvl)%dataptr(mfi)
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
                  this%rhoUint=this%rhoUint+rho*pQ(i,j,k,1)*this%amr%cell_vol(lvl)
                  this%rhoVint=this%rhoVint+rho*pQ(i,j,k,2)*this%amr%cell_vol(lvl)
                  this%rhoWint=this%rhoWint+rho*pQ(i,j,k,3)*this%amr%cell_vol(lvl)
                  ! Accumulate kinetic energy
                  this%rhoKint=this%rhoKint+0.5_WP*rho*(pQ(i,j,k,1)**2+pQ(i,j,k,2)**2+pQ(i,j,k,3)**2)*this%amr%cell_vol(lvl)
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
      ! Register face velocities, conserved variables, and VOF data via parent
      call this%amrmpflow%register_checkpoint(io)
      ! Register remaining data
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
      ! Restore face velocities, conserved variables, and VOF data via parent
      call this%amrmpflow%restore_checkpoint(io,dirname,time)
      ! Restore remaining data
      call io%read_data(dirname,this%P,'P'); call this%P%fill(time=time)
   end subroutine restore_checkpoint

end module amrmpinc_class
