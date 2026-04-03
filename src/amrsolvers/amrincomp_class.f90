!> AMR Staggered incompressible solver class
!> Provides projection and basic operations for constant-density flow
module amrincomp_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrdata_class,    only: amrdata
   use amrflow_class,    only: amrflow
   use amrmg_class,      only: amrmg
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap
   implicit none
   private

   ! Expose type
   public :: amrincomp

   !> AMR staggered incompressible solver type
   type, extends(amrflow) :: amrincomp

      ! User-configurable callbacks
      procedure(incomp_init_iface),    pointer, pass :: user_init   =>null()  !< User-defined initialization
      procedure(incomp_tagging_iface), pointer, pass :: user_tagging=>null()  !< User-defined tagging
      procedure(incomp_bc_iface),      pointer, pass :: user_bc     =>null()  !< User-defined boundary conditions

      ! Pressure
      type(amrdata) :: P
      
      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rho=1.0_WP            !< Constant density
      type(amrdata) :: visc             !< variable dynamic viscosity

      ! Monitoring quantities
      real(WP) :: Pmax=0.0_WP           !< Max pressure
      real(WP) :: rhoUint=0.0_WP        !< Integral of rho*U
      real(WP) :: rhoVint=0.0_WP        !< Integral of rho*V
      real(WP) :: rhoWint=0.0_WP        !< Integral of rho*W
      real(WP) :: rhoKint=0.0_WP        !< Integral of rho*K
      real(WP) :: CFLv_x=0.0_WP         !< Viscous CFL in x
      real(WP) :: CFLv_y=0.0_WP         !< Viscous CFL in y
      real(WP) :: CFLv_z=0.0_WP         !< Viscous CFL in z
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
      procedure :: apply_velbc=>incomp_apply_velbc
      ! Utilities
      procedure :: correct_velocity          !< Correct face velocity with pressure gradient
      ! Physics procedures
      procedure :: get_dmomdt                !< Compute momentum advection RHS
      procedure :: add_vreman                !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                   !< Compute CFL numbers
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrincomp_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrincomp

   !> Abstract interface for user-provided on_init callback
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

   !> Abstract interface for user-provided tagging callback
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
         character(len=1), intent(in) :: comp              !< Can be 'U','V','W'
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
      use amrgrid_class,    only: amrgrid
      implicit none
      class(amrincomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Initialize amrflow parent without conserved components
      call this%amrflow%initialize(amr=amr,name=name); call this%set_parent()

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
      call this%P%finalize()
      call this%visc%finalize()
      call this%psolver%finalize()
      nullify(this%user_init)
      nullify(this%user_tagging)
      nullify(this%user_bc)
      call this%amrflow%finalize()
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
      ! Parent handles face velocities and divergence
      call this%amrflow%on_init(lvl,time,ba,dm)
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
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities and divergence
      call this%amrflow%on_coarse(lvl,time,ba,dm)
      ! Pressure uses default on_coarse
      call this%P%on_coarse(lvl,time,ba,dm)
      ! Viscosity just needs to be reset
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities and divergence
      call this%amrflow%on_remake(lvl,time,ba,dm)
      ! Pressure remake
      call this%P%on_remake(lvl,time,ba,dm)
      ! Viscosity just needs to be reset
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles face velocities and divergence
      call this%amrflow%on_clear(lvl)
      ! Pressure
      call this%P%clear_level(lvl)
      ! Viscosity
      call this%visc%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles face velocities
      call this%amrflow%post_regrid(lbase,time)
      ! Average down pressure and fill ghosts
      call this%P%average_down(lbase)
      call this%P%fill(time,lbase)
      ! Rebuild pressure solver operators for new grid
      call this%psolver%setup()
   end subroutine post_regrid

   ! ============================================================================
   ! BOUNDARY CONDITIONS
   ! ============================================================================

   !> Velocity BC override: forward to user_bc with U/V/W component name
   subroutine incomp_apply_velbc(this,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrincomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp=comp,p=p)
   end subroutine incomp_apply_velbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

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
   subroutine get_dmomdt(this,drhoUdt,drhoVdt,drhoWdt)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_destroy,amrex_mfiter,amrex_box
      use amrex_interface,  only: amrmfab_average_down_cell,amrmfab_average_down_edge
      implicit none
      class(amrincomp), intent(inout) :: this
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
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
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
         ! Get mesh size
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl,mfi)
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

   !> Add Vreman SGS eddy viscosity to this%visc: assumes velocity ghosts are filled
   !> User must reset visc to molecular value before calling this routine
   subroutine add_vreman(this,dt,Cs)
      use amrsgs, only: get_vreman
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrdata) :: visc_t
      ! Create temp amrdata
      call visc_t%initialize(amr=this%amr,ncomp=1,ng=this%nover,name='visc_t'); call visc_t%reset()
      ! Compute kinematic eddy viscosity into scratch
      call get_vreman(dt=dt,visc=visc_t,U=this%U,V=this%V,W=this%W,Cs=Cs)
      ! Add rho*visc_t to dynamic viscosity
      call this%visc%saxpy(a=this%rho,src=visc_t)
      ! Destroy temp amrdata
      call visc_t%finalize()
   end subroutine add_vreman

   !> Compute CFL numbers (convective and viscous)
   subroutine get_cfl(this,dt,cfl,cflc)
      implicit none
      class(amrincomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP), intent(out), optional :: cflc
      integer :: lvl
      ! Get convective CFL from parent
      call this%amrflow%get_cflc(dt=dt)
      ! Reset CFLs
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute viscous CFL at each level (explicit stability: dt < dx^2 / (4*nu))
      do lvl=0,this%amr%clvl()
         if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*this%visc%norm0(lvl=lvl)*dt/(this%rho*this%amr%dx(lvl)**2))
         if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*this%visc%norm0(lvl=lvl)*dt/(this%rho*this%amr%dy(lvl)**2))
         if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*this%visc%norm0(lvl=lvl)*dt/(this%rho*this%amr%dz(lvl)**2))
      end do
      ! Compute max overall CFL
      this%CFL=max(this%CFLc,this%CFLv_x,this%CFLv_y,this%CFLv_z)
      ! Return max overall CFL
      cfl=this%CFL
      ! Optionally return max convective CFL
      if (present(cflc)) cflc=this%CFLc
   end subroutine get_cfl

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Get solver information: min/max velocity, min/max pressure, divergence, momentum, TKE
   subroutine get_info(this)
      implicit none
      class(amrincomp), intent(inout) :: this
      integer :: lvl

      ! Use parent's method first
      call this%amrflow%get_info()

      ! Initialize min/max values
      this%Pmax=-huge(1.0_WP)

      ! Loop over all levels for min/max
      do lvl=0,this%amr%clvl()
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
      ! Face velocities are registered with parent
      call this%amrflow%register_checkpoint(io)
      ! Register remaining data
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
      ! Restore face velocities via parent
      call this%amrflow%restore_checkpoint(io,dirname,time)
      ! Restore remaining data
      call io%read_data(dirname,this%P,'P'); call this%P%fill(time=time)
   end subroutine restore_checkpoint

end module amrincomp_class
