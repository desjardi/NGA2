!> AMR Collocated incompressible solver class
!> Provides projection and basic operations for constant-density flow
module amrcinc_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_loc,c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrdata_class,    only: amrdata
   use amrflow_class,    only: amrflow
   use amrmg_class,      only: amrmg
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter
   implicit none
   private

   ! Expose type
   public :: amrcinc

   !> AMR collocated incompressible solver type
   type, extends(amrflow) :: amrcinc

      ! User-configurable callbacks
      procedure(cinc_init_iface),    pointer, pass :: user_init   =>null()  !< User-defined initialization
      procedure(cinc_tagging_iface), pointer, pass :: user_tagging=>null()  !< User-defined tagging
      procedure(cinc_bc_iface),      pointer, pass :: user_bc     =>null()  !< User-defined boundary conditions

      ! Pressure
      type(amrdata) :: P

      ! Pressure solver
      type(amrmg) :: psolver

      ! Physical properties
      real(WP) :: rho=1.0_WP            !< Constant density
      type(amrdata) :: visc             !< Variable dynamic viscosity

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
      procedure :: apply_velbc=>cinc_apply_velbc
      procedure :: apply_Qbc=>cinc_apply_Qbc
      ! Utilities
      procedure :: get_face_velocity         !< Update face velocities from cell-centered data
      procedure :: add_pressure              !< Add pressure term consistently to face and cell-centered velocities
      ! Physics procedures
      procedure :: get_dQdt                  !< Compute rate of change of conserved variables
      procedure :: add_vreman                !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                   !< Compute CFL numbers
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrcinc_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrcinc

   !> Abstract interface for user-provided on_init callback
   abstract interface
      subroutine cinc_init_iface(solver,lvl,time,ba,dm)
         import :: amrcinc,WP,amrex_boxarray,amrex_distromap
         class(amrcinc), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine cinc_init_iface
   end interface

   !> Abstract interface for user-provided tagging callback
   abstract interface
      subroutine cinc_tagging_iface(solver,lvl,time,tags)
         import :: amrcinc,c_ptr,WP
         class(amrcinc), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine cinc_tagging_iface
   end interface

   !> Abstract interface for user-provided BC callback
   abstract interface
      subroutine cinc_bc_iface(solver,lvl,time,face,bx,comp,p)
         import :: amrcinc,WP,amrex_box
         class(amrcinc), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< Can be 'U','V','W','Q'
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      end subroutine cinc_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrcinc type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrcinc_on_init(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrcinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(lvl,time,ba,dm)
   end subroutine amrcinc_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrcinc_on_coarse(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrcinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrcinc_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrcinc_on_remake(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrcinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrcinc_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrcinc_on_clear(ctx,lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrcinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrcinc_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrcinc_tagging(ctx,lvl,time,tags)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrcinc), pointer :: this
      call c_f_pointer(ctx,this)
      if (associated(this%user_tagging)) call this%user_tagging(lvl,time,tags)
   end subroutine amrcinc_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrcinc_postregrid(ctx,lbase,time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrcinc), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrcinc_postregrid

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the incompressible solver
   subroutine initialize(this,amr,name)
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrmg_class,      only: amrmg_cstcoef
      use amrgrid_class,    only: amrgrid
      implicit none
      class(amrcinc), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Initialize amrflow parent with 3 conserved components and at least 1 ghost cell
      this%nQ=3; this%nover=max(this%nover,1)
      call this%amrflow%initialize(amr=amr,name=name); call this%set_parent()

      ! Initialize pressure with Neumann BCs
      call this%P%initialize(amr=amr,name='P',ncomp=1,ng=this%nover); this%P%parent=>this
      if (.not.amr%xper) then; this%P%lo_bc(1,1)=amrex_bc_foextrap; this%P%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%P%lo_bc(2,1)=amrex_bc_foextrap; this%P%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%P%lo_bc(3,1)=amrex_bc_foextrap; this%P%hi_bc(3,1)=amrex_bc_foextrap; end if

      ! Initialize viscosity with Neumann BCs
      call this%visc%initialize(amr=amr,name='visc',ncomp=1,ng=this%nover); this%visc%parent=>this
      if (.not.amr%xper) then; this%visc%lo_bc(1,1)=amrex_bc_foextrap; this%visc%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%visc%lo_bc(2,1)=amrex_bc_foextrap; this%visc%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%visc%lo_bc(3,1)=amrex_bc_foextrap; this%visc%hi_bc(3,1)=amrex_bc_foextrap; end if

      ! Initialize pressure solver
      call this%psolver%initialize(amr=amr,type=amrmg_cstcoef)

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      select type (this)
       type is (amrcinc)
         call this%amr%add_on_init   (amrcinc_on_init,   c_loc(this))
         call this%amr%add_on_coarse (amrcinc_on_coarse, c_loc(this))
         call this%amr%add_on_remake (amrcinc_on_remake, c_loc(this))
         call this%amr%add_on_clear  (amrcinc_on_clear,  c_loc(this))
         call this%amr%add_tagging   (amrcinc_tagging,   c_loc(this))
         call this%amr%add_postregrid(amrcinc_postregrid,c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the incompressible solver
   subroutine finalize(this)
      implicit none
      class(amrcinc), intent(inout) :: this
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
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, and conserved quantities
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
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_coarse(lvl,time,ba,dm)
      ! Pressure uses default on_coarse
      call this%P%on_coarse(lvl,time,ba,dm)
      ! Viscosity just needs to be reset
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using divergence-free interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_remake(lvl,time,ba,dm)
      ! Pressure remake
      call this%P%on_remake(lvl,time,ba,dm)
      ! Viscosity just needs to be reset
      call this%visc%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_clear(lvl)
      ! Pressure
      call this%P%clear_level(lvl)
      ! Viscosity
      call this%visc%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles face velocities and conserved quantities
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
   subroutine cinc_apply_velbc(this,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp=comp,p=p)
   end subroutine cinc_apply_velbc

   !> Q BC override: forward to user_bc with comp='Q'
   subroutine cinc_apply_Qbc(this,lvl,time,face,bx,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrcinc), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp='Q',p=p)
   end subroutine cinc_apply_Qbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Update face velocity from Q
   subroutine get_face_velocity(this)
      implicit none
      class(amrcinc), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW
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
            ! Get X-face velocity
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pU(i,j,k,1)=0.5_WP*sum(pQ(i-1:i,j,k,1))
            end do; end do; end do
            ! Get Y-face velocity
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pV(i,j,k,1)=0.5_WP*sum(pQ(i,j-1:j,k,2))
            end do; end do; end do
            ! Get Z-face velocity
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pW(i,j,k,1)=0.5_WP*sum(pQ(i,j,k-1:k,3))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_face_velocity

   !> Add (-scale*pressure gradient) to both U/V/W and Q=UVW velocities. Two flavors:
   !>   phi present -> direct path: use explicit stencil that reads phi ghost cells directly (for predictor with fs%P)
   !>   phi absent  -> MLMG path:   use psolver internal fluxes (for projection with dP)
   !> Cell-center correction averages the face gradients back to cell center
   subroutine add_pressure(this,scale,phi)
      use amrex_amr_module, only: amrex_multifab
      class(amrcinc), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrdata), intent(in), optional :: phi
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP,pFx,pFy,pFz,pQ
      real(WP) :: dxi,dyi,dzi
      integer :: lvl,i,j,k
      ! Build temp face mfabs to store pressure fluxes
      allocate(Fx(0:this%amr%clvl()),Fy(0:this%amr%clvl()),Fz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Compute -pressure gradient at faces
      if (present(phi)) then
         ! Use provided phi and its ghosts cells
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
         ! Use psolver's solution and its internal ghosts
         call this%psolver%get_fluxes(Fx,Fy,Fz)
      end if
      ! Apply to face velocities and cell-centered in one pass
      do lvl=0,this%amr%clvl()
         ! Face: use flux directly
         call this%U%mf(lvl)%saxpy(scale,Fx(lvl),1,1,1,0)
         call this%V%mf(lvl)%saxpy(scale,Fy(lvl),1,1,1,0)
         call this%W%mf(lvl)%saxpy(scale,Fz(lvl),1,1,1,0)
         ! Cell-center: average flux to cell center
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pQ(i,j,k,1)=pQ(i,j,k,1)+scale*0.5_WP*sum(pFx(i:i+1,j,k,1))
               pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.5_WP*sum(pFy(i,j:j+1,k,1))
               pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.5_WP*sum(pFz(i,j,k:k+1,1))
            end do; end do; end do
            ! Fix non-periodic boundary conditions
            if (.not.this%amr%xper.and.bx%lo(1).eq.this%amr%geom(lvl)%domain%lo(1)) then
               i=this%amr%geom(lvl)%domain%lo(1); do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
                  pQ(i,j,k,1)=pQ(i,j,k,1)+scale*0.5_WP*pFx(i+1,j,k,1)
               end do; end do
            end if
            if (.not.this%amr%xper.and.bx%hi(1).eq.this%amr%geom(lvl)%domain%hi(1)) then
               i=this%amr%geom(lvl)%domain%hi(1); do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)
                  pQ(i,j,k,1)=pQ(i,j,k,1)+scale*0.5_WP*pFx(i,  j,k,1)
               end do; end do
            end if
            if (.not.this%amr%yper.and.bx%lo(2).eq.this%amr%geom(lvl)%domain%lo(2)) then
               j=this%amr%geom(lvl)%domain%lo(2); do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
                  pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.5_WP*pFy(i,j+1,k,1)
               end do; end do
            end if
            if (.not.this%amr%yper.and.bx%hi(2).eq.this%amr%geom(lvl)%domain%hi(2)) then
               j=this%amr%geom(lvl)%domain%hi(2); do k=bx%lo(3),bx%hi(3); do i=bx%lo(1),bx%hi(1)
                  pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.5_WP*pFy(i,j,  k,1)
               end do; end do
            end if
            if (.not.this%amr%zper.and.bx%lo(3).eq.this%amr%geom(lvl)%domain%lo(3)) then
               k=this%amr%geom(lvl)%domain%lo(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.5_WP*pFz(i,j,k+1,1)
               end do; end do
            end if
            if (.not.this%amr%zper.and.bx%hi(3).eq.this%amr%geom(lvl)%domain%hi(3)) then
               k=this%amr%geom(lvl)%domain%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.5_WP*pFz(i,j,k,  1)
               end do; end do
            end if
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temps
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_destroy(Fx(lvl))
         call this%amr%mfab_destroy(Fy(lvl))
         call this%amr%mfab_destroy(Fz(lvl))
      end do
   end subroutine add_pressure

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Compute dQ/dt for all levels without pressure term (user can add it via add_pressure)
   !> Uses flux averaging at C/F interfaces for conservation
   subroutine get_dQdt(this,dQdt)
      use amrex_amr_module, only: amrex_multifab
      implicit none
      class(amrcinc), intent(inout) :: this
      class(amrdata), intent(inout) :: dQdt                           ! Output: rate of change of conserved variables
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz  ! Flux mfabs

      ! Initialize all fluxes
      define_fluxes: block
         integer :: lvl
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_build(lvl,Fx(lvl),ncomp=this%nQ,nover=1,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl,Fy(lvl),ncomp=this%nQ,nover=1,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl,Fz(lvl),ncomp=this%nQ,nover=1,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
         end do
      end block define_fluxes

      ! Compute fluxes on all levels
      compute_fluxes: block
         integer :: lvl,i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         real(WP) :: dxi,dyi,dzi
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP) :: div,visc_f
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pQ,pVisc
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            ! Loop over all tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get pointers to data
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi)
               pFy=>Fy(lvl)%dataptr(mfi)
               pFz=>Fz(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               ! Compute X-fluxes
               fbx=mfi%nodaltilebox(1)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Convective fluxes
                  pFx(i,j,k,1)=-this%rho*pU(i,j,k,1)*0.5_WP*sum(pQ(i-1:i,j,k,1))
                  pFx(i,j,k,2)=-this%rho*pU(i,j,k,1)*0.5_WP*sum(pQ(i-1:i,j,k,2))
                  pFx(i,j,k,3)=-this%rho*pU(i,j,k,1)*0.5_WP*sum(pQ(i-1:i,j,k,3))
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
               ! Compute Y-fluxes
               fbx=mfi%nodaltilebox(2)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Convective fluxes
                  pFy(i,j,k,1)=-this%rho*pV(i,j,k,1)*0.5_WP*sum(pQ(i,j-1:j,k,1))
                  pFy(i,j,k,2)=-this%rho*pV(i,j,k,1)*0.5_WP*sum(pQ(i,j-1:j,k,2))
                  pFy(i,j,k,3)=-this%rho*pV(i,j,k,1)*0.5_WP*sum(pQ(i,j-1:j,k,3))
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
               ! Compute Z-fluxes
               fbx=mfi%nodaltilebox(3)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Convective fluxes
                  pFz(i,j,k,1)=-this%rho*pW(i,j,k,1)*0.5_WP*sum(pQ(i,j,k-1:k,1))
                  pFz(i,j,k,2)=-this%rho*pW(i,j,k,1)*0.5_WP*sum(pQ(i,j,k-1:k,2))
                  pFz(i,j,k,3)=-this%rho*pW(i,j,k,1)*0.5_WP*sum(pQ(i,j,k-1:k,3))
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
      end block compute_fluxes

      ! Average down all fluxes for C/F conservation
      c_f_consistency: block
         use amrex_interface, only: amrmfab_average_down_face
         integer :: lvl
         do lvl=this%amr%clvl(),1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         end do
      end block c_f_consistency

      ! Compute divergence to get momentum RHS
      divergence_and_sources: block
         integer :: lvl,i,j,k,n
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: dxi,dyi,dzi,irho
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pdQdt
         ! Compute inverse density
         irho=1.0_WP/this%rho
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            ! Loop over all tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get pointers to data
               pdQdt=>dQdt%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi)
               pFy=>Fy(lvl)%dataptr(mfi)
               pFz=>Fz(lvl)%dataptr(mfi)
               ! Compute divergence and divide by rho
               bx=mfi%tilebox()
               do n=1,3; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdQdt(i,j,k,n)=dxi*(pFx(i+1,j,k,n)-pFx(i,j,k,n))+dyi*(pFy(i,j+1,k,n)-pFy(i,j,k,n))+dzi*(pFz(i,j,k+1,n)-pFz(i,j,k,n))
                  pdQdt(i,j,k,n)=irho*pdQdt(i,j,k,n)
               end do; end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block divergence_and_sources

      ! Cleanup flux MultiFabs
      cleanup: block
         integer :: lvl
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_destroy(Fx(lvl))
            call this%amr%mfab_destroy(Fy(lvl))
            call this%amr%mfab_destroy(Fz(lvl))
         end do
      end block cleanup

   end subroutine get_dQdt

   !> Add Vreman SGS eddy viscosity to this%visc: assumes velocity ghosts are filled
   !> User must reset visc to molecular value before calling this routine
   subroutine add_vreman(this,dt,Cs)
      use amrsgs, only: get_vreman
      implicit none
      class(amrcinc), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrdata) :: visc_t
      ! Create temp amrdata
      call visc_t%initialize(amr=this%amr,name='visc_t',ncomp=1,ng=this%nover); call visc_t%reset()
      ! Compute kinematic eddy viscosity into temp
      call get_vreman(dt=dt,visc=visc_t,U=this%Q,V=this%Q,W=this%Q,Ucomp=1,Vcomp=2,Wcomp=3,Cs=Cs)
      ! Add rho*visc_t to dynamic viscosity
      call this%visc%saxpy(a=this%rho,src=visc_t)
      ! Destroy temp amrdata
      call visc_t%finalize()
   end subroutine add_vreman

   !> Compute CFL numbers (convective and viscous)
   subroutine get_cfl(this,dt,cfl,cflc)
      implicit none
      class(amrcinc), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      real(WP), intent(out), optional :: cflc
      integer :: lvl
      ! Get convective CFL from parent
      call this%amrflow%get_cflc(dt=dt)
      ! Reset viscous CFLs
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
      class(amrcinc), intent(inout) :: this
      integer :: lvl

      ! Use parent's method first
      call this%amrflow%get_info()

      ! Initialize min/max values
      this%Pmax=-huge(1.0_WP)

      ! Loop over all levels for min/max
      do lvl=0,this%amr%clvl()
         this%Umax=max(this%Umax,this%Q%norm0(lvl=lvl,comp=1))
         this%Vmax=max(this%Vmax,this%Q%norm0(lvl=lvl,comp=2))
         this%Wmax=max(this%Wmax,this%Q%norm0(lvl=lvl,comp=3))
         this%Pmax=max(this%Pmax,this%P%norm0(lvl=lvl))
      end do

      ! Momentum integrals (rho * U * dV, summed over cells at level 0)
      this%rhoUint=this%rho*this%Q%get_sum(lvl=0,comp=1)*this%amr%cell_vol(0)
      this%rhoVint=this%rho*this%Q%get_sum(lvl=0,comp=2)*this%amr%cell_vol(0)
      this%rhoWint=this%rho*this%Q%get_sum(lvl=0,comp=3)*this%amr%cell_vol(0)

      ! Kinetic energy integral: 0.5 * rho * (Uc^2 + Vc^2 + Wc^2) * dV
      get_kinetic_energy: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         integer :: i,j,k,ierr
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
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
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over tile
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then
                     if (pMask(i,j,k,1).eq.0) cycle
                  end if
                  ! Accumulate kinetic energy
                  this%rhoKint=this%rhoKint+0.5_WP*this%rho*(pQ(i,j,k,1)**2+pQ(i,j,k,2)**2+pQ(i,j,k,3)**2)*this%amr%cell_vol(lvl)
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      end block get_kinetic_energy

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrcinc_print(this)
      use messager, only: log
      use string, only: str_long
      implicit none
      class(amrcinc), intent(in) :: this
      character(len=str_long) :: message
      call log("Incompressible collocated solver: "//trim(this%name))
      write(message,'("  rho = ",ES12.5)') this%rho
      call log(trim(message))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrcinc_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrcinc), intent(inout) :: this
      class(amrio), intent(inout) :: io
      ! Face velocities and conserved variables are registered with parent
      call this%amrflow%register_checkpoint(io)
      ! Register remaining data
      call io%add_data(this%P,'P')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrcinc), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      ! Restore face velocities and conserved variables via parent
      call this%amrflow%restore_checkpoint(io,dirname,time)
      ! Restore remaining data
      call io%read_data(dirname,this%P,'P'); call this%P%fill(time=time)
   end subroutine restore_checkpoint

end module amrcinc_class
