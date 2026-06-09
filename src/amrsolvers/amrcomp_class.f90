!> AMR Collocated compressible solver class
module amrcomp_class
   use iso_c_binding,    only: c_ptr,c_f_pointer,c_loc,c_f_pointer
   use precision,        only: WP
   use amrdata_class,    only: amrdata
   use amrflow_class,    only: amrflow
   use amrmg_class,      only: amrmg
   use material_class,   only: material
   use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter
   implicit none
   private

   ! Expose type
   public :: amrcomp

   !> AMR collocated compressible solver type
   type, extends(amrflow) :: amrcomp

      ! User-configurable callbacks
      procedure(comp_init_iface),    pointer, pass :: user_init   =>null()  !< User-defined initialization
      procedure(comp_tagging_iface), pointer, pass :: user_tagging=>null()  !< User-defined tagging
      procedure(comp_bc_iface),      pointer, pass :: user_bc     =>null()  !< User-defined boundary conditions

      ! Working material
      class(material), pointer :: mat=>null()

      ! Species index ranges in Q
      integer :: Y_lo=0,Y_hi=-1          !< Species range in  Q(:,:,:,Y_lo:Y_hi)

      ! Pressure solver for pressure projection
      logical :: use_projection=.false.
      type(amrmg) :: psolver

      ! Cell-centered primitive variables (velocities, internal energy, and pressure)
      type(amrdata) :: UVW,I,P

      ! Temperature
      type(amrdata) :: T

      ! Speed of sound
      type(amrdata) :: C

      !> Species mass fractions (ns-1 components; only allocated when ns>1)
      type(amrdata) :: Y

      ! Physical properties
      type(amrdata) :: visc              !< Dynamic viscosity
      type(amrdata) :: beta              !< Bulk viscosity
      type(amrdata) :: diff              !< Heat diffusivity

      ! CFL numbers
      real(WP) :: CFLp=0.0_WP                                !< Pressure+convection
      real(WP) :: CFLa_x=0.0_WP,CFLa_y=0.0_WP,CFLa_z=0.0_WP  !< Acoustic
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous

      ! Monitoring quantities
      real(WP) :: Imin=0.0_WP,Imax=0.0_WP
      real(WP) :: Pmin=0.0_WP,Pmax=0.0_WP
      real(WP) :: Tmin=0.0_WP,Tmax=0.0_WP
      real(WP) :: Cmin=0.0_WP,Cmax=0.0_WP
      real(WP) :: rhoKint=0.0_WP
      real(WP), dimension(:), allocatable :: Ymin,Ymax

      ! Minimum density for stability
      real(WP) :: rho_floor=1.0e-10_WP

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
      procedure :: apply_velbc=>comp_apply_velbc
      procedure :: apply_Qbc=>comp_apply_Qbc
      ! Utilities
      procedure :: get_face_velocity         !< Update face velocities from cell-centered data
      procedure :: add_pressure              !< Add pressure term to face velocities, cell-centered momentum, and internal energy
      procedure :: prepare_psolver           !< Prepare Helmholtz pressure solver
      ! Physics
      procedure :: get_primitive             !< Get primitive variables from conserved variables
      procedure :: get_dQdt                  !< Compute conserved variable time derivative
      procedure :: add_viscartif             !< Add localized artificial diffusivity
      procedure :: add_vreman                !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                   !< Compute CFL numbers
      ! Print
      procedure :: get_info
      procedure :: print=>amrcomp_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrcomp

   !> Abstract interface for user-provided init callback
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

   !> Abstract interface for user-provided tagging callback
   abstract interface
      subroutine comp_tagging_iface(solver,lvl,time,tags)
         import :: amrcomp,c_ptr,WP
         class(amrcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine comp_tagging_iface
   end interface

   !> Abstract interface for user-provided BC callback
   abstract interface
      subroutine comp_bc_iface(solver,lvl,time,face,bx,comp,p)
         use amrex_amr_module, only: amrex_box
         import :: amrcomp,WP
         class(amrcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< Can be 'U','V','W','Q'
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      end subroutine comp_bc_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrcomp_on_init(ctx,lvl,time,ba,dm)
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
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrcomp_on_clear

   !> Dispatch tagging: calls user callback if set
   subroutine amrcomp_tagging(ctx,lvl,time,tags)
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
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrmg_class,      only: amrmg_varcoef
      use amrgrid_class,    only: amrgrid
      use messager,         only: die
      implicit none
      class(amrcomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Material should have been provided by the user before this call
      if (.not.associated(this%mat)) call die('[amrcomp initialize] mat must be assigned before initialize')

      ! Q layout: 5 base components + (ns-1) extra species
      this%nQ=5+(this%mat%ns-1); this%nover=max(this%nover,2)
      this%Y_lo=6; this%Y_hi=4+this%mat%ns
      call this%amrflow%initialize(amr=amr,name=name); call this%set_parent()

      ! Initialize primitive/derived variables
      call this%UVW%initialize(amr,name='UVW',ncomp=3,ng=this%nover); this%UVW%parent=>this
      call this%I%initialize(amr,name='I',ncomp=1,ng=this%nover); this%I%parent=>this
      call this%P%initialize(amr,name='P',ncomp=1,ng=this%nover); this%P%parent=>this
      call this%T%initialize(amr,name='T',ncomp=1,ng=this%nover); this%T%parent=>this
      call this%C%initialize(amr,name='C',ncomp=1,ng=this%nover); this%C%parent=>this

      ! Species mass fractions (only if more than 1 species)
      if (this%mat%ns.gt.1) then
         call this%Y%initialize(amr,name='Y',ncomp=this%mat%ns-1,ng=this%nover); this%Y%parent=>this
         allocate(this%Ymin(this%mat%ns-1),this%Ymax(this%mat%ns-1)); this%Ymin=0.0_WP; this%Ymax=0.0_WP
      end if

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
      call this%UVW%finalize()
      call this%I%finalize()
      call this%P%finalize()
      call this%T%finalize()
      call this%C%finalize()
      call this%visc%finalize()
      call this%beta%finalize()
      call this%diff%finalize()
      if (this%mat%ns.gt.1) call this%Y%finalize()
      if (allocated(this%Ymin)) deallocate(this%Ymin)
      if (allocated(this%Ymax)) deallocate(this%Ymax)
      if (this%use_projection) call this%psolver%finalize()
      this%use_projection=.false.
      nullify(this%user_init)
      nullify(this%user_tagging)
      nullify(this%user_bc)
      nullify(this%mat)
      call this%amrflow%finalize()
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
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_init(lvl,time,ba,dm)
      ! Reset level layouts
      call this%UVW%reset_level(lvl,ba,dm)
      call this%I%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%T%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      if (this%mat%ns.gt.1) call this%Y%reset_level(lvl,ba,dm)
      ! Zero out
      call this%UVW%setval(val=0.0_WP,lvl=lvl)
      call this%I%setval(val=0.0_WP,lvl=lvl)
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%T%setval(val=0.0_WP,lvl=lvl)
      call this%C%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl)
      call this%beta%setval(val=0.0_WP,lvl=lvl)
      call this%diff%setval(val=0.0_WP,lvl=lvl)
      if (this%mat%ns.gt.1) call this%Y%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse using conservative interpolation
   subroutine on_coarse(this,lvl,time,ba,dm)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_coarse(lvl,time,ba,dm)
      ! Derived variables are just reset
      call this%UVW%reset_level(lvl,ba,dm)
      call this%I%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%T%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      if (this%mat%ns.gt.1) call this%Y%reset_level(lvl,ba,dm)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid using conservative interpolation
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_remake(lvl,time,ba,dm)
      ! Derived variables are just reset
      call this%UVW%reset_level(lvl,ba,dm)
      call this%I%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%T%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm)
      call this%beta%reset_level(lvl,ba,dm)
      call this%diff%reset_level(lvl,ba,dm)
      if (this%mat%ns.gt.1) call this%Y%reset_level(lvl,ba,dm)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles face velocities, divergence, and conserved quantities
      call this%amrflow%on_clear(lvl)
      ! Clear derived variables
      call this%UVW%clear_level(lvl)
      call this%I%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%T%clear_level(lvl)
      call this%C%clear_level(lvl)
      call this%visc%clear_level(lvl)
      call this%beta%clear_level(lvl)
      call this%diff%clear_level(lvl)
      if (this%mat%ns.gt.1) call this%Y%clear_level(lvl)
   end subroutine on_clear

   !> Override post_regrid: average down for C/F consistency
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles face velocities and conserved quantities
      call this%amrflow%post_regrid(lbase,time)
      ! Rebuild primitive variables
      call this%get_primitive(this%Q)
   end subroutine post_regrid

   ! ============================================================================
   ! BOUNDARY CONDITIONS
   ! ============================================================================

   !> Velocity BC override: forward to user_bc with U/V/W component name
   subroutine comp_apply_velbc(this,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp=comp,p=p)
   end subroutine comp_apply_velbc

   !> Q BC override: forward to user_bc with comp='Q'
   subroutine comp_apply_Qbc(this,lvl,time,face,bx,p)
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp='Q',p=p)
   end subroutine comp_apply_Qbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Update face velocity from Q
   subroutine get_face_velocity(this)
      implicit none
      class(amrcomp), intent(inout) :: this
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
               pU(i,j,k,1)=sum(pQ(i-1:i,j,k,2))/max(sum(pQ(i-1:i,j,k,1)),this%rho_floor)
            end do; end do; end do
            ! Get Y-face velocity
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pV(i,j,k,1)=sum(pQ(i,j-1:j,k,3))/max(sum(pQ(i,j-1:j,k,1)),this%rho_floor)
            end do; end do; end do
            ! Get Z-face velocity
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pW(i,j,k,1)=sum(pQ(i,j,k-1:k,4))/max(sum(pQ(i,j,k-1:k,1)),this%rho_floor)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_face_velocity

   !> Add pressure term to face velocities and cell-centered momentum and internal energy. Two flavors:
   !>   phi present -> direct path: use explicit stencil that reads phi ghost cells directly (for predictor with fs%P)
   !>   phi absent  -> MLMG path:   use psolver internal fluxes (for projection with dP)
   !> Cell-center correction averages the face gradients back to cell center
   !> Optional mask argument is used for IB masking
   !> Optional gravity(1:3) adds a constant face acceleration alongside -grad(p)/rho if phi is present
   subroutine add_pressure(this,scale,phi,mask,gravity)
      use amrex_amr_module, only: amrex_multifab,amrex_bc_reflect_odd
      use amrex_interface,  only: amrmfab_average_down_face
      use messager, only: die
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: scale
      type(amrdata), intent(in), optional :: phi
      type(amrdata), intent(in), optional :: mask
      real(WP), dimension(3), intent(in), optional :: gravity
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pP,pQ,pU,pV,pW,pMask,pPold
      real(WP) :: dxi,dyi,dzi,coeff,crossterm
      integer :: lvl,i,j,k
      ! Build temp face mfabs to store pressure fluxes
      allocate(Fx(0:this%amr%clvl()),Fy(0:this%amr%clvl()),Fz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Compute -1/rho*pressure gradient at faces
      if (present(phi)) then
         ! Use provided phi and its ghosts cells
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pP=>phi%mf(lvl)%dataptr(mfi)
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               ! Get tilebox
               bx=mfi%nodaltilebox(1)
               ! Compute pressure gradient at faces
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFx(i,j,k,1)=-2.0_WP*(pP(i,j,k,1)-pP(i-1,j,k,1))*dxi/sum(max(pQ(i-1:i,j,k,1),this%rho_floor))
                  if (present(gravity)) pFx(i,j,k,1)=pFx(i,j,k,1)+gravity(1)
               end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFy(i,j,k,1)=-2.0_WP*(pP(i,j,k,1)-pP(i,j-1,k,1))*dyi/sum(max(pQ(i,j-1:j,k,1),this%rho_floor))
                  if (present(gravity)) pFy(i,j,k,1)=pFy(i,j,k,1)+gravity(2)
               end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFz(i,j,k,1)=-2.0_WP*(pP(i,j,k,1)-pP(i,j,k-1,1))*dzi/sum(max(pQ(i,j,k-1:k,1),this%rho_floor))
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
      else
         ! Use psolver's solution and its internal ghosts
         if (this%use_projection) then
            call this%psolver%get_fluxes(Fx,Fy,Fz)
         else
            call die('[amrcomp::add_pressure] use_projection must be true to use internal phi')
         end if
      end if
      ! Apply to face velocities and cell-centered in one pass
      do lvl=0,this%amr%clvl()
         ! Get mesh size
         dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
         ! Face velocities: use flux directly
         call this%U%mf(lvl)%saxpy(scale,Fx(lvl),1,1,1,0)
         call this%V%mf(lvl)%saxpy(scale,Fy(lvl),1,1,1,0)
         call this%W%mf(lvl)%saxpy(scale,Fz(lvl),1,1,1,0)
         ! Momentum gets flux averaged to cell center, internal energy gets -P*div(U_face_updated)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            if (present(phi)) then
               pP=>phi%mf(lvl)%dataptr(mfi)
            else
               ! Use psolver's solution and its internal ghosts
               if (this%use_projection) then
                  pP=>this%psolver%sol%mf(lvl)%dataptr(mfi)
                  pPold=>this%P%mf(lvl)%dataptr(mfi)
               else
                  call die('[amrcomp::add_pressure] use_projection must be true to use internal phi')
               end if
            end if
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            if (present(mask)) pMask=>mask%mf(lvl)%dataptr(mfi)
            ! Get tilebox
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pQ(i,j,k,2)=pQ(i,j,k,2)+scale*pQ(i,j,k,1)*0.5_WP*sum(pFx(i:i+1,j,k,1))
               pQ(i,j,k,3)=pQ(i,j,k,3)+scale*pQ(i,j,k,1)*0.5_WP*sum(pFy(i,j:j+1,k,1))
               pQ(i,j,k,4)=pQ(i,j,k,4)+scale*pQ(i,j,k,1)*0.5_WP*sum(pFz(i,j,k:k+1,1))
               coeff=1.0_WP; if (present(mask)) coeff=pMask(i,j,k,1)
               crossterm=0.0_WP; if (.not.present(phi)) crossterm=-coeff*scale**2*pPold(i,j,k,1)*(dxi*(pFx(i+1,j,k,1)-pFx(i,j,k,1))+dyi*(pFy(i,j+1,k,1)-pFy(i,j,k,1))+dzi*(pFz(i,j,k+1,1)-pFz(i,j,k,1)))
               pQ(i,j,k,5)=pQ(i,j,k,5)-coeff*scale*pP(i,j,k,1)*(dxi*(pU(i+1,j,k,1)-pU(i,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k,1)))+crossterm
               ! Fix non-periodic boundary conditions
               if (.not.this%amr%xper) then
                  if (i.eq.this%amr%geom(lvl)%domain%lo(1).and.this%U%lo_bc(1,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.25_WP*sum(pQ(i:i+1,j,k,1))*pFx(i+1,j,k,1)
                  if (i.eq.this%amr%geom(lvl)%domain%hi(1).and.this%U%hi_bc(1,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,2)=pQ(i,j,k,2)+scale*0.25_WP*sum(pQ(i-1:i,j,k,1))*pFx(i  ,j,k,1)
               end if
               if (.not.this%amr%yper) then
                  if (j.eq.this%amr%geom(lvl)%domain%lo(2).and.this%V%lo_bc(2,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.25_WP*sum(pQ(i,j:j+1,k,1))*pFy(i,j+1,k,1)
                  if (j.eq.this%amr%geom(lvl)%domain%hi(2).and.this%V%hi_bc(2,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,3)=pQ(i,j,k,3)+scale*0.25_WP*sum(pQ(i,j-1:j,k,1))*pFy(i,j  ,k,1)
               end if
               if (.not.this%amr%zper) then
                  if (k.eq.this%amr%geom(lvl)%domain%lo(3).and.this%W%lo_bc(3,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,4)=pQ(i,j,k,4)+scale*0.25_WP*sum(pQ(i,j,k:k+1,1))*pFz(i,j,k+1,1)
                  if (k.eq.this%amr%geom(lvl)%domain%hi(3).and.this%W%hi_bc(3,1).ne.amrex_bc_reflect_odd) pQ(i,j,k,4)=pQ(i,j,k,4)+scale*0.25_WP*sum(pQ(i,j,k-1:k,1))*pFz(i,j,k  ,1)
               end if
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
   end subroutine add_pressure

   !> Prepare variable-coefficient pressure solver using face densities and speed of sound
   subroutine prepare_psolver(this,dt)
      use amrex_amr_module, only: amrex_multifab
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(:), allocatable :: AA,BBx,BBy,BBz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pAA,pBBx,pBBy,pBBz,pQ,pC
      ! Allocate temporary face coefficient mfabs
      allocate(AA(0:this%amr%clvl()),BBx(0:this%amr%clvl()),BBy(0:this%amr%clvl()),BBz(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,AA (lvl),ncomp=1,nover=0,atface=[.false.,.false.,.false.])
         call this%amr%mfab_build(lvl,BBx(lvl),ncomp=1,nover=0,atface=[.true., .false.,.false.])
         call this%amr%mfab_build(lvl,BBy(lvl),ncomp=1,nover=0,atface=[.false.,.true., .false.])
         call this%amr%mfab_build(lvl,BBz(lvl),ncomp=1,nover=0,atface=[.false.,.false.,.true. ])
      end do
      ! Fill Helmholtz coefficients
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ  =>this%Q%mf(lvl)%dataptr(mfi)
            pC  =>this%C%mf(lvl)%dataptr(mfi)
            pAA =>AA (lvl)%dataptr(mfi)
            pBBx=>BBx(lvl)%dataptr(mfi)
            pBBy=>BBy(lvl)%dataptr(mfi)
            pBBz=>BBz(lvl)%dataptr(mfi)
            ! Cell-centered
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pAA(i,j,k,1)=1.0_WP/(pQ(i,j,k,1)*pC(i,j,k,1)**2)
            end do; end do; end do
            ! X-faces
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBBx(i,j,k,1)=2.0_WP/max(sum(pQ(i-1:i,j,k,1)),this%rho_floor)
            end do; end do; end do
            ! Y-faces
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBBy(i,j,k,1)=2.0_WP/max(sum(pQ(i,j-1:j,k,1)),this%rho_floor)
            end do; end do; end do
            ! Z-faces
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBBz(i,j,k,1)=2.0_WP/max(sum(pQ(i,j,k-1:k,1)),this%rho_floor)
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

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Calculate primitive variables from conserved variables
   subroutine get_primitive(this,Q)
      use messager, only: die
      implicit none
      class(amrcomp), intent(inout) :: this
      type(amrdata), intent(in) :: Q
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pUVW,pI,pP,pT,pC,pY
      real(WP) :: irho
      real(WP), dimension(this%mat%ns) :: y
      ! Check passed Q is as expected
      if (Q%ncomp.ne.this%nQ) call die('[amrcomp get_primitive] Q has wrong number of components')
      if (Q%ng.lt.this%nover) call die('[amrcomp get_primitive] Q must have at least nover ghost cells')
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>Q%mf(lvl)%dataptr(mfi)
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            pI=>this%I%mf(lvl)%dataptr(mfi)
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pT=>this%T%mf(lvl)%dataptr(mfi)
            pC=>this%C%mf(lvl)%dataptr(mfi)
            if (this%mat%ns.gt.1) pY=>this%Y%mf(lvl)%dataptr(mfi)
            ! Loop over grown tiles
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Compute velocity from momentum
               irho=1.0_WP/max(pQ(i,j,k,1),this%rho_floor)
               pUVW(i,j,k,1)=pQ(i,j,k,2)*irho
               pUVW(i,j,k,2)=pQ(i,j,k,3)*irho
               pUVW(i,j,k,3)=pQ(i,j,k,4)*irho
               ! Compute internal energy per unit mass
               pI(i,j,k,1)=pQ(i,j,k,5)*irho
               ! Composition: cache first ns-1 species (clipped to [0,1]), close ns-th (clipped to [0,1])
               if (this%mat%ns.gt.1) then
                  pY(i,j,k,:)=max(0.0_WP,min(pQ(i,j,k,this%Y_lo:this%Y_hi)*irho,1.0_WP))
                  y(1:this%mat%ns-1)=pY(i,j,k,:)
               end if
               y(this%mat%ns)=max(0.0_WP,1.0_WP-sum(y(1:this%mat%ns-1)))
               ! Compute pressure via EoS: P = P(rho, I)
               pP(i,j,k,1)=this%mat%get_p_from_rho_e(rho=pQ(i,j,k,1),e=pI(i,j,k,1),y=y)
               ! Compute speed of sound via EoS: C = C(rho, P)
               pC(i,j,k,1)=this%mat%get_c_from_p_rho(p=pP(i,j,k,1),rho=pQ(i,j,k,1),y=y)
               ! Compute temperature via EoS: T = T(rho, P)
               pT(i,j,k,1)=this%mat%get_T_from_p_rho(p=pP(i,j,k,1),rho=pQ(i,j,k,1),y=y)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
   end subroutine get_primitive

   !> Calculate dQdt from passed Q without pressure term (user can add it via add_pressure)
   subroutine get_dQdt(this,dQdt)
      use amrex_amr_module, only: amrex_multifab
      implicit none
      class(amrcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: dQdt
      type(amrex_multifab), dimension(0:this%amr%maxlvl) :: Fx,Fy,Fz

      ! Initialize all fluxes
      define_fluxes: block
         integer :: lvl
         do lvl=0,this%amr%clvl()
            call this%amr%mfab_build(lvl,Fx(lvl),ncomp=this%nQ,nover=1,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl,Fy(lvl),ncomp=this%nQ,nover=1,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
            call this%amr%mfab_build(lvl,Fz(lvl),ncomp=this%nQ,nover=1,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
         end do
      end block define_fluxes
      
      ! Compute fluxes for all levels
      compute_fluxes: block
         integer :: lvl,i,j,k,n
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx
         real(WP) :: dxi,dyi,dzi,div,w
         real(WP), dimension(-2: 0) :: wenop
         real(WP), dimension(-1:+1) :: wenom
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pQ,pUVW,pI,pT,pVisc,pBeta,pDiff,pY
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz
         real(WP), parameter :: eps=1.0e-15_WP
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            ! Loop over all tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               pQ=>this%Q%mf(lvl)%dataptr(mfi)
               pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
               pI=>this%I%mf(lvl)%dataptr(mfi)
               pT=>this%T%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pDiff=>this%diff%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi)
               pFy=>Fy(lvl)%dataptr(mfi)
               pFz=>Fz(lvl)%dataptr(mfi)
               if (this%mat%ns.gt.1) pY=>this%Y%mf(lvl)%dataptr(mfi)
               ! X-fluxes
               fbx=mfi%nodaltilebox(1)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! WENO mass flux
                  w=weno_weight((abs(pQ(i-1,j,k,1)-pQ(i-2,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i+1,j,k,1)-pQ(i  ,j,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFx(i,j,k,1)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wenop*pQ(i-2:i  ,j,k,1)) &
                  &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wenom*pQ(i-1:i+1,j,k,1))
                  ! Momentum fluxes
                  pFx(i,j,k,2)=pFx(i,j,k,1)*0.5_WP*sum(pUVW(i-1:i,j,k,1))
                  pFx(i,j,k,3)=pFx(i,j,k,1)*0.5_WP*sum(pUVW(i-1:i,j,k,2))
                  pFx(i,j,k,4)=pFx(i,j,k,1)*0.5_WP*sum(pUVW(i-1:i,j,k,3))
                  ! WENO internal energy flux
                  w=weno_weight((abs(pI(i-1,j,k,1)-pI(i-2,j,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i-1,j,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pI(i+1,j,k,1)-pI(i  ,j,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFx(i,j,k,5)=0.5_WP*(pFx(i,j,k,1)-abs(pFx(i,j,k,1)))*sum(wenop*pI(i-2:i  ,j,k,1)) &
                  &           +0.5_WP*(pFx(i,j,k,1)+abs(pFx(i,j,k,1)))*sum(wenom*pI(i-1:i+1,j,k,1))
                  ! WENO species fluxes
                  do n=1,this%mat%ns-1
                     w=weno_weight((abs(pY(i-1,j,k,n)-pY(i-2,j,k,n))+eps)/(abs(pY(i,j,k,n)-pY(i-1,j,k,n))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pY(i+1,j,k,n)-pY(i  ,j,k,n))+eps)/(abs(pY(i,j,k,n)-pY(i-1,j,k,n))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                     pFx(i,j,k,this%Y_lo+n-1)=0.5_WP*(pFx(i,j,k,1)-abs(pFx(i,j,k,1)))*sum(wenop*pY(i-2:i  ,j,k,n)) &
                     &                       +0.5_WP*(pFx(i,j,k,1)+abs(pFx(i,j,k,1)))*sum(wenom*pY(i-1:i+1,j,k,n))
                  end do
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
                  ! Viscous stress at x-face
                  pFx(i,j,k,2)=pFx(i,j,k,2)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(1,1)+gradU(1,1))+0.5_WP*(sum(pBeta(i-1:i,j,k,1))-2.0_WP/3.0_WP*sum(pVisc(i-1:i,j,k,1)))*div
                  pFx(i,j,k,3)=pFx(i,j,k,3)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(2,1)+gradU(1,2))
                  pFx(i,j,k,4)=pFx(i,j,k,4)+0.5_WP*sum(pVisc(i-1:i,j,k,1))*(gradU(3,1)+gradU(1,3))
                  ! Heat diffusion flux
                  pFx(i,j,k,5)=pFx(i,j,k,5)+0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pT(i,j,k,1)-pT(i-1,j,k,1))
                  ! Species diffusion flux (Le=1)
                  do n=1,this%mat%ns-1
                     pFx(i,j,k,this%Y_lo+n-1)=pFx(i,j,k,this%Y_lo+n-1)+0.5_WP*sum(pDiff(i-1:i,j,k,1))*dxi*(pY(i,j,k,n)-pY(i-1,j,k,n))
                  end do
               end do; end do; end do
               ! Y-fluxes
               fbx=mfi%nodaltilebox(2)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! WENO mass flux
                  w=weno_weight((abs(pQ(i,j-1,k,1)-pQ(i,j-2,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i,j+1,k,1)-pQ(i,j  ,k,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFy(i,j,k,1)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wenop*pQ(i,j-2:j  ,k,1)) &
                  &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wenom*pQ(i,j-1:j+1,k,1))
                  ! Momentum fluxes
                  pFy(i,j,k,2)=pFy(i,j,k,1)*0.5_WP*sum(pUVW(i,j-1:j,k,1))
                  pFy(i,j,k,3)=pFy(i,j,k,1)*0.5_WP*sum(pUVW(i,j-1:j,k,2))
                  pFy(i,j,k,4)=pFy(i,j,k,1)*0.5_WP*sum(pUVW(i,j-1:j,k,3))
                  ! WENO internal energy flux
                  w=weno_weight((abs(pI(i,j-1,k,1)-pI(i,j-2,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j-1,k,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pI(i,j+1,k,1)-pI(i,j  ,k,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFy(i,j,k,5)=0.5_WP*(pFy(i,j,k,1)-abs(pFy(i,j,k,1)))*sum(wenop*pI(i,j-2:j  ,k,1)) &
                  &           +0.5_WP*(pFy(i,j,k,1)+abs(pFy(i,j,k,1)))*sum(wenom*pI(i,j-1:j+1,k,1))
                  ! WENO species fluxes
                  do n=1,this%mat%ns-1
                     w=weno_weight((abs(pY(i,j-1,k,n)-pY(i,j-2,k,n))+eps)/(abs(pY(i,j,k,n)-pY(i,j-1,k,n))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pY(i,j+1,k,n)-pY(i,j  ,k,n))+eps)/(abs(pY(i,j,k,n)-pY(i,j-1,k,n))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                     pFy(i,j,k,this%Y_lo+n-1)=0.5_WP*(pFy(i,j,k,1)-abs(pFy(i,j,k,1)))*sum(wenop*pY(i,j-2:j  ,k,n)) &
                     &                       +0.5_WP*(pFy(i,j,k,1)+abs(pFy(i,j,k,1)))*sum(wenom*pY(i,j-1:j+1,k,n))
                  end do
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
                  ! Viscous stress at y-face
                  pFy(i,j,k,2)=pFy(i,j,k,2)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(1,2)+gradU(2,1))
                  pFy(i,j,k,3)=pFy(i,j,k,3)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(2,2)+gradU(2,2))+0.5_WP*(sum(pBeta(i,j-1:j,k,1))-2.0_WP/3.0_WP*sum(pVisc(i,j-1:j,k,1)))*div
                  pFy(i,j,k,4)=pFy(i,j,k,4)+0.5_WP*sum(pVisc(i,j-1:j,k,1))*(gradU(3,2)+gradU(2,3))
                  ! Heat diffusion flux
                  pFy(i,j,k,5)=pFy(i,j,k,5)+0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pT(i,j,k,1)-pT(i,j-1,k,1))
                  ! Species diffusion flux (Le=1)
                  do n=1,this%mat%ns-1
                     pFy(i,j,k,this%Y_lo+n-1)=pFy(i,j,k,this%Y_lo+n-1)+0.5_WP*sum(pDiff(i,j-1:j,k,1))*dyi*(pY(i,j,k,n)-pY(i,j-1,k,n))
                  end do
               end do; end do; end do
               ! Z-fluxes
               fbx=mfi%nodaltilebox(3)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! WENO mass flux
                  w=weno_weight((abs(pQ(i,j,k-1,1)-pQ(i,j,k-2,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pQ(i,j,k+1,1)-pQ(i,j,k  ,1))+eps)/(abs(pQ(i,j,k,1)-pQ(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFz(i,j,k,1)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wenop*pQ(i,j,k-2:k,1)) &
                  &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wenom*pQ(i,j,k-1:k+1,1))
                  ! Momentum fluxes
                  pFz(i,j,k,2)=pFz(i,j,k,1)*0.5_WP*sum(pUVW(i,j,k-1:k,1))
                  pFz(i,j,k,3)=pFz(i,j,k,1)*0.5_WP*sum(pUVW(i,j,k-1:k,2))
                  pFz(i,j,k,4)=pFz(i,j,k,1)*0.5_WP*sum(pUVW(i,j,k-1:k,3))
                  ! WENO internal energy flux
                  w=weno_weight((abs(pI(i,j,k-1,1)-pI(i,j,k-2,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j,k-1,1))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                  w=weno_weight((abs(pI(i,j,k+1,1)-pI(i,j,k  ,1))+eps)/(abs(pI(i,j,k,1)-pI(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                  pFz(i,j,k,5)=0.5_WP*(pFz(i,j,k,1)-abs(pFz(i,j,k,1)))*sum(wenop*pI(i,j,k-2:k  ,1)) &
                  &           +0.5_WP*(pFz(i,j,k,1)+abs(pFz(i,j,k,1)))*sum(wenom*pI(i,j,k-1:k+1,1))
                  ! WENO species fluxes
                  do n=1,this%mat%ns-1
                     w=weno_weight((abs(pY(i,j,k-1,n)-pY(i,j,k-2,n))+eps)/(abs(pY(i,j,k,n)-pY(i,j,k-1,n))+eps)); wenop=0.5_WP*[-w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pY(i,j,k+1,n)-pY(i,j,k  ,n))+eps)/(abs(pY(i,j,k,n)-pY(i,j,k-1,n))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,-w]
                     pFz(i,j,k,this%Y_lo+n-1)=0.5_WP*(pFz(i,j,k,1)-abs(pFz(i,j,k,1)))*sum(wenop*pY(i,j,k-2:k  ,n)) &
                     &                       +0.5_WP*(pFz(i,j,k,1)+abs(pFz(i,j,k,1)))*sum(wenom*pY(i,j,k-1:k+1,n))
                  end do
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
                  ! Viscous stress at z-face
                  pFz(i,j,k,2)=pFz(i,j,k,2)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(1,3)+gradU(3,1))
                  pFz(i,j,k,3)=pFz(i,j,k,3)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(2,3)+gradU(3,2))
                  pFz(i,j,k,4)=pFz(i,j,k,4)+0.5_WP*sum(pVisc(i,j,k-1:k,1))*(gradU(3,3)+gradU(3,3))+0.5_WP*(sum(pBeta(i,j,k-1:k,1))-2.0_WP/3.0_WP*sum(pVisc(i,j,k-1:k,1)))*div
                  ! Heat diffusion flux
                  pFz(i,j,k,5)=pFz(i,j,k,5)+0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pT(i,j,k,1)-pT(i,j,k-1,1))
                  ! Species diffusion flux (Le=1)
                  do n=1,this%mat%ns-1
                     pFz(i,j,k,this%Y_lo+n-1)=pFz(i,j,k,this%Y_lo+n-1)+0.5_WP*sum(pDiff(i,j,k-1:k,1))*dzi*(pY(i,j,k,n)-pY(i,j,k-1,n))
                  end do
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
         integer :: lvl,i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: dxi,dyi,dzi,div
         real(WP), dimension(1:3,1:3) :: gradU
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW,pVisc,pBeta
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pdQdt,pFx,pFy,pFz
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            ! Loop over tiles
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pdQdt=>dQdt%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi)
               pFy=>Fy(lvl)%dataptr(mfi)
               pFz=>Fz(lvl)%dataptr(mfi)
               pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               ! Loop over interior
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Advection
                  pdQdt(i,j,k,:)=dxi*(pFx(i+1,j,k,:)-pFx(i,j,k,:))+dyi*(pFy(i,j+1,k,:)-pFy(i,j,k,:))+dzi*(pFz(i,j,k+1,:)-pFz(i,j,k,:))
                  ! Viscous heating: compute cell-centered gradU and stress tensor
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
                  ! τ:∇U = τ_ij * gradU(i,j)
                  pdQdt(i,j,k,5)=pdQdt(i,j,k,5) &
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

   !> Add artificial bulk viscosity to this%beta (and optionally this%visc)
   subroutine add_viscartif(this,dt,Cartif,Cvisc)
      use amrsgs, only: get_viscartif
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cartif
      real(WP), intent(in), optional :: Cvisc
      real(WP) :: myCvisc
      type(amrdata) :: beta_t
      ! Set shear viscosity constant
      if (present(Cvisc)) then; myCvisc=Cvisc; else; myCvisc=0.0_WP; end if
      ! Create temp amrdata
      call beta_t%initialize(this%amr,name='beta_t',ncomp=1,ng=this%nover); call beta_t%reset()
      ! Compute kinematic artificial bulk viscosity into temp
      call get_viscartif(dt=dt,visc=beta_t,U=this%UVW,V=this%UVW,W=this%UVW,Ucomp=1,Vcomp=2,Wcomp=3,C=this%C,Cartif=Cartif)
      ! Add rho*visc_t to dynamic bulk viscosity
      call beta_t%multiply(this%Q,srccomp=1,ncomp=1); call this%beta%add(beta_t)
      ! Add Cvisc*rho*visc_t to dynamic shear viscosity
      if (myCvisc.ne.0.0_WP) call this%visc%saxpy(a=myCvisc,src=beta_t)
      ! Destroy temp amrdata
      call beta_t%finalize()
   end subroutine add_viscartif

   !> Add Vreman SGS eddy viscosity to this%visc
   subroutine add_vreman(this,dt,Cs)
      use amrsgs, only: get_vreman
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrdata) :: visc_t
      ! Create temp amrdata
      call visc_t%initialize(this%amr,name='visc_t',ncomp=1,ng=this%nover); call visc_t%reset()
      ! Compute kinematic eddy viscosity into temp
      call get_vreman(dt=dt,visc=visc_t,U=this%UVW,V=this%UVW,W=this%UVW,Ucomp=1,Vcomp=2,Wcomp=3,Cs=Cs)
      ! Add rho*visc_t to dynamic viscosity
      call visc_t%multiply(src=this%Q,srccomp=1,ncomp=1); call this%visc%add(src=visc_t)
      ! Destroy temp amrdata
      call visc_t%finalize()
   end subroutine add_vreman

   !> Calculate CFL numbers
   subroutine get_cfl(this,dt,cfl)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
      implicit none
      class(amrcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: lvl,i,j,k,ierr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pP,pUVW,pVisc,pBeta,pDiff,pT,pC,pY
      real(WP) :: dxi,dyi,dzi,rho,conv,pgrad,viscmax,cv,alpha_heat
      real(WP), dimension(this%mat%ns) :: y
      ! Get convective CFL from parent
      call this%amrflow%get_cflc(dt=dt)
      ! Reset CFLs
      this%CFLp=0.0_WP
      this%CFLa_x=0.0_WP; this%CFLa_y=0.0_WP; this%CFLa_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         ! Get mesh spacing
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get data pointers
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pDiff=>this%diff%mf(lvl)%dataptr(mfi)
            pT=>this%T%mf(lvl)%dataptr(mfi)
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            pC=>this%C%mf(lvl)%dataptr(mfi)
            if (this%mat%ns.gt.1) pY=>this%Y%mf(lvl)%dataptr(mfi)
            ! Loop over interior tiles
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rho=max(pQ(i,j,k,1),this%rho_floor)
               ! Heat-diffusion CFL: thermal diffusivity alpha=lambda/(rho*cv)
               if (this%mat%ns.gt.1) y(1:this%mat%ns-1)=pY(i,j,k,:)
               y(this%mat%ns)=max(0.0_WP,1.0_WP-sum(y(1:this%mat%ns-1)))
               cv=this%mat%get_cv_from_rho_T(rho,pT(i,j,k,1),y)
               alpha_heat=pDiff(i,j,k,1)/max(rho*cv,tiny(1.0_WP))
               ! Viscous-like CFL
               viscmax=max(pVisc(i,j,k,1)/rho,pBeta(i,j,k,1)/rho,alpha_heat)
               if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt*dxi**2)
               if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt*dyi**2)
               if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt*dzi**2)
               ! Convective+pressure CFL
               conv=abs(pUVW(i,j,k,1))*dxi+abs(pUVW(i,j,k,2))*dyi+abs(pUVW(i,j,k,3))*dzi
               pgrad=(abs(pP(i+1,j,k,1)-pP(i-1,j,k,1))*0.5_WP*dxi**2+abs(pP(i,j+1,k,1)-pP(i,j-1,k,1))*0.5_WP*dyi**2+abs(pP(i,j,k+1,1)-pP(i,j,k-1,1))*0.5_WP*dzi**2)/rho
               this%CFLp=max(this%CFLp,0.5_WP*dt*(conv+sqrt(conv**2+4.0_WP*pgrad)))
               ! Acoustic CFL
               if (this%amr%nx.gt.1) this%CFLa_x=max(this%CFLa_x,(abs(pUVW(i,j,k,1))+pC(i,j,k,1))*dt*dxi)
               if (this%amr%ny.gt.1) this%CFLa_y=max(this%CFLa_y,(abs(pUVW(i,j,k,2))+pC(i,j,k,1))*dt*dyi)
               if (this%amr%nz.gt.1) this%CFLa_z=max(this%CFLa_z,(abs(pUVW(i,j,k,3))+pC(i,j,k,1))*dt*dzi)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Postprocess CFLs
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_x,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_y,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_z,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      end do
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
      if (.not.this%use_projection) cfl=max(cfl,this%CFLa_x,this%CFLa_y,this%CFLa_z)
   end subroutine get_cfl


   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Calculate monitoring info
   subroutine get_info(this)
      implicit none
      class(amrcomp), intent(inout) :: this
      integer :: lvl,n

      ! Use parent's method first
      call this%amrflow%get_info()

      ! Compute extrema across all levels
      this%Imin=huge(1.0_WP); this%Imax=-huge(1.0_WP)
      this%Pmin=huge(1.0_WP); this%Pmax=-huge(1.0_WP)
      this%Tmin=huge(1.0_WP); this%Tmax=-huge(1.0_WP)
      this%Cmin=huge(1.0_WP); this%Cmax=-huge(1.0_WP)
      if (this%mat%ns.gt.1) then; this%Ymin=huge(1.0_WP); this%Ymax=-huge(1.0_WP); end if
      do lvl=0,this%amr%clvl()
         ! Velocity norm 0
         this%Umax=max(this%Umax,this%UVW%norm0(lvl=lvl,comp=1))
         this%Vmax=max(this%Vmax,this%UVW%norm0(lvl=lvl,comp=2))
         this%Wmax=max(this%Wmax,this%UVW%norm0(lvl=lvl,comp=3))
         ! Extrema of internal energy, pressure, and temperature
         this%Imin=min(this%Imin,this%I%get_min(lvl=lvl)); this%Imax=max(this%Imax,this%I%get_max(lvl=lvl))
         this%Pmin=min(this%Pmin,this%P%get_min(lvl=lvl)); this%Pmax=max(this%Pmax,this%P%get_max(lvl=lvl))
         this%Tmin=min(this%Tmin,this%T%get_min(lvl=lvl)); this%Tmax=max(this%Tmax,this%T%get_max(lvl=lvl))
         ! Extrema of speed of sound
         this%Cmin=min(this%Cmin,this%C%get_min(lvl=lvl)); this%Cmax=max(this%Cmax,this%C%get_max(lvl=lvl))
         ! Per-species extrema
         do n=1,this%mat%ns-1
            this%Ymin(n)=min(this%Ymin(n),this%Y%get_min(lvl=lvl,comp=n))
            this%Ymax(n)=max(this%Ymax(n),this%Y%get_max(lvl=lvl,comp=n))
         end do
      end do

      ! Kinetic energy integral: 0.5 * rho * (U^2 + V^2 + W^2) * dV
      ! Uses composite integration with fine masking to avoid double-counting
      get_rhoKint: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface, only: amrmask_make_fine
         use parallel, only: MPI_REAL_WP
         use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         integer :: i,j,k,ierr
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pUVW
         integer, dimension(:,:,:,:), contiguous, pointer :: pMask
         ! Uses composite integration with fine masking to avoid double-counting
         this%rhoKint=0.0_WP
         do lvl=0,this%amr%clvl()
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
                  ! Accumulate kinetic energy
                  this%rhoKint=this%rhoKint+0.5_WP*pQ(i,j,k,1)*(pUVW(i,j,k,1)**2+pUVW(i,j,k,2)**2+pUVW(i,j,k,3)**2)*this%amr%cell_vol(lvl)
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

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      ! Face velocities and conserved variables are registered with parent
      call this%amrflow%register_checkpoint(io)
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      ! Restore face velocities and conserved variables via parent
      call this%amrflow%restore_checkpoint(io,dirname,time)
      ! Rebuild primitive variables
      call this%get_primitive(this%Q)
   end subroutine restore_checkpoint

end module amrcomp_class
