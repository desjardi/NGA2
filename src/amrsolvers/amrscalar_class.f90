!> Constant density scalar solver using AMR
!> Extends amrsolver base class
module amrscalar_class
   use iso_c_binding,    only: c_ptr, c_null_ptr, c_loc, c_f_pointer
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrflux_class,    only: amrflux
   use amrsolver_class,  only: amrsolver
   use amrio_class,      only: amrio
   use amrex_amr_module, only: amrex_multifab, amrex_boxarray, amrex_distromap, &
   &                           amrex_mfiter, amrex_box, amrex_fab
   implicit none
   private

   ! Expose type
   public :: amrscalar

   !> Constant density scalar solver object definition
   type, extends(amrsolver) :: amrscalar

      ! User-configurable callbacks
      procedure(scalar_init_iface), pointer, pass :: user_init => null()
      procedure(scalar_tagging_iface), pointer, pass :: user_tagging => null()

      ! Scalar variable definition
      integer :: nscalar
      character(len=str_medium), dimension(:), allocatable :: SCname

      ! Scalar data containers (parent pointer set to this solver)
      type(amrdata) :: SC
      type(amrdata) :: SCold

      ! Flux register for interlevel conservation (used when use_refluxing=.true.)
      type(amrflux) :: flux

      ! Conservation strategy: .true. = FluxRegister, .false. = flux averaging (default)
      logical :: use_refluxing = .false.

      ! Overlap size for advection stencil
      integer :: nover = 2

      ! Monitoring quantities
      real(WP), dimension(:), allocatable :: SCmax, SCmin, SCint

   contains
      ! Implement deferred procedures from amrsolver
      procedure :: initialize
      procedure :: finalize
      procedure :: get_info
      procedure :: register_checkpoint
      procedure :: restore_checkpoint

      ! Override internal type-bound callbacks from amrsolver
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid

      ! Physics procedures (all-level API)
      procedure :: get_dSCdt
      procedure :: reflux           ! Needs to be called before amerage_down if use_refluxing=.true.
      procedure :: average_down
   end type amrscalar

   !> Abstract interface for user-overridable on_init callback
   abstract interface
      subroutine scalar_init_iface(solver, lvl, time, ba, dm)
         import :: amrscalar, WP, amrex_boxarray, amrex_distromap
         class(amrscalar), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine scalar_init_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine scalar_tagging_iface(solver, lvl, time, tags)
         import :: amrscalar, c_ptr, WP
         class(amrscalar), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine scalar_tagging_iface
   end interface

contains

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrscalar type
   ! ============================================================================

   !> Dispatch on_init: calls user's procedure pointer with typed solver
   subroutine amrscalar_on_init(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_init(lvl, time, ba, dm)
      ! User-provided initialization
      if (associated(this%user_init)) call this%user_init(lvl, time, ba, dm)
   end subroutine amrscalar_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrscalar_on_coarse(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_coarse(lvl, time, ba, dm)
   end subroutine amrscalar_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrscalar_on_remake(ctx, lvl, time, ba, dm)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_remake(lvl, time, ba, dm)
   end subroutine amrscalar_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrscalar_on_clear(ctx, lvl)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx, this)
      call this%on_clear(lvl)
   end subroutine amrscalar_on_clear

   !> Dispatch tagging: calls user's procedure pointer with typed solver
   subroutine amrscalar_tagging(ctx, lvl, time, tags)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx, this)
      ! User-provided tagging
      if (associated(this%user_tagging)) call this%user_tagging(lvl, time, tags)
   end subroutine amrscalar_tagging

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrscalar_postregrid(ctx, lbase, time)
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrscalar), pointer :: this
      call c_f_pointer(ctx, this)
      call this%post_regrid(lbase, time)
   end subroutine amrscalar_postregrid

   ! ============================================================================
   ! LIFECYCLE METHODS
   ! ============================================================================

   !> Initialization for amrscalar solver
   subroutine initialize(this, amr, nscalar, name, use_refluxing)
      use messager, only: die
      class(amrscalar), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: nscalar
      character(len=*), intent(in), optional :: name
      logical, intent(in), optional :: use_refluxing

      ! Set the name for the solver
      if (present(name)) this%name = trim(adjustl(name))

      ! Set refluxing mode
      if (present(use_refluxing)) this%use_refluxing = use_refluxing

      ! Point to amrgrid object (stored in base class)
      this%amr => amr

      ! Set the number of scalars
      this%nscalar = nscalar
      if (this%nscalar .le. 0) call die('[amrscalar initialize] At least 1 scalar is required')

      ! Initialize scalar names
      allocate(this%SCname(1:this%nscalar))
      this%SCname = ''

      ! Initialize info storage
      allocate(this%SCmin(1:this%nscalar)); this%SCmin = +huge(1.0_WP)
      allocate(this%SCmax(1:this%nscalar)); this%SCmax = -huge(1.0_WP)
      allocate(this%SCint(1:this%nscalar)); this%SCint = 0.0_WP

      ! Initialize SC and SCold using amrdata%initialize
      call this%SC%initialize(amr, name='SC', ncomp=this%nscalar, ng=0)
      call this%SCold%initialize(amr, name='SCold', ncomp=this%nscalar, ng=0)

      ! Set parent pointers for coupled BC access
      this%SC%parent => this
      this%SCold%parent => this

      ! Initialize flux register (only used if use_refluxing=.true.)
      if (this%use_refluxing) call this%flux%initialize(amr%maxlvl, name='SC_flux', ncomp=this%nscalar)

      ! Register all 6 callbacks with amrgrid using concrete dispatchers
      ! Use select type to get non-polymorphic target for c_loc
      select type (this)
       type is (amrscalar)
         call this%amr%add_on_init   (amrscalar_on_init,    c_loc(this))
         call this%amr%add_on_coarse (amrscalar_on_coarse,  c_loc(this))
         call this%amr%add_on_remake (amrscalar_on_remake,  c_loc(this))
         call this%amr%add_on_clear  (amrscalar_on_clear,   c_loc(this))
         call this%amr%add_tagging   (amrscalar_tagging,    c_loc(this))
         call this%amr%add_postregrid(amrscalar_postregrid, c_loc(this))
      end select

   end subroutine initialize


   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this, lvl, time, ba, dm)
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Reset level layouts
      call this%SC%reset_level(lvl, ba, dm)
      call this%SCold%reset_level(lvl, ba, dm)
      ! Set to zero
      call this%SC%setval(val=0.0_WP, lvl=lvl)
      call this%SCold%setval(val=0.0_WP, lvl=lvl)
      ! Reset flux register for fine levels (if using refluxing)
      if (this%use_refluxing .and. lvl .ge. 1) call this%flux%reset_level(lvl, ba, dm, [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)])
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   subroutine on_coarse(this, lvl, time, ba, dm)
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! SC gets made from coarse
      call this%SC%on_coarse(lvl, time, ba, dm)
      ! SCold just needs geometry
      call this%SCold%reset_level(lvl, ba, dm)
      ! Reset flux register (if using refluxing)
      if (this%use_refluxing .and. lvl .ge. 1) call this%flux%reset_level(lvl, ba, dm, [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)])
   end subroutine on_coarse


   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this, lvl, time, ba, dm)
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Delegate to SC's on_remake callback
      call this%SC%on_remake(lvl, time, ba, dm)
      ! SCold just needs new geometry
      call this%SCold%reset_level(lvl, ba, dm)
      ! Rebuild flux register for fine levels (if using refluxing)
      if (this%use_refluxing .and. lvl .ge. 1) call this%flux%reset_level(lvl, ba, dm, [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)])
   end subroutine on_remake


   !> Override on_clear: delete level
   subroutine on_clear(this, lvl)
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      call this%SC%clear_level(lvl)
      call this%SCold%clear_level(lvl)
      if (this%use_refluxing .and. lvl .ge. 1)  call this%flux%clear_level(lvl)
   end subroutine on_clear


   !> Override post_regrid: average down for consistency
   subroutine post_regrid(this, lbase, time)
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      integer :: lvl
      ! Average down from finest level to lbase
      do lvl = this%amr%clvl()-1, lbase, -1
         call this%SC%average_downto(lvl)
      end do
   end subroutine post_regrid


   !> Finalization for amrscalar solver
   subroutine finalize(this)
      class(amrscalar), intent(inout) :: this
      ! Destroy containers
      call this%SC%finalize()
      call this%SCold%finalize()
      if (this%use_refluxing) call this%flux%finalize()
      ! Deallocate arrays
      if (allocated(this%SCmax)) deallocate(this%SCmax)
      if (allocated(this%SCmin)) deallocate(this%SCmin)
      if (allocated(this%SCint)) deallocate(this%SCint)
      if (allocated(this%SCname)) deallocate(this%SCname)
      ! Reset pointers
      nullify(this%amr)
      nullify(this%user_init)
      nullify(this%user_tagging)
   end subroutine finalize


   !> Calculate various information on our amrscalar object
   subroutine get_info(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl, nsc
      this%SCmin = +huge(1.0_WP)
      this%SCmax = -huge(1.0_WP)
      this%SCint = 0.0_WP
      do nsc = 1, this%nscalar
         do lvl = 0, this%amr%clvl()
            this%SCmin(nsc) = min(this%SCmin(nsc), this%SC%get_min(lvl=lvl, comp=nsc))
            this%SCmax(nsc) = max(this%SCmax(nsc), this%SC%get_max(lvl=lvl, comp=nsc))
         end do
         this%SCint(nsc) = this%SC%get_sum(lvl=0, comp=nsc) * this%amr%cell_vol(0) / this%amr%vol
      end do
   end subroutine get_info


   !> Calculate dSC/dt for all levels (all-level API)
   !> Uses flux averaging if use_refluxing=.false., FluxRegister if .true.
   subroutine get_dSCdt(this, U, V, W, SC, dSCdt)
      use amrex_interface, only: amrmfab_average_down_faces
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrdata), intent(in) :: U, V, W        ! Face-centered velocity
      class(amrdata), intent(inout) :: SC          ! Cell-centered scalar state to differentiate
      class(amrdata), intent(inout) :: dSCdt       ! Cell-centered output
      ! Local variables
      type(amrex_multifab), dimension(3,0:this%amr%maxlvl) :: flx  ! Temp flux storage
      type(amrex_multifab) :: SCfill
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx, fbx
      type(amrex_fab) :: flux_fab(3)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rhs, pSC, FX, FY, FZ, pU, pV, pW, pFlx
      integer :: lvl, i, j, k, nsc

      ! Phase 1: Compute fluxes for all levels
      do lvl = 0, this%amr%clvl()

         ! Build temp flux MultiFabs for this level
         call this%amr%mfab_build(lvl=lvl, mfab=flx(1,lvl), ncomp=this%nscalar, nover=0, atface=[.true., .false., .false.])
         call this%amr%mfab_build(lvl=lvl, mfab=flx(2,lvl), ncomp=this%nscalar, nover=0, atface=[.false., .true., .false.])
         call this%amr%mfab_build(lvl=lvl, mfab=flx(3,lvl), ncomp=this%nscalar, nover=0, atface=[.false., .false., .true.])

         ! Build SCfill with ghost cells and fill from provided SC
         call this%amr%mfab_build(lvl=lvl, mfab=SCfill, ncomp=this%nscalar, nover=this%nover, atface=[.false., .false., .false.])
         call SC%fill_mfab(SCfill, lvl, 0.0_WP)

         ! Loop over boxes/tiles
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            pSC => SCfill%dataptr(mfi)
            pU => U%mf(lvl)%dataptr(mfi)
            pV => V%mf(lvl)%dataptr(mfi)
            pW => W%mf(lvl)%dataptr(mfi)

            ! Build temporary FABs sized to this tile's nodal boxes
            fbx = bx; call fbx%nodalize(1); call flux_fab(1)%resize(fbx, this%nscalar); FX => flux_fab(1)%dataptr()
            fbx = bx; call fbx%nodalize(2); call flux_fab(2)%resize(fbx, this%nscalar); FY => flux_fab(2)%dataptr()
            fbx = bx; call fbx%nodalize(3); call flux_fab(3)%resize(fbx, this%nscalar); FZ => flux_fab(3)%dataptr()

            do nsc = 1, this%nscalar
               ! X fluxes
               do k = bx%lo(3), bx%hi(3)
                  do j = bx%lo(2), bx%hi(2)
                     do i = bx%lo(1), bx%hi(1)+1
                        FX(i,j,k,nsc) = -0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i-2,j,k,nsc)+5.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+2.0_WP/6.0_WP*pSC(i  ,j,k,nsc)) &
                        &               -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i-1,j,k,nsc)+5.0_WP/6.0_WP*pSC(i  ,j,k,nsc)-1.0_WP/6.0_WP*pSC(i+1,j,k,nsc))
                     end do
                  end do
               end do
               ! Y fluxes
               do k = bx%lo(3), bx%hi(3)
                  do j = bx%lo(2), bx%hi(2)+1
                     do i = bx%lo(1), bx%hi(1)
                        FY(i,j,k,nsc) = -0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j-2,k,nsc)+5.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+2.0_WP/6.0_WP*pSC(i,j  ,k,nsc)) &
                        &               -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j-1,k,nsc)+5.0_WP/6.0_WP*pSC(i,j  ,k,nsc)-1.0_WP/6.0_WP*pSC(i,j+1,k,nsc))
                     end do
                  end do
               end do
               ! Z fluxes
               do k = bx%lo(3), bx%hi(3)+1
                  do j = bx%lo(2), bx%hi(2)
                     do i = bx%lo(1), bx%hi(1)
                        FZ(i,j,k,nsc) = -0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*(-1.0_WP/6.0_WP*pSC(i,j,k-2,nsc)+5.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+2.0_WP/6.0_WP*pSC(i,j,k  ,nsc)) &
                        &               -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*(+2.0_WP/6.0_WP*pSC(i,j,k-1,nsc)+5.0_WP/6.0_WP*pSC(i,j,k  ,nsc)-1.0_WP/6.0_WP*pSC(i,j,k+1,nsc))
                     end do
                  end do
               end do
            end do

            ! Copy flux_fabs into flx mfab
            fbx = mfi%nodaltilebox(1); pFlx => flx(1,lvl)%dataptr(mfi)
            pFlx(fbx%lo(1):fbx%hi(1), fbx%lo(2):fbx%hi(2), fbx%lo(3):fbx%hi(3),:) = FX(fbx%lo(1):fbx%hi(1), fbx%lo(2):fbx%hi(2), fbx%lo(3):fbx%hi(3),:)
            fbx = mfi%nodaltilebox(2); pFlx => flx(2,lvl)%dataptr(mfi)
            pFlx(fbx%lo(1):fbx%hi(1), fbx%lo(2):fbx%hi(2), fbx%lo(3):fbx%hi(3),:) = FY(fbx%lo(1):fbx%hi(1), fbx%lo(2):fbx%hi(2), fbx%lo(3):fbx%hi(3),:)
            fbx = mfi%nodaltilebox(3); pFlx => flx(3,lvl)%dataptr(mfi)
            pFlx(fbx%lo(1):fbx%hi(1), fbx%lo(2):fbx%hi(2), fbx%lo(3):fbx%hi(3),:) = FZ(fbx%lo(1):fbx%hi(1), fbx%lo(2):fbx%hi(2), fbx%lo(3):fbx%hi(3),:)
         end do
         call this%amr%mfiter_destroy(mfi)
         call this%amr%mfab_destroy(SCfill)
      end do

      ! Phase 2: Handle coarse-fine interface (flux averaging only)
      if (.not.this%use_refluxing) then
         ! Flux averaging: average fine fluxes down to coarse
         do lvl = this%amr%clvl(), 1, -1
            call amrmfab_average_down_faces(flx(1,lvl)%p, flx(2,lvl)%p, flx(3,lvl)%p, &
            &   flx(1,lvl-1)%p, flx(2,lvl-1)%p, flx(3,lvl-1)%p, &
            &   this%amr%geom(lvl-1)%p, 1, this%nscalar, &
            &   [this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)])
         end do
      end if

      ! Phase 3: Compute divergence to get dSCdt
      do lvl = 0, this%amr%clvl()
         call this%amr%mfiter_build(lvl, mfi)
         do while (mfi%next())
            bx = mfi%tilebox()
            rhs => dSCdt%mf(lvl)%dataptr(mfi)
            FX => flx(1,lvl)%dataptr(mfi)
            FY => flx(2,lvl)%dataptr(mfi)
            FZ => flx(3,lvl)%dataptr(mfi)
            do nsc = 1, this%nscalar
               do k = bx%lo(3), bx%hi(3)
                  do j = bx%lo(2), bx%hi(2)
                     do i = bx%lo(1), bx%hi(1)
                        rhs(i,j,k,nsc) = (FX(i+1,j,k,nsc)-FX(i,j,k,nsc))/this%amr%dx(lvl) + &
                        &                (FY(i,j+1,k,nsc)-FY(i,j,k,nsc))/this%amr%dy(lvl) + &
                        &                (FZ(i,j,k+1,nsc)-FZ(i,j,k,nsc))/this%amr%dz(lvl)
                     end do
                  end do
               end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do

      ! Phase 4: Store scaled fluxes in FluxRegister (refluxing only)
      if (this%use_refluxing) then
         do lvl = 0, this%amr%clvl()
            ! Scale fluxes by face areas
            call flx(1,lvl)%mult(this%amr%dy(lvl) * this%amr%dz(lvl), icomp=1, ncomp=this%nscalar, nghost=0)
            call flx(2,lvl)%mult(this%amr%dz(lvl) * this%amr%dx(lvl), icomp=1, ncomp=this%nscalar, nghost=0)
            call flx(3,lvl)%mult(this%amr%dx(lvl) * this%amr%dy(lvl), icomp=1, ncomp=this%nscalar, nghost=0)
            ! Add to flux register
            if (lvl .gt. 0) call this%flux%fineadd(lvl, [flx(1,lvl), flx(2,lvl), flx(3,lvl)], -1.0_WP)
            if (lvl .lt. this%amr%clvl()) call this%flux%crseinit(lvl+1, [flx(1,lvl), flx(2,lvl), flx(3,lvl)], +1.0_WP)
         end do
      end if

      ! Cleanup temp flux storage
      do lvl = 0, this%amr%clvl()
         call this%amr%mfab_destroy(flx(1,lvl))
         call this%amr%mfab_destroy(flx(2,lvl))
         call this%amr%mfab_destroy(flx(3,lvl))
      end do


   end subroutine get_dSCdt


   !> Apply flux register correction (only if use_refluxing=.true.)
   subroutine reflux(this, dt)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: lvl
      interface
         subroutine amrex_fi_fluxregister_reflux(fr, mf, scale, geom) bind(c)
            import :: c_ptr, WP
            type(c_ptr), value :: fr, mf, geom
            real(WP), value :: scale
         end subroutine amrex_fi_fluxregister_reflux
      end interface
      if (.not. this%use_refluxing) return
      do lvl = this%amr%clvl()-1, 0, -1
         call amrex_fi_fluxregister_reflux(this%flux%fr(lvl+1)%p, this%SC%mf(lvl)%p, dt, this%amr%geom(lvl)%p)
      end do
   end subroutine reflux


   !> Average down SC from fine to coarse (all levels)
   subroutine average_down(this)
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: lvl
      do lvl = this%amr%clvl()-1, 0, -1
         call this%SC%average_downto(lvl)
      end do
   end subroutine average_down


   !> Register solver data for checkpoint
   subroutine register_checkpoint(this, io)
      class(amrscalar), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call io%add_data(this%SC, 'SC')
   end subroutine register_checkpoint


   !> Restore solver data from checkpoint
   subroutine restore_checkpoint(this, io, dirname, time)
      class(amrscalar), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      call io%read_data(dirname, this%SC, 'SC')
      call this%SC%fill(time=time)
   end subroutine restore_checkpoint

end module amrscalar_class
