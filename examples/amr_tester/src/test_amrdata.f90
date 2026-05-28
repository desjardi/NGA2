!> Stress test for amrdata callback architecture
!> Tests:
!>   1. data1: Default callbacks (baseline)
!>   2. data2: Custom on_init with context verification
!>   3. data3: Custom on_init with parent pointer pattern
!>   4. data4: Different ncomp/ng values
!>   5. data5: Gaussian bump scalar (for value-based tagging)
!>   6. data6: Face-centered data (staggered X-velocity type)
module mod_test_amrdata
   use precision,        only: WP
   use string,           only: str_medium,rtoa,itoa
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata
   use amrviz_class,     only: amrviz
   use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab,amrex_geometry
   use messager,         only: log,warn,die
   implicit none
   private
   public :: test_amrdata

   !> Mock "solver" type to test parent pointer pattern
   type :: mock_solver
      real(WP) :: fill_value = 0.0_WP
      character(len=32) :: solver_name = 'UNNAMED'
   end type mock_solver

   !> Gaussian parameters for data5
   type :: gaussian_params
      real(WP) :: x0 = 0.5_WP  !< Center x
      real(WP) :: y0 = 0.5_WP  !< Center y
      real(WP) :: z0 = 0.5_WP  !< Center z
      real(WP) :: sigma = 0.15_WP  !< Width parameter
      real(WP) :: amplitude = 1.0_WP  !< Peak amplitude
   end type gaussian_params

contains

   !> Tagger: tag bottom-left corner for refinement
   subroutine corner_tagger(ctx, lvl, time, tags_ptr)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_char
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box
      use amrgrid_class, only: SETtag
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrgrid), pointer :: amr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP) :: x,y,z,dx,dy,dz
      integer :: i,j,k
      call c_f_pointer(ctx, amr)
      tags = tags_ptr
      dx = amr%geom(lvl)%dx(1)
      dy = amr%geom(lvl)%dx(2)
      dz = amr%geom(lvl)%dx(3)
      call amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = amr%zlo + (real(k,WP)+0.5_WP)*dz
            do j = bx%lo(2), bx%hi(2)
               y = amr%ylo + (real(j,WP)+0.5_WP)*dy
               do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  if (x < 0.2_WP .and. y < 0.2_WP .and. z < 0.2_WP) then
                     tagarr(i,j,k,1) = SETtag
                  end if
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)
   end subroutine corner_tagger

   !> Tagger: tag cells where scalar value exceeds threshold
   !> Context is amrdata pointer - demonstrates value-based refinement
   subroutine value_tagger(ctx, lvl, time, tags_ptr)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_associated
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_mfiter_build, &
      &                           amrex_mfiter_destroy, amrex_box
      use amrgrid_class, only: SETtag
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrdata), pointer :: data
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), contiguous, pointer :: scalar(:,:,:,:)
      real(WP), parameter :: threshold = 0.3_WP
      integer :: i,j,k
      ! Recover amrdata from context
      call c_f_pointer(ctx, data)
      ! Skip if MultiFab not yet allocated
      if (.not.c_associated(data%mf(lvl)%p)) return
      tags = tags_ptr
      ! Iterate over the scalar field and tag based on value
      call data%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         scalar => data%mf(lvl)%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  if (scalar(i,j,k,1) > threshold) tagarr(i,j,k,1) = SETtag
               end do
            end do
         end do
      end do
      call data%amr%mfiter_destroy(mfi)
   end subroutine value_tagger

   !> Custom on_init for data2: fills MultiFab with a constant value (1.0)
   !> Demonstrates that callbacks receive typed class(amrdata)
   subroutine custom_on_init_fill_one(this, lvl, time, ba, dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Build the MultiFab
      call this%reset_level(lvl, ba, dm)
      ! Fill with 1.0 - demonstrates we can access 'this' properly
      call this%setval(val=1.0_WP, lvl=lvl)
      call log("    custom_on_init_fill_one: filled level with 1.0")
   end subroutine custom_on_init_fill_one

   !> Custom on_init for data3: uses parent pointer to get fill value
   !> Demonstrates the parent pointer pattern for solver context access
   subroutine custom_on_init_from_parent(this, lvl, time, ba, dm)
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      real(WP) :: fill_val
      ! Build the MultiFab
      call this%reset_level(lvl, ba, dm)
      ! Retrieve fill value from parent context
      fill_val = 0.0_WP
      if (associated(this%parent)) then
         select type (p => this%parent)
          type is (mock_solver)
            fill_val = p%fill_value
            call log("    custom_on_init_from_parent: recovered parent fill_value = " // trim(rtoa(fill_val)))
          class default
            call warn("    custom_on_init_from_parent: unknown parent type!")
         end select
      else
         call warn("    custom_on_init_from_parent: parent not associated!")
      end if
      call this%setval(val=fill_val, lvl=lvl)
   end subroutine custom_on_init_from_parent

   !> Custom user_init for data5: fills MultiFab with a 3D Gaussian bump
   !> Uses parent pointer to get Gaussian parameters
   subroutine custom_user_init_gaussian(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_mfiter, amrex_mfiter_build, amrex_mfiter_destroy, amrex_box
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), contiguous, pointer :: arr(:,:,:,:)
      real(WP) :: x,y,z,dx,dy,dz,r2,x0,y0,z0,sigma,amp
      integer :: i,j,k,lo(4),hi(4)
      ! Get Gaussian parameters from parent
      if (.not.associated(this%parent)) call die("custom_user_init_gaussian: parent pointer not set!")
      select type (p => this%parent)
       type is (gaussian_params)
         x0 = p%x0; y0 = p%y0; z0 = p%z0
         sigma = p%sigma; amp = p%amplitude
       class default
         call die("custom_user_init_gaussian: parent is not gaussian_params type!")
      end select
      ! Fill with Gaussian bump: amp * exp(-r^2 / sigma^2)
      dx = this%amr%geom(lvl)%dx(1)
      dy = this%amr%geom(lvl)%dx(2)
      dz = this%amr%geom(lvl)%dx(3)
      call amrex_mfiter_build(mfi, this%mf(lvl), tiling=.true.)
      do while (mfi%next())
         bx = mfi%tilebox()
         arr => this%mf(lvl)%dataPtr(mfi)
         lo = lbound(arr); hi = ubound(arr)
         do k = bx%lo(3), bx%hi(3)
            z = this%amr%zlo + (real(k,WP)+0.5_WP)*dz
            do j = bx%lo(2), bx%hi(2)
               y = this%amr%ylo + (real(j,WP)+0.5_WP)*dy
               do i = bx%lo(1), bx%hi(1)
                  x = this%amr%xlo + (real(i,WP)+0.5_WP)*dx
                  r2 = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
                  arr(i,j,k,1) = amp * exp(-r2 / sigma**2)
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
      call log("    custom_user_init_gaussian: level " // trim(itoa(lvl)) // &
      &        " min=" // trim(rtoa(this%get_min(lvl=lvl))) // &
      &        " max=" // trim(rtoa(this%get_max(lvl=lvl))))
   end subroutine custom_user_init_gaussian

   !> Main test routine - multiple amrdata with varying configurations
   subroutine test_amrdata()
      use iso_c_binding, only: c_loc, c_associated
      implicit none
      type(amrgrid), target :: amr
      type(amrdata), target :: data1, data2, data3, data4, data5, data6
      type(mock_solver), target :: solver
      type(gaussian_params), target :: gauss
      type(amrviz) :: viz
      integer :: npassed, nfailed

      npassed = 0
      nfailed = 0

      call log("===================================================")
      call log("TEST: Multi-instance amrdata stress test")
      call log("===================================================")

      ! Setup amrgrid
      amr%nx = 16; amr%ny = 16; amr%nz = 16
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper=.false.
      amr%yper=.false.
      amr%zper=.false.
      amr%maxlvl = 2
      amr%nmax = 16
      call amr%initialize("TestGrid")

      ! Setup mock solver for parent pointer test
      solver%fill_value = 42.0_WP
      solver%solver_name = 'MockSolver'

      ! -------------------------------------------------------------------
      ! TEST 1: data1 - Default callbacks (baseline)
      ! -------------------------------------------------------------------
      call log("")
      call log("--- TEST 1: data1 with default callbacks ---")
      call data1%initialize(amr, 'data1', ncomp=1, ng=0)
      call data1%register()
      call log("    data1 initialized and registered")

      ! -------------------------------------------------------------------
      ! TEST 2: data2 - Custom on_init that fills with 1.0
      ! -------------------------------------------------------------------
      call log("")
      call log("--- TEST 2: data2 with custom on_init (fill 1.0) ---")
      call data2%initialize(amr, 'data2', ncomp=1, ng=0)
      data2%on_init => custom_on_init_fill_one
      call data2%register()
      call log("    data2 initialized with custom callback and registered")

      ! -------------------------------------------------------------------
      ! TEST 3: data3 - Custom on_init with parent pointer
      ! -------------------------------------------------------------------
      call log("")
      call log("--- TEST 3: data3 with parent pointer pattern ---")
      call data3%initialize(amr, 'data3', ncomp=1, ng=0)
      data3%parent => solver  ! Set parent pointer
      data3%on_init => custom_on_init_from_parent
      call data3%register()
      call log("    data3 initialized with parent pointer and registered")

      ! -------------------------------------------------------------------
      ! TEST 4: data4 - Different ncomp and ng values
      ! -------------------------------------------------------------------
      call log("")
      call log("--- TEST 4: data4 with ncomp=3, ng=1 ---")
      call data4%initialize(amr, 'data4', ncomp=3, ng=1)
      call data4%register()
      call log("    data4 initialized with ncomp=3, ng=1 and registered")

      ! -------------------------------------------------------------------
      ! TEST 5: data5 - Gaussian bump scalar
      ! -------------------------------------------------------------------
      call log("")
      call log("--- TEST 5: data5 with Gaussian bump ---")
      gauss%x0 = 0.8_WP; gauss%y0 = 0.8_WP; gauss%z0 = 0.8_WP  ! Top-right corner
      gauss%sigma = 0.08_WP; gauss%amplitude = 1.0_WP  ! Narrow Gaussian for small region
      call data5%initialize(amr, 'data5', ncomp=1, ng=0)
      data5%parent => gauss
      data5%user_init => custom_user_init_gaussian
      call data5%register()
      call log("    data5 initialized with Gaussian at (0.8,0.8,0.8) and registered")

      ! -------------------------------------------------------------------
      ! TEST 6: data6 - Face-centered data (staggered X-velocity type)
      ! -------------------------------------------------------------------
      call log("")
      call log("--- TEST 6: data6 face-centered (nodal X) ---")
      call data6%initialize(amr, 'data6', ncomp=1, ng=1, nodal=[.true.,.false.,.false.])
      call data6%register()
      call log("    data6 initialized with nodal=[T,F,F] and registered")

      ! -------------------------------------------------------------------
      ! Initialize grid - triggers all on_init callbacks
      ! -------------------------------------------------------------------
      call log("")
      call log("--- Initializing grid (triggers all callbacks) ---")
      call log("    Adding corner_tagger (bottom-left) with amr context")
      call amr%add_tagging(corner_tagger, c_loc(amr))
      call log("    Adding value_tagger (top-right Gaussian) with data5 context")
      call amr%add_tagging(value_tagger, c_loc(data5))
      call amr%init_from_scratch(time=0.0_WP, do_postregrid=.true.)
      call log("Grid initialized with " // trim(itoa(amr%nlevels)) // " levels")

      ! -------------------------------------------------------------------
      ! VERIFICATION
      ! -------------------------------------------------------------------
      call log("")
      call log("--- Verification ---")

      ! Check data1 allocated
      if (c_associated(data1%mf(0)%p)) then
         call log("PASS: data1%mf(0) allocated")
         npassed = npassed + 1
      else
         call warn("FAIL: data1%mf(0) NOT allocated")
         nfailed = nfailed + 1
      end if

      ! Check data2 allocated and has correct value
      if (c_associated(data2%mf(0)%p)) then
         if (abs(data2%get_min(lvl=0) - 1.0_WP) < 1.0e-10_WP) then
            call log("PASS: data2%mf(0) filled with 1.0")
            npassed = npassed + 1
         else
            call warn("FAIL: data2%mf(0) has wrong value: " // trim(rtoa(data2%get_min(lvl=0))))
            nfailed = nfailed + 1
         end if
      else
         call warn("FAIL: data2%mf(0) NOT allocated")
         nfailed = nfailed + 1
      end if

      ! Check data3 allocated and has value from parent (42.0)
      if (c_associated(data3%mf(0)%p)) then
         if (abs(data3%get_min(lvl=0) - 42.0_WP) < 1.0e-10_WP) then
            call log("PASS: data3%mf(0) filled with parent value 42.0")
            npassed = npassed + 1
         else
            call warn("FAIL: data3%mf(0) has wrong value: " // trim(rtoa(data3%get_min(lvl=0))) // " (expected 42.0)")
            nfailed = nfailed + 1
         end if
      else
         call warn("FAIL: data3%mf(0) NOT allocated")
         nfailed = nfailed + 1
      end if

      ! Check data4 has correct ncomp
      if (c_associated(data4%mf(0)%p)) then
         if (data4%mf(0)%ncomp() == 3) then
            call log("PASS: data4%mf(0) has ncomp=3")
            npassed = npassed + 1
         else
            call warn("FAIL: data4%mf(0) wrong ncomp: " // trim(itoa(data4%mf(0)%ncomp())))
            nfailed = nfailed + 1
         end if
      else
         call warn("FAIL: data4%mf(0) NOT allocated")
         nfailed = nfailed + 1
      end if

      ! Check level 1 exists for all (if maxlvl >= 1)
      if (amr%nlevels > 1) then
         if (c_associated(data1%mf(1)%p) .and. c_associated(data2%mf(1)%p) .and. &
         &   c_associated(data3%mf(1)%p) .and. c_associated(data4%mf(1)%p) .and. &
         &   c_associated(data5%mf(1)%p) .and. c_associated(data6%mf(1)%p)) then
            call log("PASS: All data objects have level 1 allocated")
            npassed = npassed + 1
         else
            call warn("FAIL: Some data objects missing level 1")
            nfailed = nfailed + 1
         end if
      end if

      ! Check data5 Gaussian bump: max should be ~0.8+ at center (coarse grid), min should be small
      if (c_associated(data5%mf(0)%p)) then
         if (data5%get_max(lvl=0) > 0.8_WP .and. data5%get_min(lvl=0) < 0.1_WP) then
            call log("PASS: data5%mf(0) Gaussian bump (max=" // trim(rtoa(data5%get_max(lvl=0))) // &
            &        ", min=" // trim(rtoa(data5%get_min(lvl=0))) // ")")
            npassed = npassed + 1
         else
            call warn("FAIL: data5%mf(0) unexpected Gaussian values")
            nfailed = nfailed + 1
         end if
      else
         call warn("FAIL: data5%mf(0) NOT allocated")
         nfailed = nfailed + 1
      end if

      ! Check data6 face-centered: verify nodal flag and allocation
      if (c_associated(data6%mf(0)%p)) then
         if (data6%nodal(1) .and. .not.data6%nodal(2) .and. .not.data6%nodal(3)) then
            call log("PASS: data6%mf(0) face-centered (nodal=[T,F,F])")
            npassed = npassed + 1
         else
            call warn("FAIL: data6 has wrong nodal flags")
            nfailed = nfailed + 1
         end if
      else
         call warn("FAIL: data6%mf(0) NOT allocated")
         nfailed = nfailed + 1
      end if

      ! -------------------------------------------------------------------
      ! Summary
      ! -------------------------------------------------------------------
      call log("")
      call log("===================================================")
      call log("RESULTS: " // trim(itoa(npassed)) // " passed, " // trim(itoa(nfailed)) // " failed")
      call log("===================================================")

      ! -------------------------------------------------------------------
      ! Write to ParaView using amrviz
      ! -------------------------------------------------------------------
      call log("")
      call log("--- Writing to ParaView via amrviz ---")
      call viz%initialize(amr, 'test_amrdata')
      call viz%add_scalar(data1, 1, 'data1')
      call viz%add_scalar(data2, 1, 'data2')
      call viz%add_scalar(data3, 1, 'data3')
      call viz%add_scalar(data5, 1, 'data5')
      call viz%write(0.0_WP)
      call viz%finalize()
      call log("    Output written to amrviz/test_amrdata/")

      ! Cleanup
      call data1%finalize()
      call data2%finalize()
      call data3%finalize()
      call data4%finalize()
      call data5%finalize()
      call data6%finalize()
      call amr%finalize()

      if (nfailed > 0) then
         call warn("=== TEST FAILED ===")
      else
         call log("=== TEST COMPLETE ===")
      end if

   end subroutine test_amrdata

end module mod_test_amrdata
