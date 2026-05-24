!> Test amrincomp solver - multi-level projection test
module mod_test_project
   use precision,         only: WP
   use random,            only: random_uniform
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrincomp_class,   only: amrincomp
   use amrdata_class,     only: amrdata, amrex_interp_none
   use messager,          only: log
   use string,            only: itoa, rtoa
   use amrex_amr_module,  only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, &
   &                            amrex_mfiter_build, amrex_mfiter_destroy
   implicit none
   private
   public :: test_project

   ! Grid
   type(amrgrid), allocatable, target :: amr

   ! Solver data
   type(amrincomp), allocatable, target :: fs

   ! Visualization
   type(amrviz), allocatable :: viz

contains

   !> User-provided initialization callback for velocity - random values
   subroutine velocity_init(solver, lvl, time, ba, dm)
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      integer :: i, j, k

      ! For fine levels: fill from coarse (no direct computation)
      !if (lvl > 0) then
      !   call solver%U%fill_from_coarse(lvl, time)
      !   call solver%V%fill_from_coarse(lvl, time)
      !   call solver%W%fill_from_coarse(lvl, time)
      !   return
      !end if

      ! Build mfiter from cell-centered ba/dm (passed during regrid)
      call amrex_mfiter_build(mfi, ba, dm, tiling=solver%amr%default_tiling)
      do while (mfi%next())
         pU => solver%U%mf(lvl)%dataptr(mfi)
         pV => solver%V%mf(lvl)%dataptr(mfi)
         pW => solver%W%mf(lvl)%dataptr(mfi)
         ! U (x-faces): random
         fbx = mfi%nodaltilebox(1)
         do k = fbx%lo(3), fbx%hi(3); do j = fbx%lo(2), fbx%hi(2); do i = fbx%lo(1), fbx%hi(1)
            pU(i,j,k,1) = random_uniform(lo=-1.0_WP, hi=1.0_WP)
         end do; end do; end do
         ! V (y-faces): random
         fbx = mfi%nodaltilebox(2)
         do k = fbx%lo(3), fbx%hi(3); do j = fbx%lo(2), fbx%hi(2); do i = fbx%lo(1), fbx%hi(1)
            pV(i,j,k,1) = random_uniform(lo=-1.0_WP, hi=1.0_WP)
         end do; end do; end do
         ! W (z-faces): random
         fbx = mfi%nodaltilebox(3)
         do k = fbx%lo(3), fbx%hi(3); do j = fbx%lo(2), fbx%hi(2); do i = fbx%lo(1), fbx%hi(1)
            pW(i,j,k,1) = random_uniform(lo=-1.0_WP, hi=1.0_WP)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine velocity_init

   !> Geometric tagger: refine right half of domain (x > 0.5)
   !> This creates a C/F interface at the x=0/x=1 periodic boundary
   subroutine geometric_tagger(solver, lvl, time, tags_ptr)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP) :: x, dx
      integer :: i, j, k

      dx = solver%amr%dx(lvl)

      ! Convert c_ptr to tagboxarray
      tags = tags_ptr

      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
            x = solver%amr%xlo + (real(i,WP) + 0.5_WP) * dx
            ! Refine strip 0.5 < x < 0.94: buffer should extend right edge to x=1 (periodic boundary)
            if (x .gt. 0.5_WP .and. x .lt. 0.94_WP) tagarr(i,j,k,1) = SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine geometric_tagger

   !> Main test routine
   subroutine test_project()
      implicit none

      call log("========================================")
      call log("=TEST: amrincomp multi-level projection=")
      call log("========================================")

      ! Create amrgrid
      create_amrgrid: block
         use param, only: param_read
         allocate(amr)
         amr%name = 'project_test'
         call param_read('Base nx', amr%nx)
         call param_read('Base ny', amr%ny)
         call param_read('Base nz', amr%nz)
         amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
         amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
         amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
         amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
         call param_read('Max level', amr%maxlvl)
         call amr%initialize()
         call log("Grid initialized: "//trim(itoa(amr%nx))//"x"//trim(itoa(amr%ny))//"x"//trim(itoa(amr%nz)))
         call log("Max levels: "//trim(itoa(amr%maxlvl+1)))
      end block create_amrgrid

      ! Create flow solver
      create_flow_solver: block
         allocate(fs)
         fs%rho = 1.0_WP
         fs%user_init => velocity_init
         fs%user_tagging => geometric_tagger
         call fs%initialize(amr, name='project_fs')
         call log("Flow solver initialized")
      end block create_flow_solver

      ! Initialize grid (triggers callbacks)
      initialize: block
         call amr%init_from_scratch(time=0.0_WP)
         call log("Grid built: "//trim(itoa(amr%nlevels))//" levels")
         ! Average down velocity first to make coarse consistent with fine
         call fs%average_down_velocity()
         ! Fill velocity ghosts
         call fs%fill_velocity(time=0.0_WP)
         call log("Velocity averaged down and ghosts filled")
         ! Setup pressure solver
         fs%psolver%verbose = 2
         call fs%psolver%setup()
         call log("Pressure solver setup complete")
      end block initialize

      ! Create visualization
      create_visualization: block
         ! Prepare visualization
         allocate(viz); call viz%initialize(amr, 'test_project')
         call viz%add_scalar(fs%P, 1, 'pressure')
         call viz%add_scalar(fs%div, 1, 'divergence')
         call viz%add_scalar(fs%U, 1, 'U')
         call viz%add_scalar(fs%V, 1, 'V')
         call viz%add_scalar(fs%W, 1, 'W')
         ! Write initial visualization output
         call viz%write(time=0.0_WP)
      end block create_visualization

      ! Perform projection
      project: block

         ! 1. Compute divergence
         call fs%get_div()
         call log("Before projection: divmax = "//trim(rtoa(fs%divmax)))
         call log("RHS sum = "//trim(rtoa(fs%div%get_sum(lvl=0))))
         call viz%write(time=0.0_WP)

         ! 2. Solve pressure Poisson equation
         call fs%psolver%solve(rhs=fs%div)
         call log("Pressure solver residual = "//trim(rtoa(fs%psolver%res)))

         ! 3. Correct velocity
         call fs%add_pressure(scale=1.0_WP)

         ! 4. Average down for C/F consistency
         call fs%average_down_velocity()

         ! 5. Fill ghosts
         call fs%fill_velocity(time=0.0_WP)

         ! 6. Check divergence
         call fs%get_div()
         call log("After projection: divmax = "//trim(rtoa(fs%divmax)))

         ! Report result
         if (fs%divmax .lt. 1.0e-6_WP) then
            call log("TEST PASSED: Velocity is divergence-free (divmax < 1e-6)")
         else
            call log("TEST FAILED: Divergence not sufficiently reduced")
         end if

         ! Write visualization output
         call viz%write(time=1.0_WP)
         call log("Visualization written to amrviz/test_project/")
      end block project

      ! Cleanup
      cleanup: block
         call viz%finalize()
         call fs%finalize()
         call amr%finalize()
         deallocate(viz, fs, amr)
      end block cleanup

      call log("Test complete")

   end subroutine test_project

end module mod_test_project
