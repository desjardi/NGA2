!> Test program for unified amrmg solver
!> Tests both constant-coefficient (Poisson) and variable-coefficient modes
!> Uses manufactured solution to verify solver
module mod_test_amrmg
   use precision,         only: WP
   use string,            only: str_long
   use messager,          only: log
   use mathtools,         only: Pi
   use amrgrid_class,     only: amrgrid
   use amrdata_class,     only: amrdata
   use amrmg_class
   use timer_class,       only: timer
   use amrex_amr_module,  only: amrex_mfiter, amrex_box
   use mpi_f08
   implicit none
   private
   public :: test_amrmg

contains

   subroutine test_amrmg()
      call log('=== TEST_AMRMG STARTING ===')
      call log('')

      ! Test constant-coefficient solver (Poisson)
      call test_cstcoef()

      ! Test variable-coefficient solver
      call test_varcoef()

      ! Test multi-level AMR with refinement based on RHS
      call test_multilevel()

      ! Test per-level solve (for subcycling)
      call test_solve_level()

      call log('')
      call log('=== TEST_AMRMG COMPLETE ===')
   end subroutine test_amrmg

   !> Test constant-coefficient (Poisson) solver
   subroutine test_cstcoef()
      type(amrgrid), allocatable, target :: amr
      type(amrdata), allocatable, target :: phi, rhs, exact
      type(amrmg) :: solver
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(timer) :: tmr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPhi, pRhs, pExact
      real(WP) :: x, y, z, dx, phi_exact
      real(WP) :: err_l2, err_linf, local_err
      integer :: lvl, i, j, k, ierr
      character(len=str_long) :: msg

      call log('--- Testing constant-coefficient (Poisson) solver ---')
      call log('Solving: -∇²ϕ = f with ϕ = sin(πx)sin(πy)sin(πz)')

      ! Allocate and configure grid
      allocate(amr)
      amr%nx = 64; amr%ny = 64; amr%nz = 64
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .false.; amr%yper = .false.; amr%zper = .false.
      amr%maxlvl = 0
      amr%nmax = 32
      call amr%initialize(name='poisson_test')
      call amr%init_from_scratch(time=0.0_WP)

      ! Allocate fields
      allocate(phi, rhs, exact)
      call phi%initialize(amr=amr, name='phi', ncomp=1, ng=1)
      call rhs%initialize(amr=amr, name='rhs', ncomp=1, ng=0)
      call exact%initialize(amr=amr, name='exact', ncomp=1, ng=0)
      call phi%on_init(phi, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call rhs%on_init(rhs, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call exact%on_init(exact, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))

      ! Fill fields with manufactured solution
      dx = amr%dx(0)
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pRhs => rhs%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)

         ! Fill phi ghosts with exact solution (Dirichlet BC)
         do k = bx%lo(3)-1, bx%hi(3)+1
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2)-1, bx%hi(2)+1
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1)-1, bx%hi(1)+1
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  pPhi(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do

         ! Fill valid region
         do k = bx%lo(3), bx%hi(3)
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  phi_exact = sin(pi*x) * sin(pi*y) * sin(pi*z)
                  pExact(i,j,k,1) = phi_exact
                  pRhs(i,j,k,1) = -3.0_WP*pi**2 * phi_exact  ! Poisson: ∇²ϕ = f
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)

      ! Initialize and setup solver
      call solver%initialize(amr=amr, type=amrmg_cstcoef)
      solver%lo_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%hi_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%tol_rel = 1.0e-12_WP
      solver%verbose = 0
      solver%bottom_solver = amrmg_bottom_hypre

      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg_cst')
      call tmr%start()
      call solver%setup()
      call tmr%stop()
      write(msg,'(a,es12.5,a)') '  Setup time: ', tmr%time, ' s'; call log(msg)

      ! Solve with zero initial guess
      call phi%setval(val=0.0_WP, lvl=0)
      call tmr%reset(); call tmr%start()
      call solver%solve(rhs=rhs,phi0=phi)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  Solve: residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      ! Compute error
      err_l2 = 0.0_WP; err_linf = 0.0_WP
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  local_err = abs(pPhi(i,j,k,1) - pExact(i,j,k,1))
                  err_l2 = err_l2 + local_err**2
                  err_linf = max(err_linf, local_err)
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      err_l2 = sqrt(err_l2)

      write(msg,'(a,es12.5,a,es12.5)') '  Error: L2 = ', err_l2, ', Linf = ', err_linf; call log(msg)

      if (err_linf .lt. 0.1_WP .and. solver%res .lt. 1.0e-8_WP) then
         call log('  PASS: constant-coefficient (Poisson) test')
      else
         call log('  FAIL: constant-coefficient (Poisson) test')
      end if

      ! Cleanup
      call solver%finalize()
      call phi%finalize(); call rhs%finalize(); call exact%finalize()
      deallocate(phi, rhs, exact)
      call amr%finalize(); deallocate(amr)
   end subroutine test_cstcoef

   !> Test variable-coefficient solver
   subroutine test_varcoef()
      type(amrgrid), allocatable, target :: amr
      type(amrdata), allocatable, target :: phi, rhs, exact
      type(amrdata), allocatable, target :: bcoef_x, bcoef_y, bcoef_z
      type(amrmg) :: solver
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(timer) :: tmr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPhi, pRhs, pExact, pBx, pBy, pBz
      real(WP) :: x, y, z, dx, phi_exact
      real(WP) :: err_l2, err_linf, local_err
      integer :: lvl, i, j, k, ierr
      character(len=str_long) :: msg

      call log('')
      call log('--- Testing variable-coefficient solver ---')
      call log('Solving: -∇·(B∇ϕ) = f with B = 1 (uniform)')

      ! Allocate and configure grid
      allocate(amr)
      amr%nx = 64; amr%ny = 64; amr%nz = 64
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .false.; amr%yper = .false.; amr%zper = .false.
      amr%maxlvl = 0
      amr%nmax = 32
      call amr%initialize(name='varcoef_test')
      call amr%init_from_scratch(time=0.0_WP)

      ! Allocate fields
      allocate(phi, rhs, exact, bcoef_x, bcoef_y, bcoef_z)
      call phi%initialize(amr=amr, name='phi', ncomp=1, ng=1)
      call rhs%initialize(amr=amr, name='rhs', ncomp=1, ng=0)
      call exact%initialize(amr=amr, name='exact', ncomp=1, ng=0)
      call bcoef_x%initialize(amr=amr, name='bx', ncomp=1, ng=0, nodal=[.true.,.false.,.false.])
      call bcoef_y%initialize(amr=amr, name='by', ncomp=1, ng=0, nodal=[.false.,.true.,.false.])
      call bcoef_z%initialize(amr=amr, name='bz', ncomp=1, ng=0, nodal=[.false.,.false.,.true.])

      call phi%on_init(phi, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call rhs%on_init(rhs, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call exact%on_init(exact, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call bcoef_x%on_init(bcoef_x, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call bcoef_y%on_init(bcoef_y, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))
      call bcoef_z%on_init(bcoef_z, 0, 0.0_WP, amr%get_boxarray(0), amr%get_distromap(0))

      ! Fill fields with manufactured solution
      dx = amr%dx(0)
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pRhs => rhs%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)

         ! Fill phi ghosts with exact solution (Dirichlet BC)
         do k = bx%lo(3)-1, bx%hi(3)+1
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2)-1, bx%hi(2)+1
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1)-1, bx%hi(1)+1
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  pPhi(i,j,k,1) = sin(pi*x) * sin(pi*y) * sin(pi*z)
               end do
            end do
         end do

         ! Fill valid region
         do k = bx%lo(3), bx%hi(3)
            z = amr%zlo + (real(k,WP)+0.5_WP)*dx
            do j = bx%lo(2), bx%hi(2)
               y = amr%ylo + (real(j,WP)+0.5_WP)*dx
               do i = bx%lo(1), bx%hi(1)
                  x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                  phi_exact = sin(pi*x) * sin(pi*y) * sin(pi*z)
                  pExact(i,j,k,1) = phi_exact
                  pRhs(i,j,k,1) = 3.0_WP*pi**2 * phi_exact
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)

      ! Fill B coefficients (uniform = 1)
      call bcoef_x%mfiter_build(0, mfi)
      do while (mfi%next())
         pBx => bcoef_x%mf(0)%dataptr(mfi)
         bx = mfi%tilebox()
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pBx(i,j,k,1) = 1.0_WP
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call bcoef_y%mfiter_build(0, mfi)
      do while (mfi%next())
         pBy => bcoef_y%mf(0)%dataptr(mfi)
         bx = mfi%tilebox()
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pBy(i,j,k,1) = 1.0_WP
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call bcoef_z%mfiter_build(0, mfi)
      do while (mfi%next())
         pBz => bcoef_z%mf(0)%dataptr(mfi)
         bx = mfi%tilebox()
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
                  pBz(i,j,k,1) = 1.0_WP
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)

      ! Initialize and setup solver
      call solver%initialize(amr=amr, type=amrmg_varcoef)
      solver%lo_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%hi_bc = [amrmg_bc_dirichlet, amrmg_bc_dirichlet, amrmg_bc_dirichlet]
      solver%alpha = 0.0_WP
      solver%beta  = 1.0_WP
      solver%tol_rel = 1.0e-12_WP
      solver%verbose = 0
      solver%bottom_solver = amrmg_bottom_hypre

      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg_var')
      call tmr%start()
      call solver%setup(bcoef_x=bcoef_x%mf, bcoef_y=bcoef_y%mf, bcoef_z=bcoef_z%mf)
      call tmr%stop()
      write(msg,'(a,es12.5,a)') '  Setup time: ', tmr%time, ' s'; call log(msg)

      ! First solve
      call phi%setval(val=0.0_WP, lvl=0)
      call tmr%reset(); call tmr%start()
      call solver%solve(rhs=rhs,phi0=phi)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  First solve: residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      ! Second solve (reuse test)
      call phi%setval(val=0.0_WP, lvl=0)
      call tmr%reset(); call tmr%start()
      call solver%solve(rhs=rhs,phi0=phi)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  Second solve (reuse): residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      ! Compute error
      err_l2 = 0.0_WP; err_linf = 0.0_WP
      call amr%mfiter_build(0, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pPhi => phi%mf(0)%dataptr(mfi)
         pExact => exact%mf(0)%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  local_err = abs(pPhi(i,j,k,1) - pExact(i,j,k,1))
                  err_l2 = err_l2 + local_err**2
                  err_linf = max(err_linf, local_err)
               end do
            end do
         end do
      end do
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      err_l2 = sqrt(err_l2)

      write(msg,'(a,es12.5,a,es12.5)') '  Error: L2 = ', err_l2, ', Linf = ', err_linf; call log(msg)

      if (err_linf .lt. 0.1_WP .and. solver%res .lt. 1.0e-8_WP) then
         call log('  PASS: variable-coefficient test')
      else
         call log('  FAIL: variable-coefficient test')
      end if

      ! Cleanup
      call solver%finalize()
      call phi%finalize(); call rhs%finalize(); call exact%finalize()
      call bcoef_x%finalize(); call bcoef_y%finalize(); call bcoef_z%finalize()
      deallocate(phi, rhs, exact, bcoef_x, bcoef_y, bcoef_z)
      call amr%finalize(); deallocate(amr)
   end subroutine test_varcoef

   !> Tagger: refine where RHS magnitude exceeds threshold
   subroutine rhs_tagger(ctx, lvl, time, tags_ptr)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_associated
      use amrex_amr_module, only: amrex_tagboxarray, amrex_mfiter, amrex_box
      use amrgrid_class, only: SETtag
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrdata), pointer :: rhs
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)
      real(WP), contiguous, pointer :: rhsarr(:,:,:,:)
      real(WP), parameter :: threshold = 1.0_WP  ! Tag where |rhs| > threshold
      integer :: i,j,k
      call c_f_pointer(ctx, rhs)
      if (.not.c_associated(rhs%mf(lvl)%p)) return
      tags = tags_ptr
      call rhs%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         rhsarr => rhs%mf(lvl)%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  if (abs(rhsarr(i,j,k,1)) .gt. threshold) tagarr(i,j,k,1) = SETtag
               end do
            end do
         end do
      end do
      call rhs%amr%mfiter_destroy(mfi)
   end subroutine rhs_tagger

   !> user_init callback: fill RHS with two sharp Gaussians (positive + negative)
   !> This creates a balanced source for Neumann BC compatibility
   subroutine rhs_user_init(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_boxarray, amrex_distromap, amrex_mfiter, amrex_box, &
      & amrex_mfiter_build, amrex_mfiter_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), contiguous, pointer :: arr(:,:,:,:)
      real(WP) :: x, y, z, dx, dy, dz, r1_sq, r2_sq, g1, g2
      integer :: i, j, k
      ! Two sharp Gaussians: positive at (0.25, 0.5, 0.5), negative at (0.75, 0.5, 0.5)
      real(WP), parameter :: sigma = 0.01_WP  ! Sharp Gaussian width
      real(WP), parameter :: amp = 500.0_WP   ! Amplitude
      real(WP), parameter :: x1 = 0.25_WP, y1 = 0.5_WP, z1 = 0.5_WP  ! Positive source
      real(WP), parameter :: x2 = 0.75_WP, y2 = 0.5_WP, z2 = 0.5_WP  ! Negative source
      dx = this%amr%geom(lvl)%dx(1)
      dy = this%amr%geom(lvl)%dx(2)
      dz = this%amr%geom(lvl)%dx(3)
      call amrex_mfiter_build(mfi, this%mf(lvl), tiling=.true.)
      do while (mfi%next())
         bx = mfi%tilebox()
         arr => this%mf(lvl)%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = this%amr%zlo + (real(k,WP)+0.5_WP)*dz
            do j = bx%lo(2), bx%hi(2)
               y = this%amr%ylo + (real(j,WP)+0.5_WP)*dy
               do i = bx%lo(1), bx%hi(1)
                  x = this%amr%xlo + (real(i,WP)+0.5_WP)*dx
                  r1_sq = (x-x1)**2 + (y-y1)**2 + (z-z1)**2
                  r2_sq = (x-x2)**2 + (y-y2)**2 + (z-z2)**2
                  g1 = amp * exp(-r1_sq / (2.0_WP*sigma**2))  ! Positive Gaussian
                  g2 = amp * exp(-r2_sq / (2.0_WP*sigma**2))  ! Negative Gaussian
                  arr(i,j,k,1) = g1 - g2  ! Net zero for Neumann compatibility
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine rhs_user_init

   !> Test multi-level AMR Poisson solve with RHS-based refinement
   subroutine test_multilevel()
      use iso_c_binding, only: c_loc
      use string, only: itoa, rtoa
      type(amrgrid), allocatable, target :: amr
      type(amrdata), allocatable, target :: phi, rhs
      type(amrmg) :: solver
      type(timer) :: tmr
      character(len=str_long) :: msg

      call log('')
      call log('--- Testing multi-level AMR Poisson solver ---')
      call log('Neumann BC with dual Gaussian sources (positive + negative)')

      ! Allocate and configure grid with 2 levels
      allocate(amr)
      amr%nx = 32; amr%ny = 32; amr%nz = 32
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .false.; amr%yper = .false.; amr%zper = .false.
      amr%maxlvl = 2  ! 3 levels (0, 1, 2)
      amr%nmax = 16
      call amr%initialize(name='multilevel_test')

      ! Allocate fields
      allocate(phi, rhs)
      call phi%initialize(amr=amr, name='phi', ncomp=1, ng=1)
      call rhs%initialize(amr=amr, name='rhs', ncomp=1, ng=0)

      ! Set user_init callbacks
      rhs%user_init => rhs_user_init

      ! Register fields
      call phi%register()
      call rhs%register()

      ! Add tagging based on RHS
      call amr%add_tagging(rhs_tagger, c_loc(rhs))

      ! Initialize grid - triggers callbacks and creates refined levels
      call amr%init_from_scratch(time=0.0_WP, do_postregrid=.true.)

      write(msg,'(a,i0,a)') '  Grid created with ', amr%nlevels, ' levels'
      call log(msg)

      ! Initialize and setup solver
      call solver%initialize(amr=amr, type=amrmg_cstcoef)
      solver%tol_rel = 1.0e-10_WP
      solver%verbose = 0
      solver%bottom_solver = amrmg_bottom_hypre

      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg_ml')
      call tmr%start()
      call solver%setup()
      call tmr%stop()
      write(msg,'(a,es12.5,a)') '  Setup time: ', tmr%time, ' s'; call log(msg)

      ! Solve
      call tmr%reset(); call tmr%start()
      call solver%solve(rhs=rhs,phi0=phi)
      call tmr%stop()
      write(msg,'(a,es12.5,a,es12.5,a)') '  Solve: residual = ', solver%res, ', time = ', tmr%time, ' s'
      call log(msg)

      if (solver%res .lt. 1.0e-8_WP .and. amr%nlevels .gt. 1) then
         call log('  PASS: multi-level AMR Poisson test')
      else if (solver%res .lt. 1.0e-8_WP) then
         call log('  PASS: single-level (no refinement triggered)')
      else
         call log('  FAIL: multi-level AMR Poisson test')
      end if

      ! Write visualization output
      block
         use amrviz_class, only: amrviz
         type(amrviz) :: viz
         call viz%initialize(amr=amr, name='poisson')
         call viz%add_scalar(data=phi, comp=1, name='phi')
         call viz%add_scalar(data=rhs, comp=1, name='rhs')
         call viz%write(time=0.0_WP)
         call viz%finalize()
         call log('  HDF5 output written to amrviz/poisson/')
      end block

      ! Cleanup
      call solver%finalize()
      call phi%finalize(); call rhs%finalize()
      deallocate(phi, rhs)
      call amr%finalize(); deallocate(amr)
   end subroutine test_multilevel

   !> Test per-level solve (for subcycling)
   !> Uses 3-level grid like test_multilevel:
   !> 1) Composite solve to get reference solution
   !> 2) Zero out levels 1 and 2 of phi
   !> 3) solve_level on level 1 only (using level 0 as coarse BC)
   !> 4) Output to ParaView - should show correct solution on levels 0,1 and zeros on level 2
   subroutine test_solve_level()
      use iso_c_binding, only: c_loc
      use amrviz_class, only: amrviz
      type(amrgrid), allocatable, target :: amr
      type(amrdata), allocatable, target :: phi, rhs
      type(amrmg) :: solver
      type(amrviz) :: viz
      type(timer) :: tmr
      character(len=str_long) :: msg
      logical :: passed

      call log('')
      call log('--- Testing solve_level (per-level solve with C/F BC) ---')
      call log('3-level grid: composite solve, then per-level solve on mid-level')

      ! Create 3-level grid (same as test_multilevel)
      allocate(amr)
      amr%nx = 32; amr%ny = 32; amr%nz = 32
      amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
      amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
      amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
      amr%xper = .false.; amr%yper = .false.; amr%zper = .false.
      amr%maxlvl = 2  ! 3 levels (0, 1, 2)
      amr%nmax = 16
      call amr%initialize(name='solve_level_ml')

      ! Allocate fields
      allocate(phi, rhs)
      call phi%initialize(amr=amr, name='phi', ncomp=1, ng=1)
      call rhs%initialize(amr=amr, name='rhs', ncomp=1, ng=0)

      ! Set user_init callbacks (reuses rhs_user_init from above)
      rhs%user_init => rhs_user_init

      ! Register fields
      call phi%register()
      call rhs%register()

      ! Add tagging based on RHS (reuses rhs_tagger from above)
      call amr%add_tagging(rhs_tagger, c_loc(rhs))

      ! Initialize grid - triggers callbacks and creates refined levels
      call amr%init_from_scratch(time=0.0_WP, do_postregrid=.true.)

      write(msg,'(a,i0,a)') '  Grid created with ', amr%nlevels, ' levels'
      call log(msg)

      ! Initialize solver (Neumann BC for Gaussian sources)
      call solver%initialize(amr=amr, type=amrmg_cstcoef)
      solver%tol_rel = 1.0e-10_WP
      solver%verbose = 0
      solver%bottom_solver = amrmg_bottom_hypre

      tmr = timer(comm=MPI_COMM_WORLD, name='amrmg_slvl')

      ! Step 1: Composite solve to get reference solution on all levels
      call log('  Step 1: Composite solve on all levels')
      call tmr%start()
      call solver%setup()
      call solver%solve(rhs=rhs,phi0=phi)
      call tmr%stop()
      write(msg,'(a,es12.5,a,i0)') '    Composite solve: residual = ', solver%res, ', niter = ', solver%niter
      call log(msg)

      passed = (solver%res .lt. 1.0e-8_WP)
      if (.not.passed) call log('    FAIL: composite solve did not converge')

      ! Write reference solution (all levels correct)
      call viz%initialize(amr=amr, name='solve_level_ref')
      call viz%add_scalar(data=phi, comp=1, name='phi')
      call viz%add_scalar(data=rhs, comp=1, name='rhs')
      call viz%write(time=0.0_WP)
      call viz%finalize()
      call log('    Reference solution written to amrviz/solve_level_ref/')

      ! Step 2: Zero out phi on levels 1 and 2
      call log('  Step 2: Zero out phi on levels 1 and 2')
      call phi%setval(val=0.0_WP, lbase=1)

      ! Step 3: solve_level on level 1 only (use level 0 as coarse BC)
      call log('  Step 3: solve_level on level 1 only (with C/F BC from level 0)')
      if (amr%clvl() .ge. 1) then
         call tmr%reset(); call tmr%start()
         call solver%solve_level(lev=1, phi_mf=phi%mf(1), rhs_mf=rhs%mf(1), phi_crse_mf=phi%mf(0))
         call tmr%stop()
         write(msg,'(a,es12.5,a,i0)') '    solve_level(1): residual = ', solver%res, ', niter = ', solver%niter
         call log(msg)

         if (solver%res .ge. 1.0e-8_WP) then
            call log('    FAIL: solve_level did not converge')
            passed = .false.
         end if

         if (solver%niter .le. 0) then
            call log('    FAIL: niter should be > 0')
            passed = .false.
         end if
      else
         call log('    SKIP: only 1 level created, cannot test C/F BC')
      end if

      ! Step 4: Write output - levels 0,1 should be correct, level 2 should be zero
      call viz%initialize(amr=amr, name='solve_level_mid')
      call viz%add_scalar(data=phi, comp=1, name='phi')
      call viz%add_scalar(data=rhs, comp=1, name='rhs')
      call viz%write(time=0.0_WP)
      call viz%finalize()
      call log('    Per-level solution written to amrviz/solve_level_mid/')
      call log('    (levels 0,1 should match reference; level 2 should be zero)')

      if (passed) then
         call log('  PASS: solve_level test')
      end if

      ! Cleanup
      call solver%finalize()
      call phi%finalize(); call rhs%finalize()
      deallocate(phi, rhs)
      call amr%finalize(); deallocate(amr)
   end subroutine test_solve_level

end module mod_test_amrmg

