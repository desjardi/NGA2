!> Unified AMR Multigrid solver for elliptic problems
!> Supports both constant-coefficient (Poisson) and variable-coefficient problems.
!>   - cstcoef type: Uses amrex_poisson for ∇²ϕ = f
!>   - varcoef type: Uses amrex_abeclaplacian for (α·A - β·∇·(B∇))ϕ = f
module amrmg_class
   use precision,     only: WP
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_boxarray, amrex_distromap
   use amrex_lo_bctypes_module
   use amrex_linear_solver_module
   implicit none
   private

   ! Solver type constants (required at init)
   integer, parameter, public :: amrmg_cstcoef = 1  !< Constant coefficient (amrex_poisson)
   integer, parameter, public :: amrmg_varcoef = 2  !< Variable coefficient (amrex_abeclaplacian)

   ! Bottom solver types
   integer, parameter, public :: amrmg_bottom_smoother = 0
   integer, parameter, public :: amrmg_bottom_bicgstab = 1
   integer, parameter, public :: amrmg_bottom_cg       = 2
   integer, parameter, public :: amrmg_bottom_hypre    = 3
   integer, parameter, public :: amrmg_bottom_petsc    = 4

   ! Outer solver types
   integer, parameter, public :: amrmg_outer_mlmg        = 0  !< Standard MLMG (default)
   integer, parameter, public :: amrmg_outer_pcg_mlmg    = 1  !< PCG with MLMG as preconditioner

   ! Expose AMReX boundary condition types for user convenience
   integer, parameter, public :: amrmg_bc_interior    = amrex_lo_bogus      !< Interior (not at domain boundary)
   integer, parameter, public :: amrmg_bc_dirichlet   = amrex_lo_dirichlet  !< Dirichlet BC
   integer, parameter, public :: amrmg_bc_neumann     = amrex_lo_neumann    !< Neumann BC
   integer, parameter, public :: amrmg_bc_reflect_odd = amrex_lo_reflect_odd
   integer, parameter, public :: amrmg_bc_periodic    = amrex_lo_periodic   !< Periodic BC

   !> Unified AMR Multigrid solver
   type, public :: amrmg
      ! Solver type (required at initialize)
      integer :: type = -1  !< amrmg_cstcoef or amrmg_varcoef (no default)

      ! Equation coefficients (varcoef type only)
      real(WP) :: alpha = 0.0_WP  !< Scalar for A term
      real(WP) :: beta  = 1.0_WP  !< Scalar for ∇·(B∇) term

      ! Solver parameters
      real(WP) :: tol_rel = 1.0e-10_WP
      real(WP) :: tol_abs = 0.0_WP
      integer  :: max_iter = 200
      integer  :: verbose = 0
      integer  :: bottom_solver = amrmg_bottom_bicgstab
      integer  :: outer_solver  = amrmg_outer_mlmg      !< Outer solver type
      integer  :: maxorder = 2  !< BC/interpolation order

      ! Boundary conditions (set from grid periodicity at init)
      integer, dimension(3) :: lo_bc = amrmg_bc_interior
      integer, dimension(3) :: hi_bc = amrmg_bc_interior

      ! Output from solve
      real(WP) :: res   = 0.0_WP  !< Final residual
      integer  :: niter = 0       !< Number of iterations

      ! Semi-coarsening and hidden direction flags
      logical :: semicoarsen = .false.
      integer :: hidden_direction = -1  !< -1=none, 0=x, 1=y, 2=z

      ! Private internals
      class(amrgrid), pointer, private :: amr => null()
      logical, private :: setup_done = .false.
      ! Stored AMReX objects for reuse
      type(amrex_poisson), private :: poisson
      type(amrex_abeclaplacian), private :: abeclap
      type(amrex_multigrid), private :: multigrid
      ! Internal solution storage for incremental solves
      type(amrdata) :: sol
      ! PCG work arrays (only used when outer_solver == amrmg_outer_pcg_mlmg)
      type(amrdata) :: pcg_r   !< Residual r = b - A*x
      type(amrdata) :: pcg_z   !< Preconditioned residual z = M^{-1}*r
      type(amrdata) :: pcg_p   !< Search direction p
      type(amrdata) :: pcg_q   !< q = A*p
   contains
      procedure :: initialize
      procedure :: setup
      procedure :: solve
      procedure :: solve_level
      generic :: get_fluxes => get_fluxes_amrdata, get_fluxes_mfab
      procedure, private :: get_fluxes_amrdata
      procedure, private :: get_fluxes_mfab
      procedure :: destroy
      procedure :: print_short
      procedure :: print => print_long
      procedure :: log => log_solver
      procedure :: finalize
   end type amrmg

contains

   !> Initialize solver with grid and type
   subroutine initialize(this, amr, type)
      use messager, only: die, log
      use amrdata_class, only: interp_none
      implicit none
      class(amrmg), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      integer, intent(in) :: type  !< amrmg_cstcoef or amrmg_varcoef (required)

      ! Validate type
      if (type .ne. amrmg_cstcoef .and. type .ne. amrmg_varcoef) then
         call die('[amrmg initialize] Invalid type: must be amrmg_cstcoef or amrmg_varcoef')
      end if

      this%amr => amr
      this%type = type
      this%setup_done = .false.

      ! Set smart BC defaults based on grid periodicity
      if (amr%xper) then
         this%lo_bc(1) = amrmg_bc_periodic
         this%hi_bc(1) = amrmg_bc_periodic
      else
         this%lo_bc(1) = amrmg_bc_neumann
         this%hi_bc(1) = amrmg_bc_neumann
      end if
      if (amr%yper) then
         this%lo_bc(2) = amrmg_bc_periodic
         this%hi_bc(2) = amrmg_bc_periodic
      else
         this%lo_bc(2) = amrmg_bc_neumann
         this%hi_bc(2) = amrmg_bc_neumann
      end if
      if (amr%zper) then
         this%lo_bc(3) = amrmg_bc_periodic
         this%hi_bc(3) = amrmg_bc_periodic
      else
         this%lo_bc(3) = amrmg_bc_neumann
         this%hi_bc(3) = amrmg_bc_neumann
      end if

      ! Auto-detect hidden direction for single-cell directions
      if (amr%nx.eq.1) this%hidden_direction = 0
      if (amr%ny.eq.1) this%hidden_direction = 1
      if (amr%nz.eq.1) this%hidden_direction = 2
      if (count([amr%nx,amr%ny,amr%nz].eq.1).gt.1) call die('[amrmg] Multiple single-cell directions not supported by AMReX MLMG')

      ! Initialize internal solution storage
      call this%sol%initialize(amr,name='sol',ncomp=1,ng=1,interp=interp_none); call this%sol%register()

      ! PCG work amrdata - unregistered
      call this%pcg_r%initialize(amr,name='pcg_r',ncomp=1,ng=1,interp=interp_none)
      call this%pcg_z%initialize(amr,name='pcg_z',ncomp=1,ng=1,interp=interp_none)
      call this%pcg_p%initialize(amr,name='pcg_p',ncomp=1,ng=1,interp=interp_none)
      call this%pcg_q%initialize(amr,name='pcg_q',ncomp=1,ng=1,interp=interp_none)

      ! Log setup info
      if (type .eq. amrmg_cstcoef) then
         call log('[amrmg] Initialized constant-coefficient solver')
      else
         call log('[amrmg] Initialized variable-coefficient solver')
      end if
      if (this%hidden_direction.ge.0) then
         call log('[amrmg] Hidden direction enabled')
      end if
      if (this%semicoarsen) then
         call log('[amrmg] Semi-coarsening enabled')
      end if
   end subroutine initialize

   !> Setup solver - call after coefficients change or after regrid
   !> For varcoef type, bcoef fields must be provided
   subroutine setup(this, acoef, bcoef_x, bcoef_y, bcoef_z)
      use messager, only: die
      use amrex_interface, only: amrpoisson_build, amrabeclap_build
      implicit none
      class(amrmg), intent(inout) :: this
      type(amrex_multifab), dimension(0:), intent(in), optional :: acoef
      type(amrex_multifab), dimension(0:), intent(in), optional :: bcoef_x, bcoef_y, bcoef_z

      type(amrex_geometry), dimension(:), allocatable :: geom
      type(amrex_boxarray), dimension(:), allocatable :: ba
      type(amrex_distromap), dimension(:), allocatable :: dm
      integer :: lev, nlevs

      if (this%type .eq. -1) call die('[amrmg setup] Solver not initialized')
      if (this%setup_done) call this%destroy()

      ! Build typed arrays for operator construction
      nlevs = this%amr%clvl() + 1
      allocate(geom(0:this%amr%clvl()))
      allocate(ba(0:this%amr%clvl()))
      allocate(dm(0:this%amr%clvl()))
      do lev = 0, this%amr%clvl()
         geom(lev) = this%amr%geom(lev)
         ba(lev) = this%amr%get_boxarray(lev)
         dm(lev) = this%amr%get_distromap(lev)
      end do

      ! Build operator based on type
      select case (this%type)

       case (amrmg_cstcoef)
         call amrpoisson_build(this%poisson, nlevs, geom, ba, dm, this%semicoarsen, this%hidden_direction)
         call this%poisson%set_domain_bc(this%lo_bc, this%hi_bc)
         call amrex_multigrid_build(this%multigrid, this%poisson)
         call this%multigrid%set_verbose(this%verbose)
         call this%multigrid%set_max_iter(this%max_iter)
         call this%multigrid%set_bottom_solver(this%bottom_solver)

       case (amrmg_varcoef)
         varcoef_setup: block
            type(amrex_multifab) :: bcoef(3)
            call amrabeclap_build(this%abeclap, nlevs, geom, ba, dm, this%semicoarsen, this%hidden_direction)
            call this%abeclap%set_domain_bc(this%lo_bc, this%hi_bc)
            call this%abeclap%set_maxorder(this%maxorder)
            call this%abeclap%set_scalars(this%alpha, this%beta)
            ! Set A coefficient
            if (present(acoef)) then
               do lev = 0, this%amr%clvl()
                  call this%abeclap%set_acoeffs(lev, acoef(lev))
               end do
            end if
            ! Set B coefficients
            if (present(bcoef_x) .and. present(bcoef_y) .and. present(bcoef_z)) then
               do lev = 0, this%amr%clvl()
                  bcoef(1) = bcoef_x(lev)
                  bcoef(2) = bcoef_y(lev)
                  bcoef(3) = bcoef_z(lev)
                  call this%abeclap%set_bcoeffs(lev, bcoef)
               end do
            end if
            call amrex_multigrid_build(this%multigrid, this%abeclap)
            call this%multigrid%set_verbose(this%verbose)
            call this%multigrid%set_max_iter(this%max_iter)
            call this%multigrid%set_bottom_solver(this%bottom_solver)
         end block varcoef_setup
      end select

      ! PCG work arrays: only allocate if using PCG outer solver
      if (this%outer_solver.eq.amrmg_outer_pcg_mlmg) then
         do lev=0,this%amr%clvl()
            call this%pcg_r%reset_level(lev,ba(lev),dm(lev))
            call this%pcg_z%reset_level(lev,ba(lev),dm(lev))
            call this%pcg_p%reset_level(lev,ba(lev),dm(lev))
            call this%pcg_q%reset_level(lev,ba(lev),dm(lev))
         end do
      end if

      this%setup_done = .true.
   end subroutine setup

   !> Solve - uses stored operator and multigrid
   !> @param rhs Right-hand side field
   !> @param phi0 Optional initial guess
   subroutine solve(this, rhs, phi0)
      use messager, only: die
      class(amrmg), intent(inout) :: this
      type(amrdata), intent(in) :: rhs
      type(amrdata), intent(inout), optional :: phi0
      type(amrex_multifab), dimension(:), allocatable :: sol,rhsmf

      if (this%type.eq.-1) call die('[amrmg solve] Solver not initialized')
      if (.not.this%setup_done) call die('[amrmg solve] Solver not setup')

      ! Set initial guess
      if (present(phi0)) then
         call this%sol%copy(src=phi0)
      else
         call this%sol%setval(val=0.0_WP)
      end if

      ! Build solution/rhs arrays
      build_arrays: block
         integer :: lev
         allocate(sol(0:this%amr%clvl()))
         allocate(rhsmf(0:this%amr%clvl()))
         do lev=0,this%amr%clvl()
            sol(lev)=this%sol%mf(lev)
            rhsmf(lev)=rhs%mf(lev)
         end do
      end block build_arrays

      ! Set level BCs
      set_bcs: block
         integer :: lev
         select case (this%type)
          case (amrmg_cstcoef)
            do lev=0,this%amr%clvl()
               call this%poisson%set_level_bc(lev,sol(lev))
            end do
          case (amrmg_varcoef)
            do lev=0,this%amr%clvl()
               call this%abeclap%set_level_bc(lev,sol(lev))
            end do
         end select
      end block set_bcs

      ! Solve and get iteration count
      select case (this%outer_solver)

       case (amrmg_outer_mlmg)
         ! Standard MLMG solve
         this%res=this%multigrid%solve(sol,rhsmf,this%tol_rel,this%tol_abs)
         get_niters_mlmg: block
            use amrex_interface, only: amrmlmg_get_niters
            this%niter = amrmlmg_get_niters(this%multigrid%p)
         end block get_niters_mlmg

       case (amrmg_outer_pcg_mlmg)
         ! PCG outer loop with 1 MLMG V-cycle as preconditioner
         pcg_solve: block
            use amrex_interface, only: amrmlmg_dot_composite
            use iso_c_binding,   only: c_ptr
            real(WP) :: rho, rho_new, alpha, beta, pq, rnorm, rnorm0, rnorm_prev, pcg_dummy
            integer  :: iter, lev
            type(amrex_multifab), dimension(:), allocatable :: resmf, zmf, pmf, qmf
            type(c_ptr), dimension(:), allocatable :: mf1_ptrs, mf2_ptrs, ba_ptrs
            integer, dimension(:), allocatable :: rr_flat
            type(amrex_boxarray) :: ba_tmp
            integer :: nlevs

            nlevs = this%amr%clvl() + 1

            ! Build helper multifab arrays for pcg vectors
            allocate(resmf(0:this%amr%clvl()), zmf(0:this%amr%clvl()))
            allocate(pmf(0:this%amr%clvl()),   qmf(0:this%amr%clvl()))
            do lev=0,this%amr%clvl()
               resmf(lev)=this%pcg_r%mf(lev)
               zmf(lev)=this%pcg_z%mf(lev)
               pmf(lev)=this%pcg_p%mf(lev)
               qmf(lev)=this%pcg_q%mf(lev)
            end do

            ! Prepare pointer arrays for composite dot
            allocate(mf1_ptrs(0:this%amr%clvl()), mf2_ptrs(0:this%amr%clvl()))
            allocate(ba_ptrs(0:this%amr%clvl()))
            allocate(rr_flat(3*nlevs))
            do lev=0,this%amr%clvl()
               ba_tmp = this%amr%get_boxarray(lev); ba_ptrs(lev) = ba_tmp%p
               if (lev .lt. this%amr%clvl()) then
                  rr_flat(3*lev+1) = this%amr%rrefx(lev)
                  rr_flat(3*lev+2) = this%amr%rrefy(lev)
                  rr_flat(3*lev+3) = this%amr%rrefz(lev)
               else
                  rr_flat(3*lev+1) = 1; rr_flat(3*lev+2) = 1; rr_flat(3*lev+3) = 1
               end if
            end do

            ! --- Set MLMG to fixed 1 V-cycle for preconditioning ---
            call this%multigrid%set_fixed_iter(1)

            ! --- r = b - A*x0 (using comp_residual) ---
            call this%multigrid%comp_residual(resmf, sol, rhsmf)

            ! --- z = M^{-1} r (1 V-cycle): zero first, then set BCs, then solve ---
            call this%pcg_z%setval(0.0_WP)
            do lev=0,this%amr%clvl()
               select case (this%type)
                case (amrmg_cstcoef); call this%poisson%set_level_bc(lev, zmf(lev))
                case (amrmg_varcoef); call this%abeclap%set_level_bc(lev, zmf(lev))
               end select
            end do
            pcg_dummy = this%multigrid%solve(zmf, resmf, 0.0_WP, 0.0_WP)

            ! --- p = z ---
            call this%pcg_p%copy(src=this%pcg_z)

            ! --- rho = <r, z> ---
            do lev=0,this%amr%clvl()
               mf1_ptrs(lev) = this%pcg_r%mf(lev)%p
               mf2_ptrs(lev) = this%pcg_z%mf(lev)%p
            end do
            rho = amrmlmg_dot_composite(mf1_ptrs, mf2_ptrs, ba_ptrs, rr_flat, nlevs)

            ! --- Initial residual norm for convergence test ---
            do lev=0,this%amr%clvl()
               mf1_ptrs(lev) = this%pcg_r%mf(lev)%p
               mf2_ptrs(lev) = this%pcg_r%mf(lev)%p
            end do
            rnorm0 = sqrt(amrmlmg_dot_composite(mf1_ptrs, mf2_ptrs, ba_ptrs, rr_flat, nlevs))
            if (rnorm0.eq.0.0_WP) then
               this%res = 0.0_WP; this%niter = 0
            else
               ! --- PCG with safeguarded beta (restart when residual grows) ---
               rnorm_prev = rnorm0
               pcg_iter: do iter = 1, this%max_iter
                  ! --- q = A*p ---
                  call this%pcg_z%setval(0.0_WP)
                  call this%multigrid%comp_residual(qmf, pmf, zmf)  ! q = -A*p
                  call this%pcg_q%mult(-1.0_WP)                     ! q = A*p

                  ! --- alpha = rho / <p, q> ---
                  do lev=0,this%amr%clvl()
                     mf1_ptrs(lev) = this%pcg_p%mf(lev)%p
                     mf2_ptrs(lev) = this%pcg_q%mf(lev)%p
                  end do
                  pq = amrmlmg_dot_composite(mf1_ptrs, mf2_ptrs, ba_ptrs, rr_flat, nlevs)
                  if (abs(pq).lt.tiny(pq)) exit
                  alpha = rho / pq

                  ! --- x += alpha*p, r -= alpha*q ---
                  call this%sol%saxpy(a=alpha, src=this%pcg_p)
                  call this%pcg_r%saxpy(a=-alpha, src=this%pcg_q)

                  ! --- Convergence check ---
                  do lev=0,this%amr%clvl()
                     mf1_ptrs(lev) = this%pcg_r%mf(lev)%p
                     mf2_ptrs(lev) = this%pcg_r%mf(lev)%p
                  end do
                  rnorm = sqrt(amrmlmg_dot_composite(mf1_ptrs, mf2_ptrs, ba_ptrs, rr_flat, nlevs))
                  this%niter = iter; this%res = rnorm
                  if (this%verbose .gt. 1 .and. this%amr%amRoot) write(*,'("  [PCG-MLMG] iter ",i0," res/res0=",es12.5)') iter, rnorm/rnorm0
                  if (rnorm.le.this%tol_abs .or. rnorm/rnorm0.le.this%tol_rel) exit

                  ! --- z = M^{-1} r ---
                  call this%pcg_z%setval(0.0_WP)
                  do lev=0,this%amr%clvl()
                     select case (this%type)
                      case (amrmg_cstcoef); call this%poisson%set_level_bc(lev, zmf(lev))
                      case (amrmg_varcoef); call this%abeclap%set_level_bc(lev, zmf(lev))
                     end select
                  end do
                  pcg_dummy = this%multigrid%solve(zmf, resmf, 0.0_WP, 0.0_WP)

                  ! --- rho_new = <r, z> ---
                  do lev=0,this%amr%clvl()
                     mf1_ptrs(lev) = this%pcg_r%mf(lev)%p
                     mf2_ptrs(lev) = this%pcg_z%mf(lev)%p
                  end do
                  rho_new = amrmlmg_dot_composite(mf1_ptrs, mf2_ptrs, ba_ptrs, rr_flat, nlevs)
                  if (abs(rho).lt.tiny(rho)) exit
                  beta = rho_new / rho
                  ! Safeguard: reset beta if residual grew (prevents Krylov direction corruption)
                  if (rnorm .gt. rnorm_prev) then
                     beta = 0.0_WP
                     if (this%verbose .gt. 1 .and. this%amr%amRoot) write(*,'("  [PCG-MLMG] restart at iter ",i0)') iter
                  end if
                  rnorm_prev = rnorm
                  rho = rho_new

                  ! --- p = z + beta*p ---
                  call this%pcg_p%lincomb(a=1.0_WP, src1=this%pcg_z, b=beta, src2=this%pcg_p)
               end do pcg_iter
               ! Fill sol ghost cells: applyBC side effect of comp_residual
               ! PCG's saxpy updates leave ghost cells inconsistent; this corrects them.
               call this%multigrid%comp_residual(resmf, sol, rhsmf)

            end if  ! rnorm0 > 0

            ! Restore MLMG to normal mode
            call this%multigrid%set_fixed_iter(0)
            call this%multigrid%set_max_iter(this%max_iter)
            deallocate(resmf, zmf, pmf, qmf, mf1_ptrs, mf2_ptrs, ba_ptrs, rr_flat)
         end block pcg_solve

      end select
   end subroutine solve

   !> Level-by-level solve (for subcycling)
   !> Builds single-level operator on-the-fly (not stored)
   subroutine solve_level(this, lev, phi_mf, rhs_mf, phi_crse_mf, acoef_mf, bcoef_x_mf, bcoef_y_mf, bcoef_z_mf)
      use messager, only: die, log
      use string, only: str_long
      use amrex_interface, only: amrlinop_set_coarse_fine_bc, amrmlmg_get_niters, &
      &                          amrpoisson_build, amrabeclap_build
      class(amrmg), intent(inout) :: this
      integer, intent(in) :: lev
      type(amrex_multifab), intent(inout) :: phi_mf
      type(amrex_multifab), intent(in) :: rhs_mf
      type(amrex_multifab), intent(in), optional :: phi_crse_mf
      type(amrex_multifab), intent(in), optional :: acoef_mf
      type(amrex_multifab), intent(in), optional :: bcoef_x_mf, bcoef_y_mf, bcoef_z_mf

      ! Single-level arrays
      type(amrex_geometry) :: geom(0:0)
      type(amrex_boxarray) :: ba(0:0)
      type(amrex_distromap) :: dm(0:0)
      type(amrex_multifab) :: sol(0:0), rhsmf(0:0)
      integer, dimension(3) :: rref

      ! Validate inputs
      if (this%type .eq. -1) call die('[amrmg solve_level] Solver not initialized')
      if (lev .gt. 0 .and. .not.present(phi_crse_mf)) call die('[amrmg solve_level] phi_crse_mf required for lev>0')

      ! Build single-level arrays
      geom(0) = this%amr%geom(lev)
      ba(0) = this%amr%get_boxarray(lev)
      dm(0) = this%amr%get_distromap(lev)
      sol(0) = phi_mf
      rhsmf(0) = rhs_mf
      if (lev .gt. 0) rref = [this%amr%rrefx(lev-1),this%amr%rrefy(lev-1),this%amr%rrefz(lev-1)]

      ! Solve based on operator type
      select case (this%type)

       case (amrmg_cstcoef)
         poisson_solve: block
            type(amrex_poisson) :: linop
            type(amrex_multigrid) :: mlmg
            call amrpoisson_build(linop, 1, geom, ba, dm, this%semicoarsen, this%hidden_direction)
            call linop%set_domain_bc(this%lo_bc, this%hi_bc)
            ! Set C/F BC if on refined level
            if (lev .gt. 0) call amrlinop_set_coarse_fine_bc(linop%p, phi_crse_mf%p, rref)
            call linop%set_level_bc(0, sol(0))
            ! Solve
            call amrex_multigrid_build(mlmg, linop)
            call mlmg%set_verbose(this%verbose)
            call mlmg%set_max_iter(this%max_iter)
            call mlmg%set_bottom_solver(this%bottom_solver)
            this%res = mlmg%solve(sol, rhsmf, this%tol_rel, this%tol_abs)
            ! Get iteration count before cleanup
            this%niter = amrmlmg_get_niters(mlmg%p)
            call amrex_multigrid_destroy(mlmg)
            call amrex_poisson_destroy(linop)
         end block poisson_solve

       case (amrmg_varcoef)
         abeclap_solve: block
            type(amrex_abeclaplacian) :: linop
            type(amrex_multigrid) :: mlmg
            type(amrex_multifab) :: bcoef(3)
            call amrabeclap_build(linop, 1, geom, ba, dm, this%semicoarsen, this%hidden_direction)
            call linop%set_domain_bc(this%lo_bc, this%hi_bc)
            call linop%set_maxorder(this%maxorder)
            call linop%set_scalars(this%alpha, this%beta)
            ! Set coefficients
            if (present(acoef_mf)) call linop%set_acoeffs(0, acoef_mf)
            if (present(bcoef_x_mf) .and. present(bcoef_y_mf) .and. present(bcoef_z_mf)) then
               bcoef(1) = bcoef_x_mf
               bcoef(2) = bcoef_y_mf
               bcoef(3) = bcoef_z_mf
               call linop%set_bcoeffs(0, bcoef)
            end if
            ! Set C/F BC if on refined level
            if (lev .gt. 0) call amrlinop_set_coarse_fine_bc(linop%p, phi_crse_mf%p, rref)
            call linop%set_level_bc(0, sol(0))
            ! Solve
            call amrex_multigrid_build(mlmg, linop)
            call mlmg%set_verbose(this%verbose)
            call mlmg%set_max_iter(this%max_iter)
            call mlmg%set_bottom_solver(this%bottom_solver)
            this%res = mlmg%solve(sol, rhsmf, this%tol_rel, this%tol_abs)
            ! Get iteration count before cleanup
            this%niter = amrmlmg_get_niters(mlmg%p)
            call amrex_multigrid_destroy(mlmg)
            call amrex_abeclaplacian_destroy(linop)
         end block abeclap_solve
      end select

      ! Log result
      log_result: block
         character(len=str_long) :: message
         if (this%amr%amRoot) then
            write(message,'("[amrmg solve_level] lev=",i0," res=",es12.5)') lev, this%res
            call log(message)
         end if
      end block log_result
   end subroutine solve_level

   !> Get face-centered fluxes using solver's C/F-consistent stencils (amrdata version)
   !> Always uses this%sol (the solution from the last solve call).
   !> NOTE: Do NOT attempt to pass a different phi - MLMG boundary registers
   !>       are tied to the solution that was just solved; using a different
   !>       field will produce incorrect C/F flux values.
   subroutine get_fluxes_amrdata(this, flux_x, flux_y, flux_z)
      use iso_c_binding, only: c_ptr
      use messager, only: die
      use amrex_interface, only: amrmlmg_get_fluxes
      class(amrmg), intent(in) :: this
      type(amrdata), intent(inout) :: flux_x, flux_y, flux_z
      type(c_ptr), dimension(:), allocatable :: sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs
      integer :: lev, nlevs
      if (.not. this%setup_done) call die('[amrmg get_fluxes] Solver not setup')
      nlevs = this%amr%clvl() + 1
      ! Build pointer arrays
      allocate(sol_ptrs(0:this%amr%clvl()))
      allocate(fx_ptrs(0:this%amr%clvl()))
      allocate(fy_ptrs(0:this%amr%clvl()))
      allocate(fz_ptrs(0:this%amr%clvl()))
      ! Always use this%sol - the solution from the last solve()
      do lev = 0, this%amr%clvl()
         sol_ptrs(lev) = this%sol%mf(lev)%p
         fx_ptrs(lev) = flux_x%mf(lev)%p
         fy_ptrs(lev) = flux_y%mf(lev)%p
         fz_ptrs(lev) = flux_z%mf(lev)%p
      end do
      ! Call C++ wrapper - works for both Poisson and ABecLaplacian
      call amrmlmg_get_fluxes(this%multigrid%p, sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs, nlevs)
      deallocate(sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs)
   end subroutine get_fluxes_amrdata

   !> Get face-centered fluxes using solver's C/F-consistent stencils (multifab array version)
   !> Always uses this%sol (the solution from the last solve call).
   !> NOTE: Do NOT attempt to pass a different phi - MLMG boundary registers
   !>       are tied to the solution that was just solved; using a different
   !>       field will produce incorrect C/F flux values.
   subroutine get_fluxes_mfab(this, flux_x, flux_y, flux_z)
      use iso_c_binding, only: c_ptr
      use messager, only: die
      use amrex_interface, only: amrmlmg_get_fluxes
      class(amrmg), intent(in) :: this
      type(amrex_multifab), dimension(0:), intent(inout) :: flux_x, flux_y, flux_z
      type(c_ptr), dimension(:), allocatable :: sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs
      integer :: lev, nlevs
      if (.not. this%setup_done) call die('[amrmg get_fluxes] Solver not setup')
      nlevs = this%amr%clvl() + 1
      ! Build pointer arrays
      allocate(sol_ptrs(0:this%amr%clvl()))
      allocate(fx_ptrs(0:this%amr%clvl()))
      allocate(fy_ptrs(0:this%amr%clvl()))
      allocate(fz_ptrs(0:this%amr%clvl()))
      ! Always use this%sol - the solution from the last solve()
      do lev = 0, this%amr%clvl()
         sol_ptrs(lev) = this%sol%mf(lev)%p
         fx_ptrs(lev) = flux_x(lev)%p
         fy_ptrs(lev) = flux_y(lev)%p
         fz_ptrs(lev) = flux_z(lev)%p
      end do
      ! Call C++ wrapper - works for both Poisson and ABecLaplacian
      call amrmlmg_get_fluxes(this%multigrid%p, sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs, nlevs)
      deallocate(sol_ptrs, fx_ptrs, fy_ptrs, fz_ptrs)
   end subroutine get_fluxes_mfab

   !> Destroy solver internals - call before setup when operator changes
   subroutine destroy(this)
      class(amrmg), intent(inout) :: this
      if (this%setup_done) then
         call amrex_multigrid_destroy(this%multigrid)
         select case (this%type)
          case (amrmg_cstcoef)
            call amrex_poisson_destroy(this%poisson)
          case (amrmg_varcoef)
            call amrex_abeclaplacian_destroy(this%abeclap)
         end select
      end if
      this%setup_done = .false.
   end subroutine destroy

   !> Log solver info
   subroutine log_solver(this)
      use string,   only: str_long
      use messager, only: log
      class(amrmg), intent(in) :: this
      character(len=str_long) :: message
      if (this%amr%amRoot) then
         write(message,'("AMR Multigrid solver for grid [",a,"]")') trim(this%amr%name); call log(message)
         if (this%type .eq. amrmg_cstcoef) then
            write(message,'(" >     type = constant-coefficient (Poisson)")'); call log(message)
         else if (this%type .eq. amrmg_varcoef) then
            write(message,'(" >     type = variable-coefficient (ABecLaplacian)")'); call log(message)
            write(message,'(" >    alpha = ",es12.5)') this%alpha; call log(message)
            write(message,'(" >     beta = ",es12.5)') this%beta; call log(message)
         end if
         write(message,'(" >   it/max = ",i0,"/",i0)') this%niter,this%max_iter; call log(message)
         write(message,'(" >      res = ",es12.5)') this%res; call log(message)
         write(message,'(" >  tol_rel = ",es12.5)') this%tol_rel; call log(message)
      end if
   end subroutine log_solver

   !> Print solver info to screen
   subroutine print_long(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(amrmg), intent(in) :: this
      if (this%amr%amRoot) then
         write(output_unit,'("AMR Multigrid solver for grid [",a,"]")') trim(this%amr%name)
         if (this%type .eq. amrmg_cstcoef) then
            write(output_unit,'(" >     type = constant-coefficient (Poisson)")')
         else if (this%type .eq. amrmg_varcoef) then
            write(output_unit,'(" >     type = variable-coefficient (ABecLaplacian)")')
            write(output_unit,'(" >    alpha = ",es12.5)') this%alpha
            write(output_unit,'(" >     beta = ",es12.5)') this%beta
         end if
         write(output_unit,'(" >   it/max = ",i0,"/",i0)') this%niter,this%max_iter
         write(output_unit,'(" >      res = ",es12.5)') this%res
         write(output_unit,'(" >  tol_rel = ",es12.5)') this%tol_rel
      end if
   end subroutine print_long

   !> Short print of solver info
   subroutine print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      class(amrmg), intent(in) :: this
      if (this%amr%amRoot) then
         write(output_unit,'("AMR MG solver [",a16,"] -> it/max = ",i3,"/",i3," res = ",es12.5)') &
            trim(this%amr%name),this%niter,this%max_iter,this%res
      end if
   end subroutine print_short

   !> Finalize and release resources
   subroutine finalize(this)
      class(amrmg), intent(inout) :: this
      if (this%setup_done) call this%destroy()
      call this%sol%finalize()
      call this%pcg_r%finalize()
      call this%pcg_z%finalize()
      call this%pcg_p%finalize()
      call this%pcg_q%finalize()
      nullify(this%amr)
      this%type = -1
   end subroutine finalize

end module amrmg_class
