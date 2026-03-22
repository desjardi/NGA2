!> AMR-aware Lagrangian particle tracking solver
!> Mirrors lpt_class capabilities on an AMReX AMR hierarchy.
!> Particle communication and sorting handled by AmrParticleContainer<14,1>.
module amrlpt_class
   use precision, only: WP
   use string, only: str_medium
   use amrgrid_class, only: amrgrid
   use iso_c_binding
   implicit none
   private

   ! Public exports
   public :: amrlpt,part

   ! Particle struct layout constants (must match #define in amrlpt_wrapper.cpp)
   integer, parameter, public :: AMRLPT_NREAL=14  !< extra reals per particle
   integer, parameter, public :: AMRLPT_NINT =1   !< extra ints  per particle

   ! -----------------------------------------------------------------------
   ! Particle struct -- must match C++ Particle<14,1> memory layout exactly:
   !   pos[3]    (pos, managed by AMReX)
   !   rdata[14] (d, vel[3], angVel[3], Acol[3], Tcol[3], dt)
   !   idcpu     (packed id+cpu, private)
   !   idata[1]  (flag)
   ! -----------------------------------------------------------------------
   type, bind(C), public :: part
      !> AMReX position (physical coordinates)
      real(c_double) :: pos(3)
      !> Extra reals (NStructReal=14)
      real(c_double) :: d               !< Particle diameter
      real(c_double) :: vel(3)          !< Particle velocity
      real(c_double) :: angVel(3)       !< Angular velocity
      real(c_double) :: Acol(3)         !< Collision acceleration
      real(c_double) :: Tcol(3)         !< Collision torque
      real(c_double) :: dt              !< Particle sub-timestep
      !> Packed id (39 bits, sign=valid) + cpu (24 bits)  [AMReX internal]
      integer(c_int64_t), private :: idcpu
      !> Extra int (NStructInt=1)
      integer(c_int) :: flag            !< 0=active, 1=mark for removal
   end type part

   !> C interface: bind(C) declarations to amrlpt_wrapper.cpp
   interface

      subroutine amrlpt_new_pc(pc, amrcore) bind(c)
         import :: c_ptr
         type(c_ptr) :: pc
         type(c_ptr), value :: amrcore
      end subroutine

      subroutine amrlpt_delete_pc(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_redistribute(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine

      subroutine amrlpt_get_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_add_particle_i(pc,lev,grid,tile,p) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,p
         integer(c_int), value :: lev,grid,tile
      end subroutine

      subroutine amrlpt_get_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine

      subroutine amrlpt_set_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t), value :: id
      end subroutine

      subroutine amrlpt_get_cpu(cpu) bind(c)
         import :: c_int
         integer(c_int) :: cpu
      end subroutine

      subroutine amrlpt_set_particle_id(id,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t), value :: id
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_set_particle_cpu(cpu,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int), value :: cpu
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_particle_is_valid(valid,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int) :: valid
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_total_np(pc,np) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_write(pc,path,is_chk) bind(c)
         import :: c_ptr,c_char,c_int
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
         integer(c_int), value :: is_chk
      end subroutine

      subroutine amrlpt_read(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
      end subroutine

   end interface

   !> AMR LPT solver type
   type :: amrlpt

      !> Associated AMR grid
      class(amrgrid), pointer :: amr=>null()

      !> AMReX AmrParticleContainer<14,1> opaque handle
      type(c_ptr) :: pc=c_null_ptr

      !> Solver name
      character(len=str_medium) :: name='UNNAMED_AMRLPT'

      !> Global particle count
      integer(c_int64_t) :: np=0

      !> Physics parameters
      real(WP) :: rho                                   !< Particle material density
      real(WP), dimension(3) :: gravity=0.0_WP          !< Acceleration of gravity
      integer :: nstep=1                                !< Substeps per timestep
      character(len=str_medium) :: drag_model='Tenneti' !< Drag model

      !> Collision parameters
      real(WP) :: tau_col=1.0e-3_WP             !< Collision time scale
      real(WP) :: e_n=1.0_WP                    !< Normal restitution
      real(WP) :: e_w=1.0_WP                    !< Wall restitution
      real(WP) :: mu_f=0.0_WP                   !< Friction coefficient
      real(WP) :: clip_col=0.2_WP               !< Max overlap fraction

      !> CFL numbers
      real(WP) :: CFLp_x=0.0_WP, CFLp_y=0.0_WP, CFLp_z=0.0_WP
      real(WP) :: CFL_col=0.0_WP

      !> Monitoring data
      real(WP) :: dmin,dmax,dmean,dvar
      real(WP) :: Umin,Umax,Umean,Uvar
      real(WP) :: Vmin,Vmax,Vmean,Vvar
      real(WP) :: Wmin,Wmax,Wmean,Wvar
      real(WP) :: VFmin,VFmax,VFmean,VFvar
      integer  :: np_new=0, np_out=0
      real(WP) :: vp_new=0.0_WP, vp_out=0.0_WP, vp_tot=0.0_WP
      integer  :: ncol=0

   contains
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Physics procedures
      procedure :: advance       !< Advance particle ODEs one timestep
      procedure :: get_rhs       !< Compute drag RHS + optimal sub-dt
      procedure :: update_VF     !< Compute particle volume fraction field
      procedure :: get_cfl       !< Compute particle CFL numbers
      ! Utilities
      procedure :: redistribute  !< Call AMReX redistribute
      ! Print solver info
      procedure :: get_info
      procedure :: print
      ! Checkpoint I/O
      procedure :: read
      procedure :: write
   end type amrlpt

contains

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize amrlpt solver
   subroutine initialize(this,amr,name)
      implicit none
      class(amrlpt), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), optional :: name
      ! Assign solver name
      if (present(name)) this%name=trim(adjustl(name))
      ! Point to amr
      this%amr=>amr
      ! Create particle container
      call amrlpt_new_pc(this%pc,this%amr%amrcore)
      ! Print out info
      call this%print()
   end subroutine initialize

   !> Finalize: destroy particle container and release grid pointer
   subroutine finalize(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call amrlpt_delete_pc(this%pc)
      this%pc=c_null_ptr
      nullify(this%amr)
   end subroutine finalize

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Advance all particles by dt using RK2 with substepping
   !> U_mf, V_mf, W_mf: staggered MAC velocity multifabs (lev, with ghosts)
   !> rho_mf, visc_mf:  cell-centred multifabs (lev, with ghosts)
   !> lev:  the level particles currently live on (after redistribute)
   !> geom: geometry for that level (plo, dx)
   !> srcU/V/W_mf: optional two-way coupling momentum source (cell-centred)
   subroutine advance(this,dt,U_mf,V_mf,W_mf,rho_mf,visc_mf,lev,geom,srcU_mf,srcV_mf,srcW_mf)
      use amrex_amr_module, only: amrex_multifab,amrex_geometry,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INT,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrex_multifab), intent(in) :: U_mf,V_mf,W_mf
      type(amrex_multifab), intent(in) :: rho_mf,visc_mf
      integer, intent(in) :: lev
      type(amrex_geometry), intent(in) :: geom
      type(amrex_multifab), intent(inout), optional :: srcU_mf,srcV_mf,srcW_mf

      type(amrex_mfiter) :: mfi
      type(c_ptr) :: dp
      type(part), pointer :: p(:)
      integer(c_int64_t) :: np_tile
      integer :: i, ierr
      real(WP) :: mydt, dt_done, Ip
      real(WP), dimension(3) :: acc, dmom
      type(part) :: myp, pold

      real(WP), dimension(:,:,:,:), contiguous, pointer :: Uarr, Varr, Warr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rhoarr, viscarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: sUarr, sVarr, sWarr
      real(WP) :: plo(3), dx(3)

      plo = geom%get_physical_location([geom%domain%lo(0), geom%domain%lo(1), geom%domain%lo(2)])
      dx  = geom%dx

      ! Zero source terms
      if (present(srcU_mf)) call srcU_mf%setval(0.0_WP)
      if (present(srcV_mf)) call srcV_mf%setval(0.0_WP)
      if (present(srcW_mf)) call srcW_mf%setval(0.0_WP)

      this%np_out = 0
      this%vp_out = 0.0_WP

      ! MFIter over the level (no tiling — particles sorted by grid, not tile)
      call amrex_mfiter_build(mfi, U_mf, tiling=.false.)
      do while (mfi%valid())

         ! Bind data arrays for this FAB (includes ghost cells)
         Uarr   => U_mf%dataptr(mfi)
         Varr   => V_mf%dataptr(mfi)
         Warr   => W_mf%dataptr(mfi)
         rhoarr => rho_mf%dataptr(mfi)
         viscarr=> visc_mf%dataptr(mfi)
         if (present(srcU_mf)) sUarr => srcU_mf%dataptr(mfi)
         if (present(srcV_mf)) sVarr => srcV_mf%dataptr(mfi)
         if (present(srcW_mf)) sWarr => srcW_mf%dataptr(mfi)

         ! Get particles on this tile
         call amrlpt_get_particles_mfi(this%pc, lev, mfi%p, dp, np_tile)
         if (np_tile > 0) then
            call c_f_pointer(dp, p, [np_tile])

            do i = 1, int(np_tile)
               ! Skip dead / inactive particles
               if (p(i)%flag == 1) cycle

               myp = p(i)
               dt_done = 0.0_WP
               Ip = 0.1_WP * myp%d**2  ! moment of inertia / mass for sphere

               ! Substep loop (Euler-midpoint)
               do while (dt_done < dt)
                  mydt = min(myp%dt, dt - dt_done)
                  if (mydt <= 0.0_WP) mydt = dt - dt_done  ! first step: dt not yet set
                  pold = myp

                  ! --- Euler predictor ---
                  call this%get_rhs(plo, dx, Uarr, Varr, Warr, rhoarr, viscarr, &
                                    myp, acc, myp%dt)
                  mydt = min(myp%dt, dt - dt_done)
                  myp%pos = pold%pos + 0.5_WP*mydt*myp%vel
                  myp%vel = pold%vel + 0.5_WP*mydt*(acc + this%gravity + myp%Acol)
                  myp%angVel = pold%angVel + 0.5_WP*mydt*myp%Tcol/Ip

                  ! --- Midpoint corrector ---
                  call this%get_rhs(plo, dx, Uarr, Varr, Warr, rhoarr, viscarr, &
                                    myp, acc, myp%dt)
                  myp%pos = pold%pos + mydt*myp%vel
                  myp%vel = pold%vel + mydt*(acc + this%gravity + myp%Acol)
                  myp%angVel = pold%angVel + mydt*myp%Tcol/Ip

                  ! Two-way coupling: deposit momentum change to mesh
                  dmom = mydt * acc * this%rho * Pi/6.0_WP * myp%d**3
                  if (present(srcU_mf)) &
                     call deposit_scalar(-dmom(1), myp%pos, plo, dx, sUarr)
                  if (present(srcV_mf)) &
                     call deposit_scalar(-dmom(2), myp%pos, plo, dx, sVarr)
                  if (present(srcW_mf)) &
                     call deposit_scalar(-dmom(3), myp%pos, plo, dx, sWarr)

                  dt_done = dt_done + mydt
               end do

               ! Enforce periodicity (geometry handles actual periodic wrapping
               ! during redistribute; just clamp here for clarity)
               if (geom%is_periodic(0)) &
                  myp%pos(1) = myp%pos(1) - floor((myp%pos(1)-plo(1)) / &
                               (dx(1)*real(geom%domain%hi(0)-geom%domain%lo(0)+1,WP))) * &
                               (dx(1)*real(geom%domain%hi(0)-geom%domain%lo(0)+1,WP))

               ! Flag particles leaving domain (non-periodic faces)
               if (.not.geom%is_periodic(0)) then
                  if (myp%pos(1) < plo(1) .or. &
                      myp%pos(1) > plo(1) + dx(1)*(geom%domain%hi(0)-geom%domain%lo(0)+1)) &
                     myp%flag = 1
               end if
               if (.not.geom%is_periodic(1)) then
                  if (myp%pos(2) < plo(2) .or. &
                      myp%pos(2) > plo(2) + dx(2)*(geom%domain%hi(1)-geom%domain%lo(1)+1)) &
                     myp%flag = 1
               end if
               if (.not.geom%is_periodic(2)) then
                  if (myp%pos(3) < plo(3) .or. &
                      myp%pos(3) > plo(3) + dx(3)*(geom%domain%hi(2)-geom%domain%lo(2)+1)) &
                     myp%flag = 1
               end if

               if (myp%flag == 1) then
                  this%np_out = this%np_out + 1
                  this%vp_out = this%vp_out + Pi/6.0_WP*myp%d**3
               end if

               ! Write back (preserving flag)
               p(i) = myp
            end do
         end if

         call mfi%next()
      end do
      call amrex_mfiter_destroy(mfi)

      ! AMReX handles redistribution (moves particles to correct rank/level/tile)
      call amrlpt_redistribute(this%pc, 0, -1, 0)

      ! Remove flagged particles (AMReX does this with remove_negative during
      ! redistribute when id<0; we mark by flag, so compact manually)
      ! TODO: expose amrlpt_remove_flagged in wrapper

      ! Divide src by cell volume and sum at boundaries
      ! (user is responsible for calling mfab FillBoundary / syncsum)

      ! Update global particle count
      call amrlpt_total_np(this%pc, this%np)

   end subroutine advance

   ! -----------------------------------------------------------------------
   !> Compute RHS of particle ODE: drag acceleration + optimal sub-dt.
   !> Interpolates fluid U,V,W (MAC), rho, visc (cell-centred) to p%pos.
   ! -----------------------------------------------------------------------
   subroutine get_rhs(this, plo, dx, Uarr, Varr, Warr, rhoarr, viscarr, p, acc, opt_dt)
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in)  :: plo(3), dx(3)
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: Uarr, Varr, Warr
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: rhoarr, viscarr
      type(part), intent(in)  :: p
      real(WP), intent(out) :: acc(3)
      real(WP), intent(out) :: opt_dt

      real(WP) :: fvel(3), frho, fvisc
      real(WP) :: Re, tau, corr, b1, b2
      real(WP), parameter :: pVF = 0.0_WP  !< dilute limit (VF coupling TODO)
      real(WP), parameter :: fVF = 1.0_WP

      ! Interpolate fluid velocity (MAC staggered)
      call interp_mac(p%pos, plo, dx, Uarr, Varr, Warr, fvel)
      ! Interpolate rho and visc (cell-centred)
      frho  = interp_cc(p%pos, plo, dx, rhoarr)
      fvisc = interp_cc(p%pos, plo, dx, viscarr) + epsilon(1.0_WP)

      ! Drag correction
      select case(trim(this%drag_model))
      case('None','none')
         corr = epsilon(1.0_WP)
      case('Stokes')
         corr = 1.0_WP
      case('Schiller-Naumann','SN')
         Re   = frho*norm2(p%vel-fvel)*p%d/fvisc + epsilon(1.0_WP)
         corr = 1.0_WP + 0.15_WP*Re**(0.687_WP)
      case('Tenneti')
         Re   = fVF*frho*norm2(p%vel-fvel)*p%d/fvisc + epsilon(1.0_WP)
         b1   = 5.81_WP*pVF/fVF**3 + 0.48_WP*pVF**(1.0_WP/3.0_WP)/fVF**4
         b2   = pVF**3*Re*(0.95_WP + 0.61_WP*pVF**3/fVF**2)
         corr = fVF*((1.0_WP+0.15_WP*Re**(0.687_WP))/fVF**3+b1+b2)
      case('Beetstra')
         Re   = fVF*frho*norm2(p%vel-fvel)*p%d/fvisc + epsilon(1.0_WP)
         b1   = 10.0_WP*pVF/fVF**2 + fVF**2*(1.0_WP+1.5_WP*sqrt(pVF))
         b2   = 0.413_WP/24.0_WP*Re/fVF**2 * &
                (1.0_WP/fVF+3.0_WP*fVF*pVF+8.4_WP*Re**(-0.343_WP)) / &
                (1.0_WP+10.0_WP**(3.0_WP*pVF)*Re**(2.0_WP*fVF-2.5_WP))
         corr = b1 + b2
      case default
         corr = 1.0_WP
      end select

      tau    = this%rho*p%d**2/(18.0_WP*fvisc*corr)
      acc    = (fvel - p%vel)/tau
      opt_dt = tau/real(this%nstep, WP)

   end subroutine get_rhs

   ! -----------------------------------------------------------------------
   !> Compute particle volume fraction field by depositing Pi/6*d^3.
   !> Divides by cell volume to convert to volume fraction.
   !> The user is responsible for calling FillBoundary/syncsum on VF_mf
   !> after this call.
   ! -----------------------------------------------------------------------
   subroutine update_VF(this, VF_mf, lev, geom)
      use amrex_amr_module, only: amrex_multifab, amrex_geometry, amrex_mfiter, &
                                   amrex_mfiter_build, amrex_mfiter_destroy
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_multifab), intent(inout) :: VF_mf
      integer, intent(in) :: lev
      type(amrex_geometry), intent(in) :: geom

      type(amrex_mfiter) :: mfi
      type(c_ptr) :: dp
      type(part), pointer :: p(:)
      integer(c_int64_t) :: np_tile
      integer :: i
      real(WP) :: plo(3), dx(3)
      real(WP), dimension(:,:,:,:), contiguous, pointer :: VFarr

      plo = geom%get_physical_location([geom%domain%lo(0),geom%domain%lo(1),geom%domain%lo(2)])
      dx  = geom%dx

      call VF_mf%setval(0.0_WP)

      call amrex_mfiter_build(mfi, VF_mf, tiling=.false.)
      do while (mfi%valid())
         VFarr => VF_mf%dataptr(mfi)
         call amrlpt_get_particles_mfi(this%pc, lev, mfi%p, dp, np_tile)
         if (np_tile > 0) then
            call c_f_pointer(dp, p, [np_tile])
            do i = 1, int(np_tile)
               if (p(i)%flag == 1) cycle
               call deposit_scalar(Pi/6.0_WP*p(i)%d**3, p(i)%pos, plo, dx, VFarr)
            end do
         end if
         call mfi%next()
      end do
      call amrex_mfiter_destroy(mfi)

      ! Divide by cell volume to get volume fraction
      call VF_mf%mult(1.0_WP/(dx(1)*dx(2)*dx(3)), 0, 1, 0)

   end subroutine update_VF

   ! -----------------------------------------------------------------------
   !> CFL based on particle velocities relative to their local cell size
   ! -----------------------------------------------------------------------
   subroutine get_cfl(this, dt, cflc, cfl)
      use mpi_f08,   only: MPI_ALLREDUCE, MPI_MAX
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional    :: cfl
      ! NOTE: without per-particle dx we can only report |vel|*dt.
      ! A proper implementation queries this%amr%dx(lev) at the particle cell.
      this%CFLp_x = 0.0_WP; this%CFLp_y = 0.0_WP; this%CFLp_z = 0.0_WP
      this%CFL_col = 0.0_WP
      cflc = 0.0_WP
      if (present(cfl)) cfl = 0.0_WP
      ! TODO: loop over tiles, query geom%dx at particle position
   end subroutine get_cfl

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Redistribute particles
   subroutine redistribute(this,minlvl,maxlvl,nover)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in), optional :: minlvl,maxlvl,nover
      integer :: lmin,lmax,no
      lmin= 0; if (present(minlvl)) lmin=minlvl
      lmax=-1; if (present(maxlvl)) lmax=maxlvl
      no  = 0; if (present(nover))  no  =nover
      call amrlpt_redistribute(this%pc,lmin,lmax,no)
      call amrlpt_total_np(this%pc,this%np)
   end subroutine redistribute

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Compute particle statistics: np, d/vel min/max/mean/var.
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(c_ptr) :: dp
      type(part), pointer :: p(:)
      integer(c_int64_t) :: np_tile
      integer :: lev,n,ierr
      real(WP) :: d_sum,d_sq,vx_sum,vx_sq,vy_sum,vy_sq,vz_sum,vz_sq
      real(WP) :: inv_np
      ! Init per-rank accumulators
      this%np=0
      this%dmin=huge(1.0_WP); this%dmax=-huge(1.0_WP);  d_sum=0.0_WP;  d_sq=0.0_WP
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); vx_sum=0.0_WP; vx_sq=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); vy_sum=0.0_WP; vy_sq=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); vz_sum=0.0_WP; vz_sq=0.0_WP
      ! Loop over all AMR levels and tiles
      do lev=0,this%amr%clvl()
         call this%amr%mfiter_build(lev,mfi)
         do while (mfi%valid())
            call amrlpt_get_particles_mfi(this%pc,lev,mfi%p,dp,np_tile)
            if (np_tile.gt.0) then
               call c_f_pointer(dp,p,[np_tile])
               do n=1,int(np_tile)
                  if (p(n)%flag.eq.1) cycle
                  this%np=this%np+1
                  this%dmin=min(this%dmin,p(n)%d);      this%dmax=max(this%dmax,p(n)%d);       d_sum= d_sum+p(n)%d;       d_sq= d_sq+p(n)%d**2
                  this%Umin=min(this%Umin,p(n)%vel(1)); this%Umax=max(this%Umax,p(n)%vel(1)); vx_sum=vx_sum+p(n)%vel(1); vx_sq=vx_sq+p(n)%vel(1)**2
                  this%Vmin=min(this%Vmin,p(n)%vel(2)); this%Vmax=max(this%Vmax,p(n)%vel(2)); vy_sum=vy_sum+p(n)%vel(2); vy_sq=vy_sq+p(n)%vel(2)**2
                  this%Wmin=min(this%Wmin,p(n)%vel(3)); this%Wmax=max(this%Wmax,p(n)%vel(3)); vz_sum=vz_sum+p(n)%vel(3); vz_sq=vz_sq+p(n)%vel(3)**2
               end do
            end if
            call mfi%next()
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Global MPI reduce
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np,  1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sum,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sq,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sum,   1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sq,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sum,   1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sq,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmin,1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax,1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sum,   1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sq,    1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      ! Derive mean and variance
      if (this%np.gt.0) then
         inv_np=1.0_WP/real(this%np,WP)
         this%dmean= d_sum*inv_np; this%dvar=max(0.0_WP, d_sq*inv_np-this%dmean**2)
         this%Umean=vx_sum*inv_np; this%Uvar=max(0.0_WP,vx_sq*inv_np-this%Umean**2)
         this%Vmean=vy_sum*inv_np; this%Vvar=max(0.0_WP,vy_sq*inv_np-this%Vmean**2)
         this%Wmean=vz_sum*inv_np; this%Wvar=max(0.0_WP,vz_sq*inv_np-this%Wmean**2)
      end if
   end subroutine get_info

   !> Print solver info
   subroutine print(this)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(amrlpt), intent(in) :: this
      character(len=str_long) :: message
      call log('AMR Lagrangian particle solver: '//trim(this%name))
      write(message,'("  Particle density  : ",ES12.5)') this%rho
      call log(trim(message))
      write(message,'("  Drag model        : ",a)') trim(this%drag_model)
      call log(trim(message))
      write(message,'("  ODE substeps      : ",i0)') this%nstep
      call log(trim(message))
      call log('  Grid: '//trim(this%amr%name))
   end subroutine print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Write AMReX checkpoint for particles to dirname
   subroutine write(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_write(this%pc,trim(dirname)//c_null_char,1_c_int)
   end subroutine write

   !> Read AMReX checkpoint for particles from dirname
   subroutine read(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_read(this%pc,trim(dirname)//c_null_char)
      call amrlpt_redistribute(this%pc,0,-1,0)
      call amrlpt_total_np(this%pc,this%np)
   end subroutine read

   ! ========================================================================
   ! Interpolation / deposition helpers
   ! Verified against AMReX linear_interpolate_to_particle
   ! (AMReX_TracerParticle_mod_K.H). Formula:
   !   l  = (pos-plo)/dx - 0.5*(1 - is_nodal)
   !   i0 = floor(l),  w = {1-frac, frac},  frac = l - i0
   ! ========================================================================

   !> CIC interpolation from cell-centred scalar (is_nodal=0 in all dirs).
   pure function interp_cc(pos, plo, dx, arr) result(val)
      real(WP), intent(in) :: pos(3), plo(3), dx(3)
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: arr
      real(WP) :: val
      integer  :: i0, j0, k0
      real(WP) :: lx, ly, lz, wx, wy, wz
      lx = (pos(1)-plo(1))/dx(1) - 0.5_WP
      ly = (pos(2)-plo(2))/dx(2) - 0.5_WP
      lz = (pos(3)-plo(3))/dx(3) - 0.5_WP
      i0 = int(floor(lx)); wx = max(0.0_WP, min(1.0_WP, lx - i0))
      j0 = int(floor(ly)); wy = max(0.0_WP, min(1.0_WP, ly - j0))
      k0 = int(floor(lz)); wz = max(0.0_WP, min(1.0_WP, lz - k0))
      val = (1-wx)*(1-wy)*(1-wz)*arr(i0,  j0,  k0,  1) &
           +   wx *(1-wy)*(1-wz)*arr(i0+1,j0,  k0,  1) &
           + (1-wx)*  wy *(1-wz)*arr(i0,  j0+1,k0,  1) &
           +   wx *   wy *(1-wz)*arr(i0+1,j0+1,k0,  1) &
           + (1-wx)*(1-wy)*  wz *arr(i0,  j0,  k0+1,1) &
           +   wx *(1-wy)*  wz *arr(i0+1,j0,  k0+1,1) &
           + (1-wx)*  wy *  wz *arr(i0,  j0+1,k0+1,1) &
           +   wx *   wy *  wz *arr(i0+1,j0+1,k0+1,1)
   end function interp_cc

   !> Trilinear interpolation of a staggered MAC velocity.
   !> Each component is nodal in its own direction, cell-centred in the others.
   pure subroutine interp_mac(pos, plo, dx, Uarr, Varr, Warr, fvel)
      real(WP), intent(in) :: pos(3), plo(3), dx(3)
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: Uarr, Varr, Warr
      real(WP), intent(out) :: fvel(3)
      integer  :: iu, ju, ku, iv, jv, kv, iw, jw, kw
      real(WP) :: lx, ly, lz, wx, wy, wz

      ! U (x-nodal): index [iu, iu+1] in x; cell-centred in y,z
      lx = (pos(1)-plo(1))/dx(1)
      ly = (pos(2)-plo(2))/dx(2) - 0.5_WP
      lz = (pos(3)-plo(3))/dx(3) - 0.5_WP
      iu = int(floor(lx)); wx = max(0.0_WP, min(1.0_WP, lx - iu))
      ju = int(floor(ly)); wy = max(0.0_WP, min(1.0_WP, ly - ju))
      ku = int(floor(lz)); wz = max(0.0_WP, min(1.0_WP, lz - ku))
      fvel(1) = (1-wx)*(1-wy)*(1-wz)*Uarr(iu,  ju,  ku,  1) &
               +   wx *(1-wy)*(1-wz)*Uarr(iu+1,ju,  ku,  1) &
               + (1-wx)*  wy *(1-wz)*Uarr(iu,  ju+1,ku,  1) &
               +   wx *   wy *(1-wz)*Uarr(iu+1,ju+1,ku,  1) &
               + (1-wx)*(1-wy)*  wz *Uarr(iu,  ju,  ku+1,1) &
               +   wx *(1-wy)*  wz *Uarr(iu+1,ju,  ku+1,1) &
               + (1-wx)*  wy *  wz *Uarr(iu,  ju+1,ku+1,1) &
               +   wx *   wy *  wz *Uarr(iu+1,ju+1,ku+1,1)

      ! V (y-nodal)
      lx = (pos(1)-plo(1))/dx(1) - 0.5_WP
      ly = (pos(2)-plo(2))/dx(2)
      lz = (pos(3)-plo(3))/dx(3) - 0.5_WP
      iv = int(floor(lx)); wx = max(0.0_WP, min(1.0_WP, lx - iv))
      jv = int(floor(ly)); wy = max(0.0_WP, min(1.0_WP, ly - jv))
      kv = int(floor(lz)); wz = max(0.0_WP, min(1.0_WP, lz - kv))
      fvel(2) = (1-wx)*(1-wy)*(1-wz)*Varr(iv,  jv,  kv,  1) &
               +   wx *(1-wy)*(1-wz)*Varr(iv+1,jv,  kv,  1) &
               + (1-wx)*  wy *(1-wz)*Varr(iv,  jv+1,kv,  1) &
               +   wx *   wy *(1-wz)*Varr(iv+1,jv+1,kv,  1) &
               + (1-wx)*(1-wy)*  wz *Varr(iv,  jv,  kv+1,1) &
               +   wx *(1-wy)*  wz *Varr(iv+1,jv,  kv+1,1) &
               + (1-wx)*  wy *  wz *Varr(iv,  jv+1,kv+1,1) &
               +   wx *   wy *  wz *Varr(iv+1,jv+1,kv+1,1)

      ! W (z-nodal)
      lx = (pos(1)-plo(1))/dx(1) - 0.5_WP
      ly = (pos(2)-plo(2))/dx(2) - 0.5_WP
      lz = (pos(3)-plo(3))/dx(3)
      iw = int(floor(lx)); wx = max(0.0_WP, min(1.0_WP, lx - iw))
      jw = int(floor(ly)); wy = max(0.0_WP, min(1.0_WP, ly - jw))
      kw = int(floor(lz)); wz = max(0.0_WP, min(1.0_WP, lz - kw))
      fvel(3) = (1-wx)*(1-wy)*(1-wz)*Warr(iw,  jw,  kw,  1) &
               +   wx *(1-wy)*(1-wz)*Warr(iw+1,jw,  kw,  1) &
               + (1-wx)*  wy *(1-wz)*Warr(iw,  jw+1,kw,  1) &
               +   wx *   wy *(1-wz)*Warr(iw+1,jw+1,kw,  1) &
               + (1-wx)*(1-wy)*  wz *Warr(iw,  jw,  kw+1,1) &
               +   wx *(1-wy)*  wz *Warr(iw+1,jw,  kw+1,1) &
               + (1-wx)*  wy *  wz *Warr(iw,  jw+1,kw+1,1) &
               +   wx *   wy *  wz *Warr(iw+1,jw+1,kw+1,1)
   end subroutine interp_mac

   !> NGP deposit of scalar Sp to nearest cell in arr
   subroutine deposit_scalar(Sp, pos, plo, dx, arr)
      real(WP), intent(in) :: Sp, pos(3), plo(3), dx(3)
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(inout) :: arr
      integer :: i, j, k
      i = int(floor((pos(1)-plo(1))/dx(1)))
      j = int(floor((pos(2)-plo(2))/dx(2)))
      k = int(floor((pos(3)-plo(3))/dx(3)))
      arr(i,j,k,1) = arr(i,j,k,1) + Sp
   end subroutine deposit_scalar

end module amrlpt_class
