!> Test amrincomp solver - advection test with time loop
module mod_test_amrincomp
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrincomp_class,   only: amrincomp
   use amrdata_class,     only: amrdata, amrex_interp_none
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use string,            only: itoa, rtoa
   implicit none
   private
   public :: test_amrincomp

   ! Grid
   type(amrgrid), allocatable, target :: amr

   ! Solver data
   type(timetracker) :: time
   type(amrincomp), allocatable, target :: fs
   type(amrdata),   allocatable :: resU, resV, resW

   ! Visualization
   type(amrviz), allocatable :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: tstart_regrid
   real(WP) :: tagging_threshold=huge(1.0_WP)

   ! HIT forcing
   real(WP) :: K_target
   real(WP) :: tau_forcing

   ! Monitoring
   type(monitor) :: mfile, cflfile, consfile, gridfile

contains

   !> User-provided initialization callback for velocity - random field
   subroutine velocity_init(solver, lvl, time, ba, dm)
      use amrex_amr_module,  only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, amrex_mfiter_build, amrex_mfiter_destroy
      use random, only: random_uniform
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      integer :: i, j, k
      ! For fine levels: fill from coarse using div-free interpolation
      if (lvl .gt. 0) then
         call solver%fill_velocity_from_coarse(lvl, time)
         return
      end if
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

   !> Geometric tagger: refine center box
   subroutine geometric_tagger(solver, lvl, time, tags_ptr)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP) :: xc, yc, zc, x, y, z, dx, dy, dz, radius
      integer :: i, j, k

      xc = 0.5_WP * (solver%amr%xlo + solver%amr%xhi)
      yc = 0.5_WP * (solver%amr%ylo + solver%amr%yhi)
      zc = 0.5_WP * (solver%amr%zlo + solver%amr%zhi)
      radius = 0.25_WP * min(solver%amr%xhi - solver%amr%xlo, &
      &                      solver%amr%yhi - solver%amr%ylo, &
      &                      solver%amr%zhi - solver%amr%zlo)

      dx = solver%amr%dx(lvl); dy = solver%amr%dy(lvl); dz = solver%amr%dz(lvl)
      tags = tags_ptr

      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
            x = solver%amr%xlo + (real(i,WP) + 0.5_WP) * dx
            y = solver%amr%ylo + (real(j,WP) + 0.5_WP) * dy
            z = solver%amr%zlo + (real(k,WP) + 0.5_WP) * dz
            if (sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2) .lt. radius) tagarr(i,j,k,1) = SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine geometric_tagger

   !> Vorticity-based tagger: refine where vorticity magnitude exceeds threshold
   subroutine vorticity_tagger(solver, lvl, time, tags_ptr)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP) :: dx, dy, dz
      real(WP) :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
      real(WP) :: omega_x, omega_y, omega_z, vort_mag
      integer :: i, j, k
      dx = solver%amr%dx(lvl)
      dy = solver%amr%dy(lvl)
      dz = solver%amr%dz(lvl)
      tags = tags_ptr
      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         pU => solver%U%mf(lvl)%dataptr(mfi)
         pV => solver%V%mf(lvl)%dataptr(mfi)
         pW => solver%W%mf(lvl)%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  ! du/dy at cell center
                  dudy = 0.25_WP*(pU(i,j+1,k,1)+pU(i+1,j+1,k,1)-pU(i,j-1,k,1)-pU(i+1,j-1,k,1))/(2.0_WP*dy)
                  ! du/dz at cell center
                  dudz = 0.25_WP*(pU(i,j,k+1,1)+pU(i+1,j,k+1,1)-pU(i,j,k-1,1)-pU(i+1,j,k-1,1))/(2.0_WP*dz)
                  ! dv/dx at cell center
                  dvdx = 0.25_WP*(pV(i+1,j,k,1)+pV(i+1,j+1,k,1)-pV(i-1,j,k,1)-pV(i-1,j+1,k,1))/(2.0_WP*dx)
                  ! dv/dz at cell center
                  dvdz = 0.25_WP*(pV(i,j,k+1,1)+pV(i,j+1,k+1,1)-pV(i,j,k-1,1)-pV(i,j+1,k-1,1))/(2.0_WP*dz)
                  ! dw/dx at cell center
                  dwdx = 0.25_WP*(pW(i+1,j,k,1)+pW(i+1,j,k+1,1)-pW(i-1,j,k,1)-pW(i-1,j,k+1,1))/(2.0_WP*dx)
                  ! dw/dy at cell center
                  dwdy = 0.25_WP*(pW(i,j+1,k,1)+pW(i,j+1,k+1,1)-pW(i,j-1,k,1)-pW(i,j-1,k+1,1))/(2.0_WP*dy)
                  ! Vorticity components
                  omega_x = dwdy - dvdz
                  omega_y = dudz - dwdx
                  omega_z = dvdx - dudy
                  ! Vorticity magnitude
                  vort_mag = sqrt(omega_x**2 + omega_y**2 + omega_z**2)
                  ! Tag if above threshold
                  if (vort_mag .gt. tagging_threshold) tagarr(i,j,k,1) = SETtag
               end do
            end do
         end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine vorticity_tagger

   !> Velocity gradient magnitude tagger: refine where |∇u| > threshold
   !> |∇u| = sqrt(sum of all (dui/dxj)²) - captures all gradients equally
   subroutine gradU_tagger(solver, lvl, time, tags_ptr)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrincomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU, pV, pW
      real(WP) :: dx, dy, dz
      real(WP) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, gradU_mag
      integer :: i, j, k
      dx = solver%amr%dx(lvl); dy = solver%amr%dy(lvl); dz = solver%amr%dz(lvl)
      tags = tags_ptr
      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         tagarr => tags%dataPtr(mfi)
         pU => solver%U%mf(lvl)%dataptr(mfi)
         pV => solver%V%mf(lvl)%dataptr(mfi)
         pW => solver%W%mf(lvl)%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3); do j = bx%lo(2), bx%hi(2); do i = bx%lo(1), bx%hi(1)
            ! Diagonal gradients
            dudx = (pU(i+1,j,k,1) - pU(i,j,k,1)) / dx
            dvdy = (pV(i,j+1,k,1) - pV(i,j,k,1)) / dy
            dwdz = (pW(i,j,k+1,1) - pW(i,j,k,1)) / dz
            ! Off-diagonal gradients
            dudy = 0.25_WP*(pU(i,j+1,k,1)+pU(i+1,j+1,k,1)-pU(i,j-1,k,1)-pU(i+1,j-1,k,1))/(2.0_WP*dy)
            dudz = 0.25_WP*(pU(i,j,k+1,1)+pU(i+1,j,k+1,1)-pU(i,j,k-1,1)-pU(i+1,j,k-1,1))/(2.0_WP*dz)
            dvdx = 0.25_WP*(pV(i+1,j,k,1)+pV(i+1,j+1,k,1)-pV(i-1,j,k,1)-pV(i-1,j+1,k,1))/(2.0_WP*dx)
            dvdz = 0.25_WP*(pV(i,j,k+1,1)+pV(i,j+1,k+1,1)-pV(i,j,k-1,1)-pV(i,j+1,k-1,1))/(2.0_WP*dz)
            dwdx = 0.25_WP*(pW(i+1,j,k,1)+pW(i+1,j,k+1,1)-pW(i-1,j,k,1)-pW(i-1,j,k+1,1))/(2.0_WP*dx)
            dwdy = 0.25_WP*(pW(i,j+1,k,1)+pW(i,j+1,k+1,1)-pW(i,j-1,k,1)-pW(i,j-1,k+1,1))/(2.0_WP*dy)
            ! |∇u| = sqrt(sum of all gradients squared)
            gradU_mag = sqrt(dudx**2 + dudy**2 + dudz**2 + &
            &                dvdx**2 + dvdy**2 + dvdz**2 + &
            &                dwdx**2 + dwdy**2 + dwdz**2)
            if (gradU_mag .gt. tagging_threshold) tagarr(i,j,k,1) = SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine gradU_tagger

   !> Main test routine
   subroutine test_amrincomp()
      implicit none

      ! Create amrgrid
      create_amrgrid: block
         use param, only: param_read
         allocate(amr)
         amr%name = 'advect_test'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo = 0.0_WP; amr%xhi = 1.0_WP
         amr%ylo = 0.0_WP; amr%yhi = 1.0_WP
         amr%zlo = 0.0_WP; amr%zhi = 1.0_WP
         amr%xper = .true.; amr%yper = .true.; amr%zper = .true.
         call param_read('Max level',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         use param, only: param_read
         time = timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt = time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Create flow solver
      create_flow_solver: block
         use param, only: param_read
         allocate(fs)
         call param_read('Density',fs%rho)
         call param_read('Viscosity',fs%visc)
         fs%psolver%max_iter = 20
         fs%psolver%tol_rel = 1.0e-6_WP
         fs%user_init => velocity_init
         fs%user_tagging => gradU_tagger
         call fs%initialize(amr, name='advect_fs')
      end block create_flow_solver

      ! Create workspace arrays
      create_workspace: block
         allocate(resU); call resU%initialize(amr, name='resU', ncomp=1, ng=0, nodal=[.true., .false., .false.], interp=amrex_interp_none); call resU%register()
         allocate(resV); call resV%initialize(amr, name='resV', ncomp=1, ng=0, nodal=[.false., .true., .false.], interp=amrex_interp_none); call resV%register()
         allocate(resW); call resW%initialize(amr, name='resW', ncomp=1, ng=0, nodal=[.false., .false., .true.], interp=amrex_interp_none); call resW%register()
      end block create_workspace

      ! Initialize grid
      initialize: block
         call amr%init_from_scratch(time=time%t)
      end block initialize

      ! Initialize HIT forcing
      initialize_forcing: block
         use param, only: param_read
         call param_read('Target TKE', K_target)
         call param_read('Forcing timescale', tau_forcing)
      end block initialize_forcing

      ! Initialize visualization
      create_visualization: block
         use param, only: param_read
         ! Create visualization object
         allocate(viz); call viz%initialize(amr, 'test_advect')
         call viz%add_scalar(fs%P, 1, 'pressure')
         call viz%add_scalar(fs%div, 1, 'divergence')
         call viz%add_scalar(fs%U, 1, 'U')
         call viz%add_scalar(fs%V, 1, 'V')
         call viz%add_scalar(fs%W, 1, 'W')
         ! Create visualization output event
         viz_evt = event(time=time, name='Visualization output')
         call param_read('Output period', viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Regridding parameters
      regrid_setup: block
         use param, only: param_read
         ! Create regridding event
         regrid_evt = event(time=time, name='Regrid')
         call param_read('Regrid nsteps', regrid_evt%nper)
         ! Set regridding start time
         call param_read('Regrid start', tstart_regrid)
         ! Set tagging threshold
         call param_read('Tagging threshold', tagging_threshold)
      end block regrid_setup

      ! Create monitor
      create_monitor: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(time%dt,time%cfl)
         ! Create simulation monitor
         mfile = monitor(amRoot=amr%amRoot, name='simulation')
         call mfile%add_column(time%n, 'Timestep')
         call mfile%add_column(time%t, 'Time')
         call mfile%add_column(time%dt, 'dt')
         call mfile%add_column(fs%CFL, 'CFL')
         call mfile%add_column(fs%Umax, 'Umax')
         call mfile%add_column(fs%Vmax, 'Vmax')
         call mfile%add_column(fs%Wmax, 'Wmax')
         call mfile%add_column(fs%Pmax, 'Pmax')
         call mfile%add_column(fs%divmax, 'Divergence')
         call mfile%write()
         ! Create CFL monitor
         cflfile = monitor(amRoot=amr%amRoot, name='cfl')
         call cflfile%add_column(time%n, 'Timestep')
         call cflfile%add_column(time%t, 'Time')
         call cflfile%add_column(time%dt, 'dt')
         call cflfile%add_column(fs%CFLc_x, 'CFLc_x')
         call cflfile%add_column(fs%CFLc_y, 'CFLc_y')
         call cflfile%add_column(fs%CFLc_z, 'CFLc_z')
         call cflfile%add_column(fs%CFLv_x, 'CFLv_x')
         call cflfile%add_column(fs%CFLv_y, 'CFLv_y')
         call cflfile%add_column(fs%CFLv_z, 'CFLv_z')
         call cflfile%write()
         ! Create conservation monitor
         consfile = monitor(amRoot=amr%amRoot, name='conservation')
         call consfile%add_column(time%n, 'Timestep')
         call consfile%add_column(time%t, 'Time')
         call consfile%add_column(fs%rhoUint, 'rhoUint')
         call consfile%add_column(fs%rhoVint, 'rhoVint')
         call consfile%add_column(fs%rhoWint, 'rhoWint')
         call consfile%add_column(fs%rhoKint, 'rhoKint')
         call consfile%write()
         ! Create grid monitor
         gridfile = monitor(amRoot=amr%amRoot, name='grid')
         call gridfile%add_column(time%n, 'Timestep')
         call gridfile%add_column(time%t, 'Time')
         call gridfile%add_column(amr%nlevels, 'Nlvl')
         call gridfile%add_column(amr%nboxes, 'Nbox')
         call gridfile%add_column(amr%ncells, 'Ncell')
         call gridfile%add_column(amr%compression, 'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitor

      ! Time integration loop
      time_loop: do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old velocity
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)

         ! Sub-iterations
         do while (time%it .le. time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%U%lincomb(a=0.5_WP, src1=fs%U, b=0.5_WP, src2=fs%Uold)
            call fs%V%lincomb(a=0.5_WP, src1=fs%V, b=0.5_WP, src2=fs%Vold)
            call fs%W%lincomb(a=0.5_WP, src1=fs%W, b=0.5_WP, src2=fs%Wold)

            ! Compute advective momentum RHS
            call fs%get_dmomdt(U=fs%U, V=fs%V, W=fs%W, drhoUdt=resU, drhoVdt=resV, drhoWdt=resW, time=time%t)

            ! Add TKE-targeting linear forcing: F = A * rho * (U - U_mean)
            ! Forces only fluctuations, ensuring zero net momentum injection
            tke_target_forcing: block
               real(WP) :: A_forcing
               if (fs%rhoKint .gt. 1.0e-10_WP) then
                  A_forcing = (K_target - fs%rhoKint) / (2.0_WP * fs%rhoKint * tau_forcing)
               else
                  A_forcing = 1.0_WP       ! Strong forcing if TKE is near zero
               end if
               ! Force fluctuations: F = A * rho * (U - Umean)
               call resU%saxpy(a=A_forcing*fs%rho, src=fs%U); call resU%plus(val=-A_forcing*fs%rhoUint/fs%amr%vol)
               call resV%saxpy(a=A_forcing*fs%rho, src=fs%V); call resV%plus(val=-A_forcing*fs%rhoVint/fs%amr%vol)
               call resW%saxpy(a=A_forcing*fs%rho, src=fs%W); call resW%plus(val=-A_forcing*fs%rhoWint/fs%amr%vol)
            end block tke_target_forcing

            ! Increment velocity
            call fs%U%lincomb(a=1.0_WP, src1=fs%Uold, b=time%dt/fs%rho, src2=resU)
            call fs%V%lincomb(a=1.0_WP, src1=fs%Vold, b=time%dt/fs%rho, src2=resV)
            call fs%W%lincomb(a=1.0_WP, src1=fs%Wold, b=time%dt/fs%rho, src2=resW)

            ! Average down and fill ghosts
            call fs%average_down_velocity()
            call fs%fill_velocity(time%t)

            ! Compute divergence
            call fs%get_div()

            ! Solve pressure Poisson
            call fs%div%mult(val=fs%rho/time%dt)
            call fs%psolver%solve(rhs=fs%div)

            ! Correct velocity
            call fs%add_pressure(scale=time%dt/fs%rho)

            ! Average down and fill ghosts
            call fs%average_down_velocity()
            call fs%fill_velocity(time%t)

            ! Increment sub-iteration counter
            time%it = time%it + 1

         end do

         ! Regrid if event triggers
         if (regrid_evt%occurs().and. time%t .gt. tstart_regrid) then
            call amr%regrid(baselvl=0, time=time%t)
            call gridfile%write()
         end if

         ! Monitor output
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call consfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         
      end do time_loop

      ! Cleanup
      cleanup: block
         call viz%finalize()
         call fs%finalize()
         call amr%finalize()
         call mfile%finalize()
         call cflfile%finalize()
         call consfile%finalize()
         call gridfile%finalize()
         call time%finalize()
         call viz_evt%finalize()
         call regrid_evt%finalize()
         call resU%finalize()
         call resV%finalize()
         call resW%finalize()
         deallocate(viz, fs, amr, resU, resV, resW)
      end block cleanup

      call log("Advection test complete")

   end subroutine test_amrincomp

end module mod_test_amrincomp
