!> Test amrscalar solver with time integration (all-level API)
module mod_test_amrscalar
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrscalar_class,   only: amrscalar
   use amrdata_class,     only: amrdata, amrex_interp_reinit, amrex_interp_none
   use amrex_amr_module,  only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
   implicit none
   private
   public :: test_amrscalar

   ! Module-level objects (referenced by callbacks)
   type(amrgrid),   allocatable, target :: amr
   type(amrscalar), allocatable, target :: sc

   ! Velocity fields (amrdata, face-centered)
   type(amrdata), allocatable, target :: U, V, W

   ! dSCdt storage (amrdata, cell-centered)
   type(amrdata), allocatable, target :: dSCdt

   real(WP), parameter :: SC_REFINE_THRESH=0.01_WP  !< Refine where SC > this value

contains

   !> Custom on_init callback: initialize scalar field with Gaussian blob
   subroutine gaussian_init(solver, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_mfiter_build, amrex_mfiter_destroy
      class(amrscalar), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC
      real(WP) :: x, y, z, dx, dy, dz
      integer :: i, j, k
      ! Initialize with Gaussian blob
      dx = solver%amr%dx(lvl)
      dy = solver%amr%dy(lvl)
      dz = solver%amr%dz(lvl)
      call amrex_mfiter_build(mfi, solver%SC%mf(lvl))
      do while (mfi%next())
         bx = mfi%tilebox()
         pSC => solver%SC%mf(lvl)%dataptr(mfi)
         do k = bx%lo(3), bx%hi(3)
            z = solver%amr%zlo + (real(k,WP)+0.5_WP)*dz
            do j = bx%lo(2), bx%hi(2)
               y = solver%amr%ylo + (real(j,WP)+0.5_WP)*dy
               do i = bx%lo(1), bx%hi(1)
                  x = solver%amr%xlo + (real(i,WP)+0.5_WP)*dx
                  ! Gaussian offset from center
                  pSC(i,j,k,1) = exp(-200.0_WP*((x-0.25_WP)**2 + (y-0.5_WP)**2 + (z-0.5_WP)**2))
               end do
            end do
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine gaussian_init


   !> Custom tagging callback: refine where SC > threshold
   subroutine scalar_tagger(solver, lvl, time, tags_ptr)
      use iso_c_binding,    only: c_ptr, c_char
      use amrex_amr_module, only: amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrscalar), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSC
      integer :: i, j, k
      tags = tags_ptr
      ! Iterate over level boxes
      call solver%amr%mfiter_build(lvl, mfi)
      do while (mfi%next())
         bx = mfi%tilebox()
         pSC => solver%SC%mf(lvl)%dataptr(mfi)
         tagarr => tags%dataPtr(mfi)
         do k = bx%lo(3), bx%hi(3)
            do j = bx%lo(2), bx%hi(2)
               do i = bx%lo(1), bx%hi(1)
                  if (pSC(i,j,k,1) .gt. SC_REFINE_THRESH) then
                     tagarr(i,j,k,1) = SETtag
                  end if
               end do
            end do
         end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine scalar_tagger


   !> Velocity initialization callback (called by amrdata on init/regrid)
   !> Uses this%name to determine which component (U, V, or W)
   subroutine velocity_init(this, lvl, time, ba, dm)
      use amrex_amr_module, only: amrex_mfiter, amrex_box, amrex_boxarray, amrex_distromap, &
      &                           amrex_mfiter_build, amrex_mfiter_destroy
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: x, y, dx, dy
      integer :: i, j, k
      dx = amr%dx(lvl)
      dy = amr%dy(lvl)
      ! Initialize velocity component based on name
      call amrex_mfiter_build(mfi, this%mf(lvl), tiling=this%amr%default_tiling)
      do while (mfi%next())
         bx = mfi%tilebox()
         p => this%mf(lvl)%dataptr(mfi)
         select case (trim(this%name))
          case ('U')
            ! U on x-faces: U = -(y-0.5)
            do k = lbound(p,3), ubound(p,3)
               do j = lbound(p,2), ubound(p,2)
                  y = amr%ylo + (real(j,WP)+0.5_WP)*dy
                  do i = lbound(p,1), ubound(p,1)
                     p(i,j,k,1) = -(y-0.5_WP)
                  end do
               end do
            end do
          case ('V')
            ! V on y-faces: V = +(x-0.5)
            do k = lbound(p,3), ubound(p,3)
               do j = lbound(p,2), ubound(p,2)
                  do i = lbound(p,1), ubound(p,1)
                     x = amr%xlo + (real(i,WP)+0.5_WP)*dx
                     p(i,j,k,1) = +(x-0.5_WP)
                  end do
               end do
            end do
          case ('W')
            ! W on z-faces: W = 0
            p = 0.0_WP
         end select
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine velocity_init


   subroutine test_amrscalar()
      use string,            only: itoa,rtoa,str_medium
      use mathtools,         only: Pi
      use amrviz_class,      only: amrviz
      use amrio_class,       only: amrio
      use amrdata_class,     only: amrex_interp_none
      use monitor_class,     only: monitor
      use timetracker_class, only: timetracker
      use event_class,       only: event
      use messager,          only: log,warn
      implicit none
      type(amrviz) :: viz
      type(amrio) :: io
      type(monitor) :: mfile
      type(timetracker) :: time
      type(event) :: regrid_evt,hdf5_evt
      real(WP) :: int0,intF

      call log("--------------------------------------------------------")
      call log("Running Test: amrscalar Time Integration (All-Level API)")
      call log("--------------------------------------------------------")

      ! Allocate module objects
      allocate(amr,sc,U,V,W,dSCdt)

      ! Setup grid (periodic box)
      amr%nx  =32
      amr%ny  =32
      amr%nz  =32
      amr%xlo =0.0_WP
      amr%xhi =1.0_WP
      amr%ylo =0.0_WP
      amr%yhi =1.0_WP
      amr%zlo =0.0_WP
      amr%zhi =1.0_WP
      amr%xper=.true.
      amr%yper=.true.
      amr%zper=.true.
      amr%maxlvl=2   ! 3 levels (0,1,2)
      amr%nmax=16

      call amr%initialize("sc_amr")

      ! Initialize scalar solver (default: use_refluxing=.false. for flux averaging)
      call sc%initialize(amr, nscalar=1, name="test_scalar")
      sc%user_init => gaussian_init
      sc%user_tagging => scalar_tagger

      ! Initialize velocity fields with reinit mode (auto-recomputed on regrid)
      call U%initialize(amr, name='U', ncomp=1, ng=0, nodal=[.true., .false., .false.], interp=amrex_interp_reinit)
      call V%initialize(amr, name='V', ncomp=1, ng=0, nodal=[.false., .true., .false.], interp=amrex_interp_reinit)
      call W%initialize(amr, name='W', ncomp=1, ng=0, nodal=[.false., .false., .true.], interp=amrex_interp_reinit)
      U%user_init => velocity_init
      V%user_init => velocity_init
      W%user_init => velocity_init
      call U%register()
      call V%register()
      call W%register()

      ! Initialize dSCdt as workspace (no interpolation needed)
      call dSCdt%initialize(amr, name='dSCdt', ncomp=sc%nscalar, ng=0, interp=amrex_interp_none)
      call dSCdt%register()

      ! Initialize HDF5 viz output
      call viz%initialize(amr=amr, name='sc_advect')
      call viz%add_scalar(data=sc%SC, comp=1, name='SC')

      ! Build all levels (callbacks auto-allocate U, V, W, dSCdt)
      call amr%init_from_scratch(time=0.0_WP, do_postregrid=.true.)
      call amr%get_info()
      call log("After init_from_scratch: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")

      ! Get initial integral
      call sc%get_info()
      int0=sc%SCint(1)
      call log("Initial: SCint="//trim(rtoa(int0)))

      ! Setup timetracker
      time=timetracker(amRoot=amr%amRoot)
      time%dt=0.0025_WP
      time%dtmax=time%dt
      time%tmax=2.0_WP*Pi  ! Full circle rotation

      ! Setup regrid event (every 10 steps)
      regrid_evt=event(time=time,name='Regrid')
      regrid_evt%nper=10

      ! Setup HDF5 event
      hdf5_evt=event(time=time,name='HDF5')
      hdf5_evt%tper=0.125_WP

      ! Setup monitor file
      call sc%get_info()
      mfile=monitor(amr%amRoot,'scalar')
      call mfile%add_column(time%n,'Timestep')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(sc%SCmin(1),'SC_min')
      call mfile%add_column(sc%SCmax(1),'SC_max')
      call mfile%add_column(sc%SCint(1),'SC_int')
      call mfile%write()

      ! Write initial output
      call viz%write(time=time%t)

      call log("Advancing to t="//trim(rtoa(time%tmax))//", dt="//trim(rtoa(time%dt)))

      ! Time loop
      do while (.not.time%done())
         call time%increment()

         ! Copy to old: SCold=SC
         call sc%SCold%copy(src=sc%SC)

         ! Calculate dSC/dt for all levels
         call sc%get_dSCdt(U=U, V=V, W=W, SC=sc%SCold, dSCdt=dSCdt)

         ! Forward Euler step: SC=1.0_WP*SCold+dt*dSCdt
         call sc%SC%lincomb(a=1.0_WP, src1=sc%SCold, b=time%dt, src2=dSCdt)

         ! Reflux (only applies if use_refluxing=.true.) and average down
         call sc%reflux(dt=time%dt)
         call sc%average_down()

         ! Regrid if needed
         if (regrid_evt%occurs()) then
            call log("Regridding at step "//trim(itoa(time%n)))
            call amr%regrid(baselvl=0,time=time%t)
            call log("  Grid: "//trim(itoa(amr%nlevels))//" levels, "//trim(itoa(amr%nboxes))//" boxes")
         end if

         ! Write HDF5 output
         if (hdf5_evt%occurs()) then
            call viz%write(time=time%t)
            call log("Step "//trim(itoa(time%n))//": t="//trim(rtoa(time%t))//", SCint="//trim(rtoa(sc%SCint(1))))
         end if

         ! Update monitor file
         call sc%get_info()
         call mfile%write()
      end do

      ! Final conservation check
      call sc%get_info()
      intF=sc%SCint(1)
      call log("Final: SCint="//trim(rtoa(intF)))
      call log("Conservation error: "//trim(rtoa(abs(intF-int0)/int0*100.0_WP))//"%")

      if (abs(intF-int0)/int0<1.0e-10_WP) then
         call log("PASS: Scalar conserved to machine precision")
      else
         call warn("FAIL: Conservation error too large!")
      end if

      ! Test checkpoint round-trip
      call io%initialize(amr, nfiles=1)
      call sc%register_checkpoint(io)
      call io%add_scalar('dt', time%dt)
      call io%write('checkpoint_test', time%t, time%n)
      call log("PASS: Checkpoint written to checkpoint_test/")

      ! Save current integral, then zero the data
      call sc%get_info()
      int0 = sc%SCint(1)
      call log("  Before zero: SCint="//trim(rtoa(int0)))

      ! Zero out SC data
      call sc%SC%setval(val=0.0_WP)
      call sc%get_info()
      call log("  After zero: SCint="//trim(rtoa(sc%SCint(1))))

      ! Read checkpoint back
      block
         real(WP) :: read_time, read_dt
         integer :: read_step
         call io%read_header('checkpoint_test', read_time, read_step)
         call io%get_scalar('dt', read_dt)
         call sc%restore_checkpoint(io, 'checkpoint_test', read_time)
         call log("  Read checkpoint: time="//trim(rtoa(read_time))//" step="//trim(itoa(read_step))//" dt="//trim(rtoa(read_dt)))
      end block

      ! Verify data restored
      call sc%get_info()
      intF = sc%SCint(1)
      call log("  After read: SCint="//trim(rtoa(intF)))

      if (abs(intF-int0)/int0 .lt. 1.0e-10_WP) then
         call log("PASS: Checkpoint round-trip verified!")
      else
         call warn("FAIL: Checkpoint round-trip failed!")
      end if

      call log("PASS: HDF5 plotfiles written to amrviz/sc_advect/")

      ! Cleanup
      call io%finalize()
      call viz%finalize()
      call mfile%finalize()
      call dSCdt%finalize()
      call W%finalize()
      call V%finalize()
      call U%finalize()
      call sc%finalize()
      call amr%finalize()
      if (allocated(dSCdt)) deallocate(dSCdt)
      if (allocated(W)) deallocate(W)
      if (allocated(V)) deallocate(V)
      if (allocated(U)) deallocate(U)
      if (allocated(sc)) deallocate(sc)
      if (allocated(amr)) deallocate(amr)
      call log("PASS: amrscalar test complete!")

   end subroutine test_amrscalar

end module mod_test_amrscalar
