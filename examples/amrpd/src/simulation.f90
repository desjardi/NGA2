!> AMRPD development driver
!>
!> Long-running exploratory test for the amrpd peridynamics solver. Starts as a
!> minimal "create a cube of particles and confirm they land in the container"
!> driver; grows feature-by-feature as the solver matures (motion, BCs, bonds,
!> visualization, restart, fluid coupling).
module simulation
   use precision,         only: WP,I8
   use amrgrid_class,     only: amrgrid
   use amrpd_class,       only: amrpd,part,PART_MOVES,PART_INTEGRATES,PART_BONDS,AMRPD_WALL,AMRPD_OPEN
   use amrpdviz_class,    only: amrpdviz
   use amrviz_class,      only: amrviz
   use amrio_class,       only: amrio
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use string,            only: str_medium,str_long
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Peridynamics solver under test
   type(amrpd), target :: pd

   !> Particle visualization
   type(amrpdviz) :: pviz

   !> Eulerian-field visualization (VF and any other mesh fields we want to view)
   type(amrviz) :: viz

   !> Time tracker
   type(timetracker) :: time

   !> Regrid event (step-based; fires every Regrid nsteps)
   type(event) :: regrid_evt

   !> Visualization output event
   type(event) :: viz_evt

   !> Checkpoint I/O and event
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir=''
   logical  :: restarted=.false.
   real(WP) :: restart_time=0.0_WP

   !> Monitors (written each step)
   type(monitor) :: mfile          !< Overall simulation (time, dt, CFL, np, velocity stats)
   type(monitor) :: cflfile        !< Per-direction CFL breakdown
   type(monitor) :: gridfile       !< AMR grid stats
   type(monitor) :: balancefile    !< Load-balance metrics (np_min/max/eff)

   !> Cube initialization parameters
   integer  :: cube_n=0           !< Particles per side (derived from cube_L / elem_size)
   real(WP) :: cube_L=0.5_WP      !< Side length of the cube
   real(WP) :: elem_size=0.0_WP   !< Solid element (particle) spacing dp

contains

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none

      ! Create amrgrid
      create_amrgrid: block
         amr%name='amrpd'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-4.0_WP; amr%xhi=+4.0_WP
         amr%ylo=-4.0_WP; amr%yhi=+4.0_WP
         amr%zlo=-1.0_WP; amr%zhi=+7.0_WP
         amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
         call param_read('Max level',amr%maxlvl)
         ! Handle 2D case
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zper=.true.
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart: initialize IO object and probe for a checkpoint
      handle_restart: block
         integer :: restart_step
         call io%initialize(amr=amr,nfiles=1)
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker (use restart_time + restored dt if restarting)
      initialize_time: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',  time%dtmax)
         call param_read('Max CFL', time%cflmax,default=0.5_WP)
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_time

      ! Initialize amrpd solver (creates particle + bond containers, registers callbacks)
      init_amrpd: block
         call pd%initialize(amr,name='pd')
         ! Material parameters
         call param_read('Material density',  pd%rho)
         call param_read('Elastic modulus',   pd%elastic_modulus)
         call param_read('Poisson ratio',     pd%poisson_ratio)
         call param_read('Critical energy',   pd%crit_energy)
         ! Solid discretization: user gives the particle spacing dp (elem_size).
         ! Horizon defaults to 3.0125 * elem_size (Peridigm convention; small
         ! eps offset avoids floating-point ties with the 3rd neighbor shell).
         call param_read('Element size',      elem_size)
         call param_read('Horizon',           pd%delta, default=3.0125_WP*elem_size)
         pd%search_radius=1.5_WP*pd%delta        !< Ghost-layer radius, used by every fill_ghosts call. Placeholder
                                                 !< 1.5*delta; once damage lands this becomes (1+max_stretch)*delta+safety.
         ! Contact (soft-sphere). Leave Collision time at 0 to let amrpd
         ! auto-set tau_col = 5*dt each step (as stiff as integrable).
         ! Positive value = user override (fixed tau_col).
         call param_read('Collision time',    pd%tau_col,default=0.0_WP)
         call param_read('Restitution part.', pd%e_n,    default=0.7_WP)
         call param_read('Restitution wall',  pd%e_w,    default=0.7_WP)
         ! Particle-driven AMR refinement. Tagging VF<=0 disables tagging.
         ! (filter_width is set below from elem_size, mirroring amrlpt's convention.)
         call param_read('Tagging VF', pd%VF_tag, default=-1.0_WP)
         ! Gravity (default = freefall in -z)
         call param_read('Gravity',pd%gravity,default=[0.0_WP,0.0_WP,-9.81_WP])
         ! pd%maxlvl defaults to amr%maxlvl (particles can occupy [0, amr%maxlvl])
         ! Domain BCs: WALL on all faces so the cube bounces around inside the box.
         ! Override per face in input if desired (open exit, etc).
         pd%lo_bc=AMRPD_WALL
         pd%hi_bc=AMRPD_WALL
         ! Optional knapsack rebalance of the particle+bond containers
         call param_read('Balance particles',pd%rebalance,default=.false.)
         ! Cube initialization parameters: user gives the cube side, particle
         ! count is derived from cube_L / elem_size. Adjust cube_L to be an
         ! integer multiple of elem_size for a clean lattice.
         call param_read('Cube size', cube_L,default=0.5_WP)
         cube_n=nint(cube_L/elem_size)
         cube_L=real(cube_n,WP)*elem_size
         pd%dV =elem_size**3
         ! VF filter width: 7 particle spacings (matches the amrlpt convention
         ! of 7*inj_dmean -- smooths the deposited field over a few neighbors)
         !pd%filter_width=7.0_WP*elem_size
      end block init_amrpd

      ! Build (or restore) the AMR grid and the particle/bond populations
      build_grid_and_state: block
         if (restarted) then
            ! Restore grid hierarchy and both particle/bond containers
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            call pd%read(dirname=trim(restart_dir))
         else
            ! Fresh start: build grid, then seed a Cartesian cube of particles
            call amr%init_from_scratch(time=time%t)
            seed_cube: block
               use mathtools, only: Pi
               type(part), dimension(:), allocatable, target :: plist
               integer(I8) :: ntot,iflat
               integer :: i,j,k
               real(WP) :: dx,x0,y0,z0,x1,y1,z1,x2,y2,z2
               real(WP) :: tilt(3),cx,sx,cy,sy,cz,sz
               ! Optional Tait-Bryan tilt in degrees (rotations about x, then y, then z)
               call param_read('Cube tilt', tilt, default=[0.0_WP,0.0_WP,0.0_WP])
               tilt=tilt*Pi/180.0_WP
               cx=cos(tilt(1)); sx=sin(tilt(1))
               cy=cos(tilt(2)); sy=sin(tilt(2))
               cz=cos(tilt(3)); sz=sin(tilt(3))
               dx=cube_L/real(cube_n,WP)
               if (amr%amRoot) then
                  ntot=int(cube_n,I8)**3
                  allocate(plist(ntot))
                  iflat=0_I8
                  do k=0,cube_n-1; do j=0,cube_n-1; do i=0,cube_n-1
                     iflat=iflat+1_I8
                     x0=-0.5_WP*cube_L+(real(i,WP)+0.5_WP)*dx
                     y0=-0.5_WP*cube_L+(real(j,WP)+0.5_WP)*dx
                     z0=-0.5_WP*cube_L+(real(k,WP)+0.5_WP)*dx
                     ! Rotate about x then y then z (intrinsic, applied to lattice positions)
                     x1=x0;          y1=cx*y0-sx*z0; z1=sx*y0+cx*z0
                     x2=cy*x1+sy*z1; y2=y1;          z2=-sy*x1+cy*z1
                     plist(iflat)%pos=[cz*x2-sz*y2, sz*x2+cz*y2, z2]
                     plist(iflat)%vel    =0.0_WP
                     plist(iflat)%F_bond =0.0_WP
                     plist(iflat)%F_fluid=0.0_WP
                     plist(iflat)%mw     =0.0_WP
                     plist(iflat)%dil    =0.0_WP
                     plist(iflat)%damage =0.0_WP
                     plist(iflat)%nb0    =0.0_WP
                     plist(iflat)%flag   =PART_MOVES+PART_INTEGRATES+PART_BONDS
                  end do; end do; end do
               else
                  ntot=0_I8
                  allocate(plist(0))
               end if
               call pd%append(plist,ntot)
               deallocate(plist)
            end block seed_cube
            ! Initial VF on the (currently base-only) mesh, then trigger an
            ! AMR regrid so any fine levels exist BEFORE bond_init runs. The
            ! chicken-and-egg: init_from_scratch fires tagging with VF=0 (no
            ! particles yet) so it only creates level 0; we now have particles
            ! and a non-trivial VF, so this regrid actually refines around the
            ! cube. post_regrid handles the particle/bond redistribute and
            ! re-computes VF on the new hierarchy.
            call pd%update_VF()
            call amr%regrid(baselvl=0,time=time%t)
            call pd%get_info()
         end if
      end block build_grid_and_state

      ! Build the initial bond network and stamp reference weighted volume.
      ! Must be called after particles are in place (post seed / post restart)
      ! and on their final AMR levels (post regrid).
      if (.not.restarted) call pd%bond_init()

      ! Seed VF so the first regrid (and any pre-advance plotfile) reflects
      ! the initial particle distribution.
      call pd%update_VF()

      ! Log initial state
      log_init: block
         character(len=str_long) :: message
         if (amr%amRoot) then
            if (restarted) then
               write(message,'("[amrpd] Restarted from ",a," with ",i0," particles at t=",es12.5)') trim(restart_dir),pd%np,time%t
            else
               write(message,'("[amrpd] Seeded ",i0," particles (",i0,"^3 cube of side ",f6.3,")")') pd%np,cube_n,cube_L
            end if
            call log(message)
         end if
      end block log_init

      ! Initialize particle visualization and dump the initial state
      init_viz: block
         call pviz%initialize(pd,name='amrpd')
         call pviz%select_comp('flag',on=.true.)   ! keep flag for debugging
         call pviz%select_comp('mw',  on=.true.)   ! weighted volume (stamped at bond_init)
         call pviz%select_comp('dil', on=.true.)   ! dilatation (recomputed per step)
         call pviz%select_comp('fbx', on=.true.)   ! bond force (recomputed per step)
         call pviz%select_comp('fby', on=.true.)
         call pviz%select_comp('fbz', on=.true.)
         call pviz%select_comp('damage', on=.true.)! per-particle damage fraction in [0,1]
         ! Eulerian-field viz (particle VF on the AMR mesh). Same name as the
         ! particle viz -- amrviz uses 'plt.NNNNNN' and amrpdviz uses
         ! 'plt.part.NNNNNN', so they live side-by-side in amrviz/amrpd/.
         call viz%initialize(amr,'amrpd',use_hdf5=.false.)
         call viz%add_scalar(pd%VF,1,'pVF')
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         if (viz_evt%occurs()) then
            call pviz%write(time=time%t)
            call viz%write(time=time%t)
         end if
         if (amr%amRoot) call log('[amrpd] Initial plotfiles written under amrviz/')
      end block init_viz

      ! Create monitor files. Fields are populated by get_info / get_cfl below
      ! before the first write so that columns 1) have valid data and 2) get
      ! their headers stamped at t=0.
      create_monitors: block
         ! Populate fields so the first monitor write has real data
         call pd%get_info()
         call pd%get_cfl(dt=time%dt,cfl=time%cfl)

         ! Overall simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,    'Timestep')
         call mfile%add_column(time%t,    'Time')
         call mfile%add_column(time%dt,   'dt')
         call mfile%add_column(time%cfl,  'CFL')
         call mfile%add_column(pd%np,     'Particle count')
         call mfile%add_column(pd%nb,        'Bond count')
         call mfile%add_column(pd%nb_broken, 'Bonds broken')
         call mfile%add_column(pd%Umin,   'Umin')
         call mfile%add_column(pd%Umax,   'Umax')
         call mfile%add_column(pd%Umean,  'Umean')
         call mfile%add_column(pd%Vmin,   'Vmin')
         call mfile%add_column(pd%Vmax,   'Vmax')
         call mfile%add_column(pd%Vmean,  'Vmean')
         call mfile%add_column(pd%Wmin,   'Wmin')
         call mfile%add_column(pd%Wmax,   'Wmax')
         call mfile%add_column(pd%Wmean,  'Wmean')
         call mfile%write()

         ! Per-direction CFL breakdown
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,    'Timestep')
         call cflfile%add_column(time%t,    'Time')
         call cflfile%add_column(time%dt,   'dt')
         call cflfile%add_column(pd%CFLp, 'CFLp')
         call cflfile%add_column(pd%CFLe, 'CFLe')
         call cflfile%add_column(pd%CFLc, 'CFLc')
         call cflfile%write()

         ! AMR grid stats
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,         'Timestep')
         call gridfile%add_column(time%t,         'Time')
         call gridfile%add_column(amr%nlevels,    'Nlvl')
         call gridfile%add_column(amr%nboxes,     'Nbox')
         call gridfile%add_column(amr%ncells,     'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,     'Maximum RSS')
         call gridfile%add_column(amr%minRSS,     'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,     'Average RSS')
         call gridfile%write()

         ! Load-balance metrics across ranks
         balancefile=monitor(amRoot=amr%amRoot,name='balance')
         call balancefile%add_column(time%n,    'Timestep')
         call balancefile%add_column(time%t,    'Time')
         call balancefile%add_column(pd%np_min, 'Np min')
         call balancefile%add_column(pd%np_max, 'Np max')
         call balancefile%add_column(pd%np_eff, 'Np eff')
         call balancefile%add_column(pd%nb_min, 'Nb min')
         call balancefile%add_column(pd%nb_max, 'Nb max')
         call balancefile%add_column(pd%nb_eff, 'Nb eff')
         call balancefile%write()
      end block create_monitors

      ! Initialize checkpoint save event. amrio writes header + scalars; pd
      ! writes particle/bond data alongside.
      init_checkpoint: block
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize regridding event (step-based, like amrcomp_drop).
      init_regridding: block
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper,default=huge(1))
      end block init_regridding
   end subroutine simulation_init


   !> Time integration loop. Each step: CFL-driven dt selection, Verlet advance
   !> under gravity, monitor + viz output. No bond forces yet (F_bond stays
   !> zero); cube just falls + bounces off walls.
   subroutine simulation_run()
      implicit none

      do while (.not.time%done())
         ! CFL-driven dt selection
         call pd%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Verlet advance
         call pd%advance(time%dt)

         ! Regrid if the event triggers (rebuilds AMR hierarchy via tagging,
         ! redistributes particles+bonds via post_regrid, refreshes VF)
         if (regrid_evt%occurs()) call amr%regrid(baselvl=0,time=time%t)

         ! Refresh stats and write monitors
         call pd%get_info()
         call mfile%write()
         call cflfile%write()
         call gridfile%write()
         call balancefile%write()

         ! Periodic visualization dump (particles + Eulerian VF)
         if (viz_evt%occurs()) then
            call pviz%write(time=time%t)
            call viz%write(time=time%t)
         end if

         ! Periodic checkpoint
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               character(len=str_long) :: ckdir
               ckdir='restart/amrpd_'//trim(adjustl(rtoa(time%t)))
               call io%write(dirname=trim(ckdir),time=time%t,step=time%n)
               call pd%write(dirname=trim(ckdir))
            end block save_checkpoint
         end if
      end do

      if (amr%amRoot) call log('[amrpd] Run complete')
   end subroutine simulation_run


   !> Finalization hook
   subroutine simulation_final()
      implicit none
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
      call balancefile%finalize()
      call save_evt%finalize()
      call regrid_evt%finalize()
      call io%finalize()
      call viz_evt%finalize()
      call pviz%finalize()
      call viz%finalize()
      call pd%finalize()
      call amr%finalize()
      call time%finalize()
   end subroutine simulation_final

end module simulation
