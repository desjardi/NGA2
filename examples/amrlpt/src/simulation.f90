!> AMRLPT test
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrlptviz_class,   only: amrlptviz
   use amrgrid_class,     only: amrgrid
   use amrcinc_class,     only: amrcinc
   use amrlpt_class,      only: amrlpt
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Time integration
   type(timetracker) :: time

   ! Solver data
   type(amrcinc), target :: fs
   type(amrdata) :: resUVW,Umag

   ! IB levelset and volume fraction
   type(amrdata) :: IB

   ! LPT solver
   type(amrlpt), target :: lpt

   ! Visualization
   type(amrviz) :: viz
   type(amrlptviz) :: lptviz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: vorticity_tag=huge(1.0_WP)

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile,partfile,lbfile

   ! Physical parameters
   real(WP) :: visc_mol

   ! Injection parameters
   real(WP) :: inj_mfr=0.0_WP             !< Mass flow rate
   real(WP) :: inj_dmean=0.0_WP           !< Mean particle diameter
   real(WP) :: inj_d=0.0_WP               !< Nozzle diameter (0=full y-z domain)
   real(WP), dimension(3) :: inj_pos=0.0_WP  !< Injection center
   real(WP), dimension(3) :: inj_vel=0.0_WP  !< Injection velocity
   real(WP) :: inj_residual=0.0_WP        !< Uninjected mass from previous step

contains

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)-0.5_WP
      if (amr%nz.eq.1) G=sqrt(xyz(1)**2+xyz(2)**2)-0.5_WP ! Enable 2D case
   end function sphere_levelset

   !> Tagger for this case
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pIB,pUVW
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP), dimension(3) :: vort
      integer :: i,j,k
      ! Resolve tags
      tags=tags_ptr
      ! Get mesh spacing
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      ! Traverse level
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pIB=>IB%mf(lvl)%dataptr(mfi)
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Tag based on closeness to sphere surface
            if (pIB(i,j,k,1).lt.2.0_WP*dx.and.pIB(i,j,k,1).gt.-dx) tagarr(i,j,k,1)=SETtag
            ! Tag based on vorticity magnitude
            vort(1)=(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))*0.5_WP*dyi-(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))*0.5_WP*dzi
            vort(2)=(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))*0.5_WP*dzi-(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))*0.5_WP*dxi
            vort(3)=(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))*0.5_WP*dxi-(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))*0.5_WP*dyi
            if (norm2(vort).gt.vorticity_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Initialize immersed boundary levelset and volume fraction
   subroutine init_IB(data,lvl,time,ba,dm)
      use mms_geom, only: initialize_volume_moments
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pIB
      real(WP), dimension(3) :: BL,BG  ! Dummy barycenters
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      real(WP), parameter :: VFlo=1.0e-12_WP
      integer, parameter :: nref=3
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pIB=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get levelset value
            pIB(i,j,k,1)=sphere_levelset([data%amr%xlo+(real(i,WP)+0.5_WP)*dx,data%amr%ylo+(real(j,WP)+0.5_WP)*dy,data%amr%zlo+(real(k,WP)+0.5_WP)*dz],time)
            ! Get IB volume fraction
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=nref,VFlo=VFlo,VF=pIB(i,j,k,2),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_IB

      !> Particle injection at a prescribed MFR from a circular nozzle
   subroutine my_inject(this,dt)
      use amrlpt_class, only: amrlpt,part,PART_MOVES,PART_COLLIDES,PART_EXCHANGES
      use precision,    only: WP,I8
      use mathtools,    only: Pi,twoPi
      use random,       only: random_uniform
      use messager,     only: warn
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP) :: Mgoal,Madded,dp,r,theta
      real(WP), dimension(3) :: pos,glo,ghi
      integer(I8) :: n_inj,ncap,j,ngather,ntries
      integer(I8), parameter :: max_tries=1000
      type(part), dimension(:), allocatable :: pnew,tmp,nearby
      logical :: overlap

      ! No injection if MFR is zero
      if (inj_mfr.le.0.0_WP) return

      ! Compute current injection goal
      Mgoal=inj_mfr*dt+inj_residual; Madded=0.0_WP; n_inj=0

      ! Gather existing particles near the injector for overlap checking
      glo(1)=inj_pos(1)-1.2_WP*inj_dmean; ghi(1)=inj_pos(1)+1.2_WP*inj_dmean
      glo(2)=inj_pos(2)-0.5_WP*inj_d-1.2_WP*inj_dmean; ghi(2)=inj_pos(2)+0.5_WP*inj_d+1.2_WP*inj_dmean
      glo(3)=inj_pos(3)-0.5_WP*inj_d-1.2_WP*inj_dmean; ghi(3)=inj_pos(3)+0.5_WP*inj_d+1.2_WP*inj_dmean
      if (this%amr%nz.eq.1) then; glo(3)=this%amr%zlo; ghi(3)=this%amr%zhi; end if
      call this%gather_region(glo,ghi,nearby,ngather)

      ! Only root injects
      if (this%amr%amRoot) then
         ncap=100; allocate(pnew(ncap))
         inject_loop: do while (Madded.lt.Mgoal)
            ! Increment counter
            n_inj=n_inj+1
            ! Resize if needed
            if (n_inj.gt.ncap) then
               ncap=2*ncap; allocate(tmp(ncap))
               tmp(1:n_inj-1)=pnew(1:n_inj-1)
               call move_alloc(tmp,pnew)
            end if
            ! Diameter
            dp=inj_dmean
            ! Position based on circular nozzle (slot in 2D)
            ntries=0
            retry: do
               ntries=ntries+1
               pos(1)=inj_pos(1)
               if (this%amr%nz.eq.1) then
                  pos(2)=random_uniform(lo=inj_pos(2)-0.5_WP*inj_d,hi=inj_pos(2)+0.5_WP*inj_d)
                  pos(3)=0.5_WP*(this%amr%zlo+this%amr%zhi)
               else
                  r=0.5_WP*inj_d*sqrt(random_uniform(lo=0.0_WP,hi=1.0_WP))
                  theta=random_uniform(lo=0.0_WP,hi=twoPi)
                  pos(2)=inj_pos(2)+r*sin(theta); pos(3)=inj_pos(3)+r*cos(theta)
               end if
               ! Overlap check with 20% margin
               overlap=.false.
               ! Check against newly injected particles
               do j=1,n_inj-1
                  if (norm2(pos-pnew(j)%pos).lt.0.6_WP*(dp+pnew(j)%d)) then; overlap=.true.; exit; end if
               end do
               ! Check against existing particles in the region
               if (.not.overlap) then
                  do j=1,ngather
                     if (norm2(pos-nearby(j)%pos).lt.0.6_WP*(dp+nearby(j)%d)) then; overlap=.true.; exit; end if
                  end do
               end if
               if (.not.overlap) exit retry
               if (ntries.ge.max_tries) exit retry
            end do retry
            ! If we exhausted retries, stop injecting this step
            if (overlap) then
               call warn('[Particle injection] Injector saturated — max overlap retries reached, deferring remaining mass')
               n_inj=n_inj-1_I8; exit inject_loop
            end if
            ! Add new particle
            pnew(n_inj)%d=dp
            pnew(n_inj)%pos=pos
            pnew(n_inj)%vel=inj_vel
            pnew(n_inj)%angVel=0.0_WP
            pnew(n_inj)%Acol=0.0_WP
            pnew(n_inj)%Tcol=0.0_WP
            pnew(n_inj)%dt=0.0_WP
            pnew(n_inj)%flag=PART_MOVES+PART_COLLIDES+PART_EXCHANGES
            ! Increment mass added
            Madded=Madded+this%rho*Pi/6.0_WP*dp**3
         end do inject_loop
         ! Update injection statistics (accumulate into _loc; get_info reduces and publishes)
         this%np_new_loc=this%np_new_loc+int(n_inj); this%Vp_new_loc=this%Vp_new_loc+Madded/this%rho
         ! Adjust residual
         inj_residual=Mgoal-Madded
      end if

      ! Clean up
      if (allocated(nearby)) deallocate(nearby)

      ! Add particles to LPT
      call this%append(pnew,n_inj)
      call this%redistribute()

   end subroutine my_inject

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         amr%name='amrlpt'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-5.0_WP; amr%xhi=+5.0_WP
         amr%ylo=-5.0_WP; amr%yhi=+5.0_WP
         amr%zlo=-5.0_WP; amr%zhi=+5.0_WP
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

      ! Initialize time integration
      initialize_time: block
         ! Create time tracker and initialize
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_time

      ! Initialize LPT solver
      init_lpt: block
         use amrlpt_class, only: AMRLPT_WALL,AMRLPT_OPEN
         ! Initialize solver
         call lpt%initialize(amr)
         ! Set particle parameters
         call param_read('Particle density',lpt%rho)
         lpt%gravity=[1.0_WP,0.0_WP,0.0_WP]
         call param_read('Particle max dt' ,lpt%dtmax ,default=huge(1.0_WP))
         call param_read('Particle max CFL',lpt%cflmax,default=time%cflmax)
         lpt%dt=lpt%dtmax
         ! Injection parameters
         call param_read('Particle mfr',inj_mfr)
         call param_read('Particle diameter',inj_dmean)
         inj_d=1.0_WP
         inj_pos=[amr%xlo+0.01_WP*(amr%xhi-amr%xlo),0.0_WP,0.0_WP]
         inj_vel=[1.0_WP,0.0_WP,0.0_WP]
         lpt%inject=>my_inject
         ! Filter width
         lpt%filter_width=7.0_WP*inj_dmean
         ! Collision parameters
         lpt%tau_col=5.0_WP*lpt%dtmax
         call param_read('Restitution coeff',lpt%e_n)
         call param_read('Restitution wall' ,lpt%e_w)
         call param_read('Friction coeff'   ,lpt%mu_f)
         ! Set walls all around
         lpt%lo_bc=AMRLPT_WALL
         lpt%hi_bc=AMRLPT_WALL
         if (amr%zper) then
            lpt%lo_bc(3)=AMRLPT_OPEN
            lpt%hi_bc(3)=AMRLPT_OPEN
         end if
      end block init_lpt

      ! Create flow solver
      create_flow_solver: block
         use amrex_amr_module, only: amrex_bc_reflect_odd,amrex_bc_int_dir
         use amrdata_class, only: amrex_interp_face_linear
         ! Create flow solver with same overlap as LPT
         fs%nover=lpt%nover; call fs%initialize(amr)
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=amrex_interp_face_linear
         ! Set molecular viscosity
         call param_read('Reynolds number',visc_mol)
         visc_mol=1.0_WP/visc_mol
         ! Set pressure convergence
         fs%psolver%max_iter=20
         fs%psolver%tol_rel=1.0e-5_WP
         ! Set boundary conditions
         fs%UVW%lo_bc(:,:)=amrex_bc_reflect_odd
         fs%UVW%hi_bc(:,:)=amrex_bc_reflect_odd
         fs%U%lo_bc(:,1)  =amrex_bc_reflect_odd
         fs%V%lo_bc(:,1)  =amrex_bc_reflect_odd
         fs%W%lo_bc(:,1)  =amrex_bc_reflect_odd
         fs%U%hi_bc(:,1)  =amrex_bc_reflect_odd
         fs%V%hi_bc(:,1)  =amrex_bc_reflect_odd
         fs%W%hi_bc(:,1)  =amrex_bc_reflect_odd
         if (amr%zper) then
            fs%UVW%lo_bc(3,:)=amrex_bc_int_dir; fs%UVW%hi_bc(3,:)=amrex_bc_int_dir
            fs%U%lo_bc(3,1)  =amrex_bc_int_dir; fs%U%hi_bc(3,1)  =amrex_bc_int_dir
            fs%V%lo_bc(3,1)  =amrex_bc_int_dir; fs%V%hi_bc(3,1)  =amrex_bc_int_dir
            fs%W%lo_bc(3,1)  =amrex_bc_int_dir; fs%W%hi_bc(3,1)  =amrex_bc_int_dir
         end if
      end block create_flow_solver

      ! Create workspace array
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call resUVW%initialize(amr,name='resUVW',ncomp=3,ng=0,interp=amrex_interp_none); call resUVW%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
      end block create_workspace

      ! Create IB: IB(1)=levelset, IB(2)=volume fraction
      create_IB: block
         use amrdata_class, only: amrex_interp_reinit
         call IB%initialize(amr,name='IB',ncomp=2,ng=fs%nover,interp=amrex_interp_reinit); call IB%register()
         IB%user_init=>init_IB
      end block create_IB

      ! Initialize regridding
      init_regridding: block
         ! Could set strategy to knapsack
         !amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging vorticity',vorticity_tag)
         call param_read('Tagging VF',lpt%VF_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Initialize face velocities
         call fs%interp_vel_to_face()
         ! Set viscosity: molecular + SGS
         call fs%visc%setval(val=visc_mol)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
      end block init_regridding

      ! Initialize visualization
      create_visualization: block
         ! Create visualization object
         call viz%initialize(amr,'amrlpt',use_hdf5=.false.)
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(IB,2,'IB')
         call viz%add_scalar(lpt%VF,1,'pVF')
         ! Create LPT visualization
         call lptviz%initialize(lpt,'amrlpt')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) then
            call viz%write(time=time%t)
            call lptviz%write(time=time%t)
         end if
      end block create_visualization

      ! Create monitor
      create_monitor: block
         ! Get solver info and cfl
         call lpt%get_info()
         call lpt%get_cfl(dt=time%dt,cflc=time%cfl,cfl=time%cfl)
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create particle monitor
         partfile=monitor(amroot=amr%amRoot,name='particles')
         call partfile%add_column(time%n,'Timestep number')
         call partfile%add_column(time%t,'Time')
         call partfile%add_column(time%dt,'Timestep size')
         call partfile%add_column(lpt%np,'Particle number')
         call partfile%add_column(lpt%np_new,'Npart new')
         call partfile%add_column(lpt%np_out,'Npart removed')
         call partfile%add_column(lpt%ncol,'Particle collisions')
         call partfile%add_column(lpt%VFmax,'Max VF')
         call partfile%add_column(lpt%VFmean,'Mean VF')
         call partfile%add_column(lpt%Umin,'Particle Umin')
         call partfile%add_column(lpt%Umax,'Particle Umax')
         call partfile%add_column(lpt%Vmin,'Particle Vmin')
         call partfile%add_column(lpt%Vmax,'Particle Vmax')
         call partfile%add_column(lpt%Wmin,'Particle Wmin')
         call partfile%add_column(lpt%Wmax,'Particle Wmax')
         call partfile%add_column(lpt%dmin,'Particle dmin')
         call partfile%add_column(lpt%dmax,'Particle dmax')
         call partfile%write()
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(fs%CFL,'CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%psolver%res,'Pressure residual')
         call mfile%add_column(fs%psolver%niter,'Pressure iterations')
         call mfile%add_column(fs%divmax,'Divergence')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create load balance monitor
         lbfile=monitor(amRoot=amr%amRoot,name='balance')
         call lbfile%add_column(time%n,'Timestep')
         call lbfile%add_column(time%t,'Time')
         call lbfile%add_column(lpt%np_min,'Np min')
         call lbfile%add_column(lpt%np_max,'Np max')
         call lbfile%add_column(lpt%np_eff,'Np eff')
         call lbfile%add_column(lpt%tmr_coll_%efficiency,'Coll eff')
         call lbfile%add_column(lpt%tmr_step_%efficiency,'Step eff')
         call lbfile%add_column(lpt%tmr_coll%tmax,'Coll tmax')
         call lbfile%add_column(lpt%tmr_coll_%tmax,'Coll_ tmax')
         call lbfile%add_column(lpt%tmr_fill%tmax,'Fill tmax')
         call lbfile%add_column(lpt%tmr_nbl%tmax,'NBL tmax')
         call lbfile%add_column(lpt%tmr_step%tmax,'Step tmax')
         call lbfile%add_column(lpt%tmr_step_%tmax,'Step_ tmax')
         call lbfile%add_column(lpt%tmr_vf%tmax,'VF tmax')
         call lbfile%add_column(lpt%tmr_src%tmax,'Src tmax')
         call lbfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitor

   end subroutine simulation_init

   !> Run the simulation
   subroutine simulation_run()
      implicit none

      ! Time integration loop
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Stop injecting after time of 10
         if (time%t.gt.10.0_WP) inj_mfr=0.0_WP

         ! Advance particles to current time
         call lpt%advance_to(time=time%t,do_collide=.true.,U=fs%U,V=fs%V,W=fs%W,cst_rho=fs%rho,cst_visc=visc_mol,Gib=IB,Gibcomp=1)

         ! Store old velocities
         call fs%UVWold%copy(src=fs%UVW)
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%UVW%lincomb(a=0.5_WP,src1=fs%UVWold,b=0.5_WP,src2=fs%UVW)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Increment velocity with advection+viscous terms
            call fs%get_dmomdt(UVW=fs%UVW,U=fs%U,V=fs%V,W=fs%W,dmomdt=resUVW)
            call fs%UVW%lincomb(a=1.0_WP,src1=fs%UVWold,b=time%dt/fs%rho,src2=resUVW)
            call fs%UVW%average_down(); call fs%UVW%fill(time%t)

            ! Interpolate velocity to the faces
            call fs%interp_vel_to_face()

            ! Increment both velocities with current pressure term
            call fs%correct_both_velocities(scale=time%dt/fs%rho,phi=fs%P)

            ! Apply IB direct forcing
            call apply_ib_forcing()

            ! Average down and fill ghosts
            call fs%UVW%average_down(); call fs%UVW%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Correct outflow for mass conservation
            call fs%correct_outflow()

            ! Solve pressure Poisson and increment pressure
            call fs%get_div(); call fs%div%mult(val=fs%rho/time%dt)
            call fs%psolver%solve(rhs=fs%div)

            ! Correct both velocities with new pressure increment
            call fs%correct_both_velocities(scale=time%dt/fs%rho)

            ! Add pressure increment
            call fs%P%add(src=fs%psolver%sol)

            ! Average down and fill ghosts
            call fs%UVW%average_down(); call fs%UVW%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Update viscosity
         call fs%visc%setval(val=visc_mol)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)

         ! Monitor output
         call lpt%get_info()
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call partfile%write()
         call lbfile%write()

         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time=time%t)
            call lptviz%write(time=time%t)
         end if

      end do

   contains

      !> Apply IB direct forcing
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         implicit none
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pUVW,pIB
         integer :: i,j,k,lvl
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointer to volume fraction
               pIB=>IB%mf(lvl)%dataptr(mfi)
               ! Force cell-centered velocity
               pUVW=>fs%UVW%mf(lvl)%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pUVW(i,j,k,:)=pIB(i,j,k,2)*pUVW(i,j,k,:)
               end do; end do; end do
               ! Force face-centered velocity
               pU=>fs%U%mf(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pU(i,j,k,1)=0.5_WP*sum(pIB(i-1:i,j,k,2))*pU(i,j,k,1)
               end do; end do; end do
               pV=>fs%V%mf(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pV(i,j,k,1)=0.5_WP*sum(pIB(i,j-1:j,k,2))*pV(i,j,k,1)
               end do; end do; end do
               pW=>fs%W%mf(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pW(i,j,k,1)=0.5_WP*sum(pIB(i,j,k-1:k,2))*pW(i,j,k,1)
               end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
         end do
      end subroutine apply_ib_forcing

   end subroutine simulation_run

   !> Finalization hook
   subroutine simulation_final()
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call resUVW%finalize()
      call IB%finalize()
      call Umag%finalize()
      call lpt%finalize()
      ! Finalize visualization
      call viz%finalize()
      call lptviz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
      call partfile%finalize()
      call lbfile%finalize()
   end subroutine simulation_final

end module simulation
