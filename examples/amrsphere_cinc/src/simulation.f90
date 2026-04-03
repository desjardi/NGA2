!> AMR Sphere - Flow over a sphere with Immersed Boundary
!> Inflow/outflow in X, periodic in Y/Z
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrcinc_class,     only: amrcinc
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
   type(amrdata) :: dQdt,Umag

   ! Fluid volume fraction for IB forcing
   type(amrdata) :: VF

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile

   ! Physical parameters
   real(WP) :: visc_mol

contains

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)-0.5_WP
      if (amr%nz.eq.1) G=sqrt(xyz(1)**2+xyz(2)**2)-0.5_WP ! Enable 2D case
   end function sphere_levelset

   !> Tagger for this case based on velocity gradient magnitude and distance to sphere surface
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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dist,gradU_mag,Re_cell
      real(WP), dimension(3,3) :: gradU
      integer :: i,j,k
      tags=tags_ptr
      ! Get mesh spacing
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! No refinement 5D from the outflow
            if (solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.gt.solver%amr%xhi-5.0_WP) cycle
            ! Velocity gradient tensor
            gradU(1,1)=0.5_WP*dxi*(pQ(i+1,j,k,1)-pQ(i-1,j,k,1))
            gradU(2,1)=0.5_WP*dyi*(pQ(i,j+1,k,1)-pQ(i,j-1,k,1))
            gradU(3,1)=0.5_WP*dzi*(pQ(i,j,k+1,1)-pQ(i,j,k-1,1))
            gradU(1,2)=0.5_WP*dxi*(pQ(i+1,j,k,2)-pQ(i-1,j,k,2))
            gradU(2,2)=0.5_WP*dyi*(pQ(i,j+1,k,2)-pQ(i,j-1,k,2))
            gradU(3,2)=0.5_WP*dzi*(pQ(i,j,k+1,2)-pQ(i,j,k-1,2))
            gradU(1,3)=0.5_WP*dxi*(pQ(i+1,j,k,3)-pQ(i-1,j,k,3))
            gradU(2,3)=0.5_WP*dyi*(pQ(i,j+1,k,3)-pQ(i,j-1,k,3))
            gradU(3,3)=0.5_WP*dzi*(pQ(i,j,k+1,3)-pQ(i,j,k-1,3))
            ! |∇u| = sqrt(sum of all gradients squared)
            gradU_mag=sqrt(sum(gradU**2))
            ! Normalize into a local Reynolds number
            Re_cell=solver%rho*gradU_mag*solver%amr%min_meshsize(lvl)**2/visc_mol
            ! Tagged based on cell Re value
            if (Re_cell.gt.Re_tag) tagarr(i,j,k,1)=SETtag
            ! Also tag based on closeness to sphere surface
            dist=sphere_levelset([solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz],time)
            if (dist.lt.5.0_WP*dx.and.dist.gt.-dx) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Dirichlet BC: uniform inflow at 1 at xlo/xhi for U, 0 for V/W
   subroutine dirichlet_velocity(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrcinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), intent(inout) :: p
      integer :: i,j,k
      select case (face)
       case (1)  ! Inflow in X-
         select case (comp)
          case ('U')  ! Staggered U = 1
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=1.0_WP
            end do; end do; end do
          case ('V','W')  ! Staggered V,W = 0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')  ! Cell-centered: U=1, V=0, W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=1.0_WP
               p(i,j,k,2)=0.0_WP
               p(i,j,k,3)=0.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine dirichlet_velocity

   !> Initialize fluid volume fraction
   subroutine init_VF(data,lvl,time,ba,dm)
      use mms_geom, only: initialize_volume_moments
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP), dimension(3) :: BL,BG  ! Dummy barycenters
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=3,VFlo=1.0e-12_WP,VF=pVF(i,j,k,1),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VF

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         amr%name='amrsphere'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-05.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Handle 2D case
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
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
      
      ! Create flow solver
      create_flow_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrdata_class, only: interp_face_lin
         ! Create flow solver
         call fs%initialize(amr)
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set molecular viscosity
         call param_read('Reynolds number',visc_mol)
         visc_mol=1.0_WP/visc_mol
         ! Set pressure convergence
         fs%psolver%max_iter=20
         fs%psolver%tol_rel=1.0e-5_WP
         ! Set boundary conditions
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir
         fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%hi_bc(1,1)=amrex_bc_foextrap
         fs%user_bc=>dirichlet_velocity
      end block create_flow_solver

      ! Create workspace array
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=3,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
      end block create_workspace

      ! Create VF for IB forcing
      create_VF: block
         use amrdata_class, only: interp_reinit
         call VF%initialize(amr,name='VF',ncomp=1,ng=fs%nover,interp=interp_reinit)
         VF%user_init=>init_VF
         call VF%register()
      end block create_VF

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Reynolds',Re_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Initialize face velocities
         call fs%interp_vel_to_face()
         ! Set viscosity: molecular + SGS
         call fs%visc%setval(val=visc_mol)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)
      end block init_regridding

      ! Initialize visualization
      create_visualization: block
         ! Create visualization object
         call viz%initialize(amr,'amrsphere',use_hdf5=.false.)
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%Q,1,'U')
         call viz%add_scalar(fs%Q,2,'V')
         call viz%add_scalar(fs%Q,3,'W')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(VF,1,'VF')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_visualization

      ! Create monitor
      create_monitor: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(time%dt,time%cfl)
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
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Store old velocities
         call fs%Qold%copy(src=fs%Q)
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Increment velocity with advection+viscous terms
            call fs%get_dQdt(dQdt=dQdt)
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt/fs%rho,src2=dQdt)
            call fs%Q%average_down(); call fs%Q%fill(time%t)

            ! Interpolate velocity to the faces
            call fs%interp_vel_to_face()

            ! Increment both velocities with current pressure term
            call fs%correct_both_velocities(scale=time%dt/fs%rho,phi=fs%P)

            ! Apply IB direct forcing
            call apply_ib_forcing()

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Correct outflow for mass conservation
            call fs%correct_outflow(VF=VF)

            ! Solve pressure Poisson and increment pressure
            call fs%get_div(); call fs%div%mult(val=fs%rho/time%dt)
            call fs%psolver%solve(rhs=fs%div)

            ! Correct both velocities with new pressure increment
            call fs%correct_both_velocities(scale=time%dt/fs%rho)

            ! Add pressure increment
            call fs%P%add(src=fs%psolver%sol)

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
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
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)

         ! Monitor output
         call fs%get_info()
         call mfile%write()
         call cflfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         
      end do

   contains

      !> Apply IB direct forcing: multiply velocity by face-averaged VF
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         implicit none
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pQ,pVF
         integer :: i,j,k,lvl
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pU=>fs%U%mf(lvl)%dataptr(mfi)
               pV=>fs%V%mf(lvl)%dataptr(mfi)
               pW=>fs%W%mf(lvl)%dataptr(mfi)
               pQ=>fs%Q%mf(lvl)%dataptr(mfi)
               pVF=>VF%mf(lvl)%dataptr(mfi)
               ! Force cell-centered velocity
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pQ(i,j,k,:)=pVF(i,j,k,1)*pQ(i,j,k,:)
               end do; end do; end do
               ! Force face velocities
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pU(i,j,k,1)=0.5_WP*sum(pVF(i-1:i,j,k,1))*pU(i,j,k,1)
               end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pV(i,j,k,1)=0.5_WP*sum(pVF(i,j-1:j,k,1))*pV(i,j,k,1)
               end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pW(i,j,k,1)=0.5_WP*sum(pVF(i,j,k-1:k,1))*pW(i,j,k,1)
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
      call dQdt%finalize()
      call VF%finalize()
      call Umag%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final

end module simulation
