!> AMR LJCF - Incompressible cross-flow over a round liquid jet
!> Inflow/outflow in X, periodic in Y/Z
module simulation
   use precision,         only: WP
   use amrviz_class,      only: amrviz
   use amrgrid_class,     only: amrgrid
   use amrmpinc_class,    only: amrmpinc
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use messager,          only: log
   use amrio_class,       only: amrio
   use string,            only: str_medium
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Time integration
   type(timetracker) :: time

   ! Solver data
   type(amrmpinc), target :: fs
   type(amrdata) :: dQdt,Umag

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile

   ! Restart data
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   ! Physical parameters
   real(WP) :: Ujet
   real(WP) :: viscL_mol,viscG_mol
   real(WP), dimension(3) :: gravity

contains

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2) ! Enable 2D case
   end function sphere_levelset

   !> Levelset function for cylinder
   function cylinder_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP-abs(xyz(2)) ! Enable 2D case
   end function cylinder_levelset

   !> Compute viscosity
   subroutine get_viscosity()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Use harmonic averaging
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(viscL_mol,myeps)+(1.0_WP-pVF(i,j,k,1))/max(viscG_mol,myeps))
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosity

   !> Tagger for this case based on velocity gradient magnitude
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,dxi2,dyi2,dzi2,delta,delta2
      real(WP) :: lapU,lapV,lapW,u_sgs,Re
      integer :: i,j,k
      logical :: near_wall
      tags=tags_ptr
      ! Get mesh spacing
      dx=solver%amr%dx(lvl); dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Prevent maximum Re-driven refinement near wall
            near_wall=(solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.lt.solver%amr%xlo+2.0_WP*dx)
            if (near_wall.and.lvl.ge.solver%amr%maxlvl-1) cycle
            ! Laplacian of velocity Q=UVW
            lapU=(pQ(i+1,j,k,1)-2.0_WP*pQ(i,j,k,1)+pQ(i-1,j,k,1))*dxi2+(pQ(i,j+1,k,1)-2.0_WP*pQ(i,j,k,1)+pQ(i,j-1,k,1))*dyi2+(pQ(i,j,k+1,1)-2.0_WP*pQ(i,j,k,1)+pQ(i,j,k-1,1))*dzi2
            lapV=(pQ(i+1,j,k,2)-2.0_WP*pQ(i,j,k,2)+pQ(i-1,j,k,2))*dxi2+(pQ(i,j+1,k,2)-2.0_WP*pQ(i,j,k,2)+pQ(i,j-1,k,2))*dyi2+(pQ(i,j,k+1,2)-2.0_WP*pQ(i,j,k,2)+pQ(i,j,k-1,2))*dzi2
            lapW=(pQ(i+1,j,k,3)-2.0_WP*pQ(i,j,k,3)+pQ(i-1,j,k,3))*dxi2+(pQ(i,j+1,k,3)-2.0_WP*pQ(i,j,k,3)+pQ(i,j-1,k,3))*dyi2+(pQ(i,j,k+1,3)-2.0_WP*pQ(i,j,k,3)+pQ(i,j,k-1,3))*dzi2
            ! SGS Reynolds number
            u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
            Re=solver%rhoG*u_sgs*delta/viscG_mol
            if (Re.gt.Re_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Dirichlet BC: velocity inflow at xlo for U and ylo for V
   subroutine dirichlet_velocity(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      use mathtools, only: twoPi
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i,j,k
      real(WP) :: rad
      select case (face)
       case (1)  ! Inflow in X-
         select case (comp)
          case ('U')  ! Staggered U=Ujet
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rad=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (amr%nz.eq.1) rad=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2)
               if (rad.lt.0.5_WP) then
                  p(i,j,k,1)=Ujet
               else
                  p(i,j,k,1)=0.0_WP
               end if
            end do; end do; end do
          case ('V','W')  ! Staggered V,W = 0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')  ! Cell-centered: U=Ujet, V=0, W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rad=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (amr%nz.eq.1) rad=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2)
               if (rad.lt.0.5_WP) then
                  p(i,j,k,1)=Ujet
               else
                  p(i,j,k,1)=0.0_WP
               end if
               p(i,j,k,2)=0.0_WP
               p(i,j,k,3)=0.0_WP
            end do; end do; end do
         end select
       case (3)  ! Inflow in Y-
         select case (comp)
          case ('V')  ! Staggered V=1
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=1.0_WP
            end do; end do; end do
          case ('U','W')  ! Staggered U,W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')  ! Cell-centered: U=0, V=1, W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
               p(i,j,k,2)=1.0_WP
               p(i,j,k,3)=0.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine dirichlet_velocity

   !> User-defined VF BC - sets inlet ghost cells based on cylinder
   subroutine dirichlet_VF(solver,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
      use amrex_amr_module, only: amrex_box
      use mms_geom, only: initialize_volume_moments
      use amrmpinc_class,   only: VFlo,VFhi
      use amrvof_geometry,  only: get_plane_dist
      use mathtools,        only: normalize
      implicit none
      class(amrmpinc),  intent(inout) :: solver
      integer,          intent(in) :: lvl
      real(WP),         intent(in) :: time
      integer,          intent(in) :: face
      type(amrex_box),  intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      real(WP), dimension(3) :: BL,BG,lo,hi,nrm
      real(WP) :: dx,dy,dz,VF
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Only handle the xlo inflow
      if (face.ne.1) return
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Loop over provided ghost box
      do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
         ! Cell bounds
         lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz]
         hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz]
         ! Compute VF and barycenters from the jet column levelset
         call initialize_volume_moments(lo=lo,hi=hi,levelset=cylinder_levelset,time=time,level=nref,VFlo=VFlo,VF=VF,BL=BL,BG=BG)
         ! Store volume fraction at all levels
         pVF(i,j,k,1)=VF
         ! Store barycenters and PLIC at finest level only
         if (associated(pCL)) pCL(i,j,k,:)=BL
         if (associated(pCG)) pCG(i,j,k,:)=BG
         if (associated(pPLIC)) then
            if (VF.lt.VFlo.or.VF.gt.VFhi) then
               ! Pure cell: trivial plane
               pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0e10_WP,VF-0.5_WP)]
            else
               ! Mixed cell: normal points liquid->gas, distance from S&Z
               nrm=normalize(BG-BL)
               pPLIC(i,j,k,:)=[nrm(1),nrm(2),nrm(3),get_plane_dist(nrm,lo,hi,VF)]
            end if
         end if
      end do; end do; end do
   end subroutine dirichlet_VF

   !> User-provided initialization for jet
   subroutine jet_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpinc_class, only: VFlo
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pV,pQ
      real(WP), dimension(3) :: BL,BG  ! Dummy barycenters
      real(WP) :: dx,dy,dz,VF
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pV=>solver%V%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Set initial crossflow
            pV(i,j,k,1)=1.0_WP
            pQ(i,j,k,2)=1.0_WP
            if (i.lt.solver%amr%geom(lvl)%domain%lo(1)) then
               pV(i,j,k,1)=0.0_WP
               pQ(i,j,k,2)=0.0_WP
            end if
            ! Set interface
            if (i.lt.solver%amr%geom(lvl)%domain%lo(1)) then
               ! Compute VF and barycenters from levelset
               call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
               &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
               &                              levelset=cylinder_levelset,time=time,level=nref,VFlo=VFlo,VF=VF,BL=BL,BG=BG)
               ! Store volume fraction
               pVF(i,j,k,1)=VF
               ! Store barycenters
               if (lvl.eq.solver%amr%maxlvl) then
                  pCL(i,j,k,:)=BL
                  pCG(i,j,k,:)=BG
               end if
            else
               pVF(i,j,k,1)=0.0_WP
               if (lvl.eq.solver%amr%maxlvl) then
                  pCL(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               end if
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine jet_init

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none

      ! Create amrgrid
      create_amrgrid: block
         amr%name='LJCF'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo= 00.0_WP; amr%xhi=+20.0_WP
         amr%ylo=-05.0_WP; amr%yhi=+15.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Handle 2D case
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves here
      handle_restart: block
         integer :: restart_step
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=1)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time integration
      initialize_time: block
         ! Create time tracker and initialize
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_time

      ! Create flow solver
      create_flow_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrdata_class,    only: interp_face_lin
         use amrmpinc_class,   only: BC_GAS,BC_USER
         use amrmg_class,      only: amrmg_outer_pcg_mlmg
         ! Create flow solver
         call fs%initialize(amr,name='jet')
         ! Set initial conditions
         fs%user_init=>jet_init
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set densities
         fs%rhoG=1.0_WP; call param_read('Density ratio',fs%rhoL)
         ! Read in momentum flux ratio and set liquid velocity
         call param_read('Mom flux ratio',Ujet); Ujet=sqrt(Ujet/fs%rhoL)
         ! Set surface tension coefficient
         call param_read('Weber number',fs%sigma); fs%sigma=1.0_WP/fs%sigma
         ! Set molecular viscosities
         call param_read('Reynolds number',viscG_mol); viscG_mol=1.0_WP/viscG_mol
         call param_read('Viscosity ratio',viscL_mol); viscL_mol=viscG_mol*viscL_mol
         ! Set gravity
         gravity=0.0_WP; call param_read('Froude number',gravity(1),default=1.0e30_WP); gravity(1)=1.0_WP/gravity(1)**2
         ! Set pressure convergence
         fs%psolver%outer_solver=amrmg_outer_pcg_mlmg
         fs%psolver%tol_rel=1.0e-5_WP
         ! Dirichlet conditions for VOF at x- inlet and pure gas at y- inlet
         fs%lo_bc(1:2)=[BC_USER,BC_GAS]
         fs%user_vofbc=>dirichlet_VF
         ! Dirichlet conditions for velocities at x-/y- inlet
         fs%Q%lo_bc(1:2,:)=amrex_bc_ext_dir
         fs%U%lo_bc(1:2,1)=amrex_bc_ext_dir
         fs%V%lo_bc(1:2,1)=amrex_bc_ext_dir
         fs%W%lo_bc(1:2,1)=amrex_bc_ext_dir
         fs%user_bc=>dirichlet_velocity
         ! Neumann conditions for velocities at x+/y+ outlets
         fs%Q%hi_bc(1:2,:)=amrex_bc_foextrap
         fs%U%hi_bc(1:2,1)=amrex_bc_foextrap
         fs%V%hi_bc(1:2,1)=amrex_bc_foextrap
         fs%W%hi_bc(1:2,1)=amrex_bc_foextrap
      end block create_flow_solver

      ! Create workspace array
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=3,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Reynolds',Re_tag)
         ! Create initial grid from scratch or restore from checkpoint
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Create initial grid
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC
            call fs%build_plic(time%t)
            call fs%build_subVF()
            ! Initialize face velocities
            call fs%get_face_velocity()
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         end if
         ! Set viscosity: molecular + SGS
         call get_viscosity()
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_visualization: block
         ! Create visualization object
         call viz%initialize(amr,'jet',use_hdf5=.false.)
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%Q,1,'U')
         call viz%add_scalar(fs%Q,2,'V')
         call viz%add_scalar(fs%Q,3,'W')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_surfmesh(fs%smesh,'plic')
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
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%psolver%res,'Pressure residual')
         call mfile%add_column(fs%psolver%niter,'Pressure iterations')
         call mfile%add_column(fs%divmax,'Divergence')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLst,'CFLst')
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

         ! Store old interface and velocities
         call fs%store_old()

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity: U^{mid} = 0.5*(U + Uold)
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Increment velocity with advection+viscous terms
            call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
            call fs%Q%average_down(); call fs%Q%fill(time%t)

            ! Rebuild PLIC and sub-cell VF
            call fs%build_plic(time%t)
            call fs%build_subVF()

            ! Interpolate velocity to the faces
            call fs%get_face_velocity()

            ! Increment both velocities with current pressure term
            call fs%add_pressure(scale=time%dt,phi=fs%P,gravity=gravity)

            ! Add surface tension to both velocities
            call fs%add_surface_tension(scale=time%dt)

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Correct outflow for mass conservation
            call fs%correct_outflow()

            ! Prepare and solve pressure Poisson
            call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
            call fs%prepare_psolver()
            call fs%psolver%solve(rhs=fs%div)

            ! Correct both velocities with pressure increment
            call fs%add_pressure(scale=time%dt)

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
         call get_viscosity()
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)

         ! Monitor output
         call fs%get_info()
         call mfile%write()
         call cflfile%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/jet_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if
         
      end do

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
      call Umag%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final

end module simulation
