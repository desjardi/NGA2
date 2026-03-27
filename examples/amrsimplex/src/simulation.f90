!> AMR Simplex Atomizer
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
   use polygon_class,     only: polygon
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   ! Grid
   type(amrgrid), target :: amr

   ! Time integration
   type(timetracker) :: time

   ! Solver data
   type(amrmpinc), target :: fs
   type(amrdata) :: resUVW,Umag

   ! IB fluid volume fraction (1=fluid, 0=solid)
   type(polygon) :: poly
   type(amrdata), target :: VFib

   ! Visualization
   type(amrviz) :: viz
   type(event) :: viz_evt

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)

   ! Monitoring
   type(monitor) :: mfile,cflfile,gridfile,postproc

   ! Restart data
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   ! Physical parameters
   real(WP) :: viscL_mol,viscG_mol

   !> Inlet pipes geometry and flow rates
   real(WP) :: Rinlet=0.0023_WP
   real(WP) :: Rexit=0.00143_WP
   real(WP) :: Rpipe=0.000185_WP
   real(WP), dimension(3) :: p1=[-0.00442_WP,0.0_WP,+0.001245_WP]
   real(WP), dimension(3) :: p2=[-0.00442_WP,0.0_WP,-0.001245_WP]
   real(WP), dimension(3) :: n1=[+0.6_WP,-0.8_WP,0.0_WP]
   real(WP), dimension(3) :: n2=[+0.6_WP,+0.8_WP,0.0_WP]
   real(WP) :: mfr

   !> Post-processing info
   real(WP) :: liq_vol

contains

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
   end function sphere_levelset

   !> Levelset function for IB surface
   function simplex_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP), dimension(3) :: v,p
      ! Simplex polygon
      G=poly%get_distance([xyz(1),sqrt(xyz(2)**2+xyz(3)**2)])
      ! Add inlet pipes
      if (xyz(1).lt.p1(1).and.sqrt(xyz(2)**2+xyz(3)**2).lt.0.00203_WP) then
         v=xyz-p1; p=v-n1*dot_product(v,n1); G=max(G,Rpipe-sqrt(dot_product(p,p)))
         v=xyz-p2; p=v-n2*dot_product(v,n2); G=max(G,Rpipe-sqrt(dot_product(p,p)))
      end if
   end function simplex_levelset

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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,gradU_mag,Re_cell,dist
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
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! No refinement in the last 10% of the domain from the outflow
            if (solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.gt.solver%amr%xhi-0.1_WP*(solver%amr%xhi-solver%amr%xlo)) cycle
            ! Velocity gradient tensor
            gradU(1,1)=0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
            gradU(2,1)=0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            gradU(3,1)=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
            gradU(1,2)=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
            gradU(2,2)=0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
            gradU(3,2)=0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            gradU(1,3)=0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            gradU(2,3)=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
            gradU(3,3)=0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            ! |∇u| = sqrt(sum of all gradients squared)
            gradU_mag=sqrt(sum(gradU**2))
            ! Normalize into a local gas Reynolds number
            Re_cell=solver%rhoG*gradU_mag*solver%amr%min_meshsize(lvl)**2/viscG_mol
            ! Tagged based on cell Re value
            if (Re_cell.gt.Re_tag) tagarr(i,j,k,1)=SETtag
            ! Also tag near the IB surface
            dist=simplex_levelset([solver%amr%xlo+(real(i,WP)+0.5_WP)*dx, &
            &                      solver%amr%ylo+(real(j,WP)+0.5_WP)*dy, &
            &                      solver%amr%zlo+(real(k,WP)+0.5_WP)*dz],time)
            if (dist.lt.5.0_WP*dx.and.dist.gt.-dx) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Dirichlet BC: uniform inflow at 1 at xlo/xhi for U, 0 for V/W
   subroutine dirichlet_velocity(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      use mathtools,        only: Pi
      class(amrmpinc), intent(in) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), pointer, intent(inout) :: p
      integer :: i,j,k,ic
      real(WP), parameter :: Rin=0.00159_WP,Rout=0.00212_WP
      real(WP) :: Uin,rad
      ! Find component to modify
      if (size(p,4).eq.1) then; ic=1 ! Staggered velocity has one component
      else; ic=merge(1,merge(2,3,comp.eq.'V'),comp.eq.'U')  ! cell-centered: U→1, V→2, W→3
      end if
      ! Pick the x- face
      select case (face)
      case (1)  ! x-lo
         select case (comp)
          ! U=Uin
          case ('U')
            ! Get inflow velocity
            Uin=mfr/(solver%rhoL*Pi*(Rout**2-Rin**2))
            ! Apply
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rad=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*solver%amr%dy(lvl))**2 &
               &       +(solver%amr%zlo+(real(k,WP)+0.5_WP)*solver%amr%dz(lvl))**2)
               p(i,j,k,ic)=0.0_WP; if (rad.ge.Rin.and.rad.le.Rout) p(i,j,k,ic)=Uin
            end do; end do; end do
          ! V=W=0
          case ('V','W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,ic)=0.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine dirichlet_velocity

   !> User-provided initialization for VF
   subroutine init_VF(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrmpinc), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG
      real(WP) :: rad,dx,dy,dz
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get radial location
            rad=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
            ! Ensure the nozzle is filled with liquid up to the throat with wet walls
            if (solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.lt.-0.0015_WP.and.rad.le.Rinlet) then
               pVF(i,j,k,1)=1.0_WP
            else if (solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.ge.-0.0015_WP.and.solver%amr%xlo+(real(i,WP)+0.5_WP)*dx.lt.0.0_WP.and.rad.le.Rexit) then
               pVF(i,j,k,1)=1.0_WP
            else
               pVF(i,j,k,1)=0.0_WP
            end if
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VF

   !> User-defined VF BC - sets inlet ghost cells based on pipe geometry
   subroutine dirichlet_VF(solver,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
      use amrvof_class,     only: amrvof
      use amrex_amr_module, only: amrex_box
      implicit none
      class(amrvof),    intent(inout) :: solver
      integer,          intent(in) :: lvl
      real(WP),         intent(in) :: time
      integer,          intent(in) :: face
      type(amrex_box),  intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(inout) :: pVF,pCL,pCG,pPLIC
      real(WP) :: dx,dy,dz,rad
      integer  :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Loop over provided box
      do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
         ! Get radial location
         rad=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
         ! Set all passed variables
         if (rad.le.Rinlet) then
            pVF(i,j,k,1)=1.0_WP
            if (associated(pPLIC)) pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,+1.0e10_WP]
         else
            pVF(i,j,k,1)=0.0_WP
            if (associated(pPLIC)) pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,-1.0e10_WP]
         end if
         if (associated(pCL)) pCL(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
         if (associated(pCG)) pCG(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
      end do; end do; end do
   end subroutine dirichlet_VF

   !> Initialize IB fluid volume fraction from simplex levelset
   subroutine init_VFib(data,lvl,time,ba,dm)
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
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      real(WP), parameter :: VFlo=1.0e-12_WP
      integer, parameter :: nref=3
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=simplex_levelset,time=time,level=nref,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VFib

   !> Post-regrid dispatcher for automatic VFib filling
   subroutine vfib_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrdata), pointer :: this
      call c_f_pointer(ctx,this)
      call this%fill(time=time,lbase=lbase)
   end subroutine vfib_postregrid

   !> Initialization hook
   subroutine simulation_init()
      use param, only: param_read
      implicit none
      
      ! Create amrgrid
      create_amrgrid: block
         real(WP) :: xshift
         amr%name='amrsimplex'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('X shift',xshift)
         call param_read('Lx',amr%xhi); amr%xhi=amr%xhi-xshift; amr%xlo=-xshift
         call param_read('Ly',amr%yhi); amr%yhi=amr%yhi/2.0_WP; amr%ylo=-amr%yhi
         call param_read('Lz',amr%zhi); amr%zhi=amr%zhi/2.0_WP; amr%zlo=-amr%zhi
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
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
         use amrdata_class,    only: amrex_interp_face_linear
         use amrmpinc_class,   only: BC_USER
         use amrmg_class,      only: amrmg_outer_pcg_mlmg
         ! Create flow solver
         call fs%initialize(amr,name='simplex')
         ! Set initial conditions
         fs%user_mpinc_init=>init_VF
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=amrex_interp_face_linear
         ! Set densities
         call param_read('Liquid density',fs%rhoL)
         call param_read('Gas density'   ,fs%rhoG)
         ! Set surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Set molecular viscosities
         call param_read('Gas dynamic viscosity',viscG_mol)
         call param_read('Liquid dynamic viscosity',viscL_mol)
         ! Set pressure convergence
         fs%psolver%outer_solver=amrmg_outer_pcg_mlmg
         fs%psolver%tol_rel=1.0e-5_WP
         ! Dirichlet conditions for VOF at inlet
         fs%lo_bc(1)=BC_USER
         fs%user_vof_bc=>dirichlet_VF
         ! Dirichlet conditions for velocities at inlet
         fs%UVW%lo_bc(1,:)=amrex_bc_ext_dir
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir
         fs%user_mpinc_bc=>dirichlet_velocity
         ! Neumann conditions for velocities at outlet
         fs%UVW%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%hi_bc(1,1)=amrex_bc_foextrap
         ! Read in mass flow rate
         call param_read('Mass flow rate',mfr)
         ! Read in particle sub-stepping parameters
         call param_read('Particle dt_max', fs%dtmax, default=huge(1.0_WP))
         call param_read('Particle CFL max',fs%cflmax,default=huge(1.0_WP))
         fs%dt=fs%dtmax
      end block create_flow_solver

      ! Create workspace array
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call resUVW%initialize(amr,name='resUVW',ncomp=3,ng=0,interp=amrex_interp_none); call resUVW%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
      end block create_workspace

      ! Create IB fluid VF
      create_VFib: block
         use amrdata_class, only: amrex_interp_pc,amrex_bc_foextrap
         use iso_c_binding, only: c_loc
         ! Create polygon object
         call poly%initialize(nvert=15,name='simplex')
         poly%vert(:, 1)=[-0.10000_WP,0.00000_WP]
         poly%vert(:, 2)=[-0.00442_WP,0.00000_WP]
         poly%vert(:, 3)=[-0.00442_WP,0.00160_WP]
         poly%vert(:, 4)=[-0.00385_WP,0.00160_WP]
         poly%vert(:, 5)=[-0.00175_WP,0.00039_WP]
         poly%vert(:, 6)=[-0.00114_WP,0.00039_WP]
         poly%vert(:, 7)=[ 0.00000_WP,0.00143_WP]
         poly%vert(:, 8)=[ 0.00000_WP,0.00177_WP]
         poly%vert(:, 9)=[-0.00122_WP,0.00279_WP]
         poly%vert(:,10)=[-0.10000_WP,0.00279_WP]
         poly%vert(:,11)=[-0.10000_WP,0.00212_WP]
         poly%vert(:,12)=[-0.00543_WP,0.00212_WP]
         poly%vert(:,13)=[-0.00524_WP,0.00203_WP]
         poly%vert(:,14)=[-0.00634_WP,0.00159_WP]
         poly%vert(:,15)=[-0.10000_WP,0.00159_WP]
         ! Create VFib field with constant interpolation
         call VFib%initialize(amr,name='VFib',ncomp=1,ng=fs%nover,interp=amrex_interp_pc); call VFib%register()
         call amr%add_postregrid(vfib_postregrid,c_loc(VFib))
         VFib%user_init=>init_VFib
         VFib%lo_bc(1,1)=amrex_bc_foextrap
         VFib%hi_bc(1,1)=amrex_bc_foextrap
      end block create_VFib

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_mpinc_tagging=>my_tagger
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
            ! Initialize face velocities
            call fs%build_subVF()
            call fs%interp_vel_to_face()
            call fs%average_down_velocity()
            call fs%fill_velocity(time=time%t)
         end if
         ! Set viscosity: molecular + SGS
         call get_viscosity()
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
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
         call viz%initialize(amr,'simplex',use_hdf5=.false.)
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(VFib,1,'IB')
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
         ! Call post-processing routine
         call post_process()
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
         ! Create postproc monitor
         postproc=monitor(amRoot=amr%amRoot,name='postproc')
         call postproc%add_column(time%n,'Timestep')
         call postproc%add_column(time%t,'Time')
         call postproc%add_column(liq_vol,'Liquid volume')
         call postproc%write()
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
            call fs%UVW%lincomb(a=0.5_WP,src1=fs%UVWold,b=0.5_WP,src2=fs%UVW)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Increment velocity with advection+viscous terms
            call fs%get_dUVWdt(dUVWdt=resUVW,dt=time%dt,time=time%t)
            call fs%UVW%lincomb(a=1.0_WP,src1=fs%UVWold,b=time%dt,src2=resUVW)
            call fs%UVW%average_down(); call fs%UVW%fill(time%t)

            ! Rebuild PLIC and sub-cell VF
            call fs%build_plic(time%t)
            call fs%build_subVF()

            ! Interpolate velocity to the faces
            call fs%interp_vel_to_face()

            ! Increment both velocities with current pressure term
            call fs%correct_both_velocities(scale=time%dt,phi=fs%P)

            ! Add surface tension to both velocities
            call fs%add_surface_tension(scale=time%dt)

            ! Apply IB direct forcing
            call apply_ib_forcing()

            ! Average down and fill ghosts
            call fs%UVW%average_down(); call fs%UVW%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Correct outflow for mass conservation
            call fs%correct_outflow()

            ! Prepare and solve pressure Poisson
            call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
            call fs%prepare_psolver()
            call fs%psolver%solve(rhs=fs%div)

            ! Correct both velocities with pressure increment
            call fs%correct_both_velocities(scale=time%dt)

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
         call get_viscosity()
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)

         ! Monitor output
         call fs%get_info()
         call post_process()
         call mfile%write()
         call cflfile%write()
         call postproc%write()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/RSA_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if
         
      end do

   contains

      !> Apply IB direct forcing
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         implicit none
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pUVW,pVF
         integer :: i,j,k,lvl
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pU=>fs%U%mf(lvl)%dataptr(mfi)
               pV=>fs%V%mf(lvl)%dataptr(mfi)
               pW=>fs%W%mf(lvl)%dataptr(mfi)
               pUVW=>fs%UVW%mf(lvl)%dataptr(mfi)
               pVF=>VFib%mf(lvl)%dataptr(mfi)
               ! Force cell-centered velocity
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pUVW(i,j,k,:)=pVF(i,j,k,1)*pUVW(i,j,k,:)
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
      call resUVW%finalize()
      call Umag%finalize()
      call VFib%finalize()
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
      call postproc%finalize()
   end subroutine simulation_final

   !> Post-processing routine
   subroutine post_process()
      implicit none

      ! Get properly masked liquid volume
      get_liq_vol: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface,  only: amrmask_make_fine
         use parallel,         only: MPI_REAL_WP
         use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         integer :: lvl,i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVFib
         integer,  dimension(:,:,:,:), contiguous, pointer :: pMask
         liq_vol=0.0_WP
         do lvl=0,amr%clvl()
            if (lvl.lt.amr%clvl()) then
               call amrex_imultifab_build(mask,amr%ba(lvl),amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,amr%ba(lvl+1),[amr%rrefx(lvl),amr%rrefy(lvl),amr%rrefz(lvl)],0,1)
            end if
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
               pVFib=>VFib%mf(lvl)%dataptr(mfi)
               if (lvl.lt.amr%clvl()) pMask=>mask%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (lvl.lt.amr%clvl()) then
                     if (pMask(i,j,k,1).eq.0) cycle
                  end if
                  liq_vol=liq_vol+pVF(i,j,k,1)*pVFib(i,j,k,1)*amr%cell_vol(lvl)
               end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
            if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,liq_vol,1,MPI_REAL_WP,MPI_SUM,amr%comm)
      end block get_liq_vol

   end subroutine post_process

end module simulation
