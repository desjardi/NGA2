!> AMR compressible flow two-way coupled to a peridynamics solid (amrpd)
module simulation
   use precision,           only: WP,I8
   use amrgrid_class,       only: amrgrid
   use amrcomp_class,       only: amrcomp
   use amrpd_class,         only: amrpd,part,part_gid,PART_MOVES,PART_INTEGRATES,PART_BONDS,PART_IS_DEAD,AMRPD_OPEN
   use pdsolver_class,      only: pdsolver,pd_partition
   use amrpdviz_class,      only: amrpdviz
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use stiffened_gas_class, only: stiffened_gas
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible solver
   type(timetracker) :: time
   type(amrcomp), target :: fs
   type(amrdata) :: dQdt
   type(amrdata) :: Umag,Mach

   !> Peridynamics solid and its mesh-deposited velocity field
   type(amrpd), target :: apd

   !> The solid solver (grid-free physics); apd above is its grid-side face
   type(pdsolver) :: pd
   type(amrpdviz) :: pviz
   type(amrdata) :: Usolid          !< Solid velocity on the AMR mesh (3 comp)
   type(amrdata) :: VFf             !< Fluid volume fraction = 1 - apd%VF (amrcomp convention)
   type(amrdata) :: dStress         !< Divergence of the fluid stress tensor (3 comp; interpolated to particles as F_fluid)

   !> Coupling-direction switches (default fully coupled; for isolation tests)
   logical :: couple_s2f=.true.     !< Solid->fluid (IB forcing of the flow)
   logical :: couple_f2s=.true.     !< Fluid->solid (F_fluid reaction on particles)

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,pdfile

   !> Material
   type(stiffened_gas), target :: fluid

   !> Quiescent gas state (uniform, isothermal high-T pellet environment)
   real(WP) :: rho0,p0,e0             !< Ambient density, pressure, specific internal energy
   real(WP) :: Reynolds,Prandtl       !< Viscous parameters

   !> Isothermal degassing source (gas injected into open pores)
   real(WP) :: src_mdot               !< Mass injection rate per unit volume
   real(WP) :: src_e                  !< Specific internal energy of injected (hot) gas
   real(WP) :: src_VFthresh=0.5_WP    !< Inject only where fluid fraction exceeds this (open pore)

   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Sponge parameters (radial about the pellet, in x-y)
   real(WP) :: R_spg=3.5_WP
   real(WP) :: L_spg=1.0_WP

   !> Net body force on the IB (rough Cd-like diagnostic)
   real(WP), dimension(3) :: Fib

   !> Solid pellet initialization (porous disk)
   real(WP) :: R_disk=1.0_WP          !< Pellet radius
   integer  :: n_pore=40              !< Number of pre-seeded pores
   real(WP) :: pore_rmin=0.05_WP      !< Min pore radius
   real(WP) :: pore_rmax=0.10_WP      !< Max pore radius
   real(WP) :: elem_size=0.0_WP       !< Solid element (particle) spacing dp

   !> PD sub-steps taken per fluid step (subcycling)
   integer :: n_sub=1

contains

   !> Compute viscosity using Sutherland's law, zero bulk viscosity, and set diffusivity based on Prandtl number
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pT,pQ,pVisc,pBeta,pDiff
      real(WP) :: r_cyl,blend,nu_spg
      real(WP), parameter :: Tmax_visc=10.0_WP
      real(WP), parameter :: myeps=1.0e-15_WP
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP), parameter :: Cdiff=0.1_WP
      ! Get maximum allowable kinematic viscosity in the sponge at finest level
      nu_spg=max_cfl*amr%min_meshsize(amr%clvl())**2/(4.0_WP*time%dt)
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pT=>fs%T%mf(lvl)%dataptr(mfi)
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Sutherland's law
               pVisc(i,j,k,1)=(1.0_WP+Suth_T)*min(pT(i,j,k,1),Tmax_visc)**Suth_n/(Reynolds*(min(pT(i,j,k,1),Tmax_visc)+Suth_T))
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Heat diffusivity: k = Cp*mu/Pr = Cv*Gamma*mu/Pr
               pDiff(i,j,k,1)=fluid%gamma*fluid%cv*pVisc(i,j,k,1)/Prandtl
               ! Apply sponge layer viscosity
               r_cyl=sqrt((amr%xlo+(real(i,WP)+0.5_WP)*amr%dx(lvl))**2+(amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2)
               if (r_cyl.gt.R_spg) then
                  blend=min((r_cyl-R_spg)/L_spg,1.0_WP)**2
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*nu_spg*pQ(i,j,k,1))
                  pDiff(i,j,k,1)=max(pDiff(i,j,k,1),Cdiff*blend*nu_spg*pQ(i,j,k,1))
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback - uniform quiescent gas everywhere
   subroutine quiescent_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%growntilebox(solver%nover)
         pQ(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),1)=rho0
         pQ(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),2)=0.0_WP
         pQ(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),3)=0.0_WP
         pQ(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),4)=0.0_WP
         pQ(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),5)=rho0*e0
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine quiescent_init

   !> Tagger based on velocity and density laplacians. Refinement around the
   !> solid body is handled separately by amrpd's own VF-based tagging callback
   !> (registered in apd%initialize via apd%VF_tag).
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,dxi2,dyi2,dzi2,delta,delta2
      real(WP) ::  rho_cc, rho_xp, rho_xm, rho_yp, rho_ym, rho_zp, rho_zm
      real(WP) :: irho_cc,irho_xp,irho_xm,irho_yp,irho_ym,irho_zp,irho_zm
      real(WP) :: lapU,lapV,lapW,u_sgs,Re,lapRHO,avgRHO,r_cyl
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Sponge check
            r_cyl=sqrt((solver%amr%xlo+(real(i,WP)+0.5_WP)*dx)**2+(solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2)
            ! Get local densities and inverse
            rho_cc=max(pQ(i  ,j  ,k  ,1),solver%rho_floor); irho_cc=1.0_WP/rho_cc
            rho_xp=max(pQ(i+1,j  ,k  ,1),solver%rho_floor); irho_xp=1.0_WP/rho_xp
            rho_xm=max(pQ(i-1,j  ,k  ,1),solver%rho_floor); irho_xm=1.0_WP/rho_xm
            rho_yp=max(pQ(i  ,j+1,k  ,1),solver%rho_floor); irho_yp=1.0_WP/rho_yp
            rho_ym=max(pQ(i  ,j-1,k  ,1),solver%rho_floor); irho_ym=1.0_WP/rho_ym
            rho_zp=max(pQ(i  ,j  ,k+1,1),solver%rho_floor); irho_zp=1.0_WP/rho_zp
            rho_zm=max(pQ(i  ,j  ,k-1,1),solver%rho_floor); irho_zm=1.0_WP/rho_zm
            ! Laplacian of velocity (Q components 2,3,4 = rhoU,rhoV,rhoW)
            lapU=(pQ(i+1,j,k,2)*irho_xp-2.0_WP*pQ(i,j,k,2)*irho_cc+pQ(i-1,j,k,2)*irho_xm)*dxi2 &
            &   +(pQ(i,j+1,k,2)*irho_yp-2.0_WP*pQ(i,j,k,2)*irho_cc+pQ(i,j-1,k,2)*irho_ym)*dyi2 &
            &   +(pQ(i,j,k+1,2)*irho_zp-2.0_WP*pQ(i,j,k,2)*irho_cc+pQ(i,j,k-1,2)*irho_zm)*dzi2
            lapV=(pQ(i+1,j,k,3)*irho_xp-2.0_WP*pQ(i,j,k,3)*irho_cc+pQ(i-1,j,k,3)*irho_xm)*dxi2 &
            &   +(pQ(i,j+1,k,3)*irho_yp-2.0_WP*pQ(i,j,k,3)*irho_cc+pQ(i,j-1,k,3)*irho_ym)*dyi2 &
            &   +(pQ(i,j,k+1,3)*irho_zp-2.0_WP*pQ(i,j,k,3)*irho_cc+pQ(i,j,k-1,3)*irho_zm)*dzi2
            lapW=(pQ(i+1,j,k,4)*irho_xp-2.0_WP*pQ(i,j,k,4)*irho_cc+pQ(i-1,j,k,4)*irho_xm)*dxi2 &
            &   +(pQ(i,j+1,k,4)*irho_yp-2.0_WP*pQ(i,j,k,4)*irho_cc+pQ(i,j-1,k,4)*irho_ym)*dyi2 &
            &   +(pQ(i,j,k+1,4)*irho_zp-2.0_WP*pQ(i,j,k,4)*irho_cc+pQ(i,j,k-1,4)*irho_zm)*dzi2
            ! SGS Reynolds number
            u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
            Re=Reynolds*u_sgs*delta
            if (Re.gt.Re_tag.and.(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)) tagarr(i,j,k,1)=SETtag
            ! Normalized density Laplacian
            lapRHO=(rho_xp-2.0_WP*rho_cc+rho_xm)*dxi2+(rho_yp-2.0_WP*rho_cc+rho_ym)*dyi2+(rho_zp-2.0_WP*rho_cc+rho_zm)*dzi2
            avgRHO=(rho_cc+rho_xp+rho_xm+rho_yp+rho_ym+rho_zp+rho_zm)/7.0_WP
            lapRHO=abs(lapRHO)*delta2/avgRHO
            if (lapRHO.gt.Rho_tag.and.(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Seed a 2D porous disk of peridynamics particles, centered at the origin.
   !> A set of random circular pores is carved out of the lattice; gas later
   !> degasses into these pores. 2D: single z-layer; 3D: extruded in z.
   subroutine seed_body()
      type(part), dimension(:), allocatable, target :: plist
      integer(I8) :: ntot,iflat
      integer :: i,j,k,nr,khi,ip,pass,sz
      integer, dimension(:), allocatable :: seed
      real(WP) :: x,y,zb,u
      real(WP), dimension(:), allocatable :: pcx,pcy,pcr
      logical :: in_pore
      nr=nint(R_disk/elem_size)
      khi=0; if (amr%nz.gt.1) khi=nr
      if (amr%amRoot) then
         ! Pre-generate random pores (fixed seed -> reproducible across runs)
         allocate(pcx(n_pore),pcy(n_pore),pcr(n_pore))
         call random_seed(size=sz); allocate(seed(sz))
         do i=1,sz; seed(i)=12345+37*i; end do
         call random_seed(put=seed)
         do ip=1,n_pore
            do
               call random_number(u); pcx(ip)=(2.0_WP*u-1.0_WP)*R_disk
               call random_number(u); pcy(ip)=(2.0_WP*u-1.0_WP)*R_disk
               call random_number(u); pcr(ip)=pore_rmin+u*(pore_rmax-pore_rmin)
               ! Keep the whole pore inside the disk
               if (sqrt(pcx(ip)**2+pcy(ip)**2)+pcr(ip).lt.R_disk) exit
            end do
         end do
         ! Pass 1 counts surviving particles, pass 2 fills (pore count is data-dependent)
         do pass=1,2
            iflat=0_I8
            do k=-khi,khi; do j=-nr,nr; do i=-nr,nr
               x=real(i,WP)*elem_size; y=real(j,WP)*elem_size
               if (x*x+y*y.gt.R_disk**2) cycle             ! outside the disk
               in_pore=.false.
               do ip=1,n_pore
                  if ((x-pcx(ip))**2+(y-pcy(ip))**2.lt.pcr(ip)**2) then; in_pore=.true.; exit; end if
               end do
               if (in_pore) cycle                          ! inside a pore -> void
               iflat=iflat+1_I8
               if (pass.eq.2) then
                  zb=real(k,WP)*elem_size; if (amr%nz.eq.1) zb=0.0_WP
                  plist(iflat)%pos    =[x,y,zb]
                  plist(iflat)%vel    =0.0_WP
                  plist(iflat)%F_bond =0.0_WP
                  plist(iflat)%F_fluid=0.0_WP
                  plist(iflat)%mw     =0.0_WP
                  plist(iflat)%dil    =0.0_WP
                  plist(iflat)%damage =0.0_WP
                  plist(iflat)%nb0    =0.0_WP
                  plist(iflat)%flag   =PART_MOVES+PART_INTEGRATES+PART_BONDS
               end if
            end do; end do; end do
            if (pass.eq.1) then; ntot=iflat; allocate(plist(ntot)); end if
         end do
      else
         ntot=0_I8
         allocate(plist(0))
      end if
      call apd%append(plist,ntot)
      deallocate(plist)
   end subroutine seed_body

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS, ambient state, and degassing source
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: Gamma,Cv,T0,Tsrc
         ! EoS
         call param_read('Gamma',Gamma)
         ! Quiescent ambient gas (reference temperature T0=1 sets Cv -> e0=Cv*T0)
         call param_read('Ambient pressure',p0)
         call param_read('Ambient density', rho0)
         T0=1.0_WP
         Cv=p0/(rho0*(Gamma-1.0_WP))
         e0=Cv*T0
         ! Build material
         call fluid%initialize(gamma=Gamma,pinf=0.0_WP,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
         ! Isothermal degassing source: mass rate (per volume) of gas at hot temperature Tsrc
         call param_read('Source rate',        src_mdot)
         call param_read('Source temperature', Tsrc)
         src_e=Cv*Tsrc
         call param_read('Source VF threshold',src_VFthresh,default=0.5_WP)
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Sutherland exponent',Suth_n,default=1.5_WP)
         call param_read('Sutherland temperature',Suth_T,default=0.4042_WP)
         ! Log conditions
         write(message,'("[Ambient] rho0=",es12.5," p0=",es12.5," e0=",es12.5)') rho0,p0,e0; call log(message)
         write(message,'("[Source]  mdot=",es12.5," e_inj=",es12.5)') src_mdot,src_e; call log(message)
         call fluid%print()
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         real(WP) :: dp_slab
         amr%name='amrcomp_devol'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-05.0_WP; amr%xhi=+05.0_WP
         amr%ylo=-05.0_WP; amr%yhi=+05.0_WP
         amr%zlo=-05.0_WP; amr%zhi=+05.0_WP
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Handle 2D case: slab thickness = particle size so a single z-layer of
         ! particles fills the cell (VF -> 1 in pseudo-2D)
         if (amr%nz.eq.1) then
            call param_read('Element size',dp_slab)
            amr%zlo=-0.5_WP*dp_slab
            amr%zhi=+0.5_WP*dp_slab
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker

      ! Initialize compressible solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap
         use amrdata_class, only: interp_face_lin
         ! Assign material and create flow solver
         call param_read('Use projection',fs%use_projection)
         fs%mat=>fluid; call fs%initialize(amr=amr)
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set pressure convergence
         fs%psolver%max_iter=20
         fs%psolver%tol_rel=1.0e-5_WP
         fs%psolver%verbose=2
         ! Set initial conditions
         fs%user_init=>quiescent_init
         ! Outflow (zero-gradient) on x and y; z is periodic. No user_bc needed.
         fs%Q%lo_bc(1,:)=amrex_bc_foextrap; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%Q%lo_bc(2,:)=amrex_bc_foextrap; fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_foextrap; fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%U%lo_bc(2,1)=amrex_bc_foextrap; fs%U%hi_bc(2,1)=amrex_bc_foextrap
         fs%V%lo_bc(1,1)=amrex_bc_foextrap; fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%lo_bc(2,1)=amrex_bc_foextrap; fs%V%hi_bc(2,1)=amrex_bc_foextrap
         fs%W%lo_bc(1,1)=amrex_bc_foextrap; fs%W%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%lo_bc(2,1)=amrex_bc_foextrap; fs%W%hi_bc(2,1)=amrex_bc_foextrap
      end block create_solver

      ! Initialize the peridynamics solid (containers + AMR callbacks)
      init_solid: block
         call apd%initialize(amr,name='apd')
         ! Material parameters
         call param_read('Material density',  apd%rho)
         call param_read('Elastic modulus',   apd%elastic_modulus)
         call param_read('Poisson ratio',     apd%poisson_ratio)
         call param_read('Critical energy',   apd%crit_energy,default=huge(1.0_WP))
         ! Direct failure-stretch override (huge -> use G_c-derived s0; finite -> ductile)
         call param_read('Failure stretch',   apd%fail_stretch,default=huge(1.0_WP))
         ! Maxwell deviatoric relaxation time (huge -> elastic-brittle; small -> ductile/viscous flow)
         call param_read('Relaxation time',   apd%tau,        default=huge(1.0_WP))
         ! SLS relaxing fraction (1 = pure Maxwell/full flow; <1 keeps long-term elastic stiffness)
         call param_read('Relaxation fraction',apd%visc_lambda,default=1.0_WP)
         ! Solid spacing dp; horizon defaults to 3.0125*dp (Peridigm convention)
         call param_read('Element size',      elem_size)
         call param_read('Horizon',           apd%delta,default=3.0125_WP*elem_size)
         apd%dV=elem_size**3
         ! Porous-disk geometry
         call param_read('Disk radius',    R_disk,   default=1.0_WP)
         call param_read('Number of pores',n_pore,   default=40)
         call param_read('Pore radius min',pore_rmin,default=0.05_WP)
         call param_read('Pore radius max',pore_rmax,default=0.10_WP)
         ! Gravity off by default (body driven by the flow)
         call param_read('Gravity',apd%gravity,default=[0.0_WP,0.0_WP,0.0_WP])
         ! Open domain BCs (y/z periodicity is handled by AMReX)
         apd%lo_bc=AMRPD_OPEN; apd%hi_bc=AMRPD_OPEN
         ! Refine the AMR mesh wherever the solid volume fraction exceeds VF_tag
         call param_read('Tagging VF',apd%VF_tag,default=0.1_WP)
         ! Coupling-direction switches (absent -> fully coupled)
         call param_read('Couple solid to fluid',couple_s2f,default=.true.)
         call param_read('Couple fluid to solid',couple_f2s,default=.true.)
      end block init_solid

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         use amrex_amr_module, only: amrex_bc_foextrap
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
         ! Solid velocity on the mesh (ghosts for the IB face interpolation).
         ! interp_none: recomputed each step, so allocate-don't-fill on regrid.
         call Usolid%initialize(amr,name='Usolid',ncomp=3,ng=fs%nover,interp=interp_none); call Usolid%register()
         if (.not.amr%xper) then; Usolid%lo_bc(1,:)=amrex_bc_foextrap; Usolid%hi_bc(1,:)=amrex_bc_foextrap; end if
         if (.not.amr%yper) then; Usolid%lo_bc(2,:)=amrex_bc_foextrap; Usolid%hi_bc(2,:)=amrex_bc_foextrap; end if
         ! Fluid volume fraction = 1 - apd%VF (amrcomp's convention)
         call VFf%initialize(amr,name='VFf',ncomp=1,ng=fs%nover,interp=interp_none); call VFf%register()
         if (.not.amr%xper) then; VFf%lo_bc(1,1)=amrex_bc_foextrap; VFf%hi_bc(1,1)=amrex_bc_foextrap; end if
         if (.not.amr%yper) then; VFf%lo_bc(2,1)=amrex_bc_foextrap; VFf%hi_bc(2,1)=amrex_bc_foextrap; end if
         ! Fluid stress-tensor divergence, interpolated to particles as F_fluid
         call dStress%initialize(amr,name='dStress',ncomp=3,ng=fs%nover,interp=interp_none); call dStress%register()
         if (.not.amr%xper) then; dStress%lo_bc(1,:)=amrex_bc_foextrap; dStress%hi_bc(1,:)=amrex_bc_foextrap; end if
         if (.not.amr%yper) then; dStress%lo_bc(2,:)=amrex_bc_foextrap; dStress%hi_bc(2,:)=amrex_bc_foextrap; end if
      end block create_workspace

      ! Initialize regridding and build the initial coupled state
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging (flow features; body handled by apd's callback)
         fs%user_tagging=>my_tagger
         call param_read('Tagging Re',Re_tag)
         call param_read('Tagging Rho',Rho_tag)
         ! Create initial grid (quiescent: no flow features yet, body added next)
         call amr%init_from_scratch(time=time%t)
         ! Seed the body, deposit VF, then regrid to refine around it
         call seed_body()
         call apd%update_VF()
         call amr%regrid(baselvl=0,time=time%t)
         call apd%get_info()
         call pd%get_info()
         ! Build the initial bond network on the final AMR hierarchy
         call handoff()
         call apd%update_VF()
         ! Initial solid-velocity deposit (body at rest -> Usolid=0) and fluid VF
         call deposit_solid_velocity()
         call update_VFf()
         ! Initialize fluid primitives, face velocity, viscosities, SGS
         call fs%get_primitive(Q=fs%Q); call fs%get_face_velocity()
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize visualization
      create_viz: block
         ! Eulerian-field visualization
         call viz%initialize(amr=amr,name='amrcomp_devol',use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(apd%VF,1,'VF')
         call viz%add_scalar(Usolid,1,'Us')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         ! Particle visualization
         call pviz%initialize(apd,name='amrcomp_devol')
         call pviz%select_comp('flag',on=.true.)
         call pviz%select_comp('dil', on=.true.)
         call pviz%select_comp('damage',on=.true.)
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) then
            call viz%write(time=time%t)
            call pviz%write(time=time%t)
         end if
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call apd%get_info()
         call pd%get_info()
         call pd%get_cfl(dt=time%dt,cfl=time%cfl)
         call get_force()
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Qmin(1),'RHOmin')
         call mfile%add_column(fs%Qmax(1),'RHOmax')
         call mfile%add_column(Fib(1),'Cd')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLp,'CFLp')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Mass')
         call consfile%add_column(fs%Qint(2),'U Momentum')
         call consfile%add_column(fs%Qint(3),'V Momentum')
         call consfile%add_column(fs%Qint(4),'W Momentum')
         call consfile%add_column(fs%Qint(5),'Internal energy')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
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
         ! Create solid monitor
         pdfile=monitor(amRoot=amr%amRoot,name='solid')
         call pdfile%add_column(time%n,'Timestep')
         call pdfile%add_column(time%t,'Time')
         call pdfile%add_column(pd%np,'Particle count')
         call pdfile%add_column(pd%nb,'Bond count')
         call pdfile%add_column(pd%nb_broken,'Bonds broken')
         call pdfile%add_column(n_sub,'Subcycles')
         call pdfile%add_column(pd%CFLp,'CFLp')
         call pdfile%add_column(pd%CFLe,'CFLe')
         call pdfile%add_column(apd%CFLv,'CFLv')
         call pdfile%add_column(apd%Umin,'Umin')
         call pdfile%add_column(pd%Umax,'Umax')
         call pdfile%add_column(apd%Vmin,'Vmin')
         call pdfile%add_column(apd%Vmax,'Vmax')
         call pdfile%add_column(apd%Wmin,'Wmin')
         call pdfile%add_column(apd%Wmax,'Wmax')
         call pdfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      real(WP) :: cfl_fs,cfl_pd

      ! Perform time integration
      do while (.not.time%done())

         ! dt is set by the FLUID CFL; the stiff solid sub-cycles within it
         call fs%get_cfl(dt=time%dt,cfl=cfl_fs)
         time%cfl=cfl_fs
         call time%adjust_dt()
         call time%increment()

         ! Remember old conserved variables and face velocities
         call fs%Qold%copy(src=fs%Q)
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)

         ! Sub-iterations
         do while (time%it.le.time%itmax)

            ! Build midpoint state: Q^{mid}=0.5*(Q+Qold), U^{mid}=0.5*(U+Uold)
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Get primitive variables at midpoint
            call fs%get_primitive(Q=fs%Q)

            ! Advance Q using Q^{mid}
            call fs%get_dQdt(dQdt=dQdt)
            ! Isothermal degassing: inject gas mass+energy into open pores inside the pellet
            call add_gas_source()
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
            call fs%Q%average_down(); call fs%Q%fill(time%t)

            ! Compute face velocities and ensure C/F consistency
            call fs%get_face_velocity(); call fs%average_down_velocity()

            ! Increment both velocities with current pressure term, then force
            ! the fluid toward the solid (PD->fluid). Skipped if s2f decoupled.
            call fs%get_primitive(Q=fs%Q)
            if (couple_s2f) then
               call fs%add_pressure(scale=time%dt,phi=fs%P,mask=VFf)
               call apply_ib_forcing()
            else
               call fs%add_pressure(scale=time%dt,phi=fs%P)
            end if

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Pressure correction
            if (fs%use_projection) then
               ! Solve pressure Helmholtz equation
               call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
               call fs%prepare_psolver(dt=time%dt)
               call fs%psolver%solve(rhs=fs%div)

               ! Correct both velocities with new pressure increment, re-forcing
               ! the fluid toward the solid unless s2f is decoupled
               if (couple_s2f) then
                  call fs%add_pressure(scale=time%dt,mask=VFf)
                  call apply_ib_forcing()
               else
                  call fs%add_pressure(scale=time%dt)
               end if

               ! Average down and fill ghosts
               call fs%Q%average_down(); call fs%Q%fill(time=time%t)
               call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
            end if

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute primitive variables
         call fs%get_primitive(Q=fs%Q)

         ! Fluid->solid load: divergence of the fluid stress tensor -> F_fluid.
         ! Uses the current-grid viscosities from the previous step; get_force
         ! fills dStress and the Cd diagnostic.
         call get_force()
         if (couple_f2s) call get_fluid_force()

         ! Sub-cycle the solid over the fluid step with F_fluid held fixed.
         ! n_sub = ceil(single-step apd CFL / cflmax) keeps each sub-step stable.
         pd_subcycle: block
            real(WP) :: dt_sub
            integer :: i_sub
call exchange_solid()
            call pd%get_cfl(dt=time%dt,cfl=cfl_pd)
            n_sub=1
            if (cfl_pd.gt.time%cflmax) n_sub=ceiling(cfl_pd/time%cflmax)
            dt_sub=time%dt/real(n_sub,WP)
            do i_sub=1,n_sub
               call pd%advance(dt_sub)
            end do
            call pd%get_cfl(dt=dt_sub,cfl=cfl_pd)
            call exchange_solid()
            call apd%redistribute()
            call apd%update_VF()
         end block pd_subcycle

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Refresh solid-velocity deposit and fluid VF on the (possibly new) grid
         call deposit_solid_velocity()
         call update_VFf()

         ! Compute viscosities and add SGS models (fresh on the new grid)
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time%t)
            call pviz%write(time=time%t)
         end if

         ! Perform and output monitoring (Cd computed by get_force above)
         call fs%get_info()
         call apd%get_info()
         call pd%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call pdfile%write()

      end do

   contains

      !> Fluid->solid load: F_fluid = interp(div(sigma)), the fluid stress-tensor
      !> divergence interpolated to each particle (matches peridynamics_shock).
      !> dStress must have been filled by get_force first.
      subroutine get_fluid_force()
         use amrex_amr_module, only: amrex_mfiter
         type(amrex_mfiter) :: mfi
         type(part), dimension(:), pointer :: p
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pdS
         integer(I8) :: np_,n
         integer :: lvl
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pdS=>dStress%mf(lvl)%dataptr(mfi)
               call apd%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  p(n)%F_fluid(1)=apd%interp(lvl,p(n)%pos,pdS,1)
                  p(n)%F_fluid(2)=apd%interp(lvl,p(n)%pos,pdS,2)
                  p(n)%F_fluid(3)=apd%interp(lvl,p(n)%pos,pdS,3)
               end do
            end do
            call amr%mfiter_destroy(mfi)
         end do
      end subroutine get_fluid_force

      !> Apply IB forcing - drive the fluid in solid cells toward rho*Usolid, VF-weighted
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pVF,pUs
         real(WP), dimension(:,:,:,:), allocatable :: pQold
         real(WP) :: sum_VF,sum_VFQ1,sum_VFQ5,VFface
         integer :: i,j,k,lvl,ii,jj,kk
         ! Compressible IB scheme requires updated ghosts for Q
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Apply IB scheme in solid region
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pQ=>fs%Q%mf(lvl)%dataptr(mfi)
               pU=>fs%U%mf(lvl)%dataptr(mfi)
               pV=>fs%V%mf(lvl)%dataptr(mfi)
               pW=>fs%W%mf(lvl)%dataptr(mfi)
               pVF=>VFf%mf(lvl)%dataptr(mfi)
               pUs=>Usolid%mf(lvl)%dataptr(mfi)
               ! Get interior tilebox
               bx=mfi%tilebox()
               ! Create backup of Q
               allocate(pQold,source=pQ)
               ! Loop over tile interior
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip pure fluid cells
                  if (pVF(i,j,k,1).eq.1.0_WP) cycle
                  ! Blend momentum Q(2-4) toward rho*Usolid in the solid fraction
                  pQ(i,j,k,2)=pVF(i,j,k,1)*pQ(i,j,k,2)+(1.0_WP-pVF(i,j,k,1))*pQ(i,j,k,1)*pUs(i,j,k,1)
                  pQ(i,j,k,3)=pVF(i,j,k,1)*pQ(i,j,k,3)+(1.0_WP-pVF(i,j,k,1))*pQ(i,j,k,1)*pUs(i,j,k,2)
                  pQ(i,j,k,4)=pVF(i,j,k,1)*pQ(i,j,k,4)+(1.0_WP-pVF(i,j,k,1))*pQ(i,j,k,1)*pUs(i,j,k,3)
                  ! VF-weighted neighbor average for Q(1) and Q(5)
                  sum_VF=0.0_WP; sum_VFQ1=0.0_WP; sum_VFQ5=0.0_WP
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                     if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                     sum_VF  =sum_VF  +pVF(i+ii,j+jj,k+kk,1)
                     sum_VFQ1=sum_VFQ1+pVF(i+ii,j+jj,k+kk,1)*pQold(i+ii,j+jj,k+kk,1)
                     sum_VFQ5=sum_VFQ5+pVF(i+ii,j+jj,k+kk,1)*pQold(i+ii,j+jj,k+kk,5)
                  end do; end do; end do
                  if (sum_VF.gt.0.0_WP) then
                     pQ(i,j,k,1)=pVF(i,j,k,1)*pQold(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ1/sum_VF
                     pQ(i,j,k,5)=pVF(i,j,k,1)*pQold(i,j,k,5)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ5/sum_VF
                  end if
               end do; end do; end do
               ! Deallocate pQold
               deallocate(pQold)
               ! Force face velocities toward the face-interpolated solid velocity
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  VFface=0.5_WP*sum(pVF(i-1:i,j,k,1))
                  pU(i,j,k,1)=VFface*pU(i,j,k,1)+(1.0_WP-VFface)*0.5_WP*sum(pUs(i-1:i,j,k,1))
               end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  VFface=0.5_WP*sum(pVF(i,j-1:j,k,1))
                  pV(i,j,k,1)=VFface*pV(i,j,k,1)+(1.0_WP-VFface)*0.5_WP*sum(pUs(i,j-1:j,k,2))
               end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  VFface=0.5_WP*sum(pVF(i,j,k-1:k,1))
                  pW(i,j,k,1)=VFface*pW(i,j,k,1)+(1.0_WP-VFface)*0.5_WP*sum(pUs(i,j,k-1:k,3))
               end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
         end do
      end subroutine apply_ib_forcing

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solvers
      call fs%finalize()
      call apd%finalize()
      call pd%finalize()
      call dQdt%finalize()
      call Usolid%finalize()
      call VFf%finalize()
      call dStress%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize material
      call fluid%finalize()
      ! Finalize visualization
      call viz%finalize()
      call pviz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call pdfile%finalize()
   end subroutine simulation_final

   !> Refresh fluid volume fraction VFf = clip(1 - apd%VF, 0, 1) with ghosts filled
   subroutine update_VFf()
      call VFf%setval(1.0_WP)
      call VFf%subtract(apd%VF)
      call VFf%clip(0.0_WP,1.0_WP)
      call VFf%fill(time=time%t)
   end subroutine update_VFf

   !> Degassing source linked to SOLID SURFACE AREA: gas is generated where the
   !> solid is being converted, i.e. at solid-gas interfaces (pore walls, crack
   !> faces), not uniformly in pore volume. Rate per cell = src_mdot * |grad VFf|
   !> (coarea: int |grad VFf| dV = surface area, so this is level-consistent under
   !> AMR; src_mdot is now a per-unit-surface-area rate). Injected on the GAS side
   !> (VFf>threshold) so the IB doesn't erase it; inside the pellet (r<R_disk).
   !> Self-amplifying: a new crack -> new surface -> more local degassing -> drives
   !> progressive venting along the opening.
   subroutine add_gas_source()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pdQ,pVF
      integer :: lvl,i,j,k
      real(WP) :: x,y,gx,gy,gz,gradmag,dxi,dyi,dzi,rate
      do lvl=0,amr%clvl()
         dxi=0.5_WP/amr%dx(lvl); dyi=0.5_WP/amr%dy(lvl); dzi=0.5_WP/amr%dz(lvl)
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pdQ=>dQdt%mf(lvl)%dataptr(mfi)
            pVF=>VFf%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pVF(i,j,k,1).lt.src_VFthresh) cycle      ! gas side (IB-safe)
               x=amr%xlo+(real(i,WP)+0.5_WP)*amr%dx(lvl)
               y=amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl)
               if (x*x+y*y.gt.R_disk**2) cycle              ! inside the pellet
               ! Surface density |grad VFf| (central differences; gz=0 in 2D)
               gx=(pVF(i+1,j,k,1)-pVF(i-1,j,k,1))*dxi
               gy=(pVF(i,j+1,k,1)-pVF(i,j-1,k,1))*dyi
               gz=(pVF(i,j,k+1,1)-pVF(i,j,k-1,1))*dzi
               gradmag=sqrt(gx*gx+gy*gy+gz*gz)
               rate=src_mdot*gradmag
               pdQ(i,j,k,1)=pdQ(i,j,k,1)+rate
               pdQ(i,j,k,5)=pdQ(i,j,k,5)+rate*src_e
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine add_gas_source

   !> Deposit volume-weighted particle velocity onto Usolid, normalized by apd%VF*cell_vol
   subroutine deposit_solid_velocity()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUs,pVF
      integer(I8) :: np_,n
      integer :: lvl,ii,jj,kk,i,j,k,c
      real(WP) :: dxi,dyi,dzi,wx,wy,wz,Vp
      real(WP), parameter :: VFtiny=1.0e-12_WP
      call Usolid%setval(0.0_WP)
      ! Deposit the extensive solid momentum-volume (vel*dV); process_deposit then
      ! makes it intensive with the SAME coarse/fine reconciliation as VF.
      Vp=apd%dV
      do lvl=0,amr%clvl()
         dxi=1.0_WP/amr%dx(lvl); dyi=1.0_WP/amr%dy(lvl); dzi=1.0_WP/amr%dz(lvl)
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pUs=>Usolid%mf(lvl)%dataptr(mfi)
            call apd%get_particles(lvl,mfi,p,np_)
            do n=1,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               ii=floor((p(n)%pos(1)-amr%xlo)*dxi-0.5_WP); wx=(p(n)%pos(1)-amr%xlo)*dxi-0.5_WP-real(ii,WP)
               jj=floor((p(n)%pos(2)-amr%ylo)*dyi-0.5_WP); wy=(p(n)%pos(2)-amr%ylo)*dyi-0.5_WP-real(jj,WP)
               kk=floor((p(n)%pos(3)-amr%zlo)*dzi-0.5_WP); wz=(p(n)%pos(3)-amr%zlo)*dzi-0.5_WP-real(kk,WP)
               do c=1,3
                  pUs(ii:ii+1,jj:jj+1,kk:kk+1,c)=pUs(ii:ii+1,jj:jj+1,kk:kk+1,c)+Vp*p(n)%vel(c)*reshape([ &
                  &  (1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz), &
                  &  (1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])
               end do
            end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      ! Post-process exactly like VF: extensive->intensive with coarse/fine
      ! reconciliation, fill ghosts, then the same smoothing filter.
      call apd%process_deposit(Usolid)
      call Usolid%fill(time=time%t)
      call apd%filter(Usolid)
      ! Normalize the (filtered) momentum density by the (filtered) VF -> velocity
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pUs=>Usolid%mf(lvl)%dataptr(mfi)
            pVF=>apd%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pVF(i,j,k,1).gt.VFtiny) then
                  pUs(i,j,k,1:3)=pUs(i,j,k,1:3)/pVF(i,j,k,1)
               else
                  pUs(i,j,k,1:3)=0.0_WP
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call Usolid%average_down(); call Usolid%fill(time=time%t)
   end subroutine deposit_solid_velocity

   !> Compute force on the solid from divergence of stress tensor inside the IB
   subroutine get_force()
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab
      use mpi_f08,   only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,fbx
      type(amrex_multifab) :: Sx,Sy,Sz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW,pVisc,pBeta,pP,pVF,pdS
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSx,pSy,pSz
      real(WP) :: dxi,dyi,dzi,vol,mu_f,beta_f,div
      real(WP), dimension(3,3) :: gradU
      ! Zero out force
      Fib=0.0_WP
      ! Work at finest level only
      lvl=amr%clvl()
      ! Get grid spacing
      dxi=1.0_WP/amr%dx(lvl)
      dyi=1.0_WP/amr%dy(lvl)
      dzi=1.0_WP/amr%dz(lvl)
      ! Build face-centered stress MultiFabs (3 force components each)
      call amr%mfab_build(lvl,Sx,ncomp=3,nover=0,atface=[.true. ,.false.,.false.]); call Sx%setval(0.0_WP)
      call amr%mfab_build(lvl,Sy,ncomp=3,nover=0,atface=[.false.,.true. ,.false.]); call Sy%setval(0.0_WP)
      call amr%mfab_build(lvl,Sz,ncomp=3,nover=0,atface=[.false.,.false.,.true. ]); call Sz%setval(0.0_WP)
      ! Fill face-centered stress fluxes
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pUVW =>fs%UVW%mf(lvl)%dataptr(mfi)
         pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
         pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
         pP   =>fs%P%mf(lvl)%dataptr(mfi)
         pSx  =>Sx%dataptr(mfi)
         pSy  =>Sy%dataptr(mfi)
         pSz  =>Sz%dataptr(mfi)
         ! X-face stresses
         fbx=mfi%nodaltilebox(1)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
            ! Get velocity gradient
            gradU(1,1)=dxi*(pUVW(i,j,k,1)-pUVW(i-1,j,k,1))
            gradU(2,1)=0.25_WP*dyi*(pUVW(i-1,j+1,k,1)-pUVW(i-1,j-1,k,1)+pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            gradU(3,1)=0.25_WP*dzi*(pUVW(i-1,j,k+1,1)-pUVW(i-1,j,k-1,1)+pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
            gradU(1,2)=dxi*(pUVW(i,j,k,2)-pUVW(i-1,j,k,2))
            gradU(2,2)=0.25_WP*dyi*(pUVW(i-1,j+1,k,2)-pUVW(i-1,j-1,k,2)+pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
            gradU(3,2)=0.25_WP*dzi*(pUVW(i-1,j,k+1,2)-pUVW(i-1,j,k-1,2)+pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            gradU(1,3)=dxi*(pUVW(i,j,k,3)-pUVW(i-1,j,k,3))
            gradU(2,3)=0.25_WP*dyi*(pUVW(i-1,j+1,k,3)-pUVW(i-1,j-1,k,3)+pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
            gradU(3,3)=0.25_WP*dzi*(pUVW(i-1,j,k+1,3)-pUVW(i-1,j,k-1,3)+pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            div=gradU(1,1)+gradU(2,2)+gradU(3,3)
            ! Get face viscosities
            mu_f=0.5_WP*sum(pVisc(i-1:i,j,k,1)); beta_f=0.5_WP*sum(pBeta(i-1:i,j,k,1))
            ! Compute face stresses
            pSx(i,j,k,1)=-0.5_WP*sum(pP(i-1:i,j,k,1))+mu_f*2.0_WP*gradU(1,1)+(beta_f-2.0_WP/3.0_WP*mu_f)*div
            pSx(i,j,k,2)=mu_f*(gradU(2,1)+gradU(1,2))
            pSx(i,j,k,3)=mu_f*(gradU(3,1)+gradU(1,3))
         end do; end do; end do
         ! Y-face stresses
         fbx=mfi%nodaltilebox(2)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
            ! Get velocity gradient
            gradU(1,1)=0.25_WP*dxi*(pUVW(i+1,j-1,k,1)-pUVW(i-1,j-1,k,1)+pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
            gradU(2,1)=dyi*(pUVW(i,j,k,1)-pUVW(i,j-1,k,1))
            gradU(3,1)=0.25_WP*dzi*(pUVW(i,j-1,k+1,1)-pUVW(i,j-1,k-1,1)+pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))
            gradU(1,2)=0.25_WP*dxi*(pUVW(i+1,j-1,k,2)-pUVW(i-1,j-1,k,2)+pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
            gradU(2,2)=dyi*(pUVW(i,j,k,2)-pUVW(i,j-1,k,2))
            gradU(3,2)=0.25_WP*dzi*(pUVW(i,j-1,k+1,2)-pUVW(i,j-1,k-1,2)+pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            gradU(1,3)=0.25_WP*dxi*(pUVW(i+1,j-1,k,3)-pUVW(i-1,j-1,k,3)+pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            gradU(2,3)=dyi*(pUVW(i,j,k,3)-pUVW(i,j-1,k,3))
            gradU(3,3)=0.25_WP*dzi*(pUVW(i,j-1,k+1,3)-pUVW(i,j-1,k-1,3)+pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            div=gradU(1,1)+gradU(2,2)+gradU(3,3)
            ! Get face viscosities
            mu_f=0.5_WP*sum(pVisc(i,j-1:j,k,1)); beta_f=0.5_WP*sum(pBeta(i,j-1:j,k,1))
            ! Compute face stresses
            pSy(i,j,k,1)=mu_f*(gradU(1,2)+gradU(2,1))
            pSy(i,j,k,2)=-0.5_WP*sum(pP(i,j-1:j,k,1))+mu_f*2.0_WP*gradU(2,2)+(beta_f-2.0_WP/3.0_WP*mu_f)*div
            pSy(i,j,k,3)=mu_f*(gradU(3,2)+gradU(2,3))
         end do; end do; end do
         ! Z-face stresses
         fbx=mfi%nodaltilebox(3)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
            ! Get velocity gradient
            gradU(1,1)=0.25_WP*dxi*(pUVW(i+1,j,k-1,1)-pUVW(i-1,j,k-1,1)+pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))
            gradU(2,1)=0.25_WP*dyi*(pUVW(i,j+1,k-1,1)-pUVW(i,j-1,k-1,1)+pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            gradU(3,1)=dzi*(pUVW(i,j,k,1)-pUVW(i,j,k-1,1))
            gradU(1,2)=0.25_WP*dxi*(pUVW(i+1,j,k-1,2)-pUVW(i-1,j,k-1,2)+pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))
            gradU(2,2)=0.25_WP*dyi*(pUVW(i,j+1,k-1,2)-pUVW(i,j-1,k-1,2)+pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))
            gradU(3,2)=dzi*(pUVW(i,j,k,2)-pUVW(i,j,k-1,2))
            gradU(1,3)=0.25_WP*dxi*(pUVW(i+1,j,k-1,3)-pUVW(i-1,j,k-1,3)+pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            gradU(2,3)=0.25_WP*dyi*(pUVW(i,j+1,k-1,3)-pUVW(i,j-1,k-1,3)+pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))
            gradU(3,3)=dzi*(pUVW(i,j,k,3)-pUVW(i,j,k-1,3))
            div=gradU(1,1)+gradU(2,2)+gradU(3,3)
            ! Get face viscosities
            mu_f=0.5_WP*sum(pVisc(i,j,k-1:k,1)); beta_f=0.5_WP*sum(pBeta(i,j,k-1:k,1))
            ! Compute face stresses
            pSz(i,j,k,1)=mu_f*(gradU(1,3)+gradU(3,1))
            pSz(i,j,k,2)=mu_f*(gradU(2,3)+gradU(3,2))
            pSz(i,j,k,3)=-0.5_WP*sum(pP(i,j,k-1:k,1))+mu_f*2.0_WP*gradU(3,3)+(beta_f-2.0_WP/3.0_WP*mu_f)*div
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      ! Take divergence of the stress tensor: store the per-cell field (used as
      ! the F_fluid load on the particles) and accumulate the (1-VF)-weighted
      ! net body force for the Cd diagnostic.
      call dStress%setval(0.0_WP)
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pSx=>Sx%dataptr(mfi)
         pSy=>Sy%dataptr(mfi)
         pSz=>Sz%dataptr(mfi)
         pVF=>VFf%mf(lvl)%dataptr(mfi)
         pdS=>dStress%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Cell-centered divergence of the stress tensor (force per volume)
            pdS(i,j,k,1)=dxi*(pSx(i+1,j,k,1)-pSx(i,j,k,1))+dyi*(pSy(i,j+1,k,1)-pSy(i,j,k,1))+dzi*(pSz(i,j,k+1,1)-pSz(i,j,k,1))
            pdS(i,j,k,2)=dxi*(pSx(i+1,j,k,2)-pSx(i,j,k,2))+dyi*(pSy(i,j+1,k,2)-pSy(i,j,k,2))+dzi*(pSz(i,j,k+1,2)-pSz(i,j,k,2))
            pdS(i,j,k,3)=dxi*(pSx(i+1,j,k,3)-pSx(i,j,k,3))+dyi*(pSy(i,j+1,k,3)-pSy(i,j,k,3))+dzi*(pSz(i,j,k+1,3)-pSz(i,j,k,3))
            ! Net body force over the solid fraction (Cd diagnostic)
            if (pVF(i,j,k,1).ge.1.0_WP) cycle
            vol=(1.0_WP-pVF(i,j,k,1))*amr%cell_vol(lvl)
            Fib(1)=Fib(1)+pdS(i,j,k,1)*vol
            Fib(2)=Fib(2)+pdS(i,j,k,2)*vol
            Fib(3)=Fib(3)+pdS(i,j,k,3)*vol
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      ! Fill ghosts of the stress-divergence field for the particle interpolation
      call dStress%fill(time=time%t)
      ! Allreduce force and normalize
      call MPI_ALLREDUCE(MPI_IN_PLACE,Fib,3,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      if (amr%nz.eq.1) then
         Fib=Fib/(p0*2.0_WP*R_disk*amr%dz(lvl))
      else
         Fib=Fib/(p0*Pi*R_disk**2)
      end if
      ! Cleanup
      call amr%mfab_destroy(Sx)
      call amr%mfab_destroy(Sy)
      call amr%mfab_destroy(Sz)
   end subroutine get_force


   !> Hand the post-bond_init amrpd state off to the graph
   !> pd. Extract owned nodes and bonds, repartition the nodes by Morton
   !> order of the reference configuration (balanced, motion-invariant, uses
   !> all ranks regardless of where the solid sits), build the solver, copy
   !> the material/contact configuration, and stamp each face particle's
   !> flag with its pd owner rank (flag = 7 + 8*owner) for sync routing.
   subroutine handoff()
      use amrex_amr_module, only: amrex_mfiter
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8), allocatable :: gids(:),rgid(:)
      real(WP), allocatable :: pos(:,:),vel(:,:),voll(:),rpos(:,:),rvel(:,:),rvol(:)
      integer, allocatable :: flags(:),owner(:),rflag(:)
      integer(I8) :: np_,n
      integer :: lvl,nn,i,nr
      ! Extract this rank's owned particles
      nn=0
      do lvl=0,amr%clvl()
         call apd%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call apd%get_particles(lvl,mfi,p,np_)
            nn=nn+int(np_)
         end do
         call apd%mfiter_destroy(mfi)
      end do
      allocate(gids(max(nn,1)),pos(3,max(nn,1)),vel(3,max(nn,1)),flags(max(nn,1)),voll(max(nn,1)),owner(max(nn,1)))
      i=0
      do lvl=0,amr%clvl()
         call apd%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call apd%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               i=i+1
               gids(i) =part_gid(p(n))
               pos(:,i)=p(n)%pos
               vel(:,i)=p(n)%vel
               flags(i)=p(n)%flag
               voll(i) =apd%dV
            end do
         end do
         call apd%mfiter_destroy(mfi)
      end do
      ! Balanced static partition of the reference configuration
      call pd_partition(nn,gids,pos,vel,flags,voll,owner,nr,rgid,rpos,rvel,rflag,rvol)
      ! Configure and build the pd
      call pd%initialize(name='pd',rho=apd%rho,elastic_modulus=apd%elastic_modulus, &
      &                    poisson_ratio=apd%poisson_ratio,delta=apd%delta,dV=apd%dV,    &
      &                    gravity=apd%gravity,                                        &
      &                    collapsed=[amr%nx.eq.1,amr%ny.eq.1,amr%nz.eq.1],           &
      &                    Ldom=[amr%xhi-amr%xlo,amr%yhi-amr%ylo,amr%zhi-amr%zlo],    &
      &                    per=[amr%xper,amr%yper,amr%zper])
      pd%s0=apd%s0
      ! Viscoplastic knobs (inert at their defaults; per-side e_v in the pd)
      pd%tau=apd%tau
      pd%visc_lambda=apd%visc_lambda
      pd%yield_stretch=apd%yield_stretch
      pd%sigma_yield=apd%sigma_yield
      ! Contact + domain config (soft-sphere contact active by default;
      ! open faces mean no walls, and exits are handled by dead-node muting)
      pd%use_contact=(apd%contact_dist.gt.0.0_WP)
      pd%contact_dist=apd%contact_dist
      pd%tau_col=apd%tau_col
      pd%e_n=apd%e_n; pd%e_w=apd%e_w; pd%clip_col=apd%clip_col
      pd%lo_bc=apd%lo_bc; pd%hi_bc=apd%hi_bc
      pd%dom_lo=[amr%xlo,amr%ylo,amr%zlo]
      pd%dom_hi=[amr%xhi,amr%yhi,amr%zhi]
      call pd%set_nodes(nr,rgid,rpos,rvel,rflag,rvol)
      ! Families detected natively from the reference configuration -- the
      ! amrpd bond container is never built in pd mode
      call pd%detect_families()
      ! Demote apd to grid-face duty: stamp owner tags (same walk order as the
      ! extraction above, so owner(i) lines up)
      i=0
      do lvl=0,amr%clvl()
         call apd%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call apd%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               i=i+1
               p(n)%flag=7+8*owner(i)
            end do
         end do
         call apd%mfiter_destroy(mfi)
      end do
      deallocate(gids,pos,vel,flags,voll,owner,rgid,rpos,rvel,rflag,rvol)
   end subroutine handoff

   !> Per-fluid-step solid exchange: push each face particle's
   !> interpolated F_fluid to its pd owner; pull back the owner's current
   !> (pos, vel, damage, alive) for write-back. Dead pd nodes turn their
   !> face particle into a PART_IS_DEAD tombstone (skipped by deposits and
   !> future syncs, matching amrpd's drop-on-exit semantics).
   subroutine exchange_solid()
      use amrex_amr_module, only: amrex_mfiter
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8), allocatable :: mgid(:)
      integer, allocatable :: mown(:)
      real(WP), allocatable :: mff(:,:),mpos(:,:),mvel(:,:),mdmg(:),malive(:)
      integer(I8) :: np_,n
      integer :: lvl,nm,i
      ! Count live face particles
      nm=0
      do lvl=0,amr%clvl()
         call apd%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call apd%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               nm=nm+1
            end do
         end do
         call apd%mfiter_destroy(mfi)
      end do
      allocate(mgid(max(nm,1)),mown(max(nm,1)),mff(3,max(nm,1)))
      allocate(mpos(3,max(nm,1)),mvel(3,max(nm,1)),mdmg(max(nm,1)),malive(max(nm,1)))
      ! Pack (gid, owner tag, F_fluid)
      i=0
      do lvl=0,amr%clvl()
         call apd%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call apd%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               i=i+1
               mgid(i)=part_gid(p(n))
               mown(i)=p(n)%flag/8
               mff(:,i)=p(n)%F_fluid
            end do
         end do
         call apd%mfiter_destroy(mfi)
      end do
      ! Collective round-trip with the pd
      call pd%exchange(nm,mgid,mown,mff,mpos,mvel,mdmg,malive)
      ! Write the pd state back into the grid face (same walk order)
      i=0
      do lvl=0,amr%clvl()
         call apd%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call apd%get_particles(lvl,mfi,p,np_)
            do n=1_I8,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               i=i+1
               p(n)%pos=mpos(:,i)
               p(n)%vel=mvel(:,i)
               p(n)%damage=mdmg(i)
               if (malive(i).lt.0.5_WP) p(n)%flag=PART_IS_DEAD
            end do
         end do
         call apd%mfiter_destroy(mfi)
      end do
      deallocate(mgid,mown,mff,mpos,mvel,mdmg,malive)
   end subroutine exchange_solid


end module simulation
