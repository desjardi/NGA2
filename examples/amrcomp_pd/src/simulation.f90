!> AMR compressible flow two-way coupled to a peridynamics solid (amrpd)
module simulation
   use precision,           only: WP,I8
   use amrgrid_class,       only: amrgrid
   use amrcomp_class,       only: amrcomp
   use amrpd_class,         only: amrpd,part,PART_MOVES,PART_INTEGRATES,PART_BONDS,PART_IS_DEAD,AMRPD_OPEN
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
   type(amrpd), target :: pd
   type(amrpdviz) :: pviz
   type(amrdata) :: Usolid          !< Solid velocity on the AMR mesh (3 comp)
   type(amrdata) :: VFf             !< Fluid volume fraction = 1 - pd%VF (amrcomp convention)
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

   !> Flow parameters
   real(WP) :: M2,Xs                  !< Post-shock Mach and shock location
   real(WP) :: Ms                     !< Shock Mach number
   real(WP) :: rho1,p1,u1             !< Pre-shock state
   real(WP) :: rho2,p2,u2             !< Post-shock state
   real(WP) :: Reynolds,Prandtl       !< Viscous parameters

   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP
   real(WP) :: L_spg=1.0_WP

   !> Adimensional force on IB
   real(WP), dimension(3) :: Fib

   !> Solid body initialization
   real(WP) :: beam_L=1.0_WP          !< Beam length
   real(WP) :: beam_t=0.125_WP        !< Beam thickness
   real(WP) :: beam_angle=75.0_WP     !< Beam tilt from the flow (x) axis [deg]
   real(WP) :: elem_size=0.0_WP       !< Solid element (particle) spacing dp

   !> PD sub-steps taken per fluid step (subcycling)
   integer :: n_sub=1

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

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
               r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (amr%nz.eq.1) r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2) ! Enable quasi-2D runs
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

   !> User init callback - set normal shock profile
   subroutine shock_init(solver,lvl,time,ba,dm)
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
      real(WP) :: rho,P,U,IE,H
      integer :: i
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         ! Get pointer to data
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Get tilebox with overlap
         bx=mfi%growntilebox(solver%nover)
         do i=bx%lo(1),bx%hi(1)
            ! Evaluate Heaviside function
            H=Hshock(x=Xs-(solver%amr%xlo+(real(i,WP)+0.5_WP)*solver%amr%dx(lvl)),delta=0.5_WP*solver%amr%dx(lvl))
            ! Interpolate between post-shock (H=1, right of shock) and pre-shock (H=0, left of shock)
            rho=rho1+(rho2-rho1)*H
            U=u1+(u2-u1)*H
            P=p1+(p2-p1)*H
            IE=fluid%get_e_from_p_rho(p=P,rho=rho,y=[1.0_WP])
            ! Set conserved variables
            pQ(i,:,:,1)=rho
            pQ(i,:,:,2)=rho*U
            pQ(i,:,:,3)=0.0_WP
            pQ(i,:,:,4)=0.0_WP
            pQ(i,:,:,5)=rho*IE
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shock_init

   !> Apply inflow BC at low-x (face=1)
   subroutine shock_dirichlet(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i,j,k
      select case (face)
       case (1)  ! X-LOW: Dirichlet inflow with post-shock values
         select case (comp)
          case ('U')  ! Staggered U=u2
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=u2
            end do; end do; end do
          case ('V','W')  ! Staggered V,W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')  ! Cell-centered Q=(rho2,rho2*u2,0,0,rho2*I2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=rho2
               p(i,j,k,2)=rho2*u2
               p(i,j,k,3)=0.0_WP
               p(i,j,k,4)=0.0_WP
               p(i,j,k,5)=rho2*fluid%get_e_from_p_rho(p=p2,rho=rho2,y=[1.0_WP])
            end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

   !> Tagger based on velocity and density laplacians. Refinement around the
   !> solid body is handled separately by amrpd's own VF-based tagging callback
   !> (registered in pd%initialize via pd%VF_tag).
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
            r_cyl=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
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

   !> Seed a thin beam of peridynamics particles, centered at the origin and
   !> tilted beam_angle from the flow (x) axis. 2D: single z-layer; 3D: extruded
   !> by the thickness in z (square cross-section).
   subroutine seed_body()
      use mathtools, only: Pi
      type(part), dimension(:), allocatable, target :: plist
      integer(I8) :: ntot,iflat
      integer :: i,j,k,ns,nt,khi
      real(WP) :: s,nn,zb,ct,st
      ! Lattice counts along the length (s) and across the thickness (n)
      ns=nint(beam_L/elem_size); nt=nint(beam_t/elem_size)
      khi=0; if (amr%nz.gt.1) khi=nt
      ct=cos(beam_angle*Pi/180.0_WP); st=sin(beam_angle*Pi/180.0_WP)
      if (amr%amRoot) then
         ntot=int(ns+1,I8)*int(nt+1,I8)*int(2*khi+1,I8)
         allocate(plist(ntot))
         iflat=0_I8
         do k=-khi,khi; do j=0,nt; do i=0,ns
            ! Beam-local coords, centered: s along axis, nn across thickness
            s =-0.5_WP*beam_L+real(i,WP)*elem_size
            nn=-0.5_WP*beam_t+real(j,WP)*elem_size
            zb=real(k,WP)*elem_size; if (amr%nz.eq.1) zb=0.0_WP
            iflat=iflat+1_I8
            ! Rotate (s,nn) by beam_angle about z (from the flow/x axis)
            plist(iflat)%pos    =[s*ct-nn*st, s*st+nn*ct, zb]
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
   end subroutine seed_body

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: A,B,C
         real(WP) :: Gamma,Cv
         ! EoS parameters
         call param_read('Gamma',Gamma)
         ! Shock parameters (input is M2, post-shock lab Mach)
         call param_read('Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rho2=1, u2=1, u1=0, T2=1
         rho2=1.0_WP
         p2=1.0_WP/(Gamma*M2**2)
         ! Quadratic for rho1: A*rho1^2 - B*rho1 + C = 0
         A=2.0_WP*Gamma*p2+(Gamma-1.0_WP)
         B=4.0_WP*Gamma*p2+(Gamma+1.0_WP)
         C=2.0_WP*Gamma*p2
         rho1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)  ! smaller root for compression
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rho1)
         u2=u1-1.0_WP
         p1=p2-rho1/(1.0_WP-rho1)
         if (p1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         ! Shock Mach number
         Ms=u1/sqrt(Gamma*p1/rho1)
         ! Shift to lab frame: pre-shock stationary
         u2=1.0_WP
         u1=0.0_WP
         ! Cv from T2=1
         Cv=p2/(rho2*(Gamma-1.0_WP))
         ! Build material
         call fluid%initialize(gamma=Gamma,pinf=0.0_WP,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Sutherland exponent',Suth_n,default=1.5_WP)
         call param_read('Sutherland temperature',Suth_T,default=0.4042_WP)
         ! Log shock conditions
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rho1=",es12.5," p1=",es12.5)') rho1,p1; call log(message)
         write(message,'("[Post-shock] rho2=",es12.5," p2=",es12.5)') rho2,p2; call log(message)
         call fluid%print()
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         real(WP) :: dp_slab
         amr%name='amrcomp_pd'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-05.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
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
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
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
         fs%user_init=>shock_init
         ! Set boundary conditions
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir; fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir; fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir; fs%W%hi_bc(1,1)=amrex_bc_foextrap
         fs%user_bc=>shock_dirichlet
      end block create_solver

      ! Initialize the peridynamics solid (containers + AMR callbacks)
      init_solid: block
         call pd%initialize(amr,name='pd')
         ! Material parameters
         call param_read('Material density',  pd%rho)
         call param_read('Elastic modulus',   pd%elastic_modulus)
         call param_read('Poisson ratio',     pd%poisson_ratio)
         call param_read('Critical energy',   pd%crit_energy,default=huge(1.0_WP))
         ! Solid spacing dp; horizon defaults to 3.0125*dp (Peridigm convention)
         call param_read('Element size',      elem_size)
         call param_read('Horizon',           pd%delta,default=3.0125_WP*elem_size)
         pd%search_radius=1.5_WP*pd%delta
         pd%dV=elem_size**3
         ! Beam geometry
         call param_read('Beam length',   beam_L,    default=1.0_WP)
         call param_read('Beam thickness',beam_t,    default=0.125_WP)
         call param_read('Beam angle',    beam_angle,default=75.0_WP)
         ! Gravity off by default (body driven by the flow)
         call param_read('Gravity',pd%gravity,default=[0.0_WP,0.0_WP,0.0_WP])
         ! Open domain BCs (y/z periodicity is handled by AMReX)
         pd%lo_bc=AMRPD_OPEN; pd%hi_bc=AMRPD_OPEN
         ! Refine the AMR mesh wherever the solid volume fraction exceeds VF_tag
         call param_read('Tagging VF',pd%VF_tag,default=0.1_WP)
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
         ! Fluid volume fraction = 1 - pd%VF (amrcomp's convention)
         call VFf%initialize(amr,name='VFf',ncomp=1,ng=fs%nover,interp=interp_none); call VFf%register()
         if (.not.amr%xper) then; VFf%lo_bc(1,1)=amrex_bc_foextrap; VFf%hi_bc(1,1)=amrex_bc_foextrap; end if
         ! Fluid stress-tensor divergence, interpolated to particles as F_fluid
         call dStress%initialize(amr,name='dStress',ncomp=3,ng=fs%nover,interp=interp_none); call dStress%register()
         if (.not.amr%xper) then; dStress%lo_bc(1,:)=amrex_bc_foextrap; dStress%hi_bc(1,:)=amrex_bc_foextrap; end if
      end block create_workspace

      ! Initialize regridding and build the initial coupled state
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging (flow features; body handled by pd's callback)
         fs%user_tagging=>my_tagger
         call param_read('Tagging Re',Re_tag)
         call param_read('Tagging Rho',Rho_tag)
         ! Create initial grid (body empty here, so only the shock refines)
         call amr%init_from_scratch(time=time%t)
         ! Seed the body, deposit VF, then regrid to refine around it
         call seed_body()
         call pd%update_VF()
         call amr%regrid(baselvl=0,time=time%t)
         call pd%get_info()
         ! Build the initial bond network on the final AMR hierarchy
         call pd%bond_init()
         call pd%update_VF()
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
         call viz%initialize(amr=amr,name='amrcomp_pd',use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(pd%VF,1,'VF')
         call viz%add_scalar(Usolid,1,'Us')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         ! Particle visualization
         call pviz%initialize(pd,name='amrcomp_pd')
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
         call pd%get_info()
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
         call pdfile%add_column(pd%Umin,'Umin')
         call pdfile%add_column(pd%Umax,'Umax')
         call pdfile%add_column(pd%Vmin,'Vmin')
         call pdfile%add_column(pd%Vmax,'Vmax')
         call pdfile%add_column(pd%Wmin,'Wmin')
         call pdfile%add_column(pd%Wmax,'Wmax')
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
         ! n_sub = ceil(single-step pd CFL / cflmax) keeps each sub-step stable.
         pd_subcycle: block
            real(WP) :: dt_sub
            integer :: i_sub
            call pd%get_cfl(dt=time%dt,cfl=cfl_pd)
            n_sub=1
            if (cfl_pd.gt.time%cflmax) n_sub=ceiling(cfl_pd/time%cflmax)
            dt_sub=time%dt/real(n_sub,WP)
            do i_sub=1,n_sub
               call pd%advance(dt_sub)
            end do
            ! Refresh monitored pd CFLs to the actual (stable) sub-step values
            call pd%get_cfl(dt=dt_sub,cfl=cfl_pd)
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
               call pd%get_particles(lvl,mfi,p,np_)
               do n=1,np_
                  if (p(n)%flag.eq.PART_IS_DEAD) cycle
                  p(n)%F_fluid(1)=pd%interp(lvl,p(n)%pos,pdS,1)
                  p(n)%F_fluid(2)=pd%interp(lvl,p(n)%pos,pdS,2)
                  p(n)%F_fluid(3)=pd%interp(lvl,p(n)%pos,pdS,3)
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

   !> Refresh fluid volume fraction VFf = clip(1 - pd%VF, 0, 1) with ghosts filled
   subroutine update_VFf()
      call VFf%setval(1.0_WP)
      call VFf%subtract(pd%VF)
      call VFf%clip(0.0_WP,1.0_WP)
      call VFf%fill(time=time%t)
   end subroutine update_VFf

   !> Deposit volume-weighted particle velocity onto Usolid, normalized by pd%VF*cell_vol
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
      Vp=pd%dV
      do lvl=0,amr%clvl()
         dxi=1.0_WP/amr%dx(lvl); dyi=1.0_WP/amr%dy(lvl); dzi=1.0_WP/amr%dz(lvl)
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pUs=>Usolid%mf(lvl)%dataptr(mfi)
            call pd%get_particles(lvl,mfi,p,np_)
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
      call pd%process_deposit(Usolid)
      call Usolid%fill(time=time%t)
      call pd%filter(Usolid)
      ! Normalize the (filtered) momentum density by the (filtered) VF -> velocity
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pUs=>Usolid%mf(lvl)%dataptr(mfi)
            pVF=>pd%VF%mf(lvl)%dataptr(mfi)
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
         Fib=Fib/(0.5_WP*rho2*u2**2*1.0_WP*amr%dz(lvl))
      else
         Fib=Fib/(0.5_WP*rho2*u2**2*Pi*0.25_WP)
      end if
      ! Cleanup
      call amr%mfab_destroy(Sx)
      call amr%mfab_destroy(Sy)
      call amr%mfab_destroy(Sz)
   end subroutine get_force

end module simulation
