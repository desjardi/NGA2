!> AMR compressible sphere test case with shock initialization
module simulation
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrcomp_class,     only: amrcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
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

   !> IBs
   type(amrdata), target :: VF
   
   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   
   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile
   
   !> Stiffened gas EOS parameters
   real(WP) :: Gamma,Pinf,Cv

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
   
contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)-0.5_WP
      if (amr%nz.eq.1) G=sqrt(xyz(1)**2+xyz(2)**2)-0.5_WP ! Enable 2D case
   end function sphere_levelset

   !> P=EOS(RHO,I) - Stiffened gas
   pure real(WP) function get_P(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_P=RHO*I*(Gamma-1.0_WP)-Gamma*Pinf
   end function get_P
   
   !> T=f(RHO,P)
   pure real(WP) function get_T(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_T=(P+Pinf)/(Cv*RHO*(Gamma-1.0_WP))
   end function get_T
   
   !> C=f(RHO,P) - Speed of sound
   pure real(WP) function get_C(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_C=sqrt(Gamma*(P+Pinf)/RHO)
   end function get_C

   !> I=EOS(RHO,P)
   pure real(WP) function get_I(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_I=(P+Gamma*Pinf)/(RHO*(Gamma-1.0_WP))
   end function get_I

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
               pDiff(i,j,k,1)=Gamma*Cv*pVisc(i,j,k,1)/Prandtl
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
            IE=get_I(rho,P)
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
       case (1)  ! X-LOW: Dirichlet inflow with pre-shock (stationary) values
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
               p(i,j,k,5)=rho2*get_I(rho2,p2)
            end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

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
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      integer, parameter :: nref=3
      real(WP), parameter :: VFlo=1.0e-12_WP
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=nref,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VF

   !> Tagger based on velocity and density laplacians
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
      real(WP) :: lapU,lapV,lapW,u_sgs,Re,lapRHO,avgRHO,r_cyl,dist
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
            ! Tag near sphere surface
            dist=sphere_levelset([solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz],time)
            if (abs(dist).lt.dx) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Post-regrid dispatcher for automatic VF filling
   subroutine postregrid_VF(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrdata), pointer :: this
      call c_f_pointer(ctx,this)
      call this%fill(time=time,lbase=lbase)
   end subroutine postregrid_VF

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
         ! EoS parameters
         call param_read('Gamma',Gamma)
         Pinf=0.0_WP
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
         write(message,'("[Cv=",es12.5,"]")') Cv; call log(message)
      end block init_eos_and_flow
      
      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='amrcomp_sphere'
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
         ! Create flow solver
         call param_read('Use projection',fs%use_projection)
         call fs%initialize(amr=amr)
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set pressure convergence
         fs%psolver%max_iter=20
         fs%psolver%tol_rel=1.0e-5_WP
         fs%psolver%verbose=2
         ! Provide thermodynamic model
         fs%getP=>get_P
         fs%getC=>get_C
         fs%getT=>get_T
         ! Set initial conditions
         fs%user_init=>shock_init
         ! Set boundary conditions
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir
         fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%hi_bc(1,1)=amrex_bc_foextrap
         fs%user_bc=>shock_dirichlet
      end block create_solver

      ! Create VF for IB
      create_VF: block
         use amrex_amr_module, only: amrex_bc_foextrap
         use amrdata_class, only: interp_const
         use iso_c_binding, only: c_loc
         ! Create VF field with constant interpolation
         call VF%initialize(amr,name='VF',ncomp=1,ng=fs%nover,interp=interp_const); call VF%register()
         call amr%add_postregrid(postregrid_VF,c_loc(VF))
         VF%user_init=>init_VF
         VF%lo_bc(1,1)=amrex_bc_foextrap
         VF%hi_bc(1,1)=amrex_bc_foextrap
      end block create_VF
      
      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace
      
      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Re',Re_tag)
         call param_read('Tagging Rho',Rho_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Initialize primitive variables and face velocity
         call fs%get_primitive(Q=fs%Q); call fs%get_face_velocity()
         ! Compute viscosities and add SGS models
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding
      
      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr=amr,name='sphere',use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(VF,1,'VF')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz
      
      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
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
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
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
      end block create_monitors

   end subroutine simulation_init
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
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

            ! Interpolate velocity to the faces
            call fs%get_face_velocity()

            ! Increment both velocities with current pressure term (only 1st order in time but stable)
            call fs%get_primitive(Q=fs%Q)
            call fs%add_pressure(scale=time%dt,phi=fs%P)

            ! Apply IB direct forcing
            call apply_ibm()

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Pressure correction
            if (fs%use_projection) then
               ! Solve pressure Helmholtz equation
               call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
               call fs%prepare_psolver(dt=time%dt)
               call fs%psolver%solve(rhs=fs%div)

               ! Kill correction inside IB
               !call fs%psolver%sol%multiply(src=VF)

               ! Correct both velocities with new pressure increment
               call fs%add_pressure(scale=time%dt)

               ! Add pressure increment
               call fs%P%add(src=fs%psolver%sol)

               ! Average down and fill ghosts
               call fs%Q%average_down(); call fs%Q%fill(time=time%t)
               call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
            end if

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute primitive variables
         call fs%get_primitive(Q=fs%Q)

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities and add SGS models
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%Q,srcY=fs%Q,srcZ=fs%Q,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         
      end do

   contains

      !> Apply IB forcing - zero Q inside solid
      subroutine apply_ibm()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pVF
         real(WP) :: sum_VF,sum_VFQ1,sum_VFQ5,myVF
         integer :: i,j,k,lvl,ii,jj,kk
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pQ=>fs%Q%mf(lvl)%dataptr(mfi)
               pU=>fs%U%mf(lvl)%dataptr(mfi)
               pV=>fs%V%mf(lvl)%dataptr(mfi)
               pW=>fs%W%mf(lvl)%dataptr(mfi)
               pVF=>VF%mf(lvl)%dataptr(mfi)
               ! Get interior tilebox
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip pure fluid cells
                  if (pVF(i,j,k,1).eq.1.0_WP) cycle
                  ! Scale Q(2-4) by VF
                  pQ(i,j,k,2)=pVF(i,j,k,1)*pQ(i,j,k,2)
                  pQ(i,j,k,3)=pVF(i,j,k,1)*pQ(i,j,k,3)
                  pQ(i,j,k,4)=pVF(i,j,k,1)*pQ(i,j,k,4)
                  ! VF-weighted neighbor average for Q(1) and Q(5)
                  sum_VF=0.0_WP; sum_VFQ1=0.0_WP; sum_VFQ5=0.0_WP
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                     if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                     sum_VF  =sum_VF  +pVF(i+ii,j+jj,k+kk,1)
                     sum_VFQ1=sum_VFQ1+pVF(i+ii,j+jj,k+kk,1)*pQ(i+ii,j+jj,k+kk,1)
                     sum_VFQ5=sum_VFQ5+pVF(i+ii,j+jj,k+kk,1)*pQ(i+ii,j+jj,k+kk,5)
                  end do; end do; end do
                  if (sum_VF.gt.0.0_WP) then
                     pQ(i,j,k,1)=pVF(i,j,k,1)*pQ(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ1/sum_VF
                     pQ(i,j,k,5)=pVF(i,j,k,1)*pQ(i,j,k,5)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ5/sum_VF
                  end if
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
      end subroutine apply_ibm

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
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
      call Mach%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final
   
end module simulation
