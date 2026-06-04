!> AMR compressible drop test case
module simulation
   use precision,           only: WP
   use string,              only: str_medium
   use amrgrid_class,       only: amrgrid
   use amrmpcomp_class,     only: amrmpcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use amrio_class,         only: amrio
   use stiffened_gas_class, only: stiffened_gas
   use ideal_gas_class,     only: ideal_gas
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach
   
   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time
   
   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile
   
   !> Materials
   type(stiffened_gas), target :: water
   type(ideal_gas),     target :: air

   !> Flow parameters
   real(WP) :: rhoG1,pG1,u1           !< Pre-shock gas state
   real(WP) :: rhoG2,pG2,u2           !< Post-shock gas state
   real(WP) :: rhoL1,pL1              !< Initial liquid state
   real(WP) :: M2,Xs                  !< Post-shock Mach and shock location
   real(WP) :: Ms                     !< Shock Mach number
   real(WP) :: density_ratio          !< rhoL1/rhoG1
   real(WP) :: ML                     !< Liquid Mach number
   real(WP) :: Reynolds,visc_ratio    !< Viscosity 
   real(WP) :: Prandtl ,diff_ratio    !< Heat diffusivity
   real(WP) :: Weber                  !< Weber number
   
   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Moving wall model
   type(amrdata), target :: IBw
   real(WP) :: Xw,Uw

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP
   real(WP) :: L_spg=1.0_WP

   !> Tagging parameter
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for sphere (centered at origin)
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2) ! Enable quasi-2D runs
   end function sphere_levelset

   !> Generalized mechanical relaxation for stiffened-gas/ideal-gas pair.
   !> Solves quadratic for equilibrium pressure Peq where PL+Pjump=PG=Peq,
   !> then adjusts VF and internal energies via p*dV work exchange.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7).
   !> Enforces pressure jump provided in Pjump.
   !> Assumes gas is ideal (pinf_g = 0) — terms with pinf_g are dropped.
   subroutine P_relax_generalized(VF,Q,Pjump)
      use amrmpcomp_class, only: VFlo,VFhi
      implicit none
      real(WP), intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP), intent(in) :: Pjump
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      real(WP), parameter :: RHOGmin=1.0e-2_WP
      real(WP), parameter :: phist=1.0_WP,phi0=0.0_WP   !< Temporal weighting, phist=1 should yield best results
      ! Skip if any conserved quantity is non-positive (EOS undefined)
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Get phasic pressures
      PL=water%get_p_from_rho_e(rho=Q(1)/(       VF),e=Q(3)/Q(1),y=[1.0_WP])
      PG=air%get_p_from_rho_e  (rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*water%get_c_from_p_rho(p=PL,rho=Q(1)/(       VF),y=[1.0_WP])**2
      ZG=Q(2)/(1.0_WP-VF)*air%get_c_from_p_rho  (p=PG,rho=Q(2)/(1.0_WP-VF),y=[1.0_WP])**2
      cJ=ZL/(ZG+ZL)
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      n1=VF*phist
      n0=VF*(phi0*Pint-phist*cJ*pjump)+Q(3)
      d1=phist+1.0_WP/(water%gamma-1.0_WP)
      d0=phi0*Pint-phist*cJ*pjump+water%gamma/(water%gamma-1.0_WP)*water%pinf
      a=d1*(1.0_WP/(air%gamma-1.0_WP)+phist*VF)+n1*(-1.0_WP/(air%gamma-1.0_WP)-phist)
      b=d1*(-pjump/(air%gamma-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump))+n1*(pjump/(air%gamma-1.0_WP)-phi0*Pint+phist*cJ*pjump)+d0*(1.0_WP/(air%gamma-1.0_WP)+phist*VF)+n0*(-1.0_WP/(air%gamma-1.0_WP)-phist)
      d=d0*(-pjump/(air%gamma-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump))+n0*(pjump/(air%gamma-1.0_WP)-phi0*Pint+phist*cJ*pjump)
      ! Get equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.-water%pinf.or.Peq-Pjump.le.0.0_WP) return
      ! Get equilibrium volume fraction
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.VFlo.or.VFeq.gt.VFhi) return
      ! Adjust conserved quantities
      Q(3)=Q(3)-(phi0*Pint+phist*Peq)*(VFeq-VF)
      Q(4)=Q(4)+(phi0*Pint+phist*Peq)*(VFeq-VF)
      VF=VFeq

      ! ================ Second step: thermal relaxation ================
      !a=Q(1)*water%cv+Q(2)*air%cv
      !b=Q(1)*water%cv*water%gamma*water%pinf+Q(2)*air%cv*water%pinf-sum(Q(3:4))*(Q(1)*water%cv*(water%gamma-1.0_WP)+Q(2)*air%cv*(air%gamma-1.0_WP))
      !d=-sum(Q(3:4))*Q(2)*air%cv*(air%gamma-1.0_WP)*water%pinf
      ! Get equilibrium pressure
      !if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      !Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      !if (Peq.le.max(0.0_WP,-water%pinf)) return
      ! Get equilibrium volume fraction
      !VFeq=Q(1)*water%cv*(water%gamma-1.0_WP)*Peq/(Q(1)*water%cv*(water%gamma-1.0_WP)*Peq+Q(2)*air%cv*(air%gamma-1.0_WP)*(Peq+water%pinf))
      ! Clean up solution
      !if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-water%pinf); end if
      !if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,0.0_WP); end if
      ! Adjust conserved quantities
      !Q(3)=(       VFeq)*(Peq+water%gamma*water%pinf)/(water%gamma-1.0_WP)
      !Q(4)=(1.0_WP-VFeq)*Peq/(air%gamma-1.0_WP)
      !VF=VFeq

   end subroutine P_relax_generalized

   !> Compute viscosity: Sutherland for gas, VF-weighted blend with liquid
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pQ,pVisc,pBeta,pDiff,pRHOL,pRHOG
      real(WP) :: r_cyl,blend,nu_spg,mu_spg,mu_g,mu_l,k_g,k_l
      real(WP), parameter :: Tmax_visc=10.0_WP
      real(WP), parameter :: myeps=1.0e-15_WP
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP), parameter :: Cdiff=0.1_WP
      ! Get maximum allowable kinematic viscosity in the sponge at finest level
      nu_spg=max_cfl*amr%min_meshsize(amr%clvl())**2/(4.0_WP*time%dt)
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Get maximum allowable kinematic viscosity in the sponge at that level
         !nu_spg=max_cfl*min(amr%dx(lvl)**2,amr%dy(lvl)**2,amr%dz(lvl)**2)/(4.0_WP*time%dt)
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Gas viscosity from Sutherland
               mu_g=(1.0_WP+Suth_T)*min(pTG(i,j,k,1),Tmax_visc)**Suth_n/(Reynolds*(min(pTG(i,j,k,1),Tmax_visc)+Suth_T))
               ! Liquid viscosity from ratio
               mu_l=visc_ratio*Reynolds**(-1.0_WP)
               ! Mixture viscosity
               !pVisc(i,j,k,1)=pVF(i,j,k,1)*mu_l+(1.0_WP-pVF(i,j,k,1))*mu_g ! Arithmetic averaging
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps)) ! Harmonic averaging
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Gas heat diffusivity: k=Cv*Gamma*mu/Pr
               k_g=air%gamma*air%cv*mu_g/Prandtl
               ! Liquid heat diffusivity from ratio
               k_l=diff_ratio*air%gamma*air%cv/(Reynolds*Prandtl)
               ! Mixture diffusivity
               !pDiff(i,j,k,1)=pVF(i,j,k,1)*k_l+(1.0_WP-pVF(i,j,k,1))*k_g ! Arithmetic averaging
               pDiff(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(k_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(k_g,myeps)) ! Harmonic averaging
               ! Apply sponge layer viscosity
               r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (amr%nz.eq.1) r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2) ! Enable quasi-2D runs
               if (r_cyl.gt.R_spg) then
                  blend=min((r_cyl-R_spg)/L_spg,1.0_WP)**2
                  mu_spg=nu_spg/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),myeps)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),myeps))
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*mu_spg)
                  pDiff(i,j,k,1)=max(pDiff(i,j,k,1),Cdiff*blend*mu_spg)
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities
   
   !> User init callback - set Q and VF/barycenters for a drop at rest with a shock
   subroutine shockdrop_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,uG,H
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energy of liquid
      IEL=water%get_e_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP])
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF and barycenters from levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=myVF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Compute local gas state from shock profile
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
            rhoG=rhoG1+(rhoG2-rhoG1)*H
            pG  =pG1  +(pG2  -pG1  )*H
            uG  =u1   +(u2   -u1   )*H
            ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG, rho_mix*U, 0, 0)
            pQ(i,j,k,1)=(       myVF)*rhoL1
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*air%get_e_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
            pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*uG
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockdrop_init

   !> Apply inflow BC at low-x (face=1)
   subroutine shock_dirichlet(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i,j,k
      select case (face)
       case (1)  ! X-LOW: Dirichlet inflow with post-shock (gas only, no liquid)
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
               p(i,j,k,1)=0.0_WP                  ! No liquid
               p(i,j,k,2)=rhoG2                   ! Gas density
               p(i,j,k,3)=0.0_WP                  ! No liquid energy
               p(i,j,k,4)=rhoG2*air%get_e_from_p_rho(p=pG2,rho=rhoG2,y=[1.0_WP]) ! Gas internal energy
               p(i,j,k,5)=rhoG2*u2                ! X-momentum
               p(i,j,k,6)=0.0_WP
               p(i,j,k,7)=0.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

   !> Set wall IB
   subroutine set_IBw(data,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP) :: dx,dy,dz,xlo,xhi,xwall
      integer :: i,j,k
      xwall=Xw+Uw*time
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            xlo=data%amr%xlo+real(i  ,WP)*dx
            xhi=data%amr%xlo+real(i+1,WP)*dx
            if (xlo.ge.xwall) then
               pVF(i,j,k,1)=1.0_WP
            else if (xhi.le.xwall) then
               pVF(i,j,k,1)=0.0_WP
            else
               pVF(i,j,k,1)=(xhi-xwall)/dx
            end if
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine set_IBw

   !> Tagger based on velocity and density laplacians
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpcomp), intent(inout) :: solver
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
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      ! Recast tags
      tags=tags_ptr
      ! Compute tags
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! First compute radial location to compare with sponge
            r_cyl=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
            ! Get local densities and their inverse
            rho_cc=max(sum(pQ(i  ,j  ,k  ,1:2)),solver%rho_floor); irho_cc=1.0_WP/rho_cc
            rho_xp=max(sum(pQ(i+1,j  ,k  ,1:2)),solver%rho_floor); irho_xp=1.0_WP/rho_xp
            rho_xm=max(sum(pQ(i-1,j  ,k  ,1:2)),solver%rho_floor); irho_xm=1.0_WP/rho_xm
            rho_yp=max(sum(pQ(i  ,j+1,k  ,1:2)),solver%rho_floor); irho_yp=1.0_WP/rho_yp
            rho_ym=max(sum(pQ(i  ,j-1,k  ,1:2)),solver%rho_floor); irho_ym=1.0_WP/rho_ym
            rho_zp=max(sum(pQ(i  ,j  ,k+1,1:2)),solver%rho_floor); irho_zp=1.0_WP/rho_zp
            rho_zm=max(sum(pQ(i  ,j  ,k-1,1:2)),solver%rho_floor); irho_zm=1.0_WP/rho_zm
            ! Compute Laplacian of each velocity component
            lapU=(pQ(i+1,j,k,5)*irho_xp-2.0_WP*pQ(i,j,k,5)*irho_cc+pQ(i-1,j,k,5)*irho_xm)*dxi2+(pQ(i,j+1,k,5)*irho_yp-2.0_WP*pQ(i,j,k,5)*irho_cc+pQ(i,j-1,k,5)*irho_ym)*dyi2+(pQ(i,j,k+1,5)*irho_zp-2.0_WP*pQ(i,j,k,5)*irho_cc+pQ(i,j,k-1,5)*irho_zm)*dzi2
            lapV=(pQ(i+1,j,k,6)*irho_xp-2.0_WP*pQ(i,j,k,6)*irho_cc+pQ(i-1,j,k,6)*irho_xm)*dxi2+(pQ(i,j+1,k,6)*irho_yp-2.0_WP*pQ(i,j,k,6)*irho_cc+pQ(i,j-1,k,6)*irho_ym)*dyi2+(pQ(i,j,k+1,6)*irho_zp-2.0_WP*pQ(i,j,k,6)*irho_cc+pQ(i,j,k-1,6)*irho_zm)*dzi2
            lapW=(pQ(i+1,j,k,7)*irho_xp-2.0_WP*pQ(i,j,k,7)*irho_cc+pQ(i-1,j,k,7)*irho_xm)*dxi2+(pQ(i,j+1,k,7)*irho_yp-2.0_WP*pQ(i,j,k,7)*irho_cc+pQ(i,j-1,k,7)*irho_ym)*dyi2+(pQ(i,j,k+1,7)*irho_zp-2.0_WP*pQ(i,j,k,7)*irho_cc+pQ(i,j,k-1,7)*irho_zm)*dzi2
            ! Estimate sgs velocity from Laplacian
            u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
            ! Calculate cell Reynolds number and tag if too large
            Re=Reynolds*u_sgs*delta
            if (Re.gt.Re_tag.and.(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)) tagarr(i,j,k,1)=SETtag
            ! Compute normalized Laplacian of mixture density and tag if too large
            lapRHO=(rho_xp-2.0_WP*rho_cc+rho_xm)*dxi2+(rho_yp-2.0_WP*rho_cc+rho_ym)*dyi2+(rho_zp-2.0_WP*rho_cc+rho_zm)*dzi2
            avgRHO=(rho_cc+rho_xp+rho_xm+rho_yp+rho_ym+rho_zp+rho_zm)/7.0_WP
            lapRHO=abs(lapRHO)*delta2/avgRHO
            if (lapRHO.gt.Rho_tag.and.(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize AMR grid
      create_amrgrid: block
         ! Set name
         amr%name='amrcomp_drop'
         ! Read in base grid size
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         ! Set domain
         amr%xlo=-05.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         ! Set periodicity
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         ! Read in max level
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
         end if
         ! Initialize
         call amr%initialize()
      end block create_amrgrid

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: A,B,C
         real(WP) :: GammaL,PinfL,CvL
         real(WP) :: GammaG,CvG
         real(WP) :: T_G
         ! Gas EoS parameters (ideal gas)
         call param_read('GammaG',GammaG)
         ! Liquid EoS: gamma only, PinfL is computed below
         call param_read('GammaL',GammaL)
         ! Shock parameters (gas phase, uses GammaG)
         call param_read('Gas Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rhoG2=1, Deltau=1, T2=1
         rhoG2=1.0_WP
         pG2=1.0_WP/(GammaG*M2**2)
         ! Quadratic for rhoG1: A*rhoG1^2 - B*rhoG1 + C = 0
         A=2.0_WP*GammaG*pG2+(GammaG-1.0_WP)
         B=4.0_WP*GammaG*pG2+(GammaG+1.0_WP)
         C=2.0_WP*GammaG*pG2
         rhoG1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)  ! smaller root for compression
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rhoG1)
         u2=u1-1.0_WP
         pG1=pG2-rhoG1/(1.0_WP-rhoG1)
         if (pG1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         ! Shock Mach number
         Ms=u1/sqrt(GammaG*pG1/rhoG1)
         ! Shift to lab frame: pre-shock stationary
         u2=1.0_WP
         u1=0.0_WP
         ! CvG from T2=1
         CvG=pG2/(rhoG2*(GammaG-1.0_WP))
         ! Surface tension
         call param_read('Weber number',Weber)
         ! Liquid state from density ratio and liquid Mach number
         call param_read('Density ratio',density_ratio)
         call param_read('Liquid Mach number',ML)
         rhoL1=density_ratio
         pL1=pG1+4.0_WP/Weber                   ! Force pressure equilibrium, accounting for 3D Laplace pressure
         if (amr%nz.eq.1) pL1=pG1+2.0_WP/Weber  ! Force pressure equilibrium, accounting for 2D Laplace pressure
         PinfL=rhoL1/(GammaL*ML**2)-pL1
         ! Pre-shock gas temperature (ideal gas, T = p/((gamma-1)*Cv*rho))
         T_G=pG1/(rhoG1*(GammaG-1.0_WP)*CvG)
         CvL=(pL1+PinfL)/(rhoL1*(GammaL-1.0_WP)*T_G) ! Force thermal equilibrium
         ! Build materials
         call air%initialize  (gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP,name='air')
         call water%initialize(gamma=GammaL,pinf=PinfL,cv=CvL,q=0.0_WP,qp=0.0_WP,name='water')
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Viscosity ratio',visc_ratio)
         call param_read('Diffusivity ratio',diff_ratio)
         call param_read('Sutherland exponent',Suth_n)
         call param_read('Sutherland temperature',Suth_T)
         ! Log
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5)') rhoG1,pG1; call log(message)
         write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5)') rhoG2,pG2; call log(message)
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," ML=",es12.5)') rhoL1,pL1,ML; call log(message)
         call water%print(); call air%print()
         write(message,'("[Visc]   Re=",es12.5," mu*=",es12.5," Suth_n=",es12.5," Suth_T=",es12.5)') Reynolds,visc_ratio,Suth_n,Suth_T; call log(message)
         write(message,'("[Surface tension] We=",es12.5)') Weber; call log(message)
      end block init_eos_and_flow

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

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS
         use amrdata_class,    only: interp_face_lin
         ! Assign materials and create flow solver
         fs%liq=>water; fs%gas=>air; call fs%initialize(amr=amr,name='drop')
         ! Set surface tension coefficient
         fs%sigma=1.0_WP/Weber
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Provide pressure relaxation model
         fs%relax=>P_relax_generalized
         ! Set initial conditions
         fs%user_init=>shockdrop_init
         ! Set BCs
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_GAS
            fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%U%lo_bc(1,:)=amrex_bc_ext_dir; fs%U%hi_bc(1,:)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)=amrex_bc_ext_dir; fs%V%hi_bc(1,:)=amrex_bc_foextrap
            fs%W%lo_bc(1,:)=amrex_bc_ext_dir; fs%W%hi_bc(1,:)=amrex_bc_foextrap
            fs%user_bc=>shock_dirichlet
         end if
      end block create_solver

      ! Create IB data
      create_IB: block
         use amrdata_class, only: interp_reinit
         call param_read('Wall location',Xw,default=-10.0_WP)
         call param_read('Wall velocity',Uw,default=0.0_WP)
         call IBw%initialize(amr,name='IBw',ncomp=1,ng=fs%nover,interp=interp_reinit); call IBw%register()
         IBw%user_init=>set_IBw
      end block create_IB
      
      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=7,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Re',Re_tag)
         call param_read('Tagging Rho',Rho_tag)
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Fresh start
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC
            call fs%build_plic(time%t)
            call fs%build_subVF()
            ! Initialize primitive variables
            call fs%get_primitive(Q=fs%Q)
            ! Initialize face velocities
            call fs%get_face_velocity()
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
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
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'drop',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_surfmesh(fs%smesh,'plic')
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
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%dPmax,'dPmax')
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
         call cflfile%add_column(fs%CFLst ,'CFLst' )
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
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
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_sl,'sl_max')
         call tfile%add_column(fs%wtmin_sl,'sl_min')
         call tfile%add_column(fs%wtmax_fv,'fv_max')
         call tfile%add_column(fs%wtmin_fv,'fv_min')
         call tfile%add_column(fs%wtmax_div,'div_max')
         call tfile%add_column(fs%wtmin_div,'div_min')
         call tfile%add_column(fs%wtmax_plicnet,'plicnet_max')
         call tfile%add_column(fs%wtmin_plicnet,'plicnet_min')
         call tfile%add_column(fs%wtmax_polygon,'polygon_max')
         call tfile%add_column(fs%wtmin_polygon,'polygon_min')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%add_column(fs%nmixed_min,'mixed_min')
         call tfile%write()
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

         ! Update wall IB
         update_wall: block
            integer :: lvl
            do lvl=0,amr%clvl()
               call IBw%user_init(lvl,time%t,amr%ba(lvl),amr%dm(lvl))
            end do
         end block update_wall

         ! Remember old state
         call fs%store_old()
         
         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Get most up-to-date pressure
         call fs%apply_relax(time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities
         call fs%get_face_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt,mask=IBw)
         ! Add surface tension term
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         ! Apply IB forcing
         call apply_ib_forcing()
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Get most up-to-date pressure
         call fs%apply_relax(time=time%t)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities
         call fs%get_face_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=time%dt,mask=IBw)
         ! Add surface tension term
         call fs%add_surface_tension(scale=time%dt)
         ! Apply IB forcing
         call apply_ib_forcing()
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================================================================================

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/drop_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()
         
      end do

   contains

      !> Apply IB forcing - zero Q inside solid and apply quasi-Neumann
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pVF
         real(WP), dimension(:,:,:,:), allocatable :: pQold
         real(WP) :: sum_VF,sum_VFQ(4) 
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
               pVF=>IBw%mf(lvl)%dataptr(mfi)
               ! Get interior tilebox
               bx=mfi%tilebox()
               ! Create backup of Q
               allocate(pQold,source=pQ)
               ! Loop over tile interior
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip pure fluid cells
                  if (pVF(i,j,k,1).eq.1.0_WP) cycle
                  ! Scale Q(5-7) by VF and drive towards wall velocity
                  pQ(i,j,k,5)=pVF(i,j,k,1)*pQ(i,j,k,5)+(1.0_WP-pVF(i,j,k,1))*sum(pQ(i,j,k,1:2))*Uw
                  pQ(i,j,k,6)=pVF(i,j,k,1)*pQ(i,j,k,6)
                  pQ(i,j,k,7)=pVF(i,j,k,1)*pQ(i,j,k,7)
                  ! VF-weighted neighbor average for Q(1:4)
                  sum_VF=0.0_WP; sum_VFQ=0.0_WP
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                     if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                     sum_VF      =sum_VF      +pVF(i+ii,j+jj,k+kk,1)
                     sum_VFQ(1:4)=sum_VFQ(1:4)+pVF(i+ii,j+jj,k+kk,1)*pQold(i+ii,j+jj,k+kk,1:4)
                  end do; end do; end do
                  if (sum_VF.gt.0.0_WP) then
                     pQ(i,j,k,1:4)=pVF(i,j,k,1)*pQold(i,j,k,1:4)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ(1:4)/sum_VF
                  end if
               end do; end do; end do
               ! Deallocate pQold
               deallocate(pQold)
               ! Force face velocities
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pU(i,j,k,1)=0.5_WP*sum(pVF(i-1:i,j,k,1))*pU(i,j,k,1)+(1.0_WP-0.5_WP*sum(pVF(i-1:i,j,k,1)))*Uw
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
      call Umag%finalize()
      call Mach%finalize()
      call IBw%finalize()
      ! Finalize materials
      call water%finalize()
      call air%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
