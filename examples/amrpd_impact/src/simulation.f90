!> AMR compressible impact test case
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use amrgrid_class,     only: amrgrid
   use amrmpcomp_class,   only: amrmpcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use amrio_class,       only: amrio
   use nasg_class,        only: nasg
   use ideal_gas_class,   only: ideal_gas
   use relax_ig_nasg_class, only: relax_ig_nasg
   use amrpd_class,       only: amrpd,part,PART_MOVES,PART_INTEGRATES,PART_BONDS,PART_IS_DEAD,AMRPD_OPEN,AMRPD_WALL
   use amrpdviz_class,    only: amrpdviz
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> Peridynamics solid bodies (clamped wall + incoming projectile)
   type(amrpd), target :: pd
   type(amrpdviz) :: pviz
   real(WP) :: pd_elem                 !< particle spacing
   real(WP) :: wall_thick,wall_half    !< wall slab geometry
   real(WP) :: disk_R,disk_x,disk_vel  !< projectile geometry + impact speed
   integer  :: n_sub=1                 !< PD sub-steps per fluid step

   !> Two-way IB coupling workspaces (PD solid <-> multiphase fluid)
   type(amrdata) :: Usolid             !< Solid velocity deposited on the mesh (3 comp)
   type(amrdata) :: VFf                !< Fluid volume fraction = 1 - pd%VF
   type(amrdata) :: dStress            !< Divergence of the fluid stress tensor (3 comp) -> F_fluid
   real(WP), dimension(3) :: Fib       !< Net fluid force on the solid (diagnostic)
   logical :: couple_s2f=.true.        !< Solid->fluid (IB forcing of the flow)
   logical :: couple_f2s=.true.        !< Fluid->solid (F_fluid reaction on particles)
   
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
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,pdfile
   
   !> Materials
   type(nasg),      target :: water
   type(ideal_gas), target :: gas

   !> Relaxation model
   type(relax_ig_nasg), target :: relax_model

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

   !> Drop disturbances
   real(WP) :: dist_amp=0.005_WP      !< Disturbance amplitude
   real(WP) :: dist_num=64.0_WP       !< Disturbance wavenumber
   
   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Drop initial location
   real(WP) :: x_drop

   !> Wall BC type ('slip' or 'noslip')
   character(len=str_medium) :: wall_bc_type

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP
   real(WP) :: L_spg=1.0_WP

   !> Tagging parameter
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for drop (centered at x=x_drop)
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      !real(WP) :: G
      !G=0.5_WP-sqrt((xyz(1)-x_drop)**2+xyz(2)**2+xyz(3)**2)
      !if (amr%nz.eq.1) G=0.5_WP-sqrt((xyz(1)-x_drop)**2+xyz(2)**2) ! Enable quasi-2D runs
      real(WP) :: G,r,theta,r_perturbed
      ! Local angle in xy plane around drop center
      theta=atan2(xyz(2),xyz(1)-x_drop)
      ! Perturbed radius
      r_perturbed=0.5_WP*(1.0_WP+dist_amp*cos(real(dist_num,WP)*theta))
      ! Distance from drop center
      if (amr%nz.eq.1) then
         r=sqrt((xyz(1)-x_drop)**2+xyz(2)**2)
      else
         r=sqrt((xyz(1)-x_drop)**2+xyz(2)**2+xyz(3)**2)
      end if
      G=r_perturbed-r
   end function sphere_levelset

   !> Compute viscosity: Sutherland for gas, VF-weighted blend with liquid
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pQ,pVisc,pBeta,pDiffL,pDiffG,pRHOL,pRHOG
      real(WP) :: r_cyl,blend,nu_spg,mu_spg,mu_g,mu_l
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
            pDiffL=>fs%diffL%mf(lvl)%dataptr(mfi)
            pDiffG=>fs%diffG%mf(lvl)%dataptr(mfi)
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
               ! Phasic heat diffusivities: gas k=cp*mu/Pr, liquid from ratio
               pDiffG(i,j,k,1)=gas%cp*mu_g/Prandtl
               pDiffL(i,j,k,1)=diff_ratio*gas%cp/(Reynolds*Prandtl)
               ! Apply sponge layer viscosity
               r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (amr%nz.eq.1) r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2) ! Enable quasi-2D runs
               if (r_cyl.gt.R_spg) then
                  blend=min((r_cyl-R_spg)/L_spg,1.0_WP)**2
                  mu_spg=nu_spg/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),myeps)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),myeps))
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*mu_spg)
                  pDiffL(i,j,k,1)=max(pDiffL(i,j,k,1),Cdiff*blend*mu_spg)
                  pDiffG(i,j,k,1)=max(pDiffG(i,j,k,1),Cdiff*blend*mu_spg)
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback - set Q and VF/barycenters for a drop moving into a static shock
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
            pQ(i,j,k,4)=pQ(i,j,k,2)*gas%get_e_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
            ! Launch the drop at the impact velocity (-x) into the at-rest gas; pure-gas cells
            ! keep their local velocity uG (=0 in the post-shock region). The t=0 slip at the
            ! drop surface is an accepted initial disequilibrium.
            if (myVF.ge.VFlo) then
               pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*(-disk_vel)
            else
               pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*uG
            end if
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockdrop_init

   !> BC routine: x-HIGH is an inflow with pre-shock state in the wall frame
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
       case (2)  ! X-HIGH: Dirichlet inflow with pre-shock state (gas only, no liquid)
         select case (comp)
          case ('U')  ! Staggered U=u1 (pre-shock, shifted to wall frame)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=u1
            end do; end do; end do
          case ('V','W')  ! Staggered V,W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')  ! Cell-centered Q in pre-shock gas
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP                  ! No liquid
               p(i,j,k,2)=rhoG1                   ! Gas density
               p(i,j,k,3)=0.0_WP                  ! No liquid energy
               p(i,j,k,4)=rhoG1*gas%get_e_from_p_rho(p=pG1,rho=rhoG1,y=[1.0_WP]) ! Gas internal energy
               p(i,j,k,5)=rhoG1*u1                ! X-momentum
               p(i,j,k,6)=0.0_WP
               p(i,j,k,7)=0.0_WP
            end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

   !> Tagger based on SGS Reynolds number, density/pressure errors, Ducros sensor
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      use amrtag,           only: lap_error,grd_error
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pPL,pVF,pUVW,pC
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,delta,delta2
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
      real(WP) :: lapU,lapV,lapW,u_sgs,Re
      real(WP) :: divu,vortx,vorty,vortz,vort,Ducros,Deps
      real(WP) :: r_cyl
      logical  :: in_zone
      integer :: i,j,k
      real(WP), parameter :: Reps=1.0e-2_WP
      real(WP), parameter :: Peps=1.0e-2_WP
      real(WP), parameter :: Cduc=0.05_WP
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx; dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy; dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz; dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      ! Recast tags
      tags=tags_ptr
      ! Compute tags
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pQ  =>solver%Q%mf(lvl)%dataptr(mfi)
         pPL =>solver%PL%mf(lvl)%dataptr(mfi)
         pVF =>solver%VF%mf(lvl)%dataptr(mfi)
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         pC  =>solver%C%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Refinement zone: away from sponge unless below maxlvl-1
            r_cyl=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
            in_zone=(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)

            ! Mixture density laplacian error
            rho_cc=sum(pQ(i  ,j,  k,  1:2))
            rho_xp=sum(pQ(i+1,j,  k,  1:2)); rho_xm=sum(pQ(i-1,j,  k,  1:2))
            rho_yp=sum(pQ(i,  j+1,k,  1:2)); rho_ym=sum(pQ(i,  j-1,k,  1:2))
            rho_zp=sum(pQ(i,  j,  k+1,1:2)); rho_zm=sum(pQ(i,  j,  k-1,1:2))
            if (lap_error(rho_cc,rho_xm,rho_xp,rho_ym,rho_yp,rho_zm,rho_zp,Reps).gt.Rho_tag.and.in_zone) tagarr(i,j,k,1)=SETtag

            ! Liquid pressure gradient
            if (pVF(i,j,k,1).gt.0.0_WP) then
               if (grd_error(pPL(i,j,k,1),pPL(i-1,j,k,1),pPL(i+1,j,k,1),pPL(i,j-1,k,1),pPL(i,j+1,k,1),pPL(i,j,k-1,1),pPL(i,j,k+1,1),Peps).gt.P_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
            end if

            ! SGS cell Reynolds number
            lapU=(pUVW(i+1,j,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i-1,j,k,1))*dxi2+(pUVW(i,j+1,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j-1,k,1))*dyi2+(pUVW(i,j,k+1,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j,k-1,1))*dzi2
            lapV=(pUVW(i+1,j,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i-1,j,k,2))*dxi2+(pUVW(i,j+1,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j-1,k,2))*dyi2+(pUVW(i,j,k+1,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j,k-1,2))*dzi2
            lapW=(pUVW(i+1,j,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i-1,j,k,3))*dxi2+(pUVW(i,j+1,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j-1,k,3))*dyi2+(pUVW(i,j,k+1,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j,k-1,3))*dzi2
            u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
            Re=Reynolds*u_sgs*delta
            if (Re.gt.Re_tag.and.in_zone) tagarr(i,j,k,1)=SETtag

            ! Ducros compression switch
            divu =0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))+0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))+0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            vortx=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))-0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            vorty=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))-0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            vortz=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))-0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            vort=sqrt(vortx**2+vorty**2+vortz**2)
            Deps=(Cduc*pC(i,j,k,1)/delta)**2
            Ducros=divu**2/max(divu**2+vort**2+Deps,tiny(1.0_WP))
            if (divu.lt.0.0_WP.and.Ducros.gt.Ducros_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Seed two PD bodies: a clamped wall slab at x in [0,wall_thick] and an
   !> incoming disk projectile centered at (disk_x,0) moving -x at disk_vel.
   !> Gap (disk_x-disk_R-wall_thick) must exceed the horizon so the two bodies
   !> do NOT bond at init -- they interact only through contact.
   subroutine seed_bodies()
      use precision, only: I8
      type(part), dimension(:), allocatable, target :: plist
      integer(I8) :: ntot,iflat
      integer :: i,j,k,nx,nyh,nr,khi,pass
      real(WP) :: x,y,zb
      nyh=nint(wall_half/pd_elem); nx=nint(wall_thick/pd_elem); nr=nint(disk_R/pd_elem)
      khi=0; if (amr%nz.gt.1) khi=nr
      if (amr%amRoot) then
         do pass=1,2
            iflat=0_I8
            ! Wall slab: x in [-wall_thick,0] (front at x=0, free); back + side layers
            ! clamped (a plate held in a rigid frame; front face free to crater)
            do k=-khi,khi; do j=-nyh,nyh; do i=0,nx
               x=-wall_thick+real(i,WP)*pd_elem; y=real(j,WP)*pd_elem
               iflat=iflat+1_I8
               if (pass.eq.2) then
                  zb=real(k,WP)*pd_elem; if (amr%nz.eq.1) zb=0.0_WP
                  call set_part(plist(iflat),[x,y,zb],[0.0_WP,0.0_WP,0.0_WP],i.eq.0.or.abs(j).eq.nyh)
               end if
            end do; end do; end do
            ! Disk projectile: COMMENTED OUT -- replaced by the liquid water drop (fluid phase,
            ! seeded via shockdrop_init / Drop location). Preserved for the copper stand-in; re-enable to revert.
            !do k=-khi,khi; do j=-nr,nr; do i=-nr,nr
            !   x=real(i,WP)*pd_elem; y=real(j,WP)*pd_elem
            !   if (x*x+y*y.gt.disk_R**2) cycle
            !   iflat=iflat+1_I8
            !   if (pass.eq.2) then
            !      zb=real(k,WP)*pd_elem; if (amr%nz.eq.1) zb=0.0_WP
            !      call set_part(plist(iflat),[disk_x+x,y,zb],[-disk_vel,0.0_WP,0.0_WP],.false.)
            !   end if
            !end do; end do; end do
            if (pass.eq.1) then; ntot=iflat; allocate(plist(ntot)); end if
         end do
      else
         ntot=0_I8; allocate(plist(0))
      end if
      call pd%append(plist,ntot)
      deallocate(plist)
   contains
      subroutine set_part(p,pos,vel,clamp)
         type(part), intent(inout) :: p
         real(WP), dimension(3), intent(in) :: pos,vel
         logical, intent(in) :: clamp
         p%pos=pos; p%vel=vel
         p%F_bond=0.0_WP; p%F_fluid=0.0_WP; p%mw=0.0_WP; p%dil=0.0_WP
         p%damage=0.0_WP; p%nb0=0.0_WP
         if (clamp) then
            p%flag=PART_BONDS                              ! anchored: in network, no motion
         else
            p%flag=PART_MOVES+PART_INTEGRATES+PART_BONDS   ! free
         end if
      end subroutine set_part
   end subroutine seed_bodies

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize AMR grid
      create_amrgrid: block
         ! Set name
         amr%name='impact'
         ! Read in base grid size
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         ! Set domain: the true (reflecting) wall BC is shifted to x<0 so the PD
         ! wall can live INSIDE the domain, occupying x in [-Wall thickness, 0]
         ! (front at x=0). The original [0, Domain length] region is preserved.
         call param_read('Wall thickness',wall_thick)
         call param_read('Domain length',amr%xhi)   ! xhi temporarily holds the length
         amr%xlo=-wall_thick
         amr%xhi=amr%xhi+amr%xlo                     ! SHIFT (not extend): length unchanged -> dx,nx unchanged
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         ! Set periodicity
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         ! Read in max level
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D: slab thickness = particle spacing so one z-layer of PD
         ! particles fills the cell (VF_bulk = elem/z_slab = 1). Matches amrcomp_pd.
         ! (The finest-fluid-cell convention would give VF = elem/finest < 1.)
         if (amr%nz.eq.1) then
            call param_read('Element size',pd_elem)
            amr%zlo=-0.5_WP*pd_elem
            amr%zhi=+0.5_WP*pd_elem
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
         real(WP) :: GammaL,PinfL,bL,CvL,qpL
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
         ! Galilean shift to wall frame (wall fixed): subtract 1 from both velocities
         ! After shift: post-shock fluid at u2=0 (matches wall), pre-shock fluid at u1=-1 (toward wall)
         u2=u2-1.0_WP
         u1=u1-1.0_WP
         ! Drop initial location
         call param_read('Drop location',x_drop)
         ! CvG from T2=1
         CvG=pG2/(rhoG2*(GammaG-1.0_WP))
         ! Surface tension
         call param_read('Weber number',Weber)
         ! Liquid EoS, fit to this case's reference parameters
         call param_read('Liquid pinf',PinfL)
         call param_read('Liquid covolume',bL)
         call param_read('Liquid cv',CvL)
         call param_read('Liquid qp',qpL)
         ! Pre-shock gas temperature (ideal gas, T = p/((gamma-1)*Cv*rho))
         T_G=pG1/(rhoG1*(GammaG-1.0_WP)*CvG)
         ! Pressure equilibrium (Laplace jump): liquid pressure = gas + surface tension
         pL1=pG1+4.0_WP/Weber                   ! 3D Laplace pressure
         if (amr%nz.eq.1) pL1=pG1+2.0_WP/Weber  ! 2D Laplace pressure
         ! Liquid density from the NASG EOS at thermal+pressure equilibrium (T_L=T_G, p=pL1) [Option A]
         rhoL1=(pL1+PinfL)/((GammaL-1.0_WP)*CvL*T_G+bL*(pL1+PinfL))
         density_ratio=rhoL1/rhoG1                                        ! diagnostic (was an input under SG)
         ML=1.0_WP/sqrt(GammaL*(pL1+PinfL)/(rhoL1*(1.0_WP-bL*rhoL1)))     ! diagnostic liquid Mach (Deltau=1)
         ! Build materials: gas = ideal-gas air; liquid = NASG water
         call gas%initialize(gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP,name='gas')
         call water%initialize(gamma=GammaL,pinf=PinfL,b=bL,cv=CvL,q=0.0_WP,qp=qpL,name='water')
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
         write(message,'("[Temp]   TL=",es12.5," TG=",es12.5)') water%get_T_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP]),T_G; call log(message)
         call water%print(); call gas%print()
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
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_odd
         use amrmpcomp_class,  only: BC_GAS,BC_REFLECT
         use amrdata_class,    only: interp_face_lin
         use messager,         only: die
         ! Assign materials and create flow solver
         fs%liq=>water; fs%gas=>gas; call fs%initialize(amr=amr,name='impact')
         ! Set surface tension coefficient
         fs%sigma=1.0_WP/Weber
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Provide pressure relaxation model
         call relax_model%initialize(gas=gas,liq=water); fs%relax=>relax_model
         relax_model%model=1 ! 1=Prelax (mechanical only); 2=pT
         ! Set initial conditions
         fs%user_init=>shockdrop_init

         ! Gas inflow from x+
         !fs%hi_bc(1)=BC_GAS
         !fs%Q%hi_bc(1,:)=amrex_bc_ext_dir
         !fs%U%hi_bc(1,:)=amrex_bc_ext_dir
         !fs%V%hi_bc(1,:)=amrex_bc_ext_dir
         !fs%W%hi_bc(1,:)=amrex_bc_ext_dir
         !fs%user_bc=>shock_dirichlet

         ! Neumann at x+
         fs%hi_bc(1)=BC_REFLECT
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%hi_bc(1,:)=amrex_bc_foextrap
         fs%V%hi_bc(1,:)=amrex_bc_foextrap
         fs%W%hi_bc(1,:)=amrex_bc_foextrap

         ! Wall BC at x- (90 degree contact)
         fs%lo_bc(1)=BC_REFLECT
         fs%Q%lo_bc(1,:)=amrex_bc_foextrap ! Extrapolate everything then correct
         fs%Q%lo_bc(1,5)=amrex_bc_reflect_odd
         fs%U%lo_bc(1,:)=amrex_bc_reflect_odd
         ! Tangential momenta Q(6:7) and face velocities V, W at wall: slip vs no-slip
         call param_read('Wall BC',wall_bc_type,default='noslip')
         select case (trim(wall_bc_type))
         case ('noslip')
            fs%Q%lo_bc(1,6:7)=amrex_bc_reflect_odd
            fs%V%lo_bc(1,:)  =amrex_bc_reflect_odd
            fs%W%lo_bc(1,:)  =amrex_bc_reflect_odd
         case ('slip')
            fs%Q%lo_bc(1,6:7)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)  =amrex_bc_foextrap
            fs%W%lo_bc(1,:)  =amrex_bc_foextrap
         case default
            call die('[simulation_init] Unknown Wall BC type: must be slip or noslip')
         end select

      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class,    only: interp_none
         use amrex_amr_module, only: amrex_bc_foextrap
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1    ,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1    ,ng=0,interp=interp_none); call Mach%register()
         ! IB coupling fields (interp_none: recomputed each step, allocate-don't-fill on regrid).
         ! x is non-periodic (wall/open) -> zero-gradient extrapolation in ghosts; y,z periodic.
         call Usolid%initialize(amr,name='Usolid',ncomp=3,ng=fs%nover,interp=interp_none); call Usolid%register()
         call VFf%initialize(amr,name='VFf',ncomp=1,ng=fs%nover,interp=interp_none); call VFf%register()
         call dStress%initialize(amr,name='dStress',ncomp=3,ng=fs%nover,interp=interp_none); call dStress%register()
         Usolid%lo_bc(1,:)=amrex_bc_foextrap; Usolid%hi_bc(1,:)=amrex_bc_foextrap
         VFf%lo_bc(1,:)=amrex_bc_foextrap;    VFf%hi_bc(1,:)=amrex_bc_foextrap
         dStress%lo_bc(1,:)=amrex_bc_foextrap; dStress%hi_bc(1,:)=amrex_bc_foextrap
      end block create_workspace

      ! Initialize the peridynamics solid (registers VF-based AMR tagging)
      init_solid: block
         call pd%initialize(amr,name='pd')
         call param_read('Material density',  pd%rho)
         call param_read('Elastic modulus',   pd%elastic_modulus)
         call param_read('Poisson ratio',     pd%poisson_ratio)
         call param_read('Critical energy',   pd%crit_energy,  default=huge(1.0_WP))
         call param_read('Failure stretch',   pd%fail_stretch, default=huge(1.0_WP))
         call param_read('Relaxation time',   pd%tau,          default=huge(1.0_WP))
         call param_read('Relaxation fraction',pd%visc_lambda, default=1.0_WP)
         call param_read('Yield stretch',     pd%yield_stretch,default=0.0_WP)
         call param_read('Element size',      pd_elem)
         call param_read('Horizon',           pd%delta,        default=3.0125_WP*pd_elem)
         pd%search_radius=1.5_WP*pd%delta
         pd%dV=pd_elem**3
         ! Two-body geometry: clamped wall slab + incoming disk projectile
         call param_read('Wall thickness',  wall_thick)
         call param_read('Wall half-height',wall_half)
         call param_read('Disk radius',     disk_R)
         call param_read('Disk location',   disk_x)
         call param_read('Impact velocity', disk_vel)
         ! x-low is a rigid PD wall: backstop for the clamped wall slab (no particles
         ! dropped at the boundary, and a solid floor behind the anchor). Rest open.
         pd%lo_bc=AMRPD_OPEN; pd%hi_bc=AMRPD_OPEN
         pd%lo_bc(1)=AMRPD_WALL
         ! Refine the AMR mesh wherever the solid volume fraction exceeds VF_tag
         call param_read('Tagging VF',pd%VF_tag,default=0.1_WP)
         ! Two-way coupling switches (default on)
         call param_read('Couple solid to fluid',couple_s2f,default=.true.)
         call param_read('Couple fluid to solid',couple_f2s,default=.true.)
      end block init_solid

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tag Reynolds value',Re_tag)
         call param_read('Tag density error' ,Rho_tag)
         call param_read('Tag pressure error',P_tag)
         call param_read('Tag Ducros value',Ducros_tag)
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
         ! Seed the PD bodies, deposit VF, regrid to refine around them, build bonds
         call seed_bodies()
         call pd%update_VF()
         call amr%regrid(baselvl=0,time=time%t)
         call pd%get_info()
         call pd%bond_init()
         call pd%update_VF()
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Prime the IB coupling: deposit solid velocity, build fluid VF, and the
         ! initial fluid load (get_force needs the viscosities just computed above)
         call deposit_solid_velocity()
         call update_VFf()
         call get_force()
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
         call viz%initialize(amr,'impact',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_surfmesh(fs%smesh,'plic')
         call viz%add_scalar(pd%VF,1,'solidVF')
         call viz%add_scalar(VFf,1,'VFf')
         call viz%add_scalar(Usolid,1,'Us')
         ! Particle visualization
         call pviz%initialize(pd,name='impact')
         call pviz%select_comp('flag',on=.true.)
         call pviz%select_comp('damage',on=.true.)
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
         if (viz_evt%occurs()) call pviz%write(time=time%t)
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
         call mfile%add_column(fs%TLmin,'TLmin')
         call mfile%add_column(fs%TLmax,'TLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%TGmin,'TGmin')
         call mfile%add_column(fs%TGmax,'TGmax')
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
         call consfile%add_column(fs%VFint,'VFint')
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

      ! Solid (PD) monitor
      create_solid_monitor: block
         call pd%get_info()
         pdfile=monitor(amRoot=amr%amRoot,name='solid')
         call pdfile%add_column(time%n,'Timestep')
         call pdfile%add_column(time%t,'Time')
         call pdfile%add_column(pd%np,'Particle count')
         call pdfile%add_column(pd%nb,'Bond count')
         call pdfile%add_column(pd%nb_broken,'Bonds broken')
         call pdfile%add_column(n_sub,'Subcycles')
         call pdfile%add_column(pd%CFLe,'CFLe')
         call pdfile%add_column(pd%Umax,'Umax')
         call pdfile%add_column(Fib(1),'Ffluid_x')
         call pdfile%add_column(Fib(2),'Ffluid_y')
         call pdfile%write()
      end block create_solid_monitor

   end subroutine simulation_init
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time -- dt is set by the FLUID CFL; the stiff PD solid
         ! sub-cycles inside it (with F_fluid held fixed) further down.
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

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
         call fs%apply_relax(dt=0.5_WP*time%dt,time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities and ensure C/F consistency
         call fs%get_face_velocity(); call fs%average_down_velocity()
         ! Add pressure (masked to the fluid fraction), then force the fluid toward the solid
         if (couple_s2f) then
            call fs%add_phasic_pressure(scale=0.5_WP*time%dt,mask=VFf)
         else
            call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         end if
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         if (couple_s2f) call apply_ib_forcing()
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
         call fs%apply_relax(dt=time%dt,time=time%t)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities and ensure C/F consistency
         call fs%get_face_velocity(); call fs%average_down_velocity()
         ! Add pressure (masked to the fluid fraction), then force the fluid toward the solid
         if (couple_s2f) then
            call fs%add_phasic_pressure(scale=time%dt,mask=VFf)
         else
            call fs%add_phasic_pressure(scale=time%dt)
         end if
         call fs%add_surface_tension(scale=time%dt)
         if (couple_s2f) call apply_ib_forcing()
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================================================================================

         ! Fluid->solid load: divergence of the fluid stress tensor -> F_fluid.
         ! (get_force uses the current-grid viscosities from the previous step.)
         call get_force()
         if (couple_f2s) call get_fluid_force()

         ! Sub-cycle the PD solid over the fluid step with F_fluid held fixed.
         ! Contact handles the disk<->wall impact inside pd%advance.
         pd_subcycle: block
            real(WP) :: dt_sub,cfl_pd
            integer :: i_sub
            call pd%get_cfl(dt=time%dt,cfl=cfl_pd)
            n_sub=1
            if (cfl_pd.gt.time%cflmax) n_sub=ceiling(cfl_pd/time%cflmax)
            dt_sub=time%dt/real(n_sub,WP)
            do i_sub=1,n_sub; call pd%advance(dt_sub); end do
         end block pd_subcycle

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Refresh IB coupling fields on the (possibly new) grid for the next step
         call deposit_solid_velocity()
         call update_VFf()

         ! Compute viscosities and SGS models (fresh on the new grid)
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         if (viz_evt%occurs()) call pviz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/impact_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call pd%get_info()
         call mfile%write()
         call pdfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

   contains

      !> Solid->fluid IB forcing (multiphase Q-layout): drive the mixture momentum
      !> Q(5:7) toward rho_mix*Usolid over the solid fraction, VF-weighted neighbor
      !> average for the phasic Q(1:4), and blend face velocities toward Usolid.
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         use amrmpcomp_class,  only: VFlo,VFhi
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pVF,pUs,pVFL
         real(WP), dimension(:,:,:,:), allocatable :: pQold
         real(WP) :: sum_VF,sum_VFQ(4),rho_mix,VFface
         integer :: i,j,k,lvl,ii,jj,kk
         ! Compressible IB scheme requires updated ghosts for Q
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pQ =>fs%Q%mf(lvl)%dataptr(mfi)
               pU =>fs%U%mf(lvl)%dataptr(mfi)
               pV =>fs%V%mf(lvl)%dataptr(mfi)
               pW =>fs%W%mf(lvl)%dataptr(mfi)
               pVF=>VFf%mf(lvl)%dataptr(mfi)
               pUs=>Usolid%mf(lvl)%dataptr(mfi)
               pVFL=>fs%VF%mf(lvl)%dataptr(mfi)        ! liquid volume fraction (for the cleanup)
               bx=mfi%tilebox()
               allocate(pQold,source=pQ)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  if (pVF(i,j,k,1).eq.1.0_WP) cycle           ! pure fluid cell
                  ! VF-weighted neighbor average for the phasic mass+energy Q(1:4)
                  sum_VF=0.0_WP; sum_VFQ=0.0_WP
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                     if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                     sum_VF      =sum_VF      +pVF(i+ii,j+jj,k+kk,1)
                     sum_VFQ(1:4)=sum_VFQ(1:4)+pVF(i+ii,j+jj,k+kk,1)*pQold(i+ii,j+jj,k+kk,1:4)
                  end do; end do; end do
                  if (sum_VF.gt.0.0_WP) pQ(i,j,k,1:4)=pVF(i,j,k,1)*pQold(i,j,k,1:4)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ(1:4)/sum_VF
                  ! Light cleanup: keep the IB-extended Q strictly consistent with the liquid
                  ! VF so get_primitive/relax can't hit eL=Q(3)/rhoL=0/0 (or eG=Q(4)/rhoG=0/0).
                  if (pVFL(i,j,k,1).lt.VFlo) then
                     pQ(i,j,k,1)=0.0_WP; pQ(i,j,k,3)=0.0_WP    ! no liquid here
                  else if (pVFL(i,j,k,1).gt.VFhi) then
                     pQ(i,j,k,2)=0.0_WP; pQ(i,j,k,4)=0.0_WP    ! no gas here
                  end if
                  ! Mixture momentum toward rho_mix*Usolid in the solid fraction, using the
                  ! POST-average/cleanup density so the recovered velocity is exactly Usolid
                  ! in the solid (removes the old/new-density mismatch at the receding edge).
                  rho_mix=pQ(i,j,k,1)+pQ(i,j,k,2)
                  pQ(i,j,k,5)=pVF(i,j,k,1)*pQ(i,j,k,5)+(1.0_WP-pVF(i,j,k,1))*rho_mix*pUs(i,j,k,1)
                  pQ(i,j,k,6)=pVF(i,j,k,1)*pQ(i,j,k,6)+(1.0_WP-pVF(i,j,k,1))*rho_mix*pUs(i,j,k,2)
                  pQ(i,j,k,7)=pVF(i,j,k,1)*pQ(i,j,k,7)+(1.0_WP-pVF(i,j,k,1))*rho_mix*pUs(i,j,k,3)
               end do; end do; end do
               deallocate(pQold)
               ! Face velocities toward the face-interpolated solid velocity
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

      !> Fluid->solid load: F_fluid = interp(div(sigma)) at each particle.
      !> dStress must have been filled by get_force first.
      subroutine get_fluid_force()
         use amrex_amr_module, only: amrex_mfiter
         use precision,        only: I8
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
      call pd%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      call Usolid%finalize()
      call VFf%finalize()
      call dStress%finalize()
      ! Finalize materials
      call water%finalize()
      call gas%finalize()
      ! Finalize visualization
      call viz%finalize()
      call pviz%finalize()
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
      call pdfile%finalize()
   end subroutine simulation_final

   !> Refresh fluid volume fraction VFf = clip(1 - pd%VF, 0, 1) with ghosts filled.
   !> Additionally wall off the lateral gaps so the flow can't get around the finite-height
   !> target: x<0 AND |y|>wall_half (outside the slab footprint) -> VFf=0 (full solid). The
   !> footprint itself (|y|<=wall_half), including the crater, keeps its computed VFf, so the
   !> crater fills with fluid where copper is gone. Usolid=0 in the gaps -> static wall.
   subroutine update_VFf()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer :: lvl,i,j,k
      real(WP) :: xc,yc
      call VFf%setval(1.0_WP)
      call VFf%subtract(pd%VF)
      call VFf%clip(0.0_WP,1.0_WP)
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF=>VFf%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               xc=amr%xlo+(real(i,WP)+0.5_WP)*amr%dx(lvl)
               yc=amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl)
               if (xc.lt.0.0_WP.and.abs(yc).gt.wall_half) pVF(i,j,k,1)=0.0_WP
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call VFf%fill(time=time%t)
   end subroutine update_VFf

   !> Deposit the PD particle velocity onto the mesh (Usolid). Trilinear PIC
   !> deposit of vel*dV, then process_deposit reconciles coarse/fine exactly like
   !> VF, fill+filter, and normalize by pd%VF to recover an intensive velocity.
   subroutine deposit_solid_velocity()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use precision,        only: I8
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUs,pVF
      integer(I8) :: np_,n
      integer :: lvl,ii,jj,kk,i,j,k,c
      real(WP) :: dxi,dyi,dzi,wx,wy,wz,Vp
      real(WP), parameter :: VFtiny=1.0e-12_WP
      call Usolid%setval(0.0_WP)
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
      ! Post-process exactly like VF: extensive->intensive, C/F reconcile, fill, filter
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

   !> Fluid stress-tensor divergence inside the IB. Builds face-centered stresses
   !> sigma = -p*I + mu(grad u + grad u^T) + (beta-2/3 mu)(div u) I at the finest
   !> level (pressure = gas pressure PG; no liquid present), takes the cell-centered
   !> divergence into dStress (force/volume, interpolated to particles as F_fluid),
   !> and accumulates the net (1-VFf)-weighted force Fib for monitoring.
   subroutine get_force()
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab
      use mpi_f08,   only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,fbx
      type(amrex_multifab) :: Sx,Sy,Sz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pUVW,pVisc,pBeta,pP,pPL,pVFL,pVF,pdS
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pSx,pSy,pSz
      real(WP) :: dxi,dyi,dzi,vol,mu_f,beta_f,div
      real(WP), dimension(3,3) :: gradU
      ! Zero out net force
      Fib=0.0_WP
      ! Work at finest level only (the solid is fully refined there)
      lvl=amr%clvl()
      dxi=1.0_WP/amr%dx(lvl); dyi=1.0_WP/amr%dy(lvl); dzi=1.0_WP/amr%dz(lvl)

      ! Build face-centered stress MultiFabs (3 force components each)
      call amr%mfab_build(lvl,Sx,ncomp=3,nover=0,atface=[.true. ,.false.,.false.]); call Sx%setval(0.0_WP)
      call amr%mfab_build(lvl,Sy,ncomp=3,nover=0,atface=[.false.,.true. ,.false.]); call Sy%setval(0.0_WP)
      call amr%mfab_build(lvl,Sz,ncomp=3,nover=0,atface=[.false.,.false.,.true. ]); call Sz%setval(0.0_WP)
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pUVW =>fs%UVW%mf(lvl)%dataptr(mfi)
         pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
         pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
         pP   =>fs%PG%mf(lvl)%dataptr(mfi)   ! gas pressure
         pPL  =>fs%PL%mf(lvl)%dataptr(mfi)   ! liquid pressure
         pVFL =>fs%VF%mf(lvl)%dataptr(mfi)   ! liquid VF (for the mixture pressure VF*PL+(1-VF)*PG)
         pSx  =>Sx%dataptr(mfi)
         pSy  =>Sy%dataptr(mfi)
         pSz  =>Sz%dataptr(mfi)
         ! X-face stresses
         fbx=mfi%nodaltilebox(1)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
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
            mu_f=0.5_WP*sum(pVisc(i-1:i,j,k,1)); beta_f=0.5_WP*sum(pBeta(i-1:i,j,k,1))
            pSx(i,j,k,1)=-0.5_WP*sum(pVFL(i-1:i,j,k,1)*pPL(i-1:i,j,k,1)+(1.0_WP-pVFL(i-1:i,j,k,1))*pP(i-1:i,j,k,1))+mu_f*2.0_WP*gradU(1,1)+(beta_f-2.0_WP/3.0_WP*mu_f)*div
            pSx(i,j,k,2)=mu_f*(gradU(2,1)+gradU(1,2))
            pSx(i,j,k,3)=mu_f*(gradU(3,1)+gradU(1,3))
         end do; end do; end do
         ! Y-face stresses
         fbx=mfi%nodaltilebox(2)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
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
            mu_f=0.5_WP*sum(pVisc(i,j-1:j,k,1)); beta_f=0.5_WP*sum(pBeta(i,j-1:j,k,1))
            pSy(i,j,k,1)=mu_f*(gradU(1,2)+gradU(2,1))
            pSy(i,j,k,2)=-0.5_WP*sum(pVFL(i,j-1:j,k,1)*pPL(i,j-1:j,k,1)+(1.0_WP-pVFL(i,j-1:j,k,1))*pP(i,j-1:j,k,1))+mu_f*2.0_WP*gradU(2,2)+(beta_f-2.0_WP/3.0_WP*mu_f)*div
            pSy(i,j,k,3)=mu_f*(gradU(3,2)+gradU(2,3))
         end do; end do; end do
         ! Z-face stresses
         fbx=mfi%nodaltilebox(3)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
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
            mu_f=0.5_WP*sum(pVisc(i,j,k-1:k,1)); beta_f=0.5_WP*sum(pBeta(i,j,k-1:k,1))
            pSz(i,j,k,1)=mu_f*(gradU(1,3)+gradU(3,1))
            pSz(i,j,k,2)=mu_f*(gradU(2,3)+gradU(3,2))
            pSz(i,j,k,3)=-0.5_WP*sum(pVFL(i,j,k-1:k,1)*pPL(i,j,k-1:k,1)+(1.0_WP-pVFL(i,j,k-1:k,1))*pP(i,j,k-1:k,1))+mu_f*2.0_WP*gradU(3,3)+(beta_f-2.0_WP/3.0_WP*mu_f)*div
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      ! Cell-centered divergence of the stress tensor -> dStress (force/volume)
      call dStress%setval(0.0_WP)
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pSx=>Sx%dataptr(mfi)
         pSy=>Sy%dataptr(mfi)
         pSz=>Sz%dataptr(mfi)
         pVF=>VFf%mf(lvl)%dataptr(mfi)
         pdS=>dStress%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            pdS(i,j,k,1)=dxi*(pSx(i+1,j,k,1)-pSx(i,j,k,1))+dyi*(pSy(i,j+1,k,1)-pSy(i,j,k,1))+dzi*(pSz(i,j,k+1,1)-pSz(i,j,k,1))
            pdS(i,j,k,2)=dxi*(pSx(i+1,j,k,2)-pSx(i,j,k,2))+dyi*(pSy(i,j+1,k,2)-pSy(i,j,k,2))+dzi*(pSz(i,j,k+1,2)-pSz(i,j,k,2))
            pdS(i,j,k,3)=dxi*(pSx(i+1,j,k,3)-pSx(i,j,k,3))+dyi*(pSy(i,j+1,k,3)-pSy(i,j,k,3))+dzi*(pSz(i,j,k+1,3)-pSz(i,j,k,3))
            ! Net body force over the solid fraction (diagnostic)
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
      call MPI_ALLREDUCE(MPI_IN_PLACE,Fib,3,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call amr%mfab_destroy(Sx)
      call amr%mfab_destroy(Sy)
      call amr%mfab_destroy(Sz)
   end subroutine get_force

end module simulation
