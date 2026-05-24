!> AMR subgrid-scale (SGS) model utilities
!> Provides shared routines for computing kinematic eddy viscosity
!> and artificial viscosity from velocity fields. Output is kinematic
!> (i.e., not multiplied by density); callers handle the rho-multiply.
module amrsgs
   use precision,     only: WP
   use amrgrid_class, only: amrgrid
   implicit none
   private

   ! Expose routines
   public :: get_vreman       !< Vreman eddy viscosity model
   public :: get_viscartif    !< Artificial viscosity

   ! Shared parameters
   real(WP), parameter :: default_Cs=0.17_WP
   real(WP), parameter :: default_max_cfl=0.5_WP
   integer,  parameter :: default_nfilter=2
   real(WP), parameter :: default_Cartif=5.0_WP

contains

   !> Compute Vreman kinematic eddy viscosity from cell- or face-centered velocity into passed visc
   !> Input:  U(Ucomp), V(Vcomp), W(Wcomp) are cell-centered velocity components
   !> Output: visc is cell-centered kinematic eddy viscosity (1 component)
   subroutine get_vreman(dt,visc,U,V,W,Ucomp,Vcomp,Wcomp,Cs,max_cfl,nfilter)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      use amrdata_class,    only: amrdata
      use messager,         only: die
      implicit none
      real(WP), intent(in) :: dt
      type(amrdata), intent(inout) :: visc         !< Output kinematic eddy viscosity (ncomp=1, ng>=1)
      type(amrdata), intent(in) :: U,V,W           !< Velocities
      integer, intent(in), optional :: Ucomp,Vcomp,Wcomp  !< Component indices (default=1)
      real(WP), intent(in), optional :: Cs
      real(WP), intent(in), optional :: max_cfl
      integer,  intent(in), optional :: nfilter
      ! Locals
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pVisc
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,Cmodel,Aij,Bij,CFLmax
      real(WP), dimension(1:3,1:3) :: gradU,betaij
      integer :: lvl,i,j,k,si,sj,nf,uc,vc,wc
      logical :: is_stag

      ! Check velocity nodal locations
      check_velocity: block
         logical, dimension(3) :: nU,nV,nW
         nU=U%nodal; nV=V%nodal; nW=W%nodal
         if (all(nU.eqv.[.true.,.false.,.false.]).and. &
         &   all(nV.eqv.[.false.,.true.,.false.]).and. &
         &   all(nW.eqv.[.false.,.false.,.true.])) then
            is_stag=.true.
         else if (.not.any(nU).and..not.any(nV).and..not.any(nW)) then
            is_stag=.false.
         else
            call die('[amrsgs::get_vreman] U/V/W must be staggered (face-centered) or collocated (cell-centered)')
         end if
      end block check_velocity

      ! Check visc
      if (visc%ncomp.ne.1) call die('[amrsgs::get_vreman] visc must have exactly 1 component')
      if (visc%ng.lt.1) call die('[amrsgs::get_vreman] visc must have at least 1 ghost cell')

      ! Component indices
      if (present(Ucomp)) then; uc=Ucomp; else; uc=1; end if
      if (present(Vcomp)) then; vc=Vcomp; else; vc=1; end if
      if (present(Wcomp)) then; wc=Wcomp; else; wc=1; end if

      ! Model constant: c=2.5*Cs**2
      if (present(Cs)) then; Cmodel=2.5_WP*Cs**2; else; Cmodel=2.5_WP*default_Cs**2; end if
      ! Clipping based on maximum CFL
      if (present(max_cfl)) then; CFLmax=max_cfl; else; CFLmax=default_max_cfl; end if
      ! Number of filtering passes
      if (present(nfilter)) then; nf=nfilter; else; nf=default_nfilter; end if

      ! Loop over levels
      do lvl=0,visc%amr%clvl()
         ! Grid spacings
         dx=visc%amr%dx(lvl); dxi=1.0_WP/dx
         dy=visc%amr%dy(lvl); dyi=1.0_WP/dy
         dz=visc%amr%dz(lvl); dzi=1.0_WP/dz
         ! Max visc from CFL
         max_visc=CFLmax*visc%amr%min_meshsize(lvl)**2/(4.0_WP*dt)
         ! Zero out visc
         call visc%mf(lvl)%setval(0.0_WP)
         ! Phase 1: Compute kinematic eddy viscosity
         call visc%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pU=>U%mf(lvl)%dataptr(mfi)
            pV=>V%mf(lvl)%dataptr(mfi)
            pW=>W%mf(lvl)%dataptr(mfi)
            pVisc=>visc%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Cell-centered velocity gradient tensor...
               if (is_stag) then
                  ! ... based on staggered velocity
                  gradU(1,1)=        dxi*   (pU(i+1,j,k,uc)      -pU(i,j,k,uc)        )
                  gradU(2,1)=0.25_WP*dyi*sum(pU(i:i+1,j:j+1,k,uc)-pU(i:i+1,j-1:j,k,uc))
                  gradU(3,1)=0.25_WP*dzi*sum(pU(i:i+1,j,k:k+1,uc)-pU(i:i+1,j,k-1:k,uc))
                  gradU(1,2)=0.25_WP*dxi*sum(pV(i:i+1,j:j+1,k,vc)-pV(i-1:i,j:j+1,k,vc))
                  gradU(2,2)=        dyi*   (pV(i,j+1,k,vc)      -pV(i,j,k,vc)        )
                  gradU(3,2)=0.25_WP*dzi*sum(pV(i,j:j+1,k:k+1,vc)-pV(i,j:j+1,k-1:k,vc))
                  gradU(1,3)=0.25_WP*dxi*sum(pW(i:i+1,j,k:k+1,wc)-pW(i-1:i,j,k:k+1,wc))
                  gradU(2,3)=0.25_WP*dyi*sum(pW(i,j:j+1,k:k+1,wc)-pW(i,j-1:j,k:k+1,wc))
                  gradU(3,3)=        dzi*   (pW(i,j,k+1,wc)      -pW(i,j,k,wc)        )
               else
                  ! ... based on collocated velocity
                  gradU(1,1)=0.5_WP*dxi*(pU(i+1,j,k,uc)-pU(i-1,j,k,uc))
                  gradU(2,1)=0.5_WP*dyi*(pU(i,j+1,k,uc)-pU(i,j-1,k,uc))
                  gradU(3,1)=0.5_WP*dzi*(pU(i,j,k+1,uc)-pU(i,j,k-1,uc))
                  gradU(1,2)=0.5_WP*dxi*(pV(i+1,j,k,vc)-pV(i-1,j,k,vc))
                  gradU(2,2)=0.5_WP*dyi*(pV(i,j+1,k,vc)-pV(i,j-1,k,vc))
                  gradU(3,2)=0.5_WP*dzi*(pV(i,j,k+1,vc)-pV(i,j,k-1,vc))
                  gradU(1,3)=0.5_WP*dxi*(pW(i+1,j,k,wc)-pW(i-1,j,k,wc))
                  gradU(2,3)=0.5_WP*dyi*(pW(i,j+1,k,wc)-pW(i,j-1,k,wc))
                  gradU(3,3)=0.5_WP*dzi*(pW(i,j,k+1,wc)-pW(i,j,k-1,wc))
               end if
               ! A=gradU_ij*gradU_ij invariant
               Aij=sum(gradU**2)
               ! beta_ij=dx_m^2*gradU_mi*gradU_mj
               do sj=1,3; do si=1,3
                  betaij(si,sj)=dx**2*gradU(1,si)*gradU(1,sj)+dy**2*gradU(2,si)*gradU(2,sj)+dz**2*gradU(3,si)*gradU(3,sj)
               end do; end do
               ! B invariant
               Bij=betaij(1,1)*betaij(2,2)-betaij(1,2)**2+betaij(1,1)*betaij(3,3)-betaij(1,3)**2+betaij(2,2)*betaij(3,3)-betaij(2,3)**2
               ! Assemble eddy viscosity
               if (Bij.gt.0.0_WP) pVisc(i,j,k,1)=Cmodel*sqrt(Bij/Aij)
               ! Clip to CFL limit
               pVisc(i,j,k,1)=min(pVisc(i,j,k,1),max_visc)
            end do; end do; end do
         end do
         call visc%amr%mfiter_destroy(mfi)
         ! Phase 2: Filter
         call visc%amr%mfab_filter(lvl=lvl,mfab=visc%mf(lvl),npass=nf)
      end do

   end subroutine get_vreman

   !> Compute artificial bulk kinematic viscosity from velocity
   !> Input:  U(Ucomp), V(Vcomp), W(Wcomp) velocity, C sound speed (cell-centered)
   !> Output: visc is cell-centered kinematic artificial viscosity (ncomp=1, ng>=1)
   subroutine get_viscartif(dt,visc,U,V,W,C,Ucomp,Vcomp,Wcomp,Ccomp,Cartif,max_cfl,nfilter)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab,amrex_multifab_destroy
      use amrdata_class,    only: amrdata
      use messager,         only: die
      implicit none
      real(WP), intent(in) :: dt
      type(amrdata), intent(inout) :: visc         !< Output kinematic artificial viscosity (ncomp=1, ng>=1)
      type(amrdata), intent(in) :: U,V,W           !< Velocities
      type(amrdata), intent(in) :: C               !< Sound speed
      integer, intent(in), optional :: Ucomp,Vcomp,Wcomp   !< Component indices (default=1)
      integer, intent(in), optional :: Ccomp               !< Component index for C (default=1)
      real(WP), intent(in), optional :: Cartif
      real(WP), intent(in), optional :: max_cfl
      integer,  intent(in), optional :: nfilter
      ! Locals
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab) :: scratch
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pC,pVisc,pScratch
      real(WP) :: dxi,dyi,dzi,dx,dy,dz,max_visc,myCartif,CFLmax
      real(WP) :: dudy,dudz,dvdx,dvdz,dwdx,dwdy,vort,grad_div
      integer :: lvl,i,j,k,nf,uc,vc,wc,cc

      ! Check visc
      if (visc%ncomp.ne.1) call die('[amrsgs::get_viscartif] visc must have exactly 1 component')
      if (visc%ng.lt.1) call die('[amrsgs::get_viscartif] visc must have at least 1 ghost cell')

      ! Component indices
      if (present(Ucomp)) then; uc=Ucomp; else; uc=1; end if
      if (present(Vcomp)) then; vc=Vcomp; else; vc=1; end if
      if (present(Wcomp)) then; wc=Wcomp; else; wc=1; end if
      if (present(Ccomp)) then; cc=Ccomp; else; cc=1; end if

      ! Set model constant
      if (present(Cartif)) then; myCartif=Cartif; else; myCartif=default_Cartif; end if
      ! Clipping based on maximum CFL
      if (present(max_cfl)) then; CFLmax=max_cfl; else; CFLmax=default_max_cfl; end if
      ! Number of filtering passes
      if (present(nfilter)) then; nf=nfilter; else; nf=default_nfilter; end if

      ! Loop over levels
      do lvl=0,visc%amr%clvl()
         ! Grid spacings
         dx=visc%amr%dx(lvl); dxi=1.0_WP/dx
         dy=visc%amr%dy(lvl); dyi=1.0_WP/dy
         dz=visc%amr%dz(lvl); dzi=1.0_WP/dz
         ! Max visc from CFL
         max_visc=CFLmax*visc%amr%min_meshsize(lvl)**2/(4.0_WP*dt)
         ! Build scratch for divergence (needs 1 ghost for grad_div stencil)
         call visc%amr%mfab_build(lvl=lvl,mfab=scratch,ncomp=1,nover=1); call scratch%setval(0.0_WP)
         ! Zero out visc
         call visc%mf(lvl)%setval(0.0_WP)
         ! Phase 1: Compute divergence into scratch
         call visc%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pU=>U%mf(lvl)%dataptr(mfi)
            pV=>V%mf(lvl)%dataptr(mfi)
            pW=>W%mf(lvl)%dataptr(mfi)
            pScratch=>scratch%dataptr(mfi)
            bx=mfi%growntilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pScratch(i,j,k,1)=0.5_WP*(dxi*(pU(i+1,j,k,uc)-pU(i-1,j,k,uc))+dyi*(pV(i,j+1,k,vc)-pV(i,j-1,k,vc))+dzi*(pW(i,j,k+1,wc)-pW(i,j,k-1,wc)))
            end do; end do; end do
         end do
         call visc%amr%mfiter_destroy(mfi)
         ! Phase 2: Compute kinematic visc
         call visc%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pU=>U%mf(lvl)%dataptr(mfi)
            pV=>V%mf(lvl)%dataptr(mfi)
            pW=>W%mf(lvl)%dataptr(mfi)
            pC=>C%mf(lvl)%dataptr(mfi)
            pScratch=>scratch%dataptr(mfi)
            pVisc=>visc%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Only work in compression regions
               if (pScratch(i,j,k,1).ge.0.0_WP) cycle
               ! Compute local vorticity
               dudy=0.5_WP*dyi*(pU(i,j+1,k,uc)-pU(i,j-1,k,uc))
               dudz=0.5_WP*dzi*(pU(i,j,k+1,uc)-pU(i,j,k-1,uc))
               dvdx=0.5_WP*dxi*(pV(i+1,j,k,vc)-pV(i-1,j,k,vc))
               dvdz=0.5_WP*dzi*(pV(i,j,k+1,vc)-pV(i,j,k-1,vc))
               dwdx=0.5_WP*dxi*(pW(i+1,j,k,wc)-pW(i-1,j,k,wc))
               dwdy=0.5_WP*dyi*(pW(i,j+1,k,wc)-pW(i,j-1,k,wc))
               vort=(dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2
               ! Compute |grad(div)|
               grad_div=max(abs(pScratch(i+1,j,k,1)-pScratch(i,j,k,1)),abs(pScratch(i,j,k,1)-pScratch(i-1,j,k,1)))*dx**2 &
               &       +max(abs(pScratch(i,j+1,k,1)-pScratch(i,j,k,1)),abs(pScratch(i,j,k,1)-pScratch(i,j-1,k,1)))*dy**2 &
               &       +max(abs(pScratch(i,j,k+1,1)-pScratch(i,j,k,1)),abs(pScratch(i,j,k,1)-pScratch(i,j,k-1,1)))*dz**2
               ! Floor vorticity with sound speed
               vort=max(vort,(0.05_WP*pC(i,j,k,cc)/visc%amr%min_meshsize(lvl))**2)
               ! Compute visc
               pVisc(i,j,k,1)=myCartif*grad_div*min(4.0_WP/3.0_WP*pScratch(i,j,k,1)**2/(pScratch(i,j,k,1)**2+vort+1.0e-15_WP),1.0_WP)
               ! Clip to max
               pVisc(i,j,k,1)=min(pVisc(i,j,k,1),max_visc)
            end do; end do; end do
         end do
         call visc%amr%mfiter_destroy(mfi)
         ! Destroy scratch
         call amrex_multifab_destroy(scratch)
         ! Phase 3: Filter
         call visc%amr%mfab_filter(lvl=lvl,mfab=visc%mf(lvl),npass=nf)
      end do

   end subroutine get_viscartif

end module amrsgs
