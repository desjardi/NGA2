!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use matm_class,        only: matm
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use string,            only: str_medium
   use hypre_str_class,   only: hypre_str
   use mathtools,         only: Pi
   use iterator_class,    only: iterator
   implicit none
   private

   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   type(matm),        public :: matmod
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs
   type(iterator) :: top_layer, btm_layer 

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,dfile,cflfile,consfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Choice of relaxation model
   integer :: relax_model

   !> Problem definition
   real(WP) :: Reg,Weg,r_visc,r_rho,r_vel,delta_l,Ls,Pref_l,LP,GP,Grho,Lrho,Lyl,Lyg
   integer :: nwaveX,nwaveZ,nwave,mode_low
   real(WP), dimension(:), allocatable :: wnumbX,wshiftX,wampX,wnumbZ,wshiftZ,wampZ
   type(event) :: ppevt
   integer :: junit
   character(len=str_medium) :: filename,timestamp


contains

   !> Function that defines a level set function for a initial wavy interface
   function levelset_wavy(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      integer :: nX,nZ
      G=-xyz(2)
   end function levelset_wavy

   !> Function that localizes the top (y+) of the domain
   function top_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain

   !> Function that localizes the bottom (y-) of the domain
   function btm_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function btm_of_domain

   !> Function that localizes top sponge
   function top_sponge(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%y(pg%jmax+1)-pg%ym(j).le.Ls) isIn=.true.
   end function top_sponge

   !> Function that localizes bottom sponge
   function btm_sponge(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%ym(j)-pg%y(pg%jmin).le.Ls) isIn=.true.
   end function btm_sponge

   subroutine apply_sponges(t)
      use mathtools,  only: Pi
      implicit none
      real(WP), intent(in) :: t
      integer :: i,j,k,n,m
      logical :: in_sponge_top, in_sponge_btm
      real(WP) :: swt_top, swt_btm, psponge, rhos, usponge, vsponge, wsponge

      ! Apply sponges using iterators
      do n=1,top_layer%n_
         swt_top = min(1.0_WP,max(0.0_WP,(fs%cfg%ym(top_layer%map(2,n)) - Lyg + Ls)/Ls))**2
         ! Get sponge solution if within a sponge
         psponge = GP
         rhos = Grho
         usponge = 1.0_WP; wsponge = 0.0_WP; vsponge = 0.0_WP
         ! Apply changes to variables
         fs%Grho (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%Grho (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) & 
            - swt_top*(fs%Grho(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))-rhos)
         fs%Ui   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%Ui   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) & 
            - swt_top*(fs%Ui(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))-usponge)
         fs%Vi   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%Vi   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) & 
            - swt_top*(fs%Vi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))-vsponge)
         fs%Wi   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%Wi   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) & 
            - swt_top*(fs%Wi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))-wsponge)
         fs%GP   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%GP   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) & 
            - swt_top*(fs%GP(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))-psponge)
         fs%GrhoE(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%GrhoE(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) & 
            - swt_top*(fs%GrhoE(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))-matmod%EOS_energy(psponge,rhos,usponge,vsponge,wsponge,'gas'))
         ! Update related quantities
         fs%RHO  (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%Grho(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))
         fs%rhoUi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%RHO(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) &
            *fs%Ui(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))
         fs%rhoVi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%RHO(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) &
            *fs%Vi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))
         fs%rhoWi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = fs%RHO(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) &
            *fs%Wi(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n))
         ! Remove liquid
         vf%VF   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = 0.0_WP
         fs%Lrho (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = 0.0_WP
         fs%LP   (top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = 0.0_WP
         fs%LrhoE(top_layer%map(1,n),top_layer%map(2,n),top_layer%map(3,n)) = 0.0_WP
      end do

      do n=1,btm_layer%n_
         swt_btm = min(1.0_WP,max(0.0_WP,(-fs%cfg%ym(btm_layer%map(2,n)) - Lyl + Ls)/Ls))**2
         ! Get sponge solution if within a sponge
         psponge = GP
         rhos = Lrho
         usponge = -r_vel; wsponge = 0.0_WP; vsponge = 0.0_WP
         ! Apply changes to variables
         fs%Lrho (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%Lrho (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            - swt_btm*(fs%Lrho(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))-rhos)
         fs%Ui   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%Ui   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            - swt_btm*(fs%Ui(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))-usponge)
         fs%Vi   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%Vi   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            - swt_btm*(fs%Vi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))-vsponge)
         fs%Wi   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%Wi   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            - swt_btm*(fs%Wi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))-wsponge)
         fs%LP   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%LP   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            - swt_btm*(fs%LP(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))-psponge)
         fs%LrhoE(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%LrhoE(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            - swt_btm*(fs%LrhoE(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))-matmod%EOS_energy(psponge,rhos,usponge,vsponge,wsponge,'liquid'))
         ! Update related quantities
         fs%RHO  (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%Lrho(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))
         fs%rhoUi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%RHO(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            *fs%Ui(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))
         fs%rhoVi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%RHO(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            *fs%Vi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))
         fs%rhoWi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = fs%RHO(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) &
            *fs%Wi(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n))
         ! Remove gas
         vf%VF   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = 1.0_WP
         fs%Grho (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = 0.0_WP
         fs%GP   (btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = 0.0_WP
         fs%GrhoE(btm_layer%map(1,n),btm_layer%map(2,n),btm_layer%map(3,n)) = 0.0_WP
      end do

      ! Reset interface-based quantities due to VF changes
      call vf%remove_flotsams()
      call vf%advect_interface(0.0_WP,fs%U,fs%V,fs%W)
      call vf%build_interface()
      call vf%remove_thinstruct()
      call vf%polygonalize_interface()
      call vf%subcell_vol()
      call vf%get_curvature()
      call vf%reset_moments()

   end subroutine

   subroutine growth_rate()
      use mathtools,  only: Pi
      use mpi_f08,    only: MPI_ALLREDUCE,MPI_SUM 
      use parallel,   only: MPI_REAL_WP 
      use param,      only: param_read
      use string,     only: str_medium
      implicit none
      integer :: i,j,k,ierr,iunit
      real(WP), dimension(:), allocatable :: localrho,  totalrho
      real(WP), dimension(:), allocatable :: localrhoU, totalrhoU
      real(WP), dimension(:), allocatable :: localrhoV, totalrhoV
      real(WP), dimension(:), allocatable :: localrhoW, totalrhoW
      real(WP), dimension(:), allocatable :: localvol,  totalvol
      real(WP), dimension(:), allocatable :: localmass, totalmass
      real(WP), dimension(:), allocatable :: localrhoUU, totalrhoUU
      real(WP), dimension(:), allocatable :: localrhoVV, totalrhoVV
      real(WP), dimension(:), allocatable :: localrhoUV, totalrhoUV
      real(WP), dimension(:), allocatable :: localrhoWW, totalrhoWW
      real(WP), dimension(:), allocatable :: localU, totalU
      real(WP), dimension(:), allocatable :: localV, totalV
      real(WP), dimension(:), allocatable :: localW, totalW
      real(WP), dimension(:), allocatable :: localGvol, totalGvol
      real(WP), dimension(:), allocatable :: localLvol, totalLvol
      real(WP), dimension(:), allocatable :: localGmass, totalGmass
      real(WP), dimension(:), allocatable :: localLmass, totalLmass
      real(WP), dimension(:), allocatable :: totalGrho,totalLrho
      real(WP), dimension(:), allocatable :: localGrhoU, totalGrhoU
      real(WP), dimension(:), allocatable :: localLrhoU, totalLrhoU
      real(WP), dimension(:), allocatable :: tempGrho,tempLrho
      real(WP), dimension(:), allocatable :: tempGvol,tempLvol
      real(WP), dimension(:), allocatable :: testGrho,testLrho
      real(WP), dimension(:), allocatable :: dudy_favre,dudy_mean
      real(WP) :: localdmthick,totaldmthick,grate,totaldwthick,totaldgthick,totaldlthick
      real(WP) :: localnu, totalnu, Re_mom, Re_vort, Re_lambda, TKE, Reg_mom
      real(WP) :: uprime1, uprime2, up2, dupdx2, tmp, lambda
      character(len=str_medium) :: filename,timestamp
      ! Allocate and initialize to zero
      allocate(localvol(fs%cfg%jmin:fs%cfg%jmax));   localvol=0.0_WP
      allocate(totalvol(fs%cfg%jmin:fs%cfg%jmax));   totalvol=0.0_WP
      allocate(localrho(fs%cfg%jmin:fs%cfg%jmax));   localrho=0.0_WP
      allocate(totalrho(fs%cfg%jmin:fs%cfg%jmax));   totalrho=0.0_WP
      allocate(localrhoU(fs%cfg%jmin:fs%cfg%jmax));  localrhoU=0.0_WP
      allocate(totalrhoU(fs%cfg%jmin:fs%cfg%jmax));  totalrhoU=0.0_WP
      allocate(localrhoV(fs%cfg%jmin:fs%cfg%jmax));  localrhoV=0.0_WP
      allocate(totalrhoV(fs%cfg%jmin:fs%cfg%jmax));  totalrhoV=0.0_WP
      allocate(localrhoW(fs%cfg%jmin:fs%cfg%jmax));  localrhoW=0.0_WP
      allocate(totalrhoW(fs%cfg%jmin:fs%cfg%jmax));  totalrhoW=0.0_WP
      allocate(localrhoUU(fs%cfg%jmin:fs%cfg%jmax)); localrhoUU=0.0_WP
      allocate(totalrhoUU(fs%cfg%jmin:fs%cfg%jmax)); totalrhoUU=0.0_WP
      allocate(localrhoVV(fs%cfg%jmin:fs%cfg%jmax)); localrhoVV=0.0_WP
      allocate(totalrhoVV(fs%cfg%jmin:fs%cfg%jmax)); totalrhoVV=0.0_WP
      allocate(localrhoUV(fs%cfg%jmin:fs%cfg%jmax)); localrhoUV=0.0_WP
      allocate(totalrhoUV(fs%cfg%jmin:fs%cfg%jmax)); totalrhoUV=0.0_WP
      allocate(localrhoWW(fs%cfg%jmin:fs%cfg%jmax)); localrhoWW=0.0_WP
      allocate(totalrhoWW(fs%cfg%jmin:fs%cfg%jmax)); totalrhoWW=0.0_WP
      allocate(localmass(fs%cfg%jmin:fs%cfg%jmax));  localmass=0.0_WP
      allocate(totalmass(fs%cfg%jmin:fs%cfg%jmax));  totalmass=0.0_WP
      allocate(dudy_favre(fs%cfg%jmin:fs%cfg%jmax)); dudy_favre=0.0_WP
      allocate(dudy_mean(fs%cfg%jmin:fs%cfg%jmax));  dudy_mean=0.0_WP
      allocate(localU(fs%cfg%jmin:fs%cfg%jmax));     localU=0.0_WP
      allocate(totalU(fs%cfg%jmin:fs%cfg%jmax));     totalU=0.0_WP
      allocate(localV(fs%cfg%jmin:fs%cfg%jmax));     localV=0.0_WP
      allocate(totalV(fs%cfg%jmin:fs%cfg%jmax));     totalV=0.0_WP
      allocate(localW(fs%cfg%jmin:fs%cfg%jmax));     localW=0.0_WP
      allocate(totalW(fs%cfg%jmin:fs%cfg%jmax));     totalW=0.0_WP
      allocate(localGvol(fs%cfg%jmin:fs%cfg%jmax));   localGvol=0.0_WP
      allocate(totalGvol(fs%cfg%jmin:fs%cfg%jmax));   totalGvol=0.0_WP
      allocate(localLvol(fs%cfg%jmin:fs%cfg%jmax));   localLvol=0.0_WP
      allocate(totalLvol(fs%cfg%jmin:fs%cfg%jmax));   totalLvol=0.0_WP
      allocate(localGmass(fs%cfg%jmin:fs%cfg%jmax));   localGmass=0.0_WP
      allocate(totalGmass(fs%cfg%jmin:fs%cfg%jmax));   totalGmass=0.0_WP
      allocate(localLmass(fs%cfg%jmin:fs%cfg%jmax));   localLmass=0.0_WP
      allocate(totalLmass(fs%cfg%jmin:fs%cfg%jmax));   totalLmass=0.0_WP
      allocate(totalGrho(fs%cfg%jmin:fs%cfg%jmax));   totalGrho=0.0_WP
      allocate(totalLrho(fs%cfg%jmin:fs%cfg%jmax));   totalLrho=0.0_WP
      allocate(localGrhoU(fs%cfg%jmin:fs%cfg%jmax));  localGrhoU=0.0_WP
      allocate(totalGrhoU(fs%cfg%jmin:fs%cfg%jmax));  totalGrhoU=0.0_WP
      allocate(localLrhoU(fs%cfg%jmin:fs%cfg%jmax));  localLrhoU=0.0_WP
      allocate(totalLrhoU(fs%cfg%jmin:fs%cfg%jmax));  totalLrhoU=0.0_WP
      allocate(tempGvol(fs%cfg%jmin:fs%cfg%jmax));  tempGvol=0.0_WP
      allocate(tempLvol(fs%cfg%jmin:fs%cfg%jmax));  tempLvol=0.0_WP
      allocate(tempGrho(fs%cfg%jmin:fs%cfg%jmax));  tempGrho=0.0_WP
      allocate(tempLrho(fs%cfg%jmin:fs%cfg%jmax));  tempLrho=0.0_WP
      allocate(testGrho(fs%cfg%jmin:fs%cfg%jmax));  testGrho=0.0_WP
      allocate(testLrho(fs%cfg%jmin:fs%cfg%jmax));  testLrho=0.0_WP
      localdmthick=0.0_WP;totaldmthick=0.0_WP;grate=0.0_WP;totaldgthick=0.0_WP;totaldlthick=0.0_WP
      Re_mom=0.0_WP;Re_vort=0.0_WP;Re_lambda=0.0_WP;TKE=0.0_WP;localnu=0.0_WP;totalnu=0.0_WP
      uprime1=0.0_WP;uprime2=0.0_WP;up2=0.0_WP;dupdx2=0.0_WP;tmp=0.0_WP;lambda=0.0_WP
      ! Calculate spatial averages
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               localmass(j) = localmass(j) + fs%RHO(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localvol(j)  = localvol(j) + fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localrhoU(j) = localrhoU(j) + fs%RHO(i,j,k)*fs%Ui(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localrhoV(j) = localrhoV(j) + fs%RHO(i,j,k)*fs%Vi(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localrhoW(j) = localrhoW(j) + fs%RHO(i,j,k)*fs%Wi(i,j,k)*fs%cfg%vol(i,j,k)
               localU(j)    = localU(j) + fs%Ui(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localV(j)    = localV(j) + fs%Vi(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localW(j)    = localW(j) + fs%Wi(i,j,k)*fs%cfg%vol(i,j,k)
               ! Volume fraction specific calculations
               localGvol(j) = localGvol(j) + (1.0_WP - vf%VF(i,j,k))*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localLvol(j) = localLvol(j) + vf%VF(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localGmass(j) = localGmass(j) + (1.0_WP - vf%VF(i,j,k))*fs%Grho(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localLmass(j) = localLmass(j) + vf%VF(i,j,k)*fs%Lrho(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localGrhoU(j) = localGrhoU(j) + (1.0_WP - vf%VF(i,j,k))*fs%Grho(i,j,k)*fs%Ui(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               localLrhoU(j) = localLrhoU(j) + vf%VF(i,j,k)*fs%Lrho(i,j,k)*fs%Ui(i,j,k)*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k) !*fs%cfg%vol(i,j,k)
               ! Kinematic viscosity calculations
               localnu = localnu + fs%visc(i,j,k)*fs%cfg%vol(i,j,k)/fs%RHO(i,j,k)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE(localmass,totalmass,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localvol,totalvol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localrhoU,totalrhoU,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localrhoV,totalrhoV,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localrhoW,totalrhoW,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localU,totalU,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localV,totalV,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localW,totalW,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localGmass,totalGmass,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localGvol,totalGvol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localLmass,totalLmass,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localLvol,totalLvol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localGrhoU,totalGrhoU,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localLrhoU,totalLrhoU,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localnu,totalnu,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      ! tempGvol = totalGvol
      ! tempLvol = totalLvol
      do j=fs%cfg%jmin_,fs%cfg%jmax_
         if (totalGvol(j) .ne. 0.0_WP) then
            testGrho(j) = totalGmass(j)/totalGvol(j)
            totalGrhoU(j) = totalGrhoU(j)/totalGvol(j) 
         end if
         if (totalLvol(j) .ne. 0.0_WP) then
            testLrho(j) = totalLmass(j)/totalLvol(j)
            totalLrhoU(j) = totalLrhoU(j)/totalLvol(j)
         end if
      end do
      ! Volume-average for each plane
      totalrho  = totalmass/totalvol
      totalrhoU = totalrhoU/totalvol
      totalrhoV = totalrhoV/totalvol
      totalrhoW = totalrhoW/totalvol
      totalU    = totalU/totalvol
      totalV    = totalV/totalvol
      totalGrho = totalGmass/totalvol
      totalLrho = totalLmass/totalvol
      ! totalGrhoU = totalGrhoU/totalvol
      ! totalLrhoU = totalLrhoU/totalvol
      ! tempGrho = totalGrho
      ! tempLrho = totalLrho
      ! Average kinematic viscosity
      totalnu   = totalnu/fs%cfg%vol_total
      ! Calculate Momentum thickness
      do j=fs%cfg%jmin_,fs%cfg%jmax_
         if (totalGrho(j) .ne. 0.0_WP) then
            ! totaldgthick = totaldgthick + (1.0_WP/(Grho*(1.0_WP+r_vel)**2.0_WP))*(totalGrho(j)*(1.0_WP - (totalGrhoU(j)/testGrho(j)))*((totalGrhoU(j)/testGrho(j)) + r_vel))*vf%cfg%dy(j)
            ! totaldgthick = totaldgthick + (1.0_WP/(totalrho(j)*(1.0_WP+r_vel)**2.0_WP))*(totalGrho(j)*(1.0_WP - (totalGrhoU(j)/testGrho(j)))*((totalGrhoU(j)/testGrho(j)) + r_vel))*vf%cfg%dy(j)
            ! This is method 1
            totaldgthick = totaldgthick + (1.0_WP/(Grho*(1.0_WP+r_vel)**2.0_WP))*(totalGrho(j)*(1.0_WP - (totalGrhoU(j)/testGrho(j)))*((totalGrhoU(j)/testGrho(j)) + r_vel))*vf%cfg%dy(j)
            ! This is method 2
            ! totaldgthick = totaldgthick + (1.0_WP/(((Grho+Lrho)/2.0_WP)*(1.0_WP+r_vel)**2.0_WP))*(totalGrho(j)*(1.0_WP - (totalrhoU(j)/totalrho(j)))*((totalrhoU(j)/totalrho(j)) + r_vel))*vf%cfg%dy(j)
         end if
         if (totalLrho(j) .ne. 0.0_WP) then
            ! totaldlthick = totaldlthick + (1.0_WP/(Lrho*(1.0_WP+r_vel)**2.0_WP))*(totalLrho(j)*(1.0_WP - (totalLrhoU(j)/testLrho(j)))*((totalLrhoU(j)/testLrho(j)) + r_vel))*vf%cfg%dy(j)
            ! totaldlthick = totaldlthick + (1.0_WP/(Lrho*(1.0_WP+r_vel)**2.0_WP))*(totalLrho(j)*(1.0_WP - (totalLrhoU(j)/testLrho(j)))*((totalLrhoU(j)/testLrho(j)) + r_vel))*vf%cfg%dy(j)
            ! This is method 1
            totaldlthick = totaldlthick + (1.0_WP/(Lrho*(1.0_WP+r_vel)**2.0_WP))*(totalLrho(j)*(1.0_WP - (totalLrhoU(j)/testLrho(j)))*((totalLrhoU(j)/testLrho(j)) + r_vel))*vf%cfg%dy(j)
            ! This is method 2
            ! totaldlthick = totaldlthick + (1.0_WP/(((Grho+Lrho)/2.0_WP)*(1.0_WP+r_vel)**2.0_WP))*(totalLrho(j)*(1.0_WP - (totalrhoU(j)/totalrho(j)))*((totalrhoU(j)/totalrho(j)) + r_vel))*vf%cfg%dy(j)
         end if
         ! totaldmthick = totaldmthick + (1.0_WP/(totalrho(j)*(1.0_WP+r_vel)**2.0_WP))*(totalrho(j)*(1.0_WP - (totalrhoU(j)/totalrho(j)))*((totalrhoU(j)/totalrho(j)) + r_vel))*vf%cfg%dy(j)
         totaldmthick = totaldmthick + (1.0_WP/(((Grho+Lrho)/2.0_WP)*(1.0_WP+r_vel)**2.0_WP))*(totalrho(j)*(1.0_WP - (totalrhoU(j)/totalrho(j)))*((totalrhoU(j)/totalrho(j)) + r_vel))*vf%cfg%dy(j)
         ! totaldgthick = totaldgthick + (1.0_WP/(Grho*(1.0_WP+r_vel)**2.0_WP))*(totalGrho(j)*(1.0_WP - (totalGrhoU(j)/tempGrho(j)))*((totalGrhoU(j)/tempGrho(j)) + r_vel))*vf%cfg%dy(j)
         ! totaldlthick = totaldlthick + (1.0_WP/(Lrho*(1.0_WP+r_vel)**2.0_WP))*(totalLrho(j)*(1.0_WP - (totalLrhoU(j)/tempLrho(j)))*((totalLrhoU(j)/tempLrho(j)) + r_vel))*vf%cfg%dy(j)
      end do
      ! Calculate Reynolds Stresses
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               ! localrhoUU(j) = localrhoUU(j) + fs%RHO(i,j,k)*(fs%Ui(i,j,k) - totalrhoU(j)/totalrho(j))**2 &
               !    *fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               ! localrhoVV(j) = localrhoVV(j) + fs%RHO(i,j,k)*(fs%Vi(i,j,k) - totalrhoV(j)/totalrho(j))**2 &
               !    *fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               ! localrhoUV(j) = localrhoUV(j) + fs%RHO(i,j,k)*(fs%Ui(i,j,k) - totalrhoU(j)/totalrho(j)) &
               !    *(fs%Vi(i,j,k) - totalrhoV(j)/totalrho(j))*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               ! Alternative attempt for Reynolds Stresses
               localrhoUU(j) = localrhoUU(j) + (fs%Ui(i,j,k) - totalrhoU(j)/totalrho(j))**2 &
                  *fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               localrhoVV(j) = localrhoVV(j) + (fs%Vi(i,j,k) - totalrhoV(j)/totalrho(j))**2 &
                  *fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               localrhoUV(j) = localrhoUV(j) + (fs%Ui(i,j,k) - totalrhoU(j)/totalrho(j)) &
                  *(fs%Vi(i,j,k) - totalrhoV(j)/totalrho(j))*fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               localrhoWW(j) = localrhoWW(j) + (fs%Wi(i,j,k) - totalrhoW(j)/totalrho(j))**2 &
                  *fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
            end do
         end do
      end do
      ! Calculate vorticity thickness and favre velocity gradient component
      do j=fs%cfg%jmin_+1,fs%cfg%jmax_-1
         dudy_mean(j)  = ABS((totalU(j+1) - totalU(j-1))/(2.0_WP*fs%cfg%dy(j)))
         dudy_favre(j) = ((totalrhoU(j+1)/totalrho(j+1)) - (totalrhoU(j-1)/totalrho(j-1)))/(2.0_WP*fs%cfg%dy(j))
      end do
      totaldwthick = (1.0_WP + r_vel) / MAXVAL(dudy_mean)
      ! All-reduce the data
      call MPI_ALLREDUCE(localrhoUU,totalrhoUU,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localrhoVV,totalrhoVV,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localrhoUV,totalrhoUV,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(localrhoWW,totalrhoWW,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      totalrhoUU = (totalrho*totalrhoUU)/totalvol
      totalrhoVV = (totalrho*totalrhoVV)/totalvol
      totalrhoUV = (totalrho*totalrhoUV)/totalvol
      totalrhoWW = (totalrho*totalrhoWW)/totalvol
      ! Calculate Growth Rate
      do j=fs%cfg%jmin_,fs%cfg%jmax_
         ! grate = grate - (totalrhoUV(j)*dudy_favre(j))*vf%cfg%dy(j)
         grate = grate - (2.0_WP/(((Grho + Lrho)/2)*(1.0_WP+r_vel)**2))*(totalrhoUV(j)*dudy_favre(j+3))*vf%cfg%dy(j)
      end do
      ! Calculate Reynolds Numbers and Taylor microscale (lambda)
      Re_mom = totaldmthick*(1.0_WP + r_vel)/(totalnu+epsilon(1.0_WP))
      Re_vort = totaldwthick*(1.0_WP + r_vel)/(totalnu+epsilon(1.0_WP))
      ! Reg_mom = totaldgthick*1.0_WP/(visc_g+epsilon(1.0_WP))
      j = (fs%cfg%nyo-1)/2+1
      ! Check whether the Tayler microscale changes with compressibility 
      if (j.ge.fs%cfg%jmin_ .and. j.le.fs%cfg%jmax_) then
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               uprime1 = fs%Ui(i,j,k)-totalU(j)
               uprime2 = fs%Ui(i+1,j,k)-totalU(j)
               up2 = up2 + ((uprime1+uprime2)/2.0_WP)**2
               dupdx2 = dupdx2 + ((uprime2-uprime1)/fs%cfg%dx(i))**2
            end do
         end do
      else
         up2  = 0.0_WP
         dupdx2 = 0.0_WP
      end if
      call MPI_ALLREDUCE(up2,tmp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      up2 = tmp
      call MPI_ALLREDUCE(dupdx2,tmp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      dupdx2 = tmp
      lambda = sqrt( up2/(dupdx2+epsilon(1.0_WP)) )
      Re_lambda = lambda*up2**(0.5)/(totalnu+epsilon(1.0_WP))
      ! Calculate turbulent kinetic energy
      tmp=0.0_WP
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               tmp = tmp + fs%cfg%vol(i,j,k)*(fs%Ui(i,j,k)-(totalrhoU(j)/totalrho(j)))**2 & 
                  + (fs%Vi(i,j,k)-(totalrhoV(j)/totalrho(j)))**2 + (fs%Wi(i,j,k)-(totalrhoW(j)/totalrho(j)))**2
            end do
         end do
      end do
      call MPI_ALLREDUCE(tmp,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      TKE = TKE/(2.0_WP*fs%cfg%vol_total)
      ! If root, print it out
      if (fs%cfg%amRoot) then
         ! Momentum thickness and growth rate
         write(junit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') time%t,totaldmthick,totaldwthick,grate,totaldgthick,totaldlthick,Re_mom,TKE
         ! Reynolds Stresses
         filename='PostProc_'
         write(timestamp,'(f6.1)') time%t
         open(newunit=iunit,file='PostProc/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12)') &
            'Height','totalrho','totalrhoU','R11','R22','R12','R33','Grho','Lrho','GrhoU','LrhoU','U_favreG','U_favreL'
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') & 
               fs%cfg%ym(j),totalrho(j),totalrhoU(j),totalrhoUU(j),totalrhoVV(j),totalrhoUV(j),totalrhoWW(j),totalGrho(j),totalLrho(j),totalGrhoU(j),totalLrhoU(j),totalGrhoU(j)/totalGrho(j),totalLrhoU(j)/totalLrho(j)
         end do
         close(iunit)
      end if
      ! Deallocate work arrays
      deallocate(localrho,totalrho)
      deallocate(localrhoU,totalrhoU)
      deallocate(localvol,totalvol)
      deallocate(localrhoV,totalrhoV)
      deallocate(localrhoW,totalrhoW)
      deallocate(localrhoUU,totalrhoUU)
      deallocate(localrhoVV,totalrhoVV)
      deallocate(localrhoWW,totalrhoWW)
      deallocate(localrhoUV,totalrhoUV)
      deallocate(localmass,totalmass)
      deallocate(dudy_favre,dudy_mean)
      deallocate(localU,totalU)
      deallocate(localV,totalV)
      deallocate(localW,totalW)
      deallocate(localGvol,totalGvol)
      deallocate(localLvol,totalLvol)
      deallocate(localGmass,totalGmass)
      deallocate(localLmass,totalLmass)
      deallocate(totalGrho,totalLrho)
      deallocate(localGrhoU,totalGrhoU)
      deallocate(localLrhoU,totalLrhoU)
      deallocate(tempGrho,tempLrho)
      deallocate(tempGvol,tempLvol)
   end subroutine growth_rate


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo
         use mathtools, only: Pi,twoPi
         use random,    only: random_uniform
         use parallel,  only: MPI_REAL_WP
         use mpi_f08
         integer :: i,j,k,n,si,sj,sk,ierr
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         ! vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF') 
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_wavy,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! ! Set interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a two-phase flow solver
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use mast_class,      only: mech_egy_mech_hhz,clipped_neumann,bc_scope,bcond
         use mathtools,       only: Pi,twoPi
         use parallel,        only: amRoot
         use string,          only: str_medium
         use random,          only: random_uniform
         use parallel,        only: MPI_REAL_WP
         use mpi_f08,         only: MPI_BCAST
         use c_interface
         integer :: i,j,k,n,ierr
         real(WP) :: gamm_l,gamm_g,visc_l,visc_g,cv_g,cv_l,T
         real(WP) :: Ma_g,Ma_l,c
         real(WP) :: diffx_left,diffx_right,diffy_left,diffy_right,diffz_left,diffz_right
         real(WP), dimension(:,:,:), allocatable :: Upert,Vpert,Wpert
         real(WP), dimension(:), allocatable :: test
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         ! Get nondimensional parameters from input
         call param_read('Gas gamma',gamm_g)
         call param_read('Liquid gamma',gamm_l)
         call param_read('Gas Reynolds number',Reg); visc_g=1.0_WP/(Reg+epsilon(Reg))
         call param_read('Viscosity ratio',r_visc); visc_l=r_visc*visc_g
         call param_read('Density ratio',r_rho); Grho=1.0_WP; Lrho=r_rho*Grho
         call param_read('Velocity ratio',r_vel); delta_l=r_visc*r_vel
         call param_read('Gas Mach number',Ma_g); GP = 1.0_WP/(gamm_g*Ma_g**2.0_WP)
         call param_read('Liquid Mach number',Ma_l); Pref_l = ((r_rho*r_vel**2)/(gamm_l*Ma_l**2.0_WP)) - GP
         Pref_l = 0.0_WP
         ! Register equations of state
         call matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
         call matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g,viscconst_liquid=visc_l)
         ! Set initial fields
         fs%Lrho = Lrho; fs%Grho = Grho; fs%LP = GP; fs%GP = GP
         ! Read in necessary geometry
         call param_read('Sponge length',Ls)
         call param_read('Liquid Ly',Lyl)
         call param_read('Gas Ly',Lyg)
         call param_read('Gas Weber number',Weg); fs%sigma=1.0_WP/(Weg+epsilon(Weg))
         fs%sigma = 0.0_WP
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Set initial velocity field
         c = (1.0_WP + r_vel)**2.0_WP / ((r_vel*LOG(2.0_WP) - (LOG(2.0_WP) - 1.0_WP)))
         fs%Ui=0.0_WP; fs%Vi=0.0_WP; fs%Wi=0.0_WP
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (fs%cfg%ym(j).le.0.0_WP) then
                     ! Use the liquid profile (error function)
                     ! fs%Ui(i,j,k)=r_vel*erf(fs%cfg%ym(j)/delta_l)
                     ! Use the liquid profile (hyperbolic tangent)
                     fs%Ui(i,j,k)=r_vel*tanh(fs%cfg%ym(j)/(2.0_WP*delta_l))
                  else
                     ! Use the gas profile (error function)
                     ! fs%Ui(i,j,k)=erf(fs%cfg%ym(j))
                     ! Use the gas profile (hyperbolic tangent)
                     fs%Ui(i,j,k)=tanh(fs%cfg%ym(j)/2.0_WP)
                  end if
               end do
            end do
         end do
         ! Set velocity potential perturbation
         mode_low = 3; nwave=20
         allocate(wnumbX(nwave),wshiftX(nwave))
         if (cfg%amRoot) then
            do n=1,nwave-mode_low
               wnumbX(n)=(mode_low+n-1)*twoPi/cfg%xL
               wshiftX(n)=random_uniform(lo=0.0_WP,hi=twoPi)
            end do
         end if
         call MPI_BCAST(wshiftX,nwave-mode_low,MPI_REAL_WP,0,cfg%comm,ierr)
         call MPI_BCAST(wnumbX,nwave-mode_low,MPI_REAL_WP,0,cfg%comm,ierr)
         allocate(wnumbZ(nwave),wshiftZ(nwave))
         if (cfg%amRoot) then
            do n=1,nwave-mode_low
               wnumbZ(n)=(mode_low+n-1)*twoPi/cfg%zL
               wshiftZ(n)=random_uniform(lo=0.0_WP,hi=twoPi)
            end do
         end if
         call MPI_BCAST(wshiftZ,nwave-mode_low,MPI_REAL_WP,0,cfg%comm,ierr)
         call MPI_BCAST(wnumbZ,nwave-mode_low,MPI_REAL_WP,0,cfg%comm,ierr)
         ! Add components of velocity potential to mean flow
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  do n=1,nwave-mode_low
                  fs%Ui(i,j,k)=fs%Ui(i,j,k) + 0.15_WP/(4.0_WP*Pi**2.0_WP*wnumbX(n)*wnumbZ(n))*wnumbX(n)*sin(wnumbX(n)*fs%cfg%xm(i)+wshiftX(n)) &
                  *cos(wnumbZ(n)*fs%cfg%zm(k) + wshiftZ(n))*EXP(-5.0_WP*fs%cfg%ym(j)**2.0_WP)*(sin(fs%cfg%ym(j))+10.0_WP*fs%cfg%ym(j)*cos(fs%cfg%ym(j))) 

                  fs%Vi(i,j,k)=fs%Vi(i,j,k) + 0.15_WP/(4.0_WP*Pi**2.0_WP*wnumbX(n)*wnumbZ(n)) &
                  *cos(wnumbX(n)*fs%cfg%xm(i) + wshiftX(n))*cos(wnumbZ(n)*fs%cfg%zm(k) + wshiftZ(n)) & 
                  *(EXP(-5.0_WP*fs%cfg%ym(j)**2.0_WP)*(20.0_WP*fs%cfg%ym(j)*sin(fs%cfg%ym(j))+(100.0_WP*fs%cfg%ym(j)**2.0_WP-11.0_WP)*cos(fs%cfg%ym(j))))

                  fs%Wi(i,j,k)=fs%Wi(i,j,k) + 0.15_WP/(4.0_WP*Pi**2.0_WP*wnumbX(n)*wnumbZ(n))*wnumbZ(n)*sin(wnumbZ(n)*fs%cfg%zm(k)+wshiftZ(n)) &
                  *cos(wnumbX(n)*fs%cfg%xm(i) + wshiftX(n))*EXP(-5.0_WP*fs%cfg%ym(j)**2.0_WP)*(sin(fs%cfg%ym(j))+10.0_WP*fs%cfg%ym(j)*cos(fs%cfg%ym(j))) 
                  end do
               end do
            end do
         end do
         ! Initialize pressure and energy based on modified velocity field
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  if (fs%cfg%ym(j).le.0.0_WP) then
                     fs%LrhoE(i,j,k) = matmod%EOS_energy(GP,Lrho,fs%Ui(i,j,k),fs%Vi(i,j,k),fs%Wi(i,j,k),'liquid')
                     fs%LP(i,j,k) = matmod%EOS_liquid(i,j,k,'p')
                  else
                     fs%GrhoE(i,j,k) = matmod%EOS_energy(GP,Grho,fs%Ui(i,j,k),fs%Vi(i,j,k),fs%Wi(i,j,k),'gas')
                     fs%GP(i,j,k) = matmod%EOS_gas(i,j,k,'p')
                  end if
               end do
            end do
         end do

         ! Define BCs at top and bottom - though sponges will determine bdy behavior
         call fs%add_bcond(name='top_y'   ,type=clipped_neumann,locator=top_of_domain  ,celldir='yp')
         call fs%add_bcond(name='btm_y'   ,type=clipped_neumann,locator=btm_of_domain  ,celldir='ym')
         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         ! Use mechanical energy relaxation
         relax_model = mech_egy_mech_hhz
         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         ! Perform initial pressure relax
         call fs%pressure_relax(vf,matmod,relax_model)
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA
         ! Initialize first guess for pressure (0 works best)
         fs%psolv%sol=0.0_WP
         ! Store momentum thickness results
         if (fs%cfg%amRoot) then
            open(newunit=junit,file='thickness.csv',status='unknown')
            write(junit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'Time','d_m','d_w','Growth Rate','d_g','d_l','Re_mom','TKE'
            call execute_command_line('mkdir -p PostProc')
            ! if the execute_command_line intrinsic still gives me trouble, use the following
            ! call F_mkdir('PostProc')
         end if
      end block create_and_initialize_flow_solver


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='CompressibleMultiphaseML')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('L_Pres',fs%LP)
         call ens_out%add_scalar('G_Pres',fs%GP)
         ! call ens_out%add_scalar('L_Energy',fs%LrhoE)
         ! call ens_out%add_scalar('G_Energy',fs%GrhoE)
         ! call ens_out%add_scalar('PA',fs%PA)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_vector('Momentum',fs%rhoUi,fs%rhoVi,fs%rhoWi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('Mach',fs%Mach)
         call ens_out%add_scalar('BulkMod',fs%RHOSS2)
         ! call ens_out%add_scalar('Temp',fs%Tmptr)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call fs%get_viz()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor

      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Perform the output
         if (ppevt%occurs()) then 
            call growth_rate()
         end if
      end block create_postproc 

      ! Create an iterator for imposing the bottom and top sponges
      create_iterator_top: block
         top_layer = iterator(pg=fs%cfg,name='Top Sponge',locator=top_sponge)
      end block create_iterator_top

      create_iterator_btm: block
         btm_layer = iterator(pg=fs%cfg,name='Bottom Sponge',locator=btm_sponge)
      end block create_iterator_btm

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      integer :: i,j,k

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()

         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Diffusion and built-in source term (gravity) step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

            ! Perform sponge forcing
            call apply_sponges(time%t)

            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod) 

            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)

            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)

            ! Record convergence monitor
            call cvgfile%write()

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         fs%PA = matmod%EOS_all(vf);
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call fs%get_viz()
         call mfile%write()
         call cflfile%write()
         call dfile%write()
         if (ppevt%occurs()) then
            call growth_rate()
         end if
      end do

      close(junit)

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      ! implicit none

   end subroutine simulation_final


end module simulation
