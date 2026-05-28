!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Rcyl
   use spcomp_class,      only: spcomp
   use gp_class,          only: gpibm
   use timetracker_class, only: timetracker
   use timer_class,       only: timer
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(spcomp),      public :: fs
   type(gpibm),       public :: gp
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile,ibmfile

   !> Timing info
   type(monitor) :: timefile !< Timing monitoring
   type(timer)   :: tstep    !< Timer for step
   type(timer)   :: tibm     !< Timer for IBM
   type(timer)   :: tcom     !< Timer for solver
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc,visc_t,div

   !> Constant kinematic viscosity
   real(WP) :: visc0,T0

   !> Equations of state
   real(WP) :: Pinf,Gamma,Cv,Prandtl

   !> Flow parameters
   real(WP) :: Ms,Xs
   real(WP) :: rho1,p1,u1,M1
   real(WP) :: rho2,p2,u2,M2
   real(WP) :: Re

   !> IBM parameters
   real(WP), dimension(3) :: ibm_force
   
 contains


   !> Function that returns a smooth Heaviside of thickness delta
   real(WP) function Hshock(x,delta)
     real(WP), intent(in) :: x,delta
     ! Goes from 0 to 1 as x goes from begative to positive
     Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock


   !> P=EOS(RHO,I)
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
   !> RHO=f(T,P)
   pure real(WP) function get_RHO(T,P)
     implicit none
     real(WP), intent(in) :: T,P
     get_RHO=(P+Pinf)/(Cv*T*(Gamma-1.0_WP))
   end function get_RHO
   !> I=EOS(RHO,P)
   pure real(WP) function get_I(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_I=(P+Gamma*Pinf)/(RHO*(Gamma-1.0_WP))
   end function get_I
   !> C=f(RHO,P)
   pure real(WP) function get_C(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_C=sqrt(Gamma*(P+Pinf)/RHO)
   end function get_C
   !> S=f(RHO,P)
   pure real(WP) function get_S(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_S=Cv*log((P+Pinf)/RHO**Gamma)
   end function get_S


   !> Calculate viscosities
   subroutine prepare_viscosities()
     implicit none
     integer :: i,j,k
     real(WP) :: S
     ! Get viscosity from Sutherland's law
     S=110.4_WP/273.15_WP*T0
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           do i=fs%cfg%imino_,fs%cfg%imaxo_
              visc(i,j,k)=visc0*(T0+S)/(fs%T(i,j,k)+S)*(fs%T(i,j,k)/T0)**1.5_WP
           end do
        end do
     end do
     ! Get LAD
     call fs%get_viscartif(dt=time%dt,beta=beta); fs%BETA=fs%Q(:,:,:,1)*beta
     ! Get eddy viscosity
     call fs%get_vreman   (dt=time%dt,visc=visc_t); fs%VISC=fs%Q(:,:,:,1)*visc_t+visc
     ! Recompute thermal conductivity
     fs%diff=Gamma*Cv*fs%visc/Prandtl
   end subroutine prepare_viscosities


   !> Calculate velocity divergence
   subroutine get_div()
     implicit none
     integer :: i,j,k
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_-1; do j=fs%cfg%jmino_,fs%cfg%jmaxo_-1; do i=fs%cfg%imino_,fs%cfg%imaxo_-1
        div(i,j,k)=fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k))
     end do; end do; end do
     call fs%cfg%sync(div)
     if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) div(fs%cfg%imaxo,:,:)=div(fs%cfg%imaxo-1,:,:)
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) div(:,fs%cfg%jmaxo,:)=div(:,fs%cfg%jmaxo-1,:)
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.fs%cfg%npz) div(:,:,fs%cfg%kmaxo)=div(:,:,fs%cfg%kmaxo-1)
   end subroutine get_div


   !> Compute force on cylinder
   subroutine get_force()
     use mathtools, only: Pi
     use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
     use parallel, only: MPI_REAL_WP
     implicit none
     integer :: i,j,k,n,ierr
     real(WP) :: div,Fl,Fr
     real(WP), dimension(:,:,:,:), allocatable :: FQx,FQy,FQz

     allocate(FQx(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_,1:3)); FQx=0.0_WP
     allocate(FQy(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_,1:3)); FQy=0.0_WP
     allocate(FQz(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_,1:3)); FQz=0.0_WP

     ! Compute cell-centered momentum fluxes
     do k=fs%cfg%kmin_-1,fs%cfg%kmax_
        do j=fs%cfg%jmin_-1,fs%cfg%jmax_
           do i=fs%cfg%imin_-1,fs%cfg%imax_
              div=fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k))
              FQx(i,j,k,1)=2.0_WP*fs%VISC(i,j,k)*fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+(fs%BETA(i,j,k)-2.0_WP*fs%VISC(i,j,k)/3.0_WP)*div-fs%P(i,j,k)
              FQy(i,j,k,2)=2.0_WP*fs%VISC(i,j,k)*fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+(fs%BETA(i,j,k)-2.0_WP*fs%VISC(i,j,k)/3.0_WP)*div-fs%P(i,j,k)
              FQz(i,j,k,3)=2.0_WP*fs%VISC(i,j,k)*fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k))+(fs%BETA(i,j,k)-2.0_WP*fs%VISC(i,j,k)/3.0_WP)*div-fs%P(i,j,k)
           end do
        end do
     end do

     ! Compute edge-centered momentum viscous fluxes and corresponding viscous heating
     do k=fs%cfg%kmin_,fs%cfg%kmax_+1
        do j=fs%cfg%jmin_,fs%cfg%jmax_+1
           do i=fs%cfg%imin_,fs%cfg%imax_+1
              FQy(i,j,k,1)=0.25_WP*sum(fs%VISC(i-1:i,j-1:j,k))*(fs%dyi*(fs%U(i,j,k)-fs%U(i,j-1,k))+fs%dxi*(fs%V(i,j,k)-fs%V(i-1,j,k))); FQx(i,j,k,2)=FQy(i,j,k,1)
              FQz(i,j,k,2)=0.25_WP*sum(fs%VISC(i,j-1:j,k-1:k))*(fs%dzi*(fs%V(i,j,k)-fs%V(i,j,k-1))+fs%dyi*(fs%W(i,j,k)-fs%W(i,j-1,k))); FQy(i,j,k,3)=FQz(i,j,k,2)
              FQx(i,j,k,3)=0.25_WP*sum(fs%VISC(i-1:i,j,k-1:k))*(fs%dxi*(fs%W(i,j,k)-fs%W(i-1,j,k))+fs%dzi*(fs%U(i,j,k)-fs%U(i,j,k-1))); FQz(i,j,k,1)=FQx(i,j,k,3)
           end do
        end do
     end do

     do i=1,3
        call fs%cfg%sync(FQx(:,:,:,i))
        call fs%cfg%sync(FQy(:,:,:,i))
        call fs%cfg%sync(FQz(:,:,:,i))
     end do

     ! Sum up force
     ibm_force=0.0_WP
     do k=cfg%kmin_,cfg%kmax_
        do j=cfg%jmin_,cfg%jmax_
           do i=cfg%imin_,cfg%imax_
              if (cfg%Gib(i,j,k).lt.0.0_WP) then
                 ! Force in x
                 Fl=fs%dxi*(FQx(i  ,j,k,1)-FQx(i-1,j,k,1))+fs%dyi*(FQy(i  ,j+1,k,1)-FQy(i  ,j,k,1))+fs%dzi*(FQz(i  ,j,k+1,1)-FQz(i  ,j,k,1))*fs%vol
                 Fr=fs%dxi*(FQx(i+1,j,k,1)-FQx(i  ,j,k,1))+fs%dyi*(FQy(i+1,j+1,k,1)-FQy(i+1,j,k,1))+fs%dzi*(FQz(i+1,j,k+1,1)-FQz(i+1,j,k,1))*fs%vol
                 ibm_force(1)=ibm_force(1)+0.5_WP*(Fl+Fr)
                 ! Force in y
                 Fl=fs%dxi*(FQx(i+1,j  ,k,2)-FQx(i,j  ,k,2))+fs%dyi*(FQy(i,j  ,k,2)-FQy(i,j-1,k,2))+fs%dzi*(FQz(i,j  ,k+1,2)-FQz(i,j  ,k,2))*fs%vol
                 Fr=fs%dxi*(FQx(i+1,j+1,k,2)-FQx(i,j+1,k,2))+fs%dyi*(FQy(i,j+1,k,2)-FQy(i,j  ,k,2))+fs%dzi*(FQz(i,j+1,k+1,2)-FQz(i,j+1,k,2))*fs%vol
                 ibm_force(2)=ibm_force(2)+0.5_WP*(Fl+Fr)
                 ! Force in z
                 Fl=fs%dxi*(FQx(i+1,j,k  ,3)-FQx(i,j,k  ,3))+fs%dyi*(FQy(i,j+1,k  ,3)-FQy(i,j,k  ,3))+fs%dzi*(FQz(i,j,k  ,3)-FQz(i,j,k-1,3))*fs%vol
                 FR=fs%dxi*(FQx(i+1,j,k+1,3)-FQx(i,j,k+1,3))+fs%dyi*(FQy(i,j+1,k+1,3)-FQy(i,j,k+1,3))+fs%dzi*(FQz(i,j,k+1,3)-FQz(i,j,k  ,3))*fs%vol
                 ibm_force(3)=ibm_force(3)+0.5_WP*(Fl+Fr)
              end if
           end do
        end do
     end do
     call MPI_ALLREDUCE(MPI_IN_PLACE,ibm_force,3,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)

     ! Deallocate flux arrays
     deallocate(FQx,FQy,FQz)
   end subroutine get_force


   !> Overwrite ghostpoints to enforce BC at the cylinder
   subroutine apply_ibm()
     use gp_class, only: dirichlet,neumann
     implicit none
     integer :: i,j,k,n
     ! Recompute primitive variables
     call fs%get_primitive()
     ! Reset interior points
     do k=cfg%kmino_,cfg%kmaxo_
        do j=cfg%jmino_,cfg%jmaxo_
           do i=cfg%imino_,cfg%imaxo_
              if (cfg%Gib(i,j,k).lt.0.0_WP) then
                 fs%U(i,j,k)  =u1
                 fs%V(i,j,k)  =0.0_WP
                 fs%W(i,j,k)  =0.0_WP
                 fs%Q(i,j,k,1)=rho1
                 fs%P(i,j,k)  =p1
                 fs%I(i,j,k)  =get_I(rho1,p1)
                 fs%Q(i,j,k,2)=rho1*fs%I(i,j,k)
              end if
           end do
        end do
     end do
     ! Overwrite primitive variables
     call gp%apply_bcond(type=dirichlet,BP=0.0_WP,A=fs%U,dir='U')
     call gp%apply_bcond(type=dirichlet,BP=0.0_WP,A=fs%V,dir='V')
     call gp%apply_bcond(type=dirichlet,BP=0.0_WP,A=fs%W,dir='W')
     call gp%apply_bcond(type=neumann,  BP=0.0_WP,A=fs%P,dir='SC')
     call gp%apply_bcond(type=neumann,  BP=0.0_WP,A=fs%T,dir='SC')
     ! Rebuild conserved quantities in the ghost cells
     do n=1,gp%ngp
        i=gp%gp(n)%ind(1); j=gp%gp(n)%ind(2); k=gp%gp(n)%ind(3)
        fs%Q(i,j,k,1)=get_RHO(fs%T(i,j,k),fs%P(i,j,k))
        fs%I(i,j,k)=get_I(fs%Q(i,j,k,1),fs%P(i,j,k))
        fs%Q(i,j,k,2)=fs%Q(i,j,k,1)*fs%I(i,j,k)
     end do
     call fs%get_momentum()
   end subroutine apply_ibm


   !> Apply boundary conditions
   subroutine apply_bconds()
     implicit none
     integer :: i,j,k

     ! Apply clipped Neumann on primitive variables in x+
     if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) then
        do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           ! Copy over from imax to imax+1 and above
           do i=fs%cfg%imax+1,fs%cfg%imaxo
              ! Copy primitive variables
              fs%Q(i,j,k,1)=fs%Q(fs%cfg%imax,j,k,1)
              fs%P(i,j,k)=fs%P(fs%cfg%imax,j,k)
              fs%I(i,j,k)=fs%I(fs%cfg%imax,j,k)
              fs%U(i,j,k)=max(fs%U(fs%cfg%imax,j,k),0.0_WP)
              fs%V(i,j,k)=fs%V(fs%cfg%imax,j,k)
              fs%W(i,j,k)=fs%W(fs%cfg%imax,j,k)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in y+
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) then
        do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! Copy over from jmax to jmax+1 and above
           do j=fs%cfg%jmax+1,fs%cfg%jmaxo
              ! Copy primitive variables
              fs%Q(i,j,k,1)=fs%Q(i,fs%cfg%jmax,k,1)
              fs%P(i,j,k)=fs%P(i,fs%cfg%jmax,k)
              fs%I(i,j,k)=fs%I(i,fs%cfg%jmax,k)
              fs%U(i,j,k)=fs%U(i,fs%cfg%jmax,k)
              fs%V(i,j,k)=max(fs%V(i,fs%cfg%jmax,k),0.0_WP)
              fs%W(i,j,k)=fs%W(i,fs%cfg%jmax,k)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in y-
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) then
        do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! First copy over V from jmin+1 to jmin
           fs%V(i,fs%cfg%jmin,k)=min(fs%V(i,fs%cfg%jmin+1,k),0.0_WP)
           ! Then copy over from jmin to jmin-1 and below
           do j=fs%cfg%jmino,fs%cfg%jmin-1
              ! Copy primitive variables
              fs%Q(i,j,k,1)=fs%Q(i,fs%cfg%jmin,k,1)
              fs%P(i,j,k)=fs%P(i,fs%cfg%jmin,k)
              fs%I(i,j,k)=fs%I(i,fs%cfg%jmin,k)
              fs%U(i,j,k)=fs%U(i,fs%cfg%jmin,k)
              fs%V(i,j,k)=min(fs%V(i,fs%cfg%jmin,k),0.0_WP)
              fs%W(i,j,k)=fs%W(i,fs%cfg%jmin,k)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in z+
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.fs%cfg%npz) then
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! Copy over from kmax to kmax+1 and above
           do k=fs%cfg%kmax+1,fs%cfg%kmaxo
              ! Copy primitive variables
              fs%Q(i,j,k,1)=fs%Q(i,j,fs%cfg%kmax,1)
              fs%P(i,j,k)=fs%P(i,j,fs%cfg%kmax)
              fs%I(i,j,k)=fs%I(i,j,fs%cfg%kmax)
              fs%U(i,j,k)=fs%U(i,j,fs%cfg%kmax)
              fs%V(i,j,k)=fs%V(i,j,fs%cfg%kmax)
              fs%W(i,j,k)=max(fs%W(i,j,fs%cfg%kmax),0.0_WP)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in z-
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.1) then
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! First copy over W from kmin+1 to kmin
           fs%W(i,j,fs%cfg%kmin)=min(fs%W(i,j,fs%cfg%kmin+1),0.0_WP)
           ! Then copy over from kmin to kmin-1 and below
           do k=fs%cfg%kmino,fs%cfg%kmin-1
              ! Copy primitive variables
              fs%Q(i,j,k,1)=fs%Q(i,j,fs%cfg%kmin,1)
              fs%P(i,j,k)=fs%P(i,j,fs%cfg%kmin)
              fs%I(i,j,k)=fs%I(i,j,fs%cfg%kmin)
              fs%U(i,j,k)=fs%U(i,j,fs%cfg%kmin)
              fs%V(i,j,k)=fs%V(i,j,fs%cfg%kmin)
              fs%W(i,j,k)=min(fs%W(i,j,fs%cfg%kmin),0.0_WP)
           end do
        end do; end do
     end if

     ! Rebuild conserved quantities
     fs%Q(:,:,:,2)=fs%Q(:,:,:,1)*fs%I
     call fs%get_momentum()

   end subroutine apply_bconds


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      
      ! Create compressible flow solver
      create_flow_solver: block
        call fs%initialize(cfg=cfg,name='Compressible NS')
      end block create_flow_solver


      ! Allocate work arrays
      allocate_work_arrays: block
        allocate(dQdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
        allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Ma  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(beta(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc_t(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(div(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize eos and flow parameters
      initialize_parameters: block
        use string,   only: str_long
        use messager, only: log
        use param,    only: param_read
        character(str_long) :: message
        ! Set Pinf to zero
        Pinf=0.0_WP
        ! Read in Gamma
        call param_read('Gamma',Gamma)
        ! Read in Prandtl number
        call param_read('Prandtl number',Prandtl)
        ! Read in shock Mach number and location
        call param_read('Shock Mach number',Ms)
        call param_read('Shock location',Xs)
        ! First generate static shock with normalized pre-shock conditions
        M1=Ms
        rho1=1.0_WP
        rho2=rho1*(Gamma+1.0_WP)*M1**2/((Gamma-1.0_WP)*M1**2+2.0_WP)
        p1=0.25_WP*rho1/Gamma*((Gamma+1.0_WP)*M1/(M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
        p2=p1*(2.0_WP*Gamma/(Gamma+1.0_WP)*(M1**2-1.0_WP)+1.0_WP)
        u1=M1*sqrt(Gamma*p1/rho1)
        u2=u1*rho1/rho2
        ! Now shift frame of reference to obtain moving shock
        u2=abs(u2-u1); M2=u2/sqrt(Gamma*p2/rho2); u1=0.0_WP; M1=u1/sqrt(Gamma*p1/rho1)
        ! Set heat capacities corresponding to a normalized pre-shock
        Cv=(p1+Pinf)/(rho1*(Gamma-1.0_WP))
        ! Get reference temperature
        T0=get_T(rho1,p1)
        ! Viscous parameters
        call param_read('Reynolds number',Re); visc0=rho2*2.0_WP*Rcyl*u2/Re
        ! Output case info
        if (cfg%amRoot) then
           write(message,'("[Gas EOS]               =>  Gamma=",es12.5)')    Gamma; call log(message)
           write(message,'("[Gas EOS]               =>     Cv=",es12.5)')       Cv; call log(message)
           write(message,'("[Shock Mach number]     =>     Ms=",es12.5)')       Ms; call log(message)
           write(message,'("[Pre -shock conditions] =>   rho1=",es12.5)')     rho1; call log(message)
           write(message,'("[Pre -shock conditions] =>     p1=",es12.5)')       p1; call log(message)
           write(message,'("[Pre -shock conditions] =>     u1=",es12.5)')       u1; call log(message)
           write(message,'("[Pre -shock conditions] =>     M1=",es12.5)')       M1; call log(message)
           write(message,'("[Post-shock conditions] =>   rho2=",es12.5)')     rho2; call log(message)
           write(message,'("[Post-shock conditions] =>     p2=",es12.5)')       p2; call log(message)
           write(message,'("[Post-shock conditions] =>     u2=",es12.5)')       u2; call log(message)
           write(message,'("[Post-shock conditions] =>     M2=",es12.5)')       M2; call log(message)
           write(message,'("[Gas Reynolds]          =>     Re=",es12.5)')       Re; call log(message)
           write(message,'("[Gas viscosity]         =>     mu=",es12.5)')    visc0; call log(message)
        end if
      end block initialize_parameters


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
        time=timetracker(amRoot=cfg%amRoot)
        call param_read('Max timestep size',time%dtmax)
        call param_read('Max cfl number',time%cflmax)
        call param_read('Max time',time%tmax)
        time%dt=time%dtmax
        time%itmax=2
      end block initialize_timetracker


      ! Initialize variables
      initialize_variables: block
        integer :: i,j,k
        ! Provide thermodynamic model
        fs%getP=>get_P; fs%getC=>get_C; fs%getS=>get_S; fs%getT=>get_T
        ! Initialize primary variables to normal shock
        do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 fs%U(i,j,k)  =u2*Hshock(Xs-fs%cfg%x(i),delta=0.5_WP*fs%dx)
                 fs%V(i,j,k)  =0.0_WP
                 fs%W(i,j,k)  =0.0_WP
                 fs%Q(i,j,k,1)=rho1+(rho2-rho1)*Hshock(Xs-fs%cfg%xm(i),delta=0.5_WP*fs%dx)
                 fs%P(i,j,k)  =p1  +(p2  -p1  )*Hshock(Xs-fs%cfg%xm(i),delta=0.5_WP*fs%dx)
                 fs%I(i,j,k)  =get_I(fs%Q(i,j,k,1),fs%P(i,j,k))
              end do
           end do
        end do
        ! Initialize conserved variables
        fs%Q(:,:,:,2)=fs%Q(:,:,:,1)*fs%I
        call fs%get_momentum()
        ! Rebuild primitive variables
        call fs%get_primitive()
        ! Interpolate velocity
        call fs%interp_vel(Ui,Vi,Wi)
        ! Compute local Mach number
        Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
        ! Compute dilatation
        call get_div()
      end block initialize_variables


      ! Initialize the ghost points with 3 layers of ghost cells
      create_gp: block
        gp=gpibm(cfg=cfg,no=3)
        call gp%update()
      end block create_gp


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='cylinder')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('T',fs%T)
         call ens_out%add_scalar('Mach',Ma)
         call ens_out%add_scalar('beta',beta)
         call ens_out%add_scalar('visc',visc)
         call ens_out%add_scalar('visc_t',visc_t)
         call ens_out%add_scalar('div',div) 
         call ens_out%add_scalar('Gib',cfg%Gib)
         call ens_out%add_scalar('IBM',gp%label)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create monitor files
      create_monitor: block
        !> Perform and output monitoring
        call fs%get_info()
        call get_force()
        ! Create simulation monitor
        mfile=monitor(fs%cfg%amRoot,'simulation')
        call mfile%add_column(time%n,'Timestep number')
        call mfile%add_column(time%t,'Time')
        call mfile%add_column(time%dt,'Timestep size')
        call mfile%add_column(time%cfl,'Maximum CFL')
        call mfile%add_column(fs%Umax,'Umax')
        call mfile%add_column(fs%Vmax,'Vmax')
        call mfile%add_column(fs%Wmax,'Wmax')
        call mfile%add_column(fs%RHOmax,'max(RHO)')
        call mfile%add_column(fs%RHOmin,'min(RHO)')
        call mfile%add_column(fs%Imax  ,'max(I)'  )
        call mfile%add_column(fs%Imin  ,'min(I)'  )
        call mfile%add_column(fs%Pmax  ,'max(P)'  )
        call mfile%add_column(fs%Pmin  ,'min(P)'  )
        call mfile%add_column(fs%Tmax  ,'max(T)'  )
        call mfile%add_column(fs%Tmin  ,'min(T)'  )
        call mfile%write()
        ! Create CFL monitor
        cflfile=monitor(fs%cfg%amRoot,'cfl')
        call cflfile%add_column(time%n,'Timestep number')
        call cflfile%add_column(time%t,'Time')
        call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
        call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
        call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
        call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
        call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
        call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
        call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
        call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
        call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
        call cflfile%write()
        ! Create conservation monitor
        consfile=monitor(fs%cfg%amRoot,'conservation')
        call consfile%add_column(time%n,'Timestep number')
        call consfile%add_column(time%t,'Time')
        call consfile%add_column(fs%Qint(1),'Mass')
        call consfile%add_column(fs%Qint(2),'Energy')
        call consfile%add_column(fs%Qint(3),'U Momentum')
        call consfile%add_column(fs%Qint(4),'V Momentum')
        call consfile%add_column(fs%Qint(5),'W Momentum')
        call consfile%add_column(fs%RHOKint,'Kinetic Energy')
        call consfile%add_column(fs%RHOSint,'Entropy')
        call consfile%write()
        ! Create IBM monitor
        consfile=monitor(fs%cfg%amRoot,'ibm')
        call consfile%add_column(time%n,'Timestep number')
        call consfile%add_column(time%t,'Time')
        call consfile%add_column(ibm_force(1),'X Force')
        call consfile%add_column(ibm_force(2),'Y Force')
        call consfile%add_column(ibm_force(3),'Z Force')
        call ibmfile%write()
      end block create_monitor


      ! Create a timing monitor
      create_timing: block
        ! Create timers
        tstep=timer(comm=cfg%comm,name='Total')
        tibm =timer(comm=cfg%comm,name='IBM')
        tcom =timer(comm=cfg%comm,name='Comp')
        ! Create corresponding monitor file
        timefile=monitor(cfg%amRoot,'timing')
        call timefile%add_column(time%n,'Timestep number')
        call timefile%add_column(time%t,'Time')
        call timefile%add_column(tstep%time,trim(tstep%name))
        call timefile%add_column(tibm%time,trim(tibm%name))
        call timefile%add_column(tcom%time,trim(tcom%name))
      end block create_timing

    end subroutine simulation_init


    !> Perform an NGA2 simulation
    subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Reset all timers and start timestep timer
         call tstep%reset()
         call tibm%reset()
         call tcom%reset()
         call tstep%start()

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember conserved variables
         fs%Qold=fs%Q

         ! Prepare SGS viscosity models
         call prepare_viscosities()


         ! First RK step ====================================================================================
         ! Get non-SL RHS and increment
         call tcom%start() ! Start compressible timer
         call fs%rhs(dQdt(:,:,:,:,1))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         call tcom%stop() ! Stop compressible timer
         ! Apply IBM
         call tibm%start() ! Start IBM timer
         call apply_ibm()
         call tibm%stop() ! Stop IBM timer

         ! Second RK step ===================================================================================
         ! Get non-SL RHS and increment
         call tcom%start() ! Start compressible timer
         call fs%rhs(dQdt(:,:,:,:,2))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         call tcom%stop() ! Stop compressible timer
         ! Apply IBM
         call tibm%start() ! Start IBM timer
         call apply_ibm()
         call tibm%stop() ! Stop IBM timer

         ! Third RK step ====================================================================================
         ! Get non-SL RHS and increment
         call tcom%start() ! Start compressible timer
         call fs%rhs(dQdt(:,:,:,:,3))
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         call tcom%stop() ! Stop compressible timer
         ! Apply IBM
         call tibm%start() ! Start IBM timer
         call apply_ibm()
         call tibm%stop() ! Stop IBM timer

         ! Fourth RK step ===================================================================================
         ! Get non-SL RHS and increment
         call tcom%start() ! Start compressible timer
         call fs%rhs(dQdt(:,:,:,:,4))
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         call tcom%stop() ! Stop compressible timer
         ! Apply IBM
         call tibm%start() ! Start IBM timer
         call apply_ibm()
         call tibm%stop() ! Stop IBM timer

         ! Apply boundary conditions
         call apply_bconds()

         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)

         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C

         ! Compute dilatation
         call get_div()

         ! Stop timestep timer
         call tstep%stop()

         !> Perform and output monitoring
         call fs%get_info()
         call get_force()
         call mfile%write()
         call timefile%write()
         call cflfile%write()
         call consfile%write()
         call ibmfile%write()

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

      end do

 end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(dQdt,Ui,Vi,Wi,Ma,beta,visc,visc_t,div)
      call tstep%finalize()
      call tibm%finalize()
      call tcom%finalize()
      
   end subroutine simulation_final
   
   
end module simulation
