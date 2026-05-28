!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use spcomp_class,      only: spcomp
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(spcomp),      public :: fs
   type(lpt),         public :: lp
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile,lptfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc,visc_t,div
   real(WP), dimension(:,:,:)    , allocatable :: srcUlp,srcVlp,srcWlp,srcIlp
   real(WP), dimension(:,:,:)    , allocatable :: stressx,stressy,stressz,stressI

   !> Post-shock viscosity and temperature
   real(WP) :: visc0,T0

   !> Equations of state
   real(WP) :: Pinf,Gamma,Cv,Prandtl

   !> Flow parameters
   real(WP) :: Ms,Xs
   real(WP) :: rho1,p1,u1,M1
   real(WP) :: rho2,p2,u2,M2
   real(WP) :: Re
   
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
              lp%VF(i,j,k)=lp%VF(fs%cfg%imax,j,k)
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
              lp%VF(i,j,k)=lp%VF(i,fs%cfg%jmax,k)
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
              lp%VF(i,j,k)=lp%VF(i,fs%cfg%jmin,k)
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
              lp%VF(i,j,k)=lp%VF(i,j,fs%cfg%kmax)
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
              lp%VF(i,j,k)=lp%VF(i,j,fs%cfg%kmin)
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
        allocate(dQdt   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
        allocate(Ui     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vi     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wi     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Ma     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(beta   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc_t(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(div    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(stressI(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcUlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcVlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcWlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(srcIlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
        call param_read('Reynolds number',Re); visc0=rho2*1.0_WP*u2/Re
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


      ! Initialize our LPT solver
      initialize_lpt: block
        use random, only: random_lognormal,random_uniform
        use mathtools, only: Pi,twoPi
        use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_INTEGER
        use parallel, only: MPI_REAL_WP
        real(WP) :: VFavg,Vol_,sumVolp,dp,Wc,Xc,Tp
        integer :: i,j,k,ii,jj,kk,nn,ip,jp,kp,np,offset,ierr
        integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
        integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
        logical :: overlap
        ! Create solver
        lp=lpt(cfg=cfg,name='LPT')
        ! Get mean volume fraction from input
        call param_read('Particle volume fraction',VFavg)
        ! Get drag model from input
        ! Get particle density from input
        call param_read('Particle density',lp%rho)
        ! Get particle specific heat from input
        call param_read('Particle density',lp%Cp)
        ! Get particle diameter from input
        call param_read('Particle diameter',dp)
        ! Get particle temperature from input
        call param_read('Particle temperature',Tp)
        ! Set collision timescale
        call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
        ! Set coefficient of restitution
        call param_read('Coefficient of restitution',lp%e_n)
        call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
        ! Set filter scale to 3.5*dx
        lp%filter_width=3.5_WP*cfg%min_meshsize
        ! Get curtain properties
        call param_read('Curtain width',Wc)
        call param_read('Curtain position',Xc)
        ! Initialize particles
        ! Get volume of domain belonging to this proc
        Vol_=0.0_WP
        do k=lp%cfg%kmin_,lp%cfg%kmax_
           do j=lp%cfg%jmin_,lp%cfg%jmax_
              do i=lp%cfg%imin_,fs%cfg%imax_
                 if (lp%cfg%x(i).gt.Xc.and.lp%cfg%x(i).lt.Xc+Wc) Vol_=Vol_+lp%cfg%dx(i)*lp%cfg%dy(j)*lp%cfg%dz(k)
              end do
           end do
        end do
        ! Get particle diameters
        np=ceiling(VFavg*Vol_/(pi*dp**3/6.0_WP))
        call lp%resize(np)
        ! Allocate particle in cell arrays
        allocate(npic(     lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); npic=0
        allocate(ipic(1:40,lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); ipic=0
        ! Distribute particles
        sumVolp=0.0_WP
        do i=1,np
           ! Set the diameter
           lp%p(i)%d=dp
           ! Give position (avoid overlap)
           overlap=.true.
           do while (overlap)
              lp%p(i)%pos=[random_uniform(max(Xc,lp%cfg%x(lp%cfg%imin_)),min(Xc+Wc,lp%cfg%x(lp%cfg%imax_)-dp)),&
                   &       random_uniform(lp%cfg%y(lp%cfg%jmin_),lp%cfg%y(lp%cfg%jmax_+1)-dp),&
                   &       random_uniform(lp%cfg%z(lp%cfg%kmin_),lp%cfg%z(lp%cfg%kmax_+1)-dp)]
              if (lp%cfg%nz.eq.1) lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin_)
              lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
              overlap=.false.
              do kk=lp%p(i)%ind(3)-1,lp%p(i)%ind(3)+1
                 do jj=lp%p(i)%ind(2)-1,lp%p(i)%ind(2)+1
                    do ii=lp%p(i)%ind(1)-1,lp%p(i)%ind(1)+1
                       do nn=1,npic(ii,jj,kk)
                          j=ipic(nn,ii,jj,kk)
                          if (sqrt(sum((lp%p(i)%pos-lp%p(j)%pos)**2)).lt.0.5_WP*(lp%p(i)%d+lp%p(j)%d)) overlap=.true.
                       end do
                    end do
                 end do
              end do
           end do

           ! Activate the particle
           lp%p(i)%flag=0
           ip=lp%p(i)%ind(1); jp=lp%p(i)%ind(2); kp=lp%p(i)%ind(3)
           npic(ip,jp,kp)=npic(ip,jp,kp)+1
           ipic(npic(ip,jp,kp),ip,jp,kp)=i
           ! Give temperature
           lp%p(i)%T=Tp
           ! Give zero velocity
           lp%p(i)%vel=0.0_WP
           ! Give zero collision force
           lp%p(i)%Acol=0.0_WP
           lp%p(i)%Tcol=0.0_WP
           ! Sum up volume
           sumVolp=sumVolp+Pi/6.0_WP*lp%p(i)%d**3
        end do
        deallocate(npic,ipic)
        call lp%sync()
        ! Set ID
        offset=0
        do i=1,lp%cfg%rank
           offset=offset+lp%np_proc(i)
        end do
        do i=1,lp%np_
           lp%p(i)%id=int(i+offset,8)
        end do
        ! Get mean diameter and volume fraction
        call MPI_ALLREDUCE(sumVolp,VFavg,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); VFavg=VFavg/(Wc*lp%cfg%yL*lp%cfg%zL)
        if (lp%cfg%amRoot) then
           print*,"===== Particle Setup Description ====="
           print*,'Number of particles', lp%np
           print*,'Mean volume fraction',VFavg
        end if
        ! Get initial particle volume fraction
        call lp%update_VF()
      end block initialize_lpt


      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
        integer :: i
        pmesh=partmesh(nvar=2,nvec=1,name='lpt')
        pmesh%varname(1)='diameter'
        pmesh%varname(2)='temperature'
        pmesh%vecname(1)='velocity'
        call lp%update_partmesh(pmesh)
        do i=1,lp%np_
           pmesh%var(1,i)=lp%p(i)%d
           pmesh%var(2,i)=lp%p(i)%T
           pmesh%vec(:,1,i)=lp%p(i)%vel
        end do
      end block create_pmesh


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
        ! Multiply by volume fraction
        do i=1,fs%nQ
           fs%Q(:,:,:,i)=fs%Q(:,:,:,i)*(1.0_WP-lp%VF)
        end do
        call fs%get_momentum()
        ! Rebuild primitive variables
        call fs%get_primitive(1.0_WP-lp%VF)
        ! Interpolate velocity
        call fs%interp_vel(Ui,Vi,Wi)
        ! Compute local Mach number
        Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
        ! Compute dilatation
        call get_div()
        !> Perform and output monitoring
        call fs%get_info()
        call mfile%write()
        call cflfile%write()
        call consfile%write()
      end block initialize_variables


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='shock')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('T',fs%T)
         call ens_out%add_scalar('Mach',Ma)
         call ens_out%add_scalar('beta',beta)
         call ens_out%add_scalar('visc',visc)
         call ens_out%add_scalar('visc_t',visc_t)
         call ens_out%add_scalar('div',div) 
         call ens_out%add_scalar('epsp',lp%VF)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create monitor files
      create_monitor: block
        real(WP) :: cfl
        call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
        call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
        call fs%get_info()
        call lp%get_max()
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
        ! Create LPT monitor
        lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
        call lptfile%add_column(time%n,'Timestep number')
        call lptfile%add_column(time%t,'Time')
        call lptfile%add_column(lp%np,'Particle number')
        call lptfile%add_column(lp%VFmean,'mean(VFp)')
        call lptfile%add_column(lp%VFmax,'max(VFp)')
        call lptfile%add_column(lp%Umin,'min(U)')
        call lptfile%add_column(lp%Umax,'max(U)')
        call lptfile%add_column(lp%Vmin,'min(V)')
        call lptfile%add_column(lp%Vmax,'max(V)')
        call lptfile%add_column(lp%Wmin,'min(W)')
        call lptfile%add_column(lp%Wmax,'max(W)')
        call lptfile%add_column(lp%Remax,'max(Re)')
        call lptfile%add_column(lp%Mamax,'max(Ma)')
        call lptfile%add_column(lp%Knmax,'max(Kn)')
        call lptfile%write()
      end block create_monitor

    end subroutine simulation_init


    !> Perform an NGA2 simulation
    subroutine simulation_run
      implicit none
      real(WP) :: cfl

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember conserved variables
         fs%Qold=fs%Q

         ! Prepare SGS viscosity models
         call prepare_viscosities()

         ! First RK step ====================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=1,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,1))
         ! LPT source
         dQdt(:,:,:,2,1)=dQdt(:,:,:,2,1)+srcIlp
         dQdt(:,:,:,3,1)=dQdt(:,:,:,3,1)+srcUlp
         dQdt(:,:,:,4,1)=dQdt(:,:,:,4,1)+srcVlp
         dQdt(:,:,:,5,1)=dQdt(:,:,:,5,1)+srcWlp
         ! Advance
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Recompute primitive variables
         call fs%get_primitive(1.0_WP-lp%VF)

         ! Second RK step ===================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=2,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,2))
         ! LPT source
         dQdt(:,:,:,2,1)=dQdt(:,:,:,2,2)+srcIlp
         dQdt(:,:,:,3,1)=dQdt(:,:,:,3,2)+srcUlp
         dQdt(:,:,:,4,1)=dQdt(:,:,:,4,2)+srcVlp
         dQdt(:,:,:,5,1)=dQdt(:,:,:,5,2)+srcWlp
         ! Advance
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Recompute primitive variables
         call fs%get_primitive(1.0_WP-lp%VF)

         ! Third RK step ====================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=3,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,3))
         ! LPT source
         dQdt(:,:,:,2,1)=dQdt(:,:,:,2,3)+srcIlp
         dQdt(:,:,:,3,1)=dQdt(:,:,:,3,3)+srcUlp
         dQdt(:,:,:,4,1)=dQdt(:,:,:,4,3)+srcVlp
         dQdt(:,:,:,5,1)=dQdt(:,:,:,5,3)+srcWlp
         ! Advance
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         ! Recompute primitive variables
         call fs%get_primitive(1.0_WP-lp%VF)

         ! Fourth RK step ===================================================================================
         ! Particle increment
         call fs%get_div_stress(stressx,stressy,stressz,stressI)
         call lp%substep_rk4(stage=4,dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=fs%T,C=fs%C,&
              stress_x=stressx,stress_y=stressy,stress_z=stressz,heat_flux=stressI,srcU=srcUlp,srcV=srcVlp,srcW=srcWlp,srcI=srcIlp)
         ! Get non-SL RHS and increment
         call fs%rhs(VF=lp%VF,VFU=lp%VFU,VFV=lp%VFV,VFW=lp%VFW,dQdt=dQdt(:,:,:,:,4))
         ! LPT source
         dQdt(:,:,:,2,1)=dQdt(:,:,:,2,4)+srcIlp
         dQdt(:,:,:,3,1)=dQdt(:,:,:,3,4)+srcUlp
         dQdt(:,:,:,4,1)=dQdt(:,:,:,4,4)+srcVlp
         dQdt(:,:,:,5,1)=dQdt(:,:,:,5,4)+srcWlp
         ! Advance
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Recompute primitive variables
         call fs%get_primitive(1.0_WP-lp%VF)

         ! Apply boundary conditions
         call apply_bconds()

         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)

         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C

         ! Compute dilatation
         call get_div()

         !> Perform and output monitoring
         call fs%get_info()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         call lptfile%write()

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              integer :: i
              call lp%update_partmesh(pmesh)
              do i=1,lp%np_
                 pmesh%var(1,i)=lp%p(i)%d
                 pmesh%var(2,i)=lp%p(i)%T
                 pmesh%vec(:,1,i)=lp%p(i)%vel
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

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
      deallocate(dQdt,Ui,Vi,Wi,Ma,beta,visc,visc_t,div,srcUlp,srcVlp,srcWlp,srcIlp,stressx,stressy,stressz,stressI)
      
   end subroutine simulation_final
   
   
end module simulation
