!> Mie-Gruneisen EOS on a us-up Hugoniot reference (pure substance, ns=1)
!>    p(rho,e) = p_ref(rho) + rho*Gamma(rho)*(e - q - e_ref(rho))
!> with the standard Gruneisen closure rho*Gamma(rho) = rho0*gamma0 (constant),
!> and reference curves in eta = 1 - rho0/rho:
!>  - compression (eta>0): Hugoniot from the (up to) cubic us-up fit
!>       us = c0/D(eta), D(eta) = 1 - s1*eta - s2*eta^2 - s3*eta^3   [eta = up/us]
!>       p_ref = rho0*c0^2*eta/D^2 , e_ref = p_ref*eta/(2*rho0)
!>    (linear us = c0 + s1*up is the s2=s3=0 special case)
!>  - tension (eta<0): linear extrapolation p_ref = c0^2*(rho-rho0), e_ref = 0
!> Knee clamp (NASG brhomax analog): the reference is frozen past etamax, where D
!> falls below Dmin, or where dp_ref/deta = rho0*c0^2*(1+s1*eta+3*s2*eta^2+5*s3*eta^3)/D^3
!> loses positivity (a fitted cubic is only trustworthy on its data range; beyond it
!> it can turn non-monotone without ever reaching D=Dmin) -- bounded monotone
!> continuation: p keeps rising through the energy term, c^2 = p*Gamma/rho stays positive.
!> Temperature: entropy-consistent constant-cv completion. T_ref(eta) is integrated
!> at initialize from  dT/deta = gamma0*T + B(eta)/cv  with B = de_ref/deta - p_ref/rho0
!> (a Hugoniot reference is NOT an isentrope, so the isentrope law T0*exp(gamma0*eta)
!> would leave dS inexact; this ODE is the unique constant-cv completion for which
!> dS=(de+p dv)/T is exact), then T(rho,e) = T_ref + (e-q-e_ref)/cv and the entropy
!> has the closed form s = cv*(ln(T) + gamma0*rho0/rho) + qp.
!> Set the clamp/table fields (etamin,etamax,Dmin,ntab) BEFORE initialize -- they are
!> consumed when the T_ref table is built.
module mie_gruneisen_class
   use precision,      only: WP
   use material_class, only: material
   implicit none
   private

   public :: mie_gruneisen

   type, extends(material) :: mie_gruneisen
      real(WP) :: rho0  =0.0_WP   !< Reference density (us-up fit anchor)
      real(WP) :: c0    =0.0_WP   !< Bulk sound speed at the reference state
      real(WP) :: s1    =0.0_WP   !< us-up fit coefficients: D(eta)=1-s1*eta-s2*eta^2-s3*eta^3
      real(WP) :: s2    =0.0_WP
      real(WP) :: s3    =0.0_WP
      real(WP) :: gamma0=0.0_WP   !< Gruneisen coefficient at rho0 (Gamma(rho)=gamma0*rho0/rho)
      real(WP) :: cv    =0.0_WP   !< Constant specific heat (temperature completion)
      real(WP) :: T0    =0.0_WP   !< Temperature at the reference state (rho0, p=0, e=q)
      real(WP) :: q     =0.0_WP   !< Energy of formation
      real(WP) :: qp    =0.0_WP   !< Entropy of formation
      ! Clamp/table controls -- consumed at initialize when the T_ref table is built
      real(WP) :: etamin=-1.0_WP  !< T_ref table lower bound (tension side; frozen below)
      real(WP) :: etamax=0.95_WP  !< Compression knee: reference frozen past min(etamax, D=Dmin)
      real(WP) :: Dmin  =1.0e-3_WP!< Floor on D(eta) defining the knee (inert for sane fits)
      integer  :: ntab  =2000     !< T_ref table resolution
      ! Set at initialize
      real(WP) :: gr0m  =0.0_WP   !< rho0*gamma0 = rho*Gamma(rho) = (dp/de)_rho, constant
      real(WP) :: etaknee=0.0_WP  !< Effective knee = min(etamax, first eta with D<=Dmin)
      real(WP) :: pr_knee=0.0_WP  !< Frozen reference pressure past the knee
      real(WP) :: er_knee=0.0_WP  !< Frozen reference energy past the knee
      real(WP) :: deta   =0.0_WP  !< T_ref table spacing (etamin snapped to hit eta=0 exactly)
      real(WP), dimension(:), allocatable :: Ttab   !< T_ref(eta) on the uniform table
   contains
      procedure, private :: mg_initialize
      generic   :: initialize              => mg_initialize
      procedure, private :: get_ref        => mg_get_ref
      procedure, private :: get_ref_d      => mg_get_ref_d
      procedure, private :: get_Tref       => mg_get_Tref
      procedure, private :: ode_rhs        => mg_ode_rhs
      procedure :: get_p_from_rho_e        => mg_get_p_from_rho_e
      procedure :: get_T_from_p_rho        => mg_get_T_from_p_rho
      procedure :: get_c_from_p_rho        => mg_get_c_from_p_rho
      procedure :: get_T_from_rho_e        => mg_get_T_from_rho_e
      procedure :: get_c_from_rho_e        => mg_get_c_from_rho_e
      procedure :: get_cv_from_rho_e       => mg_get_cv_from_rho_e
      procedure :: get_cv_from_rho_T       => mg_get_cv_from_rho_T
      procedure :: get_e_from_p_rho        => mg_get_e_from_p_rho
      procedure :: get_e_from_p_T          => mg_get_e_from_p_T
      procedure :: get_p_from_rho_T        => mg_get_p_from_rho_T
      procedure :: get_rho_from_p_T        => mg_get_rho_from_p_T
      procedure :: get_h_from_p_T          => mg_get_h_from_p_T
      procedure :: get_hk_from_p_T         => mg_get_hk_from_p_T
      procedure :: get_s_from_p_T          => mg_get_s_from_p_T
      procedure :: get_g_from_p_T          => mg_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=> mg_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho     => mg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => mg_get_rhoe_from_p_T
      procedure :: print                   => mg_print
      procedure :: finalize                => mg_finalize
   end type mie_gruneisen

contains

   subroutine mg_initialize(this,rho0,c0,s1,s2,s3,gamma0,cv,T0,q,qp,name)
      use messager, only: die
      class(mie_gruneisen), intent(inout) :: this
      real(WP), intent(in) :: rho0,c0,s1,s2,s3,gamma0,cv,T0,q,qp
      character(len=*), intent(in), optional :: name
      real(WP) :: eta,D,T,k1,k2,k3,k4
      integer :: i,i0
      if (rho0.le.0.0_WP.or.c0.le.0.0_WP.or.cv.le.0.0_WP.or.T0.le.0.0_WP) call die('[mie_gruneisen initialize] rho0, c0, cv, T0 must be positive')
      if (this%etamin.ge.0.0_WP.or.this%etamax.le.0.0_WP.or.this%ntab.lt.3) call die('[mie_gruneisen initialize] need etamin<0<etamax and ntab>=3')
      this%rho0=rho0; this%c0=c0; this%s1=s1; this%s2=s2; this%s3=s3
      this%gamma0=gamma0; this%cv=cv; this%T0=T0; this%q=q; this%qp=qp
      this%gr0m=gamma0*rho0
      this%ns=1
      if (present(name)) this%name=name
      allocate(this%species_names(1)); this%species_names(1)=this%name
      ! Snap the table so eta=0 is a grid point (T_ref(0)=T0 exact), spacing from requested bounds
      this%deta=(this%etamax-this%etamin)/real(this%ntab-1,WP)
      i0=1+nint(-this%etamin/this%deta)
      this%etamin=-real(i0-1,WP)*this%deta
      this%etamax=this%etamin+real(this%ntab-1,WP)*this%deta
      ! Locate the knee: first grid eta>0 where D(eta)<=Dmin or dp_ref/deta<=0
      ! (sign of 1+s1*eta+3*s2*eta^2+5*s3*eta^3), capped at etamax
      this%etaknee=this%etamax
      do i=i0+1,this%ntab
         eta=this%etamin+real(i-1,WP)*this%deta
         D=1.0_WP-this%s1*eta-this%s2*eta**2-this%s3*eta**3
         if (D.le.this%Dmin.or.1.0_WP+this%s1*eta+3.0_WP*this%s2*eta**2+5.0_WP*this%s3*eta**3.le.0.0_WP) then
            this%etaknee=this%etamin+real(i-2,WP)*this%deta
            exit
         end if
      end do
      ! Frozen reference past the knee
      call this%get_ref(rho=this%rho0/(1.0_WP-this%etaknee),pr=this%pr_knee,er=this%er_knee,frozen=.false.)
      ! Build the T_ref table: RK4 on dT/deta = gamma0*T + B(eta)/cv from (eta=0, T=T0),
      ! upward to the knee (held constant beyond) and downward to etamin
      if (allocated(this%Ttab)) deallocate(this%Ttab)
      allocate(this%Ttab(this%ntab))
      this%Ttab(i0)=this%T0
      T=this%T0
      do i=i0+1,this%ntab
         eta=this%etamin+real(i-2,WP)*this%deta
         k1=this%ode_rhs(eta                ,T                       )
         k2=this%ode_rhs(eta+0.5_WP*this%deta,T+0.5_WP*this%deta*k1)
         k3=this%ode_rhs(eta+0.5_WP*this%deta,T+0.5_WP*this%deta*k2)
         k4=this%ode_rhs(eta+this%deta      ,T+this%deta*k3        )
         T=T+this%deta*(k1+2.0_WP*k2+2.0_WP*k3+k4)/6.0_WP
         this%Ttab(i)=T
      end do
      T=this%T0
      do i=i0-1,1,-1
         eta=this%etamin+real(i,WP)*this%deta
         k1=this%ode_rhs(eta                ,T                       )
         k2=this%ode_rhs(eta-0.5_WP*this%deta,T-0.5_WP*this%deta*k1)
         k3=this%ode_rhs(eta-0.5_WP*this%deta,T-0.5_WP*this%deta*k2)
         k4=this%ode_rhs(eta-this%deta      ,T-this%deta*k3        )
         T=T-this%deta*(k1+2.0_WP*k2+2.0_WP*k3+k4)/6.0_WP
         this%Ttab(i)=T
      end do
   end subroutine mg_initialize

   !> Completion-ODE right-hand side dT/deta = gamma0*T + B(eta)/cv,
   !> B = de_ref/deta - p_ref/rho0; held at zero past the knee (frozen reference)
   real(WP) function mg_ode_rhs(this,eta,T) result(f)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: eta,T
      real(WP) :: D,Dp,pH,dpH,B
      if (eta.gt.this%etaknee) then
         f=0.0_WP
      else if (eta.gt.0.0_WP) then
         D  =max(1.0_WP-this%s1*eta-this%s2*eta**2-this%s3*eta**3,this%Dmin)
         Dp =-(this%s1+2.0_WP*this%s2*eta+3.0_WP*this%s3*eta**2)
         pH =this%rho0*this%c0**2*eta/D**2
         dpH=this%rho0*this%c0**2*(D-2.0_WP*eta*Dp)/D**3
         B  =(dpH*eta-pH)/(2.0_WP*this%rho0)
         f  =this%gamma0*T+B/this%cv
      else
         B  =-this%c0**2*eta/(1.0_WP-eta)
         f  =this%gamma0*T+B/this%cv
      end if
   end function mg_ode_rhs

   !> Reference curve (p_ref,e_ref) at rho; frozen past the knee unless frozen=.false.
   subroutine mg_get_ref(this,rho,pr,er,frozen)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in)  :: rho
      real(WP), intent(out) :: pr,er
      logical, intent(in), optional :: frozen
      real(WP) :: eta,D
      logical :: frz
      frz=.true.; if (present(frozen)) frz=frozen
      eta=1.0_WP-this%rho0/rho
      if (frz.and.eta.gt.this%etaknee) then
         pr=this%pr_knee; er=this%er_knee
      else if (eta.gt.0.0_WP) then
         D =max(1.0_WP-this%s1*eta-this%s2*eta**2-this%s3*eta**3,this%Dmin)
         pr=this%rho0*this%c0**2*eta/D**2
         er=pr*eta/(2.0_WP*this%rho0)
      else
         pr=this%c0**2*(rho-this%rho0)
         er=0.0_WP
      end if
   end subroutine mg_get_ref

   !> Reference curve and its rho-derivatives (zero past the knee: frozen continuation)
   subroutine mg_get_ref_d(this,rho,pr,er,dpr,der)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in)  :: rho
      real(WP), intent(out) :: pr,er,dpr,der
      real(WP) :: eta,D,Dp,dpH,detadrho
      eta=1.0_WP-this%rho0/rho
      if (eta.gt.this%etaknee) then
         pr=this%pr_knee; er=this%er_knee; dpr=0.0_WP; der=0.0_WP
      else if (eta.gt.0.0_WP) then
         D  =max(1.0_WP-this%s1*eta-this%s2*eta**2-this%s3*eta**3,this%Dmin)
         Dp =-(this%s1+2.0_WP*this%s2*eta+3.0_WP*this%s3*eta**2)
         pr =this%rho0*this%c0**2*eta/D**2
         er =pr*eta/(2.0_WP*this%rho0)
         detadrho=this%rho0/rho**2
         dpH=this%rho0*this%c0**2*(D-2.0_WP*eta*Dp)/D**3
         dpr=dpH*detadrho
         der=(dpH*eta+pr)/(2.0_WP*this%rho0)*detadrho
      else
         pr =this%c0**2*(rho-this%rho0)
         er =0.0_WP
         dpr=this%c0**2
         der=0.0_WP
      end if
   end subroutine mg_get_ref_d

   !> T_ref(eta) by linear interpolation on the init-time table (clamped at both ends)
   real(WP) function mg_get_Tref(this,rho) result(Tr)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho
      real(WP) :: eta,w
      integer :: i
      eta=min(max(1.0_WP-this%rho0/rho,this%etamin),this%etamax)
      w=(eta-this%etamin)/this%deta
      i=min(int(w)+1,this%ntab-1)
      w=w-real(i-1,WP)
      Tr=(1.0_WP-w)*this%Ttab(i)+w*this%Ttab(i+1)
   end function mg_get_Tref

   real(WP) function mg_get_p_from_rho_e(this,rho,e,y) result(p)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er
      call this%get_ref(rho=rho,pr=pr,er=er)
      p=pr+this%gr0m*(e-this%q-er)
   end function mg_get_p_from_rho_e

   real(WP) function mg_get_T_from_p_rho(this,p,rho,y) result(T)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er
      call this%get_ref(rho=rho,pr=pr,er=er)
      T=this%get_Tref(rho)+(p-pr)/(this%gr0m*this%cv)
   end function mg_get_T_from_p_rho

   !> c^2 = (dp/drho)_e + p*Gamma/rho = dp_ref - Gamma*rho*de_ref + p*gr0m/rho^2
   real(WP) function mg_get_c_from_p_rho(this,p,rho,y) result(c)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er,dpr,der
      call this%get_ref_d(rho=rho,pr=pr,er=er,dpr=dpr,der=der)
      c=sqrt(max(0.0_WP,dpr-this%gr0m*der+p*this%gr0m/rho**2))
   end function mg_get_c_from_p_rho

   !> Optimal (rho,e) primitives: one reference evaluation each
   real(WP) function mg_get_T_from_rho_e(this,rho,e,y) result(T)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er
      call this%get_ref(rho=rho,pr=pr,er=er)
      T=this%get_Tref(rho)+(e-this%q-er)/this%cv
   end function mg_get_T_from_rho_e

   real(WP) function mg_get_c_from_rho_e(this,rho,e,y) result(c)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er,dpr,der,p
      call this%get_ref_d(rho=rho,pr=pr,er=er,dpr=dpr,der=der)
      p=pr+this%gr0m*(e-this%q-er)
      c=sqrt(max(0.0_WP,dpr-this%gr0m*der+p*this%gr0m/rho**2))
   end function mg_get_c_from_rho_e

   real(WP) function mg_get_cv_from_rho_e(this,rho,e,y) result(cv)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      cv=this%cv
   end function mg_get_cv_from_rho_e

   real(WP) function mg_get_cv_from_rho_T(this,rho,T,y) result(cv)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      cv=this%cv
   end function mg_get_cv_from_rho_T

   real(WP) function mg_get_e_from_p_rho(this,p,rho,y) result(e)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er
      call this%get_ref(rho=rho,pr=pr,er=er)
      e=this%q+er+(p-pr)/this%gr0m
   end function mg_get_e_from_p_rho

   real(WP) function mg_get_p_from_rho_T(this,rho,T,y) result(p)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr,er
      call this%get_ref(rho=rho,pr=pr,er=er)
      p=pr+this%gr0m*this%cv*(T-this%get_Tref(rho))
   end function mg_get_p_from_rho_T

   !> Invert p(rho,T) for rho: p is monotone increasing in rho up to the knee (then flat),
   !> so a bisection-safeguarded Newton on [near-vacuum, knee] is total; out-of-range
   !> targets return the clamped endpoint (consistent with the frozen-reference policy)
   real(WP) function mg_get_rho_from_p_T(this,p,T,y) result(rho)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rlo,rhi,flo,fhi,f,df,pr,er,dpr,der,dTr,eta
      integer :: it
      integer, parameter :: itmax=100
      real(WP), parameter :: rtol=1.0e-14_WP
      rlo=1.0e-9_WP*this%rho0
      rhi=this%rho0/(1.0_WP-this%etaknee)
      flo=this%get_p_from_rho_T(rho=rlo,T=T,y=y)-p
      fhi=this%get_p_from_rho_T(rho=rhi,T=T,y=y)-p
      if (flo.ge.0.0_WP) then; rho=rlo; return; end if
      if (fhi.le.0.0_WP) then; rho=rhi; return; end if
      rho=this%rho0
      do it=1,itmax
         call this%get_ref_d(rho=rho,pr=pr,er=er,dpr=dpr,der=der)
         f=pr+this%gr0m*this%cv*(T-this%get_Tref(rho))-p
         if (f.gt.0.0_WP) then; rhi=rho; else; rlo=rho; end if
         if (rhi-rlo.lt.rtol*this%rho0) exit
         ! dp/drho at fixed T; dT_ref/drho from the completion ODE (zero where the table is frozen)
         eta=1.0_WP-this%rho0/rho
         dTr=0.0_WP
         if (eta.gt.this%etamin.and.eta.lt.this%etaknee) dTr=this%ode_rhs(eta,this%get_Tref(rho))*this%rho0/rho**2
         df=dpr-this%gr0m*this%cv*dTr
         if (abs(df).gt.tiny(1.0_WP).and.abs(f).lt.0.5_WP*abs(df)*(rhi-rlo)) then
            rho=rho-f/df
            if (rho.le.rlo.or.rho.ge.rhi) rho=0.5_WP*(rlo+rhi)   ! Newton left the bracket -> bisect
         else
            rho=0.5_WP*(rlo+rhi)
         end if
      end do
   end function mg_get_rho_from_p_T

   real(WP) function mg_get_e_from_p_T(this,p,T,y) result(e)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      rho=this%get_rho_from_p_T(p=p,T=T,y=y)
      e=this%get_e_from_p_rho(p=p,rho=rho,y=y)
   end function mg_get_e_from_p_T

   real(WP) function mg_get_h_from_p_T(this,p,T,y) result(h)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      rho=this%get_rho_from_p_T(p=p,T=T,y=y)
      h=this%get_e_from_p_rho(p=p,rho=rho,y=y)+p/rho
   end function mg_get_h_from_p_T

   subroutine mg_get_hk_from_p_T(this,p,T,y,hk)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(:), intent(out) :: hk
      hk(1)=this%get_h_from_p_T(p=p,T=T,y=y)
   end subroutine mg_get_hk_from_p_T

   !> Closed-form entropy of the consistent completion: s = cv*(ln(T) + gamma0*rho0/rho) + qp
   real(WP) function mg_get_s_from_p_T(this,p,T,y) result(s)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      rho=this%get_rho_from_p_T(p=p,T=T,y=y)
      s=this%cv*(log(T)+this%gr0m/rho)+this%qp
   end function mg_get_s_from_p_T

   real(WP) function mg_get_g_from_p_T(this,p,T,y) result(g)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      rho=this%get_rho_from_p_T(p=p,T=T,y=y)
      g=this%get_e_from_p_rho(p=p,rho=rho,y=y)+p/rho-T*(this%cv*(log(T)+this%gr0m/rho)+this%qp)
   end function mg_get_g_from_p_T

   real(WP) function mg_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      gruneisen=this%gr0m/rho
   end function mg_get_gruneisen_from_rho_e

   real(WP) function mg_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      rhoe=rho*this%get_e_from_p_rho(p=p,rho=rho,y=y)
   end function mg_get_rhoe_from_p_rho

   real(WP) function mg_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(mie_gruneisen), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      rho=this%get_rho_from_p_T(p=p,T=T,y=y)
      rhoe=rho*this%get_e_from_p_rho(p=p,rho=rho,y=y)
   end function mg_get_rhoe_from_p_T

   subroutine mg_print(this)
      use messager, only: log
      use string,   only: str_long
      class(mie_gruneisen), intent(in) :: this
      character(len=str_long) :: msg
      integer :: n
      write(msg,'(a,a)')         '[material:mie_gruneisen] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,es12.5)') 'rho0   = ',this%rho0;   call log(msg)
      write(msg,'(2x,a,es12.5)') 'c0     = ',this%c0;     call log(msg)
      write(msg,'(2x,a,es12.5)') 's1     = ',this%s1;     call log(msg)
      write(msg,'(2x,a,es12.5)') 's2     = ',this%s2;     call log(msg)
      write(msg,'(2x,a,es12.5)') 's3     = ',this%s3;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'gamma0 = ',this%gamma0; call log(msg)
      write(msg,'(2x,a,es12.5)') 'cv     = ',this%cv;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'T0     = ',this%T0;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'q      = ',this%q;      call log(msg)
      write(msg,'(2x,a,es12.5)') 'qp     = ',this%qp;     call log(msg)
      write(msg,'(2x,a,es12.5,a,es12.5,a)') 'etaknee= ',this%etaknee,' (rho = ',this%rho0/(1.0_WP-this%etaknee),')'; call log(msg)
      write(msg,'(2x,a,es12.5,a,es12.5,a,i0,a)') 'Tref table on [',this%etamin,',',this%etamax,'] (',this%ntab,' points)'; call log(msg)
      write(msg,'(2x,a,i0)')     'ns     = ',this%ns;     call log(msg)
      do n=1,this%ns
         write(msg,'(2x,a,i0,a,a)') 'species(',n,') = ',trim(this%species_names(n)); call log(msg)
      end do
   end subroutine mg_print

   subroutine mg_finalize(this)
      class(mie_gruneisen), intent(inout) :: this
      if (allocated(this%species_names)) deallocate(this%species_names)
      if (allocated(this%Ttab)) deallocate(this%Ttab)
   end subroutine mg_finalize

end module mie_gruneisen_class
