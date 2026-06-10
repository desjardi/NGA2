!> Ideal-gas mixture: mass-fraction-weighted blend of ns calorically-perfect species
!> Species parameters: gamma, cv, q, qp + names
module igmix_class
   use precision,      only: WP
   use string,         only: str_medium
   use messager,       only: die
   use material_class, only: material
   implicit none
   private

   public :: igmix

   type, extends(material) :: igmix
      real(WP), dimension(:), allocatable :: gamma,cv,cp,R,q,qp
   contains
      procedure, private :: igmix_initialize
      generic   :: initialize              => igmix_initialize
      procedure :: get_p_from_rho_e        => igmix_get_p_from_rho_e
      procedure :: get_T_from_p_rho        => igmix_get_T_from_p_rho
      procedure :: get_c_from_p_rho        => igmix_get_c_from_p_rho
      procedure :: get_e_from_p_rho        => igmix_get_e_from_p_rho
      procedure :: get_e_from_p_T          => igmix_get_e_from_p_T
      procedure :: get_p_from_rho_T        => igmix_get_p_from_rho_T
      procedure :: get_rho_from_p_T        => igmix_get_rho_from_p_T
      procedure :: get_cv_from_rho_T       => igmix_get_cv_from_rho_T
      procedure :: get_h_from_p_T          => igmix_get_h_from_p_T
      procedure :: get_hk_from_p_T         => igmix_get_hk_from_p_T
      procedure :: get_s_from_p_T          => igmix_get_s_from_p_T
      procedure :: get_g_from_p_T          => igmix_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=> igmix_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho     => igmix_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => igmix_get_rhoe_from_p_T
      procedure :: print                   => igmix_print
      procedure :: finalize                => igmix_finalize
   end type igmix

contains

   !> Initialize from per-species arrays (length ns, inferred from size(gamma))
   subroutine igmix_initialize(this,gamma,cv,q,qp,species_names,name)
      class(igmix), intent(inout) :: this
      real(WP), dimension(:),                  intent(in) :: gamma,cv,q,qp
      character(len=*), dimension(:),          intent(in) :: species_names
      character(len=*),              optional, intent(in) :: name
      integer :: ns
      ns=size(gamma)
      if (size(cv).ne.ns.or.size(q).ne.ns.or.size(qp).ne.ns) call die('[igmix initialize] gamma/cv/q/qp array size mismatch')
      if (size(species_names).ne.ns) call die('[igmix initialize] species_names size mismatch')
      this%ns=ns
      if (present(name)) this%name=name
      if (allocated(this%gamma)) deallocate(this%gamma,this%cv,this%cp,this%R,this%q,this%qp)
      if (allocated(this%species_names)) deallocate(this%species_names)
      allocate(this%gamma(ns),this%cv(ns),this%cp(ns),this%R(ns),this%q(ns),this%qp(ns))
      allocate(this%species_names(ns))
      this%gamma=gamma
      this%cv   =cv
      this%cp   =gamma*cv
      this%R    =(gamma-1.0_WP)*cv
      this%q    =q
      this%qp   =qp
      this%species_names=species_names
   end subroutine igmix_initialize

   real(WP) function igmix_get_p_from_rho_e(this,rho,e,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cvm,qm
      Rm =sum(y(1:this%ns)*this%R (1:this%ns))
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      qm =sum(y(1:this%ns)*this%q (1:this%ns))
      p=Rm*rho*(e-qm)/cvm
   end function igmix_get_p_from_rho_e

   real(WP) function igmix_get_T_from_p_rho(this,p,rho,y) result(T)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm
      Rm=sum(y(1:this%ns)*this%R(1:this%ns))
      T=p/(rho*Rm)
   end function igmix_get_T_from_p_rho

   real(WP) function igmix_get_c_from_p_rho(this,p,rho,y) result(c)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cvm,cpm
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      cpm=sum(y(1:this%ns)*this%cp(1:this%ns))
      c=sqrt(max(0.0_WP,cpm/cvm*p/rho))
   end function igmix_get_c_from_p_rho

   real(WP) function igmix_get_e_from_p_rho(this,p,rho,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cvm,qm
      Rm =sum(y(1:this%ns)*this%R (1:this%ns))
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      qm =sum(y(1:this%ns)*this%q (1:this%ns))
      e=cvm*p/(Rm*rho)+qm
   end function igmix_get_e_from_p_rho

   real(WP) function igmix_get_e_from_p_T(this,p,T,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cvm,qm
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      qm =sum(y(1:this%ns)*this%q (1:this%ns))
      e=cvm*T+qm
   end function igmix_get_e_from_p_T

   real(WP) function igmix_get_p_from_rho_T(this,rho,T,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm
      Rm=sum(y(1:this%ns)*this%R(1:this%ns))
      p=rho*Rm*T
   end function igmix_get_p_from_rho_T

   real(WP) function igmix_get_rho_from_p_T(this,p,T,y) result(rho)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm
      Rm=sum(y(1:this%ns)*this%R(1:this%ns))
      rho=p/(Rm*T)
   end function igmix_get_rho_from_p_T

   real(WP) function igmix_get_cv_from_rho_T(this,rho,T,y) result(cv)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      cv=sum(y(1:this%ns)*this%cv(1:this%ns))
   end function igmix_get_cv_from_rho_T

   real(WP) function igmix_get_h_from_p_T(this,p,T,y) result(h)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cpm,qm
      cpm=sum(y(1:this%ns)*this%cp(1:this%ns))
      qm =sum(y(1:this%ns)*this%q (1:this%ns))
      h=cpm*T+qm
   end function igmix_get_h_from_p_T

   subroutine igmix_get_hk_from_p_T(this,p,T,y,hk)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(:), intent(out) :: hk
      hk(1:this%ns)=this%cp(1:this%ns)*T+this%q(1:this%ns)
   end subroutine igmix_get_hk_from_p_T

   !> Mixture entropy: sum_n y_n * s_n(x_n * p, T) using mole-fraction partial pressures
   real(WP) function igmix_get_s_from_p_T(this,p,T,y) result(s)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: xRsum,xn,pn,Rn
      integer :: n
      xRsum=0.0_WP
      do n=1,this%ns
         xRsum=xRsum+y(n)*this%R(n)
      end do
      if (xRsum.le.tiny(1.0_WP)) call die('[igmix s_from_p_T] non-positive species gas-constant weighted sum')
      s=0.0_WP
      do n=1,this%ns
         Rn=this%R(n)
         xn=y(n)*Rn/xRsum
         pn=max(xn*p,tiny(1.0_WP))
         s=s+y(n)*(this%cp(n)*log(T)-Rn*log(pn)+this%qp(n))
      end do
   end function igmix_get_s_from_p_T

   real(WP) function igmix_get_g_from_p_T(this,p,T,y) result(g)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g=this%get_h_from_p_T(p,T,y)-T*this%get_s_from_p_T(p,T,y)
   end function igmix_get_g_from_p_T

   real(WP) function igmix_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cvm
      Rm =sum(y(1:this%ns)*this%R (1:this%ns))
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      gruneisen=Rm/cvm
   end function igmix_get_gruneisen_from_rho_e

   real(WP) function igmix_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cvm,qm
      Rm =sum(y(1:this%ns)*this%R (1:this%ns))
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      qm =sum(y(1:this%ns)*this%q (1:this%ns))
      rhoe=cvm*p/Rm+rho*qm
   end function igmix_get_rhoe_from_p_rho

   real(WP) function igmix_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm,cvm,qm
      Rm =sum(y(1:this%ns)*this%R (1:this%ns))
      cvm=sum(y(1:this%ns)*this%cv(1:this%ns))
      qm =sum(y(1:this%ns)*this%q (1:this%ns))
      rhoe=p*(cvm*T+qm)/(Rm*T)
   end function igmix_get_rhoe_from_p_T

   subroutine igmix_print(this)
      use messager, only: log
      use string,   only: str_long
      class(igmix), intent(in) :: this
      character(len=str_long) :: msg
      integer :: n
      write(msg,'(a,a)') '[material:igmix] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,i0)') 'ns    = ',this%ns; call log(msg)
      do n=1,this%ns
         write(msg,'(2x,a,i0,a,a)') 'species(',n,') = ',trim(this%species_names(n)); call log(msg)
         write(msg,'(4x,a,es12.5,a,es12.5,a,es12.5,a,es12.5)') &
         &  'gamma=',this%gamma(n),' cv=',this%cv(n),' q=',this%q(n),' qp=',this%qp(n)
         call log(msg)
      end do
   end subroutine igmix_print

   subroutine igmix_finalize(this)
      class(igmix), intent(inout) :: this
      if (allocated(this%gamma))         deallocate(this%gamma)
      if (allocated(this%cv))            deallocate(this%cv)
      if (allocated(this%cp))            deallocate(this%cp)
      if (allocated(this%R))             deallocate(this%R)
      if (allocated(this%q))             deallocate(this%q)
      if (allocated(this%qp))            deallocate(this%qp)
      if (allocated(this%species_names)) deallocate(this%species_names)
      this%ns=0
   end subroutine igmix_finalize

end module igmix_class
