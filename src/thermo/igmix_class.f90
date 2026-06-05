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
      procedure :: get_h_from_p_T          => igmix_get_h_from_p_T
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

   !> Mixture parameters: mass-fraction-weighted averages over species
   pure subroutine mix_coeffs(this,y,cv_mix,cp_mix,R_mix,gamma_mix,q_mix,qp_mix)
      class(igmix),                     intent(in)  :: this
      real(WP), dimension(:),           intent(in)  :: y
      real(WP), optional,               intent(out) :: cv_mix,cp_mix,R_mix,gamma_mix,q_mix,qp_mix
      real(WP) :: c,p,qx,qpx
      c =sum(y(1:this%ns)*this%cv(1:this%ns))
      p =sum(y(1:this%ns)*this%cp(1:this%ns))
      qx =sum(y(1:this%ns)*this%q (1:this%ns))
      qpx=sum(y(1:this%ns)*this%qp(1:this%ns))
      if (present(cv_mix))    cv_mix   =c
      if (present(cp_mix))    cp_mix   =p
      if (present(R_mix))     R_mix    =p-c
      if (present(gamma_mix)) gamma_mix=p/c
      if (present(q_mix))     q_mix    =qx
      if (present(qp_mix))    qp_mix   =qpx
   end subroutine mix_coeffs

   real(WP) function igmix_get_p_from_rho_e(this,rho,e,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: gam,qm
      call mix_coeffs(this,y,gamma_mix=gam,q_mix=qm)
      p=(gam-1.0_WP)*rho*(e-qm)
   end function igmix_get_p_from_rho_e

   real(WP) function igmix_get_T_from_p_rho(this,p,rho,y) result(T)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm
      call mix_coeffs(this,y,R_mix=Rm)
      T=p/(rho*Rm)
   end function igmix_get_T_from_p_rho

   real(WP) function igmix_get_c_from_p_rho(this,p,rho,y) result(c)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: gam
      call mix_coeffs(this,y,gamma_mix=gam)
      c=sqrt(max(0.0_WP,gam*p/rho))
   end function igmix_get_c_from_p_rho

   real(WP) function igmix_get_e_from_p_rho(this,p,rho,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: gam,qm
      call mix_coeffs(this,y,gamma_mix=gam,q_mix=qm)
      e=p/((gam-1.0_WP)*rho)+qm
   end function igmix_get_e_from_p_rho

   real(WP) function igmix_get_e_from_p_T(this,p,T,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cvm,qm
      call mix_coeffs(this,y,cv_mix=cvm,q_mix=qm)
      e=cvm*T+qm
   end function igmix_get_e_from_p_T

   real(WP) function igmix_get_p_from_rho_T(this,rho,T,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm
      call mix_coeffs(this,y,R_mix=Rm)
      p=rho*Rm*T
   end function igmix_get_p_from_rho_T

   real(WP) function igmix_get_rho_from_p_T(this,p,T,y) result(rho)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: Rm
      call mix_coeffs(this,y,R_mix=Rm)
      rho=p/(Rm*T)
   end function igmix_get_rho_from_p_T

   real(WP) function igmix_get_h_from_p_T(this,p,T,y) result(h)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cpm,qm
      call mix_coeffs(this,y,cp_mix=cpm,q_mix=qm)
      h=cpm*T+qm
   end function igmix_get_h_from_p_T

   !> Mixture entropy: sum_n y_n * s_n(x_n * p, T) using mole-fraction partial pressures
   real(WP) function igmix_get_s_from_p_T(this,p,T,y) result(s)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: ysum,xRsum,xn,pn,Rn
      integer :: n
      ysum =sum(y(1:this%ns))
      ! Mole fractions weighted by R: x_n = (y_n*R_n) / sum_k (y_k*R_k)
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
      real(WP) :: gam
      call mix_coeffs(this,y,gamma_mix=gam)
      gruneisen=gam-1.0_WP
   end function igmix_get_gruneisen_from_rho_e

   real(WP) function igmix_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: gam,qm
      call mix_coeffs(this,y,gamma_mix=gam,q_mix=qm)
      rhoe=p/(gam-1.0_WP)+rho*qm
   end function igmix_get_rhoe_from_p_rho

   real(WP) function igmix_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cvm,Rm,qm
      call mix_coeffs(this,y,cv_mix=cvm,R_mix=Rm,q_mix=qm)
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
