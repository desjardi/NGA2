!> Stiffened-gas EOS (pure substance, ns=1)
!> Extends ideal_gas; adds stiffening parameter pinf (reduces to ideal gas when pinf=0)
module stiffened_gas_class
   use precision,       only: WP
   use ideal_gas_class, only: ideal_gas
   implicit none
   private

   public :: stiffened_gas

   type, extends(ideal_gas) :: stiffened_gas
      real(WP) :: pinf = 0.0_WP
   contains
      procedure, private :: sg_initialize
      generic   :: initialize              => sg_initialize
      procedure :: get_p_from_rho_e        => sg_get_p_from_rho_e
      procedure :: get_T_from_p_rho        => sg_get_T_from_p_rho
      procedure :: get_c_from_p_rho        => sg_get_c_from_p_rho
      procedure :: get_e_from_p_rho        => sg_get_e_from_p_rho
      procedure :: get_e_from_p_T          => sg_get_e_from_p_T
      procedure :: get_p_from_rho_T        => sg_get_p_from_rho_T
      procedure :: get_rho_from_p_T        => sg_get_rho_from_p_T
      procedure :: get_s_from_p_T          => sg_get_s_from_p_T
      procedure :: get_g_from_p_T          => sg_get_g_from_p_T
      procedure :: get_rhoe_from_p_rho     => sg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => sg_get_rhoe_from_p_T
      procedure :: print                   => sg_print
      ! Inherited unchanged from ideal_gas:
      !   get_h_from_p_T          (h = cp*T + q;     stiffening cancels in e + p/rho)
      !   get_gruneisen_from_rho_e (= gamma - 1)
      !   finalize                (deallocates species_names; no SG-specific allocatables)
   end type stiffened_gas

contains

   subroutine sg_initialize(this,gamma,pinf,cv,q,qp,name)
      class(stiffened_gas), intent(inout) :: this
      real(WP), intent(in) :: gamma,pinf,cv,q,qp
      character(len=*), intent(in), optional :: name
      ! Initialize ideal-gas parent (sets gamma, cv, cp, R, q, qp, ns=1, name, species_names)
      if (present(name)) then
         call this%ideal_gas%initialize(gamma=gamma,cv=cv,q=q,qp=qp,name=name)
      else
         call this%ideal_gas%initialize(gamma=gamma,cv=cv,q=q,qp=qp)
      end if
      ! Stiffened-gas extension
      this%pinf=pinf
   end subroutine sg_initialize

   real(WP) function sg_get_p_from_rho_e(this,rho,e,y) result(p)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      p=(this%gamma-1.0_WP)*rho*(e-this%q)-this%gamma*this%pinf
   end function sg_get_p_from_rho_e

   real(WP) function sg_get_T_from_p_rho(this,p,rho,y) result(T)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      T=(p+this%pinf)/(this%R*rho)
   end function sg_get_T_from_p_rho

   real(WP) function sg_get_c_from_p_rho(this,p,rho,y) result(c)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      c=sqrt(max(0.0_WP,this%gamma*(p+this%pinf)/rho))
   end function sg_get_c_from_p_rho

   real(WP) function sg_get_e_from_p_rho(this,p,rho,y) result(e)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      e=(p+this%gamma*this%pinf)/((this%gamma-1.0_WP)*rho)+this%q
   end function sg_get_e_from_p_rho

   real(WP) function sg_get_e_from_p_T(this,p,T,y) result(e)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      e=this%cv*T*(p+this%gamma*this%pinf)/(p+this%pinf)+this%q
   end function sg_get_e_from_p_T

   real(WP) function sg_get_p_from_rho_T(this,rho,T,y) result(p)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      p=this%R*rho*T-this%pinf
   end function sg_get_p_from_rho_T

   real(WP) function sg_get_rho_from_p_T(this,p,T,y) result(rho)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rho=(p+this%pinf)/(this%R*T)
   end function sg_get_rho_from_p_T

   real(WP) function sg_get_s_from_p_T(this,p,T,y) result(s)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      s=this%cp*log(T)-this%R*log(p+this%pinf)+this%qp
   end function sg_get_s_from_p_T

   real(WP) function sg_get_g_from_p_T(this,p,T,y) result(g)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g=(this%cp-this%qp)*T-T*(this%cp*log(T)-this%R*log(p+this%pinf))+this%q
   end function sg_get_g_from_p_T

   real(WP) function sg_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      rhoe=(p+this%gamma*this%pinf)/(this%gamma-1.0_WP)+rho*this%q
   end function sg_get_rhoe_from_p_rho

   real(WP) function sg_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(stiffened_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rhoe=((p+this%gamma*this%pinf)*this%cv*T+this%q*(p+this%pinf))/(this%R*T)
   end function sg_get_rhoe_from_p_T

   subroutine sg_print(this)
      use messager, only: log
      use string,   only: str_long
      class(stiffened_gas), intent(in) :: this
      character(len=str_long) :: msg
      integer :: n
      write(msg,'(a,a)')         '[material:stiffened_gas] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,es12.5)') 'gamma = ',this%gamma; call log(msg)
      write(msg,'(2x,a,es12.5)') 'pinf  = ',this%pinf;  call log(msg)
      write(msg,'(2x,a,es12.5)') 'cv    = ',this%cv;    call log(msg)
      write(msg,'(2x,a,es12.5)') 'cp    = ',this%cp;    call log(msg)
      write(msg,'(2x,a,es12.5)') 'R     = ',this%R;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'q     = ',this%q;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'qp    = ',this%qp;    call log(msg)
      write(msg,'(2x,a,i0)')     'ns    = ',this%ns;    call log(msg)
      do n=1,this%ns
         write(msg,'(2x,a,i0,a,a)') 'species(',n,') = ',trim(this%species_names(n)); call log(msg)
      end do
   end subroutine sg_print

end module stiffened_gas_class
