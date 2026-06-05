!> Noble-Abel stiffened-gas EOS (pure substance, ns=1)
!> Extends stiffened_gas; adds co-volume parameter b (reduces to stiffened gas when b=0)
module nasg_class
   use precision,           only: WP
   use stiffened_gas_class, only: stiffened_gas
   implicit none
   private

   public :: nasg

   type, extends(stiffened_gas) :: nasg
      real(WP) :: b = 0.0_WP   !< Co-volume
   contains
      procedure, private :: nasg_initialize
      generic   :: initialize              => nasg_initialize
      procedure :: get_p_from_rho_e        => nasg_get_p_from_rho_e
      procedure :: get_T_from_p_rho        => nasg_get_T_from_p_rho
      procedure :: get_c_from_p_rho        => nasg_get_c_from_p_rho
      procedure :: get_e_from_p_rho        => nasg_get_e_from_p_rho
      procedure :: get_p_from_rho_T        => nasg_get_p_from_rho_T
      procedure :: get_rho_from_p_T        => nasg_get_rho_from_p_T
      procedure :: get_h_from_p_T          => nasg_get_h_from_p_T
      procedure :: get_g_from_p_T          => nasg_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=> nasg_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho     => nasg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => nasg_get_rhoe_from_p_T
      procedure :: print                   => nasg_print
      ! Inherited unchanged from stiffened_gas:
      !   get_e_from_p_T   (b does not enter the calorically-perfect internal energy at given p,T)
      !   get_s_from_p_T   (same)
      !   finalize         (no NASG-specific allocatables)
   end type nasg

contains

   subroutine nasg_initialize(this,gamma,pinf,b,cv,q,qp,name)
      class(nasg), intent(inout) :: this
      real(WP), intent(in) :: gamma,pinf,b,cv,q,qp
      character(len=*), intent(in), optional :: name
      ! Initialize stiffened-gas parent (sets gamma, pinf, cv, cp, R, q, qp, ns=1, name, species_names)
      if (present(name)) then
         call this%stiffened_gas%initialize(gamma=gamma,pinf=pinf,cv=cv,q=q,qp=qp,name=name)
      else
         call this%stiffened_gas%initialize(gamma=gamma,pinf=pinf,cv=cv,q=q,qp=qp)
      end if
      this%b=b
   end subroutine nasg_initialize

   real(WP) function nasg_get_p_from_rho_e(this,rho,e,y) result(p)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      p=(this%gamma-1.0_WP)*rho*(e-this%q)/(1.0_WP-this%b*rho)-this%gamma*this%pinf
   end function nasg_get_p_from_rho_e

   real(WP) function nasg_get_T_from_p_rho(this,p,rho,y) result(T)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      T=(p+this%pinf)*(1.0_WP-this%b*rho)/(this%R*rho)
   end function nasg_get_T_from_p_rho

   real(WP) function nasg_get_c_from_p_rho(this,p,rho,y) result(c)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      c=sqrt(max(0.0_WP,this%gamma*(p+this%pinf)/(rho*(1.0_WP-this%b*rho))))
   end function nasg_get_c_from_p_rho

   real(WP) function nasg_get_e_from_p_rho(this,p,rho,y) result(e)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      e=(1.0_WP-this%b*rho)*(p+this%gamma*this%pinf)/((this%gamma-1.0_WP)*rho)+this%q
   end function nasg_get_e_from_p_rho

   real(WP) function nasg_get_p_from_rho_T(this,rho,T,y) result(p)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      p=this%R*rho*T/(1.0_WP-this%b*rho)-this%pinf
   end function nasg_get_p_from_rho_T

   real(WP) function nasg_get_rho_from_p_T(this,p,T,y) result(rho)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rho=(p+this%pinf)/(this%R*T+this%b*(p+this%pinf))
   end function nasg_get_rho_from_p_T

   real(WP) function nasg_get_h_from_p_T(this,p,T,y) result(h)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      h=this%cp*T+this%b*p+this%q
   end function nasg_get_h_from_p_T

   real(WP) function nasg_get_g_from_p_T(this,p,T,y) result(g)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g=this%cp*T+this%b*p+this%q-T*(this%cp*log(T)-this%R*log(p+this%pinf)+this%qp)
   end function nasg_get_g_from_p_T

   real(WP) function nasg_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      gruneisen=(this%gamma-1.0_WP)/(1.0_WP-this%b*rho)
   end function nasg_get_gruneisen_from_rho_e

   real(WP) function nasg_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      rhoe=(1.0_WP-this%b*rho)*(p+this%gamma*this%pinf)/(this%gamma-1.0_WP)+rho*this%q
   end function nasg_get_rhoe_from_p_rho

   real(WP) function nasg_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rhoe=((p+this%gamma*this%pinf)*this%cv*T+this%q*(p+this%pinf))/(this%R*T+this%b*(p+this%pinf))
   end function nasg_get_rhoe_from_p_T

   subroutine nasg_print(this)
      use messager, only: log
      use string,   only: str_long
      class(nasg), intent(in) :: this
      character(len=str_long) :: msg
      integer :: n
      write(msg,'(a,a)')         '[material:nasg] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,es12.5)') 'gamma = ',this%gamma; call log(msg)
      write(msg,'(2x,a,es12.5)') 'pinf  = ',this%pinf;  call log(msg)
      write(msg,'(2x,a,es12.5)') 'b     = ',this%b;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'cv    = ',this%cv;    call log(msg)
      write(msg,'(2x,a,es12.5)') 'cp    = ',this%cp;    call log(msg)
      write(msg,'(2x,a,es12.5)') 'R     = ',this%R;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'q     = ',this%q;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'qp    = ',this%qp;    call log(msg)
      write(msg,'(2x,a,i0)')     'ns    = ',this%ns;    call log(msg)
      do n=1,this%ns
         write(msg,'(2x,a,i0,a,a)') 'species(',n,') = ',trim(this%species_names(n)); call log(msg)
      end do
   end subroutine nasg_print

end module nasg_class
