!> Calorically perfect ideal-gas EOS (pure substance, ns=1)
!> Parameters: gamma, cv, q (energy of formation), qp (entropy of formation)
module ideal_gas_class
   use precision,      only: WP
   use material_class, only: material
   implicit none
   private

   public :: ideal_gas

   type, extends(material) :: ideal_gas
      real(WP) :: gamma = 0.0_WP
      real(WP) :: cv    = 0.0_WP
      real(WP) :: cp    = 0.0_WP   !< cp = gamma*cv, stored at init
      real(WP) :: R     = 0.0_WP   !< R  = (gamma-1)*cv, stored at init
      real(WP) :: q     = 0.0_WP
      real(WP) :: qp    = 0.0_WP
   contains
      procedure, private :: ig_initialize
      generic   :: initialize              => ig_initialize
      procedure :: get_p_from_rho_e        => ig_get_p_from_rho_e
      procedure :: get_T_from_p_rho        => ig_get_T_from_p_rho
      procedure :: get_c_from_p_rho        => ig_get_c_from_p_rho
      procedure :: get_e_from_p_rho        => ig_get_e_from_p_rho
      procedure :: get_e_from_p_T          => ig_get_e_from_p_T
      procedure :: get_p_from_rho_T        => ig_get_p_from_rho_T
      procedure :: get_rho_from_p_T        => ig_get_rho_from_p_T
      procedure :: get_h_from_p_T          => ig_get_h_from_p_T
      procedure :: get_s_from_p_T          => ig_get_s_from_p_T
      procedure :: get_g_from_p_T          => ig_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=> ig_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho     => ig_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       => ig_get_rhoe_from_p_T
      procedure :: print                   => ig_print
      procedure :: finalize                => ig_finalize
   end type ideal_gas

contains

   subroutine ig_initialize(this,gamma,cv,q,qp,name)
      class(ideal_gas), intent(inout) :: this
      real(WP), intent(in) :: gamma,cv,q,qp
      character(len=*), intent(in), optional :: name
      this%gamma = gamma
      this%cv    = cv
      this%cp    = gamma*cv
      this%R     = (gamma-1.0_WP)*cv
      this%q     = q
      this%qp    = qp
      this%ns    = 1
      if (present(name)) this%name=name
      allocate(this%species_names(1)); this%species_names(1)=this%name
   end subroutine ig_initialize

   real(WP) function ig_get_p_from_rho_e(this,rho,e,y) result(p)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      p = (this%gamma-1.0_WP)*rho*(e-this%q)
   end function ig_get_p_from_rho_e

   real(WP) function ig_get_T_from_p_rho(this,p,rho,y) result(T)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      T = p/(this%R*rho)
   end function ig_get_T_from_p_rho

   real(WP) function ig_get_c_from_p_rho(this,p,rho,y) result(c)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      c = sqrt(max(0.0_WP,this%gamma*p/rho))
   end function ig_get_c_from_p_rho

   real(WP) function ig_get_e_from_p_rho(this,p,rho,y) result(e)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      e = p/((this%gamma-1.0_WP)*rho)+this%q
   end function ig_get_e_from_p_rho

   real(WP) function ig_get_e_from_p_T(this,p,T,y) result(e)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      e = this%cv*T+this%q
   end function ig_get_e_from_p_T

   real(WP) function ig_get_p_from_rho_T(this,rho,T,y) result(p)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      p = this%R*rho*T
   end function ig_get_p_from_rho_T

   real(WP) function ig_get_rho_from_p_T(this,p,T,y) result(rho)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rho = p/(this%R*T)
   end function ig_get_rho_from_p_T

   real(WP) function ig_get_h_from_p_T(this,p,T,y) result(h)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      h = this%cp*T+this%q
   end function ig_get_h_from_p_T

   real(WP) function ig_get_s_from_p_T(this,p,T,y) result(s)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      s = this%cp*log(T)-this%R*log(p)+this%qp
   end function ig_get_s_from_p_T

   real(WP) function ig_get_g_from_p_T(this,p,T,y) result(g)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g = (this%cp-this%qp)*T-T*(this%cp*log(T)-this%R*log(p))+this%q
   end function ig_get_g_from_p_T

   real(WP) function ig_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      gruneisen = this%gamma-1.0_WP
   end function ig_get_gruneisen_from_rho_e

   real(WP) function ig_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      rhoe = p/(this%gamma-1.0_WP)+rho*this%q
   end function ig_get_rhoe_from_p_rho

   real(WP) function ig_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(ideal_gas), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      rhoe = p*(this%cv*T+this%q)/(this%R*T)
   end function ig_get_rhoe_from_p_T

   subroutine ig_finalize(this)
      class(ideal_gas), intent(inout) :: this
      if (allocated(this%species_names)) deallocate(this%species_names)
   end subroutine ig_finalize

   subroutine ig_print(this)
      use messager, only: log
      use string,   only: str_long
      class(ideal_gas), intent(in) :: this
      character(len=str_long) :: msg
      integer :: n
      write(msg,'(a,a)')         '[material:ideal_gas] ',trim(this%name); call log(msg)
      write(msg,'(2x,a,es12.5)') 'gamma = ',this%gamma; call log(msg)
      write(msg,'(2x,a,es12.5)') 'cv    = ',this%cv;    call log(msg)
      write(msg,'(2x,a,es12.5)') 'cp    = ',this%cp;    call log(msg)
      write(msg,'(2x,a,es12.5)') 'R     = ',this%R;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'q     = ',this%q;     call log(msg)
      write(msg,'(2x,a,es12.5)') 'qp    = ',this%qp;    call log(msg)
      write(msg,'(2x,a,i0)')     'ns    = ',this%ns;    call log(msg)
      do n=1,this%ns
         write(msg,'(2x,a,i0,a,a)') 'species(',n,') = ',trim(this%species_names(n)); call log(msg)
      end do
   end subroutine ig_print

end module ideal_gas_class
