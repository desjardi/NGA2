!> Abstract material/thermodynamic model.
!> Concrete materials (stiffened_gas, ideal_gas, ig_mixture, sg_mixture, ...) extend this.
!> All thermo methods take a composition argument y(:); pure substances ignore it.
module material_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private

   public :: material

   type, abstract :: material
      character(len=str_medium) :: name = 'UNNAMED_MATERIAL'
      integer :: ns = 1   !< Number of species (1 = pure substance)
      character(len=str_medium), dimension(:), allocatable :: species_names   !< Per-species names (length ns)
   contains
      procedure(get_p_from_rho_e_iface),    deferred :: get_p_from_rho_e
      procedure(get_T_from_p_rho_iface),    deferred :: get_T_from_p_rho
      procedure(get_c_from_p_rho_iface),    deferred :: get_c_from_p_rho
      procedure(get_e_from_p_rho_iface),    deferred :: get_e_from_p_rho
      procedure(get_e_from_p_T_iface),      deferred :: get_e_from_p_T
      procedure(get_p_from_rho_T_iface),    deferred :: get_p_from_rho_T
      procedure(get_rho_from_p_T_iface),    deferred :: get_rho_from_p_T
      procedure(get_cv_from_rho_T_iface),   deferred :: get_cv_from_rho_T
      procedure(get_h_from_p_T_iface),      deferred :: get_h_from_p_T
      procedure(get_hk_from_p_T_iface),     deferred :: get_hk_from_p_T
      procedure(get_s_from_p_T_iface),      deferred :: get_s_from_p_T
      procedure(get_g_from_p_T_iface),      deferred :: get_g_from_p_T
      procedure(get_gruneisen_iface),       deferred :: get_gruneisen_from_rho_e
      procedure(get_rhoe_from_p_rho_iface), deferred :: get_rhoe_from_p_rho
      procedure(get_rhoe_from_p_T_iface),   deferred :: get_rhoe_from_p_T
      procedure(print_iface),               deferred :: print
      procedure(finalize_iface),            deferred :: finalize
      !> Primitives from (rho,e), the well-conditioned state the solver actually holds -- prefer these
      !> over get_T/get_c_from_p_rho in any (rho,e)-driven hot loop. Deferred so every material gives an
      !> OPTIMAL impl (ideal_gas/sg/nasg/igmix closed forms; CoolProp a single (rho,e) flash).
      procedure(get_T_from_rho_e_iface),  deferred :: get_T_from_rho_e
      procedure(get_c_from_rho_e_iface),  deferred :: get_c_from_rho_e
      procedure(get_cv_from_rho_e_iface), deferred :: get_cv_from_rho_e
   end type material

   abstract interface
      real(WP) function get_p_from_rho_e_iface(this,rho,e,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,e
         real(WP), dimension(:), intent(in) :: y
      end function get_p_from_rho_e_iface

      real(WP) function get_T_from_p_rho_iface(this,p,rho,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,rho
         real(WP), dimension(:), intent(in) :: y
      end function get_T_from_p_rho_iface

      real(WP) function get_c_from_p_rho_iface(this,p,rho,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,rho
         real(WP), dimension(:), intent(in) :: y
      end function get_c_from_p_rho_iface

      real(WP) function get_e_from_p_rho_iface(this,p,rho,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,rho
         real(WP), dimension(:), intent(in) :: y
      end function get_e_from_p_rho_iface

      real(WP) function get_e_from_p_T_iface(this,p,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function get_e_from_p_T_iface

      real(WP) function get_p_from_rho_T_iface(this,rho,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,T
         real(WP), dimension(:), intent(in) :: y
      end function get_p_from_rho_T_iface

      real(WP) function get_rho_from_p_T_iface(this,p,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function get_rho_from_p_T_iface

      !> cv = (de/dT)_rho; constant for calorically-perfect EOS, composition-weighted for mixtures
      real(WP) function get_cv_from_rho_T_iface(this,rho,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,T
         real(WP), dimension(:), intent(in) :: y
      end function get_cv_from_rho_T_iface

      real(WP) function get_h_from_p_T_iface(this,p,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function get_h_from_p_T_iface

      real(WP) function get_s_from_p_T_iface(this,p,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function get_s_from_p_T_iface

      !> Partial specific enthalpies h_k = dH/dm_k at fixed (p,T,m_j) for k=1..ns
      subroutine get_hk_from_p_T_iface(this,p,T,y,hk)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
         real(WP), dimension(:), intent(out) :: hk
      end subroutine get_hk_from_p_T_iface

      real(WP) function get_g_from_p_T_iface(this,p,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function get_g_from_p_T_iface

      real(WP) function get_gruneisen_iface(this,rho,e,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,e
         real(WP), dimension(:), intent(in) :: y
      end function get_gruneisen_iface

      real(WP) function get_T_from_rho_e_iface(this,rho,e,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,e
         real(WP), dimension(:), intent(in) :: y
      end function get_T_from_rho_e_iface

      real(WP) function get_c_from_rho_e_iface(this,rho,e,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,e
         real(WP), dimension(:), intent(in) :: y
      end function get_c_from_rho_e_iface

      real(WP) function get_cv_from_rho_e_iface(this,rho,e,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: rho,e
         real(WP), dimension(:), intent(in) :: y
      end function get_cv_from_rho_e_iface

      real(WP) function get_rhoe_from_p_rho_iface(this,p,rho,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,rho
         real(WP), dimension(:), intent(in) :: y
      end function get_rhoe_from_p_rho_iface

      real(WP) function get_rhoe_from_p_T_iface(this,p,T,y)
         import :: WP,material
         class(material), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function get_rhoe_from_p_T_iface

      subroutine print_iface(this)
         import :: material
         class(material), intent(in) :: this
      end subroutine print_iface

      subroutine finalize_iface(this)
         import :: material
         class(material), intent(inout) :: this
      end subroutine finalize_iface
   end interface

end module material_class
