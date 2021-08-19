!> Two-phase material modeling class:
!> Provides equations of state and temperature-dependent properties
module matm_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: matm
   
   !> Material modeling type intended for two-phase, liquid-gas flows
   type :: matm
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the structure
      character(len=str_medium) :: name='UNNAMED_MATM'    !< class name (default=UNNAMED_MATM)

      ! Flag for if things vary with temperature or not
      logical  :: const_prop
      ! Flag for if temperature equil of multiphase cells is enforced
      logical  :: mult_iso

      ! Liquid properties
      real(WP) :: gamm_l,Pref_l,q_l,b_l
      ! Gas properties
      real(WP) :: gamm_g,Pref_g,q_g,b_g

      ! Specific heat coefficients
      real(WP), dimension(:,:,:), allocatable :: cv_l, cv_g

      ! Pointers for local thermodynamic or flow variables
      real(WP), dimension(:,:,:), pointer     :: Grho,GU,GV,GW,GrhoE,GP
      real(WP), dimension(:,:,:), pointer     :: Lrho,LU,LV,LW,LrhoE,LP

   contains
      procedure :: EOS_gas, EOS_liquid                    !< Output solver to the screen
      procedure :: EOS_energy                             !< Calculates phase total energy from pressure and kinetic energy
      procedure :: EOS_temp                               !< Calculates phase temperature from pressure directly
      procedure :: EOS_all                                !< Calculates vol-avg pressure for entire domain from conserved variables

      procedure :: register_idealgas                      !< EOS available for gas-like fluids
      procedure :: register_stiffenedgas                  !< EOS for liquid-like fluids
      procedure :: register_NobleAbelstiffenedgas         !< EOS for liquid-like fluids

      procedure :: register_thermoflow_variables          !< Creates pointers from material models to flow solver variables

      procedure :: viscosity_water,  viscosity_air        !< Empirical models for temperature dependence of dynamic viscosity
      procedure :: spec_heat_water,  spec_heat_air        !< Empirical models for temperature dependence of specific heat (cv)
      procedure :: therm_cond_water, therm_cond_air       !< Empirical models for temperature dependence of thermal conductivity
      
   end type matm
   
   
   !> Declare two-phase compressible solver constructor
   interface matm
      procedure constructor
   end interface matm
   
contains
   
   
   !> Default constructor for two-phase material modeling
   function constructor(cfg,name) result(self)
     implicit none
     type(matm) :: self
     class(config), target, intent(in) :: cfg
     character(len=*), optional :: name

     ! Set the name for the solver                                                                                   
     if (present(name)) self%name=trim(adjustl(name))
     
     ! Point to pgrid object                                                                                         
     self%cfg=>cfg
     
     ! Zero EOS parameters
     self%gamm_l = 0.0_WP; self%Pref_l = 0.0_WP; self%q_l = 0.0_WP; self%b_l = 0.0_WP
     self%gamm_g = 0.0_WP; self%Pref_g = 0.0_WP; self%q_g = 0.0_WP; self%b_g = 0.0_WP

     ! Allocate specific heat
     allocate(self%cv_l(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%cv_l=0.0_WP
     allocate(self%cv_g(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%cv_g=0.0_WP

   end function constructor

   subroutine register_idealgas(this,phase,gamm)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     real(WP), intent(in) :: gamm

     select case(trim(adjustl(phase)))
     case('liquid')
        this%gamm_l = gamm
     case('gas')
        this%gamm_g = gamm
     end select
     
   end subroutine register_idealgas

   subroutine register_stiffenedgas(this,phase,gamm,Pref)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     real(WP), intent(in) :: gamm,Pref

     select case(trim(adjustl(phase)))
     case('liquid')
        this%gamm_l = gamm
        this%Pref_l = Pref
     case('gas')
        this%gamm_g = gamm
        this%Pref_g = Pref
     end select
     
   end subroutine register_stiffenedgas

   subroutine register_NobleAbelstiffenedgas(this,phase,gamm,Pref,q,b)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     real(WP), intent(in) :: gamm,Pref,q,b

     select case(trim(adjustl(phase)))
     case('liquid')
        this%gamm_l = gamm
        this%Pref_l = Pref
        this%q_l    = q
        this%b_l    = b
     case('gas')
        this%gamm_g = gamm
        this%Pref_g = Pref
        this%q_g    = q
        this%b_g    = b
     end select
     
   end subroutine register_NobleAbelstiffenedgas

   
   subroutine register_thermoflow_variables(this,phase,rho,u,v,w,rhoe,p)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: rho,u,v,w,rhoE,p
     
     select case(trim(adjustl(phase)))
     case('liquid')
        this%Grho  => rho
        this%GU    => u
        this%GV    => v
        this%GW    => w
        this%GrhoE => rhoE
        this%GP    => p
     case('gas')
        this%Lrho  => rho
        this%LU    => u
        this%LV    => v
        this%LW    => w
        this%LrhoE => rhoE
        this%LP    => p
     end select

   end subroutine register_thermoflow_variables
     

   function EOS_liquid(this,i,j,k,flag) result(property)
     implicit none
     class(matm), intent(inout) :: this
     integer,intent(in) :: i,j,k
     character(len=1),intent(in) :: flag
     real(WP) :: property

     select case(flag)
     case('p') ! Pressure
        property = (this%gamm_l-1.0_WP)/(1.0_WP-this%Lrho(i,j,k)*this%b_l)*&
             (this%LrhoE(i,j,k)-0.5_WP*this%Lrho(i,j,k)*&
             (this%LU(i,j,k)**2+this%LV(i,j,k)**2+this%LW(i,j,k)**2) &
             -this%Lrho(i,j,k)*this%q_l)-this%gamm_l*this%Pref_l
     case('K') ! Bulk modulus
        property = this%LP(i,j,k)+this%Pref_l
        property = this%gamm_l/(1.0_WP-this%Lrho(i,j,k)*this%b_l)*property
     case('T') ! Temperature
        property = (-(1.0_WP-this%Lrho(i,j,k)*this%b_l)*this%Pref_l+this%LrhoE(i,j,k) &
             -0.5_WP*this%Lrho(i,j,k)*(this%LU(i,j,k)**2+this%LV(i,j,k)**2+this%LW(i,j,k)**2) &
             -this%Lrho(i,j,k)*this%q_l) / (this%Lrho(i,j,k)*this%cv_l(i,j,k))
     case('v') ! (1/rho)(dp/de), for viscous pressure term
        property = (this%gamm_l-1.0_WP)/(1.0_WP-this%Lrho(i,j,k)*this%b_l)
     case('M') ! Bulk modulus calculated from conserved variables
        property = (this%gamm_l-1.0_WP)/(1.0_WP-this%Lrho(i,j,k)*this%b_l)*(this%LrhoE(i,j,k) &
             -0.5_WP*this%Lrho(i,j,k)*(this%LU(i,j,k)**2+this%LV(i,j,k)**2+this%LW(i,j,k)**2) &
             -this%Lrho(i,j,k)*this%q_l)-this%gamm_l*this%Pref_l
        property = this%gamm_l/(1.0_WP-this%Lrho(i,j,k)*this%b_l)*(property+this%Pref_l)
     end select

     return
   end function EOS_liquid

   function EOS_gas(this,i,j,k,flag) result(property)
     implicit none
     class(matm), intent(inout) :: this
     integer,intent(in) :: i,j,k
     character(len=1),intent(in) :: flag
     real(WP) :: property

     select case(flag)
     case('p') ! Pressure
        property = (this%gamm_g-1.0_WP)/(1.0_WP-this%Grho(i,j,k)*this%b_g)*&
             (this%GrhoE(i,j,k)-0.5_WP*this%Grho(i,j,k)*&
             (this%GU(i,j,k)**2+this%GV(i,j,k)**2+this%GW(i,j,k)**2) &
             -this%Grho(i,j,k)*this%q_g)-this%gamm_g*this%Pref_g
     case('K') ! Bulk modulus
        property = this%GP(i,j,k)+this%Pref_g
        property = this%gamm_g/(1.0_WP-this%Grho(i,j,k)*this%b_g)*property
     case('T') ! Temperature
        property = (-(1.0_WP-this%Grho(i,j,k)*this%b_g)*this%Pref_g+this%GrhoE(i,j,k) &
             -0.5_WP*this%Grho(i,j,k)*(this%GU(i,j,k)**2+this%GV(i,j,k)**2+this%GW(i,j,k)**2) &
             -this%Grho(i,j,k)*this%q_g) / (this%Grho(i,j,k)*this%cv_g(i,j,k))
     case('v') ! (1/rho)(dp/de), for viscous pressure term
        property = (this%gamm_g-1.0_WP)/(1.0_WP-this%Grho(i,j,k)*this%b_g)
     case('M') ! Bulk modulus calculated from conserved variables
        property = (this%gamm_g-1.0_WP)/(1.0_WP-this%Grho(i,j,k)*this%b_g)*(this%GrhoE(i,j,k) &
             -0.5_WP*this%Grho(i,j,k)*(this%GU(i,j,k)**2+this%GV(i,j,k)**2+this%GW(i,j,k)**2) &
             -this%Grho(i,j,k)*this%q_g)-this%gamm_g*this%Pref_g
        property = this%gamm_g/(1.0_WP-this%Grho(i,j,k)*this%b_g)*(property+this%Pref_g)
     end select

     return
   end function EOS_gas

   function EOS_energy(this,pres,dens,uvel,vvel,wvel,phase) result(energy)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: pres,dens,uvel,vvel,wvel
     character(len=*),intent(in) :: phase
     real(WP) :: KE,energy

     KE = 0.5_WP*dens*(uvel**2+vvel**2+wvel**2)

     select case(trim(adjustl(phase)))
     case('liquid')
        energy = (pres + this%gamm_l*this%Pref_l)*(1.0_WP-dens*this%b_l)/(this%gamm_l-1.0_WP) + dens*this%q_l + KE
     case('gas')
        energy = (pres + this%gamm_g*this%Pref_g)*(1.0_WP-dens*this%b_g)/(this%gamm_g-1.0_WP) + dens*this%q_g + KE
     end select

     return
   end function EOS_energy

   function EOS_temp(this,pres,dens,cv,phase) result(temp)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: pres,dens,cv
     character(len=*),intent(in) :: phase
     real(WP) :: temp

     select case(trim(adjustl(phase)))
     case('liquid')
        temp = (1.0_WP-dens*this%b_l)*(pres+this%Pref_l)/(this%gamm_l-1.0_WP)/(dens*cv)
     case('gas')
        temp = (1.0_WP-dens*this%b_g)*(pres+this%Pref_g)/(this%gamm_g-1.0_WP)/(dens*cv)
     end select

     return
   end function EOS_temp

   function EOS_all(this,vf,buf_GrhoE,buf_LrhoE) result(buf_P)
     use vfs_class, only : vfs
     implicit none
     class(matm), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: buf_GrhoE,buf_LrhoE,buf_P

     buf_P = (1.0_WP-vf%VF)*((this%gamm_g-1.0_WP)/(1.0_WP-this%Grho*this%b_g)*(buf_GrhoE-0.5_WP* &
             this%Grho*(this%GU**2+this%GV**2+this%GW**2)-this%Grho*this%q_g)-this%gamm_g*this%Pref_g) + &
             (       vf%VF)*((this%gamm_l-1.0_WP)/(1.0_WP-this%Lrho*this%b_l)*(buf_LrhoE-0.5_WP* &
             this%Lrho*(this%LU**2+this%LV**2+this%LW**2)-this%Lrho*this%q_l)-this%gamm_l*this%Pref_l)

     return
   end function EOS_all

   function viscosity_air(this,T) result(mu)
     implicit none
     class(matm), intent(inout) :: this
     real(WP),intent(in) :: T
     real(WP) :: TC,mu

     TC = max(100.0_WP,min(T,800.0_WP))
     mu = 1.716e-5_WP*(TC/273.15_WP)**1.5_WP*(273.15+110.4_WP)/(TC+110.4_WP)
     ! mu = 2.5e-2_WP*(T/(1.488095e2/(1.4-1.0_WP)/717.6_WP))**0.75_WP

     return
   end function viscosity_air

   function viscosity_water(this,T) result(mu)
     implicit none
     class(matm), intent(inout) :: this
     real(WP),intent(in) :: T
     real(WP) :: TC,mu

     TC = max(276.15_WP,max(T,473.15_WP))
     mu = 2.414e-5_WP*10**(247.8_WP/(TC-140.0_WP))

     return
   end function viscosity_water

   function spec_heat_air(this,T) result(cv)
     implicit none
     class(matm), intent(inout) :: this
     real(WP),intent(in) :: T
     real(WP) :: TC,cv

     TC = max(173.15_WP,min(T,773.15_WP))
     cv = 717.8_WP+0.07075_WP*(TC-300_WP)+2.6125e-4_WP*(TC-300_WP)**2

     return
   end function spec_heat_air

   function spec_heat_water(this,T) result(cv)
     implicit none
     class(matm), intent(inout) :: this
     real(WP),intent(in) :: T
     real(WP) :: TC,cv

     TC = max(3.0_WP,max(T-273.15_WP,200.0_WP))
     cv = (4209.0_WP-1.31_WP*TC+0.014_WP*TC**2)/this%gamm_l

     return
   end function spec_heat_water

   function therm_cond_air(this,T) result(kappa)
     implicit none
     class(matm), intent(inout) :: this
     real(WP),intent(in) :: T
     real(WP) :: TC,kappa

     TC = max(100.0_WP,min(T,1000.0_WP))
     kappa = 2.646e-3_WP*TC**1.5_WP/(TC+245.4_WP*10.0_WP**(-12.0_WP/TC))

     return
   end function therm_cond_air

   function therm_cond_water(this,T) result(kappa)
     implicit none
     class(matm), intent(inout) :: this
     real(WP),intent(in) :: T
     real(WP) :: TC,kappa

     TC = max(1.0_WP,min(T-273.15_WP,200.0_WP))
     kappa = 0.5706_WP+1.75e-3_WP*TC-6.46e-6_WP*TC**2

     return
   end function therm_cond_water

end module matm_class

  
