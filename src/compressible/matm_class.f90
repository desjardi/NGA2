!> Single-phase material modeling class:
!> Provides equations of state and temperature-dependent properties
module matm_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: matm
   
   ! Parameters for viscosity, and heat diffusion models
   integer, parameter, public :: none     =0  !< Sets the constant to 0 (for inviscid, isothermal, etc.). Careful, overwrites input values
   integer, parameter, public :: constant =1  !< Assumes default value, available for all types, overwritten by input
   integer, parameter, public :: water    =2  !< Empirical models for water (mu, kappa, cv)
   integer, parameter, public :: air      =3  !< Empirical models for air (mu, kappa, cv)
   ! More can be added for other materials, alongside functions that feature models
   
   !> Material modeling type intended for single-phase flows (gases or liquids)
   type :: matm
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the structure
      character(len=str_medium) :: name='UNNAMED_MATM'    !< class name (default=UNNAMED_MATM)

      ! Flag for if things vary with temperature or not
      logical  :: const_prop
      ! Flag for if temperature equil of multiphase cells is enforced
      logical  :: mult_iso

      ! Fluid properties
      real(WP) :: gamm,Pref,b,q

      ! Pointers for local thermodynamic or flow variables
      real(WP), dimension(:,:,:), pointer     :: rho,U,V,W,rhoE,P
      
      ! Variables to store chosen diffusion & thermo models (defaults to constant model, M references "model")
      integer  :: M_mu    = 1
      integer  :: M_kappa = 1
      integer  :: M_cv    = 1
      
      ! Defaults for constant diffusion & thermo parameters
      real(WP) :: mu_0  = 1.81e-5_WP
      real(WP) :: kappa_0 = 0.026_WP
      real(WP) :: cv_0    = 717.6_WP
      
      ! Default temperature guess (for the sake of variable cv) upon initialization, can be overwritten
      ! Also used by EOS_gas and EOS_liquid functions if temperature guess is not directly supplied
      real(WP) :: Tdefault = 25.0_WP + 273.15_WP

   contains
      procedure :: EOS_gas, EOS_liquid                    !< Output solver to the screen
      procedure :: EOS_energy                             !< Calculates phase total energy from pressure and kinetic energy
      procedure :: EOS_temp                               !< Calculates phase temperature from pressure directly
      procedure :: EOS_density                            !< Calculates phase density from pressure and temperature directly
      procedure :: EOS_pressure                           !< Calculates phase pressure from density and temperature directly
      procedure :: EOS_all                                !< Calculates vol-avg pressure for entire domain from conserved variables

      procedure :: register_idealgas                      !< EOS available for gas-like fluids
      procedure :: register_stiffenedgas                  !< EOS for liquid-like fluids
      procedure :: register_NobleAbelstiffenedgas         !< EOS for liquid-like fluids
      
      procedure :: register_diffusion_thermo_models       !< Choose models or set parameters for viscosity, thermal conductivity, and specific heat

      procedure :: register_thermoflow_variables          !< Creates pointers from material models to flow solver variables
      
      procedure :: update_temperature                     !< Updates entire array of mixture temperature as passed in
      procedure :: get_local_temperature                  !< Calculates local mixture temperature

      procedure :: viscosity_water,  viscosity_air        !< Empirical models for temperature dependence of dynamic viscosity
      procedure :: spec_heat_water,  spec_heat_air        !< Empirical models for temperature dependence of specific heat (cv)
      procedure :: therm_cond_water, therm_cond_air       !< Empirical models for temperature dependence of thermal conductivity
      
      procedure :: viscosity, spec_heat, therm_cond  !< Provide diffusive/thermal parameters according to registered model
      
   end type matm
   
   
   !> Declare compressible solver constructor
   interface matm
      procedure constructor
   end interface matm
   
contains
   
   
   !> Default constructor for material modeling
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
     self%gamm = 0.0_WP; self%Pref = 0.0_WP; self%q = 0.0_WP; self%b = 0.0_WP

   end function constructor

   subroutine register_idealgas(this,gamm)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: gamm

     this%gamm = gamm
     
   end subroutine register_idealgas

   subroutine register_stiffenedgas(this,gamm,Pref)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: gamm,Pref

     this%gamm = gamm
     this%Pref = Pref
     
   end subroutine register_stiffenedgas

   subroutine register_NobleAbelstiffenedgas(this,gamm,Pref,q,b)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: gamm,Pref,q,b

     this%gamm = gamm
     this%Pref = Pref
     this%q    = q
     this%b    = b

   end subroutine register_NobleAbelstiffenedgas

   subroutine register_diffusion_thermo_models(this,viscmodel,hdffmodel,sphtmodel,viscconst,hdffconst,sphtconst)
     use messager,  only: die
     implicit none
     class(matm), intent(inout) :: this
     integer , intent(in), optional :: viscmodel,hdffmodel,sphtmodel
     real(WP), intent(in), optional :: viscconst,hdffconst,sphtconst

     ! If constants are supplied, replace default constants with provided constants
     if (present(viscconst)) then
       this%mu_0 = viscconst
     end if
     if (present(hdffconst)) then
       this%kappa_0 = hdffconst
     end if
     if (present(sphtconst)) then
       this%cv_0 = sphtconst
     end if
     
     ! If models are supplied, replace default models (constant) with specified models
     ! Specifying a non-constant model will override specifying a constant parameter above
     if (present(viscmodel)) then
       this%M_mu = viscmodel
       ! Check allowed values
       select case (this%M_mu)
       case (none,constant,air,water); ! do nothing
       case default; call die('[matm register_diffusion_thermo_models] Unknown viscosity model')
       end select
     end if
     if (present(hdffmodel)) then
       this%M_kappa = hdffmodel
       ! Check allowed values
       select case (this%M_kappa)
       case (none,constant,air,water); ! do nothing
       case default; call die('[matm register_diffusion_thermo_models] Unknown heat diffusivity model')
       end select
     end if
     if (present(sphtmodel)) then
       this%M_cv = sphtmodel
       ! Check allowed values
       select case (this%M_cv)
       case (constant,air,water); ! do nothing - never zero, so none is not allowed
       case default; call die('[matm register_diffusion_thermo_models] Unknown specific heat model')
       end select
     end if
     
     ! Set parameter values to zero if prescribed
     if (this%M_mu   .eq.none) this%mu_0    = 0.0_WP
     if (this%M_kappa.eq.none) this%kappa_0 = 0.0_WP
     
   end subroutine register_diffusion_thermo_models
   
   subroutine register_thermoflow_variables(this,rho,u,v,w,rhoe,p)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: rho,u,v,w,rhoE,p

     this%rho  => rho
     this%U    => u
     this%V    => v
     this%W    => w
     this%rhoE => rhoE
     this%P    => p

   end subroutine register_thermoflow_variables
   
   
   function viscosity(this,T) result(mu)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: mu
     
     select case(this%M_mu)
     case(none,constant)
        mu = this%mu_0
     case(air)
        mu = this%viscosity_air(T)
     case(water)
        mu = this%viscosity_water(T)
     end select
     
     return
   end function viscosity

   function therm_cond(this,T) result(kappa)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: kappa

     select case(this%M_kappa)
     case(none,constant)
        kappa = this%kappa_0
     case(air)
        kappa = this%therm_cond_air(T)
     case(water)
        kappa = this%therm_cond_water(T)
     end select
     
     return
   end function therm_cond

   function spec_heat(this,T) result(cv)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: cv

     select case(this%M_cv)
     case(none,constant)
        cv = this%cv_0
     case(air)
        cv = this%spec_heat_air(T)
     case(water)
        cv = this%spec_heat_water(T)
     end select

     return
   end function spec_heat

   subroutine update_temperature(this,phase,Temperature)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: Temperature
     integer :: i,j,k
     
     do k=this%cfg%kmino_,this%cfg%kmaxo_
       do j=this%cfg%jmino_,this%cfg%jmaxo_
         do i=this%cfg%imino_,this%cfg%imaxo_
            ! Get local temperature, use current temperature as first guess for cv
            Temperature(i,j,k) = this%get_local_temperature(phase,i,j,k,Temperature(i,j,k))
         end do
       end do
     end do
     
   end subroutine update_temperature
   
   function get_local_temperature(this,phase,i,j,k,Tguess) result(T)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     integer,  intent(in) :: i,j,k
     real(WP), intent(in) :: Tguess
     real(WP) :: Tlast,delta_T,temp,T
     real(WP) :: cv = 0.0_WP
     integer  :: n
     real(WP), parameter :: T_cvg = 1e-1_WP
     integer,  parameter :: n_loop = 10
     
     ! Initialize guess
     temp = Tguess
     select case(trim(adjustl(phase)))
     case('gas')
        ! Gas temperature
        n = 0; delta_T = 1.0_WP + T_cvg
        do while ((delta_T.gt.T_cvg).and.(n.lt.n_loop))
           Tlast = temp
           ! Calculate temperature with latest temperature in cv calculation
           temp = this%EOS_gas(i,j,k,'T',temp)
           ! Increment interation number
           n = n+1
           ! Convergence residual
           delta_T = abs(temp-Tlast)
           ! Exit loop if constant model
           if (this%M_cv.eq.constant) n = n_loop
           ! Exit loop if not converging
           if (temp.le.0.0_WP) then 
              n = n_loop
              temp = Tlast ! hopefully last iteration had better value
           end if
        end do
        ! Update cv value using current temperature value
        cv = this%spec_heat(temp)
     case('liquid')
        ! Liquid temperature
        n = 0; delta_T = 1.0_WP + T_cvg
        do while ((delta_T.gt.T_cvg).and.(n.lt.n_loop))
           Tlast = temp
           ! Calculate temperature with latest temperature in cv calculation
           temp = this%EOS_liquid(i,j,k,'T',temp)
           ! Increment interation number
           n = n+1
           ! Convergence residual
           delta_T = abs(temp-Tlast)
         ! Exit loop if constant model
           if (this%M_cv.eq.constant) n = n_loop
           ! Exit loop if not converging
           if (temp.le.0.0_WP) then 
              n = n_loop
              temp = Tlast ! hopefully last iteration had better value
           end if
        end do
        ! Update cv value using current temperature value
        cv = this%spec_heat(temp)
     end select
     
     ! Return temperature
     T = temp
     
     return
   end function get_local_temperature

   function EOS_liquid(this,i,j,k,flag,temp) result(property)
     implicit none
     class(matm), intent(inout) :: this
     integer,intent(in) :: i,j,k
     character(len=1),intent(in) :: flag
     real(WP), intent(in), optional :: temp ! Just for getting cv
     real(WP) :: property

     select case(flag)
     case('p') ! Pressure
        property = (this%gamm-1.0_WP)/(1.0_WP-this%rho(i,j,k)*this%b)*&
             (this%rhoE(i,j,k)-0.5_WP*this%rho(i,j,k)*&
             (this%U(i,j,k)**2+this%V(i,j,k)**2+this%W(i,j,k)**2) &
             -this%rho(i,j,k)*this%q)-this%gamm*this%Pref
     case('K') ! Bulk modulus
        property = this%P(i,j,k)+this%Pref
        property = this%gamm/(1.0_WP-this%rho(i,j,k)*this%b)*property
     case('T') ! Temperature
        if (present(temp)) then; property = temp; else; property = this%Tdefault; end if
        property = (-(1.0_WP-this%rho(i,j,k)*this%b)*this%Pref+this%rhoE(i,j,k) &
             -0.5_WP*this%rho(i,j,k)*(this%U(i,j,k)**2+this%V(i,j,k)**2+this%W(i,j,k)**2) &
             -this%rho(i,j,k)*this%q) / (this%rho(i,j,k)*this%spec_heat(property))
     case('v') ! (1/rho)(dp/de), for viscous pressure term
        property = (this%gamm-1.0_WP)/(1.0_WP-this%rho(i,j,k)*this%b)
     case('M') ! Bulk modulus calculated from conserved variables
        property = (this%gamm-1.0_WP)/(1.0_WP-this%rho(i,j,k)*this%b)*(this%rhoE(i,j,k) &
             -0.5_WP*this%rho(i,j,k)*(this%U(i,j,k)**2+this%V(i,j,k)**2+this%W(i,j,k)**2) &
             -this%rho(i,j,k)*this%q)-this%gamm*this%Pref
        property = this%gamm/(1.0_WP-this%rho(i,j,k)*this%b)*(property+this%Pref)
     end select

     return
   end function EOS_liquid

   function EOS_gas(this,i,j,k,flag,temp) result(property)
     implicit none
     class(matm), intent(inout) :: this
     integer,intent(in) :: i,j,k
     character(len=1),intent(in) :: flag
     real(WP), intent(in), optional :: temp ! Just for getting cv
     real(WP) :: property

     select case(flag)
     case('p') ! Pressure
        property = (this%gamm-1.0_WP)/(1.0_WP-this%rho(i,j,k)*this%b)*&
             (this%rhoE(i,j,k)-0.5_WP*this%rho(i,j,k)*&
             (this%U(i,j,k)**2+this%V(i,j,k)**2+this%W(i,j,k)**2) &
             -this%rho(i,j,k)*this%q)-this%gamm*this%Pref
     case('K') ! Bulk modulus
        property = this%P(i,j,k)+this%Pref
        property = this%gamm/(1.0_WP-this%rho(i,j,k)*this%b)*property
     case('T') ! Temperature
       if (present(temp)) then; property = temp; else; property = this%Tdefault; end if
        property = (-(1.0_WP-this%rho(i,j,k)*this%b)*this%Pref+this%rhoE(i,j,k) &
             -0.5_WP*this%rho(i,j,k)*(this%U(i,j,k)**2+this%V(i,j,k)**2+this%W(i,j,k)**2) &
             -this%rho(i,j,k)*this%q) / (this%rho(i,j,k)*this%spec_heat(property))
     case('v') ! (1/rho)(dp/de), for viscous pressure term
        property = (this%gamm-1.0_WP)/(1.0_WP-this%rho(i,j,k)*this%b)
     case('M') ! Bulk modulus calculated from conserved variables
        property = (this%gamm-1.0_WP)/(1.0_WP-this%rho(i,j,k)*this%b)*(this%rhoE(i,j,k) &
             -0.5_WP*this%rho(i,j,k)*(this%U(i,j,k)**2+this%V(i,j,k)**2+this%W(i,j,k)**2) &
             -this%rho(i,j,k)*this%q)-this%gamm*this%Pref
        property = this%gamm/(1.0_WP-this%rho(i,j,k)*this%b)*(property+this%Pref)
     end select

     return
   end function EOS_gas

   function EOS_energy(this,pres,dens,uvel,vvel,wvel) result(energy)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: pres,dens,uvel,vvel,wvel
     real(WP) :: KE,energy

     KE = 0.5_WP*dens*(uvel**2+vvel**2+wvel**2)
     energy = (pres + this%gamm*this%Pref)*(1.0_WP-dens*this%b)/(this%gamm-1.0_WP) + dens*this%q + KE

     return
   end function EOS_energy

   function EOS_temp(this,pres,dens,cv) result(temp)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: pres,dens,cv
     real(WP) :: temp

     temp = (1.0_WP-dens*this%b)*(pres+this%Pref)/(this%gamm-1.0_WP)/(dens*cv)

     return
   end function EOS_temp
   
   function EOS_pressure(this,temp,dens,cv) result(pres)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: temp,dens,cv
     real(WP) :: pres

     pres = temp*(dens*cv)*(this%gamm-1.0_WP)/(1.0_WP-dens*this%b)-this%Pref

     return
   end function EOS_pressure
   
   function EOS_density(this,pres,temp,cv) result(dens)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: pres,temp,cv
     real(WP) :: dens

     dens = (pres+this%Pref)/(cv*(this%gamm-1.0_WP)*temp+this%b*(pres+this%Pref))

     return
   end function EOS_density

   function EOS_all(this) result(buf_P)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: buf_P

     buf_P = (this%gamm-1.0_WP)/(1.0_WP-this%rho*this%b)*(this%rhoE-0.5_WP* &
          this%rho*(this%U**2+this%V**2+this%W**2)-this%rho*this%q)-this%gamm*this%Pref

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
     cv = (4209.0_WP-1.31_WP*TC+0.014_WP*TC**2)/this%gamm

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

  
