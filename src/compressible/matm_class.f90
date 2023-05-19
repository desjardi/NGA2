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
   
   ! Parameters for viscosity, and heat diffusion models
   integer, parameter, public :: none     =0  !< Sets the constant to 0 (for inviscid, isothermal, etc.). Careful, overwrites input values
   integer, parameter, public :: constant =1  !< Assumes default value, available for all types, overwritten by input
   integer, parameter, public :: water    =2  !< Empirical models for water (mu, kappa, cv)
   integer, parameter, public :: air      =3  !< Empirical models for air (mu, kappa, cv)
   ! More can be added for other materials, alongside functions that feature models
   
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

      ! Pointers for local thermodynamic or flow variables
      real(WP), dimension(:,:,:), pointer     :: Grho,GU,GV,GW,GrhoE,GP
      real(WP), dimension(:,:,:), pointer     :: Lrho,LU,LV,LW,LrhoE,LP
      
      ! Variables to store chosen diffusion & thermo models (defaults to constant model, M references "model")
      integer  :: M_mu_l    = 1
      integer  :: M_kappa_l = 1
      integer  :: M_cv_l    = 1
      integer  :: M_mu_g    = 1
      integer  :: M_kappa_g = 1
      integer  :: M_cv_g    = 1
      
      ! Defaults for constant diffusion & thermo parameters
      real(WP) :: mu_l0  = 8.9e-4_WP
      real(WP) :: mu_g0  = 1.81e-5_WP
      real(WP) :: kappa_l0 = 0.610_WP
      real(WP) :: kappa_g0 = 0.026_WP
      real(WP) :: cv_l0    = 4130.2_WP
      real(WP) :: cv_g0    = 717.6_WP
      
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

      procedure :: EOS_relax_quad                         !< Outputs terms for mechanical pressure relaxation step
      procedure :: EOS_thermal_relax_quad                 !< Outputs terms for thermo-mechanical pressure relaxation step
      procedure :: bulkmod_intf                           !< Calculate the bulk modulus at liquid-gas interface
      procedure :: fix_energy                             !< Can correct erroneous energy values
      
      procedure :: update_temperature                     !< Updates entire array of mixture temperature as passed in
      procedure :: get_local_temperature                  !< Calculates local mixture temperature

      procedure :: viscosity_water,  viscosity_air        !< Empirical models for temperature dependence of dynamic viscosity
      procedure :: spec_heat_water,  spec_heat_air        !< Empirical models for temperature dependence of specific heat (cv)
      procedure :: therm_cond_water, therm_cond_air       !< Empirical models for temperature dependence of thermal conductivity
      
      procedure :: viscosity_gas, spec_heat_gas, therm_cond_gas  !< Provide diffusive/thermal parameters according to registered model
      procedure :: viscosity_liquid, spec_heat_liquid, therm_cond_liquid
      
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

   subroutine register_diffusion_thermo_models(this,viscmodel_liquid,hdffmodel_liquid,sphtmodel_liquid,viscmodel_gas,hdffmodel_gas,sphtmodel_gas,viscconst_liquid,hdffconst_liquid,sphtconst_liquid,viscconst_gas,hdffconst_gas,sphtconst_gas)
     use messager,  only: die
     implicit none
     class(matm), intent(inout) :: this
     integer , intent(in), optional :: viscmodel_liquid,hdffmodel_liquid,sphtmodel_liquid,viscmodel_gas,hdffmodel_gas,sphtmodel_gas
     real(WP), intent(in), optional :: viscconst_liquid,hdffconst_liquid,sphtconst_liquid,viscconst_gas,hdffconst_gas,sphtconst_gas

     ! If constants are supplied, replace default constants with provided constants
     if (present(viscconst_liquid)) then
       this%mu_l0 = viscconst_liquid
     end if
     if (present(viscconst_gas)) then
       this%mu_g0 = viscconst_gas
     end if
     if (present(hdffconst_liquid)) then
       this%kappa_l0 = hdffconst_liquid
     end if
     if (present(hdffconst_gas)) then
       this%kappa_g0 = hdffconst_gas
     end if
     if (present(sphtconst_liquid)) then
       this%cv_l0 = sphtconst_liquid
     end if
     if (present(sphtconst_gas)) then
       this%cv_g0 = sphtconst_gas
     end if
     
     ! If models are supplied, replace default models (constant) with specified models
     ! Specifying a non-constant model will override specifying a constant parameter above
     if (present(viscmodel_liquid)) then
       this%M_mu_l = viscmodel_liquid
       ! Check allowed values
       select case (this%M_mu_l)
       case (none,constant,water); ! do nothing
       case default; call die('[matm register_diffusion_thermo_models] Unknown liquid viscosity model')
       end select
     end if
     if (present(viscmodel_gas)) then
       this%M_mu_g = viscmodel_gas
       ! Check allowed values
       select case (this%M_mu_g)
       case (none,constant,air); ! do nothing
       case default; call die('[matm register_diffusion_thermo_models] Unknown gas viscosity model')
       end select
     end if
     if (present(hdffmodel_liquid)) then
       this%M_kappa_l = hdffmodel_liquid
       ! Check allowed values
       select case (this%M_kappa_l)
       case (none,constant,water); ! do nothing
       case default; call die('[matm register_diffusion_thermo_models] Unknown liquid heat diffusivity model')
       end select
     end if
     if (present(hdffmodel_gas)) then
       this%M_kappa_g = hdffmodel_gas
       ! Check allowed values
       select case (this%M_kappa_g)
       case (none,constant,air); ! do nothing
       case default; call die('[matm register_diffusion_thermo_models] Unknown gas heat diffusivity model')
       end select
     end if
     if (present(sphtmodel_liquid)) then
       this%M_cv_l = sphtmodel_liquid
       ! Check allowed values
       select case (this%M_cv_l)
       case (constant,water); ! do nothing - never zero, so none is not allowed
       case default; call die('[matm register_diffusion_thermo_models] Unknown liquid specific heat model')
       end select
     end if
     if (present(sphtmodel_gas)) then
       this%M_cv_g = sphtmodel_gas
       ! Check allowed values
       select case (this%M_cv_g)
       case (constant,air); ! do nothing - never zero, so none is not allowed
       case default; call die('[matm register_diffusion_thermo_models] Unknown gas specific heat model')
       end select
     end if
     
     ! Set parameter values to zero if prescribed
     if (this%M_mu_l   .eq.none) this%mu_l0    = 0.0_WP
     if (this%M_kappa_l.eq.none) this%kappa_l0 = 0.0_WP
     if (this%M_mu_g   .eq.none) this%mu_g0    = 0.0_WP
     if (this%M_kappa_g.eq.none) this%kappa_g0 = 0.0_WP
     
   end subroutine register_diffusion_thermo_models
   
   subroutine register_thermoflow_variables(this,phase,rho,u,v,w,rhoe,p)
     implicit none
     class(matm), intent(inout) :: this
     character(len=*),intent(in) :: phase
     real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: rho,u,v,w,rhoE,p
     
     select case(trim(adjustl(phase)))
     case('liquid')
        this%Lrho  => rho
        this%LU    => u
        this%LV    => v
        this%LW    => w
        this%LrhoE => rhoE
        this%LP    => p
     case('gas')
        this%Grho  => rho
        this%GU    => u
        this%GV    => v
        this%GW    => w
        this%GrhoE => rhoE
        this%GP    => p
     end select

   end subroutine register_thermoflow_variables
   
   
   function viscosity_liquid(this,T) result(mu)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: mu
     
     select case(this%M_mu_l)
     case(none,constant)
       mu = this%mu_l0
     case(water)
       mu = this%viscosity_water(T)
     end select
     
     return
   end function viscosity_liquid
   
   function viscosity_gas(this,T) result(mu)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: mu
     
     select case(this%M_mu_g)
     case(none,constant)
       mu = this%mu_g0
     case(air)
       mu = this%viscosity_air(T)
     end select
     
     return
   end function viscosity_gas
   
   function therm_cond_liquid(this,T) result(kappa)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: kappa
     
     select case(this%M_kappa_l)
     case(none,constant)
       kappa = this%kappa_l0
     case(water)
       kappa = this%therm_cond_water(T)
     end select
     
     return
   end function therm_cond_liquid
   
   function therm_cond_gas(this,T) result(kappa)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: kappa
     
     select case(this%M_kappa_g)
     case(none,constant)
       kappa = this%kappa_g0
     case(air)
       kappa = this%therm_cond_air(T)
     end select
     
     return
   end function therm_cond_gas
   
   function spec_heat_liquid(this,T) result(cv)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: cv
     
     select case(this%M_cv_l)
     case(constant)
       cv = this%cv_l0
     case(water)
       cv = this%spec_heat_water(T)
     end select
     
     return
   end function spec_heat_liquid
   
   function spec_heat_gas(this,T) result(cv)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: T  ! Temperature - always supplied, not always used
     real(WP) :: cv
     
     select case(this%M_cv_g)
     case(constant)
       cv = this%cv_g0
     case(air)
       cv = this%spec_heat_air(T)
     end select
     
     return
   end function spec_heat_gas
   
   subroutine update_temperature(this,vf,Temperature)
     use vfs_class, only : vfs
     implicit none
     class(matm), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: Temperature
     integer :: i,j,k
     
     do k=this%cfg%kmino_,this%cfg%kmaxo_
       do j=this%cfg%jmino_,this%cfg%jmaxo_
         do i=this%cfg%imino_,this%cfg%imaxo_
            ! Get local temperature, use current temperature as first guess for cv
            Temperature(i,j,k) = this%get_local_temperature(vf,i,j,k,Temperature(i,j,k))
         end do
       end do
     end do
     
   end subroutine update_temperature
   
   function get_local_temperature(this,vf,i,j,k,Tguess) result(T)
     use vfs_class, only : vfs, VFhi, VFlo
     implicit none
     class(matm), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     integer,  intent(in) :: i,j,k
     real(WP), intent(in) :: Tguess
     real(WP) :: Tlast,delta_T,Gtemp,Ltemp,T
     real(WP) :: cv_g = 0.0_WP
     real(WP) :: cv_l = 0.0_WP
     integer  :: n
     real(WP), parameter :: T_cvg = 1e-1_WP
     integer,  parameter :: n_loop = 10
     
     ! Initialize guesses
     Gtemp = Tguess; Ltemp = Tguess
     ! Gas temperature
     if (vf%VF(i,j,k).le.VFhi) then
       n = 0; delta_T = 1.0_WP + T_cvg
       do while ((delta_T.gt.T_cvg).and.(n.lt.n_loop))
         Tlast = Gtemp
         ! Calculate temperature with latest temperature in cv_g calculation
         Gtemp = this%EOS_gas(i,j,k,'T',Gtemp)
         ! Increment interation number
         n = n+1
         ! Convergence residual
         delta_T = abs(Gtemp-Tlast)
         ! Exit loop if constant model
         if (this%M_cv_g.eq.constant) n = n_loop
         ! Exit loop if not converging
         if (Gtemp.le.0.0_WP) then 
           n = n_loop
           Gtemp = Tlast ! hopefully last iteration had better value
         end if
       end do
       ! Update cv value using current temperature value
       cv_g = this%spec_heat_gas(Gtemp)
     end if
     ! Liquid temperature
     if (vf%VF(i,j,k).ge.VFlo) then
       n = 0; delta_T = 1.0_WP + T_cvg
       do while ((delta_T.gt.T_cvg).and.(n.lt.n_loop))
         Tlast = Ltemp
         ! Calculate temperature with latest temperature in cv_l calculation
         Ltemp = this%EOS_liquid(i,j,k,'T',Ltemp)
         ! Increment interation number
         n = n+1
         ! Convergence residual
         delta_T = abs(Ltemp-Tlast)
         ! Exit loop if constant model
         if (this%M_cv_l.eq.constant) n = n_loop
         ! Exit loop if not converging
         if (Ltemp.le.0.0_WP) then 
           n = n_loop
           Ltemp = Tlast ! hopefully last iteration had better value
         end if
       end do
       ! Update cv value using current temperature value
       cv_l = this%spec_heat_liquid(Ltemp)
     end if
     
     ! Mixture temperature
     T = ((1.0_WP-vf%VF(i,j,k))*cv_g*this%Grho(i,j,k)*Gtemp &
        +(        vf%VF(i,j,k))*cv_l*this%Lrho(i,j,k)*Ltemp)&
        /((1.0_WP-vf%VF(i,j,k))*cv_g*this%Grho(i,j,k) &
        +(        vf%VF(i,j,k))*cv_l*this%Lrho(i,j,k) )
     
     return
   end function get_local_temperature

   function EOS_liquid(this,i,j,k,flag,Ltemp) result(property)
     implicit none
     class(matm), intent(inout) :: this
     integer,intent(in) :: i,j,k
     character(len=1),intent(in) :: flag
     real(WP), intent(in), optional :: Ltemp ! Just for getting cv
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
        if (present(Ltemp)) then; property = Ltemp; else; property = this%Tdefault; end if
        property = (-(1.0_WP-this%Lrho(i,j,k)*this%b_l)*this%Pref_l+this%LrhoE(i,j,k) &
             -0.5_WP*this%Lrho(i,j,k)*(this%LU(i,j,k)**2+this%LV(i,j,k)**2+this%LW(i,j,k)**2) &
             -this%Lrho(i,j,k)*this%q_l) / (this%Lrho(i,j,k)*this%spec_heat_liquid(property))
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

   function EOS_gas(this,i,j,k,flag,Gtemp) result(property)
     implicit none
     class(matm), intent(inout) :: this
     integer,intent(in) :: i,j,k
     character(len=1),intent(in) :: flag
     real(WP), intent(in), optional :: Gtemp ! Just for getting cv
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
       if (present(Gtemp)) then; property = Gtemp; else; property = this%Tdefault; end if
        property = (-(1.0_WP-this%Grho(i,j,k)*this%b_g)*this%Pref_g+this%GrhoE(i,j,k) &
             -0.5_WP*this%Grho(i,j,k)*(this%GU(i,j,k)**2+this%GV(i,j,k)**2+this%GW(i,j,k)**2) &
             -this%Grho(i,j,k)*this%q_g) / (this%Grho(i,j,k)*this%spec_heat_gas(property))
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
   
   function EOS_pressure(this,temp,dens,cv,phase) result(pres)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: temp,dens,cv
     character(len=*),intent(in) :: phase
     real(WP) :: pres

     select case(trim(adjustl(phase)))
     case('liquid')
        pres = temp*(dens*cv)*(this%gamm_l-1.0_WP)/(1.0_WP-dens*this%b_l)-this%Pref_l
     case('gas')
        pres = temp*(dens*cv)*(this%gamm_g-1.0_WP)/(1.0_WP-dens*this%b_g)-this%Pref_g
     end select

     return
   end function EOS_pressure
   
   function EOS_density(this,pres,temp,cv,phase) result(dens)
     implicit none
     class(matm), intent(inout) :: this
     real(WP), intent(in) :: pres,temp,cv
     character(len=*),intent(in) :: phase
     real(WP) :: dens

     select case(trim(adjustl(phase)))
     case('liquid')
        dens = (pres+this%Pref_l)/(cv*(this%gamm_l-1.0_WP)*temp+this%b_l*(pres+this%Pref_l))
     case('gas')
        dens = (pres+this%Pref_g)/(cv*(this%gamm_g-1.0_WP)*temp+this%b_g*(pres+this%Pref_g))
     end select

     return
   end function EOS_density

   function EOS_all(this,vf) result(buf_P)
     use vfs_class, only : vfs
     implicit none
     class(matm), intent(inout) :: this
     class(vfs),  intent(inout) :: vf
     real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: buf_P

     buf_P = (1.0_WP-vf%VF)*((this%gamm_g-1.0_WP)/(1.0_WP-this%Grho*this%b_g)*(this%GrhoE-0.5_WP* &
             this%Grho*(this%GU**2+this%GV**2+this%GW**2)-this%Grho*this%q_g)-this%gamm_g*this%Pref_g) + &
             (       vf%VF)*((this%gamm_l-1.0_WP)/(1.0_WP-this%Lrho*this%b_l)*(this%LrhoE-0.5_WP* &
             this%Lrho*(this%LU**2+this%LV**2+this%LW**2)-this%Lrho*this%q_l)-this%gamm_l*this%Pref_l)

     return
   end function EOS_all

   subroutine EOS_relax_quad(this,i,j,k,myVF,srcVF,my_pjump,lbm,gbm,VF_terms,quad_terms,Pint_terms)
     implicit none
     class(matm), intent(inout) :: this
     integer,  intent(in) :: i,j,k
     real(WP), intent(in) :: myVF,srcVF,my_pjump,lbm,gbm
     real(WP) :: KE, cJ, Pint
     real(WP) :: n1,n0,d1,d0
     real(WP) :: re1,r1,a1,b1,q1,g1,pr1,z1,p1
     real(WP) :: re2,r2,a2,b2,q2,g2,pr2,z2,p2
     real(WP) :: a,b,c
     real(WP), parameter :: phist = 1.0_WP
     real(WP), parameter :: phi0 = 1.0_WP-phist
     real(WP), dimension(4) :: VF_terms
     real(WP), dimension(3) :: quad_terms
     real(WP), dimension(2) :: Pint_terms

     ! Incoming quantities
     ! srcVF   : predicted change in phase volume from the advection solver
     ! myVF    : local value of volume fraction
     ! my_pjump: pressure jump from phase 1 to 2, i.e. liquid to gas
     ! lbm,gbm : current phase bulk moduli calculated in adv. step, used as backup

     ! Working quantities
     ! phist   : temporal weighting coefficient for pressure (phi^*)
     ! phi0    : corresponding weighting coefficient, at timestep 0 (unrelaxed)
     ! r1 ,r2  : phase mass for liquid, gas
     ! re1,re2 : phase internal energy for liquid, gas
     ! a1, a2  : phase volume fraction without srcVF
     ! Pint    : interface pressure (weighted with acoustic impedance) at timestep 0
     ! cJ      : part of impedance weighting applied to pressure jump

     ! Outputs
     ! a,b,c       : coefficients for quadratic equation of phase 1 eq. pressure (quad_terms)
     ! n1,n0,d1,d0 : coefficients for rational equation of phase 1 vol. frac.    (VF_terms)
     ! phist, (-phist*cJ*my_pjump + phi0*Pint):
     !   ^       ^ : coefficients for interface pressure in energy exchange      (Pint_terms)

     ! Calculate kinetic energy without the density (assuming single velocity)
     KE  = 0.5_WP*(this%LU(i,j,k)**2+this%LV(i,j,k)**2+this%LW(i,j,k)**2)

     ! Liquid is fluid 1, Gas is fluid 2
     re1 = (       myVF)*(this%LrhoE(i,j,k)-this%Lrho(i,j,k)*KE)
     re2 = (1.0_WP-myVF)*(this%GrhoE(i,j,k)-this%Grho(i,j,k)*KE)
     r1  = (       myVF)*this%Lrho (i,j,k)
     r2  = (1.0_WP-myVF)*this%Grho (i,j,k)
     a1  = myVF-srcVF; a2  = 1.0_WP-a1

     ! Calculate acoustic impedance (rho c) and initial interface pressure
     ! Adjust quantities first to get to the 0 instance (without the source term)
     this%LrhoE(i,j,k) = (       myVF)*this%LrhoE(i,j,k)/a1
     this%GrhoE(i,j,k) = (1.0_WP-myVF)*this%GrhoE(i,j,k)/a2
     this%Lrho (i,j,k) = (       myVF)*this%Lrho (i,j,k)/a1
     this%Grho (i,j,k) = (1.0_WP-myVF)*this%Grho (i,j,k)/a2
     ! Use the bulk modulus calculated from conservative variables unless it is faulty
     z1  = this%EOS_liquid(i,j,k,'M'); if (z1.lt.0.0_WP) z1 = lbm
     z2  = this%EOS_gas   (i,j,k,'M'); if (z2.lt.0.0_WP) z2 = gbm
     z1  = sqrt(this%Lrho(i,j,k)*z1)
     z2  = sqrt(this%Grho(i,j,k)*z2)
     p1  = this%EOS_liquid(i,j,k,'p')
     p2  = this%EOS_gas   (i,j,k,'p')
     Pint = (z2*p1+z1*p2)/(z1+z2)

     ! Weight based on impedance to get interface pressure after relaxation
     cJ = z1/(z1+z2)

     ! Adjust back
     this%LrhoE(i,j,k) = a1*this%LrhoE(i,j,k)/(       myVF)
     this%GrhoE(i,j,k) = a2*this%GrhoE(i,j,k)/(1.0_WP-myVF)
     this%Lrho (i,j,k) = a1*this%Lrho (i,j,k)/(       myVF)
     this%Grho (i,j,k) = a2*this%Grho (i,j,k)/(1.0_WP-myVF)

     ! Temporal weighting of pressure
     ! Rename property variables for brevity
     b1  = this%b_l; b2 = this%b_g;  g1 = this%gamm_l;  g2 = this%gamm_g
     q1  = this%q_l; q2 = this%q_g;  pr1 = this%Pref_l; pr2 = this%Pref_g
     ! Intermediate terms
     n1 = a1*phist + r1*b1/(g1-1.0_WP)
     n0 = a1*(phi0*Pint - phist*cJ*my_pjump) + re1 - r1*q1 + g1/(g1-1.0_WP)*pr1*r1*b1
     d1 = phist + 1.0_WP/(g1-1.0_WP)
     d0 = phi0*Pint - phist*cJ*my_pjump + g1/(g1-1.0_WP)*pr1
     ! Terms in quadratic equation
     a = d1*((1.0_WP - r2*b2)/(g2-1.0_WP) + phist*a1) &
       + n1*(-1.0_WP/(g2-1.0_WP) - phist)
     b = d1*((g2*pr2-my_pjump)/(g2-1.0_WP)*(1.0_WP - r2*b2) + r2*q2 - re2 + a1*(phi0*Pint-phist*cJ*my_pjump)) &
       + n1*(-(g2*pr2-my_pjump)/(g2-1.0_WP) - phi0*Pint + phist*cJ*my_pjump) &
       + d0*((1.0_WP - r2*b2)/(g2-1.0_WP) + phist*a1) &
       + n0*(-1.0_WP/(g2-1.0_WP) - phist)
     c = d0*((g2*pr2-my_pjump)/(g2-1.0_WP)*(1.0_WP - r2*b2) + r2*q2 - re2 + a1*(phi0*Pint-phist*cJ*my_pjump)) &
          + n0*(-(g2*pr2-my_pjump)/(g2-1.0_WP) - phi0*Pint + phist*cJ*my_pjump)
     ! Put into vector -> for result of Peq
     quad_terms = (/a,b,c/)
     ! Terms in new volume fraction equations -> for VF=(n1*Peq+n0)/(d1*Peq+d0)
     VF_terms   = (/n1,n0,d1,d0/)
     ! Terms in interface pressure equation -> for Pint = phist*(Peq-cJ*my_pjump)+phi0*Pint)
     Pint_terms = (/phist,-phist*cJ*my_pjump+phi0*Pint/)
     ! ^ uses acoustic impedance to calculate the interface pressure in relaxed state,
     !   and time interpolates with the value of interface pressure at initial state

     return
   end subroutine EOS_relax_quad
   
   subroutine EOS_thermal_relax_quad(this,i,j,k,myVF,my_pjump,cv_l,cv_g,VF_terms,quad_terms)
     implicit none
     class(matm), intent(inout) :: this
     integer,  intent(in) :: i,j,k
     real(WP), intent(in) :: my_pjump,myVF,cv_l,cv_g
     real(WP) :: KE
     real(WP) :: n1,n0,d1,d0
     real(WP) :: re1,r1,cv1,b1,q1,g1,pr1,temp1
     real(WP) :: re2,r2,cv2,b2,q2,g2,pr2,temp2
     real(WP) :: a,b,c
     real(WP), parameter :: phist = 1.0_WP
     real(WP), parameter :: phi0 = 1.0_WP-phist
     real(WP), dimension(4) :: VF_terms
     real(WP), dimension(3) :: quad_terms
     real(WP), dimension(2) :: Pint_terms

     ! Incoming quantities
     ! myVF    : local value of volume fraction
     ! my_pjump: pressure jump from phase 1 to 2, i.e. liquid to gas
     ! cv_l    : specific heat of liquid phase
     ! cv_g    : specific heat of gas phase

     ! Outputs
     ! a,b,c       : coefficients for quadratic equation of phase 1 eq. pressure (quad_terms)
     ! n1,n0,d1,d0 : coefficients for rational equation of phase 1 vol. frac.    (VF_terms)

     ! Calculate kinetic energy without the density (assuming single velocity)
     KE  = 0.5_WP*(this%LU(i,j,k)**2+this%LV(i,j,k)**2+this%LW(i,j,k)**2)

     ! Liquid is fluid 1, Gas is fluid 2
     re1 = (       myVF)*(this%LrhoE(i,j,k)-this%Lrho(i,j,k)*KE)
     re2 = (1.0_WP-myVF)*(this%GrhoE(i,j,k)-this%Grho(i,j,k)*KE)
     r1  = (       myVF)*this%Lrho (i,j,k)
     r2  = (1.0_WP-myVF)*this%Grho (i,j,k)
     ! Rename property variables for brevity
     b1  = this%b_l; b2 = this%b_g;  g1 = this%gamm_l;  g2 = this%gamm_g
     q1  = this%q_l; q2 = this%q_g;  pr1 = this%Pref_l; pr2 = this%Pref_g
     cv1 = cv_l; cv2 = cv_g;
     ! Intermediate terms
     temp1 = (g1-1.0_WP)*r1*cv1
     temp2 = (g2-1.0_WP)*r2*cv2
     n1 =     r1*b1*temp2 +                (1.0_WP-r2*b2)*temp1
     n0 = pr1*r1*b1*temp2 + (pr2-my_pjump)*(1.0_WP-r2*b2)*temp1
     d1 =           temp2 +                               temp1
     d0 =       pr1*temp2 +                (pr2-my_pjump)*temp1
     ! Terms in quadratic equation
     a = n1*(1.0_WP/(g1-1.0_WP) - 1.0_WP/(g2-1.0_WP)) + &
          d1*(-r1*b1/(g1-1.0_WP) + (1.0_WP-r2*b2)/(g2-1.0_WP))
     b = n1*(g1*pr1/(g1-1.0_WP) - (g2*pr2-my_pjump)/(g2-1.0_WP)) + &
          d1*(-g1*pr1*r1*b1/(g1-1.0_WP) + (g2*pr2-my_pjump)*(1.0_WP-r2*b2)/(g2-1.0_WP) &
          + r1*q1 + r2*q2 - re1 - re2) + &
          n0*(1.0_WP/(g1-1.0_WP) - 1.0_WP/(g2-1.0_WP)) + &
          d0*(-r1*b1/(g1-1.0_WP) + (1.0_WP-r2*b2)/(g2-1.0_WP))
     c = n0*(g1*pr1/(g1-1.0_WP) - (g2*pr2-my_pjump)/(g2-1.0_WP)) + &
          d0*(-g1*pr1*r1*b1/(g1-1.0_WP) + (g2*pr2-my_pjump)*(1.0_WP-r2*b2)/(g2-1.0_WP) &
          + r1*q1 + r2*q2 - re1 - re2)
     ! Put into vector -> for result of Peq
     quad_terms = (/a,b,c/)
     ! Terms in new volume fraction equations -> for VF=(n1*Peq+n0)/(d1*Peq+d0)
     VF_terms   = (/n1,n0,d1,d0/)

     return
   end subroutine EOS_thermal_relax_quad

   subroutine bulkmod_intf(this,i,j,k,rc2_l,rc2_g)
     implicit none
     class(matm), intent(in) :: this
     real(WP), intent(inout) :: rc2_l,rc2_g
     integer  :: i,j,k
     real(WP) :: p_I

     ! Calculate interface pressure (using acoustic impedance rule)
     p_I = (sqrt(this%Grho(i,j,k)*rc2_g)*this%LP(i,j,k) &
           +sqrt(this%Lrho(i,j,k)*rc2_l)*this%GP(i,j,k)) /&
           (sqrt(this%Grho(i,j,k)*rc2_g)+sqrt(this%Lrho(i,j,k)*rc2_l))

     ! Calculate each bulk modulus at the interface
     rc2_l = (p_I*(this%gamm_l-1.0_WP)+this%LP(i,j,k) &
          +this%gamm_l*this%Pref_l)/(1.0_WP-this%Lrho(i,j,k)*this%b_l)
     rc2_g = (p_I*(this%gamm_g-1.0_WP)+this%GP(i,j,k) &
          +this%gamm_g*this%Pref_g)/(1.0_WP-this%Grho(i,j,k)*this%b_g)

     return
   end subroutine bulkmod_intf

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

   subroutine fix_energy(this,vf,KE,PresH,my_pjump,i,j,k,Gflag,Lflag)
     use vfs_class, only : VFlo, VFhi
     implicit none
     class(matm), intent(inout) :: this
     real(WP),    intent(in) :: vf, KE, PresH, my_pjump
     integer,     intent(in)  :: i,j,k
     logical,     intent(inout)  :: Gflag, Lflag
     real(WP), parameter :: Lrhoe_fix = 1.0_WP
     real(WP), parameter :: Grhoe_fix = 1.0_WP
     real(WP), parameter :: P_fix = 1e-9_WP

     if (vf.lt.VFhi.and. &
          min(this%EOS_gas   (i,j,k,'p')+this%Pref_g,this%GrhoE(i,j,k)-this%Grho(i,j,k)*KE).le.0.0_WP) then
        this%GrhoE(i,j,k) = max(Grhoe_fix+this%Grho(i,j,k)*KE,this%EOS_energy(max( &
             P_fix-this%Pref_g, &
             PresH-(       vf)*my_pjump &
             ), this%Grho(i,j,k),sqrt(2.0_WP*KE),0.0_WP,0.0_WP,'gas'))
        Gflag = .true.
     end if
     if (vf.gt.VFlo.and. &
          min(this%EOS_liquid(i,j,k,'p')+this%Pref_l,this%LrhoE(i,j,k)-this%Lrho(i,j,k)*KE).le.0.0_WP) then
        this%LrhoE(i,j,k) = max(Lrhoe_fix+this%Lrho(i,j,k)*KE,this%EOS_energy(max( &
             P_fix-this%Pref_l, &
             PresH+(1.0_WP-vf)*my_pjump &
             ), this%Lrho(i,j,k),sqrt(2.0_WP*KE),0.0_WP,0.0_WP,'liquid'))
        Lflag = .true.
     end if

     return
   end subroutine fix_energy

end module matm_class

  
