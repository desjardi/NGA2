!> NASG-liquid + ideal-gas-mixture relaxation (extends SG-IG mixture with co-volume b corrections)
!> Inherits pTg_relax unchanged from parent; virtual dispatch into the overridden helpers
!> (relax_p, pT_relax, get_T_lvg, get_coeffs_lv, get_p_eq) does the right thing.
!> References Pelanti 2022 https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
module relax_igmix_nasg_class
   use precision,               only: WP
   use relax_igmix_sg_class,    only: relax_igmix_sg
   use stiffened_gas_class,     only: stiffened_gas
   use nasg_class,              only: nasg
   use igmix_class,             only: igmix
   implicit none
   private

   public :: relax_igmix_nasg

   type, extends(relax_igmix_sg) :: relax_igmix_nasg
      !> Typed pointer for direct b access (avoids select type in hot path)
      type(nasg), pointer :: liq_nasg=>null()
   contains
      procedure :: initialize
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: get_T_lvg
      procedure :: get_coeffs_lv
      procedure :: get_p_eq
   end type relax_igmix_nasg

contains

   !> Initialize: call parent (sets liq, gas, indV/A, AS-DS, ES=0), then set NASG-specific fields
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_nasg),       intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      real(WP) :: cpV,cvV,RV
      ! Parent: sets this%liq=>liq, this%gas=>gas, indV/A, saturation coeffs AS-DS, ES=0
      call this%relax_igmix_sg%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      ! Set typed pointer and override ES (select type to extract type(nasg) from class(stiffened_gas))
      select type (liq)
      type is (nasg)
         this%liq_nasg=>liq
         cpV=gas%cp(indV); cvV=gas%cv(indV); RV=cpV-cvV
         this%ES=liq%b/RV
      end select
   end subroutine initialize

   !> Mechanical relaxation (Pelanti 2022 with co-volume b)
   subroutine p_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_nasg), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: Q0,y
      real(WP) :: RHOL,RHOG,CL,CG,GL,GG,PL,PG,IL,IG,ZL,ZG,Pint
      real(WP) :: xiL,xiG,xiLinv,xiGinv
      real(WP) :: VFeq,VF0,Peq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      integer  :: iVQ
      iVQ=7+this%liq%ns+this%indV-1
      ! Store input
      allocate(Q0(size(Q))); VF0=VF; Q0=Q
      ! Vapor mass fraction and gas composition
      if (Q(2).gt.0.0_WP) then
         Yv=Q(iVQ)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Phasic thermodynamic quantities
      RHOL=Q(1)/(       VF); RHOG=Q(2)/(1.0_WP-VF)
      IL=Q(3)/Q(1);          IG=Q(4)/Q(2)
      PL=this%liq%get_p_from_rho_e(rho=RHOL,e=IL,y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=PL,rho=RHOL,y=[1.0_WP])
      CG=this%gas%get_c_from_p_rho(p=PG,rho=RHOG,y=y)
      ! Hard clipping
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(iVQ)=Yv*Q(2)
         call dealloc(); return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(iVQ)=0.0_WP
         call dealloc(); return
      end if
      ! Phasic impedances
      ZL=Q(1)/(       VF)*CL
      ZG=Q(2)/(1.0_WP-VF)*CG
      ! Interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! ODE coefficients with Gruneisen
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL,y=[1.0_WP])
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      xiL=         VF/(GL*(Pint-PL)+RHOL*CL**2)
      xiG=(1.0_WP-VF)/(GG*(Pint-PG)+RHOG*CG**2)
      xiLinv=1.0_WP/xiL; xiGinv=1.0_WP/xiG
      ! Equilibrium volume fraction
      VFeq=VF-(PG-PL)/(xiLinv+xiGinv)
      if ((VFeq.lt.0.0_WP).or.(VFeq.gt.1.0_WP)) then
         call restore(); call dealloc(); return
      end if
      ! NASG-form energy-conserving equilibrium pressure
      Peq=this%get_p_eq(VFeq,Q0,qG,gammaG)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         call restore(); call dealloc(); return
      end if
      ! Adjust densities and update conservatives (masses and velocities unchanged)
      RHOL=Q0(1)/(       VFeq)
      RHOG=Q0(2)/(1.0_WP-VFeq)
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=RHOL,y=[1.0_WP])
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=RHOG,y=y)
      call dealloc()
   contains
      subroutine restore()
         VF=VF0; Q=Q0
      end subroutine restore
      subroutine dealloc()
         if (allocated(Q0)) deallocate(Q0)
         if (allocated(y))  deallocate(y)
      end subroutine dealloc
   end subroutine p_relax

   !> Mechanical + thermal relaxation (Pelanti 2022 with co-volume b)
   subroutine pT_relax(this,dt,VF,Q,Pjump)
      implicit none
      class(relax_igmix_nasg), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: Q0,y
      real(WP) :: RHOL,RHOG,CL,CG,GL,GG,TL,TG,IL,IG
      real(WP) :: PHIL,PHIG,zetaL,zetaG,Z,D,COF
      real(WP) :: xiTL,xiTG,xiTLinv,xiTGinv
      real(WP) :: VFeq,Peq,Teq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      integer  :: iVQ
      iVQ=7+this%liq%ns+this%indV-1
      ! Store Q before relax_p so get_p_eq sees the same Q0 as the monolithic form
      allocate(Q0(size(Q))); Q0=Q
      ! Step 1: mechanical (NASG override)
      call this%p_relax(dt,VF,Q,Pjump)
      ! Step 2: thermal
      if (Q(2).gt.0.0_WP) then
         Yv=Q(iVQ)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      RHOL=Q(1)/(       VF); RHOG=Q(2)/(1.0_WP-VF)
      ! Recover p from dominant phase
      if (VF.gt.0.5_WP) then
         Peq=this%liq%get_p_from_rho_e(rho=RHOL,e=Q(3)/Q(1),y=[1.0_WP])
      else
         Peq=this%gas%get_p_from_rho_e(rho=RHOG,e=Q(4)/Q(2),y=y)
      end if
      IL=Q(3)/Q(1); IG=Q(4)/Q(2)
      TL=this%liq%get_T_from_p_rho(p=Peq,rho=RHOL,y=[1.0_WP])
      TG=this%gas%get_T_from_p_rho(p=Peq,rho=RHOG,y=y)
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL,y=[1.0_WP])
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=Peq,rho=RHOL,y=[1.0_WP])
      CG=this%gas%get_c_from_p_rho(p=Peq,rho=RHOG,y=y)
      ! ODE coefficients (with co-volume b in zetaL)
      Z=(1.0_WP-VF)*GL+VF*GG
      D=VF*RHOG*CG**2+(1.0_WP-VF)*RHOL*CL**2
      PHIL=-(this%liq%gamma-1.0_WP)*this%liq%cv*RHOL**2/(Peq+this%liq%pinf)
      PHIG=-(gammaG-1.0_WP)*cvG*RHOG**2/Peq
      zetaL=RHOL*(1.0_WP-this%liq_nasg%b*RHOL)/(Peq+this%liq%pinf)
      zetaG=RHOG/Peq
      COF=GL*RHOG*CG**2-GG*RHOL*CL**2
      xiTL=-PHIL*D/(RHOL/(       VF)*Z+zetaL*COF)
      xiTG=-PHIG*D/(RHOG/(1.0_WP-VF)*Z-zetaG*COF)
      xiTLinv=1.0_WP/xiTL; xiTGinv=1.0_WP/xiTG
      ! Equilibrium VF, T, p (NASG-form get_p_eq)
      VFeq=VF+Z/D*(TG-TL)/(xiTLinv+xiTGinv)
      Teq =(xiTL*TL+xiTG*TG)/(xiTL+xiTG)
      Peq =this%get_p_eq(VFeq,Q0,qG,gammaG)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y); return
      end if
      ! Clamp
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-this%liq%pinf); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,0.0_WP);          end if
      ! Update conservatives
      RHOL=Q(1)/(       VFeq); RHOG=Q(2)/(1.0_WP-VFeq)
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=RHOL,y=[1.0_WP])
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=RHOG,y=y)
      deallocate(Q0,y)
   end subroutine pT_relax

   !> NASG-form energy-conserving equilibrium pressure
   real(WP) function get_p_eq(this,VF_,Q0_,qG_,gammaG_) result(p_eq)
      implicit none
      class(relax_igmix_nasg), intent(in) :: this
      real(WP),                intent(in) :: VF_
      real(WP), dimension(:),  intent(in) :: Q0_
      real(WP),                intent(in) :: qG_,gammaG_
      real(WP) :: one_brho
      one_brho=1.0_WP-this%liq_nasg%b*Q0_(1)/VF_
      p_eq=(sum(Q0_(3:4))-Q0_(1)*this%liq%q-one_brho*VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0_(2)*qG_)/&
      &    (one_brho*VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG_-1.0_WP))
   end function get_p_eq

   !> Equilibrium T from energy conservation (NASG form: includes co-volume b correction)
   real(WP) function get_T_lvg(this,p_,Yv_,rho0,rhoA0) result(T)
      implicit none
      class(relax_igmix_nasg), intent(in) :: this
      real(WP),                intent(in) :: p_,Yv_,rho0,rhoA0
      T=(1.0_WP-Yv_-this%liq_nasg%b*(rho0*(1.0_WP-Yv_)-rhoA0))/                                                       &
      & ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)                          +  &
      &   rhoA0*((this%gas%gamma(this%indV)-1.0_WP)*this%gas%cv(this%indV)*Yv_                                     +  &
      &          (this%gas%gamma(this%indA)-1.0_WP)*this%gas%cv(this%indA)*(1.0_WP-Yv_))/p_)
   end function get_T_lvg

   !> Quadratic coefficients for the equilibrium-T equation (NASG form: includes co-volume terms)
   subroutine get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
      implicit none
      class(relax_igmix_nasg), intent(in)  :: this
      real(WP),                intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP),                intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      real(WP) :: cvV_,gammaV_,qV_
      cvV_   =this%gas%cv   (this%indV)
      gammaV_=this%gas%gamma(this%indV)
      qV_    =this%gas%q    (this%indV)
      ap=rho0*this%liq%cv*cvV_*((gammaV_-this%liq%gamma)*p_eq+this%liq%gamma*(gammaV_-1.0_WP)*this%liq%pinf)
      bp=(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq**2                                                    +&
      &  (this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)                                                          -&
      &   this%liq%gamma*this%liq%cv)+rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+&
      &   rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_))*p_eq                                    +&
      &   (gammaV_-1.0_WP)*cvV_*this%liq%pinf*(rho0*this%liq%q-rhoe0)
      dp=p_eq*(p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
      dapdp=rho0*this%liq%cv*cvV_*(gammaV_-this%liq%gamma)
      dbpdp=2.0_WP*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq                                             +&
      &     this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%gamma*this%liq%cv)                            +&
      &     rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)                          +&
      &     rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_)
      ddpdp=(2.0_WP*p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
   end subroutine get_coeffs_lv

end module relax_igmix_nasg_class
