!> Viscoelastic model class
!> Extends scalar class for calculation of source terms
module viscoelastic_class
   use multiscalar_class, only: multiscalar
   use config_class,      only: config
   use precision,         only: WP
   implicit none
   private
   

   ! Expose type/constructor/methods
   public :: viscoelastic
   

   ! List of available viscoelastic models
   integer, parameter, public :: fenep =0            !< FENE-P model
   integer, parameter, public :: fenecr=1            !< FENE-CR model
   integer, parameter, public :: clipped_fenecr=2    !< FENE-CR model with clipping
   integer, parameter, public :: oldroydb=3          !< Oldroyd-B model
   integer, parameter, public :: lptt=4              !< Linear Phan-Thien-Tanner model
   integer, parameter, public :: eptt=5              !< Exponential Phan-Thien-Tanner model
   

   !> Constant density viscoelastic solver object definition
   type, extends(multiscalar) :: viscoelastic
      ! Model parameters
      integer  :: model                                    !< Closure model
      real(WP) :: trelax                                   !< Polymer relaxation timescale
      real(WP) :: Lmax                                     !< Polymer maximum extensibility in FENE model
      real(WP) :: affinecoeff                              !< Parameter for affine motion in PTT model
      real(WP) :: elongvisc                                !< Extensional parameter for elognational viscosity in PTT model
   contains
      procedure :: init                                    !< Viscoelastic model initialization (different name is used because of extension)
      procedure :: get_CgradU                              !< Calculate streching and distrortion term
      procedure :: get_relax                               !< Calculate relaxation term
      procedure :: get_affine                              !< Source term in PTT equation for non-affine motion
   end type viscoelastic
   
   
contains
   
   
   !> Viscoelastic model initialization
   subroutine init(this,cfg,model,scheme,name)
      implicit none
      class(viscoelastic), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: model
      integer, intent(in) :: scheme
      character(len=*), optional :: name
      ! Create a six-scalar solver for conformation tensor in the liquid
      call this%multiscalar%initialize(cfg=cfg,scheme=scheme,nscalar=6,name=name)
      this%SCname=['Cxx','Cxy','Cxz','Cyy','Cyz','Czz']
      ! Assign closure model for viscoelastic fluid
      this%model=model
   end subroutine init
   
   
   !> Get CgradU source terms to add to multiscalar residual
   subroutine get_CgradU(this,gradU,resSC)
      implicit none
      class(viscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! xx tensor component
               resSC(i,j,k,1)=2.0_WP*(this%SC(i,j,k,1)*gradU(1,1,i,j,k)+this%SC(i,j,k,2)*gradU(2,1,i,j,k)+this%SC(i,j,k,3)*gradU(3,1,i,j,k))
               ! xy tensor component
               resSC(i,j,k,2)=(this%SC(i,j,k,2)*gradU(1,1,i,j,k)+this%SC(i,j,k,4)*gradU(2,1,i,j,k)+this%SC(i,j,k,5)*gradU(3,1,i,j,k))+&
               &              (this%SC(i,j,k,1)*gradU(1,2,i,j,k)+this%SC(i,j,k,2)*gradU(2,2,i,j,k)+this%SC(i,j,k,3)*gradU(3,2,i,j,k))
               ! xz tensor component
               resSC(i,j,k,3)=(this%SC(i,j,k,3)*gradU(1,1,i,j,k)+this%SC(i,j,k,5)*gradU(2,1,i,j,k)+this%SC(i,j,k,6)*gradU(3,1,i,j,k))+&
               &              (this%SC(i,j,k,1)*gradU(1,3,i,j,k)+this%SC(i,j,k,2)*gradU(2,3,i,j,k)+this%SC(i,j,k,3)*gradU(3,3,i,j,k))
               ! yy tensor component
               resSC(i,j,k,4)=2.0_WP*(this%SC(i,j,k,2)*gradU(1,2,i,j,k)+this%SC(i,j,k,4)*gradU(2,2,i,j,k)+this%SC(i,j,k,5)*gradU(3,2,i,j,k))
               ! yz tensor component
               resSC(i,j,k,5)=(this%SC(i,j,k,2)*gradU(1,3,i,j,k)+this%SC(i,j,k,4)*gradU(2,3,i,j,k)+this%SC(i,j,k,5)*gradU(3,3,i,j,k))+&
               &              (this%SC(i,j,k,3)*gradU(1,2,i,j,k)+this%SC(i,j,k,5)*gradU(2,2,i,j,k)+this%SC(i,j,k,6)*gradU(3,2,i,j,k))
               ! zz tensor component
               resSC(i,j,k,6)=2.0_WP*(this%SC(i,j,k,3)*gradU(1,3,i,j,k)+this%SC(i,j,k,5)*gradU(2,3,i,j,k)+this%SC(i,j,k,6)*gradU(3,3,i,j,k))
            end do
         end do
      end do
   end subroutine get_CgradU
   

   !> Get S*C terms for PTT equation
   subroutine get_affine(this,SR,resSC)
      implicit none
      class(viscoelastic), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! xx tensor component
               resSC(i,j,k,1)=-this%affinecoeff*2.0_WP*(SR(1,i,j,k)*this%SC(i,j,k,1)+SR(4,i,j,k)*this%SC(i,j,k,2)+SR(6,i,j,k)*this%SC(i,j,k,3))
               ! xy tensor component
               resSC(i,j,k,2)=-this%affinecoeff*((SR(4,i,j,k)*this%SC(i,j,k,1)+SR(2,i,j,k)*this%SC(i,j,k,2)+SR(5,i,j,k)*this%SC(i,j,k,3))+&
               &                                 (SR(1,i,j,k)*this%SC(i,j,k,2)+SR(4,i,j,k)*this%SC(i,j,k,4)+SR(6,i,j,k)*this%SC(i,j,k,5)))
               ! xz tensor component
               resSC(i,j,k,3)=-this%affinecoeff*((SR(6,i,j,k)*this%SC(i,j,k,1)+SR(5,i,j,k)*this%SC(i,j,k,2)+SR(3,i,j,k)*this%SC(i,j,k,3))+&
               &                                 (SR(1,i,j,k)*this%SC(i,j,k,3)+SR(4,i,j,k)*this%SC(i,j,k,5)+SR(6,i,j,k)*this%SC(i,j,k,6)))
               ! yy tensor component
               resSC(i,j,k,4)=-this%affinecoeff*2.0_WP*(SR(4,i,j,k)*this%SC(i,j,k,2)+SR(2,i,j,k)*this%SC(i,j,k,4)+SR(5,i,j,k)*this%SC(i,j,k,5))
               ! yz tensor component
               resSC(i,j,k,5)=-this%affinecoeff*((SR(6,i,j,k)*this%SC(i,j,k,2)+SR(5,i,j,k)*this%SC(i,j,k,4)+SR(3,i,j,k)*this%SC(i,j,k,5))+&
               &                                 (SR(4,i,j,k)*this%SC(i,j,k,3)+SR(2,i,j,k)*this%SC(i,j,k,5)+SR(5,i,j,k)*this%SC(i,j,k,6)))
               ! zz tensor component
               resSC(i,j,k,6)=-this%affinecoeff*2.0_WP*(SR(6,i,j,k)*this%SC(i,j,k,3)+SR(5,i,j,k)*this%SC(i,j,k,5)+SR(3,i,j,k)*this%SC(i,j,k,6))
            end do
         end do
      end do
   end subroutine get_affine
   

   !> Add viscoelastic relaxation source
   subroutine get_relax(this,resSC,dt)
      use messager, only: die
      implicit none
      class(viscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      real(WP), intent(in) :: dt
      integer :: i,j,k
      real(WP) :: coeff
      real(WP), parameter :: safety_margin=10.0_WP
      resSC=0.0_WP
      select case (this%model)
      case (fenep) ! Add relaxation source for FENE-P (1/lambda)(f(r)*C-I)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=(this%Lmax**2-3.00_WP)/(this%Lmax**2-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6)))
                  resSC(i,j,k,1)=-(coeff*this%SC(i,j,k,1)-1.0_WP)/this%trelax !< xx tensor component
                  resSC(i,j,k,2)=-(coeff*this%SC(i,j,k,2)-0.0_WP)/this%trelax !< xy tensor component
                  resSC(i,j,k,3)=-(coeff*this%SC(i,j,k,3)-0.0_WP)/this%trelax !< xz tensor component
                  resSC(i,j,k,4)=-(coeff*this%SC(i,j,k,4)-1.0_WP)/this%trelax !< yy tensor component
                  resSC(i,j,k,5)=-(coeff*this%SC(i,j,k,5)-0.0_WP)/this%trelax !< yz tensor component
                  resSC(i,j,k,6)=-(coeff*this%SC(i,j,k,6)-1.0_WP)/this%trelax !< zz tensor component
               end do
            end do
         end do
      case (fenecr) ! Add relaxation source for FENE-CR (f(r)/lambda*(C-I))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.0_WP-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))/this%Lmax**2
                  coeff=max(epsilon(1.0_WP),min(coeff,1.0_WP))*this%trelax      !< Build a safe adjusted relaxation time scale
                  coeff=max(coeff,safety_margin*dt)                             !< Further clip based on current time step size for stability
                  coeff=1.0_WP/coeff                                            !< Inverse coeff
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP) !> xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP) !> xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP) !> yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP) !> zz tensor component
               end do
            end do
         end do
      case (clipped_fenecr) ! Add relaxation source for FENE-CR (f(r)*(C-I))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  ! Clip conformation tensor to reasonable values
                  this%SC(i,j,k,1)=min(max(0.0_WP,this%SC(i,j,k,1)),this%Lmax**2)
                  this%SC(i,j,k,4)=min(max(0.0_WP,this%SC(i,j,k,4)),this%Lmax**2)
                  this%SC(i,j,k,6)=min(max(0.0_WP,this%SC(i,j,k,6)),this%Lmax**2)
                  this%SC(i,j,k,2)=min(max(-this%Lmax**2,this%SC(i,j,k,2)),this%Lmax**2)
                  this%SC(i,j,k,3)=min(max(-this%Lmax**2,this%SC(i,j,k,3)),this%Lmax**2)
                  this%SC(i,j,k,5)=min(max(-this%Lmax**2,this%SC(i,j,k,5)),this%Lmax**2)
                  ! Calculate forcing coefficient
                  coeff=1.0_WP-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))/this%Lmax**2
                  coeff=max(epsilon(1.0_WP),min(coeff,1.0_WP))*this%trelax      !< Build a safe adjusted relaxation time scale
                  coeff=max(coeff,safety_margin*dt)                             !< Further clip based on current time step size for stability
                  coeff=1.0_WP/coeff                                            !< Inverse coeff
                  resSC(i,j,k,1)=resSC(i,j,k,1)-coeff*(this%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                  resSC(i,j,k,2)=resSC(i,j,k,2)-coeff* this%SC(i,j,k,2)         !> xy tensor component
                  resSC(i,j,k,3)=resSC(i,j,k,3)-coeff* this%SC(i,j,k,3)         !> xz tensor component
                  resSC(i,j,k,4)=resSC(i,j,k,4)-coeff*(this%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                  resSC(i,j,k,5)=resSC(i,j,k,5)-coeff* this%SC(i,j,k,5)         !> yz tensor component
                  resSC(i,j,k,6)=resSC(i,j,k,6)-coeff*(this%SC(i,j,k,6)-1.0_WP) !> zz tensor component
               end do
            end do
         end do
      case (oldroydb) ! Add relaxation source for Oldroyd-B (1/t_relax)(C-I)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.0_WP/this%trelax                                      !< Inverse of relaxation time
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP) !> xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP) !> xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP) !> yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP) !> zz tensor component
               end do
            end do
         end do
      case (lptt) ! Add relaxation source term for lPTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.00_WP+(this%elongvisc/(1.0_WP-this%affinecoeff))*((this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))-3.0_WP)
                  coeff=coeff/this%trelax                   !< Divide by relaxation time scale
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP) !> xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP) !> xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP) !> yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP) !> zz tensor component
               end do
            end do
         end do
      case (eptt) ! Add relaxation source term for ePTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=exp(this%elongvisc/(1.0_WP-this%affinecoeff)*((this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))-3.0_WP))
                  coeff=coeff/this%trelax                   !< Divide by relaxation time scale
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP) !> xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP) !> xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP) !> yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP) !> zz tensor component
               end do
            end do
         end do
      case default
         call die('[tpviscoelastic get_relax] Unknown viscoelastic model selected')
      end select
   end subroutine get_relax
   
   
end module viscoelastic_class