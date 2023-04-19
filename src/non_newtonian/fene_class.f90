!> FENE model class
!> Extends multiscalar class for the calculation of source terms
module fene_class
   use multiscalar_class, only: multiscalar
   use config_class,      only: config
   use precision,         only: WP
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: fene
   
   ! List of available FENE models
   integer, parameter, public :: fenep =0            !< FENE-P model
   integer, parameter, public :: fenecr=1            !< FENE-CR model
   
   !> Constant density fene solver object definition
   type, extends(multiscalar) :: fene
      ! Model parameters
      integer  :: model                                    !< Closure model of FENE
      real(WP) :: lambda                                   !< Polymer relaxation timescale
      real(WP) :: visc                                     !< Polymer viscosity
      real(WP) :: Lmax                                     !< Polymer maximum extensibility
   contains
      procedure :: addsrc_CgradU                           !< Add C.gradU source term to residual
      procedure :: addsrc_relax                            !< Calculate FENE relaxation term
   end type fene
   
   !> Declare fene model constructor
   interface fene
      procedure construct_fene_from_args
   end interface fene
   
   
contains
   
   
   !> FENE model constructor from multiscalar
   function construct_fene_from_args(cfg,model,scheme,name) result(self)
      implicit none
      type(fene) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: model
      integer, intent(in) :: scheme
      character(len=*), optional :: name
      ! Create a six-scalar solver for conformation tensor
      self%multiscalar=multiscalar(cfg=cfg,scheme=scheme,nscalar=6,name=name)
      self%SCname(1)='Cxx'
      self%SCname(2)='Cxy'
      self%SCname(3)='Cxz'
      self%SCname(4)='Cyy'
      self%SCname(5)='Cyz'
      self%SCname(6)='Czz'
      ! Assign closure model for FENE
      self%model=model
   end function construct_fene_from_args
   

   !> Add C.gradU source terms to multiscalar residual
   subroutine addsrc_CgradU(this,gradU,resSC)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! xx tensor component
               resSC(i,j,k,1)=resSC(i,j,k,1)+2.0_WP*(this%SC(i,j,k,1)*gradU(1,1,i,j,k)+this%SC(i,j,k,2)*gradU(2,1,i,j,k)+this%SC(i,j,k,3)*gradU(3,1,i,j,k))
               ! xy tensor component
               resSC(i,j,k,2)=resSC(i,j,k,2)+(this%SC(i,j,k,2)*gradU(1,1,i,j,k)+this%SC(i,j,k,4)*gradU(2,1,i,j,k)+this%SC(i,j,k,5)*gradU(3,1,i,j,k))+&
               &                             (this%SC(i,j,k,1)*gradU(1,2,i,j,k)+this%SC(i,j,k,2)*gradU(2,2,i,j,k)+this%SC(i,j,k,3)*gradU(3,2,i,j,k))
               ! xz tensor component
               resSC(i,j,k,3)=resSC(i,j,k,3)+(this%SC(i,j,k,3)*gradU(1,1,i,j,k)+this%SC(i,j,k,5)*gradU(2,1,i,j,k)+this%SC(i,j,k,6)*gradU(3,1,i,j,k))+&
               &                             (this%SC(i,j,k,1)*gradU(1,3,i,j,k)+this%SC(i,j,k,2)*gradU(2,3,i,j,k)+this%SC(i,j,k,3)*gradU(3,3,i,j,k))
               ! yy tensor component
               resSC(i,j,k,4)=resSC(i,j,k,4)+2.0_WP*(this%SC(i,j,k,2)*gradU(1,2,i,j,k)+this%SC(i,j,k,4)*gradU(2,2,i,j,k)+this%SC(i,j,k,5)*gradU(3,2,i,j,k))
               ! yz tensor component
               resSC(i,j,k,5)=resSC(i,j,k,5)+(this%SC(i,j,k,2)*gradU(1,3,i,j,k)+this%SC(i,j,k,4)*gradU(2,3,i,j,k)+this%SC(i,j,k,5)*gradU(3,3,i,j,k))+&
               &                             (this%SC(i,j,k,3)*gradU(1,2,i,j,k)+this%SC(i,j,k,5)*gradU(2,2,i,j,k)+this%SC(i,j,k,6)*gradU(3,2,i,j,k))
               ! zz tensor component
               resSC(i,j,k,6)=resSC(i,j,k,6)+2.0_WP*(this%SC(i,j,k,3)*gradU(1,3,i,j,k)+this%SC(i,j,k,5)*gradU(2,3,i,j,k)+this%SC(i,j,k,6)*gradU(3,3,i,j,k))
            end do
         end do
      end do
   end subroutine addsrc_CgradU
   

   !> Add fene relaxation source
   subroutine addsrc_relax(this,resSC)
      use messager, only: die
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      real(WP) :: coeff
      select case (this%model)
      case (FENEP) ! Add relaxation source for FENE-P (f(r)*C-I)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle          !< Skip non-solved cells
                  coeff=(this%Lmax**2-3.00_WP)/(this%Lmax**2-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6)))
                  coeff=coeff/this%lambda                   !< Divide by relaxation time scale
                  resSC(i,j,k,1)=resSC(i,j,k,1)-coeff*this%SC(i,j,k,1)-1.0_WP !< xx tensor component
                  resSC(i,j,k,2)=resSC(i,j,k,2)-coeff*this%SC(i,j,k,2)        !< xy tensor component
                  resSC(i,j,k,3)=resSC(i,j,k,3)-coeff*this%SC(i,j,k,3)        !< xz tensor component
                  resSC(i,j,k,4)=resSC(i,j,k,4)-coeff*this%SC(i,j,k,4)-1.0_WP !< yy tensor component
                  resSC(i,j,k,5)=resSC(i,j,k,5)-coeff*this%SC(i,j,k,5)        !< yz tensor component
                  resSC(i,j,k,6)=resSC(i,j,k,6)-coeff*this%SC(i,j,k,6)-1.0_WP !< zz tensor component
               end do
            end do
         end do
      case (FENECR) ! Add relaxation source for FENE-CR (f(r)*(C-I))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle          !< Skip non-solved cells
                  coeff=this%Lmax**2/(this%Lmax**2-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6)))
                  coeff=coeff/this%lambda                   !< Divide by relaxation time scale
                  resSC(i,j,k,1)=resSC(i,j,k,1)-coeff*(this%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                  resSC(i,j,k,2)=resSC(i,j,k,2)-coeff* this%SC(i,j,k,2)         !> xy tensor component
                  resSC(i,j,k,3)=resSC(i,j,k,3)-coeff* this%SC(i,j,k,3)         !> xz tensor component
                  resSC(i,j,k,4)=resSC(i,j,k,4)-coeff*(this%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                  resSC(i,j,k,5)=resSC(i,j,k,5)-coeff* this%SC(i,j,k,5)         !> yz tensor component
                  resSC(i,j,k,6)=resSC(i,j,k,6)-coeff*(this%SC(i,j,k,6)-1.0_WP) !> zz tensor component
               end do
            end do
         end do
      case default
         call die('[FENE addsrc_relax] Unknown FENE model selected')
      end select
   end subroutine addsrc_relax
   

end module fene_class