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
   integer, parameter, public :: clipped_fenecr=2    !< FENE-CR model with clipping
   integer, parameter, public :: oldroydb=3          !< Oldroyd-B model
   
   !> Constant density fene solver object definition
   type, extends(multiscalar) :: fene
      ! Model parameters
      integer  :: model                                    !< Closure model of FENE
      real(WP) :: ncoeff                                   !< Carreau powerlaw coefficient
      real(WP) :: trelax                                   !< Polymer relaxation timescale
      real(WP) :: visc                                     !< Polymer viscosity at zero strain rate
      real(WP) :: Lmax                                     !< Polymer maximum extensibility
      ! Polymer viscosity
      real(WP), dimension(:,:,:), allocatable :: visc_p    !< Polymer viscosity
      ! Monitoring info
      real(WP) :: visc_pmin,visc_pmax                      !< Min and max polymer viscosity
   contains
      procedure :: update_visc_p                           !< Update visc_p given strain rate tensor using Carreau model
      procedure :: addsrc_CgradU                           !< Add C.gradU source term to residual
      procedure :: addsrc_relax                            !< Calculate FENE relaxation term
      procedure :: get_max=>fene_get_max                   !< Augment multiscalar's default monitoring
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
      ! Allocate and set polymer viscosity
      allocate(self%visc_p(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_p=0.0_WP
   end function construct_fene_from_args
   
   
   !> Compute visc_p from SR using Carreau model
   subroutine update_visc_p(this,SR)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR
      real(WP) :: SRmag
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Compute second invariant of strain rate tensor = sqrt(2*SR**2)
               SRmag=sqrt(2.0_WP*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2)))
               ! Compute polymer viscosity
               this%visc_p(i,j,k)=this%visc*(1.0_WP+(this%trelax*SRmag)**2)**(0.5_WP*this%ncoeff-0.5_WP)
            end do
         end do
      end do
   end subroutine update_visc_p


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
   subroutine addsrc_relax(this,resSC,dt)
      use messager, only: die
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      real(WP), intent(in) :: dt
      integer :: i,j,k
      real(WP) :: coeff
      real(WP), parameter :: safety_margin=10.0_WP
      select case (this%model)
      case (fenep) ! Add relaxation source for FENE-P (f(r)*C-I)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                            !< Skip non-solved cells
                  coeff=(this%Lmax**2-3.00_WP)/(this%Lmax**2-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6)))
                  coeff=coeff/this%trelax                                     !< Divide by relaxation time scale
                  resSC(i,j,k,1)=resSC(i,j,k,1)-coeff*this%SC(i,j,k,1)+1.0_WP !< xx tensor component
                  resSC(i,j,k,2)=resSC(i,j,k,2)-coeff*this%SC(i,j,k,2)        !< xy tensor component
                  resSC(i,j,k,3)=resSC(i,j,k,3)-coeff*this%SC(i,j,k,3)        !< xz tensor component
                  resSC(i,j,k,4)=resSC(i,j,k,4)-coeff*this%SC(i,j,k,4)+1.0_WP !< yy tensor component
                  resSC(i,j,k,5)=resSC(i,j,k,5)-coeff*this%SC(i,j,k,5)        !< yz tensor component
                  resSC(i,j,k,6)=resSC(i,j,k,6)-coeff*this%SC(i,j,k,6)+1.0_WP !< zz tensor component
               end do
            end do
         end do
      case (fenecr) ! Add relaxation source for FENE-CR (f(r)*(C-I))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
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
      case (oldroydb) ! Add relaxation source for FENE-CR ((C-I)/t_relax)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.0_WP/this%trelax                                      !< Inverse of relaxation time
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
   

   !> Calculate the min and max of our SC field
   subroutine fene_get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(fene), intent(inout) :: this
      integer :: ierr,nsc
      real(WP) :: my_SCmax,my_SCmin
      real(WP) :: my_visc_pmax,my_visc_pmin
      do nsc=1,this%nscalar
         my_SCmax=maxval(this%SC(:,:,:,nsc)); call MPI_ALLREDUCE(my_SCmax,this%SCmax(nsc),1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         my_SCmin=minval(this%SC(:,:,:,nsc)); call MPI_ALLREDUCE(my_SCmin,this%SCmin(nsc),1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      end do
      my_visc_pmax=maxval(this%visc_p); call MPI_ALLREDUCE(my_visc_pmax,this%visc_pmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_visc_pmin=minval(this%visc_p); call MPI_ALLREDUCE(my_visc_pmin,this%visc_pmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine fene_get_max


end module fene_class