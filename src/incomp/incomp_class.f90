!> Incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density.
module incomp_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: incomp
   
   !> Incompressible solver object definition
   type :: incomp
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_INCOMP'  !< Solver name (default=UNNAMED_INCOMP)
      
      ! Constant property fluid
      real(WP) :: rho                                     !< This is our constant fluid density
      real(WP) :: visc                                    !< These is our constant fluid dynamic viscosity
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itpu_x,itpu_y,itpu_z   !< Interpolation for U
      real(WP), dimension(:,:,:,:), allocatable :: itpv_x,itpv_y,itpv_z   !< Interpolation for V
      real(WP), dimension(:,:,:,:), allocatable :: itpw_x,itpw_y,itpw_z   !< Interpolation for W
      real(WP), dimension(:,:,:,:), allocatable :: divp_x,divp_y,divp_z   !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: divu_x,divu_y,divu_z   !< Divergence for U-cell
      real(WP), dimension(:,:,:,:), allocatable :: divv_x,divv_y,divv_z   !< Divergence for V-cell
      real(WP), dimension(:,:,:,:), allocatable :: divw_x,divw_y,divw_z   !< Divergence for W-cell
      real(WP), dimension(:,:,:,:), allocatable :: grdu_x,grdu_y,grdu_z   !< Velocity gradient for U
      real(WP), dimension(:,:,:,:), allocatable :: grdv_x,grdv_y,grdv_z   !< Velocity gradient for V
      real(WP), dimension(:,:,:,:), allocatable :: grdw_x,grdw_y,grdw_z   !< Velocity gradient for W
      
      ! Boundary condition list
      !type(bcond), pointer :: bc=NULL()                   !<
      
      ! MAYBE SOME WORK ARRAYS TOO?
      
   contains
      procedure :: print=>incomp_print                    !< Output solver to the screen
      procedure :: init_metrics                           !< Initialize metrics
   end type incomp
   
   
   !> Declare incompressible solver constructor
   interface incomp
      procedure constructor
   end interface incomp
   
contains
   
   
   !> Default constructor for incompressible flow solver
   function constructor(cfg,name) result(self)
      implicit none
      type(incomp) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Prepare metrics
      call self%init_metrics()
      
      ! Allocate flow variables
      allocate(self%U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%U=0.0_WP
      allocate(self%V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%V=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%W=0.0_WP
      allocate(self%P(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P=0.0_WP
      
   end function constructor
   
   
   !> Metric initialization that accounts for VF
   subroutine init_metrics(this)
      implicit none
      class(incomp), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference velocity interpolation coefficients
      allocate(this%itpu_x( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%itpv_y( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%itpw_z( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%itpu_y(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Edge-centered (xy)
      allocate(this%itpv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Edge-centered (xy)
      allocate(this%itpv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (yz)
      allocate(this%itpw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (yz)
      allocate(this%itpw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (zx)
      allocate(this%itpu_z(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (zx)
      ! Create velocity interpolation coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1
            do i=this%cfg%imino_+1,this%cfg%imaxo_-1
               this%itpu_x(:,i,j,k)=0.5_WP*this%cfg%VF(i,j,k)*[minval(this%cfg%VF(i-1:i,j,k)),minval(this%cfg%VF(i:i+1,j,k))] !< Linear interpolation in x of U from [x ,ym,zm]
               this%itpv_y(:,i,j,k)=0.5_WP*this%cfg%VF(i,j,k)*[minval(this%cfg%VF(i,j-1:j,k)),minval(this%cfg%VF(i,j:j+1,k))] !< Linear interpolation in y of V from [xm,y ,zm]
               this%itpw_z(:,i,j,k)=0.5_WP*this%cfg%VF(i,j,k)*[minval(this%cfg%VF(i,j,k-1:k)),minval(this%cfg%VF(i,j,k:k+1))] !< Linear interpolation in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge [x ,y ,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpu_y(:,i,j,k)=minval(this%cfg%VF(i-1:i,j-1:j,k))*this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y of U from [x ,ym,zm]
               this%itpv_x(:,i,j,k)=minval(this%cfg%VF(i-1:i,j-1:j,k))*this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x of V from [xm,y ,zm]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge [xm,y ,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%itpv_z(:,i,j,k)=minval(this%cfg%VF(i,j-1:j,k-1:k))*this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z of V from [xm,y ,zm]
               this%itpw_y(:,i,j,k)=minval(this%cfg%VF(i,j-1:j,k-1:k))*this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create velocity interpolation coefficients to cell edge [x ,ym,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%itpw_x(:,i,j,k)=minval(this%cfg%VF(i-1:i,j,k-1:k))*this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x of W from [xm,ym,z ]
               this%itpu_z(:,i,j,k)=minval(this%cfg%VF(i-1:i,j,k-1:k))*this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z of U from [x ,ym,zm]
            end do
         end do
      end do
      
      ! Allocate finite volume divergence operators
      allocate(this%divp_x( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%divp_y( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%divp_z( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%divu_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Face-centered (x)
      allocate(this%divu_y( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Face-centered (x)
      allocate(this%divu_z( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Face-centered (x)
      allocate(this%divv_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Face-centered (y)
      allocate(this%divv_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Face-centered (y)
      allocate(this%divv_z( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Face-centered (y)
      allocate(this%divw_x( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Face-centered (z)
      allocate(this%divw_y( 0:+1,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Face-centered (z)
      allocate(this%divw_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Face-centered (z)
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1
            do i=this%cfg%imino_+1,this%cfg%imaxo_-1
               this%divp_x(:,i,j,k)=this%cfg%VF(i,j,k)*this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divp_y(:,i,j,k)=this%cfg%VF(i,j,k)*this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divp_z(:,i,j,k)=this%cfg%VF(i,j,k)*this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do
      ! Create divergence operator to cell face [x ,ym,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%divu_x(:,i,j,k)=minval(this%cfg%VF(i-1:i,j,k))*this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
               this%divu_y(:,i,j,k)=minval(this%cfg%VF(i-1:i,j,k))*this%cfg%dyi (j)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,y ,zm]
               this%divu_z(:,i,j,k)=minval(this%cfg%VF(i-1:i,j,k))*this%cfg%dzi (k)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,z ]
            end do
         end do
      end do
      ! Create divergence operator to cell face [xm,y ,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divv_x(:,i,j,k)=minval(this%cfg%VF(i,j-1:j,k))*this%cfg%dxi (i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,y ,zm]
               this%divv_y(:,i,j,k)=minval(this%cfg%VF(i,j-1:j,k))*this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
               this%divv_z(:,i,j,k)=minval(this%cfg%VF(i,j-1:j,k))*this%cfg%dzi (k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,z ]
            end do
         end do
      end do
      ! Create divergence operator to cell face [xm,ym,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%divw_x(:,i,j,k)=minval(this%cfg%VF(i,j,k-1:k))*this%cfg%dxi (i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,z ]
               this%divw_y(:,i,j,k)=minval(this%cfg%VF(i,j,k-1:k))*this%cfg%dyi (j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,z ]
               this%divw_z(:,i,j,k)=minval(this%cfg%VF(i,j,k-1:k))*this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,zm]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(this%grdu_x( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%grdv_y( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%grdw_z( 0:+1,this%cfg%imino_+1:this%cfg%imaxo_-1,this%cfg%jmino_+1:this%cfg%jmaxo_-1,this%cfg%kmino_+1:this%cfg%kmaxo_-1)) !< Cell-centered
      allocate(this%grdu_y(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Edge-centered (xy)
      allocate(this%grdv_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_  :this%cfg%kmaxo_  )) !< Edge-centered (xy)
      allocate(this%grdv_z(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (yz)
      allocate(this%grdw_y(-1: 0,this%cfg%imino_  :this%cfg%imaxo_  ,this%cfg%jmino_+1:this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (yz)
      allocate(this%grdw_x(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (zx)
      allocate(this%grdu_z(-1: 0,this%cfg%imino_+1:this%cfg%imaxo_  ,this%cfg%jmino_  :this%cfg%jmaxo_  ,this%cfg%kmino_+1:this%cfg%kmaxo_  )) !< Edge-centered (zx)
      ! Create gradient coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_-1
            do i=this%cfg%imino_+1,this%cfg%imaxo_-1
               this%grdu_x(:,i,j,k)=this%cfg%VF(i,j,k)*this%cfg%dxi(i)*[-minval(this%cfg%VF(i-1:i,j,k)),+minval(this%cfg%VF(i:i+1,j,k))] !< FD gradient in x of U from [x ,ym,zm]
               this%grdv_y(:,i,j,k)=this%cfg%VF(i,j,k)*this%cfg%dyi(i)*[-minval(this%cfg%VF(i,j-1:j,k)),+minval(this%cfg%VF(i,j:j+1,k))] !< FD gradient in y of V from [xm,y ,zm]
               this%grdw_z(:,i,j,k)=this%cfg%VF(i,j,k)*this%cfg%dzi(i)*[-minval(this%cfg%VF(i,j,k-1:k)),+minval(this%cfg%VF(i,j,k:k+1))] !< FD gradient in z of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge [x ,y ,zm]
      do k=this%cfg%kmino_  ,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%grdu_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y of U from [x ,ym,zm]
               this%grdv_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x of V from [xm,y ,zm]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge [xm,y ,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_  ,this%cfg%imaxo_
               this%grdv_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,1.0_WP] !< FD gradient in z of V from [xm,y ,zm]
               this%grdw_y(:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,1.0_WP] !< FD gradient in y of W from [xm,ym,z ]
            end do
         end do
      end do
      ! Create gradient coefficients to cell edge [x ,ym,z ]
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_  ,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%grdw_x(:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,1.0_WP] !< FD gradient in x of W from [xm,ym,z ]
               this%grdu_z(:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,1.0_WP] !< FD gradient in z of U from [x ,ym,zm]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Print out info for incompressible flow solver
   subroutine incomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(incomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   density = ",es12.5)') this%rho
         write(output_unit,'(" > viscosity = ",es12.5)') this%visc
      end if
      
   end subroutine incomp_print
   
   
end module incomp_class
