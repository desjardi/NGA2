!> Fast compressible flow solver class:
!> Provides support for RHS calculation only
module fastcomp_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   implicit none
   private
   
   ! Expose type
   public :: fastcomp
   
   !> Fast compressible solver object definition
   type :: fastcomp
      
      ! This is the config around which solver is built
      class(config), pointer :: cfg
      
      ! Solver name
      character(len=str_medium) :: name='UNNAMED_FASTCOMP'
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP
      
      ! Conserved variables
      real(WP), dimension(:,:,:), allocatable :: RHO ,RHOold
      real(WP), dimension(:,:,:), allocatable :: RHOU,RHOUold
      real(WP), dimension(:,:,:), allocatable :: RHOV,RHOVold
      real(WP), dimension(:,:,:), allocatable :: RHOW,RHOWold
      real(WP), dimension(:,:,:), allocatable :: RHOE,RHOEold
      
      ! Primitive variables
      real(WP), dimension(:,:,:), allocatable :: U,V,W,E,P,T
      
      ! Viscosity, diffusivity, and specific heat
      real(WP), dimension(:,:,:), allocatable :: visc,beta,diff
      real(WP) :: Cv=1.0_WP
      
      ! Store mesh size
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      
      ! CFL numbers
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                    !< Convective CFL numbers
      real(WP) :: CFLa_x,CFLa_y,CFLa_z                    !< Acoustic CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                    !< Viscous CFL numbers
      
      ! Monitoring quantities
      real(WP) :: Umax,Vmax,Wmax,RHOUint,RHOVint,RHOWint  !< Velocity stats
      real(WP) :: RHOmin,RHOmax,RHOint                    !< Density stats
      real(WP) :: Emin,Emax,RHOEint                       !< Internal energy stats
      real(WP) :: Pmin,Pmax                               !< Pressure stats
      
   contains
      procedure :: print=>fastcomp_print                  !< Output solver to the screen
      procedure :: initialize                             !< Initialize the flow solver
      procedure :: get_RHS                                !< Calculate dQ/dt
      procedure :: get_velocity                           !< Calculate velocity from momentum
      procedure :: get_momentum                           !< Calculate momentum from velocity
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_info                               !< Calculate maximum field values
   end type fastcomp
   
   
contains
   
   
   !> Initialization for compressible flow solver
   subroutine initialize(this,cfg,name)
      use messager, only: die
      implicit none
      class(fastcomp) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to config object
      this%cfg=>cfg
      
      ! Check that config is uniform with at least 3 cells of overlap
      if (this%cfg%no.lt.3) call die('[fastcomp initialize] Fastcomp solver requires at least 2 cells of overlap')
      if (.not.all([this%cfg%uniform_x,this%cfg%uniform_y,this%cfg%uniform_z])) call die('[fastcomp initialize] Fastcomp solver requires a uniform mesh')
      
      ! Store constant cell size and its inverse
      this%dx=this%cfg%dx(this%cfg%imin_); this%dxi=1.0_WP/this%dx
      this%dy=this%cfg%dy(this%cfg%jmin_); this%dyi=1.0_WP/this%dy
      this%dz=this%cfg%dz(this%cfg%kmin_); this%dzi=1.0_WP/this%dz
      
      ! Allocate and zero out storage
      allocate(this%RHO (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHO =0.0_WP
      allocate(this%RHOU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOU=0.0_WP
      allocate(this%RHOV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOV=0.0_WP
      allocate(this%RHOW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOW=0.0_WP
      allocate(this%RHOE(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOE=0.0_WP
      allocate(this%RHOold (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOold =0.0_WP
      allocate(this%RHOUold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOUold=0.0_WP
      allocate(this%RHOVold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOVold=0.0_WP
      allocate(this%RHOWold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOWold=0.0_WP
      allocate(this%RHOEold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%RHOEold=0.0_WP
      allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      allocate(this%E(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%E=0.0_WP
      allocate(this%P(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%P=0.0_WP
      allocate(this%T(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%T=0.0_WP
      allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%visc=0.0_WP
      allocate(this%beta(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%beta=0.0_WP
      allocate(this%diff(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%diff=0.0_WP
      
   end subroutine initialize
   
   
   !> Calculate the explicit time derivative of conserved variables given primitive variables
   subroutine get_RHS(this,dRHOdt,dRHOUdt,dRHOVdt,dRHOWdt,dRHOEdt)
      implicit none
      class(fastcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: dRHOdt  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: dRHOUdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: dRHOVdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: dRHOWdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: dRHOEdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(:,:,:), allocatable :: FRX,FRY,FRZ
      real(WP), dimension(:,:,:), allocatable :: FUX,FUY,FUZ
      real(WP), dimension(:,:,:), allocatable :: FVX,FVY,FVZ
      real(WP), dimension(:,:,:), allocatable :: FWX,FWY,FWZ
      real(WP), dimension(:,:,:), allocatable :: FEX,FEY,FEZ
      integer :: i,j,k
      real(WP) :: w,div
      real(WP), parameter :: eps=1.0e-15_WP
      real(WP), dimension(-2: 0) :: wxp,wyp,wzp
      real(WP), dimension(-1:+1) :: wxm,wym,wzm
      
      ! Zero out arrays
      dRHOdt =0.0_WP
      dRHOUdt=0.0_WP
      dRHOVdt=0.0_WP
      dRHOWdt=0.0_WP
      dRHOEdt=0.0_WP
      
      ! Allocate flux arrays
      allocate(FRX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FRY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FRZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FUZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FVZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FWZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FEX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FEY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FEZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! ================================================================ !
      ! ======================== INVISID FLUXES ======================== !
      ! ================================================================ !

      ! Calculate face-centered mass fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_+1
         do j=this%cfg%jmin_-1,this%cfg%jmax_+1
            do i=this%cfg%imin_-1,this%cfg%imax_+1
               ! Prepare WENO scheme for mass transport
               w=weno_weight((abs(this%RHO(i-1,j,k)-this%RHO(i-2,j,k))+eps)/(abs(this%RHO(i,j,k)-this%RHO(i-1,j,k))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%RHO(i+1,j,k)-this%RHO(i  ,j,k))+eps)/(abs(this%RHO(i,j,k)-this%RHO(i-1,j,k))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               w=weno_weight((abs(this%RHO(i,j-1,k)-this%RHO(i,j-2,k))+eps)/(abs(this%RHO(i,j,k)-this%RHO(i,j-1,k))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%RHO(i,j+1,k)-this%RHO(i,j  ,k))+eps)/(abs(this%RHO(i,j,k)-this%RHO(i,j-1,k))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               w=weno_weight((abs(this%RHO(i,j,k-1)-this%RHO(i,j,k-2))+eps)/(abs(this%RHO(i,j,k)-this%RHO(i,j,k-1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%RHO(i,j,k+1)-this%RHO(i,j,k  ))+eps)/(abs(this%RHO(i,j,k)-this%RHO(i,j,k-1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               ! WENO-based mass fluxes
               FRX(i,j,k)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wxp*this%RHO(i-2:i  ,j,k)) &
               &          -0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wxm*this%RHO(i-1:i+1,j,k))
               FRY(i,j,k)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wyp*this%RHO(i,j-2:j  ,k)) &
               &          -0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wym*this%RHO(i,j-1:j+1,k))
               FRZ(i,j,k)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wzp*this%RHO(i,j,k-2:k  )) &
               &          -0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wzm*this%RHO(i,j,k-1:k+1))
               ! Centered mass fluxes
               !FRX(i,j,k)=-0.5_WP*this%U(i,j,k)*sum(this%RHO(i-1:i,j,k))
               !FRY(i,j,k)=-0.5_WP*this%V(i,j,k)*sum(this%RHO(i,j-1:j,k))
               !FRZ(i,j,k)=-0.5_WP*this%W(i,j,k)*sum(this%RHO(i,j,k-1:k))
            end do
         end do
      end do
      
      ! Calculate cell-centered momentum fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               FUX(i,j,k)=0.25_WP*sum(FRX(i:i+1,j,k))*sum(this%U(i:i+1,j,k))-this%P(i,j,k)
               FVY(i,j,k)=0.25_WP*sum(FRY(i,j:j+1,k))*sum(this%V(i,j:j+1,k))-this%P(i,j,k)
               FWZ(i,j,k)=0.25_WP*sum(FRZ(i,j,k:k+1))*sum(this%W(i,j,k:k+1))-this%P(i,j,k)
            end do
         end do
      end do
      
      ! Calculate edge-centered momentum fluxes and face-centered internal energy fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Momentum fluxes
               FUY(i,j,k)=0.25_WP*sum(FRY(i-1:i,j,k))*sum(this%U(i,j-1:j,k))
               FUZ(i,j,k)=0.25_WP*sum(FRZ(i-1:i,j,k))*sum(this%U(i,j,k-1:k))
               FVX(i,j,k)=0.25_WP*sum(FRX(i,j-1:j,k))*sum(this%V(i-1:i,j,k))
               FVZ(i,j,k)=0.25_WP*sum(FRZ(i,j-1:j,k))*sum(this%V(i,j,k-1:k))
               FWX(i,j,k)=0.25_WP*sum(FRX(i,j,k-1:k))*sum(this%W(i-1:i,j,k))
               FWY(i,j,k)=0.25_WP*sum(FRY(i,j,k-1:k))*sum(this%W(i,j-1:j,k))
               ! Prepare WENO scheme for internal energy transport
               !w=weno_weight((abs(this%E(i-1,j,k)-this%E(i-2,j,k))+eps)/(abs(this%E(i,j,k)-this%E(i-1,j,k))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               !w=weno_weight((abs(this%E(i+1,j,k)-this%E(i  ,j,k))+eps)/(abs(this%E(i,j,k)-this%E(i-1,j,k))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               !w=weno_weight((abs(this%E(i,j-1,k)-this%E(i,j-2,k))+eps)/(abs(this%E(i,j,k)-this%E(i,j-1,k))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               !w=weno_weight((abs(this%E(i,j+1,k)-this%E(i,j  ,k))+eps)/(abs(this%E(i,j,k)-this%E(i,j-1,k))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               !w=weno_weight((abs(this%E(i,j,k-1)-this%E(i,j,k-2))+eps)/(abs(this%E(i,j,k)-this%E(i,j,k-1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               !w=weno_weight((abs(this%E(i,j,k+1)-this%E(i,j,k  ))+eps)/(abs(this%E(i,j,k)-this%E(i,j,k-1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               ! WENO-based internal energy fluxes
               !FEX(i,j,k)=0.5_WP*(FRX(i,j,k)-abs(-FRX(i,j,k)))*sum(wxp*this%E(i-2:i  ,j,k)) &
               !&         +0.5_WP*(FRX(i,j,k)+abs(-FRX(i,j,k)))*sum(wxm*this%E(i-1:i+1,j,k))
               !FEY(i,j,k)=0.5_WP*(FRY(i,j,k)-abs(-FRY(i,j,k)))*sum(wyp*this%E(i,j-2:j  ,k)) &
               !&         +0.5_WP*(FRY(i,j,k)+abs(-FRY(i,j,k)))*sum(wym*this%E(i,j-1:j+1,k))
               !FEZ(i,j,k)=0.5_WP*(FRZ(i,j,k)-abs(-FRZ(i,j,k)))*sum(wzp*this%E(i,j,k-2:k  )) &
               !&         +0.5_WP*(FRZ(i,j,k)+abs(-FRZ(i,j,k)))*sum(wzm*this%E(i,j,k-1:k+1))
               ! Centered internal energy fluxes
               FEX(i,j,k)=0.5_WP*FRX(i,j,k)*sum(this%E(i-1:i,j,k))
               FEY(i,j,k)=0.5_WP*FRY(i,j,k)*sum(this%E(i,j-1:j,k))
               FEZ(i,j,k)=0.5_WP*FRZ(i,j,k)*sum(this%E(i,j,k-1:k))
            end do
         end do
      end do
      
      ! Assemble time derivative
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Divergence of fluxes
               dRHOdt (i,j,k)=this%dxi*(FRX(i+1,j,k)-FRX(i  ,j,k))+this%dyi*(FRY(i,j+1,k)-FRY(i,j  ,k))+this%dzi*(FRZ(i,j,k+1)-FRZ(i,j,k  ))
               dRHOUdt(i,j,k)=this%dxi*(FUX(i  ,j,k)-FUX(i-1,j,k))+this%dyi*(FUY(i,j+1,k)-FUY(i,j  ,k))+this%dzi*(FUZ(i,j,k+1)-FUZ(i,j,k  ))
               dRHOVdt(i,j,k)=this%dxi*(FVX(i+1,j,k)-FVX(i  ,j,k))+this%dyi*(FVY(i,j  ,k)-FVY(i,j-1,k))+this%dzi*(FVZ(i,j,k+1)-FVZ(i,j,k  ))
               dRHOWdt(i,j,k)=this%dxi*(FWX(i+1,j,k)-FWX(i  ,j,k))+this%dyi*(FWY(i,j+1,k)-FWY(i,j  ,k))+this%dzi*(FWZ(i,j,k  )-FWZ(i,j,k-1))
               dRHOEdt(i,j,k)=this%dxi*(FEX(i+1,j,k)-FEX(i  ,j,k))+this%dyi*(FEY(i,j+1,k)-FEY(i,j  ,k))+this%dzi*(FEZ(i,j,k+1)-FEZ(i,j,k  ))
               ! Pressure dilatation term
               dRHOEdt(i,j,k)=dRHOEdt(i,j,k)-this%P(i,j,k)*(this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k)))
            end do
         end do
      end do
      
      ! ================================================================ !
      ! ======================== VISCOUS FLUXES ======================== !
      ! ================================================================ !
      
      ! Compute cell-centered fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               ! Divergence of velocity
               div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
               ! Viscous flux of momentum
               FUX(i,j,k)=2.0_WP*this%visc(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+(this%beta(i,j,k)-2.0_WP*this%visc(i,j,k)/3.0_WP)*div
               FVY(i,j,k)=2.0_WP*this%visc(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+(this%beta(i,j,k)-2.0_WP*this%visc(i,j,k)/3.0_WP)*div
               FWZ(i,j,k)=2.0_WP*this%visc(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+(this%beta(i,j,k)-2.0_WP*this%visc(i,j,k)/3.0_WP)*div
            end do
         end do
      end do
      
      ! Calculate edge-centered momentum fluxes and face-centered internal energy fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Momentum fluxes (symmetric)
               FUY(i,j,k)=0.25_WP*sum(this%visc(i-1:i,j-1:j,k))*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
               FVZ(i,j,k)=0.25_WP*sum(this%visc(i,j-1:j,k-1:k))*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
               FUZ(i,j,k)=0.25_WP*sum(this%visc(i-1:i,j,k-1:k))*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
               FVX(i,j,k)=FUY(i,j,k)
               FWY(i,j,k)=FVZ(i,j,k)
               FWX(i,j,k)=FUZ(i,j,k)
               ! Internal energy fluxes
               FEX(i,j,k)=0.5_WP*sum(this%diff(i-1:i,j,k))*this%dxi*(this%T(i,j,k)-this%T(i-1,j,k))
               FEY(i,j,k)=0.5_WP*sum(this%diff(i,j-1:j,k))*this%dyi*(this%T(i,j,k)-this%T(i,j-1,k))
               FEZ(i,j,k)=0.5_WP*sum(this%diff(i,j,k-1:k))*this%dzi*(this%T(i,j,k)-this%T(i,j,k-1))
               ! Also prepare viscous heating term
               FRZ(i,j,k)=FUY(i,j,k)*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
               FRX(i,j,k)=FVZ(i,j,k)*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
               FRY(i,j,k)=FUZ(i,j,k)*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
            end do
         end do
      end do
      
      ! Increment time derivative
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Divergence of fluxes
               dRHOUdt(i,j,k)=dRHOUdt(i,j,k)+this%dxi*(FUX(i  ,j,k)-FUX(i-1,j,k))+this%dyi*(FUY(i,j+1,k)-FUY(i,j  ,k))+this%dzi*(FUZ(i,j,k+1)-FUZ(i,j,k  ))
               dRHOVdt(i,j,k)=dRHOVdt(i,j,k)+this%dxi*(FVX(i+1,j,k)-FVX(i  ,j,k))+this%dyi*(FVY(i,j  ,k)-FVY(i,j-1,k))+this%dzi*(FVZ(i,j,k+1)-FVZ(i,j,k  ))
               dRHOWdt(i,j,k)=dRHOWdt(i,j,k)+this%dxi*(FWX(i+1,j,k)-FWX(i  ,j,k))+this%dyi*(FWY(i,j+1,k)-FWY(i,j  ,k))+this%dzi*(FWZ(i,j,k  )-FWZ(i,j,k-1))
               dRHOEdt(i,j,k)=dRHOEdt(i,j,k)+this%dxi*(FEX(i+1,j,k)-FEX(i  ,j,k))+this%dyi*(FEY(i,j+1,k)-FEY(i,j  ,k))+this%dzi*(FEZ(i,j,k+1)-FEZ(i,j,k  ))
               ! Viscous heating term
               dRHOEdt(i,j,k)=dRHOEdt(i,j,k)+FUX(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+&
               &                             FVY(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+&
               &                             FWZ(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+&
               &                             0.25_WP*sum(FRZ(i:i+1,j:j+1,k))+0.25_WP*sum(FRX(i,j:j+1,k:k+1))+0.25_WP*sum(FRY(i:i+1,j,k:k+1))
            end do
         end do
      end do

      ! Synchronize
      call this%cfg%sync(dRHOdt )
      call this%cfg%sync(dRHOUdt)
      call this%cfg%sync(dRHOVdt)
      call this%cfg%sync(dRHOWdt)
      call this%cfg%sync(dRHOEdt)
      
      ! Deallocate flux arrays
      deallocate(FRX,FRY,FRZ,FUX,FUY,FUZ,FVX,FVY,FVZ,FWX,FWY,FWZ,FEX,FEY,FEZ)
      
   contains
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP ! Switching parameter
         real(WP), parameter :: delta=0.01_WP  ! Switching thickness
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight
   end subroutine get_RHS
   
   
   !> Calculate velocity from momentum
   subroutine get_velocity(this)
      implicit none
      class(fastcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate velocity as far as possible
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%U(i,j,k)=2.0_WP*this%RHOU(i,j,k)/(this%RHO(i-1,j,k)+this%RHO(i,j,k))
               this%V(i,j,k)=2.0_WP*this%RHOV(i,j,k)/(this%RHO(i,j-1,k)+this%RHO(i,j,k))
               this%W(i,j,k)=2.0_WP*this%RHOW(i,j,k)/(this%RHO(i,j,k-1)+this%RHO(i,j,k))
            end do
         end do
      end do
      ! Sync velocity
      call this%cfg%sync(this%U)
      call this%cfg%sync(this%V)
      call this%cfg%sync(this%W)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         this%U(this%cfg%imino,:,:)=this%RHOU(this%cfg%imino,:,:)/this%RHO(this%cfg%imino,:,:)
         this%V(this%cfg%imino,:,:)=this%RHOV(this%cfg%imino,:,:)/this%RHO(this%cfg%imino,:,:)
         this%W(this%cfg%imino,:,:)=this%RHOW(this%cfg%imino,:,:)/this%RHO(this%cfg%imino,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         this%U(:,this%cfg%jmino,:)=this%RHOU(:,this%cfg%jmino,:)/this%RHO(:,this%cfg%jmino,:)
         this%V(:,this%cfg%jmino,:)=this%RHOV(:,this%cfg%jmino,:)/this%RHO(:,this%cfg%jmino,:)
         this%W(:,this%cfg%jmino,:)=this%RHOW(:,this%cfg%jmino,:)/this%RHO(:,this%cfg%jmino,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         this%U(:,:,this%cfg%kmino)=this%RHOU(:,:,this%cfg%kmino)/this%RHO(:,:,this%cfg%kmino)
         this%V(:,:,this%cfg%kmino)=this%RHOV(:,:,this%cfg%kmino)/this%RHO(:,:,this%cfg%kmino)
         this%W(:,:,this%cfg%kmino)=this%RHOW(:,:,this%cfg%kmino)/this%RHO(:,:,this%cfg%kmino)
      end if
   end subroutine get_velocity
   
   
   !> Calculate momentum from velocity
   subroutine get_momentum(this)
      implicit none
      class(fastcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate momentum as far as possible
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%RHOU(i,j,k)=0.5_WP*(this%RHO(i-1,j,k)+this%RHO(i,j,k))*this%U(i,j,k)
               this%RHOV(i,j,k)=0.5_WP*(this%RHO(i,j-1,k)+this%RHO(i,j,k))*this%V(i,j,k)
               this%RHOW(i,j,k)=0.5_WP*(this%RHO(i,j,k-1)+this%RHO(i,j,k))*this%W(i,j,k)
            end do
         end do
      end do
      ! Sync momentum
      call this%cfg%sync(this%RHOU)
      call this%cfg%sync(this%RHOV)
      call this%cfg%sync(this%RHOW)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         this%RHOU(this%cfg%imino,:,:)=this%RHO(this%cfg%imino,:,:)*this%U(this%cfg%imino,:,:)
         this%RHOV(this%cfg%imino,:,:)=this%RHO(this%cfg%imino,:,:)*this%V(this%cfg%imino,:,:)
         this%RHOW(this%cfg%imino,:,:)=this%RHO(this%cfg%imino,:,:)*this%W(this%cfg%imino,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         this%RHOU(:,this%cfg%jmino,:)=this%RHO(:,this%cfg%jmino,:)*this%U(:,this%cfg%jmino,:)
         this%RHOV(:,this%cfg%jmino,:)=this%RHO(:,this%cfg%jmino,:)*this%V(:,this%cfg%jmino,:)
         this%RHOW(:,this%cfg%jmino,:)=this%RHO(:,this%cfg%jmino,:)*this%W(:,this%cfg%jmino,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         this%RHOU(:,:,this%cfg%kmino)=this%RHO(:,:,this%cfg%kmino)*this%U(:,:,this%cfg%kmino)
         this%RHOV(:,:,this%cfg%kmino)=this%RHO(:,:,this%cfg%kmino)*this%V(:,:,this%cfg%kmino)
         this%RHOW(:,:,this%cfg%kmino)=this%RHO(:,:,this%cfg%kmino)*this%W(:,:,this%cfg%kmino)
      end if
   end subroutine get_momentum
   
   
   !> Calculate the interpolated velocity, including overlap and ghosts
   subroutine interp_vel(this,Ui,Vi,Wi)
      implicit none
      class(fastcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate interpolated velocity as far as possible
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               Ui(i,j,k)=0.5_WP*(this%U(i,j,k)+this%U(i+1,j,k))
               Vi(i,j,k)=0.5_WP*(this%V(i,j,k)+this%V(i,j+1,k))
               Wi(i,j,k)=0.5_WP*(this%W(i,j,k)+this%W(i,j,k+1))
            end do
         end do
      end do
      ! Sync interpolated velocity
      call this%cfg%sync(Ui)
      call this%cfg%sync(Vi)
      call this%cfg%sync(Wi)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         Ui(this%cfg%imaxo,:,:)=this%U(this%cfg%imaxo,:,:)
         Vi(this%cfg%imaxo,:,:)=this%V(this%cfg%imaxo,:,:)
         Wi(this%cfg%imaxo,:,:)=this%W(this%cfg%imaxo,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         Ui(:,this%cfg%jmaxo,:)=this%U(:,this%cfg%jmaxo,:)
         Vi(:,this%cfg%jmaxo,:)=this%V(:,this%cfg%jmaxo,:)
         Wi(:,this%cfg%jmaxo,:)=this%W(:,this%cfg%jmaxo,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         Ui(:,:,this%cfg%kmaxo)=this%U(:,:,this%cfg%kmaxo)
         Vi(:,:,this%cfg%kmaxo)=this%V(:,:,this%cfg%kmaxo)
         Wi(:,:,this%cfg%kmaxo)=this%W(:,:,this%cfg%kmaxo)
      end if
   end subroutine interp_vel
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,C,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(fastcomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: C !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: cfl
      integer :: ierr
      real(WP) :: maxvisc,maxC
      ! Compute convective+acoustic CFLs
      this%CFLc_x=maxval(abs(this%U)+abs(C))*dt*this%dxi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_y=maxval(abs(this%V)+abs(C))*dt*this%dyi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_z=maxval(abs(this%W)+abs(C))*dt*this%dzi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      ! Compute acoustic CFLs
      maxC=maxval(C); call MPI_ALLREDUCE(MPI_IN_PLACE,maxC,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLa_x=maxC*dt*this%dxi
      this%CFLa_y=maxC*dt*this%dyi
      this%CFLa_z=maxC*dt*this%dzi
      ! Compute viscous CFLs
      maxvisc=maxval(this%visc/this%RHO); call MPI_ALLREDUCE(MPI_IN_PLACE,maxvisc,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLv_x=4.0_WP*maxvisc*dt*this%dxi**2
      this%CFLv_y=4.0_WP*maxvisc*dt*this%dyi**2
      this%CFLv_z=4.0_WP*maxvisc*dt*this%dzi**2
      ! Return the maximum overall CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,&
      &       this%CFLa_x,this%CFLa_y,this%CFLa_z,&
      &       this%CFLv_x,this%CFLv_y,this%CFLv_z)
   end subroutine get_cfl
   
   
   !> Calculate info about our fields
   subroutine get_info(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MIN,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(fastcomp), intent(inout) :: this
      integer :: i,j,k,ierr
      ! Compute integrals
      call this%cfg%integrate(this%RHO ,integral=this%RHOint )
      call this%cfg%integrate(this%RHOU,integral=this%RHOUint)
      call this%cfg%integrate(this%RHOV,integral=this%RHOVint)
      call this%cfg%integrate(this%RHOW,integral=this%RHOWint)
      call this%cfg%integrate(this%RHOE,integral=this%RHOEint)
      ! Calculate extrema
      this%RHOmin=+huge(1.0_WP)
      this%RHOmax=-huge(1.0_WP)
      this%Umax=0.0_WP
      this%Vmax=0.0_WP
      this%Wmax=0.0_WP
      this%Emin=+huge(1.0_WP)
      this%Emax=-huge(1.0_WP)
      this%Pmin=+huge(1.0_WP)
      this%Pmax=-huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%RHOmin=min(this%RHOmin,    this%RHO(i,j,k) )
               this%RHOmax=max(this%RHOmax,    this%RHO(i,j,k) )
               this%Umax  =max(this%Umax  ,abs(this%U  (i,j,k)))
               this%Vmax  =max(this%Vmax  ,abs(this%V  (i,j,k)))
               this%Wmax  =max(this%Wmax  ,abs(this%W  (i,j,k)))
               this%Emin  =min(this%Emin  ,    this%E  (i,j,k) )
               this%Emax  =max(this%Emax  ,    this%E  (i,j,k) )
               this%Pmin  =min(this%Pmin  ,    this%P  (i,j,k) )
               this%Pmax  =max(this%Pmax  ,    this%P  (i,j,k) )
            end do
         end do
      end do
      ! Get the parallel min/max
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Emin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Emax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Pmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   end subroutine get_info
   
   
   !> Print out info for fastcomp flow solver
   subroutine fastcomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(fastcomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("fastcomp solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine fastcomp_print
   
   
end module fastcomp_class
