!> Single phase compressible flow solver class:
!> Provides support for RHS calculation only
module spcomp_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   use timer_class,  only: timer
   implicit none
   private
   
   ! Expose type
   public :: spcomp
   
   !> Single phase compressible solver object definition
   type :: spcomp
      
      ! This is the config around which solver is built
      class(config), pointer :: cfg
      
      ! Solver name
      character(len=str_medium) :: name='UNNAMED_SPCOMP'
      
      ! Pointers to functions to evaluate P(RHO,E), T(RHO,P), and C(RHO,P)
      procedure(Pfunc_type), pointer, nopass :: getP=>NULL()
      procedure(Tfunc_type), pointer, nopass :: getT=>NULL()
      procedure(Cfunc_type), pointer, nopass :: getC=>NULL()
      procedure(Sfunc_type), pointer, nopass :: getS=>NULL()
      
      ! Conserved variables: 1=RHO, 2=RHO*I, 3=RHO*U, 4=RHO*V, 5=RHO*W
      integer :: nQ
      real(WP), dimension(:,:,:,:), allocatable :: Q,Qold
      
      ! Flow velocity
      real(WP), dimension(:,:,:), allocatable :: U,V,W
      
      ! Internal energy
      real(WP), dimension(:,:,:), allocatable :: I
      
      ! Pressure
      real(WP), dimension(:,:,:), allocatable :: P
      
      ! Temperature
      real(WP), dimension(:,:,:), allocatable :: T
      
      ! Speed of sound
      real(WP), dimension(:,:,:), allocatable :: C
      
      ! Viscosities and heat diffusivity
      real(WP), dimension(:,:,:), allocatable :: VISC,BETA,DIFF
      
      ! Gravitational acceleration
      real(WP), dimension(3) :: gravity=0.0_WP
      
      ! Store mesh info
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,vol
      
      ! CFL numbers
      real(WP) :: CFLc_x,CFLc_y,CFLc_z                    !< Convective CFL numbers
      real(WP) :: CFLa_x,CFLa_y,CFLa_z                    !< Acoustic CFL numbers
      real(WP) :: CFLv_x,CFLv_y,CFLv_z                    !< Viscous CFL numbers
      
      ! Monitoring quantities for conserved variables
      real(WP), dimension(:), allocatable :: Qmin,Qmax,Qint
      real(WP) :: RHOKint
      real(WP) :: RHOSint
      
      ! Monitoring quantities for primitive variables
      real(WP) :: Umax,Vmax,Wmax                          !< Velocity stats
      real(WP) :: RHOmin,RHOmax                           !< Density stats
      real(WP) :: Imin,Imax                               !< Internal energy stats
      real(WP) :: Pmin,Pmax                               !< Pressure stats
      real(WP) :: Tmin,Tmax                               !< Temperature stats
      
      ! Timer
      type(timer) :: trhs                                 !< Timer for RHS calculation
      
   contains
      procedure :: print=>spcomp_print                    !< Output solver to the screen
      procedure :: initialize                             !< Initialize the flow solver
      procedure :: finalize                               !< Finalize the flow solver
      procedure :: rhs                                    !< Compute rhs of our equations using standard fluxes
      procedure :: get_primitive                          !< Calculate primitive variables from conserved variables
      procedure :: get_viscartif                          !< Calculate artifical bulk kinematic viscosity
      procedure :: get_vreman                             !< Get kinematic eddy viscosity using Vreman's model
      procedure :: get_velocity                           !< Calculate velocity from momentum
      procedure :: get_ke                                 !< Calculate kinetic energy per unit mass from velocity
      procedure :: get_momentum                           !< Calculate momentum from velocity
      procedure :: interp_vel                             !< Calculate interpolated velocity
      procedure :: get_cfl                                !< Calculate maximum CFL
      procedure :: get_info                               !< Calculate maximum field values
   end type spcomp
   
   !> Interfaces for user-defined function
   abstract interface
      !> P=P(RHO,I)
      pure real(WP) function Pfunc_type(RHO,I)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: I
      end function Pfunc_type
      !> T=T(RHO,P)
      pure real(WP) function Tfunc_type(RHO,P)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: P
      end function Tfunc_type
      !> C=C(RHO,P)
      pure real(WP) function Cfunc_type(RHO,P)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: P
      end function Cfunc_type
      !> S=S(RHO,P)
      pure real(WP) function Sfunc_type(RHO,P)
         import :: WP
         implicit none
         real(WP), intent(in) :: RHO
         real(WP), intent(in) :: P
      end function Sfunc_type
   end interface
   
contains
   
   
   !> Initialization for compressible flow solver
   subroutine initialize(this,cfg,name)
      use messager, only: die
      implicit none
      class(spcomp) :: this
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to config object
      this%cfg=>cfg
      
      ! Check that config is uniform with at least 2 cells of overlap
      if (this%cfg%no.lt.2) call die('[spcomp initialize] spcomp solver requires at least 2 cells of overlap')
      if (.not.all([this%cfg%uniform_x,this%cfg%uniform_y,this%cfg%uniform_z])) call die('[spcomp initialize] spcomp solver requires a uniform mesh')
      
      ! Store constant cell size and its inverse, handle 2D conditions, store cell volume
      this%dx=this%cfg%dx(this%cfg%imin_); this%dxi=1.0_WP/this%dx; if (this%cfg%nx.eq.1) this%dxi=0.0_WP
      this%dy=this%cfg%dy(this%cfg%jmin_); this%dyi=1.0_WP/this%dy; if (this%cfg%ny.eq.1) this%dyi=0.0_WP
      this%dz=this%cfg%dz(this%cfg%kmin_); this%dzi=1.0_WP/this%dz; if (this%cfg%nz.eq.1) this%dzi=0.0_WP
      this%vol=this%dx*this%dy*this%dz
      
      ! Allocate and zero out conserved variables
      this%nQ=5
      allocate(this%Q   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); this%Q   =0.0_WP
      allocate(this%Qold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); this%Qold=0.0_WP
      
      ! Conserved variables monitoring
      allocate(this%Qmin(1:this%nQ),this%Qmax(1:this%nQ),this%Qint(1:this%nQ))
      
      ! Flow velocity
      allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
      allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%V=0.0_WP
      allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%W=0.0_WP
      
      ! Internal energy
      allocate(this%I(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%I=0.0_WP
      
      ! Pressure
      allocate(this%P(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%P=0.0_WP
      
      ! Temperature
      allocate(this%T(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%T=0.0_WP
      
      ! Speed of sound
      allocate(this%C(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%C=0.0_WP
      
      ! Viscosities and heat diffusivity
      allocate(this%VISC(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VISC=0.0_WP
      allocate(this%BETA(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%BETA=0.0_WP
      allocate(this%DIFF(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%DIFF=0.0_WP
      
      ! Create timers
      this%trhs=timer(comm=this%cfg%comm,name='RHS')
      
   end subroutine initialize
   
   
   !> Obtain RHS for all equations
   subroutine rhs(this,dQdt)
      implicit none
      class(spcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: dQdt  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nVAR)
      real(WP), dimension(:,:,:,:), allocatable :: FQx,FQy,FQz
      integer :: i,j,k,n
      real(WP) :: w,div
      real(WP), parameter :: eps=1.0e-15_WP
      real(WP), dimension(-2: 0) :: wenop
      real(WP), dimension(-1:+1) :: wenom
      
      ! Start rhs timer
      call this%trhs%start()
      
      ! Zero out RHS
      dQdt=0.0_WP
      
      ! ================================================================ !
      ! ======================== INVISID FLUXES ======================== !
      ! ================================================================ !
       
      ! Allocate fluxes of conserved variables
      allocate(FQx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); FQx=0.0_WP
      allocate(FQy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); FQy=0.0_WP
      allocate(FQz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%nQ)); FQz=0.0_WP
      
      ! Calculate standard fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X fluxes
               ! WENO mass flux
               w=weno_weight((abs(this%Q(i-1,j,k,1)-this%Q(i-2,j,k,1))+eps)/(abs(this%Q(i,j,k,1)-this%Q(i-1,j,k,1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%Q(i+1,j,k,1)-this%Q(i  ,j,k,1))+eps)/(abs(this%Q(i,j,k,1)-this%Q(i-1,j,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               FQx(i,j,k,1)=-0.5_WP*(this%U(i,j,k)+abs(this%U(i,j,k)))*sum(wenop*this%Q(i-2:i  ,j,k,1))&
               &            -0.5_WP*(this%U(i,j,k)-abs(this%U(i,j,k)))*sum(wenom*this%Q(i-1:i+1,j,k,1))
               ! Centered  mass flux
               !FQx(i,j,k,1)=-this%U(i,j,k)*0.5_WP*sum(this%Q(i-1:i,j,k,1))
               ! WENO internal energy flux
               w=weno_weight((abs(this%I(i-1,j,k)-this%I(i-2,j,k))+eps)/(abs(this%I(i,j,k)-this%I(i-1,j,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%I(i+1,j,k)-this%I(i  ,j,k))+eps)/(abs(this%I(i,j,k)-this%I(i-1,j,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               FQx(i,j,k,2)=0.5_WP*(FQx(i,j,k,1)-abs(-FQx(i,j,k,1)))*sum(wenop*this%I(i-2:i  ,j,k))&
               &           +0.5_WP*(FQx(i,j,k,1)+abs(-FQx(i,j,k,1)))*sum(wenom*this%I(i-1:i+1,j,k))
               ! Centered internal energy flux
               !FQx(i,j,k,2)=FQx(i,j,k,1)*0.5_WP*sum(this%I(i-1:i,j,k))
               ! Y fluxes
               ! WENO  mass flux
               w=weno_weight((abs(this%Q(i,j-1,k,1)-this%Q(i,j-2,k,1))+eps)/(abs(this%Q(i,j,k,1)-this%Q(i,j-1,k,1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%Q(i,j+1,k,1)-this%Q(i,j  ,k,1))+eps)/(abs(this%Q(i,j,k,1)-this%Q(i,j-1,k,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               FQy(i,j,k,1)=-0.5_WP*(this%V(i,j,k)+abs(this%V(i,j,k)))*sum(wenop*this%Q(i,j-2:j  ,k,1))&
               &            -0.5_WP*(this%V(i,j,k)-abs(this%V(i,j,k)))*sum(wenom*this%Q(i,j-1:j+1,k,1))
               ! Centered mass flux
               !FQy(i,j,k,1)=-this%V(i,j,k)*0.5_WP*sum(this%Q(i,j-1:j,k,1))
               ! WENO internal energy flux
               w=weno_weight((abs(this%I(i,j-1,k)-this%I(i,j-2,k))+eps)/(abs(this%I(i,j,k)-this%I(i,j-1,k))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%I(i,j+1,k)-this%I(i,j  ,k))+eps)/(abs(this%I(i,j,k)-this%I(i,j-1,k))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               FQy(i,j,k,2)=0.5_WP*(FQy(i,j,k,1)-abs(-FQy(i,j,k,1)))*sum(wenop*this%I(i,j-2:j  ,k))&
               &           +0.5_WP*(FQy(i,j,k,1)+abs(-FQy(i,j,k,1)))*sum(wenom*this%I(i,j-1:j+1,k))
               ! Centered internal energy flux
               !FQy(i,j,k,2)=FQy(i,j,k,1)*0.5_WP*sum(this%I(i,j-1:j,k)) 
               ! Z fluxes
               ! WENO mass flux
               w=weno_weight((abs(this%Q(i,j,k-1,1)-this%Q(i,j,k-2,1))+eps)/(abs(this%Q(i,j,k,1)-this%Q(i,j,k-1,1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%Q(i,j,k+1,1)-this%Q(i,j,k  ,1))+eps)/(abs(this%Q(i,j,k,1)-this%Q(i,j,k-1,1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               FQz(i,j,k,1)=-0.5_WP*(this%W(i,j,k)+abs(this%W(i,j,k)))*sum(wenop*this%Q(i,j,k-2:k  ,1))&
               &            -0.5_WP*(this%W(i,j,k)-abs(this%W(i,j,k)))*sum(wenom*this%Q(i,j,k-1:k+1,1))
               ! Centered mass flux
               !FQz(i,j,k,1)=-this%W(i,j,k)*0.5_WP*sum(this%Q(i,j,k-1:k,1))
               ! WENO internal energy flux
               w=weno_weight((abs(this%I(i,j,k-1)-this%I(i,j,k-2))+eps)/(abs(this%I(i,j,k)-this%I(i,j,k-1))+eps)); wenop=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
               w=weno_weight((abs(this%I(i,j,k+1)-this%I(i,j,k  ))+eps)/(abs(this%I(i,j,k)-this%I(i,j,k-1))+eps)); wenom=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
               FQz(i,j,k,2)=0.5_WP*(FQz(i,j,k,1)-abs(-FQz(i,j,k,1)))*sum(wenop*this%I(i,j,k-2:k  ))&
               &           +0.5_WP*(FQz(i,j,k,1)+abs(-FQz(i,j,k,1)))*sum(wenom*this%I(i,j,k-1:k+1))
               ! Centered internal energy flux
               !FQz(i,j,k,2)=FQz(i,j,k,1)*0.5_WP*sum(this%I(i,j,k-1:k))
            end do
         end do
      end do
      
      ! Mass fluxes will be used to build momentum fluxes, they need to be extended by one cell on the left because of staggering
      call this%cfg%sync(FQx(:,:,:,1)); if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) FQx(this%cfg%imin-1,:,:,1)=FQx(this%cfg%imin,:,:,1)
      call this%cfg%sync(FQy(:,:,:,1)); if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) FQy(:,this%cfg%jmin-1,:,1)=FQy(:,this%cfg%jmin,:,1)
      call this%cfg%sync(FQz(:,:,:,1)); if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) FQz(:,:,this%cfg%kmin-1,1)=FQz(:,:,this%cfg%kmin,1)
      
      ! Calculate cell-centered momentum fluxes with extra cell on the left due to staggering
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               FQx(i,j,k,3)=0.25_WP*sum(FQx(i:i+1,j,k,1))*sum(this%U(i:i+1,j,k))-this%P(i,j,k)
               FQy(i,j,k,4)=0.25_WP*sum(FQy(i,j:j+1,k,1))*sum(this%V(i,j:j+1,k))-this%P(i,j,k)
               FQz(i,j,k,5)=0.25_WP*sum(FQz(i,j,k:k+1,1))*sum(this%W(i,j,k:k+1))-this%P(i,j,k)
            end do
         end do
      end do
      
      ! Calculate edge-centered momentum fluxes
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FQy(i,j,k,3)=0.25_WP*sum(FQy(i-1:i,j,k,1))*sum(this%U(i,j-1:j,k))
               FQz(i,j,k,3)=0.25_WP*sum(FQz(i-1:i,j,k,1))*sum(this%U(i,j,k-1:k))
               FQx(i,j,k,4)=0.25_WP*sum(FQx(i,j-1:j,k,1))*sum(this%V(i-1:i,j,k))
               FQz(i,j,k,4)=0.25_WP*sum(FQz(i,j-1:j,k,1))*sum(this%V(i,j,k-1:k))
               FQx(i,j,k,5)=0.25_WP*sum(FQx(i,j,k-1:k,1))*sum(this%W(i-1:i,j,k))
               FQy(i,j,k,5)=0.25_WP*sum(FQy(i,j,k-1:k,1))*sum(this%W(i,j-1:j,k))
            end do
         end do
      end do
      
      ! Assemble time derivative for conserved variables
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Mass and internal energy advection
               dQdt(i,j,k,1)=this%dxi*(FQx(i+1,j,k,1)-FQx(i,j,k,1))+this%dyi*(FQy(i,j+1,k,1)-FQy(i,j,k,1))+this%dzi*(FQz(i,j,k+1,1)-FQz(i,j,k,1))
               dQdt(i,j,k,2)=this%dxi*(FQx(i+1,j,k,2)-FQx(i,j,k,2))+this%dyi*(FQy(i,j+1,k,2)-FQy(i,j,k,2))+this%dzi*(FQz(i,j,k+1,2)-FQz(i,j,k,2))
               ! Momentum advection and pressure stress
               dQdt(i,j,k,3)=this%dxi*(FQx(i  ,j,k,3)-FQx(i-1,j,k,3))+this%dyi*(FQy(i,j+1,k,3)-FQy(i,j  ,k,3))+this%dzi*(FQz(i,j,k+1,3)-FQz(i,j,k  ,3))
               dQdt(i,j,k,4)=this%dxi*(FQx(i+1,j,k,4)-FQx(i  ,j,k,4))+this%dyi*(FQy(i,j  ,k,4)-FQy(i,j-1,k,4))+this%dzi*(FQz(i,j,k+1,4)-FQz(i,j,k  ,4))
               dQdt(i,j,k,5)=this%dxi*(FQx(i+1,j,k,5)-FQx(i  ,j,k,5))+this%dyi*(FQy(i,j+1,k,5)-FQy(i,j  ,k,5))+this%dzi*(FQz(i,j,k  ,5)-FQz(i,j,k-1,5))
               ! Pressure dilatation term
               dQdt(i,j,k,2)=dQdt(i,j,k,2)-this%P(i,j,k)*(this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k)))
            end do
         end do
      end do
      
      ! ================================================================ !
      ! ======================== VISCOUS  FLUXES ======================= !
      ! ================================================================ !
      
      ! Zero out fluxes
      FQx=0.0_WP; FQy=0.0_WP; FQz=0.0_WP
      
      ! Compute cell-centered momentum viscous fluxes
      do k=this%cfg%kmin_-1,this%cfg%kmax_
         do j=this%cfg%jmin_-1,this%cfg%jmax_
            do i=this%cfg%imin_-1,this%cfg%imax_
               div=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
               FQx(i,j,k,3)=2.0_WP*this%VISC(i,j,k)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div
               FQy(i,j,k,4)=2.0_WP*this%VISC(i,j,k)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div
               FQz(i,j,k,5)=2.0_WP*this%VISC(i,j,k)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+(this%BETA(i,j,k)-2.0_WP*this%VISC(i,j,k)/3.0_WP)*div
            end do
         end do
      end do
      
      ! Compute edge-centered momentum viscous fluxes and corresponding viscous heating
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FQy(i,j,k,3)=0.25_WP*sum(this%VISC(i-1:i,j-1:j,k))*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k))); FQx(i,j,k,4)=FQy(i,j,k,3)
               FQz(i,j,k,2)=FQy(i,j,k,3)*(this%dyi*(this%U(i,j,k)-this%U(i,j-1,k))+this%dxi*(this%V(i,j,k)-this%V(i-1,j,k)))
               FQz(i,j,k,4)=0.25_WP*sum(this%VISC(i,j-1:j,k-1:k))*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k))); FQy(i,j,k,5)=FQz(i,j,k,4)
               FQx(i,j,k,2)=FQz(i,j,k,4)*(this%dzi*(this%V(i,j,k)-this%V(i,j,k-1))+this%dyi*(this%W(i,j,k)-this%W(i,j-1,k)))
               FQx(i,j,k,5)=0.25_WP*sum(this%VISC(i-1:i,j,k-1:k))*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1))); FQz(i,j,k,3)=FQx(i,j,k,5)
               FQy(i,j,k,2)=FQx(i,j,k,5)*(this%dxi*(this%W(i,j,k)-this%W(i-1,j,k))+this%dzi*(this%U(i,j,k)-this%U(i,j,k-1)))
            end do
         end do
      end do
      
      ! Assemble time derivative for conserved variables
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Viscous momentum transport
               dQdt(i,j,k,3)=dQdt(i,j,k,3)+this%dxi*(FQx(i  ,j,k,3)-FQx(i-1,j,k,3))+this%dyi*(FQy(i,j+1,k,3)-FQy(i,j  ,k,3))+this%dzi*(FQz(i,j,k+1,3)-FQz(i,j,k  ,3))
               dQdt(i,j,k,4)=dQdt(i,j,k,4)+this%dxi*(FQx(i+1,j,k,4)-FQx(i  ,j,k,4))+this%dyi*(FQy(i,j  ,k,4)-FQy(i,j-1,k,4))+this%dzi*(FQz(i,j,k+1,4)-FQz(i,j,k  ,4))
               dQdt(i,j,k,5)=dQdt(i,j,k,5)+this%dxi*(FQx(i+1,j,k,5)-FQx(i  ,j,k,5))+this%dyi*(FQy(i,j+1,k,5)-FQy(i,j  ,k,5))+this%dzi*(FQz(i,j,k  ,5)-FQz(i,j,k-1,5))
               ! Viscous heating term
               dQdt(i,j,k,2)=dQdt(i,j,k,2)+FQx(i,j,k,3)*this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+FQy(i,j,k,4)*this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+FQz(i,j,k,5)*this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))+0.25_WP*sum(FQz(i:i+1,j:j+1,k,2))+0.25_WP*sum(FQx(i,j:j+1,k:k+1,2))+0.25_WP*sum(FQy(i:i+1,j,k:k+1,2))
            end do
         end do
      end do
      
      ! Deallocate flux arrays
      deallocate(FQx,FQy,FQz)
      
      ! Synchronize all dQdt fields
      do n=1,this%nQ; call this%cfg%sync(dQdt(:,:,:,n)); end do
         
      ! Stop rhs timer
      call this%trhs%stop()
      
   contains
      !> WENO switch function
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP ! Switching parameter
         real(WP), parameter :: delta=0.01_WP  ! Switching thickness
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight
   end subroutine rhs
   
   
   !> Calculate all primitive variables from updated conserved variables
   subroutine get_primitive(this)
      implicit none
      class(spcomp), intent(inout) :: this
      integer :: i,j,k
      ! Get velocity
      call this%get_velocity()
      ! Get primitive variables
      do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
         this%I(i,j,k)=this%Q(i,j,k,2)/this%Q(i,j,k,1)
         this%P(i,j,k)=this%getP(this%Q(i,j,k,1),this%I(i,j,k))
         this%C(i,j,k)=this%getC(this%Q(i,j,k,1),this%P(i,j,k))
      end do; end do; end do
      ! Get temperature
      if (associated(this%getT)) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            this%T(i,j,k)=this%getT(this%Q(i,j,k,1),this%P(i,j,k))
         end do; end do; end do
      end if
   end subroutine get_primitive
   
   
   !> Calculate velocity from momentum and density
   subroutine get_velocity(this)
      implicit none
      class(spcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate velocity as far as possible
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%U(i,j,k)=2.0_WP*this%Q(i,j,k,3)/sum(this%Q(i-1:i,j,k,1))
               this%V(i,j,k)=2.0_WP*this%Q(i,j,k,4)/sum(this%Q(i,j-1:j,k,1))
               this%W(i,j,k)=2.0_WP*this%Q(i,j,k,5)/sum(this%Q(i,j,k-1:k,1))
            end do
         end do
      end do
      ! Sync velocity
      call this%cfg%sync(this%U)
      call this%cfg%sync(this%V)
      call this%cfg%sync(this%W)
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         this%U(this%cfg%imino,:,:)=this%Q(this%cfg%imino,:,:,3)/(this%Q(this%cfg%imino,:,:,1))
         this%V(this%cfg%imino,:,:)=this%Q(this%cfg%imino,:,:,4)/(this%Q(this%cfg%imino,:,:,1))
         this%W(this%cfg%imino,:,:)=this%Q(this%cfg%imino,:,:,5)/(this%Q(this%cfg%imino,:,:,1))
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         this%U(:,this%cfg%jmino,:)=this%Q(:,this%cfg%jmino,:,3)/(this%Q(:,this%cfg%jmino,:,1))
         this%V(:,this%cfg%jmino,:)=this%Q(:,this%cfg%jmino,:,4)/(this%Q(:,this%cfg%jmino,:,1))
         this%W(:,this%cfg%jmino,:)=this%Q(:,this%cfg%jmino,:,5)/(this%Q(:,this%cfg%jmino,:,1))
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         this%U(:,:,this%cfg%kmino)=this%Q(:,:,this%cfg%kmino,3)/(this%Q(:,:,this%cfg%kmino,1))
         this%V(:,:,this%cfg%kmino)=this%Q(:,:,this%cfg%kmino,4)/(this%Q(:,:,this%cfg%kmino,1))
         this%W(:,:,this%cfg%kmino)=this%Q(:,:,this%cfg%kmino,5)/(this%Q(:,:,this%cfg%kmino,1))
      end if
   end subroutine get_velocity
   
   
   !> Calculate kinetic energy per unit mass from pre-calculated velocity
   !> Need to redo this better
   subroutine get_ke(this,KE)
      implicit none
      class(spcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: KE !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               KE(i,j,k)=0.5_WP*sum(this%U(i:i+1,j,k)**2+this%V(i,j:j+1,k)**2+this%W(i,j,k:k+1)**2)
            end do
         end do
      end do
      call this%cfg%sync(KE)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) KE(this%cfg%imaxo,:,:)=KE(this%cfg%imaxo-1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) KE(:,this%cfg%jmaxo,:)=KE(:,this%cfg%jmaxo-1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) KE(:,:,this%cfg%kmaxo)=KE(:,:,this%cfg%kmaxo-1)
   end subroutine get_ke
   
   
   !> Calculate momentum from velocity and density
   subroutine get_momentum(this)
      implicit none
      class(spcomp), intent(inout) :: this
      integer :: i,j,k
      ! Calculate momentum as far as possible
      do k=this%cfg%kmino_+1,this%cfg%kmaxo_
         do j=this%cfg%jmino_+1,this%cfg%jmaxo_
            do i=this%cfg%imino_+1,this%cfg%imaxo_
               this%Q(i,j,k,3)=0.5_WP*sum(this%Q(i-1:i,j,k,1))*this%U(i,j,k)
               this%Q(i,j,k,4)=0.5_WP*sum(this%Q(i,j-1:j,k,1))*this%V(i,j,k)
               this%Q(i,j,k,5)=0.5_WP*sum(this%Q(i,j,k-1:k,1))*this%W(i,j,k)
            end do
         end do
      end do
      ! Sync momentum
      call this%cfg%sync(this%Q(:,:,:,3))
      call this%cfg%sync(this%Q(:,:,:,4))
      call this%cfg%sync(this%Q(:,:,:,5))
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.1) then
         this%Q(this%cfg%imino,:,:,3)=this%Q(this%cfg%imino,:,:,1)*this%U(this%cfg%imino,:,:)
         this%Q(this%cfg%imino,:,:,4)=this%Q(this%cfg%imino,:,:,1)*this%V(this%cfg%imino,:,:)
         this%Q(this%cfg%imino,:,:,5)=this%Q(this%cfg%imino,:,:,1)*this%W(this%cfg%imino,:,:)
      end if
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         this%Q(:,this%cfg%jmino,:,3)=this%Q(:,this%cfg%jmino,:,1)*this%U(:,this%cfg%jmino,:)
         this%Q(:,this%cfg%jmino,:,4)=this%Q(:,this%cfg%jmino,:,1)*this%V(:,this%cfg%jmino,:)
         this%Q(:,this%cfg%jmino,:,5)=this%Q(:,this%cfg%jmino,:,1)*this%W(:,this%cfg%jmino,:)
      end if
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         this%Q(:,:,this%cfg%kmino,3)=this%Q(:,:,this%cfg%kmino,1)*this%U(:,:,this%cfg%kmino)
         this%Q(:,:,this%cfg%kmino,4)=this%Q(:,:,this%cfg%kmino,1)*this%V(:,:,this%cfg%kmino)
         this%Q(:,:,this%cfg%kmino,5)=this%Q(:,:,this%cfg%kmino,1)*this%W(:,:,this%cfg%kmino)
      end if
   end subroutine get_momentum
   
   
   !> Interpolate velocity to cell-center, including overlap and ghosts
   subroutine interp_vel(this,Ui,Vi,Wi)
      implicit none
      class(spcomp), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Ui !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Vi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Wi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate interpolated velocity as far as possible
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               Ui(i,j,k)=0.5_WP*sum(this%U(i:i+1,j,k))
               Vi(i,j,k)=0.5_WP*sum(this%V(i,j:j+1,k))
               Wi(i,j,k)=0.5_WP*sum(this%W(i,j,k:k+1))
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
   
   
   !> Get artifical bulk kinematic viscosity
   subroutine get_viscartif(this,dt,beta)
      implicit none
      class(spcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: beta
      integer :: i,j,k,si,sj,sk,n
      integer, parameter :: nfilter=1
      real(WP) :: max_beta,dudy,dudz,dvdx,dvdz,dwdx,dwdy,vort,grad_div
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP), parameter :: Cartif=2.0_WP
      real(WP), parameter :: Cartif_vort=100.0_WP
      real(WP), dimension(:,:,:), allocatable :: div
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]
      ! Calculate max beta permissible
      max_beta=max_cfl*min(this%dx**2,this%dy**2,this%dz**2)/(4.0_WP*dt)
      ! Zero out array
      beta=0.0_WP
      ! Compute velocity divergence
      allocate(div(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1; do j=this%cfg%jmino_,this%cfg%jmaxo_-1; do i=this%cfg%imino_,this%cfg%imaxo_-1
         div(i,j,k)=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))+this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))+this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
      end do; end do; end do
      call this%cfg%sync(div)
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) div(this%cfg%imaxo,:,:)=div(this%cfg%imaxo-1,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) div(:,this%cfg%jmaxo,:)=div(:,this%cfg%jmaxo-1,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) div(:,:,this%cfg%kmaxo)=div(:,:,this%cfg%kmaxo-1)
      ! Compute artificial bulk viscosity based on gradU provided
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Only work in compression regions
         if (div(i,j,k).ge.0.0_WP) cycle
         ! Compute local vorticity
         dudy=0.25_WP*this%dyi*sum(this%U(i:i+1,j:j+1,k)-this%U(i:i+1,j-1:j,k))
         dudz=0.25_WP*this%dzi*sum(this%U(i:i+1,j,k:k+1)-this%U(i:i+1,j,k-1:k))
         dvdx=0.25_WP*this%dxi*sum(this%V(i:i+1,j:j+1,k)-this%V(i-1:i,j:j+1,k))
         dvdz=0.25_WP*this%dzi*sum(this%V(i,j:j+1,k:k+1)-this%V(i,j:j+1,k-1:k))
         dwdx=0.25_WP*this%dxi*sum(this%W(i:i+1,j,k:k+1)-this%W(i-1:i,j,k:k+1))
         dwdy=0.25_WP*this%dyi*sum(this%W(i,j:j+1,k:k+1)-this%W(i,j-1:j,k:k+1))
         vort=(dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2
         ! Compute |grad(div)|
         grad_div=max(abs(div(i+1,j,k)-div(i,j,k)),abs(div(i,j,k)-div(i-1,j,k)))*this%dx**2&
         &       +max(abs(div(i,j+1,k)-div(i,j,k)),abs(div(i,j,k)-div(i,j-1,k)))*this%dy**2&
         &       +max(abs(div(i,j,k+1)-div(i,j,k)),abs(div(i,j,k)-div(i,j,k-1)))*this%dz**2
         ! Estimate artificial kinematic viscosity using grad(div)
         beta(i,j,k)=Cartif*grad_div*div(i,j,k)**2/(div(i,j,k)**2+Cartif_vort*vort+1.0e-15_WP)
         ! Clip it so CFL<max_CFL
         beta(i,j,k)=min(beta(i,j,k),max_beta)
      end do; end do; end do
      call this%cfg%sync(beta)
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            beta(this%cfg%imin-1,:,:)=beta(this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) beta(this%cfg%imax+1,:,:)=beta(this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            beta(:,this%cfg%jmin-1,:)=beta(:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) beta(:,this%cfg%jmax+1,:)=beta(:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            beta(:,:,this%cfg%kmin-1)=beta(:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) beta(:,:,this%cfg%kmax+1)=beta(:,:,this%cfg%kmax)
      end if
      ! Filter beta
      do n=1,nfilter
         div=beta
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            beta(i,j,k)=0.0_WP
            do sk=-1,+1; do sj=-1,+1; do si=-1,+1
               beta(i,j,k)=beta(i,j,k)+filter(si)*filter(sj)*filter(sk)*div(i+si,j+sj,k+sk)
            end do; end do; end do
         end do; end do; end do
         call this%cfg%sync(beta)
         if (.not.this%cfg%xper) then
            if (this%cfg%iproc.eq.1)            beta(this%cfg%imin-1,:,:)=beta(this%cfg%imin,:,:)
            if (this%cfg%iproc.eq.this%cfg%npx) beta(this%cfg%imax+1,:,:)=beta(this%cfg%imax,:,:)
         end if
         if (.not.this%cfg%yper) then
            if (this%cfg%jproc.eq.1)            beta(:,this%cfg%jmin-1,:)=beta(:,this%cfg%jmin,:)
            if (this%cfg%jproc.eq.this%cfg%npy) beta(:,this%cfg%jmax+1,:)=beta(:,this%cfg%jmax,:)
         end if
         if (.not.this%cfg%zper) then
            if (this%cfg%kproc.eq.1)            beta(:,:,this%cfg%kmin-1)=beta(:,:,this%cfg%kmin)
            if (this%cfg%kproc.eq.this%cfg%npz) beta(:,:,this%cfg%kmax+1)=beta(:,:,this%cfg%kmax)
         end if
      end do
      ! Free up memory
      deallocate(div)
   end subroutine get_viscartif
   
   
   !> Get kinematic eddy viscosity using Vreman's model
   subroutine get_vreman(this,dt,visc)
      implicit none
      class(spcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: visc
      real(WP), parameter :: Cs_ref=0.17_WP
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP) :: max_visc,A,B,C
      real(WP), dimension(1:3,1:3) :: beta,gradU
      real(WP), dimension(:,:,:), allocatable :: tmp
      real(WP), dimension(-1:+1), parameter :: filter=[1.0_WP/6.0_WP,2.0_WP/3.0_WP,1.0_WP/6.0_WP]
      integer :: i,j,k,si,sj,sk,n
      integer, parameter :: nfilter=1
      ! Model constant is c=2.5*Cs_ref**2 - Vreman uses c=0.07 which corresponds to Cs_ref=0.17
      C=2.5_WP*Cs_ref**2
      ! Calculate max visc permissible
      max_visc=max_cfl*min(this%dx**2,this%dy**2,this%dz**2)/(4.0_WP*dt)
      ! Zero out array
      visc=0.0_WP
      ! Compute the eddy viscosity
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         ! Compute velocity gradient tensor
         gradU(1,1)=this%dxi*(this%U(i+1,j,k)-this%U(i,j,k))
         gradU(2,1)=0.25_WP*this%dyi*sum(this%U(i:i+1,j:j+1,k)-this%U(i:i+1,j-1:j,k))
         gradU(3,1)=0.25_WP*this%dzi*sum(this%U(i:i+1,j,k:k+1)-this%U(i:i+1,j,k-1:k))
         gradU(1,2)=0.25_WP*this%dxi*sum(this%V(i:i+1,j:j+1,k)-this%V(i-1:i,j:j+1,k))
         gradU(2,2)=this%dyi*(this%V(i,j+1,k)-this%V(i,j,k))
         gradU(3,2)=0.25_WP*this%dzi*sum(this%V(i,j:j+1,k:k+1)-this%V(i,j:j+1,k-1:k))
         gradU(1,3)=0.25_WP*this%dxi*sum(this%W(i:i+1,j,k:k+1)-this%W(i-1:i,j,k:k+1))
         gradU(2,3)=0.25_WP*this%dyi*sum(this%W(i,j:j+1,k:k+1)-this%W(i,j-1:j,k:k+1))
         gradU(3,3)=this%dzi*(this%W(i,j,k+1)-this%W(i,j,k))
         ! Compute A=gradu_ij*gradu_ij invariant
         A=sum(gradU**2)
         ! Compute beta_ij=dx_m*dx_m*gradu_mi*gradu_mj
         do sj=1,3; do si=1,3; beta(si,sj)=this%dx**2*gradU(1,si)*gradU(1,sj)+this%dy**2*gradU(2,si)*gradU(2,sj)+this%dz**2*gradU(3,si)*gradU(3,sj); end do; end do
         ! Compute B invariant
         B=beta(1,1)*beta(2,2)-beta(1,2)**2+beta(1,1)*beta(3,3)-beta(1,3)**2+beta(2,2)*beta(3,3)-beta(2,3)**2
         ! Assemble algebraic eddy viscosity model
         if (B.lt.1.0e-8_WP) then
            visc(i,j,k)=0.0_WP
         else
            visc(i,j,k)=C*sqrt(B/A)
         end if
         ! Clip it so CFL<max_CFL
         visc(i,j,k)=min(visc(i,j,k),max_visc)
      end do; end do; end do
      ! Synchronize visc
      call this%cfg%sync(visc)
      if (.not.this%cfg%xper) then
         if (this%cfg%iproc.eq.1)            visc(this%cfg%imin-1,:,:)=visc(this%cfg%imin,:,:)
         if (this%cfg%iproc.eq.this%cfg%npx) visc(this%cfg%imax+1,:,:)=visc(this%cfg%imax,:,:)
      end if
      if (.not.this%cfg%yper) then
         if (this%cfg%jproc.eq.1)            visc(:,this%cfg%jmin-1,:)=visc(:,this%cfg%jmin,:)
         if (this%cfg%jproc.eq.this%cfg%npy) visc(:,this%cfg%jmax+1,:)=visc(:,this%cfg%jmax,:)
      end if
      if (.not.this%cfg%zper) then
         if (this%cfg%kproc.eq.1)            visc(:,:,this%cfg%kmin-1)=visc(:,:,this%cfg%kmin)
         if (this%cfg%kproc.eq.this%cfg%npz) visc(:,:,this%cfg%kmax+1)=visc(:,:,this%cfg%kmax)
      end if
      ! Filter visc
      do n=1,nfilter
         allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); tmp=visc
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            visc(i,j,k)=0.0_WP
            do sk=-1,+1; do sj=-1,+1; do si=-1,+1
               visc(i,j,k)=visc(i,j,k)+filter(si)*filter(sj)*filter(sk)*tmp(i+si,j+sj,k+sk)
            end do; end do; end do
         end do; end do; end do
         call this%cfg%sync(visc)
         if (.not.this%cfg%xper) then
            if (this%cfg%iproc.eq.1)            visc(this%cfg%imin-1,:,:)=visc(this%cfg%imin,:,:)
            if (this%cfg%iproc.eq.this%cfg%npx) visc(this%cfg%imax+1,:,:)=visc(this%cfg%imax,:,:)
         end if
         if (.not.this%cfg%yper) then
            if (this%cfg%jproc.eq.1)            visc(:,this%cfg%jmin-1,:)=visc(:,this%cfg%jmin,:)
            if (this%cfg%jproc.eq.this%cfg%npy) visc(:,this%cfg%jmax+1,:)=visc(:,this%cfg%jmax,:)
         end if
         if (.not.this%cfg%zper) then
            if (this%cfg%kproc.eq.1)            visc(:,:,this%cfg%kmin-1)=visc(:,:,this%cfg%kmin)
            if (this%cfg%kproc.eq.this%cfg%npz) visc(:,:,this%cfg%kmax+1)=visc(:,:,this%cfg%kmax)
         end if
         deallocate(tmp)
      end do
   end subroutine get_vreman
   
   
   !> Calculate the CFL
   subroutine get_cfl(this,dt,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(spcomp), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cfl
      integer  :: ierr
      real(WP) :: maxvisc,maxC
      ! Compute convective+acoustic CFLs
      this%CFLc_x=maxval(abs(this%U)+abs(this%C))*dt*this%dxi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_y=maxval(abs(this%V)+abs(this%C))*dt*this%dyi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLc_z=maxval(abs(this%W)+abs(this%C))*dt*this%dzi; call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLc_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      ! Compute acoustic CFLs
      maxC=maxval(this%C); call MPI_ALLREDUCE(MPI_IN_PLACE,maxC,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      this%CFLa_x=maxC*dt*this%dxi
      this%CFLa_y=maxC*dt*this%dyi
      this%CFLa_z=maxC*dt*this%dzi
      ! Compute viscous CFLs
      maxvisc=maxval((this%VISC+this%BETA)/this%Q(:,:,:,1)); call MPI_ALLREDUCE(MPI_IN_PLACE,maxvisc,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
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
      class(spcomp), intent(inout) :: this
      integer :: n,i,j,k,ierr
      real(WP), dimension(:,:,:), allocatable :: tmp
      
      ! Compute integrals and extrema of conserved variables
      do n=1,this%nQ
         call this%cfg%integrate(this%Q(:,:,:,n),integral=this%Qint(n))
      end do
      this%Qmin=+huge(1.0_WP)
      this%Qmax=-huge(1.0_WP)
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         do n=1,this%nQ
            this%Qmin(n)=min(this%Qmin(n),this%Q(i,j,k,n))
            this%Qmax(n)=max(this%Qmax(n),this%Q(i,j,k,n))
         end do
      end do; end do; end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Qmin,this%nQ,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Qmax,this%nQ,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
      ! Also compute integral of KE and entropy
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      call this%get_ke(tmp); tmp=this%Q(:,:,:,1)*tmp; call this%cfg%integrate(tmp,integral=this%RHOKint)
      this%RHOSint=0.0_WP
      if (associated(this%getS)) then
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            tmp(i,j,k)=this%Q(i,j,k,1)*this%getS(this%Q(i,j,k,1),this%P(i,j,k))
         end do; end do; end do
         call this%cfg%integrate(tmp,integral=this%RHOSint)
      end if
      deallocate(tmp)
      
      ! Calculate extrema of primitive fields
      this%RHOmin=+huge(1.0_WP); this%RHOmax=-huge(1.0_WP)
      this%Imin  =+huge(1.0_WP); this%Imax  =-huge(1.0_WP)
      this%Pmin  =+huge(1.0_WP); this%Pmax  =-huge(1.0_WP)
      this%Tmin  =+huge(1.0_WP); this%Tmax  =-huge(1.0_WP)
      this%Umax=0.0_WP; this%Vmax=0.0_WP; this%Wmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
         this%RHOmin=min(this%RHOmin,this%Q(i,j,k,1)); this%RHOmax=max(this%RHOmax,this%Q(i,j,k,1))
         this%Imin  =min(this%Imin  ,this%I  (i,j,k)); this%Imax  =max(this%Imax  ,this%I  (i,j,k))
         this%Pmin  =min(this%Pmin  ,this%P  (i,j,k)); this%Pmax  =max(this%Pmax  ,this%P  (i,j,k))
         this%Tmin  =min(this%Tmin  ,this%T  (i,j,k)); this%Tmax  =max(this%Tmax  ,this%T  (i,j,k))
         this%Umax=max(this%Umax,abs(this%U(i,j,k)))
         this%Vmax=max(this%Vmax,abs(this%V(i,j,k)))
         this%Wmax=max(this%Wmax,abs(this%W(i,j,k)))
      end do; end do; end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Imin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Imax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Pmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Pmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Tmin  ,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Tmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_info
   
   
   !> Print out info for spcomp flow solver
   subroutine spcomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(spcomp), intent(in) :: this 
      if (this%cfg%amRoot) write(output_unit,'("spcomp solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
   end subroutine spcomp_print


   !> Finalize spcomp flow solver
   subroutine finalize(this)
      implicit none
      class(spcomp), intent(inout) :: this 
      nullify(this%cfg)
      this%name='UNNAMED_SPCOMP'
      nullify(this%getP)
      nullify(this%getT)
      nullify(this%getC)
      nullify(this%getS)
      this%nQ=0
      if (allocated(this%Q))    deallocate(this%Q)
      if (allocated(this%Qold)) deallocate(this%Qold)
      if (allocated(this%U))    deallocate(this%U)
      if (allocated(this%V))    deallocate(this%V)
      if (allocated(this%W))    deallocate(this%W)
      if (allocated(this%I))    deallocate(this%I)
      if (allocated(this%P))    deallocate(this%P)
      if (allocated(this%T))    deallocate(this%T)
      if (allocated(this%C))    deallocate(this%C)
      if (allocated(this%VISC)) deallocate(this%VISC)
      if (allocated(this%BETA)) deallocate(this%BETA)
      if (allocated(this%DIFF)) deallocate(this%DIFF)
      if (allocated(this%Qmin)) deallocate(this%Qmin)
      if (allocated(this%Qmax)) deallocate(this%Qmax)
      if (allocated(this%Qint)) deallocate(this%Qint)
      call this%trhs%finalize()
   end subroutine finalize
   
   
end module spcomp_class