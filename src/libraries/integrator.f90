module integrator
   use precision, only: WP
   implicit none
   private
   
   !> Function and subroutines
   public :: integrator_init,int_adapt,int_order
   public :: integrand_ftype
   
   !> Has the module been initialized?
   logical :: is_init=.false.
   
   !> Maximum quadrature order
   integer, parameter :: nmax=15
   
   !> Storage for quadrature rules
   real(WP), dimension(nmax,nmax) :: x,w

   !> Type of the integrand function used to perform integration
   interface
      real(WP) function integrand_ftype(x1,x2,x3)
         import :: WP
         real(WP), intent(in) :: x1,x2,x3
      end function integrand_ftype
   end interface
   
contains
   
   !> Initialization
   subroutine integrator_init()
      use mathtools, only: quadrature_rule
      implicit none
      integer :: n
      ! If initialize has been called already, return
      if (is_init) return
      ! Otherwise, prepare quadrature rules
      do n=1,nmax
         x(n,:)=0.0_WP; w(n,:)=0.0_WP
         call quadrature_rule(n,x(n,1:n),w(n,1:n))
      end do
      ! Set as initialized
      is_init=.true.
   end subroutine integrator_init
   
   !> Integrate spatially with adaptivity
   function int_adapt(func,tol) result(fint)
      implicit none
      procedure(integrand_ftype) :: func !< Function to integrate
      real(WP), intent(in) :: tol        !< Relative tolerance
      real(WP) :: fint,cint,err
      integer  :: i,j,k,n
      ! Evaluate n=1 estimate
      n=1; fint=func(0.5_WP,0.5_WP,0.5_WP); err=huge(1.0_WP)
      ! Refine until below tolerance
      do while (err.gt.tol.and.n.lt.nmax)
         ! Increment order
         n=n+1
         ! Remember current integral and calculate new estimate
         cint=fint; fint=0.0_WP
         do k=1,n; do j=1,n; do i=1,n
            fint=fint+w(n,i)*w(n,j)*w(n,k)*func(x(n,i),x(n,j),x(n,k))
         end do; end do; end do
         ! Calculate error
         err=abs(fint-cint)
      end do
   end function int_adapt
   
   !> Integrate spatially at a certain order
   function int_order(func,order) result(fint)
      implicit none
      procedure(integrand_ftype) :: func !< Function to integrate
      integer, intent(in) :: order       !< Order of integration
      real(WP) :: fint
      integer  :: i,j,k
      fint=0.0_WP
      do k=1,order; do j=1,order; do i=1,order
         fint=fint+w(order,i)*w(order,j)*w(order,k)*func(x(order,i),x(order,j),x(order,k))
      end do; end do; end do
   end function int_order
   
end module integrator