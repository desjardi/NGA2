module mathtools
   use precision, only: WP
   implicit none
   private
   
   ! Make public what is useful outside
   public :: Pi,twoPi
   public :: fv_itp_build,fd_itp_build
   public :: inverse_matrix
   public :: cross_product
   public :: normalize
   public :: qrotate
   public :: arctan
   
   ! Trigonometric parameters
   real(WP), parameter :: Pi   =3.1415926535897932385_WP
   real(WP), parameter :: twoPi=6.2831853071795864770_WP
   
   ! Bessel first zero
   !real(WP), parameter :: bessj1_zero=3.8317059702075123115_WP
   
   ! Blasius data points
   !real(WP), dimension(0:9) :: by0=[0.000000000000000_WP,0.165571818583440_WP,0.650024518764203_WP,1.39680822972500_WP,2.30574664618049_WP,3.28327391871370_WP,4.27962110517696_WP,5.27923901129384_WP,6.27921363832835_WP,7.27921257797747_WP]
   !real(WP), dimension(0:9) :: by1=[0.000000000000000_WP,0.329780063306651_WP,0.629765721178679_WP,0.84604458266019_WP,0.95551827831671_WP,0.99154183259084_WP,0.99897290050990_WP,0.9999216098795_WP,0.99999627301467_WP,0.99999989265063_WP]
   !real(WP), dimension(0:9) :: by2=[0.332057384255589_WP,0.323007152241930_WP,0.266751564401387_WP,0.161360240845588_WP,0.06423404047594_WP,0.01590689966410_WP,0.00240199722109_WP,0.00022016340923_WP,0.00001224984692_WP,0.00000041090325_WP]
   
contains
   
   
   !> Returns cross product in 3 dimensions: z=cross(x,y)
   pure function cross_product(x,y) result(z)
      implicit none
      real(WP), dimension(3), intent(in) :: x,y
      real(WP), dimension(3) :: z
      z(1)=x(2)*y(3)-x(3)*y(2)
      z(2)=x(3)*y(1)-x(1)*y(3)
      z(3)=x(1)*y(2)-x(2)*y(1)
   end function cross_product
   
   
   !> Finite volume interpolation metric builder
   subroutine fv_itp_build(n,x,xp,coeff)
      implicit none
      integer,  intent(in) :: n                        !< Polynomial order/number of cells
      real(WP), intent(in) :: xp                       !< Point of evaluation
      real(WP), intent(in),  dimension(n+1) :: x       !< Local mesh (of size n+1)
      real(WP), intent(out), dimension(n)   :: coeff   !< Metric coefficients (of size n)
      real(WP), dimension(:,:), allocatable :: A,B
      integer :: i,j
      ! Allocate the work arrays
      allocate(A(n,n),B(n,n))
      ! Form the matrix
      do j=1,n
         do i=1,n
            A(i,j)=(x(i+1)**(n+1-j)- x(i)**(n+1-j))/(real(n+1-j,WP)*(x(i+1)-x(i)))
         end do
      end do
      ! Invert it
      call inverse_matrix(A,B,n)
      ! Compute metrics
      coeff=0.0_WP
      do j=1,n
         do i=1,n
            coeff(i)=coeff(i)+B(j,i)*xp**(n-j)
         end do
      end do
      ! Deallocate the work arrays
      deallocate(A,B)
   end subroutine fv_itp_build
   
   
   !> Finite difference interpolation metric builder
   subroutine fd_itp_build(n,x,xp,coeff)
      implicit none
      integer,  intent(in) :: n                        !< Polynomial order/number of grid points
      real(WP), intent(in) :: xp                       !< Point of evaluation
      real(WP), intent(in),  dimension(n) :: x         !< Local mesh (of size n)
      real(WP), intent(out), dimension(n) :: coeff     !< Metric coefficients (of size n)
      real(WP), dimension(:,:), allocatable :: A,B
      integer :: i,j
      ! Allocate the work arrays
      allocate(A(n,n),B(n,n))
      ! Form the matrix
      do j=1,n
         do i=1,n
            A(i,j)=(x(i)-xp)**(j-1)
         end do
      end do
      ! Invert it
      call inverse_matrix(A,B,n)
      ! Compute metrics
      coeff=B(1,:)
      ! Deallocate the work arrays
      deallocate(A,B)
   end subroutine fd_itp_build
   
   
   !> Inverse matrix using Gauss elimination
   subroutine inverse_matrix(A,B,n)
      implicit none
      integer,  intent(in) :: n                    !< Matrix size
      real(WP), intent(inout), dimension(n,n) :: A   !< Matrix to inverse - it is destroyed
      real(WP), intent(out),   dimension(n,n) :: B   !< Matrix inverse
      integer :: i,l
      ! Zero out inverse
      B=0.0_WP
      ! Forward elimination
      do i=1,n
         B(i,i)=1.0_WP
         B(i,:)=B(i,:)/A(i,i)
         A(i,:)=A(i,:)/A(i,i)
         do l=i+1,n
            B(l,:)=B(l,:)-A(l,i)*B(i,:)
            A(l,:)=A(l,:)-A(l,i)*A(i,:)
         end do
      end do
      ! Backward substitution
      do i=n,1,-1
         do l=i+1,n
            B(i,:)=B(i,:)-A(i,l)*B(l,:)
         end do
      end do
   end subroutine inverse_matrix
   
   
   ! Returns normalized vector: w=v/|v|
   pure function normalize(v) result(w)
      implicit none
      real(WP), dimension(3), intent(in) :: v
      real(WP), dimension(3)             :: w
      w=v/(norm2(v)+tiny(1.0_WP))
   end function normalize
   
   
   ! Rotates a vector v by a specified quaternion q: w=q*v*conj(q)
   pure function qrotate(v,q) result(w)
      implicit none
      real(WP), dimension(3), intent(in) :: v    !< Vector to rotate
      real(WP), dimension(4), intent(in) :: q    !< Quaternion
      real(WP), dimension(3)             :: w    !< Rotated vector
      w(1)=(2.0_WP*(q(1)*q(1)+q(2)*q(2))-1.0_WP)*v(1)+&
      &     2.0_WP*(q(2)*q(3)-q(1)*q(4))        *v(2)+&
      &     2.0_WP*(q(2)*q(4)+q(1)*q(3))        *v(3)
      w(2)= 2.0_WP*(q(2)*q(3)+q(1)*q(4))        *v(1)+&
      &    (2.0_WP*(q(1)*q(1)+q(3)*q(3))-1.0_WP)*v(2)+&
      &     2.0_WP*(q(3)*q(4)-q(1)*q(2))        *v(3)
      w(3)= 2.0_WP*(q(2)*q(4)-q(1)*q(3))        *v(1)+&
      &     2.0_WP*(q(3)*q(4)+q(1)*q(2))        *v(2)+&
      &    (2.0_WP*(q(1)*q(1)+q(4)*q(4))-1.0_WP)*v(3)
    end function qrotate

    
    ! Safe arctan
    function arctan(dx,dy)
      implicit none
      real(WP), intent(in) :: dx,dy
      real(WP) :: arctan
      if (abs(dx)+abs(dy).lt.1.0e-9_WP) then
         arctan = 0.0_WP
      else
         arctan = atan(dy/dx)
      end if
      if (dx.le.0.0_WP) then
         arctan = Pi+arctan
      elseif (dy.le.0.0_WP .and. dx.gt.0.0_WP) then
         arctan = twoPi+arctan
      end if
  end function arctan
   
   
end module mathtools
