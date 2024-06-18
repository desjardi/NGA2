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
   public :: eigensolve3
   
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
      else if (dy.le.0.0_WP .and. dx.gt.0.0_WP) then
         arctan = twoPi+arctan
      end if
   end function arctan
   
   
   !> Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3 matrix A using the QL algorithm
   !> with implicit shifts, preceded by a Householder reduction to real tridiagonal form.
   !> The function accesses only the diagonal and upper triangular parts of A
   !> A: The symmetric input matrix
   !> Q: Storage buffer for eigenvectors
   !> W: Storage buffer for eigenvalues
   subroutine eigensolve3(A,Q,W)
      use messager, only: die
      implicit none
      real(WP), dimension(3,3), intent(in)  :: A
      real(WP), dimension(3,3), intent(out) :: Q
      real(WP), dimension(3)  , intent(out) :: W
      real(WP), dimension(3) :: E(3)
      real(WP) :: G,R,P,F,B,S,C,T
      integer :: niter,l,m,i,j,k
      ! Transform A to real tridiagonal form by the Householder method
      call householder(A,Q,W,E)
      ! Calculate eigensystem of the remaining real symmetric tridiagonal matrix with the QL method
      outer: do l=1,2
         niter=0
         do i=1,50
            ! Check for convergence and exit iteration loop if off-diagonal element E(l) is zero
            cvg: do m=l,2
               G=abs(W(m))+abs(W(m+1))
               if (abs(E(m))+G.eq.G) exit cvg
            end do cvg
            if (m.eq.l) cycle outer
            niter=niter+1
            if (niter.ge.30) call die('[DSYEVQ3] Failed to converge')
            G=(W(l+1)-W(l))/(2.0_WP*E(l))
            R=sqrt(1.0_WP+G**2)
            if (G.ge.0.0_WP) then
               G=W(m)-W(l)+E(l)/(G+R)
            else
               G=W(m)-W(l)+E(l)/(G-R)
            end if
            S=1.0_WP
            C=1.0_WP
            P=0.0_WP
            do j=m-1,l,-1
               F=S*E(j)
               B=C*E(j)
               IF (abs(F).gt.abs(G)) then
                  C=G/F
                  R=sqrt(1.0_WP+C**2)
                  E(j+1)=F*R
                  S=1.0_WP/R
                  C=C*S
               else
                  S=F/G
                  R=sqrt(1.0_WP+S**2)
                  E(j+1)=G*R
                  C=1.0_WP/R
                  S=S*C
               end if
               G=W(j+1)-P
               R=(W(j)-G)*S+2.0_WP*C*B
               P=S*R
               W(j+1)=G+P
               G=C*R-B
               ! Form eigenvectors
               do k=1,3
                  T       =  Q(k,j+1)
                  Q(k,j+1)=S*Q(k,j)+C*T
                  Q(k,j)  =C*Q(k,j)-S*T
               end do
            end do
            W(l)=W(l)-P
            E(l)=G
            E(m)=0.0_WP
         end do
      end do outer
      
   contains
      
      !> Reduces a symmetric 3x3 matrix to real tridiagonal form by applying (unitary) Householder transformations:
      !>           [ D[1]  E[1]       ]
      !>   A = Q . [ E[1]  D[2]  E[2] ] . Q^T
      !>           [       E[2]  D[3] ]
      !> The function accesses only the diagonal and upper triangular parts of A
      subroutine householder(A,Q,D,E)
         implicit none
         real(WP), dimension(3,3), intent(in)  :: A
         real(WP), dimension(3,3), intent(out) :: Q
         real(WP), dimension(3)  , intent(out) :: D
         real(WP), dimension(2)  , intent(out) :: E
         real(WP), dimension(3) :: U,P
         real(WP) :: OMEGA,F,K,H,G
         integer :: i,j
         ! Initialize Q to identity
         Q=reshape([1.0_WP,0.0_WP,0.0_WP,0.0_WP,1.0_WP,0.0_WP,0.0_WP,0.0_WP,1.0_WP],shape(Q))
         ! Bring first row and column to the desired form
         H=A(1,2)**2+A(1,3)**2
         G=sqrt(H); if (A(1,2).gt.0.0_WP) G=-G
         E(1)=G
         F   =G*A(1,2)
         U(2)=A(1,2)-G
         U(3)=A(1,3)
         OMEGA=H-F
         if (OMEGA.gt.0.0_WP) then
            OMEGA=1.0_WP/OMEGA
            K=0.0_WP
            do i=2,3
               F   =A(2,i)*U(2)+A(i,3)*U(3)
               P(i)=OMEGA*F
               K   =K+U(i)*F
            end do
            K=0.5_WP*K*OMEGA**2
            do i=2,3
               P(i)=P(i)-K*U(i)
            end do
            D(1)=A(1,1)
            D(2)=A(2,2)-2.0_WP*P(2)*U(2)
            D(3)=A(3,3)-2.0_WP*P(3)*U(3)
            ! Store inverse Householder transformation in Q
            do j=2,3
               F=OMEGA*U(j)
               do i=2,3
                  Q(i,j)=Q(i,j)-F*U(i)
               end do
            end do
            ! Calculate updated A(2,3) and store it in E(2)
            E(2)=A(2,3)-P(2)*U(3)-U(2)*P(3)
         else
            do i=1,3
               D(i)=A(i,i)
            end do
            E(2)=A(2,3)
         end if
      end subroutine householder
   
   end subroutine eigensolve3
   

end module mathtools
