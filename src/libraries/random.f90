!> This module defines several functions for random number generation.
!> @version 1.0
!> @author O. Desjardins and several others
module random
   implicit none
   private
   
   ! Function and subroutines
   public :: random_uniform,random_normal,random_lognormal,random_gamma
   public :: random_initialize
   
contains
   
   
   !> Initialization of the RNG: seeding is based on /dev/urandom
   !> or if not available on a naive RNG seeded with time and pid.
   !>
   !> This comes from the GFortran website.
   subroutine random_initialize
      use precision, only: I4,I8
      implicit none
      integer(kind=I4), allocatable, dimension(:) :: seed
      integer(kind=I4) :: i,n,un,istat,dt(8),pid
      integer(kind=I8) :: t
      ! Get seed size
      call random_seed(size=n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un,file="/dev/urandom",access="stream",form="unformatted",action="read",status="old",iostat=istat)
      if (istat.eq.0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is useful in
         ! case one launches multiple instances of the same program in parallel
         call system_clock(t)
         if (t.eq.0) then
            call date_and_time(values=dt)
            t=(dt(1)-1970)*365_I8*24*60*60*1000 &
            & +dt(2)*31_I8*24*60*60*1000+dt(3)*24_I8*60*60*1000 &
            & +dt(5)*60*60*1000+dt(6)*60*1000+dt(7)*1000+dt(8)
         end if
         pid=getpid() ! Note that this is not part of 2018 standard
         t=ieor(t,int(pid,kind(t)))
         do i=1,n
            seed(i)=lcg(t)
         end do
      end if
      ! Place the seed
      call random_seed(put=seed)
   contains
      !> Simple RNG, sufficient for seeding a better RNG
      function lcg(s) result(v)
         integer(kind=I8) :: s
         integer(kind=I4) :: v
         if (s.eq.0) then
            s=104729
         else
            s=mod(s,4294967296_I8)
         end if
         s=mod(s*279470273_I8,4294967291_I8)
         v=int(mod(s,int(huge(0),I8)),kind(0))
      end function lcg
   end subroutine random_initialize
   
   
   !> Sampling of a WP real from a uniform distribution in [lo,hi]
   function random_uniform(lo,hi) result(res)
      use precision, only: WP
      implicit none
      !> Lower and upper bounds of the distribution
      real(kind=WP), intent(in) :: lo,hi
      real(kind=WP) :: res
      ! Generate random_uniform in [0,1]
      call random_number(res)
      ! Modify to give correct bounds
      res=lo+(hi-lo)*res
   end function random_uniform
   
   
   !> Sampling of a WP real from a normal distribution
   !> Adapted from the following Fortran 77 code:
   !> ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM,
   !> IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
   !> VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
   !>
   !> The function random_normal() returns a normally distributed
   !> pseudo-random number with zero mean and unit variance.
   !> The algorithm uses the ratio of uniforms method of A.J. Kinderman
   !> and J.F. Monahan augmented with quadratic bounding curves.
   function random_normal(m,sd) result(res)
      use precision, only: WP
      implicit none
      !> Mean and standard deviation of the distribution
      real(kind=WP), intent(in), optional :: m,sd
      real(kind=WP) :: res
      real(kind=WP), parameter :: s =+0.449871_WP
      real(kind=WP), parameter :: t =-0.386595_WP
      real(kind=WP), parameter :: a =+0.19600_WP
      real(kind=WP), parameter :: b =+0.25472_WP
      real(kind=WP), parameter :: r1=+0.27597_WP
      real(kind=WP), parameter :: r2=+0.27846_WP
      real(kind=WP) :: u,v,x,y,q
      ! Generate P=(u,v) uniform in rectangle enclosing acceptance region
      try: do
         call random_number(u)
         call random_number(v)
         v=1.7156_WP*(v-0.5_WP)
         ! Evaluate the quadratic form
         x=u-s
         y=abs(v)-t
         q=x**2+y*(a*y-b*x)
         ! Accept P if inside inner ellipse
         if (q.lt.r1) exit try
         ! Reject P if outside outer ellipse
         if (q.gt.r2) cycle try
         ! Accept P if inside acceptance region
         if (v**2.lt.-4.0_WP*log(u)*u**2) exit try
      end do try
      ! Return ratio of P's coordinates as the normal deviate
      res=v/u
      ! Modify to give correct mean and standard deviation
      if (present(sd)) res=res*sd
      if (present(m))  res=res+m
   end function random_normal
   
   
   !> Sampling of a WP real from a log-normal distribution:
   !> If X has a lognormal distribution, then log(X) is normally distributed.
   !> Here the logarithm is the natural logarithm, that is to base e, sometimes
   !> denoted as ln.  To generate random variates from this distribution, generate
   !> a random deviate from the normal distribution with mean and variance equal
   !> to the mean and variance of the logarithms of X, then take its exponential.
   !>
   !> Relationship between the mean & variance of log(X) and the mean & variance
   !> of X, when X has a lognormal distribution.
   !>
   !> Let:
   !> @verbatim
   !> m = mean of log(X)
   !> s^2 = variance of log(X)
   !> @endverbatim
   !>
   !> Then:
   !> @verbatim
   !> mean of X     = exp(m + 0.5s^2)
   !> variance of X = (mean(X))^2.[exp(s^2) - 1]
   !> @endverbatim
   !>
   !> In the reverse direction:
   !> @verbatim
   !> variance of log(X) = log[1 + var(X)/(mean(X))^2]
   !> mean of log(X)     = log(mean(X) - 0.5var(log(X))
   !> @endverbatim
   function random_lognormal(m,sd) result(res)
      use precision, only: WP
      implicit none
      !> Mean and standard deviation of the log-normal distribution
      real(kind=WP), intent(in) :: m,sd
      real(kind=WP) :: res
      real(kind=WP) :: x,mlog,sdlog
      sdlog=sqrt(log(1.0_WP+(sd/m)**2))
      mlog =log(m)-0.5_WP*sdlog**2
      x    =random_normal(mlog,sdlog)
      res=exp(x)
   end function random_lognormal
   
   
   !> Function that generates a random gamma variate given a shape parameter s>0
   !> Adapted from the book:
   !     Dagpunar, J. 'Principles of random variate generation'
   !     Clarendon Press, Oxford, 1988. ISBN 0-19-852202-9
   function random_gamma(s) result(val)
      use messager,  only: die
      use precision, only: WP
      real(WP), intent(in) :: s
      real(WP) :: val
      ! Handle improper shape parameter
      if (s.lt.tiny(1.0_WP)) call die('[random_gamma] shape parameter s>=tiny(1.0) is required')
      ! Perform appropriate calculation based on s
      if (s.gt.1.0_WP) then
         gamma1: block
            real(WP) :: c,d,u,v,x
            d=s-1.0_WP/3.0_WP
            c=1.0_WP/sqrt(9.0_WP*d)
            try_gamma1: do
               inner_try_gamma1: do
                  x=random_normal(0.0_WP,1.0_WP)
                  v=(1.0_WP+c*x)**3
                  if (v.gt.0.0_WP) exit inner_try_gamma1
               end do inner_try_gamma1
               call random_number(u)
               if (u.lt.1.0_WP-0.0331_WP*x**4) then
                  val=d*v
                  exit try_gamma1
               else if (log(u).lt.0.5_WP*x**2+d*(1.0_WP-v+log(v))) then
                  val=d*v
                  exit try_gamma1
               end if
            end do try_gamma1
         end block gamma1
      else if (s.lt.1.0_WP) then
         gamma2: block
            real(WP) :: r,x,w,a,p,c,uf,vr,d
            a=1.0_WP-s
            p=a/(a+s*exp(-a))
            c=1.0_WP/s
            uf=p*(tiny(1.0_WP)/a)**s
            vr=1.0_WP-tiny(1.0_WP)
            d=a*log(a)
            try_gamma2: do
               call random_number(r)
               if (r.ge.vr) then
                  cycle try_gamma2
               else if (r.gt.p) then
                  x=a-log((1.0_WP-r)/(1.0_WP-p))
                  w=a*log(x)-d
               else if (r.gt.uf) then
                  x=a*(r/p)**c
                  w=x
               else
                  val=0.0_WP
                  return
               end if
               call random_number(r)
               if (1.0_WP-r.le.w.and.r.gt.0.0_WP) then
                  if (r*(w+1.0_WP).ge.1.0_WP) cycle try_gamma2
                  if (-log(r).le.w) cycle try_gamma2
               end if
               exit try_gamma2
            end do try_gamma2
            val=x
         end block gamma2
      else
         exponential: block
            try_exponential: do 
               call random_number(val)
               if (val.gt.0.0_WP) exit try_exponential
            end do try_exponential
            val=-log(val)
         end block exponential
      end if
   end function random_gamma
   
   
end module random
