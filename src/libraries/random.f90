!> This module defines several functions for random number generation.
!> @version 1.0
!> @author O. Desjardins and several others
module random
  use precision
  implicit none
  
  !> Whether the RNG has been initialized
  logical :: rng_init=.false.
  
contains
  
  
  !> Sampling of a WP real from a uniform distribution in [lo,hi]
  real(WP) function random_uniform(lo,hi)
    implicit none
    !> Lower and upper bounds of the distribution
    real(WP), intent(in) :: lo,hi
    ! Generate random_uniform in [0,1]
    call random_number(random_uniform)
    ! Modify to give correct bounds
    random_uniform=lo+(hi-lo)*random_uniform
    return
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
  real(WP) function random_normal(m,sd)
    implicit none
    !> Mean and standard deviation of the distribution
    real(WP), intent(in), optional :: m,sd
    real(WP) :: s = 0.449871_WP
    real(WP) :: t = -0.386595_WP
    real(WP) :: a = 0.19600_WP
    real(WP) :: b = 0.25472_WP
    real(WP) :: r1 = 0.27597_WP
    real(WP) :: r2 = 0.27846_WP
    real(WP) :: u,v,x,y,q
    ! Generate P = (u,v) uniform in rectangle enclosing acceptance region
    do    
       call random_number(u)
       call random_number(v)
       v=1.7156_WP*(v-0.5_WP)
       ! Evaluate the quadratic form
       x=u-s
       y=abs(v)-t
       q=x**2+y*(a*y-b*x)
       ! Accept P if inside inner ellipse
       if (q.lt.r1) exit
       ! Reject P if outside outer ellipse
       if (q.gt.r2) cycle
       ! Reject P if outside acceptance region
       if (v**2.lt.-4.0_WP*log(u)*u**2) exit
    end do
    ! Return ratio of P's coordinates as the normal deviate
    random_normal = v/u
    ! Modify to give correct mean and standard deviation
    if (present(sd)) random_normal = random_normal*sd
    if (present(m))  random_normal = random_normal+m
    return
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
  real(WP) function random_lognormal(m,sd)
    implicit none
    !> Mean and standard deviation of the log-normal distribution
    real(WP), intent(in) :: m,sd
    real(WP) :: x,mlog,sdlog
    sdlog = sqrt(log(1.0_WP+(sd/m)**2))
    mlog  = log(m)-0.5_WP*sdlog**2
    x     = random_normal(mlog,sdlog)
    random_lognormal = exp(x)
    return
  end function random_lognormal
  
end module random


!> Initialization of the RNG: seeding is based on /dev/urandom
!> or if not available on a naive RNG seeded with time and pid.
!>
!> This comes from the GFortran website.
subroutine random_init
  use iso_fortran_env, only: int64
  use random
  implicit none
  integer, allocatable, dimension(:) :: seed
  integer :: i,n,un,istat,dt(8),pid
  integer(int64) :: t
  ! Get seed size
  call random_seed(size=n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un,file="/dev/urandom",access="stream", &
       form="unformatted",action="read",status="old",iostat=istat)
  if (istat.eq.0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is useful in
     ! case one launches multiple instances of the same program in parallel
     call system_clock(t)
     if (t.eq.0) then
        call date_and_time(values=dt)
        t=(dt(1)-1970)*365_int64*24*60*60*1000 &
        & +dt(2)*31_int64*24*60*60*1000+dt(3)*24_int64*60*60*1000 &
        & +dt(5)*60*60*1000+dt(6)*60*1000+dt(7)*1000+dt(8)
     end if
     pid=getpid()
     t=ieor(t,int(pid,kind(t)))
     do i=1,n
        seed(i)=lcg(t)
     end do
  end if
  ! Place the seed
  call random_seed(put=seed)
  ! Adjust flag
  rng_init=.true.
contains
  !> Simple RNG, sufficient for seeding a better RNG
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s.eq.0) then
       s=104729
    else
       s=mod(s,4294967296_int64)
    end if
    s=mod(s*279470273_int64,4294967291_int64)
    lcg=int(mod(s,int(huge(0),int64)),kind(0))
  end function lcg
end subroutine random_init

!> Termination of random module.
subroutine random_final
  use random
  implicit none
  ! Don't call me a dummy - you're the dummy!
  return
end subroutine random_final
