!> This module defines various parameters related to the precision of calculations.
!> We use double precision here!
!> @version 1.0
!> @author O. Desjardins
module precision
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: SP = REAL32  !< Appropriate KIND for a single precision real.
  integer, parameter :: DP = REAL64  !< Appropriate KIND for a double precision real.
  integer, parameter :: QP = REAL128 !< Appropriate KIND for a quadruple precision real.
  integer , parameter :: WP = DP     !< Set work precision to double precision
  real(WP), private   :: sample_real_at_WP
  real(WP), parameter :: MAX_REAL_WP = HUGE(sample_real_at_WP)
  integer , private   :: sample_int
  integer , parameter :: MAX_INTEGER = HUGE(sample_int)
end module precision
