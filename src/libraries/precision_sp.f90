!> This module defines various parameters related to the precision of calculations.
!> We use single precision here!
!> @version 1.0
!> @author O. Desjardins
module precision
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  integer, public, parameter :: SP = REAL32  !< Appropriate KIND for a single precision real.
  integer, public, parameter :: DP = REAL64  !< Appropriate KIND for a double precision real.
  integer, public, parameter :: QP = REAL128 !< Appropriate KIND for a quadruple precision real.
  integer, public, parameter :: WP = SP      !< Set work precision to single precision.
  real(WP), private :: sample_real_at_WP
  real(WP), public, parameter :: MAX_REAL_WP = HUGE(sample_real_at_WP)
  integer , private :: sample_int
  integer , public, parameter :: MAX_INTEGER = HUGE(sample_int)
end module precision
