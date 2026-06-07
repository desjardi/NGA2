!> AMR refinement-criterion utilities
module amrtag
   use precision, only: WP
   implicit none
   private
   public :: lap_error      !< Normalized second-difference indicator
   public :: grd_error      !< Relative gradient-magnitude indicator
contains
   !> Normalized second-difference error indicator for a scalar field - returns error in [0,1)
   pure real(WP) function lap_error(uc,uxm,uxp,uym,uyp,uzm,uzp,eps) result(E)
      implicit none
      real(WP), intent(in) :: uc                !< Center value
      real(WP), intent(in) :: uxm,uxp           !< x-neighbors (i-1,i+1)
      real(WP), intent(in) :: uym,uyp           !< y-neighbors (j-1,j+1)
      real(WP), intent(in) :: uzm,uzp           !< z-neighbors (k-1,k+1)
      real(WP), intent(in) :: eps               !< Noise filter (~0.01)
      real(WP) :: num,den
      ! Numerator: sum of squared second differences (curvature)
      num=(uxp-2.0_WP*uc+uxm)**2+(uyp-2.0_WP*uc+uym)**2+(uzp-2.0_WP*uc+uzm)**2
      ! Denominator: sum of squared (first-difference magnitudes + noise floor)
      den=(abs(uxp-uc)+abs(uc-uxm)+eps*(abs(uxp)+2.0_WP*abs(uc)+abs(uxm)))**2 &
      &  +(abs(uyp-uc)+abs(uc-uym)+eps*(abs(uyp)+2.0_WP*abs(uc)+abs(uym)))**2 &
      &  +(abs(uzp-uc)+abs(uc-uzm)+eps*(abs(uzp)+2.0_WP*abs(uc)+abs(uzm)))**2
      ! Dimensionless indicator
      E=sqrt(num/max(den,tiny(1.0_WP)))
   end function lap_error
   !> Relative gradient-magnitude error indicator for a scalar field - returns error in [0,inf)
   pure real(WP) function grd_error(uc,uxm,uxp,uym,uyp,uzm,uzp,eps) result(E)
      implicit none
      real(WP), intent(in) :: uc                !< Center value
      real(WP), intent(in) :: uxm,uxp           !< x-neighbors (i-1,i+1)
      real(WP), intent(in) :: uym,uyp           !< y-neighbors (j-1,j+1)
      real(WP), intent(in) :: uzm,uzp           !< z-neighbors (k-1,k+1)
      real(WP), intent(in) :: eps               !< Floor on |u| (avoids divide-by-zero)
      real(WP) :: dx2,dy2,dz2
      ! Per-direction larger one-sided jump, squared
      dx2=max(abs(uxp-uc),abs(uc-uxm))**2
      dy2=max(abs(uyp-uc),abs(uc-uym))**2
      dz2=max(abs(uzp-uc),abs(uc-uzm))**2
      ! Quadrature combine (orientation-independent), normalized by local value
      E=sqrt(dx2+dy2+dz2)/(abs(uc)+eps)
   end function grd_error
end module amrtag
