!> Test amrvof geometry functions - cut_tet_vol, get_plane_dist
module mod_test_geometry
   use precision, only: WP
   use messager, only: log, die
   use string, only: rtoa
   implicit none
   private
   public :: test_geometry

contains

   !> Main test routine for geometry functions
   subroutine test_geometry()
      use amrvof_geometry, only: cut_tet_vol, get_plane_dist
      real(WP) :: tol, max_err
      
      call log("=== Testing amrvof_geometry ===")
      
      tol = 1.0e-12_WP
      max_err = 0.0_WP
      
      ! Test 1: cut_tet_vol with plane completely above tet (full tet below plane)
      test_full_tet: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, vol_exact, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         ! Unit tet at origin: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane above: n=(0,0,1), d=2 (z=2, all vertices have z<2 so all liquid)
         plane = [0.0_WP, 0.0_WP, 1.0_WP, 2.0_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         vol_exact = 1.0_WP / 6.0_WP  ! Volume of unit tet
         err = abs(vol_liq - vol_exact)
         max_err = max(max_err, err)
         
         call log("  Test 1 (full tet): vol_liq = "//rtoa(vol_liq)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 1 failed: cut_tet_vol full tet")
      end block test_full_tet
      
      ! Test 2: cut_tet_vol with plane completely below tet (zero liquid volume)
      test_empty_tet: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane below: n=(0,0,1), d=-1 (z=-1, all vertices have z>-1 so all gas)
         plane = [0.0_WP, 0.0_WP, 1.0_WP, -1.0_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         err = abs(vol_liq)
         max_err = max(max_err, err)
         
         call log("  Test 2 (empty tet): vol_liq = "//rtoa(vol_liq)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 2 failed: cut_tet_vol empty tet")
      end block test_empty_tet
      
      ! Test 3: cut_tet_vol conservation (vol_liq + vol_gas = total)
      test_conservation: block
         real(WP), dimension(3,4) :: tet
         real(WP), dimension(4) :: plane
         real(WP) :: vol_liq, vol_gas, vol_total, vol_exact, err
         real(WP), dimension(3) :: bary_liq, bary_gas
         
         tet(:,1) = [0.0_WP, 0.0_WP, 0.0_WP]
         tet(:,2) = [1.0_WP, 0.0_WP, 0.0_WP]
         tet(:,3) = [0.0_WP, 1.0_WP, 0.0_WP]
         tet(:,4) = [0.0_WP, 0.0_WP, 1.0_WP]
         
         ! Plane at z=0.3
         plane = [0.0_WP, 0.0_WP, 1.0_WP, 0.3_WP]
         
         call cut_tet_vol(tet, plane, vol_liq, vol_gas, bary_liq, bary_gas)
         vol_total = vol_liq + vol_gas
         vol_exact = 1.0_WP / 6.0_WP
         err = abs(vol_total - vol_exact)
         max_err = max(max_err, err)
         
         call log("  Test 3 (conservation): vol_liq+vol_gas = "//rtoa(vol_total)//", err = "//rtoa(err))
         if (err .gt. tol) call die("Test 3 failed: volume conservation")
      end block test_conservation
      
      ! Test 4: get_plane_dist with VF = 0.5 on unit cube
      test_plane_dist: block
         real(WP), dimension(3) :: normal, lo, hi
         real(WP) :: VF, dist, dist_expected, err
         
         ! Unit cube from (0,0,0) to (1,1,1), plane perpendicular to x-axis
         lo = [0.0_WP, 0.0_WP, 0.0_WP]
         hi = [1.0_WP, 1.0_WP, 1.0_WP]
         normal = [1.0_WP, 0.0_WP, 0.0_WP]
         VF = 0.5_WP
         
         dist = get_plane_dist(normal, lo, hi, VF)
         dist_expected = 0.5_WP
         err = abs(dist - dist_expected)
         max_err = max(max_err, err)
         
         call log("  Test 4 (plane_dist VF=0.5): dist = "//rtoa(dist)//" (expected "//rtoa(dist_expected)//")")
         if (err .gt. 1.0e-6_WP) call die("Test 4 failed: get_plane_dist")
      end block test_plane_dist
      
      ! Test 5: get_plane_dist with VF = 0 and VF = 1
      test_extreme_vf: block
         real(WP), dimension(3) :: normal, lo, hi
         real(WP) :: dist0, dist1
         
         lo = [0.0_WP, 0.0_WP, 0.0_WP]
         hi = [1.0_WP, 1.0_WP, 1.0_WP]
         normal = [1.0_WP, 0.0_WP, 0.0_WP]
         
         dist0 = get_plane_dist(normal, lo, hi, 0.0_WP)
         dist1 = get_plane_dist(normal, lo, hi, 1.0_WP)
         
         call log("  Test 5 (extreme VF): dist0 = "//rtoa(dist0)//", dist1 = "//rtoa(dist1))
         ! dist0 should be <= 0 (plane at or before origin)
         ! dist1 should be >= 1 (plane at or after far end)
         if (dist0 .gt. 1.0e-10_WP) call die("Test 5 failed: VF=0 should give dist<=0")
         if (dist1 .lt. 1.0_WP - 1.0e-10_WP) call die("Test 5 failed: VF=1 should give dist>=1")
      end block test_extreme_vf
      
      call log("=== All geometry tests passed! Max error: "//rtoa(max_err)//" ===")
      
   end subroutine test_geometry

end module mod_test_geometry

