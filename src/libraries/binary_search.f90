!> Binary search algorithms for sorted arrays
!> Provides O(log n) search and insertion position finding
module binary_search
   use precision
   implicit none
   private
   
   public :: binary_search_int, binary_search_insert_pos
   
contains
   
   
   !> Binary search to find value in sorted integer array (O(log n))
   !> @param arr Sorted integer array (assumed to be sorted in ascending order)
   !> @param val Value to search for
   !> @return Index of value if found, 0 if not found
   function binary_search_int(arr, val) result(idx)
      implicit none
      integer, dimension(:), intent(in) :: arr
      integer, intent(in) :: val
      integer :: idx
      
      integer :: left, right, mid
      
      idx = 0
      if (size(arr) .eq. 0) return
      
      left = 1
      right = size(arr)
      
      do while (left .le. right)
         mid = (left + right) / 2
         if (arr(mid) .eq. val) then
            idx = mid
            return
         else if (arr(mid) .lt. val) then
            left = mid + 1
         else
            right = mid - 1
         end if
      end do
      
   end function binary_search_int
   
   
   !> Find insertion position in sorted integer array to maintain sorted order (O(log n))
   !> @param arr Sorted integer array (assumed to be sorted in ascending order)
   !> @param val Value to find insertion position for
   !> @return Position where val should be inserted to maintain sorted order
   function binary_search_insert_pos(arr, val) result(insert_pos)
      implicit none
      integer, dimension(:), intent(in) :: arr
      integer, intent(in) :: val
      integer :: insert_pos
      
      integer :: left, right, mid
      
      if (size(arr) .eq. 0) then
         insert_pos = 1
         return
      end if
      
      ! Binary search for insertion point
      left = 1
      right = size(arr)
      
      do while (left .le. right)
         mid = (left + right) / 2
         if (arr(mid) .lt. val) then
            left = mid + 1
         else
            right = mid - 1
         end if
      end do
      
      insert_pos = left
      
   end function binary_search_insert_pos
   
   
end module binary_search


