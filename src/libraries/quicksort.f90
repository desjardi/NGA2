module quicksort
   use precision
   implicit none
   

   public  :: quick_sort,quick_sort_int
   private :: qs_partition,qs_partition_int
   
   
contains
   
   
   !> Recursive quicksort algorithm : this routine recursively sorts a real
   !> array by increasing order and sorts a integer array at the same time
   recursive subroutine quick_sort(A,B)
      implicit none
      real(WP), dimension(:), intent(inout) :: A
      integer , dimension(:), intent(inout) :: B
      integer :: imark
      if (size(A).gt.1) then
         call qs_partition(A,B,imark)
         call quick_sort(A(     :imark-1),B(     :imark-1))
         call quick_sort(A(imark:       ),B(imark:       ))
      end if
   end subroutine quick_sort
   
   
   !> Recursive quicksort algorithm, single integer array
   recursive subroutine quick_sort_int(A)
      implicit none
      integer, dimension(:), intent(inout) :: A
      integer :: imark
      if (size(A).gt.1) then
         call qs_partition_int(A,imark)
         call quick_sort_int(A(     :imark-1))
         call quick_sort_int(A(imark:       ))
      end if
   end subroutine quick_sort_int
   
   
   !> Private sorting routine for both real and integer arrays
   subroutine qs_partition(A,B,marker)
      implicit none
      real(WP), dimension(:), intent(inout) :: A
      integer , dimension(:), intent(inout) :: B
      integer , intent(out) :: marker
      integer :: i,j,itmp
      real(WP) :: dtmp,x
      x=A(1)
      i=0
      j=size(A)+1
      do
         j=j-1
         do
            if (A(j).le.x) exit
            j=j-1
         end do
         i=i+1
         do
            if (A(i).ge.x) exit
            i=i+1
         end do
         if (i.lt.j) then
            ! Exchange A(i) and A(j)
            dtmp=A(i)
            A(i)=A(j)
            A(j)=dtmp
            ! Also exchange B(i) and B(j)
            itmp=B(i)
            B(i)=B(j)
            B(j)=itmp
         else if (i.eq.j) then
            marker=i+1
            return
         else
            marker=i
            return
         endif
      end do
   end subroutine qs_partition
   
   
   ! Integer array only
   subroutine qs_partition_int(A,marker)
      implicit none
      integer, dimension(:), intent(inout) :: A
      integer, intent(out) :: marker
      integer :: i,j,itmp
      integer :: x
      x=A(1)
      i=0
      j=size(A)+1
      do
         j=j-1
         do
            if (A(j).le.x) exit
            j=j-1
         end do
         i=i+1
         do
            if (A(i).ge.x) exit
            i=i+1
         end do
         if (i.lt.j) then
            ! Exchange A(i) and A(j)
            itmp=A(i)
            A(i)=A(j)
            A(j)=itmp
         else if (i.eq.j) then
            marker=i+1
            return
         else
            marker=i
            return
         endif
      end do
   end subroutine qs_partition_int
   

end module quicksort