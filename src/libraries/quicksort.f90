module quicksort
   use precision
   implicit none
   public  :: quick_sort,quick_sort_real,quick_sort_int
contains
   
   
   !> Recursive quicksort algorithm: this routine recursively sorts a real
   !> array by increasing order and sorts a integer array at the same time
   recursive subroutine quick_sort(A,B)
      implicit none
      real(WP), dimension(:), intent(inout) :: A
      integer , dimension(:), intent(inout) :: B
      real(WP) :: x,t
      integer :: first=1,last
      integer :: i,j,n
      last=size(A,1)
      x=A((first+last)/2)
      i=first; j=last
      do
         do while (A(i).lt.x)
            i=i+1
         end do
         do while (x.lt.A(j))
            j=j-1
         end do
         if (i.ge.j) exit
         t=A(i); A(i)=A(j); A(j)=t
         n=B(i); B(i)=B(j); B(j)=n
         i=i+1; j=j-1
      end do
      if (first.lt.i-1 ) call quick_sort(A(first:i-1 ),B(first:i-1 ))
      if (j+1  .lt.last) call quick_sort(A(j+1  :last),B(j+1  :last))
   end subroutine quick_sort
   
   
   !> Recursive quicksort algorithm, single real array
   recursive subroutine quick_sort_real(A)
      implicit none
      real(WP), dimension(:), intent(inout) :: A
      real(WP) :: x,t
      integer :: first=1,last
      integer :: i,j
      last=size(A,1)
      x=A((first+last)/2)
      i=first; j=last
      do
         do while (A(i).lt.x)
            i=i+1
         end do
         do while (x.lt.A(j))
            j=j-1
         end do
         if (i.ge.j) exit
         t=A(i); A(i)=A(j); A(j)=t
         i=i+1; j=j-1
      end do
      if (first.lt.i-1 ) call quick_sort_real(A(first:i-1 ))
      if (j+1  .lt.last) call quick_sort_real(A(j+1  :last))
   end subroutine quick_sort_real
   
   
   !> Recursive quicksort algorithm, single integer array
   recursive subroutine quick_sort_int(A)
      implicit none
      integer, dimension(:), intent(inout) :: A
      integer :: x,t
      integer :: first=1,last
      integer :: i,j
      last=size(A,1)
      x=A((first+last)/2)
      i=first; j=last
      do
         do while (A(i).lt.x)
            i=i+1
         end do
         do while (x.lt.A(j))
            j=j-1
         end do
         if (i.ge.j) exit
         t=A(i); A(i)=A(j); A(j)=t
         i=i+1; j=j-1
      end do
      if (first.lt.i-1 ) call quick_sort_int(A(first:i-1 ))
      if (j+1  .lt.last) call quick_sort_int(A(j+1  :last))
   end subroutine quick_sort_int
   
   
end module quicksort