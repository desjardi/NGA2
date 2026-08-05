!> Various definitions and tools for running an NGA2 simulation
module simulation
   use avl_trees, only: avl_tree_t,avl_insert,avl_retrieve,int_cast,real_cast,avl_delete_all
   implicit none
   private
   
   ! AVL Binary Tree
   type(avl_tree_t) :: tree

   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   !> Initialization the NGA2 simulation
   subroutine simulation_init
      use, intrinsic :: iso_fortran_env, only: output_unit
      use, non_intrinsic :: avl_trees

      implicit none

      integer, parameter :: keys_count = 20

      type(avl_tree_t) :: tree
      logical :: found
      class(*), allocatable :: retval
      integer :: the_keys(1:keys_count)
      integer :: i, j

      do i = 1, keys_count
         the_keys(i) = i
      end do
      call fisher_yates_shuffle (the_keys, keys_count)

      call avl_check (tree)
      do i = 1, keys_count
         call avl_insert (lt, the_keys(i), real (the_keys(i)+10), tree)
         call avl_check (tree)
         if (avl_size (tree) /= i) error stop
         do j = 1, keys_count
            if (avl_contains (lt, the_keys(j), tree) .neqv. (j <= i)) error stop
         end do
         do j = 1, keys_count
            call avl_retrieve (lt, the_keys(j), tree, found, retval)
            if (found .neqv. (j <= i)) error stop
            if (found) then
               ! This crazy way to write ‘/=’ is to quell those tiresome
               ! warnings about using ‘==’ or ‘/=’ with floating point
               ! numbers. Floating point numbers can represent integers
               ! *exactly*.
               !if (0 < abs (real_cast (retval) - real (the_keys(j)))) error stop
               print*,real_cast(retval),the_keys(j)
            end if
            ! if (found) then
            !    block
            !       character(len = 1), parameter :: ch = '*'
            !       !
            !       ! Try replacing the data with a character and then
            !       ! restoring the number.
            !       !
            !       call avl_insert (lt, the_keys(j), ch, tree)
            !       call avl_retrieve (lt, the_keys(j), tree, found, retval)
            !       if (.not. found) error stop
            !       if (char_cast (retval) /= ch) error stop
            !       call avl_insert (lt, the_keys(j), real (the_keys(j)+10), tree)
            !       call avl_retrieve (lt, the_keys(j), tree, found, retval)
            !       if (.not. found) error stop
            !       !if (0 < abs (real_cast (retval) - real (the_keys(j)+10))) error stop
            !    end block
            ! end if
         end do
      end do

      write (output_unit, '(70("-"))')
      call avl_write (int_real_writer, output_unit, tree)
      write (output_unit, '(70("-"))')
      call print_contents (output_unit, tree)
      write (output_unit, '(70("-"))')

      call fisher_yates_shuffle (the_keys, keys_count)
      do i = 1, keys_count
         call avl_delete (lt, the_keys(i), tree)
         call avl_check (tree)
         if (avl_size (tree) /= keys_count - i) error stop
         ! Try deleting a second time.
         call avl_delete (lt, the_keys(i), tree)
         call avl_check (tree)
         if (avl_size (tree) /= keys_count - i) error stop
         do j = 1, keys_count
            if (avl_contains (lt, the_keys(j), tree) .neqv. (i < j)) error stop
         end do
         do j = 1, keys_count
            call avl_retrieve (lt, the_keys(j), tree, found, retval)
            if (found .neqv. (i < j)) error stop
            if (found) then
               !if (0 < abs (real_cast (retval) - real (the_keys(j)))) error stop
            end if
         end do
      end do

      ! Remove all nodes
      call avl_delete_all(tree)

      ! Remake the tree
      do i = 1, keys_count
         call avl_insert (lt, the_keys(i), real (the_keys(i)+100), tree)
      end do

      write (output_unit, '(70("-"))')
      call print_contents (output_unit, tree)

      contains

      subroutine fisher_yates_shuffle (keys, n)
         integer, intent(inout) :: keys(*)
         integer, intent(in) :: n

         integer :: i, j
         real :: randnum
         integer :: tmp

         do i = 1, n - 1
            call random_number (randnum)
            j = i + floor (randnum * (n - i + 1))
            tmp = keys(i)
            keys(i) = keys(j)
            keys(j) = tmp
         end do
      end subroutine fisher_yates_shuffle

      function lt (u, v) result (u_lt_v)
         class(*), intent(in) :: u, v
         logical :: u_lt_v

         select type (u)
         type is (integer)
            select type (v)
            type is (integer)
               u_lt_v = (u < v)
            class default
               ! This case is not handled.
               error stop
            end select
         class default
            ! This case is not handled.
            error stop
         end select
      end function lt

      subroutine int_real_writer (unit, key, data)
         integer, intent(in) :: unit
         class(*), intent(in) :: key, data

         write (unit, '("(", I0, ", ", F0.1, ")")', advance = 'no') &
               & int_cast(key), real_cast(data)
      end subroutine int_real_writer

      subroutine print_contents (unit, tree)
         integer, intent(in) :: unit
         class(avl_tree_t), intent(in) :: tree

         type(avl_pointer_pair_t), pointer :: ppairs, pp

         write (unit, '("tree size = ", I0)') avl_size (tree)
         ppairs => avl_pointer_pairs (tree)
         pp => ppairs
         do while (associated (pp))
            write (unit, '("(", I0, ", ", F0.1, ")")') &
                  & int_cast (pp%p_key), real_cast (pp%p_data)
            pp => pp%next
         end do
         if (associated (ppairs)) deallocate (ppairs)
      end subroutine print_contents

   end subroutine simulation_init
   
   
   !> Run the NGA2 simulation
   subroutine simulation_run
      implicit none
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
   end subroutine simulation_final
   
end module simulation
