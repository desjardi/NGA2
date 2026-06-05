!> GID -> LID hash for fast resolution of bond endpoints during force evaluation.
!>
!> Used by amrpd's compute_dilatation and compute_force kernels: after each
!> fill_neighbors_radius call, walk the tile's (owned + ghost) particles and
!> build the hash; then each bond looks up its two endpoints by their GID
!> (packed (id, cpu) as a unique int64 key via amrpd_get_particle_idcpu).
!>
!> Implementation: sorted-array + binary search. Build O(N log N), lookup
!> O(log N), where N is the tile-local particle count (owned + ghost). Cheap
!> enough that the hash is rebuilt per tile per force evaluation; no need to
!> cache across timesteps.
module amrpd_hash_class
   use iso_c_binding, only: c_int64_t
   implicit none
   private

   public :: gid_hash

   !> Sorted (key, val) pairs. Key is a unique int64 GID; val is the 1-based
   !> local index into the source particle array.
   type :: gid_hash
      integer(c_int64_t), allocatable :: keys(:)
      integer,            allocatable :: vals(:)
      integer :: n = 0
   contains
      procedure :: build
      procedure :: lookup
      procedure :: lookup_range   !< For periodic-image disambiguation: returns ALL duplicates of a key
      procedure :: finalize
   end type gid_hash

contains

   !> Build a sorted hash from an array of keys. Values are assigned 1..n
   !> (the LIDs in the source array). Caller supplies the key array; this
   !> routine copies and sorts.
   subroutine build(this,n,keys)
      implicit none
      class(gid_hash), intent(inout) :: this
      integer, intent(in) :: n
      integer(c_int64_t), intent(in) :: keys(n)
      integer :: i
      call this%finalize()
      this%n = n
      if (n.gt.0) then
         allocate(this%keys(n),this%vals(n))
         this%keys = keys
         do i = 1, n
            this%vals(i) = i
         end do
         call quicksort_pair(this%keys,this%vals,1,n)
      end if
   end subroutine build

   !> Look up a key. Returns the 1-based LID on hit, -1 on miss.
   pure function lookup(this,key) result(lid)
      implicit none
      class(gid_hash), intent(in) :: this
      integer(c_int64_t), intent(in) :: key
      integer :: lid
      integer :: lo,hi,mid
      lid = -1
      if (this%n.eq.0) return
      lo = 1; hi = this%n
      do while (lo.le.hi)
         mid = (lo + hi) / 2
         if (this%keys(mid).lt.key) then
            lo = mid + 1
         else if (this%keys(mid).gt.key) then
            hi = mid - 1
         else
            lid = this%vals(mid)
            return
         end if
      end do
   end function lookup

   !> Find the contiguous bracket of duplicates for a given key in the sorted
   !> array. Returns first_idx (1-based) and n_dup. On miss, n_dup = 0.
   !>
   !> Use case: periodic-image disambiguation. When the hash is built from a
   !> particle array that contains both an owned particle and its periodic-
   !> image ghost copy (which share the same idcpu = key), multiple entries
   !> exist. The caller walks the bracket [first_idx .. first_idx+n_dup-1]
   !> in self%vals to get all candidate LIDs, then picks the right image by
   !> minimum-image distance to an anchor position.
   !>
   !> Common case (no duplicates): n_dup = 1, self%vals(first_idx) is the LID.
   pure subroutine lookup_range(this,key,first_idx,n_dup)
      implicit none
      class(gid_hash), intent(in) :: this
      integer(c_int64_t), intent(in) :: key
      integer, intent(out) :: first_idx,n_dup
      integer :: lo,hi,mid,i,j
      first_idx = -1; n_dup = 0
      if (this%n.eq.0) return
      ! Binary search for any matching index
      lo = 1; hi = this%n
      mid = -1
      do while (lo.le.hi)
         mid = (lo + hi) / 2
         if (this%keys(mid).lt.key) then
            lo = mid + 1
         else if (this%keys(mid).gt.key) then
            hi = mid - 1
         else
            exit
         end if
      end do
      if (mid.lt.1.or.mid.gt.this%n) return
      if (this%keys(mid).ne.key) return
      ! Scan left and right for duplicates (sorted -> contiguous)
      i = mid
      do while (i.gt.1)
         if (this%keys(i-1).ne.key) exit
         i = i - 1
      end do
      j = mid
      do while (j.lt.this%n)
         if (this%keys(j+1).ne.key) exit
         j = j + 1
      end do
      first_idx = i
      n_dup     = j - i + 1
   end subroutine lookup_range

   !> Release allocated storage.
   subroutine finalize(this)
      implicit none
      class(gid_hash), intent(inout) :: this
      if (allocated(this%keys)) deallocate(this%keys)
      if (allocated(this%vals)) deallocate(this%vals)
      this%n = 0
   end subroutine finalize


   !> Recursive Hoare-partition quicksort on (key, val) pairs, sorted by key.
   !> Private module helper.
   recursive subroutine quicksort_pair(keys,vals,lo,hi)
      implicit none
      integer(c_int64_t), intent(inout) :: keys(:)
      integer,            intent(inout) :: vals(:)
      integer, intent(in) :: lo,hi
      integer :: i,j,tv
      integer(c_int64_t) :: pivot,tk
      if (lo.ge.hi) return
      pivot = keys((lo + hi) / 2)
      i = lo; j = hi
      do
         do while (keys(i).lt.pivot); i = i + 1; end do
         do while (keys(j).gt.pivot); j = j - 1; end do
         if (i.le.j) then
            tk = keys(i); keys(i) = keys(j); keys(j) = tk
            tv = vals(i); vals(i) = vals(j); vals(j) = tv
            i = i + 1; j = j - 1
         end if
         if (i.gt.j) exit
      end do
      call quicksort_pair(keys,vals,lo,j)
      call quicksort_pair(keys,vals,i,hi)
   end subroutine quicksort_pair

end module amrpd_hash_class
