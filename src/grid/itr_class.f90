!> Iterator concept is defined here.
!> The intention is to have the ability to work in a sub-region of a pgrid.
!> This is directly useful for boundary conditions, but probably for more.
module itr_class
   use precision,   only: WP
   use string,      only: str_medium
   use pgrid_class, only: pgrid
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: itr
   
   !> Bcond object definition
   !> A "bcond" is a type, direction and collection of cells indices
   type :: itr
      character(len=str_medium) :: name='UNNAMED_ITR'     !< Iterator name (default=UNNAMED_ITR)
      integer :: no_,n_                                   !< Total number of local cells in iterator, with and without overlap
      integer, dimension(:,:), allocatable :: map         !< This is our map for the iterator
      class(pgrid), pointer :: pgrid                      !< This is the pgrid the iterator is build for
   contains
      procedure :: print=>itr_print                       !< Output iterator to the screen
   end type itr
   
   
   !> Declare single-grid iterator constructor
   interface itr
      procedure construct_from_function
   end interface itr
   
contains

   !> Iterator constructor from a tester function
   function construct_from_function(pg,name,tester) result(self)
      implicit none
      type(itr) :: self
      class(pgrid), target, intent(in) :: pg
      character(len=*), intent(in) :: name
      interface
         logical function tester(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function tester
      end interface
      integer :: i,j,k
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%pgrid=>pg
      
      ! Zero out counters
      self%n_=0; self%no_=0
      ! Loop over local domain with overlap to count bcond cells
      do k=pg%kmino_,pg%kmaxo_
         do j=pg%jmino_,pg%jmaxo_
            do i=pg%imino_,pg%imaxo_
               ! Apply condition for being in the bcond
               if (tester(pg,i,j,k)) then
                  ! Cell belongs in bcond, increment counter with overlap
                  self%no_=self%no_+1
                  ! If it is not overlap, increment interior counter too
                  !if () then
                  !   self%n_=self%n_+1
                  !end if
               end if
            end do
         end do
      end do
      
   end function construct_from_function
   
   !> Basic printing of iterator
   subroutine itr_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(itr), intent(in) :: this
      integer :: i,ierr
      if (this%pgrid%amRoot) write(output_unit,'("Iterator ",a," for pgrid ",a)') trim(this%name),trim(this%pgrid%name)
      do i=0,this%pgrid%nproc-1
         ! Block for clean output
         call MPI_BARRIER(this%pgrid%comm,ierr)
         ! Output info
         if (this%pgrid%rank.eq.i) then
            write(output_unit,'(" >>> Rank ",i0," number of cells = ",i0)') this%no_
         end if
      end do
   end subroutine itr_print
   
end module itr_class
