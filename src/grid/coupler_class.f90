!> Coupler concept is defined here: it takes in two pgrid
!> objects and builds the communication and interpolation
!> layer to exchange data between them.
module coupler_class
   use precision,      only: WP
   use string,         only: str_medium
   use pgrid_class,    only: pgrid
   use mpi_f08
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: coupler
   
   !> Coupler object definition
   type :: coupler
      ! This is the name of the coupler
      character(len=str_medium) :: name='UNNAMED_CPL'     !< Coupler name (default=UNNAMED_CPL)
      ! This is our source grid
      class(pgrid), pointer :: src                        !< This is the source pgrid, from which data will be sent
      class(pgrid), pointer :: dst                        !< This is the destination pgrid, to which data will be sent
      ! Source rank on which dst point have been located
      integer, dimension(:,:,:), allocatable :: srank     !< Source-rank of each destionation point
   contains
      
   end type coupler
   
   
   !> Declare coupler constructor
   interface coupler
      procedure construct_from_two_pgrids
   end interface coupler
   
contains
   
   
   !> Coupler constructor from two pgrids
   function construct_from_two_pgrids(src,dst,name) result(self)
      use string,   only: lowercase
      use messager, only: die
      use parallel, only: comm
      implicit none
      type(coupler) :: self
      class(pgrid), target, intent(in) :: src,dst
      character(len=*), intent(in) :: name
      integer :: i,j,k
      
      ! Set name for the coupler
      self%name=trim(adjustl(name))
      
      ! Point to src and dst pgrid objects
      self%src=>src
      self%dst=>dst
      
      ! Allocate srank array
      allocate(self%srank(self%dst%imino_:self%dst%imaxo_,self%dst%jmino_:self%dst%jmaxo_,self%dst%kmino_:self%dst%kmaxo_)); self%srank=-1
      
      ! Localize each dst cell center on the src grid - *this is from the perspective of the src*
      do k=self%dst%kmin_,self%dst%kmax_
         do j=self%dst%jmin_,self%dst%jmax_
            do i=self%dst%imin_,self%dst%imax_
               !self%srank(i,j,k)=self%src%get_rank(this%p(i)%ind)
            end do
         end do
      end do
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%src%amRoot) then
            write(message,'("Coupler [",a,"] from pgrid [",a,"] to pgrid [",a,"]")') trim(self%name),trim(self%src%name),trim(self%dst%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function construct_from_two_pgrids
   
   
end module coupler_class
