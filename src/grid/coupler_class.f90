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
      ! This is our two grids
      class(pgrid), pointer :: pg1,pg2                    !< These are the pgrids that will be coupled
   contains
      
   end type coupler
   
   
   !> Declare coupler constructor
   interface coupler
      procedure construct_from_two_pgrids
   end interface coupler
   
contains
   
   
   !> Coupler constructor from two pgrids
   function construct_from_two_pgrids(pg1,pg2,name) result(self)
      use string,   only: lowercase
      use messager, only: die
      use parallel, only: comm
      implicit none
      type(coupler) :: self
      class(pgrid), target, intent(in) :: pg1,pg2
      character(len=*), intent(in) :: name
      integer :: ierr
      
      ! Set the name for the coupler
      self%name=trim(adjustl(name))
      
      ! Point to pgrid objects
      self%pg1=>pg1
      self%pg2=>pg2
      
      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%pg1%amRoot) then
            write(message,'("Coupler [",a,"] between pgrids [",a,"] and [",a,"]")') trim(self%name),trim(self%pg1%name),trim(self%pg2%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function construct_from_two_pgrids
   
   
end module coupler_class
