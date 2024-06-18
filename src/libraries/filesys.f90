module filesys
   use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
   implicit none
   private
   interface
      !> Creates a directory with all possible permitions (except those in umask)
      subroutine makedir_c(dirname) bind(C,name='c_mkdir')
         import c_char
         implicit none
         character(kind=c_char), intent(in) :: dirname
      end subroutine makedir_c
      !> Decides whether a given file name is a directory.
      function isdir_c(dirname) bind(c, name='c_isdir') result(res)
         import :: c_int,c_char
         character(kind=c_char), intent(in) :: dirname
         integer(c_int) :: res
      end function isdir_c
   end interface
   ! Routines
   public :: makedir,isdir
contains
   !> Directory creation
   subroutine makedir(dirname)
      implicit none
      character(*), intent(in) :: dirname
      character(len=len_trim(dirname)+1,kind=c_char) :: cstring
      cstring=trim(dirname)//c_null_char
      call makedir_c(cstring)
   end subroutine makedir
   !> Checks existence of directory
   function isdir(dirname) result(res)
      character(*), intent(in) :: dirname
      character(len=len_trim(dirname)+1,kind=c_char) :: cstring
      logical :: res
      integer(c_int) :: status
      cstring=trim(dirname)//c_null_char
      status=isdir_c(cstring)
      res=(status.ne.0)
    end function isdir
end module filesys