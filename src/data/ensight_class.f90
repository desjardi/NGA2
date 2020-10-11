!> Ensight class concept is defined here: given a config object,
!> it provides parallel I/O access to an ensight file
module ensight_class
   use precision,    only: WP
   use string,       only: str_short,str_medium
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ensight
   
   ! List types
   type :: scl
      type(scl), pointer :: next
      character(len=str_short) :: name
      real(WP), dimension(:,:,:), pointer :: scl_ptr
   end type scl
   type :: vct
      type(vct), pointer :: next
      character(len=str_short) :: name
      real(WP), dimension(:,:,:), pointer :: vctx_ptr
      real(WP), dimension(:,:,:), pointer :: vcty_ptr
      real(WP), dimension(:,:,:), pointer :: vctz_ptr
   end type vct
   
   !> Ensight object definition as list of pointers to arrays
   type :: ensight
      ! An ensight object has a case file
      character(len=str_medium) :: casename                           !< Name of casefile to read/write
      ! An ensight object stores time values
      integer :: ntime                                                !< Number of scalar values
      real(WP), dimension(:), allocatable :: time                     !< Time values
      ! An ensight object stores geometry data
      type(config), pointer :: cfg                                    !< Config for ensight geometry and parallel I/O
      character(len=str_medium) :: geomname                           !< Name of geometry file to write
      ! An ensight object stores lists of pointers to data
      type(scl), pointer :: first_scl                                 !< Scalar list
      type(vct), pointer :: first_vct                                 !< Vector list
   contains
      
   end type ensight
   
   
   !> Declare ensight constructor
   interface ensight
      procedure construct_ensight
   end interface ensight
   
   
contains
   
   !> Constructor for an empty ensight object
   function construct_ensight(cfg,casename) result(self)
      use monitor,  only: die
      implicit none
      type(ensight) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: casename
      
      ! Link to config
      self%cfg=>cfg
      
      ! Store casename
      self%casename=trim(adjustl(casename))
      
      
      
   end function construct_ensight
   
   
   
   
end module ensight_class
