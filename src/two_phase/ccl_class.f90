!> Connected component labeling class:
!> Provides support for identifying and manipulating Lagrangian objects from a vfs solution
module ccl_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ccl
   
   !> CCL object definition
   type :: ccl
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config object for which the CCL is build
      
      ! This is the name of the CCL
      character(len=str_medium) :: name='UNNAMED_CCL'     !< Solver name (default=UNNAMED_CCL)
      
      ! Maximum number of PLIC interfaces per cell
      integer :: max_interface_planes                     !< Number of planar interfaces per cell (0=VF-only, 1=PLIC, 2=R2P, etc)
      
      ! CCL selection parameters
      real(WP) :: VFlo                                    !< Minimum VF value considered for a structure to exist
      
      ! Output of the CCL
      integer, dimension(:,:,:), allocatable :: id        !< ID of the structure that contains the cell
      
   contains
      
   end type ccl
   
   
   !> Declare CCL algorithm constructor
   interface ccl
      procedure constructor
   end interface ccl
   
contains
   
   
   !> Default constructor for CCL algorithm
   function constructor(cfg,name) result(self)
      use messager, only: die
      implicit none
      type(ccl) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      
      ! Set the name for the object
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to cfg object
      self%cfg=>cfg
      
   end function constructor
   
   
   

end module ccl_class
