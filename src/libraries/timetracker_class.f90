module timetracker_class
   use precision, only: WP
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: timetracker
   
   !> Timetracker object
   type :: timetracker
      integer  ::   n,  nmax                    !< Current and max time step
      real(WP) ::   t,  tmax                    !< Current and max time
      real(WP) ::  dt, dtmax                    !< Current and max time step size
      real(WP) :: cfl,cflmax                    !< Current and max CFL
      real(WP) ::  wt, wtmax                    !< Current and max wallclock time
   contains
      procedure
   end type timetracker
   
   
   !> Declare timetracker constructor
   interface timetracker
      procedure constructor
   end interface timetracker
   
   
contains
   
   
   !> Constructor for timetracker object
   function constructor() result(self)
      implicit none
      type(timetracker) :: self
      
      
      
   end subroutine constructor
   
   
end module timetracker_class
