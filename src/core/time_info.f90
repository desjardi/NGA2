module time_info
   use precision, only: WP
   implicit none
   private
   
   integer , public :: ntime                 !< Current time step number
   real(WP), public :: time                  !< Current time
   real(WP), public :: dt                    !< Current time step size
   real(WP), public :: cfl                   !< Current CFL number
   real(WP), public :: wtime                 !< Current wallclock time
   
   integer , public :: max_ntime             !< Maximum number of time steps
   real(WP), public :: max_time              !< Maximum time
   real(WP), public :: max_dt                !< Maximum time step size
   real(WP), public :: max_cfl               !< Maximum CFL number
   real(WP), public :: max_wtime             !< Maximum wallclock time
   
end module time_info
