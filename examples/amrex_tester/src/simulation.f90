!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use amrcore_class, only: amrcore
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Give ourselves an amrcore object
   type(amrcore) :: amr

contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      ! Define domain and creare amrcore object
      build_amr: block
         amr%nx  =64    ; amr%ny =64    ; amr%nz =64
         amr%xlo =0.0_WP; amr%ylo=0.0_WP; amr%zlo=0.0_WP
         amr%xhi =1.0_WP; amr%yhi=1.0_WP; amr%zhi=1.0_WP
         amr%xper=1     ; amr%yper=1    ; amr%zper=1
         amr%nlvl=3
         call amr%initialize()
      end block build_amr
      

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      use amrcore_class, only: finalize_amrex
      implicit none
      call amr%finalize()
      call finalize_amrex()
   end subroutine simulation_final
   
   
end module simulation
