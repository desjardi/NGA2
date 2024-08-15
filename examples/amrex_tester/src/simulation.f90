!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use amrconfig_class, only: amrconfig,mak_lvl_stype,clr_lvl_stype,err_est_stype
   use amrscalar_class, only: amrscalar
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Give ourselves an amrconfig object
   type(amrconfig) :: amr

   !> Also create an amrscalar
   type(amrscalar) :: sc
   
contains
   
   
   !> User-provided data creation from scratch
   subroutine mak_lvl_init(lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: ba
      type(c_ptr), intent(in), value :: dm
      ! This is where data creation from scratch is done
   end subroutine mak_lvl_init
   
   
   !> User-provided data creation from coarse data
   subroutine mak_lvl_crse(lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: ba
      type(c_ptr), intent(in), value :: dm
      ! This is where data creation from coarse data is done
   end subroutine mak_lvl_crse
   
   
   !> User-provided data creation from current and coarse data
   subroutine mak_lvl_remk(lvl,time,ba,dm) bind(c)
      use iso_c_binding, only: c_ptr
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: ba
      type(c_ptr), intent(in), value :: dm
      ! This is where data creation from current and coarse data is done
   end subroutine mak_lvl_remk
   
   
   !> User-provided data deletion
   subroutine clr_lvl(lvl) bind(c)
      implicit none
      integer, intent(in) , value :: lvl
      ! This is where data deletion is done
      call sc%clr_lvl(lvl)
   end subroutine clr_lvl
   
   
   !> User-provided error estimation for regriding
   subroutine err_est(lvl,tags,time,tagval,clrval) bind(c)
      use iso_c_binding, only: c_ptr,c_char
      implicit none
      integer,                intent(in), value :: lvl
      type(c_ptr),            intent(in), value :: tags
      real(WP),               intent(in), value :: time
      character(kind=c_char), intent(in), value :: tagval
      character(kind=c_char), intent(in), value :: clrval
      ! This is where error estimation for regriding is done
   end subroutine err_est

   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      ! Define domain and creare amrcore object
      build_amr: block
         amr%nx  =64    ; amr%ny  =64    ; amr%nz  =64
         amr%xlo =0.0_WP; amr%ylo =0.0_WP; amr%zlo =0.0_WP
         amr%xhi =1.0_WP; amr%yhi =1.0_WP; amr%zhi =1.0_WP
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         amr%nlvl=3
         call amr%initialize(name='mybox')
      end block build_amr
      
      ! Create scalar solver
      create_scalar_solver: block
         call sc%initialize(amr=amr)
      end block create_scalar_solver
      
      ! Initialize regriding functions
      call amr%init_regrid_functions(mak_lvl_init,mak_lvl_crse,mak_lvl_remk,clr_lvl,err_est)
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      use amrconfig_class, only: finalize_amrex
      implicit none
      call sc%finalize()
      call amr%finalize()
      call finalize_amrex()
   end subroutine simulation_final
   
   
end module simulation
