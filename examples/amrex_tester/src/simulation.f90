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
   subroutine mak_lvl_init(lvl,time,pba,pdm) bind(c)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_multifab_build,amrex_mfiter_build,amrex_mfiter_destroy
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: pba
      type(c_ptr), intent(in), value :: pdm
      ! This is where data creation from scratch is done
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(amrex_mfiter)    :: mfi
      type(amrex_box)       :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: phi
      real(WP) :: x,y,z,r2
      integer :: i,j,k
      ! Rebuild data on new grids
      ba=pba; dm=pdm
      call sc%clr_lvl(lvl)
      call amrex_multifab_build(sc%SC(lvl),ba,dm,sc%ncomp,sc%nover)
      ! Build iterator
      call amrex_mfiter_build(mfi,sc%SC(lvl))
      ! Initialize
      do while (mfi%next())
         ! Get box and data
         bx=mfi%tilebox()
         phi=>sc%SC(lvl)%dataptr(mfi)
         ! Loop
         do k=bx%lo(3),bx%hi(3)
            do j=bx%lo(2),bx%hi(2) 
               do i=bx%lo(1),bx%hi(1)
                  ! Get position
                  x=sc%amr%xlo+(real(i,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(1)
                  y=sc%amr%ylo+(real(j,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(2)
                  z=sc%amr%zlo+(real(k,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(3)
                  ! Evaluate data
                  r2=((x-0.5_WP)**2+(y-0.75_WP)**2+(z-0.5_WP)**2)/0.01_WP
                  phi(i,j,k,1)=1.0_WP+exp(-r2)
               end do
            end do
         end do
      end do
      ! Destroy iterator
      call amrex_mfiter_destroy(mfi)
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
      print*,'please dont call me yet'
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
      print*,'please dont call me yet'
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
      
      ! Initialize user-defined functions
      call amr%register_udf(mak_lvl_init,mak_lvl_crse,mak_lvl_remk,clr_lvl,err_est)
      
      ! Initialize data
      call amr%initialize_data(0.0_WP)
      !call averagedown()
      
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
