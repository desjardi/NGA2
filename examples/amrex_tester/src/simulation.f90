!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use amrconfig_class,  only: amrconfig,mak_lvl_stype,clr_lvl_stype,err_est_stype
   use amrscalar_class,  only: amrscalar
   use amrensight_class, only: amrensight
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Give ourselves an amrconfig object
   type(amrconfig) :: amr

   !> Also create an amrscalar
   type(amrscalar) :: sc

   !> Ensight output
   type(amrensight) :: ens_out
   
contains
   
   
   !> User-provided routine to initialize level data
   subroutine initialize_lvl(lvl,time,pba,pdm) bind(c)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: pba
      type(c_ptr), intent(in), value :: pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      
      ! Convert pointers
      ba=pba; dm=pdm
      
      ! Delete level
      call delete_lvl(lvl)
      
      ! Create level for sc
      call sc%create_lvl(lvl,ba,dm)
      
      ! Initialize scalar value
      init_sc: block
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
         type(amrex_mfiter) :: mfi
         type(amrex_box)    :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: phi
         real(WP) :: x,y,z,r2
         integer :: i,j,k
         ! Loop over my boxes
         call amrex_mfiter_build(mfi,sc%SC(lvl))
         do while (mfi%next())
            bx=mfi%tilebox()
            phi=>sc%SC(lvl)%dataptr(mfi)
            ! Loop inside box
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Get position
               x=sc%amr%xlo+(real(i,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(1)
               y=sc%amr%ylo+(real(j,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(2)
               z=sc%amr%zlo+(real(k,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(3)
               ! Evaluate data
               r2=((x-0.5_WP)**2+(y-0.75_WP)**2+(z-0.5_WP)**2)/0.01_WP
               phi(i,j,k,1)=1.0_WP+exp(-r2)
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end block init_sc
      
   end subroutine initialize_lvl
   
   
   !> User-provided routine to refine level
   subroutine refine_lvl(lvl,time,pba,pdm) bind(c)
      use iso_c_binding, only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: pba
      type(c_ptr), intent(in), value :: pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      
      ! Convert pointers
      ba=pba; dm=pdm
      
      ! Delete level
      call delete_lvl(lvl)
      
      ! Create level for sc and fill it
      call sc%create_lvl(lvl,ba,dm)
      call sc%fillcoarse(lvl,time)
      
   end subroutine refine_lvl
   
   
   !> User-provided routine to remake level
   subroutine remake_lvl(lvl,time,pba,pdm) bind(c)
      use iso_c_binding, only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: time
      type(c_ptr), intent(in), value :: pba
      type(c_ptr), intent(in), value :: pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      
      ! Convert pointers
      ba=pba; dm=pdm

      ! This is where data creation from current and coarse data is done
      print*,'please dont call me yet'
   end subroutine remake_lvl
   
   
   !> User-provided routine to delete level data
   subroutine delete_lvl(lvl) bind(c)
      implicit none
      integer, intent(in), value :: lvl
      
      ! Delete level for sc
      call sc%delete_lvl(lvl)
      
   end subroutine delete_lvl
   
   
   !> User-provided routine to tag cells for refinement
   subroutine tag_lvl(lvl,ptag,time,tagval,clrval) bind(c)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_tagboxarray
      implicit none
      integer,                intent(in), value :: lvl
      type(c_ptr),            intent(in), value :: ptag
      real(WP),               intent(in), value :: time
      character(kind=c_char), intent(in), value :: tagval
      character(kind=c_char), intent(in), value :: clrval
      type(amrex_tagboxarray) :: tag
      
      ! Convert pointers
      tag=ptag
      
      ! Tag based on current sc value
      tag_sc: block
         use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
         type(amrex_mfiter) :: mfi
         type(amrex_box)    :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: phi
         character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: mytag
         integer :: i,j,k
         ! Loop over my boxes
         call amrex_mfiter_build(mfi,sc%SC(lvl))
         do while(mfi%next())
            bx=mfi%tilebox()
            phi=>sc%SC(lvl)%dataptr(mfi)
            mytag=>tag%dataptr(mfi)
            ! Loop inside box
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (phi(i,j,k,1).ge.1.15_WP) mytag(i,j,k,1)=tagval
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end block tag_sc
      
   end subroutine tag_lvl
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      implicit none
      
      ! Define domain and creare amrcore object
      build_amr: block
         amr%nx  =64    ; amr%ny  =64    ; amr%nz  =64
         amr%xlo =0.0_WP; amr%ylo =0.0_WP; amr%zlo =0.0_WP
         amr%xhi =1.0_WP; amr%yhi =1.0_WP; amr%zhi =1.0_WP
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         amr%nlvl=1
         amr%nmax=24
         amr%rref=[4]
         call amr%initialize(name='amrtest')
      end block build_amr
      
      ! Create scalar solver
      create_scalar_solver: block
         call sc%initialize(amr=amr,nscalar=1,name='scalar')

      end block create_scalar_solver
      
      ! Initialize user-defined procedures
      call amr%register_udp(initialize_lvl,refine_lvl,remake_lvl,delete_lvl,tag_lvl)
      
      ! Initialize data
      call amr%initialize_data(0.0_WP)
      call amr%average_down(sc%SC)
      
      !call amr%regrid(baselvl=0,time=0.0_WP)
      
      call amr%print()
      
      ! Prepare Ensight output
      create_ensight: block
         ! Create Ensight output from armconfig
         call ens_out%initialize(amr=amr,name='amrtest')
         ! Add variables to output
         call ens_out%add_scalar('SC',sc%SC)
         ! Output to ensight
         call ens_out%write_data(0.0_WP)
      end block create_ensight
      
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
