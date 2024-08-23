!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision, only: WP
   use amrconfig_class,   only: amrconfig,mak_lvl_stype,clr_lvl_stype,err_est_stype
   use timetracker_class, only: timetracker
   use amrscalar_class,   only: amrscalar
   use amrensight_class,  only: amrensight
   use amrdata_class,     only: amrdata
   use event_class,       only: event
   implicit none
   private
   
   !> Public declarations
   public :: simulation_init,simulation_run,simulation_final
   
   !> Give ourselves an amrconfig object
   type(amrconfig) :: amr
   
   !> Also create a timetracker and an amrscalar
   type(timetracker) :: time
   type(amrscalar)   :: sc
   
   !> Event for regriding
   type(event) :: rgd_evt
   
   !> Define a velocity field for transport
   type(amrdata) :: U,V,W
   
   !> Ensight output
   type(event)      :: ens_evt
   type(amrensight) :: ens_out
   
contains
   
   
   !> User-provided routine to initialize level data
   subroutine initialize_lvl(lvl,t,pba,pdm) bind(c)
      use iso_c_binding, only: c_ptr
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: t
      type(c_ptr), intent(in), value :: pba,pdm
      
      ! Recreate velocity at our level
      call U%create_lvl(lvl,t,pba,pdm)
      call V%create_lvl(lvl,t,pba,pdm)
      call W%create_lvl(lvl,t,pba,pdm)
      
      ! Create level for sc
      call sc%create_lvl(lvl,t,pba,pdm)
      
      ! Initialize scalar value
      init_sc: block
         use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
         type(amrex_boxarray)  :: ba
         type(amrex_distromap) :: dm
         type(amrex_mfiter)    :: mfi
         type(amrex_box)       :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: mySC
         real(WP) :: x,y,z,r2
         integer :: i,j,k
         ! Convert pointers
         ba=pba; dm=pdm
         ! Loop over my boxes
         call amrex_mfiter_build(mfi,sc%SC(lvl))
         do while (mfi%next())
            bx=mfi%tilebox()
            mySC=>sc%SC(lvl)%dataptr(mfi)
            ! Loop inside box
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Get position
               x=sc%amr%xlo+(real(i,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(1)
               y=sc%amr%ylo+(real(j,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(2)
               z=sc%amr%zlo+(real(k,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(3)
               ! Evaluate data
               r2=((x-0.5_WP)**2+(y-0.75_WP)**2+(z-0.5_WP)**2)/0.01_WP
               mySC(i,j,k,1)=1.0_WP+exp(-r2)
               mySC(i,j,k,2)=1.0_WP-exp(-r2)
            end do; end do; end do
         end do
         call amrex_mfiter_destroy(mfi)
      end block init_sc
      
   end subroutine initialize_lvl
   
   
   !> User-provided routine to refine level
   subroutine refine_lvl(lvl,t,pba,pdm) bind(c)
      use iso_c_binding, only: c_ptr
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: t
      type(c_ptr), intent(in), value :: pba,pdm
      
      ! Recreate velocity at our level
      call U%create_lvl(lvl,t,pba,pdm)
      call V%create_lvl(lvl,t,pba,pdm)
      call W%create_lvl(lvl,t,pba,pdm)
      
      ! Refine level for sc
      call sc%refine_lvl(lvl,t,pba,pdm)
      
   end subroutine refine_lvl
   
   
   !> User-provided routine to remake level
   subroutine remake_lvl(lvl,t,pba,pdm) bind(c)
      use iso_c_binding, only: c_ptr
      implicit none
      integer,     intent(in), value :: lvl
      real(WP),    intent(in), value :: t
      type(c_ptr), intent(in), value :: pba,pdm
      
      ! Recreate velocity at our level
      call U%create_lvl(lvl,t,pba,pdm)
      call V%create_lvl(lvl,t,pba,pdm)
      call W%create_lvl(lvl,t,pba,pdm)
      
      ! Remake level for sc
      call sc%remake_lvl(lvl,t,pba,pdm)
      
   end subroutine remake_lvl
   
   
   !> User-provided routine to delete level data
   subroutine delete_lvl(lvl) bind(c)
      implicit none
      integer, intent(in), value :: lvl
      
      ! Delete velocity
      call U%delete_lvl(lvl)
      call V%delete_lvl(lvl)
      call W%delete_lvl(lvl)

      ! Delete level for sc
      call sc%delete_lvl(lvl)
      
   end subroutine delete_lvl
   
   
   !> User-provided routine to tag cells for refinement
   subroutine tag_lvl(lvl,ptag,t,tagval,clrval) bind(c)
      use iso_c_binding, only: c_ptr,c_char
      implicit none
      integer,                intent(in), value :: lvl
      type(c_ptr),            intent(in), value :: ptag
      real(WP),               intent(in), value :: t
      character(kind=c_char), intent(in), value :: tagval
      character(kind=c_char), intent(in), value :: clrval
      
      ! User-provided tagging for mesh adaptation
      tag_sc: block
         use amrex_amr_module, only: amrex_tagboxarray,amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
         type(amrex_tagboxarray) :: tag
         type(amrex_mfiter)      :: mfi
         type(amrex_box)         :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: phi
         character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: mytag
         integer :: i,j,k
         real(WP) :: x,y,z,r2
         real(WP), dimension(3) :: ctr
         ! Convert pointer
         tag=ptag
         ! Loop over my boxes
         call amrex_mfiter_build(mfi,sc%SC(lvl))
         do while(mfi%next())
            bx=mfi%tilebox()
            phi=>sc%SC(lvl)%dataptr(mfi)
            mytag=>tag%dataptr(mfi)
            ! Loop inside box
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Based on position
               x=sc%amr%xlo+(real(i,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(1)
               y=sc%amr%ylo+(real(j,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(2)
               z=sc%amr%zlo+(real(k,WP)+0.5_WP)*sc%amr%geom(lvl)%dx(3)
               ctr=[0.5_WP,0.5_WP,0.5_WP]+0.25_WP*[cos(t),sin(t),0.0_WP]
               r2=(x-ctr(1))**2+(y-ctr(2))**2+(z-ctr(3))**2
               if (r2.lt.0.01_WP) then
                  mytag(i,j,k,1)=tagval
               else
                  mytag(i,j,k,1)=clrval
               end if
               ! Based on scalar value
               !if (phi(i,j,k,1).ge.1.15_WP) mytag(i,j,k,1)=tagval
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
         use param, only: param_read,param_getsize
         ! Read in mesh parameters from input
         call param_read('Base nx',amr%nx); amr%xlo=0.0_WP; call param_read('Lx',amr%xhi); amr%xper=.true.
         call param_read('Base ny',amr%ny); amr%ylo=0.0_WP; call param_read('Ly',amr%yhi); amr%yper=.true.
         call param_read('Base nz',amr%nz); amr%zlo=0.0_WP; call param_read('Lz',amr%zhi); amr%zper=.true.
         call param_read('Number of levels',amr%nlvl); amr%nlvl=amr%nlvl-1
         allocate(amr%rref(param_getsize('Refinement ratio'))); call param_read('Refinement ratio',amr%rref)
         call param_read('Max block size',amr%nmax,default=32)
         call param_read('Blocking fator',amr%nbloc,default=8)
         ! Build amrconfig
         call amr%initialize(name='amrtest')
         ! Create event for regrid
         rgd_evt=event(time=time,name='Regrid')
         call param_read('Steps between regrid',rgd_evt%nper)
      end block build_amr
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      ! Create scalar solver
      create_scalar_solver: block
         call sc%initialize(amr=amr,nscalar=2,name='scalar')
         sc%SCname=['phi ','Zmix']
      end block create_scalar_solver
      
      ! Create velocity data
      create_velocity_data: block
         call U%initialize(amr=amr,ncomp=1,nover=0,atface=[.true.,.false.,.false.])
         call V%initialize(amr=amr,ncomp=1,nover=0,atface=[.false.,.true.,.false.])
         call W%initialize(amr=amr,ncomp=1,nover=0,atface=[.false.,.false.,.true.])
      end block create_velocity_data
      
      ! Initialize user-defined procedures
      call amr%register_udp(initialize_lvl,refine_lvl,remake_lvl,delete_lvl,tag_lvl)
      
      ! Initialize grid and data
      call amr%initialize_data(time=time%t)
      call amr%average_down(sc%SC)
      
      
      call amr%print()
      
      ! Prepare Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Create Ensight output from amrconfig
         call ens_out%initialize(amr=amr,name='amrtest')
         ! Add variables to output
         call ens_out%add_scalar(name=sc%SCname(1),scalar=sc%SC,comp=1)
         call ens_out%add_scalar(name=sc%SCname(2),scalar=sc%SC,comp=2)
         ! Output to ensight
         call ens_out%write_data(time=time%t)
      end block create_ensight
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         !call fs%get_cfl(time%dt,time%cfl)
         !call time%adjust_dt()
         call time%increment()
         
         ! Time advance scalar solver
         
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time=time%t)
         
         ! Regrid if needed
         if (rgd_evt%occurs()) call amr%regrid(baselvl=0,time=time%t)
         
      end do
      
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
