!> AMR-capable volume fraction solver class:
!> Performs semi-Lagrangian geometric advancement, curvature calculation, interface reconstruction
!> This implementation solves VF transport only at the finest level
module amrvfs_class
   use precision,             only: WP
   use string,                only: str_medium
   use amrconfig_class,       only: amrconfig
   use amrex_amr_module,      only: amrex_multifab
   !use irl_fortran_interface
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrvfs,VFhi,VFlo
   
   ! List of available interface reconstructions schemes for VF
   integer, parameter, public :: plicnet=0
   
   ! List of available interface transport schemes for VF
   integer, parameter, public :: remap=0                             !< Cell-based geometric transport (faster but fluxes are not available)
   
   ! IRL cutting moment calculation method
   integer, parameter, public :: recursive_simplex=0                 !< Recursive simplex cutting
   integer, parameter, public :: half_edge=1                         !< Half-edge cutting (default)
   integer, parameter, public :: nonrecurs_simplex=2                 !< Non-recursive simplex cutting
   
   ! Default parameters for volume fraction solver
   integer,  parameter :: nband=3                                    !< Number of cells around the interfacial cells on which localized work is performed
   integer,  parameter :: advect_band=1                              !< How far we do the transport
   real(WP), parameter :: VFlo=1.0e-12_WP                            !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo                           !< Maximum VF value considered
   real(WP), parameter :: volume_epsilon_factor =1.0e-15_WP          !< Minimum volume  to consider for computational geometry (normalized by min_meshsize**3)
   real(WP), parameter :: surface_epsilon_factor=1.0e-15_WP          !< Minimum surface to consider for computational geometry (normalized by min_meshsize**2)
   real(WP), parameter :: iterative_distfind_tol=1.0e-12_WP          !< Tolerance for iterative plane distance finding
   
   !> AMR-capable volume fraction solver object definition
   type :: amrvfs
      
      ! This is our amrconfig
      class(amrconfig), pointer :: amr                               !< This is the amrconfig the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_AMRVFS'             !< Solver name (default=UNNAMED_AMRVFS)
      
      ! Phasic moment and phasic interface data
      integer :: nover                                               !< Number of overlap cells
      integer :: ncomp_volmom                                        !< Number of components for moments
      integer :: ncomp_interf                                        !< Number of components for interface
      type(amrex_multifab), dimension(:), allocatable :: volmom      !< Moments multifab array
      type(amrex_multifab), dimension(:), allocatable :: interf      !< Interface multifab array
      
      ! Interface handling methods
      integer :: reconstruction_method                               !< Interface reconstruction method
      integer :: transport_method                                    !< Interface transport method
      logical :: cons_correct=.true.                                 !< Conservative correction (true by default)
      
      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), allocatable :: lo_bc_volmom           !< Boundary condition descriptor in minus direction
      integer, dimension(:,:), allocatable :: hi_bc_volmom           !< Boundary condition descriptor in plus  direction
      integer, dimension(:,:), allocatable :: lo_bc_interf           !< Boundary condition descriptor in minus direction
      integer, dimension(:,:), allocatable :: hi_bc_interf           !< Boundary condition descriptor in plus  direction
      
      ! Monitoring quantities
      real(WP) :: VFmax,VFmin,VFint                                  !< Maximum, minimum, and integral volume fraction and surface density
      
   contains
      ! Basic procedures -------------------------------------------------------------------------------------
      procedure :: initialize        !< Initialize amrvfs object
      procedure :: finalize          !< Finalize amrvfs object
      procedure :: get_info          !< Calculate various information on amrvfs object
      !procedure :: print             !< Output solver info to the screen
      ! Regrid procedures ------------------------------------------------------------------------------------
      procedure :: delete            !< Delete data at level (lvl)
      procedure :: create            !< Create data at level (lvl) and leave it uninitialized
      procedure :: refine            !< Refine data at level (lvl) using cfill procedure
      procedure :: remake            !< Remake data at level (lvl) using  fill procedure
      ! Fill procedures --------------------------------------------------------------------------------------
      procedure ::  cfill_volmom     !< Fill provided mfab at level (lvl) from this%volmom at level (lvl-1)           - done at single time - involves boundary conditions
      procedure ::   fill_volmom     !< Fill provided mfab at level (lvl) from this%volmom at level (lvl-1) and (lvl) - done at single time - involves boundary conditions
      procedure ::  cfill_interf     !< Fill provided mfab at level (lvl) from this%interf at level (lvl-1)           - done at single time - involves boundary conditions
      procedure ::   fill_interf     !< Fill provided mfab at level (lvl) from this%interf at level (lvl-1) and (lvl) - done at single time - involves boundary conditions
      ! Time advancement procedures --------------------------------------------------------------------------
      procedure :: advance           !< Advance volume moments given an interface
      procedure :: remap_moments     !< Perform semi-Lagrangian full-cell remap given an interface
      procedure :: build_plicnet     !< Reconstruct interface given volume moments using plicnet
   end type amrvfs
   
   
contains
   
   
   !> Initialization for volume fraction solver
   subroutine initialize(this,amr,reconstruction_method,transport_method,name)
      use messager,         only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_bc_ext_dir_cc,amrex_interp_cell_cons
      implicit none
      class(amrvfs), intent(inout) :: this
      class(amrconfig), target, intent(in) :: amr
      integer, intent(in) :: reconstruction_method
      integer, intent(in) :: transport_method
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to amrconfig object
      this%amr=>amr
      
      ! Number of overlap cells - CFL<1 should require only 1
      this%nover=1
      
      ! Number of components in volmom and interf multifab
      this%ncomp_volmom=7 ! VF,liquid barycenter, gas barycenter
      this%ncomp_interf=4 ! Normal and distance
      
      ! Allocate variables
      allocate(this%volmom(0:this%amr%nlvl))
      allocate(this%interf(0:this%amr%nlvl))
      
      ! Set transport scheme
      select case (transport_method)
      case (remap)
         this%transport_method=transport_method
      case default
         call die('[amrvfs initialize] Unknown transport method.')
      end select
      
      ! Set reconstruction method
      select case (reconstruction_method)
      case (plicnet)
         this%reconstruction_method=reconstruction_method
      case default
         call die('[amrvfs initialize] Unknown interface reconstruction scheme.')
      end select
      
      ! Initialize info storage
      this%VFmin=+huge(1.0_WP)
      this%VFmax=-huge(1.0_WP)
      this%VFint= 0.0_WP
      
      ! Allocate boundary condition descriptor and initialize based on periodicity info
      allocate(this%lo_bc_volmom(1:3,1:this%ncomp_volmom))
      allocate(this%hi_bc_volmom(1:3,1:this%ncomp_volmom))
      this%lo_bc_volmom=amrex_bc_ext_dir_cc
      if (this%amr%xper) this%lo_bc_volmom(1,:)=amrex_bc_int_dir
      if (this%amr%yper) this%lo_bc_volmom(2,:)=amrex_bc_int_dir
      if (this%amr%zper) this%lo_bc_volmom(3,:)=amrex_bc_int_dir
      this%hi_bc_volmom=this%lo_bc_volmom
      allocate(this%lo_bc_interf(1:3,1:this%ncomp_interf))
      allocate(this%hi_bc_interf(1:3,1:this%ncomp_interf))
      this%lo_bc_interf=amrex_bc_ext_dir_cc
      if (this%amr%xper) this%lo_bc_interf(1,:)=amrex_bc_int_dir
      if (this%amr%yper) this%lo_bc_interf(2,:)=amrex_bc_int_dir
      if (this%amr%zper) this%lo_bc_interf(3,:)=amrex_bc_int_dir
      this%hi_bc_interf=this%lo_bc_interf
      
      ! Initialize IRL
      initialize_irl: block
         use irl_fortran_interface, only: setVFBounds,setVFTolerance_IterativeDistanceFinding,&
         &                                setMinimumVolToTrack,setMinimumSAToTrack,getMoments_setMethod
         ! Transfer small constants to IRL
         call setVFBounds(VFlo)
         call setVFTolerance_IterativeDistanceFinding(iterative_distfind_tol)
         call setMinimumVolToTrack(volume_epsilon_factor*this%amr%min_meshsize**3)
         call setMinimumSAToTrack(surface_epsilon_factor*this%amr%min_meshsize**2)
         ! Set IRL's moment calculation method
         call getMoments_setMethod(half_edge)
      end block initialize_irl
      
   end subroutine initialize
   
   
   !> Finalization for amrvfs solver
   impure elemental subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy,amrex_fluxregister_destroy
      implicit none
      class(amrvfs), intent(inout) :: this
      integer :: lvl
      do lvl=0,this%amr%nlvl
         call amrex_multifab_destroy(this%volmom(lvl))
         call amrex_multifab_destroy(this%interf(lvl))
      end do
      deallocate(this%interf,this%volmom)
      nullify(this%amr)
   end subroutine finalize
   
   
   !> Calculate various information on our amrvfs object
   subroutine get_info(this)
      implicit none
      class(amrvfs), intent(inout) :: this
      
      ! Reset info
      this%VFmin=+huge(1.0_WP)
      this%VFmax=-huge(1.0_WP)
      this%VFint= 0.0_WP
      
      ! Get min and max on the finest level
      this%VFmin=min(this%VFmin,this%volmom(this%amr%clvl())%min(comp=1))
      this%VFmax=max(this%VFmax,this%volmom(this%amr%clvl())%max(comp=1))

      ! Get int at level 0
      this%VFint=this%volmom(0)%sum(comp=1)*(this%amr%dx(0)*this%amr%dy(0)*this%amr%dz(0))/this%amr%vol
      
   end subroutine get_info


   !> Delete solver data at level lvl
   subroutine delete(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrvfs), intent(inout) :: this
      integer, intent(in) :: lvl
      call amrex_multifab_destroy(this%volmom(lvl))
      call amrex_multifab_destroy(this%interf(lvl))
   end subroutine delete
   
   
   !> Create solver data at level lvl
   subroutine create(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap
      implicit none
      class(amrvfs), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Delete level
      call this%delete(lvl)
      ! Rebuild level
      call amrex_multifab_build(this%volmom(lvl),ba,dm,this%ncomp_volmom,this%nover)
      call amrex_multifab_build(this%interf(lvl),ba,dm,this%ncomp_interf,this%nover)
   end subroutine create
   
   
   !> Refine solver data at level lvl
   subroutine refine(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrvfs), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Recreate level
      call this%create(lvl,time,ba,dm)
      ! Populate volmom from coarse level
      call this%cfill_volmom(lvl,time,this%volmom(lvl))
      ! Populate interf from coarse level
      call this%cfill_interf(lvl,time,this%interf(lvl))
   end subroutine refine
   
   
   !> Remake solver data at level lvl
   subroutine remake(this,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
      implicit none
      class(amrvfs), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab) :: volmom_new,interf_new
      ! Create volmom_new and populate it from current data
      call amrex_multifab_build(volmom_new,ba,dm,this%ncomp_volmom,this%nover)
      call this%fill_volmom(lvl,time,volmom_new)
      ! Create interf_new and populate it from current data
      call amrex_multifab_build(interf_new,ba,dm,this%ncomp_interf,this%nover)
      call this%fill_interf(lvl,time,interf_new)
      ! Recreate level
      call this%create(lvl,time,ba,dm)
      ! Copy volmom_new to volmom
      call this%volmom(lvl)%copy(volmom_new,1,1,this%ncomp_volmom,this%nover)
      ! Copy interf_new to interf
      call this%interf(lvl)%copy(interf_new,1,1,this%ncomp_interf,this%nover)
      ! Destroy volmom_new
      call amrex_multifab_destroy(volmom_new)
      ! Destroy interf_new
      call amrex_multifab_destroy(interf_new)
      ! Go back and fix interf array
      
   end subroutine remake
   
   
   !> Fill provided mfab at level (lvl) from this%volmom at level (lvl-1)
   subroutine cfill_volmom(this,lvl,time,volmom)
      use amrex_amr_module, only: amrex_multifab,amrex_fillcoarsepatch,amrex_interp_cell_cons
      implicit none
      class(amrvfs), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: volmom
      real(WP) :: old_time
      ! Set small old_time to force single-time fill
      old_time=time-1.0e200_WP
      call amrex_fillcoarsepatch(      volmom,&  !< fine data being filled...
      &           old_time,this%volmom(lvl-1),&  !< using coarse data at old time...
      &               time,this%volmom(lvl-1),&  !<   and coarse data at new time...
      &           this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
      &           this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
      &            time,1,1,this%ncomp_volmom,&  !< time when we want the data, scomp, dcomp, ncomp...
      &                  this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
      &                amrex_interp_cell_cons,&  !< interpolation strategy...
      &   this%lo_bc_volmom,this%hi_bc_volmom)   !< domain bconds
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               ! AMReX default BC tool
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc_volmom,this%hi_bc_volmom)
               
               ! We should place here our user-provided BCs
               
               ! Note also that we need to fix barycenters to account for periodicity...
               ! can it be done here, or should that be done later?
               ! Depends if periodicity has been handled yet or not.
            end if
         end do
      end subroutine fillbc
   end subroutine cfill_volmom
   
   
   !> Fill provided mfab at level (lvl) from this%volmom at level (lvl-1) and (lvl)
   subroutine fill_volmom(this,lvl,time,volmom)
      use messager,         only: die
      use amrex_amr_module, only: amrex_multifab,amrex_fillpatch,amrex_interp_cell_cons
      implicit none
      class(amrvfs), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: volmom
      real(WP) :: old_time
      ! Set small old_time to force single-time fill
      old_time=time-1.0e200_WP
      ! Fill mfab
      if (lvl.eq.0) then
         ! Fill without interpolation, just direct copy and bconds
         call amrex_fillpatch(          volmom,&  !< base data being filled...
         &           old_time,this%volmom(lvl),&  !< using base data at old time...
         &               time,this%volmom(lvl),&  !<   and base data at new time...
         &           this%amr%geom(lvl),fillbc,&  !< base geometry with function to apply bconds...
         &          time,1,1,this%ncomp_volmom)   !< time when we want the data, scomp, dcomp, ncomp
         ! Unclear why lo_bc and hi_bc aren't involved here...
      else
         ! Fill with a mix of interpolation, direct copy and bconds
         call amrex_fillpatch(          volmom,&  !< fine data being filled...
         &         old_time,this%volmom(lvl-1),&  !< using coarse data at old time...
         &             time,this%volmom(lvl-1),&  !<   and coarse data at new time...
         &         this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
         &         old_time,this%volmom(lvl  ),&  !<     and fine data at old time...
         &             time,this%volmom(lvl  ),&  !<     and fine data at new time...
         &         this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
         &          time,1,1,this%ncomp_volmom,&  !< time when we want the data, scomp, dcomp, ncomp...
         &                this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
         &              amrex_interp_cell_cons,&  !< interpolation strategy...
         & this%lo_bc_volmom,this%hi_bc_volmom)   !< domain bconds
      end if
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               ! AMReX default BC tool
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc_volmom,this%hi_bc_volmom)

               ! We should place here our user-provided BCs
               
               ! Note also that we need to fix barycenters to account for periodicity...
               ! can it be done here, or should that be done later?
               ! Depends if periodicity has been handled yet or not.
            end if
         end do
      end subroutine fillbc
   end subroutine fill_volmom
   
   
   !> Fill provided mfab at level (lvl) from this%interf at level (lvl-1)
   subroutine cfill_interf(this,lvl,time,interf)
      use amrex_amr_module, only: amrex_multifab,amrex_fillcoarsepatch,amrex_interp_cell_cons
      implicit none
      class(amrvfs), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: interf
      real(WP) :: old_time
      ! Set small old_time to force single-time fill
      old_time=time-1.0e200_WP
      call amrex_fillcoarsepatch(      interf,&  !< fine data being filled...
      &           old_time,this%interf(lvl-1),&  !< using coarse data at old time...
      &               time,this%interf(lvl-1),&  !<   and coarse data at new time...
      &           this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
      &           this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
      &            time,1,1,this%ncomp_interf,&  !< time when we want the data, scomp, dcomp, ncomp...
      &                  this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
      &                amrex_interp_cell_cons,&  !< interpolation strategy...
      &   this%lo_bc_interf,this%hi_bc_interf)   !< domain bconds
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               ! AMReX default BC tool
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc_interf,this%hi_bc_interf)
               
               ! We should place here our user-provided BCs
               
               ! Note also that we need to fix barycenters to account for periodicity...
               ! can it be done here, or should that be done later?
               ! Depends if periodicity has been handled yet or not.
            end if
         end do
      end subroutine fillbc
   end subroutine cfill_interf
   
   
   !> Fill provided mfab at level (lvl) from this%interf at level (lvl-1) and (lvl)
   subroutine fill_interf(this,lvl,time,interf)
      use messager,         only: die
      use amrex_amr_module, only: amrex_multifab,amrex_fillpatch,amrex_interp_cell_cons
      implicit none
      class(amrvfs), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: interf
      real(WP) :: old_time
      ! Set small old_time to force single-time fill
      old_time=time-1.0e200_WP
      ! Fill mfab
      if (lvl.eq.0) then
         ! Fill without interpolation, just direct copy and bconds
         call amrex_fillpatch(          interf,&  !< base data being filled...
         &           old_time,this%interf(lvl),&  !< using base data at old time...
         &               time,this%interf(lvl),&  !<   and base data at new time...
         &           this%amr%geom(lvl),fillbc,&  !< base geometry with function to apply bconds...
         &          time,1,1,this%ncomp_interf)   !< time when we want the data, scomp, dcomp, ncomp
         ! Unclear why lo_bc and hi_bc aren't involved here...
      else
         ! Fill with a mix of interpolation, direct copy and bconds
         call amrex_fillpatch(          interf,&  !< fine data being filled...
         &         old_time,this%interf(lvl-1),&  !< using coarse data at old time...
         &             time,this%interf(lvl-1),&  !<   and coarse data at new time...
         &         this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
         &         old_time,this%interf(lvl  ),&  !<     and fine data at old time...
         &             time,this%interf(lvl  ),&  !<     and fine data at new time...
         &         this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
         &          time,1,1,this%ncomp_interf,&  !< time when we want the data, scomp, dcomp, ncomp...
         &                this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
         &              amrex_interp_cell_cons,&  !< interpolation strategy...
         & this%lo_bc_interf,this%hi_bc_interf)   !< domain bconds
      end if
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr), value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP), value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter) :: mfi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
         integer, dimension(4) :: plo,phi
         ! Skip if fully periodic
         if (all([this%amr%xper,this%amr%yper,this%amr%zper])) return
         ! Convert pointers
         geom=pgeom; mf=pmf
         ! Loop over boxes
         call amrex_mfiter_build(mfi,mf)
         do while(mfi%next())
            p=>mf%dataptr(mfi)
            ! Check if part of box is outside the domain
            if (.not.geom%domain%contains(p)) then
               plo=lbound(p); phi=ubound(p)
               ! AMReX default BC tool
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc_interf,this%hi_bc_interf)

               ! We should place here our user-provided BCs
               
               ! Note also that we need to fix barycenters to account for periodicity...
               ! can it be done here, or should that be done later?
               ! Depends if periodicity has been handled yet or not.
            end if
         end do
      end subroutine fillbc
   end subroutine fill_interf
   
   
   !> Calculate the new VF based on current interface and the provided dt and U/V/W
   !! U needs to be ncomp=1,nover=2,atface=[.true. ,.false.,.false.]
   !! V needs to be ncomp=1,nover=2,atface=[.false.,.true. ,.false.]
   !! W needs to be ncomp=1,nover=2,atface=[.false.,.false.,.true. ]
   subroutine advance(this,lvl,time,dt,U,V,W)
      implicit none
      class(amrvfs), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time              !< Time at which time-dependent boundary conditions are enforced
      real(WP), intent(in) :: dt                !< Timestep size over which to advance
      type(amrex_multifab), intent(in) :: U,V,W !< Velocity to use for transport
      
      ! First perform transport
      select case (this%transport_method)
      case (remap)
         call this%remap_moments(lvl=lvl,time=time,dt=dt,U=U,V=V,W=W)
      end select
      
      ! Update the band
      !call this%update_band()
      
      ! Perform interface reconstruction from transported moments
      !select case (this%reconstruction_method)
      !case (plicnet)
      !   call this%build_plicnet(lvl=lvl,time=time)
      !end select
      
      ! Create discontinuous polygon mesh from IRL interface
      !call this%polygonalize_interface()
      
      ! Reset moments to guarantee compatibility with interface reconstruction
      !call this%reset_volume_moments()
      
   end subroutine advance
   
   
   !> Perform cell-based transport of PI based on U/V/W and dt
   subroutine remap_moments(this,lvl,time,dt,U,V,W)
      use irl_fortran_interface, only: Poly24_type,CapDod_type,SepVM_type,getPt,&
      &                                adjustCapToMatchVolume,construct,new,&
      &                                PlanarLoc_type,PlanarSep_type,LocSepLink_type,&
      &                                LocLink_type,getMoments,getVolumePtr,getCentroidPtr,&
      &                                setNumberOfPlanes,setPlane,getBoundingPts
      use amrex_amr_module,      only: amrex_mfiter,amrex_box
      implicit none
      class(amrvfs), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      real(WP), intent(in) :: dt
      type(amrex_multifab), intent(in) :: U,V,W
      ! Transport data
      real(WP) :: lvol,gvol
      real(WP), dimension(3,14) :: cell
      real(WP), dimension(3, 9) :: face
      type(Poly24_type)  :: remap_cell
      type(CapDod_type)  :: remap_face
      type( SepVM_type)  :: my_SepVM
      type(amrex_mfiter) :: mfi
      type(amrex_box)    :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: volmom,interf
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      integer :: i,j,k
      ! IRL data
      type(PlanarLoc_type),  dimension(:,:,:), allocatable :: localizer
      type(PlanarSep_type),  dimension(:,:,:), allocatable :: liquid_gas_interface
      type(LocSepLink_type), dimension(:,:,:), allocatable :: localized_separator_link
      type(LocLink_type),    dimension(:,:,:), allocatable :: localizer_link
      
      ! Allocate poly24, capdod, and SepVM objects (proper destrutors called when out of scope)
      call new(remap_cell)
      call new(remap_face)
      call new(my_SepVM)
      
      ! Loop over boxes - no tiling for now
      call this%amr%mfiter_build(lvl,mfi); do while(mfi%next())
         bx=mfi%tilebox()
         
         ! Get volmom and interf data
         volmom=>this%volmom(lvl)%dataptr(mfi)
         interf=>this%interf(lvl)%dataptr(mfi)
         
         ! Get velocity data
         pU=>U%dataptr(mfi)
         pV=>V%dataptr(mfi)
         pW=>W%dataptr(mfi)
         
         ! Outer block for ensuring proper destruction of IRL object servers
         irl_object_server: block
            use irl_fortran_interface, only: ObjServer_PlanarSep_type,ObjServer_PlanarLoc_type,&
            &                                ObjServer_LocSepLink_type,ObjServer_LocLink_type,&
            &                                IRL_LargeOffsetIndex_t,new
            type(ObjServer_PlanarSep_type)  :: allocation_planar_separator
            type(ObjServer_PlanarLoc_type)  :: allocation_planar_localizer
            type(ObjServer_LocSepLink_type) :: allocation_localized_separator_link
            type(ObjServer_LocLink_type)    :: allocation_localizer_link
            integer :: nxo_,nyo_,nzo_
            integer(IRL_LargeOffsetIndex_t) :: total_cells

            ! Calculate number of cells
            nxo_=ubound(interf,1)-lbound(interf,1)+1
            nyo_=ubound(interf,2)-lbound(interf,2)+1
            nzo_=ubound(interf,3)-lbound(interf,3)+1
            total_cells=int(nxo_,8)*int(nyo_,8)*int(nzo_,8)
            
            ! Initialize size for IRL object servers
            call new(allocation_planar_localizer,total_cells)
            call new(allocation_planar_separator,total_cells)
            call new(allocation_localized_separator_link,total_cells)
            call new(allocation_localizer_link,total_cells)
            
            ! First setup IRL on the box
            setup_irl_on_box: block
               use irl_fortran_interface, only: setFromRectangularCuboid,setId,setEdgeConnectivity
               integer :: lexico
               ! Allocate IRL arrays on box
               allocate(localizer               (lbound(interf,1):ubound(interf,1),lbound(interf,2):ubound(interf,2),lbound(interf,3):ubound(interf,3)))
               allocate(liquid_gas_interface    (lbound(interf,1):ubound(interf,1),lbound(interf,2):ubound(interf,2),lbound(interf,3):ubound(interf,3)))
               allocate(localized_separator_link(lbound(interf,1):ubound(interf,1),lbound(interf,2):ubound(interf,2),lbound(interf,3):ubound(interf,3)))
               allocate(localizer_link          (lbound(interf,1):ubound(interf,1),lbound(interf,2):ubound(interf,2),lbound(interf,3):ubound(interf,3)))
               ! Initialize arrays and setup linking
               do k=lbound(interf,3),ubound(interf,3)
                  do j=lbound(interf,2),ubound(interf,2)
                     do i=lbound(interf,1),ubound(interf,1)
                        ! Transfer cell to IRL
                        call new(localizer(i,j,k),allocation_planar_localizer)
                        call setFromRectangularCuboid(localizer(i,j,k),&
                        & [this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)],&
                        & [this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)])
                        ! PLIC interface
                        call new(liquid_gas_interface(i,j,k),allocation_planar_separator)
                        ! PLIC+mesh with connectivity (i.e., link)
                        call new(localized_separator_link(i,j,k),allocation_localized_separator_link,localizer(i,j,k),liquid_gas_interface(i,j,k))
                        ! Mesh with connectivity
                        call new(localizer_link(i,j,k),allocation_localizer_link,localizer(i,j,k))
                     end do
                  end do
               end do
               ! Give each link a unique lexicographic tag (per processor)
               do k=lbound(interf,3),ubound(interf,3)
                  do j=lbound(interf,2),ubound(interf,2)
                     do i=lbound(interf,1),ubound(interf,1)
                        lexico=(i-bx%lo(1)+this%nover)+(j-bx%lo(2)+this%nover)*nxo_+(k-bx%lo(3)+this%nover)*nxo_*nyo_
                        call setId(localized_separator_link(i,j,k),lexico)
                        call setId(localizer_link(i,j,k),lexico)
                     end do
                  end do
               end do
               ! Set the connectivity
               do k=lbound(interf,3),ubound(interf,3)
                  do j=lbound(interf,2),ubound(interf,2)
                     do i=lbound(interf,1),ubound(interf,1)
                        ! In the x- direction
                        if (i.gt.lbound(interf,1)) then
                           call setEdgeConnectivity(localized_separator_link(i,j,k),0,localized_separator_link(i-1,j,k))
                           call setEdgeConnectivity(localizer_link(i,j,k),0,localizer_link(i-1,j,k))
                        end if
                        ! In the x+ direction
                        if (i.lt.ubound(interf,1)) then
                           call setEdgeConnectivity(localized_separator_link(i,j,k),1,localized_separator_link(i+1,j,k))
                           call setEdgeConnectivity(localizer_link(i,j,k),1,localizer_link(i+1,j,k))
                        end if
                        ! In the y- direction
                        if (j.gt.lbound(interf,2)) then
                           call setEdgeConnectivity(localized_separator_link(i,j,k),2,localized_separator_link(i,j-1,k))
                           call setEdgeConnectivity(localizer_link(i,j,k),2,localizer_link(i,j-1,k))
                        end if
                        ! In the y+ direction
                        if (j.lt.ubound(interf,2)) then
                           call setEdgeConnectivity(localized_separator_link(i,j,k),3,localized_separator_link(i,j+1,k))
                           call setEdgeConnectivity(localizer_link(i,j,k),3,localizer_link(i,j+1,k))
                        end if
                        ! In the z- direction
                        if (k.gt.lbound(interf,3)) then
                           call setEdgeConnectivity(localized_separator_link(i,j,k),4,localized_separator_link(i,j,k-1))
                           call setEdgeConnectivity(localizer_link(i,j,k),4,localizer_link(i,j,k-1))
                        end if
                        ! In the z+ direction
                        if (k.lt.ubound(interf,3)) then
                           call setEdgeConnectivity(localized_separator_link(i,j,k),5,localized_separator_link(i,j,k+1))
                           call setEdgeConnectivity(localizer_link(i,j,k),5,localizer_link(i,j,k+1))
                        end if
                     end do
                  end do
               end do
            end block setup_irl_on_box
            
            ! Transfer interface to IRL
            do k=lbound(interf,3),ubound(interf,3)
               do j=lbound(interf,2),ubound(interf,2)
                  do i=lbound(interf,1),ubound(interf,1)
                     call setNumberOfPlanes(liquid_gas_interface(i,j,k),1)
                     call setPlane(liquid_gas_interface(i,j,k),0,interf(i,j,k,1:3),interf(i,j,k,4))
                  end do
               end do
            end do
            
            ! Perform cell remap
            do k=bx%lo(3),bx%hi(3)
               do j=bx%lo(2),bx%hi(2)
                  do i=bx%lo(1),bx%hi(1)

                     ! Skip cells in a fully liquid or gas neighborhood
                     if (minval(volmom(i-1:i+1,j-1:j+1,k-1:k+1,1)).eq.maxval(volmom(i-1:i+1,j-1:j+1,k-1:k+1,1))) cycle
                     
                     ! Construct the cell and project it backwards in time
                     cell(:,1)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; cell(:,1)=project(cell(:,1),-dt)
                     cell(:,2)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; cell(:,2)=project(cell(:,2),-dt)
                     cell(:,3)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; cell(:,3)=project(cell(:,3),-dt)
                     cell(:,4)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; cell(:,4)=project(cell(:,4),-dt)
                     cell(:,5)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; cell(:,5)=project(cell(:,5),-dt)
                     cell(:,6)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; cell(:,6)=project(cell(:,6),-dt)
                     cell(:,7)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; cell(:,7)=project(cell(:,7),-dt)
                     cell(:,8)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; cell(:,8)=project(cell(:,8),-dt)
                     
                     ! Correct volume of x- face
                     face(:,1)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,5)=cell(:,5)
                     face(:,2)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,6)=cell(:,6)
                     face(:,3)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,7)=cell(:,7)
                     face(:,4)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,8)=cell(:,8)
                     face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
                     call construct(remap_face,face)
                     if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*pU(i  ,j,k,1)*this%amr%dy(lvl)*this%amr%dz(lvl))
                     cell(:,14)=getPt(remap_face,8)
                     
                     ! Correct volume of x+ face
                     face(:,1)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,5)=cell(:,1)
                     face(:,2)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,6)=cell(:,2)
                     face(:,3)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,7)=cell(:,3)
                     face(:,4)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,8)=cell(:,4)
                     face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
                     call construct(remap_face,face)
                     if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*pU(i+1,j,k,1)*this%amr%dy(lvl)*this%amr%dz(lvl))
                     cell(:, 9)=getPt(remap_face,8)
                     
                     ! Correct volume of y- face
                     face(:,1)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,5)=cell(:,2)
                     face(:,2)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,6)=cell(:,6)
                     face(:,3)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,7)=cell(:,5)
                     face(:,4)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,8)=cell(:,1)
                     face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
                     call construct(remap_face,face)
                     if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*pV(i,j  ,k,1)*this%amr%dz(lvl)*this%amr%dx(lvl))
                     cell(:,10)=getPt(remap_face,8)
                     
                     ! Correct volume of y+ face
                     face(:,1)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,5)=cell(:,3)
                     face(:,2)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,6)=cell(:,7)
                     face(:,3)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,7)=cell(:,8)
                     face(:,4)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,8)=cell(:,4)
                     face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
                     call construct(remap_face,face)
                     if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*pV(i,j+1,k,1)*this%amr%dz(lvl)*this%amr%dx(lvl))
                     cell(:,12)=getPt(remap_face,8)
                     
                     ! Correct volume of z- face
                     face(:,1)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,5)=cell(:,7)
                     face(:,2)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,6)=cell(:,6)
                     face(:,3)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,7)=cell(:,2)
                     face(:,4)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)]; face(:,8)=cell(:,3)
                     face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
                     call construct(remap_face,face)
                     if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*pW(i,j,k  ,1)*this%amr%dx(lvl)*this%amr%dy(lvl))
                     cell(:,11)=getPt(remap_face,8)
                     
                     ! Correct volume of z+ face
                     face(:,1)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,5)=cell(:,8)
                     face(:,2)=[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,6)=cell(:,5)
                     face(:,3)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,7)=cell(:,1)
                     face(:,4)=[this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl),this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl),this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)]; face(:,8)=cell(:,4)
                     face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
                     call construct(remap_face,face)
                     if (this%cons_correct) call adjustCapToMatchVolume(remap_face,dt*pW(i,j,k+1,1)*this%amr%dx(lvl)*this%amr%dy(lvl))
                     cell(:,13)=getPt(remap_face,8)
                     
                     ! Form remapped cell in IRL
                     call construct(remap_cell,cell)
                     
                     ! Need full geometric flux
                     call getMoments(remap_cell,localized_separator_link(i,j,k),my_SepVM)
                     
                     ! Compute new liquid volume fraction
                     lvol=getVolumePtr(my_SepVM,0)
                     gvol=getVolumePtr(my_SepVM,1)
                     volmom(i,j,k,1)=lvol/(lvol+gvol)
                     
                     ! Only work on higher order moments if VF is in [VFlo,VFhi]
                     if (volmom(i,j,k,1).lt.VFlo) then
                        volmom(i,j,k, 1 )=0.0_WP
                        volmom(i,j,k,2:4)=[this%amr%xlo+(real(i,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(j,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(k,WP)+0.5_WP)*this%amr%dz(lvl)]
                        volmom(i,j,k,5:7)=[this%amr%xlo+(real(i,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(j,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(k,WP)+0.5_WP)*this%amr%dz(lvl)]
                     else if (volmom(i,j,k,1).gt.VFhi) then
                        volmom(i,j,k, 1 )=1.0_WP
                        volmom(i,j,k,2:4)=[this%amr%xlo+(real(i,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(j,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(k,WP)+0.5_WP)*this%amr%dz(lvl)]
                        volmom(i,j,k,5:7)=[this%amr%xlo+(real(i,WP)+0.5_WP)*this%amr%dx(lvl),this%amr%ylo+(real(j,WP)+0.5_WP)*this%amr%dy(lvl),this%amr%zlo+(real(k,WP)+0.5_WP)*this%amr%dz(lvl)]
                     else
                        ! Get old phasic barycenters and project them forward in time
                        volmom(i,j,k,2:4)=getCentroidPtr(my_SepVM,0)/lvol; volmom(i,j,k,2:4)=project(volmom(i,j,k,2:4),dt)
                        volmom(i,j,k,5:7)=getCentroidPtr(my_SepVM,1)/gvol; volmom(i,j,k,5:7)=project(volmom(i,j,k,5:7),dt)                     
                     end if
                     
                  end do
               end do
            end do
            
            ! Destroy IRL data (should call proper destructors)
            deallocate(localizer,liquid_gas_interface,localized_separator_link,localizer_link)

         end block irl_object_server
         
      end do; call this%amr%mfiter_destroy(mfi)
      
      ! Fill VF and barycenter fields
      call this%fill_volmom(lvl=lvl,time=time,volmom=this%volmom(lvl))
      
   contains
      
      !> Project function that moves a point at (pU,pV,pW) for duration mydt
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), intent(in) :: mydt
         real(WP), dimension(3) :: p2
         ! Perform implicit RK2
         real(WP), dimension(3) :: p2old,v1
         real(WP) :: tolerance
         integer :: iter
         p2=p1
         tolerance=(1.0e-3_WP*this%amr%min_meshsize)*(1.0e-3_WP*this%amr%min_meshsize)
         do iter=1,10
            v1=get_velocity(0.5_WP*(p1+p2))
            p2old=p2
            p2=p1+mydt*v1
            if (dot_product(p2-p2old,p2-p2old).lt.tolerance) exit
         end do
      end function project

      !> Function that performs trilinear interpolation of velocity to pos
      function get_velocity(pos) result(vel)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel
         integer  :: ip,jp,kp
         real(WP) :: x0,y0,z0,wx1,wy1,wz1,wx2,wy2,wz2
         ! Interpolate U velocity ----------------------------------------------------------------------------------------------------
         ! Calculate clipped U-cell index
         x0=this%amr%xlo+0.0_WP*this%amr%dx(lvl); ip=min(max(floor((pos(1)-x0)/this%amr%dx(lvl)),lbound(pU,dim=1)),ubound(pU,dim=1)-1)
         y0=this%amr%ylo+0.5_WP*this%amr%dy(lvl); jp=min(max(floor((pos(2)-y0)/this%amr%dy(lvl)),lbound(pU,dim=2)),ubound(pU,dim=2)-1)
         z0=this%amr%zlo+0.5_WP*this%amr%dz(lvl); kp=min(max(floor((pos(3)-z0)/this%amr%dz(lvl)),lbound(pU,dim=3)),ubound(pU,dim=3)-1)
         ! Prepare clipped tri-linear interpolation coefficients
         wx1=min(max(modulo(pos(1)-x0,this%amr%dx(lvl))/this%amr%dx(lvl),0.0_WP),1.0_WP); wx2=1.0_WP-wx1
         wy1=min(max(modulo(pos(2)-y0,this%amr%dy(lvl))/this%amr%dy(lvl),0.0_WP),1.0_WP); wy2=1.0_WP-wy1
         wz1=min(max(modulo(pos(3)-z0,this%amr%dz(lvl))/this%amr%dz(lvl),0.0_WP),1.0_WP); wz2=1.0_WP-wz1
         ! Tri-linear interpolation of U
         vel(1)=wz1*(wy1*(wx1*pU(i+1,j+1,k+1,1)  + &
         &                wx2*pU(i  ,j+1,k+1,1)) + &
         &           wy2*(wx1*pU(i+1,j  ,k+1,1)  + &
         &                wx2*pU(i  ,j  ,k+1,1)))+ &
         &      wz2*(wy1*(wx1*pU(i+1,j+1,k  ,1)  + &
         &                wx2*pU(i  ,j+1,k  ,1)) + &
         &           wy2*(wx1*pU(i+1,j  ,k  ,1)  + &
         &                wx2*pU(i  ,j  ,k  ,1)))
         ! Interpolate V velocity ----------------------------------------------------------------------------------------------------
         ! Calculate clipped U-cell index
         x0=this%amr%xlo+0.5_WP*this%amr%dx(lvl); ip=min(max(floor((pos(1)-x0)/this%amr%dx(lvl)),lbound(pV,dim=1)),ubound(pV,dim=1)-1)
         y0=this%amr%ylo+0.0_WP*this%amr%dy(lvl); jp=min(max(floor((pos(2)-y0)/this%amr%dy(lvl)),lbound(pV,dim=2)),ubound(pV,dim=2)-1)
         z0=this%amr%zlo+0.5_WP*this%amr%dz(lvl); kp=min(max(floor((pos(3)-z0)/this%amr%dz(lvl)),lbound(pV,dim=3)),ubound(pV,dim=3)-1)
         ! Prepare clipped tri-linear interpolation coefficients
         wx1=min(max(modulo(pos(1)-x0,this%amr%dx(lvl))/this%amr%dx(lvl),0.0_WP),1.0_WP); wx2=1.0_WP-wx1
         wy1=min(max(modulo(pos(2)-y0,this%amr%dy(lvl))/this%amr%dy(lvl),0.0_WP),1.0_WP); wy2=1.0_WP-wy1
         wz1=min(max(modulo(pos(3)-z0,this%amr%dz(lvl))/this%amr%dz(lvl),0.0_WP),1.0_WP); wz2=1.0_WP-wz1
         ! Tri-linear interpolation of V
         vel(2)=wz1*(wy1*(wx1*pV(i+1,j+1,k+1,1)  + &
         &                wx2*pV(i  ,j+1,k+1,1)) + &
         &           wy2*(wx1*pV(i+1,j  ,k+1,1)  + &
         &                wx2*pV(i  ,j  ,k+1,1)))+ &
         &      wz2*(wy1*(wx1*pV(i+1,j+1,k  ,1)  + &
         &                wx2*pV(i  ,j+1,k  ,1)) + &
         &           wy2*(wx1*pV(i+1,j  ,k  ,1)  + &
         &                wx2*pV(i  ,j  ,k  ,1)))
         ! Interpolate W velocity ----------------------------------------------------------------------------------------------------
         ! Calculate clipped U-cell index
         x0=this%amr%xlo+0.5_WP*this%amr%dx(lvl); ip=min(max(floor((pos(1)-x0)/this%amr%dx(lvl)),lbound(pW,dim=1)),ubound(pW,dim=1)-1)
         y0=this%amr%ylo+0.5_WP*this%amr%dy(lvl); jp=min(max(floor((pos(2)-y0)/this%amr%dy(lvl)),lbound(pW,dim=2)),ubound(pW,dim=2)-1)
         z0=this%amr%zlo+0.0_WP*this%amr%dz(lvl); kp=min(max(floor((pos(3)-z0)/this%amr%dz(lvl)),lbound(pW,dim=3)),ubound(pW,dim=3)-1)
         ! Prepare clipped tri-linear interpolation coefficients
         wx1=min(max(modulo(pos(1)-x0,this%amr%dx(lvl))/this%amr%dx(lvl),0.0_WP),1.0_WP); wx2=1.0_WP-wx1
         wy1=min(max(modulo(pos(2)-y0,this%amr%dy(lvl))/this%amr%dy(lvl),0.0_WP),1.0_WP); wy2=1.0_WP-wy1
         wz1=min(max(modulo(pos(3)-z0,this%amr%dz(lvl))/this%amr%dz(lvl),0.0_WP),1.0_WP); wz2=1.0_WP-wz1
         ! Tri-linear interpolation of W
         vel(3)=wz1*(wy1*(wx1*pW(i+1,j+1,k+1,1)  + &
         &                wx2*pW(i  ,j+1,k+1,1)) + &
         &           wy2*(wx1*pW(i+1,j  ,k+1,1)  + &
         &                wx2*pW(i  ,j  ,k+1,1)))+ &
         &      wz2*(wy1*(wx1*pW(i+1,j+1,k  ,1)  + &
         &                wx2*pW(i  ,j+1,k  ,1)) + &
         &           wy2*(wx1*pW(i+1,j  ,k  ,1)  + &
         &                wx2*pW(i  ,j  ,k  ,1)))
      end function get_velocity
      
   end subroutine remap_moments
   

   !> Perform interface reconstruction at lvl and time using plicnet
   subroutine build_plicnet(this,lvl,time)
      use irl_fortran_interface, only: PlanarSep_type,RectCub_type,construct_2pt,setPlane,&
      &                                setNumberOfPlanes,matchVolumeFraction,getPlane,new
      use amrex_amr_module,      only: amrex_mfiter,amrex_box
      use mathtools,             only: normalize
      use plicnet,               only: get_normal,reflect_moments
      implicit none
      class(amrvfs), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time ! Needed for time-dependent boundary conditions
      ! AMReX data
      type(amrex_mfiter) :: mfi
      type(amrex_box)    :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: volmom,interf
      integer :: i,j,k,ii,jj,kk
      ! IRL reconstruction data
      type(RectCub_type)   :: cell
      type(PlanarSep_type) :: liquid_gas_interface
      logical  :: flip
      integer  :: direction
      real(WP) :: m000,m100,m010,m001,initial_dist,xc,yc,zc
      real(WP), dimension(0:188) :: moments
      real(WP), dimension(0:2)   :: normal,center
      
      ! Get an IRL cell and interface
      call new(cell)
      call new(liquid_gas_interface)
      
      ! Loop over boxes - no tiling for now
      call this%amr%mfiter_build(lvl,mfi); do while(mfi%next())
         bx=mfi%tilebox()
         
         ! Get volmom and interf data
         volmom=>this%volmom(lvl)%dataptr(mfi)
         interf=>this%interf(lvl)%dataptr(mfi)

         ! Generate nominal interf
         do k=lbound(interf,3),ubound(interf,3)
            do j=lbound(interf,2),ubound(interf,2)
               do i=lbound(interf,1),ubound(interf,1)
                  interf(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,sign(1.0_WP,volmom(i,j,k,1)-0.5_WP)]
               end do
            end do
         end do
         
         ! Perform interface reconstruction
         do k=bx%lo(3),bx%hi(3)
            do j=bx%lo(2),bx%hi(2)
               do i=bx%lo(1),bx%hi(1)
                  
                  ! Skip full cells
                  if (volmom(i,j,k,1).lt.VFlo.or.volmom(i,j,k,1).gt.VFhi) cycle
                  
                  ! Liquid-gas symmetry
                  flip=.false.
                  if (volmom(i,j,k,1).ge.0.5_WP) flip=.true.
                  m000=0.0_WP; m100=0.0_WP; m010=0.0_WP; m001=0.0_WP

                  ! Construct neighborhood of volume moments
                  if (flip) then
                     do kk=k-1,k+1
                        do jj=j-1,j+1
                           do ii=i-1,i+1
                              xc=this%amr%xlo+(real(ii,WP)+0.5_WP)*this%amr%dx(lvl)
                              yc=this%amr%ylo+(real(jj,WP)+0.5_WP)*this%amr%dy(lvl)
                              zc=this%amr%zlo+(real(kk,WP)+0.5_WP)*this%amr%dz(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0)=1.0_WP-volmom(ii,jj,kk,1)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(volmom(ii,jj,kk,5)-xc)/this%amr%dx(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(volmom(ii,jj,kk,6)-yc)/this%amr%dy(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(volmom(ii,jj,kk,7)-zc)/this%amr%dz(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(volmom(ii,jj,kk,2)-xc)/this%amr%dx(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(volmom(ii,jj,kk,3)-yc)/this%amr%dy(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(volmom(ii,jj,kk,4)-zc)/this%amr%dz(lvl)
                              ! Calculate geometric moments of neighborhood
                              m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                              m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                              m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                              m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           end do
                        end do
                     end do
                  else
                     do kk=k-1,k+1
                        do jj=j-1,j+1
                           do ii=i-1,i+1
                              xc=this%amr%xlo+(real(ii,WP)+0.5_WP)*this%amr%dx(lvl)
                              yc=this%amr%ylo+(real(jj,WP)+0.5_WP)*this%amr%dy(lvl)
                              zc=this%amr%zlo+(real(kk,WP)+0.5_WP)*this%amr%dz(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+0)=volmom(ii,jj,kk,1)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)=(volmom(ii,jj,kk,2)-xc)/this%amr%dx(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)=(volmom(ii,jj,kk,3)-yc)/this%amr%dy(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)=(volmom(ii,jj,kk,4)-zc)/this%amr%dz(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+4)=(volmom(ii,jj,kk,5)-xc)/this%amr%dx(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+5)=(volmom(ii,jj,kk,6)-yc)/this%amr%dy(lvl)
                              moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+6)=(volmom(ii,jj,kk,7)-zc)/this%amr%dz(lvl)
                              ! Calculate geometric moments of neighborhood
                              m000=m000+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                              m100=m100+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+1)+(ii-i))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                              m010=m010+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+2)+(jj-j))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                              m001=m001+(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))+3)+(kk-k))*(moments(7*((ii+1-i)*9+(jj+1-j)*3+(kk+1-k))))
                           end do
                        end do
                     end do
                  end if
                  
                  ! Calculate geometric center of neighborhood
                  center=[m100,m010,m001]/m000
                  
                  ! Symmetry about Cartesian planes
                  call reflect_moments(moments,center,direction)
                  
                  ! Get PLIC normal vector from neural network
                  call get_normal(moments,normal)
                  normal=normalize(normal)
                  
                  ! Rotate normal vector to original octant
                  if (direction.eq.1) then
                     normal(0)=-normal(0)
                  else if (direction.eq.2) then
                     normal(1)=-normal(1)
                  else if (direction.eq.3) then
                     normal(2)=-normal(2)
                  else if (direction.eq.4) then
                     normal(0)=-normal(0)
                     normal(1)=-normal(1)
                  else if (direction.eq.5) then
                     normal(0)=-normal(0)
                     normal(2)=-normal(2)
                  else if (direction.eq.6) then
                     normal(1)=-normal(1)
                     normal(2)=-normal(2)
                  else if (direction.eq.7) then
                     normal(0)=-normal(0)
                     normal(1)=-normal(1)
                     normal(2)=-normal(2)
                  end if
                  if (.not.flip) then
                     normal(0)=-normal(0)
                     normal(1)=-normal(1)
                     normal(2)=-normal(2)
                  end if
                  
                  ! Locate PLIC plane in cell using IRL
                  call construct_2pt(cell,[this%amr%xlo+real(i  ,WP)*this%amr%dx(lvl), &
                  &                        this%amr%ylo+real(j  ,WP)*this%amr%dy(lvl), &
                  &                        this%amr%zlo+real(k  ,WP)*this%amr%dz(lvl)],&
                  &                       [this%amr%xlo+real(i+1,WP)*this%amr%dx(lvl), &
                  &                        this%amr%ylo+real(j+1,WP)*this%amr%dy(lvl), &
                  &                        this%amr%zlo+real(k+1,WP)*this%amr%dz(lvl)])
                  initial_dist=dot_product(normal,[this%amr%xlo+(real(i,WP)+0.5_WP)*this%amr%dx(lvl),&
                  &                                this%amr%ylo+(real(j,WP)+0.5_WP)*this%amr%dy(lvl),&
                  &                                this%amr%zlo+(real(k,WP)+0.5_WP)*this%amr%dz(lvl)])
                  call setNumberOfPlanes(liquid_gas_interface,1)
                  call setPlane(liquid_gas_interface,0,normal,initial_dist)
                  call matchVolumeFraction(cell,volmom(i,j,k,1),liquid_gas_interface)
                  interf(i,j,k,:)=getPlane(liquid_gas_interface,0)
                  
               end do
            end do
         end do
         
      end do; call this%amr%mfiter_destroy(mfi)
      
      ! Fill VF and barycenter fields
      call this%fill_interf(lvl=lvl,time=time,interf=this%interf(lvl))
      
   end subroutine build_plicnet

   
end module amrvfs_class
