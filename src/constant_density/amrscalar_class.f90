!> Incompressible scalar solver class built onto an amrcore
module amrscalar_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrconfig_class,  only: amrconfig
   use amrex_amr_module, only: amrex_multifab
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrscalar
   
   !> Constant density scalar solver object definition
   type :: amrscalar
      
      ! This is our amrconfig
      class(amrconfig), pointer :: amr                               !< This is the amrconfig the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_AMRSCALAR'          !< Solver name (default=UNNAMED_AMRSCALAR)
      
      ! Scalar variable definition
      integer :: nscalar                                             !< Number of scalars
      character(len=str_medium), dimension(:), allocatable :: SCname !< Names of scalars
      
      ! Scalar and old scalar data
      type(amrex_multifab), dimension(:), allocatable :: SC          !< SC multifab array
      type(amrex_multifab), dimension(:), allocatable :: SCold       !< Old SC multifab array
      
      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), allocatable :: lo_bc                  !< Boundary condition descriptor in minus direction
      integer, dimension(:,:), allocatable :: hi_bc                  !< Boundary condition descriptor in plus direction
      
      ! Interpolation method
      integer :: interp                                              !< Interpolation method

      ! Overlap size
      integer :: nover                                               !< Size of the overlap/ghost

      ! Monitoring quantities
      real(WP) :: SCmax,SCmin,SCint                                  !< Maximum and minimum, integral scalar
      
   contains
      procedure :: initialize       !< Initialize scalar solver
      procedure :: finalize         !< Finalize scalar solver
      procedure :: delete_lvl       !< Delete data at level (lvl)
      procedure :: create_lvl       !< Create data at level (lvl) and leave it uninitialized
      procedure :: refine_lvl       !< Refine data at level (lvl) using cfill procedure
      procedure :: remake_lvl       !< Remake data at level (lvl) using  fill procedure
      procedure ::  cfill_lvl       !< Fill provided mfab at level (lvl) from this%SC at level (lvl-1)
                                    !< Involves boundary conditions
                                    !< Done at a single time using this%SC
                                    !< Can be called in-place
      procedure ::   fill_lvl       !< Fill provided mfab at level (lvl) from this%SC at level (lvl-1) and (lvl)
                                    !< Involves boundary conditions
                                    !< Done at a single time using this%SC
                                    !< Cannot be called in-place

   end type amrscalar
   
   
contains
   
   
   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,nscalar,name)
      use messager,         only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_bc_reflect_even,amrex_interp_cell_cons
      implicit none
      class(amrscalar), intent(inout) :: this
      class(amrconfig), target, intent(in) :: amr
      integer, intent(in) :: nscalar
      character(len=*), optional :: name
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to amrconfig object
      this%amr=>amr
      
      ! Allocate variables
      allocate(this%SC   (0:this%amr%nlvl))
      allocate(this%SCold(0:this%amr%nlvl))
      
      ! Set the number of scalars
      this%nscalar=nscalar
      if (this%nscalar.le.0) call die('[amrscalar initialize] At least 1 scalar is required')
      
      ! Initialize scalar names
      allocate(this%SCname(1:this%nscalar))
      this%SCname='' ! User will set names

      ! Allocate boundary condition descriptor
      allocate(this%lo_bc(1:3,1:this%nscalar))
      allocate(this%hi_bc(1:3,1:this%nscalar))
      
      ! Initialize bc based on periodicity info - user can change later
      this%lo_bc=amrex_bc_reflect_even
      if (this%amr%xper) this%lo_bc(1,:)=amrex_bc_int_dir
      if (this%amr%yper) this%lo_bc(2,:)=amrex_bc_int_dir
      if (this%amr%zper) this%lo_bc(3,:)=amrex_bc_int_dir
      this%hi_bc=this%lo_bc

      ! Set overlap size for QUICK scheme
      this%nover=2
      
      ! Assume conservative interpolation - user can change later
      this%interp=amrex_interp_cell_cons
      
   end subroutine initialize


   !> Delete solver data at level lvl
   subroutine delete_lvl(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in), value :: lvl
      call amrex_multifab_destroy(this%SC   (lvl))
      call amrex_multifab_destroy(this%SCold(lvl))
   end subroutine delete_lvl
   
   
   !> Create solver data at level lvl
   subroutine create_lvl(this,lvl,time,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap
      implicit none
      class(amrscalar), intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(c_ptr),  intent(in), value :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      ! Delete data first
      call this%delete_lvl(lvl)
      ! Rebuild data
      ba=pba; dm=pdm
      call amrex_multifab_build(this%SC   (lvl),ba,dm,this%nscalar,this%nover)
      call amrex_multifab_build(this%SCold(lvl),ba,dm,this%nscalar,this%nover)
   end subroutine create_lvl
   
   
   !> Refine solver data at level lvl
   subroutine refine_lvl(this,lvl,time,pba,pdm)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrscalar), intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(c_ptr),  intent(in), value :: pba,pdm
      ! Recreate data
      call this%create_lvl(lvl,time,pba,pdm)
      ! Populate from from coarse level
      call this%cfill_lvl(lvl,time,this%SC(lvl))
   end subroutine refine_lvl


   !> Remake solver data at level lvl
   subroutine remake_lvl(this,lvl,time,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
      implicit none
      class(amrscalar), intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(c_ptr),  intent(in), value :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(amrex_multifab)  :: SCnew
      ! Create SCnew and populate it from current data
      ba=pba; dm=pdm
      call amrex_multifab_build(SCnew,ba,dm,this%nscalar,this%nover)
      call this%fill_lvl(lvl,time,SCnew)
      ! Recreate level
      call this%create_lvl(lvl,time,pba,pdm)
      ! Copy SCnew to SC
      call this%SC(lvl)%copy(SCnew,1,1,this%nscalar,this%nover)
      ! Destroy SCnew
      call amrex_multifab_destroy(SCnew)
   end subroutine remake_lvl


   !> Fill provided mfab at level (lvl) from this%SC at level (lvl-1)
   !> Can be called in-place!
   subroutine cfill_lvl(this,lvl,time,SC)
      use amrex_amr_module, only: amrex_multifab,amrex_fillcoarsepatch,amrex_interp_cell_cons
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: SC
      call amrex_fillcoarsepatch(          SC,&  !< fine data being filled...
      &                   time,this%SC(lvl-1),&  !< using coarse data at old time...
      &                   time,this%SC(lvl-1),&  !<   and coarse data at new time...
      &           this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
      &           this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
      &              time,1,1,this%SC(lvl)%nc,&  !< time when we want the data, scomp, dcomp, ncomp...
      &                  this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
      &                           this%interp,&  !< interpolation strategy...
      &                 this%lo_bc,this%hi_bc)   !< domain bconds
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
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc,this%hi_bc)
            end if
         end do
         ! This will need hooks for user-provided BCs
      end subroutine fillbc
   end subroutine cfill_lvl


   !> Fill provided mfab at level (lvl) from this%SC at level (lvl-1) and (lvl)
   !> Cannot be called in-place!
   subroutine fill_lvl(this,lvl,time,SC)
      use amrex_amr_module, only: amrex_multifab,amrex_fillpatch,amrex_interp_cell_cons
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_multifab), intent(inout) :: SC
      if (lvl.eq.0) then
         ! Fill without interpolation, just direct copy and bconds
         call amrex_fillpatch(          SC,&  !< base data being filled...
         &               time,this%SC(lvl),&  !< using base data at old time...
         &               time,this%SC(lvl),&  !<   and base data at new time...
         &       this%amr%geom(lvl),fillbc,&  !< base geometry with function to apply bconds...
         &        time,1,1,this%SC(lvl)%nc)   !< time when we want the data, scomp, dcomp, ncomp
         ! Unclear why lo_bc and hi_bc aren't involved here...
      else
         ! Fill with a mix of interpolation, direct copy and bconds
         call amrex_fillpatch(          SC,&  !< fine data being filled...
         &             time,this%SC(lvl-1),&  !< using coarse data at old time...
         &             time,this%SC(lvl-1),&  !<   and coarse data at new time...
         &     this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
         &             time,this%SC(lvl  ),&  !<     and fine data at old time...
         &             time,this%SC(lvl  ),&  !<     and fine data at new time...
         &     this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
         &        time,1,1,this%SC(lvl)%nc,&  !< time when we want the data, scomp, dcomp, ncomp...
         &            this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
         &                     this%interp,&  !< interpolation strategy...
         &           this%lo_bc,this%hi_bc)   !< domain bconds
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
               call amrex_filcc(p,plo,phi,geom%domain%lo,geom%domain%hi,geom%dx,geom%get_physical_location(plo),this%lo_bc,this%hi_bc)
            end if
         end do
         ! This will need hooks for user-provided BCs
      end subroutine fillbc
   end subroutine fill_lvl
   
   
   !> Finalization for amrscalar solver
   subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrscalar), intent(inout) :: this
      integer :: n
      do n=0,this%amr%nlvl
         call amrex_multifab_destroy(this%SC   (n))
         call amrex_multifab_destroy(this%SCold(n))
      end do
      deallocate(this%SC,this%SCold)
   end subroutine finalize
   
   
end module amrscalar_class
