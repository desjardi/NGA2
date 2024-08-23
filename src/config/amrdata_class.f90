!> Amrdata object is a wrapper for amrex's multifab
!> It is written as completely time-agnostic
module amrdata_class
   use string,           only: str_medium
   use precision,        only: WP
   use amrconfig_class,  only: amrconfig
   use amrex_amr_module, only: amrex_multifab
   implicit none
   private
   

   ! Expose type/constructor/methods
   public :: amrdata
   
   
   !> Amrdata object based on AMReX's multifab
   type :: amrdata
      ! Amrconfig on which data lives
      class(amrconfig), pointer :: amr
      ! Multifab array
      type(amrex_multifab), dimension(:), allocatable :: data
      ! Number of components and overlap cells
      integer :: ncomp,nover
      ! Location (true means at face, false means cell center)
      logical, dimension(3) :: atface
      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), allocatable :: lo_bc,hi_bc
      ! Interpolation method
      integer :: interp
      ! User-defined initialization procedure
      !procedure()
   contains
      procedure :: initialize       !< Initialize amrdata object
      procedure :: finalize         !< Finalize   amrdata object
      procedure :: delete_lvl       !< Delete    our data at level (lvl)
      procedure :: create_lvl       !< Create    our data at level (lvl) and leave it uninitialized
      !procedure ::   init_lvl       !< Init      our data at level (lvl) from udp
      procedure :: refine_lvl       !< Refine    our data at level (lvl) using cfill procedure
      procedure :: remake_lvl       !< Remake    our data at level (lvl) using  fill procedure
      procedure ::  cfill_lvl       !< Fill provided mfab at level (lvl) from our data at level (lvl-1)           - this involves boundary conditions - this can    be called in-place
      procedure ::   fill_lvl       !< Fill provided mfab at level (lvl) from our data at level (lvl-1) and (lvl) - this involves boundary conditions - this cannot be called in-place
   end type amrdata
   
   
contains
   
   
   !> Initialize an amrdata object
   subroutine initialize(this,amr,ncomp,nover,atface)
      use messager, only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_interp_cell_cons
      implicit none
      class(amrdata), intent(inout) :: this
      class(amrconfig), target, intent(in) :: amr
      integer, intent(in) :: ncomp
      integer, intent(in) :: nover
      logical, dimension(3), intent(in), optional :: atface
      ! Point to amrconfig object
      this%amr=>amr
      ! Allocate data
      allocate(this%data(0:this%amr%nlvl))
      ! Check ncomp and store
      if (ncomp.lt.1) call die('[amrdata initialize] ncomp needs to be at least 1')
      this%ncomp=ncomp
      ! Check nover and store
      if (nover.lt.0) call die('[amrdata initialize] nover needs to be at least 0')
      this%nover=nover
      ! Check if atface is present, if not default to cell-centered data
      if (present(atface)) then
         this%atface=atface
      else
         this%atface=.false.
      end if
      ! Assume periodic boundaries - user can change later
      allocate(this%lo_bc(1:3,1:this%ncomp),this%hi_bc(1:3,1:this%ncomp))
      this%lo_bc=amrex_bc_int_dir
      this%hi_bc=amrex_bc_int_dir
      ! Assume conservative interpolation - user can change later
      this%interp=amrex_interp_cell_cons
   end subroutine initialize


   !> Finalize an amrdata object
   subroutine finalize(this)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrdata), intent(inout) :: this
      integer :: n
      do n=0,this%amr%nlvl
         call amrex_multifab_destroy(this%data(n))
      end do
      deallocate(this%data)
      deallocate(this%lo_bc,this%hi_bc)
      this%amr=>NULL()
   end subroutine finalize
   
   
   !> Delete data at level lvl
   subroutine delete_lvl(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrdata), intent(inout) :: this
      integer,    intent(in), value :: lvl
      call amrex_multifab_destroy(this%data(lvl))
   end subroutine delete_lvl
   
   
   !> Create data at level lvl
   subroutine create_lvl(this,lvl,time,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(c_ptr),  intent(in), value :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      ! Delete data first
      call this%delete_lvl(lvl)
      ! Rebuild data
      ba=pba; dm=pdm; call amrex_multifab_build(this%data(lvl),ba,dm,this%ncomp,this%nover,this%atface)
   end subroutine create_lvl
   
   
   !> Refine data at level lvl
   subroutine refine_lvl(this,lvl,time,pba,pdm)
      use iso_c_binding, only: c_ptr
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(c_ptr),  intent(in), value :: pba,pdm
      ! Recreate data
      call this%create_lvl(lvl,time,pba,pdm)
      ! Fill from coarse level
      call this%cfill_lvl(lvl,time,this%data(lvl))
   end subroutine refine_lvl


   !> Remake data at level lvl
   subroutine remake_lvl(this,lvl,time,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(c_ptr),  intent(in), value :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(amrex_multifab)  :: newdata
      ! Create newdata and fill it from current data
      ba=pba; dm=pdm; call amrex_multifab_build(newdata,ba,dm,this%ncomp,this%nover,this%atface)
      call this%fill_lvl(lvl,time,newdata)
      ! Recreate level
      call this%create_lvl(lvl,time,pba,pdm)
      ! Copy newdata to data
      call this%data(lvl)%copy(newdata,1,1,this%ncomp,this%nover)
      ! Destroy newdata
      call amrex_multifab_destroy(newdata)
   end subroutine remake_lvl
   
   
   !> Fill provided mfab at level (lvl) from our data at level (lvl-1)
   !> Can be called in-place!
   subroutine cfill_lvl(this,lvl,time,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_fillcoarsepatch,amrex_interp_cell_cons
      implicit none
      class(amrdata), intent(inout) :: this
      integer,    intent(in), value :: lvl
      real(WP),   intent(in), value :: time
      type(amrex_multifab), intent(inout) :: mfab
      ! Fill with a mix of interpolation and bconds
      call amrex_fillcoarsepatch(   mfab,&  !< fine mfab being filled...
      &            time,this%data(lvl-1),&  !< using coarse data at old time...
      &            time,this%data(lvl-1),&  !<   and coarse data at new time...
      &      this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
      &      this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
      &              time,1,1,this%ncomp,&  !< time when we want the data, scomp, dcomp, ncomp...
      &             this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
      &                      this%interp,&  !< interpolation strategy...
      &            this%lo_bc,this%hi_bc)   !< domain bconds
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr),    value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP),       value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter)   :: mfi
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


   !> Fill provided mfab at level (lvl) from our data at level (lvl-1) and (lvl)
   !> Cannot be called in-place!
   subroutine fill_lvl(this,lvl,time,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_fillpatch,amrex_interp_cell_cons
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,      intent(in), value :: lvl
      real(WP),     intent(in), value :: time
      type(amrex_multifab), intent(inout) :: mfab
      if (lvl.eq.0) then
         ! Fill without interpolation, just direct copy and bconds
         call amrex_fillpatch(      mfab,&  !< base mfab being filled...
         &           time,this%data(lvl),&  !< using base data at old time...
         &           time,this%data(lvl),&  !<   and base data at new time...
         &     this%amr%geom(lvl),fillbc,&  !< base geometry with function to apply bconds...
         &           time,1,1,this%ncomp)   !< time when we want the data, scomp, dcomp, ncomp
         ! Unclear why lo_bc and hi_bc aren't involved here...
      else
         ! Fill with a mix of interpolation, direct copy and bconds
         call amrex_fillpatch(      mfab,&  !< fine mfab being filled...
         &         time,this%data(lvl-1),&  !< using coarse data at old time...
         &         time,this%data(lvl-1),&  !<   and coarse data at new time...
         &   this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
         &         time,this%data(lvl  ),&  !<     and fine data at old time...
         &         time,this%data(lvl  ),&  !<     and fine data at new time...
         &   this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
         &           time,1,1,this%ncomp,&  !< time when we want the data, scomp, dcomp, ncomp...
         &          this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
         &                   this%interp,&  !< interpolation strategy...
         &         this%lo_bc,this%hi_bc)   !< domain bconds
      end if
   contains
      subroutine fillbc(pmf,scomp,ncomp,t,pgeom) bind(c)
         use amrex_amr_module, only: amrex_filcc,amrex_geometry,amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_filcc
         use iso_c_binding,    only: c_ptr,c_int
         type(c_ptr),    value :: pmf,pgeom
         integer(c_int), value :: scomp,ncomp
         real(WP),       value :: t
         type(amrex_geometry) :: geom
         type(amrex_multifab) :: mf
         type(amrex_mfiter)   :: mfi
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

   
end module amrdata_class