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
      type(amrex_multifab), dimension(:), allocatable ::  SC         !< SC multifab array
      type(amrex_multifab), dimension(:), allocatable ::  SCold      !< Old SC multifab array
      
      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), allocatable :: lo_bc                  !< Boundary condition descriptor in minus direction
      integer, dimension(:,:), allocatable :: hi_bc                  !< Boundary condition descriptor in plus direction
      
      ! Overlap size
      integer :: nover                                               !< Size of the overlap/ghost

      ! Monitoring quantities
      real(WP) :: SCmax,SCmin,SCint                                  !< Maximum and minimum, integral scalar
      
   contains
      procedure :: initialize                                        !< Initialize scalar solver
      procedure :: finalize                                          !< Finalize scalar solver
      procedure :: delete_lvl                                        !< Delete data at level lvl
      procedure :: create_lvl                                        !< Create data at level lvl
      procedure :: fillcoarse                                        !< Populate patch of data from coarse data (single time)
   end type amrscalar
   
   
contains
   
   
   !> Initialization for amrscalar solver
   subroutine initialize(this,amr,nscalar,name)
      use messager,         only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_bc_reflect_even
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
      if (this%nscalar.le.0) call die('[amrscalar constructor] At least 1 scalar is required')
      
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

      ! Set overlap size
      this%nover=2
      
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
   subroutine create_lvl(this,lvl,ba,dm)
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap
      implicit none
      class(amrscalar), intent(inout) :: this
      type(amrex_boxarray),  intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      integer, intent(in), value :: lvl
      call amrex_multifab_build(this%SC   (lvl),ba,dm,this%nscalar,this%nover)
      call amrex_multifab_build(this%SCold(lvl),ba,dm,this%nscalar,this%nover)
   end subroutine create_lvl


   !> Fill a patch from coarse SC data and boundary conditions
   subroutine fillcoarse(this,lvl,time)
      use amrex_amr_module, only: amrex_fillcoarsepatch,amrex_interp_cell_cons
      implicit none
      class(amrscalar), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      call amrex_fillcoarsepatch(this%SC(lvl),&  !< fine data being filled...
      &                   time,this%SC(lvl-1),&  !< using coarse data at new time...
      &                   time,this%SC(lvl-1),&  !<   and coarse data at old time...
      &           this%amr%geom(lvl-1),fillbc,&  !< coarse geometry with function to apply bconds...
      &           this%amr%geom(lvl  ),fillbc,&  !<   fine geometry with function to apply bconds...
      &              time,1,1,this%SC(lvl)%nc,&  !< time when we want the data, scomp, dcomp, ncomp...
      &                  this%amr%rref(lvl-1),&  !< refinement ratio between the levels...
      &                amrex_interp_cell_cons,&  !< interpolation strategy...
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
      end subroutine fillbc
      
   end subroutine fillcoarse
   
   
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
