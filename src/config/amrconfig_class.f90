!> AMR config object is defined based on amrcore AMREX object
!> Amrconfig differs quite a bit from other configs
module amrconfig_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_char,c_int,c_funptr,c_funloc
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_geometry
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: amrconfig,finalize_amrex
   public :: mak_lvl_stype,clr_lvl_stype,err_est_stype
   
   
   !> Amrconfig object definition based on AMReX's amrcore
   type :: amrconfig
      ! Name of amrconfig
      character(len=str_medium) :: name='UNNAMED_AMRGRID'
      ! Pointer to AMReX's amrcore object
      type(c_ptr) :: amrcore=c_null_ptr
      ! Coordinate system
      integer :: coordsys=0
      ! Domain periodicity
      logical :: xper,yper,zper
      ! Domain extent
      real(WP) :: xlo,ylo,zlo
      real(WP) :: xhi,yhi,zhi
      ! Level 0 mesh dimensions
      integer :: nx,ny,nz
      ! Number of refinement levels
      integer :: nlvl
      ! Max grid size
      integer :: nmax=32
      ! Blocking factor
      integer :: nbloc=8
      ! Refinement ratio
      integer :: rref=2
      ! Geometry object at each level
      type(amrex_geometry), dimension(:), allocatable :: geom
   contains
      procedure :: initialize            ! Initialization of amrconfig object
      procedure :: finalize              ! Finalization of amrconfig object
      procedure :: register_udf          ! Register user-defined functions for regriding
      procedure :: initialize_data       ! Initialize data on armconfig according to registered function
   end type amrconfig
   
   
   !> Interfaces for regriding subroutines
   abstract interface
      subroutine mak_lvl_stype(lvl,time,ba,dm) bind(c)
         import :: c_ptr,WP
         implicit none
         integer,     intent(in), value :: lvl
         real(WP),    intent(in), value :: time
         type(c_ptr), intent(in), value :: ba
         type(c_ptr), intent(in), value :: dm
      end subroutine mak_lvl_stype
      subroutine clr_lvl_stype(lvl) bind(c)
         implicit none
         integer, intent(in) , value :: lvl
      end subroutine clr_lvl_stype
      subroutine err_est_stype(lvl,tags,time,tagval,clrval) bind(c)
         import :: c_ptr,c_char,WP
         implicit none
         integer,                intent(in), value :: lvl
         type(c_ptr),            intent(in), value :: tags
         real(WP),               intent(in), value :: time
         character(kind=c_char), intent(in), value :: tagval
         character(kind=c_char), intent(in), value :: clrval
      end subroutine err_est_stype
   end interface
   
   
contains
   
   
   !> Initialization of an amrconfig object
   subroutine initialize(this,name)
      implicit none
      class(amrconfig), intent(inout) :: this
      character(len=*), intent(in), optional :: name
      ! First of all, confirm that nga2 and amrex reals are compatible
      check_real: block
         use messager,         only: die
         use precision,        only: WP
         use amrex_amr_module, only: amrex_real
         if (amrex_real.ne.WP) call die('Incompatible real type between nga2 and amrex!')
      end block check_real
      ! Check if AMReX has been initialized - if not, initialize it here
      check_if_init: block
         use amrex_amr_module, only: amrex_initialized,amrex_init
         use parallel, only: comm
         ! Note that we are using our main communicator (not sure yet if AMReX can work with *multiple* communicators)
         ! We are not parsing the command line at all, and we're initializing without any parameters - may need to change...
         if (.not.amrex_initialized()) call amrex_init(comm=comm%MPI_VAL,arg_parmparse=.false.)
      end block check_if_init
      ! Set parameters for amrcore and geometry
      set_params: block
         use amrex_amr_module, only: amrex_parmparse,amrex_parmparse_build,amrex_parmparse,amrex_parmparse_destroy
         type(amrex_parmparse) :: pp
         integer, dimension(3) :: per
         call amrex_parmparse_build(pp,'amr')
         call pp%addarr('n_cell'         ,[this%nx,this%ny,this%nz])
         call pp%add   ('max_level'      ,this%nlvl)
         call pp%add   ('blocking_factor',this%nbloc)
         call pp%add   ('max_grid_size'  ,this%nmax)
         call pp%add   ('ref_ratio'      ,this%rref)
         call amrex_parmparse_destroy(pp)
         call amrex_parmparse_build(pp,'geometry')
         call pp%add   ('coord_sys'  ,this%coordsys)
         per=0
         if (this%xper) per(1)=1
         if (this%yper) per(2)=1
         if (this%zper) per(3)=1
         call pp%addarr('is_periodic',per)
         call pp%addarr('prob_lo'    ,[this%xlo,this%ylo,this%zlo])
         call pp%addarr('prob_hi'    ,[this%xhi,this%yhi,this%zhi])
         call amrex_parmparse_destroy(pp)
      end block set_params
      ! Create an amrcore object
      create_amrcore_obj: block
         interface
            subroutine amrex_fi_new_amrcore(core) bind(c)
               import :: c_ptr
               implicit none
               type(c_ptr) :: core
            end subroutine amrex_fi_new_amrcore
         end interface
         call amrex_fi_new_amrcore(this%amrcore)
         if (present(name)) this%name=trim(adjustl(name))
      end block create_amrcore_obj
      ! Get back geometry objects
      store_geometries: block
         use amrex_amr_module, only: amrex_geometry_init_data
         interface
            subroutine amrex_fi_get_geometry(geom,lvl,core) bind(c)
               import :: c_ptr,c_int
               implicit none
               type(c_ptr), intent(out) :: geom
               integer(c_int), value :: lvl
               type(c_ptr), value :: core
            end subroutine amrex_fi_get_geometry
         end interface
         integer :: n
         allocate(this%geom(0:this%nlvl-1))
         do n=0,this%nlvl-1
            call amrex_fi_get_geometry(this%geom(n)%p,n,this%amrcore)
            call amrex_geometry_init_data(this%geom(n))
         end do
      end block store_geometries
   end subroutine initialize
   
   
   !> Initialize data on an amrconfig object
   subroutine initialize_data(this,time)
      implicit none
      class(amrconfig), intent(inout) :: this
      real(WP), intent(in) :: time
      interface
         subroutine amrex_fi_init_from_scratch (t,core) bind(c)
            import :: c_ptr,WP
            implicit none
            real(WP), value :: t
            type(c_ptr), value :: core
         end subroutine amrex_fi_init_from_scratch
      end interface
      call amrex_fi_init_from_scratch(time,this%amrcore)
   end subroutine initialize_data
   
   
   !> Finalization of amrconfig object
   subroutine finalize(this)
      implicit none
      class(amrconfig), intent(inout) :: this
      interface
         subroutine amrex_fi_delete_amrcore(core) bind(c)
            import :: c_ptr
            implicit none
            type(c_ptr), value :: core
         end subroutine amrex_fi_delete_amrcore
      end interface
      call amrex_fi_delete_amrcore(this%amrcore)
      this%amrcore=c_null_ptr
   end subroutine finalize
   
   
   !> Register user-defined functions for regriding
   subroutine register_udf(this,mak_lvl_init,mak_lvl_crse,mak_lvl_remk,clr_lvl,err_est)
      implicit none
      class(amrconfig), intent(inout) :: this
      procedure(mak_lvl_stype) :: mak_lvl_init,mak_lvl_crse,mak_lvl_remk
      procedure(clr_lvl_stype) :: clr_lvl
      procedure(err_est_stype) :: err_est
      interface
         subroutine amrex_fi_init_virtual_functions(maklvl_init,maklvl_crse,maklvl_remk,clrlvl,errest,core) bind(c)
            import :: c_ptr,c_funptr
            implicit none
            type(c_funptr), value :: maklvl_init,maklvl_crse,maklvl_remk,clrlvl,errest
            type(c_ptr), value :: core
         end subroutine amrex_fi_init_virtual_functions
      end interface
      call amrex_fi_init_virtual_functions(c_funloc(mak_lvl_init),c_funloc(mak_lvl_crse),c_funloc(mak_lvl_remk),c_funloc(clr_lvl),c_funloc(err_est),this%amrcore)
   end subroutine register_udf
   
   
   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      call amrex_finalize()
   end subroutine finalize_amrex
   
   
end module amrconfig_class
   