!> AMR core object is defined here
module amrcore_class
   use iso_c_binding
   use precision, only: WP
   use string,    only: str_medium
   use amrex_amr_module, only: amrex_multifab,amrex_imultifab,amrex_fluxregister
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrcore,amrdata,finalize_amrex
   public :: mk_lvl_stype,clr_lvl_stype,err_est_stype

   !> List of amrdata
   type :: amrdata
      type(amrdata), pointer :: next
      character(len=str_medium) :: name
      integer :: ncomp                    ! Number of components
      integer :: nover                    ! Number of overlap cells
      logical, dimension(3) :: ctr        ! Is variable centered in direction x/y/z
      type(amrex_multifab),  dimension(:), allocatable ::  data
      type(amrex_imultifab), dimension(:), allocatable :: idata
   end type amrdata

   !> Amrcore object definition based on AMReX's amrcore
   type :: amrcore
      
      ! Name of amcore
      character(len=str_medium) :: name='default_amrcore'
      
      ! Pointer to amrcore object
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

      ! Data register
      type(amrdata), pointer :: first_amrdata=>NULL()

      ! Arrays of pointers to subroutines to call when the grid is remade
      type(mk_lvl_stype_ptr),  dimension(:), allocatable :: mk_lvl_init_array
      type(mk_lvl_stype_ptr),  dimension(:), allocatable :: mk_lvl_crse_array
      type(mk_lvl_stype_ptr),  dimension(:), allocatable :: mk_lvl_remk_array
      type(clr_lvl_stype_ptr), dimension(:), allocatable :: clr_lvl_array
      procedure(err_est_stype), pointer, nopass :: err_est=>null()
      
   contains
      procedure :: initialize ! Initialization of amrcore object
      procedure :: finalize   ! Finalization of amrcore object
      
      procedure :: register   ! Register a new real(WP) dataset
      procedure :: iregister  ! Register a new integer dataset
      procedure :: get        ! Get handle for a dataset
      
      !procedure :: register_mk_lvl_init ! Add a mk_lvl_init procedure
      !procedure :: register_mk_lvl_crse ! Add a mk_lvl_crse procedure
      !procedure :: register_mk_lvl_remk ! Add a mk_lvl_remk procedure
      !procedure :: register_clr_lvl     ! Add a clr_lvl procedure
      !procedure :: register_err_est     ! Add a err_est procedure
   end type amrcore
   
   !> Interfaces for regriding subroutines
   abstract interface
      subroutine mk_lvl_stype(lvl,time,ba,dm) bind(c)
         import
         implicit none
         integer,     intent(in), value :: lvl
         real(WP),    intent(in), value :: time
         type(c_ptr), intent(in), value :: ba
         type(c_ptr), intent(in), value :: dm
      end subroutine mk_lvl_stype
      subroutine clr_lvl_stype(lvl) bind(c)
         import
         implicit none
         integer, intent(in) , value :: lvl
      end subroutine clr_lvl_stype
      subroutine err_est_stype(lvl,tags,time,tagval,clrval) bind(c)
         import
         implicit none
         integer,                intent(in), value :: lvl
         type(c_ptr),            intent(in), value :: tags
         real(WP),               intent(in), value :: time
         character(kind=c_char), intent(in), value :: tagval
         character(kind=c_char), intent(in), value :: clrval
      end subroutine err_est_stype
   end interface

   !> Derived types for function pointers
   type mk_lvl_stype_ptr
      procedure(mk_lvl_stype), pointer, nopass :: f=>null()
   end type mk_lvl_stype_ptr
   type clr_lvl_stype_ptr
      procedure(clr_lvl_stype), pointer, nopass :: f=>null()
   end type clr_lvl_stype_ptr
   
contains
   
   
   !subroutine amrex_init_virtual_functions (mk_lev_scrtch,mk_lev_crse,mk_lev_re,clr_lev,err_est)
   !   procedure(amrex_make_level_proc)  :: mk_lev_scrtch,mk_lev_crse,mk_lev_re
   !   procedure(amrex_clear_level_proc) :: clr_lev
   !   procedure(amrex_error_est_proc)   :: err_est
   !   call amrex_fi_init_virtual_functions(c_funloc(mk_lev_scrtch),c_funloc(mk_lev_crse),c_funloc(mk_lev_re),c_funloc(clr_lev),c_funloc(err_est),amrcore)
   !end subroutine amrex_init_virtual_functions
   
   
   !> Initialization of an amrcore object
   subroutine initialize(this,name)
      implicit none
      class(amrcore), intent(inout) :: this
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
               import
               implicit none
               type(c_ptr) :: core
            end subroutine amrex_fi_new_amrcore
         end interface
         call amrex_fi_new_amrcore(this%amrcore)
         if (present(name)) this%name=trim(adjustl(name))
      end block create_amrcore_obj
      
   end subroutine initialize
   
   
   !> Register a new real(WP) dataset in our amrcore object
   subroutine register(this,name,ncomp,nover,ctr)
      use messager, only: die
      implicit none
      class(amrcore), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      integer, intent(in) :: nover
      logical, dimension(3), intent(in) :: ctr
      type(amrdata), pointer :: new_data
      ! Prepare new amrdata
      allocate(new_data)
      new_data%name=trim(adjustl(name))
      if (ncomp.lt.1) call die('[amrcore register] ncomp needs to be greater or equal to 1')
      new_data%ncomp=ncomp
      if (nover.lt.0) call die('[amrcore register] nover needs to be greater or equal to 0')
      new_data%nover=nover
      new_data%ctr=ctr
      ! Allocate the multifab array but do not create yet
      allocate(new_data%data(0:this%nlvl-1))
      ! Insert it at the front of the list
      new_data%next=>this%first_amrdata
      this%first_amrdata=>new_data
   end subroutine register


   !> Register a new integer dataset in our amrcore object
   subroutine iregister(this,name,ncomp,nover,ctr)
      use messager, only: die
      implicit none
      class(amrcore), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: ncomp
      integer, intent(in) :: nover
      logical, dimension(3), intent(in) :: ctr
      type(amrdata), pointer :: new_data
      ! Prepare new amrdata
      allocate(new_data)
      new_data%name=trim(adjustl(name))
      if (ncomp.lt.1) call die('[amrcore iregister] ncomp needs to be greater or equal to 1')
      new_data%ncomp=ncomp
      if (nover.lt.0) call die('[amrcore iregister] nover needs to be greater or equal to 0')
      new_data%nover=nover
      new_data%ctr=ctr
      ! Allocate the imultifab array but do not create yet
      allocate(new_data%idata(0:this%nlvl-1))
      ! Insert it at the front of the list
      new_data%next=>this%first_amrdata
      this%first_amrdata=>new_data
   end subroutine iregister
   
   
   !> Get handle to a dataset in our amrcore object
   function get(this,name) result(my_data)
      use messager, only: die
      implicit none
      class(amrcore), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrdata), pointer :: my_data
      ! Traverse amrdata list until name is found
      my_data=>this%first_amrdata
      do while (name.ne.my_data%name)
         if (associated(my_data%next)) then
            my_data=>my_data%next
         else
            call die('[amrcore get] Dataset '//trim(name)//' not found in amrcore '//trim(this%name))
         end if
      end do
   end function get
   
   
   !> Finalization of amrcore object
   subroutine finalize(this)
      implicit none
      class(amrcore), intent(inout) :: this

      ! Destroy our amr mesh
      destroy_mesh: block
         interface
            subroutine amrex_fi_delete_amrcore(core) bind(c)
               import
               implicit none
               type(c_ptr), value :: core
            end subroutine amrex_fi_delete_amrcore
         end interface
         call amrex_fi_delete_amrcore(this%amrcore)
      end block destroy_mesh
      this%amrcore=c_null_ptr

      ! Destroy all multifabs
      destroy_data: block
         use amrex_amr_module, only: amrex_multifab_destroy,amrex_imultifab_destroy
         type(amrdata), pointer :: my_data
         integer :: n
         my_data=>this%first_amrdata
         do while (associated(my_data))
            ! Deallocate real data
            if (allocated(my_data%data)) then
               do n=0,this%nlvl-1
                  call amrex_multifab_destroy(my_data%data(n))
               end do
               deallocate(my_data%data)
            end if
            ! Deallocate integer data
            if (allocated(my_data%idata)) then
               do n=0,this%nlvl-1
                  call amrex_imultifab_destroy(my_data%idata(n))
               end do
               deallocate(my_data%idata)
            end if
            ! Continue on to next object
            my_data=>my_data%next
         end do
      end block destroy_data

   end subroutine finalize
   
   
   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      call amrex_finalize()
   end subroutine finalize_amrex
   
   
end module amrcore_class
   