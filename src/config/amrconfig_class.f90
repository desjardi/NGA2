!> AMR config object is defined based on amrcore AMREX object
!> Amrconfig differs quite a bit from other configs
module amrconfig_class
   use iso_c_binding,    only: c_ptr,c_null_ptr,c_char,c_int,c_funptr,c_funloc
   use precision,        only: WP
   use string,           only: str_medium
   use amrex_amr_module, only: amrex_geometry,amrex_multifab,amrex_tagboxarray
   use mpi_f08,          only: MPI_Comm
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: amrconfig,amrdata,finalize_amrex
   public :: mak_lvl_stype,clr_lvl_stype,err_est_stype
   
   
   ! Types of amrdata registration available
   integer, parameter, public :: ignore=0
   integer, parameter, public :: resize=1
   integer, parameter, public :: interp=2

   
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
      integer, dimension(:), allocatable :: rref
      ! Geometry object at each level
      type(amrex_geometry), dimension(:), allocatable :: geom
      ! Shortcut for domain volume
      real(WP) :: vol
      ! Shortcut to cell size per level
      real(WP), dimension(:), allocatable :: dx,dy,dz
      ! Data register
      type(reg_list), pointer :: first_reg=>NULL()
      ! Parallel info
      type(MPI_Comm) :: comm            !< Communicator for our group
      integer        :: nproc           !< Number of processors
      integer        :: rank            !< Processor grid rank
      logical        :: amRoot          !< Am I the root?
      ! Monitoring info - call get_info() to update
      integer  :: nlevels=-1            !< Current total number of levels
      integer  :: nboxes =-1            !< Current total number of boxes
      real(WP) :: ncells =-1.0_WP       !< Current total number of cells (real!)
      real(WP) :: compression=-1.0_WP   !< Current compression ratio (ncells/total cells with uniform mesh)
   contains
      procedure :: initialize=>initialize_amrconfig   !< Initialization of amrconfig object
      procedure :: finalize=>finalize_amrconfig       !< Finalization of amrconfig object
      procedure :: register_udp                       !< Register user-defined procedures for regriding
      procedure :: define_lvl=>amr_define_lvl         !< Create and initialize all registered data on amrconfig at that lvl
      procedure :: refine_lvl=>amr_refine_lvl         !< Refine all registered data on amrconfig at that lvl
      procedure :: remake_lvl=>amr_remake_lvl         !< Remake all registered data on amrconfig at that lvl
      procedure :: delete_lvl=>amr_delete_lvl         !< Delete all registered data on amrconfig at that lvl
      procedure :: initialize_grid                    !< Initialize data on armconfig according to registered function
      procedure :: regrid                             !< Perform regriding operation on level baselvl
      procedure :: get_info=>get_info_amrconfig       !< Calculate various information on our amrconfig object
      procedure :: print=>print_amrconfig             !< Print out grid info
      procedure, private :: get_boxarray              !< Obtain box array at a given level
      procedure, private :: get_distromap             !< Obtain distromap at a given level
      procedure :: mfiter_build                       !< Build mfiter at a given level
      procedure :: mfiter_destroy                     !< Destroy mfiter
      procedure :: clvl                               !< Return current finest level
      procedure :: average_down                       !< Average down a given multifab throughout all levels
      procedure :: average_downto                     !< Average down a given multifab to level lvl
   end type amrconfig
   

   !> Amrdata object based on array of AMReX's multifab
   type :: amrdata
      ! Amrconfig on which data lives
      class(amrconfig), pointer :: amr
      ! This is the name of the amrdata
      character(len=str_medium) :: name='UNNAMED_AMRDATA'
      ! Amrdata object registration status in amrconfig
      integer :: reg
      ! Multifab array
      type(amrex_multifab), dimension(:), allocatable :: data
      ! Number of components and overlap cells
      integer :: ncomp,nover
      ! Location (true means at face, false means cell center)
      logical, dimension(3) :: atface
      ! Boundary conditions at domain boundaries
      integer, dimension(:,:), pointer :: lo_bc,hi_bc
      logical :: bc_ptr=.false.
      ! Interpolation method
      integer :: interp
      ! User-defined initialization procedure
      !procedure()
      ! Global info on our data
      real(WP), dimension(:), allocatable :: mindata,maxdata,intdata
   contains
      procedure :: initialize=>initialize_amrdata !< Initialize amrdata object
      procedure :: finalize=>finalize_amrdata     !< Finalize   amrdata object
      procedure :: get_info=>get_info_amrdata     !< Compute min/max/int values for field
      procedure :: define_lvl                     !< Define our data at level (lvl) from udp
      procedure :: refine_lvl                     !< Refine our data at level (lvl) using cfill procedure
      procedure :: remake_lvl                     !< Remake our data at level (lvl) using  fill procedure
      procedure :: delete_lvl                     !< Delete our data at level (lvl)
      procedure :: create_lvl                     !< Create our data at level (lvl) and leave it uninitialized
      procedure ::  cfill_lvl                     !< Fill provided mfab at level (lvl) from our data at level (lvl-1)           - this involves boundary conditions
      procedure ::   fill_lvl                     !< Fill provided mfab at level (lvl) from our data at level (lvl-1) and (lvl) - this involves boundary conditions
   end type amrdata
   
   
   !> List of registered amrdata objects defined on an amrconfig
   type :: reg_list
      type(reg_list), pointer :: next         !< Next list entry
      type(amrdata) , pointer :: data=>NULL() !< Amrdata object
   end type reg_list
   
   
   !> Interfaces for user-defined subroutines
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
      subroutine reftag_lvl_stype(lvl,tags,time,tagval,clrval)
         import :: c_char,WP,amrex_tagboxarray
         implicit none
         integer,                 intent(in) :: lvl
         type(amrex_tagboxarray), intent(in) :: tags
         real(WP),                intent(in) :: time
         character(kind=c_char),  intent(in) :: tagval
         character(kind=c_char),  intent(in) :: clrval
      end subroutine reftag_lvl_stype
   end interface
   
   
contains
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! AMRCONFIG specific routines !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   !> Initialization of an amrconfig object
   subroutine initialize_amrconfig(this,name)
      use messager, only: die
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
         if (this%nlvl.lt.0) call die('[amrconfig initialize] nlvl must be >= 0')
         call pp%add   ('max_level'      ,this%nlvl)
         call pp%add   ('blocking_factor',this%nbloc)
         call pp%add   ('max_grid_size'  ,this%nmax)
         if (.not.allocated(this%rref)) this%rref=[2]
         call pp%addarr('ref_ratio'      ,this%rref)
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
         allocate(this%geom(0:this%nlvl))
         do n=0,this%nlvl
            call amrex_fi_get_geometry(this%geom(n)%p,n,this%amrcore)
            call amrex_geometry_init_data(this%geom(n))
         end do
      end block store_geometries
      ! Store effective refinement ratio
      store_ref_ratio: block
         interface
            subroutine amrex_fi_get_ref_ratio(ref_ratio,core) bind(c)
               import :: c_ptr
               implicit none
               integer, dimension(*), intent(inout) :: ref_ratio
               type(c_ptr), value :: core
            end subroutine amrex_fi_get_ref_ratio
         end interface
         deallocate(this%rref); allocate(this%rref(0:this%nlvl-1))
         call amrex_fi_get_ref_ratio(this%rref,this%amrcore)
      end block store_ref_ratio
      ! Store parallel info
      store_parallel_info: block
         use parallel, only: comm,nproc,rank,amRoot
         this%comm=comm
         this%nproc=nproc
         this%rank=rank
         this%amRoot=amRoot
      end block store_parallel_info
      ! Compute shortcuts
      compute_shortcuts: block
         integer :: lvl
         ! Mesh size at each level
         allocate(this%dx(0:this%nlvl),this%dy(0:this%nlvl),this%dz(0:this%nlvl))
         do lvl=0,this%nlvl
            this%dx(lvl)=this%geom(lvl)%dx(1)
            this%dy(lvl)=this%geom(lvl)%dx(2)
            this%dz(lvl)=this%geom(lvl)%dx(3)
         end do
         ! Total domain volume
         this%vol=(this%xhi-this%xlo)*(this%yhi-this%ylo)*(this%zhi-this%zlo)
      end block compute_shortcuts
   end subroutine initialize_amrconfig
   
   
   !> Finalization of amrconfig object
   impure elemental subroutine finalize_amrconfig(this)
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
      ! Deallocate reg list
      empty_reg_list: block
         type(reg_list), pointer :: current,next
         if (associated(this%first_reg)) then
            current=>this%first_reg
            next=>current%next
            do
               deallocate(current)
               if (.not.associated(next)) exit
               current=>next
               next=>current%next
            end do
         end if
      end block empty_reg_list
   end subroutine finalize_amrconfig
   
   
   !> Register user-defined procedures for regriding
   subroutine register_udp(this,mak_lvl_init,mak_lvl_crse,mak_lvl_remk,clr_lvl,err_est)
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
   end subroutine register_udp
   

   !> Default routine to define level data
   subroutine amr_define_lvl(this,lvl,t,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrconfig), intent(inout) :: this
      integer,     intent(in) :: lvl
      real(WP),    intent(in) :: t
      type(c_ptr), intent(in) :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(reg_list), pointer :: my_reg
      ! Recast pointers
      ba=pba; dm=pdm
      ! Traverse all registered data and define each
      my_reg=>this%first_reg
      do while (associated(my_reg))
         select case (my_reg%data%reg)
         case (resize)
            call my_reg%data%create_lvl(lvl,t,ba,dm)
         case (interp)
            call my_reg%data%define_lvl(lvl,t,ba,dm)
         end select
         my_reg=>my_reg%next
      end do
   end subroutine amr_define_lvl
   
   
   !> Default routine to refine level data
   subroutine amr_refine_lvl(this,lvl,t,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrconfig), intent(inout) :: this
      integer,     intent(in) :: lvl
      real(WP),    intent(in) :: t
      type(c_ptr), intent(in) :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(reg_list), pointer :: my_reg
      ! Recast pointers
      ba=pba; dm=pdm
      ! Traverse all registered data and refine each
      my_reg=>this%first_reg
      do while (associated(my_reg))
         select case (my_reg%data%reg)
         case (resize)
            call my_reg%data%create_lvl(lvl,t,ba,dm)
         case (interp)
            call my_reg%data%refine_lvl(lvl,t,ba,dm)
         end select
         my_reg=>my_reg%next
      end do
   end subroutine amr_refine_lvl
   
   
   !> Default routine to remake level data
   subroutine amr_remake_lvl(this,lvl,t,pba,pdm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrconfig), intent(inout) :: this
      integer,     intent(in) :: lvl
      real(WP),    intent(in) :: t
      type(c_ptr), intent(in) :: pba,pdm
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(reg_list), pointer :: my_reg
      ! Recast pointers
      ba=pba; dm=pdm
      ! Traverse all registered data and remake each
      my_reg=>this%first_reg
      do while (associated(my_reg))
         select case (my_reg%data%reg)
         case (resize)
            call my_reg%data%create_lvl(lvl,t,ba,dm)
         case (interp)
            call my_reg%data%remake_lvl(lvl,t,ba,dm)
         end select
         my_reg=>my_reg%next
      end do
   end subroutine amr_remake_lvl
   
   
   !> Default routine to delete level data
   subroutine amr_delete_lvl(this,lvl)
      implicit none
      class(amrconfig), intent(inout) :: this
      integer, intent(in) :: lvl
      type(reg_list), pointer :: my_reg
      ! Traverse all registered data and delete each
      my_reg=>this%first_reg
      do while (associated(my_reg))
         call my_reg%data%delete_lvl(lvl)
         my_reg=>my_reg%next
      end do
   end subroutine amr_delete_lvl
   
   
   !> Initialize grid on an amrconfig object
   subroutine initialize_grid(this,time)
      implicit none
      class(amrconfig), intent(inout) :: this
      real(WP), intent(in) :: time
      interface
         subroutine amrex_fi_init_from_scratch(t,core) bind(c)
            import :: c_ptr,WP
            implicit none
            real(WP), value :: t
            type(c_ptr), value :: core
         end subroutine amrex_fi_init_from_scratch
      end interface
      ! Generate grid and allocate data
      call amrex_fi_init_from_scratch(time,this%amrcore)
      ! Generate info about grid
      call this%get_info()
   end subroutine initialize_grid
   
   
   !> Perform regriding operation on baselvl
   subroutine regrid(this,baselvl,time)
      implicit none
      class(amrconfig), intent(inout) :: this
      integer,  intent(in) :: baselvl
      real(WP), intent(in) :: time
      interface
         subroutine amrex_fi_regrid(blvl,t,core) bind(c)
            import :: c_ptr,c_int,WP
            implicit none
            integer(c_int), value :: blvl
            real(WP), value :: t
            type(c_ptr), value :: core
         end subroutine amrex_fi_regrid
      end interface
      ! Regenerate grid and resize/transfer data
      call amrex_fi_regrid(baselvl,time,this%amrcore)
      ! Generate info about grid
      call this%get_info()
   end subroutine regrid
   
   
   !> Get info on amrconfig object
   subroutine get_info_amrconfig(this)
      use amrex_amr_module, only: amrex_boxarray,amrex_box
      implicit none
      class(amrconfig), intent(inout) :: this
      type(amrex_boxarray) :: ba
      type(amrex_box)      :: bx
      integer :: n,m,nb
      integer, dimension(3) :: lo,hi
      ! Update number of levels
      this%nlevels=this%clvl()+1
      ! Update number of boxes and cells
      this%nboxes=0
      this%ncells=0.0_WP
      do n=0,this%clvl()
         ! Get boxarray
         ba=this%get_boxarray(lvl=n)
         ! Increment boxes
         nb=int(ba%nboxes())
         this%nboxes=this%nboxes+nb
         ! Traverse boxes and count cells
         do m=1,nb
            ! Get box
            bx=ba%get_box(m-1)
            ! Increment cells
            this%ncells=this%ncells+real((bx%hi(1)-bx%lo(1)+1),WP)*&
            &                       real((bx%hi(2)-bx%lo(2)+1),WP)*&
            &                       real((bx%hi(3)-bx%lo(3)+1),WP)
         end do
      end do
      ! Update compression ratio
      this%compression=this%geom(this%clvl())%dx(1)*&
      &                this%geom(this%clvl())%dx(2)*&
      &                this%geom(this%clvl())%dx(3)*&
      &                this%ncells/((this%xhi-this%xlo)*(this%yhi-this%ylo)*(this%zhi-this%zlo))
   end subroutine get_info_amrconfig
   
   
   !> Print amrconfig object
   subroutine print_amrconfig(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use amrex_amr_module, only: amrex_boxarray,amrex_box
      use parallel, only: amRoot
      implicit none
      class(amrconfig), intent(inout) :: this
      type(amrex_boxarray) :: ba
      type(amrex_box)      :: bx
      integer :: n,m,nb
      integer, dimension(3) :: lo,hi
      if (amRoot) then
         write(output_unit,'("AMR Cartesian grid [",a,"]")') trim(this%name)
         write(output_unit,'(" > amr level = ",i2)') this%clvl()
         write(output_unit,'(" > max level = ",i2)') this%nlvl
         write(output_unit,'(" > ref ratio = ",100(" ",i0))') this%rref
         write(output_unit,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%xlo,this%xhi,this%ylo,this%yhi,this%zlo,this%zhi
         write(output_unit,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper
         ! Loop over levels
         do n=0,this%clvl()
            ! Get boxarray at that level
            ba=this%get_boxarray(lvl=n)
            write(output_unit,'(" >>> Level ",i2,"/",i2," with [dx =",es12.5,",dy =",es12.5,",dz =",es12.5,"]")') n,this%clvl(),this%geom(n)%dx(1),this%geom(n)%dx(2),this%geom(n)%dx(3)
            ! Loop over boxes
            nb=int(ba%nboxes())
            do m=1,nb
               bx=ba%get_box(m-1)
               write(output_unit,'(" >>> Box ",i0,"/",i0," from [",i3,",",i3,",",i3,"] to [",i3,",",i3,",",i3,"]")') m,nb,bx%lo,bx%hi
            end do
         end do
      end if
   end subroutine print_amrconfig
   
   
   !> Obtain box array at a level
   function get_boxarray(this,lvl) result(ba)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrconfig), intent(inout) :: this
      integer, intent(in)  :: lvl
      type(amrex_boxarray) :: ba
      interface
         subroutine amrex_fi_get_boxarray(barray,lev,core) bind(c)
            import :: c_ptr,c_int
            implicit none
            type(c_ptr), intent(out) :: barray
            integer(c_int), value :: lev
            type(c_ptr), value :: core
         end subroutine amrex_fi_get_boxarray
      end interface
      call amrex_fi_get_boxarray(ba%p,lvl,this%amrcore)
   end function get_boxarray
   
   
   !> Obtain distromap at a level
   function get_distromap(this,lvl) result(dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrconfig), intent(inout) :: this
      integer, intent(in)  :: lvl
      type(amrex_distromap) :: dm
      interface
         subroutine amrex_fi_get_distromap(dmap,lev,core) bind(c)
            import :: c_ptr,c_int
            implicit none
            type(c_ptr), intent(out) :: dmap
            integer(c_int), value :: lev
            type(c_ptr), value :: core
         end subroutine amrex_fi_get_distromap
      end interface
      call amrex_fi_get_distromap(dm%p,lvl,this%amrcore)
   end function get_distromap
   
   
   !> Build mfiter at a level
   subroutine mfiter_build(this,lvl,mfi)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_mfiter_build
      implicit none
      class(amrconfig), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(out) :: mfi
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      ba=this%get_boxarray (lvl)
      dm=this%get_distromap(lvl)
      call amrex_mfiter_build(mfi,ba,dm)
   end subroutine mfiter_build
   
   
   !> Destroy mfiter
   subroutine mfiter_destroy(this,mfi)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_destroy
      implicit none
      class(amrconfig), intent(inout) :: this
      type(amrex_mfiter), intent(inout) :: mfi
      call amrex_mfiter_destroy(mfi)
   end subroutine mfiter_destroy
   
   
   !> Return current finest level
   function clvl(this) result(cl)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrconfig), intent(inout) :: this
      integer :: cl
      interface
         integer(c_int) function amrex_fi_get_finest_level(core) bind(c)
            import :: c_int,c_ptr
            implicit none
            type(c_ptr), value :: core
         end function amrex_fi_get_finest_level
      end interface
      cl=amrex_fi_get_finest_level(this%amrcore)
   end function clvl
   
   
   !> Average down entire multifab array
   subroutine average_down(this,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_average_down
      implicit none
      class(amrconfig), intent(inout) :: this
      type(amrex_multifab), dimension(0:) :: mfab
      integer :: n
      do n=this%clvl()-1,0,-1
         call amrex_average_down(mfab(n+1),mfab(n),this%geom(n+1),this%geom(n),1,mfab(0)%nc,this%rref(n))
      end do
   end subroutine average_down
   
   
   !> Average entire multifab array down to level lvl
   subroutine average_downto(this,mfab,lvl)
      use amrex_amr_module, only: amrex_multifab,amrex_average_down
      implicit none
      class(amrconfig), intent(inout) :: this
      type(amrex_multifab), dimension(0:) :: mfab
      integer, intent(in) :: lvl
      call amrex_average_down(mfab(lvl+1),mfab(lvl),this%geom(lvl+1),this%geom(lvl),1,mfab(0)%nc,this%rref(lvl))
   end subroutine average_downto
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! AMRDATA specific routines !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

   !> Initialize an amrdata object
   subroutine initialize_amrdata(this,amr,ncomp,nover,reg,atface,name)
      use messager, only: die
      use amrex_amr_module, only: amrex_bc_int_dir,amrex_interp_cell_cons
      implicit none
      class(amrdata),   target, intent(inout) :: this
      class(amrconfig), target, intent(in)    :: amr
      integer, intent(in) :: ncomp
      integer, intent(in) :: nover
      integer, intent(in) :: reg
      logical, dimension(3), intent(in), optional :: atface
      character(len=*),      intent(in), optional :: name
      type(reg_list), pointer :: new_reg
      ! Point to amrconfig object
      this%amr=>amr
      ! Set name
      if (present(name)) this%name=trim(adjustl(name))
      ! Allocate data
      allocate(this%data(0:this%amr%nlvl))
      ! Check ncomp and store
      if (ncomp.lt.1) call die('[amrdata initialize] ncomp needs to be at least 1')
      this%ncomp=ncomp
      ! Check nover and store
      if (nover.lt.0) call die('[amrdata initialize] nover needs to be at least 0')
      this%nover=nover
      ! Handle amrdata registration
      this%reg=reg
      select case (this%reg)
      case (ignore)
         ! Nothing needs to be done
      case (resize,interp)
         ! Add our amrdata to amrconfig's data list
         allocate(new_reg); new_reg%data=>this
         new_reg%next=>this%amr%first_reg
         this%amr%first_reg=>new_reg
      case default
         call die('[amrdata initialize] Unknown registration status specified')
      end select
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
      ! Prepare storage for field info
      allocate(this%mindata(1:this%ncomp)); this%mindata=+huge(1.0_WP)
      allocate(this%maxdata(1:this%ncomp)); this%maxdata=-huge(1.0_WP)
      allocate(this%intdata(1:this%ncomp)); this%intdata= 0.0_WP
   end subroutine initialize_amrdata
   
   
   !> Finalize an amrdata object
   impure elemental subroutine finalize_amrdata(this)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrdata), intent(inout) :: this
      integer :: n
      do n=0,this%amr%nlvl
         call amrex_multifab_destroy(this%data(n))
      end do
      deallocate(this%data)
      if (.not.this%bc_ptr) then
         if (associated(this%lo_bc)) deallocate(this%lo_bc); nullify(this%lo_bc)
         if (associated(this%hi_bc)) deallocate(this%hi_bc); nullify(this%hi_bc)
      end if
      this%amr=>NULL()
      deallocate(this%mindata,this%maxdata,this%intdata)
   end subroutine finalize_amrdata
   
   
   !> Calculate various information on our amrdata object
   subroutine get_info_amrdata(this)
      implicit none
      class(amrdata), intent(inout) :: this
      integer :: lvl,n
      
      ! Reset info
      this%mindata=+huge(1.0_WP)
      this%maxdata=-huge(1.0_WP)
      this%intdata= 0.0_WP
      
      ! Loop over components
      do n=1,this%ncomp
         ! Loop over all levels
         do lvl=0,this%amr%clvl()
            ! Get min and max at that level
            this%mindata(n)=min(this%mindata(n),this%data(lvl)%min(comp=n))
            this%maxdata(n)=max(this%maxdata(n),this%data(lvl)%max(comp=n))
         end do
         ! Get int at level 0
         this%intdata(n)=this%data(0)%sum(comp=n)*(this%amr%dx(0)*this%amr%dy(0)*this%amr%dz(0))/this%amr%vol
      end do
      
   end subroutine get_info_amrdata
   
   
   !> Create data at level lvl
   subroutine create_lvl(this,lvl,time,ba,dm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray) , intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Delete data first
      call this%delete_lvl(lvl)
      ! Rebuild data
      call amrex_multifab_build(this%data(lvl),ba,dm,this%ncomp,this%nover,this%atface)
   end subroutine create_lvl
   

   !> Define data at level lvl
   subroutine define_lvl(this,lvl,time,ba,dm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_multifab_build,amrex_boxarray,amrex_distromap
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray) , intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Delete data first
      call this%delete_lvl(lvl)
      ! Rebuild data
      call amrex_multifab_build(this%data(lvl),ba,dm,this%ncomp,this%nover,this%atface)
   end subroutine define_lvl


   !> Refine data at level lvl
   subroutine refine_lvl(this,lvl,time,ba,dm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray) , intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Recreate data
      call this%create_lvl(lvl,time,ba,dm)
      ! Fill from coarse level
      call this%cfill_lvl(lvl,time,this%data(lvl))
   end subroutine refine_lvl
   
   
   !> Remake data at level lvl
   subroutine remake_lvl(this,lvl,time,ba,dm)
      use iso_c_binding,    only: c_ptr
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_multifab_build,amrex_multifab_destroy,amrex_multifab
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray) , intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_multifab)  :: newdata
      ! Create newdata and fill it from current data
      call amrex_multifab_build(newdata,ba,dm,this%ncomp,this%nover,this%atface)
      call this%fill_lvl(lvl,time,newdata)
      ! Recreate level
      call this%create_lvl(lvl,time,ba,dm)
      ! Copy newdata to data
      call this%data(lvl)%copy(newdata,1,1,this%ncomp,this%nover)
      ! Destroy newdata
      call amrex_multifab_destroy(newdata)
   end subroutine remake_lvl
   
   
   !> Delete data at level lvl
   subroutine delete_lvl(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrdata), intent(inout) :: this
      integer, intent(in) :: lvl
      call amrex_multifab_destroy(this%data(lvl))
   end subroutine delete_lvl
   
   
   !> Fill provided mfab at level (lvl) from our data at level (lvl-1)
   subroutine cfill_lvl(this,lvl,time,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_fillcoarsepatch,amrex_interp_cell_cons
      implicit none
      class(amrdata), intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
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
   subroutine fill_lvl(this,lvl,time,mfab)
      use amrex_amr_module, only: amrex_multifab,amrex_fillpatch,amrex_interp_cell_cons
      implicit none
      class(amrdata),   intent(inout) :: this
      integer,  intent(in) :: lvl
      real(WP), intent(in) :: time
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
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Generic AMReX routines !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      call amrex_finalize()
   end subroutine finalize_amrex
   
   
end module amrconfig_class
   