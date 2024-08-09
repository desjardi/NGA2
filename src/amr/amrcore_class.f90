!> AMR core object is defined here
module amrcore_class
   use iso_c_binding
   use precision, only: WP
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrcore,finalize_amrex
   
   !> Amrcore object definition based on AMReX's amrcore
   type :: amrcore
      
      ! Pointer to amrcore object
      type(c_ptr) :: amrcore=c_null_ptr
      
      ! Coordinate system
      integer :: coordsys=0
      ! Domain periodicity (0=non-periodic,1=periodic)
      integer :: xper,yper,zper
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
      
   contains
      procedure :: initialize
      procedure :: finalize
   end type amrcore
   
   
contains
   
   
   !> Initialization of an amrcore object
   subroutine initialize(this)
      implicit none
      class(amrcore), intent(inout) :: this
      
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
         call amrex_parmparse_build(pp,'amr')
         call pp%addarr('n_cell'         ,[this%nx,this%ny,this%nz])
         call pp%add   ('max_level'      ,this%nlvl)
         call pp%add   ('blocking_factor',this%nbloc)
         call pp%add   ('max_grid_size'  ,this%nmax)
         call pp%add   ('ref_ratio'      ,this%rref)
         call amrex_parmparse_destroy(pp)
         call amrex_parmparse_build(pp,'geometry')
         call pp%add   ('coord_sys'  ,this%coordsys)
         call pp%addarr('is_periodic',[this%xper,this%yper,this%zper])
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
      end block create_amrcore_obj
      
   end subroutine initialize
   
   
   !> Finalization of amrcore object
   subroutine finalize(this)
      implicit none
      class(amrcore), intent(inout) :: this
      ! Destroy our amrcore object
      destroy_amrcore_obj: block
         interface
            subroutine amrex_fi_delete_amrcore(core) bind(c)
               import
               implicit none
               type(c_ptr), value :: core
            end subroutine amrex_fi_delete_amrcore
         end interface
         call amrex_fi_delete_amrcore(this%amrcore)
      end block destroy_amrcore_obj
      this%amrcore=c_null_ptr
   end subroutine finalize
   
   
   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      call amrex_finalize()
   end subroutine finalize_amrex
   
   
end module amrcore_class
   