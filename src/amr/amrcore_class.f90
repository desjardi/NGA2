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
      integer :: coordsys
      ! Domain periodicity
      logical :: xper,yper,zper
      ! Domain extent
      real(WP) :: xlo,ylo,zlo
      real(WP) :: xhi,yhi,zhi
      
      ! Mesh dimensions
      integer :: nx,ny,nz
      ! Max grid size
      integer :: nmax
      ! Blocking factor
      integer :: nbloc
      ! Refinement ratio
      integer, dimension(4) :: rref
      
   contains
      procedure :: initialize
      procedure :: finalize
   end type amrcore
   
   interface
      subroutine amrex_fi_new_amrcore(my_amrcore) bind(c)
         import
         implicit none
         type(c_ptr) :: my_amrcore
      end subroutine amrex_fi_new_amrcore
      subroutine amrex_fi_delete_amrcore(my_amrcore) bind(c)
         import
         implicit none
         type(c_ptr) :: my_amrcore
      end subroutine amrex_fi_delete_amrcore
   end interface
   
contains
   
   
   !> Initialization of an amrcore object
   subroutine initialize(this)
      implicit none
      class(amrcore), intent(inout) :: this
      
      ! First of all, confirm that nga2 and amrex reals are compatible
      check_real: block
         use messager,          only: die
         use precision,         only: WP
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
         call pp%addarr('n_cell',[64,64,64])
         call pp%add   ('max_level',2)
         call pp%add   ('blocking_factor',8)
         call pp%add   ('max_grid_size',32)
         call pp%add   ('ref_ratio',2)
         call amrex_parmparse_destroy(pp)
         ! Prefix 'geometry'
         call amrex_parmparse_build(pp,'geometry')
         call pp%addarr('is_periodic',[1,1,1])
         call pp%add   ('coord_sys',0) ! 0 is cartesian
         call pp%addarr('prob_lo',[0.0_WP,0.0_WP,0.0_WP])
         call pp%addarr('prob_hi',[1.0_WP,1.0_WP,1.0_WP])
         call amrex_parmparse_destroy(pp)
      end block set_params

      ! Create an amrcore object
      call amrex_fi_new_amrcore(this%amrcore)
      
   end subroutine initialize
   
   
   !> Finalization of amrcore object
   subroutine finalize(this)
      implicit none
      class(amrcore), intent(inout) :: this
      ! Destroy our amrcore object
      call amrex_fi_delete_amrcore(this%amrcore)
      this%amrcore=c_null_ptr
   end subroutine finalize


   !> Finalization of amrex
   subroutine finalize_amrex()
      use amrex_amr_module, only: amrex_finalize
      implicit none
      ! Destroy our amrcore object
      call amrex_finalize()
   end subroutine finalize_amrex
   
   
end module amrcore_class
   