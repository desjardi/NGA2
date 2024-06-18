!> Break-up model class: uses cclabel to identify thin structures and break them up into droplets
module breakup_class
   use precision,     only: WP
   use config_class,  only: config
   use cclabel_class, only: cclabel
   use lpt_class,     only: lpt
   use vfs_class,     only: vfs
   use tpns_class,    only: tpns
   use irl_fortran_interface
   implicit none
   private
   
   ! Expose type/methods
   public :: breakup
   
   !> Break-up model object
   type :: breakup
      
      !> A cclabel object is needed
      type(cclabel) :: cc
      
      !> Pointers to lpt, vfs, and tpns
      class(vfs),  pointer :: vf
      class(tpns), pointer :: fs
      class(lpt),  pointer :: lp
      
      !> Break-up model parameters
      real(WP) :: filmthickness_over_dx  =5.0e-1_WP
      real(WP) :: min_filmthickness      =1.0e-4_WP
      real(WP) :: diam_over_filmthickness=1.0e+1_WP
      real(WP) :: max_eccentricity       =5.0e-1_WP
      real(WP) :: d_threshold            =1.0e-1_WP
      
   contains
      procedure :: initialize
      !procedure :: film_breakup
      !procedure :: drop_transfer
   end type breakup
   
contains
   
   !> Initialize the breakup model
   subroutine initialize(this,vf,fs,lp)
      implicit none
      class(breakup), intent(inout) :: this
      class(vfs),  target, intent(in) :: vf
      class(tpns), target, intent(in) :: fs
      class(lpt),  target, intent(in) :: lp
      ! Store pointers to our solvers
      this%vf=>vf
      this%fs=>fs
      this%lp=>lp
      ! Create a connected-component labeling object
      !call this%cc%initialize(cfg=this%vf%cfg,name='CCL')
   end subroutine initialize
   
   
end module