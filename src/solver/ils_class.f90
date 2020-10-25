!> Iterative linear solver concept is defined here:
!> given a config, it provides parallel solvers
module ils_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ils
   
   
   ! List of known available methods
   integer, parameter, public :: rbgs=1
   
   
   !> ILS object definition
   type :: ils
      ! An iterative linear solver works for a config
      type(config), pointer :: cfg                                    !< Config for the ILS
      ! An ILS has a filename
      character(len=str_medium) :: name                               !< Name of solver
      ! An ILS has a solution method
      integer  :: method                                              !< Solution method
      ! An ILS has a stencil size
      integer  :: nst                                                 !< Stencil size in 3D
      integer, dimension(:,:), allocatable :: stc                     !< Stencil map in 3D
      ! An ILS stores the linear operator and rhs
      real(WP), dimension(:,:,:,:), allocatable :: opr                !< Linear operator
      real(WP), dimension(:,:,:),   allocatable :: rhs                !< RHS
      ! An ILS has convergence criteria
      integer  :: maxit                                               !< Maximum number of iterators allowed
      real(WP) :: rcvg                                                !< Desired relative convergence criterion
      real(WP) :: acvg                                                !< Desired absolue convergence criterion
      ! Current info
      integer  :: it                                                  !< Number of iterations
      real(WP) :: rerr                                                !< Current relative error
      real(WP) :: aerr                                                !< Current absolute error
   contains
      procedure :: print=>ils_print                                   !< Print ILS object info
      procedure :: scale_opr                                          !< Operator scaling
      procedure :: scale_rhs                                          !< RHS scaling
      final     :: destructor                                         !< Destructor for ILS
   end type ils
   
   
   !> Declare ils constructor
   interface ils
      procedure ils_from_args
   end interface ils
   
   
contains
   
   
   !> Destructor for ILS object
   subroutine destructor(this)
      implicit none
      type(ils) :: this
      if (allocated(this%stc)) deallocate(this%stc)
      if (allocated(this%opr)) deallocate(this%opr)
      if (allocated(this%rhs)) deallocate(this%rhs)
   end subroutine destructor
   
   
   !> Constructor for an ILS object
   function ils_from_args(cfg,name,nst) result(self)
      use monitor, only: die
      implicit none
      type(ils) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      integer, optional :: nst
      
      ! Link the config and store the name
      self%cfg=>cfg
      self%name=trim(adjustl(name))
      
      ! Set up stencil size and map
      self%nst=7; if (present(nst)) self%nst=nst
      allocate(self%stc(1:self%nst,1:3)); self%stc=0
      
      ! Allocate operator and RHS arrays
      allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
      allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
      
   end function ils_from_args
   
   
   !> Scaling of operator - this is simple negative volume scaling here, but could easily be generalized
   subroutine scale_opr(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%opr(:,i,j,k)=-this%opr(:,i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do
   end subroutine scale_opr
   
   
   !> Scaling of rhs - this is simple negative volume scaling here, but could easily be generalized
   subroutine scale_rhs(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               this%rhs(i,j,k)=-this%rhs(i,j,k)*this%cfg%vol(i,j,k)
            end do
         end do
      end do
   end subroutine scale_rhs
   
   
   !> Print ILS info to the screen
   subroutine ils_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(ils), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("Iterative Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         select case (this%method)
         case (rbgs)
            write(output_unit,'(" > method = ",a)') 'RBGS'
         case default
            write(output_unit,'(" > method = ",a)') 'unknown'
         end select
         write(output_unit,'(" >  maxit = ",i0)') this%maxit
         write(output_unit,'(" >   acvg = ",es12.5)') this%acvg
         write(output_unit,'(" >   rcvg = ",es12.5)') this%rcvg
      end if
   end subroutine ils_print
   
   
end module ils_class
