!> Diagonal dominant alternating direction implicit (DDADI) solver object is defined in this class.
!> Based on Pulliam and Chaussee's approximate-factorization algorithm.

module ddadi_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_short
   use diag_class,   only: diag
   use linsol_class, only: linsol
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: ddadi
   
   
   !> ddadi solver object definition
   type, extends(linsol) :: ddadi
      
      !> Diagonal object
      class(diag), allocatable :: dsol
      
   contains
      
      procedure :: print_short=>ddadi_print_short !< One-line printing of solver status
      procedure :: print=>ddadi_print             !< Long-form printing of solver status
      procedure :: log=>ddadi_log                 !< Long-form logging of solver status
      procedure :: init=>ddadi_init               !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>ddadi_setup             !< Solver setup (every time the operator changes)
      procedure :: solve=>ddadi_solve             !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>ddadi_destroy         !< Solver destruction (every time the operator changes)
      
   end type ddadi
   
   
   !> Declare ddadi constructor
   interface ddadi
      procedure ddadi_from_args
   end interface ddadi
   
   
contains
   
   
   !> Constructor for a ddadi object
   function ddadi_from_args(cfg,name,nst) result(self)
      use messager, only: die
      implicit none
      type(ddadi) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      integer, intent(in) :: nst
      
      ! Link the config and store the name
      self%cfg=>cfg
      self%name=trim(adjustl(name))
      
      ! Set solution method - not used
      self%method=0
      
      ! Set up stencil size and map
      self%nst=nst
      allocate(self%stc(1:self%nst,1:3))
      self%stc=0
      
      ! Allocate operator, rhs, and sol arrays
      allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
      allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
      allocate(self%sol(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%sol=0.0_WP
      
      ! Zero out some info
      self%it=1
      self%aerr=0.0_WP
      self%rerr=0.0_WP
      
      ! Setup is not done
      self%setup_done=.false.
      
   end function ddadi_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine ddadi_init(this)
      use messager, only: die
      implicit none
      class(ddadi), intent(inout) :: this
      integer :: st,stx1,stx2,sty1,sty2,stz1,stz2
      
      ! From the provided stencil, generate an inverse map
      stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1))
      sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2))
      stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3))
      allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
      do st=1,this%nst
         this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
      end do
      
      ! Initialize diagonal solver now that we know the stencil size
      allocate(this%dsol,source=diag(cfg=this%cfg,name=this%name,n=2*maxval(abs(this%stc))+1))
      
   end subroutine ddadi_init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine ddadi_setup(this)
      use messager, only: die
      implicit none
      class(ddadi), intent(inout) :: this
      integer :: i,j,k,st
      
      ! Update DDADI operators
      this%dsol%Ax=0.0_WP
      this%dsol%Ay=0.0_WP
      this%dsol%Az=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               do st=1,this%nst
                  this%dsol%Ax(j,k,i,this%stc(st,1))=this%dsol%Ax(j,k,i,this%stc(st,1))+this%opr(st,i,j,k)
                  this%dsol%Ay(i,k,j,this%stc(st,2))=this%dsol%Ay(i,k,j,this%stc(st,2))+this%opr(st,i,j,k)
                  this%dsol%Az(i,j,k,this%stc(st,3))=this%dsol%Az(i,j,k,this%stc(st,3))+this%opr(st,i,j,k)
               end do
            end do
         end do
      end do
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine ddadi_setup
   
   
   !> Solve the linear system iteratively
   subroutine ddadi_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(ddadi), intent(inout) :: this
      integer :: i,j,k
      
      ! Check that setup was done
      if (.not.this%setup_done) call die('[ddadi solve] Solver has not been setup.')
      
      ! Inverse in X-direction
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%dsol%Rx(j,k,i)=this%rhs(i,j,k)
            end do
         end do
      end do
      call this%dsol%linsol_x()
      
      ! Inverse in Y-direction
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%dsol%Ry(i,k,j)=this%dsol%Rx(j,k,i)*sum(this%opr(:,i,j,k))
            end do
         end do
      end do
      call this%dsol%linsol_y()
      
      ! Inverse in Z-direction
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%dsol%Rz(i,j,k)=this%dsol%Ry(i,k,j)*sum(this%opr(:,i,j,k))
            end do
         end do
      end do
      call this%dsol%linsol_z()
      
      ! Update solution vector
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%sol(i,j,k)=this%dsol%Rz(i,j,k)
            end do
         end do
      end do
      
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
      
   end subroutine ddadi_solve
   
   
   !> Log ddadi info
   subroutine ddadi_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(ddadi), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("ddadi solver [",a,"] for config [",a,"]")') trim(this%name), trim(this%cfg%name)
         call log(message)
      end if
   end subroutine ddadi_log
   
   
   !> Print ddadi info to the screen
   subroutine ddadi_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(ddadi), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("ddadi solver [",a,"] for config [",a,"]")') trim(this%name), trim(this%cfg%name)
   end subroutine ddadi_print
   
   
   !> Short print of ddadi info to the screen
   subroutine ddadi_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(ddadi), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("ddadi solver [",a16,"] for config [",a16,"]")') trim(this%name),trim(this%cfg%name)
   end subroutine ddadi_print_short
   

   !> Destroy ddadi solver
   subroutine ddadi_destroy(this)
      implicit none
      class(ddadi), intent(inout) :: this
      deallocate(this%stc,this%opr,this%rhs,this%sol,this%stmap)
   end subroutine ddadi_destroy
   
end module ddadi_class
