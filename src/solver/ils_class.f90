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
   integer, parameter, public :: hypre_smg=2
   
   
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
      real(WP), dimension(:,:,:),   allocatable :: sol                !< Solution
      ! An ILS has convergence criteria
      integer  :: maxit                                               !< Maximum number of iterators allowed
      real(WP) :: rcvg                                                !< Desired relative convergence criterion
      real(WP) :: acvg                                                !< Desired absolue convergence criterion
      ! Current info
      integer  :: it                                                  !< Number of iterations
      real(WP) :: rerr                                                !< Current relative error
      real(WP) :: aerr                                                !< Current absolute error
      
      ! Private stuff for hypre
      integer(kind=8), private :: hypre_box              !< Grid
      integer(kind=8), private :: hypre_stc              !< Stencil
      integer(kind=8), private :: hypre_mat              !< Matrix
      integer(kind=8), private :: hypre_rhs              !< Right-hand side
      integer(kind=8), private :: hypre_sol              !< Solution
      integer(kind=8), private :: hypre_solver           !< Solver
      
   contains
      procedure :: print=>ils_print                                   !< Print ILS object info
      procedure :: prep_solver                                        !< Solver preparation
      procedure :: solve                                              !< Equation solver
      procedure :: solve_rbgs                                         !< Solve equation with RBGS
      procedure :: solve_hypre_smg                                    !< Solve equation with hypre_smg
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
      if (allocated(this%sol)) deallocate(this%sol)
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
      
      ! Allocate operator, rhs, and sol arrays
      allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
      allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
      allocate(self%sol(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%sol=0.0_WP
      
   end function ils_from_args
   
   
   !> Prepare solver
   subroutine prep_solver(this)
      use monitor, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer :: ierr,st
      integer, dimension(3) :: periodicity
      integer, dimension(6) :: ghosts
      
      ! Select appropriate solver
      select case (this%method)
      case (rbgs)                    ! Initialize the Red-Black Gauss-Seidel solver here
         ! Nothing to do here
         
      case (hypre_smg)               ! Initialize the HYPRE-SMG solver here
         ! Create Hypre grid
         call HYPRE_StructGridCreate(this%cfg%comm,3,this%hypre_box,ierr)
         call HYPRE_StructGridSetExtents(this%hypre_box,[this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_],[this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_],ierr)
         periodicity=0
         if (this%cfg%xper) periodicity(1)=this%cfg%nx
         if (this%cfg%yper) periodicity(2)=this%cfg%ny
         if (this%cfg%zper) periodicity(3)=this%cfg%nz
         call HYPRE_StructGridSetPeriodic(this%hypre_box,periodicity,ierr)
         call HYPRE_StructGridAssemble(this%hypre_box,ierr)
         ! Build Hypre stencil
         call HYPRE_StructStencilCreate(3,this%nst,this%hypre_stc,ierr)
         do st=1,this%nst
            call HYPRE_StructStencilSetElement(this%hypre_stc,st-1,[this%stc(st,1),this%stc(st,2),this%stc(st,3)],ierr)
         end do
         ! Build Hypre matrix
         call HYPRE_StructMatrixCreate(this%cfg%comm,this%hypre_box,this%hypre_stc,this%hypre_mat,ierr)
         ghosts=maxval(abs(this%stc))
         call HYPRE_StructMatrixSetNumGhost(this%hypre_mat,ghosts,ierr)
         call HYPRE_StructMatrixInitialize(this%hypre_mat,ierr)
         ! Prepare Hypre RHS
         call HYPRE_StructVectorCreate(this%cfg%comm,this%hypre_box,this%hypre_rhs,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_rhs,ierr)
         ! Create Hypre solution vector
         call HYPRE_StructVectorCreate(this%cfg%comm,this%hypre_box,this%hypre_sol,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_sol,ierr)
         
      case default
         call die('[ils prep solver] Unknown solution method')
      end select
      
   end subroutine prep_solver
   
   
   !> Solve the linear system iteratively
   subroutine solve(this)
      use monitor, only: die
      implicit none
      class(ils), intent(inout) :: this
      ! Select appropriate solver
      select case (this%method)
      case (rbgs)
         call this%solve_rbgs()
      case (hypre_smg)
         call this%solve_hypre_smg()
      case default
         call die('[ils solve] Unknown solution method')
      end select
   end subroutine solve
   
   
   !> Solve the linear system using RBGS
   subroutine solve_rbgs(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,col,st
      integer :: colmin,colmax,coldir
      real(WP) :: val
      
      ! Reset iteration and error
      this%it=0
      this%rerr=huge(1.0_WP)
      this%aerr=huge(1.0_WP)
      
      ! Loop until done
      do while (this%it.lt.this%maxit.and.this%rerr.ge.this%rcvg.and.this%aerr.ge.this%acvg)
         
         ! Increment counter
         this%it=this%it+1
         
         ! Alternate sweep direction
         if (mod(this%it,2).eq.0) then
            colmin=8; colmax=1; coldir=-1 !< Negative sweep
         else
            colmin=1; colmax=8; coldir=+1 !< Positive sweep
         end if
         
         ! Loop over colors
         do col=colmin,colmax,coldir
            
            ! Loop over domain
            do k=this%cfg%kmin_+mod((col-1)/4,2),this%cfg%kmax_,2
               do j=this%cfg%jmin_+mod((col-1)/2,2),this%cfg%jmax_,2
                  do i=this%cfg%imin_+mod((col-1)/1,2),this%cfg%imax_,2
                     ! Gauss-Seidel step
                     if (abs(this%opr(1,i,j,k)).gt.0.0_WP) then
                        val=this%rhs(i,j,k)
                        do st=2,this%nst
                           val=val-this%opr(st,i,j,k)*this%sol(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
                        end do
                        this%sol(i,j,k)=val/this%opr(1,i,j,k)
                     end if
                  end do
               end do
            end do
            
            ! Communicate solution
            call this%cfg%sync(this%sol)
            
         end do
      end do
      
   end subroutine solve_rbgs
   
   
   !> Solve the linear system using hypre_smg
   subroutine solve_hypre_smg(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,st,ierr
      
      ! Build the problem
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Transfer operator
               do st=1,this%nst
                  call HYPRE_StructMatrixSetValues(this%hypre_mat,[i,j,k],1,st-1,this%opr(st,i,j,k),ierr)
               end do
               ! Tranfer RHS
               call HYPRE_StructVectorSetValues(this%hypre_rhs,[i,j,k],this%rhs(i,j,k),ierr)
               ! Tranfer initial guess
               call HYPRE_StructVectorSetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      call HYPRE_StructMatrixAssemble(this%hypre_mat,ierr)
      call HYPRE_StructVectorAssemble(this%hypre_rhs,ierr)
      call HYPRE_StructVectorAssemble(this%hypre_sol,ierr)
      
      ! Prepare the solver
      call HYPRE_StructSMGCreate    (this%cfg%comm,this%hypre_solver,ierr)
      call HYPRE_StructSMGSetMaxIter(this%hypre_solver,this%maxit,ierr)
      call HYPRE_StructSMGSetTol    (this%hypre_solver,this%rcvg ,ierr)
      call HYPRE_StructSMGSetLogging(this%hypre_solver,1         ,ierr)
      call HYPRE_StructSMGSetup     (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      print*,'here3'
      ! Solve
      call HYPRE_StructSMGSolve           (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      call HYPRE_StructSMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
      call HYPRE_StructSMGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      call HYPRE_StructSMGDestroy         (this%hypre_solver,ierr)
      print*,'here4'
      ! Get back solution
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Transfer solution
               call HYPRE_StructVectorGetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      print*,'here5'
   end subroutine solve_hypre_smg
   
   
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
