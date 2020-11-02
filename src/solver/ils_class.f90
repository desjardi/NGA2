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
   integer, parameter, public :: amg =2
   
   ! Hypre-related storage parameter
   integer(kind=8), parameter :: hypre_ParCSR=5555
   
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
      integer(kind=8), private :: hypre_mat,parse_mat    !< Matrix
      integer(kind=8), private :: hypre_rhs,parse_rhs    !< Right-hand side
      integer(kind=8), private :: hypre_sol,parse_sol    !< Solution
      integer(kind=8), private :: hypre_solver           !< Solver
      
   contains
      procedure :: print=>ils_print                                   !< Print ILS object info
      procedure :: log  =>ils_log                                     !< Log ILS object info
      procedure :: init_solver                                        !< Solver initialization (at start-up)
      procedure :: update_solver                                      !< Solver update (every time the operator changes)
      procedure :: solve                                              !< Equation solver
      procedure :: solve_rbgs                                         !< Solve equation with rbgs
      procedure :: solve_hypre_amg                                    !< Solve equation with amg
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
      use messager, only: die
      implicit none
      type(ils) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      integer, optional :: nst
      
      ! Link the config and store the name
      self%cfg=>cfg
      self%name=trim(adjustl(name))
      
      ! Zero out some info
      self%it=0; self%aerr=0.0_WP; self%rerr=0.0_WP
      
      ! Set up stencil size and map
      self%nst=7; if (present(nst)) self%nst=nst
      allocate(self%stc(1:self%nst,1:3)); self%stc=0
      
      ! Allocate operator, rhs, and sol arrays
      allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
      allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
      allocate(self%sol(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%sol=0.0_WP
      
   end function ils_from_args
   
   
   !> Initialize solver - done at start-up (probably only once?)
   subroutine init_solver(this,method)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer, intent(in) :: method
      integer :: ierr
      
      ! Set solution method
      this%method=method
      
      ! Select appropriate solver
      select case (this%method)
      case (rbgs) ! Initialize the Red-Black Gauss-Seidel solver here
         
      case (amg)  ! Initialize the HYPRE-AMG solver here
         ! Prepare an unstructured mapping for our cfg
         call this%cfg%umap_prep()
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%cfg%ind_min,this%cfg%ind_max,this%cfg%ind_min,this%cfg%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         call HYPRE_IJMatrixInitialize   (this%hypre_mat,ierr)
         ! Create a HYPRE rhs vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%cfg%ind_min,this%cfg%ind_max,this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
         ! Create a HYPRE solution vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%cfg%ind_min,this%cfg%ind_max,this%hypre_sol,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
         ! Create a HYPRE AMG solver
         call HYPRE_BoomerAMGCreate        (this%hypre_solver,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_solver,3,ierr)            ! print solve info + parameters (3 from all, 0 for none)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_solver,6,ierr)            ! interpolation default is 6 (NGA=3)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_solver,6,ierr)           ! Falgout=6 (old default, NGA); PMIS=8 and HMIS=10 (recommended)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_solver,0.15_WP,ierr)    ! 0.25 is default
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_solver,8,ierr)            ! hybrid symmetric Gauss-Seidel/SOR
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_solver,this%maxit,ierr)   ! maximum nbr of iter
         call HYPRE_BoomerAMGSetTol        (this%hypre_solver,this%rcvg ,ierr)   ! convergence tolerance
         
      case default
         call die('[ils prep solver] Unknown solution method')
      end select
      
      ! Finally, update the solver
      call this%update_solver()
      
   end subroutine init_solver
   
   
   !> Update solver - done everytime the operator changes
   subroutine update_solver(this)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count1,count2,st,ierr
      integer,  dimension(:), allocatable :: row,ncol,col
      real(WP), dimension(:), allocatable :: val
      
      ! Select appropriate solver
      select case (this%method)
      case (rbgs) ! Prepare the Red-Black Gauss-Seidel solver here
         
      case (amg)  ! Prepare the HYPRE-AMG solver here
         
         ! Allocate storage for transfer
         allocate( row(         this%cfg%ncell_))
         allocate(ncol(         this%cfg%ncell_))
         allocate( col(this%nst*this%cfg%ncell_))
         allocate( val(this%nst*this%cfg%ncell_))
         
         ! Tranfer operator to HYPRE
         count1=0; count2=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%cfg%ind(i,j,k).gt.0) then
                     count1=count1+1
                     row (count1)=this%cfg%ind(i,j,k)
                     ncol(count1)=0
                     do st=1,this%nst
                        if (this%cfg%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3)).gt.0) then
                           count2=count2+1
                           ncol(count1)=ncol(count1)+1
                           col (count2)=this%cfg%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
                           val (count2)=this%opr(st,i,j,k)
                        end if
                     end do
                  end if
               end do
            end do
         end do
         call HYPRE_IJMatrixSetValues(this%hypre_mat,this%cfg%ncell_,ncol,row,col,val,ierr)
         call HYPRE_IJMatrixAssemble (this%hypre_mat,ierr)
         call HYPRE_IJMatrixGetObject(this%hypre_mat,this%parse_mat,ierr)
         
         ! Setup AMG solver
         call HYPRE_BoomerAMGSetup(this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
         ! Deallocate
         deallocate(row,ncol,col,val)
         
      case default
         call die('[ils prep solver] Unknown solution method')
      end select
      
   end subroutine update_solver
   
   
   !> Solve the linear system iteratively
   subroutine solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(ils), intent(inout) :: this
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      ! Select appropriate solver
      select case (this%method)
      case (rbgs)
         call this%solve_rbgs()
      case (amg)
         call this%solve_hypre_amg()
      case default
         call die('[ils solve] Unknown solution method')
      end select
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print
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
   
   
   !> Solve the linear system using hypre_amg
   subroutine solve_hypre_amg(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count,ierr
      integer,  dimension(:), allocatable :: ind
      real(WP), dimension(:), allocatable :: rhs,sol
      
      ! Transfer the rhs and initial guess to hypre
      allocate(rhs(this%cfg%ncell_))
      allocate(sol(this%cfg%ncell_))
      allocate(ind(this%cfg%ncell_))
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%ind(i,j,k).gt.0) then
                  count=count+1
                  ind(count)=this%cfg%ind(i,j,k)
                  rhs(count)=this%rhs(i,j,k)
                  sol(count)=this%sol(i,j,k)
               end if
            end do
         end do
      end do
      call HYPRE_IJVectorSetValues(this%hypre_rhs,this%cfg%ncell_,ind,rhs,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_rhs,ierr)
      call HYPRE_IJVectorSetValues(this%hypre_sol,this%cfg%ncell_,ind,sol,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_sol,ierr)
      
      ! Call the solver
      call HYPRE_BoomerAMGSolve           (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      call HYPRE_BoomerAMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
      call HYPRE_BoomerAMGGetFinalReltvRes(this%hypre_solver,this%rerr,ierr)
      
      ! Get the solution back
      call HYPRE_IJVectorGetValues(this%hypre_sol,this%cfg%ncell_,ind,sol,ierr)
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%cfg%ind(i,j,k).gt.0) then
                  count=count+1
                  this%sol(i,j,k)=sol(count)
               end if
            end do
         end do
      end do
      
      ! Deallocate
      deallocate(rhs,sol,ind)
      
   end subroutine solve_hypre_amg
   
   
   !> Log ILS info
   subroutine ils_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(ils), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Iterative Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
         select case (this%method)
         case (rbgs)
            write(message,'(" > method = ",a)') 'RBGS'; call log(message)
         case (amg)
            write(message,'(" > method = ",a)') 'HYPRE AMG'; call log(message)
         case default
            write(message,'(" > method = ",a)') 'unknown'; call log(message)
         end select
         write(message,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit; call log(message)
         write(message,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg; call log(message)
         write(message,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg; call log(message)
      end if
   end subroutine ils_log
   
   
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
         case (amg)
            write(output_unit,'(" > method = ",a)') 'HYPRE AMG'
         case default
            write(output_unit,'(" > method = ",a)') 'unknown'
         end select
         write(output_unit,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit
         write(output_unit,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg
         write(output_unit,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg
      end if
   end subroutine ils_print
   
   
end module ils_class
