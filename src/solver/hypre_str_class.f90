!> Structured hypre linear solver concept is
!> defined by extension of the linsol class
module hypre_str_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   use linsol_class, only: linsol
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: hypre_str
   
   
   ! List of known available methods
   integer, parameter, public ::       pfmg=1
   integer, parameter, public ::   pcg_pfmg=2
   integer, parameter, public :: gmres_pfmg=3
   integer, parameter, public ::        smg=4
   integer, parameter, public ::    pcg_smg=5
   integer, parameter, public ::  gmres_smg=6
   integer, parameter, public ::        pcg=7
   integer, parameter, public ::      gmres=8
   integer, parameter, public ::  pcg_pfmg2=9
   
   
   ! List of key solver parameters
   integer , parameter :: sprintlvl=0                  !< Solver  printing: 0=none (default); 3=init and cvg history
   integer , parameter :: pprintlvl=0                  !< Precond printing: 0=none (default); 3=init and cvg history
   
   
   !> hypre_str object definition
   type, extends(linsol) :: hypre_str
      
      ! For multigrid solvers, maxlevel parameter is needed
      integer :: maxlevel                                !< Maximum number of multigrid levels
      
      ! Private stuff for hypre
      integer(kind=8), private :: hypre_box              !< Grid
      integer(kind=8), private :: hypre_stc              !< Stencil
      integer(kind=8), private :: hypre_mat              !< Matrix
      integer(kind=8), private :: hypre_rhs              !< Right-hand side
      integer(kind=8), private :: hypre_sol              !< Solution
      integer(kind=8), private :: hypre_solver           !< Solver
      integer(kind=8), private :: hypre_precond          !< Preconditioner
      
   contains
      
      procedure :: print_short=>hypre_str_print_short    !< One-line printing of solver status
      procedure :: print=>hypre_str_print                !< Long-form printing of solver status
      procedure :: log=>hypre_str_log                    !< Long-form logging of solver status
      procedure :: init=>hypre_str_init                  !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>hypre_str_setup                !< Solver setup (every time the operator changes)
      procedure :: solve=>hypre_str_solve                !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>hypre_str_destroy            !< Solver destruction (every time the operator changes)
      
   end type hypre_str
   
   
   !> Declare hypre_str constructor
   interface hypre_str
      procedure hypre_str_from_args
   end interface hypre_str
   
   
contains
   
   
   !> Constructor for an hypre_str object
   function hypre_str_from_args(cfg,name,method,nst) result(self)
      use messager, only: die
      implicit none
      type(hypre_str) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      integer, intent(in) :: method
      integer, intent(in) :: nst
      
      ! Link the config and store the name
      self%cfg=>cfg
      self%name=trim(adjustl(name))
      
      ! Set solution method
      self%method=method
      
      ! Set up stencil size and map
      self%nst=nst
      allocate(self%stc(1:self%nst,1:3))
      self%stc=0
      
      ! Allocate operator, rhs, and sol arrays
      allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
      allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
      allocate(self%sol(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%sol=0.0_WP
      
      ! Zero out some info
      self%it=0
      self%aerr=0.0_WP
      self%rerr=0.0_WP
      self%maxlevel=0
      
      ! Setup is not done
      self%setup_done=.false.
      
   end function hypre_str_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine hypre_str_init(this)
      use messager, only: die
      implicit none
      class(hypre_str), intent(inout) :: this
      integer :: ierr,st,stx1,stx2,sty1,sty2,stz1,stz2
      integer, dimension(3) :: periodicity,offset
      
      ! From the provided stencil, generate an inverse map
      stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1))
      sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2))
      stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3))
      allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
      do st=1,this%nst
         this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
      end do
      
      ! Initialize HYPRE
      call HYPRE_Init(ierr)
      
      ! These use HYPRE's structured environment, which requires that we create a HYPRE grid and stencil
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
         offset=this%stc(st,:)
         call HYPRE_StructStencilSetElement(this%hypre_stc,st-1,offset,ierr)
      end do
      
   end subroutine hypre_str_init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine hypre_str_setup(this)
      use messager, only: die
      implicit none
      class(hypre_str), intent(inout) :: this
      integer :: i,j,k,st,ierr
      integer,  dimension(:), allocatable :: row
      real(WP), dimension(:), allocatable :: val
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      ! Create a structured matrix
      call HYPRE_StructMatrixCreate(this%cfg%comm,this%hypre_box,this%hypre_stc,this%hypre_mat,ierr)
      call HYPRE_StructMatrixInitialize(this%hypre_mat,ierr)
      
      ! Prepare local storage
      allocate(row(1:this%nst),val(1:this%nst))
      do st=1,this%nst
         row(st)=st-1
      end do
      
      ! Transfer the operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               val(1)=1.0_WP; val(2:)=0.0_WP
               if (abs(this%opr(1,i,j,k)).gt.0.0_WP) val=this%opr(:,i,j,k)
               call HYPRE_StructMatrixSetValues(this%hypre_mat,[i,j,k],this%nst,row,val,ierr)
            end do
         end do
      end do
      call HYPRE_StructMatrixAssemble(this%hypre_mat,ierr)
      
      ! Deallocate storage
      deallocate(row,val)
      
      ! For debugging - output the operator to a file
      !call HYPRE_StructMatrixPrint('struct_mat'//char(0),this%hypre_mat,0,ierr)
      !call die('Matrix was printed out')
      
      ! Prepare structured rhs vector
      call HYPRE_StructVectorCreate(this%cfg%comm,this%hypre_box,this%hypre_rhs,ierr)
      call HYPRE_StructVectorInitialize(this%hypre_rhs,ierr)
      
      ! Prepare structured solution vector
      call HYPRE_StructVectorCreate(this%cfg%comm,this%hypre_box,this%hypre_sol,ierr)
      call HYPRE_StructVectorInitialize(this%hypre_sol,ierr)
      
      ! Initialize and setup the actual solver
      select case (this%method)
      case (smg)
         
         ! Create SMG solver
         call HYPRE_StructSMGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructSMGSetPrintLevel (this%hypre_solver,sprintlvl,ierr)
         call HYPRE_StructSMGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructSMGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_StructSMGSetLogging    (this%hypre_solver,1,ierr)
         
         ! Setup SMG solver
         call HYPRE_StructSMGSetup         (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case (pcg_smg)
         
         ! Create SMG preconditioner
         call HYPRE_StructSMGCreate      (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_StructSMGSetMaxIter  (this%hypre_precond,1,ierr)
         call HYPRE_StructSMGSetTol      (this%hypre_precond,0.0_WP,ierr)
         
         ! Create PCG solver
         call HYPRE_StructPCGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructPCGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPCGSetTol       (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructPCGSetLogging   (this%hypre_solver,1,ierr)
         
         ! Set SMG as preconditioner to PCG
         call HYPRE_StructPCGSetPrecond   (this%hypre_solver,0,this%hypre_precond,ierr)
         
         ! Setup PCG solver
         call HYPRE_StructPCGSetup        (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case (gmres_smg)
         
         ! Create SMG preconditioner
         call HYPRE_StructSMGCreate      (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_StructSMGSetMaxIter  (this%hypre_precond,1,ierr)
         call HYPRE_StructSMGSetTol      (this%hypre_precond,0.0_WP,ierr)
         
         ! Create GMRES solver
         call HYPRE_StructGMRESCreate     (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructGMRESSetMaxIter (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructGMRESSetTol     (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructGMRESSetLogging (this%hypre_solver,1,ierr)
         
         ! Set SMG as preconditioner to GMRES
         call HYPRE_StructGMRESSetPrecond (this%hypre_solver,0,this%hypre_precond,ierr)
         
         ! Setup GMRES solver
         call HYPRE_StructGMRESSetup      (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case (pfmg)
         
         ! Create PFMG solver
         call HYPRE_StructPFMGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         if (this%maxlevel.gt.0) call HYPRE_StructPFMGSetMaxLevels(this%hypre_solver,this%maxlevel,ierr)
         call HYPRE_StructPFMGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPFMGSetTol       (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_StructPFMGSetLogging   (this%hypre_solver,1,ierr)
         call HYPRE_StructPFMGSetPrintLevel(this%hypre_solver,sprintlvl,ierr)
         
         ! Setup PFMG solver
         call HYPRE_StructPFMGSetup        (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      
      case (pcg_pfmg)
         
         ! Create PFMG preconditioner
         call HYPRE_StructPFMGCreate      (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_StructPFMGSetMaxIter  (this%hypre_precond,1,ierr)
         call HYPRE_StructPFMGSetTol      (this%hypre_precond,0.0_WP,ierr)
         if (this%maxlevel.gt.0) call HYPRE_StructPFMGSetMaxLevels(this%hypre_precond,this%maxlevel,ierr)
         
         ! Create PCG solver
         call HYPRE_StructPCGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructPCGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPCGSetTol       (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructPCGSetLogging   (this%hypre_solver,1,ierr)
         
         ! Set PFMG as preconditioner to PCG
         call HYPRE_StructPCGSetPrecond   (this%hypre_solver,1,this%hypre_precond,ierr)
         
         ! Setup PCG solver
         call HYPRE_StructPCGSetup        (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)

      case (pcg_pfmg2)
         
         ! Create PFMG preconditioner
         call HYPRE_StructPFMGCreate      (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_StructPFMGSetMaxIter  (this%hypre_precond,1,ierr)
         call HYPRE_StructPFMGSetTol      (this%hypre_precond,0.0_WP,ierr)
         if (this%maxlevel.gt.0) call HYPRE_StructPFMGSetMaxLevels(this%hypre_precond,this%maxlevel,ierr)
         !call HYPRE_StructPFMGSetRelaxType(this%hypre_precond,2,ierr) ! Set symmetric RBGS as smoother
         call HYPRE_StructPFMGSetRAPType  (this%hypre_precond,1,ierr) ! Force use of 7-pt coarse stencil
         
         ! Create PCG solver
         call HYPRE_StructPCGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructPCGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPCGSetTol       (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructPCGSetLogging   (this%hypre_solver,1,ierr)
         
         ! Set PFMG as preconditioner to PCG
         call HYPRE_StructPCGSetPrecond   (this%hypre_solver,1,this%hypre_precond,ierr)
         
         ! Setup PCG solver
         call HYPRE_StructPCGSetup        (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case (pcg)
         
         ! Create PCG solver
         call HYPRE_StructPCGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructPCGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPCGSetTol       (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructPCGSetLogging   (this%hypre_solver,1,ierr)
         
         ! Setup PCG solver
         call HYPRE_StructPCGSetup        (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case (gmres_pfmg)
         
         ! Create PFMG preconditioner
         call HYPRE_StructPFMGCreate      (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_StructPFMGSetMaxIter  (this%hypre_precond,1,ierr)
         call HYPRE_StructPFMGSetTol      (this%hypre_precond,0.0_WP,ierr)
         if (this%maxlevel.gt.0) call HYPRE_StructPFMGSetMaxLevels(this%hypre_precond,this%maxlevel,ierr)
         
         ! Create GMRES solver
         call HYPRE_StructGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructGMRESSetLogging  (this%hypre_solver,1,ierr)
         
         ! Set PFMG as preconditioner to GMRES
         call HYPRE_StructGMRESSetPrecond  (this%hypre_solver,1,this%hypre_precond,ierr)
         
         ! Setup GMRES solver
         call HYPRE_StructGMRESSetup       (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case (gmres)
         
         ! Create GMRES solver
         call HYPRE_StructGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_StructGMRESSetLogging  (this%hypre_solver,1,ierr)
         
         ! Setup GMRES solver
         call HYPRE_StructGMRESSetup       (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
      case default
         
         ! Solver is unknown
         call die('[hypre_str setup] Unknown solution method')
         
      end select
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine hypre_str_setup
   
   
   !> Solve the linear system iteratively
   subroutine hypre_str_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(hypre_str), intent(inout) :: this
      integer :: i,j,k,ierr
      
      ! Check that setup was done
      if (.not.this%setup_done) call die('[hypre_str solve] Solver has not been setup.')
      
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      
      ! Transfer data to the structured vectors
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (abs(this%opr(1,i,j,k)).gt.0.0_WP) then
                  call HYPRE_StructVectorSetValues(this%hypre_rhs,[i,j,k],this%rhs(i,j,k),ierr)
                  call HYPRE_StructVectorSetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
               else
                  call HYPRE_StructVectorSetValues(this%hypre_rhs,[i,j,k],0.0_WP,ierr)
                  call HYPRE_StructVectorSetValues(this%hypre_sol,[i,j,k],0.0_WP,ierr)
               end if
            end do
         end do
      end do
      
      ! Call the appropriate solver
      select case (this%method)
      case (smg)
         call HYPRE_StructSMGSolve           (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructSMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructSMGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      case (pfmg)
         call HYPRE_StructPFMGSolve          (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructPFMGGetNumIteration(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructPFMGGetFinalRelativ(this%hypre_solver,this%rerr,ierr)
      case (pcg,pcg_smg,pcg_pfmg,pcg_pfmg2)
         call HYPRE_StructPCGSolve           (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructPCGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructPCGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      case (gmres,gmres_smg,gmres_pfmg)
         call HYPRE_StructGMRESSolve         (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructGMRESGetNumIteratio(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructGMRESGetFinalRelati(this%hypre_solver,this%rerr,ierr)
      end select
      
      ! Retrieve solution from structured vector
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call HYPRE_StructVectorGetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
      
   end subroutine hypre_str_solve
   
   
   !> Destroy solver - done everytime the operator changes
   subroutine hypre_str_destroy(this)
      use messager, only: die
      implicit none
      class(hypre_str), intent(inout) :: this
      integer :: ierr
      
      ! Destroy solver, operator, and rhs/sol vectors
      call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
      call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
      call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      select case (this%method)
      case (smg)
         call HYPRE_StructSMGDestroy(this%hypre_solver,ierr)
      case (pcg_smg)
         call HYPRE_StructSMGDestroy(this%hypre_precond,ierr)
         call HYPRE_StructPCGDestroy(this%hypre_solver,ierr)
      case (gmres_smg)
         call HYPRE_StructSMGDestroy(this%hypre_precond,ierr)
         call HYPRE_StructGMRESDestroy(this%hypre_solver,ierr)
      case (pfmg)
         call HYPRE_StructPFMGDestroy(this%hypre_solver,ierr)
      case (pcg)
         call HYPRE_StructPCGDestroy(this%hypre_solver,ierr)
      case (pcg_pfmg,pcg_pfmg2)
         call HYPRE_StructPFMGDestroy(this%hypre_precond,ierr)
         call HYPRE_StructPCGDestroy(this%hypre_solver,ierr)
      case (gmres)
         call HYPRE_StructGMRESDestroy(this%hypre_solver,ierr)
      case (gmres_pfmg)
         call HYPRE_StructPFMGDestroy(this%hypre_precond,ierr)
         call HYPRE_StructGMRESDestroy(this%hypre_solver,ierr)
      end select
      
      ! Set setup-flag to false
      this%setup_done=.false.
      
   end subroutine hypre_str_destroy
   
   
   !> Log hypre_str info
   subroutine hypre_str_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(hypre_str), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Structured HYPRE solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
         write(message,'(" >     method = ",i0)') this%method; call log(message)
         write(message,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit; call log(message)
         write(message,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg; call log(message)
         write(message,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg; call log(message)
      end if
   end subroutine hypre_str_log
   
   
   !> Print hypre_str info to the screen
   subroutine hypre_str_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(hypre_str), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("Structured HYPRE solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >     method = ",i0)') this%method
         write(output_unit,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit
         write(output_unit,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg
         write(output_unit,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg
      end if
   end subroutine hypre_str_print
   
   
   !> Short print of hypre_str info to the screen
   subroutine hypre_str_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(hypre_str), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("Structured HYPRE solver [",a16,"] for config [",a16,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   end subroutine hypre_str_print_short
   
   
end module hypre_str_class
