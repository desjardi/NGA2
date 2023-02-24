!> Unstructured hypre linear solver concept is
!> defined by extension of the linsol class
module hypre_uns_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   use linsol_class, only: linsol
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: hypre_uns
   
   
   ! List of known available methods
   integer, parameter, public ::          amg=1
   integer, parameter, public ::      pcg_amg=2
   integer, parameter, public ::    gmres_amg=3
   integer, parameter, public :: bicgstab_amg=4
   integer, parameter, public :: pcg         =5
   integer, parameter, public :: pcg_parasail=6
   integer, parameter, public ::  gmres      =7
   integer, parameter, public ::  gmres_pilut=8
   
   
   ! List of key solver parameters
   integer , parameter :: gmres_kdim=5                 !< Number of basis vectors between restarts
   real(WP), parameter :: amg_strong_threshold=0.25_WP !< Coarsening parameter (default is 0.25, 0.5 recommended in 3D)
   integer , parameter :: amg_coarsen_type=8           !< Falgout=6 (old default); PMIS=8 and HMIS=10 (recommended)
   integer , parameter :: amg_interp=6                 !< 6=extended classical modified interpolation (default); 8=standard interpolation
   integer , parameter :: amg_relax=6                  !< 6=Hybrid symmetric Gauss-Seidel (default); 8=symmetric L1-Gauss-Seidel; 0=Weighted Jacobi
   integer , parameter :: amg_relax_coarse=6           !< 9=Gauss elim; 99=GE w/ pivot; may be good to use same as above?
   integer , parameter :: sprintlvl=0                  !< Solver  printing: 0=none (default); 3=init and cvg history
   integer , parameter :: pprintlvl=0                  !< Precond printing: 0=none (default); 3=init and cvg history
   
   
   ! Hypre-related storage parameter
   integer(kind=8), parameter :: hypre_ParCSR=5555
   
   
   !> hypre_uns object definition
   type, extends(linsol) :: hypre_uns
      
      ! For multigrid solvers, maxlevel parameter is needed
      integer :: maxlevel                                !< Maximum number of multigrid levels
      
      ! Unstructured mapping
      integer,  dimension(:,:,:), allocatable :: ind     !< Unique global index
      integer :: ind_min,ind_max                         !< Local min and max indices
      integer :: ncell_,ncell                            !< Total number of local and global cells
      
      ! Private stuff for hypre
      integer(kind=8), private :: hypre_mat,parse_mat    !< Matrix
      integer(kind=8), private :: hypre_rhs,parse_rhs    !< Right-hand side
      integer(kind=8), private :: hypre_sol,parse_sol    !< Solution
      integer(kind=8), private :: hypre_solver           !< Solver
      integer(kind=8), private :: hypre_precond          !< Preconditioner
      
   contains
      
      procedure :: print_short=>hypre_uns_print_short    !< One-line printing of solver status
      procedure :: print=>hypre_uns_print                !< Long-form printing of solver status
      procedure :: log=>hypre_uns_log                    !< Long-form logging of solver status
      procedure :: init=>hypre_uns_init                  !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>hypre_uns_setup                !< Solver setup (every time the operator changes)
      procedure :: solve=>hypre_uns_solve                !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>hypre_uns_destroy            !< Solver destruction (every time the operator changes)
      
      procedure, private :: prep_umap                    !< Create unstructured mapping
      
   end type hypre_uns
   
   
   !> Declare hypre_uns constructor
   interface hypre_uns
      procedure hypre_uns_from_args
   end interface hypre_uns
   
   
contains
   
   
   !> Constructor for an hypre_uns object
   function hypre_uns_from_args(cfg,name,method,nst) result(self)
      use messager, only: die
      implicit none
      type(hypre_uns) :: self
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
      
   end function hypre_uns_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine hypre_uns_init(this)
      use messager, only: die
      implicit none
      class(hypre_uns), intent(inout) :: this
      integer :: ierr,st,stx1,stx2,sty1,sty2,stz1,stz2
      
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
      
      ! These use HYPRE's IJ environment, which requires that we allocate and prepare an unstructed mapping
      call this%prep_umap()
      
   end subroutine hypre_uns_init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine hypre_uns_setup(this)
      use messager, only: die
      implicit none
      class(hypre_uns), intent(inout) :: this
      integer :: i,j,k,st,count1,count2,ierr,nn
      integer,  dimension(:), allocatable :: row,ncol,mycol,col
      real(WP), dimension(:), allocatable :: myval,val
      integer,  dimension(:), allocatable :: sizes
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      ! Create an IJ matrix
      call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
      call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
      allocate(sizes(this%ind_min:this%ind_max)); sizes=this%nst
      call HYPRE_IJMatrixSetRowSizes  (this%hypre_mat,sizes,ierr)
      deallocate(sizes)
      ! Allocate storage for transfer
      allocate(mycol(this%nst))
      allocate(myval(this%nst))
      allocate( row(         this%ncell_))
      allocate(ncol(         this%ncell_))
      allocate( col(this%nst*this%ncell_))
      allocate( val(this%nst*this%ncell_))
      ! Transfer operator to HYPRE
      count1=0; count2=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).ge.0) then
                  ! Prepare to add the row
                  count1=count1+1
                  row (count1)=this%ind(i,j,k)
                  ncol(count1)=0
                  ! Copy all columns for processing
                  do st=1,this%nst
                     mycol(st)=this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
                     myval(st)=this%opr(st,i,j,k)
                  end do
                  ! Compress columns to avoid redundancy
                  do st=1,this%nst
                     do nn=1,st-1
                        if (mycol(st).eq.mycol(nn).and.mycol(st).ge.0) then
                           mycol(st)=0
                           myval(nn)=myval(nn)+myval(st)
                        end if
                     end do
                  end do
                  ! Finally transfer to operator
                  do st=1,this%nst
                     if (mycol(st).ge.0) then
                        count2=count2+1
                        ncol(count1)=ncol(count1)+1
                        col (count2)=mycol(st)
                        val (count2)=myval(st)
                     end if
                  end do
               end if
            end do
         end do
      end do
      call HYPRE_IJMatrixInitialize(this%hypre_mat,ierr)
      call HYPRE_IJMatrixSetValues (this%hypre_mat,this%ncell_,ncol,row,col,val,ierr)
      call HYPRE_IJMatrixAssemble  (this%hypre_mat,ierr)
      call HYPRE_IJMatrixGetObject (this%hypre_mat,this%parse_mat,ierr)
      
      ! Deallocate storage
      deallocate(row,ncol,mycol,col,myval,val)
      
      ! For debugging - output the operator to a file
      !call HYPRE_IJMatrixPrint(this%hypre_mat,'ij_mat'//char(0),ierr)
      !call die('Matrix was printed out')
      
      ! Create an IJ vector for rhs
      call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_rhs,ierr)
      call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
      call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
      call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
      
      ! Create an IJ vector for solution
      call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_sol,ierr)
      call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
      call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
      call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
      
      ! Initialize and setup the actual solver
      select case (this%method)
      case (amg)
         
         ! Create AMG solver
         call HYPRE_BoomerAMGCreate        (this%hypre_solver,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_solver,sprintlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_solver,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_solver,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_solver,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_solver,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_solver,amg_relax_coarse,3,ierr)
         !if (this%maxlevel.gt.0) call HYPRE_BoomerAMGSetMaxLevels(this%hypre_solver,this%maxlevel,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         
         ! Setup AMG solver
         call HYPRE_BoomerAMGSetup         (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case (pcg_amg)
         
         ! Create PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         
         ! Create AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,pprintlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_precond,amg_relax_coarse,3,ierr)
         !if (this%maxlevel.gt.0) call HYPRE_BoomerAMGSetMaxLevels(this%hypre_precond,this%maxlevel,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         
         ! Set AMG as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,2,this%hypre_precond,ierr)
         
         ! Setup PCG solver
         call HYPRE_ParCSRPCGSetup         (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case (pcg_parasail)
         
         ! Create PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         
         ! Create PARASAIL preconditioner
         call HYPRE_ParaSailsCreate        (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_ParaSailsSetParams     (this%hypre_precond,0.1_WP,1,ierr)
         call HYPRE_ParaSailsSetFilter     (this%hypre_precond,0.05_WP,ierr)
         call HYPRE_ParaSailsSetSym        (this%hypre_precond,0,ierr)
         call HYPRE_ParaSailsSetLogging    (this%hypre_precond,1,ierr)
         
         ! Set PARASAIL as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,4,this%hypre_precond,ierr)
         
         ! Setup PCG solver
         call HYPRE_ParCSRPCGSetup         (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case (pcg)
         
         ! Create PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         
         ! Set Diagonal Scaling as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,1,this%hypre_precond,ierr)
         
         ! Setup PCG solver
         call HYPRE_ParCSRPCGSetup         (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case (gmres)
         
         ! Create GMRES solver
         call HYPRE_ParCSRGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetPrintLevel(this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRGMRESSetKDim     (this%hypre_solver,gmres_kdim,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging  (this%hypre_solver,1,ierr)
         
         ! Set Diagonal Scaling as preconditioner to GMRES
         call HYPRE_ParCSRGMRESSetPrecond  (this%hypre_solver,1,this%hypre_precond,ierr)
         
         ! Setup GMRES solver
         call HYPRE_ParCSRGMRESSetup       (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case (gmres_amg)
         
         ! Create GMRES solver
         call HYPRE_ParCSRGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetPrintLevel(this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRGMRESSetKDim     (this%hypre_solver,gmres_kdim,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging  (this%hypre_solver,1,ierr)
         
         ! Create AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,pprintlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_precond,amg_relax_coarse,3,ierr)
         !if (this%maxlevel.gt.0) call HYPRE_BoomerAMGSetMaxLevels(this%hypre_precond,this%maxlevel,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         
         ! Set AMG as preconditioner to GMRES
         call HYPRE_ParCSRGMRESSetPrecond  (this%hypre_solver,2,this%hypre_precond,ierr)
         
         ! Setup GMRES solver
         call HYPRE_ParCSRGMRESSetup       (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
         
      case (gmres_pilut)
         
         ! Create GMRES solver
         call HYPRE_ParCSRGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetPrintLevel(this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRGMRESSetKDim     (this%hypre_solver,gmres_kdim,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging  (this%hypre_solver,1,ierr)
         
         ! Create PILUT preconditioner
         call HYPRE_ParCSRPilutCreate      (this%cfg%comm,this%hypre_precond,ierr)
         
         ! Set PILUT as preconditioner to GMRES
         call HYPRE_ParCSRGMRESSetPrecond  (this%hypre_solver,3,this%hypre_precond,ierr)
         
         ! Setup GMRES solver
         call HYPRE_ParCSRGMRESSetup       (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case (bicgstab_amg)
         
         ! Create BiCGstab solver
         call HYPRE_ParCSRBiCGSTABCreate    (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRBiCGSTABSetPrintLev(this%hypre_solver,sprintlvl,ierr)
         call HYPRE_ParCSRBiCGSTABSetMaxIter(this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRBiCGSTABSetTol    (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRBiCGSTABSetLogging(this%hypre_solver,1,ierr)
         
         ! Create AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,pprintlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_precond,amg_relax_coarse,3,ierr)
         !if (this%maxlevel.gt.0) call HYPRE_BoomerAMGSetMaxLevels(this%hypre_precond,this%maxlevel,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         
         ! Set AMG as preconditioner to BiCGstab
         call HYPRE_ParCSRBiCGSTABSetPrecond(this%hypre_solver,2,this%hypre_precond,ierr)
         
         ! Setup BiCGstab solver
         call HYPRE_ParCSRBiCGSTABSetup    (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
      case default
         
         ! Solver is unknown
         call die('[hypre_uns setup] Unknown solution method')
         
      end select
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine hypre_uns_setup
   
   
   !> Solve the linear system iteratively
   subroutine hypre_uns_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(hypre_uns), intent(inout) :: this
      integer :: i,j,k,count,ierr
      integer,  dimension(:), allocatable :: ind
      real(WP), dimension(:), allocatable :: rhs,sol
      
      ! Check that setup was done
      if (.not.this%setup_done) call die('[hypre_uns solve] Solver has not been setup.')
      
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      
      ! Transfer data to the IJ vectors
      allocate(rhs(this%ncell_))
      allocate(sol(this%ncell_))
      allocate(ind(this%ncell_))
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).ge.0) then
                  count=count+1
                  ind(count)=this%ind(i,j,k)
                  rhs(count)=this%rhs(i,j,k)
                  sol(count)=this%sol(i,j,k)
               end if
            end do
         end do
      end do
      call HYPRE_IJVectorSetValues(this%hypre_rhs,this%ncell_,ind,rhs,ierr)
      call HYPRE_IJVectorSetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      
      ! Call the appropriate solver
      select case (this%method)
      case (amg)
         call HYPRE_BoomerAMGSolve           (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         call HYPRE_BoomerAMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(this%hypre_solver,this%rerr,ierr)
      case (pcg,pcg_amg,pcg_parasail)
         call HYPRE_ParCSRPCGSolve           (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         call HYPRE_ParCSRPCGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_ParCSRPCGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      case (gmres,gmres_pilut,gmres_amg)
         call HYPRE_ParCSRGMRESSolve         (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         call HYPRE_ParCSRGMRESGetNumIteratio(this%hypre_solver,this%it  ,ierr)
         call HYPRE_ParCSRGMRESGetFinalRelati(this%hypre_solver,this%rerr,ierr)
      case (bicgstab_amg)
         call HYPRE_ParCSRBiCGSTABSolve      (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         call HYPRE_ParCSRBiCGSTABGetNumIter (this%hypre_solver,this%it  ,ierr)
         call HYPRE_ParCSRBiCGSTABGetFinalRel(this%hypre_solver,this%rerr,ierr)
      end select
      
      ! Retrieve solution from IJ vector
      call HYPRE_IJVectorGetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).ge.0) then
                  count=count+1
                  this%sol(i,j,k)=sol(count)
               end if
            end do
         end do
      end do
      deallocate(rhs,sol,ind)
      
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
      
   end subroutine hypre_uns_solve
   
   
   !> Destroy solver - done everytime the operator changes
   subroutine hypre_uns_destroy(this)
      use messager, only: die
      implicit none
      class(hypre_uns), intent(inout) :: this
      integer :: ierr
      
      ! Destroy solver, operator, and rhs/sol vectors
      call HYPRE_IJMatrixDestroy(this%hypre_mat,ierr)
      call HYPRE_IJVectorDestroy(this%hypre_rhs,ierr)
      call HYPRE_IJVectorDestroy(this%hypre_sol,ierr)
      select case (this%method)
      case (amg)
         call HYPRE_BoomerAMGDestroy(this%hypre_solver,ierr)
      case (pcg_amg)
         call HYPRE_BoomerAMGDestroy(this%hypre_precond,ierr)
         call HYPRE_ParCSRPCGDestroy(this%hypre_solver,ierr)
      case (pcg_parasail)
         call HYPRE_ParaSailsDestroy(this%hypre_precond,ierr)
         call HYPRE_ParCSRPCGDestroy(this%hypre_solver,ierr)
      case (pcg)
         call HYPRE_ParCSRPCGDestroy(this%hypre_solver,ierr)
      case (gmres)
         call HYPRE_ParCSRGMRESDestroy(this%hypre_solver,ierr)
      case (gmres_amg)
         call HYPRE_BoomerAMGDestroy(this%hypre_precond,ierr)
         call HYPRE_ParCSRGMRESDestroy(this%hypre_solver,ierr)
      case (gmres_pilut)
         call HYPRE_ParCSRPilutDestroy(this%hypre_precond,ierr)
         call HYPRE_ParCSRGMRESDestroy(this%hypre_solver,ierr)
      case (bicgstab_amg)
         call HYPRE_BoomerAMGDestroy(this%hypre_precond,ierr)
         call HYPRE_ParCSRBiCGSTABDestroy(this%hypre_solver,ierr)
      end select
      
      ! Set setup-flag to false
      this%setup_done=.false.
      
   end subroutine hypre_uns_destroy
   
   
   !> Creation of an unstructured mapping
   subroutine prep_umap(this)
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      implicit none
      class(hypre_uns), intent(inout) :: this
      integer :: i,j,k,ierr,count
      integer, dimension(:), allocatable :: ncell_per_proc
      
      ! Dump any existing mapping and recreate it
      this%ncell =0
      this%ncell_=0
      this%ind_min=0
      this%ind_max=0
      if (.not.allocated(this%ind)) allocate(this%ind(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      this%ind=-1
      
      ! Count number of active cells
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (abs(this%opr(1,i,j,k)).gt.10.0_WP*epsilon(1.0_WP)) this%ncell_=this%ncell_+1
            end do
         end do
      end do
      call MPI_ALLREDUCE(this%ncell_,this%ncell,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
      
      ! Create an array with ncell_ per cpu
      allocate(ncell_per_proc(this%cfg%nproc))
      call MPI_ALLGATHER(this%ncell_,1,MPI_INTEGER,ncell_per_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
      do i=2,this%cfg%nproc
         ncell_per_proc(i)=ncell_per_proc(i)+ncell_per_proc(i-1)
      end do
      
      ! Assign unique global index to all non-empty cells
      count=0
      if (this%cfg%rank.gt.0) count=ncell_per_proc(this%cfg%rank)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (abs(this%opr(1,i,j,k)).gt.10.0_WP*epsilon(1.0_WP)) then
                  this%ind(i,j,k)=count
                  count=count+1
               end if
            end do
         end do
      end do
      
      ! Take care of periodicity and domain decomposition
      call this%cfg%sync(this%ind)
      
      ! Get local min/max
      this%ind_min=0
      if (this%cfg%rank.gt.0) this%ind_min=ncell_per_proc(this%cfg%rank)
      this%ind_max=ncell_per_proc(this%cfg%rank+1)-1
      
      ! Deallocate
      deallocate(ncell_per_proc)
      
   end subroutine prep_umap
   
   
   !> Log hypre_uns info
   subroutine hypre_uns_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(hypre_uns), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Unstructured Hypre Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
         write(message,'(" >     method = ",i0)') this%method; call log(message)
         write(message,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit; call log(message)
         write(message,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg; call log(message)
         write(message,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg; call log(message)
      end if
   end subroutine hypre_uns_log
   
   
   !> Print hypre_uns info to the screen
   subroutine hypre_uns_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(hypre_uns), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("Unstructured Hypre Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >     method = ",i0)') this%method
         write(output_unit,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit
         write(output_unit,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg
         write(output_unit,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg
      end if
   end subroutine hypre_uns_print
   
   
   !> Short print of hypre_uns info to the screen
   subroutine hypre_uns_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(hypre_uns), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("Unstructured Hypre Linear Solver [",a16,"] for config [",a16,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   end subroutine hypre_uns_print_short
   
   
end module hypre_uns_class
