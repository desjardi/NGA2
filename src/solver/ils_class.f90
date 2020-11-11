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
   integer, parameter, public :: pcg_amg=3
   integer, parameter, public :: pcg_parasail=4
   integer, parameter, public :: gmres=5
   integer, parameter, public :: gmres_pilut=6
   integer, parameter, public :: smg=7
   integer, parameter, public :: pfmg=8
   
   
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
      
      ! Unstructured mapping
      integer,  dimension(:,:,:), allocatable :: ind     !< Unique global index
      integer :: ind_min,ind_max                         !< Local min and max indices
      integer :: ncell_,ncell                            !< Total number of local and global cells
      
      ! Private stuff for hypre
      integer(kind=8), private :: hypre_box              !< Grid
      integer(kind=8), private :: hypre_stc              !< Stencil
      integer(kind=8), private :: hypre_mat,parse_mat    !< Matrix
      integer(kind=8), private :: hypre_rhs,parse_rhs    !< Right-hand side
      integer(kind=8), private :: hypre_sol,parse_sol    !< Solution
      integer(kind=8), private :: hypre_solver           !< Solver
      integer(kind=8), private :: hypre_precond          !< Preconditioner
      
   contains
      procedure :: print=>ils_print                                   !< Print ILS object info
      procedure :: log  =>ils_log                                     !< Log ILS object info
      procedure :: init_solver                                        !< Solver initialization (at start-up)
      procedure :: update_solver                                      !< Solver update (every time the operator changes)
      procedure :: solve                                              !< Equation solver
      procedure :: solve_rbgs                                         !< Solve equation with rbgs
      procedure :: solve_hypre_amg                                    !< Solve equation with amg
      procedure :: solve_hypre_pcg                                    !< Solve equation with pcg (might be preconditioned)
      procedure :: solve_hypre_gmres                                  !< Solve equation with gmres (might be preconditioned)
      procedure :: solve_hypre_smg                                    !< Solve equation with smg
      procedure :: solve_hypre_pfmg                                   !< Solve equation with pfmg
      procedure :: prep_umap                                          !< Create unstructured mapping
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
      if (allocated(this%ind)) deallocate(this%ind)
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
   !> When creating the solver, zero values of this%opr(1,:,:,:)
   !> indicate cells that do not participate in the solver
   subroutine init_solver(this,method)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer, intent(in) :: method
      integer :: ierr,st
      integer, dimension(3) :: periodicity,offset
      
      ! Set solution method
      this%method=method
      
      ! Select appropriate solver
      select case (this%method)
      case (rbgs) ! Initialize the Red-Black Gauss-Seidel solver here
         
      case (amg)  ! Initialize the HYPRE-AMG solver here
         
         ! Allocate and prepare unstructed mapping
         call this%prep_umap()
         
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         call HYPRE_IJMatrixInitialize   (this%hypre_mat,ierr)
         ! Create a HYPRE rhs vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
         ! Create a HYPRE solution vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_sol,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
         ! Create a HYPRE AMG solver
         call HYPRE_BoomerAMGCreate        (this%hypre_solver,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_solver,3,ierr)            ! print solve info + parameters (3 from all, 0 for none)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_solver,6,ierr)            ! interpolation default is 6 (NGA=3)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_solver,6,ierr)            ! Falgout=6 (old default, NGA); PMIS=8 and HMIS=10 (recommended)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_solver,0.15_WP,ierr)    ! 0.25 is default
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_solver,8,ierr)            ! hybrid symmetric Gauss-Seidel/SOR
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_solver,this%maxit,ierr)   ! maximum nbr of iter
         call HYPRE_BoomerAMGSetTol        (this%hypre_solver,this%rcvg ,ierr)   ! convergence tolerance
         
      case (pcg_amg)  ! Initialize the HYPRE-PCG-AMG solver here
         
         ! Allocate and prepare unstructed mapping
         call this%prep_umap()
         
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         call HYPRE_IJMatrixInitialize   (this%hypre_mat,ierr)
         ! Create a HYPRE rhs vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
         ! Create a HYPRE solution vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_sol,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
         ! Create a HYPRE PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,0,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         ! Create a HYPRE AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,0,ierr)            ! print solve info + parameters (3 from all, 0 for none)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,6,ierr)            ! interpolation default is 6 (NGA=3)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,6,ierr)            ! Falgout=6 (old default, NGA); PMIS=8 and HMIS=10 (recommended)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,0.15_WP,ierr)    ! 0.25 is default
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,8,ierr)            ! hybrid symmetric Gauss-Seidel/SOR
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)            ! maximum nbr of iter
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)       ! convergence tolerance
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         ! Set AMG as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,2,this%hypre_precond,ierr)
         
      case (pcg_parasail)  ! Initialize the HYPRE-PCG-PARASAIL solver here
         
         ! Allocate and prepare unstructed mapping
         call this%prep_umap()
         
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         call HYPRE_IJMatrixInitialize   (this%hypre_mat,ierr)
         ! Create a HYPRE rhs vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
         ! Create a HYPRE solution vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_sol,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
         ! Create a HYPRE PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,0,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         ! Create a HYPRE PARASAIL preconditioner
         call HYPRE_ParaSailsCreate        (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_ParaSailsSetParams     (this%hypre_precond,0.1_WP,1,ierr)
         call HYPRE_ParaSailsSetFilter     (this%hypre_precond,0.05_WP,ierr)
         call HYPRE_ParaSailsSetSym        (this%hypre_precond,0,ierr)
         call HYPRE_ParaSailsSetLogging    (this%hypre_precond,1,ierr)
         ! Set AMG as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,4,this%hypre_precond,ierr)
         
      case (gmres)  ! Initialize the HYPRE-GMRES solver here
         
         ! Allocate and prepare unstructed mapping
         call this%prep_umap()
         
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         call HYPRE_IJMatrixInitialize   (this%hypre_mat,ierr)
         ! Create a HYPRE rhs vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
         ! Create a HYPRE solution vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_sol,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
         ! Create a HYPRE GMRES solver
         call HYPRE_ParCSRGMRESCreate    (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetKDim   (this%hypre_solver,5,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter(this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol    (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging(this%hypre_solver,1,ierr)
         
      case (gmres_pilut)  ! Initialize the HYPRE-GMRES-PILUT solver here
         
         ! Allocate and prepare unstructed mapping
         call this%prep_umap()
         
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         call HYPRE_IJMatrixInitialize   (this%hypre_mat,ierr)
         ! Create a HYPRE rhs vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_rhs,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_rhs,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_rhs,this%parse_rhs,ierr)
         ! Create a HYPRE solution vector
         call HYPRE_IJVectorCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%hypre_sol,ierr)
         call HYPRE_IJVectorSetObjectType(this%hypre_sol,hypre_ParCSR,ierr)
         call HYPRE_IJVectorInitialize   (this%hypre_sol,ierr)
         call HYPRE_IJVectorGetObject    (this%hypre_sol,this%parse_sol,ierr)
         ! Create a HYPRE GMRES solver
         call HYPRE_ParCSRGMRESCreate    (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetKDim   (this%hypre_solver,5,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter(this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol    (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging(this%hypre_solver,1,ierr)
         ! Create the PILUT preconditioner
         call HYPRE_ParCSRPilutCreate    (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_ParCSRGMRESSetPrecond(this%hypre_solver,3,this%hypre_precond,ierr)
         
      case (smg) ! Initialize the HYPRE-SMG solver here
         
         ! Create Hypre grid
         call HYPRE_StructGridCreate      (this%cfg%comm,3,this%hypre_box,ierr)
         call HYPRE_StructGridSetExtents  (this%hypre_box,[this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_],[this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_],ierr)
         periodicity=0
         if (this%cfg%xper) periodicity(1)=this%cfg%nx
         if (this%cfg%yper) periodicity(2)=this%cfg%ny
         if (this%cfg%zper) periodicity(3)=this%cfg%nz
         call HYPRE_StructGridSetPeriodic (this%hypre_box,periodicity,ierr)
         call HYPRE_StructGridAssemble    (this%hypre_box,ierr)
         ! Build Hypre stencil
         call HYPRE_StructStencilCreate   (3,this%nst,this%hypre_stc,ierr)
         do st=1,this%nst
            offset=this%stc(st,:)
            call HYPRE_StructStencilSetElement(this%hypre_stc,st-1,offset,ierr)
         end do
         ! Build Hypre matrix
         call HYPRE_StructMatrixCreate    (this%cfg%comm,this%hypre_box,this%hypre_stc,this%hypre_mat,ierr)
         call HYPRE_StructMatrixInitialize(this%hypre_mat,ierr)
         ! Prepare Hypre RHS
         call HYPRE_StructVectorCreate    (this%cfg%comm,this%hypre_box,this%hypre_rhs,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_rhs,ierr)
         ! Create Hypre solution vector
         call HYPRE_StructVectorCreate    (this%cfg%comm,this%hypre_box,this%hypre_sol,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_sol,ierr)
         ! Prepare the solver
         call HYPRE_StructSMGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructSMGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructSMGSetTol       (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_StructSMGSetLogging   (this%hypre_solver,1,ierr)
         
      case (pfmg) ! Initialize the HYPRE-PFMG solver here
         
         ! Create Hypre grid
         call HYPRE_StructGridCreate      (this%cfg%comm,3,this%hypre_box,ierr)
         call HYPRE_StructGridSetExtents  (this%hypre_box,[this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_],[this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_],ierr)
         periodicity=0
         if (this%cfg%xper) periodicity(1)=this%cfg%nx
         if (this%cfg%yper) periodicity(2)=this%cfg%ny
         if (this%cfg%zper) periodicity(3)=this%cfg%nz
         call HYPRE_StructGridSetPeriodic (this%hypre_box,periodicity,ierr)
         call HYPRE_StructGridAssemble    (this%hypre_box,ierr)
         ! Build Hypre stencil
         call HYPRE_StructStencilCreate   (3,this%nst,this%hypre_stc,ierr)
         do st=1,this%nst
            offset=this%stc(st,:)
            call HYPRE_StructStencilSetElement(this%hypre_stc,st-1,offset,ierr)
         end do
         ! Build Hypre matrix
         call HYPRE_StructMatrixCreate    (this%cfg%comm,this%hypre_box,this%hypre_stc,this%hypre_mat,ierr)
         call HYPRE_StructMatrixInitialize(this%hypre_mat,ierr)
         ! Prepare Hypre RHS
         call HYPRE_StructVectorCreate    (this%cfg%comm,this%hypre_box,this%hypre_rhs,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_rhs,ierr)
         ! Create Hypre solution vector
         call HYPRE_StructVectorCreate    (this%cfg%comm,this%hypre_box,this%hypre_sol,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_sol,ierr)
         ! Prepare the solver
         call HYPRE_StructPFMGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructPFMGSetRelaxType (this%hypre_solver,2,ierr)
         call HYPRE_StructPFMGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPFMGSetTol       (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_StructPFMGSetLogging   (this%hypre_solver,1,ierr)
         
      case default
         call die('[ils prep solver] Unknown solution method')
      end select
      
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
         allocate( row(         this%ncell_))
         allocate(ncol(         this%ncell_))
         allocate( col(this%nst*this%ncell_))
         allocate( val(this%nst*this%ncell_))
         
         ! Tranfer operator to HYPRE
         count1=0; count2=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%ind(i,j,k).gt.0) then
                     count1=count1+1
                     row (count1)=this%ind(i,j,k)
                     ncol(count1)=0
                     do st=1,this%nst
                        if (this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3)).gt.0) then
                           count2=count2+1
                           ncol(count1)=ncol(count1)+1
                           col (count2)=this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
                           val (count2)=this%opr(st,i,j,k)
                        end if
                     end do
                  end if
               end do
            end do
         end do
         call HYPRE_IJMatrixSetValues(this%hypre_mat,this%ncell_,ncol,row,col,val,ierr)
         call HYPRE_IJMatrixAssemble (this%hypre_mat,ierr)
         call HYPRE_IJMatrixGetObject(this%hypre_mat,this%parse_mat,ierr)
         
         ! Setup AMG solver
         call HYPRE_BoomerAMGSetup(this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
         ! Deallocate
         deallocate(row,ncol,col,val)
         
      case (pcg_amg,pcg_parasail)  ! Prepare the HYPRE-PCG solver here
         
         ! Allocate storage for transfer
         allocate( row(         this%ncell_))
         allocate(ncol(         this%ncell_))
         allocate( col(this%nst*this%ncell_))
         allocate( val(this%nst*this%ncell_))
         
         ! Tranfer operator to HYPRE
         count1=0; count2=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%ind(i,j,k).gt.0) then
                     count1=count1+1
                     row (count1)=this%ind(i,j,k)
                     ncol(count1)=0
                     do st=1,this%nst
                        if (this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3)).gt.0) then
                           count2=count2+1
                           ncol(count1)=ncol(count1)+1
                           col (count2)=this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
                           val (count2)=this%opr(st,i,j,k)
                        end if
                     end do
                  end if
               end do
            end do
         end do
         call HYPRE_IJMatrixSetValues(this%hypre_mat,this%ncell_,ncol,row,col,val,ierr)
         call HYPRE_IJMatrixAssemble (this%hypre_mat,ierr)
         call HYPRE_IJMatrixGetObject(this%hypre_mat,this%parse_mat,ierr)
         
         ! Setup PCG-AMG solver
         call HYPRE_ParCSRPCGSetup(this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
         ! Deallocate
         deallocate(row,ncol,col,val)
         
      case (gmres,gmres_pilut)  ! Prepare the HYPRE-GMRES solver here
         
         ! Allocate storage for transfer
         allocate( row(         this%ncell_))
         allocate(ncol(         this%ncell_))
         allocate( col(this%nst*this%ncell_))
         allocate( val(this%nst*this%ncell_))
         
         ! Tranfer operator to HYPRE
         count1=0; count2=0
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (this%ind(i,j,k).gt.0) then
                     count1=count1+1
                     row (count1)=this%ind(i,j,k)
                     ncol(count1)=0
                     do st=1,this%nst
                        if (this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3)).gt.0) then
                           count2=count2+1
                           ncol(count1)=ncol(count1)+1
                           col (count2)=this%ind(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
                           val (count2)=this%opr(st,i,j,k)
                        end if
                     end do
                  end if
               end do
            end do
         end do
         call HYPRE_IJMatrixSetValues(this%hypre_mat,this%ncell_,ncol,row,col,val,ierr)
         call HYPRE_IJMatrixAssemble (this%hypre_mat,ierr)
         call HYPRE_IJMatrixGetObject(this%hypre_mat,this%parse_mat,ierr)
         
         ! Setup GMRES solver
         call HYPRE_ParCSRGMRESSetup(this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         
         ! Deallocate
         deallocate(row,ncol,col,val)
         
      case (smg)  ! Prepare the HYPRE-SMG solver here
         
         ! Prepare local storage
         allocate(row(1:this%nst))
         do st=1,this%nst
            row(st)=st-1
         end do
         allocate(val(1:this%nst))
         
         ! Transfer the operator
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  val=this%opr(:,i,j,k)
                  call HYPRE_StructMatrixSetValues(this%hypre_mat,[i,j,k],this%nst,row,val,ierr)
               end do
            end do
         end do
         call HYPRE_StructMatrixAssemble(this%hypre_mat,ierr)
         
         ! Setup the solver
         call HYPRE_StructSMGSetup(this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
         ! Deallocate
         deallocate(row,val)
         
      case (pfmg)  ! Prepare the HYPRE-PFMG solver here
         
         ! Prepare local storage
         allocate(row(1:this%nst))
         do st=1,this%nst
            row(st)=st-1
         end do
         allocate(val(1:this%nst))
         
         ! Transfer the operator
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  val=this%opr(:,i,j,k)
                  call HYPRE_StructMatrixSetValues(this%hypre_mat,[i,j,k],this%nst,row,val,ierr)
               end do
            end do
         end do
         call HYPRE_StructMatrixAssemble(this%hypre_mat,ierr)
         
         ! Setup the solver
         call HYPRE_StructPFMGSetup(this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         
         ! Deallocate
         deallocate(row,val)
         
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
      case (pcg_amg,pcg_parasail)
         call this%solve_hypre_pcg()
      case (gmres,gmres_pilut)
         call this%solve_hypre_gmres()
      case (smg)
         call this%solve_hypre_smg()
      case (pfmg)
         call this%solve_hypre_pfmg()
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
      allocate(rhs(this%ncell_))
      allocate(sol(this%ncell_))
      allocate(ind(this%ncell_))
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).gt.0) then
                  count=count+1
                  ind(count)=this%ind(i,j,k)
                  rhs(count)=this%rhs(i,j,k)
                  sol(count)=this%sol(i,j,k)
               end if
            end do
         end do
      end do
      call HYPRE_IJVectorSetValues(this%hypre_rhs,this%ncell_,ind,rhs,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_rhs,ierr)
      call HYPRE_IJVectorSetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_sol,ierr)
      
      ! Call the solver
      call HYPRE_BoomerAMGSolve           (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      call HYPRE_BoomerAMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
      call HYPRE_BoomerAMGGetFinalReltvRes(this%hypre_solver,this%rerr,ierr)
      
      ! Get the solution back
      call HYPRE_IJVectorGetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).gt.0) then
                  count=count+1
                  this%sol(i,j,k)=sol(count)
               end if
            end do
         end do
      end do
      
      ! Deallocate
      deallocate(rhs,sol,ind)
      
   end subroutine solve_hypre_amg
   
   
   !> Solve the linear system using hypre_pcg
   subroutine solve_hypre_pcg(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count,ierr
      integer,  dimension(:), allocatable :: ind
      real(WP), dimension(:), allocatable :: rhs,sol
      
      ! Transfer the rhs and initial guess to hypre
      allocate(rhs(this%ncell_))
      allocate(sol(this%ncell_))
      allocate(ind(this%ncell_))
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).gt.0) then
                  count=count+1
                  ind(count)=this%ind(i,j,k)
                  rhs(count)=this%rhs(i,j,k)
                  sol(count)=this%sol(i,j,k)
               end if
            end do
         end do
      end do
      call HYPRE_IJVectorSetValues(this%hypre_rhs,this%ncell_,ind,rhs,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_rhs,ierr)
      call HYPRE_IJVectorSetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_sol,ierr)
      
      ! Call the solver
      call HYPRE_ParCSRPCGSolve           (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      call HYPRE_ParCSRPCGGetNumIterations(this%hypre_solver,this%it  ,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      
      ! Get the solution back
      call HYPRE_IJVectorGetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).gt.0) then
                  count=count+1
                  this%sol(i,j,k)=sol(count)
               end if
            end do
         end do
      end do
      
      ! Deallocate
      deallocate(rhs,sol,ind)
      
   end subroutine solve_hypre_pcg
   
   
   !> Solve the linear system using hypre_gmres
   subroutine solve_hypre_gmres(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count,ierr
      integer,  dimension(:), allocatable :: ind
      real(WP), dimension(:), allocatable :: rhs,sol
      
      ! Transfer the rhs and initial guess to hypre
      allocate(rhs(this%ncell_))
      allocate(sol(this%ncell_))
      allocate(ind(this%ncell_))
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).gt.0) then
                  count=count+1
                  ind(count)=this%ind(i,j,k)
                  rhs(count)=this%rhs(i,j,k)
                  sol(count)=this%sol(i,j,k)
               end if
            end do
         end do
      end do
      call HYPRE_IJVectorSetValues(this%hypre_rhs,this%ncell_,ind,rhs,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_rhs,ierr)
      call HYPRE_IJVectorSetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      call HYPRE_IJVectorAssemble (this%hypre_sol,ierr)
      
      ! Call the solver
      call HYPRE_ParCSRGMRESSolve         (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      call HYPRE_ParCSRGMRESGetNumIteratio(this%hypre_solver,this%it  ,ierr)
      call HYPRE_ParCSRGMRESGetFinalRelati(this%hypre_solver,this%rerr,ierr)
      
      ! Get the solution back
      call HYPRE_IJVectorGetValues(this%hypre_sol,this%ncell_,ind,sol,ierr)
      count=0
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%ind(i,j,k).gt.0) then
                  count=count+1
                  this%sol(i,j,k)=sol(count)
               end if
            end do
         end do
      end do
      
      ! Deallocate
      deallocate(rhs,sol,ind)
      
   end subroutine solve_hypre_gmres
   
   
   !> Solve the linear system using hypre_smg
   subroutine solve_hypre_smg(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,ierr
      
      ! Transfer the rhs and initial guess to hypre
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call HYPRE_StructVectorSetValues(this%hypre_rhs,[i,j,k],this%rhs(i,j,k),ierr)
               call HYPRE_StructVectorSetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      call HYPRE_StructVectorAssemble(this%hypre_rhs,ierr)
      call HYPRE_StructVectorAssemble(this%hypre_sol,ierr)
      
      ! Solve
      call HYPRE_StructSMGSolve           (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      call HYPRE_StructSMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
      call HYPRE_StructSMGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      
      ! Get the solution back
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call HYPRE_StructVectorGetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      
   end subroutine solve_hypre_smg
   
   
   !> Solve the linear system using hypre_pfmg
   subroutine solve_hypre_pfmg(this)
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,ierr
      
      ! Transfer the rhs and initial guess to hypre
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call HYPRE_StructVectorSetValues(this%hypre_rhs,[i,j,k],this%rhs(i,j,k),ierr)
               call HYPRE_StructVectorSetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      call HYPRE_StructVectorAssemble(this%hypre_rhs,ierr)
      call HYPRE_StructVectorAssemble(this%hypre_sol,ierr)
      
      ! Solve
      call HYPRE_StructPFMGSolve          (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      call HYPRE_StructPFMGGetNumIteration(this%hypre_solver,this%it  ,ierr)
      call HYPRE_StructPFMGGetFinalRelativ(this%hypre_solver,this%rerr,ierr)
      
      ! Get the solution back
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               call HYPRE_StructVectorGetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
            end do
         end do
      end do
      
   end subroutine solve_hypre_pfmg
   
   
   !> Creation of an unstructured mapping
   subroutine prep_umap(this)
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      implicit none
      class(ils), intent(inout) :: this
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
                  count=count+1
                  this%ind(i,j,k)=count
               end if
            end do
         end do
      end do
      
      ! Take care of periodicity and domain decomposition
      call this%cfg%sync(this%ind)
      
      ! Get local min/max
      this%ind_min=1
      if (this%cfg%rank.gt.0) this%ind_min=ncell_per_proc(this%cfg%rank)+1
      this%ind_max=ncell_per_proc(this%cfg%rank+1)
      
      ! Deallocate
      deallocate(ncell_per_proc)
      
   end subroutine prep_umap
   
   
   !> Log ILS info
   subroutine ils_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(ils), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Iterative Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
         write(message,'(" >     method = ",i0)') this%method; call log(message)
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
         write(output_unit,'(" >     method = ",i0)') this%method
         write(output_unit,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit
         write(output_unit,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg
         write(output_unit,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg
      end if
   end subroutine ils_print
   
   
end module ils_class
