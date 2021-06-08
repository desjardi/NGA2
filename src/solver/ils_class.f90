!> Iterative linear solver concept is defined here:
!> given a config, it provides parallel solvers
module ils_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   use bbmg_class,   only: bbmg
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ils
   
   
   ! List of known available methods
   integer, parameter, public ::         bbox= 1
   integer, parameter, public ::     pcg_bbox= 2
   integer, parameter, public ::          amg= 3
   integer, parameter, public ::      pcg_amg= 4
   integer, parameter, public ::    gmres_amg= 5
   integer, parameter, public :: bicgstab_amg= 6
   integer, parameter, public :: pcg         = 7
   integer, parameter, public :: pcg_parasail= 8
   integer, parameter, public ::  gmres      = 9
   integer, parameter, public ::  gmres_pilut=10
   integer, parameter, public ::         pfmg=11
   integer, parameter, public ::     pcg_pfmg=12
   integer, parameter, public ::   gmres_pfmg=13
   integer, parameter, public ::          smg=14
   integer, parameter, public ::      pcg_smg=15
   integer, parameter, public ::    gmres_smg=16
   
   
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
   
   !> ILS object definition
   type :: ils
      ! An iterative linear solver works for a config
      type(config), pointer :: cfg                                    !< Config for the ILS
      ! An ILS has a name
      character(len=str_medium) :: name                               !< Name of solver
      ! An ILS has a solution method
      integer  :: method                                              !< Solution method
      ! An ILS has a stencil size
      integer  :: nst                                                 !< Stencil size in 3D
      integer, dimension(:,:),   allocatable :: stc                   !< Stencil map in 3D (from stencil entry to i/j/k shift)
      integer, dimension(:,:,:), allocatable :: stmap                 !< Inverse stencil map in 3D (from i,j,k shift to stencil entry)
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
      
      ! Is the solver setup?
      logical :: setup_done                                           !< Check whether the solver has been setup
      
      ! For multigrid solvers, maxlevel parameter is needed
      integer :: maxlevel                                             !< Maximum number of multigrid levels
      
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
      
      ! Private BBMG solver
      type(bbmg), private :: bbox                        !< In-house Black-Box Multi-Grid solver
      
   contains
      procedure :: print_short=>ils_print_short                       !< Short print of ILS object info
      procedure :: print=>ils_print                                   !< Print ILS object info
      procedure :: log  =>ils_log                                     !< Log ILS object info
      procedure :: init                                               !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup                                              !< Solver setup (every time the operator changes)
      procedure :: solve                                              !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy                                            !< Solver destruction (every time the operator changes)
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
      self%maxlevel=0
      self%setup_done=.false.
      
      ! Set up stencil size and map
      self%nst=7; if (present(nst)) self%nst=nst
      allocate(self%stc(1:self%nst,1:3)); self%stc=0
      
      ! Allocate operator, rhs, and sol arrays
      allocate(self%opr(self%nst,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%opr=0.0_WP
      allocate(self%rhs(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rhs=0.0_WP
      allocate(self%sol(         self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%sol=0.0_WP
      
   end function ils_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine init(this,method)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer, intent(in) :: method
      integer :: ierr,st,stx1,stx2,sty1,sty2,stz1,stz2
      integer, dimension(3) :: periodicity,offset
      
      ! Set solution method
      this%method=method
      
      ! From the provided stencil, generate an inverse map
      stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1))
      sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2))
      stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3))
      allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
      do st=1,this%nst
         this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
      end do
      
      ! Initialize grid and stencil - this also catches incorrect solvers
      select case (this%method)
      case (bbox,pcg_bbox)
         ! For simplicity, directly initialize an in-house black-box multi-grid solver
         ! We'll need to have a proper destructor for bbox in order to adhere to ILS' intended use scenario...
         this%bbox=bbmg(pg=this%cfg,name=trim(this%name),nst=[3,3,3])
         ! Adjust parameters
         this%bbox%max_ite=this%maxit
         this%bbox%max_res=this%rcvg
         this%bbox%relax_pre =2
         this%bbox%relax_post=2
         this%bbox%use_direct_solve=.true.
         this%bbox%ncell_coarsest=300
         if (this%method.eq.    bbox) this%bbox%use_krylov=.false.
         if (this%method.eq.pcg_bbox) this%bbox%use_krylov=.true.
         ! Initialize solver
         call this%bbox%initialize()
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg,pcg)
         ! Initialize HYPRE
         call HYPRE_Init(ierr)
         ! These use HYPRE's IJ environment, which requires that we allocate and prepare an unstructed mapping
         call this%prep_umap()
      case (smg,pcg_smg,gmres_smg,pfmg,pcg_pfmg,gmres_pfmg)
         ! Initialize HYPRE
         call HYPRE_Init(ierr)
         ! These use HYPRE's structured environment, which requires that we create a HYPRE grid and stencil
         call HYPRE_StructGridCreate     (this%cfg%comm,3,this%hypre_box,ierr)
         call HYPRE_StructGridSetExtents (this%hypre_box,[this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_],[this%cfg%imax_,this%cfg%jmax_,this%cfg%kmax_],ierr)
         periodicity=0
         if (this%cfg%xper) periodicity(1)=this%cfg%nx
         if (this%cfg%yper) periodicity(2)=this%cfg%ny
         if (this%cfg%zper) periodicity(3)=this%cfg%nz
         call HYPRE_StructGridSetPeriodic(this%hypre_box,periodicity,ierr)
         call HYPRE_StructGridAssemble   (this%hypre_box,ierr)
         ! Build Hypre stencil
         call HYPRE_StructStencilCreate  (3,this%nst,this%hypre_stc,ierr)
         do st=1,this%nst
            offset=this%stc(st,:)
            call HYPRE_StructStencilSetElement(this%hypre_stc,st-1,offset,ierr)
         end do
      case default
         call die('[ils init] Unknown solution method')
      end select
      
   end subroutine init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine setup(this)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,st,count1,count2,ierr,nn
      integer,  dimension(:), allocatable :: row,ncol,mycol,col
      real(WP), dimension(:), allocatable :: myval,val
      integer,  dimension(:), allocatable :: sizes
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      ! Create the operator and rhs/sol vectors
      select case (this%method)
      case (bbox,pcg_bbox)
         
         ! Storage for matrix/vectors is already available within bbox, so just copy operator
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%bbox%lvl(1)%opr(:,:,:,i,j,k)=0.0_WP
                  do st=1,this%nst
                     this%bbox%lvl(1)%opr(this%stc(st,1),this%stc(st,2),this%stc(st,3),i,j,k)=this%bbox%lvl(1)%opr(this%stc(st,1),this%stc(st,2),this%stc(st,3),i,j,k)+this%opr(st,i,j,k)
                  end do
               end do
            end do
         end do
         
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg,pcg)
         
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
         
      case (smg,pcg_smg,gmres_smg,pfmg,pcg_pfmg,gmres_pfmg)
         
         ! Create a structured matrix
         call HYPRE_StructMatrixCreate    (this%cfg%comm,this%hypre_box,this%hypre_stc,this%hypre_mat,ierr)
         call HYPRE_StructMatrixInitialize(this%hypre_mat,ierr)
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
         call HYPRE_StructVectorCreate    (this%cfg%comm,this%hypre_box,this%hypre_rhs,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_rhs,ierr)
         
         ! Prepare structured solution vector
         call HYPRE_StructVectorCreate    (this%cfg%comm,this%hypre_box,this%hypre_sol,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_sol,ierr)
         
      end select
      
      
      ! Initialize and setup the actual solver
      select case (this%method)
      case (bbox,pcg_bbox)
         
         ! The bbox solver was created already
         call this%bbox%update()
         
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
         call HYPRE_StructPFMGSetRelaxType (this%hypre_solver,2,ierr)
         call HYPRE_StructPFMGSetRAPType   (this%hypre_solver,1,ierr)
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
         call HYPRE_StructPFMGSetRelaxType(this%hypre_precond,2,ierr)
         call HYPRE_StructPFMGSetRAPType  (this%hypre_precond,1,ierr)
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
         
      case (gmres_pfmg)
         
         ! Create PFMG preconditioner
         call HYPRE_StructPFMGCreate      (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_StructPFMGSetMaxIter  (this%hypre_precond,1,ierr)
         call HYPRE_StructPFMGSetTol      (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_StructPFMGSetRelaxType(this%hypre_precond,2,ierr)
         call HYPRE_StructPFMGSetRAPType  (this%hypre_precond,1,ierr)
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
         
      end select
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine setup
   
   
   !> Solve the linear system iteratively
   subroutine solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count,ierr
      integer,  dimension(:), allocatable :: ind
      real(WP), dimension(:), allocatable :: rhs,sol
      
      ! Check that setup was done
      if (.not.this%setup_done) call die('[ils solve] Solver has not been setup.')
      
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      
      ! Transfer the RHS and IG to the appropriate environment
      select case (this%method)
      case (bbox,pcg_bbox)
         ! Transfer data to bbox
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%bbox%lvl(1)%f(i,j,k)=this%rhs(i,j,k)
                  this%bbox%lvl(1)%v(i,j,k)=this%sol(i,j,k)
               end do
            end do
         end do
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg,pcg)
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
      case (smg,pcg_smg,gmres_smg,pfmg,pcg_pfmg,gmres_pfmg)
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
      end select
      
      ! Call the appropriate solver
      select case (this%method)
      case (bbox,pcg_bbox)
         call this%bbox%solve()
         this%it  =this%bbox%my_ite
         this%rerr=this%bbox%my_res
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
      case (smg)
         call HYPRE_StructSMGSolve           (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructSMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructSMGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      case (pfmg)
         call HYPRE_StructPFMGSolve          (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructPFMGGetNumIteration(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructPFMGGetFinalRelativ(this%hypre_solver,this%rerr,ierr)
      case (pcg_smg,pcg_pfmg)
         call HYPRE_StructPCGSolve           (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructPCGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructPCGGetFinalRelative(this%hypre_solver,this%rerr,ierr)
      case (gmres_smg,gmres_pfmg)
         call HYPRE_StructGMRESSolve         (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
         call HYPRE_StructGMRESGetNumIteratio(this%hypre_solver,this%it  ,ierr)
         call HYPRE_StructGMRESGetFinalRelati(this%hypre_solver,this%rerr,ierr)
      end select
      
      ! Transfer the solution out of the appropriate environment and deallocate storage
      select case (this%method)
      case (bbox,pcg_bbox)
         ! Retreive solution from bbox
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%sol(i,j,k)=this%bbox%lvl(1)%v(i,j,k)
               end do
            end do
         end do
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg,pcg)
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
      case (smg,pcg_smg,gmres_smg,pfmg,pcg_pfmg,gmres_pfmg)
         ! Retrieve solution from structured vector
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  call HYPRE_StructVectorGetValues(this%hypre_sol,[i,j,k],this%sol(i,j,k),ierr)
               end do
            end do
         end do
      end select
      
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
      
   end subroutine solve
   
   
   !> Destroy solver - done everytime the operator changes
   subroutine destroy(this)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer :: ierr
      
      ! Destroy solver, operator, and rhs/sol vectors
      select case (this%method)
      case (bbox,pcg_bbox)
         ! Do nothing here since we don't have a bbox destructor yet
      case (amg)
         call HYPRE_BoomerAMGDestroy   (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (pcg_amg)
         call HYPRE_BoomerAMGDestroy   (this%hypre_precond,ierr)
         call HYPRE_ParCSRPCGDestroy   (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (pcg_parasail)
         call HYPRE_ParaSailsDestroy   (this%hypre_precond,ierr)
         call HYPRE_ParCSRPCGDestroy   (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (pcg)
         call HYPRE_ParCSRPCGDestroy   (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (gmres)
         call HYPRE_ParCSRGMRESDestroy (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (gmres_amg)
         call HYPRE_BoomerAMGDestroy   (this%hypre_precond,ierr)
         call HYPRE_ParCSRGMRESDestroy (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (gmres_pilut)
         call HYPRE_ParCSRPilutDestroy (this%hypre_precond,ierr)
         call HYPRE_ParCSRGMRESDestroy (this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (bicgstab_amg)
         call HYPRE_BoomerAMGDestroy     (this%hypre_precond,ierr)
         call HYPRE_ParCSRBiCGSTABDestroy(this%hypre_solver,ierr)
         call HYPRE_IJMatrixDestroy    (this%hypre_mat,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_rhs,ierr)
         call HYPRE_IJVectorDestroy    (this%hypre_sol,ierr)
      case (smg)
         call HYPRE_StructSMGDestroy   (this%hypre_solver,ierr)
         call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      case (pcg_smg)
         call HYPRE_StructSMGDestroy   (this%hypre_precond,ierr)
         call HYPRE_StructPCGDestroy   (this%hypre_solver,ierr)
         call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      case (gmres_smg)
         call HYPRE_StructSMGDestroy   (this%hypre_precond,ierr)
         call HYPRE_StructGMRESDestroy (this%hypre_solver,ierr)
         call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      case (pfmg)
         call HYPRE_StructPFMGDestroy  (this%hypre_solver,ierr)
         call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      case (pcg_pfmg)
         call HYPRE_StructPFMGDestroy  (this%hypre_precond,ierr)
         call HYPRE_StructPCGDestroy   (this%hypre_solver,ierr)
         call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      case (gmres_pfmg)
         call HYPRE_StructPFMGDestroy  (this%hypre_precond,ierr)
         call HYPRE_StructGMRESDestroy (this%hypre_solver,ierr)
         call HYPRE_StructMatrixDestroy(this%hypre_mat,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_rhs,ierr)
         call HYPRE_StructVectorDestroy(this%hypre_sol,ierr)
      end select
      
      ! Set setup-flag to false
      this%setup_done=.false.
      
   end subroutine destroy
   
   
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
   
   
   !> Short print of ILS info to the screen
   subroutine ils_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(ils), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("Iterative Linear Solver [",a16,"] for config [",a16,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   end subroutine ils_print_short
   
   
end module ils_class
