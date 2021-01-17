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
   integer, parameter, public :: gmres_amg=9
   integer, parameter, public :: bicgstab_amg=10
   
   
   ! List of key solver parameters
   integer , parameter :: gmres_kdim=5                 !< Number of basis vectors between restarts
   real(WP), parameter :: amg_strong_threshold=0.25_WP !< Coarsening parameter (default is 0.25, 0.5 recommended in 3D)
   integer , parameter :: amg_coarsen_type=8           !< Falgout=6 (old default); PMIS=8 and HMIS=10 (recommended)
   integer , parameter :: amg_interp=6                 !< 6=extended classical modified interpolation (default); 8=standard interpolation
   integer , parameter :: amg_printlvl=0               !< 0=none (default); 3=init and cvg history
   integer , parameter :: amg_relax=6                  !< 6=Hybrid symmetric Gauss-Seidel (default); 8=symmetric L1-Gauss-Seidel; 0=Weighted Jacobi
   integer , parameter :: amg_relax_coarse=6           !< 9=Gauss elim; 99=GE w/ pivot; may be good to use same as above?
   
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
      procedure :: print_short=>ils_print_short                       !< Short print of ILS object info
      procedure :: print=>ils_print                                   !< Print ILS object info
      procedure :: log  =>ils_log                                     !< Log ILS object info
      procedure :: init_solver                                        !< Solver initialization (at start-up)
      procedure :: update_solver                                      !< Solver update (every time the operator changes)
      procedure :: solve                                              !< Equation solver
      procedure :: solve_rbgs                                         !< Solve equation with rbgs
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
   
   
   !> Initialize solver - done at start-up only as long as the stencil or cfg does not change
   !> When creating the solver, zero values of this%opr(1,:,:,:)
   !> indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine init_solver(this,method)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer, intent(in) :: method
      integer :: ierr,st,stx1,stx2,sty1,sty2,stz1,stz2
      integer, dimension(3) :: periodicity,offset
      integer, dimension(:), allocatable :: sizes
      
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
      
      ! Initialize matrix and vectors - this also catches incorrect solvers
      select case (this%method)
      case (rbgs)
         ! Nothing needed
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg)  ! Initialize a HYPRE IJ enviroment
         ! Allocate and prepare unstructed mapping
         call this%prep_umap()
         ! Create a HYPRE matrix
         call HYPRE_IJMatrixCreate       (this%cfg%comm,this%ind_min,this%ind_max,this%ind_min,this%ind_max,this%hypre_mat,ierr)
         call HYPRE_IJMatrixSetObjectType(this%hypre_mat,hypre_ParCSR,ierr)
         allocate(sizes(this%ind_min:this%ind_max)); sizes=this%nst
         call HYPRE_IJMatrixSetRowSizes  (this%hypre_mat,sizes,ierr)
         deallocate(sizes)
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
      case (smg,pfmg) ! Initialize a structured HYPRE environment - we will need to pay attention to 1D vs 2D vs 3D
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
      case default
         call die('[ils prep solver] Unknown solution method')
      end select
      
      ! Now initialize the actual solver
      select case (this%method)
      case (rbgs)
         ! Nothing needed
      case (amg)
         ! Create an AMG solver
         call HYPRE_BoomerAMGCreate        (this%hypre_solver,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_solver,amg_printlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_solver,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_solver,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_solver,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_solver,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_solver,amg_relax_coarse,3,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_solver,this%rcvg ,ierr)
      case (pcg_amg)
         ! Create a PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,0,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         ! Create an AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,amg_printlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_precond,amg_relax_coarse,3,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         ! Set AMG as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,2,this%hypre_precond,ierr)
      case (pcg_parasail)
         ! Create a PCG solver
         call HYPRE_ParCSRPCGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (this%hypre_solver,0,ierr)
         call HYPRE_ParCSRPCGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRPCGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (this%hypre_solver,1,ierr)
         call HYPRE_ParCSRPCGSetLogging    (this%hypre_solver,1,ierr)
         ! Create a PARASAIL preconditioner
         call HYPRE_ParaSailsCreate        (this%cfg%comm,this%hypre_precond,ierr)
         call HYPRE_ParaSailsSetParams     (this%hypre_precond,0.1_WP,1,ierr)
         call HYPRE_ParaSailsSetFilter     (this%hypre_precond,0.05_WP,ierr)
         call HYPRE_ParaSailsSetSym        (this%hypre_precond,0,ierr)
         call HYPRE_ParaSailsSetLogging    (this%hypre_precond,1,ierr)
         ! Set PARASAIL as preconditioner to PCG
         call HYPRE_ParCSRPCGSetPrecond    (this%hypre_solver,4,this%hypre_precond,ierr)
      case (gmres)
         ! Create a GMRES solver
         call HYPRE_ParCSRGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetKDim     (this%hypre_solver,gmres_kdim,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging  (this%hypre_solver,1,ierr)
      case (gmres_amg)
         ! Create a GMRES solver
         call HYPRE_ParCSRGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetKDim     (this%hypre_solver,gmres_kdim,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging  (this%hypre_solver,1,ierr)
         ! Create an AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,amg_printlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_precond,amg_relax_coarse,3,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         ! Set AMG as preconditioner to GMRES
         call HYPRE_ParCSRGMRESSetPrecond  (this%hypre_solver,2,this%hypre_precond,ierr)
      case (gmres_pilut)
         ! Create a GMRES solver
         call HYPRE_ParCSRGMRESCreate      (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRGMRESSetKDim     (this%hypre_solver,gmres_kdim,ierr)
         call HYPRE_ParCSRGMRESSetMaxIter  (this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRGMRESSetTol      (this%hypre_solver,this%rcvg,ierr)
         call HYPRE_ParCSRGMRESSetLogging  (this%hypre_solver,1,ierr)
         ! Create a PILUT preconditioner
         call HYPRE_ParCSRPilutCreate      (this%cfg%comm,this%hypre_precond,ierr)
         ! Set PILUT as preconditioner to GMRES
         call HYPRE_ParCSRGMRESSetPrecond  (this%hypre_solver,3,this%hypre_precond,ierr)
      case (bicgstab_amg)
         ! Create a BiCGstab solver
         call HYPRE_ParCSRBiCGSTABCreate    (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_ParCSRBiCGSTABSetPrintLev(this%hypre_solver,0,ierr)
         call HYPRE_ParCSRBiCGSTABSetMaxIter(this%hypre_solver,this%maxit,ierr)
         call HYPRE_ParCSRBiCGSTABSetTol    (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_ParCSRBiCGSTABSetLogging(this%hypre_solver,1,ierr)
         ! Create an AMG preconditioner
         call HYPRE_BoomerAMGCreate        (this%hypre_precond,ierr)
         call HYPRE_BoomerAMGSetPrintLevel (this%hypre_precond,amg_printlvl,ierr)
         call HYPRE_BoomerAMGSetInterpType (this%hypre_precond,amg_interp,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(this%hypre_precond,amg_coarsen_type,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(this%hypre_precond,amg_strong_threshold,ierr)
         call HYPRE_BoomerAMGSetRelaxType  (this%hypre_precond,amg_relax,ierr)
         call HYPRE_BoomerAMGSetCycleRelaxType(this%hypre_precond,amg_relax_coarse,3,ierr)
         call HYPRE_BoomerAMGSetMaxIter    (this%hypre_precond,1,ierr)
         call HYPRE_BoomerAMGSetTol        (this%hypre_precond,0.0_WP,ierr)
         call HYPRE_BoomerAMGSetNumSweeps  (this%hypre_precond,1,ierr)
         ! Set AMG as preconditioner to BiCGstab
         call HYPRE_ParCSRBiCGSTABSetPrecond(this%hypre_solver,2,this%hypre_precond,ierr)
      case (smg)
         ! Create a SMG solver
         call HYPRE_StructSMGCreate        (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructSMGSetMaxIter    (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructSMGSetTol        (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_StructSMGSetLogging    (this%hypre_solver,1,ierr)
      case (pfmg)
         ! Create a PFMG solver
         call HYPRE_StructPFMGCreate       (this%cfg%comm,this%hypre_solver,ierr)
         call HYPRE_StructPFMGSetRelaxType (this%hypre_solver,2,ierr)
         call HYPRE_StructPFMGSetRAPType   (this%hypre_solver,1,ierr)
         call HYPRE_StructPFMGSetMaxIter   (this%hypre_solver,this%maxit,ierr)
         call HYPRE_StructPFMGSetTol       (this%hypre_solver,this%rcvg ,ierr)
         call HYPRE_StructPFMGSetLogging   (this%hypre_solver,1,ierr)
         call HYPRE_StructPFMGSetPrintLevel(this%hypre_solver,0,ierr)
      end select
      
   end subroutine init_solver
   
   
   !> Update solver - done everytime the operator changes
   subroutine update_solver(this)
      use messager, only: die
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count1,count2,st,ierr,nn
      integer,  dimension(:), allocatable :: row,ncol,mycol,col
      real(WP), dimension(:), allocatable :: myval,val
      
      ! First update the matrix
      select case (this%method)
      case (rbgs)
         ! Nothing to do
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg)  ! Prepare the IJ operator
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
                  if (this%ind(i,j,k).gt.0) then
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
                           if (mycol(st).eq.mycol(nn).and.mycol(st).gt.0) then
                              mycol(st)=0
                              myval(nn)=myval(nn)+myval(st)
                           end if
                        end do
                     end do
                     ! Finally transfer to operator
                     do st=1,this%nst
                        if (mycol(st).gt.0) then
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
      case (smg,pfmg) ! Prepare the struct operator
         ! Prepare local storage
         allocate(row(1:this%nst))
         do st=1,this%nst
            row(st)=st-1
         end do
         allocate(val(1:this%nst))
         ! Transfer the operator
         call HYPRE_StructMatrixInitialize(this%hypre_mat,ierr)
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
      end select
      
      ! Setup the solver
      select case (this%method)
      case (rbgs)
         ! Nothing to do
      case (amg)
         call HYPRE_BoomerAMGSetup  (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      case (pcg_amg,pcg_parasail)
         call HYPRE_ParCSRPCGSetup  (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      case (gmres,gmres_pilut,gmres_amg)
         call HYPRE_ParCSRGMRESSetup(this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      case (bicgstab_amg)
         call HYPRE_ParCSRBiCGSTABSetup(this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
      case (smg)
         call HYPRE_StructSMGSetup  (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      case (pfmg)
         call HYPRE_StructPFMGSetup (this%hypre_solver,this%hypre_mat,this%hypre_rhs,this%hypre_sol,ierr)
      end select
      
   end subroutine update_solver
   
   
   !> Solve the linear system iteratively
   subroutine solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(ils), intent(inout) :: this
      integer :: i,j,k,count,ierr
      integer,  dimension(:), allocatable :: ind
      real(WP), dimension(:), allocatable :: rhs,sol
      
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      
      ! Transfer the RHS and IG to the appropriate environment
      select case (this%method)
      case (rbgs)
         ! Nothing to do
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg)  ! Prepare the IJ vectors
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
         call HYPRE_IJVectorInitialize(this%hypre_rhs,ierr)
         call HYPRE_IJVectorSetValues (this%hypre_rhs,this%ncell_,ind,rhs,ierr)
         call HYPRE_IJVectorAssemble  (this%hypre_rhs,ierr)
         call HYPRE_IJVectorInitialize(this%hypre_sol,ierr)
         call HYPRE_IJVectorSetValues (this%hypre_sol,this%ncell_,ind,sol,ierr)
         call HYPRE_IJVectorAssemble  (this%hypre_sol,ierr)
      case (smg,pfmg) ! Prepare the struct vectors
         call HYPRE_StructVectorInitialize(this%hypre_rhs,ierr)
         call HYPRE_StructVectorInitialize(this%hypre_sol,ierr)
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
         call HYPRE_StructVectorAssemble(this%hypre_rhs,ierr)
         call HYPRE_StructVectorAssemble(this%hypre_sol,ierr)
      end select
      
      ! Call the appropriate solver
      select case (this%method)
      case (rbgs)
         call this%solve_rbgs()
      case (amg)
         call HYPRE_BoomerAMGSolve           (this%hypre_solver,this%parse_mat,this%parse_rhs,this%parse_sol,ierr)
         call HYPRE_BoomerAMGGetNumIterations(this%hypre_solver,this%it  ,ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(this%hypre_solver,this%rerr,ierr)
      case (pcg_amg,pcg_parasail)
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
      end select
      
      ! Transfer the SOL out of the appropriate environment and deallocate storage
      select case (this%method)
      case (rbgs)
         ! Nothing to do
      case (amg,pcg_amg,pcg_parasail,gmres,gmres_pilut,gmres_amg,bicgstab_amg)  ! Retrieve the IJ vector
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
         deallocate(rhs,sol,ind)
      case (smg,pfmg) ! Retrieve the struct vector
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
   
   
   !> Short print of ILS info to the screen
   subroutine ils_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(ils), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("Iterative Linear Solver [",a20,"] for config [",a20,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   end subroutine ils_print_short
   
   
end module ils_class
