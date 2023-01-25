!> Black-box multigrid linear solver concept is
!> defined here by extension of the ilinsol class
module bbox_class
   use precision,     only: WP
   use config_class,  only: config
   use string,        only: str_medium
   use ilinsol_class, only: ilinsol
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: bbox
   
   
   ! List of known available methods
   integer, parameter, public ::     bbmg=1
   integer, parameter, public :: pcg_bbmg=2
   
   
   ! Level data object definition
   type :: lvl_data
      
      ! Global sizes
      integer :: ncell,ncello
      integer :: nx,imin,imax,nxo,imino,imaxo
      integer :: ny,jmin,jmax,nyo,jmino,jmaxo
      integer :: nz,kmin,kmax,nzo,kmino,kmaxo
      
      ! Local sizes
      integer :: ncell_,ncello_
      integer :: nx_,imin_,imax_,nxo_,imino_,imaxo_
      integer :: ny_,jmin_,jmax_,nyo_,jmino_,jmaxo_
      integer :: nz_,kmin_,kmax_,nzo_,kmino_,kmaxo_
      
      ! Operators
      real(WP), dimension(:,:,:,:,:,:), allocatable :: c2f
      real(WP), dimension(:,:,:,:,:,:), allocatable :: f2c
      real(WP), dimension(:,:,:,:,:,:), allocatable :: opr
      real(WP), dimension(:,:,:,:,:,:), allocatable :: oprc2f
      
      ! Vectors
      real(WP), dimension(:,:,:), allocatable :: v
      real(WP), dimension(:,:,:), allocatable :: f
      real(WP), dimension(:,:,:), allocatable :: r
      
      ! Parallel information
      logical, dimension(:), allocatable :: send_xm,send_xp
      logical, dimension(:), allocatable :: send_ym,send_yp
      logical, dimension(:), allocatable :: send_zm,send_zp
      integer :: nsend_xm,nsend_xp
      integer :: nsend_ym,nsend_yp
      integer :: nsend_zm,nsend_zp
      integer :: recv_xm,recv_xp
      integer :: recv_ym,recv_yp
      integer :: recv_zm,recv_zp
      
   end type lvl_data
   
   
   !> bbox object definition
   type, extends(ilinsol) :: bbox
      
      ! General operator information
      integer :: nstx                                                 !< Number of diagonal in x
      integer :: nsty                                                 !< Number of diagonal in y
      integer :: nstz                                                 !< Number of diagonal in z
      
      ! Direct solver data
      integer :: np                                                   !< Direct problem size
      integer , dimension(:),   allocatable :: piv                    !< Pivoting data
      real(WP), dimension(:),   allocatable :: myrhs,rhs              !< Local and assembled RHS
      real(WP), dimension(:,:), allocatable :: myOP,OP                !< Local and assembled operator
      
      ! Solver information
      real(WP) :: my_res0                                             !< Initial L_infty norm of the residual
      
      ! Configurable solver parameters
      integer  :: ncycle=1                                            !< Cycle type: 1=V-cycle (default), 2=W-cycle, ...
      integer  :: relax_pre =1                                        !< Number of pre-sweeps (default is 1)
      integer  :: relax_post=1                                        !< Number of post-sweeps (default is 1)
      integer  :: ncell_coarsest=0                                    !< Coarsest problem size allowed
      logical  :: use_direct_solve=.false.                            !< Use direct solve at the coarsest level (default is false)
      logical  :: use_krylov=.true.                                   !< Use BBMG as a preconditioner to a CG (default is true)
      
      ! Level data management
      integer :: nlvl                                                 !< Number of multigrid levels created - by default all will be used. This can be reduced by the user.
      type(lvl_data), dimension(:), allocatable :: lvl                !< Entire data at each level
      
      ! Krylov solver data
      real(WP), dimension(:,:,:), allocatable :: res,pp,zz,Ap,sol     !< Krylov solver work vectors
      real(WP) :: alpha,beta,rho1,rho2                                !< Krylov solver coefficients
      
   contains
      
      procedure :: print_short=>bbox_print_short         !< One-line printing of solver status
      procedure :: print=>bbox_print                     !< Long-form printing of solver status
      procedure :: log=>bbox_log                         !< Long-form logging of solver status
      procedure :: init=>bbox_init                       !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>bbox_setup                     !< Solver setup (every time the operator changes)
      procedure :: solve=>bbox_solve                     !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>bbox_destroy                 !< Solver destruction (every time the operator changes)
      
      procedure :: update                                             !< Operator update (every time the operator changes)
      procedure :: solve                                              !< Solve the linear system
      procedure :: initialize                                         !< Initialize the solver
      
      procedure, private :: mgsolve                                            !< Solve the linear system with multigrid
      procedure, private :: cgsolve                                            !< Solve the linear system with CG
      
      procedure, private :: recompute_prolongation                             !< Recompute the prolongation at a given level
      procedure, private :: recompute_restriction                              !< Recompute the restriction at a given level
      procedure, private :: recompute_operator                                 !< Recompute the operator at a given level
      procedure, private :: recompute_direct                                   !< Recompute the direct problem at the final level
      
      procedure, private :: get_residual                                       !< Compute the residual at a given level
      procedure, private :: c2f                                                !< Prolongate from a level to the next finer level
      procedure, private :: f2c                                                !< Restrict from a level to the next coarser level
      procedure, private :: relax                                              !< Relax the problem at a given level
      procedure, private :: direct_solve                                       !< Solve directly the problem at the final level
      procedure, private :: cycle                                              !< Perform a multigrid cycle
      
      generic, private :: sync=>vsync,msync                                    !< Synchronization at boundaries
      procedure, private :: vsync                                     !< Synchronize boundaries for a vector at a given level
      procedure, private :: msync                                     !< Synchronize boundaries for a matrix at a given level
      
      procedure, private :: pmodx,pmody,pmodz                                  !< Parity calculation that accounts for periodicity
      
   end type bbox
   
   
   !> Declare bbox constructor
   interface bbox
      procedure bbox_from_args
   end interface bbox
   
   
contains
   
   
   !> Division of an integer by 2
   pure integer function div(ind)
     implicit none
     integer, intent(in) :: ind
     div=ind/2+mod(ind,2)
   end function div
   
   
   !> Parity of point i shifted by n to ensure 1<=i<=n
   pure integer function pmodx(this,i,n)
      class(bbox), intent(in) :: this
      integer, intent(in) :: i,n
      if (this%cfg%xper) then
         pmodx=1-mod(mod(i+this%lvl(n)%nx-1,this%lvl(n)%nx)+1,2)
      else
         pmodx=1-mod(i,2)
      end if
   end function pmodx
   pure integer function pmody(this,i,n)
      class(bbox), intent(in) :: this
      integer, intent(in) :: i,n
      if (this%cfg%yper) then
         pmody=1-mod(mod(i+this%lvl(n)%ny-1,this%lvl(n)%ny)+1,2)
      else
         pmody=1-mod(i,2)
      end if
   end function pmody
   pure integer function pmodz(this,i,n)
      class(bbox), intent(in) :: this
      integer, intent(in) :: i,n
      if (this%cfg%zper) then
         pmodz=1-mod(mod(i+this%lvl(n)%nz-1,this%lvl(n)%nz)+1,2)
      else
         pmodz=1-mod(i,2)
      end if
   end function pmodz
   
   
   !> Constructor for a bbox object
   function bbox_from_args(cfg,name,method,nst) result(self)
      use messager, only: die
      implicit none
      type(bbox) :: self
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
      
      ! Setup is not done
      self%setup_done=.false.
      
   end function bbox_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine bbox_init(this)
      use messager, only: die
      implicit none
      class(bbox), intent(inout) :: this
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
      
      ! Initialize in-house black-box multi-grid solver
      bbox()
      
      ! Adjust parameters
      this%relax_pre =2
      this%relax_post=2
      this%use_direct_solve=.true.
      this%ncell_coarsest=300
      if (this%method.eq.    bbox) this%use_krylov=.false.
      if (this%method.eq.pcg_bbox) this%use_krylov=.true.
      
      ! Initialize solver
      call this%initialize()
      
   end subroutine bbox_init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine bbox_setup(this)
      implicit none
      class(bbox), intent(inout) :: this
      integer :: i,j,k,st
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
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
      
      ! Update bbox solver
      call this%bbox%update()
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine bbox_setup
   
   
   !> Solve the linear system iteratively
   subroutine bbox_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(bbox), intent(inout) :: this
      integer :: i,j,k
      
      ! Check that setup was done
      if (.not.this%setup_done) call die('[bbox solve] Solver has not been setup.')
      
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      
      ! Transfer the RHS and IG
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%bbox%lvl(1)%f(i,j,k)=this%rhs(i,j,k)
               this%bbox%lvl(1)%v(i,j,k)=this%sol(i,j,k)
            end do
         end do
      end do
      
      ! Call the solver
      call this%bbox%solve()
      this%it  =this%bbox%my_ite
      this%rerr=this%bbox%my_res
      
      ! Retrieve solution from bbox
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%sol(i,j,k)=this%bbox%lvl(1)%v(i,j,k)
            end do
         end do
      end do
      
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
      
   end subroutine bbox_solve
   
   
   !> Destroy solver - done everytime the operator changes
   subroutine bbox_destroy(this)
      implicit none
      class(bbox), intent(inout) :: this
      ! Set setup-flag to false
      this%setup_done=.false.
   end subroutine bbox_destroy
   
   
   !> Log bbox info
   subroutine bbox_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(bbox), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Bbox Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
         write(message,'(" >     method = ",i0)') this%method; call log(message)
         write(message,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit; call log(message)
         write(message,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg; call log(message)
         write(message,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg; call log(message)
      end if
   end subroutine bbox_log
   
   
   !> Print bbox info to the screen
   subroutine bbox_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(bbox), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("Bbox Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >     method = ",i0)') this%method
         write(output_unit,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit
         write(output_unit,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg
         write(output_unit,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg
      end if
   end subroutine bbox_print
   
   
   !> Short print of bbox info to the screen
   subroutine bbox_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(bbox), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("Bbox Linear Solver [",a16,"] for config [",a16,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   end subroutine ils_print_short
   
   
end module bbox_class
