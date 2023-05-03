!> Blackbox multigrid linear solver defined by extension of the linsol class
module bbox_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   use linsol_class, only: linsol
   use bbmg_class,   only: bbmg
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: bbox
   
   
   ! List of known available methods
   integer, parameter, public ::     mg=1
   integer, parameter, public :: pcg_mg=2
   
   
   !> Bbox object definition
   type, extends(linsol) :: bbox
      
      !> In-house bbmg solver object
      type(bbmg), allocatable :: solver
      
   contains
      
      procedure :: print_short=>bbox_print_short    !< One-line printing of solver status
      procedure :: print=>bbox_print                !< Long-form printing of solver status
      procedure :: log=>bbox_log                    !< Long-form logging of solver status
      procedure :: init=>bbox_init                  !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>bbox_setup                !< Solver setup (every time the operator changes)
      procedure :: solve=>bbox_solve                !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>bbox_destroy            !< Solver destruction (every time the operator changes)
      
   end type bbox
   
   
   !> Declare bbox constructor
   interface bbox
      procedure bbox_from_args
   end interface bbox
   
   
contains
   
   
   !> Constructor for an bbox object
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
      if (nst.gt.27) call die('[bbox constructor] bbox cannot work on operators larger than 3x3x3.')
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
      integer :: st,stx1,stx2,sty1,sty2,stz1,stz2
      
      ! From the provided stencil, generate an inverse map
      stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1))
      sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2))
      stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3))
      allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
      do st=1,this%nst
         this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
      end do
      
      ! Check stencil again
      if (2*maxval(abs(this%stc))+1.gt.3) call die('[bbox init] bbox cannot work on operators larger than 3x3x3.')
      
      ! Initialize in-house black-box multigrid solver
      this%solver=bbmg(pg=this%cfg,name=trim(this%name),nst=[3,3,3])
      
      ! Adjust parameters
      this%solver%relax_pre =3
      this%solver%relax_post=3
      this%solver%use_direct_solve=.false.
      this%solver%ncell_coarsest=10
      select case (this%method)
      case (mg)
         this%solver%use_krylov=.false.
      case (pcg_mg)
         this%solver%use_krylov=.true.
      case default
         call die('[bbox init] Unknown bbox solution method.')
      end select

      ! Initialize solver
      call this%solver%initialize()
      
   end subroutine bbox_init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine bbox_setup(this)
      use messager, only: die
      implicit none
      class(bbox), intent(inout) :: this
      integer :: i,j,k,st
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      ! Storage for matrix/vectors is already available within bbox, so just copy operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%solver%lvl(1)%opr(:,:,:,i,j,k)=0.0_WP
               do st=1,this%nst
                  this%solver%lvl(1)%opr(this%stc(st,1),this%stc(st,2),this%stc(st,3),i,j,k)=this%solver%lvl(1)%opr(this%stc(st,1),this%stc(st,2),this%stc(st,3),i,j,k)+this%opr(st,i,j,k)
               end do
            end do
         end do
      end do
      
      ! Update convergence criteria
      this%solver%max_ite=this%maxit
      this%solver%max_res=this%rcvg
      
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
      
      ! Transfer data to bbox
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%solver%lvl(1)%f(i,j,k)=this%rhs(i,j,k)
               this%solver%lvl(1)%v(i,j,k)=this%sol(i,j,k)
            end do
         end do
      end do
      
      ! Call the solver and get back convergence info
      call this%solver%solve()
      this%it  =this%solver%my_ite
      this%rerr=this%solver%my_res
      
      ! Retreive solution from bbox
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%sol(i,j,k)=this%solver%lvl(1)%v(i,j,k)
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
      use messager, only: die
      implicit none
      class(bbox), intent(inout) :: this
      
      ! Do nothing here since we don't have a bbox destructor yet
      
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
         write(message,'("bbox solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
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
         write(output_unit,'("bbox solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
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
      if (this%cfg%amRoot) write(output_unit,'("bbox solver [",a16,"] for config [",a16,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   end subroutine bbox_print_short
   
   
end module bbox_class
