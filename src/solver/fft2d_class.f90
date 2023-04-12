!> 2D FFT-tridiagonal pressure solver concept is defined by extension of the linsol class
!> This solver is specifically intended to be a FFT-tridiagonal pressure Poisson solver
!> for computational domains that have 2 periodic uniform directions and are decomposed
!> in at most 2 directions
!> Uses rfourier_class for Fourier transform
module fft2d_class
   use precision,      only: WP
   use config_class,   only: config
   use rfourier_class, only: rfourier
   use diag_class,     only: diag
   use linsol_class,   only: linsol
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: fft2d
   
   
   !> fft2d object definition
   type, extends(linsol) :: fft2d
      
      !> Fourier object
      class(rfourier), allocatable :: dft

      !> Diagonal solver
      class(diag), allocatable :: diagsolver
      logical :: xdiag,ydiag,zdiag
      
      !> Unstrided arrays
      real(WP), dimension(:,:,:), allocatable :: transformed_opr
      real(WP), dimension(:,:,:), allocatable :: transformed_rhs
      
   contains
      
      procedure :: print_short=>fft2d_print_short !< One-line printing of solver status
      procedure :: print=>fft2d_print             !< Long-form printing of solver status
      procedure :: log=>fft2d_log                 !< Long-form logging of solver status
      procedure :: init=>fft2d_init               !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>fft2d_setup             !< Solver setup (every time the operator changes)
      procedure :: solve=>fft2d_solve             !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>fft2d_destroy         !< Solver destruction (every time the operator changes)
      
   end type fft2d
   
   
   !> Declare fft2d constructor
   interface fft2d
      procedure fft2d_from_args
   end interface fft2d
   
   
contains
   
   
   !> Constructor for a fft2d object
   function fft2d_from_args(cfg,name,nst) result(self)
      use messager, only: die
      implicit none
      type(fft2d) :: self
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
      
      ! Allocate unstrided arrays
      allocate(self%transformed_opr(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
      allocate(self%transformed_rhs(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
      
      ! Create a rfourier object
      self%dft=rfourier(self%cfg)

      ! Check FFT is only unavailable in one direction
      self%xdiag=.false.; if (self%cfg%nx.gt.1.and..not.self%dft%xfft_avail) self%xdiag=.true.
      self%ydiag=.false.; if (self%cfg%ny.gt.1.and..not.self%dft%yfft_avail) self%ydiag=.true.
      self%zdiag=.false.; if (self%cfg%nz.gt.1.and..not.self%dft%zfft_avail) self%zdiag=.true.
      if (count([self%xdiag,self%ydiag,self%zdiag]).ne.1) call die('[fft2d constructor] Only one direction should rely on diagonal solver')
      
      ! Setup is not done
      self%setup_done=.false.
      
   end function fft2d_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine fft2d_init(this)
      use messager, only: die
      implicit none
      class(fft2d), intent(inout) :: this
      integer :: st,stx1,stx2,sty1,sty2,stz1,stz2
      ! From the provided stencil, generate an inverse map
      stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1))
      sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2))
      stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3))
      allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
      do st=1,this%nst
         this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
      end do
      ! Also check that the stencil is a sum of 1D-stencils, an assumption that is made when building the operator below
      do st=1,this%nst
         if (count([this%stc(st,1),this%stc(st,2),this%stc(st,3)].ne.0).gt.1) call die('[fft2d init] 3D stencil must be a sum of 1D stencils')
      end do
      ! Create a diagonal solver now that we know the stencil
      this%diagsolver=diag(this%cfg,name=this%name,n=2*maxval(abs(this%stc))+1)
   end subroutine fft2d_init
   
   
   !> Setup solver - done every time the operator changes
   subroutine fft2d_setup(this)
      use mathtools, only: twoPi
      implicit none
      class(fft2d), intent(inout) :: this
      integer :: i,j,k,st
      real(WP) :: dk,kk
      complex(WP) :: k2eff
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      ! NOTE THAT WE ARE NOT CHECKING THE VALIDITY OF THE OPERATOR HERE...
      
      ! Zero out transformed operator
      this%transformed_opr=0.0_WP
      
      ! Build the transformed operator in X
      if (this%dft%xfft_avail) then
         dk=twoPi/this%cfg%xL
         do i=this%cfg%imin_,this%cfg%imax_
            ! Get the wavenumber
            kk=real(i-this%cfg%imin,WP)*dk
            if (i-this%cfg%imin.gt.(this%cfg%nx/2)) kk=real(this%cfg%nx-i+this%cfg%imin,WP)*dk
            ! Compute effective wavenumber of X operator build by summation
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  k2eff=(0.0_WP,0.0_WP)
                  do st=1,this%nst
                     k2eff=k2eff+this%opr(st,i,j,k)*exp((0.0_WP,1.0_WP)*kk*this%cfg%xm(i+this%stc(st,1)))
                  end do
                  k2eff=k2eff/exp((0.0_WP,1.0_WP)*kk*this%cfg%xm(i))
                  this%transformed_opr(i,j,k)=this%transformed_opr(i,j,k)+realpart(k2eff)
               end do
            end do
         end do
      end if
      
      ! Build the transformed operator in Y
      if (this%dft%yfft_avail) then
         dk=twoPi/this%cfg%yL
         do j=this%cfg%jmin_,this%cfg%jmax_
            ! Get the wavenumber
            kk=real(j-this%cfg%jmin,WP)*dk
            if (j-this%cfg%jmin.gt.(this%cfg%ny/2)) kk=real(this%cfg%ny-j+this%cfg%jmin,WP)*dk
            ! Compute effective wavenumber of Y operator build by summation
            do k=this%cfg%kmin_,this%cfg%kmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  k2eff=(0.0_WP,0.0_WP)
                  do st=1,this%nst
                     k2eff=k2eff+this%opr(st,i,j,k)*exp((0.0_WP,1.0_WP)*kk*this%cfg%ym(j+this%stc(st,2)))
                  end do
                  k2eff=k2eff/exp((0.0_WP,1.0_WP)*kk*this%cfg%ym(j))
                  this%transformed_opr(i,j,k)=this%transformed_opr(i,j,k)+realpart(k2eff)
               end do
            end do
         end do
      end if
      
      ! Build the transformed operator in Z
      if (this%dft%zfft_avail) then
         dk=twoPi/this%cfg%zL
         do k=this%cfg%kmin_,this%cfg%kmax_
            ! Get the wavenumber
            kk=real(k-this%cfg%kmin,WP)*dk
            if (k-this%cfg%kmin.gt.(this%cfg%nz/2)) kk=real(this%cfg%nz-k+this%cfg%kmin,WP)*dk
            ! Compute effective wavenumber of Z operator build by summation
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  k2eff=(0.0_WP,0.0_WP)
                  do st=1,this%nst
                     k2eff=k2eff+this%opr(st,i,j,k)*exp((0.0_WP,1.0_WP)*kk*this%cfg%zm(k+this%stc(st,3)))
                  end do
                  k2eff=k2eff/exp((0.0_WP,1.0_WP)*kk*this%cfg%zm(k))
                  this%transformed_opr(i,j,k)=this%transformed_opr(i,j,k)+realpart(k2eff)
               end do
            end do
         end do
      end if
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine fft2d_setup
   
   
   !> Solve the linear system
   subroutine fft2d_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(fft2d), intent(inout) :: this
      integer :: i,j,k,st

      ! Check that setup was done
      if (.not.this%setup_done) call die('[fft2d solve] Solver has not been setup.')

      ! Copy to unstrided array
      this%transformed_rhs=this%rhs(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)
      
      ! Solve diagonal system in X
      if (this%xdiag) then
         call this%dft%ytransform_forward(this%transformed_rhs)
         call this%dft%ztransform_forward(this%transformed_rhs)
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%diagsolver%Ax(j,k,i,:)=0.0_WP
                  this%diagsolver%Ax(j,k,i,0)=this%transformed_opr(i,j,k)
                  do st=1,this%nst
                     this%diagsolver%Ax(j,k,i,this%stc(st,1))=this%diagsolver%Ax(j,k,i,this%stc(st,1))+this%opr(st,i,j,k)
                  end do
                  this%diagsolver%Rx(j,k,i)=this%transformed_rhs(i,j,k)
               end do
            end do
         end do
         call this%diagsolver%linsol_x()
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%transformed_rhs(i,j,k)=this%diagsolver%Rx(j,k,i)
               end do
            end do
         end do
         call this%dft%ytransform_backward(this%transformed_rhs)
         call this%dft%ztransform_backward(this%transformed_rhs)
      end if

      ! Solve diagonal system in Y
      if (this%ydiag) then
         call this%dft%xtransform_forward(this%transformed_rhs)
         call this%dft%ztransform_forward(this%transformed_rhs)
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%diagsolver%Ay(i,k,j,:)=0.0_WP
                  this%diagsolver%Ay(i,k,j,0)=this%transformed_opr(i,j,k)
                  do st=1,this%nst
                     this%diagsolver%Ay(i,k,j,this%stc(st,2))=this%diagsolver%Ay(i,k,j,this%stc(st,2))+this%opr(st,i,j,k)
                  end do
                  this%diagsolver%Ry(i,k,j)=this%transformed_rhs(i,j,k)
               end do
            end do
         end do
         call this%diagsolver%linsol_y()
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%transformed_rhs(i,j,k)=this%diagsolver%Ry(i,k,j)
               end do
            end do
         end do
         call this%dft%xtransform_backward(this%transformed_rhs)
         call this%dft%ztransform_backward(this%transformed_rhs)
      end if

      ! Solve diagonal system in Z
      if (this%zdiag) then
         call this%dft%xtransform_forward(this%transformed_rhs)
         call this%dft%ytransform_forward(this%transformed_rhs)
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%diagsolver%Az(i,j,k,:)=0.0_WP
                  this%diagsolver%Az(i,j,k,0)=this%transformed_opr(i,j,k)
                  do st=1,this%nst
                     this%diagsolver%Az(i,j,k,this%stc(st,3))=this%diagsolver%Az(i,j,k,this%stc(st,3))+this%opr(st,i,j,k)
                  end do
                  this%diagsolver%Rz(i,j,k)=this%transformed_rhs(i,j,k)
               end do
            end do
         end do
         call this%diagsolver%linsol_z()
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%transformed_rhs(i,j,k)=this%diagsolver%Rz(i,j,k)
               end do
            end do
         end do
         call this%dft%xtransform_backward(this%transformed_rhs)
         call this%dft%ytransform_backward(this%transformed_rhs)
      end if
      
      ! Copy to strided output
      this%sol(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=this%transformed_rhs

      ! Sync the solution vector
      call this%cfg%sync(this%sol)

      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short

   end subroutine fft2d_solve
   
   
   !> Log fft2d info
   subroutine fft2d_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(fft2d), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("fft2d solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
      end if
   end subroutine fft2d_log
   
   
   !> Print fft2d info to the screen
   subroutine fft2d_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(fft2d), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("fft2d solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
   end subroutine fft2d_print
   
   
   !> Short print of fft2d info to the screen
   subroutine fft2d_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(fft2d), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("fft2d solver [",a16,"] for config [",a16,"]")') trim(this%name),trim(this%cfg%name)
   end subroutine fft2d_print_short
   
   
   !> Solver destruction if the operator has changed
   subroutine fft2d_destroy(this)
     implicit none
     class(fft2d), intent(inout) :: this
     ! Set setup-flag to false
     this%setup_done=.false.
   end subroutine fft2d_destroy
   

end module fft2d_class