!> 3D FFT pressure solver concept is defined by extension of the linsol class
!> This solver is specifically intended to be a FFT-based pressure Poisson solver
!> for 3D periodic uniform computational domains decomposed in at most 2 directions
!> Uses rfourier_class for Fourier transform
module fft3d_class
   use precision,      only: WP
   use config_class,   only: config
   use rfourier_class, only: rfourier
   use linsol_class,   only: linsol
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: fft3d
   
   
   !> fft3d object definition
   type, extends(linsol) :: fft3d
      
      !> Fourier object
      class(rfourier), allocatable :: dft
      
      !> Unstrided arrays
      real(WP), dimension(:,:,:), allocatable :: transformed_opr
      real(WP), dimension(:,:,:), allocatable :: transformed_rhs
      
   contains
      
      procedure :: print_short=>fft3d_print_short !< One-line printing of solver status
      procedure :: print=>fft3d_print             !< Long-form printing of solver status
      procedure :: log=>fft3d_log                 !< Long-form logging of solver status
      procedure :: init=>fft3d_init               !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>fft3d_setup             !< Solver setup (every time the operator changes)
      procedure :: solve=>fft3d_solve             !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>fft3d_destroy         !< Solver destruction (every time the operator changes)
      
   end type fft3d
   
   
   !> Declare fft3d constructor
   interface fft3d
      procedure fft3d_from_args
   end interface fft3d
   
   
contains
   
   
   !> Constructor for a fft3d object
   function fft3d_from_args(cfg,name,nst) result(self)
      use messager, only: die
      implicit none
      type(fft3d) :: self
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
      allocate(self%dft,source=rfourier(self%cfg))
      
      ! Check FFT is available in all directions
      if (self%cfg%nx.gt.1.and..not.self%dft%xfft_avail) call die('[fft3d constructor] FFT is not available in x')
      if (self%cfg%ny.gt.1.and..not.self%dft%yfft_avail) call die('[fft3d constructor] FFT is not available in y')
      if (self%cfg%nz.gt.1.and..not.self%dft%zfft_avail) call die('[fft3d constructor] FFT is not available in z')
      
      ! Setup is not done
      self%setup_done=.false.
      
   end function fft3d_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine fft3d_init(this)
      use messager, only: die
      implicit none
      class(fft3d), intent(inout) :: this
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
         if (count([this%stc(st,1),this%stc(st,2),this%stc(st,3)].ne.0).gt.1) call die('[fft3d init] 3D stencil must be a sum of 1D stencils')
      end do
   end subroutine fft3d_init
   
   
   !> Setup solver - done everytime the operator changes
   subroutine fft3d_setup(this)
      use messager, only: die
      use mpi_f08,  only: MPI_BCAST,MPI_ALLREDUCE,MPI_INTEGER,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(fft3d), intent(inout) :: this
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      ! Check operator is circulant
      check_opr: block
         integer :: i,j,k,ierr
         logical :: circulant
         real(WP) :: circtol
         real(WP), dimension(1:this%nst) :: ref_opr
         if (this%cfg%amRoot) ref_opr=this%opr(:,this%cfg%imin,this%cfg%jmin,this%cfg%kmin)
         call MPI_BCAST(ref_opr,this%nst,MPI_REAL_WP,0,this%cfg%comm,ierr)
         circulant=.true.
         circtol=6.0_WP*epsilon(1.0_WP)/this%cfg%min_meshsize**4
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  if (any(abs(this%opr(:,i,j,k)-ref_opr).gt.circtol)) circulant=.false.
               end do
            end do
         end do
         if (.not.circulant) call die('[fft3d setup] operator must be uniform in xyz')
      end block check_opr
      
      ! Zero out transformed operator
      this%transformed_opr=0.0_WP
      
      ! Build the transformed operator in X
      trans_opr_x: block
         use mathtools, only: twoPi
         integer :: i,j,k,st
         real(WP) :: dk,kk
         complex(WP) :: k2eff
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
      end block trans_opr_x

      ! Build the transformed operator in Y
      trans_opr_y: block
         use mathtools, only: twoPi
         integer :: i,j,k,st
         real(WP) :: dk,kk
         complex(WP) :: k2eff
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
      end block trans_opr_y
      
      ! Build the transformed operator in Z
      trans_opr_z: block
         use mathtools, only: twoPi
         integer :: i,j,k,st
         real(WP) :: dk,kk
         complex(WP) :: k2eff
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
      end block trans_opr_z

      ! Make zero wavenumber equal to 1, therefore returning a solution with the same integral
      if (this%dft%oddball) this%transformed_opr(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_)=1.0_WP
      
      ! Make sure other wavenumbers are not close to zero
      check_trans_opr: block
         integer :: i,j,ierr
         i=count(abs(this%transformed_opr).lt.1000_WP*epsilon(1.0_WP))
         call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
         if (j.gt.0) call die('[fft3d setup] Elements of transformed operator are near zero')
      end block check_trans_opr
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine fft3d_setup
   
   
   !> Solve the linear system
   subroutine fft3d_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(fft3d), intent(inout) :: this
      ! Check that setup was done
      if (.not.this%setup_done) call die('[fft3d solve] Solver has not been setup.')
      ! Copy to unstrided array
      this%transformed_rhs=this%rhs(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)
      ! Forward transform in xyz
      call this%dft%xtransform_forward(this%transformed_rhs)
      call this%dft%ytransform_forward(this%transformed_rhs)
      call this%dft%ztransform_forward(this%transformed_rhs)
      ! Divide
      this%transformed_rhs=this%transformed_rhs/this%transformed_opr
      ! Backward transform in xyz
      call this%dft%xtransform_backward(this%transformed_rhs)
      call this%dft%ytransform_backward(this%transformed_rhs)
      call this%dft%ztransform_backward(this%transformed_rhs)
      ! Copy to strided output
      this%sol(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=this%transformed_rhs
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
   end subroutine fft3d_solve
   
   
   !> Log fft3d info
   subroutine fft3d_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(fft3d), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("fft3d solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
      end if
   end subroutine fft3d_log
   
   
   !> Print fft3d info to the screen
   subroutine fft3d_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(fft3d), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("fft3d solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
   end subroutine fft3d_print
   
   
   !> Short print of fft3d info to the screen
   subroutine fft3d_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(fft3d), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("fft3d solver [",a16,"] for config [",a16,"]")') trim(this%name),trim(this%cfg%name)
   end subroutine fft3d_print_short
   
   
   !> Solver destruction if the operator has changed
   subroutine fft3d_destroy(this)
     implicit none
     class(fft3d), intent(inout) :: this
     ! Set setup-flag to false
     this%setup_done=.false.
   end subroutine fft3d_destroy
   

end module fft3d_class