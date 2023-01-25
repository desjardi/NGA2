!> 3D FFT pressure solver concept is defined by extension of the linsol class
!> This solver is specifically intended to be a FFT-based pressure Poisson solver
!> for 3D periodic uniform computational domains decomposed in at most 2 directions
!> Makes use of the P3DFFT++
module pfft3d_class
  use p3dfft_plus_plus
  use precision,    only: WP
  use config_class, only: config
  use linsol_class, only: linsol
  implicit none
  private

  ! Expose type/constructor/methods
  public :: pfft3d

  !> pfft3d object definition
  type, extends(linsol) :: pfft3d

    logical :: oddball
    integer :: trans_f, trans_b
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: factored_operator, transformed_rhs
    real(C_DOUBLE), dimension(:,:,:), allocatable :: unstrided_rhs, unstrided_sol

  contains

    procedure :: log => pfft3d_log                  !< Long-form logging of solver status
    procedure :: print => pfft3d_print              !< Long-form printing of solver status
    procedure :: print_short => pfft3d_print_short  !< One-line printing of solver status
    procedure :: init => pfft3d_init                !< Grid and stencil initialization - done once for the grid and stencil
    procedure :: setup => pfft3d_setup              !< Solver setup (every time the operator changes)
    procedure :: solve => pfft3d_solve              !< Execute solver (assumes new RHS and initial guess at every call)
    procedure :: destroy => pfft3d_destroy          !< Solver destruction (every time the operator changes)

  end type pfft3d

  !> Declare pfft3d constructor
  interface pfft3d; procedure pfft3d_from_args; end interface pfft3d;

contains

  !> Log pfft3d info
  subroutine pfft3d_log(this)
    use string,   only: str_long
    use messager, only: log
    implicit none
    class(pfft3d), intent(in) :: this
    character(len=str_long) :: message

    if (this%cfg%amRoot) then
      write(message,'("PFFT3D solver [",a,"] for config [",a,"]")')           &
        trim(this%name), trim(this%cfg%name)
      call log(message)
    end if

  end subroutine pfft3d_log

  !> Print pfft3d info to the screen
  subroutine pfft3d_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none

    class(pfft3d), intent(in) :: this
    if (this%cfg%amRoot) then
      write(output_unit,'("PFFT3D solver [",a,"] for config [",a,"]")')       &
        trim(this%name), trim(this%cfg%name)
    end if

  end subroutine pfft3d_print

  !> Short print of pfft3d info to the screen
  subroutine pfft3d_print_short(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(pfft3d), intent(in) :: this

    if (this%cfg%amRoot) write(output_unit,'("PFFT3D solver [",a16,"] for     &
      config [",a16,"]")') trim(this%name), trim(this%cfg%name)

  end subroutine pfft3d_print_short

  !> Constructor for an pfft3d object
  function pfft3d_from_args(cfg, name, nst) result(self)
    use messager, only: die
    implicit none
    type(pfft3d) :: self
    class(config), target, intent(in) :: cfg
    character(len=*), intent(in) :: name
    integer, intent(in) :: nst
    integer :: imn_, imx_, jmn_, jmx_, kmn_, kmx_, imno_, imxo_, jmno_, jmxo_,&
      kmno_, kmxo_

    ! Link the config and store the name
    self%cfg => cfg
    self%name = trim(adjustl(name))

    imn_ = self%cfg%imin_; imx_ = self%cfg%imax_;
    jmn_ = self%cfg%jmin_; jmx_ = self%cfg%jmax_;
    kmn_ = self%cfg%kmin_; kmx_ = self%cfg%kmax_;
    imno_ = self%cfg%imino_; imxo_ = self%cfg%imaxo_;
    jmno_ = self%cfg%jmino_; jmxo_ = self%cfg%jmaxo_;
    kmno_ = self%cfg%kmino_; kmxo_ = self%cfg%kmaxo_;

    ! Set up stencil size and map
    self%nst = nst; allocate(self%stc(self%nst,3));

    ! Allocate rhs, operator, transformed operator diagonal, and solution
    allocate(self%factored_operator(imn_:imx_,jmn_:jmx_,kmn_:kmx_))
    allocate(self%transformed_rhs(imn_:imx_,jmn_:jmx_,kmn_:kmx_))
    allocate(self%unstrided_rhs(imn_:imx_,jmn_:jmx_,kmn_:kmx_))
    allocate(self%unstrided_sol(imn_:imx_,jmn_:jmx_,kmn_:kmx_))
    allocate(self%opr(self%nst,imno_:imxo_,jmno_:jmxo_,kmno_:kmxo_))
    allocate(self%rhs(imno_:imxo_,jmno_:jmxo_,kmno_:kmxo_))
    allocate(self%sol(imno_:imxo_,jmno_:jmxo_,kmno_:kmxo_))

    ! Zero out some info
    self%method = 0; self%it = 0; self%aerr = 0.0_WP; self%rerr = 0.0_WP;

    ! Setup is not done
    self%setup_done = .false.

  end function pfft3d_from_args

  !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
  !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
  !> Only the stencil needs to be defined at this point
  subroutine pfft3d_init(this)
    !use mpi_f08,  only: MPI_COMM_WORLD
    use messager, only: die
    implicit none
    include "mpif.h"
    class(pfft3d), intent(inout) :: this
    integer :: type_ccr, type_rcc, trans_f_tmp, trans_b_tmp
    integer(C_INT) :: p3dfft_grid, grid_real, grid_fourier
    integer, dimension(3) :: type_ids_f, type_ids_b,                          &
      glob_start_real, glob_start_fourier
    integer(C_INT), dimension(3) :: pdims, ldims_real, ldims_fourier,         &
      mem_order_real, mem_order_fourier, dmap_real, dmap_fourier, gdims_real, &
      gdims_fourier

    ! Set up global dimensions of the grid
    gdims_real(:) = (/ this%cfg%nx, this%cfg%ny, this%cfg%nz /)

    ! Define the initial processor grid
    pdims(:) = (/ this%cfg%npx, this%cfg%npy, this%cfg%npz /)

    ! Various checks to ensure we can use this solver
    check_solver_is_useable: block
      integer :: ndcp, ndims

      ! Periodicity and uniformity of mesh
      if (.not.(this%cfg%xper.and.this%cfg%uniform_x)) call die('[pfft3d constructor] Need x-direction needs to be periodic and uniform')
      if (.not.(this%cfg%yper.and.this%cfg%uniform_y)) call die('[pfft3d constructor] Need y-direction needs to be periodic and uniform')
      if (.not.(this%cfg%zper.and.this%cfg%uniform_z)) call die('[pfft3d constructor] Need z-direction needs to be periodic and uniform')

      ! Ensure that we have at least one non-decomposed direction
      ndims = count(gdims_real .gt. 1); ndcp = count(pdims .gt. 1);
      if (ndcp .ge. ndims) call die('[pdfft3 constructor] Need at least one NON-decomposed direction')

    end block check_solver_is_useable

    ! Set up work structures for P3DFFT
    call p3dfft_setup

    ! Set up 2 transform types for 3D transforms
    type_ids_f(:) = (/ P3DFFT_R2CFFT_D, P3DFFT_CFFT_FORWARD_D, P3DFFT_CFFT_FORWARD_D /)
    type_ids_b(:) = (/ P3DFFT_C2RFFT_D, P3DFFT_CFFT_BACKWARD_D, P3DFFT_CFFT_BACKWARD_D /)

    ! Now initialize 3D transforms (forward and backward) with these types
    call p3dfft_init_3Dtype(type_rcc, type_ids_f)
    call p3dfft_init_3Dtype(type_ccr, type_ids_b)

    ! Set up processor order and memory ordering, as well as the final global
    ! grid dimensions.  These will be different from the original dimensions in
    ! one dimension due to conjugate symmetry, since we are doing
    ! real-to-complex transform.
    gdims_fourier(:) = gdims_real(:)
    gdims_fourier(1) = gdims_fourier(1) / 2 + 1
    mem_order_real(:) = (/ 0, 1, 2 /)

    ! Set up memory order for the final grid layout (for complex array in
    ! Fourier space). It is more convenient to have the storage order of
    ! the array reversed, this helps save on memory access bandwidth, and
    ! shouldn't affect the operations in the Fourier space very much,
    ! requiring basically a change in the loop order. However it is possible
    ! to define the memory ordering the same as default (0,1,2). Note that
    ! the memory ordering is specified in C indices, i.e. starting from 0.
    mem_order_fourier(:) = (/ 1, 2, 0 /)
    dmap_real(:) = (/ 0, 1, 2 /); dmap_fourier(:) = (/ 1, 2, 0 /);

    ! Specify the default communicator for P3DFFT++. This can be different
    ! from your program default communicator.
    p3dfft_grid = p3dfft_init_proc_grid(pdims, MPI_COMM_WORLD)

    ! Initialize initial grid, no conjugate symmetry (-1)
    call p3dfft_init_data_grid(grid_real, ldims_real, glob_start_real,        &
      gdims_real, -1, p3dfft_grid, dmap_real, mem_order_real)

    ! Check ldims_real against pgrid
    if (.not. all((/ this%cfg%nx_, this%cfg%ny_, this%cfg%nz_ /) .eq.         &
      ldims_real)) then
      call die("[pfftd3] parallel fft decomposition does not match cfg")
    end if

    ! Final grid has conjugate symmetry in X dimension (0)
    call p3dfft_init_data_grid(grid_fourier, ldims_fourier, glob_start_fourier,&
      gdims_fourier, 0, p3dfft_grid, dmap_fourier, mem_order_fourier)

    ! Plan transforms
    trans_f_tmp = this%trans_f
    write(*,*) this%trans_f
    write(*,*) grid_real
    write(*,*) grid_fourier
    write(*,*) type_rcc
    call p3dfft_plan_3Dtrans(trans_f_tmp, grid_real, grid_fourier, type_rcc)
    write(*,*) "============"
    write(*,*) this%trans_f
    write(*,*) grid_real
    write(*,*) grid_fourier
    write(*,*) type_rcc
    this%trans_f = trans_f_tmp
    call p3dfft_plan_3Dtrans(this%trans_b, grid_fourier, grid_real, type_ccr)

    ! tag process with zero wavenumber
    this%oddball = all((/ this%cfg%iproc, this%cfg%jproc, this%cfg%kproc /) .eq. 1)

  end subroutine pfft3d_init

  !> Setup solver - done everytime the operator changes
  subroutine pfft3d_setup(this)
    use messager, only: die
    implicit none
    integer :: i, j, k, n, stx1, stx2, sty1, sty2, stz1, stz2
    class(pfft3d), intent(inout) :: this
    logical :: circulent
    real(WP), dimension(this%nst) :: ref_stencil
    real(C_DOUBLE), dimension(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:   &
      this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_) :: opr_col

    ! compute stencil inverse
    stx1 = minval(this%stc(:,1)); stx2 = maxval(this%stc(:,1));
    sty1 = minval(this%stc(:,2)); sty2 = maxval(this%stc(:,2));
    stz1 = minval(this%stc(:,3)); stz2 = maxval(this%stc(:,3));
    allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2))
    do n = 1, this%nst
      this%stmap(this%stc(n,1),this%stc(n,2),this%stc(n,3)) = n
    end do

    ! check circulent operator
    ref_stencil = this%opr(:,1,1,1)
    circulent = .true.
    do k = this%cfg%kmin_, this%cfg%kmax_
      do j = this%cfg%jmin_, this%cfg%jmax_
        do i = this%cfg%imin_, this%cfg%imax_
          if (any(this%opr(:,i,j,k) .ne. ref_stencil(:))) circulent = .false.
        end do
      end do
    end do
    if (.not. circulent) then
      call die("[pf3dft] stencil must be uniform in space")
    end if

    ! build
    opr_col(:,:,:) = 0.0_C_DOUBLE
    do n = 1, this%nst
      i = modulo(this%stc(n,1) - this%cfg%imin_ + 1, this%cfg%nx_) + this%cfg%imin_
      j = modulo(this%stc(n,2) - this%cfg%jmin_ + 1, this%cfg%ny_) + this%cfg%jmin_
      k = modulo(this%stc(n,3) - this%cfg%kmin_ + 1, this%cfg%nz_) + this%cfg%kmin_
      opr_col(i,j,k) = opr_col(i,j,k) + real(ref_stencil(n), C_DOUBLE)
    end do
    
    ! take transform of operator
    call p3dfft_3Dtrans_double(this%trans_f, opr_col, this%factored_operator, 0)
    this%factored_operator = 1.0_C_DOUBLE / this%factored_operator

    ! make zero wavenumber not zero
    if (this%oddball) this%factored_operator(this%cfg%imin_,this%cfg%jmin_,   &
      this%cfg%kmin_) = 1.0_WP

    ! Set setup-flag to true
    this%setup_done = .true.

  end subroutine pfft3d_setup

  !> do the solve
  subroutine pfft3d_solve(this)
    use param,    only: verbose
    implicit none
    class(pfft3d), intent(inout) :: this
    integer :: imn, imx, jmn, jmx, kmn, kmx

    imn = this%cfg%imin_; imx = this%cfg%imax_;
    jmn = this%cfg%jmin_; jmx = this%cfg%jmax_;
    kmn = this%cfg%kmin_; kmx = this%cfg%kmax_;

    ! copy to unstrided array
    this%unstrided_rhs(:,:,:) = this%rhs(imn:imx,jmn:jmx,kmn:kmx)

    ! do forward transform
    call p3dfft_3Dtrans_double(this%trans_f, this%unstrided_rhs, this%transformed_rhs, 0)

    ! divide
    this%transformed_rhs = this%transformed_rhs * this%factored_operator

    ! fix zero wavenumber
    if (this%oddball) this%transformed_rhs(imn,jmn,kmn) = (0.0_WP, 0.0_WP)

    ! do backward transform
    call p3dfft_3Dtrans_double(this%trans_b, this%transformed_rhs, this%unstrided_sol, 0)

    ! copy to strided output
    this%sol(:,:,:) = this%unstrided_sol(imn:imx,jmn:jmx,kmn:kmx)

    ! sync
    call this%cfg%sync(this%sol)

    ! If verbose run, log and or print info
    if (verbose.gt.0) call this%log
    if (verbose.gt.1) call this%print_short

  end subroutine pfft3d_solve

  subroutine pfft3d_destroy(this)
    implicit none
    class(pfft3d), intent(inout) :: this

    ! do nothing; included for consistency with other methods
    this%setup_done = .false.

  end subroutine pfft3d_destroy

end module pfft3d_class

