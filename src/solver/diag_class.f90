!> Direct diagonal solvers defined here
module diag_class
  use precision, only: WP
  use config_class,   only: config
  use string,       only: str_medium
  implicit none
  private

  ! Expose type/constructor/methods
  public :: diag

  !> diag object definition
  type :: diag
     ! Diagonalr solver works for a config
     type(config), pointer :: cfg                                    !< Config for the diag solver
     character(len=str_medium) :: name                               !< Name of solver
     integer :: ndiags                                               !< Number of diagonals
     real(WP), dimension(:,:,:,:), allocatable :: Ax,Ay,Az           !< Working arrays
     real(WP), dimension(:,:,:), allocatable :: Rx,Ry,Rz             !< Solution arrays
     real(WP), dimension(:,:), allocatable :: stackmem               !< Work arrays

   contains
     procedure :: linsol_x                                           !< Linear solver in x
     procedure :: linsol_y                                           !< Linear solver in y
     procedure :: linsol_z                                           !< Linear solver in z
     procedure :: tridiagonal
     procedure :: pentadiagonal
     procedure :: polydiagonal
     final     :: destructor                                         !< Destructor for diag
  end type diag

  !> Declare diag constructor
  interface diag
     procedure diag_from_args
  end interface diag

contains

  !> Destructor for diag object
  subroutine destructor(this)
    implicit none
    type(diag) :: this
    if (allocated(this%Ax)) deallocate(this%Ax)
    if (allocated(this%Ay)) deallocate(this%Ay)
    if (allocated(this%Az)) deallocate(this%Az)
    if (allocated(this%Rx)) deallocate(this%Rx)
    if (allocated(this%Ry)) deallocate(this%Ry)
    if (allocated(this%Rz)) deallocate(this%Rz)
    if (allocated(this%stackmem)) deallocate(this%stackmem)
  end subroutine destructor


  !> Constructor for a diag object
  function diag_from_args(cfg,name,n) result(self)
    use messager, only: die
    implicit none
    type(diag) :: self
    class(config), target, intent(in) :: cfg
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    integer :: nst
    
    ! Link the config and store the name
    self%cfg=>cfg
    self%name=trim(adjustl(name))

    ! Number of diagonals and stencil size
    self%ndiags=n
    nst=(n-1)/2

    ! Allocate arrays
    allocate(self%Ax(self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_,self%cfg%imin_:self%cfg%imax_,-nst:+nst))
    allocate(self%Ay(self%cfg%imin_:self%cfg%imax_,self%cfg%kmin_:self%cfg%kmax_,self%cfg%jmin_:self%cfg%jmax_,-nst:+nst))
    allocate(self%Az(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_,-nst:+nst))
    allocate(self%Rx(self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_,self%cfg%imin_:self%cfg%imax_))
    allocate(self%Ry(self%cfg%imin_:self%cfg%imax_,self%cfg%kmin_:self%cfg%kmax_,self%cfg%jmin_:self%cfg%jmax_))
    allocate(self%Rz(self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_))
    allocate(self%stackmem((self%cfg%imaxo_-self%cfg%imino_)*(self%cfg%jmaxo_-self%cfg%jmino_)*(self%cfg%kmaxo_-self%cfg%kmino_),nst+1))

  end function diag_from_args


  ! Real solver in x
  subroutine linsol_x(this)
    implicit none
    class(diag), intent(inout) :: this

    ! Choose based on number of diagonals
    select case(this%ndiags)
    case(3)
       call this%tridiagonal(&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,-1),&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_, 0),&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,+1),&
            this%Rx,this%cfg%nx_,this%cfg%ny_*this%cfg%nz_,'x',this%stackmem(1,1),this%stackmem(1,2))
    case(5)
       call this%pentadiagonal(&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,-2),&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,-1),&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_, 0),&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,+1),&
            this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,+2),&
            this%Rx,this%cfg%nx_,this%cfg%ny_*this%cfg%nz_,'x',this%stackmem(1,1),this%stackmem(1,3))
    case default
       call this%polydiagonal((this%ndiags-1)/2,this%Ax(this%cfg%jmin_,this%cfg%kmin_,this%cfg%imin_,&
            -(this%ndiags-1)/2),this%Rx,this%cfg%nx_,this%cfg%ny_*this%cfg%nz_,'x',this%stackmem)
    end select

    return
  end subroutine linsol_x


  ! Real solver in y
  subroutine linsol_y(this)
    implicit none
    class(diag), intent(inout) :: this

    ! Choose based on number of diagonals
    select case(this%ndiags)
    case(3)
       call this%tridiagonal(&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,-1),&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_, 0),&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,+1),&
            this%Ry,this%cfg%ny_,this%cfg%nx_*this%cfg%nz_,'y',this%stackmem(1,1),this%stackmem(1,2))
    case(5)
       call this%pentadiagonal(&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,-2),&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,-1),&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_, 0),&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,+1),&
            this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,+2),&
            this%Ry,this%cfg%ny_,this%cfg%nx_*this%cfg%nz_,'y',this%stackmem(1,1),this%stackmem(1,3))
    case default
       call this%polydiagonal((this%ndiags-1)/2,this%Ay(this%cfg%imin_,this%cfg%kmin_,this%cfg%jmin_,&
            -(this%ndiags-1)/2),this%Ry,this%cfg%ny_,this%cfg%nx_*this%cfg%nz_,'y',this%stackmem)
    end select

    return
  end subroutine linsol_y


  ! Real solver in z
  subroutine linsol_z(this)
    implicit none
    class(diag), intent(inout) :: this

    ! Choose based on number of diagonals
    select case(this%ndiags)
    case(3)
       call this%tridiagonal(&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,-1),&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_, 0),&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,+1),&
            this%Rz,this%cfg%nz_,this%cfg%ny_*this%cfg%nx_,'z',this%stackmem(1,1),this%stackmem(1,2))
    case(5)
       call this%pentadiagonal(&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,-2),&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,-1),&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_, 0),&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,+1),&
            this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,+2),&
            this%Rz,this%cfg%nz_,this%cfg%ny_*this%cfg%nx_,'z',this%stackmem(1,1),this%stackmem(1,3))
    case default
       call this%polydiagonal((this%ndiags-1)/2,this%Az(this%cfg%imin_,this%cfg%jmin_,this%cfg%kmin_,&
            -(this%ndiags-1)/2),this%Rz,this%cfg%nz_,this%cfg%ny_*this%cfg%nx_,'z',this%stackmem)
    end select

    return
  end subroutine linsol_z


  subroutine tridiagonal(this,A,B,C,R,n,lot,dir,s1,s2)
    use parallel
    use mpi_f08, only : mpi_comm
    implicit none
    class(diag), intent(inout) :: this
    ! Direction
    character(len=*), intent(in) :: dir
    ! Size of the problems
    integer, intent(in) :: n
    ! Number of problems
    integer, intent(in) :: lot
    ! Matrix
    real(WP), dimension(lot,n) :: A     ! LOWER
    real(WP), dimension(lot,n) :: B     ! DIAGONAL
    real(WP), dimension(lot,n) :: C     ! UPPER
    real(WP), dimension(lot,n) :: R     ! RHS - RESULT
    ! Local
    real(WP), dimension(lot)   :: const
    real(WP), dimension(lot)   :: r1
    real(WP), dimension(lot)   :: r2
    real(WP), dimension(lot,n) :: s1
    real(WP), dimension(lot,n) :: s2
    ! Communication
    logical :: nper
    type(MPI_Comm) :: ncom
    integer :: proc,nrank
    integer :: nremain,nlot
    real(WP), dimension(:),     allocatable :: sendbuf
    real(WP), dimension(:,:,:), allocatable :: recvbuf1
    real(WP), dimension(:,:),   allocatable :: recvbuf2
    integer,  dimension(:),     allocatable :: ngroup
    real(WP), dimension(:,:),   allocatable :: buf1,buf2,buf3,buf4,buf5,buf6
    ! Stuff
    integer :: i,igroup
    integer :: k1,L,k2,nk
    integer :: ierr

    ! Get parallel info
    select case (trim(adjustl(dir)))
    case ('x')
       proc = this%cfg%npx
       nrank = this%cfg%xrank
       ncom = this%cfg%xcomm
       nper = this%cfg%xper
    case ('y')
       proc = this%cfg%npy
       nrank = this%cfg%yrank
       ncom = this%cfg%ycomm
       nper = this%cfg%yper
    case ('z')
       proc = this%cfg%npz
       nrank = this%cfg%zrank
       ncom = this%cfg%zcomm
       nper = this%cfg%zper
    case default
       stop 'Unknown direction'
    end select

    ! If serial
    if (proc .eq. 1) then
       if (.not.nper) then
          call tridiagonal_serial(A,B,C,R,n,lot)
       else
          call tridiagonal_periodic_serial(A,B,C,R,n,lot)
       end if
       return
    end if

    ! Partition the lot
    if (lot .lt. proc) stop 'Tridiagonal solver cannot handle so many proc for such a small problem.'
    allocate(ngroup(proc))
    ngroup(:) = lot/proc
    nremain = mod(lot,proc)
    ngroup(1:nremain) = ngroup(1:nremain) + 1
    nlot = ngroup(1)
    allocate(sendbuf(nlot*12*proc))
    allocate(recvbuf1(nlot,6,2*proc))
    allocate(recvbuf2(nlot,2*proc))

    ! Initialize boundary values
    s1(:,1) = a(:,1)
    s2(:,n) = c(:,n)
    if (.not.nper) then
       if (nrank .eq. 0)      s1(:,1) = 0.0_WP
       if (nrank .eq. proc-1) s2(:,n) = 0.0_WP
    end if

    ! Forward elimination
    ! Upper boundary in s1(i)
    do i=2,n
       const(:) = a(:,i)/b(:,i-1)
       b(:,i)   = b(:,i) - c(:,i-1)*const(:)
       r(:,i)   = r(:,i) - r(:,i-1)*const(:)
       s1(:,i)  = -s1(:,i-1)*const(:)
    end do

    ! Backward elimination
    ! Lower boundary in s2(i)
    do i=n-1,1,-1
       const(:) = c(:,i)/b(:,i+1)
       r(:,i)   = r(:,i) - r(:,i+1)*const(:)
       s1(:,i)  = s1(:,i) - s1(:,i+1)*const(:)
       s2(:,i)  = -s2(:,i+1)*const(:)
    end do

    ! All dependence has been shifted to the boundary elements
    ! Communicate boundary values to root process
    ! and solve reduced pentadiagonal system
    ! Use of pentadiagonal system is more robust than the
    ! reordered (removes zeros) tridiagonal system

    ! Send rows of pentadiagonal system
    ! (0, s1, b, 0, s2; r)
    !    (s1, 0, b, s2, 0; r)

    L = 1
    k1 = 1
    do igroup=1,proc
       k2 = k1+ngroup(igroup)-1
       nk = ngroup(igroup)

       sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
       sendbuf(L:L+nk-1) = s1(k1:k2,1) ; L = L + nlot
       sendbuf(L:L+nk-1) = b(k1:k2,1)  ; L = L + nlot
       sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
       sendbuf(L:L+nk-1) = s2(k1:k2,1) ; L = L + nlot
       sendbuf(L:L+nk-1) = r(k1:k2,1)  ; L = L + nlot
       sendbuf(L:L+nk-1) = s1(k1:k2,n) ; L = L + nlot
       sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
       sendbuf(L:L+nk-1) = b(k1:k2,n)  ; L = L + nlot
       sendbuf(L:L+nk-1) = s2(k1:k2,n) ; L = L + nlot
       sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
       sendbuf(L:L+nk-1) = r(k1:k2,n)  ; L = L + nlot

       k1 = k2 + 1
    end do

    ! Gather the boundary data
    call MPI_ALLTOALL(sendbuf,nlot*12,MPI_REAL_WP,recvbuf1,nlot*12,MPI_REAL_WP,ncom,ierr)

    ! Clear unused values
    nk = ngroup(nrank+1)
    recvbuf1(nk+1:nlot, :, :) = 0.0_WP
    recvbuf1(nk+1:nlot, 3, :) = 1.0_WP

    ! Solve reduced systems
    if (.not.nper) then
       ! Store in buf to avoid creating a temporary array (there must be a better way to do this)
       allocate(buf1(nlot,2*proc),buf2(nlot,2*proc),buf3(nlot,2*proc),buf4(nlot,2*proc),buf5(nlot,2*proc),buf6(nlot,2*proc))
       buf1=recvbuf1(:,1,2:2*proc-1)
       buf2=recvbuf1(:,2,2:2*proc-1)
       buf3=recvbuf1(:,3,2:2*proc-1)
       buf4=recvbuf1(:,4,2:2*proc-1)
       buf5=recvbuf1(:,5,2:2*proc-1)
       buf6=recvbuf1(:,6,2:2*proc-1)
       call pentadiagonal_serial(buf1,buf2,buf3,buf4,buf5,buf6,2*proc-2,nlot)
       recvbuf1(:,1,2:2*proc-1)=buf1
       recvbuf1(:,2,2:2*proc-1)=buf2
       recvbuf1(:,3,2:2*proc-1)=buf3
       recvbuf1(:,4,2:2*proc-1)=buf4
       recvbuf1(:,5,2:2*proc-1)=buf5
       recvbuf1(:,6,2:2*proc-1)=buf6
       deallocate(buf1,buf2,buf3,buf4,buf5,buf6)
    else
       allocate(buf1(nlot,2*proc),buf2(nlot,2*proc),buf3(nlot,2*proc),buf4(nlot,2*proc),buf5(nlot,2*proc),buf6(nlot,2*proc))
       buf1=recvbuf1(:,1,:)
       buf2=recvbuf1(:,2,:)
       buf3=recvbuf1(:,3,:)
       buf4=recvbuf1(:,4,:)
       buf5=recvbuf1(:,5,:)
       buf6=recvbuf1(:,6,:)
       call pentadiagonal_periodic_serial(buf1,buf2,buf3,buf4,buf5,buf6,&
            2*proc,nlot,sendbuf(1),sendbuf(2*proc*nlot+1))
       recvbuf1(:,1,:)=buf1
       recvbuf1(:,2,:)=buf2
       recvbuf1(:,3,:)=buf3
       recvbuf1(:,4,:)=buf4
       recvbuf1(:,5,:)=buf5
       recvbuf1(:,6,:)=buf6
       deallocate(buf1,buf2,buf3,buf4,buf5,buf6)
    end if

    ! Move solution to first slot
    do i=1,2*proc
       recvbuf2(:,i) = recvbuf1(:,6,i)
    end do

    ! Permute the order
    do i=1,proc-1
       const(1:nlot) = recvbuf2(:,2*i)
       recvbuf2(:,2*i) = recvbuf2(:,2*i+1)
       recvbuf2(:,2*i+1) = const(1:nlot)
    end do

    ! If periodic, don't forget the end points
    if (nper) then
       const(1:nlot) = recvbuf2(:,1)
       recvbuf2(:,1) = recvbuf2(:,2*proc)
       recvbuf2(:,2*proc) = const(1:nlot)
    end if

    ! Scatter back the solution
    call MPI_ALLTOALL(recvbuf2,nlot*2,MPI_REAL_WP,sendbuf,nlot*2,MPI_REAL_WP,ncom,ierr)

    L = 1
    k1 = 1
    do igroup=1,proc
       k2 = k1+ngroup(igroup)-1
       nk = k2-k1+1

       r1(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot
       r2(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot

       k1 = k2 + 1
    end do

    ! Only if not periodic
    if (.not.nper) then
       if (nrank .eq. 0)      r1 = 0.0_WP
       if (nrank .eq. proc-1) r2 = 0.0_WP
    end if

    do i=1,n
       r(:,i) = (r(:,i) - s1(:,i)*r1(:) - s2(:,i)*r2(:))/b(:,i)
    end do

    deallocate(sendbuf)
    deallocate(recvbuf2)
    deallocate(recvbuf1)
    deallocate(ngroup)

    return
  end subroutine tridiagonal


  !> TriDiagonal Solver - Serial Case - Periodic
  subroutine tridiagonal_periodic_serial(a,b,c,r,n,lot)
    implicit none

    integer, intent(in) :: n,lot
    real(WP), intent(inout), dimension(lot,n) :: a,b,c,r
    real(WP), dimension(lot) :: const
    integer :: i

    if (n .eq. 1) then
       r = r/(a + b + c)
       return
    else if (n .eq. 2) then
       ! Solve 2x2 system
       c(:,1) = c(:,1) + a(:,1)
       a(:,2) = a(:,2) + c(:,2)
       const(:) = a(:,2)/b(:,1)
       b(:,2) = b(:,2) - c(:,1)*const(:)
       r(:,2) = r(:,2) - r(:,1)*const(:)
       r(:,2) = r(:,2)/b(:,2)
       r(:,1) = (r(:,1) - c(:,1)*r(:,2))/b(:,1)
       return
    end if

    ! Forward elimination
    do i=2,n-2
       const(:) = a(:,i)/b(:,i-1)
       b(:,i) = b(:,i) - c(:,i-1)*const(:)
       r(:,i) = r(:,i) - r(:,i-1)*const(:)
       ! Boundary is stored in a(i)
       a(:,i) = -a(:,i-1)*const(:)
    end do
    i=n-1
    const(:) = a(:,i)/b(:,i-1)
    b(:,i) = b(:,i) - c(:,i-1)*const(:)
    r(:,i) = r(:,i) - r(:,i-1)*const(:)
    a(:,i) = c(:,i) - a(:,i-1)*const(:)
    i=n
    const(:) = a(:,i)/b(:,i-1)
    r(:,i) = r(:,i) - r(:,i-1)*const(:)
    a(:,i) = b(:,i) - a(:,i-1)*const(:)

    ! Backward elimination
    do i=n-2,1,-1
       const(:) = c(:,i)/b(:,i+1)
       r(:,i) = r(:,i) - r(:,i+1)*const(:)
       a(:,i) = a(:,i) - a(:,i+1)*const(:)
    end do

    ! Eliminate oddball
    const(:) = c(:,n)/b(:,1)
    r(:,n) = r(:,n) - r(:,1)*const(:)
    a(:,n) = a(:,n) - a(:,1)*const(:)

    ! Backward substitution
    r(:,n) = r(:,n)/a(:,n)
    do i=n-1,1,-1
       r(:,i) = (r(:,i) - a(:,i)*r(:,n))/b(:,i)
    end do

    return
  end subroutine tridiagonal_periodic_serial


  !> TriDiagonal Solver - Serial Case - Not Periodic
  subroutine tridiagonal_serial(a,b,c,r,n,lot)
    implicit none

    integer, intent(in) :: n,lot
    real(WP), intent(inout), dimension(lot,n) :: a,b,c,r
    real(WP), dimension(lot) :: const
    integer :: i

    ! Forward elimination
    do i=2,n
       const(:) = a(:,i)/(b(:,i-1)+tiny(1.0_WP))
       b(:,i) = b(:,i) - c(:,i-1)*const(:)
       r(:,i) = r(:,i) - r(:,i-1)*const(:)
    end do

    ! Back-substitution
    r(:,n) = r(:,n)/b(:,n)
    do i=n-1,1,-1
       r(:,i) = (r(:,i) - c(:,i)*r(:,i+1))/(b(:,i)+tiny(1.0_WP))
    end do

    return
  end subroutine tridiagonal_serial

  subroutine pentadiagonal(this,A,B,C,D,E,R,n,lot,dir,s1,s2)
    use parallel
    use mpi_f08, only : mpi_comm
    implicit none
    class(diag), intent(inout) :: this
    ! Direction
    character(len=*), intent(in) :: dir
    ! Size of the problems
    integer, intent(in) :: n
    ! Number of problems
    integer, intent(in) :: lot
    ! Matrix
    real(WP), dimension(lot,n) :: A     ! LOWER-2
    real(WP), dimension(lot,n) :: B     ! LOWER-1
    real(WP), dimension(lot,n) :: C     ! DIAGONAL
    real(WP), dimension(lot,n) :: D     ! UPPER+1
    real(WP), dimension(lot,n) :: E     ! UPPER+2
    real(WP), dimension(lot,n) :: R     ! RHS - RESULT
    ! Local
    real(WP), dimension(lot)     :: const
    real(WP), dimension(lot,2)   :: r1
    real(WP), dimension(lot,2)   :: r2
    real(WP), dimension(lot,n,2) :: s1
    real(WP), dimension(lot,n,2) :: s2
    ! Communication
    logical :: nper
    type(MPI_Comm) :: ncom
    integer :: proc,nrank
    integer :: nremain,nlot
    real(WP), dimension(:,:),   allocatable :: sendbuf
    real(WP), dimension(:,:),   allocatable :: recvbuf
    integer,  dimension(:),     allocatable :: ngroup
    real(WP), dimension(:,:,:), allocatable :: AA
    real(WP), dimension(:,:),   allocatable :: swap

    ! Stuff
    integer :: i,igroup,j,k,m
    integer :: k1,L,k2,nk
    integer :: ierr

    ! Get parallel info
    select case (trim(adjustl(dir)))
    case ('x')
       proc = this%cfg%npx
       nrank = this%cfg%xrank
       ncom = this%cfg%xcomm
       nper = this%cfg%xper
    case ('y')
       proc = this%cfg%npy
       nrank = this%cfg%yrank
       ncom = this%cfg%ycomm
       nper = this%cfg%yper
    case ('z')
       proc = this%cfg%npz
       nrank = this%cfg%zrank
       ncom = this%cfg%zcomm
       nper = this%cfg%zper
    case default
       stop 'Unknown direction'
    end select

    ! If serial
    if (proc .eq. 1) then
       if (.not.nper) then
          call pentadiagonal_serial(A,B,C,D,E,R,n,lot)
       else
          call pentadiagonal_periodic_serial(A,B,C,D,E,R,n,lot,s1,s2)
       end if
       return
    end if

    ! Partition the lot
    if (lot .lt. proc) stop 'Pentadiagonal solver cannot handle so many proc for such a small problem.'
    allocate(ngroup(proc))
    ngroup(:) = lot/proc
    nremain = mod(lot,proc)
    ngroup(1:nremain) = ngroup(1:nremain) + 1
    nlot = ngroup(1)
    allocate(sendbuf(nlot,24*proc))
    allocate(recvbuf(nlot,24*proc))

    if (n > 1) then ! do normal stuff

       ! Initialize boundary values
       s1(:,1,1) = a(:,1)
       s1(:,1,2) = b(:,1)
       s1(:,2,1) = 0.0_WP
       s1(:,2,2) = a(:,2)
       s2(:,n-1,1) = e(:,n-1)
       s2(:,n-1,2) = 0.0_WP
       s2(:,n,1) = d(:,n)
       s2(:,n,2) = e(:,n)
       if (.not.nper) then
          if (nrank .eq. 0)      s1(:,1:2,1:2)   = 0.0_WP
          if (nrank .eq. proc-1) s2(:,n-1:n,1:2) = 0.0_WP
       end if

       ! Forward elimination
       ! Upper boundary in s1(:,i,1:2)
       do i=1,n-2
          ! Eliminate a(i+2)
          const(:) = a(:,i+2)/c(:,i)
          b(:,i+2) = b(:,i+2) - d(:,i)*const(:)
          c(:,i+2) = c(:,i+2) - e(:,i)*const(:)
          r(:,i+2) = r(:,i+2) - r(:,i)*const(:)
          s1(:,i+2,1) = -s1(:,i,1)*const(:)
          s1(:,i+2,2) = -s1(:,i,2)*const(:)

          ! Eliminate b(i+1)
          const(:) = b(:,i+1)/c(:,i)
          c(:,i+1) = c(:,i+1) - d(:,i)*const(:)
          d(:,i+1) = d(:,i+1) - e(:,i)*const(:)
          r(:,i+1) = r(:,i+1) - r(:,i)*const(:)
          s1(:,i+1,1) = s1(:,i+1,1) - s1(:,i,1)*const(:)
          s1(:,i+1,2) = s1(:,i+1,2) - s1(:,i,2)*const(:)
       end do
       ! Eliminate b(n)
       const(:) = b(:,n)/c(:,n-1)
       c(:,n) = c(:,n) - d(:,n-1)*const(:)
       r(:,n) = r(:,n) - r(:,n-1)*const(:)
       s1(:,n,1) = s1(:,n,1) - s1(:,n-1,1)*const(:)
       s1(:,n,2) = s1(:,n,2) - s1(:,n-1,2)*const(:)
       s2(:,n,1) = s2(:,n,1) - s2(:,n-1,1)*const(:)

       ! Backward elimination
       ! Lower boundary in s2(:,i,1:2)
       do i=n,3,-1
          ! Eliminate e(i-2)
          const(:) = e(:,i-2)/c(:,i)
          r(:,i-2) = r(:,i-2) - r(:,i)*const(:)
          s1(:,i-2,1) = s1(:,i-2,1) - s1(:,i,1)*const(:)
          s1(:,i-2,2) = s1(:,i-2,2) - s1(:,i,2)*const(:)
          s2(:,i-2,1) = -s2(:,i,1)*const(:)
          s2(:,i-2,2) = -s2(:,i,2)*const(:)

          ! Eliminate d(i-1)
          const(:) = d(:,i-1)/c(:,i)
          r(:,i-1) = r(:,i-1) - r(:,i)*const(:)
          s1(:,i-1,1) = s1(:,i-1,1) - s1(:,i,1)*const(:)
          s1(:,i-1,2) = s1(:,i-1,2) - s1(:,i,2)*const(:)
          s2(:,i-1,1) = s2(:,i-1,1) - s2(:,i,1)*const(:)
          s2(:,i-1,2) = s2(:,i-1,2) - s2(:,i,2)*const(:)
       end do
       ! Eliminate d(1)
       const(:) = d(:,1)/c(:,2)
       r(:,1) = r(:,1) - r(:,2)*const(:)
       s1(:,1,1) = s1(:,1,1) - s1(:,2,1)*const(:)
       s1(:,1,2) = s1(:,1,2) - s1(:,2,2)*const(:)
       s2(:,1,1) = s2(:,1,1) - s2(:,2,1)*const(:)
       s2(:,1,2) = s2(:,1,2) - s2(:,2,2)*const(:)

    end if ! n > 1

    ! All dependence has been shifted to the boundary elements
    ! Communicate boundary values to root process
    ! and solve reduced 11-diagonal system
    ! Use of 11-diagonal system is more robust than the
    ! reordered (removes zeros) 7-diagonal system

    ! Send rows of 11-diagonal system
    ! (0, 0, 0, a, b, c, 0, 0, 0, d, e; r)
    !    (0, 0, a, b, 0, c, 0, 0, d, e, 0; r)
    !       (0, a, b, 0, 0, c, 0, d, e, 0, 0; r)
    !          (a, b, 0, 0, 0, c, d, e, 0, 0, 0; r)
    ! For efficiency, only send non-zero elements

    L = 1
    k1 = 1
    do igroup=1,proc
       k2 = k1+ngroup(igroup)-1
       nk = k2-k1+1

       if (n > 1) then ! do normal stuff

          do i=1,2
             sendbuf(1:nk, L+0) = c(k1:k2, i)
             sendbuf(1:nk, L+1) = r(k1:k2, i)
             sendbuf(1:nk, L+2) = c(k1:k2, n-2+i)
             sendbuf(1:nk, L+3) = r(k1:k2, n-2+i)
             L = L + 4
          end do
          do i=1,2
             do j=1,2
                sendbuf(1:nk, L+0) = s1(k1:k2, i,j)
                sendbuf(1:nk, L+1) = s2(k1:k2, i,j)
                sendbuf(1:nk, L+2) = s1(k1:k2, n-2+i,j)
                sendbuf(1:nk, L+3) = s2(k1:k2, n-2+i,j)
                L = L + 4
             end do
          end do

       else ! n == 1 special case

          sendbuf(1:nk, L+0) = c(k1:k2, 1)
          sendbuf(1:nk, L+1) = r(k1:k2, 1)
          sendbuf(1:nk, L+2) = 1.0_WP
          sendbuf(1:nk, L+3) = 0.0_WP
          sendbuf(1:nk, L+4) = 1.0_WP
          sendbuf(1:nk, L+5) = 0.0_WP
          sendbuf(1:nk, L+6) = c(k1:k2, 1)
          sendbuf(1:nk, L+7) = r(k1:k2, 1)
          sendbuf(1:nk, L+8) = a(k1:k2, 1)
          sendbuf(1:nk, L+9) = d(k1:k2, 1)
          sendbuf(1:nk, L+10) = 0.0_WP
          sendbuf(1:nk, L+11) = 0.0_WP
          sendbuf(1:nk, L+12) = b(k1:k2, 1)
          sendbuf(1:nk, L+13) = e(k1:k2, 1)
          sendbuf(1:nk, L+14) = -1.0_WP
          sendbuf(1:nk, L+15) = 0.0_WP
          sendbuf(1:nk, L+16) = 0.0_WP
          sendbuf(1:nk, L+17) = -1.0_WP
          sendbuf(1:nk, L+18) = a(k1:k2, 1)
          sendbuf(1:nk, L+19) = d(k1:k2, 1)
          sendbuf(1:nk, L+20) = 0.0_WP
          sendbuf(1:nk, L+21) = 0.0_WP
          sendbuf(1:nk, L+22) = b(k1:k2, 1)
          sendbuf(1:nk, L+23) = e(k1:k2, 1)
          L = L + 24

       end if

       k1 = k2 + 1
    end do

    ! Gather the boundary data
    call MPI_AllToAll (sendbuf,nlot*24,MPI_REAL_WP,recvbuf,nlot*24,MPI_REAL_WP,ncom,ierr)

    ! Build reduced matrix
    allocate(AA(nlot, 4*proc, -5:5 + 1))
    AA = 0.0_WP
    L = 1
    do k=1,proc
       m = 4*(k-1)

       do i=1,2
          AA(:, m+i  , 0) = recvbuf(:, L+0)    ! c(i)
          AA(:, m+i  , 6) = recvbuf(:, L+1)    ! r(i)
          AA(:, m+i+2, 0) = recvbuf(:, L+2)    ! c(n-2+i)
          AA(:, m+i+2, 6) = recvbuf(:, L+3)    ! r(n-2+i)
          L = L + 4
       end do
       do i=1,2
          do j=1,2
             AA(:, m+i,   -2+j-i) = recvbuf(:, L+0)    ! s1(i,j)
             AA(:, m+i,    4+j-i) = recvbuf(:, L+1)    ! s2(i,j)
             AA(:, m+i+2, -4+j-i) = recvbuf(:, L+2)    ! s1(n-2+i,j)
             AA(:, m+i+2,  2+j-i) = recvbuf(:, L+3)    ! s2(n-2+i,j)
             L = L + 4
          end do
       end do

    end do

    ! Clear unused values
    nk = ngroup(nrank+1)
    AA(nk+1:nlot, :, :) = 0.0_WP
    AA(nk+1:nlot, :, 0) = 1.0_WP

    ! Solve reduced systems
    if (.not.nper) then
       call polydiagonal_serial(5, AA(:,3:4*proc-2,-5:+5), AA(:,3:4*proc-2,6), (2*proc-2)*2, nlot)
    else
       if (proc.lt.3) stop 'Pentadiagonal solver needs more proc for this problem'
       call polydiagonal_periodic_serial(5, AA(:,:,-5:+5), AA(:,:,6), (2*proc)*2, nlot, sendbuf)
    end if

    ! Move solution to beginning of recvbuf
    recvbuf(:, 1:4*proc) = AA(:, :, 6)
    deallocate(AA)

    ! Permute the order
    allocate (swap(nlot,2))
    do i=1,proc-1
       swap(:,:) = recvbuf(:, 4*i-1:4*i)
       recvbuf(:, 4*i-1:4*i) = recvbuf(:, 4*i+1:4*i+2)
       recvbuf(:, 4*i+1:4*i+2) = swap(:,:)
    end do
    ! If periodic, don't forget the end points
    if (nper) then
       swap(:,:) = recvbuf(:, 4*proc-1:4*proc)
       recvbuf(:, 4*proc-1:4*proc) = recvbuf(:, 1:2)
       recvbuf(:, 1:2) = swap(:,:)
    end if
    deallocate (swap)

    ! Scatter back the solution
    deallocate(sendbuf); allocate(sendbuf(nlot,4*proc))
    call MPI_AllToAll(recvbuf,nlot*4,MPI_REAL_WP,sendbuf,nlot*4,MPI_REAL_WP,ncom,ierr)

    L = 1
    k1 = 1
    do igroup=1,proc
       k2 = k1+ngroup(igroup)-1
       nk = k2-k1+1

       r1(k1:k2, :) = sendbuf(1:nk, L+0:L+1)
       r2(k1:k2, :) = sendbuf(1:nk, L+2:L+3)
       L = L + 4

       k1 = k2 + 1
    end do

    ! Only if not periodic
    if (.not.nper) then
       if (nrank .eq. 0)      r1 = 0.0_WP
       if (nrank .eq. proc-1) r2 = 0.0_WP
    end if

    if (n > 1) then ! do normal stuff

       do j=1,2
          do i=1,n
             r(:,i) = r(:,i) - s1(:,i,j)*r1(:,j) - s2(:,i,j)*r2(:,j)
          end do
       end do
       r = r / c

    else ! n == 1 special case

       r(:,1) = ( r(:,1) - a(:,1)*r1(:,1) - b(:,1)*r1(:,2) &
            -d(:,1)*r2(:,1) - e(:,1)*r2(:,2) )/c(:,1)

    end if

    deallocate(sendbuf)
    deallocate(recvbuf)
    deallocate(ngroup)

    return
  end subroutine pentadiagonal


  ! ================================================= !
  ! PentaDiagonal Solver - Serial Case - Not periodic !
  ! ================================================= !
  subroutine pentadiagonal_serial(A,B,C,D,E,R,n,lot)
    implicit none

    integer, intent(in) :: n,lot
    real(WP), dimension(lot,n) :: A     ! LOWER-2
    real(WP), dimension(lot,n) :: B     ! LOWER-1
    real(WP), dimension(lot,n) :: C     ! DIAGONAL
    real(WP), dimension(lot,n) :: D     ! UPPER+1
    real(WP), dimension(lot,n) :: E     ! UPPER+2
    real(WP), dimension(lot,n) :: R     ! RHS - RESULT
    real(WP), dimension(lot) :: const
    integer :: i

    if (n .eq. 1) then
       ! Solve 1x1 system
       R(:,1) = R(:,1)/C(:,1)
       return
    else if (n .eq. 2) then
       ! Solve 2x2 system
       const(:) = B(:,2)/C(:,1)
       C(:,2) = C(:,2) - D(:,1)*const(:)
       R(:,2) = R(:,2) - R(:,1)*const(:)
       R(:,2) = R(:,2)/C(:,2)
       R(:,1) = (R(:,1) - D(:,1)*R(:,2))/C(:,1)
       return
    end if

    ! Forward elimination
    do i=1,n-2
       ! Eliminate A(2,i+1)
       const(:) = B(:,i+1)/(C(:,i)+tiny(1.0_WP))
       C(:,i+1) = C(:,i+1) - D(:,i)*const(:)
       D(:,i+1) = D(:,i+1) - E(:,i)*const(:)
       R(:,i+1) = R(:,i+1) - R(:,i)*const(:)

       ! Eliminate A(1,i+2)
       const(:) = A(:,i+2)/(C(:,i)+tiny(1.0_WP))
       B(:,i+2) = B(:,i+2) - D(:,i)*const(:)
       C(:,i+2) = C(:,i+2) - E(:,i)*const(:)
       R(:,i+2) = R(:,i+2) - R(:,i)*const(:)
    end do
    ! Eliminate A(2,n)
    const(:) = B(:,n)/(C(:,n-1)+tiny(1.0_WP))
    C(:,n) = C(:,n) - D(:,n-1)*const(:)
    R(:,n) = R(:,n) - R(:,n-1)*const(:)

    ! Back-substitution
    R(:,n) = R(:,n)/(C(:,n)+tiny(1.0_WP))
    R(:,n-1) = (R(:,n-1) - D(:,n-1)*R(:,n))/(C(:,n-1)+tiny(1.0_WP))
    do i=n-2,1,-1
       R(:,i) = (R(:,i) - D(:,i)*R(:,i+1) - E(:,i)*R(:,i+2))/(C(:,i)+tiny(1.0_WP))
    end do

    return
  end subroutine pentadiagonal_serial



  ! ============================================= !
  ! PentaDiagonal Solver - Serial Case - Periodic !
  ! ============================================= !
  subroutine pentadiagonal_periodic_serial(A,B,C,D,E,R,n,lot,s1,s2)
    implicit none

    integer, intent(in) :: n,lot
    real(WP), dimension(lot,n) :: A     ! LOWER-2
    real(WP), dimension(lot,n) :: B     ! LOWER-1
    real(WP), dimension(lot,n) :: C     ! DIAGONAL
    real(WP), dimension(lot,n) :: D     ! UPPER+1
    real(WP), dimension(lot,n) :: E     ! UPPER+2
    real(WP), dimension(lot,n) :: R     ! RHS - RESULT
    real(WP), dimension(lot,n) :: s1
    real(WP), dimension(lot,n) :: s2
    real(WP), dimension(lot) :: const
    integer :: i

    if (n .eq. 1) then
       ! Solve 1x1 system
       R(:,:) = R(:,:)/(A(:,:)+B(:,:)+C(:,:)+D(:,:)+E(:,:))
       return
    else if (n .eq. 2) then
       ! Solve 2x2 system
       C(:,:) = C(:,:) + A(:,:) + E(:,:)
       D(:,1) = D(:,1) + B(:,1)
       B(:,2) = B(:,2) + D(:,2)
       const(:) = B(:,2)/C(:,1)
       C(:,2) = C(:,2) - D(:,1)*const(:)
       R(:,2) = R(:,2) - R(:,1)*const(:)
       R(:,2) = R(:,2)/C(:,2)
       R(:,1) = (R(:,1) - D(:,1)*R(:,2))/C(:,1)
       return
    else if (n .eq. 3) then
       B(:,:) = B(:,:) + E(:,:)
       D(:,:) = D(:,:) + A(:,:)
       call tridiagonal_periodic_serial(B(:,:), C(:,:), D(:,:), R(:,:), n, lot)
       return
    else if (n .eq. 4) then
       A(:,:) = A(:,:) + E(:,:)
       E(:,:) = 0.0_WP
    end if

    ! Initialize boundary data
    s1 = 0.0_WP
    s1(:,1) = A(:,1)
    s1(:,n-3) = s1(:,n-3) + E(:,n-3)
    s1(:,n-2) = D(:,n-2)
    s1(:,n-1) = C(:,n-1)
    s1(:,n) = B(:,n)
    s2 = 0.0_WP
    s2(:,1) = B(:,1)
    s2(:,2) = A(:,2)
    s2(:,n-2) = s2(:,n-2) + E(:,n-2)
    s2(:,n-1) = D(:,n-1)
    s2(:,n) = C(:,n)

    ! Forward elimination
    do i=1,n-2
       ! Eliminate b(i+1)
       const(:) = B(:,i+1)/C(:,i)
       C(:,i+1) = C(:,i+1) - D(:,i)*const(:)
       D(:,i+1) = D(:,i+1) - E(:,i)*const(:)
       R(:,i+1) = R(:,i+1) - R(:,i)*const(:)
       s1(:,i+1) = s1(:,i+1) - s1(:,i)*const(:)
       s2(:,i+1) = s2(:,i+1) - s2(:,i)*const(:)

       ! Eliminate a(i+2)
       const(:) = A(:,i+2)/C(:,i)
       B(:,i+2) = B(:,i+2) - D(:,i)*const(:)
       C(:,i+2) = C(:,i+2) - E(:,i)*const(:)
       R(:,i+2) = R(:,i+2) - R(:,i)*const(:)
       s1(:,i+2) = s1(:,i+2) - s1(:,i)*const(:)
       s2(:,i+2) = s2(:,i+2) - s2(:,i)*const(:)
    end do

    ! Backward elimination
    do i=n-2,3,-1
       ! Eliminate d(i-1)
       const(:) = D(:,i-1)/C(:,i)
       R(:,i-1) = R(:,i-1) - R(:,i)*const(:)
       s1(:,i-1) = s1(:,i-1) - s1(:,i)*const(:)
       s2(:,i-1) = s2(:,i-1) - s2(:,i)*const(:)

       ! Eliminate e(i-2)
       const(:) = E(:,i-2)/C(:,i)
       R(:,i-2) = R(:,i-2) - R(:,i)*const(:)
       s1(:,i-2) = s1(:,i-2) - s1(:,i)*const(:)
       s2(:,i-2) = s2(:,i-2) - s2(:,i)*const(:)
    end do
    i=2
    ! Eliminate d(i-1)
    const(:) = D(:,i-1)/C(:,i)
    R(:,i-1) = R(:,i-1) - R(:,i)*const(:)
    s1(:,i-1) = s1(:,i-1) - s1(:,i)*const(:)
    s2(:,i-1) = s2(:,i-1) - s2(:,i)*const(:)

    ! Eliminate oddball region
    const(:) = E(:,n-1)/C(:,1)
    R(:,n-1) = R(:,n-1) - R(:,1)*const(:)
    s1(:,n-1) = s1(:,n-1) - s1(:,1)*const(:)
    s2(:,n-1) = s2(:,n-1) - s2(:,1)*const(:)

    const(:) = D(:,n)/C(:,1)
    R(:,n) = R(:,n) - R(:,1)*const(:)
    s1(:,n) = s1(:,n) - s1(:,1)*const(:)
    s2(:,n) = s2(:,n) - s2(:,1)*const(:)

    const(:) = E(:,n)/C(:,2)
    R(:,n) = R(:,n) - R(:,2)*const(:)
    s1(:,n) = s1(:,n) - s1(:,2)*const(:)
    s2(:,n) = s2(:,n) - s2(:,2)*const(:)

    ! Eliminate corner region
    const(:) = s1(:,n)/s1(:,n-1)
    R(:,n) = R(:,n) - R(:,n-1)*const(:)
    s2(:,n) = s2(:,n) - s2(:,n-1)*const(:)

    R(:,n) = R(:,n)/s2(:,n)
    R(:,n-1) = (R(:,n-1) - s2(:,n-1)*R(:,n))/s1(:,n-1)
    do i=n-2,1,-1
       R(:,i) = (R(:,i) - s1(:,i)*R(:,n-1) - s2(:,i)*R(:,n))/C(:,i)
    end do

    return
  end subroutine pentadiagonal_periodic_serial


  subroutine polydiagonal(this,nd,A,R,n,lot,dir,b)
    implicit none
    class(diag), intent(inout) :: this
    character(len=*), intent(in) :: dir
    integer, intent(in) :: n,lot,nd
    real(WP), intent(inout), dimension(lot,n,-nd:nd) :: A
    real(WP), intent(inout), dimension(lot,n)        :: R
    real(WP), intent(inout), dimension(lot,n,nd)     :: b
    integer :: proc
    logical :: nper

    ! Get parallel info
    select case (trim(adjustl(dir)))
    case ('x')
       proc = this%cfg%npx
       nper = this%cfg%xper
    case ('y')
       proc = this%cfg%npy
       nper = this%cfg%yper
    case ('z')
       proc = this%cfg%npz
       nper = this%cfg%zper
    case default
       stop 'Unknown direction'
    end select

    ! If serial
    if (proc.eq.1) then
       if (.not.nper) then
          call polydiagonal_serial(nd,A,R,n,lot)
       else
          call polydiagonal_periodic_serial(nd,A,R,n,lot,b)
       end if
       return
    else
       stop 'polydiagonal: Implemented only in serial.'
    end if

    return
  end subroutine polydiagonal


  subroutine polydiagonal_serial(nd,A,R,n,lot)
    implicit none

    integer,  intent(in) :: n,lot,nd
    real(WP), intent(inout), dimension(lot,n,-nd:nd) :: A
    real(WP), intent(inout), dimension(lot,n)        :: R
    real(WP), dimension(lot) :: const,pivot
    integer :: i,j,k

    ! Forward elimination
    do i=1,n-1
       pivot(:) = 1.0_WP/A(:,i,0)
       do j=1,min(nd,n-i)
          const(:) = A(:,i+j,-j)*pivot(:)
          do k=1,min(nd,n-i)
             A(:,i+j,-j+k) = A(:,i+j,-j+k) - A(:,i,k)*const(:)
          end do
          R(:,i+j) = R(:,i+j) - R(:,i)*const(:)
       end do
    end do

    ! Back-substitution
    do i=n,1,-1
       do j=1,min(nd,n-i)
          R(:,i) = R(:,i) - A(:,i,j)*R(:,i+j)
       end do
       R(:,i) = R(:,i)/A(:,i,0)
    end do

    return
  end subroutine polydiagonal_serial


  subroutine polydiagonal_periodic_serial(nd,A,R,n,lot,b)
    implicit none

    integer,  intent(in) :: n,lot,nd
    real(WP), intent(inout), dimension(lot,n,-nd:nd) :: A
    real(WP), intent(inout), dimension(lot,n)        :: R
    real(WP), intent(inout), dimension(lot,n,nd)     :: b
    real(WP), dimension(lot) :: const,pivot
    integer :: i,j,k,i0

    if (n.eq.1) then
       const(:) = 0.0_WP
       do j=-nd,nd
          const(:) = const(:) + A(:,j,1)
       end do
       R(:,1) = R(:,1) / const(:)
       return
    end if

    ! Setup bridge array
    do j=1,nd
       do i=1,n
          b(:,i,j) = 0.0_WP
       end do
    end do
    do i=1,nd
       do j=i,nd
          b(:,i,j) = A(:,i,-nd+j-i)
          b(:,n-nd-i+1,j-i+1) = A(:,n-nd-i+1,j)
       end do
    end do
    do i=n-nd+1,n
       do j=1,nd
          b(:,i,j) = A(:,i,-nd+j-(i-n))
       end do
    end do

    ! Forward elimination
    do i=1,n-nd
       pivot(:) = 1.0_WP/A(:,i,0)
       do j=1,nd
          const(:) = A(:,i+j,-j)*pivot(:)
          do k=1,nd
             A(:,i+j,-j+k) = A(:,i+j,-j+k) - A(:,i,k)*const(:)
             b(:,i+j,k) = b(:,i+j,k) - b(:,i,k)*const(:)
          end do
          R(:,i+j) = R(:,i+j) - R(:,i)*const(:)
       end do
    end do

    ! Backward elimination
    do i=n-nd,1,-1
       pivot(:) = 1.0_WP/A(:,i,0)
       do j=1,min(nd,i-1)
          const(:) = A(:,i-j,j)*pivot(:)
          do k=1,nd
             b(:,i-j,k) = b(:,i-j,k) - b(:,i,k)*const(:)
          end do
          R(:,i-j) = R(:,i-j) - R(:,i)*const(:)
       end do
    end do

    ! Eliminate oddball region
    do i=1,nd
       pivot(:) = 1.0_WP/A(:,i,0)
       do j=i,nd
          const(:) = A(:,n-j+i,j)*pivot(:)
          do k=1,nd
             b(:,n-j+i,k) = b(:,n-j+i,k) - b(:,i,k)*const(:)
          end do
          R(:,n-j+i) = R(:,n-j+i) - R(:,i)*const(:)
       end do
    end do

    ! Elimination for corner matrix
    i0 = n-nd
    do i=1,nd-1
       pivot(:) = 1.0_WP/b(:,i+i0,i)
       do j=i+1,nd
          const(:) = b(:,j+i0,i)*pivot(:)
          do k=i+1,nd
             b(:,j+i0,k) = b(:,j+i0,k) - b(:,i+i0,k)*const(:)
          end do
          R(:,j+i0) = R(:,j+i0) - R(:,i+i0)*const(:)
       end do
    end do

    ! Back-substitution for corner matrix
    i0 = n-nd
    do i=nd,1,-1
       do j=i+1,nd
          R(:,i+i0) = R(:,i+i0) - b(:,i+i0,j)*R(:,j+i0)
       end do
       R(:,i+i0) = R(:,i+i0)/b(:,i+i0,i)
    end do

    ! Back-substitution for bridge
    do i=n-nd,1,-1
       do j=1,nd
          R(:,i) = R(:,i) - b(:,i,j)*R(:,n-nd+j)
       end do
       R(:,i) = R(:,i)/A(:,i,0)
    end do

    return
  end subroutine polydiagonal_periodic_serial


end module diag_class
