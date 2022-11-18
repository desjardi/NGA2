!> 3D FFT pressure solver concept is defined by extension of the linsol class
!> This solver is specifically intended to be a FFT-based pressure Poisson solver
!> for 3D periodic uniform computational domains decomposed in at most 2 directions
! Based on serial FFTW here, and manual transpose. It'd be nice to switch to a library
module pfft3d_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_short
   use linsol_class, only: linsol
   use, intrinsic :: iso_c_binding
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: pfft3d
   
   
   !> pfft3d object definition
   type, extends(linsol) :: pfft3d
      
      ! FFT's oddball
      logical :: oddball
      
      ! Data storage for FFTW plans
      real(C_DOUBLE), dimension(:), allocatable :: in_x,out_x
      real(C_DOUBLE), dimension(:), allocatable :: in_y,out_y
      real(C_DOUBLE), dimension(:), allocatable :: in_z,out_z
      
      ! FFTW plans
      type(C_PTR) :: fplan_x,bplan_x
      type(C_PTR) :: fplan_y,bplan_y
      type(C_PTR) :: fplan_z,bplan_z
      
      ! Storage for transposed data
      real(WP), dimension(:,:,:), allocatable :: xtrans
      real(WP), dimension(:,:,:), allocatable :: ytrans
      real(WP), dimension(:,:,:), allocatable :: ztrans
      
      ! Transpose partition - X
      integer, dimension(:), allocatable :: imin_x,imax_x
      integer, dimension(:), allocatable :: jmin_x,jmax_x
      integer, dimension(:), allocatable :: kmin_x,kmax_x
      integer, dimension(:), allocatable :: nx_x,ny_x,nz_x
      real(WP), dimension(:,:,:,:), allocatable :: sendbuf_x,recvbuf_x
      integer :: sendcount_x,recvcount_x
      character(len=str_short) :: xdir

      ! Transpose partition - Y
      integer, dimension(:), allocatable :: imin_y,imax_y
      integer, dimension(:), allocatable :: jmin_y,jmax_y
      integer, dimension(:), allocatable :: kmin_y,kmax_y
      integer, dimension(:), allocatable :: nx_y,ny_y,nz_y
      real(WP), dimension(:,:,:,:), allocatable :: sendbuf_y,recvbuf_y
      integer :: sendcount_y,recvcount_y
      character(len=str_short) :: ydir

      ! Transpose partition - Z
      integer, dimension(:), allocatable :: imin_z,imax_z
      integer, dimension(:), allocatable :: jmin_z,jmax_z
      integer, dimension(:), allocatable :: kmin_z,kmax_z
      integer, dimension(:), allocatable :: nx_z,ny_z,nz_z
      real(WP), dimension(:,:,:,:), allocatable :: sendbuf_z,recvbuf_z
      integer :: sendcount_z,recvcount_z
      character(len=str_short) :: zdir
      
   contains
      
      procedure :: print_short=>pfft3d_print_short    !< One-line printing of solver status
      procedure :: print=>pfft3d_print                !< Long-form printing of solver status
      procedure :: log=>pfft3d_log                    !< Long-form logging of solver status
      procedure :: init=>pfft3d_init                  !< Grid and stencil initialization - done once for the grid and stencil
      procedure :: setup=>pfft3d_setup                !< Solver setup (every time the operator changes)
      procedure :: solve=>pfft3d_solve                !< Execute solver (assumes new RHS and initial guess at every call)
      procedure :: destroy=>pfft3d_destroy            !< Solver destruction (every time the operator changes)
      
      procedure, private :: pfft3d_xtranspose_init
      procedure, private :: pfft3d_ytranspose_init
      procedure, private :: pfft3d_ztranspose_init
      
      procedure, private :: pfft3d_xtranspose_forward
      procedure, private :: pfft3d_ytranspose_forward
      procedure, private :: pfft3d_ztranspose_forward
      
      procedure, private :: pfft3d_xtranspose_backward
      procedure, private :: pfft3d_ytranspose_backward
      procedure, private :: pfft3d_ztranspose_backward
      
   end type pfft3d
   
   
   !> Declare pfft3d constructor
   interface pfft3d
      procedure pfft3d_from_args
   end interface pfft3d
   
   
contains
   
   
   !> Constructor for an pfft3d object
   function pfft3d_from_args(cfg,name) result(self)
      use messager, only: die
      implicit none
      type(pfft3d) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      
      ! Link the config and store the name
      self%cfg=>cfg
      self%name=trim(adjustl(name))
      
      ! Set solution method - not used
      self%method=0
      
      ! Set up stencil size and map
      self%nst=7
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
      
      ! Various checks to ensure we can use this solver
      check_solver_is_useable: block
         integer :: ndim,ndcp
         ! Periodicity and uniformity of mesh
         if (self%cfg%nx.gt.1.and..not.(self%cfg%xper.and.self%cfg%uniform_x)) call die('[pfft3d constructor] Need x-direction needs to be periodic and uniform')
         if (self%cfg%ny.gt.1.and..not.(self%cfg%yper.and.self%cfg%uniform_y)) call die('[pfft3d constructor] Need y-direction needs to be periodic and uniform')
         if (self%cfg%nz.gt.1.and..not.(self%cfg%zper.and.self%cfg%uniform_z)) call die('[pfft3d constructor] Need z-direction needs to be periodic and uniform')
         ! Ensure that we have at least one non-decomposed direction
         ndim=0
         if (self%cfg%nx.gt.1) ndim=ndim+1
         if (self%cfg%ny.gt.1) ndim=ndim+1
         if (self%cfg%nz.gt.1) ndim=ndim+1
         ndcp=0
         if (self%cfg%npx.gt.1) ndcp=ndcp+1
         if (self%cfg%npy.gt.1) ndcp=ndcp+1
         if (self%cfg%npz.gt.1) ndcp=ndcp+1
         if (ndcp.ge.ndim) call die('[pfft3d constructor] Need at least one NON-decomposed direction')
      end block check_solver_is_useable
      
   end function pfft3d_from_args
   
   
   !> Initialize grid and stencil - done at start-up only as long as the stencil or cfg does not change
   !> When calling, zero values of this%opr(1,:,:,:) indicate cells that do not participate in the solver
   !> Only the stencil needs to be defined at this point
   subroutine pfft3d_init(this)
      use messager, only: die
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: ierr,st,stx1,stx2,sty1,sty2,stz1,stz2
      integer, dimension(3) :: periodicity,offset
      include 'fftw3.f03'
      
      ! From the provided stencil, generate an inverse map
      stx1=minval(this%stc(:,1)); stx2=maxval(this%stc(:,1))
      sty1=minval(this%stc(:,2)); sty2=maxval(this%stc(:,2))
      stz1=minval(this%stc(:,3)); stz2=maxval(this%stc(:,3))
      allocate(this%stmap(stx1:stx2,sty1:sty2,stz1:stz2)); this%stmap=0
      do st=1,this%nst
         this%stmap(this%stc(st,1),this%stc(st,2),this%stc(st,3))=st
      end do
      
      ! Initialize transpose and FFTW plans in x
      if (this%cfg%nx.gt.1) then
         call this%pfft3d_xtranspose_init()
         allocate(this%in_x(this%cfg%nx),this%out_x(this%cfg%nx))
         this%fplan_x=fftw_plan_r2r_1d(this%cfg%nx,this%in_x,this%out_x,FFTW_R2HC,FFTW_MEASURE)
         this%bplan_x=fftw_plan_r2r_1d(this%cfg%nx,this%in_x,this%out_x,FFTW_HC2R,FFTW_MEASURE)
      end if
      
      ! Initialize transpose and FFTW plans in y
      if (this%cfg%ny.gt.1) then
         call this%pfft3d_ytranspose_init()
         allocate(this%in_y(this%cfg%ny),this%out_y(this%cfg%ny))
         this%fplan_y=fftw_plan_r2r_1d(this%cfg%ny,this%in_y,this%out_y,FFTW_R2HC,FFTW_MEASURE)
         this%bplan_y=fftw_plan_r2r_1d(this%cfg%ny,this%in_y,this%out_y,FFTW_HC2R,FFTW_MEASURE)
      end if
      
      ! Initialize transpose and FFTW plans in z
      if (this%cfg%nz.gt.1) then
         call this%pfft3d_ztranspose_init()
         allocate(this%in_z(this%cfg%nz),this%out_z(this%cfg%nz))
         this%fplan_z=fftw_plan_r2r_1d(this%cfg%nz,this%in_z,this%out_z,FFTW_R2HC,FFTW_MEASURE)
         this%bplan_z=fftw_plan_r2r_1d(this%cfg%nz,this%in_z,this%out_z,FFTW_HC2R,FFTW_MEASURE)
      end if
      
      ! Find who owns the oddball
      this%oddball=.false.
      if (this%cfg%iproc.eq.1.and.this%cfg%jproc.eq.1.and.this%cfg%kproc.eq.1) this%oddball=.true.
      
   end subroutine pfft3d_init
   
   
   !> Initialize transpose tool in x
   subroutine pfft3d_xtranspose_init(this)
      use mpi_f08
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: ierr,ip,q,r
      
      ! Determine non-decomposed direction to use for transpose
      if      (this%cfg%npx.eq.1.and.this%cfg%nx.gt.1) then
         this%xdir='x'
      else if (this%cfg%npy.eq.1.and.this%cfg%ny.gt.1) then
         this%xdir='y'
      else if (this%cfg%npz.eq.1.and.this%cfg%nz.gt.1) then
         this%xdir='z'
      end if
      
      ! Allocate global partitions
      allocate(  this%nx_x(this%cfg%npx))
      allocate(  this%ny_x(this%cfg%npx))
      allocate(  this%nz_x(this%cfg%npx))
      allocate(this%imin_x(this%cfg%npx))
      allocate(this%imax_x(this%cfg%npx))
      allocate(this%jmin_x(this%cfg%npx))
      allocate(this%jmax_x(this%cfg%npx))
      allocate(this%kmin_x(this%cfg%npx))
      allocate(this%kmax_x(this%cfg%npx))
      
      ! Partition
      select case (trim(this%xdir))
      case ('x')
         
         ! No transpose required, use local partition
         this%nx_x=this%cfg%nx_
         this%ny_x=this%cfg%ny_
         this%nz_x=this%cfg%nz_
         this%imin_x=this%cfg%imin_
         this%imax_x=this%cfg%imax_
         this%jmin_x=this%cfg%jmin_
         this%jmax_x=this%cfg%jmax_
         this%kmin_x=this%cfg%kmin_
         this%kmax_x=this%cfg%kmax_
         
      case ('y')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%cfg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
         call MPI_AllGather(this%cfg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
         this%nx_x=this%imax_x-this%imin_x+1
         
         ! Partition new local indices
         do ip=1,this%cfg%npx
            q=this%cfg%ny/this%cfg%npx
            r=mod(this%cfg%ny,this%cfg%npx)
            if (ip.le.r) then
               this%ny_x(ip)  =q+1
               this%jmin_x(ip)=this%cfg%jmin+(ip-1)*(q+1)
            else
               this%ny_x(ip)  =q
               this%jmin_x(ip)=this%cfg%jmin+r*(q+1)+(ip-r-1)*q
            end if
            this%jmax_x(ip)=this%jmin_x(ip)+this%ny_x(ip)-1
         end do
         this%nz_x=this%cfg%nz_
         this%kmin_x=this%cfg%kmin_
         this%kmax_x=this%cfg%kmax_
         
         ! Variables for AllToAll communication
         this%sendcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%cfg%nz_
         this%recvcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%cfg%nz_
         allocate(this%sendbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npx))
         allocate(this%recvbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npx))
         
         ! Zero out buffers
         this%sendbuf_x=0.0_WP
         this%recvbuf_x=0.0_WP
         
      case ('z')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%cfg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
         call MPI_AllGather(this%cfg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%cfg%xcomm,ierr)
         this%nx_x=this%imax_x-this%imin_x+1
         
         ! Partition new local indices
         do ip=1,this%cfg%npx
            q=this%cfg%nz/this%cfg%npx
            r=mod(this%cfg%nz,this%cfg%npx)
            if (ip.le.r) then
               this%nz_x(ip)  =q+1
               this%kmin_x(ip)=this%cfg%kmin+(ip-1)*(q+1)
            else
               this%nz_x(ip)  =q
               this%kmin_x(ip)=this%cfg%kmin+r*(q+1)+(ip-r-1)*q
            end if
            this%kmax_x(ip)=this%kmin_x(ip)+this%nz_x(ip)-1
         end do
         this%ny_x=this%cfg%ny_
         this%jmin_x=this%cfg%jmin_
         this%jmax_x=this%cfg%jmax_
         
         ! Variables for AllToAll communication
         this%sendcount_x=maxval(this%nx_x)*this%cfg%ny_*maxval(this%nz_x)
         this%recvcount_x=maxval(this%nx_x)*this%cfg%ny_*maxval(this%nz_x)
         allocate(this%sendbuf_x(maxval(this%nx_x),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_x),this%cfg%npx))
         allocate(this%recvbuf_x(maxval(this%nx_x),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_x),this%cfg%npx))
         
         ! Zero out buffers
         this%sendbuf_x=0.0_WP
         this%recvbuf_x=0.0_WP
         
      end select
      
      ! Allocate storage
      allocate(this%xtrans(this%cfg%imin:this%cfg%imax,this%jmin_x(this%cfg%iproc):this%jmax_x(this%cfg%iproc),this%kmin_x(this%cfg%iproc):this%kmax_x(this%cfg%iproc)))
      
   end subroutine pfft3d_xtranspose_init
   
   
   !> Initialize transpose tool in y
   subroutine pfft3d_ytranspose_init(this)
      use mpi_f08
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: ierr,jp,q,r
      
      ! Determine non-decomposed direction to use for transpose
      if      (this%cfg%npy.eq.1.and.this%cfg%ny.gt.1) then
         this%ydir='y'
      else if (this%cfg%npz.eq.1.and.this%cfg%nz.gt.1) then
         this%ydir='z'
      else if (this%cfg%npx.eq.1.and.this%cfg%nx.gt.1) then
         this%ydir='x'
      end if
      
      ! Allocate global partitions
      allocate(  this%nx_y(this%cfg%npy))
      allocate(  this%ny_y(this%cfg%npy))
      allocate(  this%nz_y(this%cfg%npy))
      allocate(this%imin_y(this%cfg%npy))
      allocate(this%imax_y(this%cfg%npy))
      allocate(this%jmin_y(this%cfg%npy))
      allocate(this%jmax_y(this%cfg%npy))
      allocate(this%kmin_y(this%cfg%npy))
      allocate(this%kmax_y(this%cfg%npy))
      
      ! Partition
      select case (trim(this%ydir))
      case ('x')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%cfg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
         call MPI_AllGather(this%cfg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
         this%ny_y=this%jmax_y-this%jmin_y+1
         
         ! Partition new local indices
         do jp=1,this%cfg%npy
            q=this%cfg%nx/this%cfg%npy
            r=mod(this%cfg%nx,this%cfg%npy)
            if (jp.le.r) then
               this%nx_y(jp)  =q+1
               this%imin_y(jp)=this%cfg%imin+(jp-1)*(q+1)
            else
               this%nx_y(jp)  =q
               this%imin_y(jp)=this%cfg%imin+r*(q+1)+(jp-r-1)*q
            end if
            this%imax_y(jp)=this%imin_y(jp)+this%nx_y(jp)-1
         end do
         this%nz_y=this%cfg%nz_
         this%kmin_y=this%cfg%kmin_
         this%kmax_y=this%cfg%kmax_
         
         ! Variables for AllToAll communication
         this%sendcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%cfg%nz_
         this%recvcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%cfg%nz_
         allocate(this%sendbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npy))
         allocate(this%recvbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%cfg%kmin_:this%cfg%kmax_,this%cfg%npy))
         
         ! Zero out buffers
         this%sendbuf_y=0.0_WP
         this%recvbuf_y=0.0_WP
         
      case ('y')
         
         ! No transpose required, use local partition
         this%nx_y=this%cfg%nx_
         this%ny_y=this%cfg%ny_
         this%nz_y=this%cfg%nz_
         this%imin_y=this%cfg%imin_
         this%imax_y=this%cfg%imax_
         this%jmin_y=this%cfg%jmin_
         this%jmax_y=this%cfg%jmax_
         this%kmin_y=this%cfg%kmin_
         this%kmax_y=this%cfg%kmax_
         
      case ('z')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%cfg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
         call MPI_AllGather(this%cfg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%cfg%ycomm,ierr)
         this%ny_y=this%jmax_y-this%jmin_y+1
         
         ! Partition new local indices
         do jp=1,this%cfg%npy
            q=this%cfg%nz/this%cfg%npy
            r=mod(this%cfg%nz,this%cfg%npy)
            if (jp.le.r) then
               this%nz_y(jp)  =q+1
               this%kmin_y(jp)=this%cfg%kmin+(jp-1)*(q+1)
            else
               this%nz_y(jp)  =q
               this%kmin_y(jp)=this%cfg%kmin+r*(q+1)+(jp-r-1)*q
            end if
            this%kmax_y(jp)=this%kmin_y(jp)+this%nz_y(jp)-1
         end do
         this%nx_y=this%cfg%nx_
         this%imin_y=this%cfg%imin_
         this%imax_y=this%cfg%imax_
         
         ! Variables for AllToAll communication
         this%sendcount_y=this%cfg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
         this%recvcount_y=this%cfg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
         allocate(this%sendbuf_y(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%cfg%npy))
         allocate(this%recvbuf_y(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%cfg%npy))
         
         ! Zero out buffers
         this%sendbuf_y=0.0_WP
         this%recvbuf_y=0.0_WP
         
      end select
      
      ! Allocate storage
      allocate(this%ytrans(this%imin_y(this%cfg%jproc):this%imax_y(this%cfg%jproc),this%cfg%jmin:this%cfg%jmax,this%kmin_y(this%cfg%jproc):this%kmax_y(this%cfg%jproc)))
      
   end subroutine pfft3d_ytranspose_init
   
   
   !> Initialize transpose tool in z
   subroutine pfft3d_ztranspose_init(this)
      use mpi_f08
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: ierr,kp,q,r
      
      ! Determine non-decomposed direction to use for transpose
      if      (this%cfg%npz.eq.1.and.this%cfg%nz.gt.1) then
         this%zdir='z'
      else if (this%cfg%npx.eq.1.and.this%cfg%nx.gt.1) then
         this%zdir='x'
      else if (this%cfg%npy.eq.1.and.this%cfg%ny.gt.1) then
         this%zdir='y'
      end if
      
      ! Allocate global partitions
      allocate(  this%nx_z(this%cfg%npz))
      allocate(  this%ny_z(this%cfg%npz))
      allocate(  this%nz_z(this%cfg%npz))
      allocate(this%imin_z(this%cfg%npz))
      allocate(this%imax_z(this%cfg%npz))
      allocate(this%jmin_z(this%cfg%npz))
      allocate(this%jmax_z(this%cfg%npz))
      allocate(this%kmin_z(this%cfg%npz))
      allocate(this%kmax_z(this%cfg%npz))
      
      ! Partition
      select case (trim(this%zdir))
      case ('x')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%cfg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
         call MPI_AllGather(this%cfg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
         this%nz_z=this%kmax_z-this%kmin_z+1
         
         ! Partition new local indices
         do kp=1,this%cfg%npz
            q=this%cfg%nx/this%cfg%npz
            r=mod(this%cfg%nx,this%cfg%npz)
            if (kp.le.r) then
               this%nx_z(kp)  =q+1
               this%imin_z(kp)=this%cfg%imin+(kp-1)*(q+1)
            else
               this%nx_z(kp)  =q
               this%imin_z(kp)=this%cfg%imin+r*(q+1)+(kp-r-1)*q
            end if
            this%imax_z(kp)=this%imin_z(kp)+this%nx_z(kp)-1
         end do
         this%ny_z=this%cfg%ny_
         this%jmin_z=this%cfg%jmin_
         this%jmax_z=this%cfg%jmax_
         
         ! Variables for AllToAll communication
         this%sendcount_z=maxval(this%nx_z)*this%cfg%ny_*maxval(this%nz_z)
         this%recvcount_z=maxval(this%nx_z)*this%cfg%ny_*maxval(this%nz_z)
         allocate(this%sendbuf_z(maxval(this%nx_z),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_z),this%cfg%npz))
         allocate(this%recvbuf_z(maxval(this%nx_z),this%cfg%jmin_:this%cfg%jmax_,maxval(this%nz_z),this%cfg%npz))
         
         ! Zero out buffers
         this%sendbuf_z=0.0_WP
         this%recvbuf_z=0.0_WP
         
      case ('y')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%cfg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
         call MPI_AllGather(this%cfg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%cfg%zcomm,ierr)
         this%nz_z=this%kmax_z-this%kmin_z+1
         
         ! Partition new local indices
         do kp=1,this%cfg%npz
            q=this%cfg%ny/this%cfg%npz
            r=mod(this%cfg%ny,this%cfg%npz)
            if (kp.le.r) then
               this%ny_z(kp)  =q+1
               this%jmin_z(kp)=this%cfg%jmin+(kp-1)*(q+1)
            else
               this%ny_z(kp)  =q
               this%jmin_z(kp)=this%cfg%jmin+r*(q+1)+(kp-r-1)*q
            end if
            this%jmax_z(kp)=this%jmin_z(kp)+this%ny_z(kp)-1
         end do
         this%nx_z=this%cfg%nx_
         this%imin_z=this%cfg%imin_
         this%imax_z=this%cfg%imax_
         
         ! Variables for AllToAll communication
         this%sendcount_z=this%cfg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
         this%recvcount_z=this%cfg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
         allocate(this%sendbuf_z(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%cfg%npz))
         allocate(this%recvbuf_z(this%cfg%imin_:this%cfg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%cfg%npz))
         
         ! Zero out buffers
         this%sendbuf_z=0.0_WP
         this%recvbuf_z=0.0_WP
         
      case ('z')
         
         ! No transpose required, use local partition
         this%nx_z=this%cfg%nx_
         this%ny_z=this%cfg%ny_
         this%nz_z=this%cfg%nz_
         this%imin_z=this%cfg%imin_
         this%imax_z=this%cfg%imax_
         this%jmin_z=this%cfg%jmin_
         this%jmax_z=this%cfg%jmax_
         this%kmin_z=this%cfg%kmin_
         this%kmax_z=this%cfg%kmax_
         
      end select
      
      ! Allocate storage
      allocate(this%ztrans(this%imin_z(this%cfg%kproc):this%imax_z(this%cfg%kproc),this%jmin_z(this%cfg%kproc):this%jmax_z(this%cfg%kproc),this%cfg%kmin:this%cfg%kmax))
      
   end subroutine pfft3d_ztranspose_init
   
   
   !> Perform forward transpose in x
   subroutine pfft3d_xtranspose_forward(this,A,At)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pfft3d), intent(inout) :: this
      real(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(in) :: A
      real(WP), dimension(this%cfg%imin :,this%jmin_x(this%cfg%iproc):,this%kmin_x(this%cfg%iproc):), intent(out) :: At
      integer :: i,j,k,ip,ii,jj,kk,ierr
      
      select case (trim(this%xdir))
      case ('x')
         ! No transpose required
         At=A
      case ('y')
         ! Transpose x=>y
         do ip=1,this%cfg%npx
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_x(ip),this%jmax_x(ip)
                  do i=this%cfg%imin_,this%cfg%imax_
                     jj=j-this%jmin_x(ip)+1
                     ii=i-this%cfg%imin_+1
                     this%sendbuf_x(ii,jj,k,ip)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_REAL_WP,this%recvbuf_x,this%recvcount_x,MPI_REAL_WP,this%cfg%xcomm,ierr)
         do ip=1,this%cfg%npx
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_x(this%cfg%iproc),this%jmax_x(this%cfg%iproc)
                  do i=this%imin_x(ip),this%imax_x(ip)
                     jj=j-this%jmin_x(this%cfg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     At(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
                  end do
               end do
            end do
         end do
      case ('z')
         ! Transpose x=>z
         do ip=1,this%cfg%npx
            do k=this%kmin_x(ip),this%kmax_x(ip)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%cfg%imin_,this%cfg%imax_
                     kk=k-this%kmin_x(ip)+1
                     ii=i-this%cfg%imin_+1
                     this%sendbuf_x(ii,j,kk,ip)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_REAL_WP,this%recvbuf_x,this%recvcount_x,MPI_REAL_WP,this%cfg%xcomm,ierr)
         do ip=1,this%cfg%npx
            do k=this%kmin_x(this%cfg%iproc),this%kmax_x(this%cfg%iproc)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_x(ip),this%imax_x(ip)
                     kk=k-this%kmin_x(this%cfg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     At(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine pfft3d_xtranspose_forward
   
   
   !> Perform forward transpose in y
   subroutine pfft3d_ytranspose_forward(this,A,At)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pfft3d), intent(inout) :: this
      real(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(in) :: A
      real(WP), dimension(this%imin_y(this%cfg%jproc):,this%cfg%jmin:,this%kmin_y(this%cfg%jproc):), intent(out) :: At
      integer :: i,j,k,jp,ii,jj,kk,ierr
      
      select case (trim(this%ydir))
      case ('x')
         ! Transpose y=>x
         do jp=1,this%cfg%npy
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_y(jp),this%imax_y(jp)
                     ii=i-this%imin_y(jp)+1
                     jj=j-this%cfg%jmin_+1
                     this%sendbuf_y(ii,jj,k,jp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_REAL_WP,this%recvbuf_y,this%recvcount_y,MPI_REAL_WP,this%cfg%ycomm,ierr)
         do jp=1,this%cfg%npy
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%imin_y(this%cfg%jproc),this%imax_y(this%cfg%jproc)
                     ii=i-this%imin_y(this%cfg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     At(i,j,k)=this%recvbuf_y(ii,jj,k,jp)
                  end do
               end do
            end do
         end do
      case ('y')
         ! No transpose required
         At=A
      case ('z')
         ! Transpose y=>z
         do jp=1,this%cfg%npy
            do k=this%kmin_y(jp),this%kmax_y(jp)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%cfg%imin_,this%cfg%imax_
                     kk=k-this%kmin_y(jp)+1
                     jj=j-this%cfg%jmin_+1
                     this%sendbuf_y(i,jj,kk,jp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_REAL_WP,this%recvbuf_y,this%recvcount_y,MPI_REAL_WP,this%cfg%ycomm,ierr)
         do jp=1,this%cfg%npy
            do k=this%kmin_y(this%cfg%jproc),this%kmax_y(this%cfg%jproc)
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%cfg%imin_,this%cfg%imax_
                     kk=k-this%kmin_y(this%cfg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     At(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine pfft3d_ytranspose_forward
   
   
   !> Perform forward transpose in z
   subroutine pfft3d_ztranspose_forward(this,A,At)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pfft3d), intent(inout) :: this
      real(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(in) :: A
      real(WP), dimension(this%imin_z(this%cfg%kproc):,this%jmin_z(this%cfg%kproc):,this%cfg%kmin:), intent(out) :: At
      integer :: i,j,k,kp,ii,jj,kk,ierr
      
      select case (trim(this%zdir))
      case ('x')
         ! Transpose z=>x
         do kp=1,this%cfg%npz
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_z(kp),this%imax_z(kp)
                     ii=i-this%imin_z(kp)+1
                     kk=k-this%cfg%kmin_+1
                     this%sendbuf_z(ii,j,kk,kp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_REAL_WP,this%recvbuf_z,this%recvcount_z,MPI_REAL_WP,this%cfg%zcomm,ierr)
         do kp=1,this%cfg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_z(this%cfg%kproc),this%imax_z(this%cfg%kproc)
                     ii=i-this%imin_z(this%cfg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     At(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
                  end do
               end do
            end do
         end do
      case ('y')
         ! Transpose z=>y
         do kp=1,this%cfg%npz
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_z(kp),this%jmax_z(kp)
                  do i=this%cfg%imin_,this%cfg%imax_
                     jj=j-this%jmin_z(kp)+1
                     kk=k-this%cfg%kmin_+1
                     this%sendbuf_z(i,jj,kk,kp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_REAL_WP,this%recvbuf_z,this%recvcount_z,MPI_REAL_WP,this%cfg%zcomm,ierr)
         do kp=1,this%cfg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%jmin_z(this%cfg%kproc),this%jmax_z(this%cfg%kproc)
                  do i=this%cfg%imin_,this%cfg%imax_
                     jj=j-this%jmin_z(this%cfg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     At(i,j,k)=this%recvbuf_z(i,jj,kk,kp)
                  end do
               end do
            end do
         end do
      case ('z')
         ! No transpose required
         At=A
      end select
      
   end subroutine pfft3d_ztranspose_forward
   
   
   !> Perform backward transpose in x
   subroutine pfft3d_xtranspose_backward(this,At,A)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pfft3d), intent(inout) :: this
      real(WP), dimension(this%cfg%imin :,this%jmin_x(this%cfg%iproc):,this%kmin_x(this%cfg%iproc):), intent(in) :: At
      real(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(out) :: A
      integer :: i,j,k,ip,ii,jj,kk,ierr
      
      select case (trim(this%xdir))
      case ('x')
         ! No transpose required
         A=At
      case ('y')
         ! Transpose y=>x
         do ip=1,this%cfg%npx
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_x(this%cfg%iproc),this%jmax_x(this%cfg%iproc)
                  do i=this%imin_x(ip),this%imax_x(ip)
                     jj=j-this%jmin_x(this%cfg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     this%sendbuf_x(ii,jj,k,ip)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_REAL_WP,this%recvbuf_x,this%recvcount_x,MPI_REAL_WP,this%cfg%xcomm,ierr)
         do ip=1,this%cfg%npx
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_x(ip),this%jmax_x(ip)
                  do i=this%imin_x(this%cfg%iproc),this%imax_x(this%cfg%iproc)
                     jj=j-this%jmin_x(ip)+1
                     ii=i-this%imin_x(this%cfg%iproc)+1
                     A(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
                  end do
               end do
            end do
         end do
      case ('z')
         ! Transpose z=>x
         do ip=1,this%cfg%npx
            do k=this%kmin_x(this%cfg%iproc),this%kmax_x(this%cfg%iproc)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_x(ip),this%imax_x(ip)
                     kk=k-this%kmin_x(this%cfg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     this%sendbuf_x(ii,j,kk,ip)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_REAL_WP,this%recvbuf_x,this%recvcount_x,MPI_REAL_WP,this%cfg%xcomm,ierr)
         do ip=1,this%cfg%npx
            do k=this%kmin_x(ip),this%kmax_x(ip)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_x(this%cfg%iproc),this%imax_x(this%cfg%iproc)
                     kk=k-this%kmin_x(ip)+1
                     ii=i-this%imin_x(this%cfg%iproc)+1
                     A(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine pfft3d_xtranspose_backward

   
   !> Perform backward transpose in y
   subroutine pfft3d_ytranspose_backward(this,At,A)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pfft3d), intent(inout) :: this
      real(WP), dimension(this%imin_y(this%cfg%jproc):,this%cfg%jmin:,this%kmin_y(this%cfg%jproc):), intent(in) :: At
      real(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(out) :: A
      integer :: i,j,k,jp,ii,jj,kk,ierr
      
      select case (trim(this%ydir))
      case ('x')
         ! Transpose x=>y
         do jp=1,this%cfg%npy
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%imin_y(this%cfg%jproc),this%imax_y(this%cfg%jproc)
                     ii=i-this%imin_y(this%cfg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     this%sendbuf_y(ii,jj,k,jp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_REAL_WP,this%recvbuf_y,this%recvcount_y,MPI_REAL_WP,this%cfg%ycomm,ierr)
         do jp=1,this%cfg%npy
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%jmin_y(this%cfg%jproc),this%jmax_y(this%cfg%jproc)
                  do i=this%imin_y(jp),this%imax_y(jp)
                     ii=i-this%imin_y(jp)+1
                     jj=j-this%jmin_y(this%cfg%jproc)+1
                     A(i,j,k)=this%recvbuf_y(ii,jj,k,jp)
                  end do
               end do
            end do
         end do
      case ('y')
         ! No transpose required
         A=At
      case ('z')
         ! Transpose z=>y
         do jp=1,this%cfg%npy
            do k=this%kmin_y(this%cfg%jproc),this%kmax_y(this%cfg%jproc)
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%cfg%imin_,this%cfg%imax_
                     kk=k-this%kmin_y(this%cfg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     this%sendbuf_y(i,jj,kk,jp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_REAL_WP,this%recvbuf_y,this%recvcount_y,MPI_REAL_WP,this%cfg%ycomm,ierr)
         do jp=1,this%cfg%npy
            do k=this%kmin_y(jp),this%kmax_y(jp)
               do j=this%jmin_y(this%cfg%jproc),this%jmax_y(this%cfg%jproc)
                  do i=this%cfg%imin_,this%cfg%imax_
                     kk=k-this%kmin_y(jp)+1
                     jj=j-this%jmin_y(this%cfg%jproc)+1
                     A(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine pfft3d_ytranspose_backward
   
   
   !> Perform backward transpose in z
   subroutine pfft3d_ztranspose_backward(this,At,A)
      use mpi_f08
      use parallel, only: MPI_REAL_WP
      implicit none
      class(pfft3d), intent(inout) :: this
      real(WP), dimension(this%imin_z(this%cfg%kproc):,this%jmin_z(this%cfg%kproc):,this%cfg%kmin:), intent(in) :: At
      real(WP), dimension(this%cfg%imin_:,this%cfg%jmin_:,this%cfg%kmin_:), intent(out) :: A
      integer :: i,j,k,kp,ii,jj,kk,ierr
      
      select case (trim(this%zdir))
      case ('x')
         ! Transpose x=>z
         do kp=1,this%cfg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_z(this%cfg%kproc),this%imax_z(this%cfg%kproc)
                     ii=i-this%imin_z(this%cfg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     this%sendbuf_z(ii,j,kk,kp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_REAL_WP,this%recvbuf_z,this%recvcount_z,MPI_REAL_WP,this%cfg%zcomm,ierr)
         do kp=1,this%cfg%npz
            do k=this%kmin_z(this%cfg%kproc),this%kmax_z(this%cfg%kproc)
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%imin_z(kp),this%imax_z(kp)
                     ii=i-this%imin_z(kp)+1
                     kk=k-this%kmin_z(this%cfg%kproc)+1
                     A(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
                  end do
               end do
            end do
         end do
      case ('y')
         ! Transpose y=>z
         do kp=1,this%cfg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%jmin_z(this%cfg%kproc),this%jmax_z(this%cfg%kproc)
                  do i=this%cfg%imin_,this%cfg%imax_
                     jj=j-this%jmin_z(this%cfg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     this%sendbuf_z(i,jj,kk,kp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_REAL_WP,this%recvbuf_z,this%recvcount_z,MPI_REAL_WP,this%cfg%zcomm,ierr)
         do kp=1,this%cfg%npz
            do k=this%kmin_z(this%cfg%kproc),this%kmax_z(this%cfg%kproc)
               do j=this%jmin_z(kp),this%jmax_z(kp)
                  do i=this%cfg%imin_,this%cfg%imax_
                     jj=j-this%jmin_z(kp)+1
                     kk=k-this%kmin_z(this%cfg%kproc)+1
                     A(i,j,k)=this%recvbuf_z(i,jj,kk,kp)
                  end do
               end do
            end do
         end do
      case ('z')
         ! No transpose required
         A=At
      end select
      
   end subroutine pfft3d_ztranspose_backward
   
   
   !> Setup solver - done everytime the operator changes
   subroutine pfft3d_setup(this)
      use messager, only: die
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: i,j,k,st,ierr
      integer,  dimension(:), allocatable :: row
      real(WP), dimension(:), allocatable :: val
      
      ! If the solver has already been setup, destroy it first
      if (this%setup_done) call this%destroy()
      
      
      
      ! Set setup-flag to true
      this%setup_done=.true.
      
   end subroutine pfft3d_setup
   
   
   !> Solve the linear system iteratively
   subroutine pfft3d_solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: i,j,k,ierr
      
      ! Check that setup was done
      if (.not.this%setup_done) call die('[pfft3d solve] Solver has not been setup.')
      
      ! Set solver it and err to standard values
      this%it=-1; this%aerr=huge(1.0_WP); this%rerr=huge(1.0_WP)
      
      
      
      ! Sync the solution vector
      call this%cfg%sync(this%sol)
      
      ! If verbose run, log and or print info
      if (verbose.gt.0) call this%log
      if (verbose.gt.1) call this%print_short
      
   end subroutine pfft3d_solve
   
   
   !> Destroy solver - done everytime the operator changes
   subroutine pfft3d_destroy(this)
      use messager, only: die
      implicit none
      class(pfft3d), intent(inout) :: this
      integer :: ierr
      
      ! Destroy solver, operator, and rhs/sol vectors
      
      
      ! Set setup-flag to false
      this%setup_done=.false.
      
   end subroutine pfft3d_destroy
   
   
   !> Log pfft3d info
   subroutine pfft3d_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(pfft3d), intent(in) :: this
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("PFFT3D solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
      end if
   end subroutine pfft3d_log
   
   
   !> Print pfft3d info to the screen
   subroutine pfft3d_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(pfft3d), intent(in) :: this
      if (this%cfg%amRoot) then
         write(output_unit,'("PFFT3D solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
   end subroutine pfft3d_print
   
   
   !> Short print of pfft3d info to the screen
   subroutine pfft3d_print_short(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(pfft3d), intent(in) :: this
      if (this%cfg%amRoot) write(output_unit,'("PFFT3D solver [",a16,"] for config [",a16,"]")') trim(this%name),trim(this%cfg%name)
   end subroutine pfft3d_print_short
   
   
end module pfft3d_class
