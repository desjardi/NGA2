!> FFT for periodic uniform computational domains decomposed in at most 2 directions.
!> Makes use of complex to complex FFTW and in-house parallel transpose operations.
!>
!> Unlike FFTW (and several other libaries), the transforms provided by this class have the correct
!> scaling and sign. Forward and backward are truly inverses no further multiplication is necessary.
!> Real-space integral and Fourier-space zero value is preserved through forward & inverse transforms.
module cfourier_class
   use precision,    only: WP
   use pgrid_class,  only: pgrid
   use string,       only: str_short
   use, intrinsic :: iso_c_binding
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: cfourier
   
   
   !> cfourier object definition
   type :: cfourier
      
      !> Pointer to our pgrid
      class(pgrid), pointer :: pg
      
      !> FFT's oddball
      logical :: oddball
      
      !> Available FFTs
      logical :: xfft_avail
      logical :: yfft_avail
      logical :: zfft_avail
      
      !> Data storage for FFTW plans
      complex(WP), dimension(:), allocatable :: in_x,out_x
      complex(WP), dimension(:), allocatable :: in_y,out_y
      complex(WP), dimension(:), allocatable :: in_z,out_z
      
      !> FFTW plans
      type(C_PTR) :: fplan_x,bplan_x
      type(C_PTR) :: fplan_y,bplan_y
      type(C_PTR) :: fplan_z,bplan_z
      
      !> Storage for transposed data
      complex(WP), dimension(:,:,:), allocatable :: xtrans,ytrans,ztrans
      
      !> Transpose partition - X
      integer, dimension(:), allocatable :: imin_x,imax_x
      integer, dimension(:), allocatable :: jmin_x,jmax_x
      integer, dimension(:), allocatable :: kmin_x,kmax_x
      integer, dimension(:), allocatable :: nx_x,ny_x,nz_x
      complex(WP), dimension(:,:,:,:), allocatable :: sendbuf_x,recvbuf_x
      integer :: sendcount_x,recvcount_x
      character :: xdir
      
      !> Transpose partition - Y
      integer, dimension(:), allocatable :: imin_y,imax_y
      integer, dimension(:), allocatable :: jmin_y,jmax_y
      integer, dimension(:), allocatable :: kmin_y,kmax_y
      integer, dimension(:), allocatable :: nx_y,ny_y,nz_y
      complex(WP), dimension(:,:,:,:), allocatable :: sendbuf_y,recvbuf_y
      integer :: sendcount_y,recvcount_y
      character :: ydir
      
      !> Transpose partition - Z
      integer, dimension(:), allocatable :: imin_z,imax_z
      integer, dimension(:), allocatable :: jmin_z,jmax_z
      integer, dimension(:), allocatable :: kmin_z,kmax_z
      integer, dimension(:), allocatable :: nx_z,ny_z,nz_z
      complex(WP), dimension(:,:,:,:), allocatable :: sendbuf_z,recvbuf_z
      integer :: sendcount_z,recvcount_z
      character :: zdir
      
   contains
      
      procedure :: xtransform_forward                 !< Forward Fourier transform in x
      procedure :: ytransform_forward                 !< Forward Fourier transform in y
      procedure :: ztransform_forward                 !< Forward Fourier transform in z
      
      procedure :: xtransform_backward                !< Backward Fourier transform in x
      procedure :: ytransform_backward                !< Backward Fourier transform in y
      procedure :: ztransform_backward                !< Backward Fourier transform in z
      
      procedure, private :: xtranspose_init           !< Transpose initialization in x
      procedure, private :: ytranspose_init           !< Transpose initialization in y
      procedure, private :: ztranspose_init           !< Transpose initialization in z
      
      procedure, private :: xtranspose_forward        !< Forward transpose in x
      procedure, private :: ytranspose_forward        !< Forward transpose in y
      procedure, private :: ztranspose_forward        !< Forward transpose in z
      
      procedure, private :: xtranspose_backward       !< Backward transpose in x
      procedure, private :: ytranspose_backward       !< Backward transpose in y
      procedure, private :: ztranspose_backward       !< Backward transpose in z
      
      procedure :: print=>cfourier_print               !< Long-form printing of transform status
      procedure :: log=>cfourier_log                   !< Long-form logging of transform status
      final :: cfourier_destroy                        !< Destructor
      
   end type cfourier
   
   
   !> Declare cfourier constructor
   interface cfourier
      procedure cfourier_from_args
   end interface cfourier
   
   
contains
   
   
   !> Constructor for a cfourier object
   function cfourier_from_args(pg) result(self)
      use messager, only: die
      use param,    only: verbose
      implicit none
      type(cfourier) :: self
      class(pgrid), target, intent(in) :: pg
      include 'fftw3.f03'
      
      ! Link the config
      self%pg=>pg
      
      ! Various checks to ensure we can use this solver
      check_solver_is_useable: block
         integer :: ndim,ndcp
         ! Check periodicity and uniformity of mesh per direction
         self%xfft_avail=.true.; if (self%pg%nx.gt.1.and..not.(self%pg%xper.and.self%pg%uniform_x)) self%xfft_avail=.false.
         self%yfft_avail=.true.; if (self%pg%ny.gt.1.and..not.(self%pg%yper.and.self%pg%uniform_y)) self%yfft_avail=.false.
         self%zfft_avail=.true.; if (self%pg%nz.gt.1.and..not.(self%pg%zper.and.self%pg%uniform_z)) self%zfft_avail=.false.
         if (.not.any([self%xfft_avail,self%yfft_avail,self%zfft_avail])) call die('[fft3d constructor] At least one direction needs to be periodic and uniform')
         ! Ensure that we have at least one non-decomposed direction
         ndim=count([self%pg%nx ,self%pg%ny ,self%pg%nz ].gt.1)
         ndcp=count([self%pg%npx,self%pg%npy,self%pg%npz].gt.1)
         if (ndcp.ge.ndim) call die('[fft3d constructor] Need at least one NON-decomposed direction')
      end block check_solver_is_useable
      
      ! Initialize transpose and FFTW plans in x
      if (self%pg%nx.gt.1.and.self%xfft_avail) then
         call self%xtranspose_init()
         allocate(self%in_x(self%pg%nx),self%out_x(self%pg%nx))
         self%fplan_x=fftw_plan_dft_1d(self%pg%nx,self%in_x,self%out_x,FFTW_FORWARD,FFTW_MEASURE)
         self%bplan_x=fftw_plan_dft_1d(self%pg%nx,self%in_x,self%out_x,FFTW_BACKWARD,FFTW_MEASURE)
      end if
      
      ! Initialize transpose and FFTW plans in y
      if (self%pg%ny.gt.1.and.self%yfft_avail) then
         call self%ytranspose_init()
         allocate(self%in_y(self%pg%ny),self%out_y(self%pg%ny))
         self%fplan_y=fftw_plan_dft_1d(self%pg%ny,self%in_y,self%out_y,FFTW_FORWARD,FFTW_MEASURE)
         self%bplan_y=fftw_plan_dft_1d(self%pg%ny,self%in_y,self%out_y,FFTW_BACKWARD,FFTW_MEASURE)
      end if
      
      ! Initialize transpose and FFTW plans in z
      if (self%pg%nz.gt.1.and.self%zfft_avail) then
         call self%ztranspose_init()
         allocate(self%in_z(self%pg%nz),self%out_z(self%pg%nz))
         self%fplan_z=fftw_plan_dft_1d(self%pg%nz,self%in_z,self%out_z,FFTW_FORWARD,FFTW_MEASURE)
         self%bplan_z=fftw_plan_dft_1d(self%pg%nz,self%in_z,self%out_z,FFTW_BACKWARD,FFTW_MEASURE)
      end if
      
      ! Find which process owns the oddball, if any
      self%oddball=.false.
      if (all([self%xfft_avail,self%yfft_avail,self%zfft_avail]).and.&
      &   all([self%pg%iproc,self%pg%jproc,self%pg%kproc].eq.1)) self%oddball=.true.
      
      ! If verbose run, log and or print grid
      if (verbose.gt.1) call self%log()
      if (verbose.gt.2) call self%print()
      
   end function cfourier_from_args
   
   
   !> Initialize transpose tool in x
   subroutine xtranspose_init(this)
      use mpi_f08, only: MPI_AllGather,MPI_INTEGER
      implicit none
      class(cfourier), intent(inout) :: this
      integer :: ierr,ip,q,r
      
      ! Determine non-decomposed direction to use for transpose
      if      (this%pg%npx.eq.1.and.this%pg%nx.gt.1) then
         this%xdir='x'
      else if (this%pg%npy.eq.1.and.this%pg%ny.gt.1) then
         this%xdir='y'
      else if (this%pg%npz.eq.1.and.this%pg%nz.gt.1) then
         this%xdir='z'
      end if
      
      ! Allocate global partitions
      allocate(  this%nx_x(this%pg%npx))
      allocate(  this%ny_x(this%pg%npx))
      allocate(  this%nz_x(this%pg%npx))
      allocate(this%imin_x(this%pg%npx))
      allocate(this%imax_x(this%pg%npx))
      allocate(this%jmin_x(this%pg%npx))
      allocate(this%jmax_x(this%pg%npx))
      allocate(this%kmin_x(this%pg%npx))
      allocate(this%kmax_x(this%pg%npx))
      
      ! Partition
      select case (trim(this%xdir))
      case ('x')
         
         ! No transpose required, use local partition
         this%nx_x=this%pg%nx_
         this%ny_x=this%pg%ny_
         this%nz_x=this%pg%nz_
         this%imin_x=this%pg%imin_
         this%imax_x=this%pg%imax_
         this%jmin_x=this%pg%jmin_
         this%jmax_x=this%pg%jmax_
         this%kmin_x=this%pg%kmin_
         this%kmax_x=this%pg%kmax_
         
      case ('y')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%pg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
         call MPI_AllGather(this%pg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
         this%nx_x=this%imax_x-this%imin_x+1
         
         ! Partition new local indices
         do ip=1,this%pg%npx
            q=this%pg%ny/this%pg%npx
            r=mod(this%pg%ny,this%pg%npx)
            if (ip.le.r) then
               this%ny_x(ip)  =q+1
               this%jmin_x(ip)=this%pg%jmin+(ip-1)*(q+1)
            else
               this%ny_x(ip)  =q
               this%jmin_x(ip)=this%pg%jmin+r*(q+1)+(ip-r-1)*q
            end if
            this%jmax_x(ip)=this%jmin_x(ip)+this%ny_x(ip)-1
         end do
         this%nz_x=this%pg%nz_
         this%kmin_x=this%pg%kmin_
         this%kmax_x=this%pg%kmax_
         
         ! Variables for AllToAll communication
         this%sendcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%pg%nz_
         this%recvcount_x=maxval(this%nx_x)*maxval(this%ny_x)*this%pg%nz_
         allocate(this%sendbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%pg%kmin_:this%pg%kmax_,this%pg%npx))
         allocate(this%recvbuf_x(maxval(this%nx_x),maxval(this%ny_x),this%pg%kmin_:this%pg%kmax_,this%pg%npx))
         
         ! Zero out buffers
         this%sendbuf_x=0.0_WP
         this%recvbuf_x=0.0_WP
         
      case ('z')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%pg%imin_,1,MPI_INTEGER,this%imin_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
         call MPI_AllGather(this%pg%imax_,1,MPI_INTEGER,this%imax_x,1,MPI_INTEGER,this%pg%xcomm,ierr)
         this%nx_x=this%imax_x-this%imin_x+1
         
         ! Partition new local indices
         do ip=1,this%pg%npx
            q=this%pg%nz/this%pg%npx
            r=mod(this%pg%nz,this%pg%npx)
            if (ip.le.r) then
               this%nz_x(ip)  =q+1
               this%kmin_x(ip)=this%pg%kmin+(ip-1)*(q+1)
            else
               this%nz_x(ip)  =q
               this%kmin_x(ip)=this%pg%kmin+r*(q+1)+(ip-r-1)*q
            end if
            this%kmax_x(ip)=this%kmin_x(ip)+this%nz_x(ip)-1
         end do
         this%ny_x=this%pg%ny_
         this%jmin_x=this%pg%jmin_
         this%jmax_x=this%pg%jmax_
         
         ! Variables for AllToAll communication
         this%sendcount_x=maxval(this%nx_x)*this%pg%ny_*maxval(this%nz_x)
         this%recvcount_x=maxval(this%nx_x)*this%pg%ny_*maxval(this%nz_x)
         allocate(this%sendbuf_x(maxval(this%nx_x),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_x),this%pg%npx))
         allocate(this%recvbuf_x(maxval(this%nx_x),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_x),this%pg%npx))
         
         ! Zero out buffers
         this%sendbuf_x=0.0_WP
         this%recvbuf_x=0.0_WP
         
      end select
      
      ! Allocate storage
      allocate(this%xtrans(this%pg%imin:this%pg%imax,this%jmin_x(this%pg%iproc):this%jmax_x(this%pg%iproc),this%kmin_x(this%pg%iproc):this%kmax_x(this%pg%iproc)))
      
   end subroutine xtranspose_init
   
   
   !> Initialize transpose tool in y
   subroutine ytranspose_init(this)
      use mpi_f08, only: MPI_AllGather,MPI_INTEGER
      implicit none
      class(cfourier), intent(inout) :: this
      integer :: ierr,jp,q,r
      
      ! Determine non-decomposed direction to use for transpose
      if      (this%pg%npy.eq.1.and.this%pg%ny.gt.1) then
         this%ydir='y'
      else if (this%pg%npz.eq.1.and.this%pg%nz.gt.1) then
         this%ydir='z'
      else if (this%pg%npx.eq.1.and.this%pg%nx.gt.1) then
         this%ydir='x'
      end if
      
      ! Allocate global partitions
      allocate(  this%nx_y(this%pg%npy))
      allocate(  this%ny_y(this%pg%npy))
      allocate(  this%nz_y(this%pg%npy))
      allocate(this%imin_y(this%pg%npy))
      allocate(this%imax_y(this%pg%npy))
      allocate(this%jmin_y(this%pg%npy))
      allocate(this%jmax_y(this%pg%npy))
      allocate(this%kmin_y(this%pg%npy))
      allocate(this%kmax_y(this%pg%npy))
      
      ! Partition
      select case (trim(this%ydir))
      case ('x')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%pg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
         call MPI_AllGather(this%pg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
         this%ny_y=this%jmax_y-this%jmin_y+1
         
         ! Partition new local indices
         do jp=1,this%pg%npy
            q=this%pg%nx/this%pg%npy
            r=mod(this%pg%nx,this%pg%npy)
            if (jp.le.r) then
               this%nx_y(jp)  =q+1
               this%imin_y(jp)=this%pg%imin+(jp-1)*(q+1)
            else
               this%nx_y(jp)  =q
               this%imin_y(jp)=this%pg%imin+r*(q+1)+(jp-r-1)*q
            end if
            this%imax_y(jp)=this%imin_y(jp)+this%nx_y(jp)-1
         end do
         this%nz_y=this%pg%nz_
         this%kmin_y=this%pg%kmin_
         this%kmax_y=this%pg%kmax_
         
         ! Variables for AllToAll communication
         this%sendcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%pg%nz_
         this%recvcount_y=maxval(this%nx_y)*maxval(this%ny_y)*this%pg%nz_
         allocate(this%sendbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%pg%kmin_:this%pg%kmax_,this%pg%npy))
         allocate(this%recvbuf_y(maxval(this%nx_y),maxval(this%ny_y),this%pg%kmin_:this%pg%kmax_,this%pg%npy))
         
         ! Zero out buffers
         this%sendbuf_y=0.0_WP
         this%recvbuf_y=0.0_WP
         
      case ('y')
         
         ! No transpose required, use local partition
         this%nx_y=this%pg%nx_
         this%ny_y=this%pg%ny_
         this%nz_y=this%pg%nz_
         this%imin_y=this%pg%imin_
         this%imax_y=this%pg%imax_
         this%jmin_y=this%pg%jmin_
         this%jmax_y=this%pg%jmax_
         this%kmin_y=this%pg%kmin_
         this%kmax_y=this%pg%kmax_
         
      case ('z')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%pg%jmin_,1,MPI_INTEGER,this%jmin_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
         call MPI_AllGather(this%pg%jmax_,1,MPI_INTEGER,this%jmax_y,1,MPI_INTEGER,this%pg%ycomm,ierr)
         this%ny_y=this%jmax_y-this%jmin_y+1
         
         ! Partition new local indices
         do jp=1,this%pg%npy
            q=this%pg%nz/this%pg%npy
            r=mod(this%pg%nz,this%pg%npy)
            if (jp.le.r) then
               this%nz_y(jp)  =q+1
               this%kmin_y(jp)=this%pg%kmin+(jp-1)*(q+1)
            else
               this%nz_y(jp)  =q
               this%kmin_y(jp)=this%pg%kmin+r*(q+1)+(jp-r-1)*q
            end if
            this%kmax_y(jp)=this%kmin_y(jp)+this%nz_y(jp)-1
         end do
         this%nx_y=this%pg%nx_
         this%imin_y=this%pg%imin_
         this%imax_y=this%pg%imax_
         
         ! Variables for AllToAll communication
         this%sendcount_y=this%pg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
         this%recvcount_y=this%pg%nx_*maxval(this%ny_y)*maxval(this%nz_y)
         allocate(this%sendbuf_y(this%pg%imin_:this%pg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%pg%npy))
         allocate(this%recvbuf_y(this%pg%imin_:this%pg%imax_,maxval(this%ny_y),maxval(this%nz_y),this%pg%npy))
         
         ! Zero out buffers
         this%sendbuf_y=0.0_WP
         this%recvbuf_y=0.0_WP
         
      end select
      
      ! Allocate storage
      allocate(this%ytrans(this%imin_y(this%pg%jproc):this%imax_y(this%pg%jproc),this%pg%jmin:this%pg%jmax,this%kmin_y(this%pg%jproc):this%kmax_y(this%pg%jproc)))
      
   end subroutine ytranspose_init
   
   
   !> Initialize transpose tool in z
   subroutine ztranspose_init(this)
      use mpi_f08, only: MPI_AllGather,MPI_INTEGER
      implicit none
      class(cfourier), intent(inout) :: this
      integer :: ierr,kp,q,r
      
      ! Determine non-decomposed direction to use for transpose
      if      (this%pg%npz.eq.1.and.this%pg%nz.gt.1) then
         this%zdir='z'
      else if (this%pg%npx.eq.1.and.this%pg%nx.gt.1) then
         this%zdir='x'
      else if (this%pg%npy.eq.1.and.this%pg%ny.gt.1) then
         this%zdir='y'
      end if
      
      ! Allocate global partitions
      allocate(  this%nx_z(this%pg%npz))
      allocate(  this%ny_z(this%pg%npz))
      allocate(  this%nz_z(this%pg%npz))
      allocate(this%imin_z(this%pg%npz))
      allocate(this%imax_z(this%pg%npz))
      allocate(this%jmin_z(this%pg%npz))
      allocate(this%jmax_z(this%pg%npz))
      allocate(this%kmin_z(this%pg%npz))
      allocate(this%kmax_z(this%pg%npz))
      
      ! Partition
      select case (trim(this%zdir))
      case ('x')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%pg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
         call MPI_AllGather(this%pg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
         this%nz_z=this%kmax_z-this%kmin_z+1
         
         ! Partition new local indices
         do kp=1,this%pg%npz
            q=this%pg%nx/this%pg%npz
            r=mod(this%pg%nx,this%pg%npz)
            if (kp.le.r) then
               this%nx_z(kp)  =q+1
               this%imin_z(kp)=this%pg%imin+(kp-1)*(q+1)
            else
               this%nx_z(kp)  =q
               this%imin_z(kp)=this%pg%imin+r*(q+1)+(kp-r-1)*q
            end if
            this%imax_z(kp)=this%imin_z(kp)+this%nx_z(kp)-1
         end do
         this%ny_z=this%pg%ny_
         this%jmin_z=this%pg%jmin_
         this%jmax_z=this%pg%jmax_
         
         ! Variables for AllToAll communication
         this%sendcount_z=maxval(this%nx_z)*this%pg%ny_*maxval(this%nz_z)
         this%recvcount_z=maxval(this%nx_z)*this%pg%ny_*maxval(this%nz_z)
         allocate(this%sendbuf_z(maxval(this%nx_z),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_z),this%pg%npz))
         allocate(this%recvbuf_z(maxval(this%nx_z),this%pg%jmin_:this%pg%jmax_,maxval(this%nz_z),this%pg%npz))
         
         ! Zero out buffers
         this%sendbuf_z=0.0_WP
         this%recvbuf_z=0.0_WP
         
      case ('y')
         
         ! Store old local indices from each processor
         call MPI_AllGather(this%pg%kmin_,1,MPI_INTEGER,this%kmin_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
         call MPI_AllGather(this%pg%kmax_,1,MPI_INTEGER,this%kmax_z,1,MPI_INTEGER,this%pg%zcomm,ierr)
         this%nz_z=this%kmax_z-this%kmin_z+1
         
         ! Partition new local indices
         do kp=1,this%pg%npz
            q=this%pg%ny/this%pg%npz
            r=mod(this%pg%ny,this%pg%npz)
            if (kp.le.r) then
               this%ny_z(kp)  =q+1
               this%jmin_z(kp)=this%pg%jmin+(kp-1)*(q+1)
            else
               this%ny_z(kp)  =q
               this%jmin_z(kp)=this%pg%jmin+r*(q+1)+(kp-r-1)*q
            end if
            this%jmax_z(kp)=this%jmin_z(kp)+this%ny_z(kp)-1
         end do
         this%nx_z=this%pg%nx_
         this%imin_z=this%pg%imin_
         this%imax_z=this%pg%imax_
         
         ! Variables for AllToAll communication
         this%sendcount_z=this%pg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
         this%recvcount_z=this%pg%nx_*maxval(this%ny_z)*maxval(this%nz_z)
         allocate(this%sendbuf_z(this%pg%imin_:this%pg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%pg%npz))
         allocate(this%recvbuf_z(this%pg%imin_:this%pg%imax_,maxval(this%ny_z),maxval(this%nz_z),this%pg%npz))
         
         ! Zero out buffers
         this%sendbuf_z=0.0_WP
         this%recvbuf_z=0.0_WP
         
      case ('z')
         
         ! No transpose required, use local partition
         this%nx_z=this%pg%nx_
         this%ny_z=this%pg%ny_
         this%nz_z=this%pg%nz_
         this%imin_z=this%pg%imin_
         this%imax_z=this%pg%imax_
         this%jmin_z=this%pg%jmin_
         this%jmax_z=this%pg%jmax_
         this%kmin_z=this%pg%kmin_
         this%kmax_z=this%pg%kmax_
         
      end select
      
      ! Allocate storage
      allocate(this%ztrans(this%imin_z(this%pg%kproc):this%imax_z(this%pg%kproc),this%jmin_z(this%pg%kproc):this%jmax_z(this%pg%kproc),this%pg%kmin:this%pg%kmax))
      
   end subroutine ztranspose_init
   
   
   !> Perform forward transpose in x
   subroutine xtranspose_forward(this,A,At)
      use mpi_f08,  only: MPI_AllToAll
      use parallel, only: MPI_COMPLEX_WP
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(in) :: A
      complex(WP), dimension(this%pg%imin :,this%jmin_x(this%pg%iproc):,this%kmin_x(this%pg%iproc):), intent(out) :: At
      integer :: i,j,k,ip,ii,jj,kk,ierr
      
      select case (trim(this%xdir))
      case ('x')
         ! No transpose required
         At=A
      case ('y')
         ! Transpose x=>y
         do ip=1,this%pg%npx
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_x(ip),this%jmax_x(ip)
                  do i=this%pg%imin_,this%pg%imax_
                     jj=j-this%jmin_x(ip)+1
                     ii=i-this%pg%imin_+1
                     this%sendbuf_x(ii,jj,k,ip)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_COMPLEX_WP,this%recvbuf_x,this%recvcount_x,MPI_COMPLEX_WP,this%pg%xcomm,ierr)
         do ip=1,this%pg%npx
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
                  do i=this%imin_x(ip),this%imax_x(ip)
                     jj=j-this%jmin_x(this%pg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     At(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
                  end do
               end do
            end do
         end do
      case ('z')
         ! Transpose x=>z
         do ip=1,this%pg%npx
            do k=this%kmin_x(ip),this%kmax_x(ip)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%pg%imin_,this%pg%imax_
                     kk=k-this%kmin_x(ip)+1
                     ii=i-this%pg%imin_+1
                     this%sendbuf_x(ii,j,kk,ip)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_COMPLEX_WP,this%recvbuf_x,this%recvcount_x,MPI_COMPLEX_WP,this%pg%xcomm,ierr)
         do ip=1,this%pg%npx
            do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_x(ip),this%imax_x(ip)
                     kk=k-this%kmin_x(this%pg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     At(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine xtranspose_forward
   
   
   !> Perform forward transpose in y
   subroutine ytranspose_forward(this,A,At)
      use mpi_f08,  only: MPI_AllToAll
      use parallel, only: MPI_COMPLEX_WP
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(in) :: A
      complex(WP), dimension(this%imin_y(this%pg%jproc):,this%pg%jmin:,this%kmin_y(this%pg%jproc):), intent(out) :: At
      integer :: i,j,k,jp,ii,jj,kk,ierr
      
      select case (trim(this%ydir))
      case ('x')
         ! Transpose y=>x
         do jp=1,this%pg%npy
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_y(jp),this%imax_y(jp)
                     ii=i-this%imin_y(jp)+1
                     jj=j-this%pg%jmin_+1
                     this%sendbuf_y(ii,jj,k,jp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_COMPLEX_WP,this%recvbuf_y,this%recvcount_y,MPI_COMPLEX_WP,this%pg%ycomm,ierr)
         do jp=1,this%pg%npy
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
                     ii=i-this%imin_y(this%pg%jproc)+1
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
         do jp=1,this%pg%npy
            do k=this%kmin_y(jp),this%kmax_y(jp)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%pg%imin_,this%pg%imax_
                     kk=k-this%kmin_y(jp)+1
                     jj=j-this%pg%jmin_+1
                     this%sendbuf_y(i,jj,kk,jp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_COMPLEX_WP,this%recvbuf_y,this%recvcount_y,MPI_COMPLEX_WP,this%pg%ycomm,ierr)
         do jp=1,this%pg%npy
            do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%pg%imin_,this%pg%imax_
                     kk=k-this%kmin_y(this%pg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     At(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine ytranspose_forward
   
   
   !> Perform forward transpose in z
   subroutine ztranspose_forward(this,A,At)
      use mpi_f08,  only: MPI_AllToAll
      use parallel, only: MPI_COMPLEX_WP
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(in) :: A
      complex(WP), dimension(this%imin_z(this%pg%kproc):,this%jmin_z(this%pg%kproc):,this%pg%kmin:), intent(out) :: At
      integer :: i,j,k,kp,ii,jj,kk,ierr
      
      select case (trim(this%zdir))
      case ('x')
         ! Transpose z=>x
         do kp=1,this%pg%npz
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_z(kp),this%imax_z(kp)
                     ii=i-this%imin_z(kp)+1
                     kk=k-this%pg%kmin_+1
                     this%sendbuf_z(ii,j,kk,kp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_COMPLEX_WP,this%recvbuf_z,this%recvcount_z,MPI_COMPLEX_WP,this%pg%zcomm,ierr)
         do kp=1,this%pg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
                     ii=i-this%imin_z(this%pg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     At(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
                  end do
               end do
            end do
         end do
      case ('y')
         ! Transpose z=>y
         do kp=1,this%pg%npz
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_z(kp),this%jmax_z(kp)
                  do i=this%pg%imin_,this%pg%imax_
                     jj=j-this%jmin_z(kp)+1
                     kk=k-this%pg%kmin_+1
                     this%sendbuf_z(i,jj,kk,kp)=A(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_COMPLEX_WP,this%recvbuf_z,this%recvcount_z,MPI_COMPLEX_WP,this%pg%zcomm,ierr)
         do kp=1,this%pg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
                  do i=this%pg%imin_,this%pg%imax_
                     jj=j-this%jmin_z(this%pg%kproc)+1
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
      
   end subroutine ztranspose_forward
   
   
   !> Perform backward transpose in x
   subroutine xtranspose_backward(this,At,A)
      use mpi_f08,  only: MPI_AllToAll
      use parallel, only: MPI_COMPLEX_WP
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin :,this%jmin_x(this%pg%iproc):,this%kmin_x(this%pg%iproc):), intent(in) :: At
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(out) :: A
      integer :: i,j,k,ip,ii,jj,kk,ierr
      
      select case (trim(this%xdir))
      case ('x')
         ! No transpose required
         A=At
      case ('y')
         ! Transpose y=>x
         do ip=1,this%pg%npx
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
                  do i=this%imin_x(ip),this%imax_x(ip)
                     jj=j-this%jmin_x(this%pg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     this%sendbuf_x(ii,jj,k,ip)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_COMPLEX_WP,this%recvbuf_x,this%recvcount_x,MPI_COMPLEX_WP,this%pg%xcomm,ierr)
         do ip=1,this%pg%npx
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_x(ip),this%jmax_x(ip)
                  do i=this%imin_x(this%pg%iproc),this%imax_x(this%pg%iproc)
                     jj=j-this%jmin_x(ip)+1
                     ii=i-this%imin_x(this%pg%iproc)+1
                     A(i,j,k)=this%recvbuf_x(ii,jj,k,ip)
                  end do
               end do
            end do
         end do
      case ('z')
         ! Transpose z=>x
         do ip=1,this%pg%npx
            do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_x(ip),this%imax_x(ip)
                     kk=k-this%kmin_x(this%pg%iproc)+1
                     ii=i-this%imin_x(ip)+1
                     this%sendbuf_x(ii,j,kk,ip)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_x,this%sendcount_x,MPI_COMPLEX_WP,this%recvbuf_x,this%recvcount_x,MPI_COMPLEX_WP,this%pg%xcomm,ierr)
         do ip=1,this%pg%npx
            do k=this%kmin_x(ip),this%kmax_x(ip)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_x(this%pg%iproc),this%imax_x(this%pg%iproc)
                     kk=k-this%kmin_x(ip)+1
                     ii=i-this%imin_x(this%pg%iproc)+1
                     A(i,j,k)=this%recvbuf_x(ii,j,kk,ip)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine xtranspose_backward
   
   
   !> Perform backward transpose in y
   subroutine ytranspose_backward(this,At,A)
      use mpi_f08,  only: MPI_AllToAll
      use parallel, only: MPI_COMPLEX_WP
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%imin_y(this%pg%jproc):,this%pg%jmin:,this%kmin_y(this%pg%jproc):), intent(in) :: At
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(out) :: A
      integer :: i,j,k,jp,ii,jj,kk,ierr
      
      select case (trim(this%ydir))
      case ('x')
         ! Transpose x=>y
         do jp=1,this%pg%npy
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
                     ii=i-this%imin_y(this%pg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     this%sendbuf_y(ii,jj,k,jp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_COMPLEX_WP,this%recvbuf_y,this%recvcount_y,MPI_COMPLEX_WP,this%pg%ycomm,ierr)
         do jp=1,this%pg%npy
            do k=this%pg%kmin_,this%pg%kmax_
               do j=this%jmin_y(this%pg%jproc),this%jmax_y(this%pg%jproc)
                  do i=this%imin_y(jp),this%imax_y(jp)
                     ii=i-this%imin_y(jp)+1
                     jj=j-this%jmin_y(this%pg%jproc)+1
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
         do jp=1,this%pg%npy
            do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
               do j=this%jmin_y(jp),this%jmax_y(jp)
                  do i=this%pg%imin_,this%pg%imax_
                     kk=k-this%kmin_y(this%pg%jproc)+1
                     jj=j-this%jmin_y(jp)+1
                     this%sendbuf_y(i,jj,kk,jp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_y,this%sendcount_y,MPI_COMPLEX_WP,this%recvbuf_y,this%recvcount_y,MPI_COMPLEX_WP,this%pg%ycomm,ierr)
         do jp=1,this%pg%npy
            do k=this%kmin_y(jp),this%kmax_y(jp)
               do j=this%jmin_y(this%pg%jproc),this%jmax_y(this%pg%jproc)
                  do i=this%pg%imin_,this%pg%imax_
                     kk=k-this%kmin_y(jp)+1
                     jj=j-this%jmin_y(this%pg%jproc)+1
                     A(i,j,k)=this%recvbuf_y(i,jj,kk,jp)
                  end do
               end do
            end do
         end do
      end select
      
   end subroutine ytranspose_backward
   
   
   !> Perform backward transpose in z
   subroutine ztranspose_backward(this,At,A)
      use mpi_f08,  only: MPI_AllToAll
      use parallel, only: MPI_COMPLEX_WP
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%imin_z(this%pg%kproc):,this%jmin_z(this%pg%kproc):,this%pg%kmin:), intent(in) :: At
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(out) :: A
      integer :: i,j,k,kp,ii,jj,kk,ierr
      
      select case (trim(this%zdir))
      case ('x')
         ! Transpose x=>z
         do kp=1,this%pg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
                     ii=i-this%imin_z(this%pg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     this%sendbuf_z(ii,j,kk,kp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_COMPLEX_WP,this%recvbuf_z,this%recvcount_z,MPI_COMPLEX_WP,this%pg%zcomm,ierr)
         do kp=1,this%pg%npz
            do k=this%kmin_z(this%pg%kproc),this%kmax_z(this%pg%kproc)
               do j=this%pg%jmin_,this%pg%jmax_
                  do i=this%imin_z(kp),this%imax_z(kp)
                     ii=i-this%imin_z(kp)+1
                     kk=k-this%kmin_z(this%pg%kproc)+1
                     A(i,j,k)=this%recvbuf_z(ii,j,kk,kp)
                  end do
               end do
            end do
         end do
      case ('y')
         ! Transpose y=>z
         do kp=1,this%pg%npz
            do k=this%kmin_z(kp),this%kmax_z(kp)
               do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
                  do i=this%pg%imin_,this%pg%imax_
                     jj=j-this%jmin_z(this%pg%kproc)+1
                     kk=k-this%kmin_z(kp)+1
                     this%sendbuf_z(i,jj,kk,kp)=At(i,j,k)
                  end do
               end do
            end do
         end do
         call MPI_AllToAll(this%sendbuf_z,this%sendcount_z,MPI_COMPLEX_WP,this%recvbuf_z,this%recvcount_z,MPI_COMPLEX_WP,this%pg%zcomm,ierr)
         do kp=1,this%pg%npz
            do k=this%kmin_z(this%pg%kproc),this%kmax_z(this%pg%kproc)
               do j=this%jmin_z(kp),this%jmax_z(kp)
                  do i=this%pg%imin_,this%pg%imax_
                     jj=j-this%jmin_z(kp)+1
                     kk=k-this%kmin_z(this%pg%kproc)+1
                     A(i,j,k)=this%recvbuf_z(i,jj,kk,kp)
                  end do
               end do
            end do
         end do
      case ('z')
         ! No transpose required
         A=At
      end select
      
   end subroutine ztranspose_backward
   
   
   !> Fourier transform A in x direction
   subroutine xtransform_forward(this,A)
      use messager, only: die
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
      integer :: j,k
      include 'fftw3.f03'
      ! Check availability of fft in x
      if (.not.this%xfft_avail) call die('[cfourier] Fourier transform is not available in x')
      ! Nothing to do if single cell
      if (this%pg%nx.eq.1) return
      ! Transpose in X
      call this%xtranspose_forward(A,this%xtrans)
      ! Forward transform - X
      do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
         do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
            this%in_x=this%xtrans(:,j,k)
            call fftw_execute_dft(this%fplan_x,this%in_x,this%out_x)
            this%xtrans(:,j,k)=this%out_x
         end do
      end do
      ! Transpose back
      call this%xtranspose_backward(this%xtrans,A)
   end subroutine xtransform_forward
   
   
   !> Fourier transform A in y direction
   subroutine ytransform_forward(this,A)
      use messager, only: die
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
      integer :: i,k
      include 'fftw3.f03'
      ! Check availability of fft in y
      if (.not.this%yfft_avail) call die('[cfourier] Fourier transform is not available in y')
      ! Nothing to do if single cell
      if (this%pg%ny.eq.1) return
      ! Transpose in Y
      call this%ytranspose_forward(A,this%ytrans)
      ! Forward transform - Y
      do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
         do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
            this%in_y=this%ytrans(i,:,k)
            call fftw_execute_dft(this%fplan_y,this%in_y,this%out_y)
            this%ytrans(i,:,k)=this%out_y
         end do
      end do
      ! Transpose back
      call this%ytranspose_backward(this%ytrans,A)
   end subroutine ytransform_forward
   
   
   !> Fourier transform A in z direction
   subroutine ztransform_forward(this,A)
      use messager, only: die
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
      integer :: i,j
      include 'fftw3.f03'
      ! Check availability of fft in z
      if (.not.this%zfft_avail) call die('[cfourier] Fourier transform is not available in z')
      ! Nothing to do if single cell
      if (this%pg%nz.eq.1) return
      ! Transpose in Z
      call this%ztranspose_forward(A,this%ztrans)
      ! Forward transform - Z
      do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
         do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
            this%in_z=this%ztrans(i,j,:)
            call fftw_execute_dft(this%fplan_z,this%in_z,this%out_z)
            this%ztrans(i,j,:)=this%out_z
         end do
      end do
      ! Transpose back
      call this%ztranspose_backward(this%ztrans,A)
   end subroutine ztransform_forward
   
   
   !> Transform A back from Fourier space in x direction
   subroutine xtransform_backward(this,A)
      use messager, only: die
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
      integer :: j,k
      include 'fftw3.f03'
      ! Check availability of fft in x
      if (.not.this%xfft_avail) call die('[cfourier] Fourier transform is not available in x')
      ! Nothing to do if single cell
      if (this%pg%nx.eq.1) return
      ! Transpose in X
      call this%xtranspose_forward(A,this%xtrans)
      ! Inverse transform
      do k=this%kmin_x(this%pg%iproc),this%kmax_x(this%pg%iproc)
         do j=this%jmin_x(this%pg%iproc),this%jmax_x(this%pg%iproc)
            this%in_x=this%xtrans(:,j,k)
            call fftw_execute_dft(this%bplan_x,this%in_x,this%out_x)
            this%xtrans(:,j,k)=this%out_x/real(this%pg%nx,WP)
         end do
      end do
      ! Transpose back
      call this%xtranspose_backward(this%xtrans,A)
   end subroutine xtransform_backward
   
   
   !> Transform A back from Fourier space in y direction
   subroutine ytransform_backward(this,A)
      use messager, only: die
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
      integer :: i,k
      include 'fftw3.f03'
      ! Check availability of fft in y
      if (.not.this%yfft_avail) call die('[cfourier] Fourier transform is not available in y')
      ! Nothing to do if single cell
      if (this%pg%ny.eq.1) return
      ! Transpose in Y
      call this%ytranspose_forward(A,this%ytrans)
      ! Inverse transform
      do k=this%kmin_y(this%pg%jproc),this%kmax_y(this%pg%jproc)
         do i=this%imin_y(this%pg%jproc),this%imax_y(this%pg%jproc)
            this%in_y=this%ytrans(i,:,k)
            call fftw_execute_dft(this%bplan_y,this%in_y,this%out_y)
            this%ytrans(i,:,k)=this%out_y/real(this%pg%ny,WP)
         end do
      end do
      ! Transpose back
      call this%ytranspose_backward(this%ytrans,A)
   end subroutine ytransform_backward
   
   
   !> Transform A back from Fourier space in z direction
   subroutine ztransform_backward(this,A)
      use messager, only: die
      implicit none
      class(cfourier), intent(inout) :: this
      complex(WP), dimension(this%pg%imin_:,this%pg%jmin_:,this%pg%kmin_:), intent(inout) :: A         !< Needs to be (imin_:imax_,jmin_:jmax_,kmin_:kmax_)
      integer :: i,j
      include 'fftw3.f03'
      ! Check availability of fft in z
      if (.not.this%zfft_avail) call die('[cfourier] Fourier transform is not available in z')
      ! Nothing to do if single cell
      if (this%pg%nz.eq.1) return
      ! Transpose in Z
      call this%ztranspose_forward(A,this%ztrans)
      ! Inverse transform
      do j=this%jmin_z(this%pg%kproc),this%jmax_z(this%pg%kproc)
         do i=this%imin_z(this%pg%kproc),this%imax_z(this%pg%kproc)
            this%in_z=this%ztrans(i,j,:)
            call fftw_execute_dft(this%bplan_z,this%in_z,this%out_z)
            this%ztrans(i,j,:)=this%out_z/real(this%pg%nz,WP)
         end do
      end do
      ! Transpose back
      call this%ztranspose_backward(this%ztrans,A)
   end subroutine ztransform_backward
   
   
   !> Destroy cfourier object
   subroutine cfourier_destroy(this)
      implicit none
      type(cfourier) :: this
      include 'fftw3.f03'
      if (this%pg%nx.gt.1.and.this%xfft_avail) then
         call fftw_destroy_plan(this%fplan_x)
         call fftw_destroy_plan(this%bplan_x)
         deallocate(this%in_x,this%out_x,this%xtrans,this%imin_x,this%imax_x,this%jmin_x,this%jmax_x,this%kmin_x,this%kmax_x,this%nx_x,this%ny_x,this%nz_x,this%sendbuf_x,this%recvbuf_x)
      end if
      if (this%pg%ny.gt.1.and.this%yfft_avail) then
         call fftw_destroy_plan(this%fplan_y)
         call fftw_destroy_plan(this%bplan_y)
         deallocate(this%in_y,this%out_y,this%ytrans,this%imin_y,this%imax_y,this%jmin_y,this%jmax_y,this%kmin_y,this%kmax_y,this%nx_y,this%ny_y,this%nz_y,this%sendbuf_y,this%recvbuf_y)
      end if
      if (this%pg%nz.gt.1.and.this%zfft_avail) then
         call fftw_destroy_plan(this%fplan_z)
         call fftw_destroy_plan(this%bplan_z)
         deallocate(this%in_z,this%out_z,this%ztrans,this%imin_z,this%imax_z,this%jmin_z,this%jmax_z,this%kmin_z,this%kmax_z,this%nx_z,this%ny_z,this%nz_z,this%sendbuf_z,this%recvbuf_z)
      end if
   end subroutine cfourier_destroy
   
   
   !> Log cfourier info
   subroutine cfourier_log(this)
      use string,   only: str_long
      use messager, only: log
      implicit none
      class(cfourier), intent(in) :: this
      character(len=str_long) :: message
      character(len=3) :: dir
      if (this%pg%amRoot) then
         dir=''
         if (this%xfft_avail) dir(1:1)='x'
         if (this%yfft_avail) dir(2:2)='y'
         if (this%zfft_avail) dir(3:3)='z'
         write(message,'("cfourier for pgrid [",a,"] available in directions [",a,"]")') trim(this%pg%name),dir; call log(message)
      end if
   end subroutine cfourier_log
   
   
   !> Print cfourier info to the screen
   subroutine cfourier_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(cfourier), intent(in) :: this
      character(len=3) :: dir
      if (this%pg%amRoot) then
         dir=''
         if (this%xfft_avail) dir(1:1)='x'
         if (this%yfft_avail) dir(2:2)='y'
         if (this%zfft_avail) dir(3:3)='z'
         write(output_unit,'("cfourier for pgrid [",a,"] available in directions [",a,"]")') trim(this%pg%name),dir
      end if
   end subroutine cfourier_print
   
   
end module cfourier_class