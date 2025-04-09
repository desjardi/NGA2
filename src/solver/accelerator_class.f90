!> Non-linear convergence accelerator
module accelerator_class
   use precision,    only: WP
   use config_class, only: config
   use string,       only: str_medium
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: accelerator
   
   !> accelerator object definition
   type :: accelerator
      type(config), pointer :: cfg                                      !< Config for the accelerator
      character(len=str_medium) :: name='UNNAMED_ACCELERATOR'           !< Name of accelerator
      real(WP), dimension(:,:,:,:), allocatable :: G,Q                  !< Data storage of the form m_max*nx*ny*nz
      real(WP), dimension(:,:)    , allocatable :: R                    !< Data storage of the form m_max*m_max
      real(WP), dimension(:,:,:)  , allocatable :: x,fcur,fold,gold,tmp !< Additional storage of intermediate fields (nx*ny*nz)
      integer :: k                                                      !< Current iteration number (first call is with k=0)
      integer :: m                                                      !< Number of fields currently stored
      integer :: m_max                                                  !< Maximum number of fields that can be stored
      integer :: k_start=0                                              !< Starting iteration (default is 0)
      real(WP) :: beta=1.0_WP                                           !< Damping parameter (should be 0<beta<=1, default is 1)
      real(WP) :: drop_tol=1.0e-5_WP                                    !< Tolerance below which data is dropped for better conditioning
      real(WP), dimension(:), allocatable :: A,b,work                   !< Storage for linear solvers
      integer , dimension(:), allocatable :: iwork                      !< Storage for linear solvers
   contains
      procedure :: initialize                                           !< Initialization of the accelerator
      procedure :: increment                                            !< Accelerated increment
      procedure :: restart                                              !< Restart accelerator
   end type accelerator
   
contains
   
   !> In-place initialization for accelerator object
   subroutine initialize(this,cfg,storage_size,name)
      use messager, only: die
      implicit none
      class(accelerator), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: storage_size
      character(len=*), optional :: name
      
      ! Set the name for the accelerator
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Link the config object
      this%cfg=>cfg
      
      ! Storage size
      if (storage_size.lt.1) call die('[Accelerator initialize] storage_size must be greater or equal to 1')
      this%m_max=storage_size
      
      ! Allocate storage
      allocate(this%R(1:this%m_max,1:this%m_max)); this%R=0.0_WP
      allocate(this%G(1:this%m_max,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%G=0.0_WP
      allocate(this%Q(1:this%m_max,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%Q=0.0_WP
      allocate(this%x             (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%x=0.0_WP
      allocate(this%fcur          (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%fcur=0.0_WP
      allocate(this%fold          (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%fold=0.0_WP
      allocate(this%gold          (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%gold=0.0_WP
      allocate(this%tmp           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%tmp=0.0_WP
      allocate(this%work(3*this%m_max),this%iwork(this%m_max),this%A(1:this%m_max*(this%m_max+1)/2),this%b(1:this%m_max))

      ! Set iteration counter to zero by default
      this%k=0
      
   end subroutine initialize
   
   !> Restart accelerator
   subroutine restart(this,ig)
      implicit none
      class(accelerator), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: ig !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      ! Reset iteration counter to zero
      this%k=0
      ! Reset number of fields stored to zero
      this%m=0
      ! Zero out G
      this%G=0.0_WP
      ! Store initial guess
      this%x=ig
   end subroutine restart
   
   !> Increment accelerator with new function evaluation provided in gcur
   subroutine increment(this,gcur)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      class(accelerator), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: gcur !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(:,:,:), allocatable :: df
      real(WP) :: tmp
      integer :: n,i,j,k,ierr,info
      
      ! Apply damping on new evaluation
      gcur=this%beta*gcur+(1.0_WP-this%beta)*this%x
      
      ! Evaluate fcur
      this%fcur=gcur-this%x
      
      ! Increment solution storage if k>k_start
      if (this%k>this%k_start) then
         ! Allocate and compute df
         allocate(df(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); df=this%fcur-this%fold
         ! If we have run out of storage, make room
         if (this%m.ge.this%m_max) then
            do n=1,this%m_max-1
               this%G(n,:,:,:)=this%G(n+1,:,:,:)
            end do
         end if
         ! Store new delta g and increment counter
         this%G(min(this%m+1,this%m_max),:,:,:)=gcur-this%gold
         this%m=this%m+1
      end if
      
      ! Increment iteration counter
      this%k=this%k+1
      
      ! Remember fields
      this%fold=this%fcur
      this%gold=     gcur
      
      ! Handle no data case
      if (this%m.eq.0) then
         this%x=gcur
         return
      end if
      
      ! Still here, so we have data for acceleration
      if (this%m.eq.1) then
         ! Calculate norm(df)
         tmp=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            tmp=tmp+df(i,j,k)**2
         end do; end do; end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); tmp=sqrt(tmp)
         ! Add new entry to Q and R
         this%Q(1,:,:,:)=df/tmp; this%R(1,1)=tmp
      else
         ! Delete data from QR decomposition if needed
         if (this%m.gt.this%m_max) then
            call qr_delete(); this%m=this%m-1
         end if
         ! Recalculate R and df
         do n=1,this%m-1
            ! Calculate Q.df
            tmp=0.0_WP
            do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
               tmp=tmp+this%Q(n,i,j,k)*df(i,j,k)
            end do; end do; end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
            ! Update R and df
            this%R(n,this%m)=tmp
            df=df-this%R(n,this%m)*this%Q(n,:,:,:)
         end do
         ! Calculate norm(df)
         tmp=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            tmp=tmp+df(i,j,k)**2
         end do; end do; end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); tmp=sqrt(tmp)
         ! Add new entry to Q and R
         this%Q(this%m,:,:,:)=df/tmp; this%R(this%m,this%m)=tmp
      end if
      
      ! Test conditioning of problem and form A for inversion
      do j=1,this%m; do i=1,this%m
         this%A(i+(j-1)*j/2)=this%R(i,j)
      end do; end do
      call dtpcon('1','U','N',this%m,this%A,tmp,this%work,this%iwork,info)
      do while (tmp.lt.this%drop_tol.and.this%m.gt.1)
         ! Drop data
         call qr_delete(); this%m=this%m-1
         ! Recalculate condition number
         do j=1,this%m; do i=1,this%m
            this%A(i+(j-1)*j/2)=this%R(i,j)
         end do; end do
         call dtpcon('1','U','N',this%m,this%A,tmp,this%work,this%iwork,info)
      end do
      
      ! Form b from Q.fcur
      do n=1,this%m
         ! Calculate Q.fcur
         tmp=0.0_WP
         do k=this%cfg%kmin_,this%cfg%kmax_; do j=this%cfg%jmin_,this%cfg%jmax_; do i=this%cfg%imin_,this%cfg%imax_
            tmp=tmp+this%Q(n,i,j,k)*this%fcur(i,j,k)
         end do; end do; end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! Store in b
         this%b(n)=tmp
      end do
      
      ! Solve with lapack
      call dtptrs('U','N','N',this%m,1,this%A,this%b,this%m,info)
      
      ! Update solution
      this%x=gcur
      do n=1,this%m
         this%x=this%x-this%b(n)*this%G(n,:,:,:)
      end do
      
      ! Standard damping method for Anderson acceleration is implemented below
      ! Tests suggest it's not as good as simply under-relaxing each new eval
      !this%x=this%x-(1.0_WP-this%beta)*this%fcur
      !do n=1,this%m
      !   tmp=sum(this%R(n,1:this%m)*this%b(1:this%m))
      !   this%x=this%x+(1.0_WP-this%beta)*this%Q(n,:,:,:)*tmp
      !end do
      
      ! Deallocate df
      deallocate(df)
      
      ! Finally, return solution in gcur for convenience
      gcur=this%x
      
   contains
      
      ! Delete data from QR
      subroutine qr_delete()
         implicit none
         integer :: ii,jj
         real(WP) :: temp,c,s
         do ii=1,this%m_max-1
            temp=sqrt(this%R(ii,ii+1)**2+this%R(ii+1,ii+1)**2)
            c=this%R(ii  ,ii+1)/temp
            s=this%R(ii+1,ii+1)/temp
            this%R(ii  ,ii+1)=temp
            this%R(ii+1,ii+1)=0.0_WP
            if (ii.lt.this%m_max-1) then
               do jj=ii+2,this%m_max
                  temp=c*this%R(ii,jj)+s*this%R(ii+1,jj)
                  this%R(ii+1,jj)=-s*this%R(ii,jj)+c*this%R(ii+1,jj)
                  this%R(ii,jj)=temp
               end do
            end if
            this%tmp=c*this%Q(ii,:,:,:)+s*this%Q(ii+1,:,:,:)
            this%Q(ii+1,:,:,:)=-s*this%Q(ii,:,:,:)+c*this%Q(ii+1,:,:,:)
            this%Q(ii  ,:,:,:)=this%tmp
         end do
         this%Q(this%m_max,:,:,:)=0.0_WP
         this%R(this%m_max,:)=0.0_WP
         do jj=1,this%m_max-1
            this%R(:,jj)=this%R(:,jj+1)
         end do
         this%R(:,this%m_max)=0.0_WP
      end subroutine qr_delete
      
   end subroutine increment
   
end module accelerator_class
