!> Volume fraction solver class:
!> Provides support for various BC, semi-Lagrangian geometric advancement,
!> curvature calculation, interface reconstruction.
module vfs_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: vfs,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   
   !> Boundary conditions for the volume fraction solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Volume fraction solver object definition
   type :: vfs
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_VFS'     !< Solver name (default=UNNAMED_VFS)
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! Volume fraction variable
      real(WP), dimension(:,:,:), allocatable :: VF       !< VF array
      
      ! Old volume fraction variable
      real(WP), dimension(:,:,:), allocatable :: VFold    !< VFold array
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for enforcing bconds
      
      ! Monitoring quantities
      real(WP) :: VFmax,VFmin,VFint                       !< Maximum, minimum, and integral volume fraction
      
   contains
      procedure :: print=>vfs_print                       !< Output solver to the screen
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: advance                                !< Advance VF to next step
      procedure :: get_max                                !< Calculate maximum field values
   end type vfs
   
   
   !> Declare volume fraction solver constructor
   interface vfs
      procedure constructor
   end interface vfs
   
contains
   
   
   !> Default constructor for volume fraction solver
   function constructor(cfg,name) result(self)
      use messager, only: die
      implicit none
      type(vfs) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate variables
      allocate(self%VF   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF   =0.0_WP
      allocate(self%VFold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFold=0.0_WP
      
      ! Prepare mask for VF
      allocate(self%mask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%mask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%mask(:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(self%cfg%imax+1:,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%mask(:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(:,self%cfg%jmax+1:,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%mask(:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(:,:,self%cfg%kmax+1:)=2
      end if
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%mask)
      
   end function constructor
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,dir,locator)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      character(len=2), intent(in) :: dir
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      select case (lowercase(dir))
      case ('+x','x+','xp','px'); new_bc%dir=1
      case ('-x','x-','xm','mx'); new_bc%dir=2
      case ('+y','y+','yp','py'); new_bc%dir=3
      case ('-y','y-','ym','my'); new_bc%dir=4
      case ('+z','z+','zp','pz'); new_bc%dir=5
      case ('-z','z-','zm','mz'); new_bc%dir=6
      case default; call die('[vfs add_bcond] Unknown bcond direction')
      end select
      new_bc%itr=iterator(this%cfg,new_bc%name,locator)
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         do n=1,new_bc%itr%n_
            i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
            this%mask(i,j,k)=2
         end do
      case (neumann)
         ! Not yet implemented
      case default
         call die('[vfs apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(vfs), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[vfs get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
            case (dirichlet)           ! Apply Dirichlet conditions
               
               ! This is done by the user directly
               ! Unclear whether we want to do this within the solver...
               
            case (neumann)             ! Apply Neumann condition
               
               ! Implement based on bcond direction
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%VF(i,j,k)=this%VF(i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[vfs apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%VF)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the new VF based on U/V/W
   subroutine advance(this,dt,U,V,W)
      implicit none
      class(vfs), intent(inout) :: this
      real(WP), intent(in) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Flux of VF
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face
      !         FX(i,j,k)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%itpsc_xp(:,i,j,k)*this%SC(i+this%stp1:i+this%stp2,j,k)) &
      !         &         -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%itpsc_xm(:,i,j,k)*this%SC(i+this%stm1:i+this%stm2,j,k)) &
      !         &         +this%diff*sum(this%grdsc_x(:,i,j,k)*this%SC(i-1:i,j,k))
               ! Fluxes on y-face
      !         FY(i,j,k)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%itpsc_yp(:,i,j,k)*this%SC(i,j+this%stp1:j+this%stp2,k)) &
      !         &         -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%itpsc_ym(:,i,j,k)*this%SC(i,j+this%stm1:j+this%stm2,k)) &
      !         &         +this%diff*sum(this%grdsc_y(:,i,j,k)*this%SC(i,j-1:j,k))
               ! Fluxes on z-face
      !         FZ(i,j,k)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%itpsc_zp(:,i,j,k)*this%SC(i,j,k+this%stp1:k+this%stp2)) &
      !         &         -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%itpsc_zm(:,i,j,k)*this%SC(i,j,k+this%stm1:k+this%stm2)) &
      !         &         +this%diff*sum(this%grdsc_z(:,i,j,k)*this%SC(i,j,k-1:k))
            end do
         end do
      end do
      ! Update VF
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
      !         drhoSCdt(i,j,k)=sum(this%divsc_x(:,i,j,k)*FX(i:i+1,j,k))+&
      !         &               sum(this%divsc_y(:,i,j,k)*FY(i,j:j+1,k))+&
      !         &               sum(this%divsc_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      ! Sync result
      call this%cfg%sync(this%VF)
   end subroutine advance
   
   
   !> Calculate the min/max/int of our VF field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(vfs), intent(inout) :: this
      integer :: ierr
      real(WP) :: my_VFmax,my_VFmin
      my_VFmax=maxval(this%VF); call MPI_ALLREDUCE(my_VFmax,this%VFmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_VFmin=minval(this%VF); call MPI_ALLREDUCE(my_VFmin,this%VFmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      call this%cfg%integrate(this%VF,integral=this%VFint)
   end subroutine get_max
   
   
   !> Print out info for vf solver
   subroutine vfs_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(vfs), intent(in) :: this
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Volume fraction solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
   end subroutine vfs_print


end module vfs_class
