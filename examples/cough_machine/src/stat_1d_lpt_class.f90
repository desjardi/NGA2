module stat_1d_lpt_class
   use precision,        only: WP
   use string,           only: str_medium,str_long
   use config_class,     only: config
   use lpt_class,        only: lpt
   use parallel,         only: comm
   use mpi_f08,          only: MPI_COMM

   implicit none
   private

   ! Expose type/constructor/methods
   public :: stat_1d_lpt

   ! A stat plane for lpt as a function of y/r for a x location
   type :: stat_1d_lpt
     ! Name of constructor
     character(len=str_medium) :: name
     ! x location to sample
     real(WP) :: xloc
     ! Number of stat variables
     integer :: nvar
     ! Number of radial bins
     integer :: nr
     ! Radial mesh spacing
     real(WP) :: dr
     ! Number of bins for diameter
     integer :: nd
     ! Drop diameter bin spacing
     real(WP) :: dd
     ! Drop bins
     real(WP), dimension(:), allocatable :: r
     ! Min and max diameter considered
     real(WP) :: dmin, dmax
     ! Tagged particle ID
     integer(kind=8) :: ID 
     ! Config file
     class(config), pointer  :: cfg
     ! Particle variable
     class(lpt), pointer :: lp
     ! Time current plane has had stats collected for
     real(WP) :: Delta_t
     ! Time current radial bin has collected for, 0 if no particle has passed
     real(WP), dimension(:,:), allocatable :: Deltat_r
     ! Number of particles that have passed through the radial bin
     real(WP), dimension(:,:), allocatable :: np_r
     ! Stat names
     character(len=str_medium), dimension(:), allocatable :: stat_name
     ! Stat variables and buffers
     real(WP), dimension(:,:,:), allocatable :: stat,buf
     ! Am I in?
     logical :: amIn
     ! Split Grid communicator
     type(MPI_Comm) :: comm
     ! Processor rank + num processors
     integer :: rank, nproc
     ! Am I the root processor?
     logical :: amRoot
  contains
    procedure :: sample => stat_1d_lpt_sample
    procedure :: write  => stat_1d_lpt_write
  end type stat_1d_lpt

  interface stat_1d_lpt
    procedure constructor
  end interface stat_1d_lpt

contains

function constructor(cfg,lp,xloc,name,dmin,dmax,dnbins,ID) result(self)
  use parallel, only : MPI_REAL_WP,comm
  use mpi_f08,  only : MPI_UNDEFINED,MPI_COMM_SPLIT
  use messager, only : die
  implicit none

  type(stat_1d_lpt) :: self
  class(config), target,  intent(in) :: cfg
  class(lpt),    target,  intent(in) :: lp
  real(WP),               intent(in) :: xloc
  character(len=*),       intent(in) :: name
  real(WP),               intent(in) :: dmin,dmax
  integer,                intent(in) :: dnbins
  integer(KIND=8),        intent(in) :: ID

  integer  :: color,key,ierr
  integer  :: j,jmid,k,kmid
  real(WP) :: dr

  ! Check if x is even in domain
  if (xloc.lt.cfg%x(cfg%imin) .or. xloc.gt.cfg%x(cfg%imax+1)) call die('[stat_1d_spray_class] error: xloc not inside domain')

  ! Name constructor and point to varibles
  self%name = adjustl(trim(name))
  self%xloc = xloc
  self%Delta_t = 0.0_WP
  self%dmin = dmin
  self%dmax = dmax
  self%nd = dnbins
  self%cfg => cfg
  self%lp  => lp
  self%ID  = ID

  ! Create directory to write to
  if (self%cfg%amRoot) call execute_command_line('mkdir -p stat')

  ! Name of stat variables
  self%nvar = 8
  allocate(self%stat_name(self%nvar))
  self%stat_name(1) = 'number'
  self%stat_name(2) = 'dd'
  self%stat_name(3) = 'dd^2'
  self%stat_name(4) = 'dd^3'
  self%stat_name(5) = 'ud'
  self%stat_name(6) = 'vd'
  self%stat_name(7) = 'wd'
  self%stat_name(8) = 'Td'

  ! Create sub-communicator
    ! Check if in local bounds
    if (self%cfg%x(self%cfg%imin_).le.self%xloc .and. self%xloc.lt.self%cfg%x(self%cfg%imax_+1)) then
      color = 0
      key = 0
      self%amIn = .true.
    else
    ! Out of bounds
      color = MPI_UNDEFINED
      key = 0
      self%amIn = .false.
    end if
    ! Create
    call MPI_COMM_SPLIT(comm,color,key,self%comm,ierr)
    ! Leave if not in
    if (.not.self%amIn) return
    ! Finish sizes and assign ranks
    call MPI_COMM_SIZE(self%comm,self%nproc,ierr)
    call MPI_COMM_RANK(self%comm,self%rank,ierr)

    ! Am I the root
    self%amRoot = (self%rank.eq.0)

    ! Basic checks
    !if (self%cfg%ny.eq.1.or.self%cfg%nz.eq.1) call die('[stat_1d_spray_class] error: Cannot be 2D in y or z')
    !if (self%cfg%ny.ne.self%cfg%nz)           call die('[stat_1d_spray_class] error: ny not equal to nz')

    ! Form radial bin
    self%nr = self%cfg%ny/2
    self%dr = self%cfg%yL/real(self%cfg%ny,WP)

    ! Form drop bin
    self%nd = dnbins
    self%dd = (dmax-dmin)/real(self%nd,WP)

    ! Form time radial bin
    allocate(self%Deltat_r(self%nd,self%nr)); self%Deltat_r = 0.0_WP

    ! Form number of particles radial bin
    allocate(self%np_r(self%nd,self%nr)); self%np_r = 0.0_WP

    ! Form buffers
    allocate(self%stat(self%nd,self%nr,self%nvar)); self%stat = 0.0_WP;
    allocate(self%buf (self%nd,self%nr,self%nvar)); self%buf  = 0.0_WP;

  return
end function constructor

subroutine stat_1d_lpt_sample(this,dt)
  use parallel, only : MPI_REAL_WP
  use mpi_f08,  only : MPI_SUM
  implicit none
  class(stat_1d_lpt), intent(inout) :: this
  real(WP), intent(in) :: dt

  integer   :: i,j,k,n,ierr
  real(WP)  :: rad, npart
  real(WP), dimension(:,:,:), allocatable :: buf, buf_sync

  ! Return if not collecting stats
  if (.not.this%amIn) return

  ! Increment time
  this%Delta_t = this%Delta_t + dt

  ! Zero out buffer
  allocate(buf     (this%nd,this%nr,this%nvar)); buf      = 0.0_WP
  allocate(buf_sync(this%nd,this%nr,this%nvar)); buf_sync = 0.0_WP

  ! Loop through all particles  in local processor
  do n = 1,this%lp%np_

    ! Return if not in particle bin
    if (this%lp%p(n)%d.lt.this%dmin .or. this%lp%p(n)%d.gt.this%dmax) cycle

    ! Return if not proper ID
    if (this%lp%p(n)%id.ne.this%ID) cycle 

    ! Test if the droplet just crossed the plane
    if ( this%lp%p(n)%pos(1).ge.this%xloc .and. this%lp%p(n)%pos(1)-dt*this%lp%p(n)%vel(1).lt.this%xloc .OR.   &
         this%lp%p(n)%pos(1).lt.this%xloc .and. this%lp%p(n)%pos(1)-dt*this%lp%p(n)%vel(1).ge.this%xloc ) then

      ! Check if it fits in a drop bin
      i = floor(  (this%lp%p(n)%d-this%dmin)  /this%dd ) + 1

      ! Check if it fits in radial bin
      rad = sqrt(this%lp%p(n)%pos(2)**2+this%lp%p(n)%pos(3)**2)
      j = floor(rad/this%dr) + 1

      ! Bin it
      if (i.le.this%nd.and.j.le.this%nr) then
        ! Increment a counter
        buf(i,j,1) = buf(i,j,1) + 1.0_WP
        ! Add diameters
        buf(i,j,2) = buf(i,j,2) + this%lp%p(n)%d**1
        buf(i,j,3) = buf(i,j,3) + this%lp%p(n)%d**2
        buf(i,j,4) = buf(i,j,4) + this%lp%p(n)%d**3
        ! Add velocities
        buf(i,j,5) = buf(i,j,5) + this%lp%p(n)%vel(1)
        buf(i,j,6) = buf(i,j,6) + this%lp%p(n)%vel(2)
        buf(i,j,7) = buf(i,j,7) + this%lp%p(n)%vel(3)
        ! Add particle time steps
        buf(i,j,8) = buf(i,j,8) + this%lp%p(n)%dt
      end if
    end if
  end do

  ! Integrate/sum over number of particles and synchronize over all processors
  call MPI_ALLREDUCE(buf,buf_sync,size(buf),MPI_REAL_WP,MPI_SUM,this%comm,ierr)

  ! Integrate over dt, particles is handled above
  !this%buf = this%buf + buf_sync*dt

  ! Increment time particles have been avg, 0 if no particles and the number of particles
  do j = 1,this%nr
    do i = 1,this%nd
      if (buf_sync(i,j,1).gt.0.0_WP) then
        this%Deltat_r(i,j) = this%Deltat_r(i,j) + dt
        this%np_r(i,j)     = this%np_r(i,j)     + buf_sync(i,j,1)
      end if
    end do
  end do

  ! Deallocate
  deallocate(buf,buf_sync)

  return
end subroutine stat_1d_lpt_sample

subroutine stat_1d_lpt_write(this)
  use parallel, only : MPI_REAL_WP, info_mpiio
  use mpi_f08
  class(stat_1d_lpt), intent(inout) :: this

  character(len=str_medium) :: filename
  integer :: i,j,ierr,var
  type(MPI_File) :: ifile
  type(MPI_Status) :: status

  ! Return if not in subcommunicator
  if (.not.this%amIn) return

  ! Get Number flux
  !this%stat(:,:,1) = this%buf(:,:,1) / this%Delta_t

  ! Normalize rest of stats by time and number of particles
  do j = 1,this%nr
    do i = 1,this%nd
      if(this%np_r(i,j).gt.0.0_WP) this%stat(i,j,2:) = this%buf(i,j,2:) / this%np_r(i,j)
    end do
  end do

  ! Open the file to write
  filename = "./stat/" // trim(this%name)
  call MPI_FILE_OPEN(this%comm,filename,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),info_mpiio,ifile,ierr)

  ! Only write if root
  if(this%amRoot) then

    ! Write the min and max particles size 
    call MPI_FILE_WRITE(ifile,this%dmin,1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_WRITE(ifile,this%dmax,1,MPI_REAL_WP,status,ierr)

    ! Write radial information first, with nd and dd, we can form droplet bin mesh completely
    call MPI_FILE_WRITE(ifile,this%nd,1,MPI_INTEGER,status,ierr)
    call MPI_FILE_WRITE(ifile,this%dd,1,MPI_REAL_WP,status,ierr)

    ! Write radial information first, with nr and dr, we can form radial bin mesh completely
    call MPI_FILE_WRITE(ifile,this%nr,1,MPI_INTEGER,status,ierr)
    call MPI_FILE_WRITE(ifile,this%dr,1,MPI_REAL_WP,status,ierr)

    ! Write out number of particle per radial bin
    call MPI_FILE_WRITE(ifile,this%np_r,this%nd*this%nr,MPI_REAL_WP,status,ierr)

    ! Write out time integrated per radial bin
    call MPI_FILE_WRITE(ifile,this%Deltat_r,this%nd*this%nr,MPI_REAL_WP,status,ierr)

    ! Write stat dimensions
    call MPI_FILE_WRITE(ifile,this%xloc,   1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_WRITE(ifile,this%Delta_t,1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_WRITE(ifile,this%nvar,   1,MPI_INTEGER,status,ierr)

    ! Write variable names
    do var=1,this%nvar
      call MPI_FILE_WRITE(ifile,this%stat_name(var),str_medium,MPI_CHARACTER,status,ierr)
    end do

    ! Write out stats
    do var=1,this%nvar
      call MPI_FILE_WRITE(ifile,this%stat(:,:,var),this%nd*this%nr,MPI_REAL_WP,status,ierr)
    end do

  end if
  ! Close file
  call MPI_FILE_CLOSE(ifile,ierr)

  return
end subroutine stat_1d_lpt_write


end module stat_1d_lpt_class
