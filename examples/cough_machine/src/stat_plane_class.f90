!> Class to record stats about lpt passing through a plane within the domain
module stat_plane_class
  use precision,        only: WP
  use string,           only: str_medium,str_long
  use config_class,     only: config
  use lpt_class,        only: lpt
  use parallel,         only: comm
  use mpi_f08,          only: MPI_COMM
  implicit none
  private

  ! Expose type/constructor/methods
  public :: stat_plane

  !> Definition of a stat plane for lpt 
  type :: stat_plane

    class(config), pointer :: cfg                                          !< Config for the domain
    
    ! Particle data
    real(WP)            :: dmin, dmax                                      !< Min and max diameter considered
    integer(kind=8)     :: ID                                              !< Tagged particle ID
    class(lpt), pointer :: lp                                              !< Particle variable
     
    ! Stat collection parameters 
    real(WP) :: dd                                                         !< Drop diameter bin spacing
    real(WP) :: Delta_t                                                    !< Time plane collected stats
    real(WP) :: location                                                   !< Location of plane along x,y or z axes
    integer  :: nvar                                                       !< Number of stat variables
    integer  :: nd                                                         !< Number of bins for diameter groupings
    real(WP), dimension(:,:), allocatable :: stat,buf                      !< Stat variables and buffers
    character(len=str_medium)                            :: plane_pos      !< Axes that the plane lies on
    character(len=str_medium)                            :: name           !< Name of constructor
    character(len=str_medium), dimension(:), allocatable :: stat_name      !< Name of recorded stats 

    ! Processors and communicators
    integer :: rank, nproc                                                 !< Processor rank + num processors
    logical :: amIn                                                        !< Check to ensure loop is in local bounds
    logical :: amRoot                                                      !< Declare root processor
    type(MPI_Comm) :: comm                                                 !< Split Grid communicator

  contains
    procedure :: sample => stat_lpt_sample                                 !< Sample and bin particles in the domain
    procedure :: write  => stat_lpt_write                                  !< Write particle stats to files
  end type stat_plane

  !> Declare lpt_stat constrcuctor 
  interface stat_plane
    procedure constructor
  end interface stat_plane

contains

!> Constructor for stat plane
function constructor(cfg,lp,location,plane_pos,name,dmin,dmax,dnbins,ID) result(self)
  use parallel, only : MPI_REAL_WP,comm
  use mpi_f08,  only : MPI_UNDEFINED,MPI_COMM_SPLIT
  use messager, only : die
  implicit none
  type(stat_plane_lpt) :: self
  integer  :: color,key,ierr
  integer  :: j,jmid,k,kmid
  ! User defined input parameters
  class(config), target, intent(in) :: cfg
  class(lpt),    target, intent(in) :: lp
  real(WP),              intent(in) :: location
  character(len=*),      intent(in) :: plane_pos
  real(WP),              intent(in) :: dmin,dmax
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: dnbins
  integer(KIND=8),       intent(in) :: ID

  ! Name constructor and point to varibles
  self%name       = adjustl(trim(name))
  self%location   = location
  self%plane_pos  = adjustl(trim(plane_pos))
  self%Delta_t    = 0.0_WP
  self%dmin       = dmin
  self%dmax       = dmax
  self%nd         = dnbins
  self%cfg        => cfg
  self%lp         => lp
  self%ID         = ID

  ! Create directory to write to
  if (self%cfg%amRoot) call execute_command_line('mkdir -p stat2')

  ! Name of stat variables
  self%nvar = 11
  allocate(self%stat_name(self%nvar))
  self%stat_name(1)  = 'number'
  self%stat_name(2)  = 'dd'
  self%stat_name(3)  = 'dd^2'
  self%stat_name(4)  = 'dd^3'
  self%stat_name(5)  = 'ud'
  self%stat_name(6)  = 'vd'
  self%stat_name(7)  = 'wd'
  self%stat_name(8)  = 'Td'
  self%stat_name(9)  = 'xpos'
  self%stat_name(10) = 'ypos'
  self%stat_name(11) = 'zpos'

  ! Create sub-communicator
    ! Check local bounds for given case
    select case(plane_pos)
    case ('x')
      if (self%cfg%x(self%cfg%imin_).le.self%location.and.self%location.lt.self%cfg%x(self%cfg%imax_+1)) then
        ! In bounds
        color=0.0_WP
        key=0.0_WP
        self%amIn=.true.
      else
      ! Out of bounds
        color=MPI_UNDEFINED
        key=0.0_WP
        self%amIn=.false.
      end if
    case ('y')
      if (self%cfg%y(self%cfg%jmin_).le.self%location.and.self%location.lt.self%cfg%y(self%cfg%jmax_+1)) then
        ! In bounds
        color=0.0_WP
        key=0.0_WP
        self%amIn=.true.
      else
      ! Out of bounds
        color=MPI_UNDEFINED
        key=0.0_WP
        self%amIn=.false.
      end if
    case ('z')
      if (self%cfg%z(self%cfg%kmin_).le.self%location.and.self%location.lt.self%cfg%z(self%cfg%kmax_+1)) then
        ! In bounds
        color=0.0_WP
        key=0.0_WP
        self%amIn=.true.
      else
      ! Out of bounds
        color=MPI_UNDEFINED
        key=0.0_WP
        self%amIn=.false.
      end if
    end select

    ! Create
    call MPI_COMM_SPLIT(comm,color,key,self%comm,ierr)
    
    ! Leave if not in
    if (.not.self%amIn) return
    
    ! Finish sizes and assign ranks
    call MPI_COMM_SIZE(self%comm,self%nproc,ierr)
    call MPI_COMM_RANK(self%comm,self%rank,ierr)

    ! Am I the root
    self%amRoot = (self%rank.eq.0)

    ! Form drop bin
    self%nd = dnbins
    self%dd = (dmax-dmin)/real(self%nd,WP)

    ! Form buffers
    allocate(self%stat(self%nd,self%nvar)); self%stat=0.0_WP;
    allocate(self%buf (self%nd,self%nvar)); self%buf =0.0_WP;

  return
end function constructor

subroutine stat_lpt_sample(this,dt)
  use parallel, only : MPI_REAL_WP
  use mpi_f08,  only : MPI_SUM
  implicit none
  class(stat_plane_lpt), intent(inout) :: this
  real(WP), intent(in) :: dt
  real(WP), dimension(:,:), allocatable :: buf,buf_sync
  integer :: i,j,k,n,ierr

  ! Return if not collecting stats
  if (.not.this%amIn) return

  ! Increment time
  this%Delta_t=this%Delta_t+dt

  ! Zero out buffer
  allocate(buf     (this%nd,this%nvar)); buf     =0.0_WP
  allocate(buf_sync(this%nd,this%nvar)); buf_sync=0.0_WP

  ! Loop through all particles  in local processor
  do n = 1,this%lp%np_

    ! Return if not in particle bin
    if (this%lp%p(n)%d.lt.this%dmin.or.this%lp%p(n)%d.gt.this%dmax) cycle

    ! Return if not proper ID
    if (this%lp%p(n)%id.ne.this%ID) cycle 

    ! Test if the droplet just crossed the plane and bin it
    select case(this%plane_pos)
    case ('x')
      if (this%lp%p(n)%pos(1).ge.this%location.and.this%lp%p(n)%pos(1)-dt*this%lp%p(n)%vel(1).lt.this%location.or. &
          this%lp%p(n)%pos(1).lt.this%location.and.this%lp%p(n)%pos(1)-dt*this%lp%p(n)%vel(1).ge.this%location) then

        ! Check if it fits in a drop bin
        i = floor((this%lp%p(n)%d-this%dmin)/this%dd)+1

        ! Bin it
        if (i.le.this%nd) then
          ! Increment a counter
          buf(i,1)=buf(i,1)+1.0_WP
          ! Add diameters
          buf(i,2)=buf(i,2)+this%lp%p(n)%d**1
          buf(i,3)=buf(i,3)+this%lp%p(n)%d**2
          buf(i,4)=buf(i,4)+this%lp%p(n)%d**3
          ! Add velocities
          buf(i,5)=buf(i,5)+this%lp%p(n)%vel(1)
          buf(i,6)=buf(i,6)+this%lp%p(n)%vel(2)
          buf(i,7)=buf(i,7)+this%lp%p(n)%vel(3)
          ! Add particle time steps
          buf(i,8)=buf(i,8)+this%lp%p(n)%dt
          ! Add posistions
          buf(i,9)=buf(i,9)+this%lp%p(n)%pos(1)
          buf(i,10)=buf(i,10)+this%lp%p(n)%pos(2)
          buf(i,11)=buf(i,11)+this%lp%p(n)%pos(3)
        end if
        
      end if
    case ('y')
      if (this%lp%p(n)%pos(2).ge.this%location.and.this%lp%p(n)%pos(2)-dt*this%lp%p(n)%vel(2).lt.this%location.or. &
          this%lp%p(n)%pos(2).lt.this%location.and.this%lp%p(n)%pos(2)-dt*this%lp%p(n)%vel(2).ge.this%location) then

        ! Check if it fits in a drop bin
        i = floor((this%lp%p(n)%d-this%dmin)/this%dd)+1

        ! Bin it
        if (i.le.this%nd) then
          ! Increment a counter
          buf(i,1)=buf(i,1)+1.0_WP
          ! Add diameters
          buf(i,2)=buf(i,2)+this%lp%p(n)%d**1
          buf(i,3)=buf(i,3)+this%lp%p(n)%d**2
          buf(i,4)=buf(i,4)+this%lp%p(n)%d**3
          ! Add velocities
          buf(i,5)=buf(i,5)+this%lp%p(n)%vel(1)
          buf(i,6)=buf(i,6)+this%lp%p(n)%vel(2)
          buf(i,7)=buf(i,7)+this%lp%p(n)%vel(3)
          ! Add particle time steps
          buf(i,8)=buf(i,8)+this%lp%p(n)%dt
          ! Add posistions
          buf(i,9)=buf(i,9)+this%lp%p(n)%pos(1)
          buf(i,10)=buf(i,10)+this%lp%p(n)%pos(2)
          buf(i,11)=buf(i,11)+this%lp%p(n)%pos(3)
        end if
        
      end if
    case ('z')
      if (this%lp%p(n)%pos(3).ge.this%location.and.this%lp%p(n)%pos(3)-dt*this%lp%p(n)%vel(3).lt.this%location.or. &
          this%lp%p(n)%pos(3).lt.this%location.and.this%lp%p(n)%pos(3)-dt*this%lp%p(n)%vel(3).ge.this%location) then

        ! Check if it fits in a drop bin
        i = floor((this%lp%p(n)%d-this%dmin)/this%dd)+1

        ! Bin it
        if (i.le.this%nd) then
          ! Increment a counter
          buf(i,1)=buf(i,1)+1.0_WP
          ! Add diameters
          buf(i,2)=buf(i,2)+this%lp%p(n)%d**1
          buf(i,3)=buf(i,3)+this%lp%p(n)%d**2
          buf(i,4)=buf(i,4)+this%lp%p(n)%d**3
          ! Add velocities
          buf(i,5)=buf(i,5)+this%lp%p(n)%vel(1)
          buf(i,6)=buf(i,6)+this%lp%p(n)%vel(2)
          buf(i,7)=buf(i,7)+this%lp%p(n)%vel(3)
          ! Add particle time steps
          buf(i,8)=buf(i,8)+this%lp%p(n)%dt
          ! Add posistions
          buf(i,9)=buf(i,9)+this%lp%p(n)%pos(1)
          buf(i,10)=buf(i,10)+this%lp%p(n)%pos(2)
          buf(i,11)=buf(i,11)+this%lp%p(n)%pos(3)
        end if
        
      end if
    end select


  end do

  ! Integrate/sum over number of particles and synchronize over all processors
  call MPI_ALLREDUCE(buf,buf_sync,size(buf),MPI_REAL_WP,MPI_SUM,this%comm,ierr)

  ! Deallocate
  deallocate(buf,buf_sync)

  return
end subroutine stat_lpt_sample

subroutine stat_lpt_write(this)
  use parallel, only : MPI_REAL_WP, info_mpiio
  use mpi_f08
  class(stat_plane_lpt), intent(inout) :: this

  character(len=str_medium) :: filename
  integer :: i,j,ierr,var
  type(MPI_File) :: ifile
  type(MPI_Status) :: status

  ! Return if not in subcommunicator
  if (.not.this%amIn) return

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

    ! Write stat dimensions
    call MPI_FILE_WRITE(ifile,this%location,1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_WRITE(ifile,this%Delta_t, 1,MPI_REAL_WP,status,ierr)
    call MPI_FILE_WRITE(ifile,this%nvar,    1,MPI_INTEGER,status,ierr)

    ! Write variable names
    do var=1,this%nvar
      call MPI_FILE_WRITE(ifile,this%stat_name(var),str_medium,MPI_CHARACTER,status,ierr)
    end do

    ! Write out stats
    do var=1,this%nvar
      call MPI_FILE_WRITE(ifile,this%stat(:,var),this%nd,MPI_REAL_WP,status,ierr)
    end do

  end if
  ! Close file
  call MPI_FILE_CLOSE(ifile,ierr)

  return
end subroutine stat_lpt_write


end module stat_plane_class
