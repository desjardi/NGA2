!> Definition of a parallel grid class in NGA.
!> @todo Optimize the information stored in the derived data type...
module pgrid_class
  use bgrid_class, only: bgrid
  use string, only: str_medium
  implicit none
  private

  !> Default parallelization strategy
  character(len=str_medium), parameter :: defstrat='fewest_dir'
  integer, parameter :: defmincell=4
  
  ! Expose type/constructor/methods
  public :: pgrid
  
  !> Parallel grid type
  type, extends(bgrid) :: pgrid
     ! Parallelization information
     integer :: group                      !< Grid group
     integer :: comm                       !< Grid communicator
     integer :: nproc                      !< Number of processors
     integer :: rank                       !< Processor grid rank
     logical :: amroot                     !< Am I grid root?
     logical :: amin                       !< Am I working in this grid?
     integer :: npx,npy,npz                !< Number of processors per direction
     integer :: iproc,jproc,kproc          !< Coordinates location of processor
     ! Local grid size
     integer :: nx_,ny_,nz_                !< Local grid size in x/y/z
     ! Local index bounds
     integer :: imin_,imax_                !< Domain-decomposed index bounds in x
     integer :: jmin_,jmax_                !< Domain-decomposed index bounds in y
     integer :: kmin_,kmax_                !< Domain-decomposed index bounds in z
   contains
     procedure :: allprint=>pgrid_allprint !< Output grid to screen - blocking and requires all procs...
     procedure :: print=>pgrid_print       !< Output grid to screen
  end type pgrid
  
  !> Declare parallel grid constructor
  interface pgrid
     procedure constructor
  end interface pgrid
    
contains
  
  !> Parallel grid constructor
  function constructor(grid,grp,strat,decomp) result(self)
    use parallel
    use string, only: lowercase
    use monitor, only: die
    implicit none
    include 'mpif.h'
    
    type(pgrid) :: self                               !< Parallel grid
    type(bgrid), intent(in) :: grid                   !< Base grid
    integer, intent(in) :: grp                        !< Processor group for parallelization
    character(len=str_medium), optional :: strat      !< Requested parallelization strategy
    integer, dimension(3), optional :: decomp         !< Requested domain decomposition

    integer, dimension(3) :: mydecomp
    character(len=str_medium) :: mystrat
    integer :: ierr,q,r
    integer, parameter :: ndims=3
    logical, parameter :: reorder=.true.
    logical, dimension(3) :: isper
    integer, dimension(3) :: coords
    
    ! Assign base grid data
    self%bgrid=grid
    
    ! Get group name, group size, and intracommunicator
    self%group=grp
    call MPI_GROUP_SIZE(self%group,self%nproc,ierr) 
    call MPI_COMM_CREATE(comm,self%group,self%comm,ierr)
    
    ! Get rank and amin info
    call MPI_GROUP_RANK(self%group,self%rank,ierr)
    self%amin=(self%rank.ne.MPI_UNDEFINED)
    
    ! Check if a decomposition was provided
    if (present(decomp)) then
       ! Check if a decomposition strategy was provided
       if (present(strat)) then
          mystrat=trim(lowercase(strat))
          if (mystrat.ne.'imposed') call die('[pgrid constructor] If decomp is provided, strat=imposed is expected')
       else
          mystrat='imposed'
       end if
       ! Set the decomposition
       mydecomp=decomp
    else
       ! Check if a decomposition strategy was provided
       if (present(strat)) then 
          mystrat=trim(lowercase(strat))
          if (mystrat.eq.'imposed') call die('[pgrid constructor] If strat=imposed, decomp is expected')
       else
          mystrat=defstrat
       end if
       ! Compute decomposition
       select case (mystrat)
       case ('imposed');    ! Handled elsewhere
       case ('fewest_dir'); mydecomp=fewest_dir_decomp([self%nx,self%ny,self%nz],self%nproc)
       case default;        call die('[pgrid constructor] Unknown parallel decomposition strategy: '//trim(mystrat))
       end select
    end if
    
    ! Perform sanity check of the decomposition
    if (mydecomp(1)*mydecomp(2)*mydecomp(3).ne.self%nproc) call die('[pgrid constructor] Parallel decomposition is improper!')

    ! Stpre decomposition on the grid
    self%npx=mydecomp(1); self%npy=mydecomp(2); self%npz=mydecomp(3)
    
    ! Handle processors that are not part of the group
    if (.not.self%amin) then
       self%amroot=.false.
       self%iproc=0; self%nx_=0; self%imin_=0; self%imax_=0
       self%jproc=0; self%ny_=0; self%jmin_=0; self%jmax_=0
       self%kproc=0; self%nz_=0; self%kmin_=0; self%kmax_=0
       return
    end if
    
    ! Give cartesian layout to intracommunicator
    isper=[.false.,.false.,.false.]
    call MPI_CART_CREATE(self%comm,ndims,mydecomp,isper,reorder,self%comm,ierr)
    call MPI_COMM_RANK(self%comm,self%rank,ierr)
    call MPI_CART_COORDS(self%comm,self%rank,ndims,coords,ierr)
    self%iproc=coords(1)+1
    self%jproc=coords(2)+1
    self%kproc=coords(3)+1
    self%amroot=(self%rank.eq.0)
    
    ! Perform decomposition in x
    q=self%nx/self%npx; r=mod(self%nx,self%npx)
    if (self%iproc.le.r) then
       self%nx_  =q+1
       self%imin_=self%imin+(self%iproc-1)*self%nx_
    else
       self%nx_  =q
       self%imin_=self%imin+(self%iproc-1)*self%nx_+r
    end if
    self%imax_=self%imin_+self%nx_-1
    
    ! Perform decomposition in y
    q=self%ny/self%npy; r=mod(self%ny,self%npy)
    if (self%jproc.le.r) then
       self%ny_  =q+1
       self%jmin_=self%jmin+(self%jproc-1)*self%ny_
    else
       self%ny_  =q
       self%jmin_=self%jmin+(self%jproc-1)*self%ny_+r
    end if
    self%jmax_=self%jmin_+self%ny_-1

    ! Perform decomposition in z
    q=self%nz/self%npz; r=mod(self%nz,self%npz)
    if (self%kproc.le.r) then
       self%nz_  =q+1
       self%kmin_=self%kmin+(self%kproc-1)*self%nz_
    else
       self%nz_  =q
       self%kmin_=self%kmin+(self%kproc-1)*self%nz_+r
    end if
    self%kmax_=self%kmin_+self%nz_-1
    
  contains
    
    !> Decomposition calculator: fewest direction strategy
    function fewest_dir_decomp(nc,nd) result(dec)
      integer, dimension(3) :: dec             !< Returned decomposition
      integer, dimension(3), intent(in) :: nc  !< Grid size in 3 directions
      integer, intent(in) :: nd                !< Number of domains
      integer :: lmin,lmed,lmax,q,r,nd_max
      logical, dimension(3) :: mask
      ! Order grid size
      mask( 1:3)=.true. ; lmax=maxloc(nc,1,mask=mask)
      mask(lmax)=.false.; lmed=maxloc(nc,1,mask=mask)
      mask(lmed)=.false.; lmin=maxloc(nc,1,mask=mask)
      ! Maximum number of domains in lmax direction
      nd_max=nc(lmax)/defmincell
      ! Eucledian division
      q=nd/nd_max; r=mod(nd,nd_max)
      if (q.eq.0) then
         ! We're done
         dec=1; dec(lmax)=r; return
      else
         ! Need zero remainder
         do while (r.ne.0)
            !q=q-1; r=nd-q*
         end do
         dec(lmax)=nd_max
         
      end if
      
    end function fewest_dir_decomp
    
  end function constructor

  !> Print out parallel grid info to screen
  !> This is a slow and blocking routine for debugging only!
  subroutine pgrid_allprint(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    use parallel
    implicit none
    class(pgrid), intent(in) :: this
    integer :: i,ierr
    ! Return all non-involved procs
    if (.not.this%amin) return
    ! Root writes header
    call pgrid_print(this)
    ! Each proc provides info on his involvement in the grid
    do i=0,this%nproc-1
       ! Block for clean output
       call MPI_BARRIER(this%comm,ierr)
       ! Output info
       if (this%rank.eq.i) then
          write(output_unit,'(" >>> Rank ",i0,"(",i0,") -> [",i0,",",i0,",",i0,"] owns [",i0,",",i0,"]x[",i0,",",i0,"]x[",i0,",",i0,"]")') this%rank,rank,this%iproc,this%jproc,this%kproc,this%imin_,this%imax_,this%jmin_,this%jmax_,this%kmin_,this%kmax_
       end if
    end do
  end subroutine pgrid_allprint

  !> Cheap print of parallel grid info to screen
  subroutine pgrid_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    use parallel, only: 
    implicit none
    class(pgrid), intent(in) :: this
    if (this%amroot) then
       write(output_unit,'("Parallel grid ",a," on ",i0," processes")') trim(this%name),this%nproc
       write(output_unit,'(" >   size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
       write(output_unit,'(" > decomp = ",i0,"x",i0,"x",i0)') this%npx,this%npy,this%npz
       write(output_unit,'(" > extent = [",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]x[",es12.6,",",es12.6,"]")') &
            this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
    end if
  end subroutine pgrid_print
  
end module pgrid_class
