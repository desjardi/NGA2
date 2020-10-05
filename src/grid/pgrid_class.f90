!> Definition of a partitioned grid class in NGA.
!> @todo Add other parallelization strategies
module pgrid_class
   use precision,   only: WP
   use sgrid_class, only: sgrid
   use string,      only: str_medium
   implicit none
   private
   
   !> Default parallelization strategy
   character(len=str_medium), parameter :: defstrat='fewest_dir'
   integer, parameter :: defmincell=4
   
   ! Expose type/constructor/methods
   public :: pgrid
   
   !> Partitioned grid type
   type, extends(sgrid) :: pgrid
      ! Parallelization information
      integer :: group                      !< Grid group
      integer :: comm                       !< Grid communicator
      integer :: nproc                      !< Number of processors
      integer :: rank                       !< Processor grid rank
      logical :: amRoot                     !< Am I grid root?
      logical :: amIn                       !< Am I working in this grid?
      integer :: npx,npy,npz                !< Number of processors per direction
      integer :: iproc,jproc,kproc          !< Coordinates location of processor
      ! Local grid size
      integer :: nx_,ny_,nz_                !< Local grid size in x/y/z
      ! Local grid size with overlap
      integer :: nxo_,nyo_,nzo_             !< Local grid size in x/y/z with overlap
      ! Local index bounds
      integer :: imin_,imax_                !< Domain-decomposed index bounds in x
      integer :: jmin_,jmax_                !< Domain-decomposed index bounds in y
      integer :: kmin_,kmax_                !< Domain-decomposed index bounds in z
      ! Local index bounds with overlap
      integer :: imino_,imaxo_              !< Domain-decomposed index bounds in x with overlap
      integer :: jmino_,jmaxo_              !< Domain-decomposed index bounds in y with overlap
      integer :: kmino_,kmaxo_              !< Domain-decomposed index bounds in z with overlap
   contains
      procedure, private :: init_mpi=>pgrid_init_mpi           !< Prepare the MPI environment for the pgrid
      procedure, private :: domain_decomp=>pgrid_domain_decomp !< Perform the domain decomposition
      procedure :: allprint=>pgrid_allprint                    !< Output grid to screen - blocking and requires all procs...
      procedure :: print=>pgrid_print                          !< Output grid to screen
      procedure :: log=>pgrid_log                          !< Output grid info to log
   end type pgrid
   
   
   !> Declare partitioned grid constructor
   interface pgrid
      procedure construct_pgrid_from_sgrid
      procedure construct_pgrid_from_file
   end interface pgrid
   
contains
   
   
   !> Partitioned grid constructor from file
   function construct_pgrid_from_file(no,file,grp,decomp) result(self)
      use string,   only: lowercase
      use monitor,  only: die
      use param,    only: verbose
      use parallel, only: MPI_REAL_WP
      use mpi
      implicit none
      
      type(pgrid) :: self                               !< Parallel grid
      integer, intent(in) :: no                         !< Overlap size
      character(len=*), intent(in) :: file              !< Grid file
      integer, intent(in) :: grp                        !< Processor group for parallelization
      integer, dimension(3), intent(in) :: decomp       !< Desired decomposition
      integer :: ierr
      character(len=str_medium) :: simu_name
      integer :: nx,ny,nz,coord
      real(WP), dimension(:), allocatable :: x
      real(WP), dimension(:), allocatable :: y
      real(WP), dimension(:), allocatable :: z
      logical :: xper,yper,zper
      
      ! Initialize MPI environment
      self%group=grp; call self%init_mpi
      
      ! Nothing more to do if the processor is not inside
      if (.not.self%amIn) return
      
      ! Root process can now read in the grid
      if (self%amRoot) then
         self%sgrid=sgrid(no,file)
         simu_name=self%name
         coord=self%coordsys
         xper=self%xper
         yper=self%yper
         zper=self%zper
         nx=self%nx
         ny=self%ny
         nz=self%nz
      end if
      
      ! Broadcast it to our group
      call MPI_BCAST(simu_name,len(simu_name),MPI_CHARACTER,0,self%comm,ierr)
      call MPI_BCAST(coord    ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      call MPI_BCAST(xper     ,1             ,MPI_LOGICAL  ,0,self%comm,ierr)
      call MPI_BCAST(yper     ,1             ,MPI_LOGICAL  ,0,self%comm,ierr)
      call MPI_BCAST(zper     ,1             ,MPI_LOGICAL  ,0,self%comm,ierr)
      call MPI_BCAST(nx       ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      call MPI_BCAST(ny       ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      call MPI_BCAST(nz       ,1             ,MPI_INTEGER  ,0,self%comm,ierr)
      
      ! Allocate x/y/z, fill it, and bc
      allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1))
      if (self%amRoot) then
         x(1:nx+1)=self%x(self%imin:self%imax+1)
         y(1:ny+1)=self%y(self%jmin:self%jmax+1)
         z(1:nz+1)=self%z(self%kmin:self%kmax+1)
      end if
      call MPI_BCAST(x,nx+1,MPI_REAL_WP,0,self%comm,ierr)
      call MPI_BCAST(y,ny+1,MPI_REAL_WP,0,self%comm,ierr)
      call MPI_BCAST(z,nz+1,MPI_REAL_WP,0,self%comm,ierr)
      
      ! Finish creating the sgrid
      if (.not.self%amRoot) self%sgrid=sgrid(coord,no,x,y,z,xper,yper,zper,trim(adjustl(simu_name)))
      
      ! Deallocate
      deallocate(x,y,z)
      
      ! Store decomposition on the grid
      self%npx=decomp(1); self%npy=decomp(2); self%npz=decomp(3)
      
      ! Perform actual domain decomposition of grid
      call self%domain_decomp()
      
      ! If verbose run, log and or print grid
      if (verbose.gt.0) call self%log
      if (verbose.gt.1) call self%print
      if (verbose.gt.2) call self%allprint
      
   end function construct_pgrid_from_file
   
   
   !> Partitioned grid constructor from sgrid
   function construct_pgrid_from_sgrid(grid,grp,strat,decomp) result(self)
      use string,  only: lowercase
      use monitor, only: die
      use param,   only: verbose
      implicit none
      include 'mpif.h'
      
      type(pgrid) :: self                               !< Parallel grid
      
      type(sgrid), intent(in) :: grid                   !< Base grid
      integer, intent(in) :: grp                        !< Processor group for parallelization
      character(len=str_medium), optional :: strat      !< Requested parallelization strategy
      integer, dimension(3), optional :: decomp         !< Requested domain decomposition
      
      integer, dimension(3) :: mydecomp
      character(len=str_medium) :: mystrat
      
      ! Initialize MPI environment
      self%group=grp; call self%init_mpi
      
      ! Nothing more to do if the processor is not inside
      if (.not.self%amIn) return
      
      ! Assign base grid data
      self%sgrid=grid
      
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
      
      ! Store decomposition on the grid
      self%npx=mydecomp(1); self%npy=mydecomp(2); self%npz=mydecomp(3)
      
      ! Perform actual domain decomposition of grid
      call self%domain_decomp()
      
      ! If verbose run, log and or print grid
      if (verbose.gt.0) call self%log
      if (verbose.gt.1) call self%print
      if (verbose.gt.2) call self%allprint
      
   contains
      
      !> Decomposition calculator: fewest direction strategy
      function fewest_dir_decomp(nc,nd) result(dec)
         use monitor, only: warn
         integer, dimension(3) :: dec             !< Returned decomposition
         integer, dimension(3), intent(in) :: nc  !< Grid size in 3 directions
         integer, intent(in) :: nd                !< Number of domains
         integer :: i1,dir1,max1,q1,r1
         integer :: i2,dir2,max2,q2,r2
         integer :: i3,dir3,max3
         logical, dimension(3) :: mask
         ! Order grid directions by size
         mask( 1:3)=.true. ; dir3=maxloc(nc,1,mask=mask)
         mask(dir3)=.false.; dir2=maxloc(nc,1,mask=mask)
         mask(dir2)=.false.; dir1=maxloc(nc,1,mask=mask)
         ! Maximum number of domains in each direction
         max1=nc(dir1)/defmincell
         max2=nc(dir2)/defmincell
         max3=nc(dir3)/defmincell
         ! Perform trial division in first direction
         loop1: do i1=1,min(nd,max1)
            ! Perform Euclidian division
            q1=nd/i1; r1=mod(nd,i1)
            ! Cycle if not an integer factor
            if (r1.ne.0) cycle loop1
            ! Perform trial division in second direction
            loop2: do i2=1,min(q1,max2)
               ! Perform Euclidian division
               q2=q1/i2; r2=mod(q1,i2)
               ! Cycle if not an integer factor
               if (r2.ne.0) cycle loop2
               ! Set decomposition in third direction
               i3=q2
               ! If it is valid, return
               if (q2.le.max3) then
                  dec(dir1)=i1
                  dec(dir2)=i2
                  dec(dir3)=i3
                  return
               end if
            end do loop2
         end do loop1
         ! If we are still here, we have explored all options and failed
         ! Return zero and a warning, as we may want to try another strategy
         dec=0; call warn('[fewest_dir_decomp] Grid parallelization failed...')
      end function fewest_dir_decomp
      
   end function construct_pgrid_from_sgrid
   
   
   !> Prepares the MPI environment for the pgrid
   subroutine pgrid_init_mpi(self)
      use parallel, only: comm
      use monitor , only: die
      implicit none
      include 'mpif.h'
      class(pgrid) :: self
      integer :: ierr
      ! Get group size, and intracommunicator
      call MPI_GROUP_SIZE(self%group,self%nproc,ierr)
      if (self%nproc.eq.0) call die('[pgrid constructor] Non-empty group is required!')
      call MPI_COMM_CREATE(comm,self%group,self%comm,ierr)
      ! Get rank and amIn info
      call MPI_GROUP_RANK(self%group,self%rank,ierr)
      self%amIn=(self%rank.ne.MPI_UNDEFINED)
      ! Handle processors that are not part of the group
      if (.not.self%amIn) then
         self%amRoot=.false.
         self%iproc=0; self%nx_=0; self%imin_=0; self%imax_=0; self%nxo_=0; self%imino_=0; self%imaxo_=0
         self%jproc=0; self%ny_=0; self%jmin_=0; self%jmax_=0; self%nyo_=0; self%jmino_=0; self%jmaxo_=0
         self%kproc=0; self%nz_=0; self%kmin_=0; self%kmax_=0; self%nzo_=0; self%kmino_=0; self%kmaxo_=0
      end if
   end subroutine pgrid_init_mpi
   
   
   !> Prepares the domain decomposition of the pgrid
   subroutine pgrid_domain_decomp(self)
      use monitor , only: die
      implicit none
      include 'mpif.h'
      class(pgrid) :: self
      integer :: ierr,q,r
      integer, parameter :: ndims=3
      logical, parameter :: reorder=.true.
      integer, dimension(3) :: coords
      
      ! Perform sanity check of the decomposition
      if (self%npx*self%npy*self%npz.ne.self%nproc) call die('[pgrid constructor] Parallel decomposition is improper')
      
      ! Give cartesian layout to intracommunicator
      call MPI_CART_CREATE(self%comm,ndims,[self%npx,self%npy,self%npz],[self%xper,self%yper,self%zper],reorder,self%comm,ierr)
      call MPI_COMM_RANK  (self%comm,self%rank,ierr)
      call MPI_CART_COORDS(self%comm,self%rank,ndims,coords,ierr)
      self%iproc=coords(1)+1; self%jproc=coords(2)+1; self%kproc=coords(3)+1
      self%amRoot=(self%rank.eq.0)
      
      ! Perform decomposition in x
      q=self%nx/self%npx; r=mod(self%nx,self%npx)
      if (self%iproc.le.r) then
         self%nx_  =q+1
         self%imin_=self%imin+(self%iproc-1)*self%nx_
      else
         self%nx_  =q
         self%imin_=self%imin+(self%iproc-1)*self%nx_+r
      end if
      self%imax_ =self%imin_+self%nx_-1
      self%nxo_  =self%nx_+2*self%no
      self%imino_=self%imin_-self%no
      self%imaxo_=self%imax_+self%no
      
      ! Perform decomposition in y
      q=self%ny/self%npy; r=mod(self%ny,self%npy)
      if (self%jproc.le.r) then
         self%ny_  =q+1
         self%jmin_=self%jmin+(self%jproc-1)*self%ny_
      else
         self%ny_  =q
         self%jmin_=self%jmin+(self%jproc-1)*self%ny_+r
      end if
      self%jmax_ =self%jmin_+self%ny_-1
      self%nyo_  =self%ny_+2*self%no
      self%jmino_=self%jmin_-self%no
      self%jmaxo_=self%jmax_+self%no
      
      ! Perform decomposition in z
      q=self%nz/self%npz; r=mod(self%nz,self%npz)
      if (self%kproc.le.r) then
         self%nz_  =q+1
         self%kmin_=self%kmin+(self%kproc-1)*self%nz_
      else
         self%nz_  =q
         self%kmin_=self%kmin+(self%kproc-1)*self%nz_+r
      end if
      self%kmax_ =self%kmin_+self%nz_-1
      self%nzo_  =self%nz_+2*self%no
      self%kmino_=self%kmin_-self%no
      self%kmaxo_=self%kmax_+self%no
      
   end subroutine pgrid_domain_decomp
   
   
   !> Print out partitioned grid info to screen
   !> This is a slow and blocking routine for debugging only!
   subroutine pgrid_allprint(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use parallel, only: rank
      implicit none
      class(pgrid), intent(in) :: this
      integer :: i,ierr
      ! Return all non-involved procs
      if (.not.this%amIn) return
      ! Root writes header
      call pgrid_print(this)
      ! Each proc provides info on his involvement in the grid
      do i=0,this%nproc-1
         ! Block for clean output
         call MPI_BARRIER(this%comm,ierr)
         ! Output info
         if (this%rank.eq.i) then
            write(output_unit,'(" >>> Rank ",i0,"(",i0,") -> [",i0,",",i0,",",i0,"] owns [",i0,",",i0,"]x[",i0,",",i0,"]x[",i0,",",i0,"]")') &
            this%rank,rank,this%iproc,this%jproc,this%kproc,this%imin_,this%imax_,this%jmin_,this%jmax_,this%kmin_,this%kmax_
         end if
      end do
   end subroutine pgrid_allprint
   
   
   !> Cheap print of partitioned grid info to screen
   subroutine pgrid_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use sgrid_class, only: cartesian,cylindrical,spherical
      implicit none
      class(pgrid), intent(in) :: this
      if (this%amRoot) then
         select case (this%coordsys)
         case (cartesian);   write(output_unit,'("Partitioned cartesian grid [",a,"] on ",i0," processes")') trim(this%name),this%nproc
         case (cylindrical); write(output_unit,'("Partitioned cylindrical grid [",a,"]")') trim(this%name)
         case (spherical);   write(output_unit,'("Partitioned spherical grid [",a,"]")') trim(this%name)
         end select
         write(output_unit,'(" >   overlap = ",i0)') this%no
         write(output_unit,'(" >      size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz
         write(output_unit,'(" >    decomp = ",i0,"x",i0,"x",i0)') this%npx,this%npy,this%npz
         write(output_unit,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1)
         write(output_unit,'(" >   uniform = ",l1,"x",l1,"x",l1)') this%uniform_x,this%uniform_y,this%uniform_z
         write(output_unit,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper
         write(output_unit,'(" >    volume = ",es12.5)') this%vol_total
      end if
   end subroutine pgrid_print
   
   
   !> Cheap print of partitioned grid info to log
   subroutine pgrid_log(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use monitor,     only: log
      use string,      only: str_long
      use sgrid_class, only: cartesian,cylindrical,spherical
      implicit none
      class(pgrid), intent(in) :: this
      character(len=str_long) :: message
      if (this%amRoot) then
         select case (this%coordsys)
         case (cartesian);   write(message,'("Partitioned cartesian grid [",a,"] on ",i0," processes")') trim(this%name),this%nproc; call log(message)
         case (cylindrical); write(message,'("Partitioned cylindrical grid [",a,"]")') trim(this%name); call log(message)
         case (spherical);   write(message,'("Partitioned spherical grid [",a,"]")') trim(this%name); call log(message)
         end select
         write(message,'(" >   overlap = ",i0)') this%no; call log(message)
         write(message,'(" >      size = ",i0,"x",i0,"x",i0)') this%nx,this%ny,this%nz; call log(message)
         write(message,'(" >    decomp = ",i0,"x",i0,"x",i0)') this%npx,this%npy,this%npz; call log(message)
         write(message,'(" >    extent = [",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]x[",es12.5,",",es12.5,"]")') this%x(this%imin),this%x(this%imax+1),this%y(this%jmin),this%y(this%jmax+1),this%z(this%kmin),this%z(this%kmax+1); call log(message)
         write(message,'(" >   uniform = ",l1,"x",l1,"x",l1)') this%uniform_x,this%uniform_y,this%uniform_z; call log(message)
         write(message,'(" >  periodic = ",l1,"x",l1,"x",l1)') this%xper,this%yper,this%zper; call log(message)
         write(message,'(" >    volume = ",es12.5)') this%vol_total; call log(message)
      end if
   end subroutine pgrid_log
   
   
end module pgrid_class
