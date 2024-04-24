!> Parallel data file concept is defined here: given a partitioned grid and
!> an I/O partition, it provides efficient parallel I/O access to a data file
module pardata_class
   use precision,      only: WP
   use string,         only: str_short,str_medium
   use pgrid_class,    only: pgrid
   use coupler_class,  only: coupler
   use datafile_class, only: datafile
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: pardata   
   

   !> pardata object definition as an array of values and variables
   !> The choice is made here to *copy* the data to simplify the process
   type :: pardata
      ! pardata has two partitioned grid
      type(pgrid), pointer     :: pg                                  !< Original partitioned grid
      type(pgrid), allocatable :: pg_io                               !< Partitioned grid for I/O
      ! I/O partition info
      logical :: ingrp_io
      integer, dimension(3) :: partition_io
      ! Couplers for transfering data for reading and writing
      type(coupler) :: wcpl,rcpl
      ! Filename for read/write operations
      character(len=str_medium) :: filename
      ! A pardata stores scalar values
      integer :: nval                                                 !< Number of scalar values
      character(len=str_short), dimension(:), allocatable :: valname  !< Name of scalar values
      real(WP), dimension(:), allocatable :: val                      !< Values
      ! A pardata stores real(WP) variables of pgrid-compatible size
      integer :: nvar                                                 !< Number of field variables
      character(len=str_short), dimension(:), allocatable :: varname  !< Name of field variables
      real(WP), dimension(:,:,:,:), allocatable :: var                !< Variables (this is with overlap!)
      ! A pardata uses a datafile object on the I/O partition
      type(datafile) :: df                                            !< Datafile object on pg_io for actual i/o
   contains
      generic :: initialize=>init_from_args,init_from_file            !< Generic initialization for pardata
      procedure, private :: init_from_args                            !< Initialize an empty pardata object from arguments
      procedure, private :: init_from_file                            !< Initialize a pardata object from an existing file
      procedure, private :: prep_iomap                                !< IO group/pgrid/map
      procedure, private :: findval                                   !< Function that returns val index if name is found, zero otherwise
      procedure, private :: findvar                                   !< Function that returns var index if name is found, zero otherwise
      generic :: push=>pushval,pushvar                                !< Generic data push
      procedure, private :: pushval                                   !< Push data to pardata
      procedure, private :: pushvar                                   !< Push data to pardata
      generic :: pull=>pullval,pullvar                                !< Generic data pull
      procedure, private :: pullval                                   !< Pull data from pardata
      procedure, private :: pullvar                                   !< Pull data from pardata
      procedure :: write=>pardata_write                               !< Parallel write a pardata object to disk
      procedure :: print=>pardata_print                               !< Print out debugging info to screen
      procedure :: log  =>pardata_log                                 !< Print out debugging info to log
      final     :: destructor                                         !< Destructor for pardata
   end type pardata
   
   
contains
   
   
   !> Destructor for pardata object
   subroutine destructor(this)
      implicit none
      type(pardata) :: this
      if (allocated(this%valname)) deallocate(this%valname)
      if (allocated(this%val    )) deallocate(this%val    )
      if (allocated(this%varname)) deallocate(this%varname)
      if (allocated(this%var    )) deallocate(this%var    )
   end subroutine destructor
   
   
   !> Prepare the io group and its pgrid, along with the couplers
   subroutine prep_iomap(this,iopartition)
      use messager, only: die
      use mpi_f08,  only: MPI_Group_range_incl,MPI_Group
      implicit none
      class(pardata) :: this
      integer, dimension(3), intent(in) :: iopartition
      type(MPI_Group) :: iogrp
      integer, dimension(3,1) :: iorange
      integer :: ierr
      
      ! Store IO partition
      this%partition_io=iopartition
      
      ! Basic checks on the partition
      if (product(this%partition_io).le.0) call die('[pardata constructor] At least one process must be assigned to IO')
      if (product(this%partition_io).gt.this%pg%nproc) call die('[pardata constructor] Number of IO processes must be less or equal to pgrid processes')
      
      ! Create IO MPI group using the first ranks (could explore other strategies here)
      iorange(:,1)=[0,product(this%partition_io)-1,1]
      call MPI_Group_range_incl(this%pg%group,1,iorange,iogrp,ierr)
      
      ! Create IO logical
      this%ingrp_io=.false.
      if (this%pg%rank.le.product(this%partition_io)-1) this%ingrp_io=.true.
      
      ! Create IO pgrid
      if (this%ingrp_io) this%pg_io=pgrid(this%pg%sgrid,iogrp,this%partition_io)
      
      ! Prepare coupler for writing
      this%wcpl=coupler(src_grp=this%pg%group,dst_grp=iogrp,name='write')
      call this%wcpl%set_src(this%pg); if (this%ingrp_io) call this%wcpl%set_dst(this%pg_io)
      call this%wcpl%initialize()
      
      ! Prepare coupler for reading
      this%rcpl=coupler(src_grp=iogrp,dst_grp=this%pg%group,name='read')
      if (this%ingrp_io) call this%rcpl%set_src(this%pg_io); call this%rcpl%set_dst(this%pg)
      call this%rcpl%initialize()
      
   end subroutine prep_iomap
   
   
   !> Initialization for an empty pardata object using user-provided arguments
   subroutine init_from_args(this,pg,iopartition,filename,nval,nvar)
      implicit none
      class(pardata), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      integer, dimension(3), intent(in) :: iopartition
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nval,nvar
      
      ! Link the partitioned grid and store the filename
      this%pg=>pg
      this%filename=trim(adjustl(filename))
      
      ! Allocate storage and store names for vals
      this%nval=nval
      allocate(this%valname(this%nval))
      this%valname=''
      allocate(this%val(this%nval))
      this%val=0.0_WP
      
      ! Allocate storage and store names for vars
      this%nvar=nvar
      allocate(this%varname(this%nvar))
      this%varname=''
      allocate(this%var(this%pg%imino_:this%pg%imaxo_,this%pg%jmino_:this%pg%jmaxo_,this%pg%kmino_:this%pg%kmaxo_,this%nvar))
      this%var=0.0_WP
      
      ! Prepare the IO-specialized partition and couplers
      call this%prep_iomap(iopartition)
      
      ! IO processes create a datafile object
      if (this%ingrp_io) then
         this%df=datafile(this%pg_io,filename,nval,nvar)
         deallocate(this%df%var) !< Reduce memory footprint
      end if
      
   end subroutine init_from_args


   !> Initialization for a pardata object based on an existing file
   subroutine init_from_file(this,pg,iopartition,fdata)
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_BARRIER,MPI_BCAST,MPI_INTEGER,MPI_CHARACTER
      implicit none
      class(pardata), intent(inout) :: this
      class(pgrid), target, intent(in) :: pg
      integer, dimension(3), intent(in) :: iopartition
      character(len=*), intent(in) :: fdata
      real(WP), dimension(:,:,:), allocatable :: temp
      integer :: ierr,n
      
      ! Link the partitioned grid and store the filename
      this%pg=>pg
      this%filename=trim(adjustl(fdata))
      
      ! Prepare the IO-specialized partition and couplers
      call this%prep_iomap(iopartition)
      
      ! IO processes create a datafile object from the file
      if (this%ingrp_io) then
         this%df=datafile(this%pg_io,fdata)
         allocate(temp(this%df%pg%imino_:this%df%pg%imaxo_,this%df%pg%jmino_:this%df%pg%jmaxo_,this%df%pg%kmino_:this%df%pg%kmaxo_))
      end if
      
      ! Root on pgio (i.e., source root on rcpl) copies this%df serial data into pardata object and broadcasts
      ! Transfer nval
      if (this%rcpl%amRoot) this%nval   =this%df%nval   ; call MPI_BCAST(this%nval,1,MPI_INTEGER,this%rcpl%sroot,this%rcpl%comm,ierr)
      ! Transfer valname
      allocate(this%valname(this%nval))
      if (this%rcpl%amRoot) this%valname=this%df%valname; call MPI_BCAST(this%valname,str_short*this%nval,MPI_CHARACTER,this%rcpl%sroot,this%rcpl%comm,ierr)
      ! Transfer val
      allocate(this%val(this%nval))
      if (this%rcpl%amRoot) this%val    =this%df%val    ; call MPI_BCAST(this%val,this%nval,MPI_REAL_WP,this%rcpl%sroot,this%rcpl%comm,ierr)
      ! Transfer nvar
      if (this%rcpl%amRoot) this%nvar   =this%df%nvar   ; call MPI_BCAST(this%nvar,1,MPI_INTEGER,this%rcpl%sroot,this%rcpl%comm,ierr)
      ! Transfer varname
      allocate(this%varname(this%nvar))
      if (this%rcpl%amRoot) this%varname=this%df%varname; call MPI_BCAST(this%varname,str_short*this%nvar,MPI_CHARACTER,this%rcpl%sroot,this%rcpl%comm,ierr)
      
      ! Transfer parallel data using the rcpl
      allocate(this%var(this%pg%imino_:this%pg%imaxo_,this%pg%jmino_:this%pg%jmaxo_,this%pg%kmino_:this%pg%kmaxo_,this%nvar))
      do n=1,this%nvar
         if (this%ingrp_io) then
            temp(this%df%pg%imin_:this%df%pg%imax_,this%df%pg%jmin_:this%df%pg%jmax_,this%df%pg%kmin_:this%df%pg%kmax_)=this%df%var(:,:,:,n)
            call this%df%pg%sync(temp)
            call this%rcpl%push(temp)
         end if
         call this%rcpl%transfer()
         call this%rcpl%pull(this%var(:,:,:,n))
         call MPI_BARRIER(this%pg%comm,ierr)
      end do
      
      ! Reduce memory footprint
      if (this%ingrp_io) deallocate(this%df%var,temp)
      
   end subroutine init_from_file
   
   
   !> Write a pardata object to a file
   subroutine pardata_write(this,fdata)
      use mpi_f08, only: MPI_BARRIER
      implicit none
      class(pardata), intent(inout) :: this
      character(len=*), optional :: fdata
      real(WP), dimension(:,:,:), allocatable :: temp
      integer :: n,ierr
      
      ! Transfer serial data to this%df
      if (this%ingrp_io) then
         this%df%valname=this%valname
         this%df%val    =this%val
         this%df%varname=this%varname
         if (.not.allocated(this%df%var)) allocate(this%df%var(this%df%pg%imin_:this%df%pg%imax_,this%df%pg%jmin_:this%df%pg%jmax_,this%df%pg%kmin_:this%df%pg%kmax_,this%df%nvar))
         allocate(temp(this%df%pg%imino_:this%df%pg%imaxo_,this%df%pg%jmino_:this%df%pg%jmaxo_,this%df%pg%kmino_:this%df%pg%kmaxo_))
         this%df%var=0.0_WP
      end if
      
      ! Transfer parallel data using the wcpl
      do n=1,this%nvar
         call this%wcpl%push(this%var(:,:,:,n))
         call this%wcpl%transfer()
         if (this%ingrp_io) then
            call this%wcpl%pull(temp)
            this%df%var(:,:,:,n)=temp(this%df%pg%imin_:this%df%pg%imax_,this%df%pg%jmin_:this%df%pg%jmax_,this%df%pg%kmin_:this%df%pg%kmax_)
         end if
         call MPI_BARRIER(this%pg%comm,ierr)
      end do
      
      ! Now perform write operation on pgio
      if (this%ingrp_io) then
         if (present(fdata)) then
            call this%df%write(fdata)
         else
            call this%df%write()
         end if
         ! Reduce memory footprint
         deallocate(this%df%var,temp)
      end if
      
   end subroutine pardata_write
   
   
   !> Print pardata content to the screen
   subroutine pardata_print(this,text)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(pardata), intent(in) :: this
      character(*), intent(in) :: text
      if (this%pg%amRoot) then
         write(output_unit,'(a," pardata [",a,"] for partitioned grid [",a,"]")') trim(text),trim(this%filename),trim(this%pg%name)
         write(output_unit,'(" >      nval = ",i0)') this%nval
         write(output_unit,'(" > val names = ",1000(a,1x))') this%valname
         write(output_unit,'(" >      nvar = ",i0)') this%nvar
         write(output_unit,'(" > var names = ",1000(a,1x))') this%varname
         write(output_unit,'(" > IO partition is ",i0,"x",i0,"x",i0)') this%pg_io%npx,this%pg_io%npy,this%pg_io%npz
      end if
   end subroutine pardata_print
   
   
   !> Print pardata content to the log
   subroutine pardata_log(this,text)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(pardata), intent(in) :: this
      character(*), intent(in) :: text
      character(len=str_long) :: message
      if (this%pg%amRoot) then
         write(message,'(a," pardata [",a,"] for partitioned grid [",a,"]")') trim(text),trim(this%filename),trim(this%pg%name); call log(message)
         write(message,'(" >      nval = ",i0)') this%nval; call log(message)
         write(message,'(" > val names = ",1000(a,1x))') this%valname; call log(message)
         write(message,'(" >      nvar = ",i0)') this%nvar; call log(message)
         write(message,'(" > var names = ",1000(a,1x))') this%varname; call log(message)
         write(message,'(" > IO partition is ",i0,"x",i0,"x",i0)') this%pg_io%npx,this%pg_io%npy,this%pg_io%npz; call log(message)
      end if
   end subroutine pardata_log
   
   
   !> Index finding for val
   function findval(this,name) result(ind)
      implicit none
      class(pardata), intent(in) :: this
      integer :: ind
      character(len=*), intent(in) :: name
      integer :: n
      ind=0
      do n=1,this%nval
         if (this%valname(n).eq.name) ind=n
      end do
   end function findval
   
   
   !> Push data to a val
   subroutine pushval(this,name,val)
      use messager, only: die
      implicit none
      class(pardata), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(in) :: val
      integer :: n
      n=this%findval(name)
      if (n.gt.0) then
         this%val(n)=val
      else
         call die('[pardata pushval] Val does not exist in the data file: '//name)
      end if
   end subroutine pushval
   
   
   !> Pull data from a val
   subroutine pullval(this,name,val)
      use messager, only: die
      implicit none
      class(pardata), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), intent(out) :: val
      integer :: n
      n=this%findval(name)
      if (n.gt.0) then
         val=this%val(n)
      else
         call die('[pardata pullval] Val does not exist in the data file: '//name)
      end if
   end subroutine pullval
   
   
   !> Index finding for var
   function findvar(this,name) result(ind)
      implicit none
      class(pardata), intent(in) :: this
      integer :: ind
      character(len=*), intent(in) :: name
      integer :: n
      ind=0
      do n=1,this%nvar
         if (this%varname(n).eq.name) ind=n
      end do
   end function findvar
   
   
   !> Push data to a var
   subroutine pushvar(this,name,var)
      use messager, only: die
      implicit none
      class(pardata), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), intent(in) :: var !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n
      n=this%findvar(name)
      if (n.gt.0) then
         this%var(:,:,:,n)=var
      else
         call die('[pardata pushvar] Var does not exist in the data file: '//name)
      end if
   end subroutine pushvar
   
   
   !> Pull data from a var and synchronize it
   subroutine pullvar(this,name,var)
      use messager, only: die
      implicit none
      class(pardata), intent(in) :: this
      character(len=*), intent(in) :: name
      real(WP), dimension(this%pg%imino_:,this%pg%jmino_:,this%pg%kmino_:), intent(out) :: var !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n
      n=this%findvar(name)
      if (n.gt.0) then
         var=this%var(:,:,:,n)
         call this%pg%sync(var)
      else
         call die('[pardata pullvar] Var does not exist in the data file: '//name)
      end if
   end subroutine pullvar
   
   
end module pardata_class
