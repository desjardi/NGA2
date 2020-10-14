!> Ensight class concept is defined here: given a config object,
!> it provides parallel I/O access to an ensight file
module ensight_class
   use precision,    only: WP
   use string,       only: str_short,str_medium
   use config_class, only: config
   use mpi_f08,      only: MPI_Datatype
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: ensight
   
   ! List types
   type :: scl
      type(scl), pointer :: next
      character(len=str_short) :: name
      real(WP), dimension(:,:,:), pointer :: scl_ptr
   end type scl
   type :: vct
      type(vct), pointer :: next
      character(len=str_short) :: name
      real(WP), dimension(:,:,:), pointer :: vctx_ptr
      real(WP), dimension(:,:,:), pointer :: vcty_ptr
      real(WP), dimension(:,:,:), pointer :: vctz_ptr
   end type vct
   
   !> Ensight object definition as list of pointers to arrays
   type :: ensight
      ! An ensight object has a case file
      character(len=str_medium) :: casename                           !< Name of casefile to read/write
      ! An ensight object stores time values
      integer :: ntime                                                !< Number of scalar values
      real(WP), dimension(:), allocatable :: time                     !< Time values
      ! An ensight object stores geometry data
      type(config), pointer :: cfg                                    !< Config for ensight geometry and parallel I/O
      ! An ensight object stores lists of pointers to data
      type(scl), pointer :: first_scl                                 !< Scalar list
      type(vct), pointer :: first_vct                                 !< Vector list
   contains
      procedure :: write_geom                                         !< Write out geometry
   end type ensight
   
   
   !> Declare ensight constructor
   interface ensight
      procedure construct_ensight
   end interface ensight
   
   
contains
   
   !> Constructor for an empty ensight object
   function construct_ensight(cfg,casename) result(self)
      use monitor, only: die
      implicit none
      type(ensight) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: casename
      integer :: iunit,ierr
      
      ! Link to config
      self%cfg=>cfg
      
      ! Initially, no time info
      self%ntime=0
      
      ! Store casename
      self%casename=trim(adjustl(casename))
      
      ! Create directory
      if (self%cfg%amRoot) call execute_command_line('mkdir -p ensight')
      
      ! Write out the geometry
      call self%write_geom()
      
      ! Empty pointer lists for now
      self%first_scl=>NULL()
      self%first_vct=>NULL()
      
      ! Output a barebone case file
      open(newunit=iunit,file='ensight/'//trim(adjustl(casename))//'.case',form='formatted',status='replace',access='stream',iostat=ierr)
      write(iunit,'(4(a,/))') 'FORMAT','type: ensight gold','GEOMETRY','model: geometry'
      close(iunit)
      
   end function construct_ensight
   
   
   !> Geometry output to a file in parallel
   subroutine write_geom(this)
      use precision, only: SP
      use monitor,   only: die
      use parallel,  only: info_mpiio
      use mpi_f08
      implicit none
      class(ensight) :: this
      integer :: iunit,ierr
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff
      type(MPI_Datatype) :: view
      type(MPI_File) :: ifile
      type(MPI_Status):: status
      integer(kind=MPI_OFFSET_KIND) :: disp
      integer, dimension(:,:,:), allocatable :: iblank
      integer :: i,i1,i2,j,j1,j2,k,k1,k2,datasize
      integer, dimension(3) :: gsizes,lsizes,lstart
      
      ! Root does most of the I/O at first
      if (this%cfg%amRoot) then
         
         ! First create a new file
         open(newunit=iunit,file='ensight/geometry',form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write geom] Could not open file: ensight/geometry')
         
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'        ; write(iunit) cbuff
         cbuff=trim(adjustl(this%cfg%name))        ; write(iunit) cbuff
         cbuff='node id off'                       ; write(iunit) cbuff
         cbuff='element id off'                    ; write(iunit) cbuff
         
         ! Extents
         cbuff='extents'                           ; write(iunit) cbuff
         rbuff=real(this%cfg%x(this%cfg%imin  ),SP); write(iunit) rbuff
         rbuff=real(this%cfg%x(this%cfg%imax+1),SP); write(iunit) rbuff
         rbuff=real(this%cfg%y(this%cfg%jmin  ),SP); write(iunit) rbuff
         rbuff=real(this%cfg%y(this%cfg%jmax+1),SP); write(iunit) rbuff
         rbuff=real(this%cfg%z(this%cfg%kmin  ),SP); write(iunit) rbuff
         rbuff=real(this%cfg%z(this%cfg%kmax+1),SP); write(iunit) rbuff
         
         ! Part header
         cbuff='part'                              ; write(iunit) cbuff
         ibuff=1                                   ; write(iunit) ibuff
         cbuff='Complete geometry'                 ; write(iunit) cbuff  ! We have a single grid-cfg here
         cbuff='block rectilinear iblanked'        ; write(iunit) cbuff
         
         ! Number of cells
         ibuff=this%cfg%nx+1                       ; write(iunit) ibuff
         ibuff=this%cfg%ny+1                       ; write(iunit) ibuff
         ibuff=this%cfg%nz+1                       ; write(iunit) ibuff
         
         ! Mesh
         write(iunit) real(this%cfg%x(this%cfg%imin:this%cfg%imax+1),SP)
         write(iunit) real(this%cfg%y(this%cfg%jmin:this%cfg%jmax+1),SP)
         write(iunit) real(this%cfg%z(this%cfg%kmin:this%cfg%kmax+1),SP)
         
         ! Close the file
         close(iunit)
         
      end if
      
      ! Set array size
      i1=this%cfg%imin_; i2=this%cfg%imax_; if (this%cfg%iproc.eq.this%cfg%npx) i2=i2+1
      j1=this%cfg%jmin_; j2=this%cfg%jmax_; if (this%cfg%jproc.eq.this%cfg%npy) j2=j2+1
      k1=this%cfg%kmin_; k2=this%cfg%kmax_; if (this%cfg%kproc.eq.this%cfg%npz) k2=k2+1
      datasize=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
      
      ! Allocate iblank array
      allocate(iblank(i1:i2,j1:j2,k1:k2))
      
      ! Build the blanking info and x/y/zbuf
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               ! Blanking from mask
               if (int(minval(this%cfg%mask(i-1:i,j-1:j,k-1:k))).eq.0) then
                  iblank(i,j,k)=1
               else
                  iblank(i,j,k)=0
               end if
            end do
         end do
      end do
      
      ! We need to define a proper MPI-I/O view
      gsizes=[this%cfg%nx+1,this%cfg%ny+1,this%cfg%nz+1]
      lsizes=[i2-i1+1,j2-j1+1,k2-k1+1]
      lstart=[this%cfg%imin_-this%cfg%imin,this%cfg%jmin_-this%cfg%jmin,this%cfg%kmin_-this%cfg%kmin]
      call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_INTEGER,view,ierr)
      call MPI_TYPE_COMMIT(view,ierr)
      
      ! Only need to parallel write masks now
      call MPI_FILE_OPEN(this%cfg%comm,'ensight/geometry',IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[ensight write geom] Problem encountered while parallel writing geometry file')
      call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_INTEGER,view,'native',info_mpiio,ierr)
      call MPI_FILE_WRITE_ALL(ifile,iblank,datasize,MPI_INTEGER,status,ierr)
      call MPI_FILE_CLOSE(ifile,ierr)
      
      ! Deallocate iblank array
      deallocate(iblank)
      
   end subroutine write_geom
   
   
end module ensight_class
