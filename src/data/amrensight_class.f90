!> Amrensight class concept is defined here: given an amr config object,
!> it provides parallel I/O access to an ensight file
module amrensight_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrconfig_class,  only: amrconfig
   use mpi_f08,          only: MPI_Datatype
   use surfmesh_class,   only: surfmesh
   use partmesh_class,   only: partmesh
   use amrex_amr_module, only: amrex_multifab,amrex_imultifab
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrensight
   
   ! List types
   type :: scl !< Scalar field
      type(scl), pointer :: next
      character(len=str_medium) :: name
      type(amrex_multifab) , dimension(:), pointer :: rptr=>NULL()  !< real(WP) data
      type(amrex_imultifab), dimension(:), pointer :: iptr=>NULL()  !< integer  data
   end type scl
   type :: vct !< Vector field
      type(vct), pointer :: next
      character(len=str_medium) :: name
      type(amrex_multifab) , dimension(:), pointer :: ptrx=>NULL()
      type(amrex_multifab) , dimension(:), pointer :: ptry=>NULL()
      type(amrex_multifab) , dimension(:), pointer :: ptrz=>NULL()
   end type vct
   type :: srf !< Surface mesh
      type(srf), pointer :: next
      character(len=str_medium) :: name
      type(surfmesh), pointer :: ptr
   end type srf
   type :: prt !< Particle mesh
      type(prt), pointer :: next
      character(len=str_medium) :: name
      type(partmesh), pointer :: ptr
   end type prt
   
   !> Amrensight object definition as list of pointers to arrays
   type :: amrensight
      ! An amrensight object has a name
      character(len=str_medium) :: name                               !< Name of amrensight directory to read/write
      ! An amrensight object stores time values
      integer :: ntime                                                !< Number of scalar values
      real(WP), dimension(:), allocatable :: time                     !< Time values
      ! An amrensight object stores geometry data
      type(amrconfig), pointer :: amr                                 !< Amrconfig for ensight geometry and parallel I/O
      ! An amrensight object stores lists of pointers to data
      type(scl), pointer :: first_scl                                 !< Scalar list
      type(vct), pointer :: first_vct                                 !< Vector list
      type(srf), pointer :: first_srf                                 !< Surface list
      type(prt), pointer :: first_prt                                 !< Particle list
   contains
      procedure :: initialize                                         !< Initialization of object
      procedure :: write_geom                                         !< Write out geometry
      procedure :: write_data                                         !< Write out data
      procedure :: write_case                                         !< Write out case file
      procedure :: write_surf                                         !< Write out surface mesh file
      procedure :: write_part                                         !< Write out particle mesh file
      generic :: add_scalar=>add_rscalar,add_iscalar                  !< Add a new scalar field
      procedure, private :: add_rscalar                               !< Add a new real(WP) scalar field
      procedure, private :: add_iscalar                               !< Add a new integer  scalar field
      procedure :: add_vector                                         !< Add a new vector field
      procedure :: add_surface                                        !< Add a new surface mesh
      procedure :: add_particle                                       !< Add a new particle mesh
   end type amrensight
   
   
contains
   
   !> Initialization-in-place for an empty amrensight object
   subroutine initialize(this,amr,name)
      use messager, only: die
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      use filesys,  only: makedir,isdir
      implicit none
      class(amrensight) :: this
      class(amrconfig), target, intent(in) :: amr
      character(len=*), intent(in) :: name
      character(len=str_medium) :: line
      integer :: iunit,ierr,stat
      logical :: file_is_there,found
      
      ! Link to amrconfig
      this%amr=>amr
      
      ! Store casename
      this%name=trim(adjustl(name))
      
      ! Start with no time stamps
      this%ntime=0
      
      ! Create directories
      if (this%amr%amRoot) then
         if (.not.isdir('ensight')) call makedir('ensight')
         if (.not.isdir('ensight/'//trim(this%name))) call makedir('ensight/'//trim(this%name))
         if (.not.isdir('ensight/'//trim(this%name)//'/geometry')) call makedir('ensight/'//trim(this%name)//'/geometry')
      end if
      
      ! Empty pointer to lists for now
      this%first_scl=>NULL()
      this%first_vct=>NULL()
      this%first_srf=>NULL()
      this%first_prt=>NULL()
      
      ! Check if a case file exists already - root only
      if (this%amr%amRoot) then
         inquire(file='ensight/'//trim(this%name)//'/nga.case',exist=file_is_there)
         if (file_is_there) then
            ! Open the case file
            open(newunit=iunit,file='ensight/'//trim(this%name)//'/nga.case',form='formatted',status='old',access='stream',iostat=ierr)
            ! Read lines until we find time values section
            stat=0; found=.false.
            do while (.not.found.and..not.is_iostat_end(stat))
               read(iunit,'(a)',iostat=stat) line
               if (line(1:16).eq.'number of steps:') found=.true.
            end do
            if (found) read(line(17:),'(i6)') this%ntime
            ! Now read the time values
            if (this%ntime.gt.0) then
               allocate(this%time(this%ntime))
               found=.false.
               do while (.not.found)
                  read(iunit,'(a)') line
                  if (line(1:12).eq.'time values:') found=.true.
               end do
               if (found) read(iunit,'(999999(es12.5,/))') this%time
            end if
            ! Close the case file
            close(iunit)
         end if
      end if
      
      ! Communicate to all processors
      call MPI_BCAST(this%ntime,1,MPI_INTEGER,0,this%amr%comm,ierr)
      if (this%ntime.gt.0) then
         if (.not.this%amr%amRoot) allocate(this%time(this%ntime))
         call MPI_BCAST(this%time,this%ntime,MPI_REAL_WP,0,this%amr%comm,ierr)
      end if
      
   end subroutine initialize
   
   
   !> Add a real scalar field for output
   subroutine add_rscalar(this,name,scalar)
      use filesys, only: makedir,isdir
      implicit none
      class(amrensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrex_multifab), dimension(:), target, intent(in) :: scalar
      type(scl), pointer :: new_scl
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%rptr=>scalar
      new_scl%iptr=>NULL()
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
      ! Also create the corresponding directory
      if (this%amr%amRoot) then
         if (.not.isdir('ensight/'//trim(this%name)//trim(new_scl%name))) call makedir('ensight/'//trim(this%name)//trim(new_scl%name))
      end if
   end subroutine add_rscalar
   
   
   !> Add an integer scalar field for output
   subroutine add_iscalar(this,name,scalar)
      use filesys, only: makedir,isdir
      implicit none
      class(amrensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrex_imultifab), dimension(:), target, intent(in) :: scalar
      type(scl), pointer :: new_scl
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%rptr=>NULL()
      new_scl%iptr=>scalar
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
      ! Also create the corresponding directory
      if (this%amr%amRoot) then
         if (.not.isdir('ensight/'//trim(this%name)//trim(new_scl%name))) call makedir('ensight/'//trim(this%name)//trim(new_scl%name))
      end if
   end subroutine add_iscalar
   
   
   !> Add a vector field for output
   subroutine add_vector(this,name,vectx,vecty,vectz)
      use filesys, only: makedir,isdir
      implicit none
      class(amrensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrex_multifab), dimension(:), target, intent(in) :: vectx
      type(amrex_multifab), dimension(:), target, intent(in) :: vecty
      type(amrex_multifab), dimension(:), target, intent(in) :: vectz
      type(vct), pointer :: new_vct
      ! Prepare new vector
      allocate(new_vct)
      new_vct%name=trim(adjustl(name))
      new_vct%ptrx=>vectx
      new_vct%ptry=>vecty
      new_vct%ptrz=>vectz
      ! Insert it up front
      new_vct%next=>this%first_vct
      ! Point list to new object
      this%first_vct=>new_vct
      ! Also create the corresponding directory
      if (this%amr%amRoot) then
         if (.not.isdir('ensight/'//trim(this%name)//trim(new_vct%name))) call makedir('ensight/'//trim(this%name)//trim(new_vct%name))
      end if
   end subroutine add_vector
   
   
   !> Add a surface mesh for output
   subroutine add_surface(this,name,surface)
      use filesys, only: makedir,isdir
      implicit none
      class(amrensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(surfmesh), target, intent(in) :: surface
      type(srf), pointer :: new_srf
      ! Prepare new surface
      allocate(new_srf)
      new_srf%name=trim(adjustl(name))
      new_srf%ptr =>surface
      ! Insert it up front
      new_srf%next=>this%first_srf
      ! Point list to new object
      this%first_srf=>new_srf
      ! Also create the corresponding directory
      if (this%amr%amRoot) then
         if (.not.isdir('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_srf%name))) call makedir('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_srf%name))
      end if
   end subroutine add_surface
   
   
   !> Add a particle mesh for output
   subroutine add_particle(this,name,particle)
      use filesys, only: makedir,isdir
      implicit none
      class(amrensight), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(partmesh), target, intent(in) :: particle
      type(prt), pointer :: new_prt
      ! Prepare new particle
      allocate(new_prt)
      new_prt%name=trim(adjustl(name))
      new_prt%ptr =>particle
      ! Insert it up front
      new_prt%next=>this%first_prt
      ! Point list to new object
      this%first_prt=>new_prt
      ! Also create the corresponding directory
      if (this%amr%amRoot) then
         if (.not.isdir('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_prt%name))) call makedir('mkdir -p ensight/'//trim(this%name)//'/'//trim(new_prt%name))
      end if
   end subroutine add_particle
   
   
   !> Output all data in the object
   subroutine write_data(this,time)
      use precision, only: SP
      use messager,  only: die
      use parallel,  only: info_mpiio,MPI_REAL_SP
      use mpi_f08
      implicit none
      class(amrensight), intent(inout) :: this
      real(WP), intent(in) :: time
      character(len=str_medium) :: filename
      integer :: iunit,ierr,n,i
      integer :: ibuff
      character(len=80) :: cbuff
      type(MPI_File) :: ifile
      integer(kind=MPI_OFFSET_KIND) :: disp
      type(MPI_Status):: status
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      type(srf), pointer :: my_srf
      type(prt), pointer :: my_prt
      real(SP), dimension(:,:,:), allocatable :: spbuff
      real(WP), dimension(:), allocatable :: temp_time
      character(len=str_medium) :: ctime
      real(WP) :: rtime

      ! Check provided time stamp and decide what to do
      if (this%ntime.eq.0) then
         ! First time stamp
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(this%ntime))
         this%time(1)=time
      else
         ! There are time stamps already, check where to insert
         n=1
         rewind: do i=this%ntime,1,-1
            ! Convert time to appropriate accuracy before comparing
            ctime=''; write(ctime,'(es12.5)') time; read(ctime,'(es12.5)') rtime
            if (this%time(i).lt.rtime) then
               n=i+1; exit rewind
            end if
         end do rewind
         this%ntime=n; allocate(temp_time(1:this%ntime))
         temp_time=[this%time(1:this%ntime-1),time]
         call move_alloc(temp_time,this%time)
      end if
      
      ! Write out the geometry
      call this%write_geom()

      !! Prepare the SP buffer
      !allocate(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))
      !
      !! Traverse all datasets and print them all out - scalars first
      !my_scl=>this%first_scl
      !do while (associated(my_scl))
      !   
      !   ! Create filename
      !   filename='ensight/'//trim(this%name)//'/'//trim(my_scl%name)//'/'//trim(my_scl%name)//'.'
      !   write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      !   
      !   ! Root process starts writing the file header
      !   if (this%cfg%amRoot) then
      !      ! Open the file
      !      open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
      !      if (ierr.ne.0) call die('[ensight write data] Could not open file '//trim(filename))
      !      ! Write the header
      !      cbuff=trim(my_scl%name); write(iunit) cbuff
      !      cbuff='part'           ; write(iunit) cbuff
      !      ibuff=1                ; write(iunit) ibuff
      !      cbuff='block'          ; write(iunit) cbuff
      !      ! Close the file
      !      close(iunit)
      !   end if
      !   
      !   ! Now parallel-write the actual data (note that we allow both real and integer fields!)
      !   call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      !   if (ierr.ne.0) call die('[ensight write data] Problem encountered while parallel writing data file '//trim(filename))
      !   call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      !   call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      !   if (associated(my_scl%rptr)) spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_scl%rptr(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
      !   if (associated(my_scl%iptr)) spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_scl%iptr(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
      !   call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      !   call MPI_FILE_CLOSE(ifile,ierr)
      !   
      !   ! Continue on to the next scalar object
      !   my_scl=>my_scl%next
      !   
      !end do
      !
      !! Traverse all datasets and print them all out - vectors second
      !my_vct=>this%first_vct
      !do while (associated(my_vct))
      !   
      !   ! Create filename
      !   filename='ensight/'//trim(this%name)//'/'//trim(my_vct%name)//'/'//trim(my_vct%name)//'.'
      !   write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      !   
      !   ! Root process starts writing the file header
      !   if (this%cfg%amRoot) then
      !      ! Open the file
      !      open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
      !      if (ierr.ne.0) call die('[ensight write data] Could not open file '//trim(filename))
      !      ! Write the header
      !      cbuff=trim(my_vct%name); write(iunit) cbuff
      !      cbuff='part'           ; write(iunit) cbuff
      !      ibuff=1                ; write(iunit) ibuff
      !      cbuff='block'          ; write(iunit) cbuff
      !      ! Close the file
      !      close(iunit)
      !   end if
      !   
      !   ! Now parallel-write the actual data
      !   call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
      !   if (ierr.ne.0) call die('[ensight write data] Problem encountered while parallel writing data file '//trim(filename))
      !   call MPI_FILE_GET_POSITION(ifile,disp,ierr)
      !   call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      !   spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptrx(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
      !   call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      !   disp=disp+int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)*int(SP,MPI_OFFSET_KIND)
      !   call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      !   spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptry(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
      !   call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      !   disp=disp+int(this%cfg%nx,MPI_OFFSET_KIND)*int(this%cfg%ny,MPI_OFFSET_KIND)*int(this%cfg%nz,MPI_OFFSET_KIND)*int(SP,MPI_OFFSET_KIND)
      !   call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      !   spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(my_vct%ptrz(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),SP)
      !   call MPI_FILE_WRITE_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      !   call MPI_FILE_CLOSE(ifile,ierr)
      !   
      !   ! Continue on to the next vector object
      !   my_vct=>my_vct%next
      !   
      !end do
      !
      !! Get rid of the SP buffer
      !deallocate(spbuff)
      
      ! Finally, re-write the case file
      call this%write_case()
      
      ! Now output all surface meshes
      my_srf=>this%first_srf
      do while (associated(my_srf))
         ! Output the surface mesh as a distinct directory
         call this%write_surf(my_srf)
         ! Continue on to the next surface mesh object
         my_srf=>my_srf%next
      end do
      
      ! Now output all particle meshes
      my_prt=>this%first_prt
      do while (associated(my_prt))
         ! Output the particle mesh as a distinct directory
         call this%write_part(my_prt)
         ! Continue on to the next surface mesh object
         my_prt=>my_prt%next
      end do
      
   end subroutine write_data
   
   
   !> Case description serial output to a text file
   subroutine write_case(this)
      implicit none
      class(amrensight), intent(in) :: this
      integer :: iunit,ierr
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      
      ! Only the root does this work
      if (.not.this%amr%amRoot) return
      
      ! Open the case file
      open(newunit=iunit,file='ensight/'//trim(this%name)//'/nga.case',form='formatted',status='replace',access='stream',iostat=ierr)
      
      ! Write all the geometry information
      write(iunit,'(a,/,a,/,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: 1 geometry/geometry.******'
      
      ! Write all the variable information
      write(iunit,'(a)') 'VARIABLE'
      my_scl=>this%first_scl
      do while (associated(my_scl))
         write(iunit,'(a)') 'scalar per element: 1 '//trim(my_scl%name)//' '//trim(my_scl%name)//'/'//trim(my_scl%name)//'.******'
         my_scl=>my_scl%next
      end do
      my_vct=>this%first_vct
      do while (associated(my_vct))
         write(iunit,'(a)') 'vector per element: 1 '//trim(my_vct%name)//' '//trim(my_vct%name)//'/'//trim(my_vct%name)//'.******'
         my_vct=>my_vct%next
      end do
      
      ! Write the time information
      write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
      write(iunit,'(999999(es12.5,/))') this%time
      
      ! Close the case file
      close(iunit)
      
   end subroutine write_case
   
   
   !> Geometry output to a file in parallel
   subroutine write_geom(this)
      use precision,        only: SP
      use messager,         only: die
      use parallel,         only: info_mpiio,MPI_REAL_SP
      use amrex_amr_module, only: amrex_box,amrex_boxarray,amrex_distromap,&
      &                           amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mpi_f08,          only: MPI_BCAST,MPI_INTEGER
      implicit none
      class(amrensight), intent(in) :: this
      integer :: iunit,ierr,n
      type(amrex_boxarray)  :: ba
      type(amrex_distromap) :: dm
      type(amrex_box)       :: bx
      type(amrex_mfiter)    :: mfi
      character(len=str_medium) :: filename
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff,rank,nbox
      
      ! Generate timestamped filename
      filename='ensight/'//trim(this%name)//'/geometry/geometry.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime

      ! Root begins geometry I/O
      if (this%amr%amRoot) then
         ! Open new file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[amrensight write geom] Could not open file '//trim(filename))
         ! General geometry header
         cbuff='C Binary'                  ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'; write(iunit) cbuff
         cbuff=trim(adjustl(this%amr%name)); write(iunit) cbuff
         cbuff='node id off'               ; write(iunit) cbuff
         cbuff='element id off'            ; write(iunit) cbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Write all the boxes in parallel to the same file
      nbox=0
      do n=0,this%amr%clvl()
         ! Build an mfiter at that level
         ba=this%amr%get_boxarray (lvl=n)
         dm=this%amr%get_distromap(lvl=n)
         call amrex_mfiter_build(mfi,ba,dm)
         ! Have all cores write out their boxes
         do rank=0,this%amr%nproc-1
            if (rank.eq.this%amr%rank) then
               ! Open the file
               open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[amrensight write geom] Could not reopen file '//trim(filename))
               ! Loop through all boxes in mfiter
               do while (mfi%next())
                  bx=mfi%tilebox()
                  nbox=nbox+1
                  ! Write box to the file
                  cbuff='part'; write(iunit) cbuff
                  ibuff=nbox  ; write(iunit) ibuff
                  write(cbuff,'("Level ",i0,"/",i0," - Box ",i0," - Rank ",i0)') n,this%amr%clvl(),nbox,rank; write(iunit) cbuff
                  cbuff='block uniform'; write(iunit) cbuff
                  ! Range of indices
                  ibuff=bx%hi(1)+1-bx%lo(1)+1; write(iunit) ibuff
                  ibuff=bx%hi(2)+1-bx%lo(2)+1; write(iunit) ibuff
                  ibuff=bx%hi(3)+1-bx%lo(3)+1; write(iunit) ibuff
                  ! Mesh origin
                  rbuff=real(this%amr%xlo+bx%lo(1)*this%amr%geom(n)%dx(1),SP); write(iunit) rbuff
                  rbuff=real(this%amr%ylo+bx%lo(2)*this%amr%geom(n)%dx(2),SP); write(iunit) rbuff
                  rbuff=real(this%amr%zlo+bx%lo(3)*this%amr%geom(n)%dx(3),SP); write(iunit) rbuff
                  ! Mesh size
                  rbuff=real(this%amr%geom(n)%dx(1),SP); write(iunit) rbuff
                  rbuff=real(this%amr%geom(n)%dx(2),SP); write(iunit) rbuff
                  rbuff=real(this%amr%geom(n)%dx(3),SP); write(iunit) rbuff
               end do
               ! Close the file
               close(iunit)
            end if
            ! Synchronize by broadcasting nbox
            call MPI_BCAST(nbox,1,MPI_INTEGER,rank,this%amr%comm,ierr)
         end do
         ! Destroy iterator
         call amrex_mfiter_destroy(mfi)
      end do
      
   end subroutine write_geom
   
   
   !> Procedure that writes out a surface mesh in Ensight format
   subroutine write_surf(this,surf)
      use precision, only: SP
      use messager,  only: die
      use mpi_f08,   only: mpi_barrier
      implicit none
      class(amrensight), intent(in) :: this
      type(srf), pointer, intent(in) :: surf
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,n
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff
      
      ! Write the case file from scratch in ASCII format
      if (this%amr%amRoot) then
         ! Open the case file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(surf%name)//'.case',form='formatted',status='replace',access='stream',iostat=ierr)
         ! Write all the geometry information
         write(iunit,'(a,/,a,/,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: 1 '//trim(surf%name)//'/'//trim(surf%name)//'.******'
         ! Write the variables
         write(iunit,'(a)') 'VARIABLE'
         do n=1,surf%ptr%nvar
            write(iunit,'(a)') 'scalar per element: 1 '//trim(surf%ptr%varname(n))//' '//trim(surf%name)//'/'//trim(surf%ptr%varname(n))//'.******'
         end do
         ! Write the time information
         write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
         write(iunit,'(999999(es12.5,/))') this%time
         ! Close the case file
         close(iunit)
      end if
      
      ! Generate the surface geometry filename
      filename='ensight/'//trim(this%name)//'/'//trim(surf%name)//'/'//trim(surf%name)//'.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      
      ! Write the file header for Ensight Gold unstructured geometry
      if (this%amr%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'        ; write(iunit) cbuff
         cbuff=trim(adjustl(surf%ptr%name))        ; write(iunit) cbuff
         cbuff='node id off'                       ; write(iunit) cbuff
         cbuff='element id off'                    ; write(iunit) cbuff
         ! Extents
         !cbuff='extents'                           ; write(iunit) cbuff
         !rbuff=real(this%amr%xlo,SP)               ; write(iunit) rbuff
         !rbuff=real(this%amr%xhi,SP)               ; write(iunit) rbuff
         !rbuff=real(this%amr%ylo,SP)               ; write(iunit) rbuff
         !rbuff=real(this%amr%yhi,SP)               ; write(iunit) rbuff
         !rbuff=real(this%amr%zlo,SP)               ; write(iunit) rbuff
         !rbuff=real(this%amr%zhi,SP)               ; write(iunit) rbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Write polygonal mesh in Ensight Gold 'nsided' format
      do rank=0,this%amr%nproc-1
         if (rank.eq.this%amr%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
            ! Part header
            cbuff='part'                              ; write(iunit) cbuff
            ibuff=rank+1                              ; write(iunit) ibuff
            cbuff='Surface geometry per processor #'
            write(cbuff(len_trim(cbuff)+1:len_trim(cbuff)+6),'(i6.6)') this%amr%rank
            write(iunit) cbuff
            ! Write part info if it exists on the processor
            if (surf%ptr%nPoly.gt.0) then
               ! Write out vertices
               cbuff='coordinates'                       ; write(iunit) cbuff
               ibuff=surf%ptr%nVert                      ; write(iunit) ibuff
               write(iunit) real(surf%ptr%xVert,SP)
               write(iunit) real(surf%ptr%yVert,SP)
               write(iunit) real(surf%ptr%zVert,SP)
               ! Write out polygons
               cbuff='nsided'                            ; write(iunit) cbuff
               ibuff=surf%ptr%nPoly                      ; write(iunit) ibuff
               write(iunit) surf%ptr%polySize
               write(iunit) surf%ptr%polyConn
            end if
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%amr%comm,ierr)
      end do
      
      ! Generate the additional variable files
      do n=1,surf%ptr%nvar
         filename='ensight/'//trim(this%name)//'/'//trim(surf%name)//'/'//trim(surf%ptr%varname(n))//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         ! Root write the header
         if (this%amr%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
            ! Write the header
            cbuff=trim(surf%name); write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         ! Write the surface variables
         do rank=0,this%amr%nproc-1
            if (rank.eq.this%amr%rank) then
               ! Open the file
               open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[ensight write surf] Could not open file: '//trim(filename))
               ! Part header
               cbuff='part'         ; write(iunit) cbuff
               ibuff=rank+1         ; write(iunit) ibuff
               ! Write surf info if it exists on the processor
               if (surf%ptr%nPoly.gt.0) then
                  cbuff='nsided'       ; write(iunit) cbuff
                  write(iunit) real(surf%ptr%var(n,:),SP)
               end if
               ! Close the file
               close(iunit)
            end if
            ! Force synchronization
            call MPI_BARRIER(this%amr%comm,ierr)
         end do
      end do
      
   end subroutine write_surf
   
   
   !> Procedure that writes out a particle mesh in Ensight format
   subroutine write_part(this,part)
      use precision, only: SP
      use messager,  only: die
      use mpi_f08,   only: MPI_BARRIER,MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      implicit none
      class(amrensight), intent(in) :: this
      type(prt), pointer, intent(in) :: part
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,n
      character(len=80) :: cbuff
      integer :: ibuff,npart
      
      ! Write the case file from scratch in ASCII format
      if (this%amr%amRoot) then
         ! Open the case file
         open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(part%name)//'.case',form='formatted',status='replace',access='stream',iostat=ierr)
         ! Write all the geometry information
         write(iunit,'(a,/,a,/,/,a,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: geometry','measured: 1 '//trim(part%name)//'/particle.******'
         ! Write the variables
         write(iunit,'(a)') 'VARIABLE'
         write(iunit,'(a)') 'scalar per element: fvf geometry.fvf'
         do n=1,part%ptr%nvar
            write(iunit,'(a)') 'scalar per measured node: 1 '//trim(part%ptr%varname(n))//' '//trim(part%name)//'/'//trim(part%ptr%varname(n))//'.******'
         end do
         do n=1,part%ptr%nvec
            write(iunit,'(a)') 'vector per measured node: 1 '//trim(part%ptr%vecname(n))//' '//trim(part%name)//'/'//trim(part%ptr%vecname(n))//'.******'
         end do
         ! Write the time information
         write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
         write(iunit,'(999999(es12.5,/))') this%time
         ! Close the case file
         close(iunit)
      end if
      
      ! Generate the particle geometry file
      filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/particle.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      ! Count the number of particles
      call MPI_ALLREDUCE(part%ptr%n,npart,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
      ! Root write the header
      if (this%amr%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
         ! General geometry header
         cbuff='C Binary'                          ; write(iunit) cbuff
         cbuff=trim(adjustl(part%ptr%name))        ; write(iunit) cbuff
         cbuff='particle coordinates'              ; write(iunit) cbuff
         ibuff=npart                               ; write(iunit) ibuff
         write(iunit) (ibuff,ibuff=1,npart)
         ! Close the file
         close(iunit)
      end if
      ! Write the particle coordinates
      do rank=0,this%amr%nproc-1
         if (rank.eq.this%amr%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
            ! Write part info if it exists on the processor
            if (part%ptr%n.gt.0) write(iunit) real(part%ptr%pos,SP)
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%amr%comm,ierr)
      end do
      
      ! Generate the particle scalar files
      do n=1,part%ptr%nvar
         filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/'//trim(part%ptr%varname(n))//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         ! Root write the header
         if (this%amr%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
            ! General header
            cbuff='particle '//trim(part%ptr%varname(n)); write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         ! Write the particle variable
         do rank=0,this%amr%nproc-1
            if (rank.eq.this%amr%rank) then
               ! Open the file
               open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
               ! Write part info if it exists on the processor
               if (part%ptr%n.gt.0) write(iunit) real(part%ptr%var(n,:),SP)
               ! Close the file
               close(iunit)
            end if
            ! Force synchronization
            call MPI_BARRIER(this%amr%comm,ierr)
         end do
      end do
      
      ! Generate the particle vector files
      do n=1,part%ptr%nvec
         filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/'//trim(part%ptr%vecname(n))//'.'
         write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
         ! Root write the header
         if (this%amr%amRoot) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
            ! General header
            cbuff='particle'//trim(adjustl(part%ptr%vecname(n)))     ; write(iunit) cbuff
            ! Close the file
            close(iunit)
         end if
         ! Write the vector at the particle location
         do rank=0,this%amr%nproc-1
            if (rank.eq.this%amr%rank) then
               ! Open the file
               open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
               if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
               ! Write part info if it exists on the processor
               if (part%ptr%n.gt.0) write(iunit) real(part%ptr%vec(:,n,:),SP)
               ! Close the file
               close(iunit)
            end if
            ! Force synchronization
            call MPI_BARRIER(this%amr%comm,ierr)
         end do
      end do
      
   end subroutine write_part
   
   
end module amrensight_class
