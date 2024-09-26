!> Amrensight2 class concept is defined here: given an amr config object and associated data,
!> it enables ensight gold format output to a directory - allows user to select levels to output
module amrensight2_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrconfig_class,  only: amrconfig
   use mpi_f08,          only: MPI_Datatype
   use surfmesh_class,   only: surfmesh
   !use partmesh_class,   only: partmesh
   use amrex_amr_module, only: amrex_multifab,amrex_imultifab
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: amrensight2
   
   ! List types
   type :: scl !< Scalar field
      type(scl), pointer :: next
      character(len=str_medium) :: name
      integer :: comp
      type(amrex_multifab) , dimension(:), pointer :: rptr=>NULL()  !< real(WP) data
      type(amrex_imultifab), dimension(:), pointer :: iptr=>NULL()  !< integer  data
   end type scl
   type :: vct !< Vector field
      type(vct), pointer :: next
      character(len=str_medium) :: name
      integer :: compx,compy,compz
      type(amrex_multifab) , dimension(:), pointer :: ptrx=>NULL()
      type(amrex_multifab) , dimension(:), pointer :: ptry=>NULL()
      type(amrex_multifab) , dimension(:), pointer :: ptrz=>NULL()
   end type vct
   type :: srf !< Surface mesh
      type(srf), pointer :: next
      character(len=str_medium) :: name
      type(surfmesh), pointer :: ptr
   end type srf
   !type :: prt !< Particle mesh
   !   type(prt), pointer :: next
   !   character(len=str_medium) :: name
   !   type(partmesh), pointer :: ptr
   !end type prt
   
   !> Amrensight2 object definition as list of pointers to arrays
   type :: amrensight2
      ! An amrensight2 object has a name
      character(len=str_medium) :: name                               !< Name of amrensight2 directory to read and write case file
      character(len=str_medium), dimension(:), allocatable :: namelvl !< Name of amrensight2 directory including level to write data
      ! An amrensight2 object stores time values
      integer :: ntime                                                !< Number of scalar values
      real(WP), dimension(:), allocatable :: time                     !< Time values
      ! An amrensight2 object stores geometry data
      type(amrconfig), pointer :: amr                                 !< Amrconfig for ensight geometry and parallel I/O
      integer :: nlvlout                                              !< Number of levels to be dumped
      integer, dimension(:), allocatable :: lvlout                    !< Levels to be dumped
      ! An amrensight2 object stores lists of pointers to data
      type(scl), pointer :: first_scl                                 !< Scalar list
      type(vct), pointer :: first_vct                                 !< Vector list
      type(srf), pointer :: first_srf                                 !< Surface list
      !type(prt), pointer :: first_prt                                 !< Particle list
   contains
      procedure :: initialize                                         !< Initialization of object
      procedure :: write_geom                                         !< Write out geometry
      procedure :: write_data                                         !< Write out data
      procedure :: write_case                                         !< Write out case file
      procedure :: write_scalar                                       !< Write out scalar file
      procedure :: write_vector                                       !< Write out vector file
      procedure :: write_surf                                         !< Write out surface mesh file
      !procedure :: write_part                                         !< Write out particle mesh file
      generic :: add_scalar=>add_rscalar,add_iscalar                  !< Add a new scalar field
      procedure, private :: add_rscalar                               !< Add a new real(WP) scalar field
      procedure, private :: add_iscalar                               !< Add a new integer  scalar field
      procedure :: add_vector                                         !< Add a new vector field
      procedure :: add_surface                                        !< Add a new surface mesh
      !procedure :: add_particle                                       !< Add a new particle mesh
   end type amrensight2
   
   
contains
   
   !> Initialization-in-place for an empty amrensight2 object
   subroutine initialize(this,amr,name,lvls)
      use messager, only: die
      use mpi_f08,  only: MPI_BCAST,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      use filesys,  only: makedir,isdir
      implicit none
      class(amrensight2) :: this
      class(amrconfig), target, intent(in) :: amr
      character(len=*), intent(in) :: name
      integer, dimension(:), intent(in), optional :: lvls
      character(len=str_medium) :: line,casename
      integer :: iunit,ierr,stat,n
      logical :: file_is_there,found
      
      ! Link to amrconfig
      this%amr=>amr

      ! Select which levels to output
      if (present(lvls)) then
         ! Check that they are all valid
         do n=1,size(lvls)
            if (lvls(n).lt.0.or.lvls(n).gt.this%amr%nlvl) call die('[amrensight2 initialize] lvls provided are improper')
         end do
         ! If so, store in our object
         this%nlvlout=size(lvls)
         allocate(this%lvlout(this%nlvlout))
         this%lvlout=lvls
      else
         ! By default, include all the levels
         this%nlvlout=this%amr%nlvl+1
         allocate(this%lvlout(this%nlvlout))
         do n=1,this%nlvlout
            this%lvlout(n)=n-1
         end do
      end if
      
      ! Store names
      this%name=trim(adjustl(name))
      allocate(this%namelvl(this%nlvlout))
      do n=1,this%nlvlout
         line=trim(adjustl(name))//'/lev'; write(line(len_trim(line)+1:),'(i0)') this%lvlout(n)
         this%namelvl(n)=line
      end do
      
      ! Start with no time stamps
      this%ntime=0
      
      ! Create ensight, ensight/name, ensight/name/lvls, and ensight/name/lvls/geometry directories
      if (this%amr%amRoot) then
         if (.not.isdir('ensight')) &
         & call makedir('ensight')
         if (.not.isdir('ensight/'//trim(this%name))) &
         & call makedir('ensight/'//trim(this%name))
         do n=1,this%nlvlout
            if (.not.isdir('ensight/'//trim(this%namelvl(n)))) &
            & call makedir('ensight/'//trim(this%namelvl(n)))
            if (.not.isdir('ensight/'//trim(this%namelvl(n))//'/geometry')) &
            & call makedir('ensight/'//trim(this%namelvl(n))//'/geometry')
         end do
      end if
      
      ! Create casefile name for out first level - we'll use it to read in time
      casename='ensight/'//trim(this%name)//'/nga.lev'
      write(casename(len_trim(casename)+1:),'(i0)') this%lvlout(1)
      casename=trim(casename)//'.case'
      
      ! Empty pointer to lists for now
      this%first_scl=>NULL()
      this%first_vct=>NULL()
      this%first_srf=>NULL()
      !this%first_prt=>NULL()
      
      ! Check if a case file exists already - root only
      if (this%amr%amRoot) then
         inquire(file=trim(casename),exist=file_is_there)
         if (file_is_there) then
            ! Open the case file
            open(newunit=iunit,file=trim(casename),form='formatted',status='old',access='stream',iostat=ierr)
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
   subroutine add_rscalar(this,name,scalar,comp)
      use filesys,  only: makedir,isdir
      use messager, only: die
      implicit none
      class(amrensight2), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrex_multifab), dimension(0:), target, intent(in) :: scalar
      integer, intent(in) :: comp
      type(scl), pointer :: new_scl
      integer :: n
      ! Check that the component is meaningful
      if (comp.le.0) call die('[amrensight2 add_rscalar] comp must be at least one')
      if (comp.gt.scalar(0)%nc) call die('[amrensight2 add_rscalar] comp is too large for provided scalar mfab')
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%comp=comp
      new_scl%rptr(0:)=>scalar
      new_scl%iptr=>NULL()
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
      ! Also create the corresponding directories
      if (this%amr%amRoot) then
         do n=1,this%nlvlout
            if (.not.isdir('ensight/'//trim(this%namelvl(n))//'/'//trim(new_scl%name))) &
            & call makedir('ensight/'//trim(this%namelvl(n))//'/'//trim(new_scl%name))
         end do
      end if
   end subroutine add_rscalar
   
   
   !> Add an integer scalar field for output
   subroutine add_iscalar(this,name,scalar,comp)
      use filesys,  only: makedir,isdir
      use messager, only: die
      implicit none
      class(amrensight2), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrex_imultifab), dimension(0:), target, intent(in) :: scalar
      integer, intent(in) :: comp
      type(scl), pointer :: new_scl
      integer :: n
      ! Check that the component is meaningful
      if (comp.le.0) call die('[amrensight2 add_iscalar] comp must be at least one')
      if (comp.gt.scalar(0)%nc) call die('[amrensight2 add_iscalar] comp is too large for provided scalar mfab')
      ! Prepare new scalar
      allocate(new_scl)
      new_scl%name=trim(adjustl(name))
      new_scl%comp=comp
      new_scl%rptr=>NULL()
      new_scl%iptr(0:)=>scalar
      ! Insert it up front
      new_scl%next=>this%first_scl
      ! Point list to new object
      this%first_scl=>new_scl
      ! Also create the corresponding directories
      if (this%amr%amRoot) then
         do n=1,this%nlvlout
            if (.not.isdir('ensight/'//trim(this%namelvl(n))//'/'//trim(new_scl%name))) &
            & call makedir('ensight/'//trim(this%namelvl(n))//'/'//trim(new_scl%name))
         end do
      end if
   end subroutine add_iscalar
   
   
   !> Add a vector field for output
   subroutine add_vector(this,name,vectx,compx,vecty,compy,vectz,compz)
      use filesys,  only: makedir,isdir
      use messager, only: die
      implicit none
      class(amrensight2), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(amrex_multifab), dimension(0:), target, intent(in) :: vectx
      type(amrex_multifab), dimension(0:), target, intent(in) :: vecty
      type(amrex_multifab), dimension(0:), target, intent(in) :: vectz
      integer, intent(in) :: compx,compy,compz
      type(vct), pointer :: new_vct
      integer :: n
      ! Check that the component is meaningful
      if (compx.le.0) call die('[amrensight2 add_vector] compx must be at least one')
      if (compx.gt.vectx(0)%nc) call die('[amrensight2 add_vector] compx is too large for provided vectx mfab')
      if (compy.le.0) call die('[amrensight2 add_vector] compy must be at least one')
      if (compy.gt.vecty(0)%nc) call die('[amrensight2 add_vector] compy is too large for provided vecty mfab')
      if (compz.le.0) call die('[amrensight2 add_vector] compz must be at least one')
      if (compz.gt.vectz(0)%nc) call die('[amrensight2 add_vector] compz is too large for provided vectz mfab')
      ! Prepare new vector
      allocate(new_vct)
      new_vct%name=trim(adjustl(name))
      new_vct%compx=compx
      new_vct%compy=compy
      new_vct%compz=compz
      new_vct%ptrx(0:)=>vectx
      new_vct%ptry(0:)=>vecty
      new_vct%ptrz(0:)=>vectz
      ! Insert it up front
      new_vct%next=>this%first_vct
      ! Point list to new object
      this%first_vct=>new_vct
      ! Also create the corresponding directories
      if (this%amr%amRoot) then
         do n=1,this%nlvlout
            if (.not.isdir('ensight/'//trim(this%namelvl(n))//'/'//trim(new_vct%name))) &
            & call makedir('ensight/'//trim(this%namelvl(n))//'/'//trim(new_vct%name))
         end do
      end if
   end subroutine add_vector
   
   
   !> Add a surface mesh for output
   subroutine add_surface(this,name,surface)
      use filesys, only: makedir,isdir
      implicit none
      class(amrensight2), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(surfmesh), target, intent(in) :: surface
      type(srf), pointer :: new_srf
      integer :: n
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
         do n=1,this%nlvlout
            if (.not.isdir('mkdir -p ensight/'//trim(this%namelvl(n))//'/'//trim(new_srf%name))) &
            & call makedir('mkdir -p ensight/'//trim(this%namelvl(n))//'/'//trim(new_srf%name))
         end do
      end if
   end subroutine add_surface
   
   
   !> Add a particle mesh for output
   !subroutine add_particle(this,name,particle)
   !   use filesys, only: makedir,isdir
   !   implicit none
   !   class(amrensight2), intent(inout) :: this
   !   character(len=*), intent(in) :: name
   !   type(partmesh), target, intent(in) :: particle
   !   type(prt), pointer :: new_prt
   !   ! Prepare new particle
   !   allocate(new_prt)
   !   new_prt%name=trim(adjustl(name))
   !   new_prt%ptr =>particle
   !   ! Insert it up front
   !   new_prt%next=>this%first_prt
   !   ! Point list to new object
   !   this%first_prt=>new_prt
   !   ! Also create the corresponding directory
   !   if (this%amr%amRoot) then
   !      if (.not.isdir('mkdir -p ensight/'//trim(this%namelvl)//'/'//trim(new_prt%name))) call makedir('mkdir -p ensight/'//trim(this%namelvl)//'/'//trim(new_prt%name))
   !   end if
   !end subroutine add_particle
   
   
   !> Output all data in the object
   subroutine write_data(this,time)
      implicit none
      class(amrensight2), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: i,n
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      type(srf), pointer :: my_srf
      !type(prt), pointer :: my_prt
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
      do n=1,this%nlvlout
         call this%write_geom(n)
      end do

      ! Traverse all datasets and print them all out - scalars first
      my_scl=>this%first_scl
      do while (associated(my_scl))
         ! Write out the geometry
         do n=1,this%nlvlout
            call this%write_scalar(my_scl,n)
         end do
         ! Continue on to the next scalar object
         my_scl=>my_scl%next
      end do
      
      ! Traverse all datasets and print them all out - vectors second
      my_vct=>this%first_vct
      do while (associated(my_vct))
         ! Write out the geometry
         do n=1,this%nlvlout
            call this%write_vector(my_vct,n)
         end do
         ! Continue on to the next vector object
         my_vct=>my_vct%next
      end do
      
      ! Finally, re-write the case files
      do n=1,this%nlvlout
         call this%write_case(n)
      end do
      
      ! Now output all surface meshes
      my_srf=>this%first_srf
      do while (associated(my_srf))
         ! Output the surface mesh as a distinct directory
         do n=1,this%nlvlout
            call this%write_surf(my_srf,n)
         end do
         ! Continue on to the next surface mesh object
         my_srf=>my_srf%next
      end do
      
      ! Now output all particle meshes
      !my_prt=>this%first_prt
      !do while (associated(my_prt))
      !   ! Output the particle mesh as a distinct directory
      !   call this%write_part(my_prt)
      !   ! Continue on to the next surface mesh object
      !   my_prt=>my_prt%next
      !end do
      
   end subroutine write_data
   
   
   !> Case description serial output to a text file
   subroutine write_case(this,lvl)
      implicit none
      class(amrensight2), intent(in) :: this
      integer, intent(in) :: lvl
      integer :: iunit,ierr
      type(scl), pointer :: my_scl
      type(vct), pointer :: my_vct
      character(len=str_medium) :: charlvl,casename
      
      ! Only the root does this work
      if (.not.this%amr%amRoot) return

      ! Create casename and charlvl
      write(charlvl,'("lev",i0,"/")') this%lvlout(lvl)
      casename='ensight/'//trim(this%name)//'/nga.lev'
      write(casename(len_trim(casename)+1:),'(i0)') this%lvlout(lvl)
      casename=trim(casename)//'.case'
      
      ! Open the case file
      open(newunit=iunit,file=trim(casename),form='formatted',status='replace',access='stream',iostat=ierr)
      
      ! Write all the geometry information
      write(iunit,'(a,/,a,/,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: 1 '//trim(charlvl)//'geometry/geometry.******'
      
      ! Write all the variable information
      write(iunit,'(a)') 'VARIABLE'
      my_scl=>this%first_scl
      do while (associated(my_scl))
         write(iunit,'(a)') 'scalar per element: 1 '//trim(my_scl%name)//' '//trim(charlvl)//trim(my_scl%name)//'/'//trim(my_scl%name)//'.******'
         my_scl=>my_scl%next
      end do
      my_vct=>this%first_vct
      do while (associated(my_vct))
         write(iunit,'(a)') 'vector per element: 1 '//trim(my_vct%name)//' '//trim(charlvl)//trim(my_vct%name)//'/'//trim(my_vct%name)//'.******'
         my_vct=>my_vct%next
      end do
      
      ! Write the time information
      write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
      write(iunit,'(999999(es12.5,/))') this%time
      
      ! Close the case file
      close(iunit)
      
   end subroutine write_case
   
   
   !> Geometry output to a file in parallel
   subroutine write_geom(this,lvl)
      use precision,        only: SP
      use messager,         only: die
      use parallel,         only: info_mpiio,MPI_REAL_SP
      use mpi_f08,          only: MPI_BCAST,MPI_INTEGER
      use amrex_amr_module, only: amrex_box,amrex_mfiter
      implicit none
      class(amrensight2), intent(in) :: this
      integer, intent(in) :: lvl
      integer :: iunit,ierr
      type(amrex_box)    :: bx
      type(amrex_mfiter) :: mfi
      character(len=str_medium) :: filename
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff,rank,nbox
      
      ! Generate timestamped filename
      filename='ensight/'//trim(this%namelvl(lvl))//'/geometry/geometry.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      
      ! Root begins geometry I/O
      if (this%amr%amRoot) then
         ! Open new file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[amrensight2 write geom] Could not open file '//trim(filename))
         ! General geometry header
         cbuff='C Binary'                  ; write(iunit) cbuff
         cbuff='Ensight Gold Geometry File'; write(iunit) cbuff
         cbuff=trim(adjustl(this%amr%name))//' - Level'; write(cbuff(len_trim(cbuff)+2:),'(i0)') this%lvlout(lvl); write(iunit) cbuff
         cbuff='node id off'               ; write(iunit) cbuff
         cbuff='element id off'            ; write(iunit) cbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Write all the boxes on our level in parallel to the same file
      nbox=0
      ! Build an mfiter at our level
      call this%amr%mfiter_build(this%lvlout(lvl),mfi)
      ! Have all cores write out their boxes
      do rank=0,this%amr%nproc-1
         if (rank.eq.this%amr%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[amrensight2 write geom] Could not reopen file '//trim(filename))
            ! Loop through all boxes in mfiter
            do while (mfi%next())
               bx=mfi%tilebox()
               nbox=nbox+1
               ! Write box to the file
               cbuff='part'; write(iunit) cbuff
               ibuff=nbox  ; write(iunit) ibuff
               write(cbuff,'("Box ",i0," - Rank ",i0)') nbox,rank; write(iunit) cbuff
               cbuff='block uniform'; write(iunit) cbuff
               ! Range of indices
               ibuff=bx%hi(1)+1-bx%lo(1)+1; write(iunit) ibuff
               ibuff=bx%hi(2)+1-bx%lo(2)+1; write(iunit) ibuff
               ibuff=bx%hi(3)+1-bx%lo(3)+1; write(iunit) ibuff
               ! Mesh origin
               rbuff=real(this%amr%xlo+bx%lo(1)*this%amr%geom(this%lvlout(lvl))%dx(1),SP); write(iunit) rbuff
               rbuff=real(this%amr%ylo+bx%lo(2)*this%amr%geom(this%lvlout(lvl))%dx(2),SP); write(iunit) rbuff
               rbuff=real(this%amr%zlo+bx%lo(3)*this%amr%geom(this%lvlout(lvl))%dx(3),SP); write(iunit) rbuff
               ! Mesh size
               rbuff=real(this%amr%geom(this%lvlout(lvl))%dx(1),SP); write(iunit) rbuff
               rbuff=real(this%amr%geom(this%lvlout(lvl))%dx(2),SP); write(iunit) rbuff
               rbuff=real(this%amr%geom(this%lvlout(lvl))%dx(3),SP); write(iunit) rbuff
            end do
            ! Close the file
            close(iunit)
         end if
         ! Synchronize by broadcasting nbox
         call MPI_BCAST(nbox,1,MPI_INTEGER,rank,this%amr%comm,ierr)
      end do
      ! Destroy iterator
      call this%amr%mfiter_destroy(mfi)
      
   end subroutine write_geom
   
   
   !> Procedure that writes out a scalar in Ensight format
   subroutine write_scalar(this,my_scl,lvl)
      use precision,        only: SP
      use messager,         only: die
      use parallel,         only: info_mpiio,MPI_REAL_SP
      use mpi_f08,          only: MPI_BCAST,MPI_INTEGER
      use amrex_amr_module, only: amrex_box,amrex_mfiter
      implicit none
      class(amrensight2), intent(inout) :: this
      type(scl), pointer, intent(in) :: my_scl
      integer, intent(in) :: lvl
      character(len=str_medium) :: filename
      integer :: iunit,ierr,nbox,rank,ibuff
      character(len=80) :: cbuff
      type(amrex_box)    :: bx
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rphi
      integer , dimension(:,:,:,:), contiguous, pointer :: iphi
      real(SP), dimension(:,:,:), allocatable :: spbuff
      
      ! Create filename
      filename='ensight/'//trim(this%namelvl(lvl))//'/'//trim(my_scl%name)//'/'//trim(my_scl%name)//'.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      
      ! Root process starts writing the file header
      if (this%amr%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[amrensight2 write scalar] Could not open file '//trim(filename))
         ! Write the header
         cbuff=trim(my_scl%name); write(iunit) cbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Write all the boxes on our level in parallel to the same file
      nbox=0
      ! Build an mfiter at our level
      call this%amr%mfiter_build(this%lvlout(lvl),mfi)
      ! Have all cores write out their data
      do rank=0,this%amr%nproc-1
         if (rank.eq.this%amr%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[amrensight2 write data] Could not reopen file '//trim(filename))
            ! Loop through all boxes in mfiter
            do while (mfi%next())
               ! Get box
               bx=mfi%tilebox()
               nbox=nbox+1
               ! Write box to the file
               cbuff='part' ; write(iunit) cbuff
               ibuff=nbox   ; write(iunit) ibuff
               cbuff='block'; write(iunit) cbuff
               ! Allocate SPbuff to right size
               allocate(spbuff(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3)))
               ! Copy data from appropriate multifab
               if (associated(my_scl%rptr)) then
                  rphi=>my_scl%rptr(this%lvlout(lvl))%dataptr(mfi)
                  spbuff=real(rphi(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),my_scl%comp),SP)
               else if (associated(my_scl%iptr)) then
                  iphi=>my_scl%iptr(this%lvlout(lvl))%dataptr(mfi)
                  spbuff=real(iphi(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),my_scl%comp),SP)
               end if
               ! Write it out
               write(iunit) spbuff
               ! Deallocate SPbuff
               deallocate(spbuff)
            end do
            ! Close the file
            close(iunit)
         end if
         ! Synchronize by broadcasting nbox
         call MPI_BCAST(nbox,1,MPI_INTEGER,rank,this%amr%comm,ierr)
      end do
      ! Destroy iterator
      call this%amr%mfiter_destroy(mfi)
      
   end subroutine write_scalar


   !> Procedure that writes out a vector in Ensight format
   subroutine write_vector(this,my_vct,lvl)
      use precision,        only: SP
      use messager,         only: die
      use parallel,         only: info_mpiio,MPI_REAL_SP
      use mpi_f08,          only: MPI_BCAST,MPI_INTEGER
      use amrex_amr_module, only: amrex_box,amrex_mfiter
      implicit none
      class(amrensight2), intent(inout) :: this
      type(vct), pointer, intent(in) :: my_vct
      integer, intent(in) :: lvl
      character(len=str_medium) :: filename
      integer :: iunit,ierr,nbox,rank,ibuff
      character(len=80) :: cbuff
      type(amrex_box)    :: bx
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: rphi
      real(SP), dimension(:,:,:), allocatable :: spbuff
      
      ! Create filename
      filename='ensight/'//trim(this%namelvl(lvl))//'/'//trim(my_vct%name)//'/'//trim(my_vct%name)//'.'
      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
      
      ! Root process starts writing the file header
      if (this%amr%amRoot) then
         ! Open the file
         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[amrensight2 write vector] Could not open file '//trim(filename))
         ! Write the header
         cbuff=trim(my_vct%name); write(iunit) cbuff
         ! Close the file
         close(iunit)
      end if
      
      ! Write all the boxes on our level in parallel to the same file
      nbox=0
      ! Build an mfiter at our level
      call this%amr%mfiter_build(this%lvlout(lvl),mfi)
      ! Have all cores write out their data
      do rank=0,this%amr%nproc-1
         if (rank.eq.this%amr%rank) then
            ! Open the file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[amrensight2 write data] Could not reopen file '//trim(filename))
            ! Loop through all boxes in mfiter
            do while (mfi%next())
               ! Get box
               bx=mfi%tilebox()
               nbox=nbox+1
               ! Write box to the file
               cbuff='part' ; write(iunit) cbuff
               ibuff=nbox   ; write(iunit) ibuff
               cbuff='block'; write(iunit) cbuff
               ! Allocate SPbuff to right size
               allocate(spbuff(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3)))
               ! Copy data from appropriate multifab and write it out
               rphi=>my_vct%ptrx(this%lvlout(lvl))%dataptr(mfi)
               spbuff=real(rphi(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),my_vct%compx),SP)
               write(iunit) spbuff
               rphi=>my_vct%ptry(this%lvlout(lvl))%dataptr(mfi)
               spbuff=real(rphi(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),my_vct%compy),SP)
               write(iunit) spbuff
               rphi=>my_vct%ptrz(this%lvlout(lvl))%dataptr(mfi)
               spbuff=real(rphi(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),my_vct%compz),SP)
               write(iunit) spbuff
               ! Deallocate SPbuff
               deallocate(spbuff)
            end do
            ! Close the file
            close(iunit)
         end if
         ! Synchronize by broadcasting nbox
         call MPI_BCAST(nbox,1,MPI_INTEGER,rank,this%amr%comm,ierr)
      end do
      ! Destroy iterator
      call this%amr%mfiter_destroy(mfi)
      
   end subroutine write_vector
   
   
   !> Procedure that writes out a surface mesh in Ensight format
   subroutine write_surf(this,surf,lvl)
      use precision, only: SP
      use messager,  only: die
      use mpi_f08,   only: mpi_barrier
      implicit none
      class(amrensight2), intent(in) :: this
      type(srf), pointer, intent(in) :: surf
      integer, intent(in) :: lvl
      character(len=str_medium) :: filename,charlvl
      integer :: iunit,ierr,rank,n
      character(len=80) :: cbuff
      real(SP) :: rbuff
      integer :: ibuff
      
      ! Write the case file from scratch in ASCII format
      if (this%amr%amRoot) then
         ! Open the case file
         filename='ensight/'//trim(this%name)//'/'//trim(surf%name)//'.lev'
         write(filename(len_trim(filename)+1:),'(i0)') this%lvlout(lvl)
         filename=trim(filename)//'.case'
         open(newunit=iunit,file=trim(filename),form='formatted',status='replace',access='stream',iostat=ierr)
         ! Create string with level
         write(charlvl,'("lev",i0,"/")') this%lvlout(lvl)
         ! Write all the geometry information
         write(iunit,'(a,/,a,/,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: 1 '//trim(charlvl)//trim(surf%name)//'/'//trim(surf%name)//'.******'
         ! Write the variables
         write(iunit,'(a)') 'VARIABLE'
         do n=1,surf%ptr%nvar
            write(iunit,'(a)') 'scalar per element: 1 '//trim(surf%ptr%varname(n))//' '//trim(charlvl)//trim(surf%name)//'/'//trim(surf%ptr%varname(n))//'.******'
         end do
         ! Write the time information
         write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
         write(iunit,'(999999(es12.5,/))') this%time
         ! Close the case file
         close(iunit)
      end if
      
      ! Generate the surface geometry filename
      filename='ensight/'//trim(this%namelvl(lvl))//'/'//trim(surf%name)//'/'//trim(surf%name)//'.'
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
         filename='ensight/'//trim(this%namelvl(lvl))//'/'//trim(surf%name)//'/'//trim(surf%ptr%varname(n))//'.'
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
   !subroutine write_part(this,part)
   !   use precision, only: SP
   !   use messager,  only: die
   !   use mpi_f08,   only: MPI_BARRIER,MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
   !   implicit none
   !   class(amrensight2), intent(in) :: this
   !   type(prt), pointer, intent(in) :: part
   !   character(len=str_medium) :: filename
   !   integer :: iunit,ierr,rank,n
   !   character(len=80) :: cbuff
   !   integer :: ibuff,npart
   !   
   !   ! Write the case file from scratch in ASCII format
   !   if (this%amr%amRoot) then
   !      ! Open the case file
   !      open(newunit=iunit,file='ensight/'//trim(this%name)//'/'//trim(part%name)//'.case',form='formatted',status='replace',access='stream',iostat=ierr)
   !      ! Write all the geometry information
   !      write(iunit,'(a,/,a,/,/,a,/,a,/,a,/)') 'FORMAT','type: ensight gold','GEOMETRY','model: geometry','measured: 1 '//trim(part%name)//'/particle.******'
   !      ! Write the variables
   !      write(iunit,'(a)') 'VARIABLE'
   !      write(iunit,'(a)') 'scalar per element: fvf geometry.fvf'
   !      do n=1,part%ptr%nvar
   !         write(iunit,'(a)') 'scalar per measured node: 1 '//trim(part%ptr%varname(n))//' '//trim(part%name)//'/'//trim(part%ptr%varname(n))//'.******'
   !      end do
   !      do n=1,part%ptr%nvec
   !         write(iunit,'(a)') 'vector per measured node: 1 '//trim(part%ptr%vecname(n))//' '//trim(part%name)//'/'//trim(part%ptr%vecname(n))//'.******'
   !      end do
   !      ! Write the time information
   !      write(iunit,'(/,a,/,a,/,a,i0,/,a,/,a,/,a)') 'TIME','time set: 1','number of steps: ',this%ntime,'filename start number: 1','filename increment: 1','time values:'
   !      write(iunit,'(999999(es12.5,/))') this%time
   !      ! Close the case file
   !      close(iunit)
   !   end if
   !   
   !   ! Generate the particle geometry file
   !   filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/particle.'
   !   write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
   !   ! Count the number of particles
   !   call MPI_ALLREDUCE(part%ptr%n,npart,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
   !   ! Root write the header
   !   if (this%amr%amRoot) then
   !      ! Open the file
   !      open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
   !      if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
   !      ! General geometry header
   !      cbuff='C Binary'                          ; write(iunit) cbuff
   !      cbuff=trim(adjustl(part%ptr%name))        ; write(iunit) cbuff
   !      cbuff='particle coordinates'              ; write(iunit) cbuff
   !      ibuff=npart                               ; write(iunit) ibuff
   !      write(iunit) (ibuff,ibuff=1,npart)
   !      ! Close the file
   !      close(iunit)
   !   end if
   !   ! Write the particle coordinates
   !   do rank=0,this%amr%nproc-1
   !      if (rank.eq.this%amr%rank) then
   !         ! Open the file
   !         open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
   !         if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
   !         ! Write part info if it exists on the processor
   !         if (part%ptr%n.gt.0) write(iunit) real(part%ptr%pos,SP)
   !         ! Close the file
   !         close(iunit)
   !      end if
   !      ! Force synchronization
   !      call MPI_BARRIER(this%amr%comm,ierr)
   !   end do
   !   
   !   ! Generate the particle scalar files
   !   do n=1,part%ptr%nvar
   !      filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/'//trim(part%ptr%varname(n))//'.'
   !      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
   !      ! Root write the header
   !      if (this%amr%amRoot) then
   !         ! Open the file
   !         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
   !         if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
   !         ! General header
   !         cbuff='particle '//trim(part%ptr%varname(n)); write(iunit) cbuff
   !         ! Close the file
   !         close(iunit)
   !      end if
   !      ! Write the particle variable
   !      do rank=0,this%amr%nproc-1
   !         if (rank.eq.this%amr%rank) then
   !            ! Open the file
   !            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
   !            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
   !            ! Write part info if it exists on the processor
   !            if (part%ptr%n.gt.0) write(iunit) real(part%ptr%var(n,:),SP)
   !            ! Close the file
   !            close(iunit)
   !         end if
   !         ! Force synchronization
   !         call MPI_BARRIER(this%amr%comm,ierr)
   !      end do
   !   end do
   !   
   !   ! Generate the particle vector files
   !   do n=1,part%ptr%nvec
   !      filename='ensight/'//trim(this%name)//'/'//trim(part%name)//'/'//trim(part%ptr%vecname(n))//'.'
   !      write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') this%ntime
   !      ! Root write the header
   !      if (this%amr%amRoot) then
   !         ! Open the file
   !         open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
   !         if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
   !         ! General header
   !         cbuff='particle'//trim(adjustl(part%ptr%vecname(n)))     ; write(iunit) cbuff
   !         ! Close the file
   !         close(iunit)
   !      end if
   !      ! Write the vector at the particle location
   !      do rank=0,this%amr%nproc-1
   !         if (rank.eq.this%amr%rank) then
   !            ! Open the file
   !            open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
   !            if (ierr.ne.0) call die('[ensight write part] Could not open file: '//trim(filename))
   !            ! Write part info if it exists on the processor
   !            if (part%ptr%n.gt.0) write(iunit) real(part%ptr%vec(:,n,:),SP)
   !            ! Close the file
   !            close(iunit)
   !         end if
   !         ! Force synchronization
   !         call MPI_BARRIER(this%amr%comm,ierr)
   !      end do
   !   end do
   !   
   !end subroutine write_part
   
   
end module amrensight2_class
