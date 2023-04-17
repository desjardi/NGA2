!> Definition for a postproc class
module postproc_class
   use precision,         only: WP,SP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use timetracker_class, only: timetracker
   use ccl_class,         only: ccl
   implicit none
   private
   
   public :: postproc
   
   !> postproc object
   type :: postproc
      !> Input file for the simulation
      type(inputfile) :: input
      !> Config
      type(config) :: cfg
      !> Time info
      type(timetracker) :: time  
      !> Ensight postprocessing
      type(partmesh) :: pmesh    !< Particle mesh for core output
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      !> CCL analysis
      type(ccl)  :: cc
      !> Data arrays
      real(WP), dimension(:,:,:), allocatable :: VOF
      real(WP), dimension(:,:,:), allocatable :: U
      real(WP), dimension(:), allocatable :: bv,by,bz
   contains
      procedure :: analyze
      procedure, private :: read_ensight_scalar
      procedure, private :: extract_core
      procedure, private :: extract_drops
   end type postproc
      
contains
   
   
   !> Read a scalar ensight file to an WP array - handle ghost cells as well
   subroutine read_ensight_scalar(this,filename,SC)
      use mpi_f08
      use parallel, only: group,info_mpiio,MPI_REAL_SP
      use string,   only: str_medium
      use messager, only: die
      class(postproc), intent(inout) :: this
      character(len=str_medium), intent(in) :: filename
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: SC
      integer(kind=MPI_OFFSET_KIND) :: disp
      real(SP), dimension(:,:,:), allocatable :: spbuff
      type(MPI_Status):: status
      type(MPI_File) :: ifile
      integer :: ierr,i
      ! Zero out SC
      SC=0.0_WP
      ! Parallel read the file
      call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[postproc read_ensight_scalar] Problem encountered while parallel reading data file '//trim(filename))
      disp=244
      call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_SP,this%cfg%SPview,'native',info_mpiio,ierr)
      allocate(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_))
      call MPI_FILE_READ_ALL(ifile,spbuff,this%cfg%nx_*this%cfg%ny_*this%cfg%nz_,MPI_REAL_SP,status,ierr)
      SC(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=real(spbuff(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_),WP)
      call MPI_FILE_CLOSE(ifile,ierr)
      deallocate(spbuff)
      ! Update ghost cells
      call this%cfg%sync(this%VOF)
      if (this%cfg%iproc.eq.1) then
         do i=this%cfg%imino_,this%cfg%imin_-1
            this%VOF(i,:,:)=this%VOF(this%cfg%imin_,:,:)
         end do
      end if
   end subroutine read_ensight_scalar
   
   
   !> Extract a pmesh skeleton of the liquid core from CCL data
   subroutine extract_core(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: twoPi
      implicit none
      class(postproc), intent(inout) :: this
      real(WP) :: maxvol
      integer :: nbig,n,l,i,j,k,ierr,np
      real(WP), dimension(:), allocatable :: mybv
      real(WP), dimension(:), allocatable :: myby
      real(WP), dimension(:), allocatable :: mybz
      ! Find index of biggest structure
      maxvol=0.0_WP
      do n=1,this%cc%n_meta_struct
         if (this%cc%meta_structures_list(n)%vol.gt.maxvol) then
            maxvol=this%cc%meta_structures_list(n)%vol
            nbig=n
         end if
      end do
      ! Find its barycenter as a function of x
      allocate(mybv(this%cfg%imino:this%cfg%imaxo)); mybv=0.0_WP
      allocate(myby(this%cfg%imino:this%cfg%imaxo)); myby=0.0_WP
      allocate(mybz(this%cfg%imino:this%cfg%imaxo)); mybz=0.0_WP
      do n=this%cc%sync_offset+1,this%cc%sync_offset+this%cc%n_struct
         if (this%cc%struct_list(this%cc%struct_map_(n))%parent.ne.this%cc%meta_structures_list(nbig)%id) cycle
         do l=1,this%cc%struct_list(this%cc%struct_map_(n))%nnode
            i=this%cc%struct_list(this%cc%struct_map_(n))%node(1,l)
            j=this%cc%struct_list(this%cc%struct_map_(n))%node(2,l)
            k=this%cc%struct_list(this%cc%struct_map_(n))%node(3,l)
            ! Integrate barycenter
            mybv(i)=mybv(i)+this%VOF(i,j,k)*this%cfg%vol(i,j,k)
            myby(i)=myby(i)+this%VOF(i,j,k)*this%cfg%vol(i,j,k)*this%cfg%ym(j)
            mybz(i)=mybz(i)+this%VOF(i,j,k)*this%cfg%vol(i,j,k)*this%cfg%zm(k)
         end do
      end do
      this%bv=0.0_WP; call MPI_ALLREDUCE(mybv,this%bv,this%cfg%nxo,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      this%by=0.0_WP; call MPI_ALLREDUCE(myby,this%by,this%cfg%nxo,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      this%bz=0.0_WP; call MPI_ALLREDUCE(mybz,this%bz,this%cfg%nxo,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      np=0
      do i=this%cfg%imino,this%cfg%imaxo
         if (this%bv(i).gt.0.0_WP) then
            this%by(i)=this%by(i)/this%bv(i)
            this%bz(i)=this%bz(i)/this%bv(i)
            np=np+1
         else
            this%by(i)=0.0_WP
            this%bz(i)=0.0_WP
         end if
      end do
      ! Now transfer to pmesh
      if (this%cfg%amRoot) then
         call this%pmesh%reset()
         call this%pmesh%set_size(size=np)
         np=0
         do i=this%cfg%imino,this%cfg%imaxo
            if (this%bv(i).gt.0.0_WP) then
               np=np+1
               this%pmesh%pos(:,np)=[this%cfg%xm(i),this%by(i),this%bz(i)]
               this%pmesh%var(1,np)=this%bv(i)*this%cfg%dxi(i)
               this%pmesh%var(2,np)=sqrt(this%pmesh%var(1,np)/twoPi)
            end if
         end do
      end if
   end subroutine extract_core


   !> Extract a liquid droplets from CCL data
   subroutine extract_drops(this)
      use mathtools, only: Pi
      implicit none
      class(postproc), intent(inout) :: this
      integer :: iunit,ierr,n
      real(WP) :: vol,diam,ecc,lmin,lmax
      ! Only root process works
      if (.not.this%cfg%amRoot) return
      ! Open a text file
      open(newunit=iunit,file='drop.stat',form='formatted',status='old',access='stream',position='append',iostat=ierr)
      ! Traverse structures and output those with x>0.025
      do n=1,this%cc%n_meta_struct
         if (this%cc%meta_structures_list(n)%x.gt.0.025_WP) then
            ! Store volume
            vol=this%cc%meta_structures_list(n)%vol
            ! Compute equivalent diameter
            diam=(6.0_WP*vol/Pi)**(1.0_WP/3.0_WP)
            ! Compute eccentricity
            lmin=this%cc%meta_structures_list(n)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=this%cc%meta_structures_list(n)%lengths(2) ! Handle 2D case
            lmax=this%cc%meta_structures_list(n)%lengths(1)
            ecc=sqrt(1.0_WP-lmin**2/lmax**2)
            ! Output structure to a text file
            write(iunit,'(999999(es12.5,x))') diam,vol,ecc
         end if
      end do
      ! Close the file
      close(iunit)
   end subroutine extract_drops
   
   
   !> Analysis of atom simulation
   subroutine analyze(this,flag)
      use parallel, only: amRoot
      use string,   only: str_medium
      use messager, only: log
      implicit none
      class(postproc), intent(inout) :: this
      logical :: flag
      character(len=str_medium) :: filename
      integer :: nfile

      ! Switch flag to true
      flag=.true.
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_postproc')
      
      ! Create the config
      create_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='postproc')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      ! Allocate work arrays
      allocate_data: block
         allocate(this%VOF(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%VOF=0.0_WP
         allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%U=0.0_WP
         !allocate(this%bv(this%cfg%imino:this%cfg%imaxo)); this%bv=0.0_WP
         !allocate(this%by(this%cfg%imino:this%cfg%imaxo)); this%by=0.0_WP
         !allocate(this%bz(this%cfg%imino:this%cfg%imaxo)); this%bz=0.0_WP
      end block allocate_data
      
      ! Create partmesh and ensight object for core output
      create_ensight: block
         !this%pmesh=partmesh(nvar=2,nvec=0,name='core')
         !this%pmesh%varname(1)='area'
         !this%pmesh%varname(2)='radius'
         !this%ens_out=ensight(cfg=this%cfg,name='postproc')
         !call this%ens_out%add_particle('core',this%pmesh)
      end block create_ensight

      ! Create CCL
      create_ccl: block
         this%cc=ccl(cfg=this%cfg,name='struct')
         this%cc%max_interface_planes=1
         this%cc%VFlo=1.0e-12_WP
      end block create_ccl
      
      ! Run on one file
      !call this%input%read('VOF file',filename)
      !call this%read_ensight_scalar(filename,this%VOF)
      !call this%cc%build_lists(VF=this%VOF,U=this%U,V=this%U,W=this%U)
      !call this%extract_core()
      !call this%cc%deallocate_lists()
      !call this%ens_out%write_data(time=0.0_WP)
      
      ! Run on all files available
      do nfile=1,188
         filename='ensight/atom/VOF/VOF.'; write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') nfile
         call log('Postprocessing file '//trim(filename)//'...')
         call this%read_ensight_scalar(filename,this%VOF)
         call log('|----> File read successfully')
         call this%cc%build_lists(VF=this%VOF,U=this%U,V=this%U,W=this%U)
         call log('|----> CCL analysis done')
         !call this%extract_core()
         !call log('|----> Core extracted')
         call this%extract_drops()
         call log('|----> Drops extracted')
         call this%cc%deallocate_lists()
         !call this%ens_out%write_data(time=real(nfile,WP)*1.0e-3_WP)
         !call log('|----> Ensight output done')
      end do
      
   end subroutine analyze
   

end module postproc_class