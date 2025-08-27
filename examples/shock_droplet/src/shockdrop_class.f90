!> Definition of a shock-drop problem
module shockdrop_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use mpcomp_class,      only: mpcomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use cclabel_class,     only: cclabel
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: shockdrop
   
   !> Shockdrop object
   type :: shockdrop
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(mpcomp) :: fs        !< Multiphase compressible solver
      type(timetracker) :: time !< Time info
      
      !> CCL for postprocessing
      type(cclabel) :: ccl
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh
      type(ensight)  :: ens_out
      
      !> Simulation monitor file
      type(monitor) :: mfile,cflfile,consfile,dropfile,meshfile
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
      real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc
      
      !> Constant phasic kinematic viscosities
      real(WP) :: cst_viscL,cst_viscG
      
      !> Various post-processing info
      real(WP) :: Vcore,Mcore,Xcore,Ycore,Zcore !< Drop core data
      real(WP), dimension(3) :: Cmin,Cmax       !< Core extent
      
      
   contains
      procedure :: initialize                      !< Initialize shock-drop simulation
      procedure :: step                            !< Advance shock-drop simulation by one time step
      procedure, private :: postproc               !< Postprocess shock-drop case
      procedure :: analyze_drops                   !< Output droplet analysis
      procedure :: output_monitor                  !< Monitoring for shock-drop case
      procedure :: output_ensight                  !< Ensight output for shock-drop case
      procedure, private :: prepare_viscosities    !< Prepare viscosities
      procedure :: apply_bconds                    !< Apply boundary conditions
      procedure :: finalize                        !< Finalize shock-drop simulation
   end type shockdrop
   
contains
   
   
   !> Various postprocessing
   subroutine postproc(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MIN,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(shockdrop), intent(inout) :: this
      integer :: i,j,k,n,ierr
      real(WP), parameter :: VFlo_extent=0.1_WP
      ! Update CCL
      call this%ccl%build(make_label,same_label)
      ! Core is id=1, skip core analysis if no core is present
      if (this%ccl%nstruct.ge.1) then
         ! Extract core volume, mass, barycenter, and extent
         this%Vcore=0.0_WP; this%Mcore=0.0_WP; this%Xcore=0.0_WP; this%Ycore=0.0_WP; this%Zcore=0.0_WP
         this%Cmin=+[huge(1.0_WP),huge(1.0_WP),huge(1.0_WP)]; this%Cmax=-[huge(1.0_WP),huge(1.0_WP),huge(1.0_WP)]
         do n=1,this%ccl%struct(1)%n_
            ! Get cell index
            i=this%ccl%struct(1)%map(1,n); j=this%ccl%struct(1)%map(2,n); k=this%ccl%struct(1)%map(3,n)
            ! Increment volume
            this%Vcore=this%Vcore+this%fs%VF(i,j,k)*this%fs%cfg%vol(i,j,k)
            ! Increment mass
            this%Mcore=this%Mcore+this%fs%Q(i,j,k,1)*this%fs%cfg%vol(i,j,k)
            ! Increment barycenter
            this%Xcore=this%Xcore+this%fs%Q(i,j,k,1)*this%fs%BL(1,i,j,k)*this%fs%cfg%vol(i,j,k)
            this%Ycore=this%Ycore+this%fs%Q(i,j,k,1)*this%fs%BL(2,i,j,k)*this%fs%cfg%vol(i,j,k)
            this%Zcore=this%Zcore+this%fs%Q(i,j,k,1)*this%fs%BL(3,i,j,k)*this%fs%cfg%vol(i,j,k)
            ! Increment extent
            this%Cmin=min(this%Cmin,[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )])
            this%Cmax=max(this%Cmax,[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k+1)])
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vcore,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Mcore,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Xcore,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%Xcore=this%Xcore/this%Mcore ! Shouldn't ever
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Ycore,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%Ycore=this%Ycore/this%Mcore ! be dividing by
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Zcore,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%Zcore=this%Zcore/this%Mcore ! zero here...
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Cmin ,3,MPI_REAL_WP,MPI_MIN,this%fs%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Cmax ,3,MPI_REAL_WP,MPI_MAX,this%fs%cfg%comm,ierr)
      end if
   contains
      !> Function that identifies cells that need a label
      logical function make_label(i1,j1,k1)
         implicit none
         integer, intent(in) :: i1,j1,k1
         if (this%fs%VF(i1,j1,k1).gt.0.0_WP) then; make_label=.true.; else; make_label=.false.; end if
      end function make_label
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         use irl_fortran_interface, only: calculateNormal,calculateCentroid
         implicit none
         integer , intent(in) :: i1,j1,k1,i2,j2,k2
         real(WP), dimension(3) :: N1,N2,O1,O2
         ! Big default, assume same label
         same_label=.true.
         ! Look more closely at the local polygon alignment to decide whether to use same labels
         if (this%fs%VF(i1,j1,k1).gt.0.0_WP.and.this%fs%VF(i1,j1,k1).lt.1.0_WP.and.this%fs%VF(i2,j2,k2).gt.0.0_WP.and.this%fs%VF(i2,j2,k2).lt.1.0_WP) then
            ! Get polygon normals
            N1=calculateNormal(this%fs%interface_polygon(i1,j1,k1))
            N2=calculateNormal(this%fs%interface_polygon(i2,j2,k2))
            ! If not at least ~75 degrees, return
            if (dot_product(N1,N2).ge.0.3_WP) return
            ! Get polygon barycenters
            O1=calculateCentroid(this%fs%interface_polygon(i1,j1,k1))
            O2=calculateCentroid(this%fs%interface_polygon(i2,j2,k2))
            ! If pointing towards one another, use different labels
            if (dot_product(O1-O2,N1).lt.0.0_WP.and.dot_product(O2-O1,N2).lt.0.0_WP) same_label=.false.
         end if
      end function same_label
   end subroutine postproc
   
   
   !> Analysis of droplets
   subroutine analyze_drops(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP,amRoot
      use mathtools, only: Pi
      use string,    only: str_medium
      use filesys,   only: makedir,isdir
      implicit none
      class(shockdrop), intent(inout) :: this
      integer :: i,j,k,n,m,ierr,iunit
      real(WP), dimension(:)  , allocatable :: Vd,Md,Pd
      real(WP), dimension(:,:), allocatable :: Bd,Ud
      character(len=str_medium) :: filename,timestamp
      
      ! Allocate volume, mass, pressure, barycenter, and velocity arrays
      allocate(Vd(    1:this%ccl%nstruct)); Vd=0.0_WP
      allocate(Md(    1:this%ccl%nstruct)); Md=0.0_WP
      allocate(Pd(    1:this%ccl%nstruct)); Pd=0.0_WP
      allocate(Bd(1:3,1:this%ccl%nstruct)); Bd=0.0_WP
      allocate(Ud(1:3,1:this%ccl%nstruct)); Ud=0.0_WP
      
      ! Loop over individual structures and compute structure properties
      do n=1,this%ccl%nstruct
         ! Loop over cells in structure and accumulate data
         do m=1,this%ccl%struct(n)%n_
            ! Get cell index
            i=this%ccl%struct(n)%map(1,m); j=this%ccl%struct(n)%map(2,m); k=this%ccl%struct(n)%map(3,m)
            ! Increment volume, mass, pressure, barycenter, and momentum
            Vd  (n)=Vd  (n)+this%fs%cfg%vol(i,j,k)*this%fs%VF(i,j,k)
            Md  (n)=Md  (n)+this%fs%cfg%vol(i,j,k)*this%fs%Q (i,j,k,1)
            Pd  (n)=Pd  (n)+this%fs%cfg%vol(i,j,k)*this%fs%VF(i,j,k)  *this%fs%PL(i,j,k)
            Bd(:,n)=Bd(:,n)+this%fs%cfg%vol(i,j,k)*this%fs%Q (i,j,k,1)*this%fs%BL(:,i,j,k)
            Ud(1,n)=Ud(1,n)+this%fs%cfg%vol(i,j,k)*this%fs%Q (i,j,k,1)*this%Ui(i,j,k)
            Ud(2,n)=Ud(2,n)+this%fs%cfg%vol(i,j,k)*this%fs%Q (i,j,k,1)*this%Vi(i,j,k)
            Ud(3,n)=Ud(3,n)+this%fs%cfg%vol(i,j,k)*this%fs%Q (i,j,k,1)*this%Wi(i,j,k)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,Vd,  this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Md,  this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Pd,  this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); Pd=Pd/Vd
      call MPI_ALLREDUCE(MPI_IN_PLACE,Bd,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); do n=1,this%ccl%nstruct; Bd(:,n)=Bd(:,n)/Md(n); end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,Ud,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); do n=1,this%ccl%nstruct; Ud(:,n)=Ud(:,n)/Md(n); end do
      
      ! Only root process outputs to a file
      if (this%cfg%amRoot) then
         ! Ensure we have a directory
         if (.not.isdir('droplets')) call makedir('droplets')
         ! Create filename
         filename='droplets_'; write(timestamp,'(es12.5)') this%time%t
         ! Open droplet file
         open(newunit=iunit,file='droplets/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         ! Write the header
         write(iunit,'(a12,10(3x,a12))') 'Label','Diameter','Volume','Mass','Pressure','X','Y','Z','U','V','W'
         ! Output the droplet data
         do n=1,this%ccl%nstruct
            write(iunit,'(i12,10(3x,es12.5))') n,(6.0_WP*Vd(n)/Pi)**(1.0_WP/3.0_WP),Vd(n),Md(n),Pd(n),Bd(1,n),Bd(2,n),Bd(3,n),Ud(1,n),Ud(2,n),Ud(3,n)
         end do
         ! Close the file
         close(iunit)
      end if
      
      ! Deallocate memory
      deallocate(Vd,Md,Pd,Bd,Ud)
      
   end subroutine analyze_drops
   
   
   !> Initialization of a shock-drop problem
   subroutine initialize(this,dx,meshsize,startloc,group,partition,continue_monitor)
      use mpi_f08, only: MPI_Group
      implicit none
      class(shockdrop), intent(inout) :: this
      real(WP), intent(in) :: dx
      integer , dimension(3), intent(in) :: meshsize
      real(WP), dimension(3), intent(in) :: startloc
      type(MPI_Group)       , intent(in) :: group
      integer , dimension(3), intent(in) :: partition
      logical , optional    , intent(in) :: continue_monitor
      logical :: monitor_continue
      
      ! Check if this shockdrop is a continuation, coming from remeshing
      is_monitor_continued: block
         monitor_continue=.false.; if (present(continue_monitor)) monitor_continue=continue_monitor
      end block is_monitor_continued
      
      ! Initialize config object
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         type(sgrid) :: grid
         integer  :: i,j,k
         logical  :: xper,yper,zper
         real(WP), dimension(:), allocatable :: x,y,z
         ! By default, the case is fully non-periodic
         xper=.false.; yper=.false.; zper=.false.
         ! Handle 2D cases
         if (meshsize(1).eq.1) xper=.true.
         if (meshsize(2).eq.1) yper=.true.
         if (meshsize(3).eq.1) zper=.true.
         ! Create simple rectilinear grid
         allocate(x(meshsize(1)+1)); do i=1,meshsize(1)+1; x(i)=real(i-1,WP)*dx+startloc(1); end do
         allocate(y(meshsize(2)+1)); do j=1,meshsize(2)+1; y(j)=real(j-1,WP)*dx+startloc(2); end do
         allocate(z(meshsize(3)+1)); do k=1,meshsize(3)+1; z(k)=real(k-1,WP)*dx+startloc(3); end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=xper,yper=yper,zper=zper,name='ShockDrop')
         ! Deallocate x/y/z
         deallocate(x,y,z)
         ! Create config
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      ! Initialize time tracker
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot,name='ShockDrop',print_info=.false.)
      end block initialize_timetracker
      
      ! Create multiphase compressible flow solver
      create_velocity_solver: block
         call this%fs%initialize(cfg=this%cfg,name='Compressible NS')
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%dQdt(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%fs%nQ,1:4))
         allocate(this%beta(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ma(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Prepare post-processing
      prep_postprocess: block
         call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      end block prep_postprocess
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='ShockDrop',time_dependent_geometry=.true.)
         ! No need to output FVF file
         this%ens_out%write_fvf=.false.
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%fs%VF)
         call this%ens_out%add_scalar('RHOL',this%fs%RHOL)
         call this%ens_out%add_scalar('RHOG',this%fs%RHOG)
         call this%ens_out%add_scalar('IL',this%fs%IL)
         call this%ens_out%add_scalar('IG',this%fs%IG)
         call this%ens_out%add_scalar('PL',this%fs%PL)
         call this%ens_out%add_scalar('PG',this%fs%PG)
         call this%ens_out%add_scalar('Mach',this%Ma)
         call this%ens_out%add_scalar('beta',this%beta)
         call this%ens_out%add_scalar('visc',this%visc)
         call this%ens_out%add_scalar('label',this%ccl%id)
         ! Create surface mesh for PLIC
         this%smesh=surfmesh(nvar=1,name='plic')
         this%smesh%varname(1)='label'
         call this%ens_out%add_surface('plic',this%smesh)
      end block create_ensight
      
      ! Create monitor files
      create_monitor: block
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'shockdrop_sim',restart=monitor_continue)
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%RHOLmax,'max(RHOL)')
         call this%mfile%add_column(this%fs%RHOLmin,'min(RHOL)')
         call this%mfile%add_column(this%fs%ILmax  ,'max(IL)'  )
         call this%mfile%add_column(this%fs%ILmin  ,'min(IL)'  )
         call this%mfile%add_column(this%fs%PLmax  ,'max(PL)'  )
         call this%mfile%add_column(this%fs%PLmin  ,'min(PL)'  )
         call this%mfile%add_column(this%fs%TLmax  ,'max(TL)'  )
         call this%mfile%add_column(this%fs%TLmin  ,'min(TL)'  )
         call this%mfile%add_column(this%fs%RHOGmax,'max(RHOG)')
         call this%mfile%add_column(this%fs%RHOGmin,'min(RHOG)')
         call this%mfile%add_column(this%fs%IGmax  ,'max(IG)'  )
         call this%mfile%add_column(this%fs%IGmin  ,'min(IG)'  )
         call this%mfile%add_column(this%fs%PGmax  ,'max(PG)'  )
         call this%mfile%add_column(this%fs%PGmin  ,'min(PG)'  )
         call this%mfile%add_column(this%fs%TGmax  ,'max(TG)'  )
         call this%mfile%add_column(this%fs%TGmin  ,'min(TG)'  )
         call this%mfile%add_column(this%fs%VFmax  ,'VFmax'    )
         call this%mfile%add_column(this%fs%VFmin  ,'VFmin'    )
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'shockdrop_cfl',restart=monitor_continue)
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLa_x,'Acoustic xCFL')
         call this%cflfile%add_column(this%fs%CFLa_y,'Acoustic yCFL')
         call this%cflfile%add_column(this%fs%CFLa_z,'Acoustic zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         ! Create conservation monitor
         this%consfile=monitor(this%fs%cfg%amRoot,'shockdrop_cons',restart=monitor_continue)
         call this%consfile%add_column(this%time%n,'Timestep number')
         call this%consfile%add_column(this%time%t,'Time')
         call this%consfile%add_column(this%fs%VFint  ,'Volume')
         call this%consfile%add_column(this%fs%Qint(1),'Liquid mass')
         call this%consfile%add_column(this%fs%Qint(2),'Gas mass')
         call this%consfile%add_column(this%fs%Qint(3),'Liquid energy')
         call this%consfile%add_column(this%fs%Qint(4),'Gas energy')
         call this%consfile%add_column(this%fs%Qint(5),'U Momentum')
         call this%consfile%add_column(this%fs%Qint(6),'V Momentum')
         call this%consfile%add_column(this%fs%Qint(7),'W Momentum')
         call this%consfile%add_column(this%fs%RHOKLint,'Liquid KE')
         call this%consfile%add_column(this%fs%RHOKGint,'Gas KE')
         call this%consfile%add_column(this%fs%RHOSLint,'Liquid entropy')
         call this%consfile%add_column(this%fs%RHOSGint,'Gas entropy')
         ! Create drop output
         this%dropfile=monitor(this%fs%cfg%amRoot,'shockdrop_drop',restart=monitor_continue)
         call this%dropfile%add_column(this%time%n,'Timestep number')
         call this%dropfile%add_column(this%time%t,'Time')
         call this%dropfile%add_column(this%fs%VFint  ,'Total volume')
         call this%dropfile%add_column(this%fs%Qint(1),'Total mass')
         call this%dropfile%add_column(this%ccl%nstruct,'N drops')
         call this%dropfile%add_column(this%Vcore,'Core volume')
         call this%dropfile%add_column(this%Mcore,'Core mass')
         call this%dropfile%add_column(this%Xcore,'Core X')
         call this%dropfile%add_column(this%Ycore,'Core Y')
         call this%dropfile%add_column(this%Zcore,'Core Z')
         call this%dropfile%add_column(this%Cmin(1),'Core Xmin')
         call this%dropfile%add_column(this%Cmin(2),'Core Ymin')
         call this%dropfile%add_column(this%Cmin(3),'Core Zmin')
         call this%dropfile%add_column(this%Cmax(1),'Core Xmax')
         call this%dropfile%add_column(this%Cmax(2),'Core Ymax')
         call this%dropfile%add_column(this%Cmax(3),'Core Zmax')
         ! Create mesh output
         this%meshfile=monitor(this%cfg%amRoot,'shockdrop_mesh',restart=monitor_continue)
         call this%meshfile%add_column(this%time%t,'Time')
         call this%meshfile%add_column(this%cfg%nx,'nx')
         call this%meshfile%add_column(this%cfg%ny,'ny')
         call this%meshfile%add_column(this%cfg%nz,'nz')
         call this%meshfile%add_column(this%cfg%x(this%cfg%imin  ),'Xmin')
         call this%meshfile%add_column(this%cfg%y(this%cfg%jmin  ),'Ymin')
         call this%meshfile%add_column(this%cfg%z(this%cfg%kmin  ),'Zmin')
         call this%meshfile%add_column(this%cfg%x(this%cfg%imax+1),'Xmax')
         call this%meshfile%add_column(this%cfg%y(this%cfg%jmax+1),'Ymax')
         call this%meshfile%add_column(this%cfg%z(this%cfg%kmax+1),'Zmax')
      end block create_monitor
      
   end subroutine initialize
   
   
   !> Take one time step
   subroutine step(this,dt)
      implicit none
      class(shockdrop), intent(inout) :: this
      real(WP) :: dt
      
      ! Increment time
      this%time%dt=dt
      call this%fs%get_cfl(dt=this%time%dt,cfl=this%time%cfl)
      call this%time%increment()
      
      ! Remember conserved variables
      this%fs%Qold=this%fs%Q
      
      ! Remember phasic quantities
      this%fs%RHOLold=this%fs%RHOL; this%fs%ILold=this%fs%IL; this%fs%PLold=this%fs%PL
      this%fs%RHOGold=this%fs%RHOG; this%fs%IGold=this%fs%IG; this%fs%PGold=this%fs%PG
      
      ! Remember volume moments and interface
      this%fs%VFold=this%fs%VF
      this%fs%BLold=this%fs%BL
      this%fs%BGold=this%fs%BG
      copy_plic_to_old: block
         use irl_fortran_interface, only: copy
         integer :: i,j,k
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            call copy(this%fs%PLICold(i,j,k),this%fs%PLIC(i,j,k))
         end do; end do; end do
      end block copy_plic_to_old
      
      ! Tag cells for semi-Lagrangian transport
      call this%fs%SLtag()
      
      ! Prepare SGS viscosity models
      call this%prepare_viscosities()
      
      ! Perform first semi-Lagrangian transport step =====================================================
      call this%fs%SLstep(dt=0.5_WP*this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%fs%build_interface()
      
      ! First RK step ====================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,1))
      this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,1)
      ! Increment Q with SL terms
      this%fs%Q=this%fs%Q+this%fs%SLdQ
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Second RK step ===================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,2))
      this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,2)
      ! Increment Q with SL terms
      this%fs%Q=this%fs%Q+this%fs%SLdQ
      ! Apply user-provided relaxation model
      !call this%fs%apply_relax()
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Perform second semi-Lagrangian transport step ====================================================
      call this%fs%SLstep(dt=1.0_WP*this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%fs%build_interface()
      
      ! Third RK step ====================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,3))
      this%fs%Q=this%fs%Qold+1.0_WP*this%time%dt*this%dQdt(:,:,:,:,3)
      ! Increment Q with SL terms
      this%fs%Q=this%fs%Q+this%fs%SLdQ
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Fourth RK step ===================================================================================
      ! Get non-SL RHS and increment
      call this%fs%rhs(this%dQdt(:,:,:,:,4))
      this%fs%Q=this%fs%Qold+this%time%dt/6.0_WP*(this%dQdt(:,:,:,:,1)+2.0_WP*this%dQdt(:,:,:,:,2)+2.0_WP*this%dQdt(:,:,:,:,3)+this%dQdt(:,:,:,:,4))
      ! Increment Q with SL terms
      this%fs%Q=this%fs%Q+this%fs%SLdQ
      ! Apply user-provided relaxation model
      call this%fs%apply_relax()
      ! Recompute primitive variables
      call this%fs%get_primitive()
      
      ! Apply boundary conditions
      !call this%apply_bconds() !< not needed with the multi-domain approach
      
      ! Interpolate velocity
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      
      ! Compute local Mach number
      this%Ma=sqrt(this%Ui**2+this%Vi**2+this%Wi**2)/this%fs%C
      
   end subroutine step
   
   
   !> Perform and output monitoring for the shockdrop problem
   subroutine output_monitor(this)
      implicit none
      class(shockdrop), intent(inout) :: this
      call this%postproc()
      call this%fs%get_info()
      call this%mfile%write()
      call this%cflfile%write()
      call this%consfile%write()
      call this%dropfile%write()
   end subroutine output_monitor
   
   
   !> Output ensight files for the shockdrop problem
   subroutine output_ensight(this,t)
      use irl_fortran_interface, only: getNumberOfVertices
      implicit none
      class(shockdrop), intent(inout) :: this
      real(WP), intent(in), optional :: t
      integer :: i,j,k,np
      ! Update surface mesh
      call this%fs%update_surfmesh(this%smesh)
      ! Update label field
      np=0
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               if (getNumberOfVertices(this%fs%interface_polygon(i,j,k)).gt.0) then
                  np=np+1; this%smesh%var(1,np)=real(this%ccl%id(i,j,k),WP)
               end if
            end do
         end do
      end do
      ! Output to ensight
      if (present(t)) then
         call this%ens_out%write_data(t)
      else
         call this%ens_out%write_data(this%time%t)
      end if
   end subroutine output_ensight
   
   
   !> Calculate viscosities
   subroutine prepare_viscosities(this)
      implicit none
      class(shockdrop), intent(inout) :: this
      real(WP) :: Lvof,Lrho,Gvof,Grho
      real(WP) :: Lvisc,Gvisc,Lbeta,Gbeta
      real(WP), parameter :: eps=1.0e-15_WP
      integer :: i,j,k
      ! Get LAD
      call this%fs%get_viscartif(dt=this%time%dt,beta=this%beta)
      ! Get eddy viscosity
      call this%fs%get_vreman   (dt=this%time%dt,visc=this%visc)
      ! Create mixture viscosities
      do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_-1; do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_-1; do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_-1
         ! Create smooth mass info distribution
         Lvof=sum(       this%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
         Gvof=sum(1.0_WP-this%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
         Lrho=sum(       this%fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,1))/(Lvof+eps)
         Grho=sum(       this%fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,2))/(Gvof+eps)
         ! Harmonic average of VISC
         Lvisc=Lrho*(this%cst_viscL+this%visc(i,j,k)); Gvisc=Grho*(this%cst_viscG+this%visc(i,j,k)); this%fs%VISC(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lvisc,eps)+Gvof/max(Gvisc,eps))
         ! Harmonic average of BETA
         Lbeta=Lrho*this%beta(i,j,k); Gbeta=Grho*this%beta(i,j,k); this%fs%BETA(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lbeta,eps)+Gvof/max(Gbeta,eps))
         ! Try adding BETA to visc
         !this%fs%VISC(i,j,k)=this%fs%VISC(i,j,k)+this%fs%BETA(i,j,k)
      end do; end do; end do
   end subroutine prepare_viscosities
   
   
   !> Apply boundary conditions
   subroutine apply_bconds(this)
      use irl_fortran_interface, only: setPlane
      implicit none
      class(shockdrop), intent(inout) :: this
      integer :: i,j,k
      
      ! Apply UNclipped Neumann on primitive variables in x+
      if (.not.this%fs%cfg%xper.and.this%fs%cfg%iproc.eq.this%fs%cfg%npx) then
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
            ! Copy over from imax to imax+1 and above
            do i=this%fs%cfg%imax+1,this%fs%cfg%imaxo
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(this%fs%cfg%imax,j,k)
               this%fs%PL  (i,j,k)=this%fs%PL  (this%fs%cfg%imax,j,k)
               this%fs%IL  (i,j,k)=this%fs%IL  (this%fs%cfg%imax,j,k)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(this%fs%cfg%imax,j,k)
               this%fs%PG  (i,j,k)=this%fs%PG  (this%fs%cfg%imax,j,k)
               this%fs%IG  (i,j,k)=this%fs%IG  (this%fs%cfg%imax,j,k)
               !this%fs%U  (i,j,k)=max(this%fs%U(this%fs%cfg%imax,j,k),0.0_WP)
               this%fs%U   (i,j,k)=this%fs%U   (this%fs%cfg%imax,j,k)
               this%fs%V   (i,j,k)=this%fs%V   (this%fs%cfg%imax,j,k)
               this%fs%W   (i,j,k)=this%fs%W   (this%fs%cfg%imax,j,k)
               this%fs%VF  (i,j,k)=this%fs%VF  (this%fs%cfg%imax,j,k)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],this%fs%cfg%x(i)+this%fs%dx*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Apply UNclipped Neumann on primitive variables in y+
      if (.not.this%fs%cfg%yper.and.this%fs%cfg%jproc.eq.this%fs%cfg%npy) then
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! Copy over from jmax to jmax+1 and above
            do j=this%fs%cfg%jmax+1,this%fs%cfg%jmaxo
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,this%fs%cfg%jmax,k)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,this%fs%cfg%jmax,k)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,this%fs%cfg%jmax,k)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,this%fs%cfg%jmax,k)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,this%fs%cfg%jmax,k)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,this%fs%cfg%jmax,k)
               this%fs%U   (i,j,k)=this%fs%U   (i,this%fs%cfg%jmax,k)
               !this%fs%V  (i,j,k)=max(this%fs%V(i,this%fs%cfg%jmax,k),0.0_WP)
               this%fs%V   (i,j,k)=this%fs%V   (i,this%fs%cfg%jmax,k)
               this%fs%W   (i,j,k)=this%fs%W   (i,this%fs%cfg%jmax,k)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,this%fs%cfg%jmax,k)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,+1.0_WP,0.0_WP],this%fs%cfg%y(j)+this%fs%dy*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Apply UNclipped Neumann on primitive variables in y-
      if (.not.this%fs%cfg%yper.and.this%fs%cfg%jproc.eq.1) then
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! First copy over V from jmin+1 to jmin
            this%fs%V(i,this%fs%cfg%jmin,k)=min(this%fs%V(i,this%fs%cfg%jmin+1,k),0.0_WP)
            ! Then copy over from jmin to jmin-1 and below
            do j=this%fs%cfg%jmino,this%fs%cfg%jmin-1
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,this%fs%cfg%jmin,k)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,this%fs%cfg%jmin,k)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,this%fs%cfg%jmin,k)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,this%fs%cfg%jmin,k)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,this%fs%cfg%jmin,k)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,this%fs%cfg%jmin,k)
               this%fs%U   (i,j,k)=this%fs%U   (i,this%fs%cfg%jmin,k)
               !this%fs%V  (i,j,k)=min(this%fs%V(i,this%fs%cfg%jmin,k),0.0_WP)
               this%fs%V   (i,j,k)=this%fs%V   (i,this%fs%cfg%jmin,k)
               this%fs%W   (i,j,k)=this%fs%W   (i,this%fs%cfg%jmin,k)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,this%fs%cfg%jmin,k)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,-1.0_WP,0.0_WP],-this%fs%cfg%y(j+1)+this%fs%dy*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Apply UNclipped Neumann on primitive variables in z+
      if (.not.this%fs%cfg%zper.and.this%fs%cfg%kproc.eq.this%fs%cfg%npz) then
         do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! Copy over from kmax to kmax+1 and above
            do k=this%fs%cfg%kmax+1,this%fs%cfg%kmaxo
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,j,this%fs%cfg%kmax)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,j,this%fs%cfg%kmax)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,j,this%fs%cfg%kmax)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,j,this%fs%cfg%kmax)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,j,this%fs%cfg%kmax)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,j,this%fs%cfg%kmax)
               this%fs%U   (i,j,k)=this%fs%U   (i,j,this%fs%cfg%kmax)
               this%fs%V   (i,j,k)=this%fs%V   (i,j,this%fs%cfg%kmax)
               !this%fs%W  (i,j,k)=max(this%fs%W(i,j,this%fs%cfg%kmax),0.0_WP)
               this%fs%W   (i,j,k)=this%fs%W   (i,j,this%fs%cfg%kmax)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,j,this%fs%cfg%kmax)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,+1.0_WP],this%fs%cfg%z(k)+this%fs%dz*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Apply UNclipped Neumann on primitive variables in z-
      if (.not.this%fs%cfg%zper.and.this%fs%cfg%kproc.eq.1) then
         do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            ! First copy over W from kmin+1 to kmin
            this%fs%W(i,j,this%fs%cfg%kmin)=min(this%fs%W(i,j,this%fs%cfg%kmin+1),0.0_WP)
            ! Then copy over from kmin to kmin-1 and below
            do k=this%fs%cfg%kmino,this%fs%cfg%kmin-1
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,j,this%fs%cfg%kmin)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,j,this%fs%cfg%kmin)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,j,this%fs%cfg%kmin)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,j,this%fs%cfg%kmin)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,j,this%fs%cfg%kmin)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,j,this%fs%cfg%kmin)
               this%fs%U   (i,j,k)=this%fs%U   (i,j,this%fs%cfg%kmin)
               this%fs%V   (i,j,k)=this%fs%V   (i,j,this%fs%cfg%kmin)
               !this%fs%W  (i,j,k)=min(this%fs%W(i,j,this%fs%cfg%kmin),0.0_WP)
               this%fs%W   (i,j,k)=this%fs%W   (i,j,this%fs%cfg%kmin)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,j,this%fs%cfg%kmin)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,-1.0_WP],-this%fs%cfg%z(k+1)+this%fs%dz*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Rebuild conserved quantities
      this%fs%Q(:,:,:,1)=        this%fs%VF *this%fs%RHOL
      this%fs%Q(:,:,:,2)=(1.0_WP-this%fs%VF)*this%fs%RHOG
      this%fs%Q(:,:,:,3)= this%fs%Q(:,:,:,1)*this%fs%IL
      this%fs%Q(:,:,:,4)= this%fs%Q(:,:,:,2)*this%fs%IG
      call this%fs%get_momentum()
      
   end subroutine apply_bconds
   
   
   !> Finalize shockdrop problem
   subroutine finalize(this)
      implicit none
      class(shockdrop), intent(inout) :: this
      ! Deallocate work arrays
      deallocate(this%dQdt,this%Ui,this%Vi,this%Wi,this%Ma,this%beta,this%visc)
      ! Finalize all objects
      call this%cfg%finalize()
      call this%fs%finalize()
      call this%time%finalize()
      call this%ccl%finalize()
      call this%smesh%finalize()
      call this%ens_out%finalize()
      call this%mfile%finalize()
      call this%cflfile%finalize()
      call this%consfile%finalize()
      call this%dropfile%finalize()
   end subroutine finalize
   
   
end module shockdrop_class
