!> Definition for a simplex class
module simplex_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use polygon_class,     only: polygon
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use cclabel_class,     only: cclabel
   use iterator_class,    only: iterator
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   implicit none
   private
   
   public :: simplex
   
   !> Simplex object
   type :: simplex
      
      !> Provide a pardata and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata)  :: df
      logical :: restarted
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB based on polygon
      type(polygon)  :: poly
      type(ibconfig) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< HYPRE linear solver for pressure
      type(ddadi)       :: vs    !< DDADI linear solver for velocity
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info
      type(cclabel)     :: ccl   !< CCLabel to transfer droplets
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight) :: ens_out   !< Ensight output for flow variables
      type(event)   :: ens_evt   !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      !real(WP), dimension(:,:,:,:), allocatable :: SR                !< Strain rate
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      integer :: nlayer=4                  !< Size of buffer layer for VOF removal
      
      !> Timing info
      type(monitor) :: timefile !< Timing monitoring
      type(timer)   :: tstep    !< Timer for step
      type(timer)   :: tsgs     !< Timer for SGS
      type(timer)   :: tvel     !< Timer for velocity
      type(timer)   :: tpres    !< Timer for pressure
      type(timer)   :: tvof     !< Timer for VOF
      type(timer)   :: ttrans   !< Timer for VOF transfer

      !> Event for flow rate analysis
      type(event) :: flowrate_evt  !< Event trigger for flow rate analysis
      
      !> Inlet pipes geometry and flow rates
      real(WP) :: Rinlet=0.002_WP
      real(WP) :: Rexit=0.00143_WP
      real(WP) :: Rpipe=0.000185_WP
      real(WP) :: Rcoflow=0.003_WP
      real(WP), dimension(3) :: p1=[-0.00442_WP,0.0_WP,+0.001245_WP]
      real(WP), dimension(3) :: p2=[-0.00442_WP,0.0_WP,-0.001245_WP]
      real(WP), dimension(3) :: n1=[+0.6_WP,-0.8_WP,0.0_WP]
      real(WP), dimension(3) :: n2=[+0.6_WP,+0.8_WP,0.0_WP]
      real(WP) :: Ucoflow,mfr,Apipe
      
   contains
      procedure :: init                            !< Initialize simplex simulation
      procedure :: step                            !< Advance simplex simulation by one time step
      procedure :: final                           !< Finalize simplex simulation
      procedure :: transfer_drops                  !< Transfer drops to a Lagrangian representation
      procedure :: analyze_flowrate                !< Compute and output flow rate through the nozzle
   end type simplex
   
   
contains
   
   
   !> Compute and output flow rate through the nozzle
   subroutine analyze_flowrate(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      use filesys,  only: makedir,isdir
      use string,   only: str_medium
      implicit none
      class(simplex), intent(inout) :: this
      real(WP), dimension(:), allocatable :: CSA_s,CSA_f,CSA_l,CSA_g !< Solid, fluid, liquid, and gas cross-sectional areas
      real(WP), dimension(:), allocatable :: VFR_s,VFR_f,VFR_l,VFR_g !< Solid, fluid, liquid, and gas volume flow rates
      character(len=str_medium) :: filename,timestamp
      integer :: i,j,k,ierr,iunit
      real(WP) :: area,rad
      ! Allocate and initialize 1D storage for areas and flow rates
      allocate(CSA_s(this%cfg%imin:this%cfg%imax)); CSA_s=0.0_WP
      allocate(CSA_f(this%cfg%imin:this%cfg%imax)); CSA_f=0.0_WP
      allocate(CSA_l(this%cfg%imin:this%cfg%imax)); CSA_l=0.0_WP
      allocate(CSA_g(this%cfg%imin:this%cfg%imax)); CSA_g=0.0_WP
      allocate(VFR_s(this%cfg%imin:this%cfg%imax)); VFR_s=0.0_WP
      allocate(VFR_f(this%cfg%imin:this%cfg%imax)); VFR_f=0.0_WP
      allocate(VFR_l(this%cfg%imin:this%cfg%imax)); VFR_l=0.0_WP
      allocate(VFR_g(this%cfg%imin:this%cfg%imax)); VFR_g=0.0_WP
      ! Integrate various areas and flow rates
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Cell cross-sectional area and radial location
               area=this%cfg%dy(j)*this%cfg%dz(k)
               rad=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
               ! Increment solid area and flow rate
               CSA_s(i)=CSA_s(i)+area*(1.0_WP-this%cfg%VF(i,j,k))
               VFR_s(i)=VFR_s(i)+area*(1.0_WP-this%cfg%VF(i,j,k))*this%Ui(i,j,k)
               ! Increment fluid areas and flow rates only inside the nozzle
               if ((this%cfg%xm(i).lt.-0.0015_WP.and.rad.le.this%Rinlet).or.(this%cfg%xm(i).ge.-0.0015_WP.and.this%cfg%xm(i).lt.0.0_WP.and.rad.le.this%Rexit).or.this%cfg%xm(i).gt.0.0_WP) then
                  CSA_f(i)=CSA_f(i)+area*this%cfg%VF(i,j,k)
                  CSA_l(i)=CSA_l(i)+area*this%cfg%VF(i,j,k)*        this%vf%VF(i,j,k)
                  CSA_g(i)=CSA_g(i)+area*this%cfg%VF(i,j,k)*(1.0_WP-this%vf%VF(i,j,k))
                  VFR_f(i)=VFR_f(i)+area*this%cfg%VF(i,j,k)                           *this%Ui(i,j,k)
                  VFR_l(i)=VFR_l(i)+area*this%cfg%VF(i,j,k)*        this%vf%VF(i,j,k) *this%Ui(i,j,k)
                  VFR_g(i)=VFR_g(i)+area*this%cfg%VF(i,j,k)*(1.0_WP-this%vf%VF(i,j,k))*this%Ui(i,j,k)
               end if
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_s,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_f,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_l,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_g,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_s,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_f,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_l,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_g,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      ! Only root process outputs to a file
      if (this%cfg%amRoot) then
         if (.not.isdir('flowrate')) call makedir('flowrate')
         filename='flowrate_'; write(timestamp,'(es12.5)') this%time%t
         open(newunit=iunit,file='flowrate/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(999999(a12,x))') 'xm','CSA_s','CSA_f','CSA_l','CSA_g','VFR_s','VFR_f','VFR_l','VFR_g'
         do i=this%cfg%imin,this%cfg%imax
            write(iunit,'(999999(es12.5,x))') this%cfg%xm(i),CSA_s(i),CSA_f(i),CSA_l(i),CSA_g(i),VFR_s(i),VFR_f(i),VFR_l(i),VFR_g(i)
         end do
         close(iunit)
      end if
   end subroutine analyze_flowrate
   
   
   !> Transfer droplet to Lagrangian representation
   subroutine transfer_drops(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      class(simplex), intent(inout) :: this
      real(WP), dimension(:), allocatable :: dvol
      integer :: n,m,ierr,nmax
      
      ! Start by performing a CCL
      call this%ccl%build(make_label,same_label)
      
      ! Allocate droplet stats arrays
      !allocate(dvol(1:this%ccl%nstruct)); dvol=0.0_WP
      
      ! First pass: loop over individual structures and accumulate stats
      !do n=1,this%ccl%nstruct
      !   ! Loop over cells in structure
      !   do m=1,this%ccl%struct(n)%n_
      !      ! Accumulate volume
      !      dvol(n)=dvol(n)+this%cfg%vol(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))*&
      !      &                 this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))
      !      ! Accumulate ...
      !   end do
      !end do
      
      ! Reduce stats
      !call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      
      
      ! Find the liquid core
      !nmax=maxloc(dvol,dim=1)
      
      ! Second pass to transfer drops
      !do n=1,this%ccl%nstruct
      !   ! Skip structures that do not meet our criteria
      !   if (n.eq.nmax) cycle
      !   ! Remove all other structures
      !   do m=1,this%ccl%struct(n)%n_
      !      this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
      !   end do
      !   ! Create a droplet instead
      !
      !end do
      !call this%vf%sync_interface()
      !call this%vf%clean_irl_and_band()
      
      ! Deallocate
      !deallocate(dvol)
      
   contains
      
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.0_WP) then
            make_label=.true.
         else
            make_label=.false.
         end if
      end function make_label
      
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
      
   end subroutine transfer_drops
   
   
   !> Initialization of simplex simulation
   subroutine init(this)
      implicit none
      class(simplex), intent(inout) :: this
      
      
      ! Setup an input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='simplex.input')
      end block read_input
      
      
      ! Initialize ibconfig object
      create_config: block
         use parallel,    only: group
         use sgrid_class, only: cartesian,sgrid
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,xshift,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni,x,y,z
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x_uni(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y_uni(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z_uni(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! Add stretching
         call this%input%read('Stretched cells in yz',ns_yz,default=0)
         if (ns_yz.gt.0) call this%input%read('Stretch ratio in yz',sratio_yz)
         call this%input%read('Stretched cells in x' ,ns_x ,default=0)
         if (ns_x .gt.0) call this%input%read('Stretch ratio in x' ,sratio_x )
         allocate(x(nx+1+1*ns_x )); x(      1:      1+nx)=x_uni
         allocate(y(ny+1+2*ns_yz)); y(ns_yz+1:ns_yz+1+ny)=y_uni
         allocate(z(nz+1+2*ns_yz)); z(ns_yz+1:ns_yz+1+nz)=z_uni
         do i=nx+2,nx+1+ns_x
            x(i)=x(i-1)+sratio_x *(x(i-1)-x(i-2))
         end do
         do j=ns_yz,1,-1
            y(j)=y(j+1)+sratio_yz*(y(j+1)-y(j+2))
         end do
         do j=ns_yz+2+ny,ny+1+2*ns_yz
            y(j)=y(j-1)+sratio_yz*(y(j-1)-y(j-2))
         end do
         do k=ns_yz,1,-1
            z(k)=z(k+1)+sratio_yz*(z(k+1)-z(k+2))
         end do
         do k=ns_yz+2+nz,nz+1+2*ns_yz
            z(k)=z(k-1)+sratio_yz*(z(k-1)-z(k-2))
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='simplex')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create ibconfig
         this%cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      
      ! Now initialize simplex nozzle geometry
      create_simplex: block
         use ibconfig_class, only: sharp
         integer :: i,j,k
         ! Create polygon
         call this%poly%initialize(nvert=10,name='simplex')
         this%poly%vert(:, 1)=[-0.01000_WP,0.00000_WP]
         this%poly%vert(:, 2)=[-0.00442_WP,0.00000_WP]
         this%poly%vert(:, 3)=[-0.00442_WP,0.00160_WP]
         this%poly%vert(:, 4)=[-0.00385_WP,0.00160_WP]
         this%poly%vert(:, 5)=[-0.00175_WP,0.00039_WP]
         this%poly%vert(:, 6)=[-0.00114_WP,0.00039_WP]
         this%poly%vert(:, 7)=[ 0.00000_WP,0.00143_WP]
         this%poly%vert(:, 8)=[ 0.00000_WP,0.00177_WP]
         this%poly%vert(:, 9)=[-0.00122_WP,0.00279_WP]
         this%poly%vert(:,10)=[-0.01000_WP,0.00279_WP]
         ! Initialize IB distance field
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Calculate distance from object obtained by revolution of polygon
                  this%cfg%Gib(i,j,k)=-this%poly%get_distance([this%cfg%xm(i),sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)])
               end do
            end do
         end do
         ! Get normal vector
         call this%cfg%calculate_normal()
         ! Get VF field
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
         ! Carve out inlet pipes
         create_inlet_pipes: block
            use mms_geom, only: cube_refine_vol
            integer :: si,sj,sk,n
            real(WP), dimension(3,8) :: cube_vertex
            real(WP), dimension(3) :: v_cent,a_cent
            real(WP) :: vol,area,contact
            integer, parameter :: amr_ref_lvl=4
            do k=this%cfg%kmino_,this%cfg%kmaxo_
               do j=this%cfg%jmino_,this%cfg%jmaxo_
                  do i=this%cfg%imino_,this%cfg%imaxo_
                     ! Only work to the left of the plenum
                     if (this%cfg%xm(i)-this%p1(1).gt.this%cfg%min_meshsize) cycle
                     ! Set cube vertices
                     n=0
                     do sk=0,1
                        do sj=0,1
                           do si=0,1
                              n=n+1; cube_vertex(:,n)=[this%cfg%x(i+si),this%cfg%y(j+sj),this%cfg%z(k+sk)]
                           end do
                        end do
                     end do
                     ! Call adaptive refinement code to get volume fraction recursively - inlet pipe 1
                     vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_inlet_pipe_1,0.0_WP,amr_ref_lvl)
                     this%cfg%VF(i,j,k)=max(this%cfg%VF(i,j,k),vol/(this%cfg%dx(i)*this%cfg%dy(j)*this%cfg%dz(k)))
                     ! Call adaptive refinement code to get volume fraction recursively - inlet pipe 2
                     vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_inlet_pipe_2,0.0_WP,amr_ref_lvl)
                     this%cfg%VF(i,j,k)=max(this%cfg%VF(i,j,k),vol/(this%cfg%dx(i)*this%cfg%dy(j)*this%cfg%dz(k)))
                  end do
               end do
            end do
         end block create_inlet_pipes
         ! Apply Neumann on VF and apply stair-stepping at entrance
         if (this%cfg%iproc.eq.1) then
            ! Stair-step entrance
            this%cfg%VF(this%cfg%imin,:,:)=max(real(nint(this%cfg%VF(this%cfg%imin,:,:)),WP),epsilon(1.0_WP))
            ! Copy into overlap layer
            do i=this%cfg%imino,this%cfg%imin-1
               this%cfg%VF(i,:,:)=this%cfg%VF(this%cfg%imin,:,:)
            end do
         end if
         ! Recompute domain volume
         call this%cfg%calc_fluid_vol()
      end block create_simplex
      
      
      ! Initialize flow rate
      set_flowrate: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         use string,   only: str_long
         use messager, only: log
         integer :: j,k,ierr
         character(len=str_long) :: message
         ! Read mass flow rate
         call this%input%read('Mass flow rate',this%mfr)
         ! Read coflow velocity
         call this%input%read('Coflow velocity',this%Ucoflow)
         ! Integrate inlet pipe surface area
         this%Apipe=0.0_WP
         if (this%cfg%iproc.eq.1) then
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).lt.this%Rinlet) then
                     this%Apipe=this%Apipe+this%cfg%VF(this%cfg%imin,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                  end if
               end do
            end do
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Apipe,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! Log this info
         if (this%cfg%amRoot) then
            write(message,'("Inlet pipe area is ",es12.5)') this%Apipe; call log(message)
         end if
      end block set_flowrate
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Handle restart/saves here
      restart_and_save: block
         use string,  only: str_medium
         use filesys, only: makedir,isdir
         character(len=str_medium) :: filename
         integer, dimension(3) :: iopartition
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',filename,default='')
         this%restarted=.false.; if (len_trim(filename).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! We are restarting, read the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/'//trim(filename))
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block restart_and_save
      
      
      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
         if (this%restarted) then
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
         end if
      end block update_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !allocate(this%SR(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: remap,plicnet,r2p,r2pnet,lvira
         use irl_fortran_interface
         integer :: i,j,k
         real(WP) :: rad
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         ! Create a VOF solver with plicnet
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=remap,name='VOF')
         !this%vf%twoplane_thld2=0.3_WP
         !this%vf%thin_thld_min=1.0e-3_WP
         ! Initialize the interface inclduing restarts
         if (this%restarted) then
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! Check if the second plane is meaningful
                     if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                     else
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                     end if
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Ensure that boundaries are correct
            if (this%vf%cfg%iproc.eq.1)               this%vf%VF(this%vf%cfg%imino:this%vf%cfg%imin-1,:,:)=0.0_WP
            if (this%vf%cfg%iproc.eq.this%vf%cfg%npx) this%vf%VF(this%vf%cfg%imax+1:this%vf%cfg%imaxo,:,:)=0.0_WP
            if (this%vf%cfg%jproc.eq.1)               this%vf%VF(:,this%vf%cfg%jmino:this%vf%cfg%jmin-1,:)=0.0_WP
            if (this%vf%cfg%jproc.eq.this%vf%cfg%npy) this%vf%VF(:,this%vf%cfg%jmax+1:this%vf%cfg%jmaxo,:)=0.0_WP
            if (this%vf%cfg%kproc.eq.1)               this%vf%VF(:,:,this%vf%cfg%kmino:this%vf%cfg%kmin-1)=0.0_WP
            if (this%vf%cfg%kproc.eq.this%vf%cfg%npz) this%vf%VF(:,:,this%vf%cfg%kmax+1:this%vf%cfg%kmaxo)=0.0_WP
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                     if (i.lt.this%vf%cfg%imin.and.rad.le.this%Rinlet) this%vf%VF(i,j,k)=1.0_WP
                  end do
               end do
            end do
            ! Update the band
            call this%vf%update_band()
            ! Set interface planes at the boundaries
            call this%vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
         else
            ! Initialize to flat interface at exit
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                     ! Ensure the nozzle is filled with liquid up to the throat with wet walls
                     if (this%vf%cfg%xm(i).lt.-0.0015_WP.and.rad.le.this%Rinlet) then
                        this%vf%VF(i,j,k)=1.0_WP
                     else if (this%vf%cfg%xm(i).ge.-0.0015_WP.and.this%vf%cfg%xm(i).lt.0.0_WP.and.rad.le.this%Rexit) then
                        this%vf%VF(i,j,k)=1.0_WP
                     else
                        this%vf%VF(i,j,k)=0.0_WP
                     end if
                     ! Initialize phasic barycenters
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end do
               end do
            end do
            ! Update the band
            call this%vf%update_band()
            ! Perform interface reconstruction from VOF field
            call this%vf%build_interface()
            ! Set interface planes at the boundaries
            call this%vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Reset moments to guarantee compatibility with interface reconstruction
            call this%vf%reset_volume_moments()
         end if
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create a two-phase flow solver with bconds
      create_flow_solver: block
         use tpns_class,      only: clipped_neumann,dirichlet,slip
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-Phase NS')
         ! Set the flow properties
         call this%input%read('Liquid dynamic viscosity',this%fs%visc_l)
         call this%input%read('Gas dynamic viscosity'   ,this%fs%visc_g)
         call this%input%read('Liquid density',this%fs%rho_l)
         call this%input%read('Gas density'   ,this%fs%rho_g)
         call this%input%read('Surface tension coefficient',this%fs%sigma)
         ! Set acceleration of gravity
         call this%input%read('Gravity',this%fs%gravity)
         ! Inlets on the left
         call this%fs%add_bcond(name='inlets',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=pipe_inlets)
         call this%fs%add_bcond(name='coflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=coflow_inlet)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         ! Zero velocity except if restarting
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         if (this%restarted) then
            ! Read data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
         end if
         ! Apply Dirichlet condition at pipe inlets
         call this%fs%get_bcond('inlets',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=this%cfg%VF(i,j,k)*this%mfr/(this%fs%rho_l*this%Apipe)
         end do
         ! Apply Dirichlet condition at coflow inlet
         call this%fs%get_bcond('coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=this%Ucoflow
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute divergence
         call this%fs%get_div()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      end block initialize_velocity
      
      
      ! Create CCL
      create_ccl: block
         call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      end block create_ccl
      
      
      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         this%smesh=surfmesh(nvar=0,name='plic')
         call this%vf%update_surfmesh_nowall(this%smesh)
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='simplex')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('divergence',this%fs%div)
         call this%ens_out%add_surface('plic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_simplex')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_simplex')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
      
      ! Create a timing monitor
      create_timing: block
         ! Create timers
         this%tstep =timer(comm=this%cfg%comm,name='Timestep')
         this%tvof  =timer(comm=this%cfg%comm,name='VOFsolve')
         this%tvel  =timer(comm=this%cfg%comm,name='Velocity')
         this%tpres =timer(comm=this%cfg%comm,name='Pressure')
         this%tsgs  =timer(comm=this%cfg%comm,name='SGSmodel')
         this%ttrans=timer(comm=this%cfg%comm,name='Transfer')
         ! Create corresponding monitor file
         this%timefile=monitor(this%fs%cfg%amRoot,'timing')
         call this%timefile%add_column(this%time%n,'Timestep number')
         call this%timefile%add_column(this%time%t,'Time')
         call this%timefile%add_column(this%tstep%time ,trim(this%tstep%name))
         call this%timefile%add_column(this%tvof%time  ,trim(this%tvof%name))
         call this%timefile%add_column(this%tvel%time  ,trim(this%tvel%name))
         call this%timefile%add_column(this%tpres%time ,trim(this%tpres%name))
         call this%timefile%add_column(this%tsgs%time  ,trim(this%tsgs%name))
         call this%timefile%add_column(this%ttrans%time,trim(this%ttrans%name))
      end block create_timing
      
      
      ! Create an event for flow rate analysis
      flowrate_analysis_prep: block
         this%flowrate_evt=event(time=this%time,name='Flow rate output')
         call this%input%read('Flow rate output period',this%flowrate_evt%tper,default=huge(1.0_WP))
      end block flowrate_analysis_prep
      
      
   contains
      
      !> Function that localizes the right domain boundary
      function right_boundary(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function right_boundary
      
      
      !> Function that localizes the pipe inlets
      function pipe_inlets(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin.and.sqrt(pg%ym(j)**2+pg%zm(k)**2).lt.this%Rinlet) isIn=.true.
      end function pipe_inlets
      
      
      !> Function that localizes the coflow inlet
      function coflow_inlet(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin.and.sqrt(pg%ym(j)**2+pg%zm(k)**2).gt.this%Rcoflow) isIn=.true.
      end function coflow_inlet
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer.or.&
         &   j.le.pg%jmin+this%nlayer.or.&
         &   j.ge.pg%jmax-this%nlayer.or.&
         &   k.le.pg%kmin+this%nlayer.or.&
         &   k.ge.pg%kmax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      
      !> Function that localizes the top (z+) of the domain
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      
      !> Function that localizes the bottom (z-) of the domain
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      
      !> Function that defines a level set function for inlet pipe 1
      function levelset_inlet_pipe_1(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP), dimension(3) :: v,p
         real(WP) :: G
         v=xyz-this%p1
         p=v-this%n1*dot_product(v,this%n1)
         G=this%Rpipe-sqrt(dot_product(p,p))
      end function levelset_inlet_pipe_1
      
      
      !> Function that defines a level set function for inlet pipe 2
      function levelset_inlet_pipe_2(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP), dimension(3) :: v,p
         real(WP) :: G
         v=xyz-this%p2
         p=v-this%n2*dot_product(v,this%n2)
         G=this%Rpipe-sqrt(dot_product(p,p))
      end function levelset_inlet_pipe_2
      
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(simplex), intent(inout) :: this
      
      ! Reset all timers and start timestep timer
      call this%tstep%reset()
      call this%tvof%reset()
      call this%tsgs%reset()
      call this%tvel%reset()
      call this%tpres%reset()
      call this%ttrans%reset()
      call this%tstep%start()
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
      
      ! VOF solver step
      call this%tvof%start() !< Start VOF timer
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%tvof%stop() !< Stop VOF timer
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
      
      ! Turbulence modeling
      call this%tsgs%start() !< Start SGS timer
      sgs_modeling: block
         use sgsmodel_class, only: vreman,dynamic_smag
         integer :: i,j,k
         this%resU=this%vf%VF*this%fs%rho_l+(1.0_WP-this%vf%VF)*this%fs%rho_g
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         !call this%fs%get_strainrate(this%SR)
         !call this%sgs%get_visc(type=dynamic_smag,dt=this%time%dtold,rho=this%resU,Ui=this%Ui,Vi=this%Vi,Wi=this%Wi,SR=this%SR)
         do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
            this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
            this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
            this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
            this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
         end do; end do; end do
      end block sgs_modeling
      call this%tsgs%stop() !< Stop SGS timer
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         !> Start velocity timer
         call this%tvel%start()
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Add momentum source terms
         call this%fs%addsrc_gravity(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply IB forcing to enforce BC at the pipe walls
         ibforcing: block
            integer :: i,j,k
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     if (this%fs%umask(i,j,k).eq.0) this%fs%U(i,j,k)=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*this%fs%U(i,j,k)
                     if (this%fs%vmask(i,j,k).eq.0) this%fs%V(i,j,k)=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*this%fs%V(i,j,k)
                     if (this%fs%wmask(i,j,k).eq.0) this%fs%W(i,j,k)=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*this%fs%W(i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibforcing
        
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         !> Stop velocity timer and start pressure timer
         call this%tvel%stop()
         call this%tpres%start()
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         !call this%fs%add_surface_tension_jump_thin(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         !call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         !> Stop pressure timer
         call this%tpres%stop()
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Transfer VOF into droplets
      call this%ttrans%start() !< Start transfer timer
      call this%transfer_drops()
      call this%ttrans%stop() !< Stop transfer timer
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            if (n.le.this%vof_removal_layer%n_) this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
      end block remove_vof
      
      ! Reset VOF in IBs
      reset_VOF_IB: block
         use ibconfig_class, only: VFhi
         integer :: i,j,k
         real(WP) :: rad
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Only work in wall/IB cells
                  if (this%vf%cfg%VF(i,j,k).gt.VFhi) cycle
                  ! Compute our local radius
                  rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                  ! Ensure the nozzle is filled with liquid up to the throat with wet walls
                  if (this%vf%cfg%xm(i).lt.-0.0015_WP.and.rad.le.this%Rinlet) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else if (this%vf%cfg%xm(i).ge.-0.0015_WP.and.this%vf%cfg%xm(i).lt.0.0_WP.and.rad.le.this%Rexit) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else
                     this%vf%VF(i,j,k)=0.0_WP
                  end if
                  ! Initialize phasic barycenters
                  this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
               end do
            end do
         end do
         call this%vf%clean_irl_and_band()
      end block reset_VOF_IB
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         call this%vf%update_surfmesh_nowall(this%smesh)
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Output flow rate
      if (this%flowrate_evt%occurs()) call this%analyze_flowrate()
      
      !> Stop timestep timer
      call this%tstep%stop()
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%timefile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! First plane
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     ! Second plane
                     plane=0.0_WP
                     if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                     P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
         end block save_restart
      end if
      
   end subroutine step
   
   
   !> Finalize simplex simulation
   subroutine final(this)
      implicit none
      class(simplex), intent(inout) :: this
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)!,this%SR)
   end subroutine final
   
   
end module simplex_class