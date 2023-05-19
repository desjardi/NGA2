!> Definition for a nozzle class
module nozzle_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use fft2d_class,       only: fft2d
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: nozzle
   
   !> Nozzle object
   type :: nozzle
      
      !> Provide a datafile and an event tracker for saving restarts
      type(event)    :: save_evt
      type(datafile) :: df
      logical :: restarted
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB
      type(ibconfig) :: cfg
      
      !> Surface mesh for IB
      type(surfmesh) :: plymesh
      
      !> Flow solver
      type(incomp)      :: fs    !< Incompressible flow solver
      type(fft2d)       :: ps    !< FFT-accelerated linear solver for pressure
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info
      
      !> Ensight postprocessing
      type(ensight) :: ens_out  !< Ensight output for flow variables
      type(event)   :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> IB velocity and mass source
      real(WP), dimension(:,:,:), allocatable :: Uib,Vib,Wib,srcM
      
      !> Fluid definition
      real(WP) :: visc
      
   contains
      procedure, private :: geometry_init          !< Initialize geometry for nozzle
      procedure, private :: simulation_init        !< Initialize simulation for nozzle
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type nozzle
   

   !> Geometric data for case
   real(WP), parameter :: Rgas1=1.5e-3_WP
   real(WP), parameter :: Rgas2=5.0e-3_WP
   real(WP), parameter :: Xaxial=-0.06394704_WP
   real(WP), parameter :: Daxial=+0.01600000_WP
   real(WP), parameter :: Xswirl=-0.06394704_WP
   real(WP), parameter :: Dswirl=+0.00711200_WP
   real(WP), parameter :: Lswirl=+0.02945130_WP
   real(WP), parameter :: sidewall=0.041_WP         !< For abs(y) and abs(z)>sidewall, set to wall
   real(WP), parameter :: nspread=2.0_WP            !< Spread IB velocity over a few cells (is it needed?)
   
   
contains
   
   
   !> Initialization of nozzle simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(nozzle), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_nozzle')
      
      ! Initialize the geometry
      call this%geometry_init()
      
      ! Initialize the simulation
      call this%simulation_init()
      
   end subroutine init
   
   
   !> Initialize geometry
   subroutine geometry_init(this)
      use sgrid_class, only: sgrid
      implicit none
      class(nozzle) :: this
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
         
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
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='nozzle')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call this%input%read('Partition',partition)
         
         ! Create partitioned grid
         this%cfg=ibconfig(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Read in the PLY geometry
      read_ply: block
         use string,   only: str_medium
         use parallel, only: MPI_REAL_WP
         use mpi_f08
         character(len=str_medium) :: plyfile
         integer :: ierr,size_conn
         
         ! Read in ply filename
         call this%input%read('PLY filename',plyfile)
         
         ! Root creates surface mesh from ply, other process create empty surface mesh
         if (this%cfg%amRoot) then
            this%plymesh=surfmesh(plyfile=plyfile,nvar=0,name='ply')
         else
            this%plymesh=surfmesh(nvar=0,name='ply')
         end if
         
         ! Go through parallel broadcast of surface mesh
         call MPI_BCAST(this%plymesh%nVert   ,1                 ,MPI_INTEGER,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%nPoly   ,1                 ,MPI_INTEGER,0,this%cfg%comm,ierr)
         if (.not.this%cfg%amRoot) call this%plymesh%set_size(nvert=this%plymesh%nVert,npoly=this%plymesh%nPoly)
         call MPI_BCAST(this%plymesh%xVert   ,this%plymesh%nVert,MPI_REAL_WP,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%yVert   ,this%plymesh%nVert,MPI_REAL_WP,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%zVert   ,this%plymesh%nVert,MPI_REAL_WP,0,this%cfg%comm,ierr)
         call MPI_BCAST(this%plymesh%polySize,this%plymesh%nPoly,MPI_INTEGER,0,this%cfg%comm,ierr)
         if (this%cfg%amRoot) size_conn=size(this%plymesh%polyConn)
         call MPI_BCAST(size_conn            ,1                 ,MPI_INTEGER,0,this%cfg%comm,ierr)
         if (.not.this%cfg%amRoot) allocate(this%plymesh%polyConn(size_conn))
         call MPI_BCAST(this%plymesh%polyConn,size_conn         ,MPI_INTEGER,0,this%cfg%comm,ierr)
         
      end block read_ply
      
      
      ! Create IB walls for this config
      create_walls: block
         use ibconfig_class, only: sharp
         use messager,       only: die
         use mathtools,      only: cross_product,normalize
         use irl_fortran_interface
         real(WP), parameter :: safe_coeff=3.0_WP
         integer :: i,j,k,np,nv,iv,ip
         real(WP) :: mydist
         real(WP), dimension(3) :: pos,nearest_pt,mynearest,mynorm
         type(Poly_type), dimension(:), allocatable :: poly
         real(WP), dimension(:,:), allocatable :: bary,vert,nvec
         
         ! Preprocess surface mesh data using IRL
         allocate(poly(1:this%plymesh%nPoly))
         allocate(vert(1:3,1:maxval(this%plymesh%polySize)))
         allocate(bary(1:3,1:this%plymesh%nPoly))
         allocate(nvec(1:3,1:this%plymesh%nPoly))
         do np=1,this%plymesh%nPoly
            ! Allocate polygon
            call new(poly(np))
            ! Fill it up
            do nv=1,this%plymesh%polySize(np)
               iv=sum(this%plymesh%polySize(1:np-1))+nv
               vert(:,nv)=[this%plymesh%xVert(this%plymesh%polyConn(iv)),this%plymesh%yVert(this%plymesh%polyConn(iv)),this%plymesh%zVert(this%plymesh%polyConn(iv))]
            end do
            call construct(poly(np),this%plymesh%polySize(np),vert(1:3,1:this%plymesh%polySize(np)))
            mynorm=normalize(cross_product(vert(:,2)-vert(:,1),vert(:,3)-vert(:,2)))
            mydist=dot_product(mynorm,vert(:,1))
            call setPlaneOfExistence(poly(np),[mynorm(1),mynorm(2),mynorm(3),mydist])
            ! Also store its barycenter
            bary(:,np)=calculateCentroid(poly(np))
            nvec(:,np)=calculateNormal  (poly(np))
         end do
         deallocate(vert)
         
         ! Create IB distance field using IRL
         this%cfg%Gib=huge(1.0_WP)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Store cell center position
                  pos=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
                  ! Traverse all polygons
                  do np=1,this%plymesh%nPoly
                     ! Calculate distance to centroid
                     nearest_pt=pos-bary(:,np)
                     mydist=dot_product(nearest_pt,nearest_pt)
                     ! If close enough, compute exact distance to the polygon instead
                     if (mydist.lt.(safe_coeff*this%cfg%min_meshsize)**2) then
                        nearest_pt=calculateNearestPtOnSurface(poly(np),pos)
                        nearest_pt=pos-nearest_pt
                        mydist=dot_product(nearest_pt,nearest_pt)
                     end if
                     ! Remember closest distance
                     if (mydist.lt.this%cfg%Gib(i,j,k)) then
                        this%cfg%Gib(i,j,k)=mydist
                        mynearest=nearest_pt
                        ip=np
                     end if
                  end do
                  ! Take the square root
                  this%cfg%Gib(i,j,k)=sqrt(this%cfg%Gib(i,j,k))
                  ! Find the sign
                  if (dot_product(mynearest,nvec(:,ip)).gt.0.0_WP) this%cfg%Gib(i,j,k)=-this%cfg%Gib(i,j,k)
               end do
            end do
         end do
         deallocate(bary,nvec,poly)
         
         ! Get normal vector
         call this%cfg%calculate_normal()
         
         ! Get VF field
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
         
         ! Add walls all around in y and z to allow periodicity
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%cfg%xm(i).lt.0.0_WP.and.abs(this%cfg%ym(j)).gt.sidewall.or.abs(this%cfg%zm(k)).gt.sidewall) this%cfg%VF(i,j,k)=epsilon(1.0_WP)
               end do
            end do
         end do
         call this%cfg%calc_fluid_vol()

         ! Redraw the x+ walls as stair-stepped to ensure conservation
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (i.ge.this%cfg%imax) then
                     mydist=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
                     if (mydist.lt.Rgas1.or.mydist.gt.Rgas2) then
                        this%cfg%VF(i,j,k)=epsilon(1.0_WP)
                     else
                        this%cfg%VF(i,j,k)=1.0_WP
                     end if
                  end if
               end do
            end do
         end do
         call this%cfg%calc_fluid_vol()
         
      end block create_walls
      
      
   end subroutine geometry_init
   
   
   !> Initialize simulation
   subroutine simulation_init(this)
      implicit none
      class(nozzle), intent(inout) :: this
      
      
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
         use string, only: str_medium
         character(len=str_medium) :: timestamp
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         if (this%restarted) then
            ! Read the datafile
            this%df=datafile(pg=this%cfg,fdata='restart/data_'//trim(adjustl(timestamp)))
         else
            ! Prepare a new directory for storing files for restart
            if (this%cfg%amRoot) call execute_command_line('mkdir -p restart')
            ! If we are not restarting, we will still need a datafile for saving restart files
            this%df=datafile(pg=this%cfg,filename=trim(this%cfg%name),nval=2,nvar=4)
            this%df%valname(1)='t'
            this%df%valname(2)='dt'
            this%df%varname(1)='U'
            this%df%varname(2)='V'
            this%df%varname(3)='W'
            this%df%varname(4)='P'
         end if
      end block restart_and_save
      
      
      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
         if (this%restarted) then
            call this%df%pullval(name='t' ,val=this%time%t )
            call this%df%pullval(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
         end if
      end block update_timetracker


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Uib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%srcM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use incomp_class, only: clipped_neumann
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call this%input%read('Density',this%fs%rho)
         call this%input%read('Dynamic viscosity',this%visc); this%fs%visc=this%visc
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary)
         ! Configure pressure solver
         this%ps=fft2d(cfg=this%cfg,name='Pressure',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)
      end block create_flow_solver
      
      
      ! Initialize our IB velocity field
      set_ib_velocity: block
         use mathtools, only: Pi
         integer  :: i,j,k
         real(WP) :: Qaxial,Uaxial,Aaxial
         real(WP) :: Qswirl,Uswirl,Aswirl
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         this%Uib=0.0_WP; this%Vib=0.0_WP; this%Wib=0.0_WP
         ! Read in axial gas flow rate, convert to SI, get velocity
         call this%input%read('Axial flow rate (SLPM)',Qaxial)
         Qaxial=Qaxial*SLPM2SI
         Aaxial=0.25_WP*Pi*Daxial**2
         Uaxial=0.25_WP*Qaxial/Aaxial ! Divided by 4 because we have 4 ports
         ! Read in swirl gas flow rate, convert to SI, get velocity
         call this%input%read('Axial flow rate (SLPM)',Qswirl)
         Qswirl=Qswirl*SLPM2SI
         Aswirl=0.25_WP*Pi*Dswirl**2
         Uswirl=0.25_WP*Qswirl/Aswirl ! Divided by 4 because we have 4 ports
         ! Set IB velocity in vicinity of injector ports
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Axial injector at ym
                  if (sqrt((this%cfg%xm(i)-Xaxial)**2+this%cfg%zm(k)**2).le.0.5_WP*Daxial+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%y(j)+sidewall).lt.nspread*this%cfg%min_meshsize) this%Vib(i,j,k)=+Uaxial
                  ! Axial injector at yp
                  if (sqrt((this%cfg%xm(i)-Xaxial)**2+this%cfg%zm(k)**2).le.0.5_WP*Daxial+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%y(j)-sidewall).lt.nspread*this%cfg%min_meshsize) this%Vib(i,j,k)=-Uaxial
                  ! Axial injector at zm
                  if (sqrt((this%cfg%xm(i)-Xaxial)**2+this%cfg%ym(j)**2).le.0.5_WP*Daxial+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%z(k)+sidewall).lt.nspread*this%cfg%min_meshsize) this%Wib(i,j,k)=+Uaxial
                  ! Axial injector at zp
                  if (sqrt((this%cfg%xm(i)-Xaxial)**2+this%cfg%ym(j)**2).le.0.5_WP*Daxial+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%z(k)-sidewall).lt.nspread*this%cfg%min_meshsize) this%Wib(i,j,k)=-Uaxial
                  ! Swirl injector at ym
                  if (sqrt((this%cfg%xm(i)-Xswirl)**2+(this%cfg%zm(k)-Lswirl)**2).le.0.5_WP*Dswirl+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%y(j)+sidewall).lt.nspread*this%cfg%min_meshsize) this%Vib(i,j,k)=+Uswirl
                  ! Swirl injector at yp
                  if (sqrt((this%cfg%xm(i)-Xswirl)**2+(this%cfg%zm(k)+Lswirl)**2).le.0.5_WP*Dswirl+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%y(j)-sidewall).lt.nspread*this%cfg%min_meshsize) this%Vib(i,j,k)=-Uswirl
                  ! Swirl injector at zm
                  if (sqrt((this%cfg%xm(i)-Xswirl)**2+(this%cfg%ym(j)+Lswirl)**2).le.0.5_WP*Dswirl+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%z(k)+sidewall).lt.nspread*this%cfg%min_meshsize) this%Wib(i,j,k)=+Uswirl
                  ! Swirl injector at zp
                  if (sqrt((this%cfg%xm(i)-Xswirl)**2+(this%cfg%ym(j)-Lswirl)**2).le.0.5_WP*Dswirl+nspread*this%cfg%min_meshsize.and.&
                  &    abs(this%cfg%z(k)-sidewall).lt.nspread*this%cfg%min_meshsize) this%Wib(i,j,k)=-Uswirl
               end do
            end do
         end do
         ! Compute IB mass source
		   do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
				do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
					do i=this%fs%cfg%imin_,this%fs%cfg%imax_
						this%srcM(i,j,k)=this%fs%rho*(1.0_WP-this%cfg%VF(i,j,k))*(sum(this%fs%divp_x(:,i,j,k)*this%Uib(i:i+1,j,k))+&
						&                                                         sum(this%fs%divp_y(:,i,j,k)*this%Vib(i,j:j+1,k))+&
						&                                                         sum(this%fs%divp_z(:,i,j,k)*this%Wib(i,j,k:k+1)))
					end do
				end do
			end do
			call this%cfg%sync(this%srcM)
      end block set_ib_velocity
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         ! Zero velocity except if restarting
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         if (this%restarted) then
            ! Read data
            call this%df%pullvar(name='U',var=this%fs%U)
            call this%df%pullvar(name='V',var=this%fs%V)
            call this%df%pullvar(name='W',var=this%fs%W)
            !call this%df%pullvar(name='P',var=this%fs%P)  !< Reset pressure upon restart because I've noticed IB is causing drift...
            ! Apply Neumann outflow condition
            call this%fs%apply_bcond(this%time%t,this%time%dt)
         end if
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr(src=this%srcM)
         ! Compute divergence
         this%resU=this%srcM/this%fs%rho      !< Careful, we need to provide
         call this%fs%get_div(src=this%resU)  !< a volume source term to div
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='nozzle')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_nozzle')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_nozzle')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   

   !> Take one time step
   subroutine step(this)
      implicit none
      class(nozzle), intent(inout) :: this
      
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         this%resU=this%fs%rho
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         this%fs%visc=this%visc+this%sgs%visc
      end block sgs_modeling
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%rho*this%fs%U-this%fs%rho*this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%rho*this%fs%V-this%fs%rho*this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%rho*this%fs%W-this%fs%rho*this%fs%Wold)+this%time%dt*this%resW
         
         ! Form implicit residuals
         !call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho
         
         ! Apply direct IB forcing
         ibforcing: block
            integer :: i,j,k
            real(WP) :: VFx,VFy,VFz
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! Compute staggered VF
                     VFx=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                     VFy=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
                     VFz=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
                     ! Enforce IB velocity
                     if (this%fs%umask(i,j,k).eq.0) this%fs%U(i,j,k)=VFx*this%fs%U(i,j,k)+(1.0_WP-VFx)*this%Uib(i,j,k)
                     if (this%fs%vmask(i,j,k).eq.0) this%fs%V(i,j,k)=VFy*this%fs%V(i,j,k)+(1.0_WP-VFy)*this%Vib(i,j,k)
                     if (this%fs%wmask(i,j,k).eq.0) this%fs%W(i,j,k)=VFz*this%fs%W(i,j,k)+(1.0_WP-VFz)*this%Wib(i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibforcing
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%correct_mfr(src=this%srcM)
         this%resU=this%srcM/this%fs%rho      !< Careful, we need to provide
         call this%fs%get_div(src=this%resU)  !< a volume source term to div
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      this%resU=this%srcM/this%fs%rho      !< Careful, we need to provide
      call this%fs%get_div(src=this%resU)  !< a volume source term to div
      
      ! Output to ensight
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%pushval(name='t' ,val=this%time%t )
            call this%df%pushval(name='dt',val=this%time%dt)
            call this%df%pushvar(name='U' ,var=this%fs%U   )
            call this%df%pushvar(name='V' ,var=this%fs%V   )
            call this%df%pushvar(name='W' ,var=this%fs%W   )
            call this%df%pushvar(name='P' ,var=this%fs%P   )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(nozzle), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      deallocate(this%Uib,this%Vib,this%Wib,this%srcM)
      
   end subroutine final
   
   
   !> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      radius=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      if (i.eq.pg%imax+1.and.radius.ge.Rgas1.and.radius.le.Rgas2) isIn=.true.
   end function right_boundary
   
   
end module nozzle_class