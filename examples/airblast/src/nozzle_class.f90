!> Definition for a nozzle class
module nozzle_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: nozzle
   
   !> Nozzle object
   type :: nozzle
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB
      type(ibconfig) :: cfg
      
		!> Surface mesh for IB
      type(surfmesh) :: plymesh
      
		!> Flow solver
		type(incomp)      :: fs    !< Incompressible flow solver
		type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
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
		
		!> Fluid definition
      real(WP) :: visc
		
   contains
	   procedure, private :: geometry_init          !< Initialize geometry for nozzle
		procedure, private :: simulation_init        !< Initialize simulation for nozzle
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type nozzle

	
	!> Hardcode nozzle inlet positions used in locator functions
	real(WP), parameter, public :: axial_xdist =-0.06394704_WP
	real(WP), parameter, public :: axial_diam  =+0.01600000_WP*1.2_WP  ! 20% larger to avoid stairstepping
	real(WP), parameter, public :: swirl_xdist =-0.06394704_WP
	real(WP), parameter, public :: swirl_diam  =+0.01000000_WP*1.2_WP  ! 20% larger to avoid stairstepping
	real(WP), parameter, public :: swirl_offset=+0.02945130_WP
	

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
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='nozzle')

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
         use string, only: str_medium
         use param,  only: param_read
         character(len=str_medium) :: plyfile
         
         ! Read in ply filename
         call this%input%read('PLY filename',plyfile)
         
         ! Create surface mesh from ply
         this%plymesh=surfmesh(comm=this%cfg%comm,plyfile=plyfile,nvar=0,name='ply')
         
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
         
      end block create_walls


      this%ens_out=ensight(cfg=this%cfg,name='nozzle')
      call this%ens_out%add_scalar('Gib',this%cfg%Gib)
      call this%ens_out%add_surface('ply',this%plymesh)
      call this%ens_out%write_data(0.0_WP)
      stop
      
      
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
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use incomp_class,    only: dirichlet,clipped_neumann
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call this%input%read('Density',this%fs%rho)
         call this%input%read('Dynamic viscosity',this%visc); this%fs%visc=this%visc
         ! Define gas port boundary conditions
         call this%fs%add_bcond(name='axial_ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=axial_ym)
         call this%fs%add_bcond(name='axial_yp',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=axial_yp)
         call this%fs%add_bcond(name='axial_zm',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=axial_zm)
         call this%fs%add_bcond(name='axial_zp',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=axial_zp)
         call this%fs%add_bcond(name='swirl_ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=swirl_ym)
         call this%fs%add_bcond(name='swirl_yp',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=swirl_yp)
         call this%fs%add_bcond(name='swirl_zm',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=swirl_zm)
         call this%fs%add_bcond(name='swirl_zp',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=swirl_zp)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         this%ps%maxlevel=20
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)
      end block create_flow_solver
      

      ! Initialize our velocity field
      initialize_velocity: block
         use mpi_f08,      only: MPI_ALLREDUCE,MPI_SUM
         use parallel,     only: MPI_REAL_WP
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k,ierr
         real(WP) :: Uaxial,myAaxial,Aaxial,Uswirl,myAswirl,Aswirl
         real(WP) :: Qaxial,Qswirl
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Read in axial gas flow rate and convert to SI
         call this%input%read('Axial flow rate (SLPM)',Qaxial)
         Qaxial=Qaxial*SLPM2SI
         ! Calculate axial flow area
         myAaxial=0.0_WP
         call this%fs%get_bcond('axial_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+this%cfg%dz(k)*this%cfg%dx(i)*sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
         end do
         call this%fs%get_bcond('axial_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+this%cfg%dz(k)*this%cfg%dx(i)*sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
         end do
         call this%fs%get_bcond('axial_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+this%cfg%dx(i)*this%cfg%dy(j)*sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
         end do
         call this%fs%get_bcond('axial_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+this%cfg%dx(i)*this%cfg%dy(j)*sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
         end do
         call MPI_ALLREDUCE(myAaxial,Aaxial,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! Calculate bulk axial velocity
         Uaxial=Qaxial/Aaxial
         ! Apply Dirichlet at 4 axial injector ports
         call this%fs%get_bcond('axial_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%V(i,j,k)=+sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*Uaxial
         end do
         call this%fs%get_bcond('axial_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%V(i,j,k)=-sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*Uaxial
         end do
         call this%fs%get_bcond('axial_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%W(i,j,k)=+sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*Uaxial
         end do
         call this%fs%get_bcond('axial_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%W(i,j,k)=-sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*Uaxial
         end do
         ! Read in swirl gas flow rate and convert to SI
         call this%input%read('Axial flow rate (SLPM)',Qswirl)
         Qswirl=Qswirl*SLPM2SI
         ! Calculate swirl flow area
         myAswirl=0.0_WP
         call this%fs%get_bcond('swirl_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+this%cfg%dz(k)*this%cfg%dx(i)*sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
         end do
         call this%fs%get_bcond('swirl_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+this%cfg%dz(k)*this%cfg%dx(i)*sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
         end do
         call this%fs%get_bcond('swirl_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+this%cfg%dx(i)*this%cfg%dy(j)*sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
         end do
         call this%fs%get_bcond('swirl_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+this%cfg%dx(i)*this%cfg%dy(j)*sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
         end do
         call MPI_ALLREDUCE(myAswirl,Aswirl,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! Calculate bulk axial velocity
         Uswirl=Qswirl/Aswirl
         ! Apply Dirichlet at 4 axial injector ports
         call this%fs%get_bcond('swirl_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%V(i,j,k)=+sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*Uswirl
         end do
         call this%fs%get_bcond('swirl_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%V(i,j,k)=-sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*Uswirl
         end do
         call this%fs%get_bcond('swirl_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%W(i,j,k)=+sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*Uswirl
         end do
         call this%fs%get_bcond('swirl_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%W(i,j,k)=-sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*Uswirl
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Compute MFR through all boundary conditions
         call this%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
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
			
			! Solve Poisson equation
			call this%fs%correct_mfr()
			call this%fs%get_div()
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
		call this%fs%get_div()
		
		! Output to ensight
		if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
		
		! Perform and output monitoring
		call this%fs%get_max()
		call this%mfile%write()
		call this%cflfile%write()
		
		
	end subroutine step
	

   !> Finalize nozzle simulation
   subroutine final(this)
		implicit none
		class(nozzle), intent(inout) :: this
		
		! Deallocate work arrays
		deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
		
	end subroutine final
   
   
	!> Function that localizes the right domain boundary
   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary
   

   !> Function that localizes axial injector at -y
   function axial_ym(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-z plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[axial_xdist,0.0_WP,0.0_WP])
      rise=abs(pg%ym(j))
      run =sqrt(hyp**2-rise**2)
      if (run.le.axial_diam/2.0_WP.and.j.eq.pg%jmin) isIn=.true.
   end function axial_ym
   
   
   !> Function that localizes axial injector at +y
   function axial_yp(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-z plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[axial_xdist,0.0_WP,0.0_WP])
      rise=abs(pg%ym(j))
      run =sqrt(hyp**2-rise**2)
      if (run.le.axial_diam/2.0_WP.and.j.eq.pg%jmax+1) isIn=.true.
   end function axial_yp
   

   !> Function that localizes axial injector at -z
   function axial_zm(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-y plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[axial_xdist,0.0_WP,0.0_WP])
      rise=abs(pg%zm(k))
      run =sqrt(hyp**2-rise**2)
      if (run.le.axial_diam/2.0_WP.and.k.eq.pg%kmin) isIn=.true.
   end function axial_zm
   

   !> Function that localizes axial injector at +z
   function axial_zp(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-y plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[axial_xdist,0.0_WP,0.0_WP])
      rise=abs(pg%zm(k))
      run =sqrt(hyp**2-rise**2)
      if (run.le.axial_diam/2.0_WP.and.k.eq.pg%kmax+1) isIn=.true.
   end function axial_zp


   !> Function that localizes swirl injector at -y
   function swirl_ym(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-z plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[swirl_xdist,0.0_WP,+swirl_offset])
      rise=abs(pg%ym(j))
      run =sqrt(hyp**2-rise**2)
      if (run.le.swirl_diam/2.0_WP.and.j.eq.pg%jmin) isIn=.true.
   end function swirl_ym
   
   
   !> Function that localizes swirl injector at +y
   function swirl_yp(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-z plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[swirl_xdist,0.0_WP,-swirl_offset])
      rise=abs(pg%ym(j))
      run =sqrt(hyp**2-rise**2)
      if (run.le.swirl_diam/2.0_WP.and.j.eq.pg%jmax+1) isIn=.true.
   end function swirl_yp
   

   !> Function that localizes swirl injector at -z
   function swirl_zm(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-y plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[swirl_xdist,-swirl_offset,0.0_WP])
      rise=abs(pg%zm(k))
      run =sqrt(hyp**2-rise**2)
      if (run.le.swirl_diam/2.0_WP.and.k.eq.pg%kmin) isIn=.true.
   end function swirl_zm
   

   !> Function that localizes swirl injector at +z
   function swirl_zp(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: hyp,rise,run
      logical :: isIn
      isIn=.false.
      ! Check injector x-y plane
      hyp =norm2([pg%xm(i),pg%ym(j),pg%zm(k)]-[swirl_xdist,+swirl_offset,0.0_WP])
      rise=abs(pg%zm(k))
      run =sqrt(hyp**2-rise**2)
      if (run.le.swirl_diam/2.0_WP.and.k.eq.pg%kmax+1) isIn=.true.
   end function swirl_zp
	

end module nozzle_class