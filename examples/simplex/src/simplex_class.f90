!> Definition for a simplex class
module simplex_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use polygon_class,     only: polygon
   use ensight_class,     only: ensight
   use fft2d_class,       only: fft2d
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
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
      procedure :: init                            !< Initialize simplex simulation
      procedure :: step                            !< Advance simplex simulation by one time step
      procedure :: final                           !< Finalize simplex simulation
   end type simplex
   
   
contains
   
   
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
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
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
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='simplex')
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
      end block create_simplex
      
      
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
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! We are restarting, read the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(adjustl(timestamp)))
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=4)
            this%df%valname=['t ','dt']
            this%df%varname=['U','V','W','P']
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
      
      
      ! Initialize our IB velocity field - our Dirichlet boundary conditions
      set_ib_velocity: block
         use mathtools, only: Pi
         integer  :: i,j,k
         real(WP) :: R,mfr,theta,rx,ry,rz,Uin
         real(WP), dimension(3) :: v1,v2,n1,n2
         real(WP), dimension(3) :: xf,yf,zf
         real(WP), dimension(3) :: xp,yp,zp
         ! Zero initial field
         this%Uib=0.0_WP; this%Vib=0.0_WP; this%Wib=0.0_WP
         ! Define geometry of inlet pipes
         R=0.000185_WP
         theta=Pi/180.0_WP*53.0_WP
         v1=[-0.00442_WP,0.0_WP,+0.001245_WP]; n1=[sin(theta),+cos(theta),0.0_WP]
         v2=[-0.00442_WP,0.0_WP,-0.001245_WP]; n2=[sin(theta),-cos(theta),0.0_WP]
         ! Read in mass flow rate
         call this%input%read('Mass flow rate',mfr)
         Uin=0.5_WP*mfr/(this%fs%rho*Pi*R**2)
         ! Apply boundary conditions within the IBs
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Inlet pipe 1
                  ! Staggered positions
                  xf=[this%cfg%x(i) ,this%cfg%ym(j),this%cfg%zm(k)]-v1
                  yf=[this%cfg%xm(i),this%cfg%y(j) ,this%cfg%zm(k)]-v1
                  zf=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%z(k) ]-v1
                  ! Project onto pipe axis
                  xp=xf-n1*dot_product(xf,n1)
                  yp=yf-n1*dot_product(yf,n1)
                  zp=zf-n1*dot_product(zf,n1)
                  ! Measure radial position in pipe 1
                  rx=sqrt(dot_product(xp,xp))/R
                  ry=sqrt(dot_product(yp,yp))/R
                  rz=sqrt(dot_product(zp,zp))/R
                  ! Apply Poiseuille profile
                  if (xf(1).le.0.0_WP) then
                     if (rx.le.1.0_WP) this%Uib(i,j,k)=Uin*n1(1)*2.0_WP*(1.0_WP-rx**2)
                     if (ry.le.1.0_WP) this%Vib(i,j,k)=Uin*n1(2)*2.0_WP*(1.0_WP-ry**2)
                     if (rz.le.1.0_WP) this%Wib(i,j,k)=Uin*n1(3)*2.0_WP*(1.0_WP-rz**2)
                  end if
                  ! Inlet pipe 2
                  ! Staggered positions
                  xf=[this%cfg%x(i) ,this%cfg%ym(j),this%cfg%zm(k)]-v2
                  yf=[this%cfg%xm(i),this%cfg%y(j) ,this%cfg%zm(k)]-v2
                  zf=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%z(k) ]-v2
                  ! Project onto pipe axis
                  xp=xf-n2*dot_product(xf,n2)
                  yp=yf-n2*dot_product(yf,n2)
                  zp=zf-n2*dot_product(zf,n2)
                  ! Measure radial position in pipe 2
                  rx=sqrt(dot_product(xp,xp))/R
                  ry=sqrt(dot_product(yp,yp))/R
                  rz=sqrt(dot_product(zp,zp))/R
                  ! Apply Poiseuille profile
                  if (xf(1).le.0.0_WP) then
                     if (rx.le.1.0_WP) this%Uib(i,j,k)=Uin*n2(1)*2.0_WP*(1.0_WP-rx**2)
                     if (ry.le.1.0_WP) this%Vib(i,j,k)=Uin*n2(2)*2.0_WP*(1.0_WP-ry**2)
                     if (rz.le.1.0_WP) this%Wib(i,j,k)=Uin*n2(3)*2.0_WP*(1.0_WP-rz**2)
                  end if
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
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            !call this%df%pullvar(name='P',var=this%fs%P)  !< Reset pressure upon restart because I've noticed IB is causing drift...
            ! Apply boundary conditions
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
         this%ens_out=ensight(cfg=this%cfg,name='simplex')
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
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_simplex')
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
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_simplex')
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
      
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(simplex), intent(inout) :: this
      
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
            call this%df%push(name='t' ,val=this%time%t )
            call this%df%push(name='dt',val=this%time%dt)
            call this%df%push(name='U' ,var=this%fs%U   )
            call this%df%push(name='V' ,var=this%fs%V   )
            call this%df%push(name='W' ,var=this%fs%W   )
            call this%df%push(name='P' ,var=this%fs%P   )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   

   !> Finalize simplex simulation
   subroutine final(this)
      implicit none
      class(simplex), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      deallocate(this%Uib,this%Vib,this%Wib,this%srcM)
      
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
   
   
end module simplex_class