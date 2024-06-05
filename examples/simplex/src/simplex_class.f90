!> Definition for a simplex class
module simplex_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use polygon_class,     only: polygon
   use ensight_class,     only: ensight
   use fft2d_class,       only: fft2d
   !use hypre_str_class,   only: hypre_str
   !use hypre_uns_class,   only: hypre_uns
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
      !type(hypre_str)   :: ps    !< HYPRE linear solver for pressure
      !type(hypre_uns)   :: ps    !< HYPRE linear solver for pressure
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
      procedure :: init                            !< Initialize simplex simulation
      procedure :: step                            !< Advance simplex simulation by one time step
      procedure :: final                           !< Finalize simplex simulation
   end type simplex
   
   !> Inlet pipes geometry
   real(WP), parameter :: Rpipe=0.000185_WP
   real(WP), dimension(3), parameter :: p1=[-0.00442_WP,0.0_WP,+0.001245_WP]
   real(WP), dimension(3), parameter :: p2=[-0.00442_WP,0.0_WP,-0.001245_WP]
   real(WP), dimension(3), parameter :: n1=[+0.6_WP,-0.8_WP,0.0_WP]
   real(WP), dimension(3), parameter :: n2=[+0.6_WP,+0.8_WP,0.0_WP]
   real(WP) :: mfr,Apipe
   
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
         real(WP), dimension(3) :: v,p
         real(WP) :: r
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
         ! Carve out stair-stepped inlet pipes
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Inlet pipe 1
                  v=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]-p1
                  p=v-n1*dot_product(v,n1)
                  r=sqrt(dot_product(p,p))/Rpipe
                  if (v(1).le.this%cfg%min_meshsize.and.r.le.1.0_WP) this%cfg%VF(i,j,k)=1.0_WP
                  ! Inlet pipe 2
                  v=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]-p2
                  p=v-n2*dot_product(v,n2)
                  r=sqrt(dot_product(p,p))/Rpipe
                  if (v(1).le.this%cfg%min_meshsize.and.r.le.1.0_WP) this%cfg%VF(i,j,k)=1.0_WP
               end do
            end do
         end do
         ! Apply Neumann on VF at entrance
         if (this%cfg%iproc.eq.1) this%cfg%VF(this%cfg%imino:this%cfg%imin-1,:,:)=this%cfg%VF(this%cfg%imino:this%cfg%imin,:,:)
         ! Recompute domain volume
         call this%cfg%calc_fluid_vol()
      end block create_simplex
      
      
      ! Initialize flow rate
      set_flow_rate: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: j,k,ierr
         ! Read mass flow rate
         call this%input%read('Mass flow rate',mfr)
         ! Integrate inlet pipe surface area
         Apipe=0.0_WP
         if (this%cfg%iproc.eq.1) then
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).lt.0.002_WP) then
                     Apipe=Apipe+this%cfg%VF(this%cfg%imin-1,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                  end if
               end do
            end do
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,Apipe,this%cfg%nproc,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      end block set_flow_rate
      
      
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
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use incomp_class,    only: clipped_neumann,dirichlet
         !use hypre_str_class, only: pcg_pfmg2
         !use hypre_uns_class, only: pcg_amg
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call this%input%read('Density',this%fs%rho)
         call this%input%read('Dynamic viscosity',this%visc); this%fs%visc=this%visc
         ! Inlets on the left
         call this%fs%add_bcond(name='inlets',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=pipe_inlets)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary)
         ! Configure pressure solver
         this%ps=fft2d(cfg=this%cfg,name='Pressure',nst=7)
         !this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         !this%ps=hypre_uns(cfg=this%cfg,name='Pressure',method=pcg_amg,nst=7)
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
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
            ! Apply boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
         end if
         ! Apply Dirichlet condition at pipe inlets
         call this%fs%get_bcond('inlets',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*mfr/(this%fs%rho*Apipe)
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
                     if (this%fs%umask(i,j,k).eq.0) this%fs%U(i,j,k)=VFx*this%fs%U(i,j,k)
                     if (this%fs%vmask(i,j,k).eq.0) this%fs%V(i,j,k)=VFy*this%fs%V(i,j,k)
                     if (this%fs%wmask(i,j,k).eq.0) this%fs%W(i,j,k)=VFz*this%fs%W(i,j,k)
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
   
   
   !> Function that localizes the pipe inlets
   function pipe_inlets(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.sqrt(pg%ym(j)**2+pg%zm(k)**2).lt.0.002_WP) isIn=.true.
   end function pipe_inlets
   
   
end module simplex_class