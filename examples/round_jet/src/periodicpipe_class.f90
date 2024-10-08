module periodicpipe_class
   use precision,         only: WP
   use ibconfig_class,    only: ibconfig
   use fft3d_class,       only: fft3d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   public :: periodicpipe
   
   !> periodicpipe object
   type :: periodicpipe
      
      !> Config
      type(ibconfig) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs     !< Incompressible flow solver
      type(fft3d)       :: ps     !< Fourier-accelerated pressure solver
      type(timetracker) :: time   !< Time info

      !> Implicit solver
      logical     :: use_implicit !< Is an implicit solver used?
      type(ddadi) :: vs           !< DDADI solver for velocity
      
      !> SGS model
      logical        :: use_sgs   !< Is an LES model used?
      type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
      
      !> Ensight postprocessing
      !type(ensight) :: ens_out    !< Event trigger for Ensight output
      !type(event)   :: ens_evt    !< Ensight output for flow variables
      
      !> Simulation monitoring files
      type(monitor) :: mfile      !< General simulation monitoring
      !type(monitor) :: cflfile    !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Flow conditions
      real(WP) :: visc,Ubulk,meanU,bforce
      
      !> Provide a pardata object for restarts
      logical       :: restarted
      type(pardata) :: df
      type(event)   :: save_evt
      
   contains
      procedure :: init                            !< Initialize pipe simulation
      procedure :: step                            !< Advance pipe simulation by one time step
      procedure :: final                           !< Finalize pipe simulation
   end type periodicpipe

contains
   
   
   !> Initialization of periodicpipe simulation
   subroutine init(this,group)
      use mpi_f08, only: MPI_Group
      use param,   only: param_read
      implicit none
      class(periodicpipe), intent(inout) :: this
      type(MPI_Group), intent(in) :: group
      

      ! Initialize the config
      initialize_config: block
         use sgrid_class,    only: sgrid,cartesian
         use ibconfig_class, only: sharp
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx,D
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         
         ! Nominal parameters
         D=1.0_WP !< Unity diameter
         no=2     !< Allow for two dead cells for IB to take effect
         
         ! Read in grid definition
         call param_read('[Pipe] Length',Lx)
         call param_read('[Pipe] nx',nx); allocate(x(nx+1))
         call param_read('[Pipe] ny',ny); allocate(y(ny+1))
         call param_read('[Pipe] nz',nz); allocate(z(nz+1))
         dx=Lx/real(nx,WP)
         
         ! Adjust domain size to account for extra cells
         if (ny.gt.1) then
            Ly=D+real(2*no,WP)*D/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=D+real(2*no,WP)*D/real(ny-2*no,WP)
         else
            Lz=dx
         end if
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='pipe')
         
         ! Read in partition
         call param_read('[Pipe] Partition',partition,short='p')
         
         ! Create partitioned grid
         this%cfg=ibconfig(grp=group,decomp=partition,grid=grid)
         
         ! Create masks for this config
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%cfg%Gib(i,j,k)=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)-0.5_WP*D
               end do
            end do
         end do
         ! Get normal vector
         call this%cfg%calculate_normal()
         ! Get VF field
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
         
      end block initialize_config
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('[Pipe] Max timestep size',this%time%dtmax)
         call param_read('[Pipe] Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         this%fs%rho=1.0_WP
         call param_read('Reynolds number',this%visc)
         this%visc=1.0_WP/this%visc
         this%fs%visc=this%visc
         ! Configure pressure solver
         this%ps=fft3d(cfg=this%cfg,name='Pressure',nst=7)
         ! Check if implicit velocity solver is used
         call param_read('[Pipe] Use implicit solver',this%use_implicit)
         if (this%use_implicit) then
            ! Configure implicit solver
            this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
            ! Finish flow solver setup
            call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         else
            ! Finish flow solver setup
            call this%fs%setup(pressure_solver=this%ps)
         end if
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: amp,VFx,VFy,VFz
         ! Initial fields
         this%Ubulk=1.0_WP
         this%fs%U=this%Ubulk; this%fs%V=0.0_WP; this%fs%W=0.0_WP; this%fs%P=0.0_WP
         this%meanU=this%Ubulk
         this%bforce=0.0_WP
         ! For faster transition
         call param_read('[Pipe] Fluctuation amp',amp,default=0.2_WP)
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  ! Add fluctuations for faster transition
                  this%fs%U(i,j,k)=this%fs%U(i,j,k)+this%Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*this%Ubulk*cos(8.0_WP*twoPi*this%fs%cfg%zm(k)/this%fs%cfg%zL)*cos(8.0_WP*twoPi*this%fs%cfg%ym(j)/this%fs%cfg%yL)
                  this%fs%V(i,j,k)=this%fs%V(i,j,k)+this%Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*this%Ubulk*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)
                  this%fs%W(i,j,k)=this%fs%W(i,j,k)+this%Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*this%Ubulk*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)
                  ! Remove values in the wall
                  VFx=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                  VFy=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
                  VFz=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
                  this%fs%U(i,j,k)=this%fs%U(i,j,k)*VFx
                  this%fs%V(i,j,k)=this%fs%V(i,j,k)*VFy
                  this%fs%W(i,j,k)=this%fs%W(i,j,k)*VFz
               end do
            end do
         end do
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%V)
         call this%fs%cfg%sync(this%fs%W)
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('[Pipe] Use SGS model',this%use_sgs)
         if (this%use_sgs) this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs


      ! Handle restart here
      perform_restart: block
         use string,  only: str_medium
         use filesys, only: makedir,isdir
         character(len=str_medium) :: filename
         integer, dimension(3) :: iopartition
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Pipe restart output')
         call param_read('[Pipe] Restart output period',this%save_evt%tper)
         ! Read in the partition for I/O
         call param_read('[Pipe] I/O partition',iopartition)
         ! Check if a restart file was provided
         call param_read('[Pipe] Restart from',filename,default='')
         this%restarted=.false.; if (len_trim(filename).gt.0) this%restarted=.true.
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata=trim(filename))
            ! Put the data at the right place
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            ! Update cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Update divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
         else
            ! Prepare a new directory for storing files for restart
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! If we are not restarting, we will still need a datafile for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=4)
            this%df%valname=['dt','t ']; this%df%varname=['U','V','W','P']
         end if
      end block perform_restart
      
      
      ! Add Ensight output
      !create_ensight: block
      !   ! Create Ensight output from cfg
      !   this%ens_out=ensight(cfg=this%cfg,name='pipe')
      !   ! Create event for Ensight output
      !   this%ens_evt=event(time=this%time,name='Ensight output')
      !   call param_read('Ensight output period',this%ens_evt%tper)
      !   ! Add variables to output
      !   call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
      !   call this%ens_out%add_scalar('Gib',this%cfg%Gib)
      !   call this%ens_out%add_scalar('pressure',this%fs%P)
      !   if (this%use_sgs) call this%ens_out%add_scalar('visc_sgs',this%sgs%visc)
      !   ! Output to ensight
      !   if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      !end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'pipe')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%meanU,'Bulk U')
         call this%mfile%add_column(this%bforce,'Body force')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         !this%cflfile=monitor(this%fs%cfg%amRoot,'pipe_cfl')
         !call this%cflfile%add_column(this%time%n,'Timestep number')
         !call this%cflfile%add_column(this%time%t,'Time')
         !all this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         !all this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         !call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         !call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         !call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         !call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         !call this%cflfile%write()
      end block create_monitor
      

   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(periodicpipe), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Turbulence modeling
      if (this%use_sgs) then
         sgs_modeling: block
            use sgsmodel_class, only: vreman
            this%resU=this%fs%rho
            call this%fs%get_gradu(this%gradU)
            call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
            where (this%cfg%Gib.gt.0.0_WP) this%sgs%visc=0.0_WP
            this%fs%visc=this%visc+this%sgs%visc
         end block sgs_modeling
      end if
      
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
         
         ! Add body forcing
         forcing: block
            use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
            use parallel, only: MPI_REAL_WP
            integer :: i,j,k,ierr
            real(WP) :: Uvol,VFx
            Uvol=0.0_WP; this%meanU=0.0_WP
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     VFx=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                     if (VFx.le.0.5_WP) cycle
                     this%meanU=this%meanU+this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)*VFx*(2.0_WP*this%fs%U(i,j,k)-this%fs%Uold(i,j,k))
                     Uvol      =Uvol      +this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)*VFx
                  end do
               end do
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,Uvol      ,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%meanU,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%meanU=this%meanU/Uvol
            this%resU=this%resU+this%fs%rho*(this%Ubulk-this%meanU)
            this%bforce=this%fs%rho*(this%Ubulk-this%meanU)/this%time%dt
         end block forcing
         
         ! Finish update
         if (this%use_implicit) then
            ! Form implicit residuals
            call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
            ! Apply these residuals
            this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
            this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
            this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         else
            ! Apply these residuals
            this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho
            this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho
            this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho
         end if
         
         ! Apply IB forcing to enforce BC at the pipe walls
         ibforcing: block
            integer :: i,j,k
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     this%fs%U(i,j,k)=this%fs%U(i,j,k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                     this%fs%V(i,j,k)=this%fs%V(i,j,k)*sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
                     this%fs%W(i,j,k)=this%fs%W(i,j,k)*sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
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
      !if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      !call this%cflfile%write()
      
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
            call this%df%write(fdata='restart/pipe_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   
   
   !> Finalize simulation
   subroutine final(this)
      implicit none
      class(periodicpipe), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      
   end subroutine final
   
   
end module periodicpipe_class