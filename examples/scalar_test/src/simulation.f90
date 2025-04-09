!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use config_class,      only: config
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   public :: simulation_init,simulation_run,simulation_final
   
   ! Objects needed
   type(config)      :: cfg       !< Config
   type(timetracker) :: time      !< Time info
   type(ensight)     :: ens_out   !< Ensight output for flow variables
   type(event)       :: ens_evt   !< Event trigger for Ensight output
   type(monitor)     :: mfile     !< General simulation monitoring
   
   ! Mesh size
   real(WP) :: dx,dy,dz
   
   ! Work arrays
   real(WP), dimension(:,:,:), allocatable :: phi,psi
   real(WP), dimension(:,:,:), allocatable :: U,V,W
   real(WP), dimension(:,:,:), allocatable :: RHO,RHOold,RHOmid
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
   real(WP), dimension(:,:,:), allocatable :: Z,Zold,Zmid
   
contains
   
   
   !> Set streamfunction and velocity potential at the current time
   subroutine set_phipsi(t)
      use mathtools, only: Pi,twoPi
      implicit none
      real(WP), intent(in) :: t
      real(WP), parameter :: Tmax=5.0_WP
      real(WP), parameter :: width=0.05_WP
      real(WP), dimension(3), parameter :: center=[0.0_WP,0.1_WP,0.0_WP]
      integer :: i,j,k
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               ! Phi is at cell center
               phi(i,j,k)=0.1_WP*width**2*exp(-((cfg%xm(i)-center(1))**2+(cfg%ym(j)-center(2))**2+(cfg%zm(k)-center(3))**2)/(2.0_WP*width**2))*cos(2.0_WP*Pi*t/Tmax)
               ! Psi is at xy edge
               psi(i,j,k)=-cos(Pi*cfg%x(i))**2*cos(Pi*cfg%y(j))**2*cos(Pi*t/Tmax)*(1.0_WP/Pi)
            end do
         end do
      end do
   end subroutine set_phipsi
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Create periodic domain with uniform mesh
      create_grid: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         ! Define grid
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1)); dx=Lx/real(nx,WP)
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1)); dy=Ly/real(ny,WP)
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1)); dz=Lz/real(nz,WP)
         ! Handle 2D case
         if (nz.eq.1) then
            dz=dx; Lz=dz
         end if
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='ScalarTest')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_grid
      
      ! Allocate work arrays
      create_work_arrays: block
         allocate(phi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(psi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(U     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(V     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(W     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(RHO   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(RHOold(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(RHOmid(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(FX    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(FY    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(FZ    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Z     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Zold  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Zmid  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block create_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         call param_read('Subiterations',time%itmax,default=2)
      end block initialize_timetracker
      
      ! Initialize the problem
      initialize_problem: block
         integer :: i,j,k
         real(WP), parameter :: width=0.1_WP
         real(WP), dimension(3), parameter :: center=[-0.1_WP,0.0_WP,0.0_WP]
         ! Initialize RHO
         RHO=1.0_WP
         ! Initialize velocity
         call set_phipsi(time%tmid)
         call calc_vel()
         call interp_vel()
         ! Initialize Z
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  Z(i,j,k)=exp(-((cfg%xm(i)-center(1))**2+(cfg%ym(j)-center(2))**2+(cfg%zm(k)-center(3))**2)/(2.0_WP*width**2))
               end do
            end do
         end do 
      end block initialize_problem
      
      ! Setup ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='ScalarTest')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('RHO',RHO)
         call ens_out%add_scalar('Z',Z)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,n,ierr
      real(WP) :: RHOerr,Zerr
      real(WP), parameter :: RHOcvg=1.0e-9_WP
      real(WP), parameter ::   Zcvg=1.0e-9_WP
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time ==================
         call time%increment()
         
         ! Recalculate velocity ============
         call set_phipsi(time%tmid)
         call calc_vel()
         
         ! Advance RHO =====================
         ! Remember RHOold
         RHOold=RHO
         ! Advance RHO
         n=0
         RHOerr=huge(1.0_WP)
         do while (RHOerr.gt.RHOcvg)
            n=n+1
            RHOmid=0.5_WP*(RHO+RHOold)
            do k=cfg%kmin_,cfg%kmax_+1
               do j=cfg%jmin_,cfg%jmax_+1
                  do i=cfg%imin_,cfg%imax_+1
                     FX(i,j,k)=0.5_WP*(RHOmid(i,j,k)+RHOmid(i-1,j,k))*U(i,j,k)
                     FY(i,j,k)=0.5_WP*(RHOmid(i,j,k)+RHOmid(i,j-1,k))*V(i,j,k)
                     FZ(i,j,k)=0.5_WP*(RHOmid(i,j,k)+RHOmid(i,j,k-1))*W(i,j,k)
                  end do
               end do
            end do
            RHOmid=RHO
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     RHO(i,j,k)=RHOold(i,j,k)-time%dt*((FX(i+1,j,k)-FX(i,j,k))/dx+(FY(i,j+1,k)-FY(i,j,k))/dy+(FZ(i,j,k+1)-FZ(i,j,k))/dz)
                  end do
               end do
            end do
            call cfg%sync(RHO)
            RHOerr=maxval(abs(RHO-RHOmid))
            call MPI_ALLREDUCE(MPI_IN_PLACE,RHOerr,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         end do
         if (cfg%amroot) print*,'RHO converged!',n,RHOerr
         
         ! Advance Z =======================
         ! Remember Zold
         Zold=Z
         ! Set RHOmid
         RHOmid=0.5_WP*(RHO+RHOold)
         ! Advance Z
         n=0
         Zerr=huge(1.0_WP)
         do while (Zerr.gt.Zcvg)
            n=n+1
            Zmid=(sqrt(RHO)*Z+sqrt(RHOold)*Zold)/(sqrt(RHO)+sqrt(RHOold))
            do k=cfg%kmin_,cfg%kmax_+1
               do j=cfg%jmin_,cfg%jmax_+1
                  do i=cfg%imin_,cfg%imax_+1
                     FX(i,j,k)=0.5_WP*(RHOmid(i,j,k)+RHOmid(i-1,j,k))*U(i,j,k)*0.5_WP*(Zmid(i,j,k)+Zmid(i-1,j,k))
                     FY(i,j,k)=0.5_WP*(RHOmid(i,j,k)+RHOmid(i,j-1,k))*V(i,j,k)*0.5_WP*(Zmid(i,j,k)+Zmid(i,j-1,k))
                     FZ(i,j,k)=0.5_WP*(RHOmid(i,j,k)+RHOmid(i,j,k-1))*W(i,j,k)*0.5_WP*(Zmid(i,j,k)+Zmid(i,j,k-1))
                  end do
               end do
            end do
            Zmid=Z
            do k=cfg%kmin_,cfg%kmax_
               do j=cfg%jmin_,cfg%jmax_
                  do i=cfg%imin_,cfg%imax_
                     Z(i,j,k)=(RHOold(i,j,k)*Zold(i,j,k)-time%dt*((FX(i+1,j,k)-FX(i,j,k))/dx+(FY(i,j+1,k)-FY(i,j,k))/dy+(FZ(i,j,k+1)-FZ(i,j,k))/dz))/RHO(i,j,k)
                  end do
               end do
            end do
            call cfg%sync(Z)
            Zerr=maxval(abs(Z-Zmid))
            call MPI_ALLREDUCE(MPI_IN_PLACE,Zerr,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         end do
         if (cfg%amroot) print*,'Z converged!',n,Zerr
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call interp_vel()
            call ens_out%write_data(time%t)
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      deallocate(U,V,W,Ui,Vi,Wi,RHO,RHOmid,RHOold,phi,psi,FX,FY,FZ,Z,Zold,Zmid)
   end subroutine simulation_final
   
   
   !> Calculate flow field from streamfunction and velocity potential
   subroutine calc_vel()
      implicit none
      integer :: i,j,k
      ! Calculate flow from velocity potential
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               ! U=+d(phi)/dx
               U(i,j,k)=(phi(i,j,k)-phi(i-1,j,k))/dx
               ! V=+d(phi)/dy
               V(i,j,k)=(phi(i,j,k)-phi(i,j-1,k))/dy
               ! W=+d(phi)/dz
               W(i,j,k)=(phi(i,j,k)-phi(i,j,k-1))/dz
            end do
         end do
      end do
      ! Calculate flow from (xy) streamfunction
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               ! U=+d(psi)/dy
               U(i,j,k)=U(i,j,k)+(psi(i,j+1,k)-psi(i,j,k))/dy
               ! V=-d(psi)/dx
               V(i,j,k)=V(i,j,k)-(psi(i+1,j,k)-psi(i,j,k))/dx
            end do
         end do
      end do
      ! Communicate
      call cfg%sync(U)
      call cfg%sync(V)
      call cfg%sync(W)
   end subroutine calc_vel
   

   !> Calculate interpolated velocity
   subroutine interp_vel()
      implicit none
      integer :: i,j,k
      do k=cfg%kmino_,cfg%kmaxo_-1
         do j=cfg%jmino_,cfg%jmaxo_-1
            do i=cfg%imino_,cfg%imaxo_-1
               Ui(i,j,k)=0.5_WP*(U(i,j,k)+U(i+1,j,k))
               Vi(i,j,k)=0.5_WP*(V(i,j,k)+V(i,j+1,k))
               Wi(i,j,k)=0.5_WP*(W(i,j,k)+W(i,j,k+1))
            end do
         end do
      end do
      ! Sync it
      call cfg%sync(Ui)
      call cfg%sync(Vi)
      call cfg%sync(Wi)
   end subroutine interp_vel


end module simulation
