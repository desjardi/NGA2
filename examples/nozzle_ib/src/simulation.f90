!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use geometry,          only: plymesh
   use geometry,          only: axial_xdist,axial_diam,swirl_xdist,swirl_diam,swirl_offset
   use hypre_str_class,   only: hypre_str
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp),      public :: fs
   type(hypre_str),   public :: ps
   type(sgsmodel),    public :: sgs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight)  :: ens_out,ply_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Fluid definition
   real(WP) :: visc
   

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
   

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         use incomp_class,    only: dirichlet,clipped_neumann
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define gas port boundary conditions
         call fs%add_bcond(name='axial_ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=axial_ym)
         call fs%add_bcond(name='axial_yp',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=axial_yp)
         call fs%add_bcond(name='axial_zm',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=axial_zm)
         call fs%add_bcond(name='axial_zp',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=axial_zp)
         call fs%add_bcond(name='swirl_ym',type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=swirl_ym)
         call fs%add_bcond(name='swirl_yp',type=dirichlet,face='y',dir=+1,canCorrect=.false.,locator=swirl_yp)
         call fs%add_bcond(name='swirl_zm',type=dirichlet,face='z',dir=-1,canCorrect=.false.,locator=swirl_zm)
         call fs%add_bcond(name='swirl_zp',type=dirichlet,face='z',dir=+1,canCorrect=.false.,locator=swirl_zp)
         ! Outflow on the right
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         !ps%maxlevel=14
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
      end block create_flow_solver
      

      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools,    only: pi
         use mpi_f08,      only: MPI_ALLREDUCE,MPI_SUM
         use parallel,     only: MPI_REAL_WP
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k,ierr
         real(WP) :: Uaxial,myAaxial,Aaxial,Uswirl,myAswirl,Aswirl
         real(WP) :: Qaxial,Qswirl
         real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Read in axial gas flow rate and convert to SI
         call param_read('Axial flow rate (SLPM)',Qaxial)
         Qaxial=Qaxial*SLPM2SI
         ! Calculate axial flow area
         myAaxial=0.0_WP
         call fs%get_bcond('axial_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+cfg%dz(k)*cfg%dx(i)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
         end do
         call fs%get_bcond('axial_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+cfg%dz(k)*cfg%dx(i)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
         end do
         call fs%get_bcond('axial_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+cfg%dx(i)*cfg%dy(j)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
         end do
         call fs%get_bcond('axial_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAaxial=myAaxial+cfg%dx(i)*cfg%dy(j)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
         end do
         call MPI_ALLREDUCE(myAaxial,Aaxial,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
         ! Calculate bulk axial velocity
         Uaxial=Qaxial/Aaxial
         ! Apply Dirichlet at 4 axial injector ports
         call fs%get_bcond('axial_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=+sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*Uaxial
         end do
         call fs%get_bcond('axial_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=-sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*Uaxial
         end do
         call fs%get_bcond('axial_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%W(i,j,k)=+sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*Uaxial
         end do
         call fs%get_bcond('axial_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%W(i,j,k)=-sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*Uaxial
         end do
         ! Read in swirl gas flow rate and convert to SI
         call param_read('Axial flow rate (SLPM)',Qswirl)
         Qswirl=Qswirl*SLPM2SI
         ! Calculate swirl flow area
         myAswirl=0.0_WP
         call fs%get_bcond('swirl_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+cfg%dz(k)*cfg%dx(i)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
         end do
         call fs%get_bcond('swirl_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+cfg%dz(k)*cfg%dx(i)*sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
         end do
         call fs%get_bcond('swirl_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+cfg%dx(i)*cfg%dy(j)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
         end do
         call fs%get_bcond('swirl_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            myAswirl=myAswirl+cfg%dx(i)*cfg%dy(j)*sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
         end do
         call MPI_ALLREDUCE(myAswirl,Aswirl,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
         ! Calculate bulk axial velocity
         Uswirl=Qswirl/Aswirl
         ! Apply Dirichlet at 4 axial injector ports
         call fs%get_bcond('swirl_ym',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=+sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*Uswirl
         end do
         call fs%get_bcond('swirl_yp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=-sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*Uswirl
         end do
         call fs%get_bcond('swirl_zm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%W(i,j,k)=+sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*Uswirl
         end do
         call fs%get_bcond('swirl_zp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%W(i,j,k)=-sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*Uswirl
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity


      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Output ply mesh
         ply_out=ensight(cfg=cfg,name='nozzle')
         call ply_out%add_surface('ply',plymesh)
         call ply_out%write_data(time%t)
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='nozzle')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('levelset',cfg%Gib)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('visc_sgs',sgs%visc)
         call ens_out%add_scalar('div',fs%div)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Turbulence modeling
         sgs_modeling: block
            use sgsmodel_class, only: vreman
            resU=fs%rho
            call fs%get_gradu(gradU)
            call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
            fs%visc=visc+sgs%visc
         end block sgs_modeling
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Form implicit residuals
            !call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            
            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                        if (fs%vmask(i,j,k).eq.0) fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                        if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
            end block ibforcing
           
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         
      end do

   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)

   end subroutine simulation_final
   
end module simulation
