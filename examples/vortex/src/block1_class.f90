!> Definition for a block 1 simulation
module block1_class
   use string,            only: str_short,str_medium
   use precision,         only: WP
   use config_class,      only: config
   use incomp_class,     only: incomp
   use hypre_str_class,   only: hypre_str
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private

   public :: block1

   !> Block 1 object
   type :: block1
      class(config), pointer :: cfg       !< Pointer to config
      type(incomp) :: fs                  !< Single-phase incompressible flow solver
      type(hypre_str) :: ps               !< Unstructured HYPRE pressure solver
      type(hypre_str) :: is               !< Unstructured HYPRE implicit solver
      type(timetracker) :: time           !< Time tracker
      type(ensight) :: ens_out            !< Ensight output
      type(event) :: ens_evt              !< Ensight output event
      type(monitor) :: mfile,cflfile      !< Monitor files
      !> Private work arrays
      real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:), allocatable :: SR
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block1

   !> Vortex parameters
   real(WP) :: tau,beta,A

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW

   !> Fluid viscosity
   real(WP) :: visc
   
   !> Convective velocity array
   real(WP), dimension(3) :: Uc
   real(WP) :: Uinf
   real(WP) :: Vinf

   real(WP) :: vheight,vint,vradius
   real(WP) :: s

   !> Vortex center location
   real(WP) :: xc,yc,r_squared
   
   !> Orientation
   character(len=str_medium) :: orientation


contains



   !> Initialization of block 1
   subroutine init(b)
      use param, only: param_read
      implicit none
      class(block1), intent(inout) :: b

      ! Allocate work arrays for cfg
      allocate_work_arrays: block
         allocate(b%resU(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resV(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resW(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Ui  (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Vi  (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Wi  (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%SR(6,b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
      end block allocate_work_arrays
      

      ! Initialize time tracker
      initialize_timetracker: block
         b%time=timetracker(b%cfg%amRoot,name='block1')
         call param_read('Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker


      ! Create a single-phase flow solver with bconds
      create_solver: block
         use hypre_str_class, only: pcg_pfmg,pcg_smg,gmres_pfmg
         use mathtools,       only: Pi,twoPi 
         integer :: i,j,k
         real(WP) :: gamma=1.40_WP
         ! Create a single-phase flow solver
         b%fs=incomp(cfg=b%cfg,name='Block 1 NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); b%fs%visc=visc
         ! Initialize our vortex feild
         call param_read('Vortex size',tau)
         call param_read('Vortex strength',beta)
         call param_read('Convective Velocity',Uc)
         call param_read('x center',xc)
         call param_read('y center',yc)
         call param_read('Vortex height',vheight)
         call param_read('Vortex radius',vradius)
         call param_read('Vortex strength',vint)
         ! Initial convection velocity field 
         ! b%fs%U=Uc(1); b%fs%V=Uc(2); b%fs%W=Uc(3)
         ! Assign constant density
         call param_read('Density',b%fs%rho)
         ! Calculate density in each cell
         ! do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
         !    do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
         !       do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
         !             r_squared=(b%fs%cfg%xm(k)-xc)**2+(b%fs%cfg%ym(k)-yc)**2
         !             b%fs%rho(i,j,k)=(1-((((gamma-1.00_WP)*beta**2)/(8.00_WP*gamma*Pi))*exp(1.00_WP-r_squared)))**(1.00_WP/(gamma-1.00_WP))
         !       end do
         !    end do
         ! end do

         ! Prepare and configure pressure solver
         b%ps=hypre_str(cfg=b%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         b%ps%maxlevel=15
         call param_read('Pressure iteration',b%ps%maxit)
         call param_read('Pressure tolerance',b%ps%rcvg)
         ! Prepare and configure implicit solver
         !b%is=hypre_str(cfg=b%cfg,name='Implicit',method=pcg_pfmg,nst=7)
         b%is=hypre_str(cfg=b%cfg,name='Implicit',method=pcg_smg,nst=7)
         call param_read('Implicit iteration',b%is%maxit)
         call param_read('Implicit tolerance',b%is%rcvg)
         ! Setup the solver
         call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%is)
         ! Initial convection velocity field 
         ! b%fs%U(i,j,k)=Uc(1); b%fs%V(i,j,k)=Uc(2); b%fs%W(i,j,k)=Uc(3)
         ! A = sqrt(2.0_WP)*Umax*exp(0.5_WP)/tau
         call param_read('U infinity',Uinf)
         call param_read('V infinity',Vinf)
         meanU=Uinf
         A=beta*exp(0.5_WP)/(twoPi)
         do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
            do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
               do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
                  ! ! Main vortex
                  ! s=sqrt((b%cfg%xm(i)-vheight)**2+(b%cfg%ym(j)-vradius)**2)+tiny(1.0_WP)
                  ! b%fs%U(i,j,k)=b%fs%U(i,j,k)+1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%ym(j)-vradius)/s
                  ! b%fs%V(i,j,k)=-1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%xm(i)-vheight)/s
                  
                  r_squared=(b%fs%cfg%xm(i)-xc)**2+(b%fs%cfg%ym(j)-yc)**2
                  b%fs%U(i,j,k)=Uinf-A*exp(1.0_WP-r_squared)*(b%fs%cfg%ym(j)-yc)
                  b%fs%V(i,j,k)=Vinf+A*exp(1.0_WP-r_squared)*(b%fs%cfg%xm(i)-xc)
                  
                  ! ! Axial
                  ! s = sqrt((b%cfg%x(i)-vheight)**2+(sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2)-vradius)**2)+tiny(1.0_WP)
                  ! b%fs%U(i,j,k) = +1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2)-vradius)/s
                  ! ! V
                  ! s = sqrt((b%cfg%xm(i)-vheight)**2+(sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)-vradius)**2)+tiny(1.0_WP)
                  ! b%fs%V(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%xm(i)-vheight)/s*b%cfg%y(j)/sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)
                  ! W
                  ! s = sqrt((b%cfg%xm(i)-vheight)**2+(sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)-vradius)**2)+tiny(1.0_WP)
                  ! b%fs%W(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%xm(i)-vheight)/s*b%cfg%z(k)/sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)
                     ! r_squared=(b%fs%cfg%xm(k)-xc)**2+(b%fs%cfg%ym(k)-yc)**2
                     ! b%fs%U(i,j,k)=b%fs%U(i,j,k)-A*exp(1.0_WP-r_squared)*(b%fs%cfg%ym(j)-yc)
                     ! b%fs%V(i,j,k)=b%fs%V(i,j,k)+A*exp(1.0_WP-r_squared)*(b%fs%cfg%xm(i)-xc)
                     ! b%fs%U(i,j,k)=b%fs%U(i,j,k)+A*exp(-(b%cfg%x(i)**2+b%cfg%zm(k)**2)/tau**2)*(-b%cfg%zm(k))
                     ! b%fs%V(i,j,k)=b%fs%V(i,j,k)+A*exp(-(b%cfg%xm(i)**2+b%cfg%y(j)**2)/tau**2)*(b%cfg%xm(i))
                     ! b%fs%W(i,j,k)=b%fs%W(i,j,k)+A*exp(-(b%cfg%xm(i)**2+b%cfg%z(k)**2)/tau**2)*(b%cfg%xm(i))    
               end do
            end do
         end do
         ! Calculate cell-centered velocities and divergence
         call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
         call b%fs%get_div()
      end block create_solver

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(b%cfg,'block1')
         ! Create event for Ensight output
         b%ens_evt=event(b%time,'Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation1')
         call b%mfile%add_column(b%time%n,'Timestep number')
         call b%mfile%add_column(b%time%t,'Time')
         call b%mfile%add_column(b%time%dt,'Timestep size')
         call b%mfile%add_column(b%time%cfl,'Maximum CFL')
         call b%mfile%add_column(b%fs%Umax,'Umax')
         call b%mfile%add_column(b%fs%Vmax,'Vmax')
         call b%mfile%add_column(b%fs%Wmax,'Wmax')
         call b%mfile%add_column(b%fs%Pmax,'Pmax')
         call b%mfile%add_column(b%fs%divmax,'Maximum divergence')
         call b%mfile%add_column(b%fs%psolv%it,'Pressure iteration')
         call b%mfile%add_column(b%fs%psolv%rerr,'Pressure error')
         call b%mfile%write()
         ! Create CFL monitor
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl1')
         call b%cflfile%add_column(b%time%n,'Timestep number')
         call b%cflfile%add_column(b%time%t,'Time')
         call b%cflfile%add_column(b%fs%CFLc_x,'Convective xCFL')
         call b%cflfile%add_column(b%fs%CFLc_y,'Convective yCFL')
         call b%cflfile%add_column(b%fs%CFLc_z,'Convective zCFL')
         call b%cflfile%add_column(b%fs%CFLv_x,'Viscous xCFL')
         call b%cflfile%add_column(b%fs%CFLv_y,'Viscous yCFL')
         call b%cflfile%add_column(b%fs%CFLv_z,'Viscous zCFL')
         call b%cflfile%write()
      end block create_monitor

   end subroutine init


   !> Take a time step with block 1
   subroutine step(b)
      implicit none
      class(block1), intent(inout) :: b

      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! ! Convect vortex through the domain
      ! convect_vortex: block
      !    use mathtools,       only: Pi,twoPi 
      !    integer :: i,j,k
      !    do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
      !       do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
      !          do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_

      !             ! r_squared=(b%fs%cfg%xm(i)-xc)**2+(b%fs%cfg%ym(j)-yc)**2
      !             ! b%fs%U(i,j,k)=b%fs%U(i,j,k)-A*exp(1.0_WP-r_squared)*(b%fs%cfg%ym(j)-yc)
      !             ! b%fs%V(i,j,k)=b%fs%V(i,j,k)+A*exp(1.0_WP-r_squared)*(b%fs%cfg%xm(i)-xc)

      !             ! r_squared=(b%fs%cfg%xm(k)-xc)**2+(b%fs%cfg%ym(k)-yc)**2
      !             ! b%fs%U(i,j,k)=b%fs%U(i,j,k)-A*exp(1.0_WP-r_squared)*(b%fs%cfg%ym(j)-yc)
      !             ! b%fs%V(i,j,k)=b%fs%V(i,j,k)+A*exp(1.0_WP-r_squared)*(b%fs%cfg%xm(i)-xc)
      !             ! ! Main vortex
      !             ! s=sqrt((b%cfg%xm(i)-vheight)**2+(b%cfg%ym(j)-vradius)**2)+tiny(1.0_WP)
      !             ! b%fs%U(i,j,k)=1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%ym(j)-vradius)/s
      !             ! b%fs%V(i,j,k)=-1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%xm(i)-vheight)/s
      !             ! ! Axial
      !             ! s = sqrt((b%cfg%x(i)-vheight)**2+(sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2)-vradius)**2)+tiny(1.0_WP)
      !             ! b%fs%U(i,j,k) = +1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2)-vradius)/s
      !             ! V
      !             ! s = sqrt((b%cfg%xm(i)-vheight)**2+(sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)-vradius)**2)+tiny(1.0_WP)
      !             ! b%fs%V(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%xm(i)-vheight)/s*b%cfg%y(j)/sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)
      !             ! ! W
      !             ! s = sqrt((b%cfg%xm(i)-vheight)**2+(sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)-vradius)**2)+tiny(1.0_WP)
      !             ! b%fs%W(i,j,k) = -1.0_WP/(twoPi*s)*(1.0_WP-exp(-s**2/vint**2))*(b%cfg%xm(i)-vheight)/s*b%cfg%z(k)/sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)
      !                ! r_squared=(b%cfg%xm(k)-xc)**2+(b%cfg%ym(k)-yc)**2
      !                ! b%fs%U(i,j,k)=b%fs%U(i,j,k)-A*exp(1.0_WP-r_squared)*(b%cfg%ym(j)-yc)
      !                ! b%fs%V(i,j,k)=b%fs%V(i,j,k)+A*exp(1.0_WP-r_squared)*(b%cfg%xm(i)-xc)
      !             ! b%fs%U(i,j,k)=b%fs%U(i,j,k)+A*exp(-(b%cfg%x(i)**2+b%cfg%zm(k)**2)/tau**2)*(-b%cfg%zm(k))
      !             ! b%fs%V(i,j,k)=b%fs%V(i,j,k)+A*exp(-(b%cfg%xm(i)**2+b%cfg%y(j)**2)/tau**2)*(b%cfg%xm(i))
      !             ! b%fs%W(i,j,k)=b%fs%W(i,j,k)+A*exp(-(b%cfg%xm(i)**2+b%cfg%z(k)**2)/tau**2)*(b%cfg%xm(i))  
      !          end do
      !       end do
      !    end do
      ! end block convect_vortex

      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W

      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)

         ! Build mid-time velocity
         b%fs%U=0.5_WP*(b%fs%U+b%fs%Uold)
         b%fs%V=0.5_WP*(b%fs%V+b%fs%Vold)
         b%fs%W=0.5_WP*(b%fs%W+b%fs%Wold)

         ! Explicit calculation of drho*u/dt from NS
         call b%fs%get_dmomdt(b%resU,b%resV,b%resW)

         ! Assemble explicit residual
         b%resU=-2.0_WP*(b%fs%rho*b%fs%U-b%fs%rho*b%fs%Uold)+b%time%dt*b%resU
         b%resV=-2.0_WP*(b%fs%rho*b%fs%V-b%fs%rho*b%fs%Vold)+b%time%dt*b%resV
         b%resW=-2.0_WP*(b%fs%rho*b%fs%W-b%fs%rho*b%fs%Wold)+b%time%dt*b%resW

         ! Add body forcing
         forcing: block
            use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
            use parallel, only: MPI_REAL_WP
            integer :: i,j,k,ierr
            real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
            myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
            do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
               do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
                  do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                     if (b%fs%umask(i,j,k).eq.0) then
                        myU   =myU   +b%fs%cfg%dxm(i)*b%fs%cfg%dy(j)*b%fs%cfg%dz(k)*(2.0_WP*b%fs%U(i,j,k)-b%fs%Uold(i,j,k))
                        myUvol=myUvol+b%fs%cfg%dxm(i)*b%fs%cfg%dy(j)*b%fs%cfg%dz(k)
                     end if
                  end do
               end do
            end do
            call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,b%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,b%fs%cfg%comm,ierr); meanU=meanU/Uvol
            where (b%fs%umask.eq.0) b%resU=b%resU+Ubulk-meanU
         end block forcing 

         ! Form implicit residuals
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)

         ! Apply these residuals
         b%fs%U=2.0_WP*b%fs%U-b%fs%Uold+b%resU
         b%fs%V=2.0_WP*b%fs%V-b%fs%Vold+b%resV
         b%fs%W=2.0_WP*b%fs%W-b%fs%Wold+b%resW

         ! Apply other boundary conditions
         call b%fs%apply_bcond(b%time%t,b%time%dt)

         ! Solve Poisson equation
         call b%fs%correct_mfr()
         call b%fs%get_div()
         b%ps%rhs=-b%fs%cfg%vol*b%fs%div*b%fs%rho/b%time%dt
         b%ps%sol=0.0_WP
         call b%ps%solve()

         ! Correct velocity
         call b%fs%get_pgrad(b%ps%sol,b%resU,b%resV,b%resW)
         b%fs%P=b%fs%P+b%ps%sol
         b%fs%U=b%fs%U-b%time%dt*b%resU/b%fs%rho
         b%fs%V=b%fs%V-b%time%dt*b%resV/b%fs%rho
         b%fs%W=b%fs%W-b%time%dt*b%resW/b%fs%rho

         ! Increment sub-iteration counter
         b%time%it=b%time%it+1

      end do

      ! Recompute interpolated velocity and divergence
      call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
      call b%fs%get_div()

      ! Output to ensight
      if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)

      ! Perform and output monitoring
      call b%fs%get_max()
      call b%mfile%write()
      call b%cflfile%write()

   end subroutine step


   !> Finalize b1 simulation
   subroutine final(b)
      implicit none
      class(block1), intent(inout) :: b

      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR)

   end subroutine final


end module block1_class
