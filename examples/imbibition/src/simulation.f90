!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Problem definition and post-processing
   real(WP), dimension(3) :: Cdrop
   real(WP) :: Rdrop,Vimb
   type(monitor) :: ppfile
   type(event) :: ppevt
   
contains
   
   
   !> Function that defines a level set function for a contacting drop problem
   function levelset_contact_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the droplet
      G=Rdrop-sqrt(sum((xyz-Cdrop)**2))
   end function levelset_contact_drop
   
   
   !> Specialized subroutine that computes the imbibed liquid volume
   subroutine get_Vimb()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr
      real(WP) :: my_Vimb
      my_Vimb=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (vf%cfg%ym(j).lt.0.0_WP) my_Vimb=my_Vimb+vf%VF(i,j,k)*vf%cfg%vol(i,j,k)*vf%cfg%VF(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_Vimb,Vimb,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
   end subroutine get_Vimb
   
   
   !> Specialized subroutine that outputs the vertical liquid distribution
   subroutine postproc_data()
      use mathtools, only: Pi
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: localvol,totalvol
      real(WP), dimension(:), allocatable :: localvel,totalvel
      real(WP), dimension(:), allocatable :: localpre,totalpre
      character(len=str_medium) :: filename,timestamp
      ! Allocate vertical line storage
      allocate(localvol(vf%cfg%jmin:vf%cfg%jmax)); localvol=0.0_WP
      allocate(localvel(vf%cfg%jmin:vf%cfg%jmax)); localvel=0.0_WP
      allocate(localpre(vf%cfg%jmin:vf%cfg%jmax)); localpre=0.0_WP
      allocate(totalvol(vf%cfg%jmin:vf%cfg%jmax)); totalvol=0.0_WP
      allocate(totalvel(vf%cfg%jmin:vf%cfg%jmax)); totalvel=0.0_WP
      allocate(totalpre(vf%cfg%jmin:vf%cfg%jmax)); totalpre=0.0_WP
      ! Initialize local data to zero
      localvol=0.0_WP; localvel=0.0_WP; localpre=0.0_WP
      ! Integrate all data over x and z
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               localvol(j)=localvol(j)+vf%VF(i,j,k)*vf%cfg%dx(i)*vf%cfg%dz(k)
               localvel(j)=localvel(j)+vf%VF(i,j,k)*vf%cfg%dx(i)*vf%cfg%dz(k)*Vi(i,j,k)
               localpre(j)=localpre(j)+vf%VF(i,j,k)*vf%cfg%dx(i)*vf%cfg%dz(k)*fs%P(i,j,k)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE(localvol,totalvol,vf%cfg%ny,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(localvel,totalvel,vf%cfg%ny,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(localpre,totalpre,vf%cfg%ny,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      ! If root, print it out
      if (vf%cfg%amRoot) then
         call execute_command_line('mkdir -p radius')
         filename='radius_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file='radius/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'Height','Liq_Radius','Liq_Area','Liq_VFR','Liq_P'
         do j=vf%cfg%jmin,vf%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') vf%cfg%ym(j),sqrt(totalvol(j)/Pi),totalvol(j),totalvel(j),totalpre(j)/(totalvol(j)+tiny(totalvol(j)))
         end do
         close(iunit)
      end if
      ! Deallocate work arrays
      deallocate(localvol,totalvol)
      deallocate(localvel,totalvel)
      deallocate(localpre,totalpre)
   end subroutine postproc_data
   
   
   !> Function that localizes the top (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function xp_locator
   
   
   !> Function that localizes the bottom (x-) of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function xm_locator
   
   
   !> Function that localizes the top (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function zp_locator
   
   
   !> Function that localizes the bottom (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin.and.pg%ym(j).ge.0.0_WP) isIn=.true.
   end function zm_locator
   
   
   !> Function that localizes the needle injection
   function needle(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1.and.sqrt(pg%xm(i)**2+pg%zm(k)**2).le.0.70e-3_WP) isIn=.true.
   end function needle
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo
         use mathtools, only: Pi
         use, intrinsic :: iso_fortran_env, only: output_unit
         use string,    only: str_long
         use messager,  only: log
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area,contact
         integer, parameter :: amr_ref_lvl=4
         character(len=str_long) :: message
         ! Create a VOF solver
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Prepare the analytical calculation of a sphere on a wall
         call param_read('Droplet radius',Rdrop)
         call param_read('Initial contact angle',contact); contact=contact*Pi/180.0_WP
         if (vf%cfg%nz.eq.1) then ! 2D analytical drop shape
            Rdrop=Rdrop*sqrt(Pi/(2.0_WP*(contact-sin(contact)*cos(contact))))
         else ! 3D analytical drop shape
            Rdrop=Rdrop*(4.0_WP/(2.0_WP-3.0_WP*cos(contact)+(cos(contact))**3))**(1.0_WP/3.0_WP)
         end if
         Cdrop=[0.0_WP,-Rdrop*cos(contact),0.0_WP]
         if (vf%cfg%amRoot) then
            write(output_unit,'("Droplet initial radius is ",es12.5)') Rdrop
            write(message    ,'("Droplet initial radius is ",es12.5)') Rdrop; call log(message)
         end if
         ! Initialize the VOF field
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Fill out part of the needle
                  if (vf%cfg%ym(j).gt.0.013_WP.and.sqrt(vf%cfg%xm(i)**2+vf%cfg%zm(k)**2).lt.0.0007_WP) then
                     vf%VF(i,j,k)=1.0_WP
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  else
                     vf%VF(i,j,k)=0.0_WP
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
                  ! ! Handle wall cells or cells below the plate surface
                  ! if (vf%mask(i,j,k).eq.1.or.vf%cfg%ym(j).lt.0.0_WP) then
                  !    vf%VF(i,j,k)=0.0_WP
                  !    vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  !    vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  !    cycle
                  ! end if
                  ! ! Fill out the needle - this is handwavy...
                  ! if (vf%cfg%ym(j).gt.0.010_WP.and.sqrt(vf%cfg%xm(i)**2+vf%cfg%zm(k)**2).lt.0.0007_WP) then
                  !    vf%VF(i,j,k)=1.0_WP
                  !    vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  !    vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  !    cycle
                  ! end if
                  ! ! Set cube vertices
                  ! n=0
                  ! do sk=0,1
                  !    do sj=0,1
                  !       do si=0,1
                  !          n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                  !       end do
                  !    end do
                  ! end do
                  ! ! Call adaptive refinement code to get volume and barycenters recursively
                  ! vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  ! call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_contact_drop,0.0_WP,amr_ref_lvl)
                  ! vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  ! if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                  !    vf%Lbary(:,i,j,k)=v_cent
                  !    vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  ! else
                  !    vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  !    vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  ! end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class, only: clipped_neumann,dirichlet,bcond
         use ils_class,  only: pcg_pfmg
         use mathtools,  only: Pi
         type(bcond), pointer :: mybc
         real(WP) :: Vneedle
         integer :: i,j,k,n
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient and contact angle
         call param_read('Surface tension coefficient',fs%sigma)
         call param_read('Static contact angle',fs%contact_angle)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Setup boundary conditions
         call fs%add_bcond(name='bc_xp' ,type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         call fs%add_bcond(name='bc_xm' ,type=clipped_neumann,face='x',dir=-1,canCorrect=.true. ,locator=xm_locator)
         call fs%add_bcond(name='bc_zp' ,type=clipped_neumann,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call fs%add_bcond(name='bc_zm' ,type=clipped_neumann,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         call fs%add_bcond(name='needle',type=dirichlet      ,face='y',dir=+1,canCorrect=.false.,locator=needle    )
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         fs%psolv%maxlevel=10
         call fs%setup(pressure_ils=pcg_pfmg,implicit_ils=pcg_pfmg)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply Dirichlet at liquid needle
         call param_read('Liquid injection velocity',Vneedle)
         call fs%get_bcond('needle',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=-Vneedle
         end do
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='twophase')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',vf%surfgrid)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
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
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
      ! Create a specialized post-processing file
      create_postproc: block
         ! Calculate imbibed volume
         call get_Vimb()
         ! Create monitor for liquid data
         ppfile=monitor(fs%cfg%amRoot,'dropinfo')
         call ppfile%add_column(time%n,'Timestep number')
         call ppfile%add_column(time%t,'Time')
         call ppfile%add_column(vf%VFmax,'VOF maximum')
         call ppfile%add_column(vf%VFmin,'VOF minimum')
         call ppfile%add_column(vf%VFint,'Total volume')
         call ppfile%add_column(Vimb,'Imbibed volume')
         call ppfile%write()
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Perform the output
         if (ppevt%occurs()) call postproc_data()
      end block create_postproc
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
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
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         
         ! Specialized post-processing
         call get_Vimb(); call ppfile%write()
         if (ppevt%occurs()) call postproc_data()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
