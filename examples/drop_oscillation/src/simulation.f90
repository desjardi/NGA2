!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use mathtools,         only: Pi
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use irl_fortran_interface
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps,vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,curvfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Problem definition
   real(WP), dimension(3) :: center,vel
   real(WP) :: radius
   integer, dimension(3) :: center_ind

   !> Max velocity magnitude
   real(WP) :: r0
   
   !> Curvature errors
   real(WP) :: l2_error,linf_error,curv_rms

contains
   
   
   !> Function that defines a level set function for a cylinder
   function levelset_circle(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      real(WP) :: epsilon, theta, costheta, R
      epsilon = 0.08_WP
      center = [0.0_WP,0.0_WP,0.0_WP]
      ! Create droplet
      theta=atan2(xyz(2)-center(2),xyz(1)-center(1))
      costheta=cos(theta)
      ! G=radius*(1.0_WP+epsilon*0.5_WP*(3.0_WP*costheta**2-1.0_WP)-epsilon**2/5.0_WP)-sqrt(sum((xyz-center)**2))
      G=radius*(1.0_WP+epsilon*0.5_WP*cos(12.0_WP*theta))-sqrt(sum((xyz-center)**2))
   end function levelset_circle

   
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
         use mms_geom, only: cube_refine_vol
         use vfs_class,only: r2p,lvira,VFhi,VFlo
         use mpi_f08,  only: MPI_WTIME
         use string,   only: str_medium,lowercase
         integer :: i,j,k,n,si,sj,sk,curvature_method,stencil_size,hf_backup_method
         character(len=str_medium) :: read_curvature_method
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=5
         real(WP) :: start, finish
         ! Create a VOF solver with r2p reconstruction
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize two droplets
         call param_read('Droplet diameter',radius); radius=0.5_WP*radius
         ! Assign a random center
         random_center: block
            use random,   only: random_uniform
            use mpi_f08,  only: MPI_BCAST
            use parallel, only: MPI_REAL_WP
            integer :: ierr
            center_ind=cfg%get_ijk_global([0.0_WP,0.0_WP,0.0_WP],[vf%cfg%imino_,vf%cfg%jmino_,vf%cfg%kmino_])
            ! if (vf%cfg%amroot) center=[random_uniform(-vf%cfg%dx(center_ind(1))/2.0_WP,vf%cfg%dx(center_ind(1))/2.0_WP),&
            ! &       random_uniform(-vf%cfg%dy(center_ind(2))/2.0_WP,vf%cfg%dy(center_ind(2))/2.0_WP),&
            ! &       random_uniform(-vf%cfg%dz(center_ind(3))/2.0_WP,vf%cfg%dz(center_ind(3))/2.0_WP)]   
            call MPI_BCAST(center,3,MPI_REAL_WP,0,vf%cfg%comm,ierr)
         end block random_center
         ! call param_read('Droplet 1 position',center)      
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_circle,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
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
         use hypre_str_class, only: pcg_pfmg
         use mathtools, only: Pi
         integer :: i,j,k
         real(WP), dimension(3) :: xyz
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initial velocity field
         call param_read('Droplet velocity',vel)
         fs%U=vel(1)
         fs%V=vel(2)
         fs%W=vel(3)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=2,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         smesh%var(1,:)=1.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
         call add_surfgrid_variable(smesh,2,vf%curv)      
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='DropOscillation')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',smesh) 
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call get_radius()
         ! call print_curvature_error()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(r0,'Radius')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         ! call mfile%add_column(fs%implicit%it,'Velocity iteration')
         ! call mfile%add_column(fs%implicit%rerr,'Velocity error')         
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%wt,'Wallclock time')
         call cflfile%add_column(time%dt,'Timestep size')
         call cflfile%add_column(time%cfl,'Maximum CFL')         
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create curvature monitor
         curvfile=monitor(fs%cfg%amroot,'curvature')
         call curvfile%add_column(time%n,'Timestep number')
         call curvfile%add_column(time%t,'Time')
         call curvfile%add_column(l2_error,'L2 error')
         call curvfile%add_column(linf_error,'Linf error')
         ! call curvfile%add_column(curv_rms,'RMS curv')
         ! call curvfile%add_column(vf%hf_failure_rate,'HF failure rate')
         ! call curvfile%add_column(vf%hf_dir1_rate,'HF direction 1 rate')
         ! call curvfile%add_column(vf%hf_dir2_rate,'HF direction 2 rate')
         ! call curvfile%add_column(vf%hf_dir3_rate,'HF direction 3 rate')
         call curvfile%write()         
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
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
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
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
         if (ens_evt%occurs()) then
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Also populate nplane variable
               smesh%var(1,:)=1.0_WP
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                           end if
                        end do
                     end do
                  end do
               end do
               ! Populate curv variable
               call add_surfgrid_variable(smesh,2,vf%curv)      
            end block update_smesh           
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call get_radius()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         call print_curvature_error()
         call curvfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      call print_curvature_error()
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   

   !> Print curvature error
   subroutine print_curvature_error
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_SUM,MPI_INTEGER
      use parallel, only: MPI_REAL_WP      
      use irl_fortran_interface, only: getNumberOfVertices
      implicit none
      real(WP) :: my_l2_error,my_linf_error,err,curv_ref,my_curv_sum,curv_sum,curv_mean,my_rms_sum
      integer :: my_count,count,ierr
      integer :: i,j,k

      my_l2_error=0.0_WP;my_linf_error=0.0_WP;my_count=0;my_curv_sum=0.0_WP
      curv_ref = 1.0_WP/radius
      if (vf%cfg%nz.gt.1) curv_ref = curv_ref*2.0_WP ! 3D
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (getNumberOfVertices(vf%interface_polygon(1,i,j,k)).eq.0) cycle
               ! print *,'ijk',i,j,k,'curv',vf%curv(i,j,k)
               err=(vf%curv(i,j,k)-curv_ref)/curv_ref
               my_l2_error=my_l2_error + (err)**2
               my_linf_error=max(my_linf_error,err)                     
               my_count=my_count+1
               my_curv_sum=my_curv_sum+vf%curv(i,j,k)
            end do
         end do
      end do

      ! Get the parallel max
      call MPI_ALLREDUCE(my_linf_error,linf_error,1,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_l2_error  ,l2_error  ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_count     ,count     ,1,MPI_INTEGER,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_curv_sum  ,curv_sum  ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)

      l2_error=sqrt(l2_error/real(count,WP))
      ! if (cfg%amroot) print *,l2_error,linf_error

      ! my_l2_error=sqrt(my_l2_error/real(my_count,WP))
      ! print *, my_l2_error,my_linf_error

      curv_mean=curv_sum/count
      my_rms_sum=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (getNumberOfVertices(vf%interface_polygon(1,i,j,k)).eq.0) cycle
               err=(vf%curv(i,j,k)-curv_mean)/curv_ref
               my_rms_sum=my_rms_sum+ (err)**2
            end do
         end do
      end do

      call MPI_ALLREDUCE(my_rms_sum  ,curv_rms  ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      curv_rms=sqrt(my_rms_sum/real(count,WP))

   end subroutine print_curvature_error


   !> Make a surface scalar variable
   subroutine add_surfgrid_variable(smesh,var_index,A)
      use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
      implicit none
      class(surfmesh), intent(inout) :: smesh
      integer, intent(in) :: var_index
      real(WP), dimension(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_), intent(in) :: A 
      integer :: i,j,k,shape,np,nplane
      
      ! Fill out arrays
      if (smesh%nPoly.gt.0) then
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(vf%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        ! Set nplane variable
                        smesh%var(var_index,np)=A(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      else
         smesh%var(var_index,1)=1
      end if      
      
   end subroutine add_surfgrid_variable

   !> Calculate the max velocity magnitude relative to the prescribed velocity
   subroutine get_radius
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr,ncells,nplane,shape,n,d
      real(WP) :: myR,myminR,mymaxR,minR,maxR
      real(WP), dimension(3) :: tmp_vert

      ! Set all to zero
      myR=0.0_WP;myminR=0.0_WP;mymaxR=0.0_WP
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         j=1+(fs%cfg%ny)/2
         i=1+(fs%cfg%nx)/2
         if (i.ge.fs%cfg%imin_.and.i.le.fs%cfg%imax_.and.j.ge.fs%cfg%jmin_.and.j.le.fs%cfg%jmax_) then
            myR=myR+(fs%cfg%z(k+1)-fs%cfg%z(k))*vf%VF(i,j,k)
            do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
               shape=getNumberOfVertices(vf%interface_polygon(nplane,i,j,k))
               if (shape.gt.0) then
                  do n=1,shape
                     tmp_vert=getPt(vf%interface_polygon(nplane,i,j,k),n-1)
                     myminR=max(myminR,-tmp_vert(3))
                     mymaxR=max(mymaxR,tmp_vert(3))
                  end do
               end if
            end do
         end if
      end do
      
      ! Get the parallel max
      call MPI_ALLREDUCE(myR,r0,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myminR,minR,1,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(mymaxR,maxR,1,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr)
      r0=maxR+minR

   end subroutine get_radius


end module simulation
