!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,bedbottom,bedheight,bedwidth
   use lpt_class,         only: lpt,part,MPI_PART_SIZE,MPI_PART
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Two-phase flow solver, volume fraction solver, LPT solver, time tracker
   type(lpt),         public :: lp
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: sfile,pfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: srcU,srcV,srcW,srcM
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Ui_old,Vi_old,Wi_old
   
   !> Problem definition
   real(WP) :: depth
   
   !> Max timestep size for LPT
   real(WP) :: lp_dt,lp_dt_max

   !> Particle bed injection
   integer :: nbed
   type(part), dimension(:), allocatable :: pbed
   real(WP) :: vbed
   logical :: freebed
   
   
contains
   
   
   !> Function that defines a level set function for a pool of water
   function levelset_pool(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Add the pool
      G=depth-xyz(2)
   end function levelset_pool
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(srcU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcM(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
      
      
      ! Initialize our LPT
      initialize_lpt: block
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Set gravity
         call param_read('Gravity',lp%gravity)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Set drag model
         lp%drag_model='Tavanashad'
         ! Initialize with zero particles
         call lp%resize(0)
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Maximum timestep size used for particles
         call param_read('Particle timestep size',lp_dt_max,default=huge(1.0_WP))
         lp_dt=lp_dt_max
         ! Set collision timescale
         call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
         ! Set coefficient of restitution
         call param_read('Coefficient of restitution',lp%e_n)
         call param_read('Wall restitution',lp%e_w)
         call param_read('Friction coefficient',lp%mu_f)
      end block initialize_lpt
      
      
      ! Initialize particles - precalculated bed entry case
      !prepare_bed: block
      !   use mpi_f08
      !   use messager,  only: die
      !   use quicksort, only: quick_sort
      !   use string,    only: str_medium
      !  use parallel,  only: comm,info_mpiio
      !   type(MPI_File) :: ifile
      !   type(MPI_Status):: status
      !   integer :: ierr,psize,i
      !   integer(kind=MPI_OFFSET_KIND) :: offset
      !   character(len=str_medium) :: bedfile
      !   real(WP), dimension(:), allocatable :: pos
      !   integer , dimension(:), allocatable :: id
      !   type(part), dimension(:), allocatable :: tmp
      !   real(WP) :: maxpos
      !   ! Injection parameters
      !   call param_read('Bed file to inject',bedfile)
      !   call param_read('Bed injection velocity',vbed)
      !   call param_read('Bed is free',freebed)
      !   ! Read the header of the particle file
      !   call MPI_FILE_OPEN(comm,trim(bedfile),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      !   if (ierr.ne.0) call die('[impinging_particles] Problem encountered while reading bed file: '//trim(bedfile))
      !   call MPI_FILE_READ_ALL(ifile,nbed,1,MPI_INTEGER,status,ierr)
      !   call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)
      !   if (psize.ne.MPI_PART_SIZE) call die('[impinging_particles] Particle type unreadable')
      !   ! Root reads in the full file and stores the particles
      !   if (lp%cfg%amRoot) then
      !      allocate(pbed(nbed))
      !      call MPI_FILE_GET_POSITION(ifile,offset,ierr)
      !      call MPI_FILE_READ_AT(ifile,offset,pbed,nbed,MPI_PART,status,ierr)
      !   end if
      !   ! Close the bed file
      !   call MPI_FILE_CLOSE(ifile,ierr)
      !   ! Manipulate the bed particles a bit
      !   if (lp%cfg%amRoot) then
      !      ! Sort particles
      !      allocate(pos(nbed),id(nbed))
      !      do i=1,nbed
      !         id(i)=i
      !         pos(i)=pbed(i)%pos(2)
      !      end do
      !      ! Flip the bed
      !      maxpos=maxval(pos)
      !      pos=maxpos-pos
      !      ! Apply quicksort
      !      call quick_sort(A=pos,B=id)
      !      ! Loop through particles and sort them
      !      allocate(tmp(nbed))
      !      do i=1,nbed
      !         tmp(i)=pbed(id(i))
      !         if (freebed) then
      !            tmp(i)%id=i
      !         else
      !            tmp(i)%id=-2
      !         end if
      !         ! Adjust position
      !         tmp(i)%pos(2)=pos(i)+lp%cfg%yL
      !         ! Adjust velocity
      !         tmp(i)%vel=[0.0_WP,vbed,0.0_WP]
      !         ! Zero out collisions
      !         tmp(i)%Acol=0.0_WP
      !         tmp(i)%Tcol=0.0_WP
      !         ! Zero out angular velocity
      !         tmp(i)%angVel=0.0_WP
      !      end do
      !      ! Store the particles
      !      pbed=tmp
      !      deallocate(tmp,pos,id)
      !   end if
      !end block prepare_bed
      
      
      ! Initialize particles using precalculated static bed
      prepare_bed: block
         use string, only: str_medium
         character(len=str_medium) :: bedfile
         integer :: i
         ! Bed file to use
         call param_read('Bed file to read',bedfile)
         ! Fixed bed
         call param_read('Bed is free',freebed,default=.true.)
         if (.not.freebed) call param_read('Bed velocity',vbed)
         ! Read it in
         call lp%read(filename=trim(bedfile))
         ! Loop through particles
         do i=1,lp%np_
            ! Shift bed position
            lp%p(i)%pos(1)=lp%p(i)%pos(1)
            lp%p(i)%pos(2)=lp%p(i)%pos(2)+bedbottom
            ! Zero out velocity
            lp%p(i)%vel=0.0_WP
            ! Handle fixed bed case
            if (.not.freebed) then
               lp%p(i)%id=-2
               lp%p(i)%vel=[0.0_WP,-vbed,0.0_WP]
            end if
            ! Clip the top
            if (lp%p(i)%pos(2).gt.bedbottom+bedheight) lp%p(i)%flag=1
            ! Relocate particles
            lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,lp%p(i)%ind)
         end do
         call lp%sync()
         ! Recalculate VF
         call lp%update_VF()
      end block prepare_bed

      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,plicnet,VFhi,VFlo,remap
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=remap,name='VOF')
         ! TEST THE REMOVAL OF WALL TREATMENT IN VOF TRANSPORT
         vf%vmask=0
         ! Initialize to a pool
         call param_read('Pool depth',depth)
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_pool,0.0_WP,amr_ref_lvl)
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
         use mathtools,       only: Pi
         use tpns_class,      only: clipped_neumann,dirichlet
         use hypre_str_class, only: pcg_pfmg2
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
         ! Transfer to lpt for surface tension modeling
         lp%sigma=fs%sigma
         lp%contact_angle=fs%contact_angle
         lp%VFst=0.25_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Outlet on the top
         call fs%add_bcond(name='outlet',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=1,nvec=1,name='lpt')
         pmesh%varname(1)='radius'
         pmesh%vecname(1)='velocity'
         call lp%update_partmesh(pmesh)
         do i=1,lp%np_
            pmesh%var(1,i)=0.5_WP*lp%p(i)%d
            pmesh%vec(:,1,i)=lp%p(i)%vel
         end do
      end block create_pmesh
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh

      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=fs%cfg,name='impinging_particles')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curv',vf%curv)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_surface('plic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call lp%get_max()
         ! Create simulation monitor
         sfile=monitor(fs%cfg%amRoot,'simulation')
         call sfile%add_column(time%n,'Timestep number')
         call sfile%add_column(time%t,'Time')
         call sfile%add_column(time%dt,'Timestep size')
         call sfile%add_column(time%cfl,'Maximum CFL')
         call sfile%add_column(fs%Umax,'Umax')
         call sfile%add_column(fs%Vmax,'Vmax')
         call sfile%add_column(fs%Wmax,'Wmax')
         call sfile%add_column(fs%Pmax,'Pmax')
         call sfile%add_column(vf%VFmax,'VOF maximum')
         call sfile%add_column(vf%VFmin,'VOF minimum')
         call sfile%add_column(vf%VFint,'VOF integral')
         call sfile%add_column(fs%divmax,'Maximum divergence')
         call sfile%add_column(fs%psolv%it,'Pressure iteration')
         call sfile%add_column(fs%psolv%rerr,'Pressure error')
         call sfile%write()
         ! Create particle monitor
         pfile=monitor(amroot=lp%cfg%amRoot,name='particles')
         call pfile%add_column(time%n,'Timestep number')
         call pfile%add_column(time%t,'Time')
         call pfile%add_column(lp_dt,'Particle dt')
         call pfile%add_column(lp%np,'Particle number')
         call pfile%add_column(lp%np_new,'Npart new')
         call pfile%add_column(lp%np_out,'Npart removed')
         call pfile%add_column(lp%ncol,'Particle collisions')
         call pfile%add_column(lp%VFmax,'Max VF')
         call pfile%add_column(lp%Umin,'Particle Umin')
         call pfile%add_column(lp%Umax,'Particle Umax')
         call pfile%add_column(lp%Vmin,'Particle Vmin')
         call pfile%add_column(lp%Vmax,'Particle Vmax')
         call pfile%add_column(lp%Wmin,'Particle Wmin')
         call pfile%add_column(lp%Wmax,'Particle Wmax')
         call pfile%add_column(lp%dmin,'Particle dmin')
         call pfile%add_column(lp%dmax,'Particle dmax')
         call pfile%write()
         ! Create CFL monitor
         cflfile=monitor(lp%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(lp%CFLp_x,'Particle xCFL')
         call cflfile%add_column(lp%CFLp_y,'Particle yCFL')
         call cflfile%add_column(lp%CFLp_z,'Particle zCFL')
         call cflfile%add_column(lp%CFL_col,'Collision CFL')
         call cflfile%write()
      end block create_monitor
      

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
      implicit none
      real(WP) :: cfl
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old Ui
         store_old_vel: block
            resU=fs%U; resV=fs%V; resW=fs%W
            fs%U=fs%Uold; fs%V=fs%Vold; fs%W=fs%Wold
            call fs%interp_vel(Ui_old,Vi_old,Wi_old)
            fs%U=resU; fs%V=resV; fs%W=resW
         end block store_old_vel
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Particle update
         lpt_step: block
            real(WP) :: dt_done,mydt
            real(WP), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3,VFold,rho
            real(WP), dimension(:,:,:), allocatable :: dVFdx,dVFdy,dVFdz
            real(WP), dimension(:,:,:,:), allocatable :: vort,acc
            integer :: i,j,k
            ! Allocate and store rho
            allocate(rho  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); rho  =vf%VF*fs%rho_l+(1.0_WP-vf%VF)*fs%rho_g
            ! Allocate and store VFold
            allocate(VFold(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFold=lp%VF
            ! Allocate src storage
            allocate(tmp1 (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp1 =0.0_WP
            allocate(tmp2 (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp2 =0.0_WP
            allocate(tmp3 (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); tmp3 =0.0_WP
            ! Allocate vorticity and ugradu storage
            allocate(vort(1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); vort=0.0_WP
            allocate(acc (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); acc =0.0_WP
            ! Allocate grad(VF) storage
            allocate(dVFdx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); dVFdx=0.0_WP
            allocate(dVFdy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); dVFdy=0.0_WP
            allocate(dVFdz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); dVFdz=0.0_WP
            do k=vf%cfg%kmin_,vf%cfg%kmax_; do j=vf%cfg%jmin_,vf%cfg%jmax_; do i=vf%cfg%imin_,vf%cfg%imax_
               dVFdx(i,j,k)=sum(fs%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k))
               dVFdy(i,j,k)=sum(fs%divu_y(:,i,j,k)*vf%VF(i,j-1:j,k))
               dVFdz(i,j,k)=sum(fs%divu_z(:,i,j,k)*vf%VF(i,j,k-1:k))
            end do; end do; end do
            call fs%cfg%sync(dVFdx)
            call fs%cfg%sync(dVFdy)
            call fs%cfg%sync(dVFdz)
            ! Get fluid stress
            resU=0.0_WP; resV=0.0_WP; resW=0.0_WP; call fs%get_div_stress(resU,resV,resW)
            ! Get vorticity
            call fs%get_vorticity(vort)
            ! Get fluid acceleration
            call fs%get_ugradu(acc)
            call fs%interp_vel(Ui,Vi,Wi)
            if (time%dtold.gt.0.0_WP) then
               acc(1,:,:,:)=acc(1,:,:,:)+(Ui-Ui_old)/time%dtold
               acc(2,:,:,:)=acc(2,:,:,:)+(Vi-Vi_old)/time%dtold
               acc(3,:,:,:)=acc(3,:,:,:)+(Wi-Wi_old)/time%dtold
            end if
            ! Zero-out LPT source terms
            srcU=0.0_WP; srcV=0.0_WP; srcW=0.0_WP
            ! Sub-iterate
            call lp%get_cfl(lp_dt,cflc=cfl,cfl=cfl)
            if (cfl.gt.0.0_WP) lp_dt=min(lp_dt*time%cflmax/cfl,lp_dt_max)
            dt_done=0.0_WP
            do while (dt_done.lt.time%dtmid)
               ! Decide the timestep size
               mydt=min(lp_dt,time%dtmid-dt_done)
               ! Inject, collide and advance particles
               !call inject_bed(dt=mydt)
               call lp%collide(dt=mydt)
               call lp%advance(dt=mydt,U=fs%U,V=fs%V,W=fs%W,rho=rho,visc=fs%visc,&
               &               stress_x=resU         ,stress_y=resV         ,stress_z=resW         ,&
               &               acc_x   =acc (1,:,:,:),acc_y   =acc (2,:,:,:),acc_z   =acc (3,:,:,:),&
               &               vort_x  =vort(1,:,:,:),vort_y  =vort(2,:,:,:),vort_z  =vort(3,:,:,:),&
               &               gradVF_x=dVFdx        ,gradVF_y=dVFdy        ,gradVF_z=dVFdz        ,&
               &               srcU=tmp1,srcV=tmp2,srcW=tmp3)
               srcU=srcU+tmp1
               srcV=srcV+tmp2
               srcW=srcW+tmp3
               ! Increment
               dt_done=dt_done+mydt
            end do
            ! Get dilatation from particles: call lp%get_dilatation(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W,VFold=VFold,dil=srcM)
            srcM=0.0_WP
            resU=1.0_WP-0.5_WP*(lp%VF+VFold)
            do k=lp%cfg%kmin_,lp%cfg%kmax_; do j=lp%cfg%jmin_,lp%cfg%jmax_; do i=lp%cfg%imin_,lp%cfg%imax_
               ! Need non-zero fluid volume fraction
               if (resU(i,j,k).le.0.0_WP) cycle
               ! Skip walls
               if (lp%cfg%VF(i,j,k).eq.0.0_WP) cycle
               ! Compute consistent divergence source
               srcM(i,j,k)=1.0_WP/(resU(i,j,k))*((lp%VF(i,j,k)-VFold(i,j,k))/time%dtmid &
               &                                +0.5_WP*(fs%U(i,j,k)*sum(lp%grd_x(:,i,j,k)*resU(i-1:i,j,k))+fs%U(i+1,j,k)*sum(lp%grd_x(:,i+1,j,k)*resU(i:i+1,j,k))) &
               &                                +0.5_WP*(fs%V(i,j,k)*sum(lp%grd_y(:,i,j,k)*resU(i,j-1:j,k))+fs%V(i,j+1,k)*sum(lp%grd_y(:,i,j+1,k)*resU(i,j:j+1,k))) &
               &                                +0.5_WP*(fs%W(i,j,k)*sum(lp%grd_z(:,i,j,k)*resU(i,j,k-1:k))+fs%W(i,j,k+1)*sum(lp%grd_z(:,i,j,k+1)*resU(i,j,k:k+1))) )
            end do; end do; end do
            call fs%cfg%sync(srcM)
            ! Deallocate
            deallocate(tmp1,tmp2,tmp3,VFold,rho,dVFdx,dVFdy,dVFdz,vort,acc)
         end block lpt_step
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

         ! Zero out surface tension in the bed
         where (lp%VF.gt.lp%VFst) vf%curv=0.0_WP
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
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
            
            ! Add momentum source term from lpt - divide by fluid volume fraction!
            add_lpt_src: block
               integer :: i,j,k
               srcU=srcU/(1.0_WP-lp%VF)
               srcV=srcV/(1.0_WP-lp%VF)
               srcW=srcW/(1.0_WP-lp%VF)
               do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
                  if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcU(i-1:i,j,k))
                  if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcV(i,j-1:j,k))
                  if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcW(i,j,k-1:k))
               end do; end do; end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
            end block add_lpt_src
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU!/fs%rho_U
            fs%V=2.0_WP*fs%V-fs%Vold+resV!/fs%rho_V
            fs%W=2.0_WP*fs%W-fs%Wold+resW!/fs%rho_W
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr(src=srcM)
            call fs%get_div(src=srcM)
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
         call fs%get_div(src=srcM)
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(1,i)=0.5_WP*lp%p(i)%d
                  pmesh%vec(:,1,i)=lp%p(i)%vel
               end do
            end block update_pmesh
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call lp%get_max()
         call sfile%write()
         call pfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Bed injection
   subroutine inject_bed(dt)
      implicit none
      real(WP), intent(in) :: dt
      integer , save :: ibed=0
      real(WP), save :: tbed=0.0_WP
      ! Root process selects particles
      if (lp%cfg%amRoot) then
         lp%np_new=0
         inject_loop: do while (ibed+1.le.nbed)
            ! Increment counter
            ibed=ibed+1
            ! Test if particle should be injected
            if (pbed(ibed)%pos(2)+vbed*(tbed+dt).lt.lp%cfg%y(lp%cfg%jmax+1)) then
               ! Add the particle
               lp%np_new=lp%np_new+1
               call lp%resize(lp%np_+lp%np_new)
               lp%p(lp%np_+lp%np_new)=pbed(ibed)
               lp%p(lp%np_+lp%np_new)%pos(2)=pbed(ibed)%pos(2)+vbed*(tbed+dt)
               lp%p(lp%np_+lp%np_new)%ind=lp%cfg%get_ijk_global(lp%p(lp%np_+lp%np_new)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               lp%p(lp%np_+lp%np_new)%vel=[0.0_WP,vbed,0.0_WP]
            else
               ! We went too far, rewind and exit
               ibed=ibed-1
               exit inject_loop
            end if
         end do inject_loop
      end if
      ! Increment injection time
      tbed=tbed+dt
      ! If bed is in free fall, update vbed
      vbed=vbed+lp%gravity(2)*dt
      ! Synchronize
      call lp%sync()
   end subroutine inject_bed
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(srcU,srcV,srcW,srcM,resU,resV,resW,Ui,Vi,Wi,Ui_old,Vi_old,Wi_old)
      
   end subroutine simulation_final
   
   
   !> Function that localizes the top (y+) of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1.and.k.ge.pg%kmin.and.k.le.pg%kmax) isIn=.true.
   end function yp_locator
   
   
end module simulation
