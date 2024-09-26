!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,nx,Lx
   use fft3d_class,       only: fft3d
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use stracker_class,    only: stracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
   type(fft3d),       public :: ps
   type(tpns),        public :: fs
   type(timetracker), public :: time
   type(vfs),         public :: vf

   !> Include structure tracker
   type(stracker) :: strack
 
   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(event)    :: ens_evt,drop_evt
   type(surfmesh) :: smesh
  
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,hitfile,cvgfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   
   !> Fluid, forcing, and particle parameters
   real(WP) :: rho,meanU,meanV,meanW
   real(WP) :: Urms0
   real(WP), dimension(3) :: center
   real(WP) :: radius
   
   !> Provide a pardata objects for restarts
   type(pardata) :: df

   !> Type for structure stats
   type :: struct_stats
      real(WP) :: vol
      real(WP) :: x_cg,y_cg,z_cg
      real(WP) :: u_avg,v_avg,w_avg
      real(WP), dimension(3,3) :: Imom
      real(WP), dimension(3) :: lengths
      real(WP), dimension(3,3) :: axes
   end type struct_stats

contains
   
   !> Function that identifies liquid cells
   logical function label_liquid(i,j,k)
      implicit none
      integer, intent(in) :: i,j,k
      if (vf%VF(i,j,k).gt.0.0_WP) then
         label_liquid=.true.
      else
         label_liquid=.false.
      end if
   end function label_liquid

   
   !> Function that defines a level set function for a sphere
   function levelset_sphere(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=radius-sqrt(sum((xyz-center)**2))
   end function levelset_sphere
   
   !> Perform droplet analysis
   subroutine analyse_drops()
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      use string,    only: str_medium
      use filesys,   only: makedir,isdir
      character(len=str_medium) :: filename,timestamp
      real(WP), dimension(:), allocatable :: dvol
      integer :: iunit,n,m,ierr
      ! Allocate droplet volume array
      allocate(dvol(1:strack%nstruct)); dvol=0.0_WP
      ! Loop over individual structures
      do n=1,strack%nstruct
         ! Loop over cells in structure and accumulate volume
         do m=1,strack%struct(n)%n_
            dvol(n)=dvol(n)+cfg%vol(strack%struct(n)%map(1,m),strack%struct(n)%map(2,m),strack%struct(n)%map(3,m))*&
            &                 vf%VF(strack%struct(n)%map(1,m),strack%struct(n)%map(2,m),strack%struct(n)%map(3,m))
         end do
      end do
      ! Reduce volume data
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,strack%nstruct,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
      ! Only root process outputs to a file
      if (cfg%amRoot) then
         if (.not.isdir('diameter')) call makedir('diameter')
         filename='diameter_'; write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file='diameter/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         do n=1,strack%nstruct
            ! Output list of diameters
            write(iunit,'(999999(es12.5,x))') (6.0_WP*dvol(n)/Pi)**(1.0_WP/3.0_WP)
         end do
         close(iunit)
      end if
   end subroutine analyse_drops


   !> Perform merge/split analysis
   subroutine analyze_merge_split
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: iunit
      logical :: file_exists
      type(struct_stats) :: stats

      ! Open the file - Created in simulation_init
      if (cfg%amRoot) open(iunit,file="merge_split.csv",form="formatted",status="old",position="append",action="write")
   
      analyze_merges: block
         integer :: n,nn

         ! Traverse merge events
         do n=1,strack%nmerge_master

            call compute_struct_stats(strack%merge_master(n)%newid,stats)
            if (cfg%amRoot) then 
               ! Write merge data to file
               strack%eventcount = strack%eventcount+1
               write(iunit,"(I0)",      advance="no")  strack%eventcount
               write(iunit,"(A)",       advance="no")  ', Merge,'
               do nn=1,strack%merge_master(n)%noldid
                  write(iunit,"(I0)",   advance="no")  strack%merge_master(n)%oldids(nn)
                  write(iunit,"(A)",    advance="no")  ';'
               end do
               write(iunit,"(A)",       advance="no")   ','
               write(iunit,"(I0)",      advance="no")  strack%merge_master(n)%newid
               write(iunit,"(A)",       advance="no")  ','
               write(iunit,"(ES12.5 )", advance="no")  time%t
               write(iunit,"(A)",       advance="no")   ','
               write(iunit,"(ES22.16)", advance="yes") stats%vol
            end if 
         end do
      end block analyze_merges

      analyze_splits: block 
      integer :: n,nn

         ! Traverse split events
         do n=1,strack%nsplit_master
            ! Write stats for each new structure after split
            do nn=1,strack%split_master(n)%nnewid
               call compute_struct_stats(strack%split_master(n)%newids(nn),stats)
               if (cfg%amRoot) then 

                  ! Write merge data to file
                  strack%eventcount = strack%eventcount+1
                  write(iunit,"(I0)",      advance="no")  strack%eventcount
                  write(iunit,"(A)",       advance="no")  ', Split,'
                  write(iunit,"(I0)",      advance="no")  strack%split_master(n)%oldid
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(I0)",      advance="no")  strack%split_master(n)%newids(nn)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES12.5 )", advance="no")  time%t
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES22.16)", advance="no")  stats%vol
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%x_cg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%y_cg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%z_cg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%u_avg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%v_avg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%w_avg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%lengths(1)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%lengths(2)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="yes")  stats%lengths(3)
               end if 
            end do
         end do
      
      end block analyze_splits

      if (cfg%amRoot) close(iunit)
   
   contains 
      subroutine compute_struct_stats(id,stats)
         implicit none 
         integer, intent(in) :: id
         type(struct_stats), intent(inout) :: stats

         integer :: n,m
         integer :: lwork,info,ierr
         integer :: ii,jj,kk
         integer  :: per_x,per_y,per_z
         real(WP) :: vol_struct
         real(WP) :: x_vol,y_vol,z_vol
         real(WP) :: u_vol,v_vol,w_vol
         real(WP), dimension(3,3) :: Imom
         real(WP) :: xtmp,ytmp,ztmp
         real(WP), dimension(3) :: lengths
         real(WP), dimension(3,3) :: axes
         
         ! Eigenvalues/eigenvectors
         real(WP), dimension(3,3) :: A
         real(WP), dimension(3) :: d
         integer , parameter :: order = 3
         real(WP), dimension(:), allocatable :: work
         real(WP), dimension(1)   :: lwork_query

         
         ! Query optimal work array size
         call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))

         ! Initialize values
         vol_struct    = 0.0_WP ! Structure volume
         x_vol = 0.0_WP; y_vol = 0.0_WP; z_vol = 0.0_WP ! Center of gravity
         u_vol = 0.0_WP; v_vol = 0.0_WP; w_vol = 0.0_WP ! Average velocity inside struct

         ! Find new structure with matching newid
         do n=1,strack%nstruct
            ! Only deal with structure matching newid
            if (strack%struct(n)%id.eq.id) then
               
               ! Periodicity
               per_x = strack%struct(n)%per(1)
               per_y = strack%struct(n)%per(2)
               per_z = strack%struct(n)%per(3)
               
               ! Loop over cells in new structure and accumulate statistics
               do m=1,strack%struct(n)%n_

                  ! Indices of cells in structure
                  ii=strack%struct(n)%map(1,m) 
                  jj=strack%struct(n)%map(2,m) 
                  kk=strack%struct(n)%map(3,m)

                  ! Location of struct node
                  xtmp = strack%vf%cfg%xm(ii)-per_x*strack%vf%cfg%xL
                  ytmp = strack%vf%cfg%ym(jj)-per_y*strack%vf%cfg%yL
                  ztmp = strack%vf%cfg%zm(kk)-per_z*strack%vf%cfg%zL

                  ! Volume
                  vol_struct = vol_struct + strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  
                  ! Center of gravity
                  x_vol = x_vol + xtmp*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  y_vol = y_vol + ytmp*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  z_vol = z_vol + ztmp*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  
                  ! Average velocity inside struct
                  u_vol = u_vol + fs%U(ii,jj,kk)*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  v_vol = v_vol + fs%V(ii,jj,kk)*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  w_vol = w_vol + fs%W(ii,jj,kk)*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
               end do
            end if
         end do

         ! Sum parallel stats
         call MPI_ALLREDUCE(MPI_IN_PLACE,vol_struct,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,x_vol,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,y_vol,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,z_vol,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,u_vol,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,v_vol,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,w_vol,1,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         
         ! Moments of inertia
         Imom=0.0_WP
         do n=1,strack%nstruct
            ! Only deal with structure matching newid
            if (strack%struct(n)%id.eq.id) then

               ! Periodicity
               per_x = strack%struct(n)%per(1)
               per_y = strack%struct(n)%per(2)
               per_z = strack%struct(n)%per(3)
                  
               ! Loop over cells in new structure and accumulate statistics
               do m=1,strack%struct(n)%n_

                  ! Indices of cells in structure
                  ii=strack%struct(n)%map(1,m) 
                  jj=strack%struct(n)%map(2,m) 
                  kk=strack%struct(n)%map(3,m)

                  ! Location of struct node
                  xtmp = strack%vf%cfg%xm(ii)-per_x*strack%vf%cfg%xL-x_vol/vol_struct
                  ytmp = strack%vf%cfg%ym(jj)-per_y*strack%vf%cfg%yL-y_vol/vol_struct
                  ztmp = strack%vf%cfg%zm(kk)-per_z*strack%vf%cfg%zL-z_vol/vol_struct

                  ! Moment of Inertia
                  Imom(1,1) = Imom(1,1) + (ytmp**2 + ztmp**2)*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  Imom(2,2) = Imom(2,2) + (xtmp**2 + ztmp**2)*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  Imom(3,3) = Imom(3,3) + (xtmp**2 + ytmp**2)*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  
                  Imom(1,2) = Imom(1,2) - xtmp*ytmp*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  Imom(1,3) = Imom(1,3) - xtmp*ztmp*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
                  Imom(2,3) = Imom(2,3) - ytmp*ztmp*strack%vf%cfg%vol(ii,jj,kk)*strack%vf%VF(ii,jj,kk)
               end do 
            end if
         end do

         ! Sum parallel stats on Imom
         do n=1,3
            call MPI_ALLREDUCE(MPI_IN_PLACE,Imom(:,n),3,MPI_REAL_WP,MPI_SUM,strack%vf%cfg%comm,ierr)
         end do

         ! Characteristic lengths and principle axes
         ! Eigenvalues/eigenvectors of moments of inertia tensor
         A = Imom
         n = 3
         call dsyev('V','U',n,Imom,n,d,work,lwork,info)
         ! Get rid of very small negative values (due to machine accuracy)
         d = max(0.0_WP,d)
         ! Store characteristic lengths
         lengths(1) = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol_struct)
         lengths(2) = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol_struct)
         lengths(3) = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol_struct)
         ! Zero out length in 3rd dimension if 2D
         if (strack%vf%cfg%nx.eq.1.or.strack%vf%cfg%ny.eq.1.or.strack%vf%cfg%nz.eq.1) lengths(3)=0.0_WP
         ! Store principal axes
         axes(:,:) = A

         ! Finish computing qantities
         stats%vol     = vol_struct
         stats%x_cg    = x_vol/vol_struct
         stats%y_cg    = y_vol/vol_struct
         stats%z_cg    = z_vol/vol_struct
         stats%u_avg   = u_vol/vol_struct
         stats%v_avg   = v_vol/vol_struct
         stats%w_avg   = w_vol/vol_struct
         stats%Imom    = Imom 
         stats%lengths = lengths
         stats%axes    = axes

      end subroutine compute_struct_stats

   end subroutine analyze_merge_split
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU         (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV         (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW         (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui           (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi           (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi           (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax,default=time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only:lvira,plicnet,remap_storage,VFhi,VFlo
         use random, only: random_uniform
         use mms_geom, only: cube_refine_vol
         use precision, only: I4
         use MPI, only: MPI_INTEGER,MPI_MAX
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         integer :: i,j,k,n,si,sj,sk
         integer :: nD,nDrop
         integer(kind=I4), allocatable, dimension(:) :: seed
         integer(kind=I4) :: nseed
         integer :: ierr
         ! Create a VOF solver with r2p reconstruction
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=remap_storage,name='VOF')
         ! Create structure tracker
         call strack%initialize(vf=vf,phase=0,make_label=label_liquid,name='stracker_test')
         ! Initialize our bubble via r2p planes
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
         call param_read('Number of droplet',nDrop);
         ! Provide seed for random number generator
         call random_seed(size=nseed)
         allocate(seed(nseed))
         seed(:)=1
         call random_seed(put=seed)
         do nD=1,nDrop
            center=[random_uniform(vf%cfg%x(vf%cfg%imin),vf%cfg%x(vf%cfg%imax+1)), &
                    random_uniform(vf%cfg%y(vf%cfg%jmin),vf%cfg%y(vf%cfg%jmax+1)), &
                    random_uniform(vf%cfg%z(vf%cfg%kmin),vf%cfg%z(vf%cfg%kmax+1))  ]
                     
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
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
                     vf%VF(i,j,k)=min(1.0_WP,vf%VF(i,j,k)+vol/vf%cfg%vol(i,j,k))
                     if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                        vf%Lbary(:,i,j,k)=v_cent
                        vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                     else
                        vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                        vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     end if
                     ! Set stracker id
                     if (vol.gt.0.0_WP) then
                        strack%id(i,j,k)=nD
                        if (nD.eq.4) then    ! Test growth (merge2)
                           strack%id(i,j,k)=0
                        end if   
                     end if
                  end do
               end do
            end do
         end do
         call vf%cfg%sync(vf%VF)
         call vf%cfg%sync(strack%id)

         ! Initialize id counter to be consistent with id's
         strack%idcount=maxval(strack%id)
         call MPI_ALLREDUCE(maxval(strack%id),strack%idcount,1,MPI_INTEGER,MPI_MAX,cfg%comm,ierr)
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
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
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
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      ! Prepare initial velocity field
      initialize_velocity: block
         use random,    only: random_normal
         use param,     only: param_exists
         use mathtools, only: Pi
         use string,    only: str_long,str_medium
         use messager,  only: die,log
         integer :: i,j,k
         ! Gaussian initial field
         Urms0=1.0
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
         fs%U=fs%U-meanU
         fs%V=fs%V-meanV
         fs%W=fs%W-meanW
         ! Project to ensure divergence-free
         call fs%get_div()
         rho=1.0_WP
         fs%psolv%rhs=-fs%cfg%vol*fs%div*rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*resU
         fs%V=fs%V-time%dt*resV
         fs%W=fs%W-time%dt*resW
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block initialize_velocity
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for structure id
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='id'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate id variable
         smesh%var(1,:)=0.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(strack%id(i,j,k),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='HIT')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('id',strack%id)
         call ens_out%add_surface('vofplic',smesh)
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
         call mfile%add_column(vf%flotsam_error,'Flotsam error')
         call mfile%add_column(vf%thinstruct_error,'Film error')
         call mfile%add_column(vf%SDint,'SD integral')
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

      create_merge_split: block 
         integer :: iunit
         logical :: file_exists
         if (cfg%amRoot) then
            ! Check if the file exists
            INQUIRE(FILE="merge_split.csv", EXIST=file_exists)
            if (.not.file_exists) then
               ! Create a new file with headers if it doesn't exist
               open(newunit=iunit, file="merge_split.csv", form="formatted", status="replace", action="write")
               write(iunit, "(A)") "Event Count, Event Type, Old IDs, New ID, Time, New Vol, X, Y, Z, U, V, W, L1, L2, L3"
               close(iunit)
            end if
         end if
      end block create_merge_split

      ! Initialize an event for drop size analysis
      drop_analysis: block
         drop_evt=event(time=time,name='Drop analysis')
         call param_read('Drop analysis period',drop_evt%tper)
         if (drop_evt%occurs()) call analyse_drops()
      end block drop_analysis
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
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

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

         ! Advance stracker
         call strack%advance(make_label=label_liquid)
         call analyze_merge_split

         ! ! Output master stracker lists
         ! print_merge_master: block
         ! integer :: n,nn,ierr
         ! do nn=0,cfg%nproc
         !    if (cfg%rank.eq.nn) then
         !       print*,'rank=',nn
         !       if (strack%nmerge_master.gt.0) print*,'Master merge list'
         !       do n=1,strack%nmerge_master
         !          print*,'noldid=',strack%merge_master(n)%noldid
         !          print*,' oldid=',strack%merge_master(n)%oldids(1:strack%merge_master(n)%noldid)
         !          print*,' newid=',strack%merge_master(n)%newid
         !       end do
         !    end if
         !    call MPI_BARRIER(cfg%comm,ierr)
         ! end do   
         ! end block print_merge_master

         ! print_split_master: block
         !    integer :: n,nn,ierr
         !    do nn=0,cfg%nproc
         !       if (cfg%rank.eq.nn) then
         !          print*,'rank=',nn   
         !          if (strack%nsplit_master.gt.0) print*,'Master split list'
         !          do n=1,strack%nsplit_master
         !             print*,' oldid=',strack%split_master(n)%oldid
         !             print*,'nnewid=',strack%split_master(n)%nnewid
         !             print*,' newid=',strack%split_master(n)%newids(1:strack%split_master(n)%nnewid)         
         !          end do
         !       end if
         !       call MPI_BARRIER(cfg%comm,ierr)
         !    end do
         ! end block print_split_master
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU!/fs%rho_U
            fs%V=2.0_WP*fs%V-fs%Vold+resV!/fs%rho_V
            fs%W=2.0_WP*fs%W-fs%Wold+resW!/fs%rho_W
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
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
            fs%U=fs%U-time%dt*resU!/fs%rho_U
            fs%V=fs%V-time%dt*resV!/fs%rho_V
            fs%W=fs%W-time%dt*resW!/fs%rho_W
            
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
               use irl_fortran_interface
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Also populate id variable
               smesh%var(1,:)=0.0_WP
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh%var(1,np)=real(strack%id(i,j,k),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            call ens_out%write_data(time%t)
         end if
         
         ! Analyse droplets
         if (drop_evt%occurs()) call analyse_drops()
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         call hitfile%write()
         call cvgfile%write()
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)
      
   end subroutine simulation_final
   
end module simulation
