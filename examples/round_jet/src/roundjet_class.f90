!> Definition for a round jet class
module roundjet_class
   use precision,         only: WP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   use stracker_class,    only: stracker
   implicit none
   private
   
   public :: roundjet
   
   !> roundjet object
   type :: roundjet
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf     !< Volume fraction solver
      type(tpns)        :: fs     !< Two-phase flow solver
      type(hypre_str)   :: ps     !< Structured Hypre linear solver for pressure
      type(timetracker) :: time   !< Time info
      
      !> Implicit solver
      logical     :: use_implicit !< Is an implicit solver used?
      type(ddadi) :: vs           !< DDADI solver for velocity
      
      !> SGS modeling
      logical        :: use_sgs   !< Is an LES model used?
      type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh     !< Surface mesh for interface
      type(ensight)  :: ens_out   !< Ensight output for flow variables
      type(event)    :: ens_evt   !< Event trigger for Ensight output
      
      !> Simulation monitoring files
      type(monitor) :: mfile      !< General simulation monitoring
      type(monitor) :: cflfile    !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      
      !> Provide a pardata object for restarts
      logical       :: restarted  !< Is the simulation restarted?
      type(pardata) :: df         !< Pardata object for restart I/O
      type(event)   :: save_evt   !< Event to trigger restart I/O

      !> Droplet diameter postprocessing
      type(event) :: drop_evt     !< Event to trigger droplet diameter output
      
      !> Include structure tracker
      logical        :: use_stracker
      type(stracker) :: strack
      
   contains
      procedure :: init                            !< Initialize round jet simulation
      procedure :: step                            !< Advance round jet simulation by one time step
      procedure :: final                           !< Finalize round jet simulation
      procedure :: analyze_drops                   !< Post-process droplet diameter using stracker
      procedure :: analyze_merge_split             !< Post-process topology change events using stracker
   end type roundjet
   
   ! Structure object
   type :: struct_stats
      real(WP) :: vol
      real(WP), dimension(3) :: pos
      real(WP), dimension(3) :: vel
      real(WP), dimension(3,3) :: Imom
      real(WP), dimension(3) :: lengths
      real(WP), dimension(3,3) :: axes
   end type struct_stats
   
   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=4
   
contains
   
   
   !> Function that defines a level set function for a cylinder
   function levelset_cylinder(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(2)**2+xyz(3)**2)
   end function levelset_cylinder
   
   
   !> Perform droplet analysis
   subroutine analyze_drops(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      use string,    only: str_medium
      use filesys,   only: makedir,isdir
      implicit none
      class(roundjet), intent(inout) :: this
      character(len=str_medium) :: filename,timestamp
      real(WP), dimension(:), allocatable :: dvol
      integer :: iunit,n,m,ierr
      ! Allocate droplet volume array
      allocate(dvol(1:this%strack%nstruct)); dvol=0.0_WP
      ! Loop over individual structures
      do n=1,this%strack%nstruct
         ! Loop over cells in structure and accumulate volume
         do m=1,this%strack%struct(n)%n_
            dvol(n)=dvol(n)+this%cfg%vol(this%strack%struct(n)%map(1,m),this%strack%struct(n)%map(2,m),this%strack%struct(n)%map(3,m))*&
            &                 this%vf%VF(this%strack%struct(n)%map(1,m),this%strack%struct(n)%map(2,m),this%strack%struct(n)%map(3,m))
         end do
      end do
      ! Reduce volume data
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,this%strack%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      ! Only root process outputs to a file
      if (this%cfg%amRoot) then
         if (.not.isdir('diameter')) call makedir('diameter')
         filename='diameter_'; write(timestamp,'(es12.5)') this%time%t
         open(newunit=iunit,file='diameter/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         do n=1,this%strack%nstruct
            ! Output list of diameters
            write(iunit,'(999999(es12.5,x))') (6.0_WP*dvol(n)/Pi)**(1.0_WP/3.0_WP)
         end do
         close(iunit)
      end if
   end subroutine analyze_drops
   
   
   !> Perform merge/split analysis
   subroutine analyze_merge_split(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(roundjet), intent(inout) :: this
      integer :: iunit
      logical :: file_exists
      type(struct_stats) :: mystat
      
      ! Open csv file
      if (this%cfg%amRoot) then
         ! Check if csv file exists, if not create a new file with headers
         inquire(file='merge_split.csv',exist=file_exists)
         if (.not.file_exists) then
            open(newunit=iunit,file='merge_split.csv',form='formatted',status='replace',action='write')
            write(iunit,'(a)') 'EventCount, EventType, OldIDs, NewID, Time, NewVol, X, Y, Z, U, V, W, L1, L2, L3'
            close(iunit)
         end if
         ! Open the csv file
         open(newunit=iunit,file='merge_split.csv',form='formatted',status='old',position='append',action='write')
      end if
      
      ! Analyze merge events
      analyze_merges: block
         integer :: n,nn
         ! Traverse merge events
         do n=1,this%strack%nmerge_master
            call compute_struct_stats(this%strack%merge_master(n)%newid,mystat)
            if (this%cfg%amRoot) then
               ! Write merge data to file
               this%strack%eventcount=this%strack%eventcount+1
               write(iunit,"(I0)",      advance="no")  this%strack%eventcount
               write(iunit,"(A)",       advance="no")  ', Merge,'
               do nn=1,this%strack%merge_master(n)%noldid
                  write(iunit,"(I0)",   advance="no")  this%strack%merge_master(n)%oldids(nn)
                  write(iunit,"(A)",    advance="no")  ';'
               end do
               write(iunit,"(A)",       advance="no")   ','
               write(iunit,"(I0)",      advance="no")  this%strack%merge_master(n)%newid
               write(iunit,"(A)",       advance="no")  ','
               write(iunit,"(ES12.5)",  advance="no")  this%time%t
               write(iunit,"(A)",       advance="no")   ','
               write(iunit,"(ES22.16)", advance="yes") mystat%vol
            end if 
         end do
      end block analyze_merges
      
      ! Analyze split events
      analyze_splits: block
         integer :: n,nn
         ! Traverse split events
         do n=1,this%strack%nsplit_master
            ! Write stats for each new structure after split
            do nn=1,this%strack%split_master(n)%nnewid
               call compute_struct_stats(this%strack%split_master(n)%newids(nn),mystat)
               if (this%cfg%amRoot) then
                  ! Write merge data to file
                  this%strack%eventcount=this%strack%eventcount+1
                  write(iunit,"(I0)",      advance="no")  this%strack%eventcount
                  write(iunit,"(A)",       advance="no")  ', Split,'
                  write(iunit,"(I0)",      advance="no")  this%strack%split_master(n)%oldid
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(I0)",      advance="no")  this%strack%split_master(n)%newids(nn)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES12.5 )", advance="no")  this%time%t
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES22.16)", advance="no")  mystat%vol
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%pos(1)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%pos(2)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%pos(3)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%vel(1)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%vel(2)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%vel(3)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%lengths(1)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  mystat%lengths(2)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="yes") mystat%lengths(3)
               end if
            end do
         end do
      end block analyze_splits
      
      ! Close file
      if (this%cfg%amRoot) close(iunit)
      
   contains
      
      ! Stats calculation for a structure
      subroutine compute_struct_stats(id,stats)
         implicit none 
         integer, intent(in) :: id
         type(struct_stats), intent(inout) :: stats
         
         integer :: n,m,nn
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
         do n=1,this%strack%nstruct
            ! Only deal with structure matching newid
            if (this%strack%struct(n)%id.eq.id) then
               
               ! Periodicity
               per_x = this%strack%struct(n)%per(1)
               per_y = this%strack%struct(n)%per(2)
               per_z = this%strack%struct(n)%per(3)
               
               ! Loop over cells in new structure and accumulate statistics
               do m=1,this%strack%struct(n)%n_

                  ! Indices of cells in structure
                  ii=this%strack%struct(n)%map(1,m) 
                  jj=this%strack%struct(n)%map(2,m) 
                  kk=this%strack%struct(n)%map(3,m)

                  ! Location of struct node
                  xtmp = this%strack%vf%cfg%xm(ii)-per_x*this%strack%vf%cfg%xL
                  ytmp = this%strack%vf%cfg%ym(jj)-per_y*this%strack%vf%cfg%yL
                  ztmp = this%strack%vf%cfg%zm(kk)-per_z*this%strack%vf%cfg%zL

                  ! Volume
                  vol_struct = vol_struct + this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  
                  ! Center of gravity
                  x_vol = x_vol + xtmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  y_vol = y_vol + ytmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  z_vol = z_vol + ztmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  
                  ! Average velocity inside struct
                  u_vol = u_vol + this%Ui(ii,jj,kk)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  v_vol = v_vol + this%Vi(ii,jj,kk)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  w_vol = w_vol + this%Wi(ii,jj,kk)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
               end do
            end if
         end do

         ! Sum parallel stats
         call MPI_ALLREDUCE(MPI_IN_PLACE,vol_struct,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,x_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,y_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,z_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,u_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,v_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,w_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         
         if (vol_struct.gt.0.0_WP) then 
            ! Moments of inertia
            Imom=0.0_WP
            do n=1,this%strack%nstruct
               do nn=1,this%strack%nmerge_master
                  ! Only deal with structure matching newid
                  if (this%strack%struct(n)%id.eq.this%strack%merge_master(nn)%newid) then
                     
                     ! Periodicity
                     per_x = this%strack%struct(n)%per(1)
                     per_y = this%strack%struct(n)%per(2)
                     per_z = this%strack%struct(n)%per(3)
                     
                     ! Loop over cells in new structure and accumulate statistics
                     do m=1,this%strack%struct(n)%n_
                        
                        ! Indices of cells in structure
                        ii=this%strack%struct(n)%map(1,m) 
                        jj=this%strack%struct(n)%map(2,m) 
                        kk=this%strack%struct(n)%map(3,m)

                        ! Location of struct node
                        xtmp = this%strack%vf%cfg%xm(ii)-per_x*this%strack%vf%cfg%xL-x_vol/vol_struct
                        ytmp = this%strack%vf%cfg%ym(jj)-per_y*this%strack%vf%cfg%yL-y_vol/vol_struct
                        ztmp = this%strack%vf%cfg%zm(kk)-per_z*this%strack%vf%cfg%zL-z_vol/vol_struct

                        ! Moment of Inertia
                        Imom(1,1) = Imom(1,1) + (ytmp**2 + ztmp**2)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                        Imom(2,2) = Imom(2,2) + (xtmp**2 + ztmp**2)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                        Imom(3,3) = Imom(3,3) + (xtmp**2 + ytmp**2)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                        
                        Imom(1,2) = Imom(1,2) - xtmp*ytmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                        Imom(1,3) = Imom(1,3) - xtmp*ztmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                        Imom(2,3) = Imom(2,3) - ytmp*ztmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                     end do 
                  end if
               end do 
            end do
            
            ! Sum parallel stats on Imom
            do n=1,3
               call MPI_ALLREDUCE(MPI_IN_PLACE,Imom(:,n),3,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
            end do
            
            ! Characteristic lengths and principle axes
            ! Eigenvalues/eigenvectors of moments of inertia tensor
            A = Imom
            n = 3
            call dsyev('V','U',n,Imom,n,d,work,lwork,info)
            ! Get rid of very small negative values (due to machine accuracy)
            d = max(0.0_WP,d)
            ! Store characteristic lengths
            lengths(1) = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol_struct) ! Need a check to prevent dividing by 0 
            lengths(2) = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol_struct)
            lengths(3) = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol_struct)
            ! Zero out length in 3rd dimension if 2D
            if (this%strack%vf%cfg%nx.eq.1.or.this%strack%vf%cfg%ny.eq.1.or.this%strack%vf%cfg%nz.eq.1) lengths(3)=0.0_WP
            ! Store principal axes
            axes(:,:) = A
            
            ! Finish computing qantities
            stats%vol     = vol_struct
            stats%pos(1)  = x_vol/vol_struct
            stats%pos(2)  = y_vol/vol_struct
            stats%pos(3)  = z_vol/vol_struct
            stats%vel(1)  = u_vol/vol_struct
            stats%vel(2)  = v_vol/vol_struct
            stats%vel(3)  = w_vol/vol_struct
            stats%Imom    = Imom 
            stats%lengths = lengths
            stats%axes    = axes
         end if 
         
      end subroutine compute_struct_stats
      
   end subroutine analyze_merge_split
   
   
   !> Initialization of roundjet simulation
   subroutine init(this)
      use param, only: param_read
      implicit none
      class(roundjet), intent(inout) :: this
      
      
      ! Initialize the config
      initialize_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,xshift,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         
         ! Read in grid definition
         call param_read('[Jet] Lx',Lx); call param_read('[Jet] nx',nx); allocate(x_uni(nx+1))
         call param_read('[Jet] Ly',Ly); call param_read('[Jet] ny',ny); allocate(y_uni(ny+1))
         call param_read('[Jet] Lz',Lz); call param_read('[Jet] nz',nz); allocate(z_uni(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! Add stretching
         call param_read('[Jet] Stretched cells in yz',ns_yz,default=0)
         if (ns_yz.gt.0) call param_read('[Jet] Stretch ratio in yz',sratio_yz)
         call param_read('[Jet] Stretched cells in x' ,ns_x ,default=0)
         if (ns_x .gt.0) call param_read('[Jet] Stretch ratio in x' ,sratio_x )
         allocate(x(nx+1+1*ns_x )); x(      1:      1+nx)=x_uni
         allocate(y(ny+1+2*ns_yz)); y(ns_yz+1:ns_yz+1+ny)=y_uni
         allocate(z(nz+1+2*ns_yz)); z(ns_yz+1:ns_yz+1+nz)=z_uni
         do i=nx+2,nx+1+ns_x
            x(i)=x(i-1)+sratio_x*(x(i-1)-x(i-2))
         end do
         do j=ns_yz,1,-1
            y(j)=y(j+1)+sratio_yz*(y(j+1)-y(j+2))
         end do
         do j=ns_yz+2+ny,ny+1+2*ns_yz
            y(j)=y(j-1)+sratio_yz*(y(j-1)-y(j-2))
         end do
         do k=ns_yz,1,-1
            z(k)=z(k+1)+sratio_yz*(z(k+1)-z(k+2))
         end do
         do k=ns_yz+2+nz,nz+1+2*ns_yz
            z(k)=z(k-1)+sratio_yz*(z(k-1)-z(k-2))
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='jet')
         
         ! Read in partition
         call param_read('[Jet] Partition',partition)
         
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         
         ! No walls in the atomization domain
         this%cfg%VF=1.0_WP

      end block initialize_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('[Jet] Max timestep size',this%time%dtmax)
         call param_read('[Jet] Max cfl number',this%time%cflmax)
         call param_read('[Jet] Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,plicnet,remap,remap_storage
         use mms_geom,  only: cube_refine_vol
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Do we use stracker
         call param_read('[Jet] Use stracker',this%use_stracker)
         ! Create a VOF solver with plicnet reconstruction
         if (this%use_stracker) then
            call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=remap_storage,name='VOF')
            call this%strack%initialize(vf=this%vf,phase=0,make_label=label_liquid,name='round_jet')
         else
            call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=remap,name='VOF')
         end if
         ! Initialize to cylindrical interface
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cylinder,0.0_WP,amr_ref_lvl)
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  else
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
                  ! Clip cylinder after two cells
                  if (i.ge.this%cfg%imin+2) then
                     this%vf%VF(i,j,k)=0.0_WP
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set simple full-liquid/full-gas interface planes in geometric overlap cells
         call this%vf%set_full_bcond()
         ! Now apply Neumann condition on interface at inlet to have proper round injection
         neumann_irl: block
            use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
            &                                setNumberOfPlanes,setPlane,matchVolumeFraction
            real(WP), dimension(1:4) :: plane
            type(RectCub_type) :: cell
            call new(cell)
            if (this%vf%cfg%iproc.eq.1) then
               do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
                  do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                     do i=this%vf%cfg%imino,this%vf%cfg%imin-1
                        ! Extract plane data and copy in overlap
                        plane=getPlane(this%vf%liquid_gas_interface(this%vf%cfg%imin,j,k),0)
                        call construct_2pt(cell,[this%vf%cfg%x(i  ),this%vf%cfg%y(j  ),this%vf%cfg%z(k  )],&
                        &                       [this%vf%cfg%x(i+1),this%vf%cfg%y(j+1),this%vf%cfg%z(k+1)])
                        plane(4)=dot_product(plane(1:3),[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)])
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        call matchVolumeFraction(cell,this%vf%VF(i,j,k),this%vf%liquid_gas_interface(i,j,k))
                     end do
                  end do
               end do
            end if
         end block neumann_irl
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
         ! Handle stracker initialization
         if (this%use_stracker) then
            init_stracker: block
               use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_INTEGER,MPI_MAX
               integer :: ierr
               where (this%vf%VF.gt.0.0_WP) this%strack%id=1
               this%strack%idcount=maxval(this%strack%id)
               call MPI_ALLREDUCE(MPI_IN_PLACE,this%strack%idcount,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
            end block init_stracker
         end if
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call param_read('Density ratio',this%fs%rho_l)
         call param_read('Reynolds number',this%fs%visc_l); this%fs%visc_l=this%fs%rho_l/this%fs%visc_l
         call param_read('Viscosity ratio',this%fs%visc_g); this%fs%visc_g=this%fs%visc_l/this%fs%visc_g
         call param_read('Weber number',this%fs%sigma); this%fs%sigma=this%fs%rho_l/this%fs%sigma
         ! Inflow on the left
         call this%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp'  ,type=slip           ,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym'  ,type=slip           ,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp'  ,type=slip           ,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm'  ,type=slip           ,face='z',dir=-1,canCorrect=.true. ,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call param_read('[Jet] Pressure iteration',this%ps%maxit)
         call param_read('[Jet] Pressure tolerance',this%ps%rcvg)
         ! Check if we want to use an implicit solver
         call param_read('[Jet] Use implicit solver',this%use_implicit)
         if (this%use_implicit) then
            ! Configure implicit velocity solver
            this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
            ! Setup the solver
            call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         else
            ! Setup the solver
            call this%fs%setup(pressure_solver=this%ps)
         end if
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply convective velocity
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).le.0.5_WP) this%fs%U(i,j,k)=1.0_WP
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Compute MFR through all boundary conditions
         call this%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('[Jet] Use SGS model',this%use_sgs)
         if (this%use_sgs) this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         character(len=str_medium) :: filename
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         integer :: i,j,k
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Jet restart output')
         call param_read('[Jet] Restart output period',this%save_evt%tper)
         ! Read in the I/O partition
         call param_read('[Jet] I/O partition',iopartition)
         ! Check if we are restarting
         call param_read('[Jet] Restart from',filename,default='')
         this%restarted=.false.; if (len_trim(filename).gt.0) this%restarted=.true.
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata=trim(filename))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     ! Check if the second plane is meaningful
                     if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                     else
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                     end if
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Update the band
            call this%vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Compute MFR through all boundary conditions
            call this%fs%get_mfr()
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block handle_restart
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         integer :: i,j,k,nplane,np
         if (this%use_stracker) then
            this%smesh=surfmesh(nvar=1,name='plic')
            this%smesh%varname(1)='id'
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate id variable
            this%smesh%var(1,:)=0.0_WP
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(this%strack%id(i,j,k),WP)
                        end if
                     end do
                  end do
               end do
            end do
         else
            this%smesh=surfmesh(nvar=0,name='plic')
            call this%vf%update_surfmesh(this%smesh)
         end if
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='jet')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('[Jet] Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         if (this%use_stracker) call this%ens_out%add_scalar('id',this%strack%id)
         call this%ens_out%add_surface('plic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'jet')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
      
      ! Create an event for drop size analysis
      if (this%use_stracker) then
         drop_analysis: block
            this%drop_evt=event(time=this%time,name='Drop analysis')
            call param_read('[Jet] Drop analysis period',this%drop_evt%tper)
            if (this%drop_evt%occurs()) call this%analyze_drops()
         end block drop_analysis
      end if
      
      
   contains
      
      !> Function that identifies liquid cells
      logical function label_liquid(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.0_WP) then
            label_liquid=.true.
         else
            label_liquid=.false.
         end if
      end function label_liquid
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(roundjet), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
      
      ! VOF solver step
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Advance stracker
      if (this%use_stracker) then
         call this%strack%advance(make_label=label_liquid)
         call this%analyze_merge_split()
      end if

      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf)
      
      ! Turbulence modeling
      if (this%use_sgs) then
         sgs_modeling: block
            use sgsmodel_class, only: vreman
            integer :: i,j,k
            this%resU=this%fs%rho_g
            call this%fs%get_gradu(this%gradU)
            call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
            do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
                     this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
                     this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
                     this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
                     this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
                  end do
               end do
            end do
         end block sgs_modeling
      end if
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Use implicit solver or not
         if (this%use_implicit) then
            ! Form implicit residuals
            call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
            ! Apply these residuals
            this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
            this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
            this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         else
            ! Apply these residuals
            this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho_U
            this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho_V
            this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho_W
         end if
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         real(WP) :: my_vof_removed
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         if (this%use_stracker) then
            update_smesh: block
               use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call this%vf%update_surfmesh(this%smesh)
               ! Also populate id variable
               this%smesh%var(1,:)=0.0_WP
               np=0
               do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
                  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                     do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; this%smesh%var(1,np)=real(this%strack%id(i,j,k),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
         else
            call this%vf%update_surfmesh(this%smesh)
         end if
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Analyse droplets
      if (this%use_stracker) then
         if (this%drop_evt%occurs()) call this%analyze_drops()
      end if

      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! First plane
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     ! Second plane
                     plane=0.0_WP
                     if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                     P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/jet_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
         end block save_restart
      end if
      
   contains
      
      !> Function that identifies liquid cells
      logical function label_liquid(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.0_WP) then
            label_liquid=.true.
         else
            label_liquid=.false.
         end if
      end function label_liquid
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(roundjet), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      
   end subroutine final
   
   
   !> Function that localizes region of VOF removal
   function vof_removal_layer_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.ge.pg%imax-nlayer.or.&
      &   j.le.pg%jmin+nlayer.or.&
      &   j.ge.pg%jmax-nlayer.or.&
      &   k.le.pg%kmin+nlayer.or.&
      &   k.ge.pg%kmax-nlayer) isIn=.true.
   end function vof_removal_layer_locator
   
   
   !> Function that localizes the right (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

   !> Function that localizes the left (x-) of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator
   
   
   !> Function that localizes the top (y+) of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator
   
   
   !> Function that localizes the bottom (y-) of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator
   
   
   !> Function that localizes the front (z+) of the domain
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator
   
   
   !> Function that localizes the back (z-) of the domain
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator
   
   
end module roundjet_class