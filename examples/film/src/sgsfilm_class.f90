!> Definition for an sgsfilm class
module sgsfilm_class
   use precision,         only: WP
   use config_class,      only: config
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: sgsfilm
   
   !> SGS film object
   type :: sgsfilm
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(timetracker) :: time  !< Time info
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      
      !> Velocity arrays
      real(WP), dimension(:,:,:), allocatable :: U,V,W      !< Face velocities
      
      !> Problem definition
      real(WP) :: dh
      integer :: nh
      real(WP), dimension(:,:), allocatable :: ph
      
   contains
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type sgsfilm
   

contains
   
   
   !> Initialization of sgsfilm simulation
   subroutine init(this)
      implicit none
      class(sgsfilm), intent(inout) :: this
      
      
      ! Create the sgsfilm mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
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
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='sgsfilm')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%U(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%V(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%W(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Prepare random holes for film perforation
      initialize_holes: block
         use param,    only: param_read
         use random,   only: random_uniform
         use parallel, only: MPI_REAL_WP
         use mpi_f08,  only: MPI_BCAST
         integer  :: n,nn,ierr
         real(WP) :: xh,zh,xxh,zzh
         logical  :: is_overlap
         real(WP), parameter :: safety_margin=1.5_WP
         ! Read in the hole size
         call param_read('Size of holes',this%dh)
         ! Read in the number of holes
         call param_read('Number of holes',this%nh)
         ! Allocate hole position array
         allocate(this%ph(2,this%nh)); this%ph=0.0_WP
         ! If only one hole, leave it in the middle
         if (this%nh.gt.1) then
            ! Root assigns random non-overlapping hole positions
            if (this%cfg%amRoot) then
               n=0
               do while (n.lt.this%nh)
                  ! Draw a random position on the film
                  xh=random_uniform(lo=this%cfg%x(this%cfg%imin),hi=this%cfg%x(this%cfg%imax+1))
                  zh=random_uniform(lo=this%cfg%z(this%cfg%kmin),hi=this%cfg%z(this%cfg%kmax+1))
                  ! Compare to all previous holes
                  is_overlap=.false.
                  do nn=1,n
                     ! Get the position of the other hole
                     xxh=this%ph(1,nn); zzh=this%ph(2,nn)
                     ! Account for periodicity
                     if (xxh-xh.gt.+0.5_WP*this%cfg%xL) xxh=xxh-this%cfg%xL
                     if (xxh-xh.lt.-0.5_WP*this%cfg%xL) xxh=xxh+this%cfg%xL
                     if (zzh-zh.gt.+0.5_WP*this%cfg%zL) zzh=zzh-this%cfg%zL
                     if (zzh-zh.lt.-0.5_WP*this%cfg%zL) zzh=zzh+this%cfg%zL
                     ! Check for overlap
                     if (norm2([xxh-xh,zzh-zh]).lt.safety_margin*this%dh) is_overlap=.true.
                  end do
                  ! If no overlap was found, add the hole to the list
                  if (.not.is_overlap) then
                     n=n+1
                     this%ph(1,n)=xh
                     this%ph(2,n)=zh
                  end if
               end do
            end if
            ! Broadcoast the hole positions
            call MPI_BCAST(this%ph,2*this%nh,MPI_REAL_WP,0,this%cfg%comm,ierr)
         end if
      end block initialize_holes
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,r2p,remap
         use irl_fortran_interface
         integer :: i,j,k,n
         ! Create a VOF solver
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2p,transport_method=remap,name='VOF')
         ! Initialize a flat film
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Setup and intact film of thickness 1
                  if (this%cfg%y(j).lt.-0.5_WP.and.this%cfg%y(j+1).gt.+0.5_WP) then
                     call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                     call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,+1.0_WP,0.0_WP],+0.5_WP)                           
                     call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[0.0_WP,-1.0_WP,0.0_WP],+0.5_WP)
                     ! Loop through the holes and perforate the film
                     do n=1,this%nh
                        if (sqrt((this%vf%cfg%xm(i)-this%ph(1,n))**2+(this%vf%cfg%zm(k)-this%ph(2,n))**2).lt.0.5_WP*this%dh) then
                           call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                           call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],-1.0_WP)
                        end if
                     end do
                  else
                     call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                     call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],-1.0_WP)
                  end if
               end do
            end do
         end do
         call this%vf%sync_interface()
         ! Reset moments from r2p
         call this%vf%reset_volume_moments()
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         !call this%vf%build_interface() !< We have directly provided the interface
         ! Set interface planes at the boundaries
         call this%vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         !call this%vf%reset_volume_moments() !< We have already done that
      end block create_and_initialize_vof
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         this%smesh=surfmesh(nvar=5,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='curv'
         this%smesh%varname(3)='edge_sensor'
         this%smesh%varname(4)='thin_sensor'
         this%smesh%varname(5)='thickness'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Also populate nplane variable
         this%smesh%var(1,:)=1.0_WP
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                        this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                        this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                        this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='sgsfilm')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_surface('vofplic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%cfg%amRoot,'simulation')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%write()
      end block create_monitor
      
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(sgsfilm), intent(inout) :: this
      
      ! Increment time
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF

      ! Set velocity
      this%U=0.0_WP
      this%V=0.0_WP
      this%W=0.0_WP
      
      ! VOF solver step
      call this%vf%advance(dt=this%time%dt,U=this%U,V=this%V,W=this%W)
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surfmesh object
         update_smesh: block
            use irl_fortran_interface
            integer :: i,j,k,nplane,np
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate nplane variable
            this%smesh%var(1,:)=1.0_WP
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                           this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                           this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                           this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Perform ensight output
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%vf%get_max()
      call this%mfile%write()
      
      
   end subroutine step
   
   
   !> Finalize film simulation
   subroutine final(this)
      implicit none
      class(sgsfilm), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%U,this%V,this%W)
      
   end subroutine final
   
   
end module sgsfilm_class