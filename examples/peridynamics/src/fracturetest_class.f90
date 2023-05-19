!> Definition of a solid simulation superclass that drives an lss_class object
module fracturetest_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use lss_class,         only: lss
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: fracturetest

   !> Solid simulation object
   type :: fracturetest
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config
      type(config) :: cfg
      
      !> Get Lagrangian solid solver and timetracker
      type(lss),         public :: ls
      type(timetracker), public :: time
      
      !> Ensight postprocessing
      type(partmesh) :: pmesh
      type(ensight)  :: ens_out
      type(event)    :: ens_evt
      
      !> Monitor file
      type(monitor) :: mfile
      real(WP) :: Xbeam1,Xbeam2
      integer :: ibeam1,ibeam2

   contains
      procedure :: init                    !< Initialize solid simulation
      procedure :: step                    !< Advance solid simulation by one time step
      procedure :: final                   !< Finalize solid simulation
   end type fracturetest
   
   
contains
   

   !> Initialization of solid simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(fracturetest), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_fracturetest')
      
      ! Create a config
      initialize_cfg: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x(nx+1))
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='solid')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block initialize_cfg
      
      ! Initialize time tracker
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
      end block initialize_timetracker
      
      ! Initialize Lagrangian solid solver
      initialize_lss: block
         use sgrid_class, only: sgrid,cartesian
         integer :: i,j,k,nx,ny,nz,np
         real(WP) :: dx
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         
         ! Create solver
         this%ls=lss(cfg=this%cfg,name='solid')
         
         ! Discretize plate
         call this%input%read('Solid dx',dx)
         Lx=0.2_WP; Ly=0.1_WP; Lz=0.025_WP
         
         ! Set material properties
         this%ls%elastic_modulus=190.0e9_WP
         this%ls%poisson_ratio=0.25_WP
         this%ls%rho=7800.0_WP
         this%ls%crit_energy=6.9e4_WP
         this%ls%dV=dx**3
         this%ls%delta=3.0_WP*dx
         
         ! Only root process initializes solid particles
         if (this%ls%cfg%amRoot) then
            ! Create simple rectilinear grid
            nx=int(Lx/dx)
            ny=int(Ly/dx)
            nz=int(Lz/dx)
            allocate(x(1:nx+1),y(1:ny+1),z(1:nz+1))
            do i=1,nx+1
               x(i)=real(i-1,WP)*dx-0.5_WP*Lx
            end do
            do j=1,ny+1
               y(j)=real(j-1,WP)*dx-Ly
            end do
            do k=1,nz+1
               z(k)=real(k-1,WP)*dx-0.5_WP*Lz
            end do
            ! General serial grid object
            grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='elements')
            ! Loop over mesh and create particles
            np=0
            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     ! Remove notches
                     if (grid%ym(j).gt.-0.05_WP.and.abs(grid%xm(i)).gt.0.0225_WP.and.abs(grid%xm(i)).lt.0.0275_WP) cycle
                     ! Increment particle
                     np=np+1
                     call this%ls%resize(np)
                     ! Set position
                     this%ls%p(np)%pos=[grid%xm(i),grid%ym(j),grid%zm(k)]
                     ! Set object id and velocity
                     this%ls%p(np)%id=1
                     this%ls%p(np)%vel=0.0_WP
                     ! Zero out force
                     this%ls%p(np)%Abond=0.0_WP
                     ! Locate the particle on the mesh
                     this%ls%p(np)%ind=this%ls%cfg%get_ijk_global(this%ls%p(np)%pos,[this%ls%cfg%imin,this%ls%cfg%jmin,this%ls%cfg%kmin])
                     ! Assign a unique integer to particle
                     this%ls%p(np)%i=np
                     ! Activate the particle
                     this%ls%p(np)%flag=0

                     ! Prescribe velocity on side
                     !if (this%ls%p(np)%pos(2).gt.-3.0_WP*dx.and.abs(this%ls%p(np)%pos(1)).lt.0.025_WP) then
                     !   this%ls%p(np)%id=-1
                     !   this%ls%p(np)%vel=[0.0_WP,-22.0_WP,0.0_WP]
                     !end if
                     
                     ! Prescribe force on side
                     if (this%ls%p(np)%pos(2).gt.-3.0_WP*dx.and.abs(this%ls%p(np)%pos(1)).lt.0.025_WP) then
                        this%ls%p(np)%id=0
                        this%ls%p(np)%Abond=[0.0_WP,-1.0e6_WP,0.0_WP]
                     end if
                  end do
               end do
            end do
         end if
         
         ! Communicate particles
         call this%ls%sync()
         
         ! Initalize bonds
         call this%ls%bond_init()
         
      end block initialize_lss
      
      ! Create partmesh object for visualizing Lagrangian particles
      create_pmesh: block
         use lss_class, only: max_bond
         integer :: i,n,nbond
         this%pmesh=partmesh(nvar=2,nvec=2,name='solid')
         this%pmesh%varname(1)='failfrac'
         this%pmesh%varname(2)='dilatation'
         this%pmesh%vecname(1)='velocity'
         this%pmesh%vecname(2)='bond_force'
         call this%ls%update_partmesh(this%pmesh)
         do i=1,this%ls%np_
            this%pmesh%var(1,i)=0.0_WP
            nbond=0
            do n=1,max_bond
               if (this%ls%p(i)%ibond(n).gt.0) nbond=nbond+1
            end do
            if (this%ls%p(i)%nbond.gt.0) then
               this%pmesh%var(1,i)=1.0_WP-real(nbond,WP)/real(this%ls%p(i)%nbond,WP)
            else
               this%pmesh%var(1,i)=0.0_WP
            end if
            this%pmesh%var(2,i)  =this%ls%p(i)%dil
            this%pmesh%vec(:,1,i)=this%ls%p(i)%vel
            this%pmesh%vec(:,2,i)=this%ls%p(i)%Abond
         end do
      end block create_pmesh
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='fracturetest')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Only output the particles
         call this%ens_out%add_particle('particles',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%ls%get_cfl(this%time%dt,this%time%cfl)
         call this%ls%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%cfg%amRoot,'fracturetest')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%ls%np,'Particle number')
         call this%mfile%add_column(this%ls%Umin,'Particle Umin')
         call this%mfile%add_column(this%ls%Umax,'Particle Umax')
         call this%mfile%add_column(this%ls%Vmin,'Particle Vmin')
         call this%mfile%add_column(this%ls%Vmax,'Particle Vmax')
         call this%mfile%add_column(this%ls%Wmin,'Particle Wmin')
         call this%mfile%add_column(this%ls%Wmax,'Particle Wmax')
         call this%mfile%write()
      end block create_monitor
      
   end subroutine init
   
   
   !> Take one time step of the solid simulation
	subroutine step(this)
		implicit none
		class(fracturetest), intent(inout) :: this
      
      ! Increment time
      call this%ls%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Advance solid solver
      call this%ls%advance(this%time%dt)
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         update_pmesh: block
            use lss_class, only: max_bond
            integer :: i,n,nbond
            call this%ls%update_partmesh(this%pmesh)
            do i=1,this%ls%np_
               nbond=0
               do n=1,max_bond
                  if (this%ls%p(i)%ibond(n).gt.0) nbond=nbond+1
               end do
               if (this%ls%p(i)%nbond.gt.0) then
                  this%pmesh%var(1,i)=1.0_WP-real(nbond,WP)/real(this%ls%p(i)%nbond,WP)
               else
                  this%pmesh%var(1,i)=0.0_WP
               end if
               this%pmesh%var(2,i)  =this%ls%p(i)%dil
               this%pmesh%vec(:,1,i)=this%ls%p(i)%vel
               this%pmesh%vec(:,2,i)=this%ls%p(i)%Abond
            end do
         end block update_pmesh
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%ls%get_max()
      call this%mfile%write()
      
   end subroutine step
   
   
   !> Finalize solid simulation
   subroutine final(this)
		implicit none
		class(fracturetest), intent(inout) :: this
		
		! Deallocate work arrays
		
	end subroutine final
   
   
end module fracturetest_class