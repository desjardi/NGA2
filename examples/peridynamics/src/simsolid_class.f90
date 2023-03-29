!> Definition of a solid simulation superclass that drives an lss_class object
module simsolid_class
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
   
   public :: simsolid

   !> Solid simulation object
   type :: simsolid
      
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
      
   contains
      procedure :: init                    !< Initialize solid simulation
      procedure :: step                    !< Advance solid simulation by one time step
      procedure :: final                   !< Finalize solid simulation
   end type simsolid
   
   
contains
   
   
   !> Initialization of solid simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(simsolid), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_solid')
      
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
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
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
         use mathtools,   only: twoPi,Pi
         use random,      only: random_uniform
         integer :: i,j,k,np,nd,nh
         real(WP) :: radius,diam,height
         real(WP), dimension(:), allocatable :: x,y
         type(sgrid) :: grid

         ! Create solver
         this%ls=lss(cfg=this%cfg,name='solid')
         
         ! Discretize cylinder
         call this%input%read('Elements per diameter',nd)
         call this%input%read('Cylinder diameter',diam)
         call this%input%read('Cylinder height',height)

         ! Set material properties
         this%ls%elastic_modulus=1.0e4_WP
         this%ls%poisson_ratio=0.333_WP
         this%ls%rho=1000.0_WP
         this%ls%crit_energy=1.0e3_WP
         this%ls%dV=(diam/real(nd,WP))**3

         ! Only root process initializes solid particles
         if (this%ls%cfg%amRoot) then
            ! Create simple rectilinear grid
            nh=ceiling(height/(diam/real(nd,WP)))
            allocate(x(1:nd+1),y(1:nh+1))
            do i=1,nd+1
               x(i)=real(i-1,WP)/real(nd,WP)*diam-0.5_WP*diam
            end do
            do j=1,nh+1
               y(j)=real(j-1,WP)/real(nh,WP)*height-0.5_WP*height
            end do
            ! General serial grid object
            grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=x,xper=.false.,yper=.false.,zper=.false.,name='elements')
            ! Loop over mesh and create particles
            np=0
            do k=1,nd
               do j=1,nh
                  do i=1,nd
                     ! Check if inside sphere
                     radius=sqrt(grid%xm(i)**2+grid%zm(k)**2)
                     if (radius.ge.0.5_WP*diam.or.abs(grid%ym(j)).ge.0.5_WP*height) cycle
                     ! Increment particle
                     np=np+1
                     call this%ls%resize(np)
                     ! Set object id
                     this%ls%p(np)%id=1
                     ! Set position
                     this%ls%p(np)%pos=[grid%xm(i),grid%ym(j),grid%zm(k)]
                     ! Assign velocity
                     this%ls%p(np)%vel=[0.0_WP,0.0_WP,0.0_WP]
                     ! Locate the particle on the mesh
                     this%ls%p(np)%ind=this%ls%cfg%get_ijk_global(this%ls%p(np)%pos,[this%ls%cfg%imin,this%ls%cfg%jmin,this%ls%cfg%kmin])
                     ! Assign a unique integer to particle
                     this%ls%p(np)%i=np
                     ! Activate the particle
                     this%ls%p(np)%flag=0
                     
                     ! Deactivate solver on top and bottom to enable boundary conditions
                     if (abs(this%ls%p(np)%pos(2)).gt.0.45_WP*height) this%ls%p(np)%id=0

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
         this%ens_out=ensight(cfg=this%cfg,name='solid')
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
         this%mfile=monitor(this%cfg%amRoot,'simulation_solid')
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
		class(simsolid), intent(inout) :: this
      
      ! Increment time
      call this%ls%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Enforce boundary conditions
      apply_bc: block
         integer :: n
         do n=1,this%ls%np_
            if (this%ls%p(n)%id.eq.0) then
               if (this%ls%p(n)%pos(2).gt.0.0_WP) this%ls%p(n)%pos(2)=this%ls%p(n)%pos(2)+1.0e-3_WP*this%time%t
               if (this%ls%p(n)%pos(2).lt.0.0_WP) this%ls%p(n)%pos(2)=this%ls%p(n)%pos(2)-1.0e-3_WP*this%time%t
            end if
         end do
      end block apply_bc
      
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
		class(simsolid), intent(inout) :: this
		
		! Deallocate work arrays
		
	end subroutine final
   
   
end module simsolid_class