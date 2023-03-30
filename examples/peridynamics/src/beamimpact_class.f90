!> Definition of a solid simulation superclass that drives an lss_class object
module beamimpact_class
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
   
   public :: beamimpact

   !> Solid simulation object
   type :: beamimpact
      
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
      procedure, private :: postproc       !< Postprocess beam position
      procedure :: init                    !< Initialize solid simulation
      procedure :: step                    !< Advance solid simulation by one time step
      procedure :: final                   !< Finalize solid simulation
   end type beamimpact
   
   
contains
   
   
   
   ! Postprocess beam position
   subroutine postproc(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      class(beamimpact), intent(inout) :: this
      integer :: n,ierr
      real(WP) :: myx1,myx2
      !integer :: myn1,myn2,n1,n2
      !myn1=0; myn2=0
      myx1=0.0_WP; myx2=0.0_WP
      do n=1,this%ls%np_
         !if (this%ls%p(n)%id.eq.1) then
         !   myn1=myn1+1; myx1=myx1+this%ls%p(n)%pos(1)
         !else if (this%ls%p(n)%id.eq.2) then
         !   myn2=myn2+1; myx2=myx2+this%ls%p(n)%pos(1)
         !end if
         if (this%ls%p(n)%i.eq.this%ibeam1) myx1=this%ls%p(n)%pos(1)
         if (this%ls%p(n)%i.eq.this%ibeam2) myx2=this%ls%p(n)%pos(1)
      end do
      !call MPI_ALLREDUCE(myn1,n1,1,MPI_INTEGER,MPI_SUM,this%ls%cfg%comm,ierr)
      call MPI_ALLREDUCE(myx1,this%Xbeam1,1,MPI_REAL_WP,MPI_SUM,this%ls%cfg%comm,ierr)
      !if (n1.gt.0) then
      !   this%Xbeam1=this%Xbeam1/real(n1,WP)
      !else
      !   this%Xbeam1=0.0_WP
      !end if
      !call MPI_ALLREDUCE(myn2,n2,1,MPI_INTEGER,MPI_SUM,this%ls%cfg%comm,ierr)
      call MPI_ALLREDUCE(myx2,this%Xbeam2,1,MPI_REAL_WP,MPI_SUM,this%ls%cfg%comm,ierr)
      !if (n2.gt.0) then
      !   this%Xbeam2=this%Xbeam2/real(n2,WP)
      !else
      !   this%Xbeam2=0.0_WP
      !end if
   end subroutine postproc
   

   !> Initialization of solid simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(beamimpact), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_beamimpact')
      
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
         use mpi_f08,     only: MPI_BCAST,MPI_INTEGER
         integer :: i,j,k,nx,ny,np,ierr
         real(WP) :: length,width,dx
         real(WP), dimension(:), allocatable :: x,y
         type(sgrid) :: grid

         ! Create solver
         this%ls=lss(cfg=this%cfg,name='solid')
         
         ! Discretize cylinder
         call this%input%read('Solid dx',dx)
         call this%input%read('Bar length',length)
         call this%input%read('Bar width',width)

         ! Set material properties
         this%ls%elastic_modulus=75.0e9_WP
         this%ls%poisson_ratio=0.25_WP
         this%ls%rho=2700.0_WP
         this%ls%crit_energy=1.0e9_WP
         this%ls%dV=dx**3
         this%ls%delta=3.0_WP*dx
         
         ! Only root process initializes solid particles
         if (this%ls%cfg%amRoot) then
            ! Create simple rectilinear grid
            nx=2*int(length/dx)
            ny=int(width/dx)
            allocate(x(1:nx+1),y(1:ny+1))
            do i=1,nx+1
               x(i)=real(i-1,WP)*dx-length
            end do
            do j=1,ny+1
               y(j)=real(j-1,WP)*dx-0.5_WP*width
            end do
            ! General serial grid object
            grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=y,xper=.false.,yper=.false.,zper=.false.,name='elements')
            ! Loop over mesh and create particles
            np=0
            do k=1,ny
               do j=1,ny
                  do i=1,nx
                     ! Increment particle
                     np=np+1
                     call this%ls%resize(np)
                     ! Set position
                     this%ls%p(np)%pos=[grid%xm(i),grid%ym(j),grid%zm(k)]
                     ! Set object id and velocity
                     if (grid%xm(i).gt.0.0_WP) then
                        this%ls%p(np)%id=1
                        this%ls%p(np)%vel=[-10.0_WP,0.0_WP,0.0_WP]
                     else
                        this%ls%p(np)%id=2
                        this%ls%p(np)%vel=[+10.0_WP,0.0_WP,0.0_WP]
                     end if
                     ! Zero out force
                     this%ls%p(np)%Abond=0.0_WP
                     ! Locate the particle on the mesh
                     this%ls%p(np)%ind=this%ls%cfg%get_ijk_global(this%ls%p(np)%pos,[this%ls%cfg%imin,this%ls%cfg%jmin,this%ls%cfg%kmin])
                     ! Assign a unique integer to particle
                     this%ls%p(np)%i=np
                     ! Activate the particle
                     this%ls%p(np)%flag=0
                     ! ID the COM of both beams
                     if (i.eq.1*nx/4.and.j.eq.ny/2.and.k.eq.ny/2) this%ibeam1=np
                     if (i.eq.3*nx/4.and.j.eq.ny/2.and.k.eq.ny/2) this%ibeam2=np
                  end do
               end do
            end do
         end if
         call MPI_BCAST(this%ibeam1,1,MPI_INTEGER,0,this%ls%cfg%comm,ierr)
         call MPI_BCAST(this%ibeam2,1,MPI_INTEGER,0,this%ls%cfg%comm,ierr)
         
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
         this%ens_out=ensight(cfg=this%cfg,name='beamimpact')
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
         call this%postproc()
         ! Create simulation monitor
         this%mfile=monitor(this%cfg%amRoot,'beamimpact')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%xbeam1,'X beam1')
         call this%mfile%add_column(this%xbeam2,'X beam2')
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
		class(beamimpact), intent(inout) :: this
      
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
      call this%postproc()
      call this%mfile%write()
      
   end subroutine step
   
   
   !> Finalize solid simulation
   subroutine final(this)
		implicit none
		class(beamimpact), intent(inout) :: this
		
		! Deallocate work arrays
		
	end subroutine final
   
   
end module beamimpact_class