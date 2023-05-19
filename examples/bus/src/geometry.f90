!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class,   only: config
   use precision,      only: WP
   use partmesh_class, only: partmesh
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
   !> Object to contain all of the bus dimensions
   type :: geom_bus
      ! Units conversion
      real(WP) :: i2m
      ! Mesh resolution factor
      integer  :: gres
      ! Dimensions of bus
      real(WP) :: l_bus,h_bus,w_bus
      ! Misc parameters
      real(WP) :: x_frontoffset,min_length
      ! Number of seat components
      integer  :: nseat_driver,nseat_curb1,nseat_curb2
      ! Size of seat components
      real(WP) :: l_seat,w_seat,th_seat,h_bseat,th_bseat
      ! Spacing of seats
      real(WP) :: l_seatrow,w_seatcol,h_seat,l_seatgap
      real(WP) :: x_seatdriver,x_seatcurb,x_seatcurb2
      ! Dimensions of parcel rack
      real(WP) :: h_rack,w_rack
      ! Knee racks and objects at boundaries
      real(WP) :: w_closet,x_closet,x_kneedriver,x_kneecurb
      real(WP) :: x_lav,w_galley,th_knee,h_knee
      ! Table
      real(WP) :: w_table,l_table,y_table,h_table
      ! Vent dimensions
      ! intakes to system (outlets of domain)
      real(WP) :: h_vvent,w_vvent,w_fventfront
      real(WP) :: l_fventfront,l_fventdriver,l_fventcurb
      real(WP) :: h_fventback,h_ventlav,w_ventlav
      real(WP) :: w_vinrack,l_vinrack
      ! outputs from system (inlets of domain)
      real(WP) :: w_ventwindow
      ! Vent locations
      real(WP) :: y_vvent,z_fventback,y_ventwindow
      real(WP) :: l_ventrow,z_ventface1,z_ventface2
      real(WP) :: x_v2facedriver,x_v2facecurb
      real(WP) :: w_ventface,l_ventface,x_vinrack
      ! Vent velocities
      real(WP) :: vel_vvent,vel_fventfront
      real(WP) :: vel_fventdriver,vel_fventcurb
      real(WP) :: vel_ventwindow
      real(WP) :: vel_ventface,vel_vinrack
      ! Passenger height
      real(WP) :: h_psg
   end type geom_bus
   
   type(geom_bus), public :: bus
   
   !> Passenger information
   type(partmesh)             , public :: psg_mesh
   real(WP)                   , public :: h_psg    !> Height of passenger mouth from the ground
   integer , parameter        , public :: npsg=32  !< Hard code 32 passengers
   integer , dimension(3,npsg), public :: ipsg     !< Index for each passenger
   real(WP), dimension(3,npsg), public :: xpsg     !< Location for each passenger
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      ! Fill bus object with desired dimensions
      populate_bus: block
         ! Units conversion
         bus%i2m = 0.0254_WP
         ! Dimensions of bus
         bus%l_bus = 45.0_WP*12.0_WP*bus%i2m
         bus%h_bus = 78.5_WP*bus%i2m ! floor to ceiling
         bus%w_bus = 2.59_WP
         ! Misc parameters
         bus%x_frontoffset = 1.2_WP
         bus%min_length = 2.0_WP*bus%i2m
         ! Number of seat components
         bus%nseat_driver = 10
         bus%nseat_curb1 = 3
         bus%nseat_curb2 = 8
         ! Size of seat components
         bus%l_seat = 0.6867_WP
         bus%w_seat = 0.5559_WP
         bus%th_seat = 7.0_WP*bus%i2m
         bus%h_bseat = 28.0_WP*bus%i2m
         bus%th_bseat = 5.0_WP*bus%i2m
         ! Spacing of seats
         bus%l_seatrow = 35.0_WP*bus%i2m
         bus%w_seatcol = 0.654_WP
         bus%h_seat = 20.0_WP*bus%i2m
         bus%l_seatgap = 0.8284_WP
         bus%x_seatdriver = 1.9838_WP
         bus%x_seatcurb = 1.4388_WP
         bus%x_seatcurb2= bus%x_seatcurb &
         & + real(bus%nseat_curb1-1,WP)*bus%l_seatrow &
         & + bus%l_seat + bus%l_seatgap
         ! Dimensions of parcel rack
         bus%h_rack = 61.0_WP*bus%i2m
         bus%w_rack = 0.78_WP
         ! Knee racks and objects at boundaries
         bus%w_closet = 1.09_WP
         bus%x_closet = 2.1473_WP
         bus%x_kneedriver = 2.6923_WP
         bus%x_kneecurb = 1.4497_WP
         bus%x_lav = 0.8502_WP
         bus%w_galley = 0.8066_WP
         bus%th_knee = 5.0_WP*bus%i2m
         bus%h_knee = bus%h_seat+0.5_WP*bus%h_bseat
         ! Table
         bus%w_table = 2.0_WP*bus%w_seatcol
         bus%l_table = 1.5_WP*bus%l_seat
         bus%y_table = bus%h_seat+0.3_WP*bus%h_bseat
         bus%h_table = 2.0_WP*bus%i2m
         ! Vent dimensions
         ! intakes to system (outlets of domain)
         bus%h_vvent = 13.0_WP*bus%i2m
         bus%w_vvent = 3.0_WP*bus%i2m
         bus%w_fventfront = 18.0_WP*bus%i2m
         bus%l_fventfront = 5.0_WP*bus%i2m
         bus%l_fventdriver = 40_WP*bus%i2m
         bus%l_fventcurb = 62.0_WP*bus%i2m
         bus%h_fventback = 4.0_WP*bus%i2m
         bus%h_ventlav = 5.0_WP*bus%i2m
         bus%w_ventlav = 2.0_WP*bus%w_seat
         bus%w_vinrack = 2.0_WP*bus%i2m
         bus%l_vinrack = bus%l_seatrow
         ! outputs from system (inlets of domain)
         bus%w_ventwindow = 2.0_WP*bus%i2m
         ! Vent locations
         bus%y_vvent = 10.0_WP*bus%i2m
         bus%z_fventback = 4.0_WP*bus%i2m
         bus%y_ventwindow = bus%h_knee
         bus%l_ventrow = 35.0_WP*bus%i2m
         bus%z_ventface1 = 19.25_WP*bus%i2m
         bus%z_ventface2 = 24.75_WP*bus%i2m
         bus%x_v2facedriver = 12.0_WP*bus%i2m
         bus%x_v2facecurb = 10.0_WP*bus%i2m
         bus%w_ventface = bus%z_ventface2-bus%z_ventface1
         bus%l_ventface = 0.5_WP*bus%w_ventface
         bus%x_vinrack = 0.5_WP*bus%l_bus
         ! Vent velocities (magnitude and sign)
         bus%vel_vvent = 1.1_WP
         bus%vel_fventfront = -0.9_WP
         bus%vel_fventdriver = -0.7_WP
         bus%vel_fventcurb = 0.5_WP
         bus%vel_ventwindow = 0.6_WP
         bus%vel_ventface = -2.0_WP
         bus%vel_vinrack = 3.0_WP ! this value gets replaced
      end block populate_bus
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nn,nx,ny,nz,nstretch
         real(WP) :: Lx,Ly,Lz
         real(WP) :: rdx,lstretch,ndiv
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Set up grid dimensions according to bus geometry
         Lx = bus%l_bus-bus%x_lav-bus%x_frontoffset
         Ly = bus%h_bus
         Lz = bus%w_bus
         
         ! Read in grid resolution multiple
         call param_read('Grid resolution',bus%gres)
         rdx = bus%min_length/real(bus%gres,WP)
         
         ! Determine required number of cells for target resolution
         nx = floor(Lx/rdx)
         ny = floor(Ly/rdx)
         nz = floor(Lz/rdx)
         !print*, nx,ny,nz
         
         ! Allocate grid size
         allocate(x(nx+1))
         allocate(y(ny+1))
         allocate(z(nz+1))
         
         ! Number of cells incorporated into stretching
         nstretch = 3
         ! number to divide by
         ndiv = 0.0_WP
         do i=1,nstretch
            ndiv = ndiv+i
         end do
         
         ! To fit dimensions, stretching is allowed in x direction at the end
         x(1) = bus%x_lav ! start with offset from lavatory
         ! uniform region
         do i=2,nx+1-nstretch
            x(i)=x(i-1)+rdx
         end do
         lstretch = Lx-real(nx,WP)*rdx ! remainder of length
         ! stretched region
         do i=nx+2-nstretch,nx+1
            nn = nx+1-nstretch
            x(i)=x(i-1)+rdx+real(i-nn,WP)/ndiv*lstretch
         end do
         
         ! Stretch at top in y direction
         y(1) = 0.0_WP ! start at floor
         ! uniform region
         do j=2,ny+1-nstretch
            y(j)=y(j-1)+rdx
         end do
         lstretch = Ly-real(ny,WP)*rdx ! remainder of length
         ! stretched region
         do j=ny+2-nstretch,ny+1
            nn = ny+1-nstretch
            y(j)=y(j-1)+rdx+real(j-nn,WP)/ndiv*lstretch
         end do
         
         ! Stretch near center in z direction
         z(1) = 0.0_WP ! start at back (driver side)
         ! uniform region 1
         do k=2,nz/2+1-nstretch
            z(k)=z(k-1)+rdx
         end do
         lstretch = Lz-real(nz,WP)*rdx ! remainder of length
         ! stretched region - increasing
         do k=nz/2+2-nstretch,nz/2+1
            nn = nz/2+1-nstretch
            z(k)=z(k-1)+rdx+0.5_WP*real(k-nn,WP)/ndiv*lstretch
         end do
         ! stretched region - decreasing
         do k=nz/2+2,nz/2+1+nstretch
            nn = nz/2+2+nstretch
            z(k)=z(k-1)+rdx+0.5_WP*real(nn-k,WP)/ndiv*lstretch
         end do
         ! uniform region 2
         do k=nz/2+2+nstretch,nz+1
            z(k)=z(k-1)+rdx
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='c2c_bus')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Adjust quantities in bus object to be mesh-aligned
      ! concerned with smaller or unevenly repeating objects
      ! not concerned with even number of inches
      adjust_bus: block
         real (WP) :: dim_factor,x1face, adist
         integer :: n
         adist = real(mod(nint(bus%l_seatrow/bus%i2m),nint(bus%min_length/bus%gres/bus%i2m)),WP)*bus%i2m
         ! Stretching occurs in the front (x), at the top (y), and in the middle (z)
         dim_factor = bus%min_length/real(bus%gres,WP)
         !! -- Make dimensions of thin/repeating objects align with cell lengths -- !!
         bus%th_seat = real(ceiling(bus%th_seat/dim_factor),WP)*dim_factor
         bus%th_bseat = real(ceiling(bus%th_bseat/dim_factor),WP)*dim_factor
         bus%l_seat = real(ceiling(bus%l_seat/dim_factor),WP)*dim_factor
         bus%w_seatcol = real(ceiling(bus%w_seatcol/dim_factor),WP)*dim_factor
         bus%th_knee = real(ceiling(bus%th_knee/dim_factor),WP)*dim_factor
         bus%w_vvent = real(ceiling(bus%w_vvent/dim_factor),WP)*dim_factor
         bus%l_fventfront = real(ceiling(bus%l_fventfront/dim_factor),WP)*dim_factor
         bus%h_ventlav = real(ceiling(bus%h_ventlav/dim_factor),WP)*dim_factor
         bus%z_ventface1 = real(ceiling(bus%z_ventface1/dim_factor),WP)*dim_factor
         bus%z_ventface2 = bus%z_ventface1+real(ceiling((bus%z_ventface2-bus%z_ventface1)/dim_factor),WP)*dim_factor
         bus%l_ventface = real(nint(bus%l_ventface/dim_factor),WP)*dim_factor
         !! -- Make locations of thin/repeating objects align with cell lengths -- !!
         bus%l_seatgap = real(ceiling(bus%l_seatgap/dim_factor),WP)*dim_factor
         bus%h_rack = real(ceiling(bus%h_rack/dim_factor),WP)*dim_factor
         bus%y_table = real(ceiling(bus%y_table/dim_factor),WP)*dim_factor
         bus%x_seatdriver = bus%x_lav+real(nint((bus%x_seatdriver-bus%x_lav)/dim_factor),WP)*dim_factor
         bus%x_seatcurb = bus%x_lav+real(nint((bus%x_seatcurb-bus%x_lav)/dim_factor),WP)*dim_factor
         ! These are measured from front, adjust according to uniform mesh from back
         bus%x_kneedriver = bus%l_bus-bus%x_lav-real(ceiling((bus%l_bus-bus%x_kneedriver-bus%x_lav)/dim_factor),WP)*dim_factor
         bus%x_kneecurb = bus%l_bus-bus%x_lav-real(ceiling((bus%l_bus-bus%x_kneecurb-bus%x_lav)/dim_factor),WP)*dim_factor
         !! -- Align parcel rack intake with mesh so that area is easily found -- !!
         bus%x_vinrack = bus%x_lav+real(nint((bus%x_vinrack-bus%x_lav)/dim_factor),WP)*dim_factor
         bus%l_vinrack = real(2*ceiling(0.5_WP*bus%l_vinrack/dim_factor),WP)*dim_factor
         !! -- Calculate resulting other dimensions -- !!
         bus%w_ventface = bus%z_ventface2-bus%z_ventface1
         bus%x_seatcurb2= bus%x_seatcurb &
         & + real(bus%nseat_curb1-1,WP)*bus%l_seatrow &
         & + bus%l_seat + bus%l_seatgap
         ! bus%l_seatrow and bus%l_ventrow are odd number of inches, but I will handle that elsewhere
         
         !! -- Enforce mass balance in the parcel rack -- !!
         ! (Forces the total flowrate of the rack intakes to match the outputs)
         ! (But does assume unrealistic symmetry)
         bus%vel_vinrack = sign(1.0_WP,bus%vel_vinrack)*abs(bus%vel_ventface)   * &
         &                 (bus%z_ventface2-bus%z_ventface1)*bus%l_ventface     * &
         &             real(bus%nseat_driver+bus%nseat_curb1+bus%nseat_curb2,WP)/ &
         &                 (2.0_WP*bus%w_vinrack*bus%l_vinrack)
         
      end block adjust_bus
      
      
      ! Create masks for this config
      create_walls: block
         integer  :: i,j,k,n
         real(WP) :: adist,x1face
         
         ! get factor to adjust seat spacing according to mesh
         adist = real(mod(nint(bus%l_seatrow/bus%i2m),nint(bus%min_length/bus%gres/bus%i2m)),WP)*bus%i2m
         
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Back of bus (left)
                  if (i.lt.cfg%imin) cfg%VF(i,j,k)=0.0_WP
                  ! Front of bus (right)
                  if (i.gt.cfg%imax) cfg%VF(i,j,k)=0.0_WP
                  ! Floor
                  if (j.lt.cfg%jmin) cfg%VF(i,j,k)=0.0_WP
                  ! Ceiling
                  if (j.gt.cfg%jmax) cfg%VF(i,j,k)=0.0_WP
                  ! Windows on driver side (back)
                  if (k.lt.cfg%kmin) cfg%VF(i,j,k)=0.0_WP
                  ! Windows on curb side (front)
                  if (k.gt.cfg%kmax) cfg%VF(i,j,k)=0.0_WP
                  ! Parcel rack - driver side
                  if (cfg%x(i)  .ge.bus%x_seatdriver-100.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%x(i+1).le.bus%x_seatdriver+bus%l_seatrow*real(bus%nseat_driver,WP).and.&
                  &   cfg%y(j)  .gt.bus%h_rack-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_rack+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Parcel rack - curb side
                  if (cfg%x(i)  .ge.bus%x_seatcurb.and.&
                  &   cfg%x(i+1).le.bus%x_seatcurb2+bus%l_seatrow*real(bus%nseat_curb2,WP).and.&
                  &   cfg%y(j)  .gt.bus%h_rack-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.bus%w_bus-bus%w_rack-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_bus+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Panel for back floor intake vent - driver side
                  if (cfg%x(i)  .gt.bus%x_seatdriver.and.&
                  &   cfg%x(i+1).lt.bus%x_seatdriver+bus%l_fventdriver.and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%h_fventback+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%z_fventback+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Panel for back floor intake vent - curb side
                  if (cfg%x(i)  .gt.bus%x_seatcurb+2.0_WP*bus%l_seatrow.and.&
                  &   cfg%x(i+1).lt.bus%x_seatcurb+2.0_WP*bus%l_seatrow+bus%l_fventcurb.and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%h_fventback+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.bus%w_bus-bus%z_fventback-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_bus+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Knee wall - driver side
                  if (cfg%x(i)  .gt.bus%l_bus-bus%x_kneedriver.and.&
                  &   cfg%x(i+1).lt.bus%l_bus-bus%x_kneedriver+bus%th_knee.and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%h_knee+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_seatcol+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Knee wall - curb side (connected to vertical vent)
                  if (cfg%x(i)  .gt.bus%l_bus-bus%x_kneecurb.and.&
                  &   cfg%x(i+1).lt.bus%l_bus-bus%x_kneecurb+bus%th_knee.and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%h_knee+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.bus%w_bus-bus%w_seatcol-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_bus+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Wall below window vents - driver side
                  if (cfg%x(i)  .gt.bus%x_seatdriver.and.&
                  &   cfg%x(i+1).lt.bus%l_bus-bus%x_closet.and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%y_ventwindow+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_ventwindow+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Wall below window vents - curb side
                  if (cfg%x(i)  .gt.bus%x_seatcurb.and.&
                  &   cfg%x(i+1).lt.bus%l_bus-bus%x_kneecurb.and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%y_ventwindow+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.bus%w_bus-bus%w_ventwindow-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_bus+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Table
                  if (cfg%x(i)  .gt.bus%x_seatcurb2-0.5_WP*bus%l_seatgap-0.5_WP*bus%l_table.and.&
                  &   cfg%x(i+1).lt.bus%x_seatcurb2-0.5_WP*bus%l_seatgap+0.5_WP*bus%l_table.and.&
                  &   cfg%y(j)  .gt.bus%y_table-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%y_table+bus%h_table+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.bus%w_bus-bus%w_table-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_bus+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Galley
                  if (cfg%x(i)  .gt.bus%x_lav-100.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%x(i+1).lt.bus%x_seatdriver+100.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_galley+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Closet
                  if (cfg%x(i)  .gt.bus%l_bus-bus%x_closet.and.&
                  &   cfg%x(i+1).lt.bus%l_bus+10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k)  .gt.-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%z(k+1).lt.bus%w_closet+10.0_WP*epsilon(1.0_WP)) cfg%VF(i,j,k)=0.0_WP
                  ! Seats
                  ! same y range
                  if (cfg%y(j)  .gt.bus%h_seat-bus%th_seat-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%h_seat+10.0_WP*epsilon(1.0_WP)) then
                     ! Driver side
                     ! set z range
                     if (cfg%z(k)  .gt.bus%w_seatcol-bus%w_seat-100.0_WP*epsilon(1.0_WP).and.&
                     &   cfg%z(k+1).lt.bus%w_seatcol+100.0_WP*epsilon(1.0_WP)) then
                        ! initial x location
                        x1face = bus%x_seatdriver
                        ! cycle through rows in x
                        do n=1,bus%nseat_driver
                           if (cfg%x(i)  .gt.x1face-100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist .and.&
                           &   cfg%x(i+1).lt.x1face+100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+bus%l_seat) then
                           cfg%VF(i,j,k)=0.0_WP
                           end if
                        end do
                     end if
                     ! Curb side
                     ! set z range for both columns
                     if ((cfg%z(k)  .gt.bus%w_bus-bus%w_seatcol-100.0_WP*epsilon(1.0_WP).and.&
                     &    cfg%z(k+1).lt.bus%w_bus-bus%w_seatcol+bus%w_seat+100.0_WP*epsilon(1.0_WP)).or.&
                     &   (cfg%z(k)  .gt.bus%w_bus-2.0_WP*bus%w_seatcol-100.0_WP*epsilon(1.0_WP).and.&
                     &    cfg%z(k+1).lt.bus%w_bus-2.0_WP*bus%w_seatcol+bus%w_seat+100.0_WP*epsilon(1.0_WP))) then
                        ! back section: initial x location
                        x1face = bus%x_seatcurb
                        do n=1,bus%nseat_curb1
                           if (cfg%x(i)  .gt.x1face-100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist.and.&
                           &   cfg%x(i+1).lt.x1face+100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+bus%l_seat) then
                              cfg%VF(i,j,k)=0.0_WP
                           end if
                        end do
                        ! front section: initial x location (includes backward seats)
                        x1face = bus%x_seatcurb2
                        do n=1,bus%nseat_curb2
                           if (cfg%x(i)  .gt.x1face-100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist.and.&
                           &   cfg%x(i+1).lt.x1face+100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+bus%l_seat) then
                              cfg%VF(i,j,k)=0.0_WP
                           end if
                        end do
                     end if
                  end if
                  ! Seat backs
                  ! same y range
                  if (cfg%y(j)  .gt.bus%h_seat-10.0_WP*epsilon(1.0_WP).and.&
                  &   cfg%y(j+1).lt.bus%h_seat+bus%h_bseat+10.0_WP*epsilon(1.0_WP)) then
                     ! Driver side
                     ! set z range
                     if (cfg%z(k)  .gt.bus%w_seatcol-bus%w_seat-100.0_WP*epsilon(1.0_WP).and.&
                     &   cfg%z(k+1).lt.bus%w_seatcol+100.0_WP*epsilon(1.0_WP)) then
                        ! initial x location
                        x1face = bus%x_seatdriver
                        ! cycle through rows in x
                        do n=1,bus%nseat_driver
                           if (cfg%x(i)  .gt.x1face-100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist .and.&
                           &   cfg%x(i+1).lt.x1face+100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+bus%th_bseat) then
                              cfg%VF(i,j,k)=0.0_WP
                           end if
                        end do
                     end if
                     ! Curb side
                     ! set z range for both columns
                     if ((cfg%z(k)  .gt.bus%w_bus-bus%w_seatcol-100.0_WP*epsilon(1.0_WP).and.&
                     &    cfg%z(k+1).lt.bus%w_bus-bus%w_seatcol+bus%w_seat+100.0_WP*epsilon(1.0_WP)).or.&
                     &   (cfg%z(k)  .gt.bus%w_bus-2.0_WP*bus%w_seatcol-100.0_WP*epsilon(1.0_WP).and.&
                     &    cfg%z(k+1).lt.bus%w_bus-2.0_WP*bus%w_seatcol+bus%w_seat+100.0_WP*epsilon(1.0_WP))) then
                        ! back section: initial x location
                        x1face = bus%x_seatcurb
                        do n=1,bus%nseat_curb1
                           if (cfg%x(i)  .gt.x1face-100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist.and.&
                           &   cfg%x(i+1).lt.x1face+100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+bus%th_bseat) then
                              cfg%VF(i,j,k)=0.0_WP
                           end if
                        end do
                        ! front section: initial x location with offset for backward seats
                        x1face = bus%x_seatcurb2+bus%l_seat-bus%th_bseat
                        do n=1,bus%nseat_curb2
                           if (cfg%x(i)  .gt.x1face-100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist.and.&
                           &   cfg%x(i+1).lt.x1face+100.0_WP*epsilon(1.0_WP) &
                           &   +real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+bus%th_bseat) then
                              cfg%VF(i,j,k)=0.0_WP
                           end if
                           ! remove offset from x1face after first step
                           x1face = bus%x_seatcurb2
                        end do
                     end if
                  end if
                  
               end do
            end do
         end do
      end block create_walls
      
      
      ! Create passengers
      add_passengers: block
         integer :: n,count
         real(WP) :: adist,x1face
         
         ! Get factor to adjust seat spacing according to mesh
         adist = real(mod(nint(bus%l_seatrow/bus%i2m),nint(bus%min_length/bus%gres/bus%i2m)),WP)*bus%i2m
         
         ! Read in source height
         call param_read('Source height',h_psg)
         
         ! Assign positions
         count=0
         ! Driver side
         x1face=bus%x_seatdriver
         do n=1,bus%nseat_driver
            count=count+1
            xpsg(:,count)=[x1face+real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+0.5_WP*bus%l_seat,h_psg,bus%w_seatcol-0.5_WP*bus%w_seat]
            ipsg(:,count)=cfg%get_ijk_global(xpsg(:,count),[cfg%imino,cfg%jmino,cfg%kmino])
         end do
         ! Curb side - back section
         x1face=bus%x_seatcurb
         do n=1,bus%nseat_curb1
            count=count+1
            xpsg(:,count)=[x1face+real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+0.5_WP*bus%l_seat,h_psg,bus%w_bus-bus%w_seatcol+0.5_WP*bus%w_seat]
            ipsg(:,count)=cfg%get_ijk_global(xpsg(:,count),[cfg%imino,cfg%jmino,cfg%kmino])
            count=count+1
            xpsg(:,count)=[x1face+real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+0.5_WP*bus%l_seat,h_psg,bus%w_bus-2.0_WP*bus%w_seatcol+0.5_WP*bus%w_seat]
            ipsg(:,count)=cfg%get_ijk_global(xpsg(:,count),[cfg%imino,cfg%jmino,cfg%kmino])
         end do
         ! Curb side - front section
         x1face=bus%x_seatcurb2+0.5_WP*bus%l_seat-2.0_WP*bus%th_bseat
         do n=1,bus%nseat_curb2
            count=count+1
            xpsg(:,count)=[x1face+real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+0.5_WP*bus%l_seat,h_psg,bus%w_bus-bus%w_seatcol+0.5_WP*bus%w_seat]
            ipsg(:,count)=cfg%get_ijk_global(xpsg(:,count),[cfg%imino,cfg%jmino,cfg%kmino])
            count=count+1
            xpsg(:,count)=[x1face+real(n-1,WP)*bus%l_seatrow+real(mod(n-1,2),WP)*adist+0.5_WP*bus%l_seat,h_psg,bus%w_bus-2.0_WP*bus%w_seatcol+0.5_WP*bus%w_seat]
            ipsg(:,count)=cfg%get_ijk_global(xpsg(:,count),[cfg%imino,cfg%jmino,cfg%kmino])
            x1face=bus%x_seatcurb2
         end do
         
         ! While we're at it, create a particle mesh for plotting the passenger positions
         psg_mesh=partmesh(nvar=1,nvec=0,name='psg'); psg_mesh%varname(1)='id'; call psg_mesh%reset()
         call psg_mesh%set_size(npsg)
         do n=1,npsg
            psg_mesh%pos(:,n)=xpsg(:,n)
            psg_mesh%var(1,n)=real(n,WP)
         end do
         
      end block add_passengers
      
      
   end subroutine geometry_init
   
   
end module geometry
