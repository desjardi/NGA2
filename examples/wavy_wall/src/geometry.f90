!> Various definitions and tools for initializing NGA2 config
module geometry
  use ibconfig_class, only: ibconfig
   use precision,     only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg

   !> Domain size
   real(WP) :: Lx,Ly,Lz
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: dy
         real(WP), dimension(:), allocatable :: x,y,z
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid in x and z, shift y by 3 cells for IB
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny-1
            y(j+2)=real(j-1,WP)/real(ny-2,WP)*Ly
         end do
         ! Add 2 cells below for the IB
         dy=y(4)-y(3)
         y(2)=y(3)-dy
         y(1)=y(2)-dy
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='channel')
      end block create_grid
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      ! Create masks for this config
      ! Taken from Almeida et al. 1993 "Wake flows behind two-dimensional model hills"
      ! Breuer et al. (2009). Flow over periodic hills–numerical and experimental study in a wide range of Reynolds numbers. Computers & Fluids, 38(2), 433-457.
      create_walls: block
        use ibconfig_class, only: bigot,sharp
        integer :: i,j,k
        real(WP) :: dist,H,h_sep,dy_min
        real(WP), parameter :: h_crest=28.0_WP
        call param_read('Hill height',H)
        h_sep = Lx/H * h_crest
        if (Lx.lt.4.0_WP*H) then
           print *, 'Lx must be > 4x hill height!'
           stop
        end if
        cfg%VF=1.0_WP
        do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 if (cfg%xm(i)*h_crest/H.le.0.5_WP*h_sep) then
                    dist = cfg%xm(i)*h_crest/H
                 else
                    dist = abs(cfg%xm(i)*h_crest/H-h_sep)
                 end if
                 if (dist.le.9.0_WP) then
                    cfg%Gib(i,j,k) = H/h_crest*min(h_crest, h_crest+6.77507096985e-3_WP*dist**2-2.1245277758e-3_WP*dist**3)-cfg%ym(j)
                 else if (dist.gt.9.0_WP .and. dist.le.14.0_WP) then
                    cfg%Gib(i,j,k) = H/h_crest*(25.07355893131_WP+9.754803562315E-01_WP*dist-1.016116352781E-01_WP*dist**2+1.889794677828E-03_WP*dist**3)-cfg%ym(j)
                 else if (dist.gt.14.0_WP .and. dist.le.20.0_WP) then 
                    cfg%Gib(i,j,k) = H/h_crest*(2.579601052357E+01_WP+8.206693007457E-01_WP*dist -9.055370274339E-02_WP*dist**2+1.626510569859E-03_WP*dist**3)-cfg%ym(j)
                 else if (dist.gt.20.0_WP .and. dist.le.30.0_WP) then
                    cfg%Gib(i,j,k) = H/h_crest*(4.046435022819E+01_WP-1.379581654948E+00_WP*dist +1.945884504128E-02_WP*dist**2-2.070318932190E-04_WP*dist**3)-cfg%ym(j)
                 else if (dist.gt.30.0_WP .and. dist.le.40.0_WP) then
                    cfg%Gib(i,j,k) = H/h_crest*( 1.792461334664E+01_WP+8.743920332081E-01_WP*dist-5.567361123058E-02_WP*dist**2+6.277731764683E-04_WP*dist**3)-cfg%ym(j)          
                 else if (dist.gt.40.0_WP .and. dist.le.54.0_WP) then 
                    cfg%Gib(i,j,k) = H/h_crest*max(0.,5.639011190988E+01_WP-2.010520359035E+00_WP*dist+1.644919857549E-02_WP*dist**2+2.674976141766E-05_WP*dist**3)-cfg%ym(j)
                 else 
                    cfg%Gib(i,j,k) = -cfg%ym(j)
                 end if
              end do
           end do
        end do
        ! Get normal vector
        call cfg%calculate_normal()
        ! Get VF field
        call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
         do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 if (cfg%ym(j).lt.cfg%y(cfg%jmin).or.cfg%ym(j).gt.cfg%y(cfg%jmax+1)) cfg%VF(i,j,k)=0.0_WP
              end do
           end do
        end do
      end block create_walls
      
   end subroutine geometry_init
   
   
end module geometry
