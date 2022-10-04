!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private

   !> Single config
   type(config), public :: cfg

   public :: geometry_init

contains


   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read, param_exists
      use parallel,    only: amRoot
      use messager,    only: die
      implicit none
      type(sgrid) :: grid


      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,jj
         integer :: cpd,n_prev,n_mid,n_post,n_wall
         real(WP) :: rdx,box_xl,box_xr,box_yl,box_yr,box_z,djet
         real(WP) :: rx1,rx2,ry1,ry2,rz
         real(WP) :: Lx,Ly,Lz
         real(WP) :: resx_1,resx_2,resy_1,resy_2,dymin
         real(WP), dimension(:), allocatable :: x,y,z
         logical :: is_sx,is_sy,is_sz,zper
         
         ! Get mesh dimensions
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)

         ! Incorporate stretching
         is_sx = param_exists('Stretching ratio in x pre')
         is_sy = param_exists('Stretching ratio in y post')
         is_sz = param_exists('Stretching ratio in z')

         ! Jet properties
         call param_read('Jet diameter',djet)
         call param_read('Wall cells in domain', n_wall, default=1)
         
         ! Ask for desired resolution
         if (is_sx.or.is_sy.or.is_sz)  then
            call param_read('Cells per jet diameter', cpd)
            rdx = djet/cpd
         end if
         ! rdx = refined dx, used for refined region when mesh stretching

         ! Initialize markers
         resx_1 = -1.0_WP; resy_1 = -1.0_WP
         ! -1 is used as an unrealistic value to be replaced

         ! Create the grid
         ! Grid in X
         if (is_sx) then
            
            ! Get stretching parameters
            call param_read('Stretching ratio in x pre',rx1)
            call param_read('Stretching ratio in x post',rx2)
            call param_read('Stretching bdy in x pre',box_xl)
            call param_read('Stretching bdy in x post',box_xr)
            ! Calculate variables for mesh generation
            ! Number of points in each region
            n_mid = ceiling((box_xr-box_xl)/rdx)
            n_prev = nint(log(-(box_xl+rdx)*(1.0_WP-rx1)/rdx+1.0_WP)/(log(rx1)+epsilon(1.0_WP)))-1
            n_post = nint(log(-(Lx-box_xr+rdx)*(1.0_WP-rx2)/rdx+1.0_WP)/(log(rx2)+epsilon(1.0_WP)))-1
            if (n_prev.lt.0) n_prev = 0
            if (n_post.lt.0) n_post = 0
            nx = n_prev+n_mid+n_post

            if (amRoot) print*, 'nx',nx,'n_prev', n_prev,'n_mid',n_mid,'n_post',n_post,'in x direction'
            !print*, 'Lx should become ', rdx*(1.0_WP-rx2**(n_post+1))/(1.0_WP-rx2+epsilon(1.0_WP))+rdx*(1.0_WP-rx1**(n_prev+1))/(1.0_WP-rx1+epsilon(1.0_WP))-box_xl+box_xr-2*rdx
            
            allocate(x(nx+1))

            ! First section (from stretched to refined region)
            do i=1,n_prev+1
               x(i) = box_xl - rdx*(rx1-rx1**(-i+n_prev+2))/(1.0_WP-rx1+epsilon(1.0_WP))
            end do
            resx_1 = x(i-1) ! store for print output, used for next section as well
            ! Second section (refined region)
            do i=n_prev+2,n_mid+n_prev+1
               x(i) = resx_1+real(i-1-n_prev,WP)*rdx
            end do
            resx_2 = x(i-1) ! store for print output, used for next section as well
            ! Third section (from refined region to stretched)
            do i=n_prev+n_mid+2,nx+1
               x(i) = resx_2 + rdx*(rx2-rx2**(i-n_mid-n_prev))/(1.0_WP-rx2+epsilon(1.0_WP))
            end do

            ! Store new Lx
            Lx = x(nx+1)-x(1)
            !print*, 'Lx, x(1), x(nx+1)', Lx, x(1), x(nx+1)
            if (amRoot) print*, 'dx as resolved as cpd from', resx_1,   'to', resx_2,  ', for', (resx_2-resx_1)/Lx*100.0_WP, '% of Lx'
         else
            ! Uniform mesh
            call param_read('nx',nx); 
            allocate(x(nx+1))
            do i=1,nx+1
               x(i) = real(i-1,WP)*Lx/real(nx,WP)
            end do
         end if

         ! Grid in Y
         if (is_sy) then
            call param_read('Cell height at wall',dymin)
            call param_read('Stretching ratio in y pre',ry1)
            call param_read('Stretching ratio in y post',ry2)
            call param_read('Stretching bdy in y post',box_yr)

            !! - Calculate number of points needed in mesh in order to allocate - !!

            ! ZEROTH REGION : wall cells

            ! FIRST REGION : from specified resolution at wall to uniform "refined box" resolution
            ! Nearest number of cells to accomplish the stretching
            ! -- this is set so that the n_prev+1 cell has dy=rdx
            n_prev = nint(log(rdx/dymin)/log(ry1))
            ! New stretching rate to correspond with number of cells
            ry1 = (rdx/dymin)**(1.0_WP/real(n_prev+1,WP))
            ! Length of new region
            box_yl = dymin*(1.0_WP-ry1**n_prev)/(1.0_WP-ry1)

            ! SECOND REGION : uniform box of a specified resolution
            n_mid = max(0,ceiling((box_yr-box_yl)/rdx))
            ! Recalculate end of refined region
            box_yr = box_yl + rdx*n_mid

            ! THIRD REGION : stretching according to a specified rate until domain length is exceeded
            n_post = nint(log(-(Ly-(box_yr-rdx))*(1.0_WP-ry2)/rdx+1.0_WP)/(log(ry2)+epsilon(1.0_WP)))-1

            if (n_post.lt.0) n_post = 0
            ny = n_wall + n_prev + n_mid + n_post

            if (amRoot) print*, 'ny',ny,'n_wall',n_wall,'n_prev',n_prev,'n_mid',n_mid,'n_post',n_post,'in y direction'
            !print*, 'Ly should become ', rdx*(1.0_WP-ry2**(n_post+1))/(1.0_WP-ry2+epsilon(1.0_WP))+box_yr-rdx+dymin*n_wall
            
            allocate(y(ny+1))
            
            !! - Calculate actual cell heights and locations - !!
            ! ZEROTH REGION
            do j=1,n_wall+1
               y(j) = real(j-1-n_wall,WP)*dymin
            end do
            ! FIRST REGION
            do j=n_wall+2,n_wall+n_prev+1
               jj = j-(n_wall+2)
               y(j) = y(j-1) + dymin*ry1**jj
            end do
            ! SECOND REGION
            do j=n_wall+n_prev+2,n_wall+n_prev+n_mid+1
               y(j) = y(j-1) + rdx
            end do
            ! THIRD REGION
            do j=n_wall+n_prev+n_mid+2,n_wall+n_prev+n_mid+n_post+1
               jj = j-(n_wall+n_prev+n_mid+1)
               y(j) = y(j-1) + rdx*ry2**jj
            end do 

            !! - Final calculations and variables for printing - !!
            Ly = y(ny+1)-y(1)
            !print*, 'Ly, y(1), y(ny+1)', Ly, y(1), y(ny+1)
            resy_1 = y(n_wall+n_prev+1); resy_2 = y(n_wall+n_prev+n_mid+1)
            if (amRoot) print*, 'dy as resolved as cpd from', resy_1,   'to', resy_2,  ', for', (resy_2-resy_1)/Ly*100.0_WP, '% of Ly'
          else
            call param_read('ny',ny); 
            allocate(y(ny+1))
            do j=1,ny+1
               y(j) = real(j-1-n_wall,WP)*Ly/real(ny,WP)
            end do
         end if

         ! Grid in Z
         if (is_sz.and.nz.ne.1) then
            call param_read('Lz',Lz)
            
            ! Get stretching parameters
            call param_read('Stretching ratio in z',rz)
            call param_read('Stretching bdy in z',box_z)
            ! Number of points in each region
            n_mid = ceiling(2.0_WP*box_z/rdx)/2*2
            box_z = n_mid*rdx
            n_prev = nint(log(-(Lz/2.0_WP-box_z/2.0_WP+rdx)*(1.0_WP-rz)/rdx+1.0_WP)/(log(rz)+epsilon(1.0_WP)))-1
            if (n_prev.lt.0) n_prev = 0
            n_post = n_prev
            nz = n_prev+n_mid+n_post

            if (amRoot) print*, 'nz',nz,'n_prev', n_prev,'n_mid',n_mid,'n_post',n_post,'in z direction'
            !print*, 'Lz should become ', box_z - 2.0_WP*rdx + 2.0_WP*rdx*(1.0_WP-rz**(n_post+1))/(1.0_WP-rz+epsilon(1.0_WP))

            allocate(z(nz+1))
            
            ! First section (from stretched to refined region)
            do k=1,n_prev+1
               z(k) = -box_z/2.0_WP - rdx*(rz-rz**(-k+n_prev+2))/(1.0_WP-rz+epsilon(1.0_WP))
            end do
            ! Second section (refined region)
            do k=n_prev+2,n_mid+n_prev+1
               z(k) = -box_z/2.0_WP + real(k-1-n_prev,WP)*rdx
            end do
            ! Third section (from refined region to stretched)
            do k=n_prev+n_mid+2,nz+1
               z(k) = box_z/2.0_WP + rdx*(rz-rz**(k-n_mid-n_prev))/(1.0_WP-rz+epsilon(1.0_WP))
            end do

            ! Store new Lz
            Lz = z(nz+1)-z(1)
            !print*, 'Lz, z(1), z(nz+1)', Lz, z(1), z(nz+1)
            if (amRoot) print*, 'dz as resolved as cpd from', -box_z, 'to', box_z, ', for', (    box_z    )/Lz*100.0_WP, '% of Lz'
         else
            ! Uniform mesh
            call param_read('nz',nz); 
            if (nz.eq.1) then
               Lz = minval(x(2:nx+1)-x(1:nx))
            else
               call param_read('Lz',Lz)
            end if
            allocate(z(nz+1))
            do k=1,nz+1
               z(k) = real(k-1,WP)*Lz/real(nz,WP)-0.5_WP*Lz
            end do
         end if
         
         if (is_sx.or.is_sy.or.is_sz) then
           if (amRoot) print*, 'total number of cells = ', nx*ny*nz, 'domain cell dimensions', nx,'x',ny,'x',nz
         end if

         ! General serial grid object
         zper = .false.; if (nz.eq.1) zper=.true.;
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=zper,name='LJSC')

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


   end subroutine geometry_init


end module geometry
