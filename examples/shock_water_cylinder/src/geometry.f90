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
      implicit none
      type(sgrid) :: grid


      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         integer :: cpd,cpl,np,nyr,nxr
         real(WP) :: rdx,dx,box_x1,box_x2,box_y1,dcyl,xcyl
         real(WP) :: Lcalc,r,rold,err,tol
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         dx = Lx/nx
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('nz',nz,default=1); allocate(z(nz+1))

         if (nz.eq.1) then
           Lz = dx
         else
           call param_read('Lz',Lz)
         end if
         ! Read in droplet information
         call param_read('Cylinder diameter',dcyl)
         call param_read('Cylinder location',xcyl)

         if (param_exists('Cells per diameter')) then
           ! Stretched grid
           call param_read('Cells per diameter',cpd)
           call param_read('Cells per left region',cpl)
           call param_read('Refined region left bdy',box_x1)
           call param_read('Refined region right bdy',box_x2)
           call param_read('Refined region top bdy',box_y1,default=0.0_WP)

           ! Refined cell size
           rdx = dcyl/cpd

           !! Left x region: 0 to x1    !!
           ! check if stretching is possible, given inputs
           if ((rdx*cpl) .ge. box_x1) then
             print*,"using uniform spacing in stretched region meets or exceeds the boundary"
             print*,"x1:",box_x1,"rdx*cpl",rdx*cpl
             print*,"geometric series is impossible. please alter inputs"
             return
           end if
           tol = 1e-10_WP ! tolerance for how close calculated Lx should be to real Lx
           ! initial values for loop to find r
           err = 1.0_WP
           r = 1.1_WP
           i = 0
           ! loop to find r
           do while (err.gt.tol)
             i = i+1
             rold = r ! store r from last iteration
             r = (1.0_WP-(1.0_WP-rold)*(box_x1+rdx)/rdx)**(1.0_WP/(real(cpl+1,WP))) ! calculate new r
             Lcalc = -rdx+rdx*(1.0_WP-r**(real(cpl+1,WP)))/(1.0_WP-r) ! use new r to calc Lx
             err = abs(Lcalc-box_x1) ! find err
           end do
           print*,'Left x stretching ratio',r
           x(1) = 0.0_WP
           x(cpl+1) = box_x1
           ! use r to populate x
           do i = cpl,1,-1
             x(i) = x(i+1)-rdx*r**(cpl+1-i)
           end do

           !! Middle x region: x1 to x2 !!
           nxr = ceiling((box_x2-box_x1)/rdx)
           do i=2,nxr+1
             x(cpl+i) = x(cpl+i-1)+rdx
           end do
           box_x2 = x(cpl+nxr+1)

           !! Right x region: x2 to Lx  !!
           np = nx-nxr-cpl
           ! check if stretching is possible, given inputs
           if ((box_x2+rdx*np) .ge. Lx) then
             print*,"using uniform spacing in stretched region meets or exceeds the boundary"
             print*,"Lx:",Lx,"x2+rdx*np",box_x2+rdx*np
             print*,"geometric series is impossible. please alter inputs"
             return
           end if
           tol = 1e-10_WP ! tolerance for how close calculated Lx should be to real Lx
           ! initial values for loop to find r
           err = 1.0_WP
           r = 1.1_WP
           i = 0
           ! loop to find r
           do while (err.gt.tol)
             i = i+1
             rold = r ! store r from last iteration
             r = (1.0_WP-(1.0_WP-rold)*(Lx-box_x2+rdx)/rdx)**(1.0_WP/real(np+1,WP)) ! calculate new r
             Lcalc = box_x2-rdx+rdx*(1.0_WP-r**real(np+1,WP))/(1.0_WP-r) ! use new r to calc Lx
             err = abs(Lcalc-Lx) ! find err
           end do
           print*,'Right x stretching ratio',r
           ! use r to populate x
           do i = nxr+cpl+2,nx+1
             x(i) = x(i-1)+rdx*r**(i-nxr-cpl-1)
           end do


           ! Check if y is meant to be stretched
           if (box_y1.eq.0.0_WP) then
             ! Use uniform spacing for y
             y(ny/2+1) = 0.0_WP
             ! Use inner resolution for all cells
             do j = 2,ny/2+1
               y(ny/2+j) = y(ny/2+j-1)+rdx
             end do
             ! All points are uniform
             nyr = ny/2-1
             print*,"No stretching in y. Ly:",Ly,"rdx*ny",rdx*ny

           else
             !! Middle y region: 0 to y1  !!
             nyr = ceiling(box_y1/rdx)
             y(ny/2+1) = 0.0_WP
             do j=2,nyr+1
               y(ny/2+j) = y(ny/2+j-1)+rdx
             end do
             box_y1 = y(ny/2+nyr+1)

             !! Top y region: y1 to Ly/2  !!
             np = ny/2-nyr
             ! check if stretching is possible, given inputs
             if ((box_y1+rdx*np) .ge. Ly/2) then
               print*,"using uniform spacing in stretched region meets or exceeds the boundary"
               print*,"Ly/2:",Ly/2,"y1+rdx*np",box_y1+rdx*np
               print*,"geometric series is impossible. please alter inputs"
               return
             end if
             tol = 1e-10_WP ! tolerance for how close calculated Ly should be to real Ly
             ! initial values for loop to find r
             err = 1.0_WP
             r = 1.1_WP
             i = 0
             ! loop to find r
             do while (err.gt.tol)
               i = i+1
               rold = r ! store r from last iteration
               r = (1.0_WP-(1.0_WP-rold)*(Ly/2.0_WP-box_y1+rdx)/rdx)**(1.0_WP/real(np+1,WP)) ! calculate new r
               Lcalc = box_y1-rdx+rdx*(1.0_WP-r**real(np+1,WP))/(1.0_WP-r) ! use new r to calc Lx
               err = abs(Lcalc-Ly/2.0_WP) ! find err
             end do
             print*,'Top y stretching ratio',r
             ! use r to populate y
             do j = nyr+2,ny/2+1
               y(ny/2+j) = y(ny/2+j-1)+rdx*r**(j-nyr-1)
             end do
           end if

           !! Mirror y for bottom half: -Ly/2 to 0  !!
           do j = 1,ny/2
             y(j) = -y(ny+2-j)
           end do

           ! Print mesh data to check
           print*,'Left stretched region'
           print*,'first dx',x(2)-x(1),'last dx',x(cpl+1)-x(cpl)
           print*,'first pt',x(1),'last point',x(cpl+1),'x1',box_x1
           print*,'Uniform x region'
           print*,'dx',x(cpl+2)-x(cpl+1),'rdx',rdx
           print*,'last pt',x(cpl+nxr+1),'number of pts',nxr
           print*,'Right stretched region'
           print*,'first dx',x(cpl+nxr+2)-x(cpl+nxr+1),'last dx',x(nx+1)-x(nx)
           print*,'last pt',x(nx+1)
           print*,'Uniform y region'
           print*,'dy',y(ny/2+2)-y(ny/2+1),'rdx',rdx
           print*,'first pt',y(ny/2+1),'number of pts',nyr
           print*,'Top stretched region'
           print*,'first dy',y(ny/2+nyr+2)-y(ny/2+nyr+1),'last dy',y(ny+1)-y(ny)
           print*,'first pt',y(ny/2+nyr+1),'last pt',y(ny+1)
           print*,'Bottom region'
           print*,'first pt',y(1),'last pt',y(ny/2)

           ! Assuming everything is spaced correctly, Make boundary points exact
           x(1) = 0.0_WP
           x(nx+1) = Lx
           y(1) = -Ly/2.0_WP
           y(ny+1) = Ly/2.0_WP
           ! This is just to eliminate round-off error, not to 'fix' any wrong calcs

         else
           ! Uniform grid
           do i=1,nx+1
             x(i) = real(i-1,WP)*dx
           end do
           do j=1,ny+1
             y(j) = real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
           end do
         end if

         ! z is always uniform
         do k=1,nz+1
           z(k) = real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='ShockWaterCylinder')

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
