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
      integer :: i,j,k,nx,ny,nz
      integer :: cpd,cpl,np,nyr,nxr,nzr
      real(WP) :: rdx,dx,box_x1,box_x2,box_y1,ddrop,box_z1
      real(WP), dimension(3) :: dctr
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
      call param_read('Droplet diameter',ddrop)
      call param_read('Droplet location',dctr)

      if (param_exists('Cells per diameter')) then
        ! Stretched grid
        if (amRoot) print*,"===== Stretched Mesh Description ====="
        call param_read('Cells per diameter',cpd)
        call param_read('Cells per left region',cpl)
        call param_read('Refined region left bdy',box_x1)
        call param_read('Refined region right bdy',box_x2)
        call param_read('Refined region top y bdy',box_y1,default=0.0_WP)
        call param_read('Refined region top z bdy',box_z1,default=0.0_WP)

        ! Refined cell size
        rdx = ddrop/cpd

        !! Left x region: 0 to x1    !!
        ! check if stretching is possible, given inputs
        if ((rdx*cpl) .ge. box_x1) then
          if (amRoot) then
            print*,"MESH ERROR: using uniform spacing in initial stretched x region meets or exceeds the boundary"
            print*,"- prescribed region length x1:",box_x1,"region length with uniform spacing:",rdx*cpl
          end if
          call die("[geometry] geometric series is impossible (0). please alter inputs")
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
        if (amRoot) print*,'Left x stretching ratio',r
        x(1) = 0.0_WP
        x(cpl+1) = box_x1
        ! use r to populate x
        do i = cpl,1,-1
          x(i) = x(i+1)-rdx*r**(cpl+1-i)
        end do

        !! Middle x region: x1 to x2 !!
        nxr = ceiling((box_x2-box_x1)/rdx)
        if (nx .le. nxr+cpl) then
          if (amRoot) then
            print*,"MESH ERROR: insufficient number of points for x direction"
            print*,"= prescribed total points",nx
            print*,"= left stretched pts",cpl,"refined pts",nxr,"right stretched pts",nx-nxr-cpl
          end if
          call die("[geometry] geometric series is impossible (1). please alter inputs")
        end if
        do i=2,nxr+1
          x(cpl+i) = x(cpl+i-1)+rdx
        end do
        box_x2 = x(cpl+nxr+1)

        !! Right x region: x2 to Lx  !!
        np = nx-nxr-cpl
        ! check if stretching is possible, given inputs
        if ((box_x2+rdx*np) .ge. Lx) then
          if (amRoot) then
            print*,"MESH ERROR: using uniform spacing in final x stretched region meets or exceeds the boundary"
            print*,"= prescribed length Lx:",Lx,"length with uniform spacing:",box_x2+rdx*np
          end if
          call die("[geometry] geometric series is impossible (2). please alter inputs")
        end if
        ! Need to make these abort statements
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
        if (amRoot) print*,'Right x stretching ratio',r
        ! use r to populate x
        do i = nxr+cpl+2,nx+1
          x(i) = x(i-1)+rdx*r**(i-nxr-cpl-1)
        end do
        ! Assuming everything is spaced correctly, Make boundary points exact
        x(1) = 0.0_WP
        x(nx+1) = Lx


        if (ny.eq.1) then
           ! Check if y is a single cell
           y(1) = -0.5 * rdx
           y(2) = 0.5 * rdx

        elseif (box_y1.eq.0.0_WP) then
          ! If y is not meant to be stretched, use uniform spacing for y
          y(ny/2+1) = 0.0_WP
          ! Use inner resolution for all cells
          do j = 2,ny/2+1
            y(ny/2+j) = y(ny/2+j-1)+rdx
          end do
          ! All points are uniform
          nyr = ny/2-1
          if (amRoot) print*,"No stretching in y. Ly:",Ly,"rdx*ny",rdx*ny

        else
          !! Middle y region: 0 to y1  !!
          nyr = ceiling(box_y1/rdx)
          if (ny/2 .le. nyr) then
            if (amRoot) then
              print*,"MESH ERROR: insufficient number of points for y direction"
              print*,"= prescribed total points",ny
              print*,"= refined pts",2*nyr,"right stretched pts",ny-2*nyr
            end if
            call die("[geometry] geometric series is impossible (3). please alter inputs")
          end if
          y(ny/2+1) = 0.0_WP
          do j=2,nyr+1
            y(ny/2+j) = y(ny/2+j-1)+rdx
          end do
          box_y1 = y(ny/2+nyr+1)

          !! Top y region: y1 to Ly/2  !!
          np = ny/2-nyr
          ! check if stretching is possible, given inputs
          if ((box_y1+rdx*np) .ge. Ly/2) then
            if (amRoot) then
              print*,"MESH ERROR: using uniform spacing in stretched y regions meets or exceeds the boundaries"
              print*,"= prescribed half height Ly/2:",Ly/2,"half height with uniform spacing",box_y1+rdx*np
            end if
            call die("[geometry] geometric series is impossible (4). please alter inputs")
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
            if (amRoot) print*,'Top y stretching ratio',r
            ! use r to populate y
            do j = nyr+2,ny/2+1
               y(ny/2+j) = y(ny/2+j-1)+rdx*r**(j-nyr-1)
            end do
            ! Assuming everything is spaced correctly, Make boundary points exact
            y(1) = -Ly/2.0_WP
            y(ny+1) = Ly/2.0_WP
        end if

        !! Mirror y for bottom half: -Ly/2 to 0  !!
        do j = 1,ny/2
           y(j) = -y(ny+2-j)
        end do

        if (nz.eq.1) then
          ! Check if y is a single cell
          z(1) = -0.5 * rdx
          z(2) = 0.5 * rdx

        elseif (box_z1.eq.0.0_WP) then
          ! If z is not meant to be stretched, use uniform spacing for z
          z(nz/2+1) = 0.0_WP
          ! Use inner resolution for all cells
          do k = 2,nz/2+1
             z(nz/2+k) = z(nz/2+k-1)+rdx
          end do
          ! All points are uniform
          nzr = nz/2-1
          if (amRoot) print*,"No stretching in z. Lz:",Lz,"rdx*nz",rdx*nz

        else
          !! Middle z region: 0 to z1  !!
          nzr = ceiling(box_z1/rdx)
          if (nz/2 .le. nzr) then
            if (amRoot) then
              print*,"MESH ERROR: insufficient number of points for z direction"
              print*,"= prescribed total points",nz
              print*,"= refined pts",2*nzr,"right stretched pts",nz-2*nzr
            end if
            call die("[geometry] geometric series is impossible (3). please alter inputs")
          end if
          z(nz/2+1) = 0.0_WP
          do k=2,nzr+1
            z(nz/2+k) = z(nz/2+k-1)+rdx
          end do
          box_z1 = z(nz/2+nzr+1)

          !! Top y region: z1 to Lz/2  !!
          np = nz/2-nzr
          ! check if stretching is possible, given inputs
          if ((box_z1+rdx*np) .ge. Lz/2) then
            if (amRoot) then
              print*,"MESH ERROR: using uniform spacing in stretched z regions meets or exceeds the boundaries"
              print*,"= prescribed half height Lz/2:",Lz/2,"half height with uniform spacing",box_z1+rdx*np
            end if
            call die("[geometry] geometric series is impossible (4). please alter inputs")
          end if
          tol = 1e-10_WP ! tolerance for how close calculated Lz should be to real Lz
          ! initial values for loop to find r
          err = 1.0_WP
          r = 1.1_WP
          i = 0
          ! loop to find r
          do while (err.gt.tol)
            i = i+1
            rold = r ! store r from last iteration
            r = (1.0_WP-(1.0_WP-rold)*(Lz/2.0_WP-box_z1+rdx)/rdx)**(1.0_WP/real(np+1,WP)) ! calculate new r
            Lcalc = box_z1-rdx+rdx*(1.0_WP-r**real(np+1,WP))/(1.0_WP-r) ! use new r to calc Lx
            err = abs(Lcalc-Lz/2.0_WP) ! find err
          end do
          if (amRoot) print*,'Top z stretching ratio',r
          ! use r to populate z
          do k = nzr+2,nz/2+1
            z(nz/2+k) = z(nz/2+k-1)+rdx*r**(k-nzr-1)
          end do
          ! Assuming everything is spaced correctly, Make boundary points exact
          z(1) = -Lz/2.0_WP
          z(nz+1) = Lz/2.0_WP
        end if

        !! Mirror y for bottom half: -Lz/2 to 0  !!
        do k = 1,nz/2
          z(k) = -z(nz+2-k)
        end do

        ! Print mesh data to check
        if (amRoot) then
          print*,'Left stretched region'
          print*,'    first dx',x(2)-x(1),'last dx',x(cpl+1)-x(cpl)
          print*,'    first pt',x(1),'last point',x(cpl+1),'x1',box_x1
          print*,'Uniform x region'
          print*,'    dx',x(cpl+2)-x(cpl+1),'rdx',rdx
          print*,'    last pt',x(cpl+nxr+1),'number of pts',nxr
          print*,'Right stretched region'
          print*,'    first dx',x(cpl+nxr+2)-x(cpl+nxr+1),'last dx',x(nx+1)-x(nx)
          print*,'    last pt',x(nx+1)
          print*,'Uniform y region'
          print*,'    dy',y(ny/2+2)-y(ny/2+1),'rdx',rdx
          print*,'    first pt',y(ny/2+1),'number of pts',nyr
          print*,'Top stretched region'
          print*,'    first dy',y(ny/2+nyr+2)-y(ny/2+nyr+1),'last dy',y(ny+1)-y(ny)
          print*,'    first pt',y(ny/2+nyr+1),'last pt',y(ny+1)
          print*,'Bottom region'
          print*,'    first pt',y(1),'last pt',y(ny/2-nyr+1)
          print*,'Uniform z region'
          print*,'    dz',z(nz/2+2)-z(nz/2+1),'rdx',rdx
          print*,'    first pt',z(nz/2+1),'number of pts',nzr
          print*,'Top stretched region'
          print*,'    first dy',z(nz/2+nzr+2)-z(nz/2+nzr+1),'last dy',z(nz+1)-z(nz)
          print*,'    first pt',z(nz/2+nzr+1),'last pt',z(nz+1)
          print*,'Bottom region'
          print*,'    first pt',z(1),'last pt',z(ny/2-nyr+1)
          print*,' '
        end if

        else
          ! Uniform grid
          do i=1,nx+1
            x(i) = real(i-1,WP)*dx
          end do
          do j=1,ny+1
            y(j) = real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
          end do
          do k=1,nz+1
            z(k) = real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
          end do
        end if


        ! General serial grid object
        grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='ShockDrop')

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
