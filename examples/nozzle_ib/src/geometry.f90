!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,      only: WP
   use surfmesh_class, only: surfmesh
   use ensight_class,  only: ensight
   implicit none
   private
   
   !> IB config
   type(ibconfig), public :: cfg
   
   !> Surface mesh
   type(surfmesh), public :: plymesh
   
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
         real(WP) :: Lx,Ly,Lz,xshift
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1)); call param_read('X shift',xshift)
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='nozzle')
         
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
      

      ! Read in the PLY geometry
      read_ply: block
         use string, only: str_medium
         use param,  only: param_read
         character(len=str_medium) :: plyfile

         ! Read in ply filename
         call param_read('ply filename',plyfile)

         ! Create surface meshg from ply
         plymesh=surfmesh(plyfile=plyfile,nvar=0,name='ply')

      end block read_ply
      
      
      ! Create IB walls for this config
      create_walls: block
         use ibconfig_class, only: sharp
         use messager,       only: die
         use mathtools,      only: cross_product,normalize
         use irl_fortran_interface
         real(WP), parameter :: safe_coeff=3.0_WP
         integer :: i,j,k,np,nv,iv,ip
         real(WP) :: mydist
         real(WP), dimension(3) :: pos,nearest_pt,mynearest,mynorm
         type(Poly_type), dimension(:), allocatable :: poly
         real(WP), dimension(:,:), allocatable :: bary,vert,nvec
         
         ! Preprocess surface mesh data using IRL
         allocate(poly(1:plymesh%nPoly))
         allocate(vert(1:3,1:maxval(plymesh%polySize)))
         allocate(bary(1:3,1:plymesh%nPoly))
         allocate(nvec(1:3,1:plymesh%nPoly))
         do np=1,plymesh%nPoly
            ! Allocate polygon
            call new(poly(np))
            ! Fill it up
            do nv=1,plymesh%polySize(np)
               iv=sum(plymesh%polySize(1:np-1))+nv
               vert(:,nv)=[plymesh%xVert(plymesh%polyConn(iv)),plymesh%yVert(plymesh%polyConn(iv)),plymesh%zVert(plymesh%polyConn(iv))]
            end do
            call construct(poly(np),plymesh%polySize(np),vert(1:3,1:plymesh%polySize(np)))
            mynorm=normalize(cross_product(vert(:,2)-vert(:,1),vert(:,3)-vert(:,2)))
            mydist=dot_product(mynorm,vert(:,1))
            call setPlaneOfExistence(poly(np),[mynorm(1),mynorm(2),mynorm(3),mydist])
            ! Also store its barycenter
            bary(:,np)=calculateCentroid(poly(np))
            nvec(:,np)=calculateNormal  (poly(np))
         end do
         deallocate(vert)

         ! Create IB distance field using IRL
         cfg%Gib=huge(1.0_WP)
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Store cell center position
                  pos=[cfg%xm(i),cfg%ym(j),cfg%zm(k)]
                  ! Traverse all polygons
                  do np=1,plymesh%nPoly
                     ! Calculate distance to centroid
                     nearest_pt=pos-bary(:,np)
                     mydist=dot_product(nearest_pt,nearest_pt)
                     ! If close enough, compute exact distance to the polygon instead
                     if (mydist.lt.(safe_coeff*cfg%min_meshsize)**2) then
                        nearest_pt=calculateNearestPtOnSurface(poly(np),pos)
                        nearest_pt=pos-nearest_pt
                        mydist=dot_product(nearest_pt,nearest_pt)
                     end if
                     ! Remember closest distance
                     if (mydist.lt.cfg%Gib(i,j,k)) then
                        cfg%Gib(i,j,k)=mydist
                        mynearest=nearest_pt
                        ip=np
                     end if
                  end do
                  ! Take the square root
                  cfg%Gib(i,j,k)=sqrt(cfg%Gib(i,j,k))
                  ! Find the sign
                  if (dot_product(mynearest,nvec(:,ip)).gt.0.0_WP) cfg%Gib(i,j,k)=-cfg%Gib(i,j,k)
               end do
            end do
         end do
         deallocate(bary,nvec,poly)
         
         ! Get normal vector
         call cfg%calculate_normal()
         
         ! Get VF field
         call cfg%calculate_vf(method=sharp,allow_zero_vf=.true.)
         
      end block create_walls
      
      
   end subroutine geometry_init
   

end module geometry
