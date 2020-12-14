module mms_geom
   use precision, only: WP
   
   interface cube_refine
      module procedure cube_refine_vol              ! compute volume and area
      module procedure cube_refine_vfunc_afunc_fd   ! compute volume and area integrals - finite difference
      module procedure cube_refine_afunc_fd         ! compute area integrals - finite difference
   end interface cube_refine
   interface rec_refine
      module procedure rec_refine_area      ! compute area
      module procedure rec_refine_afunc_fd  ! compute area integral - finite difference
   end interface rec_refine
   
   real(WP), parameter :: r16=1.0_WP/6.0_WP
   real(WP), parameter :: r13=1.0_WP/3.0_WP
   
   
   ! ==== A triangle :
   !
   !    - vertex labels (vx), edge labels (ex) and contribution to case number [1,2,4]
   !
   !              v1  -  +1
   !             /  \
   !            /    \
   !           e3     e1
   !          /        \
   !         /          \
   !  +4 - v3-----e2-----v2  - +2
   !
   ! =====
   
   ! ==== lookup table for edge to vertices
   integer, parameter, dimension(2,3) :: tri_edge_to_vertex=reshape([1,2,2,3,1,3],shape(tri_edge_to_vertex))
   
   !==== lookup table for case to intersected edges :
   !     0 is a sentinel value indicating no more edges
   !     ordering of the edges important for area calculations
   integer, parameter, dimension(2,0:7) :: tri_case_to_edge=reshape([0,0,1,3,1,2,3,2,2,3,1,2,1,3,0,0],shape(tri_case_to_edge))
   
   !==== lookup table for case to vertices with G > 0 :
   !     0 is a sentinel value
   !     ordering is important for decomposing quads into tris
   integer, parameter, dimension(2,0:7) :: tri_case_to_vertex=reshape([0,0,1,0,2,0,1,2,3,0,1,3,2,3,0,0],shape(tri_case_to_vertex))
   !==== number of vertices with G > 0 (non-zero entries above) per case
   integer, parameter, dimension(0:7) :: tri_nverts=[0,1,1,2,1,2,2,3]
   
   type tri
      real(WP), dimension(3,3) :: vertex
      real(WP), dimension(3) :: G
      ! triangle area
      real(WP) :: area
      ! intersection points (1 tuple per edge)
      real(WP), dimension(3,3) :: intp
      ! intersected area with G > 0
      real(WP) :: A_intersect
      ! centroid of intersected area
      real(WP), dimension(3) :: C_intersect
      ! centroid of triangle
      real(WP), dimension(3) :: T_centroid
      integer :: case
   end type tri
   
   !==== A Rectangle *** lookup table based on this vertex ordering ***
   !
   !       vertices
   !
   !      3---------4
   !      |         |
   !      |         |
   !      |         |
   !      1---------2
   !
   !=====
   
   !===== lookup table for rec vertices to tri vertices
   integer, parameter, dimension(3,2) :: rec_to_tri=reshape([1,4,2,1,4,3],shape(rec_to_tri))
   
   
   !==== A Tetrahedron:
   !
   !      -vertex labels (vx), edge labels (ex) and contribution to case number [1,2,4,8]
   !
   !              v1  - +1
   !             /|\
   !            / | \
   !           e4 |  e1
   !          /   e5  \
   !         /    |    \
   !  +8 - v4--e6-|-----v2  - +2
   !         \    |    /
   !          \   |   /
   !           e3 |  e2
   !            \ | /
   !             \|/
   !              v3 - + 4
   !======
   
   
   ! ==== lookup table for edge to vertices
   integer, parameter, dimension(2,6) :: edge_to_vertex=reshape([1,2,2,3,3,4,1,4,1,3,2,4],shape(edge_to_vertex))
   
   
   !==== lookup table for case to intersected edges :
   !     0 is a sentinel value indicating no more edges
   !     ordering of the edges for each case are important for area/volume calcs
   integer, parameter, dimension(4,0:15) :: case_to_edge=reshape([0,0,0,0, &                     ! case 0
   &                                                              1,4,5,0, &                     ! case 1
   &                                                              1,2,6,0, &                     ! case 2
   &                                                              2,6,5,4, &                     ! case 3
   &                                                              2,5,3,0, &                     ! case 4
   &                                                              3,2,4,1, &                     ! case 5
   &                                                              6,1,3,5, &                     ! case 6
   &                                                              6,4,3,0, &                     ! case 7
   &                                                              3,4,6,0, &                     ! case 8
   &                                                              5,1,3,6, &                     ! case 9
   &                                                              1,2,4,3, &                     ! case 10
   &                                                              2,3,5,0, &                     ! case 11
   &                                                              2,5,6,4, &                     ! case 12
   &                                                              2,1,6,0, &                     ! case 13
   &                                                              1,5,4,0, &                     ! case 14
   &                                                              0,0,0,0],shape(case_to_edge))  ! case 15
   
   !==== number of intersected edges (non-zero entries above) per case
   integer, parameter, dimension(0:15) :: nedges=[0,3,3,4,3,4,4,3,3,4,4,3,4,3,3,0]
   
   !==== lookup table for case to vertices with G > 0 :
   !     0 is a sentinel value
   !     ordering is important for decomposing polyhedra into tets
   integer, parameter, dimension(4,0:15) :: case_to_vertex=reshape([0,0,0,0, &                       ! case 0
   &                                                                1,0,0,0, &                       ! case 1
   &                                                                2,0,0,0, &                       ! case 2
   &                                                                2,1,0,0, &                       ! case 3
   &                                                                3,0,0,0, &                       ! case 4
   &                                                                3,1,0,0, &                       ! case 5
   &                                                                2,3,0,0, &                       ! case 6
   &                                                                2,1,3,0, &                       ! case 7
   &                                                                4,0,0,0, &                       ! case 8
   &                                                                1,4,0,0, &                       ! case 9
   &                                                                2,4,0,0, &                       ! case 10
   &                                                                2,4,1,0, &                       ! case 11
   &                                                                3,4,0,0, &                       ! case 12
   &                                                                3,1,4,0, &                       ! case 13
   &                                                                2,3,4,0, &                       ! case 14
   &                                                                1,2,3,4],shape(case_to_vertex))  ! case 15
   
   !==== number of vertices with G > 0 (non-zero entries above) per case
   integer, parameter, dimension(0:15) :: nverts=[0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4]
   
   type tet
      real(WP), dimension(3,4) :: vertex
      real(WP), dimension(4) :: G
      ! tet volume
      real(WP) :: Volume
      ! intersection points (1 tuple per edge)
      real(WP), dimension(3,6) :: intp
      ! Intersected volume and area with G > 0
      real(WP) :: V_intersect,A_intersect
      ! C_intersect : centroid of intersected volumne
      real(WP), dimension(3) :: C_intersect
      ! Tet and intersecting area Centroid
      real(WP), dimension(3) :: T_centroid,A_centroid
      integer :: case
   end type tet
   
   !==== A Hexahedron *** lookup table based on this vertex ordering ***
   !
   !       vertices
   !
   !      7---------8
   !     /|        /|
   !    / |       / |
   !   /  |      /  |
   !  5---3-----6---4
   !  |  /      |  /
   !  | /       | /
   !  |/        |/
   !  1---------2
   !
   !=====
   
   !===== lookup table for hex vertices to tet vertices
   integer, parameter, dimension(4,6) :: hex_to_tet=reshape([1,2,4,8,1,2,6,8,1,3,4,8,1,3,7,8,1,5,6,8,1,5,7,8],shape(hex_to_tet))
   
   !==== cube_shift is applied to cube vertices to generate 8 sub-cubes with
   !     appropriate vertices
   real(WP), parameter, dimension(3,8,8) :: cube_shift=reshape([ &
   +0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 1, vertex  :1
   -0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 1, vertex  :2
   +0.0_WP,-0.5_WP, 0.0_WP, &     ! subcube : 1, vertex  :3
   -0.5_WP,-0.5_WP, 0.0_WP, &     ! subcube : 1, vertex  :4
   +0.0_WP, 0.0_WP,-0.5_WP, &     ! subcube : 1, vertex  :5
   -0.5_WP, 0.0_WP,-0.5_WP, &     ! subcube : 1, vertex  :6
   +0.0_WP,-0.5_WP,-0.5_WP, &     ! subcube : 1, vertex  :7
   -0.5_WP,-0.5_WP,-0.5_WP, &     ! subcube : 1, vertex  :8
   +0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 2, vertex  :1
   +0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 2, vertex  :2
   +0.5_WP,-0.5_WP, 0.0_WP, &     ! subcube : 2, vertex  :3
   +0.0_WP,-0.5_WP, 0.0_WP, &     ! subcube : 2, vertex  :4
   +0.5_WP, 0.0_WP,-0.5_WP, &     ! subcube : 2, vertex  :5
   +0.0_WP, 0.0_WP,-0.5_WP, &     ! subcube : 2, vertex  :6
   +0.5_WP,-0.5_WP,-0.5_WP, &     ! subcube : 2, vertex  :7
   +0.0_WP,-0.5_WP,-0.5_WP, &     ! subcube : 2, vertex  :8
   +0.0_WP, 0.5_WP, 0.0_WP, &     ! subcube : 3, vertex  :1
   -0.5_WP, 0.5_WP, 0.0_WP, &     ! subcube : 3, vertex  :2
   +0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 3, vertex  :3
   -0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 3, vertex  :4
   +0.0_WP, 0.5_WP,-0.5_WP, &     ! subcube : 3, vertex  :5
   -0.5_WP, 0.5_WP,-0.5_WP, &     ! subcube : 3, vertex  :6
   +0.0_WP, 0.0_WP,-0.5_WP, &     ! subcube : 3, vertex  :7
   -0.5_WP, 0.0_WP,-0.5_WP, &     ! subcube : 3, vertex  :8
   +0.5_WP, 0.5_WP, 0.0_WP, &     ! subcube : 4, vertex  :1
   -0.0_WP, 0.5_WP, 0.0_WP, &     ! subcube : 4, vertex  :2
   +0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 4, vertex  :3
   -0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 4, vertex  :4
   +0.5_WP, 0.5_WP,-0.5_WP, &     ! subcube : 4, vertex  :5
   -0.0_WP, 0.5_WP,-0.5_WP, &     ! subcube : 4, vertex  :6
   +0.5_WP, 0.0_WP,-0.5_WP, &     ! subcube : 4, vertex  :7
   -0.0_WP, 0.0_WP,-0.5_WP, &     ! subcube : 4, vertex  :8
   +0.0_WP, 0.0_WP, 0.5_WP, &     ! subcube : 5, vertex  :1
   -0.5_WP, 0.0_WP, 0.5_WP, &     ! subcube : 5, vertex  :2
   +0.0_WP,-0.5_WP, 0.5_WP, &     ! subcube : 5, vertex  :3
   -0.5_WP,-0.5_WP, 0.5_WP, &     ! subcube : 5, vertex  :4
   +0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 5, vertex  :5
   -0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 5, vertex  :6
   +0.0_WP,-0.5_WP, 0.0_WP, &     ! subcube : 5, vertex  :7
   -0.5_WP,-0.5_WP, 0.0_WP, &     ! subcube : 5, vertex  :8
   +0.5_WP, 0.0_WP, 0.5_WP, &     ! subcube : 6, vertex  :1
   +0.0_WP, 0.0_WP, 0.5_WP, &     ! subcube : 6, vertex  :2
   +0.5_WP,-0.5_WP, 0.5_WP, &     ! subcube : 6, vertex  :3
   +0.0_WP,-0.5_WP, 0.5_WP, &     ! subcube : 6, vertex  :4
   +0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 6, vertex  :5
   +0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 6, vertex  :6
   +0.5_WP,-0.5_WP, 0.0_WP, &     ! subcube : 6, vertex  :7
   +0.0_WP,-0.5_WP, 0.0_WP, &     ! subcube : 6, vertex  :8
   +0.0_WP, 0.5_WP, 0.5_WP, &     ! subcube : 7, vertex  :1
   -0.5_WP, 0.5_WP, 0.5_WP, &     ! subcube : 7, vertex  :2
   +0.0_WP, 0.0_WP, 0.5_WP, &     ! subcube : 7, vertex  :3
   -0.5_WP, 0.0_WP, 0.5_WP, &     ! subcube : 7, vertex  :4
   +0.0_WP, 0.5_WP, 0.0_WP, &     ! subcube : 7, vertex  :5
   -0.5_WP, 0.5_WP, 0.0_WP, &     ! subcube : 7, vertex  :6
   +0.0_WP, 0.0_WP, 0.0_WP, &     ! subcube : 7, vertex  :7
   -0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 7, vertex  :8
   +0.5_WP, 0.5_WP, 0.5_WP, &     ! subcube : 8, vertex  :1
   -0.0_WP, 0.5_WP, 0.5_WP, &     ! subcube : 8, vertex  :2
   +0.5_WP, 0.0_WP, 0.5_WP, &     ! subcube : 8, vertex  :3
   -0.0_WP, 0.0_WP, 0.5_WP, &     ! subcube : 8, vertex  :4
   +0.5_WP, 0.5_WP, 0.0_WP, &     ! subcube : 8, vertex  :5
   -0.0_WP, 0.5_WP, 0.0_WP, &     ! subcube : 8, vertex  :6
   +0.5_WP, 0.0_WP, 0.0_WP, &     ! subcube : 8, vertex  :7
   -0.0_WP, 0.0_WP, 0.0_WP  &     ! subcube : 8, vertex  :8
   ],shape(cube_shift))
   
   type(tet) :: tet_gl
   type(tri) :: tri_gl
   
contains
   
   !##################################################
   !##################################################
   ! begin geomtry related routines
   
   function magnitude(v1) result(M)
      use mathtools, only: cross_product
      implicit none
      real(WP), dimension(3), intent(in) :: v1
      real(WP) :: M
      M = sqrt(dot_product(v1,v1))
   end function magnitude
   
   ! vertex is an array of vertices for the polygon(hedra)
   ! vertex(:,i) contains the x,y,z coordinates of vertex 'i'
   function tet_volume(vertex) result(V)
      use mathtools, only: cross_product
      implicit none
      real(WP),dimension(:,:), intent(in) :: vertex
      real(WP) :: V
      V = r16*abs(dot_product(vertex(:,1)-vertex(:,4),cross_product(vertex(:,2)-vertex(:,4),vertex(:,3)-vertex(:,4))))
   end function tet_volume
   
   ! vertex is an array of vertices for the triangle/tetrahedra
   ! vertex(:,i) contains the x,y,z coordinates of vertex 'i'
   function poly_centroid(vertex) result(centroid)
      implicit none
      real(WP), dimension(:,:), intent(in) :: vertex
      real(WP), dimension(3) :: centroid
      centroid = sum(vertex,DIM=2)/real(size(vertex,DIM=2),WP)
   end function poly_centroid
   
   ! v1,v2,v3 are the three vertices of the triangle (any order)
   function triangle_area(vertex) result(A)
      use mathtools, only: cross_product
      implicit none
      real(WP), dimension(:,:), intent(in) :: vertex
      real(WP) :: A
      A = 0.5_WP*magnitude(cross_product((vertex(:,3)-vertex(:,1)),(vertex(:,3)-vertex(:,2))))
   end function triangle_area
   
   subroutine split_polygon_area(vertex,area,centroid)
      implicit none
      real(WP), dimension(:,:), intent(in) :: vertex
      real(WP), intent(out) :: area
      real(WP), dimension(:), intent(out) :: centroid
      integer :: npts,i
      real(WP), dimension(3,3) :: tri_vertex
      real(WP) :: myArea
      npts = size(vertex,DIM=2)
      area = 0.0_WP
      centroid = 0.0_WP
      do i=3,npts
         tri_vertex = vertex(:,i-2:i)
         myArea = triangle_area(tri_vertex)
         area = area + myArea
         centroid = centroid+ myArea*poly_centroid(tri_vertex)
      end do
      if (area.gt.0.0_WP) centroid = centroid/area
   end subroutine split_polygon_area
   
   subroutine split_polygon_volume(vertex,volume,centroid)
      implicit none
      real(WP), dimension(:,:), intent(in) :: vertex
      real(WP), intent(out) :: volume
      real(WP), dimension(:), intent(out) :: centroid
      integer :: npts,i
      real(WP), dimension(3,4) :: tet_vertex
      real(WP) :: myVolume
      npts = size(vertex,DIM=2)
      volume = 0.0_WP
      centroid = 0.0_WP
      tet_vertex(:,1) = vertex(:,1)
      do i=4,npts
         tet_vertex(:,2:4) = vertex(:,i-2:i)
         myVolume = tet_volume(tet_vertex)
         volume = volume + myVolume
         centroid = centroid+ myVolume*poly_centroid(tet_vertex)
      end do
      if (volume.gt.0.0_WP) centroid = centroid/volume
   end subroutine split_polygon_volume
   
   ! end geometry related routines
   !##################################################
   !##################################################
   
   
   subroutine HexToTets(hex_vertex,hex_G,ntet)
      implicit none
      real(WP), dimension(3,8),intent(in) :: hex_vertex
      real(WP), dimension(8), intent(in) :: hex_G
      integer, intent(in) :: ntet
      tet_gl%vertex = hex_vertex(:,hex_to_tet(:,ntet))
      tet_gl%G = hex_G(hex_to_tet(:,ntet))
      tet_gl%Volume = tet_volume(tet_gl%vertex)
      tet_gl%T_Centroid = poly_centroid(tet_gl%vertex)
   end subroutine HexToTets
   
   subroutine computeCase(T)
      implicit none
      type(tet) :: T
      T%case = int(0.5_WP+sign(0.51_WP,T%G(1)))+2*int(0.5_WP+sign(0.51_WP,T%G(2)))+4*int(0.5_WP+sign(0.51_WP,T%G(3)))+8*int(0.5_WP+sign(0.51_WP,T%G(4)))
   end subroutine computeCase
   
   subroutine computeIntersectionPtsLinear(T)
      implicit none
      type(tet) :: T
      real(WP), dimension(3) :: p1,p2
      real(WP) :: m
      integer :: i,iedge
      ! loop over edges
      do i=1,nedges(T%case)
         !estimate intersection point
         iedge = case_to_edge(i,T%case)
         p1 = T%vertex(:,edge_to_vertex(1,iedge))
         p2 = T%vertex(:,edge_to_vertex(2,iedge))
         m = -T%G(edge_to_vertex(1,iedge))/(T%G(edge_to_vertex(2,iedge))-T%G(edge_to_vertex(1,iedge)))
         T%intp(1:3,iedge) = p1 + m *(p2-p1)
      end do
   end subroutine computeIntersectionPtsLinear
   
   subroutine RecToTris(rec_vertex,rec_G,ntri)
      implicit none
      real(WP), dimension(3,4),intent(in) :: rec_vertex
      real(WP), dimension(4), intent(in) :: rec_G
      integer, intent(in) :: ntri
      tri_gl%vertex = rec_vertex(:,rec_to_tri(:,ntri))
      tri_gl%G = rec_G(rec_to_tri(:,ntri))
      tri_gl%area = triangle_area(tri_gl%vertex)
      tri_gl%T_centroid = poly_centroid(tri_gl%vertex)
   end subroutine RecToTris
   
   subroutine triCase(T)
      implicit none
      type(tri) :: T
      ! sum over vertices
      T%case = int(0.5_WP+sign(0.51_WP,T%G(1)))+2*int(0.5_WP+sign(0.51_WP,T%G(2)))+4*int(0.5_WP+sign(0.51_WP,T%G(3)))
   end subroutine triCase
   
   subroutine triIntersectionPtsLinear(T)
      implicit none
      type(tri) :: T
      real(WP), dimension(3) :: p1,p2
      real(WP) :: m
      integer :: i,iedge
      ! loop over edges (two edges always!)
      do i=1,2
         !estimate intersection point
         iedge = tri_case_to_edge(i,T%case)
         p1 = T%vertex(:,tri_edge_to_vertex(1,iedge))
         p2 = T%vertex(:,tri_edge_to_vertex(2,iedge))
         m = -T%G(tri_edge_to_vertex(1,iedge))/(T%G(tri_edge_to_vertex(2,iedge))-T%G(tri_edge_to_vertex(1,iedge)))
         T%intp(1:3,iedge) = p1 + m *(p2-p1)
      end do
   end subroutine triIntersectionPtsLinear
   
   !##################################################
   !##################################################
   ! at this point case and interections points have been set
   ! so compute surface areas, volumes, and centroids
   
   subroutine computeIntersectionData(T)
      implicit none
      type(tet) :: T
      integer :: ne,nv
      real(WP) :: vertex(3,6)
      ne = nedges(T%case)
      nv = nverts(T%case)
      ! area
      vertex(:,1:ne) = T%intp(:,case_to_edge(1:ne,T%case))
      call split_polygon_area(vertex(:,1:ne),T%A_intersect,T%A_centroid)
      ! construct vertex array such that (0,i,i+1,i+2) for a tet
      vertex(:,1) = T%vertex(:,case_to_vertex(1,T%case))
      vertex(:,2:ne+1) = T%intp(:,case_to_edge(1:ne,T%case))
      vertex(:,ne+2:ne+nv) = T%vertex(:,case_to_vertex(2:nv,T%case))
      call split_polygon_volume(vertex(:,1:nv+ne),T%V_intersect,T%C_intersect)
   end subroutine computeIntersectionData
   
   subroutine triIntersectionData(T)
      implicit none
      type(tri) :: T
      integer :: nv
      real(WP) :: vertex(3,4)
      nv = tri_nverts(T%case)
      ! area
      vertex(:,1:nv) = T%vertex(:,tri_case_to_vertex(:,T%case))
      vertex(:,nv+1:nv+2) = T%intp(:,tri_case_to_edge(:,T%case))
      call split_polygon_area(vertex(:,1:nv+2),T%A_intersect,T%C_intersect)
   end subroutine triIntersectionData
   
   !##################################################
   !##################################################
   
   subroutine splitTetLinear(T)
      implicit none
      type(tet) :: T
      ! first determine case
      call computeCase(T)
      if (T%case == 0) then
         ! all nodes are negative
         T%V_intersect= 0.0_WP
         T%C_intersect=0.0_WP
         ! no surface area here
         T%A_intersect=0.0_WP
         T%A_centroid=0.0_WP
      else if (T%case == 15) then
         ! all nodes are positive
         T%V_intersect=T%volume
         T%C_intersect=T%T_Centroid
         ! no surface area here
         T%A_intersect=0.0_WP
         T%A_centroid=0.0_WP
      else
         ! determine intersection points
         call computeIntersectionPtsLinear(T)
         ! compute vol/area data for intersected tet
         call computeIntersectionData(T)
      end if
   end subroutine splitTetLinear
   
   subroutine splitTriLinear(T)
      implicit none
      type(tri) :: T
      ! first determine case
      call triCase(T)
      if (T%case == 0) then
         ! all nodes are negative
         T%A_intersect= 0.0_WP
         T%C_intersect=0.0_WP
      else if (T%case == 7) then
         ! all nodes are positive
         T%A_intersect=T%area
         T%C_intersect=T%T_Centroid
      else
         ! determine intersection points
         call triIntersectionPtsLinear(T)
         call triIntersectionData(T)
      end if
   end subroutine splitTriLinear

   !##################################################
   !##################################################
   
   subroutine marching_tets(hex_vertex,G,vol,area,v_cent,a_cent)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex_vertex
      real(WP), dimension(8  ), intent(in) :: G
      real(WP), intent(out) :: vol,area
      real(WP), intent(out),dimension(3) :: v_cent,a_cent
      
      !--------------------------------------------------
      integer :: i,cubeCase
      real(WP), dimension(3) :: c_dxyz
      
      ! initialize case
      cubeCase = 0
      ! initialize volume/area
      vol = 0.0_WP; v_cent = 0.0_WP; area = 0.0_WP; a_cent = 0.0_WP
      
      do i=1,8
         cubeCase = cubeCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if (cubeCase == 0) then
         ! no positive G values
         return
      end if
      
      if (cubeCase == 8) then
         ! all positive G values
         c_dxyz = hex_vertex(:,8)-hex_vertex(:,1)
         vol = product(c_dxyz)
         v_cent = poly_centroid(hex_vertex)
         return
      end if
      
      ! cube is intersected by G=0
      
      ! split cube into tets
      do i=1,6
         call HexToTets(hex_vertex,G,i)
         call splitTetLinear(tet_gl)
         ! sum up volume and areas
         vol = vol + tet_gl%V_intersect
         v_cent = v_cent + tet_gl%V_intersect*tet_gl%C_intersect
         
         area = area + tet_gl%A_intersect
         a_cent = a_cent + tet_gl%A_intersect*tet_gl%A_centroid
      end do
      
      if (vol .gt.0.0_WP) v_cent = v_cent/vol
      if (area.gt.0.0_WP) a_cent = a_cent/area
      
      return
   end subroutine marching_tets
   
   ! splits a cut cube recursively until level=0 - decrease level by one each call
   recursive subroutine cube_refine_vol(hex_vertex,vol,area,v_cent,a_cent,G_func,time,level)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex_vertex
      real(WP), intent(inout), dimension(3) :: v_cent,a_cent
      real(WP), intent(inout) :: vol,area
      real(WP), intent(in) :: time
      integer, intent(in) :: level
      interface
         function G_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: G_func
         end function G_func
      end interface
      !--------------------------------------------------
      real(WP), dimension(3,8) :: sub_cube
      integer :: i,j,cubeCase
      real(WP), dimension(8) :: G
      real(WP), dimension(3) :: c_dxyz
      
      ! initialize case
      cubeCase = 0
      ! dx,dy,dz for this cube
      c_dxyz = hex_vertex(:,8)-hex_vertex(:,1)
      
      do i=1,8
         G(i) = G_func(hex_vertex(:,i),time)
         cubeCase = cubeCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if(cubeCase == 0) then
         ! no positive G-values
         return
      end if
      
      if(cubeCase == 8) then
         v_cent = vol*v_cent +  product(c_dxyz)*poly_centroid(hex_vertex)
         vol = vol + product(c_dxyz)
         if(vol.gt.0.0_WP) v_cent = v_cent/vol
         return
      end if
      
      ! cube is intersected by G=0
      if(level > 0) then
         ! loop over sub_cubes
         do i=1,8
            ! loop over sub cube vertices
            do j=1,8
               sub_cube(:,j) = hex_vertex(:,j)+cube_shift(:,j,i)*c_dxyz
            end do
            call cube_refine_vol(sub_cube,vol,area,v_cent,a_cent,G_func,time,level-1)
         end do
         return
      end if
      
      ! split cube into tets
      do i=1,6
         call HexToTets(hex_vertex,G,i)
         call splitTetLinear(tet_gl)
         ! sum up volume and areas
         v_cent = vol*v_cent + tet_gl%V_intersect*tet_gl%C_intersect
         vol = vol + tet_gl%V_intersect
         
         a_cent = area*a_cent + tet_gl%A_intersect*tet_gl%A_centroid
         area = area + tet_gl%A_intersect
         if(vol > 0.0_WP) v_cent = v_cent/vol
         if(area> 0.0_WP) a_cent = a_cent/area
         
      end do
      return
   end subroutine cube_refine_vol
   
   ! splits a cut cube recursively until level=0 - decrease level by one each call
   recursive subroutine cube_refine_vfunc_afunc_fd(hex_vertex,vol,area,v_cent,a_cent,v_int,a_int,V_func,A_func,G_func,time,level)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex_vertex
      real(WP), intent(inout), dimension(3) :: v_cent,a_cent
      real(WP), intent(inout) :: vol,area,v_int,a_int
      real(WP), intent(in) :: time
      integer, intent(in) :: level
      interface
         function G_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: G_func
         end function G_func
         function V_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: V_func
         end function V_func
         function A_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: A_func
         end function A_func
      end interface
      !--------------------------------------------------
      real(WP), dimension(3,8) :: sub_cube
      real(WP), dimension(  8) :: G
      integer :: i,j,cubeCase
      real(WP), dimension(3) :: c_dxyz,c_cent
      real(WP) :: c_vol
      
      ! initialize case
      cubeCase = 0
      ! dx,dy,dz for this cube
      c_dxyz = hex_vertex(:,8)-hex_vertex(:,1)
      
      do i=1,8
         G(i) = G_func(hex_vertex(:,i),time)
         cubeCase = cubeCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if(cubeCase == 0) then
         ! no positive G-values
         return
      end if
      
      if(cubeCase == 8) then
         c_vol = product(c_dxyz); c_cent = poly_centroid(hex_vertex)
         v_int = v_int+c_vol*V_func(c_cent,time)
         v_cent = vol*v_cent +  c_vol*c_cent
         vol = vol + c_vol
         if(vol.gt.0.0_WP) v_cent = v_cent/vol
         return
      end if
      
      ! cube is intersected by G=0
      if(level > 0) then
         ! loop over sub_cubes
         do i=1,8
            ! loop over sub cube vertices
            do j=1,8
               sub_cube(:,j) = hex_vertex(:,j)+cube_shift(:,j,i)*c_dxyz
            end do
            call cube_refine_vfunc_afunc_fd(sub_cube,vol,area,v_cent,a_cent,v_int,a_int,V_func,A_func,G_func,time,level-1)
         end do
         return
      end if
      
      ! split cube into tets
      do i=1,6
         call HexToTets(hex_vertex,G,i)
         call splitTetLinear(tet_gl)
         ! sum up volume and areas
         v_int = v_int+tet_gl%V_intersect*V_func(tet_gl%C_intersect,time)
         v_cent = vol*v_cent + tet_gl%V_intersect*tet_gl%C_intersect
         vol = vol + tet_gl%V_intersect
         
         a_int = a_int+tet_gl%A_intersect*V_func(tet_gl%A_centroid,time)
         a_cent = area*a_cent + tet_gl%A_intersect*tet_gl%A_centroid
         area = area + tet_gl%A_intersect
         if(vol > 0.0_WP) v_cent = v_cent/vol
         if(area> 0.0_WP) a_cent = a_cent/area
         
      end do
      return
   end subroutine cube_refine_vfunc_afunc_fd
   
   ! splits a cut cube recursively until level=0 - decrease level by one each call
   recursive subroutine cube_refine_afunc_fd(hex_vertex,vol,area,v_cent,a_cent,a_int,A_func,G_func,time,level)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex_vertex
      real(WP), intent(inout), dimension(3) :: v_cent,a_cent
      real(WP), intent(inout) :: vol,area,a_int
      real(WP), intent(in) :: time
      integer, intent(in) :: level
      interface
         function G_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: G_func
         end function G_func
         function A_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: A_func
         end function A_func
      end interface
      !--------------------------------------------------
      real(WP), dimension(3,8) :: sub_cube
      real(WP), dimension(  8) :: G
      integer :: i,j,cubeCase
      real(WP), dimension(3) :: c_dxyz
      
      ! initialize case
      cubeCase = 0
      ! dx,dy,dz for this cube
      c_dxyz = hex_vertex(:,8)-hex_vertex(:,1)
      
      do i=1,8
         G(i) = G_func(hex_vertex(:,i),time)
         cubeCase = cubeCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if(cubeCase == 0) then
         ! no positive G-values
         return
      end if
      
      if(cubeCase == 8) then
         v_cent = vol*v_cent +  product(c_dxyz)*poly_centroid(hex_vertex)
         vol = vol + product(c_dxyz)
         if(vol.gt.0.0_WP) v_cent = v_cent/vol
         return
      end if
      
      ! cube is intersected by G=0
      if(level > 0) then
         ! loop over sub_cubes
         do i=1,8
            ! loop over sub cube vertices
            do j=1,8
               sub_cube(:,j) = hex_vertex(:,j)+cube_shift(:,j,i)*c_dxyz
            end do
            call cube_refine_afunc_fd(sub_cube,vol,area,v_cent,a_cent,a_int,A_func,G_func,time,level-1)
         end do
         return
      end if
      
      ! split cube into tets
      do i=1,6
         call HexToTets(hex_vertex,G,i)
         call splitTetLinear(tet_gl)
         ! sum up volume and areas
         v_cent = vol*v_cent + tet_gl%V_intersect*tet_gl%C_intersect
         vol = vol + tet_gl%V_intersect
         
         a_int = a_int+tet_gl%A_intersect*A_func(tet_gl%A_centroid,time)
         a_cent = area*a_cent + tet_gl%A_intersect*tet_gl%A_centroid
         area = area + tet_gl%A_intersect
         if(vol > 0.0_WP) v_cent = v_cent/vol
         if(area> 0.0_WP) a_cent = a_cent/area
         
      end do
      return
   end subroutine cube_refine_afunc_fd
   
   ! === 2D cutting
   
   subroutine marching_tris(rec_vertex,G,area,a_cent)
      implicit none
      real(WP), dimension(3,4), intent(in) :: rec_vertex
      real(WP), dimension(4  ), intent(in) :: G
      real(WP), intent(out) :: area
      real(WP), intent(out),dimension(3) :: a_cent
      
      !--------------------------------------------------
      integer :: i, recCase
      
      ! initialize case
      recCase = 0
      ! initialize area
      area = 0.0_WP; a_cent = 0.0_WP
      
      do i=1,4
         recCase = recCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if(recCase == 0) then
         ! no positive G values
         return
      end if
      
      if(recCase == 4) then
         ! all positive G values
         call split_polygon_area(rec_vertex,area,a_cent)
         return
      end if
      
      ! rec is intersected by G=0
      
      ! split rec into tris
      do i=1,2
         call RecToTris(rec_vertex,G,i)
         call splitTriLinear(tri_gl)
         ! sum up areas
         area = area + tri_gl%A_intersect
         a_cent = a_cent + tri_gl%A_intersect*tri_gl%C_intersect
      end do
      
      if(area> 0.0_WP) a_cent = a_cent/area
      
      return
   end subroutine marching_tris
   
   subroutine split_rec(rec_vertex,sub_rec)
      implicit none
      real(WP), dimension(3,4), intent(in) :: rec_vertex
      real(WP), dimension(3,4,4), intent(out) :: sub_rec
      real(WP), dimension(3,5) :: mid_rec
      
      
      ! midpoints to construct sub-rectangles
      !
      !    v3----m4----v4
      !    | III | IV  |
      !    |     |     |
      !    m3----m5----m2
      !    |     |     |
      !    | I   |  II |
      !    v1----m1----v2
      !
      mid_rec(:,1) = 0.5_WP*(rec_vertex(:,1)+rec_vertex(:,2))
      mid_rec(:,2) = 0.5_WP*(rec_vertex(:,2)+rec_vertex(:,4))
      mid_rec(:,3) = 0.5_WP*(rec_vertex(:,1)+rec_vertex(:,3))
      mid_rec(:,4) = 0.5_WP*(rec_vertex(:,3)+rec_vertex(:,4))
      mid_rec(:,5) = 0.5_WP*(rec_vertex(:,1)+rec_vertex(:,4))
      
      ! sub-rec I
      sub_rec(:,1,1)=rec_vertex(:,1)
      sub_rec(:,2,1)=mid_rec(:,1)
      sub_rec(:,3,1)=mid_rec(:,3)
      sub_rec(:,4,1)=mid_rec(:,5)
      ! sub-rec II
      sub_rec(:,1,2)=mid_rec(:,1)
      sub_rec(:,2,2)=rec_vertex(:,2)
      sub_rec(:,3,2)=mid_rec(:,5)
      sub_rec(:,4,2)=mid_rec(:,2)
      ! sub-rec III
      sub_rec(:,1,3)=mid_rec(:,3)
      sub_rec(:,2,3)=mid_rec(:,5)
      sub_rec(:,3,3)=rec_vertex(:,3)
      sub_rec(:,4,3)=mid_rec(:,4)
      ! sub-rec IV
      sub_rec(:,1,4)=mid_rec(:,5)
      sub_rec(:,2,4)=mid_rec(:,2)
      sub_rec(:,3,4)=mid_rec(:,4)
      sub_rec(:,4,4)=rec_vertex(:,4)
      
   end subroutine split_rec
   
   recursive subroutine rec_refine_area(rec_vertex,area,a_cent,G_func,time,level)
      implicit none
      real(WP), dimension(3,4), intent(in) :: rec_vertex
      real(WP), intent(inout) :: area
      real(WP), intent(inout), dimension(3) :: a_cent
      integer, intent(in) :: level
      real(WP), intent(in) :: time
      interface
         function G_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: G_func
         end function G_func
      end interface
      !--------------------------------------------------
      real(WP), dimension(3,4,4) :: sub_rec
      real(WP), dimension(4) :: G
      integer :: i,recCase
      real(WP), dimension(3) :: my_cent
      real(WP) :: my_area
      ! initialize case
      recCase = 0
      ! do not initialize area ...
      
      ! determine G/case
      do i=1,4
         G(i) = G_func(rec_vertex(:,i),time)
         recCase = recCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if(recCase == 0) then
         ! no positive G values
         return
      end if
      
      if(recCase == 4) then
         ! all positive G values
         call split_polygon_area(rec_vertex,my_area,my_cent)
         ! update centroid, area
         a_cent = area*a_cent + my_area*my_cent
         area = area + my_area
         if(area.gt.0.0_WP) a_cent = a_cent/area
         return
      end if
      
      if(level > 0 ) then
         call split_rec(rec_vertex,sub_rec)
         
         ! loop over sub-recs
         do i=1,4
            call rec_refine_area(sub_rec(:,:,i),area,a_cent,G_func,time,level-1)
         end do
         return
      end if
      
      ! rec is intersected by G=0 and at max refinement level
      do i=1,2
         call RecToTris(rec_vertex,G,i)
         call splitTriLinear(tri_gl)
         ! sum up areas
         a_cent = area*a_cent + tri_gl%A_intersect*tri_gl%C_intersect
         area = area + tri_gl%A_intersect
         if(area> 0.0_WP) a_cent = a_cent/area
      end do
      
      return
   end subroutine rec_refine_area
   
   recursive subroutine rec_refine_afunc_fd(rec_vertex,area,a_cent,a_int,A_func,G_func,time,level)
      implicit none
      real(WP), dimension(3,4), intent(in) :: rec_vertex
      real(WP), intent(inout) :: area,a_int
      real(WP), intent(inout), dimension(3) :: a_cent
      integer, intent(in) :: level
      real(WP), intent(in) :: time
      interface
         function G_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: G_func
         end function G_func
         function A_func(x,t)
            use precision, only: WP
            implicit none
            real(WP), dimension(3),intent(in) :: x
            real(WP), intent(in) :: t
            real(WP) :: A_func
         end function A_func
      end interface
      !--------------------------------------------------
      real(WP), dimension(3,4,4) :: sub_rec
      real(WP), dimension(4) :: G
      integer :: i,recCase
      real(WP), dimension(3) :: my_cent
      real(WP) :: my_area
      ! initialize case
      recCase = 0
      ! do not initialize area ...
      
      ! determine G/case
      do i=1,4
         G(i) = G_func(rec_vertex(:,i),time)
         recCase = recCase + int(0.5_WP+sign(0.51_WP,G(i)))
      end do
      
      if(recCase == 0) then
         ! no positive G values
         return
      end if
      
      if(recCase == 4) then
         ! all positive G values
         call split_polygon_area(rec_vertex,my_area,my_cent)
         ! update centroid, area
         a_int = a_int+my_area*A_func(my_cent,time)
         a_cent = area*a_cent + my_area*my_cent
         area = area + my_area
         if(area.gt.0.0_WP) a_cent = a_cent/area
         return
      end if
      
      if(level > 0 ) then
         call split_rec(rec_vertex,sub_rec)
         
         ! loop over sub-recs
         do i=1,4
            call rec_refine_afunc_fd(sub_rec(:,:,i),area,a_cent,a_int,A_func,G_func,time,level-1)
         end do
         return
      end if
      ! rec is intersected by G=0 and at max refinement level
      do i=1,2
         call RecToTris(rec_vertex,G,i)
         call splitTriLinear(tri_gl)
         ! sum up areas
         a_int = a_int+tri_gl%A_intersect*A_func(tri_gl%C_intersect,time)
         a_cent = area*a_cent + tri_gl%A_intersect*tri_gl%C_intersect
         area = area + tri_gl%A_intersect
         if(area> 0.0_WP) a_cent = a_cent/area
      end do
      
      return
   end subroutine rec_refine_afunc_fd
   
end module mms_geom
