!> AMR VOF Geometry module
!> Provides cutting tables and volume computation for native VOF geometry
module amrvof_geometry
   use precision, only: WP
   implicit none
   private

   ! Expose tables and routines
   public :: tet_map, cut_side, cut_v1, cut_v2, cut_vtet
   public :: cut_ntets, cut_nvert, cut_nntet
   public :: get_plane_dist
   public :: tet_vol, tet_sign, poly_area, cut_tet_vol, cut_hex_vol
   public :: correct_flux_poly
   public :: cut_hex_polygon, hex_poly_nvert
   public :: flux_poly_moments

   ! Cutting tables from mpcomp_class_noirl
   ! tet_map: maps a hex cell (8 vertices + center) to 8 tetrahedra
   integer, dimension(4,8), parameter :: tet_map = reshape([ &
      7, 4, 3, 6, 6, 3, 2, 4, 6, 2, 1, 4, 7, 8, 4, 6, &
      6, 5, 8, 4, 6, 5, 4, 1, 5, 6, 8, 9, 6, 7, 8, 9], shape(tet_map))

   ! cut_side: which side of plane each tet vertex is on (1=below, 2=above)
   integer, dimension(6,16), parameter :: cut_side = reshape([ &
      1,-1,-1,-1,-1,-1, 2, 1, 1, 1,-1,-1, 2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, &
      2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, &
      2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, &
      2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, 2, 2, 2, 1,-1,-1, 2,-1,-1,-1,-1,-1], shape(cut_side))

   ! cut_v1, cut_v2: edge endpoints for cut vertices
   integer, dimension(4,16), parameter :: cut_v1 = reshape([ &
      -1,-1,-1,-1, 1, 1, 1,-1, 2, 2, 2,-1, 1, 2, 1, 2, 3, 3, 3,-1, 1, 3, 1, 3, &
      2, 3, 2, 3, 4, 4, 4,-1, 4, 4, 4,-1, 1, 4, 1, 4, 2, 4, 2, 4, 3, 3, 3,-1, &
      3, 4, 3, 4, 2, 2, 2,-1, 1, 1, 1,-1,-1,-1,-1,-1], shape(cut_v1))

   integer, dimension(4,16), parameter :: cut_v2 = reshape([ &
      -1,-1,-1,-1, 2, 3, 4,-1, 3, 4, 1,-1, 4, 4, 3, 3, 4, 1, 2,-1, 4, 4, 2, 2, &
      4, 4, 1, 1, 1, 2, 3,-1, 1, 2, 3,-1, 3, 3, 2, 2, 3, 3, 1, 1, 4, 1, 2,-1, &
      2, 2, 1, 1, 3, 4, 1,-1, 2, 3, 4,-1,-1,-1,-1,-1], shape(cut_v2))

   ! cut_vtet: sub-tetrahedra for each cutting configuration
   integer, dimension(4,6,16), parameter :: cut_vtet = reshape([ &
      1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 7, 6, 1, 6, 2, 3, 4, 4, 2, 5, 6, 5, 6, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1, &
      7, 5, 6, 2, 1, 3, 4, 6, 1, 5, 3, 6, 5, 7, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 8, 6, 2, 5, 7, 8, 1, 5, 1, 8, 2, 5, 6, 8, 4, 5, 8, 7, 3, 5, 8, 3, 4, &
      6, 5, 7, 3, 2, 1, 4, 6, 6, 5, 4, 2, 6, 7, 5, 2,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 6, 8, 3, 5, 8, 7, 1, 5, 8, 1, 3, 5, 8, 6, 4, 5, 7, 8, 2, 5, 8, 4, 2, &
      8, 6, 5, 3, 5, 7, 8, 2, 8, 5, 2, 3, 8, 5, 6, 4, 5, 8, 7, 1, 5, 8, 1, 4, &
      1, 2, 3, 7, 1, 2, 7, 6, 5, 7, 6, 1, 5, 6, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 6, 7, 4, 1, 2, 3, 6, 5, 1, 3, 6, 5, 7, 6, 3,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 8, 6, 4, 5, 7, 8, 1, 5, 8, 4, 1, 5, 6, 8, 3, 5, 8, 7, 2, 5, 8, 2, 3, &
      8, 5, 6, 4, 5, 8, 7, 2, 8, 2, 5, 4, 8, 6, 5, 3, 5, 7, 8, 1, 5, 8, 3, 1, &
      1, 4, 2, 7, 4, 1, 6, 7, 6, 7, 5, 4, 6, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1,-1, &
      8, 6, 5, 4, 5, 7, 8, 3, 8, 4, 5, 3, 8, 5, 6, 2, 5, 8, 7, 1, 5, 8, 1, 2, &
      3, 4, 1, 7, 7, 6, 3, 4, 7, 6, 5, 3, 7, 5, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1, &
      7, 4, 2, 3, 2, 3, 6, 7, 5, 6, 7, 2, 5, 7, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1, &
      1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1], shape(cut_vtet))

   ! Number of resulting tetrahedra for each config
   integer, dimension(16), parameter :: cut_ntets = [1,4,4,6,4,6,6,4,4,6,6,4,6,4,4,1]
   ! Number of cut vertices for each config
   integer, dimension(16), parameter :: cut_nvert = [0,3,3,4,3,4,4,3,3,4,4,3,4,3,3,0]
   ! Number of resulting tetrahedra on one side
   integer, dimension(16), parameter :: cut_nntet = [1,2,2,4,2,4,4,4,2,4,4,4,4,4,4,2]

   !> ============================================================================
   !> Hex polygon cutting tables (from IRL lookup_tables.h)
   !> For cutting a plane by a hexahedron to get the intersection polygon
   !> Vertex ordering (x,y,z):
   !>   0:(+,-,-), 1:(+,+,-), 2:(+,+,+), 3:(+,-,+)
   !>   4:(-,-,-), 5:(-,+,-), 6:(-,+,+), 7:(-,-,+)
   !> ============================================================================
   
   ! Number of polygon vertices for each of 256 cases
   ! -1 indicates invalid/ambiguous case
   integer, dimension(0:255), parameter :: hex_poly_nvert = [ &
      0,  3,  3,  4,  3, -1,  4,  5,  3,  4, -1,  5,  4,  5,  5,  4, &
      3,  4, -1,  5, -1, -1, -1, -1, -1,  5, -1,  6, -1, -1, -1,  5, &
      3, -1,  4,  5, -1, -1,  5,  6, -1, -1, -1, -1, -1, -1, -1,  5, &
      4,  5,  5,  4, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1,  4, &
      3, -1, -1, -1,  4, -1,  5, -1, -1, -1, -1, -1,  5, -1,  6,  5, &
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      4, -1,  5, -1,  5, -1,  4,  5, -1, -1, -1, -1, -1, -1,  5,  4, &
      5, -1,  6,  5, -1, -1,  5,  4, -1, -1, -1, -1, -1, -1, -1,  3, &
      3, -1, -1, -1, -1, -1, -1, -1,  4,  5, -1, -1,  5,  6, -1,  5, &
      4,  5, -1, -1, -1, -1, -1, -1,  5,  4, -1,  5, -1,  5, -1,  4, &
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      5,  6, -1,  5, -1, -1, -1, -1, -1,  5, -1,  4, -1, -1, -1,  3, &
      4, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1,  4,  5,  5,  4, &
      5, -1, -1, -1, -1, -1, -1, -1,  6,  5, -1, -1,  5,  4, -1,  3, &
      5, -1, -1, -1,  6, -1,  5, -1, -1, -1, -1, -1,  5, -1,  4,  3, &
      4,  5,  5,  4,  5, -1,  4,  3,  5,  4, -1,  3,  4,  3,  3,  0]

   ! Edge pairs for polygon vertices: hex_poly_edges(1:2, vertex, case)
   ! Each pair (v1,v2) defines an edge; interpolate to find intersection point
   ! Stored as hex_poly_v1 and hex_poly_v2 for the two endpoints
   integer, dimension(6,0:255), parameter :: hex_poly_v1 = reshape([ &
      -1,-1,-1,-1,-1,-1,  0, 3, 0,-1,-1,-1,  0, 1, 1,-1,-1,-1,  1, 3, 0, 1,-1,-1, &
       1, 2, 2,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 1, 2, 2,-1,-1,  2, 3, 0, 1, 2,-1, &
       2, 3, 3,-1,-1,-1,  0, 2, 3, 0,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 0, 1,-1, &
       1, 2, 3, 3,-1,-1,  0, 1, 2, 3, 0,-1,  0, 1, 2, 3, 3,-1,  0, 1, 2, 3,-1,-1, &
       4, 0, 7,-1,-1,-1,  0, 3, 7, 4,-1,-1, -1,-1,-1,-1,-1,-1,  1, 3, 7, 4, 1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 2, 3, 7, 4,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 7, 4, 1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  4, 1, 2, 3, 7,-1, &
       4, 5, 1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 5, 1,-1,-1,  1, 3, 0, 4, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 5, 2, 2,-1,  2, 3, 0, 4, 5, 2, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  4, 5, 2, 3, 0,-1, &
       5, 1, 0, 7,-1,-1,  0, 3, 7, 5, 1,-1,  0, 0, 7, 5, 1,-1,  1, 3, 7, 5,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 5, 2,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 7, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 2, 3, 7,-1,-1, &
       5, 6, 2,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 6, 2,-1,-1, -1,-1,-1,-1,-1,-1,  0, 1, 5, 6, 2,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 6, 3, 3,-1, -1,-1,-1,-1,-1,-1,  0, 1, 5, 6, 3, 3,  5, 6, 3, 0, 1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       4, 6, 2, 1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 2, 1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 6, 2,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 2,-1,-1,  2, 3, 0, 4, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 3, 3,-1,  4, 6, 3, 0,-1,-1, &
       6, 2, 1, 0, 7,-1, -1,-1,-1,-1,-1,-1,  0, 0, 7, 6, 2, 1,  1, 3, 7, 6, 2,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 0, 7, 6, 2,-1,  2, 3, 7, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 3, 7,-1,-1,-1, &
       6, 7, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 3,-1,-1,  0, 2, 6, 7, 0,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 2, 6, 7, 3,-1,  0, 1, 2, 6, 7, 0, -1,-1,-1,-1,-1,-1,  6, 7, 0, 1, 2,-1, &
       4, 0, 3, 6,-1,-1,  0, 3, 3, 6, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 0, 3,-1,  0, 2, 6, 4,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 6, 4, 1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 1, 2, 6, 4,-1, -1,-1,-1,-1,-1,-1,  4, 1, 2, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       5, 1, 0, 3, 6,-1,  0, 3, 3, 6, 5, 1, -1,-1,-1,-1,-1,-1,  1, 3, 3, 6, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 2, 6, 5, 1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 6, 5,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 2, 6,-1,-1,-1, &
       5, 7, 3, 2,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 7, 3, 2,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 2, 5, 7, 3,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 7, 3,-1,-1,  0, 1, 5, 7, 0,-1,  0, 1, 5, 7, 3,-1,  5, 7, 0, 1,-1,-1, &
       4, 0, 3, 2, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 2, 5, 4, 0, 3,  0, 2, 2, 5, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 4, 0, 3,-1,  0, 1, 5, 4,-1,-1, -1,-1,-1,-1,-1,-1,  4, 1, 5,-1,-1,-1, &
       4, 7, 3, 2, 1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 7, 3, 2, -1,-1,-1,-1,-1,-1,  0, 4, 7, 3, 2,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1,  0, 4, 7, 3,-1,-1,  4, 7, 0,-1,-1,-1, &
       0, 3, 2, 1,-1,-1,  0, 3, 3, 2, 1,-1,  0, 0, 3, 2, 1,-1,  1, 3, 3, 2,-1,-1, &
       1, 1, 0, 3, 2,-1, -1,-1,-1,-1,-1,-1,  0, 0, 3, 2,-1,-1,  2, 3, 3,-1,-1,-1, &
       2, 2, 1, 0, 3,-1,  0, 2, 2, 1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 2,-1,-1,-1, &
       1, 1, 0, 3,-1,-1,  0, 1, 1,-1,-1,-1,  0, 0, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1], &
      shape=[6,256])

   integer, dimension(6,0:255), parameter :: hex_poly_v2 = reshape([ &
      -1,-1,-1,-1,-1,-1,  1, 0, 4,-1,-1,-1,  1, 5, 2,-1,-1,-1,  2, 0, 4, 5,-1,-1, &
       2, 6, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 3,-1,-1,  3, 0, 4, 5, 6,-1, &
       3, 7, 0,-1,-1,-1,  1, 3, 7, 4,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 5,-1, &
       2, 6, 7, 0,-1,-1,  1, 2, 6, 7, 4,-1,  1, 5, 6, 7, 0,-1,  4, 5, 6, 7,-1,-1, &
       5, 4, 4,-1,-1,-1,  1, 0, 4, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 0, 4, 5, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  1, 3, 7, 4, 5,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 5, 5, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6, 7, 4,-1, &
       5, 6, 5,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 2,-1,-1,  2, 0, 4, 5, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 6, 3,-1,  3, 0, 4, 5, 6, 6, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 6, 6, 7, 4,-1, &
       6, 5, 4, 4,-1,-1,  1, 0, 4, 6, 5,-1,  1, 4, 4, 6, 2,-1,  2, 0, 4, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  3, 0, 4, 6, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 6, 7, 4,-1,-1, &
       6, 7, 6,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 3,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 7, 3,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 7, 0,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 7, 7, 0,  6, 7, 7, 4, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       5, 7, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 6, 2,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 7, 3,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 3,-1,-1,  3, 0, 4, 5, 7,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 7, 0,-1,  5, 7, 7, 4,-1,-1, &
       7, 6, 5, 4, 4,-1, -1,-1,-1,-1,-1,-1,  1, 4, 4, 7, 6, 2,  2, 0, 4, 7, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 4, 4, 7, 3,-1,  3, 0, 4, 7,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  7, 7, 4,-1,-1,-1, &
       7, 4, 7,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 7, 4, 0,-1,-1,  1, 3, 7, 4, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 4, 0,-1,  1, 2, 6, 7, 4, 4, -1,-1,-1,-1,-1,-1,  7, 4, 4, 5, 6,-1, &
       5, 4, 7, 7,-1,-1,  1, 0, 7, 7, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 7, 5, 4, 0,-1,  1, 3, 7, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 5, 5,-1, &
      -1,-1,-1,-1,-1,-1,  1, 2, 6, 7, 5,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6, 7,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       6, 5, 4, 7, 7,-1,  1, 0, 7, 7, 6, 5, -1,-1,-1,-1,-1,-1,  2, 0, 7, 7, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  1, 3, 7, 6, 5,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 6, 7,-1,-1,-1, &
       6, 4, 7, 6,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 6, 6, 4, 0,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 0,-1,-1,  1, 2, 6, 4, 4,-1,  1, 5, 6, 4, 0,-1,  6, 4, 4, 5,-1,-1, &
       5, 4, 7, 6, 6,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 6, 6, 5, 4, 0,  1, 3, 6, 6, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 5, 4, 0,-1,  1, 2, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6,-1,-1,-1, &
       5, 4, 7, 6, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 4, 7, 3, -1,-1,-1,-1,-1,-1,  1, 5, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 4, 0,-1, -1,-1,-1,-1,-1,-1,  1, 5, 4, 0,-1,-1,  5, 4, 4,-1,-1,-1, &
       4, 7, 6, 5,-1,-1,  1, 0, 7, 6, 5,-1,  1, 4, 7, 6, 2,-1,  2, 0, 7, 6,-1,-1, &
       2, 5, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1,  1, 4, 7, 3,-1,-1,  3, 0, 7,-1,-1,-1, &
       3, 6, 5, 4, 0,-1,  1, 3, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 6,-1,-1,-1, &
       2, 5, 4, 0,-1,-1,  1, 2, 5,-1,-1,-1,  1, 4, 0,-1,-1,-1, -1,-1,-1,-1,-1,-1], &
      shape=[6,256])

contains

   !> Compute plane distance d such that plane (n,d) cuts cell to give target VF
   !> Cell is axis-aligned box defined by corners lo and hi
   !> Normal n must be unit length
   !> Returns d in global coordinates: n.x = d
   !> Reference: Scardovelli & Zaleski, JCP 164, 228-247 (2000)
   pure function get_plane_dist(normal, lo, hi, VF) result(d)
      implicit none
      real(WP), dimension(3), intent(in) :: normal  !< Unit normal
      real(WP), dimension(3), intent(in) :: lo      !< Lower cell corner
      real(WP), dimension(3), intent(in) :: hi      !< Upper cell corner
      real(WP), intent(in) :: VF                    !< Target volume fraction
      real(WP) :: d
      real(WP) :: m1, m2, m3, mm1, mm2, mm3, tmp
      real(WP) :: norm, factor, alpha, VOFo
      real(WP) :: m12, V1, V2, V3
      real(WP) :: a0, a1, a2, p0, q0, theta
      real(WP), parameter :: eps = 1.0e-15_WP
      real(WP), dimension(3) :: cellsize, cellctr
      
      ! Compute cell size and centroid from corners
      cellsize = hi - lo
      cellctr  = 0.5_WP * (lo + hi)
      
      ! Form scaled normal components
      mm1 = normal(1) * cellsize(1)
      mm2 = normal(2) * cellsize(2)
      mm3 = normal(3) * cellsize(3)
      
      ! Normalize to sum of absolute values = 1
      norm = abs(mm1) + abs(mm2) + abs(mm3)
      if (norm.lt.eps) then
         d = sign(1.0e10_WP, VF - 0.5_WP)
         return
      end if
      mm1 = mm1 / norm
      mm2 = mm2 / norm
      mm3 = mm3 / norm

      
      ! Track offset for negative components
      factor = 0.0_WP
      if (mm1.lt.0.0_WP) factor = factor + mm1
      if (mm2.lt.0.0_WP) factor = factor + mm2
      if (mm3.lt.0.0_WP) factor = factor + mm3
      
      ! Handle VF > 0.5 by symmetry
      VOFo = VF
      if (VF.gt.0.5_WP) VOFo = 1.0_WP - VF
      
      ! === Core S&Z algorithm ===
      ! Take absolute value and sort so that m1 <= m2 <= m3
      m1 = abs(mm1); m2 = abs(mm2); m3 = abs(mm3)
      if (m2.lt.m1) then; tmp = m2; m2 = m1; m1 = tmp; end if
      if (m3.lt.m2) then; tmp = m3; m3 = m2; m2 = tmp; end if
      if (m2.lt.m1) then; tmp = m2; m2 = m1; m1 = tmp; end if
      
      ! Form volume thresholds
      m12 = m1 + m2
      V1 = m1**2 / max(6.0_WP * m2 * m3, eps)
      V2 = V1 + 0.5_WP * (m2 - m1) / max(m3, eps)
      if (m12.le.m3) then
         V3 = 0.5_WP * m12 / max(m3, eps)
      else
         V3 = (m3**2 * (3.0_WP*m12 - m3) + m1**2 * (m1 - 3.0_WP*m3) + &
              m2**2 * (m2 - 3.0_WP*m3)) / max(6.0_WP * m1 * m2 * m3, eps)
      end if
      
      ! Calculate alpha based on regime
      if (VOFo.ge.0.0_WP.and.VOFo.lt.V1) then
         alpha = (6.0_WP * m1 * m2 * m3 * VOFo)**(1.0_WP/3.0_WP)
      else if (VOFo.lt.V2) then
         alpha = 0.5_WP * (m1 + sqrt(m1**2 + 8.0_WP * m2 * m3 * (VOFo - V1)))
      else if (VOFo.lt.V3) then
         a0 = -(m1**3 + m2**3 - 6.0_WP * m1 * m2 * m3 * VOFo)
         a1 = 3.0_WP * (m1**2 + m2**2)
         a2 = -3.0_WP * m12
         p0 = -(a1 / 3.0_WP - a2**2 / 9.0_WP)
         q0 = (a1 * a2 - 3.0_WP * a0) / 6.0_WP - a2**3 / 27.0_WP
         theta = acos(max(-1.0_WP, min(1.0_WP, q0 / sqrt(p0**3)))) / 3.0_WP
         alpha = sqrt(p0) * (sqrt(3.0_WP) * sin(theta) - cos(theta)) - a2 / 3.0_WP
      else
         if (m12.le.m3) then
            alpha = m3 * VOFo + 0.5_WP * m12
         else
            a0 = -0.5_WP * (m1**3 + m2**3 + m3**3 - 6.0_WP * m1 * m2 * m3 * VOFo)
            a1 = 1.5_WP * (m1**2 + m2**2 + m3**2)
            a2 = -1.5_WP
            p0 = -(a1 / 3.0_WP - a2**2 / 9.0_WP)
            q0 = (a1 * a2 - 3.0_WP * a0) / 6.0_WP - a2**3 / 27.0_WP
            theta = acos(max(-1.0_WP, min(1.0_WP, q0 / sqrt(p0**3)))) / 3.0_WP
            alpha = sqrt(p0) * (sqrt(3.0_WP) * sin(theta) - cos(theta)) - a2 / 3.0_WP
         end if
      end if
      ! === End core S&Z algorithm ===
      
      ! Undo VF > 0.5 symmetry
      if (VF.gt.0.5_WP) alpha = 1.0_WP - alpha
      
      ! Adjust for negative normal components
      alpha = alpha + factor
      
      ! Transform from cell-corner origin to cell-centroid origin
      alpha = alpha - 0.5_WP * (mm1 + mm2 + mm3)
      
      ! Scale back to physical coordinates and add centroid offset
      d = alpha * norm + dot_product(normal, cellctr)
      
   end function get_plane_dist

   !> Compute volume of tetrahedron given 4 vertices
   !> v(:,1:4) are the vertex coordinates
   pure function tet_vol(v) result(vol)
      implicit none
      real(WP), dimension(3,4), intent(in) :: v
      real(WP) :: vol
      real(WP), dimension(3) :: a,b,c
      a=v(:,1)-v(:,4)
      b=v(:,2)-v(:,4)
      c=v(:,3)-v(:,4)
      vol=(-a(1)*(b(2)*c(3)-c(2)*b(3))&
      &    +a(2)*(b(1)*c(3)-c(1)*b(3))&
      &    -a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
   end function tet_vol

   !> Function that calculates the sign of a tet
   function tet_sign(vert) result(s)
      implicit none
      real(WP) :: s
      real(WP), dimension(3,4), intent(in) :: vert
      real(WP), dimension(3) :: a,b,c
      a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
      s=sign(1.0_WP,-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2))))
   end function tet_sign

   !> Compute area of a convex polygon with nv vertices via triangle fan
   !> verts(:,1:nv) are the polygon vertices in cyclic order
   pure function poly_area(nv,verts) result(area)
      implicit none
      integer,  intent(in) :: nv
      real(WP), dimension(3,nv), intent(in) :: verts
      real(WP) :: area
      real(WP), dimension(3) :: a,b,c
      integer :: n
      area=0.0_WP
      do n=2,nv-1
         a=verts(:,n  )-verts(:,1)
         b=verts(:,n+1)-verts(:,1)
         c=[a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
         area=area+0.5_WP*sqrt(c(1)**2+c(2)**2+c(3)**2)
      end do
   end function poly_area
   
   !> Cut a tetrahedron by a plane and return liquid/gas volumes and barycenters
   !> Input:  v(:,1:4) = 4 tet vertices, plane(1:4) = [nx,ny,nz,d] where n.x=d defines plane
   !> Output: vol_liq, vol_gas = phase volumes
   !>         bary_liq(3), bary_gas(3) = volume-weighted barycenters
   pure subroutine cut_tet_vol(v, plane, vol_liq, vol_gas, bary_liq, bary_gas)
      implicit none
      real(WP), dimension(3,4), intent(in)  :: v
      real(WP), dimension(4),   intent(in)  :: plane
      real(WP),                 intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3),   intent(out) :: bary_liq, bary_gas
      real(WP), dimension(4)   :: d
      real(WP), dimension(3,8) :: vert
      real(WP), dimension(3)   :: a, b, c, bary
      real(WP) :: denom, mu, my_vol
      integer  :: icase, n1, v1, v2
      
      ! Compute signed distance of each vertex to plane
      d(1) = plane(1)*v(1,1) + plane(2)*v(2,1) + plane(3)*v(3,1) - plane(4)
      d(2) = plane(1)*v(1,2) + plane(2)*v(2,2) + plane(3)*v(3,2) - plane(4)
      d(3) = plane(1)*v(1,3) + plane(2)*v(2,3) + plane(3)*v(3,3) - plane(4)
      d(4) = plane(1)*v(1,4) + plane(2)*v(2,4) + plane(3)*v(3,4) - plane(4)
      
      ! Determine cutting case (1-16 based on which vertices are above plane)
      icase = 1 + int(0.5_WP+sign(0.5_WP,d(1))) + 2*int(0.5_WP+sign(0.5_WP,d(2))) &
                + 4*int(0.5_WP+sign(0.5_WP,d(3))) + 8*int(0.5_WP+sign(0.5_WP,d(4)))
      
      ! Copy original vertices
      vert(:,1:4) = v(:,1:4)
      
      ! Create interpolated vertices on cut plane
      do n1 = 1, cut_nvert(icase)
         v1 = cut_v1(n1,icase); v2 = cut_v2(n1,icase)
         denom = d(v2) - d(v1)
         if (abs(denom).ge.tiny(1.0_WP)) then
            mu = max(0.0_WP, min(1.0_WP, -d(v1)/denom))
         else
            mu = 0.0_WP
         end if
         vert(:,4+n1) = (1.0_WP-mu)*vert(:,v1) + mu*vert(:,v2)
      end do
      
      ! Initialize outputs
      vol_liq = 0.0_WP; vol_gas = 0.0_WP
      bary_liq = 0.0_WP; bary_gas = 0.0_WP
      
      ! Gas tets: from 1 to cut_nntet-1
      do n1 = 1, cut_nntet(icase)-1
         a = vert(:,cut_vtet(1,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         b = vert(:,cut_vtet(2,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         c = vert(:,cut_vtet(3,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         bary = 0.25_WP * (vert(:,cut_vtet(1,n1,icase)) + vert(:,cut_vtet(2,n1,icase)) &
                         + vert(:,cut_vtet(3,n1,icase)) + vert(:,cut_vtet(4,n1,icase)))
         vol_gas = vol_gas + my_vol
         bary_gas = bary_gas + my_vol * bary
      end do
      
      ! Liquid tets: from cut_ntets down to cut_nntet
      do n1 = cut_ntets(icase), cut_nntet(icase), -1
         a = vert(:,cut_vtet(1,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         b = vert(:,cut_vtet(2,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         c = vert(:,cut_vtet(3,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         bary = 0.25_WP * (vert(:,cut_vtet(1,n1,icase)) + vert(:,cut_vtet(2,n1,icase)) &
                         + vert(:,cut_vtet(3,n1,icase)) + vert(:,cut_vtet(4,n1,icase)))
         vol_liq = vol_liq + my_vol
         bary_liq = bary_liq + my_vol * bary
      end do
      
      ! Normalize barycenters
      if (vol_liq.gt.tiny(1.0_WP)) bary_liq = bary_liq / vol_liq
      if (vol_gas.gt.tiny(1.0_WP)) bary_gas = bary_gas / vol_gas
      
   end subroutine cut_tet_vol

   !> Adjust flux polyhedron to enforce target volume
   !> Uses IRL's direction-independent approach: moves vertex 9 along back-face normal
   subroutine correct_flux_poly(poly,target_volume)
      implicit none
      real(WP), dimension(3,9), intent(inout) :: poly
      real(WP), intent(in) :: target_volume
      real(WP) :: starting_volume,needed_change,mag,adjustment
      real(WP), dimension(3) :: cross_sum,dir
      real(WP), dimension(3) :: e1,e2,c1,c2
      integer :: n
      ! Compute starting volume
      starting_volume=0.0_WP
      do n=1,8; starting_volume=starting_volume+tet_vol([poly(:,tet_map(1,n)),poly(:,tet_map(2,n)),poly(:,tet_map(3,n)),poly(:,tet_map(4,n))]); end do
      needed_change=target_volume-starting_volume
      ! Compute volume gradient for tets 7 and 8 (the only tets using vertex 9)
      ! For tet (A,B,C,D), gradient w.r.t. D is -(1/6)*((B-A)×(C-A))
      ! Tet 7 = (5,6,8,9): gradient = -((v6-v5)×(v8-v5)). Take edges from vertex 5
      e1=poly(:,6)-poly(:,5); e2=poly(:,8)-poly(:,5)
      c1=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
      ! Tet 8 = (6,7,8,9): gradient = -((v7-v6)×(v8-v6)). Take edges from vertex 6
      e1=poly(:,7)-poly(:,6); e2=poly(:,8)-poly(:,6)
      c2=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
      ! Total gradient (negative sign absorbed into adjustment formula)
      cross_sum=c1+c2
      ! Compute adjustment along normal direction
      mag=sqrt(cross_sum(1)**2+cross_sum(2)**2+cross_sum(3)**2)
      adjustment=6.0_WP*needed_change/max(mag,tiny(1.0_WP))
      dir=cross_sum/max(mag,tiny(1.0_WP))
      ! Move vertex 9
      poly(:,9)=poly(:,9)+adjustment*dir
   end subroutine correct_flux_poly

   !> Compute total signed volume and volume-weighted barycenter of a
   !> corrected 9-vertex flux polyhedron using tet_map(4,8) decomposition
   !> No PLIC cutting, no grid-plane recursion.
   pure subroutine flux_poly_moments(poly,vol,bary)
      implicit none
      real(WP), dimension(3,9), intent(in)  :: poly
      real(WP),                 intent(out) :: vol
      real(WP), dimension(3),   intent(out) :: bary
      real(WP), dimension(3) :: a,b,c,centroid
      real(WP) :: svol
      integer :: n
      vol =0.0_WP
      bary=0.0_WP
      do n=1,8
         a=poly(:,tet_map(1,n))-poly(:,tet_map(4,n))
         b=poly(:,tet_map(2,n))-poly(:,tet_map(4,n))
         c=poly(:,tet_map(3,n))-poly(:,tet_map(4,n))
         svol=(-a(1)*(b(2)*c(3)-c(2)*b(3))+a(2)*(b(1)*c(3)-c(1)*b(3))-a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
         centroid=0.25_WP*(poly(:,tet_map(1,n))+poly(:,tet_map(2,n))+poly(:,tet_map(3,n))+poly(:,tet_map(4,n)))
         vol=vol+svol
         bary=bary+svol*centroid
      end do
   end subroutine flux_poly_moments

   !> Cut a hex cell by a plane and compute liquid/gas volumes and barycenters
   !> hex(:,1:8) = 8 vertices of hex cell (standard ordering)
   !> plane(1:3) = normal, plane(4) = distance (n·x = d)
   !> Returns vol_liq/vol_gas and bary_liq/bary_gas
   subroutine cut_hex_vol(hex, plane, vol_liq, vol_gas, bary_liq, bary_gas)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex
      real(WP), dimension(4), intent(in) :: plane
      real(WP), intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3), intent(out) :: bary_liq, bary_gas
      ! Local
      real(WP), dimension(3,4) :: tet
      real(WP), dimension(3) :: center
      real(WP) :: vl, vg
      real(WP), dimension(3) :: bl, bg
      integer :: ntet, v
      ! Hex->tet decomposition: 6 tets using center point
      ! hex_tet_map: indices into hex vertices for 6 tets (each shares center)
      integer, dimension(3,6), parameter :: hex_tet_map = reshape([ &
         1, 2, 3,  1, 3, 4,  1, 2, 6,  1, 6, 5, &
         5, 6, 7,  5, 7, 8], shape(hex_tet_map))
      integer, dimension(3,6), parameter :: hex_tet_top = reshape([ &
         5, 6, 7,  5, 7, 8,  4, 3, 7,  4, 7, 8, &
         1, 2, 3,  1, 3, 4], shape(hex_tet_top))
      ! Actually use standard 5-tet decomposition of hex
      integer, dimension(4,5), parameter :: tet5 = reshape([ &
         1, 2, 4, 5,  2, 3, 4, 7,  2, 5, 6, 7, &
         4, 5, 7, 8,  2, 4, 5, 7], shape(tet5))
      
      vol_liq = 0.0_WP; vol_gas = 0.0_WP
      bary_liq = 0.0_WP; bary_gas = 0.0_WP
      
      ! Decompose hex into 5 tetrahedra
      do ntet = 1, 5
         ! Build tetrahedron
         do v = 1, 4
            tet(:,v) = hex(:,tet5(v,ntet))
         end do
         ! Cut tet by plane
         call cut_tet_vol(tet, plane, vl, vg, bl, bg)
         ! Accumulate volumes and volume-weighted barycenters
         bary_liq = bary_liq + vl * bl
         bary_gas = bary_gas + vg * bg
         vol_liq = vol_liq + vl
         vol_gas = vol_gas + vg
      end do
      
      ! Normalize barycenters
      if (vol_liq.gt.tiny(1.0_WP)) bary_liq = bary_liq / vol_liq
      if (vol_gas.gt.tiny(1.0_WP)) bary_gas = bary_gas / vol_gas
      
   end subroutine cut_hex_vol

   !> =========================================================================
   !> Extract PLIC polygon from hex cell using IRL's algorithm
   !> Given 8 hex vertices and a plane (normal, distance), compute the 
   !> intersection polygon vertices.
   !> 
   !> Hex vertex ordering (IRL convention):
   !>   0:(+,-,-), 1:(+,+,-), 2:(+,+,+), 3:(+,-,+)
   !>   4:(-,-,-), 5:(-,+,-), 6:(-,+,+), 7:(-,-,+)
   !>
   !> For axis-aligned cell [xlo,ylo,zlo] to [xhi,yhi,zhi]:
   !>   hex(:,1) = [xhi, ylo, zlo]  ! vertex 0
   !>   hex(:,2) = [xhi, yhi, zlo]  ! vertex 1
   !>   hex(:,3) = [xhi, yhi, zhi]  ! vertex 2
   !>   hex(:,4) = [xhi, ylo, zhi]  ! vertex 3
   !>   hex(:,5) = [xlo, ylo, zlo]  ! vertex 4
   !>   hex(:,6) = [xlo, yhi, zlo]  ! vertex 5
   !>   hex(:,7) = [xlo, yhi, zhi]  ! vertex 6
   !>   hex(:,8) = [xlo, ylo, zhi]  ! vertex 7
   !>
   !> Returns:
   !>   nvert: number of polygon vertices (0-6)
   !>   poly(:,1:nvert): polygon vertices in cyclic order
   !> =========================================================================
   pure subroutine cut_hex_polygon(hex, plane, nvert, poly)
      implicit none
      real(WP), dimension(3,8), intent(in)  :: hex      !< 8 hex vertices
      real(WP), dimension(4),   intent(in)  :: plane    !< (nx, ny, nz, d)
      integer,                  intent(out) :: nvert    !< Number of polygon vertices
      real(WP), dimension(3,6), intent(out) :: poly     !< Polygon vertices (up to 6)
      
      ! Local variables
      real(WP), dimension(8) :: dist       ! Signed distance to plane for each vertex
      integer :: icase                     ! Cutting case (0-255)
      integer :: n, v1, v2
      real(WP) :: mu, denom
      
      ! Initialize output
      nvert = 0
      poly = 0.0_WP
      
      ! Compute signed distance from plane for each vertex
      ! Plane equation: n.x - d = 0, so distance = n.x - d
      do n = 1, 8
         dist(n) = plane(1)*hex(1,n) + plane(2)*hex(2,n) + plane(3)*hex(3,n) - plane(4)
      end do
      
      ! Determine cutting case from sign pattern (0-indexed vertices)
      ! Case bit i is set if vertex i is above plane (dist > 0)
      icase = 0
      if (dist(1).gt.0.0_WP) icase = icase + 1    ! vertex 0
      if (dist(2).gt.0.0_WP) icase = icase + 2    ! vertex 1
      if (dist(3).gt.0.0_WP) icase = icase + 4    ! vertex 2
      if (dist(4).gt.0.0_WP) icase = icase + 8    ! vertex 3
      if (dist(5).gt.0.0_WP) icase = icase + 16   ! vertex 4
      if (dist(6).gt.0.0_WP) icase = icase + 32   ! vertex 5
      if (dist(7).gt.0.0_WP) icase = icase + 64   ! vertex 6
      if (dist(8).gt.0.0_WP) icase = icase + 128  ! vertex 7
      
      ! Get number of polygon vertices for this case
      nvert = hex_poly_nvert(icase)
      
      ! Handle no intersection or ambiguous cases
      if (nvert.le.0) then
         nvert = 0
         return
      end if
      
      ! Compute intersection points along cut edges
      do n = 1, nvert
         ! Get edge endpoints (0-indexed in tables, convert to 1-indexed)
         v1 = hex_poly_v1(n, icase) + 1
         v2 = hex_poly_v2(n, icase) + 1
         
         ! Linear interpolation: find point where dist = 0
         ! mu = -dist(v1) / (dist(v2) - dist(v1))
         denom = dist(v2) - dist(v1)
         if (abs(denom).ge.tiny(1.0_WP)) then
            mu = max(0.0_WP, min(1.0_WP, -dist(v1)/denom))
         else
            mu = 0.5_WP
         end if
         
         ! Interpolated point
         poly(:,n) = (1.0_WP - mu)*hex(:,v1) + mu*hex(:,v2)
      end do
      
   end subroutine cut_hex_polygon

end module amrvof_geometry