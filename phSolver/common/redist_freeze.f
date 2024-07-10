      module redist_freeze
      REAL*8 numredist
      
c
c tetrahedrons
c
      integer nint_tet  ! maximum number of intersection points for a tet
      integer nedges_tet   ! number of edges on tet
      parameter (nedges_tet = 6)
      integer edge_vtx_tet(nedges_tet,2)   ! defines bounding vertices for each edge
      integer int_order_tet_4(nedges_tet,
     1  nedges_tet,nedges_tet,nedges_tet,4)  ! edge order for all possible 4-pt intersection planes
c int_order_tet_4 has the following edge assignments
c    	edge 1 connects vertices 1 and 2
c    	edge 2 connects vertices 2 and 3
c    	edge 3 connects vertices 1 and 3
c    	edge 4 connects vertices 1 and 4
c    	edge 5 connects vertices 2 and 4
c    	edge 6 connects vertices 3 and 4
      data ((edge_vtx_tet(i,j),i=1,6),j=1,2)
     &                   / 1, 2, 1, 1, 2, 3,
     &                     2, 3, 3, 4, 4, 4 /
c provides reletive order of intersection edges given the intersection edges in ascending order
c i.e. if edges 1,3,5,and 6 are intersected by the interface then the relative
c       CW/CCW order is given by
c         int_order_tet_4(1,3,5,6,:) = (1,3,4,2)
c      this means that the real order is given by 1,5,6,3
c
      data (int_order_tet_4(1,3,5,6,i),i=1,4) / 1,3,4,2 /
      data (int_order_tet_4(1,2,4,6,i),i=1,4) / 1,2,4,3 /
      data (int_order_tet_4(2,3,4,5,i),i=1,4) / 2,1,4,3 /
  
      end module

