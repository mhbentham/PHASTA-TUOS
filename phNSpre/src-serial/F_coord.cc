#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#endif
typedef double dArray[3];

void 
F_coord( pFace face, dArray *xyz ) {

  pPList  list = F_vertices(face,1);
  pEdge   edge;
  pVertex vtx;
  pPoint  pnt;
  int     i=0,j,n;
  void   *tmp=0;

  /* points of vertices */
  while (vtx=(pVertex)PList_next(list,&tmp)) {
    V_coord(vtx,xyz[i]);
    i++;
  }
  PList_delete(list);

  /* points on edges (higher order nodes) */
  list=F_edges(face,1,0);
  tmp=0;
  while (edge=(pEdge)PList_next(list, &tmp)) {
    if (n=E_numPoints(edge)) {
      for (j=0; j<n; j++) {
	pnt = E_point(edge,j);
	xyz[i][0] = P_x(pnt); xyz[i][1] = P_y(pnt); xyz[i][2] = P_z(pnt);
	i++;
      }
    }
    else
      break;
  }
  PList_delete(list);
}
