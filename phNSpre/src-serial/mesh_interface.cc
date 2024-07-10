/* missing functions when using SCOREC libraries
  Min Zhou spring 2008
*/

#ifndef SIM
#include "mesh_interface.h"
#include "MSops.h"

pEdge R_gtOppEdg(pRegion currgn,pEdge edge) {
  pEdge   opp_edge;
  pVertex v, v0,v1,oppv[2];
  pPList  vlist;
  int     i,k=0;
  void   *tmp=0;

  /* get the two vertices of the edge */
  v0=E_vertex(edge,0);
  v1=E_vertex(edge,1);

  /* get the verticies of this region */
  vlist = R_vertices(currgn,1);
  if (PList_size(vlist)!=4) {
      std::cout<<"Wrong number of vertices\n";
    return(0);
  }

  /* find the opposite edge in currgn */
  while (v=(pVertex)PList_next(vlist,&tmp)) {
    if (v!=v0 && v!=v1) {
      oppv[k]=v;
      k++;
    }
  }
  PList_delete(vlist);
  opp_edge=E_exists(oppv[0],oppv[1]);
  if (!opp_edge) {
      std::cout<<"Edge not found","R_gtOppEdg\n";
    return(0);
  }
  return opp_edge;
}

pFace R_vtOpFc(pRegion rgn,pVertex vtx)
{
  pPList vlist ;
  pFace  face ;
  int    i; 

  if ( R_numFaces(rgn)>4 )
    return((pFace)0) ;
  for ( i=0 ; i<R_numFaces(rgn) ; i++ ) {
    face = R_face(rgn,i) ;
    vlist = F_vertices(face,1) ;
    if ( !PList_inList(vlist,vtx) ) {
      PList_delete(vlist) ;
      return(face) ;
    }
    PList_delete(vlist) ;
  }
  return((pFace)0) ;
}

int F_conToFace(pFace face1, pFace face2) {
  pPList fverts, fedges;
  pEdge edge;
  pVertex vtx;
  void *temp;

  temp = 0;
  fverts = F_vertices(face1,1);
  while (vtx = (pVertex)PList_next(fverts,&temp))
    if (F_inClosure(face2,vtx)) {
      PList_delete(fverts);
      return 1;
    }
  PList_delete(fverts);

  temp = 0;
  fedges = F_edges(face1,0,0);
  while (edge = (pEdge)PList_next(fedges,&temp))
    if (F_inClosure(face2,edge)) {
      PList_delete(fedges);
      return 1;
    }
  PList_delete(fedges);

  return 0;
}

int EN_inClosure(pEntity e1, pEntity e2) {

  switch(EN_type(e1)) {
  case Tvertex:
    return (EN_type(e2) != Tvertex) ? 0 : (e1 == e2);
  case Tedge:
    return E_inClosure((pEdge)e1,e2);
  case Tface:
    return F_inClosure((pFace)e1,e2);
  case Tregion:
    return R_inClosure((pRegion)e1,e2);
  }

  return 0;
}

int C_raneql(double real1,double real2,double tolran)
{
  return ( abs(real1-real2) <= tolran ) ;
}

#endif
