#include "nspre_data.h"
#include "nspre_functions.h"

void 
UserWeight(pMesh mesh) {

  RIter rIter = M_regionIter(mesh);
  pRegion rgn;
  pFace face;
  pEdge edge;
  pVertex vertex;

  int NDOF;
  pPList vertices, edges, faces;
  void* temp=0;
 
  while(rgn = RIter_next(rIter)){
    /* the region modes */
    NDOF = Entity_getNumDOF((pEntity)rgn);
    
    /* the face modes */
    faces = R_faces(rgn, 1);
    temp =0;
    while(face = (pFace)PList_next(faces,&temp))
      NDOF += Entity_getNumDOF((pEntity)face);
    PList_delete(faces);

    /* the edge modes */
    edges = R_edges(rgn,1);
    temp =0;
    while(edge = (pEdge)PList_next(edges,&temp))
      NDOF += Entity_getNumDOF((pEntity)edge);
    PList_delete(edges);

    /* the vertex modes */
    vertices = R_vertices(rgn,1);
    temp=0;
    while(vertex = (pVertex)PList_next(vertices,&temp)) NDOF++;
    PList_delete(vertices);
    EN_attachDataI((pEntity)rgn,"WGHT",NDOF);
  }
}
