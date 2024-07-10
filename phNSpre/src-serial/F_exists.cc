#ifdef SIM
#include "MeshSim.h"
#include "MeshSimInternal.h"
#else
#include "AOMD.h"
#endif

pPList E_faces( pEdge edge );
pFace 
F_existsE(eType type, pEdge e1, pEdge e2, 
                     pEdge e3, pEdge e4) {
 pEdge ents[4] ;
 pFace face ;
 int  count;
 int i;
 pPList entities, faces;
 void *tmp = 0;

 ents[0] = e1;
 ents[1] = e2;
 ents[2] = e3;

 if (e4) {
   ents[3] = e4;
   count = 4;
 } else
   count = 3;

 /* if ents are vertices get the vertices of the face
 */

   faces = E_faces((pEdge) ents[0]);
   /* For each face connected to ents[0] see if the other edges are 
      connected to the face
   */
   
   for (;face = (pFace)PList_next(faces, &tmp);) {
     entities = F_edges(face, 1, (pVertex) 0);
     for( i = 1; i < count; i++)
         
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count) {
       PList_delete(faces);/* Face is connected to all entities */
       return face;
     }
   }
 PList_delete(faces);

 return (pFace)0 ;
}

pFace 
F_exists( eType type, pVertex e1, pVertex e2, 
                     pVertex e3, pVertex e4) {
 pVertex ents[4] ;
 pFace face ;
 int  count;
 int i;
 pPList entities, faces;
 void *tmp = 0;

 ents[0] = e1;
 ents[1] = e2;
 ents[2] = e3;

 if (e4) {
   ents[3] = e4;
   count = 4;
 } else
   count = 3;

 /* if ents are vertices get the vertices of the face
 */

   faces = V_faces((pVertex) ents[0]);
   /* For each face connected to ents[0] see if the other vertices are 
      connected to the face
   */
   
   for (;face =(pFace) PList_next(faces, &tmp);) {
     entities = F_vertices(face, 1);
     for( i = 1; i < count; i++)
       if (!PList_inList(entities, ents[i]))
	 break;
     PList_delete(entities);
     if (i == count){
       PList_delete(faces);/* Face is connected to all entities */
       return face;
     }
   }
 PList_delete(faces);

 return (pFace)0 ;
}
