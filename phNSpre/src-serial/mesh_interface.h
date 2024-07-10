/*missing functions when using SCOREC libraries
 Min Zhou spring 2008
*/

#ifndef SIM
#ifndef MESH_INTERFACE_
#define MESH_INTERFACE_
#include "AOMD.h"
#include "MeshAdjTools.h"
#include <iostream>
inline pVertex F_oppositeVertex(pFace face, pEdge edge){
    return F_edOpVt(face,edge);
}
inline void E_info(pEdge edge){
    std::cout<<"E_info() not implemented yet\n";
}
inline void F_info(pFace face) { 
    std::cout<<"F_info() not implemented yet\n";
}
inline pMesh M_new(int type, pGModel model){
    return  MS_newMesh(model);
}
inline void M_release(pMesh mesh){
    M_delete(mesh);
}
inline void GM_release(pGModel model){
    GM_delete(model);
}
pEdge R_gtOppEdg(pRegion region, pEdge edge);
pFace R_vtOpFc(pRegion region,pVertex vertex);
int F_conToFace(pFace face1, pFace face2) ;
int EN_inClosure(pEntity, pEntity);
int C_raneql(double real1,double real2,double tolran);

#endif
#endif
