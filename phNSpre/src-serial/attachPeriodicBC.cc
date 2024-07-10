/*

collect the periodic boundary condition array by finding
the degrees of freedom associated with periodic mesh entities

*/
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#endif
#include "nspre_functions.h"

/**********************************************************************/
/* compute the periodic boundary condition array. Periodic boundary   */
/* conditions which reside on other processors are treated as a       */
/* communication, and thus are not included in iper.                  */
/**********************************************************************/
void 
attachPeriodicBC( pMesh mesh, int *iper ) {

    pVertex vertex,mvert;
    pEdge edge,medge;
    pFace face,mface;
    int i,nem,nfm,sdof,mdof;

  
    /* vertex modes */
    VIter vIter = M_vertexIter(mesh);
    while (vertex = VIter_next(vIter))  {
        if (mvert = (pVertex)EN_dataP((pEntity) vertex,"PerM")){
            sdof = Entity_getDOFnumber((pEntity) vertex,0);
            mdof = Entity_getDOFnumber((pEntity) mvert,0);
            iper[sdof] = mdof;
        }
    }

    VIter_delete(vIter);
   
    /* edge modes */
    EIter eIter = M_edgeIter(mesh);
    while (edge = EIter_next(eIter))  {
        if ( isActive( (pEntity) edge ) ){
            if( medge = (pEdge)EN_dataP((pEntity) edge,"PerM") ){
                nem = Entity_getNumDOF((pEntity) edge);
                for (i=0; i < nem; i++){
                    sdof = Entity_getDOFnumber((pEntity) edge, i);
                    mdof = Entity_getDOFnumber((pEntity) medge, i);
                    iper[sdof] = mdof;
                }
            }
        }
    }
    EIter_delete(eIter);

    /* face modes */
    FIter fIter = M_faceIter(mesh);
    while (face = FIter_next(fIter))  {
        if ( isActive( (pEntity) face ) ) {
            if ( mface = (pFace)EN_dataP((pEntity) face,"PerM") ) {
                nfm = Entity_getNumDOF((pEntity) face);
                for (i=0; i < nfm; i++){
                    sdof = Entity_getDOFnumber((pEntity) face,i);
                    mdof = Entity_getDOFnumber((pEntity) mface,i);
                    iper[sdof] = mdof;
                }
            }
        }
    }
    FIter_delete(fIter);
}
