#ifndef _NSPRE_INCLUDE_
#define _NSPRE_INCLUDE_

#include "nspre_data.h"
#include "EnsaParameters.h"
#ifndef SIM
#include "MeshEntTools.h"
#include "mesh_interface.h"
#endif
#define MAXNT 100
#define MAX(a,b) ((a) < (b) ? (b) : (a) )
#define MIN(a,b) ((a) > (b) ? (b) : (a) )

#ifdef SIM
typedef class MVertex * pVertex;
typedef class MEdge * pEdge;
typedef class MFace * pFace;
typedef class MRegion * pRegion;
#endif

int topology( pRegion region );
int topology2D( pFace face );
void fixRegions( pEntity ent );
void UserWeight(pMesh mesh);
void generate_keyphrase(char* target, char* prefix, blockKey* tpblock);
void procArgs(int argc, char *argv[], char fname[], char mname[], char gname[] );

void V_reordering(pMesh mesh, globalInfo* info);
void R_reordering(pMesh mesh);

int Entity_getNumDOF(pEntity entity);
void Entity_setDOFnumber(pEntity entity, int gdof);
int Entity_getDOFnumber(pEntity entity, int ith);
void Entity_addOffset( pEntity entity, int offset );
void readBC(pGModel model,pMesh mesh);
void writeEnsaFiles( pGModel model,pMesh mesh, void* bdry);
void getX(pMesh mesh, double *x );
void getConnectivity(pMesh mesh, vector< pFace >* bdry, 
                     int** ien, int** ien_sms,  int** ienb, int** ienb_sms);
void getConnectivity2D(pMesh mesh, vector< pEdge >* bdry, 
                     int** ien, int** ien_sms,  int** ienb, int** ienb_sms);
void 
getFaceConnectivity(pMesh mesh, 
                    vector< pFace > bdry, 
                    int** ief,
                    int** iefb);

// returns number of faces of an element of a particular block
// (defined  in localInfo.cc)
int 
NEFPerBlock(blockKey key);

void R_entities(pRegion region, pVertex *vrts, pEdge *edgs, pFace *fcs, int nen);
void R_entitiesBdry(pRegion region, pFace face, pVertex *vrts, pEdge *edgs, pFace *fcs,int nen);
void attachPeriodicBC(pMesh mesh, int *iper);
void restart(double* q, int nshg, int nvr, int* lstep, char filename[] );
void F_coord( pFace face, dArray *xyz);
pFace F_exists(eType type, pVertex e1, pVertex e2, pVertex e3, pVertex e4);
pFace F_existsE(eType type, pEdge e1, pEdge e2, pEdge e3, pEdge e4);
void printPeriodicBC(pGModel model, pMesh mesh);
int getNaturalBC(pGEntity gface, pEntity face, double *BCB);
int attachEssentialBC( pGModel model, pMesh mesh, double *qtmp, int *nBC,
        		       int *iBC, double **BC );
void attachNaturalBC( pGModel model,pMesh mesh, int ***iBCB, double ***BCB, int numNBC);
void attachInitialCondition(pGModel model, pMesh mesh, double *qt, double **q);
int getBC(pGEntity gent, pEntity ent, double *BC, int i);
void setPeriodic(pGModel model,pMesh mesh);
pVertex getMasterV(VIter mvIter, pVertex vert, double theta, int trans );
pEdge getMasterE(EIter meIter, pEdge edg, double theta, int trans);
pFace getMasterF(FIter mfIter, pFace fac, double theta);
void gendual( pGModel model, pMesh mesh );
void setup_and_refine(pGModel model, pMesh mesh, globalInfo* info);
void lin2quad(pMesh mesh, pMesh lmesh, double* qtmpl, double* qtmp, int nshg);
void genblock( pMesh mesh, pGModel model, vector<pFace>& bdry);
void genblock2D( pMesh mesh, pGModel model, vector<pEdge>& bdry);
void face_coordinate_extraction(pMesh mesh, pGModel model);
void prefine_on_error( pMesh mesh, 
                       globalInfo* info, 
                       map< int, vector< int > > idshpmap,
                       int nshg_old,
                       char filename[] );

void
AllocateEnsaArrays( EnsaArrays* e, 
                    EnsaParameters *par) ;
void
DeAllocateEnsaArrays( EnsaArrays* e ) ;

void 
writeEnsaArrays( EnsaArrays* e ,
                 EnsaParameters *par, 
                 double zscale[3] );

int FakeMode( globalInfo* info ) ;
int numVertsR( pRegion region );
int numVertsF( pFace face );

bool isActive( pEntity ent );
void eqn_plane(pGModel, pMesh mesh, double *coeffs, pGFace gface);
int isCut(pPList verts, double *coeffs);
void GF_inflow_normal(pGFace gface, pGModel model, double *xyz);
#endif
