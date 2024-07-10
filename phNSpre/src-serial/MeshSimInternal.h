#ifdef SIM
#ifndef H_MeshSimInternal
#define H_MeshSimInternal

#include "MeshSim.h"

#ifndef SIMINTERNAL


/* Geometry type codes */ 
/* enum GEOM_type {*/
/*     GEOM_Unknown,*/
/*     GEOM_Point,*/
/*     GEOM_Line,*/
/*     GEOM_Circle,*/
/*     GEOM_Ellipse,*/
/*     GEOM_ParametricCurve,*/
/*     GEOM_Plane,*/
/*     GEOM_Nurb,*/
/*     GEOM_Cylinder,*/
/*     GEOM_Sphere,*/
/*     GEOM_Cone,*/
/*     GEOM_Torus,*/
/*     GEOM_ParametricSurface*/
/* };*/



typedef class MeshDataId * pMeshDataId;
typedef class AttachableData *pAttachableData;
typedef class AttachDataId *pAttachDataId;


typedef struct MeshDataId * pMeshDataId;
typedef struct AttachDataId *pAttachDataId;
typedef struct AttachableData *pAttachableData;




/* supported local mesh modification types */
typedef enum modType{
  ECOLAPS=0, 
  FSWAP, 
  ESWAP, 
  MOTION, 
  RCOLAPS, 
  ESPLIT, 
  R_REFINE,
  F_REFINE,
  E_REFINE, 
  RCOLAPS_2, 
  SPLTCLPS
} modType;

/* callback functions activated by Local mesh modifications */
typedef void (*CBFunction)(pPList, pPList, void *, modType, pEntity);
typedef void (*CBFunc_move)(pVertex, double *, void *);

/* these must match the callback types in MeshDataId.h */
/* enum { CBdelete, CBmigrateOut, CBmigrateIn };*/
typedef int (*CBfunc)(pAttachableData, pAttachDataId, int, void**, void*);

  /* these data operators are obsolete and retained for compatibility */
  /* new code should use the new operators defined below */
void EN_attachDataP(pEntity , const char *tag, void *data);
void * EN_dataP(pEntity , const char *tag);
int EN_modifyDataP(pEntity, const char *tag, void * data);

void EN_attachDataI(pEntity , const char *tag, int data);
int EN_dataI(pEntity , const char *tag);
int EN_modifyDataI(pEntity, const char *tag, int data);

void EN_removeData(pEntity , const char *tag);

  /* new data operators */
pMeshDataId MD_newMeshDataId(const char *tag);
pMeshDataId MD_lookupMeshDataId(const char *tag);
void MD_deleteMeshDataId(pMeshDataId id);
void MD_setMeshCallback(pMeshDataId id, int event, CBfunc f, void *cbdata);
void MD_removeMeshCallback(pMeshDataId id, int event);

void EN_attachDataInt(pEntity ent, pMeshDataId id, int value);
void EN_attachDataDbl(pEntity ent, pMeshDataId id, double value);
void EN_attachDataPtr(pEntity ent, pMeshDataId id, void * value);

void EN_deleteData(pEntity ent, pMeshDataId id);

int EN_getDataInt(pEntity ent, pMeshDataId id, int *value);
int EN_getDataDbl(pEntity ent, pMeshDataId id, double *value);
int EN_getDataPtr(pEntity ent, pMeshDataId id, void **value);

void EN_modifyDataInt(pEntity ent, pMeshDataId id, int value);
void EN_modifyDataDbl(pEntity ent, pMeshDataId id, double value);
void EN_modifyDataPtr(pEntity ent, pMeshDataId id, void * value);

int EN_inClosure(pEntity, pEntity);

  /* some model stuff */
void GM_name(pGModel model, char *name);
void GM_modeler(pGModel model, char *mname);

  /* mesh entity creation routines */
pRegion M_createR(pMesh mesh, int nFace, pFace *faces, int *dirs, pGEntity gent);
pFace M_createF(pMesh mesh, int nEdge, pEdge *edges, int *dirs, pGEntity gent);
pEdge M_createE(pMesh mesh, pVertex v1, pVertex v2, pGEntity ent);
pVertex M_createVP(pMesh mesh, double x, double y, double z, double *param,
		    int patch, pGEntity ent);
pVertex M_createVP2(pMesh mesh, double *xyz, double *param,
                     int patch, pGEntity ent);
pPoint M_createP(pMesh mesh, double x, double y, double z, double *param,
		   int unused, pGEntity ent);

  /* mesh entity deletion routines */
void M_removeRegion(pMesh, pRegion region);
void M_removeFace(pMesh, pFace face);
void M_removeEdge(pMesh, pEdge edge);
void M_removeVertex(pMesh, pVertex vertex);
void M_removePoint(pMesh, pPoint point);

  /* point routines */
pPoint P_new(void);
void P_delete(pPoint);

void P_setPos(pPoint  , double x, double y, double z);
void P_setParam1(pPoint, double param);
void P_setParam2(pPoint, double p1, double p2, int ptch);

  /* obsolete iterator functions */
pRegion M_nextRegion(pMesh, void **restart);
pFace M_nextFace(pMesh, void **restart);
pEdge M_nextEdge(pMesh, void **restart);
pVertex M_nextVertex(pMesh, void **restart);
pPoint M_nextPoint(pMesh, void **restart);

pRegion M_nextRegionCancel(pMesh, void **restart);
pFace M_nextFaceCancel(pMesh, void **restart);
pEdge M_nextEdgeCancel(pMesh, void **restart);
pVertex M_nextVertexCancel(pMesh, void **restart);
pPoint M_nextPointCancel(pMesh, void **restart);

  /* other functions */
typedef double dArray[3];  /* for hp compiler compatibility */

double V_size(pVertex v);
void V_setSize(pVertex v, double s);
double E_dihedral(dArray *xyz);
double E_length(pEdge edge);
double P_distToLine(double *pxyz, double *v0, double *v1, double *normal);
void crossProd(double [],double [],double []);
int R_shape(pRegion,double *);
pPList R_verticesLeft(pRegion);
int XYZ_shape(dArray *xyz,double *);
pPList E_swpCfgVerts(pEdge, int);
int XYZ_checkFlat(dArray *xyz);
pVertex F_oppositeVertex(pFace face, pEdge edge);
void M_setTolerance(pMesh);
double M_getTolerance();
void diffVt(double [],double [],double []);
double dotProd(double [],double []);
int R_checkFlat(pRegion);
pEdge F_vtOpEd(pFace face, pVertex vert);
pFace   R_vtOpFc(pRegion ,pVertex);
double XYZ_volume(dArray *xyz);
pVertex R_fcOpVt(pRegion rgn, pFace face);
int F_numPoints(pFace);
double XYZ_distance2(double [3], double [3]);
pFace F_exists(eType, pEntity, pEntity, pEntity, pEntity);

int C_closestPoint(int type, void *tag, double *inxyz, int seedflag,
                    double *seedxyz, double *seedpar, double *outxyz,
                    double *outparm);
int C_gmtype(int type, void *tag); /* Geometric entity type */

int GF_periodic(pGFace f, int dir);
int GE_periodic(pGEdge e);

int GF_geomDirection(pGFace);
int GE_geomDirection(pGEdge e);

int GF_paramDegeneracies(pGFace f, int dir, double *par);

int F_swap(pMesh,pFace, double,pPList*,pPList*);
int F_evalSwap(pFace,double);

//F_conToFace function in not include in Simmetrix 7 libs
inline int F_conToFace(pFace face1, pFace face2)
{
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


pEdge E_exists(pVertex v1, pVertex v2);

int E_numSwpCfg(pEdge);
int E_evalSwpCfg(pEdge, int);
int E_swap(pMesh, pEdge, int, pPList *);

int E_chkClpTopo(pVertex vdel, pEdge ledge);
int E_evalColaps(pMesh mesh, pEdge ledge, pVertex vertd, pVertex vertr);
int E_colaps(pMesh mesh, pEdge ledge, pVertex vertd, pVertex vertr, 
             pPList oldreg,pPList *newReg);

pVertex E_split(pMesh OM_mesh, pEdge edge, double *xyz, double *parm, 
                pPList *newReg);
pVertex F_split(pMesh OM_mesh, pFace face, double *xyz, double *parm, 
                pPList *newReg);

void V_merge(pVertex, pVertex);

void F_normalVector(pFace face, int dir, double* normal);
pEntity XYZ_onBoundary(double vxyz[3],pEntity entity,double tol);
pEdge R_gtOppEdg(pRegion region, pEdge edge);
void F_setWhatIn(pFace  , pGEntity what);
void F_chDir(pFace);
void E_setWhatIn(pEdge  , pGEntity what);
void V_setWhatIn(pVertex, pGEntity what); /*sets classification of vertex */

int M_checkAdj(pMesh mesh);

void E_setPoint(pEdge , pPoint pt);
void V_setPoint(pVertex, pPoint pt);

int E_evalColapsOnGFace(pMesh mesh, pEdge cedge, pVertex vertd, pVertex vertr,
                        double *normal);
int E_colapsOnGFace(pMesh mesh, pEdge ledge, pVertex vertd, pVertex vertr, 
		    pPList oldfaces,pPList *modfaces);

int E_evalSwpOnGFace(pEdge swpedge, pGEntity gface, double *normal);
pEdge E_swapOnGFace(pMesh mesh, pEdge swpedge, pGEntity gface, pPList *newf);


/*************************************************/
/* operators to deal with constrained operations */
/*************************************************/
typedef enum opType {
  SWAP, COLAPS, SPLIT, MOVE, DELETE
} opType;

void EN_constrain(opType ,pEntity ) ;
void EN_constrainAll(pEntity ) ;
void EN_unconstrain(opType ,pEntity ) ;
void EN_unconstrainAll(pEntity);
int  EN_okTo(int ,pEntity ) ;
int GEN_okTo(int i,pGEntity g);

void GEN_constrain(opType ,pGEntity ) ;
void GEN_constrainAll(pGEntity ) ;
void GEN_unconstrain(opType ,pGEntity ) ;
void GEN_unconstrainAll(pGEntity);



/*************************************************/
/*  Debugging utilities                          */
/*************************************************/
void M_info(pMesh);
void EN_info(pEntity);
void R_info(pRegion);
void F_info(pFace);
void E_info(pEdge);
void V_info(pVertex);
void PList_printx(pPList);

int triangulatecavity(pMesh, pPList, pGRegion*, int, pPList, pPList);
void PList_app2Lst(pPList*, pPList, int);
typedef pPList pPQueue;
pPQueue PQueue_new(void);
void PQueue_delete(pPQueue);
void PQueue_reset(pPQueue);
void PQueue_push(pPQueue, pEntity);
void PQueue_remItem(pPQueue, pEntity);
void PQueue_pushUnique(pPQueue, pEntity);
pEntity PQueue_pop(pPQueue);
int PQueue_size(pPQueue q);

struct IntDescription { double xyz[3]; pEntity IntEnt1, IntEnt2;
                        struct IntDescription *next, *prev; };
typedef struct { int dim, nbr; struct IntDescription *Description;
                 pPList boundedEntList; } IntersectionData;
IntersectionData M_intersectionData_new();
void M_intersectionData_delete(IntersectionData);
int xyz_inIntersectionData(IntersectionData*, double[3], double);
void M_updateIntersectionData(IntersectionData*, pEntity, pEntity, double[3], double);
pEdge M_boundedEdgeOnFace(pFace, pVertex, pVertex);
pEdge M_boundedEdgeOnRegion(pRegion, pVertex, pVertex);
pFace M_boundedFaceOnRegion(pRegion, pEntity, pEntity);
int M_intersectEdges(pEdge, int, pEdge, int, IntersectionData*);
int F_intersectEdge(pFace, int, pEdge, int, IntersectionData*);
int M_intersectFaces(pFace, int, pFace, int, IntersectionData*);
int R_intersectEdge(pRegion, int, pEdge, int, IntersectionData*);
int R_intersectFace(pRegion, int, pFace, int, IntersectionData*);
int M_intersectRegions(pRegion, int, pRegion, int, IntersectionData*);
int M_intersectRay3XYZ(double[3], double[3], double[3][3], double[3], int*);
int M_intersectRay2XYZ(double[3], double[3], double[2][3], double[3], int*);
int M_intersectRayXYZ(double[3], double[3], double[3][3], double[3], int*);
int M_intersectLine3XYZ(double[3], double[3], double[3][3], double[3], int*);
int inList(pPList, void*);
pPList V_bdryFaces(pVertex);
int F_checkFlat(pFace, double*);
void F_coord(pFace, dArray*);

double GEN_tolerance(pGEntity e);
void GEN_reparam(pGEntity target, pGEntity ent, double *entPar, int entDir,
		 double *targetPar);
int GEN_inClosure(pGEntity target, pGEntity ent);
int GF_periodic(pGFace f, int dir);
int GE_periodic(pGEdge e);

void GEN_attachDataP(pGEntity , const char *tag, void *data);
void * GEN_dataP(pGEntity , const char *tag);
int GEN_modifyDataP(pGEntity, const char *tag, void * data);
void GEN_attachDataI(pGEntity , const char *tag, int data);
int GEN_dataI(pGEntity , const char *tag);
int GEN_modifyDataI(pGEntity, const char *tag, int data);
void GEN_removeData(pGEntity , const char *tag);

int MS_tri2quad(pMesh mesh, int type);


double R_inscrRad(pRegion);
void XYZ_FAngles(dArray *, double *); 
void XYZ_dhdAngs(double xyz[4][3], double *cosAngs);
int C_equal(double real1,double real2);
int E_rgnLocalId(pRegion, pEdge);
void R_dihedral(pRegion reg, double *small, double *big) ;
void R_coord(pRegion rgn,dArray *xyz);
int normVt(double [], double []);
pPList E_faces(pEdge ); 
double R_volume(pRegion);
double P_distToPlane (double* pxyz, double* v0, double* v1, double* v2, double *snorm);
void P_projOnTriPlane(double fxyz[3][3], double pxyz[3], double *proxyz, 
		      double normal[3], double *distance);
		      
pPList delete_R(pMesh mesh, pRegion region);
pPList delete_F(pMesh mesh, pFace face);
pPList delete_E(pMesh mesh, pEdge edge);
void delete_V(pMesh mesh, pVertex v);


#else  /* SIMINTERNAL */
#include "MeshSimInternalPrivate.h"
#endif  /* SIMINTERNAL */

#endif /* not H_MeshSimInternal */
#endif
