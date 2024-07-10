#ifndef H_PeriodicBC
#define H_PeriodicBC

#include "BoundaryCondition.h"

class PeriodicBC : public BoundaryCondition {
public:
    PeriodicBC (pGModel , pGFace );        // constructor for face
    PeriodicBC (pGModel , pGEdge );        // constructor for edge
    PeriodicBC (pGModel , pGVertex );      // constructor for vert
    int getPerMaster (pGFace *);        // returns num jumps, also entity
    int getPerMaster (pGEdge *);        // in parentheses is the master
    int getPerMaster (pGVertex *);      // if present ent is vertex then must
                                        //   call getPerMaster(Vertex *)
    double getAngle (pGFace, pMesh);
    double getAngle (pGEdge );
    double getAngle (pGVertex );
    static int global_sanity;
    
private:
    gType gtype;
    pGFace gf;
    pGEdge ge;
    pGVertex gv;
    pGModel model;
    double myangle; /* so that I don't have to pass stuff around */
    void update_inherit ();
    double getDistance ( pGEdge );
    double getDistance ( pGVertex );
    double getDistance (double *, double *, double theta=0.0);
    double SeeAngle ( void );
    void GF_normal_flat(pGFace, pMesh, double *xyz);
};

#endif
