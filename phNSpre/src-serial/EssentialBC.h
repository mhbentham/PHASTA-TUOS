#ifndef H_EssentialBC
#define H_EssentialBC

#include <iostream>
using namespace std;
#include "BoundaryCondition.h"

#ifndef SIM
#include "nspre_functions.h"
#endif

#define numAttsE 16                       // number of essential bc attributes
enum { thermo, velo, scalar, axisym, bellf };        // to set dontinherit accordingly
enum { D, T, P, C1, C3, S1, S2, S3, S4, TBI, 
       VSC, XD, YD, XYD, XXD, YYD };//for readability of essential stuff

class EssentialBC : public BoundaryCondition {
public:
  EssentialBC (pGFace gface);
  EssentialBC (pGEdge gedge);
  EssentialBC (pGVertex gvert);
  int eval (pVertex, double *);
  int evalEdge (pEdge, double *, int);      // for hierarchic basis  
  int evalFace (pFace, double *, int);      // for hierarchic basis
  void takeBCfromIC(double *BC, double *qTot, int GDOFnum);
private:
  gType gtype;                            // Gvertex, Gedge or Gface
  pGFace gf;     // could store in a union since either face, edge or vertex
  pGEdge ge;     //   but too much bother for insignificant memory savings
  pGVertex gv;
  int zxcl;      // Flag to do Axisymmetric Centerline

  void update_inherit ();     // to update dontinherit
  void update_thermo(const double*, double*, int*);  // thermo BC's
  void update_velo(const double*, double*, int*);    // velo BC's
  void update_scalar(const double*, double*, int*);  // scalar BC's
  void update_bellf(const double*, double*, int*);  // bellf BC's (c1 quintic)
  void update_axisym(const double*, double*, int*);  // axisym centerline
  void pteval (const double*, double *, int *);
  int isZA()
  {
   GEntity* gty ;
   switch (gtype){
   case Gface:
	gty = (pGEntity) gf;
        break;
   case Gedge:
        gty = (pGEntity) ge;
        break;
   case Gvertex:
        gty = (pGEntity) gv;
        break;
   default:
        cerr << " Entity does not have type "<<endl;
        exit(1);   
   }
#ifdef SIM
   return GEN_dataI(gty,"AXCL");
#else
   int tmp;
   return GEN_dataI(gty,"AXCL",&tmp);
#endif
  }
  
};

#endif
