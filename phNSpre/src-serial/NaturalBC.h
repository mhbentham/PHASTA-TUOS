#ifndef H_NaturalBC
#define H_NaturalBC

//
// Natural Boundary condition class
//

#include "BoundaryCondition.h"
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#endif

enum { MF, NP, TV, HF, F1, F2, F3, F4, SID, TW,
	   VF, XF, YF, XYF, XXF, YYF };  // for readability
class NaturalBC : public BoundaryCondition {
public:
  NaturalBC (pGFace);
  int evalFace (pFace , double *, int *,int nflx);
  NaturalBC (pGEdge);
  int evalEdge (pEdge , double *, int *,int nflx);
};


#endif
