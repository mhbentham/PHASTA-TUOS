#ifndef H_BoundaryCondition
#define H_BoundaryCondition

//
// Boundary condition class
//

#include <stdlib.h>
#ifdef SIM
#include "MeshSim.h"
#include "MeshSimInternal.h"
#include "SimAttribute.h"
#else
#include "AOMD.h"
#include "myAttribute.h"
#ifdef MESHMODEL
#include "modeler.h"
#endif
#endif

#define MLEN 25

#define SQ(x)       ((x)*(x))
#define SMALL       1.0e-5
#define ABS(x)      ((x) < 0 ? -(x) : (x))
#define CLOSE(x,y)  (ABS((x)-(y)) < SMALL ? 1 : 0)
#define false       0
#define true        1

/* enum { false, true }; */

class BoundaryCondition {
public:
  BoundaryCondition ();
  virtual ~BoundaryCondition () ;
  int isSet ();
  int isAttSet (int i);

protected:
  pAttribute *AttList;
  int set;
  int dontinherit;            // tells you what (not) to inherit - not for face

};

#endif


