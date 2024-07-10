//
// Initial condition class
//
#ifdef SIM
#include "MeshSim.h"
#include "SimAttribute.h"
#else
#include "AOMD.h"
#include "myAttribute.h"
#endif

class InitialCondition {
public:
  InitialCondition(pGModel );
  void eval(pVertex , double *);
  void eval(pEdge , double *, int);
private:
  pAttribute attList[13]; //  two thermo + velocity vector (3) + 4 scalars + 6 for c1 quintic triangles
};
