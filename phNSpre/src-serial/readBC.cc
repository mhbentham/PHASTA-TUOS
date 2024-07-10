// 
//
// This function reads the internal boundary condition information and
// assigns it to the appropriate mesh vertices
//
// Internal boundary conditions are necessary when the user wants to 
// set boundary condtions on mesh entities which aren't classified on
// a geometric boundary. For example, a trip strip at the inflow of
// a turbulent cavity is produced by setting boundary conditions on
// a row of nodes near the base of the cavity. The user must specify
// these nodes in a separate file.
// 
#include <fstream>
#include <stdlib.h>
#include "nspre_data.h"
#include "nspre_functions.h"

using namespace std;

void readBC(pGModel model, pMesh mesh)
{
  int      nn,ibc,num,i;
  pVertex  vertex;
  pVertex  *vArray = new pVertex [M_numVertices(mesh)];
  pGEntity  g_ent;
  double   *list; // for BC values 
  void     *temp;

  // create an array of all vertices in the mesh
  temp=0;
  i=0;
  VIter vIter = M_vertexIter(mesh);
  while (vertex =  VIter_next(vIter)) vArray[i++]=vertex;
  VIter_delete(vIter);

  /******************** internal BC's ********************/
  /* To specify boundary conditions on internal nodes,   */
  /* such as a trip-strip, a file must be provided that  */
  /* contains:                                           */
  /*   nn                                                */
  /*   11 BC values                                      */
  /*   iBC-code node# (for each node)                    */
  /* which should be named: inodes.dat                   */
  /*******************************************************/
  ifstream fb( "inodes.dat" );
  list = (double *)calloc(11,sizeof(double));
  fb >> nn ;
  for (i=0; i < 11; i++) fb >> list[i];

  GRIter griter = GM_regionIter(model);
  g_ent = (GEntity *)GRIter_next(griter);
  GRIter_delete(griter);

  // attach the BC values to the first geometric region
  GEN_attachDataP(g_ent,"bc  ",(void *)list);
    
  // read the internal nodes, and attach the new ibc code 
  for (i=0; i < nn; i++){
    fb >> ibc >> num;
    vertex = vArray[num-1];
    EN_attachDataI((pEntity)vertex,"ibc ",ibc);
  }
  fb.close();
} 
