/* fortran wrapper for C function TetShapeAndDrv which returns
   the shape functions and their derivatives for tets
   */
#include <FCMangle.h>

int TetShapeAndDrv(int p,double par[3],double N[],double dN[][3]);
#define shptri FortranCInterface_GLOBAL_(shptri, SHPTRI)

shptri(int *p, double par[], double N[], double dN[][3])
{
 TriShapeAndDrv(*p,par,N,dN);

}
