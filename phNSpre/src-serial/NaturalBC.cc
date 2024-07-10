#include <iostream>
#include <math.h>
#include "nspre_data.h"
#include "nspre_functions.h"
#include "NaturalBC.h"

extern int ensa_dof;  // bring in the # of variables to help find # of scalars
extern int nscalar;
extern bool c1quintic;

#define numAttsN 16        // number of natural bc attributes

char strAttN[numAttsN][MLEN] = { "mass flux",
                                 "natural pressure", 
                                 "traction vector", 
                                 "heat flux",
                                 "scalar_1 flux",
                                 "scalar_2 flux",
                                 "scalar_3 flux",
                                 "scalar_4 flux",
                                 "surf ID",
                                 "turbulence wall",
				 "value_sc flux",
				 "x_derivative flux",
				 "y_derivative flux",
				 "xy_derivative flux",
				 "xx_derivative flux",
				 "yy_derivative flux" };

// In 3D natural boundary condition is on faces
NaturalBC::NaturalBC(pGFace gface) : BoundaryCondition()     // Natural BC's
{
  int i;
#ifdef SIM
  if (!GEN_dataI((pGEntity)gface,"NoBI"))  {
#else
  int tmp;
  if (!GEN_dataI((pGEntity)gface,"NoBI",&tmp))  {
#endif
   // model face not flagged as not in BI
    this->set = true;                            // set is true

    AttList = new pAttribute[numAttsN];          // update attlist
    for (i=0; i < numAttsN; i++)  {
      AttList[i] = GEN_attrib((pGEntity)gface, strAttN[i]);
    }
  }
}

// In 2D natural boundary condition is on edges
NaturalBC::NaturalBC(pGEdge gedge) : BoundaryCondition()     // Natural BC's
{
  int i;
#ifdef SIM
  if (!GEN_dataI((pGEntity)gedge,"NoBI"))  {
#else
  int tmp;
  if (!GEN_dataI((pGEntity)gedge,"NoBI",&tmp))  {
#endif
   // model face not flagged as not in BI
    this->set = true;                            // set is true

    AttList = new pAttribute[numAttsN];          // update attlist
    for (i=0; i < numAttsN; i++)  {
      AttList[i] = GEN_attrib((pGEntity)gedge, strAttN[i]);
    }
  }
}

// In 3D natural boundary condition is on faces
int NaturalBC::evalFace(pFace face, double *BCB, int *iBCB,int nflx)  
{
  pVertex vertex[4];
  int nenb = F_numEdges(face), i, j, ibcb,k;

  /* right now we are using nenb for interpolation of Natural Boundary
     conditions . This will have to be upgraded to bnshlb when we decide 
     to make higher order interpolation of NBCs */

  /* as of 3/13/2001 we have decided to go with piecewise constant natural 
     boundary conditions since we never seem to use anything else. we save 
     trouble moving to higher order and also avoid a potential bug on linear 
     ones when we have varying bcb for quadface wedges and triface pyramids

     so in short.. don't pay much attention to the last comment made
     above the current one but just keep in mind that we have natural
     boundary conditions which are specified as constant on a face.

     by the way this function evaluates them at the centroid. */ 


  double xyz[3];
  pPList vlist = F_vertices(face, 1);
  double centroid[3]={0.0, 0.0, 0.0};

  for (i=0; i < PList_size(vlist); i++)  {
    vertex[i] = (pVertex) PList_item(vlist, i);
  }
  PList_delete(vlist);

  for(i=0; i < nenb; i++) {
    V_coord(vertex[i], xyz);
    for(j=0; j< 3; j++) centroid[j]+= xyz[j];
  }
  for(i=0; i< 3; i++) centroid[i]/=nenb;

  if (AttList[MF])                    // mass flux
    BCB[0] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[MF], centroid);

  if (AttList[NP])                   // natural pressure
    BCB[1] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[NP], centroid);

  if (AttList[TV])                     // traction
      for( int index =0; index < 3; index++)
          BCB[2+index] = AttributeTensor1_evalDS((pAttributeTensor1)AttList[TV], index,
                                                 centroid);
  
  if (AttList[HF])                    // heat flux
      BCB[5] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[HF], centroid);
  
  for( k=0; k < nflx-5; k++)
      if (AttList[F1+k])                  // scalar flux
          BCB[6+k] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[F1+k],centroid);
  
  ibcb = 0;                            // Updating ibcb for this one
  ibcb += (AttList[MF]) ? 1 : 0;
  ibcb += (AttList[NP]) ? 2 : 0;
  ibcb += (AttList[TV]) ? 4 : 0;
  ibcb += (AttList[HF]) ? 8 : 0;
  ibcb += (AttList[TW]) ? 16 : 0;
  for( k=0; k < nflx-5; k++) 
  ibcb += (AttList[F1+k]) ? (int)(32*pow(2.0,k)) : 0;
  
  iBCB[0]=ibcb;
  iBCB[1]=0;
  if(AttList[SID]) {
    iBCB[1] = (int) AttributeTensor0_evalDS((pAttributeTensor0)AttList[SID], centroid);
    //
    //  above we seemingly misuse a real attribute to obtain an
    //  integer instead of the "right-use" shown below.  We do this
    //  because we would like to be able to vary the srfID on a face
    //  without actually cutting the model face (i.e. currently used
    //  for DtN BC's).
    //    iBCB[1] = Att_int(AttList[SID]);
  }
  
  return ibcb;
}

// In 2D natural boundary condition is on edges
int NaturalBC::evalEdge(pEdge edge, double *BCB, int *iBCB,int nflx)  
{
  pVertex vertex[2];
  int nenb = 2, i, j, ibcb,k;

  double xyz[3];
  double centroid[3]={0.0, 0.0, 0.0};

  for (i=0; i < 2; i++)  {
    vertex[i] = E_vertex(edge, i);
  }

  for(i=0; i < nenb; i++) {
    V_coord(vertex[i], xyz);
    for(j=0; j< 3; j++) centroid[j]+= xyz[j];
  }
  for(i=0; i< 3; i++) centroid[i]/=nenb;

  if (AttList[MF])                    // mass flux
    BCB[0] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[MF], centroid);

  if (AttList[NP])                   // natural pressure
    BCB[1] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[NP], centroid);

  if (AttList[TV])                     // traction
      for( int index =0; index < 3; index++)
          BCB[2+index] = AttributeTensor1_evalDS((pAttributeTensor1)AttList[TV], index,
                                                 centroid);
  
  if (AttList[HF])                    // heat flux
      BCB[5] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[HF], centroid);
  
  if (c1quintic){
     for( k=0; k < 6; k++)
        if (AttList[VF+k])                  // scalar & derivative flux for c1 
           BCB[6+k] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[VF+k],centroid);
  } else {
     for( k=0; k < nscalar; k++)
        if (AttList[F1+k])                  // scalar flux
           BCB[6+k] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[F1+k],centroid);
  }	   
  
  ibcb = 0;                            // Updating ibcb for this one
  ibcb += (AttList[MF]) ? 1 : 0;
  ibcb += (AttList[NP]) ? 2 : 0;
  ibcb += (AttList[TV]) ? 4 : 0;
  ibcb += (AttList[HF]) ? 8 : 0;
  ibcb += (AttList[TW]) ? 16 : 0;
  if (c1quintic){
     for( k=0; k < 6; k++) ibcb += (AttList[VF+k]) ? (int)(32*pow(2.0,k)) : 0;
  } else {
     for( k=0; k < nscalar; k++) ibcb += (AttList[F1+k]) ? (int)(32*pow(2.0,k)) : 0;
  }
  
  iBCB[0]=ibcb;
  iBCB[1]=0;
  if(AttList[SID]) {
    iBCB[1] = (int) AttributeTensor0_evalDS((pAttributeTensor0)AttList[SID], centroid);
  }
  
  return ibcb;
}
