#include <cstdlib>
#include <iostream>
#include "InitialCondition.h"
#include "nspre_functions.h"

extern int ensa_dof, rRead, nscalar;  
extern bool c1quintic;
extern globalInfo* info;

// Right now user must be very careful to see that what he is reading
// out of the restart file complemets what is set through the GUI

// If anything is read twice , the GUI takes precedence

InitialCondition::InitialCondition(pGModel model)
{
  //static int eflg[3] ={ 1, 1 ,1};
  static int stick=1;
  char* flag[13]={"initial pressure", "initial velocity", "initial temperature",
                 "initial scalar_1", "initial scalar_2", "initial scalar_3",
                 "initial scalar_4", "initial value_sc", "initial x_derivative",
		 "initial y_derivative", "initial xy_derivative",
		 "initial xx_derivative", "initial yy_derivative" };

  attList[0] = GM_attrib(model, "initial pressure");
  attList[1] = GM_attrib(model, "initial velocity");
  attList[2] = GM_attrib(model, "initial temperature");
  attList[3] = GM_attrib(model, "initial scalar_1");
  attList[4] = GM_attrib(model, "initial scalar_2");
  attList[5] = GM_attrib(model, "initial scalar_3");
  attList[6] = GM_attrib(model, "initial scalar_4");
  attList[7] = GM_attrib(model, "initial value_sc");
  attList[8] = GM_attrib(model, "initial x_derivative");
  attList[9] = GM_attrib(model, "initial y_derivative");
  attList[10] = GM_attrib(model, "initial xy_derivative");
  attList[11] = GM_attrib(model, "initial xx_derivative");
  attList[12] = GM_attrib(model, "initial yy_derivative");

  int nbini;
  if (info->nsd == 3 ) {
    nbini = ensa_dof-2;
  } else {
    nbini = ensa_dof-1;
  }
  if (c1quintic == true) {
    for(int i =0; i< 3; i++) if (attList[i]){
  //    if(eflg[i]) {
	cout<< flag[i] <<" Has been set in the GUI "<<endl;
  //      eflg[i] = 0;
  //    }
    }
    for(int i =7; i< 13; i++) if (attList[i]){
 //     if(eflg[i]) {
	cout<< flag[i] <<" Has been set in the GUI "<<endl;
 //       eflg[i] = 0;
 //     }
    }
  } else {
    for(int i =0; i< nbini; i++) if (attList[i]){
  //    if(eflg[i]) {
	cout<< flag[i] <<" Has been set in the GUI "<<endl;
  //      eflg[i] = 0;
   //   }
    }
  } 

  if ((!attList[0] || !attList[1] || !attList[2]) && (rRead == 0)) {
    cerr << "\nNSpre Error: you must specify an initial condition or use "
         << "the restart option" << endl;
    exit(-1);
  }
  if ( stick ){ 
     cout <<" The Rest of the Initial Conditions (if any left) "<<endl
          <<" Should be coming from the restart file " << endl; 
     stick = 0;
  }
}

void InitialCondition::eval(pVertex vertex, double *q)
{
  double x[3];
//  int nsc;
 
  V_coord(vertex,x);
  if (info->nsd == 3) {
//    nsc=ensa_dof-5;

    if(attList[0]) q[0] = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x);
    if(attList[1]) {
      q[1] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 0, x);
      q[2] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 1, x);
      q[3] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 2, x);
    }
    if(attList[2]) q[4] = AttributeTensor0_evalDS((pAttributeTensor0) attList[2], x);
    for (int i=0; i < nscalar; i++)
      if(attList[3+i]) q[5+i] = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x);
  } else {
//    nsc=ensa_dof-4;

    if(attList[0]) q[0] = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x);
    if(attList[1]) {
      q[1] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 0, x);
      q[2] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], 1, x);
    }
    if(attList[2]) q[3] = AttributeTensor0_evalDS((pAttributeTensor0) attList[2], x);
    for (int i=0; i < nscalar; i++)
      if(attList[3+i]) q[4+i] = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x);
    if (c1quintic == true){
      for (int i=0; i<6; i++)
         if(attList[7+i]) q[4+i] = AttributeTensor0_evalDS((pAttributeTensor0)attList[7+i], x);
    }
  }
}

// not corrected for 2D or c1 quintic triangles
void InitialCondition::eval(pEdge edge, double *q, int p)
{
//  int nsc = ensa_dof-5;

  if (p > 2){
    for (int i=0; i < 5+nscalar; i++) q[i]=0.0;
    return;
  }

  pVertex v1 = E_vertex(edge,0);
  pVertex v2 = E_vertex(edge,1);

  double x1[3], x2[3];
  V_coord(v1, x1);
  V_coord(v2, x2);
  double x3[3] = {0.5*(x1[0]+x2[0]), 0.5*(x1[1]+x2[1]), 0.5*(x1[2]+x2[2])};
  double d1,d2,d3,vel1[3],vel2[3],vel3[3];

  if(attList[0]) {
    d1 = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x1);
    d2 = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x2);
    d3 = AttributeTensor0_evalDS((pAttributeTensor0)attList[0], x3);
    q[0] = (d1 + d2 - 2.0 * d3);
  }
  
  if(attList[1]){
    for(int i=0; i<3; i++) {
      vel1[i] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], i, x1);
      vel2[i] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], i, x2);
      vel3[i] = AttributeTensor1_evalDS((pAttributeTensor1)attList[1], i, x3);
    }
    q[1] = (vel1[0] + vel2[0] - 2.0 * vel3[0]);
    q[2] = (vel1[1] + vel2[1] - 2.0 * vel3[1]);
    q[3] = (vel1[2] + vel2[2] - 2.0 * vel3[2]);
  }

  if(attList[2]){
    d1 = AttributeTensor0_evalDS((pAttributeTensor0)attList[2], x1);
    d2 = AttributeTensor0_evalDS((pAttributeTensor0)attList[2], x2);
    d3 = AttributeTensor0_evalDS((pAttributeTensor0)attList[2], x3);
    q[4] = (d1 + d2 - 2.0 * d3);
  }
  
  for (int i=0; i < nscalar && attList[3+i] ; i++) {
    if (attList[3+i]){
      d1 = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x1);
      d2 = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x2);
      d3 = AttributeTensor0_evalDS((pAttributeTensor0)attList[3+i], x3);
      q[5+i] =(d1 + d2 - 2.0 * d3);
    }
  }
}

