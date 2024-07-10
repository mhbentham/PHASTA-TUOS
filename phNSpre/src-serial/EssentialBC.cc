#include "EssentialBC.h"
#include <math.h>
#include "nspre_data.h"

// iBC and BCtmp array is evaluated in these routines
// BCtmp goes this  rho p T c11 c12 c13 m1 c21 c22 c23 m2 theta sc_j
// with indecies     0  1 2  3   4   5  6   7   8   9  10   11  12  13  14  15  (3D)
//                   0  1 2  3   4      5                   6    7  8  9  10    (2D)
// add 1 for phasta numbering
// iBC bits go like this rho t p u  v w sc1 sc2 sc3 sc4 perio spebc
//                        1  2 4 8 16 32 64 128 256 512 1024  2048 
// same for 3D and 2D
//
// 2D c1 quintic
// BCtmp: rho p T c11 c12 m1 theta val_sc x y xy xx yy der.
//         0  1 2  3   4  5    6    7     8 9 10 11 12
// iBC bits: rho t p u  v w val_sc x y   xy der. perio spebc xx   yy der.
//            1  2 4 8 16 32 64  128 256 512     1024  2048  4096 8192
//            0  1 2 3  4  5  6    7   8   9       10    11    12   13
//
// Modified for 2D meshes with c1 quintic triangles by Elaine Bohr March 2006

#define SQ(x)       ((x)*(x))
#define SMALL       1.0e-5
#define ABS(x)      ((x) < 0 ? -(x) : (x))
#define CLOSE(x,y)  (ABS((x)-(y)) < SMALL ? 1 : 0)

extern int nshgTot, nscalar;
extern globalInfo* info;
extern bool c1quintic;
extern int ensa_dof;  // bring in the # of variables to help find # of scalars
char strAttE[numAttsE][MLEN] = { "density", 
                                 "temperature", 
                                 "pressure", 
                                 "comp1",
                                 "comp3",
                                 "scalar_1",
                                 "scalar_2",
                                 "scalar_3",
                                 "scalar_4",
                                 "take bc from ic",
				 "value_sc",
				 "x_derivative",
				 "y_derivative",
				 "xy_derivative",
				 "xx_derivative",
				 "yy_derivative"};


int i;
EssentialBC::EssentialBC(pGFace gface) : BoundaryCondition()   // face cons.
{
  int i, loop = true;


  for (i=0; i < numAttsE && loop; i++)  {
    if (GEN_attrib((pGEntity)gface,strAttE[i]))  loop = false; // at least one
  }

  if (!loop)  {                 // at least one is set => "good"(??AK) face
    this->set = true;           // set is true
    this->dontinherit = 0;      // dontinherit nothing for face
    this->gtype = Gface;        // type of present entity
    this->gf = gface;           // pointer to entity

    AttList = new pAttribute[numAttsE];          // update attlist
    for (i=0; i < numAttsE; i++) {
        AttList[i] = GEN_attrib((pGEntity)gface,strAttE[i]);
    }
  }
}

EssentialBC::EssentialBC(pGEdge gedge) : BoundaryCondition()    // edge cons.
{
  int i;

  this->set = true;           // edges "good" by default
  this->dontinherit = 0;      // dontinherit nothing by default
  this->gtype = Gedge;        // type is edge
  this->ge = gedge;           // pointer to that edge

  AttList = new pAttribute[numAttsE];      // update attlist
                                           // many may be null at present

  for (i=0; i < numAttsE; i++){
    AttList[i] = GEN_attrib((pGEntity)gedge,strAttE[i]);
  }
  update_inherit ();          // what to inherit, what not to
}

EssentialBC::EssentialBC(pGVertex gvert) : BoundaryCondition()  // vert. cons.
{
  int i;
  this->set = true;            // model vertices "good" by default
  this->dontinherit = 0;       // dontinherit nothing by default
  this->gtype = Gvertex;       // type is vertex
  this->gv = gvert;            // pointer to this vertex

  AttList = new pAttribute[numAttsE];      // update attlist
                                           // many may be null at present

  for (i=0; i < numAttsE; i++){
    AttList[i] = GEN_attrib((pGEntity)gvert,strAttE[i]);
  }
  update_inherit ();          // what to inherit, what not to
}

void EssentialBC::update_inherit()         // private , tells what to inherit
{
  if (isAttSet(D) || isAttSet(T) || isAttSet(P))  {   // tells if attrib is set
    dontinherit = setbit (dontinherit, thermo);       // dontinherit thermo etc
  }
  if (isAttSet(C1) || isAttSet(C3))  dontinherit = setbit (dontinherit, velo);
  if (isAttSet(S1) || isAttSet(S2) || isAttSet(S3) || isAttSet(S4))  { 
    dontinherit = setbit (dontinherit, scalar);     
  }
  if (isAttSet(VSC) || isAttSet(XD) || isAttSet(YD) || 
      isAttSet(XYD) || isAttSet(XXD) || isAttSet(YYD))  { 
    dontinherit = setbit (dontinherit, bellf);     
  }
  if (isZA()) {
    dontinherit = setbit (dontinherit, axisym);
  }
}

int EssentialBC::eval(pVertex vert, double *BC)       // BC,ibc eval 
{
  double xyz[3];
  int ibc;
  V_coord (vert, xyz);
  this->pteval (xyz, BC, &ibc);
  return ibc;
}


// NOT IMPLEMENTED FOR 2D
int EssentialBC::evalEdge(pEdge edge, double *BC, int p)   // hierarchic basis
{
  pVertex v1 = E_vertex(edge, 0);
  pVertex v2 = E_vertex(edge, 1);
  double xyz1[3], xyz2[3];
  V_coord (v1, xyz1);
  V_coord (v2, xyz2);
  double xyz3[3] = { 0.5*(xyz1[0] + xyz2[0]), 
                     0.5*(xyz1[1] + xyz2[1]), 
                     0.5*(xyz1[2] + xyz2[2]) };
  double d1[11]={0.0}, d2[11]={0.0}, d3[11]={0.0};
  int i, idum, ibc;
  this->pteval(xyz1, d1, &idum);
  this->pteval(xyz2, d2, &idum);
  this->pteval(xyz3, d3, &ibc);                     // midpoint's ibc is edge's

  if (p > 2)  {                                   // all zero except cosines
    BC[0] = BC[1] = BC[2] = BC[6] = BC[10] = 0.0;

  } else {

    for (i=0; i<3; i++)  BC[i] = d1[i] + d2[i] - 2*d3[i];
    BC[6] = d1[6] + d2[6] - 2*d3[6]; 
    BC[10] = d1[10] + d2[10] - 2*d3[10]; 
  }
  for (i=3; i<6; i++)  BC[i] = d3[i];
  for (i=7; i<10; i++)  BC[i] = d3[i];            // cosines are from midpoint
  return ibc;
}

// NOT IMPLEMENTED FOR 2D
int EssentialBC::evalFace(pFace face, double *BC, int p) 
{
  // find centroid of face
  double xyz[3][3];
  F_coord(face,xyz);
  
  double xyzCen[3]={0.333333*(xyz[0][0]+xyz[1][0]+xyz[2][0]),
                    0.333333*(xyz[0][1]+xyz[1][1]+xyz[2][1]),
                    0.333333*(xyz[0][2]+xyz[1][2]+xyz[2][2])};
  double d1[11]={0.0};
  int i, ibc;

  this->pteval(xyzCen, d1, &ibc);

  BC[0] = BC[1] = BC[2] = BC[6] = BC[10] = 0.0;

  for (i=3; i<6; i++)  BC[i] = d1[i];
  for (i=7; i<10; i++)  BC[i] = d1[i];            // cosines are from centroid
  return ibc;
}

void EssentialBC::update_velo(const double *x, double *BC, int *ibcp) // private
{
  int n;
  double vecMagn;

// C1 is comp1 and C3 is comp3
  if ((static_cast<bool>(n=C1) && isAttSet(C1)) 
      || (static_cast<bool>(n=C3) && isAttSet(C3)))  {

      pAttribute adv =Attribute_childByType(AttList[n],"vector direction"); 
      for( int index=0; index< 3; index++ ) 
          BC[3+index] = AttributeTensor1_evalDS((pAttributeTensor1)adv, index, const_cast<double*>(x) );

    double mag = sqrt(SQ(BC[3]) + SQ(BC[4]) + SQ(BC[5]));
    BC[3] /= mag;
    BC[4] /= mag;
    BC[5] /= mag;
    vecMagn = AttributeDouble_evalDS ( (pAttributeDouble)Attribute_childByType( AttList[n], "vector magnitude"),
                                      const_cast<double*>(x));
    if (info->nsd == 2) BC[5] = vecMagn; 
    else BC[6] = vecMagn;	      
    
    if (n==C3)  for (i=3; i<6; i++)  *ibcp = setbit (*ibcp, i);
    else if (n==C1)  {
      int m;

      m = (ABS(BC[3]) > ABS(BC[4])) ? 3 : 4;
      //get the index of the biggest V1DC

      if (info->nsd == 3) m = (ABS(BC[m]) > ABS(BC[5])) ? m : 5;
      *ibcp = setbit (*ibcp, m);           // set the max of these anyhow
    }
  }
}

void EssentialBC::update_thermo(const double *x, double *BC, int *ibcp)
{
  if (isAttSet(D))  BC[0] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[D], const_cast<double*>(x));
  if (isAttSet(T))  BC[1] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[T], const_cast<double*>(x));
  if (isAttSet(P))  BC[2] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[P], const_cast<double*>(x));
  for (int i=0; i<3; i++)  if (isAttSet(i))  *ibcp = setbit (*ibcp, i);
}

//only works in 3D
void EssentialBC::update_axisym(const double *x, double *BC, int *ibcp)
{
  if (isZA())  {
    BC[3] = 1.0;
    BC[4] = 0.0;
    BC[5] = 0.0;
    BC[6] = 0.0;

    BC[7] = 0.0;
    BC[8] = 1.0;
    BC[9] = 0.0;
    BC[10]= 0.0;
  }
  *ibcp = setbit(*ibcp,3);
  *ibcp = setbit(*ibcp,4);
}

void EssentialBC::update_scalar(const double *x, double *BC, int *ibcp)
{
  // treat scalar variables like thermo 
  /*  Came up with a way two check the right number of attributes
      based on the value of NUMVARS  gutsy move for the first C++ code */
 
  int I3nsd, j;
  if (info->nsd == 3) {
     I3nsd = 1;
  } else {
     I3nsd = 0;
  }

  for(i=0; i<nscalar; i++){
    if (isAttSet(S1+i))
       BC[5*I3nsd+7+i] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[S1+i], const_cast<double*>(x));
  }


  for ( i=0; i<nscalar; i++)  
    if (isAttSet(S1+i))  *ibcp = setbit (*ibcp, i+6);     
}

void EssentialBC::update_bellf(const double *x, double *BC, int *ibcp)
{
  // treat scalar variables like thermo 
  /*  Came up with a way two check the right number of attributes
      based on the value of NUMVARS  gutsy move for the first C++ code */
 
  for(i=0; i<6; i++){
    if (isAttSet(VSC+i))
       BC[7+i] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[VSC+i], const_cast<double*>(x));
  }


  for ( i=0; i<4; i++)  
    if (isAttSet(VSC+i))  *ibcp = setbit (*ibcp, i+6);     
  if (isAttSet(XXD))  *ibcp = setbit (*ibcp, 12);
  if (isAttSet(YYD))  *ibcp = setbit (*ibcp, 13);
}

void EssentialBC::pteval(const double *x, double *BC, int *ibcp)
{
  int i, tbit, vbit, sbit,abit, cbit;  // whether thermo, velo or scalar is supposed to be inherited

  *ibcp = 0;

 if (info->nsd == 3){  //3D case
  if (this->gtype == Gface)  {             // this is a face essentialBC
    update_thermo (x, BC, ibcp);           // put thermo qtys. in BC
    update_velo (x, BC, ibcp);             // put velo qtys. in BC
    update_scalar (x, BC, ibcp);           // put scalar qtys. in BC

    // periodic type boundary conditions
#ifdef SIM
    if (GEN_dataI((pGEntity)gf,"iper")) {
      *ibcp = setbit (*ibcp, 10);
    } else if (GEN_dataI((pGEntity)gf,"sper")) {
      *ibcp = setbit (*ibcp, 11);
    }
#else
    int tmp;
    if (GEN_dataI((pGEntity)gf,"iper",&tmp)) {
      *ibcp = setbit (*ibcp, 10);
    } else if (GEN_dataI((pGEntity)gf,"sper",&tmp)) {
      *ibcp = setbit (*ibcp, 11);
    }    
#endif
  } else if (this->gtype == Gedge)  {        // entity is an edge
    if (tbit = getbit (dontinherit, thermo))  update_thermo (x, BC, ibcp);
    if (vbit = getbit (dontinherit, velo))    update_velo (x, BC, ibcp);
    if (sbit = getbit (dontinherit, scalar))  update_scalar (x, BC, ibcp);
    if (abit = getbit (dontinherit, axisym))  update_axisym (x, BC, ibcp);

    if (!tbit || !vbit || !sbit || !abit) { // inherit thermo or velo or scalar
      void* flast=0;
      pGFace gfi;
      pPList e_faces = GE_faces(ge);
      double ValueOnFace[16]; 
      int nf = 0, ibcf, iposBC = 0, nomore = false, none = true;

      while ((gfi = (pGFace)PList_next(e_faces,&flast)) && nf<2 )  {  
                                      // cant have >2 bdry faces per bdry edge
        EssentialBC bcf(gfi);         // go thru rigmarole -- inefficient

        for (i=0; i<ensa_dof+7; i++)  ValueOnFace[i] = 0.0;   // zero them all

        if (bcf.isSet())  {            // "good" bdry face.. user has set BC's
          bcf.pteval(x, ValueOnFace, &ibcf);  
                                   // get the BC values for face in ValueOnFace

          if (!tbit)  {              // thermo inherit
            for (i=0; i<3; i++)  {   // ?? what if D on 1 and P on other etc?
              if (getbit(ibcf, i))  {
                BC[i] = ValueOnFace[i];
                *ibcp = setbit (*ibcp, i);
              }
            }
	  }
          if (!sbit)  {              // scalar inherit
            for (i=6; i<ensa_dof+1; i++)  { 
              if (getbit(ibcf, i))  {
                BC[6+i] = ValueOnFace[6+i];
                *ibcp = setbit (*ibcp, i);
              }
            }
          }
          if (!vbit && !nomore)  {                        // inherit velo 

            if (getbit(ibcf, 3) && getbit(ibcf, 4) && getbit(ibcf, 5))  {
                                                  // if face has C3 set
              *ibcp = setbit (*ibcp, 3);
              *ibcp = setbit (*ibcp, 4);          // inherit all of it 
              *ibcp = setbit (*ibcp, 5); 

              for (i=3; i<7; i++)   BC[i] = ValueOnFace[i];
              nomore = true;                      // dont do velo updates again
              none = false;
            }
            else if (getbit(ibcf, 3) || getbit(ibcf, 4) || getbit(ibcf, 5))  {
                                                  // if face has C1 set

              // foll. takes care that nothing bad will happen if 2 faces have
              // same direction cosines -- 2nd mag will overwrite first
              // doesnt check for valid dir. cosines i.e. sum of sqs. = 1
              // which shudnt be a problem
 
              iposBC += (!CLOSE(BC[3], ValueOnFace[3]) || !CLOSE(BC[4], ValueOnFace[4]) 
                         || !CLOSE(BC[5], ValueOnFace[5]));

              for (i=0; i<4; i++)  BC[4*iposBC+i-1] = ValueOnFace[i+3];
              none = false;
            }
          }
          nf++;
        }
      }
      PList_delete(e_faces);
      // for (i=0; i<3; i++)  if (!CLOSE(BC[i],0.0))  *ibcp = setbit (*ibcp, i);

      if (!nomore && !none)  {           // not all X1, X2, X3 set for inherit
        int m, n;                        // also, at least one set

        n = (ABS(BC[3]) > ABS(BC[4])) ? 3 : 4;     // get the index of the biggest V1DC
        n = (ABS(BC[n]) > ABS(BC[5])) ? n : 5;
        *ibcp = setbit (*ibcp, n);       // set the max of these anyhow
  
        if (!CLOSE(BC[7],0.0) || !CLOSE(BC[8],0.0) || !CLOSE(BC[9],0.0))  {
  
                                         // 2 velocities need to be set
          m = (ABS(BC[7]) > ABS(BC[8])) ? 7 : 8;
          m = (ABS(BC[m]) > ABS(BC[9])) ? m-4 : 5; // get the index of biggest V2DC -4
  
          if (m == n)  {                 // both with same index
            int p, q;
  
            // find the max indices of the other BC's: p for V2DC and q for V1DC
  
            if (m==3)  p = (ABS(BC[8]) > ABS(BC[9])) ? 4 : 5;
            else if (m==4)  p = (ABS(BC[9]) > ABS(BC[7])) ? 5 : 3;
            else if (m==5)  p = (ABS(BC[7]) > ABS(BC[8])) ? 3 : 4;
            
            if (n==3)  q = (ABS(BC[4]) > ABS(BC[5])) ? 4 : 5;
            else if (n==4)  q = (ABS(BC[5]) > ABS(BC[3])) ? 5 : 3;
            else if (n==5)  q = (ABS(BC[3]) > ABS(BC[4])) ? 3 : 4;
            
            m = (ABS(BC[p+4]) > ABS(BC[q])) ? p : q; // put m to be whichever is bigger
          }
          *ibcp = setbit (*ibcp, m);     // set the mth bit - surely m!=n here
        }
      }
    }
#ifdef SIM
    if( GEN_dataI((pGEntity)ge,"iper")) {
      *ibcp = setbit (*ibcp, 10);
    }
    if( GEN_dataI((pGEntity)ge,"sper")) {
      *ibcp = setbit (*ibcp, 11);
    }
#else
    int tmp;
    if( GEN_dataI((pGEntity)ge,"iper",&tmp)) {
      *ibcp = setbit (*ibcp, 10);
    }
    if( GEN_dataI((pGEntity)ge,"sper",&tmp)) {
      *ibcp = setbit (*ibcp, 11);
    }
#endif
  } else if (this->gtype == Gvertex)  {
    if (tbit = getbit (dontinherit, thermo))  update_thermo (x, BC, ibcp);
    if (vbit = getbit (dontinherit, velo))    update_velo (x, BC, ibcp);
    if (sbit = getbit (dontinherit, scalar))  update_scalar (x, BC, ibcp);
    if (abit = getbit (dontinherit, axisym))  update_axisym (x, BC, ibcp);

    if (!tbit || !vbit || !sbit || !abit)  {                  // inherit thermo or velo
      void* elast = 0;
      pPList v_edges = GV_edges(gv);
      pGEdge gei;
      double ValueOnFace[16];  
      int ibce, iposBC = 0, nomore = false;

      while (gei = (pGEdge)PList_next(v_edges,&elast)) {    // ?? >3 edges per vert?

        for (i=0; i<ensa_dof+7; i++)  ValueOnFace[i] = 0.0;   // zero them all

        EssentialBC bce(gei);         // go thru rigmarole -- inefficient
        bce.pteval(x, ValueOnFace, &ibce);    // get the BC values for edge in ValueOnFace

        if (bce.isSet() && ibce!=0)  {   // "good" bdry edge.. 

          if (!tbit)  {                  // thermo inherit
            for (i=0; i<3; i++)  {       // ?? what if D on 1 and P on other ?
              if (getbit(ibce, i))  {    // overwrite earlier if new ones exist
                BC[i] = ValueOnFace[i];
                *ibcp = setbit (*ibcp, i);
              }
            }
          }
          if (!sbit)  {              // scalar inherit
            for (i=6; i<ensa_dof+1; i++)  {   
              if (getbit(ibce, i))  {
                BC[6+i] = ValueOnFace[6+i];
                *ibcp = setbit (*ibcp, i);
              }
            }
          }
          if (!vbit && !nomore)  {                     // inherit velo
            int uset = getbit(ibce, 3);                // ibc set for u velo
            int vset = getbit(ibce, 4);                // ibc set for v velo
            int wset = getbit(ibce, 5);                // ibc set for w velo
            int m, n, p, q;

            if (uset && vset && wset)  {               // all 3 set
              *ibcp = setbit(*ibcp, 3);
              *ibcp = setbit(*ibcp, 4);                // inherit all of it
              *ibcp = setbit(*ibcp, 5);

              for (i=3; i<7; i++)  BC[i] = ValueOnFace[i];
              nomore = true;                           // dont do velos again
            }
            else if ((uset && vset && static_cast<bool>(m=3) && (n=4) && (p=5)) ||
                     (vset && wset && static_cast<bool>(m=4) && (n=5) && (p=3)) ||
                     (wset && uset && (m=5) && (n=3) && (p=4)))  {

                                          // 2 out of 3 of X1V,X2V,X3V set
              *ibcp = setbit(*ibcp, m);
              *ibcp = setbit(*ibcp, n);

              if (getbit(*ibcp,p))  {     // ?? wont check V2 BC's if all 3 set
                double v[3];
                
                for (i=0; i<3; i++)  
                   v[i] = ValueOnFace[i+3]*ValueOnFace[6] + ValueOnFace[i+7]*ValueOnFace[10];

                // the if-else checks if 3 different velocities are set.
                // if so, finds magnitude and resultant 

                if ((!CLOSE(BC[3],ValueOnFace[3]) || !CLOSE(BC[4],ValueOnFace[4]) || 
                     !CLOSE(BC[5],ValueOnFace[5])) && (!CLOSE(BC[3],ValueOnFace[7]) || 
                     !CLOSE(BC[4],ValueOnFace[8]) || !CLOSE(BC[5],ValueOnFace[9])))  {
                  for (i=0; i<3; i++)  v[i] += BC[i+3]*BC[6];
                }
                else if ((!CLOSE(BC[7],ValueOnFace[3]) || !CLOSE(BC[8],ValueOnFace[4]) || 
                     !CLOSE(BC[9],ValueOnFace[5])) && (!CLOSE(BC[7],ValueOnFace[7]) || 
                     !CLOSE(BC[8],ValueOnFace[8]) || !CLOSE(BC[9],ValueOnFace[9])))  {
                  for (i=0; i<3; i++)  v[i] += BC[i+7]*BC[10];
                }
                  
                BC[6] = sqrt(SQ(v[0]) + SQ(v[1]) + SQ(v[2]));
                for (i=0; i<3; i++)  BC[i+3] = v[i]/BC[6];
                nomore = true;            // dont do velocities again
              }  
              else  {                     // V1 and V2 to be directly inherited
                for (i=3; i<11; i++)  BC[i] = ValueOnFace[i];
              }
            }  
            else if ((uset && (m=4) && (n=5) && (p=3)) || 
                     (vset && (m=5) && (n=3) && (p=4)) || 
                     (wset && (m=3) && (n=4) && (p=5)))  {
                                              // 1 out of 3 of X1V,X2V,X3V set

              if (getbit(*ibcp,m) && getbit(*ibcp,n))  {   // other 2 are set
                double v[3];
                
                for (i=0; i<3; i++)  v[i] = BC[i+3]*BC[6] + BC[i+7]*BC[10];

                /* getting magnitude and direction of resultant */

                if ((!CLOSE(ValueOnFace[3],BC[3]) || !CLOSE(ValueOnFace[4],BC[4]) || 
                     !CLOSE(ValueOnFace[5],BC[5])) && (!CLOSE(ValueOnFace[3],BC[7]) || 
                     !CLOSE(ValueOnFace[4],BC[8]) || !CLOSE(ValueOnFace[5],BC[9])))  {
                  for (i=0; i<3; i++)  v[i] += ValueOnFace[i+3]*ValueOnFace[6];
                }
                else if ((!CLOSE(ValueOnFace[7],BC[3]) || !CLOSE(ValueOnFace[8],BC[4]) || 
                     !CLOSE(ValueOnFace[9],BC[5])) && (!CLOSE(ValueOnFace[7],BC[7]) || 
                     !CLOSE(ValueOnFace[8],BC[8]) || !CLOSE(ValueOnFace[9],BC[9])))  {
                  for (i=0; i<3; i++)  v[i] += ValueOnFace[i+7]*ValueOnFace[10];
                }
                BC[6] = sqrt(SQ(v[0]) + SQ(v[1]) + SQ(v[2]));
                for (i=0; i<3; i++)  BC[i+3] = v[i]/BC[6];

                *ibcp = setbit(*ibcp, p);
                nomore = true;            // dont do velocities again
              }  
              else if ((getbit(*ibcp,n) && getbit(*ibcp,p) && (q=m)) || 
                       (getbit(*ibcp,p) && getbit(*ibcp,m) && (q=n)))  {

                if (!(CLOSE(BC[3],ValueOnFace[3]) && CLOSE(BC[4],ValueOnFace[4]) && 
                      CLOSE(BC[5],ValueOnFace[5])) && !(CLOSE(BC[7],ValueOnFace[3]) &&
                      CLOSE(BC[8],ValueOnFace[4]) && CLOSE(BC[9],ValueOnFace[5])))  {

                                // if new V1 is not equal to BC's V1 or V2
                  *ibcp = setbit(*ibcp,q);
                  nomore = true;
                }
              }
              else  {
                if (CLOSE(BC[7],0.0) && CLOSE(BC[8],0.0) && CLOSE(BC[9],0.0))  {
                  *ibcp = setbit(*ibcp, p);
                  iposBC += (!CLOSE(BC[3],ValueOnFace[3]) || !CLOSE(BC[4],ValueOnFace[4])
                             || !CLOSE(BC[5],ValueOnFace[5]));
                                 // ought to be valid DC's i.e. sum of sqs. = 1

                  for (i=0; i<4; i++)  BC[4*iposBC+i-1] = ValueOnFace[i+3];
                }
              }
            }         // end of loop checking for 1 out of 3 set
          }           // end of loop for velo updates
        }             // end of loop for whether parent edge is "good"
      }               // end of loop for going over edges of current vertex
      PList_delete(v_edges);
    }     
            // end of loop for whether inheritance is needed
#ifdef SIM
    if (GEN_dataI((pGEntity)gv,"iper")) {
      *ibcp = setbit (*ibcp, 10);
    } 
    if (GEN_dataI((pGEntity)gv,"sper")) {
      *ibcp = setbit (*ibcp, 11);
    }
#else
    int tmp;
    if (GEN_dataI((pGEntity)gv,"iper",&tmp)) {
      *ibcp = setbit (*ibcp, 10);
    } 
    if (GEN_dataI((pGEntity)gv,"sper",&tmp)) {
      *ibcp = setbit (*ibcp, 11);
    }
#endif  
  }                   // for gtype = Gvertex
  
  }else{  //2D case

  if (this->gtype == Gedge)  {             // this is an edge essentialBC
    update_thermo (x, BC, ibcp);           // put thermo qtys. in BC
    update_velo (x, BC, ibcp);             // put velo qtys. in BC
    update_scalar (x, BC, ibcp);           // put scalar qtys. in BC
    update_bellf (x, BC, ibcp);           // put scalar qtys. in BC

    // periodic type boundary conditions
#ifdef SIM
    if (GEN_dataI((pGEntity)ge,"iper")) {
      *ibcp = setbit (*ibcp, 10);
    } else if (GEN_dataI((pGEntity)ge,"sper")) {
      *ibcp = setbit (*ibcp, 11);
    }
#else
    int tmp;
    if (GEN_dataI((pGEntity)ge,"iper",&tmp)) {
      *ibcp = setbit (*ibcp, 10);
    } else if (GEN_dataI((pGEntity)ge,"sper",&tmp)) {
      *ibcp = setbit (*ibcp, 11);
    }
#endif
  } else if (this->gtype == Gvertex)  {
    if (tbit = getbit (dontinherit, thermo))  update_thermo (x, BC, ibcp);
    if (vbit = getbit (dontinherit, velo))    update_velo (x, BC, ibcp);
    if (sbit = getbit (dontinherit, scalar))  update_scalar (x, BC, ibcp);
    if (cbit = getbit (dontinherit, bellf))  update_bellf (x, BC, ibcp);
    if (abit = getbit (dontinherit, axisym))  update_axisym (x, BC, ibcp);

    if (!tbit || !vbit || !sbit || !abit || !cbit)  {                  // inherit thermo or velo
      void* elast = 0;
      pPList v_edges = GV_edges(gv);
      pGEdge gei;
      double ValueOnFace[16];  
      int ibce, iposBC = 0, nomore = false;

      while (gei = (pGEdge)PList_next(v_edges,&elast)) {    // ?? >3 edges per vert?

        for (i=0; i<ensa_dof+3; i++)  ValueOnFace[i] = 0.0;   // zero them all

        EssentialBC bce(gei);         // go thru rigmarole -- inefficient
        bce.pteval(x, ValueOnFace, &ibce);    // get the BC values for edge in ValueOnFace

        if (bce.isSet() && ibce!=0)  {   // "good" bdry edge.. 

          if (!tbit)  {                  // thermo inherit
            for (i=0; i<3; i++)  {       // ?? what if D on 1 and P on other ?
              if (getbit(ibce, i))  {    // overwrite earlier if new ones exist
                BC[i] = ValueOnFace[i];
                *ibcp = setbit (*ibcp, i);
              }
            }
          }
          if (!sbit)  {              // scalar inherit
            for (i=6; i<ensa_dof+2; i++)  {   
              if (getbit(ibce, i))  {
                BC[1+i] = ValueOnFace[1+i];
                *ibcp = setbit (*ibcp, i);
              }
            }
          }
          if (!cbit)  {              // c1 quintic scalar inherit
            for (i=6; i<10; i++)  {   
              if (getbit(ibce, i))  {
                BC[1+i] = ValueOnFace[1+i];
                *ibcp = setbit (*ibcp, i);
              }
            }
            for (i=11; i<13; i++)  {   
              if (getbit(ibce, i+1))  {
                BC[i] = ValueOnFace[i];
                *ibcp = setbit (*ibcp, i+1);
              }
            }
          }
	  
	  // if 2 comp1 are wanted on one model vertex (2 adjancent edges
	  // have different comp1) set a comp3 on that vertex
	  
          if (!vbit && !nomore)  {                     // inherit velo
            int uset = getbit(ibce, 3);                // ibc set for u velo
            int vset = getbit(ibce, 4);                // ibc set for v velo
            int m, n, p, q;

            if (uset && vset)  {               // both set
              *ibcp = setbit(*ibcp, 3);
              *ibcp = setbit(*ibcp, 4);                // inherit all of it
            } else if (uset)  {
              *ibcp = setbit(*ibcp, 3);
            } else if (vset)  {
              *ibcp = setbit(*ibcp, 4);
	    }         // end of loop checking for 1 out of 3 set

            for (i=3; i<7; i++)  BC[i] = ValueOnFace[i];
            nomore = true;                           // dont do velos again
          }           // end of loop for velo updates
        }             // end of loop for whether parent edge is "good"
      }               // end of loop for going over edges of current vertex
      PList_delete(v_edges);
    }                 // end of loop for whether inheritance is needed

#ifdef SIM
    if (GEN_dataI((pGEntity)gv,"iper")) {
      *ibcp = setbit (*ibcp, 10);
    } 
    if (GEN_dataI((pGEntity)gv,"sper")) {
      *ibcp = setbit (*ibcp, 11);
    }
#else
    int tmp;
    if (GEN_dataI((pGEntity)gv,"iper",&tmp)) {
      *ibcp = setbit (*ibcp, 10);
    } 
    if (GEN_dataI((pGEntity)gv,"sper",&tmp)) {
      *ibcp = setbit (*ibcp, 11);
    }
#endif
  }                   // for gtype = Gvertex for 2D

  }
}

// not implemented for c1 quintic
void EssentialBC::takeBCfromIC(double *BC, double *qTot, 
                               int GDOFnum)
{
  /* qTot here is [PPPPPPuuuuuuvvvvvvwwwwwwTTTTTTssssss]*/
  int nshg = nshgTot;
  double cosines[3];
  double mag;

  if (info->nsd == 3){  //3D meshes
    cosines[0]=qTot[1*nshgTot+GDOFnum];
    cosines[1]=qTot[2*nshgTot+GDOFnum];
    cosines[2]=qTot[3*nshgTot+GDOFnum];
    mag=sqrt(cosines[0]*cosines[0]
           +cosines[1]*cosines[1]
           +cosines[2]*cosines[2]);
    if(mag==0){
      BC[3]=1;
      BC[4]=0;
      BC[5]=0;
      BC[6]=0;
    }else{
      BC[3] = cosines[0]/mag;
      BC[4] = cosines[1]/mag;
      BC[5] = cosines[2]/mag;
      BC[6] = mag;
    }

    if (isAttSet(D))  BC[0] = qTot[GDOFnum];
    if (isAttSet(T))  BC[1] = qTot[4*nshgTot+GDOFnum];
    if (isAttSet(P))  BC[2] = qTot[GDOFnum];

    for(i=0; i<nscalar; i++){
      if (isAttSet(S1+i))  
        BC[12+i] = qTot[(5+i)*nshgTot+GDOFnum];
    }
  }else{   //2D meshes  
    cosines[0]=qTot[1*nshgTot+GDOFnum];
    cosines[1]=qTot[2*nshgTot+GDOFnum];
    mag=sqrt(cosines[0]*cosines[0]
           +cosines[1]*cosines[1]);
    if(mag==0){
      BC[3]=1;
      BC[4]=0;
      BC[5]=0;
    }else{
      BC[3] = cosines[0]/mag;
      BC[4] = cosines[1]/mag;
      BC[5] = mag;
    }

    if (isAttSet(D))  BC[0] = qTot[GDOFnum];
    if (isAttSet(T))  BC[1] = qTot[3*nshgTot+GDOFnum];
    if (isAttSet(P))  BC[2] = qTot[GDOFnum];

    for(i=0; i<nscalar; i++){
      if (isAttSet(S1+i))  
        BC[7+i] = qTot[(4+i)*nshgTot+GDOFnum];
    }
  }  
}
