#include <iostream>   
#include <math.h>  
#include "nspre_data.h"
#include "SPEBC.h"  

extern int ensa_dof;  //bring in the # of variables to help find # of scalars

#define numAttsS 1        // number of SPEBC attributes
char strAttS[numAttsS][MLEN] ={"spebc"};
enum { SP };

void 
GF_inflow_normal(pGFace gface, pGModel model, double *xyz) ;

SPEBC::SPEBC(pGFace gface) : BoundaryCondition () {
    if (GEN_attrib((pGEntity)gface, strAttS[SP]))  {
        // if slave then everything set
        this->set = true;
    }
}

SPEBC::SPEBC(pGEdge gedge) : BoundaryCondition () {

    if (GEN_attrib((pGEntity)gedge, strAttS[SP]))  { 
        // if user sets, no inheritance
        this->set = true;  
    }  else  {

        pPList gefaces = GE_faces(gedge);
        void* restart = 0;
        pGFace gfi;

        while (gfi = (pGFace)PList_next(gefaces,&restart))  { 
            // loop over faces
            SPEBC sbc(gfi);            // set bc on the face
            if (sbc.isSet())  {        // should be slave
                this->set = true;
                break;                   // break from loop => edge is slave if
            }                          // at least one connecting face is slave
        }
        PList_delete(gefaces);
    }
}

SPEBC::SPEBC(pGVertex gvert) : BoundaryCondition () {
    /* Exactly similar to edge above so look there for comments */

    if (GEN_attrib((pGEntity)gvert, strAttS[SP]))  {
        this->set = true;  
    }  else  {

        pPList gvfaces = GV_faces(gvert);
        void* restart= 0;
        pGFace gfi;
    
        while (gfi = (pGFace)PList_next(gvfaces,&restart))  {
            SPEBC sbc(gfi);
            if (sbc.isSet())  {
                this->set = true;
                break;
            }
        }
        PList_delete(gvfaces);
    }
}

int SPEBC::isSPEBC(void) 
{
    if(isSet()) return 1;
    else return 0;
}

/*******************************************************************
 * Compute the equation of the virtual plane when using SPEBC
 * and find regions that are cut by the virtual plane
 *
 * Elaine Bohr (Fall 2002)
 *******************************************************************
 */

/* Finding the equation of the virtual plane */

void 
eqn_plane( pGModel model, pMesh mesh, double *coeffs, pGFace gface ) {
    pRegion region;
    pPList verts;
    pGEntity ent;
    pAttribute att;
    int i;
    int tmp=0;
    pVertex vertex;
    int numVerts;
    double xyz[3];
    double nrml[3] ;
    double plandist;
    
    
    /* distance between the inlet & recycle planes from SPEBC attribute */    
    att = GEN_attrib((pGEntity)gface, "spebc");
    plandist = AttributeDouble_value( (pAttributeDouble) att );
    
    
    /* find the normal to the inflow plane */    
    GF_inflow_normal(gface,model,nrml);
    
    
    /* find one node and its coodrinates on inlet plane */
    RIter rIter = M_regionIter(mesh);
    while(region = RIter_next(rIter)) {
        verts = R_vertices(region,1);
        numVerts = PList_size(verts);
        for(i=0;i<numVerts;i++){
            vertex=(pVertex)PList_item(verts,i);
            ent=V_whatIn(vertex);
            if(GF_inClosure(gface,ent)){ 
                tmp=1;
                break;
            }
        }
        if(tmp == 1) break;
    }
    RIter_delete(rIter);
    
    V_coord(vertex,xyz);
    
    
    coeffs[0] = nrml[0];
    coeffs[1] = nrml[1];
    coeffs[2] = nrml[2];
    coeffs[3] = coeffs[0]*xyz[0]+coeffs[1]*xyz[1]+coeffs[2]*xyz[2]+plandist;
    
}

/* Finding if the region is cut by the virtual plane */

int 
isCut(pPList verts, double *coeffs) {
    
    pVertex vertex;
    int numVerts,i,j;
    double xyz[3];
    double* erreur;
    
    
    numVerts = PList_size(verts);
    erreur = new double [numVerts];
    for(i=0;i<numVerts;i++){
        vertex=(pVertex)PList_item(verts,i);
        V_coord(vertex,xyz);
        erreur[i] = coeffs[3] - coeffs[0]*xyz[0] - coeffs[1]*xyz[1] 
                    - coeffs[2]*xyz[2];
    }
    for(i=0;i<numVerts;i++){
        for(j=i;j<numVerts;j++){
            if (erreur[i]*erreur[j] <= 0.0) 
                return 1;
        }
    }
    return 0;
    
}

/* Calculating the inward normal to a model face */

void 
GF_inflow_normal(pGFace gface, pGModel model, double *xyz) {
    double xyztmp[3][3];
    double v[2][3];
    double fourth[3];
    pPList facevertices = GF_vertices(gface);
    pPList regions = GF_regions(gface);
    pGRegion reg = (pGRegion)PList_item(regions,0);
    pPList regionvertices = GR_vertices(reg);
    pGVertex pvtx;
    
    
    /* finding 3 model vertices on the inlet plane */
    for(int i = 0; i<3; i++) {
        pvtx = (pGVertex)PList_item(facevertices,i);
        GV_point(pvtx,xyztmp[i]);
    }
    
    /* Calculating two vectors from those three point's coodinates */
    for(int i = 0; i<2; i++) {
        v[i][0]=xyztmp[1+i][0]-xyztmp[i][0];
        v[i][1]=xyztmp[1+i][1]-xyztmp[i][1];
        v[i][2]=xyztmp[1+i][2]-xyztmp[i][2];
    }
    
    /* now take the cross-product to get a normal vector to both edges */
    xyz[0]=v[0][1]*v[1][2]-v[0][2]*v[1][1];
    xyz[1]=v[0][2]*v[1][0]-v[0][0]*v[1][2];
    xyz[2]=v[0][0]*v[1][1]-v[0][1]*v[1][0];
    
    /* and normalize */
    double mag = xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2];
    mag=sqrt(mag);
    xyz[0]=xyz[0]/mag;
    xyz[1]=xyz[1]/mag;
    xyz[2]=xyz[2]/mag;
    
    /* find one vertex from curent region that is not on inlet plane */
    int numvert = PList_size(regionvertices);
    for(int i = 0; i < numvert; i++) {
        if(PList_item(regionvertices,i) != PList_item(facevertices,i)) {
            pvtx = (pGVertex)PList_item(regionvertices,i);
            GV_point(pvtx,fourth);
            break;
        }
    }
    /* calculate the vector formed by this point and one of the inlet plane
       points */
    fourth[0] = fourth[0] - xyztmp[0][0];
    fourth[1] = fourth[1] - xyztmp[0][1];
    fourth[2] = fourth[2] - xyztmp[0][2];
    
    /* claculate scalar product between the normal and the last vector */
    double scal =  xyz[0]*fourth[0] + xyz[1]*fourth[1] + xyz[2]*fourth[2];
    
    /* if the scalar is negative inverse the normal */
    if (scal < 0.0) {
        xyz[0]= - xyz[0];
        xyz[1]= - xyz[1];
        xyz[2]= - xyz[2];
    }
}
