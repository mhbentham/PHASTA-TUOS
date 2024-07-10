// For the Vertex and Edge versions of getMaster I should remove the
// else and make, the theta case generic so that it works for theta =
// 0

//
// return the periodic master
//
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include "nspre_data.h"
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#include "MeshAdjTools.h"
#include "mesh_interface.h"
#endif

#define ABS(x)      ((x) < 0 ? -(x) : (x))

extern ofstream perr;
using namespace std;
void fixFaces(pFace face, pFace mface,double theta);

const double TOL = 1.e-5;        // for matching master/slave pairs
const double ATOL = 0.001;

extern int C_raneql(double real1,double real2,double tolran);



// find the vertex in vIter that is the periodic partner of xyz
pVertex 
getMasterV(VIter mvIter, pVertex vert,
          double theta, int trans ) {

    double xyzM[3],xyzS[3],xyzR[3];

    V_coord(vert,xyzS);

    // translate the slave coordinate by theta

    xyzR[0] = cos(theta) * xyzS[0] - sin(theta) * xyzS[1];
    xyzR[1] = sin(theta) * xyzS[0] + cos(theta) * xyzS[1];
    xyzR[2] = xyzS[2];

    pVertex mvtx;
    int neq;

    // rotated case theta != zero
    pVertex minvtx;
    double mindist=1.0e16;
    double dist;
    if (theta != 0.0) {
        while (mvtx = VIter_next(mvIter)){
            V_coord(mvtx,xyzM);
      
            // find the number of equal coordinates (x,y,z) 
            neq = (C_raneql(xyzR[0],xyzM[0],TOL))+
                (C_raneql(xyzR[1],xyzM[1],TOL))+
                (C_raneql(xyzR[2],xyzM[2],TOL));
            if (neq==(4-trans)) return mvtx;
            //
            //Below is code added to catch points that don't get caught above
            //
            dist= (xyzR[0]-xyzM[0])*(xyzR[0]-xyzM[0])+
                (xyzR[1]-xyzM[1])*(xyzR[1]-xyzM[1])+
                (xyzR[2]-xyzM[2])*(xyzR[2]-xyzM[2]);

            if(dist<mindist){
                mindist=dist;
                minvtx = mvtx;  
            }
        }
        if(mindist>TOL*TOL){

            cerr << sqrt(mindist)<< "\n";
      
        }
        return minvtx;
    } else {

        while (mvtx = VIter_next(mvIter)){
            V_coord(mvtx,xyzM);
      
            // find the number of equal coordinates (x,y,z) 
            neq = (C_raneql(xyzR[0],xyzM[0],TOL))+
                (C_raneql(xyzR[1],xyzM[1],TOL))+
                (C_raneql(xyzR[2],xyzM[2],TOL));
      
            if (neq==(3-trans)) return mvtx;
        }

    }

    cerr << "\n Error in getMaster: periodic master vertex not found...\n";
    return mvtx;

}
    
// find the edge in meIter that is the periodic partner of xyz
// checks midpoint of edge
pEdge 
getMasterE(EIter meIter, pEdge edg,
          double theta, int trans) {

    double xyzM[2][3],xyzS[2][3];
    double xcS[3],xcM[3];
    double xcSr[3];

    V_coord(E_vertex(edg,0),xyzS[0]);
    V_coord(E_vertex(edg,1),xyzS[1]);

    xcS[0] = 0.5*(xyzS[0][0]+xyzS[1][0]);
    xcS[1] = 0.5*(xyzS[0][1]+xyzS[1][1]);
    xcS[2] = 0.5*(xyzS[0][2]+xyzS[1][2]);

    // rotate the centroid of the slave by theta so that it is on the
    // master geometric face
  
    // z-restriction, so far this code assumes that the roration is on
    // the z-axis, it cannot successfully deal with anything else.
  

    xcSr[0] = cos(theta) * xcS[0] - sin(theta) * xcS[1];
    xcSr[1] = sin(theta) * xcS[0] + cos(theta) * xcS[1];
    xcSr[2] = xcS[2];

  
    pEdge medg;
    pEdge minedge;
    double mindist = 1.0e16;

    int neq;

    if (theta != 0.0) {
        while (medg = EIter_next(meIter)){
            V_coord(E_vertex(medg,0),xyzM[0]);
            V_coord(E_vertex(medg,1),xyzM[1]);

            xcM[0] = 0.5*(xyzM[0][0]+xyzM[1][0]);
            xcM[1] = 0.5*(xyzM[0][1]+xyzM[1][1]);
            xcM[2] = 0.5*(xyzM[0][2]+xyzM[1][2]);
         
            // find the number of equal coordinates (x,y,z) 
            neq = (C_raneql(xcSr[0],xcM[0],TOL))+
                (C_raneql(xcSr[1],xcM[1],TOL))+
                (C_raneql(xcSr[2],xcM[2],TOL));
       
            if ( neq == (4-trans)){
                // make sure that the edges are oriented in the same way
                xcS[0] = cos(theta) * xyzS[1][0] - sin(theta) * xyzS[1][1];
                xcS[1] = sin(theta) * xyzS[1][0] + cos(theta) * xyzS[1][1];
                xcS[2] = xyzS[1][2];
            
                neq = (C_raneql(xyzM[1][0], xcS[0], TOL))+
                    (C_raneql(xyzM[1][1], xcS[1], TOL))+
                    (C_raneql(xyzM[1][2], xcS[2], TOL));

                if( neq != 3 ) {
#ifdef DEBUG
                    perr << "Flipping Edge: "<< EN_id( (pEntity)edg ) <<  endl;
#endif 
                }
                return medg;
            } 
            double dist = (xcSr[0]-xcM[0])* (xcSr[0]-xcM[0]) +
                (xcSr[0]-xcM[1])* (xcSr[1]-xcM[1]) +
                (xcSr[0]-xcM[2])* (xcSr[2]-xcM[2]) ;
         
            if ( dist < mindist ) {
                mindist = dist ;
                minedge = medg;
            }
        }
        // Here we trap the edges which do not get a master in the usual
        // we just pick the edge with the smallest centroid distance.
        if( mindist > TOL*TOL ) 
            cerr << " Warning: slave->master distance > tolerance";

        V_coord(E_vertex(minedge,1),xyzM[1]);
        xcS[0] = cos(theta) * xyzS[1][0] - sin(theta) * xyzS[1][1];
        xcS[1] = sin(theta) * xyzS[1][0] + cos(theta) * xyzS[1][1];
        xcS[2] = xyzS[1][2];
            
        neq = (C_raneql(xyzM[1][0], xcS[0], TOL))+
            (C_raneql(xyzM[1][1], xcS[1], TOL))+
            (C_raneql(xyzM[1][2], xcS[2], TOL));

        if( neq != 3 ){
#ifdef DEBUG
            perr << "Flipping Edge: "<< EN_id( (pEntity)edg ) <<  endl;
#endif 
        }
        return medg;
       
    } else {
  
        while (medg = EIter_next(meIter)){

            V_coord(E_vertex(medg,0),xyzM[0]);
            V_coord(E_vertex(medg,1),xyzM[1]);

            xcM[0] = 0.5*(xyzM[0][0]+xyzM[1][0]);
            xcM[1] = 0.5*(xyzM[0][1]+xyzM[1][1]);
            xcM[2] = 0.5*(xyzM[0][2]+xyzM[1][2]);
        
            // find the number of equal coordinates (x,y,z) 
            neq = (C_raneql(xcS[0],xcM[0],TOL))+
                (C_raneql(xcS[1],xcM[1],TOL))+
                (C_raneql(xcS[2],xcM[2],TOL));
        
            if (neq==(3-trans)) {
          
                // reverse slave direction if necessary
                neq = (C_raneql(xyzS[0][0],xyzM[0][0],TOL))+
                    (C_raneql(xyzS[0][1],xyzM[0][1],TOL))+
                    (C_raneql(xyzS[0][2],xyzM[0][2],TOL));
          
                if ((3-trans) != neq){
#ifdef DEBUG
                    perr << "Flipping Edge: "<< EN_id( (pEntity)edg ) <<  endl;
#endif 
                }
                return medg;
            }
        }
    }    
    cerr << "\nError in getMaster: periodic master edge not found...\n";
    exit(-1);
    return 0;
}


// find the face in mfIter that is the periodic partner of xyz
// checks centroid of face
pFace 
getMasterF(FIter mfIter, pFace fac, double theta ) {

    double xyzM[4][3],xyzS[4][3];
    double xcS[3],xcM[3],xcSr[3];

    pPList sverts = F_vertices(fac,1);
    int numVerts = PList_size( sverts );
    for (int i=0; i < numVerts ; i++)  
        V_coord((pVertex)PList_item(sverts,i),xyzS[i]);
    PList_delete(sverts);
    
    for( int dir =0; dir < 3 ; dir++ ) {
        xcS[ dir ] = 0.0;
        for( int v=0; v < numVerts; v++ ) {
            xcS[ dir ] += xyzS[ v ][ dir ];
        }
        xcS[ dir ] /= numVerts;
    }
            
    xcSr[0] = cos(theta)* xcS[0] - sin(theta)* xcS[1];
    xcSr[1] = sin(theta)* xcS[0] + cos(theta)* xcS[1];
    xcSr[2] = xcS[2];

    pFace mfac;
    int neq;
    while (mfac = FIter_next(mfIter)){
      
        pPList mverts = F_vertices(mfac,1);
        int numVertsm = PList_size( mverts ) ;
        if ( numVerts != numVertsm ) continue;
        
        for (int i=0; i < numVerts ; i++)
            V_coord((pVertex)PList_item(mverts,i),xyzM[i]);
        PList_delete(mverts);
      
        for( int dir =0; dir < 3 ; dir++ ) {
            xcM[ dir ] = 0.0;
            for( int v=0; v < numVerts; v++ ) {
                xcM[ dir ] += xyzM[ v ][ dir ];
            }
            xcM[ dir ] /= numVerts;
        }
  
        // find the number of equal coordinates (x,y,z) 
        neq = (C_raneql(xcSr[0],xcM[0],TOL))+
            (C_raneql(xcSr[1],xcM[1],TOL))+
            (C_raneql(xcSr[2],xcM[2],TOL));
        
        if ((theta != 0 && neq==3) || (neq == 2)) {
            
            // since the mesh faces must be periodic, i.e. the same,
            // their edges and normals must be the same


            fixFaces(fac,mfac,theta);
            
            return mfac;
        }
    }
    
    cerr << "\nError in getMaster: periodic master face not found...\n";
    exit(-1);
    return 0;
}

//
// align two faces so their edges and normals are the same
//
void 
fixFaces(pFace face, pFace mface,double theta) {

    pPList svts,mvts;
    pVertex vS,vM;
    double xyzS[3][3],xyzM[3][3];
    double xyzSr[3];
    int i,neq;

    svts=F_vertices(face,1);
    mvts=F_vertices(mface,1);
    for (i=0; i < 3; i++){
        V_coord((pVertex)PList_item(svts,i),xyzS[i]);
        V_coord((pVertex)PList_item(mvts,i),xyzM[i]);
    }

    /* Compute the normals of the slave and master and see if they are
       in the same direction otherwise they have to be flipped this will
       after this we make sure that the vertices match exactly. */
  
    double n[3], m[3], nCm[3];
    F_normalVector(face,1, n);
    F_normalVector(mface,1, m);

    nCm[0] = n[1]*m[2]-m[1]*n[2];
    nCm[1] = n[2]*m[0]-m[2]*n[0];
    nCm[2] = n[0]*m[1]-m[0]*n[1];
 
    double normS = sqrt( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] );
    double normM = sqrt( m[0]*m[0]+m[1]*m[1]+m[2]*m[2] );
    double lthta = asin( nCm[2]/( normS * normM ) );
  
    if(lthta != theta) {
        perr << "need to flip Face for higher order: "<< EN_id( (pEntity)face ) <<  endl;
        
        // This invokes direct modification of the mesh data base.
        // From MeshSim version 6.0 onward this function is no loner available.
        // Future implementations will have to handle this differently.
        // For now, this feature is disabled.
        cout<<"F_chDir(face) for p>2 modes: from MeshSim version 6.0 onward this function is no loner available.\n";
        cout<<"Future implementations will have to handle this differently.\n";
        cout<<" For now, this feature is disabled.\n";
//        F_chDir(face);
        exit(1);
    }

    /* compare the first vertex of the first edge of slave and master */
    vS = E_vertex(F_edge(face,0),0);
    vM = E_vertex(F_edge(mface,0),0);
    V_coord(vS,xyzS[0]);
    V_coord(vM,xyzM[0]);

    /* First rotate it by theta */ 
    xyzSr[0] = xyzS[0][0] * cos(theta) - xyzS[0][1] * sin(theta);
    xyzSr[1] = xyzS[0][0] * sin(theta) + xyzS[0][1] * cos(theta);
    xyzSr[2] = xyzS[0][2];

    /* find the number of equal coordinates (x,y,z) */
    neq = (C_raneql(xyzSr[0],xyzM[0][0],TOL))+
        (C_raneql(xyzSr[1],xyzM[0][1],TOL))+
        (C_raneql(xyzSr[2],xyzM[0][2],TOL));
    if ((theta > ATOL  && neq!=3) || (theta < ATOL  && neq != 2)) {
        //F_permuteEdges(face,2);
        perr <<"No C-interface to function for permuting edges of a face\n";
    }

    vS = E_vertex(F_edge(face,0),0);
    V_coord(vS,xyzS[0]);
    xyzSr[0] = xyzS[0][0] * cos(theta) - xyzS[0][1] * sin(theta);
    xyzSr[1] = xyzS[0][0] * sin(theta) + xyzS[0][1] * cos(theta);
    xyzSr[2] = xyzS[0][2];
    /* find the number of equal coordinates (x,y,z) */
    neq = (C_raneql(xyzSr[0],xyzM[0][0],TOL))+
        (C_raneql(xyzSr[1],xyzM[0][1],TOL))+
        (C_raneql(xyzSr[2],xyzM[0][2],TOL));
    if ((theta > ATOL  && neq!=3) || (theta < ATOL  && neq != 2)) {
        //    F_permuteEdges(face,2);
        perr <<"No C-interface to function for permuting edges of a face\n";
    }

    vS = E_vertex(F_edge(face,0),0);
    V_coord(vS,xyzS[0]);
    xyzSr[0] = xyzS[0][0] * cos(theta) - xyzS[0][1] * sin(theta);
    xyzSr[1] = xyzS[0][0] * sin(theta) + xyzS[0][1] * cos(theta);
    xyzSr[2] = xyzS[0][2];
    /* find the number of equal coordinates (x,y,z) */
    neq = (C_raneql(xyzSr[0],xyzM[0][0],TOL))+
        (C_raneql(xyzSr[1],xyzM[0][1],TOL))+
        (C_raneql(xyzSr[2],xyzM[0][2],TOL));
    if ((theta > ATOL  && neq!=3) || (theta < ATOL  && neq != 2)) {
        perr << "\n master-slave edge/face, orientation/order mismatch:";
        perr << "higher order alert\n" << endl;
    }

    PList_delete(svts);
    PList_delete(mvts);
}
