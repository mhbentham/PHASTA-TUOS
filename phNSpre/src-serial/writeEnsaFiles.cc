/**********************************************************************/
/* This function collects the information from the mesh database into */
/* the arrays which will be used by ENSA.                             */
/* 								      */
/* Modified for 2D meshes by Elaine Bohr Dec 2005          	      */
/**********************************************************************/
#include <iostream>
#include <cstdlib>
#include "nspre_data.h"
#include "nspre_functions.h"
#include "EnsaParameters.h"

extern "C" void memprof( char* location );
extern time_t tstart[MAXNT];
extern double eltime[MAXNT];
extern int ICtyp,intBC,fType,zScale,ensa_dof,rStart, lStart,rRead;
int lstep=0;
extern globalInfo* info;
extern bool c1quintic;
extern bool FaceConnectivity;
extern int nscalar;

int nshgTot=1;

void 
writeEnsaFiles( pGModel model, 
                pMesh mesh,
                void* bdry ) {
                

    int I3nsd, numEBC, numNBC;
    
    nshgTot = info->nshg;

    if (info->nsd == 3) {
       I3nsd = 1;
    } else {
       I3nsd = 0;
    }
    numEBC = 5*I3nsd + nscalar + 7;
    numNBC = nscalar+6;
    if (c1quintic == true) {
       numEBC = 13;
       numNBC = 12;
    }
    
    EnsaParameters *ePar = new EnsaParameters( info->nsd, info->nenmax,
                                               0, ensa_dof, numEBC,
                                               numNBC, fType);

    /* These parameters  changed with the addition of scalar equations. 
       The formulas should key off of numVars 
       (I see that ndof is used elsewhere in the code so...). 
       numVars gets set to ensa_dof by an input variable from the command line
       numEBC  = 12   changed to ensa_dof + 7 ( 12 is for theta)
       numNBC  = 6    changed to ensa_dof + 1 ( 1 extra for turbulence wall) 
       For 2D needed to change
       numEBC = 5*I3nsd + nscalar + 7 where I3nsd =0 for 2D and 1 for 3D
       numNBC = nscalar+6 so to keep same numbering in BCB in 2D and 3D */

    double zscale[3]={ 0.0, 0.0, 0.0 };
    if ( zScale ) {
        cout << "\nEnter the scale factor for the x,y,z coordinate: " ;
        cin >> zscale[0];
        cin >> zscale[1];
        if (info->nsd == 3) cin >> zscale[2];
    }
    int nv = ensa_dof;
    int numnp= M_numVertices( mesh );
    int nshgReqst=nshgTot; // read all modes by default
    double* qTot;

    if( lStart ){
        // for now we use this to get a quadratic solution
        // load in the auxilliary linear mesh

        /* pMesh lmesh=MS_newMesh(model);
           M_load(lmesh,"geoml.sms");
           int numnpl=M_numVertices(lmesh);
           double* qtmpl = (double *)calloc(nv*numnpl,sizeof(double));
           restart(qtmpl,numnpl,rRead,&lstep,"restartc.inp");
           nshg = numnp + M_numEdges(mesh);
           qtmp = (double *)calloc(nv*nshg,sizeof(double));
           lin2quad (mesh, lmesh, qtmpl, qtmp,  nshg);*/
    } else {
        if(rStart){
            // use linear portion
            if (rStart == 3) {
                nshgReqst = numnp;
            } else if (rStart == 4) {
                // use quadratic portion
                nshgReqst = numnp + M_numEdges(mesh);
            }
      
            qTot = (double *)calloc(nv*nshgReqst,sizeof(double));
            for (int i=0; i< nv*(nshgTot); i++) qTot[i]=0;
            //
            //  We have zeroed the qTot array because restart is going to fill
            //  only the portion of this array that we request it to based on
            //  our choice of nshgReqst and rRead (below called nvReqst for
            //  consistency. Note that the array is shaped nv*nshgTot because
            //  later, higher order modes will need to index into this TOTAL
            //  solution array for all variables.  If they were not read it
            //  will find zero's there which we assume will be ok or replaced
            //  by attribute expressions. 
            //
            restart( qTot, nshgReqst, rRead, &lstep,"restartc.inp" );
        } else {
            qTot = ( double* )calloc(nv*nshgReqst,sizeof(double));
        }
    }  
    EnsaArrays* eArr = new EnsaArrays;

    AllocateEnsaArrays( eArr, ePar );

    getX( mesh, eArr->x );

    attachPeriodicBC( mesh, eArr->iper );

    info->numpbc = attachEssentialBC( model, mesh, qTot, 
                                           eArr->nBC, eArr->iBC, eArr->BC  );
                                           
    attachNaturalBC( model, mesh, eArr->iBCB, eArr->BCB, ePar->getNUMNBC() );
      
    attachInitialCondition( model, mesh, qTot, eArr->q );

    if (info->nsd == 3)  {	//3D boundary is a face
    getConnectivity( mesh, 
                     (vector<pFace>*)bdry, 
                     eArr->ien, 
                     eArr->ien_sms,
                     eArr->ienb, 
                     eArr->ienb_sms 

                     );
    } else {			//2D boundary is an edge
    getConnectivity2D( mesh, 
                     (vector<pEdge>*)bdry, 
                     eArr->ien, 
                     eArr->ien_sms,
                     eArr->ienb, 
                     eArr->ienb_sms 

                     );
    }
    if(FaceConnectivity){

        getFaceConnectivity( mesh, 
                             *((vector<pFace>*)bdry), 
                             eArr->ief, 
                             eArr->iefb 
            );

//                       eArr->ief_sms, 
//                       eArr->iefb_sms, 
        

    }//if(FaceConnectivity) 

      
    writeEnsaArrays( eArr, ePar, zscale );

    DeAllocateEnsaArrays( eArr );
         
    delete eArr;
    free( qTot );
}
