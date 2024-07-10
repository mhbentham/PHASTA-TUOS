#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "nspre_data.h"
#include "nspre_functions.h"
#include "phastaIO.h"
#include "EnsaParameters.h"
#ifndef WIN32
#include <strings.h>
#else
void  bzero(void* ptr, size_t sz) {
    int i;
    char *cptr;
    cptr = (char*) ptr;
    for (i=0; i < sz; i++) {
        cptr[i]=0;
    }
    return;
}
#endif

extern int getnumBshape(blockKey bKey);
extern int nshgTot;
extern int lstep;
extern int WRITEASC,ensa_dof;
extern int zScale;

// declared in genblock.cc
extern forwardblock Iblock;
extern forwardblock Bblock;
extern blockKey *RIblock;
extern blockKey *RBblock;
extern int* NBblock;
extern int* NIblock;
extern int* Nshape;
extern int* NshapeB;
extern char* oformat;
extern char options[];
extern char version[];
extern pMesh mesh;
extern globalInfo* info;
extern int old_format;

// command line argument
extern bool FaceConnectivity;

// in genBlock.cc 
extern int* numInteriorFacesOfBlock;
extern int* numBoundaryFacesOfBlock;

int
numVertsR( pRegion region ) {
    pPList vertx = R_vertices( region, 1 );
    int vertn = PList_size( vertx );
    PList_delete( vertx );
    return vertn;
}

int
numVertsF( pFace face ) {
    pPList vertx = F_vertices( face, 1 );
    int vertn = PList_size( vertx );
    PList_delete( vertx );
    return vertn;
}

// this function is called by writeEnsaFiles in writeEnsaFiles.cc
// allocates dynamic stack memory for ensa arrays 
// (not to get confused by usage of proc id (pid): not used in  serial
// stage  of NSPre
void
AllocateEnsaArrays( EnsaArrays* e, 
                    EnsaParameters *par) {

    int nsd      = par->getNSD();
    int numVars  = par->getNUMVARS();
    int numEBC   = par->getNUMEBC();
    int numNBC   = par->getNUMNBC();

    int num = M_numVertices( mesh );
    e->x = new double [num*nsd];
  
    num = info->nshg;
    e->q = new double* [num];
    for(int y=0; y< num; y++ ) 
        e->q[y] = new double [numVars];

    int tmpblk = Iblock.size();
    e->ien = new int* [tmpblk];
    e->ien_sms = new int* [tmpblk];
    for(int j=0; j < tmpblk; j++){
        e->ien[j] = new int [NIblock[j]*Nshape[j]];
        e->ien_sms[j] = new int [ NIblock[j] ];
    }
  
    tmpblk = Bblock.size();
    e->ienb = new int*[tmpblk];
    e->ienb_sms = new int* [tmpblk];
    e->BCB = new double**[tmpblk];  /* saving calls to size and unnecessary */
    e->iBCB = new int**[tmpblk];    /* loops by allocating some stuff here */
    for(int j=0; j<tmpblk; j++){
        num = NBblock[j];
        e->ienb[j] = new int [num*NshapeB[j]];
        e->ienb_sms[j] = new int [ num ];
        e->BCB[j] = new double* [num];
        e->iBCB[j] = new int* [num];
        for(int y=0; y< num; y++ ) {
            e->BCB[j][y] = new double [numNBC];
            e->iBCB[j][y] = new int [2];
        }
    }

    if(FaceConnectivity){
        // allocate memory for the interior face connectivity
        tmpblk = Iblock.size();
        e->ief = new int* [tmpblk];
//    e->ief_sms = new int* [tmpblk];
        
        for(int j=0; j < tmpblk; j++){
            
            // NIblock[j] : number of eles in block
            // NEFPerBlock() number of element faces characteristic to the block
            e->ief[j] =  new int [NIblock[j]*NEFPerBlock(RIblock[j])]; // DON't KNOW
            //e->ief_sms[j] = new int [ NIblock[j] ];
            
        }    
        numInteriorFacesOfBlock = new int[tmpblk];   // a one dimensional array
        // storing the total number of faces in the perticular interior block
        
        
        // allocate memory for the boundary face connectivity
        tmpblk = Bblock.size();
        e->iefb = new int*[tmpblk];
//    e->iefb_sms = new int* [tmpblk];
        
        for(int j=0; j<tmpblk; j++){
            num = NBblock[j];
//        e->iefb[j] = new int [numFaces];  wrong !
            e->iefb[j] =  new int [num*NEFPerBlock(RBblock[j])];
//        e->iefb_sms[j] = new int [ num ];
            
        }
        numBoundaryFacesOfBlock = new int[tmpblk]; // a one dimensional array
        // storing the total number of faces in the perticular boundary block
    }
    
    /* Allocation Arrays for Essential boundary conditions */

    num = info->nshg;
    e->nBC = new int[num];
    e->iBC = new int[num];
    e->BC = new double* [num];
    for( int j =0 ; j< num; j++ ) 
        e->BC[j] = new double [ numEBC ];
    e->iper = new int[num];

    // zero the boundary condition arrays
    for (int i=0; i < info->nshg; i++){
        e->iBC[i] = 0;
        e->iper[i] = -1;
        e->nBC[i] = -1;
        for (int j=0; j < numEBC; j++)
            e->BC[i][j] = 0.0;
    }
    for (int i=0; i < info->nshg; i++)
        for (int j=0; j < numVars; j++)
            e->q[i][j] = 0.0;
}

void
DeAllocateEnsaArrays( EnsaArrays* e ) {
    delete [] e->x;
    delete [] e->nBC;
    delete [] e->iBC;
    delete [] e->iper;
    for( int i =0 ; i < Iblock.size(); i++ ) { 
        delete [] e->ien[i];
        delete [] e->ien_sms[i];
    }
    delete [] e->ien;
    delete [] e->ien_sms;
    for( int i =0 ; i < Bblock.size(); i++ ) {
        delete [] e->ienb[i];
        delete [] e->ienb_sms[i];
    }
    delete [] e->ienb;

    delete [] e->ienb_sms;

    if(FaceConnectivity){

        // deallocate memory for interrior and boundary face connectivity
        for( int i =0 ; i < Iblock.size(); i++ ) {
            delete [] e->ief[i];
            //delete [] e->ief_sms[i];
        }
        
        delete [] e->ief;
        //delete [] e->ief_sms;
        
        for( int i =0 ; i < Bblock.size(); i++ ) {
            delete [] e->iefb[i];
            //  delete [] e->iefb_sms[i];
        }
        
        delete [] e->iefb;

        delete [] numInteriorFacesOfBlock;
        delete [] numBoundaryFacesOfBlock;        

            //delete [] e->iefb_sms;
    }    // if(FaceConnectivity)




}

// what is going on here ?
// re-number indices (here: numbering is 0-based, fortran is 1-based)
// called right in the beginning of WriteEnsaArrays() 
void 
fortranIndexing( EnsaArrays* e ) {
    
    int tmpblk = Iblock.size();
    for(int j=0; j < tmpblk; j++){
        int num_in_block  = NIblock[j];
        int num_shp = Nshape[j];
        for (int k=0; k < num_in_block; k++)
            for(int i =0 ; i < num_shp; i++){
                if(e->ien[j][k+i*num_in_block] >=0) e->ien[j][k+i*num_in_block]++;
                else e->ien[j][k+i*num_in_block]--;
            }
    }

    tmpblk = Bblock.size();
    for(int j=0; j < tmpblk; j++){
        int num_in_block  = NBblock[j];
        int num_shp = NshapeB[j];
        for (int k=0; k < num_in_block ; k++)
            for(int i =0; i < num_shp; i++){
                if(e->ienb[j][k+i*num_in_block] >=0) e->ienb[j][k+i*num_in_block]++;
                else e->ienb[j][k+i*num_in_block]--;
            }
    }

    if(FaceConnectivity){
        // same for ief
        tmpblk = Iblock.size();
        for(int j=0; j < tmpblk; j++){
            int num_in_block  = NIblock[j];
            int num_faces = NEFPerBlock(RIblock[j]);
            for (int k=0; k < num_in_block; k++)
                for(int i =0 ; i < num_faces; i++){
                    if(e->ief[j][k+i*num_in_block] >=0) e->ief[j][k+i*num_in_block]++;
                    else e->ief[j][k+i*num_in_block]--;
                }
        }    
        
        // same for iefb
        tmpblk = Bblock.size();
        for(int j=0; j < tmpblk; j++){
            int num_in_block  = NBblock[j];
            int num_faces = NEFPerBlock(RBblock[j]);
            for (int k=0; k < num_in_block ; k++)
                for(int i =0; i < num_faces; i++){
                    if(e->iefb[j][k+i*num_in_block] >=0) e->iefb[j][k+i*num_in_block]++;
                    else e->iefb[j][k+i*num_in_block]--;
                }
        }
    }//  if(FaceConnectivity)
    
    for (int i=0; i < info->nshg; i++){
        e->iper[i]++;
        e->nBC[i]++;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called by 
// writeEnsaFiles in writeEnsaFiles.cc
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void 
writeEnsaArrays( EnsaArrays* e ,
                 EnsaParameters *par, 
                 double zscale[3] ) {
    
    int i,j,k,ind, *iBCBf;
    double *BCf, *BCBf, *qf;
    char fname[255];
    int tmpblk,tmpblkb, nelblk, bpoly, bnen, bnsh;

    fortranIndexing( e );
    
    int numnp  = M_numVertices( mesh );
    int numel  = M_numRegions( mesh );
    if (info->nsd == 2) numel = M_numFaces( mesh );
    int nshg   = info->nshg;
    int numelb = info->numelb;
    int numpbc = info->numpbc;
    int nshapeb= 4;
    int nsd    = par->getNSD();
    int nen    = par->getNEN();

    int numflx = par->getNUMFLX();
    int numVars  = par->getNUMVARS();
    int numEBC   = par->getNUMEBC();
    int numNBC   = par->getNUMNBC();
    int fType  = par->getFTYPE();


  
    BCf   = new double[numEBC*numpbc];
    ind = 0;
    for (j=0; j < numEBC; j++){
        for (i=0; i < numpbc; i++)
            BCf[ind++] = e->BC[i][j];
    }
    for(j=0; j < nshg; j++ ) delete [] e->BC[j];
    delete [] e->BC;

    qf    = new double[numVars*nshg];
    ind = 0;
    for (j=0; j < numVars; j++){
        for (i=0; i < nshg; i++)
            qf[ind++] = e->q[i][j];
    }
    for( j=0; j <nshg; j++ ) delete []  e->q[j]; 
    delete [] e->q;
    
    if (zScale) {
        for (i=0; i < numnp; i++) {
            e->x[i+0*numnp] *= zscale[0];
            e->x[i+1*numnp] *= zscale[1];
            e->x[i+2*numnp] *= zscale[2];
        }
    }
    int magic_number = 362436;
    int* mptr = &magic_number;
    
    int fgeom, frest;
    int iarray[10];
    int size, nitems;

    sprintf(fname,"restart.%d.1",lstep);
    openfile_( fname, "write", &frest );

    bzero( (void*)fname, 255 );
    sprintf(fname,"geombc.dat.1");
    openfile_( fname, "write", &fgeom );
    
    tmpblk = Iblock.size(); /* number of interior blocks  */
    tmpblkb = Bblock.size(); /* num boundary blocks  */

    /* before anything we put in the standard headers */

    writestring_( &fgeom,"# PHASTA Input File Version 2.0\n");
    writestring_( &frest,"# PHASTA Input File Version 2.0\n");

    writestring_( &fgeom, "# Byte Order Magic Number : 362436 \n");
    writestring_( &frest, "# Byte Order Magic Number : 362436 \n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# Output generated by NSpre version: %s \n", version);
    writestring_( &fgeom, fname );
    writestring_( &frest, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"# CmdLineArgs : %s \n", options);
    writestring_( &fgeom, fname );
    writestring_( &frest, fname );

    time_t timenow = time ( &timenow);

    bzero( (void*)fname, 255 );
    sprintf(fname,"# %s\n", ctime( &timenow ));
    writestring_( &fgeom, fname );
    writestring_( &frest, fname );

    int one=1;

    size = 1;
    nitems = 1;
    iarray[ 0 ] = 1;
    writeheader_( &fgeom, "byteorder magic number ",
                  (void*)iarray, &nitems, &size, "integer", oformat );

    writedatablock_( &fgeom, "byteorder magic number ",
                     (void*)mptr, &nitems, "integer", oformat );

    writeheader_( &frest, "byteorder magic number ",
                  (void*)iarray, &nitems, &size, "integer", oformat );

    writedatablock_( &frest, "byteorder magic number ",
                     (void*)mptr, &nitems, "integer", oformat );


    /* writing the restart */
    bzero( (void*)fname, 255 );
    sprintf(fname,"number of modes : < 0 > %d\n", nshg);
    writestring_( &frest, fname );
    
    bzero( (void*)fname, 255 );
    sprintf(fname,"number of variables : < 0 > %d\n", numVars);
    writestring_( &frest, fname );
    
    size =  numVars*nshg;
    nitems = 3;
    iarray[ 0 ] = nshg;
    iarray[ 1 ] = numVars;
    iarray[ 2 ] = lstep;

    writeheader_( &frest, "solution ",
                  ( void* )iarray, &nitems, &size,"double", oformat );

    nitems = numVars*nshg;
    writedatablock_( &frest, "solution ",
                     ( void* )(qf), &nitems, "double", oformat );


    closefile_( &frest, "write" );

    /* finished writing the restart */
    bzero((void*)fname, 255 );
    sprintf(fname,"number of processors : < 0 > 1 \n");
    writestring_( &fgeom, fname );
 
    bzero( (void*)fname, 255 );
    sprintf(fname,"number of variables : < 0 > %d \n", numVars);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of spatial dimensions : < 0 > %d \n", nsd);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of nodes : < 0 > %d \n", numnp);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf( fname,"number of nodes in the mesh : < 0 > %d \n",
             M_numVertices(mesh));
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf( fname,"number of edges in the mesh : < 0 > %d \n",
             M_numEdges(mesh));
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf( fname,"number of faces in the mesh : < 0 > %d \n",
             M_numFaces(mesh));
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of modes : < 0 > %d \n", nshg);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of shapefunctions solved on processor : < 0 > %d \n", nshg);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of global modes : < 0 > %d \n", nshgTot);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of interior elements : < 0 > %d \n", numel);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of boundary elements : < 0 > %d \n", numelb);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"maximum number of element nodes  : < 0 > %d \n", nen);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of interior tpblocks : < 0 > %d \n", tmpblk);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of boundary tpblocks : < 0 > %d \n", tmpblkb);
    writestring_( &fgeom, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of nodes with Dirichlet BCs : < 0 > %d \n", numpbc);
    writestring_( &fgeom, fname );

    size = nsd*numnp;
    nitems = 2 ;
    iarray[ 0 ] = numnp;
    iarray[ 1 ] = nsd;
    writeheader_( &fgeom, "co-ordinates ", (void*)iarray, &nitems, &size,
                  "double", oformat  );

    nitems =nsd*numnp;
    writedatablock_( &fgeom, "co-ordinates ", (void*)(e->x), &nitems,
                     "double", oformat );

    /* now we enter the loop over interior blocks where we write ien */


    char keyphrase[100]; /* to store the dynamic string for each tpblock */
    int bnenbl, bnshlb, blcsyst;
    for(int i=0; i< tmpblk; i++) {     /* for interior each block */
        bnen = RIblock[i].nen;      /* nen of the block -- topology */
        bpoly = RIblock[i].maxpoly; /* polynomial order of the block */
        nelblk = NIblock[i] ;       /* numel of this block */
        bnsh  = Nshape[i] ;         /* nshape of this block */ 
        bnshlb = getnumBshape(RIblock[i]);
        bnenbl = RIblock[i].nenbl;
        blcsyst = RIblock[i].lcsyst;

        /* generate the key phrase specific to this block */
        
        generate_keyphrase( keyphrase,"connectivity interior ", &RIblock[i]);

        size = nelblk*bnsh;
        nitems = 7;
        iarray[ 0 ] = nelblk;
        iarray[ 1 ] = bnen;
        iarray[ 2 ] = bpoly;
        iarray[ 3 ] = bnsh;
        iarray[ 4 ] = bnshlb;
        iarray[ 5 ] = bnenbl;
        iarray[ 6 ] = blcsyst;

        writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                      "integer", oformat );

        nitems = nelblk*bnsh;
        writedatablock_( &fgeom, keyphrase, (void*)(e->ien[i]), &nitems,
                         "integer", oformat );

        generate_keyphrase( keyphrase,"ien to sms ", &RIblock[i]);

        size = nelblk;
        nitems = 1;
        writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                      "integer", oformat );

        nitems = nelblk;
        writedatablock_( &fgeom, keyphrase, (void*)(e->ien_sms[i]), &nitems,
                         "integer", oformat );
    }


    /* similar to the above we also write the boundary blocks */
    /* along with these  we also write the Natural boundary conditions */

    for(i =0; i< tmpblkb; i++){

        bnen = RBblock[i].nen;
        bpoly = RBblock[i].maxpoly;
        bnenbl = RBblock[i].nenbl;
        blcsyst = RBblock[i].lcsyst;
        nelblk = NBblock[i]; /* num of elements in this block */
        bnsh =  NshapeB[i]; /* num of shapefuncs for each element */
        bnshlb = getnumBshape(RBblock[i]);

        iBCBf = new int [nelblk*2];
        BCBf = new double [nelblk*numNBC];

        for(k=0; k< nelblk; k++){
            for(j=0;j<2;j++) iBCBf[k+j*nelblk] = e->iBCB[i][k][j];     
            delete [] e->iBCB[i][k];
            for(j=0;j<numNBC;j++) BCBf[k+j*nelblk] =e->BCB[i][k][j];
            delete [] e->BCB[i][k];
        } 
        delete [] e->iBCB[i];
        delete [] e->BCB[i];

        iarray[ 0 ] = nelblk;
        iarray[ 1 ] = bnen;
        iarray[ 2 ] = bpoly;
        iarray[ 3 ] = bnsh;
        iarray[ 4 ] = bnshlb;
        iarray[ 5 ] = bnenbl;
        iarray[ 6 ] = blcsyst;
        iarray[ 7 ] = numNBC;

        generate_keyphrase( keyphrase,"connectivity boundary ",
                            &RBblock[i]);
        size = nelblk*bnsh;
	if (old_format) size = nelblk*bnsh + nelblk*2 + nelblk*numNBC;
        nitems = 8;
        writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                      "integer", oformat );

        nitems = nelblk*bnsh;
        writedatablock_( &fgeom, keyphrase, (void*)(e->ienb[i]), &nitems,
                         "integer", oformat );

       if (!old_format){
	   generate_keyphrase( keyphrase,"ienb to sms ", &RBblock[i]);

           size = nelblk;
           nitems = 1;
           writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                      "integer", oformat );

           nitems = nelblk;
           writedatablock_( &fgeom, keyphrase, (void*)(e->ienb_sms[i]), &nitems,
                         "integer", oformat );

           generate_keyphrase(keyphrase,"nbc codes ", &RBblock[i]);
           size = nelblk*2;
           nitems = 8;

           writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                      "integer", oformat );
       }

        nitems = nelblk*2;
        writedatablock_( &fgeom, keyphrase, (void*)(iBCBf), &nitems,
                         "integer", oformat );


        if (!old_format){
	   generate_keyphrase(keyphrase,"nbc values ", &RBblock[i]);
           size = nelblk*numNBC;
           nitems = 8;

           writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                      "double", oformat );
	}

        nitems = nelblk*numNBC;
        writedatablock_( &fgeom, keyphrase, (void*)(BCBf), &nitems,
                         "double", oformat );

        delete [] iBCBf;
        delete [] BCBf;
   
    } // end of the ienb
    


    if(FaceConnectivity){

        // ief
        // face connectivity of interior elements 
        //          
        // now we enter the loop over interior blocks where we write ief 
        for(int i=0; i< tmpblk; i++) {     /* for interior each block */
            bnen = RIblock[i].nen;      /* nen of the block -- topology */
            bpoly = RIblock[i].maxpoly; /* polynomial order of the block */
            nelblk = NIblock[i] ;       /* numel of this block */
            bnsh  = Nshape[i] ;         /* nshape of this block */ 
            bnshlb = getnumBshape(RIblock[i]);
            bnenbl = RIblock[i].nenbl;
            blcsyst = RIblock[i].lcsyst;
            

            // get the  number of each element's faces characteristic for Block
            // (related to topology)
            int nfaces =NEFPerBlock(RIblock[i]); 
            
   
            /* generate the key phrase specific to this block */
            generate_keyphrase( keyphrase,"face connectivity interior ", 
                                &RIblock[i]);

            size = nelblk*nfaces;
            // last figure is numInteriorFacesOfBlock
            nitems = 8;
            // same order and type of variables
            // as they are read in phasta (genblk.f) (check also parallel NSpre) 
            iarray[ 0 ] = nelblk;
            iarray[ 1 ] = nfaces;
            iarray[ 2 ] = bpoly;
            iarray[ 3 ] = bnsh;
            iarray[ 4 ] = bnshlb;
            iarray[ 5 ] = bnenbl;
            iarray[ 6 ] = blcsyst;
            iarray[ 7 ] = numInteriorFacesOfBlock[i]; // redundant
            
            writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                          "integer", oformat );
            
            nitems = nelblk*bnsh;
            writedatablock_( &fgeom, keyphrase, (void*)(e->ief[i]), &nitems,
                             "integer", oformat );



//          generate_keyphrase( keyphrase,"ief to sms ", &RIblock[i]);

//          size = nelblk;
//          nitems = 1;
//          writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
//                        "integer", oformat );

//          nitems = nelblk;
//          writedatablock_( &fgeom, keyphrase, (void*)(e->ief_sms[i]), &nitems,
//                           "integer", oformat );

        }// end of ief


        // face connectivity for boundary elements 
        // 
        for(int i =0; i< tmpblkb; i++){
            bnen = RBblock[i].nen;
            bpoly = RBblock[i].maxpoly;
            bnenbl = RBblock[i].nenbl;
            blcsyst = RBblock[i].lcsyst;
            nelblk = NBblock[i]; /* num of elements in this block */
            bnsh =  NshapeB[i]; /* num of shapefuncs for each element */
            bnshlb = getnumBshape(RBblock[i]);
            
            // get the total number of element faces characteristic for Block
            // (related to topology)
            int nfaces =bnen; 


            iarray[ 0 ] = nelblk;
            iarray[ 1 ] = bnen;
            iarray[ 2 ] = bpoly;
            iarray[ 3 ] = bnsh;
            iarray[ 4 ] = bnshlb;
            iarray[ 5 ] = bnenbl;
            iarray[ 6 ] = blcsyst;
            iarray[ 7 ] = numNBC;
            iarray[ 8 ] = numBoundaryFacesOfBlock[i];
            
            generate_keyphrase( keyphrase,"face connectivity boundary ", 
                                &RBblock[i]);
            size = nelblk*nfaces; // size of unformatted array to be written
            // last figure (9th) is  numBoundaryFacesOfBlock
            nitems = 9;
            writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
                          "integer", oformat );
            
            nitems = nelblk*bnsh;
            writedatablock_( &fgeom, keyphrase, (void*)(e->iefb[i]), &nitems,
                             "integer", oformat );
            
            
//          generate_keyphrase( keyphrase,"iefb to sms ", &RBblock[i]);

//          size = nelblk;
//          nitems = 1;
//          writeheader_( &fgeom, keyphrase, (void*)iarray, &nitems, &size,
//                        "integer", oformat );

//          nitems = nelblk;
//          writedatablock_( &fgeom, keyphrase, (void*)(e->iefb_sms[i]), &nitems,
//                           "integer", oformat );


        
        } // end of iefb

    } // if(FaceConnectivity)


    /* the unblocked Essential Boundary conditions */

    size = nshg;
    nitems = 1;
    iarray[ 0 ] = nshg;
    writeheader_( &fgeom , "bc mapping array ", (void *)iarray, &nitems,
                  &size, "integer", oformat );

    nitems = nshg;
    writedatablock_( &fgeom, "bc mapping array ", (void*)(e->nBC), &nitems ,
                     "integer", oformat );

    size = numpbc;
    nitems = 1;
    iarray[ 0 ] = numpbc;
    writeheader_( &fgeom , "bc codes array ", (void *)iarray, &nitems,
                  &size, "integer", oformat );

    nitems = numpbc;
    writedatablock_( &fgeom, "bc codes array ", (void*)(e->iBC), &nitems ,
                     "integer", oformat );

    size = numpbc*numEBC;
    nitems = 3;
    iarray[ 0 ] = numpbc*numEBC;
    iarray[ 1 ] = numpbc;
    iarray[ 2 ] = numEBC;
    writeheader_( &fgeom , "boundary condition array ", (void *)iarray, &nitems,
                  &size, "double", oformat );

    nitems = numpbc*numEBC;
    writedatablock_( &fgeom, "boundary condition array ", (void*)(BCf),
                     &nitems , "double", oformat );

    size = nshg;
    nitems = 1;
    iarray[ 0 ] = nshg;
    writeheader_( &fgeom , "periodic masters array ", (void *)iarray, &nitems,
                  &size, "integer", oformat );

    nitems = nshg;
    writedatablock_( &fgeom, "periodic masters array ", (void*)(e->iper),
                     &nitems , "integer", oformat );
    closefile_( &fgeom, "write" );

    delete[] Nshape;
    delete[] NshapeB;
    delete[] NIblock;
    delete[] NBblock;
    delete[] RIblock;
    delete[] RBblock;
    delete[] BCf;
    delete[] qf;
}
