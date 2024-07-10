/********************************************************************/
/* Reading geometry from geombc.dat.nproc and reducing the nodal    */
/* coordinate array to one processor using connectivity array       */
/*                                                                  */
/* Elaine Bohr                                                      */
/* February 2004                                                    */
/********************************************************************/
#include <iostream>
#include <stdio.h>
#include <string>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "phastaIO.h"
#include "reduce.h"

using namespace std;

extern int* sms2edge;
extern int* sms2face;
extern int esize,fsize,emodes,fmodes;

void  bzero(void* ptr, size_t sz) throw() {
    int i;
    char *cptr;
    cptr = (char*) ptr;
    for (i=0; i < sz; i++) {
        cptr[i]=0;
    }
    return;
}

	    
int geometry_and_connectivity(int* array, double* &xglobal, int* &ien, 
                               int** &ncorp2d, bool RequestedVolCheck){
			       
    FILE *volcheck;
    int *ient;
    int nshgl, numnp;
    int maxnp,maxnshg, numprocs , nshgtot,i,j,k,l;
    int np, numel, nen, nelblk, nsdread;
    int neltot, nelsofar, neltp, nenl, ipordl, nshl;
    int gnum, opnum, nendx, nenlmin;
    int iarray[MAAXARRAY];
    double *xlocal;
    char gfname[40];
    char* iotype;
    int ione=1, itwo=2, ithree=3,  iseven=7;
    int igeom;
    int ixsiz, nsd, iientsiz;
    double min[MAXNUMVAR], max[MAXNUMVAR];  
    double e1[3],e2[3],e3[3];
    double vol;
    int vtxnum0,vtxnum1,vtxnum2,vtxnum3;

    igeom=1;
    iotype="binary";

    for(i=0; i< MAXNUMVAR; i++){
        max[i]=-100000000;
        min[i]= 100000000;
    }

    /* Scanning geombc.dat.1 for number of processors, total number of nodes */
    /* and number of global modes */
    sprintf(gfname,"geombc.dat.%d",1); /* geometry and bc database */
    cout << "Opening " << gfname << " to scan for number of processors and global modes..."<<endl;
    openfile(gfname, "read", &igeom);
    iarray[0] = -1;
    readheader(&igeom,"number of processors",(void*)iarray,&ione,"integer",iotype);
    numprocs=iarray[0];
    if(iarray[0] == -1){
        cout << "number of processors not found...exiting" << endl;
	return 1;
    }
    iarray[0] = -1;
    readheader(&igeom,"number of global modes",(void*)iarray,&ione,"integer",iotype); 
    nshgtot = iarray[0];
    if(iarray[0] == -1){
        cout << "number of global modes not found...exiting" << endl;
	return 1;
    }
    iarray[0] = -1;
    readheader(&igeom,"number of nodes in the mesh",(void*)iarray,&ione,"integer",iotype); 
    numnp = iarray[0];
    if(iarray[0] == -1){
        cout << "number of nodes in the mesh not found..." << endl;
	cout << "please provide it " << endl;
	cin >> numnp;
    } 
    iarray[0] = -1;
    readheader(&igeom,"number of spatial dimensions",(void*)iarray,&ione,"integer",iotype); 
    nsdread = iarray[0];
    if(iarray[0] == -1){
        cout << "number of spatial dimensions not found...only works for 3D" << endl;
    }

    iarray[0] = -1 ;
    readheader( &igeom, "edge mode mapping to sms",(void*)iarray, &itwo, "integer", iotype);
    if ( iarray[0] > 0 ) {
        esize = iarray[0]; emodes = iarray[1];
        sms2edge = new int [ esize ];
        readdatablock( &igeom, "edge mode mapping to sms", (void*)sms2edge, &esize, "integer", iotype );
    }
                                                 
    iarray[0] = -1 ;
    readheader( &igeom, "face mode mapping to sms",(void*)iarray, &itwo, "integer", iotype);
    if ( iarray[0] > 0 ) {
        fsize = iarray[0]; fmodes = iarray[1];
        sms2face = new int [ fsize ];
        readdatablock( &igeom, "face mode mapping to sms", (void*)sms2face, &fsize, "integer", iotype );
    }

    closefile(&igeom, "read");

  /* scanning geom.dat.<procnum> to add up neltot and to determine maxnshg */
    neltot=0;
    maxnshg=0;
    maxnp=0;
    for(i=0; i< numprocs; i++){
        bzero( (void*)gfname, 40);
        sprintf(gfname,"geombc.dat.%d",i+1); /* geometry and bc database */
        openfile(gfname, "read", &igeom);
        iarray[0] = -1;
	readheader(&igeom,"number of interior elements",(void*)iarray,&ione,"integer",iotype); 
        numel=iarray[0];
        if(iarray[0] == -1){
            cout << "number of interior elements not found...exiting" << endl;
	    return 1;
        }
        iarray[0] = -1;
	readheader(&igeom,"maximum number of element nodes",(void*)iarray,&ione,"integer",iotype);
        nen=iarray[0];
        if(iarray[0] == -1){
            cout << "maximum number of element nodes not found...exiting" << endl;
	    return 1;
        }
        iarray[0] = -1;
	readheader(&igeom,"number of modes",(void*)iarray,&ione,"integer",iotype); 
        nshgl=iarray[0];
        if(iarray[0] == -1){
            cout << "number of modes not found...exiting" << endl;
	    return 1;
        }
        if(nshgl>maxnshg) maxnshg=nshgl;
        iarray[0] = -1;
	readheader(&igeom,"number of nodes",(void*)iarray,&ione,"integer",iotype); 
        np=iarray[0];
        if(iarray[0] == -1){
            cout << "number of nodes not found...exiting" << endl;
	    return 1;
        }
        if(np>maxnp) maxnp=np;

        closefile(&igeom, "read" );
        neltot+=numel;
    }
    
    /* now that we have our size of array information we can allocate
     our  local and global geometry arrays (here local means on a
     given processor and global means the total array assembled 
     across all processors) */

    xglobal = (double *) malloc( 3*numnp * sizeof(double));
    ncorp2d = (int ** )malloc(sizeof(int *)*numprocs);
    nendx=nen;
    if(nen>4) nendx=8;  /* DX thinks of things as ALL tets or hexes ONLY */
    ien = (int *) malloc( nendx*neltot * sizeof(int));
    nenlmin=nendx;
    nelsofar=0;

    if(RequestedVolCheck){
        volcheck = fopen("volcheck.dat","w");
        fprintf(volcheck,"volume  row of x's, then y's then z's \n");
    }

    /*Next we loop over the processors and read each processors
    geometry database.  Using the ncorp2d array, we reconstruct the
    global geometry structures (coordinates and connectivity) */

 
    for(i=0; i< numprocs; i++){

        /* open geom file  and read header*/      
        bzero( (void*)gfname, 40);
        sprintf(gfname,"geombc.dat.%d",i+1);
        printf("Reducing : %s \n", gfname);
        openfile( gfname, "read", &igeom );

        readheader(&igeom,"number of nodes",(void*)iarray,&ione,"integer", iotype);
        np=iarray[0];
        readheader(&igeom,"number of modes",(void*)iarray,&ione,"integer", iotype);
        nshgl=iarray[0];
        iarray[0] = -1;
	readheader(&igeom,"number of interior tpblocks",(void*)iarray,&ione,"integer",iotype);
        if(iarray[0] == -1){
            cout << "number of interior tpblocks not found...exiting" << endl;
	    return 1;
        }
	nelblk=iarray[0];

        /* read coordinates and fill into global array */
        iarray[0] = -1;
	readheader(&igeom,"co-ordinates",(void*)iarray,&itwo,"double",iotype);
        np=iarray[0];
        nsd=iarray[1];
	if(iarray[0] == -1){
            cout << "co-ordinates not found...exiting" << endl;
	    return 1;
	}
        ixsiz=np*nsd;
        xlocal = (double*)malloc( ixsiz * sizeof(double) );
        readdatablock(&igeom,"co-ordinates",(void*)xlocal,&ixsiz, "double",iotype);
        /* get the map from partion numbering to global numbering */
        if ( numprocs > 1 ) {
            iarray[0] = -1;
	    readheader(&igeom,"mode number map from partition to global",(void*)iarray,&ione,"integer",iotype);
            nshgl=iarray[0];
            if(iarray[0] == -1){
                cout << "mode number map from partition to global not found...exiting" << endl;
	        return 1;
            }
	    ncorp2d[i] = (int * )malloc(sizeof(int)*nshgl);
            readdatablock(&igeom,"mode number map from partition to global",
                           (void*)ncorp2d[i],&nshgl,"integer",iotype);
        } else {
            ncorp2d[i] = (int * )malloc(sizeof(int)*nshgl);
            for(j=0; j< nshgl ; j++) 
                ncorp2d[i][j]=j+1;
        }
      
        /* map it to global numbering */
	    reduce(nsd, np, numnp, i, xlocal, xglobal, ncorp2d, max, min);  

        /*read connectivity data */
        for(k=0; k< nelblk; k++){
            /* keyphrase identifying interior connectivity element block */
            iarray[0] = -1;
	    readheader(&igeom,"connectivity interior",(void*)iarray,&iseven,"integer",iotype);
            neltp  =iarray[0];
            nenl   =iarray[1];
            ipordl =iarray[2];
            nshl   =iarray[3];
	    if(iarray[0] == -1){
                cout << "connectivity interior not found...exiting" << endl;
	        return 1;
            }

            if(nenl < nenlmin) nenlmin=nenl;
            /* allocate the array to the right size */
            ient = (int *) malloc( nshl*neltp * sizeof(int));
            iientsiz=neltp*nshl;
            /* now read the array */
            readdatablock(&igeom,"connectivity interior",(void*)ient,&iientsiz,"integer", iotype);


            /* Now we need to bring ien to the global numbering */
            for(l=0; l< neltp; l++){
                for(j=0; j<nenl ; j++){
                    gnum=nelsofar+j*neltot+l;
                    opnum=ient[j*neltp+l];
                    ien[gnum]=ncorp2d[i][opnum-1];
                }

                if(RequestedVolCheck){
                    if(nshl==4){
                        vtxnum3=ient[3*neltp+l]-1;
                        vtxnum2=ient[2*neltp+l]-1;
                        vtxnum1=ient[neltp+l]-1;
                        vtxnum0=ient[l]-1;
                        for(j=0; j<3; j++)
                            e3[j]= xlocal[j*nshgl+vtxnum3]-xlocal[j*nshgl+vtxnum0]; 
                        for(j=0; j<3; j++)
                            e2[j]= xlocal[j*nshgl+vtxnum2]-xlocal[j*nshgl+vtxnum0]; 
                        for(j=0; j<3; j++)
                            e1[j]= xlocal[j*nshgl+vtxnum1]-xlocal[j*nshgl+vtxnum0]; 

                        vol=-(e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])
                              -e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])
                              +e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))/6;

                        if(vol<1.0e-10){
                            fprintf(volcheck,"%10.8e \n",vol);
                            for (j=0; j<3; j++)
                                fprintf(volcheck,"%10.8e %10.8e %10.8e %10.8e \n",
                                        xlocal[j*nshgl+vtxnum0],
                                        xlocal[j*nshgl+vtxnum1],
                                        xlocal[j*nshgl+vtxnum2],
                                        xlocal[j*nshgl+vtxnum3]);
                        }
                    }
                }
                /* pad to largest nendx since dx has to treat pyr and wedg as
                   degenerate hex */
                opnum=ien[gnum];  /* hijack opnum to keep last real node */
                for(j=nenl; j<nendx ; j++){
                    gnum=nelsofar+j*neltot+l;
                    ien[gnum]=opnum;
                }
            }
            nelsofar+=neltp;
            free(ient);
        }
        free(xlocal);
        closefile( &igeom, "read" );
    }
    array[0]  = 0;
    array[1]  = nshgtot;
    array[2]  = 0;
    array[3]  = numnp;
    array[4]  = neltot;
    array[5]  = nendx;
    array[6]  = nelblk;
    array[7]  = nsd;
    array[8]  = ipordl;
    array[9]  = nshl;
    array[10] = nenlmin;
    array[11] = nen;
    array[12] = numprocs;
    array[13] = maxnshg;
    return 0;
}
