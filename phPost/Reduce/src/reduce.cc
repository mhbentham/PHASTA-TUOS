/********************************************************************/
/* Reduce takes solution, acceleration, boundary fluxes, turbulence */
/* statistics or errors from restart.stepnumber.nproc and reduce    */
/* them to one processor. Results can be written in :               */
/* phasta format  - restart.stepnumber.0                            */
/* ascii format   - solution.asc.out, acceleration.asc.out,         */
/*                  bflux.asc.out, stats.asc.out, errors.asc.out    */
/* dx format      - solution.bin, acceleration.bin, bflux.bin,      */
/*                  stats.bin, errors.bin                           */
/* tecplot format - plottec.dat (only solution implemented)         */
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
extern char fieldname[];

/* This routine assembles the global array from the local (processor */
/* specific) array                                                   */
void reduce(int numvar, int lsize, int gsize, int proc,
	    double* local, double* global, int** ncorp,
	    double* qmax, double* qmin){
    int i,j;
    for(i=0; i< numvar; i++){
        for(j=0; j< lsize ; j++){ 
            global[i*gsize+ncorp[proc][j]-1] = local[i*lsize+j];
            if(qmax[i] < local[i*lsize+j]) qmax[i]=local[i*lsize+j];
            if(qmin[i] > local[i*lsize+j]) qmin[i]=local[i*lsize+j];
        }
    }
    return;
}


int main(int argc, char* argv[]){
    
    int stepnumber;
    int **ncorp2d, *ien;
    int nshgl, numvar,numnp;
    int maxnshg, numprocs , nshgtot,i,k, lstep, nTimeStep;
    int nen, nelblk;
    int neltot, ipordl, nshl;
    int nendx, nenlmin;
    int iarray[MAAXARRAY];
    int requestedField = 0;
    int requestedOutput = 0;
    int clread=0, step1=0;

    double *qglobal, *qlocal, *xglobal, *aglobal;
    char rfname[40];
    char* iotype;
    char* filename="restart";
    int ione=1, itwo=2, ithree=3, ifive=5, iseven=7;
    int igeom;
    int irstin;
    int nsd, iqsiz;
    
    double qmax[MAXNUMVAR],qmin[MAXNUMVAR],min[MAXNUMVAR], max[MAXNUMVAR];  

    bool RequestedAcceleration = false;
    bool RequestedSolution = true;  /* by default solution is reduced */
    bool VolCheck = false;
    bool RequestedDX = false;
    bool RequestedASCII = false;
    bool RequestedPhasta = true;    /* by default restart.#.0 is the output */
    bool RequestedBoundaryFluxes = false;
    bool RequestedTurbStats = false;
    bool RequestedErrors = false;
    bool RequestedTecplot = false;
    bool RequestedVelbar = false;
    bool RequestedYbar = false;
    bool RequestedWSS = false;

    for(i=0; i< MAXNUMVAR; i++){
        qmax[i]=-100000000;
        qmin[i]= 100000000;
         max[i]=-100000000;
         min[i]= 100000000;
    }

    igeom=1;
    irstin=2;
    iotype="binary";

    
    /* Processing command line arguments */
    clread = processCommandLineArguments(argc, argv, &requestedField, 
                                         &requestedOutput, &VolCheck,
					 &stepnumber);
    
    if (clread == 1) return(0);
    
    /* Acceleration and solution will be reduced */
    if (requestedField == 1) 
        RequestedAcceleration = true;

    /* Boundary fluxes will be reduced */
    if (requestedField == 2){
	RequestedBoundaryFluxes = true;
	RequestedSolution = false;
    }	

    /* Trubulence statistics will be reduced */
    if (requestedField == 3){
	RequestedTurbStats = true;
	RequestedSolution = false;
    }	

    /* Errors will be reduced */
    if (requestedField == 4){
	RequestedErrors = true;
	RequestedSolution = false;
    }	

    /* Velbar will be reduced */
    if (requestedField == 5){
	RequestedVelbar = true;
	RequestedSolution = false;
    }	

    /* Ybar will be reduced */
    if (requestedField == 6){
	RequestedYbar = true;
	RequestedSolution = false;
    }	
    /* Wall shear stress will be reduced */
    if (requestedField == 7){
	RequestedWSS = true;
	RequestedSolution = false;
    }	
    

    /* DX files will be written */
    if (requestedOutput == 1){ 
        RequestedDX = true;
	RequestedPhasta = false;
    }

    /* ASCII files will be written */
    if (requestedOutput == 2){ 
        RequestedASCII = true;
	RequestedPhasta = false;
    }

    /* TecPlot files will be written */
    if (requestedOutput == 3){ 
        RequestedTecplot = true;
	RequestedPhasta = false;
    }

    /* Calculates global nodal coordinates, ien and ncorp2d arrays*/
    clread = geometry_and_connectivity(iarray, xglobal, ien, ncorp2d, VolCheck);
    
    if (clread == 1) return(0);

    nshgtot  = iarray[1];
    numnp    = iarray[3];
    neltot   = iarray[4];
    nendx    = iarray[5];
    nelblk   = iarray[6];
    nsd      = iarray[7];
    ipordl   = iarray[8];
    nshl     = iarray[9];
    nenlmin  = iarray[10];
    nen      = iarray[11];
    numprocs = iarray[12];
    maxnshg  = iarray[13];

    /* Allocation of output global and local arrays*/
    if (RequestedSolution) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
	readheader(&irstin,"solution",(void*)iarray,&ithree,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No solution in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
	numvar=iarray[1];  /* pushing the second int into numvar */
        closefile( &irstin, "read" ); //the check below is not implemented-OK?

        qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
	aglobal = (double *) malloc(1*sizeof(double));
	aglobal[0] =0;
    }

    if(RequestedAcceleration){
        free(aglobal);
	aglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
    }

    if (RequestedBoundaryFluxes) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
        readheader(&irstin,"boundary flux",(void*)iarray,&ithree,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No boundary flux in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
        numvar=iarray[1];  /* pushing the second int into numvar */
        closefile( &irstin, "read" ); 

        qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
    }

    if (RequestedTurbStats) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
        readheader(&irstin,"statistics",(void*)iarray,&ifive,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No statistics in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
        numvar=iarray[1];  /* pushing the second int into numvar */
	step1=iarray[2];
	lstep=iarray[3];
	nTimeStep = iarray[4];
        closefile( &irstin, "read" ); 

        qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
    }

    if (RequestedErrors) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
	readheader(&irstin,fieldname,(void*)iarray,&ithree,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No errors in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
	numvar=iarray[1];  /* pushing the second int into numvar */
        closefile( &irstin, "read" ); //the check below is not implemented-OK?

        qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
    }

    if (RequestedVelbar) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
	readheader(&irstin,"velbar",(void*)iarray,&ithree,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No velbar in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
        numvar=iarray[1];
        closefile( &irstin, "read" ); //the check below is not implemented-OK?

        qglobal = (double *) malloc(numvar*nshgtot *sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
    }

    if (RequestedYbar) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
	readheader(&irstin,"ybar",(void*)iarray,&ithree,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No ybar in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
	numvar=iarray[1];  /* pushing the second int into numvar */
        closefile( &irstin, "read" ); //the check below is not implemented-OK?

        qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
    }    

    if (RequestedWSS) {
        /* scanning restart.<stepnum>.1 for numvar */
        sprintf(rfname,"%s.%d.1",filename,stepnumber);
        cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
        openfile(rfname, "read", &irstin   );
        iarray[0] = -1;
	readheader(&irstin,"wall shear stresses",(void*)iarray,&ithree,"double",iotype);
	if(iarray[0]==-1){
	   cout << "No wall shear stresses in restart." << stepnumber <<".1 -- Exiting"<<endl;
	   return 1;
	}   
	numvar=iarray[1];  /* pushing the second int into numvar */
        closefile( &irstin, "read" ); //the check below is not implemented-OK?

        qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
        qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));
    }


    /* We loop over the processors and read each processors        */
    /* solution database.  Using the ncorp2d array, we reconstruct */
    /* the global solution data                                    */
    for(i=0; i<numprocs; i++){
        /* read in solution for current processor */
	if(RequestedSolution){
            sprintf(rfname,"%s.%d.%d",filename,stepnumber, i+1);
            printf("Reducing : %s \n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"solution",(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"solution",(void*)qlocal, &iqsiz, "double", iotype);
            closefile( &irstin, "read" ); 
    
            /* map solution to global */      
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }
        /* read in time derivative of solution for current processor */
        if(RequestedAcceleration){
            for(k=0;k<numvar*maxnshg;k++){ // assigning first zero acceleration
                qlocal[k]=0.0;
            }
            sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
            printf("Reducing : %s \n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"time derivative of solution",(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"time derivative of solution",(void*)qlocal,&iqsiz, "double", iotype);
            closefile(&irstin, "read");
      
            /* map solution to global */    
	    reduce(numvar, nshgl, nshgtot, i, qlocal, aglobal, ncorp2d, max, min);  
        }

        /* read in boundary fluxes for current processor */
        if(RequestedBoundaryFluxes){
            sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
            printf("Reducing : %s for boundary fluxes\n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"boundary flux",(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"boundary flux",(void*)qlocal,&iqsiz, "double", iotype);
            closefile(&irstin, "read");

            /* map solution to global */    
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }

        /* read in turbulence statistics for current processor */
        if(RequestedTurbStats){
            sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
            printf("Reducing : %s for statistics\n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"statistics",(void*)iarray,&ifive,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[3];
	    nTimeStep = iarray[4];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"statistics",(void*)qlocal,&iqsiz, "double", iotype);
            closefile(&irstin, "read");

            /* map solution to global */    
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }

        /* read in error field for current processor */
	if(RequestedErrors){
            sprintf(rfname,"%s.%d.%d",filename,stepnumber, i+1);
            printf("Reducing : %s \n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,fieldname,(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,fieldname,(void*)qlocal, &iqsiz, "double", iotype);
            closefile( &irstin, "read" ); 
    
            /* map solution to global */      
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }
        /* read in velbar for current processor */
	if(RequestedVelbar){
            sprintf(rfname,"%s.%d.%d",filename,stepnumber, i+1);
            printf("Reducing : %s \n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"velbar",(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"velbar",(void*)qlocal, &iqsiz, "double", iotype);
            closefile( &irstin, "read" ); 
    
            /* map solution to global */      
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }
        /* read in ybar field for current processor */
	if(RequestedYbar){
            sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
            printf("Reducing : %s \n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"ybar",(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"ybar",(void*)qlocal, &iqsiz, "double", iotype);
            closefile( &irstin, "read" ); 
    
            /* map solution to global */      
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }
        /* read in wall shear stresses field for current processor */
	if(RequestedWSS){
            sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
            printf("Reducing : %s \n", rfname);
            openfile(rfname, "read", &irstin   );
            readheader(&irstin,"wall shear stresses",(void*)iarray,&ithree,"double",iotype);
            nshgl=iarray[0];
            numvar=iarray[1];
            lstep=iarray[2];
            iqsiz=nshgl*numvar;
            readdatablock(&irstin,"wall shear stresses",(void*)qlocal, &iqsiz, "double", iotype);
            closefile( &irstin, "read" ); 
    
            /* map solution to global */      
	    reduce(numvar, nshgl, nshgtot, i, qlocal, qglobal, ncorp2d, qmax, qmin);  
        }
    }
    
    /* Printing minimum and maximum for each variable */
    for(k=0; k< numvar; k++)
        printf("var.=%d, min=%f, max=%f \n", k,qmin[k],qmax[k]);

    iarray[0] = numvar;
    iarray[1]  = nshgtot;
    iarray[2] = lstep;
    iarray[3]  = numnp;
    iarray[4]  = neltot;
    iarray[5]  = nendx;
    iarray[6]  = nelblk;
    iarray[7]  = nsd;
    iarray[8]  = ipordl;
    iarray[9]  = nshl;
    iarray[10] = nenlmin;
    iarray[11] = nen;
    iarray[12] = numprocs;
    iarray[13] = maxnshg;
    iarray[14] = step1;
    iarray[15] = nTimeStep;

    /* Outputing the results */
    if(RequestedPhasta){
	phastaOutput(iarray, qglobal, aglobal, requestedField);
    }

    if(requestedField == 1){
        if(RequestedDX){
	    dxOutput(iarray, qglobal, xglobal, ien, 0);
	    dxOutput(iarray, aglobal, xglobal, ien, 1);
        }
    
        if(RequestedTecplot){
	    tecplotOutput(iarray, qglobal, xglobal, ien, 0);
	    tecplotOutput(iarray, aglobal, xglobal, ien, 1);
        }

        if(RequestedASCII){
            asciiOutput(iarray, qglobal, xglobal, ien, 0);
            asciiOutput(iarray, aglobal, xglobal, ien, 1);
        }
    } else{
        if(RequestedDX){
	    dxOutput(iarray, qglobal, xglobal, ien, requestedField);
        }
    
        if(RequestedTecplot){
	    tecplotOutput(iarray, qglobal, xglobal, ien, requestedField);
        }

        if(RequestedASCII){
            asciiOutput(iarray, qglobal, xglobal, ien, requestedField);
        }
    }


    /* Desallocating arrays */
    free(qlocal);
    for(i=0; i< numprocs; i++) free(ncorp2d[i]);
    free(ncorp2d);
    free(qglobal);
    if( RequestedAcceleration){
        free(aglobal);
    }
    free(xglobal);
    free(ien);
    return 0;
}
