/********************************************************************/
/* Writes results to ascii format files :                           */
/* solution.asc.out, acceleration.asc.out, bflux.asc.out,           */
/* stats.asc.out, errors.asc.out                                    */
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

using namespace std;

void asciiOutput(int* array, double* qglobal, double* xglobal, int* ien,
	 	 int RequestedField){
    FILE* fgeombc;
    FILE* frest;
    int i,j;
    int nshgtot, numvar, lstep, nsd, neltot;
    int nen, nelblk, nendx, ipordl, nshl, step1, nTimeStep;
    
    numvar    = array[0];
    nshgtot   = array[1];
    lstep     = array[2];
    neltot    = array[4];
    nendx     = array[5];
    nelblk    = array[6];
    nsd       = array[7];
    ipordl    = array[8];
    nshl      = array[9];
    nen       = array[11];
    step1     = array[14];
    nTimeStep = array[15];
    /* echo out restart results */
      
    /* velocity field */
    if(RequestedField == 0){
        frest = fopen("solution.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    /* acceleration field */
    } else if(RequestedField == 1){	
        frest = fopen("acceleration.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    /* boundary fluxes field */
    } else if(RequestedField == 2){	
        frest = fopen("bflux.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    /* turbulence statistics field */
    } else if(RequestedField == 3){	
        frest = fopen("stats.asc.out","w");
        fprintf(frest,"%d %d %d %d %d\n", nshgtot, numvar, lstep, step1, nTimeStep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    /* error field */
    } else if(RequestedField == 4){	
        frest = fopen("errors.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    
    /* velbar field */
    } else if(RequestedField == 5){
        frest = fopen("velbar.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    

    /* ybar field */
    } else if(RequestedField == 6){	
        frest = fopen("ybar.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    
    /* wall shear stress field */
    } else if(RequestedField == 7){	
        frest = fopen("wss.asc.out","w");
        fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);
      
        for(i=0; i< nshgtot; i++){
            for(j=0; j < numvar; j++) 
                fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
            fprintf(frest,"\n");
        } 
      
        fclose(frest);
    
    }    
    if(RequestedField != 1){
        /* echo out geom results */
        fgeombc = fopen("geombc.asc.out","w");
        fprintf(fgeombc, "%d %d %d %d %d %d \n",
            nshgtot, nshgtot, nsd, neltot, nen, nelblk);
        
        for(i=0; i< nshgtot; i++){
            for(j=0; j < 3; j++) 
                fprintf(fgeombc,"%22.7e  ",xglobal[j*nshgtot+i]);
            fprintf(fgeombc,"\n");
        } 
    
        fprintf(fgeombc, "%d %d %d %d \n",
              neltot, nendx, ipordl, nshl);
        for(i=0; i< neltot; i++){
            for(j=0; j < nendx; j++) 
                fprintf(fgeombc,"%d  ",ien[j*neltot+i]);
            fprintf(fgeombc,"\n");
        } 
        fclose(fgeombc);
    }
    return;
}
