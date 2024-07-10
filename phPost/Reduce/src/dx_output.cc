/********************************************************************/
/* Formating necessary arrays for dx format                         */
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
#include "reduce.h"
#include "phastaIO.h"

using namespace std;

void dxOutput(int* array, double* qglobal, double* xglobal, int* ien,
	      int RequestedField){
    
    float *qdx, *xdx;
    int *iendx, nodes[8];
    int numvar, numnp, neltot, nendx, nenlmin,nshgtot, lstep, nsd;
    int i,j;
    
    numvar  = array[0];
    nshgtot = array[1];
    lstep   = array[2];
    numnp   = array[3];
    neltot  = array[4];
    nendx   = array[5];
    nsd     = array[7];
    nenlmin = array[10];

    xdx = (float *) malloc( nsd*numnp * sizeof(float));
    qdx = (float *) malloc( numvar*numnp * sizeof(float));
    iendx = (int *) malloc( nendx*neltot * sizeof(int));
    
    /* transpose and float data for  dx */    
    for(i=0; i< numnp; i++)
        for(j=0; j < numvar; j++) 
            qdx[i*numvar+j]=(float)qglobal[j*nshgtot+i];
    
    for(i=0; i< numnp; i++)
        for(j=0; j < nsd; j++) 
            xdx[i*nsd+j]=(float)xglobal[j*numnp+i];
    
    /* connectivity in dx is rather strange. */
      
    if (nendx==nenlmin) {   
        for(i=0; i< neltot; i++)
            for(j=0; j < nendx; j++) 
                iendx[i*nendx+j]=ien[j*neltot+i]-1;
    }
    /* multitopology */
    else {
        for(i=0; i< neltot; i++){
          
            int junique=4;
            for(j=4; j < nendx; j++)
                if(ien[j*neltot+i] != ien[(j-1)*neltot+i]) junique=j+1;

            /* is this a tet */
        
            if(junique==4 && nendx==8) {
                iendx[i*nendx+7]=ien[3*neltot+i]-1; 
                iendx[i*nendx+6]=ien[3*neltot+i]-1; 
                iendx[i*nendx+5]=ien[3*neltot+i]-1; 
                iendx[i*nendx+4]=ien[3*neltot+i]-1; 
                iendx[i*nendx+3]=ien[1*neltot+i]-1; 
                iendx[i*nendx+2]=ien[2*neltot+i]-1; 
                iendx[i*nendx+1]=ien[2*neltot+i]-1; 
                iendx[i*nendx+0]=ien[0*neltot+i]-1; 
            } 
       
            /* is this a wedge */
          
            else if(junique==6) {
                iendx[i*nendx+7]=ien[5*neltot+i]-1; 
                iendx[i*nendx+6]=ien[5*neltot+i]-1; 
                iendx[i*nendx+5]=ien[4*neltot+i]-1; 
                iendx[i*nendx+4]=ien[3*neltot+i]-1; 
                iendx[i*nendx+3]=ien[2*neltot+i]-1; 
                iendx[i*nendx+2]=ien[2*neltot+i]-1; 
                iendx[i*nendx+1]=ien[1*neltot+i]-1; 
                iendx[i*nendx+0]=ien[0*neltot+i]-1; 
            } 
       
            /*        if(junique==8) {   this is a hex */ 
   
            else{
                for(j=0; j < nendx; j++)
                    iendx[i*nendx+j]=ien[j*neltot+i]-1; 
            } 
        }
    }


    if(nendx==8) {
        for(i=0; i< neltot; i++){
            nodes[0]=iendx[i*nendx+4];
            nodes[1]=iendx[i*nendx+0];
            nodes[2]=iendx[i*nendx+7];
            nodes[3]=iendx[i*nendx+3];
            nodes[4]=iendx[i*nendx+5];
            nodes[5]=iendx[i*nendx+1];
            nodes[6]=iendx[i*nendx+6];
            nodes[7]=iendx[i*nendx+2];
            for(j=0; j < nendx; j++) 
                iendx[i*nendx+j]=nodes[j];
        }
    }
    // DX on littleEndian machines
    // does an Endian byteswap so we have to do the same to come out right

    if ( isLittleEndian() ) {
		printf("\n\nSwapping to BigEndian for DX output\n\n");
        SwapArrayByteOrder( iendx, sizeof(int), neltot*nendx );
        SwapArrayByteOrder( xdx, sizeof(float), nsd*numnp );
        SwapArrayByteOrder( qdx, sizeof(float), numvar*numnp );
    }
    
    wrtc_(numvar, neltot, iendx, lstep, nsd,
          numnp, xdx, qdx, nendx, RequestedField);
    
    free(xdx);
    free(qdx);
    free(iendx);
    return;
}


