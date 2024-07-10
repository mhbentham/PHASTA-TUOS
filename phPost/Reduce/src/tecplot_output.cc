/********************************************************************/
/* ascii tecplot format - plottec.dat (only solution implemented)   */
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

using namespace std;

void tecplotOutput(int* array, double* qglobal, double* xglobal, 
                   int* ien, int RequestedField){
    
    FILE* frest;
    int numvar, neltot, nendx, nshgtot, nen;
    int i,j, k;
    
    numvar  = array[0];
    nshgtot = array[1];
    neltot  = array[4];
    nendx   = array[5];
    nen     = array[11];

    if (RequestedField != 0) {
        cout << endl;
	cout << "Not yet implemented for this field." << endl;
	return;
    }
    frest = fopen("plottec.dat","w");

    if(nen == 8) {
        fprintf(frest,"Title = \"FE-VOLUME BRICK DATA SET\"\n");
        fprintf(frest,"VARIABLES = \"X\", \"Y\", \"Z\", \"P\", \"U\", \"V\", \"W\", \"T\"\n"); 
        fprintf(frest,"ZONE N= %d E= %d  F=FEPOINT  ET=BRICK\n", nshgtot, neltot);
        for(i=0; i< nshgtot; i++){
            for(j=0; j < 3; j++)
                fprintf(frest,"%e  ",xglobal[j*nshgtot+i]);
            for(k=0; k < numvar; k++)
                fprintf(frest,"%22.16e  ",qglobal[k*nshgtot+i]);
            fprintf(frest,"\n");
        }
 
        fprintf(frest,"\n");

        for(i=0; i< neltot; i++){
            for(j=0; j < nendx; j++)
                fprintf(frest,"%d  ",ien[j*neltot+i]);
            fprintf(frest,"\n");
        }
    } else if(nen == 4) {
        fprintf(frest,"Title = \"FE-VOLUME TETRAHEDRAL DATA SET\"\n");
        fprintf(frest,"VARIABLES = \"X\", \"Y\", \"Z\", \"P\", \"U\", \"V\", \"W\", \"T\"\n");
        fprintf(frest,"ZONE N= %d E= %d  F=FEPOINT  ET=TETRAHEDRON\n", nshgtot, neltot);
        for(i=0; i< nshgtot; i++){
            for(j=0; j < 3; j++)
                fprintf(frest,"%e  ",xglobal[j*nshgtot+i]);
            for(k=0; k < numvar; k++)
                fprintf(frest,"%22.16e  ",qglobal[k*nshgtot+i]);
            fprintf(frest,"\n");
        }

        fprintf(frest,"\n");

        for(i=0; i< neltot; i++){
            for(j=0; j < nendx; j++)
                fprintf(frest,"%d  ",ien[j*neltot+i]);
            fprintf(frest,"\n");
        }

    } else {
        cout << endl;
	cout << "Not yet implemented for this topology." << endl;
    }

    fclose(frest);

    return;
}
