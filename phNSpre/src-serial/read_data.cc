#include <iostream>
#include "phastaIO.h"

extern char* iformat;

using std::cerr;
using std::endl;

void 
restart( double* q, 
         int nshg, 
         int nvr, 
         int* lstep, 
         char filename[] ) {

	int restart;
  	openfile_( filename , "read",  &restart );

	int iarray[4];
	int isize = 3;

	readheader_( &restart, "solution", iarray,
                 &isize, "double", iformat );
    
	isize = iarray[0]*iarray[1];
    *lstep = iarray[2];
    double* qlocal = new double [ isize ];
    
    readdatablock_( &restart, "solution", qlocal, &isize,
                    "double" , iformat );

    // copy the needed part into the array passed in 

	if ( nshg > iarray[0] ) {
        cerr << "reading only " << iarray[0] << "modes from restart" << endl;
        cerr << nshg << " modes were requested " << endl;
        nshg = iarray[0];
    }
    if ( nvr  > iarray[1] ) {
        cerr << "reading only " << iarray[1] << "vars from restart" << endl;
        cerr << nvr << " modes were requested " << endl;
        nvr  = iarray[1];
    }
    for(int i=0; i < nvr; i++)
        for(int j=0; j < nshg; j++)
            q[i*nshg+j] = qlocal[i*iarray[0]+j];
    
    delete [] qlocal;
}
