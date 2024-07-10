// the functions defined in this file, generate a mapping between the
// EN_id in the sms files to the mode numbers in the geombc.dat files.
// for vertices this is always one-to-one , so we don't bother to write. 

// structure of the written array is ( 1D array , as with anything
// written by phastaIO)
// en_id nem m1 m2 m3.....mnem en_id nem m1 m2 m3.....mnem

#include <iostream>
#include <vector>
#include "nspre_data.h"
#include "nspre_functions.h"
#include "phastaIO.h"

using namespace std;
extern char* oformat;
extern bool refine_mesh;

void 
writeSMS2PHASTA(  pMesh mesh, 
                  globalInfo* info  ) {

    // if there are no higher order modes or if the refinement is
    // uniform , we can safely return out of the funtion without writing
    // anything, since the mapping is simple and can be generated on the
    // fly.
    // NOTE: only real modes are written to the file
    
    if ( 0 == info->nedgemodes || !refine_mesh )  return ;
    int gfile=0;
    openfile_( "geombc.dat.1", "append", &gfile );

    if( !gfile ) { 
        cerr << "cannot open file "<< __FILE__ <<" : "<<__LINE__ << endl;
        return ;
    }

    int nem;
    vector<int> sms2phastaMAP; 

    EIter eIter = M_edgeIter( mesh );
    while ( pEdge edge = EIter_next( eIter ) ) {
        nem = Entity_getNumDOF( (pEntity)edge );
        if ( nem > 0 ) { 
            sms2phastaMAP.push_back( EN_id( (pEntity)edge ) );
            sms2phastaMAP.push_back( nem );
            for( int em=0; em < nem ; em++ )
                sms2phastaMAP.push_back( Entity_getDOFnumber( (pEntity)edge, em) );
        }
    }
    EIter_delete( eIter );

    if ( info->nfacemodes ) {
        FIter fIter = M_faceIter( mesh );
        while ( pFace face = FIter_next( fIter ) ) {
            nem = Entity_getNumDOF( (pEntity)face );
            if ( nem > 0 ) { 
                sms2phastaMAP.push_back( EN_id( (pEntity)face ) );
                sms2phastaMAP.push_back( nem );
                for( int em=0; em < nem ; em++ )
                    sms2phastaMAP.push_back( Entity_getDOFnumber( (pEntity)face, em) );
            }
        }
        FIter_delete( fIter );
    }

    // convert the vector to an integer array, so that we can easily
    // write to file.
    int size = sms2phastaMAP.size();
    int* sms2phasta = new int [ size ];
    for( int i=0; i < size; i ++ ) sms2phasta[i] = sms2phastaMAP[i];
    sms2phastaMAP.clear();

    int ifour=4;
    int iarray[4];
                    
    iarray[0] = size ;
    iarray[1] = info->nshg ;
    iarray[2] = info->nedgemodes ;
    iarray[3] = info->nfacemodes ;

    writeheader_( &gfile, "higher order refinement mapping ",
                  (void*)iarray, &ifour,  &size, "integer", oformat );

    writedatablock_( &gfile, "higher order refinement mapping ",
                     (void*)sms2phasta,  &size, "integer", oformat );

    closefile_(&gfile,"append");
    delete [] sms2phasta;
}
