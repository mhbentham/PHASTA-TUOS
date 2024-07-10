#include <vector>
#include <iostream>
#include "nspre_data.h"
#include "MeshSimInternal.h"

extern globalInfo* info;
extern pMeshDataId poly;

pMeshDataId shpfnNum;

int 
Entity_getNumDOF( pEntity ent ) {
    // returns the active degress of freedom on each entity
    dofInfo* dnfo=NULL;
    EN_getDataPtr( ent, shpfnNum, (void**) &dnfo );
    if( !dnfo ) return 0;
    else return (dnfo->dof).size();
}

void 
Entity_setDOFnumber( pEntity ent, int gdof ) {
    dofInfo* dnfo=NULL;
    EN_getDataPtr( ent, shpfnNum, (void**)&dnfo );
    if ( !dnfo ) { // we have to create one
        dnfo = new dofInfo;
        EN_attachDataPtr( ent, shpfnNum, ( void* )dnfo );
    }
    (dnfo->dof).push_back( gdof );
}

int 
FakeMode( globalInfo* info ) {
    
    if ( info->fake_mode == false )  {
        cout << "Introducing a fake mode " << endl;
        info->fake_mode=true;
        info->nshg++;
    }
    return ( info->nshg - 1 );
}

int 
Entity_getDOFnumber( pEntity ent , int ith ) {

    dofInfo* dnfo=NULL;
    EN_getDataPtr( ent, shpfnNum, (void**)&dnfo );

    if ( !dnfo ) return FakeMode( info );
    else if ( ith < (dnfo->dof).size() ) return dnfo->dof[ith];  
    else return FakeMode( info );
}

bool 
isActive( pEntity ent ) {
    dofInfo* dnfo=NULL;
    EN_getDataPtr( ent,  shpfnNum, (void**)&dnfo );
    return ( dnfo )? true : false ;
}

void 
Entity_addOffset( pEntity ent, int offset ) {
    
    dofInfo* dnfo = NULL;
    EN_getDataPtr( ent, shpfnNum, (void**)&dnfo );
    if( dnfo ) { 
        for( std::vector<int>::iterator iter = (dnfo->dof).begin() ;
             iter != (dnfo->dof).end();
             iter++ ) 
            *iter += offset;
    }
    return ;
}

// returns number of faces of an element of a particular block
// (declared in nspre_functions.h)
int 
NEFPerBlock(blockKey key)
{
    int nfaces;
 
    // number of faces:
    // wedge: 8 nodes --> 6 faces
    // tet  : 4 nodes --> 4 faces
    // key.nen/6 gives integer multiplicity of key.nen in 6
    nfaces = key.nen == 8 ? 6: (key.nen-(key.nen/6));

    return nfaces;

}
