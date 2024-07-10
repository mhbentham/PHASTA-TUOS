// This routine is being added in place of setup.c 
// This will read the P order distribution from a file or some other
// source and then set all the entity level information. 
//
// Anilkumar Karanam  ( Summer of 2000)
//

#include <iostream>
#include <fstream>
#include <stdlib.h>

#ifndef CYGWIN
#include <unistd.h>
#include <strings.h>
#endif


#include "nspre_data.h"
#include "nspre_functions.h"

using namespace std;
extern bool refine_mesh, c1quintic;
extern globalInfo* info;
extern pMeshDataId poly;
extern bool isReorder;
extern int numVertsR( pRegion region );
extern int numVertsF( pFace face );

namespace {
    double
    Xcoord( pRegion region ) {
        double xyz[3];
        double X=0.0;
        pPList verts = R_vertices( region, 1 );
        for(int z=0; z< PList_size( verts ); z++ ) {
           V_coord( (pVertex)PList_item( verts, z ), xyz );
           X += xyz[0];
        } 
        X /= PList_size( verts );
        PList_delete( verts );

        return X;
    }
}


int 
topology(pRegion rgn) {

    // This functions determines the topology if any given region
    // 1 : Tetrahedron
    // 2 : Hexahedron
    // 3 : Wedge
    // 5 : Pyramid 
 
    switch(R_numFaces(rgn)){
    case 4:          /* tetrahedron */
        return 1;
    case 6:          /* hehaxhedron */
        return 2;
    case 5:          /* could be either a Pyramid or a Wedge */
    {
        /* To resolve this issue we check the number of vertices on
           the first face of the region ... We have reached an
           agreement with the people of Mesh Database that a wedge
           always has its tri face first and a Pyramid always has its
           quad-face first ... Or was it the other way around ????
   
           In anycase we use the first option ! 

           since numEdges == numVertices for any closed face we check
           number of edges on the first face . */


        switch(F_numEdges(R_face(rgn,0))){
        case 3:
            return 3;
        case 4:
            return 5;
        default:
            cerr << "Congratulations : You have got a Mutant Face in your Mesh"
                 << endl 
                 << "It has neither 3 nor 4 vertices "
                 << endl
                 << "This preprocessor decided to give up Here " 
                 << endl
                 << "Best of Luck ;)" << endl;
            exit(1);
        }
    }
    break;

    default:
        cerr <<" Unidentified Topology in the Input mesh "<<endl;
        cerr <<" Exiting program " << endl;
        exit(1);
    }
    return 0;
}	

int 
topology2D(pFace face) {

    // This functions determines the topology for 2D pb
    // 1 : triangle
    // 2 : rectangle
 
    switch(F_numEdges(face)){
    case 3:          /* triangle */
        if (c1quintic == true) return 3;
	else return 1;
    case 4:          /* rectangle */
        return 2;
    default:
        cerr <<" Unidentified 2D Topology in the Input mesh "<<endl;
        cerr <<" Exiting program " << endl;
        exit(1);
    }
    return 0;
}	

void 
setup_and_refine( pGModel model, 
                  pMesh mesh,
                  globalInfo* info  ) {
  
    pVertex vertex;
    pEdge edge;
    pFace face;
    pRegion region;
    map< int, vector< int > > idshpmap;
    int id, nentm, nshg_old, itmp;
    info->nedgemodes = 0;
    info->nfacemodes = 0;
    info->fake_mode = false;
    info->nshg = 0;
  
    int nvs = M_numVertices ( mesh );
    int nes = M_numEdges ( mesh );

    RIter rIter = M_regionIter( mesh ) ;
    info->nenmax = 3;
    int nenl=3;
    while ( region = RIter_next( rIter ) ) {
       nenl = numVertsR( region );
        if ( nenl > info->nenmax  ) info->nenmax = nenl;
    }
    if (info->nsd == 2) {
       FIter fIter = M_faceIter( mesh ) ;
       while ( face = FIter_next( fIter ) ) {
         nenl =  numVertsF( face );
          if ( nenl > info->nenmax  ) info->nenmax = nenl;
       }
       FIter_reset( fIter );
    }
   
    VIter vIter = M_vertexIter(mesh);

    if(!isReorder){ // 
        while( vertex = VIter_next(vIter)) {
            Entity_setDOFnumber((pEntity)vertex,(info->nshg)++);
        }
    }else{// DOF number reordering
        V_reordering(mesh, info);
    }

    if (  refine_mesh ) {
        
        ifstream idmap("idmap.dat");
        
        if ( idmap.good() ) {
            idmap >> nshg_old;
            while ( !(idmap.eof()) ) {
                vector< int > list;
                idmap >> id;
                idmap >> nentm;
                for(int row=0; row < nentm ; row++ ){
                    idmap >> itmp;
                    list.push_back( itmp); 
                }
                idshpmap[ id ] = list ;
                list.clear();
            }
        }
        idmap.close();
    }

    int PolyOrd=1;
    int nem;    /* number of modes on the entity */
    EIter eIter = M_edgeIter(mesh);
    while( edge = EIter_next(eIter)) { 
        nem = 0;
        if( refine_mesh ) {
            nem =  idshpmap[ EN_id( (pEntity)edge ) + nvs ].size();
            if ( nem > 0 ) {
                EN_modifyDataInt( (pEntity) edge, poly, nem+1 );
                fixRegions( (pEntity) edge );
            }
        }
        EN_getDataInt( (pEntity) edge, poly, &PolyOrd );
        if ( PolyOrd > 1 ) {
            nem = PolyOrd - 1 ;
            for(int z=0; z< nem ; z++) {
                Entity_setDOFnumber((pEntity)edge,(info->nedgemodes)++); 
                (info->nshg)++;
            }
        }
    }
  
    PolyOrd=1;
    FIter fIter = M_faceIter(mesh);
    while(face = FIter_next(fIter)) { 
        nem = 0;
        if( refine_mesh ) {
            nem =  idshpmap[ EN_id( (pEntity)face ) + nvs + nes ].size();
            if ( nem > 0 ) {
                PolyOrd = nem - 1 + F_numEdges( face );
                EN_modifyDataInt( (pEntity) face, poly, PolyOrd );
                fixRegions( (pEntity) face );
            }
        }
        EN_getDataInt( (pEntity) face, poly, &PolyOrd );
        if ( PolyOrd > 2 ) {
            if ( PolyOrd > 3 ) {
                cerr << "code not ready " << __FILE__ << __LINE__ << endl;
                exit( 1 );
            }
            nem = PolyOrd - F_numEdges( face ) + 1;
            for( int z=0; z< nem; z++){
                 Entity_setDOFnumber( (pEntity)face, (info->nfacemodes)++ );
                 (info->nshg)++;
            }
        }
    }

    /* this is the place to add new modes if we decide to refine */
  
     if ( refine_mesh ) 
         prefine_on_error( mesh, info, idshpmap, nshg_old , "errorc.inp" );

    /* addiing the fake equation number for all the dormant modes */
    /* now we introduce the fake degree of freedom only when needed and it is
        done automatically */
    /* info->nshg++; */

    /* now we generate the mapping between entity ids and shape fn numbers */
    /* while we do this, we also correct the shape function numbrs with the correct */
    /* entitiy offsets  */
    /* All refinement has also been done now, so this is the best place to do it  */

    VIter_reset( vIter );
    EIter_reset( eIter );
    FIter_reset( fIter );
    RIter_reset( rIter );

    unlink( "idmap.dat" );
    ofstream idmap( "idmap.dat" );
    idmap <<  ( ( info->fake_mode ) ? info->nshg - 1 : info->nshg ) << endl;
    while ( vertex = VIter_next ( vIter ) )  
        idmap << EN_id((pEntity)vertex)<<" 1 "
              << Entity_getDOFnumber( (pEntity) vertex, 0 )<< endl;

    VIter_delete(vIter);

    while ( edge = EIter_next ( eIter ) ) {
        int nem = Entity_getNumDOF( (pEntity)edge );
        if ( nem > 0 ) {
            Entity_addOffset( (pEntity)edge, nvs );
            idmap << EN_id( (pEntity)edge ) + nvs <<" " ;
            idmap << nem <<" ";
            for( int em=0; em < nem ; em++ ) 
                idmap << Entity_getDOFnumber( (pEntity)edge, em)<<" ";
            idmap << endl; 
        }
    }

    EIter_delete(eIter);

    while ( face = FIter_next ( fIter ) ) {
        int nfm = Entity_getNumDOF( (pEntity)face );
        if ( nfm > 0 ) {
            Entity_addOffset( (pEntity)face , ( nvs + info->nedgemodes) );
            idmap << EN_id( (pEntity)face ) + nvs + nes <<" " ;
            idmap << nfm <<" ";
            for( int fm=0; fm < nfm ; fm++ ) 
                idmap << Entity_getDOFnumber( (pEntity)face, fm)<<" ";
            idmap << endl; 
        }
    }
  
    FIter_delete(fIter);
    idmap.close();

    // something to visualize the refinement 

    ofstream refviz("partit.out");
    while( region = RIter_next( rIter) ) {
        PolyOrd = 1 ;
        EN_getDataInt( (pEntity) region, poly, &PolyOrd );
        refviz << PolyOrd  << endl;
    }
    refviz.close();
    RIter_delete( rIter );
}
