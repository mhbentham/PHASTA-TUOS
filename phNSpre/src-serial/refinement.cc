#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <math.h>
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#endif
#include "nspre_data.h"
#include "nspre_functions.h"

using namespace std;
extern pMeshDataId poly;
extern bool refine_mesh;
const int n_err_measures = 3 ; 
const int wgt[]={1,1,1,0,0,0};

void
fixRegions( pEntity ent ) {

    int type;
    type = EN_type( ent );

    pPList regions ;
    
    switch( type ) {
    case Tedge:
        regions = E_regions( (pEdge)ent );
        break;
    case Tface:
        regions = F_regions( (pFace)ent );
        break;
    default:
        cerr << "unknown topology" << __FILE__<<":"<< __LINE__ << endl;
        break;
    }
    pRegion reg ;
    int pe,pr;
    
    // increment the polynomial order of this ent

    EN_getDataInt( ent, poly, &pe );
        
    for( int i = 0 ; i < PList_size( regions ) ; i++ ) {
        reg = (pRegion)PList_item( regions, i ) ;
        EN_getDataInt((pEntity) reg, poly, &pr);
        if ( pr < pe ) EN_modifyDataInt( (pEntity) reg, poly, pe ); 
    }

    if ( type == Tedge ) {
        double xyz0[3];
        double xyz1[3];
        
        V_coord( E_vertex( (pEdge)ent, 0 ), xyz0 );
        V_coord( E_vertex( (pEdge)ent, 1 ), xyz1 );
        
        cout << "ent " << EN_id( ent ) << " " 
             << xyz0[0] << " "<< xyz0[1] << " " << xyz0[2] << " to "
             << xyz1[0] << " "<< xyz1[1] << " " << xyz1[2] << endl;
    }
}

void
prefine_on_error( pMesh mesh,
                  globalInfo* info,
                  map< int, vector< int > > idshpmap,
                  int nshg_old,
                  char filename[] ) { 
  

    // first we will create the id to shp map

    char inpLine[125];
    char* pch;

    ifstream chkfp( filename );
    if ( chkfp.bad()) return ;
    chkfp.close();

    double ethreshold = 100000.00;
  
    // we read in the reduced error and create an rms data structure 
    // we also combine the error based on a given set of weights. [ right now all 1 ]
  
    int junk;
    double* error = new double [ nshg_old * n_err_measures ];
    restart ( error, nshg_old, n_err_measures , &junk,  filename );

    double* error_measure = new double [ nshg_old ];

    for( int shp=0; shp < nshg_old ; shp++ ) {
        error_measure[ shp ] = 0.0;
        for( int err=0; err < n_err_measures ; err++ )
            error_measure[ shp ] += fabs(error[ err*nshg_old + shp ]) * wgt[ err ];      
        error_measure[ shp ] /= n_err_measures;
    }   

    delete [] error ;

    // now we map a single error value is mapped to each entity, this is done
    // by taking the average of all the  error modes on it.

    map< int, double > iderrormap;
    vector< double > medianlist;
    int count1 = 0;
    double minerror= 10000.0;
    double maxerror= 0.0;
    for( map<int, vector<int> >::iterator miter = idshpmap.begin() ;
         miter != idshpmap.end();
         miter++ ) {

        for( vector<int>::iterator viter = (*miter).second.begin();
             viter != (*miter).second.end();
             viter++ ) {
            iderrormap[ (*miter).first ] += error_measure[ *viter ];
        }

        //iderrormap[ (*miter).first ] /= (*miter).second.size();
        if ( iderrormap[ (*miter).first ]  > maxerror ) 
            maxerror = iderrormap[ (*miter).first ];
        if ( iderrormap[ (*miter).first ]  < minerror ) 
            minerror = iderrormap[ (*miter).first ];
        medianlist.push_back( iderrormap[ (*miter).first ] );
    }

    delete [] error_measure;
    idshpmap.clear();

    cout << "Maximum Error is : " << maxerror << endl;
    cout << "Minimum Error is : " << minerror << endl;

    sort( medianlist.begin(), medianlist.end() );
    int middle = 1+medianlist.size()/2;
    double median_error = medianlist[ middle ] ;
    if( medianlist.size() % 2 == 0 ) {
        median_error += medianlist[ middle -1 ];
        median_error /= 2;
    }
    medianlist.clear();
    cout << "Median Error is : " << median_error << endl ;

    // now we refine.

    // cout << "Enter the Error Threshold to Refine: ";
    // cin >> ethreshold;
    // we will just refine edges for now .
    
    ethreshold = median_error;

    map< pEdge, int> refined;
    VIter  vIter = M_vertexIter( mesh );
    pVertex vertex;
    pEdge medge;
    bool refinement=false;
    while ( vertex = VIter_next( vIter ) ) {
        if ( iderrormap[ EN_id( (pEntity)vertex ) ] >= ethreshold ) {
            int nedges = V_numEdges( vertex );
            for( int e=0; e< nedges; e++ ) {
                pEdge edge = V_edge( vertex, e );  
                pVertex ov = E_otherVertex( edge, vertex );
                if ( ( iderrormap[ EN_id( (pEntity)ov ) ] > ethreshold )  
                     && !refined[ edge ] ) {

                    refinement = true;
                    Entity_setDOFnumber( (pEntity)edge, (info->nedgemodes)++ ) ;
                    (info->nshg)++;
                    refined[ edge ] = 1 ;
                    /* increment the polynomial order for this edge */
                    int p;
                    EN_getDataInt( (pEntity)edge, poly, &p );
                    EN_modifyDataInt( (pEntity)edge, poly, p+1 );
                    
                    /* update the surrounding regions */
                    fixRegions( (pEntity)edge );

                    if ( ( medge = (pEdge)EN_dataP((pEntity) edge,"PerM") ) 
                         && !refined[ medge ] ) {
                        Entity_setDOFnumber( (pEntity)medge, (info->nedgemodes)++ );
                        (info->nshg)++;
                        refined[ medge ] = 1 ;
                        EN_getDataInt( (pEntity)medge, poly, &p );
                        EN_modifyDataInt( (pEntity)medge, poly, p+1 );
                        fixRegions( (pEntity)medge );
                    }
                }
            }
        }
    }
              
    if ( refinement ) {
        cerr << endl << "Refinement has taken place " << endl;
        FakeMode( info );
    }
    refined.clear();
    iderrormap.clear(); 
}
