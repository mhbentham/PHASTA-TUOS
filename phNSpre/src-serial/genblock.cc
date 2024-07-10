/* This routine generates the block structure for the mesh
   Fall 2000 AKK

   Rememmber to increase the allocation of NBblock from 12 to 18 when 
   we get either cubic wedges or pyramids working.

*/
#include <iostream>
#include <map>
#include <vector>
using namespace std;
#include <stdlib.h>
#include "nspre_functions.h"
#include "nspre_data.h"
#ifdef SIM
#include "MeshSim.h"
#else
#include "AOMD.h"
#endif

int* NIblock;
int* NBblock;
int* Nshape;
int* NshapeB;
extern globalInfo *info;
extern pMeshDataId poly;

forwardblock Iblock;
forwardblock Bblock;

// region interior block
blockKey *RIblock;
// region boundary block
blockKey *RBblock;

// taking blockid as index 
int* numInteriorFacesOfBlock;
int* numBoundaryFacesOfBlock;

namespace {
    int getnumIshape(blockKey bKey) {

        int lnen = bKey.nen;
        int poly = bKey.maxpoly;
        int nedg = lnen +(lnen+1)/2;
        int numshp=lnen;

        if( poly < 2 ){ return numshp;}
        else { numshp += nedg*(poly-1);}
 
        if (poly < 3 ) return numshp ; 
        switch(lnen){
        case 4:
        case 5:
            numshp += 4;
            break;
        case 6:
            numshp += 2;
            break;
        default:
            break; 
        }
        if ( poly > 3 ){ 
            cerr << "Code not ready for p > 3 " << endl;
            exit(-1);
        }
        return numshp;
    }
}
   
int 
getnumBshape(blockKey bkey) {

    int numEdges = bkey.nenbl;
    int poly = bkey.maxpoly;
    int numshp = numEdges;
    if( poly > 1) numshp += numEdges *( poly -1);
    if(( poly > 2) && (numEdges == 3))
        numshp ++;
    if ( poly > 3 ){ 
        cerr << "Code not ready for p > 3 " << endl;
        exit(-1);
    }
    return numshp;
}

/* 2D meshes */

void 
genblock2D( pMesh mesh, pGModel model, vector<pEdge>& bdry ) {

      int counter;

        /*  we use counter to give unique blockids to different blocks 
             but since we use 0 to identify nonexistent blocks we 
            should start counting from 1 */
        counter = 1;
        NIblock = new int [12];
        NBblock = new int [12];
        RIblock = new blockKey [ 12 ];
        RBblock = new blockKey [ 12 ];
        Nshape = new int [12];
        NshapeB = new int [12];
        for(int j=0; j< 12; j++) {
            NIblock[j] =0;
            NBblock[j] =0;
            Nshape[j] =0;
            NshapeB[j] =0;
            RIblock[j].nen = 0;
            RIblock[j].maxpoly = 0;
            RIblock[j].nenbl = 0;
            RIblock[j].lcsyst = 0;
            RBblock[j].nen = 0;
            RBblock[j].maxpoly = 0;
            RBblock[j].nenbl = 0;
            RBblock[j].lcsyst = 0;
        }

 
      FIter mfIter = M_faceIter(mesh);
      pFace face;
      blockKey bKey;
    /* First we gerenerate the interior element blocking structure */
    /* for this we loop over all the elements in the mesh and insert and
       update them into Iblock  */

      while (face = FIter_next(mfIter)){
        EN_getDataInt((pEntity)face, poly, &bKey.maxpoly);
        bKey.nen = numVertsF( face );
        bKey.nenbl = 2;
        bKey.lcsyst = topology2D(face);
    
        if (Iblock[bKey]) {        /* Already existing block */
            NIblock[Iblock[bKey]-1]++;
            /* increment the nbr of elements in */
            /* this block */
        } else {                        /* creation of a new block */      
            Iblock[bKey]=counter++;  /* assign a unique ascending id */
            NIblock[counter-2] = 1;    /* set num elements in block to 1 */
            Nshape[counter-2]  = getnumIshape(bKey);
            RIblock[counter-2] = bKey;
        }
      }

      FIter_delete(mfIter);
      counter = 1;  
      pGEdge gedge;
      pEdge edge;
      GEIter geIter = GM_edgeIter(model);
      info->numelb = 0;
  
      while( gedge = GEIter_next(geIter)){

#ifdef SIM
          if(!GEN_dataI((pGEntity)gedge,"NoBI")){    
#else
           int tmp;
          if(!GEN_dataI((pGEntity)gedge,"NoBI", &tmp)){    
#endif 
            /* If the face has boundary elements */

              EIter eIter = M_classifiedEdgeIter(mesh,(GEntity*) gedge,0);
              while (edge = EIter_next(eIter)) {
                void* rtmp = 0;
                pPList e_faces = E_faces(edge);
                face = (pFace)PList_next(e_faces,&rtmp);
                PList_delete(e_faces);
                EN_getDataInt((pEntity)face, poly, &bKey.maxpoly);
                bKey.nen = numVertsF( face );
                bKey.nenbl = 2;
                bKey.lcsyst = topology2D(face);

                if (Bblock[bKey]) {        /* Already existing block */
                    NBblock[Bblock[bKey]-1]++ ;
                    /* increment the nbr of elements in */
                    /* this block */
                } else {                       /* creation of a new block */      
                    Bblock[bKey] = counter++; /* assign a unique ascending id */
                    NBblock[counter-2] = 1; /* set num elements in block to 1 */
                    NshapeB[counter-2] = getnumIshape(bKey);
                    RBblock[counter-2] = bKey;
                }
                bdry.push_back( edge );
                (info->numelb)++;
              }
              EIter_delete(eIter);
          }
      }
      GEIter_delete(geIter);
}

/* 3D meshes */

void 
genblock( pMesh mesh, pGModel model, vector<pFace>& bdry ) {

    int counter;

        /*  we use counter to give unique blockids to different blocks 
             but since we use 0 to identify nonexistent blocks we 
            should start counting from 1 */
        counter = 1;
        NIblock = new int [12];
        NBblock = new int [12];
        RIblock = new blockKey [ 12 ];
        RBblock = new blockKey [ 12 ];
        Nshape = new int [12];
        NshapeB = new int [12];
        for(int j=0; j< 12; j++) {
            NIblock[j] =0;
            NBblock[j] =0;
            Nshape[j] =0;
            NshapeB[j] =0;
            RIblock[j].nen = 0;
            RIblock[j].maxpoly = 0;
            RIblock[j].nenbl = 0;
            RIblock[j].lcsyst = 0;
            RBblock[j].nen = 0;
            RBblock[j].maxpoly = 0;
            RBblock[j].nenbl = 0;
            RBblock[j].lcsyst = 0;
        }

 
    RIter mrIter = M_regionIter(mesh);
    pRegion region;
    blockKey bKey;
    /* First we gerenerate the interior element blocking structure */
    /* for this we loop over all the elements in the mesh and insert and
       update them into Iblock  */

    while (region = RIter_next(mrIter)){
        EN_getDataInt((pEntity)region, poly, &bKey.maxpoly);
        bKey.nen = numVertsR( region );
        bKey.nenbl = (bKey.nen == 8 ? 4 : 3);
        bKey.lcsyst = topology(region);
    
        if (Iblock[bKey]) {        /* Already existing block */
            NIblock[Iblock[bKey]-1]++;
            /* increment the nbr of elements in */
            /* this block */
        } else {                        /* creation of a new block */      
            Iblock[bKey]=counter++;  /* assign a unique ascending id */
            NIblock[counter-2] = 1;    /* set num elements in block to 1 */
            Nshape[counter-2]  = getnumIshape(bKey);
            RIblock[counter-2] = bKey;
        }
    }

    RIter_delete(mrIter);

    /* After the interior blocks, we create the boundary blocks. */
 
    counter = 1;  
    pGFace gface;
    pFace face;
    GFIter gfIter = GM_faceIter(model);
    info->numelb = 0;
  
    while( gface = GFIter_next(gfIter)){

#ifdef SIM
        if(!GEN_dataI((pGEntity)gface,"NoBI")){     
#else
        int tmp;
        if(!GEN_dataI((pGEntity)gface,"NoBI", &tmp)){     
#endif
            /* If the face has boundary elements */

              FIter fIter = M_classifiedFaceIter(mesh,(GEntity*) gface,0);
              while (face = FIter_next(fIter)) {
                void* rtmp = 0;
                pPList f_regions = F_regions(face);
                region = (pRegion)PList_next(f_regions,&rtmp);
                PList_delete(f_regions);
                EN_getDataInt((pEntity)region, poly, &bKey.maxpoly);
                bKey.nen = numVertsR( region );
                bKey.nenbl = F_numEdges(face);
                bKey.lcsyst = topology(region);

                /* here we try to seperate quad bface wedges from tri bface wedges*/

                if ( 3 == bKey.lcsyst) bKey.lcsyst = bKey.nenbl; 

                /* we also seperate quad face pyramids from triface pyramids */
                /* quad faced pyramids retain lcsyst = 5 and tri faced pyramids
                   have their lcsyst changed to 6 */

                if ( 5 == bKey.lcsyst) 
                    if ( 3 == bKey.nenbl ) bKey.lcsyst = 6;

                if (Bblock[bKey]) {        /* Already existing block */
                    NBblock[Bblock[bKey]-1]++ ;
                    /* increment the nbr of elements in */
                    /* this block */
                } else {                       /* creation of a new block */      
                    Bblock[bKey] = counter++; /* assign a unique ascending id */
                    NBblock[counter-2] = 1; /* set num elements in block to 1 */
                    NshapeB[counter-2] = getnumIshape(bKey);
                    RBblock[counter-2] = bKey;
                }
                bdry.push_back( face );
                (info->numelb)++;
              }
              FIter_delete(fIter);
        }
    }
    GFIter_delete(gfIter);
}
