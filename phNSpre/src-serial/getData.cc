/* 

This file contains functions which collect information from a mesh
into arrays needed by PHASTA

*/

#include <stdlib.h>  
#include <iostream>
#include <map>
#include "nspre_data.h"
#include "nspre_functions.h"

extern forwardblock Iblock;
extern forwardblock Bblock;
extern int* NBblock;
extern int* NIblock;
extern int* Nshape;
extern int* NshapeB;
extern globalInfo* info;
extern bool isReorder;

// taking blockid as index 
// defined in genblock.cc
extern int* numInteriorFacesOfBlock;
extern int* numBoundaryFacesOfBlock;

extern pMeshDataId poly;
extern pMeshDataId ReorderR;

int modeSign(int type, int pdof);
int Fsign(pFace face,pVertex v[3]);
int h_modeSign(int type, int pdof);

/**********************************************************************/
/* compute the coordinate array                                       */
/**********************************************************************/
void 
getX(pMesh mesh, double *x ) {

    pVertex vertex;
    dArray xloc;
    int dof;
    int numnp = M_numVertices( mesh );
    VIter vIter = M_vertexIter(mesh); 
    while( vertex = VIter_next ( vIter )  ){
        dof = Entity_getDOFnumber((pEntity )vertex, 0);
        V_coord(vertex,xloc);
        x[dof] = xloc[0];
        x[dof+numnp] = xloc[1];
        if (info->nsd == 3) x[dof+numnp*2] = xloc[2];
    }
    VIter_delete(vIter);
}

/************************************************************/
/* compute the interior and boundary element connectivity   */
/* 2D meshes: faces are elements                            */
/*                                                          */
/* only implemented for linear case                         */
/************************************************************/
void 
getConnectivity2D( pMesh mesh,
                 vector<pEdge> *bdry,
                 int** ien, 
                 int** ien_sms,
                 int** ienb ,
                 int** ienb_sms ) {


    pVertex *Fvertices;
    pFace face;
    pEdge edge;
    pPList ents;
    int i,ldof,nen;
    int blockid, porder;
    blockKey BLOCK;
    std::map<int,int> iel ;
    std::map<int,int> ielb; 

    /******************** interior elements ********************/

    Fvertices = (pVertex*)malloc(4*sizeof(pVertex)); 
    FIter fIter = M_faceIter(mesh);
    while (face = FIter_next(fIter) ){
        EN_getDataInt((pEntity) face, poly, &porder);
        nen =  numVertsF( face );
        BLOCK.nen = nen;
        BLOCK.maxpoly = porder;
        BLOCK.nenbl = 2;
        BLOCK.lcsyst = topology2D(face);
        blockid = Iblock[BLOCK]-1;
	
	ents = F_vertices(face,1);
        /* collect vertex mode numbers */
        for (i=0; i < nen; i++) {
	    Fvertices[i] = (pVertex)PList_item(ents,i);
            ien[blockid][iel[blockid]+i*NIblock[blockid]] = 
                Entity_getDOFnumber((pEntity )Fvertices[i],0);
	}
	PList_delete(ents);

        ldof = nen;
        ien_sms[ blockid ][ iel[ blockid ] ] = EN_id( (pEntity)face );
        iel[blockid]++; /* increment the element counter for this block on 
                           this processor */

    }
    FIter_delete(fIter);

    /******************** boundary elements ********************/

    vector< pEdge >::iterator bedgeIter = bdry->begin();
    while ( bedgeIter != bdry->end() ) {
        edge = *bedgeIter;
        /* find the list of faces connected to this edge. Since this
           is a boundary edge, there should only be one such face. */
        ents = E_faces(edge);
        face = (pFace)PList_item(ents,0); /* get the region */
        PList_delete(ents);
        
        nen = numVertsF( face );
        EN_getDataInt((pEntity) face, poly, &porder);
        BLOCK.nen = nen;
        BLOCK.maxpoly = porder;
        BLOCK.nenbl = 2;
        BLOCK.lcsyst = topology2D(face);
        
        blockid = Bblock[BLOCK]-1;
      
        /* get the vertices */
	/* only works for triangles */
        ents = F_vertices(face,1);
	for (i=0; i < 2; i++)
	    Fvertices[i] = E_vertex(edge,i);

        for (i=0; i < nen; i++){
	    if ( ((pVertex)PList_item(ents,i) != Fvertices[0]) &&
	         ((pVertex)PList_item(ents,i) != Fvertices[1]) ) {
	      Fvertices[2] = (pVertex)PList_item(ents,i);
	    }  
	}

        for (i=0; i < nen; i++){
            ienb[blockid][ielb[blockid]+NBblock[blockid]*i] = 
                Entity_getDOFnumber((pEntity )Fvertices[i],0);
        }

        ldof = nen;
        ienb_sms[ blockid ][ ielb[ blockid ] ] = EN_id( (pEntity)face);
        ielb[blockid]++;
        bedgeIter++;
    }
    ielb.clear();

    free(Fvertices);
}
/************************************************************/
/* compute the interior and boundary element connectivity   */
/* 3D meshes: regions are elements                          */
/************************************************************/
void 
getConnectivity( pMesh mesh,
                 vector<pFace> *bdry,
                 int** ien, 
                 int** ien_sms,
                 int** ienb ,
                 int** ienb_sms ) {

    /* Connectivity Templates Begin */ 
    
    // Needed for polynomial order greather then 1
  
    // ordering of edges and faces for tets
    int TetEMap[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
    int TetFMap[4][4]={{0,1,2,-1},{0,3,1,-1},{1,3,2,-1},{0,2,3,-1}};
  
    
    // ordering of edges and faces for hexes
    int HexEMap[12][2]={{0,1},{1,2},{3,2},{0,3},{4,5},{5,6},{7,6},{4,7},
                        {0,4},{1,5},{2,6},{3,7}};
    int HexFMap[6][4]={{0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}, 
                       {4,5,6,7}}; 
    int HexfaceDir[6]={0,0,1,0,0,1};
    
    // ordering of edges and faces for wedges
    int WedEMap[9][2]={{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5}};
    int WedFMap[5][4]={{0,2,1,-1},{0,1,4,3},{1,2,5,4},{2,0,3,5},{3,4,5,-1}};
    int WedfaceDir[6]={0,0,0,0,0,-1};
  
    // ordering of edges and faces for pyramids
    int PyrEMap[8][2]={{0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4}};
    int PyrFMap[5][4]={{0,1,2,3},{1,0,4,-1},{2,1,4,-1},{3,2,4,-1},{0,3,4,-1}};
    int PyrfaceDir[6]={1,1,1,1,1,-1};

    /* Connectivity Templates End */ 

    int** GenericEdgeMap;
    int** GenericFaceMap;
    int*  GenericfaceDir;
 
    pVertex *Rvertices;
    pEdge   *Redges; 
    pFace   *Rfaces; 
  
    pRegion region;
    pPList ents;
    pVertex v1,v[3];
    pFace face;
    int i,j,nem;
    int ldof,esign,fsign,msign,k,nsh,is,ip;
    int PolyOrd;
    int nen,nfaces, nedges;
    int nvf ; /* number of vertices on this face */
    int blockid, porder;

    blockKey BLOCK;
    std::map<int,int> iel ;
    std::map<int,int> ielb; 
    
    /******************** interior elements ********************/
  
    GenericEdgeMap = (int **)malloc(12*sizeof(int *));
    for(i=0; i< 12; i++) GenericEdgeMap[i]=(int *)malloc(2*sizeof(int));
	
    GenericFaceMap = (int **)malloc(6*sizeof(int *));
    for(i=0; i< 6; i++) GenericFaceMap[i]=(int *)malloc(4*sizeof(int));

    GenericfaceDir = (int *)malloc(6*sizeof(int));

    /* Using the maximum values possible for verts, edges and faces */

    Rvertices = (pVertex*)malloc(8*sizeof(pVertex)); 
    Redges = (pEdge*)malloc(12*sizeof(pEdge));     
    Rfaces = (pFace*)malloc(6*sizeof(pFace));     

    RIter rIter = M_regionIter(mesh);
    while (region = RIter_next(rIter) ){
        /* processor for this region */
        EN_getDataInt((pEntity) region, poly, &porder);
        nen =  numVertsR( region );
        BLOCK.nen = nen;
        BLOCK.maxpoly = porder;
        BLOCK.nenbl = nen == 8 ? 4 : 3;
        BLOCK.lcsyst = topology(region);
        blockid = Iblock[BLOCK]-1;

        nedges = nen+(nen+1)/2;
        nfaces = nen == 8 ? 6: (nen-(nen/6));

        /* get all mesh entities connected to this region */
        R_entities(region, Rvertices, Redges, Rfaces,nen);

        /* collect vertex mode numbers */

        if(isReorder){
            int Rlable;
            EN_getDataInt(region,ReorderR,&Rlable);
            iel[blockid] = Rlable;
        }

        for (i=0; i < nen; i++)
            ien[blockid][iel[blockid]+i*NIblock[blockid]] = 
                Entity_getDOFnumber((pEntity )Rvertices[i],0);

        ldof = nen;
   
        if( porder > 1 ) { /* we will never have edges on a linear region */ 
            /* Before we start on the higher order modes we should point the
               topology connectivity templates to the right arrays based on
               the nen ( uniquely defines topology ) */


            switch(nen){
            case 4:   /* Tetrahedra */

                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = TetEMap[i][0];
                    GenericEdgeMap[i][1] = TetEMap[i][1];
                }
         
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = TetFMap[i][j];

                break;
            case 5:   /* Pyramid */
        
                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = PyrEMap[i][0];
                    GenericEdgeMap[i][1] = PyrEMap[i][1];
                }
         
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = PyrFMap[i][j];

                for(i =0; i< nfaces; i++)
                    GenericfaceDir[i] = PyrfaceDir[i];
        
                break;
        
            case 6:  /* Wedge Element */
        
                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = WedEMap[i][0];
                    GenericEdgeMap[i][1] = WedEMap[i][1];
                }
         
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = WedFMap[i][j];

                for(i =0; i< nfaces; i++)
                    GenericfaceDir[i] = WedfaceDir[i];
        
                break;
        
            case 8:  /* Hexahedron */

                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = HexEMap[i][0];
                    GenericEdgeMap[i][1] = HexEMap[i][1];
                }
         
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = HexFMap[i][j];

                for(i =0; i< nfaces; i++)
                    GenericfaceDir[i] = HexfaceDir[i];
        
                break;
            }
      
            /* collect edge mode numbers */
            /* find the direction the edge is being used by the region.
               If the first vertex in the edge map for this edge coincides
               with the first vertex of the edge (as it was defined),
               then the sign is 1, otherwise it is 0 */

            nem = porder-1;
      
            for ( i=0; i < nedges; i++){
                v1 = Rvertices[GenericEdgeMap[i][0]];
                esign = (v1 == E_vertex(Redges[i],0));
                for ( j=0; j < nem; j++){
                    msign = h_modeSign(1,j+2);
                    ien[blockid][iel[blockid]+NIblock[blockid]*ldof++] = 
                        Entity_getDOFnumber((pEntity )Redges[i],j);
                    if((esign == 0 ) &&(msign ==0 ))
                        ien[blockid][iel[blockid]+NIblock[blockid]*(ldof-1)] *= -1;
                }
            }

            if ( porder > 2 ) { // I have to figure something out for cubic hexes
    
                /* collecting  face mode numbers */
                /* Note : Code Reuse
   
                Reagardless of element topology only 2 types of faces exist ,
                tri/quad so we can write just 2 cases and use them for all
                topologies. Ofcourse the connectivity template is different for
                each template but that is taken care of using the switch
                statement above .

                In any case it is useless to add edge+triface+quadface
                specifically for each topology we encounter */

     

                for(i =0; i < nfaces; i++){

                    j =0;
                    nvf = F_numEdges(Rfaces[i]);
                    PolyOrd =porder;

                    /* Triangular Face */
                    if (nvf == 3){ 
                        /* find the direction the face is being used */
                        for (k=0; k < 3; k++)
                            v[k] = Rvertices[GenericFaceMap[i][k]];
                        fsign = Fsign(Rfaces[i],v);
                        for (ip = 3; ip <= PolyOrd; ip++){
                            msign = modeSign(2,ip);
                            /* get the number of shape functions of each order */
                            nsh = ip-2;
                            for (is=0; is < nsh; is++){
                                ien[blockid][iel[blockid]+NIblock[blockid]*ldof++] = 
                                    Entity_getDOFnumber((pEntity )Rfaces[i],j++);
                                if ((fsign == 0) && (msign == 0))
                                    ien[blockid][iel[blockid]+NIblock[blockid]*(ldof-1)] *= -1;
                            }
                        }
                        /* Quadrilateral Face */
                    } else if ( nvf == 4 ) {
                        /* find the direction in which the face is being used */
                        fsign = R_dirUsingFace(region,Rfaces[0]);
                        /* A zero here usually means that the face is defined pointing
                           into this particualar region */
                        for(ip = 4; ip <= PolyOrd; ip++) {
                            msign = h_modeSign(2,ip);
                            nsh = ip-3; /* number of shapefunctions of each order */
                            for( is =0; is< nsh; is ++) {
                                ien[blockid][iel[blockid]+NIblock[blockid]*ldof++] = 
                                    Entity_getDOFnumber((pEntity )Rfaces[0],j++);
                                if(fsign != GenericfaceDir[0] && msign == 0){
                                    ien[blockid][iel[blockid]+NIblock[blockid]*(ldof-1)] *= -1;
                                }
                            }
                        }
                    }
                }
            }  // p > 2 
        }  // p > 1
    
        /* ien[blockid][iel[blockid]+NIblock[blockid]*ldof] = EN_id( (pEntity)region); */
        ien_sms[ blockid ][ iel[ blockid ] ] = EN_id( (pEntity)region );
        iel[blockid]++; /* increment the element counter for this block on 
                           this processor */

    }
    RIter_delete(rIter);

    /******************** boundary elements ********************/

    vector< pFace >::iterator bfaceIter = bdry->begin();
    while ( bfaceIter != bdry->end() ) {
        face = *bfaceIter;
        /* find the list of regions connected to this face. Since this
           is a boundary face, there should only be one such region. */
        ents = F_regions(face);
        region = (pRegion)PList_item(ents,0); /* get the region */
        PList_delete(ents);
        
        nen = numVertsR( region );
        EN_getDataInt((pEntity) region, poly, &porder);
        BLOCK.nen = nen;
        BLOCK.maxpoly = porder;
        BLOCK.nenbl = F_numEdges(face);
        BLOCK.lcsyst = topology(region);
        if ( BLOCK.lcsyst == 3 ) BLOCK.lcsyst = BLOCK.nenbl;
        if ( BLOCK.lcsyst == 5 )
            if ( BLOCK.nenbl == 3 ) BLOCK.lcsyst = 6;
        
        blockid = Bblock[BLOCK]-1;
        nedges = nen+(nen+1)/2; 
        nfaces = nen == 8 ? 6: (nen-(nen/6));
      
        /* get the mesh bounding mesh entities */
        R_entitiesBdry(region,face,Rvertices,Redges,Rfaces,nen);
      
        /*********** gather equation numbers ***********/
        /* vertices */

        for (i=0; i < nen; i++){
            ienb[blockid][ielb[blockid]+NBblock[blockid]*i] = 
                Entity_getDOFnumber((pEntity )Rvertices[i],0);
        }

        ldof = nen;

        if ( porder > 1 ) {

            /* All the stuff below is same as for interior elements , just
               look the same comments... someday should find elegant way of
               reusing code for this */

   
            switch(nen){
            case 4:   /* Tetrahedra */
          
                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = TetEMap[i][0];
                    GenericEdgeMap[i][1] = TetEMap[i][1];
                }
          
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = TetFMap[i][j];

                break;
            case 5:   /* Pyramid */
          
                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = PyrEMap[i][0];
                    GenericEdgeMap[i][1] = PyrEMap[i][1];
                }
          
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = PyrFMap[i][j];
          
                for(i =0; i< nfaces; i++)
                    GenericfaceDir[i] = PyrfaceDir[i];
          
                break;
          
            case 6:  /* Wedge Element */
          
                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = WedEMap[i][0];
                    GenericEdgeMap[i][1] = WedEMap[i][1];
                }
          
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = WedFMap[i][j];
          
                for(i =0; i< nfaces; i++)
                    GenericfaceDir[i] = WedfaceDir[i];
          
                break;
          
            case 8:  /* Hexahedron */
          
                for(i=0; i< nedges; i++){
                    GenericEdgeMap[i][0] = HexEMap[i][0];
                    GenericEdgeMap[i][1] = HexEMap[i][1];
                }
          
                for(i=0; i< nfaces; i++) 
                    for(j=0; j<4;  j++)
                        GenericFaceMap[i][j] = HexFMap[i][j];
          
                for(i =0; i< nfaces; i++)
                    GenericfaceDir[i] = HexfaceDir[i];
          
                break;
            }

            nem = porder-1;
            for ( i=0; i < nedges; i++){
                v1 = Rvertices[GenericEdgeMap[i][0]];
                esign = (v1 == E_vertex(Redges[i],0));
                for ( j=0; j < nem; j++){
                    msign = h_modeSign(1,j+2);
                    ienb[blockid][ielb[blockid]+NBblock[blockid]*ldof++] =
                        Entity_getDOFnumber((pEntity )Redges[i],j);
                    if((esign == 0 ) &&(msign ==0 ))
                        ienb[blockid][ielb[blockid]+NBblock[blockid]*(ldof-1)] *= -1;
                }
            }
      
            /* look at the detailed comment for interior elements */

            if ( porder > 2 ) {
        
                for(i =0; i < nfaces; i++){
                    j =0;
                    nvf = F_numEdges(Rfaces[i]);
                    PolyOrd = porder;
                    if (nvf == 3){
                        for (k=0; k < 3; k++)
                            v[k] = Rvertices[GenericFaceMap[i][k]];
                        fsign = Fsign(Rfaces[i],v);
                        for (ip = 3; ip <= PolyOrd; ip++){
                            msign = modeSign(2,ip);
                            nsh = ip-2;
                            for (is=0; is < nsh; is++){
                                ienb[blockid][ielb[blockid]+NBblock[blockid]*ldof++] = 
                                    Entity_getDOFnumber((pEntity )Rfaces[i],j++);
                                if ((fsign == 0) && (msign == 0))
                                    ienb[blockid][ielb[blockid]+NBblock[blockid]*(ldof-1)] *= -1;
                            }
                        }
                    } else if( nvf == 4 ){
                        fsign = R_dirUsingFace(region,Rfaces[0]);
                        for(ip = 4; ip <= PolyOrd; ip++) {
                            msign = h_modeSign(2,ip);
                            nsh = ip-3; /* number of shapefunctions of each order */
                            for( is =0; is< nsh; is ++) {
                                ienb[blockid][ielb[blockid]+NBblock[blockid]*ldof++] = 
                                    Entity_getDOFnumber((pEntity )Rfaces[0],j++);
                                if(fsign != GenericfaceDir[0] && msign == 0){
                                    ienb[blockid][ielb[blockid]+NBblock[blockid]*(ldof-1)] *= -1;
                                }
                            }
                        }
                    }
                }
            } // p > 2
        } // p > 1 
      
        /* ienb[blockid][ielb[blockid]+NBblock[blockid]*ldof] = EN_id( (pEntity)region); */
        ienb_sms[ blockid ][ ielb[ blockid ] ] = EN_id( (pEntity)region);
        ielb[blockid]++;
        bfaceIter++;
    }
    ielb.clear();
  
    /* free the counter array */
    for(i=0; i< 6; i++) free(GenericFaceMap[i]);
    for(i=0; i< 12; i++) free(GenericEdgeMap[i]);
    free(GenericfaceDir);
    free(GenericEdgeMap);
    free(GenericFaceMap);
    free(Rvertices);
    free(Redges);
    free(Rfaces);
}
// ONLY WORKS FOR 3D
////////////////////////////////////////////////////////////////////////////////////////////////////
// obtain face connectivity i.e a map from local to global face numbering
// later -- to be used as phasta array ief, iefb (interior, and boundary
// elements respectively)
// called by: 
// writeEnsaFiles in writeEnsaFiles.cc
//
// each block is uniquiely identified by a blockid and that is used to sort
// ief,  ien, too.
// iel[blockid] is the current running counter on the number of interior elements is a particular block.
// ielb[blockid] is the current running counter on the number of boundary elements is a particular block.
///////////////////////////////////////////////////////////////////////////////////////////////////
void 
getFaceConnectivity(pMesh mesh, 
                    vector<pFace> bdry, 
                    int **ief,
                    int **iefb)

{
//                    int **ief_sms,                    
//                    int **iefb_sms)    
    int i,k;
    int nen,nfaces, nedges;
    int blockid, porder;
    void* Iter;

    int numberOfBlocks = Iblock.size();
    int numberOfBoundaryBlocks = Bblock.size(); 

    pPList* faceLists = new pPList[numberOfBlocks];


    for (i=0;i<numberOfBlocks;i++){

        faceLists[i]= PList_new();
        numInteriorFacesOfBlock[i]=0;

    }

    pPList* boundaryFaceLists = new pPList[numberOfBoundaryBlocks];
    for (i=0;i<numberOfBoundaryBlocks;i++){

        boundaryFaceLists[i]= PList_new();
        numBoundaryFacesOfBlock[i]=0;

    }




    // blockKey  defined in nspre_functions.h
    // uniquely identifies each block of elements 
    blockKey BLOCK;

    // iel, ielb assigned where ???
    std::map<int,int> iel ;
    std::map<int,int> ielb; 
 
    // declare region MeshSim datastructure
    pRegion region;

    // loop over regions (=elements)
    RIter rIter = M_regionIter(mesh);

    // index of regions
    i=0;

    // temporary face objects used in loops below
    pFace face;      
    pPList facesOfRegion;
    pPList facesOfBoundaryRegion;
    while (region = RIter_next(rIter) ){


        EN_getDataInt((pEntity) region, poly, &porder);
        /* processor for this region */

        
        nen =  numVertsR( region );
        
        
        ////////////////////////////////////////////////////////////////////////////////////
        // BLOCK definition (blockkey defined in nspre_data.h))
        // nen =  EN_dataI((pEntity)region,"RNEN"); 
        // used as in getConnectivity (above)
        // block indexing works how?
        ////////////////////////////////////////////////////////////////////////////////////
        BLOCK.nen = nen;
        BLOCK.maxpoly = porder;
        BLOCK.nenbl = nen == 8 ? 4 : 3;
        BLOCK.lcsyst = topology(region);
        blockid = Iblock[BLOCK]-1;

        nedges = nen+(nen+1)/2;
            
        // number of element faces for particular block
        // each block consists only of ONE type of topological elements 
        nfaces = nen == 8 ? 6: (nen-(nen/6));
        
        // parameter 1 is: faces in the order of MeshSim's R_face 
        facesOfRegion = R_faces(region,1);
      

        ////////////////////////////////////////////////////////////////////////////////
        // debug output control: dump vertices
        ////////////////////////////////////////////////////////////////////////////////
        void* itera=0;
        void* ite=0;
        pVertex vert;
        pPList vertList;
        
        // 1 is MeshSim cyclic order: 4th vertex in direction of normal
        // 0 is oppposite
//        vertList = R_vertices(region, 0); 
        pPList vertList1 = R_vertices(region, 1); 
        vertList = PList_new();
        int mapVerts[4] = {0,2,1,3};
        for (int iVert=0; iVert<4; iVert++)
            PList_append(vertList,PList_item(vertList1,mapVerts[iVert]));
        PList_delete(vertList1);
        

        int vrts[4];
            
        int l=0;
        while (vert =(pVertex)PList_next(vertList, &ite) ){
            
            vrts[l]=EN_id((pEntity)vert);
            l++;
        }          
           
        // assign faces correctly according to the re-orientation convention
        // of pasta (negative Volume) 
        int i0,i1,i2,i3;
        i0=0;
        i1=1;
        i2=2;
        i3=3;
        // flip nodes if 0th face normal points outside the region 
        // see R_entities(...)
        if ( R_faceDir( region, 0 )) {
            int tmpVertex = vrts[0];
            vrts[0] = vrts[2];
            vrts[2] = vrts[1];
            vrts[1] = tmpVertex;
            
            // re-ordering of local face numbers
            i1=3;
            i2=1;
            i3=2;

#ifdef DEBUG            
            cout<<"\ninverted!\n";
#endif // DEBUG            
        }            

#ifdef DEBUG
        cout<<"\nRegion "<<EN_id((pEntity)region)+1<<" has nodes\n";
        for(i = 0; i < 4; i++)
        {
            
            cout<<vrts[i]+1<<" ";
        }
#endif // DEBUG

        PList_delete(vertList);
        
              
        // loop over each region's faces
        
        i=0;// index for faces of region
        

        // face order has to be changed according to the renumbering of the
        // region's nodes
        Iter=0;
        
        k=0;
        while (face =(pFace)PList_next(facesOfRegion, &Iter) ){
            
            ////////////////////////////////////////////////////////////////////////////////
            // assign global face number to face array
            // first index is region Id of element in particular block
            // second index should be local face number for block ???
            ////////////////////////////////////////////////////////////////////////////////
            int ithindex=blockid;
            int jthindex=iel[blockid];
            int kthindex=NIblock[blockid];
            

            // adjust the face order according to the phasta convention
            if(i==0) k = i0;
            else if (i==1) k = i1;
            else if (i==2) k = i2;
            else if (i==3) k = i3;

            ief[blockid][iel[blockid]+k*NIblock[blockid]] = 
                EN_id((pEntity)face);
            
                

            PList_appUnique(faceLists[blockid],face);
                
                

            // increase total number of interior faces of this block
#ifdef DEBUG
            cout<<"\nface number (local "<<k+1<<"): "<<
                EN_id((pEntity)face)+1
                <<" for region "<<EN_id((pEntity)region)+1;
            cout<<" - face has nodes :\n";
            void* iterator=0;
            pVertex vertex;
            pPList vertexList; // list of nodes on face
            vertexList = F_vertices(face,1); // 1 is counterclockwise around
            // normal 
                
            while (vertex =(pVertex)PList_next(vertexList, &iterator) ){
                
                cout<<"  "<<EN_id((pEntity)vertex)+1<<"\n";
            } 
#endif // DEBUG            

            i=i+1;
        }// end loop over faces
        
        iel[blockid]++; /* increment the element counter for this block on 
                           this processor */

        PList_delete(facesOfRegion);


    
    }// end while-loop over regions

    RIter_delete(rIter);

    // extract each block's face number
    for(i=0;i<numberOfBlocks;i++){
        numInteriorFacesOfBlock[i]=PList_size(faceLists[i]);
    }

    


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // instantiate boundary face conectivity
    ////////////////////////////////////////////////////////////////////////////////////////////////
    void* temp = 0;
    pPList ents1;
    // loop over Mesh faces on model boundary 

    vector< pFace >::iterator bfaceIter = bdry.begin();
    while ( bfaceIter != bdry.end() ) {
        face = *bfaceIter;    


        // find the list of regions connected to this face. Since this
        // is a boundary face, there should only be one such region. 
        ents1 = F_regions(face);
        region = (pRegion)PList_item(ents1,0); /* get the region */
        PList_delete(ents1);
        
        EN_getDataInt((pEntity) region, poly, &porder);

        //nen = EN_dataI((pEntity)region,"RNEN");
        nen = numVertsR( region );
        
        BLOCK.nen = nen;
        BLOCK.maxpoly = porder;
        BLOCK.nenbl = F_numEdges(face);
        BLOCK.lcsyst = topology(region);
        if ( BLOCK.lcsyst == 3 ) BLOCK.lcsyst = BLOCK.nenbl;
        if ( BLOCK.lcsyst == 5 )
            if ( BLOCK.nenbl == 3 ) BLOCK.lcsyst = 6;
        
        blockid = Bblock[BLOCK]-1;
        nedges = nen+(nen+1)/2; 
        nfaces = nen == 8 ? 6: (nen-(nen/6));     

        // looping over the faces of the boundary element of the particular
        // block
        // Pick the first face: it is the face on the boundary by phasta
        // convention
        // parameter 1 is: faces in the order of MeshSim's R_face 
        facesOfBoundaryRegion = R_faces(region,1);
        // Assign the first face of connectivity as 
        // iefb(element#,1)=global id of face on the model boundary
        i=0;
            
        iefb[blockid][ielb[blockid]+NBblock[blockid]*i] =
            EN_id((pEntity)face);
        // cout<<"BOUNDARY face number: "<< EN_id((pEntity)face)<<" for region "<<EN_id((pEntity)region)<<"\n";

        // 
        PList_appUnique(boundaryFaceLists[blockid],face);


        

        // extract each block's face number
        for(i=0;i<numberOfBlocks;i++){
            numBoundaryFacesOfBlock[i]=PList_size(boundaryFaceLists[i]);
        }



        //  loop over each boundary region's faces
        // index for faces
        i=1;
        pFace otherFace; 
        Iter=0;
        while (otherFace = (pFace)PList_next(facesOfBoundaryRegion, &Iter) ){
            if(EN_id((pEntity)face)==EN_id((pEntity)otherFace)){
            } // already assigned above, so we do nothing here
            else{
                
                // assign global face number to face array
                // first index is region Id of element in particular block
                // second index should be local face number for block

                // assigning of the global BOUNDARY face number 
                iefb[blockid][ielb[blockid]+NBblock[blockid]*i] = 
                    EN_id((pEntity)otherFace);

                i=i+1;
            } // close of if loop to check for the first face
        }// end loop over faces
        PList_delete(facesOfBoundaryRegion); 
        ielb[blockid]++;
        bfaceIter++;
    } // end of the loop of the mesh faces defined on the model boundary       

}




/* ---------------------------------------------------------------
   return if the shape function associated with the given mode of 
   a finite element entity is ODD/EVEN
   Makes sense only for edge/face modes, for regions/vertex its
   always even.

   return 1 for EVEN mode
   0 for ODD  mode

   type = 1 ==> edge mode
   2 ==> tri-face mode

   pdof is the spectral order for the entity.
   ----------------------------------------------------------------- */
int 
modeSign(int type, int pdof) {
    if( type == 1 )         /* edge */
        return ( pdof%2 ? 0 : 1 );    
    else if( type == 2 )    /* tri face */
        return ( (pdof%3)%2 ? 0 : 1 ) ;
    else
        return 1;
}

int 
h_modeSign(int type, int pdof) {
    /* read comments above , pretty much the same thing but for hexes */
    if(type == 1) /* edge */
        return ((pdof%2) ? 0 : 1);
    else if(type == 2 )  /* quad face */
        return ((pdof%2) ? 0 : 1) ;
    else 
        return 1;
}

/*
  return the direction that a face is being used with respect to 
  how it was defined
*/
int 
Fsign(pFace face,pVertex v[3]) {

    pPList verts;
    pVertex vl[3];
    int j;

    /* get the vertices in the order they were defined */
    verts = F_vertices(face,1);
    for (j=0; j < 3; j++)
        vl[j] = (pVertex )PList_item(verts,j);
    PList_delete(verts);

    /* compare the vertices */
    if ( (v[0]==vl[0] && v[1]==vl[1])
         || (v[0]==vl[1] && v[1]==vl[2])
         || (v[0]==vl[2] && v[1]==vl[0]) ) return 1;
    else return 0; 
}
