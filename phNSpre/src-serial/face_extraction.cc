/* Face co-ordinate and connectivity extraction on a given model face */
/* Anil Karanam  , Apr 2001  */

#include <cstdlib>
#include <iostream>		
#include <iomanip>
#include <fstream>
#include <vector>               
#include <map>

#ifdef SIM
#include "MeshSimInternal.h"
#include "MeshSim.h"
#else
#include "AOMD.h"
#include "MeshAdjTools.h"
#include "mesh_interface.h"
#endif
using namespace std;


extern int NSFTAG;
int* facetaglist;

/* The map class is an associative data structure which deals with
   pairs of objets <key , value>.  The pairs are internally kept
   sorted by the criterion we provide interms of the function object
   above

   In the next 3 instances of the map class we use mesh entities as
   keys and integers as their values, in this case we are interested
   in only two specific values of the 'value' either 0 (meaining that
   the given key is not a member of the list ) or 1 which means
   otherwise.
   
   Map has the important property that inserting a new element into a
   map does not invalidate iterators that point to existing
   elements. Erasing an element from a map also does not invalidate
   any iterators, except, of course, for iterators that actually point
   to the element that is being erased. This gives us tremendous
   flexibility in adding and removing elements to the vertex and edge
   lists we are currently working on.

   More information on the specifics of the map class can be found at.
   
   http://www.sgi.com/tech/stl/Map.html 

*/

/* List of counted faces  */
map<pFace , int > counted;    
/* List of vertices on the front bdry */
map<pEntity , int > vertex_list;
/* List of edges on the front bdry */
map<pEntity , int > edge_list; 

/* List of deleted entities */
map<pEntity, int > dvertex;

pGFace* gface;

/* vector is a simple extension of an array which does its own memory 
   management and can expand dynamically 

   http://www.sgi.com/tech/stl/Vector.html 
*/

/* This is the final destination for out newly counted faces */
/* we also maintain the counted facemap so that we don't have to completely */
/* traverse this everytime we need to check if a face has been counted.*/ 
/* like in the next function which is going to be called many times */

vector<pFace> facelist;

/* This function, give an mesh edge, model face */
/* retuns the other mesh face touching the edge and classified on the */
/* given model face, which hasn't been marked yet */

/*short for 'other_face_onthis_gface_and_not_marked_yet' */

pFace other_face(pEdge edge )
{
  pFace face;
  int numFaces = E_numFaces(edge);
 
  while ( numFaces-- > 0) {
    face = E_face(edge,numFaces); 
    if ( gface[0] == (pGFace)F_whatIn(face) )
      if ( !counted[face] ) return face;
  }

  /* If this function returns a NULL, 0,  then the current edge has no 
     business staying in the free boundary list */

  cerr<<"Problem in other_face: "<< EN_id((pEntity)edge) <<endl;
  exit(0);
  return 0;
}

void Add_Vertex(pEntity vert ){
  
  vertex_list[vert] = 1 ;
}
    
void Add_Edge(pEntity edge ){
  
  if ( E_whatInType((pEdge)edge) == Gedge )  return;

  /* sanity check */
  if ( E_whatIn((pEdge)edge) != (pGEntity)gface[0] ) {
    cerr<<"Error: New Edge not on model face \n";
    E_info((pEdge)edge);
    exit(0);
  }

  /* Also need to check if both atleast one face is not marked yet */

  pFace face, freeface=0;
  int numFaces = E_numFaces((pEdge)edge);
  
  while( numFaces-- > 0) {
    face = E_face((pEdge)edge,numFaces);
    if ( gface[0] == (pGFace)F_whatIn(face) )
      if ( !counted[face] )
	{
	  freeface=face;
	  break;
	}
  }

  if ( !freeface ){
    cerr<<"No Free faces available on:  "<< EN_id((pEntity)edge)<<endl;
    //    E_info((pEdge)edge);
    numFaces = E_numFaces((pEdge)edge);
    while ( numFaces-- > 0) F_info(E_face((pEdge)edge,numFaces));
    exit(0);
  }
  edge_list[edge] = 1 ;
}
    
  
/* main function to do face co-ordinate extraction */
void face_coordinate_extraction(pMesh mesh, pGModel model)
{
  ofstream coord("mesh.crd");
  ofstream fbct("faceBCT.nod");
  pVertex vertex;
  void* temp=0;
  double xyzC[3];
  int vid;
  gface = new pGFace [NSFTAG];
  facetaglist = new int [NSFTAG];

  for(int u=0; u< NSFTAG; u++)
    {
      cout<< "Enter Face tag: ";
      cin >> facetaglist[u];
    }

  for(int f=0; f<NSFTAG; f++) 
    gface[f] = (pGFace)GM_entityByTag(model,Gface,facetaglist[f]);

  /* First the coordinate extration */

  /* EDListIter<MVertex> gvIter = mesh->classifiedVertexIter(gface);
     while( vertex = gvIter.next() ) { */

  VIter vIter = M_vertexIter(mesh);
  while( vertex = VIter_next(vIter)) {
    /* incase you change your mind, oncomment the 2 lines above the while */
    /* comment out the existing while */
    /* also, comment out the line  "if ( V_whatIn(vertex) == gface)" */

    V_coord(vertex,xyzC);  
    xyzC[0]=xyzC[0]*0.1;
    xyzC[1]=xyzC[1]*0.1;
    xyzC[2]=xyzC[2]*0.1;   
    vid = EN_id((pEntity)vertex); 
    coord.setf(ios::scientific);
    coord << setprecision(17) ;
    coord <<vid<<"\t"<<xyzC[0]<<" "<<xyzC[1]<<" "<<xyzC[2]<< endl;
    int match =0;
    for(int h=0; h<NSFTAG; h++) 
      if( V_whatIn(vertex) == (pGEntity)gface[h] )  match = 1;

      if (match) fbct << vid << endl;
  }
  VIter_delete(vIter);
  coord.close();
  fbct.close();
  
  /* Then the connectivity extraction 
     we do this by the propagating front algorithm , where by we make,
     what we believe is a good start, by starting on a model edge and 
     keep traversing the free edge list while any free edges which are 
     not on model boundaries exist .

     In each pass of our traversal, we dynamically try to advance the edge 
     front, without voilating the condition of ``no loops on the boundary '' */

  /*  The first step is to make the good start we were talking about */ 
  /*  Get the first model edge of the modelface and the first mesh edge */

  /* first remove the intersection of the given list of face */

  if( NSFTAG > 1) {

      cout<<"face_coordinate_extraction for  NSFTAG > 1 is currently disabled.\n";
      cout<<"The algorithm relies on reclassifying mesh entities which is  a direct intervention\n";
      cout<<"into the mesh data base MDB\n";
      cout<<"from MeshSim version 6.0 on, no more direct manipulation of the MDB  is publicly available.\n";
       cout<<"Future implementations will have to handle this differently.\n";
      exit(1);

    pGEdge edge;
    pGFace face;
    void* elast = 0;
    pPList e_faces; 
  
    GEIter geIter = GM_edgeIter(model);
    while( edge = GEIter_next(geIter) ) {
      e_faces = GE_faces(edge);
      int matched = 0;
      void* etmp = 0;
      while ( face =(pGFace)PList_next(e_faces,&etmp))
	for(int f=0; f<NSFTAG;f++) if ( face == gface[f] ) matched++;

      PList_delete(e_faces); 

      if ( matched == 2 ) { /* seam --> unclassify */
        cout<<"Seam Edge : "<< GEN_tag((pGEntity)edge) << endl;
	EIter eIter = M_classifiedEdgeIter(mesh, (pGEntity)edge, 0);
	pEdge edgy ;
    while(edgy= EIter_next(eIter)){
        //E_setWhatIn(edgy, (pGEntity)gface[0]);
	}
	EIter_delete(eIter);
      }
    }
    GEIter_delete(geIter);

    for(int f=1; f<NSFTAG;f++) {
      /* reclassifying the edges */
      EIter eIter = M_classifiedEdgeIter(mesh,(pGEntity)gface[f],0);
      pEdge edgy;
      while( edgy = EIter_next(eIter)) 
          //E_setWhatIn(edgy,(pGEntity)gface[0]);
      EIter_delete(eIter);

      /* reclassifying  NSFTAG > 1the faces */
      FIter fIter = M_classifiedFaceIter(mesh,(pGEntity)gface[f],0);
      pFace facey;
      while ( facey = FIter_next(fIter)) 
          //F_setWhatIn(facey,(pGEntity)gface[0]);
      FIter_delete(fIter);
    }
  }
  pGEdge model_edge ;
  pPList gedgeList = GF_edges(gface[0]);
  void*  edge_restart=0 ;
  EIter  mesh_edge_iter;
  pEdge  meshEdge;

  do {

    model_edge = (pGEdge)PList_next(gedgeList,&edge_restart);
    mesh_edge_iter = M_classifiedEdgeIter(mesh,(pGEntity)model_edge,0);
    meshEdge = EIter_next(mesh_edge_iter);
    EIter_delete(mesh_edge_iter);

  } while ( !meshEdge );

  pFace firstMFace  = other_face(meshEdge);

  /* Now that we have the firstMeshFace, lets start building our list */

  counted[firstMFace]++;
  facelist.push_back(firstMFace);
  
  /* Add all the vertives */

  pPList vlist = F_vertices(firstMFace,1);
  void* fvtmp =0;
  while( vertex = (pVertex)PList_next(vlist,&fvtmp) ) 
    Add_Vertex((pEntity)vertex);
  PList_delete(vlist);

  /* Add all edges on the model face */

  pPList elist = F_edges(firstMFace,1,0);
  void* fetmp =0;
  pEdge edge;
  while( edge = (pEdge)PList_next(elist,&fetmp) ) 
    Add_Edge((pEntity)edge);
  PList_delete(elist);
  
//    {
//      EIter eIter = M_edgeIter(mesh);
//      pEdge edgy;
//      while( edgy = EIter_next(eIter)) {
//        if ( 12364 == EN_id((pEntity)edgy)){
//  	E_info(edgy);
//          int numFaces = E_numFaces(edgy);
//  	while ( numFaces-- > 0) F_info(E_face(edgy,numFaces));
//        }
//      }
//      EIter_delete(eIter);
//      exit(0);
//    }

  /* Now we start our edgelist traversal */
  
  map<pEntity, int >::iterator edgeiter;
  pFace candidate; /* candidate for inclusion  */
  pEntity nuke[3];
  int numchanges = 1;
  while ( numchanges ) { /* there still exists a growable  boundary */
    numchanges = 0;
    edgeiter = edge_list.begin();

    while ( edgeiter != edge_list.end() ) {
      /* we will try to pop each on of these edges and see if it voilates
         out no loop condition */
      if ( (*edgeiter).second ) {
	/* The current edge is not a deleted edge */

	candidate = other_face((pEdge)((*edgeiter).first));
	vertex = F_oppositeVertex( candidate ,(pEdge)((*edgeiter).first ));

	if ( !dvertex[(pEntity)vertex] ) {
	  if ( !vertex_list[(pEntity)vertex]){ /* not in list, no conflict */
	    /* we have a pop case where we add 2 edges and delete 1 */
	    numchanges++;
	    nuke[0] = (*edgeiter).first;
	    Add_Vertex((pEntity) vertex );
	    
	    pPList cedges = F_edges(candidate,1,0);
	    void* cetemp=0;
	    while( meshEdge = (pEdge)PList_next(cedges,&cetemp))
	      Add_Edge((pEntity) meshEdge);
	    PList_delete(cedges);
	    
	    edgeiter++;
	    edge_list[nuke[0]]= 0;
	    counted[candidate] = 1 ;
	    facelist.push_back(candidate);
	    
	  } else { /* the vertex is already in the list */
	    /* Two cases exist here, one acceptable and the other, not */
	    /* The acceptable case is the case of closing a valley */
	    /* for this we need to check if any of the other edges of the */
	    /* candidate face are in the current freeboundary list */
	    
	    pEdge medge_candidate;
	    
	    nuke[0] = NULL;
	    nuke[1] = NULL;
	    nuke[2] = NULL;
	    
	    int mctr =0 ;
	    
	    pPList elist1 = F_edges(candidate,1,0);
	    void* el1 = 0;
	    while( meshEdge = (pEdge)PList_next(elist1, &el1)) {
	      if ( edge_list[(pEntity)meshEdge] ) {
		nuke[mctr++] = (pEntity) meshEdge;
	      } else  {
		medge_candidate = meshEdge;
	      }
	    }
	    
	    if (nuke[1]) { /*  we have found a valley case */
	      
	      numchanges++;
	      nuke[2] = (pEntity)F_oppositeVertex(candidate,medge_candidate);
	      
	      counted[candidate] = 1;
	      facelist.push_back(candidate);
	      Add_Edge((pEntity)medge_candidate);
	      
	      edgeiter++;

	      edge_list[nuke[0]] = 0;
	      edge_list[nuke[1]] = 0 ;
	      dvertex[nuke[2]] = 1 ;
	      vertex_list[nuke[2]] = 0 ;
	      
	    } else { 
	      /* we are faced with an unacceptable case , loop formation*/
	      /* we don't do anything, but just move on */
	      
	      edgeiter++;
	    }
	  }
	}
      } else { edgeiter++ ;}
    }
  } 
  
  /* by the time the control reaches here, the mesh should have been traversed
     in the intended fashion and the final facelist should be ready to dump.
     following the rules of good housekeeping, free unneeded resources first */

  counted.clear();
  vertex_list.clear();
  edge_list.clear();
  dvertex.clear();

  ofstream bcnn("meshb.cnn");
  vector<pFace>::iterator cvIter;
  pRegion region;
  int dir;

  for( cvIter = facelist.begin(); cvIter != facelist.end(); cvIter++) {

    pPList regions = F_regions(*cvIter);
    void* frtmp=0;
    region = (pRegion)PList_next(regions, &frtmp);
    dir = R_dirUsingFace(region, *cvIter );
    PList_delete(regions);

    pPList vertices = F_vertices(*cvIter,dir);
    frtmp =0 ;
    
    while( vertex = (pVertex)PList_next(vertices,&frtmp) ) 
      bcnn << EN_id((pEntity)vertex) + 1 <<" ";

    bcnn << EN_id((pEntity)R_fcOpVt(region, *cvIter)) + 1 ;

    bcnn << endl;
  }

  cout <<"faces are in a valid reordering \n";

  facelist.clear();
  bcnn.close();
}
