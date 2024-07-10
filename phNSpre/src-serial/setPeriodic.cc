////////////////////////////////////////////////////////////////////////
//
// This function attaches the periodic mesh partner to each mesh      
// entity classified on a model entity which has the attribute
// "periodic slave"
//
// Each of these periodic slave model entities will find periodic
// masters for all of the mesh entities classified on it.
//
// A special case arises when a model entity is multiply periodic
// (periodic in more than one spatial dimension). In this case, an
// integer representing the number of jumps to the master (equivilently,
// the number of periodic dimensions) is also assigned as an attribute,
// designated "number of jumps". This information is necessary to
// correctly find a mesh entity's periodic partner.
// 
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <math.h>
#include "nspre_data.h"
#include "PeriodicBC.h"
#include "NaturalBC.h"
#include "SPEBC.h"
#include "nspre_functions.h"

#define ATOL 0.01

extern int prCd;
extern int globalP;
extern globalInfo* info;
ofstream perr;

void 
setPeriodic ( pGModel model, pMesh mesh ) {
    pVertex vertex, mvert;
    pEdge medge;
    pFace mface;
    pGEntity gment;
    pGVertex gvert;
    pGEdge gedge;
    pGFace gface;
    void *temp;
    int njumps;
    perr.open("periodic_errors.log");
    ofstream fout;

    if (prCd)  {
        fout.open("BC.out");
        fout << "Slave entity" << "\tMaster entity" << "\tNum jumps\n";
        fout << "------------" << "\t-------------" << "\t---------\n";
    }

    temp = 0;
    GFIter gftr1 = GM_faceIter(model);
    while (gface = GFIter_next(gftr1))  { // over model faces
        pPList gf_regions = GF_regions(gface);
        if (PList_size(gf_regions) == 2) { // interior face

            GEN_attachDataI((pGEntity)gface,"intr",1);
            GEN_attachDataI((pGEntity)gface,"NoBI",1);

            pPList gf_edges = GF_edges(gface);
            void* etmp = 0;
            while( pGEdge gedge = (pGEdge)PList_next(gf_edges,&etmp)) 
                GEN_attachDataI((pGEntity)gedge,"intr",1);
            PList_delete(gf_edges);

            pPList gf_vertices = GF_vertices(gface);
            void* vtmp = 0;
            while( pGVertex gvertex = (pGVertex)PList_next(gf_vertices,&vtmp)) 
                GEN_attachDataI((pGEntity)gvertex,"intr",1);
            PList_delete(gf_vertices);
      
        }
        // 2D mesh: the face is the model and does not have any regions
	// attached to it
	if (PList_size(gf_regions) == 0) { 

            GEN_attachDataI((pGEntity)gface,"NoBI",1);
        }
    }
  
    // loop over all model entities (faces, edges, and vertices) and
    // if they are a "real" periodic slave( as in not an SPEBC), find 
    // the master of all the  mesh entities classified on them.
    temp=0;
    GFIter_reset(gftr1);
    while (gface = GFIter_next(gftr1)) {// over model faces
        pGEntity gent = (pGEntity) gface;    
#ifdef SIM  
        if (!GEN_dataI(gent,"intr"))  {
#else
         int tmp;
        if (!GEN_dataI(gent,"intr",&tmp))  {
#endif

            SPEBC spbc(gface);	
      
            if (spbc.isSPEBC()) {
                GEN_attachDataI(gent, "sper", 1);
            } 
      
            PeriodicBC bc(model, gface);
            njumps = bc.getPerMaster ((pGFace*) &gment);
      
            if (njumps != 0)  {                

                // check for axisymmetric boundary condition
                // If the slave face makes a nonzero angle with the master
                // then we assume that we have an Axisymmetric BC case
  
                // Here theta represents the  rotation for the slave face to
                // the master face and hence -theta from the masterface to the
                // slave face.
        
                double theta = bc.getAngle((pGFace) gment, mesh);
                if (fabs(theta) < ATOL ) theta = 0.0;
                double *ttmp = new double;
                double *ttmpn = new double;
                *ttmp = theta;
                *ttmpn = -theta;
                double *doubletmp;
                // If there's no "thta" yet for slave, attach it and set to
                //  wherever ttmp points
                // If "thta" has been initialized already for slave, change
                //  it to point where ttmp points
#ifdef SIM
                if (!GEN_dataP(gent,"thta")){
#else
                if (!GEN_dataP(gent,"thta",(void**)&doubletmp)){
#endif
                    GEN_attachDataP(gent,"thta",ttmp);
                } else {
                    GEN_modifyDataP(gent,"thta",ttmp);
                }
                // Similar procedure for master face
#ifdef SIM
                if (!GEN_dataP(gment,"thta")){
#else
                if (!GEN_dataP(gment,"thta",(void**)&doubletmp)){
#endif
                    GEN_attachDataP(gment,"thta",ttmpn);
                } else {
                    GEN_modifyDataP(gment,"thta",ttmpn);
                }

                if (prCd)  fout << "F-" 
                                << GEN_tag((GEntity *) gface) 
                                << "\t\tF-"
                                << GEN_tag(gment) 
                                << "\t\t" << njumps << "\n";
	
       
                NaturalBC nbc(gface);
                if(!(nbc.isAttSet(SID))){       // for faces with Surf ID
                    GEN_attachDataI (gent, "NoBI", 1);    // flag as noBI
                } 
	
                // wont count in bdry integral
                GEN_attachDataI (gment, "NoBI", 1);   

                GEN_attachDataI(gent, "iper", 1);
	
                VIter vIter = M_classifiedVertexIter(mesh,gent,0);
                VIter mvIter = M_classifiedVertexIter(mesh,gment,0);

                pVertex vtx;
                while (vtx = VIter_next(vIter))  {                 // mesh verts
                    mvert = getMasterV (mvIter, vtx, theta, 1);
                    EN_attachDataP ((pEntity)vtx, "PerM", mvert);
                    VIter_reset( mvIter );
                }

                VIter_delete(vIter);
                VIter_delete(mvIter);

                if ( globalP > 1 ) {
                    // mesh edges
                    EIter eIter  = M_classifiedEdgeIter(mesh,gent,0);
                    EIter meIter = M_classifiedEdgeIter(mesh,gment,0);
                    pEdge edg;
                    while (edg = EIter_next(eIter))  { 
                        medge = getMasterE (meIter, edg, theta, 1);
                        EN_attachDataP ((pEntity)edg, "PerM", medge);
                        EIter_reset( meIter );
                    }
                    EIter_delete(eIter);
                    EIter_delete(meIter);
        
                    if ( globalP > 2 ) {
                        // mesh faces
                        FIter fIter  = M_classifiedFaceIter(mesh,gent,0);
                        FIter mfIter = M_classifiedFaceIter(mesh,gment,0);
                        pFace fac;
                        while (fac = FIter_next(fIter))  {
                            mface = getMasterF (mfIter, fac, theta);
                            EN_attachDataP ((pEntity)fac, "PerM", mface);
                            FIter_reset( mfIter );
                        }
                        FIter_delete(fIter);
                        FIter_delete(mfIter);
                    }
                }
            }
        }
#ifdef SIM    
        if (!GEN_dataP(gent,"thta")){
#else
        double doubletmp;
        if (!GEN_dataP(gent,"thta",(void**)&doubletmp)){
#endif
            double* thetat = new double;
            *thetat = 0.0;
            GEN_attachDataP(gent,"thta",thetat);
        }
    }
    GFIter_delete(gftr1);
    temp = 0;
    GEIter getr1 = GM_edgeIter(model);
    while (gedge = GEIter_next(getr1)){          // all model edges
    
        SPEBC spbc(gedge);	
        pGEntity gent = (pGEntity) gedge;    
        if (spbc.isSPEBC()) {
            GEN_attachDataI(gent, "sper", 1);
        } 
        PeriodicBC bc(model, gedge );
        njumps = bc.getPerMaster((pGEdge *) &gment);
    
        // Recognizing the axisymmetric centerline
        if ( njumps == 2475 ) { // magic number
      
            GEN_attachDataI(gent,"AXCL", 1);
      
        } else if (njumps != 0)  {
      
            double theta = bc.getAngle((pGEdge) gment);
            if( fabs(theta) < ATOL ) theta = 0.0;
            double *ttmp = new double;
            *ttmp = theta;
            GEN_attachDataP(gent,"thta",ttmp);
      
      
            if (prCd)  fout << "E-" << GEN_tag((GEntity *) gedge) << "\t\tE-"
                            << GEN_tag(gment) << "\t\t" << njumps << "\n";
      
            if (info->nsd  == 2) {
	        // periodic master and slave do not count in bdry integral
                GEN_attachDataI (gment, "NoBI", 1);   
                GEN_attachDataI (gent, "NoBI", 1);   
	    }

            GEN_attachDataI(gent, "iper", 1);
      
            VIter vIter = M_classifiedVertexIter(mesh,gent,0);
            VIter mvIter = M_classifiedVertexIter(mesh,gment,0);
            pVertex vtx;
      
            if (info->nsd == 3) {
	      while (vtx = VIter_next(vIter)){                  // mesh vertices
                mvert = getMasterV (mvIter, vtx, theta, njumps);
                VIter_reset( mvIter);
                if (mvert) {
                    EN_attachDataP ((pEntity)vtx, "PerM", mvert);
                } else {
                    cerr << " special case: periodic circular faces where model";
                    cerr << " verts don't line up" << endl;
                    cerr << " Warning: taking bounding vert as periodic master";
                    cerr << endl;
                    VIter mvertiter  = 
                        M_classifiedVertexIter( mesh,
                                                (pGEntity)GE_vertex((pGEdge)gment,0),0);
                    mvert = VIter_next(mvertiter);
                    VIter_delete(mvertiter);
                }
	      }
            } else {  //2D mesh
                while (vtx = VIter_next(vIter))  {                 // mesh verts
                    mvert = getMasterV (mvIter, vtx, theta, 1);
                    EN_attachDataP ((pEntity)vtx, "PerM", mvert);
                    VIter_reset( mvIter );
                }
	    }
            VIter_delete(vIter);
            VIter_delete(mvIter);
      
            if ( globalP > 1 ) {
                // mesh edges
                EIter eIter = M_classifiedEdgeIter(mesh,gent,0);
                EIter meIter = M_classifiedEdgeIter(mesh,gment,0);
                pEdge edg;
                while (edg = EIter_next(eIter)){ 
                    medge = getMasterE (meIter, edg, theta, njumps);
                    EN_attachDataP ((pEntity)edg, "PerM", medge);
                    EIter_reset( meIter );
                }
                EIter_delete(eIter);
                EIter_delete(meIter);
            }
        }

#ifdef SIM    
        if (!(GEN_dataP(gent,"thta"))){
#else
        double doubletmp;
        if (!GEN_dataP(gent,"thta",(void**)&doubletmp)){
#endif
            double* thetat = new double;
            *thetat = 0.0;
            GEN_attachDataP(gent,"thta",thetat);
        }
    }
    GEIter_delete(getr1);
  
    temp = 0;
    GVIter gvtr1 = GM_vertexIter(model);
    while (gvert = GVIter_next(gvtr1)){       // all model vertices
    
        SPEBC spbc(gvert);	
    
        pGEntity gent = (pGEntity) gvert;    
    
        if (spbc.isSPEBC()) {
            GEN_attachDataI(gent, "sper", 1);
        } 
    
        PeriodicBC bc(model, gvert );
        njumps = bc.getPerMaster((pGVertex *) &gment);
    
        // Recognizing the axisymmetric centerline
    
        if ( njumps == 2475 ) { // magic number
      
            GEN_attachDataI(gent,"AXCL", 1);
      
        } else if (njumps != 0)  {              // means the entity is a slave
      
            double theta = bc.getAngle((pGVertex) gment);
            if (fabs(theta) < ATOL) theta = 0.0;
            double *ttmp = new double;
            *ttmp = theta;
            GEN_attachDataP(gent,"thta",ttmp);
      
            if (prCd)  fout << "V-" << GEN_tag((GEntity *) gvert) << "\t\tV-"
                            << GEN_tag(gment) << "\t\t" << njumps << "\n";
      
            GEN_attachDataI(gent, "iper", 1);
      
            VIter vIter = M_classifiedVertexIter(mesh,gent,0);
            VIter mvIter = M_classifiedVertexIter(mesh,gment,0);

            vertex = VIter_next(vIter);
            mvert  = VIter_next(mvIter);

            VIter_delete(vIter);
            VIter_delete(mvIter);
      
            // attach the periodic master (pointer) to the slave mesh entity
            EN_attachDataP ((pEntity)vertex, "PerM", mvert);
        }
#ifdef SIM
        if (!(GEN_dataP(gent,"thta"))){
#else
        double doubletmp;
        if (!GEN_dataP(gent,"thta",(void**)&doubletmp)){
#endif
            double* thetat = new double;
            *thetat = 0.0;
            GEN_attachDataP(gent,"thta",thetat);
        }
    }
    GVIter_delete(gvtr1);
    perr.close();
}
