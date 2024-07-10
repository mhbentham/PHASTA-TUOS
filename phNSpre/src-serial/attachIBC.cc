////////////////////////////////////////////////////////////////////////
//
// These functions set up the boundary condition information.
// 
////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <map>
#ifdef SIM
#include "MeshSim.h"
#else 
#include "AOMD.h"
#endif
#include "nspre_data.h"
#include "nspre_functions.h"
#include "EssentialBC.h"
#include "NaturalBC.h"
#include "InitialCondition.h"

#define RGAS 288.294382801664369

// global variables 
extern int rStart,nshgTot,intBC,prCd,ensa_dof,rRead;
extern globalInfo* info;
extern forwardblock Bblock;
extern pMeshDataId poly;

int 
attachEssentialBC( pGModel model, pMesh mesh, 
                   double *qTot, int *nBC,
                   int *iBC, double **BC ) {
    pGEntity gent;
    double   *intlist;
    int      code,numpbc=0,sign;
    int GDOFnum; // DOF number in qTot 
    ofstream fout;

    int done;

    // if there are internal boundary nodes (on a trip strip)
    // get their values now
    // not updated for 2D
    if (intBC == 1)  {
        GRIter grIter = GM_regionIter(model);
        gent    = (pGEntity)GRIter_next(grIter);
        GRIter_delete(grIter);
#ifdef SIM
        intlist = (double *)GEN_dataP(gent,"bc  ");
#else
        GEN_dataP(gent,"bc  ", (void **)&intlist);
#endif
        VIter vIter = M_vertexIter(mesh);
        pVertex vtx;
        while (vtx = VIter_next(vIter)) {
            code = EN_dataI((pEntity)vtx,"ibc ");
            if (code > 0) {
                nBC[Entity_getDOFnumber((pEntity)vtx,0)] = numpbc;
                iBC[numpbc] = code;
                for (int i=0; i < 11; i++) {
                    BC[numpbc][i] = intlist[i];
                }
                numpbc++;
            }
        }
        VIter_delete(vIter);
    }

    if (prCd)  {
        fout.open("BC.out", ios::app); 
        fout << "\n\nEssential boundary condition codes\n";
        fout << "----------------------------------\n\n";
        fout << "Entity" << "\tiBC-Code\n";
        fout << "------" << "\t--------\n";
    }
  
    // model faces
    GFIter gfIter = GM_faceIter(model);
    pGFace gface;
    while(gface = GFIter_next(gfIter))  {

        EssentialBC bc(gface);

#ifdef SIM
        double theta =  *(double*)(GEN_dataP((pGEntity)gface,"thta"));
#else
        double *ptheta;
        GEN_dataP((pGEntity)gface,"thta",(void **)&ptheta);
        double theta = *ptheta;
#endif
        if (prCd)  done = 0;
#ifdef SIM
        if (bc.isSet() || GEN_dataI((pGEntity)gface,"iper") 
            || GEN_dataI((pGEntity)gface,"sper") )  {
#else
            int tmp;
            if (bc.isSet() || GEN_dataI((pGEntity)gface,"iper", &tmp)
            || GEN_dataI((pGEntity)gface,"sper",&tmp))  {
#endif
            // mesh vertices on this face
            VIter vIter  = M_classifiedVertexIter(mesh, (pGEntity)gface, 0);
            pVertex vert;
            while (vert = VIter_next(vIter))  {
                if (EN_dataI((pEntity)vert,"ibc ") > 0) continue; 
                // nodal bc code attached
        
                if (bc.isSet())  {
                    if ((code = bc.eval(vert,BC[numpbc])) != 0)  {

                        // here is a function that will override the BC values
                        // obtained by evaluating the expressions (in the
                        // previous line) with those given in the initial
                        // condition (typically read in solution from a previous mesh)
                        if(bc.isAttSet(TBI)) bc.takeBCfromIC(BC[numpbc],
                                                             qTot,EN_id((pEntity)vert));

                        nBC[Entity_getDOFnumber((pEntity)vert,0)] = numpbc;
                        iBC[numpbc] = code;

                        sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
                        if (info->nsd == 3) BC[numpbc][11] = sign*theta;
  
                        if (prCd && !done)  {
                            fout << "F-" << GEN_tag((pGEntity)gface)<< "\t\t" ;
                            fout << code << "\n";
                            done = 1;
                        }
                        numpbc++;
                    }
                } else {  
                    code =0;
#ifdef SIM
                    if (GEN_dataI((pGEntity)gface,"iper")) {
#else
                    int tmp;
                    if (GEN_dataI((pGEntity)gface,"iper",&tmp)) {
#endif
                        code = 1024;
                    }
                    sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
                    if (info->nsd == 3) BC[numpbc][11] = sign*theta;
#ifdef SIM
                    if (GEN_dataI((pGEntity)gface,"sper"))  code += 2048;
#else
                    if (GEN_dataI((pGEntity)gface,"sper",&tmp))  code += 2048;
#endif
                    nBC[Entity_getDOFnumber((pEntity)vert,0)] = numpbc;
                    iBC[numpbc] = code;
  
                    if (prCd && !done)  {
                        fout << "F-" << GEN_tag((pGEntity)gface)  << "\t\t" ;
                        fout << code << "\n";
                        done = 1;
                    }
                    numpbc++;
                }
            }
            VIter_delete(vIter);

            // mesh edges

            EIter eIter = M_classifiedEdgeIter(mesh,(pGEntity)gface, 0);
            pEdge edge;
            int nem;
            while (edge = EIter_next(eIter))  {
                if ( isActive( (pEntity) edge ) ) {
                    nem = Entity_getNumDOF((pEntity)edge);
                    for (int i=0; i < nem; i++)  {
                        if (bc.isSet())  {
                            if ((code = bc.evalEdge(edge,BC[numpbc],i+2)) != 0)  {
                                nBC[Entity_getDOFnumber((pEntity)edge,i)] = numpbc;
                                iBC[numpbc] = code;
                                sign=((EN_dataI((pEntity)edge,"swit"))? -1 : 1);
                                if (info->nsd == 3) BC[numpbc][11] = sign*theta;
                                numpbc++;
                            }
                        }
                        else  {
                            nBC[Entity_getDOFnumber((pEntity)edge,i)] = numpbc;
                            code=0;
                            sign=((EN_dataI((pEntity)edge,"swit"))? -1 : 1);
                            if (info->nsd == 3) BC[numpbc][11] = sign*theta;
#ifdef SIM
                            if (GEN_dataI((pGEntity)gface,"iper")) {
                                code = 1024;
                            }
                            if (GEN_dataI((pGEntity)gface,"sper"))  code += 2048;
#else
                            if (GEN_dataI((pGEntity)gface,"iper",&tmp)) {
                                code = 1024;
                            }
                            if (GEN_dataI((pGEntity)gface,"sper",&tmp))  code += 2048;

#endif

                            iBC[numpbc] = code;
                            numpbc++;
                        }
                    }
                }
            }
            EIter_delete(eIter);

            // mesh faces
            FIter fIter  = M_classifiedFaceIter(mesh,(pGEntity)gface, 0);
            pFace face;
            int nfm;
            while (face = FIter_next(fIter))  {
                if (isActive((pEntity)face) ) {
                    nfm = Entity_getNumDOF((pEntity)face);
                    for (int i=0; i < nfm; i++)  {
                        if (bc.isSet())  {
                            if ((code = bc.evalFace(face,BC[numpbc],i+3)) != 0)  {
                                nBC[Entity_getDOFnumber((pEntity)face,i)] = numpbc;
                                iBC[numpbc] = code;
                                numpbc++;
                            }
                        }
                        else  {
                            nBC[Entity_getDOFnumber((pEntity)face,i)] = numpbc;
                            code=0;
                            sign=((EN_dataI((pEntity)face,"swit"))? -1 : 1);
                            if (info->nsd == 3) BC[numpbc][11] = sign*theta;
#ifdef SIM
                            if (GEN_dataI((pGEntity)gface,"iper")) {
#else
                            if (GEN_dataI((pGEntity)gface,"iper",&tmp)) { 
#endif
                                code = 1024;
                            }
                            iBC[numpbc] = code;
		
                            numpbc++;
                        }
                    }
                }
            }
        }
    }
    GFIter_delete(gfIter);
  
    // model edges

    GEIter geIter = GM_edgeIter(model);
    pGEdge gedge;
    while( gedge = GEIter_next(geIter) )  {
 
        // axisymmetric bc's
#ifdef SIM
        double theta = *(double*)(GEN_dataP((pGEntity)gedge,"thta"));
#else
        double *ptheta;
        GEN_dataP((pGEntity)gedge,"thta",(void **)&ptheta);
        double theta = *ptheta;
#endif
  
        EssentialBC bc(gedge);
    
        if (prCd)  done = 0;
        if (bc.isSet())  {
            // mesh vertices

            VIter vIter  = M_classifiedVertexIter(mesh,(pGEntity)gedge, 0);
            pVertex vert;

            while (vert = VIter_next(vIter))  {
                if (EN_dataI((pEntity)vert,"ibc ") > 0) continue;
                if (isActive((pEntity)vert)) {
                    if ((code = bc.eval(vert,BC[numpbc])) != 0)  {
                        // here is a function that will override the BC values
                        // obtained by evaluating the expressions (in the
                        // previous line) with those given in the initial
                        // condition (typically read in solution from a previous mesh)
                        if(bc.isAttSet(TBI)) bc.takeBCfromIC(BC[numpbc],
                                                             qTot,EN_id((pEntity)vert));

                        nBC[Entity_getDOFnumber((pEntity)vert,0)] = numpbc;
                        iBC[numpbc] = code;
                        sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
                        if (info->nsd == 3) BC[numpbc][11] = sign*theta;
                        if (prCd && !done)  {
                            fout << "E-" << GEN_tag((pGEntity)gedge) << "\t\t" << code << "\n";
                            done = 1;
                        }
                        numpbc++;
                    }
                }
            }
            VIter_delete(vIter);
            // mesh edges
            EIter eIter  = M_classifiedEdgeIter(mesh, (pGEntity) gedge, 0 );
            pEdge edge;
            int nem;
            while (edge = EIter_next(eIter))  {
                if (isActive((pEntity)edge))  {
                    nem = Entity_getNumDOF((pEntity)edge) ;
                    for (int i=0; i < nem; i++)  {
                        if ((code = bc.evalEdge(edge,BC[numpbc],i+2)) != 0)  {
                            nBC[Entity_getDOFnumber((pEntity)edge,i)] = numpbc;
                            iBC[numpbc] = code;
                            sign=((EN_dataI((pEntity)edge,"swit"))? -1 : 1);
                            if (info->nsd == 3) BC[numpbc][11] = sign*theta;
                            numpbc++;
                        }
                    }
                }
            }
            EIter_delete(eIter);
        }
    }
    GEIter_delete(geIter);

    // model vertices
    GVIter gvIter = GM_vertexIter(model);
    pGVertex gvert;
    while(gvert = GVIter_next(gvIter)) {
        EssentialBC bc(gvert);
        
        // axisymmetric bc's
#ifdef SIM
        double theta =  *(double*)(GEN_dataP((pGEntity)gvert,"thta"));
#else
        double *ptheta;
        GEN_dataP((pGEntity)gvert,"thta",(void **)&ptheta);
        double theta = *ptheta;
#endif
     
        if (prCd)  done = 0;
        if (bc.isSet())  {

            // mesh vertex
            pVertex vert;
            vert = VIter_next(M_classifiedVertexIter(mesh, (pGEntity)gvert, 0));
             if (EN_dataI((pEntity)vert,"ibc ") > 0 ) continue;
      
            if (isActive((pEntity)vert)) {
                if ((code = bc.eval(vert,BC[numpbc])) != 0)  {
                    // here is a function that will override the BC values
                    // obtained by evaluating the expressions (in the
                    // previous line) with those given in the initial
                    // condition (typically read in solution from a previous mesh)
                    if(bc.isAttSet(TBI)) bc.takeBCfromIC(BC[numpbc],
                                                         qTot,EN_id((pEntity)vert));

                    nBC[Entity_getDOFnumber((pEntity)vert,0)] = numpbc;
                    iBC[numpbc] = code;
                    sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
                    if (info->nsd == 3) BC[numpbc][11] = sign*theta;
	  
                    if (prCd && !done)  {
                        fout << "V-" << GEN_tag((pGEntity)gvert) << "\t\t" << code << "\n";
                        done = 1;
                    }
                    numpbc++;
                }
            }
        }
    }
    GVIter_delete(gvIter);
    if (prCd ) fout.close();
    /* before we quit we want to add the constraining bcs for the fake node  */
    if ( info->fake_mode ) { 
        code =0; 
        for( int u=0; u< 10; u++ ) code = setbit( code, u );
        nBC[ info->nshg - 1 ] = numpbc;
        iBC[ numpbc ] = code ;
        int numEBC=ensa_dof+7;
        for( int u=0; u < numEBC; u++ ) BC[numpbc][u] =0.0;
        BC[numpbc][3] = BC[numpbc][7] = 1.0;
        numpbc++;
    }
    return numpbc;
}

/**********************************************************************/
/* Natural boundary conditions                                        */
/**********************************************************************/
void 
attachNaturalBC( pGModel model, pMesh mesh, int ***iBCB, double ***BCB, int numNBC) {
    blockKey BLOCK;
    int blockid;
    std::map<int, int> iel;
    pPList ents;
    pRegion region;
    pFace face;

    int code;
    int done;
    ofstream fout("BC.out", ios::app); 

    if (prCd) {
        fout << "\n\nNatural boundary condition codes\n";
        fout <<     "--------------------------------\n\n";
        fout << "Entity" << "\tiBCB-Code\n";
        fout << "------" << "\t--------\n";
    }

    // For 3D loop over model faces and evaluate the natural BC over meshes
    // faces, then block the corresponding boundary element
    if (info->nsd == 3){
       GFIter gfIter = GM_faceIter(model);
       pGFace gface;
       while(gface = GFIter_next(gfIter))  {
    
          NaturalBC bc(gface);

          if (prCd)  done = 0;
          if (bc.isSet())  {                // if code is nonzero 
             FIter fIter  = M_classifiedFaceIter(mesh,(pGEntity)gface,0);
             while (face = FIter_next(fIter))  {
                ents = F_regions(face);
                region = (pRegion)PList_item(ents,0);
                PList_delete(ents);
                BLOCK.nen = numVertsR( region );
                EN_getDataInt( (pEntity) region, poly, &BLOCK.maxpoly );
                BLOCK.nenbl = F_numEdges(face);
                BLOCK.lcsyst = topology(region);
                if( 3 == BLOCK.lcsyst ) BLOCK.lcsyst = BLOCK.nenbl;
                if( 5 == BLOCK.lcsyst )
                    if ( 3 == BLOCK.nenbl ) BLOCK.lcsyst = 6;
                blockid = Bblock[BLOCK]-1;
                bc.evalFace(face,BCB[blockid][iel[blockid]],iBCB[blockid][iel[blockid]], numNBC);
                code=iBCB[blockid][iel[blockid]][1]; 
                if (prCd && !done)  {
                    fout << "F-" << GEN_tag((pGEntity)gface) << "\t\t" << code << "\n";
                    done = 1;
                }
                iel[blockid]++;
             }
             FIter_delete(fIter);
          }
       }
       GFIter_delete(gfIter);

    // For 2D loop over model edges and evaluate the natural BC over meshes
    // edges, then block the corresponding boundary element
    } else {
       GEIter geIter = GM_edgeIter(model);
       pGEdge gedge;
       while(gedge = GEIter_next(geIter))  {
    
          NaturalBC bc(gedge);

          if (prCd)  done = 0;
          if (bc.isSet())  {                // if code is nonzero 
             EIter eIter  = M_classifiedEdgeIter(mesh,(pGEntity)gedge,0);
             pEdge edge;
             while (edge = EIter_next(eIter))  {
                ents = E_faces(edge);
                face = (pFace)PList_item(ents,0);
                PList_delete(ents);
                BLOCK.nen = numVertsF( face );
                EN_getDataInt( (pEntity) face, poly, &BLOCK.maxpoly );
                BLOCK.nenbl = 2;
                BLOCK.lcsyst = topology2D(face);
                blockid = Bblock[BLOCK]-1;
                bc.evalEdge(edge,BCB[blockid][iel[blockid]],iBCB[blockid][iel[blockid]], numNBC);
                code=iBCB[blockid][iel[blockid]][1]; 
                if (prCd && !done)  {
                    fout << "F-" << GEN_tag((pGEntity)gedge) << "\t\t" << code << "\n";
                    done = 1;
                }
                iel[blockid]++;
             }
             EIter_delete(eIter);
          }
       }
       GEIter_delete(geIter);
    }
    iel.clear();
}

    
/**********************************************************************/
/* Initial condition                                                  */
/**********************************************************************/
void 
attachInitialCondition( pGModel model, pMesh mesh,
                        double *qTot, double **q ) {
    pVertex  vertex;
    pEdge    edge;
    pFace    face;
    pRegion  region;
    int      i,j,numnp=M_numVertices(mesh),count,nv,nem=0,nfm=0,nrm=0;
    int nshg = info->nshg;
    nv=ensa_dof;
  
    // user is supplying a restart file to use as
    // the initial condition
    if (rStart)  {
        // vertex modes
        VIter vIter = M_vertexIter(mesh);
        while (vertex = VIter_next(vIter))  {
            for (i=0; i < nv; i++)
                if (rStart == 2 && i==0)
                    /* convert from density to pressure */
                    q[Entity_getDOFnumber((pEntity)vertex,0)][0] =
                        qTot[0*nshg+EN_id((pEntity)vertex)]
                        *RGAS
                        *qTot[4*nshg+EN_id((pEntity)vertex)];
                else
                    q[Entity_getDOFnumber((pEntity)vertex,0)][i] =
                        qTot[i*nshg+EN_id((pEntity)vertex)];
        }
        VIter_delete(vIter);
        // restart file contains higher order modes
        if (rStart == 4) { // only quadratic in restart file
            cout << "RESTART OPTION WITH C-WRITES NOT READY FOR rStart==4"<<endl;
            count = M_numVertices(mesh);
            EIter eIter = M_edgeIter(mesh);
            while (edge = EIter_next(eIter))  {
                for (i=0; i < nv; i++)  {
                        q[Entity_getDOFnumber((pEntity)edge,0)][i] =
                            qTot[count+i*nshg];
                    }
                    count++;
            }
            EIter_delete(eIter);
        } else if (rStart != 3) { 
            count = M_numVertices(mesh);
            // edge modes 
            EIter eIter = M_edgeIter(mesh);
            while (edge = EIter_next(eIter))  {
                nem = Entity_getNumDOF((pEntity)edge);
                    for (j=0; j < nem; j++)  {
                        for (i=0; i < nv; i++)  {
                            q[Entity_getDOFnumber((pEntity)edge,j)][i] = qTot[i*nshg+count];
                        }
                        count++;
                    }
            }
            EIter_delete(eIter);
            // face modes 
            FIter fIter = M_faceIter(mesh);
            while (face = FIter_next(fIter))  {
                nfm = Entity_getNumDOF((pEntity)face);
                    for (j=0; j < nfm; j++)  {
                        for (i=0; i < nv; i++)  {
                            q[Entity_getDOFnumber((pEntity)face,j)][i] = qTot[count+i*nshg];
                        }
                        count++;
                    }
            }
            FIter_delete(fIter);
        }

        // intitial condition is generated from the attribute
        // manager
    } 
    
    // find the initial condition attribute expression
    InitialCondition ic(model);
    
    VIter vIter = M_vertexIter(mesh);
    while (vertex = VIter_next(vIter))  {
        if(isActive((pEntity)vertex))  {
            ic.eval(vertex,q[Entity_getDOFnumber((pEntity)vertex,0)]);
        }
    }
    VIter_delete(vIter);
  
    // edge modes
    EIter eIter = M_edgeIter(mesh);
    while (edge = EIter_next(eIter))  {
        if(isActive((pEntity)edge))  {
            nem = Entity_getNumDOF((pEntity)edge);
            if (nem > 0)  {
                for (j=0; j < nem; j++)  {
                    for (i=0; i < nv; i++)  {
                        ic.eval(edge,q[Entity_getDOFnumber((pEntity)edge,j)],j+2);
                    }
                }
            }
        }
    }
    EIter_delete(eIter);
    
    // face modes 
    FIter fIter = M_faceIter(mesh);
    while (face = FIter_next(fIter))  {
        if(isActive((pEntity)face))  {
            nfm = Entity_getNumDOF((pEntity)face);
            for(j=0;j<nfm;j++) {
                for (i=0; i < nv; i++)  {
                    q[Entity_getDOFnumber((pEntity)face,j)][i]=0.0;
                }
            }
        } 
    }
    FIter_delete(fIter);
    if ( info->fake_mode ) {
        for (i=0; i < nv; i++)  {
            q[ info->nshg - 1 ][i]=0.0;
        }
    }
}
