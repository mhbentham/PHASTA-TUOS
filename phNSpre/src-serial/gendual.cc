#include <fstream>
#include <stdio.h>
#include "nspre_data.h"
#include "nspre_functions.h"
#include "phastaIO.h"

extern char* oformat;

void gendual( pGModel model, pMesh mesh ) {

    pRegion other;
    pRegion region;
    pPList verts,faces;
    pEntity ent;
    pGFace gface;
    void *temp=0;
    int numel=M_numRegions(mesh);
    int nIfaces=0;
    int tmp;  
    int iel,i;

    int  *xadj, *adjncy, *vwgt;
    int weight, sspebc;
    double coeffs[4];

    RIter rIter = M_regionIter(mesh);

    /* counting interior faces */

    FIter fIter = M_faceIter(mesh);
    pFace face;
      
    while (face = FIter_next(fIter)) {

        pPList regions = F_regions(face);
        if(PList_size(regions)==2) nIfaces++;
        PList_delete(regions);
    }
    FIter_reset(fIter);

    /* Allocating the necessary memory */

    xadj = new int [numel+1];
    adjncy = new int [2*nIfaces];
    vwgt = new int [numel];

    pVertex vertex;
    pGEntity gent;
    int Sregion=0;
    
    /* Finding out if SPEBC is used */
    GFIter gfiter = GM_faceIter(model);
    sspebc = 0;
    while(gface = GFIter_next(gfiter)){
#ifdef SIM
        if((GEN_dataI((pGEntity)gface,"sper"))) {
#else
        int tmp;
        if((GEN_dataI((pGEntity)gface,"sper",&tmp))) {
#endif
	        sspebc = 1;
	        break;
	    }
    }

        
    GFIter_delete(gfiter);

    if(sspebc) { 
	
        /* do this constrained partitioning only if the SPEBC is set  
         * here I set a tag "Locn" to each region, Locn = 0 for the regions 
         * on the master slave face pair , and Locn =1 for all other regions 
         */
	
        eqn_plane(model, mesh, coeffs, gface);
	    
    	int loc;	
        while(region = RIter_next(rIter)) {
            verts = R_vertices(region,1);
    	    int numVert = PList_size(verts);
    	    if (isCut(verts, coeffs)) loc = 0;
    	    else {
        		loc=1;
                for(i=0;i<numVert;i++){
                    tmp=1;
                    vertex=(pVertex)PList_item(verts,i);
                    gent=V_whatIn(vertex);
                    if(GF_inClosure(gface,gent)) tmp=0;
                    else tmp=1;
                    loc=loc*tmp;
		        }
            }
	    
            if(loc) {
                EN_attachDataI((pEntity)region,"Locn",1);
            } else {
                EN_attachDataI((pEntity)region,"Locn",0);
                Sregion++;
            }
            PList_delete(verts);
        }
        
        /* All the regions now have a location tag */

        /* SPEBC case requires weights to do constrained partition */
      
        int subnotset = 1;
        iel =0;            
        RIter_reset(rIter);
        while(region = RIter_next(rIter)){
            if(EN_dataI((pEntity)region,"Locn")) weight = 1;
            else if(subnotset) { subnotset = 0; weight=Sregion;}
            else weight=0;
            vwgt[iel] = weight;
            EN_attachDataI((pEntity)region, "WGHT", weight);
            iel++;
        }
    }    
    iel =0;
    int adj=0;
    RIter_reset(rIter);
    xadj[0] = 0;
    
    while(region=RIter_next(rIter)){
        if(!sspebc) vwgt[iel] = EN_dataI((pEntity)region,"WGHT");
        faces = R_faces(region,1);
        for (i =0; i<PList_size(faces);i++){
            face = (pFace)PList_item(faces,i);
            pPList fregions = F_regions(face);
            if (PList_size(fregions) == 2) {
                other = (pRegion)PList_item(fregions,0);
                if (other == region) other = (pRegion)PList_item(fregions,1);
                adjncy[adj++]=EN_id((pEntity)other);
            }
            PList_delete(fregions);
        }
        xadj[++iel]=adj;
        PList_delete(faces);
    }

//   ofstream dfile( "dual.dat" );
//
//   dfile << numel <<" "<< nIfaces <<" "<< sspebc << endl;
//   dfile.write( (char*)xadj , (numel+1)* sizeof( int ) );
//   dfile << endl;
//
//   dfile.write( (char*)adjncy , 2*nIfaces* sizeof( int ) );
//   dfile << endl;
//
//   dfile.write( (char*)vwgt , numel* sizeof( int ) );
//   dfile << endl;
//   dfile.close();
    
    char fname[255];
    int iarray[5];
    int fgeom, isize, nitems;

    sprintf(fname,"geombc.dat.1");
    openfile_( fname, "append", &fgeom );

    isize = numel+1;
    nitems = 2;
    iarray[0] = numel;
    iarray[1] = sspebc;
    writeheader_( &fgeom, "keyword xadj ", (void*)iarray, &nitems,
                  &isize, "integer", oformat );

    nitems = numel+1;
    writedatablock_( &fgeom, "keyword xadj ", (void*)xadj, &nitems,
                     "integer", oformat );

    isize = 2*nIfaces;
    nitems = 1;
    iarray[0] = nIfaces;
    writeheader_( &fgeom, "keyword adjncy ", (void*)iarray, &nitems,
                  &isize, "integer", oformat );

    nitems = 2*nIfaces;
    writedatablock_( &fgeom, "keyword adjncy ", (void*)adjncy, &nitems,
                     "integer", oformat );

    isize = numel;
    nitems = 1;
    iarray[0] = numel;
    writeheader_( &fgeom, "keyword vwgt ", (void*)iarray, &nitems,
                  &isize, "integer", oformat );

    nitems = numel;
    writedatablock_( &fgeom, "keyword vwgt ", (void*)vwgt, &nitems,
                     "integer", oformat );

    closefile_( &fgeom,"append");

    delete[] vwgt;
    delete[] xadj;
    delete[] adjncy;
    
    RIter_delete(rIter);
    FIter_delete(fIter);
}
