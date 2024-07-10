#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cstring>
#include <map>
#include <utility>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "phastaIO.h"
#include <FCMangle.h>

#define sonfath FortranCInterface_GLOBAL_(sonfath,SONFATH)
#define for if (0) {} else for
#ifndef WIN32
#include <unistd.h>
#include <strings.h>
#else
#include <direct.h>
void  bzero(void* ptr, size_t sz) {
    int i;
    char *cptr;
    cptr = (char*) ptr;
    for (i=0; i < sz; i++) {
        cptr[i]=0;
    }
    return;
}
#endif

char keyphrase[100];
char filename[255];
char dirname[25];
char dirname2[50];

extern "C" {
    void
    METIS_PartGraphKway(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
//    METIS_PartGraphKway(int*,long int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
    void sonfath( int*, int*, int*, int*, int*, int*, int*, int* );
}

void 
Gather_Headers( int* fileDescriptor, std::vector< std::string >& headers );
    
using namespace std;

typedef vector< int > block ;
map< int, map< int, int >  > ParallelData;
 
// ParallelData is the most important data structure of this code it maitains 
// the  global to local mapping of shapefuntion numbers 
// ParallelData[ global_number][ processor_id ] = processor_local_no 
// all shape function numbers, global and local are stored in fortran indexing

struct lessKey{
    bool operator()(block b1, block b2) const
        {
            return ( 100*b1[1] + 10*b1[5] < 100*b2[1] + 10*b2[5] );
        }
};

typedef struct lessKey lessKey;
typedef map< block, vector< vector< int > >, lessKey > Block_Map;

struct CommuTask {
    int tag;
    int type;
    int other_pid;
    vector<int> shpFn ; 
};

typedef struct CommuTask CommuTask;

// since ParallelData[x] is a sorted data structure, the lowest pid 
// adjacency always  becomes the master.

int 
getMaster( int x ) {
    map< int , int >::iterator iter  =  ParallelData[x].begin();
    return (*iter).first;

}
void 
generate_keyphrase(char* target, const char* prefix, block& tpblock) {

    strcpy(target,prefix);  /* interior or boundary */

    switch( tpblock[1] ){
    case 1:
        strcat(target,"linear ");
        break;
    case 2:
        strcat(target,"quadratic ");
        break;
    case 3:
        strcat(target,"cubic ");
        break;
    case 4:
        strcat(target,"quartic ");
        break;
    }

    switch( tpblock[5] ){
    case 1:
        strcat(target,"tetrahedron ");
        break;
    case 2:
        strcat(target,"hexahedron ");
        break;
    case 3:
        if(!strcmp(prefix,"boundary "))
            strcat(target,"wedge triface ");
        else
            strcat(target,"wedge ");
        break;
    case 4:
        strcat(target,"wedge quadface ");
        break;
    case 5:
        if(!strcmp(prefix,"boundary "))
            strcat(target,"pyramid quadface ");
        else
            strcat(target,"pyramid ");
        break;
    case 6:
        strcat(target,"pyramid triface ");
        break;
    }
}


int
main( int argc,
      char* argv[] ) { 

    /* load the ien and ienb and X into arrays */ 

    int stepno;
    int igeombc;  /* file handle for geombc */
    int irestart; /* file handle for restart */
    int iarray[10];
    int ione=1;
    int itwo=2;
    int ithree=3;
    int iseven=7;
    int ieight=8;
    int isize;
    int nitems;
    int numProcs = 1;
    char* iformat="binary";
    char* oformat="binary";
    int SONFATH_VAR = 0;
    int old_format = 0;

    // process the command line arguments 
    int c;
    while (--argc > 0 && (*++argv)[0] == '-')  {
        while ((argc > 0) && (c = *++argv[0]))  {
            switch( c ) {
            case 'n':
                numProcs = atoi( argv[1] );
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'i':
                iformat = argv[1] ;
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'o':
                oformat = argv[1] ;
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'S':
                SONFATH_VAR = atoi( argv[1] ) ;
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'B':
                old_format = 1;
                break;
            }
        }
    }
    /* Lets us first load the dual and partition it */

   
    sprintf(dirname,"%d-procs_case/",numProcs);
#if !(defined WIN32)
    struct stat buf;
    if ( !stat( dirname, &buf ) &&  S_ISDIR( buf.st_mode ) )
#else
    struct _stat buf;
    char probdir[25];
    sprintf(probdir,"%d-procs_case",numProcs);
    int result = _stat( probdir, &buf );
    if ( (result == 0 ) &&  ( 01 & ( buf.st_mode >> _S_IFDIR ) ) ) 
#endif
        {
        cout << "Directory " << dirname <<" already exists " << endl;
        cout << "using the existing inputfiles" << endl;
        return 0 ;
    } else {
        cout << dirname << " does not exist or is unusable" << endl;
        cout << "creating a new one" << endl;
#if !(defined WIN32)
        unlink( dirname );
#else 
        _unlink( dirname );
#endif
    }
    // else if the case does not already exist

#if !(defined WIN32)
    if ( mkdir( dirname, 00755 ) ) {
#else 
    if ( _mkdir( dirname ) ) {
#endif
        cerr << "cannot create directory " << dirname << endl;
        exit(0);
    }

// Igor, 08/2012: directory structure for restarts (important for NumProcs > ~20000 ?)

    int idirstep = 512;    // Number of files in each directory
    int idirtrigger = 10;  // Number of procs to trigger separate directories

    if (numProcs > idirtrigger) {
      for (int i=0; i<(int)(numProcs/idirstep)+1; i++) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,i*idirstep);
#if !(defined WIN32)
    if ( mkdir( dirname2, 00755 ) ) {
#else
    if ( _mkdir( dirname2 ) ) {
#endif
        cerr << "cannot create directory " << dirname << endl;
        exit(0);
    }

      }
    }
    sprintf(filename,"%snumpe.in",dirname);
    ofstream npe(filename);
    npe << numProcs << endl ;
    npe.close();
    
    openfile( "geombc.dat.1", "read", &igeombc );
    vector< string > headers;
    Gather_Headers( &igeombc, headers );
    
    readheader( &igeombc, "number of interior elements", (void*)iarray, 
                 &ione, "integer", iformat );

    int numel = iarray[0];
    int* epart = new int [ numel];
    if ( numProcs > 1 ) {
        int nIfaces, spebc;
        readheader( &igeombc, "keyword xadj", (void*)iarray,
                     &itwo, "integer", iformat );
        numel = iarray[0];
        spebc = iarray[1];
        isize = numel + 1;
//        long int* xadj = new long int [ isize ];   // July 2012 Igor, see the kmetis.c and struct.h in phMetis for idxtype !
        int* xadj = new int [ isize ];
        readdatablock( &igeombc, "keyword xadj", (void*)xadj,
                        &isize, "integer", iformat );
    
        readheader( &igeombc, "keyword adjncy", (void*)iarray,
                     &ione, "integer", iformat );
        nIfaces = iarray[0];
        isize = 2*nIfaces;
        int* adjncy = new int [ isize ];
        readdatablock( &igeombc, "keyword adjncy", (void*)adjncy,
                        &isize, "integer", iformat );
    
        isize = numel;
        int* vwgt = new int [ numel ];
        readheader( &igeombc, "keyword vwgt", (void*)iarray,
                     &ione, "integer", iformat );
        readdatablock( &igeombc, "keyword vwgt", (void*)vwgt,
                        &isize, "integer", iformat );
    
        int wgtflag = 2;
        int numflag = 0;
        int edgecut;
        int options[10];
        options[0] = 0; 
        METIS_PartGraphKway( &numel, xadj, adjncy, vwgt, NULL, &wgtflag,
                             &numflag, &numProcs, options, &edgecut, epart);
    
        /* putting SPEBC elements together on the processor 0 */
        if ( spebc > 0 ) {
            int swap;
            for (int i=0; i<numel; i++)
                if ( vwgt[i] > 1 ) { 
                    swap = epart[i];
                    break;
                }
            for (int i=0; i<numel; i++) {
                if ( epart[i] == swap ) epart[i] = 0;
                else if ( epart[i] == 0 ) epart[i] = swap;
                if ( vwgt[i] == 0 ) epart[i] = 0;
            }
        }
    
        delete [] xadj;
        delete [] adjncy;
        delete [] vwgt;

    } else {

        for(int nel=0; nel < numel; nel++ ) epart[nel]=0;

    }

    
    /* now we have the epart array which is all the partition Info we need */

    int ndof=5; 
    ifstream numstart("numstart.dat");
    numstart >> stepno;
    numstart.close();

    readheader( &igeombc, "number of global modes", (void*)iarray, 
                 &ione, "integer", iformat );
    int nshgTot = iarray[0];
    int nshg = iarray[0] ;

    readheader( &igeombc, "number of variables", (void*)iarray, 
                 &ione, "integer", iformat );

    ndof = iarray[0] ;

    readheader( &igeombc, "maximum number of element nodes", (void*)iarray, 
                 &ione, "integer", iformat );
    int nenmax = iarray[0];

    readheader( &igeombc, "number of interior tpblocks", (void*)iarray, 
                 &ione, "integer", iformat  );
    int nblock = iarray[0];

    readheader( &igeombc, "number of nodes in the mesh", (void*)iarray, 
                 &ione, "integer", iformat );
    int numnp = iarray[0];

    readheader( &igeombc, "number of edges in the mesh", (void*)iarray, 
                 &ione, "integer", iformat );
    int numedges = iarray[0];

    vector< Block_Map >  ien( numProcs );

    int* vcount= new int [numProcs];
    int* ecount= new int [numProcs];
    int* fcount= new int [numProcs];

    for( int e=0; e < numProcs; e++ ) 
        vcount[e] = ecount[e] = fcount[e]= 1;

    for( int b = 0; b < nblock ; b++ ) {

        readheader( &igeombc, "connectivity interior?", (void*)iarray,
                     &iseven, "integer", iformat );
        block CurrentBlock ;
        for( int w=1; w< iseven; w++ ) 
            CurrentBlock.push_back( iarray[w] );
        isize = iarray[0]*iarray[3];
        int* ient = new int [ isize ];

        readdatablock( &igeombc, "connectivity interior?", (void*)ient, &isize,
                        "integer", iformat );

        readheader( &igeombc, "ien to sms?", (void*)iarray, &ione,
                     "integer", iformat );

        isize=iarray[0];
        int* ient_sms = new int [ isize ];
        readdatablock( &igeombc, "ien to sms?", (void*)ient_sms,
                        &isize, "integer", iformat );

        vector< int > element;
        for( int c=0; c < iarray[0] ;  c++ ) {
            for( int d=0; d < iarray[3] ; d++ ) {
                element.push_back( ient[d*iarray[0]+c] );
            }
            int gid = ient_sms[c];
            int pid = epart[ gid ];

            ien[pid][CurrentBlock].push_back( element );

            for(  vector<int>::const_iterator iterv = element.begin();
                  iterv != element.end();
                  iterv++ ) {
                int id = abs ( *iterv );
                bool test = false;

                map<int,map<int, int> >::const_iterator iterm=ParallelData.find( id );
                if ( iterm == ParallelData.end() ) test = true;
                else  if ( 0 == ((*iterm).second).count(pid)) test = true;
                if( test ) {
                    ParallelData[id][pid] = (id>numnp) ? 
                                                (id>numnp+numedges)? 
                                                    fcount[pid]++:ecount[pid]++ 
                                                : vcount[pid]++;
                }
            }
            element.clear();
        }
        delete [] ient;
        delete [] ient_sms;
    }

    /* At this point the values stored in vcount ecount and fcount are 1 more 
     * than  the number of vertices edges and faces respectively , sp we need 
     * to decrement  all of them by one so that we can directly added them to 
     * the lower order counts  in the next  step. 
     */
 
    for( int e=0; e < numProcs; e++ ) {
        --vcount[e] ; 
        --ecount[e] ; 
        --fcount[e];
    }

    /* our parallel data structure has indepenent numbering for vertices, faces
     * and edges, now we need to correct it and add the correct offsets , we
     * cannot write out ien before we make this correction 
     */

    for( map< int, map< int, int >  >::iterator iter = ParallelData.begin();
         iter != ParallelData.end();
         iter++ ) {
        if ( (*iter).first > numnp ) {         // it is not a vertex 
            if ( (*iter).first > ( numnp + numedges ) ) { // it is not a edge, 
                // must be a face
                for( map< int, int> :: iterator modeIter = (*iter).second.begin();
                     modeIter != (*iter).second.end();
                     modeIter++ ) 
                    (*modeIter).second = (*modeIter).second 
                                         + vcount[ (*modeIter).first ] 
                                         + ecount[ (*modeIter).first ] ;
            } else { // is an edge
                for( map< int, int> :: iterator modeIter = (*iter).second.begin();
                     modeIter != (*iter).second.end();
                     modeIter++ ) 
                    (*modeIter).second = (*modeIter).second 
                                         + vcount[ (*modeIter).first ];
            }
        }
    }

    delete [] vcount;
    delete [] ecount;
    delete [] fcount;

    /* after correction all local mode numbers are accurate ( hopefully ) and
     * are in fortran style ( starting 1 ) */
                        
    /* at this point, ien still has global numbering */
    /* we now use the parallel data structure to correct ien */

    for( int p = 0 ; p < numProcs; p++ ) {
        for( Block_Map::iterator iter = ien[p].begin();
             iter != ien[p].end();
             iter++ ) {
            for( vector<vector<int> >::iterator iter1 = ((*iter).second).begin();
                 iter1 != ((*iter).second).end();
                 iter1++ ) {
                for( vector<int>::iterator iter2 = (*iter1).begin();
                     iter2 != (*iter1).end();
                     iter2++ ) {
                    int id = *iter2 ;
                    int absid = abs( id );
                    map<int,int> ::const_iterator pIter = ParallelData[absid].find( p );
                    *iter2 = (absid/id)*((*pIter).second);
                }
            }
        }
    }

    /* now we can write out ien and delete that data structure. */
    /* ien is stored as ien[pid][block][iel][nshl] */
    /* we need to write this block by block for each pid */

    int magic_number = 362436;
    int* mptr = &magic_number;
    int fgeom;
    int frest;

    for( int p=0; p< numProcs; p++ ) {
        int numel=0;
        
        bzero( (void*)filename,255);

// Make sure you include the longer dir name here   ! Igor, 08/2012

    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);
        sprintf(filename,"%sgeombc.dat.%d",dirname2,p+1 );
      }
	else
	{
        sprintf(filename,"%sgeombc.dat.%d",dirname,p+1 );
	}
        openfile( filename, "write", &fgeom );

        for( vector< string >::iterator iter = headers.begin();
             iter != headers.end();
             iter++ ) {
            writestring( &fgeom, (*iter).c_str() );
            writestring( &fgeom, "\n" );
        }
        
#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);
        sprintf(filename,"%sascii_out.%d",dirname2,p+1);
	}
	else
	{
        sprintf(filename,"%sascii_out.%d",dirname,p+1);
	}
        ofstream fascii( filename );
        fascii <<"Interior Element Modal Connectivity:"<<endl;
        fascii <<"------------------------------------"<<endl;
#endif

        isize = 1;
        nitems = 1;
        iarray[ 0 ] = 1;
        writeheader( &fgeom, "byteorder magic number ", 
                      (void*)iarray, &nitems, &isize, "integer", oformat );

        nitems = 1;
        writedatablock( &fgeom, "byteorder magic number ",
                         (void*)mptr, &nitems, "integer", oformat );

        bzero( (void*)filename,255 );
        sprintf( filename, "number of interior tpblocks : < 0 > %d \n",
                 ien[p].size()  );
        writestring( &fgeom, filename );

        bzero( (void*)filename,255 );
        sprintf( filename, "number of global modes : < 0 > %d \n" ,nshgTot );
        writestring( &fgeom, filename );

        bzero( (void*)filename, 255 );
        sprintf( filename,"number of nodes in the mesh : < 0 > %d \n",numnp);
        writestring( &fgeom, filename );
    
            for( Block_Map::iterator bIter = ien[p].begin();
             bIter != ien[p].end();
             bIter++ ) {
            block CurrentBlock = (*bIter).first;
            vector< vector< int > > blockIEN = (*bIter).second;
            numel += blockIEN.size();
            isize = blockIEN.size()*CurrentBlock[2];
            iarray[ 0 ] = blockIEN.size();        

#if defined ( DEBUG )
            fascii << endl;
            fascii << "********************"<< endl;
#endif

            int xct = 1;
            for( block::iterator ibtr = CurrentBlock.begin();
                 ibtr != CurrentBlock.end();
                 ibtr++ ) {
                iarray[ xct++ ] = *ibtr;

#if defined ( DEBUG )
                fascii << *ibtr <<" ";
#endif
            }

            generate_keyphrase( keyphrase, "connectivity interior ", 
                                CurrentBlock );
            writeheader( &fgeom, keyphrase, (void*)iarray, &xct, &isize,
                          "integer", oformat );
   
#if defined ( DEBUG )
            fascii << endl;
            fascii << "********************"<< endl;
#endif
            /* now we need to create the fortran style array for this block of
             * ien, so need to invert it */

            int* ient = new int [ blockIEN.size()*CurrentBlock[2] ];
            for( int h=0; h < blockIEN.size(); h++ ) {
                for( int y=0; y< CurrentBlock[2]; y++ ) {
                    ient[ y*blockIEN.size() + h] = blockIEN[h][y];
#if defined ( DEBUG )
                    fascii << ient[ y*blockIEN.size() + h ] <<" ";
#endif
                }
                ((*bIter).second)[h].clear();

#if defined ( DEBUG )
                fascii << endl;
#endif
            }
            ((*bIter).second).clear();
            nitems = blockIEN.size()*CurrentBlock[2];
            writedatablock( &fgeom, keyphrase, (void*)( ient ),
                             &nitems, "integer", oformat ); 
            
            delete [] ient;
            blockIEN.clear();
        }// loop over blocks with in a processor
        ien[p].clear();

        bzero( (void*)filename, 255 );
        sprintf(filename ,"number of interior elements : < 0 > %d \n",numel);
        writestring( &fgeom, filename );

        closefile( &fgeom, "write" );
#if defined ( DEBUG )
        fascii.close();
#endif
    }// loop over procs 
    //ien.clear();
   
    /* co-ords */

    readheader( &igeombc, "co-ordinates", (void*)iarray, &itwo, "double",
                 iformat );
    isize =  iarray[0]*iarray[1];
    double* xloc = new double [ isize ];
    readdatablock( &igeombc, "co-ordinates", (void*)xloc, &isize,
                    "double", iformat );
    vector< map<int, vector< double > > > Xpart( numProcs );

    for( int x=1; x < iarray[0]+1; x++ ){
        for( map< int, int>::const_iterator iter = ParallelData[x].begin();
             iter != ParallelData[x].end();
             iter++ ) {
            for( int s=0; s < 3 ; s++ ) 
                Xpart[(*iter).first][(*iter).second].push_back( xloc[x-1+s*iarray[0]] );
        }
    }
    delete [] xloc;

    // if we need to calculate sonfath later we will need to have vnumnp[i]
    int* vnumnp;
    if ( SONFATH_VAR ) vnumnp = new int [ numProcs ];

    for( int p=0; p< numProcs; p++ ) {

        bzero((void*)filename, 255 );
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sgeombc.dat.%d",dirname2,p+1 );
	}
	else
        sprintf(filename,"%sgeombc.dat.%d",dirname,p+1 );

        openfile( filename, "append", &fgeom );
        
        if ( SONFATH_VAR ) vnumnp[p] = Xpart[p].size();

        bzero((void*)filename, 255 );
        sprintf(filename,"number of nodes : < 0 > %d \n", Xpart[p].size() );
        writestring( &fgeom, filename );

        bzero((void*)filename, 255 );
        sprintf(filename,"number of processors : < 0 > %d \n", numProcs);
        writestring( &fgeom, filename );
        
        double* xf = new double [ Xpart[p].size() * 3 ];
        for(   map< int, vector< double > >::iterator mIter = Xpart[p].begin();
               mIter != Xpart[p].end();
               mIter++ )
            for( int f=0; f<3; f++) 
                xf[ f*Xpart[p].size()+ (*mIter).first - 1] = (*mIter).second[f];
            
        isize = Xpart[p].size() * 3 ;
        nitems = 2;
        iarray[ 0 ] = Xpart[p].size();
        iarray[ 1 ] = 3;
        writeheader( &fgeom, "co-ordinates", (void*)iarray, &nitems, &isize,
                      "double", oformat );

        nitems = Xpart[p].size() * 3;
        writedatablock( &fgeom, "co-ordinates", (void*)( xf ), &nitems,
                         "double", oformat );

#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);   // Igor, August 2012
        sprintf( filename, "%sascii_out.%d",dirname2, p+1 );
	}
	else
        sprintf( filename, "%sascii_out.%d",dirname, p+1 );
        FILE* cfascii=fopen(filename,"a");
        fprintf(cfascii,"\n");
        fprintf(cfascii, "Nodal coordinates:\n");
        fprintf(cfascii, "------------------\n");
        int numnp_tmp = Xpart[p].size();
        for( int c=0; c< numnp_tmp ; c++ )
            fprintf( cfascii,"%f %f %f \n", 
                     xf[ c ], xf[ c + numnp_tmp ],xf[c+2*numnp_tmp] );
        fclose( cfascii );
#endif          

        delete [] xf ;
        Xpart[p].clear();
        closefile( &fgeom, "append" );
    }
    Xpart.clear();

   
    /* let us take care of Essential BCs.*/
    // BCs are in the indirect numbering of nBC so only nBC needs to be sorted 
    // according to nshg
    // BC abd iBC are numbered as numpbc which us given by nBC.

    int numEBC = ndof + 7;
    readheader( &igeombc, "bc mapping array", (void*)iarray, &ione, "integer",
                 iformat );
    nshg = iarray[0];
    int* nBC = new int [ nshg ];
    readdatablock( &igeombc, "bc mapping array", (void*)nBC, &nshg, "integer",
                    iformat );

    readheader( &igeombc, "bc codes array", (void*)iarray, &ione,"integer",
                 iformat );
    int numpbc = iarray[0];
    int* iBC = new int [ numpbc ];
    readdatablock( &igeombc, "bc codes array", (void*)iBC, &numpbc,
                    "integer", iformat  );

    readheader( &igeombc, "boundary condition array", (void*)iarray,
                 &ione, "double", iformat );
    double* BC = new double  [ numpbc*numEBC ];
    readdatablock( &igeombc, "boundary condition array", (void*)BC,
                    &iarray[0], "double", iformat );

    int* inumpbc = new int [ numProcs ];
    for( int c=0; c < numProcs; c++ ) inumpbc[ c ] = 0;
    vector< map< int, int> > nBCpart( numProcs );
    vector< vector< int> > iBCpart ( numProcs );
    vector< vector< double> > BCpart( numProcs );
    
    for( int x=1; x < nshg+1 ; x++ ) {
        if ( nBC[x-1] > 0 ) {
            for( map< int, int>::const_iterator iter = ParallelData[x].begin();
                 iter != ParallelData[x].end();
                 iter++ ) {
                nBCpart[ (*iter).first ][ (*iter).second ] 
                    = inumpbc[ (*iter).first ]++; 
                iBCpart[ (*iter).first ].push_back( iBC[ nBC[x-1] - 1 ] );
                for( int g=0; g< numEBC; g++ )
                    BCpart[ (*iter).first ].push_back(BC[nBC[x-1]-1+g*numpbc]); 
            }
        } else {
            for( map< int, int>::const_iterator iter = ParallelData[x].begin();
                 iter != ParallelData[x].end();
                 iter++ ) 
                nBCpart[ (*iter).first ][ (*iter).second ] = -1;
        }
    }
   
    delete [] nBC;
    delete [] iBC;
    delete [] BC;

    // we will calculate a maxnshg here for later use with sonfath.
    
    int maxnshg =0;

    for( int p=0; p< numProcs; p++ ) {

        bzero( (void*)filename, 255 );
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sgeombc.dat.%d",dirname2,p+1 );
	}
	else
        sprintf(filename,"%sgeombc.dat.%d",dirname,p+1 );

        openfile( filename, "append", &fgeom );

#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);   // Igor, August 2012
        sprintf( filename, "%sascii_out.%d",dirname2, p+1 );
	}
	else
        sprintf( filename, "%sascii_out.%d",dirname, p+1 );
        FILE* cfascii = fopen( filename, "a" );
        fprintf(cfascii,"\nBoundary Condition Mapping Array (nBC):\n" ); 
        fprintf(cfascii,"---------------------------------------\n" );
#endif
        if ( (SONFATH_VAR) && (maxnshg < nBCpart[p].size()) ) 
            maxnshg = nBCpart[p].size();
        bzero( (void*)filename, 255 );
        sprintf( filename,  "maximum number of element nodes : < 0 > %d\n",
                 nenmax );
        writestring( &fgeom, filename );
        
        bzero( (void*)filename, 255 );
        sprintf(filename,  "number of modes : < 0 > %d\n", nBCpart[p].size() );
        writestring( &fgeom, filename );
        
        bzero( (void*)filename, 255 );
        sprintf( filename,
                 "number of shapefunctions soved on processor : < 0 > %d\n", 
                 nBCpart[p].size() );
        writestring( &fgeom, filename );
        
        bzero( (void*)filename, 255 );
        sprintf( filename, 
                 "number of nodes with Dirichlet BCs : < 0 > %d\n",
                 inumpbc[ p ] );
        writestring( &fgeom, filename );
        
        isize =  nBCpart[p].size();
        nitems = 1 ;
        iarray[ 0 ] = nBCpart[p].size();
        writeheader( &fgeom, "bc mapping array", (void*)iarray, &nitems,
                      &isize, "integer", oformat );
        
        int* inBC= new int [ nBCpart[p].size() ];
        int count = 0;
        for( map< int, int>::iterator iter = nBCpart[p].begin();
             iter != nBCpart[p].end();
             iter++ ) {
            inBC[ count++ ] = ((*iter).second + 1) ;

#if defined ( DEBUG )
            fprintf( cfascii, "%d\n", inBC[ count - 1 ] );
#endif

        }
        writedatablock( &fgeom, "bc mapping array", (void*)( inBC ), &count,
                         "integer", oformat );
        delete [] inBC;
        nBCpart[p].clear();

#if defined ( DEBUG )
        fprintf(cfascii,"\nBoundary Condition Codes Array (iBC):\n" ); 
        fprintf(cfascii,"-------------------------------------\n" );
#endif
        isize = inumpbc[p];
        nitems = 1;
        iarray[ 0 ] = inumpbc[ p ];

        writeheader( &fgeom, "bc codes array ", (void*) iarray, &nitems,
                      &isize, "integer", oformat );

        int* iiBC = new int [ inumpbc[p] ];
        if ( inumpbc[p] != iBCpart[p].size() ) 
            cerr<<"Error 1" << __LINE__ <<endl;
        for( int i = 0; i< inumpbc[p]; i++ ) {
            iiBC[i] = iBCpart[p][i];
#if defined ( DEBUG )
            fprintf( cfascii, "%d\n", iiBC[ i ] );
#endif
        }
        iBCpart[p].clear();
        nitems = inumpbc[ p ];
        writedatablock( &fgeom, "bc codes array ", (void*)( iiBC ), &nitems,
                         "integer", oformat );
        delete [] iiBC;
        
#if defined ( DEBUG )
        fprintf(cfascii,"\nBoundary Condition Array (BC):\n" ); 
        fprintf(cfascii,"--------------------------------\n" );
#endif
        isize = inumpbc[p]*numEBC;
        nitems = 1;
        iarray[ 0 ] = inumpbc[p]*numEBC;
        writeheader( &fgeom, "boundary condition array ", (void*) iarray, 
                      &nitems, &isize, "double", oformat );

        double* BCf = new double [inumpbc[p]*numEBC];
        for( int a = 0 ; a < inumpbc[p]; a++){
            for( int b =0; b < numEBC; b++ ){
                BCf[b*inumpbc[p]+a] = BCpart[p][a*numEBC+b];

#if defined ( DEBUG )
                fprintf( cfascii,"%.3f ",BCpart[p][a*numEBC+b]);
#endif

            }

#if defined ( DEBUG )
            fprintf( cfascii, "\n" );
#endif

        }
        BCpart[p].clear();
        nitems = inumpbc[p]*numEBC; 
        writedatablock( &fgeom, "boundary condition array ", (void*)( BCf ), 
                         &nitems, "double", oformat );

        delete [] BCf;
        closefile( &fgeom, "append" );

#if defined ( DEBUG )
        fclose( cfascii );
#endif
    }

    iBCpart.clear();
    BCpart.clear();
        

    /* done writing ebcs */

    /* ienb and BCB and iBCB */

    
    int numNBC = ndof + 1;


    readheader( &igeombc, "number of boundary tpblocks", (void*)iarray, 
                 &ione, "integer", iformat );
    nblock = iarray[0];

    //ienb -- Boundary Element Nodal Connectivity.
    //ibcb -- Natural Boundary Condition Codes
    //bcb  -- Natural Boundary Condition Values

    vector< map< block, vector< int >, lessKey > > iBCBpart( numProcs );
    vector< map< block, vector< double >, lessKey > > BCBpart( numProcs );

    for( int b = 0; b < nblock ; b++ ) {
        
        readheader( &igeombc, "connectivity boundary?", (void*)iarray, &ieight,
                     "integer", iformat );

        block CurrentBlock;
        for( int w=1; w< ieight; w++ ) CurrentBlock.push_back( iarray[w] );
        
        isize = iarray[0]*iarray[3];
        int* ient = new int [ isize ];
        readdatablock( &igeombc, "connectivity boundary?", (void*)ient,
                        &isize, "integer", iformat );
        
        readheader( &igeombc, "ienb to sms?", (void*)iarray, &ione,
                     "integer", iformat );

        isize=iarray[0];
        int* ient_sms = new int [ isize ];
        readdatablock( &igeombc, "ienb to sms?", (void*)ient_sms,
                        &isize, "integer", iformat );


        readheader( &igeombc, "nbc codes?", (void*)iarray, &ieight,
                     "integer", iformat );

        isize = iarray[0]*2;
        int* iBCB = new int [ isize ];
        readdatablock( &igeombc, "nbc codes?", (void*)iBCB,
                        &isize, "integer", iformat );

        readheader( &igeombc, "nbc values?", (void*)iarray, &ieight,
                     "double", iformat );

        isize = iarray[0]*numNBC;
        double* BCB = new double [ isize ];
        readdatablock( &igeombc,"nbc values?",(void*)BCB, &isize,
                        "double", iformat );


        for( int c=0; c < iarray[0] ;  c++ ) {
            vector< int > element;
            vector< int > element2;
            for( int d=0; d < iarray[3] ; d++ )
                element.push_back( ient[d*iarray[0]+c] );

            int gid = ient_sms[c];
            int pid = epart[ gid ];
            int _local_id ;
            for( vector< int >::iterator iter = element.begin();
                 iter != element.end();
                 iter++ ) {
                int id = abs( *iter );
                map<int,int> ::const_iterator pIter = ParallelData[id].find( pid );
                _local_id = (*pIter).second;
                if ( *iter < 0 ) _local_id = -1*_local_id;
                element2.push_back( _local_id );
            }
            element.clear();
            ien[pid][CurrentBlock].push_back( element2 );
            element2.clear();

            iBCBpart[pid][CurrentBlock].push_back( iBCB[c] );
            iBCBpart[pid][CurrentBlock].push_back( iBCB[c + iarray[0]] );

            for( int o=0; o< numNBC; o++) 
                BCBpart[pid][CurrentBlock].push_back( BCB[ c + iarray[0]*o ] ); 

        }
          
        delete [] ient;
        delete [] ient_sms;
        delete [] iBCB;
        delete [] BCB;
    }

    delete [] epart;

    /* format and write out */

    for( int p=0; p< numProcs; p++ ) {
        
        bzero( (void*)filename,255);
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sgeombc.dat.%d",dirname2,p+1 );
	}
	else
        sprintf(filename,"%sgeombc.dat.%d",dirname,p+1 );
        openfile( filename , "append", &fgeom );
        
        bzero( (void*)filename,255);
        sprintf( filename,"number of boundary tpblocks : < 0 > %d\n", 
                 ien[p].size() );
        writestring( &fgeom, filename );

#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(p/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sascii_out.%d",dirname2,p+1);
	}
	else
        sprintf(filename,"%sascii_out.%d",dirname,p+1);
        ofstream fascii( filename, ios::app );
        fascii << endl;
        fascii <<"Boundary Element Modal Connectivity:"<<endl;
        fascii <<"------------------------------------"<<endl;
#endif
        int numelb = 0;
        for( Block_Map::iterator bIter = ien[p].begin();
             bIter != ien[p].end();
             bIter++ ) {
            block CurrentBlock = (*bIter).first;
            vector< vector< int > > blockIEN = (*bIter).second;
            numelb += blockIEN.size();
             
            isize = blockIEN.size()*CurrentBlock[2];
            if (old_format) isize = blockIEN.size()*CurrentBlock[2]
                                + blockIEN.size()*2
                                +blockIEN.size()*numNBC;
            iarray[ 0 ] = blockIEN.size();

#if defined ( DEBUG )
            fascii << endl;
            fascii <<"************************"<< endl;    
#endif
            int xct = 1;
            for( block::iterator ibtr = CurrentBlock.begin();
                 ibtr != CurrentBlock.end();
                 ibtr++ ) {
                iarray[ xct++ ] = *ibtr;
#if defined ( DEBUG )
                fascii << *ibtr <<" ";
#endif
            }
            generate_keyphrase( keyphrase, "connectivity boundary ", 
                                CurrentBlock );
            writeheader( &fgeom, keyphrase, (void*)iarray, &xct, &isize,
                          "integer", oformat );

#if defined ( DEBUG )
            fascii << endl;
            fascii <<"************************"<< endl;    
#endif

            /* now we need to create the fortran style array for this block of
             * ien, so need to invert it */

            int* ient = new int [ blockIEN.size()*CurrentBlock[2] ];
            for( int h=0; h< blockIEN.size(); h++ ) {
                for( int y=0; y< CurrentBlock[2]; y++ ) {
                    ient[ y*blockIEN.size() + h] = blockIEN[h][y];
#if defined ( DEBUG )
                    fascii << ient[ y*blockIEN.size() + h ] <<" ";
#endif
                }
                ((*bIter).second)[h].clear();
#if defined ( DEBUG )
                fascii << endl;
#endif
            }
          
            ((*bIter).second).clear();
            nitems = blockIEN.size()*CurrentBlock[2];
            writedatablock( &fgeom, keyphrase, (void*)( ient ),
                             &nitems, "integer", oformat ); 

            delete [] ient;

            /* need to also write IBCB and BCB now */
#if defined ( DEBUG )
            fascii << "------------------------" << endl;
            fascii << "iBCB for the above block" << endl;
            fascii << "------------------------" << endl;
#endif

            int* iBCBf = new int [ blockIEN.size() * 2 ];
            for( int u=0; u< blockIEN.size(); u++ )  {
                iBCBf[ u ] = iBCBpart[p][CurrentBlock][u*2];
                iBCBf[ u + blockIEN.size() ] = iBCBpart[p][CurrentBlock][u*2+1];

#if defined ( DEBUG )
                fascii << iBCBf[u] <<"  "<<iBCBf[ u + blockIEN.size() ]<< endl;
#endif

            }

            iBCBpart[p][CurrentBlock].clear();

            if (!old_format) {
	            isize = blockIEN.size()*2;
                nitems = xct;
                generate_keyphrase( keyphrase, "nbc codes ", CurrentBlock );
                writeheader( &fgeom, keyphrase, (void*)iarray, &xct, &isize,
                              "integer", oformat );
	        }

            nitems = blockIEN.size()*2;
            writedatablock( &fgeom, keyphrase, (void*)( iBCBf ), &nitems,
                             "integer", oformat );

            delete [] iBCBf;

#if defined ( DEBUG )
            fascii.precision( 8 );
            fascii << "------------------------" << endl;
            fascii << "BCB for the above block" << endl;
            fascii << "------------------------" << endl;
            fascii.precision( 8 );
#endif
            double* BCBf = new double [ blockIEN.size() * numNBC ];
            for( int u=0; u< blockIEN.size(); u++ ) { 
                for( int v=0; v< numNBC; v++ ) { 
                    BCBf[v* blockIEN.size()+u] 
                        = BCBpart[p][CurrentBlock][u*numNBC+v];

#if defined ( DEBUG )
                    fascii << BCBf[v* blockIEN.size()+u] << " ";
#endif
                }

#if defined ( DEBUG )
                fascii << endl ;
#endif
            }

            BCBpart[p][CurrentBlock].clear();

            if (!old_format) {
	            isize = blockIEN.size()*numNBC;
                nitems = xct;
                generate_keyphrase( keyphrase, "nbc values ", CurrentBlock );
                writeheader( &fgeom, keyphrase, (void*)iarray, &xct, &isize,
                              "double", oformat );
	        }

            nitems = blockIEN.size()*numNBC;        
            writedatablock( &fgeom, keyphrase, (void*)( BCBf ), &nitems,
                             "double", oformat );
            delete [] BCBf;
                    
            blockIEN.clear();
        }// loop over blocks with in a processor
        ien[p].clear();

        bzero((void*)filename , 255 );
        sprintf(filename,"number of boundary elements : < 0 > %d \n", numelb );
        writestring( &fgeom, filename );

        closefile( &fgeom, "append" );

#if defined ( DEBUG )
        fascii.close();
#endif
    }// loop over procs 
    ien.clear();


    /* now starts commu and periodicity */
    /* We need to generate ilwork and ncorp and the iper array for each
     * processor . To achieve this we first read in the complete iper array and
     * then go through "ParallelData" iterating all the global shapefunctions
     * and then generate a consistant communication trace for all 
     * for on processor periodiciy we just dump an array */

    readheader( &igeombc, "periodic masters array?", (void*)iarray, &ione, 
                 "integer", iformat );
    int* periodic = new int [ nshg ];
    readdatablock( &igeombc,"periodic masters array?", (void*)periodic, &nshg,
                    "integer", iformat );
    vector< map< int, int > > PeriodicPart( numProcs );
    
    CommuTask*** rtask ;
    CommuTask*** stask;

    rtask = new CommuTask** [numProcs];   
    stask = new CommuTask** [numProcs];   
    for( int b=0; b< numProcs; b++ ){
        rtask[b] = new CommuTask* [numProcs];   
        stask[b] = new CommuTask* [numProcs];   
        for( int c=0; c< numProcs; c++ ){
            rtask[b][c] = NULL;
            stask[b][c] = NULL;
        }
    }
    int newtag=1;
    bool isPERIODIC = false ;
    for( int x=1; x < nshg+1; x++ ) {

        // first we need to fill a 0 in all the iper entries of x and its images 
        // since this will be the case even if the node has no periodicity at all 

        isPERIODIC = false ;

        for( map< int, int >::const_iterator iter = ParallelData[x].begin();
             iter!= ParallelData[x].end();
             iter++ ) 
            PeriodicPart[ (*iter).first ][ (*iter).second ] = 0;

             
        int master_pid = getMaster( x ); // assuming no periodicity
        map<int,int>::const_iterator pIter = ParallelData[x].find(master_pid);
        int Global_Master = (*pIter).second;

        if ( periodic[ x-1 ] > 0 ){  

            // if there is periodicity, then the real master node and processor are 
            // not that of the current mode but of the master mode,  and need 
            // updating...also now the image which shares this master_pid should get 
            // an iper and not a commu.

            master_pid = getMaster( periodic[x-1] );

            // master_pid is the processor on which the master is real ( present ) 

            pIter = ParallelData[periodic[x-1]].find(master_pid);
            Global_Master = (*pIter).second;  

            // Global_Master is the processor local shape function number of the master
            // mode on the master processor

            isPERIODIC = true;
        }
        // periodicity has been taken care of 
        // now we do the simple commu to the master_pid 
         
        for( map< int, int>::const_iterator iter = ParallelData[x].begin();
             iter != ParallelData[x].end();
             iter++ )             
            if ( (*iter).first != master_pid ) {

                int sender = (*iter).first;
                int receiver = master_pid;
                 
                // if we ever add a send task from A to B then we also immediately add a 
                // receive task from B to A so its enuf if we check one of them 
                if ( !stask[sender][receiver] ) {
                    stask[sender][receiver] = new CommuTask;
                    rtask[receiver][sender] = new CommuTask;
                    stask[sender][receiver]->tag = newtag;
                    rtask[receiver][sender]->tag = newtag++;
                    stask[sender][receiver]->type = 0;
                    rtask[receiver][sender]->type = 1;
                    stask[sender][receiver]->other_pid = receiver;
                    rtask[receiver][sender]->other_pid = sender;
                }
                 
                stask[sender][receiver]->shpFn.push_back( (*iter).second );
                rtask[receiver][sender]->shpFn.push_back( Global_Master );

            } else {
                // we  have an image of the original slave which is on the master
                // processor of the master node for this image we need to fill iper
                // instead of commu. ( ofcourse only for periodic modes )

                if ( isPERIODIC ) {
                    PeriodicPart[ (*iter).first ][ (*iter).second ] 
                        = Global_Master;
                }           
            }
    }

    delete [] periodic;

    // now to generate ILwork and Periodicity

    vector<int>  ilwork;

#if defined ( DEBUG)
    sprintf(filename,"%silwork.info",dirname );
    ofstream ilwf( filename );
#endif

    for( int a=0; a< numProcs; a++ ){ //outer loop over procs 
        ilwork.push_back( 0 );

#if defined ( DEBUG )
        ilwf << a+1 <<" receives  ";
#endif

        for( int b=0; b< numProcs; b++ )  { // inner loop 1 over procs
            if( rtask[a][b] ) {
                ilwork[0]++;
                ilwork.push_back( rtask[a][b]->tag );
                ilwork.push_back( rtask[a][b]->type );
                ilwork.push_back( rtask[a][b]->other_pid+1 );
                ilwork.push_back( rtask[a][b]->shpFn.size() );

#if defined ( DEBUG )
                ilwf << " from "   << rtask[a][b]->other_pid+1 
                     << " tag "    << rtask[a][b]->tag 
                     << " numSeg " << rtask[a][b]->shpFn.size() << endl ; 
#endif

                for( vector<int>::iterator iter = rtask[a][b]->shpFn.begin();
                     iter != rtask[a][b]->shpFn.end();
                     iter++ ) {
                    ilwork.push_back( *iter );
                    ilwork.push_back( 1 );

#if defined ( DEBUG )
                    ilwf << *iter << " ,";
#endif

                }

#if defined ( DEBUG )
                ilwf << endl;
#endif

                delete rtask[a][b];

            } 

        } // end inner loop 1 over procs

#if defined ( DEBUG )

        ilwf << endl;
        
        ilwf << a+1 << " sends  ";
#endif

        for( int b=0; b< numProcs; b++ ) { // inner loop 2 over procs
            if( stask[a][b] ) {
                ilwork[0]++;
                ilwork.push_back( stask[a][b]->tag );
                ilwork.push_back( stask[a][b]->type );
                ilwork.push_back( stask[a][b]->other_pid+1 );
                ilwork.push_back( stask[a][b]->shpFn.size() );

#if defined ( DEBUG )
                ilwf << " to "   << stask[a][b]->other_pid+1 
                     << " tag "    << stask[a][b]->tag
                     << " numSeg " << stask[a][b]->shpFn.size() << endl;
#endif

                for( vector<int>::iterator iter = stask[a][b]->shpFn.begin();
                     iter != stask[a][b]->shpFn.end();
                     iter++ ) {
                    ilwork.push_back( *iter );
                    ilwork.push_back( 1 );

#if defined ( DEBUG )
                    ilwf << *iter << " ,";
#endif

                }

#if defined ( DEBUG )
                ilwf << endl;
#endif

                delete stask[a][b];
            } 
        } // end inner 2loop over procs

#if defined ( DEBUG )
        ilwf << endl;
        ilwf << "Size of ILwork for "<< a+1 << " : " << ilwork.size() << endl;
        if ( ilwork.size() > 0 ) {
            for( vector<int>::iterator iter = ilwork.begin();
                 iter != ilwork.end();
                 iter++ ) 
                ilwf << *iter <<" ";

            ilwf << endl;
            ilwf << endl;
        }
#endif

        delete [] rtask[a];
        delete [] stask[a];

        bzero( (void*)filename, 255 );
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sgeombc.dat.%d",dirname2,a+1 );
	}
	else
        sprintf(filename,"%sgeombc.dat.%d",dirname,a+1 );
        openfile( filename, "append", &fgeom );

        bzero( (void*)filename, 255 );
        sprintf(filename,"size of ilwork array : < 0 > %d\n",ilwork.size() );
        writestring( &fgeom, filename );

        isize = ilwork.size();
        nitems = 1;
        iarray[ 0 ] = ilwork.size();
        writeheader( &fgeom, "ilwork ", (void*)iarray,
                      &nitems, &isize, "integer", oformat );

        int* filwork = new int [ilwork.size()];
        for( int r=0; r< ilwork.size(); r++ ) filwork[r] = ilwork[r];
        nitems = ilwork.size();
        writedatablock( &fgeom, "ilwork ", (void*)(filwork), &nitems,
                         "integer", oformat );
        ilwork.clear();
        delete [] filwork;

        isize = PeriodicPart[a].size();
        nitems = 1 ;
        iarray[ 0 ] = PeriodicPart[a].size();
        writeheader( &fgeom, "periodic masters array ", (void* )iarray,
                      &nitems, &isize, "integer", oformat );

        int* fper= new int [ PeriodicPart[a].size() ];
        int count = 0;
        for( map< int, int>::iterator iter = PeriodicPart[a].begin();
             iter != PeriodicPart[a].end();
             iter++ ) 
            fper[ count++ ] = (*iter).second;

#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf( filename, "%sascii_out.%d",dirname2, a+1 );
	}
	else
        sprintf( filename, "%sascii_out.%d",dirname, a+1 );
        ofstream fascii( filename, ios::app );
        fascii <<"\n------------------------"<<endl;
        fascii <<"Periodic Partners Array:"<<endl;
        fascii <<"------------------------"<<endl;
        for( int shg=0; shg < count ; shg++ ) 
            fascii << fper[ shg ] << endl;
        fascii.close();
#endif
        
        nitems = count;
        writedatablock( &fgeom, "periodic masters array ", (void*)( fper ),
                         &nitems, "integer", oformat );

        delete [] fper;
        PeriodicPart[a].clear();
        closefile( &fgeom , "append" );

    } // end outer loop over procs

#if defined ( DEBUG )
    ilwf.close();
#endif

    PeriodicPart.clear();

    delete [] rtask;
    delete [] stask;

    // write ncorp
    // generating ncorp 
    // ncorp is our map between the "partition local" and the global (sms) 
    // numbering of modes  ncorp[ proc ][ local number ] = sms_number

    vector< map< int, int > > ncorp( numProcs ) ;
    for( int x=1; x < nshg+1; x++ ) {
        for( map< int, int>::const_iterator iter=ParallelData[x].begin();
             iter != ParallelData[x].end();
             iter++ )
            ncorp[(*iter).first][(*iter).second] = x ;
    }

// memLS output

       for( int a=0; a < numProcs ; a++ ) {
         bzero( (void*)filename, 255 );
         ofstream myfileltg;
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf( filename, "%sltg.dat.%d",dirname2, a+1 );
        }
        else
        sprintf( filename, "%sltg.dat.%d",dirname, a+1 );
         myfileltg.open (filename, ios::out);
         myfileltg << nshg << endl;
         myfileltg << ncorp[a].size() << endl;
         for( map<int, int>::iterator niter = ncorp[a].begin(); niter != ncorp[a].end(); niter++ ) {
            myfileltg << (*niter).second << endl;
         }
         myfileltg.close();
      }

// -------------------------

    // SONFATH stuff

    if ( SONFATH_VAR > 0 ) {
        int* ncvec = new int [ maxnshg*numProcs ];
        for( int x=0; x < numProcs*maxnshg ; x++ ) ncvec[ x ] =0;

        for ( int p = 0; p < numProcs; p++ )
            for( map<int, int>::iterator niter = ncorp[p].begin();
                 niter != ncorp[p].end();
                 niter++ )
                ncvec[ numProcs * ((*niter).first-1) + p ] = (*niter).second;
            
        int* nsons = new int [ SONFATH_VAR ];
        int* ifath = new int [ maxnshg*numProcs ];

        sonfath( ncvec, vnumnp, ifath, nsons, 
                  &SONFATH_VAR, &numProcs, &nshgTot, &maxnshg );

        delete [] vnumnp ;
        delete [] ncvec;

        // now we have ifath and nsons for the whole mesh, need to split it 
        // and write out.

        int* ifath1d;

        for( int a=0; a < numProcs ; a++ ) {

            bzero( (void*)filename, 255 );
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sgeombc.dat.%d",dirname2,a+1 );
	}
	else
            sprintf( filename, "%sgeombc.dat.%d",dirname,a+1 );
            openfile( filename, "append", &fgeom );

            bzero( (void*)filename, 255 );
            sprintf( filename,"number of father-nodes : < 0 > %d\n", 
                     SONFATH_VAR );
            writestring( &fgeom, filename );

            isize = SONFATH_VAR;
            nitems = 1;
            iarray[ 0 ] = SONFATH_VAR;
            writeheader( &fgeom, "number of son-nodes for each father ",
                          (void*)iarray, &nitems, &isize, "integer", 
                          oformat );
            nitems = SONFATH_VAR;
            writedatablock( &fgeom, 
                             "number of son-nodes for each father ", 
                             (void*)( nsons ), &nitems, "integer", 
                             oformat );
                
            ifath1d = new int [ ncorp[a].size() ];
            for(int j=0; j< ncorp[a].size() ; j++) 
                ifath1d[j]=ifath[j*numProcs+a];

            isize = ncorp[a].size();
            nitems = 1;
            iarray[ 0 ] = ncorp[a].size();
            writeheader( &fgeom, "keyword ifath ", (void*)iarray, &nitems,
                          &isize, "integer", oformat );
            nitems = ncorp[a].size();
            writedatablock( &fgeom, "keyword ifath ",
                             (void*)( ifath1d ), &nitems, "integer", 
                             oformat );

#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
            sprintf( filename, "%sascii_out.%d",dirname2, a+1 );
	}
	else
            sprintf( filename, "%sascii_out.%d",dirname, a+1 );
            ofstream fascii( filename, ios::app );
            fascii <<"\n-------------------------------"<< endl;
            fascii <<"IFATH array "<<endl ;
            fascii <<"-------------------------------"<< endl;
            for( int j=0; j < ncorp[a].size() ; j++ ) 
                fascii << ifath1d[j] << endl;
            fascii <<"-------------------------------"<< endl;
            fascii.close();
#endif
                
            closefile( &fgeom, "append" );
            delete [] ifath1d;
        }

        delete [] nsons;
        delete [] ifath;
    }

    for( int a=0; a< numProcs; a++ ){

        bzero( (void*)filename, 255 );
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%sgeombc.dat.%d",dirname2,a+1 );
	}
	else
        sprintf( filename, "%sgeombc.dat.%d",dirname,a+1 );
        openfile( filename, "append", &fgeom );
        
        isize = ncorp[a].size() ;
        nitems = 1;
        iarray[ 0 ] = ncorp[a].size() ;
        writeheader( &fgeom, " mode number map from partition to global",
                      (void*)iarray, &nitems, &isize, "integer", oformat );


        int* fncorp = new int [  ncorp[a].size() ];
        int count = 0;
        for( map< int, int>::iterator iter = ncorp[a].begin();
             iter != ncorp[a].end();
             iter++ )
            fncorp[count++] = (*iter).second ;

        nitems = count;
        writedatablock( &fgeom, " mode number map from partition to global",
                         (void*)(fncorp), &nitems, "integer", oformat );

        delete [] fncorp ;
        closefile( &fgeom, "append" );
        ncorp[a].clear();
    } 

    //write restart

    sprintf( filename, "restart.%d.1", stepno );  
    openfile( filename, "read", &irestart );
    readheader( &irestart, "solution?", (void*)iarray, &ithree, 
                 "double", iformat );

    nshg  = iarray[0];
    ndof  = iarray[1];
    isize = nshg * ndof ;

    double* solution = new double [ isize ];

    readdatablock( &irestart,"solution?", (void*)solution, &isize, 
                    "double", iformat  );

    //time derivative of solution
    int accelexist;
    double* acceleration;
    iarray[0] = -1;
    readheader(&irestart,"time derivative of solution?",(void*)iarray,&ithree,
                "double",iformat);
    if(iarray[0]==-1) accelexist = 0;
    else{
        accelexist = 1;
        nshg  = iarray[0];
        ndof  = iarray[1];
        isize = nshg * ndof ;

        acceleration = new double [ isize ];

        readdatablock( &irestart,"time derivative of solution?", (void*)acceleration, &isize, 
                        "double", iformat  ); 
    }
       
    //displacement
    int dispexist, ndisp;
    double* displacement;
    iarray[0] = -1;
    readheader( &irestart, "displacement?", (void*)iarray, &ithree, 
                 "double", iformat );
    if(iarray[0] == -1) dispexist = 0; //displacement doesnot exist
    else{
        dispexist = 1;
        nshg  = iarray[0];
        ndisp  = iarray[1];
        isize = nshg * ndisp ;

        displacement = new double [ isize ];

        readdatablock( &irestart,"displacement?", (void*)displacement, &isize, 
                        "double", iformat  ); 
    }

    // now we need to create the partitioned data structure.
    // keeping in mind that unlike BCs this has a direct corelation 
    // to each nshg so we need to use a 
    // sorted data structure and we choose a std::map for simplicity

    vector< map< int, vector< double > > > solPart( numProcs );

    for( int x = 1; x < nshg+1; x++ ){
        for( map<int,int>::const_iterator pIter = ParallelData[x].begin();
             pIter != ParallelData[x].end();
             pIter++ ) {
            for( int v=0; v< ndof ; v++ )
                solPart[ (*pIter).first ][ (*pIter).second ].push_back( 
                    solution[ v*nshg + x -1 ] );
        }
    }
      
    delete [] solution;   

   // partition time derivative of solution
    vector< map< int, vector< double > > > accelPart( numProcs );

    if(accelexist){
        for( int x = 1; x < nshg+1; x++ ){
            for( map<int,int>::const_iterator pIter = ParallelData[x].begin();
                 pIter != ParallelData[x].end();
                 pIter++ ) {
                for( int v=0; v< ndof ; v++ )
                    accelPart[ (*pIter).first ][ (*pIter).second ].push_back( 
                        acceleration[ v*nshg + x -1 ] );
            }
        }
        
        delete [] acceleration;
    }

    //partition displacement
    vector< map< int, vector< double > > > dispPart( numProcs );

    if(dispexist){
        for( int x = 1; x < nshg+1; x++ ){
            for( map<int,int>::const_iterator pIter = ParallelData[x].begin();
                 pIter != ParallelData[x].end();
                 pIter++ ) {
                for( int v=0; v< ndisp ; v++ )
                    dispPart[ (*pIter).first ][ (*pIter).second ].push_back( 
                        displacement[ v*nshg + x -1 ] );
            }
        }
        
        delete [] displacement;
    }

    // now we have the partitioned data structure we need to write it out.

    for( int a=0; a< numProcs; a++ ){ 
        
        int nshgLocal = solPart[a].size();
        double* fSolution = new double [ nshgLocal * ndof ];
        double* fAcceleration;
        double* fDisplacement;
        if(acceleration)fAcceleration = new double[nshgLocal * ndof];
        if(dispexist) fDisplacement = new double[nshgLocal * ndisp];
#if defined ( DEBUG )
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf( filename,"%srestart.asc.%d",dirname2, a+1 );
	}
	else
        sprintf( filename,"%srestart.asc.%d",dirname, a+1 );
        FILE* rascii = fopen( filename,"w");
        fprintf( rascii, "nshgLocal: %d\n", nshgLocal);
        fprintf( rascii, "numVars: %d\n", ndof );
        fprintf( rascii, "Step Number: %d\n", stepno);
        fprintf( rascii, "---------------------------------------------\n");
        for( int y=1; y < nshgLocal+1 ; y++ ) {
            fprintf( rascii, "%d :", y ) ;
            for( int w=0; w< ndof; w++ ) { 
                fprintf(rascii," %f ", solPart[a][y][w] );
            }
            fprintf( rascii, "\n");
        }
        fclose( rascii );
#endif
        for( int w=0; w< ndof; w++ ) 
            for( int y=1; y < nshgLocal+1 ; y++ ) {
                fSolution[ w*nshgLocal+(y-1) ] = solPart[a][y][w];
            }

        solPart[a].clear();

        if(accelexist){
            for( int w=0; w< ndof; w++ ) 
                for( int y=1; y < nshgLocal+1 ; y++ ) {
                    fAcceleration[ w*nshgLocal+(y-1) ] = accelPart[a][y][w];
                }
            
            accelPart[a].clear();        
        }
 
        if(dispexist){
            for( int w=0; w< ndisp; w++ ) 
                for( int y=1; y < nshgLocal+1 ; y++ ) {
                    fDisplacement[ w*nshgLocal+(y-1) ] = dispPart[a][y][w];
                }
            
            dispPart[a].clear();        
        }
        
        bzero( (void*)filename, 255 );
    if (numProcs > idirtrigger) {
        sprintf(dirname2,"%d-procs_case/%d-set/",numProcs,((int)(a/idirstep))*idirstep);   // Igor, August 2012
        sprintf(filename,"%srestart.%d.%d",dirname2, stepno, a+1 );
	}
	else
        sprintf(filename,"%srestart.%d.%d",dirname, stepno, a+1 );
        openfile( filename, "write", &frest );
        
        for( vector< string >::iterator iter = headers.begin();
             iter != headers.end();
             iter++ ) {
            writestring( &frest, (*iter).c_str() );
            writestring( &frest, "\n" );
        }
        headers.clear();

        bzero( (void*)filename, 255 );
        sprintf(filename,"number of modes : < 0 > %d\n", nshgLocal);
        writestring( &frest, filename );
        
        bzero( (void*)filename, 255 );
        sprintf(filename,"number of variables : < 0 > %d\n", ndof);
        writestring( &frest, filename );
        
        isize = 1;
        nitems = 1;
        iarray[ 0 ] = 1;
        writeheader( &frest, "byteorder magic number ", 
                      (void*)iarray, &nitems, &isize, "integer", oformat );
        
        nitems = 1;
        writedatablock( &frest, "byteorder magic number ",
                         (void*)mptr, &nitems, "integer", oformat );
        
        
        isize = nshgLocal*ndof;
        nitems = 3;
        iarray[ 0 ] = nshgLocal;
        iarray[ 1 ] = ndof;
        iarray[ 2 ] = stepno;
        writeheader( &frest, "solution ", 
                      (void*)iarray, &nitems, &isize, "double", oformat );
        
        
        nitems = nshgLocal*ndof;
        writedatablock( &frest, "solution ",
                         (void*)(fSolution), &nitems, "double", oformat );
        
        
        delete [] fSolution;

        if(accelexist) { //has time derivative of solution
            isize = nshgLocal*ndof;
            nitems = 3;
            iarray[ 0 ] = nshgLocal;
            iarray[ 1 ] = ndof;
            iarray[ 2 ] = stepno;
            writeheader( &frest, "time derivative of solution ", 
                          (void*)iarray, &nitems, &isize, "double", oformat );
        
        
            nitems = nshgLocal*ndof;
            writedatablock( &frest, "time derivative of solution ",
                             (void*)(fAcceleration), &nitems, "double", oformat );
        
            delete [] fAcceleration;
        }

        if(dispexist) { //has displacement
            isize = nshgLocal*ndisp;
            nitems = 3;
            iarray[ 0 ] = nshgLocal;
            iarray[ 1 ] = ndisp;
            iarray[ 2 ] = stepno;
            writeheader( &frest, "displacement ", 
                          (void*)iarray, &nitems, &isize, "double", oformat );
        
        
            nitems = nshgLocal*ndisp;
            writedatablock( &frest, "displacement ",
                             (void*)(fDisplacement), &nitems, "double", oformat );
        
            delete [] fDisplacement;
        }

        closefile( &frest, "write" );

    }        
    sprintf(filename,"%snumstart.dat", dirname );
    ofstream nstart( filename );
    nstart << stepno << endl;
    nstart.close();
    closefile( &igeombc, "read" );
    closefile( &irestart, "read" );
    ParallelData.clear();
    return 0;
}
