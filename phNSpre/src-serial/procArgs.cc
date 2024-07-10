#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


using namespace std;

extern int globalP,rStart,intBC,prCd,
           zScale,ensa_dof=0,rRead=5,
           NSFTAG=-1, lStart=0, old_format;

extern bool FaceConnectivity;
extern bool isReorder;

bool refine_mesh=false, c1quintic=false;
char options[255];
char* iformat = "binary";
char* oformat = "binary";
#ifdef VERSION
  char version[]=VERSION;
#elif
  char version[]="UNKNOWN";
#endif

void 
procArgs( int argc, char *argv[], char fname[], char mname[], char gname[] ) {
    int c;
   
    char* executable = argv[0];
    /* first copy the entire argv into out dumpstring */
    if ( argc > 1 ){
        strcpy(options,argv[1]);
        for(c=2; c< argc; c++){
            strcat(options," ");
            strcat(options,argv[c]);
        }
    }

    /* look for command line arguments */
    while (--argc > 0 && (*++argv)[0] == '-')  {
        while ((argc > 0) && (c = *++argv[0]))  {
            switch (c)  {
            case 'c':                        /* input redirection*/
                if(!freopen( argv[1],"r",stdin)){
                    cerr << "could not redirect input" << endl ;
                    exit(1);
                }
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'o':                        /* output format*/
                oformat = argv[1];
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
             case 'i':                        /* input format*/
                iformat = argv[1];
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'p':                        /* polynomial order */
                globalP = atoi(argv[1]);
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'T':                   
                NSFTAG = atoi(argv[1]);
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'v':
                ensa_dof=atoi(argv[1]);
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                break;
            case 'f':
                strcpy(fname,argv[1]);  /* attribute file name */
                argc--;
                argc--;
                (*++argv); 
                (*++argv); 
                break;
            case 'm':
                strcpy(mname,argv[1]);  /* Mesh file name */
                argc--;
                argc--;
                (*++argv); 
                (*++argv); 
                break;
            case 'G':
                strcpy(gname,argv[1]);  /* Model file name */
                argc--;
                argc--;
                (*++argv); 
                (*++argv); 
                break;

            case 'h':
                cout <<"\nUsage: " << executable 
                     <<" [-r <numdof to read>][-i ascii][ -o ascii ][-c stdinfile]" 
                     <<" [-tVTsIlLRb][-p <order>][-f <filename>] [-m <mesh file name>] \n" 
                     <<"\t-p: sets the global polynomial order.\n"
                     <<"\t-f: sets the attribute file name. Default is 'geom.atdb'.\n"
                     <<"\t-m: sets the mesh file name. Default is 'geom.sms'.\n"
                     <<"\t-G: sets the model file name. Default is geom.xmt_txt (PARASOL) - geom.sdm (DISCR)'.\n"
                     <<"\t-o: write output in ascii.\n"
                     <<"\t-i: read restart in ascii.\n"
                     <<"\t-[r numVars]: use file 'restartc.inp' for initial condition  .\n"
                     <<"\t-b: print boundary condition codes for geometric model.\n"
                     <<"\t-g: write graph.out, the partition graph\n"
                     <<"\t-s: prompt user for a z-coordinate scaling factor\n"
		     <<"\t-B: prints out old format of geombc.dat.1. Use it only for\n"
		     <<"\t    1 processor case (for more proc. use -B in parallel part.\n"   
                     <<"\t-F: get the face connectivity\n"
                     <<"\t-R: use file 'restar.inp' for initial condition and\n"
                     <<"\t    convert first variable from density to pressure.\n"
                     <<"\t-l: use linear 'restartc.inp' for initial condition\n"
                     <<"\t    for p > 1.\n"
                     <<"\t-t: refine based on error.\n"
                     <<"\t-C: c1 quintic triangular 2D mesh.\n"
                     <<"\t-N: No DOF and region reOrdering by RCM algorithm.\n"
                     <<"\n"
                     <<"\n\n";
          
                exit(1);
                break;
            case 't':
                refine_mesh = true;
                break;
            case 'C':			/* c1 quintic triagular mesh */
                c1quintic = true;
                break;
            case 's':
                zScale = 1;
                break;
            case 'I':                        /* lin-2-quad*/
                lStart=1;
                break;
            case 'r':                        /* use restart file */
                rRead = atoi(argv[1]);
                argc--;
                argc--;
                (*++argv);
                (*++argv);
                rStart = 1;
                break;
            case 'R':
                rStart = 2; /* use restart file and convert from density
                               to pressure */
                break;
            case 'l':
                rStart = 3; /* use linear portion of restart file */
                break;
            case 'L':
                rStart = 4; /* use quadratic portion of restart file */
                break;
            case 'b':
                prCd = 1;
                break;
	    case 'B':
	        old_format = 1;
		break;
            case 'F':
                /* get the face connectivity */
                FaceConnectivity = 1;
                break;
            case 'N': /* dof and mesh region reorder by RCM algorithm */
                isReorder = 0;
                break;
            case 'V':
                cout << "You are using NSpre Version : "<< version << endl;
		        exit(0);
            default:
                cout <<" Unknown option  "<< c 
                     << " Please use NSpre -h for the list of available options\n";
            }
        }
    }
}
