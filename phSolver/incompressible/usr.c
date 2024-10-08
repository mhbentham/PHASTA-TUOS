/*===========================================================================
 *
 * "usr.c":  user's function
 *
 *===========================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include "les.h"
#include "usr.h"
#include "common_c.h"
#include "phastaIO.h"
#include "new_interface.h"
#include <FCMangle.h>

extern char phasta_iotype[80];
extern int field_flag; //new_interface.c //TODO: figure out where these
extern int f_descriptor; //new_interface.c // should be

/*===========================================================================
 *
 * "usrNew":  Put all the values in usrHd
 *
 * From FORTRAN
 *
 *	integer		usr(100)
 *	dimension	aperm(numnp,nperm)
 *	...
 *	call usrnew( usr, aperm, ..., numnp, ...)
 *
 *
 *===========================================================================
 */
#include "mpi.h"
static int lmNum = 0;
static LesHd lesArray[8];
void   usrNew(	UsrHd	  usrHd,
                        int*      eqnType,
                        double*	  aperm,
                        double*	  atemp,
                        double*   resf,
                        double*   solinc,
                        double*   flowDiag,
                        double*   sclrDiag,
                        double*   lesP,
                        double*   lesQ,
                        int*      iBC,
                        double*   BC,
                        int*      iper,
                        int*      ilwork,
                        int*      numpe,
                        int*      nNodes,
                        int*      nenl,
                        int*	  nPermDims,
                        int*	  nTmpDims,
                        int*	  rowp,
                        int*	  colm,
                        double*   lhsK,
                        double*   lhsP,
                        double*   lhsS,
                        int*      nnz_tot,
                        double*   CGsol
    )
{
    char*	funcName = "usrNew" ;	/* function name		*/

/*---------------------------------------------------------------------------
 * Stick the parameters
 *---------------------------------------------------------------------------
 */
    usrHd->eqnType      = *eqnType ;
    usrHd->aperm	= aperm ;
    usrHd->atemp	= atemp ;
    usrHd->resf         = resf ;
    usrHd->solinc       = solinc ;
    usrHd->flowDiag     = flowDiag ;
    usrHd->sclrDiag     = sclrDiag ;
    usrHd->lesP         = lesP ;
    usrHd->lesQ         = lesQ ;
    usrHd->iBC          = iBC  ;
    usrHd->BC           = BC   ;
    usrHd->iper         = iper ;
    usrHd->ilwork       = ilwork ;
    usrHd->numpe        = *numpe ;
    usrHd->nNodes	= *nNodes ;
    usrHd->nenl         = *nenl ;
    usrHd->nPermDims	= *nPermDims ;
    usrHd->nTmpDims	= *nTmpDims ;
    usrHd->rowp	        = rowp ;
    usrHd->colm	        = colm ;
    usrHd->lhsK	        = lhsK ;
    usrHd->lhsP	        = lhsP ;
    usrHd->lhsS         = lhsS ;
    usrHd->nnz_tot      = nnz_tot ;
    usrHd->CGsol        = CGsol;
} /* end of usrNew() */

/*===========================================================================
 *
 * "usrPointer":  Get the pointer
 *
 *===========================================================================
 */
Real*
usrPointer(	UsrHd	usrHd,
            Integer	id,
            Integer	offset,
            Integer	nDims )
{
    char*	funcName = "usrPointer";/* function name		*/
    Real*	pnt ;			/* pointer			*/

/*---------------------------------------------------------------------------
 * Get the head of the memory
 *---------------------------------------------------------------------------
 */
    if ( id == LES_RES_PNT ) {

        pnt	= usrHd->resf ;
        id	= 0 ;

    } else if ( id == LES_SOL_PNT ) {

        pnt	= usrHd->solinc ;
        id	= 0 ;

    } else if ( id < 0 ) {

        pnt	= usrHd->aperm ;
        id	= id + usrHd->nPermDims ;

    } else {

        pnt	= usrHd->atemp ;
        id	= id ;

    }
/*---------------------------------------------------------------------------
 * Get the offset
 *---------------------------------------------------------------------------
 */
    pnt		= pnt + (id + offset) * usrHd->nNodes ;

/*---------------------------------------------------------------------------
 * Return the pointer
 *---------------------------------------------------------------------------
 */
    return( pnt ) ;

} /* end of usrPointer() */

#define myflesnew FortranCInterface_GLOBAL_(myflesnew,MYFLESNEW)
#define myflessolve FortranCInterface_GLOBAL_(myflessolve,MYFLESSOLVE)
#define savelesrestart FortranCInterface_GLOBAL_(savelesrestart,SAVELESRESTART)
#define readlesrestart FortranCInterface_GLOBAL_(readlesrestart,READLESRESTART)
#define solverlicenseserver FortranCInterface_GLOBAL_(solverlicenseserver,SOLVERLICENSESERVER)



#ifdef intel
        lesArray[ *lesId ] = lesNew( fileName, *lmport, &lmNum, *eqnType,
                                     *nDofs, *minIters, *maxIters, *nKvecs,
                                     *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,presPrec,
                                     *tol, *presTol, *verbose, stats, nPermDims,
                                     nTmpDims );
    return ;}
/* the following is a fake function that was required when we moved to
   a C++ main on in the MS Visual Studio environment.  It fails to
   link because it is looking for this function
*/
void  _CrtDbgReport() {
    return ;}

double __vcos_(double fg) { fflush(stdout); printf(" vcos got called \n"); fflush(stdout);}
double __vlog_(double fg)  { fflush(stdout); printf(" vlog got called \n"); fflush(stdout);}


#endif /* we are in unix land... whew.  secretly we have equivalenced fileName and  */

/* #ifdef LINUX*/
/* void flush_(int* junk ){ return; }*/
/* #endif*/
void    myflesnew(	     Integer*	lesId,
                         Integer*	lmport,
                         Integer*	eqnType,
                         Integer*	nDofs,
                         Integer*	minIters,
                         Integer*	maxIters,
                         Integer*	nKvecs,
                         Integer*	prjFlag,
                         Integer*	nPrjs,
                         Integer*	presPrjFlag,
                         Integer*	nPresPrjs,
                         Real*	    tol,
                         Real*     	presTol,
                         Integer*	verbose,
                         Real*     	stats,
                         Integer*	nPermDims,
                         Integer*	nTmpDims,
                         char*      lmhost          ) {
    int procId;
#ifdef AMG
    int presPrec=1;
#else
    int presPrec=0;
#endif
    MPI_Comm_rank( MPI_COMM_WORLD, &procId ) ;
//    printf("\n USR.C:  line 232 \n");    
    if(lmNum==0){
        if(procId==0){
//    printf("\n USR.C:  line 235 \n");
            lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                         *nDofs, *minIters, *maxIters, *nKvecs,
                                         *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,presPrec,
                                         *tol, *presTol, *verbose, stats, nPermDims,
                                         nTmpDims );
            MPI_Bcast( &lmNum, 1, MPI_INT, 0, MPI_COMM_WORLD ) ;
        } else {
//    printf("\n USR.C:  line 243 \n");
            MPI_Bcast( &lmNum, 1, MPI_INT, 0, MPI_COMM_WORLD ) ;
            lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                         *nDofs, *minIters, *maxIters, *nKvecs,
                                         *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,presPrec,
                                         *tol, *presTol, *verbose, stats, nPermDims,
                                         nTmpDims );
        }
    } else {
//    printf("\n USR.C:  line 252 \n");
        lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                     *nDofs, *minIters, *maxIters, *nKvecs,
                                     *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,presPrec,
                                     *tol, *presTol, *verbose, stats, nPermDims,
                                     nTmpDims );
    }
//    printf("\n USR.C:  line 259 \n");
    return ;
}


void
savelesrestart( Integer* lesId,
                 Real*    aperm,
                 Integer* nshg,
                 Integer* myrank,
                 Integer* lstep,
                 Integer* nPermDims ) {

    int nPrjs, PrjSrcId;
    int nPresPrjs, PresPrjSrcId;
    char filename[255];
    int fileHandle=0;
    int iarray[3];
    int size, nitems;
    double* projVec;
    int i, j, count;

    nPrjs = (Integer) lesGetPar( lesArray[ *lesId ], LES_ACT_PRJS );
    PrjSrcId = (Integer) lesGetPar( lesArray[ *lesId ], LES_PRJ_VEC_ID );

    if ( PrjSrcId < 0 ) PrjSrcId += *nPermDims;

    projVec = (double*)malloc( nPrjs * ( *nshg ) * sizeof( double ) );

    count = 0;
    for( i = PrjSrcId; i < PrjSrcId+nPrjs; i ++ ) {
        for( j = 0 ; j < *nshg; j++ ) {
            projVec[ count++ ] = aperm[ (*nshg) * i + j ];
        }
    }

    iarray[ 0 ] = *nshg;
    iarray[ 1 ] = nPrjs;
    nitems = 2;
    size = (*nshg)*nPrjs;

    int name_length;
    name_length = 18;
    Write_Field(myrank,"a","projection vectors",&name_length, (void *)projVec,"d", nshg, &nPrjs, lstep);
    free(projVec);

    /*************************************************************************/

    nPresPrjs = (Integer) lesGetPar( lesArray[ *lesId ], LES_ACT_PRES_PRJS );
    PresPrjSrcId =(Integer)lesGetPar( lesArray[ *lesId ], LES_PRES_PRJ_VEC_ID );
    if ( PresPrjSrcId < 0 ) PresPrjSrcId += *nPermDims;

    projVec = (double*)malloc( nPresPrjs * ( *nshg ) * sizeof( double ) );

    count = 0;
    for( i = PresPrjSrcId; i < (PresPrjSrcId + nPresPrjs) ; i ++ ) {
        for( j = 0 ; j < *nshg; j++ ) {
            projVec[ count++ ] = aperm[ (*nshg) * i + j ];
        }
    }

    iarray[ 0 ] = *nshg;
    iarray[ 1 ] = nPresPrjs;
    nitems = 2;
    size = (*nshg)*nPresPrjs;

    name_length = 27;
    Write_Field(myrank,"a","pressure projection vectors",&name_length, projVec,"d", nshg, &nPresPrjs, lstep);

    free( projVec);
}

void
readlesrestart( Integer* lesId,
                 Real*    aperm,
                 Integer* nshg,
                 Integer* myrank,
                 Integer* lstep ,
                 Integer* nPermDims ) {

    int nPrjs, PrjSrcId;
    int nPresPrjs, PresPrjSrcId;
    char filename[255];
    int fileHandle=0;
    int iarray[3]={-1,-1,-1};
    int size, nitems;
    int itwo=2;
    int lnshg;
    double* projVec;
    int i,j,count;

//MR CHANGE

    int nfiles;
    int nfields;
    int numParts;
    int nprocs;
    int nppf;

    nfiles = outpar.nsynciofiles;
    numParts = workfc.numpe;
    nprocs = workfc.numpe;
//MR CHANGE END

    char fieldname[255];

    // Calculate number of parts each proc deal with and where it start and end ...
    int nppp = numParts/nprocs;        // nppp : Number of parts per proc ...
    int startpart = *myrank * nppp +1;    // Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;  // Part id to which I (myrank) end ...

    sprintf( filename,"restart-dat.%d.%d",*lstep,((int)(*myrank/(numParts/nfiles))+1));
    queryphmpiio(filename, &nfields, &nppf);
    initphmpiio(&nfields, &nppf, &nfiles,&fileHandle,"read");
//MR CHANGE END

    openfile( filename, "read", &fileHandle );



    if ( fileHandle < 0 ) return; // See phastaIO.cc for error fileHandle
    bzero((void*)fieldname,255);
    sprintf(fieldname,"projection vectors@%d",startpart);
    readheader( &fileHandle, fieldname, (void*)iarray,
                &itwo, "integer", phasta_iotype );
//MR CHANGE END

    if ( iarray[0] != *nshg ) {
        closefile( &fileHandle, "read" );
        if(workfc.myrank==workfc.master)
          printf("projection vectors are being initialized to zero (SAFE)\n");
        return;
    }

    lnshg = iarray[ 0 ] ;
    nPrjs = iarray[ 1 ] ;

    size = (*nshg)*nPrjs;

    if (nPrjs > 0) {
      projVec = (double*)malloc( size * sizeof( double ));

      readdatablock( &fileHandle, fieldname, (void*)projVec,
                    &size, "double", phasta_iotype );
    
      lesSetPar( lesArray[ *lesId ], LES_ACT_PRJS, (Real) nPrjs );
      PrjSrcId = (Integer) lesGetPar( lesArray[ *lesId ], LES_PRJ_VEC_ID );
      if ( PrjSrcId < 0 ) PrjSrcId += *nPermDims;

      count = 0;
      for( i = PrjSrcId; i < PrjSrcId+nPrjs; i ++ ) {
        for( j = 0 ; j < *nshg; j++ ) {
            aperm[ (*nshg) * i + j ] = projVec[ count++ ] ;
        }
      }

      free( projVec );
    }

    /************************************************************************/

    iarray[0] = -1; iarray[1] = -1; iarray[2] = -1;

//MR CHANGE
    bzero((void*)fieldname,255);
    sprintf(fieldname,"pressure projection vectors@%d",startpart);
//MR CHANGE END


    readheader( &fileHandle, fieldname, (void*)iarray,
                 &itwo, "integer", phasta_iotype );

    lnshg = iarray[ 0 ] ;
    nPresPrjs = iarray[ 1 ] ;

    if ( lnshg != *nshg )  {
        closefile( &fileHandle, "read" );
        if(workfc.myrank==workfc.master)
          printf("pressure projection vectors are being initialized to zero (SAFE)\n");
        return;
    }

    size = (*nshg)*nPresPrjs;
    if (nPresPrjs > 0) {
      projVec = (double*)malloc( size * sizeof( double ));
    
      readdatablock( &fileHandle, fieldname, (void*)projVec,
                    &size, "double", phasta_iotype );
    
      lesSetPar( lesArray[ *lesId ], LES_ACT_PRES_PRJS, (Real) nPresPrjs );
      PresPrjSrcId=(Integer)lesGetPar( lesArray[ *lesId ], LES_PRES_PRJ_VEC_ID );
      if ( PresPrjSrcId < 0 ) PresPrjSrcId += *nPermDims;
    
      count = 0;
      for( i = PresPrjSrcId; i < PresPrjSrcId+nPresPrjs; i ++ ) {
        for( j = 0 ; j < *nshg; j++ ) {
            aperm[ (*nshg) * i + j ] = projVec[ count++ ] ;
        }
      }

      free( projVec );
    }

    closefile( &fileHandle, "read" );

//MR CHANGE
    finalizephmpiio(&fileHandle);
//MR CHANGE END

}

void  myflessolve( Integer* lesId,
                    UsrHd    usrHd){
    lesSolve( lesArray[ *lesId ], usrHd );
}


int solverlicenseserver(char key[]){
#ifdef intel
    strcpy(key,"C:\\cygwin\\license.dat");
#else
    char* env_server_name;
    env_server_name = getenv("LES_LICENSE_SERVER");
    if(env_server_name) strcpy(key, env_server_name);
    else {
        if(workfc.myrank==workfc.master) {
          fprintf(stderr,"environment variable LES_LICENSE_SERVER not defined \n");
//          fprintf(stderr,"using wesley as default \n");
          fprintf(stderr,"using SCOREC license server as default \n");
        }
//        strcpy(key, "wesley.scorec.rpi.edu");
        strcpy(key, "acusim.license.scorec.rpi.edu");
    }
#endif
     if (workfc.myrank==workfc.master) { fprintf(stdout, "Acusim server is set.  \n");
     	fprintf(stdout, key); 
        fprintf(stdout, "\n");
       }
    return 1;
}
