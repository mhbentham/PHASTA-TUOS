/* This file provides interface functions for 'partial ' random 
   access into the PHASTA input files 

   Anil Karanam March 2001 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "phastaIO.h"
#include <FCMangle.h>
#include "new_interface.h"

#include "common_c.h"

#ifdef intel
#include <winsock2.h>
#else
#include <unistd.h>
#include <strings.h>
#endif

void igetMinMaxAvg(int *ivalue, double *stats, int *statRanks) {
  int isThisRank;

  double *value = (double*)malloc(sizeof(double));
  *value = 1.0*(*ivalue);

  rgetMinMaxAvg(value,stats,statRanks);

  /* MPI_Allreduce(value,&stats[0],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[0])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[0],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[1],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[1])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[1],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[2],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  stats[2] /= workfc.numpe; */

  free(value);
}

void rgetMinMaxAvg(double *value, double *stats, int *statRanks) {
  int isThisRank;

  MPI_Allreduce(value,&stats[0],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[0])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[0],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[1],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[1])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[1],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[2],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  stats[2] /= workfc.numpe;

  double sqValue = (*value)*(*value), sqValueAvg = 0.;
  MPI_Allreduce(&sqValue,&sqValueAvg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  sqValueAvg /= workfc.numpe;
  // stats[3] = sqValueAvg;

  stats[3] = sqrt(sqValueAvg-stats[2]*stats[2]);
}

void print_mesh_stats(void) {
  int statRanks[2];
  double iStats[4], rStats[4];

  igetMinMaxAvg(&conpar.nshg,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("nshg    : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&conpar.numel,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("numel   : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&conpar.numelb,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("numelb  : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&conpar.nnz_tot,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("nnz_tot : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
}

void print_mpi_stats(void) {
  int statRanks[2];
  double iStats[4], rStats[4];

  igetMinMaxAvg(&mpistats.iISend,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iISend : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iIRecv,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iIRecv : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iWaitAll,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iWtAll : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iAllR,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iAllR  : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);

  rgetMinMaxAvg(&mpistats.rISend,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rISend : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rIRecv,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rIRecv : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rWaitAll,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rWtAll : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rCommu,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rCommu : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rAllR,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rAllR  : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
}

//void print_system_stats(double tcorecp[2]) {
void print_system_stats(double *tcorecp) {
  int statRanks[2];
  double iStats[4], rStats[4];
  double syst_assembly, syst_solve;

  syst_assembly = tcorecp[0];
  syst_solve = tcorecp[1];

  rgetMinMaxAvg(&syst_assembly,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("Elm. form. : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);

  rgetMinMaxAvg(&syst_solve,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("Lin. alg. sol : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);

  //printf("rank %d - syst_assembly %f - syst_solve %f\n",workfc.myrank,syst_assembly,syst_solve);
}
//MR CHANGE END



void countfieldstowriterestart()
{
  int nfields;

//     printf("TEST: %d %d %d %d %d\n",timdat.istep,timdat.itseq,inpdat.nstep[0],inpdat.nstep[1],timdat.lstep);

  nfields = 3; //time stamp, solution, time derivatives

  if(outpar.ivort == 1){
    nfields++; //vorticity
  }

  if(abs(turbvar.itwmod) != 1 && outpar.iowflux == 1) { 
    nfields++; //instantaneous wss in bflux.f
  }

//   if(ideformwall.eq.1) not handled yet

  if(timdat.istep == inpdat.nstep[timdat.itseq-1]){ //Last time step of the computation

    //projection vectors and pressure projection vectors (call saveLesRestart in itrdrv)
    nfields = nfields +2;

    //if Print Error Indicators = true (call write_error in itrdrv)
    if(turbvar.ierrcalc == 1){
      nfields++;
    }

    //if Print ybar = True (call write_field(myrank,'a','ybar',4,... in itrdrv)
    if(outpar.ioybar == 1){
      nfields++;  //ybar

      //phase average fields
      //if(outpar.nphasesincycle >0) {
      //  nfields = nfields + outpar.nphasesincycle;
      //}

      if(abs(turbvar.itwmod) != 1 && outpar.iowflux == 1) { 
        nfields++; //wssbar
      }

    }

//    if(turbvari.irans < 0) {
//    printf("irans: %d - idns: %d - itwmod: %d - id2w: %d\n",turbvari.irans,turbvari.idns,turbvar.itwmod,levlset.id2w);
    if (turbvari.irans < 0 || (turbvari.idns!=0 && abs(turbvar.itwmod)==1) || levlset.id2w==1) {
      nfields++; //dwal
    }

  }

  outpar.nsynciofieldswriterestart = nfields;

  if(workfc.myrank == 0) {
    printf("Number of fields to write in restart files: %d\n", nfields);
  }
}


void
Write_Restart(  int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
		        double* soltime,
                double* array1,
                double* array2 ) {

    char fname[255];
    char rfile[60];
    char existingfile[30], linkfile[30];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];

//MR CHANGE
    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    // First, count the number of fields to write and store the result in
    countfieldstowriterestart();

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    nfields = outpar.nsynciofieldswriterestart;
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
    nprocs = workfc.numpe;
//MR CHANGE END
    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    int descriptor;
    char filename[255],path[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    sprintf(filename,"restart-dat.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");

    if (*pid==0) {
      printf("Filename is %s \n",filename);
    }

    openfile(filename, "write", &f_descriptor);

    field_flag=0;

     int i;
     for ( i = 0; i < nppp; i++) { 
     // This loop is useful only if several parts per processor during computation
     // In practice, nppp = 1
     // Keep only this loop as an example in case of future implementation change
     // GPID : global part id, corresponds to rank ...
        // e.g : (in this example)
        // proc 0 : 1--4
        // proc 1 : 5--8 ...
        GPID = startpart + i;

        sprintf(fieldtag_s,"TimeStamp@%d",GPID);
        isize = 1;
        nitems = 1;
        iarray[ 0 ] = 1;

        writeheader( &f_descriptor, fieldtag_s, 
                      (void*)iarray, &nitems, &isize,"double", phasta_iotype);

        writedatablock( &f_descriptor, fieldtag_s, 
                        (void*)soltime, &isize, "double", phasta_iotype );
     }
     field_flag++;


    // Write solution field ...
    sprintf(fieldtag_s,"solution@%d",startpart);

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
    writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );
    field_flag++;


    // Write time derivative of solution field ...
    sprintf(fieldtag_s,"time derivative of solution@%d",startpart);

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
    writedatablock( &f_descriptor, fieldtag_s, (void*)(array2), &isize, "double", phasta_iotype );
    field_flag++;


    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("\n");
      }
    }

}


void
Write_Error(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 ) {


    char fieldtag[255];
    int isize, nitems;
    int iarray[10];

    int nfields;
    int numparts;
    int irank;
    int startpart;

    nfields = outpar.nsynciofieldswriterestart;
    irank = *pid; //workfc.myrank;

    // Calculate number of parts each  proc deal with and where it start and end ...
    startpart = irank +1;// Part id from which I (myrank) start ...

    field_flag++;

    sprintf(fieldtag,"errors@%d",startpart);

    if(*pid==0) {
      printf("\n");
      printf("The %d/%d th field to be written is '%s'\n",field_flag,nfields,fieldtag);
    }

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &f_descriptor, fieldtag, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
    writedatablock( &f_descriptor, fieldtag, (void*)array1, &isize, "double", phasta_iotype );

    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
      }
    }

}


void
Write_Displ(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 ) {
// Needs to be updated for SyncIO

    char fname[255];
    char rfile[60];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    time_t timenow = time ( &timenow);
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    sprintf(rfile,"restart.%d.%d",*stepno,*pid+1);
    openfile(rfile,"append", &irstou);

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &irstou, "displacement", (void*)iarray, &nitems, &isize, "double", phasta_iotype );

    nitems = (*nshg)*(*numVars);
    writedatablock( &irstou, "displacement", (void*)(array1), &nitems, "double", phasta_iotype );

    closefile( &irstou, "append" );

}

void
Write_Field(  int *pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno) {

    char *fieldlabel = (char *)malloc((*tagsize+1)*sizeof(char));
    strncpy(fieldlabel, fieldtag, *tagsize);
    fieldlabel[*tagsize] = '\0';

    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    char fmode[10];
    if(!strncmp(filemode,"w",1))
      strcpy(fmode,"write");
    else // default is append
      strcpy(fmode,"append");

    char datatype[10];
    if(!strncmp(arraytype,"i",1))
      strcpy(datatype,"int");
    else // default is double
      strcpy(datatype,"double");

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    nfiles = outpar.nsynciofiles;
    nfields = outpar.nsynciofieldswriterestart;
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
    nprocs = workfc.numpe;

    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each  proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    char filename[255],path[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    strncpy(fieldlabel, fieldtag, *tagsize);

    field_flag++;
    if(*pid==0) {
      printf("\n");
      printf("The %d/%d th field to be written is '%s'\n",field_flag,nfields,fieldlabel);
    }

    sprintf(filename,"restart-dat.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));

    sprintf(fieldtag_s,"%s@%d",fieldlabel,startpart);

    isize = (*nshg)*(*numvars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numvars);
    iarray[ 2 ] = (*stepno);
    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, datatype, phasta_iotype);
    writedatablock( &f_descriptor, fieldtag_s, array, &isize, datatype, phasta_iotype );


    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
      }
    }

    free(fieldlabel);
}


void
Write_d2wall(   int* pid,
                int* numnp,
                int* numVars,
                double* array1 ) {


    int isize, nitems;
    int iarray[10];


    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    //  First, count the number of fields to write and store the result in
    //countfieldstowriterestart();

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    nfields = 1; //outpar.nsynciofieldswriterestart;  // Only the distance to the walls in d2wall
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
    nprocs = workfc.numpe;
    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each proc deal with and where it start and end ...
    int startpart = irank+1;// Part id from which I (myrank) start ...

    int descriptor;
    char filename[255],path[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    sprintf(filename,"d2wall.%d",((int)(irank/(nprocs/nfiles))+1));

    if (irank==0) {
      printf("Filename is %s \n",filename);
    }

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");

    openfile(filename, "write", &f_descriptor);

    field_flag=0;

    // Write solution field ...
    sprintf(fieldtag_s,"d2wall@%d",startpart);

    isize = (*numnp)*(*numVars);
    nitems = 2;
    iarray[ 0 ] = (*numnp);
    iarray[ 1 ] = (*numVars);

    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);

    writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );

    field_flag++;

    if (field_flag==nfields){

      closefile(&f_descriptor, "write");

      finalizephmpiio(&f_descriptor);

      if (irank==0) {
        printf("\n");
      }
    }
}

void
Write_GradPhi(  int* pid, 
                int* stepno, 
                int* nshg, 
                int* numVars,
		        double* x,
                double* y, 
                double* gradphi,
	            double* mag ) {

    printf("Write_GradPhi not updated to SyncIO yet\n");
    return;

    char fname[255];
    char rfile[60];
    FILE* file=NULL ;
    int j;
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    time_t timenow = time ( &timenow);
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    sprintf(rfile,"gradphi.%d.%d",*stepno,*pid+1);
    openfile(rfile,"write", &irstou);

    /* writing the top ascii header for the gradphi file */

    writestring( &irstou,"# PHASTA GradPhi File Version 1.0\n");
    writestring( &irstou,
                  "# format \"keyphrase : sizeofnextblock usual headers\"\n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# Output generated by phasta version (NOT YET CURRENT): %lf \n", version);
    writestring( &irstou, fname );

    bzero( (void*)fname, 255 );
    gethostname(fname,255);
    writestring( &irstou,"# This result was produced on: ");
    writestring( &irstou, fname );
    writestring( &irstou,"\n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# %s\n", ctime( &timenow ));
    writestring( &irstou, fname );
        
    isize = 1;
    nitems = 1;
    iarray[ 0 ] = 1;
    writeheader( &irstou, "byteorder magic number ", 
                  (void*)iarray, &nitems, &isize, "integer", phasta_iotype );
    
    nitems = 1;
    writedatablock( &irstou, "byteorder magic number ",
                     (void*)mptr, &nitems, "integer", phasta_iotype );
    
    
    bzero( (void*)fname, 255 );
    sprintf(fname,"number of modes : < 0 > %d\n", *nshg);
    writestring( &irstou, fname );
    
    bzero( (void*)fname, 255 );
    sprintf(fname,"number of variables : < 0 > %d\n", *numVars);
    writestring( &irstou, fname );
    
    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &irstou, "gradphi ", 
                  (void*)iarray, &nitems, &isize, "double", phasta_iotype );
    
        
    nitems = (*nshg)*(*numVars);
    writedatablock( &irstou, "gradphi ",
                     (void*)(gradphi), &nitems, "double", phasta_iotype );
        

    closefile( &irstou, "write" );
}

void
Write_PrimVert( int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
                double* primvertval ) {


    printf("Write_PrimVert not updated to SyncIO yet\n");
    return;

    char fname[255];
    char rfile[60];
    FILE* file=NULL ;
    int j;
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    time_t timenow = time ( &timenow);
    double version=0.0;
    int isize, nitems;
    int iarray[10];

// Binary version
//

    sprintf(rfile,"primvert.%d.%d",*stepno,*pid+1);
    openfile(rfile,"write", &irstou);

    /* writing the top ascii header for the gradphi file */

    writestring( &irstou,"# PHASTA Primary Vertex File Version 1.0\n");
    writestring( &irstou,
                  "# format \"keyphrase : sizeofnextblock usual headers\"\n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# Output generated by phasta version (NOT YET CURRENT): %lf \n", version);
    writestring( &irstou, fname );

    bzero( (void*)fname, 255 );
    gethostname(fname,255);
    writestring( &irstou,"# This result was produced on: ");
    writestring( &irstou, fname );
    writestring( &irstou,"\n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# %s\n", ctime( &timenow ));
    writestring( &irstou, fname );

    isize = 1;
    nitems = 1;
    iarray[ 0 ] = 1;
    writeheader( &irstou, "byteorder magic number ",
                  (void*)iarray, &nitems, &isize, "integer", phasta_iotype );

    nitems = 1;
    writedatablock( &irstou, "byteorder magic number ",
                     (void*)mptr, &nitems, "integer", phasta_iotype );


    bzero( (void*)fname, 255 );
    sprintf(fname,"number of modes : < 0 > %d\n", *nshg);
    writestring( &irstou, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of variables : < 0 > %d\n", *numVars);
    writestring( &irstou, fname );

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &irstou, "primvert",
                  (void*)iarray, &nitems, &isize, "double", phasta_iotype );


    nitems = (*nshg)*(*numVars);
    writedatablock( &irstou, "primvert",
                     (void*)(primvertval), &nitems, "double", phasta_iotype );


    closefile( &irstou, "write" );
}

