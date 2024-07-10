/********************************************************************/
/* Writes dx format files:                                          */
/* solution.bin, acceleration.bin, bflux.bin, stats.bin, errors.bin */
/*                                                                  */
/* Elaine Bohr                                                      */
/* February 2004                                                    */
/********************************************************************/
#include <stdio.h>
#ifndef intel
#include <unistd.h>
#endif
#ifdef intel
#include <direct.h>
#define getcwd _getcwd
#endif
/* 
   The code calling this has glued together the solution and geometry
   files and has passed in the following:

   nv is the number of degrees of freedom per node
   numel is the total number of elements
   ien is the TRUE global connectivity  
       (index local element node number and element
        number to TRUE global node number)
   numnp is the total number of nodes                          
   x is the coordinates (3*numnp)
   q is the solution (nv*numnp)
   nen is the number of nodes in an element
   RequestedField is the flag for differentiating different fields:
   RequestedField = 0 solution
   RequestedField = 1 acceleration
   RequestedField = 2 boundary fluxes
   RequestedField = 3 turbulent statistics
   RequestedField = 4 errors
   RequestedField = 5 velbar
*/

void wrtc_(int nv,    int numel, int* ien, int lstep, int nsd,
           int numnp, float* x, float* q, int nen, int RequestedField)
{
  int i,j,count;
  char *eltype, fname1[100], fname2[100];			/* element type */
  sprintf(fname1,"header%d.dx",lstep);
  if (RequestedField == 0) sprintf(fname1,"header%d-sol.dx",lstep);
  if (RequestedField == 5) sprintf(fname1,"header%d-vel.dx",lstep);
  FILE *f0 = fopen(fname1,"w");
  FILE *f1 = fopen("connect.bin","wb");
  FILE *f2 = fopen("points.bin","wb");
  FILE *f3;
  
  if (RequestedField == 0) sprintf(fname2,"solution%d.bin",lstep);
  if (RequestedField == 1) sprintf(fname2,"acceleration%d.bin",lstep);
  if (RequestedField == 2) sprintf(fname2,"bflux%d.bin",lstep);
  if (RequestedField == 3) sprintf(fname2,"stats%d.bin",lstep);
  if (RequestedField == 4) sprintf(fname2,"errors%d.bin",lstep);
  if (RequestedField == 5) sprintf(fname2,"velbar%d.bin",lstep);
  if (RequestedField == 6) sprintf(fname2,"ybar%d.bin",lstep);
  if (RequestedField == 7) sprintf(fname2,"wss%d.bin",lstep);
  f3 = fopen(fname2,"wb");

  /* element type */
  if (nen == 4) eltype = "tetrahedra";
  else if (nen == 3) eltype = "triangles";
  else eltype = "cubes";
     
  /* write the positions array */
  count=0;
  for(i=0; i < numnp; i++){
    for(j=0; j < nsd; j++) {
      fwrite((void *)&(x[count]),sizeof(float),1,f2);
      count++;
    }
  }

  /* write the connections array */
  count =0;
  for(i=0; i < numel; i++) {
    for(j=0; j < nen; j++) {
      fwrite((void *)&(ien[count]),sizeof(int),1,f1);
      count++;
    }
  }

/* write the data array */
  count = 0;
  for(i=0; i< numnp; i++){
    for(j=0; j < nv; j++){
      fwrite((void *)&(q[count]),sizeof(float),1,f3);
      count++;
    }
  }
  
  /*  write the header file */
  fprintf(f0,"object 1 class array type float rank 1 shape %d items %d msb binary\n",
	  nsd,numnp);
  fprintf(f0," data file %s/points.bin,0 \n",getcwd(NULL,128));
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 2 class array type int rank 1 shape %d items %d msb binary\n",
	  nen,numel);
  fprintf(f0," data file %s/connect.bin,0 \n",getcwd(NULL,128));
  fprintf(f0," attribute \"element type\" string \"%s\" \n",eltype);
  fprintf(f0," attribute \"ref\" string \"positions\" \n\n");

  fprintf(f0,"object 3 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/solution%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 4 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/acceleration%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 5 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/bflux%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 6 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/stats%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 7 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/errors%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 8 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/velbar%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 9 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/ybar%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f0,"object 10 class array type float rank 1 shape %d items %d msb binary\n",
	  nv,numnp);
  fprintf(f0," data file %s/wss%d.bin,0 \n",getcwd(NULL,128),lstep);
  fprintf(f0," attribute \"dep\" string \"positions\" \n\n");



  fprintf(f0,"object 11 class array type int rank 0 items %d \n",
	  numel);
  fprintf(f0," data file %s/partit.out \n",getcwd(NULL,128));
  fprintf(f0," attribute \"dep\" string \"connections\" \n\n");

  fprintf(f0,"object \"irregular positions irregular connections\" class field\n");
  fprintf(f0," component \"positions\" value 1\n");
  fprintf(f0," component \"connections\" value 2\n");
  fprintf(f0," component \"data\" value %d \n",(RequestedField + 3));
  fprintf(f0,"\n end \n");

  printf("\n Files successfully written... \n");
  printf("\n Number of elements = %d \n",numel);
  printf(" Number of nodes = %d \n\n",numnp);


  fclose(f0);
  fclose(f1);
  fclose(f2);
  fclose(f3);
}
