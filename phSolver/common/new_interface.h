#ifndef __NEW_INTERFACE_H__
#define __NEW_INTERFACE_H__

#include <FCMangle.h>

#define igetMinMaxAvg FortranCInterface_GLOBAL_(igetminmaxavg,IGETMINMAXAVG)
#define rgetMinMaxAvg FortranCInterface_GLOBAL_(rgetminmaxavg,RGETMINMAXAVG)
#define print_mesh_stats FortranCInterface_GLOBAL_(print_mesh_stats,PRINT_MESH_STATS)
#define print_mpi_stats FortranCInterface_GLOBAL_(print_mpi_stats,PRINT_MPI_STATS)
#define print_system_stats FortranCInterface_GLOBAL_(print_system_stats,PRINT_SYSTEM_STATS)
#define Write_Restart FortranCInterface_GLOBAL_(write_restart,WRITE_RESTART)
#define Write_Error   FortranCInterface_GLOBAL_(write_error,WRITE_ERROR)
#define Write_Displ   FortranCInterface_GLOBAL_(write_displ,WRITE_DISPL)
#define Write_Field   FortranCInterface_GLOBAL_(write_field,WRITE_FIELD)
#define Write_d2wall  FortranCInterface_GLOBAL_(write_d2wall,WRITE_D2WALL)
#define read_d2wall FortranCInterface_GLOBAL_(read_d2wall,READ_D2WALL)
#define Write_GradPhi   FortranCInterface_GLOBAL_(write_gradphi,WRITE_GRADPHI)
#define Write_PrimVert   FortranCInterface_GLOBAL_(write_primvert,WRITE_PRIMVERT)

extern char phasta_iotype[80];
extern int field_flag;
extern int f_descriptor;

void igetMinMaxAvg(int *ivalue, double *stats, int *statRanks);
void rgetMinMaxAvg(double *value, double *stats, int *statRanks);
void print_mesh_stats(void);
void print_mpi_stats(void);
void print_system_stats(double *tcorecp);

void countfieldstowriterestart();
void
Write_Restart(  int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
	        double* soltime,
                double* array1,
                double* array2 );

void
Write_Error(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 );

void
Write_Displ(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 );

void
Write_Field(  int *pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno );

void
Write_d2wall(   int* pid,
                int* numnp,
                int* numVars,
                double* array1 );	


void
read_d2wall(  int* pid,
              int* numnp,
              int* numVars,
              double* array1,
              int* foundd2wall );

void
Write_GradPhi(  int* pid, 
                int* stepno, 
                int* nshg, 
                int* numVars,
	        double* x,
                double* y, 
                double* gradphi,
                double* mag );

void
Write_PrimVert( int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
                double* primvertval );

#endif //header guard
