/* MAXNUMVAR is the maximum limit set on the number of variables   */
/* MAAXARRAY is the maximum size of the working array in reduce.cc */
/* and geom_connectivity.cc                                        */
#define MAXNUMVAR 20
#define MAAXARRAY 20

int processCommandLineArguments(int argc, char* argv[], int* RequestedField,
				int* RequestedOutput, bool* RequestedVolCheck,
				int* stepnumber);

int geometry_and_connectivity(int* array, double* &xglobal, int* &ien, 
                               int** &ncorp2d, bool RequestedVolCheck);
			       
void asciiOutput(int* array, double* qglobal, double* xglobal, 
                 int* ien, int RequestedField);

void dxOutput(int* array, double* qglobal, double* xglobal, 
              int* ien, int RequestedField);

void tecplotOutput(int* array, double* qglobal, double* xglobal, 
                   int* ien, int RequestedField);

void phastaOutput(int* array, double* qglobal, double* aglobal, 
                  int RequestedField);
		  
void reduce(int numvar, int lsize, int gsize, int proc,
	    double* local, double* global, int** ncorp,
	    double* qmax, double* qmin);

void wrtc_(int nv, int numel, int* ien, int lstep, int nsd, int numnp, 
           float* x, float* q, int nen, int RequestedField);


void  bzero(void* ptr, size_t sz) ;
