#ifndef __NSPRE_DATA__
#define __NSPRE_DATA__

#include <map>
#include <vector>
#ifdef SIM
#include "MeshSim.h"
#include "MeshSimInternal.h"
#include "SimAttribute.h"
#else
#include "AOMD.h"
#include "AOMDInternals.h"
#include "myAttribute.h"
#ifdef MESHMODEL
#include "modeler.h"
#include "ModelTypes.h"
#endif
#endif
using namespace std;

class dofInfo {
public:
    std::vector< int >  dof;   
    inline dofInfo( ) { dof.clear(); }
} ;


/* This structure which uniquely identifies each block of elements 
   this includes the number of vertices to identify the topology and
   the polynomial order for the block. */

struct blockKey {
    int nen;
    int maxpoly;
    int nenbl;
    int lcsyst;
};
typedef struct blockKey blockKey;

struct lessKey{
    bool operator()(blockKey b1, blockKey b2) const {
            return ( 1000*b1.nen+b1.maxpoly+10*b1.nenbl+100*b1.lcsyst <
                     1000*b2.nen+b2.maxpoly+10*b2.nenbl+100*b2.lcsyst );
        }
};

typedef struct lessKey lessKey;
typedef struct pidKey pidKey;
typedef map<blockKey, int, lessKey> forwardblock;


struct EnsaArrays {
    int** ien; 
    int** ienb;
    int** ien_sms;
    int** ienb_sms;
    int** ief;
    int** iefb;
/*      int** ief_sms; */
/*      int** iefb_sms; */
    
    int* nBC; 
    int* iBC;
    int*** iBCB;
    int* iper;
    double* x; 
    double** q; 
    double*** BCB; 
    double** BC;
};
typedef struct EnsaArrays EnsaArrays;


struct globalInfo{

    int nenmax;
    int nedgemodes;
    int nfacemodes;
    int numpbc;			               /* # of prescribed essential bc's */
    int nshg;			               /* # of global modes */
    int numelb;			               /* # of boundary elements */
    int nsd;
    bool fake_mode;                    /* # false unless we need a fake mode */
};
typedef struct globalInfo globalInfo;


inline int 
setbit (int n, int i) { return n | (01 << i); }

inline int 
getbit (int n, int i) { return 01 & (n >> i); }

#endif
