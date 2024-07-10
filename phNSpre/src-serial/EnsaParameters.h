//
// class for ENSA parameters
//
#include "nspre_data.h"

#ifndef _EnsaParameters
#define _EnsaParameters

class EnsaParameters {
public:
    EnsaParameters ( int nsd, int nen,
		           int numflx,  int numVars, int numEBC,
		           int numNBC,  int fType ) {
  _nsd = nsd;
  _nen = nen;

  _numflx = numflx;
  _numVars = numVars;
  _numEBC = numEBC;
  _numNBC = numNBC;
  _fType = fType;

}  
  int getNSD() { return _nsd; }
  int getNEN() { return _nen; }

  int getNUMFLX() { return _numflx; }
  int getNUMVARS() { return _numVars; }
  int getNUMEBC() { return _numEBC; }
  int getNUMNBC() { return _numNBC; }
  int getFTYPE() { return _fType; }

private:
  int _nsd, _nen,  _numflx, _numVars, _numEBC, _numNBC, _fType;
};


#endif
