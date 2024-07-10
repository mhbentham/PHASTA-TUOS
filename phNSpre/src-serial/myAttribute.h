/* create our own attribute system
  Min Zhou spring 2008
*/

#ifndef SIM
#ifndef _MYATTRIBUTE_
#define _MYATTRIBUTE_
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "AOMD.h"
#include "AOMDInternals.h"
#ifdef MESHMODEL
#include "modeler.h"
#endif

using namespace std;

class Attribute {
 public:
    int dim;
    int tag;
    char BCtype[50];
    Attribute(int dimension, int entitytag, char* conditiontype){
        dim = dimension;
        tag = entitytag;
        sprintf(BCtype,"%s",conditiontype);
    }
};

typedef class Attribute * pAttribute;
typedef class Attribute * pAttributeTensor0;
typedef class Attribute * pAttributeTensor1;
typedef class Attribute * pAttributeDouble;
typedef class Attribute * pAttributeInt;
typedef pGModel pACase;
typedef pGModel pAManager;

pGModel AttCase_model(pACase acase);

pAManager AMAN_load(char* fname);
void AttCase_associate(pACase acase);
pACase AMAN_findCase(pAManager attmngr, const char* casename);
void AMAN_release(pAManager attmngr);

pAttribute GM_attrib(pGModel gmodel, char* Conditiontype);
pAttribute  GEN_attrib(pGEntity gent, char* strAtt);

double AttributeTensor0_evalDS(pAttributeTensor0 attList, double* x);
double AttributeTensor1_evalDS(pAttributeTensor1 attList, int i, double* x);
double AttributeDouble_value(pAttributeDouble att);
double AttributeDouble_evalDS(pAttributeDouble att, double* x);
pAttribute Attribute_childByType(pAttribute AttList,char* request); 
int AttributeInt_value(pAttributeInt AttList);

void IBCinput(char* fname,pGModel model,vector<pAttribute> &AttList);
void get_input_lines(vector<string> *input_text, ifstream &infile);
void buildattlist(vector<string> *input_text, pGModel model, vector<pAttribute> &AttList);


#endif
#endif
