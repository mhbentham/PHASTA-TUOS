/* create our attribute system. Only useful when using SCOREC libraries

   Min Zhou spring 2008
*/

#ifndef SIM
#include "myAttribute.h"
#ifdef MESHMODEL
#include "modelerDiscrete.h"
#endif
#include <vector>

using namespace std;
vector<pAttribute> AttList;
pGModel meshmodel;

pGModel AttCase_model(pACase acase){
    return acase;
}

pAManager AMAN_load(char* fname){

    meshmodel =GM_createFromDmgFile("geom.dmg");
    IBCinput(fname, meshmodel, AttList);    
    return meshmodel;
}

void AttCase_associate(pACase acase){
    cout<<"AttCase_associate not implemented yet\n";
}

pACase AMAN_findCase(pAManager attmngr, const char* casename){
    return (pACase)attmngr;
}
void AMAN_release(pAManager attmngr){
    
    int EntityTag, dimension;
    double *value;
    meshmodel = attmngr;
    for(unsigned int i=0;i<AttList.size();i++){
        EntityTag = AttList[i]->tag;
        dimension = AttList[i]->dim;
        pGEntity gent = GM_entityByTag(meshmodel,dimension,EntityTag);
        if(GEN_dataP(gent,AttList[i]->BCtype,(void**)&value)!=0){
            delete [] value;
            GEN_removeData(gent,AttList[i]->BCtype);
        }
        delete AttList[i];
    }        
}

pAttribute GM_attrib(pGModel gmodel, char* Conditiontype){
    for(unsigned int i=0;i<AttList.size();i++)
        if(strcmp(AttList[i]->BCtype,Conditiontype)==0 && AttList[i]->dim==3)
            return AttList[i];
    return NULL;

}

pAttribute GEN_attrib(pGEntity gent, char* Conditiontype){
    int Entitytag = GEN_tag((pGEntity)gent);
    for(unsigned int i=0;i<AttList.size();i++)
        if(strcmp(AttList[i]->BCtype,Conditiontype)==0 && AttList[i]->dim==GEN_type(gent) && AttList[i]->tag == Entitytag)
            return AttList[i];
    return NULL;
}

double AttributeTensor0_evalDS(pAttributeTensor0 attList, double* x){
    double *value;
    pGEntity gent = GM_entityByTag(meshmodel,attList->dim,attList->tag);

    GEN_dataP(gent,attList->BCtype,(void**)&value);
    return *value;
}

double AttributeTensor1_evalDS(pAttributeTensor1 attList, int i, double* x){
    double *value;
    pGEntity gent = GM_entityByTag(meshmodel,attList->dim,attList->tag);

    GEN_dataP(gent,attList->BCtype,(void**)&value);
    return value[i];
}

double AttributeDouble_value(pAttributeDouble att){
    cout<<"AttributeDouble_value is not implemented yet \n";
    return 0.0;
}

double AttributeDouble_evalDS(pAttributeDouble attList, double* x){
    double *value;
    pGEntity gent = GM_entityByTag(meshmodel,attList->dim,attList->tag);
    GEN_dataP(gent,attList->BCtype,(void**)&value);
    return (*value)*(1-x[0]*x[0]-x[1]*x[1]);
}

pAttribute Attribute_childByType(pAttribute attList,char* request){

    double *value;
    Attribute *att;
    pGEntity gent = GM_entityByTag(meshmodel,attList->dim,attList->tag);

    GEN_dataP(gent,attList->BCtype,(void**)&value);
    
    if(strcmp(request,"vector direction")==0){
        double *dir = new double[3];
        for(int i=0;i<3;i++)
            dir[i] = value[i+1];
        GEN_attachDataP(gent,"vector direction", dir);
        att = new Attribute(attList->dim,attList->tag,"vector direction");
        AttList.push_back(att);
    }
    else if(strcmp(request,"vector magnitude")==0){        
        double *mag = new double;
        *mag = value[0];
        GEN_attachDataP(gent,"vector magnitude", mag);
        att = new Attribute(attList->dim,attList->tag,"vector magnitude");
        AttList.push_back(att);      
    }
    return att;
        
}
int AttributeInt_value(pAttributeInt att){
    return att->tag;
}
#endif
