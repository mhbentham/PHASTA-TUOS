#ifndef SIM
#include "myAttribute.h"
#include <map>
#include <strstream>

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

void IBCinput(char* filename, pGModel model, vector<pAttribute> &AttList)
{
    //open the input file
    string fname = string(filename);
    ifstream infile(fname.c_str(), ios::in);
    
    if(!infile || infile.eof()){
        cerr<<"Attribute file doesnot exist! \n";
        exit(-1);
    }
    //allocate memory for input text
    vector<string> *input_text = new vector<string>;

    // get the lines of text from the input file
    // and build the attlist from the input text
    get_input_lines(input_text, infile);
    buildattlist(input_text, model, AttList);

    delete input_text;
}

void get_input_lines(vector<string> *input_text, ifstream &infile)
{
    string textline;
    while(getline(infile,textline,'\n')){
        //ignore everything on a comment line
        if(textline[0]!='#'){
            input_text->push_back(textline);
        }
    }
}

void buildattlist(vector<string> *input_text, pGModel model, vector<pAttribute> &AttList)
{
    map<const char*, int, ltstr> numvalues;

    numvalues["initial velocity"] = 3;
    numvalues["comp1"] = 4;
    numvalues["comp3"] = 4;
    numvalues["traction vector"]=3;


    for(unsigned int i=0;i<input_text->size();i++) {
        string textlineAll = (*input_text)[i];
        string textline;
        string bctype, taginfo, diminfo, valueinfo;
        char BCtype[50];
        int numvalue, dim, tag;
        double *value;
        unsigned int pos = 0;

        if( (pos = textlineAll.find_first_of('#',pos))!=string::npos){
            textline = textlineAll.substr(0,pos);
        }else{
            textline = textlineAll;
        }

        pos = 0;

        if((pos = textline.find_first_of(':',pos)) !=string::npos){
            bctype = textline.substr(0,pos);
            strcpy(BCtype,bctype.c_str());
            if(numvalues.find(BCtype)!=numvalues.end())
                numvalue = numvalues[BCtype];
            else numvalue = 1;

            pos+=2;
            int pos2 = pos;
            if((pos2 = textline.find_first_of(' ',pos2))!=string::npos){
                taginfo = textline.substr(pos,pos2-pos);
                char Taginfo[50];
                strcpy(Taginfo,taginfo.c_str());
                tag = atoi(Taginfo);
            }else{
                cout<<"cannot find out the tag\n";
                exit(-1);
            }
            pos=pos2+1;
            pos2 = pos;
            if((pos2 = textline.find_first_of(" ",pos2))!=string::npos){
                diminfo= textline.substr(pos,pos2-pos);
                char Diminfo[50];
                strcpy(Diminfo,diminfo.c_str());
                dim = atoi(Diminfo);
            }else{
                cout<<"cannot find out the dimension\n";
                exit(-1);
            }
            value = new double[numvalue];
            for(int i = 0;i<numvalue;i++){
                pos=pos2+1;
                pos2 = pos;
                if(i==numvalue-1){
                    valueinfo = textline.substr(pos,textline.length()-pos);
                    char Valueinfo[50];
                    strcpy(Valueinfo,valueinfo.c_str());
                    value[i] = atof(Valueinfo);
                }else{
                    if((pos2 = textline.find_first_of(" ",pos2))!=string::npos){
                        valueinfo = textline.substr(pos,pos2-pos);
                        char Valueinfo[50];
                        strcpy(Valueinfo,valueinfo.c_str());
                        value[i] = atof(Valueinfo);
                    }else{
                        cout<<"cannot find out values \n";
                        exit(-1);
                    }
                }
            }

            Attribute *att = new Attribute(dim,tag,BCtype);
            AttList.push_back(att);
            pGEntity gent = GM_entityByTag(model, dim, tag);
            GEN_attachDataP(gent,BCtype,value);

        }            
    }
}
#endif
