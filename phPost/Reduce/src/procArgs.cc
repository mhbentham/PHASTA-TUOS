/********************************************************************/
/* Processing command line arguments                                */
/*                                                                  */
/* Elaine Bohr                                                      */
/* February 2004                                                    */
/********************************************************************/
#include <iostream>
#include <stdio.h>
#include <string>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

using namespace std;
char fieldname[255]="errors";

int processCommandLineArguments(int argc, char* argv[], int* RequestedField,
				int* RequestedOutput, bool* RequestedVolCheck,
				int* stepnumber){
    FILE *nstart;
    int iarg, arglength;
    string tmpstr;
    bool StepNumberAvailable;
    bool BogusCmdLine = false;
    bool RequestedHelp = false;
    bool RequestedStepNumber = false;

    /* argc is the number of strings on the command-line */
    /*  starting with the program name */
    for(iarg=1; iarg<argc; iarg++){
        arglength = strlen(argv[iarg]);
        /* replace 0..arglength-1 with argv[iarg] */
        tmpstr.replace(0,arglength,argv[iarg],0,arglength);
        if(tmpstr=="-h"){
            RequestedHelp = true;
            cout << endl;
            cout << "usage:" <<endl;
            cout << "  Reduce [-sn stepnumber][-(output format)][-(requested field)]" << endl;
            cout << endl;
            cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
            cout << "  -h                : Display usage and command-line argument summary"<< endl;
            cout << "  -sn stepnumber    : Specify step number to reduce"<< endl;
            cout << "  -vol              : Calculate element volumes and write volcheck.dat" <<endl;
	    cout << endl;
	    cout << "Possible output formats" << endl;
            cout << "  -ph               : Write phasta-format file restart.<stepnumber>.0"<<endl;
            cout << "  -dx               : Write files for visualization in Data Explorer" <<endl;
            cout << "  -ASCII            : Write ASCII-format restart and geometry files" <<endl;
	    cout << "  -tec              : Write files for visualization in Tecplot" << endl;
	    cout << endl;
	    cout << "Possible fields to be reduced" << endl;
            cout << "  -y_ac             : Reduce velocity and time-derivative fields"<< endl;
            cout << "  -bflux            : Reduce boundary fluxes"<<endl;
            cout << "  -stats            : Reduce turbulent statistics"<<endl;
            cout << "  -errors           : Reduce error field"<<endl;	 
            cout << "  -velbar           : Reduce velbar field"<<endl;	    
            cout << "  -ybar             : Reduce ybar field"<<endl;	
            cout << "  -wss              : Reduce wall shear stress field"<<endl;	 
            cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
            cout << endl;
        }
        else if(tmpstr=="-sn"){
            RequestedStepNumber = true;
            iarg++;
            *stepnumber = atoi(argv[iarg]);
            StepNumberAvailable = true;
        } 
        else if(tmpstr=="-tag"){
            iarg++;
            strcpy( fieldname,argv[iarg]);
        }
	/* Here are the possible fields from restart files             */
	/* that can be reduced from multiproc case to one proc format. */
	/* As default the solution is reduced.                         */
	else if(tmpstr=="-y_ac"){
	    *RequestedField = 1;
        }
        else if(tmpstr=="-bflux"){
	    *RequestedField = 2;
        }
	else if(tmpstr=="-stats"){
	    *RequestedField = 3;
        }
	else if(tmpstr=="-errors"){
	    *RequestedField = 4;
        }
	else if(tmpstr=="-velbar"){
	    *RequestedField = 5;
        }
	else if(tmpstr=="-ybar"){
	    *RequestedField = 6;
        }
	else if(tmpstr=="-wss"){
	    *RequestedField = 7;
        }        
	else if(tmpstr=="-vol"){
            *RequestedVolCheck = true;
        }
	
	/* Here are the different output possible.*/
        else if(tmpstr=="-dx"){
            *RequestedOutput = 1;
        }
        else if(tmpstr=="-ph"){
            *RequestedOutput = 0;
        }
        else if(tmpstr=="-ASCII"){
            *RequestedOutput = 2;
        }
	else if(tmpstr=="-tec"){
            *RequestedOutput = 3;
        }
        
	else {
            BogusCmdLine = true;
        }
        /* reset tmpstr for next argument */
        tmpstr.erase(0,arglength);
    }
    if(BogusCmdLine){
        cout << "Use -h option for help" << endl;
	return(1);
    }
    if(RequestedHelp){
        cout << endl;
        cout << "Exiting due to -h flag";
        cout << endl;
        return(1);
    }
    if(RequestedStepNumber){
        cout << "Will use requested step number " << *stepnumber << endl;
    } else {
        /* User did not make a request */
        /* Try to use the number in numstart.dat */
        nstart = fopen("numstart.dat","r");
        if(nstart){
            /* The file numstart.dat is present*/
            fscanf(nstart,"%d", stepnumber);
            fclose(nstart);
            StepNumberAvailable=true;
            cout << "Will use numstart.dat step number " << *stepnumber << endl;
        }
    }
    if(*RequestedField == 1){
        cout << "Will reduce time-derivative field as requested" << endl;
    }
    if(*RequestedField == 0){
        cout << "Will reduce solution field" << endl;
    }
    if(*RequestedOutput == 1){
        cout << "Will write files for visualization in Data Explorer as requested" << endl;
    }
    if(*RequestedVolCheck){
        cout << "Will calculate element volumes and write volcheck.dat as requested" << endl;
    }
    if(*RequestedOutput == 2){
        cout << "Will write ASCII restart and geometry as requested" << endl;
    }
    if(*RequestedField == 2){
        cout << "Will reduce boundary flux field as requested" << endl;
    }
    if(*RequestedField == 3){
        cout << "Will reduce turbulent statistics as requested" << endl;
    }
    if(*RequestedField == 4){
        cout << "Will reduce error field as requested" << endl;
    }
    if(*RequestedField == 5){
        cout << "Will reduce velbar field as requested" << endl;
    }
    if(*RequestedField == 6){
        cout << "Will reduce ybar field as requested" << endl;
    }
    if(*RequestedField == 7){
        cout << "Will reduce wss field as requested" << endl;
    }

    /* Keep these last */
    if(!StepNumberAvailable){
        cout << "No step number or range of steps given, so exiting." << endl;
        return(1);
    }
    return(0);
}
