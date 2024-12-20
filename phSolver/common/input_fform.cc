#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cstring>

#include "Input.h"
#include "common_c.h"

using namespace std; //::cout;
void print_error_code(int ierr);
int input_fform(char inpfname[]);
int SONFATH=0;
extern "C" char phasta_iotype[80];

extern "C" {

  int input_fform_(char inpfname[]) {
    return input_fform(inpfname);
  }

}

int input_fform(char inpfname[])
{

  int ierr = 0 ;
  int i,j;
  char* path_to_config = 0 ;
  char complete_filename[256];

  try {
    // get the input file stream
    path_to_config = getenv("PHASTA_CONFIG");
    if(path_to_config) strcpy(complete_filename, path_to_config);
    else strcpy(complete_filename,".");
    strcat(complete_filename, "/input.config");
    if(workfc.myrank==workfc.master) {
      printf("\n Complete Filename: %s \n", complete_filename);
      printf("\n Local Config: %s \n\n", inpfname);
    }
    string def(complete_filename);
    Input inp(inpfname,def);


	// AMG PARAMETERS 
	
	if ((string)inp.GetValue("Employ AMG") == "True" ) {
	  
	    amgvari.irun_amg = 1;

        amgvari.irun_amg_sa = inp.GetValue("Run AMG Stand-alone");
        amgvari.irun_amg_prec = inp.GetValue("Run AMG As CG-preconditioner");

        amgvarr.strong_eps       = inp.GetValue("Strong Criterion Eps");
        
        amgvarr.ramg_eps         = inp.GetValue("AMG Convergence Eps");
        
        amgvarr.ramg_relax       = inp.GetValue("AMG Relaxation Omega");        
        amgvarr.ramg_trunc       = inp.GetValue("AMG Truncation Set");
        
        amgvari.iamg_verb        = inp.GetValue("AMG Verbosity");
        
        amgvari.iamg_neg_sten    = inp.GetValue("AMG Neg_Sten");           

        amgvari.iamg_nlevel      = inp.GetValue("AMG Nlevel"); 
        
        amgvari.iamg_c_solver  = inp.GetValue("AMG Coarsest Solver");
        
        amgvari.iamg_prescale     = inp.GetValue("AMG Prescaling");

        amgvari.iamg_init = 0;
        amgvari.iamg_ppe_frez = inp.GetValue("AMG Freeze PPE");
        amgvari.iamg_setup_frez = inp.GetValue("AMG Freeze Setup");
        amgvari.iamg_iai_frez = inp.GetValue("AMG Freeze Construct Coarse");
        if ((string)inp.GetValue("AMG Interpolation Type")=="Standard")
            amgvari.iamg_interp = 1;
        else
            amgvari.iamg_interp = 0;
        amgvari.maxnev = inp.GetValue("AMG GGB nev");
        amgvari.maxncv = inp.GetValue("AMG GGB ncv");
        amgvari.mlsdeg = inp.GetValue("AMG MLS Degree");
        amgvari.iamg_scale = inp.GetValue("AMG PPE Scale");
        amgvari.iamg_reduce = inp.GetValue("AMG Run Reduced Serial");
	}	


    // Disabled Features 

    conpar.iALE = inp.GetValue("iALE");
    conpar.icoord = inp.GetValue("icoord");
    conpar.irs = inp.GetValue("irs");
    conpar.iexec = inp.GetValue("iexec");
    timpar.ntseq = inp.GetValue("ntseq");
    solpar.imap = inp.GetValue("imap");

//-------------------------------------------------------------------
// Magnus, need the follow block for the bubble controller 
    //Control Force
    bubstudy.iClrLiq = inp.GetValue("Color local liquid");
    bubstudy.iBK = inp.GetValue("Break up tracking");
    bubstudy.iCForz = inp.GetValue("Matts Control Force");
    bubstudy.iCForz_where = inp.GetValue("Apply to Whole Domain");
//    std::cout << "Value of iCForz is: " << bubstudy.iCForz << std::endl;

    if (bubstudy.iCForz == 1) {
       vector<double> xcfcoeffvec;
       vector<double> ycfcoeffvec;
       vector<double> zcfcoeffvec;
       xcfcoeffvec = inp.GetValue("X Control Force Coefficients");
       ycfcoeffvec = inp.GetValue("Y Control Force Coefficients");
       zcfcoeffvec = inp.GetValue("Z Control Force Coefficients");
       for(i=0; i<10; i++){
          bubstudy.xcfcoeff[i]=xcfcoeffvec[i];
       }
       for(i=0; i<9; i++){
          bubstudy.ycfcoeff[i]=ycfcoeffvec[i];
          bubstudy.zcfcoeff[i]=zcfcoeffvec[i];
       }
       xcfcoeffvec.erase(xcfcoeffvec.begin(),xcfcoeffvec.end());
       ycfcoeffvec.erase(ycfcoeffvec.begin(),ycfcoeffvec.end());
       zcfcoeffvec.erase(zcfcoeffvec.begin(),zcfcoeffvec.end()); 
       
       bubstudy.numts_histyavg = inp.GetValue("Y History Term Averaging Range");
       bubstudy.shear_rate = inp.GetValue("Shear Rate");
       bubstudy.vel_centre = inp.GetValue("Velocity at Centerline"); }
    else {
       for(i=0; i<10; i++){
          bubstudy.xcfcoeff[i]=0.0;
       }
       for(i=0; i<9; i++){
          bubstudy.ycfcoeff[i]=0.0;
          bubstudy.zcfcoeff[i]=0.0;
       }
    }
//------------------------------------------------------------------------------

    // bubble study flags
    //bubstudy.iClrLiq = inp.GetValue("Local Liquid Tracking");
    //bubstudy.iBK = inp.GetValue("Break-up Tracking");
    bubstudy.phi_inner = inp.GetValue("Liquid Shell Inner Front Location");
    bubstudy.phi_outer = inp.GetValue("Liquid Shell Outer Front Location");

    // Solution Control Keywords

    if((string)inp.GetValue("Equation of State") == "Incompressible") matdat.matflg[0][0] =-1 ;
    if((string)inp.GetValue("Equation of State") == "Compressible") matdat.matflg[0][0] =0;
    inpdat.Delt[0] = inp.GetValue("Time Step Size");
    inpdat.timestart = inp.GetValue("Beginning Time");
    inpdat.nstep[0] = inp.GetValue("Number of Timesteps");
    if((string)inp.GetValue("Time Step Based on Maximum CFL") == "True") {
      inpdat.iflag_cfl_dt = 1;
      inpdat.CFLfl_limit = inp.GetValue("Maximum CFL Number Allowed");
    }
    else {
      inpdat.iflag_cfl_dt = 0;
    }
    if((string)inp.GetValue("Viscous Control")=="Viscous") conpar.navier=1 ; else conpar.navier=0;
   
      turbvari.idns = 0; 
    if ((string)inp.GetValue("Turbulence Model") == "No-Model" ) {
      turbvari.irans = 0;
      turbvari.iles  = 0;
    } else if ((string)inp.GetValue("Turbulence Model") == "DNS-WallFunc" ) {
      turbvari.irans = 0;
      turbvari.iles  = 0;
      turbvari.idns = 1;
    } else if ((string)inp.GetValue("Turbulence Model") == "LES" ) {
      turbvari.iles  = 1;
      turbvari.irans = 0;
    } else if ((string)inp.GetValue("Turbulence Model") == "RANS-SA" ) {
      turbvari.iles  = 0;
      turbvari.irans = -1;
    } else if ((string)inp.GetValue("Turbulence Model") == "RANS" ) {
      turbvari.iles  = 0;
      turbvari.irans = -1; // assume S-A for backward compatibility
    } else if ((string)inp.GetValue("Turbulence Model") == "RANS-KE" ) {
      turbvari.iles  = 0;
      turbvari.irans = -2;
    } else if ((string)inp.GetValue("Turbulence Model") == "DES" ) {
      turbvari.iles  = 1;
      turbvari.irans = -1;
    } else {
    if (workfc.myrank==workfc.master) {
      cout << " Turbulence Model: Only Legal Values ( No-Model, LES, RANS-SA, RANS-KE, DES, DNS-WallFunc )";
      cout << endl;
      cout << " Turbulence Model:  the code got: " << (string)inp.GetValue("Turbulence Model");
      cout << endl;
     }
      exit(1);
    }

    if (turbvari.iles*turbvari.irans!=0) turbvar.eles=
                                           inp.GetValue("DES Edge Length");

    int solflow, solheat , solscalr, ilset, coalcon, coaltimtrak, coalest, iBT, id2w, icoalCtrl;
    ((string)inp.GetValue("Solve Flow") == "True")? solflow=1:solflow=0;
    ((string)inp.GetValue("Solve Heat") == "True")? 
    pcboiling.solheat=1:pcboiling.solheat=0;
    levlset.coalcon=0;
    levlset.coaltimtrak=1;
    levlset.coalest=2;
    levlset.coalbubrad=5.0E-3;
    levlset.epsilonBT=3.0E-5;
    levlset.iBT=0;
    ((string)inp.GetValue("d2wall Calculation") == "Enable")? levlset.id2w=1:levlset.id2w=0;
    levlset.icoalCtrl=0;
    //for compressible solheat= False so
    if((string)inp.GetValue("Equation of State") == "Compressible") pcboiling.solheat=0;
    ilset = (int)inp.GetValue("Solve Level Set");
    solscalr = (int)inp.GetValue("Solve Scalars");
    solscalr += ilset;
    if(turbvari.irans == -1) solscalr++;
    if(turbvari.irans == -2) solscalr=solscalr+2;
    if ( solscalr > 4 ) {
      cout << " Only Four Scalars are supported \n";
      cout <<" Please reduce number of scalars \n";
      exit(1);
    }
    inpdat.impl[0] = 10*solflow+solscalr*100+pcboiling.solheat;

    levlset.iLSet = ilset;
    if( ilset > 0) {
    levlset.epsilon_ls = inp.GetValue("Number of Elements Across Interface");
    levlset.epsilon_lsd = inp.GetValue("Number of Elements Across Interface for Redistancing");
    levlset.dtlset = inp.GetValue("Pseudo Time step for Redistancing");
    levlset.dtlset_cfl = inp.GetValue("Base pseudo time step for redistancing on CFL number");
    ((string)inp.GetValue("Bubble Tracking") == "Enable")?  levlset.iBT=1:levlset.iBT=0;
    if(levlset.iBT == 1){
       levlset.id2w = 1;
    }
    ((string)inp.GetValue("Jun Fang Coalescence Control") == "Enable")?  levlset.icoalCtrl=1:levlset.icoalCtrl=0;
    levlset.icoalcon_verbo = inp.GetValue("Jun Fang Coalescence Control Verbosity");
    levlset.coalest = inp.GetValue("Number of Estimated Coalescence Events");
    levlset.coalcon_dist = inp.GetValue("Jun Fang Coalescence Control Interface Distance");
    levlset.epsilonBT = inp.GetValue("Bubble Tracking Epsilon");
    levlset.coalbubrad = inp.GetValue("Bubble Radius for Coalescence Control");
    ((string)inp.GetValue("Coalescence Control") == "Active")? levlset.coalcon=1:levlset.coalcon=0;
    ((string)inp.GetValue("Coalescence Time Control") == "Active")? levlset.coaltimtrak=0:levlset.coaltimtrak=1;
       if (levlset.dtlset_cfl > 0.0) {
       levlset.i_dtlset_cfl = 1; }
    else {
       levlset.i_dtlset_cfl = 0;
    }
    levlset.AdjRedistVelCFL = inp.GetValue("Adjust Redistance Velocity to Satisfy CFL Limit");
    if (levlset.AdjRedistVelCFL > 0.0) {
       levlset.i_AdjRedistVel = 1; }
    else { 
       levlset.i_AdjRedistVel = 0;
    }
    levlset.iSolvLSSclr2 = inp.GetValue("Solve for Redistance Field");
    levlset.iSolvLSSclr1 = inp.GetValue("Solve for Scalar 1 Field");
    levlset.i_focusredist = inp.GetValue("Focus redistancing about interface");
    if ((string)inp.GetValue("Apply Volume Constraint") == "True" ) {
      levlset.ivconstraint = 1; } 
    else if((string)inp.GetValue("Apply Volume Constraint") == "False" ) {
      levlset.ivconstraint = 0; }
    else {
      cout << "Apply Volume Constraint: Only Legal Values (True, False) ";
      cout << endl;
      exit(1);
    }   
//*********************
     if ((string)inp.GetValue("Redistance loop") == "True" ) {
       levlset.i_redist_loop_flag = 1; }
     else if((string)inp.GetValue("Redistance loop") == "False" ) {
       levlset.i_redist_loop_flag = 0; }
     else {
       cout << "Redistance loop: Only Legal Values (True, False) ";
       cout << endl;
       exit(1);
     }
     levlset.redist_toler = inp.GetValue("Tolerance for redistance loop");
     levlset.i_redist_max_iter = inp.GetValue("Maximum number of redistance iterations");
     levlset.i_spat_var_eps_flag = inp.GetValue("Use spatial varying epsilon_ls");
//*******************
// For checking proximity of the interface to large elements
     if ((string)inp.GetValue("Check proximity of interface to large elements") == "True" ) {
       levlset.i_check_prox = 1; }
     else if((string)inp.GetValue("Check proximity of interface to large elements") == "False" ) {
       levlset.i_check_prox = 0; }
     else {
       cout << "Check proximity of interface to large elements: Only Legal Values (True, False) ";
       cout << endl;
       exit(1);
     }
     levlset.r_int_buffer = inp.GetValue("Check proximity interface buffer thickness");
     levlset.r_int_elem_size = inp.GetValue("Check proximity maximum element size");
// For output of gradient of level set function
     if ((string)inp.GetValue("Output level set gradient") == "True" ) {
       levlset.i_gradphi = 1; }
     else if((string)inp.GetValue("Output level set gradient") == "False" ) {
       levlset.i_gradphi = 0; }
     else {
       cout << "Output level set gradient: Only Legal Values (True, False) ";
       cout << endl;
       exit(1);
     }
    }

    vector<double> vec;

    // OUTPUT CONTROL KEY WORDS.
    outpar.iMeshingTool = inp.GetValue("Meshing Tools Used");  //JF, Jan 2015

    conpar.necho = inp.GetValue("Verbosity Level");
    outpar.ntout = inp.GetValue("Number of Timesteps between Restarts");
    if((string)inp.GetValue("Print Statistics") == "True") outpar.ioform = 2;
    else outpar.ioform = 1;
  
    if((string)inp.GetValue("Print Wall Fluxes") == "True") outpar.iowflux = 1;
    else outpar.iowflux = 0;

    if((string)inp.GetValue("Print FieldView") == "True") outpar.iofieldv = 1;
    else outpar.iofieldv = 0;

    if((string)inp.GetValue("Print ybar") == "True") outpar.ioybar = 1;
    else outpar.ioybar = 0;

//MR CHANGE
    if((string)inp.GetValue("Print vorticity") == "True") outpar.ivort = 1;
    else outpar.ivort = 0;
//MR CHANGE END
    strcpy( outpar.iotype , ((string)inp.GetValue("Data Block Format")).c_str());
    strcpy( phasta_iotype , ((string)inp.GetValue("Data Block Format")).c_str());
    SONFATH = inp.GetValue("Number of Father Nodes");
  
    if((string)inp.GetValue("Print Residual at End of Step") == "True") genpar.lstres = 1;
    else genpar.lstres = 0;
  
    if((string)inp.GetValue("Print Error Indicators") == "True") turbvar.ierrcalc = 1;
    else turbvar.ierrcalc = 0;

    if((string)inp.GetValue("Print Velocity Hessian") == "True") turbvar.ihessian = 1;
    else turbvar.ihessian = 0;

    if ( turbvar.ierrcalc == 1 )
        turbvari.ierrsmooth = inp.GetValue("Number of Error Smoothing Iterations");

    int nsrfCM = inp.GetValue("Number of Force Surfaces");
    if (nsrfCM > 0) {
      vector<int> ivec = inp.GetValue("Surface ID's for Force Calculation");
      for(i=0;i<MAXSURF+1; i++) aerfrc.nsrflist[i] = 0;
      for(i=0; i< nsrfCM; i++){
        aerfrc.nsrflist[ivec[i]] = 1;
        //        cout <<"surface in force list "<< ivec[i] << endl;
      }
      ivec.erase(ivec.begin(),ivec.end());
    }

    aerfrc.isrfIM = inp.GetValue("Surface ID for Integrated Mass");
    //Limiting
    vec = inp.GetValue("Limit u1");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[0][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit u2");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[1][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit u3");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[2][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit Pressure");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[3][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit Temperature");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[4][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    //Material Properties Keywords 
    matdat.nummat = levlset.iLSet+1;
    if((string)inp.GetValue("Shear Law") == "Constant Viscosity") 
      for(i=0; i < levlset.iLSet+1; i++) matdat.matflg[i][1] = 0;

    if((string)inp.GetValue("Bulk Viscosity Law") == "Constant Bulk Viscosity") 
      for(i=0; i < levlset.iLSet+1; i++) matdat.matflg[i][2] = 0;

    mmatpar.pr = inp.GetValue("Prandtl Number"); 

    mmatpar.Texp = inp.GetValue("Thermal Expansion Coefficient");

    if((string)inp.GetValue("Conductivity Law") == "Constant Conductivity") 
      for(i=0; i < levlset.iLSet+1; i++) matdat.matflg[i][3] = 0;

    vec = inp.GetValue("Density");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][0][0] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Viscosity");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][1][0] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

// This parameter set the Bubble Radius for the wall repellant force:
    if (levlset.iLSet > 0) {
    vec = inp.GetValue("Repellant Force Bubble Radius");
    levlset.BubRad = vec[0];
    levlset.vf_target = inp.GetValue("Target Constant Void Fraction");
    levlset.vfcontrcoeff = inp.GetValue("Void Fraction Control Coefficient");
    levlset.iuse_vfcont_cap = inp.GetValue("Use Void Fraction Control Coefficient Cap");
    levlset.C_int_cap = inp.GetValue("Void Fraction Control Coefficient Cap");
    if(workfc.myrank==workfc.master) {
      printf("\n Target Void Fraction is Set: %f \n", levlset.vf_target);
    }
     }

// Those parameters are needed for ramping up/down the density and viscosity (initially being developed for bubble-in-a-turbulent-channel-flow):

    matramp.iramp = 0;
    if((string)inp.GetValue("Ramp Properties") == "Yes")
      matramp.iramp = 1;

    if (matramp.iramp == 1) {
// Read the additional parameters in case we ramping the properties:
      matramp.tmu = inp.GetValue("Target Viscosity");
      matramp.trho = inp.GetValue("Target Density");
      matramp.qrts0 = inp.GetValue("Ramp Start Time");
      matramp.qrts1 = inp.GetValue("Ramp Stop Time");
// save the original rho & mu:
      matramp.omu = matdat.datmat[1][1][0];
      matramp.orho = matdat.datmat[1][0][0];
    }
// ****************************************

// Bubble evaporation and condensation model

    if ( (string)inp.GetValue("Bubbly boiling Mode") == "True"){
        pcboiling.bubboil = 1.0;  levlset.vfcontrcoeff = 0.0; 
//	else pcboiling.bubboil = 0;
    }
    if ((string)inp.GetValue("Bubble growth Mode") == "True"){
        pcboiling.bubgrow = 1.0; levlset.vfcontrcoeff = 0.0;
//        else pcboiling.bubgrow = 0;
    }
    if ( pcboiling.bubgrow == 1||pcboiling.solheat == 1){
    vec = inp.GetValue("Specific Heat");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][2][0] = vec[i];
    }
      vec.erase(vec.begin(),vec.end());
    }
     if ( pcboiling.solheat == 1)
        pcboiling.h_fg =  inp.GetValue("The latent heat");
      if (pcboiling.bubgrow == 1||pcboiling.bubboil == 1){
       pcboiling.h_fg =  inp.GetValue("The latent heat");

       vec = inp.GetValue("The superheated rate");
       for(i=0; i< 4 ; i++){
       pcboiling.delt_T[i]= vec[i];
        }
       vec.erase(vec.begin(),vec.end());
//     Mengnan
//     //       pcboiling.delt_T = inp.GetValue("The superheated rate");
            pcboiling.epsilon_lst = inp.GetValue("Number of Elements Across Interface for heat transfer");
      }
       vec = inp.GetValue("Thermal Conductivity");
      for(i=0; i< levlset.iLSet +1 ; i++){
        matdat.datmat[i][3][0] = vec[i];
        }
        vec.erase(vec.begin(),vec.end());
  	
    vec = inp.GetValue("Scalar Diffusivity");
    for(i=0; i< solscalr ; i++){
      sclrs.scdiff[i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    if((string)inp.GetValue("Zero Mean Pressure") == "True")
      turbvar.pzero=1;

    turbvar.rmutarget = inp.GetValue("Target Viscosity For Step NSTEP");

    if ( (string)inp.GetValue("Body Force Option") == "None" ) {
      for( i=0; i< levlset.iLSet +1 ; i++)  matdat.matflg[i][4] = 0;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Vector" ) {
      for( i=0; i< levlset.iLSet +1 ; i++)  matdat.matflg[i][4] = 1;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "User e3source.f" ) {
      for( i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 3;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Boussinesq" ) {
      for(i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 2;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Cooling Analytic" ) {
      for(i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 4;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Cooling Initial Condition" ) {
      for(i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 5;
    }

    // the following block of stuff is common to all cooling type sponges. 
    // Specific stuff belongs in the conditionals above

    if(matdat.matflg[0][4] >=4) {
      spongevar.betamax = inp.GetValue("Maximum Value of Sponge Parameter");
      spongevar.zinsponge = inp.GetValue("Inflow Cooling Sponge Ends at z");
      spongevar.zoutsponge= inp.GetValue("Outflow Cooling Sponge Begins at z");
      spongevar.radsponge = inp.GetValue("Radial Cooling Sponge Begins at r");
      spongevar.grthosponge = inp.GetValue("Sponge Growth Coefficient Outflow");
      spongevar.grthisponge = inp.GetValue("Sponge Growth Coefficient Inflow");


      spongevar.spongecontinuity = 0;
      spongevar.spongemomentum1 = 0;
      spongevar.spongemomentum2 = 0;
      spongevar.spongemomentum3 = 0;
      spongevar.spongeenergy = 0;
 
      if((string)inp.GetValue("Sponge for Continuity Equation") == "True")
	spongevar.spongecontinuity = 1;
      if((string)inp.GetValue("Sponge for x Momentum Equation") == "True")
	spongevar.spongemomentum1 = 1;
      if((string)inp.GetValue("Sponge for y Momentum Equation") == "True")
	spongevar.spongemomentum2 = 1;
      if((string)inp.GetValue("Sponge for z Momentum Equation") == "True")
	spongevar.spongemomentum3 = 1;
      if((string)inp.GetValue("Sponge for Energy Equation") == "True")
	spongevar.spongeenergy = 1;
      
    }

    vec = inp.GetValue("Body Force");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][4][0] = vec[0+i*3];
      matdat.datmat[i][4][1] = vec[1+i*3];
      matdat.datmat[i][4][2] = vec[2+i*3];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Body Force Pressure Gradient");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][6][0] = vec[0+i*3];
      matdat.datmat[i][6][1] = vec[1+i*3];
      matdat.datmat[i][6][2] = vec[2+i*3];
    }
    vec.erase(vec.begin(),vec.end());

    if ( (string)inp.GetValue("Surface Tension Option") == "No" ){
        genpar.isurf = 0;
    }
    else if ((string)inp.GetValue("Surface Tension Option") == "Yes" ){
        genpar.isurf = 1;
    }
    else {
      cout << " Surface Tension: Only Legal Values (Yes, No) ";
      cout << endl;
      exit(1);
    }
    if( genpar.isurf > 0) {
      genpar.Bo = inp.GetValue("Inverse Surface Tension");
     genpar.CoalInvSigma = inp.GetValue("Coalescence Control Inverse Surface Tension");
    }

    genpar.EntropyPressure = inp.GetValue("Entropy Form of Pressure Constraint on Weight Space");

    
    if ( (string)inp.GetValue("Rotating Frame of Reference") == "True" ) {
      matdat.matflg[0][5] = 1;
      vec = inp.GetValue("Rotating Frame of Reference Rotation Rate");
      matdat.datmat[0][5][0] = vec[0];
      matdat.datmat[0][5][1] = vec[1];
      matdat.datmat[0][5][2] = vec[2];
      vec.erase(vec.begin(),vec.end());
    }
    else {
      matdat.matflg[0][5] = 0;
      matdat.datmat[0][5][0] = 0.;
      matdat.datmat[0][5][1] = 0.;
      matdat.datmat[0][5][2] = 0.;
    }

    //Linear Solver parameters
    //inpdat.svLSFlag=0;    // svLS
    //inpdat.memLSFlag=0;    // memLS
    inpdat.svLSFlag=0;    // memLS
    if( (string)inp.GetValue("Solver Type") =="ACUSIM with P Projection" ){
      incomp.iprjFlag = 0; incomp.ipresPrjFlag=1;}
    else if ( (string)inp.GetValue("Solver Type") =="ACUSIM" ){
      incomp.iprjFlag = 0; incomp.ipresPrjFlag=0;}
    else if( (string)inp.GetValue("Solver Type") =="ACUSIM with Velocity Projection" ){
      incomp.iprjFlag = 1; incomp.ipresPrjFlag=0;}
    else if( (string)inp.GetValue("Solver Type") =="ACUSIM with Full Projection" ){
      incomp.iprjFlag = 1; incomp.ipresPrjFlag=1;}
    else if( (string)inp.GetValue("Solver Type") =="GMRES Matrix Free"){ 
      inpdat.impl[0] += 10*solflow;}
    else if( (string)inp.GetValue("Solver Type") =="GMRES EBE"){ 
      inpdat.impl[0] += 20*solflow;}
    else if( (string)inp.GetValue("Solver Type") =="svLS"){    // memLS
      inpdat.svLSFlag=1;}
    //else if( (string)inp.GetValue("Solver Type") =="memLS"){    // memLS
    //  inpdat.memLSFlag=1;}
    // MB, added flag options for svLS
    //GMRES sparse is assumed default and has the value of 10, MFG 20,
    // EBE 30


    //    inpdat.niter[0] = inp.GetValue("Number of Solves per Time Step");
    solpar.nGMRES = inp.GetValue("Number of GMRES Sweeps per Solve");
    solpar.Kspace = inp.GetValue("Number of Krylov Vectors per GMRES Sweep");
    inpdat.LHSupd[0] = inp.GetValue("Number of Solves per Left-hand-side Formation");
    inpdat.epstol[0] = inp.GetValue("Tolerance on Momentum Equations");
    inpdat.epstol[6] = inp.GetValue("Tolerance on Continuity Equations"); //MB, uncommented line
    inpdat.epstol[7] = inp.GetValue("Tolerance on svLS NS Solver");  //MB, added svLS tolerance
    incomp.prestol = inp.GetValue("Tolerance on ACUSIM Pressure Projection"); 
    incomp.minIters = inp.GetValue("Minimum Number of Iterations per Nonlinear Iteration");
    incomp.maxIters = inp.GetValue("Maximum Number of Iterations per Nonlinear Iteration");
    inpdat.deltol[0][0]=inp.GetValue("Velocity Delta Ratio"); 
    inpdat.deltol[1][0]=inp.GetValue("Pressure Delta Ratio"); 
    incomp.nPrjs = inp.GetValue("Number of Velocity Projection Vectors");
    incomp.nPresPrjs = inp.GetValue("Number of Pressure Projection Vectors");
    incomp.iverbose = inp.GetValue("ACUSIM Verbosity Level");

    if(pcboiling.solheat==1){ 
      inpdat.epstol[1]=inp.GetValue("Temperature Solver Tolerance");
      inpdat.LHSupd[1]=inp.GetValue("Number of Solves of Temperature per Left-hand-side Formation");
    }

    // The following is where you should put any inputs that are able to 
    // input differently for each scalar.  It is a little tedious in the code 
    // but it should make the solver.inp easier to understand. Note this will 
    // require some care with regression tests.


    if(solscalr>0){
      inpdat.epstol[2]=inp.GetValue("Scalar 1 Solver Tolerance");
      inpdat.LHSupd[2]=inp.GetValue("Number of Solves of Scalar 1 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 1");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[5][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    } 

    if(solscalr>1){
      inpdat.epstol[3]=inp.GetValue("Scalar 2 Solver Tolerance");
      inpdat.LHSupd[3]=inp.GetValue("Number of Solves of Scalar 2 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 2");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[6][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    } 

    if(solscalr>2){
      inpdat.epstol[4]=inp.GetValue("Scalar 3 Solver Tolerance");
      inpdat.LHSupd[4]=inp.GetValue("Number of Solves of Scalar 3 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 3");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[7][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    } 

    if(solscalr>3){
      inpdat.epstol[5]=inp.GetValue("Scalar 4 Solver Tolerance");
      inpdat.LHSupd[5]=inp.GetValue("Number of Solves of Scalar 4 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 4");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[8][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    }
 
    // DISCRETIZATION CONTROL
    
    genpar.ipord = inp.GetValue("Basis Function Order");
    if((string)inp.GetValue("Time Integration Rule") == "First Order")
      inpdat.rhoinf[0] = -1 ;
    else inpdat.rhoinf[0] = (double)inp.GetValue("Time Integration Rho Infinity");
    if((string)inp.GetValue("Predictor at Start of Step")=="Same Velocity")
      genpar.ipred = 1;
    if((string)inp.GetValue("Predictor at Start of Step")=="Zero Acceleration")
      genpar.ipred = 2;
    if((string)inp.GetValue("Predictor at Start of Step")=="Same Acceleration")
      genpar.ipred = 3;
    if((string)inp.GetValue("Predictor at Start of Step")=="Same Delta")
      genpar.ipred = 4;
    
    if((string)inp.GetValue("Weak Form") == "Galerkin")
      solpar.ivart = 1;
    if((string)inp.GetValue("Weak Form") == "SUPG")
      solpar.ivart = 2;

    if((string)inp.GetValue("Flow Advection Form") == "Convective")
      solpar.iconvflow = 2;
    else if((string)inp.GetValue("Flow Advection Form") == "Conservative")
      solpar.iconvflow = 1;
    if((string)inp.GetValue("Scalar Advection Form") == "Convective")
      solpar.iconvsclr = 2;
    else if((string)inp.GetValue("Scalar Advection Form") == "Conservative")
      solpar.iconvsclr = 1;
    if((string)inp.GetValue("Use Conservative Scalar Convection Velocity") == "True")
      sclrs.consrv_sclr_conv_vel = 1;
    else if((string)inp.GetValue("Use Conservative Scalar Convection Velocity") == "False")
      sclrs.consrv_sclr_conv_vel = 0;
    // TAU INPUT 
    if((string)inp.GetValue("Tau Matrix") == "Diagonal-Shakib")
      genpar.itau = 0;
    else  if((string)inp.GetValue("Tau Matrix") == "Diagonal-Franca")
      genpar.itau =1;
    else if((string)inp.GetValue("Tau Matrix") == "Diagonal-Jansen(dev)") 
      genpar.itau = 2;
    else if((string)inp.GetValue("Tau Matrix") == "Diagonal-Compressible")
      genpar.itau = 3;
    else if((string)inp.GetValue("Tau Matrix") == "Matrix-Mallet") 
      genpar.itau = 10;
    else if((string)inp.GetValue("Tau Matrix") == "Matrix-Modal")
      genpar.itau = 11;

    genpar.dtsfct = inp.GetValue("Tau Time Constant");
    genpar.dtsfctsclr = inp.GetValue("Tau Time Constant for Scalars");
    genpar.taucfct = inp.GetValue("Tau C Scale Factor");

    // FLOW DISCONTINUITY CAPTURING

      if((string)inp.GetValue("Discontinuity Capturing") == "Off") solpar.iDC = 0;
    else if((string)inp.GetValue("Discontinuity Capturing") == "DC-mallet") solpar.iDC = 1;
    else if((string)inp.GetValue("Discontinuity Capturing") == "DC-quadratic") solpar.iDC = 2;
   else if((string)inp.GetValue("Discontinuity Capturing") == "DC-minimum") solpar.iDC = 3;    
    else {
      cout<< "Condition not defined for Discontinuity Capturing \n ";
      exit(1);
    }

    // SCALAR DISCONTINUITY CAPTURING

      vector<int> ivec = inp.GetValue("Scalar Discontinuity Capturing");
      for(i=0; i< 2; i++)  solpar.idcsclr[i] = ivec[i];
      ivec.erase(ivec.begin(),ivec.end());
 

//        if((string)inp.GetValue("Scalar Discontinuity Capturing") == "No") solpar.idcsclr = 0;
//      else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "1") solpar.idcsclr = 1; 
//   else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "2") solpar.idcsclr = 2; 
//   else {
//        cout<< "Condition not defined for Scalar Discontinuity Capturing \n ";
//        exit(1);
//      }     
    if((string)inp.GetValue("Include Viscous Correction in Stabilization") == "True")
      {  
        if(genpar.ipord == 1 ) genpar.idiff = 1;
        else genpar.idiff = 2;
      }
    else { genpar.idiff = 0;}

    timdat.flmpl = inp.GetValue("Lumped Mass Fraction on Left-hand-side");
    timdat.flmpr = inp.GetValue("Lumped Mass Fraction on Right-hand-side");

    timdat.iCFLworst = 0;
    if((string)inp.GetValue("Dump CFL") == "True")
      timdat.iCFLworst = 1;

    intdat.intg[0][0]=inp.GetValue("Quadrature Rule on Interior");
    intdat.intg[0][1]=inp.GetValue("Quadrature Rule on Boundary");
    genpar.ibksiz = inp.GetValue("Number of Elements Per Block");

    ((string)inp.GetValue("Turn Off Source Terms for Scalars") 
     == "True")? sclrs.nosource=1:sclrs.nosource=0;
    sclrs.tdecay=inp.GetValue("Decay Multiplier for Scalars");

    // TURBULENCE MODELING PARAMETER
    int tpturb = turbvari.iles-turbvari.irans;
    int ifrule;
    
    if( tpturb != 0 ){


      turbvari.nohomog =inp.GetValue("Number of Homogenous Directions");

      if((string)inp.GetValue("Turbulence Wall Model Type") == "Slip Velocity") turbvar.itwmod = 1;
      else if((string)inp.GetValue("Turbulence Wall Model Type") == "Effective Viscosity") turbvar.itwmod = 2; 
      else  turbvar.itwmod = 0;
      if (turbvari.irans < 0) turbvar.itwmod = turbvar.itwmod*(-1);
      ifrule  = inp.GetValue("Velocity Averaging Steps");
      turbvar.wtavei =(ifrule >0)? 1.0/ifrule : -1.0/ifrule;
 
      if(turbvari.iles == 1){
        
        if((string)inp.GetValue("Dynamic Model Type") == "Bardina") turbvari.iles += 10;
        else if((string)inp.GetValue("Dynamic Model Type") == "Projection") turbvari.iles += 20;

        ifrule = inp.GetValue("Filter Integration Rule");
        turbvari.iles += ifrule-1;
        ifrule = inp.GetValue("Dynamic Model Averaging Steps");
        turbvar.dtavei = (ifrule >0)? 1.0/ifrule : -1.0/ifrule;
        turbvar.fwr1 = inp.GetValue("Filter Width Ratio");
        turbvar.flump = inp.GetValue("Lumping Factor for Filter");


        if ((string)inp.GetValue("Model Statistics") == "True" ) {
          turbvari.modlstats = 1; } 
        else {
          turbvari.modlstats = 0; }   
 
        if ((string)inp.GetValue("Double Filter") == "True" ) {
          turbvari.i2filt = 1; } 
        else {
          turbvari.i2filt = 0; }  

        if ((string)inp.GetValue("Model/SUPG Dissipation") == "True" ) {
          turbvari.idis = 1; } 
        else {
          turbvari.idis = 0; }


        if((string)inp.GetValue("Dynamic Model Type") == "Standard") {

          if((string)inp.GetValue("Dynamic Sub-Model Type") == "None") 
            turbvari.isubmod = 0;
          else if((string)inp.GetValue("Dynamic Sub-Model Type") =="DFWR")
            turbvari.isubmod = 1;
          else if((string)inp.GetValue("Dynamic Sub-Model Type") =="SUPG") 
            turbvari.isubmod = 2;
        }
        else if((string)inp.GetValue("Dynamic Model Type") == "Projection") {

          if((string)inp.GetValue("Projection Filter Type") == "Linear")
            turbvari.ifproj = 0;
          else if((string)inp.GetValue("Projection Filter Type") =="Quadratic") 
            turbvari.ifproj = 1;          

          if((string)inp.GetValue("Dynamic Sub-Model Type") == "None")
            turbvari.isubmod = 0;
          else if((string)inp.GetValue("Dynamic Sub-Model Type") =="ConsistentProj") 
            turbvari.isubmod = 1;          
        }

      }
//
// could be turbulence wall model in DNS (i.e. no turb. model)
//
// will assume we are using RANS models here
//
    } else if (turbvari.idns > 0) {
      if((string)inp.GetValue("Turbulence Wall Model Type") == "Slip Velocity") turbvar.itwmod = -1;
      else if((string)inp.GetValue("Turbulence Wall Model Type") == "Effective Viscosity") turbvar.itwmod = -2;
      else  turbvar.itwmod = 0;
    }
  
    // SPEBC MODELING PARAMETERS

    if ( (spebcvr.irscale = inp.GetValue("SPEBC Model Active")) >= 0 ){

      ifrule  = inp.GetValue("Velocity Averaging Steps");
      turbvar.wtavei =(ifrule >0)? 1.0/ifrule : 1.0/inpdat.nstep[0];
      spebcvr.intpres = inp.GetValue("Interpolate Pressure");
      spebcvr.plandist = inp.GetValue("Distance between Planes");
      spebcvr.thetag  = inp.GetValue("Theta Angle of Arc");
      spebcvr.ds = inp.GetValue("Distance for Velocity Averaging");
      spebcvr.tolerence = inp.GetValue("SPEBC Cylindrical Tolerance");
      spebcvr.radcyl = inp.GetValue("Radius of recycle plane");
      spebcvr.rbltin  = inp.GetValue("Inlet Boundary Layer Thickness");
      spebcvr.rvscal  = inp.GetValue("Vertical Velocity Scale Factor");
    } 
        // Contact Angle Control
     if ( (string)inp.GetValue("Contact Angle Control") == "True")
        contactangle.CA_flag = 1.0; else contactangle.CA_flag = 0.0;
    
     if (contactangle.CA_flag == 1.0){
        contactangle.theta_adv = inp.GetValue("Advancing Contact Angle");
        contactangle.theta_rec = inp.GetValue("Receding Contact Angle");
        contactangle.Forcecont = inp.GetValue("Constant Multiplied In The Force");
	contactangle.stretch = inp.GetValue("Constant To Adjust The Tangent Value");
	contactangle.Fapp_thick = inp.GetValue("Relative Thickness Of The Force Application Region");
	contactangle.Fapp_heigh = inp.GetValue("Relative Height Of The Force Application Region");
	}
      
    // CARDIOVASCULAR MODELING PARAMETERS
    if ( (string)inp.GetValue("Time Varying Boundary Conditions From File") == "True") 
      nomodule.itvn = 1; else nomodule.itvn = 0;
    if ( nomodule.itvn ==1)
      nomodule.bcttimescale = inp.GetValue("BCT Time Scale Factor");

    nomodule.ipvsq=0;
    if(nomodule.icardio = inp.GetValue("Number of Coupled Surfaces")){
      if ( nomodule.icardio > MAXSURF ) {
        cout << "Number of Coupled Surfaces > MAXSURF \n";
        exit(1);
      } 
      if ( (string)inp.GetValue("Pressure Coupling") == "None") 
        nomodule.ipvsq=0;
      if ( (string)inp.GetValue("Pressure Coupling") == "Explicit") 
        nomodule.ipvsq=1;
      if ( (string)inp.GetValue("Pressure Coupling") == "Implicit") 
        nomodule.ipvsq=2;
      if ( (string)inp.GetValue("Pressure Coupling") == "P-Implicit") 
        nomodule.ipvsq=3;

      if(nomodule.numResistSrfs=inp.GetValue("Number of Resistance Surfaces")){
          ivec = inp.GetValue("List of Resistance Surfaces");          
          for(i=0;i<MAXSURF+1; i++) nomodule.nsrflistResist[i] = 0;
          for(i=0; i< nomodule.numResistSrfs; i++){
              nomodule.nsrflistResist[i+1] = ivec[i];
          }
          vec = inp.GetValue("Resistance Values");
          for(i =0; i< MAXSURF+1 ; i++) nomodule.ValueListResist[i] = 0;
          for(i =0; i< nomodule.numResistSrfs ; i++) nomodule.ValueListResist[i+1] = vec[i];
          vec.erase(vec.begin(),vec.end());
      }
      if(nomodule.numImpSrfs=inp.GetValue("Number of Impedance Surfaces")){
          ivec = inp.GetValue("List of Impedance Surfaces");
          for(i=0;i<MAXSURF+1; i++) nomodule.nsrflistImp[i] = 0;
          for(i=0; i< nomodule.numImpSrfs; i++){
              nomodule.nsrflistImp[i+1] = ivec[i];
          }
          if ( (string)inp.GetValue("Impedance From File") == "True")
              nomodule.impfile = 1; else nomodule.impfile = 0;
      }

    }
    if ( (string)inp.GetValue("Time Varying Boundary Conditions From File") == "True") 
      nomodule.itvn = 1; else nomodule.itvn = 0;
    if ( (string)inp.GetValue("Time Varying Boundary Conditions Switch") == "True")
      nomodule.tvbcswitch = 1; else nomodule.tvbcswitch = 0;
    if ( (string)inp.GetValue("Time Varying Boundary Conditions From Routine") == "True") 
      nomodule.itvbc = 1; else nomodule.itvbc = 0;
    if ( (string)inp.GetValue("Convective Pressure Boundary") == "True") 
        nomodule.ibcb_conv_p = 1; else nomodule.ibcb_conv_p = 0;
    if ( (string)inp.GetValue("Convective Pressure Boundary Normalization On") == "True")
        nomodule.ibcb_conv_p_norm = 1; else nomodule.ibcb_conv_p_norm = 0;
    nomodule.ideformwall = 0;
    if((string)inp.GetValue("Deformable Wall")=="True"){
        nomodule.ideformwall = 1;
        nomodule.rhovw = inp.GetValue("Density of Vessel Wall");
        nomodule.thicknessvw = inp.GetValue("Thickness of Vessel Wall");
        nomodule.evw = inp.GetValue("Young Mod of Vessel Wall");
        nomodule.rnuvw = inp.GetValue("Poisson Ratio of Vessel Wall");
        nomodule.rshearconstantvw = inp.GetValue("Shear Constant of Vessel Wall");
        if((string)inp.GetValue("Wall Mass Matrix for LHS") == "True") nomodule.iwallmassfactor = 1;
        else nomodule.iwallmassfactor = 0;
        if((string)inp.GetValue("Wall Stiffness Matrix for LHS") == "True") nomodule.iwallstiffactor = 1;
        else nomodule.iwallstiffactor = 0; 
    } 

   
    // Scaling Parameters Keywords

    outpar.ro = inp.GetValue("Density");
    outpar.vel = inp.GetValue("Velocity");
    outpar.press = inp.GetValue("Pressure");
    outpar.temper = inp.GetValue("Temperature");
    outpar.entrop = inp.GetValue("Entropy");

    // Step Sequencing
 

    ivec = inp.GetValue("Step Construction");
    sequence.seqsize = ivec.size();
    if( sequence.seqsize > 200 || sequence.seqsize < 2 )
     cerr<<"Sequence size must be between 2 and 200 "<<endl;
   
    for(i=0; i< sequence.seqsize; i++)
      sequence.stepseq[i] = ivec[i];
  }
  catch ( exception &e ) {
    cout << endl << "Input exception: " << e.what() << endl << endl;
    ierr = 001;
    print_error_code(ierr);
    return ierr;
  }

  return ierr;
  
}

void print_error_code(int ierr) {
  /*
    Return Error codes:
    0xx         Input error
    1xx         Solution Control
    105         Turbulence Model not supported

    2xx         Material Properties

    3xx         Output Control

    4xx         Discretization Control

    5xx         Scaling Parameters

    6xx         Linear Solver Control
    601         linear solver type not supported
  */
  cout << endl << endl << "Input error detected: " << endl << endl;
  if ( ierr == 001 ) {
    cout << endl << "Input Directive not understood" << endl << endl;
  }
  if ( ierr == 105 ) {
    cout << endl << "Turbulence Model Not Supported" << endl << endl;
  }
  if ( ierr == 601 ) {
    cout << endl << "Linear Solver Type Not Supported" << endl << endl;
  }

}
