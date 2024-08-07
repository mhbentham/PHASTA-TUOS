// Routine contains the structures for reading the user input through
// input_fform.cc. The default values for all these variables are defined in
// input.config.
//
// Input variables that have been previously declared in common.h have to be
// re-declared here, in a consistant structure.   

#include <FCMangle.h>

#define workfc FortranCInterface_GLOBAL_(workfc,WORKFC)
#define fronts FortranCInterface_GLOBAL_(fronts,FRONTS)
#define newdim FortranCInterface_GLOBAL_(newdim,NEWDIM)
#define timer4 FortranCInterface_GLOBAL_(timer4,TIMER4)
#define extrat FortranCInterface_GLOBAL_(extrat,EXTRAT)
#define spongevar FortranCInterface_GLOBAL_(spongevar,SPONGEVAR)
#define turbvar FortranCInterface_GLOBAL_(turbvar,TURBVAR)
#define turbvari FortranCInterface_GLOBAL_(turbvari,TURBVARI)
#define mpistats FortranCInterface_GLOBAL_(mpistats,MPISTATS)
#define spebcvr FortranCInterface_GLOBAL_(spebcvr,SPEBCVR)
#define aerfrc FortranCInterface_GLOBAL_(aerfrc,AERFRC)
#define astore FortranCInterface_GLOBAL_(astore,ASTORE)
#define conpar FortranCInterface_GLOBAL_(conpar,CONPAR)
#define shpdat FortranCInterface_GLOBAL_(shpdat,SHPDAT)
#define datpnt FortranCInterface_GLOBAL_(datpnt,DATPNT)
#define errpar FortranCInterface_GLOBAL_(errpar,ERRPAR)
#define elmpar FortranCInterface_GLOBAL_(elmpar,ELMPAR)
#define genpar FortranCInterface_GLOBAL_(genpar,GENPAR)
#define inpdat FortranCInterface_GLOBAL_(inpdat,INPDAT)
#define intdat FortranCInterface_GLOBAL_(intdat,INTDAT)
#define mio FortranCInterface_GLOBAL_(mio,MIO)
#define mioname FortranCInterface_GLOBAL_(mioname,MIONAME)
#define itrpar FortranCInterface_GLOBAL_(itrpar,ITRPAR)
#define itrpnt FortranCInterface_GLOBAL_(itrpnt,ITRPNT)
#define matdat FortranCInterface_GLOBAL_(matdat,MATDAT)
#define matramp FortranCInterface_GLOBAL_(matramp,MATRAMP)
#define mmatpar FortranCInterface_GLOBAL_(mmatpar,MMATPAR)
#define outpar FortranCInterface_GLOBAL_(outpar,OUTPAR)
#define point FortranCInterface_GLOBAL_(point,POINT)
#define precis FortranCInterface_GLOBAL_(precis,PRECIS)
#define propar FortranCInterface_GLOBAL_(propar,PROPAR)
#define resdat FortranCInterface_GLOBAL_(resdat,RESDAT)
#define solpar FortranCInterface_GLOBAL_(solpar,SOLPAR)
#define timdat FortranCInterface_GLOBAL_(timdat,TIMDAT)
#define timpar FortranCInterface_GLOBAL_(timpar,TIMPAR)
#define incomp FortranCInterface_GLOBAL_(incomp,INCOMP)
#define mtimer1 FortranCInterface_GLOBAL_(mtimer1,MTIMER1)
#define mtimer2 FortranCInterface_GLOBAL_(mtimer2,MTIMER2)
#define timer3 FortranCInterface_GLOBAL_(timer3,TIMER3)
#define title FortranCInterface_GLOBAL_(title,TITLE)
#define sclrs FortranCInterface_GLOBAL_(sclrs,SCLRS)
#define levlset FortranCInterface_GLOBAL_(levlset,LEVLSET)
#define bubstudy FortranCInterface_GLOBAL_(bubstudy,BUBSTUDY)
#define nomodule FortranCInterface_GLOBAL_(nomodule,NOMODULE)
#define sequence FortranCInterface_GLOBAL_(sequence,SEQUENCE)
#define amgvarr FortranCInterface_GLOBAL_(amgvarr,AMGVARR)
#define amgvari FortranCInterface_GLOBAL_(amgvari,AMGVARI)
#define pcboiling FortranCInterface_GLOBAL_(pcboiling,PCBOILING)
#define contactangle FortranCInterface_GLOBAL_(contactangle,CONTACTANGLE)


#define MAXBLK   5000
#define MAXSURF  20  
#define MAXTS   100
#define MAXTOP   5
#define MAXQPT   125
#define MAXSH    125
#define NSD      3
#define NSDSQ    9
#define machin   'RS/6000'
#define machfl   4
#define zero   0.0000000000000000000000000000000d0
#define pt125   0.1250000000000000000000000000000d0
#define pt25   0.2500000000000000000000000000000d0
#define pt33   0.3333333333333333333333333333333d0
#define pt39   0.3968502629920498686879264098181d0
#define pt5   0.5000000000000000000000000000000d0
#define pt57   0.5773502691896257645091487805020d0
#define pt66   0.6666666666666666666666666666667d0
#define pt75   0.7500000000000000000000000000000d0
#define one   1.0000000000000000000000000000000d0
#define sqt2   1.4142135623730950488016887242097d0
#define onept5   1.5000000000000000000000000000000d0
#define two   2.0000000000000000000000000000000d0
#define three   3.0000000000000000000000000000000d0
#define four   4.0000000000000000000000000000000d0
#define five   5.0000000000000000000000000000000d0
#define pi   3.1415926535897932384626433832795d0

#ifdef __cplusplus
extern "C" {
#endif
  extern struct { 
    int master;
    int numpe;
    int myrank;
  } workfc ;

  extern struct { 
    int maxfront;
    int nlwork;
    int idirstep;
    int idirtrigger;
  } fronts ;

  extern struct { 
    long long int numper;
    long long int nshgt;
    long long int nshg0;
  } newdim ;

  extern struct { 
    double birth;
    double death;
    double comtim;
  } timer4 ;

  extern struct { 
    double ttim[100];
  } extrat ;

  extern struct {
    double zoutsponge, radsponge, zinsponge, grthosponge, grthisponge;
    double betamax;
    int spongecontinuity, spongemomentum1, spongemomentum2;
    int spongeenergy, spongemomentum3;
  } spongevar ;

  extern struct {
    double eles;
    double ylimit[9][3]; /* 9 = 5 + 4 = puvwT + 4Scalars */
    double rmutarget;
    double pzero;
    double wtavei;
    double dtavei;
    double dke;
    double fwr1;
    double flump;
    int ierrcalc;
    int ihessian;
    int itwmod;
    int ngaussf;
    int idim;
    int nlist;
    int nintf[MAXTOP];
  } turbvar ;

  extern struct {
    int irans, iles, idns, isubmod;
    int ifproj;
    int i2filt;
    int modlstats;
    int idis;
    int nohomog;
    int ierrsmooth;
    /* wonder if we should include nintf(MAXTOP) and MAXTOP since its
       in common.h */

/*      int itwmod; */
/*      double rtavei; */
/*      int ierrcalc; */
  } turbvari ;

   extern struct { 
    int iISend;
    int iIRecv;
    int iWaitAll;
    int iAllR;
    int impistat;
    int impistat2;
    double rmpitmr;
    double rISend;
    double rIRecv;
    double rWaitAll;
    double rAllR;
    double rCommu;
  } mpistats ;

 extern struct { 
    int irscale;
    int intpres;
    double plandist;
    double thetag;
    double ds;
    double tolerence;
    double radcyl;
    double rbltin;
    double rvscal;
  } spebcvr ;

  extern struct {
    double scdiff[5];
    double tdecay;
    int nsclr, isclr,nsolt, nosource;
    int consrv_sclr_conv_vel;
  } sclrs;

  extern struct { 
    double flxID[MAXSURF+1][10] ;
    double Force[3];
    double HFlux;
    int nsrflist[MAXSURF+1];
    int isrfIM;
    double flxIDsclr[MAXSURF][4];
  } aerfrc ;

  extern struct { 
    double a[100000];
  } astore ;

  extern struct { 
    int numnp;
    int numel;
    int numelb;
    int numpbc;
    int nen;
    int nfaces;
    int numflx;
    int ndof;
    int iALE;
    int icoord;
    int navier;
    int iblk;
    int irs;
    int iexec;
    int necho;
    int ichem;
    int iRK;
    int nedof;
    int nshg;
    int nnz;
    int istop;
    int nflow;
    int nnz_tot;
    int idtn;
  } conpar ;
  
  extern struct { 
    double epsilon_ls;
    double epsilon_lsd;
    double dtlset;
    double dtlset_cfl;
    double redist_toler;
    double redist_toler_curr;
    double r_int_buffer;
    double r_int_elem_size;
    double phvol[2];
    double AdjRedistVelCFL;
    double BubRad;
    double vf_target;
    double C_int_adjust;
    double vf_now;
    double vf_obj;
    double vfcontrcoeff;
    double C_int_cap;
    double epsilonBT;

    double coalbubrad;
    double coalcon_dist;    
    double avgxcoordold[100];
    double avgycoordold[100];
    double avgzcoordold[100];

    int iLSet;
    int iuse_vfcont_cap;
    int i_num_bubbles;
    int ivconstraint;
    int iSolvLSSclr1;
    int iSolvLSSclr2;
    int i_redist_loop_flag;
    int i_redist_max_iter;
    int i_spat_var_eps_flag;
    int i_dtlset_cfl;
    int i_check_prox;
    int i_gradphi;
    int i_focusredist;
    int i_AdjRedistVel;
    int iBT;
    int id2w;
    int icoalCtrl;
    int icoalcon_verbo;
    int coalcon;
    int update_coalcon;
    int coaltimtrak;
    int coalest;
    int coalcon_rem[100];
  } levlset;

  extern struct {
    double DomainSize[6];
    double phi_inner;
    double phi_outer;
    int iClrLiq;
    int iBK;
    int Nbubtot;
    int Nghost;
  } bubstudy;

  extern struct { 
    int nshape;
    int nshapeb;
    int maxshb;
    int nshl;
    int nshlb;
    int nfath;
    int ntopsh;
    int nsonmax;
  } shpdat ;

  extern struct { 
    int mshp;
    int mshgl;
    int mwght;
    int mshpb;
    int mshglb;
    int mwghtb;
    int mmut;
    int mrhot;
    int mxst;
  } datpnt ;

  extern struct { 
    int lelCat;
    int lcsyst;
    int iorder;
    int nenb;
    int nelblk;
    int nelblb;
    int ndofl;
    int nsymdl;
    int nenl;
    int nfacel;
    int nenbl;
    int intind;
    int mattyp;
  } elmpar ;

  extern struct {
    int numerr;
  } errpar ;

  extern struct { 
    double E3nsd;
    int I3nsd;
    int nsymdf;:
    int ndofBC;
    int ndiBCB;
    int ndBCB;
    int Jactyp;
    int jump;
    int ires;
    int iprec;
    int iprev;
    int ibound;
    int idiff;
    int lhs;
    int itau;
    int ipord;
    int ipred;
    int lstres;
    int iepstm;
    double dtsfct;
    double dtsfctsclr;
    double taucfct;
    int ibksiz;
    int iabc;
    int isurf;
    int idflx;
    double Bo;
    double CoalInvSigma;
    double presavg;
    int EntropyPressure;
  } genpar ;

  extern struct { 
    // MAGNUS, changed epstol[6] to epstol[8]
    double epstol[8];  /* 1+ max number of scalars  (beginning of the
                          end of time sequences) */
    double Delt[MAXTS];
    double CFLfl[MAXTS];
    double CFLsl[MAXTS];
    int nstep[MAXTS];
    int niter[MAXTS];
    int impl[MAXTS];
    double rhoinf[MAXTS];
    int LHSupd[6];
    int loctim[MAXTS];
    double deltol[2][MAXTS];
    double CFLfl_max;
    int iCFLfl_maxelem;
    int iflag_cfl_dt;
    double CFLfl_limit;
    double timestart; 
    double CFLls_max;
    int iCFLls_maxelem;
    int svLSFlag; //MAGNUS, added the svLSFlag
  } inpdat ;

  extern struct { 
    int iin;
    int igeom;
    int ipar;
    int ibndc;
    int imat;
    int iecho;
    int iout;
    int ichmou;
    int irstin;
    int irstou;
    int ihist;
    int iflux;
    int ierror;
    int itable;
    int iforce;
    int igraph;
    int itime;
    int ivol;
    int istat;
    int ivhist;
  } mio ;

  extern struct { 
    double fin;
    double fgeom;
    double fpar;
    double fbndc;
    double fmat;
    double fecho;
    double frstin;
    double frstou;
    double fhist;
    double ferror;
    double ftable;
    double fforce;
    double fgraph;
    double ftime;
    double fvol;
    double fstat;
    double fvhist;
  } mioname ;

  extern struct { 
    double eGMRES;
    int lGMRES;
    int iKs;
    int ntotGM;
  } itrpar ;

  extern struct { 
    int mHBrg;
    int meBrg;
    int myBrg;
    int mRcos;
    int mRsin;
  } itrpnt ;

  extern struct { 
    double datmat[MAXTS][7][3];
    int matflg[MAXTS][6];
    int nummat;
    int mexist;
  } matdat ;

  extern struct {
    double tmu ;
    double trho ;
    double omu ;
    double orho ;
    double qrts0 ;
    double qrts1 ;
    int iramp ;
    int nrts0 ;
    int nrts1 ;
  } matramp ;

  extern struct { 
    double pr,Texp, Planck, Stephan, Nh, Rh, Rgas;
    double gamma, gamma1, s0;
    //, const, xN2, xO2;
    //double yN2,    yO2,    Msh[5], cpsh[5],s0sh[5],h0sh[5];
    //double Rs[5],  cps[5], cvs[5], h0s[5], Trot[5],sigs[5];
    //double Tvib[5],g0s[5], dofs[5],ithm;
  } mmatpar ;

  extern struct { 
    double ro;
    double vel;
    double temper;
    double press;
    double entrop;
    int ntout;
    int ioform;
    int iowflux;
    int iofieldv;
    char iotype[80];
    int ioybar;
    int ivort;
    int icomputevort;
    int nsynciofiles;
    int nsynciofieldswriterestart;
    int iMeshingTool;
    /*  int iostats; */
/*      int ipresref; */
  } outpar ;

  extern struct { 
    int mbeg;
    int mend;
    int mprec;
  } point ;

  extern struct { 
    double epsM;
    int iabres;
  } precis ;

  extern struct { 
    int npro;
  } propar ;

  extern struct { 
    double resfrt;
  } resdat ;

  extern struct { 
    int imap;
    int ivart;
    int iDC;
    int iPcond;
    int Kspace;
    int nGMRES;
    int iconvflow;
    int iconvsclr;
    int idcsclr[2];
  } solpar ;

  extern struct { 
    double time;
    double CFLfld;
    double CFLsld;
    double Dtgl;
    double Dtmax;
    double alpha;
    double etol;
    int lstep;
    int ifunc;
    int itseq;
    int istep;
    int iter;
    int nitr;
    double almi;
    double alfi;
    double gami;
    double flmpl;
    double flmpr;
    double dtol[2];
    int iCFLworst;
  } timdat ;

  extern struct { 
    int LCtime;
    int ntseq;
  } timpar ;

  extern struct { 
    int numeqns[100];
    int minIters;
    int maxIters;
    int iprjFlag;
    int nPrjs;
    int ipresPrjFlag;
    int nPresPrjs;
    double prestol;
    double statsflow[6];
    double statssclr[6];
    int iverbose;
  } incomp ;

  extern struct { 
    double ccode[13];
  } mtimer1 ;

  extern struct { 
    double flops;
    double gbytes;
    double sbytes;
    int iclock;
    int icd;
    int icode;
    int icode2;
    int icode3;
  } mtimer2 ;

  extern struct { 
    double cpu[11];
    double cpu0[11];
    int nacess[11];
  } timer3 ;

  extern struct { 
    double title;
    int ititle;
  } title ;

  extern struct {
    int intg[MAXTS][2];
  }intdat;

  extern struct {
    double bcttimescale;    
    double ValueListResist[MAXSURF+1];
    double rhovw;
    double thicknessvw;
    double evw;
    double rnuvw;
    double rshearconstantvw;
    double betai;
    int icardio;
    int itvn;
    int ipvsq;
    int numResistSrfs;
    int nsrflistResist[MAXSURF+1];
    int numImpSrfs;
    int nsrflistImp[MAXSURF+1];
    int impfile;
    int ideformwall;  
    int iwallmassfactor;
    int iwallstiffactor;      
    int tvbcswitch;
    int ibcb_conv_p;
    int ibcb_conv_p_norm;
    int itvbc;
 } nomodule;

  extern struct {
    double bubboil;
    double solheat;
    double bubgrow;
    double h_fg;
    double T_sat;
    double epsilon_lst;
    double delt_T[4];
    double numshell[2][4];
    double bubble_vol[4];
    double bubble_tempG[4];
    double bubbleID[4];
   } pcboiling;

  extern struct {
    double CA_flag;
    double theta_adv;
    double theta_rec;
    double Forcecont;
    double stretch;
    double Fapp_thick;
    double Fapp_heigh;
  } contactangle;

  extern struct {
    int seqsize;
    int stepseq[100];
  } sequence;

  extern struct {
    double strong_eps;      /* strong criterion Stuben factor    */
    double ramg_eps;        /* AMG convergence eps               */
    double ramg_relax;       /* relaxation factor Gauss-Seidel/Jac*/
    double ramg_trunc;      /* truncation select */
 } amgvarr ;
  
  extern struct {
    int irun_amg;           /* Employ AMG feature solfar.f      */
    int irun_amg_sa;        /* Run AMG stand alone ,tamg or SAMG */
    int irun_amg_prec;      /* Run AMG as preconditioner to CG */
    int iamg_verb;          /* amg verbosity flag                */
    int iamg_neg_sten;      /* neg only stencil or neg and pos   */
    int iamg_nlevel;        /* number of levels 2-V etc.         */
    int iamg_c_solver;     /* solve fine level iter. method     */
    int iamg_prescale;       /* diagonal scale AMG             */
    int iamg_init;           /* setup flag */
    int iamg_ppe_frez;       /* how many solfars to re extract ppe */
    int iamg_setup_frez;    /* how many solfars to re setup amg */
    int iamg_iai_frez;      /* how many solfars to re iai amg */
    int iamg_interp;        /* interpolation select */
    int maxnev;             /* total eigenvectors used for ggb*/
    int maxncv;             /* total iterative vectors for ggb*/
    int mlsdeg;             /* Polynomial Smoothing (MLS) degree */
    int iamg_scale;          /* control the scaling of PPE */
    int iamg_reduce;        /* Run a reduced case */
 } amgvari ;

#ifdef __cplusplus
}
#endif
