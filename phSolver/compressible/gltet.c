int GaussLegendreTet(int,int,int,double [][4],double[]);

#ifdef sun4_5
gltet_(int *n1, int *n2,int *n3, double pt[][4], double wt[], int *npt)
#endif
#ifdef ibm6000
gltet(int *n1, int *n2,int *n3, double pt[][4], double wt[], int *npt)
#endif
#ifdef sgi
void gltet_(int *n1, int *n2,int *n3, double pt[][4], double wt[], int *npt)
#endif
#ifdef decalp
void gltet_(int *n1, int *n2,int *n3, double pt[][4], double wt[], int *npt)
#endif
#if defined(intel) || defined(__GNUC__)
void GLTET(int *n1, int *n2,int *n3, double pt[][4], double wt[], int *npt)
#endif
{ 
  *npt = GaussLegendreTet(*n1,*n2,*n3,pt,wt);
}
