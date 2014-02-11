#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __channelpot_h
#define __channelpot_h

typedef struct {
  double Swid;
  int Nlev;
  int Nb;
  int Nt;
  int pmax;
  double H;
  double sigma;
  double *z0;
  double *d0;
  int *levellist;
  int *ixptr;
  int *jxptr;
} channelparam;


void setchannelparams(channelparam* cparams,double width,double tol);
void channelsort(void *This0,channelparam* cparams,const double* zr,
                 const double* er,int* ix,int *jx,int* ixptr,int* jxptr,
                 double xmin, int N, int NE);
void cleanupchannelparam(channelparam* cparams);
void replicatezd(void *This0,channelparam* cparams);
void replicatepositions(void *This0,channelparam* cparams, void *tmpvar0,int N);
void replicatecoeffs(void *This0,channelparam* cparams);
void MPexp_connect_CPU_channelpot_(void *This0,void *C0,channelparam* cparams,
                                   int *levm2p,int *maxm2p,double cutoff);
void addmirrorcorrection(void *This0,channelparam* cparams);
void streaminteraction(void *This0,channelparam* cparams,double *pr,double* pi,double* qr,
                       double* qi,const double* zr,const double* zi,
                       const double* er,const double* ei,const double* mr,
                       const double* mi,int N,int NE,const void *panels0);
#endif
