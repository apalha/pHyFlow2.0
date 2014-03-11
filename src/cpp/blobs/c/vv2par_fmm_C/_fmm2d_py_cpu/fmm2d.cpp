/* fmm2d.cpp -- Fast multipole method for 2-D potentials in free space
   (C++/CUDA-version) */

/* S. Engblom 2011-06-28 (Panel+Timing I/O) */
/* S. Engblom and A. Goude 2011-04-12 (Major revision) */
/* A. Goude 2010-01-01 (Panels, port to CUDA/GPU) */
/* S. Engblom 2009-05-08 (Revision) */
/* S. Engblom 2007-07-08 (Major revision) */
/* S. Engblom 2006-10-11 (Port to Mex and a major revision) */
/* S. Engblom 2005-01-05 (original code 'fmm2dlp' in C99) */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#else
#include "mex.h"
#include "matrix.h"
#include "expint.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>


typedef enum {NDIRECT = 0,PANEL = 1,ETA = 2,SORTLIMITS = 3,CUTOFF = 4,
              PRINTTIME = 6,CONT = 7,TOL = 9,
              OUT = 10,CHANNELHEIGHT = 11,SMOOTH = 12,POT = 13,XOPT = 14,
              ERR_PROPERTY = -1} PROPERTY;

#include "fmm.h"
#include "panel.h"
#include "mexpanel.h"

#ifndef CUDASUPPORT
#include "direct.h"
#else
#include "cudaeval.h"
#endif
#ifdef CHANNELPOT
#include "directchannelpot.h"
#include "channelpotpanel.h"
#endif

// forward declarations
typedef void (*getvalFun)(void *,const mxArray *);

void getNDIRECT(void *,const mxArray *);
void getCUTOFF(void *,const mxArray *);
void getTOL(void *,const mxArray *);
void getOUT(void *,const mxArray *);
void getPRINTTIME(void *,const mxArray *);
void getCONT(void *,const mxArray *);
void getXOPT(void *,const mxArray *);
void getPOT(void *,const mxArray *);
void getSMOOTH(void *,const mxArray *);
void getSORTLIMITS(void *,const mxArray*);
void getETA(void *,const mxArray*);
void getCHANNELHEIGHT(void *,const mxArray*);

PROPERTY getprop(const char *str);
SMOOTHER getsmooth(const char *str);

#define ISDOUBLEMATRIX(A) (mxIsDouble(A) && !mxIsSparse(A) && \
                           mxGetNumberOfDimensions(A) == 2)
#define ISREALMATRIX(A) (ISDOUBLEMATRIX(A) && !mxIsComplex(A))
#define ISLOGICALMATRIX(A) (mxIsLogical(A) && !mxIsSparse(A) && \
                            mxGetNumberOfDimensions(A) == 2)
#define ISREALSCALAR(A) (mxIsDouble(A) && !mxIsSparse(A) && \
                         !mxIsComplex(A) && mxGetNumberOfElements(A) == 1)
#define ISLOGICALSCALAR(A) (mxIsLogical(A) && \
                            mxGetNumberOfElements(A) == 1)
#ifdef CUDASUPPORT
static int initialized=0;
#endif
/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  // main input
  const double *zr,*zi,*mr,*mi,*er,*ei;
  panel *panels = NULL; // non-const since input is interpreted into struct
  int N,NE = 0,Npanel = 0;

  // values of FMM properties with defaults
  double xopt = 1e-3,cutoff = 2e-3,tol = 1e-6,channelheight=-1;
  int pot = 0,out = 0,Ndirect = 5;
  SMOOTHER smooth = DIRAC;
  bool printtime = false,cont = true;
  float sortlimits[]={2,2,1,1,0,0};
  double eta[]={2,2};

  // output
  double *pr = NULL,*pi = NULL,*qr = NULL,*qi = NULL,*timing = NULL;
  mxArray *Timing = NULL;
  int eval_syntax = 2 < nrhs && ISDOUBLEMATRIX(prhs[2]);
  mxComplexity kind;

  // parsing functions
  static const getvalFun getval[] = {&getNDIRECT,NULL,&getETA,&getSORTLIMITS,getCUTOFF,
                                     NULL,&getPRINTTIME,&getCONT,NULL,&getTOL,
                                     &getOUT,&getCHANNELHEIGHT,&getSMOOTH,&getPOT,&getXOPT};
  void *val[] = {&Ndirect,NULL,eta,sortlimits,&cutoff,
                 NULL,&printtime,&cont,NULL,&tol,
                 &out,&channelheight,&smooth,&pot,&xopt};


  #ifdef CUDASUPPORT
  if(initialized==0) {
    mexAtExit(cudareset);
    initialized=1;
  }
  #endif

  // hidden output: to facilitate testing of expint()
  if (nrhs == 1) {
    double *p = mxGetPr(plhs[0] = mxDuplicateArray(prhs[0]));
    for (int i = 0; i < mxGetNumberOfElements(prhs[0]); i++)
      p[i] = expint(p[i]);
    return;
  }

  // syntax
  if (2 < nlhs || nrhs < 2 || (nrhs-eval_syntax)&1)
    mexErrMsgTxt("Expecting one or two outputs and two or three inputs "
                 "followed by property/value-pairs.");

  // inputs Z and M
  if (!ISDOUBLEMATRIX(prhs[0]) || !ISDOUBLEMATRIX(prhs[1]) ||
      mxGetM(prhs[0]) != mxGetM(prhs[1]) ||
      mxGetN(prhs[0]) != mxGetN(prhs[1]))
    mexErrMsgTxt("First two arguments must be double matrices of "
                 "matching sizes.");
  N = mxGetNumberOfElements(prhs[0]);
  if (N != 0 && !mxIsComplex(prhs[0]))
    mexErrMsgTxt("First argument must be complex.");

  zr = mxGetPr(prhs[0]); zi = mxGetPi(prhs[0]);
  mr = mxGetPr(prhs[1]); mi = mxGetPi(prhs[1]);

  // input ZE -- or just (er,ei) for short
  if (eval_syntax) {
    NE = mxGetNumberOfElements(prhs[2]);
    if (NE != 0 && !mxIsComplex(prhs[2]))
      mexErrMsgTxt("Evaluation points must be complex.");
    er = mxGetPr(prhs[2]); ei = mxGetPi(prhs[2]);
  }

  // remaining checks
  if (nlhs > 1 && !eval_syntax && out != 1)
    mexErrMsgTxt("Three input matrices required for two outputs.");

  if (cutoff <= xopt)
    mexErrMsgTxt("Parameter 'cutoff' must be larger than 'xopt'.");

  if(Npanel && pot != 1)
    mexErrMsgTxt("Panels only supported for 'pot' = 1.");

#ifdef CUDASUPPORT
  if(Npanel>0)
    mexWarnMsgTxt("Panels not supported on GPU version, all panels ignored");
#endif

  // complexity
  kind = pot == 1 || mi != NULL ? mxCOMPLEX : mxREAL;

  // allocate and produce output
  if (out == 1) {
    // special syntax: return FMM-mesh only
    fmm2dGetMesh(nlhs,plhs,
                 N,zr,zi,mr,mi,NE,er,ei,
                 panels,Npanel,pot,tol,Ndirect,cutoff,NULL,
                 printtime,sortlimits,eta,channelheight);
    return;
  }
  else if (nlhs > 1) {
    // evaluate both at potentials Z and at coordinates ZE
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),kind);
    pr = mxGetPr(plhs[0]); pi = mxGetPi(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[2]),mxGetN(prhs[2]),kind);
    qr = mxGetPr(plhs[1]); qi = mxGetPi(plhs[1]);
  }
  else if (!eval_syntax || zr == er && zi == ei && N == NE) {
    // only at potentials Z
    NE = 0;
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]),kind);
    pr = mxGetPr(plhs[0]); pi = mxGetPi(plhs[0]);
  }
  else {
    // at coordinates ZE
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[2]),mxGetN(prhs[2]),kind);
    qr = mxGetPr(plhs[0]); qi = mxGetPi(plhs[0]);
  }

  // special syntax: measure time
  if(out == 2) {
    Timing = mxCreateDoubleMatrix(31,1,mxREAL);
    timing = mxGetPr(Timing);
  }

  // evaluate
  if (tol == 0.0) {
#ifdef CHANNELPOT
    if(channelheight>0)
      directInteractchannelpot(N,zr,zi,mr,mi,NE,er,ei,pr,pi,qr,qi,panels,Npanel,
                   channelheight,smooth,xopt,cutoff,cont,timing,printtime);
    else
#endif
#ifdef CUDASUPPORT
    direct_eval_cuda(N,zr,zi,mr,mi,NE,er,ei,pr,pi,qr,qi,panels,Npanel,
                     pot,smooth,xopt,cutoff,cont,timing,printtime);
#else
    directInteract(N,zr,zi,mr,mi,NE,er,ei,pr,pi,qr,qi,panels,Npanel,
                   pot,smooth,xopt,cutoff,cont,timing,printtime);
#endif /* CUDASUPPORT */
  }
  else
    fmm2dInteract(N,zr,zi,mr,mi,NE,er,ei,pr,pi,qr,qi,panels,Npanel,
                  pot,tol,Ndirect,smooth,xopt,cutoff,cont,timing,
                  printtime,sortlimits,eta,channelheight);

  // throw away the result and return the measured time instead
  if(out == 2) {
    mxDestroyArray(plhs[0]);
    plhs[0] = Timing;
  }

  // deallocate panels
  mxFree(panels);
}
/*------------------------------------------------------------------------*/
PROPERTY getprop(const char *str)
/* Small perfect hash for the supported properties. Note that
   str[0..2] must be readable for this to work. Returns ERR_PROPERTY
   on failure. */
{
  static const char *tab[] = {"ndirect","panel","eta","sortlimits","cutoff",
                              "","printtime","cont","","tol",
                              "out","channelheight","smooth","pot","xopt"};
  const int hash = (((int)str[0] * 3)+((int)str[2] * 2)) % 15;

  return str[0] != '\0' && strcmp(str,tab[hash]) == 0 ?
    (PROPERTY)hash : ERR_PROPERTY;
}
/*------------------------------------------------------------------------*/
void getNDIRECT(void *Ndirect_,const mxArray *rhs)
{
  int *Ndirect = (int *)Ndirect_;

  if (!ISREALSCALAR(rhs))
    mexErrMsgTxt("Expecting a real and scalar option 'ndirect'.");
  *Ndirect = (int)*mxGetPr(rhs);
  if (*Ndirect <= 0)
    mexErrMsgTxt("Property 'ndirect' must be strictly positive.");
}
/*------------------------------------------------------------------------*/
void getCUTOFF(void *cutoff_,const mxArray *rhs)
{
  double *cutoff = (double *)cutoff_;

  if (!ISREALSCALAR(rhs))
      mexErrMsgTxt("Expecting a real and scalar option 'cutoff'.");
  *cutoff = *mxGetPr(rhs);
  if (*cutoff <= 0.0)
    mexErrMsgTxt("Property 'cutoff' must be strictly positive.");
}
/*------------------------------------------------------------------------*/
void getTOL(void *tol_,const mxArray *rhs)
{
  double *tol = (double *)tol_;

  if (!ISREALSCALAR(rhs))
      mexErrMsgTxt("Expecting a real and scalar option 'tol'.");
  *tol = *mxGetPr(rhs);
  if (*tol < 0)
    mexErrMsgTxt("Property 'tol' must not be negative.");
}
/*------------------------------------------------------------------------*/
void getOUT(void *out_,const mxArray *rhs)
{
  int *out = (int *)out_,err = 0;
  char str[10];

  if (!mxIsChar(rhs))
    mexErrMsgTxt("Expecting a character array as option 'out'.");
  if (mxGetString(rhs,str,10) != 0)
    err = 1;
  else if (strcmp(str,"sol") == 0)
    *out = 0;
  else if (strcmp(str,"mesh") == 0)
    *out = 1;
  else if (strcmp(str,"time") == 0)
    *out = 2;
  else
    err = 1;

  if (err)
    mexErrMsgTxt("Property 'out' must be 'sol', 'mesh', or 'time'.");
}
/*------------------------------------------------------------------------*/
void getPRINTTIME(void *printtime_,const mxArray *rhs)
{
  bool *printtime = (bool *)printtime_;

  if (!ISLOGICALSCALAR(rhs))
    mexErrMsgTxt("Expecting a logical scalar option 'printtime'.");
  *printtime = mxIsLogicalScalarTrue(rhs);
}
/*------------------------------------------------------------------------*/
void getCONT(void *cont_,const mxArray *rhs)
{
  bool *cont = (bool *)cont_;

  if (!ISLOGICALSCALAR(rhs))
    mexErrMsgTxt("Expecting a logical scalar option 'cont'.");
  *cont = mxIsLogicalScalarTrue(rhs);
}
/*------------------------------------------------------------------------*/
void getXOPT(void *xopt_,const mxArray *rhs)
{
  double *xopt = (double *)xopt_;

  if (!ISREALSCALAR(rhs))
      mexErrMsgTxt("Expecting a real and scalar option 'xopt'.");
  *xopt = *mxGetPr(rhs);
  if (*xopt <= 0)
    mexErrMsgTxt("Property 'xopt' must be strictly positive.");
}
/*------------------------------------------------------------------------*/
void getPOT(void *pot_,const mxArray *rhs)
{
  int *pot = (int *)pot_;

  if (!ISREALSCALAR(rhs))
    mexErrMsgTxt("Expecting a real and scalar option 'pot'.");
  *pot = (int)*mxGetPr(rhs);
  if (*pot != 0 && *pot != 1)
    mexErrMsgTxt("Property 'pot' must be 0 or 1.");
}
/*------------------------------------------------------------------------*/
void getSMOOTH(void *smooth_,const mxArray *rhs)
{
  SMOOTHER *smooth = (SMOOTHER *)smooth_;
  char str[10];

  if (!mxIsChar(rhs))
    mexErrMsgTxt("Expecting a character array as option 'smooth'.");
  if (mxGetString(rhs,str,10) != 0 ||
      (*smooth = getsmooth(str)) == ERR_SMOOTHER)
    mexErrMsgTxt("Unknown smoother.");
}
/*------------------------------------------------------------------------*/
SMOOTHER getsmooth(const char *str)
/* Small perfect hash for the supported smoothers. Note that str[4]
   must be readable for this to work. Returns ERR_SMOOTHER on
   failure. */
{
  static const char *tab[] = {"scully","rankine","oseen","dirac"};
  const int hash = ((int)str[4])%4;

  return strcmp(str,tab[hash]) == 0 ? (SMOOTHER)hash : ERR_SMOOTHER;
}
/*------------------------------------------------------------------------*/
void getSORTLIMITS(void *sortlimits,const mxArray *rhs)
{
  if(!ISDOUBLEMATRIX(rhs))
    mexErrMsgTxt("Expecting array as sortlimits");
  int Nrhs=mxGetNumberOfElements(rhs);
  if(Nrhs!=4&&Nrhs!=6)
    mexErrMsgTxt("sortlimits must contain 4 or 6 elements");
  double* limits=mxGetPr(rhs);
  float* sortlimits2=(float*)sortlimits;
  sortlimits2[0]=limits[0];
  sortlimits2[1]=limits[1];
  sortlimits2[2]=limits[2];
  sortlimits2[3]=limits[3];
  if(Nrhs==4) {
    sortlimits2[4]=limits[2];
    sortlimits2[5]=limits[3];
  }
  else {
    sortlimits2[4]=limits[4];
    sortlimits2[5]=limits[5];
  }
}
/*------------------------------------------------------------------------*/
void getETA(void *eta,const mxArray *rhs)
{

  if(!ISDOUBLEMATRIX(rhs))
    mexErrMsgTxt("Expecting array as sortlimits");
  int Nrhs=mxGetNumberOfElements(rhs);
  if(Nrhs!=1&&Nrhs!=2)
    mexErrMsgTxt("eta must contain 1 or 2 elements");
  double* limits=mxGetPr(rhs);
  double* etalimits=(double*)eta;
  etalimits[0]=limits[0];
  if(Nrhs==1)
    etalimits[1]=etalimits[0];
  else
    etalimits[1]=limits[1];
}
/*------------------------------------------------------------------------*/
void getCHANNELHEIGHT(void *channelheight_,const mxArray *rhs)
{
  double *channelheight = (double *)channelheight_;

  if (!ISREALSCALAR(rhs))
      mexErrMsgTxt("Expecting a real and scalar option 'channelheight'.");
  *channelheight = *mxGetPr(rhs);
}
/*------------------------------------------------------------------------*/
