/* fmmsort.h */

/* S. Engblom and A. Goude 2011-10-21 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __fmmsort_h
#define __fmmsort_h

#include "fmm.h"
#include "panel.h"

// main routine
void MPexp_sort_(MPexp *This,
                 int N,const double *zr,const double *zi,
                 int NE,const double *er,const double *ei,
                 const panel *panels,int Npanel,int **panelptrlist,
                 int *maxm2p,int *levm2p,int Ndirect,double cutoff,
                 cudavariables *GPUvars,channelparam* cparams,
                 double* timing,bool printtime,bool meshonly);

// CPU/CUDA versions
void MPexp_partition_CPU_(MPexp *This,
                          int N,const double *zr,const double *zi,
                          int NE,const double *er,const double *ei,
                          const panel *panels,int Npanel,
                          int **panelptrlist,
                          int Ndirect,double cutoff,
                          int *ix,int *ixptr,int *ixptr0,
                          int *jx,int *jxptr,int *jxptr0,
                          channelparam* cparams);
void MPexp_connect_CPU_(MPexp *This,mpSparse *C,
                        int *levm2p,int *maxm2p,double cutoff);
void MPexp_partition_post_CUDA_(MPexp *This,cudavariables *GPUvars,
                                int *ix,int *ixptr,int *jx,int *jxptr,
                                int N,int NE,double *timing,bool printtime);
void MPexp_connect_CUDA_(MPexp *This,mpSparse *C,
                         cudavariables *GPUvars,
                         int *levm2p,int* maxm2p,double cutoff,
                         double *timing,bool printtime,bool meshonly);

// CPU helper functions
void MPexp_minmax_(double *xmin,double *xmax,double *ymin,double *ymax,
                   int N,const double *x,const double *y, int* validinput);
void MPexp_box_(MPexp *This,int N,int NE,
                const double *zr,const double *zi,
                const double *er,const double *ei, int* validinput);
void MPexp_partition_(int begin,int *im,int end,double *z0,
                      int *ix,const double *z,int Nmax);
void MPexp_split_(int begin,int *im,int end,double z0,
                  int *ix,const double *z);

#endif /* __fmmsort_h */
