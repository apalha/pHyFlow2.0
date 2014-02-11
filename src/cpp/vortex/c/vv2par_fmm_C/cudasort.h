/* cudasort.h */

/* S. Engblom and A. Goude 2011-10-21 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __cudasort_h
#define __cudasort_h

#include "cudaeval.h"

#define maxblockcount 4096*2
#define threadsperblock 256

// always use integer power of two for the number of threads
#define ZDMAXTHREADS 256
#define PREPARESPLITTHREADCOUNT 128
#define PREPARESPLIT32THREADCOUNT 32
#define SINGLETHREADTHREADCOUNT 128
#define SINGLETHREADMAXTHREADS 16384
#define SPLITSHIFT 4096
#define MAXBLOCKSZDMAXMULTI 16384
#define MAXCONNECTIVITYTHREADS 256
#define MAXCONNECTIVITYBLOCKS 16384
#define CONVERTTOFLOATMAXTHREADS 256
#define CONVERTTOFLOATMAXBLOCKS 16384
#define CUMSUMSHIFTSTEP 8
#define CUMSUMTHREADS 256
#define CUMSUMMAXBLOCKS 16384

#define MEDIAN_OF_32

// work-around for Visual Studio; for gcc, the definitions in fmm.h work fine
#define meshratio FMM_MESHRATIO
#define theta FMM_THETA

#define INDEXCACHELENGTH 2

#ifdef SORTLIMIT
#define SORTLIMITSTRING ,SORT_REAL *leftlimitvalues,SORT_REAL *rightlimitvalues,float distlimit,float disttarget
#define SORTLIMITCALLINGSTRING ,leftlimitvalues,rightlimitvalues,distlimit,disttarget
#define SORTLIMITSTRING2 ,SORT_REAL *leftlimitvalues,SORT_REAL *rightlimitvalues,float sortlimit
#define SORTLIMITCALLINGSTRING2 ,leftlimitvalues,rightlimitvalues,sortlimit
#else
#define SORTLIMITSTRING
#define SORTLIMITCALLINGSTRING
#define SORTLIMITSTRING2
#define SORTLIMITCALLINGSTRING2
#endif


//VALIDATEPARTITIONING checks if the partitioning was correct for zr and zi. Will slow down the code significantly
// #define VALIDATEPARTITIONING
//CHECKPARTITIONING debugs multiblockpartition. Will take some time
// #define CHECKPARTITIONING
// #define DEBUGSINGLEBLOCKPARTITION

#ifdef VALIDATEPARTITIONING
#define VALIDATEPARTITIONINGSTRING1 ,const double* zr,const double* zi
#define VALIDATEPARTITIONINGSTRING2 ,zr,zi
#else
#define VALIDATEPARTITIONINGSTRING1
#define VALIDATEPARTITIONINGSTRING2
#endif

void CUDA_perform_partitioning(cudavariables *GPUvars,
                               int N,int NE,int nlevels VALIDATEPARTITIONINGSTRING1);
void CUDA_copy_vectors(cudavariables *GPUvars,
                               int N,int NE,
                               const double* zr,const double *zi,
                               const double *er,const double *ei);

// wrapper functions to allow for C++-compilation:
void calcdabs(const void *d,SORT_REAL *dabs,int count);
void cudaCreateConnectivity(int *jcptr,int *kcptr,int *ir,
                            int *oldjcptr,int *oldkcptr,int *oldir,
                            int count,int maxm2p,void *z,
                            SORT_REAL *dabs,SORT_REAL cutoff,
                            int lastlevel,int *outputvector
                            DEBUGVECTORSTRING);
void cumsumlist(int* oldjcptr,int* oldkcptr,int* jcptr,size_t count,cudavariables* GPUvars,int evalshift);

// corresponding CUDA functions
#ifdef __CUDACC__
__global__ void cudacreateconnectivity(int *jcptr,int *kcptr,int *ir,
                                       int *oldjcptr,int *oldkcptr,int *oldir,
                                       int count,int maxm2p,
                                       SORT_DCMPLX *z,
                                       SORT_REAL *dabs,SORT_REAL cutoff,
                                       int lastlevel,
                                       int *outputvector DEBUGVECTORSTRING);
__global__ void calculatedabs(const SORT_DCMPLX *d,SORT_REAL *dabs,int count);
#endif /* __CUDACC__ */

#endif /* __cudasort_h */
