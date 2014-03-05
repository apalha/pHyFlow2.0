/* cudaeval.h -- CUDA host code for direct evaluation. */

/* S. Engblom and A. Goude 2011-10-19 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __cudaeval_h
#define __cudaeval_h

#include "fmm.h"
#define INF __longlong_as_double(0x7ff0000000000000)
#define NEGINF __longlong_as_double(0xfff0000000000000)

#if !defined(CUDASORT) || !defined(FULLCUDAFMM) //FLOATSORT currently only implemented for CUDASORT
#undef FLOATSORT
#endif

#ifdef FLOATSORT
#define SORT_REAL float
#define SORT_DCMPLX fcmplx
#else
#define SORT_REAL double
#define SORT_DCMPLX dcmplx
#endif

#ifdef CUDASUPPORT
#include "cuda.h"
#include "driver_types.h"
#include "cudatiming.h"

#ifdef __CUDACC__
static __inline__ __device__ double fetch_double(texture<int2,1> t,int i)
{
  int2 v = tex1Dfetch(t,i);
  return __hiloint2double(v.y,v.x);
}
static __inline__ __device__ void fetch_double2(double2 *dest,texture<int4,1> t,int i)
{
  int4 v = tex1Dfetch(t,i);
  *dest= make_double2(__hiloint2double(v.y,v.x),__hiloint2double(v.w,v.z));
}
#endif
#define TEXTUREFETCH(dest,texref,index) dest=fetch_double(texref,index)
// #define TEXTUREFETCH(dest,texref,index)                   \
//   ((float *)&(dest))[0] = tex1Dfetch(texref,(index)*2);   \
//   ((float *)&(dest))[1] = tex1Dfetch(texref,(index)*2+1);

// #define CUDADEBUGVECTOR // sends an extra vector to the GPU kernel for debug purposes

// uses old slower direct summation instead (compiles faster)
// #define SYNCKERNELOLD
// #define CHECKCUDACALLS //synchronizes after each cuda call and checks for errors

#ifdef CUDADEBUGVECTOR
#define DEBUGVECTORSTRING ,double* debugvector
#define DEBUGVECTORSTRING2 ,GPUvars->debugvector
#define DEBUGVECTORSTRING3 ,debugvector
#else
#define DEBUGVECTORSTRING
#define DEBUGVECTORSTRING2
#define DEBUGVECTORSTRING3
#endif

#ifdef CUDADEBUGCHECKTIME
/* This is to check the time between two different positions in the
   code. Add GPUCHECKSTART at the starting position and GPUCHECKEND at
   the end position. Not used by the standard compilation */

#define GPUCHECKSTART                                                   \
  cudasafe(cudaEventCreate(&((cudavariables*)GPUvars)->checkstart),     \
           "cudaEventCreate checkstart");                               \
  cudasafe(cudaEventCreate(&((cudavariables*)GPUvars)->checkstop),      \
           "cudaEventCreate checkstop");                                \
  cudasafe(cudaEventRecord(((cudavariables*)GPUvars)->checkstart,0),    \
           "cudaEventRecord checkstart");

#define GPUCHECKEND                                                     \
  cudasafe(cudaEventRecord(((cudavariables*)GPUvars)->checkstop),       \
           "cudaEventRecord checkstop");                                \
#else
#define GPUCHECKSTART
#define GPUCHECKEND
#endif

// *** cleanup here...
#define CUDAHORNERM2PS
#define DYNAMICHORNER

#ifdef DYNAMICHORNER
#define CUDAHORNERM2PS
#endif
// #define CUDAM2PSTIMING

#if !defined(CUDASUPPORT) || !defined(MULTIPOLEINITCUDA)
#undef CHECKINITCUDA
#endif

// struct collecting all data to send to the GPU
struct cudavariables {
  int Nf;
  double *er,*ei,*zr,*zi,*mr,*mi,*qr,*qi,*pr,*pi;
  int *ixptr,*jxptr;
  cudaEvent_t start,stop;
  cudaEvent_t mpstart,mpstop;
  cudaEvent_t initstart,initstop;
  cudaEvent_t stepm2mstart,stepm2mstop;
  cudaEvent_t stepp2pstart,stepp2pstop;
  cudaEvent_t stepm2psstart,stepm2psstop;
  cudaEvent_t fulltimestart,fulltimestop;
  cudaEvent_t reorder1start,reorder1stop;
  cudaEvent_t reorder2start,reorder2stop;
  cudaEvent_t sortstart,sortstop;
  cudaEvent_t checkstart,checkstop;
#ifdef CUDADEBUGVECTOR
  double *debugvector;
  int N;
  int NE;
#endif
#ifdef MULTIPOLESHIFTCUDA
  dcmplx *tmpcoeff;
  int *interactionlist;
  int *sortboxstart;
  double *binomial;
  int** jcptr2,**kcptr2,**ir2;
#endif
#if defined(MULTIPOLEEVALCUDA) || defined(MULTIPOLEINITCUDA)

  double *coeff1;
  double *coeff2;
#endif
  SORT_REAL *z0;
  SORT_REAL *d0;
  int *ix,*jx;
  // allocated on host memory, but internal variables on GPU memory
  mpSparse* connect;
  SORT_REAL *dabs;
  cudaDeviceProp info;
  int evalonly; //set if all evaluation is to be performed on eval points
#ifdef SORTLIMIT
  float *sortlimits; //stored on CPU
#endif
  double *eta;
};

//Keeps track of all memory allocations making it easy to spot memory leaks
//Make all allocations through these functions
cudaError_t cudaMallocDebug(void **addr,size_t size);
cudaError_t cudaFreeDebug(void *addr);
void printcudaalloccount(char *message);
int getcudaalloccount();
void resetalloccount();

//Checks the error message of a cura call and aborts the code if a cuda call fails
#define cudasafe(command,message) cudasafe2(command,message,#command,__FILE__,__LINE__)
void cudasafe2(cudaError_t error,const char *message,const char* command,const char *file,int line);
void cudasafeMalloc(void **addr,size_t size);
void checkcudaerror(const char *message);
void cudareset();
void cudastart();

#ifdef CHECKCUDACALLS
#define CHECKCUDAERROR cudasafe(cudaThreadSynchronize(),"check");
void checkcudaerrorwrapper(int line);
#define CHECKCUDAERROR2 checkcudaerrorwrapper(__LINE__);
#else
#define CHECKCUDAERROR
#define CHECKCUDAERROR2
#endif
// CudaMemcpy and CudaMemset with automatic check of the return value. Aborts if call fails
void cudaSafeMemcpy(void *dst,const void *src,size_t count,
                    enum cudaMemcpyKind kind,const char *str);
void cudaSafeMemset(void *dst,int value,size_t count,const char* str);
cudaError_t cudaGetLastErrorDummy();

// Prototypes for functions called by the GPU initialization of multipole coefficients
void MPexp_init_cuda(const MPexp *This,cudavariables *GPUvars,
                     int pot,int complexpoint);
void cuda_check_init(MPexp *This,int Nt,cudavariables* GPUvars);
void bucketsort_init(int *cumulativelist,int *tmpindices,
                     int *buckets,int *listtosort,int *additionallist,
                     int *outputlist1,int *outputlist2,int count,int bc);
#ifdef __CUDACC__
__host__ __device__
#endif
inline int imin(int x,int y)
{
  return x < y ? x : y;
}
#ifdef  __CUDACC__
__host__ __device__
#endif
inline int imax(int x,int y)
{
  return x < y ? y : x;
}

//Used to sort the interaction lists, to be removed when assymetric lists are working
void quicksort(int *A,int q,int r,int *dummy1,int *dummy2,int *dummy3);
void quicksort2(int *A,int q,int r,int *dummy1);

void MPexp_eval_cuda(const MPexp *This,
		     const double *er,const double *ei,
		     const double *zr,const double *zi,
		     const double *mr,const double *mi,
		     const double *pr,cudavariables *GPUvars,int N,int NE,
		     int pot,SMOOTHER smooth,double xopt,
		     double cutoff,bool cont,double *timing,int printtime);

//collect results from GPU
void MPexp_eval_cuda_result(double *qr,double *qi,double *pr,double *pi,
                            int *ix,int *jx,int N,int NE,int pot,
                            cudavariables *GPUvars,
                            const MPexp *This,
                            double *timing,int printtime);
void getGPUinfo(cudavariables* GPUvars);
#ifdef MULTIPOLEEVALCUDA
void mpexp_eval_cuda(const MPexp *This,cudavariables *GPUvars,
                     const double *pr,int N,int NE,int pot,
                     double *timing,int printtime);
#endif

#ifdef MULTIPOLEINITCUDA
void MPexp_shiftm2m_cuda(const MPexp *This,cudavariables *GPUvars,
                         int pot,double *timing,int printtime);
void MPexp_init_cuda(const MPexp *This,cudavariables *GPUvars,
                     int pot,int complexpoint,double *timing,int printtime);
#endif

#ifdef MULTIPOLESHIFTCUDA
void MPexp_shiftp2p_cuda(const MPexp *This,cudavariables *GPUvars,
                         int levm2p,double *timing,int printtime);
void MPexp_shiftp2p_cuda_scaled(const MPexp *This,cudavariables *GPUvars,
                                int levm2p,double *timing,int printtime);

void MPexp_shiftm2ps_cuda(const MPexp *This,cudavariables *GPUvars,
                          int levm2p,int pot,double *timing,int printtime);
// #define M2PSBINOMIALSIZE 64
#endif /* MULTIPOLESHIFTCUDA */

void direct_eval_cuda(int N,
		      const double *zr,const double *zi,
		      const double *mr,const double *mi,
		      int NE,
		      const double *er,const double *ei,
		      double *pr,double *pi,
		      double *qr,double *qi,
		      const panel *panels,int panelcount,int pot,
		      SMOOTHER smooth,double xopt,double cutoff,
                      bool cont,double* timing,int printtime);

// #ifdef CUDASORT
void cleanupsort(cudavariables* GPUvars,const MPexp *This);
// #endif /* CUDASORT */

#else
//define as empty strings
#define DEBUGVECTORSTRING
#define DEBUGVECTORSTRING2
#define DEBUGVECTORSTRING3
#endif /* CUDASUPPORT */
#endif /* __cudaeval_h */
