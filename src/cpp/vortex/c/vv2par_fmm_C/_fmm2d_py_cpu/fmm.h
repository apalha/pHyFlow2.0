/* fmm.h  -- Fast multipole method for 2-D potentials in free space */

/* S. Engblom and A. Goude 2011-10-19 (Revision, platforms+CUDA) */
/* S. Engblom 2009-10-08 (Revision) */
/* S. Engblom 2009-04-15 (Major revision, adaptivity) */
/* S. Engblom 2007-07-08 (Major revision, performance) */
/* S. Engblom 2007-06-19 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __fmm_h
#define __fmm_h

#include "stdlib.h" //for VS compatibility, restrict has to be defined here
#include <math.h>
#include <float.h>

#ifdef __cplusplus
#define restrict __restrict
#endif

#ifdef _MSC_VER
#ifndef gamma_
#define gamma_ 0.5772156649015328606065121
#endif
#ifndef VC16ALIGN
#define VC16ALIGN __declspec(align(16))
#endif
#ifndef GCC16ALIGN
#define GCC16ALIGN
#endif
#else
const static double gamma_ = 0.5772156649015328606065121;
#ifndef VC16ALIGN
#define VC16ALIGN
#endif
#ifndef GCC16ALIGN
#define GCC16ALIGN __attribute__((aligned(16)))
#endif
#endif

#ifndef M_PI
#define M_PI 3.1415926535897932384626434
#endif

#ifndef C_CODE          /*if C_CODE is not defined, we are in MEX file mode*/
#include "mex.h"
#else                   /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mexPrintf printf
#define mxLogical int
#define mxArray double
#endif

#include "channelpot.h"

typedef enum {SCULLY = 0,RANKINE = 1,OSEEN = 2,DIRAC = 3,
              ERR_SMOOTHER = -1} SMOOTHER;

// handle support for CUDA
#ifdef CUDASUPPORT
#undef CHANNELPOT
// necessary requirements, do not change unless you know what you are
// doing:
#ifndef CUDATIME
#define CUDATIME
#endif

#ifndef INLINECOMPLEX
#define INLINECOMPLEX
#endif

#ifndef INLINECOMPLEX_DIRECT
#define INLINECOMPLEX_DIRECT
#endif

#if defined(CUDASORT) && !defined(REORDERCUDA)
#define REORDERCUDA
#endif

#if defined(MULTIPOLESHIFTCUDA) && !defined(MULTIPOLEINITCUDA)
#define MULTIPOLEINITCUDA
#endif

#ifdef CUDASORT
#undef PANELSORT
#endif

#else /* CUDASUPPORT */

#undef MULTIPOLEEVALCUDA
#undef MULTIPOLEINITCUDA
#undef MULTIPOLESHIFTCUDA
#undef REORDERCUDA
#undef CUDASORT

#endif /* CUDASUPPORT */

// #define CHECKM2PS // check coefficients in m2ps
// #define CHECKP2P // check coefficients in p2p (will give errors if m2ps fails as well as initmp)
// #define CHECKINITCUDA // check coefficients in initialization (init & m2m)

#define CUDABLOCKSPERMP 8 //blocks per multiprocessor. Used for determining number of threads in some algorithms

#if !defined(CUDASUPPORT) ||                                            \
  (defined(VALIDATECUDACONNECTIVITY) ||                                 \
   defined(CHECKINITCUDA) ||                                            \
   defined(CHECKM2PS) ||                                                \
   defined(CHECKP2P) ||                                                 \
   !defined(MULTIPOLEEVALCUDA) ||                                       \
   !defined(MULTIPOLEINITCUDA) ||                                       \
   !defined(MULTIPOLESHIFTCUDA))
#undef FULLCUDAFMM //if some part of the algorithm is on the CPU (sorting not included)
#else
#define FULLCUDAFMM //Full FMM algorithm (except sorting) is on the GPU
#endif

#ifdef CUDASUPPORT

/* must avoid external (host) calls on the GPU (device/global) and
   ensure that these operations are inlined */
#if defined(__CUDACC__)
typedef struct { double re,im; } dcmplx;
typedef struct { float re,im; } fcmplx;
#define COMPLEXADD(z,x,y)    {(z).re += x; (z).im += y;}
#define COMPLEXSUB(z,x,y)    {(z).re -= x; (z).im -= y;}
#define COMPLEXASSIGN(z,x,y) {(z).re = x;  (z).im = y;}
#define creal(z) ((z).re)
#define cimag(z) ((z).im)
#else /* __CUDACC__*/
#include <complex>
using namespace std;
typedef complex<double> dcmplx;
typedef complex<float> fcmplx;
#define COMPLEXADD(z,x,y)    {creal(z) += x; cimag(z) += y;}
#define COMPLEXSUB(z,x,y)    {creal(z) -= x; cimag(z) -= y;}
#define COMPLEXASSIGN(z,x,y) ((z) = dcmplx(x,y))
/* "hard" but reasonable assumption that dcmplx is stored as two
   doubles in the order [real imag]: */
#define creal(z) (*(double *)&(z))
#define cimag(z) (*((double *)&(z)+1))
// redirect C99-syntax
#define clog(x) (log(x))
#define cabs(x) (abs(x))
const dcmplx I = dcmplx(0.0,1.0);
#endif /* CUDACC */

// incomplete typedef here, finalized in cudaeval.h
typedef struct cudavariables cudavariables;

#else /* !CUDASUPPORT */

#include <complex>
using namespace std;
typedef complex<double> dcmplx;
// here we can rely on the dcmplx class
#define COMPLEXADD(z,x,y)    ((z) += dcmplx(x,y))
#define COMPLEXSUB(z,x,y)    ((z) -= dcmplx(x,y))
#define COMPLEXASSIGN(z,x,y) ((z) = dcmplx(x,y))
#define creal(z) ((z).real())
#define cimag(z) ((z).imag())
// redirect C99-syntax
#define clog(x) (log(x))
#define cabs(x) (abs(x))
const dcmplx I = dcmplx(0.0,1.0);

// dummy typedef, saves some code
typedef int cudavariables;

#endif /* CUDASUPPORT */

// non-C99 math under VS
#ifdef _MSC_VER
#define snprintf _snprintf
#if !defined(__CUDACC__)
#define isfinite(x) _finite(x)
#define isnan(x) _isnan(x)
inline double fmax(double x,double y)
{
  return x > y ? x : y;
}
inline double exp2(double x)
{
  return exp(x*0.69314718055994533);
}
inline double log2(double x)
{
  return log(x)/0.69314718055994530942;
}
inline double expm1(double x)
{
  if (fabs(x) < 1e-4)
    return x*(1.0+0.5*x*(0.33333333333333333333*x+1.0));
  else
    return exp(x)-1.0;
}
#endif /*_MSC_VER */
#endif /* !__CUDACC__ */

// timing using clock()
#include <time.h>
static clock_t TIME_before, TIME_after;
#define StartTime(printtime,timing) \
  if(printtime || timing)           \
    TIME_before = clock()
#define StopTime(printtime,timing) \
  if(printtime||timing)            \
    TIME_after = clock()
#define PrintTime(s,timing,index,printtime)                             \
  if(printtime)                                                         \
    mexPrintf(s,(double)(TIME_after-TIME_before)/CLOCKS_PER_SEC);       \
  if(timing != NULL)                                                    \
    timing[index] = (double)(TIME_after-TIME_before)/CLOCKS_PER_SEC;

// timing using the GPU
#ifdef CUDATIME
#include "cuda.h"
#include "driver_types.h"
#include "cudatiming.h"
static cudaEvent_t cudaStart,cudaStop;
#define CudaStartTime(printtime,timing)                                 \
  if (printtime || timing)                                              \
    cudaTimingCreateAndStart(&cudaStart,&cudaStop,"cudaStart", "cudaStop");
#define CudaStopTime(printtime,timing)                                  \
  if (printtime || timing)                                              \
    cudaSafeEventRecord(cudaStop,"cudaStop");
#define CudaPrintTime(s,timing,index,printtime)                         \
  if (printtime || timing)                                              \
    cudaTimingSyncPrintAndDestroy(cudaStart,cudaStop,timing,index,      \
                                  printtime,s, "cudaStart", "cudaStop");
#else /* !CUDATIME */
#define CudaStartTime(printtime,timing)
#define CudaStopTime(printtime,timing)
#define CudaPrintTime(s,timing,index,printtime)
#endif /* CUDATIME */

/* Criterion for near/far approximation: two boxes interact through
   cluster-cluster approximations provided that the function
   mpexp_theta() below returns true. Generally, small theta (say, less
   than 0.5) means that less coefficients are needed but that the
   interaction lists get larger, and vice versa. */
#ifndef FMM_THETA
#define FMM_THETA 0.5
#endif
const double theta = FMM_THETA; // 0 < theta < 1

/* Criterion for dividing a box. If one side is more than meshratio
   times the other, then the longer side is split into two. Otherwise
   it is split according to an ordering which yields good cache
   effects. */
#ifndef FMM_MESHRATIO
#define FMM_MESHRATIO 2.0
#endif
const double meshratio = FMM_MESHRATIO; // 1.0 <= meshratio

/* Multipole expansion structure. The coefficients are to be
   understood as
  mpexp(z) = coeff1[0]*log(z-z0)+
     coeff1[1]/(z-z0)+coeff1[2]/(z-z0)^2+coeff1[3]/(z-z0)^3+...
  or, depending on the stage of the algorithm,
  mpexp(z) = coeff2[0]+coeff2[1]*(z-z0)+coeff2[2]*(z-z0)^2+... */
typedef struct {
  dcmplx z0,d0;            // midpoint and radius (z0 to upper right corner)
  dcmplx *coeff1,*coeff2;  // multipole coefficients
  int npanel,*panelptr;    // indices of panels in the box
#ifdef RADIALSHRINK
  double absd0;            // distance to the vortex furthest away (<= |d0|)
#endif
} mpexp;

/* Sparse connectivity matrix (special CCS format integer matrix): the
   row indices in ir of column j are found in [jcptr[j],jcptr[j+1]),
   the near connections (think of these as the value 1) in
   [jcptr[j],kcptr[j]) and the far connections (i.e. value -1) in
   [kcptr[j],jcptr[j+1]). The row indices of the near part is sorted
   in ascending- and the far part in descending order,
   respectively. Only the strictly lower triangular part is stored. At
   the finest level a special syntax is used in which the near part is
   grouped into a 'truly' near part (value 1) in [jcptr[j],kcptr[j+N])
   and a 'less strongly' connected part (value 2) in
   [kcptr[j+N],kcptr[j]). In this case the near part is not sorted and
   also, the rowindices are signed (positive for cluster to particle,
   negative for particle to cluster connections). */
typedef struct {
  int *jcptr,*kcptr; // column pointer and pointer to the beginning of
                     // the 'far'-part
  int *ir;           // row indices
} mpSparse;

/* Multipole expansion tree in vectorized form. First comes the root,
   then the 4 children of the root and so on */
typedef struct {
  double xmin,xmax,
    ymin,ymax;       // bounding box
  int pcoeff;        // number of coefficients offset by one; coeff[0]
                     // through coeff[pcoeff]
  int nlevel;        // number of levels offset by one; level =
                     // 0..nlevel
  int *lptr;         // pointer to levels in tree; level l can be
                     // found in root[i] where lptr[l] <= i <
                     // lptr[l+1] and lptr[nlevel+1] is the total
                     // number of boxes
  mpexp *root;       // box 0, the root of the tree
  mpSparse *connect; // connectivity information for each level
  int *ix,*ixptr;    // permutation of potentials
  int *jx,*jxptr;    // permutation of evaluation points
} MPexp;

#include "fmmsort.h"
#include "panel.h"

typedef struct {
  double *zrtmp,*zitmp,*mrtmp,*mitmp,*ertmp,
         *eitmp,*prtmp,*pitmp,*qrtmp,*qitmp,
         *zrnew,*zinew,*mrnew,*minew,*ernew,
         *einew,*prnew,*pinew,*qrnew,*qinew;
} tmpvariables;

// main interface routines
void fmm2dInteract(int N,
                   const double *zr,const double *zi,
                   const double *mr,const double *mi,
                   int NE,
                   const double *er,const double *ei,
                   double *pr,double *pi,
                   double *qr,double *qi,
                   const panel *panels,int Npanel,
                   int pot,double tol,int Ndirect,
                   SMOOTHER smooth,double xopt,double cutoff,bool cont,
                   double* timing,bool printtime,
                   float *sortlimits,double *eta,double channelheight);
void fmm2dGetMesh(int nlhs,mxArray *plhs[],
                  int N,
                  const double *zr,const double *zi,
                  const double *mr,const double *mi,
                  int NE,
                  const double *er,const double *ei,
                  const panel *panels,int Npanel,
                  int pot,double tol,int Ndirect,double cutoff,
                  double *timing,bool printtime,
                  float *sortlimits,double* eta,double channelheight);

// routines at the MPexp-level
void MPexp_setup(MPexp *This,int N,int NE,
                 const double *zr,const double *zi,
                 const double *mr,const double *mi,
                 const double *er,const double *ei,const double *pr,
                 const panel *panels,int Npanel,int ***panelptrlist,
                 int pot, double tol,int Ndirect,bool meshonly,double cutoff,
                 cudavariables *GPUvars,channelparam *cparam,SMOOTHER smooth,double xopt,bool cont,
                 tmpvariables *tmpvar,double* timing,bool printtime);
void MPexp_free(MPexp *This);
void MPexp_eval2(const MPexp *This,
                 double *pr,double *pi,
                 const double *zr,const double *zi,
                 const double *mr,const double *mi,const panel *panels,
                 int pot,SMOOTHER smooth,double xopt,double cutoff,bool cont,channelparam* cparams);
void MPexp_eval(const MPexp *This,
                double *qr,double *qi,
                const double *er,const double *ei,
                const double *zr,const double *zi,
                const double *mr,const double *mi,const panel *panels,
                int pot,SMOOTHER smooth,double xopt,double cutoff,bool cont,channelparam* cparams);

// private helper functions at the MPexp-level
void MPexp_smooth_(int pot,SMOOTHER smooth,double xopt,bool cont,
                   double *cutoff,double *shape,double *scale);
void MPexp_init_CPU_(MPexp *restrict This,int Nf,int Nt,
                     const double *zr,const double *zi,
                     const double *mr,const double *mi,
                     const panel *panels,channelparam* cparams,int pot,
                     double *timing,bool printtime);
void MPexp_m2m_CPU_(MPexp *This,int Nf,int pot,
                    double *timing,bool printtime);
void MPexp_m2ps_CPU_(MPexp *This,int levm2p,
                     double *timing,bool printtime);
void MPexp_p2p_CPU_(MPexp *This,int levm2p,
                    double *timing,bool printtime);

// routines at the mpexp-level
//On the win32 platform, code this code is approximately 15 % faster if inline. No difference on a64
inline void mpexp_init(mpexp *This,int p,
                const double zr,const double zi,
                const double mr,const double mi,
                int pot);
inline void mpexp_init_sse(mpexp *This,int p,
                const double zr,const double zi,const double zr2,const double zi2,
                const double mr,const double mi,const double mr2,const double mi2,
                int pot);
inline void mpexp_initp(mpexp *This,int p,
                 const double zr,const double zi,
                 const double mr,const double mi,
                 int pot);
inline void mpexp_initp_sse(mpexp *This,int p,
                 const double zr,const double zi,const double zr2,const double zi2,
                 const double mr,const double mi,const double mr2,const double mi2,
                 int pot);
bool mpexp_theta(const mpexp *This,const mpexp *that,double cutoff);
int mpexp_theta2(const mpexp *This,const mpexp *that,double cutoff);
void mpexp_eval(const mpexp *This,int p,
                double *vr,double *vi,
                const double *zr,const double *zi,
                const int *ix,int begin,int end);
void mpexp_eval_sse(const mpexp *This,int p,
                double *vr,double *vi,
                const double *zr,const double *zi,
                const int *ix,int begin,int end);
void mpexp_evalmp(const mpexp *This,int p, int pot,
                  double *vr,double *vi,
                  const double *zr,const double *zi,
                  const int *ix,int begin,int end);
void mpexp_evalmp_sse(const mpexp *This,int p, int pot,
                  double *vr,double *vi,
                  const double *zr,const double *zi,
                  const int *ix,int begin,int end);
void mpexp_directInteract2(double *pr,double *pi,
                           const double *zr,const double *zi,
                           const double *mr,const double *mi,
                           const int *ix,
                           int begin1,int end1,int begin2,int end2,
                           int pot,SMOOTHER smooth,
                           double cutoff,double shape,double scale);
void mpexp_directInteract(double *qr,double *qi,
                          const double *er,const double *ei,
                          const int *jx,int begin1,int end1,
                          const double *zr,const double *zi,
                          const double *mr,const double *mi,
                          const int *ix,int begin2,int end2,
                          int pot,SMOOTHER smooth,
                          double cutoff,double shape,double scale);

/*------------------------------------------------------------------------*/
#endif /* __fmm_h */
