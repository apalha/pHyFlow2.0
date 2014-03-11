/* fmm.cpp  -- Fast multipole method for 2-D potentials in free space */

/* S. Engblom and A. Goude 2011-10-19 (Revision, CUDA-support) */
/* S. Engblom 2009-10-08 (Revision) */
/* S. Engblom 2009-04-15 (Major revision, adaptivity) */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#include <stdio.h>
#include <string.h>

#ifdef CUDASUPPORT
#include "cuda.h"
#include "cudaeval.h"
#endif

#include "fmm.h"
#include "fmmshift.h"
#include "expint.h"
#include "channelpot.h"
#include "channelpotpanel.h"

#if defined(SSE2DIRECT) || defined (SSE2INIT)
#ifdef __SSE3__
#include <pmmintrin.h>
#else
#include <emmintrin.h>
#endif
#endif

/* #defines: (preferably handled via make)
   ---------------------------------------
   INLINECOMPLEX
                 Avoid complex arithmetic in various places in the
                 FMM-algorithm and use inline versions instead. Care
                 has been taken to do this in such a way that the
                 result does not differ from C++'s native complex
                 aritmetic. When in doubt, try #undefining
                 INLINECOMPLEX.

   INLINECOMPLEX_DIRECT
                 Inline complex aritmetic in direct evaluations. On
                 some platforms, it is faster by a factor of up to 8
                 not to rely on certain math-functions and/or complex
                 arithmetics, but to write the complex arithmetic
                 inline. On the one hand this implies possible
                 overflow as the evaluations of complex division and
                 absolute values demands properly scaling; on the
                 other hand this will only be observed with very
                 extreme scaling of the input since the operations are
                 only used in near-field computations. Again, when in
                 doubt, try #undefining INLINECOMPLEX_DIRECT.

   THETACUTOFFCHECK
                 Respect cutoff in theta criterion on CPU (always
                 checked on the GPU).

   It is recommended to #define all these.

   RADIALSHRINK  Use radial distance to the point furthest away in a
                 box instead of the diagonal of the box as a measure
                 of its size.
*/

#ifdef CUDAVELOCITYDEBUG
int NEglobal; // debug purposes
#endif

#ifdef FULLCUDAFMM
const int fullcudafmm=1;
#else
const int fullcudafmm=0;
#endif

#ifdef CUDASORT
const int cudasort=1;
#else
const int cudasort=0;
#endif
void checkboxes(MPexp *This,const double* zr,const double* zi,const double *er,const double *ei,int Nf,int Nt,int NE);
/*------------------------------------------------------------------------*/
void fmm2dInteract(int N,
                   const double *zr,const double *zi,
                   const double *mr,const double *mi,
                   int NE,
                   const double *er,const double *ei,
                   double *pr,double *pi,
                   double *qr,double *qi,
                   const panel *panels,int Npanel,
                   int pot,double tol,int Ndirect,
                   SMOOTHER smooth,double xopt,double cutoff,
                   bool cont,
                   double *timing,bool printtime,
                   float *sortlimits,double* eta,double channelheight)
/* Evaluation by the fast multipole method to accuracy TOL of the
   potential from N pointmasses (mr,mi) (with mi possibly NULL) at the
   positions (zr,zi), and Npanel panels. The evaluation is done at the
   points (zr,zi) if (pr,pi) is not NULL and/or at the points (er,ei)
   if (qr,qi) is not NULL (it is required that NE is zero if (qr,qi)
   are NULL). The result is added to the vectors (pr,pi) and (qr,qi),
   respectively, which must be allocated prior to call. Logarithmic
   potential for pot == 0, harmonic potential otherwise. The smoothing
   is done using smoother smooth with parameters
   (xopt,cutoff,cont). */
{
#ifdef CUDAVELOCITYDEBUG
  NEglobal = NE;
#endif

  // early return: cannot allocate anything
#ifdef PANELSORT
  if (N+Npanel <= 0) return;
#else
  if (N <= 0) return;
#endif

  // multipole structure
  MPexp MP;
  int **panelptr = NULL; /* Has do be declared here to be deallocated
                            at the end. Allocation managed by
                            MPexp_setup(). */

  cudavariables GPUvars; // note: dummy int #ifndef CUDASUPPORT (see fmm.h)
  channelparam cparams;
  int multiplier=1;
#ifdef CHANNELPOT
  cparams.H=channelheight;
  if(channelheight>0)
    multiplier=3;
#endif
#ifdef CUDASUPPORT
  /* Create the structure GPUvars, which contains pointers to all
     global variables on the GPU. Set all pointers to zero since they
     may not always be allocated in all cases. */
  memset(&GPUvars,0,sizeof(GPUvars));
  //check that no errors exist from prevoius cuda functions
  cudastart();
  checkcudaerror("Check before first cuda call");
  getGPUinfo(&GPUvars);
  if (timing || printtime)
    cudaTimingCreateAndStart(&GPUvars.sortstart,&GPUvars.sortstop,
                             "sortstart","sortstop");
  if(pr==NULL&&NE>0)
    GPUvars.evalonly=1;
#ifdef SORTLIMIT
  GPUvars.sortlimits=sortlimits;
#endif
  GPUvars.eta=eta;
#endif /*CUDASOPPORT*/
#ifdef CUDATIME
  cudaEvent_t totaltimestart,totaltimestop;
  if (timing || printtime)
    cudaTimingCreateAndStart(&totaltimestart,&totaltimestop,
                             "totaltimestart","totaltimestop");
#endif
  tmpvariables tmpvar;
  memset(&tmpvar,0,sizeof(tmpvar));
  // allocate and setup representation
  MPexp_setup(&MP,N,NE,zr,zi,mr,mi,er,ei,pr,
              panels,Npanel,&panelptr,pot,tol,Ndirect,false,cutoff,
              &GPUvars,&cparams,smooth,xopt,cont,&tmpvar,timing,printtime);

  // evaluation...
  StartTime(printtime,timing);


#if !defined(CUDASUPPORT) || !defined(MULTIPOLEEVALCUDA)
  //use reordered arrays aligned to 16 byte
  if(pr!=NULL) { //should always be true...
    tmpvar.prtmp=(double*)mxCalloc((N*multiplier+3),sizeof(double)); //for channelpot, symmetric interactions will cause writes to virtual points, therefore, allocate extra memory
    char *tmp=(char*)tmpvar.prtmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    tmpvar.prnew=(double*)tmp;
  }

  if(pi!=NULL) {
    tmpvar.pitmp=(double*)mxCalloc((N*multiplier+3),sizeof(double));
    char *tmp=(char*)tmpvar.pitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    tmpvar.pinew=(double*)tmp;
  }

  if(qr!=NULL) {
    tmpvar.qrtmp=(double*)mxCalloc((NE+3),sizeof(double));
    char *tmp=(char*)tmpvar.qrtmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    tmpvar.qrnew=(double*)tmp;
  }

  if(qi!=NULL) {
    tmpvar.qitmp=(double*)mxCalloc((NE+3),sizeof(double));
    char *tmp=(char*)tmpvar.qitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    tmpvar.qinew=(double*)tmp;
  }
  CudaStartTime(printtime,timing); //for a good comparison with gpu, only take the time over the actual evaluation
  // evaluation at potentials (zr,zi)
  if (pr != NULL) MPexp_eval2(&MP,tmpvar.prnew,tmpvar.pinew,tmpvar.zrnew,tmpvar.zinew,tmpvar.mrnew,tmpvar.minew,panels,
                             pot,smooth,xopt,cutoff,cont,&cparams);

  // evaluation at coordinates (er,ei)
  if (qr != NULL) MPexp_eval(&MP,tmpvar.qrnew,tmpvar.qinew,tmpvar.ernew,tmpvar.einew,tmpvar.zrnew,tmpvar.zinew,tmpvar.mrnew,tmpvar.minew,panels,
                             pot,smooth,xopt,cutoff,cont,&cparams);
#ifdef CHANNELPOT
  if(cparams.H>0)
    streaminteraction((void*)&MP,&cparams,tmpvar.prnew,tmpvar.pinew,tmpvar.qrnew,tmpvar.qinew,
                      tmpvar.zrnew,tmpvar.zinew,tmpvar.ernew,tmpvar.einew,tmpvar.mrnew,tmpvar.minew,N,NE,panels);
#endif
#ifndef CUDASUPPORT
 CudaStopTime(printtime,timing);
#endif
  if(pr!=NULL) {
    if(pi!=NULL)
      for(int i=0;i<N;i++) {
        pr[MP.ix[i]]=tmpvar.prnew[i];
        pi[MP.ix[i]]=tmpvar.pinew[i];
      }
    else
      for(int i=0;i<N;i++)
        pr[MP.ix[i]]=tmpvar.prnew[i];
  }
  if(qr != NULL) {
    if(qi!=NULL)
      for(int i=0;i<NE;i++) {
        qr[MP.jx[i]]=tmpvar.qrnew[i];
        qi[MP.jx[i]]=tmpvar.qinew[i];
      }
    else
      for(int i=0;i<NE;i++)
        qr[MP.jx[i]]=tmpvar.qrnew[i];
  }
  mxFree(tmpvar.zrtmp);
  mxFree(tmpvar.zitmp);
  mxFree(tmpvar.mrtmp);
  mxFree(tmpvar.mitmp);
  mxFree(tmpvar.ertmp);
  mxFree(tmpvar.eitmp);
  mxFree(tmpvar.prtmp);
  mxFree(tmpvar.pitmp);
  mxFree(tmpvar.qrtmp);
  mxFree(tmpvar.qitmp);
#else
  CudaStartTime(printtime,timing);
#endif /*!defined(CUDASUPPORT) || !defined(MULTIPOLEEVALCUDA)*/

#if defined(CUDASUPPORT) && defined(MULTIPOLEEVALCUDA)
  // perform the evaluation of multipole coefficients on the GPU
  mpexp_eval_cuda(&MP,&GPUvars,pr,N,NE,pot,timing,printtime);
#endif

  StopTime(printtime,timing);
#if !defined(SSE2DIRECT) || defined(CUDASUPPORT)
  CudaStopTime(printtime,timing);
#endif
  PrintTime("Evaluation: %f\n",timing,18,printtime);
  CudaPrintTime("Evaluation CUDA time",timing,8,printtime);

#ifdef CUDASUPPORT
  // copy all results from the GPU back to the CPU
  MPexp_eval_cuda_result(qr,qi,pr,pi,MP.ix,MP.jx,N,NE,pot,
                         &GPUvars,&MP,timing,printtime);
#endif
#ifdef CUDATIME
  TimeEventRecord(totaltimestop);
  if (timing || printtime) {
    TimeSyncPrintAndDestroy(totaltimestart,totaltimestop,9,"Full time");
  }
#endif
  // dellocate
  if (panelptr != NULL) {
    for (int i = 0; i <= MP.nlevel; i++) {
      if (panelptr[i] != NULL)
        mxFree(panelptr[i]);
    }
    mxFree(panelptr);
  }

  MPexp_free(&MP);

#ifdef CUDASUPPORT
  if (getcudaalloccount() != 0)
    mexPrintf("Warning: %d memory allocations on GPU not deallocated.\n",
              getcudaalloccount());
  CHECKCUDAERROR2
#endif
#ifdef CHANNELPOT
  if(cparams.H>0)
    cleanupchannelparam(&cparams);
#endif
}

#ifndef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
/*------------------------------------------------------------------------*/
void fmm2dGetMesh(int nlhs,mxArray *plhs[],
                  int N,
                  const double *zr,const double *zi,
                  const double *mr,const double *mi,
                  int NE,
                  const double *er,const double *ei,
                  const panel *panels,int Npanel,
                  int pot,double tol,int Ndirect,double cutoff,
                  double *timing,bool printtime,
                  float *sortlimits,double* eta,double channelheight)
/* Returns the multipole mesh in a cell vector of complex
   mxArrays. Additionally, all sparse connection matrices are
   optionally returned. For debugging and pedagogical purposes. */
{
#ifdef PANELSORT
  if (N+Npanel <= 0) return;
#else
  if (N <= 0) return;
#endif

  // setup multipole tree
  MPexp MP;
  int **panelptr = NULL;

  cudavariables GPUvars;
  channelparam cparams;
#ifdef CHANNELPOT
  cparams.H=channelheight;
#endif
#ifdef CUDASUPPORT
  memset(&GPUvars,0,sizeof(GPUvars));
  cudastart();
  checkcudaerror("Check before first cuda call");
  if (timing || printtime)
    cudaTimingCreateAndStart(&GPUvars.sortstart,&GPUvars.sortstop,
                             "sortstart","sortstop");
  checkcudaerror("Check after first cuda call");
#ifdef SORTLIMIT
  GPUvars.sortlimits=sortlimits;
#endif
  GPUvars.eta=eta;
#endif

  MPexp_setup(&MP,N,NE,zr,zi,mr,mi,er,ei,NULL,
              panels,Npanel,&panelptr,pot,tol,Ndirect,true,cutoff,
              &GPUvars,&cparams,DIRAC,0,false,NULL,timing,printtime);

  int multiplier=1;
#ifdef CHANNELPOT
  if(cparams.H>0) {
    cleanupchannelparam(&cparams);
    multiplier=3;
  }
  const int Nf = 1 << (MP.nlevel << 1);
  const int Nt = (Nf-1)/3+Nf;
#endif
  // create cell vector and copy the FMM mesh
  plhs[0] = mxCreateCellMatrix(1,MP.nlevel+1);
  for (int i = 0; i <= MP.nlevel; i++) {
    mxArray *cell = mxCreateDoubleMatrix(2,(MP.lptr[i+1]-MP.lptr[i])*multiplier,
                                         mxCOMPLEX);
    double *mshr = mxGetPr(cell),*mshi = mxGetPi(cell);
    mpexp *offset = &MP.root[MP.lptr[i]];

    for (int j = 0; j < MP.lptr[i+1]-MP.lptr[i]; j++) {
      mshr[2*j] = creal(offset[j].z0);
      mshi[2*j] = cimag(offset[j].z0);
      mshr[2*j+1] = creal(offset[j].d0);
      mshi[2*j+1] = cimag(offset[j].d0);
    }
#ifdef CHANNELPOT
    if(cparams.H>0) {
      int base=MP.lptr[i+1]-MP.lptr[i];
      for (int j = 0; j < MP.lptr[i+1]-MP.lptr[i]; j++) {
        mshr[2*(j+base)] = creal(offset[j+Nt].z0);
        mshi[2*(j+base)] = cimag(offset[j+Nt].z0);
        mshr[2*(j+base)+1] = creal(offset[j+Nt].d0);
        mshi[2*(j+base)+1] = cimag(offset[j+Nt].d0);
        mshr[2*(j+2*base)] = creal(offset[j+2*Nt].z0);
        mshi[2*(j+2*base)] = cimag(offset[j+2*Nt].z0);
        mshr[2*(j+2*base)+1] = creal(offset[j+2*Nt].d0);
        mshi[2*(j+2*base)+1] = cimag(offset[j+2*Nt].d0);
      }
    }
#endif
    mxSetCell(plhs[0], i, cell);
  }

  // additionally, the connection matrices
  if (nlhs > 1) {
    plhs[1] = mxCreateCellMatrix(1,MP.nlevel+1);
    for (int i = 0; i <= MP.nlevel; i++) {
      const mpSparse *S = &MP.connect[i];
      const int N = MP.lptr[i+1]-MP.lptr[i];
      mxArray *cell = mxCreateSparse(N*multiplier,N,S->jcptr[N],mxREAL);
      mwSize *jc = mxGetJc(cell);
      mwIndex *ir = mxGetIr(cell);
      double *pr = mxGetPr(cell);
#ifdef CHANNELPOT
      int Nfl=1 << (i << 1);
#endif

      // mwIndex not defined? use '-DmwIndex=int' in the mex-command!

      // note: interaction through parents as in the C99-code is not
      // implemented here (very little efficiency gain but more
      // complicated memory access pattern)
      jc[0] = 0;
      for (size_t j = 0,k; j < N; j++) {
        jc[j+1] = S->jcptr[j+1];

        // near connections
        if (i < MP.nlevel)
          for (k = jc[j]; k < S->kcptr[j]; k++) {
            ir[k] = S->ir[k];
            pr[k] = 1.0;
#ifdef CHANNELPOT
            if(cparams.H>0) {
              if(ir[k]>=N)
                ir[k]-=(Nt-Nfl);
              if(ir[k]>=2*N)
                ir[k]-=(Nt-Nfl);
            }
#endif
          }
        else {
          // final level contains more information:
          for (k = jc[j]; k < S->kcptr[j+N]; k++) {
            ir[k] = S->ir[k];
            pr[k] = 1.0;
#ifdef CHANNELPOT
            if(cparams.H>0) {
              if(ir[k]>=N)
                ir[k]-=(Nt-Nfl);
              if(ir[k]>=2*N)
                ir[k]-=(Nt-Nfl);
            }
#endif
          }
          for ( ; k < S->kcptr[j]; k++) {
            // two cases here; when debugging it is occasionally
            // convenient to output -2 in the first case, but this is
            // not particularly logical (1 and -1 are truly different,
            // but 2 and -2 are just two sides of the same coin...)
            if (S->ir[k] < 0) {
              ir[k] = -S->ir[k]; // must clear sign
              pr[k] = 2.0; // particle-to-cluster
            }
            else {
              ir[k] = S->ir[k];
              pr[k] = 2.0; // cluster-to-particle
            }
#ifdef CHANNELPOT
            if(cparams.H>0) {
              if(ir[k]>=N)
                ir[k]-=(Nt-Nfl);
              if(ir[k]>=2*N)
                ir[k]-=(Nt-Nfl);
            }
#endif
          }
        }

        // far connections
        for ( ; k < S->jcptr[j+1]; k++) {
          ir[k] = S->ir[k];
#ifdef CHANNELPOT
          if(cparams.H>0) {
            if(ir[k]>=N)
              ir[k]-=(Nt-Nfl);
            if(ir[k]>=2*N)
              ir[k]-=(Nt-Nfl);
          }
#endif
          pr[k] = -1.0;
        }
      }
      mxSetCell(plhs[1],i,cell);
    }
  }

  // deallocate and return
  if (panelptr != NULL) {
    for (int i = 0; i <= MP.nlevel; i++) {
      if (panelptr[i] != NULL)
        mxFree(panelptr[i]);
    }
    mxFree(panelptr);
  }
  MPexp_free(&MP);

#ifdef CUDASUPPORT
  cleanupsort(&GPUvars,&MP);
  if (getcudaalloccount() != 0)
    mexPrintf("Warning: %d memory allocations on GPU not deallocated.\n",
              getcudaalloccount());
#endif
}
#endif
/*------------------------------------------------------------------------*/
void MPexp_setup(MPexp *This,int N,int NE,
                 const double *zr,const double *zi,
                 const double *mr,const double *mi,
                 const double *er,const double *ei,const double *pr,
                 const panel *panels,int Npanel,int ***panelptrlist,
                 int pot, double tol,int Ndirect,bool meshonly,double cutoff,
                 cudavariables *GPUvars,channelparam *cparams,SMOOTHER smooth,double xopt,bool cont,
                 tmpvariables *tmpvar,double* timing,bool printtime)
/* Setups and allocates a multipole tree for N potentials at (zr,zi)
   to be evaluated at the NE points (er,ei) to tolerance tol. The
   discretization is determined so that the number of potentials in
   each box at the finest level is <= Ndirect. If the boolean variable
   meshonly is true, then only the multipole tree and the connectivity
   information is determined. */
{
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);

  // size of the domain
  int validinput=1;
#ifndef CUDASORT
  MPexp_box_(This,N,NE,zr,zi,er,ei,&validinput);
  MPexp_box_panel(This,panels,Npanel,&validinput);
#endif
#ifdef CHECKNANINPUT
  if(!validinput)
    mexErrMsgTxt("NaN detected in input vectors, aborting");
#endif

  // number of multipole terms
  This->pcoeff = (int)floor(log2(4.0*tol)/log2(theta)-1.0)+(pot != 0);
  /* the magical constants 4.0 is a (not overly optimistic) estimate
     obtained from numerous experiments with the model error ~
     constant * theta^(p+1) and theta = 0.5 */
#ifdef SSE2M2PSSCALING
  if (This->pcoeff < 2) This->pcoeff = 2;//due to other ordering of elements, SSE2M2PSSCALING does not work for pcoeff=1
#else
  if (This->pcoeff < 1) This->pcoeff = 1;
#endif

  // number of levels
#ifdef PANELSORT // additional levels may be necessary for many panels
  This->nlevel = (int)ceil(0.5*log2(5.0/8.0*
                               (double)(N+PANELDUMMYFACTOR*Npanel)/Ndirect));
#else
  This->nlevel = (int)ceil(0.5*log2(5.0/8.0*(double)N/Ndirect));
#endif
  if (This->nlevel < 0) This->nlevel = 0;
  /* the magical constant 5/8 comes from solving for the number of
     levels the equation E(N/Nf) = Ndirect, where E is expectation,
     N is considered a uniform random number on [1,M] (with M
     tending to infinity) and Nf is defined below */

  /* used frequently below: total number of coefficients, number of
     boxes at the finest level and total number of boxes */

  int multiplier=1;
#ifdef CHANNELPOT
  if(cparams->H>0) {
    setchannelparams(cparams,This->xmax-This->xmin,tol);
    if(cparams->Nlev/2>This->nlevel) //Number of levels required by stream expansions, note: only one split/level for stream expansions, therefore division by 2
      This->nlevel=cparams->Nlev/2;
    multiplier=3;
    if(This->ymax>cparams->H||This->ymin<0)
      mexPrintf("All points are not within the channel");
  }
#endif
  const int P = This->pcoeff+1;
  const int Nf = 1 << (This->nlevel << 1);
  int Nt = ((Nf-1)/3+Nf)*multiplier;//two dummy levels as well, put after the rest of the elements

  // pointer to levels
  This->lptr = (int *)mxMalloc((This->nlevel+2)*sizeof(int));
  This->lptr[0] = 0;
  for (int level = 0, nb = 1; level <= This->nlevel; level++, nb <<= 2)
    This->lptr[level+1] = This->lptr[level]+nb;

  if(!fullcudafmm||!cudasort||meshonly) { //If GPU sorts and evaluates, no need to allocate on CPU
    // multipole tree
    This->root = (mpexp *)mxMalloc(Nt*sizeof(mpexp));

    // multipole coefficients
    dcmplx *coeffs1 = (dcmplx *)mxCalloc(P*Nt, sizeof(dcmplx));
    dcmplx *coeffs2 = (dcmplx *)mxCalloc(P*Nt, sizeof(dcmplx));
    for (int j = 0; j < Nt; j++) {
      This->root[j].coeff1 = coeffs1; coeffs1 += P;
      This->root[j].coeff2 = coeffs2; coeffs2 += P;
    }

    // permutation of input
    This->ixptr = (int *)mxMalloc((Nf*multiplier+1)*sizeof(int));
    This->ix = (int *)mxMalloc(N*sizeof(int));
    if (NE) {
      This->jxptr = (int *)mxMalloc((Nf+1)*sizeof(int));
      This->jx = (int *)mxMalloc(NE*sizeof(int));
    }
    else
      This->jxptr = This->jx = NULL;
  }
  else {
    This->root=NULL;
    This->ixptr=NULL;
    This->ix=NULL;
    This->jxptr=NULL;
    This->jx=NULL;
  }
  // connections between boxes at each level
  This->connect = (mpSparse *)mxMalloc((This->nlevel+1)*sizeof(mpSparse));
  *panelptrlist = (int **)mxCalloc((This->nlevel+1),sizeof(int *));

  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Alloc: %f\n",timing,10,printtime);
  CudaPrintTime("Alloc CUDA time",timing,0,printtime);

  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);
  /* setup the multipole mesh, sort all potentials and evaluation
     points and determine all connections */
  int maxm2p,levm2p;
#ifdef CUDASUPPORT
  checkcudaerror("Check before MPexp_sort_.");
#endif
#if defined(CUDASUPPORT)/* && defined(CUDASORT)*/
  GPUvars->connect =
    (mpSparse *)mxMalloc((This->nlevel+1)*sizeof(mpSparse));
#endif
  MPexp_sort_(This,N,zr,zi,NE,er,ei,panels,Npanel,*panelptrlist,
              &maxm2p,&levm2p,Ndirect,cutoff,GPUvars,cparams,timing,printtime,meshonly);
  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Sort: %f\n",timing,13,printtime);
  CudaPrintTime("Sort CUDA time",timing,3,printtime);

  if (meshonly) return;
#if !defined(CUDASUPPORT) || !defined(MULTIPOLEEVALCUDA) //use reordered arrays aligned to 16 byte
    if(zr!=NULL) { //should always be true...
      tmpvar->zrtmp=(double*)mxMalloc((N*multiplier+3)*sizeof(double));
      char *tmp=(char*)tmpvar->zrtmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->zrnew=(double*)tmp;
    }

    if(zi!=NULL) { //should always be true...
      tmpvar->zitmp=(double*)mxMalloc((N*multiplier+3)*sizeof(double));
      char *tmp=(char*)tmpvar->zitmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->zinew=(double*)tmp;
    }

    if(mr!=NULL) { //should always be true...
      tmpvar->mrtmp=(double*)mxMalloc((N*multiplier+3)*sizeof(double));
      char *tmp=(char*)tmpvar->mrtmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->mrnew=(double*)tmp;
    }

    if(mi!=NULL) {
      tmpvar->mitmp=(double*)mxMalloc((N*multiplier+3)*sizeof(double));
      char *tmp=(char*)tmpvar->mitmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->minew=(double*)tmp;
    }

    if(er!=NULL) {
      tmpvar->ertmp=(double*)mxMalloc((NE+3)*sizeof(double));
      char *tmp=(char*)tmpvar->ertmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->ernew=(double*)tmp;
    }

    if(ei!=NULL) {
      tmpvar->eitmp=(double*)mxMalloc((NE+3)*sizeof(double));
      char *tmp=(char*)tmpvar->eitmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->einew=(double*)tmp;
    }
    if(mi==NULL) {
      for(int i=0;i<N;i++) {
        tmpvar->zrnew[i]=zr[This->ix[i]];
        tmpvar->zinew[i]=zi[This->ix[i]];
        tmpvar->mrnew[i]=mr[This->ix[i]];
      }
    }
    else {
      tmpvar->mitmp=(double*)mxMalloc((N*multiplier+3)*sizeof(double));
      char *tmp=(char*)tmpvar->mitmp;
      if(((size_t)tmp&15)!=0) {
        tmp+=16-((size_t)tmp&15);
      }
      tmpvar->mitmp=(double*)tmp;
      for(int i=0;i<N;i++) {
        tmpvar->zrnew[i]=zr[This->ix[i]];
        tmpvar->zinew[i]=zi[This->ix[i]];
        tmpvar->mrnew[i]=mr[This->ix[i]];
        tmpvar->minew[i]=mi[This->ix[i]];
      }
    }
    for(int i=0;i<NE;i++) {
      tmpvar->ernew[i]=er[This->jx[i]];
      tmpvar->einew[i]=ei[This->jx[i]];
    }
#endif /*defined(CUDASUPPORT) || !defined(MULTIPOLEEVALCUDA)*/
#ifdef CUDASUPPORT
  checkcudaerror("Check before MPexp_eval_cuda");
  // as soon as the sorting is completed, the GPU direct summation can begin
  MPexp_eval_cuda(This,er,ei,zr,zi,mr,mi,pr,GPUvars,N,NE,pot,
                  smooth,xopt,cutoff,cont,timing,printtime);
  // (this function also calls cuda_init and cuda_m2m)
#endif

#if !defined(CUDASUPPORT) || !defined(MULTIPOLEINITCUDA) ||  \
  !defined (MULTIPOLESHIFTCUDA) || defined(CHECKINITCUDA) || \
  defined(CHECKM2PS) || defined(CHECKP2P)
  // if CPU-shift will be used, prepare for this
  shift_alloc(This->pcoeff,pot,maxm2p);
#endif

#ifdef CHANNELPOT
  if(cparams->H>0) {
    replicatepositions((void *)This,cparams,(void*)tmpvar,N);
    replicatepanelindices(This,cparams,Npanel);
  }
#endif
  /*******************************************/
  /* Start of the actual multipole algorithm */
  /*******************************************/

#if !defined(CUDASUPPORT) || !defined(MULTIPOLEINITCUDA) || \
  defined(CHECKINITCUDA)

  /* Step #1: Initialization at the finest level. */
    MPexp_init_CPU_(This,Nf,Nt,tmpvar->zrnew,tmpvar->zinew,tmpvar->mrnew,tmpvar->minew,panels,cparams,pot,timing,printtime);

  /* Step #2: Upward pass of the fast multipole algorithm. After the
     finest level has been initialized, This loop performs the
     'upward' pass where all multipole contributions are shifted and
     added onto all coarser levels. */
  MPexp_m2m_CPU_(This,Nf,pot,timing,printtime);

//   checkboxes(This,zr,zi,er,ei,Nf,Nt,NE);
#ifdef CHANNELPOT
  if(cparams->H>0) {
    replicatecoeffs((void *)This,cparams);
    addmirrorcorrection((void *)This,cparams);
  }
#endif
#ifdef CHECKINITCUDA
  cuda_check_init(This,Nt,GPUvars);
#endif
#endif

  /* Step #3: Downward pass of the fast multipole method. This
     shifts expansions from the interaction list of each box to the
     center of the box and adds up the result. */
#if defined(CHECKM2PS) || defined(CHECKP2P) || \
  !defined(MULTIPOLESHIFTCUDA)
  // CPU code if not MULTIPOLESHIFTCUDA, or if shift validation should
  // be performed
  MPexp_m2ps_CPU_(This,levm2p,timing,printtime);
#endif
#ifdef MULTIPOLESHIFTCUDA
  MPexp_shiftm2ps_cuda(This,GPUvars,levm2p,pot,timing,printtime);
#endif
#ifdef CHANNELPOT
  if(cparams->H>0) {
    levm2p=cparams->Nlev/2; //for p2p, there will be coeff2 values at this level due to mirror correction
  }
#endif
  /* Before continuing to finer levels, each parent shifts and adds
     its expansion to each of its children. */
#if defined(CHECKP2P) || !defined(MULTIPOLESHIFTCUDA)
  MPexp_p2p_CPU_(This,levm2p,timing,printtime);
#endif
#ifdef MULTIPOLESHIFTCUDA
  MPexp_shiftp2p_cuda_scaled(This,GPUvars,levm2p,timing,printtime);
#endif
}
/*------------------------------------------------------------------------*/
void MPexp_init_CPU_(MPexp *restrict This,int Nf,int Nt,
                     const double *restrict zr,const double *restrict zi,
                     const double *restrict mr,const double *restrict mi,
                     const panel *restrict panels,channelparam* cparams,int pot,
                     double *restrict timing,bool printtime)
/* Step #1: CPU-initialization at the finest level of the tree. */
{
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);

  mpexp *offset = &This->root[This->lptr[This->nlevel]];
  mpSparse *C = &This->connect[This->nlevel];
#if defined(CHECKINITCUDA) && defined(MULTIPOLEINITCUDA)
  //cuda debug purpose
  memset(This->root->coeff2,0,Nt*(This->pcoeff+1)*sizeof(dcmplx));
  memset(This->root->coeff1,0,Nt*(This->pcoeff+1)*sizeof(dcmplx));
#endif
#ifdef CHANNELPOT
  const int Ntl = (Nf-1)/3+Nf; //original Nt multiplied by 3
#endif
  if (mi == NULL) // real mass
    for (int j = 0; j < Nf; j++) {
      // initialize multipole expansion in box j
      int i = This->ixptr[j];
#ifdef SSE2INIT
        for (; i < This->ixptr[j+1]-1; i+=2)
          mpexp_init_sse(&offset[j],This->pcoeff,
                        zr[i],zi[i],zr[i+1],zi[i+1],
                        mr[i],0.0,mr[i+1],0.0,pot);
#endif
        for (; i < This->ixptr[j+1]; i++)
          mpexp_init(&offset[j],This->pcoeff,
                     zr[i],zi[i],
                     mr[i],0.0,pot);

      // panels in box j
      panel smallpanel;
      for (int i = 0; i < offset[j].npanel; i++)
        if (panelinbox(&smallpanel,panels+offset[j].panelptr[i],
                       offset[j].z0,offset[j].d0))
          expandpanel(&smallpanel,offset[j].z0,offset[j].coeff1,This->pcoeff);

      // 'less near' connections: shift into polynomial directly
      for (int jj = C->kcptr[j+Nf]; jj < C->kcptr[j]; jj++) {
        int k = C->ir[jj];
        if (k < 0) {
#ifdef CHANNELPOT
          if(cparams->H>0) {
            if(-k>=Nf)
              continue; //no use initializing at a point where no evaluation will be performed
          }
#endif
          int i = This->ixptr[j];
#ifdef SSE2INIT
            for (; i < This->ixptr[j+1]-1; i+=2)
              mpexp_initp_sse(&offset[-k],This->pcoeff,
                              zr[i],zi[i],zr[i+1],zi[i+1],
                              mr[i],0.0,mr[i+1],0.0,pot);
#endif
            for (; i < This->ixptr[j+1]; i++)
              mpexp_initp(&offset[-k],This->pcoeff,
                          zr[i],zi[i],
                          mr[i],0.0,pot);
          // ditto panels
          panel smallpanel;
          for (int i = 0; i < offset[j].npanel; i++)
            if (panelinbox(&smallpanel,panels+offset[j].panelptr[i],
                           offset[j].z0,offset[j].d0))
              farexpandpanel(&smallpanel,offset[-k].z0,offset[-k].coeff2,
                             This->pcoeff);
        }
        else {
#ifdef CHANNELPOT
          if(cparams->H>0) {
            if(k>=Nf)
              k-=(Ntl-Nf);
            if(k>=2*Nf)
              k-=(Ntl-Nf);
          }
#endif
          int i = This->ixptr[k];
#ifdef SSE2INIT
            for (; i < This->ixptr[k+1]-1; i+=2)
              mpexp_initp_sse(&offset[j],This->pcoeff,
                          zr[i],zi[i],zr[i+1],zi[i+1],
                          mr[i],0.0,mr[i+1],0.0,pot);
#endif
            for (; i < This->ixptr[k+1]; i++)
              mpexp_initp(&offset[j],This->pcoeff,
                          zr[i],zi[i],
                          mr[i],0.0,pot);
          panel smallpanel;
          for (int i = 0; i < offset[k].npanel; i++)
            if (panelinbox(&smallpanel,panels+offset[k].panelptr[i],
                           offset[k].z0,offset[k].d0))
              farexpandpanel(&smallpanel,offset[j].z0,offset[j].coeff2,
                             This->pcoeff);
        }
      }
    }
  else // complex mass, mi != NULL
    for (int j = 0; j < Nf; j++) {
      int i = This->ixptr[j];
#ifdef SSE2INIT
        for (; i < This->ixptr[j+1]-1; i+=2)
          mpexp_init_sse(&offset[j], This->pcoeff,
                         zr[i], zi[i],zr[i+1], zi[i+1],
                         mr[i], mi[i],mr[i+1], mi[i+1], pot);
#endif
        for (; i < This->ixptr[j+1]; i++)
          mpexp_init(&offset[j], This->pcoeff,
                     zr[i], zi[i],
                     mr[i], mi[i], pot);
      panel smallpanel;
      for (int i = 0; i < offset[j].npanel; i++)
        if (panelinbox(&smallpanel,panels+offset[j].panelptr[i],
                       offset[j].z0,offset[j].d0))
          expandpanel(&smallpanel,offset[j].z0,offset[j].coeff1,This->pcoeff);

      for (int jj = C->kcptr[j+Nf]; jj < C->kcptr[j]; jj++) {
        int k = C->ir[jj];
        if (k < 0) {
#ifdef CHANNELPOT
          if(cparams->H>0) {
            if(-k>=Nf)
              continue;
          }
#endif
          int i = This->ixptr[j];
#ifdef SSE2INIT
            for (; i < This->ixptr[j+1]-1; i+=2)
              mpexp_initp_sse(&offset[-k],This->pcoeff,
                              zr[i],zi[i],zr[i+1],zi[i+1],
                              mr[i],mi[i],mr[i+1],mi[i+1],pot);
#endif
            for (; i < This->ixptr[j+1]; i++)
              mpexp_initp(&offset[-k],This->pcoeff,
                          zr[i],zi[i],
                          mr[i],mi[i],pot);
          panel smallpanel;
          for (int i = 0; i < offset[j].npanel; i++)
            if (panelinbox(&smallpanel,panels+offset[j].panelptr[i],
                           offset[j].z0,offset[j].d0))
              farexpandpanel(&smallpanel,offset[-k].z0,offset[-k].coeff2,
                             This->pcoeff);
        }
        else {
#ifdef CHANNELPOT
          if(cparams->H>0) {
            if(k>=Nf)
              k-=(Ntl-Nf);
            if(k>=2*Nf)
              k-=(Ntl-Nf);
          }
#endif
          int i = This->ixptr[k];
#ifdef SSE2INIT
            for (; i < This->ixptr[k+1]-1; i+=2)
              mpexp_initp_sse(&offset[j],This->pcoeff,
                          zr[i],zi[i],zr[i+1],zi[i+1],
                          mr[i],mi[i],mr[i+1],mi[i+1],pot);
#endif
            for (; i < This->ixptr[k+1]; i++)
              mpexp_initp(&offset[j],This->pcoeff,
                          zr[i],zi[i],
                          mr[i],mi[i],pot);
          panel smallpanel;
          for (int i = 0; i < offset[k].npanel; i++)
            if (panelinbox(&smallpanel,panels+offset[k].panelptr[i],
                           offset[k].z0,offset[k].d0))
              farexpandpanel(&smallpanel,offset[j].z0,offset[j].coeff2,
                             This->pcoeff);
        }
      }
    }
  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Step #1: %f\n",timing,14,printtime);
  CudaPrintTime("Step #1 CUDA time",timing,4,printtime);
}
/*------------------------------------------------------------------------*/
void MPexp_m2m_CPU_(MPexp *restrict This,int Nf,int pot,
                    double *restrict timing,bool printtime)
/* Step #2: CPU-upward pass of the fast multipole algorithm. After the
   finest level has been initialized, this loop performs the 'upward'
   pass where all multipole contributions are shifted and added onto
   all coarser levels. */
{
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);

  for (int l = This->nlevel, nb = Nf; l > 0; l--, nb >>= 2) {
    mpexp *parent = &This->root[This->lptr[l-1]];
    mpexp *child = &This->root[This->lptr[l]];

    // loop over all children at This level
    for (int j = 0; j < nb; j += 4) {
      // shift and add result from 4 children at once
      shift_m2m(child[j].z0-parent->z0,child[j+1].z0-parent->z0,
                child[j+2].z0-parent->z0,child[j+3].z0-parent->z0,
                parent->coeff1,child[j].coeff1);
      parent++;
    }
  }

  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Step #2: %f\n",timing,15,printtime);
  CudaPrintTime("Step #2 CUDA time",timing,5,printtime);
}
/*------------------------------------------------------------------------*/
void MPexp_m2ps_CPU_(MPexp *restrict This,int levm2p,
                     double *restrict timing,bool printtime)
/* Step #3(1/2): CPU-downward pass of the fast multipole method. This
   shifts expansions from the interaction list of each box to the
   center of the box and adds up the result. */
{
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);

  for (int l = levm2p-1, nb = 1 << 2*levm2p; l < This->nlevel; l++, nb <<= 2) {
    mpexp *child = &This->root[This->lptr[l+1]];
    mpSparse *C = &This->connect[l+1];

    /* tricky part: loop over all boxes (parents) at this level and
       symmetrically shift results to and from the interaction list of
       each child */
    for (int j = 0; j < nb; )
      // loop through the interaction list of each child
      for (int c = 0; c < 4; c++, j++)
        shift_m2ps(child[j].z0,(void *)&child[0].z0,sizeof(mpexp),
                   C->ir,C->kcptr[j],C->jcptr[j+1],
                   child[j].coeff1,child[j].coeff2,
                   child[0].coeff1,child[0].coeff2);
  }

  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Step #3(1/2): %f\n",timing,16,printtime);
  CudaPrintTime("Step #3(1/2) CUDA time",timing,6,printtime);
}
/*------------------------------------------------------------------------*/
void MPexp_p2p_CPU_(MPexp *restrict This,int levm2p,
                    double *restrict timing,bool printtime)
/* Step #3(2/2): Before continuing to finer levels, each parent shifts
   and adds its expansion to each of its children. */
{
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);

  for (int l = levm2p-1, nb = 1 << 2*levm2p; l < This->nlevel; l++, nb <<= 2) {
    mpexp *parent = &This->root[This->lptr[l]];
    mpexp *child = &This->root[This->lptr[l+1]];

    for (int j = 0; j < nb; j += 4) {
      // contribution from parent to 4 children at once
      if (l >= levm2p)
        // need not perform p2p before there has been a m2ps-call
        // (note: since the C99-version allows for interactions
        // through parents such a shift could have occurred at the
        // previous level)
        shift_p2p(parent->z0-child[j].z0,parent->z0-child[j+1].z0,
                  parent->z0-child[j+2].z0,parent->z0-child[j+3].z0,
                  child[j].coeff2,parent->coeff2);
      parent++;
    }
  }

  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Step #3(2/2): %f\n",timing,17,printtime);
  CudaPrintTime("Step #3(2/2) CUDA time",timing,7,printtime);
}
/*------------------------------------------------------------------------*/
void MPexp_free(MPexp *This)
/* Deallocates the whole tree. */
{
  for (int l = This->nlevel; l >= 0; l--)
    mxFree(This->connect[l].jcptr);
  mxFree(This->connect);
  mxFree(This->jx);
  mxFree(This->jxptr);
  mxFree(This->ix);
  mxFree(This->ixptr);
  if(This->root!=NULL) {
    mxFree(This->root->coeff2);
    mxFree(This->root->coeff1);
    mxFree(This->root);
  }
  mxFree(This->lptr);
}
/*------------------------------------------------------------------------*/
void MPexp_smooth_(int pot,SMOOTHER smooth,double xopt,bool cont,
                   double *cutoff,double *shape,double *scale)
/* Precomputes the three parameters (cutoff,shape,scale) given
   (pot,smooth,xopt,cutoff,cont). This yields a slicker code for
   mpexp_directInteract(). */
{
  *scale = 1.0;

  if (pot == 0)
    switch (smooth) {
    case DIRAC:
      break;

    case RANKINE:
#ifdef INLINECOMPLEX_DIRECT
      *cutoff = xopt*xopt;
      *shape = 0.5*log(*cutoff)-0.5;
#else
      *cutoff = xopt;
      *shape = log(*cutoff)-0.5;
#endif
      break;

    case SCULLY:
#ifdef INLINECOMPLEX_DIRECT
      *cutoff = *cutoff**cutoff;
      *shape = xopt*xopt;
      if (cont)
        *scale = log(*cutoff)/log(*shape+*cutoff);
#else
      *shape = xopt*xopt;
      if (cont)
        *scale = 2.0*log(*cutoff)/log(*shape+*cutoff**cutoff);
#endif
      break;

    case OSEEN:
#ifdef INLINECOMPLEX_DIRECT
    *cutoff = *cutoff**cutoff;
    *shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      *scale = log(*cutoff)/(log(*cutoff)+expint(*shape**cutoff));
#else
    *shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      *scale = log(*cutoff)/(log(*cutoff)+0.5*expint(*shape**cutoff**cutoff));
#endif
      break;
    }
  else
    switch (smooth) {
    case DIRAC:
      break;

    case RANKINE:
#ifdef INLINECOMPLEX_DIRECT
      *cutoff = xopt*xopt;
      *shape = 1.0/(*cutoff);
#else
      *cutoff = xopt;
      *shape = 1.0/(*cutoff**cutoff);
#endif
      break;

    case SCULLY:
#ifdef INLINECOMPLEX_DIRECT
      *cutoff = *cutoff**cutoff;
      *shape = xopt*xopt;
      if (cont)
        *scale = 1.0+*shape/(*cutoff);
#else
      *shape = xopt*xopt;
      if (cont)
        *scale = 1.0+*shape/(*cutoff**cutoff);
#endif
      break;

    case OSEEN:
#ifdef INLINECOMPLEX_DIRECT
      *cutoff = *cutoff**cutoff;
      *shape = 1.2564312086261696770/(xopt*xopt);
      if (cont)
        *scale = -1.0/expm1(-*shape**cutoff);
#else
      *shape = 1.2564312086261696770/(xopt*xopt);
      if (cont)
        *scale = -1.0/expm1(-*shape**cutoff**cutoff);
#endif
      break;
    }
}
/*------------------------------------------------------------------------*/
void MPexp_eval2(const MPexp *This,
                 double *restrict pr,double *restrict pi,
                 const double *restrict zr,const double *restrict zi,
                 const double *restrict mr,const double *restrict mi,const panel *restrict panels,
                 int pot,SMOOTHER smooth,double xopt,double cutoff,bool cont,channelparam* cparams)
/* Evaluation of the potential at the potentials themselves. Function
   is critical to performance. */
{
  // hopefully, some optimizing compiler can use this:
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  const mpexp *restrict offset = &This->root[This->lptr[This->nlevel]];
  const int *restrict ix = This->ix,*restrict ixptr = This->ixptr;
  const mpSparse *C = &This->connect[This->nlevel];
  const int *restrict jcptr = C->jcptr,*restrict j2cptr = &C->kcptr[Nf],
    *restrict kcptr = C->kcptr,*restrict ir = C->ir;
  double shape,scale;
#ifdef CHANNELPOT
  const int Nt = (Nf-1)/3+Nf;
#endif
  MPexp_smooth_(pot,smooth,xopt,cont,&cutoff,&shape,&scale);
  for (int j = 0, jj; j < Nf; j++) {
    // far-field: contribution from distant potentials
#ifndef MULTIPOLEEVALCUDA
#ifndef NOMULTIPOLEEVAL //Only for timings, will give wrong results
#ifdef SSE2DIRECT
    mpexp_eval_sse(&offset[j],This->pcoeff,pr,pi,zr,zi,ix,ixptr[j],ixptr[j+1]);
#else
    mpexp_eval(&offset[j],This->pcoeff,pr,pi,zr,zi,ix,ixptr[j],ixptr[j+1]);
#endif
#endif
#endif

    // near-field: symmetric interaction within list of near connections
#ifndef NODIRECTEVAL //Only for timings, will give wrong results
    for (jj = jcptr[j]; jj < j2cptr[j]; jj++) {
      int k = ir[jj];
      int kbase=k;
#ifdef CHANNELPOT
      if(cparams->H>0) {
        if(k>=Nf)
          k-=(Nt-Nf);
        if(k>=2*Nf)
          k-=(Nt-Nf);
      }
#endif
#ifdef ASYMMETRIC_CPU
      if (k != j) { // cannot interact with self here (division by zero)
#endif
#ifndef CUDASUPPORT
      mpexp_directInteract2(pr,pi,zr,zi,mr,mi,ix,
                            ixptr[j],ixptr[j+1],
                            ixptr[k],ixptr[k+1],
                            pot,smooth,cutoff,shape,scale);
#endif
      mpexp_directInteractPanel(panels,&offset[kbase],pr,pi,zr,zi,ix,
                                ixptr[j],ixptr[j+1]);
#ifdef CHANNELPOT
      if(k<Nf) //k should only be larger than Nf if channelpot is active
#endif
      mpexp_directInteractPanel(panels,&offset[j],pr,pi,zr,zi,ix,
                                ixptr[k],ixptr[k+1]);
#ifdef ASYMMETRIC_CPU
      }
#endif
    }
#else
    jj=j2cptr[j];
#endif /*NODIRECTEVAL*/

#ifndef MULTIPOLEEVALCUDA
#ifndef NOMULTIPOLEEVAL //Only for timings, will give wrong results
    // cluster to particle interactions
    for ( ; jj < kcptr[j]; jj++) {
      int k = ir[jj];
      if (k > 0) {
#ifdef CHANNELPOT
        if(cparams->H>0)
          if(k>=Nf)
            continue;
#endif
        // multipole expansion in box j evaluated at all points in box k
        #ifdef SSE2DIRECT
        mpexp_evalmp_sse(&offset[j],This->pcoeff,pot,pr,pi,zr,zi,
                     ix,ixptr[k],ixptr[k+1]);
        #else
        mpexp_evalmp(&offset[j],This->pcoeff,pot,pr,pi,zr,zi,
                     ix,ixptr[k],ixptr[k+1]);
        #endif
      }
      else {
        #ifdef SSE2DIRECT
        mpexp_evalmp_sse(&offset[-k],This->pcoeff,pot,pr,pi,zr,zi,
                     ix,ixptr[j],ixptr[j+1]);
        #else
        mpexp_evalmp(&offset[-k],This->pcoeff,pot,pr,pi,zr,zi,
                     ix,ixptr[j],ixptr[j+1]);
        #endif
      }
    }
#endif
#endif

    // near-field: direct interaction within box
#ifndef NODIRECTEVAL //Only for timings, will give wrong results
#ifndef CUDASUPPORT
    mpexp_directInteract2(pr,pi,zr,zi,mr,mi,ix,
                          ixptr[j],ixptr[j+1],-1,ixptr[j+1],
                          pot,smooth,cutoff,shape,scale);
#endif
    mpexp_directInteractPanel(panels,&offset[j],pr,pi,zr,zi,ix,
                              ixptr[j],ixptr[j+1]);
#ifdef ASYMMETRIC_CPU
    // repeat interaction with self twice (see below)
#ifndef CUDASUPPORT
    mpexp_directInteract2(pr,pi,zr,zi,mr,mi,ix,
                          ixptr[j],ixptr[j+1],-1,ixptr[j+1],
                          pot,smooth,cutoff,shape,scale);
#endif
    mpexp_directInteractPanel(panels,&offset[j],pr,pi,zr,zi,ix,
                              ixptr[j],ixptr[j+1]);
#endif
#endif /*NODIRECTEVAL*/
  }

#ifdef ASYMMETRIC_CPU
  // since asymmetric connections means everything has been done twice
  // (including the interaction with self)
  for (int i = ixptr[0]; i < ixptr[Nf]; i++)
    pr[i] *= 0.5;

  if (pi != NULL)
    for (int i = ixptr[0]; i < ixptr[Nf]; i++)
      pi[i] *= 0.5;
#endif /* ASYMMETRIC_CPU */
}
/*------------------------------------------------------------------------*/
void MPexp_eval(const MPexp *restrict This,
                double *restrict qr,double *restrict qi,
                const double *restrict er,const double *restrict ei,
                const double *restrict zr,const double *restrict zi,
                const double *restrict mr,const double *restrict mi,const panel *restrict panels,
                int pot,SMOOTHER smooth,double xopt,double cutoff,bool cont,channelparam* cparams)
/* Evaluation of the potential at the points (er,ei). Function is
   critical to performance. */
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  const mpexp *restrict offset = &This->root[This->lptr[This->nlevel]];
#ifndef CUDASUPPORT
  const int *restrict ix = This->ix,*restrict ixptr = This->ixptr;
#endif
  const int *restrict jx = This->jx,*restrict jxptr = This->jxptr;
  const mpSparse *C = &This->connect[This->nlevel];
  const int *restrict jcptr = C->jcptr,*restrict j2cptr = &C->kcptr[Nf],
    *restrict kcptr = C->kcptr,*restrict ir = C->ir;
  double shape,scale;
#ifdef CHANNELPOT
  const int Nt = (Nf-1)/3+Nf;
#endif
  MPexp_smooth_(pot,smooth,xopt,cont,&cutoff,&shape,&scale);
  for (int j = 0, jj; j < Nf; j++) {
    // far-field: contribution from distant potentials
#ifndef MULTIPOLEEVALCUDA
#ifndef NOMULTIPOLEEVAL //Only for timings, will give wrong results
#ifdef SSE2DIRECT
    mpexp_eval_sse(&offset[j],This->pcoeff,qr,qi,er,ei,jx,jxptr[j],jxptr[j+1]);
#else
    mpexp_eval(&offset[j],This->pcoeff,qr,qi,er,ei,jx,jxptr[j],jxptr[j+1]);
#endif
#endif
#endif

#ifndef NODIRECTEVAL //Only for timings, will give wrong results
    // near-field: symmetric interaction within list of near connections
    for (jj = jcptr[j]; jj < j2cptr[j]; jj++) {
      int k = ir[jj];
      int kbase=k;
#ifdef CHANNELPOT
      if(cparams->H>0) {
        if(k>=Nf)
          k-=(Nt-Nf);
        if(k>=2*Nf)
          k-=(Nt-Nf);
      }
#endif
#ifndef CUDASUPPORT
      mpexp_directInteract(qr,qi,er,ei,
                           jx,jxptr[j],jxptr[j+1],
                           zr,zi,mr,mi,
                           ix,ixptr[k],ixptr[k+1],
                           pot,smooth,cutoff,shape,scale);
#ifdef CHANNELPOT
      if(k<Nf) //k should only be larger than Nf if channelpot is active
#endif
      mpexp_directInteract(qr,qi,er,ei,
                           jx,jxptr[k],jxptr[k+1],
                           zr,zi,mr,mi,
                           ix,ixptr[j],ixptr[j+1],
                           pot,smooth,cutoff,shape,scale);
#endif
      mpexp_directInteractPanel(panels,&offset[kbase],qr,qi,er,ei,
                                jx,jxptr[j],jxptr[j+1]);
#ifdef CHANNELPOT
      if(k<Nf) //k should only be larger than Nf if channelpot is active
#endif
      mpexp_directInteractPanel(panels,&offset[j],qr,qi,er,ei,
                                jx,jxptr[k],jxptr[k+1]);
    }
#endif /*!NODIRECTEVAL*/

#ifndef MULTIPOLEEVALCUDA
#ifndef NOMULTIPOLEEVAL //Only for timings, will give wrong results
    // cluster to particle interactions
    for ( ; jj < kcptr[j]; jj++) {
      int k = ir[jj];
      if (k > 0) {
#ifdef CHANNELPOT
        if(cparams->H>0)
          if(k>=Nf)
            continue;
#endif
        // multipole expansion in box j evaluated at all points in box k
        #ifdef SSE2DIRECT
        mpexp_evalmp_sse(&offset[j],This->pcoeff,pot,qr,qi,er,ei,
                     jx,jxptr[k],jxptr[k+1]);
        #else
        mpexp_evalmp(&offset[j],This->pcoeff,pot,qr,qi,er,ei,
                     jx,jxptr[k],jxptr[k+1]);
        #endif
      }
      else {
        #ifdef SSE2DIRECT
        mpexp_evalmp_sse(&offset[-k],This->pcoeff,pot,qr,qi,er,ei,
                     jx,jxptr[j],jxptr[j+1]);
        #else
        mpexp_evalmp(&offset[-k],This->pcoeff,pot,qr,qi,er,ei,
                     jx,jxptr[j],jxptr[j+1]);
        #endif
      }
    }
#endif
#endif

#ifndef NODIRECTEVAL //Only for timings, will give wrong results
#ifndef ASYMMETRIC_CPU // has been handled above
#ifndef CUDASUPPORT
    // near-field: direct interaction within box
    mpexp_directInteract(qr,qi,er,ei,
                         jx,jxptr[j],jxptr[j+1],
                         zr,zi,mr,mi,
                         ix,ixptr[j],ixptr[j+1],
                         pot,smooth,cutoff,shape,scale);
#endif
    mpexp_directInteractPanel(panels,&offset[j],qr,qi,er,ei,jx,
                              jxptr[j],jxptr[j+1]);
#endif /* !ASYMMETRIC_CPU */
#endif /*NODIRECTEVAL*/
  }


#ifdef ASYMMETRIC_CPU // account for twice the work
  for (int j = jxptr[0]; j < jxptr[Nf]; j++)
    qr[j] *= 0.5;

  if (qi != NULL)
    for (int j = jxptr[0]; j < jxptr[Nf]; j++)
      qi[j] *= 0.5;
#endif /* ASYMMETRIC_CPU */
}
/*------------------------------------------------------------------------*/
inline void mpexp_init(mpexp *This,int p,
                const double zr,const double zi,
                const double mr,const double mi,int pot)
/* Initializes a multipole expansion at the lowest level with a single
   pointmass (mr,mi) at (zr,zi). */
{
  /* We have that the potential is given by (pot = 0)
     mp(z) = -m*log(z-zi)
           = -m*log(z-z0-(zi-z0))
           = -m*log((z-z0)*(1-(zi-z0)/(z-z0)))
           = -m*[log(z-z0)+log(1-(zi-z0)/(z-z0))]
           = -m*[log(z-z0)-sum_{j >= 1} (zi-z0)^j/j*(z-z0)^{-j},
     or othwerwise by (pot = 1)
     mp(z) = -m*1/(z-zi)
           = -m*1/(z-z0-(zi-z0))
           = -m/(z-z0)*1/(1-(zi-z0)/(z-z0))
           = -m*sum_{j >= 1} (zi-z0)^{j-1}*(z-z0)^{-j}.
  */
#ifndef INLINECOMPLEX
  const dcmplx r = (zr+I*zi)-This->z0;
  const dcmplx m = mr+I*mi;
  dcmplx *coeff = This->coeff1;

  if (pot == 0) {
    dcmplx z = r*m;
    coeff[0] -= m;
    for (int j = 1; j <= p; j++,z *= r) coeff[j] += z/j;
  }
  else {
    dcmplx z = m;
    for (int j = 1; j <= p; j++,z *= r) coeff[j] -= z;
  }
#else
  const double re = zr-creal(This->z0),im = zi-cimag(This->z0);
  dcmplx *coeff = This->coeff1;

  if (pot == 0) {
    double zre = mr,zim = mi;
    COMPLEXSUB(coeff[0],mr,mi);
    for (int j = 1; j <= p; j++) {
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
      COMPLEXADD(coeff[j],zre/j,zim/j);
    }
  }
  else {
    double zre = mr,zim = mi;
    for (int j = 1; j < p; j++) {
      COMPLEXSUB(coeff[j],zre,zim);
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
    }
    COMPLEXSUB(coeff[p],zre,zim);
  }
#endif /* INLINECOMPLEX */
}
/*------------------------------------------------------------------------*/
#ifdef SSE2INIT
inline void mpexp_init_sse(mpexp *This,int p,
                 const double zr,const double zi,const double zr2,const double zi2,
                 const double mr,const double mi,const double mr2,const double mi2,int pot)
/* Initializes a multipole expansion at the lowest level with a single
   pointmass (mr,mi) at (zr,zi). */
{
  /* We have that the potential is given by (pot = 0)
     mp(z) = -m*log(z-zi)
           = -m*log(z-z0-(zi-z0))
           = -m*log((z-z0)*(1-(zi-z0)/(z-z0)))
           = -m*[log(z-z0)+log(1-(zi-z0)/(z-z0))]
           = -m*[log(z-z0)-sum_{j >= 1} (zi-z0)^j/j*(z-z0)^{-j},
     or othwerwise by (pot = 1)
     mp(z) = -m*1/(z-zi)
           = -m*1/(z-z0-(zi-z0))
           = -m/(z-z0)*1/(1-(zi-z0)/(z-z0))
           = -m*sum_{j >= 1} (zi-z0)^{j-1}*(z-z0)^{-j}.
  */

  dcmplx *coeff = This->coeff1;

  if (pot == 0) {
    double re = zr-creal(This->z0),im = zi-cimag(This->z0);
    double zre = mr,zim = mi;
    COMPLEXSUB(coeff[0],mr,mi);
    for (int j = 1; j <= p; j++) {
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
      COMPLEXADD(coeff[j],zre/j,zim/j);
    }
    re = zr2-creal(This->z0),im = zi2-cimag(This->z0);
    zre = mr2;
    zim = mi2;
    COMPLEXSUB(coeff[0],mr2,mi2);
    for (int j = 1; j <= p; j++) {
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
      COMPLEXADD(coeff[j],zre/j,zim/j);
    }
  }
  else {
    __m128d mre=_mm_set_pd(zr-creal(This->z0),zr2-creal(This->z0));
    __m128d mim=_mm_set_pd(zi-cimag(This->z0),zi2-cimag(This->z0));
    __m128d mzre=_mm_set_pd(mr,mr2);
    __m128d mzim=_mm_set_pd(mi,mi2);
    for (int j = 1; j < p; j++) {
#ifdef __SSE3__
      *((__m128d*)&coeff[j])=_mm_sub_pd(*((__m128d*)&coeff[j]),_mm_hadd_pd(mzre,mzim));
#else
      double VC16ALIGN dre[2] GCC16ALIGN;
      double VC16ALIGN dim[2] GCC16ALIGN;
      _mm_store_pd(dre, mzre);
      _mm_store_pd(dim, mzim);
      COMPLEXSUB(coeff[j],dre[0]+dre[1],dim[0]+dim[1]);
#endif
      __m128d mzre0=_mm_sub_pd(_mm_mul_pd(mzre,mre),_mm_mul_pd(mzim,mim));
      mzim=_mm_add_pd(_mm_mul_pd(mzim,mre),_mm_mul_pd(mzre,mim));
      mzre=mzre0;
    }
#ifdef __SSE3__
    *((__m128d*)&coeff[p])=_mm_sub_pd(*((__m128d*)&coeff[p]), _mm_hadd_pd(mzre, mzim));
#else
    double VC16ALIGN dre[2] GCC16ALIGN;
    double VC16ALIGN dim[2] GCC16ALIGN;
    _mm_store_pd(dre, mzre);
    _mm_store_pd(dim, mzim);
    COMPLEXSUB(coeff[p], dre[0]+dre[1], dim[0]+dim[1]);
#endif
  }
}
#endif /*SSE2INIT*/
/*------------------------------------------------------------------------*/
inline void mpexp_initp(mpexp *This,int p,
                 const double zr,const double zi,
                 const double mr,const double mi,int pot)
/* Initializes a multipole expansion at the lowest level with a single
   pointmass (mr,mi) at (zr,zi). The pointmass is assumed to be
   sufficiently far away from the center of the box that a polynomial
   expansion may be used directly. */
{
  /* We have that the potential is given by (pot = 0)
     mp(z) = -m*log(z-zi)
           = -m*log(z0-zi+(z-z0))
           = -m*[log(z0-zi)-sum_{j >= 1} (z-z0)^j/j*(zi-z0)^{-j},
     or othwerwise by (pot = 1)
     mp(z) = -m*1/(z-zi)
           = -m*1/(z0-zi+(z-z0))
           = m*sum_{j >= 1} (z-z0)^{j-1}*(zi-z0)^{-j}.
  */
#ifndef INLINECOMPLEX
  const dcmplx r = 1.0/((zr+I*zi)-This->z0);
  const dcmplx m = mr+I*mi;
  dcmplx *coeff = This->coeff2;

  if (pot == 0) {
    dcmplx z = m*r;
    coeff[0] += m*clog(-r);
    for (int j = 1; j <= p; j++,z *= r) coeff[j] += z/j;
  }
  else {
    dcmplx z = m*r;
    for (int j = 0; j <= p; j++,z *= r) coeff[j] += z;
  }
#else
  const dcmplx r = 1.0/((zr+I*zi)-This->z0);
  const double re = creal(r),im = cimag(r);
  dcmplx *coeff = This->coeff2;

  if (pot == 0) {
    double zre = mr*re-mi*im,zim = mi*re+mr*im;
    {
      const dcmplx logr = clog(-r);
      coeff[0] += (mr*creal(logr)-mi*cimag(logr)+
                   I*(mi*creal(logr)+mr*cimag(logr)));
    }
    for (int j = 1; j <= p; j++) {
      COMPLEXADD(coeff[j],zre/j,zim/j);
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
    }
  }
  else {
    double zre = mr*re-mi*im,zim = mi*re+mr*im;
    for (int j = 0; j < p; j++) {
      COMPLEXADD(coeff[j],zre,zim);
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
    }
    COMPLEXADD(coeff[p],zre,zim);
  }
#endif
}
/*------------------------------------------------------------------------*/
#ifdef SSE2INIT
inline void mpexp_initp_sse(mpexp *This,int p,
                 const double zr,const double zi,const double zr2,const double zi2,
                 const double mr,const double mi,const double mr2,const double mi2,int pot)
/* Initializes a multipole expansion at the lowest level with a single
   pointmass (mr,mi) at (zr,zi). The pointmass is assumed to be
   sufficiently far away from the center of the box that a polynomial
   expansion may be used directly. */
{
  /* We have that the potential is given by (pot = 0)
     mp(z) = -m*log(z-zi)
           = -m*log(z0-zi+(z-z0))
           = -m*[log(z0-zi)-sum_{j >= 1} (z-z0)^j/j*(zi-z0)^{-j},
     or othwerwise by (pot = 1)
     mp(z) = -m*1/(z-zi)
           = -m*1/(z0-zi+(z-z0))
           = m*sum_{j >= 1} (z-z0)^{j-1}*(zi-z0)^{-j}.
  */
  const dcmplx r = 1.0/((zr+I*zi)-This->z0);
  const dcmplx r2 = 1.0/((zr2+I*zi2)-This->z0);
  dcmplx *coeff = This->coeff2;

  if (pot == 0) {
    double re = creal(r),im = cimag(r);
    double zre = mr*re-mi*im,zim = mi*re+mr*im;
    {
      const dcmplx logr = clog(-r);
      coeff[0] += (mr*creal(logr)-mi*cimag(logr)+
                   I*(mi*creal(logr)+mr*cimag(logr)));
    }
    for (int j = 1; j <= p; j++) {
      COMPLEXADD(coeff[j],zre/j,zim/j);
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
    }
    re = creal(r2);
    im = cimag(r2);
    zre = mr2*re-mi2*im;
    zim = mi2*re+mr2*im;
    {
      const dcmplx logr = clog(-r2);
      coeff[0] += (mr2*creal(logr)-mi2*cimag(logr)+
                   I*(mi2*creal(logr)+mr2*cimag(logr)));
    }
    for (int j = 1; j <= p; j++) {
      COMPLEXADD(coeff[j],zre/j,zim/j);
      double zre0 = zre*re-zim*im;
      zim = zim*re+zre*im;
      zre = zre0;
    }
  }
  else {
    __m128d mre=_mm_set_pd(creal(r),creal(r2));
    __m128d mim=_mm_set_pd(cimag(r),cimag(r2));
    __m128d mzre0=_mm_set_pd(mr,mr2);
    __m128d mzim=_mm_set_pd(mi,mi2);
    __m128d mzre=_mm_sub_pd(_mm_mul_pd(mre,mzre0),_mm_mul_pd(mzim,mim));
    mzim=_mm_add_pd(_mm_mul_pd(mzim,mre),_mm_mul_pd(mzre0,mim));
    for (int j = 0; j < p; j++) {
      #ifdef __SSE3__
      *((__m128d*)&coeff[j])=_mm_add_pd(*((__m128d*)&coeff[j]),_mm_hadd_pd(mzre,mzim));
      #else
      double VC16ALIGN dre[2] GCC16ALIGN;
      double VC16ALIGN dim[2] GCC16ALIGN;
      _mm_store_pd(dre, mzre);
      _mm_store_pd(dim, mzim);
      COMPLEXADD(coeff[j],dre[0]+dre[1],dim[0]+dim[1]);
      #endif
      mzre0=_mm_sub_pd(_mm_mul_pd(mzre,mre),_mm_mul_pd(mzim,mim));
      mzim=_mm_add_pd(_mm_mul_pd(mzim,mre),_mm_mul_pd(mzre,mim));
      mzre=mzre0;
    }
    #ifdef __SSE3__
      *((__m128d*)&coeff[p])=_mm_add_pd(*((__m128d*)&coeff[p]),_mm_hadd_pd(mzre,mzim));
    #else
    double VC16ALIGN dre[2] GCC16ALIGN;
    double VC16ALIGN dim[2] GCC16ALIGN;
    _mm_store_pd(dre, mzre);
    _mm_store_pd(dim, mzim);
    COMPLEXADD(coeff[p], dre[0]+dre[1], dim[0]+dim[1]);
    #endif
  }
}
#endif /*SSE2INIT*/
/*------------------------------------------------------------------------*/
bool mpexp_theta(const mpexp *This, const mpexp *that,double cutoff)
/* Theta criterion for accurate cluster to cluster interaction between
   two boxes This and that.

   This routine is a workhorse and needs to be fast.

   The maximum relative error in using a multipole expansion from a
   box 'this' to evaluate the field in a box 'that' is error <=
   constant * (R/|R+minimum distance between the boxes|)^(p+1) for
   points |z| > R. Using that the minimum distance is >= d-R-r (where
   d is the center-to-center distance and where r is the radius of the
   smaller box), we are led to the requirement that R/(d-r) <= theta
   and this is what is implemented below (see [4] for the full
   analysis).

   This also means that theta = 1/2 is equivalent to the boxes being
   classically well-separated as in [1,2]. With a general theta, any
   of the two sets may be expanded by a factor of 1/theta and
   arbitrarily rotated about its center point without touching the
   other set.

   However, in that case one should note that the theta used in the
   classical uniform multipole method [2] is actually 1/(2*sqrt(2)-1)
   = 0.5469181 and not 1/2 -- perhaps for convenience, given the
   quadratic boxes.

   Also, in [3] a different and less stringent criterion is used; R <=
   c*d. For uniform boxes, their c = 1/2 would mean theta = sqrt(2)/3
   = 0.4714045 here. Note that trivially, R <= theta/(1+theta)*d
   implies R+theta*r <= theta*d, but the former criterion is less
   sharp when the two boxes have different sizes.

   A slightly sharper result stems from squaring and using Young's
   weighted inequality; then R^2+theta*r^2 <= theta^2/(1+theta)*d^2
   implies R+theta*r <= theta*d. This is used below as a first check
   in order to avoid computing the square-root. Note: not used in the
   current implementation.

   Also, note that the GKZ-criterion in [3] yields an error
   proportional to (c/(1-2*c))^(p+1) and not to c^(p+1) as they
   claim. Hence c = theta/(1+2*theta) would be appropriate -- this is
   even less sharp compared to the trivial interpretation mentioned
   above. Unfortunately, following their method of proof, the
   currently implemented criterion yields an error proportional to
   (theta/(1-theta))^(p+1), which would suggest that the choice theta
   = 1/2 should be avoided.

   Notes: the GKZ-criterion is found in [3], Eq. (8.31) at p. 379 and
   the error estimate is derived at p. 386. The first estimate [3] in
   (8.40) is incorrect (change (1+2*theta) to (1-2*theta) and the
   second is unsharp (replace (1-2*theta) with 1/(1+2*theta)). In the
   proof that follows, the problematic factor (1-2*theta) is
   inadvertently omitted.

   References:
     [1] J. Carrier, L. Greengard, and V. Rokhlin: "A fast adaptive
         multipole algorithm for particle simulations". SIAM
         J. Sci. Stat. Comput., 9(4):669--686, 1988.
     [2] L. Greengard and V. Rokhlin: "A fast algorithm for particle
         simulations". J. Comput. Phys., 73:325--348, 1987.
     [3] M. Griebel, S. Knapek, and G. Zumbusch: "Numerical Simulation
         in Molecular Dynamics", volume 5 of Texts in Computational
         Science and Engineering. Springer Verlag, Berlin, 2007.
     [4] S. Engblom: "On well-separated sets and fast mutlipole
         methods", Appl. Numer. Math. 61(10):1096--1102, 2011.
*/
{
  // *** (almost) sqrt-free implementation in the C99-version
#ifdef RADIALSHRINK
//  const double d = cabs(This->z0-that->z0);
//  const double d = hypot(creal(This->z0)-creal(that->z0),cimag(This->z0)-cimag(that->z0));
  const double d = sqrt((creal(This->z0)-creal(that->z0))*(creal(This->z0)-creal(that->z0))+(cimag(This->z0)-cimag(that->z0))*(cimag(This->z0)-cimag(that->z0)));

  /* note: criterion uses in fact a strict inequality to get rid of
     the singular case when both boxes are the same points (This can
     happen!) */

#ifdef THETACUTOFFCHECK
  if (d-This->absd0-that->absd0 < cutoff)
    return false;
#endif
  if (This->absd0 >= that->absd0)
    return This->absd0+theta*that->absd0 < theta*d;
  else
    return that->absd0+theta*This->absd0 < theta*d;

#else /* !RADIALSHRINK */

  const double Thisrad = cabs(This->d0),thatrad = cabs(that->d0);
  const double d = cabs(This->z0-that->z0);

#ifdef THETACUTOFFCHECK
  if (d-Thisrad-thatrad < cutoff)
    return false;
#endif
  if (Thisrad >= thatrad)
    return Thisrad+theta*thatrad < theta*d;
  else
    return thatrad+theta*Thisrad < theta*d;
#endif /* RADIALSHRINK */
}
/*------------------------------------------------------------------------*/
int mpexp_theta2(const mpexp *This,const mpexp *that,double cutoff)
/* 2nd theta criterion: if mpexp_theta() above is false, is it due to
   only one of the boxes? Then the multipoles of the larger box can be
   shifted into a polynomial expansion in the smaller, and the
   multipole expansion in the smaller can be evaluated at each point
   in the larger box. Returns -1 if this is the smaller box (cluster
   to particle) and 1 otherwise (particle to cluster). */
{
  // derived as a kind of adjoint to mpexp_theta() above
#ifdef RADIALSHRINK
//  const double d = cabs(This->z0-that->z0);
//  const double d = hypot(creal(This->z0)-creal(that->z0),cimag(This->z0)-cimag(that->z0));
  const double d = sqrt((creal(This->z0)-creal(that->z0))*(creal(This->z0)-creal(that->z0))+(cimag(This->z0)-cimag(that->z0))*(cimag(This->z0)-cimag(that->z0)));

#ifdef THETACUTOFFCHECK
  if (d-This->absd0-that->absd0 < cutoff)
    return 0;
#endif
  if (This->absd0 >= that->absd0)
    return that->absd0+theta*This->absd0 < theta*d ? -1 : 0;
  else
    return This->absd0+theta*that->absd0 < theta*d ? 1 : 0;

#else /* !RADIALSHRINK */

  const double Thisrad = cabs(This->d0),thatrad = cabs(that->d0);
  const double d = cabs(This->z0-that->z0);
#ifdef THETACUTOFFCHECK
  if (d-Thisrad-thatrad < cutoff)
    return 0;
#endif
  if (Thisrad >= thatrad)
    return thatrad+theta*Thisrad < theta*d ? -1 : 0;
  else
    return Thisrad+theta*thatrad < theta*d ? 1 : 0;
#endif /* RADIALSHRINK */
}
/*------------------------------------------------------------------------*/
void mpexp_eval(const mpexp *This,int p,
                double *vr,double *vi,
                const double *zr,const double *zi,
                const int *ix,int begin,int end)
/* Evaluation of near field at coordinates (zr,zi). */
{
  const dcmplx *coeff = This->coeff2;

  if (vi != NULL)
    for (int i = begin; i < end; i++) {
#ifndef INLINECOMPLEX
      dcmplx polyval = coeff[p];
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);

      // Horner's method
      for (int j = p-1; j >= 0; j--) polyval = r*polyval+coeff[j];

      // add the result
      vr[i] += creal(polyval);
      vi[i] += cimag(polyval);
#else
      const double re = zr[i]-creal(This->z0),im = zi[i]-cimag(This->z0);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method
      for (int j = p-1; j >= 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }

      // add the result
      vr[i] += pre;
      vi[i] += pim;
#endif
    }
  else
    for (int i = begin; i < end; i++) {
#ifndef INLINECOMPLEX
      dcmplx polyval = coeff[p];
      const dcmplx r = (zr[i]+I*zi[i])-This->z0;
      for (int j = p-1; j >= 0; j--) polyval = r*polyval+coeff[j];
      vr[i] += creal(polyval);
#else
      const double re = zr[i]-creal(This->z0),im = zi[i]-cimag(This->z0);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);
      for (int j = p-1; j >= 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      vr[i] += pre;
#endif
    }
}
/*------------------------------------------------------------------------*/
void mpexp_evalmp(const mpexp *restrict This,int p,int pot,
                  double *restrict vr,double *restrict vi,
                  const double *restrict zr,const double *restrict zi,
                  const int *restrict ix,int begin,int end)
/* Evaluation of multipole expansion at coordinates (zr,zi). */
{
  const dcmplx *coeff = This->coeff1;

  if (vi != NULL) {
    for (int i = begin; i < end; i++) {
#ifndef INLINECOMPLEX
      dcmplx polyval = coeff[p];
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) polyval = r*polyval+coeff[j];
      if (pot == 0)
        polyval = r*polyval-coeff[0]*clog(r);
      else
        polyval *= r;

      vr[i] += creal(polyval);
      vi[i] += cimag(polyval);

#else
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const double re = creal(r),im = cimag(r);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      {
        double pre0;
        pre0 = re*pre-im*pim;
        pim = re*pim+im*pre;
        pre = pre0;
        if (pot == 0) {
          const dcmplx logr = clog(r);
          pre -= creal(coeff[0])*creal(logr)-cimag(coeff[0])*cimag(logr);
          pim -= creal(coeff[0])*cimag(logr)+cimag(coeff[0])*creal(logr);
        }
      }

      vr[i] += pre;
      vi[i] += pim;
#endif
    }
  }
  else {
    for (int i = begin; i < end; i++) {
#ifndef INLINECOMPLEX
      dcmplx polyval = coeff[p];
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      for (int j = p-1; j > 0; j--) polyval = r*polyval+coeff[j];
      if (pot == 0)
        polyval = r*polyval-coeff[0]*clog(r);
      else
        polyval *= r;
      vr[i] += creal(polyval);
#else
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const double re = creal(r),im = cimag(r);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      {
        double pre0;
        pre0 = re*pre-im*pim;
        pim = re*pim+im*pre;
        pre = pre0;
        if (pot == 0) {
          const dcmplx logr = clog(r);
          pre -= creal(coeff[0])*creal(logr)-cimag(coeff[0])*cimag(logr);
          pim -= creal(coeff[0])*cimag(logr)+cimag(coeff[0])*creal(logr);
        }
      }

      vr[i] += pre;
#endif
    }
  }
}
#ifdef SSE2DIRECT
/*------------------------------------------------------------------------*/
void mpexp_evalmp_sse(const mpexp *restrict This,int p,int pot,
                  double *restrict vr,double *restrict vi,
                  const double *restrict zr,const double *restrict zi,
                  const int *restrict ix,int begin,int end)
/* Evaluation of multipole expansion at coordinates (zr,zi). */
{
  const dcmplx *coeff = This->coeff1;

  if (vi != NULL) {
    int i = begin;
    if (i < end && i&1) {
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const double re = creal(r),im = cimag(r);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      {
        double pre0;
        pre0 = re*pre-im*pim;
        pim = re*pim+im*pre;
        pre = pre0;
      }
      vr[i] += pre;
      vi[i] += pim;
      i++;
    }
    #ifdef SSEEVALUNROLL
    for (; i < end-3; i+=4) {
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const dcmplx r2 = 1.0/((zr[i+1]+I*zi[i+1])-This->z0);
      __m128d mre=_mm_set_pd(creal(r2),creal(r));
      __m128d mim=_mm_set_pd(cimag(r2),cimag(r));
      const dcmplx r3 = 1.0/((zr[i+2]+I*zi[i+2])-This->z0);
      const dcmplx r4 = 1.0/((zr[i+3]+I*zi[i+3])-This->z0);
      __m128d mre2=_mm_set_pd(creal(r4),creal(r3));
      __m128d mim2=_mm_set_pd(cimag(r4),cimag(r3));
      __m128d mpre=_mm_set1_pd(creal(coeff[p]));
      __m128d mpim=_mm_set1_pd(cimag(coeff[p]));
      __m128d mpre2=mpre;
      __m128d mpim2=mpim;
//      const double re = creal(r),im = cimag(r);
//      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        __m128d rcoeff=_mm_set1_pd(creal(coeff[j]));
        __m128d icoeff=_mm_set1_pd(cimag(coeff[j]));
        __m128d mpre0=_mm_add_pd(rcoeff,_mm_sub_pd(_mm_mul_pd(mpre,mre),_mm_mul_pd(mpim,mim)));
        mpim=_mm_add_pd(icoeff,_mm_add_pd(_mm_mul_pd(mpim,mre),_mm_mul_pd(mpre,mim)));
        mpre=mpre0;
        mpre0=_mm_add_pd(rcoeff,_mm_sub_pd(_mm_mul_pd(mpre2,mre2),_mm_mul_pd(mpim2,mim2)));
        mpim2=_mm_add_pd(icoeff,_mm_add_pd(_mm_mul_pd(mpim2,mre2),_mm_mul_pd(mpre2,mim2)));
        mpre2=mpre0;
//        double pre0;
//        pre0 = re*pre-im*pim+creal(coeff[j]);
//        pim = re*pim+im*pre+cimag(coeff[j]);
//        pre = pre0;
      }
      {
        __m128d mpre0=_mm_sub_pd(_mm_mul_pd(mpre,mre),_mm_mul_pd(mpim,mim));
        mpim=_mm_add_pd(_mm_mul_pd(mpim,mre),_mm_mul_pd(mpre,mim));
        mpre=mpre0;
        mpre0=_mm_sub_pd(_mm_mul_pd(mpre2,mre2),_mm_mul_pd(mpim2,mim2));
        mpim2=_mm_add_pd(_mm_mul_pd(mpim2,mre2),_mm_mul_pd(mpre2,mim2));
        mpre2=mpre0;
      }
      __m128d* mvr=(__m128d*)&vr[i];
      __m128d* mvi=(__m128d*)&vi[i];
      mvr[0]=_mm_add_pd(mvr[0],mpre);
      mvi[0]=_mm_add_pd(mvi[0],mpim);
      mvr[1]=_mm_add_pd(mvr[1],mpre2);
      mvi[1]=_mm_add_pd(mvi[1],mpim2);
    }
    #endif /*SSEEVALUNROLL*/
    for (; i < end-1; i+=2) {
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const dcmplx r2 = 1.0/((zr[i+1]+I*zi[i+1])-This->z0);
      __m128d mre=_mm_set_pd(creal(r2),creal(r));
      __m128d mim=_mm_set_pd(cimag(r2),cimag(r));
      __m128d mpre=_mm_set1_pd(creal(coeff[p]));
      __m128d mpim=_mm_set1_pd(cimag(coeff[p]));
//      const double re = creal(r),im = cimag(r);
//      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        __m128d rcoeff=_mm_set1_pd(creal(coeff[j]));
        __m128d icoeff=_mm_set1_pd(cimag(coeff[j]));
        __m128d mpre0=_mm_add_pd(rcoeff,_mm_sub_pd(_mm_mul_pd(mpre,mre),_mm_mul_pd(mpim,mim)));
        mpim=_mm_add_pd(icoeff,_mm_add_pd(_mm_mul_pd(mpim,mre),_mm_mul_pd(mpre,mim)));
        mpre=mpre0;
//        double pre0;
//        pre0 = re*pre-im*pim+creal(coeff[j]);
//        pim = re*pim+im*pre+cimag(coeff[j]);
//        pre = pre0;
      }
      {
        __m128d mpre0=_mm_sub_pd(_mm_mul_pd(mpre,mre),_mm_mul_pd(mpim,mim));
        mpim=_mm_add_pd(_mm_mul_pd(mpim,mre),_mm_mul_pd(mpre,mim));
        mpre=mpre0;
      }
      __m128d* mvr=(__m128d*)&vr[i];
      __m128d* mvi=(__m128d*)&vi[i];
      mvr[0]=_mm_add_pd(mvr[0],mpre);
      mvi[0]=_mm_add_pd(mvi[0],mpim);
    }
    if (i < end) {
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const double re = creal(r),im = cimag(r);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      {
        double pre0;
        pre0 = re*pre-im*pim;
        pim = re*pim+im*pre;
        pre = pre0;
      }
      vr[i] += pre;
      vi[i] += pim;
    }
  }
  else {
    for (int i = begin; i < end; i++) {
#ifndef INLINECOMPLEX
      dcmplx polyval = coeff[p];
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      for (int j = p-1; j > 0; j--) polyval = r*polyval+coeff[j];
      if (pot == 0)
        polyval = r*polyval-coeff[0]*clog(r);
      else
        polyval *= r;
//       vr[ix[i]] += creal(polyval);
      vr[i] += creal(polyval);
#else
      const dcmplx r = 1.0/((zr[i]+I*zi[i])-This->z0);
      const double re = creal(r),im = cimag(r);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method for inverse polynomial
      for (int j = p-1; j > 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      {
        double pre0;
        pre0 = re*pre-im*pim;
        pim = re*pim+im*pre;
        pre = pre0;
        if (pot == 0) {
          const dcmplx logr = clog(r);
          pre -= creal(coeff[0])*creal(logr)-cimag(coeff[0])*cimag(logr);
          pim -= creal(coeff[0])*cimag(logr)+cimag(coeff[0])*creal(logr);
        }
      }

      vr[i] += pre;
#endif
    }
  }
}
/*------------------------------------------------------------------------*/
void mpexp_eval_sse(const mpexp *This,int p,
                double *vr,double *vi,
                const double *zr,const double *zi,
                const int *ix,int begin,int end)
/* Evaluation of near field at coordinates (zr,zi). */
{
  const dcmplx *coeff = This->coeff2;

  if (vi != NULL) {
    int i=begin;
    if (i < end && i&1) {
      const double re = zr[i]-creal(This->z0),im = zi[i]-cimag(This->z0);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method
      for (int j = p-1; j >= 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }

      // add the result
      vr[i] += pre;
      vi[i] += pim;
      i++;
    }
    __m128d z0re=_mm_set1_pd(creal(This->z0));
    __m128d z0im=_mm_set1_pd(cimag(This->z0));
    #ifdef SSEEVALUNROLL
    for (;i < end-3; i+=4) {
      __m128d mre=_mm_load_pd(&zr[i]);
      __m128d mim=_mm_load_pd(&zi[i]);
      mre=_mm_sub_pd(mre,z0re);
      mim=_mm_sub_pd(mim,z0im);
      __m128d mre2=_mm_load_pd(&zr[i+2]);
      __m128d mim2=_mm_load_pd(&zi[i+2]);
      mre2=_mm_sub_pd(mre2,z0re);
      mim2=_mm_sub_pd(mim2,z0im);
//      const double re = zr[i]-creal(This->z0),im = zi[i]-cimag(This->z0);
      __m128d mpre=_mm_set1_pd(creal(coeff[p]));
      __m128d mpim=_mm_set1_pd(cimag(coeff[p]));
      __m128d mpre2=mpre;
      __m128d mpim2=mpim;
//      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method
      for (int j = p-1; j >= 0; j--) {
        __m128d rcoeff=_mm_set1_pd(creal(coeff[j]));
        __m128d icoeff=_mm_set1_pd(cimag(coeff[j]));
        __m128d mpre0=_mm_add_pd(rcoeff,_mm_sub_pd(_mm_mul_pd(mpre,mre),_mm_mul_pd(mpim,mim)));
        mpim=_mm_add_pd(icoeff,_mm_add_pd(_mm_mul_pd(mpim,mre),_mm_mul_pd(mpre,mim)));
        mpre=mpre0;
        mpre0=_mm_add_pd(rcoeff,_mm_sub_pd(_mm_mul_pd(mpre2,mre2),_mm_mul_pd(mpim2,mim2)));
        mpim2=_mm_add_pd(icoeff,_mm_add_pd(_mm_mul_pd(mpim2,mre2),_mm_mul_pd(mpre2,mim2)));
        mpre2=mpre0;
      }
      __m128d* mvr=(__m128d*)&vr[i];
      __m128d* mvi=(__m128d*)&vi[i];
      mvr[0]=_mm_add_pd(mvr[0],mpre);
      mvi[0]=_mm_add_pd(mvi[0],mpim);
      mvr[1]=_mm_add_pd(mvr[1],mpre2);
      mvi[1]=_mm_add_pd(mvi[1],mpim2);
    }
    #endif /*SSEEVALUNROLL*/
    for (;i < end-1; i+=2) {
      __m128d mre=_mm_load_pd(&zr[i]);
      __m128d mim=_mm_load_pd(&zi[i]);
      mre=_mm_sub_pd(mre,z0re);
      mim=_mm_sub_pd(mim,z0im);
      __m128d mpre=_mm_set1_pd(creal(coeff[p]));
      __m128d mpim=_mm_set1_pd(cimag(coeff[p]));

      // Horner's method
      for (int j = p-1; j >= 0; j--) {
        __m128d rcoeff=_mm_set1_pd(creal(coeff[j]));
        __m128d icoeff=_mm_set1_pd(cimag(coeff[j]));
        __m128d mpre0=_mm_add_pd(rcoeff,_mm_sub_pd(_mm_mul_pd(mpre,mre),_mm_mul_pd(mpim,mim)));
        mpim=_mm_add_pd(icoeff,_mm_add_pd(_mm_mul_pd(mpim,mre),_mm_mul_pd(mpre,mim)));
        mpre=mpre0;
      }
      __m128d* mvr=(__m128d*)&vr[i];
      __m128d* mvi=(__m128d*)&vi[i];
      mvr[0]=_mm_add_pd(mvr[0],mpre);
      mvi[0]=_mm_add_pd(mvi[0],mpim);
    }
    if (i < end) {
      const double re = zr[i]-creal(This->z0),im = zi[i]-cimag(This->z0);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);

      // Horner's method
      for (int j = p-1; j >= 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }

      // add the result
      vr[i] += pre;
      vi[i] += pim;
      i++;
    }
  }
  else {
    for (int i = begin; i < end; i++) {
      const double re = zr[i]-creal(This->z0),im = zi[i]-cimag(This->z0);
      double pre = creal(coeff[p]),pim = cimag(coeff[p]);
      for (int j = p-1; j >= 0; j--) {
        double pre0;
        pre0 = re*pre-im*pim+creal(coeff[j]);
        pim = re*pim+im*pre+cimag(coeff[j]);
        pre = pre0;
      }
      vr[i] += pre;
    }
  }
}
#endif
/*------------------------------------------------------------------------*/
#define BEGIN(begin, i) (begin == -1 ? i+1 : begin)
//note: SSE support only implemented for harmonic potential with real masses
void mpexp_directInteract2(double *restrict pr,double *restrict pi,
                           const double *restrict zr,const double *restrict zi,
                           const double *restrict mr,const double *restrict mi,
                           const int *restrict ix,
                           int begin1,int end1,int begin2,int end2,
                           int pot,SMOOTHER smooth,
                           double cutoff,double shape,double scale)
/* Direct interaction between potentials inside a box (begin2 == -1)
   or between two boxes. */
{
  if (pot == 0) { // logaritmic potential
    if (mi == NULL) { // rlog2()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++)
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 0.5*log(re*re+im*im);
#else
            const double P = log(hypot(zr[i]-zr[j],zi[i]-zi[j]));
#endif
            pr[i] -= mr[j]*P;
            pr[j] -= mr[i]*P;
          }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++) {
          if (begin2 == -1)
            pr[i] -= mr[i]*shape;
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*(ad/cutoff)+shape;
#else
            const double ad = hypot(zr[i]-zr[j],zi[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*(ad/cutoff)*(ad/cutoff)+shape;
#endif
            pr[i] -= mr[j]*P;
            pr[j] -= mr[i]*P;
          }
        }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++) {
          if (begin2 == -1)
            pr[i] -= mr[i]*0.5*log(shape)*scale;
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*log(shape+ad)*scale;
#else
            const double ad = hypot(zr[i]-zr[j],zi[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*log(shape+ad*ad)*scale;
#endif
            pr[i] -= mr[j]*P;
            pr[j] -= mr[i]*P;
          }
        }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++) {
          if (begin2 == -1)
            pr[i] += mr[i]*0.5*(log(shape)+gamma_)*scale;
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else if (ad != 0.0)
              P = 0.5*(log(ad)+expint(shape*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#else
            const double ad = hypot(zr[i]-zr[j],zi[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else if (ad != 0.0)
              P = (log(ad)+0.5*expint(shape*ad*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#endif
            pr[i] -= mr[j]*P;
            pr[j] -= mr[i]*P;
          }
        }
        break;
      }
    }
    else { // zlog2()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++)
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 0.5*log(re*re+im*im);
#else
            const double P = log(hypot(zr[i]-zr[j],
                                       zi[i]-zi[j]));
#endif
            pr[i] -= mr[j]*P;
            pi[i] -= mi[j]*P;
            pr[j] -= mr[i]*P;
            pi[j] -= mi[i]*P;
          }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++) {
          if (begin2 == -1) {
            pr[i] -= mr[i]*shape;
            pi[i] -= mi[i]*shape;
          }
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*(ad/cutoff)+shape;
#else
            const double ad = hypot(zr[i]-zr[j],zi[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*(ad/cutoff)*(ad/cutoff)+shape;
#endif
            pr[i] -= mr[j]*P;
            pi[i] -= mi[j]*P;
            pr[j] -= mr[i]*P;
            pi[j] -= mi[i]*P;
          }
        }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++) {
          if (begin2 == -1) {
            pr[i] -= mr[i]*0.5*log(shape)*scale;
            pi[i] -= mi[i]*0.5*log(shape)*scale;
          }
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*log(shape+ad)*scale;
#else
            const double ad = hypot(zr[i]-zr[j],zi[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*log(shape+ad*ad)*scale;
#endif
            pr[i] -= mr[j]*P;
            pi[i] -= mi[j]*P;
            pr[j] -= mr[i]*P;
            pi[j] -= mi[i]*P;
          }
        }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++) {
          if (begin2 == -1) {
            pr[i] += mr[i]*0.5*(log(shape)+gamma_)*scale;
            pi[i] += mi[i]*0.5*(log(shape)+gamma_)*scale;
          }
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else if (ad != 0.0)
              P = 0.5*(log(ad)+expint(shape*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#else
            const double ad = hypot(zr[i]-zr[j],zi[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else if (ad != 0.0)
              P = (log(ad)+0.5*expint(shape*ad*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#endif
            pr[i] -= mr[j]*P;
            pi[i] -= mi[j]*P;
            pr[j] -= mr[i]*P;
            pi[j] -= mi[i]*P;
          }
        }
        break;
      }
    }
  }
  else { // harmonic potential
    if (mi == NULL) { // rinv2()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 1.0/(re*re+im*im);
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d P=_mm_add_pd(_mm_mul_pd(mre, mre), _mm_mul_pd(mim, mim));
            P=_mm_div_pd(_mm_set1_pd(1.0), P);
            __m128d mmr=_mm_load_pd(&mr[j]);
            mmr=_mm_mul_pd(mmr, P);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 1.0/(re*re+im*im);
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 1.0/(re*re+im*im);
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
#else
            const dcmplx P = 1.0/(zr[i]-zr[j]+I*(zi[i]-zi[j]));
            pr[i] -= mr[j]*creal(P);
            pi[i] -= mr[j]*cimag(P);
            pr[j] += mr[i]*creal(P);
            pi[j] += mr[i]*cimag(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadl_pd(P, &shape);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadh_pd(P, &shape);
            }
            mmr=_mm_load_pd(&mr[j]);
            mmr=_mm_mul_pd(mmr, P);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
#else
            const dcmplx d = zr[i]-zr[j]+I*(zi[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)*shape;
            pr[i] -= mr[j]*creal(P);
            pi[i] -= mr[j]*cimag(P);
            pr[j] += mr[i]*creal(P);
            pi[j] += mr[i]*cimag(P);
#endif
          }
#endif
        }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++) {
          #ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadl_pd(P, &Pd);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadh_pd(P, &Pd);
            }
            mmr=_mm_load_pd(&mr[j]);
            mmr=_mm_mul_pd(mmr, P);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            pr[i] -= mr[j]*P*re;
            pi[i] += mr[j]*P*im;
            pr[j] += mr[i]*P*re;
            pi[j] -= mr[i]*P*im;
#else
            const dcmplx d = zr[i]-zr[j]+I*(zi[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)/(shape+ad*ad)*scale;
            pr[i] -= mr[j]*creal(P);
            pi[i] -= mr[j]*cimag(P);
            pr[j] += mr[i]*creal(P);
            pi[j] += mr[i]*cimag(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              pr[i] -= mr[j]*P*re;
              pi[i] += mr[j]*P*im;
              pr[j] += mr[i]*P*re;
              pi[j] -= mr[i]*P*im;
            }
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
#ifdef OSEENSSEALTERNATIVE
            double VC16ALIGN ad[2],Pd[2] GCC16ALIGN;
            _mm_store_pd(ad,mmr);
            if(ad[0]<cutoff||ad[1]<cutoff) {
              _mm_store_pd(Pd,P);
              if(ad[0]<cutoff) {
                if(ad[0]!=0.0)
                  Pd[0]*=-expm1(-shape*ad[0])*scale;
                else
                  Pd[0]=0;
              }
              if(ad[1]<cutoff) {
                if(ad[1]!=0.0)
                  Pd[1]*=-expm1(-shape*ad[1])*scale;
                else
                  Pd[1]=0;
              }
              P=_mm_load_pd(Pd);
            }
#else
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storel_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadl_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadl_pd(P, &Pd);
              }
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storeh_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadh_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadh_pd(P, &Pd);
              }
            }
#endif
            mmr=_mm_load_pd(&mr[j]);
            mmr=_mm_mul_pd(mmr, P);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              pr[i] -= mr[j]*P*re;
              pi[i] += mr[j]*P*im;
              pr[j] += mr[i]*P*re;
              pi[j] -= mr[i]*P*im;
            }
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              pr[i] -= mr[j]*P*re;
              pi[i] += mr[j]*P*im;
              pr[j] += mr[i]*P*re;
              pi[j] -= mr[i]*P*im;
            }
#else
            const dcmplx d = zr[i]-zr[j]+I*(zi[i]-zi[j]);
            const double ad = cabs(d);
            if (ad != 0.0) {
              dcmplx P;
              if (ad >= cutoff)
                P = 1.0/d;
              else
                P = -expm1(-shape*ad*ad)/d*scale;
              pr[i] -= mr[j]*creal(P);
              pi[i] -= mr[j]*cimag(P);
              pr[j] += mr[i]*creal(P);
              pi[j] += mr[i]*cimag(P);
            }
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;
      }
    }
    else { // zinv2
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 1.0/(re*re+im*im);
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d mmii=_mm_set1_pd(mi[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d P=_mm_add_pd(_mm_mul_pd(mre, mre), _mm_mul_pd(mim, mim));
            P=_mm_div_pd(_mm_set1_pd(1.0), P);
            __m128d mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            __m128d P2=_mm_mul_pd(P, mmii);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P2, mim));
            pij[0]=_mm_sub_pd(pij[0], _mm_mul_pd(P2, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 1.0/(re*re+im*im);
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double P = 1.0/(re*re+im*im);
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
#else
            const dcmplx P = 1.0/(zr[i]-zr[j]+I*(zi[i]-zi[j]));
            pr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
            pi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
            pr[j] += mr[i]*creal(P)-mi[i]*cimag(P);
            pi[j] += mr[i]*cimag(P)+mi[i]*creal(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d mmii=_mm_set1_pd(mi[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadl_pd(P, &shape);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadh_pd(P, &shape);
            }
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            __m128d P2=_mm_mul_pd(P, mmii);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P2, mim));
            pij[0]=_mm_sub_pd(pij[0], _mm_mul_pd(P2, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
#else
            const dcmplx d = zr[i]-zr[j]+I*(zi[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)*shape;
            pr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
            pi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
            pr[j] += mr[i]*creal(P)-mi[i]*cimag(P);
            pi[j] += mr[i]*cimag(P)+mi[i]*creal(P);
#endif
          }
#endif
        }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d mmii=_mm_set1_pd(mi[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadl_pd(P, &Pd);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadh_pd(P, &Pd);
            }
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            __m128d P2=_mm_mul_pd(P, mmii);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P2, mim));
            pij[0]=_mm_sub_pd(pij[0], _mm_mul_pd(P2, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            pr[i] -= (mr[j]*re+mi[j]*im)*P;
            pi[i] += (mr[j]*im-mi[j]*re)*P;
            pr[j] += (mr[i]*re+mi[i]*im)*P;
            pi[j] -= (mr[i]*im-mi[i]*re)*P;
#else
            const dcmplx d = zr[i]-zr[j]+I*(zi[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)/(shape+ad*ad)*scale;
            pr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
            pi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
            pr[j] += mr[i]*creal(P)-mi[i]*cimag(P);
            pi[j] += mr[i]*cimag(P)+mi[i]*creal(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = BEGIN(begin2,i);
          if (j < end2&&j&1) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              pr[i] -= (mr[j]*re+mi[j]*im)*P;
              pi[i] += (mr[j]*im-mi[j]*re)*P;
              pr[j] += (mr[i]*re+mi[i]*im)*P;
              pi[j] -= (mr[i]*im-mi[i]*re)*P;
            }
            j++;
          }
          __m128d mzr=_mm_set1_pd(zr[i]);
          __m128d mzi=_mm_set1_pd(zi[i]);
          __m128d mmi=_mm_set1_pd(mr[i]);
          __m128d mmii=_mm_set1_pd(mi[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            __m128d* prj=(__m128d*)&pr[j];
            __m128d* pij=(__m128d*)&pi[j];
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
#ifdef OSEENSSEALTERNATIVE
            double VC16ALIGN ad[2],Pd[2] GCC16ALIGN;
            _mm_store_pd(ad,mmr);
            if(ad[0]<cutoff||ad[1]<cutoff) {
              _mm_store_pd(Pd,P);
              if(ad[0]<cutoff) {
                if(ad[0]!=0.0)
                  Pd[0]*=-expm1(-shape*ad[0])*scale;
                else
                  Pd[0]=0;
              }
              if(ad[1]<cutoff) {
                if(ad[1]!=0.0)
                  Pd[1]*=-expm1(-shape*ad[1])*scale;
                else
                  Pd[1]=0;
              }
              P=_mm_load_pd(Pd);
            }
#else
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storel_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadl_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadl_pd(P, &Pd);
              }
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storeh_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadh_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadh_pd(P, &Pd);
              }
            }
#endif
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            __m128d P2=_mm_mul_pd(P, mmii);
            P=_mm_mul_pd(P, mmi);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P, mre));
            pij[0]=_mm_add_pd(pij[0], _mm_mul_pd(P, mim));
            prj[0]=_mm_sub_pd(prj[0], _mm_mul_pd(P2, mim));
            pij[0]=_mm_sub_pd(pij[0], _mm_mul_pd(P2, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          pr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          pi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              pr[i] -= (mr[j]*re+mi[j]*im)*P;
              pi[i] += (mr[j]*im-mi[j]*re)*P;
              pr[j] += (mr[i]*re+mi[i]*im)*P;
              pi[j] -= (mr[i]*im-mi[i]*re)*P;
            }
          }
#else /*SSE2DIRECT*/
          for (int j = BEGIN(begin2,i); j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = zr[i]-zr[j],im = zi[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              pr[i] -= (mr[j]*re+mi[j]*im)*P;
              pi[i] += (mr[j]*im-mi[j]*re)*P;
              pr[j] += (mr[i]*re+mi[i]*im)*P;
              pi[j] -= (mr[i]*im-mi[i]*re)*P;
            }
#else
            const dcmplx d = zr[i]-zr[j]+I*(zi[i]-zi[j]);
            const double ad = cabs(d);
            if (ad != 0.0) {
              dcmplx P;
              if (ad >= cutoff)
                P = 1.0/d;
              else
                P = -expm1(-shape*ad*ad)/d*scale;
              pr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
              pi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
              pr[j] += mr[i]*creal(P)-mi[i]*cimag(P);
              pi[j] += mr[i]*cimag(P)+mi[i]*creal(P);
            }
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;
      }
    }
  }
}
/*------------------------------------------------------------------------*/
void mpexp_directInteract(double *restrict qr,double *restrict qi,
                          const double *restrict er,const double *restrict ei,
                          const int *restrict jx,int begin1,int end1,
                          const double *restrict zr,const double *restrict zi,
                          const double *restrict mr,const double *restrict mi,
                          const int *restrict ix,int begin2,int end2,
                          int pot,SMOOTHER smooth,
                          double cutoff,double shape,double scale)
/* Direct evaluation of field from potentials inside one box computed
   at points (er,ei) inside another box. */
{
  if (pot == 0) { // logarithmic potential
    if (mi == NULL) { // rlog()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++)
            if (er[i] != zr[j] || ei[i] != zi[j]) {
#ifdef INLINECOMPLEX_DIRECT
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 0.5*log(re*re+im*im);
#else
              const double P = log(hypot(er[i]-zr[j],ei[i]-zi[j]));
#endif
              qr[i] -= mr[j]*P;
            }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*(ad/cutoff)+shape;
#else
            const double ad = hypot(er[i]-zr[j],ei[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*(ad/cutoff)*(ad/cutoff)+shape;
#endif
            qr[i] -= mr[j]*P;
          }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*log(shape+ad)*scale;
#else
            const double ad = hypot(er[i]-zr[j],ei[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*log(shape+ad*ad)*scale;
#endif
            qr[i] -= mr[j]*P;
          }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else if (ad != 0.0)
              P = 0.5*(log(ad)+expint(shape*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#else
            const double ad = hypot(er[i]-zr[j],ei[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else if (ad != 0.0)
              P = (log(ad)+0.5*expint(shape*ad*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#endif
            qr[i] -= mr[j]*P;
          }
        break;
      }
    }
    else { // zlog()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++)
            if (er[i] != zr[j] || ei[i] != zi[j]) {
#ifdef INLINECOMPLEX_DIRECT
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 0.5*log(re*re+im*im);
#else
              const double P = log(hypot(er[i]-zr[j],ei[i]-zi[j]));
#endif
              qr[i] -= mr[j]*P;
              qi[i] -= mi[j]*P;
            }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*(ad/cutoff)+shape;
#else
            const double ad = hypot(er[i]-zr[j],ei[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*(ad/cutoff)*(ad/cutoff)+shape;
#endif
            qr[i] -= mr[j]*P;
            qi[i] -= mi[j]*P;
          }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else
              P = 0.5*log(shape+ad)*scale;
#else
            const double ad = hypot(er[i]-zr[j],ei[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else
              P = 0.5*log(shape+ad*ad)*scale;
#endif
            qr[i] -= mr[j]*P;
            qi[i] -= mi[j]*P;
          }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++)
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 0.5*log(ad);
            else if (ad != 0.0)
              P = 0.5*(log(ad)+expint(shape*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#else
            const double ad = hypot(er[i]-zr[j],ei[i]-zi[j]);
            double P;
            if (ad >= cutoff)
              P = log(ad);
            else if (ad != 0.0)
              P = (log(ad)+0.5*expint(shape*ad*ad))*scale;
            else
              P = -0.5*(log(shape)+gamma_)*scale;
#endif
            qr[i] -= mr[j]*P;
            qi[i] -= mi[j]*P;
          }
        break;
      }
    }
  }
  else { // harmonic potential
    if (mi == NULL) { // rinv()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2 && j&1) {
            if (er[i] != zr[j] || ei[i] != zi[j]) {
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 1.0/(re*re+im*im);
              qr[i] -= mr[j]*P*re;
              qi[i] += mr[j]*P*im;
            }
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d P=_mm_load_pd(&mr[j]);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            P=_mm_div_pd(P,mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad==0) {
              double Pd=0;
              P=_mm_loadl_pd(P, &Pd);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad==0) {
              double Pd=0;
              P=_mm_loadh_pd(P, &Pd);
            }
            pri=_mm_add_pd(pri, _mm_mul_pd(P, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(P, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            if (er[i] != zr[j] || ei[i] != zi[j]) {
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 1.0/(re*re+im*im);
              qr[i] -= mr[j]*P*re;
              qi[i] += mr[j]*P*im;
            }
          }
#else /*SSE2DIRECT*/
          for (int j = begin2; j < end2; j++)
            if (er[i] != zr[j] || ei[i] != zi[j]) {
#ifdef INLINECOMPLEX_DIRECT
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 1.0/(re*re+im*im);
              qr[i] -= mr[j]*P*re;
              qi[i] += mr[j]*P*im;
#else
              const dcmplx P = 1.0/(er[i]-zr[j]+I*(ei[i]-zi[j]));
              qr[i] -= mr[j]*creal(P);
              qi[i] -= mr[j]*cimag(P);
#endif
            }
#endif /*SSE2DIRECT*/
        }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2&&j&1) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            qr[i] -= mr[j]*P*re;
            qi[i] += mr[j]*P*im;
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadl_pd(P, &shape);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadh_pd(P, &shape);
            }
            mmr=_mm_load_pd(&mr[j]);
            mmr=_mm_mul_pd(mmr,P);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            qr[i] -= mr[j]*P*re;
            qi[i] += mr[j]*P*im;
          }
#else /*SSE2DIRECT*/
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            qr[i] -= mr[j]*P*re;
            qi[i] += mr[j]*P*im;
#else
            const dcmplx d = er[i]-zr[j]+I*(ei[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)*shape;
            qr[i] -= mr[j]*creal(P);
            qi[i] -= mr[j]*cimag(P);
#endif
          }
#endif
        }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2&&j&1) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            qr[i] -= mr[j]*P*re;
            qi[i] += mr[j]*P*im;
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadl_pd(P, &Pd);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadh_pd(P, &Pd);
            }
            mmr=_mm_load_pd(&mr[j]);
            mmr=_mm_mul_pd(mmr,P);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            qr[i] -= mr[j]*P*re;
            qi[i] += mr[j]*P*im;
          }
#else /*SSE2DIRECT*/
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            qr[i] -= mr[j]*P*re;
            qi[i] += mr[j]*P*im;
#else
            const dcmplx d = er[i]-zr[j]+I*(ei[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)/(shape+ad*ad)*scale;
            qr[i] -= mr[j]*creal(P);
            qi[i] -= mr[j]*cimag(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2&&j&1) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              qr[i] -= mr[j]*P*re;
              qi[i] += mr[j]*P*im;
            }
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d P=_mm_load_pd(&mr[j]);
            __m128d mad=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            P=_mm_div_pd(P,mad);
#ifdef OSEENSSEALTERNATIVE
            double VC16ALIGN ad[2],Pd[2] GCC16ALIGN;
            _mm_store_pd(ad,mad);
            if(ad[0]<cutoff||ad[1]<cutoff) {
              _mm_store_pd(Pd,P);
              if(ad[0]<cutoff) {
                if(ad[0]!=0.0)
                  Pd[0]*=-expm1(-shape*ad[0])*scale;
                else
                  Pd[0]=0;
              }
              if(ad[1]<cutoff) {
                if(ad[1]!=0.0)
                  Pd[1]*=-expm1(-shape*ad[1])*scale;
                else
                  Pd[1]=0;
              }
              P=_mm_load_pd(Pd);
            }
#else
            double ad;
            _mm_storel_pd(&ad, mad);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storel_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadl_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadl_pd(P, &Pd);
              }
            }
            _mm_storeh_pd(&ad, mad);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storeh_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadh_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadh_pd(P, &Pd);
              }
            }
#endif
            pri=_mm_add_pd(pri, _mm_mul_pd(P, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(P, mim));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              qr[i] -= mr[j]*P*re;
              qi[i] += mr[j]*P*im;
            }
          }
#else /*SSE2DIRECT*/
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              qr[i] -= mr[j]*P*re;
              qi[i] += mr[j]*P*im;
            }
#else
            const dcmplx d = er[i]-zr[j]+I*(ei[i]-zi[j]);
            const double ad = cabs(d);
            if (ad != 0.0) {
              dcmplx P;
              if (ad >= cutoff)
                P = 1.0/d;
              else
                P = -expm1(-shape*ad*ad)/d*scale;
              qr[i] -= mr[j]*creal(P);
              qi[i] -= mr[j]*cimag(P);
            }
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;
      }
    }
    else { // zinv()
      switch (smooth) {
      case DIRAC:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2 && j&1) {
            if (er[i] != zr[j] || ei[i] != zi[j]) {
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 1.0/(re*re+im*im);
              qr[i] -= (mr[j]*re+mi[j]*im)*P;
              qi[i] += (mr[j]*im-mi[j]*re)*P;
            }
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad==0) {
              double Pd=0;
              P=_mm_loadl_pd(P, &Pd);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad==0) {
              double Pd=0;
              P=_mm_loadh_pd(P, &Pd);
            }
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            if (er[i] != zr[j] || ei[i] != zi[j]) {
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 1.0/(re*re+im*im);
              qr[i] -= (mr[j]*re+mi[j]*im)*P;
              qi[i] += (mr[j]*im-mi[j]*re)*P;
            }
          }
#else
          for (int j = begin2; j < end2; j++)
            if (er[i] != zr[j] || ei[i] != zi[j]) {
#ifdef INLINECOMPLEX_DIRECT
              const double re = er[i]-zr[j],im = ei[i]-zi[j];
              const double P = 1.0/(re*re+im*im);
              qr[i] -= (mr[j]*re+mi[j]*im)*P;
              qi[i] += (mr[j]*im-mi[j]*re)*P;
#else
              const dcmplx P = 1.0/(er[i]-zr[j]+I*(ei[i]-zi[j]));
              qr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
              qi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
#endif
            }
#endif /*SSE2DIRECT*/
        }
        break;

      case RANKINE:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2&&j&1) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            qr[i] -= (mr[j]*re+mi[j]*im)*P;
            qi[i] += (mr[j]*im-mi[j]*re)*P;
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadl_pd(P, &shape);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              P=_mm_loadh_pd(P, &shape);
            }
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            qr[i] -= (mr[j]*re+mi[j]*im)*P;
            qi[i] += (mr[j]*im-mi[j]*re)*P;
          }
#else
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = shape;
            qr[i] -= (mr[j]*re+mi[j]*im)*P;
            qi[i] += (mr[j]*im-mi[j]*re)*P;
#else
            const dcmplx d = er[i]-zr[j]+I*(ei[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)*shape;
            qr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
            qi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case SCULLY:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2&&j&1) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            qr[i] -= (mr[j]*re+mi[j]*im)*P;
            qi[i] += (mr[j]*im-mi[j]*re)*P;
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadl_pd(P, &Pd);
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              double Pd = 1.0/(shape+ad)*scale;
              P=_mm_loadh_pd(P, &Pd);
            }
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            qr[i] -= (mr[j]*re+mi[j]*im)*P;
            qi[i] += (mr[j]*im-mi[j]*re)*P;
          }
#else
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad >= cutoff)
              P = 1.0/ad;
            else
              P = 1.0/(shape+ad)*scale;
            qr[i] -= (mr[j]*re+mi[j]*im)*P;
            qi[i] += (mr[j]*im-mi[j]*re)*P;
#else
            const dcmplx d = er[i]-zr[j]+I*(ei[i]-zi[j]);
            const double ad = cabs(d);
            dcmplx P;
            if (ad >= cutoff)
              P = 1.0/d;
            else
              P = conj(d)/(shape+ad*ad)*scale;
            qr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
            qi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;

      case OSEEN:
        for (int i = begin1; i < end1; i++) {
#ifdef SSE2DIRECT
          int j = begin2;
          if (j < end2&&j&1) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              qr[i] -= (mr[j]*re+mi[j]*im)*P;
              qi[i] += (mr[j]*im-mi[j]*re)*P;
            }
            j++;
          }
          __m128d mzr=_mm_set1_pd(er[i]);
          __m128d mzi=_mm_set1_pd(ei[i]);
          __m128d pri=_mm_set1_pd(0);
          __m128d pii=_mm_set1_pd(0);
          for (; j < end2-1; j+=2) {
            __m128d mre=_mm_load_pd(&zr[j]);
            __m128d mim=_mm_load_pd(&zi[j]);
            mre=_mm_sub_pd(mre, mzr);
            mim=_mm_sub_pd(mim, mzi);
            __m128d mmr=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
            __m128d P=_mm_div_pd(_mm_set1_pd(1.0),mmr);
#ifdef OSEENSSEALTERNATIVE
            double VC16ALIGN ad[2],Pd[2] GCC16ALIGN;
            _mm_store_pd(ad,mmr);
            if(ad[0]<cutoff||ad[1]<cutoff) {
              _mm_store_pd(Pd,P);
              if(ad[0]<cutoff) {
                if(ad[0]!=0.0)
                  Pd[0]*=-expm1(-shape*ad[0])*scale;
                else
                  Pd[0]=0;
              }
              if(ad[1]<cutoff) {
                if(ad[1]!=0.0)
                  Pd[1]*=-expm1(-shape*ad[1])*scale;
                else
                  Pd[1]=0;
              }
              P=_mm_load_pd(Pd);
            }
#else
            double ad;
            _mm_storel_pd(&ad, mmr);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storel_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadl_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadl_pd(P, &Pd);
              }
            }
            _mm_storeh_pd(&ad, mmr);
            if(ad<cutoff) {
              if(ad!=0.0) {
                double Pd;
                _mm_storeh_pd(&Pd, P);
                Pd*= -expm1(-shape*ad)*scale;
                P=_mm_loadh_pd(P, &Pd);
              }
              else {
                double Pd=0;
                P=_mm_loadh_pd(P, &Pd);
              }
            }
#endif
            mmr=_mm_load_pd(&mr[j]);
            __m128d mmir=_mm_load_pd(&mi[j]);
            mmr=_mm_mul_pd(mmr, P);
            mmir=_mm_mul_pd(mmir, P);
            pri=_mm_add_pd(pri, _mm_mul_pd(mmr, mre));
            pii=_mm_sub_pd(pii, _mm_mul_pd(mmr, mim));
            pri=_mm_add_pd(pri, _mm_mul_pd(mmir, mim));
            pii=_mm_add_pd(pii, _mm_mul_pd(mmir, mre));
          }
          double VC16ALIGN a2[2] GCC16ALIGN;
          _mm_store_pd(a2, pri);
          qr[i]+=a2[0]+a2[1];
          _mm_store_pd(a2, pii);
          qi[i]+=a2[0]+a2[1];
          if (j < end2) {
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              qr[i] -= (mr[j]*re+mi[j]*im)*P;
              qi[i] += (mr[j]*im-mi[j]*re)*P;
            }
          }
#else
          for (int j = begin2; j < end2; j++) {
#ifdef INLINECOMPLEX_DIRECT
            const double re = er[i]-zr[j],im = ei[i]-zi[j];
            const double ad = re*re+im*im;
            double P;
            if (ad != 0.0) {
              if (ad >= cutoff)
                P = 1.0/ad;
              else
                P = -expm1(-shape*ad)/ad*scale;
              qr[i] -= (mr[j]*re+mi[j]*im)*P;
              qi[i] += (mr[j]*im-mi[j]*re)*P;
            }
#else
            const dcmplx d = er[i]-zr[j]+I*(ei[i]-zi[j]);
            const double ad = cabs(d);
            if (ad != 0.0) {
              dcmplx P;
              if (ad >= cutoff)
                P = 1.0/d;
              else
                P = -expm1(-shape*ad*ad)/d*scale;
              qr[i] -= mr[j]*creal(P)-mi[j]*cimag(P);
              qi[i] -= mr[j]*cimag(P)+mi[j]*creal(P);
            }
#endif
          }
#endif /*SSE2DIRECT*/
        }
        break;
      }
    }
  }
}
/*------------------------------------------------------------------------*/
//debug function to check if points are within box limits
void checkboxes(MPexp *This,const double* zr,const double* zi,const double *er,const double *ei,int Nf,int Nt,int NE)
{
  double *rmaxbox=(double*)mxCalloc(Nt,sizeof(double));
  double *rminbox=(double*)mxCalloc(Nt,sizeof(double));
  double *imaxbox=(double*)mxCalloc(Nt,sizeof(double));
  double *iminbox=(double*)mxCalloc(Nt,sizeof(double));
  for(int k=0;k<Nf;k++) {
    if(This->ixptr[k+1]>This->ixptr[k]) {
      int ind=This->lptr[This->nlevel]+k;
      rmaxbox[ind]=zr[This->ix[This->ixptr[k]]];
      rminbox[ind]=zr[This->ix[This->ixptr[k]]];
      imaxbox[ind]=zi[This->ix[This->ixptr[k]]];
      iminbox[ind]=zi[This->ix[This->ixptr[k]]];
      for(int m=This->ixptr[k]+1;m<This->ixptr[k+1];m++) {
        if(zr[This->ix[m]]>rmaxbox[ind])
          rmaxbox[ind]=zr[This->ix[m]];
        if(zr[This->ix[m]]<rminbox[ind])
          rminbox[ind]=zr[This->ix[m]];
        if(zi[This->ix[m]]>imaxbox[ind])
          imaxbox[ind]=zi[This->ix[m]];
        if(zi[This->ix[m]]<iminbox[ind])
          iminbox[ind]=zi[This->ix[m]];
      }
      if(NE!=0) {
        for(int m=This->jxptr[k]+1;m<This->jxptr[k+1];m++) {
          if(er[This->jx[m]]>rmaxbox[ind])
            rmaxbox[ind]=er[This->jx[m]];
          if(er[This->jx[m]]<rminbox[ind])
            rminbox[ind]=er[This->jx[m]];
          if(ei[This->jx[m]]>imaxbox[ind])
            imaxbox[ind]=ei[This->jx[m]];
          if(ei[This->jx[m]]<iminbox[ind])
            iminbox[ind]=ei[This->jx[m]];
        }
      }
    }
  }
  for(int k=This->nlevel-1;k>=0;k--) {
    for(int m=0;m<This->lptr[k+1]-This->lptr[k];m++) {
      int ind=This->lptr[k]+m;
      int ind2=This->lptr[k+1]+m*4;
      rmaxbox[ind]=rmaxbox[ind2];
      rminbox[ind]=rminbox[ind2];
      imaxbox[ind]=imaxbox[ind2];
      iminbox[ind]=iminbox[ind2];
      for(int n=1;n<4;n++) {
        if(rmaxbox[ind2+n]>rmaxbox[ind])
          rmaxbox[ind]=rmaxbox[ind2+n];
        if(rminbox[ind2+n]<rminbox[ind])
          rminbox[ind]=rminbox[ind2+n];
        if(imaxbox[ind2+n]>imaxbox[ind])
          imaxbox[ind]=imaxbox[ind2+n];
        if(iminbox[ind2+n]<iminbox[ind])
          iminbox[ind]=iminbox[ind2+n];
      }
    }
  }
  for(int k=0;k<=This->nlevel;k++) {
    for(int m=0;m<This->lptr[k+1]-This->lptr[k];m++) {
      int ind=This->lptr[k]+m;
      if(rmaxbox[ind]>creal(This->root[ind].z0)+creal(This->root[ind].d0)*1.000001)
        mexPrintf("Level %d box %d real value %e larger than box %e\n",k,m,rmaxbox[ind],creal(This->root[ind].z0)+creal(This->root[ind].d0));
      if(rminbox[ind]<creal(This->root[ind].z0)-creal(This->root[ind].d0)*1.000001)
        mexPrintf("Level %d box %d real value %e smaller than box %e\n",k,m,rminbox[ind],creal(This->root[ind].z0)-creal(This->root[ind].d0));
      if(imaxbox[ind]>cimag(This->root[ind].z0)+cimag(This->root[ind].d0)*1.000001)
        mexPrintf("Level %d box %d imag value %e larger than box %e\n",k,m,imaxbox[ind],cimag(This->root[ind].z0)+cimag(This->root[ind].d0));
      if(iminbox[ind]<cimag(This->root[ind].z0)-cimag(This->root[ind].d0)*1.000001)
        mexPrintf("Level %d box %d imag value %e smaller than box %e\n",k,m,iminbox[ind],cimag(This->root[ind].z0)-cimag(This->root[ind].d0));
      mexPrintf("%d: min=%.16e+%.16ei,max=%.16e+%.16ei,z0=%.16e+%.16ei,d0=%.16e+%.16ei\n",ind,rminbox[ind],iminbox[ind],rmaxbox[ind],imaxbox[ind],creal(This->root[ind].z0),cimag(This->root[ind].z0),creal(This->root[ind].d0),cimag(This->root[ind].d0));
    }
  }
  mxFree(rmaxbox);
  mxFree(rminbox);
  mxFree(imaxbox);
  mxFree(iminbox);
}
