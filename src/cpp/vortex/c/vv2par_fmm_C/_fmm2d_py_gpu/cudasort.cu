/* cudasort.cu */

/* S. Engblom and A. Goude 2011-10-21 */


#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#ifndef SWAP
#define SWAP(x,y,tmp) (tmp)=(x);(x)=(y);(y)=(tmp);
#endif

#ifndef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#include "mex.h"
#include "matrix.h"
#endif
#include "cudaeval.h"
#include "cudasort.h"

__host__ __device__ int imin(int x, int y);
__host__ __device__ int imax(int x, int y);

// include all CUDA kernels
#include "cudasortkernels.h"

#if defined(CUDASUPPORT) && defined(CUDASORT)
// prototypes
int checkpartitioning(const SORT_REAL* originalpositions,
                      const SORT_REAL *originalpositions2,
                      SORT_REAL *cudapositions,
                      SORT_REAL *cudapositions2,
                      SORT_REAL *splitpoints,
                      int *cudaindices,int *oldcudaindices,
                      int count,int *cudallimits,int *cudarlimits,
                      int splitcount,
                      int *cudaysplit,
                      int printsuccess,int printfailure,
                      SORT_DCMPLX* z,SORT_DCMPLX *d);
void multiblockpartition(SORT_REAL *positions,SORT_REAL *positions2,
                         int *indices,int *newindices,
                         SORT_REAL *splitpoints,
                         int *llimits,int *rlimits,int *newllimits,
                         int *newrlimits,int *tmpllimits,int *tmprlimits,
                         int splitcount,int count,
                         SORT_REAL *newpositions,SORT_REAL *newpositions2,
                         int *ysplit,int *lrcount,
                         int *threesplitvector,int *outputvector,int N
                         SORTLIMITSTRING DEBUGVECTORSTRING);
#ifdef SORTLIMIT
void calculatesortlimits(float *distlimit,float *disttarget,float *sortlimit,float* input,float currentlevel,float maxlevel);
#endif
double calculateetalimits(double *eta,double currentlevel,double maxlevel);
/*------------------------------------------------------------------------*/
//this function copies data for CUDA_perform_partitioning
//put in its own function for timing issues only
void CUDA_copy_vectors(cudavariables *GPUvars,
                               int N,int NE,
                               const double* zr,const double *zi,
                               const double *er,const double *ei)
{
  if(NE) {
    cudasafeMalloc((void**)&GPUvars->er,NE*sizeof(double));
    cudasafeMalloc((void**)&GPUvars->ei,NE*sizeof(double));
    cudasafeMalloc((void**)&GPUvars->jx,NE*sizeof(int));
    cudasafe( cudaMemcpy(GPUvars->er, er, NE*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy ertmp" );
    cudasafe( cudaMemcpy(GPUvars->ei, ei, NE*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy eitmp" );
  }
  else {
    GPUvars->er=NULL;
    GPUvars->ei=NULL;
  }
  cudasafeMalloc((void**)&GPUvars->zr,N*sizeof(double));
  cudasafeMalloc((void**)&GPUvars->zi,N*sizeof(double));
  cudasafeMalloc((void**)&GPUvars->ix,N*sizeof(int));
  cudasafe( cudaMemcpy(GPUvars->zr, zr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zrtmp" );
  cudasafe( cudaMemcpy(GPUvars->zi, zi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zitmp" );
}
/*------------------------------------------------------------------------*/
// this is the main function that performes partitioning. Will
// allocate GPUvariables, that will be left on the GPU
// GPUvars-variables for GPU N,NE number of potential/evaluation
// points nlevels, the number of times to split the boxes into four
// smaller ones zr,zr,er,ei coordinates of potential and evaluation
// points eta, for special split with eta-criterion
void CUDA_perform_partitioning(cudavariables *GPUvars,
                               int N,int NE,int nlevels VALIDATEPARTITIONINGSTRING1)
{
  SORT_REAL *splitpoints;
  int* rlimits;
  int* newllimits;
  int* newrlimits;
  int* tmpllimits;
  int* tmprlimits;
  int* newindices;
  int* ysplit;
  int* lrcount;
  int* outputvector;
  int* splitside;
  SORT_REAL* xpositionstmp;
  SORT_REAL* ypositionstmp;
  SORT_REAL* zdmaxpositions;
  SORT_REAL* zdminpositions;
  SORT_REAL* zdmaxpositions2;
  SORT_REAL* zdminpositions2;
#ifdef SORTLIMIT
  SORT_REAL* leftlimitvalues;
  SORT_REAL* rightlimitvalues;
  float distlimit;
  float disttarget;
  float sortlimit;
#endif
  double eta;
  int* zdllimits;
  int* zdrlimits;
  int* threesplit;
  SORT_DCMPLX *ztmp;
  SORT_DCMPLX *dtmp;
  SORT_DCMPLX *z,*d;
  SORT_REAL *dabs;

  int* rlimitsNE;
  int* newllimitsNE;
  int* newrlimitsNE;
  int* lrcountNE;
  int* newindicesNE;
  SORT_REAL* xpositionstmpNE;
  SORT_REAL* ypositionstmpNE;

  int ptrinitvector[2];
  ptrinitvector[0]=0;
  int threadcount=ZDMAXTHREADS;
#ifdef FLOATSORT
  float *fer=NULL,*fei=NULL,*fzr=NULL,*fzi=NULL;
#else
  double *fer,*fei,*fzr,*fzi;
#endif
  checkcudaerror("Partitioning start\n");
  //init
#ifdef CUDADEBUGVECTOR
  cudasafe(cudaMallocDebug((void**)&GPUvars->debugvector,imax((NE<N?N:NE),2000)*sizeof(double)),"cudaMalloc");
  cudasafe(cudaMemset(GPUvars->debugvector,0,imax((NE<N?N:NE),2000)*sizeof(double)),"cudaMemset");
#endif

  //allocate positions and move them to gpu (moved to function above)
  #ifdef FLOATSORT
  if(NE) {
    cudasafeMalloc((void**)&fer,NE*sizeof(float));
    cudasafeMalloc((void**)&fei,NE*sizeof(float));
    int blockcountcf=imin((NE+4*CONVERTTOFLOATMAXTHREADS-1)/(4*CONVERTTOFLOATMAXTHREADS),CONVERTTOFLOATMAXBLOCKS);
    checkcudaerror("before converttofloatNE\n");
    converttofloat<<<blockcountcf,CONVERTTOFLOATMAXTHREADS>>>(fer,fei,GPUvars->er,GPUvars->ei,NE);
    CHECKCUDAERROR
    checkcudaerror("converttofloatNE\n");
  }
  cudasafeMalloc((void**)&fzr,N*sizeof(float));
  cudasafeMalloc((void**)&fzi,N*sizeof(float));
  int blockcountcf=imin((N+4*CONVERTTOFLOATMAXTHREADS-1)/(4*CONVERTTOFLOATMAXTHREADS),CONVERTTOFLOATMAXBLOCKS);
  converttofloat<<<blockcountcf,CONVERTTOFLOATMAXTHREADS>>>(fzr,fzi,GPUvars->zr,GPUvars->zi,N);
  CHECKCUDAERROR
  checkcudaerror("converttofloat\n");
#else
  fer=GPUvars->er;
  fei=GPUvars->ei;
  fzr=GPUvars->zr;
  fzi=GPUvars->zi;
#endif
  //calculate the number of boxes etc.
  int Nf=1,Nt=1;
  for(int i=1;i<=nlevels;i++) {
    Nf<<=2;
    Nt+=Nf;
  }
  if(NE)
    cudasafeMalloc((void**)&GPUvars->jxptr, (Nf+1)*sizeof(int));
  cudasafeMalloc((void**)&GPUvars->ixptr,(Nf+1)*sizeof(int));
  cudasafeMalloc((void**)&GPUvars->z0,Nt*sizeof(SORT_DCMPLX));
  cudasafeMalloc((void**)&GPUvars->d0,Nt*sizeof(SORT_DCMPLX));
  cudasafeMalloc((void **)&GPUvars->dabs,Nt*sizeof(SORT_REAL));
  z=(SORT_DCMPLX*)GPUvars->z0;
  d=(SORT_DCMPLX*)GPUvars->d0;
  dabs=GPUvars->dabs;

#ifdef CUDATIMESORT
  cudaEvent_t start;
  cudaEvent_t stop;
  cudaEvent_t start2;
  cudaEvent_t stop2;
  cudasafe(cudaEventCreate(&start),"cudaEventCreate sort1start");
  cudasafe(cudaEventCreate(&stop),"cudaEventCreate sort1stop");
  cudasafe(cudaEventCreate(&start2),"cudaEventCreate sort1start");
  cudasafe(cudaEventCreate(&stop2),"cudaEventCreate sort1stop");
  cudasafe(cudaEventRecord(start,0),"cudaEventRecord sort1start");
  float elapsedtime,elapsedtime2;
#endif

  //allocate temporary variables. These should be cleaned up afterwards
  cudasafe(cudaMallocDebug((void**)&xpositionstmp, N*sizeof(SORT_REAL)), "cudaMalloc xpositionstmp");
  cudasafe(cudaMallocDebug((void**)&ypositionstmp, N*sizeof(SORT_REAL)), "cudaMalloc ypositionstmp");
  cudasafe(cudaMallocDebug((void**)&newindices, N*sizeof(int)), "cudaMalloc newindices");
  cudasafe(cudaMallocDebug((void**)&lrcount, Nf*sizeof(int)), "cudaMalloc lrcount");
  cudasafe(cudaMallocDebug((void**)&rlimits, Nf*sizeof(int)), "cudaMalloc rlimits");
  cudasafe(cudaMallocDebug((void**)&newllimits, (Nf/2+1)*sizeof(int)), "cudaMalloc newllimits");
  cudasafe(cudaMallocDebug((void**)&newrlimits, Nf/2*sizeof(int)), "cudaMalloc newrlimits");
  cudasafe(cudaMallocDebug((void**)&tmpllimits, Nf/2*sizeof(int)), "cudaMalloc tmpllimits");
  cudasafe(cudaMallocDebug((void**)&tmprlimits, Nf/2*sizeof(int)), "cudaMalloc tmprlimits");
  cudasafe(cudaMallocDebug((void**)&threesplit, Nf/2*sizeof(int)), "cudaMalloc threesplit");
#ifdef SORTLIMIT
  cudasafe(cudaMallocDebug((void**)&leftlimitvalues, Nf/2*sizeof(double)), "cudaMalloc leftlimitvalues");
  cudasafe(cudaMallocDebug((void**)&rightlimitvalues, Nf/2*sizeof(double)), "cudaMalloc rightlimitvalues");
#endif
  cudasafe(cudaMallocDebug((void**)&outputvector, 2*sizeof(int)), "cudaMalloc outputvector");
  cudasafe(cudaMallocDebug((void**)&ysplit, Nf/2*sizeof(int)), "cudaMalloc ysplit");
  cudasafe(cudaMallocDebug((void**)&splitside, Nf/2*sizeof(int)), "cudaMalloc splitside");
  cudasafe(cudaMallocDebug((void**)&splitpoints, Nf/2*sizeof(SORT_REAL)), "cudaMalloc cxpositions");
  cudasafe(cudaMallocDebug((void**)&zdllimits, MAXBLOCKSZDMAXMULTI*sizeof(int)), "cudaMalloc zdllimits");
  cudasafe(cudaMallocDebug((void**)&zdrlimits, MAXBLOCKSZDMAXMULTI*sizeof(int)), "cudaMalloc zdrlimits");
  cudasafe(cudaMallocDebug((void**)&zdmaxpositions, MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)), "cudaMalloc zdmaxpositions");
  cudasafe(cudaMallocDebug((void**)&zdminpositions, MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)), "cudaMalloc zdminpositions");
  cudasafe(cudaMallocDebug((void**)&zdmaxpositions2, MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)), "cudaMalloc zdmaxpositions2");
  cudasafe(cudaMallocDebug((void**)&zdminpositions2, MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)), "cudaMalloc zdminpositions2");
  cudasafe(cudaMallocDebug((void**)&ztmp, Nf/2*sizeof(SORT_DCMPLX)), "cudaMalloc cxpositions");
  cudasafe(cudaMallocDebug((void**)&dtmp, Nf/2*sizeof(SORT_DCMPLX)), "cudaMalloc cxpositions");

  if(NE) {
    cudasafe(cudaMallocDebug((void**)&xpositionstmpNE, NE*sizeof(SORT_REAL)), "cudaMalloc xpositionstmpNE");
    cudasafe(cudaMallocDebug((void**)&ypositionstmpNE, NE*sizeof(SORT_REAL)), "cudaMalloc ypositionstmpNE");
    cudasafe(cudaMallocDebug((void**)&newindicesNE, NE*sizeof(int)), "cudaMalloc newindicesNE");
    cudasafe(cudaMallocDebug((void**)&newllimitsNE, (Nf/2+1)*sizeof(int)), "cudaMalloc newllimitsNE");
    cudasafe(cudaMallocDebug((void**)&newrlimitsNE, Nf/2*sizeof(int)), "cudaMalloc newrlimitsNE");
    cudasafe(cudaMallocDebug((void**)&rlimitsNE, Nf*sizeof(int)), "cudaMalloc rlimitsNE");
    cudasafe(cudaMallocDebug((void**)&lrcountNE, Nf*sizeof(int)), "cudaMalloc lrcountNE");
    ptrinitvector[1]=NE; //make sure the last element is there in case of only 1 level
    cudasafe(cudaMemcpy(GPUvars->jxptr,ptrinitvector, 2*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy GPUvars->jxptr");
    cudasafe( cudaMemcpy(rlimitsNE, &NE, sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy rlimitsNE" );
  }
  checkcudaerror("allocation\n");
  //initiation
  ptrinitvector[1]=N; //this is necessary in case only one level is used. Otherwise, it has no point
  cudasafe(cudaMemcpy(GPUvars->ixptr,ptrinitvector, 2*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy GPUvars->ixptr");
  cudasafe( cudaMemcpy(rlimits, &N, sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy rlimits" );

  //useful number of blocks for calls. Could probably be optimized
  int blockcount=imin((N+4*threadsperblock-1)/(4*threadsperblock),maxblockcount);
  int blockcountNE=imin((NE+4*threadsperblock-1)/(4*threadsperblock),maxblockcount);
  int zdblockcount=imin((imax(N,NE)+4*ZDMAXTHREADS-1)/(4*ZDMAXTHREADS),MAXBLOCKSZDMAXMULTI);

#ifdef CUDATIMESORT
  cudasafe(cudaEventRecord(start2,0),"cudaEventRecord stop2");
#endif
  //initiate all indices to [1,2,3,....]
  checkcudaerror("before initiateindices\n");
  initiateindices<<<blockcount,threadsperblock>>>(GPUvars->ix,N);
  CHECKCUDAERROR
  checkcudaerror("initiateindices\n");
  if(NE) {
    initiateindices<<<blockcountNE,threadsperblock>>>(GPUvars->jx,NE);
    CHECKCUDAERROR
    checkcudaerror("initiateindicesNE\n");
  }

#ifdef CUDATIMESORT
  cudasafe(cudaEventRecord(stop2),"cudaEventRecord fulltimestop");
  cudasafe(cudaEventSynchronize(stop2),"cudaEventSynchronize fulltimestop");
  cudasafe(cudaEventElapsedTime(&elapsedtime,start2,stop2),"cudaEventElapsedTime");
  mexPrintf("Initiateindices, init time: %f\n",elapsedtime/1000);
  cudasafe(cudaEventRecord(start2,0),"cudaEventRecord sort1start");
#endif

  //determine base block size
  cudasafe(cudaMemset(outputvector,0,sizeof(int)),"Memset outputvector");
  findzdmulti<<<zdblockcount,ZDMAXTHREADS>>>(GPUvars->ixptr,rlimits,fzr,fzi,GPUvars->jxptr,rlimitsNE,fer,fei,zdllimits,zdrlimits,zdmaxpositions,zdmaxpositions2,zdminpositions,zdminpositions2,1,1,outputvector DEBUGVECTORSTRING2);
  CHECKCUDAERROR
  checkcudaerror("findzdmulti\n");
  int hasnan;
  cudasafe(cudaMemcpy(&hasnan,outputvector, sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy coutputvector");
  if(hasnan) {
    cudaThreadExit();
    resetalloccount();
    mexErrMsgTxt("NaN detected in input vectors, aborting");
  }
#ifdef CUDATIMESORT
  cudasafe(cudaEventRecord(stop2),"cudaEventRecord fulltimestop");
  cudasafe(cudaEventSynchronize(stop2),"cudaEventSynchronize fulltimestop");
  cudasafe(cudaEventElapsedTime(&elapsedtime,start2,stop2),"cudaEventElapsedTime");
  mexPrintf("findzdmulti, init time: %f\n",elapsedtime/1000);
  cudasafe(cudaEventRecord(start2,0),"cudaEventRecord sort1start");
#endif
  findzdmultistep2<<<1,ZDMAXTHREADS>>>(zdllimits,zdrlimits,zdmaxpositions,zdmaxpositions2,zdminpositions,zdminpositions2,z,d,dabs,1);
  CHECKCUDAERROR
  checkcudaerror("findzdmultistep2\n");
#ifdef CUDATIMESORT
  cudasafe(cudaEventRecord(stop2),"cudaEventRecord fulltimestop");
  cudasafe(cudaEventSynchronize(stop2),"cudaEventSynchronize fulltimestop");
  cudasafe(cudaEventElapsedTime(&elapsedtime,start2,stop2),"cudaEventElapsedTime");
  mexPrintf("findzdmultistep2, init time: %f\n",elapsedtime/1000);
#endif

#ifdef CUDATIMESORT
  cudasafe(cudaEventRecord(stop),"cudaEventRecord fulltimestop");
  cudasafe(cudaEventSynchronize(stop),"cudaEventSynchronize fulltimestop");
  cudasafe(cudaEventElapsedTime(&elapsedtime,start,stop),"cudaEventElapsedTime");
  mexPrintf("Partition, init time: %f\n",elapsedtime/1000);
  cudasafe(cudaEventRecord(start,0),"cudaEventRecord sort1start");
#endif

          //the main partitioning loop
  for(int i=0,Nb=1;i<nlevels;i++) {
//     mexPrintf("partitioning step %d of %d\n",i,nlevels);
#ifdef SORTLIMIT
    calculatesortlimits(&distlimit,&disttarget,&sortlimit,GPUvars->sortlimits,i,nlevels-0.5);
#endif
    eta=calculateetalimits(GPUvars->eta,i,nlevels-0.5);
    while(threadcount>imax(N,NE)/Nb&&threadcount>32)
      threadcount>>=1;
    if(threadcount<32) threadcount=32; //should not happen
    //setup
    int singleblockwork=(Nb+SINGLETHREADTHREADCOUNT-1)/SINGLETHREADTHREADCOUNT;
    checkcudaerror("before setuppartition\n");
    setuppartition<<<imin(singleblockwork,SINGLETHREADMAXTHREADS),SINGLETHREADTHREADCOUNT>>>(z,d,splitpoints,ysplit,1,Nb SORTLIMITCALLINGSTRING2);
    CHECKCUDAERROR
    cudasafe(cudaMemset(lrcount,0, Nb*2*sizeof(int)), "cudaMemset lrcount");
    cudasafe(cudaMemset(threesplit,0, Nb*sizeof(int)), "cudaMemset threesplit");
#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(start2,0),"cudaEventRecord sort1start");
#endif

    //make first partitioning
    checkcudaerror("before multiblockpartition\n");
    multiblockpartition(fzr,fzi,GPUvars->ix,newindices,splitpoints,GPUvars->ixptr,rlimits,newllimits,newrlimits,tmpllimits,tmprlimits,Nb,N,xpositionstmp,ypositionstmp,ysplit,lrcount,threesplit,outputvector,N SORTLIMITCALLINGSTRING DEBUGVECTORSTRING2/*,zr,zi*/);
    if(eta<1) { //if eta<1, make one additional split
      setupetasplit<<<imin(singleblockwork,SINGLETHREADMAXTHREADS),SINGLETHREADTHREADCOUNT>>>(z,splitpoints,ysplit,eta,Nb,splitside,newllimits,newrlimits,tmpllimits,tmprlimits);
      CHECKCUDAERROR
      cudasafe(cudaMemset(lrcount,0, Nb*2*sizeof(int)), "cudaMemset lrcount");

      //make partitioning, and copy data back to original array
      if(blockcount>Nb) { //single or multiblock mode?
        partitionsplit<<<blockcount,threadsperblock>>>(xpositionstmp,ypositionstmp,newindices,GPUvars->ix,splitpoints,tmpllimits,tmprlimits,lrcount,Nb,fzr,fzi,ysplit,NULL DEBUGVECTORSTRING2);
        CHECKCUDAERROR
        splitetacopymultithread<<<blockcount,threadsperblock>>>(tmpllimits,tmprlimits,fzr,fzi,GPUvars->ix,xpositionstmp,ypositionstmp,newindices,Nb);
        CHECKCUDAERROR
      }
      else {
        partitionsplitsinglethread<<<blockcount,threadsperblock>>>(xpositionstmp,ypositionstmp,newindices,GPUvars->ix,splitpoints,tmpllimits,tmprlimits,lrcount,Nb,fzr,fzi,ysplit DEBUGVECTORSTRING2);
        CHECKCUDAERROR
        splitetacopysinglethread<<<blockcount,threadsperblock>>>(tmpllimits,tmprlimits,fzr,fzi,GPUvars->ix,xpositionstmp,ypositionstmp,newindices,Nb);
        CHECKCUDAERROR
      }
      correctetalimits<<<imin(singleblockwork,SINGLETHREADMAXTHREADS),SINGLETHREADTHREADCOUNT>>>(newllimits,newrlimits,splitside,lrcount,Nb);
      CHECKCUDAERROR
    }
#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(stop2),"cudaEventRecord stop2");
#endif
    if(NE) { //if evaluation points, split them
      if(blockcountNE>Nb) { //single or multi block mode
        cudasafe(cudaMemset(lrcountNE, 0, Nb*2*sizeof(int)), "cudaMemset lrcount");
        partitionsplit<<<blockcountNE, threadsperblock>>>(fer, fei, GPUvars->jx, newindicesNE, splitpoints, GPUvars->jxptr, rlimitsNE, lrcountNE, Nb, xpositionstmpNE, ypositionstmpNE, ysplit, NULL DEBUGVECTORSTRING2);
        CHECKCUDAERROR
      }
      else {
        partitionsplitsinglethread<<<blockcountNE, threadsperblock>>>(fer, fei, GPUvars->jx, newindicesNE, splitpoints, GPUvars->jxptr, rlimitsNE, lrcountNE, Nb, xpositionstmpNE, ypositionstmpNE, ysplit DEBUGVECTORSTRING2);
        CHECKCUDAERROR
      }
      setNElimits<<<imin(singleblockwork, SINGLETHREADMAXTHREADS), SINGLETHREADTHREADCOUNT>>>(GPUvars->jxptr, rlimitsNE, newllimitsNE, newrlimitsNE, lrcountNE, Nb);
      CHECKCUDAERROR
    }
    //in the middle step, only approximate z and d, since these values will no be used anymore (except for the eta split, and as starting values in the next partitioning)
    approximatezd<<<imin(singleblockwork,SINGLETHREADMAXTHREADS),SINGLETHREADTHREADCOUNT>>>(z,d,splitpoints,ztmp,dtmp,ysplit,Nb);
    CHECKCUDAERROR

#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(stop),"cudaEventRecord stop 287");
    cudasafe(cudaEventSynchronize(stop),"cudaEventSynchronize stop");
    cudasafe(cudaEventSynchronize(stop2),"cudaEventSynchronize stop2");
    cudasafe(cudaEventElapsedTime(&elapsedtime,start,stop),"cudaEventElapsedTime");
    cudasafe(cudaEventElapsedTime(&elapsedtime2,start2,stop2),"cudaEventElapsedTime");
    mexPrintf("Partition, loop %d: %f partitiontime: %f\n",i,elapsedtime/1000,elapsedtime2/1000);
#endif

#ifdef VALIDATEPARTITIONING
    checkpartitioning(zr,zi,xpositionstmp,ypositionstmp,splitpoints,newindices,GPUvars->ix,N,newllimits,newrlimits,Nb,ysplit,0,1,ztmp,dtmp);
#endif

    z+=Nb;
    d+=Nb;
    dabs+=Nb;

#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(start,0),"cudaEventRecord sort1start");
#endif

    Nb<<=1;
#ifdef SORTLIMITS
    calculatesortlimits(&distlimit,&disttarget,&sortlimit,GPUvars->sortlimits,i+0.5,nlevels-0.5);
#endif
    eta=calculateetalimits(GPUvars->eta,i+0.5,nlevels-0.5);
    singleblockwork=(Nb+SINGLETHREADTHREADCOUNT-1)/SINGLETHREADTHREADCOUNT;
    cudasafe(cudaMemset(lrcount,0, Nb*2*sizeof(int)), "cudaMemset lrcount");
    cudasafe(cudaMemset(threesplit,0, Nb*sizeof(int)), "cudaMemset threesplit");
    setuppartition<<<imin(singleblockwork,SINGLETHREADMAXTHREADS),SINGLETHREADTHREADCOUNT>>>(ztmp,dtmp,splitpoints,ysplit,0,Nb SORTLIMITCALLINGSTRING2);
    CHECKCUDAERROR
#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(start2,0),"cudaEventRecord sort1start");
#endif

    //second partitioning on this level
    checkcudaerror("before second multiblockpartition\n");
    multiblockpartition(xpositionstmp,ypositionstmp,newindices,GPUvars->ix,splitpoints,newllimits,newrlimits,GPUvars->ixptr,rlimits,tmpllimits,tmprlimits,Nb,N,fzr,fzi,ysplit,lrcount,threesplit,outputvector,N SORTLIMITCALLINGSTRING DEBUGVECTORSTRING2/*,zr,zi*/);


#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(stop2),"cudaEventRecord fulltimestop");
#endif

    if(eta<1) { //eta split again
      setupetasplit<<<imin(singleblockwork, SINGLETHREADMAXTHREADS), SINGLETHREADTHREADCOUNT>>>(ztmp, splitpoints, ysplit, eta, Nb, splitside, GPUvars->ixptr, rlimits, tmpllimits, tmprlimits);
      CHECKCUDAERROR
      cudasafe(cudaMemset(lrcount, 0, Nb*2*sizeof(int)), "cudaMemset lrcount");
      if(blockcount>Nb) {
        partitionsplit<<<blockcount, threadsperblock>>>(fzr, fzi, GPUvars->ix, newindices, splitpoints, tmpllimits, tmprlimits, lrcount, Nb, xpositionstmp, ypositionstmp, ysplit, NULL DEBUGVECTORSTRING2);
        CHECKCUDAERROR
        splitetacopymultithread<<<blockcount, threadsperblock>>>(tmpllimits, tmprlimits, xpositionstmp, ypositionstmp, newindices, fzr, fzi, GPUvars->ix, Nb);
        CHECKCUDAERROR
      }
      else {
        partitionsplitsinglethread<<<blockcount, threadsperblock>>>(fzr, fzi, GPUvars->ix, newindices, splitpoints, tmpllimits, tmprlimits, lrcount, Nb, xpositionstmp, ypositionstmp, ysplit DEBUGVECTORSTRING2);
        CHECKCUDAERROR
        splitetacopysinglethread<<<blockcount, threadsperblock>>>(tmpllimits, tmprlimits, xpositionstmp, ypositionstmp, newindices, fzr, fzi, GPUvars->ix, Nb);
        CHECKCUDAERROR
      }
      correctetalimits<<<imin(singleblockwork, SINGLETHREADMAXTHREADS), SINGLETHREADTHREADCOUNT>>>(GPUvars->ixptr, rlimits, splitside, lrcount, Nb);
      CHECKCUDAERROR
    }
    if(NE) { //evaluation point split
      if(blockcountNE>Nb) {
        cudasafe(cudaMemset(lrcountNE, 0, Nb*2*sizeof(int)), "cudaMemset lrcount");
        partitionsplit<<<blockcountNE, threadsperblock>>>(xpositionstmpNE, ypositionstmpNE, newindicesNE, GPUvars->jx, splitpoints, newllimitsNE, newrlimitsNE, lrcountNE, Nb, fer, fei, ysplit, NULL DEBUGVECTORSTRING2);
        CHECKCUDAERROR
      }
      else {
        partitionsplitsinglethread<<<blockcountNE, threadsperblock>>>(xpositionstmpNE, ypositionstmpNE, newindicesNE, GPUvars->jx, splitpoints, newllimitsNE, newrlimitsNE, lrcountNE, Nb, fer, fei, ysplit DEBUGVECTORSTRING2);
        CHECKCUDAERROR
      }
      setNElimits<<<imin(singleblockwork, SINGLETHREADMAXTHREADS), SINGLETHREADTHREADCOUNT>>>(newllimitsNE, newrlimitsNE, GPUvars->jxptr, rlimitsNE, lrcountNE, Nb);
      CHECKCUDAERROR

    }
#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(stop),"cudaEventRecord stop 350");
    cudasafe(cudaEventSynchronize(stop),"cudaEventSynchronize stop");
    cudasafe(cudaEventSynchronize(stop2),"cudaEventSynchronize stop2");
    cudasafe(cudaEventElapsedTime(&elapsedtime,start,stop),"cudaEventElapsedTime");
    cudasafe(cudaEventElapsedTime(&elapsedtime2,start2,stop2),"cudaEventElapsedTime");
    mexPrintf("Partition, loop %d: %f partitiontime: %f\n",i,elapsedtime/1000,elapsedtime2/1000);
    cudasafe(cudaEventRecord(start,0),"cudaEventRecord sort1start");
#endif
    Nb<<=1;
    singleblockwork=(Nb+SINGLETHREADTHREADCOUNT-1)/SINGLETHREADTHREADCOUNT;

    //now, this will be the actual fmm level. Here, calculate the real values of z and d to use by the theta-criterion etc.
    if(Nb>zdblockcount) { //single or multi block mode
      findzd<<<imin(MAXBLOCKSZDMAXMULTI, Nb), threadcount,threadcount*4*sizeof(SORT_REAL)>>>(GPUvars->ixptr, rlimits, fzr, fzi, GPUvars->jxptr, rlimitsNE, fer, fei, z, d, dabs, Nb);
      CHECKCUDAERROR
    }
    else {
//       mexPrintf("findzdmulti\n");
//       SORT_REAL *hzdmaxpositions=(SORT_REAL*)mxMalloc(MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL));
//       SORT_REAL *hzdminpositions=(SORT_REAL*)mxMalloc(MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL));
//       SORT_REAL *hzdmaxpositions2=(SORT_REAL*)mxMalloc(MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL));
//       SORT_REAL *hzdminpositions2=(SORT_REAL*)mxMalloc(MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL));
//       int *hzdllimits=(int*)mxMalloc(Nb*sizeof(int));
//       int *hzdrlimits=(int*)mxMalloc(Nb*sizeof(int));
//       SORT_DCMPLX *hd0=(SORT_DCMPLX*)mxMalloc(Nb*sizeof(SORT_DCMPLX));
//       SORT_DCMPLX *hz0=(SORT_DCMPLX*)mxMalloc(Nb*sizeof(SORT_DCMPLX));
//       double* hdebugvector=(double*)mxMalloc(MAXBLOCKSZDMAXMULTI*sizeof(double));
//       cudasafe(cudaMemset(GPUvars->debugvector,0,MAXBLOCKSZDMAXMULTI*sizeof(double)),"cudaMeMset");
//       cudasafe(cudaMemset(zdmaxpositions,0,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)),"cudaMemset");
//       cudasafe(cudaMemset(zdminpositions,0,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)),"cudaMemset");
//       cudasafe(cudaMemset(zdmaxpositions2,0,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)),"cudaMemset");
//       cudasafe(cudaMemset(zdminpositions2,0,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL)),"cudaMemset");
      findzdmulti<<<zdblockcount,ZDMAXTHREADS>>>(GPUvars->ixptr,rlimits,fzr,fzi,GPUvars->jxptr,rlimitsNE,fer,fei,zdllimits,zdrlimits,zdmaxpositions,zdmaxpositions2,zdminpositions,zdminpositions2,Nb,0,NULL DEBUGVECTORSTRING2);
      CHECKCUDAERROR
//       cudasafe(cudaMemcpy(hzdmaxpositions,zdmaxpositions,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hzdminpositions,zdminpositions,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hzdmaxpositions2,zdmaxpositions2,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hzdminpositions2,zdminpositions2,MAXBLOCKSZDMAXMULTI*sizeof(SORT_REAL),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hdebugvector,GPUvars->debugvector,MAXBLOCKSZDMAXMULTI*sizeof(double),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hzdllimits,zdllimits,Nb*sizeof(int),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hzdrlimits,zdrlimits,Nb*sizeof(int),cudaMemcpyDeviceToHost),"Memcpy");
      findzdmultistep2<<<Nb,ZDMAXTHREADS>>>(zdllimits,zdrlimits,zdmaxpositions,zdmaxpositions2,zdminpositions,zdminpositions2,z,d,dabs,Nb);
      CHECKCUDAERROR
//       cudasafe(cudaMemcpy(hd0,d,Nb*sizeof(SORT_DCMPLX),cudaMemcpyDeviceToHost),"Memcpy");
//       cudasafe(cudaMemcpy(hz0,z,Nb*sizeof(SORT_DCMPLX),cudaMemcpyDeviceToHost),"Memcpy");
//       if(Nb==1024) {
//       for(int k=0;k<zdblockcount;k++)
//         mexPrintf("limits[%d]: %14e %14e %14e %14e, dv=%e\n",k,hzdminpositions[k],hzdmaxpositions[k],hzdminpositions2[k],hzdmaxpositions2[k],hdebugvector[k]);
//       for(int k=0;k<Nb;k++)
//         mexPrintf("inner size[%d]: [%d %d) z0=%e+%ei d0=%e+%ei\n",k,hzdllimits[k],hzdrlimits[k],creal(hz0[k]),cimag(hz0[k]),creal(hd0[k]),cimag(hd0[k]));
//       }
//       mexPrintf("zdblockcount=%d\n",zdblockcount);
//       mxFree(hzdmaxpositions);
//       mxFree(hzdminpositions);
//       mxFree(hzdmaxpositions2);
//       mxFree(hzdminpositions2);
//       mxFree(hzdllimits);
//       mxFree(hzdrlimits);
//       mxFree(hd0);
//       mxFree(hz0);
//       mxFree(hdebugvector);
    }
#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(stop),"cudaEventRecord stop 370");
    cudasafe(cudaEventSynchronize(stop),"cudaEventSynchronize stop");
    cudasafe(cudaEventElapsedTime(&elapsedtime,start,stop),"cudaEventElapsedTime");
    mexPrintf("find zd, loop %d: %f\n",i,elapsedtime/1000);
#endif

#ifdef VALIDATEPARTITIONING
    checkpartitioning(zr,zi,fzr,fzi,splitpoints,GPUvars->ix,newindices,N,GPUvars->ixptr,rlimits,Nb/2,ysplit,0,1,z,d);
#endif
#ifdef CUDATIMESORT
    cudasafe(cudaEventRecord(start,0),"cudaEventRecord sort1start");
#endif
  }

  //cleanup
#ifdef CUDATIMESORT
  cudasafe(cudaEventDestroy(stop),"cudaEventDestroy");
  cudasafe(cudaEventDestroy(start),"cudaEventDestroy");
  cudasafe(cudaEventDestroy(stop2),"cudaEventDestroy");
  cudasafe(cudaEventDestroy(start2),"cudaEventDestroy");
#endif
  cudaFreeDebug(xpositionstmp);
  cudaFreeDebug(ypositionstmp);
  cudaFreeDebug(newindices);
  cudaFreeDebug(rlimits);
  cudaFreeDebug(newllimits);
  cudaFreeDebug(newrlimits);
  cudaFreeDebug(tmpllimits);
  cudaFreeDebug(tmprlimits);
  cudaFreeDebug(threesplit);
  cudaFreeDebug(outputvector);
  cudaFreeDebug(ysplit);
  cudaFreeDebug(splitside);
  cudaFreeDebug(lrcount);
  cudaFreeDebug(splitpoints);
  cudaFreeDebug(zdllimits);
  cudaFreeDebug(zdrlimits);
  cudaFreeDebug(zdmaxpositions);
  cudaFreeDebug(zdminpositions);
  cudaFreeDebug(zdmaxpositions2);
  cudaFreeDebug(zdminpositions2);
  cudaFreeDebug(ztmp);
  cudaFreeDebug(dtmp);
#ifdef FLOATSORT
  cudaFreeDebug(fzr);
  cudaFreeDebug(fzi);
#endif
#ifdef SORTLIMIT
  cudaFreeDebug(leftlimitvalues);
  cudaFreeDebug(rightlimitvalues);
#endif
  if(NE) {
    cudaFreeDebug(xpositionstmpNE);
    cudaFreeDebug(ypositionstmpNE);
    cudaFreeDebug(newindicesNE);
    cudaFreeDebug(newllimitsNE);
    cudaFreeDebug(newrlimitsNE);
    cudaFreeDebug(rlimitsNE);
    cudaFreeDebug(lrcountNE);
#ifdef FLOATSORT
    cudaFreeDebug(fer);
    cudaFreeDebug(fei);
#endif
  }
}
/*------------------------------------------------------------------------*/
#ifdef CHECKPARTITIONING
#include "cudasortdebugdefs.h"
#endif


/*------------------------------------------------------------------------*/
//multiblockpartition is the equivalent to singleblockpartition, but uses multiple blocks per partitioning. Good in the beginning with few boxes
//positions,positions2 is positions in x and y direction
//indices,newindices is input and output values for the permutation array
//splitpoints is where the split has been performed
//llimits,rlimits indicates where the array to be split are
//newllimits,newrlimits output for llimits/rlimits
//tmpllimits,tmprlimits temporary storage for limits during the algorithm
//splitcount number of partitions to be performed
//count number of elements
//newpositions,newpositions2 output values for positions,positions2
//ysplit split with respect to y or x coordinates
//lrcount for a split, number of elements on each sidet
//threesplitvector if threesplit mode should be used for the box
//outputvector vector on GPU to commuticate results from split
//debugvector debug purposes only
void multiblockpartition(SORT_REAL* positions,SORT_REAL *positions2,int* indices,int* newindices,SORT_REAL *splitpoints,int *llimits,int* rlimits,int *newllimits,int* newrlimits,int* tmpllimits,int* tmprlimits,int splitcount,int count,SORT_REAL* newpositions,SORT_REAL *newpositions2,int* ysplit,int* lrcount,int* threesplitvector,int* outputvector,int N SORTLIMITSTRING DEBUGVECTORSTRING)
{
  int *itmp,i;
  checkcudaerror("start multiblockpartition\n");
#ifdef CHECKPARTITIONING
  static int printcount=0;
  SORT_REAL outsplitpoints[4];
#endif
  SORT_REAL *ctmp;
  int outputvectorlocal[2];
  int threadcount;
  int blockcount=imin((count+4*threadsperblock-1)/(4*threadsperblock),maxblockcount); //determine a reasonable number of blocks to use
  if(splitcount>blockcount) { //if more splits than blocks, use singleblockpartition instead
#ifdef CUDATIMESORT
    cudaEvent_t start;
    cudaEvent_t stop;
    cudasafe(cudaEventCreate(&start), "cudaEventCreate start");
    cudasafe(cudaEventCreate(&stop), "cudaEventCreate stop");
    cudasafe(cudaEventRecord(start, 0), "cudaEventRecord start");
#endif
    threadcount=threadsperblock;
    while(threadcount>N/splitcount&&threadcount>32)
      threadcount>>=1;
    if(threadcount<32) threadcount=32; //should not happen
//     int* hllimits=(int*)mxMalloc((splitcount+1)*sizeof(int));
//     int* hrlimits=(int*)mxMalloc(splitcount*sizeof(int));
//     int* hllimitsnew=(int*)mxMalloc((splitcount+1)*sizeof(int));
//     int* hrlimitsnew=(int*)mxMalloc(splitcount*sizeof(int));
//     cudasafe( cudaMemcpy(hllimits, llimits, splitcount*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(hrlimits, rlimits, splitcount*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
// //     for(int k=0;k<splitcount;k++)
// //       mexPrintf("limits[%d]=[%d %d)\n",k,hllimits[k],hrlimits[k]);
//     for(int k=0;k<splitcount;k++) {
//       if(hrlimits[k]<hllimits[k]||hllimits[k]<0||hrlimits[k]>N)
//         mexPrintf("ERROR: limits[%d]=[%d %d)\n",k,hllimits[k],hrlimits[k]);
//       if(k<splitcount-1&&hrlimits[k]!=hllimits[k+1])
//         mexPrintf("ERROR: not connected %d\n",k);
//     }
//     double *positionsbak,*positions2bak,*splitpointsbak;
//     int* indicesbak;
//     cudasafeMalloc((void **)&positionsbak,N*sizeof(double));
//     cudasafeMalloc((void **)&positions2bak,N*sizeof(double));
//     cudasafeMalloc((void **)&indicesbak,N*sizeof(int));
//     cudasafeMalloc((void **)&splitpointsbak,splitcount*sizeof(int));
//     cudasafe( cudaMemcpy(positionsbak, positions, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(positions2bak, positions2, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(indicesbak, indices, N*sizeof(int), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(splitpointsbak, splitpoints, splitcount*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//     double *hdebugvector=(double*)mxMalloc(imax(2000,splitcount)*sizeof(double));
//     mexPrintf("distlimit=%f disttarget=%f\n",distlimit,disttarget);
//     for(int k=0;k<splitcount;k++) {
//       for(int m=0;m<splitcount;m++) {
//         if(m==k) {
//           hllimitsnew[m]=hllimits[k];
//           hrlimitsnew[m]=hrlimits[k];
//         }
//         else {
//           hllimitsnew[m]=0;
//           hrlimitsnew[m]=0;
//         }
//       }
//       cudasafe( cudaMemset(debugvector,0,imax(2000,splitcount)*sizeof(double)),"dd");
//       cudasafe( cudaMemcpy(llimits, hllimitsnew, splitcount*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(rlimits, hrlimitsnew, splitcount*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy outputvector" );
//       mexPrintf("testing box %d limits = [%d,%d) (%d elements) maxblockcount=%d splitcount=%d\n",k,hllimits[k],hrlimits[k],hrlimits[k]-hllimits[k],maxblockcount,splitcount);
//       mexEvalString("drawnow");
// //       singleblockpartition<<<imin(splitcount, maxblockcount), threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions, positions2, indices, newindices, splitpoints, llimits, rlimits, newllimits, newrlimits, NULL, NULL, splitcount, newpositions, newpositions2, ysplit, 0 SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
//       singleblockpartitiondebug<<<imin(splitcount, maxblockcount), threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions, positions2, indices, newindices, splitpoints, llimits, rlimits, newllimits, newrlimits, NULL, NULL, splitcount, newpositions, newpositions2, ysplit, 0 ,leftlimitvaluesdebug,rightlimitvaluesdebug,distlimit,disttarget DEBUGVECTORSTRING3);
//       cudasafe( cudaMemcpy(positions, positionsbak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(positions2, positions2bak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(indices, indicesbak, N*sizeof(int), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(splitpoints, splitpointsbak, splitcount*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       CHECKCUDAERROR
//       cudasafe( cudaMemcpy(hdebugvector, debugvector, imax(2000,splitcount)*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
//       for(int m=0;m<imax(2000,splitcount);m++)
//         if(hdebugvector[m]!=0.0)
//           mexPrintf("debugvector[%d]=%e\n",m,hdebugvector[m]);
//       mexPrintf("test2\n");
//       mexEvalString("drawnow");
//
// //       if(k==0)
//
//       singleblockpartition<<<imin(splitcount, maxblockcount), threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions, positions2, indices, newindices, splitpoints, llimits, rlimits, newllimits, newrlimits, NULL, NULL, splitcount, newpositions, newpositions2, ysplit, 0 ,leftlimitvalues,rightlimitvalues,distlimit,disttarget DEBUGVECTORSTRING3);
//       cudasafe( cudaMemcpy(positions, positionsbak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(positions2, positions2bak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(indices, indicesbak, N*sizeof(int), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(splitpoints, splitpointsbak, splitcount*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       CHECKCUDAERROR
//       mexPrintf("test3\n");
//       mexEvalString("drawnow");
//       singleblockpartition<<<imin(splitcount, maxblockcount), threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions, positions2, indices, newindices, splitpoints, llimits, rlimits, newllimits, newrlimits, NULL, NULL, splitcount, newpositions, newpositions2, ysplit, 0 ,leftlimitvaluesdebug,rightlimitvaluesdebug,distlimit,disttarget DEBUGVECTORSTRING3);
//       cudasafe( cudaMemcpy(positions, positionsbak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(positions2, positions2bak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(indices, indicesbak, N*sizeof(int), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(splitpoints, splitpointsbak, splitcount*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       CHECKCUDAERROR
//
//     }
//     cudasafe( cudaMemcpy(llimits, hllimits, splitcount*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(rlimits, hrlimits, splitcount*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy outputvector" );
//
//     mxFree(hllimits);
//     mxFree(hrlimits);
//     mxFree(hllimitsnew);
//     mxFree(hrlimitsnew);
//     mxFree(hdebugvector);
//     double* leftlimits=(double*)mxMalloc(splitcount*sizeof(double));
//     double* rightlimits=(double*)mxMalloc(splitcount*sizeof(double));
//     double* hsplitpoints=(double*)mxMalloc(splitcount*sizeof(double));
//     cudasafe( cudaMemcpy(leftlimits, leftlimitvaluesdebug, splitcount*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(rightlimits, rightlimitvaluesdebug, splitcount*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
//     cudasafe( cudaMemcpy(hsplitpoints, splitpoints, splitcount*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
// //     for(int k=0;k<splitcount;k++)
// //       mexPrintf("sortlimits[%d]=[%e %e %e)\n",k,leftlimits[k],hsplitpoints[k],rightlimits[k]);
//     for(int k=0;k<splitcount;k++) {
//       if(hsplitpoints[k]<leftlimits[k]||hsplitpoints[k]>rightlimits[k])
//         mexPrintf("ERROR sortlimits[%d]=[%e %e %e)\n",k,leftlimits[k],hsplitpoints[k],rightlimits[k]);
//     }
//     mxFree(leftlimits);
//     mxFree(rightlimits);
//     mxFree(hsplitpoints);
// //     int* hindices=(int*)mxMalloc(N*sizeof(int));
// //     cudasafe( cudaMemcpy(hindices, indices, N*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
// //     for(int k=0;k<N;k++)
// //       if(hindices[k]<0||hindices[k]>=N)
// //         mexPrintf("Error index[%d]=%d\n",k,hindices[k]);
// //
//     mexPrintf(".");
//     mexEvalString("drawnow");
//     singleblockpartition<<<imin(splitcount, maxblockcount), threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions, positions2, indices, newindices, splitpoints, llimits, rlimits, newllimits, newrlimits, NULL, NULL, splitcount, newpositions, newpositions2, ysplit, 0 ,leftlimitvaluesdebug,rightlimitvaluesdebug,distlimit,disttarget DEBUGVECTORSTRING3);
//     cudasafe( cudaMemcpy(positions, positionsbak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(positions2, positions2bak, N*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(indices, indicesbak, N*sizeof(int), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//       cudasafe( cudaMemcpy(splitpoints, splitpointsbak, splitcount*sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy outputvector" );
//     CHECKCUDAERROR
//     mexPrintf(",");
//     mexEvalString("drawnow");
    singleblockpartition<<<imin(splitcount, maxblockcount), threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions, positions2, indices, newindices, splitpoints, llimits, rlimits, newllimits, newrlimits, NULL, NULL, splitcount, newpositions, newpositions2, ysplit, 0 SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
    CHECKCUDAERROR
//     mexPrintf(":");
//     mexEvalString("drawnow");
//     cudaFreeDebug(positionsbak);
//     cudaFreeDebug(positions2bak);
//     cudaFreeDebug(indicesbak);
//     cudaFreeDebug(splitpointsbak);

#ifdef CUDATIMESORT
    float elapsedtime;
    cudasafe(cudaEventRecord(stop), "cudaEventRecord stop 722");
    cudasafe(cudaEventSynchronize(stop), "cudaEventSynchronize stop");
    cudasafe(cudaEventElapsedTime(&elapsedtime, start, stop), "cudaEventElapsedTime");
    mexPrintf("singleblockpartition (%d blocks): %f\n", imin(splitcount, maxblockcount), elapsedtime/1000);
    cudasafe(cudaEventDestroy(stop), "cudaEventDestroy");
    cudasafe(cudaEventDestroy(start), "cudaEventDestroy");
#endif
    return;
  }
  //now start the real multiblock partition
  blockcount=imax(blockcount,splitcount); //at least as many blocks as splits, but in this case, use singleblockpartition instead
  cudasafe(cudaMemset(lrcount,0, 2*splitcount*sizeof(int)), "cudaMemset lrcount");
  cudasafe(cudaMemset(threesplitvector,0, splitcount*sizeof(int)), "cudaMemset threesplitvector");

  //initial split using the given value in splitpoints (center of box)
  partitionsplit<<<blockcount,threadsperblock>>>(positions,positions2,indices,newindices,splitpoints,llimits,rlimits,lrcount,splitcount,newpositions,newpositions2,ysplit,threesplitvector DEBUGVECTORSTRING3);
  CHECKCUDAERROR

  checkcudaerror("partitionsplit\n");
  SWAP(indices,newindices,itmp);
  SWAP(positions,newpositions,ctmp);
  SWAP(positions2,newpositions2,ctmp);
  cudasafe(cudaMemset(outputvector,0, 2*sizeof(int)), "cudaMemset outputvector");

  //complete the split and set up the next one
#ifdef MEDIAN_OF_32
  preparesplit32<<<splitcount,PREPARESPLIT32THREADCOUNT>>>(positions, positions2,indices,splitpoints,llimits,rlimits,tmpllimits,tmprlimits,llimits,rlimits,lrcount,splitcount,ysplit,threesplitvector,outputvector SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
  CHECKCUDAERROR
#else
  preparesplit<<<(splitcount+PREPARESPLITTHREADCOUNT-1)/PREPARESPLITTHREADCOUNT,PREPARESPLITTHREADCOUNT>>>(positions, positions2,indices,splitpoints,llimits,rlimits,tmpllimits,tmprlimits,llimits,rlimits,lrcount,splitcount,ysplit,threesplitvector,outputvector SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
  CHECKCUDAERROR
#endif
  cudasafe( cudaMemcpy(outputvectorlocal, outputvector, 2*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
  i=0;
#ifdef CHECKPARTITIONING
  CHECKPARTITIONINGSTRING1
  CHECKPARTITIONINGSTRINGSORT1
#endif
  //outputvectorlocal contains maximum number of elements in the coming split.
  //If this is small enough, move to singleblockpartition.
  //If outputvector[1]==0, one split did not move any elements at all, move to singleblockpartition
  //in this case as well, since singleblockpartition handles this more properly
  while(outputvectorlocal[0]>SPLITSHIFT&&outputvectorlocal[1]==0||i%2==1) { //always make an even number of splits from this point. This is to keep the results in the correct vector when moving to singleblockpartition
    cudasafe(cudaMemset(lrcount,0, 2*splitcount*sizeof(int)), "cudaMemset lrcount");
#ifdef CHECKPARTITIONING
    CHECKPARTITIONINGSTRING2
#endif

    //partitioning
    partitionsplit<<<blockcount,threadsperblock>>>(positions,positions2,indices,newindices,splitpoints,tmpllimits,tmprlimits,lrcount,splitcount,newpositions,newpositions2,ysplit,threesplitvector DEBUGVECTORSTRING3);
    CHECKCUDAERROR
#ifdef CHECKPARTITIONING
    CHECKPARTITIONINGSTRING3
#endif
    SWAP(indices,newindices,itmp);
    SWAP(positions,newpositions,ctmp);
    SWAP(positions2,newpositions2,ctmp);
    cudasafe(cudaMemset(outputvector,0, 2*sizeof(int)), "cudaMemset outputvector");

    //set up next one
#ifdef MEDIAN_OF_32
    preparesplit32<<<splitcount,PREPARESPLIT32THREADCOUNT>>>(positions, positions2,indices,splitpoints,tmpllimits,tmprlimits,newllimits,newrlimits,llimits,rlimits,lrcount,splitcount,ysplit,threesplitvector,outputvector SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
    CHECKCUDAERROR
#else
    preparesplit<<<(splitcount+PREPARESPLITTHREADCOUNT-1)/PREPARESPLITTHREADCOUNT,PREPARESPLITTHREADCOUNT>>>(positions, positions2,indices,splitpoints,tmpllimits,tmprlimits,newllimits,newrlimits,llimits,rlimits,lrcount,splitcount,ysplit,threesplitvector,outputvector SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
    CHECKCUDAERROR
#endif

    SWAP(newllimits,tmpllimits,itmp);
    SWAP(newrlimits,tmprlimits,itmp);
#ifdef CHECKPARTITIONING
    CHECKPARTITIONINGSTRING4
    CHECKPARTITIONINGSTRINGSORT1
#endif
    if(i%2==0) {
#ifdef CHECKPARTITIONING
      CHECKPARTITIONINGSTRING5
#endif

      //in this direction, copy the elements outside the split region to the output vectors
      splitoutsidecopy<<<blockcount,threadsperblock>>>(newllimits,newrlimits,tmpllimits,tmprlimits,positions,positions2,indices,newpositions,newpositions2,newindices,splitcount DEBUGVECTORSTRING3);
      CHECKCUDAERROR
#ifdef CHECKPARTITIONING
      CHECKPARTITIONINGSTRING6
#endif
    }
    cudasafe( cudaMemcpy(outputvectorlocal, outputvector, 2*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy outputvector" );
#ifdef CHECKPARTITIONING
    CHECKPARTITIONINGSTRING7
#endif
    i++;

  }
#ifdef CHECKPARTITIONING
  CHECKPARTITIONINGSTRING8
#endif
#ifdef CUDATIMESORT
  cudaEvent_t start;
  cudaEvent_t stop;
  cudasafe(cudaEventCreate(&start),"cudaEventCreate sort1start");
  cudasafe(cudaEventCreate(&stop),"cudaEventCreate sort1stop");
  cudasafe(cudaEventRecord(start,0),"cudaEventRecord sort1start");
#endif


  //finish by using singleblockpartition
  threadcount=threadsperblock;
  while(threadcount>N/splitcount&&threadcount>32)
    threadcount>>=1;
  if(threadcount<32) threadcount=32; //should not happen
  singleblockpartition<<<imin(splitcount,maxblockcount),threadcount,INDEXCACHELENGTH*threadcount*(2*sizeof(int)+sizeof(SORT_REAL))>>>(positions,positions2,indices,newindices,splitpoints,tmpllimits,tmprlimits,newllimits,newrlimits,llimits,rlimits,splitcount,newpositions,newpositions2,ysplit,1 SORTLIMITCALLINGSTRING DEBUGVECTORSTRING3);
  CHECKCUDAERROR
#ifdef CUDATIMESORT
  float elapsedtime;
  cudasafe(cudaEventRecord(stop),"cudaEventRecord fulltimestop");
  cudasafe(cudaEventSynchronize(stop),"cudaEventSynchronize fulltimestop");
  cudasafe(cudaEventElapsedTime(&elapsedtime,start,stop),"cudaEventElapsedTime");
  mexPrintf("singleblockpartition (%d blocks): (%d previous loops) %f\n",imin(splitcount,maxblockcount),i,elapsedtime/1000);
  cudasafe(cudaEventDestroy(stop),"cudaEventDestroy");
  cudasafe(cudaEventDestroy(start),"cudaEventDestroy");
#endif

#ifdef CHECKPARTITIONING
  CHECKPARTITIONINGSTRING9
#endif
}
/*------------------------------------------------------------------------*/
//debug function, checks if partitioning is valid or not
int checkpartitioning(const SORT_REAL* originalpositions,const SORT_REAL *originalpositions2,SORT_REAL *cudapositions,SORT_REAL *cudapositions2,SORT_REAL *splitpoints,int* cudaindices,int* oldcudaindices,int count,int* cudallimits,int* cudarlimits,int splitcount,int* cudaysplit,int printsuccess,int printfailure,SORT_DCMPLX* z,SORT_DCMPLX *d)
{
  SORT_REAL *positions=(SORT_REAL*)mxMalloc(count*sizeof(SORT_REAL));
  SORT_REAL *positions2=(SORT_REAL*)mxMalloc(count*sizeof(SORT_REAL));
  int *indices=(int*)mxMalloc(count*sizeof(int));
  int *oldindices=(int*)mxMalloc(count*sizeof(int));
  int *llimits=(int*)mxMalloc(2*splitcount*sizeof(int));
  int *rlimits=(int*)mxMalloc(2*splitcount*sizeof(int));
  int *ysplit=(int*)mxMalloc(splitcount*sizeof(int));
  int* buckets=(int*)mxCalloc(count, sizeof(int));
  SORT_DCMPLX *hz=(SORT_DCMPLX*)mxMalloc(2*splitcount*sizeof(SORT_DCMPLX));
  SORT_DCMPLX *hd=(SORT_DCMPLX*)mxMalloc(2*splitcount*sizeof(SORT_DCMPLX));
  SORT_REAL* hsplitpoints=(SORT_REAL*)mxMalloc(splitcount*sizeof(SORT_REAL));
  cudasafe( cudaMemcpy(positions, cudapositions, count*sizeof(SORT_REAL), cudaMemcpyDeviceToHost), "cudaMemcpy cudapositions" );
  cudasafe( cudaMemcpy(positions2, cudapositions2, count*sizeof(SORT_REAL), cudaMemcpyDeviceToHost), "cudaMemcpy cudapositions2" );
  cudasafe( cudaMemcpy(indices, cudaindices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy cudaindices" );
  cudasafe( cudaMemcpy(oldindices, oldcudaindices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy oldcudaindices" );
  cudasafe( cudaMemcpy(llimits, cudallimits, 2*splitcount*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy cudallimits" );
  cudasafe( cudaMemcpy(rlimits, cudarlimits, 2*splitcount*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy cudarlimits" );
  cudasafe( cudaMemcpy(ysplit, cudaysplit, splitcount*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy cudaysplit" );
  cudasafe( cudaMemcpy(hz, z, 2*splitcount*sizeof(SORT_DCMPLX), cudaMemcpyDeviceToHost), "cudaMemcpy z" );
  cudasafe( cudaMemcpy(hd, d, 2*splitcount*sizeof(SORT_DCMPLX), cudaMemcpyDeviceToHost), "cudaMemcpy d" );
  cudasafe( cudaMemcpy(hsplitpoints, splitpoints, splitcount*sizeof(SORT_REAL), cudaMemcpyDeviceToHost), "cudaMemcpy splitpoints" );
//     mexPrintf("running checkpartitioning\n");
  int returnvalue=0;
  int outside=0, failures=0, failurecount;
  SORT_REAL max, min, max2, min2;

  //start by checking that all elements exists in the vector
  for(int j=0;j<count;j++) {
    if(indices[j]>=count||indices[j]<0)
      outside++;
    else
      buckets[indices[j]]++;
  }
  for(int j=0;j<count;j++) {
    if(buckets[j]!=1)
      failures++;
  }
  if(failures==0) {
    if(printsuccess)
      mexPrintf("All elements accounted for\n");
  }
  else {
    returnvalue&=1;
    if(printfailure)
      mexPrintf("Invalid split, %d elements do not occur 1 time\n", failures);
    failurecount=0;
    //find where the elements were lost
    for(int i=0;i<splitcount;i++) {
      failures=0;
      memset(buckets, 0, count*sizeof(int));
      for(int j=llimits[2*i];j<rlimits[2*i+1];j++) {
        buckets[indices[j]]++;
        buckets[oldindices[j]]--;
      }
      for(int j=0;j<count;j++) {
        if(buckets[j]!=0) { //if element is missing, determine which split it belongs to
          failures++;
          if(printfailure&&failures<10) {
            if(buckets[j]<0) {
              mexPrintf("E%d missing. ", j);
              for(int k=llimits[2*i];k<rlimits[2*i+1];k++) {
                if(oldindices[k]==j)
                  mexPrintf("pos %d. ", k);
              }
            }
            if(buckets[j]>0) {
              mexPrintf("E%d extra. ", j);
              for(int k=llimits[2*i];k<rlimits[2*i+1];k++) {
                if(oldindices[k]==j)
                  mexPrintf("from %d. ", k);
                if(indices[k]==j)
                  mexPrintf("to %d. ", k);
              }
            }
          }
        }
      }
      if(failures!=0&&failurecount<20&&printfailure)
        mexPrintf("split %d invalid, %d elements do not occur 1 time\n", i, failures);
    }
    if(printfailure)
      mexPrintf("Invalid split, %d splits incorrect\n", failurecount);
  }
  if(outside!=0) {
    returnvalue&=2;
    if(printfailure)
      mexPrintf("Invalid split, %d indices outside interval\n", outside);
    mxFree(positions);
    mxFree(positions2);
    mxFree(indices);
    mxFree(llimits);
    mxFree(rlimits);
    mxFree(ysplit);
    mxFree(buckets);
    mxFree(hz);
    mxFree(hd);
    return returnvalue;  //abort in this case, considering that successive tests could cause segmentation fault
  }

  //check permutation of positions
  failures=0;
  for(int j=0;j<count;j++) {
    if(originalpositions[indices[j]]!=positions[j]) {
      if(failures<10&&printfailure) {
        int k=0;
        for(;j>llimits[k];k++);
        k--;
        mexPrintf("element %d (indices=%d) in split %d (%d to %d) (ysplit=%d) not premuted correctly for positions, permuted value: %e, original value: %e\n", j, indices[j], k, llimits[k], rlimits[k], ysplit[k>>1], positions[j], originalpositions[indices[j]]);
      }
      failures++;
    }
  }
  if(failures==0) {
    if(printsuccess)
      mexPrintf("Positions permuted correctly\n");
  }
  else {
    returnvalue&=4;
    if(printfailure)
      mexPrintf("Invalid permutation of positions, %d elements incorrect\n", failures);
  }

  //check permutation of positions2
  failures=0;
  for(int j=0;j<count;j++) {
    if(originalpositions2[indices[j]]!=positions2[j]) {
      if(failures<10&&printfailure) {
        int k=0;
        for(;j>llimits[k];k++);
        k--;
        mexPrintf("element %d (indices=%d) in split %d (%d to %d) (ysplit=%d) not premuted correctly for positions 2, permuted value: %e, original value: %e\n", j, indices[j], k, llimits[k], rlimits[k], ysplit[k>>1], positions2[j], originalpositions2[indices[j]]);
      }
      failures++;
    }
  }
  if(failures==0) {
    if(printsuccess)
      mexPrintf("Positions2 permuted correctly\n");
  }
  else {
    returnvalue&=8;
    if(printfailure)
      mexPrintf("Invalid permutation of positions2, %d elements incorrect\n", failures);
  }

  //check that all elements on left side are smaller than all on right side
  failures=0;
  int splitpointfailures=0;
  for(int k=0;k<splitcount;k++) {
    if(rlimits[2*k]>llimits[2*k]) {
      if(ysplit[k]) {
        max=positions2[llimits[2*k]];
        for(int j=llimits[2*k];j<rlimits[2*k];j++) {
          if(positions2[j]>max)
            max=positions2[j];
        }
      }
      else {
        max=positions[llimits[2*k]];
        for(int j=llimits[2*k];j<rlimits[2*k];j++) {
          if(positions[j]>max)
            max=positions[j];
        }
      }
    }
    if(rlimits[2*k+1]>llimits[2*k+1]) {
      if(ysplit[k]) {
        min=positions2[llimits[2*k+1]];
        for(int j=llimits[2*k+1];j<rlimits[2*k+1];j++) {
          if(positions2[j]<min)
            min=positions2[j];
        }
      }
      else {
        min=positions[llimits[2*k+1]];
        for(int j=llimits[2*k+1];j<rlimits[2*k+1];j++) {
          if(positions[j]<min)
            min=positions[j];
        }
      }
    }
    if(max>min&&rlimits[2*k]>llimits[2*k]&&rlimits[2*k+1]>llimits[2*k+1]) {
      if(printfailure&&failures<10)
        mexPrintf("split %d: min[%d,%d)=%.16e max[%d,%d)=%.16e\n", k, llimits[2*k+1], rlimits[2*k+1], min, llimits[2*k], rlimits[2*k], max);
      failures++;
    }
    //check that the value of the splitpoints are in between the two boxes
    if(max>hsplitpoints[k]&&rlimits[2*k]>llimits[2*k]||min<hsplitpoints[k]&&rlimits[2*k+1]>llimits[2*k+1]) {
      if(printfailure&&splitpointfailures<10)
        mexPrintf("split %d: min[%d,%d)=%.16e max[%d,%d)=%.16e splitpoint=%16e\n", k, llimits[2*k+1], rlimits[2*k+1], min, llimits[2*k], rlimits[2*k], max, hsplitpoints[k]);
      splitpointfailures++;
    }
  }
  //split properly. All elements on the correct side of the split
  if(failures==0) {
    if(printsuccess)
      mexPrintf("Split values ok\n");
  }
  else {
    returnvalue&=8;
    if(printfailure)
      mexPrintf("Invalid values in vector, split not correct, %d elements incorrect\n", failures);
  }
  if(splitpointfailures==0) {
    if(printsuccess)
      mexPrintf("Splitpoints ok\n");
  }
  else {
    returnvalue&=16;
    if(printfailure)
      mexPrintf("Splitpoints incorrect, %d elements incorrect\n", splitpointfailures);
  }
  for(int k=0;k<2*splitcount;k++) {
    if(rlimits[k]>llimits[k]) {
      max=positions[llimits[k]];
      max2=positions2[llimits[k]];
      min=positions[llimits[k]];
      min2=positions2[llimits[k]];
      for(int j=llimits[k];j<rlimits[k];j++) {
        if(positions[j]>max)
          max=positions[j];
        if(positions2[j]>max2)
          max2=positions2[j];
        if(positions[j]<min)
          min=positions[j];
        if(positions[j]<min2)
          min2=positions2[j];
      }

      if(creal(hz[k])-creal(hd[k])*1.0000001>min) {
        failures++;
        if(printfailure)
          mexPrintf("Box %d of %d size incorrect, real value %e smaller than box %e\n",k,2*splitcount,min, creal(hz[k])-creal(hd[k]));
      }
      if(creal(hz[k])+creal(hd[k])*1.0000001<max) {
        failures++;
        if(printfailure)
          mexPrintf("Box %d of %d size incorrect, real value %e larger than box %e\n",k,2*splitcount,max, creal(hz[k])+creal(hd[k]));
      }
      if(cimag(hz[k])-cimag(hd[k])*1.0000001>min2) {
        failures++;
        if(printfailure)
          mexPrintf("Box %d of %d size incorrect, imag value %e smaller than box %e\n",k,2*splitcount,min2, cimag(hz[k])-cimag(hd[k]));
      }
      if(cimag(hz[k])+cimag(hd[k])*1.0000001<max2) {
        failures++;
        if(printfailure)
          mexPrintf("Box %d of %d size incorrect, imag value %e larger than box %e\n",k,2*splitcount,max2, cimag(hz[k])+cimag(hd[k]));
      }
//       mexPrintf("B %d limits: [%.16e %.16e %.16e %.16e ] z0=%.16e + %.16e d0=%.16e + %.16e\n",k,min,max,min2,max2,creal(hz[k]),cimag(hz[k]),creal(hd[k]),cimag(hd[k]));
    }
  }

  mxFree(positions);
  mxFree(positions2);
  mxFree(indices);
  mxFree(oldindices);
  mxFree(llimits);
  mxFree(rlimits);
  mxFree(ysplit);
  mxFree(buckets);
  mxFree(hsplitpoints);
  mxFree(hz);
  mxFree(hd);
  return returnvalue;
}
/*------------------------------------------------------------------------*/
#endif /* defined(CUDASUPPORT) && defined(CUDASORT) */
//last two functions used by fmmsort in assymetric case even without CUDASORT
//as connectivity always is built on the GPU
#ifdef CUDASUPPORT
//the use of void* pointers instead of dcmplx is because of visual studio compatibility, where the built in dcmplx class does not work in cuda, and cuda files are compiled with different implementation of dcmplx compared to cpu files
void cudaCreateConnectivity(int *jcptr,int *kcptr,int *ir,
                            int *oldjcptr,int *oldkcptr,int *oldir,
                            int count,int maxm2p,
                            void *z,SORT_REAL *dabs,
                            SORT_REAL cutoff,int lastlevel,
                            int *outputvector DEBUGVECTORSTRING)
//wrapper since sort.cpp is not a cuda-file
{
  cudacreateconnectivity<<<imin(MAXCONNECTIVITYBLOCKS,(count+MAXCONNECTIVITYTHREADS-1)/MAXCONNECTIVITYTHREADS),MAXCONNECTIVITYTHREADS>>>(jcptr,kcptr,ir,oldjcptr,oldkcptr,oldir,count,maxm2p,(SORT_DCMPLX*)z,dabs,cutoff,lastlevel,outputvector DEBUGVECTORSTRING3);
  CHECKCUDAERROR
}
/*------------------------------------------------------------------------*/
void calcdabs(const void *d,SORT_REAL *dabs,int count)
{
  calculatedabs<<<imin(MAXCONNECTIVITYBLOCKS,(count+MAXCONNECTIVITYTHREADS-1)/MAXCONNECTIVITYTHREADS),MAXCONNECTIVITYTHREADS>>>((SORT_DCMPLX*)d,dabs,count);
  CHECKCUDAERROR
}
/*------------------------------------------------------------------------*/
void cumsumlist(int* oldjcptr,int* oldkcptr,int* jcptr,size_t count,cudavariables* GPUvars,int evalshift)
{
  size_t shift=CUMSUMSHIFTSTEP;
  size_t blockcount=(count+CUMSUMTHREADS-1)/CUMSUMTHREADS;
  size_t blockcounts[64/CUMSUMSHIFTSTEP+1]; //size=the limit where a 64 bit int would overflow anyway
  size_t shiftfactor=1<<CUMSUMSHIFTSTEP;
  size_t fullshiftfactor=1<<shift;
  size_t i=1;
  blockcounts[0]=blockcount;
  if(GPUvars->evalonly) {
    cumsuminitevalonly<<<imin(blockcount,CUMSUMMAXBLOCKS),CUMSUMTHREADS>>>(oldjcptr,oldkcptr,jcptr,count,blockcount,GPUvars->jxptr,evalshift);
    CHECKCUDAERROR
  }
  else {
    cumsuminit<<<imin(blockcount,CUMSUMMAXBLOCKS),CUMSUMTHREADS>>>(oldjcptr,oldkcptr,jcptr,count,blockcount,GPUvars->ixptr,GPUvars->jxptr,evalshift);
    CHECKCUDAERROR
  }
  blockcount=(blockcount+shiftfactor-2)/shiftfactor;
  while(blockcount>=1) {
    blockcounts[i]=blockcount;
//     mexPrintf("blockcount=%d shift=%d shiftfactor=%d fullshiftfactor=%d\n",blockcount,shift,1<<shiftfactor,fullshiftfactor);
    cumsumpass1<<<imin(blockcount,CUMSUMMAXBLOCKS),CUMSUMTHREADS>>>(jcptr+fullshiftfactor-1,count-fullshiftfactor+1,blockcount,shift);
    CHECKCUDAERROR
//     cudasafe(cudaThreadSynchronize(),"cudaThreadSynchronize");
    shift+=CUMSUMSHIFTSTEP;
    fullshiftfactor=1<<shift;
    blockcount=(blockcount+shiftfactor-2)/shiftfactor;
    i++;
  }
  i--;
  shift-=CUMSUMSHIFTSTEP*2;
  while(i>0) {
    fullshiftfactor=1<<shift;
    cumsumpass2<<<imin(blockcounts[i-1],CUMSUMMAXBLOCKS),CUMSUMTHREADS>>>(jcptr+fullshiftfactor-1,count-fullshiftfactor+1,blockcounts[i-1],shift/*,debugvector*/);
    CHECKCUDAERROR
//     cudasafe(cudaThreadSynchronize(),"cudaThreadSynchronize");
    shift-=CUMSUMSHIFTSTEP;
    i--;
  }
}
/*------------------------------------------------------------------------*/
#ifdef SORTLIMIT
void calculatesortlimits(float *distlimit,float *disttarget,float *sortlimit,float* input,float currentlevel,float maxlevel)
{
  float tmpdistlimit,tmpdisttarget,tmpsortlimit;
  float interppos=currentlevel/(float)maxlevel;
  tmpsortlimit=input[0]*(1-interppos)+input[1]*interppos;
  tmpdistlimit=input[2]*(1-interppos)+input[3]*interppos;
  tmpdisttarget=input[4]*(1-interppos)+input[5]*interppos;
  //tmpsortlimit needs to be atleast 0, larger than 1 gives no effect, but should work.
  //put all variables in the region [0 1] as it is here they have effect.
  if(tmpsortlimit<0)
    tmpsortlimit=0;
  if(tmpdisttarget<0)
    tmpdisttarget=0;
  if(tmpdisttarget>1)
    tmpdisttarget=1;
  if(tmpdistlimit<=0) { //split in position centre
    tmpdistlimit=0;
    tmpsortlimit=2; //disable this for speed issues
  }
  if(tmpdistlimit>=1)
    tmpdistlimit=1;
  if(tmpdisttarget>tmpdistlimit) //not reasonable, and not implemented in the code, as sorting could fail otherwise
    tmpdisttarget=tmpdistlimit;
  *distlimit=tmpdistlimit;
  *disttarget=tmpdisttarget;
  *sortlimit=tmpsortlimit;
//   mexPrintf("interppos=%f distlimit=%f disttarget=%f sortlimit=%f\n",interppos,tmpdistlimit,tmpdisttarget,tmpsortlimit);
}
#endif
double calculateetalimits(double *eta,double currentlevel,double maxlevel)
{
  double ret;
  double interppos=currentlevel/maxlevel;
  ret=eta[0]*(1-interppos)+eta[1]*interppos;
  if(ret<0) //eta should be larger than 0.
    ret=0;
  return ret;
}
#endif
