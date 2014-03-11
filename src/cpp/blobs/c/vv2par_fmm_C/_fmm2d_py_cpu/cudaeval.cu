/* cudaeval.cu -- CUDA host code for direct evaluation. */

/* S. Engblom and A. Goude 2011-10-19 */

#ifdef CUDASUPPORT

#include "cuda.h"
#include "cudaeval.h"

/*------------------------------------------------------------------------*/
int cudaalloccount; // remembers the number of cuda allocations

cudaError_t cudaMallocDebug(void **addr,size_t size) 
// GPU memory allocation debug, can track memory leaks
{
  cudaalloccount++;
  cudaError_t ret = cudaMalloc(addr,size);
  if (*addr == NULL) cudaalloccount--;
  return ret;
}
/*------------------------------------------------------------------------*/
cudaError_t cudaFreeDebug(void* addr)
// corresponding to CudaMallocDebug, all allocations should be
// performed through these calls
{
  if (addr != NULL) cudaalloccount--;
  return cudaFree(addr);
}
/*------------------------------------------------------------------------*/
void printcudaalloccount(char *message)
// for tracking memory leaks
{
  mexPrintf("%s %d memory allocations on GPU.\n",message,cudaalloccount);
}
/*------------------------------------------------------------------------*/
int getcudaalloccount()
// since the variable currently only is declared in this file.
{
  return cudaalloccount;
}
/*------------------------------------------------------------------------*/
void resetalloccount()
// since the variable currently only is declared in this file.
{
  cudaalloccount=0;;
}
/*------------------------------------------------------------------------*/
void cudasafe2(cudaError_t error,const char *message,const char* command,const char *file,int line)
// checks if a cuda function succeded, if not, print error message and
// abort (Cuda functions fail a lot more often than normal cpu
// functions. it is important to check the return values)
{
  char errormessage[1000];
  if (error != cudaSuccess) {
    snprintf(errormessage,999,"%s : %s (%i) \nCommand: %s\nFile %s at line %d",message,
             cudaGetErrorString(error),error,command,file,line);
    cudaGetLastError(); // remove the error message
    cudaThreadExit();
    resetalloccount();//cudaThreadExit() should deallocate variables
    mexErrMsgTxt(errormessage); 
  }
}
/*------------------------------------------------------------------------*/
void cudasafeMalloc(void** addr,size_t size)
// safe memory allocation, cudasafe(cudaMalloc(...),"message") can
// also be used
{
  char errormessage[1000];
  cudaError_t error = cudaMallocDebug(addr,size);
  if(error != cudaSuccess) {
    sprintf(errormessage,"cudaMalloc (%u bytes): %s (%i) \n",
            (unsigned)size,cudaGetErrorString(error),error);
    cudaGetLastError(); //remove the error message
    cudaThreadExit();
    mexErrMsgTxt(errormessage); 
  }
}
/*------------------------------------------------------------------------*/
void checkcudaerror(const char *message)
{
  //only check for error and abort if error is found
  char errormessage[1000];
  cudaError_t error=cudaGetLastError();
  if(error!=cudaSuccess) {
    sprintf(errormessage,"%s : %s (%i) \n",message,
            cudaGetErrorString(error),error);
    cudaThreadExit();
    mexErrMsgTxt(errormessage); 
  }
}
/*------------------------------------------------------------------------*/
void cudareset()
{
  cudaThreadExit();
}
/*------------------------------------------------------------------------*/
void cudastart()
{
  char errormessage[1000];
  cudaError_t error=cudaGetLastError();
  if(error!=cudaSuccess) {
    sprintf(errormessage,"Cuda initialization error: %s (%i) \n",
            cudaGetErrorString(error),error);
    cudaThreadExit();
    mexWarnMsgTxt(errormessage); 
  }
}
/*------------------------------------------------------------------------*/

//textures have to be global
texture<int2,1> ter;
texture<int2,1> tei;
texture<int2,1> tzr;
texture<int2,1> tzi;
texture<int2,1> tmr;
texture<int2,1> tmi;
texture<int> tixptr;
texture<int> tjxptr;
#if defined(MULTIPOLEEVALCUDA) || defined(MULTIPOLEINITCUDA)
texture<float> tz0;
// texture<int2,1> tcoeff1;
// texture<int2,1> tcoeff2;
texture<int4,1> tcoeff1d;
texture<int4,1> tcoeff2d;
#endif

#define MAXBLOCKS 16384
#define MAXTHREADSEVALUATION 64 //max threads for fmm direct evaluation and multipole evaluation
#define DIRECTMAXTHREADS 256 //max threads for tol=0, i.e. direct summation
#define REORDERMAXTHREADS 256 //max threads for CUDA reordering
#define M2PSMAXTHREADS 32 //number of threads in each block for M2PS interaction
#define P2PMAXTHREADS 64 //number of threads in each block for P2P interaction (needs to be atleast P2PMINTHREADS)
#define P2PMINTHREADS 32 //number of threads in each block for P2P interaction (lowest value)
#define SHIFT_M2M_MAXTHREADS 64 //number of threads in each block for P2P interaction (needs to be atleast SHIFT_M2M_MINTHREADS)
#define SHIFT_M2M_MINTHREADS 32 //number of threads in each block for P2P interaction (lowest value)

//cuda kernels has to be declared in same file as they are calles
#include "cudaevalkernels.h"

//function prototype
void launchdirectkernel(const MPexp *This,int pot,const double *mi,
                        cudavariables *GPUvars,
                        double *GPUqr,double *GPUqi,
                        SMOOTHER smooth,double shape,
                        double scale,double cutoff);


//Allocates GPU memory and copies data to the GPU. For the full GPU algorithm
//all GPU data is allocated and copied inside this function.
void MPexp_setup_cuda(const MPexp *This,
                     const double *er,const double *ei,
                     const double *zr,const double *zi,
                     const double *mr,const double *mi,
                     const double *pr,cudavariables *GPUvars,int N,int NE,
                     int pot,SMOOTHER smooth,
                     double xopt,double cutoff,bool cont,
                     double *timing,int printtime)
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  int Nt=(Nf-1)/3+Nf;
#ifndef CUDASORT
  const int *restrict ix = This->ix, *restrict jx = This->jx;
#endif
  checkcudaerror("Check eval_cuda start");
    
  TimeEventRecord(GPUvars->sortstop);
  TimeCreateAndStart(GPUvars->fulltimestart,GPUvars->fulltimestop);
            
  //start by allocating variables on GPU. If CUDASORT is used, some of them are already allocated
  if(NE) {
#ifndef CUDASORT
    cudasafeMalloc((void**)&GPUvars->er,NE*sizeof(double));
    cudasafeMalloc((void**)&GPUvars->ei,NE*sizeof(double));
#endif
    cudasafeMalloc((void**)&GPUvars->qr,NE*sizeof(double));
    if(pot!=0||mi!=NULL)
      cudasafeMalloc((void**)&GPUvars->qi,NE*sizeof(double));
  }
  else {
    GPUvars->er=NULL;
    GPUvars->ei=NULL;
  }
  if(pr!=NULL) {
    cudasafeMalloc((void**)&GPUvars->pr,N*sizeof(double));
    if(pot!=0||mi!=NULL)
      cudasafeMalloc((void**)&GPUvars->pi,N*sizeof(double));
  }
#ifndef CUDASORT
  cudasafeMalloc((void**)&GPUvars->zr,N*sizeof(double));
  cudasafeMalloc((void**)&GPUvars->zi,N*sizeof(double));
#endif
  cudasafeMalloc((void**)&GPUvars->mr,N*sizeof(double));
    
  if(mi!=NULL) {
    cudasafeMalloc((void**)&GPUvars->mi,N*sizeof(double));
  }
  else
    GPUvars->mi=NULL;
    
  //reorder the elements on the GPU according to ix/jx to avoid random accesses to device memory. If cudasort, some variables are already reordered
  TimeCreateAndStart(GPUvars->reorder1start,GPUvars->reorder1stop);
            
  double* mrtmp,*mitmp;
  //If cudasort is active, the elements are reordered during the sorting, otherwise, do it here
#ifndef CUDASORT
  double* zrtmp,*zitmp;
  int* cudaix,*cudajx;
  if(NE) { //reorder the evaluation points
    if(pot==0&&mi==NULL)
      cudasafeMalloc((void**)&GPUvars->qi, NE*sizeof(double));
    cudasafeMalloc((void**)&cudajx, NE*sizeof(int));
    cudasafe( cudaMemcpy(cudajx, jx, NE*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy jx" );
    cudasafe( cudaMemcpy(GPUvars->qr, er, NE*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy ertmp" );
    cudasafe( cudaMemcpy(GPUvars->qi, ei, NE*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy eitmp" );
        
    reorder_cuda_2<<<imin(MAXBLOCKS, (NE+REORDERMAXTHREADS-1)/REORDERMAXTHREADS), REORDERMAXTHREADS>>>(GPUvars->er, GPUvars->ei, GPUvars->qr, GPUvars->qi, cudajx, NE);
    CHECKCUDAERROR
    if(pot==0&&mi==NULL)
      cudaFreeDebug(GPUvars->qi); //used as temporary variable, if not used anymore, cleanup
    GPUvars->jx=cudajx;
  }
  //reorder the potential points
  if(pr==NULL) { //the case no evaluation should happen on the potential points, in this case, do a proper cleanup afterwards
    cudasafeMalloc((void**)&zrtmp,N*sizeof(double));
    cudasafeMalloc((void**)&zitmp,N*sizeof(double));
    cudasafeMalloc((void**)&mrtmp,N*sizeof(double));
    cudasafeMalloc((void**)&cudaix,N*sizeof(int));
        
    if(mi!=NULL) {
      cudasafeMalloc((void**)&mitmp,N*sizeof(double));
      cudasafe( cudaMemcpy(mitmp, mi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mitmp" );
    }
    cudasafe( cudaMemcpy(zrtmp, zr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zrtmp" );
    cudasafe( cudaMemcpy(zitmp, zi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zitmp" );
    cudasafe( cudaMemcpy(mrtmp, mr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mrtmp" );
    cudasafe( cudaMemcpy(cudaix, ix, N*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy ix" );
        
    if(mi!=NULL) {
      reorder_cuda_4<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->zr,GPUvars->zi,GPUvars->mr,GPUvars->mi,zrtmp,zitmp,mrtmp,mitmp,cudaix,N);
      CHECKCUDAERROR
    }
    else {
      reorder_cuda_3<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->zr,GPUvars->zi,GPUvars->mr,zrtmp,zitmp,mrtmp,cudaix,N);
      CHECKCUDAERROR
    }
    cudaFreeDebug(zrtmp);
    cudaFreeDebug(zitmp);
    cudaFreeDebug(mrtmp);
        
    if(mi!=NULL)
      cudaFreeDebug(mitmp);
    cudaFreeDebug(cudaix);
  }
  else { //evaluation will happen at the potential points. Now, the arrays can be allocated and kept on the GPU for future use
    if(pot==0&&mi==NULL)
      cudasafeMalloc((void**)&GPUvars->pi,N*sizeof(double));
    cudasafeMalloc((void**)&mrtmp,N*sizeof(double));
    cudasafeMalloc((void**)&cudaix,N*sizeof(int));
    if(mi!=NULL) {
      cudasafeMalloc((void**)&mitmp,N*sizeof(double));
      cudasafe( cudaMemcpy(mitmp, mi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mitmp" );
    }
    cudasafe( cudaMemcpy(GPUvars->pr, zr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zrtmp" );
    cudasafe( cudaMemcpy(GPUvars->pi, zi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zitmp" );
    cudasafe( cudaMemcpy(mrtmp, mr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mrtmp" );
    cudasafe( cudaMemcpy(cudaix, ix, N*sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy ix" );
    if(mi!=NULL) {
      reorder_cuda_4<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->zr,GPUvars->zi,GPUvars->mr,GPUvars->mi,GPUvars->pr,GPUvars->pi,mrtmp,mitmp,cudaix,N);
      CHECKCUDAERROR
    }
    else {
      reorder_cuda_3<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->zr,GPUvars->zi,GPUvars->mr,GPUvars->pr,GPUvars->pi,mrtmp,cudaix,N);
      CHECKCUDAERROR
    }
    cudaFreeDebug(mrtmp);
    if(pot==0&&mi==NULL)
      cudaFreeDebug(GPUvars->pi);
    if(mi!=NULL)
      cudaFreeDebug(mitmp);
    GPUvars->ix=cudaix;
  }
#else 
#ifdef FLOATSORT //data transferred to CPU, but not ordered
  double *zrtmp,*zitmp;
  cudasafeMalloc((void**)&mrtmp,N*sizeof(double));
  cudasafeMalloc((void**)&zrtmp,N*sizeof(double));
  cudasafeMalloc((void**)&zitmp,N*sizeof(double));
  cudasafe( cudaMemcpy(mrtmp, mr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mrtmp" );
  if(mi!=NULL) { //imaginary strength?
    cudasafeMalloc((void**)&mitmp, N*sizeof(double));
    cudasafe( cudaMemcpy(mitmp, mi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mitmp" );
    reorder_cuda_4<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->mr,GPUvars->mi,zrtmp,zitmp,mrtmp,mitmp,GPUvars->zr,GPUvars->zi,GPUvars->ix,N);
    CHECKCUDAERROR
    cudaFreeDebug(mitmp);
  }
  else {
    reorder_cuda_3<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->mr,zrtmp,zitmp,mrtmp,GPUvars->zr,GPUvars->zi,GPUvars->ix,N);
    CHECKCUDAERROR
  }
  double* zrbackup=GPUvars->zr;
  double* zibackup=GPUvars->zi;
  GPUvars->zr=zrtmp;
  GPUvars->zi=zitmp;
  if(NE) {
    if(NE!=N) {
      cudaFreeDebug(zrbackup);
      cudaFreeDebug(zibackup);
      cudasafeMalloc((void**)&zrbackup,NE*sizeof(double));
      cudasafeMalloc((void**)&zibackup,NE*sizeof(double));
    }
    reorder_cuda_2<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(zrbackup,zibackup,GPUvars->er,GPUvars->ei,GPUvars->jx,NE);
    CHECKCUDAERROR
    cudaFreeDebug(GPUvars->er);
    GPUvars->er=zrbackup;
    cudaFreeDebug(GPUvars->ei);
    GPUvars->ei=zibackup;
  }
  else {
    cudaFreeDebug(zrbackup);
    cudaFreeDebug(zibackup);
  }
  cudaFreeDebug(mrtmp);
  if(pr==NULL) { //no evaluation on potential points, cleanup
    cudaFreeDebug(GPUvars->ix);
    GPUvars->ix=0;
  }
#else//CUDASORT, Now, only the strengths has to be reordered, the positions are reordered by the sorting code
  cudasafeMalloc((void**)&mrtmp,N*sizeof(double));
  cudasafe( cudaMemcpy(mrtmp, mr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mrtmp" );
  if(mi!=NULL) { //imaginary strength?
    cudasafeMalloc((void**)&mitmp, N*sizeof(double));
    cudasafe( cudaMemcpy(mitmp, mi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mitmp" );
    reorder_cuda_2<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->mr,GPUvars->mi,mrtmp,mitmp,GPUvars->ix,N);
    CHECKCUDAERROR
    cudaFreeDebug(mitmp);
  }
  else {
    reorder_cuda_1<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->mr,mrtmp,GPUvars->ix,N);
    CHECKCUDAERROR
  }
  cudaFreeDebug(mrtmp);
  if(pr==NULL) { //no evaluation on potential points, cleanup
    cudaFreeDebug(GPUvars->ix);
    GPUvars->ix=0;
  }
#endif
#endif
  TimeEventRecord(GPUvars->reorder1stop);
  
#ifdef CUDADEBUGVECTOR
  GPUvars->N=N;//debugging purposes, to always have the values available
  GPUvars->NE=NE;
#endif
  GPUvars->Nf=Nf;
  
  //Allocate output
  if(NE) {
    cudasafe( cudaMemset(GPUvars->qr, 0, NE*sizeof(double)), "cudaMemset qr" );
    if(pot!=0||mi!=NULL)
      cudasafe( cudaMemset(GPUvars->qi, 0, NE*sizeof(double)), "cudaMemset qi" );
  }
  if(pr!=NULL) {
    cudasafe( cudaMemset(GPUvars->pr, 0, N*sizeof(double)), "cudaMemset qr" );
    if(pot!=0||mi!=NULL)
      cudasafe( cudaMemset(GPUvars->pi, 0, N*sizeof(double)), "cudaMemset qi" );
  }
  
  //use textures for faster memory access
  if(mi!=NULL) {
    cudasafe(cudaBindTexture(NULL,tmi,GPUvars->mi,N*sizeof(double)),"cudaBindTexture tmi");
  }
  cudasafe(cudaBindTexture(NULL,tzr,GPUvars->zr,N*sizeof(double)),"cudaBindTexture tzr");
  cudasafe(cudaBindTexture(NULL,tzi,GPUvars->zi,N*sizeof(double)),"cudaBindTexture tzi");
  cudasafe(cudaBindTexture(NULL,tmr,GPUvars->mr,N*sizeof(double)),"cudaBindTexture tmr");
    
#if defined(CUDADEBUGVECTOR) &&!defined(CUDASORT) //debug code
  cudasafe(cudaMallocDebug((void**)&GPUvars->debugvector,imax((NE<N?N:NE),2000)*sizeof(double)),"cudaMalloc");
  cudasafe(cudaMemset(GPUvars->debugvector,0,imax((NE<N?N:NE),2000)*sizeof(double)),"cudaMemset");
#endif
  
  cudasafe(cudaBindTexture(NULL,tixptr,GPUvars->ixptr,(Nf+1)*sizeof(int)),"cudaBindTexture ixptr");
#if defined(MULTIPOLEEVALCUDA) || defined(MULTIPOLEINITCUDA)
  cudasafeMalloc((void**)&GPUvars->coeff1,(This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudasafeMalloc((void**)&GPUvars->coeff2,(This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudasafe(cudaBindTexture(NULL,tz0,GPUvars->z0,Nt*sizeof(SORT_DCMPLX)),"cudaBindTexture z0");
//   cudasafe(cudaBindTexture(NULL,tcoeff1,GPUvars->coeff1,(This->pcoeff+1)*Nt*sizeof(dcmplx)),"cudaBindTexture coeff1");
//   cudasafe(cudaBindTexture(NULL,tcoeff2,GPUvars->coeff2,(This->pcoeff+1)*Nt*sizeof(dcmplx)),"cudaBindTexture coeff2");
  cudasafe(cudaBindTexture(NULL,tcoeff1d,GPUvars->coeff1,(This->pcoeff+1)*Nt*sizeof(dcmplx)),"cudaBindTexture coeff1d");
  cudasafe(cudaBindTexture(NULL,tcoeff2d,GPUvars->coeff2,(This->pcoeff+1)*Nt*sizeof(dcmplx)),"cudaBindTexture coeff2d");
#ifdef MULTIPOLEINITCUDA
  cudasafe( cudaMemset(GPUvars->coeff1, 0, (This->pcoeff+1)*Nt*sizeof(dcmplx)), "cudaMemset coeff1" );
  cudasafe( cudaMemset(GPUvars->coeff2, 0, (This->pcoeff+1)*Nt*sizeof(dcmplx)), "cudaMemset coeff2" );
#endif /*MULTIPOLEINITCUDA*/
#endif /*defined(MULTIPOLEEVALCUDA) || defined(MULTIPOLEINITCUDA)*/
  
#ifdef MULTIPOLESHIFTCUDA
//for m2ps, save an array of the connecticity pointers on the GPU memory to
//allow for that one call to m2ps can handle all levels
  int** jcptr2,**kcptr2,**ir2;
  jcptr2=(int**)mxMalloc((This->nlevel+1)*sizeof(int*));
  kcptr2=(int**)mxMalloc((This->nlevel+1)*sizeof(int*));
  ir2=(int**)mxMalloc((This->nlevel+1)*sizeof(int*));
  cudasafe(cudaMallocDebug((void**)&GPUvars->jcptr2,(This->nlevel+1)*sizeof(int*)),"cudaMalloc jcptr2");
  cudasafe(cudaMallocDebug((void**)&GPUvars->kcptr2,(This->nlevel+1)*sizeof(int*)),"cudaMalloc kcptr2");
  cudasafe(cudaMallocDebug((void**)&GPUvars->ir2,(This->nlevel+1)*sizeof(int*)),"cudaMalloc ir2");
  for(int j=0;j<=This->nlevel;j++) {
    jcptr2[j]=GPUvars->connect[j].jcptr;
    kcptr2[j]=GPUvars->connect[j].kcptr;
    ir2[j]=GPUvars->connect[j].ir;
  }
  cudasafe( cudaMemcpy(GPUvars->jcptr2, jcptr2, (This->nlevel+1)*sizeof(int*), cudaMemcpyHostToDevice), "cudaMemcpy jcptr2" );
  cudasafe( cudaMemcpy(GPUvars->kcptr2, kcptr2, (This->nlevel+1)*sizeof(int*), cudaMemcpyHostToDevice), "cudaMemcpy kcptr2" );
  cudasafe( cudaMemcpy(GPUvars->ir2, ir2, (This->nlevel+1)*sizeof(int*), cudaMemcpyHostToDevice), "cudaMemcpy ir2" );
#endif
}

/*------------------------------------------------------------------------*/
//function to start the GPU direct summation calculations, should be
//called as soon as the sorting is completed
void MPexp_eval_cuda(const MPexp *This,
                     const double *er,const double *ei,
                     const double *zr,const double *zi,
                     const double *mr,const double *mi,
                     const double *pr,cudavariables *GPUvars,int N,int NE,
                     int pot,SMOOTHER smooth,
                     double xopt,double cutoff,bool cont,
                     double *timing,int printtime)
/* Evaluation of the potential at the points (er,ei). Function is
   critical to performance. */
{
  double shape, scale;
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  MPexp_smooth_(pot, smooth, xopt, cont, &cutoff, &shape, &scale);
  
  MPexp_setup_cuda(This,er,ei,zr,zi,mr,mi,pr,GPUvars,
                   N,NE,pot,smooth,xopt,cutoff,cont,
                   timing,printtime);
  if(pot==0&&mi==NULL) { //mark in the GPUvars structure if evaluations should be done or not
    if(NE)
      GPUvars->qi=NULL;
    if(pr!=NULL)
      GPUvars->pi=NULL;
  }
  //     checkcudaerror("Check before Init_cuda");  
    
  //Allocations and memory copies to cuda done. Start the algorithm
#if defined( MULTIPOLEINITCUDA) //make multipole init first, to enable CPU to continue
  MPexp_init_cuda(This,GPUvars,pot,mi!=NULL,timing,printtime);
#endif

  //launch kernel for direct evaluation
  checkcudaerror("Check before direct interact");
  //start the timing of the direct evaluation algorithm
  TimeCreateAndStart(GPUvars->start,GPUvars->stop);
  //one case for evaluation in potential points and one for evaluation points
  if(pr!=NULL) {
    cudasafe(cudaBindTexture(NULL,ter,GPUvars->zr,N*sizeof(double)),"cudaBindTexture tzr");
    cudasafe(cudaBindTexture(NULL,tei,GPUvars->zi,N*sizeof(double)),"cudaBindTexture tzi");
    cudasafe(cudaBindTexture(NULL,tjxptr,GPUvars->ixptr,(Nf+1)*sizeof(int)),"cudaBindTexture jxptr");
    launchdirectkernel(This,pot,mi,GPUvars,GPUvars->pr,GPUvars->pi,smooth,shape,scale,cutoff);
  }
  if(NE) {
    //move the textures to the evaluation points before calling the kernel again
    if(pr!=NULL) {
      cudaUnbindTexture(ter);
      cudaUnbindTexture(tei);
      cudaUnbindTexture(tjxptr);
    }
    cudasafe(cudaBindTexture(NULL,ter,GPUvars->er,NE*sizeof(double)),"cudaBindTexture ter");
    cudasafe(cudaBindTexture(NULL,tei,GPUvars->ei,NE*sizeof(double)),"cudaBindTexture tei");
    cudasafe(cudaBindTexture(NULL,tjxptr,GPUvars->jxptr,(Nf+1)*sizeof(int)),"cudaBindTexture jxptr");
    launchdirectkernel(This,pot,mi,GPUvars,GPUvars->qr,GPUvars->qi,smooth,shape,scale,cutoff);
  }

  checkcudaerror("Kernel launch error direct interact mpexp_eval_cuda");
  TimeEventRecord(GPUvars->stop);

}
/*------------------------------------------------------------------------*/
//chooses the correct direct evaluation kernel and launches GPU calculation.
void launchdirectkernel(const MPexp *This,int pot,const double *mi,cudavariables *GPUvars,
                        double *GPUqr,double *GPUqi,
                        SMOOTHER smooth,double shape,double scale,double cutoff)
{
  if(pot==1) {
    if(mi==NULL) {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        pot1_nomi_dirac_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi,GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir DEBUGVECTORSTRING2);
        break;
      case RANKINE:

        pot1_nomi_rankine_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape DEBUGVECTORSTRING2);
        break;
      case SCULLY:
        pot1_nomi_scully_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      case OSEEN:
        pot1_nomi_oseen_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      }
    }
    else {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        pot1_mi_dirac_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi,GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir DEBUGVECTORSTRING2);
        break;
      case RANKINE:

        pot1_mi_rankine_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape DEBUGVECTORSTRING2);
        break;
      case SCULLY:
        pot1_mi_scully_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      case OSEEN:
        pot1_mi_oseen_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      }
    }
  }
  else {
    if(mi==NULL) {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        pot0_nomi_dirac_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr,GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir DEBUGVECTORSTRING2);
        break;
      case RANKINE:
        pot0_nomi_rankine_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, cutoff DEBUGVECTORSTRING2);
        break;
      case SCULLY:
        pot0_nomi_scully_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      case OSEEN:
        pot0_nomi_oseen_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      }
    }
    else {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        pot0_mi_dirac_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi,GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir DEBUGVECTORSTRING2);
        break;
      case RANKINE:

        pot0_mi_rankine_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, cutoff DEBUGVECTORSTRING2);
        break;
      case SCULLY:
        pot0_mi_scully_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      case OSEEN:
        pot0_mi_oseen_synchronized<<<imin(MAXBLOCKS, GPUvars->Nf), MAXTHREADSEVALUATION>>>(GPUqr, GPUqi, GPUvars->Nf,GPUvars->connect[This->nlevel].jcptr,&GPUvars->connect[This->nlevel].kcptr[GPUvars->Nf],GPUvars->connect[This->nlevel].ir, shape, scale, cutoff DEBUGVECTORSTRING2);
        break;
      }
    }
  }
  CHECKCUDAERROR
  checkcudaerror("Kernel launch error direct interact direct_interact_kernel");
}
/*------------------------------------------------------------------------*/
//function to collect results and performes cleanup, should be called at the end of the algorithm
void MPexp_eval_cuda_result(double *qr,double *qi,
                            double *pr,double *pi,
                            int *ix,int *jx,int N,int NE,
                            int pot,cudavariables *GPUvars,
                            const MPexp *This,
                            double *timing,int printtime)
{
  //copy results back to normal memory
#ifdef MULTIPOLEEVALCUDA //reordercuda only implemented for multipoleevalcuda, where the cpu does not change the velocities
  TimeCreateAndStart(GPUvars->reorder2start,GPUvars->reorder2stop);
  //reorder the values on the CPU and copy them back to the CPU
  if(NE) { //For evaluation points
    checkcudaerror("Kernel launch error");
    if(qi!=NULL) {
      reorder_cuda_2_inv<<<imin(MAXBLOCKS, (NE+REORDERMAXTHREADS-1)/REORDERMAXTHREADS), REORDERMAXTHREADS>>>(GPUvars->er, GPUvars->ei, GPUvars->qr, GPUvars->qi, GPUvars->jx, NE);
      CHECKCUDAERROR
      cudasafe( cudaMemcpy(qr, GPUvars->er, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qr" );
      cudasafe( cudaMemcpy(qi, GPUvars->ei, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qi" );
    }
    else {
      reorder_cuda_1_inv<<<imin(MAXBLOCKS,(NE+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->er,GPUvars->qr,GPUvars->jx,NE);
      CHECKCUDAERROR
      cudasafe( cudaMemcpy(qr, GPUvars->er, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qr" );
    }
    checkcudaerror("Kernel launch error reorder");
  }
  if(pr!=NULL) { //for potential points
    checkcudaerror("Kernel launch error");
    if(pi!=NULL) {
      reorder_cuda_2_inv<<<imin(MAXBLOCKS, (N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS), REORDERMAXTHREADS>>>(GPUvars->zr, GPUvars->zi, GPUvars->pr, GPUvars->pi, GPUvars->ix, N);
      CHECKCUDAERROR
      cudasafe( cudaMemcpy(pr, GPUvars->zr, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qr" );
      cudasafe( cudaMemcpy(pi, GPUvars->zi, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qi" );
    }
    else {
      reorder_cuda_1_inv<<<imin(MAXBLOCKS,(N+REORDERMAXTHREADS-1)/REORDERMAXTHREADS),REORDERMAXTHREADS>>>(GPUvars->zr,GPUvars->pr,GPUvars->ix,N);
      CHECKCUDAERROR
      cudasafe( cudaMemcpy(pr, GPUvars->zr, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qr" );
    }
    checkcudaerror("Kernel launch error reorder");
  }
    
  TimeEventRecord(GPUvars->reorder2stop);
#else //special case if MULTIPOLEEVALCUDA is false when results from CPU and GPU calculations has to be compined.
 double *fqr,*fqi,*fpr,*fpi;
  if(NE) {
    fqr=(double*)mxMalloc(NE*sizeof(double));
    cudasafe( cudaMemcpy(fqr, GPUvars->qr, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qr" );
    if(qi!=NULL) {
      fqi=(double*)mxMalloc(NE*sizeof(double));
      cudasafe( cudaMemcpy(fqi, GPUvars->qi, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qi" );
    }
  }
  if(pr!=NULL) {
    fpr=(double*)mxMalloc(N*sizeof(double));
    cudasafe( cudaMemcpy(fpr, GPUvars->pr, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy pr" );
    if(pi!=NULL) {
      fpi=(double*)mxMalloc(N*sizeof(double));
      cudasafe( cudaMemcpy(fpi, GPUvars->pi, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy pi" );
    }
  }
  if(NE) {
    if(qi!=NULL) {
      for(int i=0;i<NE;i++) {
        qr[jx[i]]+=(double)fqr[i];
        qi[jx[i]]+=(double)fqi[i];
      }
      mxFree(fqi);
    }
    else {
      for(int i=0;i<NE;i++) {
        qr[jx[i]]+=(double)fqr[i];
      }
    }
    mxFree(fqr);
  }
  if(pr!=NULL) {
    if(pi!=NULL) {
      for(int i=0;i<N;i++) {
        pr[ix[i]]+=(double)fpr[i];
        pi[ix[i]]+=(double)fpi[i];
      }
      mxFree(fpi);
    }
    else {
      for(int i=0;i<N;i++) {
        pr[ix[i]]+=(double)fpr[i];
      }
    }
    mxFree(fpr);
  }
#endif //MULTIPOLEEVALCUDA
    
  //convert to double and add to the original vectors
#ifdef CUDADEBUGVECTOR
  double *debugvector0=(double*)mxMalloc(imax(2000,(GPUvars->NE<GPUvars->N?GPUvars->N:GPUvars->NE))*sizeof(double));
  cudasafe( cudaMemcpy(debugvector0, GPUvars->debugvector, imax(2000,(GPUvars->NE<GPUvars->N?GPUvars->N:GPUvars->NE))*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy debugvector" );
  for(int i=0;i<imax(2000,(GPUvars->NE<GPUvars->N?GPUvars->N:GPUvars->NE));i++) {
    if(debugvector0[i]!=0)
      mexPrintf("%d: debugv=%.14e\n"/* zr=%f diff=%f\n"*/,i,debugvector0[i]/*,temvector[i],debugvector0[i]-temvector[i]*/);
  }
  cudaFreeDebug(GPUvars->debugvector);
  mxFree(debugvector0);
#endif //CUDADEBUGVECTOR
            
  //unbind all textures and deallocate all GPU variables. (Note that GPU variables remain in GPU memory until "clear fmm2d" is called, i.e. they are not deallocated when mex file has completed)
    
  cudaUnbindTexture(tzr);
  cudaUnbindTexture(tzi);
  cudaUnbindTexture(tmr);
  if(GPUvars->mi!=NULL){
    cudaUnbindTexture(tmi);
    cudaFreeDebug(GPUvars->mi);
  }
  cudaFreeDebug(GPUvars->zr);
  cudaFreeDebug(GPUvars->zi);
  if(NE) {
    cudaFreeDebug(GPUvars->qr);
    if(qi!=NULL)
      cudaFreeDebug(GPUvars->qi);
  }
  if(pr!=NULL) {
    cudaFreeDebug(GPUvars->pr);
    if(pi!=NULL)
      cudaFreeDebug(GPUvars->pi);
  }
  cudaFreeDebug(GPUvars->mr);
    
  //Print timing information
  TimeSyncPrintAndDestroy(GPUvars->sortstart,GPUvars->sortstop,20,"sort");
  TimeSyncPrintAndDestroy(GPUvars->reorder1start,GPUvars->reorder1stop,21,"GPU reorder");
#ifdef MULTIPOLEINITCUDA
  TimeSyncPrintAndDestroy(GPUvars->initstart,GPUvars->initstop,22,"GPU init");
  TimeSyncPrintAndDestroy(GPUvars->stepm2mstart,GPUvars->stepm2mstop,23,"GPU m2m");
#endif
#ifdef MULTIPOLESHIFTCUDA
  TimeSyncPrintAndDestroy(GPUvars->stepm2psstart,GPUvars->stepm2psstop,24,"GPU m2ps");
  TimeSyncPrintAndDestroy(GPUvars->stepp2pstart,GPUvars->stepp2pstop,25,"GPU p2p");
#endif
#ifdef MULTIPOLEEVALCUDA
  TimeSyncPrintAndDestroy(GPUvars->mpstart,GPUvars->mpstop,26,"GPU mpeval");
  TimeSyncPrintAndDestroy(GPUvars->reorder2start,GPUvars->reorder2stop,28,"GPU reverse reorder");
#endif
  TimeSyncPrintAndDestroy(GPUvars->start,GPUvars->stop,27,"GPU direct");
  TimeEventRecord(GPUvars->fulltimestop);
  TimeSyncPrintAndDestroy(GPUvars->fulltimestart,GPUvars->fulltimestop,29,"GPU total");
#ifdef CUDADEBUGCHECKTIME
  TimeSyncPrintAndDestroy(GPUvars->checkstart,GPUvars->checkstop,30,"GPU check");
#endif
            
  //continue CLEANUP
  //                     checkcudaerror("Check before cleanup");
  cudaUnbindTexture(tixptr);
  cudaFreeDebug(GPUvars->ixptr);
  cudaUnbindTexture(tjxptr);
  if(GPUvars->er!=NULL) {
    cudaFreeDebug(GPUvars->jxptr);
  }
#if defined(MULTIPOLEEVALCUDA) || defined(MULTIPOLEINITCUDA)
  cudaUnbindTexture(tz0);
//   cudaUnbindTexture(tcoeff1);
//   cudaUnbindTexture(tcoeff2);
  cudaUnbindTexture(tcoeff1d);
  cudaUnbindTexture(tcoeff2d);
    
  cudaFreeDebug(GPUvars->z0);
  cudaFreeDebug(GPUvars->coeff1);
  cudaFreeDebug(GPUvars->coeff2);
#endif

#if !defined(MULTIPOLEEVALCUDA) && !defined(MULTIPOLEINITCUDA)
  cudaFreeDebug(GPUvars->z0);
#endif
  cudaUnbindTexture(ter);
  cudaUnbindTexture(tei);
  if(GPUvars->er!=NULL) {
        
    cudaFreeDebug(GPUvars->er);
    cudaFreeDebug(GPUvars->ei);
  }
  if(NE)
    cudaFreeDebug(GPUvars->jx);
  if(pr!=NULL)
    cudaFreeDebug(GPUvars->ix);
  for (int l = 0; l <= This->nlevel; l++) {
    cudaFreeDebug(GPUvars->connect[l].jcptr);
    GPUvars->connect[l].jcptr = NULL;
  }
#ifdef MULTIPOLESHIFTCUDA
  cudaFreeDebug(GPUvars->jcptr2);
  cudaFreeDebug(GPUvars->kcptr2);
  cudaFreeDebug(GPUvars->ir2);
#endif
  mxFree(GPUvars->connect);
  CHECKCUDAERROR
  checkcudaerror("Check after final cuda call");
}
/*------------------------------------------------------------------------*/
//multipole evaluation on gpu
#ifdef MULTIPOLEEVALCUDA
void mpexp_eval_cuda(const MPexp *This,cudavariables *GPUvars,
                     const double *pr,int N,int NE,
                     int pot,double *timing,int printtime)
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
#if !defined(MULTIPOLEINITCUDA) || !defined(MULTIPOLESHIFTCUDA)
  int Nt=(Nf-1)/3+Nf;
#endif
  TimeCreateAndStart(GPUvars->mpstart,GPUvars->mpstop);
  
//special cases if some part of the algorithm is done on the CPU
#ifndef MULTIPOLEINITCUDA
  cudasafe( cudaMemcpy(GPUvars->coeff1, This->root->coeff1, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyHostToDevice), "cudaMemcpy coeff1" );
#endif
#ifndef MULTIPOLESHIFTCUDA //if the correct one is not already on cuda
  cudasafe( cudaMemcpy(GPUvars->coeff2, This->root->coeff2, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyHostToDevice), "cudaMemcpy coeff2" );
#endif
  //call the kernel two times if evaluation points are used as well as source point evaluations
  if(NE) { //the textures should already be left in this position after direct summation
    //         cudaUnbindTexture(ter);
    //         cudaUnbindTexture(tei);
    //         cudaUnbindTexture(tjxptr);
    //         cudasafe(cudaBindTexture(NULL,ter,GPUvars->er,NE*sizeof(double)),"cudaBindTexture ter (zr)");
    //         cudasafe(cudaBindTexture(NULL,tei,GPUvars->ei,NE*sizeof(double)),"cudaBindTexture ter (zr)");
    //         cudasafe(cudaBindTexture(NULL,tjxptr,GPUvars->jxptr,(Nf+1)*sizeof(double)),"cudaBindTexture ter (zr)");
    cuda_mpexp_eval<<<imin(MAXBLOCKS, Nf), MAXTHREADSEVALUATION>>>(This->pcoeff,GPUvars->qr,GPUvars->qi,Nf,This->lptr[This->nlevel],&GPUvars->connect[This->nlevel].kcptr[Nf],GPUvars->connect[This->nlevel].kcptr,GPUvars->connect[This->nlevel].ir,pot DEBUGVECTORSTRING2);
    CHECKCUDAERROR
  }
  if(pr!=NULL) {
    //move the textures to the source points and make the same call to cuda
    cudaUnbindTexture(ter);
    cudaUnbindTexture(tei);
    cudaUnbindTexture(tjxptr);
    cudasafe(cudaBindTexture(NULL,ter,GPUvars->zr,N*sizeof(double)),"cudaBindTexture ter (zr)");
    cudasafe(cudaBindTexture(NULL,tei,GPUvars->zi,N*sizeof(double)),"cudaBindTexture ter (zr)");
    cudasafe(cudaBindTexture(NULL,tjxptr,GPUvars->ixptr,(Nf+1)*sizeof(double)),"cudaBindTexture ter (zr)");
    cuda_mpexp_eval<<<imin(MAXBLOCKS, Nf), MAXTHREADSEVALUATION>>>(This->pcoeff,GPUvars->pr,GPUvars->pi,Nf,This->lptr[This->nlevel],&GPUvars->connect[This->nlevel].kcptr[Nf],GPUvars->connect[This->nlevel].kcptr,GPUvars->connect[This->nlevel].ir,pot DEBUGVECTORSTRING2);
    CHECKCUDAERROR
  }
  TimeEventRecord(GPUvars->mpstop);
  checkcudaerror("Kernel launch error mpeval");
}
#endif /* MULTIPOLEEVALCUDA */
/*------------------------------------------------------------------------*/
//m2m evaluation. Simply make the CUDA call
#ifdef MULTIPOLEINITCUDA
//two possiblies for m2m shift, first one without pre/post scaling
void MPexp_shiftm2m_cuda_horner(const MPexp *This,cudavariables *GPUvars,
                                int pot,double* timing,int printtime)
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  TimeCreateAndStart(GPUvars->stepm2mstart,GPUvars->stepm2mstop);
  for (int l = This->nlevel, nb = Nf>>2; l > 0; l--, nb >>= 2) { //loop over all levels
    int dynalloc=(SHIFT_M2M_MAXTHREADS*(This->pcoeff+1)*sizeof(dcmplx)); //calculate how much shared memory that has to be allocated
    shift_m2m_cuda_horner<<<imin((nb*4+SHIFT_M2M_MAXTHREADS-1)/SHIFT_M2M_MAXTHREADS,MAXBLOCKS),SHIFT_M2M_MAXTHREADS,dynalloc>>>(l-1,This->lptr[l-1],This->pcoeff,
                                                                                                                                ((dcmplx *)GPUvars->coeff1),(SORT_DCMPLX*)GPUvars->z0,pot  DEBUGVECTORSTRING2);
    CHECKCUDAERROR
  }
  checkcudaerror("Kernel launch error shift_m2m");
  TimeEventRecord(GPUvars->stepm2mstop);
}
/*------------------------------------------------------------------------*/
//Calls the scaled version instead. Is usually faster
void MPexp_shiftm2m_cuda_horner_scaled(const MPexp *This,
                                       cudavariables *GPUvars,
                                       int pot,double *timing,int printtime)
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  TimeCreateAndStart(GPUvars->stepm2mstart,GPUvars->stepm2mstop);
  double rminshift = exp2(-1000.0/(This->pcoeff+1));
  for (int l = This->nlevel, nb = Nf>>2; l > 0; l--, nb >>= 2) { //loop over all levels
    int dynalloc=(SHIFT_M2M_MINTHREADS*(This->pcoeff+1)*sizeof(dcmplx))/2; //calculate how much shared memory that has to be allocated
    if((dynalloc+SHIFT_M2M_MAXTHREADS/4*sizeof(dcmplx))*CUDABLOCKSPERMP<GPUvars->info.sharedMemPerBlock) {
      dynalloc=(SHIFT_M2M_MAXTHREADS*(This->pcoeff+1)*sizeof(dcmplx))/2; //calculate how much shared memory that has to be allocated
      shift_m2m_cuda_horner_scaled_t<<<imin((nb*8+SHIFT_M2M_MAXTHREADS-1)/SHIFT_M2M_MAXTHREADS,MAXBLOCKS),SHIFT_M2M_MAXTHREADS,dynalloc>>>(l-1,This->lptr[l-1],This->pcoeff,
                                                                                                                                          ((dcmplx *)GPUvars->coeff1),(SORT_DCMPLX*)GPUvars->z0,pot,rminshift  DEBUGVECTORSTRING2);
      CHECKCUDAERROR
    }
    else {
      shift_m2m_cuda_horner_scaled_t<<<imin((nb*8+SHIFT_M2M_MINTHREADS-1)/SHIFT_M2M_MINTHREADS,MAXBLOCKS),SHIFT_M2M_MINTHREADS,dynalloc>>>(l-1,This->lptr[l-1],This->pcoeff,
                                                                                                                                          ((dcmplx *)GPUvars->coeff1),(SORT_DCMPLX*)GPUvars->z0,pot,rminshift  DEBUGVECTORSTRING2);
      CHECKCUDAERROR
    }
  }
  checkcudaerror("Kernel launch error shift_m2m");
  TimeEventRecord(GPUvars->stepm2mstop);
}
/*------------------------------------------------------------------------*/
//write initial coefficients on gpu. This code is similar to the multipole evaluation
void MPexp_init_cuda(const MPexp *This,cudavariables *GPUvars,
                     int pot,int complexpoint,double *timing,int printtime)
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  checkcudaerror("Check Init start");
  TimeCreateAndStart(GPUvars->initstart,GPUvars->initstop);
  //launch the GPU kernel
  checkcudaerror("Check before Init");
  cuda_mpexp_init<<<imin(MAXBLOCKS, GPUvars->Nf), initthreadperblock,(This->pcoeff+1)*sizeof(dcmplx)>>>(This->pcoeff,GPUvars->coeff1,GPUvars->coeff2,Nf,This->lptr[This->nlevel],&GPUvars->connect[This->nlevel].kcptr[Nf],GPUvars->connect[This->nlevel].kcptr,GPUvars->connect[This->nlevel].ir,complexpoint,pot DEBUGVECTORSTRING2);
  CHECKCUDAERROR
  checkcudaerror("Kernel launch error Init");
  TimeEventRecord(GPUvars->initstop);
            
  //Continue by calling m2m call shiftm2m
  MPexp_shiftm2m_cuda_horner_scaled(This,GPUvars,pot,timing,printtime);
#if !defined( MULTIPOLESHIFTCUDA) || defined(CHECKP2P) || defined(CHECKM2PS)
  int Nt=(Nf-1)/3+Nf;
  cudasafe( cudaMemcpy(This->root->coeff1,GPUvars->coeff1, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff1" );
  cudasafe( cudaMemcpy(This->root->coeff2,GPUvars->coeff2, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff1" );
#endif
    
// #ifndef MULTIPOLEEVALCUDA //if multipole evaluation, keep the vectors on the gpu, note however that distbox and boxstart has to be reordered
//   cudaUnbindTexture(tz0);
//   cudaUnbindTexture(tcoeff1d);
//   cudaUnbindTexture(tcoeff2d);
//   cudaFreeDebug(GPUvars->z0);
//   cudaFreeDebug(GPUvars->coeff1);
//   cudaFreeDebug(GPUvars->coeff2);
//   GPUvars->z0=NULL;
//   GPUvars->coeff1=NULL;
//   GPUvars->coeff2=NULL;
// #endif
}
/*------------------------------------------------------------------------*/
//validates the coefficient initialization against the values calculated by the CPU
//this is debug code, activated by defining CHECKINITCUDA in fmm.h
void cuda_check_init(MPexp *This,int Nt,cudavariables* GPUvars)
{
  mexPrintf("check init cuda, pcoeff=%d Nt=%d\n",This->pcoeff,Nt);
  dcmplx *backup2=(dcmplx*)mxMalloc((This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudaSafeMemcpy(backup2, ((cudavariables*)GPUvars)->coeff1, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost, "cudaMemcpy coeff1" );
  for(int i=0;i<(This->pcoeff+1)*Nt;i++) {
    if((hypot(creal(This->root->coeff1[i]),cimag(This->root->coeff1[i]))!=0&&(fabs(hypot(creal(This->root->coeff1[i]),cimag(This->root->coeff1[i]))-hypot(creal(backup2[i]),cimag(backup2[i])))/hypot(creal(This->root->coeff1[i]),cimag(This->root->coeff1[i]))>1e-12||atan2(creal(backup2[i]),cimag(backup2[i]))-atan2(creal(This->root->coeff1[i]),cimag(This->root->coeff1[i]))>1e-8))||
       hypot(creal(This->root->coeff1[i]),cimag(This->root->coeff1[i]))==0&&hypot(creal(backup2[i]),cimag(backup2[i]))!=0||
       !isfinite(hypot(creal(backup2[i]),cimag(backup2[i])))||
       !isfinite(hypot(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))))
      mexPrintf("Difference: Box[%d] coeff1[%d]=%.14e+%.14ei, GPU=%.14e+%.14ei\n",i/(This->pcoeff+1),i%(This->pcoeff+1),creal(This->root->coeff1[i]),cimag(This->root->coeff1[i]),creal(backup2[i]),cimag(backup2[i]));
  }
  memcpy(This->root->coeff1,backup2,(This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudaSafeMemcpy(backup2, ((cudavariables*)GPUvars)->coeff2, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost, "cudaMemcpy coeff1" );
  for(int i=0;i<(This->pcoeff+1)*Nt;i++) {
    if((hypot(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))!=0&&(fabs(hypot(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))-hypot(creal(backup2[i]),cimag(backup2[i])))/hypot(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))>1e-12||atan2(creal(backup2[i]),cimag(backup2[i]))-atan2(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))>1e-8))||
       hypot(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))==0&&hypot(creal(backup2[i]),cimag(backup2[i]))!=0||
       !isfinite(hypot(creal(backup2[i]),cimag(backup2[i])))||
       !isfinite(hypot(creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]))))
      mexPrintf("Difference: Box[%d] coeff2[%d]=%.14e+%.14ei, GPU=%.14e+%.14ei\n",i/(This->pcoeff+1),i%(This->pcoeff+1),creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]),creal(backup2[i]),cimag(backup2[i]));
  }
  memcpy(This->root->coeff2,backup2,(This->pcoeff+1)*Nt*sizeof(dcmplx));
  mxFree(backup2);
}
#endif /* MULTIPOLEINITCUDA */
/*------------------------------------------------------------------------*/
//shift p2p, traditional horner scheme.
#ifdef MULTIPOLESHIFTCUDA
void MPexp_shiftp2p_cuda(const MPexp *This,cudavariables *GPUvars,
                         int levm2p,double* timing,int printtime)
{
  TimeCreateAndStart(GPUvars->stepp2pstart,GPUvars->stepp2pstop);
    
    
  for (int l = levm2p, nb = 1 << 2*(levm2p+1); l < This->nlevel; l++, nb <<= 2) {
    int p2pthreadcount=imin((nb+P2PMAXTHREADS-1)/P2PMAXTHREADS, MAXBLOCKS);
    shift_p2p_cuda<<<p2pthreadcount, P2PMAXTHREADS, P2PMAXTHREADS*(This->pcoeff+1)*sizeof(dcmplx)>>>(l, This->lptr[l], This->pcoeff,
                                                                                                     ((dcmplx *)GPUvars->coeff2)  DEBUGVECTORSTRING2);
    CHECKCUDAERROR
  }
  TimeEventRecord(GPUvars->stepp2pstop);
            
  //if error checkin is active, activated from fmm.h, but should be turned of for production codes
#ifdef CHECKP2P
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  int Nt=(Nf-1)/3+Nf;
  mexPrintf("checking p2p interaction\n");
  dcmplx *backup=(dcmplx*)mxMalloc((This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudasafe( cudaMemcpy(backup, GPUvars->coeff2, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff2" );
  for(int i=0;i<(This->pcoeff+1)*Nt;i++) {
    if((creal(backup[i])!=0&&fabs((creal((This->root->coeff2[i]))-creal(backup[i]))/creal(backup[i]))>1e-9)||(cimag(backup[i])!=0&&fabs((cimag((This->root->coeff2[i]))-cimag(backup[i]))/cimag(backup[i]))>1e-9))
      mexPrintf("Difference: Box[%d] coeff[%d]=%.14e+%.14ei, GPU=%.14e+%.14ei\n",i/(This->pcoeff+1),i%(This->pcoeff+1),creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]),creal(backup[i]),cimag(backup[i]));
  }
  mxFree(backup);
#endif
  checkcudaerror("Kernel launch error p2p");
#if !defined(MULTIPOLEEVALCUDA)
  const int Nf2 = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  const int Nt2=(Nf2-1)/3+Nf2;
  cudasafe( cudaMemcpy(This->root->coeff1,GPUvars->coeff1, (This->pcoeff+1)*Nt2*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff1" );
  cudasafe( cudaMemcpy(This->root->coeff2,GPUvars->coeff2, (This->pcoeff+1)*Nt2*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff1" );
#endif
    
}
/*------------------------------------------------------------------------*/
//the scaled version, is usually faster
void MPexp_shiftp2p_cuda_scaled(const MPexp *This,cudavariables *GPUvars,
                                int levm2p,double *timing,int printtime)
{
  TimeCreateAndStart(GPUvars->stepp2pstart,GPUvars->stepp2pstop);
    
  double rminshift = exp2(-1000.0/(This->pcoeff+1));
  for (int l = levm2p, nb = 1 << 2*(levm2p+1); l < This->nlevel; l++, nb <<= 2) {
    int p2pthreadcount=imin((nb*2+P2PMAXTHREADS-1)/P2PMAXTHREADS, MAXBLOCKS);
    int dynalloc=P2PMINTHREADS/2*(This->pcoeff+1)*sizeof(dcmplx);
    if((dynalloc+P2PMAXTHREADS/4*sizeof(dcmplx))*CUDABLOCKSPERMP<GPUvars->info.sharedMemPerBlock) {
      dynalloc=P2PMAXTHREADS/2*(This->pcoeff+1)*sizeof(dcmplx);
      int p2pthreadcount=imin((nb*2+P2PMAXTHREADS-1)/P2PMAXTHREADS, MAXBLOCKS);
      shift_p2p_cuda_horner_scaled_t<<<p2pthreadcount, P2PMAXTHREADS, dynalloc>>>(l, This->lptr[l], This->pcoeff,
                                                                                 ((dcmplx *)GPUvars->coeff2),((SORT_DCMPLX *)GPUvars->z0),rminshift  DEBUGVECTORSTRING2);
      CHECKCUDAERROR
    }
    else {
//       dynalloc=P2PMINTHREADS/2*(This->pcoeff+1)*sizeof(dcmplx);
      int p2pthreadcount=imin((nb*2+P2PMINTHREADS-1)/P2PMINTHREADS, MAXBLOCKS);
      shift_p2p_cuda_horner_scaled_t<<<p2pthreadcount, P2PMINTHREADS, dynalloc>>>(l, This->lptr[l], This->pcoeff,
                                                                                 ((dcmplx *)GPUvars->coeff2),((SORT_DCMPLX *)GPUvars->z0),rminshift  DEBUGVECTORSTRING2);
      CHECKCUDAERROR
    }
  }
  TimeEventRecord(GPUvars->stepp2pstop);
            
  //if error checkin is active, activated from fmm.h, but should be turned of for production codes
#ifdef CHECKP2P
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  int Nt=(Nf-1)/3+Nf;
  mexPrintf("checking p2p interaction p = %d\n",This->pcoeff);
  dcmplx *backup=(dcmplx*)mxMalloc((This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudasafe( cudaMemcpy(backup, GPUvars->coeff2, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff2" );
  for(int i=0;i<(This->pcoeff+1)*Nt;i++) {
    if((creal(backup[i])!=0&&fabs((creal((This->root->coeff2[i]))-creal(backup[i]))/creal(backup[i]))>1e-9)||(cimag(backup[i])!=0&&fabs((cimag((This->root->coeff2[i]))-cimag(backup[i]))/cimag(backup[i]))>1e-9))
      mexPrintf("Difference: Box[%d] coeff[%d]=%.14e+%.14ei, GPU=%.14e+%.14ei\n",i/(This->pcoeff+1),i%(This->pcoeff+1),creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]),creal(backup[i]),cimag(backup[i]));
  }
  mxFree(backup);
#endif
  checkcudaerror("Kernel launch error p2p");
#if !defined(MULTIPOLEEVALCUDA)
  const int Nf2 = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  const int Nt2=(Nf2-1)/3+Nf2;
  cudasafe( cudaMemcpy(This->root->coeff1,GPUvars->coeff1, (This->pcoeff+1)*Nt2*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff1" );
  cudasafe( cudaMemcpy(This->root->coeff2,GPUvars->coeff2, (This->pcoeff+1)*Nt2*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff1" );
#endif
}
/*------------------------------------------------------------------------*/
//The m2ps interaction. Since there is no connection between the levels in
//this interaction, one single cuda call is enough for all levels
//only scaled version exists for this shift operator
void MPexp_shiftm2ps_cuda(const MPexp *This,cudavariables *GPUvars,
                          int levm2p,int pot,double* timing,int printtime)
{
  int count=This->lptr[This->nlevel+1];
  TimeCreateAndStart(GPUvars->stepm2psstart,GPUvars->stepm2psstop);
  if(pot==0) {
    shift_m2ps_cuda_assym_log<<<imin(count, MAXBLOCKS), M2PSMAXTHREADS, (This->pcoeff+1)*(M2PSMAXTHREADS+2)*sizeof(double)>>>(1,This->nlevel, This->pcoeff,
                                                                                                                             (dcmplx*)GPUvars->coeff1,(dcmplx*)GPUvars->coeff2, (SORT_DCMPLX*)GPUvars->z0, GPUvars->jcptr2,GPUvars->kcptr2,GPUvars->ir2  DEBUGVECTORSTRING2);
    CHECKCUDAERROR
  }
  else {
    shift_m2ps_cuda_assym<<<imin(count, MAXBLOCKS), M2PSMAXTHREADS, (This->pcoeff+1)*(M2PSMAXTHREADS+2)*sizeof(double)>>>(1,This->nlevel, This->pcoeff,
                                                                                                                         (dcmplx*)GPUvars->coeff1,(dcmplx*)GPUvars->coeff2, (SORT_DCMPLX*)GPUvars->z0, GPUvars->jcptr2,GPUvars->kcptr2,GPUvars->ir2  DEBUGVECTORSTRING2);
    CHECKCUDAERROR
  }
  TimeEventRecord(GPUvars->stepm2psstop);
  
  //if error checkin is active, activated from fmm.h, but should be turned of for production codes
#ifdef CHECKM2PS
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  int Nt=(Nf-1)/3+Nf;
  mexPrintf("checking m2ps interaction p = %d\n",This->pcoeff);
  dcmplx *backup=(dcmplx*)mxMalloc((This->pcoeff+1)*Nt*sizeof(dcmplx));
  cudasafe( cudaMemcpy(backup, GPUvars->coeff2, (This->pcoeff+1)*Nt*sizeof(dcmplx), cudaMemcpyDeviceToHost), "cudaMemcpy coeff2" );
  for(int i=0;i<(This->pcoeff+1)*Nt;i++) {
    if((creal(backup[i])!=0&&fabs((creal((This->root->coeff2[i]))-creal(backup[i]))/creal(backup[i]))>1e-9)||(cimag(backup[i])!=0&&fabs((cimag((This->root->coeff2[i]))-cimag(backup[i]))/cimag(backup[i]))>1e-9)||creal(backup[i])==0&&creal(This->root->coeff2[i])!=0||cimag(backup[i])==0&&cimag(This->root->coeff2[i])!=0)
      mexPrintf("Difference: Box[%d] coeff2[%d]=%.14e+%.14ei, GPU=%.14e+%.14ei\n",i/(This->pcoeff+1),i%(This->pcoeff+1),creal(This->root->coeff2[i]),cimag(This->root->coeff2[i]),creal(backup[i]),cimag(backup[i]));
  }
  mxFree(backup);
#endif
}
#endif /*MULTIPOLESHIFTCUDA
/*------------------------------------------------------------------------*/
//prototype
void directsumkernellaunch(const double *mi,double *GPUqr,double *GPUqi,
                           int N,int NE,int pot,SMOOTHER smooth,double cutoff,
                           double shape,double scale);

//direct evaluation (tol==0). This function is independent of the others and
//performs its own memory allocation, cuda calls cleanup etc.
void direct_eval_cuda(int N,
                      const double *zr,const double *zi,
                      const double *mr,const double *mi,
                      int NE,
                      const double *er,const double *ei,
                      double *pr,double *pi,
                      double *qr,double *qi,
                      const panel *panels,int panelcount,int pot,
                      SMOOTHER smooth,double xopt,double cutoff,
                      bool cont,double* timing,int printtime)
//function to start the GPU direct summation calculations
{
  cudavariables ctmp;
  cudavariables *GPUvars=&ctmp; //simply to always use GPUvars as a pointer
  double shape, scale;
  if(N==0) /*no potential points, code will give Unspecified Driver Error otherwise*/
    return;
  cudastart();
  checkcudaerror("Check eval_cuda start");
  TimeCreateAndStart(GPUvars->fulltimestart,GPUvars->fulltimestop);
            
  //start by making necessary allocations and copy memory to GPU
  if(NE) {
    cudasafeMalloc((void**)&GPUvars->er,NE*sizeof(double));
    cudasafeMalloc((void**)&GPUvars->ei,NE*sizeof(double));
    cudasafeMalloc((void**)&GPUvars->qr,NE*sizeof(double));
    cudasafe( cudaMemcpy(GPUvars->er, er, NE*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy er" );
    cudasafe( cudaMemcpy(GPUvars->ei, ei, NE*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy ei" );
    if(pot!=0||mi!=NULL)
      cudasafeMalloc((void**)&GPUvars->qi,NE*sizeof(double));
  }
  else {
    GPUvars->er=NULL;
    GPUvars->ei=NULL;
  }
  if(pr!=NULL) {
    cudasafeMalloc((void**)&GPUvars->pr,N*sizeof(double));
    if(pot!=0||mi!=NULL) {
      cudasafeMalloc((void**)&GPUvars->pi,N*sizeof(double));
    }
  }
  cudasafeMalloc((void**)&GPUvars->zr,N*sizeof(double));
  cudasafeMalloc((void**)&GPUvars->zi,N*sizeof(double));
  cudasafeMalloc((void**)&GPUvars->mr,N*sizeof(double));
    
  if(mi!=NULL) {
    cudasafeMalloc((void**)&GPUvars->mi,N*sizeof(double));
    cudasafe(cudaBindTexture(NULL,tmi,GPUvars->mi,N*sizeof(double)),"cudaBindTexture tmi");
    cudasafe( cudaMemcpy(GPUvars->mi, mi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mi" );
  }
  else
    GPUvars->mi=NULL;
  MPexp_smooth_(pot, smooth, xopt, cont, &cutoff, &shape, &scale);

  cudasafe( cudaMemcpy(GPUvars->zr, zr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zr" );
  cudasafe( cudaMemcpy(GPUvars->zi, zi, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy zi" );
  cudasafe( cudaMemcpy(GPUvars->mr, mr, N*sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy mr" );

  cudasafe(cudaBindTexture(NULL,tzr,GPUvars->zr,N*sizeof(double)),"cudaBindTexture tzr");
  cudasafe(cudaBindTexture(NULL,tzi,GPUvars->zi,N*sizeof(double)),"cudaBindTexture tzi");
  cudasafe(cudaBindTexture(NULL,tmr,GPUvars->mr,N*sizeof(double)),"cudaBindTexture tmr");

  if(timing!=NULL||printtime) {
    cudaTimingCreateAndStart(&GPUvars->start,&GPUvars->stop,"start","stop");
  }
  //move textures to make it possible to use same code in both cases.
  if(pr!=NULL) {
    cudasafe(cudaBindTexture(NULL,ter,GPUvars->zr,N*sizeof(double)),"cudaBindTexture ter");
    cudasafe(cudaBindTexture(NULL,tei,GPUvars->zi,N*sizeof(double)),"cudaBindTexture tei");
    
    //launch GPU kernel
    directsumkernellaunch(mi,GPUvars->pr,GPUvars->pi,N,N,pot,smooth,cutoff,shape,scale);
  }
  if(NE) {
    if(pr!=NULL) {
      cudaUnbindTexture(ter);
      cudaUnbindTexture(tei);
    }
    cudasafe(cudaBindTexture(NULL,ter,GPUvars->er,NE*sizeof(double)),"cudaBindTexture ter");
    cudasafe(cudaBindTexture(NULL,tei,GPUvars->ei,NE*sizeof(double)),"cudaBindTexture tei");
    
    //launch GPU kernel
    directsumkernellaunch(mi,GPUvars->qr,GPUvars->qi,N,NE,pot,smooth,cutoff,shape,scale);
  }
  if(timing!=NULL||printtime) {
    cudasafe(cudaEventRecord(GPUvars->stop), "cudaEventRecord");
  }

  //start the cleanup    
  if(NE) {
    cudasafe( cudaMemcpy(qr,GPUvars->qr, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qr" );
    cudaFreeDebug(GPUvars->er);
    cudaFreeDebug(GPUvars->ei);
    cudaFreeDebug(GPUvars->qr);
    if(pot!=0||mi!=NULL) {
      cudasafe( cudaMemcpy(qi,GPUvars->qi, NE*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy qi" );
      cudaFreeDebug(GPUvars->qi);
    }
  }
  if(pr!=NULL) {
    cudasafe( cudaMemcpy(pr,GPUvars->pr, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy pr" );
    cudaFreeDebug(GPUvars->pr);
    if(pot!=0||mi!=NULL) {
      cudasafe( cudaMemcpy(pi,GPUvars->pi, N*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy pr" );
      cudaFreeDebug(GPUvars->pi);
    }
  }
  cudaUnbindTexture(tzr);
  cudaUnbindTexture(tzi);
  cudaUnbindTexture(tmr);
  cudaUnbindTexture(ter);
  cudaUnbindTexture(tei);
    
  cudaFreeDebug(GPUvars->zr);
  cudaFreeDebug(GPUvars->zi);
  cudaFreeDebug(GPUvars->mr);
    
  if(mi!=NULL) {
    cudaUnbindTexture(tmi);
    cudaFreeDebug(GPUvars->mi);
  }  
  TimeEventRecord(GPUvars->fulltimestop);
  //do not synchronize until last possible moment
  TimeSyncPrintAndDestroy(GPUvars->start, GPUvars->stop,0,"GPU direct sum");
  TimeSyncPrintAndDestroy(GPUvars->fulltimestart, GPUvars->fulltimestop,1,"GPU total");
}
/*------------------------------------------------------------------------*/
//selects correct direct summation kernel and starts it
void directsumkernellaunch(const double *mi,double *GPUqr,double *GPUqi,
                           int N,int NE,int pot,SMOOTHER smooth,double cutoff,
                           double shape,double scale)
{
  cudaDeviceProp info;
  int device;
  cudasafe(cudaGetDevice(&device),"cudaGetDevice");
  cudasafe(cudaGetDeviceProperties(&info,device),"cudaGetDeviceProperties");
  //load balancing. Each multiprocessor can run up to 8 blocks (current models, this information is sadly not available in the info structure)
  //mexPrintf("DIRECTMAXTHREADS: %d NE: %d multiprocessorcount: %d CUDABLOCKSPERMP: %d warpsize: %d", DIRECTMAXTHREADS,NE,info.multiProcessorCount, CUDABLOCKSPERMP,info.warpSize);
  int nthreads=imin(DIRECTMAXTHREADS,ceil((double)NE/(double)info.multiProcessorCount/(double)CUDABLOCKSPERMP/(double)info.warpSize)*info.warpSize);
  if(nthreads<0)
    nthreads=DIRECTMAXTHREADS;
  if(pot==1) {
    if(mi==NULL) {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        direct_pot1_nomi_dirac_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr, GPUqi,N,NE);
        break;
      case RANKINE:
        direct_pot1_nomi_rankine_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape);
        break;
      case SCULLY:
        direct_pot1_nomi_scully_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, scale, cutoff);
        break;
      case OSEEN:
        direct_pot1_nomi_oseen_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, scale, cutoff);
        break;
      }
    }
    else {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        direct_pot1_mi_dirac_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE);
        break;
      case RANKINE:

        direct_pot1_mi_rankine_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape);
        break;
      case SCULLY:
        direct_pot1_mi_scully_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, scale, cutoff);
        break;
      case OSEEN:
        direct_pot1_mi_oseen_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, scale, cutoff);
        break;
      }
    }
  }
  else {
    if(mi==NULL) {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        direct_pot0_nomi_dirac_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr,N,NE);
        break;
      case RANKINE:
        direct_pot0_nomi_rankine_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr,N,NE, shape, cutoff);
        break;
      case SCULLY:
        direct_pot0_nomi_scully_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr,N,NE, shape, scale, cutoff);
        break;
      case OSEEN:
        direct_pot0_nomi_oseen_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*3*sizeof(double)>>>(GPUqr,N,NE, shape, scale, cutoff);
        break;
      }
    }
    else {
      switch(smooth) { //one kernel for each type. This is to reduce the number of local variables for each kernel
      case DIRAC:
        direct_pot0_mi_dirac_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE);
        break;
      case RANKINE:
        direct_pot0_mi_rankine_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, cutoff);
        break;
      case SCULLY:
        direct_pot0_mi_scully_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, scale, cutoff);
        break;
      case OSEEN:
        direct_pot0_mi_oseen_synchronized<<<imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double)>>>(GPUqr, GPUqi,N,NE, shape, scale, cutoff);
        break;
      }
    }
  }
  //mexPrintf("blocks: %d threads: %d sharedmem: %d\n", imin(MAXBLOCKS, (NE+nthreads-1)/nthreads), nthreads,nthreads*4*sizeof(double));
  checkcudaerror("Kernel launch error direct interact directsumkernellaunch");
  CHECKCUDAERROR
}
/*------------------------------------------------------------------------*/
//function to cleanup CUDASORT/MPexp_connect_CUDA_, used if only sorting is done, like out=mesh
void cleanupsort(cudavariables* GPUvars,const MPexp *This)
{
#ifdef CUDADEBUGVECTOR
  cudaFreeDebug(GPUvars->debugvector);
#endif
    
  if(GPUvars->er!=NULL)
    cudaFreeDebug(GPUvars->er);
  if(GPUvars->ei!=NULL)
    cudaFreeDebug(GPUvars->ei);
  if(GPUvars->jxptr!=NULL)
    cudaFreeDebug(GPUvars->jxptr);
  if(GPUvars->zr!=NULL)
    cudaFreeDebug(GPUvars->zr);
  if(GPUvars->zi!=NULL)
    cudaFreeDebug(GPUvars->zi);
  if(GPUvars->ixptr!=NULL)
    cudaFreeDebug(GPUvars->ixptr);
  if(GPUvars->d0!=NULL)
    cudaFreeDebug(GPUvars->d0);
  if(GPUvars->z0!=NULL)
    cudaFreeDebug(GPUvars->z0);
  if(GPUvars->ix!=NULL)
    cudaFreeDebug(GPUvars->ix);
  if(GPUvars->jx!=NULL)
    cudaFreeDebug(GPUvars->jx);
  for (int l = 0; l <= This->nlevel; l++)
    if(GPUvars->connect[l].jcptr!=NULL)
      cudaFreeDebug(GPUvars->connect[l].jcptr);
}
/*------------------------------------------------------------------------*/
void cudaSafeMemcpy(void *dst,const void *src,size_t count,
                    enum cudaMemcpyKind kind,const char *str)
// wrapper
{
  cudasafe(cudaMemcpy(dst,src,count,kind),str);
}
/*------------------------------------------------------------------------*/
void cudaSafeMemset(void *dst,int value,size_t count,const char *str)
// wrapper
{
  cudasafe(cudaMemset(dst,value,count),str);
}
/*------------------------------------------------------------------------*/
cudaError_t cudaGetLastErrorDummy()
// wrapper
{
  return cudaGetLastError();
}
/*------------------------------------------------------------------------*/
void getGPUinfo(cudavariables* GPUvars)
{
  int device;
  cudasafe(cudaGetDevice(&device),"cudaGetDevice");
  cudasafe(cudaGetDeviceProperties(&GPUvars->info,device),"cudaGetDeviceProperties");
}
/*------------------------------------------------------------------------*/
#ifdef CHECKCUDACALLS
void checkcudaerrorwrapper(int line)
{
  char number[51];
  snprintf(number,50,"Line %d",line);
  cudasafe(cudaThreadSynchronize(),number);
}
#endif
#endif /* CUDASUPPORT */
