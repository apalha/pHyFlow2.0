/* fmmsort.cpp */

/* S. Engblom and A. Goude 2011-10-21 */

//#define VALIDATECUDACONNECTIVITY
#ifndef CUDASORT
#undef VALIDATECUDACONNECTIVITY
#endif

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#include <stdlib.h>
#include "fmmsort.h"
#include "cudasort.h"
#include "cudatiming.h"
#include "channelpot.h"
#include "channelpotpanel.h"

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

/*------------------------------------------------------------------------*/
void MPexp_sort_(MPexp *This,
                 int N,const double *zr,const double *zi,
                 int NE,const double *er,const double *ei,
                 const panel *panels,int Npanel,int **panelptrlist,
                 int *maxm2p,int *levm2p,int Ndirect,double cutoff,
                 cudavariables *GPUvars,channelparam* cparams,
                 double* timing,bool printtime,bool meshonly)
/* Sorts the N potentials (zr,zi) in the multipole tree This which is
   simultaneously constructed and, for NE > 0, also assigns each
   evaluation point (er,ei) a multipole box. After each level has been
   thus constructed the connectivity information is also
   determined. The integer maxm2p is set to the maximum number of
   cluster to cluster interactions in the whole tree (the maximum
   number of -1's in any column of all the connectivity matrices) and
   levm2p is the first level where a cluster to cluster interaction
   occurs. This function is critical to performance. */
{
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];
  int *ix,*ixptr,*ixptr0;
  int *jx,*jxptr,*jxptr0;
  mpSparse *C = This->connect;

#ifdef CUDASUPPORT
  cudaGetLastErrorDummy();
#endif
#ifdef CUDASORT
  CUDA_copy_vectors(GPUvars,N,NE,zr,zi,er,ei); //copy data before timing starts
#endif
#ifdef CUDATIME
  cudaEvent_t sort1start,sort1stop,sort2start,sort2stop;
  TimeCreateAndStart(sort1start,sort1stop);
#endif

  if(!fullcudafmm||!cudasort||meshonly) {
    // (ixptr,ix) to be built
    ix = This->ix;
    ixptr = This->ixptr;
    ixptr0 = &ixptr[Nf-1];
    ixptr0[0] = 0;
    ixptr0[1] = N;

    // ditto (jxptr,jx)
    if (NE) {
      jx = This->jx;
      jxptr = This->jxptr;
      jxptr0 = &jxptr[Nf-1];
      jxptr0[0] = 0;
      jxptr0[1] = NE;
    }
    else
      jx = jxptr = NULL;
  }

  /* connectivity: each matrix is strictly lower triangular since a
     box is always in near connection with itself and since
     connectivity is a symmetric relation */
  C[0].jcptr = (int *)mxMalloc(5*sizeof(int));
  C[0].jcptr[0] = 0;
  // convenient in fmm2dGetMesh():
  C[0].kcptr = &C[0].jcptr[2];
  C[0].ir = &C[0].kcptr[2];
  C[0].ir[0] = 0;

  *maxm2p = 0;
  *levm2p = This->nlevel+1;

  // partitioning on the GPU/CPU, respectively
#ifdef CUDASORT
  cudaGetLastErrorDummy();
  // found in cudasort.h:
  CUDA_perform_partitioning(GPUvars,N,NE,This->nlevel VALIDATEPARTITIONINGSTRING2);
  // postprocessing on the CPU
  if(!fullcudafmm||meshonly)
    MPexp_partition_post_CUDA_(This,GPUvars,ix,ixptr,jx,jxptr,N,NE,
                               timing,printtime);
#else
  MPexp_partition_CPU_(This,N,zr,zi,NE,er,ei,
                       panels,Npanel,panelptrlist,
                       Ndirect,cutoff,ix,ixptr,ixptr0,
                       jx,jxptr,jxptr0,cparams);
#endif

#ifdef CUDATIME
  TimeEventRecord(sort1stop);
  TimeSyncPrintAndDestroy(sort1start,sort1stop,1,"GPU sort part 1");
  TimeCreateAndStart(sort2start,sort2stop);
#endif

  // connectivity on the GPU/CPU, respectively
#ifdef CUDASUPPORT
  C[0].jcptr[1] = 1; // direct interaction to self is included
  C[0].kcptr[0] = 1;
  C[0].kcptr[1] = 1;

#ifndef CUDASORT
  // copy z0 and d0 to GPU for connectivity build
  int Nt=(Nf-1)/3+Nf;
  cudasafe(cudaMallocDebug((void **)&GPUvars->dabs,
                           Nt*sizeof(double)),
           "cudaMalloc dabs");
  dcmplx *z0=(dcmplx*)mxMalloc(Nt*sizeof(dcmplx));
  dcmplx *d0=(dcmplx*)mxMalloc(Nt*sizeof(dcmplx));
#ifdef RADIALSHRINK
  double *dabs=(double*)mxMalloc(Nt*sizeof(double));
#endif
  for (int j = 0; j < Nt; j++) {
    COMPLEXASSIGN(z0[j],creal(This->root[j].z0),cimag(This->root[j].z0));
    COMPLEXASSIGN(d0[j],creal(This->root[j].d0),cimag(This->root[j].d0));
#ifdef RADIALSHRINK
    dabs[j]=This->root[j].absd0;
#endif
  }
  cudasafeMalloc((void**)&GPUvars->ixptr,(Nf+1)*sizeof(int));
  cudaSafeMemcpy(GPUvars->ixptr, This->ixptr, (Nf+1)*sizeof(int), cudaMemcpyHostToDevice, "cudaMemcpy ixptr" );
  if(NE) {
    cudasafeMalloc((void**)&GPUvars->jxptr,(Nf+1)*sizeof(int));
    cudaSafeMemcpy(GPUvars->jxptr, This->jxptr, (Nf+1)*sizeof(int), cudaMemcpyHostToDevice, "cudaMemcpy jxptr" );
  }
  cudasafeMalloc((void**)&GPUvars->z0,Nt*sizeof(dcmplx));
  cudasafeMalloc((void**)&GPUvars->d0,Nt*sizeof(dcmplx));
  cudaSafeMemcpy(GPUvars->z0, z0, Nt*sizeof(dcmplx), cudaMemcpyHostToDevice, "cudaMemcpy z0" );
  cudaSafeMemcpy(GPUvars->d0, d0, Nt*sizeof(dcmplx), cudaMemcpyHostToDevice, "cudaMemcpy d0" );
#ifdef RADIALSHRINK
  cudaSafeMemcpy(GPUvars->dabs, dabs, Nt*sizeof(double), cudaMemcpyHostToDevice, "cudaMemcpy dabs" );
  mxFree(dabs);
#else
  // faster (?) to calculate abs(d) for all boxes once, do it here
  calcdabs((void *)GPUvars->d0,GPUvars->dabs,Nt);
  // *** one GPUvars!
  /* wrapper to allow for C++ compilation, the intent is

   calculatedabs<<<imin(MAXCONNECTIVITYBLOCKS,
                   (This->lptr[This->nlevel+1]+MAXCONNECTIVITYTHREADS-1)/
                     MAXCONNECTIVITYTHREADS),
                   MAXCONNECTIVITYTHREADS>>>
   ((dcmplx*)GPUvars->d0,GPUvars->dabs,This->lptr[This->nlevel+1]);
  */
#endif /*RADIALSHRINK*/
  mxFree(z0);
  mxFree(d0);

#endif /* !CUDASORT */

//   calcdabs((void *)GPUvars->d0,GPUvars->dabs,(Nf-1)/3+Nf);
  MPexp_connect_CUDA_(This,C,GPUvars,levm2p,maxm2p,cutoff,timing,printtime,meshonly);
#endif /* CUDASUPPORT */

  // conditions for when to sort on the CPU
#ifndef FULLCUDAFMM

#ifndef ASYMMETRIC_CPU
  C[0].jcptr[1] = 0;
  C[0].kcptr[0] = 0;
  C[0].kcptr[1] = 0; // note: 0 could be the final level!
#else
  C[0].jcptr[1] = 1; // direct interaction to self is included
  C[0].kcptr[0] = 1;
  C[0].kcptr[1] = 1;
#endif
  // if the GPU built the lists and the CPU also should do it, clean
  // the old values
#if defined(CUDASUPPORT) && !defined(VALIDATECUDACONNECTIVITY)
  for (int l = 1; l <= This->nlevel; l++) {
    mxFree(C[l].jcptr);
    C[l].jcptr=NULL;
  }
#endif
#if defined(CUDASORT) && defined(RADIALSHRINK)
  int Nt=(Nf-1)/3+Nf;
  SORT_REAL *dabs=(SORT_REAL*)mxMalloc(Nt*sizeof(SORT_REAL));
  cudaSafeMemcpy(dabs,GPUvars->dabs, Nt*sizeof(SORT_REAL), cudaMemcpyDeviceToHost, "cudaMemcpy dabs" );
  for (int j = 0; j < Nt; j++) {
    This->root[j].absd0=dabs[j];
  }
  mxFree(dabs);
#endif /*defined(CUDASORT) && defined(RADIALSHRINK)*/
#ifdef CHANNELPOT
  if(cparams->H>0)
    MPexp_connect_CPU_channelpot_((void*)This,(void*)C,cparams,levm2p,maxm2p,cutoff);
  else
    MPexp_connect_CPU_(This,C,levm2p,maxm2p,cutoff);
#else
  MPexp_connect_CPU_(This,C,levm2p,maxm2p,cutoff);
#endif /*CHANNELPOT*/
#endif /*!FULLCUDAFMM*/

#ifdef CUDATIME
  TimeEventRecord(sort2stop);
  TimeSyncPrintAndDestroy(sort2start,sort2stop,2,"GPU sort part 2");
#endif
#ifdef CUDASUPPORT
  cudaFreeDebug(GPUvars->dabs);
  cudaFreeDebug(GPUvars->d0);
  GPUvars->d0 = 0;
  GPUvars->dabs = 0;
#endif
}
/*------------------------------------------------------------------------*/
void MPexp_minmax_(double *xmin,double *xmax,double *ymin,double *ymax,
                   int N,const double *x,const double *y, int *validinput)
/* Finds min/max of two vectors. Note: it is assumed that suitable
   initial values in xmin/xmax are already set. */
{
  // loop unrolled to get rid of a few N/2 comparisions
  for (int i = 0; i < N-1; i += 2) {
    double min, max;

    if (x[i] < x[i+1]) {
      min = x[i];
      max = x[i+1];
    }
    else {
      min = x[i+1];
      max = x[i];
    }
    if (min < *xmin) *xmin = min;
    if (max > *xmax) *xmax = max;

    if (y[i] < y[i+1]) {
      min = y[i];
      max = y[i+1];
    }
    else {
      min = y[i+1];
      max = y[i];
    }
    if (min < *ymin) *ymin = min;
    if (max > *ymax) *ymax = max;
#ifdef CHECKNANINPUT
    if(!isfinite(y[i])||!isfinite(y[i+1])||!isfinite(x[i])||!isfinite(x[i+1]))
        *validinput=0;
#endif
  }
  if (N&1) {
    if (x[N-1] < *xmin)
      *xmin = x[N-1];
    else if (x[N-1] > *xmax)
      *xmax = x[N-1];

    if (y[N-1] < *ymin)
      *ymin = y[N-1];
    else if (y[N-1] > *ymax)
      *ymax = y[N-1];
#ifdef CHECKNANINPUT
    if(!isfinite(y[N-1])||!isfinite(x[N-1]))
        *validinput=0;
#endif
  }
}
/*------------------------------------------------------------------------*/
void MPexp_box_(MPexp *This,int N,int NE,
                const double *zr,const double *zi,
                const double *er,const double *ei,int* validinput)
/* Determines and sets the dimensions of the computational domain. */
{
  if (N) {
    This->xmin = zr[0];
    This->xmax = zr[0];
    This->ymin = zi[0];
    This->ymax = zi[0];
    MPexp_minmax_(&This->xmin,&This->xmax,
                  &This->ymin,&This->ymax,
                  N-1,zr+1,zi+1,validinput);

    // bounding box might become larger:
    if (NE)
      MPexp_minmax_(&This->xmin,&This->xmax,
                    &This->ymin,&This->ymax,
                    NE,er,ei,validinput);
  }
  else if (NE) {
    This->xmin = er[0];
    This->xmax = er[0];
    This->ymin = ei[0];
    This->ymax = ei[0];
    MPexp_minmax_(&This->xmin, &This->xmax,
                  &This->ymin, &This->ymax,
                  NE-1, er+1, ei+1,validinput);
  }
  else {
    // provide for the case of sorting only panels and no potentials
    This->xmin = 0.0;
    This->xmax = 0.0;
    This->ymin = 0.0;
    This->ymax = 0.0;
  }
}
/*------------------------------------------------------------------------*/
#ifndef SWAPI
#define SWAPI(X, Y) { int tmpi = X; X = Y; Y = tmpi; }
#endif

void MPexp_partition_(int begin,int *im,int end,double *z0,
                      int *ix,const double *z,int Nmax)
/* Partitions the set of points z[ix[begin]] through z[ix[end]]
  (exclusive) in two sets indicated by [begin,im) and [im,end) such
  that all points are <= or >= than z0 (which is both input and
  output). The division is done such that Nlarge <= Nmax, where Nlarge
  is the larger number of particles in the two lists. It is assumed
  that this split is possible and hence that Nmax is not too small.

  The routine is a workhorse and performance is important; for each
  box which is not a leaf it will be called three times. */
{
  int il,ir,left = begin,right = end;
  int Nsmall,Nlarge;

  // first split according to input z0 (geometric midpoint)
  {
    const double z0_ = *z0; // compiler: This is really a constant!
    il = left-1; ir = right;
    do il++; while (il < right && z[ix[il]] < z0_);
    do ir--; while (left <= ir && z0_ < z[ix[ir]]);
    if (il < ir) {
      SWAPI(ix[il],ix[ir])
      /* after the first swap, there is no need to check for index
         out of bounds */
      for ( ; ; ) {
        do il++; while (z[ix[il]] < z0_);
        do ir--; while (z0_ < z[ix[ir]]);
        if (il >= ir) break;
        SWAPI(ix[il],ix[ir])
      }
    }
  }

  // proceed with largest interval [left,il) or [il,right)
  Nsmall = il-begin;
  Nlarge = end-il;
  if (Nsmall > Nlarge) {
    Nlarge = Nsmall;
    *im = right = il;
  }
  else
    *im = left = il;

  while (Nlarge > Nmax) {
    if (right-left >= 3) {
      /* "median of 3": find pivot and create boundaries such that
         left <= left+1 <= right-1 */
      const int mid = (left+right)>>1;
      SWAPI(ix[left+1],ix[mid])
      if (z[ix[left]] > z[ix[right-1]])
        SWAPI(ix[left], ix[right-1])
      if (z[ix[left+1]] > z[ix[right-1]])
        SWAPI(ix[left+1],ix[right-1])
      else if (z[ix[left]] > z[ix[left+1]])
        SWAPI(ix[left],ix[left+1])
      const double z0_ = *z0 = z[ix[left+1]];

      for (il = left+1, ir = right-1; ; ) {
        do il++; while (z[ix[il]] < z0_);
        do ir--; while (z0_ < z[ix[ir]]);
        if (il >= ir) break;

        SWAPI(ix[il],ix[ir])
      }
      SWAPI(ix[left+1],ix[ir]) // ir is the correct place for the pivot

      // the pivot belongs to the shortest interval:
      Nsmall = ir-begin; // [begin,ir)
      Nlarge = end-ir-1; // [ir+1,end)
      if (Nsmall > Nlarge) {
        Nlarge = Nsmall;
        *im = right = ir;
      }
      else
        /* note: if Nsmall == Nlarge (can only happen if right-left is
           odd), then the right array will be shorter than the left by
           one */
        *im = left = ir+1;
    }
    else {
      // code above does not work properly for less than three elements
      if (z[ix[left]] > z[ix[right-1]])
        SWAPI(ix[left], ix[right-1])
      *im = right-1;
      *z0 = z[ix[right-1]];
      break;
    }
  }
}
/*------------------------------------------------------------------------*/
void MPexp_split_(int begin,int *im,int end,double z0,
                  int *ix,const double *z)
/* Split as above, but the coordinate z0 around which to split is
   already known. */
{
  int il = begin-1, ir = end;
  do il++; while (il < end && z[ix[il]] < z0);
  do ir--; while (begin <= ir && z0 < z[ix[ir]]);
  if (il < ir) {
    SWAPI(ix[il], ix[ir])
    for ( ; ; ) {
      do il++; while (z[ix[il]] < z0);
      do ir--; while (z0 < z[ix[ir]]);
      if (il >= ir) break;
      SWAPI(ix[il], ix[ir])
    }
  }
  *im = il;
}
/*-----------------------------------------------------------------------*/
void MPexp_partition_CPU_(MPexp *This,
                          int N,const double *zr,const double *zi,
                          int NE,const double *er,const double *ei,
                          const panel *panels,int Npanel,
                          int **panelptrlist,
                          int Ndirect,double cutoff,
                          int *ix,int *ixptr,int *ixptr0,
                          int *jx,int *jxptr,int *jxptr0,
                          channelparam* cparams)
/* Partitioning of potentials, evaluation points, and panels. This
   function creates the multipole mesh on the CPU. */
{
  int *ixptr1,*jxptr1;
  const int Nf = This->lptr[This->nlevel+1]-This->lptr[This->nlevel];

  // (ixptr,ix) and (jxptr,jx) to be built
  for (int i = 0; i < N; i++) ix[i] = i;
  if(NE) for (int i = 0; i < NE; i++) jx[i] = i;

  // panel variables
  int *tmpptr1 = (int *)mxMalloc(Npanel*sizeof(int));
  int *tmpptr2 = (int *)mxMalloc(Npanel*sizeof(int));
  int panelsum = Npanel,panelptrindex;
  int outcount1,outcount2;

  panelptrlist[0] = (int *)mxMalloc(Npanel*sizeof(int));
  for(int l = 0; l < Npanel; l++)
    panelptrlist[0][l] = l; // root level contains all panels

  This->root[0].panelptr = panelptrlist[0];
  This->root[0].npanel = Npanel;

  // root box
  This->root[0].z0 = 0.5*(This->xmin+This->xmax)+I*0.5*(This->ymin+This->ymax);
  This->root[0].d0 = 0.5*(This->xmax-This->xmin)+I*0.5*(This->ymax-This->ymin);

  // dummies unless used
  panel *smallerpanels = NULL;
  double *dummyvortices = NULL;

#ifdef RADIALSHRINK
  // shrink never performed on the root; explicitly set absd0
  This->root[0].absd0 = cabs(This->root[0].d0);

  // variable smallerpanels is used to temporary save the panel array
  // (this avoids calling panelinbox() twice)
  smallerpanels = (panel *)mxMalloc(Npanel*sizeof(panel));
#endif
#ifdef PANELSORT
  dummyvortices = (double *)mxMalloc(Npanel*2*sizeof(double));
#endif

  int startlevel=1;
  int panelmultiplier;
#ifdef CHANNELPOT
  if(cparams->H>0) {
    channelsort((void*)This,cparams,zr,er,ix,jx,&ixptr[Nf-cparams->Nb],
                &jxptr[Nf-cparams->Nb],This->xmin, N, NE);
    ixptr0=&ixptr[Nf-cparams->Nb];
    if (NE) jxptr0=&jxptr[Nf-cparams->Nb];
    startlevel=cparams->Nlev/2+1;
    channelsortpanel(This,panels,tmpptr1,tmpptr2,panelptrlist,startlevel);
    panelmultiplier=3;
//     mexPrintf("startlevel=%d nlevel=%d\n",startlevel,This->nlevel);
//     for (int i = 0; i < This->lptr[startlevel+1]; i++)
//       mexPrintf("Box %d npanels=%d\n",i,This->root[i].npanel);

  }
#endif
  // loop over levels
  for (int l = startlevel,Nmax = N,Nmax0; l <= This->nlevel; l++) {
    mpexp *parent = &This->root[This->lptr[l-1]];
    mpexp *child = &This->root[This->lptr[l]];
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];
// mexPrintf("N0=%d N1=%d Nf=%d 4(N0-1)+3=%d\n",N0,N1,Nf,4*(N0-1)+3);
    // maximum number of panels this round
    panelptrlist[l] = (int *)mxMalloc(panelsum*panelmultiplier*4*sizeof(int));
    // the number of panels in the lower level to save memory in the
    // next allocation
    panelsum = 0;
    panelptrindex = 0; // the index in the array

    /* divide each box on level l into 4 children such that about
       1/4th of the potentials are located in each child box. */

    /* expanding arrays; they are written at in the lower part of the
       _same_ memory area (parents use the "0-part" and children the
       "1-part") */
    ixptr1 = &ixptr[Nf-N1];
    if (NE) jxptr1 = &jxptr[Nf-N1];

    /* determine the maximum number of coordinates after each split
       (the map is x --> ceil(0.5*x)) */
    if ((Nmax0 = (Nmax >> 1)+(Nmax&1)) <= Ndirect)
      Nmax0 = Nmax = Ndirect;
    else if ((Nmax = (Nmax0 >> 1)+(Nmax0&1)) <= Ndirect)
      Nmax = Ndirect;
    /* the intent with using Ndirect as a cutoff is to save some work
       at the final level whenever the number of particles is already
       sufficiently low */

    // loop over all parent boxes
    for (int i = 0, j = 0; i < N0; i++, j += 4) {
      dcmplx z,d;
      dcmplx z0,d0,z1,d1;
      dcmplx z00,d00,z01,d01,z10,d10,z11,d11;
      double split;

      // box to be split
      z = parent[i].z0;
      d = parent[i].d0;
      /* partition the coordinates [ixptr0[i],ixptr[i+1]) in box i
         into 4 sections indexed by [ixptr1[4*i+n],ixptr1[4*i+n+1])
         for n = 0..3 */
      ixptr1[4*i] = ixptr0[i];
      if (NE) jxptr1[4*i] = jxptr0[i];

      /* cut across the longest side if it is a certain factor longer,
         otherwise take a y-cut followed by two x-cuts; this is the
         Lebesgue ordering */
      if (creal(d) > meshratio*cimag(d)) {
        split = creal(z);
#ifdef PANELSORT
        // create the dummy list and call the other partition function
        createdummylist(dummyvortices,
                        panels,parent[i].panelptr,parent[i].npanel,1,z,d);
        MPexp_partition_dummy(ixptr0[i],&ixptr1[4*i+2],ixptr0[i+1],
                              &split,ix,zr,(ixptr0[i+1]-ixptr0[i])/2,
                              dummyvortices,parent[i].npanel*2);
#else
        MPexp_partition_(ixptr0[i],&ixptr1[4*i+2],ixptr0[i+1],
                         &split,ix,zr,Nmax0);
#endif
        COMPLEXASSIGN(z0,0.5*(creal(z)-creal(d)+split),cimag(z));
        d0 = z0-(z-d);
        COMPLEXASSIGN(z1,0.5*(creal(z)+creal(d)+split),cimag(z));
        d1 = (z+d)-z1;
        if (NE)
          MPexp_split_(jxptr0[i],&jxptr1[4*i+2],jxptr0[i+1],split,jx,er);

        outcount1 = outcount2 = 0;
        if (parent[i].npanel)
          MPexp_split_panel(z,d,panels,parent[i].panelptr,
                            tmpptr1,tmpptr2,parent[i].npanel,
                            &outcount1,&outcount2,split,1);
      }
      else {
        split = cimag(z);
#ifdef PANELSORT
        createdummylist(dummyvortices,
                        panels,parent[i].panelptr,parent[i].npanel,0,z,d);
        MPexp_partition_dummy(ixptr0[i],&ixptr1[4*i+2],ixptr0[i+1],
                              &split,ix,zi,(ixptr0[i+1]-ixptr0[i])/2,
                              dummyvortices,parent[i].npanel*2);
#else
        MPexp_partition_(ixptr0[i],&ixptr1[4*i+2],ixptr0[i+1],
                         &split,ix,zi,Nmax0);
#endif
        COMPLEXASSIGN(z0,creal(z),0.5*(cimag(z)-cimag(d)+split));
        d0 = z0-(z-d);
        COMPLEXASSIGN(z1,creal(z),0.5*(cimag(z)+cimag(d)+split));
        d1 = (z+d)-z1;
        if (NE)
          MPexp_split_(jxptr0[i],&jxptr1[4*i+2],jxptr0[i+1],split,jx,ei);

        outcount1 = outcount2 = 0;
        if (parent[i].npanel)
          MPexp_split_panel(z,d,panels,parent[i].panelptr,
                            tmpptr1,tmpptr2,parent[i].npanel,
                            &outcount1,&outcount2,split,0);
      }

      // set up the panel pointers
      child[j].panelptr = panelptrlist[l]+panelptrindex;
      child[j+1].panelptr = panelptrlist[l]+panelptrindex+outcount1;
      panelptrindex += outcount1*2;
      child[j+2].panelptr = panelptrlist[l]+panelptrindex;
      child[j+3].panelptr = panelptrlist[l]+panelptrindex+outcount2;
      panelptrindex += outcount2*2;

      // same procedure for the resulting two boxes...
      if (meshratio*creal(d0) > cimag(d0)) {
        split = creal(z0);
#ifdef PANELSORT
        createdummylist(dummyvortices,panels,tmpptr1,outcount1,1,z0,d0);
        MPexp_partition_dummy(ixptr0[i],&ixptr1[4*i+1],ixptr1[4*i+2],
                              &split, ix,zr,(ixptr1[4*i+2]-ixptr0[i])/2,
                              dummyvortices,outcount1*2);
#else
        MPexp_partition_(ixptr0[i],&ixptr1[4*i+1],ixptr1[4*i+2],
                         &split,ix,zr,Nmax);
#endif
        COMPLEXASSIGN(z00,0.5*(creal(z0)-creal(d0)+split),cimag(z0));
        d00 = z00-(z0-d0);
        COMPLEXASSIGN(z01,0.5*(creal(z0)+creal(d0)+split),cimag(z0));
        d01 = (z0+d0)-z01;
        if (NE)
          MPexp_split_(jxptr0[i],&jxptr1[4*i+1],jxptr1[4*i+2],split,jx,er);

        if (outcount1)
          MPexp_split_panel(z0,d0,panels,tmpptr1,
                            child[j].panelptr,child[j+1].panelptr,
                            outcount1,
                            &(child[j].npanel),&(child[j+1].npanel),split,1);
        else
          child[j].npanel = child[j+1].npanel = 0;
      }
      else {
        split = cimag(z0);
#ifdef PANELSORT
        createdummylist(dummyvortices,panels,tmpptr1,outcount1,0,z0,d0);
        MPexp_partition_dummy(ixptr0[i],&ixptr1[4*i+1],ixptr1[4*i+2],
                              &split,ix,zi,(ixptr1[4*i+2]-ixptr0[i])/2,
                              dummyvortices,outcount1*2);
#else
        MPexp_partition_(ixptr0[i],&ixptr1[4*i+1],ixptr1[4*i+2],
                         &split,ix,zi,Nmax);
#endif
        COMPLEXASSIGN(z00,creal(z0),0.5*(cimag(z0)-cimag(d0)+split));
        d00 = z00-(z0-d0);
        COMPLEXASSIGN(z01,creal(z0),0.5*(cimag(z0)+cimag(d0)+split));
        d01 = (z0+d0)-z01;
        if (NE)
          MPexp_split_(jxptr0[i],&jxptr1[4*i+1],jxptr1[4*i+2],split,jx, ei);

        if (outcount1)
          MPexp_split_panel(z0,d0,panels,tmpptr1,child[j].panelptr,
                            child[j+1].panelptr,outcount1,
                            &(child[j].npanel),&(child[j+1].npanel),split,0);
        else
          child[j].npanel = child[j+1].npanel = 0;
      }
      // and again...
      if (meshratio*creal(d1) > cimag(d1)) {
        split = creal(z1);
#ifdef PANELSORT
        createdummylist(dummyvortices,panels,tmpptr2,outcount2,1,z1,d1);
        MPexp_partition_dummy(ixptr1[4*i+2],&ixptr1[4*i+3],ixptr0[i+1],
                              &split,ix,zr,(ixptr0[i+1]-ixptr1[4*i+2])/2,
                              dummyvortices,outcount1*2);
#else
        MPexp_partition_(ixptr1[4*i+2],&ixptr1[4*i+3],ixptr0[i+1],
                         &split,ix,zr,Nmax);
#endif
        COMPLEXASSIGN(z10,0.5*(creal(z1)-creal(d1)+split),cimag(z1));
        d10 = z10-(z1-d1);
        COMPLEXASSIGN(z11,0.5*(creal(z1)+creal(d1)+split),cimag(z1));
        d11 = (z1+d1)-z11;
        if (NE)
          MPexp_split_(jxptr1[4*i+2],&jxptr1[4*i+3],jxptr0[i+1],split,jx,er);

        if (outcount2)
          MPexp_split_panel(z1,d1,panels,tmpptr2,
                            child[j+2].panelptr,child[j+3].panelptr,
                            outcount2,
                            &(child[j+2].npanel),&(child[j+3].npanel),split,1);
        else
          child[j+2].npanel = child[j+3].npanel=0;
      }
      else {
        split = cimag(z1);
#ifdef PANELSORT
        createdummylist(dummyvortices,panels,tmpptr2,outcount2,0,z1,d1);
        MPexp_partition_dummy(ixptr1[4*i+2],&ixptr1[4*i+3],ixptr0[i+1],
                              &split,ix,zi,(ixptr0[i+1]-ixptr1[4*i+2])/2,
                              dummyvortices,outcount1*2);
#else
        MPexp_partition_(ixptr1[4*i+2],&ixptr1[4*i+3],ixptr0[i+1],
                         &split,ix,zi,Nmax);
#endif
        COMPLEXASSIGN(z10,creal(z1),0.5*(cimag(z1)-cimag(d1)+split));
        d10 = z10-(z1-d1);
        COMPLEXASSIGN(z11,creal(z1),0.5*(cimag(z1)+cimag(d1)+split));
        d11 = (z1+d1)-z11;
        if (NE)
          MPexp_split_(jxptr1[4*i+2],&jxptr1[4*i+3],jxptr0[i+1],split,jx,ei);

        if (outcount2)
          MPexp_split_panel(z1,d1,panels,tmpptr2,
                            child[j+2].panelptr,child[j+3].panelptr,
                            outcount2,&(child[j+2].npanel),
                            &(child[j+3].npanel),split,0);
        else
          child[j+2].npanel = child[j+3].npanel = 0;
      }
      panelsum += child[j].npanel+child[j+1].npanel+
        child[j+2].npanel+child[j+3].npanel;

      // produced boxes
      child[j].z0 = z00;
      child[j].d0 = d00;
      child[j+1].z0 = z01;
      child[j+1].d0 = d01;
      child[j+2].z0 = z10;
      child[j+2].d0 = d10;
      child[j+3].z0 = z11;
      child[j+3].d0 = d11;

#ifdef PANELSHRINKBOX
      if(NE) {
        panelshrinkbox(child+j,zr,zi,er,ei,ix,jx,
                       ixptr1[4*i],ixptr1[4*i+1],jxptr1[4*i],jxptr1[4*i+1],
                       panels,smallerpanels,cutoff);
        panelshrinkbox(child+j+1,zr,zi,er,ei,ix,jx,
                       ixptr1[4*i+1],ixptr1[4*i+2],jxptr1[4*i+1],jxptr1[4*i+2],
                       panels,smallerpanels,cutoff);
        panelshrinkbox(child+j+2,zr,zi,er,ei,ix,jx,
                       ixptr1[4*i+2],ixptr1[4*i+3],jxptr1[4*i+2],jxptr1[4*i+3],
                       panels,smallerpanels,cutoff);
        panelshrinkbox(child+j+3,zr,zi,er,ei,ix,jx,
                       ixptr1[4*i+3],ixptr0[i+1],jxptr1[4*i+3],jxptr0[i+1],
                       panels,smallerpanels,cutoff);
      }
      else {
        panelshrinkbox(child+j,zr,zi,NULL,NULL,ix,jx,
                       ixptr1[4*i],ixptr1[4*i+1],0,0,
                       panels,smallerpanels,cutoff);
        panelshrinkbox(child+j+1,zr,zi,NULL,NULL,ix,jx,
                       ixptr1[4*i+1],ixptr1[4*i+2],0,0,
                       panels,smallerpanels,cutoff);
        panelshrinkbox(child+j+2,zr,zi,NULL,NULL,ix,jx,
                       ixptr1[4*i+2],ixptr1[4*i+3],0,0,
                       panels,smallerpanels,cutoff);
        panelshrinkbox(child+j+3,zr,zi,NULL,NULL,ix,jx,
                       ixptr1[4*i+3],ixptr0[i+1],0,0,
                       panels,smallerpanels,cutoff);
      }
#endif /* PANELSHRINKBOX */
    }

    // arrays grow to the left
    ixptr0 = ixptr1;
    if (NE) jxptr0 = jxptr1;
  }
#ifdef CHANNELPOT
  if(cparams->H>0)
    replicatezd((void*)This,cparams); //create z and d for mirror boxes
#endif
  // deallocate
  mxFree(dummyvortices);
  mxFree(smallerpanels);
  mxFree(tmpptr2);
  mxFree(tmpptr1);
}
/*-----------------------------------------------------------------------*/
#ifndef SWAPI
#define SWAPI(X,Y) { int tmpi = X; X = Y; Y = tmpi; }
#endif

#ifndef ASYMMETRIC_CPU
void MPexp_connect_CPU_(MPexp *This,mpSparse *C,
                        int *levm2p,int *maxm2p,double cutoff)
/* Once the mesh has been obtained, this function determines the type
   of connection between all multipole boxes. Executed on the CPU. */
{
  // loop over all levels
  for (int l = 1,nc = 0; l <= This->nlevel; l++) {
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];
    mpexp *child = &This->root[This->lptr[l]];

    /* allocate connectivity matrix using that all near connections to
       higher index boxes at the previous level implies 4 connections
       for each of 4 children at this level (and 6 additional
       connections to higher indices {0->[1 2 3],1->[2 3],2->[3]}
       between the children of each parent) */
    C[l].jcptr = (int *)mxMalloc(((2+(l == This->nlevel))*N1+1+
                                  (6*N0+4*4*nc))*sizeof(int));
    C[l].jcptr[0] = 0;
    C[l].kcptr = &C[l].jcptr[N1+1];
    C[l].ir = &C[l].kcptr[(1+(l == This->nlevel))*N1];


    // setup connections
    nc = 0;  // counter for near connections
    for (int i = 0, j = 0; i < N0; i++)  // parent i
      for (int c = 0; c < 4; c++, j++) { // children j
        // forward- and backward counters in each column
        int fwd = C[l].jcptr[j];
        int bwd = fwd+3-c+4*(C[l-1].kcptr[i]-C[l-1].jcptr[i]);
        C[l].jcptr[j+1] = bwd;

        /* near connections to siblings > j: the connection could in
           principle be of the far type (except to the neighbour
           sibling), but this should occur only very rarely */

        for (int d = 1; d < 4-c; d++)
#if !defined(CUDASUPPORT) && !defined(M2PSINBOX)
          C[l].ir[fwd++] = j+d;
#else
        if (!mpexp_theta(&child[j], &child[j+d], cutoff))
          // near connection
          C[l].ir[fwd++] = j+d;
        else
          // far connection
          C[l].ir[--bwd] = j+d;
#endif /*CUDASUPPORT*/

        // parents (near) connections
        for (int ii = C[l-1].jcptr[i]; ii < C[l-1].kcptr[i]; ii++) {
          /* the children of the near connections of box i become
             either near connections or far connections of box j */
          int k = C[l-1].ir[ii] << 2; // 0-child of connection
          // (note: minor differences to the C99-code in that
          // interactions through parents are not supported here)

          for (int d = 0; d < 4; d++, k++)
            // connection according to theta criterion
            if (!mpexp_theta(&child[j],&child[k],cutoff))
              // near connection
              C[l].ir[fwd++] = k;
            else
              // far connection
              C[l].ir[--bwd] = k;
        }
        C[l].kcptr[j] = bwd;
        /* this splits the column into near connections
           [jcptr[j],kcptr[j]) and far connections
           [kcptr[j],jcptr[j+1]) */
        nc += bwd-C[l].jcptr[j];
        if (C[l].jcptr[j+1]-bwd > *maxm2p)
          *maxm2p = C[l].jcptr[j+1]-bwd;
        if (l < *levm2p && *maxm2p > 0)
          *levm2p = l;

        // last level: sort the near connections one step further
        if (l == This->nlevel) {
          int dir1,dir2; // 1: cluster to particle, -1: particle to cluster
          fwd = C[l].jcptr[j]-1;
          bwd = C[l].kcptr[j];
#if defined(CUDASUPPORT) || defined(M2PSINBOX)
          if (j == 0)
            bwd--; // skip box 0
          else {
            if (C[l].ir[fwd+1] == 0)
              fwd++;
#endif
            for ( ; ; ) {
              do
                fwd++;
              while (fwd < bwd &&
                      !(dir1 = mpexp_theta2(&child[j], &child[C[l].ir[fwd]],
                      cutoff)));

              do {
                bwd--;
                if (fwd <= bwd &&
                        (dir2 = mpexp_theta2(&child[j], &child[C[l].ir[bwd]],
                        cutoff)))
                  C[l].ir[bwd] *= dir2;
                else
                  break;
              } while (1);

              if (fwd >= bwd) break;
              SWAPI(C[l].ir[fwd], C[l].ir[bwd])
              C[l].ir[bwd] *= dir1;
            }
#if defined(CUDASUPPORT) || defined(M2PSINBOX)
          }
#endif
          C[l].kcptr[j+N1] = bwd+1;
          /* this splits the near connections into 'truly near'
             [jcptr[j],kcptr[j+N]) and 'less near'
             [kcptr[j+N],kcptr[j]) connections */
        }
      }
  }
}
/*-----------------------------------------------------------------------*/
#else /* ASYMMETRIC_CPU */

void MPexp_connect_CPU_(MPexp *This,mpSparse *C,
                        int *levm2p,int *maxm2p,double cutoff)
/* Once the mesh has been obtained, this function determines the type
   of connection between all multipole boxes. Executed on the CPU. */
{
  // loop over all levels (nc = 1 since the root box is strongly
  // connected to itself)
  for (int l = 1,nc = 1; l <= This->nlevel; l++) {
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];
    mpexp *child = &This->root[This->lptr[l]];

    /* allocate connectivity matrix using that all near connections at
       the previous level implies 4 connections for each of 4 children
       at this level */

#ifdef VALIDATECUDACONNECTIVITY
    int jcptrlen=((2+(l == This->nlevel))*N1+1+(4*4*nc));
    int *jcptrbackup=(int *)mxMalloc(jcptrlen*sizeof(int));
    memcpy(jcptrbackup,C[l].jcptr,jcptrlen*sizeof(int));
#endif

    C[l].jcptr = (int *)mxMalloc(((2+(l == This->nlevel))*N1+1+
                                  (4*4*nc))*sizeof(int));
    C[l].jcptr[0] = 0;
    C[l].kcptr = &C[l].jcptr[N1+1];
    C[l].ir = &C[l].kcptr[(1+(l == This->nlevel))*N1];

    // setup connections
    nc = 0;  // counter for near connections
    for (int i = 0, j = 0; i < N0; i++)  // parent i
      for (int c = 0; c < 4; c++, j++) { // children j
        // forward- and backward counters in each column
        int fwd = C[l].jcptr[j];
        int bwd = fwd+4*(C[l-1].kcptr[i]-C[l-1].jcptr[i]);
        C[l].jcptr[j+1] = bwd;

        // parents (near) connections (includes the parent itself)
        for (int ii = C[l-1].jcptr[i]; ii < C[l-1].kcptr[i]; ii++) {
          /* the children of the near connections of box i become
             either near connections or far connections of box j */
          int k = C[l-1].ir[ii] << 2; // 0-child of connection
          // (note: minor differences to the C99-code in that
          // interactions through parents are not supported here)

          for (int d = 0; d < 4; d++, k++)
            // connection according to theta criterion
            if (!mpexp_theta(&child[j],&child[k],cutoff))
              // near connection
              C[l].ir[fwd++] = k;
            else
              // far connection
              C[l].ir[--bwd] = k;
        }
        C[l].kcptr[j] = bwd;
        /* this splits the column into near connections
           [jcptr[j],kcptr[j]) and far connections
           [kcptr[j],jcptr[j+1]) */
        nc += bwd-C[l].jcptr[j];
        if (C[l].jcptr[j+1]-bwd > *maxm2p)
          *maxm2p = C[l].jcptr[j+1]-bwd;
        if (l < *levm2p && *maxm2p > 0)
          *levm2p = l;

        // last level: sort the near connections one step further
        if (l == This->nlevel) {
          int dir1,dir2; // 1: cluster to particle, -1: particle to cluster
          fwd = C[l].jcptr[j]-1;
          bwd = C[l].kcptr[j];
          // since the format uses the sign of the index to denote the
          // type of the near connection, box 0 must be skipped
          if (j == 0)
            bwd--; // skip box 0
          else {
            // again, skip box 0 (must be first since column is sorted)
            if (C[l].ir[fwd+1] == 0)
              fwd++;
            for ( ; ; ) {
              do
                fwd++;
              while (fwd < bwd &&
                     !(dir1 = mpexp_theta2(&child[j],&child[C[l].ir[fwd]],
                                           cutoff)));

              do {
                bwd--;
                if (fwd <= bwd &&
                    (dir2 = mpexp_theta2(&child[j],&child[C[l].ir[bwd]],
                                         cutoff)))
                  C[l].ir[bwd] *= dir2;
                else
                  break;
              } while (1);

              if (fwd >= bwd) break;
              SWAPI(C[l].ir[fwd],C[l].ir[bwd])
              C[l].ir[bwd] *= dir1;
            }
          }
          C[l].kcptr[j+N1] = bwd+1;
          /* this splits the near connections into 'truly near'
             [jcptr[j],kcptr[j+N]) and 'less near'
             [kcptr[j+N],kcptr[j]) connections */
        }
      }
#ifdef VALIDATECUDACONNECTIVITY
    if(memcmp(C[l].jcptr,jcptrbackup,jcptrlen*sizeof(int)))
      mexPrintf("Warning, connectivity different between CPU and GPU\n");
    mxFree(jcptrbackup);
#endif

  }
}
#endif /* ASYMMETRIC_CPU */
/*-----------------------------------------------------------------------*/
#ifdef CUDASUPPORT
void MPexp_partition_post_CUDA_(MPexp* This,cudavariables* GPUvars,
                                int *ix,int *ixptr,int *jx,int *jxptr,
                                int N,int NE,double* timing,bool printtime)
/* After CUDA_perform_partitioning has finished, this function copies
   some data back to the CPU and performs some cleanup. */
{
  CudaStartTime(printtime,NULL)

  // note: all data are copied back to the CPU here without regards to
  // whether it is actually used or not

  // copy variables back to the CPU after the GPU sorting
  SORT_DCMPLX *ztmp = (SORT_DCMPLX *)mxMalloc(This->lptr[This->nlevel+1]*sizeof(SORT_DCMPLX));
  SORT_DCMPLX *dtmp = (SORT_DCMPLX *)mxMalloc(This->lptr[This->nlevel+1]*sizeof(SORT_DCMPLX));
  cudaSafeMemcpy(ixptr,GPUvars->ixptr,
                 (This->lptr[This->nlevel+1]-This->lptr[This->nlevel]+1)*
                 sizeof(int),cudaMemcpyDeviceToHost,"cudaMemcpy ixptr" );
  cudaSafeMemcpy(ix,GPUvars->ix,N*sizeof(int),
                 cudaMemcpyDeviceToHost,"cudaMemcpy ix");
  if(NE) {
    cudaSafeMemcpy(jxptr,GPUvars->jxptr,
                   (This->lptr[This->nlevel+1]-This->lptr[This->nlevel]+1)*
                   sizeof(int),cudaMemcpyDeviceToHost,"cudaMemcpy jxptr");
    cudaSafeMemcpy(jx,GPUvars->jx,NE*sizeof(int),
                   cudaMemcpyDeviceToHost,"cudaMemcpy ix");
  }
  cudaSafeMemcpy(ztmp,GPUvars->z0,This->lptr[This->nlevel+1]*sizeof(SORT_DCMPLX),
                 cudaMemcpyDeviceToHost,"cudaMemcpy z0");
  cudaSafeMemcpy(dtmp,GPUvars->d0,This->lptr[This->nlevel+1]*sizeof(SORT_DCMPLX),
                 cudaMemcpyDeviceToHost,"cudaMemcpy d0");
  for (int j = 0; j < This->lptr[This->nlevel+1]; j++) {
    COMPLEXASSIGN(This->root[j].z0,((SORT_REAL*)&ztmp[j])[0],((SORT_REAL*)&ztmp[j])[1]);
    COMPLEXASSIGN(This->root[j].d0,((SORT_REAL*)&dtmp[j])[0],((SORT_REAL*)&dtmp[j])[1]);
    This->root[j].npanel = 0; // panels not handled by the GPU
  }
  mxFree(dtmp);
  mxFree(ztmp);

  CudaStopTime(printtime,NULL)
  CudaPrintTime("Time taken to copy elements back to CPU after sort",
                NULL,0,printtime)


}
/*-----------------------------------------------------------------------*/
void MPexp_connect_CUDA_(MPexp *This,mpSparse *C,
                         cudavariables *GPUvars,
                         int *levm2p,int *maxm2p,double cutoff,
                         double *timing,bool printtime,bool meshonly)
/* Like MPexp_connect_CPU_, but with the main work executed on the
   GPU. */
{
  // copy the first level of connectivity information to the GPU
  cudasafe(cudaMallocDebug((void **)&GPUvars->connect[0].jcptr,
                           5*sizeof(int)),"cudaMalloc GPUvars->connect");
  cudaSafeMemcpy(GPUvars->connect[0].jcptr,C[0].jcptr,
                 5*sizeof(int),cudaMemcpyHostToDevice,
                 "cudaMemcpy GPUvars->connect" );
  GPUvars->connect->kcptr = &GPUvars->connect[0].jcptr[2];
  GPUvars->connect->ir = &GPUvars->connect[0].kcptr[2];

  int *coutput; // = [maxm2p nc]
  cudasafe(cudaMallocDebug((void **)&coutput,2*sizeof(int)),
           "cudaMalloc coutput");

  // loop over levels (nc = 1 since the root box is strongly
  // connected to itself)
//   int jxsize=This->lptr[This->nlevel+1]-This->lptr[This->nlevel]+1;
//   int* testlist2= (int *)mxMalloc(jxsize*sizeof(int));
//   cudaSafeMemcpy(testlist2,GPUvars->jxptr,jxsize*sizeof(int),
//           cudaMemcpyDeviceToHost,"cudaMemcpy GPUvars->connect");
//   for(int i=0;i<jxsize;i++)
// //       if(testlist[i]!=C[l].jcptr[i])
// //         mexPrintf("level %d list differs at position %d, cpu: %d gpu: %d\n",l,i,C[l].jcptr[i],testlist[i]);
//     mexPrintf("jxptr[%d]=%d\n",i,testlist2[i]);
//
//   mxFree(testlist2);

  for (int l = 1,nc = 1; l <= This->nlevel; l++) {
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];

    /* allocate connectivity matrix using that all near connections at
       the previous level implies 4 connections for each of 4 children
       at this level */
    int jcptrlen = ((2+(l == This->nlevel))*N1+1+(4*4*nc));
    if(!fullcudafmm||meshonly) {
      C[l].jcptr = (int *)mxMalloc(jcptrlen*sizeof(int));
      C[l].jcptr[0] = 0;
      C[l].kcptr = &C[l].jcptr[N1+1];
      C[l].ir = &C[l].kcptr[(1+(l == This->nlevel))*N1];
    }
    else {
      C[l].jcptr=NULL;
      C[l].kcptr=NULL;
      C[l].ir=NULL;
    }
#ifdef CUDATIMESORT
    cudaEvent_t start,stop;
    cudaTimingCreateAndStart(&start,&stop,"start","stop");
//     TimeCreateAndStart(start,stop);
#endif

    // this step of the connectivity search cannot easily be done on
    // the GPU, since each loop depends on previous steps
//     for (int i = 0, j = 0; i < N0; i++) // parent i
//       for (int c = 0; c < 4; c++, j++) { // children j
//         // forward- and backward counters in each column
//         int fwd = C[l].jcptr[j];
//         int bwd = fwd+4*(C[l-1].kcptr[i]-C[l-1].jcptr[i]);
//         C[l].jcptr[j+1] = bwd;
//       }

    // perform allocations for the GPU connecticity search
    cudasafe(cudaMallocDebug((void **)&GPUvars->connect[l].jcptr,
                             jcptrlen*sizeof(int)),
             "cudaMalloc GPUvars->connect");

    cudaSafeMemset(GPUvars->connect[l].jcptr,0,sizeof(int),"cudaMemset jcptr");
    cumsumlist(GPUvars->connect[l-1].jcptr,GPUvars->connect[l-1].kcptr,GPUvars->connect[l].jcptr+1,N0*4,GPUvars,(This->nlevel-l)*2);

    //testing code
//     int* testlist= (int *)mxMalloc(jcptrlen*sizeof(int));
//     cudaSafeMemcpy(testlist,GPUvars->connect[l].jcptr,(N1+1)*sizeof(int),
//                    cudaMemcpyDeviceToHost,"cudaMemcpy GPUvars->connect");
//     for(int i=0;i<N1+1;i++)
// //       if(testlist[i]!=C[l].jcptr[i])
// //         mexPrintf("level %d list differs at position %d, cpu: %d gpu: %d\n",l,i,C[l].jcptr[i],testlist[i]);
//       mexPrintf("level %d at position %d gpu: %d\n",l,i,testlist[i]);
//
//     mxFree(testlist);


    //end testcode

//     cudaSafeMemcpy(GPUvars->connect[l].jcptr,C[l].jcptr,(N1+1)*sizeof(int),
//                    cudaMemcpyHostToDevice,"cudaMemcpy GPUvars->connect");
    GPUvars->connect[l].kcptr = &GPUvars->connect[l].jcptr[N1+1];
    GPUvars->connect[l].ir =
      &GPUvars->connect[l].kcptr[(1+(l == This->nlevel))*N1];
    cudaSafeMemset(coutput,0,2*sizeof(int),"cudaMemset coutput");

    // the CUDA call for this level
    cudaCreateConnectivity(GPUvars->connect[l].jcptr, // *** one GPUvars!
                           GPUvars->connect[l].kcptr,
                           GPUvars->connect[l].ir,
                           GPUvars->connect[l-1].jcptr,
                           GPUvars->connect[l-1].kcptr,
                           GPUvars->connect[l-1].ir,N1,*maxm2p,
                           (void *)(((SORT_DCMPLX *)GPUvars->z0)+This->lptr[l]),
                           GPUvars->dabs+This->lptr[l],cutoff,
                           l == This->nlevel,coutput
                           DEBUGVECTORSTRING2);
    /* wrapper to allow for C++ compilation, the intent is

    cudacreateconnectivity<<<imin(MAXCONNECTIVITYBLOCKS,
                             (N1+MAXCONNECTIVITYTHREADS-1)/
                               MAXCONNECTIVITYTHREADS),
                             MAXCONNECTIVITYTHREADS>>>
      (GPUvars->connect[l].jcptr,
      GPUvars->connect[l].kcptr,
      GPUvars->connect[l].ir,
      GPUvars->connect[l-1].jcptr,
      GPUvars->connect[l-1].kcptr,
      GPUvars->connect[l-1].ir,N1,*maxm2p,
      ((dcmplx *)GPUvars->z0)+This->lptr[l],G
      PUvars->dabs+This->lptr[l],cutoff,
      l == This->nlevel,coutput,
      GPUvars->debugvector);
    */

    // copy of coutput at the CPU
    int output[2]; // = [maxm2p nc]
    cudaSafeMemcpy(output,coutput,2*sizeof(int),
                   cudaMemcpyDeviceToHost,"cudaMemcpy coutput");

#ifdef CUDATIMESORT
    cudaSafeEventRecord(stop,"stop 1116");
    cudaTimingSyncPrintAndDestroy(start,stop,NULL,0,printtime,
                                  "connectivity matrix","start","stop");
#endif
    if(!fullcudafmm||meshonly)
//       jcptrlen=2*N1+1;
      cudaSafeMemcpy(C[l].jcptr,GPUvars->connect[l].jcptr,
                     jcptrlen*sizeof(int),cudaMemcpyDeviceToHost,
                     "cudaMemcpy jcptr");

    if (l < *levm2p && output[0] > 0)
      *levm2p = l;
    *maxm2p = output[0];
    nc = output[1];
  }
  cudaFreeDebug(coutput);
}
#endif /* CUDASUPPORT */
/*------------------------------------------------------------------------*/
