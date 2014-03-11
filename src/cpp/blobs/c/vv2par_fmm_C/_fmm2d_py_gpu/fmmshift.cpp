/* fmmshift.cpp */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#include <stdio.h>
#endif

/* S. Engblom 2010-02-18 */
#include <string.h>
#include "fmmshift.h"
#include "hornershift.h"
#ifdef SSE3INTRINSIC
#include "pmmintrin.h"
#elif defined(SSE2INTRINSIC)
#include "emmintrin.h"
#endif
// deallocate persistent data
void shift_atExit(void);

// code responds to #define INLINECOMPLEX

// workspace arrays
static dcmplx *wksp,*rpow;

// scalars defining shift-type and sizes
static int pshift,potshift;
static double rminshift;
#ifdef SSE2INTRINSIC
void printm128d(__m128d min,const char* str)
{
  double tmp[2];
  _mm_store_pd(tmp,min);
  mexPrintf("%s = [%e,%e]\n",str,tmp[0],tmp[1]);
}
#define PRINTM128D(in) printm128d(in,#in)
#endif
/*------------------------------------------------------------------------*/
void shift_alloc(int p,int pot,int maxm2p)
/* Allocates and initializes static data given a polynomial order p,
   type of potential pot and maximum number of cluster to cluster
   interactions maxm2p. */
{
  if (pshift != p) {
    // workspace arrays
    mxFree(wksp);
    wksp = (dcmplx *)mxMalloc(2*(p+1)*sizeof(dcmplx));
    #ifndef C_CODE
    mexMakeMemoryPersistent(wksp);
    #endif
    // powers of distances (used by m2ps)
    mxFree(rpow);
    rpow = (dcmplx *)mxMalloc((p+1)*sizeof(dcmplx));
    rpow[0] = 1.0;
    #ifndef C_CODE
    mexMakeMemoryPersistent(rpow);
    #endif
    pshift = p;
    potshift = pot;
    rminshift = exp2(-1000.0/(pshift+1));
    #ifndef C_CODE
    mexAtExit(shift_atExit);
    #endif
  }
}
/*------------------------------------------------------------------------*/
#ifndef C_CODE
void shift_atExit(void)
/* Deallocates persistent data. */
{
  mxFree(rpow);
  mxFree(wksp);
  pshift = 0;
}
#endif
/*------------------------------------------------------------------------*/
void shift_m2m(dcmplx r0,dcmplx r1,dcmplx r2,dcmplx r3,
	       dcmplx *This,const dcmplx *that)
/* Multipole to multipole shift from 4 children to 1 parent; This +=
   shift(that[0...p]) {child 0} + shift(that[p+1...2*p]) {child 1} +
   ...  and so on. */
{
  dcmplx rpow[4],mul,invr;
  rpow[0] = r0;
  rpow[1] = r1;
  rpow[2] = r2;
  rpow[3] = r3;
  for(int i = 0; i < 4; i++, that += pshift+1) {
    //if shift distance is too small, run safe code
    if(fabs(creal(rpow[i]))<rminshift&&fabs(cimag(rpow[i]))<rminshift) { //if overflow is possible, run safe code
      memcpy(wksp, that, (pshift+1)*sizeof(dcmplx));

      for(int k = pshift; k >= 2; k--)
        for(int j = k; j <= pshift; j++)
          wksp[j] += rpow[i]*wksp[j-1];

      if (potshift == 0) {
        mul = rpow[i];
        This[0] += that[0];
        for (int k = 1; k <= pshift; k++) {
          This[k] += wksp[k]-(that[0]/(double)k)*mul;
          mul *= rpow[i];
        }
      }
      else
        for(int k = 0; k <= pshift; k++)
          This[k] += wksp[k];
    }
    else { //otherwise, no overflow, run faster code
      mul=invr=1.0/rpow[i];
      wksp[0]=that[0];

      #ifdef SSE2INTRINSIC
      {
        double rsqr=creal(mul)*creal(mul)-cimag(mul)*cimag(mul);
        double isqr=2*creal(mul)*cimag(mul);
        __m128d mrtmp=_mm_set_pd(rsqr,creal(mul));
        __m128d mitmp=_mm_set_pd(isqr,cimag(mul));
        __m128d mrsqr=_mm_set1_pd(rsqr);
        __m128d misqr=_mm_set1_pd(isqr);
        __m128d mtmp=_mm_load_pd((double*)&that[1]);
        __m128d mtmp2=_mm_load_pd((double*)&that[2]);
        __m128d mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
        __m128d mithat=_mm_unpackhi_pd(mtmp,mtmp2);
        mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
        mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
        _mm_store_pd((double*)&wksp[1],_mm_unpacklo_pd(mtmp,mtmp2));
        _mm_store_pd((double*)&wksp[2],_mm_unpackhi_pd(mtmp,mtmp2));
        int j=3;
        for(;j<=pshift-1;j+=2) {
          mtmp=_mm_load_pd((double*)&that[j]);
          mtmp2=_mm_load_pd((double*)&that[j+1]);
          mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
          mithat=_mm_unpackhi_pd(mtmp,mtmp2);
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
          mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
          mrtmp=mtmp;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
          _mm_store_pd((double*)&wksp[j],_mm_unpacklo_pd(mtmp,mtmp2));
          _mm_store_pd((double*)&wksp[j+1],_mm_unpackhi_pd(mtmp,mtmp2));
        }
        if(j<=pshift) {
          mtmp=_mm_load_pd((double*)&that[j]);
          mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
          mithat=_mm_unpackhi_pd(mtmp,mtmp2);
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
          mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
          mrtmp=mtmp;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
          _mm_store_pd((double*)&wksp[j],_mm_unpacklo_pd(mtmp,mtmp2));
        }
      }
      #else
      wksp[1]=that[1]*mul; //pshift should never be smaller than 1
      for(int k=2;k<=pshift;k++) {
        mul*=invr;
        wksp[k]=that[k]*mul;
      }
      #endif
//    for(int k=pshift;k>=2;k--) {
//      for(int j=k;j<=pshift;j++) {
//        wksp[j]+=wksp[j-1];
//      }
//    }
      downwardshift((double*)wksp, 1,pshift);
      mul=rpow[i];
      This[0]+=that[0];
      if (potshift == 0) {
        for (int k = 1; k <= pshift; k++) {
          This[k] += (wksp[k] -that[0]*(1.0/k))*mul;
          mul*=rpow[i];
        }
      }
      else {
        #ifdef SSE2INTRINSIC
        {
          double rsqr=creal(mul)*creal(mul)-cimag(mul)*cimag(mul);
          double isqr=2*creal(mul)*cimag(mul);
          __m128d mrtmp=_mm_set_pd(rsqr,creal(mul));
          __m128d mitmp=_mm_set_pd(isqr,cimag(mul));
          __m128d mrsqr=_mm_set1_pd(rsqr);
          __m128d misqr=_mm_set1_pd(isqr);
          __m128d mtmp=_mm_load_pd((double*)&wksp[1]);
          __m128d mtmp2=_mm_load_pd((double*)&wksp[2]); //should not really be a problem even if pshift=1, as wksp should be larger than 2
          __m128d mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
          __m128d miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
          __m128d* mThis=(__m128d*)This;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrwksp),_mm_mul_pd(mitmp,miwksp));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,miwksp),_mm_mul_pd(mitmp,mrwksp));
          mThis[1]=_mm_add_pd(mThis[1],_mm_unpacklo_pd(mtmp,mtmp2));
          if(pshift>=2)
            mThis[2]=_mm_add_pd(mThis[2],_mm_unpackhi_pd(mtmp,mtmp2));

          int j=3;
          for(;j<=pshift-1;j+=2) {
            mtmp=_mm_load_pd((double*)&wksp[j]);
            mtmp2=_mm_load_pd((double*)&wksp[j+1]);
            mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
            miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
            mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
            mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
            mrtmp=mtmp;
            mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrwksp),_mm_mul_pd(mitmp,miwksp));
            mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,miwksp),_mm_mul_pd(mitmp,mrwksp));
            mThis[j]=_mm_add_pd(mThis[j],_mm_unpacklo_pd(mtmp,mtmp2));
            mThis[j+1]=_mm_add_pd(mThis[j+1],_mm_unpackhi_pd(mtmp,mtmp2));
          }
          if(j<=pshift) {
            mtmp=_mm_load_pd((double*)&wksp[j]);
            mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
            miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
            mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
            mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
            mrtmp=mtmp;
            mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrwksp),_mm_mul_pd(mitmp,miwksp));
            mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,miwksp),_mm_mul_pd(mitmp,mrwksp));
            mThis[j]=_mm_add_pd(mThis[j],_mm_unpacklo_pd(mtmp,mtmp2));
          }
        }
        #else
        This[1]+=wksp[1]*mul;
        for (int k = 2; k <= pshift; k++) {
          mul*=rpow[i];
          This[k] += wksp[k]*mul;
        }
        #endif
      }
    }
  }
}
/*------------------------------------------------------------------------*/
void shift_m2ps(dcmplx z0,const void *zi,size_t siz,
                const int *ix,int begin,int end,
                const dcmplx *This1,dcmplx *This2,
                const dcmplx *that1,dcmplx *that2)
/* Symmetric multipole to polynomial shift; This2 +=
   shift(that1[ix[i]]) followed by that2[ix[i]] += shift(This1) for
   index i in [begin,end). Center coordinates of the boxes are given
   in z0 and zi (in which each entry is separated by siz bytes).

   A somewhat singular feature is that a negative index in ix,
   pointing to a box at a previous level, is detected and implies that
   the 3 following indices are ignored. This ensures that connections
   between levels are properly handled and that the connectivity
   indices are convenient to assemble. Note: this interaction through
   parents as in the C99-code is not implemented by the main routines,
   but is supported here. */
{
  for (int i = begin; i < end; i++) {
#ifdef INLINECOMPLEX
    // typcast (char *) to avoid compiler warning:
    const dcmplx r0 = 1.0/(*(dcmplx *)((char *)zi+ix[i]*siz)-z0);
    const double re = creal(r0),im = cimag(r0);
    double *ri = (double *)rpow;
    double *wksp1 = (double *)&wksp[0],*wksp2 = (double *)&wksp[pshift+1];
    const dcmplx *that1i = &that1[ix[i]*(pshift+1)];
    double mre = 1.0,mim = 0.0,pre,pim;
    {

#ifdef SSE2M2PSSCALING
#ifdef _MSC_VER
    __m128d SIGNMASK=_mm_castsi128_pd(_mm_set_epi32(0x80000000,0,0x80000000,0));
#else
    __m128d SIGNMASK=_mm_castsi128_pd(_mm_set1_epi64((__m64)0x8000000000000000));
#endif
    double rsqr=re*re-im*im;
    double isqr=2*re*im;
    ri[2]=re;
    ri[4]=im;
    ri[3]=rsqr;
    ri[5]=isqr;
    __m128d mrtmp=_mm_set_pd(rsqr,re);
    __m128d mitmp=_mm_set_pd(isqr,im);
    __m128d mrsqr=_mm_set1_pd(rsqr);
    __m128d misqr=_mm_set1_pd(isqr);
    __m128d mtmp=_mm_load_pd((double*)&that1i[1]);
    __m128d mtmp2=_mm_load_pd((double*)&that1i[2]);
    __m128d mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
    __m128d mithat=_mm_unpackhi_pd(mtmp,mtmp2);
    mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
    mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
    _mm_store_pd(&wksp1[0],_mm_xor_pd(_mm_unpacklo_pd(mtmp,mtmp2),SIGNMASK));
    _mm_store_pd(&wksp1[2],_mm_unpackhi_pd(mtmp,mtmp2));
    mtmp=_mm_load_pd((double*)&This1[1]);
    mtmp2=_mm_load_pd((double*)&This1[2]);
    mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
    mithat=_mm_unpackhi_pd(mtmp,mtmp2);
    mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
    mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
    _mm_store_pd(&wksp2[0],_mm_unpacklo_pd(mtmp,mtmp2));
    _mm_store_pd(&wksp2[2],_mm_unpackhi_pd(mtmp,mtmp2));
    #endif
    // prescale
    #ifdef SSE2M2PSSCALING
    for (int j = 3; j < pshift; )
    #else
    for (int j = 1; j < pshift; )
    #endif
    {
      #ifdef SSE2M2PSSCALING
      mtmp=_mm_load_pd((double*)&that1i[j]);
      mtmp2=_mm_load_pd((double*)&that1i[j+1]);
      mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
      mithat=_mm_unpackhi_pd(mtmp,mtmp2);
      mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
      mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
      mrtmp=mtmp;
      _mm_store_pd(&ri[2*j],mrtmp);
      _mm_store_pd(&ri[2*j+2],mitmp);
      mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
      mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
      _mm_store_pd(&wksp1[2*j-2],_mm_xor_pd(_mm_unpacklo_pd(mtmp,mtmp2),SIGNMASK));
      _mm_store_pd(&wksp1[2*j],_mm_unpackhi_pd(mtmp,mtmp2));
      mtmp=_mm_load_pd((double*)&This1[j]);
      mtmp2=_mm_load_pd((double*)&This1[j+1]);
      mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
      mithat=_mm_unpackhi_pd(mtmp,mtmp2);
      mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
      mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
      _mm_store_pd(&wksp2[2*j-2],_mm_unpacklo_pd(mtmp,mtmp2));
      _mm_store_pd(&wksp2[2*j],_mm_unpackhi_pd(mtmp,mtmp2));
      j+=2;
      #else
      ri[2*j] = mre*re-mim*im;
      ri[2*j+1] = mre*im+mim*re;
      mre = ri[2*j];
      mim = ri[2*j+1];
      pre = creal(that1i[j])*mre-cimag(that1i[j])*mim;
      pim = creal(that1i[j])*mim+cimag(that1i[j])*mre;
      wksp1[2*j-2] = -pre;
      wksp1[2*j-1] = -pim;
      pre = creal(This1[j])*mre-cimag(This1[j])*mim;
      pim = creal(This1[j])*mim+cimag(This1[j])*mre;
      wksp2[2*j-2] = pre;
      wksp2[2*j-1] = pim;
      j++;

      ri[2*j] = mre*re-mim*im;
      ri[2*j+1] = mre*im+mim*re;
      mre = ri[2*j];
      mim = ri[2*j+1];
      pre = creal(that1i[j])*mre-cimag(that1i[j])*mim;
      pim = creal(that1i[j])*mim+cimag(that1i[j])*mre;
      wksp1[2*j-2] = pre;
      wksp1[2*j-1] = pim;
      pre = creal(This1[j])*mre-cimag(This1[j])*mim;
      pim = creal(This1[j])*mim+cimag(This1[j])*mre;
      wksp2[2*j-2] = pre;
      wksp2[2*j-1] = pim;
      j++;
      #endif
    }
    if (pshift&1) {
      #ifdef SSE2M2PSSCALING
      mre=ri[2*pshift-3];
      mim=ri[2*pshift-1];
      #endif
      ri[2*pshift] = mre*re-mim*im;
      ri[2*pshift+1] = mre*im+mim*re;
      mre = ri[2*pshift];
      mim = ri[2*pshift+1];
      pre = creal(that1i[pshift])*mre-cimag(that1i[pshift])*mim;
      pim = creal(that1i[pshift])*mim+cimag(that1i[pshift])*mre;
      wksp1[2*pshift-2] = -pre;
      wksp1[2*pshift-1] = -pim;
      pre = creal(This1[pshift])*mre-cimag(This1[pshift])*mim;
      pim = creal(This1[pshift])*mim+cimag(This1[pshift])*mre;
      wksp2[2*pshift-2] = pre;
      wksp2[2*pshift-1] = pim;
    }
    }
    wksp1[2*pshift] = wksp1[2*pshift+1] = 0.0;
    wksp2[2*pshift] = wksp2[2*pshift+1] = 0.0;

    // (1) "upward and downward" phase
    upwardshift2(wksp1,wksp2,pshift);
    downwardshift2(wksp1,wksp2,pshift);


    // postscale
    if (potshift == 0) {
      dcmplx *that2i = &that2[ix[i]*(pshift+1)];
      const double rec0 = creal(This1[0]),imc0 = cimag(This1[0]);
      const double red0 = creal(that1[ix[i]*(pshift+1)]),
        imd0 = cimag(that1[ix[i]*(pshift+1)]);

      {
#ifdef SSE2M2PSSCALING
        const dcmplx logr = clog(ri[2]+I*ri[4]);
#else
        const dcmplx logr = clog(ri[2]+I*ri[3]);
#endif
        This2[0] += wksp1[0]+I*wksp1[1]
          -(red0*creal(logr)-imd0*cimag(logr)+
            I*(imd0*creal(logr)+red0*cimag(logr)));
        that2i[0] += wksp2[0]+I*wksp2[1]
          -(rec0*creal(logr)-imc0*(cimag(logr)-M_PI)+
            I*(imc0*creal(logr)+rec0*(cimag(logr)-M_PI)));
      }
      for (int j = 1; j < pshift; ) {
#ifdef SSE2M2PSSCALING
        pre = (wksp1[2*j]-red0/j)*ri[2*j]-
          (wksp1[2*j+1]-imd0/j)*ri[2*j+2];
        pim = (wksp1[2*j]-red0/j)*ri[2*j+2]+
          (wksp1[2*j+1]-imd0/j)*ri[2*j];
        COMPLEXADD(This2[j],pre,pim); //does not have any significant effect for log pot compared to complex version above (macro is much faster than built in operator+)
        pre = (wksp2[2*j]-rec0/j)*ri[2*j]-
          (wksp2[2*j+1]-imc0/j)*ri[2*j+2];
        pim = (wksp2[2*j]-rec0/j)*ri[2*j+2]+
          (wksp2[2*j+1]-imc0/j)*ri[2*j];
        COMPLEXSUB(that2i[j],pre,pim);
        j++;

        pre = (wksp1[2*j]-red0/j)*ri[2*j-1]-
          (wksp1[2*j+1]-imd0/j)*ri[2*j+1];
        pim = (wksp1[2*j]-red0/j)*ri[2*j+1]+
          (wksp1[2*j+1]-imd0/j)*ri[2*j-1];
        COMPLEXADD(This2[j],pre,pim);
        pre = (wksp2[2*j]-rec0/j)*ri[2*j-1]-
          (wksp2[2*j+1]-imc0/j)*ri[2*j+1];
        pim = (wksp2[2*j]-rec0/j)*ri[2*j+1]+
          (wksp2[2*j+1]-imc0/j)*ri[2*j-1];
        COMPLEXADD(that2i[j],pre,pim);
        j++;
#else
        pre = (wksp1[2*j]-red0/j)*ri[2*j]-
          (wksp1[2*j+1]-imd0/j)*ri[2*j+1];
        pim = (wksp1[2*j]-red0/j)*ri[2*j+1]+
          (wksp1[2*j+1]-imd0/j)*ri[2*j];
        COMPLEXADD(This2[j],pre,pim);
        pre = (wksp2[2*j]-rec0/j)*ri[2*j]-
          (wksp2[2*j+1]-imc0/j)*ri[2*j+1];
        pim = (wksp2[2*j]-rec0/j)*ri[2*j+1]+
          (wksp2[2*j+1]-imc0/j)*ri[2*j];
        COMPLEXSUB(that2i[j],pre,pim);
        j++;

        pre = (wksp1[2*j]-red0/j)*ri[2*j]-
          (wksp1[2*j+1]-imd0/j)*ri[2*j+1];
        pim = (wksp1[2*j]-red0/j)*ri[2*j+1]+
          (wksp1[2*j+1]-imd0/j)*ri[2*j];
        COMPLEXADD(This2[j],pre,pim);
        pre = (wksp2[2*j]-rec0/j)*ri[2*j]-
          (wksp2[2*j+1]-imc0/j)*ri[2*j+1];
        pim = (wksp2[2*j]-rec0/j)*ri[2*j+1]+
          (wksp2[2*j+1]-imc0/j)*ri[2*j];
        COMPLEXADD(that2i[j],pre,pim);
        j++;
#endif
      }
      if (pshift&1) {
        pre = (wksp1[2*pshift]-red0/pshift)*ri[2*pshift]-
          (wksp1[2*pshift+1]-imd0/pshift)*ri[2*pshift+1];
        pim = (wksp1[2*pshift]-red0/pshift)*ri[2*pshift+1]+
          (wksp1[2*pshift+1]-imd0/pshift)*ri[2*pshift];
        COMPLEXADD(This2[pshift],pre,pim);
        pre = (wksp2[2*pshift]-rec0/pshift)*ri[2*pshift]-
          (wksp2[2*pshift+1]-imc0/pshift)*ri[2*pshift+1];
        pim = (wksp2[2*pshift]-rec0/pshift)*ri[2*pshift+1]+
          (wksp2[2*pshift+1]-imc0/pshift)*ri[2*pshift];
        COMPLEXSUB(that2i[pshift],pre,pim);
      }
    }
    else {
      dcmplx *that2i = &that2[ix[i]*(pshift+1)];

#ifdef SSE2M2PSSCALING
      __m128d* mThis=(__m128d*)This2;
      __m128d* mthat2i=(__m128d*)that2i;
      mThis[0]=_mm_add_pd(mThis[0],((__m128d*)wksp1)[0]);
      mthat2i[0]=_mm_add_pd(mthat2i[0],((__m128d*)wksp2)[0]);
      for (int j = 1; j < pshift; ) { //same loop indices as above, to enable different storage of ri than in the other code
        __m128d mrtmp=_mm_load_pd(&ri[2*j]);
        __m128d mitmp=_mm_load_pd(&ri[2*j+2]);
        __m128d mtmp=_mm_load_pd(&wksp1[2*j]);
        __m128d mtmp2=_mm_load_pd(&wksp1[2*j+2]);
        __m128d mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
        __m128d miwksp=_mm_unpackhi_pd(mtmp,mtmp2);

        mtmp=_mm_sub_pd(_mm_mul_pd(mrwksp,mrtmp),_mm_mul_pd(miwksp,mitmp));
        mtmp2=_mm_add_pd(_mm_mul_pd(mrwksp,mitmp),_mm_mul_pd(miwksp,mrtmp));
        mThis[j]=_mm_add_pd(mThis[j],_mm_unpacklo_pd(mtmp,mtmp2));
        mThis[j+1]=_mm_add_pd(mThis[j+1],_mm_unpackhi_pd(mtmp,mtmp2));
        mtmp=_mm_load_pd(&wksp2[2*j]);
        mtmp2=_mm_load_pd(&wksp2[2*j+2]);
        mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
        miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
        mtmp=_mm_sub_pd(_mm_mul_pd(mrwksp,mrtmp),_mm_mul_pd(miwksp,mitmp));
        mtmp2=_mm_add_pd(_mm_mul_pd(mrwksp,mitmp),_mm_mul_pd(miwksp,mrtmp));
        mthat2i[j]=_mm_sub_pd(mthat2i[j],_mm_unpacklo_pd(mtmp,mtmp2));
        mthat2i[j+1]=_mm_add_pd(mthat2i[j+1],_mm_unpackhi_pd(mtmp,mtmp2));
        j+=2;
      }
      if (pshift&1) {
        pre = wksp1[2*pshift]*ri[2*pshift]-wksp1[2*pshift+1]*ri[2*pshift+1];
        pim = wksp1[2*pshift]*ri[2*pshift+1]+wksp1[2*pshift+1]*ri[2*pshift];
        COMPLEXADD(This2[pshift],pre,pim);
        pre = wksp2[2*pshift]*ri[2*pshift]-wksp2[2*pshift+1]*ri[2*pshift+1];
        pim = wksp2[2*pshift]*ri[2*pshift+1]+wksp2[2*pshift+1]*ri[2*pshift];
        COMPLEXSUB(that2i[pshift],pre,pim);
      }
#else
      for (int j = 0; j < pshift; ) {
        pre = wksp1[2*j]*ri[2*j]-wksp1[2*j+1]*ri[2*j+1];
        pim = wksp1[2*j]*ri[2*j+1]+wksp1[2*j+1]*ri[2*j];
        COMPLEXADD(This2[j],pre,pim); //significant speed improvement compared to using operator+
        pre = wksp2[2*j]*ri[2*j]-wksp2[2*j+1]*ri[2*j+1];
        pim = wksp2[2*j]*ri[2*j+1]+wksp2[2*j+1]*ri[2*j];
        COMPLEXADD(that2i[j],pre,pim);
        j++;

        pre = wksp1[2*j]*ri[2*j]-wksp1[2*j+1]*ri[2*j+1];
        pim = wksp1[2*j]*ri[2*j+1]+wksp1[2*j+1]*ri[2*j];
        COMPLEXADD(This2[j],pre,pim);
        pre = wksp2[2*j]*ri[2*j]-wksp2[2*j+1]*ri[2*j+1];
        pim = wksp2[2*j]*ri[2*j+1]+wksp2[2*j+1]*ri[2*j];
        COMPLEXSUB(that2i[j],pre,pim);
        j++;
      }
      if (~pshift&1) {
        pre = wksp1[2*pshift]*ri[2*pshift]-wksp1[2*pshift+1]*ri[2*pshift+1];
        pim = wksp1[2*pshift]*ri[2*pshift+1]+wksp1[2*pshift+1]*ri[2*pshift];
        COMPLEXADD(This2[pshift],pre,pim);
        pre = wksp2[2*pshift]*ri[2*pshift]-wksp2[2*pshift+1]*ri[2*pshift+1];
        pim = wksp2[2*pshift]*ri[2*pshift+1]+wksp2[2*pshift+1]*ri[2*pshift];
        COMPLEXADD(that2i[pshift],pre,pim);
      }
#endif
    }
#else // !INLINECOMPLEX
    // typcast (char *) to avoid compiler warning:
    const dcmplx r0 = 1.0/(*(dcmplx *)((char *)zi+ix[i]*siz)-z0); //
    dcmplx *wksp1 = &wksp[0],*wksp2 = &wksp[pshift+1];
    const dcmplx *that1i = &that1[ix[i]*(pshift+1)];
    dcmplx mul = 1.0;

    // prescale
    for (int j = 1; j < pshift; ) {
      mul *= r0;
      rpow[j] = mul;
      wksp1[j-1] = -that1i[j]*mul;
      wksp2[j-1] = This1[j]*mul;
      j++;
      mul *= r0;
      rpow[j] = mul;
      wksp1[j-1] = that1i[j]*mul;
      wksp2[j-1] = This1[j]*mul;
      j++;
    }
    if (pshift&1) {
      mul *= r0;
      rpow[pshift] = mul;
      wksp1[pshift-1] = -that1i[pshift]*mul;
      wksp2[pshift-1] = This1[pshift]*mul;
    }
    wksp1[pshift] = wksp2[pshift] = 0.0;

    // (1) "upward" phase
    for (int k = 2; k <= pshift; k++)
      for (int j = pshift-k; j < pshift; j++) {
        wksp1[j] += wksp1[j+1];
        wksp2[j] += wksp2[j+1];
      }

    // (2) "downward" phase
    for (int k = pshift; k > 0; k--)
      for (int j = k; j <= pshift; j++) {
        wksp1[j] += wksp1[j-1];
        wksp2[j] += wksp2[j-1];
      }

    // postscale
    if (potshift == 0) {
      dcmplx *that2i = &that2[ix[i]*(pshift+1)];
      const dcmplx This0 = This1[0];
      const dcmplx that0 = that1[ix[i]*(pshift+1)];

      This2[0] += wksp1[0]-that0*clog(rpow[1]);
      that2i[0] += wksp2[0]-This0*clog(-rpow[1]); // *** can avoid the 2nd call
      for (int j = 1; j < pshift; ) {
        This2[j] += (wksp1[j]-that0/(double)j)*rpow[j];
        that2i[j] -= (wksp2[j]-This0/(double)j)*rpow[j];
        j++;
        This2[j] += (wksp1[j]-that0/(double)j)*rpow[j];
        that2i[j] += (wksp2[j]-This0/(double)j)*rpow[j];
        j++;
      }
      if (pshift&1) {
        This2[pshift] += (wksp1[pshift]-that0/(double)pshift)*rpow[pshift];
        that2i[pshift] -= (wksp2[pshift]-This0/(double)pshift)*rpow[pshift];
      }
    }
    else {
      dcmplx *that2i = &that2[ix[i]*(pshift+1)];

      for (int j = 0; j < pshift; ) {
        This2[j] += wksp1[j]*rpow[j];
        that2i[j] += wksp2[j]*rpow[j];
        j++;
        This2[j] += wksp1[j]*rpow[j];
        that2i[j] -= wksp2[j]*rpow[j];
        j++;
      }
      if (~pshift&1) {
        This2[pshift] += wksp1[pshift]*rpow[pshift];
        that2i[pshift] += wksp2[pshift]*rpow[pshift];
      }
    }
#endif // INLINECOMPLEX

    // in case the shift was managed via parents
    if (ix[i] < 0) i += 3;
  }
}

/*------------------------------------------------------------------------*/
void shift_p2p(dcmplx r0,dcmplx r1,dcmplx r2,dcmplx r3,
	       dcmplx *This,const dcmplx *that)
{
  dcmplx tmp,invr;
  dcmplx rpow[4];
  rpow[0] = r0;
  rpow[1] = r1;
  rpow[2] = r2;
  rpow[3] = r3;
  for(int i = 0; i < 4; i++, This += pshift+1) {
    if(fabs(creal(rpow[i]))<rminshift&&fabs(cimag(rpow[i]))<rminshift) {
      memcpy(wksp, that, (pshift+1)*sizeof(dcmplx));

      for (int k = 1; k <= pshift; k++)
        for(int j = pshift-k; j < pshift; j++)
          wksp[j] -= rpow[i]*wksp[j+1];

      for(int k = 0; k <= pshift; k++)
        This[k] += wksp[k];
    }
    else {
      rpow[i]=-rpow[i];
      #ifndef SSE2INTRINSIC
      #ifndef INLINECOMPLEX
      tmp=rpow[i];
      wksp[0]=that[0]*tmp;
      for(int k=1;k<=pshift;k++) {
        tmp*=rpow[i];
        wksp[k]=that[k]*tmp;
      }
      #else
      double rtmp=creal(rpow[i]);
      double itmp=cimag(rpow[i]);
      COMPLEXASSIGN(wksp[0],creal(that[0])*rtmp-cimag(that[0])*itmp,creal(that[0])*itmp+cimag(that[0])*rtmp);
      for(int k=1;k<=pshift;k++) {
        double tmp=rtmp*creal(rpow[i])-itmp*cimag(rpow[i]);
        itmp=itmp*creal(rpow[i])+rtmp*cimag(rpow[i]);
        rtmp=tmp;
        COMPLEXASSIGN(wksp[k],creal(that[k])*rtmp-cimag(that[k])*itmp,creal(that[k])*itmp+cimag(that[k])*rtmp);
      }
      #endif
      #else /*SSE2INTRINSIC*/
      {
        double rsqr=creal(rpow[i])*creal(rpow[i])-cimag(rpow[i])*cimag(rpow[i]);
        double isqr=2*creal(rpow[i])*cimag(rpow[i]);
        __m128d mrtmp=_mm_set_pd(rsqr,creal(rpow[i]));
        __m128d mitmp=_mm_set_pd(isqr,cimag(rpow[i]));
        __m128d mrsqr=_mm_set1_pd(rsqr);
        __m128d misqr=_mm_set1_pd(isqr);
        __m128d mtmp=_mm_load_pd((double*)&that[0]);
        __m128d mtmp2=_mm_load_pd((double*)&that[1]);
        __m128d mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
        __m128d mithat=_mm_unpackhi_pd(mtmp,mtmp2);
        mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
        mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
        _mm_store_pd((double*)&wksp[0],_mm_unpacklo_pd(mtmp,mtmp2));
        _mm_store_pd((double*)&wksp[1],_mm_unpackhi_pd(mtmp,mtmp2));
        int j=2;
        for(;j<=pshift-1;j+=2) {
          mtmp=_mm_load_pd((double*)&that[j]);
          mtmp2=_mm_load_pd((double*)&that[j+1]);
          mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
          mithat=_mm_unpackhi_pd(mtmp,mtmp2);
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
          mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
          mrtmp=mtmp;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
          _mm_store_pd((double*)&wksp[j],_mm_unpacklo_pd(mtmp,mtmp2));
          _mm_store_pd((double*)&wksp[j+1],_mm_unpackhi_pd(mtmp,mtmp2));
        }
        if(j<=pshift) {
          mtmp=_mm_load_pd((double*)&that[j]);
          mrthat=_mm_unpacklo_pd(mtmp,mtmp2);
          mithat=_mm_unpackhi_pd(mtmp,mtmp2);
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
          mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
          mrtmp=mtmp;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrthat),_mm_mul_pd(mitmp,mithat));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,mithat),_mm_mul_pd(mitmp,mrthat));
          _mm_store_pd((double*)&wksp[j],_mm_unpacklo_pd(mtmp,mtmp2));
        }
      }
      #endif
//      for(int j=0;j<=pshift;j++) {
//        for(int k=pshift-j;k<pshift;k++) {
//          wksp[k]+=wksp[k+1];
//        }
//      }
      upwardshift((double*)wksp,1,pshift);
      tmp=invr=1.0/rpow[i];
      #ifdef SSE2INTRINSIC
      {
        double rsqr=creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp);
        double isqr=2*creal(tmp)*cimag(tmp);
        __m128d mrtmp=_mm_set_pd(rsqr,creal(tmp));
        __m128d mitmp=_mm_set_pd(isqr,cimag(tmp));
        __m128d mrsqr=_mm_set1_pd(rsqr);
        __m128d misqr=_mm_set1_pd(isqr);
        __m128d mtmp=_mm_load_pd((double*)&wksp[0]);
        __m128d mtmp2=_mm_load_pd((double*)&wksp[1]);
        __m128d mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
        __m128d miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
        __m128d* mThis=(__m128d*)This;
        mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrwksp),_mm_mul_pd(mitmp,miwksp));
        mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,miwksp),_mm_mul_pd(mitmp,mrwksp));
        mThis[0]=_mm_add_pd(mThis[0],_mm_unpacklo_pd(mtmp,mtmp2));
        mThis[1]=_mm_add_pd(mThis[1],_mm_unpackhi_pd(mtmp,mtmp2));

        int j=2;
        for(;j<=pshift-1;j+=2) {
          mtmp=_mm_load_pd((double*)&wksp[j]);
          mtmp2=_mm_load_pd((double*)&wksp[j+1]);
          mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
          miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
          mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
          mrtmp=mtmp;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrwksp),_mm_mul_pd(mitmp,miwksp));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,miwksp),_mm_mul_pd(mitmp,mrwksp));
          mThis[j]=_mm_add_pd(mThis[j],_mm_unpacklo_pd(mtmp,mtmp2));
          mThis[j+1]=_mm_add_pd(mThis[j+1],_mm_unpackhi_pd(mtmp,mtmp2));
        }
        if(j<=pshift) {
          mtmp=_mm_load_pd((double*)&wksp[j]);
          mrwksp=_mm_unpacklo_pd(mtmp,mtmp2);
          miwksp=_mm_unpackhi_pd(mtmp,mtmp2);
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrsqr),_mm_mul_pd(mitmp,misqr));
          mitmp=_mm_add_pd(_mm_mul_pd(mrtmp,misqr),_mm_mul_pd(mitmp,mrsqr));
          mrtmp=mtmp;
          mtmp=_mm_sub_pd(_mm_mul_pd(mrtmp,mrwksp),_mm_mul_pd(mitmp,miwksp));
          mtmp2=_mm_add_pd(_mm_mul_pd(mrtmp,miwksp),_mm_mul_pd(mitmp,mrwksp));
          mThis[j]=_mm_add_pd(mThis[j],_mm_unpacklo_pd(mtmp,mtmp2));
        }
      }
      #else
      This[0]+=wksp[0]*tmp;
      for(int k=1;k<=pshift;k++) {
        tmp*=invr;
        This[k]+=wksp[k]*tmp;
      }
      #endif
    }
  }
}
/*------------------------------------------------------------------------*/
