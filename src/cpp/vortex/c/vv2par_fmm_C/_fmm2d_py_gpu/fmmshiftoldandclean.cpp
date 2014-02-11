/* fmmshift.cpp */

/* S. Engblom 2010-02-18 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#include <string.h>
#include "fmmshift.h"

// deallocate persistent data
void shift_atExit(void);

// code responds to #define INLINECOMPLEX

// workspace arrays
static dcmplx *wksp,*rpow;

// scalars defining shift-type and sizes
static int pshift,potshift;

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
    mexMakeMemoryPersistent(wksp);

    // powers of distances (used by m2ps)
    mxFree(rpow);
    rpow = (dcmplx *)mxMalloc((p+1)*sizeof(dcmplx));
    rpow[0] = 1.0;
    mexMakeMemoryPersistent(rpow);

    pshift = p;
    potshift = pot;
    mexAtExit(shift_atExit);
  }
}
/*------------------------------------------------------------------------*/
void shift_atExit(void)
/* Deallocates persistent data. */
{
  mxFree(rpow);
  mxFree(wksp);
  pshift = 0;
}
/*------------------------------------------------------------------------*/
void shift_m2m(dcmplx r0,dcmplx r1,dcmplx r2,dcmplx r3,
               dcmplx *This,const dcmplx *that)
/* Multipole to multipole shift from 4 children to 1 parent; This +=
   shift(that[0...p]) {child 0} + shift(that[p+1...2*p]) {child 1} +
   ...  and so on. */
{
  dcmplx rpow[4],mul;
  rpow[0] = r0;
  rpow[1] = r1;
  rpow[2] = r2;
  rpow[3] = r3;

  for(int i = 0; i < 4; i++, that += pshift+1) {
    memcpy(wksp,that,(pshift+1)*sizeof(dcmplx));

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

    // prescale
    for (int j = 1; j < pshift; ) {
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
    }
    if (pshift&1) {
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
    wksp1[2*pshift] = wksp1[2*pshift+1] = 0.0;
    wksp2[2*pshift] = wksp2[2*pshift+1] = 0.0;

    // (1) "upward" phase
    for (int k = 2; k <= pshift; k++)
      for (int j = pshift-k; j < pshift; j++) {
        wksp1[2*j] += wksp1[2*(j+1)];
        wksp1[2*j+1] += wksp1[2*(j+1)+1];
        wksp2[2*j] += wksp2[2*(j+1)];
        wksp2[2*j+1] += wksp2[2*(j+1)+1];
      }

    // (2) "downward" phase
    for (int k = pshift; k > 0; k--)
      for (int j = k; j <= pshift; j++) {
        wksp1[2*j] += wksp1[2*(j-1)];
        wksp1[2*j+1] += wksp1[2*(j-1)+1];
        wksp2[2*j] += wksp2[2*(j-1)];
        wksp2[2*j+1] += wksp2[2*(j-1)+1];
      }

    // postscale
    if (potshift == 0) {
      dcmplx *that2i = &that2[ix[i]*(pshift+1)];
      const double rec0 = creal(This1[0]),imc0 = cimag(This1[0]);
      const double red0 = creal(that1[ix[i]*(pshift+1)]),
        imd0 = cimag(that1[ix[i]*(pshift+1)]);

      {
        const dcmplx logr = clog(ri[2]+I*ri[3]);
        This2[0] += wksp1[0]+I*wksp1[1]
          -(red0*creal(logr)-imd0*cimag(logr)+
            I*(imd0*creal(logr)+red0*cimag(logr)));
        that2i[0] += wksp2[0]+I*wksp2[1]
          -(rec0*creal(logr)-imc0*(cimag(logr)-M_PI)+
            I*(imc0*creal(logr)+rec0*(cimag(logr)-M_PI)));
      }
      for (int j = 1; j < pshift; ) {
        pre = (wksp1[2*j]-red0/j)*ri[2*j]-
          (wksp1[2*j+1]-imd0/j)*ri[2*j+1];
        pim = (wksp1[2*j]-red0/j)*ri[2*j+1]+
          (wksp1[2*j+1]-imd0/j)*ri[2*j];
        This2[j] += pre+I*pim;
        pre = (wksp2[2*j]-rec0/j)*ri[2*j]-
          (wksp2[2*j+1]-imc0/j)*ri[2*j+1];
        pim = (wksp2[2*j]-rec0/j)*ri[2*j+1]+
          (wksp2[2*j+1]-imc0/j)*ri[2*j];
        that2i[j] -= pre+I*pim;
        j++;

        pre = (wksp1[2*j]-red0/j)*ri[2*j]-
          (wksp1[2*j+1]-imd0/j)*ri[2*j+1];
        pim = (wksp1[2*j]-red0/j)*ri[2*j+1]+
          (wksp1[2*j+1]-imd0/j)*ri[2*j];
        This2[j] += pre+I*pim;
        pre = (wksp2[2*j]-rec0/j)*ri[2*j]-
          (wksp2[2*j+1]-imc0/j)*ri[2*j+1];
        pim = (wksp2[2*j]-rec0/j)*ri[2*j+1]+
          (wksp2[2*j+1]-imc0/j)*ri[2*j];
        that2i[j] += pre+I*pim;
        j++;
      }
      if (pshift&1) {
        pre = (wksp1[2*pshift]-red0/pshift)*ri[2*pshift]-
          (wksp1[2*pshift+1]-imd0/pshift)*ri[2*pshift+1];
        pim = (wksp1[2*pshift]-red0/pshift)*ri[2*pshift+1]+
          (wksp1[2*pshift+1]-imd0/pshift)*ri[2*pshift];
        This2[pshift] += pre+I*pim;
        pre = (wksp2[2*pshift]-rec0/pshift)*ri[2*pshift]-
          (wksp2[2*pshift+1]-imc0/pshift)*ri[2*pshift+1];
        pim = (wksp2[2*pshift]-rec0/pshift)*ri[2*pshift+1]+
          (wksp2[2*pshift+1]-imc0/pshift)*ri[2*pshift];
        that2i[pshift] -= pre+I*pim;
      }
    }
    else {
      dcmplx *that2i = &that2[ix[i]*(pshift+1)];

      for (int j = 0; j < pshift; ) {
        pre = wksp1[2*j]*ri[2*j]-wksp1[2*j+1]*ri[2*j+1];
        pim = wksp1[2*j]*ri[2*j+1]+wksp1[2*j+1]*ri[2*j];
        This2[j] += pre+I*pim;
        pre = wksp2[2*j]*ri[2*j]-wksp2[2*j+1]*ri[2*j+1];
        pim = wksp2[2*j]*ri[2*j+1]+wksp2[2*j+1]*ri[2*j];
        that2i[j] += pre+I*pim;
        j++;

        pre = wksp1[2*j]*ri[2*j]-wksp1[2*j+1]*ri[2*j+1];
        pim = wksp1[2*j]*ri[2*j+1]+wksp1[2*j+1]*ri[2*j];
        This2[j] += pre+I*pim;
        pre = wksp2[2*j]*ri[2*j]-wksp2[2*j+1]*ri[2*j+1];
        pim = wksp2[2*j]*ri[2*j+1]+wksp2[2*j+1]*ri[2*j];
        that2i[j] -= pre+I*pim;
        j++;
      }
      if (~pshift&1) {
        pre = wksp1[2*pshift]*ri[2*pshift]-wksp1[2*pshift+1]*ri[2*pshift+1];
        pim = wksp1[2*pshift]*ri[2*pshift+1]+wksp1[2*pshift+1]*ri[2*pshift];
        This2[pshift] += pre+I*pim;
        pre = wksp2[2*pshift]*ri[2*pshift]-wksp2[2*pshift+1]*ri[2*pshift+1];
        pim = wksp2[2*pshift]*ri[2*pshift+1]+wksp2[2*pshift+1]*ri[2*pshift];
        that2i[pshift] += pre+I*pim;
      }
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
/* Polynomial to polynomial shift from 1 parent to 4 children;
   This[0..p] += shift(that) {child 0}, This[p+1..2*p] += shift(that)
   {child 1}, and so on. */
{
  dcmplx rpow[4];
  rpow[0] = r0;
  rpow[1] = r1;
  rpow[2] = r2;
  rpow[3] = r3;

  for(int i = 0; i < 4; i++, This += pshift+1) {
    memcpy(wksp,that,(pshift+1)*sizeof(dcmplx));

    for (int k = 1; k <= pshift; k++)
      for(int j = pshift-k; j < pshift; j++)
        wksp[j] -= rpow[i]*wksp[j+1];

    for(int k = 0; k <= pshift; k++)
      This[k] += wksp[k];
  }
}
/*------------------------------------------------------------------------*/
