/* direct.cpp */

/* S. Engblom 2010-02-18 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#include <stdio.h>
#include "expint.h"
#include "direct.h"

#include <emmintrin.h>
#include <string.h>

// code responds to #define INLINECOMPLEX_DIRECT
void printm128(__m128d in,char* intro) {
  double a[2];
  _mm_store_pd (a, in);
  mexPrintf("%s = %e + %ei\n",intro,a[0],a[1]);
}
#define PRINTM128(a) printm128(a,#a)
void rlog2(int N,
           const double *zr,const double *zi,
           const double *mr,
           double *pr,double *pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont);
void zlog2(int N,
           const double *zr,const double *zi,
           const double *mr,const double *mi,
           double *pr,double *pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont);
void rlog(int N,
          const double *zr,const double *zi,
          const double *mr,
          int NE,
          const double *er,const double *ei,
          double *qr,double *qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont);
void zlog(int N,
          const double *zr,const double *zi,
          const double *mr,const double *mi,
          int NE,
          const double *er,const double *ei,
          double *qr,double *qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont);

void rinv2(int N,
           const double *zr,const double *zi,
           const double *mr,
           double *pr,double *pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont);
void zinv2(int N,
           const double *zr,const double *zi,
           const double *mr,const double *mi,
           double *pr,double *pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont);
void rinv(int N,
          const double *zr,const double *zi,
          const double *mr,
          int NE,
          const double *er,const double *ei,
          double *qr,double *qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont);
void zinv(int N,
          const double *zr,const double *zi,
          const double *mr,const double *mi,
          int NE,
          const double *er,const double *ei,
          double *qr,double *qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont);

/*------------------------------------------------------------------------*/
void directInteract(int N,
                    const double *zr,const double *zi,
                    const double *mr,const double *mi,
                    int NE,
                    const double *er,const double *ei,
                    double *pr,double *pi,
                    double *qr,double *qi,
                    const panel *panels,int Npanel,int pot,
                    SMOOTHER smooth,double xopt,double cutoff,
                    bool cont,double* timing,int printtime)
/* Direct evaluation of the potential from N pointmasses (mr,mi) (with
   mi possibly NULL) at the positions (zr,zi). The evaluation is done
   at the points (zr,zi) if (pr,pi) is not NULL and/or at the points
   (er,ei) if (qr,qi) is not NULL. The result is stored in the vectors
   (pr,pi) and (qr,qi), respectively, which must be allocated and
   cleared prior to call. Logarithmic potential for pot == 0, harmonic
   potential otherwise. The smoothing is done using smoother smooth
   with parameters (xopt,cutoff,cont). */
{
  // evaluation at (zr,zi)
#ifdef SSE2DIRECT
  //SSE2 requires 16 byte alignment for proper memory loads. This is usually the case, except for some very small arrays
  //This is not included in the timing, as a normal code would use aligned vectors, but this cannot be controlled if the input comes from matlab
  double *zrtmp=NULL,*zitmp=NULL,*mrtmp=NULL,*mitmp=NULL,*ertmp=NULL,
         *eitmp=NULL,*prtmp=NULL,*pitmp=NULL,*qrtmp=NULL,*qitmp=NULL,
         *prold=NULL,*piold=NULL,*qrold=NULL,*qiold=NULL;
  if(zr!=NULL && ((size_t)zr&15)!=0) {
    zrtmp=(double*)mxMalloc((N+3)*sizeof(double));
    char *tmp=(char*)zrtmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    memcpy(tmp,zr,N*sizeof(double));
    zr=(double*)tmp;
  }
  if(zi!=NULL && ((size_t)zi&15)!=0) {
    zitmp=(double*)mxMalloc((N+3)*sizeof(double));
    char *tmp=(char*)zitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    memcpy(tmp,zi,N*sizeof(double));
    zi=(double*)tmp;
  }
  if(mr!=NULL && ((size_t)mr&15)!=0) {
    mrtmp=(double*)mxMalloc((N+3)*sizeof(double));
    char *tmp=(char*)mrtmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    memcpy(tmp,mr,N*sizeof(double));
    mr=(double*)tmp;
  }
  if(mi!=NULL && ((size_t)mi&15)!=0) {
    mitmp=(double*)mxMalloc((N+3)*sizeof(double));
    char *tmp=(char*)mitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    memcpy(tmp,mi,N*sizeof(double));
    mi=(double*)tmp;
  }
  if(pr!=NULL && ((size_t)pr&15)!=0) {
    prtmp=(double*)mxCalloc((N+3),sizeof(double));
    char *tmp=(char*)prtmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    prold=pr;
    pr=(double*)tmp;
  }
  if(pi!=NULL && ((size_t)pi&15)!=0) {
    pitmp=(double*)mxCalloc((N+3),sizeof(double));
    char *tmp=(char*)pitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    piold=pi;
    pi=(double*)tmp;
  }
  if(er!=NULL && ((size_t)er&15)!=0) {
    ertmp=(double*)mxMalloc((NE+3)*sizeof(double));
    char *tmp=(char*)ertmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    memcpy(tmp,er,NE*sizeof(double));
    er=(double*)tmp;
  }
  if(ei!=NULL && ((size_t)ei&15)!=0) {
    eitmp=(double*)mxMalloc((NE+3)*sizeof(double));
    char *tmp=(char*)eitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    memcpy(tmp,ei,NE*sizeof(double));
    ei=(double*)tmp;
  }
  if(qr!=NULL && ((size_t)qr&15)!=0) {
    qrtmp=(double*)mxCalloc((NE+3),sizeof(double));
    char *tmp=(char*)qrtmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    qrold=qi;
    qr=(double*)tmp;
  }
  if(qi!=NULL && ((size_t)qi&15)!=0) {
    qitmp=(double*)mxCalloc((NE+3),sizeof(double));
    char *tmp=(char*)qitmp;
    if(((size_t)tmp&15)!=0) {
      tmp+=16-((size_t)tmp&15);
    }
    qiold=qi;
    qi=(double*)tmp;
  }
#endif /*SSE2DIRECT*/
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);
  if (pr != NULL) {
    // mass is real
    if (mi == NULL)
      pot ?
        rinv2(N,zr,zi,mr,pr,pi,smooth,xopt,cutoff,cont) :
        rlog2(N,zr,zi,mr,pr,pi,smooth,xopt,cutoff,cont);
    // mass is complex
    else
      pot ?
        zinv2(N,zr,zi,mr,mi,pr,pi,smooth,xopt,cutoff,cont) :
        zlog2(N,zr,zi,mr,mi,pr,pi,smooth,xopt,cutoff,cont);
    directInteractPanel(panels,pr,pi,zr,zi,Npanel,N);
  }

  // evaluation at (er,ei)
  if (qr != NULL) {
    // mass is real
    if (mi == NULL)
      pot ?
        rinv(N,zr,zi,mr,NE,er,ei,qr,qi,smooth,xopt,cutoff,cont) :
        rlog(N,zr,zi,mr,NE,er,ei,qr,qi,smooth,xopt,cutoff,cont);
    // mass is complex
    else
      pot ?
        zinv(N,zr,zi,mr,mi,NE,er,ei,qr,qi,smooth,xopt,cutoff,cont) :
        zlog(N,zr,zi,mr,mi,NE,er,ei,qr,qi,smooth,xopt,cutoff,cont);
    directInteractPanel(panels,qr,qi,er,ei,Npanel,NE);
  }
  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Direct summation System time: %f\n",timing,0,printtime);
  CudaPrintTime("Direct summation Cuda time",timing,1,printtime);

#ifdef SSE2DIRECT
  if(qiold!=NULL)
    memcpy(qiold,qi,NE*sizeof(double));
  if(qrold!=NULL)
    memcpy(qrold,qr,NE*sizeof(double));
  if(piold!=NULL)
    memcpy(qiold,pi,N*sizeof(double));
  if(prold!=NULL)
    memcpy(prold,pr,N*sizeof(double));
  mxFree(zrtmp);
  mxFree(zitmp);
  mxFree(mrtmp);
  mxFree(mitmp);
  mxFree(ertmp);
  mxFree(eitmp);
  mxFree(prtmp);
  mxFree(pitmp);
  mxFree(qrtmp);
  mxFree(qitmp);
#endif /*SSE2DIRECT*/
}
/*------------------------------------------------------------------------*/
/* Four routines for the logarithmic potential. */
void rlog2(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,
           double *restrict pr,double *restrict pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont)
/* Logarithmic potential for N potentials located at (zr,zi) with real
   mass mr. The contribution is added to (pr,pi). */
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 0.5*log(cutoff)-0.5;
#else
    cutoff = xopt;
    shape = log(cutoff)-0.5;
#endif
    for (int i = 0; i < N; i++) {
      pr[i] -= mr[i]*shape;
      for (int j = i+1; j < N; j++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = log(cutoff)/log(shape+cutoff);
#else
    shape = xopt*xopt;
    if (cont)
      scale = 2.0*log(cutoff)/log(shape+cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++) {
      pr[i] -= mr[i]*0.5*log(shape)*scale;
      for (int j = i+1; j < N; j++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+expint(shape*cutoff));
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+0.5*expint(shape*cutoff*cutoff));
#endif
    for (int i = 0; i < N; i++) {
      pr[i] += mr[i]*0.5*(log(shape)+gamma_)*scale;
      for (int j = i+1; j < N; j++) {
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
/*------------------------------------------------------------------------*/
void zlog2(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,const double *restrict mi,
           double *restrict pr,double *restrict pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont)
/* Same as above but with complex mass (mr,mi). */
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
#ifdef INLINECOMPLEX_DIRECT
        const double re = zr[i]-zr[j],im = zi[i]-zi[j];
        const double P = 0.5*log(re*re+im*im);
#else
        const double P = log(hypot(zr[i]-zr[j],zi[i]-zi[j]));
#endif
        pr[i] -= mr[j]*P;
        pi[i] -= mi[j]*P;
        pr[j] -= mr[i]*P;
        pi[j] -= mi[i]*P;
      }
    break;

  case RANKINE:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 0.5*log(cutoff)-0.5;
#else
    cutoff = xopt;
    shape = log(cutoff)-0.5;
#endif
    for (int i = 0; i < N; i++) {
      pr[i] -= mr[i]*shape;
      pi[i] -= mi[i]*shape;
      for (int j = i+1; j < N; j++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = log(cutoff)/log(shape+cutoff);
#else
    shape = xopt*xopt;
    if (cont)
      scale = 2.0*log(cutoff)/log(shape+cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++) {
      pr[i] -= mr[i]*0.5*log(shape)*scale;
      pi[i] -= mi[i]*0.5*log(shape)*scale;
      for (int j = i+1; j < N; j++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+expint(shape*cutoff));
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+0.5*expint(shape*cutoff*cutoff));
#endif
    for (int i = 0; i < N; i++) {
      pr[i] += mr[i]*0.5*(log(shape)+gamma_)*scale;
      pi[i] += mi[i]*0.5*(log(shape)+gamma_)*scale;
      for (int j = i+1; j < N; j++) {
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
/*------------------------------------------------------------------------*/
void rlog(int N,
          const double *restrict zr,const double *restrict zi,
          const double *restrict mr,
          int NE,
          const double *restrict er,const double *restrict ei,
          double *restrict qr,double *restrict qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont)
/* Logarithmic potential evaluated at NE points (er,ei) with
   potentials at the N points (zr,zi) with real mass mr. The
   contribution is added to (qr,qi). */
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++)
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 0.5*log(cutoff)-0.5;
#else
    cutoff = xopt;
    shape = log(cutoff)-0.5;
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = log(cutoff)/log(shape+cutoff);
#else
    shape = xopt*xopt;
    if (cont)
      scale = 2.0*log(cutoff)/log(shape+cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+expint(shape*cutoff));
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+0.5*expint(shape*cutoff*cutoff));
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
/*------------------------------------------------------------------------*/
void zlog(int N,
          const double *restrict zr,const double *restrict zi,
          const double *restrict mr,const double *restrict mi,
          int NE,
          const double *restrict er,const double *restrict ei,
          double *restrict qr,double *restrict qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont)
/* Same as above but with complex mass (mr,mi). */
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++)
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 0.5*log(cutoff)-0.5;
#else
    cutoff = xopt;
    shape = log(cutoff)-0.5;
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = log(cutoff)/log(shape+cutoff);
#else
    shape = xopt*xopt;
    if (cont)
      scale = 2.0*log(cutoff)/log(shape+cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+expint(shape*cutoff));
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = log(cutoff)/(log(cutoff)+0.5*expint(shape*cutoff*cutoff));
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
/*------------------------------------------------------------------------*/
/* Now follows the 4 corresponding functions for the Harmonic
   potential. */
void rinv2(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,
           double *restrict pr,double *restrict pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont)
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int i = 0; i < N; i++) {
#ifdef SSE2DIRECT
      int j=i+1;
      __m128d mzr=_mm_set1_pd(zr[i]);
      __m128d mzi=_mm_set1_pd(zi[i]);
      __m128d mmi=_mm_set1_pd(mr[i]);
      __m128d pri=_mm_set1_pd(0);
      __m128d pii=_mm_set1_pd(0);
      if (j<N&&j&1) { //SSE need 16 byte alignment, if j is odd, make one normal step
        const double re = zr[i]-zr[j],im = zi[i]-zi[j];
        const double P = 1.0/(re*re+im*im);
        pr[i] -= mr[j]*P*re;
        pi[i] += mr[j]*P*im;
        pr[j] += mr[i]*P*re;
        pi[j] -= mr[i]*P*im;
        j++;
      }
      for (; j < N-1; j+=2) {
        __m128d mre=_mm_load_pd(&zr[j]);
        __m128d mim=_mm_load_pd(&zi[j]);
        __m128d* prj=(__m128d*)&pr[j];
        __m128d* pij=(__m128d*)&pi[j];
        mre=_mm_sub_pd(mre,mzr);
        mim=_mm_sub_pd(mim,mzi);
        __m128d P=_mm_add_pd(_mm_mul_pd(mre,mre),_mm_mul_pd(mim,mim));
        P=_mm_div_pd(_mm_set1_pd(1.0),P);
        __m128d mmr=_mm_load_pd(&mr[j]);
        mmr=_mm_mul_pd(mmr,P);
        P=_mm_mul_pd(P,mmi);
        pri=_mm_add_pd(pri,_mm_mul_pd(mmr,mre));
        pii=_mm_sub_pd(pii,_mm_mul_pd(mmr,mim));
        prj[0]=_mm_sub_pd(prj[0],_mm_mul_pd(P,mre));
        pij[0]=_mm_add_pd(pij[0],_mm_mul_pd(P,mim));
      }
      for (;j < N; j++) {
        const double re = zr[i]-zr[j],im = zi[i]-zi[j];
        const double P = 1.0/(re*re+im*im);
        pr[i] -= mr[j]*P*re;
        pi[i] += mr[j]*P*im;
        pr[j] += mr[i]*P*re;
        pi[j] -= mr[i]*P*im;
      }
      double VC16ALIGN a[2] GCC16ALIGN;
      _mm_store_pd (a, pri);
      pr[i]+=a[0]+a[1];
      _mm_store_pd (a, pii);
      pi[i]+=a[0]+a[1];

#else /*SSE2DIRECT*/
      for (int j = i+1; j < N; j++) {
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
#else
    cutoff = xopt;
    shape = 1.0/(cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
    break;

  case SCULLY:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
#else
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/(cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
    break;

  case OSEEN:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++) {
#ifdef SSE2DIRECT
      int j=i+1;
      __m128d mzr=_mm_set1_pd(zr[i]);
      __m128d mzi=_mm_set1_pd(zi[i]);
      __m128d mmi=_mm_set1_pd(mr[i]);
      __m128d pri=_mm_set1_pd(0);
      __m128d pii=_mm_set1_pd(0);
      if (j<N&&j&1) {
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
      for (; j < N-1; j+=2) {
        __m128d mre=_mm_load_pd(&zr[j]);
        __m128d mim=_mm_load_pd(&zi[j]);
        __m128d* prj=(__m128d*)&pr[j];
        __m128d* pij=(__m128d*)&pi[j];
        mre=_mm_sub_pd(mre,mzr);
        mim=_mm_sub_pd(mim,mzi);
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
        _mm_storel_pd(&ad,mmr);
        if(ad<cutoff) {
          if(ad!=0.0) {
            double Pd;
            _mm_storel_pd(&Pd,P);
            Pd*= -expm1(-shape*ad)*scale;
            P=_mm_loadl_pd(P,&Pd);
          }
          else {
            double Pd=0;
            P=_mm_loadl_pd(P,&Pd);
          }
        }
        _mm_storeh_pd(&ad,mmr);
        if(ad<cutoff) {
          if(ad!=0.0) {
            double Pd;
            _mm_storeh_pd(&Pd,P);
            Pd*= -expm1(-shape*ad)*scale;
            P=_mm_loadh_pd(P,&Pd);
          }
          else {
            double Pd=0;
            P=_mm_loadh_pd(P,&Pd);
          }
        }
#endif
        mmr=_mm_load_pd(&mr[j]);
        mmr=_mm_mul_pd(mmr,P);
        P=_mm_mul_pd(P,mmi);
        pri=_mm_add_pd(pri,_mm_mul_pd(mmr,mre));
        pii=_mm_sub_pd(pii,_mm_mul_pd(mmr,mim));
        prj[0]=_mm_sub_pd(prj[0],_mm_mul_pd(P,mre));
        pij[0]=_mm_add_pd(pij[0],_mm_mul_pd(P,mim));
      }
      for (;j < N; j++) {
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
      double VC16ALIGN a[2] GCC16ALIGN;
      _mm_store_pd (a, pri);
      pr[i]+=a[0]+a[1];
      _mm_store_pd (a, pii);
      pi[i]+=a[0]+a[1];

#else /*SSE2DIRECT*/
      for (int j = i+1; j < N; j++) {
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
/*------------------------------------------------------------------------*/
void zinv2(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,const double *restrict mi,
           double *restrict pr,double *restrict pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont)
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
    break;

  case RANKINE:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
#else
    cutoff = xopt;
    shape = 1.0/(cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
    break;

  case SCULLY:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
#else
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/(cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
    break;

  case OSEEN:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff*cutoff);
#endif
    for (int i = 0; i < N; i++)
      for (int j = i+1; j < N; j++) {
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
    break;
  }
}
/*------------------------------------------------------------------------*/
void rinv(int N,
          const double *restrict zr,const double *restrict zi,
          const double *restrict mr,
          int NE,
          const double *restrict er,const double *restrict ei,
          double *restrict qr,double *restrict qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont)
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int j = 0; j < N; j++) {
#ifdef SSE2DIRECT
      int i=0;
      __m128d mzr=_mm_set1_pd(zr[j]);
      __m128d mzi=_mm_set1_pd(zi[j]);
      __m128d mmi=_mm_set1_pd(mr[j]);
      for (; i < NE-1; i+=2) {
        __m128d mre=_mm_load_pd(&er[i]);
        __m128d mim=_mm_load_pd(&ei[i]);
        __m128d* qrj=(__m128d*)&qr[i];
        __m128d* qij=(__m128d*)&qi[i];
        mre=_mm_sub_pd(mre, mzr);
        mim=_mm_sub_pd(mim, mzi);
        __m128d mad=_mm_add_pd(_mm_mul_pd(mre, mre), _mm_mul_pd(mim, mim));
        __m128d P=_mm_div_pd(mmi, mad);
        double ad;
        _mm_storel_pd(&ad, mad);
        if(ad==0) {
          double Pd=0;
          P=_mm_loadl_pd(P, &Pd);
        }
        _mm_storeh_pd(&ad, mad);
        if(ad==0) {
          double Pd=0;
          P=_mm_loadh_pd(P, &Pd);
        }
        qrj[0]=_mm_sub_pd(qrj[0], _mm_mul_pd(P, mre));
        qij[0]=_mm_add_pd(qij[0], _mm_mul_pd(P, mim));
      }
      if (i < NE) {
        if (er[i] != zr[j] || ei[i] != zi[j]) {
          const double re = er[i]-zr[j],im = ei[i]-zi[j];
          const double P = 1.0/(re*re+im*im);
          qr[i] -= mr[j]*P*re;
          qi[i] += mr[j]*P*im;
        }
      }
#else /*SSE2DIRECT*/
      for (int i = 0; i < NE; i++)
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
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
#else
    cutoff = xopt;
    shape = 1.0/(cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
    break;

  case SCULLY:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
#else
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/(cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
    break;

  case OSEEN:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++) {
#ifdef SSE2DIRECT
      int i=0;
      __m128d mzr=_mm_set1_pd(zr[j]);
      __m128d mzi=_mm_set1_pd(zi[j]);
      __m128d mmi=_mm_set1_pd(mr[j]);
      for (; i < NE-1; i+=2) {
        __m128d mre=_mm_load_pd(&er[i]);
        __m128d mim=_mm_load_pd(&ei[i]);
        __m128d qrj=_mm_load_pd(&qr[i]);
        __m128d qij=_mm_load_pd(&qi[i]);
        mre=_mm_sub_pd(mre, mzr);
        mim=_mm_sub_pd(mim, mzi);
        __m128d mad=_mm_add_pd(_mm_mul_pd(mre, mre), _mm_mul_pd(mim, mim));
        __m128d P=_mm_div_pd(mmi, mad);
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
        qrj=_mm_sub_pd(qrj, _mm_mul_pd(P, mre));
        qij=_mm_add_pd(qij, _mm_mul_pd(P, mim));
        _mm_store_pd(&qr[i],qrj);
        _mm_store_pd(&qi[i],qij);
      }
      if (i < NE) {
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
      for (int i = 0; i < NE; i++) {
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
/*------------------------------------------------------------------------*/
void zinv(int N,
          const double *restrict zr,const double *restrict zi,
          const double *restrict mr,const double *restrict mi,
          int NE,
          const double *restrict er,const double *restrict ei,
          double *restrict qr,double *restrict qi,
          SMOOTHER smooth,double xopt,double cutoff,bool cont)
{
  double shape,scale = 1.0;

  switch (smooth) {
  case DIRAC:
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++)
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
    break;

  case RANKINE:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
#else
    cutoff = xopt;
    shape = 1.0/(cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
    break;

  case SCULLY:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
#else
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/(cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
    break;

  case OSEEN:
#ifdef INLINECOMPLEX_DIRECT
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
#else
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff*cutoff);
#endif
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
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
    break;
  }
}
/*------------------------------------------------------------------------*/
