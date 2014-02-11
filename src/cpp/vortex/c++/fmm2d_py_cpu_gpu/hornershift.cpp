//#define SSE2INTRINSIC
//#define UNROLLLEVEL 7
//#define CROSSORDER

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#include <stdio.h>
#endif

#ifdef UNROLLHORNER //new code may get problems with gcc 3.4 for 32 bit systems due to non-16 byte aligned stack
#include "hornershiftdefs.h"
#elif defined(UNROLLHORNEROLD) //old code works with gcc 3.4 but does not optimize properly with gcc 4.5
#include "hornershiftdefsold.h"
#endif

#ifndef NOINCLUDE //Dummy to remove headers if only preprocessor is being used (gcc -E -DNOINCLUDE)
#ifdef SSE2INTRINSIC
#include <emmintrin.h>
#endif
#endif

#ifdef __cplusplus
#define restrict __restrict
#endif

#ifdef UNROLLHORNER
void upwardshift(double* wksp1,int startk,int pshift) //note: only works for startk>=1
{
  int k=startk;
  #ifdef SSE2INTRINSIC //SSE2 support, use the 128 bit variable instead
  __m128d a[UNROLLLEVEL+1];
  #else
  double a[UNROLLLEVEL+1];
  #ifdef CROSSORDER
  double a2[UNROLLLEVEL+1];
  #endif
  #endif
  #if UNROLLLEVEL>=7 //loop is performed on different levels depending on UNROLLLEVEL
  for (; k <= pshift-6; k+=7) {
    int j;
    INNERUNROLLBASE8(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k) //here, imaginary part is taken directly after real part for each operation
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC) //If not crossorder, take the imaginary part afterwards
    INNERUNROLLBASE8(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
  }
  if (k <= pshift-5)
  #elif UNROLLLEVEL>=6 //The pattern repeats itself for each unroll level
  for (;k <= pshift-5;)
  #endif
  #if UNROLLLEVEL>=6
  {
    int j;
    INNERUNROLLBASE7(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE7(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    k+=6;
  }
  if (k <= pshift-4)
  #elif UNROLLLEVEL>=5
  for (;k <= pshift-4;)
  #endif
  #if UNROLLLEVEL>=5
  {
    int j;
    INNERUNROLLBASE6(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE6(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    k+=5;
  }
  if (k <= pshift-3)
  #elif UNROLLLEVEL>=4
  for (;k <= pshift-3;)
  #endif
  #if UNROLLLEVEL>=4
  {
    int j;
    INNERUNROLLBASE5(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE5(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    k+=4;
  }
  if (k <= pshift-2)
  #elif UNROLLLEVEL>=3
  for (;k <= pshift-2;)
  #endif
  #if UNROLLLEVEL>=3
  {
    int j;
    INNERUNROLLBASE4(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE4(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    k+=3;
  }
  if (k <= pshift-1)
  #elif UNROLLLEVEL>=2
  for (;k <= pshift-1;)
  #endif
  #if UNROLLLEVEL>=2
  {
    int j;
    INNERUNROLLBASE3(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE3(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    k+=2;
  }
  if (k <= pshift)
  #else
  for (;k <= pshift;k++)
  #endif
  {
    int j;
    INNERUNROLLBASE2(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE2(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
  }
}
/*------------------------------------------------------------------------*/
//The m2ps shift
void upwardshift2(double *restrict wksp1,double *restrict wksp2,int pshift)
{
  int k=2;
  #ifdef SSE2INTRINSIC
  __m128d a[UNROLLLEVEL+1];
  #ifdef FULLCROSSORDER //full crossorder requires twice as many variables. Useful for 64 bit systems where there are twice as many xmm registers
  __m128d a2[UNROLLLEVEL+1];
  #endif
  #else
  double a[UNROLLLEVEL+1];
  #ifdef FULLCROSSORDER
  double a2[UNROLLLEVEL+1];
  double a3[UNROLLLEVEL+1];
  double a4[UNROLLLEVEL+1];
  #elif defined(CROSSORDER)
  double a2[UNROLLLEVEL+1];
  #endif
  #endif
  #if UNROLLLEVEL>=7
  for (; k <= pshift-6; k+=7) {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE8(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k) //all four additions in one call
    #else
    INNERUNROLLBASE8(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k) //shift wksp1 first
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE8(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE8(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)//shift wksp2 after
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE8(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
  }
  if (k <= pshift-5)
  #elif UNROLLLEVEL>=6
  for (;k <= pshift-5;)
  #endif
  #if UNROLLLEVEL>=6
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE7(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k)
    #else
    INNERUNROLLBASE7(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE7(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE7(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE7(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
    k+=6;
  }
  if (k <= pshift-4)
  #elif UNROLLLEVEL>=5
  for (;k <= pshift-4;)
  #endif
  #if UNROLLLEVEL>=5
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE6(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k)
    #else
    INNERUNROLLBASE6(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE6(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE6(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE6(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
    k+=5;
  }
  if (k <= pshift-3)
  #elif UNROLLLEVEL>=4
  for (;k <= pshift-3;)
  #endif
  #if UNROLLLEVEL>=4
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE5(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k)
    #else
    INNERUNROLLBASE5(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE5(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE5(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE5(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
    k+=4;
  }
  if (k <= pshift-2)
  #elif UNROLLLEVEL>=3
  for (;k <= pshift-2;)
  #endif
  #if UNROLLLEVEL>=3
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE4(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k)
    #else
    INNERUNROLLBASE4(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE4(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE4(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE4(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
    k+=3;
  }
  if (k <= pshift-1)
  #elif UNROLLLEVEL>=2
  for (;k <= pshift-1;)
  #endif
  #if UNROLLLEVEL>=2
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE3(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k)
    #else
    INNERUNROLLBASE3(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE3(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE3(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE3(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
    k+=2;
  }
  if (k <= pshift)
  #else
  for (;k <= pshift;k++)
  #endif
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE2(NULL,0,STOREMUP2,LOADM2,ADDEUP2,0,pshift-k)
    #else
    INNERUNROLLBASE2(wksp1,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE2(wksp1,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    INNERUNROLLBASE2(wksp2,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE2(wksp2,1,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    #endif
    #endif
  }
}
/*------------------------------------------------------------------------*/
//same pattern for downward as upward. Only different outer loop and different starting j
void downwardshift(double* wksp1,int endk,int pshift)
{
  int k=pshift;
  #ifdef SSE2INTRINSIC
  __m128d a[UNROLLLEVEL+1];
  #else
  double a[UNROLLLEVEL+1];
  #ifdef CROSSORDER
  double a2[UNROLLLEVEL+1];
  #endif
  #endif
  #if UNROLLLEVEL>=7
  for (; k>endk+6; k-=7) {
    int j;
    INNERUNROLLBASE8(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE8(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
  }
  if (k>endk+5)
  #elif UNROLLLEVEL>=6
  for (;k>endk+5;)
  #endif
  #if UNROLLLEVEL>=6
  {
    int j;
    INNERUNROLLBASE7(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE7(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    k-=6;
  }
  if (k>endk+4)
  #elif UNROLLLEVEL>=5
  for (;k>endk+4;)
  #endif
  #if UNROLLLEVEL>=5
  {
    int j;
    INNERUNROLLBASE6(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE6(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    k-=5;
  }
  if (k>endk+3)
  #elif UNROLLLEVEL>=4
  for (;k>endk+3;)
  #endif
  #if UNROLLLEVEL>=4
  {
    int j;
    INNERUNROLLBASE5(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE5(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    k-=4;
  }
  if (k>endk+2)
  #elif UNROLLLEVEL>=3
  for (;k>endk+2;)
  #endif
  #if UNROLLLEVEL>=3
  {
    int j;
    INNERUNROLLBASE4(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE4(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    k-=3;
  }
  if (k>endk+1)
  #elif UNROLLLEVEL>=2
  for (;k>endk+1;)
  #endif
  #if UNROLLLEVEL>=2
  {
    int j;
    INNERUNROLLBASE3(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE3(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    k-=2;
  }
  if (k > endk)
  #else
  for (;k > endk;k--)
  #endif
  {
    int j;
    INNERUNROLLBASE2(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE2(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
  }
}
/*------------------------------------------------------------------------*/
void downwardshift2(double *restrict wksp1,double *restrict wksp2,int pshift)
{
  int k=pshift;
  #ifdef SSE2INTRINSIC
  __m128d a[UNROLLLEVEL+1];
  #ifdef FULLCROSSORDER
  __m128d a2[UNROLLLEVEL+1];
  #endif
  #else
  double a[UNROLLLEVEL+1];
  #ifdef FULLCROSSORDER
  double a2[UNROLLLEVEL+1];
  double a3[UNROLLLEVEL+1];
  double a4[UNROLLLEVEL+1];
  #elif defined(CROSSORDER)
  double a2[UNROLLLEVEL+1];
  #endif
  #endif
  #if UNROLLLEVEL>=7
  for (; k>6; k-=7) {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE8(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE8(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE8(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE8(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE8(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
  }
  if (k>5)
  #elif UNROLLLEVEL>=6
  for (;k>5;)
  #endif
  #if UNROLLLEVEL>=6
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE7(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE7(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE7(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE7(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE7(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
    k-=6;
  }
  if (k>4)
  #elif UNROLLLEVEL>=5
  for (;k>4;)
  #endif
  #if UNROLLLEVEL>=5
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE6(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE6(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE6(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE6(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE6(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
    k-=5;
  }
  if (k>3)
  #elif UNROLLLEVEL>=4
  for (;k>3;)
  #endif
  #if UNROLLLEVEL>=4
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE5(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE5(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE5(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE5(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE5(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
    k-=4;
  }
  if (k>2)
  #elif UNROLLLEVEL>=3
  for (;k>2;)
  #endif
  #if UNROLLLEVEL>=3
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE4(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE4(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE4(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE4(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE4(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
    k-=3;
  }
  if (k>1)
  #elif UNROLLLEVEL>=2
  for (;k>1;)
  #endif
  #if UNROLLLEVEL>=2
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE3(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE3(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE3(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE3(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE3(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
    k-=2;
  }
  if (k > 0)
  #else
  for (;k > 0;k--)
  #endif
  {
    int j;
    #ifdef FULLCROSSORDER
    INNERUNROLLBASE2(NULL,0,STOREMDN2,LOADM2,ADDEDN2,1,k)
    #else
    INNERUNROLLBASE2(wksp1,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE2(wksp1,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    INNERUNROLLBASE2(wksp2,0,STOREMDN,LOADM1,ADDEDN,1,k)
    #if !defined(CROSSORDER) && !defined(SSE2INTRINSIC)
    INNERUNROLLBASE2(wksp2,1,STOREMDN,LOADM1,ADDEDN,1,k)
    #endif
    #endif
  }
}
/*------------------------------------------------------------------------*/
#elif defined(UNROLLHORNEROLD)
//The older version of the shift. More primitive unroll, otherwise, basically the same as the first one

void upwardshift(double* wksp1, int startk, int pshift)
 {
  int k=startk;
  #ifdef SSE2INTRINSIC
  __m128d *wkspm1=(__m128d*)wksp1;
  #endif
  //note: unrolling one additional step did decrease performance
  #if UNROLLLEVEL>=6
  for (; k <= pshift-5; k+=6) {
    INNERUNROLLUP(UPREALIMAG,UPREALIMAGEND,UPWARDADDU6,UPWARDEND6,wksp1);
  }
  if (k <= pshift-4) {
  #elif UNROLLLEVEL>=5
  for(;k<= pshift-4;) {
  #endif
  #if UNROLLLEVEL>=5
    INNERUNROLLUP(UPREALIMAG,UPREALIMAGEND,UPWARDADDU5,UPWARDEND5,wksp1);
    k+=5;
  }
  if (k <= pshift-3) {
  #elif UNROLLLEVEL>=4
  for(;k<= pshift-3;) {
  #endif
  #if UNROLLLEVEL>=4
    INNERUNROLLUP(UPREALIMAG,UPREALIMAGEND,UPWARDADDU4,UPWARDEND4,wksp1);
    k+=4;
  }
  if (k <= pshift-2) {
  #elif UNROLLLEVEL>=3
  for(;k<= pshift-2;) {
  #endif
  #if UNROLLLEVEL >=3
    INNERUNROLLUP(UPREALIMAG,UPREALIMAGEND,UPWARDADDU3,UPWARDEND3,wksp1);
    k+=3;
  }
  if (k <= pshift-1) {
  #elif UNROLLLEVEL>=2
  for(;k<= pshift-1;) {
  #endif
  #if UNROLLLEVEL >=2
    INNERUNROLLUP(UPREALIMAG,UPREALIMAGEND,UPWARDADDU2,UPWARDEND2,wksp1);
    k+=2;
  }
  if (k <= pshift) {
  #else
  for (;k <= pshift;k++) {
  #endif
    INNERUNROLL1UP(UPREALIMAG,wksp1);
  }
}
/*------------------------------------------------------------------------*/
void upwardshift2(double *restrict wksp1,double *restrict wksp2, int pshift)
 {
  int k=2;
  #ifdef SSE2INTRINSIC
  __m128d *wkspm1=(__m128d*)wksp1;
  __m128d *wkspm2=(__m128d*)wksp2;
  #endif
  //note: unrolling one additional step did decrease performance
  #if UNROLLLEVEL>=6
  for (; k <= pshift-5; k+=6) {
    INNERUNROLLUP(UPREALIMAGD,UPREALIMAGENDD,UPWARDADDU6,UPWARDEND6,wksp1);
  }
  if (k <= pshift-4) {
  #elif UNROLLLEVEL>=5
  for(;k<= pshift-4;) {
  #endif
  #if UNROLLLEVEL>=5
    INNERUNROLLUP(UPREALIMAGD,UPREALIMAGENDD,UPWARDADDU5,UPWARDEND5,wksp1);
    k+=5;
  }
  if (k <= pshift-3) {
  #elif UNROLLLEVEL>=4
  for(;k<= pshift-3;) {
  #endif
  #if UNROLLLEVEL>=4
    INNERUNROLLUP(UPREALIMAGD,UPREALIMAGENDD,UPWARDADDU4,UPWARDEND4,wksp1);
    k+=4;
  }
  if (k <= pshift-2) {
  #elif UNROLLLEVEL>=3
  for(;k<= pshift-2;) {
  #endif
  #if UNROLLLEVEL >=3
    INNERUNROLLUP(UPREALIMAGD,UPREALIMAGENDD,UPWARDADDU3,UPWARDEND3,wksp1);
    k+=3;
  }
  if (k <= pshift-1) {
  #elif UNROLLLEVEL>=2
  for(;k<= pshift-1;) {
  #endif
  #if UNROLLLEVEL >=2
    INNERUNROLLUP(UPREALIMAGD,UPREALIMAGENDD,UPWARDADDU2,UPWARDEND2,wksp1);
    k+=2;
  }
  if (k <= pshift) {
  #else
  for (;k <= pshift;k++) {
  #endif
    INNERUNROLL1UP(UPREALIMAGD,wksp1);
  }
}
/*------------------------------------------------------------------------*/
void downwardshift(double* wksp1,int endk, int pshift)
{
  #ifdef SSE2INTRINSIC
  __m128d *wkspm1=(__m128d*)wksp1;
  #endif
  int k=pshift;
  //note: unrolling one additional step did decrease performance
  #if UNROLLLEVEL>=6
  for (; k > 5+endk; k-=6) {
    INNERUNROLLDN(UPREALIMAG,UPREALIMAGEND,DNWARDADDU6,DNWARDEND6,wksp1);
  }
  if (k > 4+endk) {
  #elif UNROLLLEVEL>=5
  for(;k > 4+endk;) {
  #endif
  #if UNROLLLEVEL>=5
    INNERUNROLLDN(UPREALIMAG,UPREALIMAGEND,DNWARDADDU5,DNWARDEND5,wksp1);
    k-=5;
  }
  if (k > 3+endk) {
  #elif UNROLLLEVEL>=4
  for(;k > 3+endk;) {
  #endif
  #if UNROLLLEVEL>=4
    INNERUNROLLDN(UPREALIMAG,UPREALIMAGEND,DNWARDADDU4,DNWARDEND4,wksp1);
    k-=4;
  }
  if (k > 2+endk) {
  #elif UNROLLLEVEL>=3
  for(;k > 2+endk;) {
  #endif
  #if UNROLLLEVEL >=3
    INNERUNROLLDN(UPREALIMAG,UPREALIMAGEND,DNWARDADDU3,DNWARDEND3,wksp1);
    k-=3;
  }
  if (k > 1+endk) {
  #elif UNROLLLEVEL>=2
  for(;k > 1+endk;) {
  #endif
  #if UNROLLLEVEL >=2
    INNERUNROLLDN(UPREALIMAG,UPREALIMAGEND,DNWARDADDU2,DNWARDEND2,wksp1);
    k-=2;
  }
  if (k > endk) {
  #else
  for (;k > endk;k--) {
  #endif
    INNERUNROLL1DN(UPREALIMAG,wksp1);
  }
}
/*------------------------------------------------------------------------*/
void downwardshift2(double *restrict wksp1,double *restrict wksp2, int pshift)
{
  int k=pshift;
  #ifdef SSE2INTRINSIC
  __m128d *wkspm1=(__m128d*)wksp1;
  __m128d *wkspm2=(__m128d*)wksp2;
  #endif
  //note: unrolling one additional step did decrease performance
  #if UNROLLLEVEL>=6
  for (; k > 5; k-=6) {
    INNERUNROLLDN(UPREALIMAGD,UPREALIMAGENDD,DNWARDADDU6,DNWARDEND6,wksp1);
  }
  if (k > 4) {
  #elif UNROLLLEVEL>=5
  for(;k > 4;) {
  #endif
  #if UNROLLLEVEL>=5
    INNERUNROLLDN(UPREALIMAGD,UPREALIMAGENDD,DNWARDADDU5,DNWARDEND5,wksp1);
    k-=5;
  }
  if (k > 3) {
  #elif UNROLLLEVEL>=4
  for(;k > 3;) {
  #endif
  #if UNROLLLEVEL>=4
    INNERUNROLLDN(UPREALIMAGD,UPREALIMAGENDD,DNWARDADDU4,DNWARDEND4,wksp1);
    k-=4;
  }
  if (k > 2) {
  #elif UNROLLLEVEL>=3
  for(;k > 2;) {
  #endif
  #if UNROLLLEVEL >=3
    INNERUNROLLDN(UPREALIMAGD,UPREALIMAGENDD,DNWARDADDU3,DNWARDEND3,wksp1);
    k-=3;
  }
  if (k > 1) {
  #elif UNROLLLEVEL>=2
  for(;k > 1;) {
  #endif
  #if UNROLLLEVEL >=2
    INNERUNROLLDN(UPREALIMAGD,UPREALIMAGENDD,DNWARDADDU2,DNWARDEND2,wksp1);
    k-=2;
  }
  if (k > 0) {
  #else
  for (;k > 0;k--) {
  #endif
    INNERUNROLL1DN(UPREALIMAGD,wksp1);
  }
}
#else
/*------------------------------------------------------------------------*/
//No unroll at all, but support for SSE2
#ifdef SSE2INTRINSIC
void downwardshift(double* wksp1,int endk,int pshift)
{
  __m128d *wksp=(__m128d*)wksp1;
  for (int k=pshift; k > endk; k--) {
    for (int j=k; j <= pshift; j++) {
      wksp[j]=_mm_add_pd(wksp[j],wksp[j-1]);
    }
  }
}
/*------------------------------------------------------------------------*/
void upwardshift(double* wksp1, int startk,int pshift)
{
  __m128d *wksp=(__m128d*)wksp1;
  for (int k=startk; k <= pshift; k++) {
    for (int j = pshift-k; j < pshift; j++) {
      wksp[j]=_mm_add_pd(wksp[j],wksp[j+1]);
    }
  }
}
/*------------------------------------------------------------------------*/
void downwardshift2(double *restrict wksp1,double *restrict wksp2,int pshift)
{
  __m128d *wkspm1=(__m128d*)wksp1;
  __m128d *wkspm2=(__m128d*)wksp2;
  for (int k=pshift; k > 0; k--) {
    for (int j=k; j <= pshift; j++) {
      wkspm1[j]=_mm_add_pd(wkspm1[j],wkspm1[j-1]);
      wkspm2[j]=_mm_add_pd(wkspm2[j],wkspm2[j-1]);
    }
  }
}
/*------------------------------------------------------------------------*/
void upwardshift2(double *restrict wksp1,double *restrict wksp2,int pshift)
{
  __m128d *wkspm1=(__m128d*)wksp1;
  __m128d *wkspm2=(__m128d*)wksp2;
  for (int k=2; k <= pshift; k++) {
    for (int j = pshift-k; j < pshift; j++) {
      wkspm1[j]=_mm_add_pd(wkspm1[j],wkspm1[j+1]);
      wkspm2[j]=_mm_add_pd(wkspm2[j],wkspm2[j+1]);
    }
  }
}
#else //NO SSE2
void downwardshift(double* wksp1,int endk,int pshift)
{
  for (int k=pshift; k > endk; k--) {
    for (int j=k; j <= pshift; j++) {
      wksp1[2*j] += wksp1[2*(j-1)];
      wksp1[2*j+1] += wksp1[2*(j-1)+1];
    }
  }
}
/*------------------------------------------------------------------------*/
void upwardshift(double* wksp1, int startk,int pshift)
{
  for (int k=startk; k <= pshift; k++) {
    for (int j = pshift-k; j < pshift; j++) {
      wksp1[2*j] += wksp1[2*(j+1)];
      wksp1[2*j+1] += wksp1[2*(j+1)+1];
    }
  }
}
/*------------------------------------------------------------------------*/
void downwardshift2(double *restrict wksp1,double *restrict wksp2,int pshift)
{
  for (int k=pshift; k > 0; k--) {
    for (int j=k; j <= pshift; j++) {
      wksp1[2*j] += wksp1[2*(j-1)];
      wksp1[2*j+1] += wksp1[2*(j-1)+1];
      wksp2[2*j] += wksp2[2*(j-1)];
      wksp2[2*j+1] += wksp2[2*(j-1)+1];
    }
  }
}
/*------------------------------------------------------------------------*/
void upwardshift2(double *restrict wksp1,double *restrict wksp2,int pshift)
{
  for (int k=2; k <= pshift; k++) {
    for (int j = pshift-k; j < pshift; j++) {
      wksp1[2*j] += wksp1[2*(j+1)];
      wksp1[2*j+1] += wksp1[2*(j+1)+1];
      wksp2[2*j] += wksp2[2*(j+1)];
      wksp2[2*j+1] += wksp2[2*(j+1)+1];
    }
  }
}
#endif /*SSE2INTRINSIC*/
#endif
