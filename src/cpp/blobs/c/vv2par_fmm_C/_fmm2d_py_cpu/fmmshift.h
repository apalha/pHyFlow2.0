/* fmmshift.h */

/* S. Engblom 2009-10-08 (Revision) */
/* S. Engblom 2009-04-15 (Major revision) */
/* S. Engblom 2007-05-29 */

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __fmmshift_h
#define __fmmshift_h

#include "fmm.h"

/* Definition of the three shift-operators used by the
   FMM. Performance is critical in general and in particular for the
   m2p-shift. */

void shift_alloc(int p,int pot,int maxm2p);

void shift_m2m(dcmplx r0,dcmplx r1,dcmplx r2,dcmplx r3,
               dcmplx *This,const dcmplx *that);
void shift_m2ps(dcmplx z0,const void *zi,size_t siz,
                const int *ix,int begin,int end,
                const dcmplx *This1,dcmplx *This2,
                const dcmplx *that1,dcmplx *that2);
void shift_p2p(dcmplx r0,dcmplx r1,dcmplx r2,dcmplx r3,
               dcmplx *This,const dcmplx *that);

#endif /* __fmmshift_h */
