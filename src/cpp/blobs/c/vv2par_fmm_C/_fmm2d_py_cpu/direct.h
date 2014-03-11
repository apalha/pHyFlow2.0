/* direct.h */

/* S. Engblom 2007-07-08 (Revision) */
/* S. Engblom 2007-06-26 */

#ifndef __direct_h
#define __direct_h

#include "fmm.h"

/* Direct evaluation of the potential from N pointmasses (mr,mi) at
   the positions (zr,zi). */
void directInteract(int N,
                    const double *zr,const double *zi,
                    const double *mr,const double *mi,
                    int NE,
                    const double *er,const double *ei,
                    double *pr,double *pi,
                    double *qr,double *qi,
                    const panel *panels,int Npanel,int pot,
                    SMOOTHER smooth,double xopt,double cutoff,
                    bool cont,double* timing,int printtime);

#endif /* __direct_h */
