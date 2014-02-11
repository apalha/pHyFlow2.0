#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __directchannelpot_h
#define __directchannelpot_h

void directInteractchannelpot(int N,
                    const double *zr,const double *zi,
                    const double *mr,const double *mi,
                    int NE,
                    const double *er,const double *ei,
                    double *pr,double *pi,
                    double *qr,double *qi,
                    const panel *panels,int Npanel,double channelheight,
                    SMOOTHER smooth,double xopt,double cutoff,
                    bool cont,double* timing,int printtime);
#endif
