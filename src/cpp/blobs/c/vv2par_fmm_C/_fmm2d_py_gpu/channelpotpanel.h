#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __channelpotpanel_h
#define __channelpotpanel_h

#include "fmm.h"
#include "panel.h"
void channelsortpanel(MPexp* This,const panel *panels,int* tmpptr1,int* tmpptr2,
                      int** panelptrlist,int endlevel);
void replicatepanelpositions(panel* panels, int Npanel,double channelheight);
void replicatepanelindices(MPexp *This,channelparam* cparams,int Npanel);
void initstreamcoeffpanel(double* ucoeff, double *dcoeff, const panel* panels,
                          mpexp *box, double z0, double d0,
                          int pmaxs,double sigma,double channelheight);
void directInteractChannelpotPanel(const panel *panels, double *qr, double *qi,const double* er,
                         const double* ei, const int Npanel,const int N,const double channelheight);

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif
#ifndef INV2PI
#define INV2PI 0.159154943091895
#endif
#endif
