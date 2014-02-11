/* panel.h */
/* A. Goude 2010-01-01 */

#ifndef C_CODE          /*if C_CODE is not defined, we are in MEX file mode*/
#include "mex.h"
#else                   /*if in C file mode, redefine all mex functions to c functions*/
#define mxLogical int
#define mxFree free
#define mxMalloc malloc
#define mexPrintf printf
#endif


#ifndef __PANEL_H
#define __PANEL_H

#include <complex>

typedef struct {
    double x1;        /*first x coordinate of the panel */
    double y1;        /*first y coordinate of the panel */
    double x2;        /*second x coordinate of the panel*/
    double y2;        /*second y coordinate of the panel*/
    double k1;        /*slope of panel (y2-y2)/(x2-x1)  */
    double k2;        /*inverse slope of panel (x2-x2)/(y2-y1)*/
    double theta;     /*Angle of the panel */
    double costheta;  /*cosine of the angle*/
    double sintheta;  /*sinus of the angle */
    double lambda;    /*panel length*/
    double strength1; /*panel circulation at start*/
    double strength2; /*panel circulation at end  */
    double istrength1;/*imaginary part of panel strength*/
    double istrength2;/*imaginary part of panel strength*/
    int horizontal;   /*if k2>k1 (which one to use if one may be close to inf*/
    int smoother;     /*smoother to be used for the panel*/
    double cutoff;    /*cutoff for given smoother*/
    mxLogical side;   /*which side to evaluate on, only applies for points on the panel itself*/
} panel;
#include "fmm.h"
//#define DEBUGpanelinbox
#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef PANELDUMMYFACTOR
#define PANELDUMMYFACTOR 2
#endif

#ifndef DIRACSMOOTHER
#define DIRACSMOOTHER 1
#define RANKINESMOOTHER 2
#define HATSMOOTHER 3
#endif
void printpanels(const panel* panels,int count);
int panelinbox(panel* outpanel,const panel* inpanel,dcmplx z0,dcmplx d0);
void readpanels(panel* panels,const double* realpart,const double* imagpart,const double* strength,const double* istrength,mxLogical* side,int* smoother,double* cutoff,int count);
void writepanels(const panel* panels,double* realpart,double* imagpart,double* strength,double* istrength,mxLogical* side,int count);
void expandpanel(const panel* panels,dcmplx z0,dcmplx *coeff,int p);
void farexpandpanel(const panel* panels,dcmplx z0,dcmplx *coeff,int p);
void MPexp_box_panel(MPexp *This,const panel *panels,int panelcount,int* validinput);
void MPexp_split_panel(dcmplx z0,dcmplx d0,const panel *panels,const int *inptr, int *outptr1,int* outptr2,int incnt,int* outcnt1,int* outcnt2,double splitpoint,int xsplit);
void printindexlist(const double *z,const int *ix,int begin, int end);
void printdummylist(const double *dummy,int begin,int end);
void MPexp_partition_dummy(int begin, int *im, int end, double *z0,
        int *ix, const double *z, int Nmax, double* dummy, int Ndummy);
void mpexp_directInteractPanel(const panel *panels,const mpexp *This,double *qr, double *qi,const double* er,const double* ei,const int *jx,
        int begin, int end);
void directInteractPanel(const panel *panels, double *qr, double *qi,const double* er,const double* ei, const int panelcount,const int N);
void createdummylist(double* dummy,const panel* panels,const int* indices,int N,int xval,dcmplx z,dcmplx d);

void panelshrinkbox(mpexp* This,const double *zr,const double *zi,
        const double *er, const double *ei,
        const int *ix, const int *jx,
        int begin1, int end1,
        int begin2, int end2,
        const panel* panels,panel* smallerpanels,double cutoff);
#endif /*__PANEL_H*/
