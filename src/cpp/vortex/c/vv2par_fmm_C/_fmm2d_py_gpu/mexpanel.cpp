/* mexpanel.cpp */

/* S. Engblom and A. Goude 2011-06-15 */

#include <string.h>
#include "matrix.h"
#include "mex.h"
#include "panel.h"

// macros to be independent of fmm2d.cpp
#define ISDOUBLEMATRIX(A) (mxIsDouble(A) && !mxIsSparse(A) && \
                           mxGetNumberOfDimensions(A) == 2)
#define ISREALMATRIX(A) (ISDOUBLEMATRIX(A) && !mxIsComplex(A))
#define ISLOGICALMATRIX(A) (mxIsLogical(A) && !mxIsSparse(A) && \
                            mxGetNumberOfDimensions(A) == 2)
#define ERR_PROPERTY -1

/*------------------------------------------------------------------------*/
int getpanelprop(const char *str)
/* Small perfect hash for the supported panel properties. */
{
  static const char *tab[] = {"cutoff","dir","smooth"};
  const int hash = ((int)str[0]+(int)str[1])%3;

  return strcmp(str,tab[hash]) == 0 ? hash : ERR_PROPERTY;
}
/*------------------------------------------------------------------------*/
void getpanelCUTOFF(int Npanel,void *cutoff_,const mxArray *rhs)
{
  double *cutoff = (double *)cutoff_;

  // accepts double vector or scalar as input
  if (ISREALMATRIX(rhs)) {
    if (mxGetNumberOfElements(rhs) == Npanel) {
      const double *cutoffin = mxGetPr(rhs);
      for (int k = 0; k < Npanel; k++) {
        cutoff[k] = cutoffin[k];
        if(cutoff[k] <= 0.0)
          mexErrMsgTxt("Panel property 'cutoff' must be positive.");
      }
    }
    else if (mxGetNumberOfElements(rhs) == 1) {
      const double cutoffin = *mxGetPr(rhs);
      if (cutoffin < 0.0)
        mexErrMsgTxt("Panel property 'cutoff' must be positive.");
      for (int k = 0; k < Npanel; k++)
        cutoff[k] = cutoffin;
    }
    else
      mexErrMsgTxt("Expecting a matching vector as panel property "
                   "'cutoff'.");
  }
  else
    mexErrMsgTxt("Expecting a real vector or scalar as panel property "
                 "'cutoff'.");
}
/*------------------------------------------------------------------------*/
void getpanelDIR(int Npanel,void *dir_,const mxArray *rhs)
{
  mxLogical *dir = (mxLogical *)dir_;

  // accepts logical vector or scalar as input
  if (ISLOGICALMATRIX(rhs)) {
    if (mxGetNumberOfElements(rhs) == Npanel) {
      const mxLogical *dirin = mxGetLogicals(rhs);
      for (int k = 0; k < Npanel; k++) {
        dir[k] = dirin[k];
      }
    }
    else if (mxGetNumberOfElements(rhs) == 1) {
      const mxLogical dirin = *mxGetLogicals(rhs);
      for (int k = 0; k < Npanel; k++)
        dir[k] = dirin;
    }
    else
      mexErrMsgTxt("Expecting a matching vector as panel property "
                   "'dir'.");
  }
  else
    mexErrMsgTxt("Expecting a logical vector/scalar as panel property "
                 "'dir'.");
}
/*------------------------------------------------------------------------*/
int getpanelsmooth(const char *str)
/* Small perfect hash for the supported panel smoothers. */
{
  static const char *tab[] = {"rankine","dirac","hat"};
  const int hash = ((int)str[0])%3;

  return strcmp(str,tab[hash]) == 0 ? hash : ERR_SMOOTHER;
}
/*------------------------------------------------------------------------*/
void getpanelSMOOTH(int Npanel,void *smooth_,const mxArray *rhs)
{
  int *smooth = (int *)smooth_;

  // accepts double vector/scalar or a string as input
  if (ISREALMATRIX(rhs)) {
    if (mxGetNumberOfElements(rhs) == Npanel) {
      const double *smoothin = mxGetPr(rhs);
      for (int k = 0; k < Npanel; k++) {
        smooth[k] = (int)smoothin[k];
        if (smooth[k] < 0 || 2 < smooth[k])
          mexErrMsgTxt("Panel property 'smooth' must be an integer 0..2.");
      }
    }
    else if (mxGetNumberOfElements(rhs) == 1) {
      const int smoothin = (int)(*mxGetPr(rhs));
      if (smoothin < 0 || smoothin > 2)
        mexErrMsgTxt("Panel property 'smooth' must be an integer 0..2.");
      for (int k = 0; k < Npanel; k++)
        smooth[k] = smoothin;
    }
    else
      mexErrMsgTxt("Expecting a matching vector as panel property "
                   "'smooth'.");
  }
  else if (mxIsChar(rhs)) {
    char str[10];
    int ismooth;

    if (mxGetString(rhs,str,10) != 0 || 
        (ismooth = getpanelsmooth(str)) == ERR_SMOOTHER)
      mexErrMsgTxt("Unknown panel smoother.");
    /* the perfect hash above yields a different order compared to the
       one in panel.h */
    ismooth = ismooth == 0 ? RANKINESMOOTHER :
      ismooth == 1 ? DIRACSMOOTHER : HATSMOOTHER;
    for (int k = 0; k < Npanel; k++)
      smooth[k] = ismooth;
  }
  else
    mexErrMsgTxt("Expecting a real vector/scalar or string as panel property "
                 "'smooth'.");
}
/*------------------------------------------------------------------------*/
