/* mexpanel.h -- Property/value inputs for handling panels. */

/* S. Engblom and A. Goude 2011-06-15 */

#ifndef __mexpanel_h
#define __mexpanel_h

#include "mex.h"

typedef void (*getpanelvalFun)(int,void *,const mxArray *);

void getpanelCUTOFF(int,void *,const mxArray *);
void getpanelDIR(int,void *,const mxArray *);
void getpanelSMOOTH(int,void *,const mxArray *);

int getpanelprop(const char *str);
int getpanelsmooth(const char *str);

#endif /* __mexpanel_h */
