#include "mex.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "hornershift.h"

//reference implementations
void downwardshiftref(double* wksp1,int endk,int pshift)
{
  for (int k=pshift; k > endk; k--) {
    for (int j=k; j <= pshift; j++) {
      wksp1[2*j] += wksp1[2*(j-1)];
      wksp1[2*j+1] += wksp1[2*(j-1)+1];
    }
  }
}
/*------------------------------------------------------------------------*/
void upwardshiftref(double* wksp1, int startk,int pshift)
{
  for (int k=startk; k <= pshift; k++) {
    for (int j = pshift-k; j < pshift; j++) {
      wksp1[2*j] += wksp1[2*(j+1)];
      wksp1[2*j+1] += wksp1[2*(j+1)+1];
    }
  }
}
void downwardshiftref2(double* wksp1,double* wksp2,int pshift)
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
void upwardshiftref2(double* wksp1, double* wksp2,int pshift)
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
/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double teststring1[500],teststring2[500],teststring3[500],teststring4[500],
         teststring5[500],teststring6[500],teststring7[500],teststring8[500];
  int pshift,i,j,k,l;
  srand(time(NULL));
  int stoponfail=0,fail=0;
  if(nrhs>=1) {
    if(mxIsDouble(prhs[0]))
      stoponfail=(int)mxGetScalar(prhs[0]);
    else if(mxIsLogicalScalarTrue(prhs[0]))
      stoponfail=1;
  }
  k=0;
  for(int pshift=1;pshift<60;pshift++) {
    for(l=0;l<400;l++) {
      teststring1[l]=teststring2[l]=1.0;
      teststring3[l]=teststring4[l]=((double)rand()/(double)RAND_MAX);
    }
    upwardshiftref(teststring1,1,pshift);
    upwardshiftref(teststring3,1,pshift);
    upwardshift(teststring2,1,pshift);
    upwardshift(teststring4,1,pshift);
    j=0;
    for(i=0;i<=pshift*2+1;i++) {
      if(fabs(teststring3[i]-teststring4[i])>1e-13||!isfinite(teststring4[i]))
        j++;
    }
    if(j!=0) {
      mexPrintf("Error for pshift=%d, %d elements differ upwards\nReference: ",pshift,j);
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring1[i]);
      mexPrintf("\nResult: ");
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring2[i]);
      mexPrintf("\n");
      k++;
    }
  }
  if(k==0)
    mexPrintf("Upward test successful\n");
  else {
    mexPrintf("Upward test failed\n");
    fail=1;
  }
  
  k=0;
  for(int pshift=1;pshift<60;pshift++) {
    for(l=0;l<400;l++) {
      teststring1[l]=teststring2[l]=teststring3[l]=teststring4[l]=1.0;
      teststring5[l]=teststring6[l]=teststring7[l]=teststring8[l]=((double)rand()/(double)RAND_MAX);
    }
    upwardshiftref2(teststring1,teststring2,pshift);
    upwardshiftref2(teststring5,teststring6,pshift);
    upwardshift2(teststring3,teststring4,pshift);
    upwardshift2(teststring7,teststring8,pshift);
    j=0;
    for(i=0;i<=pshift*2+1;i++) {
      if(fabs(teststring5[i]-teststring7[i])>1e-13||!isfinite(teststring7[i]))
        j++;
    }
    if(j!=0) {
      mexPrintf("Error for pshift=%d string1, %d elements differ upwards2\nReference: ",pshift,j);
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring1[i]);
      mexPrintf("\nResult: ");
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring3[i]);
      mexPrintf("\n");
      k++;
    }
    j=0;
    for(i=0;i<=pshift*2+1;i++) {
      if(fabs(teststring6[i]-teststring8[i])>1e-13||!isfinite(teststring8[i]))
        j++;
    }
    if(j!=0) {
      mexPrintf("Error for pshift=%d string2, %d elements differ upwards2\nReference: ",pshift,j);
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring2[i]);
      mexPrintf("\nResult: ");
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring4[i]);
      mexPrintf("\n");
      k++;
    }
  }
  if(k==0)
    mexPrintf("Upward2 test successful\n");
  else {
    mexPrintf("Upward2 test failed\n");
    fail=1;
  }
  
  k=0;
  for(int pshift=1;pshift<60;pshift++) {
    for(l=0;l<400;l++) {
      teststring1[l]=teststring2[l]=1.0;
      teststring3[l]=teststring4[l]=((double)rand()/(double)RAND_MAX);
    }
    downwardshiftref(teststring1,1,pshift);
    downwardshiftref(teststring3,1,pshift);
    downwardshift(teststring2,1,pshift);
    downwardshift(teststring4,1,pshift);
    j=0;
    for(i=0;i<=pshift*2+1;i++) {
      if(fabs(teststring3[i]-teststring4[i])>1e-13||!isfinite(teststring4[i]))
        j++;
    }
    if(j!=0) {
      mexPrintf("Error for pshift=%d, %d elements differ downwards\nReference: ",pshift,j);
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring1[i]);
      mexPrintf("\nResult: ");
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring2[i]);
      mexPrintf("\n");
      k++;
    }
  }
  if(k==0)
    mexPrintf("Downward test successful\n");
  else {
    mexPrintf("Downward test failed\n");
    fail=1;
  }
  
  k=0;
  for(int pshift=1;pshift<60;pshift++) {
    for(l=0;l<400;l++) {
      teststring1[l]=teststring2[l]=teststring3[l]=teststring4[l]=1.0;
      teststring5[l]=teststring6[l]=teststring7[l]=teststring8[l]=((double)rand()/(double)RAND_MAX);
    }
    downwardshiftref2(teststring1,teststring2,pshift);
    downwardshiftref2(teststring5,teststring6,pshift);
    downwardshift2(teststring3,teststring4,pshift);
    downwardshift2(teststring7,teststring8,pshift);
    j=0;
    for(i=0;i<=pshift*2+1;i++) {
      if(fabs(teststring5[i]-teststring7[i])>1e-13||!isfinite(teststring7[i]))
        j++;
    }
    if(j!=0) {
      mexPrintf("Error for pshift=%d string1, %d elements differ downwards2\nReference: ",pshift,j);
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring1[i]);
      mexPrintf("\nResult: ");
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring3[i]);
      mexPrintf("\n");
      k++;
    }
    j=0;
    for(i=0;i<=pshift*2+1;i++) {
      if(fabs(teststring6[i]-teststring8[i])>1e-13||!isfinite(teststring8[i]))
        j++;
    }
    if(j!=0) {
      mexPrintf("Error for pshift=%d string2, %d elements differ downwards2\nReference: ",pshift,j);
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring2[i]);
      mexPrintf("\nResult: ");
      for(i=0;i<=pshift*2+1;i++)
        mexPrintf("%.1f ",teststring4[i]);
      mexPrintf("\n");
      k++;
    }
  }
  if(k==0)
    mexPrintf("Downward2 test successful\n");
  else {
    mexPrintf("Downward2 test failed\n");
    fail=1;
  }
  if(stoponfail&&fail)
    mexErrMsgTxt("Test failed");
}