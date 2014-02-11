#include <stdlib.h>
#include<stdio.h>
#include "cuda.h" //why does the code compile without this?
#ifndef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#include "mex.h"
#endif
#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

//Like cudasafe, but takes two messages
void cudatimingsafe( cudaError_t error,const char* basemessage,const char* message)
{
    char errormessage[1000];
    if(error!=cudaSuccess) {
        snprintf(errormessage,999,"%s %s : %s (%i) \n",basemessage,message,cudaGetErrorString(error),error);
        cudaGetLastError(); //remove the error message
        //mexErrMsgTxt(errormessage);
        printf("Error message goes here!");
    }
}

//create events
void cudaSafeEventCreate(cudaEvent_t *event, const char* str) //only a wrapper to be able to call it from cpp files
{
  cudatimingsafe(cudaEventCreate(event),"cudaEventCreate", str);
}
//start and stop timer
void cudaSafeEventRecord(cudaEvent_t event, const char* str) //only a wrapper to be able to call it from cpp files
{
  cudatimingsafe(cudaEventRecord(event),"cudaEventRecord", str);
}
//synchronizes with GPU
void cudaSafeEventSynchronize(cudaEvent_t event, const char* str) //only a wrapper to be able to call it from cpp files
{
  cudatimingsafe(cudaEventSynchronize(event),"cudaEventSynchronize", str);
}
//get elapsed time
void cudaSafeEventElapsedTime(float* elapsedtime, cudaEvent_t startevent, cudaEvent_t stopevent, const char* str) //only a wrapper to be able to call it from cpp files
{
  cudatimingsafe(cudaEventElapsedTime(elapsedtime, startevent, startevent),"cudaEventElapsedTime", str);
}
//remove event
void cudaSafeEventDestroy(cudaEvent_t event, const char* str) //only a wrapper to be able to call it from cpp files
{
  cudatimingsafe(cudaEventDestroy(event),"cudaEventDestroy", str);
}
//create events and start the first timer
void cudaTimingCreateAndStart(cudaEvent_t *startevent, cudaEvent_t *stopevent, const char* startstr, const char *stopstr)
{
  cudatimingsafe(cudaEventCreate(startevent),"cudaEventCreate", startstr);
  cudatimingsafe(cudaEventCreate(stopevent),"cudaEventCreate", stopstr);
  cudatimingsafe(cudaEventRecord(*startevent),"cudaEventRecord", startstr);
}
//synchronizes with GPU, and prints the time between start and stop event and removes the two events
void cudaTimingSyncPrintAndDestroy(cudaEvent_t startevent, cudaEvent_t stopevent,double* timing,int index,int printtime,const char* outputstr, const char* startstr, const char *stopstr)
{
  float elapsedtime;
  cudatimingsafe(cudaEventSynchronize(stopevent), "cudaEventSynchronize", stopstr);
  cudatimingsafe(cudaEventElapsedTime(&elapsedtime, startevent, stopevent), "cudaEventElapsedTime",startstr);
  cudatimingsafe(cudaEventDestroy(stopevent), "cudaEventDestroy",stopstr);
  cudatimingsafe(cudaEventDestroy(startevent), "cudaEventDestroy",startstr);
  if(printtime)
    mexPrintf("%s: %f\n",outputstr, elapsedtime/1000);
  if(timing!=NULL)
    timing[index]=elapsedtime/1000;
}
