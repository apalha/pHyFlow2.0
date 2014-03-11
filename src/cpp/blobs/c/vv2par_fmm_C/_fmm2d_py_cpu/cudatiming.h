#ifdef C_CODE          /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mxCalloc calloc
#define mexPrintf printf
#define mxLogical int
#endif

#ifndef __cudatiming_h
#define __cudatiming_h

#ifdef CUDATIME
void cudaSafeEventCreate(cudaEvent_t *event,const char* str);//wrapper
void cudaSafeEventRecord(cudaEvent_t event,const char* str);//wrapper
void cudaSafeEventSynchronize(cudaEvent_t event,const char* str); //wrapper
void cudaSafeEventElapsedTime(float* elapsedtime,cudaEvent_t startevent,cudaEvent_t stopevent,const char* str); //wrapper
void cudaSafeEventDestroy(cudaEvent_t event,const char* str); //wrapper
void cudaTimingCreateAndStart(cudaEvent_t *startevent, cudaEvent_t *stopevent, const char* startstr, const char *stopstr);
void cudaTimingSyncPrintAndDestroy(cudaEvent_t startevent, cudaEvent_t stopevent,double* timing,int index,int printtime,const char* outputstr, const char* startstr, const char *stopstr);
void cudatimingsafe( cudaError_t error,const char* basemessage,const char* message);

//useful macros
#define TimeCreateAndStart(startevent,stopevent) if (timing||printtime) cudaTimingCreateAndStart(&startevent,&stopevent,#startevent,#stopevent);
#define TimeSyncPrintAndDestroy(startevent,stopevent,index,outputstr) if (timing||printtime) cudaTimingSyncPrintAndDestroy(startevent,stopevent,timing,index,printtime,outputstr,#startevent,#stopevent);
#define TimeEventRecord(stopevent) if (timing||printtime) cudaSafeEventRecord(stopevent, #stopevent);

#endif /* CUDATIME */
#endif /* __cudatiming_h */
