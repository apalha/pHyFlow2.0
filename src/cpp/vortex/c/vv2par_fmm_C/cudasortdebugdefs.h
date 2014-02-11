#ifdef SORTLIMIT
#define CHECKPARTITIONINGSTRINGSORT1 double hleftlimitvalues[4],hrightlimitvalues[4];\
cudasafe( cudaMemcpy(hleftlimitvalues, leftlimitvalues, imin(splitcount,4)*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy hleftlimitvalues" );\
cudasafe( cudaMemcpy(hrightlimitvalues, rightlimitvalues, imin(splitcount,4)*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy hrightlimitvalues" );\
mexPrintf("limits are [0]: [%f %f], [1]: [%f %f], [2]: [%f %f], [3]: [%f %f]\n",hleftlimitvalues[0],hrightlimitvalues[0],hleftlimitvalues[1],hrightlimitvalues[1],hleftlimitvalues[2],hrightlimitvalues[2],hleftlimitvalues[3],hrightlimitvalues[3]);
#else
#define CHECKPARTITIONINGSTRINGSORT1 
#endif

//Now follows a lot of debugging definitions, that should be removed in the future. They are placed outside to make the code more readable
#define CHECKPARTITIONINGSTRING1 int outllimits[4],outrlimits[4],outlrcount[8];\
cudasafe( cudaMemcpy(outllimits, tmpllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(outrlimits, tmprlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmprlimits" );\
cudasafe( cudaMemcpy(outlrcount, lrcount, imin(splitcount*2,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy lrcount" );\
cudasafe( cudaMemcpy(outsplitpoints, splitpoints, imin(splitcount,4)*sizeof(SORT_REAL), cudaMemcpyDeviceToHost), "cudaMemcpy splitpoints" );\
mexPrintf("left limits are %d %d, right limits are %d %d\n",outllimits[0],outllimits[1],outrlimits[0],outrlimits[1]);\
mexPrintf("lrcount[0]=%d lrcount[1]=%d lrcount[0]=%d lrcount[1]=%d\n",outlrcount[0],outlrcount[1],outlrcount[2],outlrcount[3]);\
mexPrintf("splitpoint=%e outputvector[0]=%d outputvector[1]=%d\n",outsplitpoints[0],outputvectorlocal[0],outputvectorlocal[1]);\
int* dummyindices=(int*)mxMalloc(count*sizeof(int));\
int* buckets=(int*)mxCalloc(count, sizeof(int));\
cudasafe( cudaMemcpy(dummyindices, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
for(int j=0;j<count;j++)\
  buckets[dummyindices[j]]++;\
int failures=0;\
for(int j=0;j<count;j++) {\
  if(buckets[j]!=1)\
    failures++;\
}\
if(failures==0)\
  mexPrintf("All elements accounted for\n");\
else\
  mexPrintf("Invalid split, %d elements do not occur 1 time\n", failures);\
mxFree(dummyindices);\
mxFree(buckets);\
int startllimits[4],startrlimits[4];\
cudasafe( cudaMemcpy(startllimits, llimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(startrlimits, rlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
for(int k=0;k<imin(4, splitcount); k++) {\
  mexPrintf("Limits are %d: [%d,%d)\n",k,startllimits[k],startrlimits[k]);\
}\
mexEvalString("drawnow;");

#define CHECKPARTITIONINGSTRING2 int *indicesbackup=(int*)mxMalloc(count*sizeof(int)); \
int *passivebackup=(int*)mxMalloc(count*sizeof(int));\
cudasafe( cudaMemcpy(indicesbackup, newindices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(passivebackup, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );
        
#define CHECKPARTITIONINGSTRING3 dummyindices=(int*)mxMalloc(count*sizeof(int));\
buckets=(int*)mxCalloc(count,sizeof(int));\
int* dummyindices2=(int*)mxMalloc(count*sizeof(int));\
int* buckets2=(int*)mxCalloc(count,sizeof(int));\
int limitvector[4],limitvector2[4],baselimitvector[4],baselimitvector2[4];\
cudasafe( cudaMemcpy(limitvector, tmpllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(limitvector2, tmprlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(baselimitvector, llimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(baselimitvector2, rlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(dummyindices, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(dummyindices2, newindices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
for(int k=0;k<imin(4, splitcount); k++) {\
  memset(buckets,0,count*sizeof(int));\
  memset(buckets2,0,count*sizeof(int));\
  for(int j=limitvector[k];j<limitvector2[k];j++) {\
    buckets[dummyindices[j]]++;\
    buckets2[dummyindices2[j]]++;\
  }\
  failures=0;\
  for(int j=0;j<count;j++) {\
    if(buckets[j]!=buckets2[j])\
      failures++;\
  }\
  if(failures==0)\
    mexPrintf("All elements accounted for in split between %d and %d\n", limitvector[k], limitvector2[k]);\
  else {\
    mexPrintf("Invalid internal split, %d elements do not occur 1 time\n", failures);\
    failures=0;\
    for(int j=limitvector[k];j<limitvector2[k];j++) {\
      if(dummyindices[j]==dummyindices2[j])\
        failures++;\
    }\
    mexPrintf("%d elements equal before and after split\n", failures);\
    failures=0;\
    for(int j=limitvector[k];j<limitvector2[k];j++) {\
      if(dummyindices[j]==indicesbackup[j])\
        failures++;\
    }\
    mexPrintf("%d elements unchanged before and after split\n", failures);\
  }\
  failures=0;\
  for(int j=baselimitvector[k];j<limitvector[k];j++) {\
    if(dummyindices2[j]!=indicesbackup[j])\
      failures++;\
  }\
  if(failures==0)\
    mexPrintf("Left outside region ok (%d to %d)\n",baselimitvector[k], limitvector[k]);\
  else {\
    mexPrintf("Left outside region invalid, %d elements not equal after split\n", failures);\
  }\
  failures=0;\
  for(int j=limitvector2[k];j<baselimitvector2[k];j++) {\
    if(dummyindices2[j]!=indicesbackup[j])\
      failures++;\
  }\
  if(failures==0)\
    mexPrintf("Right outside region ok (%d to %d)\n", limitvector2[k], baselimitvector2[k]);\
  else {\
    mexPrintf("Rights outside region invalid, %d elements not equal after split\n", failures);\
  }\
}\
failures=0;\
for(int j=0;j<count;j++) {\
  if(dummyindices[j]!=passivebackup[j])\
    failures++;\
}\
if(failures==0)\
  mexPrintf("Passive region ok (%d to %d)\n", 0, count);\
else {\
  mexPrintf("Passive region invalid, %d elements not equal after split\n", failures);\
}\
mxFree(dummyindices);\
mxFree(dummyindices2);\
mxFree(buckets);\
mxFree(buckets2);\
mxFree(indicesbackup);\
mxFree(passivebackup);\
mexEvalString("drawnow;");
        
#define CHECKPARTITIONINGSTRING4 if(i%2==1) {\
  dummyindices=(int*)mxMalloc(count*sizeof(int));\
  buckets=(int*)mxCalloc(count, sizeof(int));\
  cudasafe( cudaMemcpy(dummyindices, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
  for(int j=0;j<count;j++)\
    buckets[dummyindices[j]]++;\
  failures=0;\
  for(int j=0;j<count;j++) {\
    if(buckets[j]!=1)\
      failures++;\
  }\
  if(failures==0)\
    mexPrintf("All elements accounted for\n");\
  else\
    mexPrintf("Invalid split, %d elements do not occur 1 time\n", failures);\
  mxFree(dummyindices);\
  mxFree(buckets);\
}\
mexEvalString("drawnow;");
#define CHECKPARTITIONINGSTRING5 int limitvector[4];\
int *splitbackup=(int*)mxMalloc(count*sizeof(int));\
cudasafe( cudaMemcpy(splitbackup, newindices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(limitvector, tmpllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
mexPrintf("tmpllimits[0]=%d tmpllimits[1]=%d\n",limitvector[0],limitvector[1]);\
cudasafe( cudaMemcpy(limitvector, tmprlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
mexPrintf("tmprlimits[0]=%d tmprlimits[1]=%d\n",limitvector[0],limitvector[1]);\
cudasafe( cudaMemcpy(limitvector, newllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
mexPrintf("newllimits[0]=%d newllimits[1]=%d\n",limitvector[0],limitvector[1]);\
cudasafe( cudaMemcpy(limitvector, newrlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
mexPrintf("newrlimits[0]=%d newrlimits[1]=%d\n",limitvector[0],limitvector[1]);\
mexPrintf("count=%d debugvector=%p\n",count,debugvector);\
cudasafe(cudaMemset(debugvector,0, count*sizeof(double)), "cudaMemset debugvector");\
mexEvalString("drawnow;");
            
#define CHECKPARTITIONINGSTRING6 dummyindices=(int*)mxMalloc(count*sizeof(int));\
int* dummyindices2=(int*)mxMalloc(count*sizeof(int));\
int limitvector2[4], limitvector3[4], limitvector4[4],baselimits[4],baselimits2[4];\
cudasafe( cudaMemcpy(limitvector, newllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy newllimits" );\
cudasafe( cudaMemcpy(limitvector2, newrlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy newrlimits" );\
cudasafe( cudaMemcpy(limitvector3, tmpllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(limitvector4, tmprlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmprlimits" );\
cudasafe( cudaMemcpy(baselimits, llimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy llimits" );\
cudasafe( cudaMemcpy(baselimits2, rlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy rlimits" );\
cudasafe( cudaMemcpy(dummyindices, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy indices" );\
cudasafe( cudaMemcpy(dummyindices2, newindices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy newindices" );\
for(int k=0;k<imin(splitcount, 4);k++) {\
  failures=0;\
  int rights=0;\
  for(int j=limitvector[k];j<limitvector3[k];j++) {\
    if(dummyindices[j]!=dummyindices2[j])\
      failures++;\
    else\
      rights++;\
  }\
  for(int j=limitvector4[k];j<limitvector2[k];j++) {\
    if(dummyindices[j]!=dummyindices2[j])\
      failures++;\
    else\
      rights++;\
  }\
  if(failures==0)\
    mexPrintf("All %d elements correct in copy between %d and %d, and %d and %d\n", rights,limitvector4[k],limitvector2[k],limitvector[k],limitvector3[k]);\
  else {\
    mexPrintf("Invalid copy between %d and %d, and %d and %d, %d elements not copied\n",limitvector4[k],limitvector2[k],limitvector[k],limitvector3[k], failures);\
    if(printcount==0) {\
      double *hdebugvector=(double*)mxMalloc(count*sizeof(double));\
      cudasafe( cudaMemcpy(hdebugvector, debugvector, count*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy debugvector" );\
      for(int m=0;m<count;m++) {\
        if(hdebugvector[m]!=0.0)\
          mexPrintf("dv[%d]=%e\n",m,hdebugvector[m]);\
      }\
      mxFree(hdebugvector);\
      printcount++;\
    }\
  }\
  failures=0;\
  for(int j=baselimits[k];j<limitvector[k];j++) {\
    if(splitbackup[j]!=dummyindices2[j])\
      failures++;\
  }\
  for(int j=limitvector2[k];j<baselimits2[k];j++) {\
    if(splitbackup[j]!=dummyindices2[j])\
      failures++;\
  }\
  if(failures==0)\
    mexPrintf("copy outside check ok\n", rights);\
  else\
    mexPrintf("Invalid copy, %d outside elements changed\n", failures);\
}\
mxFree(dummyindices);\
mxFree(dummyindices2);\
mxFree(splitbackup);
#define CHECKPARTITIONINGSTRING7 cudasafe( cudaMemcpy(outllimits, tmpllimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(outrlimits, tmprlimits, imin(splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmprlimits" );\
cudasafe( cudaMemcpy(outlrcount, lrcount, imin(splitcount*2,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy lrcount" );\
cudasafe( cudaMemcpy(outsplitpoints, splitpoints, imin(splitcount,4)*sizeof(SORT_REAL), cudaMemcpyDeviceToHost), "cudaMemcpy splitpoints" );\
mexPrintf("left limits are %d %d, right limits are %d %d\n", outllimits[0], outllimits[1], outrlimits[0], outrlimits[1]);\
mexPrintf("lrcount[0]=%d lrcount[1]=%d\n", outlrcount[0], outlrcount[1]);\
mexPrintf("splitpoint=%e outputvector[0]=%d outputvector[1]=%d\n", outsplitpoints[0], outputvectorlocal[0], outputvectorlocal[1]);\
if(i==1000000)\
  break;\
mexEvalString("drawnow;");
        
#define CHECKPARTITIONINGSTRING8 dummyindices=(int*)mxMalloc(count*sizeof(int));\
buckets=(int*)mxCalloc(count, sizeof(int));\
cudasafe( cudaMemcpy(dummyindices, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
for(int j=0;j<count;j++)\
  buckets[dummyindices[j]]++;\
failures=0;\
for(int j=0;j<count;j++) {\
  if(buckets[j]!=1)\
    failures++;\
}\
if(failures==0)\
  mexPrintf("All elements accounted for before singlethreadpartition\n");\
else\
  mexPrintf("Invalid split, %d elements do not occur 1 time\n", failures);\
mxFree(dummyindices);\
mxFree(buckets);\
mexEvalString("drawnow;");
    
#define CHECKPARTITIONINGSTRING9 dummyindices=(int*)mxMalloc(count*sizeof(int));\
buckets=(int*)mxCalloc(count, sizeof(int));\
cudasafe( cudaMemcpy(dummyindices, indices, count*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
for(int j=0;j<count;j++)\
  buckets[dummyindices[j]]++;\
failures=0;\
for(int j=0;j<count;j++) {\
  if(buckets[j]!=1)\
    failures++;\
}\
if(failures==0)\
  mexPrintf("All elements accounted for after singlethreadpartition\n");\
else\
  mexPrintf("Invalid split, %d elements do not occur 1 time\n", failures);\
mxFree(dummyindices);\
mxFree(buckets);\
cudasafe( cudaMemcpy(startllimits, newllimits, imin(2*splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
cudasafe( cudaMemcpy(startrlimits, newrlimits, imin(2*splitcount,4)*sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy tmpllimits" );\
for(int k=0;k<imin(4, 2*splitcount); k++) {\
  mexPrintf("Output limits are %d: [%d,%d)\n",k,startllimits[k],startrlimits[k]);\
}\
mexEvalString("drawnow;");
