/* getgpuinfo.cu -- Prints information of the installed CUDA GPU-card(s). */

/* A. Goude 2011-04-12 */

#include "cuda.h"
#include "mex.h"

/*------------------------------------------------------------------------*/
void cudasafe(cudaError_t error,char* message)
/* Function-call wrapper. */
{
    if (error != cudaSuccess) {
      mexPrintf("ERROR: %s : %i\n",message,error);
      exit(-1);
    }
}
/*------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  cudaDeviceProp info;
  int count,i;

  cudasafe(cudaGetDeviceCount(&count),"cudaGetDeviceCount");
  for(i = 0; i < count; i++) {
    cudasafe(cudaGetDeviceProperties(&info,i),"cudaGetDeviceProperties");
    mexPrintf("\nDevice #%d\n",i);
    mexPrintf("---------\n");
    mexPrintf("Device name:\t\t%s\n",info.name);
    mexPrintf("totalGlobalMem:\t\t%d bytes\n",info.totalGlobalMem);
    mexPrintf("sharedMemPerBlock:\t%d bytes\n",info.sharedMemPerBlock);
    mexPrintf("regsPerBlock:\t\t%d\n",info.regsPerBlock);
    mexPrintf("warpSize:\t\t%d threads\n",info.warpSize);

    mexPrintf("memPitch:\t\t%d bytes\n",info.memPitch);
    mexPrintf("maxThreadsPerBlock:\t%d\n",info.maxThreadsPerBlock);
    mexPrintf("maxThreadsDim:\t\t[%d %d %d]\n",
              info.maxThreadsDim[0],
              info.maxThreadsDim[1],
              info.maxThreadsDim[2]);
    mexPrintf("maxGridSize:\t\t[%d %d %d]\n",
              info.maxGridSize[0],
              info.maxGridSize[1],
              info.maxGridSize[2]);
    mexPrintf("totalConstMem:\t\t%d bytes\n\n",info.totalConstMem);

    mexPrintf("Compute Capability:\t%d.%d\n",info.major,info.minor);
    mexPrintf("clockRate:\t\t%d kHz\n",info.clockRate);
    mexPrintf("textureAlignment:\t%d\n",info.textureAlignment);
    mexPrintf("deviceOverlap:\t\t%d\n",info.deviceOverlap);
    mexPrintf("multiProcessorCount:\t%d\n",info.multiProcessorCount);
    if (info.kernelExecTimeoutEnabled)
      mexPrintf("kernelExecTimeout:\tEnabled\n");
    else
      mexPrintf("kernelExecTimeout:\tDisabled\n");
    if (info.integrated)
      mexPrintf("integrated:\t\tmotherboard GPU\n");
    else
      mexPrintf("integrated:\t\tcomponent\n");
    mexPrintf("canMapHostMemory:\t%d\n",info.canMapHostMemory);
    switch (info.computeMode) {
    case cudaComputeModeDefault:
      mexPrintf("computeMode:\t\tcudaComputeModeDefault\n"); break;
    case cudaComputeModeExclusive:
      mexPrintf("computeMode:\t\tcudaComputeModeExclusive\n"); break;
    case cudaComputeModeProhibited:
      mexPrintf("computeMode:\t\tcudaComputeModeProhibited\n"); break;
    default:
      mexPrintf("computeMode:\t\tUNKNOWN\n"); break;
    }
    mexPrintf("maxTexture1D:\t\t%d\n",info.maxTexture1D);
    mexPrintf("maxTexture2D:\t\t[%d %d]\n\n",
              info.maxTexture2D[0],
              info.maxTexture2D[1]);

    mexPrintf("maxTexture3D:\t\t[%d %d %d]\n",
              info.maxTexture3D[0],
              info.maxTexture3D[1],
              info.maxTexture3D[2]);
/*
    mexPrintf("maxTexture2DArray:\t[%d %d %d]\n",
              info.maxTexture2DArray[0],
              info.maxTexture2DArray[1],
              info.maxTexture2DArray[2]);
*/
    mexPrintf("concurrentKernels:\t%d\n",info.concurrentKernels);
    mexPrintf("ECCEnabled:\t\t%d\n",info.ECCEnabled);
    mexPrintf("pciBusID:\t\t%d\n",info.pciBusID);
    mexPrintf("pciDeviceID:\t\t%d\n",info.pciDeviceID);
    mexPrintf("tccDriver:\t\t%d\n",info.tccDriver);
  }
}
/*------------------------------------------------------------------------*/
