/*_________________________________________________________________________
* ww2parCC_python_066DP.cu
*
* Python wrapping to manage data passing and call the CUDA kernel code.
* Pass the input data from python to C environment;
* Copy input data from host to device;
* Set parallel computation;
* Call device code;
* Copy output data back from device to host;
* Pass the output data from C back to Matlab environment.
*
* DUWIND- Delft University Wind Energy Research Institute
* Developer: Artur Palha
*            Giuseppe Tescione
*
* Version: 0.8DP (alpha) - 20130508
* Gaussian kernel of order 1
* double precision (for GPUs of computing capability 2.x)
*________________________________________________________________________*/

#include "Python.h"
#include "arrayobject.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>   //C standard basic mathematical operations
#include <cstdlib>  //C standard general utilities library
#include <stdio.h>
#include <string.h> //C standard constant and function declaration library

/*--------------------- CUDA Kernel Code - Device -----------------------*/
#include "ww2parGauss_device.cu"

/*--------------- Python callable function defintions -------------------*/
extern "C" static PyObject *ww2parGauss_gpu(PyObject *self, PyObject *args);
extern "C" void init_ww2parGauss_gpu();

/*----------------   C utility function definitions ---------------------*/


/*---------------- cuda utility function definitions --------------------*/
__global__ void ww2par_kernel(void *cxBlob_gpu_ondevice, void *cyBlob_gpu_ondevice, void *cwBlob_gpu_ondevice, void *cxTarget_gpu_ondevice, void *cyTarget_gpu_ondevice,void *cvorticity_gpu_ondevice);
__device__ double ww2par_thread(double THR_vorticity, double THR_xTarget, double THR_yTarget, double THR_xBlob, double THR_yBlob, double THR_wBlob);
__device__ double ww2par_block(double BLK_xTarget, double BLK_yTarget, double BLK_vorticity);

// define the doctrings for the functions so that help is available in Python
static char ww2parGauss_gpu_docs[] =
"ww2parGauss_gpu computes induced velocities of particle distribution for a\n"
"Gaussian kernel using the GPU.\n"
"\n"
"Usage\n"
"-----\n"
"	ww2parGauss_gpu(xBlob,yBlob,wBlob,xTarget,yTarget,ksigmasqr,blocksize)\n"
"\n"
"Parameters\n"
"----------\n"
"	xBlob :: the x coordinates of the vortex blobs (nBlobs=N*blocksize)\n"
"        -----    (type: numpy.ndarray (float64); shape: (nBlobs,))\n"
"\n"
"	yBlob :: the y coordinates of the vortex blobs (nBlobs=N*blocksize)\n"
"        -----    (type: numpy.ndarray (float64); shape: (nBlobs,))\n"
"\n"
"	wBlob :: the circulations associated to each of the vortex blobs\n"
"        -----    (nBlobs=N*blocksize)\n" 
"                (type: numpy.ndarray (float64); shape: (nBlobs,))\n"
"\n"
"	xTarget :: the x coordinates of the target points where to compute\n"
"        -------    velocity\n"
"                  (type: numpy.ndarray (float64); shape: (nTargets,))\n"
"\n"
"	yTarget :: the y coordinates of the target points where to compute\n"
"        -------    velocity\n"
"                  (type: numpy.ndarray (float64); shape: (nTargets,))\n"
"\n"
"       ksigmasqr :: the square of the core size of all the vortex blobs\n"
"       ---------    (type: float (64bits); shape: single value)\n"
"\n"
"       blocksize :: the blocksize of the memory used in the gpu\n"
"       ---------    (type: int (any); shape: single value)\n"
"\n"
"Returns\n"
"-------\n"
"       vorticity :: the induced velocities in each of the (xTarget,yTarget)\n"
"                     points\n"
"                     (type: numpy.ndarray; shape: (nEval,2))\n"
"\n"
"First added:	2013-06-24\n"
"\n"
"Copyright (C) 2013 Artur Palha\n"
"                   pHyFlow";

/*------------------- Set up the methods table for Python ---------------*/
static PyMethodDef _ww2parGauss_gpu[] = {
   {"ww2parGauss_gpu", ww2parGauss_gpu, METH_VARARGS,ww2parGauss_gpu_docs},
   {NULL} /*Sentiner - marks the end of this Python structure*/
};

/*--------------Initialize the ww2par functions -------------------------*/
void init_ww2parGauss_gpu(){
   (void) Py_InitModule("_ww2parGauss_gpu", _ww2parGauss_gpu);
   import_array(); /* Must be present for Numpy */
}

/*-------------------- Python wrapping function - HOST -----------------------*/
static PyObject *ww2parGauss_gpu(PyObject *self, PyObject *args){
// Python wrapping function, pass data between Python and C and call CUDA Host Code.
  // declare variables
  PyArrayObject *xBlob;
  PyArrayObject *yBlob;
  PyArrayObject *wBlob;
  PyArrayObject *xTarget;
  PyArrayObject *yTarget;
  PyArrayObject *w;
  double *cxBlob,*cyBlob,*cwBlob,*cxTarget,*cyTarget, *cw;
  void *cxBlob_gpu, *cyBlob_gpu, *cwBlob_gpu, *cxTarget_gpu, *cyTarget_gpu, *cw_gpu;
  double myeps = 1e-12;
  double ksigmasqr;
  double inv_pi = 1.0/3.141592653589793;
  
  int blocksize, nParticles, nTargets, nTargetBlocks, nParticleBlocks;
  int dims[2];

  /* Parse tuples separately since args will differ between C functions*/
  if(!PyArg_ParseTuple(args, "O!O!O!O!O!di", &PyArray_Type, &xBlob, &PyArray_Type, &yBlob, &PyArray_Type, &wBlob, &PyArray_Type, &xTarget, &PyArray_Type, &yTarget, &ksigmasqr, &blocksize)) return NULL;
  if(NULL == xBlob) return NULL; // if something went wrong and xBlob matrix is null, exit
  if(NULL == yBlob) return NULL; // if something went wrong and yBlob matrix is null, exit
  if(NULL == wBlob) return NULL; // if something went wrong and wBlob matrix is null, exit
  if(NULL == xTarget) return NULL; // if something went wrong and xTarget matrix is null, exit
  if(NULL == yTarget) return NULL; // if something went wrong and yTarget matrix is null, exit

  // get the number of particles
  nParticles = xBlob->dimensions[0];
  nTargets = dims[0] = xTarget->dimensions[0];
  nTargetBlocks = nTargets/blocksize;
  nParticleBlocks = nParticles/blocksize;
  dims[1] = 1; 

  // allocate memory space for the induced vorticity matrix
  w = (PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);

  // pointers to arrays of double precision values of Python matrices
  cxBlob = (double*)xBlob->data;
  cyBlob = (double*)yBlob->data;
  cwBlob = (double*)wBlob->data;
  cxTarget = (double*)xTarget->data;
  cyTarget = (double*)yTarget->data;
  cw = (double*)w->data;

    // Copy constants to Constant memory of device. (CUDA syntax)
    cudaMemcpyToSymbol(ksigmasqr_gpu, &ksigmasqr,      sizeof(ksigmasqr),      0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(inv_pi_gpu, &inv_pi,      sizeof(inv_pi),      0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(blocksize_gpu,  &blocksize,  sizeof(blocksize),  0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(nParticles_gpu,  &nParticles, sizeof(nParticles), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(nTargets_gpu,  &nTargets, sizeof(nTargets), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(nParticleBlocks_gpu,  &nParticleBlocks, sizeof(nParticleBlocks), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(nTargetBlocks_gpu,  &nTargetBlocks, sizeof(nTargetBlocks), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(myeps_gpu,  &myeps, sizeof(myeps), 0, cudaMemcpyHostToDevice);


    // Allocate space on device for the input (cxBlob_gpu,cyBlob_gpu,cwBlob_gpu,cxTarget_gpu,cyTarget_gpu,) and output (cvelocities_gpu) arrays.
    size_t SIZE_POS = sizeof(double) * nParticles; //double
    size_t SIZE_IND = sizeof(double) * nTargets; //double

    cudaMalloc(&cxBlob_gpu, SIZE_POS);
    cudaMalloc(&cyBlob_gpu, SIZE_POS);
    cudaMalloc(&cwBlob_gpu, SIZE_POS);
    cudaMalloc(&cxTarget_gpu, SIZE_IND);
    cudaMalloc(&cyTarget_gpu, SIZE_IND);
    cudaMalloc(&cw_gpu, SIZE_IND);

    // Copy input array from host to device memory.
    cudaMemcpy(cxBlob_gpu, cxBlob, SIZE_POS, cudaMemcpyHostToDevice);
    cudaMemcpy(cyBlob_gpu, cyBlob, SIZE_POS, cudaMemcpyHostToDevice);
    cudaMemcpy(cwBlob_gpu, cwBlob, SIZE_POS, cudaMemcpyHostToDevice);
    cudaMemcpy(cxTarget_gpu, cxTarget, SIZE_IND, cudaMemcpyHostToDevice);
    cudaMemcpy(cyTarget_gpu, cyTarget, SIZE_IND, cudaMemcpyHostToDevice);

    // RUN KERNEL
    size_t Sharedmemsize = sizeof(double) * 3 * blocksize;
    dim3 threads(blocksize, 1, 1);
    dim3 grid(nTargets/blocksize, 1, 1);
    ww2par_kernel <<<grid, threads, Sharedmemsize>>> (cxBlob_gpu, cyBlob_gpu, cwBlob_gpu, cxTarget_gpu, cyTarget_gpu, cw_gpu);

    // Copy induction array from device to host memory
    cudaMemcpy(cw, cw_gpu, SIZE_IND, cudaMemcpyDeviceToHost);

    // Free memory space of particle position and induction arrays
    cudaFree(cxBlob_gpu);
    cudaFree(cyBlob_gpu);
    cudaFree(cwBlob_gpu);
    cudaFree(cxTarget_gpu);
    cudaFree(cyTarget_gpu);
    cudaFree(cw_gpu);

    // return the Python array with induced vorticities
    return PyArray_Return(w);

}
