/*_________________________________________________________________________
* vv2parCC_python_066DP.cu
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
* Developer: Giuseppe Tescione
*
* Version: 0.6.6DP (alpha) - 20120814
* simple cut-off constant for desingularization
* double precision (for GPUs of computing capability 2.x)
*________________________________________________________________________*/

#include "Python.h"
#include "arrayobject.h"
#include <math.h>   //C standard basic mathematical operations
#include <omp.h>
#include <stdlib.h>  //C standard general utilities library
#include <stdio.h>
#include <string.h> //C standard constant and function declaration library

/*--------------- Python callable function defintions -------------------*/
static PyObject *vv2parGauss_C_cpu(PyObject *self, PyObject *args);
void init_vv2parGauss_C_cpu();

/*----------------   C utility function definitions ---------------------*/

// define the doctrings for the functions so that help is available in Python
static char vv2parGauss_C_cpu_docs[] =
"vv2parGauss_C_cpu computes induced velocities of particle distribution for a Gaussian kernel.\n"
"\n"
"Usage\n"
"-----\n"
"	vv2parGauss_C_cpu(xBlob,yBlob,wBlob,xTarget,yTarget,ksigmasqr)\n"
"\n"
"Parameters\n"
"----------\n"
"	xBlob :: the x coordinates of the vortex blobs\n"
"        -----    (type: numpy.ndarray (float64); shape: (nBlobs,))\n"
"\n"
"	yBlob :: the y coordinates of the vortex blobs\n"
"        -----    (type: numpy.ndarray (float64); shape: (nBlobs,))\n"
"\n"
"	wBlob :: the circulations associated to each of the vortex blobs\n"
"        -----    (type: numpy.ndarray (float64); shape: (nBlobs,))\n"
"\n"
"	xTarget :: the x coordinates of the target points where to compute\n"
"        -------    velocity\n"
"                  (type: numpy.ndarray (float64); shape: (nTargets,))\n"
"\n"
"	yTarget :: the y coordinates of the target points where to compute\n"
"        -------    velocity\n"
"                  (type: numpy.ndarray (float64); shape: (nTargets,))\n"
"\n"
"       ksigmasqr :: k times the square of the core size of all the vortex blobs\n"
"       ---------    (type: float (64bits); shape: single value)\n"
"\n"
"Returns\n"
"-------\n"
"       velocities :: the induced velocities in each of the (xTarget,yTarget)\n"
"                     points\n"
"                     (type: numpy.ndarray; shape: (nEval,2))\n"
"\n"
"First added:	2013-06-07\n"
"\n"
"Copyright (C) 2013 Artur Palha\n"
"                   pHyFlow";


/*------------------- Set up the methods table for Python ---------------*/
static PyMethodDef _vv2parGauss_C_cpu[] = {
   {"vv2parGauss_C_cpu", vv2parGauss_C_cpu, METH_VARARGS,vv2parGauss_C_cpu_docs},
   {NULL} /*Sentiner - marks the end of this Python structure*/
};

/*--------------Initialize the vv2par functions -------------------------*/
void init_vv2parGauss_C_cpu(){
   (void) Py_InitModule("_vv2parGauss_C_cpu", _vv2parGauss_C_cpu);
   import_array(); /* Must be present for Numpy */
}

/*-------------------- Python wrapping function - HOST -----------------------*/
static PyObject *vv2parGauss_C_cpu(PyObject *self, PyObject *args){
// Python wrapping function, pass data between Python and C and call CUDA Host Code.
    // declare variables
    PyArrayObject *xB;
    PyArrayObject *yB;
    PyArrayObject *wB;
    PyArrayObject *xT;
    PyArrayObject *yT;
    PyArrayObject *vx;
    PyArrayObject *vy;
    double *cxB,*cyB,*cwB, *cxT, *cyT, *cvx,*cvy;
    double EPS = 1e-12;
    double ksigmasqr;
    double inv_2pi = 0.5/3.141592653589793;
    double xk, yk;
    double s, dx, dy, radius_squared;
    double vxk,vyk;
    int nBlobs;
    int nTargets;
    int dims[2];
    int k, m;
    int max_threads;

    /* Parse tuples separately since args will differ between C functions*/
    if(!PyArg_ParseTuple(args, "O!O!O!O!O!d", &PyArray_Type, &xB, &PyArray_Type, &yB, &PyArray_Type, &wB,  &PyArray_Type, &xT, &PyArray_Type, &yT, &ksigmasqr)) return NULL;
    if(NULL == xB) return NULL; // if something went wrong and blobs matrix is null, exit

    // get the number of particles
    nBlobs = xB->dimensions[0];
    nTargets = dims[0] = xT->dimensions[0];
    dims[1] = 1; // the 2 is the number of components of the velocity

    // allocate memory space for the induced velocities matrix
    vx = (PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);
    vy = (PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);

    // pointers to arrays of double precision values of Python matrices
    cxB = (double*)xB->data;
    cyB = (double*)yB->data;
    cwB = (double*)wB->data;
    cxT = (double*)xT->data;
    cyT = (double*)yT->data;
    cvx = (double*)vx->data;
    cvy = (double*)vy->data;

    #define PRIVATE_VARS vxk, vyk, dx, dy, s, radius_squared

    #pragma omp parallel for shared(cvx,cvy,cxB,cyB,cwB,EPS, inv_2pi,ksigmasqr) private(m, PRIVATE_VARS)
    for(k=0; k<nTargets; k++){
        //xk = cxB[k];
        //yk = cyB[k];
        vxk = 0.0;
        vyk = 0.0;
        for(m=0; m<nBlobs; m++){
            dx = cxT[k]-cxB[m];
            dy = cyT[k]-cyB[m];
            radius_squared = dx*dx + dy*dy + EPS;
            s = (cwB[m]/radius_squared)*(1.0 - exp(-radius_squared/ksigmasqr));
            //exp_decay = 1.0 - exp(-radius_squared/core_size);
            vxk -= dy * s;// * exp_decay;
            vyk += dx * s;// * exp_decay;
        }
        cvx[k] = inv_2pi*vxk;
        cvy[k] = inv_2pi*vyk;
    }

    PyObject *tupleresult = PyTuple_New(2);
    PyTuple_SetItem(tupleresult, 0, PyArray_Return(vx));
    PyTuple_SetItem(tupleresult, 1, PyArray_Return(vy));
    return tupleresult;

}













