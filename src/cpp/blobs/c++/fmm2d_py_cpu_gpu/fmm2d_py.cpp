/* fmm2d_py.cpp -- Fast multipole method for 2-D potentials in free space
   (Python/C++/CUDA-version) */

/* A. Palha 2013-10-08 (Python/C++ interface, no use of Matlab)*/
/* S. Engblom 2011-06-28 (Panel+Timing I/O) */
/* S. Engblom and A. Goude 2011-04-12 (Major revision) */
/* A. Goude 2010-01-01 (Panels, port to CUDA/GPU) */
/* S. Engblom 2009-05-08 (Revision) */
/* S. Engblom 2007-07-08 (Major revision) */
/* S. Engblom 2006-10-11 (Port to Mex and a major revision) */
/* S. Engblom 2005-01-05 (original code 'fmm2dlp' in C99) */

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <math.h>   //C standard basic mathematical operations
#include <stdlib.h>  //C standard general utilities library
#include <stdio.h>
#include <string.h> //C standard constant and function declaration library
#include <float.h>
using namespace std;

typedef enum {NDIRECT = 0,PANEL = 1,ETA = 2,SORTLIMITS = 3,CUTOFF = 4,
              PRINTTIME = 6,CONT = 7,TOL = 9,
              OUT = 10,CHANNELHEIGHT = 11,SMOOTH = 12,POT = 13,XOPT = 14,
              ERR_PROPERTY = -1} PROPERTY;

#include "fmm.h"
#include "panel.h"

#ifndef CUDASUPPORT
#include "direct.h"
#else
#include "cudaeval.h"
#endif
#ifdef CHANNELPOT
#include "directchannelpot.h"
#include "channelpotpanel.h"
#endif

/* PYTHON */
#include <boost/python.hpp>
#include <boost/foreach.hpp>
#include <boost/range/combine.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <pyublas/numpy.hpp>
namespace python = boost::python;

#ifdef CUDASUPPORT
static int initialized=0;
#endif

// the gpu fmm
#ifdef CUDASUPPORT
boost::python::tuple fmm2d_py(pyublas::numpy_vector<double> xBlob, pyublas::numpy_vector<double> yBlob, pyublas::numpy_vector<double> wBlob,
                                           pyublas::numpy_vector<double> xTarget, pyublas::numpy_vector<double> yTarget,
                                           int Ndirect, double xopt, double cutoff, double tol){
    /*
        FMM calculator of induced velocities given a set of vortex blobs.

        Usage
        -----
            tuple(vx,vy) = fmm2d_py_cpu(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,tol)

        Parameters
        ----------
            xBlob :: the x coordinates of the vortex blobs.
                     (type: pyublas::numpy_vector<double>; shape: (nBlobs,1))

            yBlob :: the y coordinates of the vortex blobs.
                     (type: pyublas::numpy_vector<double>; shape: (nBlobs,1))

            wBlob :: the circulation associated to the vortex blobs.
                     (type: pyublas::numpy_vector<double>; shape: (nBlobs,1))

            xTarget :: the x coordinates of the target points, where to compute induced velocities.
                       (type: pyublas::numpy_vector<double>; shape: (nTargets,1))

            yTarget :: the y coordinates of the target points, where to compute induced velocities.
                       (type: pyublas::numpy_vector<double>; shape: (nTargets,1))

            Ndirect :: the number of neighbor blobs where to use direct calculation of induced velocities.
                       (type: int; shape: single value)

            xopt :: the core size of the Oseen vortex. It is related to the sigma by the following expression
                        xopt = numpy.sqrt(alpha*ksigmasqr)
                    with:
                        alpha = 1.2564312086261696770
                        ksigmasqr = k*sigma*sigma
                        k = 2.0
                    (type: double; shape: single value)

            cutoff :: the threshold after which to consider 1/r^2 decay of vorticity.
                      typically :: cutoff = 5.0*xopt
                      (type: double; shape: single value)

            tol :: the tolerance used to compute the velocities with FMM, if tol == 0 then direct calculation is performed

        Returns
        -------
            tuple(vx,vy) :: the induced velocities at the points (xTarget,yTarget).
                            (type: boost::python::tuple of two pyublas::numpy_vector<double>; shape: (2,1) and each vx and vy have size nTargets)


        First added:     2013-09-25

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    */

    /*
        Reviews:
    */

    // c arrays to pass input c++ arrays to c functions
    const double *xBlob_c,*yBlob_c,*wBlob_c,*xTarget_c,*yTarget_c;

    // default values of FMM properties
    bool cont = true;
    double channelheight=-1;
    SMOOTHER smooth = OSEEN; // Oseen vortex
    float sortlimits[]={2,2,1,1,0,0};
    double eta[]={2,2};
    int pot = 1;
    double *timing = NULL;
    double *wBLob_c_i = NULL;
    double *vx_blob_c = NULL;
    double *vy_blob_c = NULL;
    bool printtime = false;
    panel *panels = NULL;
    int nPanels = 0;

    // initialize the number of blobs and targets
    int nBlobs, nTargets;

    // compute the number of targets, it is assumed that yTarget and xTarget have the same size
    nTargets = xTarget.size();

    // compute the number of blobs, it is assumed that yBlob and xBlob have the same size
    nBlobs = xBlob.size();

    // declare and initialize x and y components of velocity as a pyublas::numpy_vector<double> vector
    // of size of nTargets
    pyublas::numpy_vector<double> vx(nTargets);
    pyublas::numpy_vector<double> vy(nTargets);

    // c arrays to pass output c++ arrays to c functions
    double *vx_c, *vy_c;

    // extract pointers to the inputs and outputs

    // blobs x and y coordinates and circulation
    xBlob_c = &xBlob[0];
    yBlob_c = &yBlob[0];
    wBlob_c = &wBlob[0];

    // targets x and y coordinates
    xTarget_c = &xTarget[0];
    yTarget_c = &yTarget[0];

    // x and y components of induced velocities
    vx_c = &vx[0];
    vy_c = &vy[0];

    // the velocity vectors have to be initialized to zero
    std::fill(vx.begin(),vx.end(),0.0f);
    std::fill(vy.begin(),vy.end(),0.0f);

    // compute the fmm2d velocity

    if (tol == 0.0){ // if tol == 0.0 then we must use direct method to compute the induced velocities
        direct_eval_cuda(nBlobs,xBlob_c,yBlob_c,wBlob_c,wBLob_c_i,nTargets,xTarget_c,yTarget_c,vx_blob_c,vy_blob_c,vx_c,vy_c,panels,nPanels,
                         pot,smooth,xopt,cutoff,cont,timing,printtime);
    }
    else { // if tol > 0.0 then we use fmm to compute the induced velocities
        fmm2dInteract(nBlobs,xBlob_c,yBlob_c,wBlob_c,wBLob_c_i,nTargets,xTarget_c,yTarget_c,vx_blob_c,vy_blob_c,vx_c,vy_c,panels,nPanels,
                  pot,tol,Ndirect,smooth,xopt,cutoff,cont,timing,
                  printtime,sortlimits,eta,channelheight);
    }
    // return the velocities as a python tuple
    return boost::python::make_tuple(vx, vy);
}
#else
// the cpu fmm version
boost::python::tuple fmm2d_py(pyublas::numpy_vector<double> xBlob, pyublas::numpy_vector<double> yBlob, pyublas::numpy_vector<double> wBlob,
                                           pyublas::numpy_vector<double> xTarget, pyublas::numpy_vector<double> yTarget,
                                           int Ndirect, double xopt, double cutoff, double tol){
    /*
        FMM calculator of induced velocities given a set of vortex blobs.

        Usage
        -----
            tuple(vx,vy) = fmm2d_py_cpu(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,tol)

        Parameters
        ----------
            xBlob :: the x coordinates of the vortex blobs.
                     (type: pyublas::numpy_vector<double>; shape: (nBlobs,1))

            yBlob :: the y coordinates of the vortex blobs.
                     (type: pyublas::numpy_vector<double>; shape: (nBlobs,1))

            wBlob :: the circulation associated to the vortex blobs.
                     (type: pyublas::numpy_vector<double>; shape: (nBlobs,1))

            xTarget :: the x coordinates of the target points, where to compute induced velocities.
                       (type: pyublas::numpy_vector<double>; shape: (nTargets,1))

            yTarget :: the y coordinates of the target points, where to compute induced velocities.
                       (type: pyublas::numpy_vector<double>; shape: (nTargets,1))

            Ndirect :: the number of neighbor blobs where to use direct calculation of induced velocities.
                       (type: int; shape: single value)

            xopt :: the core size of the Oseen vortex. It is related to the sigma by the following expression
                        xopt = numpy.sqrt(alpha*ksigmasqr)
                    with:
                        alpha = 1.2564312086261696770
                        ksigmasqr = k*sigma*sigma
                        k = 2.0
                    (type: double; shape: single value)

            cutoff :: the threshold after which to consider 1/r^2 decay of vorticity.
                      typically :: cutoff = 5.0*xopt
                      (type: double; shape: single value)

            tol :: the tolerance used to compute the velocities with FMM, if tol == 0 then direct calculation is performed

        Returns
        -------
            tuple(vx,vy) :: the induced velocities at the points (xTarget,yTarget).
                            (type: boost::python::tuple of two pyublas::numpy_vector<double>; shape: (2,1) and each vx and vy have size nTargets)


        First added:     2013-09-25

        Copyright (C) 2013 Artur Palha
                           pHyFlow
    */

    /*
        Reviews:
    */

    // c arrays to pass input c++ arrays to c functions
    const double *xBlob_c,*yBlob_c,*wBlob_c,*xTarget_c,*yTarget_c;

    // default values of FMM properties
    bool cont = true;
    double channelheight=-1;
    SMOOTHER smooth = OSEEN; // Oseen vortex
    float sortlimits[]={2,2,1,1,0,0};
    double eta[]={2,2};
    int pot = 1;
    double *timing = NULL;
    double *wBLob_c_i = NULL;
    double *vx_blob_c = NULL;
    double *vy_blob_c = NULL;
    bool printtime = false;
    panel *panels = NULL;
    int nPanels = 0;

    // initialize the number of blobs and targets
    int nBlobs, nTargets;

    // compute the number of targets, it is assumed that yTarget and xTarget have the same size
    nTargets = xTarget.size();

    // compute the number of blobs, it is assumed that yBlob and xBlob have the same size
    nBlobs = xBlob.size();

    // declare and initialize x and y components of velocity as a pyublas::numpy_vector<double> vector
    // of size of nTargets
    pyublas::numpy_vector<double> vx(nTargets);
    pyublas::numpy_vector<double> vy(nTargets);

    // c arrays to pass output c++ arrays to c functions
    double *vx_c, *vy_c;

    // extract pointers to the inputs and outputs

    // blobs x and y coordinates and circulation
    xBlob_c = &xBlob[0];
    yBlob_c = &yBlob[0];
    wBlob_c = &wBlob[0];

    // targets x and y coordinates
    xTarget_c = &xTarget[0];
    yTarget_c = &yTarget[0];

    // x and y components of induced velocities
    vx_c = &vx[0];
    vy_c = &vy[0];

    // the velocity vectors have to be initialized to zero
    std::fill(vx.begin(),vx.end(),0.0f);
    std::fill(vy.begin(),vy.end(),0.0f);

    // compute the fmm2d velocity

    if (tol == 0.0){ // if tol == 0.0 then we must use direct method to compute the induced velocities
        directInteract(nBlobs,xBlob_c,yBlob_c,wBlob_c,wBLob_c_i,nTargets,xTarget_c,yTarget_c,vx_blob_c,vy_blob_c,vx_c,vy_c,panels,nPanels,
                   pot,smooth,xopt,cutoff,cont,timing,printtime);
    }
    else { // if tol > 0.0 then we use fmm to compute the induced velocities
        fmm2dInteract(nBlobs,xBlob_c,yBlob_c,wBlob_c,wBLob_c_i,nTargets,xTarget_c,yTarget_c,vx_blob_c,vy_blob_c,vx_c,vy_c,panels,nPanels,
                  pot,tol,Ndirect,smooth,xopt,cutoff,cont,timing,
                  printtime,sortlimits,eta,channelheight);
    }
    // return the velocities as a python tuple
    return boost::python::make_tuple(vx, vy);
}
#endif

// initialize the fmm_cpu function, such that python can see it
void initfmm2d_py() {;}

#ifdef CUDASUPPORT
static char fmm2d_py_docs[] =
"FMM calculator of induced velocities given a set of vortex blobs (GPU version).\n"
"\n"
"Usage\n"
"-----\n"
"    tuple(vx,vy) = fmm2d_py(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,tol)\n"
"\n"
"Parameters\n"
"----------\n"
"\n    xBlob :: the x coordinates of the vortex blobs.\n"
"             (type: numpy.ndarray<double>; shape: (nBlobs,1))\n"
"\n"
"    yBlob :: the y coordinates of the vortex blobs.\n"
"             (type: numpy.ndarray<double>; shape: (nBlobs,1))\n"
"\n"
"    wBlob :: the circulation associated to the vortex blobs.\n"
"             (type: numpy.ndarray<double>; shape: (nBlobs,1))\n"
"\n"
"    xTarget :: the x coordinates of the target points, where to compute induced velocities.\n"
"               (type: numpy.ndarray<double>; shape: (nTargets,1))\n"
"\n"
"    yTarget :: the y coordinates of the target points, where to compute induced velocities.\n"
"               (type: numpy.ndarray<double>; shape: (nTargets,1))\n"
"\n"
"    Ndirect :: the number of neighbor blobs where to use direct calculation of induced velocities.\n"
"               (type: int; shape: single value)\n"
"\n"
"    xopt :: the core size of the Oseen vortex. It is related to the sigma by the following expression\n"
"                xopt = numpy.sqrt(alpha*ksigmasqr)\n"
"            with:\n"
"                alpha = 1.2564312086261696770\n"
"                ksigmasqr = k*sigma*sigma\n"
"                k = 2.0\n"
"            (type: double; shape: single value)\n"
"\n"
"    cutoff :: the threshold after which to consider 1/r^2 decay of vorticity.\n"
"              typically :: cutoff = 5.0*xopt\n"
"              (type: double; shape: single value)\n"
"\n"
"    tol :: the tolerance used to compute the velocities with FMM, if tol == 0 then direct calculation is performed\n"
"\n"
"Returns\n"
"-------\n"
"    (vx,vy) :: the induced velocities at the points (xTarget,yTarget).\n"
"               (type: (numpy.ndarray,numpy.ndarray); shape: (2,1) and each vx and vy have size nTargets)\n"
"\n"
"\n"
"First added:     2013-09-25\n"
"\n"
"Copyright (C) 2013 Artur Palha\n"
"                   pHyFlow\n";
#else
static char fmm2d_py_docs[] =
"FMM calculator of induced velocities given a set of vortex blobs (CPU version).\n"
"\n"
"Usage\n"
"-----\n"
"    tuple(vx,vy) = fmm2d_py(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,tol)\n"
"\n"
"Parameters\n"
"----------\n"
"\n    xBlob :: the x coordinates of the vortex blobs.\n"
"             (type: numpy.ndarray<double>; shape: (nBlobs,1))\n"
"\n"
"    yBlob :: the y coordinates of the vortex blobs.\n"
"             (type: numpy.ndarray<double>; shape: (nBlobs,1))\n"
"\n"
"    wBlob :: the circulation associated to the vortex blobs.\n"
"             (type: numpy.ndarray<double>; shape: (nBlobs,1))\n"
"\n"
"    xTarget :: the x coordinates of the target points, where to compute induced velocities.\n"
"               (type: numpy.ndarray<double>; shape: (nTargets,1))\n"
"\n"
"    yTarget :: the y coordinates of the target points, where to compute induced velocities.\n"
"               (type: numpy.ndarray<double>; shape: (nTargets,1))\n"
"\n"
"    Ndirect :: the number of neighbor blobs where to use direct calculation of induced velocities.\n"
"               (type: int; shape: single value)\n"
"\n"
"    xopt :: the core size of the Oseen vortex. It is related to the sigma by the following expression\n"
"                xopt = numpy.sqrt(alpha*ksigmasqr)\n"
"            with:\n"
"                alpha = 1.2564312086261696770\n"
"                ksigmasqr = k*sigma*sigma\n"
"                k = 2.0\n"
"            (type: double; shape: single value)\n"
"\n"
"    cutoff :: the threshold after which to consider 1/r^2 decay of vorticity.\n"
"              typically :: cutoff = 5.0*xopt\n"
"              (type: double; shape: single value)\n"
"\n"
"    tol :: the tolerance used to compute the velocities with FMM, if tol == 0 then direct calculation is performed\n"
"\n"
"Returns\n"
"-------\n"
"    (vx,vy) :: the induced velocities at the points (xTarget,yTarget).\n"
"               (type: (numpy.ndarray,numpy.ndarray); shape: (2,1) and each vx and vy have size nTargets)\n"
"\n"
"\n"
"First added:     2013-09-25\n"
"\n"
"Copyright (C) 2013 Artur Palha\n"
"                   pHyFlow\n";
#endif


// generate all python definitions for the kernel function
#ifdef CUDASUPPORT
BOOST_PYTHON_MODULE(_fmm2d_py_gpu)
{
    boost::python::def("fmm2d_py", fmm2d_py, boost::python::args("xBlob","yBlob","wBlob","xTarget","yTarget","Ndirect","xopt","cutoff","tol"), fmm2d_py_docs);
}
#else
BOOST_PYTHON_MODULE(_fmm2d_py_cpu)
{
    boost::python::def("fmm2d_py", fmm2d_py, boost::python::args("xBlob","yBlob","wBlob","xTarget","yTarget","Ndirect","xopt","cutoff","tol"), fmm2d_py_docs);
}
#endif
