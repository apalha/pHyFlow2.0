// Copyright Artur Palha 2013 and Lento Manickathan
// numpy_m4.cpp
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
using namespace std;

/* PYTHON */
#include <boost/python.hpp>
#include <boost/foreach.hpp>
#include <boost/range/combine.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <pyublas/numpy.hpp>
namespace python = boost::python;


pyublas::numpy_vector<double> m4kernel_1d(pyublas::numpy_vector<double> x, double c){
    /*
        The one dimensional M4' smooth interpolating kernel given by:

                  /
                 | 0                                        if |x| > 2
                 |
        M4'(x) =<  0.5*((2-|x|)^2)*(1-|x|)     
                 |      - (c^2)*((2-|x|)^2)*(1-2*|x|)       if |x| \in [1,2]
                 |
                 | 1 - 2.5*(x^2) + 1.5*(|x|^3)
                 \      - (c^2)*(2 - 9*(x^2) + 6*(|x|^3))   if |x| \in [0,1]
                  \

        Usage
        -----
            interpolationWeights = _M4prime(x,c)

        Parameters
        ----------
            x :: the points where to compute the interpolation weights
            -    (type: pyublas::numpy_vector<double>; shape: (n,1))

            c :: the coefficient that models diffusion, given by:
            _    c = \frac{\sqrt(\nu * delta)}{\delta x^{2}}
                 (default value 0.0 (no diffusion))
                 (type: double (64bits); shape: single value)


        Returns
        -------
            weights :: the interpolation weight with diffusion
                       (type: pyublas::numpy_vector<double>; shape: x.shape())


        First added:     2013-07-15

        Copyright (C) 2013 Artur Palha and Lento Manickathan
                           pHyFlow
    */

    /*
        Reviews: 
                (1) - [Lento Manickathan] interpolation equation modified to 
                        account diffusion.
    */

    // declare and initialize weights as a pyublas::numpy_vector<double> vector
    // of size x.size()
    pyublas::numpy_vector<double> weights(x.size());

    // loop over all the x coordinates and compute the interpolating kernel
    // weights given by the M4' kernel
    #pragma omp parallel for // for loop in openMP
    for(size_t k=0;k<x.size();k++){
        // |x| \in [1,2]
        if ( (std::abs(x[k]) > 1.0) && (std::abs(x[k]) < 2.0) ){
            weights[k] = 0.5*pow(2.0-std::abs(x[k]),2.0)*(1.0-std::abs(x[k])) - pow(c,2.0)*pow(2.0-std::abs(x[k]),2.0)*(1.0-2.0*std::abs(x[k]));
        }
        // |x| \in [0,1]
        else if ( (std::abs(x[k]) >= 0.0) && (std::abs(x[k]) < 1.0) ){
            weights[k] = 1.0 - 2.5*pow(x[k],2.0) + 1.5*pow(std::abs(x[k]),3.0) - pow(c,2.0)*(2.0 - 9.0*pow(x[k],2) + 6.0*pow(std::abs(x[k]),3.0));
        }
        // |x| > 2
        else {
            weights[k] = 0.0;
        }

    }

    // return the weights
    return weights;
}

// initialize the the interpolating kernel function, such that python can see it
void initm4kernel_1d() {;}

static char m4kernel_1d_docs[] =
"The one dimensional M4' smooth interpolating kernel given by:\n"
"\n"
"              /\n"
"             | 0                                       if |x| > 2\n"
"             |\n"
"   M4'(x) = <  0.5*((2-|x|)^2)*(1-|x|)\n"
"             |      - (c^2)*((2-|x|)^2)*(1-2*|x|)      if |x| \in [1,2]\n"
"             |\n"
"             | 1 - 2.5*(x^2) + 1.5*(|x|^3)\n"
"             \      - (c^2)*(2 - 9*(x^2) + 6*(|x|^3))  if |x| \in [0,1]\n"
"              \\ \n"
"\n"
"   Usage \n"
"   ----- \n"
"       interpolationWeights = _M4prime(x,c)\n"
"\n"
"   Parameters\n"
"\n   ----------\n"
"       x :: the points where to compute the interpolation weights\n"
"       -    (type: pyublas::numpy_vector<double>; shape: (n,1))\n"
"\n"
"       c :: the coefficient that models diffusion, given by:\n"
"       -    c = \frac{\sqrt(\nu * delta)}{\delta x^{2}}\n"
"            (default value 0.0 (no diffusion))\n"
"            (type: double (64bits); shape: single value)\n"
"\n"
"\n"
"   Returns\n"
"   -------\n"
"       weights :: the interpolation weight\n"
"                  (type: pyublas::numpy_vector<double>; shape: x.size())\n"
"\n"
"\n"
"   First added:     2013-07-11\n"
"   Last Modified:   2013-08-27\n"
"\n"
"\n   Copyright (C) 2013 Artur Palha and Lento Manickathan\n"
"\n                      pHyFlow\n";


// generate all python definitions for the kernel function
BOOST_PYTHON_MODULE(_interpolation_kernels_1d)
{
    boost::python::def("m4kernel_1d", m4kernel_1d, boost::python::args("x","c"), m4kernel_1d_docs);
}
