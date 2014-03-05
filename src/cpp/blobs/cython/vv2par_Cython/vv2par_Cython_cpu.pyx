# -*- coding: utf-8 -*-
"""
Created on Thu May  2 15:27:01 2013

@author: gorkiana
"""

from __future__ import division
from cython.parallel import prange
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.float64
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.float64_t DTYPE_t
# The builtin min and max functions works with Python objects, and are
# so very slow. So we create our own.
#  - "cdef" declares a function which has much less overhead than a normal
#    def function (but it is not Python-callable)
#  - "inline" is passed on to the C compiler which may inline the functions
#  - The C type "int" is chosen as return type and argument types
#  - Cython allows some newer Python constructs like "a if x else b", but
#    the resulting C file compiles with Python 2.3 through to Python 3.0 beta.
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
# "def" can type its arguments but not have a return type. The type of the
# arguments for a "def" function is checked at run-time when entering the
# function.
#
# The arrays blobs and velocities are typed as "np.ndarray" instances. The only effect
# this has is to a) insert checks that the function arguments really are
# NumPy arrays, and b) make some attribute access like particles.shape[0] much
# more efficient. (In this example this doesn't matter though.)
cimport cython
cimport openmp
@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.cdivision(True)

def vv2par_Cython_cpu(np.ndarray[DTYPE_t, ndim=2] blobs, DTYPE_t sigmasqr):
    assert blobs.dtype == DTYPE
    # The "cdef" keyword is also used within functions to type variables. It
    # can only be used at the top indendation level (there are non-trivial
    # problems with allowing them in other places, though we'd love to see
    # good and thought out proposals for it).
    #
    # For the indices, the "int" type is used. This corresponds to a C int,
    # other C types (like "unsigned int") could have been used instead.
    # Purists could use "Py_ssize_t" which is the proper Python type for
    # array indices.
    cdef int nBlobs = blobs.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=2] velocities = np.zeros([nBlobs, 2], dtype=DTYPE)
    cdef int k, m
    # It is very important to type ALL your variables. You do not get any
    # warnings if not, only much slower code (they are implicitly typed as
    # Python objects).
    cdef int s_from, s_to, t_from, t_to
    # For the value variable, we want to use the same data type as is
    # stored in the array, so we use "DTYPE_t" as defined above.
    # NB! An important side-effect of this is that if "value" overflows its
    # datatype size, it will simply wrap around like in C, rather than raise
    # an error like in Python.
    cdef DTYPE_t EPS    
    cdef DTYPE_t value
    cdef DTYPE_t inv_2pi
    cdef DTYPE_t radiusx
    cdef DTYPE_t radiusy
    cdef DTYPE_t radius_squared
    cdef DTYPE_t s

    cdef DTYPE_t xk    
    cdef DTYPE_t yk
    
    # define floating point precision
    EPS = 1e-12
    # pre-compute 1/2*pi
    inv_2pi = 0.5/np.pi    
    
    with nogil:
        for k in prange(nBlobs):
            xk = blobs[k,0]
            yk = blobs[k,1]
            for m in range(nBlobs):
                #radiusx = blobs[k,0]-blobs[m,0]
                #radiusy = blobs[k,1]-blobs[m,1]
                radiusx = xk-blobs[m,0]
                radiusy = yk-blobs[m,1]
                radius_squared = radiusx*radiusx + radiusy*radiusy + sigmasqr
                s = blobs[m,2]/radius_squared
                velocities[k,0] -= radiusy * s
                velocities[k,1] += radiusx * s
            
            velocities[k,0] *= inv_2pi
            velocities[k,1] *= inv_2pi
         
    return velocities