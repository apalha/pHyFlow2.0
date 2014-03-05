""" _sourcePanel_Cython_cpu

Source panel cython cpu implementation. 

Description
-----------
Source panel cython cpu (multi-threaded) solver for solving the problem of
*source panels*. Contains two functions: one to assemble the self-induction
matrix for the source panels. The functions is constructed such that the data
of multiple geometry can simply be concatenated into a 1-D array. The second
function evaluates the induced velocity of all the panels on a given list of
x,y evaluation points.

The outer-looping is parallelized using prange.

Reference
---------
Katz, J., & Plotkin, A. (2001). Low-speed aerodynamics. Cambridge 
University Press, [Chapter 11, 11.2.1]

Info
----
By      : Lento Manickathan
Template: Artur Palha

:First Added:   2013-09-31
:Last Modified: 2013-10-08
:Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
:Licence:       GNU GPL version 3 or any later version


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
#cdef inline int int_max(int a, int b): return a if a >= b else b
#cdef inline int int_min(int a, int b): return a if a <= b else b

# "def" can type its arguments but not have a return type. The type of the
# arguments for a "def" function is checked at run-time when entering the
# function.


# The arrays blobs and velocities are typed as "np.ndarray" instances. The only effect
# this has is to a) insert checks that the function arguments really are
# NumPy arrays, and b) make some attribute access like particles.shape[0] much
# more efficient. (In this example this doesn't matter though.)


cimport cython
cimport openmp
from libc.math cimport log, atan2, M_PI, sqrt

@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
@cython.cdivision(True)
def assemble_influenceMatrix(np.ndarray[DTYPE_t, ndim=1] xCP, 
                             np.ndarray[DTYPE_t, ndim=1] yCP,
                             np.ndarray[DTYPE_t, ndim=1] xPanelStart,
                             np.ndarray[DTYPE_t, ndim=1] yPanelStart,
                             np.ndarray[DTYPE_t, ndim=1] xPanelEnd,
                             np.ndarray[DTYPE_t, ndim=1] yPanelEnd,
                             np.ndarray[DTYPE_t, ndim=1] cosAlpha,
                             np.ndarray[DTYPE_t, ndim=1] sinAlpha):
    r"""
    Assemble the self-induction influence matrix :math:`\mathbf{A}` for the 
    source panels. The operation is done using parallel CPU operations.
  
    The source panels are represented by a *constant-strength source element*
    singularity element as given by:
    
        Katz, J., & Plotkin, A. (2001). Low-speed aerodynamics. Cambridge 
        University Press, [Chapter 11, 11.2.1]
  
    The boundary intergral equations are discretized into :math:`\mathbf{M}`
    source panels, where the strength of those panels are the solution to the
    system of :math:`\mathbf{M}` linear equations of form:
    
            ---            ---       ---  ---       --  --
           | a11 a12 ...  a1M |     | sigma1 |     | RHS1 |
           | a21  .           |     |   .    |     |  .   |
           |  .      .        |  o  |   .    |  =  |  .   |
           |  .         .     |     |   .    |     |  .   |
           | aM1    ...   aMM |     | sigmaM |     | RHSM |
            ---             ---      ---  ---       --  --
    
                :math:`\mathbf{A} \cdot \sigma=\mathbf{RHS}`
        
    This function will return the self-induction influence matrix 
    :math:`\mathbf{A}` which can be used to solve a given :math:`\mathbf{RHS}`
    problem.
    

    Usage
    -----
    .. code-block:: python
    
        A = assemble_influenceMatrix(xCP,yCP,xPanelStart,yPanelStart,
                                     xPanelEnd,yPanelEnd,cosAlpha,sinAlpha)


    Parameters 
    ---------- 
    xCP : ndarray (float64), shape (M,)
          the :math:`x`-coordinates of :math:`\mathbf{M}`panel collocation 
          points in `global coordinate system`. For *source panels*, the
          collocation points should be at the center of the panel and slightly
          OUTSIDE.

    yCP : ndarray (float64), shape (M,)
          the :math:`y`-coordinates of :math:`\mathbf{M}`panel collocation 
          points in `global coordinate system`. For *source panels*, the
          collocation points should be at the center of the panel and slightly
          OUTSIDE.        
        
    xPanelStart : ndarray (float64), shape (M,)
                  the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
                  *starting point* in `global coordinate system`. This
                  variable is split so that it can be used for multiple 
                  geometries.

    yPanelStart : ndarray (float64), shape (M,)
                  the :math:`y`-coordinates of :math:`\mathbf{M}`panel corner
                  *starting point* in `global coordinate system`. This
                  variable is split so that it can be used for multiple 
                  geometries.

    xPanelEnd : ndarray (float64), shape (M,)
                the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
                *starting point* in `global coordinate system`. This
                variable is split so that it can be used for multiple 
                geometries.

    yPanelEnd : ndarray (float64), shape (M,)
                the :math:`y`-coordinates of :math:`\mathbf{M}`panel corner
                *starting point* in `global coordinate system`. This
                variable is split so that it can be used for multiple 
                geometries.
                
    cosAlpha : ndarray (float64), shape (M,)
               the :math:`\cos\alpha` of the panel :math:`\mathbf{M}`,
               where the angle :math:`\alpha` is w.r.t to the panel 
               :math:`\mathbf{x}^\prime`-axis and the global
               :math:`\mathbf{x}`-axis.

    sinAlpha : ndarray (float64), shape (M,)
               the :math:`\sin\alpha` of the panel :math:`\mathbf{M}`,
               where the angle :math:`\alpha` is w.r.t to the panel 
               :math:`\mathbf{x}^\prime`-axis and the global
               :math:`\mathbf{x}`-axis.

                    
    Returns
    -------
    A : ndarray (float64), shape (M,M)
        the assembled self-induction influence matrix :math:`\mathbf{A}` for 
        the source panels. *Note:* when :math:`i == j, a_{ii} = 1/2`.


    :First Added:   2013-09-31
    :Last Modified: 2013-10-08
    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version
    
    """                       

    # Panel Parameters
    cdef int M = xCP.shape[0] # Number of panels
    
    # Constants
    cdef DTYPE_t inv_2pi = 1.0/(2.0*M_PI)
    cdef DTYPE_t inv_4pi = inv_2pi/2.0
        
    # Variables of panels
    cdef DTYPE_t Lx,Ly,x2,xP,yP,vxP,vyP,vx,vy

    # Allocate memory for self-induction matrix
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.empty((M,M), dtype=DTYPE) # Influence matrix

    # Iteration parameters
    cdef int i,j
    
    # Iterating through all the collocation points
    for i in prange(M, nogil=True, schedule='static'):
        
        # Iterating through all the panels 
        for j in range(M):
            
            # Calculate the length of panel
            Lx = xPanelEnd[j] - xPanelStart[j]
            Ly = yPanelEnd[j] - yPanelStart[j]
            
            # Calculate the panel end coordinate in the panel coordinate system
            x2 = cosAlpha[j]*Lx + sinAlpha[j]*Ly
    
            # Locating the control point in panel coordinate system
            xP =  cosAlpha[j] * (xCP[i]-xPanelStart[j])   +   sinAlpha[j] * (yCP[i]-yPanelStart[j])
            yP = -sinAlpha[j] * (xCP[i]-xPanelStart[j])   +   cosAlpha[j] * (yCP[i]-yPanelStart[j])
    
            # Calculate the induced velocity in panel coordinates
            vxP = inv_4pi * log( (xP*xP + yP*yP) / ((xP-x2)*(xP-x2) + yP*yP) )
            vyP = inv_2pi * ( atan2(yP,(xP-x2)) - atan2(yP, xP) )
    
            # Transform to global coordinates
            vx = cosAlpha[j]*vxP - sinAlpha[j]*vyP
            vy = sinAlpha[j]*vxP + cosAlpha[j]*vyP
            
            # Assemble the A matrix [influence Matrix]
            A[i,j] = vx*-sinAlpha[i] + vy*cosAlpha[i]
                
    # return the influence matrix                
    return A
    
    

@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
@cython.cdivision(True)    
def inducedVelocity(np.ndarray[DTYPE_t, ndim=1] sigma,
                    np.ndarray[DTYPE_t, ndim=1] xEval,
                    np.ndarray[DTYPE_t, ndim=1] yEval,
                    np.ndarray[DTYPE_t, ndim=1] xPanelStart,
                    np.ndarray[DTYPE_t, ndim=1] yPanelStart,
                    np.ndarray[DTYPE_t, ndim=1] xPanelEnd,
                    np.ndarray[DTYPE_t, ndim=1] yPanelEnd,
                    np.ndarray[DTYPE_t, ndim=1] cosAlpha,
                    np.ndarray[DTYPE_t, ndim=1] sinAlpha):
    r"""
    Calculates the induced velocity from :math:`\mathbf{M}` source panels on 
    :math:`\mathbf{n}` evaluations points.
    
    The source panels are represented by a *constant-strength source element*
    singularity element as given by:
    
        Katz, J., & Plotkin, A. (2001). Low-speed aerodynamics. Cambridge 
        University Press, [Chapter 11, 11.2.3]
        
    The velocity potential of the constant stength source element is given as:
    
    .. math::
    
        u_p = \frac{\sigma}{4\pi}\ln\frac{\left({x-x_1}\right)^2 + \left(y-y_1\right)^2}{\left(x-x_2\right)^2+ \left(y-y_1\right)^2},
        
        w_p = \frac{\sigma}{2\pi}\left[\tan^{-1}\frac{y-y_1}{x-x_2} - \tan^{-1}\frac{y-y_1}{x-x_1}\right]
        
    
    Usage
    -----
    .. code-block:: python
    
        vx, vy = inducedVelocity(gamma,xEval,yEval,xPanelStart,yPanelStart,
                                 xPanelEnd,yPanelEnd,cosAlpha,sinAlpha)


    Parameters 
    ---------- 
    gamma : ndarray (float64), shape (M,)
            the source panel strengths :math:`\gamma` of :math:`\mathbf{M}`
            panels.

    xEval : ndarray (float64), shape (nEval,)
            the :math:`x`-coordinates of the points where to evaluate the 
            induced velocities.
  
    yEval : ndarray (float64), shape (nEval,)
            the :math:`y`-coordinates of the points where to evaluate the 
            induced velocities
            
    xPanelStart : ndarray (float64), shape (M,)
                  the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
                  *starting point* in `global coordinate system`. This
                  variable is split so that it can be used for multiple 
                  geometries.

    yPanelStart : ndarray (float64), shape (M,)
                  the :math:`y`-coordinates of :math:`\mathbf{M}`panel corner
                  *starting point* in `global coordinate system`. This
                  variable is split so that it can be used for multiple 
                  geometries.

    xPanelEnd : ndarray (float64), shape (M,)
                the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
                *starting point* in `global coordinate system`. This
                variable is split so that it can be used for multiple 
                geometries.

    yPanelEnd : ndarray (float64), shape (M,)
                the :math:`y`-coordinates of :math:`\mathbf{M}`panel corner
                *starting point* in `global coordinate system`. This
                variable is split so that it can be used for multiple 
                geometries.

    cosAlpha : ndarray (float64), shape (M,)
               the :math:`\cos\alpha` of the panel :math:`\mathbf{M}`,
               where the angle :math:`\alpha` is w.r.t to the panel 
               :math:`\mathbf{x}^\prime`-axis and the global
               :math:`\mathbf{x}`-axis.

    sinAlpha : ndarray (float64), shape (M,)
               the :math:`\sin\alpha` of the panel :math:`\mathbf{M}`,
               where the angle :math:`\alpha` is w.r.t to the panel 
               :math:`\mathbf{x}^\prime`-axis and the global
               :math:`\mathbf{x}`-axis.
               
                    
    Returns
    -------
    vx : ndarray (float64), shape (nEval,)
         the :math:`x`-component of the induced velocities in each of the
         :math:`(\mathbf{xEval},\mathbf{yEval})` points.

    vy : ndarray (float64), shape (nEval,)
         the :math:`y`-component of the induced velocities in each of the
         :math:`(\mathbf{xEval},\mathbf{yEval})` points.
         

    :First Added:   2013-09-31
    :Last Modified: 2013-10-08
    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version
    
    """                        
    
    # Panel Parameters
    cdef int M      = sigma.shape[0] # number of panels
    cdef int nEval  = xEval.shape[0] # number of evaluation points

    # Constants
    cdef DTYPE_t inv_2pi = 1.0/(2.0*M_PI)
    cdef DTYPE_t inv_4pi = inv_2pi/2.0
    
    # panel parameters
    cdef DTYPE_t Lx,Ly,x2,xP,yP,vxP,vyP
    
    # Return variables
    cdef np.ndarray[DTYPE_t, ndim=1] vx = np.zeros(nEval, dtype=DTYPE) # induced velocity vx, vy
    cdef np.ndarray[DTYPE_t, ndim=1] vy = np.zeros(nEval, dtype=DTYPE) # induced velocity vx, vy
    
    # iteration parameters
    cdef int i,j    
    
    # Iterating through all the evaluation points
    for i in prange(nEval, nogil=True, schedule='static'):
        
        # Summing through all the panels
        for j in range(M):
            
            # Calculate the length of panel
            Lx = xPanelEnd[j] - xPanelStart[j]
            Ly = yPanelEnd[j] - yPanelStart[j]
            
            # Calculate the panel end coordinate in the panel coordinate system
            x2 = cosAlpha[j]*Lx + sinAlpha[j]*Ly
            
            # Locating the evaluation point in panel coordinate system
            xP =  cosAlpha[j]*(xEval[i]-xPanelStart[j]) + sinAlpha[j]*(yEval[i]-yPanelStart[j])
            yP = -sinAlpha[j]*(xEval[i]-xPanelStart[j]) + cosAlpha[j]*(yEval[i]-yPanelStart[j]) 

            # Calculate the induced velocity in panel coordinates
            vxP = (sigma[j]*inv_4pi)* log( (xP*xP + yP*yP) / ((xP-x2)*(xP-x2) + yP*yP) )            
            vyP = (sigma[j]*inv_2pi)* ( atan2(yP,(xP-x2)) - atan2(yP, xP) )
            
            # Summing the induced velocity
            vx[i] += cosAlpha[j]*vxP - sinAlpha[j]*vyP
            vy[i] += sinAlpha[j]*vxP + cosAlpha[j]*vyP
            
    # return the induced velocities            
    return vx,vy
