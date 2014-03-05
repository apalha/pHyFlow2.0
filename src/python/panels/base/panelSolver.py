"""
This module handles all panel method related operations in python
"""
# Copyright (C) 2013 Lento Manickathan                                                                          
#                                                                                                                                   
# This file is part of pHyFlow.                                                                                                      
#                                                                                                                                   
# pHyFlow is free software: you can redistribute it and/or modify                                                                    
# it under the terms of the GNU Lesser General Public License as published by                                                       
# the Free Software Foundation, either version 3 of the License, or                                                                 
# (at your option) any later version.                                                                                               
#                                                                                                                                   
# pHyFlow is distributed in the hope that it will be useful,                                                                         
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                    
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                      
# GNU Lesser General Public License for more details.                                                                               
#                                                                                                                                   
# You should have received a copy of the GNU Lesser General Public License                                                          
# along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.                                                                    
#                                                                                                                                   
# First added:  2013-09-09                                                                                                         
# Last changed: 2013-10-08


from pHyFlow import options
from pHyFlow.cpp.panels import kernels #,sourcePanel
#from pHyFlow.vortex.induced import velocity

# External packages
#import numpy
#from scipy.sparse.linalg import gmres, bicgstab


__all__ = ['influenceMatrix','inducedVelocity']


def influenceMatrix(xCP,yCP,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,
                    cosAlpha,sinAlpha,hardware=0,method=1):     
    r"""
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    Assemble the self-induction influence matrix for constant strength
    vortex panels. The operations are done using the appropriate hardware 
    and method.

    
    Usage
    -----
    .. code-block:: python
        
        A = influenceMatrix(xCP,yCP,xPanelStart,yPanelStart,xPanelEnd,
                            yPanelEnd,cosAlpha,sinAlpha,hardware,method)
                            
    
    Parameters
    ----------
    xCP : ndarray (float64), shape (M,)
          the :math:`x`-coordinates of :math:`\mathbf{M}`panel collocation 
          points in `global coordinate system`. For *vortex panels*, the
          collocation points should be at the cetner of the panel and slightly
          inside. 
           
    yCP : ndarray (float64), shape (M,)
          the :math:`y`-coordinates of :math:`\mathbf{M}`panel collocation 
          points in `global coordinate system`. For *vortex panels*, the
          collocation points should be at the cetner of the panel and slightly
          inside.                 
        
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
                  
    #    `panelKernel` : int (0 = constant source, 1 = constant vortex), optional, default = 1
    #                    the type of panel kernel that should be used for the 
    #                    computation.
                  
    `hardware` : int (0 = CPU, 1 = GPU), optional, default = 0
                 the hardware to use to compute the induced velocities
                  
    `method` : int (0 = FMM, 1 = Direct), optional, default = 1
               the method used to compute the induced velocities


    Returns
    -------
    A : ndarray (float64), shape (M,M)
        the assembled self-induction influence matrix :math:`\mathbf{A}` for 
        the panels.
        
        
    :First Added:   2013-09-31
    :Last Modified: 2013-11-26
    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version      
              
    """

    # compute optimized fullOption
    # the fullOtion is obtained as:
    #    100*method + 10*hardware
    # this means that if method = 2, hardware = 0we get as
    # full option:
    #    100*2 + 10*0 = 200 this reduces the options
    fullOption = 100*method + 10*hardware

    # Choose the computation type

    if fullOption == 100*options.DIRECT_METHOD + 10*options.CPU_HARDWARE:
        A = kernels.vortex._vortexPanel_Cython_cpu.assemble_influenceMatrix(xCP,yCP,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,cosAlpha,sinAlpha)

    elif fullOption == 100*options.DIRECT_METHOD + 10*options.GPU_HARDWARE:
       print "DIRECT + GPU + Constant strength VORTEX :: Not Implemented!"                      
                         
    elif fullOption < 100*options.DIRECT_METHOD:
        print "FMM :: Not Implemented!"

    # return the influence Matrix
    return A



def inducedVelocity(strength,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,
                    cosAlpha,sinAlpha,xEval,yEval,hardware=0,method=1):
                    
    r"""
    Compute the induced velocities of the panels on the collocation points.
    
    
    Usage
    -----
    .. code-block:: python
        
        A = inducedVelocity(strength,xPanelStart,yPanelStart,xPanelEnd,
                            yPanelEnd,cosAlpha,sinAlpha,xEval,yEval,
                            hardware,method)
                            

    Parameters
    ----------
    strength : ndarray (float64), shape (M,)
               the strengths :math:`\gamma` of :math:`\mathbf{M}` panels.
            
    xPanelStart : ndarray (float64), shape (M,)
                  the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
                  *starting point* in `local coordinate system`. This
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

    xEval : ndarray (float64), shape (nEval,)
            the :math:`x`-coordinates of the points where to evaluate the 
            induced velocities.
  
    yEval : ndarray (float64), shape (nEval,)
            the :math:`y`-coordinates of the points where to evaluate the 
            induced velocities
            
    #    `panelKernel` : int (0 = constant source, 1 = constant vortex), optional, default = 1
    #                    the type of panel kernel that should be used for the 
    #                    computation.
                  
    `hardware` : int (0 = CPU, 1 = GPU), optional, default = 0
                 the hardware to use to compute the induced velocities
                  
    `method` : int (0 = FMM, 1 = Direct), optional, default = 1
               the method used to compute the induced velocities           
            
                
    Returns
    -------
    vx : ndarray (float64), shape (nEval,)
         the :math:`x`-component of the induced velocities in each of the
         :math:`(\mathbf{xEval},\mathbf{yEval})` points.

    vy : ndarray (float64), shape (nEval,)
         the :math:`y`-component of the induced velocities in each of the
         :math:`(\mathbf{xEval},\mathbf{yEval})` points.    
    
    
    :First Added:   2013-09-31
    :Last Modified: 2013-11-26
    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version                            
                            
    """
        
    # compute optimized fullOption
    # the fullOtion is obtained as:
    #    100*method + 10*hardware + kernel
    # this means that if method = 2, hardware = 0  we get as
    # full option:
    #    100*2 + 10*0 = 200 this reduces the options
    fullOption = 100*method + 10*hardware

    # Choose the computation type
    
    if fullOption == 100*options.DIRECT_METHOD + 10*options.CPU_HARDWARE:
        vx,vy = kernels.vortex._vortexPanel_Cython_cpu.inducedVelocity(strength,xEval,yEval,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,cosAlpha,sinAlpha)

    elif fullOption == 100*options.DIRECT_METHOD + 10*options.GPU_HARDWARE:
       print "DIRECT + GPU + Constant strength VORTEX :: Not Implemented!"                      
                         
    elif fullOption < 100*options.DIRECT_METHOD:
        print "FMM :: Not Implemented!"        

    # return the induced velocities
    return vx,vy        
  
  

#def solve(xCP,yCP,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,cosAlpha,
#          sinAlpha,vxInduced,vyInduced,A=None,panelKernel=1,hardware=0,
#          method=1):
#    r"""
#    Solve for the panel strengths to ensure no-through flow condition.
#    Solved using a BIConjugate Gradient STABilized iteration to solve 
#    :math:`A\mathbf{x} = \mathbf{RHS}`
#    
#    
#    Usage
#    -----
#    .. code-block:: python
#        
#        A = solve(xCP,yCP,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,cosAlpha,
#                  sinAlpha,vxInduced,vyInduced,A,panelKernel,hardware,method)    
#    
#    
#    Parameters
#    ----------
#    xCP : ndarray (float64), shape (M,)
#          the :math:`x`-coordinates of :math:`\mathbf{M}`panel collocation 
#          points in `global coordinate system`. For *vortex panels*, the
#          collocation points should be at the cetner of the panel and slightly
#          inside.
#
#    yCP : ndarray (float64), shape (M,)
#          the :math:`y`-coordinates of :math:`\mathbf{M}`panel collocation 
#          points in `global coordinate system`. For *vortex panels*, the
#          collocation points should be at the cetner of the panel and slightly
#          inside.        
#        
#    xPanelStart : ndarray (float64), shape (M,)
#                  the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
#                  *starting point* in `global coordinate system`. This
#                  variable is split so that it can be used for multiple 
#                  geometries.
#
#    yPanelStart : ndarray (float64), shape (M,)
#                  the :math:`y`-coordinates of :math:`\mathbf{M}`panel corner
#                  *starting point* in `global coordinate system`. This
#                  variable is split so that it can be used for multiple 
#                  geometries.
#
#    xPanelEnd : ndarray (float64), shape (M,)
#                the :math:`x`-coordinates of :math:`\mathbf{M}`panel corner
#                *starting point* in `global coordinate system`. This
#                variable is split so that it can be used for multiple 
#                geometries.
#
#    yPanelEnd : ndarray (float64), shape (M,)
#                the :math:`y`-coordinates of :math:`\mathbf{M}`panel corner
#                *starting point* in `global coordinate system`. This
#                variable is split so that it can be used for multiple 
#                geometries.
#                
#    cosAlpha : ndarray (float64), shape (M,)
#               the :math:`\cos\alpha` of the panel :math:`\mathbf{M}`,
#               where the angle :math:`\alpha` is w.r.t to the panel 
#               :math:`\mathbf{x}^\prime`-axis and the global
#               :math:`\mathbf{x}`-axis.
#
#    sinAlpha : ndarray (float64), shape (M,)
#               the :math:`\sin\alpha` of the panel :math:`\mathbf{M}`,
#               where the angle :math:`\alpha` is w.r.t to the panel 
#               :math:`\mathbf{x}^\prime`-axis and the global
#               :math:`\mathbf{x}`-axis.
#               
#    vxInduced : ndarray (float64), shape (M,)
#                the :math:`x`-component of the induced velocities on each of 
#                the panel :math:`(\mathbf{x_{CP}},\mathbf{y_{CP}})` 
#                collocation points.
#         
#    vyInduced : ndarray (float64), shape (M,)
#                the :math:`y`-component of the induced velocities on each of 
#                the panel :math:`(\mathbf{x_{CP}},\mathbf{y_{CP}})` 
#                collocation points.
#
#    `A` : ndarray (float64), shape (M,M), optional, default = None
#          the assembled self-induction influence matrix :math:`\mathbf{A}` for 
#          the panels. If ``A`` is not given, it will be *recomputed*.
#          
#    `panelKernel` : int (0 = constant source, 1 = constant vortex), optional, default = 1
#                    the type of panel kernel that should be used for the 
#                    computation.
#                  
#    `hardware` : int (0 = CPU, 1 = GPU), optional, default = 0
#                 the hardware to use to compute the induced velocities
#                  
#    `method` : int (0 = FMM, 1 = Direct), optional, default = 1
#               the method used to compute the induced velocities            
#    
#    
#    Returns
#    -------
#    strength : ndarray (float64), shape (M,)
#               the strengths :math:`\gamma` of :math:`\mathbf{M}` panels.   
#    
#    
#    :First Added:   2013-09-31
#    :Last Modified: 2013-10-10
#    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
#    :Licence:       GNU GPL version 3 or any later version       
#    
#    """
#
#    # References:
#    # normX, normY = -sinAlpha, cosAlpha
#    # tangX, tangY = cosAlpha, sinAlpha
#
#    # Calculate the RHS
#
#    # Constant source panels
#    if panelKernel == options.PANELKERNEL_CONSTANTSOURCE:    
#        
#        # Solve for no-through condition - NORMAL Velocity at the panel
#        RHS = - (vxInduced*-sinAlpha + vyInduced*cosAlpha)
#
#    # Constant vortex panel        
#    elif panelKernel == options.PANELKERNEL_CONSTANTVORTEX:
#        
#        # Solve for no-slip velocity - TANGENT velocity at the panel
#        RHS = - (vxInduced*cosAlpha + vyInduced*sinAlpha)       
#
#    # others        
#    else:
#        raise RuntimeError('Not Implemented!')
#    
#    # Assemble the influence Matrix A
#    if A == None:
#        A = influenceMatrix(xCP,yCP,xPanelStart,yPanelStart,xPanelEnd,yPanelEnd,
#                            cosAlpha,sinAlpha,panelKernel,hardware,method)
#
#    # Solve for the panel strenths
#    # Note: dense 'A' matrix as an input is also valid (as in our case)
#
#    # using LAPACK routines, compute the "exact" solution of well-determined full-rank lin. matrix
#    #strength = numpy.linalg.solve(A,RHS) 
#    
#    # Use Generalized Minimal RESidual iteration to solve A x = b.
#    #strength = gmres(A,RHS,tol=1e-12)[0]  # tol = 1e-05 
#
#    # Use BIConjugate Gradient STABilized iteration to solve A x = b
#    strength = bicgstab(A,RHS,tol=1e-12)[0]
#    
#    # return the panel strength
#    return strength
    
