#-*- coding: utf-8 -*-
__doc__ = r"""

PANELS
======
The main class for the panel method. 

Description
-----------
This is the class for solving the panel method problem. The panel method
problem is an important part of the viscous vortex particle method. In order
to define the vorticity generation at the solid boundary, the panel method is 
need. 

As Cottet and Koumoutsakos have explained [1]_, the vortex sheet is used
to enforce the no-through-flow boundary condtion. The problem is solved using
panel method, where the body is discretized into :math:`\mathbf{M}` vortex 
panels:

.. math::
    
    \mathbf{A} \cdot \mathbf{\gamma} = \mathbf{RHS}
    
The vortex panel is based on the formulation of Katz [2]_. 


Methodology
-----------
The inviscid boundary problem is solved as follows

[1]. Setup the panel method problem by defining the body panels. To define
     the panel body (or multiple bodies), we require a list of parameters:
     
         xPanel,yPanel  : the panel end points.
         cmGlobal       : the global position of the panel geometry
         thetaLocal     : the local rotation of the panel geometry
         dPanel         : the off-set of the panel collocation points from
                          the mid-point of the panel. For the vortex panels,
                          the collocation points has to be placed slightly
                          inside.
                          
[2]. Assemble the inter-induction matrix :math:`\mathbf{A}`. This will have to
     be re-assembled if the body is in motion.
     
[3]. Solve for the panel vortex sheet strength :math:`\gamma`, such that no-slip
     boundary condition is satisfied at the panel collocation points
     :math:`\mathbf{x}_{cp},\mathbf{y}_{cp}`.
     
[4]. (Optional) If the body moves, update the global position and local rotation
     of the panel body. With this, the inter-induction matrix will also be 
     re-assembled.     
     
With the new panel strengths, the induced velocity of the panels can be
calculated. This is important in order to couple the vortex particles and
the no-through b.c.

Implemented Algorithms
----------------------

If the problem is fixed:
    - Pre-compute LU decomposition of the influence matrix :math:`\mathbf{A}`
    
If the problem is moving:
    - GMRES : Generalized Minimal RESidual iteration method
    - BICGSTAG : BIConjugate Gradient STABilized iteration method
    - direct : direct solver.

Reference
---------
.. [1] Cottet, G. H., & Koumoutsakos, P. D. (2000). Vortex Methods: Theory and 
       Practice. Measurement Science and Technology (Vol. 12, pp. 354–354). 
       Cambridge University Press. doi:10.1088/0957-0233/12/3/704
       
.. [2] Katz, J., & Plotkin, A. (2001). Low-speed aerodynamics.

:First Added:   2013-11-19
:Last Modified: 2014-02-20
:Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
:License:       GNU GPL version 3 or any later version
"""

#   Copyright (C) 2014 Lento Manickathan                                                                         
#   
#   This file is part of pHyFlow.                                                                                                      
#   
#   pHyFlow is free software: you can redistribute it and/or modify                                                                    
#   it under the terms of the GNU Lesser General Public License as published by                                                       
#   the Free Software Foundation, either version 3 of the License, or                                                                 
#   (at your option) any later version.                                                                                               
#   
#   pHyFlow is distributed in the hope that it will be useful,                                                                         
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                    
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                                      
#   GNU Lesser General Public License for more details.                                                                               
#      
#   You should have received a copy of the GNU Lesser General Public License                                                          
#   along with pHyFlow. If not, see <http://www.gnu.org/licenses/>.       

__all__ = ['Panels']


# External packages
import numpy as _numpy
import scipy.linalg as _splin
from scipy.sparse.linalg import gmres as _gmres
from scipy.sparse.linalg import bicgstab as _bicgstab
#import types as _types

# Import pHyFlow options
from pHyFlow.aux.customDecorators import simpleGetProperty
from pHyFlow.panels.base import panelSolver as _panelSolver
from pHyFlow.panels import panelOptions as _panelOptions

class Panels(object):
    r"""
    Panel method solver for defining the no-through-flow boundary condition of
    the vortex method.
    
    Usage
    -----
    .. code-block:: python
    
        panelBody = Panels(geometries, panelKernel='csv', problemType='moving',
                           velCompParams={'method':'direct', 'hardware':'cpu'},
                           solverCompParams={'method':'bicgstab', 'tol':1e-12,
                                             'assemble':'all'})
       
    Parameters
    ----------
    geometries : dict of dict, nBodies of dict
                 the dictionary containing dictionary of body name and its
                 parameters.
                 
                 Format : {<'geometryKey'>: <geometryDataName>}
                 
                 'geometryKey' : str, nBodies of key
                                the name (or key) of the body.
                               
                 geometryData : dict of data, {'xPanel','yPanel','cmGlobal','thetaLocal','dPanel'}
                                a dictionary containing all the parameters of 
                                the geometry (i.d. by geometryKey)
                                
                                'xPanel' : numpy.ndarray(float64), (nPanels+1,)
                                           the :math:`x`-coordinate of the :math:`\mathbf{M}` panel
                                           corners. The coordinates are defined in the local
                                           coordinates system with the reference point being (0,0).
                     
                                           * Note: Close the loop 
                     
                                'yPanel' : numpy.ndarray(float64), (nPanels+1,)
                                           the :math:`y`-coordinate of the :math:`\mathbf{M}` panel
                                           corners. The coordinates are defined in the local
                                           coordinates system with the reference point being (0,0).
                                             
                                           * Note: Close the loop
                     
                                'cmGlobal' : numpy.ndarray(float64), shape (2,) 
                                             the :math:`x,y` global coordinates of the local panel
                                             body reference point. For solving the problem, the panel
                                             body will be displaced according the global displacement
                                             vector **cmGlobal**.
                       
                                'thetaLocal' : float, unit (rad)
                                               the local rotation angle :math:`\theta` w.r.t to 
                                               the local coordinate system. The rotational will be 
                                               performed around the local reference point (0,0), i.e
                                               around the global cm point **cmGlobal**. 
                                               * Note: Positive rotation is anti-clockwise.
                                               
                                 'dPanel' : float
                                            the off-set of the panel collocation point from the panel 
                                            mid-point.                                               
                      
                  example:    
                  .. code-block:: python                        
                  
                      geometries = {'airfoil': {'xPanel' : np.array(-1.,1.,1.,1.),
                                                'yPanel' : np.array(-1.,1.,1.,1.),
                                                'cmGlobal': np.array(-1.,1.),
                                                'thetaLocal': 0.,
                                                'dPanel': np.spacing(100)},
                                    'cylinder': {'xPanel' : np.array(-1.,1.,1.,1.),
                                                 'yPanel' : np.array(-1.,1.,1.,1.),
                                                 'cmGlobal': np.array(-1.,1.),
                                                 'thetaLocal': 0.,
                                                 'dPanel': np.spacing(100)}}
                
                
    panelKernel : str, optional
                  A string defining panel kernel type:
                  'csv' : [default] Constant-strength vortex panel
                  'css' : Constant-strength source panel
                      
                  * Note: options are available in .. py:module:: pHyFlow.panel.panelOptions
              
   problemType : str, optional
                 A string defining the panel problem is of a moving type or of 
                 a fixed type.
                 'moving' : [default] The panel geometries change during each
                             iteration. Therefore, new coordinates and new
                             inter-induction matrix has to be recalculated.
                             * Note: LU Factorization algorithm is not possible.
                 'fixed' : The panel body is not moving and therefore the 
                           inter-induction matrix only needed to be computed
                           once. The LU factorization is very efficient for this
                           type of problem.
   
   velCompParams : dict, optional
                   A dictionary containing the velocity computation parameters,
                   method and hardware.
                   
                   'method' : str
                              the method to compute the induced velocities.
                              'direct' : [default] direct calculation with
                                         :math:`\mathcal{O}\left(n^2\right)`
                                         
                   'hardware' : str
                                the hardware for computing the induced velocity
                                of the panel.
                                'cpu' : [default] use the cpu.
                                      
                    * Note: options are available in .. py:module:: pHyFlow.panel.panelOptions

    solverCompParams : dict, optional
                       the dictionary containing solver computation parameters
                       
                       'method' : str
                                  the method to solve the problem.
                                  'bicgstab' : [default] BIConjugate Gradient 
                                               STABilized iteration method
                                  'gmres' : Generalized Minimal RESidual 
                                            iteration method
                                  'direct' : Direct solving.
                                  
                       'tol': str
                               the tolerance of the absolute and relative
                               residual for the iterative solvers: ('bicgstab','gmres')
                               'tol' : [default] 1e-12
                               
                       'assemble': str
                                   the type of assembling for the inter-
                                   induction matrix. 
                                   'all' : [default] all of the inter-induction
                                           matrix will be recalculated and
                                           assembled.
                                   * Note: Smart assembling is not implemented.
                           
                       * Note: options are available in .. py:module:: pHyFlow.panel.panelOptions           
   
    Attributes
    ----------
    A
    cmGlobal
    deltaT
    geometryKeys
    nBodies    
    norm
    normCat    
    nPanels
    nPanelsTotal
    panelKernel
    problemType
    solverCompParams
    sPanel
    t        
    tang
    tangCat    
    thetaLocal
    tStep
    velCompParams
    xyCPGlobal
    xyCPGlobalCat
    xyPanelGlobal
    xyPanelGlobalCat
    xyPanelLocal
    
    __A : numpy.ndarray(float64), shape (nPanelsTotal, nPanelsTotal)
          the inter-induction matrix :math:`\mathbf{A}`, the LHS of the problem.
                
    __cosSinAlpha : numpy.ndarray(float64), shape(2, nPanelsTotal)
                    the :math:`\cos\alpha` and :math:`\sin\alpha` of the panel
                    :math:`\mathbf{M}`, where the angle :math:`\alpha` is w.r.t 
                    to the panel :math:`\mathbf{x}^\prime`-axis and the local
                    :math:`\mathbf{x}`-axis.
                    
    __deltaT : float
               the current time step size
                    
    __geometries : dict of dict, nBodies of dict
                   the dictionary containing dictionary of body name and its
                   parameters.  
                   
                   {'xPanel','yPanel','cmGlobal','thetaLocal','dPanel'}
                                   
    __index : numpy.ndarray(float64), shape(nBodies+1,)
              the index of panel start point and end point in concatenated
              numpy array containing all the panel coordinates.
              
    __LU_FACTOR : tuple of (numpy.array(float64), shape(nPanelsTotal, nPanelsTotal) (1) and
                  numpy.array(float64), shape (nPanelsTotal,) (2).
                  The output of the computed pivoted LU decomposition of a 
                  inter-induction matrix. It is used to solver the system using
                  LU factorization.
                  1. The Matrix containing the U in its upper triange, L in its
                     lower triangle. 
                  2. The pivot indices represeting the permutation matrix P.
                  
    __nBodies : float
                the number of panel bodies.
    
    __norm : numpy.ndarray(float64), shape (2, nPanelsTotal)
             the :math:`x,y` normal vector of each panel.
             
    __nPanels : numpy.ndarray(float64), shape (nBodies,)
                the number of panels in each body.
                
    __nPanelsTotal: float
                    the total number of panels
                    
    __panelKernel : str
                    A string defining panel kernel type:
                    'csv' : [default] Constant-strength vortex panel
                    'css' : Constant-strength source panel
                      
    __problemType : str
                    A string defining the panel problem is of a moving type or
                    of a fixed type.
                    'moving' : [default] The panel geometries change during each
                               iteration. Therefore, new coordinates and new
                               inter-induction matrix has to be recalculated.
                               * Note: LU Factorization algorithm is not possible.
                    'fixed' : The panel body is not moving and therefore the 
                              inter-induction matrix only needed to be computed
                              once. The LU factorization is very efficient for this
                              type of problem.
                              
    __rotMat : list of numpy.ndarray(float64), N list of shape(2,2)
               the rotation matrix for rotating the local panel coordiantes
               with **thetaLocal** in the anti-clockwise direction.
               
    __solverCompParams : dict, optional
                         the dictionary containing solver computation parameters
                         'method' : str
                                    the method to solve the problem.
                                    'bicgstab' : [default] BIConjugate Gradient 
                                                 STABilized iteration method
                                    'gmres' : Generalized Minimal RESidual 
                                              iteration method
                                    'direct' : Direct solving.
                         'tol': str
                                the tolerance of the absolute and relative
                                residual for the iterative solvers: ('bicgstab','gmres')
                                'tol' : [default] 1e-12
                         'assemble': str
                                     the type of assembling for the inter-
                                     induction matrix. 
                                     'all' : [default] all of the inter-induction
                                             matrix will be recalculated and
                                             assembled.
                                     * Note: Smart assembling is not implemented.
                                     
    __sPanel : numpy.ndarray(float64), shape (nPanelsTotal,)
               the vortex sheet strengths :math:`\gamma` of :math:`\mathbf{M}`
               panels.  
               
     _t : float
          the current time of the simulation
          
    __tang : numpy.ndarray(float64), shape (2, nPanelsTotal)
             the :math:`x,y` tangent vector of each panel.
             
    __thetaLocal : list of float, N list of float
                   the local rotation angle :math:`\theta` w.r.t to 
                   the local coordinate system. The rotational will be 
                   performed around the local reference point (0,0), i.e
                   around the global cm point **cmGlobal**. 
                   * Note: Positive rotation is anti-clockwise.
                   
    __tStep : float
              the current step of the simulation
              
    __velCompParams : dict, optional
                      A dictionary containing the velocity computation parameters,
                      method and hardware.
                      'method' : str
                                 the method to compute the induced velocities.
                                 'direct' : [default] direct calculation with
                                            :math:`\mathcal{O}\left(n^2\right)`
                      'hardware' : str
                                   the hardware for computing the induced velocity
                                   of the panel.
                                   'cpu' : [default] use the cpu.
                      * Note: options are available in .. py:module:: pHyFlow.panel.panelOptions
               
    __xyCP_global : numpy.ndarray(float64), shape (2, nPanelsTotal)
                    the :math:`x,y` global coordinate of the collocation point.
                    
    __xyPanelEnd__global : numpy.ndarray(float64), shape (2, nPanelsTotal)
                           the :math:`x,y` global coordinate of the panel
                           starting points.
                           
    __xyPanelStart_global: numpy.ndarray(float64), shape (2, nPanelsTotal)
                           the :math:`x,y` global coordinate of the panel
                           end points.

    Methods
    -------
    solve(vxExternel,vyExternel)
        solve the panel strength to satisfy no-slip b.c.
        
    updateBody(thetaLocals, cmGlobals)
        update all the panel body coordinates.

    evaluateVelocity(xTarget,yTarget)
        evaluate the induced velocity due to the panels.
        
    save(dir)
        save all the ..py:class::`panels` data to a file.
        
    :First Added:   2013-11-21
    :Last Modified: 2014-02-25
    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version               
                  
    """
    """
    Revisions
    ---------
    2013-11-22, Lento Manickathan
        - Init using data file fixed.
        
    2014-02-20, Lento Manickathan
        - Restructured
        - Added attributes.

    2014-02-25, Lento Manickathan
        - Restructured
        - implemented geometry dictionary storage and retrival
        - Added attributes.

    """
    
    def __init__(self,geometries,panelKernel=_panelOptions.PANEL_KERNEL['default'],
                 problemType=_panelOptions.PROBLEM_TYPE['default'],
                 velCompParams={'method'  : _panelOptions.VELOCITY_COMPUTATION_METHOD['default'],
                                'hardware': _panelOptions.VELOCITY_COMPUTATION_HARDWARE['default']},
                 solverCompParams={'method' : _panelOptions.SOLVER_COMPUTATION_METHOD['default'],
                                   'tol'    : _panelOptions.SOLVER_COMPUTATION_TOL['default'],
                                   'assemble' :_panelOptions.SOLVER_COMPUTATION_ASSEMBLE['default']}):

        #----------------------------------------------------------------------
        # Check input parameters

        # Set the geometry(ies) datas.
        self.__set('geometries',geometries)
        
        # Set velocity computation parameters
        self.__set('velCompParams', velCompParams)
        
        # Set panel kernel type
        self.__set('panelKernel', panelKernel)
        
        # Set problem type
        self.__set('problemType', problemType)
        
        # Set solver computation parameters
        self.__set('solverCompParams',solverCompParams)
        #----------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------
        # Calculate parameters
        
        # Total number of panels
        self.__nPanelsTotal = _numpy.int64(self.__nPanels.sum())
        
        # Panel indices
        self.__index = _numpy.cumsum(_numpy.hstack((0,self.__nPanels)))
        
        # Update the rotation matrix
        self.__updateRotMat()
        
        # Initialize the global panel parameters       
        self.__norm = _numpy.zeros((2,self.__nPanelsTotal))
        self.__tang = _numpy.zeros((2,self.__nPanelsTotal))               
        self.__xyCP_global = _numpy.zeros((2,self.__nPanelsTotal))
        self.__cosSinAlpha = _numpy.zeros((2,self.__nPanelsTotal))
        self.__xyPanelStart_global = _numpy.zeros((2,self.__nPanelsTotal))
        self.__xyPanelEnd_global   = _numpy.zeros((2,self.__nPanelsTotal))
        
        # Determine the global coordinates (also, concatenate the data)
        self.__updateCoordinates()
        
        
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        # Assemble influence Matrix
        self.__A = _numpy.zeros((self.__nPanelsTotal,self.__nPanelsTotal))
        self.__assembleInfluenceMatrix()

        #----------------------------------------------------------------------
        # Define and initialize the time step counter and absolute time

        self.__tStep = 0 # Current Step
        self.__t = 0.0 # Current t
        
        #----------------------------------------------------------------------
   

    def updateBody(self,cmGlobalNew,thetaLocalNew):
        """
        Update all the panel body coordinates. This internally calculate the
        new global panel coordinates of `xPanel, yPanel, xCP, yCP` and will
        also rebuild the inter-induction influence matrix :math:`\mathbf{A}`
        
        Usage
        -----
        .. code-block :: python
        
            updateBody(cmGlobalNew,thetaLocalNew)
        
        Parameters
        ----------
        cmGlobalNew : dict of cmGlobalNew of each body, nBodies of keys
                      The dictionary containing the geometry key and it's new
                      cmGlobal.
                      Format : {'<geometryKey>': cmGlobalNew}
                   
                      'geometryKey' : str, nBodies of key
                                      the geometryKey is used to identify the body
                                      inside the panel class. The geometryKey
                                      can be found in **geometryKeys**.
                                 
                      cmGlobalNew : numpy.ndarray(float64), shape (2,)
                                    the new :math:`x,y` global coordinates of the
                                    local panel body reference point. For solving
                                    the problem, the panel body will be displaced
                                    according the global displacement vector **cmGlobalNew**.
                                   
        thetaLocalNew : dict of thetaLocalNew of each body, nBodies of keys
                        The dictionary containing the geometry key and it's new
                        thetaLocal.
                        Format: {'<geometryKey>': thetaLocalNew}
                        
                        'geometryKey' : str, nBodies of key
                                        the geometryKey is used to identify the body
                                        inside the panel class. The geometryKey
                                        can be found in **geometryKeys**.
                       
                                        
                        thetaLocalNew : float
                                        the new local rotation angle :math:`\theta`
                                        w.r.t to the local coordinate system. The
                                        rotational will be performed around the 
                                        local reference point (0,0), i.e around
                                        the global cm point **cmGlobal**. 
                                        * Note: Positive rotation is anti-clockwise.
        
        Returns
        -------
        None returned.
        
        Attribute
        ---------
        __A : numpy.ndarray(float64), shape (nPanelsTotal, nPanelsTotal)
              the inter-induction matrix :math:`\mathbf{A}`, the LHS of the problem.
                           
        __cosSinAlpha : numpy.ndarray(float64), shape(2, nPanelsTotal)
                        the :math:`\cos\alpha` and :math:`\sin\alpha` of the panel
                        :math:`\mathbf{M}`, where the angle :math:`\alpha` is w.r.t 
                        to the panel :math:`\mathbf{x}^\prime`-axis and the local
                        :math:`\mathbf{x}`-axis.
              
        __geometries : dict of dict, nBodies of dict
                       the dictionary containing dictionary of body name and its
                       parameters.  
                       
                       Updated parameters: {'cmGlobal','thetaLocal'}    
                       
        __norm : numpy.ndarray(float64), shape (2, nPanelsTotal)
                 the :math:`x,y` normal vector of each panel.
             
        __rotMat : list of numpy.ndarray(float64), nBodies list of shape(2,2)
                   the rotation matrix for rotating the local panel coordiantes
                   with **thetaLocal** in the anti-clockwise direction.
               
        __tang : numpy.ndarray(float64), shape (2, nPanelsTotal)
                 the :math:`x,y` tangent vector of each panel.
                
        __xyCP_global : numpy.ndarray(float64), shape (2, nPanelsTotal)
                        the :math:`x,y` global coordinate of the collocation point.
                    
        __xyPanelEnd__global : numpy.ndarray(float64), shape (2, nPanelsTotal)
                               the :math:`x,y` global coordinate of the panel
                               starting points.
                           
        __xyPanelStart_global: numpy.ndarray(float64), shape (2, nPanelsTotal)
                               the :math:`x,y` global coordinate of the panel
                               end points.
    
        :First Added:   2013-11-22
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version            
        
        """
        # Assign the new body location
        self.__set('thetaLocal', thetaLocalNew)
        self.__set('cmGlobal', cmGlobalNew)

        # Update the rotational matrix
        self.__updateRotMat()

        # Update the panel coordinates and angles
        self.__updateCoordinates()

        # Re-assemble the influence matrix A
        self.__assembleInfluenceMatrix()     


    def evaluateVelocity(self,xTarget, yTarget):
        """
        Evaluate the velocity field on the given target locations. The induced
        velocity of the vortex panel is defaulty calculed using CPU hardware 
        with direct computation. The CPU computation is parallelized using
        Cython prange.
        
        * Note: Can only be evaluated after solving the no-slip vortex sheets.
        
        Usage
        -----
        .. code-block :: python
            
            vx,vy = evaluateVelocity(xTarget,yTarget)
            
        Parameters
        ----------
        xTarget : numpy.ndarray(float64), shape (nTargets,)
                  the :math:`x`-coordinate of the target location, where the 
                  total velocity should be evaluated.
                  
        yTarget : numpy.ndarray(float64), shape (nTargets,)
                  the :math:`y`-coordinate of the target location, where the 
                  total velocity should be evaluated.
                  
        Returns
        -------
        vx : numpy.ndarray(float64), shape (nTargets,)
             the :math:`x`-component of the induced velocity at each of the
             **(xTarget,yTarget)** points.
             
        vy : numpy.ndarray(float64), shape (nTargets,)
             the :math:`x`-component of the induced velocity at each of the
             **(xTarget,yTarget)** points.
             
        Attributes
        ----------
        None changed.

        :First Added:   2013-11-22
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version   
        """
        
        # Determine the computation hardware
        if self.__velCompParams['hardware'] == 'cpu':
            hardware = _panelOptions._VELOCITY_COMPUTATION_HARDWARE_CPU
        else:
            raise NotImplementedError('GPU computation not implemented !')
        
        # Determine the computation method
        if self.__velCompParams['method'] == 'direct':
            method = _panelOptions._VELOCITY_COMPUTATION_METHOD_DIRECT
        else:
            raise NotImplementedError('FMM computation not implemented !')
            
        # Calculate the induced velocity (internal cython solver)
        vx,vy = _panelSolver.inducedVelocity(self.__sPanel,self.__xyPanelStart_global[0],
                                             self.__xyPanelStart_global[1],self.__xyPanelEnd_global[0],
                                             self.__xyPanelEnd_global[1],self.__cosSinAlpha[0],self.__cosSinAlpha[1],
                                             xTarget,yTarget,hardware,method)

        # return induced velocity
        return vx,vy
        
        
    def solve(self,vxExternel,vyExternel):
        """
        Solve the panel strengths satisfying the no-slip boundary condition
        for the vortex panel at the collocation point. With the parameters:
        **vxExternel**, **vyExternel**, new panel strenght is calculated
        such that the no-slip velocity is found at the collocation points.
        This ensure the inviscid boundary condition for the vortex method:
        no-through-flow at the body.
        
        Fixed:
            - the fixed problem is solved using LU factorization for efficiency.
            
        Moving:
            - the moving problem required constant recomputation of the LHS, so
              LU factorization is inefficient.
              
              - direct : first option is to solve the system of equation directly.
              - iterative: 'gmres' or 'bicgstab' can be used with 'tol' for
                           faster iterative solving.
        
        Usage
        -----
        .. code-block :: python
        
            solve(vxExternel,vyExternel)
            
        Parameters
        ----------
        vxExternel : numpy.ndarray(float64), shape (nPanelsTotal,)
                     the :math:`x`-component of the external velocity acting
                     of the :math:`\sum_i^N M_i` panel collocation points,
                     defining by **xCPGlobalCat** and **yCPGlobalCat**.
  
        vyExternel : numpy.ndarray(float64), shape (nPanelsTotal,)
                     the :math:`y-component of the external velocity acting
                     of the :math:`\sum_i^N M_i` panel collocation points,
                     defining by **xCPGlobalCat** and **yCPGlobalCat**.
          
        Returns
        -------
        None

        Attributes
        ----------
        __sPanel : numpy.ndarray(float64), shape (nPanelsTotal,)
                   the vortex sheet strengths :math:`\gamma` of :math:`\mathbf{M}`
                   panels.                   
               
        :First Added:   2013-11-21
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version                
        """
        
        # Determine the RHS
        RHS = (-vxExternel)*self.__tang[0] + (-vyExternel)*self.__tang[1] # vortex panel

        # Solve for panel strength
        if self.__problemType == 'fixed':
            # For the fixed problem, utilize LU factor algorithm
            self.__sPanel = _splin.lu_solve(self.__LU_FACTOR, RHS)
        
        # For moving problem
        elif self.__problemType == 'moving':
            
            # Use the iterative GMRES solver
            if self.__solverCompParams['method'] == 'gmres':   
                self.__sPanel = _gmres(self.__A,RHS,tol=self.__solverCompParams['tol'])[0]

            # Use the iterative bicgstab method                
            elif self.__solverCompParams['method'] == 'bicgstab':
                self.__sPanel = _bicgstab(self.__A,RHS,tol=self.__solverCompParams['tol'])[0]
            
            # Use the direct matrix solver.
            elif self.__solverCompParams['method'] == 'direct':
                self.__sPanel = _numpy.linalg.solve(self.__A,RHS)
     

    def _advanceTime(self,deltaTNew):
        r"""
        Function to advance the internal time clock.

        Usage
        -----
        .. code-block :: python

            _advanceTime(deltaT)

        Parameters
        ----------
        deltaT : float
                 the current time step size
                 
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        __deltaT : float
                   the current time step size
    
         _t : float
              the current time of the simulation
    
        __tStep : float
                  the current step of the simulation
        
        :First Added:   2014-02-25
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        # store the current time step size
        self.__deltaT = deltaTNew

        # advance tStep
        self.__tStep += 1

        # advance t
        self.__t = self.__deltaT*self.__tStep


    def __assembleInfluenceMatrix(self):
        r"""
        Assemble the inter-influence matrix :math:`\mathbf{A}`, i.e the LHS
        of the panel problem defined by Katz [1]_.
        
        Currently, the complete influence matrix is recalculated. However, in
        future, smart recalculation can be utilized where the section of the 
        LHS that changes needs to be recomputed.
        
        Fixed:
            - If the problem is fixed, LU_FACTOR is also calculated which
              can be used to solve the matrix faster.
        
        Moving:
            - LU Factorization cannot be used.
            
        .. [1] Katz, J., & Plotkin, A. (2001). Low-speed aerodynamics.
                
        Usage
        -----
        .. code-block:: python
        
            __assembleInfluenceMatrix()
            
        Parameters
        ----------
        None 
        
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        __A : numpy.ndarray(float64), shape (nPanelsTotal, nPanelsTotal)
              the inter-induction matrix :math:`\mathbf{A}`, the LHS of the problem.
              
        __LU_FACTOR : tuple of (numpy.array(float64), shape(nPanelsTotal, nPanelsTotal) (1) and
                      numpy.array(float64), shape (nPanelsTotal,) (2).
                      The output of the computed pivoted LU decomposition of a 
                      inter-induction matrix. It is used to solver the system using
                      LU factorization.
                      1. The Matrix containing the U in its upper triange, L in its
                         lower triangle. 
                      2. The pivot indices represeting the permutation matrix P.
        
        :First Added:   2013-11-25
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version          
        """
        
        # Determine the computation hardware
        if self.__velCompParams['hardware'] == 'cpu':
            hardware = _panelOptions._VELOCITY_COMPUTATION_HARDWARE_CPU
        else:
            raise NotImplementedError('GPU computation not implemented !')
        
        # Determine the computation method
        if self.__velCompParams['method'] == 'direct':
            method = _panelOptions._VELOCITY_COMPUTATION_METHOD_DIRECT
        else:
            raise NotImplementedError('FMM computation not implemented !')
            
        # Compute the complete A matrix
        if self.__solverCompParams['assemble'] == 'all':
            
            # Call the fast cythonized function to calculate A matrix
            self.__A = _panelSolver.influenceMatrix(self.__xyCP_global[0],
                                                   self.__xyCP_global[1],
                                                   self.__xyPanelStart_global[0],
                                                   self.__xyPanelStart_global[1],
                                                   self.__xyPanelEnd_global[0],
                                                   self.__xyPanelEnd_global[1],
                                                   self.__cosSinAlpha[0],
                                                   self.__cosSinAlpha[1],
                                                   hardware, method)
        
        else:
            raise NotImplementedError("Smart Assemble not implemented")

        # Compute pivoted LU decomposition of the Influence Matrix
        # only if problem is fixed.
        if self.__problemType == 'fixed':
            self.__LU_FACTOR = _splin.lu_factor(self.__A)
            
            
    def __updateCoordinates(self):
        """
        Function to update the coordinates. Calculates the global panel
        end, start and collocation points. After that, we can also calculate 
        the global normal and tangent vector of each panel.
        
        Usage
        -----
        .. code-block:: python
            
            __updateCoordinates()
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        __cosSinAlpha : numpy.ndarray(float64), shape(2, nPanelsTotal)
                        the :math:`\cos\alpha` and :math:`\sin\alpha` of the panel
                        :math:`\mathbf{M}`, where the angle :math:`\alpha` is w.r.t 
                        to the panel :math:`\mathbf{x}^\prime`-axis and the local
                        :math:`\mathbf{x}`-axis.
                    
        __norm : numpy.ndarray(float64), shape (2, nPanelsTotal)
                 the :math:`x,y` normal vector of each panel.
               
        __tang : numpy.ndarray(float64), shape (2, nPanelsTotal)
                 the :math:`x,y` tangent vector of each panel.
                
        __xyCP_global : numpy.ndarray(float64), shape (2, nPanelsTotal)
                        the :math:`x,y` global coordinate of the collocation point.
                    
        __xyPanelEnd__global : numpy.ndarray(float64), shape (2, nPanelsTotal)
                               the :math:`x,y` global coordinate of the panel
                               starting points.
                           
        __xyPanelStart_global: numpy.ndarray(float64), shape (2, nPanelsTotal)
                               the :math:`x,y` global coordinate of the panel
                               end points.
    
      
        :First Added:   2013-11-25
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version  
        
        """
            
        # Update panel coordinates
        for i, geometryData in enumerate(self.__geometries.itervalues()):
            
            # Start and end index
            iS, iE = self.__index[i], self.__index[i+1]
            
            # New global panel coordinates
            xyPanelNew = _numpy.dot(self.__rotMat[i], 
                                    _numpy.vstack((geometryData['xPanel'],
                                                   geometryData['yPanel']))) \
                               + geometryData['cmGlobal'].reshape(2,1)
            
            # Separate and concatenate the (start) and (end) points of the panels
            self.__xyPanelStart_global[:,iS:iE]  = xyPanelNew[:,:-1]
            self.__xyPanelEnd_global[:,iS:iE]    = xyPanelNew[:,1:]
            
            #------------------------------------------------------------------
            # Calculate the angle functions
            
            # Determine the length of the panels
            Lxy = self.__xyPanelEnd_global[:,iS:iE] - self.__xyPanelStart_global[:,iS:iE]
            L   = _numpy.sqrt(_numpy.sum(Lxy*Lxy,axis=0))
            
            # Panel Angles
            self.__cosSinAlpha[:,iS:iE] = Lxy / L # cosAlpha, sinAlpha
            #------------------------------------------------------------------
            
            # Unit Vectors (Tangent and Normal vector)
            self.__norm[:,iS:iE] = _numpy.vstack((-self.__cosSinAlpha[1,iS:iE], 
                                                  self.__cosSinAlpha[0,iS:iE])) 
                                                  
            self.__tang[:,iS:iE] = _numpy.vstack((self.__cosSinAlpha[0,iS:iE],
                                                  self.__cosSinAlpha[1,iS:iE]))
            
            # Determine the collocation points (at the new position)
            self.__xyCP_global[:,iS:iE] = 0.5 * (self.__xyPanelEnd_global[:,iS:iE] + \
                                                 self.__xyPanelStart_global[:,iS:iE])\
                                           - self.__norm[:,iS:iE]*geometryData['dPanel'] # shift collocation point by dPanel, in the noraml dir.
        
            

    def __updateRotMat(self):
        """
        Calculates the global transformation matrix for rotating the local
        panel coordinates with **thetaLocal** around **cmGlobal** of each
        body.
        
        Usage
        -----
        .. code-block:: python
        
            __updateRotMat()
            
        Parameters
        ----------
        None
        
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        __rotMat : list of numpy.ndarray(float64), nBodies list of shape(2,2)
                   the rotation matrix for rotating the local panel coordiantes
                   with **thetaLocal** in the anti-clockwise direction.
                   
        :First Added:   2013-11-25
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version                     
        """ 
        
        # The global rotational matrix
        self.__rotMat = [_numpy.array([[_numpy.cos(geometryData['thetaLocal']), -_numpy.sin(geometryData['thetaLocal'])],
                                       [_numpy.sin(geometryData['thetaLocal']), _numpy.cos(geometryData['thetaLocal'])]])
                         for geometryData in self.__geometries.itervalues()]                                         


    def __set(self,varName,var):
        """
        Function to check the shape, type of the input parameters and to
        set them to local attributes.
        
        varNames:
            - 'geometries'       
            - 'panelKernel'
            - 'problemType'
            - 'solverCompParams'
            - 'cmGlobal'
            - 'thetaLocal' 
            - 'velCompParams'
        
        Usage
        -----
        .. code-block:: python
        
            __set('varName',var)
            
        Paramaters
        ----------
        varName : str
                  the name of the local attribute to be checked and set.
                  
        var : any
              the data of the attribute that needs to be set.  
              
        Returns
        -------
        None returned.
        
        Attribute
        ---------
        ...
        
        :First Added:   2013-11-25
        :Last Modified: 2014-02-20
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version        
        """
        
        #----------------------------------------------------------------------
        # Load the panel data
        
        if varName == 'geometries':
            
            # Check geometries data type
            if type(var) != dict:
                raise TypeError("'geometries' must either be a dictionary of "\
                                "geometries, containing all the parameters of "\
                                "its body: ['xPanel','yPanel','cmGlobal','thetaLocal', 'dPanel'] "\
                                "in a dictionary format.")

            # Determine number of bodies
            self.__nBodies = len(var)

            # Parameters of a given geometry that needs to loaded and checked
            panelParmsList = ['xPanel','yPanel','cmGlobal','thetaLocal','dPanel']
                        
            # Allocate nPanel variable
            # variable storing the number of panels in each body
            self.__nPanels = _numpy.zeros(self.__nBodies,dtype=int)
            
            for i,(geometryName, geometryData) in enumerate(var.iteritems()):
                
                # Raise error if it not a dictionary                
                if type(geometryData) != dict:
                    raise TypeError("%geometries['%s'] must either be a dictionary, "\
                                    "containing the geometry's panel datas." % geometryName)
                                     
                # Check if all the keys are found inside 'geometryData'                                 
                for panelParameter in panelParmsList:
                    if panelParameter not in geometryData.keys():
                        raise ValueError("%s was not found in geometries['%s']."\
                                          % (panelParameter, geometryName))
            
                # Determine the number of panels of the current body
                # with M panels, we will have M+1 panel end points
                self.__nPanels[i] = geometryData['xPanel'].shape[0] - 1
                
                #--------------------------------------------------------------
                # Check 'xPanel'
                
                # Check data type
                if type(geometryData['xPanel']) != _numpy.ndarray or geometryData['xPanel'].dtype != float:
                    raise TypeError("'geometries['%s']['xPanel']' must be a "\
                                    "float64 numpy ndarray. It is %s."\
                                    % (geometryName,str(type(geometryData['xPanel']))))
                
                # Check shape of data
                if (geometryData['xPanel'].shape[0] - 1) != self.__nPanels[i]:
                    raise ValueError("geometries['%s']['xPanel'] must have shape " \
                                     "%s. It has shape %s." % (geometryName,
                                                               str((self.__nPanels[i],)),
                                                               str(((geometryData['xPanel'].shape[0] - 1),))))
                
                #--------------------------------------------------------------
                
                #--------------------------------------------------------------
                # Check 'yPanel'
                
                # Check data type
                if type(geometryData['yPanel']) != _numpy.ndarray or geometryData['yPanel'].dtype != float:
                    raise TypeError("'geometries['%s']['yPanel']' must be a "
                                    "float64 numpy ndarray. It is %s."\
                                     % (geometryName,str(type(geometryData['yPanel']))))
                    
                # Check shape of data    
                if (geometryData['yPanel'].shape[0] - 1) != self.__nPanels[i]:
                    raise ValueError("geometries['%s']['yPanel'] must have shape %s. "
                                     "It has shape %s." % (geometryName,
                                                           str((self.__nPanels[i],)),
                                                           str(((geometryData['yPanel'].shape[0] - 1),))))
                    
                #--------------------------------------------------------------                    
                    
                #--------------------------------------------------------------
                # Check 'cmGlobal'
                
                # Check data type
                if type(geometryData['cmGlobal']) != _numpy.ndarray or geometryData['cmGlobal'].dtype != float:
                    raise TypeError("'geometries['%s']['cmGlobal'] must be a float64 "\
                                    "numpy ndarray. It is a %s" % (geometryName,
                                                                   str(type(geometryData['cmGlobal']))))
                # Check data shape    
                if geometryData['cmGlobal'].shape[0] != 2:
                    raise ValueError("geometries['%s']['cmGlobal'] must have "\
                                     "shape (2,). It has shape %s."% (geometryName,
                                                                     str(geometryData['cmGlobal'].shape)))
                                                                     
                #--------------------------------------------------------------
                                                                     
                #--------------------------------------------------------------                                                                     
                # Check 'thetaLocal'
                                                                     
                # Check data type
                if type(geometryData['thetaLocal']) != float and type(geometryData['thetaLocal']) != _numpy.float64:
                    raise TypeError("geometries['%s']['thetaLocal'] must be a "\
                                    "float. It is a %s" % (geometryName,str(type(geometryData['thetaLocal']))))
                                    
                #--------------------------------------------------------------
                        
                #--------------------------------------------------------------                        
                # Check 'dPanel'
                
                # Check data type
                if type(geometryData['dPanel']) != float and type(geometryData['dPanel']) != _numpy.float64:
                    raise TypeError("geometries['%s']['dPanel'] must be a "\
                                    "float. It is a %s" % (geometryName,str(type(geometryData['dPanel']))))
                
                #--------------------------------------------------------------                        
                
            # Everything checked, simply save the data.
            self.__geometries = var
                

        #----------------------------------------------------------------------                 
        # Set the velocity computation parameters
        elif varName == 'velCompParams':
            
            # Check if the input method is an available options
            if var['method'] not in _panelOptions.VELOCITY_COMPUTATION_METHOD['available']:
                raise ValueError(r"velCompParams['method'] must be in [%s]. It is %s" % (_panelOptions.VELOCITY_COMPUTATION_METHOD['available'],var['method']))
            
            # Check if input harware is an available option
            if var['hardware'] not in _panelOptions.VELOCITY_COMPUTATION_HARDWARE['available']:
                raise ValueError(r"velCompParams['hardware'] must be in [%s]. It is %s" % (_panelOptions.VELOCITY_COMPUTATION_HARDWARE['available'],var['hardware']))   
            self.__velCompParams = var
        #----------------------------------------------------------------------            
            
        #----------------------------------------------------------------------            
        # Set solver computation parameters
        elif varName == 'solverCompParams':
            
            # Check if input method is an available option
            if var['method'] not in _panelOptions.SOLVER_COMPUTATION_METHOD['available']:
                raise ValueError(r"solverCompParams['method'] must be in [%s]. It is %s" % (_panelOptions.SOLVER_COMPUTATION_METHOD['available'],var['method']))
            
            # Check if the tolerance is float
            if type(var['tol']) != float and type(var['tol']) != _numpy.float64:
                    raise TypeError("solverCompParams['tol'] must be a numpy float. It is a %s" % (str(type(var['tol']))))

            # Check if input assemble is an available option.                    
            if var['assemble'] not in _panelOptions.SOLVER_COMPUTATION_ASSEMBLE['available']:
                raise ValueError(r"solverCompParams['assemble'] must be in [%s]. It is %s" % (_panelOptions.SOLVER_COMPUTATION_ASSEMBLE['available'],var['assemble']))
            self.__solverCompParams = var
            
        #----------------------------------------------------------------------            
        # Set the panel kernel
        elif varName == 'panelKernel':
            if var == 'csv':
                self.__panelKernel = var
            else:
                raise NotImplementedError('Use Constant Strength Vortex kernel: csv')

        #----------------------------------------------------------------------                
        # Set the problem type
        elif varName == 'problemType':
            if var not in _panelOptions.PROBLEM_TYPE['available']:
                raise ValueError('%s: %s unknown! Should be in [%s]' % (varName,var,str(_panelOptions.PROBLEM_TYPE['available'])))
            self.__problemType = var

        #----------------------------------------------------------------------                
        # Set thetaLocal
        elif varName == 'thetaLocal':
            
            # Check the parent type
            if type(var) != dict:
                raise TypeError("'thetaLocal' must either be a dictionary containing "\
                                "the following keys: %s" % str(self.geometryKeys))

            # Iterate through all the bodies                                
            for key in self.geometryKeys:
                if key not in var:
                    raise TypeError("'thetaLocal' must either be a dictionary containing "\
                                    "the following keys: %s" % str(self.geometryKeys))

                # Check the child data type
                if type(var[key]) != float and type(var[key]) != _numpy.float64:
                    raise TypeError("thetaLocal['%s'] must be a float. It is a %s"\
                                     % (key,type(var[key])))
                
                # Replace the thetaLocal
                self.__geometries[key]['thetaLocal'] = var[key]
     
        elif varName == 'cmGlobal':
            
            # Check the parent type
            if type(var) != dict:
                raise TypeError("'cmGlobal' must either be a dictionary containing "\
                                "the following keys: %s" % str(self.geometryKeys))
            
            # Iterate through all the bodies                                
            for key in self.geometryKeys:
                if key not in var:
                    raise TypeError("'thetaLocal' must either be a dictionary containing "\
                                    "the following keys: %s" % str(self.geometryKeys))

                # Check the child data type
                if type(var[key]) != _numpy.ndarray and var[key].dtype != float:
                    raise TypeError("cmGlobal['%s'] must be a float64 numpy ndarray."\
                                    "It is a %s" % (key,type(var[key])))
                                    
                # Check the child data shape
                if var[key].shape[0] != 2:
                    raise ValueError("cmGlobal['%s'] must have shape (2,). "\
                                     "It has shape %s."% (key,str(var[key].shape)))
                                     
                # Replace the thetaLocal
                self.__geometries[key]['cmGlobal'] = var[key]     
                
        
    #--------------------------------------------------------------------------        
    # All the attributes
        
    # Inter induction matrix
    @simpleGetProperty        
    def A(self):
        r"""
        A : numpy.ndarray(float64), shape (nPanelsTotal, nPanelsTotal)
            the inter-induction matrix :math:`\mathbf{A}`, the LHS of the problem.
        """
        return self.__A
        
    @simpleGetProperty
    def cmGlobal(self):
        r"""
        cmGlobal : list of numpy.ndarray(float64), nBodies list of shape (2,)
                   the global position vector for each of the :math:`\mathbf{N}`
                   body, refining the position of the local panel (0,0) in the
                   global coordinate system.
        """
        # Extract cmGlobal
        return [geometryData['cmGlobal'] for geometryData in self.__geometries.itervalues()]

    @simpleGetProperty
    def deltaT(self):
        r"""
        deltaT : float
                 the simulation time step size :math:`\Delta T`
        """
        return self.__deltaT
        
    @simpleGetProperty
    def geometryKeys(self):
        r"""
        geometryKeys : list of str
                       the key/names of the geometries in the panel problem. 
        """
        return self.__geometries.keys()
    
    @simpleGetProperty
    def nBodies(self):
        r"""
        nBodies : int
                  The total number of panel bodies/geometries
        """
        return self.__nBodies

    @simpleGetProperty
    def norm(self):
        r"""
        norm : list of numpy.ndarray(float64), nBodies of shape (2,nPanels_i)
               the global :math:`x,y`-component of the panel normal vector
               at each collocation points.
        """
        # Split norm concatenated array
        return [self.__norm[0,self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)],\
               [self.__norm[1,self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)]    
  
    @simpleGetProperty
    def normCat(self):
        r"""
        normCat : numpy.ndarray(float64), nBodies of shape (2,nPanelsTotal)
                  the global concatenated :math:`x,y`-component of the panel
                  normal vector at each collocation points.
        """
        return self.__norm
  
    @simpleGetProperty
    def nPanels(self):
        r"""
        nPanels : numpy.ndarray(int64), shape (nBodies,)
                  the number of panels in each body/geometry
        """
        return self.__nPanels
        
    @simpleGetProperty
    def nPanelsTotal(self):
        r"""
        nPanelsTotal : int
                       the total number of panels
        """
        return self.__nPanelsTotal

    @simpleGetProperty        
    def panelKernel(self):
        r"""
        panelKernel : str
                      A string defining panel kernel type:
                      
                      'csv' : [default] Constant-strength vortex panel
                      'css' : Constant-strength source panel
        """
        return self.__panelKernel
        
    @simpleGetProperty
    def problemType(self):
        r"""
        problemType : str
                      A string defining the panel problem is of a moving type or
                      of a fixed type.
                    
                      'moving' : [default] The panel geometries change during each
                                 iteration. Therefore, new coordinates and new
                                 inter-induction matrix has to be recalculated.
                               
                                 * Note: LU Factorization algorithm is not possible.
                      
                      'fixed' : The panel body is not moving and therefore the 
                                inter-induction matrix only needed to be computed
                                once. The LU factorization is very efficient for this
                                type of problem.
        """        
        return self.__problemType
        
    @simpleGetProperty
    def solverCompParams(self):
        r"""
        solverCompParams : dict
                           the dictionary containing solver computation parameters
                           
                           'method' : str
                                      the method to solve the problem.
                                      
                                      'bicgstab' : [default] BIConjugate Gradient 
                                                   STABilized iteration method
                                      'gmres' : Generalized Minimal RESidual 
                                                iteration method
                                      'direct' : Direct solving.
                           'tol' : str
                                   the tolerance of the absolute and relative
                                   residual for the iterative solvers: ('bicgstab','gmres')
                                   
                                   'tol' : [default] 1e-12
                         
                           'assemble' : str
                                        the type of assembling for the inter-
                                        induction matrix. 
                                     
                                        'all' : [default] all of the inter-induction
                                                matrix will be recalculated and
                                                assembled.
                                      
                                        * Note: Smart assembling is not implemented.        
        """
        return self.__solverCompParams
                            
    @simpleGetProperty
    def sPanel(self):
        r"""
        sPanel : numpy.ndarray(float64), shape (nPanelsTotal,)
        """
        # Split the sPanel data
        return [self.__sPanel[self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)]
     
    @simpleGetProperty
    def t(self):
        r"""
        t : float
            the current time of the simulation
        """
        return self.__t
        
    @simpleGetProperty
    def tang(self):
        r"""
        tang : list of numpy.ndarray(float64), nBodies of shape (2,nPanels_i)
               the global :math:`x,y`-component of the panel tangent vector
               at each collocation points.
        """
        # Split norm concatenated array
        return [self.__tang[0,self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)],\
               [self.__tang[1,self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)]    
  
    @simpleGetProperty
    def tangCat(self):
        r"""
        tangCat : numpy.ndarray(float64), nBodies of shape (2,nPanelsTotal)
                  the global concatenated :math:`x,y`-component of the panel
                  normal vector at each collocation points.
        """
        return self.__tang
        
    @simpleGetProperty
    def thetaLocal(self):
        r"""
        thetaLocal : list of float, nBodies of float
                     the local rotation angle :math:`\theta` w.r.t to the local
                     coordinate system. The rotational will be performed around
                     the local reference point (0,0), i.e around the global cm 
                     point **cmGlobal**. 
                                                   
                     * Note: Positive rotation is anti-clockwise.                              
        """
        # Extract data
        return [geometryData['thetaLocal'] for geometryData in self.__geometries.itervalues()]
        
    @simpleGetProperty
    def tStep(self):
        r"""
        tStep : float
                the current step of the simulation
        """
        return self.__tStep

    @simpleGetProperty
    def velCompParams(self):
        r"""
        velCompParams : dict, optional
                        A dictionary containing the velocity computation parameters,
                        method and hardware.
                      
                        'method' : str
                                   the method to compute the induced velocities.
                                   
                                   'direct' : [default] direct calculation with
                                              :math:`\mathcal{O}\left(n^2\right)`
                        
                        'hardware' : str
                                     the hardware for computing the induced velocity
                                     of the panel.
                                     
                                     'cpu' : [default] use the cpu.
        """
        return self.__velCompParams        
        
    @simpleGetProperty
    def xyCPGlobal(self):
        r"""
        xyCPGlobal : list of numpy.ndarray(float64), nBodies of shape (2,nPanels_i)
                     the global :math:`x,y`-coordinate of the panel collocation
                     points.
        """
        # Split global coordinate
        return [self.__xyCP_global[0,self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)],\
               [self.__xyCP_global[1,self.__index[i]:self.__index[i+1]] for i in range(self.__nBodies)]
                        
    @simpleGetProperty
    def xyCPGlobalCat(self):
        r"""
        xyCPGlobalCat : numpy.ndarray(float64), shape (2,nPanelsTotal)
                        the global concatenated :math:`x,y`-coordinate of the
                        panel collocation points.
        """                 
        return self.__xyCP_global  
        
    @simpleGetProperty
    def xyPanelGlobal(self):
        r"""
        xyPanelGlobal : list of numpy.ndarray(float64), nBodies of shape (2,nPanels_i+1)
                        the global :math:`x,y`-coordinate of the panel bodies.
                       
                        * Note: closed loop
        """
        return [_numpy.append(self.__xyPanelStart_global[0,self.__index[i]:self.__index[i+1]],
                              self.__xyPanelStart_global[0,self.__index[i]]) for i in range(self.__nBodies)],\
               [_numpy.append(self.__xyPanelStart_global[1,self.__index[i]:self.__index[i+1]],
                              self.__xyPanelStart_global[1,self.__index[i]]) for i in range(self.__nBodies)]               
                 
    
    @simpleGetProperty
    def xyPanelGlobalCat(self):
        r"""
        xyPanelGlobalCat : numpy.ndarray(float64), shape (2,nPanelsTotal+nBodies)
                           the global concatenated :math:`x,y`-coordinate of the
                           panel bodies.
                          
                           * Note: closed loop
        """
        # Initialize the variable
        xyPanelGlobalCat = _numpy.zeros((2,self.__nPanelsTotal + self.__nBodies))
        
        # Close the panel body (append start point)
        for i in range(self.__nBodies):
            iS,iE = self.__index[i] + i, self.__index[i+1] + (i+1)
            xyPanelGlobalCat[:,iS:iE] = _numpy.hstack((self.__xyPanelStart_global[:,iS:iE],
                                                       self.__xyPanelStart_global[:,iS].reshape(2,1)))
        # return xPanel Global concatenated data                                                   
        return xyPanelGlobalCat
  
    @simpleGetProperty
    def xyPanelLocal(self):
        r"""
        xyPanelLocal : list of numpy.ndarray(float64), nBodies of shape (2,nPanels_i+1) 
                       the local :math:`x,y`-coordinate of the panel bodies.
                      
                      * Note: closed loop.
        """
        return [geometryData['xPanel'] for geometryData in self.__geometries.itervalues()],\
               [geometryData['yPanel'] for geometryData in self.__geometries.itervalues()]
 
 
    #--------------------------------------------------------------------------


#    def save(self,fileName=None):
#        """
#        Save all the data of the ..py:class::`panels` class to a file. This
#        file can be used later to re-initalize the class.
#        
#        Usage
#        -----
#        ...
#        
#
#        Parameters
#        ----------
#        ...
#        
#
#        Assigns
#        -------
#        ...
#        
#        
#        Returns
#        -------
#        ...
#        
#        
#        :First Added:   2013-11-21
#        :Last Modified: 2013-11-22
#        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
#        :Licence:       GNU GPL version 3 or any later version                     
#        
#        """
#        
#        # Define the file store location            
#        if fileName is None:
#            fileName = './panelData'
#            
#        # Save data
#        _numpy.savez_compressed(fileName,xCP=self.__xCP,yCP=self.__yCP,
#                            xPanel=self.__xPanel,yPanel=self.__yPanel,
#                            cmGlobal=self.__cmGlobal,thetaLocal=self.__thetaLocal,
#                            nPanels=self.__nPanels,nBody=self.__nBodies)