__doc__ = """

Vortex particle solver master class. 

* Note: No Turbulence schemes and vorticity generation at the solid boundaries implemented! *

Description
-----------
Vortex particle solver [1]_ for the solution of 2D incompressible laminar fluid. As no
turbulence scheme is implemented, the problems should be **LAMINAR**. Solves 
the problem for a given initial vorticity distribution. The time stepping algorithm
is a viscous splitting one:

    1- advection step with no viscosity
    2- viscosity step
    
Vortex particles are redistributed into a regular grid at given time steps. This
redistribution can include viscosity modelling, [2]_.

Induced velocities are computed using a Fast Multipole Biot-Savart solver on a GPU [3]_.
    
    
Implemented Algorithms
----------------------

Time stepping:
    - Euler time stepping ['euler']
    - Runge-Kutta 4th order time stepping ['rk4']
    
Diffusion (viscosity):
    - Modified interpolation kernel, [2]_

Induced velocity:
    - 2D FMM, [3]_

References
----------
.. [1] Cottet, G-H., Koumoutsakos, P.D., Vortex Methods: Theory and practice,
        Cambridge University Press, 2000.
        
.. [2] Wee, D., Ghoniem, A., Modified interpolation kernels, for treating diffusion
        and remeshing in vortex methods, Journal of Computational Physics, 213(1),
        p. 239-263, 2006.
        
    [3] Engblom S., On well-separated sets and fast multipole methods,
        Applied Numerical Mathematics 61(10):1096--1102, 2011.


:First added:   2013-12-18  
:Last updated:  2013-12-21     
:Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
:License:       GNU GPL version 3 or any later version
      
"""

"""
Reviews:
-------
        1- Added method VortexBlobs.regrid                 (apalha, 2014-01-15)
        2- Added method VortexBlobs.populationControl      (apalha, 2014-01-20)
        3- Added method VortexBlobs.evaluateVelocity       (apalha, 2014-01-21)

"""


__all__ = ['VortexBlobs']

import numpy as _numpy
import types as __types
import inspect as __inspect

from pHyFlow import options # import options definitions for pHyFlow

# import required base functions from where to build up the class
from pHyFlow.vortex.base.induced import velocity as _base_velocity
from pHyFlow.vortex.base.induced import vorticity_blobs as _base_vorticity_blobs
from pHyFlow.vortex.base.induced import vorticity as _base_vorticity
from pHyFlow.vortex.base.regrid import Regrid as _base_redistribute
from pHyFlow.vortex.base.regrid import PopulationControl as _base_populationControl



default_params = {'stepRedistribution':1, 'integrationMethod':'rk4',\
                  'computationMethod':('fmm','gpu'), 'stepPopulationControl':1,\
                  'gThreshold':(1e-8,1e-8),'diffusion':{'method':'regrid_wee','c2':0.2}}

default_diffusion_params = {'method':'regrid_wee','c2':0.2}
default_time_integration_params = {'method':'rk4','deltaTc':0.1,\
                                   'stepRedistribution':1}
default_population_control_params = {'stepPopulationControl':1,'gThresholdLocal':1e-8,\
                                     'gThresholdGlobal':1e-8}
default_velocity_computation_params = {'method':'fmm','hardware':'gpu'}


class VortexBlobs(object):
    r"""
        Vortex blob solver for the Navier-Stokes equation. No solid boundaries
        are implemented.
        
        Usage
        -----
        .. code-block:: python
        
            VortexBlobs(wField,vInf,overlap,h,deltaTc,nu,parameters)
            
        Parameters
        ----------
        wField : tuple of three numpy.array(float64) (1) or tuple of a function
                 and two numpy.array(float64) (2)
                 (1) (`xBlobs`, `yBlobs`, `gBlob`) defines the initial distribution of
                     vortex blobs. `xBlobs` contains the x coordinates, `yBlobs`
                     contains the y coordinates and `gBlobs` contains the circulations.
                     Blobs are recursively redistributed into a regular mesh. This
                     mesh is assumed to have a particle located at (0,0). All
                     arrays must have shape (1,nBlobs).
                 (2) (`wExactFunction`, `xBounds`, `yBounds`) defines the initial
                     distribution of vorticity and the region of interest, in x
                     is the region between `xBounds` [0] and `xBounds` [1] and in y
                     is the region between `yBounds` [0] and `yBounds` [1]. It is
                     assumed that `xBounds` [1] > `xBounds` [0] and `yBounds` [1] > `yBounds` [0].
                     Vortex blobs are then generated inside this region with the
                     overlap `overlap` and spacing `h`.
        
        vInf : numpy.array(float)
               The free-stream velocity,  vx = `vInf` [0] and vy = `vInf` [1].
                        
        overlap : float64
                  The overlap ratio between neighboring blobs. It is related to
                  the core size of the blob, `sigma`, and to the spacing 'h' by
                  the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.
                  
        h : float64
            The size of the cell associated to the vortex blobs. Corresponds to
            the minimum spacing between the core of two neighboring cells. It is related to
            the core size of the blob, `sigma`, and to the spacing 'h' by
            the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`
            
        deltaTc : float64
                  The size of the convective time step :math:`\Delta t_{c}`.
                  
        nu : float64
             The fluid kinematic viscosity, used to calculate the diffusion coefficient
             :math:`c^{2}` and diffusion time step `deltaTd`, :math:`\Delta t_{d}`.
        
        parameters : dict, optional
                     A dictionary containing all additional parameters of the
                     simulation: `stepRedistribution`, `integrationMethod`,
                     `computationMethod`, `stepPopulationControl`, `gThreshold`.
                     
                     `stepRedistribution` (int) the redistribution step frequency.
                     
                     `integrationMethod` ('fe', 'rk4') the type of time integrator
                     used: 'fe' forward Euler, 'rk4' Runge-Kutta 4th order.
                     
                     'computationMethod' (tuple) refering with the type of Biot-Savart
                     solver ('direct', 'fmm') and the type of hardware to use
                     ('cpu','gpu').
                     
                     'stepPopulationControl' (int) the population control step frequency.
                     
                     'gThreshold' (tuple) with minimum value of circulation to
                     consider for each vortex blob and maximum value of total
                     circulation change by discarding particles.
                     
                     'diffusion' (dict) defining the diffusion method to use and
                     the parameters of the diffusion method. It has the key 'method'
                     that defines the method to use: 'regrid_wee' (based in [1]_).
                     Parameters for each different method:
                     'regrid_wee':
                         `c2` (float64), free parameter that specifies the stability
                                         of the diffusion scheme. Best values: [1/6,1/2].
                     
                     default: {`stepRedistribution`:1, `integrationMethod`:'rk4',
                     `computationMethod`:('fmm','gpu'), `stepPopulationControl`:1,
                     `gThreshold`:(1e-8,1e-8),'diffusion':{'method':'regrid_wee','c2',0.2}}
        
        Attributes
        ----------
        computationMethod
        deltaTc
        deltaTd
        diffusionParameters        
        g
        gThreshold
        h
        integrationMethod        
        nu
        num_blobs
        overlap
        parameters        
        sigma
        stepPopulationControl
        stepRedistribution
        vInf        
        x
        y
        __diffusionParameters : dict
                                Contains all parameters specific to the diffusion
                                method selected.
        __deltaTc : float64
                    The size of the convective time step :math:`\Delta t_{c}`.
        __deltaTd : float64
                    The size of the diffusive time step :math:`\Delta t_{d}`. It
                    is a integer multiple of the convective time step deltaTc.
        __g : numpy.array(float64) (1,num_blobs)
              The g circulation of the vortex blobs.
        __h : float64
              The size of the cell associated to the vortex blobs. Corresponds to
              the minimum spacing between the core of two neighboring cells. It is related to
              the core size of the blob, `__sigma`, and to the overlap '__overlap' by
              the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.
        __nu : float64
               The fluid kinematic viscosity, used to calculate the diffusion coefficient
               :math:`c^{2}` and diffusion time step `deltaTd`, :math:`\Delta t_{d}`.
        __overlap : float64
                    The overlap ratio between neighboring blobs. It is related to
                    the core size of the blob, `__sigma`, and to the spacing '__h' by
                    the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.
        __parameters : dict, optional
                     A dictionary containing all additional parameters of the
                     simulation: `stepRedistribution`, `integrationMethod`,
                     `computationMethod`, `stepPopulationControl`, `gThreshold`.
                     
                     `stepRedistribution` (int) the redistribution step frequency.
                     
                     `integrationMethod` ('fe', 'rk4') the type of time integrator
                     used: 'fe' forward Euler, 'rk4' Runge-Kutta 4th order.
                     
                     'computationMethod' (tuple) refering with the type of Biot-Savart
                     solver ('direct', 'fmm') and the type of hardware to use
                     ('cpu','gpu').
                     
                     'stepPopulationControl' (int) the population control step frequency.
                     
                     'gThreshold' (tuple) with minimum value of circulation to
                     consider for each vortex blob and maximum value of total
                     circulation change by discarding particles.
                     
                     'diffusion' (dict) defining the diffusion method to use and
                     the parameters of the diffusion method. It has the key 'method'
                     that defines the method to use: 'regrid_wee' (based in [1]_).
                     Parameters for each different method:
                     'regrid_wee':
                         `c2` (float64), free parameter that specifies the stability
                                         of the diffusion scheme. Best values: [1/6,1/2].
                                         
        __sigma : float64
                  The core size of the vortex blobs.
        __vInf : float64
                 The free stream velocity.          
        __x : numpy.array(float64) (1,num_blobs)
              The x coordinates of the vortex blobs.
        __y : numpy.array(float64) (1,num_blobs)
              The y coordinates of the vortex blobs.
        
        
        .. [1] Wee, D., Ahmed, F., Ghoniem, F., Modified interpolation
              kernels for treating diffusion and remeshing in vortex
              methods, Journal of Computational Physics, 213, 239--263, 2006.
              
              
        :First Added:   2013-12-18
        :Last Modified: 2014-01-15
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
    """
        
    """    
    Reviews:  Added method self.regrid (apalha, 2014-01-15)
              Added method self.populationControl (apalha, 2014-01-20)
              Added method self.evaluateVelocity (apalha, 2014-01-21)
        
    """ 
        
    def __init__(self,wField,vInf,overlap,h,deltaTc,nu,parameters=default_params):
        """    
        Reviews:  
        
        """ 
        
        # read the inflow input parameters
        self.__read_input_flow_parameters(vInf,overlap,h,deltaTc,nu)
            
        # read the computation input parameters
        self.__read_input_computation_parameters(parameters)    
        
        # read wField input
        self.__read_input_wField(wField)
        
        # define the time step counter
        self.__tStep = 0

    
    def evaluateVelocity(self,xEval,yEval,Ndirect=35,tol=1.0e-6,cutoff=None,hardware=0):
        r"""
        Computes the induced velocities generated by the vortex blobs. It can use
        a FMM approximation (faster, but less accurate) or direct calculation
        (slower, but machine precision accurate).
        
        Usage
        -----
        .. code-block self.evaluateVelocity(xEval,yEval,Ndirect=35,tol=1.0e-6,cutoff=None)
    
        Parameters
        ----------
            xEval : numpy.array(float64)
                    the x coordinates of the points where to evaluate the
                    induced velocities
                    (if value is None, xEval = self.x)
                    shape: (nEval,)
            yEval : numpy.array(float64)
                    the y coordinates of the points where to evaluate the
                    induced velocities
                    (if value is None, yEval = self.y)
                    shape: (nEval,)
            Ndirect : int (optional)
                      the number of neighbor blobs where to perform direct calculation
                      (default value is 35)
                      shape: single value
            tol : float64 (optional)
                  the tolerance (error) with which the induced velocities are computed
                  (default value is 1.0e-6, if tol=0.0 then direct calculation
                  is performed instead of FMM)
                  shape: single value
            cutoff : float64 (optional)
                     the radius over which the velocity field is approximately 1/r^2
                     (default value is None, which corresponds to 5.0*xopt)
                     shape: single value
            hardware : int (optional)
                       the hardware to use to compute the velocity: 0 for CPU,
                       1 for GPU.
                       (default value is 0, CPU)
                       shape: single value
                    
        Returns
        -------
            vx :: numpy.ndarray(float64)
                  the x component of the induced velocities in each of the 
                  (xEval,yEval) points
                  shape: (nEval,)
                  
            vy :: numpy.ndarray(float64)
                  the y component of the induced velocities in each of the 
                  (xEval,yEval) points
                  shape: (nEval,)
            
            
        Attributes
        ----------
            None changed            
                  
        :First Added:   2014-01-21
        :Last Modified: 2014-01-21
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """
    
        """    
        Reviews:  
        
        """
        
        # compute the velocity at the evaluation points
        # computes de velocity with fmm or direct depending on tol
        vx,vy = _base_velocity(self.x,self.y,self.g,self.sigma,\
                               xEval=xEval,yEval=yEval,hardware=hardware,method=0,\
                               Ndirect=Ndirect,tol=tol,cutoff=cutoff)

        return vx,vy
        
        
    def evaluateVorticity(self,xEval,yEval,hardware=0,blocksize=128,method=1):
        r"""
        Computes the induced velocities generated by the vortex blobs. It can use
        a FMM approximation (faster, but less accurate) or direct calculation
        (slower, but machine precision accurate).
        
        Usage
        -----
        .. code-block self.evaluateVelocity(xEval,yEval,Ndirect=35,tol=1.0e-6,cutoff=None)
    
        Parameters
        ----------
            xEval : numpy.array(float64)
                    the x coordinates of the points where to evaluate the
                    induced velocities
                    (if value is None, xEval = self.x)
                    shape: (nEval,)
            yEval : numpy.array(float64)
                    the y coordinates of the points where to evaluate the
                    induced velocities
                    (if value is None, yEval = self.y)
                    shape: (nEval,)
                    
            Ndirect :: int (optional)
                       the number of neighbor blobs where to perform direct calculation
                       (default value is 35)
                       shape: single value
            
            tol :: float64
                   the tolerance (error) with which the induced velocities are computed
                   (default value is 1.0e-6, if tol=0.0 then direct calculation
                   is performed instead of FMM)
                   shape: single value
            
            cutoff :: float64
                      the radius over which the velocity field is approximately 1/r^2
                      (default value is None, which corresponds to 5.0*xopt)
                      shape: single value
                    
        
        Returns
        -------
            vx :: numpy.ndarray(float64)
                  the x component of the induced velocities in each of the 
                  (xEval,yEval) points
                  shape: (nEval,)
                  
            vy :: numpy.ndarray(float64)
                  the y component of the induced velocities in each of the 
                  (xEval,yEval) points
                  shape: (nEval,)
            
            
        Attributes
        ----------
            None changed            
                  
        :First Added:   2014-01-21
        :Last Modified: 2014-01-21
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """
    
        """    
        Reviews:  
        
        """
        
        pass
    
    
    def evolution(self):
        pass
    
        
    def populationControl(self):
        r"""
        PopulationControl controls the number of particles by discarding
        particles with circulation higher than self.gThreshold[0] (the minimum value
        of circulation to consider). If total circulation of the flagged blobs
        (flaggedCirculation) is more than self.gThreshold[1] (the maximum variation
        in total circulation due to population control) then self.gThreshold[1] in
        this population control step is multiplied by 10 and the test is
        performed again.
        
        In the end it is insured that the total circulation did not change more
        than self.gThreshold[1].
        
        Usage
        -----
        .. code-block self.populationControl()
    
        Parameters
        ----------
    
        Attributes
        ----------
            x
            y
            g
                  
        :First Added:   2014-01-20
        :Last Modified: 2014-01-20
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """
    
        """    
        Reviews:  
        
        """
        
        # perform the population control and update the blobs witht the new blobs
        self.__x,self.__y,self.__g = _base_populationControl(self.x,self.y,self.g,self.gThreshold[0],self.gThreshold[1])
        
        
    def redistribute(self,interpKernel=0):
        r"""
        Redistribute vortex blobs over evenly spaced grid points with separation
        h = sigma * overlap.
        The interpolation kernel can be one of the following:
           - M4' (0)
        
        Usage
        -----
        .. code-block self.redistribute(interpKernel=0)
            
        Parameters
        ----------
        interpKernel :: int64, optional
                        the interpolating kernel using for remeshing
                        0 : M4' (default)
             
        Attributes
        ----------
            x
            y
            g
                  
                  
        :First Added:   2014-01-15
        :Last Modified: 2014-01-15
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """
    
        """    
        Reviews:  
        
        """    
    
        # base_regrid is a general function the regrid blobs into a grid that
        # does not need to have a node at the point [0.0, 0.0], for that it
        # uses the original grid placement of particles to redistribute them.
        # in here it is assumed that the redistribution grid always contains the
        # point [0.0, 0.0]. In order to re-use base_regrid, the xBounds and
        # yBounds variables are computes in order to redistribute particles
        # into a grid that contains the point [0.0, 0.0]. Essentially, the
        # redistribution is done into the following initial grid:
        #
        #       -------------------------
        #       |       |       |       |
        #       |   X   |   X   |   X   |
        #       |  (3)  |  (7)  |  (11) |
        #       -------------------------
        #       |       |       |      O|
        #       |   X   |   X   |   X   |
        #       |  (2)  |  (6)  |  (10) |
        #       -------------------------
        #       |       |       |       |
        #       |   X   |   X   |   X   |
        #       |  (1)  |  (5)  |  (9)  |
        #       -------------------------
        #
        # Where (6) is the point: [0.0, 0.0]. For this reason the self.h*1.5
        # terms appear.
        xBounds = _numpy.array([-self.h*1.5, self.h*1.5])
        yBounds = _numpy.array([-self.h*1.5, self.h*1.5])
        
        # Note that here c=0.0, because we only want to redistribute the blobs
        # no diffusion is to be performed.
        self.__x, self.__y, self.__g = _base_redistribute(self.__x,self.__y,self.__g,self.__sigma,self.__overlap,\
                                                                   xBounds,yBounds,interpKernel=interpKernel,c=0.0)
        
        
    def __read_input_wField(self,wField):
        """
            A function that reads the input wField given to self.__init__.
            
            Usage
            -----
            .. code-block:: python
            
                self.__read_input_wField(wField)
                
            Parameters
            ----------
            wField : tuple of three numpy.array(float64) (1) or tuple of a function
                     and two numpy.array(float64) (2)
                     (1) (`xBlobs`, `yBlobs`, `gBlob`) defines the initial distribution of
                         vortex blobs. `xBlobs` contains the x coordinates, `yBlobs`
                         contains the y coordinates and `gBlobs` contains the circulations.
                         Blobs are recursively redistributed into a regular mesh. This
                         mesh is assumed to have a particle located at (0,0). All
                         arrays must have shape (1,nBlobs).
                     (2) (`wExactFunction`, `xBounds`, `yBounds`) defines the initial
                         distribution of vorticity and the region of interest, in x
                         is the region between `xBounds` [0] and `xBounds` [1] and in y
                         is the region between `yBounds` [0] and `yBounds` [1]. It is
                         assumed that `xBounds` [1] > `xBounds` [0] and `yBounds` [1] > `yBounds` [0].
                         Vortex blobs are then generated inside this region with the
                         overlap `overlap` and spacing `h`.
            
            
            :First Added:   2013-12-19
            :Last Modified: 2013-12-19
            :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
            :License:       GNU GPL version 3 or any later version
            
        """
        
        """    
        Reviews:  
        
        """
        
        # check if wField contains the blobs or the vorticity function
        if type(wField[0]) == _numpy.ndarray: # wField has the blobs
            self.__x,self.__y,self.__g = self.__read_input_wField_blobs(wField)
            
        elif type(wField[0]) == __types.FunctionType:   # wField has the vorticity function
            self.__x,self.__y,self.__g = self.__read_input_wField_function(wField)
        
        else:   # wField is not a numpy.ndarray and not a function, then raise error
            raise TypeError('wField is not a tuple with x,y coordinates and g circulation of the blobs or a function with vorticity.')
        
        
    def __read_input_wField_blobs(self,wField):
        """
            A function that reads the input wField when given as three
            numpy.array with x and y coordinates of the blobs and circulation
            of the blobs. Returns x,y,g if correctly given in wField.
            
            Usage
            -----
            .. code-block:: python
            
                self.__read_input_wField_blobs(wField)
                
            Parameters
            ----------
            wField : tuple of three numpy.array(float64) 
                     (`xBlobs`, `yBlobs`, `gBlob`) defines the initial distribution of
                     vortex blobs. `xBlobs` contains the x coordinates, `yBlobs`
                     contains the y coordinates and `gBlobs` contains the circulations.
                     Blobs are recursively redistributed into a regular mesh. This
                     mesh is assumed to have a particle located at (0,0). All
                     arrays must have shape (1,nBlobs).
                     
            
            :First Added:   2013-12-19
            :Last Modified: 2013-12-19
            :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
            :License:       GNU GPL version 3 or any later version
            
        """
        
        """    
        Reviews:  
        
        """
        
        # check if there are three elements (xBlobs,yBlobs,gBlobs)
        if len(wField) != 3: # if there are not exaclty three arraysm something is wrong so raise error
            raise ValueError('wField must have exactly 3 numpy.array\'s. It has %d.' % len(wField))
        
        # check if all elements of wField are numpy.array and that sizes
        # are all the same
        for k,temp_array in enumerate(wField):
            # first check if it is a numpy.array
            if type(temp_array) != _numpy.ndarray: # raise error if it is not
                raise ValueError('wField[%d] must be numpy.array.' % k)
                
            # check the size and compare with previous sizes
            if k==0: # for the first array get the size
                array_size = temp_array.shape
            
            if array_size != temp_array.shape:
                raise ValueError('wField[%d] has a different size. All arrays must have the same size' % k)
                
        # flatten all arrays so that they are one dimensional
        # x coordinate, y coordinate and g circulation
        return wField[0].flatten().copy(),\
               wField[1].flatten().copy(),\
               wField[2].flatten().copy()   
        
    
    def __read_input_wField_function(self,wField):
        """
            A function that reads the input wField when given as a function
            and the bounds. Returns x,y,g if correctly given in wField.
            
            Usage
            -----
            .. code-block:: python
            
                self.__read_input_wField_function(wField)
                
            Parameters
            ----------
            wField : tuple of a function and two numpy.array(float64) (2)
                     (`wExactFunction`, `xBounds`, `yBounds`) defines the initial
                     distribution of vorticity and the region of interest, in x
                     is the region between `xBounds` [0] and `xBounds` [1] and in y
                     is the region between `yBounds` [0] and `yBounds` [1]. It is
                     assumed that `xBounds` [1] > `xBounds` [0] and `yBounds` [1] > `yBounds` [0].
                     Vortex blobs are then generated inside this region with the
                     overlap `overlap` and spacing `h`.
                     
            
            :First Added:   2013-12-19
            :Last Modified: 2013-12-19
            :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
            :License:       GNU GPL version 3 or any later version
            
        """
        
        """    
        Reviews:  
        
        """
        
        # check if there are three elements (xBlobs,yBlobs,gBlobs)
        if len(wField) != 3: # if there are not exaclty three arrays something is wrong so raise error
            raise ValueError('wField must have exactly 3 entries. It has %d.' % len(wField))
        
        # check if the first element is a function
        if type(wField[0]) != __types.FunctionType: # if is not a function raise error
            raise TypeError('wField[0] must be a function. It is %s.' % type(wField[0]))
        
        # check if first element is a function with two arguments
        if len(__inspect.getargspec(wField[0]).args) != 2: # if it is not 2 (x,y) raise error
            raise ValueError('wField[0] must be a function that takes two inputs (x,y). It takes %s.' % str(len(__inspect.getargspec(wField[0]).args)))
            
        # check if the second and third elements of wField are numpy.array and that 
        # it has size (2,)
        for k,temp_array in enumerate(wField):
            if k != 0:
                # first check if it is a numpy.array
                if type(temp_array) != _numpy.ndarray: # raise error if it is not
                    raise ValueError('wField[%d] must be numpy.array.' % k)
                    
                # check the size
                if temp_array.shape != (2,):
                    raise ValueError('wField[%d] must have size (2,). It has size %s.' % (k,str(wField[k].shape)))
        
        # compute the number of blobs in x and y directions inside the box defined by wField[1] (xBounds)
        # wField[2] (yBounds)
        nBlobs_x = _numpy.floor((wField[1][1]-wField[1][0])/self.__h).astype(long)
        nBlobs_y = _numpy.floor((wField[2][1]-wField[2][0])/self.__h).astype(long)
        
        # generate the coordinates of the vortex blobs
        x = wField[1][0]+0.5*self.__h + _numpy.arange(0,nBlobs_x)*self.__h      # nBlobs_x evenly spaced in x \in [xBoundsDomain[0], xBoundsDomain[1]]
        y = wField[2][0]+0.5*self.__h + _numpy.arange(0,nBlobs_y)*self.__h      # nBlobs_y evenly spaced in y \in [yBoundsDomain[0], yBoundsDomain[1]]
        x,y = _numpy.meshgrid(x,y)                                           # generate the 2d grid of blobs in [xBoundsDomain[0], xBoundsDomain[1]] x [yBoundsDomain[0], yBoundsDomain[1]]
        
        # flatten the grids into vectors
        x = x.flatten()                           
        y = y.flatten()
        
        # compute the circulations of the blobs
        g = wField[0](x,y)*self.__h*self.__h
        
        return x,y,g
        
        
    def __read_input_flow_parameters(self,vInf,overlap,h,deltaTc,nu):
        """
        A function that reads the input flow parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_input_flow_parameters(vInf,overlap,h,deltaTc,nu)
            
        Parameters
        ----------
        vInf : numpy.array(float)
               The free-stream velocity: vx = `vInf` [0] and vy = `vInf` [1].
                        
        overlap : float64
                  The overlap ratio between neighboring blobs. It is related to
                  the core size of the blob, `sigma`, and to the spacing 'h' by
                  the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.
                  
        h : float64
            The size of the cell associated to the vortex blobs. Corresponds to
            the minimum spacing between the core of two neighboring cells. It is related to
            the core size of the blob, `sigma`, and to the spacing 'h' by
            the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`
            
        deltaTc : float64
                  The size of the convective time step :math:`\Delta t_{c}`.
                  
        nu : float64
             The fluid kinematic viscosity, used to calculate the diffusion coefficient
             :math:`c^{2}` and diffusion time step `deltaTd`, :math:`\Delta t_{d}`.
        
        
        :First Added:   2013-12-19
        :Last Modified: 2013-12-19
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """
        
        # check if vInf (free stream velocity) is a numpy.array with two float values
        if type(vInf) != _numpy.ndarray:
            raise TypeError('vInf must be a numpy.ndarray. It is %s.' % str(type(vInf)))
        # check now if the type is float64
        elif vInf.dtype != 'float64':
            raise TypeError('vInf must have dtype float64. It has type %s.' % str(vInf.dtype))
        # check if it has only 2 values
        elif vInf.shape != (2,):
            raise ValueError('vInf must have shape (2,0). It has shape %s.' % str(vInf.shape))
        
        self.__vInf = vInf
        
        # check if overlap is positive
        if overlap <= 0.0:
            raise ValueError('overlap must be > 0.0. It is %f.' % overlap)
        
        self.__overlap = overlap

        # check if cell size is positive
        if h <= 0.0:
            raise ValueError('h must be > 0.0. It is %f.' % h)

        self.__h = h

        # compute the core spreading, sigma
        self.__sigma = h/overlap                
        
        # check if deltaTc is positive
        if deltaTc <= 0.0:
            raise ValueError('deltaTc must be > 0.0. It is %f.' % deltaTc)
        
        self.__deltaTc = deltaTc
        
        # check if nu is positive
        if nu <= 0.0:
            raise ValueError('nu must be > 0.0. It is %f.' % nu)
        
        self.__nu = nu
        
    
    def __read_input_computation_parameters(self,parameters):
        """
        A function that reads the input computation parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_input_computation_parameters(parameters)
            
        Parameters
        ----------
        parameters : dict, optional
                     A dictionary containing all additional parameters of the
                     simulation: `stepRedistribution`, `integrationMethod`,
                     `computationMethod`, `stepPopulationControl`, `gThreshold`.
                     
                     `stepRedistribution` (int) the redistribution step frequency.
                     
                     `integrationMethod` ('fe', 'rk4') the type of time integrator
                     used: 'fe' forward Euler, 'rk4' Runge-Kutta 4th order.
                     
                     'computationMethod' (tuple) refering with the type of Biot-Savart
                     solver ('direct', 'fmm') and the type of hardware to use
                     ('cpu','gpu').
                     
                     'stepPopulationControl' (int) the population control step frequency.
                     
                     'gThreshold' (tuple) with minimum value of circulation to
                     consider for each vortex blob and maximum value of total
                     circulation change by discarding particles.
                     
                     'diffusion' (dict) defining the diffusion method to use and
                     the parameters of the diffusion method. It has the key 'method'
                     that defines the method to use: 'regrid_wee' (based in [1]_).
                     Parameters for each different method:
                     'regrid_wee':
                         `c2` (float64), free parameter that specifies the stability
                                         of the diffusion scheme. Best values: [1/6,1/2].        
        
        :First Added:   2013-12-19
        :Last Modified: 2014-01-23
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        """
        
        """    
        Reviews:  1- Added diffusion parameters. (2014-01-23)
        
        """
        
        # stepRedistribution
        if type(parameters['stepRedistribution']) not in (int,long): # stepRedistribution must be an integer
            raise TypeError('parameters[\'stepRedistribution\'] must be an int. It is %s.' % str(type(parameters['stepRedistribution'])))
        elif parameters['stepRedistribution']<0: # stepRedistribution must be positive
            raise ValueError('parameters[\'stepRedistribution\'] must be positive. It is %d.' % parameters['stepRedistribution'])
            
        # integrationMethod
        if parameters['integrationMethod'] not in options.time_integrator_options:
            raise ValueError(('parameters[\'integrationMethod\'] must be one of ' +\
                             '\'%s\','*(len(options.time_integrator_options)-1) +\
                             '\'%s\'.' + ' It is %s.') % (options.time_integrator_options + (str(parameters['integrationMethod']),)))

        # computationMethod
        # check first the biot-savart option
        if parameters['computationMethod'][0] not in options.biot_savart_options:
            raise ValueError(('parameters[\'computationMethod\'][0] must be one of ' +\
                             '\'%s\','*(len(options.biot_savart_options)-1) +\
                             '\'%s\'.' + ' It is %s.') % (options.biot_savart_options + (str(parameters['computationMethod'][0]),)))
        # check the hardware option
        if parameters['computationMethod'][1] not in options.hardware_options:
            raise ValueError(('parameters[\'computationMethod\'][1] must be one of ' +\
                             '\'%s\','*(len(options.hardware_options)-1) +\
                             '\'%s\'.' + ' It is %s.') % (options.hardware_options + (str(parameters['computationMethod'][1]),)))
        
        # stepPopulationControl
        if type(parameters['stepPopulationControl']) not in (int,long): # stepPopulationControl must be an integer
            raise ValueError('parameters[\'stepPopulationControl\'] must be an int. It is %s.' % str(type(parameters['stepPopulationControl'])))
        elif parameters['stepPopulationControl']<0: # stepPopulationControl must be positive
            raise ValueError('parameters[\'stepPopulationControl\'] must be positive. It is %d.' % parameters['stepPopulationControl'])
        
        # gThreshold
        # check first for the minimum value of circulation to consider
        if parameters['gThreshold'][0] < 0: # gThreshold[0] must be positive or 0
            raise ValueError('parameters[\'gThreshold\'][0] must be positive. It is %f.' % parameters['gThreshold'][0])
        
        # check for the maximum value of total circulation allowed to change in population control
        if parameters['gThreshold'][1] < 0: # gThreshold[1] must be positive or 0
            raise ValueError('parameters[\'gThreshold\'][1] must be positive. It is %f.' % parameters['gThreshold'][1])

        # check if diffusion parameters are correctly given
        if type(parameters['diffusion']) == dict:
            self.__read_input_diffusion_parameters(parameters['diffusion'])
        else: # in case is not dict raise error
            raise TypeError('parameters[\'diffusion\'] must be of type dict. It is %s.' % type(parameters['diffusion']))
                
    
    def __read_input_diffusion_parameters(self,parameters):
        """
        A function that reads the input diffusion parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_input_diffusion_parameters(parameters)
            
        Parameters
        ----------
        parameters : dict, optional
                     A dictionary containing all additional parameters of the diffusion
                     method. The parameters depend on the diffusion method chosen.

                     `method` (string) (non-optional) defines the diffusion method used. Can be: 
                     
                     'regrid_wee' :: implements the diffusion-regrid method described in [1]_
                     
                     Parameters for each different method:
                     
                     'regrid_wee':
                         `c2` (float64), free parameter that specifies the stability
                                         of the diffusion scheme. Best values: [1/6,1/2].
        
        
        .. [1] Wee, D., Ahmed, F., Ghoniem, F., Modified interpolation
              kernels for treating diffusion and remeshing in vortex
              methods, Journal of Computational Physics, 213, 239--263, 2006.
        
        
        :First Added:   2013-12-22
        :Last Modified: 2013-12-22
        :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        """
        
        """    
        Reviews:  
        
        """
        
        # method
        if parameters['method'] == 'regrid_wee':
            # the parameters of the 'regrid_wee' method
            if 'c2' in parameters.keys(): # c2 parameter is required
                if type(parameters['c2']) == float:
                    if parameters['c2'] < 0.0: # c2 values must be positive
                        raise ValueError('parameters[\'c2\'] must be positive. It is %f.' % parameters['c2'])
                else:
                    raise TypeError('parameters[\'c2\'] must be a float. It is %s.' % type(parameters['c2']))
            else:
                raise ValueError('parameters[\'c2\'] must exist for diffusion method \'regrid_wee\'.')
            
            # define the diffusive time step
            # for stability reasons, c should be such that c^2 \in [1/6,1/2] 
            self.__deltaTd = (parameters['c2']*self.__h*self.__h)/self.__nu
            
        else:
            raise ValueError('parameters[\'method\'] must be one of \'regrid_wee\'. It is %s.' % str(parameters['method']))

        
        
    @property
    def stepRedistribution(self):
        """
        `stepRedistribution` (int) the redistribution step frequency.
            
        Usage
        -----
        .. code-block:: python
        
        self.stepRedistribution
        """
        return self.__parameters['stepRedistribution']

    @stepRedistribution.setter
    def stepRedistribution(self,stepRedistributionnew):
        raise AttributeError('stepRedistribution cannot be manually set by the user, user can only get the value of stepRedistribution.')
            
    @stepRedistribution.deleter
    def stepRedistribution(self):
        raise AttributeError('stepRedistribution cannot be manually deleted by the user, user can only get the value of stepRedistribution.')
    
    @property
    def integrationMethod(self):
        """
        `integrationMethod` ('fe','rk4') the type of time integrator
        used: 'fe' forward Euler, 'rk4' Runge-Kutta 4th order.
            
        Usage
        -----
        .. code-block:: python
    
            self.integrationMethod
        """
        return self.__parameters['integrationMethod']
    
    @integrationMethod.setter
    def integrationMethod(self,integrationMethodnew):
        raise AttributeError('integrationMethod cannot be manually set by the user, user can only get the value of integrationMethod.')
            
    @integrationMethod.deleter
    def integrationMethod(self):
        raise AttributeError('integrationMethod cannot be manually deleted by the user, user can only get the value of integrationMethod.')

    @property
    def computationMethod(self):
        """
        'computationMethod' (tuple) with the type of Biot-Savart
        solver ('direct', 'fmm') and the type of hardware to use
        ('cpu','gpu').
        
        Usage
        -----
        .. code-block:: python
        
            self.computationMethod
        """
        return self.__parameters['computationMethod']
    
    @computationMethod.setter
    def computationMethod(self,computationMethodnew):
        raise AttributeError('computationMethod cannot be manually set by the user, user can only get the value of computationMethod.')
            
    @computationMethod.deleter
    def computationMethod(self):
        raise AttributeError('computationMethod cannot be manually deleted by the user, user can only get the value of computationMethod.')

    @property
    def stepPopulationControl(self):
        """
            'computationMethod' (tuple) with the type of Biot-Savart
            solver ('direct', 'fmm') and the type of hardware to use
            ('cpu','gpu').
            
            Usage
            -----
            .. code-block:: python
        
                self.stepPopulationControl
        """
        return self.__parameters['stepPopulationControl']
    
    @stepPopulationControl.setter
    def stepPopulationControl(self,stepPopulationControlnew):
        raise AttributeError('stepPopulationControl cannot be manually set by the user, user can only get the value of stepPopulationControl.')
            
    @stepPopulationControl.deleter
    def stepPopulationControl(self):
        raise AttributeError('stepPopulationControl cannot be manually deleted by the user, user can only get the value of stepPopulationControl.')

    @property
    def gThreshold(self):
        """
            'gThreshold' (tuple) with minimum value of circulation to
            consider for each vortex blob and maximum value of total
            circulation change by discarding particles.
            
            Usage
            -----
            .. code-block:: python
        
                self.gThreshold
        """
        return self.__parameters['gThreshold']
    
    @gThreshold.setter
    def gThreshold(self,gThreshold):
        raise AttributeError('gThreshold cannot be manually set by the user, user can only get the value of gThreshold.')
            
    @gThreshold.deleter
    def gThreshold(self):
        raise AttributeError('gThreshold cannot be manually deleted by the user, user can only get the value of gThreshold.')


    # set the property h in order to access self.__h in a safe way
    @property
    def h(self):
        """
            The size of the cell associated to the vortex blobs. Corresponds to
            the minimum spacing between the core of two neighboring cells. It is related to
            the core size of the blob, `sigma`, and to the overlap 'overlap' by
            the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.
            
            Usage
            -----
            .. code-block:: python
        
                self.h

        """
        return self.__h
        
    @h.setter
    def h(self,hnew):
        raise AttributeError('h cannot be manually set by the user, user can only get the value of h.')
            
    @h.deleter
    def h(self):
        raise AttributeError('h cannot be manually deleted by the user, user can only get the value of h.')
        
    # set the property deltaTc in order to access self.__sigma in a safe way
    @property
    def sigma(self):
        """
            The core size of the vortex blobs.
            
            Usage
            -----
            .. code-block:: python
        
                self.sigma

        """
        return self.__sigma
        
    @sigma.setter
    def sigma(self,sigmanew):
        raise AttributeError('sigma cannot be manually set by the user, user can only get the value of sigma.')
            
    @sigma.deleter
    def sigma(self):
        raise AttributeError('sigma cannot be manually deleted by the user, user can only get the value of sigma.')
        
        
    # set the property deltaTc in order to access self.__overlap in a safe way
    @property
    def overlap(self):
        """
            The overlap ratio between neighboring blobs. It is related to
            the core size of the blob, `sigma`, and to the spacing 'h' by
            the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.
            
            Usage
            -----
            .. code-block:: python
        
                self.overlap

        """
        return self.__overlap
        
    @overlap.setter
    def overlap(self,overlapnew):
        raise AttributeError('overlap cannot be manually set by the user, user can only get the value of overlap.')
            
    @overlap.deleter
    def overlap(self):
        raise AttributeError('overlap cannot be manually deleted by the user, user can only get the value of overlap.')    
        
        
    # set the property deltaTc in order to access self.__deltaTc in a safe way
    @property
    def deltaTc(self):
        """
            The size of the convective time step :math:`\Delta t_{c}`.
            
            Usage
            -----
            .. code-block:: python
        
                self.deltaTc

        """
        return self.__deltaTc
        
    @deltaTc.setter
    def deltaTc(self,deltaTcnew):
        raise AttributeError('deltaTc cannot be manually set by the user, user can only get the value of deltaTc.')
            
    @deltaTc.deleter
    def deltaTc(self):
        raise AttributeError('deltaTc cannot be manually deleted by the user, user can only get the value of deltaTc.')
        
    # set the property deltaTd in order to access self.__deltaTd in a safe way
    @property
    def deltaTd(self):
        """
            The size of the diffusive time step :math:`\Delta t_{d}`.
            
            Usage
            -----
            .. code-block:: python
        
                self.deltaTd

        """
        return self.__deltaTd
        
    @deltaTd.setter
    def deltaTd(self,deltaTdnew):
        raise AttributeError('deltaTd cannot be manually set by the user, user can only get the value of deltaTd.')
            
    @deltaTd.deleter
    def deltaTd(self):
        raise AttributeError('deltaTd cannot be manually deleted by the user, user can only get the value of deltaTd.')

    
    # set the property vInf in order to access self.__vInf in a safe way
    @property
    def vInf(self):
        """
            The free stream velocity.
            
            Usage
            -----
            .. code-block:: python
        
                self.vInf

        """
        return self.__vInf
        
    @vInf.setter
    def vInf(self,vInfnew):
        raise AttributeError('vInf cannot be manually set by the user, user can only get the value of vInf.')
            
    @vInf.deleter
    def vInf(self):
        raise AttributeError('vInf cannot be manually deleted by the user, user can only get the value of vInf.')
        
        
    # set the property vInf in order to access self.__vInf in a safe way
    @property
    def nu(self):
        """
            The fluid kinematic viscosity, used to calculate the diffusion coefficient
             :math:`c^{2}` and diffusion time step `deltaTd`, :math:`\Delta t_{d}`.
            
            Usage
            -----
            .. code-block:: python
        
                self.nu

        """
        return self.__nu
        
    @nu.setter
    def nu(self,nunew):
        raise AttributeError('nu cannot be manually set by the user, user can only get the value of nu.')
            
    @nu.deleter
    def nu(self):
        raise AttributeError('nu cannot be manually deleted by the user, user can only get the value of nu.')
        
    # set the property x in order to access self.__x in a safe way
    @property
    def x(self):
        """
            The x coordinates of the vortex blobs.
            
            Usage
            -----
            .. code-block:: python
        
                self.x

        """
        return self.__x
        
    @x.setter
    def x(self,xnew):
        raise AttributeError('x cannot be manually set by the user, user can only get the value of x.')
            
    @x.deleter
    def x(self):
        raise AttributeError('x cannot be manually deleted by the user, user can only get the value of x.')
        
    # set the property y in order to access self.__y in a safe way
    @property
    def y(self):
        """
            The y coordinates of the vortex blobs.
            
            Usage
            -----
            .. code-block:: python
        
                self.y

        """
        return self.__y
        
    @y.setter
    def y(self,ynew):
        raise AttributeError('y cannot be manually set by the user, user can only get the value of y.')
            
    @y.deleter
    def y(self):
        raise AttributeError('y cannot be manually deleted by the user, user can only get the value of y.')
        
        
    # set the property y in order to access self.__y in a safe way
    @property
    def g(self):
        """
            The g coordinates of the vortex blobs.
            
            Usage
            -----
            .. code-block:: python
        
                self.g

        """
        return self.__g
        
    @g.setter
    def g(self,gnew):
        raise AttributeError('g cannot be manually set by the user, user can only get the value of g.')
            
    @g.deleter
    def g(self):
        raise AttributeError('g cannot be manually deleted by the user, user can only get the value of g.')
        
    # set the property num_blobs in order to get the number of blobs
    @property
    def num_blobs(self):
        """
            The number of blobs.
            
            Usage
            -----
            .. code-block:: python
        
                self.num_blobs

        """
        return self.__x.shape[0]
        
    @num_blobs.setter
    def num_blobs(self,gnew):
        raise AttributeError('num_blobs cannot be set by the user, user can only get the value of num_blobs.')
            
    @num_blobs.deleter
    def num_blobs(self):
        raise AttributeError('num_blobs cannot be deleted by the user, user can only get the value of num_blobs.')
        
    # set the property diffusionParameters in order to get the diffusion parameters
    # of the simulation
    @property      
    def diffusionParameters(self):
        """
            The diffusion parameters. It is a dictionary containing all the
            parameters of the diffusion method used for the simulation.
            
            Usage
            -----
            .. code-block:: python
                
                self.diffusionParameters
                
        """
        
        return self.__diffusionParameters
        
    @diffusionParameters.setter
    def diffusionParameters(self,diffusionParametersNew):
        raise AttributeError('diffusionParameters cannot be set by the user, user can only get the value of diffusionParameters.')
            
    @diffusionParameters.deleter
    def diffusionParameters(self):
        raise AttributeError('diffusionParameters cannot be deleted by the user, user can only get the value of diffusionParameters.')
    
























