__doc__ = """

Vortex particle (blobs) solver master class. 

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
:Last updated:  2014-03-05
:Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
:License:       GNU GPL version 3 or any later version
      
"""

"""
Reviews:
-------
        1- Added method VortexBlobs.regrid                 (apalha, 2014-01-15)
        2- Added method VortexBlobs.populationControl      (apalha, 2014-01-20)
        3- Added method VortexBlobs.evaluateVelocity       (apalha, 2014-01-21)
        4- Added method VortexBlobs.evolve, restructured
           input parameters, correct
           VortexBlobs.evaluateVelocity, added properties
           for most of the attributes for easy and secure
           access to simulation parameters. Comments of the
           input parameters have been greatly revised and
           corrected.                                      (apalha, 2014-02-04)
           
        5- Class renamed/restructured (VortexBlobs to Blobs) (lmanickathan, 2014-03-05)

"""


__all__ = ['Blobs']

# External packages
import numpy as _numpy
import types as __types
import inspect as __inspect

# pHyFlow packages
from pHyFlow.aux import ErrorOutput as _ErrorOutput # import error output functions

# import blob options
from pHyFlow.blobs import blobOptions # import options definitions for pHyFlow

# import required base functions from where to build up the class
from pHyFlow.blobs.base.induced import velocity as _base_velocity
#from pHyFlow.blobs.base.induced import vorticity_blobs as _base_vorticity_blobs
#from pHyFlow.blobs.base.induced import vorticity as _base_vorticity
from pHyFlow.blobs.base.regrid import Regrid as _base_redistribute
from pHyFlow.blobs.base.regrid import PopulationControl as _base_populationControl
from pHyFlow.blobs.base.evolution import convection as _base_convection
from pHyFlow.blobs.base.evolution import diffusion_wee as _base_diffusion_wee
from pHyFlow.blobs.base.evolution import diffusion_wee as _base_diffusion_wee
from pHyFlow.blobs.base.evolution import diffusion_tutty as _base_diffusion_tutty


# default parameters

default_diffusion_params = {'method':'regrid_wee','c2':'optimized'}
default_time_integration_params = {'method':'rk4'}
default_blob_control_params = {'stepRedistribution':1,'stepPopulationControl':1,\
                               'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}
default_velocity_computation_params = {'method':'fmm','hardware':'gpu'}


class Blobs(object):
    r"""
    Vortex blob solver for the Navier-Stokes equation. No solid boundaries
    are implemented.
        
    Usage
    -----
    .. code-block:: python
    
        Blobs(wField,vInf,nu,deltaTc,h,overlap,
              timeIntegrationParams={'method':'rk4'},
              blobControlParams={'stepRedistribution':1,'stepPopulationControl':1,
              'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8},
              velocityComputationParams={'method':'fmm','hardware':'gpu'},
              diffusionParams={'method':'regrid_wee','c2':'optimized'})
        
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

    nu : float64
         The fluid kinematic viscosity, used to calculate the diffusion coefficient
         :math:`c^{2}` and diffusion time step `deltaTd`, :math:`\Delta t_{d}`.

    deltaTc : float64
              The size of the convective time step :math:`\Delta t_{c}`.

    h : float64
        The size of the cell associated to the vortex blobs. Corresponds to
        the minimum spacing between the core of two neighboring cells. It is related to
        the core size of the blob, `sigma`, and to the spacing 'h' by
        the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`

    overlap : float64
              The overlap ratio between neighboring blobs. It is related to
              the core size of the blob, `sigma`, and to the spacing 'h' by
              the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`.

    timeIntegrationParams : dict, optional
                            A dictionary containing all time integration parameters of the simulation. Contains
                            the definition of the time integration scheme possibly additional parameters
                            specific to the scheme.

                            Contains always the following keys:

                            'method' : string
                                       Specifies the time integration scheme.
                                       'rk4' : Runge-Kutta 4th order.
                                       'euler' : Forward Euler scheme (avoid using this).

                            Other keys are optional and specific to the time integration scheme.

                            default:
                                {'method':'rk4'}

    blobControlParams : dict, optional
                        A dictionary containing all parameters related to blob control, from population
                        control to regularity (redistribution).

                        Contains the following keys:

                        'stepRedistribution' : int
                                               The redistribution step frequency.

                        'stepPopulationControl' : int
                                                  The population control step frequency. The frequency at which
                                                  blobs with too small circulation are removed, in order to
                                                  keep the number of blobs from growing too fast.

                        'gThresholdLocal' : float64
                                            Minimum value of circulation to consider for each vortex blob when
                                            selecting blobs to remove during population control.

                        'gThresholdGlobal' : float64
                                             Maximum value of variation of total vorticity due to the removal of
                                             blobs during population control.

                        default:
                            {'stepRedistribution':1,'stepPopulationControl':1, 'gThresholdLocal':1e-8,'gThresholdGlobal':1e-8}

    velocityComputationParams : dict, optional
                                A dictionary containing all the parameters related to the computation of induced
                                velocities. Specifies computation scheme (direct or fmm) and hardware to use (cpu
                                or gpu).

                                Contains the following keys:

                                'method' : string
                                           Specifies the induced velocity computation scheme.
                                           'direct' : direct calculation O(n^2) time, but exact.
                                           'fmm' : fmm (fast multipole method) calculation O(n log n) time,
                                                   but approximated.

                                'hardware' : string
                                             Specifies the hardware to use in the computation of velocity.
                                             'cpu' : use the cpu.
                                             'gpu' : use the gpu (CUDA enabled).

                                default:
                                    {'method':'fmm','hardware':'gpu'}
              
    diffusionParams : dict, optional
                      A dictionary containing all the parameters related to the computation of the diffusion step.
                      Specifies the diffusion scheme and other specific parameters.

                      Contains the following keys:

                      'method' : string
                                 Specifies the diffusion scheme used.
                                 'regrid_wee' : diffusion scheme presented by Wee et al. in [1]_, fusing VRM method
                                                with regriding in one step.

                      Specific keys:

                      For 'regrid_wee':

                            'c2' : float64 or string
                                   the c^2 parameter introduced by Wee in his scheme, see [1]_ . Should be
                                   such that c2 is in the interval [1/6, 1/2].

                                   As an alternative, 'optimized' can be given as input, which then chooses 1/3.

                      default:
                        {'method':'regrid_wee','c2':'optimized'}


    Attributes
    ----------
    blobControlParams
    computationMethod
    deltaTc
    deltaTd
    diffusionParams
    g
    gThresholdGlobal
    gThresholdLocal
    h
    integrationMethod        
    nu
    numBlobs
    overlap
    plotVelocity
    sigma
    stepDiffusion
    stepPopulationControl
    stepRedistribution
    timeIntegrationParams
    t
    tStep
    velocityComputationParams
    vInf        
    x
    y
    __blobControlParams : dict, optional
                        A dictionary containing all parameters related to blob control, from population
                        control to regularity (redistribution).

                        Contains the following keys:

                        'stepRedistribution' : int
                                               The redistribution step frequency.

                        'stepPopulationControl' : int
                                                  The population control step frequency. The frequency at which
                                                  blobs with too small circulation are removed, in order to
                                                  keep the number of blobs from growing too fast.

                        'gThresholdLocal' : float64
                                            Minimum value of circulation to consider for each vortex blob when
                                            selecting blobs to remove during population control.

                        'gThresholdGlobal' : float64
                                             Maximum value of variation of total vorticity due to the removal of
                                             blobs during population control.
    __diffusionParams : dict, optional
                      A dictionary containing all the parameters related to the computation of the diffusion step.
                      Specifies the diffusion scheme and other specific parameters.

                      Contains the following keys:

                      'method' : string
                                 Specifies the diffusion scheme used.
                                 'regrid_wee' : diffusion scheme presented by Wee et al. in [1]_, fusing VRM method
                                                with regriding in one step.

                      Specific keys:

                      For 'regrid_wee':

                            'c2' : float64 or string
                                   the c^2 parameter introduced by Wee in his scheme, see [1]_ . Should be
                                   such that c2 is in the interval [1/6, 1/2].

                                   As an alternative, 'optimized' can be given as input, which then chooses 1/3.
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
    __plotVelocity : bool
                     A flag that defines if velocity is to be plotted or not.
    __sigma : float64
              The core size of the vortex blobs.
    __stepDiffusion : int
                      The frequency of diffusion steps.
    __stepPopulationControl : int
                              The frequency of population control.
    __stepRedistribution : int
                           The frequency of redistribution of blobs.
    __t : float64
          The current time of the simulation.
    __tStep : int
              The current time step of the simulation.
    __timeIntegrationParams : dict, optional
                            A dictionary containing all time integration parameters of the simulation. Contains
                            the definition of the time integration scheme possibly additional parameters
                            specific to the scheme.

                            Contains always the following keys:

                            'method' : string
                                       Specifies the time integration scheme.
                                       'rk4' : Runge-Kutta 4th order.
                                       'euler' : Forward Euler scheme (avoid using this).

                            Other keys are optional and specific to the time integration scheme.
    __velocityComputationParams : dict, optional
                                A dictionary containing all the parameters related to the computation of induced
                                velocities. Specifies computation scheme (direct or fmm) and hardware to use (cpu
                                or gpu).

                                Contains the following keys:

                                'method' : string
                                           Specifies the induced velocity computation scheme.
                                           'direct' : direct calculation O(n^2) time, but exact.
                                           'fmm' : fmm (fast multipole method) calculation O(n log n) time,
                                                   but approximated.

                                'hardware' : string
                                             Specifies the hardware to use in the computation of velocity.
                                             'cpu' : use the cpu.
                                             'gpu' : use the gpu (CUDA enabled).
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
    Reviews:  (apalha, 2014-01-15) Added method self.regrid
              (apalha, 2014-01-20) Added method self.populationControl
              (apalha, 2014-01-21) Added method self.evaluateVelocity
              (apalha, 2014-02-10) Added method self.evolve, introduced possibility to change hardware
                                   used for computation, added method self.addBlobs,added method self.removeBlobs
              (apalha, 2014-04-14) Added Tutty diffusion scheme.

        
    """ 
        
    def __init__(self,wField,vInf,nu,deltaTc,h,overlap,\
                      timeIntegrationParams=default_time_integration_params,\
                      blobControlParams=default_blob_control_params,\
                      velocityComputationParams=default_velocity_computation_params,\
                      diffusionParams=default_diffusion_params):
        """    
        Reviews:  
        
        """ 

        #--------------------------
        # Check if input parameters all satisfy the required ranges and types
        self.__check_wField(wField)
        self.__check_vInf(vInf)
        self.__check_nu(nu)
        self.__check_deltaTc(deltaTc)
        self.__check_h(h)
        self.__check_overlap(overlap)
        self.__check_timeIntegrationParams(timeIntegrationParams)
        self.__check_blobControlParams(blobControlParams)
        self.__check_velocityComputationParams(velocityComputationParams)
        self.__check_diffusionParams(diffusionParams)
        #--------------------------

        #--------------------------
        # Read the input parameters

        self.__read_wField(wField)
        self.__read_vInf(vInf)
        self.__read_nu(nu)
        self.__read_deltaTc(deltaTc)
        self.__read_h(h)
        self.__read_overlap(overlap)
        self.__read_timeIntegrationParams(timeIntegrationParams)
        self.__read_blobControlParams(blobControlParams)
        self.__read_velocityComputationParams(velocityComputationParams)
        self.__read_diffusionParams(diffusionParams)
        #--------------------------


        # compute deltaTd, sigma, stepDiffusion, etc..

        #--------------------------
        # compute derived flow and blob parameters

        self.__sigma = self.__h/self.__overlap # the core spreading constant

        # define the diffusive time step
        # depending on the diffusion scheme chosen deltaTd will be computed
        self.__deltaTd, self.__stepDiffusion = self.__compute_deltaTd()

        #--------------------------

        #--------------------------
        # define and initialize time step counters and absolute time

        # define the time step counter
        self.__tStep = 0

        # define the current time instant
        self.__t = 0.0

        #--------------------------

        #--------------------------
        # define flags

        self.plotVelocity = False

        #--------------------------

        #--------------------------
        # perform redistribution and population control

        #if _numpy.abs(self.__g).sum() > _numpy.spacing(1):# TODO:Check
        if self.stepPopulationControl > 0: # only if population control is to be performed during time stepping
            self.populationControl()         

        if self.stepRedistribution > 0: # only if redistribution is to be performed during time stepping
            self.redistribute()      

        #--------------------------


    def addBlobs(self,xNew,yNew,gNew):
        r"""
        Adds the blobs defined by coordinates xNew,yNew and circulation gNew,
        to the existing blobs.

        Usage
        -----
        .. code-block :: python

            self.addBlobs(xNew,yNew,gNew)


        Parameters
        ----------
        xNew : numpy.array(float64), (1,nBlobsNew)
               the x coordinates of the new blobs.
        yNew : numpy.array(float64), (1,nBlobsNew)
               the y coordinates of the new blobs.
        gNew : numpy.array(float64), (1,nBlobsNew)
               the circulation g of the new blobs.

        Returns
        -------
        None returned.

        Attributes
        ----------
        x
        y
        g

        :First Added:   2014-02-10
        :Last Modified: 2014-02-10
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        # first check if the inputs are valid
        self.__check_input_wField_blobs((xNew,yNew,gNew))

        # add the new blobs to the existing blobs
        self.__x = _numpy.concatenate((self.__x,xNew))
        self.__y = _numpy.concatenate((self.__y,yNew))
        self.__g = _numpy.concatenate((self.__g,gNew))


    def evaluateVelocity(self,xEval,yEval,Ndirect=35,tol=1.0e-6,cutoff=None):
        r"""
        Computes the induced velocities generated by the vortex blobs. It can use
        a FMM approximation (faster, but less accurate) or direct calculation
        (slower, but machine precision accurate).
        
        Usage
        -----
        .. code-block :: python

            self.evaluateVelocity(xEval,yEval,Ndirect=35,tol=1.0e-6,cutoff=None)
    
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

        # convert the hardware flag into an int to use in _base_velocity
        if self.__velocityComputationParams['hardware'] == 'gpu': hardware = blobOptions.GPU_HARDWARE
        else: hardware = blobOptions.CPU_HARDWARE

        # convert the method flag into an int to use in _base_velocity
        if self.__velocityComputationParams['method'] == 'fmm': method = blobOptions.FMM_METHOD
        else: method = blobOptions.DIRECT_METHOD

        # compute the velocity at the evaluation points
        # computes de velocity with fmm or direct depending on tol
        vx,vy = _base_velocity(self.__x,self.__y,self.__g,self.__sigma,\
                               xEval=xEval,yEval=yEval,hardware=hardware,method=method,\
                               Ndirect=Ndirect,tol=tol,cutoff=cutoff)

        return vx,vy
        
        
    def evaluateVorticity(self,xEval,yEval):
        r"""
        Computes the induced velocities generated by the vortex blobs. It can use
        a FMM approximation (faster, but less accurate) or direct calculation
        (slower, but machine precision accurate).
        
        Usage
        -----
        .. code-block :: python

            self.evaluateVorticity(xEval,yEval,Ndirect=35,tol=1.0e-6,cutoff=None)
    
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
    
    
    def evolve(self):
        r"""
        Compute the new blobs as given by the vorticity formulation of
        Navier-Stokes equation. The method used is a splitting method which is
        constituted by a convection step followed by a diffusion step. Depending
        on the diffusion scheme selected, the diffusion step can be together
        with a convection step or at a prescribed frequency, depending on the
        relation between deltaTc (convection time step) and deltaTd (diffusion
        time step).

        Usage
        -----
        .. code-block :: python

            self.evolve()

        Parameters
        ----------
        None
        
        Attributes
        ----------
        x
        y
        g


        :First Added:   2014-02-10
        :Last Modified: 2014-02-10
        :Copyright:     Copyright (C) 2013,2014 Artur Palha, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """

        """
            Reviews: (apalha, 2014-04-14) Added diffusion as a callable function.
                     (apalha, 2014-04-15) Put redistribution before diffusion.
        """

        # convert the hardware flag into an int to use in _base_convection
        if self.__velocityComputationParams['hardware'] == 'gpu': hardware = blobOptions.GPU_HARDWARE
        else: hardware = blobOptions.CPU_HARDWARE

        # convert the method flag into an int to use in _base_convection
        if self.__velocityComputationParams['method'] == 'fmm': method = blobOptions.FMM_METHOD
        else: method = blobOptions.DIRECT_METHOD

        # convert the time integration method into an int to use in _base_convection
        if self.__timeIntegrationParams['method'] == 'rk4': integrator = blobOptions.RK4_INTEGRATOR
        elif self.__timeIntegrationParams['method'] == 'euler': integrator = blobOptions.FE_INTEGRATOR

        # convection step
        xBlobTemp,yBlobTemp,gBlobTemp = _base_convection(self.__deltaTc,self.__x,self.__y,self.__g,self.__sigma,
                                                         hardware=hardware,method=method,integrator=integrator,vInf=self.__vInf)

        # update the blobs
        self.__x = xBlobTemp
        self.__y = yBlobTemp
        self.__g = gBlobTemp

        # # diffusion step
        # if self.__stepDiffusion != 0: # corresponds to nu = 0.0, therefore no diffusion is to be computed
        #     if (self.__tStep % self.__stepDiffusion) == 0: # if the time step is a multiple of the stepDiffusion perform diffusion
        #
        #         # select the diffusion method to use
        #         if self.__diffusionParams['method'] == 'regrid_wee':
        #             # _base_diffusion_wee is a general function that implements the wee method for diffusion which is based
        #             # on regrid of the blobs into a grid that
        #             # does not need to have a node at the point [0.0, 0.0], for that it
        #             # uses the original grid placement of particles to redistribute them.
        #             # in here it is assumed that the redistribution grid always contains the
        #             # point [0.0, 0.0]. In order to re-use base_redistribute, the xBounds and
        #             # yBounds variables are computes in order to redistribute particles
        #             # into a grid that contains the point [0.0, 0.0]. Essentially, the
        #             # redistribution is done into the following initial grid:
        #             #
        #             #       -------------------------
        #             #       |       |       |       |
        #             #       |   X   |   X   |   X   |
        #             #       |  (3)  |  (7)  |  (11) |
        #             #       -------------------------
        #             #       |       |       |      O|
        #             #       |   X   |   X   |   X   |
        #             #       |  (2)  |  (6)  |  (10) |
        #             #       -------------------------
        #             #       |       |       |       |
        #             #       |   X   |   X   |   X   |
        #             #       |  (1)  |  (5)  |  (9)  |
        #             #       -------------------------
        #             #
        #             # Where (6) is the point: [0.0, 0.0]. For this reason the self.__h*1.5
        #             # terms appear.
        #             xBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])
        #             yBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])
        #
        #             # perform diffusion and store the results directly in the object
        #             self.__x,self.__y,self.__g = _base_diffusion_wee(self.__deltaTd,self.__nu,self.__h,
        #                                                              xBlobTemp,yBlobTemp,gBlobTemp,
        #                                                              self.__sigma,xBounds,yBounds,
        #                                                              overlap=self.__overlap)
        #
        #    else: # if the time step is not a multiple of stepDiffusion then don't perform diffusion and just update blobs
        #        self.__x = xBlobTemp
        #        self.__y = yBlobTemp
        #        self.__g = gBlobTemp
        #
        #else: # if the stepDiffusion is 0 then just update the blobs
        #    self.__x = xBlobTemp
        #    self.__y = yBlobTemp
        #    self.__g = gBlobTemp

        # redistribution step
        if self.__stepRedistribution != 0: # meaning that redistribution should be done
            if (self.diffusionParams['method'] == 'regrid_tutty'): # if the time step is a multiple of the stepRedistribution perform redistribution
                self.redistribute()

        self._diffusion()

        # redistribution step
        if self.__stepRedistribution != 0: # meaning that redistribution should be done
            if ((self.__tStep % self.__stepRedistribution) == 0) and (self.diffusionParams['method'] != 'regrid_tutty'): # if the time step is a multiple of the stepRedistribution perform redistribution
                self.redistribute()

        # update the time counter
        self._advanceTime()



        # population control step
        if self.__stepPopulationControl != 0: # meaning that population control should be done
            if (self.__tStep % self.__stepPopulationControl) == 0: # if the time step is a multiple of the stepPopulationControl perform population control
                self.populationControl()

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
        .. code-block :: python

            self.populationControl()
    
        Parameters
        ----------
        None
        
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
        self.__x,self.__y,self.__g = _base_populationControl(self.x,self.y,self.g,self.gThresholdLocal,self.gThresholdGlobal)
        
        
    def redistribute(self,interpKernel=0):
        r"""
        Redistribute vortex blobs over evenly spaced grid points with separation
        h = sigma * overlap.
        The interpolation kernel can be one of the following:
           - M4' (0)
        
        Usage
        -----
        .. code-block :: python

            self.redistribute(interpKernel=0)
            
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
    
        # base_redistribute is a general function the regrid blobs into a grid that
        # does not need to have a node at the point [0.0, 0.0], for that it
        # uses the original grid placement of particles to redistribute them.
        # in here it is assumed that the redistribution grid always contains the
        # point [0.0, 0.0]. In order to re-use base_redistribute, the xBounds and
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
        xBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])
        yBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])
        
        # Note that here c=0.0, because we only want to redistribute the blobs
        # no diffusion is to be performed.
        self.__x, self.__y, self.__g = _base_redistribute(self.__x,self.__y,self.__g,self.__sigma,self.__overlap,\
                                                                   xBounds,yBounds,interpKernel=interpKernel,c=0.0)


    def removeBlobs(self,iBlobs):
        r"""
        Remove the blobs indexed by the indices iBlobs or if iBlobs if of type bool, remove the blobs
        that have a value False in iBlobs (in this case iBlobs must have as many elements as self.numBlobs).

        Usage
        -----
        .. code-block :: python

            self.removeBlobs(iBlobs)

        Parameters
        ----------
        iBlobs :: numpy.array(int64) or numpy.array(bool) (in this case the dimensions have to be self.numBlobs)
                  Contains the indices of the blobs to remove.


        Attributes
        ----------
        x
        y
        g


        :First Added:   2014-02-10
        :Last Modified: 2014-02-11
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:  (apalha, 2014-02-11) Added the possibility of iBlobs to be of type bool. In this case
                                       iBlobs has to have as many elements as self.numBlobs. True means
                                       that the blob is kept and False that it is deleted.

        """

        # first check if iBlobs is a valid value
        if type(iBlobs) != _numpy.ndarray: # check if it is numpy.ndarray
            raise TypeError('iBlobs must be a numpy.array, it is a %s.' % str(type(iBlobs)))

        if iBlobs.dtype == _numpy.dtype('int64'): # check if it is int64
            if iBlobs.max() >= self.numBlobs:
                raise TypeError('iBlobs values must be smaller than the numBlobs: %d.' % self.numBlobs)

            # remove the blobs by removing the blobs whose index is in iBlobs
            self.__x = _numpy.delete(self.__x,iBlobs)
            self.__y = _numpy.delete(self.__y,iBlobs)
            self.__g = _numpy.delete(self.__g,iBlobs)

        elif iBlobs.dtype == _numpy.dtype('bool'): # check if it is bool
            if iBlobs.shape != (self.numBlobs,): # check if it has as many elements as the number of blobs
                raise TypeError('iBlobs of type bool must have as many elements as numBlobs: %d.' % self.numBlobs)

            # remove the blobs by removing the blobs whose value of iBlobs is False
            self.__x = self.__x[iBlobs]
            self.__y = self.__y[iBlobs]
            self.__g = self.__g[iBlobs]

        else:
            raise TypeError('iBlobs must be of type int64 or bool, it is a %s.' % str(iBlobs.dtype))

    def _advanceTime(self):
        r"""
        _advanceTime updates all time related variables: self.__tStep and self.__t.

        Usage
        -----
        .. code-block :: python

            self.__advanceTime()

        Parameters
        ----------

        Attributes
        ----------
        __tStep
        __t

        :First Added:   2014-02-10
        :Last Modified: 2014-03-06
        :Copyright:     Copyright (C) 2014 Artur Palha, Lento Manickathan **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # advance tStep
        self.__tStep += 1

        # advance t
        #self.__t += self.__deltaTc
        self.__t = self.__deltaTc * self.__tStep


    
    def _diffusion(self):
        """
        Function to perform the diffusion of the vortex blobs
        
        Usage
        -----
        .. code-block:: python
        
            self._diffusion()
            
        Parameters
        ----------
        None
        
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        x
        y
        g
        
        :First Added:   2014-02-10
        :Last Modified: 2014-03-06
        :Copyright:     Copyright (C) 2014 Artur Palha, Lento Manickathan **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """        
        # select the diffusion method to use
        if self.__diffusionParams['method'] == 'regrid_wee':
            # _base_diffusion_wee is a general function that implements the wee method for diffusion which is based
            # on regrid of the blobs into a grid that
            # does not need to have a node at the point [0.0, 0.0], for that it
            # uses the original grid placement of particles to redistribute them.
            # in here it is assumed that the redistribution grid always contains the
            # point [0.0, 0.0]. In order to re-use base_redistribute, the xBounds and
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
            # Where (6) is the point: [0.0, 0.0]. For this reason the self.__h*1.5
            # terms appear.
            xBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])
            yBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])

            # perform diffusion and store the results directly in the object
            self.__x,self.__y,self.__g = _base_diffusion_wee(self.__deltaTd,self.__nu,self.__h,
                                                             self.__x,self.__y,self.__g,
                                                             self.__sigma,xBounds,yBounds,
                                                             overlap=self.__overlap)

        if self.__diffusionParams['method'] == 'regrid_tutty':
            # _base_diffusion_wee is a general function that implements the wee method for diffusion which is based
            # on regrid of the blobs into a grid that
            # does not need to have a node at the point [0.0, 0.0], for that it
            # uses the original grid placement of particles to redistribute them.
            # in here it is assumed that the redistribution grid always contains the
            # point [0.0, 0.0]. In order to re-use base_redistribute, the xBounds and
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
            # Where (6) is the point: [0.0, 0.0]. For this reason the self.__h*1.5
            # terms appear.
            xBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])
            yBounds = _numpy.array([-self.__h*1.5, self.__h*1.5])

            # perform diffusion and store the results directly in the object
            self.__x,self.__y,self.__g = _base_diffusion_tutty(self.__deltaTd,self.__nu,self.__h,
                                                             self.__x,self.__y,self.__g,
                                                             self.__sigma,xBounds,yBounds,
                                                             overlap=self.__overlap)


    def __check_wField(self,wField):
        """
        A function that reads the input wField given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_wField(wField)
            
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
            self.__check_input_wField_blobs(wField)
            
        elif type(wField[0]) == __types.FunctionType:   # wField has the vorticity function
            self.__check_input_wField_function(wField)
        
        else:   # wField is not a numpy.ndarray and not a function, then raise error
            raise TypeError('wField is not a tuple with x,y coordinates and g circulation of the blobs or a function with vorticity.')
        
    
    def __read_wField(self,wField):
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
        
        
    def __check_input_wField_blobs(self,wField):
        """
        A function that checks the input wField when given as three
        numpy.array with x and y coordinates of the blobs and circulation
        of the blobs.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_input_wField_blobs(wField)
            
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
        

    def __read_input_wField_blobs(self,wField):
        """
        A function that reads the input wField when given as three
        numpy.array with x and y coordinates of the blobs and circulation
        of the blobs.

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

        # return a copy of the input blobs
        # x, y and g (circulation)
        return wField[0].copy(), wField[1].copy(), wField[2].copy()


    def __check_input_wField_function(self,wField):
        """
        A function that checks the input wField when given as a function
        and the bounds.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_input_wField_function(wField)
            
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
        
        # check if there are three elements (wFunction,xBounds,yBounds)
        if len(wField) != 3: # if there are not exaclty three values something is wrong so raise error
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
        
        
    def __check_vInf(self,vInf):
        """
        A function that checks the input free stream velocity given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_vInf(vInf)
            
        Parameters
        ----------
        vInf : numpy.array(float)
               The free-stream velocity: vx = `vInf` [0] and vy = `vInf` [1].
                        
        
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
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
            
    
    def __read_vInf(self,vInf):
        """
        A function that reads the input free stream velocity given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_vInf(vInf)
            
        Parameters
        ----------
        vInf : numpy.array(float)
               The free-stream velocity: vx = `vInf` [0] and vy = `vInf` [1].
                        
        
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """ 
        
        self.__vInf = vInf
        

    def __check_nu(self,nu):
        """
        A function that checks the input viscosity given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__check_nu(nu)

        Parameters
        ----------
        nu : float64
             The dynamic viscosity nu.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # check if nu is a float
        if type(nu) != float:
            raise TypeError('nu must be of type float. It is %s.' % str(type(nu)))

        # check if nu is negative:
        if nu < 0.0:
            raise ValueError('nu must be positive. It is %f.' % nu)


    def __read_nu(self,nu):
        """
        A function that reads the input viscosity given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__read_nu(nu)

        Parameters
        ----------
        nu : float64
             The dynamic viscosity nu.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # save nu
        self.__nu = nu


    def __check_deltaTc(self,deltaTc):
        """
        A function that checks the input convective time step given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__check_deltaTc(deltaTc)

        Parameters
        ----------
        deltaTc : float64
                  The time step to use for the convective steps.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # check if deltaTc is a float
        if type(deltaTc) != float:
            raise TypeError('deltaTc must be of type float. It is %s.' % str(type(deltaTc)))

        # check if deltaTc is negative or zero:
        if deltaTc <= 0.0:
            raise ValueError('deltaTc must be positive and larger than zero. It is %f.' % deltaTc)


    def __read_deltaTc(self,deltaTc):
        """
        A function that reads the input convective time step given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__read_deltaTc(deltaTc)

        Parameters
        ----------
        deltaTc : float64
                  The time step to use for the convective steps.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # save deltaTc
        self.__deltaTc = deltaTc


    def __check_h(self,h):
        """
        A function that checks the blob cell size given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__check_h(h)

        Parameters
        ----------
        h : float64
            The blob cell size.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # check if h is a float
        if type(h) not in [float, _numpy.float64]:
            raise TypeError('h must be of type float. It is %s.' % str(type(h)))

        # check if h is negative or zero:
        if h <= 0.0:
            raise ValueError('h must be positive and larger than zero. It is %f.' % h)


    def __read_h(self,h):
        """
        A function that reads the blob cel size given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__read_h(h)

        Parameters
        ----------
        h : float64
            The blob cell size.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # save deltaTc
        self.__h = h


    def __check_overlap(self,overlap):
        """
        A function that checks the blob overlap ratio between neighboring blobs given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__check_overlap(overlap)

        Parameters
        ----------
        h : float64
            The overlap ratio.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # check if h is a float
        if type(overlap) != float:
            raise TypeError('overlap must be of type float. It is %s.' % str(type(overlap)))

        # check if h is negative or zero:
        if overlap <= 0.0:
            raise ValueError('overlap must be positive and larger than zero. It is %f.' % overlap)


    def __read_overlap(self,overlap):
        """
        A function that reads the overlap ratio between neighboring blobs given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__read_overlap(overlap)

        Parameters
        ----------
        overlap : float64
                  The overlap ratio.


        :First Added:   2014-01-29
        :Last Modified: 2014-01-29
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # save deltaTc
        self.__overlap = overlap


    def __check_blobParams(self,blobParams):
        """
        A function that checks the input blob parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_blobParams(blobParams)
            
        Parameters
        ----------
        blobParams : dict
                     Contains the blob paramaters:
                     'overlap': float64
                                The overlap ratio between neighboring blobs. It is related to
                                the core size of the blob, `sigma`, and to the spacing 'h' by
                                the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`. 
                     'h': float64
                          The size of the cell associated to the vortex blobs. Corresponds to
                          the minimum spacing between the core of two neighboring cells. It is related to
                          the core size of the blob, `sigma`, and to the spacing 'h' by
                          the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`
        
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """ 
        
        # check if overlap is positive
        if blobParams['overlap'] <= 0.0:
            raise ValueError('overlap must be > 0.0. It is %f.' % blobParams['overlap'])
        
        # check if cell size is positive
        if blobParams['h'] <= 0.0:
            raise ValueError('h must be > 0.0. It is %f.' % blobParams['h'])


    def __read_blobParams(self,blobParams):
        """
        A function that reads the input blob parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_blobParams(blobParams)
            
        Parameters
        ----------
        blobParams : dict
                     Contains the blob paramaters:
                     'overlap': float64
                                The overlap ratio between neighboring blobs. It is related to
                                the core size of the blob, `sigma`, and to the spacing 'h' by
                                the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`. 
                     'h': float64
                          The size of the cell associated to the vortex blobs. Corresponds to
                          the minimum spacing between the core of two neighboring cells. It is related to
                          the core size of the blob, `sigma`, and to the spacing 'h' by
                          the expression :math:`\mathtt{overlap} = \frac{\mathtt{h}}{\mathtt{sigma}}`
        
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """ 
        
        # read overlap
        self.__overlap = blobParams['overlap']
        
        # read blob cell size
        self.__h = blobParams['h']
        
        # compute the core spreading, sigma
        self.__sigma = self.__h/self.__overlap
    

    def __check_timeIntegrationParams(self,timeIntegrationParams):
        """
        A function that checks the input time integrations parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_timeIntegrationParams(timeIntegrationParams)
            
        Parameters
        ----------
        timeIntegrationParams : dict
                                Contains the time integration paramaters:
                                'method': string
                                          ('fe', 'rk4') the type of time integrator
                                          used: 'fe' forward Euler, 'rk4' Runge-Kutta 4th order.

        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """
        
        # integrationMethod
        if timeIntegrationParams['method'] not in blobOptions.time_integrator_options:
            raise ValueError(('timeIntegrationParams[\'method\'] must be one of ' +\
                             '\'%s\','*(len(blobOptions.time_integrator_options)-1) +\
                             '\'%s\'.' + ' It is %s.') % (blobOptions.time_integrator_options + (str(timeIntegrationParams['method']),)))


    def __read_timeIntegrationParams(self,timeIntegrationParams):
        """
        A function that reads the input time integrations parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_timeIntegrationParams(timeIntegrationParams)
            
        Parameters
        ----------
        timeIntegrationParams : dict
                                Contains the time integration paramaters:
                                'method': string
                                          ('fe', 'rk4') the type of time integrator
                                          used: 'fe' forward Euler, 'rk4' Runge-Kutta 4th order.

                                
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """
        
        # save the timeIntegrationParams
        self.__timeIntegrationParams = timeIntegrationParams

    
    def __check_blobControlParams(self,blobControlParams):
        """
        A function that checks the input blob control parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_blobControlParams(blobControlParams)
            
        Parameters
        ----------
        blobControlParams : dict
                            Contains the time integration paramaters:
                            'stepRedistribution': int
                                                  the blob redistribution step frequency.
                            'stepPopulationControl': int
                                                     the population control step frequency.
                            'gThresholdLocal': float64
                                               minimum value of circulation to consider for each vortex blob.
                            'gThresholdGlobal': float64
                                                maximum value of total circulation change allowed by discarding blobs.
                                           
                                
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """

        # stepRedistribution
        if type(blobControlParams['stepRedistribution']) not in (int,long): # stepRedistribution must be an integer
            raise ValueError('blobControlParams[\'stepRedistribution\'] must be an int. It is %s.' % str(type(blobControlParams['stepRedistribution'])))
        elif blobControlParams['stepRedistribution']<0: # stepRedistribution must be positive
            raise ValueError('blobControlParams[\'stepRedistribution\'] must be positive. It is %d.' % blobControlParams['stepRedistribution'])

        # stepPopulationControl
        if type(blobControlParams['stepPopulationControl']) not in (int,long): # stepPopulationControl must be an integer
            raise ValueError('blobControlParams[\'stepPopulationControl\'] must be an int. It is %s.' % str(type(blobControlParams['stepPopulationControl'])))
        elif blobControlParams['stepPopulationControl']<0: # stepPopulationControl must be positive
            raise ValueError('blobControlParams[\'stepPopulationControl\'] must be positive. It is %d.' % blobControlParams['stepPopulationControl'])
        
        # gThresholdLocal
        # check first for the minimum value of circulation to consider
        if blobControlParams['gThresholdLocal'] < 0: # gThresholdLocal must be positive or 0
            raise ValueError('blobControlParams[\'gThresholdLocal\'] must be positive. It is %f.' % blobControlParams['gThresholdLocal'])
        
        # gThresholdGlobal
        # check for the maximum value of total circulation allowed to change in population control
        if blobControlParams['gThresholdGlobal'] < 0: # gThresholdGlobal must be positive or 0
            raise ValueError('blobControlParams[\'gThresholdGlobal\'] must be positive. It is %f.' % blobControlParams['gThresholdGlobal'])

    
    def __read_blobControlParams(self,blobControlParams):
        """
        A function that reads the input blob control parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__read_blobControlParams(blobControlParams)
            
        Parameters
        ----------
        blobControlParams : dict
                           Contains the time integration paramaters:
                           'stepRedistribution': int
                                                 the blob redistribution step frequency.
                           'stepPopulationControl': int
                                                    the population control step frequency.
                           'gThresholdLocal': float64
                                              minimum value of circulation to consider for each vortex blob.
                           'gThresholdGlobal': float64
                                               maximum value of total circulation change allowed by discarding blobs.
                                           
                                
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """
        
        # read popControlParams
        self.__blobControlParams = blobControlParams

        # read step redistribution
        self.__stepRedistribution = blobControlParams['stepRedistribution']
        
        # read step population control
        self.__stepPopulationControl = blobControlParams['stepPopulationControl']
        
        # read local circulation control parameter
        self.__gThresholdLocal = blobControlParams['gThresholdLocal']
        
        # read global circulation control parameter
        self.__gThresholdGlobal = blobControlParams['gThresholdGlobal']
        
    
    def __check_velocityComputationParams(self,velocityComputationParams):
        """
        A function that checks the input velocity computation parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_velocityComputationParams(popControlParams)
            
        Parameters
        ----------
        velocityComputationParams : dict
                                    Contains the time integration paramaters:
                                    'method': string
                                              the type of Biot-Savart solver
                                              'direct' or 'fmm'
                                    'hardware': string
                                                the type of hardware to use
                                                'cpu' or 'gpu'
                                                
                                                              
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """

        """    
        Reviews:  
    
        """
        
        # check first the biot-savart computation option
        if velocityComputationParams['method'] not in blobOptions.biot_savart_options:
            raise ValueError(('velocityComputationParams[\'mehtod\'] must be one of ' +\
                             '\'%s\','*(len(blobOptions.biot_savart_options)-1) +\
                             '\'%s\'.' + ' It is %s.') % (blobOptions.biot_savart_options + (str(velocityComputationParams['method']),)))
        # check the hardware option
        if velocityComputationParams['hardware'] not in blobOptions.hardware_options:
            raise ValueError(('velocityComputationParams[\'hardware\'] must be one of ' +\
                             '\'%s\','*(len(blobOptions.hardware_options)-1) +\
                             '\'%s\'.' + ' It is %s.') % (blobOptions.hardware_options + (str(velocityComputationParams['hardware']),)))

    
    def __read_velocityComputationParams(self,velocityComputationParams):
        """
        A function that checks the input velocity computation parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_velocityComputationParams(popControlParams)
            
        Parameters
        ----------
        velocityComputationParams : dict
                                    Contains the time integration paramaters:
                                    'method': string
                                              the type of Biot-Savart solver
                                              'direct' or 'fmm'
                                    'hardware': string
                                                the type of hardware to use
                                                'cpu' or 'gpu'
                                                
                                                              
        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  
    
        """
        
        # save velocityComputationParams
        self.__velocityComputationParams = velocityComputationParams
    
    
    def __check_diffusionParams(self,diffusionParams):
        """
        A function that checks the input diffusion parameters given to self.__init__.
        
        Usage
        -----
        .. code-block:: python
        
            self.__check_diffusionParams(diffusionParams)
            
        Parameters
        ----------
        diffusionParams : dict
                          Contains the time integration parameters:
                          'method': string
                                    the type of diffusion scheme to use, can be
                                    one of:
                                    defining the diffusion method to use and
                                    'regrid_wee' (based in [1]_)
                                    'regrid_tutty' (based in [2]_).
                          Parameters for each different method:
                              'regrid_wee':
                                  `c2` : (float64 or 'optimized')
                                         free parameter that specifies the stability
                                         of the diffusion scheme. Best values: [1/6,1/2]
                                         if 'optimized' is given, 1/3 is selected.

                              'regrid_tutty'
        
        .. [1] Wee, D., Ahmed, F., Ghoniem, F., Modified interpolation
              kernels for treating diffusion and remeshing in vortex
              methods, Journal of Computational Physics, 213, 239--263, 2006.
        .. [2] Tutty, O.R., A simple redistribution vortex method (with accurate body forces),
                   arXiv:1009.0166v1, 2010
                                         
        :First Added:   2014-01-27
        :Last Modified: 2014-04-14
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version
        
        """
    
        """    
        Reviews:  1- Added regrid_tutty diffusion scheme.
    
        """
        
        # check first the diffusion method used
        if diffusionParams['method'] not in blobOptions.diffusion_method_options:
            raise ValueError(('diffusionParams[\'method\'] must be one of ' +\
                              '\'%s\','*(len(blobOptions.diffusion_method_options)-1) +\
                              '\'%s\'.' + ' It is %s.') % (blobOptions.diffusion_method_options + (str(diffusionParams['method']),)))

        # check the specific parameters for each diffusion method

        # regrid_wee
        if diffusionParams['method'] == 'regrid_wee':
            if 'c2' not in diffusionParams.keys():
                raise ValueError('c2 parameter must be defined for diffusion method regrid_wee in diffusionParams')
            # c2 parameter must be a float between 1/6 and 1/2 or take the string value 'optimized'
            if type(diffusionParams['c2']) != float:
                if diffusionParams['c2'] != 'optimized':
                    raise ValueError('c2 parameter must be a float or \'optimized\'. You introduced %s.' % str(diffusionParams['c2']))
            elif (diffusionParams['c2'] < 1.0/6.0) or (diffusionParams['c2'] > 1.0/2.0):
                    raise ValueError('c2 parameter must be between 1/6 and 1/2. You introduced %f.' % diffusionParams['c2'])

        if diffusionParams['method'] == 'regrid_tutty':
            pass


    def __read_diffusionParams(self,diffusionParams):
        """
        A function that reads the input diffusion parameters given to self.__init__.

        Usage
        -----
        .. code-block:: python

            self.__read_diffusionParams(diffusionParams)

        Parameters
        ----------
        diffusionParams : dict
                          Contains the time integration paramaters:
                          'method': string
                                    the type of diffusion scheme to use, can be
                                    one of:
                                    defining the diffusion method to use and
                                    'regrid_wee' (based in [1]_).
                          Parameters for each different method:
                              'regrid_wee':
                                  `c2` : (float64)
                                         free parameter that specifies the stability
                                         of the diffusion scheme. Best values: [1/6,1/2].

        .. [1] Wee, D., Ahmed, F., Ghoniem, F., Modified interpolation
              kernels for treating diffusion and remeshing in vortex
              methods, Journal of Computational Physics, 213, 239--263, 2006.

        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        """
        Reviews:

        """

        # check the diffusion method used and then read its specific parameters

        # regrid_wee
        if diffusionParams['method'] in blobOptions.diffusion_method_options:
            self.__diffusionParams = diffusionParams
        else:
            raise ValueError('No other diffusion method have been implemented yet! Use only \'regrid_wee\'.')


    def __compute_deltaTd(self):
        """
        A function that computes the diffusive time step (deltaTd) and the frequency of diffusion computation.
        Diffusion step will be a multiple of the convective step.

        The ideal diffusion time step will be computed and then readjusted to the closest multiple of the
        convective time step. Additionally, the multiplicative constant between the convective and the diffusive
        steps is also computed (stepDiffusion).

        Usage
        -----
        .. code-block:: python

            deltaTd,stepDiffusion = self.__compute_deltaTd()

        Parameters
        ----------

        .. [1] Wee, D., Ahmed, F., Ghoniem, F., Modified interpolation
              kernels for treating diffusion and remeshing in vortex
              methods, Journal of Computational Physics, 213, 239--263, 2006.

        :First Added:   2014-01-27
        :Last Modified: 2014-01-27
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        # compute the ideal deltaTd for the diffusion method selected
        if self.__diffusionParams['method'] == 'regrid_wee': # regrid_wee diffusion scheme of [1]
            if self.__nu == 0.0: # if there is no viscosity then deltaTd is None
                deltaTd = 0.0
            else:
                if self.__diffusionParams['c2'] == 'optimized':
                    deltaTd = ((1.0/3.0)*self.__h*self.__h)/self.__nu
                else:
                    deltaTd = (self.__diffusionParams['c2']*self.__h*self.__h)/self.__nu

            # readjust deltaTd to be a multiple of self.__deltaTc and compute the multiplicity constant between convection
            # and diffusion: stepDiffusion
            if deltaTd >= self.deltaTc:
                stepDiffusion = round(deltaTd/self.deltaTc)
                deltaTd = stepDiffusion*self.deltaTc
            elif deltaTd == 0.0:
                stepDiffusion = 0
                deltaTd = 0.0
            else:
                raise ValueError('deltaTc is too large. Reduce deltaTc to a value smaller than %f.' % deltaTd)


        if self.__diffusionParams['method'] == 'regrid_tutty': #regrid_tutty diffusion scheme of [2]
            if self.__nu == 0.0:
                deltaTd = 0.0
                stepDiffusion = 0
            else:
                deltaTd = self.deltaTc
                stepDiffusion = 1

                # characteristic diffusion scale
                hNu = _numpy.sqrt(deltaTd*self.__nu)

                # check if deltaTd satisfies the stability condition
                if ((hNu/self.__h) >= (1.0/_numpy.sqrt(2.0))):

                    deltaTd = (((1.0/_numpy.sqrt(2.0))*self.__h)**2)/self.__nu

                    raise ValueError('deltaTc is toooo large. Reduce deltaTc to a value smaller than %f.' % deltaTd)



        return deltaTd, stepDiffusion

    #----------------------------------------------------
    # Definition of properties

    @property
    def stepRedistribution(self):
        """
        `stepRedistribution` (int) the redistribution step frequency.
            
        Usage
        -----
        .. code-block:: python
        
        self.stepRedistribution
        """
        return self.__stepRedistribution

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
        return self.__timeIntegrationParams['method']
    
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
        return self.__velocityComputationParams['method']
    
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
        return self.__stepPopulationControl
    
    @stepPopulationControl.setter
    def stepPopulationControl(self,stepPopulationControlnew):
        raise AttributeError('stepPopulationControl cannot be manually set by the user, user can only get the value of stepPopulationControl.')
            
    @stepPopulationControl.deleter
    def stepPopulationControl(self):
        raise AttributeError('stepPopulationControl cannot be manually deleted by the user, user can only get the value of stepPopulationControl.')

    @property
    def gThresholdGlobal(self):
        """
            'gThresholdGlobal' with maximum value of total circulation change allowed when performing population control.
            
            Usage
            -----
            .. code-block:: python
        
                self.gThresholdGlobal
        """
        return self.__gThresholdGlobal
    
    @gThresholdGlobal.setter
    def gThresholdGlobal(self,gThresholdGlobal):
        raise AttributeError('gThresholdGlobal cannot be manually set by the user, user can only get the value of gThresholdGlobal.')
            
    @gThresholdGlobal.deleter
    def gThresholdGlobal(self):
        raise AttributeError('gThresholdGlobal cannot be manually deleted by the user, user can only get the value of gThresholdGlobal.')


    @property
    def gThresholdLocal(self):
        """
            'gThresholdLocal' with minimum value of circulation to
            consider for each vortex blob.

            Usage
            -----
            .. code-block:: python

                self.gThresholdLocal
        """
        return self.__gThresholdLocal

    @gThresholdLocal.setter
    def gThresholdLocal(self,gThresholdLocal):
        raise AttributeError('gThresholdLocal cannot be manually set by the user, user can only get the value of gThresholdLocal.')

    @gThresholdLocal.deleter
    def gThresholdLocal(self):
        raise AttributeError('gThresholdLocal cannot be manually deleted by the user, user can only get the value of gThresholdLocal.')


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
        """
            Set a new free stream velocity

            Usage
            -----
            .. code-block:: python

                self.vInf = vInfnew

            Parameters
            ----------
            vInfnew : numpy.array(float)
                      The free-stream velocity: vx = `vInfnew` [0] and vy = `vInfnew` [1].


        """

        # first check if vInf is a good value
        self.__check_vInf(vInfnew)

        # update the free stream velocity
        self.__read_vInf(vInfnew)

            
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
        
        
    # set the property g in order to access self.__y in a safe way
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
        
    # set the property numBlobs in order to get the number of blobs
    @property
    def numBlobs(self):
        """
            The number of blobs.
            
            Usage
            -----
            .. code-block:: python
        
                self.num_blobs

        """
        return self.__x.shape[0]
        
    @numBlobs.setter
    def numBlobs(self,gnew):
        raise AttributeError('numBlobs cannot be set by the user, user can only get the value of numBlobs.')
            
    @numBlobs.deleter
    def numBlobs(self):
        raise AttributeError('numBlobs cannot be deleted by the user, user can only get the value of numBlobs.')
        
    # set the property diffusionParams in order to get the diffusion parameters
    # of the simulation
    @property      
    def diffusionParams(self):
        """
            The diffusion parameters. It is a dictionary containing all the
            parameters of the diffusion method used for the simulation.
            
            Usage
            -----
            .. code-block:: python
                
                self.diffusionParameters
                
        """
        
        return self.__diffusionParams
        
    @diffusionParams.setter
    def diffusionParams(self,diffusionParamNew):
        raise AttributeError('diffusionParam cannot be set by the user, user can only get the value of diffusionParam.')
            
    @diffusionParams.deleter
    def diffusionParams(self):
        raise AttributeError('diffusionParam cannot be deleted by the user, user can only get the value of diffusionParam.')


    # set the property timeIntegrationParams in order to get the time integration parameters
    # of the simulation
    @property
    def timeIntegrationParams(self):
        """
            The diffusion parameters. It is a dictionary containing all the
            parameters of the diffusion method used for the simulation.

            Usage
            -----
            .. code-block:: python

                self.diffusionParameters

        """

        return self.__timeIntegrationParams

    @timeIntegrationParams.setter
    def timeIntegrationParams(self,timeIntegrationParamsNew):
        raise AttributeError('timeIntegrationParams cannot be set by the user, user can only get the value of timeIntegrationParams.')

    @timeIntegrationParams.deleter
    def timeIntegrationParams(self):
        raise AttributeError('timeIntegrationParams cannot be deleted by the user, user can only get the value of timeIntegrationParams.')


    # set the property blobControlParams in order to get the blob control parameters
    # of the simulation
    @property
    def blobControlParams(self):
        """
            The diffusion parameters. It is a dictionary containing all the
            parameters of the diffusion method used for the simulation.

            Usage
            -----
            .. code-block:: python

                self.diffusionParameters

        """

        return self.__blobControlParams

    @blobControlParams.setter
    def blobControlParams(self,blobControlParamsNew):
        raise AttributeError('blobControlParams cannot be set by the user, user can only get the value of blobControlParams.')

    @blobControlParams.deleter
    def blobControlParams(self):
        raise AttributeError('blobControlParams cannot be deleted by the user, user can only get the value of blobControlParams.')


    # set the property velocityComputationParams in order to get the velocity computation parameters
    # of the simulation that defines if direct or fmm calculation should be performed using the gpu or cpu
    @property
    def velocityComputationParams(self):
        """
            The diffusion parameters. It is a dictionary containing all the
            parameters of the diffusion method used for the simulation.

            Usage
            -----
            .. code-block:: python

                self.diffusionParameters

        """

        return self.__velocityComputationParams

    @velocityComputationParams.setter
    def velocityComputationParams(self,velocityComputationParamsNew):
        raise AttributeError('velocityComputationParams cannot be set by the user, user can only get the value of velocityComputationParams.')

    @velocityComputationParams.deleter
    def velocityComputationParams(self):
        raise AttributeError('velocityComputationParams cannot be deleted by the user, user can only get the value of velocityComputationParams.')


    # set the property tStep in order to access self.__tStep in a safe way
    @property
    def tStep(self):
        """
            The current time step of the simulation.

            Usage
            -----
            .. code-block:: python

                self.tStep

        """
        return self.__tStep

    @tStep.setter
    def tStep(self,tStepnew):
        raise AttributeError('tStep cannot be manually set by the user, user can only get the value of tStep.')

    @tStep.deleter
    def tStep(self):
        raise AttributeError('tStep cannot be manually deleted by the user, user can only get the value of tStep.')


    # set the property t in order to access self.__t in a safe way
    @property
    def t(self):
        """
            The current time of the simulation.

            Usage
            -----
            .. code-block:: python

                self.t

        """
        return self.__t

    @t.setter
    def t(self,tnew):
        raise AttributeError('t cannot be manually set by the user, user can only get the value of t.')

    @t.deleter
    def t(self):
        raise AttributeError('t cannot be manually deleted by the user, user can only get the value of t.')


    # set the property stepDiffusion in order to access self.__stepDiffusion in a safe way
    @property
    def stepDiffusion(self):
        """
            The frequency of the diffusion step.

            Usage
            -----
            .. code-block:: python

                self.stepDiffusion

        """
        return self.__stepDiffusion

    @stepDiffusion.setter
    def stepDiffusion(self,stepDiffusionnew):
        raise AttributeError('stepDiffusion cannot be manually set by the user, user can only get the value of stepDiffusion.')

    @stepDiffusion.deleter
    def stepDiffusion(self):
        raise AttributeError('stepDiffusion cannot be manually deleted by the user, user can only get the value of stepDiffusion.')


    # set the property hardware that defines the hardware to use to compute the induced velocities
    @property
    def hardware(self):
        """
            The hardware parameter. It is a string that defines if the CPU or GPU should be used
            for computing the induced velocities.

            Usage
            -----
            .. code-block:: python

                self.hardware

        """

        return self.__velocityComputationParams['hardware']

    @hardware.setter
    def hardware(self,hardwarenew):
        """
            Set a new hardware parameter

            Usage
            -----
            .. code-block:: python

                self.hardware = hardwarenew

            Parameters
            ----------

            hardwarenew : string
                          The new value for the hardware. It can be one of:
                          'cpu' for using the CPU
                          'gpu' for using the GPU

        """

        # set the new parameter
        self.__velocityComputationParams['hardware'] = hardwarenew

        # check if the new velocityComputationParams is allowed
        self.__check_velocityComputationParams(self.__velocityComputationParams)


    @hardware.deleter
    def hardware(self):
        raise AttributeError('hardware cannot be deleted by the user, user can only get and set the value of hardware.')


    # set the property plotVelocity that defines if velocity is to be plotted at the blobs
    @property
    def plotVelocity(self):
        """
            The flag that defines if the induced velocity field at the blobs is to be plotted.

            Usage
            -----
            .. code-block:: python

                self.plotVelocity

        """

        return self.__plotVelocity

    @plotVelocity.setter
    def plotVelocity(self,plotVelocitynew):
        """
            Set a the plotVelocity flag

            Usage
            -----
            .. code-block:: python

                self.plotVelocity = plotVelocitynew

            Parameters
            ----------

            plotVelocityNew : bool
                              The new value for the plotVelocity flag.

        """

        if type(plotVelocitynew) != bool:
            raise TypeError('plotVelocity must be of type bool. It is of type: %s' % str(type(plotVelocitynew)))

        # update the flag
        self.__plotVelocity = plotVelocitynew

    @plotVelocity.deleter
    def plotVelocity(self):
        raise AttributeError('plotVelocity cannot be deleted by the user, user can only get and set the value of plotVelocity.')