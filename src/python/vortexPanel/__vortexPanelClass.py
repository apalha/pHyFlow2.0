#-*- coding: utf-8 -*-
__doc__ = """

VORTEXPANELCLASS
================

The main class for the vortex-panel coupled.

* Note: No Turbulence scheme implemented (only laminar flow)

Description
-----------
This class wraps the vortex blobs and panel together.

Methodology
-----------


Implemented NS Algorithms
-------------------------


:First Added:   2014-02-21
:Last Modified: 2014-02-21                         
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

__all__ = ['VortexPanel']

# External packages
import numpy as _numpy
import sys

# Import pHyFlow
from pHyFlow import options
from pHyFlow.vortexPanel import base as _base
from pHyFlow.vortexPanel import vpOptions as _vpOptions

from pHyFlow.vortex import VortexBlobs as _VortexBlobs
from pHyFlow.vortex.base.induced import velocity as _blobs_velocity
from pHyFlow.panel import Panels as _Panels

   

class VortexPanel(object):
    r"""
    
    Usage
    -----
    .. code-block:: python
    
        VortexPanel(blobs, panels, 
                    couplingParams={})
        
    Parameters
    ----------
    blobs : pHyFlow.vortex.VortexBlobs
            the vortex blob class containing all the parameters that define
            the vortex-blob problem
    
    panels : pHyFlow.panels.Panels
             the panel class containing all the parameters that define the
             panel problem.
    
    couplingParams : dict, optional
                     Dictionary containing parameters for coupling the
                     panels and the blobs together.
                     
                     'panelStrengthUpdate' : str
                                             String parameters specifying if
                                             panel strength should be updated
                                             during the mult-step time integration
                                             scheme (such as RK4).
                                             
                                             'constant' : solve for panel strength
                                                          only in the beginning.
                                              'varying' : constantly update
                                                          the panel strength depending
                                                          on the sub-step position
                                                          of the particles.
                          
    Attribute
    ---------
    blobs
    panels
    __couplingParams     
    
    Methods
    -------
    evolve
    __coupled_convection
   

    :First Added:   2014-02-21
    :Last Modified: 2014-02-24                         
    :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
    :License:       GNU GPL version 3 or any later version   
   
    """
    """
    Revisions
    ---------


    """
    
    def __init__(self,blobs,panels,couplingParams={
                 'panelStrengthUpdate': _vpOptions.PANEL_STRENGTH_UPDATE['default']}):
        
        
        #---------------------------------------------------------------------
        # Check/Set input parameters
            
        # Initialize the parameters
        self.__set('blobs',blobs) # vortex-blobs class
        self.__set('panels', panels) # panel class
        
        # Initialize couplingParams
        self.__set('couplingParams',couplingParams)
            
        #---------------------------------------------------------------------            
        # Make references

        # Modify the parameters
        self.blobs.evolve = self.__mute_blobs_evolve # remove blobs.evolve option
       
    def evaluateVelocity(self,xTarget,yTarget):
        """
        Function to evaluate the total induced velocity due to vortex blobs,
        panels and the free-stream flow.
        .. math::
        
            \mathbf{u} = \mathbf{u}_{\omega} + \mathbf{u}_{\gamma} + \mathbf{u}_\infty
            
        Usage
        -----
        .. code-block :: python
        
            vx, vy = evaluteVelocity(xTarget, yTarget)
            
        Parameters
        ----------
        xTarget : numpy.ndarray(float64), shape (nTargets,)
                  the :math:`x`-coordinate of the target location, where the
                  total velocity is to be evaluated
                  
        yTarget : numpy.ndarray(float64), shape (nTargets,)
                  the :math:`y`-coordinate of the target location, where the
                  total velocity is to be evaluated
          
        Returns
        -------
        vx : numpy.ndarray(float64), shape (nTarget,)
             the :math:`x`-component of the total induced velocity on each
             **(xTarget, yTarget)** point.

        vy : numpy.ndarray(float64), shape (nTarget,)
             the :math:`y`-component of the total induced velocity on each
             **(xTarget, yTarget)** point.

        Attributes
        ----------
        None changed.
        
        :First Added:   2014-02-24
        :Last Modified: 2014-02-24                         
        :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
        :License:       GNU GPL version 3 or any later version   
        """
        # Determine the induced velocity due to blobs + panels + free-stream.
        vx,vy =self.blobs.evaluateVelocity(xTarget,yTarget) \
                        + self.panels.evaluateVelocity(xTarget,yTarget) \
                        + self.vInf.reshape(2,-1)
        
        # return the induced velocity
        return vx,vy
                
    def evolve(self):
        r"""
        Function to evolve the coupled vortex-panel problem. During the evolution
        of the vortex blobs, the no-through boundary condition is taken into 
        account by solving the no-slip panel problem.
        
        The problem can be solved in two manners: 
                
            panelStrengthUpdate = 'constant'
            
                In this case, the no-slip b.c is only solved at the beginning
                of the multi-stepping process. The solution panel problem is 
                then assumed to stay constant during the sub-stepping of the
                vortex blobs. 
                
            panelStrengthUpdate = 'varying'
            
                In this case, the panel problem is solved in every sub-step.
                This ensures that that no-slip b.c is also satisfied during
                the sub-step.
                        
        Usage
        -----
        .. code-block:: python

            evolve()

        Parameters
        ----------
        None.
        
        Returned
        --------
        None returned.

        Attributes
        ----------
        blobs : pHyFlow.vortex.VortexBlobs
                the position and the strength of the blobs are updated
                
                x : numpy.ndarray(float64), shape (blobs.numBlobs, )
                    the :math:`x`-position of the particle.
                y : numpy.ndarray(float64), shape (blobs.numBlobs, )
                    the :math:`y`-postion of the particles.
                g : numpy.ndarray(float64), shape (blobs.numBlobs, )
                    the circulation :math:`\Gamma` of the particles.
                    
        panels : pHyFlow.panels.Panels
                 the strength of the panels are updated.
                 
                 sPanel : numpy.ndarray(float64), shape (panels.nPanelTotal, )
                          the strength :math:`\gamma` of the panels.
                  
                 * Note: position update not implemented yet.
         
        t : float
            the current time
            
        tStep : int
                the current time step
         
        :First Added:   2014-02-21
        :Last Modified: 2014-02-24
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """        
        #----------------------------------------------------------------------
        # Evolve : Convect + Diffuse the blobs with regard of panel geometry
      
        # convect blobs with coupled panels
        self.__coupled_convection() # Algorithm for convection
        
        # Diffusion Step : diffuse vortex blobs
        if self.blobs.stepDiffusion != 0: # corresponds to nu = 0.0, therefore no diffusion is to be computed
            if (self.tStep % self.blobs.stepDiffusion) == 0: # if the time step is a multiple of the stepDiffusion perform diffusion
                self.blobs._diffusion()

        # update the time counter
        self._advanceTime()

        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        # Redistribution step
        if self.blobs.stepRedistribution != 0: # meaning that redistribution should be done
            if (self.tStep % self.blobs.stepRedistribution) == 0: # if the time step is a multiple of the stepRedistribution perform redistribution
                self.blobs.redistribute()
                
        #----------------------------------------------------------------------                

        #---------------------------------------------------------------------- 
        # population control step
        if self.blobs.stepPopulationControl != 0: # meaning that population control should be done
            if (self.tStep % self.blobs.stepPopulationControl) == 0: # if the time step is a multiple of the stepPopulationControl perform population control
                self.blobs.populationControl()
                
        #----------------------------------------------------------------------                
            
            
            
    def __coupled_convection(self):
        """
        Internal function to convect the blobs with taking account of the 
        inviscid geometry (panel body).
        
        During the time-stepping, the no-through boundary condition is taken  
        into account by solving the no-slip panel problem.
        
        The problem can be solved in two manners: 
                
            panelStrengthUpdate = 'constant'
            
                In this case, the no-slip b.c is only solved at the beginning
                of the multi-stepping process. The solution panel problem is 
                then assumed to stay constant during the sub-stepping of the
                vortex blobs. 
                
            panelStrengthUpdate = 'varying'
            
                In this case, the panel problem is solved in every sub-step.
                This ensures that that no-slip b.c is also satisfied during
                the sub-steps.
                
        Usage
        -----
        .. code-block :: python
        
            __coupled_convection()
            
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        
        :First Added:   2014-02-21
        :Last Modified: 2014-02-24
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version        
        
        """      

        #----------------------------------------------------------------------
        # Determine parameters

        # convert the hardware flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['hardware'] == 'gpu': 
            blobs_hardware = options.GPU_HARDWARE
        else: 
            blobs_hardware = options.CPU_HARDWARE

        # convert the method flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['method'] == 'fmm': 
            blobs_method = options.FMM_METHOD
        else: 
            blobs_method = options.DIRECT_METHOD
    
        # convert the time integration method into an int to use in _base_convection
        if self.blobs.timeIntegrationParams['method'] == 'rk4': 
            blobs_integrator = options.RK4_INTEGRATOR
        elif self.blobs.timeIntegrationParams['method'] == 'euler':
            blobs_integrator = options.FE_INTEGRATOR
            
        #----------------------------------------------------------------------


        #----------------------------------------------------------------------
        # Time integrate the vortex blobs

        # If integration method is RK4
        if blobs_integrator == options.RK4_INTEGRATOR:

            # Make references to vortex-blobs
            xBlob, yBlob, gBlob = self.blobs.x, self.blobs.y, self.blobs.g 
            
            # Make references to panel collocation points (where no-slip b.c. is enforced.)
            xCP, yCP = self.panels.xCPGlobalCat, self.panels.yCPGlobalCat
        
            # Runge-Kutta Butcher Tablaeu:
            #  
            #   - Parameters according to the Butcher tableau of RK4: 
            #   
            #  (c_i)
            #    |   0  |       (a_ij)
            #    v  1/2 | 1/2            
            #       1/2 |  0    1/2 
            #        1  |  0     0    1
            #      ------------------------- 
            #           | 1/6   1/3  1/3   1     (b_j) ->
            #
            
            # Step-up for time-integration.
            kxTemp = 0
            kyTemp = 0
        
            # For-loop for the RK4: Step 1 to 4
            for stage, RK4_b,RK4_c in zip([1,2,3,4],
                                          [1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0], # b coefficients
                                          [0.0, 0.5, 0.5, 1.0]):                # c coefficients
            
                # Determine the temporary blob position
                xBlobTemp, yBlobTemp = xBlob + RK4_c*self.deltaT*kxTemp, yBlob + RK4_c*self.deltaT*kyTemp
            
                #--------------------------------------------------------------
                # Solve panel strength
            
                # Check if panel strength should be updated.
                if self.__couplingParams['panelStrengthUpdate'] == 'varying':
                    # In the case of varying panel strength, solve for panel
                    # strength at each sub-steps
                
                    # Determine induced velocity of blobs on panels [blobs --> panels]
                    # * Note: don't forget the free-stream flow.
                    vxSlip, vySlip = _blobs_velocity(xBlobTemp,yBlobTemp,gBlob,self.blobs.sigma,
                                                     xEval=xCP,yEval=yCP,hardware=blobs_hardware,   
                                                     method=blobs_method) \
                                             + self.vInf.reshape(2,-1)
                
                    # Solve for no-slip panel strengths
                    self.panels.solve(vxSlip, vySlip)
                    
                # In the case of constant panels strength, evaluate only once.
                elif self.__couplingParams['panelStrengthUpdate'] == 'constant':
                    # If constant, only calculate panel strengt at first substep                        
                    if stage==1:
                        
                        # Determine induced velocity of blobs on panels [blobs --> panels]
                        # * Note: don't forget the free-stream flow.
                        vxSlip, vySlip = _blobs_velocity(xBlobTemp, yBlobTemp, gBlob, self.blobs.sigma,
                                                         xEval=xCP, yEval=yCP, hardware=blobs_hardware,   
                                                         method=blobs_method) \
                                                 + self.vInf.reshape(2,-1)
                    
                        # Solve for no-slip panel strengths
                        self.panels.solve(vxSlip, vySlip) 

                #--------------------------------------------------------------

                #--------------------------------------------------------------
                # Calculate the induced velocities on the blobs [ --> blobs] 
                   
                # Determine the induced velocity of panels on blobs [panels -- > blobs]
                vxPanel, vyPanel = self.panels.evaluateVelocity(xBlobTemp,yBlobTemp)
                    
                # Determine the self-induction velocity of blobs [blobs --> blobs]
                vxBlob, vyBlob = _blobs_velocity(xBlobTemp,yBlobTemp,gBlob,self.blobs.sigma,\
                                                 xEval=None,yEval=None,
                                                 hardware=blobs_hardware,
                                                 method=blobs_method)
                                                 
                #--------------------------------------------------------------                                                 
                                                 
                #--------------------------------------------------------------
                # Perform the time-integration
                                                 
                # compute kx and ky                              
                kxTemp = vxBlob + vxPanel + self.vInf[0] # compute kx1
                kyTemp = vyBlob + vyPanel + self.vInf[1] # compute ky1
                
                # update the positions using
                self.blobs._VortexBlobs__x += RK4_b*kxTemp*self.deltaT
                self.blobs._VortexBlobs__y += RK4_b*kyTemp*self.deltaT                                                      
                #--------------------------------------------------------------
            
            
        elif blobs_integrator == options.FE_INTEGRATOR:
            raise NotImplementedError('FE time integration not available !')
        
        #----------------------------------------------------------------------           

    
    def _advanceTime(self):
        """
        Function to advance time.
        
        Usage
        -----
        .. code-block :: python
        
            _advanceTime()
            
        Parameters
        ----------
        None.
        
        Returns
        -------
        None returned.
        
        Attributes
        ----------
        t : float
            the current time
            
        tStep : int
                the current time step
        
        :First Added:   2014-02-24
        :Last Modified: 2014-02-24
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
        """
        # The time-stepping algorithm is found inside blobs.
        self.blobs._advanceTime()
    
            
    def __set(self,varName,var):
        """
        Function to set all the parameters.
        
        Usage
        -----
        .. code-block :: python
        
            __set('varName', var)
            
        Parameters
        ----------
        varName : str
                  the variable name that is to be assigned. Parameters that
                  are set:

        var : <var>
              the variable that is to be set.

        Returns              
        -------
        None returned.
        
        Attributes
        ----------
        blobs : pHyFlow.vortex.VortexBlobs
                the vortex-blobs class is assigned
                
        panels : pHyFlow.panels.Panels
                 the panel blass
                 
        couplingParams : dict
                         Dictionary containing parameters for coupling the
                         panels and the blobs together.
                         
        :First Added:   2014-02-24
        :Last Modified: 2014-02-24
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
        """
        # Vortex-blobs
        if varName == 'blobs':
            if type(var) != _VortexBlobs:
                raise ValueError("'blobs' should be of type pHyFlow.vortex.VortexBlobs. It is %s" % type(var))
                
            self.blobs = var
            self.blobs.evolve = self.__mute_blobs_evolve
        
        # panels
        elif varName == 'panels':
            if type(var) != _Panels:
                raise ValueError("'blobs' should be of type pHyFlow.vortex.VortexBlobs. It is %s" % type(var))
                
            self.panels = var
            
        # coupling parameters
        elif varName == 'couplingParams':
            if type(var) != dict:
                raise ValueError("'couplingParams' should be of type dict. It is %s" % type(var))
            if var['panelStrengthUpdate'] not in _vpOptions.PANEL_STRENGTH_UPDATE['available']:
                raise ValueError("'couplingParams['panelStrengthUpdate']' is unknown. It should be in [%s]. It is %s" % (str(_vpOptions.PANEL_STRENGTH_UPDATE['available']),var['panelStrengthUpdate']))
            self.__couplingParams = var
            
            
    def __mute_blobs_evolve(self):
        raise AttributeError('Not possible with vortexPanel coupled class !')
            
    def __setError(self,setVar):
        raise AttributeError('Cannot be manually set !')
        
    def __delError(self):
        raise AttributeError('Cannot be manually deleted !')         
    #--------------------------------------------------------------------------
    # Define properties
            
    # Time-step
    #    test = my_property(getter = lambda self: self.blobs.tStep,
    #                       doc = "Documentation...")
        
    tStep = property(fget = lambda self: self.blobs.tStep,
                     fset = __setError,
                     fdel = __delError,
                     doc  = r"""tStep : int
                                        the current time step 
                     """)
            
    # Current time t
    t = property(fget = lambda self: self.blobs.t)
    
    # Time-step size
    deltaT= property(fget = lambda self: self.blobs.deltaTc)
    
    # Free-stream velocity
    vInf  = property(fget = lambda self: self.blobs.vInf)