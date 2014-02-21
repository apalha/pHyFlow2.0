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
        
    Parameters
    ----------

    Attribute
    ---------
   
    Methods
    -------
   

    :First Added:   2014-02-21
    :Last Modified: 2014-02-21                         
    :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
    :License:       GNU GPL version 3 or any later version   
   
    """
    """
    Revisions
    ---------


    """
    
    def __init__(self,blobs,panels):
        
        
        #---------------------------------------------------------------------
        # Check/Set input parameters
            
        # Initialize the parameters
        self.__set('blobs',blobs)
        self.__set('panels', panels)
            
        

    def evolve(self):
        r"""
        

        Usage
        -----
        .. code-block :: python

            evolve()

        Parameters
        ----------

        Attributes
        ----------
            

        :First Added:   2014-02-21
        :Last Modified: 2014-02-21
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version

        """
        
        #------------------------------------------------------------------
        # convert the hardware flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['hardware'] == 'gpu': blobs_hardware = options.GPU_HARDWARE
        else: blobs_hardware = options.CPU_HARDWARE

        # convert the method flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['method'] == 'fmm': blobs_method = options.FMM_METHOD
        else: blobs_method = options.DIRECT_METHOD
        
        
        # convert the time integration method into an int to use in _base_convection
        if self.blobs.timeIntegrationParams['method'] == 'rk4': blobs_integrator = options.RK4_INTEGRATOR
        elif self.blobs.timeIntegrationParams['method'] == 'euler': blobs_integrator = options.FE_INTEGRATOR

        #------------------------------------------------------------------        
        
        # TODO:
        
        # convect blobs with coupled panels
        self.__coupled_convect()
        
        # diffuse vortex blobs
        self.blobs._diffuse()

        # update the time counter
        self.blobs._advanceTime()

        # redistribution step
        if self.blobs.stepRedistribution != 0: # meaning that redistribution should be done
            if (self.blobs.tStep % self.blobs.stepRedistribution) == 0: # if the time step is a multiple of the stepRedistribution perform redistribution
                self.blobs.redistribute()

        # population control step
        if self.blobs.stepPopulationControl != 0: # meaning that population control should be done
            if (self.blobs.tStep % self.blobs.stepPopulationControl) == 0: # if the time step is a multiple of the stepPopulationControl perform population control
                self.populationControl()
            
            
            
        def __varyingPanelStrength_coupled_convection_rk4(self):
            """
            coupled convection.
            """      

    
    
            #            # convection step
            #            xBlobTemp,yBlobTemp,gBlobTemp = _base_coupled_convection(self.__deltaTc,self.__x,self.__y,self.__g,self.__sigma,
            #                                                                     hardware=hardware,method=method,integrator=integrator,vInf=self.__vInf)

            # Make references

            # blobs
            xBlob, yBlob, gBlob = self.blobs.x, self.blobs.y, self.blobs.g    
        
            # panels
            xCP, yCP = self.panels.xCPGlobalCat, self.panels.yCPGlobalCat
            
            # free-stream
            vInf = self.blobs.vInf
            
            # time-step
            dt = self.blobs.deltaTc
            
            # allocate memory space for new blobs
            xBlobNew = _numpy.zeros(xBlob.shape)
            yBlobNew = _numpy.zeros(yBlob.shape)
            
            kxTemp = 0
            kyTemp = 0
            
            for k in [1.0, 2.0, 2.0, 1.0]:
                
                # Determine the temporary blob position
                xBlobTemp, yBlobTemp = xBlob + 0.5*dt*kxTemp, yBlob + 0.5*dt*kyTemp
            
                #TODO: we are here
            #------------------------------------------------------------------
            # STEP 1 :: x_new_1 = k1 (sub-step 1)
    
            # evaluate the velocities at/from (xBlob,yBlob)
    
            # Evalute velocity [ blobs --> ]

            # Determine induced velocity of blobs on panels [blobs --> panels]
            # * Note: don't forget the free-stream flow.
            vxSlip, vySlip = _blobs_velocity(xBlob,yBlob,gBlob,self.blobs.sigma,\
                                             xEval=xCP,yEval=yCP,
                                             hardware=blobs_hardware,
                                             method=blobs_method) + vInf.reshape(2,-1)
            
            # Solve for no-slip panel strengths
            self.panels.solve(vxSlip, vySlip)
            
            # Evaluate velocity [  --> blobs]
            
            # Determine the induced velocity of panels on blobs [panels -- > blobs]
            vxPanel, vyPanel = self.panels.evaluteVelocity(xBlob,yBlob)
            
            
            # Determine the self-induction velocity of blobs [blobs --> blobs]
            vxBlob,vyBlob = _blobs_velocity(xBlob,yBlob,gBlob,self.blobs.sigma,\
                                            xEval=None,yEval=None,
                                            hardware=blobs_hardware,
                                            method=blobs_method)
                                    
            # compute kx1 and ky1 using (2)                                    
            kxTemp = vxBlob + vxPanel + vInf[0] # compute kx1 using (2)
            kyTemp = vyBlob + vyPanel + vInf[1] # compute ky1 using (2)
        
            # update the positions using (step 1)
            xBlobNew += kxTemp # kxTemp is now kx1
            yBlobNew += kyTemp # kyTemp is now ky1
            
            #------------------------------------------------------------------
            
            
            #------------------------------------------------------------------
            # STEP 2 :: x_new_2 = x_new_1 + 2*k2 (sub-step 2)

            # evaluate the velocities at/from (xBlob+0.5*dt*kx1,yBlob+0.5*dt*ky1)            
            xBlobTemp, yBlobTemp = xBlob + 0.5*dt*kxTemp, yBlob + 0.5*dt*kyTemp
            
            # Determine induced velocity of blobs on panels [blobs --> panels]
            # * Note: don't forget the free-stream flow.
            vxSlip, vySlip = _blobs_velocity(xBlobTemp,yBlobTemp,gBlob,self.blobs.sigma,\
                                             xEval=xCP,yEval=yCP,
                                             hardware=blobs_hardware,
                                             method=blobs_method) + vInf.reshape(2,-1) 
                                             
            # Solve for no-slip panel strengths
            self.panels.solve(vxSlip, vySlip)     

            # Evaluate velocity [  --> blobs]
            
            # Determine the induced velocity of panels on blobs [panels -- > blobs]
            vxPanel, vyPanel = self.panels.evaluteVelocity(xBlobTemp,yBlobTemp)                                        
            
            # Determine the self-induction velocity of blobs [blobs --> blobs]
            vxBlob,vyBlob = _blobs_velocity(xBlobTemp,yBlobTemp,gBlob,self.blobs.sigma,\
                                            xEval=None,yEval=None,
                                            hardware=blobs_hardware,
                                            method=blobs_method)
            # compute kx2, kx2 using (3)
            kxTemp = vxBlob + vxPanel + vInf[0] # compute kx2 using (3)
            kyTemp = vyBlob + vxPanel + vInf[1] # compute ky2 using (3)
            
            # update the positions using (step 2)
            xBlobNew += 2.0*kxTemp # kxTemp is now kx2
            yBlobNew += 2.0*kyTemp # kyTemp is now ky2
            
            #------------------------------------------------------------------

        
            #------------------------------------------------------------------        
            # STEP 3 :: x_new_3 = x_new_2 + 2*k3
        
            # evaluate the velocities at (xBlob+0.5*dt*kx2,yBlob+0.5*dt*ky2)
            xBlobTemp, yBlobTemp = xBlob + 0.5*dt*kxTemp, yBlob + 0.5*dt*kyTemp
            
            # Determine induced velocity of blobs on panels [blobs --> panels]
            # * Note: don't forget the free-stream flow.
            vxSlip, vySlip = _blobs_velocity(xBlobTemp,yBlobTemp,gBlob,self.blobs.sigma,\
                                             xEval=xCP,yEval=yCP,
                                             hardware=blobs_hardware,
                                             method=blobs_method) + vInf.reshape(2,-1)             

            # Solve for no-slip panel strengths
            self.panels.solve(vxSlip, vySlip)     

            # Evaluate velocity [  --> blobs]
            
            # Determine the induced velocity of panels on blobs [panels -- > blobs]
            vxPanel, vyPanel = self.panels.evaluteVelocity(xBlobTemp,yBlobTemp)                                        
            
            # Determine the self-induction velocity of blobs [blobs --> blobs]
            vxBlob,vyBlob = _blobs_velocity(xBlobTemp,yBlobTemp,gBlob,self.blobs.sigma,\
                                            xEval=None,yEval=None,
                                            hardware=blobs_hardware,
                                            method=blobs_method)
            
            # compute kx3 and kx3 using (4)
            kxTemp = vxBlob + vxPanel + vInf[0] # compute kx3 using (4)
            kyTemp = vyBlob + vxPanel + vInf[1] # compute ky3 using (4)
            
            # uptade the positions using (step 3)
            xBlobNew += 2.0*kxTemp # kxTemp is now kx3
            yBlobNew += 2.0*kyTemp # kyTemp is now ky3
        
            #------------------------------------------------------------------
        
        
            #------------------------------------------------------------------
            # STEP 4 :: x_new_4 = x_new_3 + k4
            
            # compute kx4 and kx4 using (5)
            vx,vy = velocity(xBlob + dt*kxTemp,yBlob + dt*kyTemp,wBlob,sigma,k=k,kernel=kernel,\
                     xEval=None,yEval=None,hardware=hardware,blocksize=blocksize,method=method) # evaluate the velocities at (xBlob+0.5*dt*kx1,yBlob+0.5*dt*ky1)
                                                                                                # note that kx3 = kxTemp and ky3 = kyTemp at this stage
            kxTemp = vx + vInf[0] # compute kx3 using (5)
            kyTemp = vy + vInf[1] # compute ky3 using (5)
            
            # uptade the positions using (step 3)
            xBlobNew += kxTemp # kxTemp is now kx4
            yBlobNew += kyTemp # kyTemp is now ky4
        
        
            # STEP 5 :: x_new_5 = (1/6)*dt*x_new_4 + x_old
            
            # x_new_5a = (1/6)*dt*x_new_4
            xBlobNew *= (1.0/6.0) * dt
            yBlobNew *= (1.0/6.0) * dt
            
            # x_new_5b = x_new_5a + x_old
            xBlobNew += xBlob 
            yBlobNew += yBlob
                
                
            # circulations are conserved in the convection step
            wBlobNew = wBlob.copy()
            
            # Return new blob positions
            return xBlobNew,yBlobNew,wBlobNew    
    
    
    
            
    def __set(self,varName,var):
        """
        Function to set parameters.
        """
        
        if varName == 'blobs':
            if type(var) != _VortexBlobs:
                raise ValueError("'blobs' should be of type pHyFlow.vortex.VortexBlobs. It is %s" % type(var))
                
            self.blobs = var
            self.blobs.evolve = self.__mute_blobs_evolve
            
        elif varName == 'panels':
            if type(var) != _Panels:
                raise ValueError("'blobs' should be of type pHyFlow.vortex.VortexBlobs. It is %s" % type(var))
                
            self.panels = var
            
            
    def __mute_blobs_evolve(self):
        raise AttributeError('Not possible with vortexPanel coupled class !')
            
          
    tStep = property(fget = lambda self: self.blobs.tStep)
    t = property(fget = lambda self: self.blobs.t)
