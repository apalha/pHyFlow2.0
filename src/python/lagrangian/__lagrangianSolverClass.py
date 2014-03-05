#-*- coding: utf-8 -*-
__doc__ = """

LAGRANGIANSOLVER
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

__all__ = ['LagrangianSolver']

# External packages
import numpy as _numpy

# pHyFlow modules
from pHyFlow import options as _pHyFlowOptions
from pHyFlow.aux.customDecorators import simpleGetProperty

# pHyFlow blob modules
from pHyFlow.blobs import Blobs as _BlobsClass
from pHyFlow.blobs.base.induced import velocity as _blobs_velocity

# phyFlow panel modules
from pHyFlow.panels import Panels as _PanelsClass

# pHyFlow lagrangian options
from pHyFlow.lagrangian import lagrangianOptions

class LagrangianSolver(object):
    r"""
    
    Usage
    -----
    .. code-block:: python
    
        LagrangianSolver(blobs, panels, 
                         couplingParams={'panelStrengthUpdate': 'constant'})
        
    Parameters
    ----------
    blobs : pHyFlow.blobs.Blobs
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
    blobs : pHyFlow.blobs.Blobs
            the vortex-blobs class is assigned
            
    deltaT
    
    panels : pHyFlow.panels.Panels
             the panel class    
             
    t
    tStep
    vInf    

    __couplingParams : dict
                       Dictionary containing parameters for coupling the
                       panels and the blobs together.  
    
    Methods
    -------
    evolve
    evaluateVelocity
    __coupled_convection
    __advanceTime
    __set
    __mute_blobs_evolve
    
    :First Added:   2014-02-21
    :Last Modified: 2014-02-25                         
    :Copyright:     Copyright (C) 2014 Lento Manickathan **pHyFlow**
    :License:       GNU GPL version 3 or any later version   
   
    """
    """
    Revisions
    ---------


    """
    
    def __init__(self,blobs,panels=None,couplingParams={
                 'panelStrengthUpdate': lagrangianOptions.PANEL_STRENGTH_UPDATE['default']}):
        
        
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
        if self.panels is not None:
            self.blobs.evolve = self.__mute_blobs_evolve # remove blobs.evolve option
        
            # Solve for the initial panel strengths
            self.__solve_panelStrength()
       
       
    def evaluateVelocity(self,xTarget,yTarget,Ndirect=35,tol=1.0e-6,cutoff=None):
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
        # Determine the induced velocity of blobs and free-stream
        if self.panels is None:
            # Induced velocity
            vx,vy = self.blobs.evaluateVelocity(xTarget,yTarget,Ndirect,tol,cutoff)
            
            # Total induced velocity
            vx += self.blobs.vInf[0]
            vy += self.blobs.vInf[1]
        
        # Determine the induced velocity of blobs, panels and free-stream                
        else:
            #            vx,vy = self.blobs.evaluateVelocity(xTarget,yTarget,Ndirect,tol,cutoff) \
            #                            + self.panels.evaluateVelocity(xTarget,yTarget) \
            #                            + self.blobs.vInf.reshape(2,-1)
            # Induced velocity from blobs
            vxBlob, vyBlob = self.blobs.evaluateVelocity(xTarget,yTarget,Ndirect,tol,cutoff)
            # Induced velocity from panels
            vxPanel, vyPanel = self.panels.evaluateVelocity(xTarget,yTarget)
            
            # Total induced velocity
            vx = vxBlob + vxPanel + self.blobs.vInf[0]
            vy = vyBlob + vyPanel + self.blobs.vInf[1]
        
        # return the induced velocity
        return vx,vy
                
                
    def evolve(self):
        r"""
        Function to evolve the coupled vortex-panel problem. During the evolution
        of the vortex blobs, the no-through boundary condition is taken into 
        account by solving the no-slip panel problem.
        
        * Note: moving panels not implemented yet !
        
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
        blobs : pHyFlow.blobs.Blobs
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
        # Only convect the blobs
        if self.panels is None:
            
            # Evolve : convection + diffusion + redistribution + population control
            self.blobs.evolve()
            
        # Coupled convection : blobs with coupled panels            
        else:
            
            #----------------------------------------------------------------------
            # Evolve: Convect + Diffuse the blobs with regard of panel geometry
      
             # Algorithm for convecting the coupled problem     
            self.__coupled_convection()
            
            # Diffusion Step : diffuse vortex blobs
            if self.blobs.stepDiffusion != 0: # corresponds to nu = 0.0, therefore no diffusion is to be computed
                if (self.tStep % self.blobs.stepDiffusion) == 0: # if the time step is a multiple of the stepDiffusion perform diffusion
                    self.blobs._Blobs__diffusion()
    
            # update the time counter
            self.__advanceTime() # update both blobs+panels internal time.

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
        blobs : pHyFlow.blobs.Blobs
                the position (only) of the blobs are updated
                
                x : numpy.ndarray(float64), shape (blobs.numBlobs, )
                    the :math:`x`-position of the particle.
                y : numpy.ndarray(float64), shape (blobs.numBlobs, )
                    the :math:`y`-postion of the particles.
                    
        panels : pHyFlow.panels.Panels
                 the strength of the panels are updated.
                 
                 sPanel : numpy.ndarray(float64), shape (panels.nPanelTotal, )
                          the strength :math:`\gamma` of the panels.
                  
                 * Note: position update not implemented yet.        
        
        :First Added:   2014-02-21
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version        
        
        """      

        #----------------------------------------------------------------------
        # Determine parameters

        # convert the hardware flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['hardware'] == 'gpu': 
            blobs_hardware = _pHyFlowOptions.GPU_HARDWARE
        else: 
            blobs_hardware = _pHyFlowOptions.CPU_HARDWARE

        # convert the method flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['method'] == 'fmm': 
            blobs_method = _pHyFlowOptions.FMM_METHOD
        else: 
            blobs_method = _pHyFlowOptions.DIRECT_METHOD
    
        # convert the time integration method into an int to use in _base_convection
        if self.blobs.timeIntegrationParams['method'] == 'rk4': 
            blobs_integrator = _pHyFlowOptions.RK4_INTEGRATOR
        elif self.blobs.timeIntegrationParams['method'] == 'euler':
            blobs_integrator = _pHyFlowOptions.FE_INTEGRATOR
            
        #----------------------------------------------------------------------


        #----------------------------------------------------------------------
        # Time integrate the vortex blobs

        # If integration method is RK4
        if blobs_integrator == _pHyFlowOptions.RK4_INTEGRATOR:

            # Make references to vortex-blobs
            xBlob, yBlob, gBlob = self.blobs.x, self.blobs.y, self.blobs.g 
            
            # Make references to panel collocation points (where no-slip b.c. is enforced.)
            xCP, yCP = self.panels.xyCPGlobalCat
        
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
                self.blobs._Blobs__x += RK4_b*kxTemp*self.deltaT
                self.blobs._Blobs__y += RK4_b*kyTemp*self.deltaT                                                      
                #--------------------------------------------------------------
            
            
        elif blobs_integrator == _pHyFlowOptions.FE_INTEGRATOR:
            raise NotImplementedError('FE time integration not available !')
        
        #----------------------------------------------------------------------           

    
    def __solve_panelStrength(self):
        r"""
        Function to solve for the initial panel strength
        """
        
        #----------------------------------------------------------------------
        # Determine parameters

        # convert the hardware flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['hardware'] == 'gpu': 
            blobs_hardware = _pHyFlowOptions.GPU_HARDWARE
        else: 
            blobs_hardware = _pHyFlowOptions.CPU_HARDWARE

        # convert the method flag into an int to use in _base_convection
        if self.blobs.velocityComputationParams['method'] == 'fmm': 
            blobs_method = _pHyFlowOptions.FMM_METHOD
        else: 
            blobs_method = _pHyFlowOptions.DIRECT_METHOD
            
        #----------------------------------------------------------------------


        #----------------------------------------------------------------------
        # Solve panel strength

        # Make references to vortex-blobs
        xBlob, yBlob, gBlob = self.blobs.x, self.blobs.y, self.blobs.g 
            
        # Make references to panel collocation points (where no-slip b.c. is enforced.)
        xCP, yCP = self.panels.xyCPGlobalCat
        
        # Determine the initial slip velocity
        vxSlip, vySlip = _blobs_velocity(xBlob,yBlob,gBlob,self.blobs.sigma,
                                         xEval=xCP,yEval=yCP,hardware=blobs_hardware,   
                                         method=blobs_method) \
                                 + self.vInf.reshape(2,-1)
                
        # Solve for no-slip panel strengths
        self.panels.solve(vxSlip, vySlip)
        
                            
    def __advanceTime(self):
        """
        Function to advance time. To advance time, we advance both the internal
        time of Blob and Panels. Then we assert that both :math:`t`
        are the same. 
        
        The LagrangianSolver uses internal time of Blobs.
        
        Usage
        -----
        .. code-block :: python
        
            __advanceTime()
            
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
        :Last Modified: 2014-02-25
        :Copyright:     Copyright (C) 2014 Lento Manickathan, **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version    
        """
        # Advance the internal time of blobs
        self.blobs._advanceTime()
    
        # Advance the internal time of panels
        self.panels._advanceTime(self.deltaT)
        
        # Check if both are sync
        if _numpy.abs(self.blobs.t - self.panels.t) > _numpy.spacing(100):
            raise ValueError('The internal time of Blobs and Panels'\
                             ' are not synchronized. blobs = %g, panels = %g'\
                             % (self.blobs.t, self.panels.t)  )
                             
            
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
        blobs : pHyFlow.blobs.Blobs
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
            if type(var) != _BlobsClass:
                raise ValueError("'blobs' should be of type pHyFlow.blobs.Blobs. It is %s" % type(var))
                
            self.blobs = var
        
        # panels
        elif varName == 'panels':
            if var is not None:
                if type(var) != _PanelsClass:
                    raise ValueError("'blobs' should be of type pHyFlow.blobs.Blobs. It is %s" % type(var))
                
            self.panels = var
            
        # coupling parameters
        elif varName == 'couplingParams':
            if type(var) != dict:
                raise ValueError("'couplingParams' should be of type dict. It is %s" % type(var))
            if var['panelStrengthUpdate'] not in lagrangianOptions.PANEL_STRENGTH_UPDATE['available']:
                raise ValueError("'couplingParams['panelStrengthUpdate']' is unknown. It should be in [%s]. It is %s" % (str(lagrangianOptions.PANEL_STRENGTH_UPDATE['available']),var['panelStrengthUpdate']))
            self.__couplingParams = var
            
            
    def __mute_blobs_evolve(self):
        raise AttributeError('Not possible with LagrangianSolver coupled class !')
            

    #--------------------------------------------------------------------------
    # All attributes
            
    @simpleGetProperty
    def tStep(self):
        r"""
        tStep : int
                the current time step of the simulation
        """
        return self.blobs.tStep
            
    # Current time t
    @simpleGetProperty
    def t(self):
        r"""
        t : float
            the current time of the simulation
        """
        return self.blobs.t
    
    # Time-step size
    @simpleGetProperty
    def deltaT(self):
        r"""
        deltaT : float
                 the time step size of the simulation :math:`|Delta t`
        """
        return self.blobs.deltaTc
    
    # Free-stream velocity
    @simpleGetProperty
    def vInf(self):
        """
        vInf : numpy.ndarray(float64), shape (2,)
               the :math:`x,y` component of the free-stream velocity
        """
        return self.blobs.vInf
        
        