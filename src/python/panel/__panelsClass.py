"""

PANELS
======

The main class for the panel method. This class contains of the functions
related to the calculation of the panel method.


Class structure
---------------
Panels:

#. __init__
#. solve
#. updateBody
#. evaluateVelocity
#. plotPanels
#. plotVelocity
#. save
#. savePanels
#. saveVelocity 

:First Added:   2013-11-19
:Last Modified: 2013-11-21                                                             
:Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
:License:       GNU GPL version 3 or any later version
"""

#   Copyright (C) 2013 Lento Manickathan                                                                         
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
import scipy.linalg as splin
#from scipy.sparse.linalg import gmres, bicgstab
import types as _types

# Import pHyFlow options
from pHyFlow import options
from pHyFlow.panel.base import panelSolver as _panelSolver


# default parameters
default_velocity_computation_params = {'method':'direct', 'hardware':'cpu'}

class Panels(object):
    r"""
    The class containing all the functions related to the calculation of the
    panel method.
    
    Usage
    -----
    .. code-block:: python
    
        panelBody = panel(fileName)
        
        or
        
        panelBody = panel(xCP,yCP,xPanel,yPanel,cmGlobal,thetaLocal,externVel)
       
    Parameters
    ----------
    initFile : str, fileType (.txt)
               the file containing all the parameters necessary for the
               reinitialization of the ..py:class::panel class.

    or

    xCP : list of ndarray (float64), shape(M,)
          the :math:`x`-coordinates of the :math:`\mathbf{M}` panel collocation
          points in the `local` coordinates. In the case of multiple geometry,
          the input should be a list coordinates. `Note`: Should have a closed 
          loop. The origin of the body should be at (0,0).          

    yCP : list of ndarray (float64), shape(M,)    
    
    xPanel : list of ndarray (float64), shape(M+1,)
    
    yPanel : list of ndarray (float64), shape(M+1,)
    
    cmGlobal : ndarray (float64), shape(n,2)
               the `global` position of the panel body. `Note`: The reference
               point of the panel body is assumed to be at its (0,0).

    thetaLocal : ndarray (float64), shape(n,)
                 the `local` **anti-clockwise** rotation angle w.r.t its
                 :math:`x`-axis. The rotation is done about its `local` origin
                 at (0,0).
                 
    externVel : py:function
                reference to an external velocity fuction action on the panels.

    
    Attributes
    ----------
    ...    
    

    Methods
    -------
    solve()
        solve the panel strength to satisfy no-slip b.c.
        
    updateBody(thetaLocals, cmGlobals)
        update all the panel body coordinates.

    evaluateVelocity(xTarget,yTarget,totalVelocity=False)
        evaluate the induced velocity due to the panels.
    
    plotPanels(plot=True, save=False)
        plot and/or save panel data
    
    plotVelocity(**kwargs)
        plot and/or save induced velocity
        
    save(dir)
        save all the ..py:class::`panels` data to a file.
        
    savePanels(dir)
        save all the panel data
        
    saveVelocity(**kwargs)
        save velocity field data
    
                
    Returns
    -------
    panels : py:class
             the single or multi-body class for the calculation of the panel
             bodies.
             

    :First Added:   2013-11-21
    :Last Modified: 2013-11-22
    :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
    :Licence:       GNU GPL version 3 or any later version               
                  
    """
    """
    Revisions
    ---------
    2013-11-22, Lento Manickathan
        - Init using data file fixed.


    """
    
    def __init__(self,externVel,initFile=None,
                 velocityComputationParams=default_velocity_computation_params,
                 **kwargs):

        # Check input parameters

        # The file containing all the data for the init
        self.__set('initFile',initFile)
        
        # Assign externel velocity
        self.__set('externVel',externVel) # Freestream + arbitrary
        
        # Assign velocity computation parameters
        self.__set('velocityComputationParams',velocityComputationParams)
        
        # Computation parameters
        #self._hardware  = options.CPU_HARDWARE  # CPU
        #self._method    = options.DIRECT_METHOD # DIRECT METHOD
                
        # The important parameters (the parameters that needed to be loaded)       
        inputParameters = ['xCP','yCP','xPanel','yPanel','cmGlobal','thetaLocal']
        
        # If init file is not given, explicitly define the parameters
        if self.__initFile is None:
            
            # Number of panel bodies
            self.__nBodies = len(kwargs['xCP'])
            
            # Number of panels per body
            self.__nPanels = _numpy.zeros(self.__nBodies)
            for i,xCP in enumerate(kwargs['xCP']):
                self.__nPanels[i] =xCP.shape[0] # Determine number of panels per body
            
            # Initialize the parameters
            for parameter in inputParameters:
                self.__set(parameter, kwargs[parameter])
                #setattr(self, '_' + parameter, np.array(kwargs[parameter]))
               
        # If init file is given
        else:
            raise NotImplementedError('Not implemented !')
            #            # Load the file
            #            dataFile = np.load(self._initFile)            
            #            
            #            # Iterative load the parameters
            #            for parameter in self._inputParameters:
            #                setattr(self, '_' + parameter, dataFile[parameter]) # store and hide parameters
        
        
        # -------------------
        # Setup panel bodies


        # Panel parameters
        self.__nPanelsTotal = _numpy.int64(self.__nPanels.sum()) # Total number of panels
        self.__index = _numpy.cumsum(_numpy.hstack((0,self.__nPanels))) # Panel index
        
        # Update the rotation matrix
        self.__updateRotMat() # function of thetaLocal
              
        # Initialize the global panel parameters
        self.__norm = _numpy.zeros((2,self.__nPanelsTotal))
        self.__tang = _numpy.zeros((2,self.__nPanelsTotal))               
        self.__xyCP_global = _numpy.zeros((2,self.__nPanelsTotal))
        self.__cosSinAlpha = _numpy.zeros((2,self.__nPanelsTotal))
        self.__xyPanelStart_global = _numpy.zeros((2,self.__nPanelsTotal))
        self.__xyPanelEnd_global   = _numpy.zeros((2,self.__nPanelsTotal))
        
        # Determine the global coordinates (also, concatenate the data)
        self.__updateCoordinates()        
        
        # Assemble influence Matrix
        self.__A = _numpy.zeros((self.__nPanelsTotal,self.__nPanelsTotal))
        self.__assembleInfluenceMatrix(assembleAll=True)



    def solve(self):
        """
        Solve the panel strengths
        """
        
        # Determine the external induced velocity on the panel collocation points
        vx,vy = self.__externVel(self.__xyCP_global[0],self.__xyCP_global[1]) #TODO: not implemented
        
        # Determine the RHS
        RHS = (-vx)*self.__tang[0] + (-vy)*self.__tang[1] # vortex panel
        
        # Solve for panel strength
        self.__sPanel = splin.lu_solve(self.__LU_FACTOR, RHS)
                
        
        
    def updateBody(self,cmGlobal,thetaLocal):
        """
        Update all the panel body coordinates. This internally calculate the
        new global panel coordinates of `xPanel, yPanel, xCP, yCP` and will
        also rebuild the inter-induction influence matrix `A`
        
        Usage
        -----
        ...
        
        Parameters
        ----------
        ...
        
        Assigns
        -------
        ...
        
        Returns
        -------
        ...
        
         
        :First Added:   2013-11-22
        :Last Modified: 2013-11-22
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version            
        
        """
        # Assign the new body location
        self.__set('thetaLocal', thetaLocal)
        self.__set('cmGlobal', cmGlobal)
        #self.__thetaLocal = _numpy.copy(thetaLocal)
        #self.__cmGlobal   = _numpy.copy(cmGlobal)
        
        # Update the rotational matrix
        self.__updateRotMat()
        
        # Update the panel coordinates and angles
        self.__updateCoordinates()
        
        # Re-assemble the influence matrix A
        self.__assembleInfluenceMatrix(assembleAll=True) 
        # Implement the smarted method later, i.e partial assembly

        
        
    def evaluateVelocity(self,xTarget, yTarget,addExternVel=False):
        """
        Evaluate the velocity field on the given target locations.
        """
        
        if self.__velocityComputationParams['hardware'] == 'gpu':
            hardware = options.GPU_HARDWARE
        else:
            hardware = options.CPU_HARDWARE
        
        if self.__velocityComputationParams['method'] == 'fmm':
            method = options.FMM_METHOD
        else:
            method = options.DIRECT_METHOD
            
        # Calculate the induced velocity
        vx,vy = _panelSolver.inducedVelocity(self.__sPanel,self.__xyPanelStart_global[0],
                                             self.__xyPanelStart_global[1],self.__xyPanelEnd_global[0],
                                             self.__xyPanelEnd_global[1],self.__cosSinAlpha[0],self.__cosSinAlpha[1],
                                             xTarget,yTarget,hardware,method)

        # Add the external velocity                                             
        if addExternVel is True:
            # Determine the external induced velocity
            vxExtern,vyExtern = self.__externVel(xTarget,yTarget) #TODO: not implemented
            # Total velocity field
            vx += vxExtern
            vy += vyExtern 
        
        return vx,vy

                                             
        
    def save(self,fileName=None):
        """
        Save all the data of the ..py:class::`panels` class to a file. This
        file can be used later to re-initalize the class.
        
        Usage
        -----
        ...
        

        Parameters
        ----------
        ...
        

        Assigns
        -------
        ...
        
        
        Returns
        -------
        ...
        
        
        :First Added:   2013-11-21
        :Last Modified: 2013-11-22
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version                     
        
        """
        
        # Define the file store location            
        if fileName is None:
            fileName = './panelData'
            
        # Save data
        _numpy.savez_compressed(fileName,xCP=self.__xCP,yCP=self.__yCP,
                            xPanel=self.__xPanel,yPanel=self.__yPanel,
                            cmGlobal=self.__cmGlobal,thetaLocal=self.__thetaLocal,
                            nPanels=self.__nPanels,nBody=self.__nBodies)

        

    def __updateCoordinates(self):
        """
        Function to update the coordinates. Calculates the global
        collocation, panel coordinates. Also calculate the normal, tangent and 
        the panel angle functions.
        
        Usage
        -----
        ...
        
        Parameters
        ----------
        ...
        
        Assigns
        -------
        ...
        
        Returns
        -------
        -
      
        :First Added:   2013-11-25
        :Last Modified: 2013-11-25
        :Copyright:     Copyright (C) 2013 Lento Manickathan **pHyFlow**
        :Licence:       GNU GPL version 3 or any later version  
        
        """
        
        # Update collocation points
        for i,(x,y) in enumerate(zip(self.__xCP, self.__yCP)):
            
            # Start and end index
            iS, iE = self.__index[i], self.__index[i+1] 
            # Concatenate the global collocation coordinates
            self.__xyCP_global[:,iS:iE] = _numpy.dot(self.__rotMat[i], _numpy.vstack((x,y))) + self.__cmGlobal[i].reshape(2,1)
        
        # Update panel coordinates
        for i,(x,y) in enumerate(zip(self.__xPanel, self.__yPanel)):
            
            # Start and end index
            iS, iE = self.__index[i], self.__index[i+1]
            
            # New global panel coordinates
            xyPanelNew = _numpy.dot(self.__rotMat[i], _numpy.vstack((x,y))) + self.__cmGlobal[i].reshape(2,1)
            
            # Separate and concatenate the (start) and (end) points of the panels
            self.__xyPanelStart_global[:,iS:iE]  = xyPanelNew[:,:-1]
            self.__xyPanelEnd_global[:,iS:iE]    = xyPanelNew[:,1:]
            
            # Calculate the angle functions
            
            # Determine the length of the panels
            Lxy = xyPanelNew[:,1:] - xyPanelNew[:,:-1]
            L  = _numpy.sqrt(_numpy.sum(Lxy*Lxy,axis=0))
            
            # Panel Angles
            self.__cosSinAlpha[:,iS:iE] = Lxy / L # cosAlpha, sinAlpha
            
            # Unit Vectors (Tangent and Normal vector)
            self.__norm[:,iS:iE] = _numpy.vstack((-self.__cosSinAlpha[1,iS:iE], self.__cosSinAlpha[0,iS:iE])) 
            self.__tang[:,iS:iE] = _numpy.vstack((self.__cosSinAlpha[0,iS:iE], self.__cosSinAlpha[1,iS:iE]))
        
        
        
    def __updateRotMat(self):
        """
        Calculates the global tranformation matrix.
        """        
        # The global rotational matrix
        self.__rotMat = [_numpy.array([[_numpy.cos(theta), -_numpy.sin(theta)],
                                  [_numpy.sin(theta), _numpy.cos(theta)]])
                        for theta in self.__thetaLocal.flatten()]
                                         


    def __assembleInfluenceMatrix(self,assembleAll=True):
        """
        Assemble the influence matrix.       
        """
        
        if self.__velocityComputationParams['hardware'] == 'gpu':
            hardware = options.GPU_HARDWARE
        else:
            hardware = options.CPU_HARDWARE
        
        if self.__velocityComputationParams['method'] == 'fmm':
            method = options.FMM_METHOD
        else:
            method = options.DIRECT_METHOD
            
        
        if assembleAll:
            
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
        self.__LU_FACTOR = splin.lu_factor(self.__A)
        
  


    def __set(self,varName,var):
        """
        Function to check the input parameter
        """
        
        # Check/Set Init File
        if varName == 'initFile':
            if type(var) == str or var == None:
                self.__initFile = var
            else:
                raise TypeError('initFile must be a str. It is %s.' % str(type(var)))
        
        # Check/Set externVel
        elif varName == 'externVel':
            if type(var) != _types.FunctionType:
                raise TypeError('externVel must be a function. It is %s.' % str(type(var)))
            self.__externVel = var
                
        elif varName == 'velocityComputationParams':
            
            # check first the biot-savart computation option
            if var['method'] not in options.biot_savart_options:
                raise ValueError(('velocityComputationParams[\'mehtod\'] must be one of ' +\
                                 '\'%s\','*(len(options.biot_savart_options)-1) +\
                                 '\'%s\'.' + ' It is %s.') % (options.biot_savart_options + (str(var['method']),)))
            # check the hardware option
            if var['hardware'] not in options.hardware_options:
                raise ValueError(('velocityComputationParams[\'hardware\'] must be one of ' +\
                                 '\'%s\','*(len(options.hardware_options)-1) +\
                                 '\'%s\'.' + ' It is %s.') % (options.hardware_options + (str(var['hardware']),)))
            self.__velocityComputationParams = var         
             
                
        # Panel variables                
        elif varName in ['xCP','yCP','xPanel','yPanel','cmGlobal','thetaLocal']:
            if len(var) != self.__nBodies:
                raise ValueError('%s must have exactly %g numpy arrays. It has %g.' % (varName, self._nBodies, len(var)))
                
            
            for k, temp_value in enumerate(var):
                if varName == 'thetaLocal':
                    if type(temp_value) != float:
                        raise TypeError('%s[%g] must be a float. It is %s.') % (varName, k, str(type(temp_value)))     
                    else:
                        self.__thetaLocal = _numpy.array(var)
                        
                else:
                    if type(temp_value) != _numpy.ndarray:
                        raise TypeError('%s[%g] must be a numpy ndarray. It is %s.') % (varName, k, str(type(temp_value)))
                    if temp_value.dtype != 'float64':
                            raise TypeError('%s[%g] must have dtype float64. It is %s.') % (varName, k, str(temp_value.dtype))
                    
                    if varName == 'cmGlobal':
                        if temp_value.shape != (2,):
                            raise ValueError('%s[%g] must have shape (2,0). It has shape %s.') % (varName, k, str(temp_value.shape))
                        else:
                            self.__cmGlobal = _numpy.array(var)
                            
                    elif varName == 'xPanel' or varName == 'yPanel':
                        if temp_value.shape !=(self.__nPanels[k]+1,):
                            raise ValueError('%s[%g] must have shape (%g,0). It has shape %s.') % (varName, k, str((self.__nPanels[k]+1,)), str(temp_value.shape))
                        if varName == 'xPanel':
                            self.__xPanel = _numpy.array(var)
                        else:
                            self.__yPanel = _numpy.array(var)
                            
                    elif varName == 'xCP' or varName == 'yCP':
                        if temp_value.shape !=(self.__nPanels[k],):
                            raise ValueError('%s[%g] must have shape (%g,0). It has shape %s.') % (varName, k, str((self.__nPanels[k],)), str(temp_value.shape))
                        if varName == 'xCP':
                            self.__xCP = _numpy.array(var)
                        else:
                            self.__yCP = _numpy.array(var)
                
                
    @property  
    def sPanel(self):
        """
        The new strength of the panel, satisfying the no-through boundary
        condition fo body. 
        """
        
        __sPanel = _numpy.copy(self.__sPanel)
        
        # Separate panel strengths according to the body
        sPanel = [] # Init list
        
        # Split the data
        for i in range(self.__nBodies):
            iS,iE = self.__index[i], self.__index[i+1]
            sPanel.append(__sPanel[iS:iE])
        
        # Return list of panel strengths
        return sPanel
        
    @property
    def xCPLocal(self):
        """
        The local :math:`x`-coordinates of the panel collocation points,
        split into sub-arrays according to the bodies.  
        """
        return self.__xCP
        
        
    @property
    def xCPGlobal(self):
        """
        The global :math:`x`-coordinates of the panel collocation points.
        
        Note: 
        
           *Slower* function. This function contains a for-loop for the 
            splitting.
        """
        
        xCPGlobal = []
        # Split the data
        for i in range(self.__nBodies):
            iS,iE = self.__index[i], self.__index[i+1]
            xCPGlobal.append(self.__xyCP_global[0,iS:iE])
        
        # return the list of xCPGlobal
        return xCPGlobal
        
    @property
    def xCPGlobalCat(self):
        """
        The global :math:`x`-coordinates of the panel collocation points.
        """
        return self.__xyCP_global[0]        
        
    @property
    def yCPLocal(self):
        """
        The local :math:`y`-coordinates of the panel collocation points,
        split into sub-arrays according to the bodies.  
        """
        return self.__yCP
        
    @property
    def yCPGlobal(self):
        """
        The global :math:`y`-coordinates of the panel collocation points.
        
        Note: 
        
           *Slower* function. This function contains a for-loop for the 
            splitting.
        """
        
        yCPGlobal = []
        # Split the data
        for i in range(self.__nBodies):
            iS,iE = self.__index[i], self.__index[i+1]
            yCPGlobal.append(self.__xyCP_global[1,iS:iE])
        
        # return the list of xCPGlobal
        return yCPGlobal
        
    @property
    def yCPGlobalCat(self):
        """
        The global :math:`y`-coordinates of the panel collocation points.
        """
        return self.__xyCP_global[1] 
        
    @property
    def xPanelLocal(self):
        """
        The local :math:`x`-coordinate of the panel edge points.
        
        Note: It is a closed loop.
        """
        return _numpy.copy(self.__xPanel)
        
    @property
    def xPanelGlobal(self):
        """
        The global :math:`x`-coordinate of the panel edge points.
        
            Note: It is a open loop
        """
        __xPanelGlobal = _numpy.copy(self.__xyPanelStart_global[0])
        
        xPanelGlobal = []
        # Split the data
        for i in range(self.__nBodies):
            iS,iE = self.__index[i], self.__index[i+1]
            xPanelGlobal.append(__xPanelGlobal[iS:iE])
            
        return xPanelGlobal
        
    @property
    def yPanel(self):
        """
        The local :math:`x`-coordinate of the panel edge points.
        
            Note: It is a closed loop.
        """
        return _numpy.copy(self.__yPanel)        

    @property
    def yPanelGlobal(self):
        """
        The global :math:`x`-coordinate of the panel edge points.
        
            Note: It is a open loop
        """
        __yPanelGlobal = self.__xyPanelStart_global[1]
        
        yPanelGlobal = []
        # Split the data
        for i in range(self.__nBodies):
            iS,iE = self.__index[i], self.__index[i+1]
            yPanelGlobal.append(__yPanelGlobal[iS:iE])
            
        return yPanelGlobal
        
    @property
    def cmGlobal(self):
        """
        The position of the reference ppoint of the given panel bodies.
        """
        return self.__cmGlobal
        
    @property
    def thetaLocal(self):
        """
        The rotational angles of the apnel bodies w.r.t to the global
        :math:`x`-axis.
        """
        return self.__thetaLocal
        
    @property
    def nBodies(self):
        """
        Number of panel bodies.
        """
        return self.__nBodies
    
    @property
    def nPanels(self):
        """
        Number of panels in each panel body.
        """
        return self.__nPanels
        
    @property
    def nPanelsTotal(self):
        """
        Number of panels in total.
        """
        return self.__nPanelsTotal   

    @property
    def A(self):
        """
        The panel self-induction matrix.
        """
        return self.__A
        